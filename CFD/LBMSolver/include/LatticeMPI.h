#pragma once

#include "GlobalData.h"
#include "D2Q9.h"
#include "Range.h"
#include "LBMGeometry.h"
#include "Collision.h"
#include "GenBcMap.h"
#include <mpi.h>

class DirectionMask
{
private:
  static const uint32_t directionInBinary[32];
  uint32_t directionMask;
public:
  explicit DirectionMask(const uint32_t dm = 0) : directionMask(dm) {}

  DirectionMask(const DirectionMask &dm) : directionMask(dm) {};

  operator uint32_t () const { return directionMask; }

  uint32_t *data() { return &directionMask; }
  const uint32_t *data() const { return &directionMask; }

  inline void add_direction(
    const size_t dId
  ) {
    directionMask |= directionInBinary[dId];
  }

  inline void del_direction(
    const size_t dId
  ) {
    directionMask &= ~(directionInBinary[dId]);
  }

  inline bool next_direction(
    size_t &dIdPlusOne, 
    const size_t &maxId = 32
  ) const {
    for(; dIdPlusOne < maxId; dIdPlusOne++) {
      if((directionMask & directionInBinary[dIdPlusOne]) > 0) {
        dIdPlusOne++; 
        return true;
      }
    }
    return false;
  }

  inline size_t size(
    const size_t &maxId = 32
  ) const {
    size_t counter = 0;
    for(size_t i; i < maxId; i++) {
      if((directionMask & directionInBinary[i]) > 0) {
        counter++; 
      }
    }
    return counter;
  }
};


template<const size_t dim>
nnwlbm::Range<dim> get_send_range(
  const nnwlbm::Range<dim> &sender,
  const nnwlbm::Range<dim> &receiver,
  const int &margin = 1
) {
  auto temp = nnwlbm::Range<dim>::enlarge(receiver, margin);
  return nnwlbm::Range<dim>::CalIntersection(sender, temp);
}

template<const size_t dim>
nnwlbm::Range<dim> get_recv_range(
  const nnwlbm::Range<dim> &receiver,
  const nnwlbm::Range<dim> &sender,
  const int &margin = 1
) {
  auto temp = nnwlbm::Range<dim>::enlarge(receiver, margin);
  return nnwlbm::Range<dim>::CalIntersection(sender, temp);
}

template<const size_t dim>
class RangeNeighbors
{
public:
  typedef nnwlbm::Range<dim> Neighbor;
  typedef nnwlbm::Range<dim> ReceiverBuffer;
  typedef nnwlbm::Range<dim> SenderBuffer;
  typedef std::pair<ReceiverBuffer, SenderBuffer> Buffer;

  explicit RangeNeighbors(const nnwlbm::Range<dim> &rg) : hostRange(rg) {}
  RangeNeighbors(const RangeNeighbors &rn) : hostRange(rn.hostRange), neighborAndBuffer(rn.neighborAndBuffer) {}
  bool AddNeiborCandidate(const nnwlbm::Range<dim> &neighborCandidate, const int &bufferSize = 1) 
  {
    if(hostRange.IsNeighbor(neighborCandidate, bufferSize)) {
      ReceiverBuffer rbf = get_recv_range<dim>(hostRange, neighborCandidate);
      SenderBuffer   sbf = get_send_range<dim>(hostRange, neighborCandidate);
      neighborAndBuffer.emplace(neighborCandidate, Buffer(rbf, sbf));
      return true;
    }
    return false;
  }
  const std::map<Neighbor, Buffer> &get_neighbors() const { return neighborAndBuffer; }
  const nnwlbm::Range<dim> &GetHostRange() const { return hostRange; }
private:
  std::map<Neighbor, Buffer> neighborAndBuffer;
  const nnwlbm::Range<dim> &hostRange;
};

class RequestHeap
{
protected:
  virtual bool expand();
  virtual bool shrink() { return true; }
protected:
  int maxSize;
  int increment;
  int freeIndex;
  int completeIndex;
  MPI_Request *heap;
  //public:
  int* flag;
  MPI_Status *status;
public:
  RequestHeap();
  ~RequestHeap();
  virtual void reset();
public:
  MPI_Request *get_request(int i);
  MPI_Request *GetNew_request();
  void set_increment(int n) { increment = n; }
  const MPI_Status *get_complete_status();
  int GetSize() { return freeIndex; }
  int get_completeNum() { return completeIndex; }
  int is_some_complete();
  bool test_all_complete();
  void WaitAllComplete();
};

//class EnvMpiET 
template<typename LC, template<typename C, const size_t dim, typename Compare = nnwlbm::Color, const int ndef = -99999>class ColorMaskData>
class ColorMapCommunicator 
{
  //typedef nnwlbm::MaskData<double, LC::GetDim(), nnwlbm::Color> ColorMaskData;
public:
  ColorMapCommunicator(
    ColorMaskData<double, LC::GetDim(), nnwlbm::Color> &colorIn, 
    const RangeNeighbors<LC::GetDim()> &nbs, 
    const std::map<nnwlbm::Range<LC::GetDim()>, int> &rrMap,
    MPI_Comm commIn
  ) : colorMd(colorIn), comm(commIn), neighbors(nbs), rangeToRank(rrMap)
  { allocate_buffer(); }
  bool SendAndRecv();
  inline bool is_completed() { return requestHeap.test_all_complete(); }
  void wait() { requestHeap.WaitAllComplete(); }
  bool update();
private:
  bool allocate_buffer();
private:
  RequestHeap requestHeap;
  const RangeNeighbors<LC::GetDim()> &neighbors;
  std::vector<std::vector<byte>> sendBuff;
  std::vector<std::vector<byte>> recvBuff;
  ColorMaskData<double, LC::GetDim(), nnwlbm::Color> &colorMd;
  MPI_Comm comm;
  const std::map<nnwlbm::Range<LC::GetDim()>, int> &rangeToRank;
};

template<typename LC, template<typename C, const size_t dim, typename Compare = nnwlbm::Color, const int ndef = -99999>class ColorMaskData>
bool ColorMapCommunicator<LC, ColorMaskData>::allocate_buffer() 
{
  size_t neighborNum = neighbors.get_neighbors().size();
  sendBuff.resize(neighborNum);
  recvBuff.resize(neighborNum);
  // index of buffer;
  size_t ib = 0;
  for(const auto& nb : neighbors.get_neighbors()) 
  {
    const nnwlbm::Range<LC::GetDim()> &recvRange = nb.second.first;
    const nnwlbm::Range<LC::GetDim()> &sendRange = nb.second.second;
    int buffSize;
    MPI_Pack_size(recvRange.GetSize(), MPI_INT, comm, &buffSize);
    recvBuff[ib].resize(buffSize+MPI_BSEND_OVERHEAD);
    MPI_Pack_size(sendRange.GetSize(), MPI_INT, comm, &buffSize);
    sendBuff[ib].resize(buffSize+MPI_BSEND_OVERHEAD);
  }
}

template<typename LC, template<typename C, const size_t dim, typename Compare = nnwlbm::Color, const int ndef = -99999>class ColorMaskData>
bool ColorMapCommunicator<LC, ColorMaskData>::SendAndRecv() 
{
  nnwlbm::RangeIndex<LC::GetDim()> rgIndex = colorMd.GetNodeIndex();
  // index of buffer;
  size_t ib = 0;
  for(const auto &nb : neighbors.get_neighbors()) 
  {
    const nnwlbm::Range<LC::GetDim()> &nbRange = nb.first;
    const nnwlbm::Range<LC::GetDim()> &recvRange = nb.second.first;
    const nnwlbm::Range<LC::GetDim()> &sendRange = nb.second.second;
    int loc[LC::GetDim()];
    int pos = 0;
    for(int carry = sendRange.IniCarry(loc); 
      carry == 0; carry = sendRange.AddOne(loc))
    {
      const size_t nodeId = rgIndex.CalNodeID(loc);
      int color = colorMd.GetColor(nodeId);
      MPI_Pack(&color, 1, MPI_INT, (void*)sendBuff[ib].data(), 
        sendBuff[ib].size(), &pos, comm);
    }

    auto reqSend = requestHeap.GetNew_request();
    int tag = 0 + 100 * rangeToRank.at(nnwlbm::Range<LC::GetDim()>::enlarge(colorMd.GetRange(), -1));
    MPI_Isend((void*)sendBuff[ib].data(), sendRange.GetSize(),
      MPI_INT, rangeToRank.at(nbRange), tag, comm, reqSend);

    auto reqRecv = requestHeap.GetNew_request();
    tag = 0 + 100 * rangeToRank.at(nbRange);
    MPI_Irecv((void*)recvBuff[ib].data(), recvRange.GetSize(),
      MPI_INT, rangeToRank.at(nbRange), tag, comm, reqRecv);

    ib++;
  }
  return true;
}

template<typename LC, template<typename C, const size_t dim, typename Compare = nnwlbm::Color, const int ndef = -99999>class ColorMaskData>
bool ColorMapCommunicator<LC, ColorMaskData>::update() 
{
  if(!is_completed()) return false;
  nnwlbm::RangeIndex<LC::GetDim()> rgIndex = colorMd.GetNodeIndex();
  // index of buffer;
  size_t ib = 0;
  for(const auto &nb : neighbors.get_neighbors()) 
  {
    const nnwlbm::Range<LC::GetDim()> &nbRange = nb.first;
    const nnwlbm::Range<LC::GetDim()> &recvRange = nb.second.first;
    //const nnwlbm::Range<LC::GetDim()>& sendRange = nb.second.second;
    int loc[LC::GetDim()];
    int pos = 0;
    auto buff = (const int*)recvBuff[ib].data();
    for(int carry = recvRange.IniCarry(loc); 
      carry == 0; carry = recvRange.AddOne(loc))
    {
      const size_t nodeId = rgIndex.CalNodeID(loc);
      colorMd.SetColor(nodeId, buff[pos]);
      pos++;
    }
    ib++;
  }
  return true;
}

template<typename LC, template<typename C, const size_t dim, typename Compare = nnwlbm::Color, const int ndef = -99999>class PdfMaskData>
class PdfCommunicator 
{
public:
  PdfCommunicator(
    PdfMaskData<double, LC::GetDim(), nnwlbm::Color> *pdfIn,
    const RangeNeighbors<LC::GetDim()>& nbs, 
    const std::map<nnwlbm::Range<LC::GetDim()>, int> &rrMap,
    MPI_Comm commIn
  ) : pdfData(pdfIn), comm(commIn), neighbors(nbs), rangeToRank(rrMap) 
  { allocate_buffer(); }
  void SetPdfDataWithIdenticalMask(PdfMaskData<double, LC::GetDim(), nnwlbm::Color> *pdfIn) { pdfData = pdfIn; }
  bool SendAndRecv();
  inline bool is_completed() { return requestHeap.test_all_complete(); }
  void wait() { requestHeap.WaitAllComplete(); }
  bool update();
private:
  bool allocate_buffer();
private:
  RequestHeap requestHeap;
  RangeNeighbors<LC::GetDim()> neighbors;
  std::vector<std::vector<double>> sendBuff;
  std::vector<std::vector<double>> recvBuff;
  PdfMaskData<double, LC::GetDim(), nnwlbm::Color>* pdfData;
  MPI_Comm comm;
  std::vector<std::vector<std::pair<size_t, DirectionMask>>> recvIndex;
  std::vector<std::vector<std::pair<size_t, DirectionMask>>> sendIndex;
  const std::map<nnwlbm::Range<LC::GetDim()>, int> &rangeToRank;
};

template<typename LC, template<typename C, const size_t dim, typename Compare = nnwlbm::Color, const int ndef = -99999>class PdfMaskData>
bool PdfCommunicator<LC, PdfMaskData>::allocate_buffer() 
{
  auto hr = neighbors.GetHostRange();
  size_t neighborNum = neighbors.get_neighbors().size();
  sendBuff.resize(neighborNum);
  recvBuff.resize(neighborNum);
  sendIndex.resize(neighborNum);
  recvIndex.resize(neighborNum);

  size_t ib = 0;

  auto hostRange = neighbors.GetHostRange();
  nnwlbm::RangeIndex<LC::GetDim()> rgIndex = pdfData->GetNodeIndex();
  for(const auto &nb : neighbors.get_neighbors()) 
  {
    const nnwlbm::Range<LC::GetDim()> &nbRange = nb.first;
    const nnwlbm::Range<LC::GetDim()> &recvRange = nb.second.first;
    const nnwlbm::Range<LC::GetDim()> &sendRange = nb.second.second;

    int sendVelNum = 0;
    int recvVelNum = 0;

    int loc[LC::GetDim()];

    // recvBuff
    for(int carry = recvRange.IniCarry(loc); carry == 0;
      carry = recvRange.AddOne(loc))
    {
      const size_t nodeId = rgIndex.CalNodeID(loc);
      if(pdfData->IsValid(nodeId)) {
        int locN[LC::GetDim()];
        DirectionMask dm;
        for(size_t d = 0; d < LC::size(); d++) {
          for(size_t iDim = 0; iDim < LC::GetDim(); iDim++) {
            locN[iDim] = loc[iDim] + LC::direction(d)[iDim];
          }
          if(hostRange.IsContaining(locN)) {
            const size_t nodeIdN = rgIndex.CalNodeID(locN);
            if(pdfData->IsValid(nodeIdN)) {
              recvVelNum++;
              dm.add_direction(d);
            }
          }
        }
        if(0 != dm ) {
          recvIndex[ib].push_back(std::make_pair(nodeId, dm));
        }
      }
    }
    recvBuff[ib].resize(recvVelNum);

    // sendBuff
    for(int carry = sendRange.IniCarry(loc); carry == 0;
      carry = sendRange.AddOne(loc))
    {
      const size_t nodeId = rgIndex.CalNodeID(loc);
      if(pdfData->IsValid(nodeId)) {
        int locN[LC::GetDim()];
        DirectionMask dm;
        for(size_t d = 0; d < LC::size(); d++) {
          for(size_t iDim = 0; iDim < LC::GetDim(); iDim++) {
            locN[iDim] = loc[iDim] + LC::direction(d)[iDim];
          }
          if(nbRange.IsContaining(locN)) {
            const size_t nodeIdN = rgIndex.CalNodeID(locN);
            if(pdfData->IsValid(nodeIdN)) {
              sendVelNum++;
              dm.add_direction(d);
            }
          }
        }
        if(0 != dm) sendIndex[ib].push_back(std::make_pair(nodeId, dm));
      }
    }
    sendBuff[ib].resize(sendVelNum);
    ib++;
  }
  return true;
}

template<typename LC, template<typename C, const size_t dim, typename Compare = nnwlbm::Color, const int ndef = -99999>class PdfMaskData>
bool PdfCommunicator<LC, PdfMaskData>::SendAndRecv() 
{
  auto hr = neighbors.GetHostRange();
  auto rgIndex = pdfData->GetNodeIndex();

  size_t ib = 0;
  for(const auto &nb : neighbors.get_neighbors()) 
  {
    const nnwlbm::Range<LC::GetDim()> &nbRange = nb.first;
    const nnwlbm::Range<LC::GetDim()> &recvRange = nb.second.first;
    const nnwlbm::Range<LC::GetDim()> &sendRange = nb.second.second;

    double *buf = (sendBuff[ib].data());
    int pos = 0;
    for(const auto &nodeIdAndDirectionMask : sendIndex[ib])
    {
      auto nodeId = nodeIdAndDirectionMask.first;
      auto dm = nodeIdAndDirectionMask.second;
      
      const double* pdf = pdfData->GetVar(nodeId);
      for(size_t dPos = 0; dm.next_direction(dPos, LC::size()); )
      {
        //MPI_Pack((void*)(pdf+dPos-1), 1, MPI_DOUBLE, buf, (int)sendBuff[ib].size(), &pos, comm);
        buf[pos] = pdf[dPos - 1];
        pos++;
      }
    }
    auto req = requestHeap.GetNew_request();
    int tag = 1 + 100 * rangeToRank.at(neighbors.GetHostRange());
    MPI_Isend((void*)sendBuff[ib].data(), (int)sendBuff[ib].size(), 
      MPI_DOUBLE, rangeToRank.at(nbRange), tag, comm, req);
    MPI_Request_free(req);

    tag = 1 + 100 * rangeToRank.at(nbRange);
    MPI_Irecv((void*)recvBuff[ib].data(), (int)recvBuff[ib].size(),
      MPI_DOUBLE, rangeToRank.at(nbRange), tag, comm, req);

    ib++;
  }
  return true;
}

template<typename LC, template<typename C, const size_t dim, typename Compare = nnwlbm::Color, const int ndef = -99999>class PdfMaskData>
bool PdfCommunicator<LC, PdfMaskData>::update() 
{
  if(!is_completed()) return false;
  nnwlbm::RangeIndex<LC::GetDim()> rgIndex = pdfData->GetNodeIndex();
  size_t ib = 0;
  for(const auto &nb : neighbors.get_neighbors()) 
  {
    int pos = 0;
    auto buff = (const double*)recvBuff[ib].data();
    for(const auto &nodeIdAndDirectionMask : recvIndex[ib])
    {
      auto nodeId = nodeIdAndDirectionMask.first;
      auto dm = nodeIdAndDirectionMask.second;
      double *pdf = pdfData->GetVar(nodeId);
      for(size_t dPos = 0; dm.next_direction(dPos, LC::size()); )
      {
        pdf[dPos - 1] = buff[pos];
        pos++;
      }
    }
    ib++;
  }
  return true;
}

template<typename Output, typename MainLoop>
void time_marching_mpi(
  Output &output,
  MainLoop &mainLoop,
  const size_t &maxStep,
  const size_t &outputInterval,
  const size_t &reportInterval
) {
  int rank, rankSize;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &rankSize);
  for(int iTime = 0; iTime < maxStep; iTime++) {
    if (iTime % outputInterval == 0) {
      //output.SetFileName(std::to_string(iTime)+"rank"+std::to_string(rank));
      output.SetFileName("rank"+std::to_string(rank));
      output();
    }
    mainLoop();
    if (iTime % reportInterval == 0) {
      if((iTime/reportInterval) % rankSize == rank) std::cout << "rank = " << rank << ", \ti = " << iTime << std::endl;
    }
  }
  output.SetFileName(std::to_string(maxStep)+"rank"+std::to_string(rank));
  output();
  std::cout << "rank = " << rank << " finish." << std::endl;
}

template<typename LC, typename CLS, typename EQ>
class CollisionAndStream {
  const LC &lc;
  const CLS &kernel;
  const EQ &equilibrium;
  nnwlbm::LinkedMaskData<double, LC::GetDim()> &pdfData;
  nnwlbm::LinkedMaskData<double, LC::GetDim()> pdfData2;
  nnwlbm::LinkedMaskData<double, LC::GetDim()> &macroData;
  nnwlbm::LinkedMaskData<nnwlbm::SingleNodeBc<LC>*, LC::GetDim(), nnwlbm::BcIdentifier > &bcMapSingle;
  const nnwlbm::UniformGrid<LC::GetDim()> &gridC;
  nnwlbm::LinkedMaskData<double, LC::GetDim()> *ptrCurrent;
  nnwlbm::LinkedMaskData<double, LC::GetDim()> *ptrNext;
  nnwlbm::Range<LC::GetDim()> loopRange;
  nnwlbm::RangeIndex<LC::GetDim()> rIndex;
  PdfCommunicator<LC, nnwlbm::LinkedMaskData > &pdfComm;
public:
  CollisionAndStream(
    const LC &lcIn,
    const CLS &kernelIn,
    const EQ &equilibriumIn,
    nnwlbm::LinkedMaskData<double, LC::GetDim()> &pdfDataIn,
    nnwlbm::LinkedMaskData<double, LC::GetDim()> &macroDataIn,
    nnwlbm::LinkedMaskData<nnwlbm::SingleNodeBc<LC>*, LC::GetDim(), nnwlbm::BcIdentifier > &bcMapSingleIn,
    const nnwlbm::UniformGrid<LC::GetDim()> &gridCIn,
    PdfCommunicator<LC, nnwlbm::LinkedMaskData> &pdfCommIn
  ) : lc(lcIn), kernel(kernelIn), equilibrium(equilibriumIn), pdfData(pdfDataIn), macroData(macroDataIn),
    bcMapSingle(bcMapSingleIn), gridC(gridCIn), pdfData2(pdfDataIn), ptrCurrent(&pdfData), ptrNext(&pdfData2),
    loopRange(pdfData.GetRange()), rIndex(pdfData.GetNodeIndex()), pdfComm(pdfCommIn)
  { 
    loopRange.AddBuffer(-1);
  }
  bool operator() (){
    //return true;
    int loc[LC::GetDim()];
    int locN[LC::GetDim()];
    // collision
    for (int carry = loopRange.IniCarry(loc); carry == 0; carry = loopRange.AddOne(loc))
    {
      size_t nodeId = rIndex.CalNodeID(loc);
      if(macroData.IsValid(nodeId)) {
        kernel.collide(ptrCurrent->GetVar(nodeId), macroData.GetVar(nodeId));
      }
    }
    pdfComm.SendAndRecv();
    pdfComm.wait();
    pdfComm.update();
    for (int carry = loopRange.IniCarry(loc); carry == 0; carry = loopRange.AddOne(loc))
    {
      size_t nodeId = rIndex.CalNodeID(loc);
      if( !pdfData.IsValid(nodeId) ) continue;
      size_t nodeIdN = std::string::npos;
      double* pdfNext = ptrNext->GetVar(nodeId);
      switch(macroData.GetColor(nodeId)) {
      case nnwlbm::geo::inside :
        for (size_t i = 0; i < lc.size(); i++)
        {
          const int* displacement = lc.direction(i);
          for (size_t iDim = 0; iDim < lc.GetDim(); iDim++) {
            locN[iDim] = loc[iDim] - displacement[iDim];
          }
          nodeIdN = rIndex.CalNodeID(locN);
          pdfNext[i] = ptrCurrent->GetVar(nodeIdN)[i];
        }
        break;
      case nnwlbm::geo::boundary :
        for (size_t i = 0; i < lc.size(); i++)
        {
          const int* displacement = lc.direction(i);
          for (size_t iDim = 0; iDim < lc.GetDim(); iDim++) {
            locN[iDim] = loc[iDim] - displacement[iDim];
          }
          nodeIdN = rIndex.CalNodeID(locN);
          if(macroData.IsValid(nodeIdN)) {
            pdfNext[i] = ptrCurrent->GetVar(nodeIdN)[i];
          }
          else {
            auto bcNode = bcMapSingle.GetVar(nodeId);
            (*(bcNode[i]))(pdfNext, ptrCurrent->GetVar(nodeId), i, macroData.GetVar(nodeId));
          }
        }
        break;
      default:
        continue;
      }
    }
    std::swap(ptrCurrent, ptrNext);
    pdfComm.SetPdfDataWithIdenticalMask(ptrCurrent);
    return true;
  }
};

template<typename LC, typename CLS, typename EQ>
class CollisionAndStream1 {
  const LC& lc;
  const CLS& kernel;
  const EQ& equilibrium;
  nnwlbm::LinkedMaskData<double, LC::GetDim()> &pdfData;
  nnwlbm::LinkedMaskData<double, LC::GetDim()> pdfData2;
  nnwlbm::LinkedMaskData<double, LC::GetDim()> &macroData;
  nnwlbm::LinkedMaskData<std::shared_ptr<nnwlbm::SimpleBC>, LC::GetDim(), nnwlbm::BcIdentifier > &bcMapSingle;
  const nnwlbm::UniformGrid<LC::GetDim()> &gridC;
  nnwlbm::LinkedMaskData<double, LC::GetDim()> *ptrCurrent;
  nnwlbm::LinkedMaskData<double, LC::GetDim()> *ptrNext;
  nnwlbm::Range<LC::GetDim()> loopRange;
  nnwlbm::RangeIndex<LC::GetDim()> rIndex;
  PdfCommunicator<LC, nnwlbm::LinkedMaskData > &pdfComm;
public:
  CollisionAndStream1(
    const LC &lcIn,
    const CLS &kernelIn,
    const EQ &equilibriumIn,
    nnwlbm::LinkedMaskData<double, LC::GetDim()> &pdfDataIn,
    nnwlbm::LinkedMaskData<double, LC::GetDim()> &macroDataIn,
    nnwlbm::LinkedMaskData<std::shared_ptr<nnwlbm::SimpleBC>, LC::GetDim(), nnwlbm::BcIdentifier > &bcMapSingleIn,
    const nnwlbm::UniformGrid<LC::GetDim()> &gridCIn,
    PdfCommunicator<LC, nnwlbm::LinkedMaskData> &pdfCommIn
  ) : lc(lcIn), kernel(kernelIn), equilibrium(equilibriumIn), pdfData(pdfDataIn), macroData(macroDataIn),
    bcMapSingle(bcMapSingleIn), gridC(gridCIn), pdfData2(pdfDataIn), ptrCurrent(&pdfData), ptrNext(&pdfData2),
    loopRange(pdfData.GetRange()), rIndex(pdfData.GetNodeIndex()), pdfComm(pdfCommIn)
  {
    loopRange.AddBuffer(-1);
  }
  bool operator() () {
    //return true;
    int loc[LC::GetDim()];
    int locN[LC::GetDim()];
    // collision
    for (int carry = loopRange.IniCarry(loc); carry == 0; carry = loopRange.AddOne(loc))
    {
      size_t nodeId = rIndex.CalNodeID(loc);
      if (macroData.IsValid(nodeId)) {
        kernel.collide(ptrCurrent->GetVar(nodeId), macroData.GetVar(nodeId));
      }
    }
    pdfComm.SendAndRecv();
    pdfComm.wait();
    pdfComm.update();
    for (int carry = loopRange.IniCarry(loc); carry == 0; carry = loopRange.AddOne(loc))
    {
      size_t nodeId = rIndex.CalNodeID(loc);
      if (!pdfData.IsValid(nodeId)) continue;
      size_t nodeIdN = std::string::npos;
      double* pdfNext = ptrNext->GetVar(nodeId);
      switch (macroData.GetColor(nodeId)) {
      case nnwlbm::geo::inside:
        for (size_t i = 0; i < lc.size(); i++)
        {
          const int* displacement = lc.direction(i);
          for (size_t iDim = 0; iDim < lc.GetDim(); iDim++) {
            locN[iDim] = loc[iDim] - displacement[iDim];
          }
          nodeIdN = rIndex.CalNodeID(locN);
          pdfNext[i] = ptrCurrent->GetVar(nodeIdN)[i];
        }
        break;
      case nnwlbm::geo::boundary:
        for (size_t i = 0; i < lc.size(); i++)
        {
          const int* displacement = lc.direction(i);
          for (size_t iDim = 0; iDim < lc.GetDim(); iDim++) {
            locN[iDim] = loc[iDim] - displacement[iDim];
          }
          nodeIdN = rIndex.CalNodeID(locN);
          if (macroData.IsValid(nodeIdN)) {
            pdfNext[i] = ptrCurrent->GetVar(nodeIdN)[i];
          }
          else {
            auto bcNode = bcMapSingle.GetVar(nodeId);
            (*(bcNode[i]))(pdfNext, ptrCurrent->GetVar(nodeId), i);
          }
        }
        break;
      default:
        continue;
      }
    }
    std::swap(ptrCurrent, ptrNext);
    pdfComm.SetPdfDataWithIdenticalMask(ptrCurrent);
    return true;
  }
};
