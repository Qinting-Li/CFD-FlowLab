#include "Geo_Interface.h"
#include "Geo_Interpoint.h"
#include "Geo_UnstructGrid.h"
#include "PHMpi.h"
using namespace std;

namespace PHSPACE
{
vector< vector< Grid *> * > ggrids;
ZoneConnectivity *zoneConnectivity;
ZoneConnectivityForPoint *zoneConnectivityForPoint;
RDouble globalGridMinBox[3], globalGridMaxBox[3];

Grid *GetGrid(int izone, int level)
{
    if (!ggrids[izone])
    {
        return 0;
    }

    return (*ggrids[izone])[level];
}

void DeleteGGrids()
{
    ggrids.clear();
}

InterfaceInfo * GetInterfaceInfo(int izone, int level)
{
    Grid *grid = GetGrid(izone, level);
    if (!grid) return 0;
    return grid->GetInterfaceInfo();
}

InterpointInformation *GetInterpointInfo(int izone, int level)
{
    Grid *grid = GetGrid(izone, level);
    if (! grid)
    {
        return 0;
    }
    //! Temporarily debug to obtain the interpoint information of unstructured grid.
    return (static_cast<UnstructGrid*>(grid))->GetInterpointInfo();
}

OversetInfoProxy * GetOversetInfoProxy(int izone, int level)
{
    Grid *grid = GetGrid(izone, level);
    if (!grid) return 0;
    return grid->GetOversetInfoProxy();
}

uint_t GetNumberOfMultiGrid(int izone)
{
    return ggrids[izone]->size();
}

uint_t GetNumberOfMultiGrid()
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    for (int izone = 0; izone < nZones; ++ izone)
    {
        if (ggrids[izone])
        {
            return ggrids[izone]->size();
        }
    }
    return 0;
}

void AddGridToGlobal(int izone, vector < Grid   * > *grids)
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    if (ggrids.size() == 0)
    {
        ggrids.resize(nZones,0);
    }
    ggrids[izone] = grids;
}

void GetGlobalMinMaxBox(RDouble * globalMinBox, RDouble * globalMaxBox)
{
    for (int iDim = 0; iDim < 3; ++ iDim)
    {
        globalMinBox[iDim] = globalGridMinBox[iDim];
        globalMaxBox[iDim] = globalGridMaxBox[iDim];
    }
}

void SetGlobalMinMaxBox(RDouble * globalMinBox, RDouble * globalMaxBox)
{
    for (int iDim = 0; iDim < 3; ++ iDim)
    {
        globalGridMinBox[iDim] = globalMinBox[iDim];
        globalGridMaxBox[iDim] = globalMaxBox[iDim];
    }
}

ZoneConnectivityForPoint::ZoneConnectivityForPoint()
{
    numberOfZonesForPoint = 0;
    zoneNeighborForPoint  = 0;
};

ZoneConnectivityForPoint::~ZoneConnectivityForPoint()
{
    for (int izone = 0; izone < numberOfZonesForPoint; ++ izone)
    {
        delete zoneNeighborForPoint[izone];
    }
    delete [] zoneNeighborForPoint;
}

int  ZoneConnectivityForPoint::GetNumberOfZonesForPoint() const
{ 
    return numberOfZonesForPoint;
};

void ZoneConnectivityForPoint::SetNumberOfZonesForPoint(int numberOfZonesForPoint)
{
    this->numberOfZonesForPoint = numberOfZonesForPoint;
};

void ZoneConnectivityForPoint::CreateAllZones(int numberOfZonesForPoint)
{
    SetNumberOfZonesForPoint(numberOfZonesForPoint);
    zoneNeighborForPoint = new ZoneNeighbor *[numberOfZonesForPoint];
    for (int izone = 0; izone < numberOfZonesForPoint; ++ izone)
    {
        zoneNeighborForPoint[izone] = new ZoneNeighbor();
    }
}

ZoneNeighbor * ZoneConnectivityForPoint::GetZoneNeighborForPoint(int izone) const 
{
    return zoneNeighborForPoint[izone];
}

GridID::GridID(int index, int local_index)
{
    this->index       = index;
    this->local_index = local_index;
}

GridID::~GridID()
{
}

OversetInfo::OversetInfo()
{
    cellIndexForTransfer = 0;
    Data_ParamFieldSuite *dataStorageForTransfer = new Data_ParamFieldSuite();
}

OversetInfo::~OversetInfo()
{
    delete dataStorageForTransfer;
    delete [] cellIndexForTransfer;
}

void OversetInfo::AllocateData()
{
    int numberOfCellForTransfer = this->GetNumberOfCellForTransfer();
    int * cellIndexForTransfer = new int [ numberOfCellForTransfer ];
    this->SetCellIndexForTransfer(cellIndexForTransfer);
}

void OversetInfo::Encode(DataContainer *cdata)
{
    ;
}

void OversetInfo::Decode(DataContainer *cdata)
{
    int zoneIndexOfNeighbor = 0;
    cdata->read(reinterpret_cast<char *>(&zoneIndexOfNeighbor),sizeof(int));
    SetZoneIndexOfNeighbor(zoneIndexOfNeighbor);

    int numberOfCellForTransfer;
    cdata->read(reinterpret_cast<char *>(&numberOfCellForTransfer),sizeof(int));

    int * cellIndexForTransfer = new int [numberOfCellForTransfer];
    cdata->read(reinterpret_cast<char *>(cellIndexForTransfer), numberOfCellForTransfer * sizeof(int));
    SetCellIndexForTransfer(cellIndexForTransfer);
}

CompositeOversetInfo::CompositeOversetInfo()
{
}

CompositeOversetInfo::~CompositeOversetInfo()
{
}


void CompositeOversetInfo::Encode(DataContainer *cdata)
{
    int numberOfZone = this->GetNumberOfZone();
    cdata->write(reinterpret_cast<char *>(&numberOfZone),sizeof(int));

    int numberOfCell = this->GetNumberOfCell();
    cdata->write(reinterpret_cast<char *>(&numberOfCell),sizeof(int));

    for (int izone = 0; izone < numberOfZone; ++ izone)
    {
        OversetInfo *oversetInfo = GetOversetInfo(izone);
        oversetInfo->Encode(cdata);
    }
}

void CompositeOversetInfo::Decode(DataContainer *cdata)
{
    int numberOfZone;
    cdata->read(reinterpret_cast<char *>(&numberOfZone),sizeof(int));

    SetNumberOfZone(numberOfZone);

    int numberOfCell;
    cdata->read(reinterpret_cast<char *>(&numberOfCell),sizeof(int));
    SetNumberOfCell(numberOfCell);

    for (int izone = 0; izone < numberOfZone; ++ izone)
    {
        OversetInfo *oversetInfo = new OversetInfo();
        AddOversetInfo(oversetInfo);
        oversetInfo->Decode(cdata);
    }
}

OversetInfoProxy::OversetInfoProxy()
{
    sendIndexProxy = 0;
    recvIndexProxy = 0;
}

OversetInfoProxy::~OversetInfoProxy()
{
}

void OversetInfoProxy::AllocateNeighborInfo()
{
}

void OversetInfoProxy::Encode(DataContainer *cdata)
{
    sendInfo->Encode(cdata);
    recvInfo->Encode(cdata);
}

void OversetInfoProxy::Decode(DataContainer *cdata)
{
    sendInfo->Decode(cdata);
    recvInfo->Decode(cdata);
}

int OversetInfoProxy::GetNIFaceOfNeighbor(int ineighbor) const
{
    return -1;
}

void OversetInfoProxy::InitReceivingInformation(vector<int> &cellIndexRef1, vector<int> &zoneIndexRef2, vector<int> &cellIndexRef2)
{
    this->cellIndexRef1 = cellIndexRef1;
    this->zoneIndexRef2 = zoneIndexRef2;
    this->cellIndexRef2 = cellIndexRef2;

    set<int> neighborZoneIndexSet;
    GetNeighborZoneIndexSet(neighborZoneIndexSet);

    IndexProxy * recvIndexProxy      = GetRecvIndexProxy();
    IndexProxy * temporaryIndexProxy = GetTemporaryIndexProxy();

    for (set<int>::iterator iter = neighborZoneIndexSet.begin(); iter != neighborZoneIndexSet.end(); ++ iter)
    {
        //! jzone is the neighbor,which is the interpolated zone overlapped with izone.
        int jzone = *iter;

        vector<int> cellIndexI, cellIndexJ;

        GetPartOfArrayByCondition(zoneIndexRef2, jzone, cellIndexRef1, cellIndexI);
        GetPartOfArrayByCondition(zoneIndexRef2, jzone, cellIndexRef2, cellIndexJ);

        IndexSlice * indexSlice = new IndexSlice();

        indexSlice->SetTag(jzone);
        indexSlice->SetIndex(cellIndexI);
        recvIndexProxy->AddIndexSlice(indexSlice);

        IndexSlice * temporaryIndexSlice = new IndexSlice();

        temporaryIndexSlice->SetTag(jzone);
        temporaryIndexSlice->SetIndex(cellIndexJ);
        temporaryIndexProxy->AddIndexSlice(temporaryIndexSlice);
    }
}

void OversetInfoProxy::PrepareSendingInformation(Geo_OversetGrid **oversetGrids)
{
    IndexProxy * temporaryIndexProxy = GetTemporaryIndexProxy();

    uint_t numberOfNeighbor = temporaryIndexProxy->size();
    int zoneIndex = this->GetGrid()->GetZoneID();
    
    for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
    {
        IndexSlice * temporaryIndexSlice = temporaryIndexProxy->GetIndexSlice(iNeighbor);
        int jzone = temporaryIndexSlice->GetTag();
        OversetInfoProxy * jOversetInfoProxy = GetOversetInfoProxy(jzone);
        IndexProxy *sendIndexProxy = jOversetInfoProxy->GetSendIndexProxy();

        IndexSlice * sendIndexSlice = new IndexSlice();
        sendIndexSlice->SetTag(zoneIndex);
        sendIndexSlice->SetIndex(temporaryIndexSlice->GetIndex());
        sendIndexProxy->AddIndexSlice(sendIndexSlice);
    }
}

void OversetInfoProxy::PostProcessRecv()
{
    IndexProxy * recvIndexProxy = GetRecvIndexProxy();
    int numberOfCellForRecv = 0;

    uint_t numberOfNeighbor = recvIndexProxy->size();
    
    for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
    {
        int iLength = recvIndexProxy->GetLength(iNeighbor);
        recvIndexProxy->SetIndexMap(iNeighbor, iLength);
        numberOfCellForRecv += iLength;
    }

    this->SetNumberOfCellForRecv(numberOfCellForRecv);
}

void OversetInfoProxy::PostProcessSend()
{
    IndexProxy * sendIndexProxy = GetSendIndexProxy();

    int numberOfCellForSend = 0;
    uint_t numberOfNeighbor = sendIndexProxy->size();
    
    for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
    {
        int iLength = sendIndexProxy->GetLength(iNeighbor);
        sendIndexProxy->SetIndexMap(iNeighbor, iLength);
        numberOfCellForSend += iLength;
    }

    this->SetNumberOfCellForSend(numberOfCellForSend);
}

void OversetInfoProxy::PostProcess(Geo_OversetGrid **oversetGrids)
{
    PostProcessRecv();
    PostProcessSend();
}

void OversetInfoProxy::GetNeighborZoneIndexSet(set<int> &neighborZoneIndexSet)
{
    Array2Set(zoneIndexRef2, neighborZoneIndexSet);
}

void OversetInfoProxy::Process()
{
    set<int> neighborZoneIndexSet;
    GetNeighborZoneIndexSet(neighborZoneIndexSet);

    IndexProxy * recvIndexProxy      = GetRecvIndexProxy();
    IndexProxy * temporaryIndexProxy = GetTemporaryIndexProxy();

    vector<int> & cellIndexRef1 = GetCellIndexRef1();
    vector<int> & cellIndexRef2 = GetCellIndexRef2();
    vector<int> & zoneIndexRef2 = GetZoneIndexRef2();

    for (set<int>::iterator iter = neighborZoneIndexSet.begin(); iter != neighborZoneIndexSet.end(); ++ iter)
    {
        //! jzone is the neighbor,which is the interpolated zone overlapped with izone.
        int jzone = *iter;

        vector<int> cellIndexI, cellIndexJ;

        GetPartOfArrayByCondition(zoneIndexRef2, jzone, cellIndexRef1, cellIndexI);
        GetPartOfArrayByCondition(zoneIndexRef2, jzone, cellIndexRef2, cellIndexJ);

        IndexSlice * indexSlice = new IndexSlice();

        indexSlice->SetTag(jzone);
        indexSlice->SetIndex(cellIndexI);
        recvIndexProxy->AddIndexSlice(indexSlice);

        IndexSlice * temporaryIndexSlice = new IndexSlice();

        temporaryIndexSlice->SetTag(jzone);
        temporaryIndexSlice->SetIndex(cellIndexJ);
        temporaryIndexProxy->AddIndexSlice(temporaryIndexSlice);
    }
}


}
