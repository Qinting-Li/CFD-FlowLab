#pragma once

#include "BoundaryCondition.h"
#include "Collision.h"
#include "GenBcMap.h"
#include "GenColorMap.h"
#include "InitField.h"

namespace nnwlbm {
  template< typename MD, typename Grid >
  class Output
  {
    const MD &maskData;
    const Grid &grid;
    std::string fileName;
  public:
    Output(const MD &md, const Grid &gridIn) : maskData(md), grid(gridIn) {}
    void SetFileName(const std::string &fname = "") { fileName = fname; }
    bool operator()() {
      return nnwlbm::OutputCellCenteredcellData(maskData, grid, fileName);
    }
  };

  template<typename Output, typename MainLoop>
  void time_marching(
    Output &output,
    MainLoop &mainLoop,
    const size_t &maxStep,
    const size_t &outputInterval,
    const size_t &reportInterval
  ) {
    for(int iTime = 0; iTime < maxStep; iTime++) {
      if (iTime % outputInterval == 0) {
        output.SetFileName(std::to_string(iTime));
        output();
      }
      mainLoop();
      if (iTime % reportInterval == 0) std::cout << "i = " << iTime << std::endl;
    }
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
  public:
    CollisionAndStream(
      const LC &lcIn,
      const CLS &kernelIn,
      const EQ &equilibriumIn,
      nnwlbm::LinkedMaskData<double, LC::GetDim()> &pdfDataIn,
      nnwlbm::LinkedMaskData<double, LC::GetDim()> &macroDataIn,
      nnwlbm::LinkedMaskData<nnwlbm::SingleNodeBc<LC>*, LC::GetDim(), nnwlbm::BcIdentifier > &bcMapSingleIn,
      const nnwlbm::UniformGrid<LC::GetDim()> &gridCIn
    ) : lc(lcIn), kernel(kernelIn), equilibrium(equilibriumIn), pdfData(pdfDataIn), macroData(macroDataIn),
      bcMapSingle(bcMapSingleIn), gridC(gridCIn), pdfData2(pdfDataIn), ptrCurrent(&pdfData), ptrNext(&pdfData2),
      loopRange(pdfData.GetRange()), rIndex(pdfData.GetNodeIndex())
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
      for (int carry = loopRange.IniCarry(loc); carry == 0; carry = loopRange.AddOne(loc))
      {
        size_t nodeId = rIndex.CalNodeID(loc);
        if( !pdfData.IsValid(nodeId) ) continue;
        size_t nodeIdN = std::string::npos;
        double *pdfNext = ptrNext->GetVar(nodeId);
        switch(macroData.GetColor(nodeId)) {
        case nnwlbm::geo::inside :
          for (size_t i = 0; i < lc.size(); i++)
          {
            const int *displacement = lc.direction(i);
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
            const int *displacement = lc.direction(i);
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
      return true;
    }
  };
}


