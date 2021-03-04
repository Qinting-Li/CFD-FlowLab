//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//          PPPPP  H   H  EEEEE  N    N  GGGGG  L      EEEEE  III         +
//          P   P  H   H  E      NN   N  G      L      E       I          +
//          PPPPP  HHHHH  EEEEE  N N  N  G  GG  L      EEEEE   I          +
//          P      H   H  E      N  N N  G   G  L      E       I          +
//          P      H   H  EEEEE  N    N  GGGGG  LLLLL  EEEEE  III         +
//------------------------------------------------------------------------+
//          Platform for Hybrid Engineering Simulation of Flows           +
//          China Aerodynamics Research and Development Center            +
//                     (C) Copyright, Since 2010                          +
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      TurbSolverStructFD.h
//! @brief     turbulence solver for struct grid.
//! @author    Ma Yankai, Min Yaobing.

#pragma once
#include "TurbSolverStruct.h"

namespace PHSPACE
{
class TurbSolverStrFD : public TurbSolverStr
{
public:
    TurbSolverStrFD();
    ~TurbSolverStrFD();

private:
    void Boundary(Grid *gridIn);
    void ComputeViscousCoeff(Grid *grid);
    void ComputeGradientCellCenter(Grid *gridIn);
    void GetDependentVariablesforStructHighOrder(Grid *gridIn);
    void ObtainViscousCoefWithGhost(Grid *gridIn);
    void ObtainqlqrLowOrder(Grid *gridIn);
    void GetFaceVariable(Grid *gridIn, FieldProxy *qlProxy, FieldProxy *qrProxy, int nsurf);
    RDouble Limiter_ThirdSmooth(const RDouble &x, const RDouble &y);
    void GetFaceVariableCorrectionAtPhysicalBoundary(Grid *gridIn, FieldProxy *qlProxy, FieldProxy *qrProxy, int nsurf);
    void ObtainGradientCellCenter_Laminar(Grid *gridIn);
    void ObtainBoundaryValue(Grid *gridIn);
    void ObtainViscosity(Grid *gridIn);
    void ObtainGradientCellCenter(Grid *gridIn);
    void ObtainCrossing(Grid *gridIn);
    void ObtainBlending(Grid *gridIn);
    void InFlowBC(Grid *grid, BCRegion *BCRegion);
    void OutFlowBC(Grid *gridIn, BCRegion *BCRegion);
    void VisWall(Grid *gridIn, BCRegion *BCRegion);
    void VisWallWithWallFunctionStandard(Grid* gridIn, BCRegion* BCRegion);
    void VisWallWithWallFunctionPAB3D(Grid* gridIn, BCRegion* BCRegion);
    void SymmetryBC(Grid *gridIn, BCRegion *BCRegion);
    void FarFieldBC(Grid *gridIn, BCRegion *BCRegion);
    void CornerPoint(Grid *gridIn);
    void UpdateFlowField(Grid *gridIn, FieldProxy *qProxy, FieldProxy *dqProxy);
    void Spectrum(Grid *gridIn);
    void SpectrumRadiusOfOneEquation(Grid *gridIn);
    void SpectrumRadiusOfTwoEquation(Grid *gridIn);
    void LUSGSInitializationStructHighOrder(Grid *gridIn, FieldProxy *dqProxy, FieldProxy *LUplusDQ);
    void CompressDQ(DataContainer *&dataContainer, FieldProxy *fieldProxy, Grid *gridIn, const int &zoneIndex, const int &neighborZoneIndex, const int &nEquation);
    void DecompressDQ(DataContainer *&dataContainer, FieldProxy *fieldProxy, Grid *gridIn, const int &zoneIndex, const int &neighborZoneIndex, const int &nEquation);
    void SolveLUSGSForward(Grid *gridIn, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal);
    void SolveLUSGSBackward(Grid *gridIn, FieldProxy * dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal);
};

}
