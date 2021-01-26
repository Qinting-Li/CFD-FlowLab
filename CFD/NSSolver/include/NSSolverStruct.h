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
//! @file      NSSolverStruct.h
//! @brief     NS solver for struct grid.
//! @author    He Xin, Bell, He Kun, He Xianyao, Wan Yunbo, Liu Jian, Xu Qingxin, Zhang Zipei, Li peng, Ma Yankai, 
//!            Min Yaobing, Meng Liyuan, Xu Gang, Guo Yongyan, Zhang Yang.

#pragma once
#include "NSSolver.h"
using namespace std;

namespace PHSPACE
{
class GeomProxy;
class BCRegion;
class Param_NSSolverStruct;
class GradientCellCenter;
class Gradient;

class StructIter
{
public:
    int i0,j0,k0;
    int i1,j1,k1;
    int il1,jl1,kl1;
    int ist, ied, jst, jed, kst, ked;
    int nsurf;
    int nsize;
};

class NSSolverStruct : public NSSolver
{
private:
    GradientCellCenter *gradientCellCenterUVWT[5];
    int systemgridtype;
public:

    NSSolverStruct();
    ~NSSolverStruct();
private:
    void ReadParameter();
    void AllocateGlobalVar  (Grid *gridIn);
    void DeAllocateGlobalVar(Grid *gridIn);
    
    //! Judge if the flow field variables file already exists in NS solver.
    bool JudgeIfRestart();
    bool JudgeIfReadAverage();

    void InitFlowAsRestart();
    void InitGrad(Grid *gridIn);
    void InitCGrid(Grid *fineGridIn, Grid *coarseGridIn);
    void ZeroResiduals(Grid *grid);

    void CreateStandardADT(Grid *grid, DataStruct_AdtTree<int,RDouble> *&adtTree, RDouble &tol);

    void CompNodeVar(Grid *gridIn, FieldProxy *qnProxy, FieldProxy *qProxy);
    void InterfaceValueFix(Grid *gridIn, FieldProxy *qnProxy, FieldProxy *qProxy);
    void WallValueFix(Grid *gridIn, FieldProxy *qnProxy, FieldProxy *qProxy);
    
    void CompNodeVar(Grid *gridIn, RDouble4D &qn, int m, RDouble4D &q, int n);
    void CompNodeVarForVelocity(Grid *gridIn, RDouble4D &qn, int m, RDouble4D &q, int n);
    void CompNodeVar(Grid *gridIn, RDouble4D &qn, int m, RDouble3D &q);
    
    void LoadQ(Grid *grid, FieldProxy *qProxy);
    void LoadResiduals(Grid *gridIn, FieldProxy *rightHandSideProxy);
    void ComputeMinTimeStep(Grid *grid, RDouble &minDt, RDouble &maxDt);
    void ReduceMaxTimeStep(Grid *grid, RDouble globalMinDt);
    void LocalTimeStep(Grid *gridIn);
    void GlobalTimeStep(Grid *gridIn);
    void LocalGlobalTimeStep(Grid *gridIn);
    void SpectrumRadius(Grid *grid);
    void SpectrumRadiusInviscid(Grid *gridIn);
    void PreconSpectrumRadiusInviscid(Grid *gridIn);
    void SpectrumRadiusViscous(Grid *gridIn);

    void Spectrum(Grid *gridIn);

    void Boundary(Grid *grid);
    void CalPressureFactor(Grid *gridIn);
    void InviscidFlux(Grid *gridIn);
    void ZeroResidualOfSpecialCells(Grid *gridIn);    //! Added By Guo Yongheng 20160825.
    void ViscousFlux(Grid *gridIn);
    void DualTimeSource(Grid *gridIn);

    //! Compute the chemical source terms.
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    void ChemicalSource(Grid *gridIn);

    void UpdateUnsteadyFlow(Grid *gridIn);
    RDouble UnsteadyConvergence(Grid *gridIn);
    
    void ComputeViscousCoeff(Grid *gridIn);
    void ComputeViscousCoefficientWithoutChemical(Grid *gridIn);

    void ComputeGamaAndTemperature(Grid *gridIn);
    void RungeKuttaResidual(Grid *gridIn, FieldProxy *dqProxy, RDouble coef);
    
    void SolutionFix(Grid *gridIn, RDouble *prim, int i, int j, int k);
    void UpdateFlowField(Grid *gridIn, FieldProxy *qProxy, FieldProxy *dqProxy);
    
    void CopyInviscidfluxProxytoGlobalfacefluxOnlyforMixGrid(Grid *gridIn, FieldProxy *fluxProxy, int iSurface);
    void CopyViscousfluxProxytoGlobalfacefluxOnlyforMixGrid(Grid *gridIn, FieldProxy *fluxProxy, int iSurface);
    void LoadGlobalfacefluxtoResOnlyforMixGrid(Grid *gridIn);

    void LoadFlux(Grid *gridIn, FieldProxy *fluxProxy, int iSurface);
    void CompVisFlux(Grid *gridIn, FieldProxy *fluxProxy, int iSurface);
    //! To correct the viscous fluxes of solid boundaries under the non-equilibrium condition.
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    //! @param[out]: fluxProxy is an array that stores the inviscid fluxes of grid faces.
    //! @param[in ]: iSurface denotes the direction index of the grid faces, where 1 indicates the I-direction, 2 for J-direction, and 3 is K-direction.
    void CorrectViscousFlux(Grid *gridIn, FieldProxy *fluxProxy, int iSurface);

    void CompVisFluxnew(Grid *gridIn, Gradient *gradientUVWT, FieldProxy *fluxProxy, int iSurface);
    void CompVisFluxnew(Grid *gridIn, FieldProxy *fluxProxy, int iSurface);
    void ComputeInviscidFlux(Grid *gridIn, FieldProxy *qlProxy, FieldProxy *qrProxy, int iSurface);
    void InviscidFluxWrap(Grid *gridIn, FieldProxy *qlProxy, FieldProxy *qrProxy, FieldProxy *fluxProxy, INVScheme inviscidScheme, int iSurface);
    //! To correct the inviscid fluxes of solid boundaries under the non-equilibrium condition.
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    //! @param[out]: fluxProxy is an array that stores the inviscid fluxes of grid faces.
    //! @param[in ]: iSurface denotes the direction index of the grid faces, where 1 indicates the I-direction, 2 for J-direction, and 3 is K-direction.
    void CorrectInviscidFlux(Grid *gridIn, FieldProxy *fluxProxy, int iSurface);

    //! To correct the primitive variables of faces nearby the wall.
    //! @param[in ]: gridIn denotes the current computational regions of the grids.
    //! @param[out]: inviscidFaceProxy is an array that stores the primitive variables of grid faces.
    //! @param[in ]: iSurface denotes the direction index of the grid faces, where 1 indicates the I-direction, 2 for J-direction, and 3 is K-direction.
    void CorrectFaceVar(Grid *gridIn, FieldProxy *qlProxy, FieldProxy *qrProxy, int iSurface);
    void CorrectChemicalFaceVar(Grid *gridIn, FieldProxy *qlProxy, FieldProxy *qrProxy, int iSurface);

    void GetInvFaceVar(Grid *gridIn, FieldProxy *qlProxy, FieldProxy *qrProxy, int iSurface);
    
    //! Compute gradient of primitive variables & temperature for viscous flux computation.
    void ComputeGradient(Grid *grid);
    void ComputeGradientCellCenter(Grid *grid);
    void ComputeGradientCellCenterOfruvwptOnlyForMixGrid(Grid *grid);
    void ReconGrad(Grid *gridIn, int iSurface);

    void RegisterCFDSolverInterfaceField();
    void CompressSpecificArrayToInterface(DataContainer *&dataContainer, const string &fieldName, Grid *grid, const int &neighborZoneIndex, const int &neqn);
    void DecompressArrayFromInterface(DataContainer *&dataContainer, const string &fieldName, Grid *grid, const int &neighborZoneIndex, const int &neqn);

    void CompressGradientToInterfaceOnlyForMixGrid(DataContainer *&dataContainer, Grid *grid, const int &neighborZoneIndex);
    void DecompressGradientToInterfaceOnlyForMixGrid(DataContainer *&dataContainer, Grid *grid, const int &neighborZoneIndex);

    void CompressFacefluxToInterfaceOnlyForMixGrid(DataContainer *&dataContainer, Grid *grid, const int &neighborZoneIndex);
    void DecompressFacefluxToInterfaceOnlyForMixGrid(DataContainer *&dataContainer, Grid *grid, const int &neighborZoneIndex);

    void CompressDQ(DataContainer *&dataContainer, FieldProxy *fieldProxy, Grid *grid, const int &zoneIndex, const int &neighborZoneIndex, const int &neqn);
    void DecompressDQ(DataContainer *&dataContainer, FieldProxy *fieldProxy, Grid *grid, const int &zoneIndex, const int &neighborZoneIndex, const int &neqn);

    void RestrictAllQ(Grid *fineGridIn, Grid *coarseGridIn);
    void RestrictDefect(Grid *fineGridIn, Grid *coarseGridIn);
    void PutCorrectionBack(Grid *grid, FieldProxy *qProxy);
    void CorrectFineGrid(Grid *fineGridIn, Grid *coarseGridIn);
    void InterpolatFineGrid(Grid *fineGridIn, Grid *coarseGridIn);
    void prol   (Grid *fineGridIn, Grid *coarseGridIn, FieldProxy *workProxy);
    void prolhic_eric(Grid *fineGridIn, Grid *coarseGridIn, FieldProxy *workProxy);
    void PositiveUpdate(Grid *grid, FieldProxy *qProxy, FieldProxy *dqProxy);

    void SolveLUSGSForward(Grid *grid, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal);
    void SolveLUSGSBackward(Grid *grid, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal);


    void SetGhostDQLUSGS(Grid *gridIn, RDouble4D &dq);
    void OutFlowBC3D (Grid *gridIn, BCRegion *BCRegion);

    void SymmetryBC3D(Grid *gridIn, BCRegion *BCRegion);
    void VisWallADI3D(Grid *gridIn, BCRegion *BCRegion);
    void VisWallISO3D(Grid *gridIn, BCRegion *BCRegion);

    void InFlowBC3D  (Grid *gridIn, BCRegion *BCRegion);
    void FarFieldRiemannInvariants(Grid *gridIn, BCRegion *BCRegion);
    void PreconFarfieldBC(Grid *gridIn, BCRegion *BCRegion);
    void PressureInletBCRiemann(Grid *gridIn, BCRegion *BCRegion);
    void PressureOutletBC(Grid *gridIn, BCRegion *BCRegion);
    void CornerPoint(Grid *gridIn);
    
    void FillField(Grid *gridIn, FieldProxy *fieldProxy, RDouble value);
    void FillField(Grid *gridIn, FieldProxy *field1Proxy, FieldProxy *field2Proxy);
    
    FieldProxy * GetFieldProxy(Grid *gridIn, const string &field_name);
    FieldProxy * CreateFieldProxy(Grid *gridIn);
    FieldProxy * GetResidualProxy(Grid *gridIn);

    void ReconGradVolWeightedwithCorrection(Grid *gridIn, int nSurface);
    void GetGradientAtFace_CorrectionAtPhysicalBoundary(Grid *gridin, int iSurface);

    void StoreRhsByResidual(Grid *gridIn, FieldProxy *rightHandSideProxy);
    void InitResidual(Grid *gridIn, FieldProxy *rightHandSideProxy);
    void RecoverResidual(Grid *gridIn, FieldProxy *rightHandSideProxy);
    
    FaceProxy    *CreateFaceProxy   (Grid *gridIn);
    GeomProxy    *CreateGeomProxy   (Grid *gridIn);
    
    StructIter *CreateInvStructIter(Grid *gridIn, int iSurface);
    
    bool GetInviscidFaceValue(Grid *gridIn, FaceProxy *faceProxy, StructIter *structIter, FieldProxy *qlProxy, FieldProxy *qrProxy, int nlen);
    bool LoadFaceFlux(Grid *gridIn, FaceProxy *faceProxy, StructIter *structIter, FieldProxy *fluxProxy, int nlen);

    //! Restart from perfect gas flow field for chemical reaction simulation.
    void ReadPrefectGasRestartH5(ActionKey *actkey);

    void DumpRestart(ActionKey *actkey);
    void ReadRestart(ActionKey *actkey);

    void DumpRestartH5(ActionKey *actkey);
    void ReadRestartH5(ActionKey *actkey);

    void AirForceCoef(ActionKey *actkey);

private:
    void AirForceCoefParts(ActionKey *actkey, int partID, int Coordinate = GlobalCoordinate);
    
    //! Write the variables of Cp, Cf, Qw, y+ etc. on surface to the file.
    void CpDistriCoef(ActionKey *actkey);
    
    //! To write variables on surface such as p, rho, T etc. to the file. added by Li Peng on Mar.5, 2019.
    //! Add this function in order to verify the chemical codes and the related algorithm of computing thermal distribution.
    void ExportSurfaceVariables(ActionKey *actkey);

    //! To obtain the derivatives of the designate primary variable whose sequence number in the array is marked with the valued of nIndex.
    void GetDerivatives(int nIndex, int iSurface, int ndim, int nlr, int i, int j, int k, int ni,int nj, int nk, RDouble4D &q,
                        RDouble4D &xvn, RDouble4D &yvn, RDouble4D &zvn, RDouble3D &vol, RDouble &dvdx, RDouble &dvdy, RDouble &dvdz);
    
    void GetResidual(ActionKey *actkey);

    void VisualizationAverageFlow(ActionKey *actkey);

public:
	RDouble4D *CompMachNumber(Grid *gridIn);
	RDouble4D *ComputeMachNumberField(Grid *gridIn);

	void ComputeVorticitybyQCriteria(Grid *gridIn, RDouble3D *&vorticity_x, RDouble3D *&vorticity_y, RDouble3D *&vorticity_z, RDouble3D *&vorticityMagnitude,  RDouble3D *&strain_rate, RDouble3D *&Q_criteria);
	void ComputeCp(Grid *gridIn, RDouble4D *&cp);

public:
    void VisualizationVariables(Grid *gridIn, FieldProxy *workProxy, const string &filename);

    void ComputeInterfaceGrad(int level);
    void UploadInterfaceData  (ActionKey *actkey);

    void DownloadInterfaceData(ActionKey *actkey);

	void UploadOversetData  (ActionKey *actkey);

	void DownloadOversetData(ActionKey *actkey);

	GradientCellCenter *GetGradientCellCenterUVWT(Grid *gridIn, FieldProxy *q);
	GradientCellCenter *GetGradientCellCenterUVWT(Grid *gridIn);

    void CalLimiters(Grid *gridIn, int ig, RDouble *limit, RDouble4D *q, RDouble **dqdx, RDouble **dqdy, RDouble **dqdz);

	//! Prepare post visualization variables, store it in to grid's postVisualization database.
	LIB_EXPORT void ComputePostVisualVariables(Post_Visual * postVisualization);

	//! Create and init control parameters.
	LIB_EXPORT void InitControlParameters();

 	//! Get controlParamters( a Para_NSSolver type pointer).
	LIB_EXPORT Param_NSSolverStruct *GetControlParameters() const;

    void ComputeHybridInterface();
};

}

