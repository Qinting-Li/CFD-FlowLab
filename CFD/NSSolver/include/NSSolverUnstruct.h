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
//! @file      NSSolverUnstruct.h
//! @brief     NS solver for unstructured grid.
//! @author    Bell, He Xin, Zhang Laiping, Dr. Wang, Zhang Yang, He Kun, 
//!            Wan Yunbo, Xu Gang, Meng Liyuan, Zhang Yaobing.

#pragma once
#include "NSSolver.h"
#include "PHMatrix.h"
#include "Limiter.h"

namespace PHSPACE
{
class GeomProxy;
class Gradient;
class Param_NSSolverUnstruct;

class NSSolverUnstruct : public NSSolver
{
private:
    Gradient * gradientNSField;
    Gradient * gradientTemperature;
    int systemgridtype;
    Limiter *limiterOnlyForMixGrid;
    FaceProxy *InviscidfaceProxyOnlyForMixGrid;
    FaceProxy *ViscousfaceProxyOnlyForMixGrid;
public:
    NSSolverUnstruct();
    ~NSSolverUnstruct();
public:
    void ReadParameter();

    void AllocateGlobalVar(Grid *gridIn);
    void DeAllocateGlobalVar(Grid *gridIn);

    //! Judge if the flow field variables file already exists in NS solver.
    bool JudgeIfRestart();
    bool JudgeIfReadAverage();

    void DumpRestart(ActionKey *actkey);
    void ReadRestart(ActionKey *actkey);

    void DumpRestartH5(ActionKey *actkey);
    void ReadRestartH5(ActionKey *actkey);

    void InitFlowAsRestart();
    void LoadQ(Grid *gridIn, FieldProxy *q_proxy);
    void RestrictDefect(Grid *fgrid_in, Grid *cgrid_in);
    void PutCorrectionBack(Grid *gridIn, FieldProxy *q_proxy);

    //! Prolongation operator, from coarse mesh to fine mesh.
    void CorrectFineGrid(Grid *fgrid_in, Grid *cgrid_in);

    //! Prolongation by the current coarse cell directly.
    void CorrectFineGridZeroOrder(Grid *fgrid_in, Grid *cgrid_in);

    //! Prolongation by the current and the neighbor coarse cells.
    void CorrectFineGridHighOrder(Grid *fgrid_in, Grid *cgrid_in);

    void InterpolatFineGrid(Grid *fgrid_in, Grid *cgrid_in);

#ifdef USE_INHOUSE_VERSION
    //! Compute the Matrix of left hand side for time integration on unstructured grid:\n
    //! @param[in] gridIn    computational mesh.
    void LHSMatrixs(Grid *gridIn);

    //! Compute the Matrix of left hand side for time integration on unstructured grid:\n
    //! @param[in] gridIn    computational mesh.
    void LineLUSGSMatrixs(Grid* gridIn);
#endif

    void LoadResiduals(Grid *gridIn, FieldProxy *rhsProxy);   

    void Boundary(Grid *gridIn);
    void ComputeMassFlow(Grid *gridIn);
    void ComputeMassInFlow(Grid *gridIn, UnstructBC *bcRegionUnstruct);
    void ComputeMassOutFlow(Grid *gridIn, UnstructBC *bcRegionUnstruct);
    void ComputeGamaAndTemperature(Grid *gridIn);

    //! Compute the node value in the flow.
    //! @param[in] gridIn         the mesh in the computing.
    void ComputeNodeValue(Grid *gridIn);

    //! Modify the node value for the interpoint.
    //! @param[in] gridIn         the mesh in the computing.
    void ModifyNodeValue(Grid *gridIn);

    void GetQlQr  (Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd);
    void GetGamaLR(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd);
    void ReConstructFaceValue(Grid *gridIn, FaceProxy *faceProxy, Limiter * limiter, int localStart, int localEnd);
    void ZeroResiduals(Grid *gridIn);
    void SpectrumRadiusInviscid(Grid *grid);
    void PreconInvSpectrumRadius(Grid *gridIn);
    void SpectrumRadiusViscous(Grid *grid);
    void Diagonal(Grid *gridIn);
    void ComputeMinTimeStep(Grid *grid, RDouble & minDt, RDouble & maxDt);
    void ReduceMaxTimeStep(Grid *grid, RDouble globalMinDt);
    void LocalTimeStep(Grid *grid);

    //! calculate time step by another methods different with PHengLEI code.
    void LocalTimeStep_Prince(Grid *grid);
    
    void GlobalTimeStep(Grid *grid);
    void LocalGlobalTimeStep(Grid *grid);
    void RungeKuttaResidual(Grid *gridIn, FieldProxy *dqProxy, RDouble coef);
    void LoadFlux(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd);

    void ZeroResidualOfSpecialCells(Grid *gridIn);

    void UpdateQlQrOnlyforMixGrid(Grid *grid);
    void UpdateInviscidfluxOnlyforMixGrid(Grid *grid);
    void UpdateViscousfluxOnlyforMixGrid(Grid *grid);
    void LoadGlobalfacefluxtoResOnlyforMixGrid(Grid *gridIn);
    void RestrictAllQ(Grid *fgrid, Grid *cgrid);
    void ComputeViscousCoeff(Grid *gridIn);
    void ComputeViscousCoefficientWithoutChemical(Grid *gridIn);

    void BoundaryQlQrFix(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd);
    void InitCGrid(Grid *fgrid_in, Grid *cgrid_in);
    RDouble * CompSoundSpeed(Grid *gridIn);
    RDouble * CompMachNumber(Grid *gridIn);
    RDouble** ComputeDimensionalVariables(Grid *gridIn);

    RDouble * ComputeVorticitybyQCriteria(Grid *gridIn);
    void ComputeVorticitybyQCriteria(Grid *gridIn, RDouble * vorticity_x, RDouble * vorticity_y, RDouble * vorticity_z, RDouble * vorticityMagnitude,  RDouble * strain_rate,  RDouble * Q_criteria);
    RDouble * ComputeCp(Grid *gridIn);
    RDouble * ComputeCpAverage(Grid *gridIn);

    void ComputeGradient(Grid *grid);
    void ComputelimiterOnlyForMixGrid(Grid *grid);
    void InitiallimiterOnlyForMixGrid(Grid *grid);
    void InitialInviscidfaceproxyOnlyForMixGrid(Grid *grid);
    void InitialViscousfaceproxyOnlyForMixGrid(Grid *grid);

    void FillField(Grid *gridIn, FieldProxy *field_proxy, RDouble value);
    void FillField(Grid *gridIn, FieldProxy *fieldTarget, FieldProxy *fieldSource);
    FieldProxy * CreateFieldProxy(Grid *gridIn);
    FieldProxy * GetFieldProxy(Grid *gridIn, const string &field_name);
    FieldProxy * GetResidualProxy(Grid *gridIn);

    void RecoverResidual(Grid *gridIn, FieldProxy *rhsProxy);
    void InitResidual(Grid *gridIn, FieldProxy *rhsProxy);
    void StoreRhsByResidual(Grid *gridIn, FieldProxy *rhsProxy);

    FaceProxy     * CreateFaceProxy    (Grid *gridIn);
    GeomProxy     * CreateGeomProxy    (Grid *gridIn);
    NSFaceValue   * CreateNSFaceValue  (Grid *gridIn);

    void FillGeomProxy(Grid *gridIn, GeomProxy *geom_proxy, int localStart, int localEnd);

    void CalPressureFactor(Grid *gridIn);
    void GetHeatTransferCoeff(Grid *gridIn);

    void FreeGradientProxy(Grid *grid);
    void InviscidFlux  (Grid *gridIn);
    void ViscousFlux   (Grid *gridIn);
    void ChemicalSource(Grid *gridIn);
    void DualTimeSource(Grid *gridIn);
    void StoreHistoricalResidual(Grid *gridIn);
    void UpdateUnsteadyFlow(Grid *gridIn);
    RDouble UnsteadyConvergence(Grid *gridIn);

    void Turb_Sengy(Grid *gridIn);
    void GetInviscidFaceValue(Grid *gridIn, FaceProxy *faceProxy, Limiter * limiter, int localStart, int localEnd);
    void GetVisFaceValue(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd);
    void ComputeVisflux(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd);

    void ComputeFaceWeight(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd);
    void ComputeInviscidFlux(Grid *grid, FaceProxy *faceProxy, int localStart, int localEnd);
    void SolutionFix(Grid *gridIn, RDouble *prim, int icell);

    void UpdateFlowField(Grid *gridIn, FieldProxy *q_proxy, FieldProxy *dqProxy);
    void SetGhostDQLUSGS(Grid *gridIn, RDouble **dq, RDouble **q);
    void SetSymmBCGhostDQLUSGS(Grid *gridIn, UnstructBC *bcRegionUnstruct, RDouble **dq, RDouble **q);
    void SetWallBCGhostDQLUSGS(Grid *gridIn, UnstructBC *bcRegionUnstruct, RDouble **dq, RDouble **q);
    void SetExtraBCGhostDQLUSGS(Grid *gridIn, UnstructBC *bcRegionUnstruct, RDouble **dq, RDouble **q);
    
    void SwapNeighborsDQNonblocking(Grid *gridIn, FieldProxy *dqProxy);
    void CompressDQ(DataContainer *&dataContainer, FieldProxy *fieldProxy, Grid *grid, const int &zoneIndex, const int &neighborZoneIndex, const int &neqn);
    void DecompressDQ(DataContainer *&dataContainer, FieldProxy *fieldProxy, Grid *grid, const int &zoneIndex, const int &neighborZoneIndex, const int &neqn);
    RDouble GetRad(RDouble *prim, RDouble *nxyz, RDouble gama);
    void SolveLUSGSForward(Grid *grid, FieldProxy * dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal);
    void SolveLUSGSBackward(Grid *grid, FieldProxy * dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal);

#ifdef USE_INHOUSE_VERSION
    //! Forward sweep in line implicit LUSGS time methods.
    void SolveLineLUSGSForward(Grid* grid, FieldProxy* dqProxy, FieldProxy* LUplusDQ, RDouble& sweepNormal);

    //! Backward sweep in line implicit LUSGS time methods.
    void SolveLineLUSGSBackward(Grid* grid, FieldProxy* dqProxy, FieldProxy* LUplusDQ, RDouble& sweepNormal);

    //! This is for Block Lower-Upper Symmetric Gauss¨CSeidel Scheme according wzj reference in 2000 year.
    //! Forward sweep in BLUSGS time methods.
    //! @param[in] grid    computational mesh.
    //! @param[in] dqProxy        a proxy to store dq.
    //! @param[in] LUplusDQ       a proxy to store dq.
    //! @param[in] sweepNormal    a value depends the number of inner iterations.
    void SolveBLUSGSForward(Grid *grid, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal);
    
    //! Backward sweep in BLUSGS time methods.
    void SolveBLUSGSBackward(Grid *grid, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal);
#endif

    void AirForceCoef(ActionKey *actkey);
    void AirForceCoefParts(ActionKey *actkey, int partID, int Coordinate = GlobalCoordinate);
    void AirForceCoefComponents(ActionKey *actkey, int partID, int Coordinate = GlobalCoordinate);
    void CpDistriCoef(ActionKey *actkey);
    void GetResidual(ActionKey *actkey);

    //! Dump out typical standard case: turbulent flat plate.
    void Turbulence_Flat_Plate_Output_NASA(Grid *gridIn);
    void CreateNodeVariable(Grid *gridIn, RDouble **nodeVariable);

    void ComputeResidual(ActionKey *actkey);
    void VisualizationAverageFlow(ActionKey *actkey);

    void ComputeAverageVisualNodeField(Grid *gridIn, RDouble **qn);
    void BoundaryVisualization(Grid *gridIn, ActionKey *actkey, RDouble **qn, int nvarplot);

    void BoundaryVisualization(Grid *gridIn, std::ostringstream &oss, vector<string> &title, RDouble **qn, int nvarplot);
    void SliceVisualization(Grid *gridIn, ActionKey *actkey, RDouble **qn, int nvarplot);

    void SpectrumRadius(Grid *gridIn);

    void RestrictDQ(RDouble *dqtry, int nm, RDouble lim_max);

    void UploadInterfaceData(ActionKey *actkey);
    void DownloadInterfaceData(ActionKey *actkey);

    //! Upload the interpoint variables for send.
    //! @param[in] actkey         the actkey store the information of the interpoint.
    void UploadInterpointData(ActionKey *actkey);

    //! Download the interpoint variables from the receive data.
    //! @param[in] actkey         the actkey store the information of the interpoint.
    void DownloadInterpointData(ActionKey *actkey);

    void CommunicationInterpointWeight(ActionKey *actkey);
    void DownloadInterpointWeight(ActionKey *actkey);

    void UploadOversetData  (ActionKey *actkey);
    void DownloadOversetData(ActionKey *actkey);

    //Gradient * GetGradientNSField(Grid *gridIn);
    //Gradient * GetGradientTemperature(Grid *gridIn);

    void GetGradientField(Grid *gridIn);

    //! Prepare post visualization variables, store it in to postVisualization database.
    LIB_EXPORT void ComputePostVisualVariables(Post_Visual * postVisualization);
    //! Prepare monitored probes variables for post-processing , store it into postProbesVar database.
    LIB_EXPORT void ComputePostProbesVariables(Post_Probes * postProbesVar);
    //! Create and init ControlParameters
    LIB_EXPORT void InitControlParameters();

    //! Get controlParamters(a Para_NSSolver type pointer)
    LIB_EXPORT Param_NSSolverUnstruct * GetControlParameters() const;

    void CompressSpecificArrayToInterface(DataContainer *&dataContainer, const string &fieldName, Grid *grid, const int &neighborZoneIndex, const int &neqn);
    void DecompressArrayFromInterface(DataContainer *&dataContainer, const string &fieldName, Grid *grid, const int &neighborZoneIndex, const int &neqn);

    void CompressGradientandLimitToInterfaceOnlyForMixGrid(DataContainer *&dataContainer, Grid *grid, const int &neighborZoneIndex);
    void DecompressGradientandLimitToInterfaceOnlyForMixGrid(DataContainer *&dataContainer, Grid *grid, const int &neighborZoneIndex);

    void CompressQlQrToInterfaceOnlyForMixGrid(DataContainer *&dataContainer, Grid *grid, const int &neighborZoneIndex);
    void DecompressQlQrToInterfaceOnlyForMixGrid(DataContainer *&dataContainer, Grid *grid, const int &neighborZoneIndex);

    void CompressFacefluxToInterfaceOnlyForMixGrid(DataContainer *&dataContainer, Grid *grid, const int &neighborZoneIndex);
    void DecompressFacefluxToInterfaceOnlyForMixGrid(DataContainer *&dataContainer, Grid *grid, const int &neighborZoneIndex);

    //! Out flow boundary condition for one Region of a bounday.
    //! A BC Region Contains lots of faces with the same boundary name.
    //! @param[in] gridIn              an unstruct grid pointer.
    //! @param[in] bcRegionUnstruct    an object of some BCRegionUnstruct class, contain the information of this BC.
    void OutflowBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct);

    //! Symmetry boundary condition for one Region of a bounday.
    //! @param[in] gridIn              an unstruct grid pointer.
    //! @param[in] bcRegionUnstruct    an object of some BCRegionUnstruct class, contain the information of this BC.
    void SymmetryBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct);
    
    //! Viscous Adiabatic wall boundary condition for one Region of a bounday.
    //! @param[in] gridIn              an unstruct grid pointer.
    //! @param[in] bcRegionUnstruct    an object of some BCRegionUnstruct class, contain the information of this BC.
    void ViscousAdiabaticWallBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct);
    
    //! Viscous Isotropic wall boundary condition for one Region of a bounday.
    //! @param[in] gridIn              an unstruct grid pointer.
    //! @param[in] bcRegionUnstruct    an object of some BCRegionUnstruct class, contain the information of this BC.
    void ViscousIsotropicWallBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct);
    
    //! Farfield boundary condition for one Region of a bounday.
    //! @param[in] gridIn              an unstruct grid pointer.
    //! @param[in] bcRegionUnstruct    an object of some BCRegionUnstruct class, contain the information of this BC.
    void FarfieldBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct); 
    
    //! Farfield  boundary condition for one Region of a bounday with Precon for subsonic flows.
    //! @param[in] gridIn              an unstruct grid pointer.
    //! @param[in] bcRegionUnstruct    an object of some BCRegionUnstruct class, contain the information of this BC.
    void PreconFarfieldBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct);   
    
    //! Inflow boundary condition for one Region of a bounday.
    //! @param[in] gridIn              an unstruct grid pointer.
    //! @param[in] bcRegionUnstruct    an object of some BCRegionUnstruct class, contain the information of this BC.
    void InflowBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct);   
    
    //! PressureOutlet boundary condition for one Region of a bounday.
    //! @param[in] gridIn              an unstruct grid pointer.
    //! @param[in] bcRegionUnstruct    an object of some BCRegionUnstruct class, contain the information of this BC.
    void PressureOutletBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct);   
    
    //! PressureInlet boundary condition for one Region of a bounday.
    //! @param[in] gridIn              an unstruct grid pointer.
    //! @param[in] bcRegionUnstruct    an object of some BCRegionUnstruct class, contain the information of this BC.
    void PressureInletBCRiemannRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct); 
    
    //calculate the total area of all BC regions with the same BC name.
    void calculateBCArea(Grid *gridIn, UnstructBC *bcRegionUnstruct); 

    //calculate the mass flux ratio of all mass outlet BC regions with the same BC name.
    void calculateMassFluxRatio(Grid *gridIn, UnstructBC *bcRegionUnstruct); 

    // mass flow in boundary condition.
    void MassFlowInletBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct); 
    
    // mass flow out boundary condition.
    void MassFlowOutletBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct);
private:
    void InviscidFluxWrap(Grid *gridIn, FaceProxy *faceProxy, INVScheme inv_scheme, int localStart, int localEnd);
    void RegisterCFDSolverInterfaceField();
    //! Register the interpoint field in the NS unstructured solver.
    void RegisterCFDSolverInterpointField();
    void RegisterCFDSolverOversetField();
    RDouble P0Wall(RDouble p1, RDouble gama, RDouble mach);
    RDouble R0Wall(RDouble r1, RDouble gama, RDouble mach);
    Limiter * CreateLimiter(Grid *gridIn);
};

}
