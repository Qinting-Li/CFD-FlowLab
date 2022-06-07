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
//! @file      PBSolverFlow.h
//! @brief     Pressure Based flow solver.
//! @author    Ming Pingjian, Sun Huawen, Song Yang.

#pragma once

#include <map>
#include <iostream>

#include "PBObject.h"
#include "PBSolver.h"

namespace PHSPACE
{
    class PBGrad;
	class PBBoundary;
	class PBTransient;
	class PBConvection;
	class PBDiffusion;
	class PBSource;
	class UnstructGrid;

#define SOLVER_FLOW FLOW

class PBSolverFlow : public PBSolver
{
public:
    //! Pointer for discrete scheme.
    map<string, PBGrad *> gradCalc; //! Gradient calculation
    vector<PBBoundary *> boundCalc; //! Boundary condition calculation
    PBTransient *tranCalc; //! Transient calculation
    PBConvection *convCalc; //! Convection calculation
    PBDiffusion *diffCalc; //! Diffusion term calculation
    vector<PBSource *> sourceCalc; //! Source term calculation

protected:
    //! From configuration parameters
    RDouble initRho; //! density.
    RDouble initMu; //! viscosity.
    RDouble initU; //! Init velocity in U direction.
    RDouble initV; //! Init velocity in V direction.
    RDouble initW; //! Init velocity in W direction.
    RDouble initP; //! Init P.
    RDouble urfU; //! Relaxing factor of velocity in U direction.
    RDouble urfV; //! Relaxing factor of velocity in V direction. 
    RDouble urfW; //! Relaxing factor of velocity in W direction. 
    RDouble urfP; //! Relaxing factor of pressure.
    RDouble urfFlux; //! Relaxing factor of flux.
    RDouble resU; //! Residual of velocity in U direction.
    RDouble resV; //! Residual of velocity in V direction.
    RDouble resW; //! Residual of velocity in W direction.
    RDouble resP; //! Residual of pressure.
    int bodyForceFlag; //! body force flag

protected:
    //! Global register
    RDouble *u; //! velocity in U direction.
    RDouble *v; //! velocity in V direction.
    RDouble *w; //! velocity in W direction.
    RDouble *p; //! pressure.
    RDouble *uc; //! correction of velocity in U direction.
    RDouble *vc; //! correction of velocity in V direction.
    RDouble *wc; //! correction of velocity in W direction.
    RDouble *pc; //! correction of pressure.
    RDouble *apU; //! diagonal coefficient in U direction.
    RDouble *apV; //! diagonal coefficient in V direction.
    RDouble *apW; //! diagonal coefficient in W direction.
    RDouble *ap; //! diagonal coefficient of general transport equation.
    RDouble *ajLeft; //! Upper triangular coefficient.
    RDouble *ajRight; //! Lower triangular coefficient.
    RDouble *scU; //! right-hand side in U direction.
    RDouble *scV; //! right-hand side in V direction.
    RDouble *scW; //! right-hand side in W direction.
    RDouble *sc; //! right-hand side of general transport equation.
    RDouble *bfU; //! bodyForce in U direction.
    RDouble *bfV; //! bodyForce in V direction.
    RDouble *bfW; //! bodyForce in W direction.

protected:
    //! Temp value in calculation.
    RDouble *uOld; //! velocity in U direction at last time step.
    RDouble *vOld; //! velocity in V direction at last time step.
    RDouble *wOld; //! velocity in W direction at last time step.
    RDouble *pOld; //! pressure at last time step.
    RDouble *uOldOld;
    RDouble *vOldOld;
    RDouble *wOldOld;
    RDouble *pOldOld;
    RDouble *uf;
    RDouble *ub; //! boundary velocity in U direction.
    RDouble *vb; //! boundary velocity in V direction.
    RDouble *wb; //! boundary velocity in W direction.
    RDouble *pb; //! boundary pressure.
    RDouble *pcb; //! boundary pressure correction.
    RDouble *faceFlux; //! Pointer to initial position of flux on internal cell surface.
    RDouble *flux; //! flux on face.
    RDouble *dun; //! relaxation of flux on face.
    RDouble *rho; //! density.
    RDouble *rhof; //! density on face.
    RDouble *mu; //! viscosity.
    RDouble *hu; //! for pressure correction equation.
    RDouble *hv;
    RDouble *hw;
    RDouble *dudx; //! gradient in X direction.
    RDouble *dudy;
    RDouble *dudz;
    RDouble *dvdx; //! gradient in Y direction.
    RDouble *dvdy;
    RDouble *dvdz;
    RDouble *dwdx; //! gradient in Z direction.
    RDouble *dwdy;
    RDouble *dwdz;
    RDouble *dpdx; //! gradient of pressure.
    RDouble *dpdy;
    RDouble *dpdz;
    RDouble *dpcdx; //! gradient of pressure correction.
    RDouble *dpcdy;
    RDouble *dpcdz;
    RDouble uResNow; //! residual in U direction.
    RDouble vResNow; //! residual in V direction.
    RDouble wResNow; //! residual in W direction.
    RDouble pResNow; //! residual of pressure
    RDouble dt;
protected:
    //! Solver iter Pointer
    PBIterSolv *iterSolvU; //! in U direction.
    PBIterSolv *iterSolvV; //! in V direction.
    PBIterSolv *iterSolvW; //! in W direction.
    PBIterSolv *iterSolvP; //! for pressure
    int maxSweepU; //! max iter steps in U direction.
    int maxSweepV; //! max iter steps in V direction.
    int maxSweepW; //! max iter steps in W direction.
    int maxSweepP; //! max iter steps of pressure
    RDouble iterSolvTolU; //! tolerance of iteration in U direction.
    RDouble iterSolvTolV; //! tolerance of iteration in V direction.
    RDouble iterSolvTolW; //! tolerance of iteration in W direction.
    RDouble iterSolvTolP; //! tolerance of iteration for pressure

protected:
    UnstructGrid *grid;
    RDouble *deltL; //! weighting coefficient of face. 
    RDouble *deltR; //! weighting coefficient of face. 
    int nTotalInFaces; //! number of interior faces.

    inline const int &NTotalInFace()
    {
        return this->nTotalInFaces;
    }

    vector<string> phiNameList; //! name set of variables.

public:
    PBSolverFlow(const string &name);
    ~PBSolverFlow();
public:
    virtual void AllocateGlobalVar(Grid *grid);
    virtual void DeAllocateGlobalVar(Grid *grid);
    virtual void ZeroResiduals(Grid *grid);
    virtual void TimeStep(Grid *grid);
    virtual void ComputeResidual(Grid *grid, FieldProxy *rhs_proxy);
    virtual void SourceFlux(Grid *grid);
    virtual void InviscidFlux(Grid *grid);
    virtual void ViscousFlux(Grid *grid);
    virtual void RungeKutta(Grid *grid, FieldProxy *rhs_proxy);
    virtual void Post();
    virtual void PostSolve(Grid *grid, int stage, int level);
    virtual void InitIncompressibleResultFile(Grid *grid, string fileName);
    virtual void DumpIncompressibleResultFile(Grid *grid, string fileName);
    virtual RDouble UnsteadyConvergence(Grid *grid);

public:
    virtual void AllocateMemory();
    virtual void InitField();
    virtual void ReinitField();
    virtual void UpdateUnsteadyVariable(Grid *grid);
    virtual void SolutionIsConverged(int iZone, bool *flag);
    virtual void OutputResidual(ofstream &file, int flag);

    virtual void InitMemory();

public:
    void UpdatePhi(string phiName, int cellBegin, int cellEnd, RDouble *phic, RDouble *phi);

public:
    virtual void DiscrtMomentumEquations();
    virtual void SolveMomentumEquations(Region &region);
    virtual void DiscrtPressureCorrectionEquation();
    virtual void SolvePressureCorrectionEquation(Region &region);
    virtual void UpdateVelocity();
    virtual void UpdatePressureAndVelocity();
    virtual void UpdateBoundary();
    virtual void ComputePostVisualVariables(Post_Visual *postVisualization);
    virtual void ComputePostProbesVariables(Post_Probes *postProbesVar);
    virtual void UpdateProperties();

protected: //!SIMPLEC
    void UpdateBoundaryPressureAndVelocity();
    void UpdateCorrectionBoundaryPressureAndVelocity();
    void CalcPressureGrad();
    void CalcVelocityGrad();
    void CalcVelocityBoundaryConditions();
    void CalcVelocityTransientTerm();
    void CalcVelocityConvectionTerm();
    void CalcVelocityDiffusionTerm();
    void CalcVelocitySourceTerm();
    void PrimeVelocityEquations();
    void CalcVelocityRes();
    void RelaxVelocity();
    void CoefPrepForPres();
    void PrimePressureEquations();
    void CalcPressureRes();
    void CalcPressureCorrectionGrad();
    void UpdateVelocityAndFlux();
    void CheckFaceCells();
};

DECLARE_PBSOLVER_REGISTER_FUNCTION(PBSolver, PBSolverFlow, SOLVER_FLOW);
}
