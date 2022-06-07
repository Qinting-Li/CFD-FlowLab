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
//! @file      PBS.h
//! @brief     register all PBSolver Objects.
//! @author    Ming Pingjian, Sun Huawen, Song Yang, Huang Jin.
 
#pragma once

//! For boundary conditions.
#include "PBBoundary.h"
#include "PBBoundaryFlowPressureOutlet.h"
#include "PBBoundaryFlowSolidSurface.h"
#include "PBBoundaryFlowVelocityInlet.h"
#include "PBBoundaryScalar.h"
#include "PBBoundaryScalarConvection.h"
#include "PBBoundaryScalarDefault.h"
#include "PBBoundaryScalarDiffusion.h"
#include "PBBoundaryTurb.h"
#include "PBBoundaryTurbDefault.h"
#include "PBBoundaryTurbEpsilonInlet.h"
#include "PBBoundaryTurbEpsilonPressureOutlet.h"
#include "PBBoundaryTurbEpsilonWall.h"
#include "PBBoundaryTurbKInlet.h"
#include "PBBoundaryTurbKPressureOutlet.h"
#include "PBBoundaryTurbKWall.h"
#include "PBBoundaryTurbSAInlet.h"
#include "PBBoundaryTurbSAWall.h"
#include "PBBoundaryTurbSAPressureOutlet.h"
#include "PBBoundaryEnergy.h"
#include "PBBoundaryEnergyDefault.h"
#include "PBBoundaryEnergyInlet.h"
#include "PBBoundaryEnergyOutlet.h"
#include "PBBoundaryEnergyWall.h"

//! For convection term.
#include "PBConvection.h"
#include "PBConvectionUDS.h"
#include "PBConvectionSUDS.h"
#include "PBConvectionCDS.h"
#include "PBConvectionQUICK.h"

//! For diffusion term.
#include "PBDiffusion.h"
#include "PBDiffusionNonorthogonal.h"

//! For transient term.
#include "PBTransient.h"
#include "PBTransientImplicitEuler.h"
#include "PBTransientImplicit2ndOrder.h"
#include "PBTransientCrankNicolson.h"

//! For source term.
#include "PBSource.h"
#include "PBSourceFlowDefault.h"
#include "PBSourceFlowGravity.h"
#include "PBSourceEnergyDefault.h"
#include "PBSourceEnergyDissipation.h"
#include "PBSourceTurbEpsilonDefault.h"
#include "PBSourceTurbKDefault.h"
#include "PBSourceTurbSADefault.h"

//! For gradient calculation.
#include "PBGrad.h"
#include "PBGradGauss.h"
#include "PBGradLsq.h"

//! For equation solvers.
#include "PBIterSolv.h"
#include "PBIterSolvBiCGSTAB.h"
#include "PBIterSolvGMRES.h"

//! For preconditioner.
#include "PBPreconditioner.h"
#include "PBPreconditionerDILU.h"

//! For solvers of governing equations.
#include "PBSolverFlow.h"
#include "PBSolverScalar.h"
#include "PBSolverEnergy.h"
#include "PBSolverTurbEpsilon.h"
#include "PBSolverTurbK.h"
#include "PBSolverTurbSA.h"
#include "PBSolverKEpsilon.h"

//! For post.
#include "PBOutput.h"
#include "PBOutputTecplot.h"


namespace PHSPACE
{
    //! Boundary conditions.
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBBoundary, PBBoundary, BOUNDARY_EMPTY);
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBBoundaryFlow, PBBoundaryFlowPressureOutlet, FLOW_PRESSURE_OUTLET);
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBBoundaryFlow, PBBoundaryFlowSolidSurface, FLOW_SOLID_SURFACE);
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBBoundaryFlow, PBBoundaryFlowVelocityInlet, FLOW_VELOCITY_INLET);
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBBoundary, PBBoundaryScalar, SCALAR_EMPTY);
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBBoundary, PBBoundaryScalarConvection, SCALAR_CONVECTION);
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBBoundary, PBBoundaryScalarDefault, SCALAR_DEFAULT);
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBBoundary, PBBoundaryScalarDiffusion, SCALAR_DIFFUSION);
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBBoundary, PBBoundaryTurb, TURB_EMPTY);
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBBoundary, PBBoundaryTurbDefault, TURB_DEFAULT);
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBBoundary, PBBoundaryTurbEpsilonInlet, TURB_EPSILON_INLET);
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBBoundary, PBBoundaryTurbEpsilonPressureOutlet, TURB_EPSILON_PRESSUREOUTLET);
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBBoundary, PBBoundaryTurbEpsilonWall, TURB_EPSILON_WALL);
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBBoundary, PBBoundaryTurbSAWall, TURB_SA_WALL);
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBBoundary, PBBoundaryTurbSAInlet, TURB_SA_INLET);
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBBoundary, PBBoundaryTurbSAPressureOutlet, TURB_SA_PRESSUREOUTLET);
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBBoundary, PBBoundaryTurbKInlet, TURB_K_INLET);
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBBoundary, PBBoundaryTurbKPressureOutlet, TURB_K_PRESSUREOUTLET);
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBBoundary, PBBoundaryTurbKWall, TURB_K_WALL);
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBBoundary, PBBoundaryEnergy, ENERGY_EMPTY);
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBBoundary, PBBoundaryEnergyDefault, ENERGY_DEFAULT);
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBBoundary, PBBoundaryEnergyInlet, ENERGY_INLET);
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBBoundary, PBBoundaryEnergyOutlet, ENERGY_OUTLET);
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBBoundary, PBBoundaryEnergyWall, ENERGY_WALL);

    //! Calculations for Convection.
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBConvection, PBConvection, CONVECTION_EMPTY);
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBConvection, PBConvectionUDS, CONVECTION_UPWIND);
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBConvection, PBConvectionCDS, CONVECTION_SUDS);
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBConvection, PBConvectionCDS, CONVECTION_CDS);
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBConvection, PBConvectionQUICK, CONVECTION_QUICK);

    //! Calculations for Diffusion.
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBDiffusion, PBDiffusion, DIFFUSION_EMPTY);
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBDiffusion, PBDiffusionNonorthogonal, DIFFUSION_NON_ORTHOGONAL);

    //! Calculations for Transient.
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBTransient, PBTransient, TRANSIENT_EMPTY);
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBTransient, PBTransientImplicitEuler, TRANSIENT_IMPLICIT_EULER);
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBTransient, PBTransientImplicit2ndOrder, TRANSIENT_IMPLICIT_2ND_ORDER);
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBTransient, PBTransientCrankNicolson, TRANSIENT_CRANK_NICOLSON);

    //! Calculations for Source.
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBSource, PBSource, SOURCE_EMPTY);
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBSource, PBSourceFlowDefault, SOURCE_FLOW_DEFAULT);
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBSource, PBSourceFlowGravity, SOURCE_FLOW_GRAVITY);
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBSource, PBSourceEnergyDefault, SOURCE_ENERGY_DEFAULT);
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBSource, PBSourceEnergyDissipation, SOURCE_ENERGY_DISSIPATION);
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBSource, PBSourceTurbEpsilonDefault, SOURCE_TURB_EPSILON_DEFAULT);
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBSource, PBSourceTurbKDefault, SOURCE_TURB_K_DEFAULT);
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBSource, PBSourceTurbKDefault, SOURCE_TURB_SA_DEFAULT);

    //! Calculations for Grad.
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBGrad, PBGrad, GRAD_EMPTY);
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBGrad, PBGradGauss, GRAD_GAUSS);
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBGrad, PBGradLsq, GRAD_LSQ);

    //! Solver for systems of linear equations.
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBIterSolv, PBIterSolv, ITERSOLV_EMPTY);
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBIterSolv, PBIterSolvBiCGSTAB, ITERSOLV_BiCGSTAB);
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBIterSolv, PBIterSolvGMRES, ITERSOLV_GMRES);

    //! Solver for Equation.
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBSolver, PBSolverFlow, SOLVER_FLOW);
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBSolver, PBSolverScalar, SOLVER_SCALAR);
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBSolver, PBSolverTurbEpsilon, SOLVER_TURB_EPSILON);
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBSolver, PBSolverTurbK, SOLVER_TURB_K);
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBSolver, PBSolverKEpsilon, SOLVER_TURB_KEPSILON);
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBSolver, PBSolverTurbSA, SOLVER_TURB_SA);
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBSolver, PBSolverEnergy, SOLVER_ENERGY);

    //! Solver for Preconditioner.
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBPreconditioner, PBPreconditioner, ITERSOLVPRECOND_EMPTY);
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBPreconditioner, PBPreconditionerDILU, ITERSOLVPRECOND_DILU);

    //! Solver for Output.
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBOutput, PBOutput, OUTPUT_EMPTY);
    DEFINE_PBSOLVER_REGISTER_OBJECT(PBOutput, PBOutputTecplot, OUTPUT_TECPLOT);

}
