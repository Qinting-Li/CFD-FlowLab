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
//! @file      Simulation.h
//! @brief     Explain this file briefly.
//! @author    PHengLEI team.

#pragma once
#include "Precision.h"

namespace PHSPACE
{
class Zone;
class TimeSpan;
class Region;

class Simulation
{
private:
    Region *region;
public:
    Simulation(int dimension);
    ~Simulation();
public:
    //! select task according to taskType
    void Start();
    //! main part of simulation
    void SolveSimulation();

    //! The numerical iteration process similar to the function named SolveSimulation.
    void SelfAdaptionSolve();

    //! post process of simulation
    void PostSimulation();

    //! initionalize flow field in the first flowInitStep steps with coarse grid when multigrid simulation
    void MultiGridInitFlow();
private:
    //! read the grid such as grid coordinate, boundary information, wall distance for turbulent flow and interpoint information
    void ReadGeometry();
    void InitGeometry();
    void SwapGeommetry();
    void InitOversetGrid();
    void OversetGridConfig();

    //! Define solver names for all zones and create the newly defined solvers.
    //! IMPORTANT for solver developers:
    //! Define your solver name which is same with the solver class name here!
    void DefineAndCreateSolvers();

    //! Initialize each solver on each zone.
    void InitializeSolvers();

    //! Release the dynamic memories related to the solvers.
    void RefreshSolvers();

    //main work flow of CFD simulation process
    void RunCFD();

    //! Check whether the current computation is restart process.
    bool IsRestartProcess();

    //run CFD using the self-adaption method.
    void SelfAdaptionRunCFD();
    void FirstStepRunCFD(bool isStart = true);
    void SecondStepRunCFD(bool isStart = false);
    void ThirdStepRunCFD(bool isStart = false);
    void FinalStepRunCFD();

    void RunIntegrativeSolver();
    void RunGridCheckSolver();

    #ifdef USE_INHOUSE_VERSION
    //! work flow of some specific simulation process
    void RunDSMC();
    void RunUGKS();
    void RunHelicopterSolver();

     //! Identify the motive derivatives from the unsteady result file "./results/aircoef.dat"
    void IdentifyDerivative();
    #endif

    #ifdef USE_SpecDiffHybSolver
    void RunSpecDiffHybSolver();
    void RunSpecSolver();
    #endif

    //! compute wall distance for turbulence flow while there was no wdt files.
    void ComputeWallDistance();

    void IsNeedConvert();
    void UpdateAllZonesPara();

    //! run post processing without flow field simulation
    void RunPostProcessing();

    void SolveOneOuterStep();
    void UpdateOuterIterationStep(int &outerIterStep);

    void PreSolve();
    void MainSolve();
    void PostSolve();

    void InitGlobalBoundaryCondition();
    
    void SolveIncomFlow();
    void SolveUnsteady();
    void SolveSteady();
    
    bool IsSubIterationConvergence(const int &innerIterStep, const int &minSubIteration, const int &maxSubIteration, const RDouble &tolalSubIteration, const RDouble &averageQuotientMaxInSolvers);
};
}
