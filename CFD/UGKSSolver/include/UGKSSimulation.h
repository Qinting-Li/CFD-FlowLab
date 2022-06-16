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
//! @file      UGKSSimulation.h
//! @brief     It defines UGKS simulation's frame function.
//! @author    

#pragma once

namespace PHSPACE
{

class Zone;
class TimeSpan;
class Region;
//! @brief UGKSSimulation defines UGKS simulation's frame function, it calls functions of class UGKSSolver.
class UGKSSimulation
{
private:
    Region * region;
    TimeSpan * simu_ts;
public:
    UGKSSimulation();
    ~UGKSSimulation();
public:
    //! UGKSSolver's main function.
    void Run();
private:
    //! Preprocess.
    void Preprocess();

    //! Initialize memory.
    void InitGeometry();

    //! Read geometry information.
    void ReadGeometry();

    //! Read grid information.
    void ReadGrid();

    //! Update parameters of all zones.
    void UpdateAllZonesPara();

    //! Define and creat solver.
    void DefineAndCreateSolvers();

    //! Initialize some common parameters.
    void InitCommonParameter();

    //! Initialize solvers.
    void InitializeSolvers();

    //! Initialize global boundary condition.
    void InitGlobalBoundaryCondition();

    //! Solve simulation.
    void SolveSimulation();

    //! Solve simulation of one step.
    void SolveOneStep();

    //! Solve main flow field evolution.
    void MainSolve();
   
    //! Computer minimum time step.
    void ComputeMinTimeStep();

    //! Solve main flow field evolution of UGKS.
    void UGKSSolve();

    //! Solve proprocess of field in one step.
    void PostSolve();

    //! Visualize flow field.
    void Visualize();

    //! Postprocess when flow converges to an ideal state.
    void Postprocess();
};

}