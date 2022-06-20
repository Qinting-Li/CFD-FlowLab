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
//! @file      PBSolverTurbK.h
//! @brief     TurbK.
//! @author    Wan Yunbo, Ming Pingjian, Meng Liyuan, Sun Huawen, Song Yang, Huang Jin, Bell..

#pragma once

#include "PBSolverScalar.h"
#include "PBSolverTurbK.h"
#include "PBSolverTurbEpsilon.h"

using namespace std;

namespace PHSPACE
{
#define SOLVER_TURB_KEPSILON TURB_KEPSILON

class PBSolverKEpsilon : public PBSolverScalar
{
public:
    PBSolverKEpsilon(const std::string &name);
    virtual ~PBSolverKEpsilon();

    void InitMemory();
    void ReinitField();
    void UpdateUnsteadyVariable(Grid *grid);
    void SolveKEquation(Region &region);
    void SolveEpsilonEquation(Region &region);
    void UpdateAfterIterloop();
    void ComputePostVisualVariables(Post_Visual *postVisualization);
    void InitIncompressibleResultFile(Grid *grid, string fileName);
    void DumpIncompressibleResultFile(Grid *grid, string fileName);
    void SolutionIsConverged(int iZone, bool *flag);
    void SetGeometry(PHGeometry* geometry) { this->geometry = geometry; }

protected:
    vector<PHSolver *> solverList;
    CFDSolver *solverKEqn;
    CFDSolver *solverEpsilonEqn;
    PHGeometry *geometry;
};

DECLARE_PBSOLVER_REGISTER_FUNCTION(PBSolver, PBSolverKEpsilon, SOLVER_TURB_KEPSILON);
}