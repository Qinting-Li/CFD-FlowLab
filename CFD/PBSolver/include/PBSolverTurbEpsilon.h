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
//! @file      PBSolverTurbEpsilon.h
//! @brief     TurbEpsilon of Pressure Based Solver (PBSolver).
//! @author    Wan Yunbo, Ming Pingjian, Meng Liyuan, Sun Huawen, Song Yang, Huang Jin, Bell..

#pragma once

#include <vector>
#include "PBSolverScalar.h"

using namespace std;

namespace PHSPACE
{
#define SOLVER_TURB_EPSILON TURB_EPSILON

class PBSolverTurbEpsilon : public PBSolverScalar
{
public:
    PBSolverTurbEpsilon(const string &name);
    virtual ~PBSolverTurbEpsilon();

    void InitMemory();
    void Discrt();
    virtual void ComputePostVisualVariables(Post_Visual *postVisualization);
    void OutputResidual(ofstream &file, int flag);

protected:
    RDouble prandtl;//! the Prandtl number of Turbulent dissipation rate.

};

DECLARE_PBSOLVER_REGISTER_FUNCTION(PBSolver, PBSolverTurbEpsilon, SOLVER_TURB_EPSILON);
}
