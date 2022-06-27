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
//! @file      PBSolverEnergy.h
//! @brief     Energy of Pressure Based Solver (PBSolver).
//! @author    Wan Yunbo, Ming Pingjian, Meng Liyuan, Sun Huawen, Huang Jin, Bell..

#pragma once

#include "PBSolverScalar.h"

namespace PHSPACE
{
#define SOLVER_ENERGY ENERGY

class PBSolverEnergy : public PBSolverScalar
{
public:
    PBSolverEnergy(const string &name);
    ~PBSolverEnergy();

public:
    void InitMemory();
    void AllocateMemory();
    void IncompressibleInitial();
    void CalcDiffusionTerm();
    void PrimeEquations();
    void Relax();
    void Update();
    void UpdateProperties();
    void SetOutputs(vector<string> &phiNameList);
    bool SolutionIsConverged();
    void ComputePostVisualVariables(Post_Visual *postVisualization);
    void OutputResidual(ofstream &file, int flag);

protected:
    RDouble initT;
    RDouble initH;
    RDouble urfT;
    RDouble urfH;
    RDouble initk;
    RDouble initCPg;
    RDouble tMax;
    RDouble tMin;
    RDouble *cpg;
    RDouble *cpgb;
    RDouble *k;
    RDouble *kb;
    RDouble *muT;
    RDouble *muTb;
    RDouble *h;
    RDouble *hb;
    RDouble *t;
    RDouble *tb;
    std::string *energyPrintName;
};

DECLARE_PBSOLVER_REGISTER_FUNCTION(PBSolver, PBSolverEnergy, SOLVER_ENERGY);
}