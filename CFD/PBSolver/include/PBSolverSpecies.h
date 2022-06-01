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
//! @file      PBSolverSpecies.h
//! @brief     Species.
//! @author    Huang Jin.

#pragma once

#include "PBSolverScalar.h"

namespace PHSPACE
{
#define SOLVER_SPECIES SPECIES

    class PBSolverSpecies : public PBSolverScalar
    {
    public:
        PBSolverSpecies(const std::string& name);
        ~PBSolverSpecies();

    public:
        virtual void Init();
        virtual void AllocateMemory();
        virtual void ReinitField();
        virtual void UpdateDt();
        virtual bool SolutionIsConverged();

    public:
        virtual void UpdateProperties();
        virtual void UpdateAfterIterloop();

    protected:
        int isSolveTurb;
        std::string* speciesName; /// 控制面板输出变量控制；
        RDouble* speciesMolarMass; /// 各组分摩尔质量；
        RDouble* speciesRho; /// 各组分密度；
        RDouble* speciesMu; /// 各组分粘性系数；

        RDouble dilute; /// 湍流开关，用于组分方程扩散系数的求解；
        RDouble sct; /// Schmidt数，用于求组分方程扩散系数；
        RDouble gasR; /// 理想气体常数;
        RDouble speciesThreshold;
        RDouble initMuSp;
        RDouble* muSp;
        RDouble* muSpb;

    };

    DECLARE_PBSOLVER_REGISTER_FUNCTION(PBSolver, PBSolverSpecies, SOLVER_SPECIES);
}