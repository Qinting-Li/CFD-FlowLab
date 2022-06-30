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
//! @file      PBSourceTurbKDefault.h
//! @brief     TurbK of Pressure Based Solver (PBSolver).
//! @author    Wan Yunbo, Ming Pingjian, Meng Liyuan, Sun Huawen, Song Yang, Bell.

#pragma once

#include "PBSource.h"
#include "PBObject.h"

namespace PHSPACE
{
#define SOURCE_TURB_K_DEFAULT TURB_K_DEFAULT

class PBSourceTurbKDefault : public PBSource
{
public:
    PBSourceTurbKDefault(const string &name);
    virtual ~PBSourceTurbKDefault();

    void InitMemory(UnstructGrid *grid, string phiName);
    void Discrt(UnstructGrid *grid, string phiName);
};

DECLARE_PBSOLVER_REGISTER_FUNCTION(PBSource, PBSourceTurbKDefault, SOURCE_TURB_K_DEFAULT);
}
