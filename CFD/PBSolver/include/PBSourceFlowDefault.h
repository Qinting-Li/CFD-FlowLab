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
//! @file      PBSourceFlowDefault.h
//! @brief     Flow of Pressure Based Solver (PBSolver).
//! @author    Wan Yunbo, Ming Pingjian, Meng Liyuan, Sun Huawen, Bell.

#pragma once

#include "PBObject.h"
#include "PBSource.h"

namespace PHSPACE
{
#define SOURCE_FLOW_DEFAULT FLOW_DEFAULT

class PBSourceFlowDefault : public PBSource
{
public:
    PBSourceFlowDefault(const string &name);
    virtual ~PBSourceFlowDefault();

    void Discrt(UnstructGrid *grid, string phiName);

protected:
    void DiscrtProxy(UnstructGrid *grid);
};

DECLARE_PBSOLVER_REGISTER_FUNCTION(PBSource, PBSourceFlowDefault, SOURCE_FLOW_DEFAULT);
}
