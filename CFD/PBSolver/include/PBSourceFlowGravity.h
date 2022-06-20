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
//! @file      PBSourceFlowGravity.h
//! @brief     Flow.
//! @author    Wan Yunbo, Ming Pingjian, Meng Liyuan, Sun Huawen, Bell.

#pragma once

#include "PBSource.h"

namespace PHSPACE
{
#define SOURCE_FLOW_GRAVITY FLOW_GRAVITY

class PBSourceFlowGravity : public PBSource
{
public:
    PBSourceFlowGravity(const string &name);
    virtual ~PBSourceFlowGravity();

    void InitMemory(UnstructGrid *grid, string phiName);
    void Discrt(UnstructGrid *grid, string phiName);

protected:
    void InitProxy(UnstructGrid *grid);

    void DiscrtProxy(UnstructGrid *grid);

    RDouble gx_;
    RDouble gy_;
    RDouble gz_;
};

DECLARE_PBSOLVER_REGISTER_FUNCTION(PBSource, PBSourceFlowGravity, SOURCE_FLOW_GRAVITY);
}
