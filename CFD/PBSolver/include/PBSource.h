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
//! @file      PBSource.h
//! @brief     Discretization of the Source Term in Pressure Based Solver (PBSolver).
//! @author    Wan Yunbo, Ming Pingjian, Meng Liyuan, Sun Huawen, Bell.

#pragma once

#include "PBObject.h"
#include "Geo_UnstructGrid.h"

namespace PHSPACE
{
#define SOURCE_EMPTY EMPTY

class PBSource : public PBObject
{
public:
    PBSource(const string &name);
    virtual ~PBSource();

public:
    virtual void InitMemory(UnstructGrid *grid, string phiName);

    virtual void Discrt(UnstructGrid *grid, string phiName);

    virtual void UpdateBodyForce(UnstructGrid *grid, string phiName);
};

DEFINE_PBOBJECT_FACTORY(PBSource);

DECLARE_PBSOLVER_REGISTER_FUNCTION(PBSource, PBSource, SOURCE_EMPTY);
}
