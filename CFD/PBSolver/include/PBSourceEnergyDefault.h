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
//! @file      PBSourceEnergyDefault.h
//! @brief     Energy.
//! @author    Wan Yunbo, Ming Pingjian, Meng Liyuan, Sun Huawen, Huang Jin, Bell.

#pragma once

#include "PBObject.h"
#include "PBSource.h"

namespace PHSPACE
{
#define SOURCE_ENERGY_DEFAULT SOURCE_ENERGY_DEFAULT

class PBSourceEnergyDefault : public PBSource
{
public:
    PBSourceEnergyDefault(const string &name);
    virtual ~PBSourceEnergyDefault();

    void Discrt(UnstructGrid *grid, string &phiName);

protected:
    void DiscrtProxy(UnstructGrid *grid);
};

DECLARE_PBSOLVER_REGISTER_FUNCTION(PBSource, PBSourceEnergyDefault, SOURCE_ENERGY_DEFAULT);
}

