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
//! @file      PBBoundaryEnergyDefault.h
//! @brief     Boundary Conditions.
//! @author    Ming Pingjian, Sun Huawen, Huang Jin.

#pragma once

#include "PBBoundaryScalarDefault.h"

namespace PHSPACE
{
#define ENERGY_DEFAULT ENERGY_DEFAULT 

    class PBBoundaryEnergyDefault : public PBBoundaryScalarDefault
    {
    public:
        PBBoundaryEnergyDefault(const string &name);
        virtual ~PBBoundaryEnergyDefault();
    };

    DECLARE_PBSOLVER_REGISTER_FUNCTION(PBBoundary, PBBoundaryEnergyDefault, ENERGY_DEFAULT);
}
