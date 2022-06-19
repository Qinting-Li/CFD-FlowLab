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
//! @file      PBBoundaryEnergyInlet.h
//! @brief     Boundary Conditions.
//! @author    Wan Yunbo, Ming Pingjian, Meng Liyuan, Sun Huawen, Huang Jin, Bell.

#pragma once

#include "PBBoundaryEnergyDefault.h"

namespace PHSPACE
{
#define ENERGY_INLET ENERGY_INLET 

class PBBoundaryEnergyInlet : public PBBoundaryEnergyDefault
{
public:
    PBBoundaryEnergyInlet(const string &name);
    virtual ~PBBoundaryEnergyInlet();

    void Discrt(UnstructGrid *grid, string phiName, vector<int> *faceIndex, Data_Param *bcData, RDouble *diffCoe);

    void Update(UnstructGrid *grid, string phiName, vector<int> *faceIndex, Data_Param *bcData);

};


DECLARE_PBSOLVER_REGISTER_FUNCTION(PBBoundary, PBBoundaryEnergyInlet, ENERGY_INLET);
}