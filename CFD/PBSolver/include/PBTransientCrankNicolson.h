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
//! @file      PBTransientCrankNicolson.h
//! @brief     CrankNicolson.
//! @author    Wan Yunbo, Ming Pingjian, Meng Liyuan, Sun Huawen, Bell.

#pragma once

#include "PBTransient.h"

namespace PHSPACE
{
#define TRANSIENT_CRANK_NICOLSON CRANK_NICOLSON

class PBTransientCrankNicolson : public PBTransient
{
public:
    PBTransientCrankNicolson(const string &name);
    virtual ~PBTransientCrankNicolson();

    void AllocateMemory(UnstructGrid *grid, vector<string> &phiName, vector<RDouble> &initPhi);

    void Discrt(UnstructGrid *grid, string phiName, string phiOldName, string phiOldOldName, string rhoOldName = NULL, string rhoOldOldName = NULL);

    void SaveOldTimeValue(UnstructGrid *grid, vector<string> &phiName);

};

DECLARE_PBSOLVER_REGISTER_FUNCTION(PBTransient, PBTransientCrankNicolson, TRANSIENT_CRANK_NICOLSON);
}
