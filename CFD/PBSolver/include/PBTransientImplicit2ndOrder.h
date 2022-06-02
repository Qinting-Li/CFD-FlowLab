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
//! @file      PBTransientImplicit2ndOrder.h
//! @brief     Implicit2ndOrder.
//! @author    Ming Pingjian, Sun Huawen.

#pragma once

#include "PBTransient.h"

namespace PHSPACE
{
#define TRANSIENT_IMPLICIT_2ND_ORDER IMPLICIT_2ND_ORDER

    class PBTransientImplicit2ndOrder : public PBTransient
    {
    public:
        PBTransientImplicit2ndOrder(const string &name);
        virtual ~PBTransientImplicit2ndOrder();

        virtual void AllocateMemory(UnstructGrid *grid, vector<string> &phiName, vector<RDouble> &initPhi);

        virtual void Discrt(UnstructGrid *grid, string phiName, string phiOldName, string phiOldOldName, string rhoOldName = NULL, string rhoOldOldName = NULL);

        virtual void SaveOldTimeValue(UnstructGrid *grid, vector<string> &phiName);
    };

    DECLARE_PBSOLVER_REGISTER_FUNCTION(PBTransient, PBTransientImplicit2ndOrder, TRANSIENT_IMPLICIT_2ND_ORDER);
}
