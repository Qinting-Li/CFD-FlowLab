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
//! @file      PBGradLsq.h
//! @brief     Lsq.
//! @author    Wan Yunbo, Ming Pingjian, Meng Liyuan, Sun Huawen, Song Yang, Bell.
 
#pragma once

#include "PBGrad.h"

namespace PHSPACE
{
#define GRAD_LSQ LSQ

class PBGradLsq : public PBGrad
{
public:
    PBGradLsq(const string &name);
    virtual ~PBGradLsq();

    void Calc(UnstructGrid *grid, string phiName, string dphidxName, string dphidyName, string dphidzName, bool isComm = false);
};

DECLARE_PBSOLVER_REGISTER_FUNCTION(PBGrad, PBGradLsq, GRAD_LSQ);
}
