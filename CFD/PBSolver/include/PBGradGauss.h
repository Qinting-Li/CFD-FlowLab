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
//! @file      PBGradGauss.h
//! @brief     Gauss of Pressure Based Solver (PBSolver).
//! @author    Wan Yunbo, Ming Pingjian, Meng Liyuan, Sun Huawen, Song Yang, Bell.

#pragma once

#include "PBGrad.h"

namespace PHSPACE
{
#define GRAD_GAUSS GAUSS

class PBGradGauss : public PBGrad
{
public:
    PBGradGauss(const string &name);
    virtual ~PBGradGauss();

    void Calc(UnstructGrid *grid, string phiName, string dphidxName, string dphidyName, string dphidzName, bool isComm = false);
};

DECLARE_PBSOLVER_REGISTER_FUNCTION(PBGrad, PBGradGauss, GRAD_GAUSS);
}
