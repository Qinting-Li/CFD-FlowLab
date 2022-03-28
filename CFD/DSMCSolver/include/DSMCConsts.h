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
//! @file      DSMCConsts.h
//! @brief     Consts for DSMC.
//! @author    Li Jin.

#pragma once
#include "Precision.h"
using namespace std;

namespace PHSPACE
{
    static const int MaxMoleculeNumber = 1000000;
    static const int BufferLen = 500;
    static const int MaxCPUNum = 9999;
    static const int MaxPartnerNum = 20000;
    static const int MaxSurfNum = 6;
    static const int MaxPntsNum = 8;

    static const RDouble BOLTZ = 1.380658e-23;
    static const RDouble fMOV_INT_BUFF_factor = 10.0;
    static const RDouble ZERO_V_ParticleInCellChk = 1.0e-7;
    static const RDouble ZERO_XYZC_SlightMove = 2.0e-3;
    static const RDouble ZERO_ParticleInSurf = 9.0e-4;
    static const RDouble ZERO_StPara_Start = 5.0e-6;
    static const RDouble ZERO_StPara_End = 5.0e-4;
    static const RDouble DOUB_XYZC_SlightMove = 3.0;
    static const RDouble ZERO_ENTER2_SlightMove = 1.0e-5;

    static const int METER      = 0;
    static const int MILLIMETER = 1;
}