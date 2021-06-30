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
//! @file      Flux_RoeEntropyCorrection.h
//! @brief     This file defines the entropy correction (fix) method in Roe inviscid flux method.
//! @author    Bell, Liu Jian.

#pragma once
#include "LIB_Macro.h"
#include "Math_BasisFunction.h"

namespace PHSPACE
{
//! The entropy correction (fix) in Roe scheme of inviscid flux, which is used to limit the magnitude of the three eigenvalue.\n
//!    @param[in] eigv1,2,3                    The three eigenvalue, which are V/V+c/V-c on the face. They will be corrected and returned.\n
//!    @param[in] vn                           The normal velocity on the face.\n
//!    @param[in] RoeEntropyCorrectionMethod   The Roe entropy correction method.\n
//!    @param[in] entropyCorrectionCoef        The Roe entropy correction coefficient.\n
LIB_EXPORT void Flux_RoeEntropyCorrection(RDouble &eigv1, RDouble &eigv2, RDouble &eigv3, RDouble lrtL, RDouble lrtR, 
                                          RDouble vn, RDouble cm, int RoeEntropyCorrectionMethod, RDouble acousticEigenLimit, RDouble convectEigenLimit);

//! Harten type of entropy correction (fix).\n
//!    @param[in] eigenValue        The original eigenvalue, which may be invalid.
//!    @param[in] eigenValueLimit   The epsilon in the Harten formula.
//!    @return                      The new eigenvalue which has been limited. 
inline RDouble Harten(RDouble eigenValue, RDouble eigenValueLimit)
{
    RDouble eigenValueNew = ABS(eigenValue);
    if(eigenValueNew < eigenValueLimit) 
    {
        eigenValueNew = 0.5 * (eigenValue * eigenValue + eigenValueLimit * eigenValueLimit) / eigenValueLimit;
    }
    return eigenValueNew;
}

}