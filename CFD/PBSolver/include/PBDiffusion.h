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
//! @file      PBDiffusion.h
//! @brief     The Diffusion Term.
//! @author    Wan Yunbo, Ming Pingjian, Meng Liyuan, Sun Huawen, Bell.
 
#pragma once

#include "PBObject.h"

namespace PHSPACE
{
#define DIFFUSION_EMPTY EMPTY
class UnstructGrid;
class PBDiffusion : public PBObject
{
public:
    PBDiffusion(const string &name);
    virtual ~PBDiffusion();

public:
    virtual void Discrt(UnstructGrid *grid, string phiName, RDouble *diffCoe, RDouble *dphidx, RDouble *dphidy, RDouble *dphidz,
        RDouble *ajLeft, RDouble *ajRight, RDouble *scPhi, RDouble *dvdx = NULL, RDouble *dvdy = NULL, RDouble *dvdz = NULL,
        RDouble *scV = NULL, RDouble *dwdx = NULL, RDouble *dwdy = NULL, RDouble *dwdz = NULL, RDouble *scW = NULL);
};

DEFINE_PBOBJECT_FACTORY(PBDiffusion);

DECLARE_PBSOLVER_REGISTER_FUNCTION(PBDiffusion, PBDiffusion, DIFFUSION_EMPTY);
}
