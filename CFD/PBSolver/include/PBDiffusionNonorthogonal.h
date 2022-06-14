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
//! @file      PBDiffusionNonorthogonal.h
//! @brief     The Diffusion Term.
//! @author    Ming Pingjian, Sun Huawen.
 
#pragma once

#include "PBDiffusion.h"
#include "Geo_UnstructGrid.h"

namespace PHSPACE
{
#define DIFFUSION_NON_ORTHOGONAL NON_ORTHOGONAL

    class UnstructGrid;

    class PBDiffusionNonorthogonal : public PBDiffusion
    {
    public:
        PBDiffusionNonorthogonal(const string &name);
        virtual ~PBDiffusionNonorthogonal();

        virtual void Discrt(UnstructGrid *grid, string phiName, RDouble *diffCoe, RDouble *dphidx, RDouble *dphidy, RDouble *dphidz, RDouble *ajLeft,
            RDouble *ajRight, RDouble *scPhi, RDouble *dvdx = NULL, RDouble *dvdy = NULL, RDouble *dvdz = NULL,
            RDouble *scV = NULL, RDouble *dwdx = NULL, RDouble *dwdy = NULL, RDouble *dwdz = NULL, RDouble *scW = NULL);
    };

    DECLARE_PBSOLVER_REGISTER_FUNCTION(PBDiffusion, PBDiffusionNonorthogonal, DIFFUSION_NON_ORTHOGONAL);
}
