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
//! @file      PBBoundaryScalarDefault.h
//! @brief     Boundary Conditions.
//! @author    Ming Pingjian, Sun Huawen.

#pragma once

#include "PBBoundaryScalar.h"
#include "PBBoundaryScalarConvection.h"
#include "PBBoundaryScalarDiffusion.h"

namespace PHSPACE
{
#define SCALAR_DEFAULT SCALAR_DEFAULT

    class PBBoundaryScalarDefault : public PBBoundaryScalar
    {
    public:
        PBBoundaryScalarDefault(const string &name);
        virtual ~PBBoundaryScalarDefault();

        virtual void Discrt(UnstructGrid *grid, string phiName, vector<int> *faceIndex, Data_Param *bcData, RDouble *diffCoe);

    protected:
        PBBoundaryScalarConvection *convCalc;
        PBBoundaryScalarDiffusion *diffCalc;
    };

    DECLARE_PBSOLVER_REGISTER_FUNCTION(PBBoundary, PBBoundaryScalarDefault, SCALAR_DEFAULT);
}
