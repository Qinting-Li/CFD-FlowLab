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
//! @file      PBBoundaryTurbEpsilonInlet.h
//! @brief     Boundary Conditions.
//! @author    Wan Yunbo, Ming Pingjian, Meng Liyuan, Sun Huawen, Song Yang, Bell.

#pragma once

#include "Geo_UnstructGrid.h"
#include "Data_Param.h"
#include "PBBoundaryTurbDefault.h"

namespace PHSPACE
{
#define TURB_EPSILON_INLET TURB_EPSILON_INLET 

class PBBoundaryTurbEpsilonInlet : public PBBoundaryTurbDefault
{
public:
    PBBoundaryTurbEpsilonInlet (const string &name);
    virtual ~PBBoundaryTurbEpsilonInlet ();

    void Update(UnstructGrid *grid, string phiName, vector<int> *faceIndex, Data_Param *bcData);

    void Discrt(UnstructGrid *grid, string phiName, vector<int> *faceIndex, Data_Param *bcData, RDouble *diffCoe);

protected:
    RDouble cmu;
};


DECLARE_PBSOLVER_REGISTER_FUNCTION(PBBoundary, PBBoundaryTurbEpsilonInlet, TURB_EPSILON_INLET);
} // PB_SPACE
