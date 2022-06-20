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
//! @file      PBBoundaryTurbEpsilonWall.h
//! @brief     Boundary Conditions.
//! @author    Wan Yunbo, Ming Pingjian, Meng Liyuan, Sun Huawen, Song Yang, Bell.

#pragma once

#include "Geo_UnstructGrid.h"
#include "Data_Param.h"
#include "PBBoundaryScalar.h"

namespace PHSPACE
{
#define TURB_EPSILON_WALL TURB_EPSILON_WALL

class PBBoundaryTurbEpsilonWall : public PBBoundaryScalar
{
public:
    PBBoundaryTurbEpsilonWall(const string &name);
    virtual ~PBBoundaryTurbEpsilonWall();

    void Discrt(UnstructGrid *grid, string phiName, vector<int> *faceIndex, Data_Param *bcData, RDouble *diffCoe);
    void Update(UnstructGrid *grid, string phiName, vector<int> *faceIndex, Data_Param *bcData);

protected:
    RDouble cmu;
    RDouble cmu25_;
    RDouble cmu75_;
    RDouble cappa_;
    RDouble yptr_;
};

DECLARE_PBSOLVER_REGISTER_FUNCTION(PBBoundary, PBBoundaryTurbEpsilonWall, TURB_EPSILON_WALL);
} // PB_SPACE
