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
//! @file      PBBoundary.h
//! @brief     Boundary Conditions of Pressure Based Solver (PBSolver).
//! @author    Wan Yunbo, Ming Pingjian, Meng Liyuan, Sun Huawen, Bell.
 
#pragma once

#include "PBObject.h"
#include "Geo_UnstructGrid.h"
#include "Precision.h"

namespace PHSPACE
{
#define BOUNDARY_EMPTY EMPTY

class Data_Param;

class PBBoundary : public PBObject
{
public:
    PBBoundary(const string &name);
    virtual ~PBBoundary();

public:
    virtual void Init(string phiName);

    virtual void Discrt(UnstructGrid *grid, string phiName, vector<int> *faceIndex, Data_Param *bcData, RDouble *diffCoe);

    virtual void Update(UnstructGrid *grid, string phiName, vector<int> *faceIndex, Data_Param *bcData);

    virtual void WallModification(UnstructGrid *grid, string phiName, vector<int> *faceIndex, Data_Param *bcData);

    virtual void UpdateVistb(UnstructGrid *grid, string phiName, vector<int> *faceIndex, Data_Param *bcData);
};

DEFINE_PBOBJECT_FACTORY(PBBoundary);

DECLARE_PBSOLVER_REGISTER_FUNCTION(PBBoundary, PBBoundary, BOUNDARY_EMPTY);
}
