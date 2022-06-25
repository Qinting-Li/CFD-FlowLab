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
//! @file      PBBoundaryFlowVelocityInlet.h
//! @brief     Boundary Conditions of Pressure Based Solver (PBSolver).
//! @author    Wan Yunbo, Ming Pingjian, Meng Liyuan, Sun Huawen, Song Yang, Bell.

#pragma once

#include "PBBoundaryFlow.h"

namespace PHSPACE
{
#define FLOW_VELOCITY_INLET FLOW_VELOCITY_INLET

class PBBoundaryFlowVelocityInlet : public PBBoundaryFlow
{
public:
    PBBoundaryFlowVelocityInlet(const string& name);
    virtual ~PBBoundaryFlowVelocityInlet();

    void CalcVelocityBoundaryCondition(UnstructGrid *grid, vector<int> *faceIndex);

    void UpdateVelocityBoundaryCondition(UnstructGrid *grid, vector<int> *faceIndex, Data_Param *bcData);

    void UpdatePressureBoundaryCondition(UnstructGrid *grid, vector<int> *faceIndex, Data_Param *bcData);

    void CalcPressureCorrectionBoundaryCondition(UnstructGrid *grid, vector<int> *faceIndex, Data_Param *bcData);

    void UpdatePressureCorrectionBoundaryCondition(UnstructGrid *grid, vector<int> *faceIndex, Data_Param *bcData);

    void UpdateVelocityCorrectionBoundary(UnstructGrid *grid, vector<int> *faceIndex, Data_Param *bcData);
};

DECLARE_PBSOLVER_REGISTER_FUNCTION(PBBoundaryFlow, PBBoundaryFlowVelocityInlet, FLOW_VELOCITY_INLET);
}
