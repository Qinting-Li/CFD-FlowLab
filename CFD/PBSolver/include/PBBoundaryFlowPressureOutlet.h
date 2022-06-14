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
//! @file      PBBoundaryFlowPressureOutlet.h
//! @brief     Boundary Conditions.
//! @author    Ming Pingjian, Sun Huawen, Song Yang.

#pragma once

#include "PBBoundaryFlow.h"

namespace PHSPACE
{
#define FLOW_PRESSURE_OUTLET FLOW_PRESSURE_OUTLET

    class PBBoundaryFlowPressureOutlet : public PBBoundaryFlow
    {
    public:
        PBBoundaryFlowPressureOutlet(const string &name);
        virtual ~PBBoundaryFlowPressureOutlet();

        virtual void CalcVelocityBoundaryCondition(UnstructGrid *grid, vector<int> *faceIndex);

        virtual void UpdatePressureBoundaryCondition(UnstructGrid *grid, vector<int> *faceIndex, Data_Param *bcData);

        virtual void CalcPressureCorrectionBoundaryCondition(UnstructGrid *grid, vector<int> *faceIndex, Data_Param *bcData);

        virtual void UpdateVelocityCorrectionBoundary(UnstructGrid *grid, vector<int> *faceIndex, Data_Param *bcData);
    };

    DECLARE_PBSOLVER_REGISTER_FUNCTION(PBBoundaryFlow, PBBoundaryFlowPressureOutlet, FLOW_PRESSURE_OUTLET);
}
