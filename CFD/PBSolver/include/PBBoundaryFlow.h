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
//! @file      PBBoundaryFlow.h
//! @brief     Boundary Conditions of Pressure Based Solver (PBSolver).
//! @author    Wan Yunbo, Ming Pingjian, Meng Liyuan, Sun Huawen, Song Yang, Bell.

#pragma once

#include "PBBoundary.h"

namespace PHSPACE
{
class Data_Param;

class PBBoundaryFlow : public PBBoundary
{
public:
    PBBoundaryFlow(const string &name);
    virtual ~PBBoundaryFlow();

    //! @brief Calculate the velocity boundary conditions.
    //! @param[in]: gridIn denotes the current computational regions of the grids.
    //! @param[in]: faceIndex denotes the faceIndex set storing the all face numbers in this boundary region.
    virtual void CalcVelocityBoundaryCondition(UnstructGrid *grid, vector<int> *faceIndex);

    //! @brief Update the velocity boundary conditions.
    //! @param[in]: gridIn denotes the current computational regions of the grids.
    //! @param[in]: faceIndex denotes the faceIndex set storing the all face numbers in this boundary region.
    //! @param[in]: bcData denotes the information of boundaries.
    virtual void UpdateVelocityBoundaryCondition(UnstructGrid *grid, vector<int> *faceIndex, Data_Param *bcData);

    //! @brief Update the pressure boundary conditions.
    //! @param[in]: gridIn denotes the current computational regions of the grids.
    //! @param[in]: faceIndex denotes the faceIndex set storing the all face numbers in this boundary region.
    //! @param[in]: bcData denotes the information of boundaries.
    virtual void UpdatePressureBoundaryCondition(UnstructGrid *grid, vector<int> *faceIndex, Data_Param *bcData);

    //! @brief Calculate the pressure boundary conditions.
    //! @param[in]: gridIn denotes the current computational regions of the grids.
    //! @param[in]: faceIndex denotes the faceIndex set storing the all face numbers in this boundary region.
    //! @param[in]: bcData denotes the information of boundaries.
    virtual void CalcPressureCorrectionBoundaryCondition(UnstructGrid *grid, vector<int> *faceIndex, Data_Param *bcData);



    //! @brief Update the pressure correction boundary conditions.
    //! @param[in]: gridIn denotes the current computational regions of the grids.
    //! @param[in]: faceIndex denotes the faceIndex set storing the all face numbers in this boundary region.
    //! @param[in]: bcData denotes the information of boundaries.
    virtual void UpdatePressureCorrectionBoundaryCondition(UnstructGrid *grid, vector<int> *faceIndex, Data_Param *bcData);

    //! @brief Update the velocity correction boundary conditions.
    //! @param[in]: gridIn denotes the current computational regions of the grids.
    //! @param[in]: faceIndex denotes the faceIndex set storing the all face numbers in this boundary region.
    //! @param[in]: bcData denotes the information of boundaries.
    virtual void UpdateVelocityCorrectionBoundary(UnstructGrid *grid, vector<int> *faceIndex, Data_Param *bcData);
};

DEFINE_PBOBJECT_FACTORY(PBBoundaryFlow);
}
