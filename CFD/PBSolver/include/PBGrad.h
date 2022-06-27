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
//! @file      PBGrad.h
//! @brief     Compute Grad of Pressure Based Solver (PBSolver).
//! @author    Wan Yunbo, Ming Pingjian, Meng Liyuan, Sun Huawen, Bell.
 
#pragma once

#include "PBObject.h"

namespace PHSPACE
{
#define GRAD_EMPTY EMPTY

class UnstructGrid;

class PBGrad : public PBObject
{
public:
    PBGrad(const string &name);
    virtual ~PBGrad();

public:
    virtual void Calc(UnstructGrid *grid, string phiName, string dphidxName, string dphidyName, string dphidzName, bool isComm = false);

    void SetSolverID(int solverID) { this->solverID = solverID; }
    int GetSolverID() { return this->solverID; }

protected:
    UnstructGrid *grid;
    int solverID;
};

DEFINE_PBOBJECT_FACTORY(PBGrad);

DECLARE_PBSOLVER_REGISTER_FUNCTION(PBGrad, PBGrad, GRAD_EMPTY);
}
