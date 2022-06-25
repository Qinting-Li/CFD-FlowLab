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
//! @file      PBOutput.h
//! @brief     Post-processing output of Pressure Based Solver (PBSolver).
//! @author    Wan Yunbo, Ming Pingjian, Meng Liyuan, Sun Huawen, Bell.
 
#pragma once

#include "PBObject.h"

namespace PHSPACE
{
#define OUTPUT_EMPTY EMPTY

class Zone;

class PBOutput : public PBObject
{
public:
    PBOutput(const string &name);
    virtual ~PBOutput();

public:
    virtual void Run(Zone *zone, vector<string> &varNameList);
};

DEFINE_PBOBJECT_FACTORY(PBOutput);

DECLARE_PBSOLVER_REGISTER_FUNCTION(PBOutput, PBOutput, OUTPUT_EMPTY);
}
