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
//! @file      PBOutputTecplot.h
//! @brief     Tecplot.
//! @author    Ming Pingjian, Sun Huawen.
 
#pragma once

#include <vector>

#include "PBOutput.h"

namespace PHSPACE
{
#define OUTPUT_TECPLOT TECPLOT

    class PBOutputTecplot : public PBOutput
    {
    public:
        PBOutputTecplot(const string &name);
        virtual ~PBOutputTecplot();

        virtual void Run(Zone *zone, vector<string> &varNameList);
    };

    DECLARE_PBSOLVER_REGISTER_FUNCTION(PBOutput, PBOutputTecplot, OUTPUT_TECPLOT);
}
