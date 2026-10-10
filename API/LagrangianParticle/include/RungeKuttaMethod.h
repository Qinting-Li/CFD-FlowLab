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
//! @file      RungeKuttaMethod.h
//! @brief     The Particle point solver.
//! @author    Lei Yinghaonan (Lanzhou University).

#pragma once
#include "Precision.h"
namespace PHSPACE
{
class RungeKuttaIndex
{
private:
    RDouble *kindex;
    RDouble *findex;
    int index;
    int nRK;
public:
    RungeKuttaIndex(int nRK);
    ~RungeKuttaIndex();
    void AddIndex();

    int GetIndex();
    int GetNumRK();

    RDouble GetFindex();
    RDouble GetKindex();

    void ResetIndex();
};
}