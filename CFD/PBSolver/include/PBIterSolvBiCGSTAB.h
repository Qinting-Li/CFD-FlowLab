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
//! @file      PBIterSolvBiCGSTAB.h
//! @brief     Solving the System of Algebraic Equations.
//! @author    Ming Pingjian, Sun Huawen, Huang Jin.
 
#pragma once
#include "PBIterSolv.h"

namespace PHSPACE
{
#define ITERSOLV_BiCGSTAB BiCGSTAB

    class PBIterSolvBiCGSTAB : public PBIterSolv
    {
    public:
        PBIterSolvBiCGSTAB(const string &name);
        ~PBIterSolvBiCGSTAB();

        virtual void Run(Region &region, vector<int> &zoneIdList, string varStr, string apStr, string ajUpperStr, string ajLowerStr,
            string ajbUpperStr, string ajbLowerStr, string scStr, string rsStr, int maxSweep = 50, RDouble tol = 1e-4);

    protected:
        virtual void Calc(int nFaces, int nCells, int *leftCellOfFace, int *rightCellOfFace, RDouble *ap, RDouble *ajLeft, RDouble *ajRight, RDouble *sc,
            RDouble *rs, string varStr, int maxSweep = 50, RDouble tol = 1e-4);
    };

    DECLARE_PBSOLVER_REGISTER_FUNCTION(PBIterSolv, PBIterSolvBiCGSTAB, ITERSOLV_BiCGSTAB);
}