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
//! @file      PBIterSolvAMG.h
//! @brief     Solving the System of Algebraic Equations in Pressure Based Solver (PBSolver).
//! @author    Wan Yunbo, Ming Pingjian, Meng Liyuan, Sun Huawen, Huang Jin, Bell.

#pragma once
#include "PBIterSolv.h"

namespace PHSPACE
{
#define ITERSOLV_AMG AMG

#ifndef _DLLPGAMG_H
#define _DLLPGAMG_H
#pragma comment(lib, "Dll1.lib")
extern "C" {void _stdcall PGAMG(int *nrows, int *ncols, int *nnonzeros, int *IA, int *JA, double *d, double *A, double *u, double *b, int *groupsize, int *maxSweep, double *eps, int *sweep, double *res0, double *res1); }
#endif

class PBIterSolvAMG : public PBIterSolv
{
public:
    PBIterSolvAMG(const string &name);
    ~PBIterSolvAMG();

    void CSRsp(RDouble *ajLeft, RDouble *ajRight, int ncells, int nfaces, int *leftCellOfFace, int *rightCellOfFace, 
        int *faceNumberOfEachCell, int **faceIndexOfEachCell, int boundFaces);

    void Run(Region &region, vector<int> &zoneIdList, string varStr, string apStr, string ajUpperStr, string ajLowerStr,
        string ajbUpperStr, string ajbLowerStr, string scStr, string rsStr, int maxSweep = 50, RDouble tol = 1e-4);

protected:
    void Calc(int nfaces, int ncells, int *leftCellOfFace, int *rightCellOfFace, RDouble *ap, RDouble *ajLeft, RDouble *ajRight, RDouble *sc, RDouble *rs,
        string varStr, int maxSweep = 50, RDouble tol = 1e-4);

private:
    int entryNum;
    int groupsize;
    int *IA;
    int *JA;
    int *faceNumberOfEachCell;
    int **faceIndexOfEachCell;
    int boundFaces;
    RDouble *AMG_A;
};

DECLARE_PBSOLVER_REGISTER_FUNCTION(PBIterSolv, PBIterSolvAMG, ITERSOLV_AMG);
}