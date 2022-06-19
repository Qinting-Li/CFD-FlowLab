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
//! @file      PBIterSolv.h
//! @brief     Solving the System of Algebraic Equations.
//! @author    Wan Yunbo, Ming Pingjian, Meng Liyuan, Sun Huawen, Bell.
 
#pragma once

#include "PBObject.h"
#include "Region.h"
#include "Zone.h"
#include "Geometry.h"
#include <vector>

namespace PHSPACE
{
#define ITERSOLV_EMPTY EMPTY

class PBPreconditioner;

class PBIterSolv : public PBObject
{
public:
    PBIterSolv(const string &name);
    virtual ~PBIterSolv();

public:
    virtual void Run(Region &region, vector<int> &zoneIdList, string varStr, string apStr, string ajUpperStr, string ajLowerStr,
        string ajbUpperStr, string ajbLowerStr, string scStr, string rsStr, int maxSweep = 50, RDouble tol = 1e-4);

protected:
    void Matrix(int nFaces, int nCells, int *leftCellOfFace, int *rightCellOfFace, RDouble *ap,RDouble *ajLeft,
        RDouble *ajRight, RDouble *p, RDouble *rs, string *phiName = NULL);

    RDouble Dot(RDouble *b, RDouble *c, int rows);

    void Rotg(RDouble *X, RDouble *Y, RDouble *c, RDouble *s);

    void Rot(RDouble *X, RDouble *Y, RDouble *c, RDouble *s);

    PBPreconditioner *preconditioner_;
    UnstructGrid *grid;
};

DEFINE_PBOBJECT_FACTORY(PBIterSolv);

DECLARE_PBSOLVER_REGISTER_FUNCTION(PBIterSolv, PBIterSolv, ITERSOLV_EMPTY);
}
