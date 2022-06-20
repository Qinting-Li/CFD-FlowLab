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
//! @file      PBPreconditioner.h
//! @brief     Preconditioner.
//! @author    Wan Yunbo, Ming Pingjian, Meng Liyuan, Sun Huawen, Bell.
 
#pragma once

#include "PBObject.h"
#include "Region.h"

namespace PHSPACE
{
#define ITERSOLVPRECOND_EMPTY EMPTY

class PBPreconditioner : public PBObject
{
public:
    PBPreconditioner(const string &name);
    virtual ~PBPreconditioner();

    virtual void Preconditioning(Region &region, vector<int> &zoneId, string &varStr, string &apStr, string &ajUpperStr,
        string &ajLowerStr, string &aiiStr, string &scStr, string &rsStr);

public: 
    virtual void PreconditioningProxy(int nFaces, int nFells, int *leftCellOfFace, int *rightCellOfFace,
        RDouble *ap, RDouble *ajLeft, RDouble *ajRight, RDouble *rs, string *phiName = NULL);

    virtual void ForwardAndBackward(int nFaces, int nCells, int *leftCellOfFace, int *rightCellOfFace, RDouble *ap0,
        RDouble *ajLeft, RDouble *ajRight, RDouble *p, RDouble *rs, string *phiName = NULL);

    UnstructGrid *grid;
};

DEFINE_PBOBJECT_FACTORY(PBPreconditioner);

DECLARE_PBSOLVER_REGISTER_FUNCTION(PBPreconditioner, PBPreconditioner, ITERSOLVPRECOND_EMPTY);
}