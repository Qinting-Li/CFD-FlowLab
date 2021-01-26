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
//! @file      Force.h
//! @brief     Explain this file briefly.
//! @author    He Xin.


#pragma once
#include "DataContainer.h"
using namespace std;

namespace PHSPACE
{

class ActionKey;
class Force
{
public:
    static const int ALL_SOLID_WALL = -1;

    //! Number of force variables.
    static const int nForceVar = 16;

private:
    //! Number of solid wall part.
    int numberofSolidWallPart;

    //! Number of total part, is numberofSolidWallPart plus ONE.
    //! 'ONE' represent the total solid wall force, stored on the first position.
    int numberofTotalPart;

    //! The corresponding Boundary Condition ID of each part.
    int *partBCID;

    RDouble *cpx, *cpy, *cpz;
    RDouble *cl_tot, *cd_tot, *cd_pr, *cd_sf, *cd_cl2pa, *xcp, *ycp, *zcp, *side;
    RDouble *CA_f, *CA_p, *CN_f, *CN_p, *CZ_f, *CZ_p, *Cl_f, *Cl_p, *Cn_f, *Cn_p, *Cm_f, *Cm_p;

    RDouble *hingeMoment; 
    RDouble attack, sideslip;
    RDouble forceRefenenceArea, forceRefenenceLength, forceRefenenceLengthSpanWise;

public:
    Force();
    ~Force();

public:
    void Init(RDouble attack, RDouble sideslip, RDouble forceRefenenceArea, RDouble forceRefenenceLength, RDouble forceRefenenceLengthSpanWise);
    void CollectionForce(DataContainer *data_in);
    void CompForce();
    void DumpForce(int outnstep,std::ostringstream &oss);
    void DumpPartForce(int outnstep, int partID, const string &partName, std::ostringstream &oss);
    void DumpForce(ActionKey *actkey, int outnstep);
    void DumpData(DataContainer *data);

    //! Dump all force data in forcefile.
    void DumpAllForceData(ActionKey *actkey, int outnstep);

    //! Dump global force data.
    void DumpAllGlobalForceData(const string &partName, int partID, std::ostringstream &oss);

    //! Dump local force data.
    void DumpAllLocalForceData(const string &partName, int partID, std::ostringstream &oss);


private:
    bool isPartExist;
};

void InitForce(ActionKey *actkey);
void CollectionForce(ActionKey *actkey);
void DumpForce(ActionKey *actkey);

}
