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
//! @file      ParticleBoundaryCondition.h
//! @brief     It is the class for the Boundary Condition of particle.
//! @author    Lei Yinghaonan (Lanzhou University).

using namespace std;

#pragma once
#include "Precision.h"
#include "TypeDefine.h"
#include "Constants.h"
#include "Data_Param.h"

namespace PHSPACE
{
namespace PARTICLEBCSPACE
{
//! Interface as flow.
const int INTERFACE = -1;
//! No partcile boundary.
const int NO_PARTICLE_BC = 0;
//! Wall
const int FULLY_ELASTIC_COLLISION_WALL = 1;
//! For out particle remove.
const int REMOVE_PARTICLE = 2;
//! For out particle remove but reset as initFile.
const int REMOVE_PARTICLE_AS_INIT = 3;
//! For out particle remove but reset as initFile only coordinate (but others don't change).
const int REMOVE_PARTICLE_AS_INIT_ONLY_COORD = 4;
//! Periodic boundary condition (with link). (To do later).
const int PERIODIC = 5;
};

class GlobalParticleBoundaryLink
{
private:
    //! a list of link for particle boundary condition id.
    static vector< vector<string>* >  *globalLinkName;

public:
    GlobalParticleBoundaryLink();

    ~GlobalParticleBoundaryLink();
public:
    static void InitGlobalParticleBoundaryLink(int nBC);

    static void AddBCLinkName(int iBC,string bcLinkName);

    static int GetNumOfGlobalParticleBoundaryLink();

    static int GetNumOfLinkBCByID(int iBC);

    static bool IsBCLink(int iBC);

    static vector< vector<string>* > * GetGlobalLinkName();

};

class ParticleBoundaryCondition
{
private:
    //! BC name ,which should equal flow bc name.
    string particleBCName;

    //! Type of particle boundary condition in PARTICLEBCSPACE for each bc.
    int particleBCType;

    //! param
    Data_Param *particleBCParamDataBase;

    //! particle boundary condition link.
    vector<ParticleBoundaryCondition* > *particleBCLink;

public:
    ParticleBoundaryCondition();

    ~ParticleBoundaryCondition();

public:
    void SetParticleBCName(string &particleBCName);
    void SetParticleBCTye(int &particleBCType);
    void SetParticleBCParam(Data_Param* particleBCParamDataBase);

    const string &GetParticleBCName();
    int GetParticleBCType();
    Data_Param *GetParticleBCParam();

    bool IsPeriodicBC();
};



}