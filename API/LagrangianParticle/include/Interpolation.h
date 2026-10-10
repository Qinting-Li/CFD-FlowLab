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
//! @file      Interpolation.h
//! @brief     It is the class for particle interpolation.
//! @author    Lei Yinghaonan (Lanzhou University).

#pragma once
#include "TypeDefine.h"
#include "Math_BasisFunction.h"
#include "Geo_StructGrid.h"
#include "Geo_UnstructGrid.h"
#include "ParticlePointGroup.h"
#include "LagrangianVariable.h"

using namespace std;

namespace PHSPACE
{

//! The Base class of interpolation method.
class ParticleInterpolation
{
protected:
    Grid *grid;

public:
    ParticleInterpolation();
    ParticleInterpolation(Grid *gridIn);

public:
    void SetGrid(Grid *gridIn);

    virtual void ComputeNodeVar();

    virtual void ComputePointVar();

    Grid *GetGrid();
};

//! Trilinear PointInterpolation
//! At present, it is only applicable to rectangular orthogonal mesh of struct grid.
class ParticleTrilinearInterpolation : public ParticleInterpolation
{
private:
    int nPoint;
    vector<RDouble*> *pointToInterpolation;

public:
    ParticleTrilinearInterpolation();
    ParticleTrilinearInterpolation(Grid *gridIn);
    ~ParticleTrilinearInterpolation();

public:
    void SetPointToInterpolation(vector<RDouble*> *pointToInterpolationIn);
    void InitPointToInterpolation();
    void DeletePointToInterpolation();

    //! Compute node var for flow.
    void ComputeNodeVar(RDouble4D &q, RDouble4D &qNode, int &nVar);
    //! Compute node var for viscous.
    void ComputePointVar(RDouble4D &qNode, ParticlePointGroup *particlePointGroup,LDouble *lagrangianVariable ,int *varIndex,int *inIndex,int &nVar);

private:
    RDouble LinearInterpolation(RDouble &x0, RDouble &f0, RDouble &x1, RDouble &f1, RDouble &x);

};

//! The Local List for Interpolation of particle(or point) on current region.
//! Note here, PHMPI::nZones is equal to the sum of zones in all processors,
//! which get bt PHMPI::GetNumberofGlobalZones().
//! But , PHMPI::numberofZoneInThisProcess is equal to the sum of zones in current processors.
//! So this class,LocalParticleInterpolation is a local list on current region(current processors).
class LocalParticleInterpolation
{
private:
    //! The local Interpolation list.
    //! The size of localInterpolationList is equal to PHMPI::numberofZoneInThisProcess.
    //! localInterpolationList will be inited by InitMemory of Controller.
    static vector< ParticleInterpolation *> *localInterpolationList;

    //! initInterpolation: 
    //!     -- 0 Not initialized yet.
    //!     -- 1 Has been initialized.
    static vector<int> *initInterpolation;

public:
    //! Return the method of local iZone-th zone from current region zone.
    //! @param[in] iZone : the Local zone index on current region.
    static ParticleInterpolation *GetInterpolation(int iZone);

    //! If the interpolation method for each zone has been initialized,
    //! this function will return ture, else return false.
    static bool JudgeIfInitInterpolation();

    //! Init the interpolation method and push to list.
    //! @param[in] grid : The grid of iZone on current region.
    static void InitParticleInterpolation(Grid *grid);
};

}