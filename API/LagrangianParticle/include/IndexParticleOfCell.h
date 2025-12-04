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
//! @file      IndexParticleOfCell.h
//! @brief     The index of particle in one grid cell.
//! @author    Lei Yinghaonan (Lanzhou University).

#pragma once
#include "TypeDefine.h"
#include <iostream>
#include <fstream>

namespace PHSPACE
{

namespace INDEX_CELLTYPE
{

const int NO_PAR_IN = -1;

const int NO_BC = 0;

const int GHOST_IN = 1;
const int CORNER_IN = 2;

const int GHOST_OUT = 3;
const int CORNER_OUT = 4;

const int CENTER_ZONE = 5;

}

class IndexParticleOfCell
{
private:
    //! The set of particle id in one grid cell.
    set<int> *particleCellID;

    //! cellType : location of cell.
    //!            -1 -- the particle will not in this cell.
    //!            0 -- not init.
    //!            1 -- ghost inner
    //!                 cell on the sides inside the zone 
    //!                (except for the inner corner)
    //!            2 -- corner inner :
    //!                 cell in the corner inside the zone.
    //!            3 -- ghost ourner :
    //!                 cell in the ghost area at ournner
    //!                 the boundary of the zone.
    //!            4 -- corner ourner
    //!                 cell in corner outside the zone.
    //!            5 -- cell in zone (except for 
    //!                 the inner corner and inner side)
    int cellType;

    //! bc index.
    //! The first is the index of face (from 0 to nFace-1 of cell  *2).
    //! The seonde is the bc Type.
    map<int,int> bcIndex;

public:
    //! Constructor function, new set<int>
    IndexParticleOfCell();

    //! Copy Constructor, new set<int>
    //! and insert IndexParticleOfCell.
    //! Note here pass by value.
    IndexParticleOfCell(const IndexParticleOfCell &rhs);

    //! Destructor function , delete the set.
    ~IndexParticleOfCell();

public:
    //! Judge is the existence of particles.
    bool IsExistParticleInCell(int particleID);
    //! Judge is the set empty. 
    bool IsEmpty();

    //! Print particle ID in current set on windows.
    void PrintParticleIDInCell();

    //! Add particle ID in current set.
    void AddParticleID(int particleID);
    //! Remove particle ID in current set.
    void RemoveParticleID(int particleID);

    //! Return the size of set.
    int GetNumOfParticleOfCell();

    //! Return the particle ID array,
    //! note here pass by value.
    //! If the set is empty, return 0.
    int *GetParticleIDOfCell();

    //! Clear set.
    void ClearParticle();
    //! Delete set.
    void DeleteSet();

    //! Init map<int,int> bcIndex by nFace on 0.
    void InitBCIndex(int nFace);

    void SetBCType(int iFace,int iBC);

    void SetCellType(int cellType);

    int GetCellType();

    map<int, int> *GetBCType();

    int GetBCType(int iFace);
};

}