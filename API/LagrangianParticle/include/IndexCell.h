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
//! @file      IndexCell.h
//! @brief     It is the class for each cell id both on struct grid 
//!            and unstruct grid for each zone.
//! @author    Lei Yinghaonan (Lanzhou University).

#pragma once
#include "IndexParticleOfCell.h"
#include "OnePointVariable.h"

namespace PHSPACE
{
class Grid;
class StructGrid;
class UnstructGrid;

class IndexCell
{
private:
    //! Typedef for map iterator.
    typedef map<int, IndexParticleOfCell*>::iterator iterCell;

    //! The mapping relationship of cell local id and particle global id.
    //! The first varibal is th cell id on this zone.
    //! The second varibal is the global particle id on each cell.
    //! Initially, we kept the cell information for each grid (including Ghost). 
    //! However, this is not necessary. 
    //! Here, we only need to keep Ghost, Corner, 
    //! and the grid cells close to avoid.
    map<int, IndexParticleOfCell*> *cellParticleID;

    //! The number of cell (we need may not all) in current zone.
    int nCell;

    int nLayersBC;

public:
    //! Constructor function, new map
    //! The nCell will be init to -1.
    IndexCell();

    //! Destructor function , delete the map.
    ~IndexCell();

public:
    //! Init the nCell
    void InitNumOfCell(Grid *gridIn);

    //! Add empty set in each cell.
    void InitCellID(Grid *grid);
    void InitCellID(StructGrid *grid);
    void InitCellID(UnstructGrid *grid);

    //! Add particle ID whih search particle.
    //! if return -1, the particle is not in current zone.
    int InitIndexOfParticleInCell(Grid *grid, int idParticle,OnePointVariable *onePointVariable, bool &inZone);
    //! Add particle ID with search particle.
    void InitIndexOfParticleInCell(Grid *grid, map<int, OnePointVariable*> *particlePointVariable);

    //! Search Cell ID of one particle.
    int SearchParticleCell(Grid *grid, SimplePointer<RDouble> *onePointCoordinate,bool &inZone);
    int SearchParticleCell(StructGrid *grid, SimplePointer<RDouble> *onePointCoordinate, bool &inZone);
    int SearchParticleCell(UnstructGrid *grid, SimplePointer<RDouble> *onePointCoordinate,bool &inZone);

    //! Set Boundary index.
    void SetBoundaryCellIndex(Grid *grid);
    void SetBoundaryCellIndex(StructGrid *grid);
    void SetBoundaryCellIndex(UnstructGrid *grid);

    //! Add and remove particleID.
    void AddParticleID(int idCell,int idParticle);
    void RemoveParticleID(int idCell, int idParticle);

    //! Judge whether the cell's set is empty.
    bool IsCellEmpty(int idCell);
    //! Judge whether the particle is in idCell.
    bool IsParticleOnCell(int idCell,int idParticle);
    //! Judge whether the idCell is in current map.
    bool IsCellInMap(int idCell);

    map<int, IndexParticleOfCell*> *GetCellParticleID();

    //! Get the number of cell in current zone.
    int GetNumOfCell();
    //! Get the number of particle in idCell.
    int GetNumOfParticleOnCell(int idCell);
    //! Get the particle ID in idCell.
    int *GetIDOfParticleForCell(int idCell);

    //! Delete Map and 
    void DeleteMap();

    //! If particle in IndexCell.
    bool SetParticleInCellInfo(Grid *grid,OnePointVariable *onePointVariable);
    bool SetParticleInCellInfo(StructGrid *grid, OnePointVariable *onePointVariable);
    bool SetParticleInCellInfo(UnstructGrid *grid, OnePointVariable *onePointVariable);

private:

};

void GetIndexCellIJKStruct(int &nI, int &nJ, int &nK, int &iP, int &jP, int &kP, int nGhost, int &cellID);

void GetIDCellOfParticleByIJKSturct(int &nI, int &nJ, int &nK, int &iP, int &jP, int &kP, int nGhost, int &cellID);

}