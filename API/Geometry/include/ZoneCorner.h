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
//! @file      ZoneCorner.h
//! @brief     This module is used to store information about corner zones 
//!            of zone in Region. This is because particles may run to 
//!            ghost zone and corner zone when calculating in multi-block 
//!            grids, but this is not needed in traditional Euler continuum.
//! @author    Lei Yinghaonan (Lanzhou University).

#pragma once
#include "Precision.h"
#include "TypeDefine.h"
#include "TK_Exit.h"
#include "DataContainer.h"

namespace PHSPACE
{
class Grid;
class StructGrid;
class UnstructGrid;
class ZoneNeighbor;

class CornerGridIndex
{
private:

    //! The Corner point index for grid. 
    //! For struct grid, corner point is whitin corner zone
    int indexStr[3];

    //! For unstruct grid, corner point is every point.
    int indexUnstr[3];

public:
    CornerGridIndex();
    ~CornerGridIndex();
public:
    void SetCornerGirdIndex(Grid *grid, int &zoneCornerIndex, set<int> &zoneCornerNeighbor);
    void SetCornerGirdIndex(StructGrid *structGrid, int &zoneCornerIndex, set<int> &zoneCornerNeighbor);
    void SetCornerGirdIndex(UnstructGrid *unstructGrid, int &zoneCornerIndex, set<int> &zoneCornerNeighbor);
    
    int *GetIndex(Grid *grid);
    int *GetIndex(int gridType);

private:

};

//! The index for zoneCorner
//! This class include the information for only one corner.
class ZoneCornerIndex
{
private:

    //! The Global zone id of corner zone for current zone.
    int zoneCornerIndex;

    //! The corner index.
    CornerGridIndex cornerGridIndex;

    //! The neighbor of corner around center zone.    
    //! In 3d case, these may be ghost zone, so we need this
    //! variable to check.
    //! The first variable is global zone index of corner.
    set<int> zoneCornerNeighbor;

public:
    ZoneCornerIndex();
    ~ZoneCornerIndex();
public:

    //! Set the corner index.
    void SetZoneCornerIndex(int iCornerZoneIndex);

    void AddZoneCornerNeighborIndex(set<int> indexOfCornerNeighbor);

    void SetCornerGirdIndex(Grid *grid);

    set<int> &GetZoneCornerNeighbor();

    int *GetCornerGridIndex(int gridType);

    int *GetCornerGridIndex(Grid *grid);

    int GetCornerZoneIndex();

    //! Init the index of corner zone.
    //void InitZoneCornerIndex(int gridType, InterfaceInfo *iinfo, InterpointInformation *ipnfo);
private:

};

//! This class on ZoneConective.
//! Note that ZoneCorner is on only one zone.
class ZoneCorner
{
private:
    typedef set<int>::iterator zoneIndexIter;
    typedef map<int, set<int> >::iterator NeighborIndexIter;
    typedef map<int, ZoneCornerIndex >::iterator CornerIndexIter;
    //! The neighbor index of neighbor zone index on global.
    //! The first variable represents the global 
    //! neighbor zone index of the current zone.
    //! The second variable set<int> represents 
    //! the global number of the zone around each 
    //! adjacent zone of each neighbor zone.
    map<int, set<int> > neighborOfNeighbor;

    //! The corner information;
    //! The first variable represents the Corner Point zone index.
    //! The second variable represents the information of each corner.
    map<int, ZoneCornerIndex> cornerInfo;

    //! Global zone index.
    int globalZoneIndex;

    //! grid type for current zone.
    int gridType;

public:
    ZoneCorner();
    ~ZoneCorner();
public:
    //! Add the global zone index of neighbor's neighbor.
    void AddNeighborOfNeighbor(int iNeighbor, set<int> indexOfNeighbor);

    //! Add for unstruct grid.
    void AddCornerIndexForPoint(ZoneNeighbor *zoneNeighbor);

    //! Find the corner index by neighborOfNeighbor.
    void CalcCornerIndex();

    bool isExistCornerIndex(int iCornerIndex);

    void SetZoneCornerIndexGridType(int type);

    void SetZoneCornerIndexByGrid(Grid *grid);

    //! Init the index of corner zone.
    // void InitZoneCornerIndex(int gridType, InterfaceInfo *iinfo, InterpointInformation *ipnfo);

    void SetGlobalZoneIndex(int iZone);

    map<int, set<int> > &GetSetIndexNeighborOfNeighbor();

    map<int, ZoneCornerIndex> &GetCornerInfo();

    int GetNumOfCorner();

    //! Note that ,the interface variable iCorner 
    //! is a index from 0 to nCorner - 1, 
    //! which is not the corner of zone index.
    int GetZoneIndexOfCorner(int iCorner);

    map<int, ZoneCornerIndex >::iterator GetCornerIterBegin();
    map<int, ZoneCornerIndex >::iterator GetCornerIterEnd();

private:

};

void MainTaskBlockingCorner(ActionKey *actkey);

namespace FILLCORNER
{
void FillGridCorner(Grid *grid);

void FillGridCorner(StructGrid *grid);

void FillGridCorner(UnstructGrid *grid);

template < typename T >
void GhostCell3D(PHArray<T, 4> &w, int ni, int nj, int nk, int mst, int med)
{
    //! w(-1:ni+1, -1:nj+1, -1:nk+1, mst:med).
    if (nk != 1)
    {
        for (int m = mst; m <= med; ++m)
        {
            for (int k = 1; k <= nk - 1; ++k)
            {
                for (int j = 1; j <= nj - 1; ++j)
                {
                    w(0, j, k, m) = w(1, j, k, m);
                    w(ni, j, k, m) = w(ni - 1, j, k, m);
                }
                for (int i = 0; i <= ni; ++i)
                {
                    w(i, 0, k, m) = w(i, 1, k, m);
                    w(i, nj, k, m) = w(i, nj - 1, k, m);
                }
                //! Initialize volume in second row of ghost cells.
                for (int j = 0; j <= nj; ++j)
                {
                    w(-1, j, k, m) = w(0, j, k, m);
                    w(ni + 1, j, k, m) = w(ni, j, k, m);
                }

                for (int i = -1; i <= ni + 1; ++i)
                {
                    w(i, -1, k, m) = w(i, 0, k, m);
                    w(i, nj + 1, k, m) = w(i, nj, k, m);
                }
            }

            for (int j = -1; j <= nj + 1; ++j)
            {
                for (int i = -1; i <= ni + 1; ++i)
                {
                    w(i, j, 0, m) = w(i, j, 1, m);
                    w(i, j, -1, m) = w(i, j, 1, m);

                    w(i, j, nk, m) = w(i, j, nk - 1, m);
                    w(i, j, nk + 1, m) = w(i, j, nk - 1, m);
                }
            }
        }
    }
    else
    {
        int k = 1;
        for (int m = mst; m <= med; ++m)
        {
            for (int j = 1; j <= nj - 1; ++j)
            {
                w(0, j, k, m) = w(1, j, k, m);
                w(ni, j, k, m) = w(ni - 1, j, k, m);
            }

            for (int i = 0; i <= ni; ++i)
            {
                w(i, 0, k, m) = w(i, 1, k, m);
                w(i, nj, k, m) = w(i, nj - 1, k, m);
            }

            //! Initialize volume in second row of ghost cells.
            for (int j = 0; j <= nj; ++j)
            {
                w(-1, j, k, m) = w(0, j, k, m);
                w(ni + 1, j, k, m) = w(ni, j, k, m);
            }

            for (int i = -1; i <= ni + 1; ++i)
            {
                w(i, -1, k, m) = w(i, 0, k, m);
                w(i, nj + 1, k, m) = w(i, nj, k, m);
            }
        }
    }
}

}

}