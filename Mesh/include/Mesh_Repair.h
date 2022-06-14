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
//! @file      Mesh_Repair.h
//! @brief     Mesh repair for unstructured grid: for negative grid, repairing to valid quality.
//! @author    Bell.

#pragma once
#include "Geo_UnstructGrid.h"
#include "Zone.h"
using namespace std;

namespace PHSPACE
{

//! @brief Mesh repair for unstructured grid: for negative grid, repairing to valid quality.
class Mesh_Repair
{
public:
    //! @param[in] originalGridIn   Original grid, with negative volume cells, which need to be repaired..\n
    LIB_EXPORT Mesh_Repair(UnstructGrid *originalGridIn);

    LIB_EXPORT ~Mesh_Repair();

    //! Repair the original grid with negative to be valid grid.
    LIB_EXPORT void Repair();

    //! Return the agglomerated coarse grid.
    LIB_EXPORT UnstructGrid * GetRepairedGrid() { return this->repairedGrid; }

private:
    //! The original grid, with negative volume cells, which is need to be repaired.
    UnstructGrid *originalGrid;

    //! The repaired grid, which is valid for simulation.
    UnstructGrid *repairedGrid;
    
};

int PickOutBadSkewnessFaces(UnstructGrid *grid, bool *isBadFace, bool *isBadCell);

//! Repair grid to make sure sutiable mesh quality.
LIB_EXPORT Grid * RepairGrid(Grid *originalGridIn);

LIB_EXPORT void RepairGrid(Zone *zone, ActionKey *actkey);
}