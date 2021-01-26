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
//! @file      GridFactory.h
//! @brief     Some functions about grid.
//! @author    He Xin, Bell, Baka.

#pragma once
#include "MixGrid.h"

namespace PHSPACE
{
const int GRID_CREATION   = 0;
const int GRID_CONVERSION = 1;
const int GRID_ADAPTATION = 2;
const int GRID_MERGING    = 3;
const int GRID_BIRD       = 4;
const int GRID_REPAIR     = 5;
const int GRID_MIRROR     = 6;
const int GRID_FLAPPING   = 7;

void CreateGrid     ();
void ConvertGrid    ();
void GeneGridFromDataBase();
void ReconstructGrid();
void CombinationGrid();
void GeneMixGrid    ();
void OversetGridView();
void RepairGrid();
void MirrorGrid();
void CGNS2UnstructGrid();
void CGNS2StructGrid();
void Fantasy2Gridgen();
void Ustar2PHengLEI();
void PHengLEI2CGNS();

void ReadGMSH(CGNSFactory *&factory_cgns, int &nBlocks, const string &filename);
void GMSH2PHengLEI();
void GridgenToMIX();
void Gridgen2Fantasy();
void GridgenToSingleUns();
void GridgenToMBlockUns();
void GridgenToStr();

void GenDataBaseGrid(const string &gridfile, GENE_SIMPLE_GRID gene_simle_grid, MESH2D_GENE mesh2dgen, int gridtype);
}
