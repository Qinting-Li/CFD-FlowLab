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
//! @file      GridgenIn.h
//! @brief     Explain this file briefly.
//! @author    He Xin.

#pragma once
using namespace std;

namespace PHSPACE
{
class Grid;

//! To read the grid such as coordinate and boundary condition information from plot3d format file.
//! @param[in ]: structGrids denotes the current computational regions of the structure grids.
//! @param[in ]: numberOfBlocks denotes the number of total grid blocks.
//! @param[in ]: gridgenFile denotes the source file in which the grid stored.
//! @param[in ]: BCFile denotes the source file in which the boundary information stored.
void ReadGridgen(Grid **& structGrids, int &numberOfBlocks, const string &gridgenFile, const string &BCFile);

//! To read the grid coordinate from plot3d binary format file.
//! @param[in ]: structGrids denotes the current computational regions of the structure grids.
//! @param[in ]: numberOfBlocks denotes the number of total grid blocks.
//! @param[in ]: coordinateFile denotes the source file in which the grid coordinate stored.
void ReadGridgenCoor(const string &coordinateFile, int &numberOfBlocks, Grid **&structGrids);

//! To read the grid coordinate from plot3d ASCII format file.
//! @param[in ]: structGrids denotes the current computational regions of the structure grids.
//! @param[in ]: numberOfBlocks denotes the number of total grid blocks.
//! @param[in ]: coordinateFile denotes the source file in which the grid coordinate stored.
int  CheckGridgenCoorASCIIFileType(const string &coordinateFile);
void ReadGridgenCoorASCIIPlot3D(const string &coordinateFile, int &numberOfBlocks, Grid **&structGrids);
void ReadGridgenCoorASCIIGridgen(const string &coordinateFile, int &numberOfBlocks, Grid **&structGrids);

//! To read the grid boundary information from ".inp" format file.
//! @param[in ]: structGrids denotes the current computational regions of the structure grids.
//! @param[in ]: numberOfBlocks denotes the number of total grid blocks.
//! @param[in ]: BCFile denotes the source file in which the grid boundary condition information stored.
void ReadGridgenBC(Grid **& structGrids, int &numberOfBlocks, const string &BCFile);

int CheckGridgenFileIfASCII(const string &coordinateFile);
}
