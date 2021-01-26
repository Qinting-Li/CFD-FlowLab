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
//! @file      Pre_GridRefine.h
//! @brief     Refining grid.
//! @author    He Xianyao.
#pragma once
#include "Region.h"
#include "Geometry.h"

using namespace std;

namespace PHSPACE
{
class Pre_GridRefine
{
public:
	string gridFileName;

public:
	LIB_EXPORT Pre_GridRefine(const string & from_gfile);
	LIB_EXPORT ~Pre_GridRefine();

public:
	LIB_EXPORT void Run();

private:
	void RefineGrid();
	void RefineGridStr(int izone, StructGrid * grid_in);
	void RefineBCRegion();

private:
	int numberOfBlocks;
	vector<Grid *> OrdinaryGrid;
	Grid ** NewGrids;
};

}