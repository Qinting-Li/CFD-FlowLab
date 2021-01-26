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
//! @file      Geo_CellMetrics_Unstruct.h
//! @brief     It defines the cell metrics of the unstructured grid,
//!			   such as the cell center and cell volume.
//! @author    Bell, He Xin.
#pragma once
#include "LIB_Macro.h"
#include "Precision.h"
namespace PHSPACE
{
//! @brief Geo_CellMetrics_Unstruct class defines the cell metrics of unstructured grid.\n
//! 1: cell center.\n
//! 2: cell volume.
class Geo_CellMetrics_Unstruct
{
private:
	//! Cell center data, including ghost cells.
	RDouble * xcc, * ycc, * zcc;

	//! Cell volume data, excluding ghosts.
	RDouble * vol; 

	//! Used for DES simulation ONLY: largest local grid length, distance of neighbor cell centers.
	RDouble * largestLocalGridLength;

	//! Used for DES simulation ONLY: sub-grid length.
	RDouble * subgridLength;
public:
	LIB_EXPORT Geo_CellMetrics_Unstruct();
	LIB_EXPORT ~Geo_CellMetrics_Unstruct();

public:
	//! Assigns the given pointer to the X Cell center data.
	void SetCellCenterX(RDouble * xcc);

	//! Assigns the given pointer to the Y Cell center data.
	void SetCellCenterY(RDouble * ycc);

	//! Assigns the given pointer to the Z Cell center data.
	void SetCellCenterZ(RDouble * zcc);

	//! Assigns the given pointer to the Cell volume data.
	void SetCellVolume(RDouble * vol);

	//! Assigns the largest local grid spacing.
	void SetLargestLocalGridLength(RDouble * length);

	//! Assigns the sub-grid length, used for DES simulation ONLY.
	void SetSubgridLength(RDouble * length);

	//! Returns the X pointer of cell center.
	RDouble * GetCellCenterX() const;

	//! Returns the X pointer of cell center.
	RDouble * GetCellCenterY() const;

	//! Returns the X pointer of cell center.
	RDouble * GetCellCenterZ() const;

	//! Returns the cell volume pointer.
	RDouble * GetCellVolume()  const;

	//! Returns the largest local grid spacing.
	RDouble * GetLargestLocalGridLength();

	//! Returns the sub-grid length, used for DES simulation ONLY.
	RDouble * GetSubgridLength();
};

#include "Geo_CellMetrics_Unstruct.hxx"

}