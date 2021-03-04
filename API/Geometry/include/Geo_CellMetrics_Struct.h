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
//! @file      Geo_CellMetrics_Struct.h
//! @brief     It defines the cell metrics of the structured grid,
//!            such as the cell center and cell volume.
//! @author    Bell, He Xin.
#pragma once
#include "TypeDefine.h"

namespace PHSPACE
{
//! @brief Geo_CellMetrics_Struct class defines the cell metrics of structured grid.\n
//! 1: cell center. \n
//! 2: cell volume. \n
class Geo_CellMetrics_Struct
{
private:
	//! Cell center
	RDouble3D * xcc, * ycc, * zcc;

	//! Cell volume
	RDouble3D * vol;

	//! Cell length
	RDouble3D * xlen, * ylen, * zlen;

	//! Add by Min Yaobing.
	RDouble3D * jacobian;
	//! Add by Min Yaobing.

    //! Used for DES simulation ONLY: largest local grid length, distance of neighbor cell centers.
	RDouble3D * largestLocalGridLength;

	//! Used for DES simulation ONLY: sub-grid length.
	RDouble3D * subgridLength;

public:
	Geo_CellMetrics_Struct();
	~Geo_CellMetrics_Struct();

public:
	//! Returns the X component of cell center (xcc).
	RDouble3D * GetCellCenterX() const;

	//! Returns the Y component of cell center (ycc).
	RDouble3D * GetCellCenterY() const;

	//! Returns the Z component of cell center (zcc).
	RDouble3D * GetCellCenterZ() const;

	//! Returns the cell volume (vol).
	RDouble3D * GetCellVolume() const;

	//! Return the X direction cell length.
	RDouble3D * GetCellLengthX() const;

	//! Return the Y direction cell length.
	RDouble3D * GetCellLengthY() const;

	//! Return the Z direction cell length.
	RDouble3D * GetCellLengthZ() const;

	//! Add by Min Yaobing.
	RDouble3D * GetCellJacobian() const;
	//! Add by Min Yaobing.

    //! Returns the largest local grid spacing.
    RDouble3D * GetLargestLocalGridLength();

	//! Returns the sub-grid length, used for DES simulation ONLY.
	RDouble3D * GetSubgridLength();

    //! Assigns the largest local grid spacing.
    void SetLargestLocalGridLength(RDouble3D * length);

	//! Assigns the sub-grid length, used for DES simulation ONLY.
	void SetSubgridLength(RDouble3D * length);

	//! Init all the metrics pointers ( allocate memory ). 
	//! @param[in] I
	//! @param[in] J
	//! @param[in] K
	//! @note construct the three dimensional array with Range I, J, K.
	void AllocateMetrics(Range &I, Range &J, Range &K);

private:
	//! Deallocate memory of all the metrics array pointers.
	void DeallocateAll();
};
#include "Geo_CellMetrics_Struct.hxx"
}