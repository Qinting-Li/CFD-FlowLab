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
//! @file      Geo_DynamicGridMetrics_Unstruct.h
//! @brief     It defines the dynamic grid metrics of the unstructured grid,
//!            such as the face normal velocity.
//! @author    Bell (Modified).

#pragma once
#include "LIB_Macro.h"
#include "Precision.h"
namespace PHSPACE
{
//! @brief Geo_DynamicGridMetrics_Unstruct class defines dynamic mesh metrics of unstructured grid.\n
//! 1: face velocity.\n
//! 2: node velocity.
class Geo_DynamicGridMetrics_Unstruct
{
private:
	//! Cell volume of the last time step.
	RDouble * voln;

	//! Face normal velocity.
	RDouble * vgn;

	//! Node velocity.
	RDouble * xtn, * ytn, * ztn;

public:
	LIB_EXPORT Geo_DynamicGridMetrics_Unstruct();
	LIB_EXPORT ~Geo_DynamicGridMetrics_Unstruct();

public:
	//! Assign the given pointer to the X node velocity.
	void SetNodeVelocityX(RDouble * xtn);

	//! Assign the given pointer to the Y node velocity.
	void SetNodeVelocityY(RDouble * ytn);

	//! Assign the given pointer to the Z node velocity.	
	void SetNodeVelocityZ(RDouble * ztn);
	
	//! Assign the given pointer to the face normal velocity.
	void SetFaceNormalVelocity(RDouble * vgn);

	//! Assign the given pointer to the Cell volume of the last time step.
	void SetCellVolumeOld(RDouble * voln);

	//! Return the X node velocity pointer.
	RDouble * GetNodeVelocityX() const;

	//! Return the Y node velocity pointer.
	RDouble * GetNodeVelocityY() const;

	//! Return the Z node velocity pointer.
	RDouble * GetNodeVelocityZ() const;

	//! Return the face normal velocity pointer.
	RDouble * GetFaceNormalVelocity() const;

	//! Return the cell volume of the last time step.
	RDouble * GetCellVolumeOld() const;
};

#include "Geo_DynamicGridMetrics_Unstruct.hxx"





}