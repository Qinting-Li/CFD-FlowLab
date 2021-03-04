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
//! @file      Geo_DynamicGridMetrics_Struct.h
//! @brief     It defines the dynamic grid metrics of the structured grid,
//!            such as the face normal velocity.
//! @author    Bell, Chang Xinghua, Zhang Jian, He Xin.

#pragma once
#include "TypeDefine.h"
namespace PHSPACE
{
//! @brief Geo_DynamicGridMetrics_Struct class defines the dynamic grid metrics of the structured grid.\n
//! 1: Cell volume of the last time step. \n
//! 2: Face normal velocity. \n
//! 3: Node velocity. \n
class Geo_DynamicGridMetrics_Struct
{
private:
	//! Cell volume of the last time step.
	RDouble3D * voln;

	//! Face normal velocity.
	RDouble4D * vgn;

	//! Node velocity.
	RDouble4D * xtn, * ytn, * ztn;
	
public:
	Geo_DynamicGridMetrics_Struct();
	~Geo_DynamicGridMetrics_Struct();

public:
	//! Returns the cell volume of the last time step (voln).
	RDouble3D * GetCellVolumeOld() const;

	//! Returns the Face normal velocity (vgn).
	RDouble4D * GetFaceNormalVelocity() const;

	//! Returns the x component of face velocity (xtn).
	RDouble4D * GetFaceVelocityX() const;

	//! Returns the y component of face velocity (ytn).
	RDouble4D * GetFaceVelocityY() const;

	//! Returns the z component of face velocity (ztn).
	RDouble4D * GetFaceVelocityZ() const;

	//! Assigns the cell volume of the last time step (voln).
	void SetCellVolumeOld(RDouble3D * voln);

	//! Assigns the Face normal velocity (vgn).
	void SetFaceNormalVelocity(RDouble4D * vgn);

	//! Assigns the x component of face velocity (xtn).
	void SetFaceVelocityX(RDouble4D * xtn);

	//! Assigns the y component of face velocity (ytn).
	void SetFaceVelocityY(RDouble4D * ytn);

	//! Assigns the z component of face velocity (ztn).
	void SetFaceVelocityZ(RDouble4D * ztn);
private:
	void DeallocateAll();

};
#include "Geo_DynamicGridMetrics_Struct.hxx"
}