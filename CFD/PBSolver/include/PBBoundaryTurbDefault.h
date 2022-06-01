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
//! @file      PBBoundaryTurbDefault.h
//! @brief     Boundary Conditions.
//! @author    Ming Pingjian, Sun Huawen, Song Yang.

#pragma once

#include "Geo_UnstructGrid.h"
#include "Data_Param.h"
#include "PBBoundaryScalarDefault.h"

namespace PHSPACE
{
#define TURB_DEFAULT TURB_DEFAULT

	class PBBoundaryTurbDefault : public PBBoundaryScalarDefault
	{
	public:
		PBBoundaryTurbDefault(const string& name);
		virtual ~PBBoundaryTurbDefault();
	};


	DECLARE_PBSOLVER_REGISTER_FUNCTION(PBBoundary, PBBoundaryTurbDefault, TURB_DEFAULT);
} // PB_SPACE
