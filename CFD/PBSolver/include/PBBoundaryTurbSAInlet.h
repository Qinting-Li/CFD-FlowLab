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
//! @file      PBBoundaryTurbSAInlet.h
//! @brief     Boundary Conditions.
//! @author    Ming Pingjian, Sun Huawen, Song Yang.

#pragma once

#include "Geo_UnstructGrid.h"
#include "Data_Param.h"
#include "PBBoundaryTurbDefault.h"

namespace PHSPACE
{
#define TURB_SA_INLET TURB_SA_INLET 

	class PBBoundaryTurbSAInlet : public PBBoundaryTurbDefault
	{
	public:
		PBBoundaryTurbSAInlet(const string& name);
		virtual ~PBBoundaryTurbSAInlet();

		virtual void Update(UnstructGrid* grid, string phiName, vector<int>* faceIndex, Data_Param* bcData);

		virtual void Discrt(UnstructGrid* grid, string phiName, vector<int>* faceIndex, Data_Param* bcData, RDouble* diffCoe);
	};

	DECLARE_PBSOLVER_REGISTER_FUNCTION(PBBoundary, PBBoundaryTurbSAInlet, TURB_SA_INLET);
} // PB_SPACE

