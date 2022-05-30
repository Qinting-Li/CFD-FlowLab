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
//! @file      PBBoundaryTurbEpsilonPressureOutlet.h
//! @brief     Boundary Conditions.
//! @author    Ming Pingjian, Sun Huawen, Song Yang.

#pragma once

#include "Geo_UnstructGrid.h"
#include "Data_Param.h"
#include "PBBoundaryTurbDefault.h"

namespace PHSPACE
{
#define TURB_EPSILON_PRESSUREOUTLET TURB_EPSILON_PRESSUREOUTLET 

	class PBBoundaryTurbEpsilonPressureOutlet : public PBBoundaryTurbDefault
	{
	public:
		PBBoundaryTurbEpsilonPressureOutlet(const string& name);
		virtual ~PBBoundaryTurbEpsilonPressureOutlet();

		virtual void Update(UnstructGrid* grid, string phiName, vector<int>* faceIndex, Data_Param* bcData);

		virtual void Discrt(UnstructGrid* grid, string phiName, vector<int>* faceIndex, Data_Param* bcData, RDouble* diffCoe);

	protected:
		RDouble cmu;
	};


	DECLARE_PBSOLVER_REGISTER_FUNCTION(PBBoundary, PBBoundaryTurbEpsilonPressureOutlet, TURB_EPSILON_PRESSUREOUTLET);
} // PB_SPACE
