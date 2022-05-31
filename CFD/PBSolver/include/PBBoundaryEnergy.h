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
//! @file      PBBoundaryEnergy.h
//! @brief     Boundary Conditions.
//! @author    Ming Pingjian, Sun Huawen, Huang Jin.

#pragma once

#include "PBBoundaryScalar.h"

namespace PHSPACE
{
#define ENERGY_EMPTY ENERGY_EMPTY 

	class PBBoundaryEnergy : public PBBoundaryScalar
	{
	public:
		PBBoundaryEnergy(const string& name);
		virtual ~PBBoundaryEnergy();

	public:
		virtual void Init(string phiName) {};

		virtual void Discrt(UnstructGrid* grid, string phiName, vector<int>* faceIndex, Data_Param* bcData, RDouble* diffCoe) {};

		virtual void Update(UnstructGrid* grid, string phiName, vector<int>* faceIndex, Data_Param* bcData) {};
	};

	DECLARE_PBSOLVER_REGISTER_FUNCTION(PBBoundary, PBBoundaryEnergy, ENERGY_EMPTY);
}
