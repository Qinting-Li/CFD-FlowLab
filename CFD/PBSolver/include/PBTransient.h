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
//! @file      PBTransient.h
//! @brief     Discretization of the Transient Term.
//! @author    Ming Pingjian, Sun Huawen.

#pragma once

#include "PBObject.h"

namespace PHSPACE
{
#define TRANSIENT_EMPTY EMPTY

	class UnstructGrid;

	class PBTransient : public PBObject
	{
	public:
		PBTransient(const string& name);
		virtual ~PBTransient();

	public:
		virtual void AllocateMemory(UnstructGrid* grid, vector<string>& phiName, vector<RDouble>& initPhi) {};

		virtual void Discrt(UnstructGrid* grid, string phiName, string phiOldName, string phiOldOldName, string rhoOldName = NULL, string rhoOldOldName = NULL) {};

		virtual void SaveOldTimeValue(UnstructGrid* grid, vector<string>& phiName) {};
	};

	DEFINE_PBOBJECT_FACTORY(PBTransient);

	DECLARE_PBSOLVER_REGISTER_FUNCTION(PBTransient, PBTransient, TRANSIENT_EMPTY);
}
