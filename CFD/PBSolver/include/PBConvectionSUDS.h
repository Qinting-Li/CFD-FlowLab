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
//! @file      PBConvectionSUDS.h
//! @brief     SUDS.
//! @author    Ming Pingjian, Sun Huawen, Song Yang.

#pragma once
#include "PBConvection.h"

namespace PHSPACE
{
#define CONVECTION_SUDS SUDS

	class PBConvectionSUDS : public PBConvection
	{
	public:
		PBConvectionSUDS(const string& name);
		virtual ~PBConvectionSUDS();

		//! @brief SOU scheme.
		//! @param[in]: gridIn denotes the current computational regions of the grids.
		//! @param[in]: phiName denotes the name of parameter. 
		//! @param[out]: apName denotes the leading diagonal.
		//! @param[out]: scName denotes the right-hand side. 
		virtual void Discrt(UnstructGrid* grid, string phiName, string apName, string scName, RDouble* v = nullptr,
			                RDouble* apV = nullptr, RDouble* scV = nullptr, RDouble* w = nullptr, RDouble* apW = nullptr, RDouble* scW = nullptr);
	};

	DECLARE_PBSOLVER_REGISTER_FUNCTION(PBConvection, PBConvectionSUDS, CONVECTION_SUDS);
};