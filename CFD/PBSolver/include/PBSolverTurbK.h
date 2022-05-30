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
//! @file      PBSolverTurbK.h
//! @brief     TurbK.
//! @author    Ming Pingjian, Sun Huawen, Song Yang, Huang Jin.

#pragma once

#include <vector>
#include "PBSolverScalar.h"

using namespace std;

namespace PHSPACE
{
#define SOLVER_TURB_K TURB_K

	class PBSolverTurbK : public PBSolverScalar
	{
	public:
		PBSolverTurbK(const string& name);
		virtual ~PBSolverTurbK();

		virtual void InitMemory();
		virtual void AllocateMemory();
		virtual void ReinitField();
		virtual void Discrt();
		virtual void UpdateAfterIterloop();
		virtual void ComputePostVisualVariables(Post_Visual* postVisualization);
		virtual void OutputResidual(ofstream& file, int flag);

	protected:
		virtual void CalcGenTerm();

	protected:
		RDouble cmu; 
		RDouble prandtl;//! the Prandtl number of Turbulent kinetic energy
		RDouble urfMu; //! Viscosity relaxation factor
		RDouble* k;//! Turbulent kinetic energy
		RDouble* kb;//! Turbulent kinetic energy
		RDouble* epsilon;//! Turbulent dissipation rate
		RDouble* epsilonb;//! Turbulent dissipation rate on wall
		RDouble* visl;//! Laminar viscosity
		RDouble* vist;//! Turbulent viscosity
		RDouble* gen;//! Turbulent kinetic energy generation term
	};

	DECLARE_PBSOLVER_REGISTER_FUNCTION(PBSolver, PBSolverTurbK, SOLVER_TURB_K);
}
