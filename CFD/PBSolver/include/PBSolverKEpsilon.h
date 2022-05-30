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

#include "PBSolverScalar.h"
#include "PBSolverTurbK.h"
#include "PBSolverTurbEpsilon.h"

using namespace std;

namespace PHSPACE
{
#define SOLVER_TURB_KEPSILON TURB_KEPSILON

	class PBSolverKEpsilon : public PBSolverScalar
	{
	public:
		PBSolverKEpsilon(const std::string& name);
		virtual ~PBSolverKEpsilon();

		virtual void InitMemory();
		virtual void ReinitField();
		virtual void UpdateUnsteadyVariable(Grid* grid);
		virtual void SolveKEquation(Region& region);
		virtual void SolveEpsilonEquation(Region& region);
		virtual void UpdateAfterIterloop();
		virtual void ComputePostVisualVariables(Post_Visual* postVisualization);
		virtual void InitIncompressibleResultFile(Grid* grid, string fileName);
		virtual void DumpIncompressibleResultFile(Grid* grid, string fileName);
		virtual void SolutionIsConverged(int iZone, bool* flag);
		virtual void SetGeometry(PHGeometry* geometry) { this->geometry = geometry; }

	protected:
		vector<PHSolver*> solverList;
		CFDSolver* solverKEqn;
		CFDSolver* solverEpsilonEqn;
		PHGeometry* geometry;
	};

	DECLARE_PBSOLVER_REGISTER_FUNCTION(PBSolver, PBSolverKEpsilon, SOLVER_TURB_KEPSILON);
}