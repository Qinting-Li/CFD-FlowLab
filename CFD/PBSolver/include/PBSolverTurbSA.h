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
//! @brief     TurbSA.
//! @author    Ming Pingjian, Sun Huawen, Song Yang, Huang Jin.

#pragma once

#include <vector>
#include "PBSolverScalar.h"

using namespace std;

namespace PHSPACE
{
#define SOLVER_TURB_SA TURB_SA

	class PBSolverTurbSA : public PBSolverScalar
	{
	public:
		PBSolverTurbSA(const string& name);
		virtual ~PBSolverTurbSA();

		virtual void InitMemory();
		virtual void AllocateMemory();
		virtual void ReinitField();
		virtual void Discrt();
		virtual void UpdateAfterIterloop();//! Updates the viscosity of the momentum equation
		virtual void ComputePostVisualVariables(Post_Visual* postVisualization);
		virtual void OutputResidual(ofstream& file, int flag);

	protected:
		virtual void CalcGenAndDesTerm(UnstructGrid* grid);
		
		virtual void CalcGraTerm(UnstructGrid* grid);



	protected:
		RDouble delatv;//! 2/3
		RDouble Cb2;//! 0.622
		RDouble Cv1;//! 7.1
		RDouble Cb1;//! 0.1355
		RDouble kappa;//! 0.4187
		RDouble Cw2;//! 0.3
		RDouble Cw3;//! 2
		RDouble urfMu;//! Viscosity relaxation factor
		RDouble* k;//! Kinematic viscosity
		RDouble* kb;//! Kinematic viscosity on the wall
		RDouble* visl;//! Laminar viscosity
		RDouble* vist;//! Turbulent viscosity
		RDouble* des;//! Reduced turbulence viscosity
		RDouble* gra;//! Gradient term
		RDouble* gen;//! Turbulent viscosity produces terms
		RDouble* fv1;//! Viscous damping coefficient
		RDouble* dkdx;//! The gradient of Kinematic viscosity
		RDouble* dkdy;
		RDouble* dkdz;
		RDouble* WD;//! Wall distance

	};

	DECLARE_PBSOLVER_REGISTER_FUNCTION(PBSolver, PBSolverTurbSA, SOLVER_TURB_SA);
}

