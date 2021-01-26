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
//! @file      Boundary.h
//! @brief     boundary condition.
//! @author    Bell, He Xianyao, Wan Yunbo, He Xin.

#pragma once
#include "Precision.h"

namespace PHSPACE
{
void Pressure_InletBC (RDouble *prims, RDouble *primt, RDouble gama, RDouble mach_inf, RDouble totalP, RDouble totalT, RDouble* angle, int nl, int nchem);
void PressureInletBCRiemann (RDouble *prims, RDouble *primt, RDouble gama, RDouble mach_inf, 
							 RDouble nx, RDouble ny, RDouble nz, RDouble totalP, RDouble totalT, RDouble* angle, int nl, int nchem);
void PressureInletBCPlate (RDouble *prims, RDouble *primt, RDouble gama, RDouble mach_inf, RDouble totalP, RDouble totalT, RDouble* angle, int nl, int nchem);
void PressureOutletBC(RDouble *prims, RDouble *primt, RDouble pres, int nl, int nchem);

void InflowBC(RDouble *prims, RDouble *primt, int nl, int nchem);
void Isentropic_Inflow_BC(RDouble *prims, RDouble *primt, RDouble gama, int nl, int nchem);
 
void OutflowBC(RDouble *prims, RDouble *primt, int nl, int nchem);
void MassOutflow_BC(RDouble *prims, RDouble *primt, RDouble mass_alpha, int nl, int nchem);

void SymmetryBC(RDouble *prims, RDouble *primt, RDouble nxs, RDouble nys, RDouble nzs, RDouble vgn, int nl, int nchem);
void ViscousAdiabaticWallBC(RDouble *prims, RDouble *primt, RDouble nxs, RDouble nys, RDouble nzs, RDouble xtn, RDouble ytn, RDouble ztn, RDouble refMachNumber, int nl, int nchem);
void ViscousIsotropicWallBC(RDouble *prims, RDouble *primt, RDouble nxs, RDouble nys, RDouble nzs, RDouble xtn, RDouble ytn, RDouble ztn, RDouble refMachNumber, RDouble tw, int nl, int nchem);
void FarfieldBC(RDouble *prims, RDouble *prim_inf, RDouble *primt, RDouble gama0, RDouble gama, RDouble nxs, RDouble nys, RDouble nzs, RDouble vgn, int nl, int nchem);
void PreconFarfieldBC(RDouble *prims, RDouble *prim_inf, RDouble *primt, RDouble gama0, RDouble gama, RDouble nxs, RDouble nys, RDouble nzs, RDouble vgn, RDouble refMachNumber, int nl, int nchem);
void FarFieldRiemannInvariants( RDouble *prims, RDouble *prim_inf, RDouble *primt, RDouble gama0, RDouble gama, RDouble nxs, RDouble nys, RDouble nzs, RDouble vgn, int nl, int nchem );
}
