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
//! @file      WallDist.h
//! @brief     Explain this file briefly.
//! @author    He Xin.

#pragma once
#include "ActionKey.h"

namespace PHSPACE
{
class Grid;

//! Process of generate wall distance.
void GenerateWallDist();

//! Compute wall distance.
void ComputeWallDist();

//!
void AllocateWalldist();

//! Write wall distance.
void WriteWallDist();

//! Dump out the wall distance to file.
void DumpWallDistToFile(ActionKey *actkey, int iZone);

#ifdef USE_TecplotLib
//! Visualization of wall distance.
void VisalWallDist();
#endif

//!
void WallDistStatics();

//! If read wall distance file or not.
bool IsReadWallDist();

}