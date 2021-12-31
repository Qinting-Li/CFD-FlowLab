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
//! @file      GridDeformationSPRING.h
//! @brief     Grid deformation by SPRING method.
//! @author    Bell, Baka.

//! Zhao Zhong, Duan Xupeng, Chang Xinghua, etc. Numerical simulation of unsteady 
//! flows on bird-like flapping wing[J]. Transactions of Nanjing University of 
//! Aeronautics & Astronautics, 2013, 30(sup):93-100.

#pragma once
#include "GridDeformation.h"

namespace PHSPACE
{

//! @brief GridDeformationSPRING class achieve unstructured grid deform by SPRING method.
class GridDeformationSPRING : public GridDeformation
{
public:
    //! @param[in] nZonesIn         Number of grid.
    //! @param[in] stationalGridIn  Original grid.
    GridDeformationSPRING(int nZonesIn, Grid **stationalGridIn);
    ~GridDeformationSPRING();

private:
    void Deforming();
    void SurfaceGridMove(int iStep);
    void TransformGrid();
    void MatchControlPoints();
    void PostTreat();
};

}
