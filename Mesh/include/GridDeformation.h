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
//! @file      GridDeformation.h
//! @brief     Grid deformation.
//! @author    Bell, Baka.

//! Zhao Zhong, Ma Rong, He Lei, Chang Xinghua, Zhang Laiping. An efficient large-scale 
//! mesh deformation method based on MPI/OpenMP hybrid parallel radial basis function 
//! interpolation[J]. Chinese Journal of Aeronautics, DOI: 10.1016/j.cja.2019.12.025.

#pragma once
#include "DynamicGrid.h"

namespace PHSPACE
{

//! @brief DeformParameter class storage the parameter of grid deform.
class DeformParameter
{
private:
    //! Number of deform step.
    int nDeformStep;

    //! The max flap angle
    RDouble flapAngle;

    //! Rotate postion.
    RDouble rotatePostionZ, rotatePostionY;

    //! If dump slice grid.
    int gridSlice;

    //! Grid slice axis
    int sliceAxis;

    //! Grid slice position.
    RDouble slicePosition;

public:
    void SetNDeformStep(int nDeformStepIn);
    void SetFlapAngle(RDouble flapAngleIn);
    void SetRotatePostionZ(RDouble rotatePostionZIn);
    void SetRotatePostionY(RDouble rotatePostionYIn);
    void SetGridSlice(int gridSliceIn);
    void SetSliceAxis(int sliceAxisIn);
    void SetSlicePosition(RDouble slicePositionIn);

    int GetNDeformStep();
    RDouble GetFlapAngle();
    RDouble GetRotatePostionZ();
    RDouble GetRotatePostionY();
    int GetGridSlice();
    int GetSliceAxis();
    RDouble GetSlicePosition();
};

//! @brief GridDeformation achieve unstructured grid deform.
class GridDeformation
{
public:
    //! @param[in] nZonesIn         Number of grid.
    //! @param[in] stationalGridIn  Original grid.
    GridDeformation(int nZonesIn, Grid **stationalGridIn);
    ~GridDeformation();

public:
    void Run();

    //! @param[in] deformParameterIn   Grid deform parameter.
    void SetDeformParameter(DeformParameter *deformParameterIn);

    virtual void SetNumberOfReferenceCP(int numberOfReferenceCPIn) {};
    virtual void SetInfluencePara(int influenceParaIn) {};

private:
    //! Transform stationalGrid to deformedGrid.
    virtual void TransformGrid() = 0;

    //! Grid deformed.
    virtual void Deforming() = 0;

    //! Set vector of surface node.
    virtual void SurfaceGridMove(int iStep) = 0;

    virtual void MatchControlPoints() = 0;
    virtual void PostTreat() = 0;

protected:
    void SetNodeAttribute(DynamicGrid *grid);

    //! Set z coordinate of node(in symmetry face) to zero.
    void SetSymmetryToZero();

    //! Dump surface and slice grid into tecpolt file.
    void PostVisual(string visualFileName, int iStep);

protected:
    //! Number of grid.
    int nZones;

    //! Original grid.
    Grid **stationalGrid;

    //! Deformed grid.
    DynamicGrid **deformedGrid;

    //! Grid deform parameter.
    DeformParameter *deformParameter;
};

void DeformGrid();

const int SPRING = 1;
const int RBF    = 2;
}
