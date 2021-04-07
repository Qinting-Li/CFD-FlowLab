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
//! @file      LESSolverStruct.h
//! @brief     Large-eddy simulation solver for struct grid.
//! @author    Zhang Zipei.

#pragma once
#include "LESSolver.h"

namespace PHSPACE
{
class Grid;
class BCRegion;
class Param_LESSolverStruct;

class LESSolverStruct : public LESSolver
{
private:

public:
    LESSolverStruct();
    ~LESSolverStruct();

public:
    //! Create and initialize controlparameters.
    LIB_EXPORT void InitControlParameters();

    //! Get control paramters.
    LIB_EXPORT Param_LESSolverStruct * GetControlParameters();
private:
    void AllocateGlobalVar(Grid *gridIn);
    void DeAllocateGlobalVar(Grid *gridIn);

    void Init(Grid *grid_in);

    void ComputeCellLengthScale(Grid *gridIn);
    void ComputeCellLengthScaleforStructHighOrder(Grid *gridIn);
    void ComputeFaceLengthScale(Grid *gridIn);

    void ComputeWallFunction(Grid *gridIn);
    void WallFunctionofVanDriestOld(Grid *gridIn);
    void WallFunctionofVanDriest(Grid *gridIn);
    void WallFunctionofDXB(Grid *gridIn);
    void WallFunctionofPiomelli(Grid *gridIn);

    void ComputeViscousCoeff(Grid *gridIn);

    void Boundary(Grid *gridIn);

    //void ComputeSubgridScaleEnergy(Grid *gridIn);

    void SmagorinskyViscosity(Grid *gridIn);
    void Smagorinsky(Grid *gridIn);
    void SmagorinskyOld(Grid *gridIn);

    void DynamicSmagViscosity(Grid *gridIn);
    void DynamicSmagViscosityCompressible(Grid *gridIn);
    //! CHANG-YUE XU, LI-WEI CHEN AND XI-YUN LU, 2010, Large-eddy simulation of the compressible flow past a wavy cylinder.

    void DynamicSmagViscosityCompressibleOld(Grid *gridIn);
    void DynamicSmagViscosityIncompressible(Grid *gridIn);
    void DynamicSmagViscosityIncompressibleOld(Grid *gridIn);
    void ComputeTestFilteredQ(Grid *gridIn);
    void ComputeTestFilteredQfd(Grid *gridIn);
    //void ComputeTestFilteredDensity(Grid *gridIn);
    void ComputeTestFilteredTemperature(Grid *gridIn);
    void ComputeStrainRateTensor(Grid *gridIn);
    void ComputeStrainRateTensorforStructHighOrder(Grid *gridIn);
    void ComputeGradT(Grid *gridIn);
    void ComputeTestFilteredStrainRateTensor(Grid *gridIn);
    void ComputeTestFilteredRhoUiUj(Grid *gridIn);
    void ComputeLeonardStress(Grid *gridIn);
    void ComputeModeledStress(Grid *gridIn);
    void ComputeTestFilteredAlphaIJ(Grid *gridIn);
    void ComputeBetaIJ(Grid *gridIn);
    void ComputeAnisotropicConstant(Grid *gridIn);
    void ComputeIsotropicConstant(Grid *gridIn);
    void ComputePrandtlNumber(Grid *gridIn);
    void ComputeLeonardTemperatureTerm(Grid *gridIn);
    void ComputeModeledTemperatureTerm(Grid *gridIn);
    void ComputeTestFilteredRhoUT(Grid *gridIn);
    void ComputeTestFilteredGradT(Grid *gridIn);
    void ComputeTestFilteredRhoStrainRateMagnitudeGradT(Grid *gridIn);

    void WALE(Grid *grid);
    //! F. Nicoud and F. Ducros, 1999, Subgrid-Scale Stress Modelling Based on the Square of the Velocity Gradient Tensor.
    //! Christian Wollblad et.al., 2006, Large Eddy Simulation of Transonic Flow with Shock Wave/Turbulent Boundary Layer Interaction.

    void filter(const int & nn, RFloat1D * fu1, RFloat1D * fu2, RFloat1D * fu3
                ,RFloat1D * fuifuj11, RFloat1D * fuifuj22, RFloat1D * fuifuj33, RFloat1D * fuifuj12, RFloat1D * fuifuj13, RFloat1D * fuifuj23
                ,RFloat1D * fs11    , RFloat1D * fs22    , RFloat1D * fs33    , RFloat1D * fs12    , RFloat1D * fs13    , RFloat1D * fs23    
                ,RFloat1D * fsfs11  , RFloat1D * fsfs22  , RFloat1D * fsfs33  , RFloat1D * fsfs12  , RFloat1D * fsfs13  , RFloat1D * fsfs23);
    void TestFilter(const int &nn,RFloat1D *fu);
    void TestFilter(const int &n0, const int &nn, RDouble1D *&tf, RDouble1D *&f);
    void TestFilter(const int &n0, const int &nn, RDouble *tf, RDouble *f);

    //! To obtain the grid scale using in sub-grid model.
    void GetDeltaBar(const int iDeltaBar, const RDouble deltai, const RDouble deltaj, const RDouble deltak, RDouble &deltaBar);

    //! To calculate the non-isotriopic fix of filter width of Scotti&Meneveau (scotti1993, scotti1997, physics of fluid).
    void GetFilterWidthOfScotti(const RDouble deltai, const RDouble deltaj, const RDouble deltak, RDouble &deltaBar);
};

}