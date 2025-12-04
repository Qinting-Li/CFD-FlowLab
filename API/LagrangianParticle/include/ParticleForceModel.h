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
//! @file      FoceModel.h
//! @brief     It is the class for only one particle force model.
//! @author    Lei Sanhui, Lei Yinghaonan (Lanzhou University).

#include "OnePointVariable.h"
#include "ParticleAnalysis.h"
#include "SimplePointer.h"  

# pragma warning (disable:4819)

namespace PHSPACE
{

namespace PARTICLE_FORCETYPE
{
//! Do not calculate this force
const int NO_FORCE = 0;

//! =================================================================
//! ===                  Gravity force model                      ===
//! =================================================================
const int GRAVITY = 1;

//! =================================================================
//! ===                    Drag force model                       ===
//! =================================================================
//! A Drag Coefficient from Stokes 1850.
const int DRAG_STOKES = 1;
//! A Drag Coefficient from SCHILLER,1933.
const int DRAG_SN_1933 = 2;
//! A Drag Coefficient from XingLi 2015.
const int DRAG_Xing_2015 = 3;

//! =================================================================
//! ===                   Lift Grad force model                   ===
//! =================================================================
//! Saffman force from XingLi 2015.
const int SAFFMAX_XING_2015 = 1;

//! =================================================================
//! ===                  Lift Rotate force model                  ===
//! =================================================================
//! Magnus force from RuyangLI 2020.
const int MAGNUS_RUYANLI_2020 = 1;

//! =================================================================
//! ===                    Addmass  force model                   ===
//! =================================================================
//! Addmass force from Reza Barati 2018.
const int ADDMASS_REZA_2018 = 1;

//! =================================================================
//! ===              Fluid acceleration force model              ===
//! =================================================================
//! Fluid Acceleration force from Reza Barati 2018.
const int FlLUID_ACCELERATION_REZA_2018 = 1;

//! =================================================================
//! ===                   Brownian force model                    ===
//! =================================================================
//! Brownian fore from XingLi 2015.
const int BROWMIAN_XING_2015 = 1;

//! =================================================================
//! ===                 Thermophoretic force model                ===
//! =================================================================
//!  Thermophoretic force from XingLi 2015.
const int THERMOPHORETIC_XING_2015 = 1;

//! Besset force

//! =================================================================
//! ===           particle TemperatureModel  model                ===
//! =================================================================
const int TEMPERATURE_CHING_2020 = 1;

}

class ParticleForceModel
{
private:

public:
    ParticleForceModel();
    ~ParticleForceModel();
public:
    static void CalcParticleForce(OnePointVariable *onePointVariable, Data_Param *parParam);

private:

    //! =================================================================
    //! ===                  Gravity force model                      ===
    //! =================================================================
    static void CalcGravityForce(OnePointVariable *onePointVariable, Data_Param *parParam);
    static void CalcGravityForceOnDirection(OnePointVariable *onePointVariable, Data_Param *parParam);

    //! =================================================================
    //! ===                    Drag force model                       ===
    //! =================================================================
    static void CalcDragForce(OnePointVariable *onePointVariable, Data_Param *parParam);

    //! A Drag Coefficient from Stokes 1850.
    static void CalcDragForceStokes(OnePointVariable *onePointVariable, Data_Param *parParam);

    //! A Drag Coefficient from SCHILLER,1933.
    static void CalcDragForceSN(OnePointVariable *onePointVariable, Data_Param *parParam);

    //! A Drag Coefficient from XingLi 2015.
    static void CalcDragForceXingLi(OnePointVariable *onePointVariable, Data_Param *parParam);
    //! A Drag Coefficient from XingLi 2015.
    //! Solve nonlinear equations about slip coefficients on XingLi 2015.
    //! slipCoefficentK is the variable which to calculate.
    //! This return the result of function on slipCoefficentK.
    static RDouble atknf(RDouble slipCoefficentK, RDouble particleRe, RDouble knNumber, RDouble maNumber);

    //! =================================================================
    //! ===                   Lift Grad force model                   ===
    //! =================================================================
    static void CalcSaffmanForce(OnePointVariable *onePointVariable, Data_Param *parParam);

    //! Saffman force from XingLi 2015.
    static void CalcSaffmanForceXingLi(OnePointVariable *onePointVariable, Data_Param *parParam);

    //! =================================================================
    //! ===                  Lift Rotate force model                  ===
    //! =================================================================
    static void CalcMagnusForce(OnePointVariable *onePointVariable, Data_Param *parParam);

    //! Magnus force from RuyangLI 2020.
    static void CalcMagnusForceRuyangLI(OnePointVariable *onePointVariable, Data_Param *parParam);

    //! =================================================================
    //! ===                    Addmass  force model                   ===
    //! =================================================================
    static void CalcAddmassForce(OnePointVariable *onePointVariable, Data_Param *parParam);

    //! Addmass force from Reza Barati 2018.
    static void CalcAddmassForceRezaBarati(OnePointVariable *onePointVariable, Data_Param *parParam);

    //! =================================================================
    //! ===              Fluid acceleration  force model              ===
    //! =================================================================
    //! fluid acceleration foece model;
    static void CalcFluidAccelerationForce(OnePointVariable *onePointVariable, Data_Param *parParam);

    //! Fluid Acceleration force from Reza Barati 2018.
    static void CalcFluidAccelerationForceRezaBarati(OnePointVariable *onePointVariable, Data_Param *parParam);

    //! =================================================================
    //! ===                   Brownian force model                    ===
    //! =================================================================
    static void CalcBrownianForce(OnePointVariable *onePointVariable, Data_Param *parParam);

    //! Brownian fore from XingLi 2015.
    static void CalcBrownianForceXingLi(OnePointVariable *onePointVariable, Data_Param *parParam);
    //! Brownian fore from XingLi 2015.
    //! Calculate the spectral intensity.
    static void CalculatespectralIntensity(RDouble *sI, SPDouble &flowOnParticle, RDouble diameter, RDouble density);
    //! Generates a Gaussian random number
    static RDouble CalculateKeci();
    //! Calculate the uniform. Used by CalculateKeci().
    static void Uniform(double *p);

    //! =================================================================
    //! ===                 Thermophoretic force model                ===
    //! =================================================================
    //! Thermophoretic force model.
    static void CalcThermophoreticForce(OnePointVariable *onePointVariable, Data_Param *parParam);

    //! Thermophoretic force from XingLi 2015.
    static void CalcThermophoreticForceXingLi(OnePointVariable *onePointVariable, Data_Param *parParam);
    //! Thermophoretic force from XingLi 2015.
    //！Calculate the coefficient by Thermophoretic.
    static void CalcThermodynamicFactor(SPDouble &factorAH, RDouble knNumber);
    //！ Using the five Lagrange extrapolation.Used by CalculatethermodynamicFactor.
    static RDouble CalcLagrangeInterpolation(RDouble x[], RDouble y[], RDouble xx);
    //! Using cubic spline interpolation which used by CalculatethermodynamicFactor.
    static RDouble CalcSpline(RDouble x[], RDouble y[], int s, RDouble x1, RDouble xn, RDouble t);
    //! Using shooting method to solve linear equations which used by CalculateSpline.
    static RDouble capway(RDouble nu[], RDouble d[], RDouble lan[], RDouble h[], RDouble x[], RDouble y[], RDouble t);

    //! =================================================================
    //! ===               Particle temperature model                  ===
    //! =================================================================
    static void CalcParticleTemperature(OnePointVariable *onePointVariable, Data_Param *parParam);
    
    static void CalcParticleTemperatureChing2020(OnePointVariable *onePointVariable, Data_Param *parParam);

};

}