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
//! @file      Param_NSSolver.h
//! @brief     Record paramters of NS Solver.
//! @author    Bell, Zhang Jian, Wan Yunbo, Meng Liyuan.

#pragma once
#include "Param_CFDSolver.h"

namespace PHSPACE
{

class Param_NSSolver : public Param_CFDSolver
{
public:
	LIB_EXPORT Param_NSSolver();

	LIB_EXPORT ~Param_NSSolver();

public:
	//! Init all parameters.
	LIB_EXPORT void Init();

	//! Get interval steps of dumping force.
	int GetIntervalStepForce() const;

	//! Get wall temperature.
	RDouble GetWallTemperature() const; 

    //! Get the flag of the gas model.
    int GetChemicalFlag() const;

    //! Get the number of N-S equations without the species transport equations.
    int GetNSEquationNumber() const;

    //! Get the number of the equations including the species transport equations.
    int GetLaminarNumber() const;

    //! Get the flag of thermodynamic temperature model.
    int GetTemperatureModel() const;

    //! Get if carry out precondition or not.
	int GetIfLowSpeedPrecon() const;

	//! Get attacked of angle.
	RDouble GetAoA() const;

    //! Get the Entropy fix (correction) method.
    int GetRoeEntropyFixMethod() const;

    //! Get the Entropy coefficient scale, default is 1.0.
    RDouble GetRoeEntropyScale() const;

	//! Get the entropy fix coefficient 1 of Roe method.
	RDouble GetRoeEntropyFixCoef1() const; 
	
	//! Get the entropy fix coefficient 2 of Roe method.
	RDouble GetRoeEntropyFixCoef2() const; 

    //! Get angle of slide.
	RDouble GetAngleOfSlide () const;

    //! Set the entropy fix coefficients by the inflow mach number.
    void SetEntropyFixCoefficients();

    //! Get prandtlLaminar.
    RDouble GetPrandtlLaminar() const;

    //! Get prandtlTurbulence.
    RDouble GetPrandtlTurbulence() const;

    //! Get oPrandtlLaminar.
    RDouble GetoPrandtlLaminar() const;

    //! Get oPrandtlTurbulence.
    RDouble GetoPrandtlTurbulence() const;

    int GetNChemicalRadius() const;

    int GetNChemicalSource() const;

    //! Get timeIntegration.
    int GetTimeIntegration() const;

    //! Get nonDimensionalSutherlandTemperature.
    RDouble  GetNonDimensionalSutherlandTemperature() const;
    
    int GetMgProlongationType() const;

    RDouble GetMassFractionO2() const;

    RDouble GetMassFractionN2() const;

    RDouble * GetPrimitiveVarFarfield() const;

    RDouble GetSubsonicInletTotalPressure() const;

    RDouble GetSubsonicInletTotalTemperature() const;

    RDouble GetSubsonicOutletTotalPressure() const;

    RDouble GetSubsonicOutletTotalTemperature() const;

    int GetSubsonicDirectionType() const;

    int GetNumberOfSpecies() const;
    
private:
	int intervalStepForce;

	RDouble wallTemperature;

    //! The number of N-S equations without the species transport equations.
    int nNSEquation;
    //! The number of the equations including the species transport equations.
    int nLaminar;
    //! The flag of the gas model, 1 indicates the chemical non-equilibrium gas, 0 indicates the perfect gas.
    int nChemical;
    //! The flag of thermodynamic temperature model, 1, 2 and 3 denote one/two/three temperature model for chemical non-equilibrium gas.
    int nTemperatureModel;

    //! For incompressible flow (mach < 0.3), using precondition method or not.
	int ifLowSpeedPrecon;

    //! Angle of attacked.
    double AoA;

    //! Angle of slide.
    RDouble angleSlide;

    //! Entropy fix (correction) method.\n
    //!   -# 1: direct fix, which limits the minimum eigenvalue directly.\n
    //!   -# 2: multi-dimensional fix, which is derived from structured solver and now is only valid for struct solver.\n
    //!   -# 3: Harten type, which is default used.
    int roeEntropyFixMethod;

    //! Entropy coefficient scale, default is 1.0.
    //! It is used to scale the default RoeEntropyFixCoef1 and RoeEntropyFixCoef2.
    RDouble roeEntropyScale;

    //! Entropy fix for Acoustic eigenvalue.
	RDouble RoeEntropyFixCoef1;

    //! Entropy fix for convective eigenvalue.
	RDouble RoeEntropyFixCoef2;

    RDouble prandtlLaminar;

    RDouble prandtlTurbulence;

    RDouble oPrandtlLaminar;

    RDouble oPrandtlTurbulence;

    int nChemicalRadius;

    int nChemicalSource;

    int timeIntegration;

    RDouble nonDimensionalSutherlandTemperature;

    int mgProlongationType;
    RDouble massFractionO2;

    RDouble massFractionN2;

    RDouble * primitiveVarFarfield;

    RDouble subsonicInletTotalPressure;

    RDouble subsonicInletTotalTemperature;

    RDouble subsonicOutletTotalPressure;

    RDouble subsonicOutletTotalTemperature;

    int subsonicDirectionType;

    int numberOfSpecies; //! The number of species.
};

#include "Param_NSSolver.hxx"
}