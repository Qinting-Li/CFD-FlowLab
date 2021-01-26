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
//! @file      Gas.h
//! @brief     Explain this file briefly.
//! @author    Wan YunBo.

#pragma once
#include "Precision.h"
#include "DataContainer.h"
#include <string>

using namespace std;

namespace PHSPACE
{
namespace GAS_SPACE
{
//! General gas constant,the unit is J/(Mol*K).
const RDouble rjmk = 8.31434;

class FluidParameter
{
private:
    RDouble machnumber;
    RDouble density, pressure, temperature;
    RDouble sound_velocity, enthalpy;       //! The speed of sound, enthapy.
    RDouble velocity, velocity_x, velocity_y, velocity_z;
    RDouble viscosity;
    RDouble average_molecular_weight;
    RDouble reynoldsnumber;
    RDouble average_gamma;
    RDouble average_Rg;
public:
    RDouble GetDensity() const { return density; }
    void SetDensity(RDouble density) { this->density = density; }

    RDouble GetPressure() const { return pressure; }
    void SetPressure(RDouble pressure) { this->pressure = pressure; }

    RDouble GetTemperature() const { return temperature; }
    void SetTemperature(RDouble temperature) { this->temperature = temperature; }

    RDouble GetVelocity() const { return velocity; }
    void SetVelocity(RDouble velocity) { this->velocity = velocity; }

    RDouble GetMachNumber() const { return machnumber; }
    void SetMachNumber(RDouble machnumber) { this->machnumber = machnumber; }

    RDouble GetViscosity() const { return viscosity; }
    void SetViscosity(RDouble viscosity) { this->viscosity = viscosity; }

    RDouble GetSoundVelocity() const { return sound_velocity; }
    void SetSoundVelocity(RDouble sound_velocity) { this->sound_velocity = sound_velocity; }

    RDouble GetAverageMolecularWeight() const { return average_molecular_weight; }
    void SetAverageMolecularWeight(RDouble average_molecular_weight) { this->average_molecular_weight = average_molecular_weight; }

    RDouble GetReynoldsNumber() const { return reynoldsnumber; }
    void SetReynoldsNumber(RDouble reynoldsnumber) { this->reynoldsnumber = reynoldsnumber; }

    RDouble GetAverageSpecificHeatRatio() const { return average_gamma; }
    void SetAverageSpecificHeatRatio(RDouble average_gamma) { this->average_gamma = average_gamma; }

    RDouble GetAverageGeneralGasConstant() const { return average_Rg; }
    void SetAverageGeneralGasConstant(RDouble average_Rg) { this->average_Rg = average_Rg; }
};

class Gas
{
public:
    Gas(){};
    ~Gas(){};
public:

public:
    virtual int GetSpeciesEquationNumber() const {return 0;}
    virtual int GetNitrogenIndex() const {return 0;}
    virtual int GetElectronIndex() const {return 0;}
    virtual void SetSpeciesEquationNumber(int ns) {}
    virtual void SetNitrogenIndex(int nn) {}
    virtual void SetElectronIndex(int ne) {}
    virtual RDouble GetElectronMolecularWeight() const {return 0;}
    virtual void SetElectronMolecularWeight(RDouble me) {}
    virtual RDouble GetDimensionalElectronMolecularWeight() const {return 0;}
    virtual int GetNSEquationNumber() const{return 0;}; 
    virtual int GetLaminarNumber() const{return 0;};
    virtual int GetTemperatureModel() const {return 0;};
    virtual int GetChemicalType() const {return 0;};

    public:
    virtual RDouble GetTemperatureMin(){return 0;};
    virtual RDouble GetTemperatureMax(){return 0;};
    public:
    virtual void AllocateWorkingSpace(){};
    virtual void DeAllocateWorkingSpace(){};
    virtual void InitCommonParameter(){};
    virtual void InitOtherParameter(){};
    virtual void ComputeReferenceParameter(){};
    virtual void ComputeReferenceTsuth(){};

    virtual void ComputeReferenceViscosityDimensionalChemical(){};
    virtual void ComputeReferenceSpecificHeatRatio(){};
    virtual void ComputeReferenceGeneralGasConstant(){};

    virtual void ComputeReferenceSpecificHeatRatioWithChemical(){};
    virtual void ComputeReferenceViscosityDimensionalNoChemical(){};
    virtual void ComputeReferenceViscosityDimensional(){};
    virtual void ComputeReferenceGasInformation(){};
    virtual void ComputeReferenceSoundVelocity(){};
    virtual void ComputeReferenceMolecularInformation(){};
    virtual void ComputeReferenceMolecularInformationWithChemical(){};
    virtual void ComputeReferenceMolecularInformationWithoutChemical(){};
    virtual void ComputeReferenceVelocity(){};
    virtual void ComputeReferenceReynoldsNumber(){};
    virtual void ComputeProperReynoldsNumberForGrid(){};
    virtual void ComputeOtherProperParameterForGrid(){};
    virtual RDouble ComputeReferenceReynoldsNumberFromTotalPT( void ){return 0;};
    virtual void ComputeDensityWithReynoldsNumber(){};
    virtual void ComputePressureInGasStateEquation(){};
    virtual void InitReferenceParameter(){};
    virtual void InitNumberOfSpecies(){};
    virtual void ComputeCoefficientOfStateEquation(){};

    //! Obtain the T, p, rho and c by the height (Rg=287.053, gamma=1.4 are used). However, only the T and p is useful.
    //! The other two variables will be recomputed to be compatible with the condition of chemical reactions existing.
    virtual void SetAirInformationByDataBase(){};

    //! Obtain the T, p, rho and c by the T0, P0, Ma (Rg=287.053, gamma=1.4 are used){}. However, only the T and p is useful.
    //! The other two variables will be recomputed to be compatible with the condition of chemical reactions existing.
    virtual void SetAirInformationByExpData(){};

    //! This function is mainly for the average density of the mixture in chemical reactions.
    //! The value of mixture density has little difference with that of no chemical reaction, which can be ignored.
    virtual void NormalizeAirInformation(){};

    virtual void ComputeReferencePrimitive(){};

    virtual void CreateAllDataClass(int numberOfSpecies, int numberOfReaction){};
    virtual void Read(fstream &file){};
    virtual void Read(DataContainer *cdata){};
    virtual void Write(DataContainer *cdata){};
public:
    virtual void CompressData  (DataContainer *&cdata){};
    virtual void DecompressData(DataContainer *cdata){};
    virtual RDouble * GetLaminarSchmidtNumberReciprocal  () const {return 0;};
    virtual RDouble * GetTurbulentSchmidtNumberReciprocal() const {return 0;};
public:
    //! To obtain the laminar Schmidt number of each component.
    virtual void GetLaminarSchmidtNumber  (RDouble *Scl) const {};
    //! To obtain the turbulent Schmidt number of each component.
    virtual void GetTurbulentSchmidtNumber(RDouble *Sct) const {};
    //! To obtain the molecular weight of each component.
    virtual RDouble * GetMolecularWeight() const {return 0;}

    virtual int  GetNumberOfSpecies() const{return 0;};
    virtual void SetNumberOfSpecies(int numberOfSpecies) {}
    virtual void InitGlobalParameterOfNSEquation(){};
    virtual void ReadGasModel(fstream &file){};
    //! To Generate the gas model.
    //! Param[in]: chemicalModel indicates the name of the chemical model.
    //! Param[in]: ns indicates the number of species.
    //! Param[in]: nr indicates the number of reactions.
    virtual void CreateGasModel(string chemicalModel, int ns, int nr){};

    virtual void InitGasModel(){};
    //! To import the data of gas model from file.
    virtual void ReadGasModel(){};

    virtual void GetTemperatureRangeIndex(const RDouble &temperature_dimensional, int &indexOfTemperatureRange){};
    virtual RDouble GetCoefficientOfStateEquation(){return 0;};
    virtual string * GetNameOfSpecies() const{return 0;};

    virtual void GetElectronMolecularWeightAndMassFraction(RDouble *massFraction, RDouble &ce, RDouble &Me){};

    //! Compute pressure when the total internal energy e is known.
    //! @param[in ]: primitiveVars is the array that stores the primitive variables.
    //! @param[in ]: gama is the specific heat ratio.
    //! @param[in ]: internalEnergy denotes the total internal energy that equals to the density multiplies the specific internal energy.
    //! @param[in ]: temperature denotes the translation-rotation temperature in the one-temperature model.
    //! @param[out]: pressure is the returned value that denotes the pressure of the mixed gas.
    virtual void GetPressure(RDouble *primitiveVars, RDouble gama, RDouble internalEnergy, RDouble temperature, RDouble &pressure){};

    //! To obtain the non-dimensional enthalpy of each species in accordance with the temperature of mixture gas.
    //! @param[in ]: temperature denotes the non-dimensional temperature of the mixture gas.
    //! @param[out]: speciesEnthaly saves the non-dimensional static enthalpy of each species.
    virtual void ComputeSpeciesEnthalpy(RDouble temperature, RDouble *speciesEnthaly){};

    //! To obtain the dimensional total enthalpy of each species in accordance with the temperature of mixture gas.
    //! @param[in ]: temperature denotes the dimensional temperature of the mixture gas.
    //! @param[out]: speciesEnthaly saves the dimensional total enthalpy of each species.
    virtual void ComputeSpeciesEnthalpyDimensional(RDouble temperature, RDouble *speciesEnthaly){};

    //! To obtain the non-dimensional molecular weight of the mixture gas in accordance with the mass fraction of each species.
    //! @param[in ]: massFraction is an array saving the mass fraction of each species.
    virtual RDouble GetMixedGasMolecularWeight(RDouble *massFraction){return 0;};

    //! To obtain the non-dimensional molecular weight of the mixture gas eliminating electron.
    //! @param[in ]: massFraction is an array saving the mass fraction of each species.
    //! @param[out]: ceDivideMe is the value of the electron mass fraction dividing the electron molecular weight.
    virtual RDouble GetMixedGasMolecularWeightReciprocalWithoutElectron(RDouble *massFraction, RDouble &ceDivideMe){return 0;};

    //! To obtain the non-dimensional universal gas constant called R = 8.314J/(mol*K).
    virtual RDouble GetUniversalGasConstant(){return 0;};

    //! To obtain the static enthalpy of the mixed gas.
    //! @param[in ]: massFraction is an array saving the mass fraction of each species.
    //! @param[in ]: temperature denotes the non-dimensional temperature of the mixture gas.
    //! @param[out]: return the non-dimensional value of static enthalpy of mixed gas.
    virtual RDouble GetMixtureGasEnthalpy(RDouble *massFraction, RDouble temperature){return 0;};

    //! To obtain the non-dimensional constant pressure specific heat called cp according to the mass fraction of each species.
    //! @param[in ]: massFraction is an array saving the mass fraction of each species.
    //! @param[in ]: temperature is the non-dimensional temperature of the mixture gas.
    //! @param[out]: mixedcp denotes the non-dimensional constant pressure specific heat of the mixture gas.
    virtual void ComputeConstantPressureSpecificHeatByMassFraction(RDouble *massFraction, RDouble temperature, RDouble &mixedcp){};

    //! To obtain the non-dimensional constant pressure specific heat called cp of each species by fitting formula.
    //! @param[in ]: temperature is the non-dimensional temperature of the mixture gas.
    //! @param[out]: cpSpecies is an array saving the non-dimensional constant pressure specific heat of each species.
    virtual void ComputeSpeciesConstantPressureSpecificHeat(RDouble temperature, RDouble *cpSpecies){};

    //! To obtain the dimensional constant pressure specific heat called cp[J/(kg*K)] of each species by fitting formula.
    //! @param[in ]: temperature is the non-dimensional temperature of the mixture gas.
    //! @param[out]: cpSpeciesDimensional is an array saving the dimensional constant pressure specific heat of each species.
    virtual void ComputeSpeciesConstantPressureSpecificHeatDimensional(RDouble temperature, RDouble *cpSpeciesDimensional){};

    //! To obtain the non-dimensional constant volume specific heat called cv according to the mass fraction of each species.
    //! @param[in ]: massFraction is an array saving the mass fraction of each species.
    //! @param[in ]: temperature is the non-dimensional temperature of the mixture gas.
    //! @param[out]: mixedGasCv denotes the non-dimensional constant volume specific heat of the mixture gas.
    virtual void ComputeConstantVolumeSpecificHeatByMassFraction(RDouble *massFraction, RDouble temperature, RDouble &mixedGasCv){};

    //! To obtain the dimensional constant volume specific heat called cv[J/(kg*K)] of each species by fitting formula.
    //! @param[in ]: temperature is the non-dimensional temperature of the mixture gas.
    //! @param[out]: cvSpeciesDimensional is an array saving the dimensional constant volume specific heat of each species.
    virtual void ComputeSpeciesConstantVolumeSpecificHeatDimensional(RDouble temperature, RDouble *cvSpeciesDimensional){};

    //! To obtain the non-dimensional constant volume specific heat called cv of each species by fitting formula.
    //! @param[in ]: temperature is the non-dimensional temperature of the mixture gas.
    //! @param[out]: cvSpecies is an array saving the non-dimensional constant volume specific heat of each species.
    virtual void ComputeSpeciesConstantVolumeSpecificHeat(RDouble temperature, RDouble *cvSpecies){};

    //! To obtain the non-dimensional viscosity of each species according to the Blotter fitting formula.
    //! @param[in ]: temperature is the non-dimensional temperature of the mixture gas.
    //! @param[out]: viscosityOfSpecies is an array saving the non-dimensional viscosity of each species.
    virtual void ComputeOneTemperatureModelSpeciesViscosity(RDouble temperature, RDouble *viscosityOfSpecies){};

    //! To obtain the dimensional viscosity of each species according to the Blotter fitting formula.
    //! @param[in ]: temperature is the non-dimensional temperature of the mixture gas.
    //! @param[out]: viscosityOfSpeciesDimensional is an array saving the dimensional viscosity of each species [N*s/m2]
    virtual void ComputeOneTemperatureModelSpeciesViscosityDimensional(RDouble temperature, RDouble *viscosityOfSpeciesDimensional){};

    virtual void ComputeMixtureByMassFraction(RDouble *massFraction, RDouble *mixtureOfSpecies, RDouble &mixture){};
    virtual void ComputeMixtureByPrimitive(RDouble *prim, RDouble *mixtureOfSpecies, RDouble &mixture){};
    virtual void ComputeMixtureCoefficientByWilkeFormula(RDouble *moleFractionOfSpecies, RDouble *chem_var, RDouble *chem_phi){};

    virtual void ComputeMixtureByWilkeFormula(RDouble *moleFractionOfSpecies, RDouble *mixtureOfSpecies, RDouble *phiOfSpecies, RDouble &mixture){};

    //! To obtain the mole fraction of each species according to the mass fraction of each species.
    //! @param[in ]: massFractionOfSpecies is an array of saving the mass fraction of each species.
    //! @param[out]: moleFractionOfSpecies is an array of saving the mole fraction of each species.
    virtual void ComputeMoleFractionByMassFraction(RDouble *massFractionOfSpecies, RDouble *moleFractionOfSpecies){};

    //! To compute the non-dimensional temperature of the mixture gas in accordance with primary variables saving in the array of prim.
    //! @param[in ]: primitiveVariables is an array of saving the primitive variables.
    //! @param[out]: temperature is the non-dimensional value of temperature.
    virtual void ComputeTemperature(RDouble *primitiveVariables, RDouble &temperature){};

    //! To compute the non-dimensional temperature of the mixture gas in accordance with primary variables saving in the array of prim.
    //! @param[in ]: primitiveVariables is an array of saving the primitive variables.
    //! @param[out]: transRotationTemperature is the non-dimensional value of translation-rotation temperature.
    //! @param[out]: vibrationTemperature is the non-dimensional value of vibration temperature.
    //! @param[out]: electronTemperature is the non-dimensional value of electron temperature.
    virtual void GetTemperature(RDouble *primitiveVariables, RDouble &transRotationTemperature, RDouble &vibrationTemperature, RDouble &electronTemperature){};

    //! To obtain the mass fraction of the last species with the label ns-1.
    //! @param[in ]: primitiveVariables is an array of saving the primitive variables.
    virtual void NormalizePrimitive(RDouble *primitiveVariables){};

    //! Compute the molecular weight reciprocal of the mixed gas.
    //! @param[in ]: primitiveVariables is an array of saving the primitive variables.
    //! @param[out]: massReciprocal denotes the reciprocal of the molecular weight.
    virtual void ComputeMolecularWeightReciprocal(RDouble *primitiveVariables, RDouble &massReciprocal){};

    virtual void ComputeMolecularWeightReciprocalDimensional(RDouble *prim, RDouble &omavDimensional){};

    //! Compute the specific heat at constant pressure using the Eucken method.
    //! @param[out]: cpSpecies is an array to saving the specific heat at constant pressure of each species.
    virtual void ComputeSpeciesConstantPressureSpecificHeatByEuckenFormula(RDouble *cpSpecies){};

    //! Compute the static enthalpy with the primitive variables.
    //! @param[in ]: primitiveVariables is an array of saving the primitive variables.
    //! @param[in ]: gama is the specific heat ratio.
    //! @param[out]: enthalpy denotes the static enthalpy.
    //! @param[in ]: temperatures stores the translation-rotation temperature, vibration temperature and the electron temperature.
    virtual void ComputeEnthalpyByPrimitive(RDouble *primitiveVariables, RDouble gama, RDouble &enthalpy, RDouble* temperatures){};

    //! To obtain the total energy E = e + 1/2*V2.
    //! @param[in ]: primitiveVariables is an array of saving the primitive variables, gama is the specific heat ratio of the mixture gas.
    //! @param[in ]: vibrationTemperature is non-dimensional value of vibration temperature.
    //! @param[in ]: electronTemperature is non-dimensional value of electron  temperature.
    //! @param[OUT]: totalEnergy denotes the total internal energy of the mixed gas.
    virtual void ComputeInternalEnergy(RDouble *primitiveVariables, RDouble gama, RDouble vibrationTemperature, RDouble electronTemperature, RDouble &totalEnergy){};

    //! To compute the total enthalpy and the difference of the static enthalpy.
    //! @param[in ]: primitiveVariables is an array of saving the primitive variables.
    //! @param[in ]: gama is the specific heat ratio.
    //! @param[in ]: deltaQ is an array to saving the differences of the primitive variables.
    //! @param[out]: deltaEnthalpy denotes the difference of static enthalpy.
    //! @param[out]: totalEnthalpy denotes the total enthalpy.
    virtual void ComputeDHAndTotalEnthalpy(RDouble *primitiveVariables, RDouble &gama, RDouble *deltaQ, RDouble &deltaEnthalpy, RDouble &totalEnthalpy){};

    //! This function is used for perfect gas and single temperature model.
    //! To compute the total enthalpy and the variable dh in the function MXDQ_STD(). dh=b2, b2 denotes the coefficient \n
    //! of the vector M*dQ which can be referred to the forumla (A.7) and (A.8) in the appendix A of the PHengLEI Theory manual.
    //! @param[in ]: primitiveVariables is an array of saving the primitive variables.
    //! @param[in ]: gama is the specific heat ratio.
    //! @param[in ]: deltaQ is an array to saving the differences of the primitive variables.
    //! @param[out]: totalEnthalpy denotes the total enthalpy.
    //! @param[out]: deltaEnthalpy denotes the difference of static enthalpy.
    virtual void ComputeTotalEnthalpyAndDH(RDouble *primitiveVariables, RDouble &gama, RDouble *deltaQ, RDouble &totalEnthalpy, RDouble &deltaEnthalpy){};

    //! To obtain the non-dimensional viscosity of the mixture gas according to the primary variables.
    //! @param[in ]: primitiveVariables is an array of saving the primitive variables.
    //! @param[in ]: transRotationTemperature is the non-dimensional translation-rotation temperature of the mixture gas.
    //! @param[in ]: electronTemperature is the non-dimensional electron temperature of the mixture gas.
    //! @param[in ]: suthTemperature denotes the temperature value in the Sutherland formula.
    //! @param[out]: viscosity is the non-dimensional viscosity of the mixture gas [N*s/m2].
    virtual void ComputeViscosityByPrimitive(RDouble *primitiveVariables, RDouble &transRotationTemperature, RDouble &electronTemperature, RDouble &suthTemperature, RDouble &viscosity){};

    //£¡To Compute the viscosity, heat conductivity and mass diffusion coefficient.
    //! @param[in ]: primitiveVariables is an array of saving the primitive variables.
    //! @param[in ]: temperature is the array of temperature for each temperature model.
    //! @param[out]: viscosity is the non-dimensional viscosity of the mixture gas.
    //! @param[out]: heatConductivity denotes the array of heat conductivity for each mode.
    //! @param[out]: rhoDs denotes the array of mass diffusion coefficients for each species.
    virtual void ComputeTransportCoefficients(RDouble *primitiveVariables, RDouble *temperature, RDouble visTurbulence, RDouble &viscosity, RDouble *heatConductivity, RDouble *rhoDs){};

    //! To obtain the non-dimensional viscosity of the mixture gas according to the mass fraction of each species.
    //! @param[in ]: density denotes the non-dimensional value of density of the mixed gas.
    //! @param[in ]: massFractionOfSpecies is an array of saving the mass fraction of each species.
    //! @param[in ]: transRotationTemperature is the non-dimensional translation-rotation temperature of the mixture gas.
    //! @param[in ]: electronTemperature is the non-dimensional electron temperature of the mixture gas.
    //! @param[in ]: suthTemperature denotes the temperature value in the Sutherland formula.
    //! @param[out]: viscosity is the non-dimensional viscosity of the mixture gas [N*s/m2].  
    virtual void ComputeViscosity(RDouble density, RDouble *massFractionOfSpecies, RDouble &transRotationTemperature, RDouble &electronTemperature, RDouble &suthTemperature, RDouble &viscosity){};

    //! To compute the heat conductivity of each species using the Eucken formula.
    //! @param[in ]: viscosity is array of saving the species viscosity.
    //! @param[in ]: cps is the array of saving the species values of specific heat at constant pressure.
    //! @param[out]: conductivity is an array of saving the heat conductivity of each species.
    virtual void ComputeSpeciesHeatConductivityByEuckenFormula(RDouble *viscosity, RDouble *cps, RDouble *conductivity){};

    //! To compute heat conductivity of mixture gas using Eucken formula and Wassilewa relation.
    //! @param[in ]: primitiveVariables is an array of saving the primitive variables.
    //! @param[in ]: temperature is the non-dimensional temperature of the mixture gas.
    //! The function returns the heat conductivity of mixture gas.
    virtual RDouble ComputeMixtureGasHeatConductivityByWassilewaFormula(RDouble *primitiveVariables, RDouble temperature){return 0;};

    //! To obtain the partition function values based on the targets, e.g conductivity or viscosity.
    //! @param[in ]: target is array of viscosity or heat conductivity of each species.
    //! @param[in ]: speciesMass is array of molecular weight.
    //! @param[in ]: moleFraction is array of mole fraction.
    //! @param[out]: phi is an array of saving the partition function value of each species.
    virtual void GetPartitionFunctionValues(RDouble *target, RDouble *speciesMass, RDouble *moleFraction, RDouble *phi){};

    //! Compute the source terms of chemical reaction equations.
    //! @param[in ]: primitiveVariables stores the primitive variables of the specified control cells.
    //! @param[in ]: temperatures stores the temperature values of the specified control cells.
    //! @param[in ]: cellVolume is an array that stores the volume values of the specified control cells.
    //! @param[out]: chemicalSourceTerms returns the source terms of chemical reaction equations in the specified control cells.
    //! @param[in ]: nCellNumber denotes the number of the specified control cells.
    virtual void ChemicalSource(RDouble **primitiveVariables, RDouble **temperatures, RDouble *cellVolume, RDouble **chemicalSourceTerms, int nCellNumber){};

    //! Compute the source terms of vibration-electron energy equations.
    //! @param[in ]: primitiveVariables stores the primitive variables of the specified control cells.
    //! @param[in ]: temperatures stores the temperature values of the specified control cells.
    //! @param[in ]: cellVolume is an array that stores the volume values of the specified control cells.
    //! @param[in ]: chemicalSourceTerms stores the source terms of chemical reaction equations in the specified control cells.
    //! @param[out]: energySourceTerms returns the source terms of vibration-electron energy equations in the specified control cells.
    //! @param[in ]: nCellNumber denotes the number of the specified control cells.
    virtual void EnergySource(RDouble **primitiveVariables, RDouble **temperatures, RDouble *cellVolume, RDouble **chemicalSourceTerms, RDouble **energySourceTerms, int nCellNumber){};

    //! Compute the source term of vibration and electron energy equation.
    //! @param[in ]: primitiveVariables is an array of saving the primitive variables.
    //! @param[in ]: omega is the source term of chemical reaction equation.
    //! @param[in ]: temperature is non-dimensional value of temperatures.
    //! @param[out]: return the dimensional value of source term of vibration and electron energy equation.
    virtual RDouble ComputeVibrationAndElectronEnergySourceTerm(RDouble *primitiveVariables, RDouble *omega, RDouble *temperature){return 0;};

    //! Compute the source term of vibration energy equation.
    //! @param[in ]: primitiveVariables is an array of saving the primitive variables.
    //! @param[in ]: omega is the source term of chemical reaction equation.
    //! @param[in ]: temperature is non-dimensional value of temperatures.
    //! @param[out]: return the dimensional value of source term of vibration energy equation.
    virtual RDouble ComputeVibrationEnergySourceTerm(RDouble *primitiveVariables, RDouble *omega, RDouble *temperature){return 0.0;};

    //! Compute the source term of electron energy equation.
    //! @param[in ]: primitiveVariables is an array of saving the primitive variables.
    //! @param[in ]: omega is the source term of chemical reaction equation.
    //! @param[in ]: temperature is non-dimensional value of temperatures.
    //! @param[out]: return the dimensional value of source term of vibration energy equation.
    virtual RDouble ComputeElectronEnergySourceTerm(RDouble *primitiveVariables, RDouble *omega, RDouble *temperature){return 0.0;};

    virtual void ChemicalSpectrumRadius(RDouble **q, RDouble **t, RDouble *vol, RDouble **src, int nlen){};

    //! To obtain the approximating polynomial coefficients from the polynomial curve fits in five temperature ranges. \n
    //! The approximating polynomial is input for evaluating the species specific heats and enthalpies, and the temperature\n
    //! boundary is smoothed by linearly averaging the polynomial coefficients.
    //! @param[in ]: nSpecies denotes the species index in the collection.
    //! @param[in ]: temperature is the dimensional temperature value.
    //! @param[out]: polynomialCoef is an array of saving the polynomial coefficients.
    virtual void GetLinearlyAveragedPolynomialCoefficients(int nSpecies, RDouble temperature, RDouble *polynomialCoef){};

    //! Compute chemical source term and its Jacobian.
    virtual void ComputeChemicalSourceDerivatives(RDouble **primitiveVariables, RDouble **temperatures, RDouble *cellVolume, RDouble **chemicalSourceTerms, RDouble ***sourceDerivatives, int nCellNumber){};
    virtual void ComputeChemicalSourceJacobian(RDouble *primitiveVariables, RDouble *temperatures, RDouble cellVolume, RDouble *chemicalSourceTerms, RDouble **sourceDerivatives){};
    virtual void ComputeChemicalSourceAndJacobian(RDouble **primitiveVariables, RDouble **temperatures, RDouble *cellVolume, RDouble **chemicalSourceTerms, RDouble ***sourceDerivatives, int nCellNumber){};

    //! Compute the partial derivatives of source terms of chemical reaction equations without electron.
    //! @param[in ]: primitiveVariables stores the primitive variables of the specified control cells.
    //! @param[in ]: temperatures stores the temperature values of the specified control cells.
    //! @param[in ]: cellVolume is an array that stores the volume values of the specified control cells.
    //! @param[out]: chemicalSourceTerms returns the source terms of chemical reaction equations in the specified control cells.
    //! @param[out]: sourceDerivatives returns the elements of the Jacobian matrix in the source terms of chemical reaction equations in the specified control cells.
    //! @param[in ]: nCellNumber denotes the number of the specified control cells.
    virtual void ComputeChemicalSourceDerivatives1(RDouble **primitiveVariables, RDouble **temperatures, RDouble *cellVolume, RDouble **chemicalSourceTerms, RDouble ***sourceDerivatives, int nCellNumber){};

    //! Compute the partial derivatives of source terms of chemical reaction equations, it is used for the reactions including electron.
    //! @param[in ]: primitiveVariables stores the primitive variables of the specified control cells.
    //! @param[in ]: temperatures stores the temperature values of the specified control cells.
    //! @param[in ]: cellVolume is an array that stores the volume values of the specified control cells.
    //! @param[out]: chemicalSourceTerms returns the source terms of chemical reaction equations in the specified control cells.
    //! @param[out]: sourceDerivatives returns the elements of the Jacobian matrix in the source terms of chemical reaction equations in the specified control cells.
    //! @param[in ]: nCellNumber denotes the number of the specified control cells.
    virtual void ComputeChemicalSourceDerivatives2(RDouble **primitiveVariables, RDouble **temperatures, RDouble *cellVolume, RDouble **chemicalSourceTerms, RDouble ***sourceDerivatives, int nCellNumber){};

    //! Compute the partial derivatives of source terms of chemical reaction equations. The function is used for two-temperature model.
    //! @param[in ]: primitiveVariables stores the primitive variables of the specified control cells.
    //! @param[in ]: temperatures stores the temperature values of the specified control cells.
    //! @param[in ]: cellVolume is an array that stores the volume values of the specified control cells.
    //! @param[in ]: chemicalSourceTerms stores the source terms of chemical reaction equations in the specified control cells.
    //! @param[in ]: chemicalSourceDerivatives stores the partial derivatives of the source terms of chemical reaction equations in the specified control cells.
    //! @param[out]: energySourceDerivatives returns the partial derivatives of the source terms of energy equations in the specified control cells.
    //! @param[in ]: nCellNumber denotes the number of the specified control cells.
    virtual void ComputeEnergySourceDerivatives1(RDouble **primitiveVariables, RDouble **temperatures, RDouble *cellVolume, RDouble **chemicalSourceTerms,
        RDouble ***chemicalSourceDerivatives, RDouble ***energySourceDerivatives, int nCellNumber){};

    //! Compute the partial derivatives of source terms of vibration-electron energy equations. The function is used for two-temperature model.
    //! @param[in ]: primitiveVariables stores the primitive variables of the specified control cells.
    //! @param[in ]: temperatures stores the temperature values of the specified control cells.
    //! @param[in ]: cellVolume is an array that stores the volume values of the specified control cells.
    //! @param[in ]: chemicalSourceTerms stores the source terms of chemical reaction equations in the specified control cells.
    //! @param[in ]: chemicalSourceDerivatives stores the partial derivatives of the source terms of chemical reaction equations in the specified control cells.
    //! @param[out]: energySourceDerivatives returns the diagonal element of partial derivatives to the source terms of energy equations in the specified control cells.
    //! @param[in ]: nCellNumber denotes the number of the specified control cells.
    virtual void ComputeEnergySourceDiagonalDerivatives1(RDouble **primitiveVariables, RDouble **temperatures, RDouble *cellVolume, RDouble **chemicalSourceTerms,
        RDouble ***chemicalSourceDerivatives, RDouble **energySourceDerivatives, int nCellNumber){};

    //! Compute the partial derivatives of source terms of vibration-electron energy equations. The function is used for three-temperature model.
    //! @param[in ]: primitiveVariables stores the primitive variables of the specified control cells.
    //! @param[in ]: temperatures stores the temperature values of the specified control cells.
    //! @param[in ]: cellVolume is an array that stores the volume values of the specified control cells.
    //! @param[in ]: chemicalSourceTerms stores the source terms of chemical reaction equations in the specified control cells.
    //! @param[in ]: chemicalSourceDerivatives stores the partial derivatives of the source terms of chemical reaction equations in the specified control cells.
    //! @param[out]: energySourceDerivatives returns the diagonal element of partial derivatives to the source terms of energy equations in the specified control cells.
    //! @param[in ]: nCellNumber denotes the number of the specified control cells.
    virtual void ComputeEnergySourceDiagonalDerivatives2(RDouble **primitiveVariables, RDouble **temperatures, RDouble *cellVolume, RDouble **chemicalSourceTerms,
            RDouble ***chemicalSourceDerivatives, RDouble **energySourceDerivatives, int nCellNumber){};

public:
    virtual void Primitive2Conservative(RDouble *prim, RDouble gama, RDouble Tv, RDouble Te, RDouble *q){};
    virtual void Conservative2Primitive(RDouble *q, RDouble gama, RDouble *prim, RDouble temperature){};
    virtual void GetSpecificHeatRatio(RDouble *prim, RDouble &gama){};

    //! Compute the specific heat ratio and the temperatures.
    //! @param[in ]: primitiveVariables stores the primitive variables of the specified control cells.
    //! @param[out]: gama denotes the specific heat ratio of the mixed gas.
    //! @param[out]: transRotationTemperature denotes the translation-rotation temperature of the mixed gas.
    //! @param[out]: vibrationTemperature denotes the vibration temperature of the mixed gas.
    //! @param[out]: electronTemperature denotes the electron temperature of the mixed gas.
    virtual void GetSpecificHeatRatioAndTemperatute(RDouble *primitiveVariables, RDouble &gama, RDouble &transRotationTemperature, RDouble &vibrationTemperature, RDouble &electronTemperature){};

    virtual void ComputeDensityDiffusionAndKCP(RDouble temperature, RDouble *primface, RDouble mul, RDouble mut, RDouble *rho_ds_face, RDouble *hintSpeciesOfFace, RDouble &kcp){};
    virtual void GetAirInfo(const RDouble &height, RDouble &temperature, RDouble &pressure, RDouble &density, RDouble &soundspeed){};

    //! To obtain the coefficients of computing the M*dQ in the function MXDQ_STD_LP(), added by LiPeng on Jan. 10, 2019.
    //! @param[in ]: primitiveVars is the array that stores the primitive variables of the current cell.
    //! @param[out]: alpha, beta and speciesBeta are referred to the formula (A.5) in the appendix A of the PHengLEI Theory manual.
    virtual void ComputeCoefficientInMXDQ(RDouble *primitiveVars, RDouble &alpha, RDouble &beta, RDouble *speciesBeta){};

    //! To obtain the coefficients of computing the M*dQ in the function MVXDQ_STD_LP(), added by LiPeng on Jan. 21, 2019.
    //! K=visl/Prl + vist/Prt denotes the thermal conductivity coefficient, rD=visl/Scl + vist/Sct is the species diffusion coefficient,\n
    //! Ds is an array of saving the diffusion coefficient of each species.
    virtual void ComputeCoefficientInMVXDQ(RDouble *prim, RDouble K, RDouble rD, RDouble *theta_s, RDouble *phi_s, RDouble *Ds){};

public:
    //! Compute the mass fraction of the Oxygen and Nitrogen which is used for the fully catalytic wall.
    //! @param[in ]: fs is array of mass fractions of the components in the multi-species.
    //! @param[out]: co2 stores the mass fraction of Oxygen.
    //! @param[out]: cn2 stores the mass fraction of Nitrogen.
    virtual void GetOxygenAndNitrogenMassFraction(RDouble *fs, RDouble &co2, RDouble &cn2){};

    //! To find the chemical element in the species, if found, return the coefficient of the element.
    //! @param[in ]: species_name denotes the name of the species.
    //! @param[out]: element_name denotes the name of the specified element.
    //! @param[out]: coef is the stoichiometric coefficient of the element.
    virtual bool FindChemicalElement(string species_name, char element_name, int &coef){return 0;};

    //! To obtain the index of the specified species in the sequence in accordance with its name.
    //! @param[in ]: species_name denotes the name of the species.
    //! @param[out]: the function returns the index of the species if found, otherwise, returns value of -1.
    virtual int GetSpeciesIndex(string species_name){return 0;};

public:     //! The function library of multi-temperature model.

    //! To compute the translation and rotation specific heat at constant volume.
    //! @param[out]: speciesTransRotationCv is an array saving the translation and rotation specific heat at constant volume of each species.
    virtual void ComputeTranslationAndRotationSpecificHeatAtConstantVolume(RDouble *speciesTransRotationCv){};
    
    //! To compute the translation and rotation specific heat at constant volume.
    //! @param[out]: speciesTransRotationCv is an array saving the dimensional value of translation and rotation specific heat at constant volume of each species.
    virtual void ComputeDimensionalTranslationAndRotationSpecificHeatAtConstantVolume(RDouble *speciesTransRotationCv){};
    
    //! To compute the the translation and rotation energy.
    //! @param[in ]: transRotationTemperature is the non-dimensional translation-rotation temperature of the mixture gas.
    //! @param[out]: speciesTransRotationEnergy is an array saving the translation and rotation energy of each species.
    virtual void ComputeTranslationAndRotationEnergy(RDouble transRotationTemperature, RDouble *speciesTransRotationEnergy){};
    
    //! To compute the vibration specific heat at constant volume.
    //! @param[in ]: vibrationTemperature is the non-dimensional vibration temperature of the mixture gas.
    //! @param[out]: speciesVibrationCv is an array saving the vibration specific heat at constant volume of each species.
    virtual void ComputeVibrationSpecificHeatAtConstantVolume(RDouble vibrationTemperature, RDouble *speciesVibrationCv){};
    
    //! To compute the vibration specific heat at constant volume.
    //! @param[in ]: vibrationTemperature is the non-dimensional vibration temperature of the mixture gas.
    //! @param[out]: speciesVibrationCv is an array saving the dimensional value of vibration specific heat at constant volume of each species.
    virtual void ComputeDimensionalVibrationSpecificHeatAtConstantVolume(RDouble vibrationTemperature, RDouble *speciesVibrationCv){};
    
    //! To compute the dimensional value of the vibration energy (J/kg).
    //! @param[in ]: vibrationTemperature is the non-dimensional vibration temperature of the mixture gas.
    //! @param[out]: speciesVibrationEnergy is an array saving the vibration energy of each species.
    virtual void ComputeDimensionalVibrationEnergy(RDouble vibrationTemperature, RDouble *speciesVibrationEnergy){};
    
    //! To compute the vibration energy.
    //! @param[in ]: vibrationTemperature is the non-dimensional vibration temperature of the mixture gas.
    //! @param[out]: speciesVibrationEnergy is an array saving the vibration energy of each species.
    virtual void ComputeVibrationEnergy(RDouble vibrationTemperature, RDouble *speciesVibrationEnergy){};
    
    //! To compute the electron specific heat at constant volume.
    //! @param[in ]: electronTemperature is the non-dimensional electron temperature of the mixture gas.
    //! @param[out]: speciesElectronCv is an array saving the electron specific heat at constant volume of each species.
    virtual void ComputeElectronSpecificHeatAtConstantVolume(RDouble electronTemperature, RDouble *speciesElectronCv){};
    
    //! To compute the electron specific heat at constant volume.
    //! @param[in ]: electronTemperature is the non-dimensional electron temperature of the mixture gas.
    //! @param[out]: speciesElectronCv is an array saving the dimensional value of electron specific heat at constant volume of each species.
    virtual void ComputeDimensionalElectronSpecificHeatAtConstantVolume(RDouble electronTemperature, RDouble *speciesElectronCv){};
    
    //! To compute the dimensional value of electron energy (J/kg).
    //! @param[in ]: electronTemperature is the non-dimensional electron temperature of the mixture gas.
    //! @param[out]: speciesElectronEnergy is an array saving the electron energy of each species.
    virtual void ComputeDimensionalElectronEnergy(RDouble electronTemperature, RDouble *speciesElectronEnergy){};
    
    //! To compute the electron energy.
    //! @param[in ]: electronTemperature is the non-dimensional electron temperature of the mixture gas.
    //! @param[out]: speciesElectronEnergy is an array saving the electron energy of each species.
    virtual void ComputeElectronEnergy(RDouble electronTemperature, RDouble *speciesElectronEnergy){};
    
    //! To compute the enthalpies of species.
    //! @param[in ]: transRotationTemperature is the non-dimensional translation-rotation temperature of the mixture gas.
    //! @param[in ]: vibrationTemperature is the non-dimensional vibration temperature of the mixture gas.
    //! @param[in ]: electronTemperature is the non-dimensional electron temperature of the mixture gas.
    //! @param[out]: speciesEnthalpy is an array saving the enthalpies of each species.
    virtual void ComputeSpeciesEnthalpy(RDouble transRotationTemperature, RDouble vibrationTemperature, RDouble electronTemperature, RDouble *speciesEnthalpy){};
    
    //! To compute the translation-rotation specific heat at constant volume of the mixed gas.
    //! @param[in ]: massFraction is an array to storing the mass fractions of species. 
    //! @param[out]: return the total translation-rotation specific heat at constant volume.
    virtual RDouble GetMixedGasTranslationAndRotationSpecificHeatAtConstantVolume(RDouble *massFraction){return 0;};
    
    //! To compute the translation-rotation specific heat at constant volume of the mixed gas.
    //! @param[in ]: massFraction is an array to storing the mass fractions of species. 
    //! @param[out]: return the dimensional value of total translation-rotation specific heat at constant volume.
    virtual RDouble GetMixedGasDimensionalTranslationAndRotationSpecificHeatAtConstantVolume(RDouble *massFraction){return 0;};
    
    //! To compute the vibration specific heat at constant volume of the mixed gas.
    //! @param[in ]: vibrationTemperature is the non-dimensional vibration temperature of the mixture gas.
    //! @param[in ]: massFraction is an array to storing the mass fractions of species. 
    //! @param[out]: return the total vibration specific heat at constant volume.
    virtual RDouble GetMixedGasVibrationSpecificHeatAtConstantVolume(RDouble vibrationTemperature, RDouble *massFraction){return 0;};
    
    //! To compute the vibration specific heat at constant volume of the mixed gas.
    //! @param[in ]: vibrationTemperature is the non-dimensional vibration temperature of the mixture gas.
    //! @param[in ]: massFraction is an array to storing the mass fractions of species. 
    //! @param[out]: return the dimensional value of total vibration specific heat at constant volume.
    virtual RDouble GetMixedGasDimensionalVibrationSpecificHeatAtConstantVolume(RDouble vibrationTemperature, RDouble *massFraction){return 0;};
    
    //! To compute the electron specific heat at constant volume of the mixed gas.
    //! @param[in ]: electronTemperature is the non-dimensional electron temperature of the mixture gas.
    //! @param[in ]: massFraction is an array to storing the mass fractions of species. 
    //! @param[out]: return the total electron specific heat at constant volume.
    virtual RDouble GetMixedGasElectronSpecificHeatAtConstantVolume(RDouble electronTemperature, RDouble *massFraction){return 0;};
    
    //! To compute the electron specific heat at constant volume of the mixed gas.
    //! @param[in ]: electronTemperature is the non-dimensional electron temperature of the mixture gas.
    //! @param[in ]: massFraction is an array to storing the mass fractions of species. 
    //! @param[out]: return the dimensional value of total electron specific heat at constant volume.
    virtual RDouble GetMixedGasDimensionalElectronSpecificHeatAtConstantVolume(RDouble electronTemperature, RDouble *massFraction){return 0;};
    
    //! To compute the translation-rotation energy of the mixed gas.
    //! @param[in ]: transRotationTemperature is the non-dimensional translation-rotation temperature of the mixture gas.
    //! @param[in ]: massFraction is an array to storing the mass fractions of species. 
    //! @param[out]: return the total translation-rotation energy.
    virtual RDouble GetMixedGasTranslationAndRotationEnergy(RDouble transRotationTemperature, RDouble *massFraction){return 0;};
    
    //! To compute the dimensional value of vibration energy of the mixed gas.
    //! @param[in ]: vibrationTemperature is the non-dimensional vibration temperature of the mixture gas.
    //! @param[in ]: massFraction is an array to storing the mass fractions of species. 
    //! @param[out]: return the total vibration energy.
    virtual RDouble GetMixedGasDimensionalVibrationEnergy(RDouble vibrationTemperature, RDouble *massFraction){return 0;};
    
    //! To compute the vibration energy of the mixed gas.
    //! @param[in ]: vibrationTemperature is the non-dimensional vibration temperature of the mixture gas.
    //! @param[in ]: massFraction is an array to storing the mass fractions of species. 
    //! @param[out]: return the total vibration energy.
    virtual RDouble GetMixedGasVibrationEnergy(RDouble vibrationTemperature, RDouble *massFraction){return 0;};
    
    //! To compute the dimensional value of electron energy of the mixed gas.
    //! @param[in ]: electronTemperature is the non-dimensional electron temperature of the mixture gas.
    //! @param[in ]: massFraction is an array to storing the mass fractions of species. 
    //! @param[out]: return the total electron energy.
    virtual RDouble GetMixedGasDimensionalElectronEnergy(RDouble electronTemperature, RDouble *massFraction){return 0;};
    
    //! To compute the electron energy of the mixed gas.
    //! @param[in ]: electronTemperature is the non-dimensional electron temperature of the mixture gas.
    //! @param[in ]: massFraction is an array to storing the mass fractions of species. 
    //! @param[out]: return the total electron energy.
    virtual RDouble GetMixedGasElectronEnergy(RDouble electronTemperature, RDouble *massFraction){return 0;};
    
    //! To compute the total enthalpy of the mixed gas.
    //! @param[in ]: transRotationTemperature is the non-dimensional translation-rotation temperature of the mixture gas.
    //! @param[in ]: vibrationTemperature is the non-dimensional vibration temperature of the mixture gas.
    //! @param[in ]: electronTemperature is the non-dimensional electron temperature of the mixture gas.
    //! @param[out]: return the total enthalpy of the mixed gas.
    virtual RDouble ComputeMixedGasEnthalpy(RDouble transRotationTemperature, RDouble vibrationTemperature, RDouble electronTemperature, RDouble *massFraction){return 0;};
    
    //! To compute the total enthalpy of the mixed gas.
    //! @param[in ]: primitiveVars indicates the array of primitive variables.
    //! @param[in ]: temperatures indicates the array of temperatures.
    virtual RDouble ComputeMixedGasEnthalpy(RDouble *primitiveVars, RDouble *temperatures){return 0;};
    
    //! To compute the total enthalpy of formation.
    //! @param[in ]: massFraction is an array to storing the mass fractions of species.
    //! @param[out]: return the non-dimensional value of the formation enthalpy.
    virtual RDouble GetTotalFormationEnthalpy(RDouble *massFraction){return 0;};
    
    //! To compute translation-rotation temperature where the pressure is known.
    //! @param[in ]: primitiveVariables is an array of saving the primary variables.
    //! @param[in ]: electronPressure is the pressure of electron.
    //! @param[out]: return the non-dimensional value of translation-rotation temperature.
    virtual RDouble ComputeTranslationAndRotationTemperatureByPrimitiveVariables(RDouble *primitiveVariables, RDouble electronPressure = 0.0){return 0;};
    
    //! To compute translation-rotation temperature.
    //! @param[in ]: transRotationEnergy is the total internal energy of translation and rotation mode.
    //! @param[in ]: massFraction is an array of saving the mass fractions of species.
    //! @param[out]: return the non-dimensional value of translation-rotation temperature.
    virtual RDouble ComputeTranslationAndRotationTemperature(RDouble *massFraction, RDouble transRotationEnergy){return 0;};
    
    //! To compute vibration temperature.
    //! @param[in ]: massFraction is an array of saving the mass fractions of species.
    //! @param[in ]: vibrationEnergy is the total internal energy of vibration mode.
    //! @param[in ]: initTemperature is the initial temperature of the newton iteration method.
    //! @param[out]: return the non-dimensional value of vibration temperature.
    virtual RDouble ComputeVibrationTemperature(RDouble *massFraction, RDouble vibrationEnergy, RDouble initTemperature){return 0;};
    
    //! To compute electron temperature.
    //! @param[in ]: massFraction is an array of saving the mass fractions of species.
    //! @param[in ]: electronEnergy is the total internal energy of electron mode.
    //! @param[in ]: initTemperature is the initial temperature of the newton iteration method.
    //! @param[out]: return the non-dimensional value of electron temperature.
    virtual RDouble ComputeElectronTemperature(RDouble *massFraction, RDouble electronEnergy, RDouble initTemperature){return 0;};
    
    //! To compute vibration-electron temperature.
    //! @param[in ]: massFraction is an array of saving the mass fractions of species.
    //! @param[in ]: vibrationElectronEnergy is the total internal energy of vibration and electron mode.
    //! @param[in ]: initTemperature is the initial temperature of the newton iteration method.
    //! @param[out]: return the non-dimensional value of vibration-electron temperature.
    virtual RDouble ComputeVibrationAndElectronTemperature(RDouble *massFraction, RDouble vibrationElectronEnergy, RDouble initTemperature){return 0;};
    
    //! To compute the temperature of one-temperature model.
    //! @param[in ]: massFraction is an array of mass fractions of species.
    //! @param[in ]: internalEnergy is total internal energy of the mixed gas.
    //! @param[out]: return the non-dimensional value of temperature.
    virtual RDouble ComputeOneTemperatureModelTemperature(RDouble *massFraction, RDouble internalEnergy){return 0;};
    
    //! To compute non-equilibrium temperature of chemical reaction flow for one-temperature model using the curve fitting method.
    //! @param[in ]: massFraction is an array of mass fractions of species.
    //! @param[in ]: internalEnergy is total internal energy of the mixed gas.
    //! @param[out]: return the non-dimensional value of temperature.
    virtual RDouble ComputeNonequilibrumTemperatureViaBisectionMethod(RDouble *massFraction, RDouble internalEnergy){return 0;};
    
    //! To compute non-equilibrium temperature of chemical reaction flow for one-temperature model using the newton method.
    //! @param[in ]: massFraction is an array of mass fractions of species.
    //! @param[in ]: internalEnergy is total internal energy of the mixed gas.
    //! @param[in ]: initTemperature is the initial temperature of the newton iteration method.
    //! @param[out]: return the non-dimensional value of temperature.
    virtual RDouble ComputeNonequilibrumTemperature(RDouble *massFraction, RDouble internalEnergy, RDouble initTemperature){return 0;};
    
    //! To obtain the electron pressure.
    //! @param[in ]: primitiveVariables is an array of saving the primary variables.
    //! @param[in ]: electronTemperature is non-dimensional value of electron temperature.
    //! @param[out]: return the non-dimensional value of electron pressure.
    virtual RDouble GetElectronPressure(RDouble *primitiveVariables, RDouble electronTemperature){return 0;};
    
    //! To compute the heat conductivity of translation and rotation using the Eucken formula.
    //! @param[in ]: viscosity is array of saving the species viscosity.
    //! @param[out]: conductivity is an array of saving the translation-rotation heat conductivity of each species.
    virtual void ComputeSpeciesTranslationAndRotationHeatConductivity(RDouble *viscosity, RDouble *conductivity){};
    
    //! To compute the heat conductivity of vibration using the Eucken formula.
    //! @param[in ]: viscosity is array of saving the species viscosity.
    //! @param[in ]: speciesVibrationCv is the array of saving the vibration specific heat at constant pressure.
    //! @param[out]: conductivity is an array of saving the vibration heat conductivity of each species.
    virtual void ComputeSpeciesVibrationHeatConductivity(RDouble *viscosity, RDouble *speciesVibrationCv, RDouble *conductivity){};
    
    //! To compute the heat conductivity of electron using the Eucken formula.
    //! @param[in ]: viscosity is array of saving the species viscosity.
    //! @param[in ]: speciesElectronCv is the array of saving the electron specific heat at constant pressure.
    //! @param[out]: conductivity is an array of saving the electron heat conductivity of each species.
    virtual void ComputeSpeciesElectronHeatConductivity(RDouble *viscosity, RDouble *speciesElectronCv, RDouble *conductivity){};
    
    //! To compute the translation-rotation heat conductivity of the mixed gas using Wilke method.
    //! @param[in ]: primitiveVariables is an array of saving the primary variables.
    //! @param[in ]: transRotationTemperature denotes the non-dimensional value of translation-rotation temperature.
    //! @param[in ]: electronTemperature denotes the non-dimensional value of electron temperature.
    //! @param[out]: return the translation-rotation heat conductivity of the mixed gas.
    virtual RDouble GetMixedGasTranslationAndRotationHeatConductivity(RDouble *primitiveVariables, RDouble transRotationTemperature, RDouble electronTemperature){return 0;};
    
    //! To compute the vibration heat conductivity of the mixed gas using Wilke method.
    //! @param[in ]: primitiveVariables is an array of saving the primary variables.
    //! @param[in ]: transRotationTemperature denotes the non-dimensional value of translation-rotation temperature.
    //! @param[in ]: vibrationTemperature denotes the non-dimensional value of vibration temperature.
    //! @param[in ]: electronTemperature denotes the non-dimensional value of electron temperature.
    //! @param[out]: return the vibration heat conductivity of the mixed gas.
    virtual RDouble GetMixedGasVibrationHeatConductivity(RDouble *primitiveVariables, RDouble transRotationTemperature, RDouble vibrationTemperature, RDouble electronTemperature){return 0;};
    
    //! To compute the electron heat conductivity of the mixed gas using Wilke method.
    //! @param[in ]: primitiveVariables is an array of saving the primary variables.
    //! @param[in ]: transRotationTemperature denotes the non-dimensional value of translation-rotation temperature.
    //! @param[in ]: electronTemperature denotes the non-dimensional value of electron temperature.
    //! @param[out]: return the electron heat conductivity of the mixed gas.
    virtual RDouble GetMixedGasElectronHeatConductivity(RDouble *primitiveVariables, RDouble transRotationTemperature, RDouble electronTemperature){return 0;};
    
    virtual void ComputeMixedGasTranslationVibrationElectronHeatConductivity(RDouble *primitiveVariables, RDouble transRotationTemperature, RDouble vibrationTemperature, RDouble electronTemperature, RDouble &transRotationConductivity, RDouble &vibrationConductivity, RDouble &electronConductivity){};
    
    //! To obtain the heat conductivity of the mixed gas using Wilke method, the value is identified by the type of the heat conductivity.
    //! @param[in ]: moleFraction is an array of saving the mole fractions of species.
    //! @param[in ]: conductivity is an array of saving the heat conductivity of species. The heat conductivity could be translation-rotation\n
    //! heat conductivity, vibration heat conductivity or electron heat conductivity.
    //! @param[in ]: Phi denotes the array of the partition function.
    //! @param[out]: return the heat conductivity of the mixed gas.
    virtual RDouble GetMixedGasHeatConductivityWithWilkeFormula(RDouble *moleFraction, RDouble *conductivity, RDouble *phi){return 0;};
    
    //! To compute the non-dimensional values of viscosities of species using Gupta or Blotter method. 
    //! @param[in ]: transRotationTemperature denotes the non-dimensional value of translation-rotation temperature.
    //! @param[in ]: electronTemperature denotes the non-dimensional value of electron temperature.
    //! @param[out]: speciesViscosity is an array to storing the non-dimensional values of viscosities of species.
    virtual void ComputeSpeciesViscosityWithCurveFitMethod(RDouble transRotationTemperature, RDouble electronTemperature, RDouble *speciesViscosity){};
    
    //! To compute the dimensional values of viscosities of species using Gupta or Blotter method. 
    //! @param[in ]: transRotationTemperature denotes the non-dimensional value of translation-rotation temperature.
    //! @param[in ]: electronTemperature denotes the non-dimensional value of electron temperature.
    //! @param[out]: speciesViscosity is an array to storing the dimensional values of viscosities of species(Unit: kg/(m*s)).
    virtual void ComputeSpeciesDimensionalViscosityWithCurveFitMethod(RDouble transRotationTemperature, RDouble electronTemperature, RDouble *speciesViscosity){};
    
    //! To compute the non-dimensional values of viscosities of species using Lennard-Jones method. 
    //! @param[in ]: transRotationTemperature denotes the non-dimensional value of translation-rotation temperature.
    //! @param[in ]: electronPressure denotes the non-dimensional value of electron pressure.
    //! @param[out]: speciesViscosity is an array to storing the non-dimensional values of viscosities of species.
    virtual void ComputeSpeciesViscosityWithLennardJonesMethod(RDouble transRotationTemperature, RDouble electronPressure, RDouble *speciesViscosity){};
    
    //! To compute the dimensional values of viscosities of species using using Lennard-Jones method. 
    //! @param[in ]: transRotationTemperature denotes the non-dimensional value of translation-rotation temperature.
    //! @param[in ]: electronPressure denotes the non-dimensional value of electron pressure.
    //! @param[out]: speciesViscosity is an array to storing the dimensional values of viscosities of species(Unit: kg/(m*s)).
    virtual void ComputeSpeciesDimensionalViscosityWithLennardJonesMethod(RDouble transRotationTemperature, RDouble electronPressure, RDouble *speciesViscosity){};
    
    //! To compute the collision area between the same species using the curve fit method.
    //! @param[in ]: transRotationTemperature denotes the non-dimensional value of translation-rotation temperature.
    //! @param[in ]: electronPressure denotes the non-dimensional value of electron pressure.
    //! @param[out]: piOmega is an array to storing the dimensional values of collision area of species(Unit: m2).
    virtual void ComputeAverageCollisionAreaOmega22(RDouble transRotationTemperature, RDouble electronPressure, RDouble *piOmega){};
    
    //! To obtain the viscosity of the mixed gas using Wilke method.
    //! @param[in ]: moleFraction is an array of saving the mole fractions of species.
    //! @param[in ]: viscosity is an array of saving the viscosity of species.
    //! @param[in ]: phi denotes the array of the partition function.
    //! @param[out]: return the viscosity of the mixed gas.
    virtual RDouble GetMixedGasViscosityWithWilkeFormula(RDouble *moleFraction, RDouble *viscosity, RDouble *phi){return 0;};
    
    //! To obtain the diffusion coefficient of each species. added by LiPeng in Mar. 11,2019.
    //! @param[in ]: massFraction is an array of saving the primary variables.
    //! @param[in ]: viscosityLaminar is the laminar viscosity, and viscosityTurbulence is the turbulent viscosity.
    //! @param[out]: speciesDiffusionCoef stores the species mass diffusion coefficient of each component.
    //! The function returns the diffusion coefficients of each species which saved in rDs = (1 - Cs)/(1- Xs) * (mul/Scl + mut/Sct).
    virtual void ComputeSpeciesMassDiffusionCoefficient(RDouble *massFraction, RDouble viscosityLaminar, RDouble viscosityTurbulence, RDouble *speciesDiffusionCoef){};
    
    //! To obtain the mass diffusion coefficient with CLN method, the obtained value is enjoyed by all of the species.
    //! @param[in ]: primitiveVariables is an array of saving the primitive variables.
    //! @param[in ]: transRotationTemperature denotes the non-dimensional value of translation-rotation temperature.
    //! @param[in ]: electronTemperature denotes the non-dimensional value of electron temperature.
    //! @param[out]: return the mass diffusion coefficient.
    virtual RDouble GetMassDiffusionCoefficientWithCLNModel(RDouble *primitiveVariables, RDouble transRotationTemperature, RDouble electronTemperature){return 0;};
    
    //! Compute the coefficient called gama system in the multi-temperature model.
    //! @param[in ]: primitiveVars is an array of saving the primitive variables.
    //! @param[in ]: electronTemperature is the non-dimensional value of vibration or electron temperature.
    //! @param[out]: gama are partial derivative of electron pressure to vibration or electron energy.
    //! @param[out]: speciesGama is an array of partial derivatives of electron pressure to species density.
    virtual void GetElectronPressurePartialDerivatives(RDouble *primitiveVars, RDouble electronTemperature, RDouble &gama, RDouble *speciesGama){};
    
    //! Compute the coefficient called alpha system in the multi-temperature model.
    //! @param[in ]: primitiveVars is an array of saving the primitive variables.
    //! @param[in ]: transRotationTemperature is the non-dimensional value of translation-rotation temperature.
    //! @param[in ]: electronTemperature is the non-dimensional value of vibration or electron temperature.
    //! @param[out]: alpha is partial derivative of pressure to internal energy.
    //! @param[out]: beta is partial derivative of pressure to total density.
    //! @param[out]: speciesBeta is an array of partial derivatives of pressure to species density..
    virtual void GetPressurePartialDerivatives(RDouble *primitiveVars, RDouble transRotationTemperature, RDouble electronTemperature,
        RDouble &alpha, RDouble &beta, RDouble *speciesBeta){};
    
    //! Compute the main partial derivatives in the multi-temperature model.
    //! @param[in ]: primitiveVars is an array of saving the primitive variables.
    //! @param[in ]: transRotationTemperature is the non-dimensional value of translation-rotation temperature.
    //! @param[in ]: electronTemperature is the non-dimensional value of vibration or electron temperature.
    //! @param[out]: alpha is the partial derivative of pressure to total energy.
    //! @param[out]: beta is the partial derivative of pressure to total density.
    //! @param[out]: gama is the partial derivative of electron pressure to vibration-electron enegy.
    //! @param[out]: speciesBeta¡¢speciesGama are the partial derivatives of pressure and electron pressure to species density.
    virtual void GetMultiTemperatureModelPartialDerivatives(RDouble *primitiveVars, RDouble transRotationTemperature, RDouble electronTemperature,
        RDouble &alpha, RDouble &beta, RDouble &gama, RDouble *speciesBeta, RDouble *speciesGama){};
    
    //! Compute partial derivative of pressure to total density.
    //! @param[in ]: primitiveVars is an array of saving the primitive variables.
    //! @param[out]: return the partial derivative called alpha.
    virtual RDouble ComputePressureDerivativeToTotalDensity(RDouble *primitiveVars){return 0;};
    
    //! To compute the static enthalpy of the mixed gas.
    //! @param[in ]: primitiveVars is an array of saving the primitive variables.
    //! @param[in ]: transRotationTemperature denotes the temperature of translation-rotation mode.
    //! @param[in ]: vibrationTemperature denotes the temperature of vibration mode.
    //! @param[in ]: electronTemperature denotes the temperature of electron mode.
    virtual RDouble GetMixedGasEnthalpy(RDouble *primitiveVars, RDouble transRotationTemperature, RDouble vibrationTemperature, RDouble electronTemperature){return 0;};
    
    virtual RDouble ComputeMolecularWeightReciprocal(RDouble *prim){return 0;};
    
    virtual RDouble ComputeMolecularWeight(RDouble *prim){return 0;};
    
    //! Compute chemical energy of species stored in the variable of qx, qy and qz. The variable of fvis stores\n
    //! the mass diffusion flux namely ro*Ds * (dfsdx * nx + dfsdy * ny + dfsdz * nz).
    //! The function is the optimized format of the previous function, added by LiPeng on April 12, 2019.
    //! @param[in ]: primitiveVars is an array of saving the primitive variables.
    //! @param[in ]: dfsdx, dfsdy, dfsdz are arrays of saving the partial derivatives of the mass fractions.
    //! @param[in ]: transRotationTemperature, vibrationTemperatuare and electronTemperature denote the temperatures of translation-rotation mode, vibration mode and electron mode, respectively.
    //! @param[in ]: viscosityLaminar and viscosityTurbulence are the laminar viscosity and turbulent viscosity, respectively.
    //! @param[in ]: nx, ny and nz are the components of the normal vector on the face.
    //! @param[out]: qx, qy and qz are chemical energy of species in each coordinates.
    //! @param[out]: speciesFluxes is array of saving the mass diffusion flux namely ro*Ds * (dfsdx * nx + dfsdy * ny + dfsdz * nz).
    virtual void ComputeSpeciesMassDiffusionFlux(RDouble *primitiveVars, RDouble *dfsdx, RDouble *dfsdy, RDouble *dfsdz, RDouble transRotationTemperature, RDouble vibrationTemperatuare, RDouble electronTemperature,
        RDouble viscosityLaminar, RDouble viscosityTurbulence, RDouble nx, RDouble ny, RDouble nz, RDouble &qx, RDouble &qy, RDouble &qz, RDouble *speciesFluxes){};
    
    //! To compute heat conductivity of mixture gas using Eucken formula and Wassilewa relation.
    //! @param[in ]: primitiveVariables is an array of saving the primitive variables.
    //! @param[in ]: transRotationTemperature is the non-dimensional value of translation-rotation temperature of the mixture gas.
    //! @param[in ]: vibrationTemperature is the non-dimensional value of vibration temperature of the mixture gas.
    //! @param[in ]: electronTemperature is the non-dimensional value of electron temperature of the mixture gas.
    //! @param[out]: transRotationConductivity is heat conductivity of translation-rotation of mixed gas.
    //! @param[out]: vibrationConductivity is heat conductivity of vibration of mixed gas.
    //! @param[out]: electronConductivity is heat conductivity of electron of mixed gas.
    virtual void ComputeMixtureGasHeatConductivity(RDouble *primitiveVariables, RDouble transRotationTemperature, RDouble vibrationTemperature, RDouble electronTemperature,
        RDouble &transRotationConductivity, RDouble &vibrationConductivity, RDouble &electronConductivity){};
    
    //! Compute the partial derivatives of source terms of chemical reaction equations.
    //! @param[in ]: primitiveVariables stores the primitive variables of the specified control cells.
    //! @param[in ]: temperatures stores the temperature values of the specified control cells.
    //! @param[in ]: cellVolume is an array that stores the volume values of the specified control cells.
    //! @param[in ]: chemicalSourceTerms stores the source terms of chemical reaction equations in the specified control cells.
    //! @param[in ]: chemicalSourceDerivatives stores the partial derivatives of the source terms of chemical reaction equations in the specified control cells.
    //! @param[out]: energySourceDerivatives returns the partial derivatives of the source terms of energy equations in the specified control cells.
    //! @param[in ]: nCellNumber denotes the number of the specified control cells.
    virtual void ComputeEnergySourceDerivatives(RDouble **primitiveVariables, RDouble **temperatures, RDouble *cellVolume, RDouble **chemicalSourceTerms,
        RDouble ***chemicalSourceDerivatives, RDouble ***energySourceDerivatives, int nCellNumber){};
    
    //! Compute the partial derivatives of source terms of chemical reaction equations.
    //! @param[in ]: primitiveVariables stores the primitive variables of the specified control cells.
    //! @param[in ]: temperatures stores the temperature values of the specified control cells.
    //! @param[in ]: cellVolume is an array that stores the volume values of the specified control cells.
    //! @param[in ]: chemicalSourceTerms stores the source terms of chemical reaction equations in the specified control cells.
    //! @param[in ]: chemicalSourceDerivatives stores the partial derivatives of the source terms of chemical reaction equations in the specified control cells.
    //! @param[out]: energySourceDerivatives returns the diagonal element of partial derivatives to the source terms of energy equations in the specified control cells.
    //! @param[in ]: nCellNumber denotes the number of the specified control cells.
    virtual void ComputeEnergySourceDiagonalDerivatives(RDouble **primitiveVariables, RDouble **temperatures, RDouble *cellVolume, RDouble **chemicalSourceTerms,
        RDouble ***chemicalSourceDerivatives, RDouble **energySourceDerivatives, int nCellNumber){};

    virtual void ComputeEnergySourceJacobianT2(RDouble *primitiveVariables, RDouble *temperatures, RDouble cellVolume, RDouble *chemicalSourceTerms, RDouble **sourceDerivatives){};
    virtual void ComputeEnergySourceJacobianT3(RDouble *primitiveVariables, RDouble *temperatures, RDouble cellVolume, RDouble *chemicalSourceTerms, RDouble **sourceDerivatives){};

    virtual void ComputeEnergySourceAndJacobianT2(RDouble **primitiveVariables, RDouble **temperatures, RDouble *cellVolume, RDouble **chemicalSourceTerms, RDouble ***sourceDerivatives, int nCellNumber){};
    virtual void ComputeEnergySourceAndJacobianT3(RDouble **primitiveVariables, RDouble **temperatures, RDouble *cellVolume, RDouble **chemicalSourceTerms, RDouble ***sourceDerivatives, int nCellNumber){};

    //! To obtain the maximum mass diffusion coefficient among the species. added by LiPeng in Jan. 2019.
    //! @param[in ]: primitiveVars is an array of saving the primitive variables.
    //! @param[in ]: viscosityLaminar is the laminar viscosity, and viscosityTurbulence is the turbulent viscosity.
    //! The function returns the maximum mass diffusion coefficient called (1 - Cs)/(1- Xs) * (mul/Scl + mut/Sct).
    virtual RDouble GetMaximumSpeciesMassDiffusionCoef(RDouble *primitiveVars, RDouble viscosityLaminar, RDouble viscosityTurbulence){return 0;};
    
    virtual RDouble ComputeMolecularWeightReciprocalWithoutElectron(RDouble *prim, RDouble &ce_div_me){return 0;};
    
    //! To obtain the diffusion coefficient of each species. added by LiPeng in Mar. 11,2019.
    //! @param[in ]: primitiveVars is an array of saving the primitive variables.
    //! @param[in ]: viscosityLaminar is the laminar viscosity, and viscosityTurbulence is the turbulent viscosity.
    //! @param[out]: speciesDiffusionCoef stores the species mass diffusion coefficient of each component.
    //! The function returns the diffusion coefficients of each species which saved in rDs = (1 - Cs)/(1- Xs) * (mul/Scl + mut/Sct).
    virtual void GetSpeciesMassDiffusionCoef(RDouble *primitiveVars, RDouble viscosityLaminar, RDouble viscosityTurbulence, RDouble *speciesDiffusionCoef){};
    
    //! To obtain the non-dimensional constant pressure specific heat in accordance with the primary variables.
    //! prim is an array of saving the primary variables, the non-dimensional temperature called Tm and the non-dimensional \n
    //! constant pressure specific heat called cp are returned.
    virtual void GetTemperatureAndCP(RDouble *prim, RDouble &Tm, RDouble &cp){};
    
    //! To obtain the non-dimensional enthalpy of each species, the variable Tm is the non-dimensional temperature of the mixture gas,\n
    //! the enthalpy of each species saved in the array of hs is returned.
    virtual void GetEverySpeciesEnthalpy(RDouble trTemperature, RDouble vTemperature, RDouble eTemperature, RDouble *hs){};
    
    //! To compute the electron energy of the mixed gas.
    //! @param[in ]: prim is an array of saving the primary variables.
    //! @param[in ]: Te is the non-dimensional electron temperature of the mixture gas.
    //! @param[out]: return the total electron energy.
    virtual RDouble GetMixedGasElectronEnergy(RDouble *prim, RDouble Te){return 0;};
    
    //! To compute the vibration specific heat at constant volume of the mixed gas.
    //! @param[in ]: prim is an array of saving the primary variables.
    //! @param[in ]: Tv is the non-dimensional vibration temperature of the mixture gas.
    //! @param[out]: return the total vibration specific heat at constant volume.
    virtual RDouble GetMixedGasVibrationSpecificHeatAtConstantVolume(RDouble *prim, RDouble Tv){return 0;};
    
    //! To compute the vibration specific heat at constant volume of the mixed gas.
    //! @param[in ]: prim is an array of saving the primary variables.
    //! @param[in ]: Te is the non-dimensional electron temperature of the mixture gas.
    //! @param[out]: return the dimensional value of total vibration specific heat at constant volume.
    virtual RDouble GetMixedGasElectronSpecificHeatAtConstantVolume(RDouble *prim, RDouble Te){return 0;};
    
    //! To compute the translation-rotation specific heat at constant volume of the mixed gas.
    //! @param[in ]: prim is an array of saving the primary variables.
    //! @param[out]: return the dimensional value of total translation-rotation specific heat at constant volume.
    virtual RDouble GetMixedGasTransAndRotatSpecHeatAtConstVolume(RDouble *prim){return 0;};
    
    //! To compute the translation-rotation energy of the mixed gas.
    //! @param[in ]: prim is an array of saving the primary variables.
    //! @param[in ]: Ttr is the non-dimensional translation-rotation temperature of the mixture gas.
    //! @param[out]: return the total translation-rotation energy.
    virtual RDouble GetMixedGasTranslationAndRotationEnergy(RDouble *prim, RDouble Ttr){return 0;};
    
    //! To compute the vibration energy of the mixed gas.
    //! @param[in ]: prim is an array of saving the primary variables.
    //! @param[in ]: Tv is the non-dimensional vibration temperature of the mixture gas.
    //! @param[out]: return the total vibration energy.
    virtual RDouble GetMixedGasVibrationEnergy(RDouble *prim, RDouble Tv){return 0;};
    
    //! Compute the frozen sound speed in the multi-temperature model.
    //! @param[in ]: primitiveVars is an array of saving the primitive variables.
    //! @param[out]: return the frozen sound speed.
    virtual RDouble ComputeFrozenSoundSpeed(RDouble *primitiveVars){return 0;};
    
    //! To obtain the non-dimension temperature according to the primary variables,\n
    //! prim is an array of saving the primary variables, the non-dimensional temperature is returned.
    virtual RDouble GetTemperature(RDouble *prim){return 0;};

};

extern Gas *gas;

}
}
