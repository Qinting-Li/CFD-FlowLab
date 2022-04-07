//==============================================
//! @file      ChemicalModel.h
//! @brief     define chemical reactions.
//! @author    Li Peng.
//==============================================
#pragma once
#include <string>
#include "Precision.h"
using namespace std;

namespace PHSPACE
{
    namespace GAS_SPACE
    {
        //! Define the class of chemical model that is used to construct the chemical reactions and their reaction rates.
        class ChemicalModel{
        private:
            int nReactionType; //nType = 0 indicates the Gupta Model, nType = 1 indicates the Park Model, and nType = 2 indicates the Dunn-Kang Model.
            //nType = 3 indicates the Park Model of Mars, and nType = 4 indicates the McKenzie Model of Mars.
        private:
            int nSpecies, nReaction;
            int **stoichForward, **stoichBackward;
            RDouble **collisionCoef;
            RDouble **reactionData; //curve fitting data of reaction rate.
            int **reactionType;
        public:
            ChemicalModel();
            ChemicalModel(int ns, int nr);
            ~ChemicalModel();
            void AllocateMemories();
            void DeAllocateMemories();
            int GetChemicalReactionType() const {return nReactionType;};
            void SetChemicalReactionType(int ntype) {this->nReactionType = ntype;};
            void GetChemicalModelInfo(int &ns, int &nr) {ns = nSpecies; nr = nReaction;};
            void SetChemicalModelInfo(int ns, int nr);
            int** GetReactionType() const {return reactionType;}
            int** GetForwardStoichiometricCoefficient() const {return stoichForward;};
            int** GetBackwardStoichiometricCoefficient() const {return stoichBackward;};
            RDouble** GetThirdBodyCollisionCoefficient() const {return collisionCoef;};
            RDouble** GetReactionRateCurveFittingData() const {return reactionData;};
        public:
            //! Create the Gupta Reaction Model.
            //! Param[in]: ns indicates the number of species, 11 is set as the default value.
            //! Param[in]: nr indicates the number of reactions, 20 is set as the default value.
            //! the order of species is O, O2, NO, N, O+, O2+, NO+, N+, N2+, N2, e-.
            void CreateGuptaModel(int ns = 11, int nr = 20);

            //! Create the Park Reaction Model.
            //! Param[in]: ns indicates the number of species, 11 is set as the default value.
            //! Param[in]: nr indicates the number of reactions, 20 is set as the default value.
            //! the order of species is O, O2, NO, N, O+, O2+, NO+, N+, N2+, N2, e-.
            void CreateParkModel(int ns = 11, int nr = 20);

            //! Create the Dunn-Kang Reaction Model.
            //! Param[in]: ns indicates the number of species, 11 is set as the default value.
            //! Param[in]: nr indicates the number of reactions, 26 is set as the default value.
            //! the order of species is O, O2, NO, N, O+, O2+, NO+, N+, N2+, N2, e-.
            void CreateDunnKangModel(int ns = 11, int nr = 26);

            //! Create the Mars Reaction Model.
            //! Param[in]: ns indicates the number of species, 8 is set as the default value.
            //! Param[in]: nr indicates the number of reactions, 12 is set as the default value.
            //! Param[in]: nType = 0 indicates the Park model, nType = 1 indicates the McKenzie model.
            //! the order of species is O, O2, NO, N, N2, C, CO, CO2.
            void CreateMarsGasModel(int ns = 8, int nr = 12, int nType = 0);

            //! To obtain the data of the Gupta 5s6r reaction model.
            //! Param[in/out]: reactionCoef denotes the stoichiometric coefficients of forward and backward reaction equation.
            //! Param[in/out]: collisionCoef indicates the collision coefficients of third body in the reaction equation.
            //! Param[in/out]: rateCoef indicates the curve fitting data of reaction rates.
            void GenerateGupta5S6RData(int *reactionCoef, RDouble *collisionCoef, RDouble *rateCoef);

            //! To obtain the data of the Gupta 7s9r reaction model.
            //! Param[in/out]: reactionCoef denotes the stoichiometric coefficients of forward and backward reaction equation.
            //! Param[in/out]: collisionCoef indicates the collision coefficients of third body in the reaction equation.
            //! Param[in/out]: rateCoef indicates the curve fitting data of reaction rates.
            void GenerateGupta7S9RData(int *reactionCoef, RDouble *collisionCoef, RDouble *rateCoef);

            //! To obtain the data of the Gupta 11s20r reaction model.
            //! Param[in/out]: reactionCoef denotes the stoichiometric coefficients of forward and backward reaction equation.
            //! Param[in/out]: collisionCoef indicates the collision coefficients of third body in the reaction equation.
            //! Param[in/out]: rateCoef indicates the curve fitting data of reaction rates.
            void GenerateGupta11S20RData(int *reactionCoef, RDouble *collisionCoef, RDouble *rateCoef);

            //! To obtain the data of the Park 5s17r reaction model.
            //! Param[in/out]: reactionCoef denotes the stoichiometric coefficients of forward and backward reaction equation.
            //! Param[in/out]: collisionCoef indicates the collision coefficients of third body in the reaction equation.
            //! Param[in/out]: rateCoef indicates the curve fitting data of reaction rates.
            void GeneratePark5S17RData(int *reactionCoef, RDouble *collisionCoef, RDouble *rateCoef);

            //! To obtain the data of the Park 7s22r reaction model.
            //! Param[in/out]: reactionCoef denotes the stoichiometric coefficients of forward and backward reaction equation.
            //! Param[in/out]: collisionCoef indicates the collision coefficients of third body in the reaction equation.
            //! Param[in/out]: rateCoef indicates the curve fitting data of reaction rates.
            void GeneratePark7S22RData(int *reactionCoef, RDouble *collisionCoef, RDouble *rateCoef);

            //! To obtain the data of the Park 11s48r reaction model.
            //! Param[in/out]: reactionCoef denotes the stoichiometric coefficients of forward and backward reaction equation.
            //! Param[in/out]: collisionCoef indicates the collision coefficients of third body in the reaction equation.
            //! Param[in/out]: rateCoef indicates the curve fitting data of reaction rates.
            void GeneratePark11S48RData(int *reactionCoef, RDouble *collisionCoef, RDouble *rateCoef);

            //! To obtain the data of the Dunn-Kang 5s11r reaction model.
            //! Param[in/out]: reactionCoef denotes the stoichiometric coefficients of forward and backward reaction equation.
            //! Param[in/out]: collisionCoef indicates the collision coefficients of third body in the reaction equation.
            //! Param[in/out]: rateCoef indicates the curve fitting data of reaction rates.
            void GenerateDunnKang5S11RData(int *reactionCoef, RDouble *collisionCoef, RDouble *rateCoef);

            //! To obtain the data of the Dunn-Kang 7s15r reaction model.
            //! Param[in/out]: reactionCoef denotes the stoichiometric coefficients of forward and backward reaction equation.
            //! Param[in/out]: collisionCoef indicates the collision coefficients of third body in the reaction equation.
            //! Param[in/out]: rateCoef indicates the curve fitting data of reaction rates.
            void GenerateDunnKang7S15RData(int *reactionCoef, RDouble *collisionCoef, RDouble *rateCoef);

            //! To obtain the data of the Dunn-Kang 11s26r reaction model.
            //! Param[in/out]: reactionCoef denotes the stoichiometric coefficients of forward and backward reaction equation.
            //! Param[in/out]: collisionCoef indicates the collision coefficients of third body in the reaction equation.
            //! Param[in/out]: rateCoef indicates the curve fitting data of reaction rates.
            void GenerateDunnKang11S26RData(int *reactionCoef, RDouble *collisionCoef, RDouble *rateCoef);

            //! To obtain the data of the Park 12s12r reaction model of Mars.
            //! Param[in/out]: reactionCoef denotes the stoichiometric coefficients of forward and backward reaction equation.
            //! Param[in/out]: collisionCoef indicates the collision coefficients of third body in the reaction equation.
            //! Param[in/out]: rateCoef indicates the curve fitting data of reaction rates.
            void GenerateMarsPark8S12RData(int *reactionCoef, RDouble *collisionCoef, RDouble *rateCoef);

            //! To obtain the data of the McKenzie 12s12r reaction model of Mars.
            //! Param[in/out]: reactionCoef denotes the stoichiometric coefficients of forward and backward reaction equation.
            //! Param[in/out]: collisionCoef indicates the collision coefficients of third body in the reaction equation.
            //! Param[in/out]: rateCoef indicates the curve fitting data of reaction rates.
            void GenerateMarsMcKenzie8S12RData(int *reactionCoef, RDouble *collisionCoef, RDouble *rateCoef);

        };

        //! Define the class of ChemkinData that is used to compute the specific heat of species via curve fitting method.
        class ChemkinData{
        private:
            int nSpecies;//number of species;
            RDouble ***curveFittingData;
            int nInterval, nCoefNum;
            RDouble *temperatureRange;
        public:
            ChemkinData();
            ~ChemkinData();
            void AllocateMemories();
            void DeAllocateMemories();
            void GetInterpolationInfo(int &ns, int &m, int &n) {ns = nSpecies; m = nInterval; n = nCoefNum;}
            void SetInterpolationInfo(int ns, int m = 5, int n = 6);
            RDouble*** GetCurveFittingData() const {return curveFittingData;}
            RDouble* GetTemperatureRange() const {return temperatureRange;}
        public:
            //! To generate the data of chemkin method.
            //! Param[in]: ns denotes the number of species, the default value is set 11.
            //! the order of species is O, O2, NO, N, O+, O2+, NO+, N+, N2+, N2, e-, C, C2, CN, CO, CO2.
            void CreateChemkinData(int nGasModel = 0, int ns = 11);

            //! To Generate the curve fitting data of 11 species in the chemkin method.
            //! Param[out]: fittingData indicates the three-dimensional array of saving the curve fitting data.
            void GenerateFullChemkinData(RDouble ***fittingData);
        };

        //! Define the class of SpeciesData that is used to save the physical and chemical data.
        //! This is the atmosphere of Earth.
        class SpeciesData{
        private:
            int nSpecies; //number of species;
            string *speciesName; // the name of species.
            int *nSpeciesType; // the type of species.
            int *nSpeciesCharge; // the charge of species.
            int *nSpeciesG0; // the degeneracy of species.
            int *nSpeciesG1; // the degeneracy of species.
            RDouble *speciesMass; // the molecular weight of species, (g/mol).
            RDouble *speciesTtr; // the characteristic temperature of species, (K).
            RDouble *speciesTvs; // the characteristic vibration temperature of species, (K).
            RDouble *speciesTes; // the characteristic electron temperature of species, (K).
            RDouble *speciesEnthalpy; // the formation enthalpy of species, (J/kg).
            RDouble *speciesSigma; // the collision area of species, (m2).
        public:
            SpeciesData();
            ~SpeciesData();
            void AllocateMemories();
            void DeAllocateMemories();
            void SetSpeciesNumber(int ns);
            string* GetSpeciesName() const {return speciesName;}
            int GetSpeciesNumber() const {return nSpecies;}
            int* GetSpeciesType() const {return nSpeciesType;}
            int* GetSpeciesCharge() const {return nSpeciesCharge;}
            RDouble* GetSpeciesMolecularWeight() const{ return speciesMass; }
            int* GetSpeciesDegeneracy0() const {return nSpeciesG0;}
            int* GetSpeciesDegeneracy1() const {return nSpeciesG1;}
            RDouble* GetSpeciesTemperature() const {return speciesTtr;}
            RDouble* GetSpeciesVibrationTemperature() const {return speciesTvs;}
            RDouble* GetSpeciesElectronTemperature() const {return speciesTes;}
            RDouble* GetSpeciesFormationEnthalpy() const {return speciesEnthalpy;}
            RDouble* GetSpeciesCollisionArea() const {return speciesSigma;}
        public:
            //! To generate the physical and chemical data of species.
            //! Param[in]: ns denotes the number of species, the default value is set 11.
            //! the order of species is O, O2, NO, N, O+, O2+, NO+, N+, N2+, N2, e-.
            void CreatePhysicalChemicalData(int ns = 11);

            //! To Generate the physical and chemicla data of 11 species.
            //! Param[out]: nameArray indicates the one-dimensional array of saving the species name.
            //! Param[out]: dataArray1 indicates the two-dimensional array of saving the type and charge.
            //! Param[out]: dataArray2 indicates the two-dimensional array of saving the other data.
            void GenerateFullPhysicalChemicalData(string *nameArray, int **dataArray1, RDouble **dataArray2);

            //! To obtain the index of species in the sequence..
            //! Param[in]: speciesName indicates the name of species.
            static int GetSpeciesIndex(string speciesName);

            //! To obtain the name of species.
            //! Param[in]: nIndex indicates the index of species in the sequence.
            static string GetSpeciesName(int nIndex);
        };

        //! Define the class of TransportData that is used to compute the transport coefficients such as viscosity.
        class TransportData{
        private:
            int nSpecies; //number of species;
            RDouble **viscosityFitData; //the curve fitting data of viscosity.
            RDouble **collisionAreaFitData; // the curve fitting data of the average collision area.
        public:
            TransportData();
            ~TransportData();
            void AllocateMemories();
            void DeAllocateMemories();
            void SetSpeciesNumber(int ns);
            int GetSpeciesNumber() const {return nSpecies;}
            RDouble** GetViscosityCurveFittingData() const {return viscosityFitData;}
            RDouble** GetCollisionAreaCurveFittingData() const {return collisionAreaFitData;}
        public:
            //! To generate the curve fitting data of viscosity and average collision area.
            //! Param[in]: nGasModel = 0 indicates the Earth atmosphere, while nGasModel = 1 indicates the Mars atmopsphere.
            //! Param[in]: ns denotes the number of species, the default value is set 11.
            //! the order of species for earth atmopsphere is O, O2, NO, N, O+, O2+, NO+, N+, N2+, N2, e-.
            //! the order of species for mars atmopsphere is O, O2, NO, N, N2, C, CO, CO2, C2, CN, NO+, e-.
            void CreateCurveFittingData(int nGasModel = 0, int ns = 11);

            //! To Generate the curve fitting data of 11 species.
            //! Param[out]: viscosityData indicates the two-dimensional array of saving the curve fitting data of viscosity.
            //! Param[out]: collisionData indicates the two-dimensional array of saving the curve fitting data of average collision area.
            //! Param[in]: nGasModel = 0 indicates the Earth atmosphere, while nGasModel = 1 indicates the Mars atmopsphere.
            void GenerateFullCurveFittingData(RDouble **viscosityData, RDouble **collisionData, int nGasModel = 0);
        };

        //! Define the class of MarsSpeciesData that is used to save the physical and chemical data.
        //! This is the atmosphere of Earth.
        class MarsSpeciesData{
        private:
            int nSpecies; //number of species;
            string *speciesName; // the name of species.
            int *nSpeciesType; // the type of species.
            int *nSpeciesCharge; // the charge of species.
            int *nSpeciesG0; // the degeneracy of species.
            int *nSpeciesG1; // the degeneracy of species.
            RDouble *speciesMass; // the molecular weight of species, (g/mol).
            RDouble *speciesTtr; // the characteristic temperature of species, (K).
            RDouble *speciesTvs; // the characteristic vibration temperature of species, (K).
            RDouble *speciesTes; // the characteristic electron temperature of species, (K).
            RDouble *speciesEnthalpy; // the formation enthalpy of species, (J/kg).
            RDouble *speciesSigma; // the collision area of species, (m2).
        public:
            MarsSpeciesData();
            ~MarsSpeciesData();
            void AllocateMemories();
            void DeAllocateMemories();
            void SetSpeciesNumber(int ns);
            string* GetSpeciesName() const {return speciesName;}
            int GetSpeciesNumber() const {return nSpecies;}
            int* GetSpeciesType() const {return nSpeciesType;}
            int* GetSpeciesCharge() const {return nSpeciesCharge;}
            RDouble* GetSpeciesMolecularWeight() const{ return speciesMass; }
            int* GetSpeciesDegeneracy0() const {return nSpeciesG0;}
            int* GetSpeciesDegeneracy1() const {return nSpeciesG1;}
            RDouble* GetSpeciesTemperature() const {return speciesTtr;}
            RDouble* GetSpeciesVibrationTemperature() const {return speciesTvs;}
            RDouble* GetSpeciesElectronTemperature() const {return speciesTes;}
            RDouble* GetSpeciesFormationEnthalpy() const {return speciesEnthalpy;}
            RDouble* GetSpeciesCollisionArea() const {return speciesSigma;}
        public:
            //! To generate the physical and chemical data of species.
            //! Param[in]: ns denotes the number of species, the default value is set 12.
            //! the order of species is O, O2, C, CO, CO2, NO, N, N2, C2, CN, NO+, e-.
            void CreatePhysicalChemicalData(int ns = 12);

            //! To Generate the physical and chemicla data of 11 species.
            //! Param[out]: nameArray indicates the one-dimensional array of saving the species name.
            //! Param[out]: dataArray1 indicates the two-dimensional array of saving the type and charge.
            //! Param[out]: dataArray2 indicates the two-dimensional array of saving the other data.
            void GenerateFullPhysicalChemicalData(string *nameArray, int **dataArray1, RDouble **dataArray2);

            //! To obtain the index of species in the sequence..
            //! Param[in]: speciesName indicates the name of species.
            static int GetSpeciesIndex(string speciesName);

            //! To obtain the name of species.
            //! Param[in]: nIndex indicates the index of species in the sequence.
            static string GetSpeciesName(int nIndex);
        };

    }// end of GAS_SPACE

}// end of PHSPACE
