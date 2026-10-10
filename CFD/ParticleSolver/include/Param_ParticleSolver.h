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
//! @file      Param_ParticleSolver.h
//! @brief     The param of Particle solver.
//! @author    Lei Yinghaonan (Lanzhou University).

#pragma once
#include "Param_CFDSolver.h"
#include "Data_Param.h"

namespace PHSPACE
{

class Param_ParticleSolver : public Param_CFDSolver
{
private:
    //! Particle solver model: The type of particle solver.
    //!                        1 -- Particle Point Solver.
    //!                        2 -- Particle Resolved(IBM) Solver.
    int iParticleModel;

    //! Init or restart particle file.
    string initParticleFile;

    //! ifReadRestartParticle : If read restart particle field.
    //!                         0 -- Init particle field by "initParticleFile".
    //!                         1 -- Read particle field by "initParticleFile".
    int ifReadRestartParticle;
    //! ifWirteRestartParticle : If write restart particle file.
    //!                          0 -- Do not write restart particle file : par_init.h5
    //!                          1 -- Write restart particle file. file : par_xxx.h5
    int ifWirteRestartParticle;

    //! The step intervals for particle file.
    int intervalStepRestartParticle;

    //! If init particle info as init file.
    //! ifInitParticleBC : If use particle BC as REMOVE_PARTICLE_AS_INIT.
    //!                    0 -- not use.
    //!                    1 -- when remove particle bc as REMOVE_PARTICLE_AS_INIT.
    //!                         reset particle as "initParticleBCFile".
    //!                    2 -- reset particle only coordinate and cellID as "initParticleBCFile".
    //!                         but other variable but others don't change.
    //!                    3 -- reset particle as "initParticleBCFile" with with periodic time and
    //!                         random with velocity.
    int ifInitParticleBC;

    //! initParticleBCPeriodicStep : The step interval with periodic time.
    //!                              The Only use when ifInitParticleBC is 3.
    int initParticleBCPeriodicStep;
    //! initParticleBCRandomDirectionalRatio: The proportion that changes the direction 
    //!                                       of the particle velocity during initial release. 
    //!                                       Notice that this is proportional to [-PI,PI]. 
    //!                                       For example, the ratio 0.5 corresponds to [-0.5*PI,0.5*PI].
    //!                                       The Only use when ifInitParticleBC is 3.
    int particleRandomReleaseMode;
    int particleRandomDiameterMode;

    RDouble randomDiameterMaxRatio;
    RDouble randomDiameterMinRatio;
    RDouble initParticleBCRandomDirectionalRatio;
    RDouble initParticleBCRandomValueRatio;
    
    string initParticleBCFile;

    //! gridTypeForPreAndPost : The grid type for pre - and post-processing
    //!                         0 : Default arbitrary mesh (structured and unstructured)
    //!                         1 : Orthogonal grid (only structured).
    int gridTypeForPreAndPost = 0;

    //! The version of particle file.
    RDouble version;

    //! Num of particle sum in all zone.
    int nParticleTotal;
    //! Num of particle in current zone.
    int nParticleLocal;

    //! Currently, there is only one Group on each zone. 
    //! But in the future, 
    //! it might be possible to have multiple groups on a zone, 
    //! such as the CMP-PIC method.
    int nParticleGroup;

    //! Here we explain what the grid Factor does.
    //! For example, CGNS grids mesh.cgns do not contain units in their output.
    //! So when we do the grid transformation, 
    //! the coordinate X (in cgns) is multiplied by the grid factor to become a length unit in meters.
    //! For example, in CGNS the length is 1.2, which is dimensionless.
    //! When we use grid factor 0.01, it means that the length unit is mm, 
    //! so the length variable in the FTS file after conversion becomes 0.012 (m).
    //! 
    //reynoldsReferenceLengthDimensional and gridScaleFactor

    //! DimensionlnalParam -- flow ref var.
    RDouble refReynoldsLengthDimensional;
    RDouble refVelocityDimensional;
    RDouble refTemperatureDimensional;
    RDouble refPressureDimensional;
    RDouble refDensityDimensional;
    RDouble refAverageGeneralGasConstantDimensional;
    RDouble refDynamicViscosityDimensional;
    RDouble refWallTemperatureDimensional;
    //! DimensionlnalParam -- flow const.
    RDouble refMachNumber;
    RDouble refReNumber;
    RDouble refGama;
    
    //! DimensionlnalParam -- 
    RDouble refMassDimensional;
    RDouble refTimeDimensional;
    RDouble refAngularVelocityDimensional;
    RDouble refAccelerationyDimensional;
    RDouble refForceDimensional;

    //! particleHeatCD : the specific heat of the disperse phase
    RDouble particleHeatCD;
    //! particle radiative temperature
    RDouble radiativeTemperature;

    //! the specific heat of the carrier phase
    //! Specific heat at constant pressure
    RDouble flowHeatCPConstPressureDim;
    RDouble flowHeatCPConstPressureDimless;

    //! Prandtl number of flow.
    RDouble prl;

    int enableRadiativeHeatTransfer;
    int enableVariableDiameter;

    //! ifInitParticleAsFlow : If init particle velocity, temperature as local flow.
    //!                        0 -- Do note set.
    //!                        1 -- Set as local flow.
    int ifInitParticleAsFlow;

    int initOutStep;

    //! The param from NSSolver.
    int nNSEquation;
    int nLaminar;
    int nChemical;
    int nTemperatureModel;

public:
    LIB_EXPORT Param_ParticleSolver();
    LIB_EXPORT ~Param_ParticleSolver();

public:
    LIB_EXPORT void Init();

    int GetiParticleModel();

    string GetInitParticleFile();
    string GetInitParticleBCFile();

    RDouble GetVersion();
    
    bool IfInitParticle();
    bool IfReadRestartParticle();
    bool IfWirteRestartParticle();
    bool ifInitParticleBCInfo();

    int GetIntervalStepRestartParticle();

    int GetNSEquationNumber();
    int GetLaminarNumber();
    int GetChemicalFlag();
    int GetTemperatureModel();

    RDouble GetRefDimLength();
    RDouble GetRefDimVelocity();
    RDouble GetRefDimTemperature();
    RDouble GetRefDimPressure();
    RDouble GetRefDimDensity();
    RDouble GetRefDimR();

    RDouble GetRefDimMass();
    RDouble GetRefDimTime();
    RDouble GetRefDimAngularVelocity();
    RDouble GetRefDimAcceleration();

    RDouble GetRefMach();
    RDouble GetRefReIn();
    RDouble GetRefGamma();
    RDouble GetRefDimDynamicViscosity();
    RDouble GetRefDimForce();

    void SetNumParticleTotal(int nParticleTotal);
    void SetNumParticleLocal(int nParticleLocal);
    void SetNumParticleGroup(int &nGroup);

    int GetNumParticleTotal();
    int GetNumParticleLocal();
    int GetNumParticleGroup();

    int GetInitParticleBCInfoType();

    void GetRefDataParam(Data_Param *parameter);

    bool ifInitParticleInFlow();

    void SetInitOutStep(int initOutStep);

    int GetInitOutStep();

    bool ifIsOrthogonalForPreAndPost();

    bool ifAddParticleInitBC();

    int GetInitParticleBCPeriodicStep();

    int GetParticleRandomReleaseMode();

    int GetParticleRandomParticleDiameterMode();

    RDouble GetInitParticleBCRandomDirectionalRatio();
    RDouble GetInitParticleBCRandomValueRatio();
    RDouble GetRandomDiameterMaxRatio();
    RDouble GetRandomDiameterMinRatio();
private:
    //! DimensionlessParam
    void InitDimensionlessParam();
};

}