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
//! @file      Param_ParticlePointSolverParallelDomain.h
//! @brief     The param of Particle point solver.
//! @author    Lei Yinghaonan (Lanzhou University).

#pragma once
#include "Param_ParticleSolver.h"
#include "Data_Param.h"

using namespace std;

namespace PHSPACE
{
class Param_ParticlePointSolverParallelDomain : public Param_ParticleSolver
{
private:

    //! simulationSplashFunction: The type of SplashFunction.
    //!                           0 -- No Splash.
    //!                           1 -- Use splash function on wall.
    //! 
    int simulationSplashFunction;

    //! simulationPointRotate: If use Rotate.
    //!                        0 -- No Rotate.
    //!                        1 -- Caculate the Rotate of point.
    int simulationPointRotate;

    //! simulationBackCouping: If use back-couping.
    //!                        0 -- One way couping.
    //!                        1 -- Caculate back force on flow..
    int simulationParticleTemperature;

    //! simulationParticleAcceleration : If use particle acceleration.
    //!                                  0 -- No acceleration.
    //!                                  1 -- Caculate particle acceleration.
    int simulationParticleAcceleration;

    //! simulationBackCouping: If use back-couping.
    //!                        0 -- One way couping.
    //!                        1 -- Caculate back force on flow.
    int simulationBackCouping;

    //! particleInterpolation_struct: The method of particle interpolation on struct grid.
    //!                               TrilinearInterpolation -- TrilinearInterpolation.
    string particleInterpolation_struct;
    string particleInterpolation_unstruct;

    //! nInitGlobalParticle: The number of particles on all zones when program start.
    int nInitGlobalParticle;
    int nParticleGroup;

    //! Gravity type.
    int forceGravityType;

    //! The quasi-steady force (steady-state resistance) force
    int forceDragType;

    //! The lift force
    int forceSaffmanType;
    int forceMagnusType;
    
    //! The added-mass force
    int forceAddedMassType;

    //! The stress-gradient force
    int forceFluidAccelerationType;

    //! Brownain force
    int forceBrownianType;

    //! Thermophoretic force.
    int forceThermophoreticType;

    //! The Basset history force;
    int forceBassetType;

    //! Statistic for particle.
    //! Particle statistic file.
    string fileParticleOutput;

    //! Particle statistic file.
    int fileParticleOutputType;

    //! particle temperature model.
    int particleTemperatureModel;

    //! Runge Kutte method.
    int nRK;

    //! Number of equation for particle.
    int nEquationParticle;

    //! Number of temperatura model for partocle.
    int nDimTem;

    double domainMesh[6];

public:
    LIB_EXPORT Param_ParticlePointSolverParallelDomain();
    LIB_EXPORT ~Param_ParticlePointSolverParallelDomain();

public:
    LIB_EXPORT void Init();
    void ReadParameter();

    void GetNumOfParticleGroup();
    int GetNumOfInitGlobalParticle();

    int GetForceGravityType();
    int GetForceDragType();
    int GetForceMagnusType();
    int GetForceSaffmanType();
    int GetForceBrownianType();
    int GetForceThermophoreticType();
    int GetForceAddedMassType();
    int GetForceBassetType();
    int GetFluidAccelerationType();

    int GetNumRK();
    int GetNumDimTemperature();
    RDouble *GetDomainMesh();

    string GetParticleInterpolationStruct();
    string GetParticleInterpolationUnstruct();

    int GetNumOfParticleEquation();

    void GetForceTypeDataParam(Data_Param *parameter);


    bool ifUseParticleSplash();
    bool ifUseParticleRotate();
    bool ifUseParticleTemperature();
    bool ifUseParticleBackCouping();
    bool ifUseParticleAcceleration();

    void SetNumDimTemperature(int nDimTemp);

private:

};

}