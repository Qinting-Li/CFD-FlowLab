//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//          PPPPP  H   H  EEEEE  N    N  GGGGG  L      EEEEE  III         +
//          P   P  H   H  E      NN   N  G      L      E       I          +
//          PPPPP  HHHHH  EEEEE  N N  N  G  GG  L      EEEEE   I          +
//          P      H   H  E      N  N N  G   G  L      E       I          +
//          P      H   H  EEEEE  N    N  GGGGG  LLLLL  EEEEE  III         +
//------------------------------------------------------------------------+
//          Platform for Hybrid Engineering Simulation of Flows           +
//          China Aerodynamics Research and Development Center             +
//                     (C) Copyright, Since 2010                          +
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      OnePointVariable.h
//! @brief     It is the variabl of one point particle.
//! @author    Lei Yinghaonan (Lanzhou University).

#pragma once
#include "SimplePointer.h"
#include "Data_Param.h"

namespace PHSPACE
{

namespace INDEX_FLOWONPARTICLE
{

//! Index of flowOnParticleInEular
const int nVarOnFlow = 19;

const int PFU = 0;
const int PFV = 1;
const int PFW = 2;

//! flow density.
const int PFR = 3;

//! flow pressure.
const int PFP = 4;

//! flow dynamic viscosity
const int PFMU = 5;

//! flow temperature.
const int PFT = 6;
const int PFDTDX = 16;
const int PFDTDY = 17;
const int PFDTDZ = 18;

//! flow vorticity.
const int PFDUDX = 7;
const int PFDUDY = 8;
const int PFDUDZ = 9;

const int PFDVDX = 10;
const int PFDVDY = 11;
const int PFDVDZ = 12;

const int PFDWDX = 13;
const int PFDWDY = 14;
const int PFDWDZ = 15;

}

namespace INDEX_PARTICLE_FORCE
{
const int nForceDir = 24;

//! 0 -- Gravity force.
const int FG = 0;
const int FGX = 0;
const int FGY = 1;
const int FGZ = 2;

//! 1 -- Drag force.
const int FD = 3;
const int FDX = 3;
const int FDY = 4;
const int FDZ = 5;

//! 2 -- Saffman force.
const int FSAFF = 6;
const int FSAFFX = 6;
const int FSAFFY = 7;
const int FSAFFZ = 8;

//! 3 -- Magnus force.
const int FMAG = 9;
const int FMAGX = 9;
const int FMAGY = 10;
const int FMAGZ = 11;

//! 4 -- Addmass force.
const int FAM = 12;
const int FAMX = 12;
const int FAMY = 13;
const int FAMZ = 14;

//! 5 -- FluidAcceleration force.
const int FFA = 15;
const int FFAX = 15;
const int FFAY = 16;
const int FFAZ = 17;

//! 6 -- Brownian force.
const int FB = 18;
const int FBX = 18;
const int FBY = 19;
const int FBZ = 20;

//! 7 -- Thermophoretic force.
const int FT = 21;
const int FTX = 21;
const int FTY = 22;
const int FTZ = 23;
}

namespace INDEX_RK
{

//! dx/dt = v
const int RK_COORD = 0;
const int RK_COORD_X = 0;
const int RK_COORD_Y = 1;
const int RK_COORD_Z = 2;

//! dv/dt = f/m
const int RK_VEL = 3;
const int RK_VEL_X = 3;
const int RK_VEL_Y = 4;
const int RK_VEL_Z = 5;

//! d\omega/dt = T_c/m
const int RK_ANGULAR_VEL = 6;
const int RK_ANGULAR_VEL_X = 6;
const int RK_ANGULAR_VEL_Y = 7;
const int RK_ANGULAR_VEL_Z = 8;

//! dT/dt = F_T
const int RK_TEMP = 9;

}

class OnePointVariable
{
private:
    //! The global particle ID.
    int particleID;

    //! The local cell ID of these point.
    int cellID;

    //! To Get the bc of particle for current cell.
    int particleBC;

    //! isZoneParticle : 1 particle is on zone, 0 out of zone.
    //! This not include ghost and corner.
    bool isZoneParticle;

    //! Global zone index for current zone.
    int zoneIndex;

    //! cell Type as IndexParticleOfCell's INDEX_CELLTYPE.
    int cellType;

    //! ifOneRK : if only rk = 1 when particle touch wall.
    //!   0 -- rk4.
    //!   1 -- particle touch wall.
    int ifOneRK;

    //! particle diameter of sphere
    RDouble onePointDiameter;

    //! particle density
    RDouble onePointDensity;

    //! particle specificHeatCapacity
    RDouble onePointSpecificHeatCapacity;

    //! particle Emissivity
    RDouble onePointEmissivity;

    //! particle Re Number
    RDouble onePointParticleReNum;

    //! particle Ma Number
    RDouble onePointParticleMaNum;

    //! particle Nu Number
    RDouble onePointParticleNuseltNum;

    //! particle Stokes Number
    RDouble onePointTotalForceSt;
    RDouble onePointStokesSt;
    RDouble onePointModifiedSt;

    //! particle ConvectiveHeatTransfer
    RDouble onePointParticleConvectiveHeatTransfer;

    //! particle RadiativeHeatTransfer
    RDouble onePointParticleRadiativeHeatTransfer;

    //! particle Work
    RDouble onePointParticleWork;

    RDouble onePointParticleDissipation;

    //! Coordinate of point center
    SPDouble onePointCoordinate;

    //! Velocity of point center
    SPDouble onePointVelocity;

    //! Relative Velocity of point center
    SPDouble onePointRelativeVelocity;

    //! Velocity of point center
    SPDouble onePointAngularVelocity;

    //! Acceleration of point center for unsteady force.
    SPDouble onePointAcceleration;

    //! Particle Temperature.
    int nDimTem;
    SPDouble onePointTemperature;

    //! Flow variable on particle location.
    //! The index is INDEX_FLOWONPARTICLE.
    //! 0 -- flow velocity. 3
    //! 1 -- flow density. 1
    //! 2 -- flow dynamic viscosity. 1
    //! 3 -- flow temperature. 1
    //! 4 -- flow velocity grad. 3*3
    //! 5 -- flow pressure grad. 3
    //! 6 -- flow temperature grad. 3
    SPDouble flowOnParticle;

    //! The force from flow to particle.
    //! The index is INDEX_PARTICLE_FORCE.
    //! 0 -- Gravity force.
    //! 1 -- Drag force.
    //! 2 -- Saffman force.
    //! 3 -- Magnus force.
    //! 4 -- Addmass force.
    //! 5 -- FluidAcceleration force.
    //! 6 -- Brownian force.
    //! 7 -- Thermophoretic force.
    SPDouble particleForce;

    //! Value for RK.
    //! eg : RK4.
    //! The equation y = f(x) and f(x_0) = y_0
    //! k1 = f(x_0)
    //! k2 = f(x_0 + 0.5*h)
    //! k3 = f(x_0 + 0.5*h)
    //! k4 = f(x_0 + 1.0*h)
    //! f(x_1) = y_0 + (h/6) *( k1 + 2*k2 + 2*k3+k4)
    //! dRKValue is k.
    //! initRKValue is y_0.
    //! sumDRK is ( k1 + 2*k2 + 2*k3+k4).
    int nEquation;
    SPDouble initRKValue;
    SPDouble dRKValue;
    SPDouble sumDRK;

    //! The old value in RK and particle Track.
    //! initCellID : The init value of Cell id in RK.
    int oldCellID;
    //! onePointCoordinateOld : The coordinates when 
    //!   the particle position was last changed. 
    //!   This variable is used for the 
    //!   old coordinates when the particle is tracked
    SPDouble onePointCoordinateOld;

    RDouble particelCFL;

public:
    //! In this class, there is no new space, 
    //! while new of simplepoint is called automatically at the beginning of this constructor.
    OnePointVariable();

    OnePointVariable(int simulationRotate, int simulationTemperature,int simulationAcceleration);

    //! Since we do not have any new in this class, 
    //! the destructor of SimplePointer will be called automatically 
    //! at the end of the destructor. 
    //! There is no need to delete it manually.
    ~OnePointVariable();

public:

    void SetIsParticleInZone(bool &inZone);

    void SetParticleID(int &particleID);
    //! Set the local cell ID of this point.
    void SetCellID(int cellID);
    void SetOldCellID(int cellID);

    void SetParticleBCType(int parBCType);

    void SetIfOneRK(int value);

    void SetDiameter(RDouble diameter);
    void SetDensity(RDouble density);
    void SetSpecificHeatCapacity(RDouble specificHeatCapacity);
    void SetEmissivity(RDouble particleEmissivity);

    void SetCoordinate(RDouble *coordinateX,RDouble *coordinateY,RDouble *coordinateZ);
    void SetVelocity(RDouble *velocityX, RDouble *velocityY, RDouble *velocityZ);
    void SetRelativeVelocity(RDouble *relativeVelocityX, RDouble *relativeVelocityY, RDouble *relativeVelocityZ);
    void SetAngularVelocity(RDouble *angularVelocityX, RDouble *angularVelocityY, RDouble *angularVelocityZ);
    void SetAcceleration(RDouble *accelerationX, RDouble *accelerationY, RDouble *accelerationZ);
    void SetTemperature(RDouble *temperature);
    void SetForce(RDouble *particleForceX, RDouble *particleForceY, RDouble *particleForceZ, int &forceType);

    void SetParticleReNum(RDouble particleReNum);
    void SetParticleMaNum(RDouble particleMaNum);
    void SetParticleNuseltNum(RDouble particleNuseltNum);

    void SetParticleTotalForceSt(RDouble particleNuseltNum);
    void SetParticleStokesSt(RDouble particleNuseltNum);
    void SetParticleModifiedSt(RDouble particleNuseltNum);

    void SetParticleConvectiveHeatTransfer(RDouble particleConvectiveHeatTransfer);
    void SetParticleRadiativeHeatTransfer(RDouble particleRadiativeHeatTransfer);
    void SetParticleWork(RDouble particleWork);
    void SetParticleDissipation(RDouble particleDissipation);

    //! Set the old value as current value.
    void SetCoordinateOld();

    void SetParticleCFL(RDouble CFL);

    //! Dimensionless
    void SetDimensionlessLength(RDouble refLengthDimensional);
    void SetDimensionlessVelocity(RDouble refVelocityDimensional);
    void SetDimensionlessRelativeVelocity(RDouble refVelocityDimensional);
    void SetDimensionlessAngularVelocity(RDouble refAngularVelocity);
    void SetDimensionlessAcceleration(RDouble refAcceleration);
    void SetDimensionlessTemperature(RDouble refTemperature);
    void SetDimensionlessSpecificHeatCapacity(RDouble refVelocityDimensional, RDouble refTemperature);
    void SetDimensionlessDensity(RDouble refDensity);
    void SetDimensionlessForce(RDouble refForce);

    SPDouble *GetOnePointCoordinate();
    SPDouble *GetOnePointVelocity();
    SPDouble *GetOnePointRelativeVelocity();
    SPDouble *GetOnePointAcceleration();
    SPDouble *GetOnePointAngularVelocity();
    SPDouble *GetOnePointTemperature();

    SPDouble *GetOnePointFlowOnParticle();
    SPDouble *GetOnePointForce();

    SPDouble *GetOnePointCoordinateOld();

    RDouble GetOnePointDiameter();

    RDouble GetOnePointDensity();
    RDouble GetOnePointSpecificHeatCapacity();
    RDouble GetOnePointEmissivity();

    RDouble &GetOnePointCoordinateOnDirction(int iDim);
    RDouble &GetOnePointVelocityOnDirction(int iDim);
    RDouble &GetOnePointRelativeVelocityOnDirction(int iDim);
    RDouble &GetOnePointAngularVelocityOnDirction(int iDim);
    RDouble &GetOnePointAccelerationOnDirction(int iDim);
    RDouble &GetOnePointTemperatureOnDirction(int iDim);
    RDouble &GetOnePointForceOnDirction(int iDim,int forceType);

    SPDouble *GetOnePointInitRKValue();
    SPDouble *GetOnePointDRKValue();
    SPDouble *GetOnePointSumDRK();

    int GetParticleBCType();

    int GetIfOneRK();

    int GetNumDimTemperature();

    int GetParticleID();
    RDouble GetOnePointParticleReNum();
    RDouble GetOnePointParticleNuseltNum();
    RDouble GetOnePointParticleMaNum();

    RDouble GetOnePointTotalForceSt();
    RDouble GetOnePointStokesSt();
    RDouble GetOnePointModifiedSt();

    RDouble GetOnePointParticleConvectiveHeatTransfer();
    RDouble GetOnePointParticleRadiativeHeatTransfer();
    RDouble GetOnePointParticleWork();
    RDouble GetOnePointParticleDissipation();

    //! Get the local cell ID of this point.
    int GetCellID();
    int GetOldCellID();

    RDouble GetParticleCFL();

    bool IsParticleInZone();

    void Print2Window();
    string Print2WindowOS();

};

}