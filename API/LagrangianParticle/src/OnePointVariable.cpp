#include "OnePointVariable.h"
#include "GlobalDataBase.h"
#include "ParticleBoundaryCondition.h"

namespace PHSPACE
{

OnePointVariable::OnePointVariable()
{
    //! Note here that when entering the constructor at the beginning, 
    //! the constructor of SimplePointer will be called.
    this->cellID = 0;

    using namespace PARTICLEBCSPACE;
    this->particleBC = NO_PARTICLE_BC;

    this->isZoneParticle = false;
    this->zoneIndex = 0;
    this->ifOneRK = 0;

    this->onePointDiameter = 0.0;
    this->onePointDensity = 0.0;

    this->onePointSpecificHeatCapacity = 0.0;
    this->onePointParticleReNum = 0.0;
    this->onePointParticleMaNum = 0.0;
    this->onePointParticleNuseltNum = 2.0;
    this->onePointTotalForceSt=0.0;
    this->onePointStokesSt=0.0;
    this->onePointModifiedSt=0.0;
    this->onePointParticleConvectiveHeatTransfer = 0.0;

    //! particle RadiativeHeatTransfer
    this->onePointParticleRadiativeHeatTransfer = 0.0;

    //! particle Work
    this->onePointParticleWork = 0.0;
    this->onePointParticleDissipation = 0.0;

    int nDimVar = 3;
    RDouble initValue = 0.0;

    //! Set size of data and new data.
    this->onePointCoordinate.SetSize(nDimVar, initValue);
    this->onePointVelocity.SetSize(nDimVar, initValue);
    this->onePointRelativeVelocity.SetSize(nDimVar, initValue);
    this->onePointAngularVelocity.SetSize(nDimVar, initValue);
    this->onePointAcceleration.SetSize(nDimVar, initValue);

    this->nDimTem = 1;
    nDimVar = nDimTem;
    this->onePointTemperature.SetSize(nDimVar, initValue);

    nDimVar = INDEX_FLOWONPARTICLE::nVarOnFlow;
    this->flowOnParticle.SetSize(nDimVar, initValue);

    nDimVar = INDEX_PARTICLE_FORCE::nForceDir;
    this->particleForce.SetSize(nDimVar, initValue);

    //! dx/dt,dv/dt,d\omega/dt,dT/dt.
    nDimVar = 3+3+3+1;
    this->initRKValue.SetSize(nDimVar, initValue);
    this->dRKValue.SetSize(nDimVar, initValue);
    this->sumDRK.SetSize(nDimVar, initValue);

    this->oldCellID = 0;

    nDimVar = 3;
    this->onePointCoordinateOld.SetSize(nDimVar, initValue);

    this->particelCFL = 0.0;
}

OnePointVariable::OnePointVariable(int simulationRotate, int simulationTemperature, int simulationAcceleration)
{

}

OnePointVariable::~OnePointVariable()
{
    //! Destructors can cause many problems of deletion by mistake.
    //! We must be careful.
    //! Note that since we do not have any new in this class, 
    //! the destructor of simpleponter will be called automatically 
    //! at the end of the destructor. 
    //! There is no need to delete it manually.
}

void OnePointVariable::SetParticleID(int &particleID)
{
    this->particleID = particleID;
}

void OnePointVariable::SetCellID(int cellID)
{
    this->cellID = cellID;
}

void OnePointVariable::SetOldCellID(int cellID)
{
    this->oldCellID = cellID;
}

void OnePointVariable::SetCoordinate(RDouble *coordinateX, RDouble *coordinateY, RDouble *coordinateZ)
{
    int indexData = 0;
    int xDir = 0;
    int yDir = 1;
    int zDir = 2;
    if (coordinateX != nullptr && coordinateY != nullptr && coordinateZ != nullptr)
    {
        this->onePointCoordinate.CopyDataOnIndex(coordinateX, indexData, xDir);
        this->onePointCoordinate.CopyDataOnIndex(coordinateY, indexData, yDir);
        this->onePointCoordinate.CopyDataOnIndex(coordinateZ, indexData, zDir);
    }
}

void OnePointVariable::SetVelocity(RDouble *velocityX, RDouble *velocityY, RDouble *velocityZ)
{
    int indexData = 0;
    int xDir = 0;
    int yDir = 1;
    int zDir = 2;
    if (velocityX != nullptr && velocityY != nullptr && velocityZ != nullptr)
    {
        this->onePointVelocity.CopyDataOnIndex(velocityX, indexData, xDir);
        this->onePointVelocity.CopyDataOnIndex(velocityY, indexData, yDir);
        this->onePointVelocity.CopyDataOnIndex(velocityZ, indexData, zDir);
    }
}

void OnePointVariable::SetRelativeVelocity(RDouble *relativeVelocityX, RDouble *relativeVelocityY, RDouble *relativeVelocityZ)
{
    int indexData = 0;
    int xDir = 0;
    int yDir = 1;
    int zDir = 2;
    if (relativeVelocityX != nullptr && relativeVelocityY != nullptr && relativeVelocityZ != nullptr)
    {
        this->onePointRelativeVelocity.CopyDataOnIndex(relativeVelocityX, indexData, xDir);
        this->onePointRelativeVelocity.CopyDataOnIndex(relativeVelocityY, indexData, yDir);
        this->onePointRelativeVelocity.CopyDataOnIndex(relativeVelocityZ, indexData, zDir);
    }
}

void OnePointVariable::SetAngularVelocity(RDouble *angularVelocityX, RDouble *angularVelocityY, RDouble *angularVelocityZ)
{
    int indexData = 0;
    int xDir = 0;
    int yDir = 1;
    int zDir = 2;
    int nDimVar = 3;
    if (angularVelocityX != nullptr && angularVelocityY != nullptr && angularVelocityZ != nullptr)
    {
        this->onePointAngularVelocity.CopyDataOnIndex(angularVelocityX, indexData, xDir);
        this->onePointAngularVelocity.CopyDataOnIndex(angularVelocityY, indexData, yDir);
        this->onePointAngularVelocity.CopyDataOnIndex(angularVelocityZ, indexData, zDir);
    }
}

void OnePointVariable::SetAcceleration(RDouble *accelerationX, RDouble *accelerationY, RDouble *accelerationZ)
{
    int indexData = 0;
    int xDir = 0;
    int yDir = 1;
    int zDir = 2;
    this->onePointAcceleration.CopyDataOnIndex(accelerationX, indexData, xDir);
    this->onePointAcceleration.CopyDataOnIndex(accelerationY, indexData, yDir);
    this->onePointAcceleration.CopyDataOnIndex(accelerationZ, indexData, zDir);
    if ((accelerationX != nullptr) && (accelerationY != nullptr) && (accelerationZ != nullptr))
    {
        this->onePointAcceleration.CopyDataOnIndex(accelerationX, indexData, xDir);
        this->onePointAcceleration.CopyDataOnIndex(accelerationY, indexData, yDir);
        this->onePointAcceleration.CopyDataOnIndex(accelerationZ, indexData, zDir);
    }
}

void OnePointVariable::SetTemperature(RDouble *temperature)
{
    int indexData = 0;
    int xDir = 0;

    if (temperature != nullptr)
    {
        this->onePointTemperature.CopyDataOnIndex(temperature, indexData, xDir);
    }
}

void OnePointVariable::SetForce(RDouble *particleForceX, RDouble *particleForceY, RDouble *particleForceZ, int &forceType)
{
    int indexData = 0;

    int nDimVar = 3;
    int xDir = forceType * nDimVar + 0;
    int yDir = forceType * nDimVar + 1;
    int zDir = forceType * nDimVar + 2;

    this->particleForce.CopyDataOnIndex(particleForceX, indexData, xDir);
    this->particleForce.CopyDataOnIndex(particleForceY, indexData, yDir);
    this->particleForce.CopyDataOnIndex(particleForceZ, indexData, zDir);
}

void OnePointVariable::SetCoordinateOld()
{
    int nDimVar = 3;
    this->onePointCoordinateOld.CopyData(this->onePointCoordinate, nDimVar);
}

void OnePointVariable::SetParticleCFL(RDouble CFL)
{
    this->particelCFL = CFL;
}

void OnePointVariable::SetDiameter(RDouble diameter)
{
    this->onePointDiameter = diameter;
}


void OnePointVariable::SetDensity(RDouble density)
{
    this->onePointDensity = density;
}

void OnePointVariable::SetSpecificHeatCapacity(RDouble specificHeatCapacity)
{
    this->onePointSpecificHeatCapacity = specificHeatCapacity;
}

void OnePointVariable::SetEmissivity(RDouble particleEmissivity)
{
    this->onePointEmissivity = particleEmissivity;
}

void OnePointVariable::SetParticleReNum(RDouble particleReNum)
{
    this->onePointParticleReNum = particleReNum;
}
void OnePointVariable::SetParticleMaNum(RDouble particleMaNum)
{
    this->onePointParticleMaNum = particleMaNum;
}
void OnePointVariable::SetParticleNuseltNum(RDouble particleNuseltNum)
{
    this->onePointParticleNuseltNum = particleNuseltNum;
}
void OnePointVariable::SetParticleTotalForceSt(RDouble particleTotalForceSt)
{
    this->onePointTotalForceSt = particleTotalForceSt;
}
void OnePointVariable::SetParticleStokesSt(RDouble particleStokesSt)
{
    this->onePointStokesSt = particleStokesSt;
}
void OnePointVariable::SetParticleModifiedSt(RDouble particleModifiedSt)
{
    this->onePointModifiedSt = particleModifiedSt;
}
void OnePointVariable::SetParticleConvectiveHeatTransfer(RDouble particleConvectiveHeatTransfer)
{
    this->onePointParticleConvectiveHeatTransfer = particleConvectiveHeatTransfer;
}
void OnePointVariable::SetParticleRadiativeHeatTransfer(RDouble particleRadiativeHeatTransfer)
{
    this->onePointParticleRadiativeHeatTransfer = particleRadiativeHeatTransfer;
}
void OnePointVariable::SetParticleWork(RDouble particleWork)
{
    this->onePointParticleWork = particleWork;
}
void OnePointVariable::SetParticleDissipation(RDouble particleDissipation)
{
    this->onePointParticleDissipation = particleDissipation;
}

void OnePointVariable::SetDimensionlessLength(RDouble refLengthDimensional)
{
    int nDim = 3;
    RDouble unit = 1.0 / refLengthDimensional;
    this->onePointCoordinate.DotValue(unit, nDim);
    this->onePointDiameter /= refLengthDimensional;
}

void OnePointVariable::SetDimensionlessVelocity(RDouble refVelocityDimensional)
{
    int nDim = 3;
    RDouble unit = 1.0 / refVelocityDimensional;
    this->onePointVelocity.DotValue(unit, nDim);
}

void OnePointVariable::SetDimensionlessRelativeVelocity(RDouble refVelocityDimensional)
{
    int nDim = 3;
    RDouble unit = 1.0 / refVelocityDimensional;
    this->onePointRelativeVelocity.DotValue(unit, nDim);
}


void OnePointVariable::SetDimensionlessAngularVelocity(RDouble refAngularVelocity)
{
    int nDim = 3;
    RDouble unit = 1.0 / refAngularVelocity;
    this->onePointAngularVelocity.DotValue(unit, nDim);
}

void OnePointVariable::SetDimensionlessAcceleration(RDouble refAcceleration)
{
    int nDim = 3;
    RDouble unit = 1.0 / refAcceleration;
    this->onePointAcceleration.DotValue(unit, nDim);
}

void OnePointVariable::SetDimensionlessTemperature(RDouble refTemperature)
{
    int nDim = 1;
    RDouble unit = 1.0 / refTemperature;
    this->onePointTemperature.DotValue(unit, nDim);
}

void OnePointVariable::SetDimensionlessDensity(RDouble refDensity)
{
    this->onePointDensity /= refDensity;
}

void OnePointVariable::SetDimensionlessSpecificHeatCapacity(RDouble refVelocityDimensional,RDouble refTemperature)
{
    this->onePointSpecificHeatCapacity /= (refVelocityDimensional* refVelocityDimensional/ refTemperature);
}

void OnePointVariable::SetDimensionlessForce(RDouble refForce)
{
    int nDim = 8 * 3;
    RDouble unit = 1.0 / refForce;
    this->particleForce.DotValue(unit, nDim);
}

SPDouble *OnePointVariable::GetOnePointCoordinate()
{
    return &(this->onePointCoordinate);
}

SPDouble *OnePointVariable::GetOnePointCoordinateOld()
{
    return &(this->onePointCoordinateOld);
}

SPDouble *OnePointVariable::GetOnePointVelocity()
{
    return &(this->onePointVelocity);
}

SPDouble *OnePointVariable::GetOnePointRelativeVelocity()
{
    return &(this->onePointRelativeVelocity);
}

SPDouble *OnePointVariable::GetOnePointAcceleration()
{
    return &(this->onePointAcceleration);
}

SPDouble *OnePointVariable::GetOnePointAngularVelocity()
{
    return &(this->onePointAngularVelocity);
}

SPDouble *OnePointVariable::GetOnePointTemperature()
{
    return &(this->onePointTemperature);
}

RDouble OnePointVariable::GetOnePointSpecificHeatCapacity()
{
    return this->onePointSpecificHeatCapacity;
}

RDouble OnePointVariable::GetOnePointEmissivity()
{
    return this->onePointEmissivity;
}

SPDouble *OnePointVariable::GetOnePointFlowOnParticle()
{
    return &(this->flowOnParticle);
}

SPDouble *OnePointVariable::GetOnePointForce()
{
    return &(this->particleForce);
}

RDouble OnePointVariable::GetOnePointDiameter()
{
    return this->onePointDiameter;
}

RDouble OnePointVariable::GetOnePointParticleReNum()
{
    return this->onePointParticleReNum;
}

RDouble OnePointVariable::GetOnePointParticleMaNum()
{
    return this->onePointParticleMaNum;
}

RDouble OnePointVariable::GetOnePointParticleNuseltNum()
{
    return this->onePointParticleNuseltNum;
}

RDouble OnePointVariable::GetOnePointTotalForceSt()
{
    return this->onePointTotalForceSt;
}

RDouble OnePointVariable::GetOnePointStokesSt()
{
    return this->onePointStokesSt;
}

RDouble OnePointVariable::GetOnePointModifiedSt()
{
    return this->onePointModifiedSt;
}

RDouble OnePointVariable::GetOnePointParticleConvectiveHeatTransfer()
{
    return this->onePointParticleConvectiveHeatTransfer;
}

RDouble OnePointVariable::GetOnePointParticleRadiativeHeatTransfer()
{
    return this->onePointParticleRadiativeHeatTransfer;
}

RDouble OnePointVariable::GetOnePointParticleWork()
{
    return this->onePointParticleWork;
}

RDouble OnePointVariable::GetOnePointParticleDissipation()
{
    return this->onePointParticleDissipation;
}

RDouble OnePointVariable::GetOnePointDensity()
{
    return this->onePointDensity;
}

SPDouble *OnePointVariable::GetOnePointInitRKValue()
{
    return &(this->initRKValue);
}

SPDouble *OnePointVariable::GetOnePointDRKValue()
{
    return &(this->dRKValue);
}

SPDouble *OnePointVariable::GetOnePointSumDRK()
{
    return &(this->sumDRK);
}

RDouble &OnePointVariable::GetOnePointCoordinateOnDirction(int iDim)
{
    return onePointCoordinate[iDim];
}

RDouble &OnePointVariable::GetOnePointVelocityOnDirction(int iDim)
{
    return onePointVelocity[iDim];
}

RDouble &OnePointVariable::GetOnePointRelativeVelocityOnDirction(int iDim)
{
    return onePointRelativeVelocity[iDim];
}

RDouble &OnePointVariable::GetOnePointAngularVelocityOnDirction(int iDim)
{
    return onePointAngularVelocity[iDim];
}

RDouble &OnePointVariable::GetOnePointAccelerationOnDirction(int iDim)
{
    return onePointAcceleration[iDim];
}

RDouble &OnePointVariable::GetOnePointTemperatureOnDirction(int iDim)
{
    return onePointTemperature[iDim];
}

RDouble &OnePointVariable::GetOnePointForceOnDirction(int iDim, int forceType)
{
    return particleForce[forceType + iDim];
}

int OnePointVariable::GetParticleID()
{
    return this->particleID;
}



int OnePointVariable::GetCellID()
{
    return this->cellID;
}

int OnePointVariable::GetOldCellID()
{
    return this->oldCellID;
}


RDouble OnePointVariable::GetParticleCFL()
{
    return this->particelCFL;
}

void OnePointVariable::Print2Window()
{
    int nDim = 3;
    this->onePointCoordinate.PrintData2Window(nDim);
    this->onePointVelocity.PrintData2Window(nDim);
    this->onePointAngularVelocity.PrintData2Window(nDim);
    this->onePointAcceleration.PrintData2Window(nDim);
}

string OnePointVariable::Print2WindowOS()
{
    int nDimVar;
    ostringstream oss;
    oss << "  Coordinate : ";
    nDimVar = 3;
    for (int iDim = 0; iDim < nDimVar; ++iDim)
    {
        oss << " " << onePointCoordinate[iDim];
    }
    oss << "\n";

    oss << "  Velocity : ";
    nDimVar = 3;
    for (int iDim = 0; iDim < nDimVar; ++iDim)
    {
        oss << " " << onePointVelocity[iDim];
    }
    oss << "\n";

    oss << "  AngularVelocity : ";
    nDimVar = 3;
    for (int iDim = 0; iDim < nDimVar; ++iDim)
    {
        oss << " " << onePointAngularVelocity[iDim];
    }
    oss << "\n";

    oss << "  Acceleration : ";
    nDimVar = 3;
    for (int iDim = 0; iDim < nDimVar; ++iDim)
    {
        oss << " " << onePointAcceleration[iDim];
    }
    oss << "\n";

    oss << "  Temperature : ";
    nDimVar = this->nDimTem;
    for (int iDim = 0; iDim < nDimVar; ++iDim)
    {
        oss << " " << onePointTemperature[iDim];
    }
    oss << "\n";

    oss << "  Force : " << "\n";
    using namespace INDEX_PARTICLE_FORCE;
    nDimVar = 3;
    int nForce = nForceDir / nDimVar;
    for (int iForce = 0; iForce < nForce; ++iForce)
    {
        oss << "   iForce " << iForce << " : ";
        oss << particleForce[iForce * nDimVar + 0] << " ";
        oss << particleForce[iForce * nDimVar + 1] << " ";
        oss << particleForce[iForce * nDimVar + 2] << " ";
        oss << "\n";
    }

    oss << endl;
    return oss.str();
}

void OnePointVariable::SetIfOneRK(int value)
{
    this->ifOneRK = value;
}

void OnePointVariable::SetIsParticleInZone(bool &inZone)
{
    this->isZoneParticle = inZone;
}

bool OnePointVariable::IsParticleInZone()
{
    return this->isZoneParticle;
}

void OnePointVariable::SetParticleBCType(int parBCType)
{
    this->particleBC = parBCType;
}

int OnePointVariable::GetParticleBCType()
{
    return this->particleBC;
}

int OnePointVariable::GetIfOneRK()
{
    return this->ifOneRK;
}

int OnePointVariable::GetNumDimTemperature()
{
    return this->nDimTem;
}

};