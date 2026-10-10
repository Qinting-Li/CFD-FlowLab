#include "Param_ParticlePointSolverParallelDomain.h"
#include "Data_Param.h"
#include "TK_Parse.h"
#include "TK_Exit.h"
#include "Constants.h"
#include "GlobalDataBase.h"
#include "ParticleAnalysis.h"
using namespace std;

namespace PHSPACE
{

LIB_EXPORT Param_ParticlePointSolverParallelDomain::Param_ParticlePointSolverParallelDomain()
{
    simulationSplashFunction = 0;
    simulationPointRotate = 0;
    simulationParticleTemperature = 0;
    simulationBackCouping = 0;

    particleInterpolation_struct = "";

    nInitGlobalParticle = 0;
    nParticleGroup = 0;

    forceDragType = 0;
    forceMagnusType = 0;
    forceBrownianType = 0;
    forceThermophoreticType = 0;
    forceSaffmanType = 0;
    forceFluidAccelerationType = 0;
    forceAddedMassType = 0;
    forceBassetType = 0;

    nDimTem = 1;
    particleTemperatureModel = 0;
}

LIB_EXPORT Param_ParticlePointSolverParallelDomain::~Param_ParticlePointSolverParallelDomain()
{
    //! Destructors can cause many problems of deletion by mistake.
    //! We must be careful. You can comment out this sentence 
    //! after the program debugging is perfect.
}

LIB_EXPORT void Param_ParticlePointSolverParallelDomain::Init()
{
    Param_ParticleSolver::Init();
}

void Param_ParticlePointSolverParallelDomain::ReadParameter()
{
    Data_Param *parameter = new Data_Param;
    string paraFileName = "./bin/particlepoint_para.hypara";

    using namespace PARTICLE_PARAM;
    ReadParticleParamFile(parameter, paraFileName);

    //! basic param of particle.
    simulationSplashFunction = GetParticleIntPara(parameter, "simulationSplashFunction");
    simulationPointRotate = GetParticleIntPara(parameter, "simulationPointRotate");
    simulationParticleTemperature = GetParticleIntPara(parameter, "simulationParticleTemperature");
    simulationBackCouping = GetParticleIntPara(parameter, "simulationBackCouping");

    particleInterpolation_struct = GetParticleStrPara(parameter,"particleInterpolation_struct");
    particleInterpolation_unstruct = GetParticleStrPara(parameter,"particleInterpolation_unstruct");

    //! Force type.
    forceGravityType = GetParticleIntPara(parameter, "forceGravityType");
    forceDragType = GetParticleIntPara(parameter, "forceDragType");
    forceMagnusType = GetParticleIntPara(parameter, "forceMagnusType");
    forceSaffmanType = GetParticleIntPara(parameter, "forceSaffmanType");
    forceFluidAccelerationType = GetParticleIntPara(parameter, "forceFluidAccelerationType");
    forceAddedMassType = GetParticleIntPara(parameter, "forceAddedMassType");
    forceBassetType = GetParticleIntPara(parameter, "forceBassetType");
    forceBrownianType = GetParticleIntPara(parameter, "forceBrownianType");
    forceThermophoreticType = GetParticleIntPara(parameter, "forceThermophoreticType");

    particleTemperatureModel = GetParticleIntPara(parameter, "particleTemperatureModel");

    fileParticleOutputType = GetParticleIntPara(parameter, "fileParticleOutputType");

    nRK = GetParticleIntPara(parameter, "nRK");

    this->nEquationParticle = 3 + 3 + simulationPointRotate * 3 + simulationParticleTemperature * 1;

    parameter->GetData("domainMesh", &(this->domainMesh), PHDOUBLE, 6);

    delete parameter;
}


void Param_ParticlePointSolverParallelDomain::GetNumOfParticleGroup()
{

}

int Param_ParticlePointSolverParallelDomain::GetNumOfInitGlobalParticle()
{
    return nInitGlobalParticle;
}

int Param_ParticlePointSolverParallelDomain::GetForceGravityType()
{
    return this->forceGravityType;
}

int Param_ParticlePointSolverParallelDomain::GetForceDragType()
{
    return forceDragType;
}

int Param_ParticlePointSolverParallelDomain::GetForceMagnusType()
{
    return forceMagnusType;
}
int Param_ParticlePointSolverParallelDomain::GetForceBrownianType()
{
    return forceBrownianType;
}

int Param_ParticlePointSolverParallelDomain::GetForceThermophoreticType()
{
    return forceThermophoreticType;
}

int Param_ParticlePointSolverParallelDomain::GetForceSaffmanType()
{
    return forceSaffmanType;
}

int Param_ParticlePointSolverParallelDomain::GetForceAddedMassType()
{
    return forceAddedMassType;
}

int Param_ParticlePointSolverParallelDomain::GetForceBassetType()
{
    return forceBassetType;
}
int Param_ParticlePointSolverParallelDomain::GetFluidAccelerationType()
{
    return forceFluidAccelerationType;
}

int Param_ParticlePointSolverParallelDomain::GetNumRK()
{
    return this->nRK;
}

int Param_ParticlePointSolverParallelDomain::GetNumDimTemperature()
{
    return this->nDimTem;
}

RDouble *Param_ParticlePointSolverParallelDomain::GetDomainMesh()
{
    return this->domainMesh;
}

string Param_ParticlePointSolverParallelDomain::GetParticleInterpolationStruct()
{
    return this->particleInterpolation_struct;
}

string Param_ParticlePointSolverParallelDomain::GetParticleInterpolationUnstruct()
{
    return this->particleInterpolation_unstruct;
}

int Param_ParticlePointSolverParallelDomain::GetNumOfParticleEquation()
{
    return this->nEquationParticle;
}

void Param_ParticlePointSolverParallelDomain::GetForceTypeDataParam(Data_Param *parameter)
{
    parameter->UpdateData("forceGravityType", &forceGravityType, PHINT, 1);
    parameter->UpdateData("forceDragType", &forceDragType, PHINT, 1);

    parameter->UpdateData("forceSaffmanType", &forceSaffmanType, PHINT, 1);
    parameter->UpdateData("forceMagnusType", &forceMagnusType, PHINT, 1);

    parameter->UpdateData("forceAddedMassType", &forceAddedMassType, PHINT, 1);
    parameter->UpdateData("forceFluidAccelerationType", &forceFluidAccelerationType, PHINT, 1);

    parameter->UpdateData("forceBrownianType", &forceBrownianType, PHINT, 1);
    parameter->UpdateData("forceThermophoreticType", &forceThermophoreticType, PHINT, 1);

    parameter->UpdateData("particleTemperatureModel", &particleTemperatureModel, PHINT, 1);
}

bool Param_ParticlePointSolverParallelDomain::ifUseParticleSplash()
{
    bool useParticleSplash = false;
    if (1 == this->simulationSplashFunction)
    {
        useParticleSplash = true;
    }
    return useParticleSplash;
}

bool Param_ParticlePointSolverParallelDomain::ifUseParticleRotate()
{
    bool useParticleRotate = false;
    if (1 == this->simulationPointRotate)
    {
        useParticleRotate = true;
    }
    return useParticleRotate;
}

bool Param_ParticlePointSolverParallelDomain::ifUseParticleTemperature()
{
    bool useParticleTemperature = false;
    if (1 == this->simulationParticleTemperature)
    {
        useParticleTemperature = true;
    }
    return useParticleTemperature;
}

bool Param_ParticlePointSolverParallelDomain::ifUseParticleBackCouping()
{
    bool useParticleBackCouping = false;
    if (1 == this->simulationBackCouping)
    {
        useParticleBackCouping = true;
    }
    return useParticleBackCouping;
}

bool Param_ParticlePointSolverParallelDomain::ifUseParticleAcceleration()
{
    bool useParticleAcceleration = false;
    if (1 == this->simulationParticleAcceleration)
    {
        useParticleAcceleration = true;
    }
    return useParticleAcceleration;
}

void Param_ParticlePointSolverParallelDomain::SetNumDimTemperature(int nDimTemp)
{
    this->nDimTem = nDimTemp;
}


}