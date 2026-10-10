#include "Param_ParticleSolver.h"
#include "GlobalDataBase.h"
#include "TK_Log.h"
#include "TK_Exit.h"
#include "TK_Parse.h"
#include "Constants.h"

namespace PHSPACE
{

LIB_EXPORT Param_ParticleSolver::Param_ParticleSolver()
{
    ifInitParticleAsFlow = 0;

    initOutStep = 0;

    iParticleModel = 0;

    initParticleFile = "";
    initParticleBCFile = "";

    ifReadRestartParticle = 0;
    ifWirteRestartParticle = 0;
    ifInitParticleBC = 0;

    gridTypeForPreAndPost = 0;

    version = 0.0;
}

LIB_EXPORT Param_ParticleSolver::~Param_ParticleSolver()
{
    //! Destructors can cause many problems of deletion by mistake.
    //! We must be careful. You can comment out this sentence 
    //! after the program debugging is perfect.
    //cout << "Call destructor function for Param_ParticleSolver" << endl;
}

LIB_EXPORT void Param_ParticleSolver::Init()
{
    Param_CFDSolver::Init();

    //! Param in particle_para.hypara in GlobalDataBase.
    if (GlobalDataBase::IsExist("iParticleModel", PHINT, 1))
    {
        iParticleModel = GlobalDataBase::GetIntParaFromDB("iParticleModel");
    }
    else
    {
        iParticleModel = 0;
    }

    initParticleFile = GlobalDataBase::GetStrParaFromDB("initParticleFile");
    ifReadRestartParticle = GlobalDataBase::GetIntParaFromDB("ifReadRestartParticle");
    ifWirteRestartParticle = GlobalDataBase::GetIntParaFromDB("ifWirteRestartParticle");

    if (GlobalDataBase::IsExist("iParticleModel", PHINT, 1))
    {
        ifInitParticleAsFlow = GlobalDataBase::GetIntParaFromDB("ifInitParticleAsFlow");
    }
    else
    {
        ifInitParticleAsFlow = 0;
    }

    if (GlobalDataBase::IsExist("gridTypeForPreAndPost", PHINT, 1))
    {
        gridTypeForPreAndPost = GlobalDataBase::GetIntParaFromDB("gridTypeForPreAndPost");
    }
    else
    {
        gridTypeForPreAndPost = 0;
    }

    if (GlobalDataBase::IsExist("particleHeatCD", PHDOUBLE, 1))
    {
        particleHeatCD = GlobalDataBase::GetDoubleParaFromDB("particleHeatCD");
    }
    else
    {
        particleHeatCD = 0.0;
    }

    if (GlobalDataBase::IsExist("radiativeTemperature", PHDOUBLE, 1))
    {
        radiativeTemperature = GlobalDataBase::GetDoubleParaFromDB("radiativeTemperature");
    }
    else
    {
        radiativeTemperature = 0.0;
    }

    if (GlobalDataBase::IsExist("enableRadiativeHeatTransfer", PHINT, 1))
    {
        enableRadiativeHeatTransfer = GlobalDataBase::GetIntParaFromDB("enableRadiativeHeatTransfer");
    }
    else
    {
        enableRadiativeHeatTransfer = 0;
    }


    if (GlobalDataBase::IsExist("enableVariableDiameter", PHINT, 1))
    {
        enableVariableDiameter = GlobalDataBase::GetIntParaFromDB("enableVariableDiameter");
    }
    else
    {
        enableVariableDiameter = 0;
    }

    if (GlobalDataBase::IsExist("ifInitParticleBC", PHDOUBLE, 1))
    {
        ifInitParticleBC = GlobalDataBase::GetIntParaFromDB("ifInitParticleBC");

        if (ifInitParticleBC != 0)
        {
            initParticleBCFile = GlobalDataBase::GetStrParaFromDB("initParticleBCFile");
        }

        if (3 == ifInitParticleBC)
        {
            initParticleBCPeriodicStep = GlobalDataBase::GetIntParaFromDB("initParticleBCPeriodicStep");
            particleRandomReleaseMode = GlobalDataBase::GetIntParaFromDB("particleRandomReleaseMode");
            particleRandomDiameterMode = GlobalDataBase::GetIntParaFromDB("particleRandomDiameterMode");

            randomDiameterMaxRatio = GlobalDataBase::GetDoubleParaFromDB("randomDiameterMaxRatio");
            randomDiameterMinRatio = GlobalDataBase::GetDoubleParaFromDB("randomDiameterMinRatio");


            initParticleBCRandomDirectionalRatio = GlobalDataBase::GetDoubleParaFromDB("initParticleBCRandomDirectionalRatio");
            initParticleBCRandomValueRatio = GlobalDataBase::GetDoubleParaFromDB("initParticleBCRandomValueRatio");


            if (initParticleBCRandomDirectionalRatio < 0)
            {
                TK_Exit::UnexpectedVarValue("initParticleBCRandomDirectionalRatio", initParticleBCRandomDirectionalRatio);
            }

            if (initParticleBCRandomValueRatio < 0)
            {
                TK_Exit::UnexpectedVarValue("initParticleBCRandomValueRatio", initParticleBCRandomValueRatio);
            }
        }

    }
    else
    {
        ifInitParticleBC = 0;
    }

    version = GlobalDataBase::GetDoubleParaFromDB("version");

    intervalStepRestartParticle = GlobalDataBase::GetIntParaFromDB("intervalStepRestartParticle");

    this->InitDimensionlessParam();

    flowHeatCPConstPressureDim = (1.0 / (this->refGama - 1.0)) * this->refAverageGeneralGasConstantDimensional;
    flowHeatCPConstPressureDimless = 1.0 / (this->refGama - 1.0) / refMachNumber / refMachNumber;

    prl = GlobalDataBase::GetDoubleParaFromDB("prl");

    //! The param from NSSolver.
    //! nm: Equation number of the physics, but is out of commision now.
    //! 4 -- for 2D,5 -- for 3D,but usually for both 2D and 3D,nm = 5.
    nNSEquation = GlobalDataBase::GetIntParaFromDB("nm");
    //! This param is not read from cfd_para.hypara,which is init on PerfectGas::InitNumberOfSpecies(),
    //! and in that function ,nl is equal to nm.
    nLaminar = GlobalDataBase::GetIntParaFromDB("nl");
    //! nchem:0 -- without chemical reaction flow.1 -- the chemical reaction flow is considered.
    nChemical = GlobalDataBase::GetIntParaFromDB("nchem");
    //! ntmodel: The thermodynamic temperature model.
    //!              1 -- One-temperature model.
    //!              2 -- Two-temperature model.
    //!              3 -- Three-temperature model.
    nTemperatureModel = GlobalDataBase::GetIntParaFromDB("ntmodel");
    if (nChemical == 0)
    {
        if (nTemperatureModel > 1)
        {
            nTemperatureModel = 1;
            GlobalDataBase::UpdateData("ntmodel", &nTemperatureModel, PHINT, 1);
        }
    }
}

bool Param_ParticleSolver::ifInitParticleInFlow()
{
    bool returnBool = false;
    if (1 == ifInitParticleAsFlow)
    {
        returnBool = true;
    }
    return returnBool;
}

void Param_ParticleSolver::SetInitOutStep(int initOutStep)
{
    this->initOutStep = initOutStep;
}

int Param_ParticleSolver::GetInitOutStep()
{
    return initOutStep;
}

bool Param_ParticleSolver::ifIsOrthogonalForPreAndPost()
{
    bool isOrthogonalGrid = false;
    if (1 == this->gridTypeForPreAndPost)
    {
        isOrthogonalGrid = true;
    }
    return isOrthogonalGrid;
}

bool Param_ParticleSolver::ifAddParticleInitBC()
{
    bool addParticle = false;
    if (3 == this-> ifInitParticleBC)
    {
        addParticle = true;
    }
    return addParticle;
}

int Param_ParticleSolver::GetInitParticleBCPeriodicStep()
{
    return this->initParticleBCPeriodicStep;
}

int Param_ParticleSolver::GetParticleRandomReleaseMode()
{
    return this->particleRandomReleaseMode;
}

int Param_ParticleSolver::GetParticleRandomParticleDiameterMode()
{
    return this->particleRandomDiameterMode;
}

RDouble Param_ParticleSolver::GetInitParticleBCRandomDirectionalRatio()
{
    return this->initParticleBCRandomDirectionalRatio;
}


RDouble Param_ParticleSolver::GetRandomDiameterMaxRatio()
{
    return this->randomDiameterMaxRatio;
}

RDouble Param_ParticleSolver::GetRandomDiameterMinRatio()
{
    return this->randomDiameterMinRatio;
}


RDouble Param_ParticleSolver::GetInitParticleBCRandomValueRatio()
{
    return this->initParticleBCRandomValueRatio;
}

int Param_ParticleSolver::GetiParticleModel()
{
    return this->iParticleModel;
}

string Param_ParticleSolver::GetInitParticleFile()
{
    return this->initParticleFile;
}

string Param_ParticleSolver::GetInitParticleBCFile()
{
    return this->initParticleBCFile;
}

RDouble Param_ParticleSolver::GetVersion()
{
    return this->version;
}

int Param_ParticleSolver::GetIntervalStepRestartParticle()
{
    return this->intervalStepRestartParticle;
}

int Param_ParticleSolver::GetNSEquationNumber()
{
    return this->nNSEquation;
}

int Param_ParticleSolver::GetLaminarNumber()
{
    return this->nLaminar;
}

int Param_ParticleSolver::GetChemicalFlag()
{
    return this->nChemical;
}

int Param_ParticleSolver::GetTemperatureModel()
{
    return this->nTemperatureModel;
}

bool Param_ParticleSolver::IfInitParticle()
{
    bool ifInit = false;
    if (0 == ifReadRestartParticle)
    {
        return ifInit = true;
    }
    return ifInit;
}

bool Param_ParticleSolver::IfReadRestartParticle()
{
    bool ifRestart = false;
    if (1 == ifReadRestartParticle)
    {
        ifRestart = true;
    }
    return ifRestart;
}

bool Param_ParticleSolver::IfWirteRestartParticle()
{
    bool ifWirte = false;
    if (1 == ifWirteRestartParticle)
    {
        ifWirte =  true;
    }
    return ifWirte;
}

bool Param_ParticleSolver::ifInitParticleBCInfo()
{
    bool ifInit = false;
    if (this->ifInitParticleBC != 0)
    {
        ifInit = true;
    }
    return ifInit;
}

void Param_ParticleSolver::InitDimensionlessParam()
{
    //! ================== ref dimensional param for flow =====================
    //! Flow field parameters are used to dimensionless particle parameters.
    //! length
    this->refReynoldsLengthDimensional = GlobalDataBase::GetDoubleParaFromDB("reynoldsReferenceLengthDimensional");

    //! velocity
    this->refVelocityDimensional = GlobalDataBase::GetDoubleParaFromDB("refDimensionalVelocity");

    //! temperature
    this->refTemperatureDimensional = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");

    //! pressure
    this->refPressureDimensional = GlobalDataBase::GetDoubleParaFromDB("refDimensionalPressure");

    //! density
    this->refDensityDimensional = GlobalDataBase::GetDoubleParaFromDB("refDimensionalDensity");

    //! R
    this->refAverageGeneralGasConstantDimensional = GlobalDataBase::GetDoubleParaFromDB("reference_average_general_gas_constant");

    this->refMachNumber = GlobalDataBase::GetDoubleParaFromDB("refMachNumber");

    this->refReNumber = GlobalDataBase::GetDoubleParaFromDB("refReNumber");

    this->refGama = GlobalDataBase::GetDoubleParaFromDB("refGama");

    this->refDynamicViscosityDimensional = GlobalDataBase::GetDoubleParaFromDB("refDynamicViscosityDimensional");

    this->refMassDimensional = refReynoldsLengthDimensional * refReynoldsLengthDimensional * refReynoldsLengthDimensional * refDensityDimensional;

    this->refTimeDimensional = refReynoldsLengthDimensional / refVelocityDimensional;

    //! T^{-1}
    this->refAngularVelocityDimensional = 1.0 / refTimeDimensional;

    this->refAccelerationyDimensional = refVelocityDimensional / refTimeDimensional;

    this->refForceDimensional = refMassDimensional * refAccelerationyDimensional;
}

RDouble Param_ParticleSolver::GetRefDimLength()
{
    return this->refReynoldsLengthDimensional;
}

RDouble Param_ParticleSolver::GetRefDimVelocity()
{
    return this->refVelocityDimensional;
}

RDouble Param_ParticleSolver::GetRefDimTemperature()
{
    return this->refTemperatureDimensional;
}

RDouble Param_ParticleSolver::GetRefDimPressure()
{
    return this->refPressureDimensional;
}

RDouble Param_ParticleSolver::GetRefDimDensity()
{
    return this->refDensityDimensional;
}

RDouble Param_ParticleSolver::GetRefDimR()
{
    return this->refAverageGeneralGasConstantDimensional;
}

RDouble Param_ParticleSolver::GetRefDimMass()
{
    return this->refMassDimensional;
}

RDouble Param_ParticleSolver::GetRefDimTime()
{
    return this->refTimeDimensional;
}

RDouble Param_ParticleSolver::GetRefDimAngularVelocity()
{
    return this->refAccelerationyDimensional;
}

RDouble Param_ParticleSolver::GetRefDimAcceleration()
{
    return this->refAccelerationyDimensional;
}

RDouble Param_ParticleSolver::GetRefMach()
{
    return this->refMachNumber;
}

RDouble Param_ParticleSolver::GetRefReIn()
{
    return this->refReNumber;
}

RDouble Param_ParticleSolver::GetRefGamma()
{
    return this->refGama;
}

RDouble Param_ParticleSolver::GetRefDimDynamicViscosity()
{
    return this->refDynamicViscosityDimensional;
}

RDouble Param_ParticleSolver::GetRefDimForce()
{
    return this->refForceDimensional;
}

void Param_ParticleSolver::SetNumParticleTotal(int nParticleTotal)
{
    this->nParticleTotal = nParticleTotal;
}

void Param_ParticleSolver::SetNumParticleLocal(int nParticleLocal)
{
    this->nParticleLocal = nParticleLocal;
}

void Param_ParticleSolver::SetNumParticleGroup(int &nGroup)
{
    this->nParticleGroup = nGroup;
}

int Param_ParticleSolver::GetNumParticleTotal()
{
    return this->nParticleTotal;
}

int Param_ParticleSolver::GetNumParticleLocal()
{
    return this->nParticleLocal;
}

int Param_ParticleSolver::GetNumParticleGroup()
{
    return this->nParticleGroup;
}

int Param_ParticleSolver::GetInitParticleBCInfoType()
{
    return this->ifInitParticleBC;
}

void Param_ParticleSolver::GetRefDataParam(Data_Param *parameter)
{
    parameter->UpdateData("refReynoldsLengthDimensional", &refReynoldsLengthDimensional,PHDOUBLE,1);
    parameter->UpdateData("refVelocityDimensional", &refVelocityDimensional, PHDOUBLE, 1);
    parameter->UpdateData("refTemperatureDimensional", &refTemperatureDimensional, PHDOUBLE, 1);
    parameter->UpdateData("refPressureDimensional", &refPressureDimensional, PHDOUBLE, 1);
    parameter->UpdateData("refDensityDimensional", &refDensityDimensional, PHDOUBLE, 1);
    parameter->UpdateData("refAverageGeneralGasConstantDimensional", &refAverageGeneralGasConstantDimensional, PHDOUBLE, 1);
    parameter->UpdateData("refDynamicViscosityDimensional", &refDynamicViscosityDimensional, PHDOUBLE, 1);

    parameter->UpdateData("refMachNumber", &refMachNumber, PHDOUBLE, 1);
    parameter->UpdateData("refReNumber", &refReNumber, PHDOUBLE, 1);
    parameter->UpdateData("refGama", &refGama, PHDOUBLE, 1);

    parameter->UpdateData("refMassDimensional", &refMassDimensional, PHDOUBLE, 1);
    parameter->UpdateData("refTimeDimensional", &refTimeDimensional, PHDOUBLE, 1);
    parameter->UpdateData("refAngularVelocityDimensional", &refAngularVelocityDimensional, PHDOUBLE, 1);
    parameter->UpdateData("refAccelerationyDimensional", &refAccelerationyDimensional, PHDOUBLE, 1);
    parameter->UpdateData("refForceDimensional", &refForceDimensional, PHDOUBLE, 1);

    parameter->UpdateData("particleHeatCD", &particleHeatCD, PHDOUBLE, 1);
    parameter->UpdateData("radiativeTemperature", &radiativeTemperature, PHDOUBLE, 1);
    parameter->UpdateData("enableRadiativeHeatTransfer", &enableRadiativeHeatTransfer, PHINT, 1);
    parameter->UpdateData("enableVariableDiameter", &enableVariableDiameter, PHINT, 1);

    parameter->UpdateData("flowHeatCPConstPressureDim", &flowHeatCPConstPressureDim, PHDOUBLE, 1);
    parameter->UpdateData("flowHeatCPConstPressureDimless", &flowHeatCPConstPressureDimless, PHDOUBLE, 1);

    parameter->UpdateData("prl", &prl, PHDOUBLE, 1);
}

}