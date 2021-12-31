#include "Param_NSSolver.h"
#include "GlobalDataBase.h"
#include "TK_Log.h"
#include "Constants.h"

namespace PHSPACE
{

LIB_EXPORT Param_NSSolver::Param_NSSolver()
{
    
}

LIB_EXPORT Param_NSSolver::~Param_NSSolver()
{

}

LIB_EXPORT void Param_NSSolver::Init()
{
    Param_CFDSolver::Init();
    intervalStepForce   = GlobalDataBase::GetIntParaFromDB("intervalStepForce");

    wallTemperature     = GlobalDataBase::GetDoubleParaFromDB("wallTemperature");
    AoA                 = GlobalDataBase::GetDoubleParaFromDB("attackd");
    angleSlide          = GlobalDataBase::GetDoubleParaFromDB("angleSlide");

    ifLowSpeedPrecon    = GlobalDataBase::GetIntParaFromDB("ifLowSpeedPrecon");

    roeEntropyFixMethod = GlobalDataBase::GetIntParaFromDB("roeEntropyFixMethod");
    roeEntropyScale     = GlobalDataBase::GetDoubleParaFromDB("roeEntropyScale");

    nNSEquation         = GlobalDataBase::GetIntParaFromDB("nm");
    nLaminar            = GlobalDataBase::GetIntParaFromDB("nl");
    nChemical           = GlobalDataBase::GetIntParaFromDB("nchem");
    nTemperatureModel   = GlobalDataBase::GetIntParaFromDB("ntmodel");
    numberOfSpecies     = GlobalDataBase::GetIntParaFromDB("numberOfSpecies");

    nDensityModify      = GlobalDataBase::GetIntParaFromDB("nDensityModify");
    nDebug              = GlobalDataBase::GetIntParaFromDB("nDebug");
    nTEnergyModel       = GlobalDataBase::GetIntParaFromDB("nTEnergyModel");

    nEnergyRecycle      = GlobalDataBase::GetIntParaFromDB("nEnergyRecycle") * nChemical;

    nTurblenceForChemical = GlobalDataBase::GetIntParaFromDB("nTurblenceForChemical");

    nViscosityFluxSublevelModified = GlobalDataBase::GetIntParaFromDB("nViscosityFluxSublevelModified");

    if (nChemical == 0)
    {
        if (nTemperatureModel > 1)
        {
            nTemperatureModel = 1;
            GlobalDataBase::UpdateData("ntmodel", &nTemperatureModel, PHINT, 1);
        }
        if (nLaminar > nNSEquation)
        {
            nLaminar = nNSEquation;
            GlobalDataBase::UpdateData("nl", &nLaminar, PHINT, 1);
        }
        if (numberOfSpecies > 0)
        {
            numberOfSpecies = 0;
            GlobalDataBase::UpdateData("numberOfSpecies", &numberOfSpecies, PHINT, 1);
        }
    }

    prandtlLaminar    = GlobalDataBase::GetDoubleParaFromDB("prl");
    prandtlTurbulence = GlobalDataBase::GetDoubleParaFromDB("prt");

    oPrandtlLaminar    = GlobalDataBase::GetDoubleParaFromDB("oprl");
    oPrandtlTurbulence = GlobalDataBase::GetDoubleParaFromDB("oprt");

    nChemicalRadius = GlobalDataBase::GetIntParaFromDB("nchemrad");
    nChemicalSource = GlobalDataBase::GetIntParaFromDB("nchemsrc");

    timeIntegration = GlobalDataBase::GetIntParaFromDB("tscheme");

    nonDimensionalSutherlandTemperature = GlobalDataBase::GetDoubleParaFromDB("tsuth");

    mgProlongationType = GlobalDataBase::GetIntParaFromDB("mprol");

    if (nChemical > 0)
    {
        indexOfNitrogen = GlobalDataBase::GetIntParaFromDB("nNitrogenIndex");
        indexOfElectron = GlobalDataBase::GetIntParaFromDB("nElectronIndex");
        catalyticMassFraction = reinterpret_cast<RDouble *>(GlobalDataBase::GetDataPtr("catalyticMassFraction"));
    }
    else
    {
        indexOfNitrogen = -1;
        indexOfElectron = -1;
        catalyticMassFraction = 0;
    }
    nEquilibriumGas = GlobalDataBase::GetIntParaFromDB("nEquilibriumGas");

    primitiveVarFarfield = reinterpret_cast <RDouble *> (GlobalDataBase::GetDataPtr("prim_inf"));

    mTT = 0;
    if(nTemperatureModel == 3)
    {
        mTV = 1;
        mTE = 2;
    }
    else if(nTemperatureModel ==2)
    {
        mTV = 1;
        mTE = 1;
    }
    else
    {
        mTV = 0;
        mTE = 0;
    }
}

void Param_NSSolver::SetEntropyFixCoefficients()
{
    RDouble entropyAcoustic = 0.01, entropyConvective = 0.01;
    double referenceMachNumber = this->GetRefMachNumber();

    int ifLowSpeedPrecon = this->GetIfLowSpeedPrecon();
    if (ifLowSpeedPrecon)
    {
        //! Incompressible, precondition.
        entropyAcoustic   = 1.0e-8;
    }
    else
    {
        //! Compressible.
        RDouble x0, x1, y0, y1;
        if(referenceMachNumber <= 0.6)
        {
            entropyAcoustic = 0.0001;
        }
        else if(referenceMachNumber <= 2.0)
        {
            x0 = 0.6;    x1 = 2.0;
            y0 = 0.0001; y1 = 0.01;
            entropyAcoustic = y0 * (referenceMachNumber - x1) / (x0 - x1) + y1 * (referenceMachNumber - x0) / (x1 - x0);
        }
        else if(referenceMachNumber <= 5.0)
        {
            x0 = 2.0;  x1 = 5.0;
            y0 = 0.01; y1 = 0.5;
            entropyAcoustic = y0 * (referenceMachNumber - x1) / (x0 - x1) + y1 * (referenceMachNumber - x0) / (x1 - x0);
        }
        else
        {
            entropyAcoustic = 1.0;
        }  

        //! Do not scale the entropy for incompressible precondition.
        entropyAcoustic *= this->roeEntropyScale;
    }

    entropyConvective = entropyAcoustic;

    this->RoeEntropyFixCoef1 = entropyAcoustic;
    this->RoeEntropyFixCoef2 = entropyConvective;
}

}
