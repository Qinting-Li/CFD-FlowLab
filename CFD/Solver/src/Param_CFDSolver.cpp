#include "Param_CFDSolver.h"
#include "GlobalDataBase.h"
#include "Constants.h"

namespace PHSPACE
{

LIB_EXPORT Param_CFDSolver::Param_CFDSolver()
{

}

LIB_EXPORT Param_CFDSolver::~Param_CFDSolver()
{

}

LIB_EXPORT void Param_CFDSolver::Init()
{
    restartNSFile            = GlobalDataBase::GetStrParaFromDB("restartNSFile");
    resSaveFile              = GlobalDataBase::GetStrParaFromDB("resSaveFile");
    viscousName              = GlobalDataBase::GetStrParaFromDB("viscousName");

    plotFieldType            = GlobalDataBase::GetIntParaFromDB("plotFieldType");

    viscousType              = GlobalDataBase::GetIntParaFromDB("viscousType");
    isUnsteady               = GlobalDataBase::GetIntParaFromDB("iunsteady");

    intervalStepRes          = GlobalDataBase::GetIntParaFromDB("intervalStepRes");
    intervalStepPlot         = GlobalDataBase::GetIntParaFromDB("intervalStepPlot");
    intervalStepFlow         = GlobalDataBase::GetIntParaFromDB("intervalStepFlow");

    ifLocalTimeStep          = GlobalDataBase::GetIntParaFromDB("ifLocalTimeStep");

    CFLVaryStep              = GlobalDataBase::GetIntParaFromDB("CFLVaryStep");
    nLUSGSSweeps             = GlobalDataBase::GetIntParaFromDB("nLUSGSSweeps");
    nMGLevel                 = GlobalDataBase::GetIntParaFromDB("nMGLevel");

    ifStartFromSteadyResults = GlobalDataBase::GetIntParaFromDB("ifStartFromSteadyResults");
    ifStaticsFlowField       = GlobalDataBase::GetIntParaFromDB("ifStaticsFlowField");
    ifStaticsReynoldsStress  = 0;
    if (GlobalDataBase::IsExist("ifStaticsReynoldsStress", PHINT, 1))
    {
        ifStaticsReynoldsStress = GlobalDataBase::GetIntParaFromDB("ifStaticsReynoldsStress");
    }
    statisticMethod = 0;
    if (GlobalDataBase::IsExist("statisticMethod", PHINT, 1))
    {
        statisticMethod = GlobalDataBase::GetIntParaFromDB("statisticMethod");
    }

    CFLStart                  = GlobalDataBase::GetDoubleParaFromDB("CFLStart");
    CFLEnd                    = GlobalDataBase::GetDoubleParaFromDB("CFLEnd");

    LUSGSTolerance            = GlobalDataBase::GetDoubleParaFromDB("LUSGSTolerance");

    refMachNumber             = GlobalDataBase::GetDoubleParaFromDB("refMachNumber");
    refReNumber               = GlobalDataBase::GetDoubleParaFromDB("refReNumber");
    oRefReNumber              = GlobalDataBase::GetDoubleParaFromDB("oreynolds");

    refGama                   = GlobalDataBase::GetDoubleParaFromDB("refGama");

    TorqueRefX                = GlobalDataBase::GetDoubleParaFromDB("TorqueRefX");
    TorqueRefY                = GlobalDataBase::GetDoubleParaFromDB("TorqueRefY");
    TorqueRefZ                = GlobalDataBase::GetDoubleParaFromDB("TorqueRefZ");

    refDimensionalPressure    = GlobalDataBase::GetDoubleParaFromDB("refDimensionalPressure");
    refDimensionalDensity     = GlobalDataBase::GetDoubleParaFromDB("refDimensionalDensity");
    refDimensionalVelocity    = GlobalDataBase::GetDoubleParaFromDB("refDimensionalVelocity");
    refDimensionalSonicSpeed  = GlobalDataBase::GetDoubleParaFromDB("refDimensionalSonicSpeed");
    refDimensionalTemperature = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");

    GlobalDataBase::GetData("lowerPlotFieldBox", &lowerPlotFieldBox, PHDOUBLE, 3);
    GlobalDataBase::GetData("upperPlotFieldBox", &upperPlotFieldBox, PHDOUBLE, 3);

    isOverset = GlobalDataBase::GetIntParaFromDB("codeOfOversetGrid");

    pMax = 0.5, pMin = 0.2, deltaMax = 0.5, magnifyFactor = 2.0;
    if (GlobalDataBase::IsExist("pMaxForCFL", PHDOUBLE, 1))
    {
        pMax = GlobalDataBase::GetDoubleParaFromDB("pMaxForCFL");
    }
    if (GlobalDataBase::IsExist("pMinForCFL", PHDOUBLE, 1))
    {
        pMin = GlobalDataBase::GetDoubleParaFromDB("pMinForCFL");
    }
    if (GlobalDataBase::IsExist("deltaMaxForCFL", PHDOUBLE, 1))
    {
        deltaMax = GlobalDataBase::GetDoubleParaFromDB("deltaMaxForCFL");
    }
    if (GlobalDataBase::IsExist("magnifyFactorForCFL", PHDOUBLE, 1))
    {
        magnifyFactor = GlobalDataBase::GetDoubleParaFromDB("magnifyFactorForCFL");
    }
    reduceFactor = 1.0 / magnifyFactor;
    if (GlobalDataBase::IsExist("reduceFactorForCFL", PHDOUBLE, 1))
    {
        reduceFactor = GlobalDataBase::GetDoubleParaFromDB("reduceFactorForCFL"); //reduceFactorForCFL
    }
    fineCFL = 0.5;
    if (GlobalDataBase::IsExist("fineCFL", PHDOUBLE, 1))
    {
        fineCFL = GlobalDataBase::GetDoubleParaFromDB("fineCFL"); //fineCFL
    }
    AusmpwPlusLimiter = 1.0;
    if (GlobalDataBase::IsExist("AusmpwPlusLimiter", PHDOUBLE, 1))
    {
        AusmpwPlusLimiter = GlobalDataBase::GetDoubleParaFromDB("AusmpwPlusLimiter"); //reduceFactorForCFL
    }


}

}