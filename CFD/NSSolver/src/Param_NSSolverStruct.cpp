#include "Param_NSSolverStruct.h"
#include "GlobalDataBase.h"

namespace PHSPACE
{

LIB_EXPORT Param_NSSolverStruct::Param_NSSolverStruct()
{

}

LIB_EXPORT Param_NSSolverStruct::~Param_NSSolverStruct()
{

}

LIB_EXPORT void Param_NSSolverStruct::Init()
{
    Param_NSSolver::Init();
    MUSCLCoefXk = GlobalDataBase::GetDoubleParaFromDB("MUSCLCoefXk");
    MUSCLCoefXb = GlobalDataBase::GetDoubleParaFromDB("MUSCLCoefXb");
    catalyticCoef = GlobalDataBase::GetDoubleParaFromDB("catalyticCoef");
    inviscidSchemeName = GlobalDataBase::GetStrParaFromDB("inviscidSchemeName");
}

}