#include "PHMpi.h"
#include "LESSolver.h"
#include "Constants.h"
#include "Precision.h"
#include "GlobalDataBase.h"

namespace PHSPACE
{

bool LESSolver::init_param = false;
bool LESSolver::free_param = false;
int  LESSolver::nTurbEquations = 0;

LESSolver::LESSolver()
{
    //InitParameter();
}

LESSolver::~LESSolver()
{
    FreeParameter();
}

int LESSolver::GetNumberOfEquations()
{
    int nTurbulence = GlobalDataBase::GetIntParaFromDB("n_turb");
    return nTurbulence;
}

void LESSolver::FreeParameter()
{
    if (free_param) return;
    free_param = true;

    RDouble * turboo = reinterpret_cast<RDouble *> (GlobalDataBase::GetDataPtr("turboo"));
    delete [] turboo;
}

void LESSolver::InitParameter()
{
    /*
    string sgsmodel;
    GlobalDataBase::GetData("sgsmodel", &sgsmodel, PHSTRING, 1);
    int isgsmodel = 0;
    int nwtf = 0;
    RDouble coef_testfilter[3];

    if (sgsmodel.substr(0, 11) == "smagorinsky")
    {
        isgsmodel = 1;
    }
    else if (sgsmodel.substr(0, 18) == "dynamicsmagorinsky")
    {
        isgsmodel = 2;
        coef_testfilter[0] = 0.16666667;
        coef_testfilter[1] = 0.66666667;
        coef_testfilter[2] = 0.16666667;
        nwtf = 1;
        GlobalDataBase::UpdateData("coef_testfilter" ,&coef_testfilter , PHDOUBLE, 3);
        GlobalDataBase::UpdateData("nwtf" ,&nwtf , PHINT, 1);
    }
    else if (sgsmodel.substr(0, 4) == "wale")
    {
        isgsmodel = 3;
    }

    GlobalDataBase::UpdateData("isgsmodel", &isgsmodel , PHINT, 1);

    return;
    */
}

void LESSolver::Action(ActionKey *actkey)
{
    CFDSolver::Action(actkey);
}

void LESSolver::FillActionKey(ActionKey *actkey, int action, int level)
{
    actkey->action   = action;
    actkey->solver   = TURBULENCE;
    actkey->solverID = this->GetIndex();
    actkey->kind     = SOLVER_BASED;
    actkey->level    = level;
}

int LESSolver::GetNPostSolve()
{
    return 3;
}

void LESSolver::PostSolve(Grid *grid_in, int stage, int level)
{
    return;
}

void LESSolver::Post()
{
    return;
}

LIB_EXPORT Param_LESSolver *LESSolver::GetControlParameters()
{
    return static_cast <Param_LESSolver *> (controlParameters);
}

}