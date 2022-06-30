#include "PHMpi.h"
#include "LESSolver.h"
#include "Constants.h"
#include "Precision.h"
#include "GlobalDataBase.h"
#include "TK_Log.h"
#include "TK_Exit.h"

namespace PHSPACE
{

bool LESSolver::init_param = false;
bool LESSolver::free_param = false;
int  LESSolver::nTurbEquations = 0;

LESSolver::LESSolver()
{
}

LESSolver::~LESSolver()
{
    FreeParameter();
}

void LESSolver::InitDependentVariables()
{
    Grid *fineGrid = GetGrid(0);
    //! Don't be used for STRUCTGRID because of turb viscosity being stored.
    if (fineGrid->Type() != STRUCTGRID)
    {
        ComputeGradient(fineGrid);

        ObtainViscosity(fineGrid);
    }

    ComputeQTurbNodeValue(fineGrid);
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

void LESSolver::PostSolve(Grid *grid_in, int stage, int level)
{
    CommunicationInterfaceData();
}

void LESSolver::DumpResultFile(Grid *grid, int level)
{
    int outnstep = GlobalDataBase::GetIntParaFromDB("outnstep");
    int intervalStepFlow = GlobalDataBase::GetIntParaFromDB("intervalStepFlow");

    bool isSubIterationDump = true; //! It is true default in steady simulations.

    Param_LESSolver *parameters = GetControlParameters();
    int isUnsteady = parameters->GetIsUnsteady();
    if (isUnsteady)
    {
        GlobalDataBase::GetData("isSubIterationDump", &isSubIterationDump, PHBOOL, 1);
    }

    if (outnstep % intervalStepFlow == 0 && level == 0 && isSubIterationDump)
    {
        //! To dump restart data for continual simulation.
        ActionKey *actkeyDumpRestartData = new ActionKey();
        FillActionKey(actkeyDumpRestartData, DUMP_RESTART, 0);
        DumpRestartData(actkeyDumpRestartData);
        delete actkeyDumpRestartData;
    }
}

void LESSolver::Post()
{
    //! To dump restart data for continual simulation.
    ActionKey *actkeyDumpRestartData = new ActionKey();
    FillActionKey(actkeyDumpRestartData, DUMP_RESTART, 0);
    DumpRestartData(actkeyDumpRestartData);
    delete actkeyDumpRestartData;
}

void LESSolver::RegisterCFDSolverInterfaceField()
{
    CFDSolver::RegisterInterfaceField("vist", PHDOUBLE, 1);

    CFDSolver::RegisterInterfaceField("subgridScaleEnergy", PHDOUBLE, 1);

    CFDSolver::RegisterInterfaceField("turbulentPrandtlNumber", PHDOUBLE, 1);
}

//! Compare the outnstep in the flowfield file between the NS and turbulence.
void LESSolver::CompareOutStepOfFlowfieldFile(int outnstepofNS, int outnstepofTurb) const
{
    if (outnstepofNS == outnstepofTurb) return;

    ostringstream oss;
    oss << "outnstepofNS   = " << outnstepofNS   << endl;
    oss << "outnstepofTurb = " << outnstepofTurb << endl;
    oss << "The NS flowfield file is not matched with the turbulence." << endl;

    int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
    if (!isUnsteady)
    {
        oss << "Warning: Cluster writes files error! It needs more iteration step." << endl;
        PrintToWindow(oss.str());
        WriteLogFile (oss.str());

		GlobalDataBase::UpdateData("outnstep", &outnstepofTurb, PHINT, 1);
    }
    else
    {
        TK_Exit::ExceptionExit("Error: Cluster writes files error! Please restart the calculation. \n");
    }
}

//! Get the file name for Residual dumping.
const string LESSolver::GetRestartFileName()
{
    string turbfile = "turbres.dat";
    GlobalDataBase::GetData("turbfile", &turbfile, PHSTRING, 1);
    return turbfile;
}

LIB_EXPORT Param_LESSolver *LESSolver::GetControlParameters()
{
    return static_cast <Param_LESSolver *> (controlParameters);
}

}