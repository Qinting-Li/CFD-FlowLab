#include "TurbSolver.h"
#include "Constants.h"
#include "TK_Exit.h"
#include "TK_Log.h"

namespace PHSPACE
{

void Turb_MxDQ(RDouble *mat,RDouble *dq,int numberOfTurbulentEquation, RDouble *df)
{
    for (int m = 0; m < numberOfTurbulentEquation; ++ m)
    {
        df[m] = mat[m] * dq[m];
    }
}

TurbSolver::TurbSolver()
{
}

TurbSolver::~TurbSolver()
{
}

string TurbSolver::CastName(const string &name)
{
    if (name == "res")
    {
        return "res_turb";
    }
    else if (name == "q")
    {
        return "q_turb";
    }

    return name;
}

int TurbSolver::GetNumberOfEquations()
{
    Param_TurbSolver *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();
    return nTurbulentEquation;
}

void TurbSolver::InitDependentVariables()
{
    Grid *fineGrid = GetGrid(0);
    //! Don't be used for STRUCTGRID because of turb viscosity being stored.
    if (fineGrid->Type() != STRUCTGRID)
    {
        ComputeViscousCoeff(fineGrid);
    }

    ComputeQTurbNodeValue(fineGrid);
}

void TurbSolver::InitCoarseGridsFlow()
{
}

void TurbSolver::Action(ActionKey *actkey)
{
    CFDSolver::Action(actkey);
}

void TurbSolver::FillActionKey(ActionKey *actkey, int action, int level)
{
    actkey->action   = action;
    actkey->solver   = TURBULENCE;
    actkey->solverID = this->GetIndex();
    actkey->kind     = SOLVER_BASED;
    actkey->level    = level;
}

void TurbSolver::Solve()
{
    int outnstep = 0;
    //! Get iteration and zonal number.
    GlobalDataBase::GetData("outnstep", &outnstep, PHINT, 1);

    Grid *grid = GetGrid(0);

    //! Zero the residuals of the finest grid.
    ZeroResiduals(grid);

    //! Solve NS for all grids starting from the finest grid.
    SolveOnGrid(grid);

    CommunicationInterfaceData();

    InterpointInformation *interpointInformation = grid->GetInterpointInfo();
    if (interpointInformation)
    {
        CommunicationInterpointData();
    }
    GlobalDataBase::UpdateData("outnstep", &outnstep, PHINT, 1);

    int intervalStepFlow = 100;
    GlobalDataBase::GetData("intervalStepFlow", &intervalStepFlow, PHINT, 1);

    if (outnstep % intervalStepFlow == 0)
    {
        //! To dump restart data for continual simulation.
        ActionKey *actkeyDumpRestartData = new ActionKey();
        FillActionKey(actkeyDumpRestartData, DUMP_RESTART, 0);
        DumpRestartData(actkeyDumpRestartData);
        delete actkeyDumpRestartData;
    }

    int intervalStepRes;
    GlobalDataBase::GetData("intervalStepRes", &intervalStepRes, PHINT, 1);
    if (outnstep % intervalStepRes == 0)
    {
        //! To dump residual data.
        ActionKey *actkeyDumpResidual = new ActionKey();
        FillActionKey(actkeyDumpResidual, DUMP_RESIDUAL, 0);
        DumpResidual(actkeyDumpResidual);
        delete actkeyDumpResidual;
    }
}

void TurbSolver::PostSolve(Grid *grid_in, int stage, int level)
{
    CommunicationInterfaceData();
}

void TurbSolver::DumpResultFile(Grid *grid, int level)
{
    int outnstep = GlobalDataBase::GetIntParaFromDB("outnstep");
    int intervalStepFlow = GlobalDataBase::GetIntParaFromDB("intervalStepFlow");

    bool isSubIterationDump = true; //! It is true default in steady simulations.

    Param_TurbSolver *parameters = GetControlParameters();
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

    int intervalStepRes = GlobalDataBase::GetIntParaFromDB("intervalStepRes");
    if (outnstep == 1 || outnstep % intervalStepRes == 0 && (!isUnsteady || !isSubIterationDump))
    {
        //! To dump residual data.
        ActionKey * actkeyDumpResidual = new ActionKey();
        FillActionKey(actkeyDumpResidual, DUMP_RESIDUAL, level);
        DumpResidual(actkeyDumpResidual);
        delete actkeyDumpResidual;
    }
}

void TurbSolver::Post()
{
    //! To dump restart data for continual simulation.
    ActionKey *actkeyDumpRestartData = new ActionKey();
    FillActionKey(actkeyDumpRestartData, DUMP_RESTART, 0);
    DumpRestartData(actkeyDumpRestartData);
    delete actkeyDumpRestartData;
}

void TurbSolver::RegisterCFDSolverInterfaceField()
{
    Param_TurbSolver *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();
    
    CFDSolver::RegisterInterfaceField("turb::q", PHDOUBLE, nTurbulentEquation);

    CFDSolver::RegisterInterfaceField("turb::dq", PHDOUBLE, nTurbulentEquation);

    CFDSolver::RegisterInterfaceField("turb::dqdx", PHDOUBLE, nTurbulentEquation);
    CFDSolver::RegisterInterfaceField("turb::dqdy", PHDOUBLE, nTurbulentEquation);
    CFDSolver::RegisterInterfaceField("turb::dqdz", PHDOUBLE, nTurbulentEquation);

    CFDSolver::RegisterInterfaceField("vist", PHDOUBLE, 1);
    CFDSolver::RegisterInterfaceField("blend", PHDOUBLE, 1);
}

void TurbSolver::RegisterCFDSolverInterpointField()
{
    Param_TurbSolver *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    CFDSolver::RegisterInterpointField("qTurbNode", PHDOUBLE, nTurbulentEquation);
}

void TurbSolver::RegisterCFDSolverOversetField()
{
    Param_TurbSolver *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

	CFDSolver::RegisterOversetField( "turb::q", PHDOUBLE, nTurbulentEquation );

	CFDSolver::RegisterOversetField( "turb::dqdx", PHDOUBLE, nTurbulentEquation );
	CFDSolver::RegisterOversetField( "turb::dqdy", PHDOUBLE, nTurbulentEquation );
	CFDSolver::RegisterOversetField( "turb::dqdz", PHDOUBLE, nTurbulentEquation );
 
	CFDSolver::RegisterOversetField( "turb::limit", PHDOUBLE, 1 );

	CFDSolver::RegisterOversetField( "turb::dq", PHDOUBLE, nTurbulentEquation );
}

//! Drive functions to calculate residuals for different orders and dimensions.
void TurbSolver::UpdateResiduals(Grid *grid)
{
    //! Set values for flow variables in ghost cells.
    Boundary(grid);

    //! Modify the variable value in the cell center with wall boundary face.
    //ResetWallScalar(grid);

    //! Zero the other variables of the grid.
    InitSpectrum(grid);

    //! Compute the source term.
    SourceFlux(grid);

    //! Compute the spectrum radius.
    //Spectrum(grid);

    //! Compute the inviscid flux.
    InviscidFlux(grid);

    //! Compute the viscous flux.
    ViscousFlux(grid);

    DualTimeSource(grid);

    ZeroResidualOfSpecialCells( grid );       //! Added By Guo Yongheng 20160825.

    ModifyResiduals(grid);

    FreeGradientProxy(grid);
}

void TurbSolver::SourceFlux(Grid *grid)
{
    Param_TurbSolver *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();
    if (nTurbulentEquation == 1)
    {
        //! One equation turbulence model.
        SourceFluxOneEquation(grid);
    }
    else if (nTurbulentEquation == 2)
    {
		//! Two equation turbulence model.
        SourceFluxTwoEquation(grid);
    }
}

//! Compare the outnstep in the flowfield file between the NS and turbulence.
void TurbSolver::CompareOutStepOfFlowfieldFile(int outnstepofNS, int outnstepofTurb) const
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

		GlobalDataBase::UpdateData( "outnstep", &outnstepofTurb, PHINT, 1 );
    }
    else
    {
        TK_Exit::ExceptionExit("Error: Cluster writes files error! Please restart the calculation. \n");
    }
}

//! Get the file name for Residual dumping.
const string TurbSolver::GetResidualFileName()
{
    string turbresfile = "turbres.dat";
    GlobalDataBase::GetData("turbresfile", &turbresfile, PHSTRING, 1);
    return turbresfile;
}

//! Get the file name for Residual dumping.
const string TurbSolver::GetRestartFileName()
{
    string turbfile = "turbres.dat";
    GlobalDataBase::GetData("turbfile", &turbfile, PHSTRING, 1);
    return turbfile;
}

LIB_EXPORT Param_TurbSolver *TurbSolver::GetControlParameters()
{
    return static_cast <Param_TurbSolver *> (controlParameters);
}

}
