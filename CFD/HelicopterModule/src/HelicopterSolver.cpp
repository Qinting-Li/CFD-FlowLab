#include "HelicopterSolver.h"
#include "MBReactor.h"
#include "Controller.h"
#include "AleModel.h"
#include "Post_WriteTecplot.h"
#include "NSSolverUnstruct.h"
#include "TurbSolverUnstr.h"
#include "Zone.h"
#include "Task.h"
using namespace std;

namespace PHSPACE
{
MYHelicopter::MYHelicopter(Region * region)
{
    this->region = region;

    this->codeOfTurbulenceFlow = (GlobalDataBase::GetIntParaFromDB("viscousType") > LAMINAR);
    this->isOverset            = (GlobalDataBase::GetIntParaFromDB("codeOfOversetGrid") != 0);
    this->codeOfUnsteadyFlow   = (GlobalDataBase::GetIntParaFromDB("iunsteady") != 0);
    this->codeOfStaticDes      = (GlobalDataBase::GetIntParaFromDB("ifStaticsFlowField") != 0);

    this->numberOfExternalIterations = GlobalDataBase::GetIntParaFromDB("maxSimuStep");
    this->intervalOfPlotFlow         = GlobalDataBase::GetIntParaFromDB("intervalStepPlot");
}

MYHelicopter::~MYHelicopter()
{

}

void MYHelicopter::Run()
{
    InitEvolutionModel();

    DoTimeEvolution();

    OutputNumericalField();
}

void MYHelicopter::ReadGeometry()
{
    region->ReadGrid();
    region->UpdateAllZonesPara();
    region->InitGeometry();
    region->SwapGeommetry();
}

void MYHelicopter::InitEvolutionModel()
{
    this->numberOfZones = PHMPI::GetNumberofGlobalZones();
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        Zone * zone = region->GetZone(iZone);
        if (! zone) continue;

        vector<PHSolver *> solverList;
        PHSolver *solver = new NSSolverUnstruct();
        solverList.push_back(solver);
        
        if (codeOfTurbulenceFlow)
        {
            PHSolver *solver = new TurbSolverUnstr();
            solverList.push_back(solver);
        }

        //! Finally, push the solvers into the list.
        PHGeometry *geometry = zone->GetGeometry();
        for (std::size_t iSolver = 0; iSolver < solverList.size(); ++ iSolver)
        {
            PHSolver *solver = solverList[iSolver];

            solver->SetGeometry(geometry);
            solver->SetIndex(static_cast<int>(iSolver));
            zone->AddSolver(solver);
        }
    }

    this->numberOfSolvers = Zone::GetNumberOfSolvers();
    for (int iSolver = 0; iSolver < numberOfSolvers; ++ iSolver)
    {
        MYReactor * reactor = new MYReactor(this->region, iSolver);

        reactor->InitFlowField(0);

        delete reactor; reactor = NULL;
    }
    region->ReSetBoundaryConditionByGlobalBC();
    return;
}

void MYHelicopter::UpdateIterationStep(int & outnstep)
{
    outnstep += 1;
    GlobalDataBase::UpdateData("outnstep", & outnstep, PHINT, 1);

    int newnstep = GlobalDataBase::GetIntParaFromDB("newnstep");
    newnstep += 1;
    GlobalDataBase::UpdateData("newnstep", & newnstep, PHINT, 1);

    if (codeOfUnsteadyFlow && codeOfStaticDes)
    {
        int startStatisticStep = GlobalDataBase::GetIntParaFromDB("startStatisticStep");
        if(outnstep >= startStatisticStep)
        {
            int nStatisticalStep = GlobalDataBase::GetIntParaFromDB("nStatisticalStep");
            nStatisticalStep += 1;
            GlobalDataBase::UpdateData("nStatisticalStep", & nStatisticalStep, PHINT, 1);
        }
    }
    return;
}

void MYHelicopter::PlotFlowField()
{
    Post_WriteTecplot * plotter = new Post_WriteTecplot();

    plotter->Run();

    delete plotter; plotter = NULL;

    return;
}

void MYHelicopter::DoTimeEvolution()
{
    int outnstep = GlobalDataBase::GetIntParaFromDB("outnstep");

    bool active = true;
    while (active)
    {
        UpdateIterationStep(outnstep);

        AleModel::UnsteadyPreSolve();

        if (codeOfUnsteadyFlow)
        {
            SolveUnsteadyFlowField();
        }
        else
        {
            SolveSteadyFlowField();
        }

        if (outnstep % intervalOfPlotFlow == 0) PlotFlowField();

        if (outnstep >= numberOfExternalIterations) active = false;
    }
    return;
}

void MYHelicopter::SolveSteadyFlowField()
{
    for (int iSolver = 0; iSolver < numberOfSolvers; ++ iSolver)
    {
        MYReactor * reactor = new MYReactor(this->region, iSolver);

        reactor->SolveSteadyFlowField();

        reactor->PostCurrentStatus(0);

        delete reactor; reactor = NULL;
    }
    return;
}

bool MYHelicopter::IsLocallyConvergent(const int & innstep, const int & min_sub_iter, const int & max_sub_iter, const double & tol_sub_iter, const double & valueOfErrorNorm)
{
    bool convergent = false;
    if (innstep < min_sub_iter)
    {
        convergent = false;
    }
    else if (valueOfErrorNorm <= tol_sub_iter)
    {
        convergent = true;
    }
    else if (innstep >= max_sub_iter)
    {
        convergent = true;
    }
    else
    {
        convergent = false;
    }
    return convergent;
}

void MYHelicopter::SolveUnsteadyFlowField()
{
    int     min_sub_iter = GlobalDataBase::GetIntParaFromDB("min_sub_iter");
    int     max_sub_iter = GlobalDataBase::GetIntParaFromDB("max_sub_iter");
    RDouble tol_sub_iter = GlobalDataBase::GetDoubleParaFromDB("tol_sub_iter");

    int innstep = 0;
    int nSolver = Zone::GetNumberOfSolvers();

    bool isSubIterationDump = false;
    GlobalDataBase::UpdateData("isSubIterationDump", &isSubIterationDump, PHBOOL, 1);

    bool isSubIterationConvergence = false;
    while(!isSubIterationConvergence)
    {
        innstep += 1;
        GlobalDataBase::UpdateData("innstep", &innstep, PHINT, 1);

        RDouble subIterationNormOneSolver     =  0.0;
        RDouble maxSubIterationNormAllSolvers = -LARGE;

        for (int iSolver = 0; iSolver < nSolver; ++ iSolver)
        {
            Controller * controller = new Controller(iSolver, this->region);

            controller->SolveSteadyField();

            //! Here is in the sub-iteration, PostSolve is only for the CommunicationInterfaceData and dumps the "res.dat", which contains index of the sub-iteration residual.
            controller->PostSolve();

            subIterationNormOneSolver     = controller->ComputeSubIterationNorm(iSolver, 0);
            maxSubIterationNormAllSolvers = max(maxSubIterationNormAllSolvers, subIterationNormOneSolver);

            delete controller;
        }

        PH_CompareMaxMin(maxSubIterationNormAllSolvers, 1);

        //! Sub-iteration continue or not judgement.
        isSubIterationConvergence = IsLocallyConvergent(innstep, min_sub_iter, max_sub_iter, tol_sub_iter, maxSubIterationNormAllSolvers);
    }

    //Here is already (innstep == max_sub_iter || isSubIterationConvergence).
    isSubIterationDump = true;
    GlobalDataBase::UpdateData("isSubIterationDump", &isSubIterationDump, PHBOOL, 1);

    //! Update physical time step flowfield after sub-iteration convergence.
    for (int iSolver = 0; iSolver < nSolver; ++ iSolver)
    {
        Controller * controller = new Controller(iSolver, this->region);

        controller->UpdateUnsteadyFlow(iSolver, 0);

        //! Here is out of the sub-iteration, PostSolve is to dump various files according to the parameters in the "cfd_para.hypara", such as intervalStepFlow, intervalStepPlot, intervalStepForce, intervalStepRes, et al.
        controller->PostSolve();

        delete controller;
    }
}

void MYHelicopter::OutputNumericalField()
{
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        Zone * zone = region->GetZone(iZone);
        if (! zone) continue;
        for (int iSolver = 0; iSolver < numberOfSolvers; ++ iSolver)
        {
            PHSolver * solver = PHSPACE::GlobalSolvers::GetSolver(iZone, iSolver);
            solver->Post();
        }
    }
    this->RunGlobalTaskList();
    this->PlotFlowField();
    return;
}

void MYHelicopter::RunGlobalTaskList()
{
    int numberOfTasks = static_cast< int >(Task::ootasklist->size());
    for (int iTask = 0; iTask < numberOfTasks; ++ iTask)
    {
        (* Task::ootasklist)[iTask]->RunTask();
        delete (* Task::ootasklist)[iTask]; (* Task::ootasklist)[iTask] = NULL;
    }
    Task::ootasklist->clear();
    return;
}

}