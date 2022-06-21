#include <iostream>
#include <algorithm>
#include "PBSolver.h"
#include "TK_Exit.h"
#include "TK_Time.h"
#include "PBOutput.h"
#include "Post_WriteTecplot.h"
#include "Post_ProbesWrite.h"
#include "Task_ServerUpdateInterface.h"

namespace PHSPACE
{

PBSolver::PBSolver(const string &name):PBObject(name)
{

}

PBSolver::~PBSolver()
{

}

void PBSolver::AllocateGlobalVar(Grid *grid)
{
    int nTotalNode = grid->GetNTotalNode();

    if (grid->GetLevel() == 0)
    {
        int *nodeBCType = new int[nTotalNode];
        PHSPACE::SetField(nodeBCType, 0, nTotalNode);
        grid->UpdateDataPtr("nodeBCType", nodeBCType);
    }
}

void PBSolver::DeAllocateGlobalVar(Grid *grid)
{
    PrintToWindow("PBSolver::DeAllocateGlobalVar\n");
}

void PBSolver::ZeroResiduals(Grid *grid)
{
    PrintToWindow("PBSolver::ZeroResiduals\n");
}

void PBSolver::TimeStep(Grid *grid)
{
    PrintToWindow("PBSolver::TimeStep\n");
}

void PBSolver::ComputeResidual(Grid *grid_in, FieldProxy *rhsProxy)
{
    PrintToWindow("PBSolver::ComputeResidual\n");
}

void PBSolver::SourceFlux(Grid *grid)
{
    PrintToWindow("PBSolver::SourceFlux\n");
}

void PBSolver::InviscidFlux(Grid *grid)
{
    PrintToWindow("PBSolver::InviscidFlux\n");
}

void PBSolver::ViscousFlux(Grid *grid)
{
    PrintToWindow("PBSolver::ViscousFlux\n");
}

void PBSolver::RungeKutta(Grid *grid, FieldProxy *rhsProxy)
{
    PrintToWindow("PBSolver::RungeKutta\n");
}

void PBSolver::Post()
{
    PrintToWindow("PBSolver::Post\n");
}

void PBSolver::PostSolve(Grid *grid, int stage, int level)
{
    PrintToWindow("PBSolver::PostSolve\n");
}

RDouble PBSolver::UnsteadyConvergence(Grid *grid)
{
    return 0;
}

void PBSolver::InitMemory()
{

}

void PBSolver::CommunicateInterfaceData(vector<string> &nameList)
{
    UploadInterfaceData(nameList);     
    UpdateInterfaceData();             //!Including data exchange,CFDSolver
    DownloadInterfaceData(nameList);
}

void PBSolver::UploadInterfaceData(vector<string> &nameList)
{
    int zoneIDofThisSolver = GetZoneLocalID();
    if (zoneIDofThisSolver != 0)
    {
        //! Very Important: only the first zone on a processor need to implement the task.
        //! other wise, the tasks would be implemented several times.
        return;
    }

    using namespace PHMPI;

    int nZones = GetNumberofGlobalZones();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int ZoneProcessorID = GetZoneProcessorIDSepMode(iZone);
        int CurrentProcessorID = GetCurrentProcessorID();

        if (CurrentProcessorID == ZoneProcessorID)
        {
            PBSolver *solver = static_cast <PBSolver *> (GlobalSolvers::GetSolver(iZone, GetIndex()));
            solver->UploadInterfaceValue(nameList);
        }
    }
}

void PBSolver::DownloadInterfaceData(vector<string> &nameList)
{
    int zoneIDofThisSolver = GetZoneLocalID();
    if (zoneIDofThisSolver != 0)
    {
        //! Very Important: only the first zone on a processor need to implement the task.
        //! other wise, the tasks would be implemented several times.
        return;
    }

    using namespace PHMPI;

    int nZones = GetNumberofGlobalZones();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int ZoneProcessorID = GetZoneProcessorIDSepMode(iZone);
        int CurrentProcessorID = GetCurrentProcessorID();

        if (CurrentProcessorID == ZoneProcessorID)
        {
            PBSolver *solver = static_cast <PBSolver *> (GlobalSolvers::GetSolver(iZone, GetIndex()));
            solver->DownloadInterfaceValue(nameList);
        }
    }
}

void PBSolver::UploadInterfaceValue(vector<string> &nameList)
{
    UnstructGrid *grid = UnstructGridCast(GetGrid());
    InterfaceInfo *iinfo = grid->GetInterfaceInfo();
    if (iinfo == 0) return;

    for (int i = 0; i < nameList.size(); ++ i)
    {
        RDouble *variable = reinterpret_cast <RDouble *> (grid->GetDataPtr(nameList[i]));
        PHSPACE::UploadInterfaceValue(grid, variable, nameList[i]);
    }
}

void PBSolver::DownloadInterfaceValue(vector<string> &nameList)
{
    UnstructGrid *grid = UnstructGridCast(GetGrid());
    InterfaceInfo *iinfo = grid->GetInterfaceInfo();
    if (iinfo == 0) return;

    for (int i = 0; i < nameList.size(); ++ i)
    {
        RDouble *variable = reinterpret_cast <RDouble *> (grid->GetDataPtr(nameList[i]));
        PHSPACE::DownloadInterfaceValue(grid, variable, nameList[i]);
    }
}

void PBSolver::UpdateInterfaceData()
{
    int zoneIDofThisSolver = GetZoneLocalID();
    if (zoneIDofThisSolver != 0)
    {
        //! Very Important: only the first zone on a processor need to implement the task.
        //! other wise, the tasks would be implemented several times.
        return;
    }

    using namespace PHMPI;
    int nLevel = GetNumberOfMultiGrid();
    for (int iLevel = 0; iLevel < nLevel; ++ iLevel)
    {
        ActionKey *actkey = new ActionKey();
        FillActionKey(actkey, UPDATE_INTERFACE_DATA, iLevel);

        //Task_ServerUpdateInterface *task = new Task_ServerUpdateInterface();
        if (IsNeedNonBlockingCommunication(actkey))
        {
            MainTaskNonBlocking(actkey);
        }
        else
        {
            MainTaskBlocking(actkey);
        }
        delete actkey;
    }
}

void PBSolver::CommunicationInterfaceData()
{
    if (PHMPI::IsParallelRun())
    {
        UploadOversetData();
        UpdateOversetData();
        DownloadOversetData();
    }
}

void PBSolver::FillActionKey(ActionKey *actkey, int action, int level)
{
    actkey->action = action;
    actkey->solver = NS_EQUATION;
    actkey->solverID = GetIndex();
    actkey->kind = SOLVER_BASED;
    actkey->level = level;
}
}