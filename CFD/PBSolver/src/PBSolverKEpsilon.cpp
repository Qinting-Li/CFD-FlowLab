#include "PBSolverKEpsilon.h"

namespace PHSPACE
{
PBSolverKEpsilon::PBSolverKEpsilon(const string &name):PBSolverScalar(name)
{

}

PBSolverKEpsilon::~PBSolverKEpsilon()
{

}

void PBSolverKEpsilon::InitMemory()
{
    geometry = GetGeometry();
    vector<PHSolver *> solverList;
    PHSolver *solver = new PBSolverTurbK("TURB_K");
    solverList.push_back(solver);
    solver->SetGeometry(geometry);
    solverKEqn = static_cast<CFDSolver *>(solver);
    solver = new PBSolverKEpsilon("TURB_EPSILON");
    solverList.push_back(solver);
    solver->SetGeometry(geometry);
    solverEpsilonEqn = static_cast<CFDSolver *>(solver);

    solverKEqn->InitMemory();
    solverEpsilonEqn->InitMemory();
}

void PBSolverKEpsilon::ReinitField()
{
    solverKEqn->ReinitField();
    solverEpsilonEqn->ReinitField();
}

void PBSolverKEpsilon::UpdateUnsteadyVariable(Grid *grid)
{
    solverKEqn->UpdateUnsteadyVariable(grid);
    solverEpsilonEqn->UpdateUnsteadyVariable(grid);
}

void PBSolverKEpsilon::SolveKEquation(Region &region)
{
    solverKEqn->Discrt();

    solverKEqn->SolveScalarEquation(region);

    solverKEqn->Update();

    solverKEqn->UpdateBoundary();

    solverKEqn->UpdateProperties();
}

void PBSolverKEpsilon::SolveEpsilonEquation(Region &region)
{
    solverEpsilonEqn->Discrt();

    solverEpsilonEqn->SolveScalarEquation(region);

    solverEpsilonEqn->Update();

    solverEpsilonEqn->UpdateBoundary();

    solverEpsilonEqn->UpdateProperties();
}

void PBSolverKEpsilon::UpdateAfterIterloop()
{
    solverKEqn->UpdateAfterIterloop();

    solverEpsilonEqn->UpdateAfterIterloop();
}

void PBSolverKEpsilon::ComputePostVisualVariables(Post_Visual *postVisualization)
{
    solverKEqn->ComputePostVisualVariables(postVisualization);

    solverEpsilonEqn->ComputePostVisualVariables(postVisualization);
}

void PBSolverKEpsilon::InitIncompressibleResultFile(Grid *grid, string fileName)
{
    solverKEqn->InitIncompressibleResultFile(grid, fileName);
    solverEpsilonEqn->InitIncompressibleResultFile(grid, fileName);
}

void PBSolverKEpsilon::DumpIncompressibleResultFile(Grid *grid, string fileName)
{
    solverKEqn->DumpIncompressibleResultFile(grid, fileName);
    solverEpsilonEqn->DumpIncompressibleResultFile(grid, fileName);
}

void PBSolverKEpsilon::SolutionIsConverged(int iZone, bool *flag)
{
    bool flagK = false;
    bool flagEpsilon = false;

    solverKEqn->SolutionIsConverged(iZone, &flagK);
    solverEpsilonEqn->SolutionIsConverged(iZone, &flagEpsilon);

    if (flagK && flagEpsilon)
    {
        *flag = true;
    }
    else
    {
        *flag = false;
    }
}

DEFINE_PBSOLVER_REGISTER_FUNCTION(PBSolver, PBSolverKEpsilon, SOLVER_TURB_KEPSILON);
}