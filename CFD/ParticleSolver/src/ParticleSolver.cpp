#include "ParticleSolver.h"
#include "LIB_Macro.h"
#include "Param_ParticleSolver.h"
using namespace std;

namespace PHSPACE
{

ParticleSolver::ParticleSolver()
{

}

ParticleSolver::~ParticleSolver()
{
    //! Destructors can cause many problems of deletion by mistake.
    //! We must be careful. You can comment out this sentence 
    //! after the program debugging is perfect.
    DeAllocateGlobalVariables();
    FreeControlParameters();
}

LIB_EXPORT Param_ParticleSolver *ParticleSolver::GetControlParameters() const
{
    return static_cast <Param_ParticleSolver*> (controlParameters);
}

void ParticleSolver::InitMemory()
{

}

void ParticleSolver::InitFlow()
{
    Param_CFDSolver *parameters = GetControlParameters();
}

void ParticleSolver::CommunicationInterfaceData()
{

}

void ParticleSolver::CommunicationInterpointData()
{

}

void ParticleSolver::Action(ActionKey *actkey)
{

}

void ParticleSolver::TranslateAction(ActionKey *actkey)
{

}

streamsize ParticleSolver::TranslateActionLength(ActionKey *actkey)
{
    return 0;
}

void ParticleSolver::PrepareInterfaceData(Grid *grid, Data_ParamFieldSuite *datastore, 
                                          InterfaceDataProxy *interfaceDataProxy)
{

}

void ParticleSolver::Solve()
{

}

void ParticleSolver::ComputePostVisualVariables(Post_Visual *postVisualization)
{

}

void ParticleSolver::InitCoarseGridsFlow()
{

}

void ParticleSolver::ActionReflect(ActionTag *acttag)
{

}

void ParticleSolver::ComputePostProbesVariables(Post_Probes *postProbesVar)
{

}

int ParticleSolver::GetNumberOfEquations()
{
    return 0;
}

void ParticleSolver::Post()
{

}

void ParticleSolver::DumpResultFile(Grid *grid, int level)
{

}

void ParticleSolver::RungeKuttaResidual(Grid *grid, FieldProxy *dqProxy, RDouble coef)
{

}

void ParticleSolver::ZeroResiduals(Grid *grid)
{

}

}