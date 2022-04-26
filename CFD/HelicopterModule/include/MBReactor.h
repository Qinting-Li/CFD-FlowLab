#pragma once
#include "Region.h"
#include "CFDSolver.h"
using namespace std;

namespace PHSPACE
{
class ReginalFieldProxy
{
private:
    FieldProxy ** fieldProxy;
    int numberOfZones;
public:
    ReginalFieldProxy(int nZones);
    ~ReginalFieldProxy();
public:
    FieldProxy ** GetFieldProxy() { return fieldProxy; }
};

class MultiFieldProxy
{
private:
    int numberOfLevels;
    ReginalFieldProxy ** reginalFieldProxy;
public:
    MultiFieldProxy(int numberOfLevels, int nZones);
    ~MultiFieldProxy();
public:
    FieldProxy ** GetFieldProxy(int level) { return reginalFieldProxy[level]->GetFieldProxy(); }
};

class MYReactor
{
private:
    Region * region;
    MultiFieldProxy * q_proxy, * rhs_proxy;
    int solverIndex, timeMethod;
    int numberOfZones, numberOfLevels, numberOfEquations;
    int numberOfSmoothSteps, numberOfFilterSteps, numberOfFasSteps, numberOfSubSteps;
    bool codeOfPrecondition;
public:
    MYReactor(Region * region, int solverIndex);
    ~MYReactor();
public:
    void InitFlowField(int level);
    void SolveSteadyFlowField();
    void SolveUnsteadyFlowField();
    RDouble ComputeErrorNorm(int level);
    void PostCurrentStatus(int level);
public:
    bool IsCoarsestLevel(int level)      { return (level == numberOfLevels - 1); }
    FieldProxy ** GetQProxy(int level)   { return q_proxy->GetFieldProxy(level);   }
    FieldProxy ** GetRHSProxy(int level) { return rhs_proxy->GetFieldProxy(level); }
    Zone * GetZone(int izone)            { return region->GetZone(izone); }
    CFDSolver * GetCFDSolver(int izone)  { return static_cast< CFDSolver * >(GlobalSolvers::GetSolver(izone,this->solverIndex)); }
private:
    void SolveLinearSystem(int level, int numberOfSteps);
    void UpdateSteadyFlowField(int level, FieldProxy ** dq_proxy);
    void UpdateUnsteadyFlowField(int level);
    void CalculateTimeStep(int level);
    RDouble CalculateMinimumTimeStep(int level);
    void ReduceMaximumTimeStep(int level);
    void ZeroResidualOfModel(int level);
    void IntegrateResidualValue(int level);
    void RunGlobalTaskList();
    void UpdateResidualValue(int level);
    void TradeContactFlowField(int level);
    void RunModernLusgsModel(int level);
    void InitResidualByNegativeRhs(int level);
    void InitRhsByPositiveResidual(int level);
    void RestrictCoarseFlowField(int level);
    void RestrictDefectOfModel(int level);
    void PutNegativeCorrection(int level);
    void PutPositiveCorrection(int level);
    void CorrectFlowFieldOnFineGrid(int level);
    void RecoverResidualOnFineGrid(int level);
    int  GetPostNumberOfSolver();
    bool IsLocallyConvergent(int level);
    void LoadFlowField(int level);
    void ScanEquationForwardly(int level, FieldProxy ** dq_proxy);
    void ScanEquationBackwardly(int level, FieldProxy ** dq_proxy);
    void RunMultiGridSolver(int level);
    int  FindZoneSeriesOnNeighbor(int jZone, int iZone);
    void TradeContactDeltaField(int level, FieldProxy ** dq_proxy);
};

}