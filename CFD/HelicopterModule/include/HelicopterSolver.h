#pragma once
using namespace std;

namespace PHSPACE
{
    class Region;
    class MYHelicopter
    {
    private:
        Region * region;
        int numberOfZones, numberOfSolvers, numberOfLevels;
        int numberOfExternalIterations, intervalOfPlotFlow;
        bool codeOfTurbulenceFlow, isOverset, codeOfUnsteadyFlow, codeOfStaticDes;
    public:
        MYHelicopter(Region * region);
        ~MYHelicopter();
    public:
        void ReadGeometry();
        void Run();
    private:
        void InitEvolutionModel();
        void DoTimeEvolution();
        void UpdateIterationStep(int & outnstep);
        void SolveSteadyFlowField();
        void SolveUnsteadyFlowField();
        bool IsLocallyConvergent(const int & innstep, const int & min_sub_iter, const int & max_sub_iter, const double & tol_sub_iter, const double & valueOfErrorNorm);
        void PlotFlowField();
        void OutputNumericalField();
        void RunGlobalTaskList();
    };

}