#include "MBReactor.h"
#include "MBConstTable.h"
#include "Task.h"
#include "TK_Log.h"
using namespace std;

namespace PHSPACE
{
    ReginalFieldProxy::ReginalFieldProxy(int numberOfZones)
    {
        this->numberOfZones = numberOfZones;
        fieldProxy = new FieldProxy * [numberOfZones];
        for (int iZone = 0; iZone < numberOfZones; ++ iZone)
        {
            fieldProxy[iZone] = NULL;
        }
    }

    ReginalFieldProxy::~ReginalFieldProxy()
    {
        for (int iZone = 0; iZone < numberOfZones; ++ iZone)
        {
            delete fieldProxy[iZone]; fieldProxy[iZone] = NULL;
        }
        delete [] fieldProxy; fieldProxy = NULL;
    }

    MultiFieldProxy::MultiFieldProxy(int numberOfLevels, int numberOfZones)
    {
        this->numberOfLevels = numberOfLevels;
        reginalFieldProxy = new ReginalFieldProxy * [numberOfLevels];
        for (int iLevel = 0; iLevel < numberOfLevels; ++ iLevel)
        {
            reginalFieldProxy[iLevel] = new ReginalFieldProxy(numberOfZones);
        }
    }

    MultiFieldProxy::~MultiFieldProxy()
    {
        for (int iLevel = 0; iLevel < numberOfLevels; ++ iLevel)
        {
            delete reginalFieldProxy[iLevel]; reginalFieldProxy[iLevel] = NULL;
        }
        delete [] reginalFieldProxy; reginalFieldProxy = NULL;
    }

    MYReactor::MYReactor(Region * region, int solverIndex)
    {
        this->region      = region;
        this->solverIndex = solverIndex;

        if (solverIndex == 0)
        {
            this->numberOfEquations = 5;
        }
        else
        {
            this->numberOfEquations = 0;
            int viscousType = GlobalDataBase::GetIntParaFromDB("viscousType");
            if (viscousType == ALGEBRAIC)
            {
                numberOfEquations = 0;
            }
            else if (viscousType == ONE_EQU)
            {
                numberOfEquations = 1;
            }
            else if (viscousType == TWO_EQU)
            {
                numberOfEquations = 2;
            }
            GlobalDataBase::UpdateData("n_turb", & numberOfEquations, PHINT, 1);
        }

        this->numberOfZones       = PHMPI::GetNumberofGlobalZones();
        this->numberOfLevels      = GlobalDataBase::GetIntParaFromDB("nMGLevel");
        this->numberOfSmoothSteps = GlobalDataBase::GetIntParaFromDB("MGPreIteration");
        this->numberOfFilterSteps = GlobalDataBase::GetIntParaFromDB("n_post");
        if (numberOfLevels > 1)
        {
            this->numberOfSubSteps = numberOfSmoothSteps + numberOfFilterSteps;
            this->q_proxy = new MultiFieldProxy(numberOfLevels, numberOfZones);
        }
        else
        {
            this->numberOfSubSteps = GlobalDataBase::GetIntParaFromDB("MGCoarsestIteration");
            this->q_proxy = NULL;
        }
        this->rhs_proxy = new MultiFieldProxy(numberOfLevels, numberOfZones);

        this->numberOfFasSteps = GlobalDataBase::GetIntParaFromDB("MGFasType");
        this->timeMethod       = GlobalDataBase::GetIntParaFromDB("ifLocalTimeStep");

        int ifLowSpeedPrecon = GlobalDataBase::GetIntParaFromDB("ifLowSpeedPrecon");
        this->codeOfPrecondition = (ifLowSpeedPrecon == 1 && this->solverIndex == 0);
    }

    MYReactor::~MYReactor()
    {
        delete q_proxy;   q_proxy   = NULL;
        delete rhs_proxy; rhs_proxy = NULL;
    }

    void MYReactor::InitFlowField(int level)
    {
        for (int iZone = 0; iZone < numberOfZones; ++ iZone)
        {
            Zone * zone = region->GetZone(iZone);
            if (! zone) continue;
            PHSolver * solver = static_cast< CFDSolver * >(PHSPACE::GlobalSolvers::GetSolver(iZone, this->solverIndex));
            solver->InitMemory();
        }

        for (int iZone = 0; iZone < numberOfZones; ++ iZone)
        {
            Zone * zone = region->GetZone(iZone);
            if (! zone) continue;
            PHSolver * solver = static_cast< CFDSolver * >(PHSPACE::GlobalSolvers::GetSolver(iZone, this->solverIndex));
            solver->InitFlow();
        }

        for (int iZone = 0; iZone < numberOfZones; ++ iZone)
        {
            Zone * zone = region->GetZone(iZone);
            if (! zone) continue;
            PHSolver * solver = static_cast< CFDSolver * >(PHSPACE::GlobalSolvers::GetSolver(iZone, this->solverIndex));
            solver->InitCoarseGridsFlow();
        }

        for (int iZone = 0; iZone < numberOfZones; ++ iZone)
        {
            Zone * zone = region->GetZone(iZone);
            if (! zone) continue;
            PHSolver * solver = static_cast< CFDSolver * >(PHSPACE::GlobalSolvers::GetSolver(iZone, this->solverIndex));
            solver->CommunicationInterfaceData();
        }
        RunGlobalTaskList();
        return;
    }

    void MYReactor::SolveSteadyFlowField()
    {
        ZeroResidualOfModel(0);

        RunMultiGridSolver(0);
    }

    void MYReactor::SolveUnsteadyFlowField()
    {
        int innstep = GlobalDataBase::GetIntParaFromDB("innstep");
        while(true)
        {
            innstep += 1;
            GlobalDataBase::UpdateData("innstep", & innstep, PHINT, 1);
            SolveSteadyFlowField();
            if (IsLocallyConvergent(0))
            {
                UpdateUnsteadyFlowField(0);
                break;
            }
            TradeContactFlowField(0);
        }
        return;
    }

    void MYReactor::UpdateUnsteadyFlowField(int level)
    {
        for (int iZone = 0; iZone < numberOfZones; ++ iZone)
        {
            Zone * zone = GetZone(iZone);
            if (! zone) continue;
            Grid *grid = PHSPACE::GetGrid(iZone, level);
            CFDSolver * solver = GetCFDSolver(iZone);
            solver->UpdateUnsteadyFlow(grid);
        }
        return;
    }

    void MYReactor::RunMultiGridSolver(int level)
    {
        InitRhsByPositiveResidual(level);

        if (IsCoarsestLevel(level) || this->solverIndex == 1)
        {
            SolveLinearSystem(level, numberOfSubSteps);
        }
        else
        {
            int coarseLevel = level + 1;
            SolveLinearSystem(level, numberOfSmoothSteps);
            RestrictCoarseFlowField(level);

            LoadFlowField(coarseLevel);
            InitResidualByNegativeRhs(level);
            UpdateResidualValue(level);
            RestrictDefectOfModel(level);
            UpdateResidualValue(coarseLevel);

            for (int iFasStep = 0; iFasStep < numberOfFasSteps; ++ iFasStep)
            {
                RunMultiGridSolver(coarseLevel);
            }
            PutNegativeCorrection(coarseLevel);
            CorrectFlowFieldOnFineGrid(level);

            SolveLinearSystem(level, numberOfFilterSteps);
            PutPositiveCorrection(coarseLevel);
        }
        RecoverResidualOnFineGrid(level);
        return;
    }

    void MYReactor::SolveLinearSystem(int level, int numberOfSteps)
    {
        int nLUSGSSweeps = GlobalDataBase::GetIntParaFromDB("nLUSGSSweeps");
        int turb_sweeps  = GlobalDataBase::GetIntParaFromDB("turb_sweeps");
        bool runSingleSweep = (this->solverIndex == 0 && nLUSGSSweeps == 1) || (this->solverIndex == 1 && turb_sweeps == 1);

        for (int iStep = 0; iStep < numberOfSteps; ++ iStep)
        {
            CalculateTimeStep(level);

            if (runSingleSweep || level != 0 || this->codeOfPrecondition)
            {
                IntegrateResidualValue(level);
                RunModernLusgsModel(level);
            }
            else
            {
            }
        }
        return;
    }

    void MYReactor::CalculateTimeStep(int level)
    {
        for (int iZone = 0; iZone < numberOfZones; ++ iZone)
        {
            Zone * zone = GetZone(iZone);
            if (! zone) continue;
            Grid *grid = PHSPACE::GetGrid(iZone, level);
            CFDSolver * solver = GetCFDSolver(iZone);
            solver->TimeStep(grid);
        }
        ReduceMaximumTimeStep(level);
        return;
    }

    RDouble MYReactor::CalculateMinimumTimeStep(int level)
    {
        using namespace PHMPI;

        int outnstep             = GlobalDataBase::GetIntParaFromDB("outnstep");
        int allReduceStep        = GlobalDataBase::GetIntParaFromDB("allReduceStep");
        RDouble globalMinTimeStep = GlobalDataBase::GetDoubleParaFromDB("globalMinTimeStep");
        if (outnstep % allReduceStep != 0 && globalMinTimeStep > 0.0)
        {
            return globalMinTimeStep;
        }

        RDouble globalMinDt = LARGE, localMaxDt = -LARGE;
        for (int iZone = 0; iZone < numberOfZones; ++ iZone)
        {
            Zone * zone = GetZone(iZone);
            if (! zone) continue;
            Grid *grid = PHSPACE::GetGrid(iZone, level);
            CFDSolver * solver = GetCFDSolver(iZone);
            RDouble minDt, maxDt;
            solver->ComputeMinTimeStep(grid, minDt, maxDt);
            globalMinDt = min(globalMinDt, minDt);
            localMaxDt  = max(localMaxDt, maxDt);
        }

        PH_CompareMaxMin(globalMinDt, 2);
        GlobalDataBase::UpdateData("globalMinTimeStep", & globalMinDt , PHDOUBLE, 1);

        int intervalStepForce = GlobalDataBase::GetIntParaFromDB("intervalStepForce");
        intervalStepForce = max(intervalStepForce, 200);
        if(outnstep == 1 || outnstep % intervalStepForce == 0)
        {
            RDouble globalMaxDt;
            int numberOfGridProcessor = PHMPI::GetNumberOfGridProcessor();
            if (numberOfGridProcessor == 0)
            {
                PH_Reduce(& localMaxDt, & globalMaxDt, 1, MPI_MAX);
            }
            else
            {
                PH_AllreduceSepMode(& localMaxDt, & globalMaxDt, 1, MPI_DOUBLE, MPI_MAX);
            }

            if (PHMPI::GetCurrentProcessorID() == PHMPI::GetServerProcessorID())
            {
                ostringstream oss;
                oss.setf(ios::scientific);
                oss.precision(2);
                oss << "    iter = " << outnstep 
                    << ", level = " << level
                    << ", max && min time step = " << globalMaxDt << ", " << globalMinDt
                    << ", max/min ratio = " << globalMaxDt / globalMinDt;
                WriteLogFile(oss);
            }
        }
        return globalMinDt;
    }

    void MYReactor::ReduceMaximumTimeStep(int level)
    {
        using namespace PHMPI;
        if (GlobalDataBase::GetIntParaFromDB("ifLocalTimeStep") != 0) return;
        if (this->solverIndex != 0) return;

        RDouble globalMinDt = CalculateMinimumTimeStep(level);
        for (int iZone = 0; iZone < numberOfZones; ++ iZone)
        {
            Zone * zone = GetZone(iZone);
            if (! zone) continue;
            Grid *grid = PHSPACE::GetGrid(iZone, level);
            CFDSolver * solver = GetCFDSolver(iZone);
            solver->ReduceMaxTimeStep(grid, globalMinDt);
        }
        return;
    }

    void MYReactor::ZeroResidualOfModel(int level)
    {
        for (int iZone = 0; iZone < numberOfZones; ++ iZone)
        {
            Zone * zone = GetZone(iZone);
            if (! zone) continue;
            Grid *grid = PHSPACE::GetGrid(iZone, level);
            CFDSolver * solver = GetCFDSolver(iZone);
            solver->ZeroResiduals(grid);
        }
        return;
    }

    void MYReactor::IntegrateResidualValue(int level)
    {
        FieldProxy ** rhs_proxy = GetRHSProxy(level);
        for (int iZone = 0; iZone < numberOfZones; ++ iZone)
        {
            Zone * zone = GetZone(iZone);
            if (! zone) continue;
            Grid *grid = PHSPACE::GetGrid(iZone, level);
            CFDSolver * solver = GetCFDSolver(iZone);
            solver->ComputeResidual(grid, rhs_proxy[iZone]);
        }
        return;
    }

    void MYReactor::UpdateResidualValue(int level)
    {
        for (int iZone = 0; iZone < numberOfZones; ++ iZone)
        {
            Zone * zone = GetZone(iZone);
            if (! zone) continue;
            Grid *grid = PHSPACE::GetGrid(iZone, level);
            CFDSolver * solver = GetCFDSolver(iZone);
            solver->UpdateResiduals(grid);
        }
        return;
    }

    void MYReactor::PostCurrentStatus(int level)
    {
        int npostsolve = GetPostNumberOfSolver();
        for (int n = 0; n < npostsolve; ++ n)
        {
            for (int iZone = 0; iZone < numberOfZones; ++ iZone)
            {
                Zone * zone = GetZone(iZone);
                if (! zone) continue;
                Grid *grid = PHSPACE::GetGrid(iZone, level);
                CFDSolver * solver = GetCFDSolver(iZone);
                solver->PostSolve(grid, n, level);
            }
            RunGlobalTaskList();
        }
        return;
    }

    RDouble MYReactor::ComputeErrorNorm(int level)
    {
        RDouble subIterationNorm    =  0.0;
        RDouble maxSubIterationNorm = -LARGE;
        for (int iZone = 0; iZone < numberOfZones; ++ iZone)
        {
            Zone * zone = GetZone(iZone);
            if (! zone) continue;
            Grid      * grid    = PHSPACE::GetGrid(iZone, level);
            CFDSolver * solver  = GetCFDSolver(iZone);
            subIterationNorm    = solver->UnsteadyConvergence(grid);
            maxSubIterationNorm = max(maxSubIterationNorm, subIterationNorm);
        }
        return maxSubIterationNorm;
    }

    void MYReactor::RunGlobalTaskList()
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

    int MYReactor::GetPostNumberOfSolver()
    {
        for (int iZone = 0; iZone < numberOfZones; ++ iZone)
        {
            Zone * zone = GetZone(iZone);
            if (! zone) continue;
            CFDSolver * solver = GetCFDSolver(iZone);
            return solver->GetNPostSolve();
        }
        return 1;
    }

    void MYReactor::RunModernLusgsModel(int level)
    {
        FieldProxy ** dq_proxy = new FieldProxy * [numberOfZones];
        for (int iZone = 0; iZone < numberOfZones; ++ iZone)
        {
            dq_proxy[iZone] = NULL;
        }

        ScanEquationForwardly(level, dq_proxy);

        TradeContactDeltaField(level, dq_proxy);

        ScanEquationBackwardly(level, dq_proxy);

        UpdateSteadyFlowField(level, dq_proxy);

        for (int iZone = 0; iZone < numberOfZones; ++ iZone)
        {
            delete dq_proxy[iZone]; dq_proxy[iZone] = NULL;
        }
        delete [] dq_proxy; dq_proxy = NULL;
        return;
    }

    void MYReactor::ScanEquationForwardly(int level, FieldProxy ** dq_proxy)
    {
        RDouble sweepNormal;
        for (int iZone = 0; iZone < numberOfZones; ++ iZone)
        {
            Zone * zone = GetZone(iZone);
            if (! zone) continue;
            Grid *grid = PHSPACE::GetGrid(iZone, level);
            CFDSolver * solver = GetCFDSolver(iZone);
            dq_proxy[iZone] = solver->CreateFieldProxy(grid);

            TK_Exit::UnexpectedVarValue("Warning: please check the right parameter of SolveLUSGSForward function!(Bell note)", 0);
            solver->SolveLUSGSForward(grid, dq_proxy[iZone], 0, sweepNormal);
        }
        return;
    }

    int MYReactor::FindZoneSeriesOnNeighbor(int jZone, int iZone)
    {
        ZoneNeighbor * neighbor = zoneConnectivity->GetZoneNeighbor(jZone);
        int numberOfNeighbors = neighbor->GetNumberOfNeighbor();
        for (int iNeighbor = 0; iNeighbor < numberOfNeighbors; ++ iNeighbor)
        {
            if (neighbor->GetZoneIndexOfNeighbor(iNeighbor) == iZone) return iNeighbor;
        }
        return MB_SPACE::GHOST_INDEX;
    }

    void MYReactor::TradeContactDeltaField(int level, FieldProxy ** dq_proxy)
    {
        if (this->solverIndex != 0) return;

        int rank = PHMPI::GetCurrentProcessorID();
        vector< MPI_Request > request;
        vector< DataContainer * > * dataForSend = new vector< DataContainer * > [numberOfZones];
        vector< DataContainer * > * dataForRecv = new vector< DataContainer * > [numberOfZones];

        for (int iZone = 0; iZone < numberOfZones; ++ iZone)
        {
            int ip = PHMPI::GetZoneProcessorID(iZone);

            ZoneNeighbor * zone_neighbor = zoneConnectivity->GetZoneNeighbor(iZone);
            int numberOfNeighbors = zone_neighbor->GetNumberOfNeighbor();
            if(rank == ip)
            {
                CFDSolver * solver = GetCFDSolver(iZone);
                Grid *grid = PHSPACE::GetGrid(iZone, level);

                dataForSend[iZone].resize(numberOfNeighbors);
                for (int iNeighbor = 0; iNeighbor < numberOfNeighbors; ++ iNeighbor)
                {
                    int jZone = zone_neighbor->GetZoneIndexOfNeighbor(iNeighbor);
                    dataForSend[iZone][iNeighbor] = new DataContainer();
                    solver->CompressDQ(dataForSend[iZone][iNeighbor], dq_proxy[iZone], grid, iZone, jZone, 5);
                }

                dataForRecv[iZone].resize(numberOfNeighbors);
                for (int iNeighbor = 0; iNeighbor < numberOfNeighbors; ++ iNeighbor)
                {
                    dataForRecv[iZone][iNeighbor] = new DataContainer();
                    dataForRecv[iZone][iNeighbor]->resize(dataForSend[iZone][iNeighbor]->size());
                }
            }
        }

        for (int iZone = 0; iZone < numberOfZones; ++ iZone)
        {
            int tag = iZone;
            int ip = PHMPI::GetZoneProcessorID(iZone);

            ZoneNeighbor * zone_neighbor = zoneConnectivity->GetZoneNeighbor(iZone);
            int numberOfNeighbors = zone_neighbor->GetNumberOfNeighbor();
            for (int iNeighbor = 0; iNeighbor < numberOfNeighbors; ++ iNeighbor)
            {
                int jZone = zone_neighbor->GetZoneIndexOfNeighbor(iNeighbor);
                int jp    = PHMPI::GetZoneProcessorID(jZone);

                if (rank == ip && rank == jp)
                {
                    int series = FindZoneSeriesOnNeighbor(jZone, iZone);

                    DataContainer * address      = dataForSend[iZone][iNeighbor];
                    dataForSend[iZone][iNeighbor] = dataForRecv[jZone][series];
                    dataForRecv[jZone][series]    = address;

                    continue;
                }

                if (rank == ip)
                {
                    streamsize length = dataForSend[iZone][iNeighbor]->size();
                    send(dataForSend[iZone][iNeighbor], jp, length, request, tag);
                }
                else if (rank == jp)
                {
                    int series = FindZoneSeriesOnNeighbor(jZone, iZone);
                    streamsize length = dataForRecv[jZone][series]->size();
                    receive(dataForRecv[jZone][series], ip, length, request, tag);
                }
            }
        }

        if(PHMPI::GetNumberOfProcessor() > 1)
        {
            PH_Wait(static_cast< int >(request.size()), & request[0]);
        }

        for (int iZone = 0; iZone < numberOfZones; ++ iZone)
        {
            int ip = PHMPI::GetZoneProcessorID(iZone);

            ZoneNeighbor * zone_neighbor = zoneConnectivity->GetZoneNeighbor(iZone);
            int numberOfNeighbors = zone_neighbor->GetNumberOfNeighbor();
            if (rank == ip)
            {
                CFDSolver * solver = GetCFDSolver(iZone);
                Grid *grid = PHSPACE::GetGrid(iZone, level);
                for (int iNeighbor = 0; iNeighbor < numberOfNeighbors; ++ iNeighbor)
                {
                    int jZone = zone_neighbor->GetZoneIndexOfNeighbor(iNeighbor);
                    solver->DecompressDQ(dataForRecv[iZone][iNeighbor], dq_proxy[iZone], grid, iZone, jZone, 5);
                }
            }
        }

        for (int iZone = 0; iZone < numberOfZones; ++ iZone)
        {
            int numberOfNeighbors = static_cast< int >(dataForSend[iZone].size());
            for (int iNeighbor = 0; iNeighbor < numberOfNeighbors; ++ iNeighbor)
            {
                delete dataForSend[iZone][iNeighbor]; dataForSend[iZone][iNeighbor] = NULL;
                delete dataForRecv[iZone][iNeighbor]; dataForRecv[iZone][iNeighbor] = NULL;
            }
        }
        delete [] dataForSend; dataForSend = NULL;
        delete [] dataForRecv; dataForRecv = NULL;

        return;
    }

    void MYReactor::ScanEquationBackwardly(int level, FieldProxy ** dq_proxy)
    {
        RDouble sweepNormal;
        for (int iZone = 0; iZone < numberOfZones; ++ iZone)
        {
            Zone * zone = GetZone(iZone);
            if (! zone) continue;
            Grid *grid = PHSPACE::GetGrid(iZone, level);
            CFDSolver * solver = GetCFDSolver(iZone);
            solver->SolveLUSGSBackward(grid, dq_proxy[iZone], 0, sweepNormal);
        }
        return;
    }

    void MYReactor::UpdateSteadyFlowField(int level, FieldProxy ** dq_proxy)
    {
        for (int iZone = 0; iZone < numberOfZones; ++ iZone)
        {
            Zone * zone = GetZone(iZone);
            if (! zone) continue;
            Grid *grid = PHSPACE::GetGrid(iZone, level);
            CFDSolver * solver = GetCFDSolver(iZone);
            FieldProxy * q_proxy = solver->GetFieldProxy(grid, solver->CastName("q"));
            solver->UpdateFlowField(grid, q_proxy, dq_proxy[iZone]);
            delete q_proxy; q_proxy = NULL;
        }
        return;
    }

    void MYReactor::InitResidualByNegativeRhs(int level)
    {
        FieldProxy ** rhs_proxy = GetRHSProxy(level);
        for (int iZone = 0; iZone < numberOfZones; ++ iZone)
        {
            Zone * zone = GetZone(iZone);
            if (! zone) continue;
            Grid *grid = PHSPACE::GetGrid(iZone, level);
            CFDSolver * solver = GetCFDSolver(iZone);
            solver->InitResidual(grid, rhs_proxy[iZone]);
        }
        return;
    }

    void MYReactor::LoadFlowField(int level)
    {
        FieldProxy ** q_proxy = GetQProxy(level);
        for (int iZone = 0; iZone < numberOfZones; ++ iZone)
        {
            Zone * zone = GetZone(iZone);
            if (! zone) continue;
            Grid *grid = PHSPACE::GetGrid(iZone, level);
            CFDSolver * solver = GetCFDSolver(iZone);
            q_proxy[iZone] = solver->CreateFieldProxy(grid);
            solver->LoadQ(grid, q_proxy[iZone]);
        }
        return;
    }

    void MYReactor::RestrictCoarseFlowField(int level)
    {
        for (int iZone = 0; iZone < numberOfZones; ++ iZone)
        {
            Zone * zone = GetZone(iZone);
            if (! zone) continue;
            Grid * fine_grid = PHSPACE::GetGrid(iZone, level);
            Grid * coarse_grid = fine_grid->GetCoarseGrid();
            CFDSolver * solver = GetCFDSolver(iZone);
            solver->RestrictAllQ(fine_grid, coarse_grid);
        }
        return;
    }

    void MYReactor::RestrictDefectOfModel(int level)
    {
        for (int iZone = 0; iZone < numberOfZones; ++ iZone)
        {
            Zone * zone = GetZone(iZone);
            if (! zone) continue;
            Grid * fine_grid   = PHSPACE::GetGrid(iZone, level);
            Grid * coarse_grid = fine_grid->GetCoarseGrid();
            CFDSolver * solver = GetCFDSolver(iZone);
            solver->RestrictDefect(fine_grid, coarse_grid);
        }
        return;
    }

    void MYReactor::PutNegativeCorrection(int level)
    {
        FieldProxy ** q_proxy = GetQProxy(level);
        for (int iZone = 0; iZone < numberOfZones; ++ iZone)
        {
            Zone * zone = GetZone(iZone);
            if (! zone) continue;
            Grid *grid = PHSPACE::GetGrid(iZone, level);
            CFDSolver * solver = GetCFDSolver(iZone);
            solver->PutCorrection(grid, q_proxy[iZone]);
        }
        return;
    }

    void MYReactor::PutPositiveCorrection(int level)
    {
        FieldProxy ** q_proxy = GetQProxy(level);
        for (int iZone = 0; iZone < numberOfZones; ++ iZone)
        {
            Zone * zone = GetZone(iZone);
            if (! zone) continue;
            Grid *grid = PHSPACE::GetGrid(iZone, level);
            CFDSolver * solver = GetCFDSolver(iZone);
            solver->PutCorrectionBack(grid, q_proxy[iZone]);
            delete q_proxy[iZone]; q_proxy[iZone] = NULL;
        }
        return;
    }

    void MYReactor::CorrectFlowFieldOnFineGrid(int level)
    {
        for (int iZone = 0; iZone < numberOfZones; ++ iZone)
        {
            Zone * zone = GetZone(iZone);
            if (! zone) continue;
            Grid *grid = PHSPACE::GetGrid(iZone, level);
            CFDSolver * solver = GetCFDSolver(iZone);
            solver->CorrectFineGrid(grid, grid->GetCoarseGrid());
        }
        return;
    }

    void MYReactor::InitRhsByPositiveResidual(int level)
    {
        FieldProxy ** rhs_proxy = GetRHSProxy(level);
        for (int iZone = 0; iZone < numberOfZones; ++ iZone)
        {
            Zone * zone = GetZone(iZone);
            if (! zone) continue;
            Grid *grid = PHSPACE::GetGrid(iZone, level);
            CFDSolver * solver = GetCFDSolver(iZone);
            rhs_proxy[iZone] = solver->CreateFieldProxy(grid);
            solver->StoreRhsByResidual(grid, rhs_proxy[iZone]);
        }
        return;
    }

    void MYReactor::RecoverResidualOnFineGrid(int level)
    {
        FieldProxy ** rhs_proxy = GetRHSProxy(level);
        for (int iZone = 0; iZone < numberOfZones; ++ iZone)
        {
            Zone * zone = GetZone(iZone);
            if (! zone) continue;
            Grid *grid = PHSPACE::GetGrid(iZone, level);
            CFDSolver * solver = GetCFDSolver(iZone);
            solver->RecoverResidual(grid, rhs_proxy[iZone]);
            delete rhs_proxy[iZone]; rhs_proxy[iZone] = NULL;
        }
        return;
    }

    bool MYReactor::IsLocallyConvergent(int level)
    {
        int innstep = GlobalDataBase::GetIntParaFromDB("innstep");
        RDouble tol_sub_iter = GlobalDataBase::GetDoubleParaFromDB("tol_sub_iter");
        int min_sub_iter = GlobalDataBase::GetIntParaFromDB("min_sub_iter");
        int max_sub_iter = GlobalDataBase::GetIntParaFromDB("max_sub_iter");

        RDouble cvgmax = 0.0, cvg;
        for (int iZone = 0; iZone < numberOfZones; ++ iZone)
        {
            Zone * zone = GetZone(iZone);
            if (!zone) continue;
            Grid *grid = PHSPACE::GetGrid(iZone, level);
            CFDSolver * solver = GetCFDSolver(iZone);
            cvg = solver->UnsteadyConvergence(grid);
            cvgmax = max(cvgmax, cvg);
        }
        RunGlobalTaskList();
        PH_CompareMaxMin(cvgmax, 1);

        if (innstep <  min_sub_iter) return false;
        if (cvgmax  <= tol_sub_iter) return true ;
        if (innstep >= max_sub_iter) return true ;

        return false;
    }

    void MYReactor::TradeContactFlowField(int level)
    {
        for (int iZone = 0; iZone < numberOfZones; ++ iZone)
        {
            Zone * zone = GetZone(iZone);
            if (! zone) continue;
            Grid *grid = PHSPACE::GetGrid(iZone,level);
            CFDSolver * solver = GetCFDSolver(iZone);
            solver->PostSolve(grid, 0);
        }
        RunGlobalTaskList();
        return;
    }

}