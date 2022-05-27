#include "UGKSSimulation.h"
#include "Post_WriteTecplot.h"
#include "TK_Time.h"
#include "UGKSSolver.h"
#include "Region.h"
#include "Math_BasisFunction.h"
#include "Task.h"
#include "Zone.h"
#include "omp.h"

namespace PHSPACE
{
    UGKSSimulation::UGKSSimulation()
    {
        region = new Region();
        simu_ts = NULL;
    }
    
    UGKSSimulation::~UGKSSimulation()
    {
        using namespace PHMPI;
        delete region;
        region = NULL;
    }
    
    void UGKSSimulation::InitGeometry()
    {
        region->InitGeometry();
    }
    
    void UGKSSimulation::ReadGeometry()
    {
        TimeSpan* testTime = new TimeSpan();
        ReadGrid();
        testTime->ShowTimeSpanToLogFile("Grid Reading");
        delete testTime;
    }
    
    void UGKSSimulation::ReadGrid()
    {
        region->ReadGrid();
    }
    
    void UGKSSimulation::UpdateAllZonesPara()
    {
        this->region->UpdateAllZonesPara();
    }
    
    void UGKSSimulation::Run()
    {
        //cout << PHMPI::GetCurrentProcessorID() << endl;
        //omp_set_num_threads(1);
        //RDouble start = omp_get_wtime();
        //TimeSpan* testTime = new TimeSpan();
        Preprocess();
        //testTime->ShowTimeSpanToLogFile("Process");
        //delete testTime;

        SolveSimulation();
    }
    
    void UGKSSimulation::Preprocess()
    {
        ReadGeometry();

        UpdateAllZonesPara();

        InitGeometry();

        DefineAndCreateSolvers();

        InitializeSolvers();

        //InitGlobalBoundaryCondition();
    }
    
    void UGKSSimulation::DefineAndCreateSolvers()
    {
        Region * this_region = this->region;

        int nZones = PHMPI::GetNumberofGlobalZones();

        vector<vector<string> > allSolverNames;
        allSolverNames.resize(nZones);

        for (int iZone = 0; iZone < nZones; ++ iZone)
        {
            Zone * zone = this_region->GetZone(iZone); 
            if (!zone) continue;
            PHGeometry * geometry = zone->GetGeometry();

            vector<PHSolver *> solverList; 
            int gridType = geometry->GetGrid()->Type();
            if (gridType == STRUCTGRID)
            {
                PHSolver *solver = new UGKSSolver();
                solverList.push_back(solver);
            }
            else if (gridType == UNSTRUCTGRID)
            {
                cout << "unstructured grid of UGKS is not defined";
            }
            //! Finally, push the solvers into the list.
            for (std::size_t iSolver = 0; iSolver < solverList.size(); ++ iSolver)
            {
                PHSolver *solver = solverList[iSolver];

                solver->SetGeometry(geometry);
                solver->SetIndex(static_cast<int>(iSolver));
                zone->AddSolver(solver);
            }
        }
        InitCommonParameter();
        //GlobalBoundaryCondition::ReadGlobalBoundaryCondition();
        //GlobalBoundaryCondition::SetGlobalDataBaseByGlobalBC();
    }
    
    void UGKSSimulation::InitCommonParameter()
    {
        RDouble reynoldsReferenceLengthDimensional = 1.0;

        RDouble gridScaleFactor = GlobalDataBase::GetDoubleParaFromDB("gridScaleFactor");

        reynoldsReferenceLengthDimensional = gridScaleFactor;

        GlobalDataBase::UpdateData("reynoldsReferenceLengthDimensional", &reynoldsReferenceLengthDimensional, PHDOUBLE, 1);

    }
    
    void UGKSSimulation::InitializeSolvers()
    {
        WriteLogFile("Creating and initialization solvers start\n");

        using namespace PHMPI;
        int nLocalZones = GetNumberofLocalZones();
        int iSolver = 0; //! Only one solver.

        //! First: initialize memory.
        for (int iZone = 0; iZone < nLocalZones; ++ iZone)
        {
            int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

            Zone *zone = region->GetZone(zoneID);
            if (!zone) continue;
            PHSolver *solver = static_cast<UGKSSolver *>(GlobalSolvers::GetSolver(zoneID, iSolver));
            solver->InitMemory();
        }
        //! Second: initialize flow variables on the finest mesh.
        for (int iZone = 0; iZone < nLocalZones; ++ iZone)
        {
            int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

            Zone *zone = region->GetZone(zoneID);
            if (!zone) continue;
            PHSolver *solver = static_cast<UGKSSolver *>( GlobalSolvers::GetSolver(zoneID, iSolver));

            solver->InitFlow();
        }
        //! Finally: interfaces communications.
        for (int iZone = 0; iZone < nLocalZones; ++ iZone)
        {
            int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

            Zone *zone = region->GetZone(zoneID);
            if (!zone) continue;
            UGKSSolver *solver = static_cast<UGKSSolver *>( GlobalSolvers::GetSolver(zoneID, iSolver));

            solver->CommunicationInterfaceDataInit();
        }
        //GlobalSolve();

        WriteLogFile("Creating and initialization solvers over\n");
    }
    
    void UGKSSimulation::InitGlobalBoundaryCondition()
    {
        this->region->ReSetBoundaryConditionByGlobalBC();
    }
    
    void UGKSSimulation::SolveSimulation()
    {
        WriteLogFile("Start iterating ...");
        int maxSimuStep = PHSPACE::GlobalDataBase::GetIntParaFromDB("maxSimuStep");
        bool continueSolve = true;
        int outnstep = PHSPACE::GlobalDataBase::GetIntParaFromDB("outnstep");
        while (continueSolve)
        {
            outnstep = outnstep + 1;
            GlobalDataBase::UpdateData("outnstep", &outnstep, PHINT, 1);
            //RDouble start = omp_get_wtime();
            SolveOneStep();
            //RDouble end = omp_get_wtime();
            //RDouble diff = end - start;
            if (outnstep >= maxSimuStep && maxSimuStep > 0) continueSolve = false;
        }
    }
    
    void UGKSSimulation::SolveOneStep()
    {
        MainSolve();

        //TimeSpan* testTime = new TimeSpan();
        //RDouble start, end, diff;
        //start = omp_get_wtime();
        PostSolve();
        //end = omp_get_wtime();
        //diff = end - start;
        //cout <<"postsolve:"<<diff << endl;
        //testTime->ShowTimeSpanToLogFile("Process Time");
        //delete testTime;
    }
    
    void UGKSSimulation::MainSolve()
    {
        ComputeMinTimeStep();
        UGKSSolve();
    }
    
    void UGKSSimulation::ComputeMinTimeStep()
    {
        using namespace PHMPI;
        int nLocalZones = GetNumberofLocalZones();
        int level = 0;
        int iSolver = 0; //only one solver

        RDouble globalMinDt = LARGE,minDt;
        for (int iZone = 0; iZone < nLocalZones; ++ iZone)
        {
            int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

            Zone *zone = region->GetZone(zoneID);
            if (!zone) continue;
            Grid *grid = GetGrid(zoneID, level);
            UGKSSolver *solver = static_cast<UGKSSolver *>(GlobalSolvers::GetSolver(zoneID, iSolver));
            solver->ComputeMinTimeStep(grid, minDt);
            globalMinDt = MIN(globalMinDt, minDt);
        }
        PH_CompareMaxMin(globalMinDt, 2);
        GlobalDataBase::UpdateData("globalMinTimeStep", &globalMinDt, PHDOUBLE, 1);

        if (PHMPI::GetCurrentProcessorID() == PHMPI::GetServerProcessorID())
        {
            int outnstep = GlobalDataBase::GetIntParaFromDB("outnstep");
            ostringstream oss;
            oss.setf(ios::scientific);
            oss.precision(2);
            oss << "    iter = " << outnstep 
                << ", level = " << level
                << ", min time step = " << globalMinDt;
            WriteLogFile(oss);
        }
    }
    
    void UGKSSimulation::UGKSSolve()
    {
        using namespace PHMPI;
        int nLocalZones = GetNumberofLocalZones();
        int level = 0;
        int iSolver = 0; //only one solver

        for (int iZone = 0; iZone < nLocalZones; ++ iZone)
        {
            int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

            Zone *zone = region->GetZone(zoneID);
            if (!zone) continue;
            Grid *grid = GetGrid(zoneID, level);
            UGKSSolver *solver = static_cast<UGKSSolver *>(GlobalSolvers::GetSolver(zoneID, iSolver));

            solver->UGKSSolve(grid);
        }
    }
    
    void UGKSSimulation::PostSolve()
    {
        using namespace PHMPI;
        int nLocalZones = GetNumberofLocalZones();
        int level = 0;
        int iSolver = 0; //only one solver

        for (int iZone = 0; iZone < nLocalZones; ++ iZone)
        {
            int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

            Zone *zone = region->GetZone(zoneID);
            if (!zone) continue;
            Grid *grid = GetGrid(zoneID, level);
            UGKSSolver *solver = static_cast<UGKSSolver *>(GlobalSolvers::GetSolver(zoneID, iSolver));

            solver->PostSolve(grid);
        }
        GlobalSolve();

        Visualize();
    }
    
    void UGKSSimulation::Visualize()
    {
        int intervalStepPlot = GlobalDataBase::GetIntParaFromDB("intervalStepPlot");
        int outnstep  = GlobalDataBase::GetIntParaFromDB("outnstep");
        if (outnstep % intervalStepPlot == 0)
        {
            using namespace PHMPI;
            int nLocalZones = GetNumberofLocalZones();
            int level = 0;
            int iSolver = 0; //only one solver
            for (int iZone = 0; iZone < nLocalZones; ++ iZone)
            {
                int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

                Zone *zone = region->GetZone(zoneID);
                if (!zone) continue;
                Grid *grid = GetGrid(zoneID, level);
                UGKSSolver *solver = static_cast<UGKSSolver *>(GlobalSolvers::GetSolver(zoneID, iSolver));

                solver->VisulizeGhostCell(grid);
            }  

            Post_WriteTecplot generateVisualFile;
            generateVisualFile.Run();
        }
    }
    
    void UGKSSimulation::Postprocess()
    {
        this->region->PostSimulation();
    }
}