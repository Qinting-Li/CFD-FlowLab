#include "Simulation.h"
#include "IO_FileName.h"
#include "TK_Parse.h"
#include "GridFactory.h"
#include "AleModel.h"
#include "GridPartition.h"
#include "Controller.h"
#include "TK_Time.h"
#include "PHIO.h"
#include "Geo_Probes.h"
#include "HOSimulation.h"

#ifdef BUILD_LBM
#include "cavity.h"
#endif

#include "OversetInformation.h"
#include "SolversList.h"
#include "PostProcess.h"
#include "Glb_Dimension.h"
#include "Zone.h"
#include "WallDist.h"
#include "Post_WriteTecplot.h"
#include "Post_ProbesWrite.h"
using namespace std;

namespace PHSPACE
{
Simulation::Simulation(int dimension)
{
    using namespace PHMPI;

    SetDim(dimension);

    region = new Region();
}

Simulation::~Simulation()
{
    using namespace PHMPI;

    delete region;

    TK_Exit::ExitPHengLEI();
}

void Simulation::Start()
{
    int taskType = GetTaskCode();

    switch (taskType)
    {
    case SOLVE_FIELD:
        RunCFD();
        break;
    case CAL_WALL_DIST:
        ComputeWallDistance();
        break;
    case CREATE_GRID:
        PHSPACE::CreateGrid();
        break;
    case PARTITION_GRID:
        IsNeedConvert();
        PHSPACE::PartitionGrid();
        break;
    case OVERSETGRID_VIEW:
        PHSPACE::OversetGridView();
        break;
    case HO_SOLVER:
        {
            HOUnstruct::HOSimulation * uho = new HOUnstruct::HOSimulation();
            uho->Run();
            delete uho;
            break;
        }
    case CONVERT_PARTITION:
        //ConvertAndPartition();
        break;
	#ifdef BUILD_LBM
	case LBM_SOLVER:
        RunLBMSolver();
        break;
	#endif
    case INTEGRATIVE_SOLVER:
        RunIntegrativeSolver();
        break;
    case POST_PROCESSING:
        RunPostProcessing();
        break;
    default:
        TK_Exit::UnexpectedVarValue("taskType = ", taskType);
        break;
    }
}

void Simulation::IsNeedConvert()
{
    string originalGridFile = GlobalDataBase::GetStrParaFromDB("original_grid_file");

    string fileMame, fileSuffix;
    GetNameExt(originalGridFile, fileMame, fileSuffix, ".");

    set < pair<string, int> > gridTypeMap;
    gridTypeMap.insert(pair <string, int> ("fts",  FANTASY_TYPE  ));
    gridTypeMap.insert(pair <string, int> ("cgns", CGNS_TYPE     ));
    gridTypeMap.insert(pair <string, int> ("grd",  GRIDGEN_TYPE  ));
    gridTypeMap.insert(pair <string, int> ("uns",  FIELDVIEW_TYPE));
    gridTypeMap.insert(pair <string, int> ("cas",  FLUENT_TYPE   ));

    int sourceGridType = 0;
    set< pair<string, int> >::iterator iter;
    for (iter = gridTypeMap.begin(); iter != gridTypeMap.end(); ++ iter)
    {
        //if ((*iter).first != fext) continue;
        //if (strcmp(fext, (*iter).first)
        if (fileSuffix.compare((*iter).first) != 0)
        {
            continue;
        }

        sourceGridType = (*iter).second;
        break;
    }

    if (sourceGridType == 0)
    {
        ostringstream oss;
        oss << "  Error: this situation has not been considered, for filetype = ." << fileSuffix << endl;
        TK_Exit::ExceptionExit(oss.str());
    }

    if (sourceGridType != FANTASY_TYPE)
    {
        GlobalDataBase::UpdateData("from_gtype", &sourceGridType, PHINT, 1);
        GlobalDataBase::UpdateData("from_gfile", &originalGridFile, PHSTRING, 1);

        int targetGridType = 1;
        GlobalDataBase::UpdateData("to_gtype", &targetGridType, PHINT, 1);

        string outputGridFile = ChangeExtensionOfFileName(originalGridFile, "fts");
        GlobalDataBase::UpdateData("out_gfile", &outputGridFile, PHSTRING, 1);

        int partitionedGridType = GlobalDataBase::GetIntParaFromDB("pgridtype");
        GlobalDataBase::UpdateData("gridtype", &partitionedGridType, PHINT, 1);

        int gridObjection = 1;
        GlobalDataBase::UpdateData("gridobj", &gridObjection, PHINT, 1);

        int numberOfSimulationTask = PHSPACE::CREATE_GRID;
        GlobalDataBase::UpdateData("nsimutask", &numberOfSimulationTask, PHINT, 1);

        int multiBlock = 1;
        GlobalDataBase::UpdateData("multiblock", &multiBlock, PHINT, 1);

        int axisUp = 1;
        GlobalDataBase::UpdateData("axisup", &axisUp, PHINT, 1);

        int gridDatabaseIndex = 3;
        GlobalDataBase::UpdateData("grid_database_index", &gridDatabaseIndex, PHINT, 1);

        int iadapt = 0;
        GlobalDataBase::UpdateData("iadapt", &iadapt, PHINT, 1);

        int iovrlap = 0;
        GlobalDataBase::UpdateData("iovrlap", &iovrlap, PHINT, 1);

        CompleteParamInfo();

        PHSPACE::CreateGrid();

        GlobalDataBase::UpdateData("original_grid_file", &outputGridFile, PHSTRING, 1);

        numberOfSimulationTask = PHSPACE::PARTITION_GRID;
        GlobalDataBase::UpdateData("nsimutask", &numberOfSimulationTask, PHINT, 1);
    }
}

void Simulation::RunCFD()
{
    double last_time[3];
    TIME_SPACE::StartTime(last_time);

    ReadGeometry();

    UpdateAllZonesPara();

    InitGeometry();

    SwapGeommetry();

    //TODO: To develop in the future
    InitOversetGrid();

    if (PHMPI::CurrentProcessorIsGridProcessor())
    {
        return;
    }

    //! Define solver names for all zones and create the new defined solvers.
    //! IMPORTANT for solver developers:
    //! Define your solver name which is same with the solver class name here!
    DefineAndCreateSolvers();

    InitGlobalBoundaryCondition();

    GenerateWallDist();

    //! Initialize each solver on each zone.
    InitializeSolvers();

    MultiGridInitFlow();

    //! The main simulation entrance.
	SolveSimulation();

    //! Do some work after simulation.
	PostSimulation();
}

void Simulation::RunIntegrativeSolver()
{
    int numberOfProcessor = PHMPI::GetNumberOfProcessor();
    int currentProcessor  = PHMPI::GetCurrentProcessorID();
    int serverProcessor   = PHMPI::GetServerProcessorID();
    int tasktype;

    string fromGridName = GlobalDataBase::GetStrParaFromDB("from_gfile");
    string outGridName = ChangeExtensionOfFileName(fromGridName, "fts");

    tasktype = PHSPACE::CREATE_GRID;
    GlobalDataBase::UpdateData("out_gfile", &outGridName, PHSTRING, 1);
    GlobalDataBase::UpdateData("nsimutask", &tasktype, PHINT, 1);

    int isConvertOver = 0;
    PHMPI::SetNumberOfProcessor(1);

    if (currentProcessor == serverProcessor)
    {
        CreateGrid();
        isConvertOver = 1;
    }

    PHMPI::SetNumberOfProcessor(numberOfProcessor);
    PHMPI::PH_Bcast(&isConvertOver, 1, serverProcessor);
    WriteLogFile("Convert grid complete !\n");
    PrintToWindow("Convert grid complete !\n\n");

    //! If need PartitionGrid.
    if (numberOfProcessor > 1)
    {
        WriteLogFile("Begin partition grid !\n");
        PrintToWindow("Begin partition grid !\n\n");
        PHMPI::SetNumberOfProcessor(1);

        tasktype = PHSPACE::PARTITION_GRID;
        int maxproc = numberOfProcessor;
        int numberOfMultigrid = 1;
        int gridtype = GlobalDataBase::GetIntParaFromDB("gridtype");
        string partition_grid_file = AddSymbolToFileName(outGridName, "__", numberOfProcessor);

        GlobalDataBase::UpdateData("nsimutask", &tasktype, PHINT, 1);
        GlobalDataBase::UpdateData("maxproc", &maxproc, PHINT, 1);
        GlobalDataBase::UpdateData("pgridtype", &gridtype, PHINT, 1);
        GlobalDataBase::UpdateData("numberOfMultigrid", &numberOfMultigrid, PHINT, 1);
        GlobalDataBase::UpdateData("partition_grid_file", &partition_grid_file, PHSTRING, 1);
        GlobalDataBase::UpdateData("original_grid_file", &outGridName, PHSTRING, 1);

        int isPartitionOver = 0;
        if (currentProcessor == serverProcessor)
        {
            PartitionGrid();
            isPartitionOver = 1;

            DeleteGGrids();
            PHMPI::SetNumberofLocalZones(0);
        }
        else
        {
            bool fileExist = true;
            bool fileCheck = true;
            PH_AllReduce(&fileExist, &fileCheck, 1, MPI_MIN);
        }

        PHMPI::SetNumberOfProcessor(numberOfProcessor);
        PHMPI::PH_Bcast(&isPartitionOver, 1, serverProcessor);
        WriteLogFile("Partition grid complete !\n");
        PrintToWindow("Partition grid complete !\n\n");

        string walldistfile = ChangeExtensionOfFileName(partition_grid_file, "wdt");
        GlobalDataBase::UpdateData("gridfile", &partition_grid_file, PHSTRING, 1);
        GlobalDataBase::UpdateData("walldistfile", &walldistfile, PHSTRING, 1);
    }
    else
    {
        string walldistfile = ChangeExtensionOfFileName(outGridName, "wdt");
        GlobalDataBase::UpdateData("gridfile", &outGridName, PHSTRING, 1);
        GlobalDataBase::UpdateData("walldistfile", &walldistfile, PHSTRING, 1);
    }

    tasktype = PHSPACE::SOLVE_FIELD;
    GlobalDataBase::UpdateData("nsimutask", &tasktype, PHINT, 1);

    RunCFD();
}

void Simulation::RunPostProcessing()
{
    PostProcess *postProcess = new PostProcess(this->region);
    postProcess->Run();

    delete postProcess;

    postProcess = NULL;

}


void Simulation::ComputeWallDistance()
{
    double last_time[3];
    TIME_SPACE::StartTime(last_time);

    ReadGeometry();

    UpdateAllZonesPara();

    InitGeometry();

    GlobalBoundaryCondition::ReadGlobalBoundaryCondition();

    GlobalBoundaryCondition::ChangeBCTypeByGlobalBC();

    InitGlobalBoundaryCondition();

    GenerateWallDist();

    WriteLogFile("GenerateWallDist ...");

}

void Simulation::InitOversetGrid()
{
    //! The overset config process must be taken after reading in ordinary grid and initializing!!!
    //! Read overset information files or config.

    int iOverlap = PHSPACE::GlobalDataBase::GetIntParaFromDB("iovrlap");

    if (!iOverlap)
    {
        return;
    }

    //! Structured overset mesh! need to be revised!
    int sys_gridtype;
    GlobalDataBase::GetData("sys_gridtype",&sys_gridtype,PHINT,1);

    if (sys_gridtype == STRUCTGRID)      //! Overset system configuration????
    {
        region->ReadOversetGrid();

        PHSPACE::GenerateStructGridTable();

        PHSPACE::CreateOversetStructGrid();

        OversetStructGrid * oversetStructGrid = PHSPACE::GetOversetStructGrid();

        oversetStructGrid->RunGeometricAnalysis();

        return;
    }

    int readOversetFileOrNot = PHSPACE::GlobalDataBase::GetIntParaFromDB("readOversetFileOrNot");

    //! if the processor is grid processor, and overset grid need to be assembled.
    if (PHMPI::CurrentProcessorIsGridProcessor() && readOversetFileOrNot != 1)
    {
        OversetGridConfig();
        WriteLogFile("Run overset configuration\n");
    }

    //! waite for the grid processors.
    PH_Barrier();

    //! for the original processors
    if (!PHMPI::CurrentProcessorIsGridProcessor())
    {
        //! if the overset files is aviable or has been prepared by grid processors, just read in it!
        if (readOversetFileOrNot == 1 || PHMPI::GetNumberOfGridProcessor() != 0)
        {
            region->ReadOversetGrid();
            WriteLogFile("Read in overset information\n");
        }
        //! else it should be done by these processors by themselves.
        else
        {
            OversetGridConfig();
            WriteLogFile("Run overset configuration\n");
        }

        InitializeOversetInterfaceTopology();

        WriteLogFile("Initialize OversetInterface Topology\n");
    }
}

void Simulation::OversetGridConfig()
{
    int sysGridType;
    GlobalDataBase::GetData("sys_gridtype",&sysGridType,PHINT,1);

    if (sysGridType == STRUCTGRID)      //! Overset system configuration????
    {
        PHSPACE::GenerateStructGridTable();

        PHSPACE::CreateOversetStructGrid();

        OversetStructGrid *oversetStructGrid = PHSPACE::GetOversetStructGrid();

        oversetStructGrid->RunGeometricAnalysis();
    }
    else if (sysGridType == UNSTRUCTGRID)
    {
    }
}

void Simulation::ReadGeometry()
{
    TimeSpan testTime;

    this->region->ReadGrid();

    testTime.ShowTimeSpanToLogfile("Grid Reading");
}

void Simulation::InitGeometry()
{
    this->region->InitGeometry();
}

void Simulation::SwapGeommetry()
{
    this->region->SwapGeommetry();
}

void Simulation::InitializeSolvers()
{
    WriteLogFile("Creating and initialization solvers start\n");

    int nSolver = Zone::GetNumberOfSolvers();
    for (int iSolver = 0; iSolver < nSolver; ++ iSolver)
    {
        Controller *controller = new Controller(iSolver, this->region);

        controller->Initialize(iSolver);

        delete controller;
    }

    WriteLogFile("Creating and initialization solvers over\n");
}

void Simulation::DefineAndCreateSolvers()
{
    Region *this_region = this->region;
    Zone::Init();

    Region *region = this->region;

    int nZones = PHMPI::GetNumberofGlobalZones();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        Zone *zone = region->GetZone(iZone); 
        if (!zone)
        {
            continue;
        }

        vector<PHSolver *> solverList;
#ifdef USE_NSDEMO
        solver = new NSSolverDemo();
        solverList.push_back(solver);
#else  
        PHGeometry *geometry = zone->GetGeometry();        

        //! For each zone, define the solver names on it.
        //! If you are solver developer, define your solve name here,
        //! according your solver logic.
        int dgHighOrder = 0;
        zone->GetData("dg_high_order",&dgHighOrder, PHINT, 1);  

        //! First: N-S equations solver.
        if (!dgHighOrder)
        {
            int gridType = geometry->GetGrid()->Type();

            int compressible  = 1;
            zone->GetData("compressible",&compressible, PHINT, 1);

            if (compressible)
            {
                //! Add compressible solvers.
                if (gridType == UNSTRUCTGRID)
                {
                    //! Add unstructured NS solver.
                    PHSolver *solver = new NSSolverUnstruct();
                    solverList.push_back(solver);
                }
                else if (gridType == STRUCTGRID)
                {
                    const int FV_METHOD = 0;
                    const int FD_METHOD = 1;

                    int FVorFD = 0;
                    zone->GetData("ifvfd",&FVorFD, PHINT, 1);

                    if (FVorFD == FV_METHOD)
                    {
                        //! Add structured FV NS solver.
                        PHSolver *solver = new NSSolverStructFV();
                        solverList.push_back(solver);
                    }
                    else
                    {
                        //! Add structured FD NS solver.
                        PHSolver *solver = new NSSolverStructFD();
                        solverList.push_back(solver);
                    }
                }			    
            }
            else
            {
            }
        }
        else
        {
            //
        }

        //! Second: turbulent equations solver.
        int viscousType = ONE_EQU;
        zone->GetData("viscousType", &viscousType, PHINT, 1);
        if (viscousType > 1)
        {
            int gridType = geometry->GetGrid()->Type();
            if (gridType == UNSTRUCTGRID)
            {
                //! Add unstructured turbulent solver.
                PHSolver *solver = new TurbSolverUnstr();
                solverList.push_back(solver);
            }
            else if (gridType == STRUCTGRID)
            {
                const int FV_METHOD = 0;
                const int FD_METHOD = 1;
                
                int FVorFD = 0;
                zone->GetData("ifvfd",&FVorFD, PHINT, 1);
                
                if (FVorFD == FV_METHOD)
                {
                    //! Add structured turbulent solver.
                    PHSolver *solver = new TurbSolverStr();
                    solverList.push_back(solver);
                }
                else
                {
                    //! Add structured FD turbulent solver.
                    PHSolver *solver = new TurbSolverStrFD();
                    solverList.push_back(solver);
                }
            }
        }
#endif

        //! Finally, push the solvers into the list.
        for (std::size_t iSolver = 0; iSolver < solverList.size(); ++ iSolver)
        {
            PHSolver *solver = solverList[iSolver];

            solver->SetGeometry(geometry);
            solver->SetIndex(static_cast<int> (iSolver));
            zone->AddSolver(solver);
        }

    }// iZone
}

void Simulation::PreSolve()
{
    int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
    if (isUnsteady != 0)
    {
        AleModel::UnsteadyPreSolve();
    }

    int isInIniting = GlobalDataBase::GetIntParaFromDB("isInIniting");
    if (isInIniting == 1)
    {
        int nMGLevel = GlobalDataBase::GetIntParaFromDB("nMGLevel");
        if (nMGLevel > 1)
        {
            //! For MG case, only initialize on coarse grid.
            isInIniting = 0;
            GlobalDataBase::UpdateData("isInIniting", &isInIniting, PHINT, 1);
        }
        else
        {
            //! For single grid case, check if in the state of initialization.
            int outIterStep = PHSPACE::GlobalDataBase::GetIntParaFromDB("outnstep");
            int flowInitStep = PHSPACE::GlobalDataBase::GetIntParaFromDB("flowInitStep");
            if ( (flowInitStep > 0 && outIterStep > flowInitStep) ||
                  flowInitStep == 0 )
            {
                isInIniting = 0;
                GlobalDataBase::UpdateData("isInIniting", &isInIniting, PHINT, 1);

                PrintToWindow("End init flow by first order ...\n");
            }
        }
    }
}

void Simulation::PostSolve()
{
    int intervalStepPlot = GlobalDataBase::GetIntParaFromDB("intervalStepPlot");
    int intervalStepSample = GlobalDataBase::GetIntParaFromDB("intervalStepSample");
    int outIterStep  = GlobalDataBase::GetIntParaFromDB("outnstep");
    if (outIterStep % intervalStepPlot == 0)
    {
        Post_WriteTecplot generateVisualFile;
        generateVisualFile.Run();
    }

    if (outIterStep % intervalStepSample == 0)
    {
        Post_ProbesWrite generateProbesVarFile;
        generateProbesVarFile.Run();
    }
}

void Simulation::MainSolve()
{
    int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");

    if (isUnsteady == 0)
    {
        SolveSteady();
    }
    else
    {
        SolveUnsteady();
    }
}

void Simulation::SolveOneOuterStep()
{
    PreSolve();

    MainSolve();

    PostSolve();
}

void Simulation::SolveSimulation()
{

    WriteLogFile("Start iterating ...");
    bool continueSolve = true;

    int maxSimuStep = PHSPACE::GlobalDataBase::GetIntParaFromDB("maxSimuStep");
    int iOverlap    = PHSPACE::GlobalDataBase::GetIntParaFromDB("iovrlap");
    int sysGridType = PHSPACE::GlobalDataBase::GetIntParaFromDB("sys_gridtype");

    if (iOverlap != 0 && sysGridType == STRUCTGRID)
    {
        GetOversetStructGrid()->FillSolverContainer();
        GetOversetStructGrid()->InitField();
    }

    //! Wait for all processors here!.
    if (iOverlap == 0)
    {
        PrintToWindow("Wait for all processors before iteration ...\n");
        PH_Barrier();
        PrintToWindow("Wait is over,start iteration ...\n");
    }

    int outIterStep = PHSPACE::GlobalDataBase::GetIntParaFromDB("outnstep");
    if (maxSimuStep == 0)
    {
        continueSolve = false;
    }

    while (continueSolve)
    {
        //TODO: To modify the function name
        UpdateOuterIterationStep(outIterStep);

        SolveOneOuterStep();

        if (outIterStep >= maxSimuStep && maxSimuStep > 0)
        {
            continueSolve = false;
        }
    }
}

void Simulation::PostSimulation()
{
	this->region->PostSimulation();

    int iOverlap = PHSPACE::GlobalDataBase::GetIntParaFromDB("iovrlap");

    int sysGridType = PHSPACE::GlobalDataBase::GetIntParaFromDB("sys_gridtype");

    if (iOverlap != 0 && sysGridType == STRUCTGRID)
    {
        PHSPACE::DeleteOversetStructGrid();

        PHSPACE::DeleteStructGridTable();
    }

    return;
}

void Simulation::UpdateOuterIterationStep(int &outIterStep)
{
    outIterStep = outIterStep + 1;
    GlobalDataBase::UpdateData("outnstep", &outIterStep, PHINT, 1);

    int newIterStep = GlobalDataBase::GetIntParaFromDB("newnstep");
    newIterStep ++;
    GlobalDataBase::UpdateData("newnstep", &newIterStep, PHINT, 1);

    int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
    if (isUnsteady == 0)
    {
        return;
    }

    int ifStaticsFlowField = GlobalDataBase::GetIntParaFromDB("ifStaticsFlowField");
    if (ifStaticsFlowField == 1)
    {
        int startStatisticStep = GlobalDataBase::GetIntParaFromDB("startStatisticStep");
        if (outIterStep >= startStatisticStep)
        {
            int numberOfstatisticsDESInterval = GlobalDataBase::GetIntParaFromDB("nStatisticalStep");
            ++ numberOfstatisticsDESInterval;
            GlobalDataBase::UpdateData("nStatisticalStep", &numberOfstatisticsDESInterval, PHINT, 1);
        }
    }
}

void Simulation::UpdateAllZonesPara()
{
    this->region->UpdateAllZonesPara();
}

void Simulation::InitGlobalBoundaryCondition()
{
    this->region->ReSetBoundaryConditionByGlobalBC();
}

void Simulation::MultiGridInitFlow()
{
    WriteLogFile("Multi-Grid Init Flow ...");
    int nMGLevel = GlobalDataBase::GetIntParaFromDB("nMGLevel");
    if (nMGLevel <= 1)
    {
        return;
    }

    int outIterStep = PHSPACE::GlobalDataBase::GetIntParaFromDB("outnstep");
    if (outIterStep > 0)
    {
//      int newnstep = outnstep;
//      GlobalDataBase::UpdateData("newnstep",&newnstep, PHINT, 1);
        return; //! Read from flow file.
    }

    int flowInitStep = PHSPACE::GlobalDataBase::GetIntParaFromDB("flowInitStep");
    if (flowInitStep <= 0) return;

    int nSolver = 1;        //! Consider NS solver only.

    int isInMGIniting = 1;
    GlobalDataBase::UpdateData("isInMGIniting", &isInMGIniting, PHINT, 1);

    for (int iLevel = nMGLevel - 1; iLevel > 0; -- iLevel)
    {
        //TODO: to modify the levelStep for alterative using
        int levelStep = 0;
        int totalLevelStep = flowInitStep * iLevel;

        ostringstream oss;
        oss << "Start init flow on level " << iLevel << ", " << totalLevelStep << " iterations ...";
        WriteLogFile(oss.str()); oss << endl;
        PrintToWindow(oss.str());

        while (levelStep < totalLevelStep)
        {
            ++ levelStep;

            UpdateOuterIterationStep(outIterStep);

            //using namespace PHMPI;
            for (int isolver = 0; isolver < nSolver; ++ isolver)
            {
                Controller * controller = new Controller(isolver, this->region);
                controller->MultiGridInitFlow(iLevel);
                controller->PostSolve(iLevel);
                delete controller;
            }
        }

        //! Interpolate from coarse grid to fine grid.
        for (int isolver = 0; isolver < nSolver; ++ isolver)
        {
            Controller *controller = new Controller(isolver, this->region);
            controller->InterpolatFineGrid(isolver, iLevel - 1);
            delete controller;
        }
    }

    isInMGIniting = 0;
    GlobalDataBase::UpdateData("isInMGIniting", &isInMGIniting, PHINT, 1);

    int newIterStep = 0;
    GlobalDataBase::UpdateData("newnstep",&newIterStep, PHINT, 1);

    PrintToWindow("End init flow by coarse grid ...\n");
    WriteLogFile("End init flow by coarse grid ...");
}

void Simulation::SolveSteady()
{
    int nSolver = Zone::GetNumberOfSolvers();
    // Get iteration number
    int outerIterStep = GlobalDataBase::GetIntParaFromDB("outnstep");
    int turbInterval = GlobalDataBase::GetIntParaFromDB("turbInterval");

    for (int iSolver = 0; iSolver < nSolver; ++ iSolver)
    {
        // isolver=1 or 2, is supposed as turbulence computing
        if ((iSolver != 0) && (outerIterStep % turbInterval != 0))
        {
            return;
        }

        Controller *controller = new Controller(iSolver, this->region);

        controller->SolveSteadyField();

        controller->PostSolve();

        delete controller;
    }
}

void Simulation::SolveUnsteady()
{
    int     minSubItertation    = GlobalDataBase::GetIntParaFromDB("min_sub_iter");
    int     maxSubItertation    = GlobalDataBase::GetIntParaFromDB("max_sub_iter");
    RDouble toltalSubItertation = GlobalDataBase::GetDoubleParaFromDB("tol_sub_iter");

    int innerIterStep = 0;
    int nSolver = Zone::GetNumberOfSolvers();

    bool isSubIterationDump = false;
    GlobalDataBase::UpdateData("isSubIterationDump", &isSubIterationDump, PHBOOL, 1);

    bool isSubIterationConvergence = false;
    while (!isSubIterationConvergence)
    {
        innerIterStep += 1;
        GlobalDataBase::UpdateData("innstep", &innerIterStep, PHINT, 1);

        RDouble subIterationNormOneSolver     =  zero;
        RDouble maxSubIterationNormAllSolvers = -LARGE;

        for (int iSolver = 0; iSolver < nSolver; ++ iSolver)
        {
            Controller *controller = new Controller(iSolver, this->region);

            controller->SolveSteadyField();

            //! Here is in the sub-iteration, PostSolve is only for the CommunicationInterfaceData and dumps the "res.dat", which contains index of the sub-iteration residual.
            controller->PostSolve();

            subIterationNormOneSolver     = controller->ComputeSubIterationNorm(iSolver, 0);
            maxSubIterationNormAllSolvers = MAX(maxSubIterationNormAllSolvers, subIterationNormOneSolver);

            delete controller;
        }

        PH_CompareMaxMin(maxSubIterationNormAllSolvers, 1);

        //! Sub-iteration continue or not judgement.
        isSubIterationConvergence = IsSubIterationConvergence(innerIterStep, minSubItertation, maxSubItertation, toltalSubItertation, maxSubIterationNormAllSolvers);
    }

    //Here is already (innstep == max_sub_iter || isSubIterationConvergence).
    isSubIterationDump = true;
    GlobalDataBase::UpdateData("isSubIterationDump", &isSubIterationDump, PHBOOL, 1);

    //! Update physical time step flowfield after sub-iteration convergence.
    for (int iSolver = 0; iSolver < nSolver; ++ iSolver)
    {
        Controller *controller = new Controller(iSolver, this->region);

        controller->UpdateUnsteadyFlow(iSolver, 0);

        //! Here is out of the sub-iteration, PostSolve is to dump various files according to the parameters in the "cfd_para.hypara", such as intervalStepFlow, intervalStepPlot, intervalStepForce, intervalStepRes, et al.
        controller->PostSolve();

        delete controller;
    }
}

bool Simulation::IsSubIterationConvergence(const int &innerIterStep, const int &minSubIteration, const int &maxSubIteration, const RDouble &tolalSubIteration, const RDouble &averageQuotientMaxInSolvers)
{
    bool isSubIterationConvergence = false;

    if (innerIterStep < minSubIteration)
    {
        isSubIterationConvergence = false;
    }
    else if (averageQuotientMaxInSolvers <= tolalSubIteration)
    {
        isSubIterationConvergence = true;
    }
    else if (innerIterStep >= maxSubIteration)
    {
        isSubIterationConvergence = true;
    }
    else
    {
        isSubIterationConvergence = false;
    }

    return isSubIterationConvergence;
}

}

