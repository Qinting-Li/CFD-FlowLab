#include "Simulation.h"
#include "IO_FileName.h"
#include "TK_Parse.h"
#include "GridFactory.h"
#include "AleManager.h"
#include "GridPartition.h"
#include "Controller.h"
#include "TK_Time.h"
#include "PHIO.h"
#include "Geo_Sample.h"
#include "UnstructuredOversetConfig.h"
#include "HOSimulation.h"
#include "OversetInformation.h"
#include "SolversList.h"
#include "PostProcess.h"
#include "Glb_Dimension.h"
#include "Zone.h"
#include "WallDist.h"
#include "Post_WriteTecplot.h"
#include "LESSolverStruct.h"
#include "LESSolverUnstr.h"
#ifdef USE_SpecDiffHybSolver
#include "SpecDiffHybSolver.h"
#include "SpecSolver.h"
#endif
#include "Post_Identify.h"
#include "Post_ProbesWrite.h"
#ifdef USE_DEMOSOLVER
#include "DemoSolver.h"
#endif
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
    case CREATE_GRID:
        PHSPACE::CreateGrid();
        break;
    case CAL_WALL_DIST:
        ComputeWallDistance();
        break;
    case PARTITION_GRID:
        IsNeedConvert();
        PHSPACE::PartitionGrid();
        break;
    case OVERSETGRID_VIEW:
        PHSPACE::OversetGridView();
        break;
    case OVERSET_CONFIG:
        PHSPACE::RunOversetGridConfig();
        break;
    case CONVERT_PARTITION:
        //ConvertAndPartition();
        break;
    case HO_SOLVER:
        {
            HOUnstruct::HOSimulation * uho = new HOUnstruct::HOSimulation();
            uho->Run();
            delete uho;
            break;
        }
    case INTEGRATIVE_SOLVER:
        RunIntegrativeSolver();
        break;
    case POST_PROCESSING:
        RunPostProcessing();
        break;
    #ifdef USE_SpecDiffHybSolver
    case SPECDIFFHYB_SOLVER:
        RunSpecDiffHybSolver();
        break;
    case SPEC_SOLVER:
        RunSpecSolver();
        break;
    #endif
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

        int iadapt = 0;
        GlobalDataBase::UpdateData("iadapt", &iadapt, PHINT, 1);

        CompleteParamInfo();

        PHSPACE::CreateGrid();

        GlobalDataBase::UpdateData("original_grid_file", &outputGridFile, PHSTRING, 1);

        numberOfSimulationTask = PHSPACE::PARTITION_GRID;
        GlobalDataBase::UpdateData("nsimutask", &numberOfSimulationTask, PHINT, 1);
    }
}

void Simulation::RunCFD()
{
    RDouble last_time[3] = {0.0};
    TIME_SPACE::StartTime(last_time);

    ReadGeometry();

    UpdateAllZonesPara();

    InitGeometry();

    SwapGeommetry();

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
    int tasktype = 0;

    string fromGridName = GlobalDataBase::GetStrParaFromDB("from_gfile");
    string outGridName = ChangeExtensionOfFileName(fromGridName, "fts");
    string partition_grid_file = AddSymbolToFileName(outGridName, "__", numberOfProcessor);
    GlobalDataBase::UpdateData("out_gfile", &outGridName, PHSTRING, 1);

    bool partitionFtsExist = false;
    if (numberOfProcessor > 1)
    {
        if (SerialFileExist(partition_grid_file))
        {
            partitionFtsExist = true;
        }
    }

    //! fts file not exist, need convert grid.
    if (!SerialFileExist(outGridName))
    {
        if (!partitionFtsExist)
        {
            WriteLogFile("Begin convert grid !\n");
            PrintToWindow("Begin convert grid !\n\n");

            tasktype = PHSPACE::CREATE_GRID;
            GlobalDataBase::UpdateData("nsimutask", &tasktype, PHINT, 1);

            int isConvertOver = 0;
            PHMPI::SetNumberOfProcessor(1);

            if (currentProcessor == serverProcessor)
            {
                CreateGrid();
                isConvertOver = 1;
            }

            PHMPI::SetNumberOfProcessor(numberOfProcessor);
            PHMPI::PH_Bcast(&isConvertOver, 1*sizeof(int), serverProcessor);
            WriteLogFile("Convert grid complete !\n");
            PrintToWindow("Convert grid complete !\n\n");
        }
    }

    //! If need PartitionGrid.
    if (numberOfProcessor > 1)
    {
        if (!partitionFtsExist)
        {
            WriteLogFile("Begin partition grid !\n");
            PrintToWindow("Begin partition grid !\n\n");
            PHMPI::SetNumberOfProcessor(1);

            tasktype = PHSPACE::PARTITION_GRID;
            GlobalDataBase::UpdateData("nsimutask", &tasktype, PHINT, 1);

            int numberOfMultigrid = GlobalDataBase::GetIntParaFromDB("nMGLevel");
            int maxproc = numberOfProcessor;
            int pgridtype = GlobalDataBase::GetIntParaFromDB("gridtype");

            GlobalDataBase::UpdateData("pgridtype", &pgridtype, PHINT, 1);
            GlobalDataBase::UpdateData("maxproc", &maxproc, PHINT, 1);
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
            PHMPI::PH_Bcast(&isPartitionOver, 1 * sizeof(int), serverProcessor);
            WriteLogFile("Partition grid complete !\n");
            PrintToWindow("Partition grid complete !\n\n");
        }

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

    WriteLogFile("Begin CFD simulation !\n");
    PrintToWindow("Begin CFD simulation !\n\n");

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


#ifdef USE_SpecDiffHybSolver
void Simulation::RunSpecDiffHybSolver()
{
    double last_time[3];
    TIME_SPACE::StartTime(last_time);

    SpecDiffHybSolver * solver = new SpecDiffHybSolver();

    solver->Run();

    delete solver;

    WriteLogFile("SpecDiffHybSolver Finished!\n");
}

void Simulation::RunSpecSolver()
{
    double last_time[3];
    TIME_SPACE::StartTime(last_time);

    SpecSolver * solver = new SpecSolver();

    solver->Run();

    delete solver;
}
#endif

void Simulation::ComputeWallDistance()
{
    RDouble last_time[3] = {0.0};
    TIME_SPACE::StartTime(last_time);

    ReadGeometry();

    UpdateAllZonesPara();

    InitGeometry();

    GlobalBoundaryCondition::ReadGlobalBoundaryCondition();

    GlobalBoundaryCondition::ChangeBCTypeByGlobalBC();

    InitGlobalBoundaryCondition();

    GenerateWallDist();

    WriteLogFile("Generate wall distance ...");
}

void Simulation::InitOversetGrid()
{
    //! The overset config process must be taken after reading in ordinary grid and initializing!
    //! Read overset information files or config.

    int codeOfOversetGrid = PHSPACE::GlobalDataBase::GetIntParaFromDB("codeOfOversetGrid");

    if (!codeOfOversetGrid)
    {
        return;
    }

    int readOversetFileOrNot = PHSPACE::GlobalDataBase::GetIntParaFromDB("readOversetFileOrNot");

    //! If the processor is grid processor, and overset grid need to be assembled.
    if (PHMPI::CurrentProcessorIsGridProcessor() && readOversetFileOrNot != 1)
    {
        OversetGridConfig();
        WriteLogFile("Run overset configuration ...");
    }

    //! Wait for the grid processors.
    PH_Barrier();

    //! For the original processors.
    if (!PHMPI::CurrentProcessorIsGridProcessor())
    {
        //! If the overset files is aviable or has been prepared by grid processors, just read it!
        if (readOversetFileOrNot == 1 || PHMPI::GetNumberOfGridProcessor() != 0)
        {
            region->ReadOversetGrid();
            WriteLogFile("Read overset information \n");
        }
        //! Else it should be done by these processors by themselves.
        else
        {
            OversetGridConfig();
            WriteLogFile("Run overset configuration \n");
        }

        InitializeOversetInterfaceTopology();

        WriteLogFile("Initialize oversetInterface topology \n");
    }
}

void Simulation::OversetGridConfig()
{
    OversetConfigFactoryMulti *oversetConfigFactory = new OversetConfigFactoryMulti();

    oversetConfigFactory->Initialize();

    oversetConfigFactory->OversetConfigMulti();

    oversetConfigFactory->OutPutOversetInformation();

    delete oversetConfigFactory;
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
    WriteLogFile("Creating and initialization solvers start \n");

    CreateAleManager();

    int nSolver = Zone::GetNumberOfSolvers();
    for (int iSolver = 0; iSolver < nSolver; ++ iSolver)
    {
        Controller *controller = new Controller(iSolver, this->region);

        controller->Initialize(iSolver);

        delete controller;
    }

    WriteLogFile("Creating and initialization solvers over \n");
}

void Simulation::DefineAndCreateSolvers()
{
    InitGlobalValuesOfAllZones();

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
        PHGeometry *geometry = zone->GetGeometry();

#ifdef USE_DEMOSOLVER
        PHSolver *solver = new DemoSolver();
        solverList.push_back(solver);
#else  

        //! For each zone, define the solver names on it.
        //! If you are solver developer, define your solve name here,
        //! according your solver logic.
        int dgHighOrder = 0;
        zone->GetData("dg_high_order", &dgHighOrder, PHINT, 1);

        //! First: NS equations solver.
        if (!dgHighOrder)
        {
            int gridType = geometry->GetGrid()->Type();

            int compressible = 1;
            zone->GetData("compressible", &compressible, PHINT, 1);

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
                    int FVorFD = 0;
                    zone->GetData("isFVMOrFDM", &FVorFD, PHINT, 1);

                    int iCodeBranch = 0;
                    zone->GetData("iCodeBranch", &iCodeBranch, PHINT, 1);

                    if (FVorFD == FV_METHOD) //! Add structured FV NS solver.
                    {
                        if (iCodeBranch == 0)
                        {
                            PHSolver *solver = new NSSolverStructFV();
                            solverList.push_back(solver);
                        }
                        else //! A new code branch is active for Navier-Stokes solver.
                        {
                        }
                    }
                    else //! Add structured FD NS solver.
                    {
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
        if (viscousType > LAMINAR)
        {
            int compressible = 1;
            zone->GetData("compressible", &compressible, PHINT, 1);
            if (compressible)
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
                    const int NS_METHOD = 2;

                    int FVorFD = 0;
                    zone->GetData("isFVMOrFDM", &FVorFD, PHINT, 1);

                    if (FVorFD == FV_METHOD || FVorFD == NS_METHOD)
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
        }

        int iLES = 0;
        zone->GetData("iLES", &iLES, PHINT, 1);
        if (iLES == LES_SOLVER)
        {
            int gridType = geometry->GetGrid()->Type();
            if (gridType == UNSTRUCTGRID)
            {
                //! Add unstructured LES solver.
                PHSolver *solver = new LESSolverUnstr();
                solverList.push_back(solver);
            }
            else if (gridType == STRUCTGRID)
            {
                //! Add structured LES solver.
                PHSolver *solver = new LESSolverStruct();
                solverList.push_back(solver);
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
    } //! iZone
}

void Simulation::PreSolve()
{
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
            if ((flowInitStep > 0 && outIterStep > flowInitStep) || flowInitStep == 0)
            {
                isInIniting = 0;
                GlobalDataBase::UpdateData("isInIniting", &isInIniting, PHINT, 1);

                PrintToWindow("End init flow by first order ... \n");
            }
        }
    }
}

void Simulation::PostSolve()
{
    using namespace PHMPI;

    //! Dump tecflow file.
    int intervalStepPlot = GlobalDataBase::GetIntParaFromDB("intervalStepPlot");
    int intervalStepSample = GlobalDataBase::GetIntParaFromDB("intervalStepSample");
    int outIterStep = GlobalDataBase::GetIntParaFromDB("outnstep");
    if (outIterStep % intervalStepPlot == 0)
    {
        Post_WriteTecplot generateVisualFile;
        generateVisualFile.Run();
    }

    if (!GlobalDataBase::IsExist("ifSetDataMonitor", PHINT, 1))
    {
        return;
    }
    else
    {
        int ifSetDataMonitor = GlobalDataBase::GetIntParaFromDB("ifSetDataMonitor");
        if (!ifSetDataMonitor) return;

        if (outIterStep % intervalStepSample == 0)
        {
            Post_ProbesWrite generateProbesVarFile;
            generateProbesVarFile.Run();
        }
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
    int codeOfOversetGrid = PHSPACE::GlobalDataBase::GetIntParaFromDB("codeOfOversetGrid");
    int sysGridType = PHSPACE::GlobalDataBase::GetIntParaFromDB("sys_gridtype");

    //! Wait for all processors here!
    if (codeOfOversetGrid == 0 || sysGridType == STRUCTGRID)
    {
        PrintToWindow("Wait for all processors before iteration ... \n");
        PH_Barrier();
        PrintToWindow("Wait is over,start iteration ... \n");
    }

    int outIterStep = PHSPACE::GlobalDataBase::GetIntParaFromDB("outnstep");
    if (maxSimuStep == 0)
    {
        continueSolve = false;
    }

    while (continueSolve)
    {
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

    int codeOfOversetGrid = PHSPACE::GlobalDataBase::GetIntParaFromDB("codeOfOversetGrid");

    int sysGridType = PHSPACE::GlobalDataBase::GetIntParaFromDB("sys_gridtype");

    int IsUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
    int IsAleSimulation = GlobalDataBase::GetIntParaFromDB("codeOfAleModel");
    if (IsUnsteady != 0  && IsAleSimulation !=0)
    {
        FreeAleManager();
    }

    return;
}

void Simulation::UpdateOuterIterationStep(int &outIterStep)
{
    outIterStep ++;
    GlobalDataBase::UpdateData("outnstep", &outIterStep, PHINT, 1);

    int newIterStep = GlobalDataBase::GetIntParaFromDB("newnstep");
    newIterStep ++;
    GlobalDataBase::UpdateData("newnstep", &newIterStep, PHINT, 1);

    int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
    if (isUnsteady == 0)
    {
        return;
    }

    RDouble physicalTime = GlobalDataBase::GetDoubleParaFromDB("physicalTime");
    RDouble physicalTimeStep = GlobalDataBase::GetDoubleParaFromDB("physicalTimeStep");

    physicalTime += physicalTimeStep;

    GlobalDataBase::UpdateData("physicalTime", &physicalTime, PHDOUBLE, 1);

    int ifStaticsFlowField = GlobalDataBase::GetIntParaFromDB("ifStaticsFlowField");
    if (ifStaticsFlowField == 1)
    {
        int startStatisticStep = GlobalDataBase::GetIntParaFromDB("startStatisticStep");
        if (outIterStep >= startStatisticStep)
        {
            int nStatisticalStep = GlobalDataBase::GetIntParaFromDB("nStatisticalStep");
            ++ nStatisticalStep;
            GlobalDataBase::UpdateData("nStatisticalStep", &nStatisticalStep, PHINT, 1);
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
        return; //! Read from flow file.
    }

    int flowInitStep = PHSPACE::GlobalDataBase::GetIntParaFromDB("flowInitStep");
    if (flowInitStep <= 0) return;

    int nSolver = 1;    //! Consider NS solver only.

    int isInMGIniting = 1;
    GlobalDataBase::UpdateData("isInMGIniting", &isInMGIniting, PHINT, 1);

    for (int iLevel = nMGLevel - 1; iLevel > 0; -- iLevel)
    {
        int levelStep = 0;
        int totalLevelStep = flowInitStep * iLevel;

        ostringstream oss;
        oss << "Start init flow on level " << iLevel << ", " << totalLevelStep << " iterations ...";
        WriteLogFile(oss.str());
        oss << endl;
        PrintToWindow(oss.str());

        while (levelStep < totalLevelStep)
        {
            ++ levelStep;

            UpdateOuterIterationStep(outIterStep);

            //using namespace PHMPI;
            for (int isolver = 0; isolver < nSolver; ++ isolver)
            {
                Controller *controller = new Controller(isolver, this->region);
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
    GlobalDataBase::UpdateData("newnstep", &newIterStep, PHINT, 1);

    PrintToWindow("End init flow by coarse grid ... \n");
    WriteLogFile("End init flow by coarse grid ...");
}

void Simulation::SolveSteady()
{
    int nSolver = Zone::GetNumberOfSolvers();
    //! Get iteration number.
    int outerIterStep = GlobalDataBase::GetIntParaFromDB("outnstep");
    int turbInterval = GlobalDataBase::GetIntParaFromDB("turbInterval");

    for (int iSolver = 0; iSolver < nSolver; ++ iSolver)
    {
        //! isolver = 1 or 2, is supposed as turbulence computing.
        if ((iSolver == TURB_SOLVER) && (outerIterStep % turbInterval != 0))
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
    int minSubItertation = GlobalDataBase::GetIntParaFromDB("min_sub_iter");
    int maxSubItertation = GlobalDataBase::GetIntParaFromDB("max_sub_iter");
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

        SolveSingleInnerIterationStepForAleSolvers();

        RDouble subIterationNormOneSolver = zero;
        RDouble maxSubIterationNormAllSolvers = -LARGE;

        for (int iSolver = 0; iSolver < nSolver; ++ iSolver)
        {
            Controller *controller = new Controller(iSolver, this->region);

            controller->SolveSteadyField();

            //! Here is in the sub-iteration, PostSolve is only for the CommunicationInterfaceData and dumps the "res.dat", which contains index of the sub-iteration residual.
            controller->PostSolve();

            subIterationNormOneSolver = controller->ComputeSubIterationNorm(iSolver, 0);
            maxSubIterationNormAllSolvers = MAX(maxSubIterationNormAllSolvers, subIterationNormOneSolver);

            delete controller;
        }

        PH_CompareMaxMin(maxSubIterationNormAllSolvers, 1);

        //! Sub-iteration continue or not judgement.
        isSubIterationConvergence = IsSubIterationConvergence(innerIterStep, minSubItertation, maxSubItertation, toltalSubItertation, maxSubIterationNormAllSolvers);
    }

    PostprocessAfterInnerIterationStepForAleSolvers();

    //! Here is already (innstep == max_sub_iter || isSubIterationConvergence).
    isSubIterationDump = true;
    GlobalDataBase::UpdateData("isSubIterationDump", &isSubIterationDump, PHBOOL, 1);

    //! Update physical time step flow field after sub-iteration convergence.
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

