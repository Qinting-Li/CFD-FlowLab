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
#ifdef USE_INHOUSE_VERSION
#include "Post_Identify.h"
#endif
#ifdef USE_LBM_CODE
#include "cavity.h"
#endif

#include "OversetInformation.h"
#include "SolversList.h"
#include "PostProcess.h"
#include "Glb_Dimension.h"
#include "Zone.h"
#include "Pre_WalldistCompute.h"
#include "Post_WriteTecplot.h"
#include "LESSolverStruct.h"
#include "LESSolverUnstr.h"
#ifdef USE_SpecDiffHybSolver
#include "SpecDiffHybSolver.h"
#include "SpecSolver.h"
#endif

#include "Post_ProbesWrite.h"
#include "HyFlowSolverStruct.h"
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
    #ifdef USE_INHOUSE_VERSION
    case DSMC_SOLVER:
        RunDSMC();
        break;
    case HELICOPTER:
        RunHelicopterSolver();
        break;
    case UGKS_SOLVER:
        RunUGKS();
        break;
    case DERIVATIVE_IDENITFY:
        IdentifyDerivative();
        break;
    #endif
    #ifdef USE_LBM_CODE
    case LBM_SOLVER:
        RunLBMSolver();
        break;
    #endif
    case INTEGRATIVE_SOLVER:
        RunIntegrativeSolver();
        break;
    #ifdef USE_INHOUSE_VERSION
    case GRIDCHECK_SOLVER:
        RunGridCheckSolver();
        break;
    #endif
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
    gridTypeMap.insert(pair <string, int> ("grd",  PLOT3D_TYPE   ));
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

    int iCodeBranch = GlobalDataBase::GetIntParaFromDB("iCodeBranch");
    if (iCodeBranch > 0)
    {
        GenerateWallDist();

        SelfAdaptionRunCFD(); //run CFD using the adaptive method.
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

bool Simulation::IsRestartProcess()
{
    using namespace PHMPI;

    string restartFile = ".\results\flow.dat";
    GlobalDataBase::GetData("restartNSFile", &restartFile, PHSTRING, 1);
    if(PHMPI::IsParallelRun())
    {
        restartFile = PHSPACE::AddSymbolToFileName(restartFile, "_", 0);
    }

    ifstream infile(restartFile.c_str(), ios::in);
    bool isRestart = false;
    if (infile)
    {
        isRestart = true;
        infile.close();
        infile.clear();
    }

    if (isRestart)
    {
        hid_t file;
        file = OpenHDF5File(restartFile);
        int gridID = 0;
        hid_t grploc;
        string grpName;
        ostringstream oss;
        oss<<"Group"<<gridID;
        grpName = oss.str();

        grploc = OpenGroup(file, grpName);

        int preChemical = 0, preTemperatureModel = 1;
        if (CheckDataExist(grploc, "nIsChemicalFlow"))
        {
            ReadData(grploc, &preChemical, "nIsChemicalFlow");
        }
        if (CheckDataExist(grploc, "nTemperatureModel"))
        {
            ReadData(grploc, &preTemperatureModel, "nTemperatureModel");
        }

        GlobalDataBase::UpdateData("previousChemical", &preChemical, PHINT, 1);
        GlobalDataBase::UpdateData("previousTModel", &preTemperatureModel, PHINT, 1);

        H5Gclose(grploc);
        H5Fclose(file);
    }

    return isRestart;
}

void Simulation::SelfAdaptionRunCFD()
{
    int nFinalChemical = GlobalDataBase::GetIntParaFromDB("nchem");
    GlobalDataBase::UpdateData("nFinalChemical", &nFinalChemical, PHINT, 1);
    int nFinalTModel = GlobalDataBase::GetIntParaFromDB("ntmodel");
    GlobalDataBase::UpdateData("nFinalTModel", &nFinalTModel, PHINT, 1);
    string nFinalNameList = GlobalDataBase::GetStrParaFromDB("speciesName");
    string nFinalMassFraction = GlobalDataBase::GetStrParaFromDB("initMassFraction");
    GlobalDataBase::UpdateData("nFinalNameList", &nFinalNameList, PHSTRING, 1);
    GlobalDataBase::UpdateData("nFinalMassFraction", &nFinalMassFraction, PHSTRING, 1);

    if (IsRestartProcess())
    {
        int preChemical = GlobalDataBase::GetIntParaFromDB("previousChemical");
        int preTModel = GlobalDataBase::GetIntParaFromDB("previousTModel");

        if (preChemical == 1 && preTModel == 1)
        {
            //! restart from the one-temperature model.
            SecondStepRunCFD(true);
            if (nFinalTModel == 1)
            {
                //The third step: the final gas model.
                FinalStepRunCFD();
            }
            else
            {
                //The second step: two-temperature model.
                ThirdStepRunCFD();

                //The third step: the final gas model.
                FinalStepRunCFD();
            }

            return;
        }
        else if (preChemical == 1 && preTModel == 2)
        {
            //! restart from the two-temperature model.
            ThirdStepRunCFD(true);

            //The third step: the final gas model.
            FinalStepRunCFD();

            return;
        }
    }

    //The first step: the perfect gas.
    FirstStepRunCFD();
    if (nFinalChemical == 0)
    {
        return;
    }
    //The second step: one-temperature model.
    SecondStepRunCFD();
    if (nFinalTModel == 1)
    {
        //The third step: the final gas model.
        FinalStepRunCFD();
    }
    else
    {
        //The second step: two-temperature model.
        ThirdStepRunCFD();

        //The third step: the final gas model.
        FinalStepRunCFD();
    }
}

void Simulation::FirstStepRunCFD(bool isStart/* = true*/)
{
    //Record the serial number of current iteration step.
    int nCurrentStep = 1;
    GlobalDataBase::UpdateData("nCurrentStep", &nCurrentStep, PHINT, 1);
    int nCurrentProcessorID = PHMPI::GetCurrentProcessorID();
    if (nCurrentProcessorID == PHMPI::server)
    {
        cout<<endl;
        cout<<"nCurrentStep="<<nCurrentStep<<endl;
    }

    //The first step: run CFD using the perfect gas model.
    int nCurChemical = 0;
    GlobalDataBase::UpdateData("nchem", &nCurChemical, PHINT, 1);
    int nCurTModel = 1;
    GlobalDataBase::UpdateData("ntmodel", &nCurTModel, PHINT, 1);

    //The index of the variable pressure.
    int nKeyVariableIndex = 4;
    GlobalDataBase::UpdateData("nKeyVariableIndex", &nKeyVariableIndex, PHINT, 1);
    string keyVariable = "pressure";
    GlobalDataBase::UpdateData("keyVariable", &keyVariable, PHSTRING, 1);

    DefineAndCreateSolvers();

    InitGlobalBoundaryCondition();

    GenerateWallDist();

    InitializeSolvers();

    MultiGridInitFlow();

    SelfAdaptionSolve();

    PostSimulation();

    //Release the dynamic memories at first.
    CleanGlobalValuesSolver();

    RefreshSolvers();
}

void Simulation::SecondStepRunCFD(bool isStart/* = false*/)
{
    //Record the serial number of current iteration step.
    int nCurrentStep = 2;
    GlobalDataBase::UpdateData("nCurrentStep", &nCurrentStep, PHINT, 1);
    int nCurrentProcessorID = PHMPI::GetCurrentProcessorID();
    if (nCurrentProcessorID == PHMPI::server)
    {
        cout<<endl;
        cout<<"nCurrentStep="<<nCurrentStep<<endl;
    }

    //The second step: run CFD using the one-temperature model.
    int nCurChemical = 1;
    GlobalDataBase::UpdateData("nchem", &nCurChemical, PHINT, 1);
    int nCurTModel = 1;
    GlobalDataBase::UpdateData("ntmodel", &nCurTModel, PHINT, 1);
    string speciesName = "O, O2, NO, N, N2";
    GlobalDataBase::UpdateData("speciesName", &speciesName, PHSTRING, 1);
    string initMassFraction = "0.0, 0.233, 0.0, 0.0, 0.767";
    GlobalDataBase::UpdateData("initMassFraction", &initMassFraction, PHSTRING, 1);

    int nControlVariable = 1;
    if (GlobalDataBase::IsExist("nControlVariable",PHINT, 1))
    {
        nControlVariable = GlobalDataBase::GetIntParaFromDB("nControlVariable");
    }
    int nKeyVariableIndex = nControlVariable;
    if (nControlVariable == 5)
    {
        nKeyVariableIndex = 6;
    }
    else if (nControlVariable == 6)
    {
        nKeyVariableIndex = 9;
    }
    else if (nControlVariable > 6) //default
    {
        nKeyVariableIndex = 1; //The index of the variable temperature.
    }
    GlobalDataBase::UpdateData("nKeyVariableIndex", &nKeyVariableIndex, PHINT, 1);

    string keyVariable = "";
    switch (nControlVariable)
    {
    case 0:
        keyVariable = "density";
    case 1:
        keyVariable = "translation temperature";
    case 2:
        keyVariable = "vibration temperature";
    case 3:
        keyVariable = "electron temperature";
    case 4:
        keyVariable = "pressure";
    default:
        keyVariable = "translation temperature";
        break;
    }
    GlobalDataBase::UpdateData("keyVariable", &keyVariable, PHSTRING, 1);

    if (!isStart)
    {
        InitGlobalValuesSolver();
    }
    else
    {
        DefineAndCreateSolvers();

        GenerateWallDist();
    }

    InitGlobalBoundaryCondition();

    InitializeSolvers();

    MultiGridInitFlow();

    SelfAdaptionSolve();

    PostSimulation();

    //Release the dynamic memories at first.
    CleanGlobalValuesSolver();

    RefreshSolvers();
}

void Simulation::ThirdStepRunCFD(bool isStart/* = false*/)
{
    //Record the serial number of current iteration step.
    int nCurrentStep = 3;
    GlobalDataBase::UpdateData("nCurrentStep", &nCurrentStep, PHINT, 1);
    int nCurrentProcessorID = PHMPI::GetCurrentProcessorID();
    if (nCurrentProcessorID == PHMPI::server)
    {
        cout<<endl;
        cout<<"nCurrentStep="<<nCurrentStep<<endl;
    }

    //The third step: run CFD using the two-temperature model.
    int nCurChemical = 1;
    GlobalDataBase::UpdateData("nchem", &nCurChemical, PHINT, 1);
    int nCurTModel = 2;
    GlobalDataBase::UpdateData("ntmodel", &nCurTModel, PHINT, 1);
    int nEnergyAssembly = 0;
    if (GlobalDataBase::IsExist("nEnergyAssembly", PHINT, 1))
    {
        nEnergyAssembly = GlobalDataBase::GetIntParaFromDB("nEnergyAssembly");
    }
    if (nEnergyAssembly == 1)
    {
        //the vibration energy is computed using the polynomial fitting method.
        int nTEnergyModel = 1;
        GlobalDataBase::UpdateData("nTEnergyModel", &nTEnergyModel, PHINT, 1);
    }

    string speciesName = "O, O2, NO, N, N2";
    GlobalDataBase::UpdateData("speciesName", &speciesName, PHSTRING, 1);
    string initMassFraction = "0.0, 0.233, 0.0, 0.0, 0.767";
    GlobalDataBase::UpdateData("initMassFraction", &initMassFraction, PHSTRING, 1);

    //The index of the variable vibrational temperature.
    int nKeyVariableIndex = 2;
    GlobalDataBase::UpdateData("nKeyVariableIndex", &nKeyVariableIndex, PHINT, 1);
    string keyVariable = "vibration temperature";
    GlobalDataBase::UpdateData("keyVariable", &keyVariable, PHSTRING, 1);

    if (!isStart)
    {
        InitGlobalValuesSolver();
    }
    else
    {
        DefineAndCreateSolvers();

        GenerateWallDist();
    }

    InitGlobalBoundaryCondition();

    InitializeSolvers();

    MultiGridInitFlow();

    SelfAdaptionSolve();

    PostSimulation();

    //Release the dynamic memories at first.
    CleanGlobalValuesSolver();

    RefreshSolvers();
}

void Simulation::FinalStepRunCFD()
{
    //Record the serial number of current iteration step.
    int nCurrentStep = 4;
    GlobalDataBase::UpdateData("nCurrentStep", &nCurrentStep, PHINT, 1);
    int nCurrentProcessorID = PHMPI::GetCurrentProcessorID();
    if (nCurrentProcessorID == PHMPI::server)
    {
        cout<<endl;
        cout<<"nCurrentStep="<<nCurrentStep<<endl;
    }

    //The third step: run CFD using the final gas model.
    int nCurChemical = GlobalDataBase::GetIntParaFromDB("nFinalChemical");
    GlobalDataBase::UpdateData("nchem", &nCurChemical, PHINT, 1);
    int nCurTModel = GlobalDataBase::GetIntParaFromDB("nFinalTModel");
    GlobalDataBase::UpdateData("ntmodel", &nCurTModel, PHINT, 1);
    string speciesName = GlobalDataBase::GetStrParaFromDB("nFinalNameList");
    GlobalDataBase::UpdateData("speciesName", &speciesName, PHSTRING, 1);
    string initMassFraction = GlobalDataBase::GetStrParaFromDB("nFinalMassFraction");
    GlobalDataBase::UpdateData("initMassFraction", &initMassFraction, PHSTRING, 1);
    int nEnergyAssembly = 0;
    if (GlobalDataBase::IsExist("nEnergyAssembly", PHINT, 1))
    {
        nEnergyAssembly = GlobalDataBase::GetIntParaFromDB("nEnergyAssembly");
    }
    if (nEnergyAssembly == 1)
    {
        //the vibration energy is computed using the default method.
        int nTEnergyModel = 0;
        GlobalDataBase::UpdateData("nTEnergyModel", &nTEnergyModel, PHINT, 1);
    }

    //default: the index of the variable vibrational temperature.
    int nKeyVariableIndex = 2;
    string keyVariable = "vibration temperature";
    if (nCurTModel == 1)
    {
        int nControlVariable = 1;
        if (GlobalDataBase::IsExist("nControlVariable",PHINT, 1))
        {
            nControlVariable = GlobalDataBase::GetIntParaFromDB("nControlVariable");
        }

        switch (nControlVariable)
        {
        case 0:
            keyVariable = "density";
        case 1:
            keyVariable = "translation temperature";
        case 2:
            keyVariable = "vibration temperature";
        case 3:
            keyVariable = "electron temperature";
        case 4:
            keyVariable = "pressure";
        default:
            keyVariable = "translation temperature";
            break;
        }

        nKeyVariableIndex = nControlVariable; //The index of the variable temperature.
        if (nControlVariable == 5)
        {
            nKeyVariableIndex = 6;
        }
        else if (nControlVariable == 6)
        {
            nKeyVariableIndex = 9;
        }
        else if (nControlVariable > 6)
        {
            nKeyVariableIndex = 1;
        }
    }
    GlobalDataBase::UpdateData("nKeyVariableIndex", &nKeyVariableIndex, PHINT, 1);
    GlobalDataBase::UpdateData("keyVariable", &keyVariable, PHSTRING, 1);

    InitGlobalValuesSolver();

    InitGlobalBoundaryCondition();

    InitializeSolvers();

    MultiGridInitFlow();

    SelfAdaptionSolve();

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

void Simulation::RunGridCheckSolver()
{
    string inputGridName = GlobalDataBase::GetStrParaFromDB("input_gfile");
    GlobalDataBase::UpdateData("gridfile", &inputGridName, PHSTRING, 1);

    Region *region = new Region();
    region->ReadGrid();

    int nZones = PHMPI::GetNumberofGlobalZones();
    Grid **gridArray = new Grid * [nZones];

    fstream file;
    string FantasyCheckFileName = ChangeExtensionOfFileName(inputGridName, "fts_check");
    PHSPACE::OpenSerialFile(file, FantasyCheckFileName, ios_base::out);

    int *globalBlockProcInp = PHMPI::GetZoneProcessorID_INP();
    int *globalBlockType    = PHMPI::GetZoneGridType();
    int *globalBlockIdx     = PHMPI::GetZoneGridID();

    file << "TotalInfo : Global information of current grid." << endl;
    file << "block_idxNumber  ZoneID" << endl;
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        file << iZone << setw(10) << globalBlockIdx[iZone] << endl;
    }

    file << "block_procNumber  ZoneProcessorID" << endl;
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        file << iZone << setw(10) << globalBlockProcInp[iZone] << endl;
    }

    file << "block_typeNumber  ZoneGridType" << endl;
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        file << iZone << setw(10) << globalBlockType[iZone] << endl;
    }

    file << "nBlocks" << setw(10) << nZones << endl;
    file << endl;

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        gridArray[iZone] = GetGrid(iZone);
        int gridType = gridArray[iZone]->Grid::Type();

        if(gridType == UNSTRUCTGRID)
        {
            UnstructGrid *grid = UnstructGridCast(gridArray[iZone]);

            int nTotalNode = grid->GetNTotalNode();
            int nTotalFace = grid->GetNTotalFace();
            int nTotalCell = grid->GetNTotalCell();

            file << "ZoneID : " << iZone << endl;
            file << "GridInformation  Value" << endl;
            file << "nTotalNode" << setw(10) << nTotalNode << endl;
            file << "nTotalFace" << setw(10) << nTotalFace << endl;
            file << "nTotalCell" << setw(10) << nTotalCell << endl;

            file << "CellTopology." << endl;
            int *nodeNumberOfEachCell = grid->GetNodeNumberOfEachCell();
            int *cell2Node = grid->GetCell2Node();
            int nodeCount = 0;

            file << "nodeNumberOfEachCell :" << endl;
            file << "cellNumber nodeNumber " << endl;
            for (int iCell =0; iCell < nTotalCell; iCell++)
            {
                file << iCell << setw(8) << nodeNumberOfEachCell[iCell] << endl;
                nodeCount += nodeNumberOfEachCell[iCell];
            }

            file << "cell2Node :" << endl;
            file << "nodeCountNumber nodeNumber " << endl;
            for (int iCount =0; iCount < nodeCount; iCount++)
            {
                file << iCount << setw(8) << cell2Node[iCount] << endl;
            }

            file << "FaceTopology." << endl;

            int nBoundFace = grid->GetNBoundFace();
            CompositeUnstructBC **bcRecord = grid->GetBCRecord();

            file << "nBoundFace : " << nBoundFace << endl;

            file << "faceNumber  BcType  BCName" << endl;
            for (int iFace = 0; iFace < nBoundFace; iFace++)
            {
                int bcType = bcRecord[iFace]->GetKey();
                string bcName = bcRecord[iFace]->GetBCName();
                file << iFace << setw(10) << bcType << setw(20) << bcName <<endl;
            }

            int *nodeNumberOfEachFace = grid->GetNodeNumberOfEachFace();
            int *face2Node = grid->GetFace2Node();
            int faceCount = 0;

            file << "nodeNumberOfEachFace :" << endl;
            file << "faceNumber nodeNumber " << endl;
            for (int iFace = 0; iFace < nTotalFace; iFace++)
            {
                file << iFace << setw(8) << nodeNumberOfEachFace[iFace] << endl;
                faceCount += nodeNumberOfEachFace[iFace];
            }

            file << "face2Node :" << endl;
            file << "faceCountNumber nodeNumber " << endl;
            for (int iCount =0; iCount < faceCount; iCount++)
            {
                file << iCount << setw(8) << face2Node[iCount] << endl;
            }

            int *leftCellOfFace  = grid->GetLeftCellOfFace();
            int *rightCellOfFace = grid->GetRightCellOfFace();
            
            file << "leftCellOfFace :" << endl;
            file << "faceNumber leftCellNumber " << endl;
            for (int iFace = 0; iFace < nTotalFace; iFace++)
            {
                file << iFace << setw(8) << leftCellOfFace[iFace] << endl;;
            }

            file << "rightCellOfFace :" << endl;
            file << "faceNumber rightCellNumber " << endl;
            for (int iFace = 0; iFace < nTotalFace; iFace++)
            {
                file << iFace << setw(8) << rightCellOfFace[iFace] << endl;;
            }

            InterfaceInfo *interfaceInfo = grid->GetInterfaceInfo();
            if (interfaceInfo)
            {
                int nIFace = interfaceInfo->GetNIFace();
                file << "nIFace = " << nIFace << endl;

                int *interFace2ZoneID       = interfaceInfo->GetInterFace2ZoneID();
                int *interFace2InterFaceID  = interfaceInfo->GetInterFace2InterFaceID();
                int *interFace2BoundaryFace = interfaceInfo->GetInterFace2BoundaryFace();

                file << "interFace2ZoneID" << endl;
                file << "interFaceNumber  ZoneID" << endl;
                for (int iFace = 0; iFace < nIFace; iFace++)
                {
                    file << iFace << setw(8) << interFace2ZoneID[iFace] << endl;
                }

                file << "interFace2InterFaceID" << endl;
                file << "interFaceNumber  InterFaceID" << endl;
                for (int iFace = 0; iFace < nIFace; iFace++)
                {
                    file << iFace << setw(8) << interFace2InterFaceID[iFace] << endl;
                }

                file << "interFace2BoundaryFace" << endl;
                file << "interFaceNumber  BoundaryFaceID" << endl;
                for (int iFace = 0; iFace < nIFace; iFace++)
                {
                    file << iFace << setw(8) << interFace2BoundaryFace[iFace] << endl;
                }
                
            }
            else
            {
                file << "There is no interface exist, nIFace = 0." << endl;
            }

            InterpointInformation *interpointInformation = grid->GetInterpointInfo();
            if (interpointInformation)
            {
                int nIPoint = interpointInformation->GetNumberOfInterpoints();
                file << "nIPoint = " << nIPoint << endl;

                int *cellNumberOfInterPoint = interpointInformation->GetCellNumberOfInterPoint();
                int *interPoint2GlobalPoint = interpointInformation->GetInterPoint2GlobalPoint();
                int *interPoint2InterPointID = interpointInformation->GetInterPoint2InterPointID();
                int *interPoint2ZoneID = interpointInformation->GetInterPoint2ZoneID();
                int *labelOfInterPoint = interpointInformation->GetLabelOfInterPoint();
                int *totalZonesOfInterPoint = interpointInformation->GetTotalZonesOfInterPoint();

                file << "cellNumberOfInterPoint" << endl;
                file << "interPointNumber  totalCellNumbers" << endl;
                for (int iPoint = 0; iPoint < nIPoint; iPoint++)
                {
                    file << iPoint << setw(8) << cellNumberOfInterPoint[iPoint] << endl;
                }

                file << "interPoint2GlobalPoint" << endl;
                file << "interPointNumber  gobalPointNumber" << endl;
                for (int iPoint = 0; iPoint < nIPoint; iPoint++)
                {
                    file << iPoint << setw(8) << interPoint2GlobalPoint[iPoint] << endl;
                }

                file << "interPoint2InterPointID" << endl;
                file << "interPointNumber  interPointID" << endl;
                for (int iPoint = 0; iPoint < nIPoint; iPoint++)
                {
                    file << iPoint << setw(8) << interPoint2InterPointID[iPoint] << endl;
                }

                file << "interPoint2ZoneID" << endl;
                file << "interPointNumber  ZoneID" << endl;
                for (int iPoint = 0; iPoint < nIPoint; iPoint++)
                {
                    file << iPoint << setw(8) << interPoint2ZoneID[iPoint] << endl;
                }

                file << "labelOfInterPoint" << endl;
                file << "interPointNumber  label" << endl;
                for (int iPoint = 0; iPoint < nIPoint; iPoint++)
                {
                    file << iPoint << setw(8) << labelOfInterPoint[iPoint] << endl;
                }

                file << "totalZonesOfInterPoint" << endl;
                file << "interPointNumber  totalZones" << endl;
                for (int iPoint = 0; iPoint < nIPoint; iPoint++)
                {
                    file << iPoint << setw(8) << totalZonesOfInterPoint[iPoint] << endl;
                }
                
            }
            else
            {
                file << "There is no interpoint exist, nIPoint = 0." << endl;
            }

            file << endl;

        }
        else if (gridType == STRUCTGRID)
        {
            int iMin, iMax, jMin, jMax, kMin, kMax, bcType;
            int sLr, sNd;
            StructGrid *grid = StructGridCast(gridArray[iZone]);

            int iDimensions[3];
            iDimensions[0] = grid->GetNI();
            iDimensions[1] = grid->GetNJ();
            iDimensions[2] = grid->GetNK();

            file << "ZoneID : " << iZone << endl;
            file << "iDimensions : Dimensions of different directions in current zone." << endl;
            file << "iDimensions_i" << setw(10) << iDimensions[0] << endl;
            file << "iDimensions_j" << setw(10) << iDimensions[1] << endl;
            file << "iDimensions_k" << setw(10) << iDimensions[2] << endl;

            CompositeStructBC *compositeBCRegion = grid->GetCompositeBCRegion();
            int nBCRegion = compositeBCRegion->GetnBCRegion();
            file << "nBCRegion : " << nBCRegion << endl;

            file << "BCInfor : Boundary information of current zone ." << endl;
            file << "bcNumber  iMin  iMax  jMin  jMax  kMin  kMax  bcType  sLr  sNd  BCName" << endl;
            int nBCTotal = nBCRegion;
            for (int ibcRegion = 0; ibcRegion < nBCRegion; ++ ibcRegion)
            {
                StructBC *bcRegion = compositeBCRegion->GetBCRegion(ibcRegion);
                bcType = bcRegion->GetBCType();
                if (bcType == -1)
                {
                    nBCTotal ++;
                }
            }

            int **BCInfor = NewPointer2<int>(nBCTotal, 9);
            int iBCinfo = 0;

            for (int ibcRegion = 0; ibcRegion < nBCRegion; ++ ibcRegion)
            {
                int index = 0;
                StructBC *bcRegion = compositeBCRegion->GetBCRegion(ibcRegion);
                bcRegion->GetIJKRegion(iMin, iMax, jMin, jMax, kMin, kMax);
                bcType = bcRegion->GetBCType();
                sLr = bcRegion->GetFaceLeftOrRightIndex();
                sNd = bcRegion->GetFaceDirection();
                string bcName = bcRegion->GetBCName();

                BCInfor[iBCinfo][index++] = iMin;
                BCInfor[iBCinfo][index++] = iMax;
                BCInfor[iBCinfo][index++] = jMin;
                BCInfor[iBCinfo][index++] = jMax;
                BCInfor[iBCinfo][index++] = kMin;
                BCInfor[iBCinfo][index++] = kMax;
                BCInfor[iBCinfo][index++] = bcType;
                BCInfor[iBCinfo][index++] = sLr;
                BCInfor[iBCinfo][index++] = sNd;

                file << iBCinfo << setw(8) << iMin << setw(8) << iMax 
                    << setw(8) << jMin << setw(8) << jMax << setw(8) << kMin 
                    << setw(8) << kMax << setw(8) << bcType << setw(8) 
                    << sLr << setw(10) << sNd << setw(20) << bcName <<endl;

                iBCinfo ++;

                if (bcType == -1)
                {
                    index = 0;
                    int *t_st = bcRegion->GetTargetStart();
                    int *t_ed = bcRegion->GetTargetEnd();
                    int neighborZoneIndex = bcRegion->GetTargetRegionBlock();

                    BCInfor[iBCinfo][index++] = t_st[0];
                    BCInfor[iBCinfo][index++] = t_ed[0];
                    BCInfor[iBCinfo][index++] = t_st[1];
                    BCInfor[iBCinfo][index++] = t_ed[1];
                    BCInfor[iBCinfo][index++] = t_st[2];
                    BCInfor[iBCinfo][index++] = t_ed[2];
                    BCInfor[iBCinfo][index++] = neighborZoneIndex;
                    BCInfor[iBCinfo][index++] = 0;
                    BCInfor[iBCinfo][index++] = 0;

                    file << iBCinfo << setw(8) << t_st[0] << setw(8) << t_ed[0] 
                    << setw(8) << t_st[1] << setw(8) << t_ed[1] << setw(8) << t_st[2] 
                    << setw(8) << t_ed[2] << setw(8) << neighborZoneIndex << setw(8) 
                        << 0 << setw(8) << 0 << endl;

                    iBCinfo ++;
                }
            }

            InterfaceInfo *interfaceInfo = grid->GetInterfaceInfo();
            if (interfaceInfo)
            {
                int nIFace = interfaceInfo->GetNIFace();
                file << "nIFace = " << nIFace << endl;

                int *interFace2ZoneID       = interfaceInfo->GetInterFace2ZoneID();
                int *interFace2InterFaceID  = interfaceInfo->GetInterFace2InterFaceID();
                int *interFace2BoundaryFace = interfaceInfo->GetInterFace2BoundaryFace();

                file << "interFace2ZoneID" << endl;
                file << "interFaceNumber  ZoneID" << endl;
                for (int iFace = 0; iFace < nIFace; iFace++)
                {
                    file << iFace << setw(8) << interFace2ZoneID[iFace] << endl;
                }

                file << "interFace2InterFaceID" << endl;
                file << "interFaceNumber  InterFaceID" << endl;
                for (int iFace = 0; iFace < nIFace; iFace++)
                {
                    file << iFace << setw(8) << interFace2InterFaceID[iFace] << endl;
                }

                file << "interFace2BoundaryFace" << endl;
                file << "interFaceNumber  BoundaryFaceID" << endl;
                for (int iFace = 0; iFace < nIFace; iFace++)
                {
                    file << iFace << setw(8) << interFace2BoundaryFace[iFace] << endl;
                }

            }
            else
            {
                file << "There is no interface exist, nIFace = 0." << endl;
            }

            file << endl;

        }
        else
        {
            TK_Exit :: UnexpectedVarValue("gridType", gridType);
        }

    }

    PHSPACE::CloseFile(file);

    delete region;
    delete [] gridArray;
}

void Simulation::RunPostProcessing()
{
    PostProcess *postProcess = new PostProcess(this->region);
    postProcess->Run();

    delete postProcess;

    postProcess = NULL;
}

#ifdef USE_INHOUSE_VERSION
void Simulation::IdentifyDerivative()
{
    Identify *identify = new Identify();
    identify->Run();

    delete identify;

    identify = NULL;
}

void Simulation::RunHelicopterSolver()
{
    MYHelicopter *helicopter = new MYHelicopter(this->region);

    helicopter->ReadGeometry();

    helicopter->Run();

    delete helicopter;

    helicopter = NULL;
}

void Simulation::RunDSMC()
{
    DSMC *dsmc = new DSMC();

    dsmc->Run();

    delete dsmc;

    dsmc = NULL;
}

void Simulation::RunUGKS()
{
    UGKSSimulation *ugks = new UGKSSimulation();

    ugks->Run();

    delete ugks;

    ugks = NULL;
}
#endif

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

    int isOverset = PHSPACE::GlobalDataBase::GetIntParaFromDB("codeOfOversetGrid");
    int isAle = PHSPACE::GlobalDataBase::GetIntParaFromDB("codeOfAleModel");

    if (!isOverset || isAle)
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
            TimeSpan *timeSpan = new TimeSpan();
            WriteLogFile("Run overset configuration ...\n");
            OversetGridConfig();
            timeSpan->ShowTimeSpanToLogfile("Run overset configuration");
            delete timeSpan;
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

void Simulation::RefreshSolvers()
{
    FreeAleManager();

    int nSolver = Zone::GetNumberOfSolvers();
    for (int iSolver = 0; iSolver < nSolver; ++ iSolver)
    {
        Controller *controller = new Controller(iSolver, this->region);

        controller->CleanUp(iSolver);

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
                    int isFVMOrFDM = 0;
                    zone->GetData("isFVMOrFDM", &isFVMOrFDM, PHINT, 1);

                    int iCodeBranch = 0;
                    zone->GetData("iCodeBranch", &iCodeBranch, PHINT, 1);

                    if (isFVMOrFDM == FV_METHOD) //! Add structured FV NS solver.
                    {
#ifdef USE_OLD_CODE
                        if (iCodeBranch == 0)
                        {
                            PHSolver *solver = new NSSolverStructFV();
                            solverList.push_back(solver);
                        }
                        else //! A new code branch is active for Navier-Stokes solver.
                        {
                            PHSolver *solver = new HyFlowSolverStruct();
                            solverList.push_back(solver);
                        }
#endif
                        PHSolver *solver = new NSSolverStructFV();
                        solverList.push_back(solver);
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
                #ifdef USE_IncompressibleSolver
                //! Add incompressible solvers.
                if (gridType == UNSTRUCTGRID)
                {
                    //! Add unstructured incompressible solver.
                    PHSolver *solver = new INCompSolverUnstruct();
                    solverList.push_back(solver);
                }
                else if (gridType == STRUCTGRID)
                {
                    //! Add structured incompressible solver.
                    PHSolver *solver = new INCompSolverStruct();
                    solverList.push_back(solver);
                }
                #endif
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

                    int isFVMOrFDM = 0;
                    zone->GetData("isFVMOrFDM", &isFVMOrFDM, PHINT, 1);

                    if (isFVMOrFDM == FV_METHOD || isFVMOrFDM == NS_METHOD)
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

    if (!isUnsteady)
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
    #ifdef USE_IncompressibleSolver
    //! If incompressible, use incompressible solver.
    int compressible = GlobalDataBase::GetIntParaFromDB("compressible");
    if (compressible == 0)
    {
        using namespace PHMPI;
        int nLocalZones = GetNumberofLocalZones();

        for (int iZone = 0; iZone < nLocalZones; ++ iZone)
        {
            int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);
            Zone *zone = region->GetZone(iZone);
            PHGeometry *geometry = zone->GetGeometry();

            int gridType = geometry->GetGrid()->Type();
            //! Add compressible solvers.
            if (gridType == UNSTRUCTGRID)
            {
                INCompSolverUnstruct *solver = static_cast<INCompSolverUnstruct *>(GlobalSolvers::GetSolver(iZone, 0));
                solver->IncompSloverUnstructMain();
            }
        }

        INCompSolver *incompressibleSolver = new INCompSolver(this->region);
        incompressibleSolver->IncompSloverMain();
        delete incompressibleSolver;
        return;
    }
    #endif

    WriteLogFile("Start iterating ...");
    bool continueSolve = true;

    int maxSimuStep = PHSPACE::GlobalDataBase::GetIntParaFromDB("maxSimuStep");
    int isOverset = PHSPACE::GlobalDataBase::GetIntParaFromDB("codeOfOversetGrid");
    int sysGridType = PHSPACE::GlobalDataBase::GetIntParaFromDB("sys_gridtype");

    //! Wait for all processors here!
    if (!isOverset || sysGridType == STRUCTGRID)
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

void Simulation::SelfAdaptionSolve()
{
#ifdef USE_IncompressibleSolver
    //! If incompressible, use incompressible solver.
    int compressible = GlobalDataBase::GetIntParaFromDB("compressible");
    if (compressible == 0)
    {
        using namespace PHMPI;
        int nLocalZones = GetNumberofLocalZones();

        for (int iZone = 0; iZone < nLocalZones; ++ iZone)
        {
            int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);
            Zone *zone = region->GetZone(iZone);
            PHGeometry *geometry = zone->GetGeometry();

            int gridType = geometry->GetGrid()->Type();
            //! Add compressible solvers.
            if (gridType == UNSTRUCTGRID)
            {
                INCompSolverUnstruct *solver = static_cast<INCompSolverUnstruct *>(GlobalSolvers::GetSolver(iZone, 0));
                solver->IncompSloverUnstructMain();
            }
        }

        INCompSolver *incompressibleSolver = new INCompSolver(this->region);
        incompressibleSolver->IncompSloverMain();
        delete incompressibleSolver;
        return;
    }
#endif

    WriteLogFile("Start iterating ...");
    bool continueSolve = true;

    int maxSimuStep = PHSPACE::GlobalDataBase::GetIntParaFromDB("maxSimuStep");
    int isOverset = PHSPACE::GlobalDataBase::GetIntParaFromDB("codeOfOversetGrid");
    int sysGridType = PHSPACE::GlobalDataBase::GetIntParaFromDB("sys_gridtype");

    //! Wait for all processors here!
    if (!isOverset || sysGridType == STRUCTGRID)
    {
        PrintToWindow("Wait for all processors before iteration ... \n");
        PH_Barrier();
        PrintToWindow("Wait is over,start iteration ... \n");
    }

    //Set the flag to pause the numerical iteration.
    RDouble monitorResidual = 0.0;
    PHSPACE::GlobalDataBase::UpdateData("monitorResidual", &monitorResidual, PHDOUBLE, 1);
    int nPauseKey = 0;
    PHSPACE::GlobalDataBase::UpdateData("nPauseKey", &nPauseKey, PHINT, 1);

    int outIterStep = PHSPACE::GlobalDataBase::GetIntParaFromDB("outnstep");
    if (maxSimuStep == 0)
    {
        continueSolve = false;
    }

    while (continueSolve)
    {
        UpdateOuterIterationStep(outIterStep);

        SolveOneOuterStep();

        int nPauseKey = PHSPACE::GlobalDataBase::GetIntParaFromDB("nPauseKey");

        if ((outIterStep >= maxSimuStep && maxSimuStep > 0) || nPauseKey == 1)
        {
            continueSolve = false;
        }
    }
}

void Simulation::PostSimulation()
{
    this->region->PostSimulation();

    int IsUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
    int isAle = GlobalDataBase::GetIntParaFromDB("codeOfAleModel");
    if (IsUnsteady  && isAle)
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
    if (!isUnsteady)
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
