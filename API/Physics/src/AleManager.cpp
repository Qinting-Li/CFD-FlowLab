#include <iostream>
#include "PHMpi.h"
#include "Precision.h"
#include "Constants.h"
#include "AleManager.h"
#include "Force.h"
#include "Solver.h"
#include "Geometry.h"
#include "SixDofManager.h"
#include "DeformingSolverManager.h"
#include "MovingMesh.h"
#include "Zone.h"
#include "Region.h"
#include "OversetInformation.h"
#include "GlobalDataBase.h"
#include "PHIO.h"
#include "IO_FileName.h"
#include "IO_FileReader.h"
#include "TK_Exit.h"
#include "Task_ServerUpdateInterface.h"
#include "Geo_UnstructGrid.h"
#include "Glb_Dimension.h"

namespace PHSPACE
{

AleManager * aleManager = NULL;

void CreateAleManager()
{
    int IsUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
    if (IsUnsteady == 0)
    {
        return;
    }
    int IsAleSimulation = GlobalDataBase::GetIntParaFromDB("codeOfAleModel");
    if( IsAleSimulation == 0 )
    {
        return;
    }

    aleManager = new AleManager;

    aleManager->Initialize();
}

AleManager * GetAleManager()
{
    return aleManager;
}

void FreeAleManager()
{
    if( aleManager )
    {
        delete aleManager;
        aleManager = NULL;
    }
}

void SolveSingleInnerIterationStepForAleSolvers()
{
    int IsAleSimulation = GlobalDataBase::GetIntParaFromDB("codeOfAleModel");
    if (IsAleSimulation == 0)
    {
        return;
    }

    //准定常问题不需要求解Ale
    int aleStartStrategy = GlobalDataBase::GetIntParaFromDB("aleStartStrategy");
    if (aleStartStrategy == -1)
    {
        return;
    }

    //内迭代第一步备份旧网格（n时刻）信息
    int innStep = GlobalDataBase::GetIntParaFromDB("innstep");
    if (innStep == 1)
    {
        BackupAllZonesOldGrid();
    }

    //如果是内迭代第一步，或者是隐式（紧耦合）方法，需要更新网格几何信息
    if (innStep == 1 || IsImplicitAleMethod() )
    {
        aleManager->Run();

        //SetCurrentGridLevelIndex( 0 );
        //更新各个zone的网格几何信息
        //DoMysteriousGridTask( "COMPUTE_METRICS" );
        aleManager->ComputeMetrics();
        //ComputeMetrics();
        //通讯交界面左右信息。切记！
        aleManager->CommCellCenterData();
        //CommunicateCellCenterData();

        //计算faceNormalVelocityContainer等。对于准定常问题（本版本程序，非定常问题一定），不计算运动速度。
        int methodOfDualTime = GlobalDataBase::GetIntParaFromDB("methodOfDualTime");
        if(methodOfDualTime != 0)
        {
            aleManager->ComputeGridFaceVelocity();
        }

        //如果是重叠网格，则执行挖洞操作
        int codeOfOversetGrid = PHSPACE::GlobalDataBase::GetIntParaFromDB("codeOfOversetGrid");
        int outIterStep = PHSPACE::GlobalDataBase::GetIntParaFromDB("outnstep");
        if (codeOfOversetGrid == 1 && outIterStep != 1)
        {
            //TODO: overset config
            //需要在装配之前将重叠区域的变量存储空间删除
            DeAllocateOversetStorage();

            ResetAllZonesOversetGrid();

            InitializeOversetInterfaceTopology();

            //由于重叠网格信息改变，重叠区域的变量存储空间需要重新申请
            ReAllocateOversetStorage();
        }
        
    }
}

void PostprocessAfterInnerIterationStepForAleSolvers()
{
    int IsAleSimulation = GlobalDataBase::GetIntParaFromDB("codeOfAleModel");
    if (IsAleSimulation == 0)
    {
        return;
    }

    aleManager->UpdateAleForce();

    aleManager->Post();

    aleManager->DumpRestart();
}

AleManager::AleManager()
{
    aleForceManager      = NULL;
    sixDofSolverManager  = NULL;
    deformingManager     = NULL;
    movingMeshManager    = NULL;
    currentAleFile       = NULL;
}

AleManager::~AleManager()
{
    if (aleForceManager     ) delete aleForceManager;
    if (sixDofSolverManager ) delete sixDofSolverManager;
    if (deformingManager    ) delete deformingManager;
    if (movingMeshManager   ) delete movingMeshManager;
}

void AleManager::Initialize()
{
    Allocate();

    int aleStartStrategy = GlobalDataBase::GetIntParaFromDB("aleStartStrategy");
    if(aleStartStrategy ==  1 )
    {
        Restart();
    }

    if(aleStartStrategy == 0 )
    {
        aleManager->Run();

        //SetCurrentGridLevelIndex( 0 );
        //DoMysteriousGridTask( "COMPUTE_METRICS" );
        aleManager->ComputeMetrics();

        aleManager->CommCellCenterData();

        //ExecuteMessageTask( "DEALLOCATEOVERSETSTORAGE" );
        int codeOfOversetGrid = PHSPACE::GlobalDataBase::GetIntParaFromDB("codeOfOversetGrid");
        if (codeOfOversetGrid == 1)
        {
            ResetAllZonesOversetGrid();
            InitializeOversetInterfaceTopology();
        }
    }
}

void AleManager::ComputeMetrics()
{
    using namespace PHMPI;

    int nZones = GetNumberofGlobalZones();
    ActionKey* actkey = new ActionKey();
    actkey->action = COMPUTEMETRICS;

    for (int izone = 0; izone < nZones; ++izone)
    {
        int send_proc = GetZoneProcessorIDSepMode(izone);
        int recv_proc = GetServerSepMode();

        int tag = GetSendRecvTag(GRID_BASED, 0, izone);
        int proc = GetZoneProcessorID(izone);
        int myid = PHMPI::GetCurrentProcessorID();
        if (myid == send_proc)
        {
            //! If the ProcessorID of i-th zone equal to the running ProcessorID, compute grid metrics.
            Grid* grid = GetGrid(izone, 0);
            grid->ComputeMetrics(actkey);
        }

        // collect all info out by ComputeMetrics function to cdata.
        PH_Trade(actkey, send_proc, recv_proc, tag);

        if (myid == recv_proc)
        {
            //! If ProcessorID ProcessorID is serve ProcessorID, output the info in cdata.
            WriteScreen(actkey);
        }
    }
    delete actkey;
}

void AleManager::CommCellCenterData()
{
    using namespace PHMPI;
    int nMGLevel = GlobalDataBase::GetIntParaFromDB("nMGLevel");
    for (int level = 0; level < nMGLevel; ++level)
    {
        ActionKey* actkey = new ActionKey();
        actkey->action = COMMCELLCENTERDATA;
        actkey->kind = GRID_BASED;
        actkey->level = level;

        //Task_ServerUpdateInterface *task = new Task_ServerUpdateInterface();
        if (IsNeedNonBlockingCommunication(actkey))
        {
            //task->MainTaskNonBlocking(actkey);
            MainTaskNonBlocking(actkey);
        }
        else
        {
            //task->MainTaskBlocking(actkey);
            MainTaskBlocking(actkey);
        }

        //delete task;
        delete actkey;
    }
}

void AleManager::ComputeGridFaceVelocity()
{
    int numberOfZones = PHMPI::GetNumberofGlobalZones();

    int gridLevelIndex = 0;
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        Grid* grid = PHSPACE::GetGrid(iZone, gridLevelIndex);
        if (!grid) continue;

        grid->ComputeGridFaceVelocity();
    }
}

void AleManager::Allocate()
{
    //需要首先求一下气动力，用于sixdofSolver的初始化
    //PHSPACE::ExecuteMessageTask( "INITIALIZE_ALLZONE_AERODYNAMICFORCES" );
    InitializeAerodynamicForce();

    aleForceManager         = new AleForceManager       ;
    deformingManager        = new DeformingSolverManager;
    sixDofSolverManager     = new SixDofSolverManager   ;
    movingMeshManager       = new MovingMeshManager     ;

    //初始化参考时间
    InitializeDimensionalReferenceTime();

    //初始化有量纲时间
    RDouble dimensionalAleSimulationTime = 0.0;
    GlobalDataBase::UpdateData("dimensionalAleSimulationTime", &dimensionalAleSimulationTime, PHDOUBLE, 1);

    //初始化无量纲时间 & 初始化内/外迭代步数
    int aleStartStrategy = GlobalDataBase::GetIntParaFromDB("aleStartStrategy");
    if(aleStartStrategy == 0 )
    {
        RDouble physicalTime = 0.0;
        GlobalDataBase::UpdateData("physicalTime", &physicalTime, PHDOUBLE, 1);

        int outnstep = 0;
        GlobalDataBase::UpdateData("outnstep", &outnstep, PHINT, 1);

        int innstep = 0;
        GlobalDataBase::UpdateData("innstep", &innstep, PHINT, 1);
    }

    //一般初始化
    aleForceManager     ->Initialize();
    deformingManager    ->Initialize();
    sixDofSolverManager ->Initialize();
    movingMeshManager   ->Initialize();
}

void AleManager::Restart()
{
    fstream aleFile;
    string aleFileName = "./results/aleRestartFile.dat";

    int aleStartStrategy = GlobalDataBase::GetIntParaFromDB("aleStartStrategy");
    ifstream infile(aleFileName.c_str(), ios::in);
    if (!infile)
    {
        int aleStartStrategy = 0;
        GlobalDataBase::UpdateData("aleStartStrategy", &aleStartStrategy, PHINT, 1);
        return;
    }

    ios_base::openmode openMode = ios_base::in | ios_base::binary;
    PHSPACE::OpenFile( aleFile, aleFileName, openMode );
    currentAleFile = & aleFile;
    SetIsReadAleFile( true );

    RDouble dimensionalAleSimulationTime = 0.0;
    PHRead( aleFile, dimensionalAleSimulationTime);
    GlobalDataBase::UpdateData("dimensionalAleSimulationTime", &dimensionalAleSimulationTime, PHDOUBLE, 1);
    //SetDimensionalAleSimulationTime( aleTime );

    //如果是重启，动网格需要从0时刻转换到n时刻
    aleForceManager       ->Restart();
    deformingManager      ->Restart();
    sixDofSolverManager   ->Restart();
    movingMeshManager     ->Run();

    //网格运动了，需要重新求几何信息
    aleManager->ComputeMetrics();

    aleManager->CommCellCenterData();

    PHSPACE::CloseFile( aleFile );
}

void AleManager::InitializeDimensionalReferenceTime()
{
    RDouble referenceLength    = GlobalDataBase::GetDoubleParaFromDB("referenceLength");
    RDouble referenceVelocity  = GlobalDataBase::GetDoubleParaFromDB("referenceVelocity");

    RDouble dimensionalReferenceTime = referenceLength / referenceVelocity;

    GlobalDataBase::UpdateData("dimensionalReferenceTime", &dimensionalReferenceTime, PHDOUBLE, 1);
}

void AleManager::UpdateAleForce()
{
    //气动力还是放在迭代后进行更新，以防止重启的时候得到错误的能耗
    aleForceManager->Run();
}

void AleManager::Run()
{
    RDouble dimensionalAleSimulationTime = GlobalDataBase::GetDoubleParaFromDB("dimensionalAleSimulationTime");
    RDouble dimensionalAleTimeStep = GetDimensionalAleTimeStep();

    dimensionalAleSimulationTime += dimensionalAleTimeStep;
    GlobalDataBase::UpdateData("dimensionalAleSimulationTime", &dimensionalAleSimulationTime, PHDOUBLE, 1);

    deformingManager    ->Run();
    sixDofSolverManager ->Run();
    movingMeshManager   ->Run();
}

void AleManager::Post()
{
    aleForceManager     ->Post();
    deformingManager    ->Post();
    sixDofSolverManager ->Post();
    movingMeshManager   ->Post();
}

void AleManager::Dump()
{
    deformingManager      ->Dump();
    sixDofSolverManager   ->Dump();
}

void AleManager::DumpRestart()
{
    using namespace PHMPI;
    int currentProcessor = GetCurrentProcessorID();
    int serverProcessor = GetServerProcessorID();
    if (currentProcessor != serverProcessor)
    {
        return;
    }

    fstream aleFile;

    int outerStep = GlobalDataBase::GetIntParaFromDB("outnstep");
    int numberOfFieldSaveSteps = GlobalDataBase::GetIntParaFromDB("intervalStepFlow");

    if( outerStep % numberOfFieldSaveSteps == 0 )
    {
        string aleFileName = "./results/aleRestartFile.dat";
        ios_base::openmode openMode = ios_base::out | ios_base::binary;
        PHSPACE::OpenFile( aleFile, aleFileName, openMode );

        RDouble aleTime = GlobalDataBase::GetDoubleParaFromDB("dimensionalAleSimulationTime");

        PHWrite( aleFile, aleTime );

        aleForceManager->DumpRestartAleForceOfBodyContainer( aleFile );
        deformingManager->DumpRestart(aleFile);
        sixDofSolverManager->DumpRestart( aleFile );

        PHSPACE::CloseFile(aleFile);
    }
}


void BackupAllZonesOldGrid()
{
    int numberOfZones = PHMPI::GetNumberofGlobalZones();

    int gridLevelIndex = 0;
    for ( int iZone = 0; iZone < numberOfZones; ++ iZone )
    {
        Grid* grid = GetGrid( iZone, gridLevelIndex );
        if ( ! grid ) continue;
        grid->BackUpOldGrid();
    }
}

void ComputeMetrics()
{
    int numberOfZones = PHMPI::GetNumberofGlobalZones();
    ActionKey* actkey = new ActionKey();
    actkey->action = COMPUTEMETRICS;

    int gridLevelIndex = 0;
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        Grid* grid = GetGrid(iZone, gridLevelIndex);
        if (!grid) continue;
        grid->ComputeMetrics(actkey);
    }
}


void DeAllocateOversetStorage()
{
    int nSolver = Zone::GetNumberOfSolvers();
    for (int iSolver = 0; iSolver < nSolver; ++iSolver)
    {
        int numberOfZones = PHMPI::GetNumberofGlobalZones();
        int gridLevelIndex = 0;
        for (int iZone = 0; iZone < numberOfZones; ++iZone)
        {
            Grid* grid = GetGrid(iZone, gridLevelIndex);
            //CFDSolver* solver = GetCFDSolver(iZone, iSolver);
            if (!grid) continue;
            if (iSolver == NS_SOLVER)
            {
                int nNSEquation = GlobalDataBase::GetIntParaFromDB("nm");
                int nTemperatureModel = GlobalDataBase::GetIntParaFromDB("ntmodel");
                InterfaceFieldNameProperty< PHNsTag >::AddFieldNameProperty("q", nNSEquation, TWO_D);
                InterfaceFieldNameProperty< PHNsTag >::AddFieldNameProperty("t", nTemperatureModel, TWO_D);
                GreatDeAllocateOversetInterfaceVariables< PHNsTag >(grid);
            }
            else if(iSolver == TURB_SOLVER)
            {
                int nEquation = GlobalDataBase::GetIntParaFromDB("n_turb");
                InterfaceFieldNameProperty< PHNsTag >::AddFieldNameProperty("turb::q", nEquation, TWO_D);
                GreatDeAllocateOversetInterfaceVariables< PHTurbulentTag >(grid);
            }
        }

    }
}

void ReAllocateOversetStorage()
{
    int nSolver = Zone::GetNumberOfSolvers();
    for (int iSolver = 0; iSolver < nSolver; ++iSolver)
    {
        int numberOfZones = PHMPI::GetNumberofGlobalZones();
        int gridLevelIndex = 0;
        for (int iZone = 0; iZone < numberOfZones; ++iZone)
        {
            Grid* grid = GetGrid(iZone, gridLevelIndex);
            //CFDSolver* solver = GetCFDSolver(iZone, iSolver);
            if (!grid) continue;
            if (iSolver == NS_SOLVER)
            {
                int nNSEquation = GlobalDataBase::GetIntParaFromDB("nm");
                int nTemperatureModel = GlobalDataBase::GetIntParaFromDB("ntmodel");
                InterfaceFieldNameProperty< PHNsTag >::AddFieldNameProperty("q", nNSEquation, TWO_D);
                InterfaceFieldNameProperty< PHNsTag >::AddFieldNameProperty("t", nTemperatureModel, TWO_D);
                GreatAllocateOversetInterfaceVariables< PHNsTag >(grid);
            }
            else if (iSolver == TURB_SOLVER)
            {
                int nEquation = GlobalDataBase::GetIntParaFromDB("n_turb");
                InterfaceFieldNameProperty< PHNsTag >::AddFieldNameProperty("turb::q", nEquation, TWO_D);
                GreatAllocateOversetInterfaceVariables< PHTurbulentTag >(grid);
            }
        }

    }
}

template< typename T >
void InterfaceFieldNameProperty< T >::AllocateInterfaceVariables(int numberOfInterfaces, Data_ParamFieldSuite* dataStorage)
{
    if (numberOfInterfaces <= 0) return;

    for (std::map< string, AllocateParameter >::iterator iter = InterfaceFieldNameProperty< T >::fieldNameParameterData.begin(); iter != InterfaceFieldNameProperty< T >::fieldNameParameterData.end(); ++iter)
    {
        int numberOfEquations = iter->second.GetNumberOfEquations();

        RDouble** field = NewPointer2< RDouble >(numberOfEquations, numberOfInterfaces);
        ZeroField(field, numberOfEquations, numberOfInterfaces);
        dataStorage->UpdateDataPtr(iter->first, field);
    }
}

template< typename T >
void InterfaceFieldNameProperty< T >::DeAllocateInterfaceVariables(Data_ParamFieldSuite* dataStorage)
{
    for (std::map< string, AllocateParameter >::iterator iter = InterfaceFieldNameProperty< T >::fieldNameParameterData.begin(); iter != InterfaceFieldNameProperty< T >::fieldNameParameterData.end(); ++iter)
    {
        RDouble** q = reinterpret_cast <RDouble**> (dataStorage->GetDataPtr(iter->first));
        DelPointer2(q);
    }
}

template< typename T >
std::map< string, AllocateParameter > InterfaceFieldNameProperty< T >::fieldNameParameterData;

template< typename T >
InterfaceFieldNameProperty< T >::InterfaceFieldNameProperty()
{
}

template< typename T >
InterfaceFieldNameProperty< T >::~InterfaceFieldNameProperty()
{
}

template< typename T >
void InterfaceFieldNameProperty< T >::AddFieldNameProperty(const string& fieldName, int numberOfEquations)
{
    AllocateParameter allocateParameter;
    allocateParameter.SetNumberOfEquations(numberOfEquations);
    allocateParameter.SetDimensions(TWO_D);

    InterfaceFieldNameProperty< T >::fieldNameParameterData[fieldName] = allocateParameter;
}

template< typename T >
void InterfaceFieldNameProperty< T >::AddFieldNameProperty(const string& fieldName, int numberOfEquations, int dimensions)
{
    AllocateParameter allocateParameter;
    allocateParameter.SetNumberOfEquations(numberOfEquations);
    allocateParameter.SetDimensions(dimensions);

    InterfaceFieldNameProperty< T >::fieldNameParameterData[fieldName] = allocateParameter;
}

template< typename T >
int InterfaceFieldNameProperty< T >::GetNumberOfEquations(const string& fieldName)
{
    std::map< string, AllocateParameter >::iterator iter;
    iter = InterfaceFieldNameProperty< T >::fieldNameParameterData.find(fieldName);
    if (iter != InterfaceFieldNameProperty< T >::fieldNameParameterData.end())
    {
        return iter->second.GetNumberOfEquations();
    }
    return -1;
}

template< typename T >
int InterfaceFieldNameProperty< T >::GetNumberOfDimensions(const string& fieldName)
{
    std::map< string, AllocateParameter >::iterator iter;
    iter = InterfaceFieldNameProperty< T >::fieldNameParameterData.find(fieldName);
    if (iter != SuperFieldNameProperty< T >::fieldNameParameterData.end())
    {
        return iter->second.GetDimensions();
    }
    return -1;
}

template< typename T >
void GreatAllocateOversetInterfaceVariables(Grid* gridIn)
{
    OversetInformationProxy* oversetInformationProxy = gridIn->GetOversetInformationProxy();
    if (!oversetInformationProxy) return;

    OversetDonorCellManager* oversetDonorCellManager = oversetInformationProxy->GetOversetDonorCellManager();
    OversetInterplateCellManager* oversetInterplateCellManager = oversetInformationProxy->GetOversetInterplateCellManager();

    InterfaceFieldNameProperty< T >::AllocateInterfaceVariables(oversetDonorCellManager->GetTotalNumberOfOversetDonorCells(), oversetDonorCellManager->GetDataStorage());
    InterfaceFieldNameProperty< T >::AllocateInterfaceVariables(oversetInterplateCellManager->GetTotalNumberOfOversetInterplateCells(), oversetInterplateCellManager->GetDataStorage());
}

template< typename T >
void GreatDeAllocateOversetInterfaceVariables(Grid* gridIn)
{
    OversetInformationProxy* oversetInformationProxy = gridIn->GetOversetInformationProxy();
    if (!oversetInformationProxy) return;

    OversetDonorCellManager* oversetDonorCellManager = oversetInformationProxy->GetOversetDonorCellManager();
    OversetInterplateCellManager* oversetInterplateCellManager = oversetInformationProxy->GetOversetInterplateCellManager();

    InterfaceFieldNameProperty< T >::DeAllocateInterfaceVariables(oversetDonorCellManager->GetDataStorage());
    InterfaceFieldNameProperty< T >::DeAllocateInterfaceVariables(oversetInterplateCellManager->GetDataStorage());
}

bool isReadAleFile = false;

bool IsReadAleFile()
{
    return isReadAleFile;
}

void SetIsReadAleFile( bool isReadAleFileIn )
{
    PHSPACE::isReadAleFile = isReadAleFileIn;
}

fstream & GetCurrentAleFile()
{
    return aleManager->GetCurrentAleFile();
}

bool IsImplicitAleMethod()
{
    int keyMethod = PHSPACE::GlobalDataBase::GetIntParaFromDB( "methodForKineticEquation" );

    if ( keyMethod == KINETIC_SCHEME::ADMAS_BASHFORTH_FIRST_ORDER )   return false;
    if ( keyMethod == KINETIC_SCHEME::ADMAS_BASHFORTH_SECOND_ORDER )  return false;

    return true;
}

bool IsExplicitAleMethod()
{
    return ! IsImplicitAleMethod();
}

int GetIntegerParameterFromDataBase(int bodyIndex, const string& parameterName)
{
    ostringstream oss;
    oss << parameterName << "_" << bodyIndex;
    return PHSPACE::GlobalDataBase::GetIntParaFromDB(oss.str());
}

RDouble GetRDoubleParameterFromDataBase(int bodyIndex, const string& parameterName)
{
    ostringstream oss;
    oss << parameterName << "_" << bodyIndex;
    return PHSPACE::GlobalDataBase::GetDoubleParaFromDB(oss.str());
}

string GetStringParameterFromDataBase(int bodyIndex, const string& parameterName)
{
    ostringstream oss;
    oss << parameterName << "_" << bodyIndex;
    return PHSPACE::GlobalDataBase::GetStrParaFromDB(oss.str());
}

int* GetIntegerArrayFromDataBase(const string& nameOfIntegerArray, int numberOfElements)
{
    int* integerArray = new int[numberOfElements];

    GlobalDataBase::GetData(nameOfIntegerArray, integerArray, PHINT, numberOfElements);

    return integerArray;

}

void GetRDoubleVectorFromDataBase(vector< RDouble >& parameter, int bodyIndex, const string& parameterName, int numberOfElements)
{
    ostringstream oss;
    oss << parameterName << "_" << bodyIndex;
    PHSPACE::GlobalDataBase::GetRDoubleVectorFromDB(parameter, oss.str(), numberOfElements);
}

RDouble GetDimensionalAleTimeStep()
{
    //RDouble reynoldsReferenceLengthDimensional = GlobalDataBase::GetDoubleParaFromDB("reynoldsReferenceLengthDimensional");

    //RDouble refDimensionalSonicSpeed = GlobalDataBase::GetDoubleParaFromDB("refDimensionalSonicSpeed");
    //RDouble refMachNumber = GlobalDataBase::GetDoubleParaFromDB("refMachNumber");
    //RDouble refVelocityDimensional = refMachNumber * refDimensionalSonicSpeed;
    //RDouble referenceVelocity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalVelocity");
    RDouble referenceLength = GlobalDataBase::GetDoubleParaFromDB("referenceLength");

    RDouble refVelocityDimensional = GlobalDataBase::GetDoubleParaFromDB("refDimensionalVelocity");
    RDouble unit = referenceLength / refVelocityDimensional;
    RDouble physicalTimeStep = PHSPACE::GlobalDataBase::GetDoubleParaFromDB("physicalTimeStep");
    RDouble dimensionalAleTimeStep = physicalTimeStep * unit;
    //RDouble dimensionalAleTimeStep = physicalTimeStep;
    return dimensionalAleTimeStep;
}

}