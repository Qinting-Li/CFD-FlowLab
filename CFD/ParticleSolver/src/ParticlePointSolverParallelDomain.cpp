#include "ParticlePointSolverParallelDomain.h"
#include "Param_ParticleSolver.h"
#include "Param_ParticlePointSolverParallelDomain.h"
#include "LIB_Macro.h"
#include "ParticlePointGroup.h"
#include "SimplePointer.h"
#include "IO_FileName.h"
#include "TK_Log.h"
#include "Glb_Dimension.h"
#include "Parallel.h"
#include "TK_Exit.h"
#include "Geo_Grid.h"
#include "Interpolation.h"
#include "LagrangianVariable.h"
#include "Geo_StructBC.h"
#include "ParticleForceModel.h"
#include "Math_BasisFunction.h"
#include <algorithm>
#include <random>
#include <fstream>
using namespace std;

namespace PHSPACE
{

ParticlePointSolverParallelDomain::ParticlePointSolverParallelDomain()
{

}

ParticlePointSolverParallelDomain::~ParticlePointSolverParallelDomain()
{
    //! Destructors can cause many problems of deletion by mistake.
    //! We must be careful. You can comment out this sentence 
    //! after the program debugging is perfect.
    DeAllocateGlobalVariables();
    FreeControlParameters();
}

LIB_EXPORT Param_ParticlePointSolverParallelDomain* ParticlePointSolverParallelDomain::GetControlParameters() const
{
    return static_cast <Param_ParticlePointSolverParallelDomain*> (controlParameters);
}

void ParticlePointSolverParallelDomain::InitMemory()
{
    //! new param on each zone.
    InitControlParameters();

    //! read param file one each zone.
    //! can change to parallel reading in future.
    ReadParameter();
 
    //! particle array is carried out on each zone.
    AllocateGlobalVariables();
}

void ParticlePointSolverParallelDomain::InitControlParameters()
{
    //! Note here, Param_CFDSolver *controlParameters is Class member of CFDSolver.
    FreeControlParameters();
    //! Note that new param_* one each zone.
    controlParameters = new Param_ParticlePointSolverParallelDomain();
    controlParameters->Init();
}

void ParticlePointSolverParallelDomain::ReadParameter()
{
    Param_ParticlePointSolverParallelDomain *param_ParticlePointSolverParallelDomain = GetControlParameters();
    param_ParticlePointSolverParallelDomain->ReadParameter();
}

void ParticlePointSolverParallelDomain::AllocateGlobalVariables()
{
    Grid *grid = GetGrid();
    int gridType = grid->Type();

    CFDSolver::AllocateGlobalVar(grid);
    Param_ParticlePointSolverParallelDomain *parameter = GetControlParameters();

    //! Note that here we only declare the pointer address of the variable, 
    //! but do not initialize the size of the variable.

    ParticlePointGroup *particlePointGroup = new ParticlePointGroup();
    grid->UpdateDataPtr("particlePointGroup", particlePointGroup);

    //! If init particle BC file.
    ParticlePointGroup *particleBCInit = new ParticlePointGroup();
    grid->UpdateDataPtr("particleBCInit", particleBCInit);

    //! Inie the IndexCell.
    particlePointGroup->InitCellIndex(grid);

    //! Variable for flow
    //! Init flow angular. 
    if (gridType == PHSPACE::STRUCTGRID)
    {
        int nLayers = GetNumberOfGhostCellLayers();

        StructGrid *structgrid = StructGridCast(grid);
        int ni = structgrid->GetNI();
        int nj = structgrid->GetNJ();
        int nk = structgrid->GetNK();
        Range I, J, K;
        //! eg,nlayers = 0,start from 1 to ni - 1.
        //! eg,nlayers = 2,start from 1-nlayers to ni - 1 + nlayers.
        GetRange(ni, nj, nk, -nLayers, nLayers - 1, I, J, K);

        if ("Gimenez2019JCP" == parameter->GetParticleInterpolationStruct())
        {
            Range M(0, INDEX_FLOWONPARTICLE::nVarOnFlow - 1);
            RDouble4D *flowOnParticleInEular = new RDouble4D(I, J, K, M, fortranArray);
            grid->UpdateDataPtr("flowOnParticleInEular", flowOnParticleInEular);

            Range MGrad(0, INDEX_FLOWONPARTICLE::nVarOnFlow * 3 - 1);
            RDouble4D *flowOnParticleInEularGrad = new RDouble4D(I, J, K, MGrad, fortranArray);
            grid->UpdateDataPtr("flowOnParticleInEularGrad", flowOnParticleInEularGrad);

            //! Hessian matrix
            Range MGradSecond(0, INDEX_FLOWONPARTICLE::nVarOnFlow * 3 * 3 - 1);
            RDouble4D *flowOnParticleInEularGradSecond = new RDouble4D(I, J, K, MGradSecond, fortranArray);
            grid->UpdateDataPtr("flowOnParticleInEularGradSecond", flowOnParticleInEularGradSecond);

            *flowOnParticleInEular = 0.0;
            *flowOnParticleInEularGrad = 0.0;
            *flowOnParticleInEularGradSecond = 0.0;
        }

        //! Back coupling, where 0,1,2 are force feedback and 3 is temperature feedback.
        Range P(0, 3);
        RDouble4D *sourceParticle2Flow = new RDouble4D(I, J, K, P, fortranArray);
        *sourceParticle2Flow = 0.0;
        grid->UpdateDataPtr("sourceParticle2Flow", sourceParticle2Flow);
    }
    else if (gridType == PHSPACE::UNSTRUCTGRID)
    {
        cout << "Particle Solver Warning: Unstructured grid is not supported yet!" << endl;
    }
    else
    {
        cout << "Particle Solver Warning: Only structured grid is supported yet!" << endl;
    }
}

void ParticlePointSolverParallelDomain::InitFlow()
{
    //! Read particle HDF5 file
    ActionKey *actkeyReadParticleFile = new ActionKey();
    FillActionKey(actkeyReadParticleFile, READ_PARTICLE, 0);

    actkeyReadParticleFile->subact = 0;
    ReadParticleH5(actkeyReadParticleFile);

    actkeyReadParticleFile->subact = 1;
    ReadParticleH5(actkeyReadParticleFile);

    delete actkeyReadParticleFile;

    //! Calculate flow on particle interpolation.
    if (GetControlParameters()->ifInitParticleInFlow())
    {
        ActionKey *actkeyCommFlowOnParticle = new ActionKey();
        FillActionKey(actkeyCommFlowOnParticle, FLOW_ON_PARTICLE_MPI, 0);

        //! Init the flow on flowOnParticleInEular.
        actkeyCommFlowOnParticle->subact = 0;
        CalculateVariableForParticleOnEular(actkeyCommFlowOnParticle);
        CalcFlowOnLocationOfParticle(actkeyCommFlowOnParticle);
        delete actkeyCommFlowOnParticle;
    }
    
    //! Sum all particle number on each group.
    CheckNumOfParticle();

    //! initialize LDouble variables. 
    InitLagrangianVariable();

    //! Check the ghost and corner variable on grid.
    CheckCornerAndGhost();

    //! Calculate CFL number of particle.
    CalcCFLNumber();

    //! Source for restart.
    GetSourceForRestart();
}

bool ParticlePointSolverParallelDomain::JudgeIfRestart()
{
    Param_ParticlePointSolverParallelDomain *param_ParticlePointSolverParallelDomain = GetControlParameters();
    bool readRestartFlag = param_ParticlePointSolverParallelDomain->IfReadRestartParticle();
    if (!readRestartFlag)
    {
        return false;
    }
    else
    {
        string restartParticleFile = GetRestartFileName();
        if (restartParticleFile == "")
        {
            return false;
        }

        ifstream infile(restartParticleFile.c_str(), ios::in);
        if (infile)
        {
            infile.close();
            infile.clear();
            return true;
        }
        else
        {
            return false;
        }
    }
}

bool ParticlePointSolverParallelDomain::JudgeIfWriteRestart()
{
    Param_ParticlePointSolverParallelDomain *param_ParticlePointSolverParallelDomain = GetControlParameters();
    bool writeRestartFlag = param_ParticlePointSolverParallelDomain->IfWirteRestartParticle();
    if (!writeRestartFlag)
    {
        return false;
    }
    else
    {
        string restartParticleFile = GetRestartFileName();
        if (restartParticleFile == "")
        {
            return false;
        }

        return true;
    }
}

bool ParticlePointSolverParallelDomain::JudgeIfInit()
{
    Param_ParticlePointSolverParallelDomain *param_ParticlePointSolverParallelDomain = GetControlParameters();
    bool initFlag = param_ParticlePointSolverParallelDomain->IfInitParticle();
    if (!initFlag)
    {
        return false;
    }
    else
    {
        string startParticleFile = GetInitFileName();
        if (startParticleFile == "")
        {
            return false;
        }
        ifstream infile(startParticleFile.c_str(), ios::in);
        if (infile)
        {
            infile.close();
            infile.clear();
            return true;
        }
        else
        {
            ostringstream oss;
            oss << "Need init particle file :" << startParticleFile << "\n";
            oss << "But there is no file." << endl;
            TK_Exit::ExceptionExit(oss);
            return false;
        }
    }
}

void ParticlePointSolverParallelDomain::FillActionKey(ActionKey *actkey, int action, int level)
{
    actkey->action = action;
    actkey->solver = PARTICLE;
    actkey->solverID = this->GetIndex();
    actkey->kind = SOLVER_BASED;
    actkey->level = level;
}

void ParticlePointSolverParallelDomain::ReadParticleH5(ActionKey *actkey)
{
    int zoneIDofThisSolver = GetZoneLocalID();
    if (zoneIDofThisSolver != 0)
    {
        //! Very Important: only the first zone on a processor need to implement the task.
        //! other wise, the tasks would be implemented several times.
        return;
    }

    Param_ParticlePointSolverParallelDomain *parameter = GetControlParameters();

    if (0 == actkey->subact)
    {
        //! init particle file.
        actkey->filename = GetInitFileName();
    }
    else if (1 == actkey->subact)
    {
        //! init particle bc file.
        if (parameter->ifInitParticleBCInfo())
        {
            actkey->filename = parameter->GetInitParticleBCFile();
        }
        else
        {
            return;
        }
    } 

    if (actkey->filename == "")
    {
        TK_Exit::UnexpectedVarValue("No init particle file :", actkey->filename);
    }

    //! Print messege for read particle.
    ostringstream ossLog; 
    ossLog << " Read particle file : " << actkey->filename ;
    ossLog << " on process :" << PHMPI::GetServerProcessorID() << "\n";
    ossLog << endl;
    bool writeEachProcessor = true;
    
    PrintToWindow(ossLog);

    //! Open HDF5 file.
    //! Note that each process can only be opened once when local zone = 0, 
    //! and then distributed to other zones of the current region.
    hid_t file;
    file = OpenHDF5File(actkey->filename);
    actkey->filepos = file;

    using namespace PHMPI;
    //! The number of zones on all processor global.
    int nZones = GetNumberofGlobalZones();
    //! Get current processor ID.
    int currentProcessorID = GetCurrentProcessorID();
    for (int iZone = 0; iZone < nZones; ++iZone)
    {
        //! Judge the zones on current processor.
        int recvProcess = GetZoneProcessorID(iZone);

        //! The loop for zone on current processor.
        if (currentProcessorID == recvProcess)
        {
            //! For each zone'solver.
            ParticlePointSolverParallelDomain *solver = static_cast <ParticlePointSolverParallelDomain*> (GlobalSolvers::GetSolver(iZone, actkey->solverID));
            solver->ReadParticleHDF5File(actkey);
        }
    }

    //! Close HDF5 file.
    H5Fclose(file);
    actkey->filepos = 0;
}

const string ParticlePointSolverParallelDomain::GetInitFileName()
{
    Param_ParticlePointSolverParallelDomain *param_ParticlePointSolverParallelDomain = GetControlParameters();
    return param_ParticlePointSolverParallelDomain->GetInitParticleFile();
}

void ParticlePointSolverParallelDomain::ReadParticleHDF5File(ActionKey *actkey)
{
    using namespace PHMPI;
    //! Get grid, parameter and particlePointGroup.
    int level = actkey->level;
    Grid *grid = GetGrid(level);

    int nZones = GetNumberofGlobalZones();
    int currentProcessorID = GetCurrentProcessorID();

    //! The global zone id.
    int zoneIDGlobal = grid->GetZoneID();
    //! The local zone id on current processer.
    int zoneIDLocal = grid->GetZoneLocalID();

    ostringstream ossReadParticle;
    ossReadParticle << " --- Read Particle HDF5---" << "\n";
    ossReadParticle << " zone ID Global :" << zoneIDGlobal << "\n";
    ossReadParticle << " zone ID Local :" << zoneIDLocal << endl;
    ossReadParticle << " currentProcessorID : " << currentProcessorID << endl;
    bool writeEachProcessor = true;
    WriteLogFile(ossReadParticle, writeEachProcessor);
    ossReadParticle.clear();
    ossReadParticle.str("");

    Param_ParticlePointSolverParallelDomain *parameter = GetControlParameters();

    ParticlePointGroup *particlePointGroup;
    if (actkey->subact == 0)
    {
        particlePointGroup = GetParticlePointGroup("particlePointGroup");
    }
    else
    {
        particlePointGroup = GetParticlePointGroup("particleBCInit");
    }

    //! ref dimensional length
    RDouble refReynoldsLengthDimensional = parameter->GetRefDimLength();
    //! ref dimensional velocity
    RDouble refVelocityDimensional = parameter->GetRefDimVelocity();
    //! ref dimensional angularVelocity
    RDouble refAngularVelocityDimensional = parameter->GetRefDimAngularVelocity();
    //! ref dimensional acceleration
    RDouble refAccelerationyDimensional = parameter->GetRefDimAcceleration();
    //! ref dimensional temperature
    RDouble refTemperatureDimensional = parameter->GetRefDimTemperature();
    //! ref dimensional density
    RDouble refDensityDimensional = parameter->GetRefDimDensity();
    //! ref dimensional force
    RDouble refForceDimensional = parameter->GetRefDimForce();

    hid_t grploc;
    string grpName;

    //! When use ReadDoubleOneRow, it will new array.

    //! Read the version of particle file.
    grpName = "Version";
    grploc = OpenGroup(actkey->filepos, grpName);
    RDouble version = 0.0;
    ReadData(grploc, &version, "version");
    RDouble versionParam = parameter->GetVersion();
    if (version != versionParam)
    {
        ostringstream ossExit;
        ossExit << "Error on particle file version" << "\n";
        ossExit << "version in hdf5 : " << version << "\n";
        ossExit << "version in param : " << versionParam << endl;
        TK_Exit::ExceptionExit(ossExit);
    }
    H5Gclose(grploc);

    //! Read the number of particles.
    grpName = "nGroup";
    grploc = OpenGroup(actkey->filepos, grpName);
    int nParticleTotal = 0;
    int nGroup = 0;
    ReadData(grploc, &nParticleTotal, "nParticleTotal");
    ReadData(grploc, &nGroup, "nGroup");
    H5Gclose(grploc);

    int fileType = 0;
    int iGroup = 0;

    if (actkey->subact == 0)
    {
        //! fileType : type of init particle.
        //! 0, init read only one group.
        //! 1, init read for each group.
        //! 2, restart file for each group.
        if (JudgeIfRestart())
        {
            //! For restart particle file.
            iGroup = zoneIDGlobal;
            fileType = 2;
        }
        else if (JudgeIfInit())
        {
            if (nGroup > 1)
            {
                //! For parallel init particle on each zone.
                iGroup = zoneIDGlobal;
                fileType = 1;
            }
            else
            {
                //! The init file is by Matlab one Group
                iGroup = 0;
                nGroup = nZones;
                fileType = 0;
            }
        }
        else
        {
            ostringstream ossError;
            ossError << "Error on particle file" << "\n";
            ossError << "JudgeIfRestart : " << JudgeIfRestart() << "\n";
            ossError << "Restart particle file : " << GetRestartFileName() << "\n";
            ossError << "JudgeIfInit : " << JudgeIfInit() << "\n";
            ossError << "Init particle file : " << GetInitFileName() << endl;
            TK_Exit::ExceptionExit(ossError);
        }

        int restartOutnstep = 0;
        if (2 == fileType)
        {
            grpName = "Time";
            grploc = OpenGroup(actkey->filepos, grpName);
            ReadData(grploc, &restartOutnstep, "outnstep");
            H5Gclose(grploc);
        }
        parameter->SetInitOutStep(restartOutnstep);
    }
    else
    {
        //! Note that, the init BC particle file must be initfile.
        //! That mean, the group of BC particle file is only one.
        iGroup = 0;
        nGroup = nZones;
        fileType = 0;
    }

    ostringstream particleGroupPath;
    ostringstream particleVarPath;
    string groupNameString;

    particleGroupPath << "Particle_Group" << iGroup;

    //! Read data on each group.
    particleVarPath << particleGroupPath.str() << "/Info";
    groupNameString = particleVarPath.str();
    grpName = groupNameString;
    grploc = OpenGroup(actkey->filepos, grpName);
    int nDimTem = 1;
    int nParticleLocal = 0;
    ReadData(grploc, &nDimTem, "nDimTem");
    ReadData(grploc, &nParticleLocal, "nParticle");

    int *particleID = ReadIntOneRow(grploc, "ParticleID");
    H5Gclose(grploc);
    particleVarPath.clear();
    particleVarPath.str("");

    if (0 == actkey->subact)
    {
        //! Add param to param_solver.
        if (0 == nParticleLocal)
        {
            if (fileType != 0)
            {
                //! For multi-group read file.
                parameter->SetNumParticleLocal(nParticleLocal);
                return;
            }
            else
            {
                //! For init only one group.
                nParticleLocal = nParticleTotal;
            }
        }
        else
        {
            if (fileType != 0)
            {
                parameter->SetNumParticleLocal(nParticleLocal);
            }
        }
        parameter->SetNumParticleTotal(nParticleTotal);
        parameter->SetNumParticleGroup(nGroup);
    }
    else
    {
        //! For init BC file .
        if (0 == nParticleLocal)
        {
            nParticleLocal = nParticleTotal;
        }
    }

    //! Basic param. 3unit
    particleVarPath << particleGroupPath.str() << "/BasicParam";
    groupNameString = particleVarPath.str();
    grpName = groupNameString;
    grploc = OpenGroup(actkey->filepos, grpName);
    RDouble *diameter = ReadDoubleOneRow(grploc, "Diameter");
    RDouble *density = ReadDoubleOneRow(grploc, "Density");
    RDouble *specificHeatCapacity = ReadDoubleOneRow(grploc, "SpecificHeatCapacity");
    RDouble *particleEmissivity = ReadDoubleOneRow(grploc, "ParticleEmissivity");
    H5Gclose(grploc);
    particleVarPath.clear();
    particleVarPath.str("");

    //! Read particle coodinate.
    particleVarPath << particleGroupPath.str() << "/Coordinate";
    groupNameString = particleVarPath.str();
    grpName = groupNameString;
    grploc = OpenGroup(actkey->filepos, grpName);
    RDouble *particleCoordinateX = ReadDoubleOneRow(grploc, "CoordinateX");
    RDouble *particleCoordinateY = ReadDoubleOneRow(grploc, "CoordinateY");
    RDouble *particleCoordinateZ = ReadDoubleOneRow(grploc, "CoordinateZ");
    H5Gclose(grploc);
    particleVarPath.clear();
    particleVarPath.str("");

    //! Read particle velocity.
    particleVarPath << particleGroupPath.str() << "/Velocity";
    groupNameString = particleVarPath.str();
    grpName = groupNameString;
    grploc = OpenGroup(actkey->filepos, grpName);
    RDouble *particleVelocityX = ReadDoubleOneRow(grploc, "VelocityX");
    RDouble *particleVelocityY = ReadDoubleOneRow(grploc, "VelocityY");
    RDouble *particleVelocityZ = ReadDoubleOneRow(grploc, "VelocityZ");
    H5Gclose(grploc);
    particleVarPath.clear();
    particleVarPath.str("");

    RDouble *particleAngularVelocityX = 0,  *particleAngularVelocityY = 0,  *particleAngularVelocityZ = 0;
    if (parameter->ifUseParticleRotate())
    {
        //! Read particle Angular velocity.
        particleVarPath << particleGroupPath.str() << "/AngularVelocity";
        groupNameString = particleVarPath.str();
        grpName = groupNameString;
        grploc = OpenGroup(actkey->filepos, grpName);
        RDouble *particleAngularVelocityX = ReadDoubleOneRow(grploc, "AngularVelocityX");
        RDouble *particleAngularVelocityY = ReadDoubleOneRow(grploc, "AngularVelocityY");
        RDouble *particleAngularVelocityZ = ReadDoubleOneRow(grploc, "AngularVelocityZ");
        H5Gclose(grploc);
        particleVarPath.clear();
        particleVarPath.str("");
    }

    RDouble *particleAccelerationX = 0,  *particleAccelerationY = 0,  *particleAccelerationZ = 0;
    if (parameter->ifUseParticleAcceleration())
    {
        //! Read particle acceleration.
        particleVarPath << particleGroupPath.str() << "/Acceleration";
        groupNameString = particleVarPath.str();
        grpName = groupNameString;
        grploc = OpenGroup(actkey->filepos, grpName);
        particleAccelerationX = ReadDoubleOneRow(grploc, "AccelerationX");
        particleAccelerationY = ReadDoubleOneRow(grploc, "AccelerationY");
        particleAccelerationZ = ReadDoubleOneRow(grploc, "AccelerationZ");
        H5Gclose(grploc);
        particleVarPath.clear();
        particleVarPath.str("");
    }

    RDouble **particleTemperature = 0;
    if (parameter->ifUseParticleTemperature())
    {
        parameter->SetNumDimTemperature(nDimTem);

        //! Read particle temperature.
        particleVarPath << particleGroupPath.str() << "/Temperature";
        groupNameString = particleVarPath.str();
        grpName = groupNameString;
        grploc = OpenGroup(actkey->filepos, grpName);
        particleTemperature = new RDouble * [nDimTem];
        for (int iTem = 0; iTem < nDimTem; ++iTem)
        {
            ostringstream temperatureOss;
            temperatureOss << "Temperature" << iTem;
            particleTemperature[iTem] = ReadDoubleOneRow(grploc, temperatureOss.str());
        }
        H5Gclose(grploc);
        particleVarPath.clear();
        particleVarPath.str("");
    }

    //! For particle force.
    RDouble *particleGravityX = 0,  *particleGravityY = 0,  *particleGravityZ = 0;
    RDouble *particleDragX = 0,  *particleDragY = 0,  *particleDragZ = 0;
    RDouble *particleSaffmanX = 0,  *particleSaffmanY = 0,  *particleSaffmanZ = 0;
    RDouble *particleMagnusX = 0,  *particleMagnusY = 0,  *particleMagnusZ = 0;
    RDouble *particleAddedMassX = 0,  *particleAddedMassY = 0,  *particleAddedMassZ = 0;
    RDouble *particleFluidAccelerationX = 0,  *particleFluidAccelerationY = 0,  *particleFluidAccelerationZ = 0;
    RDouble *particleBrownianX = 0,  *particleBrownianY = 0,  *particleBrownianZ = 0;
    RDouble *particleThermophoreticX = 0,  *particleThermophoreticY = 0,  *particleThermophoreticZ = 0;
    using namespace PARTICLE_FORCETYPE;
    if (2 == fileType)
    {
        if (parameter->GetForceGravityType() != NO_FORCE)
        {
            particleVarPath << particleGroupPath.str() << "/Force";
            groupNameString = particleVarPath.str();
            grpName = groupNameString;
            grploc = OpenGroup(actkey->filepos, grpName);
            particleGravityX = ReadDoubleOneRow(grploc, "GravityX");
            particleGravityY = ReadDoubleOneRow(grploc, "GravityY");
            particleGravityZ = ReadDoubleOneRow(grploc, "GravityZ");
            H5Gclose(grploc);
            particleVarPath.clear();
            particleVarPath.str("");
        }

        if (parameter->GetForceDragType() != NO_FORCE)
        {
            //! Drag
            particleVarPath << particleGroupPath.str() << "/Force";
            groupNameString = particleVarPath.str();
            grpName = groupNameString;
            grploc = OpenGroup(actkey->filepos, grpName);
            particleDragX = ReadDoubleOneRow(grploc, "DragX");
            particleDragY = ReadDoubleOneRow(grploc, "DragY");
            particleDragZ = ReadDoubleOneRow(grploc, "DragZ");
            H5Gclose(grploc);
            particleVarPath.clear();
            particleVarPath.str("");
        }

        if (parameter->GetForceSaffmanType() != NO_FORCE)
        {
            //! Saffman
            particleVarPath << particleGroupPath.str() << "/Force";
            groupNameString = particleVarPath.str();
            grpName = groupNameString;
            grploc = OpenGroup(actkey->filepos, grpName);
            particleSaffmanX = ReadDoubleOneRow(grploc, "SaffmanX");
            particleSaffmanY = ReadDoubleOneRow(grploc, "SaffmanY");
            particleSaffmanZ = ReadDoubleOneRow(grploc, "SaffmanZ");
            H5Gclose(grploc);
            particleVarPath.clear();
            particleVarPath.str("");
        }

        if (parameter->GetForceMagnusType() != NO_FORCE)
        {
            //! Magnus
            particleVarPath << particleGroupPath.str() << "/Force";
            groupNameString = particleVarPath.str();
            grpName = groupNameString;
            grploc = OpenGroup(actkey->filepos, grpName);
            particleMagnusX = ReadDoubleOneRow(grploc, "MagnusX");
            particleMagnusY = ReadDoubleOneRow(grploc, "MagnusY");
            particleMagnusZ = ReadDoubleOneRow(grploc, "MagnusZ");
            H5Gclose(grploc);
            particleVarPath.clear();
            particleVarPath.str("");
        }

        if (parameter->GetForceAddedMassType() != NO_FORCE)
        {
            //! Added Mass
            particleVarPath << particleGroupPath.str() << "/Force";
            groupNameString = particleVarPath.str();
            grpName = groupNameString;
            grploc = OpenGroup(actkey->filepos, grpName);
            particleAddedMassX = ReadDoubleOneRow(grploc, "AddedMassX");
            particleAddedMassY = ReadDoubleOneRow(grploc, "AddedMassY");
            particleAddedMassZ = ReadDoubleOneRow(grploc, "AddedMassZ");
            H5Gclose(grploc);
            particleVarPath.clear();
            particleVarPath.str("");
        }

        if (parameter->GetFluidAccelerationType() != NO_FORCE)
        {
            //! Fluid acceleration force
            particleVarPath << particleGroupPath.str() << "/Force";
            groupNameString = particleVarPath.str();
            grpName = groupNameString;
            grploc = OpenGroup(actkey->filepos, grpName);
            particleFluidAccelerationX = ReadDoubleOneRow(grploc, "FluidAccelerationX");
            particleFluidAccelerationY = ReadDoubleOneRow(grploc, "FluidAccelerationY");
            particleFluidAccelerationZ = ReadDoubleOneRow(grploc, "FluidAccelerationZ");
            H5Gclose(grploc);
            particleVarPath.clear();
            particleVarPath.str("");
        }

        if (parameter->GetForceBrownianType() != NO_FORCE)
        {
            particleVarPath << particleGroupPath.str() << "/Force";
            groupNameString = particleVarPath.str();
            grpName = groupNameString;
            grploc = OpenGroup(actkey->filepos, grpName);
            particleBrownianX = ReadDoubleOneRow(grploc, "BrownianX");
            particleBrownianY = ReadDoubleOneRow(grploc, "BrownianY");
            particleBrownianZ = ReadDoubleOneRow(grploc, "BrownianZ");
            H5Gclose(grploc);
            particleVarPath.clear();
            particleVarPath.str("");
        }

        if (parameter->GetForceThermophoreticType() != NO_FORCE)
        {
            particleVarPath << particleGroupPath.str() << "/Force";
            groupNameString = particleVarPath.str();
            grpName = groupNameString;
            grploc = OpenGroup(actkey->filepos, grpName);
            particleThermophoreticX = ReadDoubleOneRow(grploc, "ThermophoreticX");
            particleThermophoreticY = ReadDoubleOneRow(grploc, "ThermophoreticY");
            particleThermophoreticZ = ReadDoubleOneRow(grploc, "ThermophoreticZ");
            H5Gclose(grploc);
            particleVarPath.clear();
            particleVarPath.str("");
        }
    }

    //! The number of particle on current zone.
    int nParticleOnCurrentZoneCheck = 0;

    //!Loop for each particle to add.
    for (int iParticle = 0; iParticle < nParticleLocal; ++iParticle)
    {
        int ndim = 3;
        SPDouble onePointCoordinate(ndim);
        onePointCoordinate[0] = particleCoordinateX[iParticle];
        onePointCoordinate[1] = particleCoordinateY[iParticle];
        onePointCoordinate[2] = particleCoordinateZ[iParticle];

        //! if the particle is in current zone.
        //! For example, when particles are at Ghost or Corner 
        //! and need to be transported or removed by MPI or OutBC.
        bool inZone = false;
        //! Note that the search does not include particles 
        //! from the Ghostand corner areas.
        //! Particles in these areas will be accepted into the required zone.
        int idCellOfParticle = 0;
        if (parameter->ifIsOrthogonalForPreAndPost())
        {
            idCellOfParticle = this->SearchParticelForOrthogonalStrGrid(grid, onePointCoordinate, inZone);
        }
        else
        {
            idCellOfParticle = this->QuickSearchParticleCell(grid, onePointCoordinate, inZone);
        }

        if (inZone)
        {
            //! Init the point variable.
            OnePointVariable *onePointVariable = new OnePointVariable();
            onePointVariable->SetParticleID(particleID[iParticle]);
            onePointVariable->SetDiameter(diameter[iParticle]);
            onePointVariable->SetDensity(density[iParticle]);
            onePointVariable->SetSpecificHeatCapacity(specificHeatCapacity[iParticle]);
            onePointVariable->SetEmissivity(particleEmissivity[iParticle]);

            onePointVariable->SetCoordinate(&particleCoordinateX[iParticle], &particleCoordinateY[iParticle], &particleCoordinateZ[iParticle]);
            onePointVariable->SetVelocity(&particleVelocityX[iParticle], &particleVelocityY[iParticle], &particleVelocityZ[iParticle]);

            //! Dimensionless initialization variables.
            onePointVariable->SetDimensionlessLength(refReynoldsLengthDimensional);
            onePointVariable->SetDimensionlessVelocity(refVelocityDimensional);
            onePointVariable->SetDimensionlessDensity(refDensityDimensional);
            onePointVariable->SetDimensionlessSpecificHeatCapacity(refVelocityDimensional, refTemperatureDimensional);

            if (parameter->ifUseParticleRotate())
            {
                onePointVariable->SetAngularVelocity(&particleAngularVelocityX[iParticle], &particleAngularVelocityY[iParticle], &particleAngularVelocityZ[iParticle]);
                onePointVariable->SetDimensionlessAngularVelocity(refAngularVelocityDimensional);
            }

            if (parameter->ifUseParticleAcceleration())
            {
                onePointVariable->SetAcceleration(&particleVelocityX[iParticle], &particleAccelerationY[iParticle], &particleAccelerationZ[iParticle]);
                onePointVariable->SetDimensionlessAcceleration(refAccelerationyDimensional);
            }
            
            if (parameter->ifUseParticleTemperature())
            {
                onePointVariable->SetTemperature(&particleTemperature[nDimTem - 1][iParticle]);
                onePointVariable->SetDimensionlessTemperature(refTemperatureDimensional);
            }
            
            //! Add particle force by restart.
            if (2 == fileType)
            {
                int forceType;

                if (parameter->GetForceGravityType() != NO_FORCE)
                {
                    forceType = 0;
                    onePointVariable->SetForce(&particleGravityX[iParticle], &particleGravityY[iParticle], &particleGravityZ[iParticle], forceType);
                }
                if (parameter->GetForceDragType() != NO_FORCE)
                {
                    forceType = 1;
                    onePointVariable->SetForce(&particleDragX[iParticle], &particleDragY[iParticle], &particleDragZ[iParticle], forceType);
                }

                if (parameter->GetForceSaffmanType() != NO_FORCE)
                {
                    forceType = 2;
                    onePointVariable->SetForce(&particleSaffmanX[iParticle], &particleSaffmanY[iParticle], &particleSaffmanZ[iParticle], forceType);
                }

                if (parameter->GetForceMagnusType() != NO_FORCE)
                {
                    forceType = 3;
                    onePointVariable->SetForce(&particleMagnusX[iParticle], &particleMagnusY[iParticle], &particleMagnusZ[iParticle], forceType);
                }

                if (parameter->GetForceAddedMassType() != NO_FORCE)
                {
                    forceType = 4;
                    onePointVariable->SetForce(&particleAddedMassX[iParticle], &particleAddedMassY[iParticle], &particleAddedMassZ[iParticle], forceType);
                }

                if (parameter->GetFluidAccelerationType() != NO_FORCE)
                {
                    forceType = 5;
                    onePointVariable->SetForce(&particleFluidAccelerationX[iParticle], &particleFluidAccelerationY[iParticle], &particleFluidAccelerationZ[iParticle], forceType);
                }

                if (parameter->GetForceBrownianType() != NO_FORCE)
                {
                    forceType = 6;
                    onePointVariable->SetForce(&particleBrownianX[iParticle], &particleBrownianY[iParticle], &particleBrownianZ[iParticle], forceType);
                }

                if (parameter->GetForceThermophoreticType() != NO_FORCE)
                {
                    forceType = 7;
                    onePointVariable->SetForce(&particleThermophoreticX[iParticle], &particleThermophoreticY[iParticle], &particleThermophoreticZ[iParticle], forceType);
                }

                onePointVariable->SetDimensionlessForce(refForceDimensional);
            }

            //! Set particle cellID and add particle to particleParticleGroup (only inZone).
            particlePointGroup->QuickAddParticle(particleID[iParticle], idCellOfParticle, onePointVariable, inZone);

            //! Is particle to IndexCell.
            //! Here it means that the particles 
            //! are in one layer(inner) of cells close to the wall, 
            //! so it is necessary to decide whether to use boundary conditions.
            bool ifParticleInIndexCell = (particlePointGroup->GetIndexCell())->SetParticleInCellInfo(grid, onePointVariable);

            if (ifParticleInIndexCell)
            {
                bool ifRemoveParticle = ParticleBoundary(grid,onePointVariable, particlePointGroup->GetIndexCell());

                if (ifRemoveParticle)
                {
                    ostringstream ossError;
                    ossError << " Error remove particle when init but inZone is true. " << "\n";
                    ossError << " subact : " << actkey->subact << "\n";
                    ossError << " zoneIDGlobal : " << zoneIDGlobal << " zoneIDLocal : " << zoneIDLocal << "\n";
                    ossError << " Particle ID: " << particleID[iParticle] << "\n";
                    ossError << (onePointVariable->Print2WindowOS()) << endl;
                    PrintToWindow(ossError);
                }
            }
            nParticleOnCurrentZoneCheck += 1;
        }
    }

    if (0 == actkey->subact)
    {
        if (0 == fileType)
        {
            nParticleLocal = nParticleOnCurrentZoneCheck;
            parameter->SetNumParticleLocal(nParticleLocal);
        }
    }

    ossReadParticle << " --- Check the particle on each zone ---" << "\n";
    ossReadParticle << " subact : " << actkey->subact << "\n";
    ossReadParticle << " zoneIDGlobal : " << zoneIDGlobal << " zoneIDLocal : " << zoneIDLocal << "\n";
    ossReadParticle << " currentProcessorID : " << currentProcessorID << "\n";
    ossReadParticle << " nParticle local number by Sum : " << nParticleOnCurrentZoneCheck << "\n";
    ossReadParticle << " nParticle local number by HDF5  : " << nParticleLocal << "\n";
    ossReadParticle << " nParticle total number by HDF5  : " << nParticleTotal << endl;
    WriteLogFile(ossReadParticle, writeEachProcessor);
    ossReadParticle.clear();
    ossReadParticle.str("");

    if (nParticleOnCurrentZoneCheck != nParticleLocal)
    {
        ostringstream ossError;
        ossError << " Error in local particle number in read particle file " << "\n";
        ossError << " nParticle local number by Sum : " << nParticleOnCurrentZoneCheck << "\n";
        ossError << " nParticle local number by HDF5  : " << nParticleLocal << "\n";
        ossError << endl;
        TK_Exit::ExceptionExit(ossError);
    }

    DelPointer(particleID);

    //! 3unit
    DelPointer(diameter);
    DelPointer(density);
    DelPointer(specificHeatCapacity);

    //! Read particle coodinate.
    DelPointer(particleCoordinateX);
    DelPointer(particleCoordinateY);
    DelPointer(particleCoordinateZ);
    //! Read particle velocity.
    DelPointer(particleVelocityX);
    DelPointer(particleVelocityY);
    DelPointer(particleVelocityZ);

     //! Read particle Angular velocity.
    DelPointer(particleAngularVelocityX);
    DelPointer(particleAngularVelocityY);
    DelPointer(particleAngularVelocityZ);

    //! Read particle acceleration.
    DelPointer(particleAccelerationX);
    DelPointer(particleAccelerationY);
    DelPointer(particleAccelerationZ);

    //! Read particle acceleration.
    DelPointer2(particleTemperature);
    //! For particle force.
    DelPointer(particleGravityX);
    DelPointer(particleGravityY);
    DelPointer(particleGravityZ);

    DelPointer(particleDragX);
    DelPointer(particleDragY);
    DelPointer(particleDragZ);

    DelPointer(particleSaffmanX);
    DelPointer(particleSaffmanY);
    DelPointer(particleSaffmanZ);

    DelPointer(particleMagnusX);
    DelPointer(particleMagnusX);
    DelPointer(particleMagnusZ);

    DelPointer(particleAddedMassX);
    DelPointer(particleAddedMassY);
    DelPointer(particleAddedMassZ);

    DelPointer(particleFluidAccelerationX);
    DelPointer(particleFluidAccelerationY);
    DelPointer(particleFluidAccelerationZ);

    DelPointer(particleBrownianX);
    DelPointer(particleBrownianY);
    DelPointer(particleBrownianZ);

    DelPointer(particleThermophoreticX);
    DelPointer(particleThermophoreticY);
    DelPointer(particleThermophoreticZ);
}

int ParticlePointSolverParallelDomain::QuickSearchParticleCell(Grid *grid, SPDouble &onePointCoordinate, bool &inZone)
{
    IndexCell indexCell;
    int IDCellOfParticle = indexCell.SearchParticleCell(grid, &onePointCoordinate, inZone);
    return IDCellOfParticle;
}

int ParticlePointSolverParallelDomain::SearchParticelForOrthogonalStrGrid(Grid *gridIn, SPDouble &onePointCoordinate, bool &inZone)
{
    int IDCellOfParticle = 0;

    int gridType = gridIn->Type();
    if (gridType == PHSPACE::STRUCTGRID)
    {
        StructGrid *grid = StructGridCast(gridIn);

        int nLayers = 0;
        //! Index start and end of cell in current zone.
        int ist, ied;
        int jst, jed;
        int kst, ked;
        grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, nLayers);

        //! Index start and end of surface in .
        //! From 1 to nDim.
        int nDim = GetDim();

        RDouble3D &structX = *(grid->GetStructX());
        RDouble3D &structY = *(grid->GetStructY());
        RDouble3D &structZ = *(grid->GetStructZ());

        //! The bool for particle in current cell.
        //! 1 -- the particle is in current zone.
        //! 0 -- else.
        bool isParticleInCell = false;
        bool isParticleInCellx = false;
        bool isParticleInCelly = false;
        bool isParticleInCellz = false;

        int iP=0, jP=0, kP=0;

        //! First, loop for index of cell in struct grid.
        for (int k = kst; k <= ked; ++k)
        {
            for (int j = jst; j <= jed; ++j)
            {
                if (!isParticleInCellx)
                {
                    for (int i = ist; i <= ied; ++i)
                    {
                        if (structX(i, j, k) <= onePointCoordinate[0] && structX(i + 1, j, k) >= onePointCoordinate[0])
                        {
                            iP = i;
                            isParticleInCellx = true;
                            break;
                        }
                    }
                }

                if (!isParticleInCelly)
                {
                    if (structY(iP, j, k) <= onePointCoordinate[1] && structY(iP, j + 1, k) >= onePointCoordinate[1])
                    {
                        jP = j;
                        isParticleInCelly = true;
                        break;
                    }
                }
            }

            if (!isParticleInCellz)
            {
                if (GetDim() == TWO_D)
                {
                    kP = k;
                    isParticleInCellz = true;
                    break;
                }
                else if (GetDim() == THREE_D)
                {
                    if (structZ(iP, jP, k) <= onePointCoordinate[2] && structZ(iP, jP, k + 1) >= onePointCoordinate[2])
                    {
                        kP = k;
                        isParticleInCellz = true;
                        break;
                    }
                }
            }
        }

        if (isParticleInCellx && isParticleInCelly && isParticleInCellz)
        {
            isParticleInCell = true;
        }

        //! Get index IJK of particle position in current zone.
        int nI, nJ, nK;
        grid->GetND(nI, nJ, nK);
        nI -= 1;
        nJ -= 1;
        nK -= 1;
        if (GetDim() == TWO_D)
        {
            nK = 1;
        }

        GetIDCellOfParticleByIJKSturct(nI, nJ, nK, iP, jP, kP, GetNumberOfGhostCellLayers(), IDCellOfParticle);

        if (isParticleInCell)
        {
            inZone = true;
        }
        else
        {
            inZone = false;
        }

    }
    else
    {
        ostringstream ossError;
        ossError << " Error on gridTypeForPreAndPost. " << "\n";
        ossError << " the gridType is not struct grid " << "\n";
        ossError << " but gridTypeForPreAndPost = 1" << "\n";
        ossError << endl;
    }
    return IDCellOfParticle;
}

void ParticlePointSolverParallelDomain::CheckNumOfParticle()
{
    using namespace PHMPI;
    int level = 0;
    Grid *grid = GetGrid(level);

    //! The global zone id on current zone on current processer.
    int zoneIDGlobal = grid->GetZoneID();
    //! The number of zones on global all processer.
    int nZonesGlobal = GetNumberofGlobalZones();

    //! The local zone id on current zone on current processer.
    int zoneIDLocal = grid->GetZoneLocalID();
    //! The number of zones on on current processer.
    int nZonesLocal = GetNumberofLocalZones();

    if (zoneIDLocal != 0)
    {
        return;
    }

    //! Get current processor ID on current zones..
    int currentProcessorID = GetCurrentProcessorID();

    //! The ProcessorID of current zone.
    int currentZoneProcessorID = GetZoneProcessorID(zoneIDGlobal);

    //! The number of particle on current Region(ProcessorID) by sum.
    int nPartucleLocalRegion = 0;
    //! The number of particle on all Region(ProcessorID) by sum.
    int nParticleTotalSum = 0;

    Param_ParticlePointSolverParallelDomain *parameter = GetControlParameters();
    //! The number of particle on all Region(ProcessorID) in Param.
    int nParticleTotal = parameter->GetNumParticleTotal();

    //! Sum the particle on on current Region(ProcessorID).
    for (int iZoneGlobal = 0; iZoneGlobal < nZonesGlobal; ++iZoneGlobal)
    {
        //! Judge the zones on current processor.
        int recvProcess = GetZoneProcessorID(iZoneGlobal);

        //! The loop for zone on current processor.
        if (currentProcessorID == recvProcess)
        {
            //! For each zone'solver.
            int solverID = this->GetIndex();
            ParticlePointSolverParallelDomain *solver = static_cast <ParticlePointSolverParallelDomain*> (GlobalSolvers::GetSolver(iZoneGlobal, solverID));
            Param_ParticlePointSolverParallelDomain *parameter = solver->GetControlParameters();
            ParticlePointGroup *particlePointGroup = solver->GetParticlePointGroup("particlePointGroup");

            //! The number of particle on each zone.
            int nParticleLocalZone = particlePointGroup->GetNumOfLocalParticle();

            if (parameter->GetNumParticleLocal() != nParticleLocalZone)
            {
                ostringstream ossError;
                ossError << " Error in local particle number in check number" << "\n";
                ossError << " nParticle local number by Sum : " << nParticleLocalZone << "\n";
                ossError << " nParticle local number by HDF5  : " << parameter->GetNumParticleLocal() << "\n";
                ossError << endl;
                TK_Exit::ExceptionExit(ossError);
            }

            nPartucleLocalRegion += nParticleLocalZone;
        }
    }
;
    PH_Barrier();
    PH_AllReduce(&nPartucleLocalRegion, &nParticleTotalSum, 1, MPI_SUM);

    //! Print messege for read particle.
    ostringstream ossLog;
    ossLog << "--- Check particle number on all region ---" << "\n";
    ossLog << " zoneIDGlobal : " << zoneIDGlobal << " zoneIDLocal : " << zoneIDLocal << "\n";
    ossLog << " nParticle by HDF5 :" << nParticleTotal << "\n";
    ossLog << " nParticle by sum  :" << nParticleTotalSum << endl;
    bool writeEachProcessor = true;

    WriteLogFile(ossLog, writeEachProcessor);
    PrintToWindow(ossLog);

    if (nParticleTotal != nParticleTotalSum)
    {
        ostringstream oss;
        oss << " zoneIDGlobal : " << zoneIDGlobal << " zoneIDLocal : " << zoneIDLocal << "\n";
        oss << " Error nParticle by parameter : " << " nParticle = " << nParticleTotal << "\n";
        oss << " nParticleCheck by sum : " << " nParticleCheck = " << nParticleTotalSum << endl;
        TK_Exit::ExceptionExit(oss);
    }
}

void ParticlePointSolverParallelDomain::UpdateNumOfParticle()
{
    using namespace PHMPI;
    int level = 0;
    Grid *grid = GetGrid(level);

    //! The global zone id on current zone on current processer.
    int zoneIDGlobal = grid->GetZoneID();
    //! The number of zones on global all processer.
    int nZonesGlobal = GetNumberofGlobalZones();

    //! The local zone id on current zone on current processer.
    int zoneIDLocal = grid->GetZoneLocalID();
    //! The number of zones on on current processer.
    int nZonesLocal = GetNumberofLocalZones();

    if (zoneIDLocal != 0)
    {
        return;
    }

    //! Get current processor ID on current zones..
    int currentProcessorID = GetCurrentProcessorID();

    //! The ProcessorID of current zone.
    int currentZoneProcessorID = GetZoneProcessorID(zoneIDGlobal);

    //! The number of particle on current Region(ProcessorID) by sum.
    int nPartucleLocalRegion = 0;
    //! The number of particle on all Region(ProcessorID) by sum.
    int nParticleTotalSum = 0;

    Param_ParticlePointSolverParallelDomain *parameter = GetControlParameters();
    //! The number of particle on all Region(ProcessorID) in Param.
    int nParticleTotal = parameter->GetNumParticleTotal();

    //! Sum the particle on on current Region(ProcessorID).
    for (int iZoneGlobal = 0; iZoneGlobal < nZonesGlobal; ++iZoneGlobal)
    {
        //! Judge the zones on current processor.
        int recvProcess = GetZoneProcessorID(iZoneGlobal);

        //! The loop for zone on current processor.
        if (currentProcessorID == recvProcess)
        {
            //! For each zone'solver.
            int solverID = this->GetIndex();
            ParticlePointSolverParallelDomain *solver = static_cast <ParticlePointSolverParallelDomain*> (GlobalSolvers::GetSolver(iZoneGlobal, solverID));
            Param_ParticlePointSolverParallelDomain *parameter = solver->GetControlParameters();
            ParticlePointGroup *particlePointGroup = solver->GetParticlePointGroup("particlePointGroup");

            //! The number of particle on each zone.
            int nParticleLocalZone = particlePointGroup->GetNumOfLocalParticle();

            nPartucleLocalRegion += nParticleLocalZone;
        }
    }

    PH_Barrier();
    PH_AllReduce(&nPartucleLocalRegion, &nParticleTotalSum, 1, MPI_SUM);

    //! Print messege for read particle.
    ostringstream ossLog;
    ossLog << "--- Update particle number on all region ---" << "\n";
    ossLog << " zoneIDGlobal : " << zoneIDGlobal << " zoneIDLocal : " << zoneIDLocal << "\n";
    ossLog << " old nParticle Total :" << nParticleTotal << "\n";
    ossLog << " new nParticle by sum  :" << nParticleTotalSum << endl;
    bool writeEachProcessor = true;

    if (nParticleTotal != nParticleTotalSum)
    {
        WriteLogFile(ossLog, writeEachProcessor);
        PrintToWindow(ossLog);
    }

    nParticleTotal = nParticleTotalSum;
    parameter->SetNumParticleTotal(nParticleTotal);
}

void ParticlePointSolverParallelDomain::CheckParticleCell()
{
    int level = 0;
    Grid *grid = GetGrid(level);
    Param_ParticlePointSolverParallelDomain *param_ParticlePointSolverParallelDomain = GetControlParameters();
    ParticlePointGroup *particlePointGroup = reinterpret_cast<ParticlePointGroup*> (grid->GetDataPtr("particlePointGroup"));

    //! The number of particle on current zone.
    int nLocalParticle = particlePointGroup->GetNumOfLocalParticle();

    if (nLocalParticle == 0)
    {
        return;
    }

    //! The global zone id.
    int zoneIDGlobal = grid->GetZoneID();

    for (int iParticle = 0; iParticle < nLocalParticle; ++iParticle)
    {
        int particleID = particlePointGroup->GetGlobalParticleIDByLocalIParticle(iParticle);
        int cellID = particlePointGroup->GetCellIDOfParticle(particleID);
        bool inZone = particlePointGroup->isParticleInZone(particleID);
        if (!inZone)
        {
            cout << "particleID " << particleID << endl;
            TK_Exit::UnexpectedVarValue(" cellID : ", cellID);
        }
        else
        {
            ostringstream oss;
            oss << " iParticle = " << iParticle << "\n";
            oss << " Particle ID = " << particleID << " \n";
            oss << " cell ID = " << cellID << endl;
            cout << oss.str();
        }
    }

    int nCell = particlePointGroup->GetNumOfLocalCell();
    int checkParticleOnCell = 0;
    for (int iCell = 0; iCell < nCell; ++iCell)
    {
        int nParticleOfCell = particlePointGroup->GetNumOfParticleOfCell(iCell);
        if (nParticleOfCell != 0)
        {
            int *particleIDList = particlePointGroup->GetIDOfParticleOfCell(iCell);
            ostringstream oss;
            oss << " iCell = " << iCell << "\n";
            for (int iParticleOfCell = 0; iParticleOfCell < nParticleOfCell; ++iParticleOfCell)
            {
                int particleID = particleIDList[iParticleOfCell];
                oss << "---particleID = " << particleID << "\n";
            }
        }
        checkParticleOnCell += nParticleOfCell;
    }

    if (checkParticleOnCell != nLocalParticle)
    {
        TK_Exit::UnexpectedVarValue(" checkParticleOnCell isn't equal to nParticleOfCell :", checkParticleOnCell);
    }
    else
    {
        ostringstream oss;
        oss << " nLocalParticle = " << nLocalParticle << " checkParticleOnCell = " << checkParticleOnCell << endl;
        cout << oss.str();
    }
}

void ParticlePointSolverParallelDomain::DumpResultFile(Grid *grid, int level)
{
    int innerIterStep = PHSPACE::GlobalDataBase::GetIntParaFromDB("innstep");
    if (innerIterStep != 1)
    {
        return;
    }

    Param_ParticlePointSolverParallelDomain *param_ParticlePointSolverParallelDomain = GetControlParameters();

    //! The step intervals for flow.
    int intervalStepFlow = param_ParticlePointSolverParallelDomain->GetIntervalStepFlow();
    //! The step intervals for particle.
    int intervalStepRestartParticle = param_ParticlePointSolverParallelDomain->GetIntervalStepRestartParticle();

    int outnstep = GlobalDataBase::GetIntParaFromDB("outnstep");

    if ((1 == outnstep || 0 == outnstep % intervalStepRestartParticle))
    {
        //! To dump particle data.
        ActionKey *actkeyDumpRestartData = new ActionKey();
        FillActionKey(actkeyDumpRestartData, DUMP_RESTART, level);
        DumpRestartData(actkeyDumpRestartData);
        delete actkeyDumpRestartData;
    }
}

void ParticlePointSolverParallelDomain::DumpRestartData(ActionKey *actkey)
{
    using namespace PHMPI;
    Param_ParticlePointSolverParallelDomain *parameter = GetControlParameters();

    int zoneIDofThisSolver = GetZoneLocalID();
    if (zoneIDofThisSolver != 0)
    {
        //! Very Important: only the first zone on a processor need to implement the task.
        //! other wise, the tasks would be implemented several times.
        return;
    }

    if (!JudgeIfWriteRestart())
    {
        return;
    }

    actkey->filename = "./results/particle.h5";

    int outnstep = GlobalDataBase::GetIntParaFromDB("outnstep");
    int initOutStep = parameter->GetInitOutStep();
    outnstep += initOutStep;
    actkey->filename = PHSPACE::AddSymbolToFileName(actkey->filename, "_", outnstep);

    hid_t file;
    file = CreateHDF5File(actkey->filename);
    actkey->filepos = file;

    CreateH5ParticleFile(actkey);

    int nZones = GetNumberofGlobalZones();
    int currentProcessorID = GetCurrentProcessorID();

    for (int iZone = 0; iZone < nZones; ++iZone)
    {
        int sendProcessID = GetZoneProcessorID(iZone);
        if (currentProcessorID == sendProcessID)
        {
            ParticlePointSolverParallelDomain *solver = static_cast <ParticlePointSolverParallelDomain*> (GlobalSolvers::GetSolver(iZone, actkey->solverID));
            solver->DumpRestartH5(actkey);
        }
    }

    H5Fclose(file);
    actkey->filepos = 0;
}

void ParticlePointSolverParallelDomain::CreateH5ParticleFile(ActionKey *actkey)
{
    //! This function is similar to CFDSolver::CreateH5RestartFile.
    //! But CFDSolver::CreateH5RestartFile is not a virtual function.
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();
    int currentProcessorID = GetCurrentProcessorID();

    hid_t grploc,grpVar;
    string grpName;
    ostringstream ossData;

    //! Basic information of file.
    //! Version
    grploc = GetGroup(actkey->filepos, "Version");
    CreateEmptyData(grploc, "version", 1, 1, PHDOUBLE);
    H5Gclose(grploc);
    //! nGroup.
    grploc = GetGroup(actkey->filepos, "nGroup");
    CreateEmptyData(grploc, "nGroup", 1, 1, PHINT);
    CreateEmptyData(grploc, "nParticleTotal", 1, 1, PHINT);
    H5Gclose(grploc);
    //! Time
    grploc = GetGroup(actkey->filepos, "Time");
    CreateEmptyData(grploc, "outnstep", 1, 1, PHINT);
    H5Gclose(grploc);

    //! Flow param.
    grploc = GetGroup(actkey->filepos, "FlowParam");
    CreateEmptyData(grploc, "refReynoldsLengthDimensional", 1, 1, PHDOUBLE);
    CreateEmptyData(grploc, "refVelocityDimensional", 1, 1, PHDOUBLE);
    CreateEmptyData(grploc, "refTemperatureDimensional", 1, 1, PHDOUBLE);
    CreateEmptyData(grploc, "refPressureDimensional", 1, 1, PHDOUBLE);
    CreateEmptyData(grploc, "refDensityDimensional", 1, 1, PHDOUBLE);
    CreateEmptyData(grploc, "refAverageGeneralGasConstantDimensional", 1, 1, PHDOUBLE);
    CreateEmptyData(grploc, "refDynamicViscosityDimensional", 1, 1, PHDOUBLE);

    CreateEmptyData(grploc, "refMachNumber", 1, 1, PHDOUBLE);
    CreateEmptyData(grploc, "refReNumber", 1, 1, PHDOUBLE);
    CreateEmptyData(grploc, "refGama", 1, 1, PHDOUBLE);

    CreateEmptyData(grploc, "refMassDimensional", 1, 1, PHDOUBLE);
    CreateEmptyData(grploc, "refTimeDimensional", 1, 1, PHDOUBLE);
    CreateEmptyData(grploc, "refAngularVelocityDimensional", 1, 1, PHDOUBLE);
    CreateEmptyData(grploc, "refAccelerationyDimensional", 1, 1, PHDOUBLE);
    CreateEmptyData(grploc, "refForceDimensional", 1, 1, PHDOUBLE);

    CreateEmptyData(grploc, "physicalTimeStep", 1, 1, PHDOUBLE);

    H5Gclose(grploc);

    for (int iZone = 0; iZone < nZones; ++iZone)
    {
        if (currentProcessorID == GetZoneProcessorID(iZone))
        {
            //! Group name.
            ostringstream ossFolder;
            ossFolder << "Particle_Group" << iZone;
            grpName = ossFolder.str();
            ossFolder.clear();
            ossFolder.str("");

            CreateGroup(actkey->filepos, grpName);
            grploc = OpenGroup(actkey->filepos, grpName);

            ParticlePointSolverParallelDomain *solver = static_cast <ParticlePointSolverParallelDomain*> (GlobalSolvers::GetSolver(iZone, actkey->solverID));
            Param_ParticlePointSolverParallelDomain *parameter = solver->GetControlParameters();

            int nLocalParticle = solver->GetNumOfLocalParticle();

            //! Info
            ossData << "Info";
            grpName = ossData.str();
            grpVar = GetGroup(grploc, grpName);
            CreateEmptyData(grpVar, "ParticleID", 1, nLocalParticle, PHINT);
            CreateEmptyData(grpVar, "nDimTem", 1, 1, PHINT);
            CreateEmptyData(grpVar, "nParticle", 1, 1, PHINT);
            H5Gclose(grpVar);
            ossData.clear();
            ossData.str("");

            if (0 == nLocalParticle)
            {
                H5Gclose(grploc);
                continue;
            }

            //! BasicParam.
            ossData << "BasicParam";
            grpName = ossData.str();
            grpVar = GetGroup(grploc, grpName);
            CreateEmptyData(grpVar, "Density", 1, nLocalParticle, PHDOUBLE);
            CreateEmptyData(grpVar, "Diameter", 1, nLocalParticle, PHDOUBLE);
            CreateEmptyData(grpVar, "particleSpecificHeatCapacity", 1, nLocalParticle, PHDOUBLE);
            CreateEmptyData(grpVar, "ParticleReNumber", 1, nLocalParticle, PHDOUBLE);
            H5Gclose(grpVar);
            ossData.clear();
            ossData.str("");

            //! Heat transfer.
            ossData << "HeatTransfer";
            grpName = ossData.str();
            grpVar = GetGroup(grploc, grpName);
            CreateEmptyData(grpVar, "ConvectiveHeatTransfer", 1, nLocalParticle, PHDOUBLE);
            CreateEmptyData(grpVar, "RadiativeHeatTransfer", 1, nLocalParticle, PHDOUBLE);
            CreateEmptyData(grpVar, "ParticleWork", 1, nLocalParticle, PHDOUBLE);
            CreateEmptyData(grpVar, "ParticleDissipation", 1, nLocalParticle, PHDOUBLE);
            CreateEmptyData(grpVar, "NuseltNumber", 1, nLocalParticle, PHDOUBLE);
            H5Gclose(grpVar);
            ossData.clear();
            ossData.str("");

            //! Coordinate
            ossData << "Coordinate";
            grpName = ossData.str();
            grpVar = GetGroup(grploc, grpName);
            CreateEmptyData(grpVar, "CoordinateX", 1, nLocalParticle, PHDOUBLE);
            CreateEmptyData(grpVar, "CoordinateY", 1, nLocalParticle, PHDOUBLE);
            CreateEmptyData(grpVar, "CoordinateZ", 1, nLocalParticle, PHDOUBLE);
            H5Gclose(grpVar);
            ossData.clear();
            ossData.str("");

            //! Stokes Number
            ossData << "StokesNumber";
            grpName = ossData.str();
            grpVar = GetGroup(grploc, grpName);
            CreateEmptyData(grpVar, "TotalForceSt", 1, nLocalParticle, PHDOUBLE);
            CreateEmptyData(grpVar, "StokesSt", 1, nLocalParticle, PHDOUBLE);
            CreateEmptyData(grpVar, "ModifiedSt", 1, nLocalParticle, PHDOUBLE);
            H5Gclose(grpVar);
            ossData.clear();
            ossData.str("");

            //! Velocity
            ossData << "Velocity";
            grpName = ossData.str();
            grpVar = GetGroup(grploc, grpName);
            CreateEmptyData(grpVar, "VelocityX", 1, nLocalParticle, PHDOUBLE);
            CreateEmptyData(grpVar, "VelocityY", 1, nLocalParticle, PHDOUBLE);
            CreateEmptyData(grpVar, "VelocityZ", 1, nLocalParticle, PHDOUBLE);
            H5Gclose(grpVar);
            ossData.clear();
            ossData.str("");

            //! Relative Velocity
            ossData << "RelativeVelocity";
            grpName = ossData.str();
            grpVar = GetGroup(grploc, grpName);
            CreateEmptyData(grpVar, "ParticleMaNumber", 1, nLocalParticle, PHDOUBLE);
            CreateEmptyData(grpVar, "RelativeVelocityX", 1, nLocalParticle, PHDOUBLE);
            CreateEmptyData(grpVar, "RelativeVelocityY", 1, nLocalParticle, PHDOUBLE);
            CreateEmptyData(grpVar, "RelativeVelocityZ", 1, nLocalParticle, PHDOUBLE);
            H5Gclose(grpVar);
            ossData.clear();
            ossData.str("");

            //! AngularVelocity
            ossData << "AngularVelocity";
            grpName = ossData.str();
            grpVar = GetGroup(grploc, grpName);
            CreateEmptyData(grpVar, "AngularVelocityX", 1, nLocalParticle, PHDOUBLE);
            CreateEmptyData(grpVar, "AngularVelocityY", 1, nLocalParticle, PHDOUBLE);
            CreateEmptyData(grpVar, "AngularVelocityZ", 1, nLocalParticle, PHDOUBLE);
            H5Gclose(grpVar);
            ossData.clear();
            ossData.str("");

            //! Acceleration
            ossData << "Acceleration";
            grpName = ossData.str();
            grpVar = GetGroup(grploc, grpName);
            CreateEmptyData(grpVar, "AccelerationX", 1, nLocalParticle, PHDOUBLE);
            CreateEmptyData(grpVar, "AccelerationY", 1, nLocalParticle, PHDOUBLE);
            CreateEmptyData(grpVar, "AccelerationZ", 1, nLocalParticle, PHDOUBLE);
            H5Gclose(grpVar);
            ossData.clear();
            ossData.str("");

            //! Temperature
            ossData << "Temperature";
            grpName = ossData.str();
            grpVar = GetGroup(grploc, grpName);
            CreateEmptyData(grpVar, "Temperature0", 1, nLocalParticle, PHDOUBLE);
            H5Gclose(grpVar);
            ossData.clear();
            ossData.str("");

            //! Force.
            ossData << "Force";
            grpName = ossData.str();
            grpVar = GetGroup(grploc, grpName);
            //! Gravity
            CreateEmptyData(grpVar, "GravityX", 1, nLocalParticle, PHDOUBLE);
            CreateEmptyData(grpVar, "GravityY", 1, nLocalParticle, PHDOUBLE);
            CreateEmptyData(grpVar, "GravityZ", 1, nLocalParticle, PHDOUBLE);
            //! Drag
            CreateEmptyData(grpVar, "DragX", 1, nLocalParticle, PHDOUBLE);
            CreateEmptyData(grpVar, "DragY", 1, nLocalParticle, PHDOUBLE);
            CreateEmptyData(grpVar, "DragZ", 1, nLocalParticle, PHDOUBLE);
            //! Saffman
            CreateEmptyData(grpVar, "SaffmanX", 1, nLocalParticle, PHDOUBLE);
            CreateEmptyData(grpVar, "SaffmanY", 1, nLocalParticle, PHDOUBLE);
            CreateEmptyData(grpVar, "SaffmanZ", 1, nLocalParticle, PHDOUBLE);
            //! Magnus
            CreateEmptyData(grpVar, "MagnusX", 1, nLocalParticle, PHDOUBLE);
            CreateEmptyData(grpVar, "MagnusY", 1, nLocalParticle, PHDOUBLE);
            CreateEmptyData(grpVar, "MagnusZ", 1, nLocalParticle, PHDOUBLE);
            //! AddedMassX
            CreateEmptyData(grpVar, "AddedMassX", 1, nLocalParticle, PHDOUBLE);
            CreateEmptyData(grpVar, "AddedMassY", 1, nLocalParticle, PHDOUBLE);
            CreateEmptyData(grpVar, "AddedMassZ", 1, nLocalParticle, PHDOUBLE);
            //! FluidAccelerationX
            CreateEmptyData(grpVar, "FluidAccelerationX", 1, nLocalParticle, PHDOUBLE);
            CreateEmptyData(grpVar, "FluidAccelerationY", 1, nLocalParticle, PHDOUBLE);
            CreateEmptyData(grpVar, "FluidAccelerationZ", 1, nLocalParticle, PHDOUBLE);
            //! BrownianX
            CreateEmptyData(grpVar, "BrownianX", 1, nLocalParticle, PHDOUBLE);
            CreateEmptyData(grpVar, "BrownianY", 1, nLocalParticle, PHDOUBLE);
            CreateEmptyData(grpVar, "BrownianZ", 1, nLocalParticle, PHDOUBLE);
            //! ThermophoreticX
            CreateEmptyData(grpVar, "ThermophoreticX", 1, nLocalParticle, PHDOUBLE);
            CreateEmptyData(grpVar, "ThermophoreticY", 1, nLocalParticle, PHDOUBLE);
            CreateEmptyData(grpVar, "ThermophoreticZ", 1, nLocalParticle, PHDOUBLE);
            H5Gclose(grpVar);
            ossData.clear();
            ossData.str("");

            H5Gclose(grploc);
        }
    }
}

void ParticlePointSolverParallelDomain::DumpRestartH5(ActionKey *actkey)
{
    Grid *grid = GetGrid(actkey->level);
    int gridID = grid->GetZoneID();

    Param_ParticlePointSolverParallelDomain *parameter = GetControlParameters();
    ParticlePointGroup *particlePointGroup = GetParticlePointGroup("particlePointGroup");

    hid_t grploc, grpVar;
    string grpName;

    RDouble refReynoldsLengthDimensional = parameter->GetRefDimLength();
    RDouble refVelocityDimensional = parameter->GetRefDimVelocity();
    RDouble refTemperatureDimensional = parameter->GetRefDimensionalTemperature();
    RDouble refPressureDimensional = parameter->GetRefDimPressure();
    RDouble refDensityDimensional = parameter->GetRefDimDensity();
    RDouble refAverageGeneralGasConstantDimensional = parameter->GetRefDimR();
    RDouble refDynamicViscosityDimensional = parameter->GetRefDimDynamicViscosity();

    RDouble refMachNumber = parameter->GetRefMach();
    RDouble refReNumber = parameter->GetRefReIn();
    RDouble refGama = parameter->GetRefGamma();

    RDouble refMassDimensional = parameter->GetRefDimMass();
    RDouble refTimeDimensional = parameter->GetRefDimTime();
    RDouble refAngularVelocityDimensional = parameter->GetRefDimAngularVelocity();
    RDouble refAccelerationyDimensional = parameter->GetRefDimAcceleration();
    RDouble refForceDimensional = parameter->GetRefDimForce();

    RDouble physicalTimeStep = GlobalDataBase::GetDoubleParaFromDB("physicalTimeStep");

    using namespace PHMPI;
    int currentProcessorID = GetCurrentProcessorID();
    if (currentProcessorID == GetServerProcessorID())
    {
        //! Version.
        grpVar = OpenGroup(actkey->filepos, "Version");
        RDouble version = parameter->GetVersion();
        WriteData(grpVar, &version, "version");
        H5Gclose(grpVar);

        //! Group.
        grpVar = OpenGroup(actkey->filepos, "nGroup");
        int nGroup = parameter->GetNumParticleGroup();
        int nParticleTotal = parameter->GetNumParticleTotal();
        WriteData(grpVar, &nGroup, "nGroup");
        WriteData(grpVar, &nParticleTotal, "nParticleTotal");
        H5Gclose(grpVar);

        //! Time
        grpVar = OpenGroup(actkey->filepos, "Time");
        int outIterStep = GlobalDataBase::GetIntParaFromDB("outnstep");
        int initOutStep = parameter->GetInitOutStep();
        outIterStep += initOutStep;
        WriteData(grpVar, &outIterStep, "outnstep");
        H5Gclose(grpVar);

        grpVar = OpenGroup(actkey->filepos, "FlowParam");

        WriteData(grpVar, &refReynoldsLengthDimensional, "refReynoldsLengthDimensional");
        WriteData(grpVar, &refVelocityDimensional, "refVelocityDimensional");
        WriteData(grpVar, &refTemperatureDimensional, "refTemperatureDimensional");
        WriteData(grpVar, &refPressureDimensional, "refPressureDimensional");
        WriteData(grpVar, &refDensityDimensional, "refDensityDimensional");
        WriteData(grpVar, &refAverageGeneralGasConstantDimensional, "refAverageGeneralGasConstantDimensional");
        WriteData(grpVar, &refDynamicViscosityDimensional, "refDynamicViscosityDimensional");

        WriteData(grpVar, &refMachNumber, "refMachNumber");
        WriteData(grpVar, &refReNumber, "refReNumber");
        WriteData(grpVar, &refGama, "refGama");

        WriteData(grpVar, &refMassDimensional, "refMassDimensional");
        WriteData(grpVar, &refTimeDimensional, "refTimeDimensional");
        WriteData(grpVar, &refAngularVelocityDimensional, "refAngularVelocityDimensional");
        WriteData(grpVar, &refAccelerationyDimensional, "refAccelerationyDimensional");
        WriteData(grpVar, &refForceDimensional, "refForceDimensional");

        WriteData(grpVar, &physicalTimeStep, "physicalTimeStep");

        H5Gclose(grpVar);
    }

    int nLocalParticle = GetNumOfLocalParticle();

    ostringstream ossFolder;
    ossFolder << "Particle_Group" << gridID;
    grpName = ossFolder.str();
    grploc = OpenGroup(actkey->filepos, grpName);
    ossFolder.clear();
    ossFolder.str("");

    ostringstream ossData;

    //! Info
    ossData << "Info";
    grpName = ossData.str();
    grpVar = OpenGroup(grploc, grpName);
    WriteData(grpVar, &nLocalParticle, "nParticle");
    int *prticleID = particlePointGroup->GetGlobalParticleIDOnLocalZone();
    int nDimTem = parameter->GetNumDimTemperature();

    WriteData(grpVar, prticleID, "ParticleID");
    WriteData(grpVar, &nDimTem, "nDimTem");

    DelPointer(prticleID);
    H5Gclose(grpVar);
    ossData.clear();
    ossData.str("");

    if (0 == nLocalParticle)
    {
        H5Gclose(grploc);
        return;
    }

    //! BasicParam
    ossData << "BasicParam";
    grpName = ossData.str();
    grpVar = OpenGroup(grploc, grpName); 
    RDouble *particleDensity = particlePointGroup->GetParticleDensity(refDensityDimensional);
    RDouble *particleDiameter = particlePointGroup->GetParticleDiameter(refReynoldsLengthDimensional);
    RDouble *particleReNum = particlePointGroup->GetParticleReNum();
    RDouble *particleSpecificHeatCapacity = particlePointGroup->GetParticleSpecificHeatCapacity(refVelocityDimensional, refTemperatureDimensional);
    WriteData(grpVar, particleDensity, "Density");
    WriteData(grpVar, particleDiameter, "Diameter");
    WriteData(grpVar, particleReNum, "ParticleReNumber");
    WriteData(grpVar, particleSpecificHeatCapacity, "particleSpecificHeatCapacity");

    DelPointer(particleDiameter);
    DelPointer(particleDensity);
    DelPointer(particleReNum);
    DelPointer(particleSpecificHeatCapacity);
    H5Gclose(grpVar);
    ossData.clear();
    ossData.str("");

    //! Heat transfer
    ossData << "HeatTransfer";
    grpName = ossData.str();
    grpVar = OpenGroup(grploc, grpName);
    RDouble *particleConvectiveHeatTransfer = particlePointGroup->GetParticleConvectiveHeatTransfer(refDensityDimensional, refVelocityDimensional);
    RDouble *particleRadiativeHeatTransfer = particlePointGroup->GetParticleRadiativeHeatTransfer(refDensityDimensional, refVelocityDimensional);
    RDouble *particleNuseltNum = particlePointGroup->GetParticleNuseltNum();
    RDouble *particleWork = particlePointGroup->GetParticleWork(refDensityDimensional, refVelocityDimensional);
    RDouble *particleDissipation = particlePointGroup->GetParticleDissipation(refDensityDimensional, refVelocityDimensional);

    WriteData(grpVar, particleConvectiveHeatTransfer, "ConvectiveHeatTransfer");
    WriteData(grpVar, particleRadiativeHeatTransfer, "RadiativeHeatTransfer");
    WriteData(grpVar, particleNuseltNum, "NuseltNumber");
    WriteData(grpVar, particleWork, "ParticleWork");
    WriteData(grpVar, particleDissipation, "ParticleDissipation");

    DelPointer(particleConvectiveHeatTransfer);
    DelPointer(particleRadiativeHeatTransfer);
    DelPointer(particleNuseltNum);
    DelPointer(particleWork);
    DelPointer(particleDissipation);

    H5Gclose(grpVar);
    ossData.clear();
    ossData.str("");

    //! Coordinate.
    ossData << "Coordinate";
    grpName = ossData.str();
    grpVar = OpenGroup(grploc, grpName);
    RDouble *particleCoordinateX = particlePointGroup->GetParticleCoordinateDim(refReynoldsLengthDimensional, 0);
    RDouble *particleCoordinateY = particlePointGroup->GetParticleCoordinateDim(refReynoldsLengthDimensional, 1);
    RDouble *particleCoordinateZ = particlePointGroup->GetParticleCoordinateDim(refReynoldsLengthDimensional, 2);
    WriteData(grpVar, particleCoordinateX, "CoordinateX");
    WriteData(grpVar, particleCoordinateY, "CoordinateY");
    WriteData(grpVar, particleCoordinateZ, "CoordinateZ");
    DelPointer(particleCoordinateX);
    DelPointer(particleCoordinateY);
    DelPointer(particleCoordinateZ);
    H5Gclose(grpVar);
    ossData.clear();
    ossData.str("");

    //! Velocity.
    ossData << "Velocity";
    grpName = ossData.str();
    grpVar = OpenGroup(grploc, grpName);
    RDouble *particleVelocityX = particlePointGroup->GetParticleVelocityDim(refVelocityDimensional, 0);
    RDouble *particleVelocityY = particlePointGroup->GetParticleVelocityDim(refVelocityDimensional, 1);
    RDouble *particleVelocityZ = particlePointGroup->GetParticleVelocityDim(refVelocityDimensional, 2);
    WriteData(grpVar, particleVelocityX, "VelocityX");
    WriteData(grpVar, particleVelocityY, "VelocityY");
    WriteData(grpVar, particleVelocityZ, "VelocityZ");
    DelPointer(particleVelocityX);
    DelPointer(particleVelocityY);
    DelPointer(particleVelocityZ);
    H5Gclose(grpVar);
    ossData.clear();
    ossData.str("");

    //! Relative Velocity.
    ossData << "RelativeVelocity";
    grpName = ossData.str();
    grpVar = OpenGroup(grploc, grpName);
    RDouble *particleMaNum = particlePointGroup->GetParticleMaNum();
    RDouble *particleRelativeVelocityX = particlePointGroup->GetParticleRelativeVelocityDim(refVelocityDimensional, 0);
    RDouble *particleRelativeVelocityY = particlePointGroup->GetParticleRelativeVelocityDim(refVelocityDimensional, 1);
    RDouble *particleRelativeVelocityZ = particlePointGroup->GetParticleRelativeVelocityDim(refVelocityDimensional, 2);
    WriteData(grpVar, particleMaNum, "ParticleMaNumber");
    WriteData(grpVar, particleRelativeVelocityX, "RelativeVelocityX");
    WriteData(grpVar, particleRelativeVelocityY, "RelativeVelocityY");
    WriteData(grpVar, particleRelativeVelocityZ, "RelativeVelocityZ");
    DelPointer(particleMaNum);
    DelPointer(particleRelativeVelocityX);
    DelPointer(particleRelativeVelocityY);
    DelPointer(particleRelativeVelocityZ);
    H5Gclose(grpVar);
    ossData.clear();
    ossData.str("");

    //! Stokes number.
    ossData << "StokesNumber";
    grpName = ossData.str();
    grpVar = OpenGroup(grploc, grpName);
    RDouble *particleTotalForceSt = particlePointGroup->GetParticleTotalForceSt();
    RDouble *particleStokesSt = particlePointGroup->GetParticleStokesSt();
    RDouble *particleModifiedSt = particlePointGroup->GetParticleModifiedSt();
    WriteData(grpVar, particleTotalForceSt, "TotalForceSt");
    WriteData(grpVar, particleStokesSt, "StokesSt");
    WriteData(grpVar, particleModifiedSt, "ModifiedSt");
    DelPointer(particleTotalForceSt);
    DelPointer(particleStokesSt);
    DelPointer(particleModifiedSt);
    H5Gclose(grpVar);
    ossData.clear();
    ossData.str("");

    //! AngularVelocity.
    ossData << "AngularVelocity";
    grpName = ossData.str();
    grpVar = OpenGroup(grploc, grpName);
    RDouble *particleAngularVelocityX = particlePointGroup->GetParticleAngularVelocityDim(refAngularVelocityDimensional, 0);
    RDouble *particleAngularVelocityY = particlePointGroup->GetParticleAngularVelocityDim(refAngularVelocityDimensional, 1);
    RDouble *particleAngularVelocityZ = particlePointGroup->GetParticleAngularVelocityDim(refAngularVelocityDimensional, 2);
    WriteData(grpVar, particleAngularVelocityX, "AngularVelocityX");
    WriteData(grpVar, particleAngularVelocityX, "AngularVelocityY");
    WriteData(grpVar, particleAngularVelocityX, "AngularVelocityZ");
    DelPointer(particleAngularVelocityX);
    DelPointer(particleAngularVelocityY);
    DelPointer(particleAngularVelocityZ);
    H5Gclose(grpVar);
    ossData.clear();
    ossData.str("");

    //! Acceleration.
    ossData << "Acceleration";
    grpName = ossData.str();
    grpVar = OpenGroup(grploc, grpName);
    RDouble *particleAccelerationX = particlePointGroup->GetParticleAccelerationDim(refAccelerationyDimensional, 0);
    RDouble *particleAccelerationY = particlePointGroup->GetParticleAccelerationDim(refAccelerationyDimensional, 1);
    RDouble *particleAccelerationZ = particlePointGroup->GetParticleAccelerationDim(refAccelerationyDimensional, 2);
    WriteData(grpVar, particleAccelerationX, "AccelerationX");
    WriteData(grpVar, particleAccelerationY, "AccelerationY");
    WriteData(grpVar, particleAccelerationZ, "AccelerationZ");
    DelPointer(particleAccelerationX);
    DelPointer(particleAccelerationY);
    DelPointer(particleAccelerationZ);

    H5Gclose(grpVar);
    ossData.clear();
    ossData.str("");

    //! Temperature
    ossData << "Temperature";
    grpName = ossData.str();
    grpVar = OpenGroup(grploc, grpName);
    RDouble *particleTemperature0 = particlePointGroup->GetParticleTemperatureDim(refTemperatureDimensional,0);
    WriteData(grpVar, particleTemperature0, "Temperature0");
    DelPointer(particleTemperature0);
    H5Gclose(grpVar);
    ossData.clear();
    ossData.str("");

    //! Force.
    using namespace INDEX_PARTICLE_FORCE;
    RDouble *particleForceX,  *particleForceY,  *particleForceZ;
    int forceType = 0;
    ossData << "Force";
    grpName = ossData.str();
    grpVar = OpenGroup(grploc, grpName);
    //! Gravity
    forceType = FG;
    particleForceX = particlePointGroup->GetParticleForceDim(refForceDimensional, 0, forceType);
    particleForceY = particlePointGroup->GetParticleForceDim(refForceDimensional, 1, forceType);
    particleForceZ = particlePointGroup->GetParticleForceDim(refForceDimensional, 2, forceType);
    WriteData(grpVar, particleForceX, "GravityX");
    WriteData(grpVar, particleForceY, "GravityY");
    WriteData(grpVar, particleForceZ, "GravityZ");
    DelPointer(particleForceX);
    DelPointer(particleForceY);
    DelPointer(particleForceZ);
    //! Drag
    forceType = FD;
    particleForceX = particlePointGroup->GetParticleForceDim(refForceDimensional, 0, forceType);
    particleForceY = particlePointGroup->GetParticleForceDim(refForceDimensional, 1, forceType);
    particleForceZ = particlePointGroup->GetParticleForceDim(refForceDimensional, 2, forceType);
    WriteData(grpVar, particleForceX, "DragX");
    WriteData(grpVar, particleForceY, "DragY");
    WriteData(grpVar, particleForceZ, "DragZ");
    DelPointer(particleForceX);
    DelPointer(particleForceY);
    DelPointer(particleForceZ);
    //! Saffman
    forceType = FSAFF;
    particleForceX = particlePointGroup->GetParticleForceDim(refForceDimensional, 0, forceType);
    particleForceY = particlePointGroup->GetParticleForceDim(refForceDimensional, 1, forceType);
    particleForceZ = particlePointGroup->GetParticleForceDim(refForceDimensional, 2, forceType);
    WriteData(grpVar, particleForceX, "SaffmanX");
    WriteData(grpVar, particleForceY, "SaffmanY");
    WriteData(grpVar, particleForceZ, "SaffmanZ");
    DelPointer(particleForceX);
    DelPointer(particleForceY);
    DelPointer(particleForceZ);
    //! Magnus
    forceType = FMAG;
    particleForceX = particlePointGroup->GetParticleForceDim(refForceDimensional, 0, forceType);
    particleForceY = particlePointGroup->GetParticleForceDim(refForceDimensional, 1, forceType);
    particleForceZ = particlePointGroup->GetParticleForceDim(refForceDimensional, 2, forceType);
    WriteData(grpVar, particleForceX, "MagnusX");
    WriteData(grpVar, particleForceY, "MagnusY");
    WriteData(grpVar, particleForceZ, "MagnusZ");
    DelPointer(particleForceX);
    DelPointer(particleForceY);
    DelPointer(particleForceZ);
    //! Addmass
    forceType = FAM;
    particleForceX = particlePointGroup->GetParticleForceDim(refForceDimensional, 0, forceType);
    particleForceY = particlePointGroup->GetParticleForceDim(refForceDimensional, 1, forceType);
    particleForceZ = particlePointGroup->GetParticleForceDim(refForceDimensional, 2, forceType);
    WriteData(grpVar, particleForceX, "AddedMassX");
    WriteData(grpVar, particleForceY, "AddedMassY");
    WriteData(grpVar, particleForceZ, "AddedMassZ");
    DelPointer(particleForceX);
    DelPointer(particleForceY);
    DelPointer(particleForceZ);
    //! FluidAcceleration
    forceType = FFA;
    particleForceX = particlePointGroup->GetParticleForceDim(refForceDimensional, 0, forceType);
    particleForceY = particlePointGroup->GetParticleForceDim(refForceDimensional, 1, forceType);
    particleForceZ = particlePointGroup->GetParticleForceDim(refForceDimensional, 2, forceType);
    WriteData(grpVar, particleForceX, "FluidAccelerationX");
    WriteData(grpVar, particleForceY, "FluidAccelerationY");
    WriteData(grpVar, particleForceZ, "FluidAccelerationZ");
    DelPointer(particleForceX);
    DelPointer(particleForceY);
    DelPointer(particleForceZ);
    //! Brownian
    forceType = FB;
    particleForceX = particlePointGroup->GetParticleForceDim(refForceDimensional, 0, forceType);
    particleForceY = particlePointGroup->GetParticleForceDim(refForceDimensional, 1, forceType);
    particleForceZ = particlePointGroup->GetParticleForceDim(refForceDimensional, 2, forceType);
    WriteData(grpVar, particleForceX, "BrownianX");
    WriteData(grpVar, particleForceY, "BrownianY");
    WriteData(grpVar, particleForceZ, "BrownianZ");
    DelPointer(particleForceX);
    DelPointer(particleForceY);
    DelPointer(particleForceZ);
    //! Thermophoretic
    forceType = FT;
    particleForceX = particlePointGroup->GetParticleForceDim(refForceDimensional, 0, forceType);
    particleForceY = particlePointGroup->GetParticleForceDim(refForceDimensional, 1, forceType);
    particleForceZ = particlePointGroup->GetParticleForceDim(refForceDimensional, 2, forceType);
    WriteData(grpVar, particleForceX, "ThermophoreticX");
    WriteData(grpVar, particleForceY, "ThermophoreticY");
    WriteData(grpVar, particleForceZ, "ThermophoreticZ");
    DelPointer(particleForceX);
    DelPointer(particleForceY);
    DelPointer(particleForceZ);

    H5Gclose(grpVar);
    ossData.clear();
    ossData.str("");

    H5Gclose(grploc);

}

int ParticlePointSolverParallelDomain::GetNumOfLocalParticle()
{
    Param_ParticlePointSolverParallelDomain *parameter = GetControlParameters();
    return parameter->GetNumParticleLocal();
}

int ParticlePointSolverParallelDomain::GetNumOfTotalParticle()
{
    Param_ParticlePointSolverParallelDomain *parameter = GetControlParameters();
    return parameter->GetNumParticleTotal();
}

ParticlePointGroup *ParticlePointSolverParallelDomain::GetParticlePointGroup(string varibaleName)
{
    int level = 0;
    Grid *grid = GetGrid(level);
    ParticlePointGroup *particlePointGroup = reinterpret_cast<ParticlePointGroup*> (grid->GetDataPtr(varibaleName));
    return particlePointGroup;
}

LagrangianVariable<RDouble> *ParticlePointSolverParallelDomain::GetLagVariableDouble(string varibaleName)
{
    int level = 0;
    Grid *grid = GetGrid(level);
    LDouble *lagrangianVariable = reinterpret_cast<LDouble*> (grid->GetDataPtr(varibaleName));
    return lagrangianVariable;
}

void ParticlePointSolverParallelDomain::CalcFlowOnLocationOfParticle(ActionKey *actkey)
{
    Grid *grid = GetGrid();
    int gridType = grid->Type();

    Param_ParticlePointSolverParallelDomain *parameter = GetControlParameters();

    //! Init for particlePointGroup.
    ParticlePointGroup *particlePointGroup = GetParticlePointGroup("particlePointGroup");
    for (map<int, OnePointVariable*>::iterator iter = particlePointGroup->GetIterBegin(); iter != particlePointGroup->GetIterEnd(); ++iter)
    {
        //! The global particleID.
        int particleID = iter->first;

        //! OnePointVariable
        OnePointVariable *onePointVariable = iter->second;
        ParticleInterpolation(actkey, grid, onePointVariable);
    }

    //! Init for particleBCInit.
    if (parameter->ifInitParticleBCInfo())
    {
        ParticlePointGroup *particleBCInit = GetParticlePointGroup("particleBCInit");
        for (map<int, OnePointVariable*>::iterator iter = particleBCInit->GetIterBegin(); iter != particleBCInit->GetIterEnd(); ++iter)
        {
            //! The global particleID.
            int particleID = iter->first;

            //! OnePointVariable
            OnePointVariable *onePointVariable = iter->second;
            ParticleInterpolation(actkey, grid, onePointVariable);
        }
    }
}

void ParticlePointSolverParallelDomain::InitLagrangianVariable()
{
    int level = 0;
    Grid *grid = GetGrid(level);
    Param_ParticlePointSolverParallelDomain *param_ParticlePointSolverParallelDomain = GetControlParameters();
    ParticlePointGroup *particlePointGroup = GetParticlePointGroup("particlePointGroup");

    //! The number of particle on current zone.
    int nLocalParticle = particlePointGroup->GetNumOfLocalParticle();

    //! The list of particleID on current zone.
    int *particleIDList = particlePointGroup->GetGlobalParticleIDOnLocalZone();;
}

void ParticlePointSolverParallelDomain::SolverMultiphaseFlow(FieldProxy *rhsProxy)
{
    int innerIterStep = PHSPACE::GlobalDataBase::GetIntParaFromDB("innstep");
    int outnstep = PHSPACE::GlobalDataBase::GetIntParaFromDB("outnstep");

    if (innerIterStep != 1)
    {
        return;
    }

    if (0 == GetNumOfLocalParticle())
    {
        return;
    }

    int level = 0;
    Grid *grid = GetGrid(level);

    //! Get particle point group.
    ParticlePointGroup *particlePointGroup = GetParticlePointGroup("particlePointGroup");
    Param_ParticlePointSolverParallelDomain *parameter = GetControlParameters();

    //! Add Particle.
    if (parameter->ifAddParticleInitBC())
    {
        if (outnstep % parameter->GetInitParticleBCPeriodicStep() == 0)
        {
            this->AddParticleAsInitBC(grid);
        }
    }

    //! Calculate flowOnParticle variable in Eular grid.
    ActionKey *actkeyCommFlowOnParticle = new ActionKey();
    FillActionKey(actkeyCommFlowOnParticle, FLOW_ON_PARTICLE_MPI, 0);

    //! Init the flow on flowOnParticleInEular.
    actkeyCommFlowOnParticle->subact = 1;
    CalculateVariableForParticleOnEular(actkeyCommFlowOnParticle);

    //! MPI flowOnParticleInEular
    actkeyCommFlowOnParticle->kind = 0;
    CommunicationFlowOnParticleInEular(actkeyCommFlowOnParticle, grid);

    //! Calc the flow on flowOnParticleInEularGrad and grad second.
    actkeyCommFlowOnParticle->subact = 2;
    CalculateVariableForParticleOnEular(actkeyCommFlowOnParticle);

    //! MPI flowOnParticleInEularGrad and grad2.
    actkeyCommFlowOnParticle->kind = 2;
    CommunicationFlowOnParticleInEular(actkeyCommFlowOnParticle, grid);

    delete actkeyCommFlowOnParticle;    actkeyCommFlowOnParticle = nullptr;

    //! Init source.
    InitSource(grid);

    //! Runge Kutta Coef.
    RungeKuttaIndex *rungeKutta = new RungeKuttaIndex(GetControlParameters()->GetNumRK());

    //! Add ref param from solver to particleParam.
    Data_Param *particleParam = new Data_Param;
    parameter->GetRefDataParam(particleParam);
    parameter->GetForceTypeDataParam(particleParam);

    ActionKey *actkeyParticleInterpolation = new ActionKey();
    FillActionKey(actkeyParticleInterpolation, FLOW_ON_PARTICLE_MPI, 0);
    actkeyParticleInterpolation->subact = 1;
    IndexCell *indexCell = particlePointGroup->GetIndexCell();
    //! Loop for each particle.
    map<int, OnePointVariable*>::iterator iter = particlePointGroup->GetIterBegin();
    while (iter != particlePointGroup->GetIterEnd())
    {
        //! The global particleID.
        int particleID = iter->first;

        //! OnePointVariable
        OnePointVariable *onePointVariable = iter->second;

        //! Interpolation for partcle loaction.
        int cellID = onePointVariable->GetCellID();
        int *cellIndexID = this->GetCellIndex(grid, cellID);

        //! Stores the start values xi,yi of RK iteration.
        StoresParticleStartValue(onePointVariable);

        int ifOneRK = 0;

        rungeKutta->ResetIndex();

        //! Runge Kutta Loop.
        for (int iRK = 0; iRK < rungeKutta->GetNumRK(); ++iRK)
        {
            if (0 == ifOneRK)
            {
                //! Particle interpolation.
                ParticleInterpolation(actkeyParticleInterpolation, grid, onePointVariable);

                //! Calculate the force on the particle from the fluid.
                ParticleForceModel::CalcParticleForce(onePointVariable, particleParam);

                //! Update the RK iteration to advance the particle governing equation.
                UpdateRK4(grid, onePointVariable, rungeKutta, particleParam);

                //! Particle tracking in Lagrange to update the cell ID.
                ParticleTrack(grid, onePointVariable);

                //! Determine whether the particle touches the boundary 
                //! and add particle boundary information.
                //! Here we need to determine whether the particle is touching the boundary. 
                //! If the particle is near the critical region of the boundary, 
                //! we need to record it and add it to the IndexCell class.
                indexCell->SetParticleInCellInfo(grid, onePointVariable);

                ifOneRK = onePointVariable->GetIfOneRK();

                //! Push the RK coefficient.
                rungeKutta->AddIndex();
            }
            else
            {
                break;
            }
        }

        if (ifOneRK == 1)
        {
            UpdateRK1(grid, onePointVariable, rungeKutta, particleParam);
            //! Particle tracking in Lagrange to update the cell ID.
            ParticleTrack(grid, onePointVariable);
            indexCell->SetParticleInCellInfo(grid, onePointVariable);
        }
        
        bool ifRemoveParticle = ParticleBoundary(grid, onePointVariable,indexCell);

        if (ifRemoveParticle)
        {
            delete iter->second;
            iter->second = nullptr;
            (particlePointGroup->GetMap())->erase(iter++);
        }
        else
        {
            //! Calculate Two-way coupling.
            GetSourceFromParticle2Flow(grid, onePointVariable, particleParam);

            iter++;
        }

    }

    delete rungeKutta;    rungeKutta = nullptr;

    delete actkeyParticleInterpolation;    actkeyParticleInterpolation = nullptr;

    delete particleParam;    particleParam = nullptr;

    //! Update particle number.
    parameter->SetNumParticleLocal(particlePointGroup->GetNumOfLocalParticle());

    //! Calculate particle CFL on current zone..
    CalcCFLNumber();

}

void ParticlePointSolverParallelDomain::CommunicateMultiphaseFlow(FieldProxy *rhsProxy)
{
    int innerIterStep = PHSPACE::GlobalDataBase::GetIntParaFromDB("innstep");
    if (innerIterStep != 1)
    {
        return;
    }

    int level = 0;
    Grid *grid = GetGrid(level);

    UpdateNumOfParticle();

    //! This part is communicate particle information between zones and processors, we acheived this
    //! once and we are trying a better way.
    //! CommunicationParticle(grid, rhsProxy);
}

void ParticlePointSolverParallelDomain::ParticleTrack(Grid *grid, OnePointVariable *onePointVariable)
{
    int gridType = grid->Type();

    if (gridType == PHSPACE::STRUCTGRID)
    {
        StructGrid *structgrid = StructGridCast(grid);
        ParticleTrack(structgrid, onePointVariable);
    }
    else if (gridType == PHSPACE::UNSTRUCTGRID)
    {
        UnstructGrid *unstructgrid = UnstructGridCast(grid);
        ParticleTrack(unstructgrid, onePointVariable);
    }
    else if (gridType == PHSPACE::MIXGRID)
    {
        //! Mixing
        TK_Exit::UnexpectedVarValue(" gridType = ", PHSPACE::MIXGRID);
    }
}

void ParticlePointSolverParallelDomain::ParticleTrack(StructGrid *grid, OnePointVariable *onePointVariable)
{
    //! Grid bc
    StructBCSet *structBCSet = grid->GetStructBCSet();
    int nBCRegion = structBCSet->GetnBCRegion();
    //! Loop for each bc on current zone.
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++iBCRegion)
    {
        StructBC *structBC = structBCSet->GetBCRegion(iBCRegion);
        int bcType = structBC->GetBCType();
        int particleBCType = structBC->GetParticleBCType();
    }
    //! Index start and end of surface in .From 1 to nDim.
    int nDimGlobal = GetDim();
    int iSurfaceStart = 1;
    int iSurfaceEnd = nDimGlobal;
    //! Get face normal.
    RDouble4D &xfn = *(grid->GetFaceNormalX());
    RDouble4D &yfn = *(grid->GetFaceNormalY());
    RDouble4D &zfn = *(grid->GetFaceNormalZ());
    //! Cell center coordinate.
    RDouble3D &xcc = *(grid->GetCellCenterX());
    RDouble3D &ycc = *(grid->GetCellCenterY());
    RDouble3D &zcc = *(grid->GetCellCenterZ());
    //! Cell node coordinate.
    RDouble3D &xNode = *(grid->GetStructX());
    RDouble3D &yNode = *(grid->GetStructY());
    RDouble3D &zNode = *(grid->GetStructZ());
    //! Index start and end of cell in current zone.
    int ist, ied;
    int jst, jed;
    int kst, ked;
    int nLayers = GetNumberOfGhostCellLayers();
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    //! Particle Information 
    int cellID = onePointVariable->GetCellID();
    RDouble particleDiameter = onePointVariable->GetOnePointDiameter();
    SPDouble &particleVelocity = *(onePointVariable->GetOnePointVelocity());

    //! New particle coordinate.
    SPDouble &newCoordinate = *(onePointVariable->GetOnePointCoordinate());
    //! Old particle coordinate.
    int nDimVar = 3;
    SPDouble &oldCoordinate = *(onePointVariable->GetOnePointCoordinateOld());

    //! Get index IJK of particle cell position in current zone.
    int nI, nJ, nK;
    int iP, jP, kP;
    grid->GetND(nI, nJ, nK);
    nI -= 1;
    nJ -= 1;
    nK -= 1;
    if (GetDim() == TWO_D)
    {
        nK = 1;
    }
    GetIndexCellIJKStruct(nI, nJ, nK, iP, jP, kP, nLayers, cellID);

    //! Particle Track by
    //! "Particle tracking in unstructured, 
    //! arbitrary polyhedral meshes for use in CFD and molecular dynamics"
    //! Macpherson,2009
    SPDouble lambda_a(nDimGlobal * 2);
    SPDouble lambda_a1(nDimVar);
    SPDouble lambda_a2(nDimVar);

    SPDouble lambda_c(nDimGlobal * 2);
    SPDouble lambda_c1(nDimVar);
    SPDouble lambda_c2(nDimVar);

    //! The local index of surface in one cell.
    //! The index start from 0 to nDim-1;
    int iSurLocal, jSurLocal, kSurLocal;
    //! The local index for parallel surface of local cell.
    int iSurLocalParallel, jSurLocalParallel, kSurLocalParallel;
    //! The cell id when while loop (After loop, it is the new cell).
    int idCellLoop = -1;
    //! If in loop.
    int ifWhile = 1;
    int nWhile = 0;

    //! Loop to get the new index iP,jP,kP,cellID.
    while ( 1 == ifWhile)
    {
        //! Store the lowest lambda_a.
        RDouble lambda_value = 1.0;
        //! The index which surface is the lowest lambda_a.
        int indexSurface = 0;

        //! How many surfaces are used to compare the lowest lambda_a.
        int countSurface = 0;

        //! Loop for surface.
        //! This loop will find lowest lambda_a when 0<= lambda_c <=1.
        //! Otherwise, cell new is cell old.
        for (int iSurface = iSurfaceStart; iSurface <= iSurfaceEnd; ++iSurface)
        {
            //! The center coordinate of surface.
            RDouble xfc1, yfc1, zfc1;
            //! Get the center coodinate of surface.
            grid->FaceCoor(iP, jP, kP, iSurface, xfc1, yfc1, zfc1);
            //! Get the normal vector of surface
            RDouble xfn1 = xfn(iP, jP, kP, iSurface);
            RDouble yfn1 = yfn(iP, jP, kP, iSurface);
            RDouble zfn1 = zfn(iP, jP, kP, iSurface);

            //! a particle is located at position A and is required to move to B.
            //! B is newCoordinate and A is oldCoordinate.
            //! lambda_a = \frac{(C_f - a)\cdot S}{(b-a)\codt S} in Macpherson,2009.
            //! where S is face normal vector (outward) of surface of current cell,
            //! Cf is a coordinate vector at the face centre of surface of cell,
            //! a and b are the coordinate of point A and B.
            //! Here,  (C_f - a)\cdot S is lambda_a1[0] + lambda_a1[1] + lambda_a1[2].
            //! which is a scalar.
            //! (b-a)\codt S is lambda_a2[0] + lambda_a2[1] + lambda_a2[2].
            //! Note here, in PHengLEI, xfn and xfc is inwards only at (i,j,k).
            lambda_a1[0] = - (xfc1 - oldCoordinate[0]) * xfn1;
            lambda_a1[1] = - (yfc1 - oldCoordinate[1]) * yfn1;
            lambda_a1[2] = - (zfc1 - oldCoordinate[2]) * zfn1;

            lambda_a2[0] = - (newCoordinate[0] - oldCoordinate[0]) * xfn1;
            lambda_a2[1] = - (newCoordinate[1] - oldCoordinate[1]) * yfn1;
            lambda_a2[2] = - (newCoordinate[2] - oldCoordinate[2]) * zfn1;

            lambda_a[iSurface - iSurfaceStart] = (lambda_a1[0] + lambda_a1[1] + lambda_a1[2]) / (lambda_a2[0] + lambda_a2[1] + lambda_a2[2]);

            //! lambda_c = \frac{(C_f - C_c)\cdot S}{(b-C_c)\cdot S}
            //! C_c is a coordinate vector at cell center.
            //! (C_f - C_c)\cdot S is lambda_c1[0] + lambda_c1[1] + lambda_c1[2],
            //! and (b-C_c)\cdot S is lambda_c2[i] + lambda_c2[1] + lambda_c2[2].
            RDouble xcc1 = xcc(iP, jP, kP);
            RDouble ycc1 = ycc(iP, jP, kP);
            RDouble zcc1 = zcc(iP, jP, kP);

            lambda_c1[0] = - (xfc1 - xcc1) * xfn1;
            lambda_c1[1] = - (yfc1 - ycc1) * yfn1;
            lambda_c1[2] = - (zfc1 - zcc1) * zfn1;

            lambda_c2[0] = - (newCoordinate[0] - xcc1) * xfn1;
            lambda_c2[1] = - (newCoordinate[1] - ycc1) * yfn1;
            lambda_c2[2] = - (newCoordinate[2] - zcc1) * zfn1;

            lambda_c[iSurface - iSurfaceStart] = (lambda_c1[0] + lambda_c1[1] + lambda_c1[2]) / (lambda_c2[0] + lambda_c2[1] + lambda_c2[2]);

            //! Using lambda_c to find which face the particle will hit result in 0 <= lambda_c <= 1.0.
            //! The lowest value of lambda_a determines which face was actually hit.
            //! If lambda_c < 0 or lambda_c > 1, B must be inside the old cell.
            if (lambda_c[iSurface - iSurfaceStart] >= 0.0 && lambda_c[iSurface - iSurfaceStart] <= 1.0)
            {
                if (lambda_a[iSurface - iSurfaceStart] < lambda_value)
                {
                    lambda_value = lambda_a[iSurface - iSurfaceStart];
                    indexSurface = iSurface;
                    countSurface += 1;
                }
            }
        
            //! Get the local index of parallel surface.
            //! Get the local index of surface in current one cell.
            GetNsurfIndex(iSurface, iSurLocal, jSurLocal, kSurLocal);
            iSurLocalParallel = iP + iSurLocal;
            jSurLocalParallel = jP + jSurLocal;
            kSurLocalParallel = kP + kSurLocal;

            //! Get the normal vector of parallel surface.
            RDouble xfn2 = xfn(iSurLocalParallel, jSurLocalParallel, kSurLocalParallel, iSurface);
            RDouble yfn2 = yfn(iSurLocalParallel, jSurLocalParallel, kSurLocalParallel, iSurface);
            RDouble zfn2 = zfn(iSurLocalParallel, jSurLocalParallel, kSurLocalParallel, iSurface);

            //! The center coordinate of parallel surface.
            RDouble xfc2, yfc2, zfc2;
            grid->FaceCoor(iSurLocalParallel, jSurLocalParallel, kSurLocalParallel, iSurface, xfc2, yfc2, zfc2);

            //! a particle is located at position A and is required to move to B.
            //! B is newCoordinate and A is oldCoordinate.
            //! lambda_a = \frac{(C_f - a)\cdot S}{(b-a)\codt S} in Macpherson,2009.
            //! where S is face normal vector (outward) of surface of current cell,
            //! Cf is a coordinate vector at the face centre of surface of cell,
            //! a and b are the coordinate of point A and B.
            //! Here,  (C_f - a)\cdot S is lambda_a1[0] + lambda_a1[1] + lambda_a1[2].
            //! (b-a)\codt S is lambda_a2[0] + lambda_a2[1] + lambda_a2[2].
            //! Note here, in PHengLEI, xfn and xfc is inwards only at (i,j,k).
            lambda_a1[0] = - (xfc2 - oldCoordinate[0]) * xfn2;
            lambda_a1[1] = - (yfc2 - oldCoordinate[1]) * yfn2;
            lambda_a1[2] = - (zfc2 - oldCoordinate[2]) * zfn2;

            lambda_a2[0] = - (newCoordinate[0] - oldCoordinate[0]) * xfn2;
            lambda_a2[1] = - (newCoordinate[1] - oldCoordinate[1]) * yfn2;
            lambda_a2[2] = - (newCoordinate[2] - oldCoordinate[2]) * yfn2;

            lambda_a[iSurface - iSurfaceStart + nDimGlobal] = (lambda_a1[0] + lambda_a1[1] + lambda_a1[2]) / (lambda_a2[0] + lambda_a2[1] + lambda_a2[2]);

            //! Note here xcc is the coordinate of Cell A,
            //! so xcc2 is still be xcc1 ,it only loop for face but not cell.
            //! lambda_c = \frac{(C_f - C_c)\cdot S}{(b-C_c)\cdot S}
            //! C_c is a coordinate vector at cell center.
            //! (C_f - C_c)\cdot S is lambda_c1[0] + lambda_c1[1] + lambda_c1[2],
            //! and (b-C_c)\cdot S is lambda_c2[i] + lambda_c2[1] + lambda_c2[2].
            lambda_c1[0] = - (xfc2 - xcc1) * xfn2;
            lambda_c1[1] = - (yfc2 - ycc1) * yfn2;
            lambda_c1[2] = - (zfc2 - zcc1) * zfn2;

            lambda_c2[0] = - (newCoordinate[0] - xcc1) * xfn2;
            lambda_c2[1] = - (newCoordinate[1] - ycc1) * yfn2;
            lambda_c2[2] = - (newCoordinate[2] - zcc1) * zfn2;

            lambda_c[iSurface - iSurfaceStart + nDimGlobal] = (lambda_c1[0] + lambda_c1[1] + lambda_c1[2]) / (lambda_c2[0] + lambda_c2[1] + lambda_c2[2]);

            if (lambda_c[iSurface - iSurfaceStart + nDimGlobal] >= 0.0 && lambda_c[iSurface - iSurfaceStart + nDimGlobal] <= 1.0)
            {
                if (lambda_a[iSurface - iSurfaceStart + nDimGlobal] < lambda_value)
                {
                    lambda_value = lambda_a[iSurface - iSurfaceStart + nDimGlobal];
                    indexSurface = iSurface + nDimGlobal;
                    countSurface += 1;
                }
            }
        }
    
        //! Move the point.
        if (0 == countSurface)
        {
            //! cell new is cell old. The cell ID stays the same.
            ifWhile = 0;
        }
        else
        {
            //! move the particle according to the hit surface (hit point P).
            //! P = A + \lambda_m \times(B-A).
            //! \lambda_m = min(1,max(0,\lambda_a)).
            RDouble lambda_m = min(1.0, max(0.0, lambda_value));
            for (int iDim = 0; iDim < nDimGlobal; ++iDim)
            {
                oldCoordinate[iDim] = oldCoordinate[iDim] + lambda_m * (newCoordinate[iDim] - oldCoordinate[iDim]);
                if (std::isnan(newCoordinate[iDim]))
                {
                    cout << "error" << endl;
                }
            }
        }
        
        //! Get the new cell index.
        int iM, jM, kM;
        GetIndexIJKStructOnTrackSurface(indexSurface, iM, jM, kM);
        //! The new index iP,jP,kP.
        iP = iP + iM;
        jP = jP + jM;
        kP = kP + kM;
        GetIDCellOfParticleByIJKSturct(nI, nJ, nK, iP, jP, kP, nLayers, idCellLoop);
    }

    onePointVariable->SetCellID(idCellLoop);
}

void ParticlePointSolverParallelDomain::GetIndexIJKStructOnTrackSurface(int &indexSurface, int &iM, int &jM, int &kM)
{
    int nDimGlobal = GetDim();

    iM = 0;
    jM = 0;
    kM = 0;

    if (nDimGlobal == TWO_D)
    {
        switch (indexSurface)
        {
        case 1:
            iM = -1;
            break;
        case 2:
            jM = -1;
            break;
        case 3:
            iM = +1;
            break;
        case 4:
            jM = +1;
            break;
        default:
            break;
        }
    }
    else
    {
        switch (indexSurface)
        {
        case 1:
            iM = -1;
            break;
        case 2:
            jM = -1;
            break;
        case 3:
            kM = -1;
            break;
        case 4:
            iM = +1;
            break;
        case 5:
            jM = +1;
            break;
        case 6:
            kM = +1;
        default:
            break;
        }
    }
}

void ParticlePointSolverParallelDomain::ParticleTrack(UnstructGrid *grid, OnePointVariable *onePointVariable)
{

}

int ParticlePointSolverParallelDomain::CreatDataContanier(int &nParticle, DataContainer *cdata)
{
    int isize = 0;

    for (int iParticle = 0; iParticle < nParticle; ++iParticle)
    {
        isize = 0;

        int particleID = -1;
        cdata->Write(&particleID, sizeof(int));
        isize += sizeof(int);

        int iGlobalzone = -1;
        cdata->Write(&iGlobalzone, sizeof(int));
        isize += sizeof(int);

        int iP, jP, kP;
        iP = 0;
        jP = 0;
        kP = 0;

        cdata->Write(&iP, sizeof(int));
        cdata->Write(&jP, sizeof(int));
        cdata->Write(&kP, sizeof(int));
        isize += sizeof(int);
        isize += sizeof(int);
        isize += sizeof(int);

        int nDimVar = 12;
        RDouble *data = new RDouble[nDimVar];

        for (int idata = 0; idata < nDimVar; ++idata)
        {
            data[idata] = 1.0;
            cdata->Write(&data[idata], sizeof(RDouble));
            isize += sizeof(RDouble);
        }

        delete [] data;    data = nullptr;

    }

    return isize;
}

void ParticlePointSolverParallelDomain::StoresParticleStartValue(OnePointVariable *onePointVariable)
{
    //! Get the velocity of particle.
    RDouble particleDiameter = onePointVariable->GetOnePointDiameter();
    RDouble particleDensity = onePointVariable->GetOnePointDensity();
    SPDouble &particleVelocity = *(onePointVariable->GetOnePointVelocity());
    SPDouble &particleCoordinate = *(onePointVariable->GetOnePointCoordinate());
    SPDouble &particleAngularVelocity = *(onePointVariable->GetOnePointAngularVelocity());
    SPDouble &particleTemperature = *(onePointVariable->GetOnePointTemperature());
    SPDouble &initRKValue = *(onePointVariable->GetOnePointInitRKValue());
    SPDouble &dRKValue = *(onePointVariable->GetOnePointDRKValue());
    SPDouble &sumDRK = *(onePointVariable->GetOnePointSumDRK());

    RDouble initValue = 0.0;
    int nRKVar = 10;

    initRKValue.CopyData(initValue, nRKVar);
    dRKValue.CopyData(initValue, nRKVar);
    sumDRK.CopyData(initValue, nRKVar);

    using namespace INDEX_RK;
    int nDimVar = 3;
    for (int iDim = 0; iDim < nDimVar; ++iDim)
    {
        initRKValue[RK_COORD + iDim] = particleCoordinate[iDim];
        initRKValue[RK_VEL + iDim] = particleVelocity[iDim];
        initRKValue[RK_ANGULAR_VEL + iDim] = particleAngularVelocity[iDim];
    }
    initRKValue[RK_TEMP] = particleTemperature[0];

    onePointVariable->SetIfOneRK(0);
    onePointVariable->SetOldCellID(onePointVariable->GetCellID());
}

void ParticlePointSolverParallelDomain::UpdateRK4(Grid *grid, OnePointVariable *onePointVariable, RungeKuttaIndex *indexRK, Data_Param *particleParam)
{
    int nDimVar;

    //! Get the variable of particle.
    RDouble particleDiameter = onePointVariable->GetOnePointDiameter();
    RDouble particleDensity = onePointVariable->GetOnePointDensity();
    SPDouble &particleVelocity = *(onePointVariable->GetOnePointVelocity());
    SPDouble &particleCoordinate = *(onePointVariable->GetOnePointCoordinate());
    SPDouble &particleAngularVelocity = *(onePointVariable->GetOnePointAngularVelocity());
    SPDouble &particleTemperature = *(onePointVariable->GetOnePointTemperature());
    SPDouble &particleForce = *(onePointVariable->GetOnePointForce());

    RDouble particleMass = (4.0 / 3.0) * PI * pow(particleDiameter * 0.5, 3.0) * particleDensity;
    RDouble physicalTimeStep = GlobalDataBase::GetDoubleParaFromDB("physicalTimeStep");

    //! RK variable.
    using namespace INDEX_RK;
    SPDouble &initRKValue = *(onePointVariable->GetOnePointInitRKValue());
    SPDouble &dRKValue = *(onePointVariable->GetOnePointDRKValue());
    SPDouble &sumDRK = *(onePointVariable->GetOnePointSumDRK());

    //! Reset rk variable without initRK.
    nDimVar = 10;
    dRKValue.CopyData(0.0, nDimVar);

    //! Coordinate.
    dRKValue[RK_COORD + 0] = particleVelocity[0];
    dRKValue[RK_COORD + 1] = particleVelocity[1];
    dRKValue[RK_COORD + 2] = particleVelocity[2];

    //! velocity.
    using namespace INDEX_PARTICLE_FORCE;
    nDimVar = 3;
    for (int iForce = 0; iForce < nForceDir/ nDimVar; ++iForce)
    {
        dRKValue[RK_VEL + 0] += particleForce[iForce * nDimVar + 0] / (particleMass+SMALL);
        dRKValue[RK_VEL + 1] += particleForce[iForce * nDimVar + 1] / (particleMass + SMALL);
        dRKValue[RK_VEL + 2] += particleForce[iForce * nDimVar + 2] / (particleMass + SMALL);
    }

    dRKValue[RK_TEMP] = particleTemperature[0];
    nDimVar = 10;
    for (int iDim = 0; iDim < nDimVar; ++iDim)
    {
        sumDRK[iDim] += (indexRK->GetKindex()) * dRKValue[iDim];
    }

    ostringstream ossLog;
    ossLog << " --- init rk 4 ---" << "\n";
    ossLog << " outnstep " << GlobalDataBase::GetIntParaFromDB("outnstep") << "\n";
    ossLog << " particle ID : " << onePointVariable->GetParticleID() << "\n";
    ossLog << " iRK : " << indexRK->GetIndex() << "\n";
    ossLog << " initRKValue coord : ";
    ossLog << initRKValue[RK_COORD + 0] << " " << initRKValue[RK_COORD + 1] << " " << initRKValue[RK_COORD + 2] << "\n";
    ossLog << " initRKValue vel : ";
    ossLog << initRKValue[RK_VEL + 0] << " " << initRKValue[RK_VEL + 1] << " " << initRKValue[RK_VEL + 2] << "\n";
    ossLog << " dRKValue coord :";
    ossLog << dRKValue[RK_COORD + 0] << " " << dRKValue[RK_COORD + 1] << " " << dRKValue[RK_COORD + 2] << "\n";
    ossLog << " dRKValue vel : ";
    ossLog << dRKValue[RK_VEL + 0] << " " << dRKValue[RK_VEL + 1] << " " << dRKValue[RK_VEL + 2] << "\n";
    ossLog << "\n";

    nDimVar = 3;
    if (0 == onePointVariable->GetIfOneRK())
    {
        //! RK4.
        //! The equation y = f(x) and f(x_0) = y_0
        //! k1 = f(x_0)
        //! k2 = f(x_0 + 0.5*h)
        //! k3 = f(x_0 + 0.5*h)
        //! k4 = f(x_0 + 1.0*h)
        //! f(x_1) = y_0 + (h/6) *(k1 + 2*k2 + 2*k3+k4)
        //! dRKValue is k.
        //! initRKValue is y_0.
        //! sumDRK is (k1 + 2*k2 + 2*k3+k4)

        if (indexRK->GetIndex() == indexRK->GetNumRK() - 1)
        {
            onePointVariable->SetCoordinateOld();
            for (int iDim = 0; iDim < nDimVar; ++iDim)
            {
                particleCoordinate[iDim] = initRKValue[RK_COORD + iDim] + physicalTimeStep * (indexRK->GetFindex()) * sumDRK[RK_COORD + iDim];

                particleVelocity[iDim] = initRKValue[RK_VEL + iDim] + physicalTimeStep * (indexRK->GetFindex()) * sumDRK[RK_VEL + iDim];
            }
            particleTemperature[0] = initRKValue[RK_TEMP] + physicalTimeStep * (indexRK->GetFindex()) * dRKValue[RK_TEMP];
        }
        else
        {
            onePointVariable->SetCoordinateOld();
            for (int iDim = 0; iDim < nDimVar; ++iDim)
            {
                particleCoordinate[iDim] = initRKValue[RK_COORD + iDim] + physicalTimeStep * (indexRK->GetFindex()) * dRKValue[RK_COORD + iDim];
                particleVelocity[iDim] = initRKValue[RK_VEL + iDim] + physicalTimeStep * (indexRK->GetFindex()) * dRKValue[RK_VEL + iDim];
            }
            particleTemperature[0] = initRKValue[RK_TEMP] + physicalTimeStep * (indexRK->GetFindex()) * dRKValue[RK_TEMP];
        }

    }
    else
    {
        ostringstream ossError;
        ossError << " error on rk4 " << onePointVariable->GetIfOneRK() << endl;
        TK_Exit::ExceptionExit(ossError);
    }

    ossLog << " dRKValue coord :" << dRKValue[RK_COORD + 0] << " " << dRKValue[RK_COORD + 1] << " " << dRKValue[RK_COORD + 2] << "\n";
    ossLog << " dRKValue vel : " << dRKValue[RK_VEL + 0] << " " << dRKValue[RK_VEL + 1] << " " << dRKValue[RK_VEL + 2] << "\n";
    ossLog << " particleCoordinate : ";
    ossLog << particleCoordinate[0] << " " << particleCoordinate[1] << " " << particleCoordinate[2] << "\n";
    ossLog << " particleVelocity : ";
    ossLog << particleVelocity[0] << " " << particleVelocity[1] << " " << particleVelocity[2] << "\n";
    ossLog << " cellID : " << onePointVariable->GetCellID() << "\n";
    ossLog << " (indexRK->GetFindex()) " << (indexRK->GetFindex()) << "\n";
    ossLog << " (indexRK->GetKindex()) " << (indexRK->GetKindex()) << "\n";
    ossLog << " physicalTimeStep : " << physicalTimeStep << "\n";
    ossLog << " sumDRK coord : ";
    ossLog << sumDRK[RK_COORD + 0] << " " << sumDRK[RK_COORD + 1] << " " << sumDRK[RK_COORD + 2] << "\n";
    ossLog << " sumDRK vel : ";
    ossLog << sumDRK[RK_VEL + 0] << " " << sumDRK[RK_VEL + 1] << " " << sumDRK[RK_VEL + 2] << "\n";
    ossLog << "\n" << endl;

    bool writeEachProcessor = false;

}

    //! Particle interpolation.
void ParticlePointSolverParallelDomain::UpdateRK1(Grid *grid, OnePointVariable *onePointVariable, RungeKuttaIndex *indexRK, Data_Param *particleParam)
{
    int nDimVar;

    //! Get the variable of particle.
    RDouble particleDiameter = onePointVariable->GetOnePointDiameter();
    RDouble particleDensity = onePointVariable->GetOnePointDensity();
    SPDouble &particleVelocity = *(onePointVariable->GetOnePointVelocity());
    SPDouble &particleCoordinate = *(onePointVariable->GetOnePointCoordinate());
    SPDouble &particleAngularVelocity = *(onePointVariable->GetOnePointAngularVelocity());
    SPDouble &particleTemperature = *(onePointVariable->GetOnePointTemperature());
    SPDouble &particleForce = *(onePointVariable->GetOnePointForce());

    RDouble particleMass = (4.0 / 3.0) * PI * pow(particleDiameter * 0.5, 3.0) * particleDensity;
    RDouble physicalTimeStep = GlobalDataBase::GetDoubleParaFromDB("physicalTimeStep");

    //! RK variable.
    using namespace INDEX_RK;
    SPDouble &initRKValue = *(onePointVariable->GetOnePointInitRKValue());
    SPDouble &dRKValue = *(onePointVariable->GetOnePointDRKValue());
    SPDouble &sumDRK = *(onePointVariable->GetOnePointSumDRK());

    //! Reset rk variable without initRK.
    nDimVar = 10;
    dRKValue.CopyData(0.0, nDimVar);

    //! Coordinate.
    dRKValue[RK_COORD + 0] = particleVelocity[0];
    dRKValue[RK_COORD + 1] = particleVelocity[1];
    dRKValue[RK_COORD + 2] = particleVelocity[2];

    //! velocity.
    using namespace INDEX_PARTICLE_FORCE;
    nDimVar = 3;
    for (int iForce = 0; iForce < nForceDir / nDimVar; ++iForce)
    {
        dRKValue[RK_VEL + 0] += particleForce[iForce * nDimVar + 0] / particleMass;
        dRKValue[RK_VEL + 1] += particleForce[iForce * nDimVar + 1] / particleMass;
        dRKValue[RK_VEL + 2] += particleForce[iForce * nDimVar + 2] / particleMass;
    }

    nDimVar = 10;
    for (int iDim = 0; iDim < nDimVar; ++iDim)
    {
        sumDRK[iDim] += (indexRK->GetKindex()) * dRKValue[iDim];
    }

    ostringstream ossLog;
    ossLog << " --- init rk 1 ---" << "\n";
    ossLog << " outnstep " << GlobalDataBase::GetIntParaFromDB("outnstep") << "\n";
    ossLog << " particle ID : " << onePointVariable->GetParticleID() << "\n";
    ossLog << " iRK : " << indexRK->GetIndex() << "\n";
    ossLog << " initRKValue coord : ";
    ossLog << initRKValue[RK_COORD + 0] << " " << initRKValue[RK_COORD + 1] << " " << initRKValue[RK_COORD + 2] << "\n";
    ossLog << " initRKValue vel : ";
    ossLog << initRKValue[RK_VEL + 0] << " " << initRKValue[RK_VEL + 1] << " " << initRKValue[RK_VEL + 2] << "\n";
    ossLog << " dRKValue coord :";
    ossLog << dRKValue[RK_COORD + 0] << " " << dRKValue[RK_COORD + 1] << " " << dRKValue[RK_COORD + 2] << "\n";
    ossLog << " dRKValue vel : ";
    ossLog << dRKValue[RK_VEL + 0] << " " << dRKValue[RK_VEL + 1] << " " << dRKValue[RK_VEL + 2] << "\n";
    ossLog << "\n";

    nDimVar = 3;

    if (1 == onePointVariable->GetIfOneRK())
    {
        //! RK1 
        for (int iDim = 0; iDim < nDimVar; ++iDim)
        {
            particleCoordinate[iDim] = initRKValue[RK_COORD + iDim];
            particleVelocity[iDim] = initRKValue[RK_VEL + iDim];
            particleAngularVelocity[iDim] = initRKValue[RK_ANGULAR_VEL + iDim];
        }
        particleTemperature[0] = initRKValue[RK_TEMP];

        //! Old particle cell ID.
        int oldCellID = onePointVariable->GetOldCellID();

        //! Cunrretn index.
        int cellID = onePointVariable->GetCellID();

        onePointVariable->SetCellID(oldCellID);

        onePointVariable->SetCoordinateOld();

        //! Particle interpolation.
        ActionKey *actkeyParticleInterpolation = new ActionKey();
        FillActionKey(actkeyParticleInterpolation, FLOW_ON_PARTICLE_MPI, 0);
        actkeyParticleInterpolation->subact = 1;
        ParticleInterpolation(actkeyParticleInterpolation, grid, onePointVariable);
        delete actkeyParticleInterpolation;

        ParticleForceModel::CalcParticleForce(onePointVariable, particleParam);

        for (int iDim = 0; iDim < nDimVar; ++iDim)
        {
            particleCoordinate[iDim] = initRKValue[RK_COORD + iDim] + physicalTimeStep * (indexRK->GetFindex()) * dRKValue[RK_COORD + iDim];
            particleVelocity[iDim] = initRKValue[RK_VEL + iDim] + physicalTimeStep * (indexRK->GetFindex()) * dRKValue[RK_VEL + iDim];
            particleAngularVelocity[iDim] = initRKValue[RK_ANGULAR_VEL + iDim] + physicalTimeStep * (indexRK->GetFindex()) * dRKValue[RK_ANGULAR_VEL + iDim];
        }
        particleTemperature[0] = initRKValue[RK_TEMP] + physicalTimeStep * (indexRK->GetFindex()) * dRKValue[RK_TEMP];
    }
    else
    {
        ostringstream ossError;
        ossError << " error on rk1 " << onePointVariable->GetIfOneRK() << endl;
        TK_Exit::ExceptionExit(ossError);
    }

    ossLog << " dRKValue coord :" << dRKValue[RK_COORD + 0] << " " << dRKValue[RK_COORD + 1] << " " << dRKValue[RK_COORD + 2] << "\n";
    ossLog << " dRKValue vel : " << dRKValue[RK_VEL + 0] << " " << dRKValue[RK_VEL + 1] << " " << dRKValue[RK_VEL + 2] << "\n";
    ossLog << " particleCoordinate : ";
    ossLog << particleCoordinate[0] << " " << particleCoordinate[1] << " " << particleCoordinate[2] << "\n";
    ossLog << " particleVelocity : ";
    ossLog << particleVelocity[0] << " " << particleVelocity[1] << " " << particleVelocity[2] << "\n";
    ossLog << " cellID : " << onePointVariable->GetCellID() << "\n";
    ossLog << " (indexRK->GetFindex()) " << (indexRK->GetFindex()) << "\n";
    ossLog << " (indexRK->GetKindex()) " << (indexRK->GetKindex()) << "\n";
    ossLog << " physicalTimeStep : " << physicalTimeStep << "\n";
    ossLog << " sumDRK coord : ";
    ossLog << sumDRK[RK_COORD + 0] << " " << sumDRK[RK_COORD + 1] << " " << sumDRK[RK_COORD + 2] << "\n";
    ossLog << " sumDRK vel : ";
    ossLog << sumDRK[RK_VEL + 0] << " " << sumDRK[RK_VEL + 1] << " " << sumDRK[RK_VEL + 2] << "\n";
    ossLog << "\n" << endl;

    bool writeEachProcessor = false;
}

void ParticlePointSolverParallelDomain::CheckParticleBCForOne(SPDouble &oneParticleCoordinate)
{
    int nDim = GetDim();
    Param_ParticlePointSolverParallelDomain *parameters = GetControlParameters();

    double *domainMesh;
    domainMesh = parameters->GetDomainMesh();

    RDouble lengthPar2Domain;
    RDouble lengthDomain;
    int nLength;

    for (int iDim = 0; iDim < nDim; ++iDim)
    {
        lengthDomain = abs(domainMesh[iDim * 2 + 1] - domainMesh[iDim * 2]);

        if (oneParticleCoordinate[iDim] <= domainMesh[iDim * 2])
        {
            lengthPar2Domain = abs(domainMesh[iDim * 2] - oneParticleCoordinate[iDim]);
            nLength = floor(lengthPar2Domain / lengthDomain);
            lengthPar2Domain -= nLength * lengthDomain;

            oneParticleCoordinate[iDim] = domainMesh[iDim * 2 + 1] - lengthPar2Domain;
        }
        else if (oneParticleCoordinate[iDim] > domainMesh[iDim * 2 + 1])
        {
            lengthPar2Domain = abs(oneParticleCoordinate[iDim] - domainMesh[iDim * 2 + 1]);
            nLength = floor(lengthPar2Domain / lengthDomain);
            lengthPar2Domain -= nLength * lengthDomain;

            oneParticleCoordinate[iDim] = domainMesh[iDim * 2] + lengthPar2Domain;
        }
    }
}

void ParticlePointSolverParallelDomain::CommunicationParticle(Grid *grid, FieldProxy *rhsProxy)
{
    using namespace PHMPI;
    ParticlePointGroup *particlePointGroup = GetParticlePointGroup("particlePointGroup");

    //! The number of particle on current zone.
    int nLocalParticle = particlePointGroup->GetNumOfLocalParticle();

    //! The list of particleID on current zone.
    int *particleIDList = particlePointGroup->GetGlobalParticleIDOnLocalZone();

    //! The global zone id.
    int zoneIDGlobal = grid->GetZoneID();
    int nZonesGlobal = GetNumberofGlobalZones();

    //! The local zone id on current processer.
    int zoneIDLocal = grid->GetZoneLocalID();
    int nZonesLocal = GetNumberofLocalZones();

    int gridType = grid->Type();

    if (gridType == PHSPACE::STRUCTGRID)
    {
        StructGrid *structgrid = StructGridCast(grid);
        CommunicationParticle(structgrid, particlePointGroup, rhsProxy);
    }
    else if (gridType == PHSPACE::UNSTRUCTGRID)
    {
        UnstructGrid *unstructgrid = UnstructGridCast(grid);
        CommunicationParticle(unstructgrid, particlePointGroup, rhsProxy);
    }
    else if (gridType == PHSPACE::MIXGRID)
    {
        //! Mixing
        TK_Exit::UnexpectedVarValue(" gridType = ", PHSPACE::MIXGRID);
    }
}

void ParticlePointSolverParallelDomain::CalcNumOfParticleToSendStruct(FieldProxy *rhsProxy, ActionKey *action, DataContainer *cdataIn, int &nParticleIn, int &size_in, DataContainer *cdataOut, int &nParticleToSend, int &size_send, set<int> &sendIDSet)
{

}

void ParticlePointSolverParallelDomain::DecompressParticleDataContainerStruct(DataContainer *cdata, int &iGlobalzone, int &particleID, int &iP, int &jP, int &kP, OnePointVariable *onePointVariable)
{

    cdata->Read(&particleID, sizeof(int));

    cdata->Read(&iGlobalzone, sizeof(int));

    cdata->Read(&iP, sizeof(int));
    cdata->Read(&jP, sizeof(int));
    cdata->Read(&kP, sizeof(int));

    int nDimVar = 12;
    RDouble *data = new RDouble[nDimVar];

    for (int idata = 0; idata < nDimVar; ++idata)
    {
        cdata->Read(&data[idata], sizeof(RDouble));
    }

    SPDouble &onePointCoordinate = *(onePointVariable->GetOnePointCoordinate());
    SPDouble &onePointVelocity = *(onePointVariable->GetOnePointVelocity());
    SPDouble &onePointAcceleration = *(onePointVariable->GetOnePointAcceleration());
    SPDouble &onePointAngularVelocity = *(onePointVariable->GetOnePointAngularVelocity());

    int index = -1;
    int indexDim = -1;
    onePointCoordinate[++indexDim] = data[++index];
    onePointCoordinate[++indexDim] = data[++index];
    onePointCoordinate[++indexDim] = data[++index];

    indexDim = -1;
    onePointVelocity[++indexDim] = data[++index];
    onePointVelocity[++indexDim] = data[++index];
    onePointVelocity[++indexDim] = data[++index];

    indexDim = -1;
    onePointAcceleration[++indexDim] = data[++index];
    onePointAcceleration[++indexDim] = data[++index];
    onePointAcceleration[++indexDim] = data[++index];

    indexDim = -1;
    onePointAngularVelocity[++indexDim] = data[++index];
    onePointAngularVelocity[++indexDim] = data[++index];
    onePointAngularVelocity[++indexDim] = data[++index];

}

void ParticlePointSolverParallelDomain::CommunicationParticle(StructGrid *grid, ParticlePointGroup *particlePointGroup, FieldProxy *rhsProxy)
{
 
}

void ParticlePointSolverParallelDomain::CommunicationParticle(UnstructGrid *grid, ParticlePointGroup *particlePointGroup, FieldProxy *rhsProxy)
{

}

void ParticlePointSolverParallelDomain::AddParticle(FieldProxy *rhsProxy,StructGrid *grid, int &particleID, int &iP, int &jP, int &kP, OnePointVariable *onePointVariable)
{
    int nI, nJ, nK;
    // int iP, jP, kP;
    grid->GetND(nI, nJ, nK);
    nI -= 1;
    nJ -= 1;
    nK -= 1;
    if (TWO_D == GetDim())
    {
        nK = 1;
    }

    if (-1 == iP)
    {
        iP = nI + 1;
    }

    if (-1 == jP)
    {
        jP = nJ + 1;
    }

    if (-1 == kP)
    {
        kP = nK + 1;
    }

    ParticlePointGroup *particlePointGroup = GetParticlePointGroup("particlePointGroup");
    IndexCell *indexCell = particlePointGroup->GetIndexCell();

    int IDCellOfParticle;
    GetIDCellOfParticleByIJKSturct(nI,nJ,nK,iP,jP,kP,2, IDCellOfParticle);

    onePointVariable->SetCellID(IDCellOfParticle);
    indexCell->AddParticleID(IDCellOfParticle,particleID);

    particlePointGroup->AddParticle2Map(particleID, onePointVariable);

    //! residual.
    LDouble *particleRes = GetLagVariableDouble("particleRes");
    particleRes->AddParticle(particleID);

    LDouble *dParticleRes = GetLagVariableDouble("dParticleRes");
    dParticleRes->AddParticle(particleID);

    LDouble *ifOneRK = GetLagVariableDouble("ifOneRK");
    ifOneRK->AddParticle(particleID);

    //! flow on particle.
    LDouble *flowOnParticle = GetLagVariableDouble("flowOnParticle");
    flowOnParticle->AddParticle(particleID);

    //! The Reynolds number of particle.
    LDouble *particleRe = GetLagVariableDouble("particleRe");
    particleRe->AddParticle(particleID);
}

bool ParticlePointSolverParallelDomain::ParticleBoundary(Grid *gridIn, OnePointVariable *onePointVariable, IndexCell *indexCell)
{
    using namespace PARTICLEBCSPACE;
    int gridType = gridIn->Type();

    int particleID = onePointVariable->GetParticleID();
    int cellID = onePointVariable->GetCellID();
    map<int, IndexParticleOfCell*> *cellParticleID = (indexCell->GetCellParticleID());

    //! wall condition.
    map<int, IndexParticleOfCell*>::iterator iterCell = cellParticleID->find(cellID);
    if (iterCell != cellParticleID->end())
    {
        IndexParticleOfCell *indexParticleOfCell = iterCell->second;
        map<int, int> *bcIndex = indexParticleOfCell->GetBCType();
        bool changCellID;
        for (map<int, int>::iterator iterFace = bcIndex->begin(); iterFace != bcIndex->end(); ++iterFace)
        {
            int iBCFace = iterFace->first;
            int bcType = iterFace->second;
            if (bcType == FULLY_ELASTIC_COLLISION_WALL)
            {
                if (gridType == PHSPACE::STRUCTGRID)
                {
                    StructGrid *structgrid = StructGridCast(gridIn);
                    changCellID = FullyElasticCollisionWallStrcut(structgrid, onePointVariable, indexParticleOfCell, iBCFace);
                }
                else if (gridType == PHSPACE::UNSTRUCTGRID)
                {
                    UnstructGrid *unstructgrid = UnstructGridCast(gridIn);
                }
            }
        }
    }

    //! Here we should decide whether to remove 
    //! the particles and then determine the other boundaries.
    bool inZone = onePointVariable->IsParticleInZone();
    bool ifRemoveParticle = false;
    if (!inZone)
    {
        if (iterCell != cellParticleID->end())
        {
            IndexParticleOfCell *indexParticleOfCell = iterCell->second;
            map<int, int> *bcIndex = indexParticleOfCell->GetBCType();
            int cellType = indexParticleOfCell->GetCellType();
            bool changCellID;
            for (map<int, int>::iterator iterFace = bcIndex->begin(); iterFace != bcIndex->end(); ++iterFace)
            {
                int iBCFace = iterFace->first;
                int bcType = iterFace->second;
                if (REMOVE_PARTICLE == bcType)
                {
                    using namespace INDEX_CELLTYPE;
                    if ( GHOST_OUT == cellType || CORNER_OUT == cellType)
                    {
                        ifRemoveParticle = true;
                    }
                }
            }
        }
    }

    if (ifRemoveParticle)
    {
        return ifRemoveParticle;
    }

    //! Others condition,as remove_as_init_random
    bool ifResetParticle = false;
    bool checkCornerBC1 = false;
    bool checkCornerBC2 = false;
    if (!inZone)
    {
        if (iterCell != cellParticleID->end())
        {
            IndexParticleOfCell *indexParticleOfCell = iterCell->second;
            map<int, int> *bcIndex = indexParticleOfCell->GetBCType();
            int cellType = indexParticleOfCell->GetCellType();
            bool changCellID;
            for (map<int, int>::iterator iterFace = bcIndex->begin(); iterFace != bcIndex->end(); ++iterFace)
            {
                int iBCFace = iterFace->first;
                int bcType = iterFace->second;
                if (REMOVE_PARTICLE_AS_INIT == bcType ||  REMOVE_PARTICLE_AS_INIT_ONLY_COORD == bcType)
                {
                    //! Note that for a Ghost_out cell, the boundary conditions cannot 
                    //! be both REMOVE_PARTICLE_AS_INIT and REMOVE_PARTICLE_AS_INIT_ONLY_COORD.
                    //! But for corner,it seem to be complicated, 
                    //! and we need to save the process so try 
                    //! not to use both of these bounds at the adjacent corner locations.
                    using namespace INDEX_CELLTYPE;
                    if ( GHOST_OUT == cellType || CORNER_OUT == cellType)
                    {
                        ifResetParticle = true;

                        if (CORNER_OUT == cellType &&  REMOVE_PARTICLE_AS_INIT == bcType)
                        {
                            checkCornerBC1 = true;
                        }
                        if (CORNER_OUT == cellType && REMOVE_PARTICLE_AS_INIT_ONLY_COORD == bcType)
                        {
                            checkCornerBC2 = true;
                        }
                    }
                }
            }
        }
    }

    if (ifResetParticle)
    {

        if (true == checkCornerBC1  && true == checkCornerBC2)
        {
            ostringstream ossError;
            ossError << " error both REMOVE_PARTICLE_AS_INIT and REMOVE_PARTICLE_AS_INIT_ONLY_COORD for corner" << endl;
            TK_Exit::ExceptionExit(ossError);
        }

        ResetParticleWhenOut(gridIn, onePointVariable,iterCell->second);

    }

    return ifRemoveParticle;
}


bool ParticlePointSolverParallelDomain::FullyElasticCollisionWallStrcut(StructGrid *grid, OnePointVariable *onePointVariable, IndexParticleOfCell *indexParticleOfCell, int &iBCFace)
{

    int cellID = onePointVariable->GetCellID();

    //! Get index IJK of particle position in current zone.
    int nI, nJ, nK;
    int iP, jP, kP;
    grid->GetND(nI, nJ, nK);
    nI -= 1;
    nJ -= 1;
    nK -= 1;
    if (TWO_D == GetDim())
    {
        nK = 1;
    }
    GetIndexCellIJKStruct(nI, nJ, nK, iP, jP, kP, GetNumberOfGhostCellLayers(), cellID);

    //! Index start and end of cell in current zone.
    int ist, ied;
    int jst, jed;
    int kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    //! Index start and end of surface in .
    //! From 1 to nDim.
    int nDim = GetDim();
    int iSurfaceStart = 1;
    int iSurfaceEnd = nDim;

    //! Get cell center.
    RDouble3D &xcc = *(grid->GetCellCenterX());
    RDouble3D &ycc = *(grid->GetCellCenterY());
    RDouble3D &zcc = *(grid->GetCellCenterZ());

    RDouble3D &xNode = *(grid->GetStructX());
    RDouble3D &yNode = *(grid->GetStructY());
    RDouble3D &zNode = *(grid->GetStructZ());

    //! Get face normal.
    RDouble4D &xfn = *(grid->GetFaceNormalX());
    RDouble4D &yfn = *(grid->GetFaceNormalY());
    RDouble4D &zfn = *(grid->GetFaceNormalZ());

    int iSurLocal, jSurLocal, kSurLocal;

    using namespace INDEX_CELLTYPE;
    int cellType = indexParticleOfCell->GetCellType();

    //! iBCFace in 2D : i:0,1, j:2,3.
    //! iBCFace in 3D : i:0,1, j:2,3, k:4,5
    //! iDir : thr direction of face in i,j,k.
    //! ilr : Left or right of boundary condition face (0 or 1).
    //! ilrDir : seem to ilr, 1 left, -1 right, 
    //!          which is to make the face normal always point into the cell.
    int iDir = 0;
    int ilr = 0;
    int ilrDir = 0;
    for (int iDim = 1; iDim <= GetDim(); ++iDim)
    {
        //! The direction of boundary condition face.
        if (iBCFace <= 2 * iDim - 1 && iBCFace > 2 * (iDim-1) - 1)
        {
            iDir = iDim;
            ilr = iBCFace - (2 * (iDim - 1) - 1);
            //! ilr : 0,left, 1,right.
            ilr = ilr - 1;
            if (ilr < 0 || ilr > 1)
            {
                ostringstream ossError;
                ossError << " error on ilr = " << ilr << endl;
                TK_Exit::ExceptionExit(ossError);
            }
            else
            {
                if (0 == ilr)
                {
                    ilrDir = 1;
                }
                else
                {
                    ilrDir = -1;
                }
            }
        }
    }

    //! Get the surface (or parallel) index.
    GetNsurfIndex(iDir, iSurLocal, jSurLocal, kSurLocal);
    iSurLocal = iP + ilr * iSurLocal;
    jSurLocal = jP + ilr * jSurLocal;
    kSurLocal = kP + ilr * kSurLocal;

    const int nDimVar = 3;
    RDouble fn[nDimVar];
    fn[0] = xfn(iSurLocal, jSurLocal, kSurLocal, iDir) * ilrDir;
    fn[1] = yfn(iSurLocal, jSurLocal, kSurLocal, iDir) * ilrDir;
    fn[2] = zfn(iSurLocal, jSurLocal, kSurLocal, iDir) * ilrDir;
    double lengthOfNormalVector = DISTANCE(fn[0], fn[1], fn[2]);

    for (int iDim = 0; iDim < nDimVar; ++iDim)
    {
        fn[iDim]= fn[iDim] / lengthOfNormalVector;
    }
    //! The center coordinate of surface.
    RDouble fc[nDimVar];
    //! Get the center coodinate of surface.
    grid->FaceCoor(iSurLocal, jSurLocal, kSurLocal, iDir, fc[0], fc[1], fc[2]);

    RDouble particleDiameter = onePointVariable->GetOnePointDiameter();
    SPDouble &particleCoordinate = *(onePointVariable->GetOnePointCoordinate());
    SPDouble &particleVelocity = *(onePointVariable->GetOnePointVelocity());

    RDouble parVel[nDimVar];
    RDouble distanceParticleToFace = 0.0;
    RDouble surEqD = -fn[0] * fc[0] - fn[1] * fc[1] - fn[2] * fc[2];
    for (int iDim = 0; iDim < nDimVar; ++iDim)
    {
        parVel[iDim] = particleVelocity[iDim];
        distanceParticleToFace += fn[iDim] * particleCoordinate[iDim];
    }
    distanceParticleToFace += surEqD;
    distanceParticleToFace = ABS(distanceParticleToFace);
    distanceParticleToFace /= DISTANCE(fn[0], fn[1], fn[2]);

    int outIterStep = GlobalDataBase::GetIntParaFromDB("outnstep");
    fstream f;

    //! Notice, this should change the velocity first,
    //! and then should change the coordinates.
    if ( GHOST_IN == cellType || CORNER_IN == cellType)
    {

        RDouble physicalTimeStep = GlobalDataBase::GetDoubleParaFromDB("physicalTimeStep");
        RDouble refv = GlobalDataBase::GetDoubleParaFromDB("refDimensionalVelocity");

        RDouble velNDot = 0;
        for (int iDim = 0; iDim < nDimVar; ++iDim)
        {
            velNDot += fn[iDim] * particleVelocity[iDim];
        }

        RDouble displacementN = 0;
        RDouble checkDistance = particleDiameter*0.5;

        if (GHOST_IN == cellType)
        {
            if (distanceParticleToFace > checkDistance)
            {
                if (velNDot < 0)
                {
                    displacementN = abs(velNDot);
                    displacementN *= physicalTimeStep;
                    if (displacementN > (distanceParticleToFace + checkDistance))
                    {
                        RDouble parVelNew[nDimVar] = { 0 };

                        KernelBoundaryTransformation(fn, parVel, parVelNew);
                        f.open("./results/ParticleForce.txt", ios::out | ios::app);
                        f << outIterStep <<" " <<particleDiameter << " " << parVelNew[0]-parVel[0] << " " << parVelNew[1] - parVel[1] <<" "<< parVelNew[2] - parVel[2] << endl;
                        f.close();
                        for (int iDim = 0; iDim < nDimVar; ++iDim)
                        {
                            particleVelocity[iDim] = parVelNew[iDim];
                        }
                    }
                }
            }
            else
            {
                RDouble parVelNew[nDimVar] = { 0 };
                KernelBoundaryTransformation(fn, parVel, parVelNew);
                f.open("./results/ParticleForce.txt", ios::out | ios::app);
                f << outIterStep << " " << particleDiameter << " " << parVelNew[0] - parVel[0] << " " << parVelNew[1] - parVel[1] << " " << parVelNew[2] - parVel[2] << endl;
                f.close();
                for (int iDim = 0; iDim < nDimVar; ++iDim)
                {
                    particleVelocity[iDim] = parVelNew[iDim];
                }
            }
        }
        if (CORNER_IN == cellType)
        {
            if (distanceParticleToFace > checkDistance)
            {

                if (velNDot < 0)
                {
                    displacementN = abs(velNDot);
                    displacementN *= physicalTimeStep;
                    if (displacementN > (distanceParticleToFace+ checkDistance))
                    {
                        RDouble parVelNew[nDimVar] = { 0 };

                        KernelBoundaryTransformation(fn, parVel, parVelNew);
                        f.open("./results/ParticleForce.txt", ios::out | ios::app);
                        f << outIterStep << " " << particleDiameter << " " << parVelNew[0] - parVel[0] << " " << parVelNew[1] - parVel[1] << " " << parVelNew[2] - parVel[2] << endl;
                        f.close();
                        for (int iDim = 0; iDim < nDimVar; ++iDim)
                        {
                            particleVelocity[iDim] = parVelNew[iDim];
                        }
                    }
                }
            }
            else
            {
                RDouble parVelNew[nDimVar] = { 0 };
                KernelBoundaryTransformation(fn, parVel, parVelNew);

                f.open("./results/ParticleForce.txt", ios::out | ios::app);
                f << outIterStep << " " << particleDiameter << " " << parVelNew[0] - parVel[0] << " " << parVelNew[1] - parVel[1] << " " << parVelNew[2] - parVel[2] << endl;
                f.close();
                for (int iDim = 0; iDim < nDimVar; ++iDim)
                {
                    particleVelocity[iDim] = parVelNew[iDim];
                }
            }
        }
    }
    else if (GHOST_OUT == cellType || CORNER_OUT == cellType)
    {

        RDouble checkDistance = particleDiameter * 0.5;
      
        cout << "###WARNING:particle ejected into GHOST GRID!### " << endl;
        SPDouble &particleCoordinateOld = *(onePointVariable->GetOnePointCoordinateOld());
       
        for (int iDim = 0; iDim < nDimVar; ++iDim)
        {
            particleCoordinate[iDim] = particleCoordinateOld[iDim];
        }
        onePointVariable->SetCellID(onePointVariable->GetOldCellID());
}

    bool changeCellID = false;

    return changeCellID;
}

void ParticlePointSolverParallelDomain::ResetParticleWhenOut(Grid *grid, OnePointVariable *onePointVariable, IndexParticleOfCell *indexParticleOfCell)
{
    using namespace PARTICLEBCSPACE;
    int bcType = 0;

    map<int, int> *bcIndex = indexParticleOfCell->GetBCType();
    int cellType = indexParticleOfCell->GetCellType();
    bool changCellID;
    for (map<int, int>::iterator iterFace = bcIndex->begin(); iterFace != bcIndex->end(); ++iterFace)
    {
        int iBCFace = iterFace->first;
        int bcFaceType = iterFace->second;
        if (REMOVE_PARTICLE_AS_INIT == bcFaceType)
        {
            bcType = bcFaceType;
        }
        else if(REMOVE_PARTICLE_AS_INIT_ONLY_COORD == bcFaceType)
        {
            bcType = bcFaceType;
        }
    }

    //! Get particleBCInit.
    ParticlePointGroup *particleBCInit = GetParticlePointGroup("particleBCInit");

    //! Fine the init variable.
    OnePointVariable *onePointBCInit = particleBCInit->GetOnePoint(onePointVariable->GetParticleID());

    SPDouble &particleCoordinate = *(onePointVariable->GetOnePointCoordinate());
    SPDouble &particleVelocity = *(onePointVariable->GetOnePointVelocity());
    SPDouble &particleAngularVelocity = *(onePointVariable->GetOnePointAngularVelocity());
    SPDouble &particleAcceleration = *(onePointVariable->GetOnePointAcceleration());
    SPDouble &particleTemperature = *(onePointVariable->GetOnePointTemperature());
    SPDouble &flowOnParticle = *(onePointVariable->GetOnePointFlowOnParticle());
    SPDouble &particleForce = *(onePointVariable->GetOnePointForce());

    int nDimVar;

    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();
    int currentProcessorID = GetCurrentProcessorID();
    //! The global zone id.
    int zoneIDGlobal = grid->GetZoneID();
    //! The local zone id on current processer.
    int zoneIDLocal = grid->GetZoneLocalID();
    ostringstream ossLog;
    ossLog << " --- check reset particle when out --- " << "\n";
    ossLog << " zone ID Global :" << zoneIDGlobal << "\n";
    ossLog << " zone ID Local :" << zoneIDLocal << endl;
    ossLog << " currentProcessorID : " << currentProcessorID << "\n";
    ossLog << " bcType : " << bcType << "\n";
    ossLog << " particle ID : " << onePointVariable->GetParticleID() << "\n"; 
    ossLog << " old value of particle : " << "\n";
    ossLog << "  cellID : " << onePointVariable->GetCellID() << "\n";
    ossLog << (onePointVariable->Print2WindowOS())<< "\n";

    if (REMOVE_PARTICLE_AS_INIT == bcType)
    {
        onePointVariable->SetCellID(onePointBCInit->GetCellID());

        nDimVar = 3;
        particleCoordinate.CopyData(*(onePointBCInit->GetOnePointCoordinate()), nDimVar);
        particleVelocity.CopyData(*(onePointBCInit->GetOnePointVelocity()), nDimVar);
        particleAngularVelocity.CopyData(*(onePointBCInit->GetOnePointAngularVelocity()), nDimVar);
        particleAcceleration.CopyData(*(onePointBCInit->GetOnePointAcceleration()), nDimVar);

        nDimVar = 1;
        particleTemperature.CopyData(*(onePointBCInit->GetOnePointTemperature()), nDimVar);

    }
    else if (REMOVE_PARTICLE_AS_INIT_ONLY_COORD == bcType)
    {
        onePointVariable->SetCellID(onePointBCInit->GetCellID());

        nDimVar = 3;
        particleCoordinate.CopyData(*(onePointBCInit->GetOnePointCoordinate()), nDimVar);

        ActionKey *actkeyParticleInterpolation = new ActionKey();
        FillActionKey(actkeyParticleInterpolation, FLOW_ON_PARTICLE_MPI, 0);
        actkeyParticleInterpolation->subact = 0;
        ParticleInterpolation(actkeyParticleInterpolation, grid, onePointVariable);
        delete actkeyParticleInterpolation;
    }

    ossLog << " new value of particle : " << "\n";
    ossLog << "  cellID : " << onePointVariable->GetCellID() << "\n";
    ossLog << (onePointVariable->Print2WindowOS()) << "\n";
    WriteLogFile(ossLog);
}

RDouble ParticlePointSolverParallelDomain::GetMaxParIDLocal()
{
    ParticlePointGroup *particlePointGroup = GetParticlePointGroup("particlePointGroup");
    set<int> parIDSet;
    for (ParticleIDIter iterPar = particlePointGroup->GetIterBegin();iterPar != particlePointGroup->GetIterEnd(); ++iterPar)
    {
        parIDSet.insert(iterPar->first);
    }
    int maxParIDLocal = *(parIDSet.rbegin());
    return maxParIDLocal;
}

RDouble ParticlePointSolverParallelDomain::GetRandomNumber(RDouble ValueOfLeftOrMean, RDouble ValueOfRightOrStandardDeviation,int mode)
{
    std::random_device Gen;
        RDouble randomValue;
    std::default_random_engine Number;
    Number.seed(Gen());
    if (1 == mode)
    {
        std::uniform_real_distribution<double> uniformNumber(ValueOfLeftOrMean, ValueOfRightOrStandardDeviation);
        randomValue = uniformNumber(Number);
    }

    if (2 == mode)
    {
        std::normal_distribution<double> normalDistribution(ValueOfLeftOrMean, ValueOfRightOrStandardDeviation);
        randomValue = normalDistribution(Number);
    }   

    return randomValue;
}

void ParticlePointSolverParallelDomain::AddParticleAsInitBC(Grid *grid)
{
    Param_ParticlePointSolverParallelDomain *parameters = GetControlParameters();

    //! Get particleBCInit.
    ParticlePointGroup *particleBCInit = GetParticlePointGroup("particleBCInit");
    ParticlePointGroup *particlePointGroup = GetParticlePointGroup("particlePointGroup");

    using namespace PHMPI;
    //! The number of zones on global all processer.
    int nZonesGlobal = GetNumberofGlobalZones();
    int zoneIDLocal = grid->GetZoneLocalID();
    //! Get current processor ID on current zones..
    int currentProcessorID = GetCurrentProcessorID();
    int maxParIDLocal = 0;
    int maxParID = 0;
    if (0 == zoneIDLocal)
    {
        for (int iZoneGlobal = 0; iZoneGlobal < nZonesGlobal; ++iZoneGlobal)
        {
            //! Judge the zones on current processor.
            int recvProcess = GetZoneProcessorID(iZoneGlobal);
            //! The loop for zone on current processor.
            if (currentProcessorID == recvProcess)
            {
                //! For each zone'solver.
                int solverID = this->GetIndex();
                ParticlePointSolverParallelDomain *solver = static_cast <ParticlePointSolverParallelDomain*> (GlobalSolvers::GetSolver(iZoneGlobal, solverID));
                
                maxParIDLocal = MAX(solver->GetMaxParIDLocal(), maxParIDLocal);
            }
        }
    }
    PH_Barrier();
    PH_AllReduce(&maxParIDLocal, &maxParID, 1, MPI_MAX);

    //! Random velocity direction and value.
    if (3 == parameters->GetInitParticleBCInfoType())
    {
        int randomMode = parameters->GetParticleRandomReleaseMode();
        int RandomDp = parameters->GetParticleRandomParticleDiameterMode();
        if (1 == randomMode)
        {
            cout << "Particle random release mode is ACTIVATED! Releasing particle ...Dp mode is " << RandomDp << endl;
        }
        else
        {
            cout << "Particle random release mode is SUPPRESSED! Releasing particle ..." << endl;
        }
        RDouble randAngle = parameters->GetInitParticleBCRandomDirectionalRatio();
        randAngle = randAngle * PI/180;
        RDouble valueRatio = parameters->GetInitParticleBCRandomValueRatio();

        int countPar = 0;
        map<int, OnePointVariable*>::iterator iter = particleBCInit->GetIterBegin();
        while (iter != particleBCInit->GetIterEnd())
        {
            OnePointVariable *onePointBCInit = iter->second;
            SPDouble &particleVelocity = *(onePointBCInit->GetOnePointVelocity());
            if (1 == randomMode)
            {
                OnePointVariable *onePointVariable = new OnePointVariable();

                RDouble Dp = onePointBCInit->GetOnePointDiameter();

                if (1 == RandomDp) 
                {
                    RDouble MaxDp = parameters->GetRandomDiameterMaxRatio();
                    RDouble MinDp = parameters->GetRandomDiameterMinRatio();
                    MaxDp *= Dp;
                    MinDp *= Dp;
                    Dp = GetRandomNumber(MinDp, MaxDp, 1);
                 
                 }
                else if (2 == RandomDp)
                {
                    RDouble VarianceDp = parameters->GetRandomDiameterMaxRatio();
                    RDouble MeanDp = parameters->GetRandomDiameterMinRatio();

                    Dp = GetRandomNumber(MeanDp, VarianceDp, 2);
                    Dp*= 1e-6;

                }
                onePointVariable->SetDiameter(abs(Dp));
                onePointVariable->SetDensity(onePointBCInit->GetOnePointDensity());
                onePointVariable->SetSpecificHeatCapacity(onePointBCInit->GetOnePointSpecificHeatCapacity());
                onePointVariable->SetEmissivity(onePointBCInit->GetOnePointEmissivity());

                int nDimVar = 3;
                onePointVariable->GetOnePointCoordinate()->CopyData(*(onePointBCInit->GetOnePointCoordinate()), nDimVar);
                onePointVariable->GetOnePointVelocity()->CopyData(*(onePointBCInit->GetOnePointVelocity()), nDimVar);
                onePointVariable->GetOnePointAngularVelocity()->CopyData(*(onePointBCInit->GetOnePointAngularVelocity()), nDimVar);
                onePointVariable->GetOnePointAcceleration()->CopyData(*(onePointBCInit->GetOnePointAcceleration()), nDimVar);

                nDimVar = 1;
                onePointVariable->GetOnePointTemperature()->CopyData(*(onePointBCInit->GetOnePointTemperature()), nDimVar);

                //! ##################################################
                //! #    Region Release
                //! ##################################################
                SPDouble &particleCoordinateNew = *(onePointVariable->GetOnePointCoordinate());
                //! Release domin
                //! ##################################################
                RDouble XRange_up = -0.005;
                RDouble XRange_down = -0.005;
                RDouble YRange_up = 1.0*Dp;
                RDouble YRange_down = 0.01;

                particleCoordinateNew[0] = GetRandomNumber(XRange_up, XRange_down, 1);
                particleCoordinateNew[1] = GetRandomNumber(YRange_up, YRange_down, 1);

                if (THREE_D == GetDim())
                {
                    RDouble ZRange_up = -0.01;
                    RDouble ZRange_down = -1.0 * Dp;
                    particleCoordinateNew[2] = GetRandomNumber(ZRange_up, ZRange_down, 1);
                }
                else
                {

                }

                //! if the particle is in current zone.
                //! For example, when particles are at Ghost or Corner 
                //! and need to be transported or removed by MPI or OutBC.
                bool inZone = true;
                //! Note that the search does not include particles 
                //! from the Ghostand corner areas.
                //! Particles in these areas will be accepted into the required zone.
                int idCellOfParticle = 0;

                idCellOfParticle = QuickSearchParticleCell(grid, particleCoordinateNew, inZone);

                onePointBCInit->SetCellID(idCellOfParticle);
                onePointVariable->SetParticleBCType(onePointBCInit->GetParticleBCType());
                onePointVariable->SetCellID(idCellOfParticle);
                onePointVariable->SetIsParticleInZone(inZone);
                countPar += 1;
                int newIndex = maxParID + countPar;
                onePointVariable->SetParticleID(newIndex);
                particlePointGroup->QuickAddParticle(onePointVariable->GetParticleID(), onePointVariable->GetCellID(), onePointVariable, inZone);

                iter++;
            }
            else
            {
                OnePointVariable *onePointVariable = new OnePointVariable();
                onePointVariable->SetCellID(onePointBCInit->GetCellID());
                onePointVariable->SetParticleBCType(onePointBCInit->GetParticleBCType());
                bool isZoneParticle = onePointBCInit->IsParticleInZone();

                onePointVariable->SetIsParticleInZone(isZoneParticle);

                onePointVariable->SetDiameter(onePointBCInit->GetOnePointDiameter());
                onePointVariable->SetDensity(onePointBCInit->GetOnePointDensity());
                onePointVariable->SetSpecificHeatCapacity(onePointBCInit->GetOnePointSpecificHeatCapacity());
                onePointVariable->SetEmissivity(onePointBCInit->GetOnePointEmissivity());

                int nDimVar = 3;
                onePointVariable->GetOnePointCoordinate()->CopyData(*(onePointBCInit->GetOnePointCoordinate()), nDimVar);
                onePointVariable->GetOnePointVelocity()->CopyData(*(onePointBCInit->GetOnePointVelocity()), nDimVar);
                onePointVariable->GetOnePointAngularVelocity()->CopyData(*(onePointBCInit->GetOnePointAngularVelocity()), nDimVar);
                onePointVariable->GetOnePointAcceleration()->CopyData(*(onePointBCInit->GetOnePointAcceleration()), nDimVar);

                nDimVar = 1;
                onePointVariable->GetOnePointTemperature()->CopyData(*(onePointBCInit->GetOnePointTemperature()), nDimVar);

                countPar += 1;
                int newIndex = maxParID + countPar;
                onePointVariable->SetParticleID(newIndex);

                particlePointGroup->QuickAddParticle(onePointVariable->GetParticleID(), onePointVariable->GetCellID(), onePointVariable, isZoneParticle);

                iter++;
            }
        }
    }
    UpdateNumOfParticle();
}

int *ParticlePointSolverParallelDomain::GetCellIndex(Grid *grid, int cellID)
{
    int *cellIndexID = 0;

    int gridType = grid->Type();

    if (PHSPACE::STRUCTGRID == gridType)
    {
        StructGrid *structgrid = StructGridCast(grid);

        int nI, nJ, nK;
        int iP, jP, kP;
        structgrid->GetND(nI, nJ, nK);
        nI -= 1;
        nJ -= 1;
        nK -= 1;
        if (TWO_D == GetDim())
        {
            nK = 1;
        }

        GetIndexCellIJKStruct(nI, nJ, nK, iP, jP, kP, GetNumberOfGhostCellLayers(), cellID);

        cellIndexID = new int[3];
        cellIndexID[0] = iP;
        cellIndexID[1] = jP;
        cellIndexID[2] = kP;
    }
    else if (PHSPACE::UNSTRUCTGRID == gridType)
    {
        UnstructGrid *unstructgrid = UnstructGridCast(grid);
        cellIndexID = new int(1);
        cellIndexID[0] = -1;
    }
    else if (PHSPACE::MIXGRID == gridType)
    {
        //! Mixing
        TK_Exit::UnexpectedVarValue(" gridType = ", PHSPACE::MIXGRID);
    }

    return cellIndexID;
}

void ParticlePointSolverParallelDomain::CalculateVariableForParticleOnEular(ActionKey *actkey)
{
    Grid *grid = GetGrid();
    int gridType = grid->Type();

    if (PHSPACE::STRUCTGRID == gridType)
    {
        StructGrid *structgrid = StructGridCast(grid);
        CalculateVariableForParticleOnEular(structgrid, actkey);
    }
    else if ( PHSPACE::UNSTRUCTGRID == gridType)
    {
        UnstructGrid *unstructgrid = UnstructGridCast(grid);
        CalculateVariableForParticleOnEular(unstructgrid, actkey);
    }
    else if (PHSPACE::MIXGRID == gridType)
    {
        //! Mixing
        TK_Exit::UnexpectedVarValue(" gridType = ", PHSPACE::MIXGRID);
    }
}

void ParticlePointSolverParallelDomain::CalculateVariableForParticleOnEular(StructGrid *grid, ActionKey *actkey)
{
    Param_ParticlePointSolverParallelDomain *parameters = GetControlParameters();

    //! Get the interpolation method of particle on struct grid.
    string particleInterpolation_struct = parameters->GetParticleInterpolationStruct();

    if (particleInterpolation_struct == "Gimenez2019JCP")
    {
        using namespace PHSPACE::INDEX_FLOWONPARTICLE;
        //! Flow variable on particle location.
        //! 1 -- flow velocity.
        //! 2 -- flow density.
        //! 3 -- flow dynamic viscosity.
        //! 3 -- flow temperature.
        //! 5 -- flow velocity grad.

        //! use Taylor series by cell-centered gradient of the field.
        RDouble4D &flowOnParticleInEular = *reinterpret_cast <RDouble4D*> (grid->GetDataPtr("flowOnParticleInEular"));
        RDouble4D &flowOnParticleInEularGrad = *reinterpret_cast <RDouble4D*> (grid->GetDataPtr("flowOnParticleInEularGrad"));
        RDouble4D &flowOnParticleInEularGradSecond = *reinterpret_cast <RDouble4D*> (grid->GetDataPtr("flowOnParticleInEularGradSecond"));

        //! volmu of cell.
        RDouble3D &vol = *(grid->GetCellVolume());

        //! Get face vecttor of cell.
        //! The length of the vector is the magnitude of the area.
        RDouble4D &xfv = *(grid->GetFaceVectorX());
        RDouble4D &yfv = *(grid->GetFaceVectorY());
        RDouble4D &zfv = *(grid->GetFaceVectorZ());

        //! Get dynamic viscous of Laminar.
        RDouble3D &viscousLaminar = *reinterpret_cast <RDouble3D*> (grid->GetDataPtr("visl"));
        RDouble3D &viscousTurbulence = *reinterpret_cast <RDouble3D*> (grid->GetDataPtr("vist"));
        int viscousType = parameters->GetViscousType();    //! Viscous model.

        //! Get temperature
        RDouble4D &t = *reinterpret_cast <RDouble4D*> (grid->GetDataPtr("t"));

        //! Get the variable of flow.
        RDouble4D &q = *reinterpret_cast <RDouble4D*> (grid->GetDataPtr("q"));

        //! Grad of flow.
        RDouble4D &gradUVWTCellCenterX = *reinterpret_cast <RDouble4D*> (grid->GetDataPtr("gradUVWTCellCenterX"));
        RDouble4D &gradUVWTCellCenterY = *reinterpret_cast <RDouble4D*> (grid->GetDataPtr("gradUVWTCellCenterY"));
        RDouble4D &gradUVWTCellCenterZ = *reinterpret_cast <RDouble4D*> (grid->GetDataPtr("gradUVWTCellCenterZ"));

        //! Index start and end of surface in .
        //! From 1 to nDim.
        int nDim = GetDim();
        int iSurfaceStart = 1;
        int iSurfaceEnd = nDim;

        //! The node of grid.
        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();

        //! Index for cell center.
        //! Index start and end of cell in current zone, layer is 2.
        int ist, ied;
        int jst, jed;
        int kst, ked;
        int nLayers = GetNumberOfGhostCellLayers();
        grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);

        int nDimVar = 3;

        //! Note that the scope of our subact three times is not the same. 
        //! The first one contains Ghost and corner, 
        //! the second one contains Ghost and corner on the periphery, 
        //! and the third one only contains the central zone area.

        //! Here, we take a corner in two dimensions 
        //! as an example to illustrate the influence of Corner.
        //! (-1,-1),(-1,0),(0,-1) and (0,0) is the corner cell.
        //! and (-1,1),(0,1),(1,0) and (1,-1) is the ghost cell.
        //! (0,0) is the zone center cell.

        //! Let's say we fill the corner values with 
        //! the program's own fill and decide 
        //! that the corner values are inappropriate (with some error).
        //! Take the variable q, for example.
        //! problem is (-1,-1),(-1,0),(0,-1) and (0,0).
        //! The corner of grad q,
        //! here's the problem is (0,0),(0,1) and (1,0).
        //! For grad second q,
        //! the problem is (0,0).
        //! So, as long as we don't pass the value after the computation, 
        //! we're fine, for example.
        //! We mpi communicates with q after subact1, 
        //! and then communicates with Grad and grad 2 at last.

        if (0 == actkey->subact || 1 == actkey->subact)
        {
            //! Step 1: init value.
            flowOnParticleInEular = 0.0;
            //! Note that this loop index include none layer ghost cell
            grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, nLayers);
            for (int k = kst; k <= ked; ++k)
            {
                for (int j = jst; j <= jed; ++j)
                {
                    for (int i = ist; i <= ied; ++i)
                    {
                        //! flow density.
                        flowOnParticleInEular(i, j, k, PFU) = q(i, j, k, IDX::IU);
                        flowOnParticleInEular(i, j, k, PFV) = q(i, j, k, IDX::IV);
                        flowOnParticleInEular(i, j, k, PFW) = q(i, j, k, IDX::IW);

                        //! flow density.
                        flowOnParticleInEular(i, j, k, PFR) = q(i, j, k, IDX::IR);

                        //! flow dynamic viscosity
                        if (viscousType > LAMINAR)
                        {
                            flowOnParticleInEular(i, j, k, PFMU) = viscousLaminar(i, j, k) + viscousTurbulence(i, j, k);

                        }
                        else
                        {
                            flowOnParticleInEular(i, j, k, PFMU) = viscousLaminar(i, j, k);
                        }

                        //! flow temperature.
                        flowOnParticleInEular(i, j, k, PFT) = t(i, j, k, 0);

                        flowOnParticleInEular(i, j, k, PFDTDX) = gradUVWTCellCenterX(i, j, k, 3);
                        flowOnParticleInEular(i, j, k, PFDTDY) = gradUVWTCellCenterY(i, j, k, 3);
                        flowOnParticleInEular(i, j, k, PFDTDZ) = gradUVWTCellCenterZ(i, j, k, 3);

                        //! Flow grad.
                        flowOnParticleInEular(i, j, k, PFDUDX ) = gradUVWTCellCenterX(i, j, k, 0);
                        flowOnParticleInEular(i, j, k, PFDUDY ) = gradUVWTCellCenterY(i, j, k, 0);
                        flowOnParticleInEular(i, j, k, PFDUDZ ) = gradUVWTCellCenterZ(i, j, k, 0);

                        flowOnParticleInEular(i, j, k, PFDVDX) = gradUVWTCellCenterX(i, j, k, 1);
                        flowOnParticleInEular(i, j, k, PFDVDY) = gradUVWTCellCenterY(i, j, k, 1);
                        flowOnParticleInEular(i, j, k, PFDVDZ) = gradUVWTCellCenterZ(i, j, k, 1);

                        flowOnParticleInEular(i, j, k, PFDWDX) = gradUVWTCellCenterX(i, j, k, 2);
                        flowOnParticleInEular(i, j, k, PFDWDY) = gradUVWTCellCenterY(i, j, k, 2);
                        flowOnParticleInEular(i, j, k, PFDWDZ) = gradUVWTCellCenterZ(i, j, k, 2);
                    }
                }
            }
        }

        if (0 == actkey->subact|| 2 == actkey->subact)
        {
            //! Step 2 : get grad.
            flowOnParticleInEularGrad = 0.0;
            grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, nLayers-1);
            for (int nsurf = 1; nsurf <= nDim; ++nsurf)
            {
                //! For 2D
                //! nsurf = 1, index is (1,0,0)
                //! nsurf = 2, index is (0,1,0)
                //! So the index is the surf index for each cell.
                int il1, jl1, kl1;
                GetNsurfIndex(nsurf, il1, jl1, kl1);

                int ist = 1;
                int ied = ni - 1 + il1;
                int jst = 1;
                int jed = nj - 1 + jl1;
                int kst = 1;
                int ked = nk - 1 + kl1;

                if (TWO_D == nDim) ked = 1;

                int il, jl, kl;

                //! The gradient of u,v,w.
                //! Loop for cell with no layer.
                for (int k = kst; k <= ked; ++k)
                {
                    for (int j = jst; j <= jed; ++j)
                    {
                        for (int i = ist; i <= ied; ++i)
                        {
                            //! i,j,k and nsurf is the nsurf-st surface for cell (i,j,k)
                            //! il,jl,kl is the cell which shard the surf whit cell(i,j,k)

                            il = i - il1;
                            jl = j - jl1;
                            kl = k - kl1;

                            //! Loop for variable.
                            for (int m = 0; m < nVarOnFlow; ++m)
                            {
                                //! Note that xfv does not have value on (0,0,1) and (-1,1,1)
                                //! but is has value on (0,1,1)

                                //! considers a linear variation of the field between 
                                //! the centroids of the cells which share the face.

                                RDouble phis = flowOnParticleInEular(i, j, k, m) + flowOnParticleInEular(il, jl, kl, m);

                                RDouble ddx = phis * xfv(i, j, k, nsurf);
                                RDouble ddy = phis * yfv(i, j, k, nsurf);
                                RDouble ddz = phis * zfv(i, j, k, nsurf);

                                //! For the surf in current cell, the surf is inwards
                                //! but in Gauss formula we need the outwards pointing face area vector.
                                //! The index,eg : 
                                //! 0 - dudx, 1 - dudy, 2 - dudz, PFU*3 + 0,1,2
                                //! 3 - dvdx, 4 - dvdy, 5 - dvdz, PFV*3 + 0,1,2
                                //! 6 - dwdx, 7 - dwdy, 8 - dwdz, PFW*3 + 0,1,2
                                //! So that, drhody is PFR*3+1.
                                int mGrad;
                                mGrad = m * 3 - 1;
                                flowOnParticleInEularGrad(i, j, k, ++mGrad) -= ddx;
                                flowOnParticleInEularGrad(i, j, k, ++mGrad) -= ddy;
                                flowOnParticleInEularGrad(i, j, k, ++mGrad) -= ddz;

                                mGrad = m * 3 - 1;
                                //! These face area vector is outwards pointing.
                                flowOnParticleInEularGrad(il, jl, kl, ++mGrad) += ddx;
                                flowOnParticleInEularGrad(il, jl, kl, ++mGrad) += ddy;
                                flowOnParticleInEularGrad(il, jl, kl, ++mGrad) += ddz;
                            }
                        }
                    }
                }
            }

            for (int k = kst; k <= ked; ++k)
            {
                for (int j = jst; j <= jed; ++j)
                {
                    for (int i = ist; i <= ied; ++i)
                    {
                        //! half is a linear variation.
                        RDouble oov = half / vol(i, j, k);

                        //! Loop for variable.
                        for (int mGrad = 0; mGrad < nVarOnFlow * 3; ++mGrad)
                        {
                            flowOnParticleInEularGrad(i, j, k, mGrad) *= oov;
                        }
                    }
                }
            }

            grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);

            //! Step 3: calc value grad 2
            flowOnParticleInEularGradSecond = 0.0;
            for (int nsurf = 1; nsurf <= nDim; ++nsurf)
            {
                //! For 2D
                //! nsurf = 1, index is (1,0,0)
                //! nsurf = 2, index is (0,1,0)
                //! So the index is the surf index for each cell.
                int il1, jl1, kl1;
                GetNsurfIndex(nsurf, il1, jl1, kl1);

                int ist = 1;
                int ied = ni - 1 + il1;
                int jst = 1;
                int jed = nj - 1 + jl1;
                int kst = 1;
                int ked = nk - 1 + kl1;

                if (TWO_D == nDim) ked = 1;

                int il, jl, kl;

                //! The gradient of u,v,w.
                //! Loop for cell with 1 layer.
                for (int k = kst; k <= ked; ++k)
                {
                    for (int j = jst; j <= jed; ++j)
                    {
                        for (int i = ist; i <= ied; ++i)
                        {
                            //! i,j,k and nsurf is the nsurf-st surface for cell (i,j,k)
                            //! il,jl,kl is the cell which shard the surf whit cell(i,j,k)

                            il = i - il1;
                            jl = j - jl1;
                            kl = k - kl1;

                            //! Loop for variable.
                            for (int m = 0; m < nVarOnFlow; ++m)
                            {
                                //! Note that xfv does not have value on (0,0,1) and (-1,1,1)
                                //! but is has value on (0,1,1)

                                //! considers a linear variation of the field between 
                                //! the centroids of the cells which share the face.

                                int mGrad = m * 3 - 1;

                                RDouble phisx = flowOnParticleInEularGrad(i, j, k, mGrad + 1) + flowOnParticleInEularGrad(il, jl, kl, mGrad + 1);
                                RDouble phisy = flowOnParticleInEularGrad(i, j, k, mGrad + 2) + flowOnParticleInEularGrad(il, jl, kl, mGrad + 2);
                                RDouble phisz = flowOnParticleInEularGrad(i, j, k, mGrad + 3) + flowOnParticleInEularGrad(il, jl, kl, mGrad + 3);

                                RDouble ddxx = xfv(i, j, k, nsurf) * phisx;
                                RDouble ddxy = xfv(i, j, k, nsurf) * phisy;
                                RDouble ddxz = xfv(i, j, k, nsurf) * phisz;

                                RDouble ddyx = yfv(i, j, k, nsurf) * phisx;
                                RDouble ddyy = yfv(i, j, k, nsurf) * phisy;
                                RDouble ddyz = yfv(i, j, k, nsurf) * phisz;

                                RDouble ddzx = zfv(i, j, k, nsurf) * phisx;
                                RDouble ddzy = zfv(i, j, k, nsurf) * phisy;
                                RDouble ddzz = zfv(i, j, k, nsurf) * phisz;

                                //! For the surf in current cell, the surf is inwards
                                //! but in Gauss formula we need the outwards pointing face area vector.
                                //! 3 * 3 for a Hessian matrix.
                                int mGradSecond;
                                mGradSecond = m * 3 * 3 - 1;

                                flowOnParticleInEularGradSecond(i, j, k, ++mGradSecond) -= ddxx;
                                flowOnParticleInEularGradSecond(i, j, k, ++mGradSecond) -= ddxy;
                                flowOnParticleInEularGradSecond(i, j, k, ++mGradSecond) -= ddxz;

                                flowOnParticleInEularGradSecond(i, j, k, ++mGradSecond) -= ddyx;
                                flowOnParticleInEularGradSecond(i, j, k, ++mGradSecond) -= ddyy;
                                flowOnParticleInEularGradSecond(i, j, k, ++mGradSecond) -= ddyz;

                                flowOnParticleInEularGradSecond(i, j, k, ++mGradSecond) -= ddzx;
                                flowOnParticleInEularGradSecond(i, j, k, ++mGradSecond) -= ddzy;
                                flowOnParticleInEularGradSecond(i, j, k, ++mGradSecond) -= ddzz;

                                //! These face area vector is outwards pointing.
                                mGradSecond = m * 3 * 3 - 1;
                                flowOnParticleInEularGradSecond(il, jl, kl, ++mGradSecond) += ddxx;
                                flowOnParticleInEularGradSecond(il, jl, kl, ++mGradSecond) += ddxy;
                                flowOnParticleInEularGradSecond(il, jl, kl, ++mGradSecond) += ddxz;

                                flowOnParticleInEularGradSecond(il, jl, kl, ++mGradSecond) += ddyx;
                                flowOnParticleInEularGradSecond(il, jl, kl, ++mGradSecond) += ddyy;
                                flowOnParticleInEularGradSecond(il, jl, kl, ++mGradSecond) += ddyz;

                                flowOnParticleInEularGradSecond(il, jl, kl, ++mGradSecond) += ddzx;
                                flowOnParticleInEularGradSecond(il, jl, kl, ++mGradSecond) += ddzy;
                                flowOnParticleInEularGradSecond(il, jl, kl, ++mGradSecond) += ddzz;
                            }
                        }
                    }
                }
            }

            for (int k = kst; k <= ked; ++k)
            {
                for (int j = jst; j <= jed; ++j)
                {
                    for (int i = ist; i <= ied; ++i)
                    {
                        //! half is a linear variation.
                        RDouble oov = half / vol(i, j, k);

                        //! Loop for variable.
                        for (int mGradSecond = 0; mGradSecond < nVarOnFlow * 3 * 3; ++mGradSecond)
                        {
                            flowOnParticleInEularGradSecond(i, j, k, mGradSecond) *= oov;
                        }
                    }
                }
            }
        }
    }
}

void ParticlePointSolverParallelDomain::CalculateVariableForParticleOnEular(UnstructGrid *grid, ActionKey *actkey)
{

}

void ParticlePointSolverParallelDomain::ParticleInterpolation(ActionKey *actkey, Grid *grid, OnePointVariable *onePointVariable)
{
    int gridType = grid->Type();

    Param_ParticlePointSolverParallelDomain *parameter = GetControlParameters();

    if (PHSPACE::STRUCTGRID == gridType)
    {
        StructGrid *structgrid = StructGridCast(grid);
        //! Get the interpolation method of particle on struct grid.
        if (parameter->GetParticleInterpolationStruct() == "Gimenez2019JCP")
        {
            RDouble4D &flowOnParticleInEular = *reinterpret_cast <RDouble4D*> (structgrid->GetDataPtr("flowOnParticleInEular"));
            RDouble4D &flowOnParticleInEularGrad = *reinterpret_cast <RDouble4D*> (structgrid->GetDataPtr("flowOnParticleInEularGrad"));
            RDouble4D &flowOnParticleInEularGradSecond = *reinterpret_cast <RDouble4D*> (structgrid->GetDataPtr("flowOnParticleInEularGradSecond"));

            InterpolationGimenez2019JCP(actkey, structgrid, onePointVariable, flowOnParticleInEular, flowOnParticleInEularGrad, flowOnParticleInEularGradSecond);
        }
     
    }
    else if (PHSPACE::UNSTRUCTGRID == gridType)
    {
        UnstructGrid *unstructgrid = UnstructGridCast(grid);
        if (parameter->GetParticleInterpolationUnstruct() == "Gimenez2019JCP")
        {

        }
    }
    else if (PHSPACE::MIXGRID == gridType)
    {
        //! Mixing
        TK_Exit::UnexpectedVarValue(" gridType = ", PHSPACE::MIXGRID);
    }
}

void ParticlePointSolverParallelDomain::InterpolationGimenez2019JCP(ActionKey *actkey, Grid *grid, OnePointVariable *onePointVariable)
{
    int gridType = grid->Type();

    if (PHSPACE::STRUCTGRID == gridType)
    {
        StructGrid *structgrid = StructGridCast(grid);
        RDouble4D &flowOnParticleInEular = *reinterpret_cast <RDouble4D*> (structgrid->GetDataPtr("flowOnParticleInEular"));
        RDouble4D &flowOnParticleInEularGrad = *reinterpret_cast <RDouble4D*> (structgrid->GetDataPtr("flowOnParticleInEularGrad"));
        RDouble4D &flowOnParticleInEularGradSecond = *reinterpret_cast <RDouble4D*> (structgrid->GetDataPtr("flowOnParticleInEularGradSecond"));

        InterpolationGimenez2019JCP(actkey, structgrid, onePointVariable,flowOnParticleInEular, flowOnParticleInEularGrad, flowOnParticleInEularGradSecond);
    }
    else if (PHSPACE::UNSTRUCTGRID == gridType)
    {
        UnstructGrid *unstructgrid = UnstructGridCast(grid);
    }
    else if (PHSPACE::MIXGRID == gridType)
    {
        //! Mixing
        TK_Exit::UnexpectedVarValue(" gridType = ", PHSPACE::MIXGRID);
    }
}

void ParticlePointSolverParallelDomain::InterpolationGimenez2019JCP(ActionKey *actkey, UnstructGrid *grid, OnePointVariable *onePointVariable,
    RDouble **flowOnParticleInEular, RDouble **flowOnParticleInEularGrad, RDouble **flowOnParticleInEularGradSecond)
{

}

void ParticlePointSolverParallelDomain::InterpolationGimenez2019JCP(ActionKey *actkey, StructGrid *grid, OnePointVariable *onePointVariable,
    RDouble4D &flowOnParticleInEular, RDouble4D &flowOnParticleInEularGrad, RDouble4D &flowOnParticleInEularGradSecond)
{
    Param_ParticlePointSolverParallelDomain *parameter = GetControlParameters();

    //! Get the velocity of particle.
    RDouble particleDiameter = onePointVariable->GetOnePointDiameter();
    RDouble particleDensity = onePointVariable->GetOnePointDensity();
    SPDouble &particleVelocity = *(onePointVariable->GetOnePointVelocity());
    SPDouble &particleCoordinate = *(onePointVariable->GetOnePointCoordinate());
    SPDouble &particleAngularVelocity = *(onePointVariable->GetOnePointAngularVelocity());
    SPDouble &particleAcceleration = *(onePointVariable->GetOnePointAcceleration());
    SPDouble &particleTemperature = *(onePointVariable->GetOnePointTemperature());

    SPDouble &flowOnParticle = *(onePointVariable->GetOnePointFlowOnParticle());

    RDouble3D &xcc = *(grid->GetCellCenterX());
    RDouble3D &ycc = *(grid->GetCellCenterY());
    RDouble3D &zcc = *(grid->GetCellCenterZ());

    int cellID = onePointVariable->GetCellID();
    int nLayers = GetNumberOfGhostCellLayers();

    int nI, nJ, nK;
    int iP, jP, kP;
    grid->GetND(nI, nJ, nK);
    nI -= 1;
    nJ -= 1;
    nK -= 1;
    if (TWO_D == GetDim())
    {
        nK = 1;
    }

    GetIndexCellIJKStruct(nI, nJ, nK, iP, jP, kP, nLayers, cellID);

    int nDimVar = 3;
    SPDouble relativeCoord(nDimVar);
    relativeCoord[0] = particleCoordinate[0] - xcc(iP, jP, kP);
    relativeCoord[1] = particleCoordinate[1] - ycc(iP, jP, kP);
    relativeCoord[2] = particleCoordinate[2] - zcc(iP, jP, kP);

    using namespace INDEX_FLOWONPARTICLE;

    //! Init particle as flow.
    if (0 == actkey->subact)
    {
        RDouble dotCoordGrad, dotCoordGrad2;
        int mGrad;

        mGrad = PFU;
        dotCoordGrad = GradDot(relativeCoord, flowOnParticleInEularGrad,iP,jP,kP, mGrad);
        dotCoordGrad2 = GradDoubleDot(relativeCoord, flowOnParticleInEularGradSecond, iP, jP, kP, mGrad);
        particleVelocity[0] = flowOnParticleInEular(iP, jP, kP, mGrad) + dotCoordGrad + 0.5 * dotCoordGrad2;

        mGrad = PFV;
        dotCoordGrad = GradDot(relativeCoord, flowOnParticleInEularGrad, iP, jP, kP, mGrad);
        dotCoordGrad2 = GradDoubleDot(relativeCoord, flowOnParticleInEularGradSecond, iP, jP, kP, mGrad);
        particleVelocity[1] = flowOnParticleInEular(iP, jP, kP, mGrad) + dotCoordGrad + 0.5 * dotCoordGrad2;

        mGrad = PFW;
        dotCoordGrad = GradDot(relativeCoord, flowOnParticleInEularGrad, iP, jP, kP, mGrad);
        dotCoordGrad2 = GradDoubleDot(relativeCoord, flowOnParticleInEularGradSecond, iP, jP, kP, mGrad);
        particleVelocity[2] = flowOnParticleInEular(iP, jP, kP, mGrad) + dotCoordGrad + 0.5 * dotCoordGrad2;

        if (parameter->ifUseParticleTemperature())
        {
            mGrad = PFT;
            dotCoordGrad = GradDot(relativeCoord, flowOnParticleInEularGrad, iP, jP, kP, mGrad);
            dotCoordGrad2 = GradDoubleDot(relativeCoord, flowOnParticleInEularGradSecond, iP, jP, kP, mGrad);
            particleTemperature[0] = flowOnParticleInEular(iP, jP, kP, mGrad) + dotCoordGrad + 0.5 * dotCoordGrad2;
        }
       
    }
    else if (1 == actkey->subact)
    {
        RDouble dotCoordGrad, dotCoordGrad2;
        int mGrad;

        for (int iGrad = 0; iGrad < nVarOnFlow; ++iGrad)
        {
            mGrad = iGrad;
            dotCoordGrad = GradDot(relativeCoord, flowOnParticleInEularGrad, iP, jP, kP, mGrad);
            dotCoordGrad2 = GradDoubleDot(relativeCoord, flowOnParticleInEularGradSecond, iP, jP, kP, mGrad);
            flowOnParticle[mGrad] = flowOnParticleInEular(iP, jP, kP, mGrad) + dotCoordGrad + 0.5 * dotCoordGrad2;
        }
    }

    ostringstream ossLog;
    ossLog << " --- Check particle interpolation --- " << "\n";
    ossLog << " ip jp kp : " << iP << " " << jP << " " << kP << "\n";
    ossLog << " flowOnParticleInEularGrad : " << "\n";

    ossLog << " 0 flow  u v w rho mu t: " << "\n";
    ossLog << "  " << flowOnParticleInEular(iP, jP, kP, PFU) << " ";
    ossLog << "  " << flowOnParticleInEular(iP, jP, kP, PFV) << " ";
    ossLog << "  " << flowOnParticleInEular(iP, jP, kP, PFW) << " ";
    ossLog << "\n";
    ossLog << "  " << flowOnParticleInEular(iP, jP, kP, PFR) << " ";
    ossLog << "  " << flowOnParticleInEular(iP, jP, kP, PFMU) << " ";
    ossLog << "  " << flowOnParticleInEular(iP, jP, kP, PFT) << " ";
    ossLog << "\n";

    ossLog << " 1 grad u v w rho mu t: " << "\n";
    ossLog << "  " << flowOnParticleInEularGrad(iP, jP, kP, PFU * 3 + 0) << " ";
    ossLog << "  " << flowOnParticleInEularGrad(iP, jP, kP, PFV * 3 + 1) << " ";
    ossLog << "  " << flowOnParticleInEularGrad(iP, jP, kP, PFW * 3 + 2) << " ";
    ossLog << "\n";
    ossLog << "  " << flowOnParticleInEularGrad(iP, jP, kP, PFR) << " ";
    ossLog << "  " << flowOnParticleInEularGrad(iP, jP, kP, PFMU) << " ";
    ossLog << "  " << flowOnParticleInEularGrad(iP, jP, kP, PFT) << " ";
    ossLog << "\n";

    ossLog << " 2 grad u v w rho mu t : " << "\n";
    ossLog << "  " << flowOnParticleInEularGradSecond(iP, jP, kP, PFU * 3 * 3 + 0) << " ";
    ossLog << "  " << flowOnParticleInEularGradSecond(iP, jP, kP, PFV * 3 * 3 + 1) << " ";
    ossLog << "  " << flowOnParticleInEularGradSecond(iP, jP, kP, PFW * 3 * 3 + 2) << " ";
    ossLog << "\n";
    ossLog << "  " << flowOnParticleInEularGradSecond(iP, jP, kP, PFR) << " ";
    ossLog << "  " << flowOnParticleInEularGradSecond(iP, jP, kP, PFMU) << " ";
    ossLog << "  " << flowOnParticleInEularGradSecond(iP, jP, kP, PFT) << " ";
    ossLog << "\n";
    if(std::isnan(particleCoordinate[0]))
    {
        cout << "error" << endl;
    }
    ossLog << " particle var : " << "\n";
    ossLog << " par x y z : " << "\n";
    ossLog << "  " << particleCoordinate[0] << " ";
    ossLog << "  " << particleCoordinate[1] << " ";
    ossLog << "  " << particleCoordinate[2] << " ";
    ossLog << "\n";
    ossLog << " par u v w : " << "\n";
    ossLog << "  " << particleVelocity[0] << " ";
    ossLog << "  " << particleVelocity[1] << " ";
    ossLog << "  " << particleVelocity[2] << " ";
    ossLog << "\n";

    ossLog << " interpolation : " << "\n";
    ossLog << " inter u v w rho mu t: " << "\n";
    ossLog << flowOnParticle[PFU] << " " << flowOnParticle[PFV] << " " << flowOnParticle[PFW] << "\n";
    ossLog << flowOnParticle[PFR] << " " << flowOnParticle[PFMU] << " " << flowOnParticle[PFT] << "\n";

    //! Check interpolation by ParticleTrilinearInterpolation.
    RDouble trilinearVar;
}

RDouble ParticlePointSolverParallelDomain::GradDot(SPDouble &relativeCoord, RDouble4D &gradient, int &iP, int &jP, int &kP, int index)
{
    RDouble dotCoordGrad =
        relativeCoord[0] * gradient(iP, jP, kP, index * 3 + 0) +
        relativeCoord[1] * gradient(iP, jP, kP, index * 3 + 1) +
        relativeCoord[2] * gradient(iP, jP, kP, index * 3 + 2);
    return dotCoordGrad;
}

RDouble ParticlePointSolverParallelDomain::GradDoubleDot(SPDouble &relativeCoord, RDouble4D &gradient2ord, int &iP, int &jP, int &kP, int index)
{
    RDouble dotCoordGrad2 =
        relativeCoord[0] * gradient2ord(iP, jP, kP, index * 3 * 3 + 0) * relativeCoord[0] +
        relativeCoord[0] * gradient2ord(iP, jP, kP, index * 3 * 3 + 1) * relativeCoord[1] +
        relativeCoord[0] * gradient2ord(iP, jP, kP, index * 3 * 3 + 2) * relativeCoord[2] +
        relativeCoord[1] * gradient2ord(iP, jP, kP, index * 3 * 3 + 3) * relativeCoord[0] +
        relativeCoord[1] * gradient2ord(iP, jP, kP, index * 3 * 3 + 4) * relativeCoord[1] +
        relativeCoord[1] * gradient2ord(iP, jP, kP, index * 3 * 3 + 5) * relativeCoord[2] +
        relativeCoord[2] * gradient2ord(iP, jP, kP, index * 3 * 3 + 6) * relativeCoord[0] +
        relativeCoord[2] * gradient2ord(iP, jP, kP, index * 3 * 3 + 7) * relativeCoord[1] +
        relativeCoord[2] * gradient2ord(iP, jP, kP, index * 3 * 3 + 8) * relativeCoord[2];
    return dotCoordGrad2;
}

RDouble ParticlePointSolverParallelDomain::LinearInterpolation(RDouble &x0, RDouble &f0, RDouble &x1, RDouble &f1, RDouble &x)
{
    RDouble value = 0.0;
    if (x == x0)
    {
        value = f0;
    }
    else if (x == x1)
    {
        value = f1;
    }
    else if ((x - x0) * (x - x1) > 0)
    {
        ostringstream ossError;
        ossError << " error in linear interpolation " << "\n";
        ossError << " x0 : " << x0 << " x1 : " << x1 << "\n";
        ossError << " but for x : " << x << endl;
        TK_Exit::ExceptionExit(ossError);
    }
    else
    {
        value = f0 * (x1 - x) / (x1 - x0) + f1 * (x - x0) / (x1 - x0);
    }
    return value;
}

void ParticlePointSolverParallelDomain::CheckCornerAndGhost()
{
    Grid *grid = GetGrid();
    if (UNSTRUCTGRID == grid->Type())
    {
        return;
    }

    StructGrid *structGrid = StructGridCast(grid);

    //! Grid variable.
    RDouble3D &structX = *(structGrid->GetStructX());
    RDouble3D &structY = *(structGrid->GetStructY());
    RDouble3D &structZ = *(structGrid->GetStructZ());

    RDouble3D &xcc = *(structGrid->GetCellCenterX());
    RDouble3D &ycc = *(structGrid->GetCellCenterY());
    RDouble3D &zcc = *(structGrid->GetCellCenterZ());

    RDouble4D &xfv = *(structGrid->GetFaceVectorX());
    RDouble4D &yfv = *(structGrid->GetFaceVectorY());
    RDouble4D &zfv = *(structGrid->GetFaceVectorZ());

    RDouble4D &xfn = *(structGrid->GetFaceNormalX());
    RDouble4D &yfn = *(structGrid->GetFaceNormalY());
    RDouble4D &zfn = *(structGrid->GetFaceNormalZ());

    RDouble3D &vol = *structGrid->GetCellVolume();

    //! Get dynamic viscous of Laminar.
    RDouble3D &visl = *reinterpret_cast <RDouble3D*> (structGrid->GetDataPtr("visl"));

    //! Get temperature
    RDouble4D &t = *reinterpret_cast <RDouble4D*> (structGrid->GetDataPtr("t"));

    //! Get the variable of flow.
    RDouble4D &q = *reinterpret_cast <RDouble4D*> (structGrid->GetDataPtr("q"));

    using namespace PHMPI;
    int zoneGloablID = structGrid->GetZoneID();
    int ip, jp, kp;
    ostringstream ossLog;
    ossLog << " --- Check grid corner and ghost --- " << "\n";
    ossLog << " zoneGloablID : " << zoneGloablID << "\n";
    ossLog << " structX xcc xfn" << "\n";
    ossLog << " variable : structX,xcc,xfv,visl,q,vol : " << "\n";
    ossLog << " zone center : \n";
    ip = 2; jp = 1; kp = 1;
    ossLog << ip << " " << jp << " " << kp << " \n";
    ossLog << structX(ip, jp, kp) << " \n";
    ossLog << xcc(ip, jp, kp) << " \n";
    ossLog << xfv(ip, jp, kp, 1) << " \n";
    ossLog << visl(ip, jp, kp) << " \n";
    ossLog << q(ip, jp, kp, 2) << " \n";
    ossLog << vol(ip, jp, kp) << " \n";
    ossLog << "\n";
    ip = 1; jp = 1; kp = 1;
    ossLog << ip << " " << jp << " " << kp << " \n";
    ossLog << structX(ip, jp, kp) << " \n";
    ossLog << xcc(ip, jp, kp) << " \n";
    ossLog << xfv(ip, jp, kp, 1) << " \n";
    ossLog << visl(ip, jp, kp) << " \n";
    ossLog << q(ip, jp, kp, 2) << " \n";
    ossLog << vol(ip, jp, kp) << " \n";
    ossLog << "\n";

    ossLog << " zone ghost 1 layer : \n";
    ip = 0; jp = 1; kp = 1;
    ossLog << ip << " " << jp << " " << kp << " \n";
    ossLog << structX(ip, jp, kp) << " \n";
    ossLog << xcc(ip, jp, kp) << " \n";
    ossLog << xfv(ip, jp, kp, 1) << " \n";
    ossLog << visl(ip, jp, kp) << " \n";
    ossLog << q(ip, jp, kp, 2) << " \n";
    ossLog << vol(ip, jp, kp) << " \n";
    ossLog << "\n";
    ip = 1; jp = 0; kp = 1;
    ossLog << ip << " " << jp << " " << kp << " \n";
    ossLog << structX(ip, jp, kp) << " \n";
    ossLog << xcc(ip, jp, kp) << " \n";
    ossLog << xfv(ip, jp, kp, 1) << " \n";
    ossLog << visl(ip, jp, kp) << " \n";
    ossLog << q(ip, jp, kp, 2) << " \n";
    ossLog << vol(ip, jp, kp) << " \n";
    ossLog << "\n";

    ossLog << " zone ghost 2 layer : \n";
    ip = -1; jp = 1; kp = 1;
    ossLog << ip << " " << jp << " " << kp << " \n";
    ossLog << structX(ip, jp, kp) << " \n";
    ossLog << xcc(ip, jp, kp) << " \n";
    ossLog << xfv(ip, jp, kp, 1) << " \n";
    ossLog << visl(ip, jp, kp) << " \n";
    ossLog << q(ip, jp, kp, 2) << " \n";
    ossLog << vol(ip, jp, kp) << " \n";
    ossLog << "\n";
    ip = 1; jp = -1; kp = 1;
    ossLog << ip << " " << jp << " " << kp << " \n";
    ossLog << structX(ip, jp, kp) << " \n";
    ossLog << xcc(ip, jp, kp) << " \n";
    ossLog << xfv(ip, jp, kp, 1) << " \n";
    ossLog << visl(ip, jp, kp) << " \n";
    ossLog << q(ip, jp, kp, 2) << " \n";
    ossLog << vol(ip, jp, kp) << " \n";
    ossLog << "\n";

    ossLog << " zone corner 1 layer : \n";
    ip = 0; jp = 0; kp = 1;
    ossLog << ip << " " << jp << " " << kp << " \n";
    ossLog << structX(ip, jp, kp) << " \n";
    ossLog << xcc(ip, jp, kp) << " \n";
    ossLog << xfv(ip, jp, kp, 1) << " \n";
    ossLog << visl(ip, jp, kp) << " \n";
    ossLog << q(ip, jp, kp, 2) << " \n";
    ossLog << vol(ip, jp, kp) << " \n";
    ossLog << "\n";

    ossLog << " zone corner 2 layer : \n";
    ip = 0; jp = -1; kp = 1;
    ossLog << ip << " " << jp << " " << kp << " \n";
    ossLog << structX(ip, jp, kp) << " \n";
    ossLog << xcc(ip, jp, kp) << " \n";
    ossLog << xfv(ip, jp, kp, 1) << " \n";
    ossLog << visl(ip, jp, kp) << " \n";
    ossLog << q(ip, jp, kp, 2) << " \n";
    ossLog << vol(ip, jp, kp) << " \n";
    ossLog << "\n";

    ip = -1; jp = 0; kp = 1;
    ossLog << ip << " " << jp << " " << kp << " \n";
    ossLog << structX(ip, jp, kp) << " \n";
    ossLog << xcc(ip, jp, kp) << " \n";
    ossLog << xfv(ip, jp, kp, 1) << " \n";
    ossLog << visl(ip, jp, kp) << " \n";
    ossLog << q(ip, jp, kp, 2) << " \n";
    ossLog << vol(ip, jp, kp) << " \n";
    ossLog << "\n";

    ip = -1; jp = -1; kp = 1;
    ossLog << ip << " " << jp << " " << kp << " \n";
    ossLog << structX(ip, jp, kp) << " \n";
    ossLog << xcc(ip, jp, kp) << " \n";
    ossLog << xfv(ip, jp, kp, 1) << " \n";
    ossLog << visl(ip, jp, kp) << " \n";
    ossLog << q(ip, jp, kp, 2) << " \n";
    ossLog << vol(ip, jp, kp) << " \n";
    ossLog << "\n";

    bool writeEachProcessor = true;
    WriteLogFile(ossLog, writeEachProcessor);
    
}

RDouble ParticlePointSolverParallelDomain::CheckInterpolation(StructGrid *grid, SPDouble &particleCoordinate, int cellID, RDouble4D &q, int mVar)
{
    RDouble interpolationVar = 0.0;

    //! Check particle interpolation by ParticleTrilinearInterpolation.
    //! Compute node var.
    Range INode(0, 1), JNode(0, 1), KNode(0, 1);
    RDouble3D nodeVar(INode, JNode, KNode);

    //! Get the cell index of particle.
    int nI, nJ, nK;
    int iP, jP, kP;
    grid->GetND(nI, nJ, nK);
    nI -= 1;
    nJ -= 1;
    nK -= 1;
    if (TWO_D == GetDim())
    {
        nK = 1;
    }
    GetIndexCellIJKStruct(nI, nJ, nK, iP, jP, kP, GetNumberOfGhostCellLayers(), cellID);

    RDouble3D &xcc = *(grid->GetCellCenterX());
    RDouble3D &ycc = *(grid->GetCellCenterY());
    RDouble3D &zcc = *(grid->GetCellCenterZ());

    RDouble3D &structX = *(grid->GetStructX());
    RDouble3D &structY = *(grid->GetStructY());
    RDouble3D &structZ = *(grid->GetStructZ());

    int indexLoopIter[3];
    indexLoopIter[0] = 1;
    indexLoopIter[1] = 1;
    indexLoopIter[2] = 1;
    if (TWO_D == GetDim())
    {
        indexLoopIter[2] = 0;
    }

    int dI, dJ, dK;
    dI = 1;
    dJ = 1;
    dK = 1;
    if (TWO_D == GetDim())
    {
        dK = 0;
    }

    ostringstream ossLogCC;
    ossLogCC << " --- Check linear interpolation xcc -- " << "\n";

    for (int k = kP; k <= kP + indexLoopIter[2]; ++k)
    {
        for (int j = jP; j <= jP + indexLoopIter[1]; ++j)
        {
            for (int i = iP; i <= iP + indexLoopIter[0]; ++i)
            {
                //! The index of node.
                RDouble f0, f1, f2, f3;

                ossLogCC.clear();
                ossLogCC.str("");
                ossLogCC << " i j k : " << i << " " << j << " " << k << "\n";
                ossLogCC << " xcc : " << "\n";

                ossLogCC << i - dI << " " << j - dJ << " " << k - dK << " " << xcc(i - dI, j - dJ, k - dK) << "\n";
                ossLogCC << i - dI << " " << j - dJ << " " << k << " " << xcc(i - dI, j - dJ, k) << "\n";

                ossLogCC << i - dI << " " << j << " " << k - dK << " " << xcc(i - dI, j, k - dK) << "\n";
                ossLogCC << i - dI << " " << j << " " << k << " " << xcc(i - dI, j, k) << "\n";

                ossLogCC << i << " " << j - dJ << " " << k - dK << " " << xcc(i, j - dJ, k - dK) << "\n";
                ossLogCC << i << " " << j - dJ << " " << k << " " << xcc(i, j - dJ, k) << "\n";

                ossLogCC << i << " " << j << " " << k - dK << " " << xcc(i, j, k - dK) << "\n";
                ossLogCC << i << " " << j << " " << k << " " << xcc(i, j, k) << "\n";

                f0 = LinearInterpolation(xcc(i - dI, j - dJ, k - dK), q(i - dI, j - dJ, k - dK, mVar), xcc(i, j - dJ, k - dK), q(i, j - dJ, k - dK, mVar), structX(i, j, k));
                f1 = LinearInterpolation(xcc(i - dI, j, k - dK), q(i - dI, j, k - dK, mVar), xcc(i, j, k - dK), q(i, j, k - dK, mVar), structX(i, j, k));
                f2 = LinearInterpolation(xcc(i - dI, j - dJ, k), q(i - dI, j - dJ, k, mVar), xcc(i, j - dJ, k), q(i, j - dJ, k, mVar), structX(i, j, k));
                f3 = LinearInterpolation(xcc(i - dI, j, k), q(i - dI, j, k, mVar), xcc(i, j, k), q(i, j, k, mVar), structX(i, j, k));

                //! For y direction.
                f0 = LinearInterpolation(ycc(i, j - dJ, k - dK), f0, ycc(i, j, k - dK), f1, structY(i, j, k));
                f1 = LinearInterpolation(ycc(i, j - dJ, k), f2, ycc(i, j, k), f3, structY(i, j, k));

                //! For k direction.
                if (TWO_D == GetDim())
                {
                    nodeVar(i - iP, j - jP, 0) = f0;
                }
                else
                {
                    f0 = LinearInterpolation(zcc(i, j, k - dK), f0, zcc(i, j, k), f1, structZ(i, j, k));
                    nodeVar(i - iP, j - jP, k - kP) = f0;
                }
            }
        }
    }
    ossLogCC.clear();
    ossLogCC.str("");

    //! The index of node.
    RDouble f0, f1, f2, f3;
    int i, j, k;
    i = iP;
    j = jP;
    k = kP;

    f0 = LinearInterpolation(structX(i, j, k), nodeVar(0, 0, 0), structX(i + dI, j, k), nodeVar(0 + dI, 0, 0), particleCoordinate[0]);
    f1 = LinearInterpolation(structX(i, j + dJ, k), nodeVar(0, 0 + dJ, 0), structX(i + dI, j + dJ, k), nodeVar(0 + dI, 0 + dJ, 0), particleCoordinate[0]);
    f2 = LinearInterpolation(structX(i, j, k + dK), nodeVar(0, 0, 0 + dK), structX(i + dI, j, k + dK), nodeVar(0 + dI, 0, 0 + dK), particleCoordinate[0]);
    f3 = LinearInterpolation(structX(i, j + dJ, k + dK), nodeVar(0, 0 + dJ, 0 + dK), structX(i + dI, j + dJ, k + dK), nodeVar(0 + dI, 0 + dJ, 0 + dK), particleCoordinate[0]);

    //! For y direction.
    f0 = LinearInterpolation(structY(i, j, k), f0, structY(i, j + dJ, k), f1, particleCoordinate[1]);
    f1 = LinearInterpolation(structY(i, j, k + dK), f2, structY(i, j + dJ, k + dK), f3, particleCoordinate[1]);

    //! For k direction.
    if (TWO_D == GetDim())
    {
        interpolationVar = f0;
    }
    else
    {
        f0 = LinearInterpolation(structZ(i, j, k), f0, structZ(i, j, k + dK), f1, particleCoordinate[2]);
        interpolationVar = f0;
    }

    return interpolationVar;
}

void ParticlePointSolverParallelDomain::CommunicationFlowOnParticleInEular(ActionKey *actkey, Grid *grid)
{
    int gridType = grid->Type();
    Param_ParticlePointSolverParallelDomain *parameter = GetControlParameters();

    if (PHSPACE::STRUCTGRID == gridType)
    {
        StructGrid *structgrid = StructGridCast(grid);
        if (parameter->GetParticleInterpolationStruct() == "Gimenez2019JCP")
        {
            if (0 == actkey->kind)
            {
                //! Only corner for flowOnParticleInEular.
                actkey->format = 1;
                //! MPI corner.
                actkey->subact = 1;
                actkey->filename = "flowOnParticleInEular";
                CommunicateSpecificVariable(actkey);
            }
            else if (1 == actkey->kind)
            {
                actkey->format = 1;
                //! MPI neighbor.
                actkey->subact = 0;
                actkey->filename = "flowOnParticleInEular";
                CommunicateSpecificVariable(actkey);
                //! MPI corner.
                actkey->subact = 1;
                actkey->filename = "flowOnParticleInEular";
                CommunicateSpecificVariable(actkey);
            }
            else if(2 == actkey->kind)
            {
                actkey->format = 1;
                //! MPI neighbor.
                actkey->subact = 0;
                actkey->filename = "flowOnParticleInEularGrad";
                CommunicateSpecificVariable(actkey);
                actkey->filename = "flowOnParticleInEularGradSecond";
                CommunicateSpecificVariable(actkey);
                //! MPI corner.
                actkey->subact = 1;
                actkey->filename = "flowOnParticleInEularGrad";
                CommunicateSpecificVariable(actkey);
                actkey->filename = "flowOnParticleInEularGradSecond";
                CommunicateSpecificVariable(actkey);
            }
        } 
    }
    else if (PHSPACE::UNSTRUCTGRID == gridType)
    {

    }
    else if (PHSPACE::MIXGRID == gridType)
    {
        //! Mixing
        TK_Exit::UnexpectedVarValue(" gridType = ", PHSPACE::MIXGRID);
    }
}

void ParticlePointSolverParallelDomain::CommunicateSpecificVariable(ActionKey *actkey)
{
    int zoneIDofThisSolver = GetZoneLocalID();
    if (zoneIDofThisSolver != 0)
    {
        //! Very Important: only the first zone on a processor need to implement the task.
        //! other wise, the tasks would be implemented several times.
        return;
    }

    //! Same as Controller::CommunicateSpecificArray().
    using namespace PHMPI;
    //! The number of zones on all processor global.
    int nZones = GetNumberofGlobalZones();
    //! Get current processor ID.
    int currentProcessorID = GetCurrentProcessorID();
    
    vector <PH_Request > requestContainer;

    vector <vector <DataContainer*> > receivedDataBuffer;
    receivedDataBuffer.resize(nZones);

    vector <vector <DataContainer*> > sendDataBuffer;
    sendDataBuffer.resize(nZones);

    //! Step 0: Compressing data firstly.
    for (int iZone = 0; iZone < nZones; ++iZone)
    {
        //! Judge the zones on current processor.
        int sendProcessor = GetZoneProcessorID(iZone);

        ZoneNeighbor *zoneNeighbor = zoneConnectivity->GetZoneNeighbor(iZone);
        int nNeighbor = zoneNeighbor->GetNumberOfNeighbor();

        ZoneCorner *zoneCorner = zoneConnectivity->GetZoneCorner(iZone);
        int nCorner = zoneCorner->GetNumOfCorner();

        int numSendBuffer = 0;
        if (0 == actkey->subact)
        {
            //! Neighbor variable.
            numSendBuffer = nNeighbor;
        }
        else if (1 == actkey->subact)
        {
            //! Corner variable.
            numSendBuffer = nCorner;
        }
        else
        {
            ostringstream ossError;
            ossError << " Error on actkey->subact at CompressSpecificArray :";
            ossError << actkey->subact << endl;
            TK_Exit::ExceptionExit(ossError);
        }

        //! The zone on current processor.
        if (currentProcessorID == sendProcessor)
        {
            sendDataBuffer[iZone].resize(numSendBuffer);
        }

        for (int iSendBuffer = 0; iSendBuffer < numSendBuffer; ++iSendBuffer)
        {
            int recvZone = 0;
            if (0 == actkey->subact)
            {
                recvZone = zoneNeighbor->GetZoneIndexOfNeighbor(iSendBuffer);
            }
            else if (1 == actkey->subact)
            {
                recvZone = zoneCorner->GetZoneIndexOfCorner(iSendBuffer);
            }
            actkey->ipos = recvZone;

            int receiveProcessor = PHMPI::GetZoneProcessorID(recvZone);

            if (currentProcessorID != sendProcessor && currentProcessorID != receiveProcessor)
            {
                continue;
            }

            //! Allocate and compress the buffers for sending.
            DataContainer *sendBuffer = 0;
            if (currentProcessorID == sendProcessor)
            {
                sendBuffer = new DataContainer;
                sendDataBuffer[iZone][iSendBuffer] = sendBuffer;

                ParticlePointSolverParallelDomain *solver = static_cast <ParticlePointSolverParallelDomain*> (GlobalSolvers::GetSolver(iZone, actkey->solverID));
                //! Compress the send information into the actkey.
                solver->CompressSpecificArray(sendBuffer, actkey);
            }

        } 
    }

    //! Step 1: Communication.
    for (int iZone = 0; iZone < nZones; ++iZone)
    {
        ZoneNeighbor *zoneNeighbor = zoneConnectivity->GetZoneNeighbor(iZone);
        int nNeighbor = zoneNeighbor->GetNumberOfNeighbor();

        ZoneCorner *zoneCorner = zoneConnectivity->GetZoneCorner(iZone);
        int nCorner = zoneCorner->GetNumOfCorner();

        int numSendBuffer = 0;
        if (0 == actkey->subact)
        {
            //! Neighbor variable.
            numSendBuffer = nNeighbor;
        }
        else if (1 == actkey->subact)
        {
            //! Corner variable.
            numSendBuffer = nCorner;
        }

        //! Communicating.
        for (int iSendBuffer = 0; iSendBuffer < numSendBuffer; ++iSendBuffer)
        {
            int recvZone = 0;
            if (0 == actkey->subact)
            {
                recvZone = zoneNeighbor->GetZoneIndexOfNeighbor(iSendBuffer);
            }
            else if (1 == actkey->subact)
            {
                recvZone = zoneCorner->GetZoneIndexOfCorner(iSendBuffer);
            }
            actkey->ipos = recvZone;

            int sendProcessor = PHMPI::GetZoneProcessorID(iZone);
            int receiveProcessor = PHMPI::GetZoneProcessorID(recvZone);

            if (currentProcessorID != sendProcessor && currentProcessorID != receiveProcessor)
            {
                continue;
            }

            //! Allocate the buffers for sending and receiving.
            DataContainer *receivedBuffer = 0;
            DataContainer *sendBuffer = 0;
            if (currentProcessorID == receiveProcessor)
            {
                receivedBuffer = new DataContainer;
                receivedDataBuffer[iZone].push_back(receivedBuffer);
            }
            if (currentProcessorID == sendProcessor)
            {
                sendBuffer = sendDataBuffer[iZone][iSendBuffer];
            }

            int tag = iZone;

            //! Communication.
            if (currentProcessorID == sendProcessor)
            {
                if (sendProcessor != receiveProcessor)
                {
                    streamsize nlen = sendBuffer->Size();

                    //! Send the data to neighbors.
                    send(sendBuffer, receiveProcessor, nlen, requestContainer, tag);
                }
            }

            if (currentProcessorID == receiveProcessor)
            {
                //! Get the neighbor(Corner) order: 
                //! which order that 'iZone' on neighbors(Corner) of neighbor/Corner Zone(recvZone).
                int recvOrer = -1;
                int numberOfRecvTemp = 0;
                if (0 == actkey->subact)
                {
                    //! The neighbor of recvZone.
                    ZoneNeighbor *zoneNeighborTemp = zoneConnectivity->GetZoneNeighbor(recvZone);
                    //! Number of neighbor of recvZone.
                    numberOfRecvTemp = zoneNeighborTemp->GetNumberOfNeighbor();

                    //! Loop for the neighbor of neighbor(recvZone).
                    for (int neighborID = 0; neighborID < numberOfRecvTemp; ++neighborID)
                    {
                        int neighborZoneTemp = zoneNeighborTemp->GetZoneIndexOfNeighbor(neighborID);
                        if (neighborZoneTemp == iZone)
                        {
                            recvOrer = neighborID;
                            break;
                        }
                    }
                }
                else if(1 == actkey->subact)
                {
                    //! The corner of recvZone.
                    ZoneCorner *zoneCornerTemp = zoneConnectivity->GetZoneCorner(recvZone);

                    //! Number of neighbor of recvZone.
                    numberOfRecvTemp = zoneCornerTemp->GetNumOfCorner();

                    for (map<int, ZoneCornerIndex >::iterator iterCorner = zoneCornerTemp->GetCornerIterBegin();
                        iterCorner!=zoneCornerTemp->GetCornerIterEnd();++iterCorner)
                    {
                        int cornerZoneTemp = iterCorner->first;
                        if (cornerZoneTemp == iZone)
                        {
                            recvOrer = cornerZoneTemp;
                            break;
                        }
                    }
                }

                ASSERT(recvOrer != -1);

                if (sendProcessor != receiveProcessor)
                {
                    //! The data length of the received data is same to the length that send to the neighbor.
                    CharVecSizeType nlen = sendDataBuffer[recvZone][recvOrer]->Size();

                    //! Receive data from neighbors.
                    receive(receivedBuffer, sendProcessor, nlen, requestContainer, tag);
                }
                else
                {
                    SWAP(receivedDataBuffer[iZone].back(), sendDataBuffer[iZone][iSendBuffer]);
                }
            }
        }
    }

    //! Step2: MPI waiting for Non-blocking communications.
    if (PHMPI::GetNumberOfProcessor() > 1)
    {
        PH_Wait(static_cast<int>(requestContainer.size()), &(requestContainer[0]));
    }

    //! Step3: Translate the data container.
    for (int iZone = 0; iZone < nZones; ++iZone)
    {
        ZoneNeighbor *zoneNeighbor = zoneConnectivity->GetZoneNeighbor(iZone);
        int nNeighbor = zoneNeighbor->GetNumberOfNeighbor();

        ZoneCorner *zoneCorner = zoneConnectivity->GetZoneCorner(iZone);
        int nCorner = zoneCorner->GetNumOfCorner();

        int numSendBuffer = 0;
        if (0 == actkey->subact)
        {
            //! Neighbor variable.
            numSendBuffer = nNeighbor;
        }
        else if (1 == actkey->subact)
        {
            //! Corner variable.
            numSendBuffer = nCorner;
        }

        int count = 0;
        for (int iSendBuffer = 0; iSendBuffer < numSendBuffer; ++iSendBuffer)
        {
            int recvZone = 0;
            if (0 == actkey->subact)
            {
                recvZone = zoneNeighbor->GetZoneIndexOfNeighbor(iSendBuffer);
            }
            else if (1 == actkey->subact)
            {
                recvZone = zoneCorner->GetZoneIndexOfCorner(iSendBuffer);
            }
            actkey->ipos = recvZone;

            int sendProcessor = PHMPI::GetZoneProcessorID(iZone);
            int receiveProcessor = PHMPI::GetZoneProcessorID(recvZone);

            if (currentProcessorID == receiveProcessor)
            {
                DataContainer *receiveData = receivedDataBuffer[iZone][count];

                //! Because get into here when currentProcessor==receiveProcessor,
                //! so, use zoneGridSolver of neighborZone, this zone is on the current processor!

                ParticlePointSolverParallelDomain *solver = static_cast <ParticlePointSolverParallelDomain*> (GlobalSolvers::GetSolver(iZone, actkey->solverID));
                //! Compress the send information into the actkey.
                
                actkey->ipos = iZone;
                solver->DecompressSpecificArray(receiveData, actkey);
                ++count;
            }
        }

    }

    //! Step4: Free the buffers.
    for (int iDim = 0; iDim < receivedDataBuffer.size(); ++iDim)
    {
        for (int jDim = 0; jDim < receivedDataBuffer[iDim].size(); ++jDim)
        {
            delete receivedDataBuffer[iDim][jDim];
        }
    }
    receivedDataBuffer.clear();
    for (int iDim = 0; iDim < sendDataBuffer.size(); ++iDim)
    {
        for (int jDim = 0; jDim < sendDataBuffer[iDim].size(); ++jDim)
        {
            delete sendDataBuffer[iDim][jDim];
        }
    }
    sendDataBuffer.clear();
    requestContainer.clear();

}

void ParticlePointSolverParallelDomain::CompressSpecificArray(DataContainer *&dataContainer, ActionKey *actkey)
{
    Grid *grid = GetGrid();
    int gridType = grid->Type();

    if (PHSPACE::STRUCTGRID == gridType)
    {
        StructGrid *structgrid = StructGridCast(grid);
        CompressSpecificArray(structgrid, dataContainer, actkey);
    }
    else if (PHSPACE::UNSTRUCTGRID == gridType)
    {
        UnstructGrid *unstructgrid = UnstructGridCast(grid);
        CompressSpecificArray(unstructgrid, dataContainer, actkey);
    }
    else if (PHSPACE::MIXGRID == gridType)
    {
        //! Mixing
        TK_Exit::UnexpectedVarValue(" gridType = ", PHSPACE::MIXGRID);
    }
}

void ParticlePointSolverParallelDomain::CompressSpecificArray(UnstructGrid *grid, DataContainer *&dataContainer, ActionKey *actkey)
{

}

void ParticlePointSolverParallelDomain::CompressSpecificArray(StructGrid *grid, DataContainer *&dataContainer, ActionKey *actkey)
{
    dataContainer->MoveToBegin();

    if (0 == actkey->format)
    {
        RDouble3D &fieldSend = *reinterpret_cast <RDouble3D*> (grid->GetDataPtr(actkey->filename));
        if (0 == actkey->subact)
        {
            InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
            int iNeighborZone = interfaceInformation->FindIthNeighbor(actkey->ipos);
            int interfaceNumberBetweenTwoNeighboringZone = interfaceInformation->GetNIFaceOfNeighbor(iNeighborZone);
            int *interfaceIndexContainerForSend = interfaceInformation->GetFaceIndexForSend(iNeighborZone);
            for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; --iGhostLayer)
            {
                for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++iLocalFace)
                {
                    int iFace = interfaceIndexContainerForSend[iLocalFace];
                    int is, js, ks;
                    grid->GetSourceIndexIJK(iFace, iGhostLayer + 1, is, js, ks);
                    grid->RemapMultigridIJK(actkey->level, is, js, ks);
                    PHWrite(dataContainer, fieldSend(is, js, ks));
                }
            }
        }
        else if(1 == actkey->subact)
        {
            ZoneCorner *zoneCorner = zoneConnectivity->GetZoneCorner(grid->GetZoneID());
            map<int, ZoneCornerIndex> cornerInfo = zoneCorner->GetCornerInfo();

            int *indexCornerGrid = 0;
            map<int, ZoneCornerIndex>::iterator iterZoneCorner = cornerInfo.find(actkey->ipos);
            if (iterZoneCorner != cornerInfo.end())
            {
                ZoneCornerIndex zoneCornerIndex = iterZoneCorner->second;

                if (iterZoneCorner->second.GetCornerZoneIndex() != actkey->ipos)
                {
                    ostringstream ossError;
                    ossError << " Error on iZoneToRecv " << "\n";
                    ossError << "iZoneToRecv in Action : " << actkey->ipos << "\n";
                    ossError << "iZoneToRecv in cornerInfo : " << iterZoneCorner->second.GetCornerZoneIndex() << "\n";
                    ossError << endl;
                    TK_Exit::ExceptionExit(ossError);
                }
                indexCornerGrid = zoneCornerIndex.GetCornerGridIndex(PHSPACE::STRUCTGRID);
            }
            else
            {
                ostringstream ossError;
                ossError << " Error on cornerInfo " << "\n";
                ossError << " There is no iZoneToRecv  " << actkey->ipos << "\n";
                ossError << " the Size of Corner : " << cornerInfo.size() << "\n";
                for (map<int, ZoneCornerIndex>::iterator iter = cornerInfo.begin(); iter != cornerInfo.end(); ++iter)
                {
                    ossError << " corner Zone index : " << iter->first << "\n";
                }
                ossError << endl;
                TK_Exit::ExceptionExit(ossError);
            }

            int iCorner, jCorner, kCorner;
            int nLayers = GetNumberOfGhostCellLayers();
            for (int iLayer = 1; iLayer <= nLayers; ++iLayer)
            {
                set<int*> indexCorner;
                grid->GetCornerSourceIndexIJK(indexCornerGrid, indexCorner, iLayer);
                int nCornerCell = indexCorner.size();
                for (set<int*>::iterator iterCorner = indexCorner.begin(); iterCorner != indexCorner.end(); ++iterCorner)
                {
                    iCorner = (*iterCorner)[0];
                    jCorner = (*iterCorner)[1];
                    kCorner = (*iterCorner)[2];
                    PHWrite(dataContainer, fieldSend(iCorner, jCorner, kCorner));
                }

                indexCorner.clear();
            }
        }
    }
    else if(1 == actkey->format)
    {
        RDouble4D &fieldSend = *reinterpret_cast <RDouble4D*> (grid->GetDataPtr(actkey->filename));
        if (0 == actkey->subact)
        {
            InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
            int iNeighborZone = interfaceInformation->FindIthNeighbor(actkey->ipos);
            int interfaceNumberBetweenTwoNeighboringZone = interfaceInformation->GetNIFaceOfNeighbor(iNeighborZone);
            int *interfaceIndexContainerForSend = interfaceInformation->GetFaceIndexForSend(iNeighborZone);
            for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; --iGhostLayer)
            {
                for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++iLocalFace)
                {
                    int iFace = interfaceIndexContainerForSend[iLocalFace];
                    int is, js, ks;
                    grid->GetSourceIndexIJK(iFace, iGhostLayer + 1, is, js, ks);
                    grid->RemapMultigridIJK(actkey->level, is, js, ks);
                    for (int m = 0; m < actkey->filepos; ++m)
                    {
                        PHWrite(dataContainer, fieldSend(is, js, ks, m));
                    }
                }
            }
        }
        else if(1 == actkey->subact)
        {
            ZoneCorner *zoneCorner = zoneConnectivity->GetZoneCorner(grid->GetZoneID());
            map<int, ZoneCornerIndex> cornerInfo = zoneCorner->GetCornerInfo();

            int *indexCornerGrid = 0;
            map<int, ZoneCornerIndex>::iterator iterZoneCorner = cornerInfo.find(actkey->ipos);
            if (iterZoneCorner != cornerInfo.end())
            {
                ZoneCornerIndex zoneCornerIndex = iterZoneCorner->second;

                if (iterZoneCorner->second.GetCornerZoneIndex() != actkey->ipos)
                {
                    ostringstream ossError;
                    ossError << " Error on iZoneToRecv " << "\n";
                    ossError << "iZoneToRecv in Action : " << actkey->ipos << "\n";
                    ossError << "iZoneToRecv in cornerInfo : " << iterZoneCorner->second.GetCornerZoneIndex() << "\n";
                    ossError << endl;
                    TK_Exit::ExceptionExit(ossError);
                }
                indexCornerGrid = zoneCornerIndex.GetCornerGridIndex(PHSPACE::STRUCTGRID);
            }
            else
            {
                ostringstream ossError;
                ossError << " Error on cornerInfo " << "\n";
                ossError << " There is no iZoneToRecv  " << actkey->ipos << "\n";
                ossError << " the Size of Corner : " << cornerInfo.size() << "\n";
                for (map<int, ZoneCornerIndex>::iterator iter = cornerInfo.begin(); iter != cornerInfo.end(); ++iter)
                {
                    ossError << " corner Zone index : " << iter->first << "\n";
                }
                ossError << endl;
                TK_Exit::ExceptionExit(ossError);
            }

            int iCorner, jCorner, kCorner;
            int nLayers = GetNumberOfGhostCellLayers();
            for (int iLayer = 1; iLayer <= nLayers; ++iLayer)
            {
                set<int*> indexCorner;
                grid->GetCornerSourceIndexIJK(indexCornerGrid, indexCorner, iLayer);
                int nCornerCell = indexCorner.size();
                for (set<int*>::iterator iterCorner = indexCorner.begin(); iterCorner != indexCorner.end(); ++iterCorner)
                {
                    iCorner = (*iterCorner)[0];
                    jCorner = (*iterCorner)[1];
                    kCorner = (*iterCorner)[2];

                    for (int m = 0; m < actkey->filepos; ++m)
                    {
                        PHWrite(dataContainer, fieldSend(iCorner, jCorner, kCorner, m));
                    }
                }

                indexCorner.clear();
            }
        }
    }
    else if (1 == actkey->format)
    {

    }
    else
    {
        ostringstream ossError;
        ossError << " Error on actkey->format at CompressSpecificArray :";
        ossError << actkey->format << endl;
        TK_Exit::ExceptionExit(ossError);
    }
}

void ParticlePointSolverParallelDomain::DecompressSpecificArray(DataContainer *&dataContainer, ActionKey *actkey)
{
    Grid *grid = GetGrid();
    int gridType = grid->Type();

    if (PHSPACE::STRUCTGRID == gridType)
    {
        StructGrid *structgrid = StructGridCast(grid);
        DecompressSpecificArray(structgrid, dataContainer, actkey);
    }
    else if (PHSPACE::UNSTRUCTGRID == gridType)
    {
        UnstructGrid *unstructgrid = UnstructGridCast(grid);
        DecompressSpecificArray(unstructgrid, dataContainer, actkey);
    }
    else if (PHSPACE::MIXGRID == gridType)
    {
        //! Mixing
        TK_Exit::UnexpectedVarValue(" gridType = ", PHSPACE::MIXGRID);
    }
}

void ParticlePointSolverParallelDomain::DecompressSpecificArray(UnstructGrid *grid, DataContainer *&dataContainer, ActionKey *actkey)
{

}

void ParticlePointSolverParallelDomain::DecompressSpecificArray(StructGrid *grid, DataContainer *&dataContainer, ActionKey *actkey)
{
    dataContainer->MoveToBegin();

    if (1 == actkey->format)
    {
        RDouble3D &fieldRecv = *reinterpret_cast <RDouble3D*> (grid->GetDataPtr(actkey->filename));
        if (0 == actkey->subact)
        {
            InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
            int iNeighbor = interfaceInformation->FindIthNeighbor(actkey->ipos);
            int interfaceNumberBetweenTwoNeighboringZone = interfaceInformation->GetNIFaceOfNeighbor(iNeighbor);
            int *interfaceIndexContainerForReceive = interfaceInformation->GetFaceIndexForRecv(iNeighbor);
            for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; --iGhostLayer)
            {
                for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++iLocalFace)
                {
                    int iFace = interfaceIndexContainerForReceive[iLocalFace];
                    int is, js, ks;
                    grid->GetSourceIndexIJK(iFace, iGhostLayer + 1, is, js, ks);
                    grid->RemapMultigridIJK(actkey->level, is, js, ks);
                    PHRead(dataContainer, fieldRecv(is, js, ks));
                }
            }
        }
        else if (1 == actkey->subact)
        {
            ZoneCorner *zoneCorner = zoneConnectivity->GetZoneCorner(grid->GetZoneID());
            map<int, ZoneCornerIndex> cornerInfo = zoneCorner->GetCornerInfo();

            int *indexCornerGrid = 0;
            map<int, ZoneCornerIndex>::iterator iterZoneCorner = cornerInfo.find(actkey->ipos);
            if (iterZoneCorner != cornerInfo.end())
            {
                ZoneCornerIndex zoneCornerIndex = iterZoneCorner->second;

                if (iterZoneCorner->second.GetCornerZoneIndex() != actkey->ipos)
                {
                    ostringstream ossError;
                    ossError << " Error on iZoneToRecv " << "\n";
                    ossError << "iZoneToSend in Action : " << actkey->ipos << "\n";
                    ossError << "iZoneToSend in cornerInfo : " << iterZoneCorner->second.GetCornerZoneIndex() << "\n";
                    ossError << endl;
                    TK_Exit::ExceptionExit(ossError);
                }
                indexCornerGrid = zoneCornerIndex.GetCornerGridIndex(PHSPACE::STRUCTGRID);
            }
            else
            {
                ostringstream ossError;
                ossError << " Error on cornerInfo " << "\n";
                ossError << " There is no iZoneToRecv  " << actkey->ipos << "\n";
                ossError << " the Size of Corner : " << cornerInfo.size() << "\n";
                for (map<int, ZoneCornerIndex>::iterator iter = cornerInfo.begin(); iter != cornerInfo.end(); ++iter)
                {
                    ossError << " corner Zone index : " << iter->first << "\n";
                }
                ossError << endl;
                TK_Exit::ExceptionExit(ossError);
            }

            int iCorner, jCorner, kCorner;
            int nLayers = GetNumberOfGhostCellLayers();
            for (int iLayer = 1; iLayer <= nLayers; ++iLayer)
            {
                set<int*> indexCorner;
                grid->GetCornerTargetIndexIJK(indexCornerGrid, indexCorner, iLayer);
                int nCornerCell = indexCorner.size();
                for (set<int*>::iterator iterCorner = indexCorner.begin(); iterCorner != indexCorner.end(); ++iterCorner)
                {
                    iCorner = (*iterCorner)[0];
                    jCorner = (*iterCorner)[1];
                    kCorner = (*iterCorner)[2];
                    PHRead(dataContainer, fieldRecv(iCorner, jCorner, kCorner));
                }

                indexCorner.clear();
            }
        }
    }
    else if (1 == actkey->format)
    {
        RDouble4D &fieldRecv = *reinterpret_cast <RDouble4D*> (grid->GetDataPtr(actkey->filename));
        if (0 == actkey->subact)
        {
            InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
            int iNeighbor = interfaceInformation->FindIthNeighbor(actkey->ipos);
            int interfaceNumberBetweenTwoNeighboringZone = interfaceInformation->GetNIFaceOfNeighbor(iNeighbor);
            int *interfaceIndexContainerForReceive = interfaceInformation->GetFaceIndexForRecv(iNeighbor);
            for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; --iGhostLayer)
            {
                for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++iLocalFace)
                {
                    int iFace = interfaceIndexContainerForReceive[iLocalFace];
                    int is, js, ks;
                    grid->GetSourceIndexIJK(iFace, iGhostLayer + 1, is, js, ks);
                    grid->RemapMultigridIJK(actkey->level, is, js, ks);
                    for (int m = 0; m < actkey->filepos; ++m)
                    {
                        PHRead(dataContainer, fieldRecv(is, js, ks, m));
                    }
                }
            }
        }
        else if (1 == actkey->subact)
        {
            ZoneCorner *zoneCorner = zoneConnectivity->GetZoneCorner(grid->GetZoneID());
            map<int, ZoneCornerIndex> cornerInfo = zoneCorner->GetCornerInfo();

            int *indexCornerGrid = 0;
            map<int, ZoneCornerIndex>::iterator iterZoneCorner = cornerInfo.find(actkey->ipos);
            if (iterZoneCorner != cornerInfo.end())
            {
                ZoneCornerIndex zoneCornerIndex = iterZoneCorner->second;

                if (iterZoneCorner->second.GetCornerZoneIndex() != actkey->ipos)
                {
                    ostringstream ossError;
                    ossError << " Error on iZoneToRecv " << "\n";
                    ossError << "iZoneToSend in Action : " << actkey->ipos << "\n";
                    ossError << "iZoneToSend in cornerInfo : " << iterZoneCorner->second.GetCornerZoneIndex() << "\n";
                    ossError << endl;
                    TK_Exit::ExceptionExit(ossError);
                }
                indexCornerGrid = zoneCornerIndex.GetCornerGridIndex(PHSPACE::STRUCTGRID);
            }
            else
            {
                ostringstream ossError;
                ossError << " Error on cornerInfo " << "\n";
                ossError << " There is no iZoneToSend " << actkey->ipos << "\n";
                ossError << " the Size of Corner : " << cornerInfo.size() << "\n";
                for (map<int, ZoneCornerIndex>::iterator iter = cornerInfo.begin(); iter != cornerInfo.end(); ++iter)
                {
                    ossError << " corner Zone index : " << iter->first << "\n";
                }
                ossError << endl;
                TK_Exit::ExceptionExit(ossError);
            }

            int iCorner, jCorner, kCorner;
            int nLayers = GetNumberOfGhostCellLayers();
            for (int iLayer = 1; iLayer <= nLayers; ++iLayer)
            {
                set<int*> indexCorner;
                grid->GetCornerTargetIndexIJK(indexCornerGrid, indexCorner, iLayer);
                int nCornerCell = indexCorner.size();
                for (set<int*>::iterator iterCorner = indexCorner.begin(); iterCorner != indexCorner.end(); ++iterCorner)
                {
                    iCorner = (*iterCorner)[0];
                    jCorner = (*iterCorner)[1];
                    kCorner = (*iterCorner)[2];

                    for (int m = 0; m < actkey->filepos; ++m)
                    {
                        PHRead(dataContainer, fieldRecv(iCorner, jCorner, kCorner, m));
                    }
                }

                indexCorner.clear();
            }
        }
    }
    else if (1 == actkey->format)
    {

    }
    else
    {
        ostringstream ossError;
        ossError << " Error on actkey->format at CompressSpecificArray :";
        ossError << actkey->format << endl;
        TK_Exit::ExceptionExit(ossError);
    }
}

void ParticlePointSolverParallelDomain::CalcCFLNumber()
{
    RDouble particleCFLThreshold = 1.0;
    RDouble particleCFLMax = SMALL;
    RDouble particleCFL = 0.0;
    RDouble timeCFL = 0.0;
    RDouble timeCFLMin = LARGE;

    Grid *grid = GetGrid();

    using namespace PHMPI;
    //! The global zone id on current zone on current processer.
    int zoneIDGlobal = grid->GetZoneID();
    //! The number of zones on global all processer.
    int nZonesGlobal = GetNumberofGlobalZones();

    //! The local zone id on current zone on current processer.
    int zoneIDLocal = grid->GetZoneLocalID();
    //! The number of zones on on current processer.
    int nZonesLocal = GetNumberofLocalZones();

    if (zoneIDLocal != 0)
    {
        return;
    }

    //! Get current processor ID on current zones..
    int currentProcessorID = GetCurrentProcessorID();

    //! The ProcessorID of current zone.
    int currentZoneProcessorID = GetZoneProcessorID(zoneIDGlobal);

    Param_ParticlePointSolverParallelDomain *parameter = GetControlParameters();
    //! The number of particle on all Region(ProcessorID) in Param.
    int nParticleTotal = parameter->GetNumParticleTotal();

    //! Sum the particle on on current Region(ProcessorID).
    for (int iZoneGlobal = 0; iZoneGlobal < nZonesGlobal; ++iZoneGlobal)
    {
        //! Judge the zones on current processor.
        int recvProcess = GetZoneProcessorID(iZoneGlobal);

        //! The loop for zone on current processor.
        if (currentProcessorID == recvProcess)
        {
            //! For each zone'solver.
            int solverID = this->GetIndex();
            ParticlePointSolverParallelDomain *solver = static_cast <ParticlePointSolverParallelDomain*> (GlobalSolvers::GetSolver(iZoneGlobal, solverID));
            Param_ParticlePointSolverParallelDomain *parameter = solver->GetControlParameters();
            ParticlePointGroup *particlePointGroup = solver->GetParticlePointGroup("particlePointGroup");

            //! Calculater the CFL on each zone.
            map<int, OnePointVariable*>::iterator iter = particlePointGroup->GetIterBegin();
            while (iter != particlePointGroup->GetIterEnd())
            {
                OnePointVariable *onePointVariable = iter->second;

                //! See for 'A Second Order Time Accurate Finite Volume Method 
                //! for Unsteady Incompressible Flow on Hybrid Unstructured Grids',
                //! on JCP 2000.
                timeCFL = CalcCFLNumber(grid, onePointVariable, particleCFLThreshold);

                particleCFL = onePointVariable->GetParticleCFL();

                particleCFLMax = MAX(particleCFL, particleCFLMax);
                timeCFLMin = MIN(timeCFLMin, timeCFL);
                iter++;
            }
        }
    }

    //! MPI.
    RDouble particelCFLMaxGlobal = 0.0;
    PH_Barrier();
    PH_AllReduce(&particleCFLMax, &particelCFLMaxGlobal, 1, MPI_MAX);

    int outIterStep = GlobalDataBase::GetIntParaFromDB("outnstep");
    int intervalStepRes = GlobalDataBase::GetIntParaFromDB("intervalStepRes");
    int initOutStep = parameter->GetInitOutStep();
    outIterStep += initOutStep;

    if (particleCFLMax >= particleCFLThreshold)
    {
        ostringstream ossError;
        ossError << " Error on particle CFL number : " << particleCFLMax << "\n";;
        ossError << " outIterStep :" << outIterStep << "\n";
        ossError << " the min dt of time : " << timeCFLMin << "\n";
        ossError << endl;
        TK_Exit::ExceptionExit(ossError);
    }

    ostringstream ossLog;
    ossLog << " particle CFL number (n,dtMin,CFL) : ";
    ossLog << outIterStep << " ";
    ossLog << timeCFLMin << " ";
    ossLog << particleCFLMax << " ";
    ossLog << "\n" << endl;

    //WriteLogFile(ossLog);
    if (0 == outIterStep % intervalStepRes)
    {
        //PrintToWindow(ossLog);
    }
}

RDouble ParticlePointSolverParallelDomain::CalcCFLNumber(Grid *grid, OnePointVariable *onePointVariable, RDouble &particleCFLThreshold)
{
    int gridType = grid->Type();
    RDouble timeCFL = 0.0;
    if (PHSPACE::STRUCTGRID == gridType)
    {
        StructGrid *structgrid = StructGridCast(grid);
        timeCFL = CalcCFLNumber(structgrid, onePointVariable, particleCFLThreshold);
    }
    else if (PHSPACE::UNSTRUCTGRID == gridType)
    {
        UnstructGrid *unstructgrid = UnstructGridCast(grid);
        timeCFL = CalcCFLNumber(unstructgrid, onePointVariable, particleCFLThreshold);
    }
    else if (PHSPACE::MIXGRID == gridType)
    {
        //! Mixing
        TK_Exit::UnexpectedVarValue(" gridType = ", PHSPACE::MIXGRID);
    }

    return timeCFL;
}

RDouble ParticlePointSolverParallelDomain::CalcCFLNumber(StructGrid *grid, OnePointVariable *onePointVariable, RDouble &particleCFLThreshold)
{
    RDouble timeCFL = 0.0;

    RDouble physicalTimeStep = GlobalDataBase::GetDoubleParaFromDB("physicalTimeStep");

    RDouble4D &xfn = *(grid->GetFaceNormalX());
    RDouble4D &yfn = *(grid->GetFaceNormalY());
    RDouble4D &zfn = *(grid->GetFaceNormalZ());

    RDouble3D &vol = *(grid->GetCellVolume());
    RDouble4D &area = *(grid->GetFaceArea());

    int nI, nJ, nK;
    int iP, jP, kP;
    grid->GetND(nI, nJ, nK);
    nI -= 1;
    nJ -= 1;
    nK -= 1;
    if (GetDim() == TWO_D)
    {
        nK = 1;
    }
    int cellID = onePointVariable->GetCellID();
    GetIndexCellIJKStruct(nI, nJ, nK, iP, jP, kP, GetNumberOfGhostCellLayers(), cellID);

    SPDouble &particleVelocity = *(onePointVariable->GetOnePointVelocity());
    SPDouble &particleCoordinate = *(onePointVariable->GetOnePointCoordinate());

    const int nDimVar = 3;
    const int nFace = GetDim() * 2;

    int iSurLocal, jSurLocal, kSurLocal;

    RDouble particleCFL = 0.0;

    //! iFace in 2D : i:0,1, j:2,3.
    //! iFace in 3D : i:0,1, j:2,3, k:4,5.
    //! iDir : thr direction of face in i,j,k.
    //! ilr : Left or right of boundary condition face (0 or 1).
    //! ilrDir : seem to ilr, 1 left, -1 right, 
    //!          which is to make the face normal always point into the cell.
    int iDir = 0;
    int ilr = 0;
    int ilrDir = 0;
    for (int iFace = 0; iFace < nFace; ++iFace)
    {
        for (int iDim = 1; iDim <= GetDim(); ++iDim)
        {
            //! The direction of boundary condition face.
            if (iFace <= 2 * iDim - 1 && iFace > 2 * (iDim - 1) - 1)
            {
                iDir = iDim;
                ilr = iFace - (2 * (iDim - 1) - 1);
                //! ilr : 0,left, 1,right.
                ilr = ilr - 1;
                if (ilr < 0 || ilr > 1)
                {
                    ostringstream ossError;
                    ossError << " error on ilr = " << ilr << endl;
                    TK_Exit::ExceptionExit(ossError);
                }
                else
                {
                    if (0 == ilr)
                    {
                        ilrDir = 1;
                    }
                    else
                    {
                        ilrDir = -1;
                    }
                }
            }
        }

        //! Get the surface (or parallel) index.
        GetNsurfIndex(iDir, iSurLocal, jSurLocal, kSurLocal);
        iSurLocal = iP + ilr * iSurLocal;
        jSurLocal = jP + ilr * jSurLocal;
        kSurLocal = kP + ilr * kSurLocal;

        RDouble fn[nDimVar];
        //! Out side.
        fn[0] = xfn(iP, jP, kP, iDir) * ilrDir * (-1);
        fn[1] = yfn(iP, jP, kP, iDir) * ilrDir * (-1);
        fn[2] = zfn(iP, jP, kP, iDir) * ilrDir * (-1);

        RDouble uFaceNormal = 0.0;

        for (int iDim = 0; iDim < nDimVar; ++iDim)
        {
            uFaceNormal += fn[iDim] * particleVelocity[iDim];
        }

        uFaceNormal = ABS(uFaceNormal);

        particleCFL += uFaceNormal * area(iP,jP,kP, iDir);
    }

    particleCFL = particleCFL * (0.5) / vol(iP, jP, kP);

    timeCFL = particleCFLThreshold / particleCFL;

    particleCFL *= physicalTimeStep;

    onePointVariable->SetParticleCFL(particleCFL);

    return timeCFL;
}

RDouble ParticlePointSolverParallelDomain::CalcCFLNumber(UnstructGrid *grid, OnePointVariable *onePointVariable, RDouble &particleCFLThreshold)
{
    RDouble timeCFL = 0.0;
    return timeCFL;
}

void ParticlePointSolverParallelDomain::InitSource(Grid *grid)
{
    int gridType = grid->Type();
    if (PHSPACE::STRUCTGRID == gridType)
    {
        StructGrid *structgrid = StructGridCast(grid);

        RDouble4D &sourceParticle2Flow = *reinterpret_cast <RDouble4D*> (grid->GetDataPtr("sourceParticle2Flow"));
        sourceParticle2Flow = 0.0;
    }
    else if (PHSPACE::UNSTRUCTGRID == gridType)
    {
        UnstructGrid *unstructgrid = UnstructGridCast(grid);
    }
    else if (PHSPACE::MIXGRID == gridType)
    {
        //! Mixing
        TK_Exit::UnexpectedVarValue(" gridType = ", PHSPACE::MIXGRID);
    }
}

void ParticlePointSolverParallelDomain::GetSourceForRestart()
{
    Param_ParticlePointSolverParallelDomain *parameter = GetControlParameters();

    if (parameter->ifUseParticleBackCouping())
    {
        ParticlePointGroup *particlePointGroup = GetParticlePointGroup("particlePointGroup");

        Grid *grid = GetGrid();

        Data_Param *particleParam = new Data_Param;
        parameter->GetRefDataParam(particleParam);
        parameter->GetForceTypeDataParam(particleParam);

        map<int, OnePointVariable*>::iterator iter = particlePointGroup->GetIterBegin();
        while (iter != particlePointGroup->GetIterEnd())
        {
            OnePointVariable *onePointVariable = iter->second;

            GetSourceFromParticle2Flow(grid, onePointVariable, particleParam);
            iter++;
        }

        delete particleParam;    particleParam = nullptr;
    }
}

void ParticlePointSolverParallelDomain::GetSourceFromParticle2Flow(Grid *grid, OnePointVariable *onePointVariable, Data_Param *particleParam)
{
    Param_ParticlePointSolverParallelDomain *parameter = GetControlParameters();

    if (parameter->ifUseParticleBackCouping())
    {
        int gridType = grid->Type();
        if (PHSPACE::STRUCTGRID == gridType)
        {
            StructGrid *structgrid = StructGridCast(grid);
            GetSourceFromParticle2Flow(structgrid, onePointVariable, particleParam);
        }
        else if (PHSPACE::UNSTRUCTGRID == gridType)
        {
            UnstructGrid *unstructgrid = UnstructGridCast(grid);
            GetSourceFromParticle2Flow(unstructgrid, onePointVariable, particleParam);
        }
        else if (PHSPACE::MIXGRID == gridType)
        {
            //! Mixing
            TK_Exit::UnexpectedVarValue(" gridType = ", PHSPACE::MIXGRID);
        }
    }
}

void ParticlePointSolverParallelDomain::GetSourceFromParticle2Flow(StructGrid *grid, OnePointVariable *onePointVariable, Data_Param *particleParam)
{
    int nLayers = 0;
    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, nLayers);

    RDouble4D &sourceParticle2Flow = *reinterpret_cast <RDouble4D*> (grid->GetDataPtr("sourceParticle2Flow"));

    RDouble4D &q = *reinterpret_cast <RDouble4D*> (grid->GetDataPtr("q"));

    RDouble particleDiameter = onePointVariable-> GetOnePointDiameter();
    RDouble particleDensity = onePointVariable->GetOnePointDensity();
    RDouble particleSpecificHeatCapacity = onePointVariable->GetOnePointSpecificHeatCapacity();
    RDouble particleEmissivity = onePointVariable->GetOnePointEmissivity();

    SPDouble &particleVelocity = *(onePointVariable->GetOnePointVelocity());
    SPDouble &particleCoordinate = *(onePointVariable->GetOnePointCoordinate());
    SPDouble &particleAngularVelocity = *(onePointVariable->GetOnePointAngularVelocity());
    SPDouble &particleAcceleration = *(onePointVariable->GetOnePointAcceleration());
    SPDouble &particleTemperature = *(onePointVariable->GetOnePointTemperature());
    SPDouble &flowOnParticle = *(onePointVariable->GetOnePointFlowOnParticle());
    SPDouble &particleForce = *(onePointVariable->GetOnePointForce());

    int nI, nJ, nK;
    int iP, jP, kP;
    grid->GetND(nI, nJ, nK);
    nI -= 1;
    nJ -= 1;
    nK -= 1;
    if (TWO_D == GetDim())
    {
        nK = 1;
    }
    int cellID = onePointVariable->GetCellID();
    GetIndexCellIJKStruct(nI, nJ, nK, iP, jP, kP, GetNumberOfGhostCellLayers(), cellID);

    int nDimVar;

    //! force.
    using namespace INDEX_PARTICLE_FORCE;
    nDimVar = 3;
    RDouble totalForce[3];
    totalForce[0] = 0;
    totalForce[1] = 0;
    totalForce[2] = 0;

    for (int iForce = 0; iForce < nForceDir / nDimVar; ++iForce)
    {
        totalForce[0] += particleForce[iForce * nDimVar + 0];
        totalForce[1] += particleForce[iForce * nDimVar + 1];
        totalForce[2] += particleForce[iForce * nDimVar + 2];
    }

    for (int iDim = 0; iDim < nDimVar; ++iDim)
    {
        sourceParticle2Flow(iP, jP, kP, iDim) -= totalForce[iDim];
    }

    //! temperature.
    using namespace PARTICLE_PARAM;
    using namespace INDEX_FLOWONPARTICLE;

    RDouble R = GetParticleDoublePara(particleParam, "refAverageGeneralGasConstantDimensional");
    RDouble refGama = GetParticleDoublePara(particleParam, "refGama");
    RDouble refReNumber = GetParticleDoublePara(particleParam, "refReNumber");
    RDouble refTemperatureDimensional = GetParticleDoublePara(particleParam, "refTemperatureDimensional");
    RDouble refVelocityDimensional = GetParticleDoublePara(particleParam, "refVelocityDimensional");
    RDouble refDynamicViscosityDimensional = GetParticleDoublePara(particleParam, "refDynamicViscosityDimensional");
    RDouble refDensityDimensional = GetParticleDoublePara(particleParam, "refDensityDimensional");

    RDouble prl = GetParticleDoublePara(particleParam, "prl");

    RDouble radiativeTemperature = GetParticleDoublePara(particleParam, "radiativeTemperature");

    int enableRadiativeHeatTransfer = GetParticleIntPara(particleParam, "enableRadiativeHeatTransfer");
    int enableVariableDiameter = GetParticleIntPara(particleParam, "enableVariableDiameter");

    RDouble flowHeatCPConstPressureDim = GetParticleDoublePara(particleParam, "flowHeatCPConstPressureDim");
    RDouble flowHeatCPConstPressureDimless = GetParticleDoublePara(particleParam, "flowHeatCPConstPressureDimless");

    RDouble cd = particleSpecificHeatCapacity;
    RDouble cp = flowHeatCPConstPressureDim;

    //! Init and Calculate the relative velocity of particle.
    int nDim = 3;
    SPDouble relativeVelocity(nDim);
    relativeVelocity[0] = flowOnParticle[PFU] - particleVelocity[0];
    relativeVelocity[1] = flowOnParticle[PFV] - particleVelocity[1];
    relativeVelocity[2] = flowOnParticle[PFW] - particleVelocity[2];

    onePointVariable->SetRelativeVelocity(&relativeVelocity[0], &relativeVelocity[1], &relativeVelocity[2]);

    //! Calculate the Reynolds number of particle.
    RDouble relativeVelocityNorm = sqrt(SQR(relativeVelocity[0], relativeVelocity[1], relativeVelocity[2]));
    if (std::isnan(flowOnParticle[PFMU]))
    {
        cout << "flow on particle Error!!" << endl;
    }
    RDouble particleRe = refReNumber * flowOnParticle[PFR] * relativeVelocityNorm * particleDiameter / (flowOnParticle[PFMU]+SMALL);

    RDouble Nu = 2.0;

    //Total force stokes number!
    RDouble totalForceNorm = sqrt(SQR(totalForce[0], totalForce[1], totalForce[2]));
    RDouble Cdcoeff = totalForceNorm / (flowOnParticle[PFR] * relativeVelocityNorm * relativeVelocityNorm * 3.1416 * particleDiameter * particleDiameter / 8);
    RDouble Cds = 24 / particleRe;

    RDouble taup = refReNumber * particleDensity * pow(particleDiameter, 2.0) / 18.0 / flowOnParticle[PFMU];
    RDouble TotalForceSt = taup * relativeVelocityNorm / particleDiameter * (Cds / Cdcoeff);

    Nu = 2.0 + 0.6 * pow(prl, 1.0 / 3.0) * pow(particleRe, 0.5);

    onePointVariable->SetParticleTotalForceSt(TotalForceSt);

    //! the work part of the particle.
    RDouble particleWork = particleVelocity[0] * totalForce[0] + particleVelocity[1] * totalForce[1] + particleVelocity[2] * totalForce[2];
    RDouble particleDissipation= relativeVelocity[0] * totalForce[0] + relativeVelocity[1] * totalForce[1] + relativeVelocity[2] * totalForce[2];

    //! kc is the coefficient of thermal conductivity of  flow.
    //! kc is dim.
    RDouble kc = cp * refDynamicViscosityDimensional * flowOnParticle[PFMU] / prl;
    kc = kc / (refDensityDimensional * pow(refVelocityDimensional, 3));
    RDouble particleMass = (4.0 / 3.0) * PI * pow(particleDiameter * 0.5, 3.0) * particleDensity;

    RDouble dparticleTemperature = 6.0 * kc *
        (flowOnParticle[PFT] - particleTemperature[0]) * Nu
        / (particleDensity * pow(particleDiameter, 2.0) * cd);

    //! the ConvectiveHeatTransfer part of the particle.
    RDouble convectiveHeatTransfer = dparticleTemperature * particleMass * cd;

    onePointVariable->SetParticleWork(particleWork);
    onePointVariable->SetParticleDissipation(particleDissipation);

    RDouble energyTotal = convectiveHeatTransfer + particleWork;
    if (std::isnan(energyTotal))
    {
        RDouble ID =onePointVariable->GetParticleID();
        cout <<"The NO. "<< ID << "  EnergyTotal back to flow field error !    ##############" << endl;
        sourceParticle2Flow(iP, jP, kP, 3) -= 0;
    }
    else
    {
        sourceParticle2Flow(iP, jP, kP, 3) -= energyTotal;
    }
    

    ostringstream ossLog;
    nDimVar = 4;

    if (1 == enableVariableDiameter)
    {
        cout << "Variable particle diameter is activated!" << endl;
        particleDiameter -= 1e-7;
        onePointVariable->SetDiameter(particleDiameter);
    }
}

void ParticlePointSolverParallelDomain::GetSourceFromParticle2Flow(UnstructGrid *grid, OnePointVariable *onePointVariable, Data_Param *particleParam)
{

}

//! @author ZHANGSHH
void ParticlePointSolverParallelDomain::KernelBoundaryTransformation(RDouble fn[3], RDouble particleVelocity[3],RDouble parVelNew[3])
{
    //! ###############################
    //! #
    //! #Boundary kernel transformation
    //! #
    //! ###############################

    const int nDimVar = 3;
    RDouble particleVel[nDimVar];
    for (int iDim = 0; iDim < nDimVar; ++iDim)
    {
        particleVel[iDim] = particleVelocity[iDim];
    }
    RDouble homoFn[nDimVar + 1];

    ConvertVectorToHomogeneous(fn, homoFn);

    RDouble homoVel[nDimVar + 1];

    ConvertVectorToHomogeneous(particleVel, homoVel);

    RDouble reflectMatrix[nDimVar + 1][nDimVar + 1] = { 0 };
    reflectMatrix[0][0] = 1;
    for (int iDim = 1; iDim < 4; ++iDim)
    {
        for (int jDim = 1; jDim < 4; ++jDim)
        {
            reflectMatrix[iDim][jDim] = homoFn[iDim] * homoFn[jDim];
        }
    }

    for (int iDim = 1; iDim < 4; ++iDim)
    {
        for (int jDim = 1; jDim < 4; ++jDim)
        {
            if (iDim == jDim)
            {
                reflectMatrix[iDim][jDim] = 1 - 2 * reflectMatrix[iDim][jDim];
            }
            else
            {
                reflectMatrix[iDim][jDim] *= -2;
            }
        }
    }

    RDouble homogeneousReflectVel[nDimVar + 1] = { 0 };
    for (int iDim = 0; iDim < nDimVar + 1; ++iDim)
    {
        for (int jDim = 0; jDim < nDimVar + 1; ++jDim)
        {
            homogeneousReflectVel[iDim] += reflectMatrix[iDim][jDim] * homoVel[jDim];
        }
    }

    for (int iDim = 0; iDim < nDimVar; ++iDim)
    {
        parVelNew[iDim]= homogeneousReflectVel[iDim+1];
    }
    //! ###########################################
    //! #                                         #
    //! #The end of Boundary kernel transformation#
    //! #                                         #
    //! ###########################################
}
void ParticlePointSolverParallelDomain::ConvertVectorToHomogeneous(RDouble originVector[3], RDouble homoVector[4])
{
    homoVector[0] = 1;
    for (int iDim = 0; iDim < 3; ++iDim)
    {
        homoVector[iDim + 1] = originVector[iDim];
    }
}
void ParticlePointSolverParallelDomain::MirrorCoordinate(RDouble originPoint[3], RDouble wallDistance, RDouble fn[3], RDouble mirrorPoint[3])
{
    RDouble homoPoint[4];
    RDouble homoPoint1[4] = { 0 };
    ConvertVectorToHomogeneous(originPoint, homoPoint);
    RDouble mirrorMatrix[4][4] = { 0 };

    for (int iDim = 0; iDim < 4; ++iDim)
    {
        for (int jDim = 0; jDim < 4; ++jDim)
        {
            if (iDim == jDim)
            {
                mirrorMatrix[iDim][jDim] = 1.0;
            }
        }
    }

    for (int iDim = 1; iDim < 4; ++iDim)
    {
        mirrorMatrix[iDim][0] = -2 * wallDistance * fn[iDim - 1];
    }

    for (int iDim = 0; iDim < 4; ++iDim)
    {
        for (int jDim = 0; jDim < 4; ++jDim)
        {
            homoPoint1[iDim] += mirrorMatrix[iDim][jDim] * homoPoint[jDim];
        }
    }
    for (int iDim = 0; iDim < 3; ++iDim)
    {
        mirrorPoint[iDim] = homoPoint1[iDim + 1];
    }
}
RDouble ParticlePointSolverParallelDomain::DistanceBetweenPointAndSurface(RDouble point[3], RDouble surfacePoint[3], RDouble surfn[3])
{
    RDouble surEqD = -surfn[0] * surfacePoint[0] - surfn[1] * surfacePoint[1] - surfn[2] * surfacePoint[2];
    RDouble distance = 0;
    for (int iDim = 0; iDim < 3; ++iDim)
    {
        distance += surfn[iDim] * point[iDim];
    }
    distance += surEqD;
    distance = ABS(distance);
    distance /= DISTANCE(surfn[0], surfn[1], surfn[2]);
    return distance;
}

}