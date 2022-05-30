#include "Controller.h"
#include "CFDSolver.h"
#include "Region.h"
#include "Math_BasisFunction.h"
#include "Task.h"
#include "TK_Log.h"
using namespace std;

namespace PHSPACE
{
Controller::Controller(int iSolver, Region *region)
{
    this->solverID = iSolver;
    this->region = region;

    q_proxy = 0;

    using namespace PHMPI;

    int nZones = GetNumberofGlobalZones();

    rhs_proxy = new ZoneFieldProxy(iSolver, GetNumberOfMultiGrid(), nZones);

    if (GetNumberOfMultiGrid() > 1)
    {
        q_proxy = new ZoneFieldProxy(iSolver, GetNumberOfMultiGrid(), nZones);
    }
}

Controller::~Controller()
{
    if (GetNumberOfMultiGrid() > 1)
    {
        delete q_proxy;
    }

    delete rhs_proxy;
}

LIB_EXPORT void Controller::InitCoarseGridsFlow(int iSolver, int level)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Zone *zone = GetZone(zoneID);
        if (!zone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        Grid *coarseGrid = grid->GetCoarseGrid();
        if (coarseGrid)
        {
            if (iSolver == TURB_SOLVER) continue;
            CFDSolver *solver = GetCFDSolver(zoneID, iSolver);
            //solver->InitCGrid(grid, coarseGrid);
            solver->ComputeGamaAndTemperature(coarseGrid);
            solver->Boundary(coarseGrid);
            solver->ComputeGamaAndTemperature(coarseGrid);
            solver->ComputeNodeValue(coarseGrid);
        }
    }
}

LIB_EXPORT void Controller::Initialize(int iSolver)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();
    int nProtectData = 0;
    if(GlobalDataBase::IsExist("nProtectData", PHINT, 1))
    {
        nProtectData = GlobalDataBase::GetIntParaFromDB("nProtectData");
    }

    int level = 0;
    //! First: Initialize memory.
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Zone *zone = GetZone(zoneID);
        if (!zone)
        {
            continue;
        }

        PHSolver *solver = GetCFDSolver(zoneID, iSolver);

        solver->InitMemory();
    }

    //! Second: Initialize flow variables on the finest mesh.
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Zone *zone = GetZone(zoneID);
        if (!zone)
        {
            continue;
        }

        PHSolver *solver = GetCFDSolver(zoneID, iSolver);

        if(nProtectData == 0)
        {
            solver->InitFlow();
        }

        else
        {
            solver->InitProtectedFlow();
        }
    }

    //! Third: Initialize flow variables on coarser grids, if multi-grid is used.
    int nMGLevel = PHSPACE::GlobalDataBase::GetIntParaFromDB("nMGLevel");
    for (int iLevel = 0; iLevel < nMGLevel - 1; ++ iLevel)
    {
        RestrictAllQ(iSolver, iLevel);
        InitCoarseGridsFlow(iSolver, iLevel);
    }

    //! Finally: Interfaces communications.
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Zone *zone = GetZone(zoneID);
        if (!zone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        PHSolver *solver = GetCFDSolver(zoneID, iSolver);

        solver->CommunicationInterfaceData();
        InterpointInformation *interPointInfor = grid->GetInterpointInfo();
        if (interPointInfor)
        {
            solver->CommunicationInterpointData();
        }
    }

    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);
        Zone *zone = GetZone(zoneID);
        if (!zone) continue;

        for (int iLevel = 0; iLevel < GetNumberOfMultiGrid(); iLevel ++)
        {
            Grid *grid = GetGrid(zoneID, iLevel);
            CFDSolver *solver = GetCFDSolver(zoneID, iSolver);
            solver->ModifyNodeValue(grid);
            solver->ModifyQTurbNodeValue(grid);
        }
    }

    int compressible = GlobalDataBase::GetIntParaFromDB("compressible");
    if (compressible)
    {
    ComputeNodeBCType();
}
}

void Controller::CleanUp(int iSolver)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

    int level = 0;
    //! First: Initialize memory.
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Zone *zone = GetZone(zoneID);
        if (!zone)
        {
            continue;
        }

        PHSolver *solver = GetCFDSolver(zoneID, iSolver);

        solver->ReleaseMemory();
    }
}

void Controller::ReinitField(int iSolver)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

    int level = 0;
    //! First: Initialize memory.
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Zone *zone = GetZone(zoneID);
        if (!zone)
        {
            continue;
        }

        PHSolver *solver = GetCFDSolver(zoneID, iSolver);

        solver->ReinitField();
    }
}

void Controller::ComputeNodeBCType()
{
    int sysGridType = GlobalDataBase::GetIntParaFromDB("sys_gridtype");
    if (sysGridType != UNSTRUCTGRID)
    {
        return;
    }

    using namespace PHMPI;

    int needCommunication = 1;
    int nZones = GetNumberofGlobalZones();
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int iZoneProcessorID = GetZoneProcessorID(iZone);
        int myid = GetCurrentProcessorID();
        if (myid == iZoneProcessorID)
        {
            //! Compute node bctype in each zone.
            Grid *grid = GetGrid(iZone, 0);
            grid->ComputeNodeBCType();
            if (!(grid->GetInterpointInfo()))
            {
                needCommunication = 0;
            }
        }
    }

    PH_CompareMaxMin(needCommunication, 2);
    if (!needCommunication)
    {
        return;
    }
    CommunicationNodeBCType();
}

void Controller::CommunicationNodeBCType()
{
    using namespace PHMPI;
    using namespace PHENGLEI;

    map<int, int> bcTypeMap;
    bcTypeMap.insert(pair<int,int>(SOLID_SURFACE   , 4));
    bcTypeMap.insert(pair<int,int>(FARFIELD        , 3));
    bcTypeMap.insert(pair<int,int>(SYMMETRY        , 2));
    bcTypeMap.insert(pair<int,int>(INTERFACE       , 1));
    bcTypeMap.insert(pair<int,int>(INFLOW          , 1));
    bcTypeMap.insert(pair<int,int>(OUTFLOW         , 1));
    bcTypeMap.insert(pair<int,int>(PRESSURE_INLET  , 1));
    bcTypeMap.insert(pair<int,int>(PRESSURE_OUTLET , 1));
    bcTypeMap.insert(pair<int,int>(OUTFLOW_CONFINED, 1));

    int nZones = GetNumberofGlobalZones();
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int myid = GetCurrentProcessorID();
        int tag = iZone;
        int sendProc = GetZoneProcessorID(iZone);

        ZoneNeighbor *zoneNeighborForPoint = zoneConnectivityForPoint->GetZoneNeighborForPoint(iZone);
        std::size_t numberOfNeighbor = zoneNeighborForPoint->GetNumberOfNeighbor();

        for (std::size_t iNeighborZone = 0; iNeighborZone < numberOfNeighbor; ++ iNeighborZone)
        {
            int neighborZoneIndex = zoneNeighborForPoint->GetZoneIndexOfNeighbor(iNeighborZone);
            int recv_proc = GetZoneProcessorID(neighborZoneIndex);

            DataContainer *cdata = new DataContainer();

            if (myid == sendProc)
            {
                UnstructGrid *grid = UnstructGridCast(PHSPACE::GetGrid(iZone, 0));
                InterpointInformation *interpointInformation = grid->GetInterpointInfo();
                int *nodeBCType = reinterpret_cast <int *> (grid->GetDataPtr("nodeBCType"));

                cdata->moveto_begin();

                int iNeighbor = interpointInformation->FindIthNeighbor(neighborZoneIndex);
                int numberOfInterpointsForNeighbor = interpointInformation->GetNumberOfInterpointsForNeighbor(iNeighbor);
                int *pointIndexForSend = interpointInformation->GetPointIndexForSend(iNeighbor);
                int *interPoint2GlobalPoint = interpointInformation->GetInterPoint2GlobalPoint();

                int iPoint;
                int globalPoint;
                int BCType;
                for (int ipointLocal = 0; ipointLocal < numberOfInterpointsForNeighbor; ++ ipointLocal)
                {
                    iPoint = pointIndexForSend[ipointLocal];
                    globalPoint = interPoint2GlobalPoint[iPoint];
                    BCType = nodeBCType[globalPoint];

                    cdata->write(&BCType, sizeof(int));
                }
            }

            PH_Trade(cdata, sendProc, recv_proc, tag);

            if (myid == recv_proc)
            {
                UnstructGrid *gridNeighbor = UnstructGridCast(PHSPACE::GetGrid(neighborZoneIndex, 0));
                InterpointInformation *interpointInformationNeighbor = gridNeighbor->GetInterpointInfo();
                int *nodeBCType = reinterpret_cast<int *>(gridNeighbor->GetDataPtr("nodeBCType"));

                cdata->moveto_begin();

                int iNeighbor = interpointInformationNeighbor->FindIthNeighbor(iZone);
                int numberOfInterpointsForNeighbor = interpointInformationNeighbor->GetNumberOfInterpointsForNeighbor(iNeighbor);
                int *pointIndexForReceive = interpointInformationNeighbor->GetPointIndexForReceive(iNeighbor);
                int *interPoint2GlobalPoint = interpointInformationNeighbor->GetInterPoint2GlobalPoint();

                int iPoint;
                int globalPoint;
                int recvBCType;
                for (int ipointLocal = 0; ipointLocal < numberOfInterpointsForNeighbor; ++ ipointLocal)
                {
                    iPoint = pointIndexForReceive[ipointLocal];
                    globalPoint = interPoint2GlobalPoint[iPoint];
                    cdata->read(&recvBCType, sizeof(int));

                    int localBCIndex = bcTypeMap[nodeBCType[globalPoint]];
                    int recvBCIndex  = bcTypeMap[recvBCType];
                    if (recvBCIndex > localBCIndex)
                    {
                        nodeBCType[globalPoint] = recvBCType;
                    }

                }
            }

            delete cdata;
        }
    }
}

LIB_EXPORT void Controller::SolveIncomFlow()
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();
    int iSolver = this->GetSolverIndex();
#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Zone *zone = GetZone(zoneID);
        if (!zone)
        {
            continue;
        }
        
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        solver->Solve();
    }

}

LIB_EXPORT void Controller::SolveSteadyField()
{
    int iSolver = this->GetSolverIndex();

    //! Zero the residuals of the finest grid.
    ZeroResiduals(iSolver);

    MultiGrid(iSolver);

}

LIB_EXPORT void Controller::SolveUnsteadyField()
{
    int iSolver = this->GetSolverIndex();
    int inerStepInterval = GlobalDataBase::GetIntParaFromDB("innstep");

    while (true)
    {
        inerStepInterval += 1;
        GlobalDataBase::UpdateData("innstep", &inerStepInterval, PHINT, 1);

        SolveSteadyField();

        if (CVGNorm(iSolver))
        {
            UpdateUnsteadyFlow(iSolver);
            break;
        }

        PHSolver *solver = GetCFDSolver(0, iSolver);
        solver->CommunicationInterpointData();
    }
}

//! Incompressible Flow
void Controller::SolveIncompressibleSteadyField()
{
    using namespace PHMPI;

    int iSolver = this->GetSolverIndex();
    int nLocalZones = GetNumberofLocalZones();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Zone *zone = GetZone(zoneID);
        if (!zone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, 0);


        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);
        const string solverName = solver->GetName();
        if (solverName == std::to_string(iSolver) + "_" + "FLOW")
        {
            this->SolveFlowField(zoneID);
        }
        else if (solverName == "TURB")
        {
            this->SolveTurbField(zoneID);
        }
        else if (solverName == "ENERGY")
        {
            this->SolveEnergyField(zoneID);
        }
    }
}

void Controller::InitResultFile(string &fileName)
{
    using namespace PHMPI;

    int iSolver = this->GetSolverIndex();
    int nLocalZones = GetNumberofLocalZones();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Zone *zone = GetZone(zoneID);
        if (!zone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, 0);

        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        solver->InitIncompressibleResultFile(grid, fileName);
    }
}

void Controller::DumpResultFile(string &fileName)
{
    using namespace PHMPI;

    int iSolver = this->GetSolverIndex();
    int nLocalZones = GetNumberofLocalZones();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Zone *zone = GetZone(zoneID);
        if (!zone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, 0);

        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        solver->DumpIncompressibleResultFile(grid, fileName);
    }
}

void Controller::MultiGrid(int iSolver, int level)
{
    CreateRHSProxy(iSolver, level);
    StoreRhsByResidual(iSolver, level);
    int outerStepInterval = GlobalDataBase::GetIntParaFromDB("outnstep");
    int MGPreIteration    = GlobalDataBase::GetIntParaFromDB("MGPreIteration");
    bool computeOnFinestGrid = (outerStepInterval % MGPreIteration != 0);

    //! If grid is in the coarsest level.
    if (ISCoarsestGrid(level) || iSolver == TURB_SOLVER || computeOnFinestGrid)
    {
        Relaxation(iSolver, level);
    }
    else
    {
        int cLevel = level + 1;

        //! Smooth the solutions before corrections.
        Relaxation(iSolver, level);

        //! Restrict from fine to coarse grid for all q.
        //! Corresponding to nsmb: wsav = restr(w).
        //! Call the following function to update the Q value(cq) on the coarse grid.
        RestrictAllQ(iSolver, level);

        //! Store coarse grid solution in temporary array cqsav.
        //! cqsav represents the variables on the coarse grid.
        CreateQProxy(iSolver, cLevel);

        //! After LoadQ(cgrid, cqsav), cqsav is actually wsav(Corresponding to NSMB).
        //! LoadQ  gets out the value of q(cq) on coarse grid and assigns it to cqsav.
        LoadQ(iSolver, cLevel);

        //! Defect.
        //! The followings corresponding to nsmb d = RL-1(wsav) - restr(RL(w)-f).
        //! in which (*res) = - (*rhs) correspond to - f
        //! Note that the value of res is on fine grid, rhs ought be f£¬and this should be checked.
        //! Here should be ni and nj,rather than cni and cnj!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        InitResidual(iSolver, level);

        //! After UpdateResiduals,(*res) = RL(w) - f.
        UpdateResiduals(iSolver, level);

        //! After RestrictDefect,(*cres) = - restr(RL(w) - f).
        RestrictDefect(iSolver, level);

        //! After UpdateResiduals,(*cres) = RL-1(wsav) - restr(RL(w) - f).
        UpdateResiduals(iSolver, cLevel);

        //! Solve on coarse grid.
        int MGFasType = GlobalDataBase::GetIntParaFromDB("MGFasType");
        int ifLowSpeedPrecon = GlobalDataBase::GetIntParaFromDB("ifLowSpeedPrecon");
        if (ifLowSpeedPrecon == 1)
        {
            MGFasType = 1;    //! Precondition method is only available for V cycle.
        }

        //! After MultiGrid(cgrid), the value of q on coarse grid has been updated.
        //! At this time,the value of q on coarse grid equals to w0 of NSMB.
        //! MGFasType = 1 V cycle
        //! MGFasType = 2 W cycle
        for (int i = 0; i < MGFasType; ++ i)
        {
            MultiGrid(iSolver, cLevel);
        }


        //! Put correction back to the coarse grid.
        //! Actually, PutCorrection has done the operation of (w0 - wsav).
        PutCorrection(iSolver, cLevel);


        //! Correct variables in fine grid.
        //! Actually CorrectFineGrid has done the operation of w = w + prol(w0 - wsav).
        CorrectFineGrid(iSolver, level);

        //! Smooth the solutions after corrections.
        int n_post = GlobalDataBase::GetIntParaFromDB("n_post");
        for (int iter = 0; iter < n_post; ++ iter)
        {
            Relaxation(iSolver, level);
        }

        //! The following function actually recover the value of q on coarse grid.
        PutCorrectionBack(iSolver, cLevel);
        //! Delete the temporary variables of QProxy.
        DestroyQProxy(iSolver, cLevel);
    }
    //! Put RHS back into its places.
    //! Here, we just recover the initial value of res.
    RecoverResidual(iSolver, level);

    //! Delete the temporary variables of rhs_proxy.
    DestroyRHSProxy(iSolver, level);
}

void Controller::UpdateUnsteadyFlow(int iSolver, int level)
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        Zone *zone = GetZone(iZone);
        if (!zone)
        {
            continue;
        }

        Grid *grid = GetGrid(iZone, level);
        CFDSolver *solver = GetCFDSolver(iZone, iSolver);

        solver->UpdateUnsteadyFlow(grid);
    }
}

void Controller::UpdateUnsteadyVariable(int iSolver, int level)
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        Zone *zone = GetZone(iZone);
        if (!zone)
        {
            continue;
        }

        Grid *grid = GetGrid(iZone, level);
        CFDSolver *solver = GetCFDSolver(iZone, iSolver);

        solver->UpdateUnsteadyVariable(grid);
    }
}

void Controller::UpdateDt(int iSolver, int level)
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        Zone *zone = GetZone(iZone);
        if (!zone)
        {
            continue;
        }

        Grid *grid = GetGrid(iZone, level);
        CFDSolver *solver = GetCFDSolver(iZone, iSolver);

        solver->UpdateDt();
    }
}

void Controller::Relaxation(int iSolver, int level)
{
    using namespace PHSPACE;
    int isFVMOrFDM = GlobalDataBase::GetIntParaFromDB("isFVMOrFDM");
    if (isFVMOrFDM == FD_METHOD)
    {
        PreSolveforStructHighOrder(iSolver, level);
    }

    int isUseLocalCFL = 0;
    if (GlobalDataBase::IsExist("isUseLocalCFL", PHINT, 1))
    {
        isUseLocalCFL = GlobalDataBase::GetIntParaFromDB("isUseLocalCFL");
    }

    int tscheme = GlobalDataBase::GetIntParaFromDB("tscheme");
    if (tscheme == LU_SGS || level != 0)
    {
        //! LU-SGS implicit time integration method.
        //if(!runOriginalSingleSweep)
        {
            LoadResiduals(iSolver, level);

            UpdateResiduals(iSolver, level);

            SpectrumRadius(iSolver, level);

            TimeStep(iSolver, level);

            if (isUseLocalCFL == 1)
            {
                DualTimeStepsLUSGS(iSolver, level);
            }
            else
            {
                LUSGS(iSolver, level);
            }
        }
    }
    else if (tscheme == MULTI_STAGE)
    {
        //! Explicit multi-stage time marching method.
        SetBoundaryValueOfField(iSolver, level);

        // To obtain the viscous coff in spectrum radium term
        ComputeViscousCoeff(iSolver, level);

        SpectrumRadius(iSolver, level);

        TimeStep(iSolver, level);

        RungeKutta(iSolver, level);
    }
    else if (tscheme == Block_LU_SGS)
    {
        LoadResiduals(iSolver, level);

        UpdateResiduals(iSolver, level);

        SpectrumRadius(iSolver, level);

        TimeStep(iSolver, level);

        BLUSGS(iSolver, level);
    }
    else if (tscheme == Line_LU_SGS)
    {
        LoadResiduals(iSolver, level);

        UpdateResiduals(iSolver, level);

        SpectrumRadius(iSolver, level);

        TimeStep(iSolver, level);

        LineLUSGS(iSolver, level);
    }
	else if (tscheme == Matrix_LU_SGS)
    {
        LoadResiduals(iSolver, level);

        UpdateResiduals(iSolver, level);

        SpectrumRadius(iSolver, level);

        TimeStep(iSolver, level);

        MLUSGS(iSolver, level);
    }
    else
    {
        //! default
        LoadResiduals(iSolver, level);
        
        UpdateResiduals(iSolver, level);
        
        SpectrumRadius(iSolver, level);
        
        TimeStep(iSolver, level);
        
        if (isUseLocalCFL == 1)
        {
            DualTimeStepsLUSGS(iSolver, level);
        }
        else
        {
            LUSGS(iSolver, level);
        }
    }
}

void Controller::TimeStep(int iSolver, int level)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

    if (iSolver != NS_SOLVER)
    {
        //! Only get time step for NS solver.
        return;
    }

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Zone *zone = GetZone(zoneID);
        if (!zone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        solver->TimeStep(grid);
    }

    ReduceMaxTimeStep(iSolver, level);
}

RDouble Controller::ComputeMinTimeStep(int iSolver, int level)
{
    using namespace PHMPI;
    int outnstep      = GlobalDataBase::GetIntParaFromDB("outnstep");
    int allReduceStep = GlobalDataBase::GetIntParaFromDB("allReduceStep");
    RDouble globalMinTimeStep = GlobalDataBase::GetDoubleParaFromDB("globalMinTimeStep");

    if (outnstep % allReduceStep != 0 && globalMinTimeStep > 0.0)
    {
        //! Compute global minimum time step each 'allReduceStep' iterations.
        return globalMinTimeStep;
    }

    int nLocalZones = GetNumberofLocalZones();
    RDouble globalMinDt =  LARGE;
    RDouble localMaxDt  = -LARGE;

    //! Important: Can not use OMP directly, due to the MIN/MAX comparison!
    // #ifdef USE_OMP
    // #pragma omp parallel for
    // #endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Zone *zone = GetZone(zoneID);
        if (!zone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        RDouble minDt, maxDt;
        solver->ComputeMinTimeStep(grid, minDt, maxDt);

        globalMinDt = MIN(globalMinDt, minDt);
        localMaxDt  = MAX(localMaxDt, maxDt);
    }


    PH_CompareMaxMin(globalMinDt, 2);
    GlobalDataBase::UpdateData("globalMinTimeStep", &globalMinDt, PHDOUBLE, 1);

    return globalMinDt;
}

void Controller::ReduceMaxTimeStep(int iSolver, int level)
{
    using namespace PHMPI;
    int ifLocalTimeStep = GlobalDataBase::GetIntParaFromDB("ifLocalTimeStep");
    if (ifLocalTimeStep != 0)
    {
        //! NON local time step method, return.
        return;
    }
    if (iSolver != NS_SOLVER)
    {
        //! Only reduce time step for NS solver.
        return;
    }

    RDouble globalMinDt = ComputeMinTimeStep(iSolver, level);

    int nLocalZones = GetNumberofLocalZones();
#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Zone *zone = GetZone(zoneID);
        if (!zone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        solver->ReduceMaxTimeStep(grid, globalMinDt);
    }
}

void Controller::Diagonal(int iSolver, int level)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Zone *zone = GetZone(zoneID);
        if (!zone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        solver->Diagonal(grid);
    }
}

void Controller::SpectrumRadius(int iSolver, int level)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();
    if (iSolver != NS_SOLVER)
    {
        return;
    }

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Zone *zone = GetZone(zoneID);
        if (!zone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        solver->SpectrumRadius(grid);
    }
}

void Controller::ZeroResiduals(int iSolver, int level)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Zone *zone = GetZone(zoneID);
        if (!zone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        solver->ZeroResiduals(grid);
    }
}

void Controller::UpdateResiduals(int iSolver, int level)
{
    int systemGridType = PHSPACE::GlobalDataBase::GetIntParaFromDB("sys_gridtype");
    int isFVMOrFDM = GlobalDataBase::GetIntParaFromDB("isFVMOrFDM");
    string highordersolvername = GlobalDataBase::GetStrParaFromDB("str_highorder_solver");

    if((systemGridType == MIXGRID && iSolver == NS_SOLVER && isFVMOrFDM == FV_METHOD) || (iSolver == NS_SOLVER && isFVMOrFDM == FD_METHOD && highordersolvername.substr(0, 4) == "HDCS"))
    {
        //! ugly codes, modify later.
        //! for mixgrid, more variables need communication for conservation.
        //! InvscidFlux() and ViscousFlux() must be breaked into a lot of small functions.
        UpdateResidualsOnlyforMixGrid(iSolver, level);
        return;
    }

    using namespace PHMPI;

    SetBoundaryValueOfField(iSolver, level);

    UpdateIndependentField(iSolver, level);

    int nLocalZones = GetNumberofLocalZones();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Zone *zone = GetZone(zoneID);
        if (!zone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        solver->UpdateResiduals(grid);
    }
}

void Controller::UpdateResidualsOnlyforMixGrid(int iSolver, int level)
{
    using namespace PHMPI;

    SetBoundaryValueOfField(iSolver, level);

    UpdateIndependentField(iSolver, level);

    int nLocalZones = GetNumberofLocalZones();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Zone *zone = GetZone(zoneID);
        if (!zone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        solver->UpdateQlQrOnlyforMixGrid(grid);
    }
    CommunicateQlQrOnlyforMixGrid();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Zone *zone = GetZone(zoneID);
        if (!zone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        solver->UpdateInviscidfluxOnlyforMixGrid(grid);

        solver->UpdateViscousfluxOnlyforMixGrid(grid);
    }
    CommunicateFacefluxOnlyforMixGrid();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);
        
        Zone *zone = GetZone(zoneID);
        if (!zone)
        {
            continue;
        }
        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);
        solver->LoadGlobalfacefluxtoResOnlyforMixGrid(grid);
    }

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Zone *zone = GetZone(zoneID);
        if (!zone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        solver->SourceFlux(grid);
    }
}


void Controller::LoadResiduals(int iSolver, int level)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();
    FieldProxy **rhsProxy = GetRHSProxy(level);

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Zone *zone = GetZone(zoneID);
        if (!zone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID,level);
        CFDSolver *solver = GetCFDSolver(zoneID,iSolver);

        solver->LoadResiduals(grid, rhsProxy[zoneID]);
    }
}

void Controller::PressureFactor(int iSolver, int level)
{
    int nEntropyFixMethod = GlobalDataBase::GetIntParaFromDB("roeEntropyFixMethod");

    if (NS_SOLVER != iSolver || (nEntropyFixMethod != 2 && nEntropyFixMethod != 6))
    {
        return;
    }

    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Zone *zone = GetZone(zoneID);
        if (!zone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        solver->CalPressureFactor(grid);
    }

    CommunicateSpecificArray(iSolver, level, "rtem", 0);
}

void Controller::ComputeGradient(int iSolver, int level)
{
    int systemGridType = PHSPACE::GlobalDataBase::GetIntParaFromDB("sys_gridtype");
    int viscousType = GlobalDataBase::GetIntParaFromDB("viscousType");
    int isFVMOrFDM = GlobalDataBase::GetIntParaFromDB("isFVMOrFDM");
    if ((viscousType == INVISCID && STRUCTGRID == systemGridType))
    {
        return;
    }

    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Zone *zone = GetZone(zoneID);
        if (!zone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        int gridTypeofLocalZone = GetZoneGridType(zoneID);

        if (gridTypeofLocalZone == STRUCTGRID)
        {
            string sgsmodel;
            GlobalDataBase::GetData("sgsmodel", &sgsmodel, PHSTRING, 1);
            solver->ComputeGradientCellCenter(grid);

            if (MIXGRID == systemGridType && isFVMOrFDM == FV_METHOD)
            {
                solver->ComputeGradientCellCenterOfruvwptOnlyForMixGrid(grid);
            }
         }
         else if(gridTypeofLocalZone == UNSTRUCTGRID)
         {
            solver->ComputeGradient(grid);
            if (MIXGRID == systemGridType && isFVMOrFDM == FV_METHOD)
            {
                solver->ComputelimiterOnlyForMixGrid(grid);
            }
        }
    }

    CommunicateGradientArray(iSolver, level);
}

void Controller::ComputeGamaAndTField(int iSolver, int level)
{
    if (NS_SOLVER != iSolver)
    {
        return;
    }

    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Zone *zone = GetZone(zoneID);
        if (!zone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        solver->ComputeGamaAndTemperature(grid);
    }
}

void Controller::ComputeViscousCoeff(int iSolver, int level)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();
    int viscousType = GlobalDataBase::GetIntParaFromDB("viscousType");
    if (viscousType == INVISCID)
    {
        return;
    }

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Zone *zone = GetZone(zoneID);
        if (!zone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        solver->ComputeViscousCoeff(grid);
    }
}

void Controller::ComputeBlendFunctionSST(int iSolver, int level)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();
    int viscousType = GlobalDataBase::GetIntParaFromDB("viscousType");
    int isFVMOrFDM = GlobalDataBase::GetIntParaFromDB("isFVMOrFDM");
    if (iSolver == NS_SOLVER || viscousType != TWO_EQU || 0 != level || isFVMOrFDM == FD_METHOD)
    {
        return;
    }

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Zone *zone = GetZone(zoneID);
        if (!zone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        solver->Crossing(grid);
        solver->Blending(grid);
    }
    CommunicateSpecificArray(iSolver, level, "blend", 0);
}

void Controller::SetBoundaryValueOfField(int iSolver, int level)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Zone *zone = GetZone(zoneID);
        if (!zone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        solver->Boundary(grid);
    }
}

void Controller::PreSolveforStructHighOrder(int iSolver, int level)
{
    //int mission_Comm(q_FD                                    ) =  5;
    //int mission_Comm(rtem                                    ) =  6;
    //int mission_Comm(rtem, dqdxyznokxi,dqdxyznoeta,dqdxyznocta) =  7;
    //int mission_Comm(rtem, gradUVWTCellCenterXYZ             ) =  8;
    //int mission_Comm(vist, gradqturb                         ) =  9;
    //int mission_Comm(vist, gradqturb, blending               ) = 10;

    int viscousType = GlobalDataBase::GetIntParaFromDB("viscousType");
    int iLES = GlobalDataBase::GetIntParaFromDB("iLES");

    if(iSolver == NS_SOLVER)
    {
        GetIndependentVariablesforStructHighOrder(iSolver, level);
        this->region->CommGridInfoOnlyforStructHighOrder(5);

        GetDependentVariablesforStructHighOrder(iSolver, level);
        if(viscousType == INVISCID)
        {
            //! Euler, no gradient communicated.
            this->region->CommGridInfoOnlyforStructHighOrder(6);
        }
        else
        {
            string gradient_method = GlobalDataBase::GetStrParaFromDB("structhighordergradient");
            if (gradient_method.substr(0,10) == "chain_rule")
            {
                this->region->CommGridInfoOnlyforStructHighOrder(7);
            }
            else
            {
                this->region->CommGridInfoOnlyforStructHighOrder(8);
            }
        }
    }
    else if(iSolver == TURB_SOLVER)
    {
        GetDependentVariablesforStructHighOrder(iSolver, level);
        if(viscousType == ONE_EQU)
        {
            this->region->CommGridInfoOnlyforStructHighOrder(9);
        }
        else if(viscousType == TWO_EQU)
        {
            this->region->CommGridInfoOnlyforStructHighOrder(10);
        }
        if(iLES == LES_SOLVER)
        {
            this->region->CommGridInfoOnlyforStructHighOrder(8);

            this->region->CommGridInfoOnlyforStructHighOrder(11);
        }
    }
}

void Controller::GetIndependentVariablesforStructHighOrder(int iSolver, int level)
{    
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Zone *zone = GetZone(zoneID);
        if (!zone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        solver->GetIndependentVariablesforStructHighOrder(grid);
    }
}

void Controller::GetDependentVariablesforStructHighOrder(int iSolver, int level)
{    
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Zone *zone = GetZone(zoneID);
        if (!zone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        solver->GetDependentVariablesforStructHighOrder(grid);
    }   
}

void Controller::RungeKuttaThirdTVD(int iSolver, int level)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

    FieldProxy **q0Proxy = new FieldProxy *[nLocalZones];

    CreatqProxyfor3rdRK(iSolver, level, q0Proxy);

    Stagefor3rdRK(iSolver, level, q0Proxy, 0);

    Stagefor3rdRK(iSolver, level, q0Proxy, 1);

    Stagefor3rdRK(iSolver, level, q0Proxy, 2);
    
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        delete q0Proxy[iZone];
    }
    delete []q0Proxy;
}

void Controller::CreatqProxyfor3rdRK(int iSolver, int level, FieldProxy **q0Proxy)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

#ifdef USE_OMP
    #pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Zone *zone = GetZone(zoneID);
        if (!zone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);
        
        RDouble4D &q = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));
        q0Proxy[iZone] = solver->CreateFieldProxy(grid);
        RDouble4D &q0 = q0Proxy[iZone]->GetField_STR();
        q0 = q;
    }
}

void Controller::Stagefor3rdRK(int iSolver, int level, FieldProxy **q0Proxy, int nstage)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();    
    
    if (nstage != 0)
    {
        //! already solved before timestep;
        PreSolveforStructHighOrder(iSolver, level);
    }

    LoadResiduals(iSolver, level);
    UpdateResiduals(iSolver, level);      

    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Zone *zone = GetZone(zoneID);
        if (!zone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);

        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);
        FieldProxy *resRKProxy = solver->CreateFieldProxy(grid);

        if (nstage == 0)
        {
            solver->RungeKuttaResidual(grid, resRKProxy, 1.0);
            solver->UpdateFlowField(grid, q0Proxy[iZone], resRKProxy);            
        }
        else if (nstage == 1)
        {
            FieldProxy *qtempProxy = solver->CreateFieldProxy(grid);
            
            RDouble4D &qtemp   = qtempProxy->GetField_STR();
            RDouble4D &q0   = q0Proxy[iZone]->GetField_STR();
            RDouble4D &q = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));

            qtemp = 0.75 * q0 + 0.25 * q;
            solver->RungeKuttaResidual(grid, resRKProxy, 0.25);
            solver->UpdateFlowField(grid, qtempProxy, resRKProxy);

            delete qtempProxy;
        }
        else if (nstage == 2)
        {
            FieldProxy *qtempProxy = solver->CreateFieldProxy(grid);
            
            RDouble4D &qtemp   = qtempProxy->GetField_STR();
            RDouble4D &q0   = q0Proxy[iZone]->GetField_STR();
            RDouble4D &q = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));

            qtemp = (q0 + 2.0 * q) / 3.0;
            solver->RungeKuttaResidual(grid, resRKProxy, 2.0 / 3.0);
            solver->UpdateFlowField(grid, qtempProxy, resRKProxy);

            delete qtempProxy;
        }
        delete resRKProxy;
    }
}

void Controller::UpdateIndependentField(int iSolver, int level)
{
    this->ComputeGamaAndTField(iSolver, level);

    this->PressureFactor(iSolver, level);

    this->ComputeViscousCoeff(iSolver, level);

    this->ComputeGradient(iSolver, level);

    this->ComputeBlendFunctionSST(iSolver, level);
}

void Controller::LUSGSInitialization(int iSolver, int level, int nSweep, FieldProxy **dqProxy, FieldProxy **LUplusDQ)
{
    int isFVMOrFDM = GlobalDataBase::GetIntParaFromDB("isFVMOrFDM");

    int nLocalZones = PHMPI::GetNumberofLocalZones();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = PHMPI::GetLocalZoneIDToGlobalZoneID(iZone);

        Zone *zone = GetZone(zoneID);
        if (!zone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        //! Memory allocating.
        dqProxy[iZone] = solver->CreateFieldProxy(grid);
        LUplusDQ[iZone]= solver->CreateFieldProxy(grid);
        
        int gridtype = grid->Type();
        if (isFVMOrFDM == FD_METHOD && gridtype == STRUCTGRID)
        {
            solver->LUSGSInitializationStructHighOrder(grid, dqProxy[iZone], LUplusDQ[iZone]);
        }
        else
        {
            //! The delta Flux is initialized by zero, which will be used in the forward sweep.
            solver->FillField(grid, LUplusDQ[iZone], zero);
            
            //* warning: it is important to research the effect of the initialization in future!
            //! e.g., using rhs to init when nsweep > 1.
            //! Initialization.
            if(nSweep == 1)
            {
                //! Single sweep, init by res.
                FieldProxy *rhsProxy = solver->GetFieldProxy(grid, solver->CastName("res")); 
                solver->FillField(grid, dqProxy[iZone], rhsProxy);
                delete rhsProxy;
            }
            else
            {
                //! Multi-sweep, init by zero.
                solver->FillField(grid, dqProxy[iZone], zero);
            }
        }
    }
}

void Controller::BLUSGSInitialization(int iSolver, int level, int nSweep, FieldProxy **dqProxy, FieldProxy **LUplusDQ)
{
    int nLocalZones = PHMPI::GetNumberofLocalZones();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = PHMPI::GetLocalZoneIDToGlobalZoneID(iZone);

        Zone *zone = GetZone(zoneID);
        if (!zone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        //! Memory allocating.
        dqProxy[iZone] = solver->CreateFieldProxy(grid);
        LUplusDQ[iZone]= solver->CreateFieldProxy(grid); 

        //! The delta Flux is initialized by zero, which will be used in the forward sweep.
        solver->FillField(grid, LUplusDQ[iZone], zero);

        //* warning: it is important to research the effect of the initialization in future!
        //! e.g., using rhs to init when nsweep > 1.
        //! Initialization.
        if(nSweep == 1)
        {
            //! Single sweep, init by res.
            FieldProxy *rhsProxy = solver->GetFieldProxy(grid, solver->CastName("res")); 
            solver->FillField(grid, dqProxy[iZone], rhsProxy);
            delete rhsProxy;
        }
        else
        {
            //! Multi-sweep, init by zero.
            solver->FillField(grid, dqProxy[iZone], zero);
        }
    }
}

void Controller::MLUSGSInitialization(int iSolver, int level, int nSweep, FieldProxy **dqProxy, FieldProxy **LUplusDQ)
{
    int nLocalZones = PHMPI::GetNumberofLocalZones();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = PHMPI::GetLocalZoneIDToGlobalZoneID(iZone);

        Zone *zone = GetZone(zoneID);
        if (!zone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        //! Memory allocating.
        dqProxy[iZone] = solver->CreateFieldProxy(grid);
        LUplusDQ[iZone]= solver->CreateFieldProxy(grid);

        //! The delta Flux is initialized by zero, which will be used in the forward sweep.
        solver->FillField(grid, LUplusDQ[iZone], zero);

        //* warning: it is important to research the effect of the initialization in future!
        //! e.g., using rhs to init when nsweep > 1.
        //! Initialization.
        if(nSweep == 1)
        {
            //! Single sweep, init by res.
            FieldProxy *rhsProxy = solver->GetFieldProxy(grid, solver->CastName("res")); 
            solver->FillField(grid, dqProxy[iZone], rhsProxy);
            delete rhsProxy;
        }
        else
        {
            //! Multi-sweep, init by zero.
            solver->FillField(grid, dqProxy[iZone], zero);
        }
        //*/    
    }
}

void Controller::LUSGS(int iSolver, int level)
{
    Diagonal(iSolver, level);

    int nLocalZones = PHMPI::GetNumberofLocalZones();
    int nLUSGSSweeps = GlobalDataBase::GetIntParaFromDB("nLUSGSSweeps");

    if (level != 0)
    {
        nLUSGSSweeps = 1;
    }

    RDouble LUSGSTolerance = GlobalDataBase::GetDoubleParaFromDB("LUSGSTolerance");
    int systemGridType = PHSPACE::GlobalDataBase::GetIntParaFromDB("sys_gridtype");
    int isFVMOrFDM = GlobalDataBase::GetIntParaFromDB("isFVMOrFDM");

    FieldProxy **dqProxy = new FieldProxy *[nLocalZones];

    //! L*Dq or U*Dq, which is equal to DeltaFlux between beighbors.
    FieldProxy **LUplusDQ= new FieldProxy *[nLocalZones];

    LUSGSInitialization(iSolver, level, nLUSGSSweeps, dqProxy, LUplusDQ);

    RDouble sweepNormalInit = 0.0, sweepEpsilon = 1.1, globleNormal = 0.0;
    for (int iSweep = 0; iSweep < nLUSGSSweeps; ++ iSweep)
    {
        RDouble sweepNormal = 0.0;

        if (systemGridType == MIXGRID || isFVMOrFDM == FD_METHOD)
        {
            if (iSolver == NS_SOLVER)
            {
                CommunicationDQ(iSolver, level, dqProxy);
            }

            ForwardSweep(iSolver, level, dqProxy, LUplusDQ, sweepNormal, false);

            BackwardSweep(iSolver, level, dqProxy, LUplusDQ, sweepNormal, false);
        }
        else
        {
            ForwardSweep(iSolver, level, dqProxy, LUplusDQ, sweepNormal, false);

            CommunicationDQ(iSolver, level, dqProxy);

            BackwardSweep(iSolver, level, dqProxy, LUplusDQ, sweepNormal, false);
        }

        if (nLUSGSSweeps > 1)
        {
            //! Get the global normal.
            PH_AllReduce(&sweepNormal, &globleNormal, 1, PH_SUM);
            sweepNormal = globleNormal;
        }

        //! Check the convergence.
        //! This criterion will lead to at least 2 sweep.
        //! It may be more considered, such as using the maximum normal of dq/dqOld.
        if(iSweep == 0)
        {
            sweepNormalInit = sweepNormal;
        }
        else
        {
            sweepEpsilon = sweepNormal / sweepNormalInit;
        }

        if (sweepEpsilon < LUSGSTolerance)
        {
            //! The equations have been converged.
            break;
        }
    }

    UpdateFlowField(iSolver, level, dqProxy);

    for (int iLocalZone = 0; iLocalZone < nLocalZones; ++ iLocalZone)
    {
        delete dqProxy[iLocalZone];
        delete LUplusDQ[iLocalZone];
    }
    delete []dqProxy;
    delete []LUplusDQ;
}

void Controller::DualTimeStepsLUSGS(int iSolver, int level)
{
    Diagonal(iSolver, level);

    int nLocalZones = PHMPI::GetNumberofLocalZones();
    int nLUSGSSweeps = GlobalDataBase::GetIntParaFromDB("nLUSGSSweeps");

    if (level != 0)
    {
        nLUSGSSweeps = 1;
    }

    RDouble LUSGSTolerance = GlobalDataBase::GetDoubleParaFromDB("LUSGSTolerance");
    int systemGridType = PHSPACE::GlobalDataBase::GetIntParaFromDB("sys_gridtype");
    int isFVMOrFDM = GlobalDataBase::GetIntParaFromDB("isFVMOrFDM");

    FieldProxy **dqProxy = new FieldProxy *[nLocalZones];
    //! L*Dq or U*Dq, which is equal to DeltaFlux between beighbors.
    FieldProxy **LUplusDQ= new FieldProxy *[nLocalZones];

    FieldProxy **dqProxy1 = new FieldProxy *[nLocalZones];
    //! L*Dq or U*Dq, which is equal to DeltaFlux between beighbors.
    FieldProxy **LUplusDQ1= new FieldProxy *[nLocalZones];

    LUSGSInitialization(iSolver, level, nLUSGSSweeps, dqProxy, LUplusDQ);
    LUSGSInitialization(iSolver, level, nLUSGSSweeps, dqProxy1, LUplusDQ1);

    //The basic iteration step.
    RDouble sweepNormalInit = 0.0, sweepEpsilon = 1.1, globleNormal = 0.0;
    for (int iSweep = 0; iSweep < nLUSGSSweeps; ++ iSweep)
    {
        RDouble sweepNormal = 0.0;

        if (systemGridType == MIXGRID || isFVMOrFDM == FD_METHOD)
        {
            if (iSolver == NS_SOLVER)
            {
                CommunicationDQ(iSolver, level, dqProxy);
            }

            ForwardSweep(iSolver, level, dqProxy, LUplusDQ, sweepNormal, false);

            BackwardSweep(iSolver, level, dqProxy, LUplusDQ, sweepNormal, false);
        }
        else
        {
            ForwardSweep(iSolver, level, dqProxy, LUplusDQ, sweepNormal, false);

            CommunicationDQ(iSolver, level, dqProxy);

            BackwardSweep(iSolver, level, dqProxy, LUplusDQ, sweepNormal, false);
        }

        if (nLUSGSSweeps > 1)
        {
            //! Get the global normal.
            PH_AllReduce(&sweepNormal, &globleNormal, 1, PH_SUM);
            sweepNormal = globleNormal;
        }

        //! Check the convergence.
        //! This criterion will lead to at least 2 sweep.
        //! It may be more considered, such as using the maximum normal of dq/dqOld.
        if(iSweep == 0)
        {
            sweepNormalInit = sweepNormal;
        }
        else
        {
            sweepEpsilon = sweepNormal / sweepNormalInit;
        }

        if (sweepEpsilon < LUSGSTolerance)
        {
            //! The equations have been converged.
            break;
        }
    }

    //The advancing time iteration step.
    for (int iSweep = 0; iSweep < nLUSGSSweeps; ++ iSweep)
    {
        RDouble sweepNormal = 0.0;

        if (systemGridType == MIXGRID || isFVMOrFDM == FD_METHOD)
        {
            if (iSolver == NS_SOLVER)
            {
                CommunicationDQ(iSolver, level, dqProxy1);
            }

            ForwardSweep(iSolver, level, dqProxy1, LUplusDQ1, sweepNormal, true);

            BackwardSweep(iSolver, level, dqProxy1, LUplusDQ1, sweepNormal, true);
        }
        else
        {
            ForwardSweep(iSolver, level, dqProxy1, LUplusDQ1, sweepNormal, true);

            CommunicationDQ(iSolver, level, dqProxy1);

            BackwardSweep(iSolver, level, dqProxy1, LUplusDQ1, sweepNormal, true);
        }

        if (nLUSGSSweeps > 1)
        {
            //! Get the global normal.
            PH_AllReduce(&sweepNormal, &globleNormal, 1, PH_SUM);
            sweepNormal = globleNormal;
        }

        //! Check the convergence.
        //! This criterion will lead to at least 2 sweep.
        //! It may be more considered, such as using the maximum normal of dq/dqOld.
        if(iSweep == 0)
        {
            sweepNormalInit = sweepNormal;
        }
        else
        {
            sweepEpsilon = sweepNormal / sweepNormalInit;
        }

        if (sweepEpsilon < LUSGSTolerance)
        {
            //! The equations have been converged.
            break;
        }
    }

    //Determine the CFL number of the next iteration step.
    DetermineCFLNumber(iSolver, level, dqProxy, dqProxy1);

    UpdateFlowField(iSolver, level, dqProxy);

    for (int iLocalZone = 0; iLocalZone < nLocalZones; ++ iLocalZone)
    {
        delete dqProxy[iLocalZone];
        delete LUplusDQ[iLocalZone];
        delete dqProxy1[iLocalZone];
        delete LUplusDQ1[iLocalZone];
    }
    delete []dqProxy;
    delete []LUplusDQ;
    delete []dqProxy1;
    delete []LUplusDQ1;
}

void Controller::ForwardSweep(int iSolver, int level, FieldProxy **dqProxy, FieldProxy **LUplusDQ, RDouble &sweepNormal, bool iAdvanceStep)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();
    RDouble *zoneNormal = new RDouble [nLocalZones]();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Zone *zone = GetZone(zoneID);
        if (!zone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID,level);
        CFDSolver *solver = GetCFDSolver(zoneID,iSolver);       

        zoneNormal[iZone] = 0.0;
        solver->SolveLUSGSForward(grid, dqProxy[iZone], LUplusDQ[iZone], zoneNormal[iZone], iAdvanceStep);
    }

    //! Sum the zone normal, which is to void OMP data competition.
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        sweepNormal += zoneNormal[iZone];
    }
    delete []zoneNormal;    zoneNormal = NULL;
}

void Controller::BackwardSweep(int iSolver, int level, FieldProxy **dqProxy, FieldProxy **LUplusDQ, RDouble &sweepNormal, bool iAdvanceStep)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();
    RDouble *zoneNormal = new RDouble [nLocalZones]();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Zone *zone = GetZone(zoneID);
        if (!zone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID,level);
        CFDSolver *solver = GetCFDSolver(zoneID,iSolver);

        zoneNormal[iZone] = 0.0;
        solver->SolveLUSGSBackward(grid, dqProxy[iZone], LUplusDQ[iZone], zoneNormal[iZone], iAdvanceStep);
    }

    //! Sum the zone normal, which is to void OMP data competition.
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        sweepNormal += zoneNormal[iZone];
    }
    delete []zoneNormal;    zoneNormal = NULL;
}

void Controller::UpdateFlowField(int iSolver, int level, FieldProxy **dqProxy)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();
    int isFVMOrFDM = GlobalDataBase::GetIntParaFromDB("isFVMOrFDM");

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Zone *zone = GetZone(zoneID);
        if (!zone)
        {
            continue;
        }
        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        int gridtype = grid->Type();
        if (iSolver == NS_SOLVER && isFVMOrFDM == FD_METHOD && gridtype == STRUCTGRID)
        {
            FieldProxy *qtemp_Proxy = solver->CreateFieldProxy(grid);
            FieldProxy *qfield      = solver->GetFieldProxy(grid, solver->CastName("q"));
            solver->FillField(grid, qtemp_Proxy, qfield);

            solver->UpdateFlowField(grid, qtemp_Proxy, dqProxy[iZone]);
            delete qtemp_Proxy;
        }
        else
        {
            FieldProxy *qProxy = solver->GetFieldProxy(grid, solver->CastName("q"));
            solver->UpdateFlowField(grid, qProxy, dqProxy[iZone]);
            delete qProxy;
        }
    }
}

void Controller::DetermineCFLNumber(int iSolver, int level, FieldProxy **dqProxy, FieldProxy **dqProxy1)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Zone *zone = GetZone(zoneID);
        if (!zone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        solver->DetermineCFLNumber(grid, dqProxy[iZone], dqProxy1[iZone]);
    }
}

void Controller::CommunicateGradientandLimitOnlyForMixGrid()
{
    int currentProcessor = PHMPI::GetCurrentProcessorID();
    int nZones = PHMPI::GetNumberofGlobalZones();
    vector <PH_Request > requestContainer;
    vector <vector <DataContainer *> > receivedDataBuffer;
    vector <vector <DataContainer *> > sendDataBuffer;
    receivedDataBuffer.resize(nZones);
    sendDataBuffer.resize(nZones);

    //! Step 0: Compressing data firstly.
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int sendProcessor = PHMPI::GetZoneProcessorID(iZone);

        ZoneNeighbor *globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
        int numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();
        if (currentProcessor == sendProcessor)
        {
            sendDataBuffer[iZone].resize(numberOfNeighbor);
        }

        for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
        {
            int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);

            int receiveProcessor = PHMPI::GetZoneProcessorID(neighborZone);

            if (currentProcessor != sendProcessor && currentProcessor != receiveProcessor)
            {
                continue;
            }

            //! Allocate and compress the buffers for sending.
            DataContainer *sendBuffer = 0;
            if (currentProcessor == sendProcessor)
            {
                sendBuffer = new DataContainer;
                sendDataBuffer[iZone][iNeighbor] = sendBuffer;

                //CFDSolver *solver = GetCFDSolver(iZone, iSolver);
                CFDSolver *solver = GetCFDSolver(iZone, 0);
                //! Compress the send information into the actkey.
                Grid *grid = PHSPACE::GetGrid(iZone, 0);

                solver->CompressGradientandLimitToInterfaceOnlyForMixGrid(sendBuffer, grid, neighborZone);
            }
        }
    }

    //! Step 1: Communication.
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        ZoneNeighbor *globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
        int numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();

        //! Communicating.
        for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
        {
            int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);

            int sendProcessor = PHMPI::GetZoneProcessorID(iZone);
            int receiveProcessor = PHMPI::GetZoneProcessorID(neighborZone);

            if (currentProcessor != sendProcessor && currentProcessor != receiveProcessor)
            {
                continue;
            }

            //! Allocate the buffers for sending and receiving.
            DataContainer *receivedBuffer = 0;
            DataContainer *sendBuffer = 0;
            if (currentProcessor == receiveProcessor)
            {
                receivedBuffer = new DataContainer;
                receivedDataBuffer[iZone].push_back(receivedBuffer);
            }
            if (currentProcessor == sendProcessor)
            {
                sendBuffer = sendDataBuffer[iZone][iNeighbor];
            }

            int tag = iZone;

            //if (sendProcessor == receiveProcessor)
            //{
            //    continue;
            //}

            //! Communication.
            if (currentProcessor == sendProcessor)
            {
                if (sendProcessor != receiveProcessor)
                {
                    streamsize nlen = sendBuffer->size();

                    //! Send the data to neighbors.
                    send(sendBuffer, receiveProcessor, nlen, requestContainer, tag);
                }
            }
            if (currentProcessor == receiveProcessor)
            {
                //! Get the neighbor order: which order that 'iZone' on neighbors of neighborZone.
                int neighborOrer = -1;
                ZoneNeighbor *globalNeighborZonesTemp = zoneConnectivity->GetZoneNeighbor(neighborZone);
                int numberOfNeighborTemp = globalNeighborZonesTemp->GetNumberOfNeighbor();

                for (int neighborID = 0; neighborID < numberOfNeighborTemp; ++ neighborID)
                {
                    int neighborZoneTemp = globalNeighborZonesTemp->GetZoneIndexOfNeighbor(neighborID);
                    if (neighborZoneTemp == iZone)
                    {
                        neighborOrer = neighborID;
                        break;
                    }
                }

                ASSERT(neighborOrer != -1);

                if (sendProcessor != receiveProcessor)
                {
                    //! The data length of the received data is same to the length that send to the neighbor.
                    CharVecSizeType nlen = sendDataBuffer[neighborZone][neighborOrer]->size();

                    //! Receive data from neighbors.
                    receive(receivedBuffer, sendProcessor, nlen, requestContainer, tag);
                }
                else
                {
                    SWAP(receivedDataBuffer[iZone].back(), sendDataBuffer[iZone][iNeighbor]);
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
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        ZoneNeighbor *globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
        int numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();

        int count = 0;
        for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
        {
            int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);

            int receiveProcessor = PHMPI::GetZoneProcessorID(neighborZone);

            if (currentProcessor == receiveProcessor)
            {
                DataContainer *receiveData = receivedDataBuffer[iZone][count];

                //! Because get into here when currentProcessor==receiveProcessor,
                //! so, use zoneGridSolver of neighborZone, this zone is on the current processor!

                //! Decompress the interface data from data-container.
                //CFDSolver *solver = GetCFDSolver(neighborZone, iSolver);
                CFDSolver *solver = GetCFDSolver(neighborZone, 0);

                Grid *grid = PHSPACE::GetGrid(neighborZone, 0);

                solver->DecompressGradientandLimitToInterfaceOnlyForMixGrid(receiveData, grid, iZone);

                ++count;
            }
        }
    }

    //! Step4: Free the buffers.
    for (int iDim = 0; iDim < receivedDataBuffer.size(); ++ iDim)
    {
        for (int jDim = 0; jDim < receivedDataBuffer[iDim].size(); ++ jDim)
        {
            delete receivedDataBuffer[iDim][jDim];
        }
    }
    receivedDataBuffer.clear();
    for (int iDim = 0; iDim < sendDataBuffer.size(); ++ iDim)
    {
        for (int jDim = 0; jDim < sendDataBuffer[iDim].size(); ++ jDim)
        {
            delete sendDataBuffer[iDim][jDim];
        }
    }
    sendDataBuffer.clear();
    requestContainer.clear();
}

void Controller::CommunicateQlQrOnlyforMixGrid()
{
    int systemGridType = PHSPACE::GlobalDataBase::GetIntParaFromDB("sys_gridtype");
    if(systemGridType != MIXGRID)
    {
        return;
    }
    int currentProcessor = PHMPI::GetCurrentProcessorID();
    int nZones = PHMPI::GetNumberofGlobalZones();
    vector <PH_Request > requestContainer;
    vector <vector <DataContainer *> > receivedDataBuffer;
    vector <vector <DataContainer *> > sendDataBuffer;
    receivedDataBuffer.resize(nZones);
    sendDataBuffer.resize(nZones);

    //! Step 0: Compressing data firstly.
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int gridtypeofiZone = PHMPI::GetZoneGridType(iZone);

        int sendProcessor = PHMPI::GetZoneProcessorID(iZone);

        ZoneNeighbor *globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
        int numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();
        if (currentProcessor == sendProcessor)
        {
            sendDataBuffer[iZone].resize(numberOfNeighbor);
        }

        for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
        {
            int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);

            int gridtypeofneighborZone = PHMPI::GetZoneGridType(neighborZone);

            int receiveProcessor = PHMPI::GetZoneProcessorID(neighborZone);

            if (currentProcessor != sendProcessor && currentProcessor != receiveProcessor)
            {
                continue;
            }

            //! Allocate and compress the buffers for sending.
            DataContainer *sendBuffer = 0;
            if (currentProcessor == sendProcessor)
            {
                sendBuffer = new DataContainer;
                sendDataBuffer[iZone][iNeighbor] = sendBuffer;

                //CFDSolver *solver = GetCFDSolver(iZone, iSolver);
                CFDSolver *solver = GetCFDSolver(iZone, 0);
                //! Compress the send information into the actkey.
                Grid *grid = PHSPACE::GetGrid(iZone, 0);

                InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
                if (!interfaceInformation)
                {
                }
                else
                {
                    if(gridtypeofiZone == gridtypeofneighborZone)
                    {
                        RDouble uploadData = 0.0;
                        sendBuffer->moveto_begin();
                        PHWrite(sendBuffer, uploadData);
                    }
                    else
                    {
                        solver->CompressQlQrToInterfaceOnlyForMixGrid(sendBuffer, grid, neighborZone);
                    }
                }
            }
        }
    }

    //! Step 1: Communication.
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        ZoneNeighbor *globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
        int numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();

        //! Communicating.
        for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
        {
            int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);

            int sendProcessor = PHMPI::GetZoneProcessorID(iZone);
            int receiveProcessor = PHMPI::GetZoneProcessorID(neighborZone);

            if (currentProcessor != sendProcessor && currentProcessor != receiveProcessor)
            {
                continue;
            }

            //! Allocate the buffers for sending and receiving.
            DataContainer *receivedBuffer = 0;
            DataContainer *sendBuffer = 0;
            if (currentProcessor == receiveProcessor)
            {
                receivedBuffer = new DataContainer;
                receivedDataBuffer[iZone].push_back(receivedBuffer);
            }
            if (currentProcessor == sendProcessor)
            {
                sendBuffer = sendDataBuffer[iZone][iNeighbor];
            }

            int tag = iZone;

            //if (sendProcessor == receiveProcessor)
            //{
            //    continue;
            //}

            //! Communication.
            if (currentProcessor == sendProcessor)
            {
                if (sendProcessor != receiveProcessor)
                {
                    streamsize nlen = sendBuffer->size();

                    //! Send the data to neighbors.
                    send(sendBuffer, receiveProcessor, nlen, requestContainer, tag);
                }
            }
            if (currentProcessor == receiveProcessor)
            {
                //! Get the neighbor order: which order that 'iZone' on neighbors of neighborZone.
                int neighborOrer = -1;
                ZoneNeighbor *globalNeighborZonesTemp = zoneConnectivity->GetZoneNeighbor(neighborZone);
                int numberOfNeighborTemp = globalNeighborZonesTemp->GetNumberOfNeighbor();

                for (int neighborID = 0; neighborID < numberOfNeighborTemp; ++ neighborID)
                {
                    int neighborZoneTemp = globalNeighborZonesTemp->GetZoneIndexOfNeighbor(neighborID);
                    if (neighborZoneTemp == iZone)
                    {
                        neighborOrer = neighborID;
                        break;
                    }
                }

                ASSERT(neighborOrer != -1);

                if (sendProcessor != receiveProcessor)
                {
                    //! The data length of the received data is same to the length that send to the neighbor.
                    CharVecSizeType nlen = sendDataBuffer[neighborZone][neighborOrer]->size();

                    //! Receive data from neighbors.
                    receive(receivedBuffer, sendProcessor, nlen, requestContainer, tag);
                }
                else
                {
                    SWAP(receivedDataBuffer[iZone].back(), sendDataBuffer[iZone][iNeighbor]);
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
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int gridtypeofiZone = PHMPI::GetZoneGridType(iZone);

        ZoneNeighbor *globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
        int numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();

        int count = 0;
        for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
        {
            int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);

            int gridtypeofneighborZone = PHMPI::GetZoneGridType(neighborZone);

            int receiveProcessor = PHMPI::GetZoneProcessorID(neighborZone);

            if (currentProcessor == receiveProcessor)
            {
                DataContainer *receiveData = receivedDataBuffer[iZone][count];

                //! Because get into here when currentProcessor==receiveProcessor,
                //! so, use zoneGridSolver of neighborZone, this zone is on the current processor!

                //! Decompress the interface data from data-container.
                //CFDSolver *solver = GetCFDSolver(neighborZone, iSolver);
                CFDSolver *solver = GetCFDSolver(neighborZone, 0);

                Grid *grid = PHSPACE::GetGrid(neighborZone, 0);

                InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
                if (!interfaceInformation)
                {
                    //! can not use return here, because ++count in the after
                }
                else
                {
                    if(gridtypeofiZone == gridtypeofneighborZone)
                    {
                        RDouble downdData = 0.0;
                        receiveData->moveto_begin();
                        PHRead(receiveData, downdData);
                    }
                    else
                    {
                        solver->DecompressQlQrToInterfaceOnlyForMixGrid(receiveData, grid, iZone);
                    }
                }

                ++count;
            }
        }
    }

    //! Step4: Free the buffers.
    for (int iDim = 0; iDim < receivedDataBuffer.size(); ++ iDim)
    {
        for (int jDim = 0; jDim < receivedDataBuffer[iDim].size(); ++ jDim)
        {
            delete receivedDataBuffer[iDim][jDim];
        }
    }
    receivedDataBuffer.clear();
    for (int iDim = 0; iDim < sendDataBuffer.size(); ++ iDim)
    {
        for (int jDim = 0; jDim < sendDataBuffer[iDim].size(); ++ jDim)
        {
            delete sendDataBuffer[iDim][jDim];
        }
    }
    sendDataBuffer.clear();
    requestContainer.clear();
}

void Controller::CommunicateFacefluxOnlyforMixGrid()
{
    int systemGridType = PHSPACE::GlobalDataBase::GetIntParaFromDB("sys_gridtype");
    if(systemGridType != MIXGRID)
    {
        return;
    }

    int currentProcessor = PHMPI::GetCurrentProcessorID();
    int nZones = PHMPI::GetNumberofGlobalZones();
    vector <PH_Request > requestContainer;
    vector <vector <DataContainer *> > receivedDataBuffer;
    vector <vector <DataContainer *> > sendDataBuffer;
    receivedDataBuffer.resize(nZones);
    sendDataBuffer.resize(nZones);

    //! Step 0: Compressing data firstly.
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int sendProcessor = PHMPI::GetZoneProcessorID(iZone);

        ZoneNeighbor *globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
        int numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();
        if (currentProcessor == sendProcessor)
        {
            sendDataBuffer[iZone].resize(numberOfNeighbor);
        }

        for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
        {
            int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);

            int receiveProcessor = PHMPI::GetZoneProcessorID(neighborZone);

            if (currentProcessor != sendProcessor && currentProcessor != receiveProcessor)
            {
                continue;
            }

            //! Allocate and compress the buffers for sending.
            DataContainer *sendBuffer = 0;
            if (currentProcessor == sendProcessor)
            {
                sendBuffer = new DataContainer;
                sendDataBuffer[iZone][iNeighbor] = sendBuffer;

                //CFDSolver *solver = GetCFDSolver(iZone, iSolver);
                CFDSolver *solver = GetCFDSolver(iZone, 0);
                //! Compress the send information into the actkey.
                Grid *grid = PHSPACE::GetGrid(iZone, 0);

                solver->CompressFacefluxToInterfaceOnlyForMixGrid(sendBuffer, grid, neighborZone);
            }
        }
    }

    //! Step 1: Communication.
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        ZoneNeighbor *globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
        int numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();

        //! Communicating.
        for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
        {
            int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);

            int sendProcessor = PHMPI::GetZoneProcessorID(iZone);
            int receiveProcessor = PHMPI::GetZoneProcessorID(neighborZone);

            if (currentProcessor != sendProcessor && currentProcessor != receiveProcessor)
            {
                continue;
            }

            //! Allocate the buffers for sending and receiving.
            DataContainer *receivedBuffer = 0;
            DataContainer *sendBuffer = 0;
            if (currentProcessor == receiveProcessor)
            {
                receivedBuffer = new DataContainer;
                receivedDataBuffer[iZone].push_back(receivedBuffer);
            }
            if (currentProcessor == sendProcessor)
            {
                sendBuffer = sendDataBuffer[iZone][iNeighbor];
            }

            int tag = iZone;

            //! Communication.
            if (currentProcessor == sendProcessor)
            {
                if (sendProcessor != receiveProcessor)
                {
                    streamsize nlen = sendBuffer->size();

                    //! Send the data to neighbors.
                    send(sendBuffer, receiveProcessor, nlen, requestContainer, tag);
                }
            }
            if (currentProcessor == receiveProcessor)
            {
                //! Get the neighbor order: which order that 'iZone' on neighbors of neighborZone.
                int neighborOrer = -1;
                ZoneNeighbor *globalNeighborZonesTemp = zoneConnectivity->GetZoneNeighbor(neighborZone);
                int numberOfNeighborTemp = globalNeighborZonesTemp->GetNumberOfNeighbor();

                for (int neighborID = 0; neighborID < numberOfNeighborTemp; ++ neighborID)
                {
                    int neighborZoneTemp = globalNeighborZonesTemp->GetZoneIndexOfNeighbor(neighborID);
                    if (neighborZoneTemp == iZone)
                    {
                        neighborOrer = neighborID;
                        break;
                    }
                }

                ASSERT(neighborOrer != -1);

                if (sendProcessor != receiveProcessor)
                {
                    //! The data length of the received data is same to the length that send to the neighbor.
                    CharVecSizeType nlen = sendDataBuffer[neighborZone][neighborOrer]->size();

                    //! Receive data from neighbors.
                    receive(receivedBuffer, sendProcessor, nlen, requestContainer, tag);
                }
                else
                {
                    SWAP(receivedDataBuffer[iZone].back(), sendDataBuffer[iZone][iNeighbor]);
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
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        ZoneNeighbor *globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
        int numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();

        int count = 0;
        for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
        {
            int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);

            int receiveProcessor = PHMPI::GetZoneProcessorID(neighborZone);

            if (currentProcessor == receiveProcessor)
            {
                DataContainer *receiveData = receivedDataBuffer[iZone][count];

                //! Because get into here when currentProcessor==receiveProcessor,
                //! so, use zoneGridSolver of neighborZone, this zone is on the current processor!

                //! Decompress the interface data from data-container.
                //CFDSolver *solver = GetCFDSolver(neighborZone, iSolver);
                CFDSolver *solver = GetCFDSolver(neighborZone, 0);

                Grid *grid = PHSPACE::GetGrid(neighborZone, 0);

                solver->DecompressFacefluxToInterfaceOnlyForMixGrid(receiveData, grid, iZone);

                ++count;
            }
        }
    }

    //! Step4: Free the buffers.
    for (int iDim = 0; iDim < receivedDataBuffer.size(); ++ iDim)
    {
        for (int jDim = 0; jDim < receivedDataBuffer[iDim].size(); ++ jDim)
        {
            delete receivedDataBuffer[iDim][jDim];
        }
    }
    receivedDataBuffer.clear();
    for (int iDim = 0; iDim < sendDataBuffer.size(); ++ iDim)
    {
        for (int jDim = 0; jDim < sendDataBuffer[iDim].size(); ++ jDim)
        {
            delete sendDataBuffer[iDim][jDim];
        }
    }
    sendDataBuffer.clear();
    requestContainer.clear();
}

void Controller::CommunicationDQ(int iSolver, int level, FieldProxy **dqProxy)
{
    //if (iSolver != NS_SOLVER)
    //{
    //    //! Only consider the communication of NS-Solver.
    //    //return;
    //}

    int swapDq = GlobalDataBase::GetIntParaFromDB("swapDq");
    if (swapDq == 0)
    {
        return;
    }

    int zoneID = PHMPI::GetLocalZoneIDToGlobalZoneID(0);
    PHSolver *solver = GetCFDSolver(zoneID, iSolver);
    int nEquation  = solver->GetNumberOfEquations();

    int currentProcessor = PHMPI::GetCurrentProcessorID();
    int nZones = PHMPI::GetNumberofGlobalZones();
    vector <PH_Request > requestContainer;
    vector <vector <DataContainer *> > receivedDataBuffer;
    vector <vector <DataContainer *> > sendDataBuffer;
    receivedDataBuffer.resize(nZones);
    sendDataBuffer.resize(nZones);

    //! Step 0: Compressing data firstly.
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int sendProcessor = PHMPI::GetZoneProcessorID(iZone);

        ZoneNeighbor *globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
        int numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();
        if (currentProcessor == sendProcessor)
        {
            sendDataBuffer[iZone].resize(numberOfNeighbor);
        }

        for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
        {
            int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);

            int receiveProcessor = PHMPI::GetZoneProcessorID(neighborZone);

            if (currentProcessor != sendProcessor && currentProcessor != receiveProcessor)
            {
                continue;
            }

            //! Allocate and compress the buffers for sending.
            DataContainer *sendBuffer = 0;
            if (currentProcessor == sendProcessor)
            {
                sendBuffer = new DataContainer;
                sendDataBuffer[iZone][iNeighbor] = sendBuffer;

                int localZoneID = region->GetProcessGlobalZoneIndexToLocalZoneIndex(iZone);
                CFDSolver *solver = GetCFDSolver(iZone, iSolver);

                //! Compress the send information into the actkey.
                Grid *grid = PHSPACE::GetGrid(iZone, level);
                solver->CompressDQ(sendBuffer, dqProxy[localZoneID], grid, iZone, neighborZone, nEquation);
            }
        }
    }

    //! Step 1: Communication.
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        ZoneNeighbor *globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
        int numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();

        //! Communicating.
        for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
        {
            int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);

            int sendProcessor    = PHMPI::GetZoneProcessorID(iZone);
            int receiveProcessor = PHMPI::GetZoneProcessorID(neighborZone);

            if (currentProcessor != sendProcessor && currentProcessor != receiveProcessor)
            {
                continue;
            }

            //! Allocate the buffers for sending and receiving.
            DataContainer *receivedBuffer = 0;
            DataContainer *sendBuffer = 0;
            if (currentProcessor == receiveProcessor)
            {
                receivedBuffer = new DataContainer;
                receivedDataBuffer[iZone].push_back(receivedBuffer);
            }
            if (currentProcessor == sendProcessor)
            {
                sendBuffer = sendDataBuffer[iZone][iNeighbor];
            }

            int tag = iZone;

            //! Communication.
            if (currentProcessor == sendProcessor)
            {
                if (sendProcessor != receiveProcessor)
                {
                    streamsize nlen = sendBuffer->size();

                    // Send the data to neighbors.
                    send(sendBuffer, receiveProcessor, nlen, requestContainer, tag);
                }
            }
            if (currentProcessor == receiveProcessor)
            {
                //! Get the neighbor order: which order that 'iZone' on neighbors of neighborZone.
                int neighborOrer = -1;
                ZoneNeighbor *globalNeighborZonesTemp = zoneConnectivity->GetZoneNeighbor(neighborZone);
                int numberOfNeighborTemp = globalNeighborZonesTemp->GetNumberOfNeighbor();

                for (int neighborID = 0; neighborID < numberOfNeighborTemp; ++ neighborID)
                {
                    int neighborZoneTemp = globalNeighborZonesTemp->GetZoneIndexOfNeighbor(neighborID);
                    if (neighborZoneTemp == iZone)
                    {
                        neighborOrer = neighborID;
                        break;
                    }
                }

                ASSERT(neighborOrer != -1);

                if (sendProcessor != receiveProcessor)
                {
                    //! The data length of the received data is same to the length that send to the neighbor.
                    CharVecSizeType nlen = sendDataBuffer[neighborZone][neighborOrer]->size();

                    //! Receive data from neighbors.
                    receive(receivedBuffer, sendProcessor, nlen, requestContainer, tag);
                }
                else
                {
                    SWAP(receivedDataBuffer[iZone].back(), sendDataBuffer[iZone][iNeighbor]);
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
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        ZoneNeighbor *globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
        int numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();

        int count = 0;
        for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
        {
            int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);

            int receiveProcessor = PHMPI::GetZoneProcessorID(neighborZone);

            if (currentProcessor == receiveProcessor)
            {
                DataContainer *receiveData = receivedDataBuffer[iZone][count];

                //! Because get into here when currentProcessor==receiveProcessor,
                //! so, use zoneGridSolver of neighborZone, this zone is on the current processor!

                //! Decompress the interface data from data-container.
                int localZoneID = region->GetProcessGlobalZoneIndexToLocalZoneIndex(neighborZone);
                CFDSolver *solver = GetCFDSolver(neighborZone, iSolver);

                Grid *grid = PHSPACE::GetGrid(neighborZone, level);
                solver->DecompressDQ(receiveData, dqProxy[localZoneID], grid, neighborZone, iZone, nEquation);

                ++ count;
            }
        }
    }

    //! Step4: Free the buffers.
    for (int iDim = 0; iDim < receivedDataBuffer.size(); ++ iDim)
    {
        for (int jDim = 0; jDim < receivedDataBuffer[iDim].size(); ++ jDim)
        {
            delete receivedDataBuffer[iDim][jDim];
        }
    }
    receivedDataBuffer.clear();
    for (int iDim = 0; iDim < sendDataBuffer.size(); ++ iDim)
    {
        for (int jDim = 0; jDim < sendDataBuffer[iDim].size(); ++ jDim)
        {
            delete sendDataBuffer[iDim][jDim];
        }
    }
    sendDataBuffer.clear();
    requestContainer.clear();
}

void Controller::CommunicateSpecificArray(int iSolver, int level, string arrayName, int nEquation)
{
    int currentProcessor = PHMPI::GetCurrentProcessorID();
    int nZones = PHMPI::GetNumberofGlobalZones();
    vector <PH_Request > requestContainer;
    vector <vector <DataContainer *> > receivedDataBuffer;
    vector <vector <DataContainer *> > sendDataBuffer;
    receivedDataBuffer.resize(nZones);
    sendDataBuffer.resize(nZones);

    //! Step 0: Compressing data firstly.
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int sendProcessor = PHMPI::GetZoneProcessorID(iZone);

        ZoneNeighbor *globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
        int numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();
        if (currentProcessor == sendProcessor)
        {
            sendDataBuffer[iZone].resize(numberOfNeighbor);
        }

        for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
        {
            int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);

            int receiveProcessor = PHMPI::GetZoneProcessorID(neighborZone);

            if (currentProcessor != sendProcessor && currentProcessor != receiveProcessor)
            {
                continue;
            }

            //! Allocate and compress the buffers for sending.
            DataContainer *sendBuffer = 0;
            if (currentProcessor == sendProcessor)
            {
                sendBuffer = new DataContainer;
                sendDataBuffer[iZone][iNeighbor] = sendBuffer;

                //CFDSolver *solver = GetCFDSolver(iZone, iSolver);
                CFDSolver *solver = GetCFDSolver(iZone, iSolver);
                //! Compress the send information into the actkey.
                Grid *grid = PHSPACE::GetGrid(iZone, level);

                solver->CompressSpecificArrayToInterface(sendBuffer, arrayName, grid, neighborZone, nEquation);
            }
        }
    }

    //! Step 1: Communication.
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        ZoneNeighbor *globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
        int numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();

        //! Communicating.
        for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
        {
            int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);

            int sendProcessor = PHMPI::GetZoneProcessorID(iZone);
            int receiveProcessor = PHMPI::GetZoneProcessorID(neighborZone);

            if (currentProcessor != sendProcessor && currentProcessor != receiveProcessor)
            {
                continue;
            }

            //! Allocate the buffers for sending and receiving.
            DataContainer *receivedBuffer = 0;
            DataContainer *sendBuffer = 0;
            if (currentProcessor == receiveProcessor)
            {
                receivedBuffer = new DataContainer;
                receivedDataBuffer[iZone].push_back(receivedBuffer);
            }
            if (currentProcessor == sendProcessor)
            {
                sendBuffer = sendDataBuffer[iZone][iNeighbor];
            }

            int tag = iZone;

            //! Communication.
            if (currentProcessor == sendProcessor)
            {
                if (sendProcessor != receiveProcessor)
                {
                    streamsize nlen = sendBuffer->size();

                    //! Send the data to neighbors.
                    send(sendBuffer, receiveProcessor, nlen, requestContainer, tag);
                }
            }
            if (currentProcessor == receiveProcessor)
            {
                //! Get the neighbor order: which order that 'iZone' on neighbors of neighborZone.
                int neighborOrer = -1;
                ZoneNeighbor *globalNeighborZonesTemp = zoneConnectivity->GetZoneNeighbor(neighborZone);
                int numberOfNeighborTemp = globalNeighborZonesTemp->GetNumberOfNeighbor();

                for (int neighborID = 0; neighborID < numberOfNeighborTemp; ++ neighborID)
                {
                    int neighborZoneTemp = globalNeighborZonesTemp->GetZoneIndexOfNeighbor(neighborID);
                    if (neighborZoneTemp == iZone)
                    {
                        neighborOrer = neighborID;
                        break;
                    }
                }

                ASSERT(neighborOrer != -1);

                if (sendProcessor != receiveProcessor)
                {
                    //! The data length of the received data is same to the length that send to the neighbor.
                    CharVecSizeType nlen = sendDataBuffer[neighborZone][neighborOrer]->size();

                    //! Receive data from neighbors.
                    receive(receivedBuffer, sendProcessor, nlen, requestContainer, tag);
                }
                else
                {
                    SWAP(receivedDataBuffer[iZone].back(), sendDataBuffer[iZone][iNeighbor]);
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
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        ZoneNeighbor *globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
        int numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();

        int count = 0;
        for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
        {
            int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);

            int receiveProcessor = PHMPI::GetZoneProcessorID(neighborZone);

            if (currentProcessor == receiveProcessor)
            {
                DataContainer *receiveData = receivedDataBuffer[iZone][count];

                //! Because get into here when currentProcessor==receiveProcessor,
                //! so, use zoneGridSolver of neighborZone, this zone is on the current processor!

                //! Decompress the interface data from data-container.
                //CFDSolver *solver = GetCFDSolver(neighborZone, iSolver);
                CFDSolver *solver = GetCFDSolver(neighborZone, iSolver);

                Grid *grid = PHSPACE::GetGrid(neighborZone, level);
                solver->DecompressArrayFromInterface(receiveData, arrayName, grid, iZone, nEquation);

                ++count;
            }
        }
    }

    //! Step4: Free the buffers.
    for (int iDim = 0; iDim < receivedDataBuffer.size(); ++ iDim)
    {
        for (int jDim = 0; jDim < receivedDataBuffer[iDim].size(); ++ jDim)
        {
            delete receivedDataBuffer[iDim][jDim];
        }
    }
    receivedDataBuffer.clear();
    for (int iDim = 0; iDim < sendDataBuffer.size(); ++ iDim)
    {
        for (int jDim = 0; jDim < sendDataBuffer[iDim].size(); ++ jDim)
        {
            delete sendDataBuffer[iDim][jDim];
        }
    }
    sendDataBuffer.clear();
    requestContainer.clear();
}

void Controller::CommunicateGenericArray(int iSolver, int level)
{
    using namespace PHMPI;
    int zoneID = GetLocalZoneIDToGlobalZoneID(0);
    CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

    if (iSolver == NS_SOLVER)
    {
        int nEquation = solver->GetNumberOfEquations();
        CommunicateSpecificArray(iSolver, level, "q", nEquation);

        int nTemperatureModel = GlobalDataBase::GetIntParaFromDB("ntmodel");
        CommunicateSpecificArray(iSolver, level, "t", nTemperatureModel);
    }
    else if (iSolver == TURB_SOLVER)
    {
        int nEquation = solver->GetNumberOfEquations();
        CommunicateSpecificArray(iSolver, level, "q_turb", nEquation);
    }
}

void Controller::CommunicateGradientArray(int iSolver, int level)
{
    using namespace PHMPI;
    int zoneID = GetLocalZoneIDToGlobalZoneID(0);
    CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

    int systemGridType = PHSPACE::GlobalDataBase::GetIntParaFromDB("sys_gridtype");
    int isFVMOrFDM = PHSPACE::GlobalDataBase::GetIntParaFromDB("isFVMOrFDM");
    int periodicType = PHSPACE::GlobalDataBase::GetIntParaFromDB("periodicType");
    int iLES = PHSPACE::GlobalDataBase::GetIntParaFromDB("iLES");

    //! gradient has been communicated in another place for isFVMOrFDM == FDM.
    if (systemGridType == STRUCTGRID && isFVMOrFDM == FV_METHOD)
    {
        int nZones = PHMPI::GetNumberofGlobalZones();
        int currentProcessor = PHMPI::GetCurrentProcessorID();

        if (iSolver == NS_SOLVER)
        {
            //! NSSolver, GradientCommunication for struct high order has been done at another place, with another function.
            int nChemical = GlobalDataBase::GetIntParaFromDB("nchem");
            int nLaminar = GlobalDataBase::GetIntParaFromDB("nl");
            int nNSEquation = GlobalDataBase::GetIntParaFromDB("nm");
            int nTemperatureModel = GlobalDataBase::GetIntParaFromDB("ntmodel");
            int nSpeciesNumber = nLaminar + nChemical - nNSEquation;

            int nTotalVariable = nTemperatureModel + 3;
            CommunicateSpecificArray(iSolver, level, "gradUVWTCellCenterX", nTotalVariable);
            CommunicateSpecificArray(iSolver, level, "gradUVWTCellCenterY", nTotalVariable);
            CommunicateSpecificArray(iSolver, level, "gradUVWTCellCenterZ", nTotalVariable);

            if (nChemical > 0 && nSpeciesNumber > 0)
            {
                CommunicateSpecificArray(iSolver, level, "dcdx", nSpeciesNumber);
                CommunicateSpecificArray(iSolver, level, "dcdy", nSpeciesNumber);
                CommunicateSpecificArray(iSolver, level, "dcdz", nSpeciesNumber);
            }
            if(periodicType == 2)
            {
                for (int iZone = 0; iZone < nZones; ++ iZone)
                {
                    int sendProcessor = PHMPI::GetZoneProcessorID(iZone);
                    ZoneNeighbor *globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
                    int numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();
                    for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
                    {
                        int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);
                        int receiveProcessor = PHMPI::GetZoneProcessorID(neighborZone);
                        if (currentProcessor != sendProcessor && currentProcessor != receiveProcessor)
                        {
                            continue;
                        }
                        if (currentProcessor == sendProcessor)
                        {
                            CFDSolver *solver = GetCFDSolver(iZone, iSolver);
                            Grid *grid = PHSPACE::GetGrid(iZone, level);
                            solver->RotateVectorFromInterface(grid, iZone, nTotalVariable);
                        }
                    }
                }
            }
        }
        else if (iSolver == TURB_SOLVER && iLES == NOLES_SOLVER)
        {
            int nEquation = solver->GetNumberOfEquations();
            CommunicateSpecificArray(iSolver, level, "gradTurbulenceCellCenterX", nEquation);
            CommunicateSpecificArray(iSolver, level, "gradTurbulenceCellCenterY", nEquation);
            CommunicateSpecificArray(iSolver, level, "gradTurbulenceCellCenterZ", nEquation);
            if(periodicType == 2)
            {
                for (int iZone = 0; iZone < nZones; ++ iZone)
                {
                    int sendProcessor = PHMPI::GetZoneProcessorID(iZone);
                    ZoneNeighbor *globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
                    int numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();
                    for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
                    {
                        int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);
                        int receiveProcessor = PHMPI::GetZoneProcessorID(neighborZone);
                        if (currentProcessor != sendProcessor && currentProcessor != receiveProcessor)
                        {
                            continue;
                        }
                        if (currentProcessor == sendProcessor)
                        {
                            CFDSolver *solver = GetCFDSolver(iZone, iSolver);
                            Grid *grid = PHSPACE::GetGrid(iZone, level);
                            solver->RotateVectorFromInterface(grid, iZone, nEquation);
                        }
                    }
                }
            }
        }
        else if(iSolver == TURB_SOLVER && iLES == LES_SOLVER)
        {
            //for LES, Gradient Communication should be done another time after NSSolver ends.
            int nTemperatureModel = GlobalDataBase::GetIntParaFromDB("ntmodel");
            int nTotalVariable = nTemperatureModel + 3;

            CommunicateSpecificArray(0, level, "gradUVWTCellCenterX", nTotalVariable);
            CommunicateSpecificArray(0, level, "gradUVWTCellCenterY", nTotalVariable);
            CommunicateSpecificArray(0, level, "gradUVWTCellCenterZ", nTotalVariable);

            if(periodicType == 2)
            {
                for (int iZone = 0; iZone < nZones; ++ iZone)
                {
                    int sendProcessor = PHMPI::GetZoneProcessorID(iZone);
                    ZoneNeighbor *globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
                    int numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();
                    for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
                    {
                        int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);
                        int receiveProcessor = PHMPI::GetZoneProcessorID(neighborZone);
                        if (currentProcessor != sendProcessor && currentProcessor != receiveProcessor)
                        {
                            continue;
                        }
                        if (currentProcessor == sendProcessor)
                        {
                            CFDSolver *solver = GetCFDSolver(iZone, iSolver);
                            Grid *grid = PHSPACE::GetGrid(iZone, level);
                            solver->RotateVectorFromInterface(grid, iZone, nTotalVariable);
                        }
                    }
                }
            }
        }
    }
    else if (systemGridType == MIXGRID && isFVMOrFDM == FV_METHOD)
    {
        if (iSolver == NS_SOLVER)
        {
            CommunicateGradientandLimitOnlyForMixGrid();
        }
        else if (iSolver == TURB_SOLVER)
        {
        }
    }
}

void Controller::BLUSGS(int iSolver, int level)
{
    //! If the solver is Turb Solver, it will take the LUSGS methods.
    if (iSolver == TURB_SOLVER)
    {
        LUSGS(iSolver, level);
        return;
    }
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

    int nLUSGSSweeps = GlobalDataBase::GetIntParaFromDB("nLUSGSSweeps");
    RDouble LUSGSTolerance = GlobalDataBase::GetDoubleParaFromDB("LUSGSTolerance");

    FieldProxy **dqProxy = new FieldProxy *[nLocalZones];

    //! L*Dq or U*Dq, which is equal to DeltaFlux between beighbors.
    FieldProxy **LUplusDQ= new FieldProxy *[nLocalZones];

    BLUSGSInitialization(iSolver, level, nLUSGSSweeps, dqProxy, LUplusDQ);

    RDouble sweepNormalInit = 0.0, sweepEpsilon = 1.1, globleNormal = 0.0;
    
    for (int iSweep = 0; iSweep < nLUSGSSweeps; ++ iSweep)
    {
        RDouble sweepNormal = 0.0;
        ForwardSweepBLUSGS(iSolver, level, dqProxy, LUplusDQ, sweepNormal,iSweep);

        CommunicationDQ(iSolver, level, dqProxy);

        BackwardSweepBLUSGS(iSolver, level, dqProxy, LUplusDQ, sweepNormal);

        if (nLUSGSSweeps > 1)
        {
            //! Get the global normal.
            PH_AllReduce(&sweepNormal, &globleNormal, 1, PH_SUM);
            sweepNormal = globleNormal;
        }

        //! Check the convergence.
        //! This criterion will lead to at least 2 sweep.
        //! It may be more considered, such as using the maximum normal of dq/dqOld.
        if(iSweep == 0)
        {
            sweepNormalInit = sweepNormal;
        }else
        {
            sweepEpsilon = sqrt(sweepNormal / sweepNormalInit);
        }

        if (sweepEpsilon < LUSGSTolerance)
        {
            //! The equations have been converged.
            break;
        }
    }

    UpdateFlowField(iSolver, level, dqProxy);

    for (int iLocalZone = 0; iLocalZone < nLocalZones; ++ iLocalZone)
    {
        delete dqProxy[iLocalZone];
        delete LUplusDQ[iLocalZone];
    }
    delete []dqProxy;
    delete []LUplusDQ;
}

void Controller::ForwardSweepBLUSGS(int iSolver, int level, FieldProxy **dqProxy, FieldProxy **LUplusDQ, RDouble &sweepNormal , int &iSweep)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();
    RDouble *zoneNormal = new RDouble [nLocalZones]();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Zone *zone = GetZone(zoneID);
        if (!zone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        if (iSweep == 0)
        {
            solver->LHSMatrixs(grid);
        }

        zoneNormal[iZone] = 0.0;
        solver->SolveBLUSGSForward(grid, dqProxy[iZone], LUplusDQ[iZone], zoneNormal[iZone]);//! The Block LU-SGS method.
    }

    //! Sum the zone normal, which is to void OMP data competition.
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        sweepNormal += zoneNormal[iZone];
    }
    delete []zoneNormal;    zoneNormal = NULL;
}

void Controller::BackwardSweepBLUSGS(int iSolver, int level, FieldProxy **dqProxy, FieldProxy **LUplusDQ, RDouble &sweepNormal)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();
    RDouble *zoneNormal = new RDouble [nLocalZones]();
#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Zone *zone = GetZone(zoneID);
        if (!zone)
        {
            continue;
        }
        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        zoneNormal[iZone] = 0.0;
        solver->SolveBLUSGSBackward(grid, dqProxy[iZone], LUplusDQ[iZone], zoneNormal[iZone]);//The Block LU-SGS method
    }

    //! Sum the zone normal, which is to void OMP data competition.
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        sweepNormal += zoneNormal[iZone];
    }
    delete []zoneNormal;    zoneNormal = NULL;
}

void Controller::MLUSGS(int iSolver, int level)
{
    //! If the solver is Turb Solver, it will take the LUSGS methods.
    if (iSolver == TURB_SOLVER)
    {
        LUSGS(iSolver, level);
        return;
    }
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

    int nLUSGSSweeps = GlobalDataBase::GetIntParaFromDB("nLUSGSSweeps");
    RDouble LUSGSTolerance = GlobalDataBase::GetDoubleParaFromDB("LUSGSTolerance");

    FieldProxy **dqProxy = new FieldProxy *[nLocalZones];

    //! L*Dq or U*Dq, which is equal to DeltaFlux between beighbors.
    FieldProxy **LUplusDQ= new FieldProxy *[nLocalZones];

    MLUSGSInitialization(iSolver, level, nLUSGSSweeps, dqProxy, LUplusDQ);

    RDouble sweepNormalInit, sweepEpsilon = 1.1, globleNormal;
    
    for (int iSweep = 0; iSweep < nLUSGSSweeps; ++ iSweep)
    {
        RDouble sweepNormal = 0.0;
        ForwardSweepMLUSGS(iSolver, level, dqProxy, LUplusDQ, sweepNormal,iSweep);

		CommunicationDQ(iSolver, level, dqProxy);

        BackwardSweepMLUSGS(iSolver, level, dqProxy, LUplusDQ, sweepNormal);		

        if (nLUSGSSweeps > 1)
        {
            //! Get the global normal.
            PH_AllReduce(&sweepNormal, &globleNormal, 1, PH_SUM);
            sweepNormal = globleNormal;
			
        }

        //! Check the convergence.
        //! This criterion will lead to at least 2 sweep.
        //! It may be more considered, such as using the maximum normal of dq/dqOld.
        if(iSweep == 0)
        {
            sweepNormalInit = sweepNormal;
        }else
        {
            sweepEpsilon = sqrt(sweepNormal / sweepNormalInit);
        }

        if (sweepEpsilon < LUSGSTolerance)
        {
            //! The equations have been converged.
            break;
        }
    }

	UpdateFlowField(iSolver, level, dqProxy);

    for (int iLocalZone = 0; iLocalZone < nLocalZones; ++ iLocalZone)
    {
        delete dqProxy[iLocalZone];
        delete LUplusDQ[iLocalZone];
    }
    delete []dqProxy;
    delete []LUplusDQ;
}

void Controller::ForwardSweepMLUSGS(int iSolver, int level, FieldProxy **dqProxy, FieldProxy **LUplusDQ, RDouble &sweepNormal , int &iSweep)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();
    RDouble *zoneNormal = new RDouble [nLocalZones];

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int izone = 0; izone < nLocalZones; ++ izone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(izone);

        Zone *zone = GetZone(zoneID);
        if (!zone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        if (iSweep == 0)
        {
            solver->LHSMatrixs(grid);
        }

        zoneNormal[izone] = 0.0;
		solver->SolveMatrixLUSGSForward(grid, dqProxy[izone], LUplusDQ[izone], zoneNormal[izone]); //! The Matrix LU-SGS method.
    }

    //! Sum the zone normal, which is to void OMP data competition.
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        sweepNormal += zoneNormal[iZone];
    }
    delete []zoneNormal;
}

void Controller::BackwardSweepMLUSGS(int iSolver, int level, FieldProxy **dqProxy, FieldProxy **LUplusDQ, RDouble &sweepNormal)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();
    RDouble *zoneNormal = new RDouble [nLocalZones];
#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int izone = 0; izone < nLocalZones; ++ izone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(izone);

        Zone *zone = GetZone(zoneID);
        if (!zone)
        {
            continue;
        }
        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        zoneNormal[izone] = 0.0;
        solver->SolveMatrixLUSGSBackward(grid, dqProxy[izone], LUplusDQ[izone], zoneNormal[izone]);//The Matrix LU-SGS method
    }

    //! Sum the zone normal, which is to void OMP data competition.
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        sweepNormal += zoneNormal[iZone];
    }
    delete []zoneNormal;
}

void Controller::LineLUSGS(int iSolver, int level)
{
    //! If the solver is Turb Solver, it will take the LUSGS methods.
    if (iSolver == TURB_SOLVER)
    {
        LUSGS(iSolver, level);
        return;
    }
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

    int nLUSGSSweeps = GlobalDataBase::GetIntParaFromDB("nLUSGSSweeps");
    RDouble LUSGSTolerance = GlobalDataBase::GetDoubleParaFromDB("LUSGSTolerance");

    FieldProxy** dqProxy = new FieldProxy * [nLocalZones];

    //! L*Dq or U*Dq, which is equal to DeltaFlux between beighbors.
    FieldProxy** LUplusDQ = new FieldProxy * [nLocalZones];

    LineLUSGSInitialization(iSolver, level, nLUSGSSweeps, dqProxy, LUplusDQ);

    RDouble sweepNormalInit = 0.0, sweepEpsilon = 1.1, globleNormal = 0.0;

    for (int iSweep = 0; iSweep < nLUSGSSweeps; ++iSweep)
    {
        RDouble sweepNormal = 0.0;
        ForwardSweepLineLUSGS(iSolver, level, dqProxy, LUplusDQ, sweepNormal, iSweep);

        CommunicationDQ(iSolver, level, dqProxy);

        BackwardSweepLineLUSGS(iSolver, level, dqProxy, LUplusDQ, sweepNormal);

        if (nLUSGSSweeps > 1)
        {
            //! Get the global normal.
            PH_AllReduce(&sweepNormal, &globleNormal, 1, PH_SUM);
            sweepNormal = globleNormal;
        }

        //! Check the convergence.
        //! This criterion will lead to at least 2 sweep.
        //! It may be more considered, such as using the maximum normal of dq/dqOld.
        if (iSweep == 0)
        {
            sweepNormalInit = sweepNormal;
        }
        else
        {
            sweepEpsilon = sweepNormal / sweepNormalInit;
        }

        if (sweepEpsilon < LUSGSTolerance)
        {
            //! The equations have been converged.
            break;
        }
    }

    UpdateFlowField(iSolver, level, dqProxy);

    for (int iLocalZone = 0; iLocalZone < nLocalZones; ++iLocalZone)
    {
        delete dqProxy[iLocalZone];
        delete LUplusDQ[iLocalZone];
    }
    delete [] dqProxy; 
    delete [] LUplusDQ;
}

void Controller::LineLUSGSInitialization(int iSolver, int level, int nSweep, FieldProxy** dqProxy, FieldProxy** LUplusDQ)
{
    int nLocalZones = PHMPI::GetNumberofLocalZones();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++iZone)
    {
        int zoneID = PHMPI::GetLocalZoneIDToGlobalZoneID(iZone);

        Zone* zone = GetZone(zoneID);
        if (!zone)
        {
            continue;
        }

        Grid* grid = GetGrid(zoneID, level);
        CFDSolver* solver = GetCFDSolver(zoneID, iSolver);

        //! Memory allocating.
        dqProxy[iZone] = solver->CreateFieldProxy(grid);
        LUplusDQ[iZone] = solver->CreateFieldProxy(grid);

        //! The delta Flux is initialized by zero, which will be used in the forward sweep.
        solver->FillField(grid, LUplusDQ[iZone], zero);

        //* warning: it is important to research the effect of the initialization in future!
        //! e.g., using rhs to init when nsweep > 1.
        //! Initialization.
        if (nSweep == 1)
        {
            //! Single sweep, init by res.
            FieldProxy* rhsProxy = solver->GetFieldProxy(grid, solver->CastName("res"));
            solver->FillField(grid, dqProxy[iZone], rhsProxy);
            delete rhsProxy;
        }
        else
        {
            //! Multi-sweep, init by zero.
            solver->FillField(grid, dqProxy[iZone], zero);
        }
    }
}

void Controller::ForwardSweepLineLUSGS(int iSolver, int level, FieldProxy** dqProxy, FieldProxy** LUplusDQ, RDouble& sweepNormal, int& iSweep)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();
    RDouble* zoneNormal = new RDouble[nLocalZones]();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Zone* zone = GetZone(zoneID);
        if (!zone)
        {
            continue;
        }

        Grid* grid = GetGrid(zoneID, level);
        CFDSolver* solver = GetCFDSolver(zoneID, iSolver);

        if (iSweep == 0)
        {
            solver->LineLUSGSMatrixs(grid);
        }

        zoneNormal[iZone] = 0.0;
        solver->SolveLineLUSGSForward(grid, dqProxy[iZone], LUplusDQ[iZone], zoneNormal[iZone]);//! The Block LU-SGS method.
    }

    //! Sum the zone normal, which is to void OMP data competition.
    for (int iZone = 0; iZone < nLocalZones; ++iZone)
    {
        sweepNormal += zoneNormal[iZone];
    }
    delete [] zoneNormal;    zoneNormal = NULL;
}

void Controller::BackwardSweepLineLUSGS(int iSolver, int level, FieldProxy** dqProxy, FieldProxy** LUplusDQ, RDouble& sweepNormal)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();
    RDouble* zoneNormal = new RDouble[nLocalZones]();
#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Zone* zone = GetZone(zoneID);
        if (!zone)
        {
            continue;
        }
        Grid* grid = GetGrid(zoneID, level);
        CFDSolver* solver = GetCFDSolver(zoneID, iSolver);

        zoneNormal[iZone] = 0.0;
        solver->SolveLineLUSGSBackward(grid, dqProxy[iZone], LUplusDQ[iZone], zoneNormal[iZone]);//The Block LU-SGS method
    }

    //! Sum the zone normal, which is to void OMP data competition.
    for (int iZone = 0; iZone < nLocalZones; ++iZone)
    {
        sweepNormal += zoneNormal[iZone];
    }
    delete [] zoneNormal;    zoneNormal = NULL;
}

void Controller::RestrictAllQ(int iSolver, int level)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Zone *zone = GetZone(zoneID);
        if (!zone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        Grid *coarseGrid = grid->GetCoarseGrid();
        if (coarseGrid)
        {
            CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

            solver->RestrictAllQ(grid, coarseGrid);
        }
    }
}

void Controller::RestrictDefect(int iSolver, int level)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Zone *zone = GetZone(zoneID);
        if (!zone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID,level);
        Grid * coarseGrid = grid->GetCoarseGrid();
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        solver->RestrictDefect(grid, coarseGrid);
    }
}

void Controller::PutCorrection(int iSolver, int level)
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    FieldProxy **qProxy = GetQProxy(level);
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        Zone *zone = GetZone(iZone);
        if (!zone)
        {
            continue;
        }

        Grid *grid = GetGrid(iZone,level);
        CFDSolver *solver = GetCFDSolver(iZone,iSolver);

        solver->PutCorrection(grid, qProxy[iZone]);
    }
}

void Controller::PutCorrectionBack(int iSolver, int level)
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    FieldProxy **qProxy = GetQProxy(level);
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        Zone *zone = GetZone(iZone);
        if (!zone)
        {
            continue;
        }

        Grid *grid = GetGrid(iZone,level);
        CFDSolver *solver = GetCFDSolver(iZone,iSolver);

        solver->PutCorrectionBack(grid, qProxy[iZone]);
    }
}

void Controller::CorrectFineGrid(int iSolver, int level)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Zone *zone = GetZone(zoneID);
        if (!zone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        solver->CorrectFineGrid(grid, grid->GetCoarseGrid());
    }
}

Zone * Controller::GetZone(int iZone)
{
    return region->GetZone(iZone);
}

CFDSolver * Controller::GetCFDSolver(int iZone, int iSolver)
{
    return static_cast<CFDSolver *>(GlobalSolvers::GetSolver(iZone, iSolver));
}

bool Controller::CVGNorm(int iSolver, int level)
{
    int innstep = 0;
    GlobalDataBase::GetData("innstep", &innstep, PHINT, 1);

    RDouble tol_sub_iter = 0.01;
    GlobalDataBase::GetData("tol_sub_iter", &tol_sub_iter, PHDOUBLE, 1);

    int min_sub_iter = 1;
    GlobalDataBase::GetData("min_sub_iter", &min_sub_iter, PHINT, 1);

    int max_sub_iter = 1;
    GlobalDataBase::GetData("max_sub_iter", &max_sub_iter, PHINT, 1);

    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    RDouble cvgmax = zero, cvg;
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        Zone *zone = GetZone(iZone);
        if (!zone)
        {
            continue;
        }

        Grid *grid = GetGrid(iZone,level);
        CFDSolver *solver = GetCFDSolver(iZone,iSolver);

        cvg = solver->UnsteadyConvergence(grid);
        cvgmax = MAX(cvgmax, cvg);
    }

    //cout << "iner_iter = " << iner_iter << " cvg = " << cvg << "\n";

    GlobalSolve();

    PH_CompareMaxMin(cvgmax, 1);

    if (innstep <  min_sub_iter) return false;
    if (cvgmax  <= tol_sub_iter) return true ;
    if (innstep >= max_sub_iter) return true ;
    return false;
}

RDouble Controller::ComputeSubIterationNorm(int iSolver, int level)
{
    RDouble subIterationNorm    =  zero;
    RDouble maxSubIterationNorm = -LARGE;

    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        Zone *zone = GetZone(iZone);
        if (!zone)
        {
            continue;
        }

        Grid *grid    = GetGrid(iZone, level);
        CFDSolver *solver  = GetCFDSolver(iZone, iSolver);

        subIterationNorm    = solver->UnsteadyConvergence(grid);

        maxSubIterationNorm = MAX(maxSubIterationNorm, subIterationNorm);
    }

    //PH_CompareMaxMin(averageQuotientMAX, 1);

    return maxSubIterationNorm;
}

int Controller::GetNPostSolve(int iSolver, int level)
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        Zone *zone = GetZone(iZone);
        if (!zone)
        {
            continue;
        }

        CFDSolver *solver = GetCFDSolver(iZone,iSolver);

        return solver->GetNPostSolve();
    }

    return 1;
}

void Controller::PostSolve(int level)
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();
    int iSolver = this->GetSolverIndex();
    int nPostSolve = GetNPostSolve(iSolver, level);

    //! Face Communication.
    for (int n = 0; n < nPostSolve; ++ n)
    {
        for (int iZone = 0; iZone < nZones; ++ iZone)
        {
            Zone *zone = GetZone(iZone);
            if (!zone)
            {
                continue;
            }

            Grid *grid = GetGrid(iZone,level);
            CFDSolver *solver = GetCFDSolver(iZone, iSolver);

            solver->PostSolve(grid, n, level);
        }
    }

    SetBoundaryValueOfField(iSolver, level);

    //! Post tasks, such as residual computing.
    int zoneIndex = GetLocalZoneIDToGlobalZoneID(0);
    Grid *grid = GetGrid(zoneIndex, level);
    CFDSolver *solver = static_cast<CFDSolver *>(GlobalSolvers::GetSolver(zoneIndex, iSolver));
    solver->DumpResultFile(grid, level);

     //exam the change of key variables such as the pressure, the translational temperature and the vibrational temperature, etc.
    int isAdaptiveSolver = GlobalDataBase::GetIntParaFromDB("isAdaptiveSolver");
    if (isAdaptiveSolver > 0)
    {
        solver->CheckResult(grid, level);
    }

    //dump the surface heating change.
    int nSurfHeatMonitor = 0;
    if (GlobalDataBase::IsExist("nSurfHeatMonitor", PHINT, 1))
    {
        nSurfHeatMonitor = GlobalDataBase::GetIntParaFromDB("nSurfHeatMonitor");
    }
    if (nSurfHeatMonitor > 0)
    {
        solver->CheckSurfaceHeatingChange(grid, level);
    }

    //dump the CFL number to the file.
    int nDumpCFLNumber = 0;
    if (GlobalDataBase::IsExist("nDumpCFLNumber", PHINT, 1))
    {
        nDumpCFLNumber = GlobalDataBase::GetIntParaFromDB("nDumpCFLNumber");
    }
    if (nDumpCFLNumber > 0)
    {
        solver->DumpCFLNumber(grid, level);
    }

    //! Node communication: part1.
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        Zone *zone = GetZone(iZone);
        if (!zone)
        {
            continue;
        }

        CFDSolver *solver = GetCFDSolver(iZone, iSolver);
        for (int iLevel = 0; iLevel < GetNumberOfMultiGrid(); ++ iLevel)
        {
            Grid *grid = GetGrid(iZone, iLevel);
            solver->ComputeNodeValue(grid);
            solver->ComputeQTurbNodeValue(grid);
        }
    }

    //! Node communication: part2.
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        Zone *zone = GetZone(iZone);
        if (!zone)
        {
            continue;
        }
        Grid *grid = GetGrid(iZone,level);
        CFDSolver *solver = GetCFDSolver(iZone, iSolver);
        InterpointInformation *inpterPointInfor = grid->GetInterpointInfo();
        if (inpterPointInfor)
        {
            solver->CommunicationInterpointData();
        }
    }

    //! Node communication: part3.
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        Zone *zone = GetZone(iZone);
        if (! zone)
        {
            continue;
        }

        CFDSolver *solver = GetCFDSolver(iZone, iSolver);
        for (int iLevel = 0; iLevel < GetNumberOfMultiGrid(); ++ iLevel)
        {
            Grid *grid = GetGrid(iZone, iLevel);
            solver->ModifyNodeValue(grid);
            solver->ModifyQTurbNodeValue(grid);
        }
    }
}

void Controller::FillField(int iSolver, int level, FieldProxy **targetProxy, FieldProxy **sourceProxy)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);
        Zone *zone = GetZone(zoneID);
        if (!zone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID,level);
        CFDSolver *solver = GetCFDSolver(zoneID,iSolver);

        solver->FillField(grid, targetProxy[iZone], sourceProxy[zoneID]);
    }
}

void Controller::RungeKutta(int iSolver, int level)
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    FieldProxy **rightHandSideProxy = GetRHSProxy(level);

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        Zone *zone = GetZone(iZone);
        if (!zone)
        {
            continue;
        }

        Grid *grid = GetGrid(iZone,level);
        CFDSolver *solver = GetCFDSolver(iZone,iSolver);

        solver->RungeKutta(grid, rightHandSideProxy[iZone]);
    }
}

void Controller::InitResidual(int iSolver, int level)
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    FieldProxy **rightHandSideProxy = GetRHSProxy(level);
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        Zone *zone = GetZone(iZone);
        if (!zone)
        {
            continue;
        }

        Grid *grid = GetGrid(iZone,level);
        CFDSolver *solver = GetCFDSolver(iZone,iSolver);

        solver->InitResidual(grid, rightHandSideProxy[iZone]);
    }
}

void Controller::CreateQProxy(int iSolver, int level)
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    FieldProxy **qProxy = GetQProxy(level);
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        Zone *zone = GetZone(iZone);
        if (!zone)
        {
            continue;
        }

        Grid *grid = GetGrid(iZone, level);
        CFDSolver *solver = GetCFDSolver(iZone, iSolver);

        qProxy[iZone] = solver->CreateFieldProxy(grid);
    }
}

void Controller::DestroyQProxy(int iSolver, int level)
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    FieldProxy **qProxy = GetQProxy(level);
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        Zone *zone = GetZone(iZone);
        if (!zone)
        {
            continue;
        }
        delete qProxy[iZone];
        qProxy[iZone] = NULL;
    }
}

void Controller::CreateRHSProxy(int iSolver, int level)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

    FieldProxy **rightHandSideProxy = GetRHSProxy(level);

#ifdef USE_OMP
    #pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Zone *zone = GetZone(zoneID);
        if (!zone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        rightHandSideProxy[zoneID] = solver->CreateFieldProxy(grid);
    }
}

void Controller::DestroyRHSProxy(int iSolver, int level)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

    FieldProxy **rightHandSideProxy = GetRHSProxy(level);
#ifdef USE_OMP
    #pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Zone *zone = GetZone(zoneID);
        if (!zone)
        {
            continue;
        }

        delete rightHandSideProxy[zoneID];
        rightHandSideProxy[zoneID] = NULL;
    }
}

void Controller::StoreRhsByResidual(int iSolver, int level)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

    FieldProxy **rightHandSideProxy = GetRHSProxy(level);

#ifdef USE_OMP
    #pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Zone *zone = GetZone(zoneID);
        if (!zone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID, level);
        CFDSolver *solver = GetCFDSolver(zoneID, iSolver);

        solver->StoreRhsByResidual(grid, rightHandSideProxy[zoneID]);
    }
}

void Controller::RecoverResidual(int iSolver, int level)
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

    FieldProxy **rightHandSideProxy = GetRHSProxy(level);
#ifdef USE_OMP
    #pragma omp parallel for
#endif
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Zone *zone = GetZone(zoneID);
        if (!zone)
        {
            continue;
        }

        Grid *grid = GetGrid(zoneID,level);
        CFDSolver *solver = GetCFDSolver(zoneID,iSolver);

        solver->RecoverResidual(grid,rightHandSideProxy[zoneID]);
    }
}

FieldProxy **Controller::GetQProxy(int level)
{
    return q_proxy->GetFieldProxy(level);
}

FieldProxy **Controller::GetRHSProxy(int level)
{
    return rhs_proxy->GetFieldProxy(level);
}

void Controller::LoadQ(int iSolver, int level)
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    FieldProxy **qProxy = GetQProxy(level);
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        Zone *zone = GetZone(iZone);
        if (!zone)
        {
            continue;
        }

        Grid *grid = GetGrid(iZone, level);
        CFDSolver *solver = GetCFDSolver(iZone, iSolver);

        solver->LoadQ(grid, qProxy[iZone]);
    }
}

bool Controller::ISCoarsestGrid(int level)
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        Zone *zone = GetZone(iZone);
        if (!zone)
        {
            continue;
        }

        Grid *grid = GetGrid(iZone,level);
        if (grid->GetCoarseGrid())
        {
            return false;
        }
        else
        {
            return true;
        }
    }
    return false;
}

LIB_EXPORT void Controller::InterpolatFineGrid(int iSolver, int level)
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        Zone *zone = GetZone(iZone);
        if (!zone)
        {
            continue;
        }

        Grid *grid = GetGrid(iZone, level);
        CFDSolver *solver = GetCFDSolver(iZone, iSolver);

        Grid *coarseGrid = grid->GetCoarseGrid();

        solver->InterpolatFineGrid(grid, coarseGrid);
        solver->CommunicationInterfaceData();
    }
}

LIB_EXPORT void Controller::MultiGridInitFlow(int level)
{
    int iSolver = this->GetSolverIndex();

    ZeroResiduals(iSolver, level);

    //! Create Right Hand Side term.
    CreateRHSProxy(iSolver, level);
    StoreRhsByResidual(iSolver, level);

    //! Relaxation.
    Relaxation(iSolver, level);

    //! Recover the residual.
    DestroyRHSProxy(iSolver, level);
}

//! Incompressible Flow
void Controller::SolveFlowField(int iZone)
{

    //! Discret momentum equation.
    CFDSolver *solver = GetCFDSolver(iZone, this->solverID);
    solver->DiscrtMomentumEquations();
    solver->SolveMomentumEquations(*this->region);

    //! Update the velocity and form the coefficient of pressure equation.
    solver->UpdateVelocity();
    solver->DiscrtPressureCorrectionEquation();
    solver->SolvePressureCorrectionEquation(*this->region);

    //! Update the pressure and velocity.
    solver->UpdatePressureAndVelocity();

    //! Update the properties.
    solver->UpdateProperties();
}

//! Incompressible Flow
void Controller::SolveScalarField(int iZone)
{
    CFDSolver *solver = GetCFDSolver(iZone, this->solverID);

    solver->Discrt();

    solver->SolveScalarEquation(*this->region);

    solver->Update();

    solver->UpdateBoundary();

    solver->UpdateProperties();

    solver->UpdateAfterIterloop();
}

//! Incompressible Flow
void Controller::SolveTurbField(int iZone)
{
    int viscousType = GlobalDataBase::GetIntParaFromDB("viscousType");
    if (viscousType == SPALART_ALLMARAS)
    {
        SolveScalarField(iZone);
    }
    else if (viscousType == KEPSILON)
    {
        CFDSolver *solver = GetCFDSolver(iZone, this->solverID);

        solver->SolveKEquation(*this->region);

        solver->SolveEpsilonEquation(*this->region);

        solver->UpdateAfterIterloop();
    }
}

//! Incompressible Flow
void Controller::SolveEnergyField(int iZone)
{
    SolveScalarField(iZone);
}

bool Controller::IsConverged(int iZone)
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();
    bool isConverged;

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        Zone *zone = GetZone(iZone);
        if (!zone)
        {
            continue;
        }

        Grid *grid = GetGrid(iZone, 0);

        CFDSolver *solver = GetCFDSolver(iZone, this->solverID);

        solver->SolutionIsConverged(iZone, &isConverged);
    }

    return isConverged;
}

MultiGridFieldProxy::MultiGridFieldProxy(int iSolver, int level, int nZone)
{
    this->iSolver = iSolver;
    this->level   = level;
    this->nzone   = nZone;

    field_proxy = new FieldProxy *[nZone];
    for (int iZone = 0; iZone < nZone; ++ iZone)
    {
        field_proxy[iZone] = 0;
    }
}

MultiGridFieldProxy::~MultiGridFieldProxy()
{
    for (int iZone = 0; iZone < nzone; ++ iZone)
    {
        delete field_proxy[iZone];
    }
    delete [] field_proxy;
}

ZoneFieldProxy::ZoneFieldProxy(int iSolver, uint_t nLevel, int nZone)
{
    this->iSolver = iSolver;
    this->nzone   = nZone;
    this->nlevel  = nLevel;

    field_proxy = new MultiGridFieldProxy *[nLevel];

    for (int level = 0; level < nLevel; ++ level)
    {
        field_proxy[level] = new MultiGridFieldProxy(iSolver, level, nZone);
    }
}

ZoneFieldProxy::~ZoneFieldProxy()
{
    for (int level = 0; level < nlevel; ++ level)
    {
        delete field_proxy[level];
    }
    delete []field_proxy;
}

FieldProxy **ZoneFieldProxy::GetFieldProxy(int level)
{
    return field_proxy[level]->GetFieldProxy();
}

}
