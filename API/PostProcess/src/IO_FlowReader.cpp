#include "IO_FlowReader.h"
#include "Glb_Dimension.h"
#include "IO_FileName.h"
#include "GridType.h"
#include "PHIO.h"

using namespace std;
namespace PHSPACE
{

IO_FlowReader::IO_FlowReader(Region *region_in)
{
    this->region = region_in;
    IfReadTurb = false;
}

IO_FlowReader::~IO_FlowReader()
{

}

void IO_FlowReader::Init()
{
    using namespace PHMPI;

    restartNSFile = GlobalDataBase::GetStrParaFromDB("restartNSFile");
    turbfile = GlobalDataBase::GetStrParaFromDB("turbfile");

    if (GetNumberOfProcessor() > 1)
    {
        restartNSFile = AddSymbolToFileName(restartNSFile, "_", 0);
        turbfile = AddSymbolToFileName(turbfile, "_", 0);
    }

    int viscousType = GlobalDataBase::GetIntParaFromDB("viscousType");

    if (viscousType > 1)
    {
        IfReadTurb = true;
    }

    int n_turb = 0;

    if (viscousType == 2)
    {
        n_turb = 0;
    }
    else if (viscousType == 3)
    {
        n_turb = 1;
    }
    else if (viscousType == 4)
    {
        n_turb = 2;
    }

    GlobalDataBase::UpdateData("n_turb", &n_turb, PHINT, 1);


    int nLocalZones = GetNumberofLocalZones();
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Grid *grid = GetGrid(zoneID);
        int gridType = grid->Type();

        if (gridType == STRUCTGRID)
        {
            InitMemoryStr(grid);
        }
        else
        {
            InitMemoryUnStr(grid);
        }
    }
}

void IO_FlowReader::InitMemoryStr(Grid *grid_in)
{
    StructGrid *grid = StructGridCast(grid_in);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();
    Range I(-1, ni+1);
    Range J(-1, nj+1);
    Range K(-1, nk+1);
    if (nk == 1)
    {
        K.setRange(1, 1);
    }

    int nm    = GlobalDataBase::GetIntParaFromDB("nm");
    int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
    int neqn  = nm + nchem;

    int nTemperatureModel = GlobalDataBase::GetIntParaFromDB("ntmodel");
    int nDim = GetDim();

    Range M(0, neqn - 1);

    RFloat4D *q = new RFloat4D(I, J, K, M, fortranArray);
    grid->UpdateDataPtr("q", q);

    RFloat4D *q_turb = new RFloat4D(I, J, K, M, fortranArray);
    grid->UpdateDataPtr("q_turb", q_turb);
}

void IO_FlowReader::InitMemoryUnStr(Grid *grid_in)
{
    UnstructGrid *grid = UnstructGridCast(grid_in);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    int nm = 5;
    GlobalDataBase::GetData("nm", &nm, PHINT, 1);

    int nchem = 0;
    GlobalDataBase::GetData("nchem", &nchem, PHINT, 1);

    int neqn = nm + nchem;

    RFloat **q = NewPointer2<RFloat>(neqn, nTotal);
    RFloat **q_turb = NewPointer2<RFloat>(neqn, nTotal);

    grid->UpdateDataPtr("q", q);
    grid->UpdateDataPtr("q_turb", q_turb);
}

void IO_FlowReader::Run()
{
    Init();
    ReadFlowFile();
    ReadTurbFile();
}

void IO_FlowReader::ReadFlowFile()
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    fstream file;
    ParallelOpenFile(file, restartNSFile, ios_base::in|ios_base::binary);

    for (int izone = 0; izone < nZones; ++ izone)
    {
        int send_proc = server;
        int recv_proc = GetZoneProcessorID(izone);

        int myid = GetCurrentProcessorID();
        int tag = izone;

        DataContainer *flowData = new DataContainer();
        
        if (myid == send_proc)
        {
            ReadFile(file, flowData);
        }

        PH_Trade(flowData, send_proc, recv_proc);

        if (myid == recv_proc)
        {
            Grid *grid = GetGrid(izone);
            int gridType = grid->Type();

            if (gridType == STRUCTGRID)
            {
                ReadFlowFileStr(grid, flowData);
            }
            else
            {
                ReadFlowFileUnStr(grid, flowData);
            }
        }

        delete flowData;
    }

    ParallelCloseFile(file);
}

void IO_FlowReader::ReadFlowFileUnStr(Grid *grid_in, DataContainer *flowData)
{
    UnstructGrid *grid = UnstructGridCast(grid_in);

    int nm    = GlobalDataBase::GetIntParaFromDB("nm");
    int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
    int neqn  = nm + nchem;

    flowData->moveto_begin();
    int outnstep = 0;
    flowData->read(&outnstep, sizeof(int));
    GlobalDataBase::UpdateData("outnstep", &outnstep, PHINT, 1);

    int nBoundFace = grid->GetNBoundFace();
    int nTotalCell = grid->GetNTotalCell();
    int nTotal = nTotalCell + nBoundFace;

    int nIFace = grid->GetNIFace();

    RFloat **q = reinterpret_cast<RFloat **>(grid->GetDataPtr("q"));

    for (int m = 0; m < neqn; ++ m)
    {
        flowData->read(q[m], sizeof(RFloat) * nTotal);
    }
}

void IO_FlowReader::ReadFlowFileStr(Grid *grid_in, DataContainer *flowData)
{
    StructGrid *grid = StructGridCast(grid_in);

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 2);

    int nm    = GlobalDataBase::GetIntParaFromDB("nm");
    int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
    int neqn  = nm + nchem;

    Range M(0, neqn - 1);
    int mst = M.first();
    int med = M.last();

    flowData->moveto_begin();

    int outnstep = 0;
    flowData->read(&outnstep, sizeof(int));
    GlobalDataBase::UpdateData("outnstep", &outnstep, PHINT, 1);

    RFloat4D &q = *reinterpret_cast<RFloat4D *>(grid->GetDataPtr("q"));

    for (int m = mst; m <= med; ++ m)
    {
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    flowData->read(&q(i, j, k, m), sizeof(RFloat));
                }
            }
        }
    }
}

void IO_FlowReader::ReadTurbFile()
{
    if (IfReadTurb == false) return;

    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    fstream file;

    ParallelOpenFile(file, turbfile, ios_base::in|ios_base::binary);

    for (int izone = 0; izone < nZones; ++ izone)
    {
        int send_proc = server;
        int recv_proc = GetZoneProcessorID(izone);

        int myid = GetCurrentProcessorID();
        int tag = izone;

        DataContainer *flowData = new DataContainer();
        if (myid == send_proc)
        {
            ReadFile(file, flowData);
        }

        PH_Trade(flowData, send_proc, recv_proc);

        if (myid == recv_proc)
        {
            Grid *grid = GetGrid(izone, 0);
            int gridType = grid->Type();

            if (gridType == STRUCTGRID)
            {
                ReadTurbFileStr(grid, flowData);
            }
            else
            {
                ReadTurbFileUnStr(grid, flowData);
            }
        }

        delete flowData;
    }

    ParallelCloseFile(file);
}

void IO_FlowReader::ReadTurbFileStr(Grid *grid_in, DataContainer *flowData)
{
    StructGrid *grid = StructGridCast(grid_in);

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 2);

    int n_turb = GlobalDataBase::GetIntParaFromDB("n_turb");

    RFloat4D &q = *reinterpret_cast<RFloat4D *>(grid->GetDataPtr("q_turb"));

    Range M(0, n_turb - 1);
    int mst = M.first();
    int med = M.last();

    flowData->moveto_begin();
    int outnstepofTurb = 0;
    flowData->read(&outnstepofTurb, sizeof(int));

    int outnstepofNS = GlobalDataBase::GetIntParaFromDB("outnstep");

    for (int m = mst; m <= med; ++ m)
    {
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    flowData->read(&q(i, j, k, m), sizeof(RFloat));
                }
            }
        }
    }
}

void IO_FlowReader::ReadTurbFileUnStr(Grid *grid_in, DataContainer *flowData)
{
    UnstructGrid *grid = UnstructGridCast(grid_in);

    int nBoundFace = grid->GetNBoundFace();
    int nTotalCell = grid->GetNTotalCell();
    int nTotal = nTotalCell + nBoundFace;

    int n_turb = GlobalDataBase::GetIntParaFromDB("n_turb");

    RFloat **q = reinterpret_cast<RFloat **>(grid->GetDataPtr("q_turb"));

    flowData->moveto_begin();
    int outnstepofTurb = 0;
    flowData->read(&outnstepofTurb, sizeof(int));

    int outnstepofNS = GlobalDataBase::GetIntParaFromDB("outnstep");

    for (int m = 0; m < n_turb; ++ m)
    {
        for (int icell = 0; icell < nTotal; ++ icell)
        {
            flowData->read(&q[m][icell], sizeof(RFloat));
        }
    }
}

}