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

    this->nEquation           = 0;
    this->nTemperatureModel   = 0;
    this->nTurbulenceEquation = 0;
}

IO_FlowReader::~IO_FlowReader()
{

}

void IO_FlowReader::Run()
{
    Init();
    ReadFlowFile();
    ReadTurbFile();
}

void IO_FlowReader::Init()
{
    using namespace PHMPI;

    restartNSFile = GlobalDataBase::GetStrParaFromDB("restartNSFile");
    turbfile = GlobalDataBase::GetStrParaFromDB("turbfile");

    int nChemical           = GlobalDataBase::GetIntParaFromDB("nchem");
    int nNSEquation         = GlobalDataBase::GetIntParaFromDB("nm");
    this->nTemperatureModel = GlobalDataBase::GetIntParaFromDB("ntmodel");
    this->nEquation         = nNSEquation + nChemical + nTemperatureModel - 1;

    int viscousType = GlobalDataBase::GetIntParaFromDB("viscousType");
    if (viscousType == 3)
    {
        this->nTurbulenceEquation = 1;
    }
    else if (viscousType == 4)
    {
        this->nTurbulenceEquation = 2;
        int transitionType = GlobalDataBase::GetIntParaFromDB("transitionType");
        this->nTurbulenceEquation += transitionType;
    }

    int nLocalZones = GetNumberofLocalZones();
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Grid *grid = GetGrid(zoneID);
        int gridType = grid->Type();
        if (gridType == STRUCTGRID)
        {
            InitMemoryStruct(grid);
        }
        else
        {
            InitMemoryUnstruct(grid);
        }
    }
}

void IO_FlowReader::InitMemoryStruct(Grid *grid_in)
{
    StructGrid *grid = StructGridCast(grid_in);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    Range M(0, nEquation - 1);
    Range T(0, nTemperatureModel - 1);

    RDouble4D *q = new RDouble4D(I, J, K, M, fortranArray);
    grid->UpdateDataPtr("q", q);

    RDouble4D *t = new RDouble4D(I, J, K, T, fortranArray);
    grid->UpdateDataPtr("t", t);

    if (nTurbulenceEquation > 0)
    {
        Range Mturb(0, nTurbulenceEquation - 1);
        RDouble4D *q_turb = new RDouble4D(I, J, K, Mturb, fortranArray);
        grid->UpdateDataPtr("q_turb", q_turb);

        RDouble3D *viscousTurbulent = new RDouble3D(I, J, K, fortranArray);
        grid->UpdateDataPtr("vist", viscousTurbulent);
    }
}

void IO_FlowReader::InitMemoryUnstruct(Grid *grid_in)
{

}

void IO_FlowReader::ReadFlowFile()
{
    hid_t file;
    file = OpenHDF5File(restartNSFile);

    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Grid *grid = GetGrid(zoneID);
        int gridType = grid->Type();
        if (gridType == STRUCTGRID)
        {
            ReadFlowFileStr(file, grid);
        }
        else
        {

        }
    }

    H5Fclose(file);
}

void IO_FlowReader::ReadFlowFileUnStr(Grid *grid_in, DataContainer *flowData)
{

}

void IO_FlowReader::ReadFlowFileStr(hid_t loc_id, Grid *grid_in)
{
    RDouble4D &temperatures       = *reinterpret_cast <RDouble4D *> (grid_in->GetDataPtr("t"));
    RDouble4D &primitiveVariables = *reinterpret_cast <RDouble4D *> (grid_in->GetDataPtr("q"));

    StructGrid *strGrid = StructGridCast(grid_in);
    int iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd;
    strGrid->GetCellIterationIndex(iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd, TWO_GHOST_LAYER);

    hid_t grploc;
    string grpName;

    int outIterStep = 0;
    ReadData(loc_id, &outIterStep, "outnstep");
    GlobalDataBase::UpdateData("outnstep", &outIterStep, PHINT, 1);

    ostringstream oss;
    oss << "Group" << grid_in->GetZoneID();
    grpName = oss.str();

    grploc = OpenGroup(loc_id, grpName);
    ReadData(grploc, &primitiveVariables(iCellStart, jCellStart, kCellStart, 0), "q");
    ReadData(grploc, &temperatures(iCellStart, jCellStart, kCellStart, 0), "t");

    H5Gclose(grploc);
}

void IO_FlowReader::ReadTurbFile()
{
    if (!nTurbulenceEquation)
    {
        return;
    }

    hid_t file;
    file = OpenHDF5File(turbfile);

    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);

        Grid *grid = GetGrid(zoneID);
        int gridType = grid->Type();
        if (gridType == STRUCTGRID)
        {
            ReadTurbFileStr(file, grid);
        }
        else
        {

        }
    }

    H5Fclose(file);
}

void IO_FlowReader::ReadTurbFileStr(hid_t loc_id, Grid *grid_in)
{
    StructGrid *strGrid = StructGridCast(grid_in);
    int iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd;
    strGrid->GetCellIterationIndex(iCellStart, iCellEnd, jCellStart, jCellEnd, kCellStart, kCellEnd, TWO_GHOST_LAYER);

    RDouble4D &qTurbulent = *reinterpret_cast<RDouble4D *> (strGrid->GetDataPtr("q_turb"));
    RDouble3D &viscousTurbulent = *reinterpret_cast<RDouble3D *> (strGrid->GetDataPtr("vist"));

    hid_t grploc;
    string grpName;

    ostringstream oss;
    oss << "Group" << grid_in->GetZoneID();
    grpName = oss.str();

    grploc = OpenGroup(loc_id, grpName);
    ReadData(grploc, &qTurbulent(iCellStart, jCellStart, kCellStart, 0), "q_turb");
    ReadData(grploc, &viscousTurbulent(iCellStart, jCellStart, kCellStart), "visturb");

    H5Gclose(grploc);
}

void IO_FlowReader::ReadTurbFileUnStr(Grid *grid_in, DataContainer *flowData)
{
    
}

}