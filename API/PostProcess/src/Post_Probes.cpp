#include "Post_Probes.h"
#include "TK_Exit.h"
#include "Pointer.h"
#include "IO_FileName.h"
#include "PHIO.h"
#include "Constants.h"
#include "Glb_Dimension.h"
#include "TK_Log.h"
#include "TK_Time.h"
#include "Gas.h"

using namespace std;
namespace PHSPACE
{
LIB_EXPORT Post_Probes::Post_Probes(int nProbeVariables_in, int *probeVariables_in, Grid *grid_in)
{
    probesVariablesMap.insert(pair<int, string>(PROBE_DENSITY    , "density"));
    probesVariablesMap.insert(pair<int, string>(PROBE_U          , "u"));
    probesVariablesMap.insert(pair<int, string>(PROBE_V          , "v"));
    probesVariablesMap.insert(pair<int, string>(PROBE_W          , "w"));
    probesVariablesMap.insert(pair<int, string>(PROBE_PRESSURE   , "pressure"));
    probesVariablesMap.insert(pair<int, string>(PROBE_TEMPERATURE, "temperature"));

    nProbeVariables = nProbeVariables_in;
    for (int iVar = 0; iVar < nProbeVariables; ++ iVar)
    {
        int varTemp = probeVariables_in[iVar];
        probeVariables.insert(varTemp);
    }

    probeVariablesPtr = new Data_Field();

    this->grid = grid_in;

    if (grid_in)
    {
        gridType = grid_in->Type();
    }
    else
    {
        gridType = -1;
    }
}

LIB_EXPORT Post_Probes::~Post_Probes()
{
    FreeAllProbeVarPtr();
    delete probeVariablesPtr;
}

LIB_EXPORT void Post_Probes::GetASCIIFileHeader(vector<string> &title_tecplot)
{
    if (nProbeVariables <= 0) return;
    int dataMonitorType = GlobalDataBase::GetIntParaFromDB("dataMonitorType");

    if (dataMonitorType == PROBESMONITOR)
    {
        GetProbesASCIIFileHeader(title_tecplot);
    }
    else if (dataMonitorType == LINESMONITOR || dataMonitorType == SURFACESMONITOR)
    {
        GetLinesOrSurfacesASCIIFileHeader(title_tecplot);
    }
    else
    {
        TK_Exit::UnexpectedVarValue("dataMonitorType", dataMonitorType);
    }
}

LIB_EXPORT void Post_Probes::GetProbesASCIIFileHeader(vector<string> &title_tecplot)
{
    ostringstream oss;
    oss << "title = \"Probes Flow Fields Data of PHengLEI\"" << endl;
    oss << "Variables =" << endl;
    oss << "\"iter\"" << endl;
    set<int>::iterator varIter;
    for (varIter = probeVariables.begin(); varIter != probeVariables.end(); ++ varIter)
    {
        int variableType = *varIter;
        string varName = probesVariablesMap[variableType];
        string varNameTemp = "\"" + varName + "\"";
        oss << varNameTemp << endl;
    }
    oss << "\"WallTime\"" << endl;
    title_tecplot.push_back(oss.str());
    return;
}

LIB_EXPORT void Post_Probes::GetLinesOrSurfacesASCIIFileHeader(vector<string> &title_tecplot)
{
    ostringstream oss;
    oss << "title = \"Probes Flow Fields Data of PHengLEI\"" << endl;
    oss << "Variables =" << endl;
    oss << "\"iter\"" << endl;
    oss << "\"x\"" << endl;
    oss << "\"y\"" << endl;
    oss << "\"z\"" << endl;
    set<int>::iterator varIter;
    for (varIter = probeVariables.begin(); varIter != probeVariables.end(); ++ varIter)
    {
        int variableType = *varIter;
        string varName = probesVariablesMap[variableType];
        string varNameTemp = "\"" + varName + "\"";
        oss << varNameTemp << endl;
    }
    oss << "\"WallTime\"" << endl;
    title_tecplot.push_back(oss.str());
    return;
}

//LIB_EXPORT void Post_Probes::UpdateProbesVarPtr(const string &name, RDouble *CellCenterData)
//{
//    if (!grid)
//    {
//        TK_Exit::ExceptionExit("Post_Probe::UpdateProbesVarPtr needs Geometry infomation.");
//    }
//    UnstructGrid * unstructGrid = UnstructGridCast(grid);
//    int ZoneProbesNumber = unstructGrid->GetZoneProbesNumber();
//    if (ZoneProbesNumber == 0) return;
//
//    vector<vector<RDouble> > ProbesCoordinates = unstructGrid->GetProbesCoordinates();
//    vector<int> ProbesCellID = unstructGrid->GetProbesCellID();
//
//    RDouble *q_probes = new RDouble [ZoneProbesNumber];
//    SetField(q_probes, 0.0, ZoneProbesNumber);
//
//    for (int izoneProbe = 0; izoneProbe < ZoneProbesNumber; ++ izoneProbe)
//    {
//        int probeCell = ProbesCellID[izoneProbe];
//        q_probes[izoneProbe] = CellCenterData[probeCell];
//    }
//    probeVariablesPtr->UpdateDataPtr(name, q_probes);
//
//}

//LIB_EXPORT void Post_Probes::UpdateProbesVarPtr(const string &name, RDouble4D &CellCenterData, int index)
//{
//    if (!grid)
//    {
//        TK_Exit::ExceptionExit("Post_Probe::UpdateProbesVarPtr needs Geometry infomation.");
//    }
//
//    StructGrid * structGrid = StructGridCast(grid);
//
//    int ZoneProbesNumber = structGrid->GetZoneProbesNumber();
//    if (ZoneProbesNumber == 0) return;
//
//    vector<vector<RDouble> > ProbesCoordinates = structGrid->GetProbesCoordinates();
//    vector<int> ProbesCellNI = structGrid->GetProbesCellNI();
//    vector<int> ProbesCellNJ = structGrid->GetProbesCellNJ();
//    vector<int> ProbesCellNK = structGrid->GetProbesCellNK();
//
//    RDouble *q_probes = new RDouble [ZoneProbesNumber];
//    SetField(q_probes, 0.0, ZoneProbesNumber);
//
//    for (int izoneProbe = 0; izoneProbe < ZoneProbesNumber; ++ izoneProbe)
//    {
//        int probeNI = ProbesCellNI[izoneProbe];
//        int probeNJ = ProbesCellNJ[izoneProbe];
//        int probeNK = ProbesCellNK[izoneProbe];
//        q_probes[izoneProbe] = CellCenterData(probeNI, probeNJ, probeNK, index);
//    }
//
//    probeVariablesPtr->UpdateDataPtr(name, q_probes);
//
//}

LIB_EXPORT void Post_Probes::UpdateProbesVarPtr(const string &name, RDouble *cellCenterData)
{
    if (!grid)
    {
        TK_Exit::ExceptionExit("Post_Probe::UpdateProbesVarPtr needs Geometry infomation.");
    }

    int probeVariablesInterpolationMethod = GlobalDataBase::GetIntParaFromDB("probeVariablesInterpolationMethod");

    if (probeVariablesInterpolationMethod == CELLVALUE)
    {
        UpdateProbesVarPtrM1(name, cellCenterData);
    }
    else if (probeVariablesInterpolationMethod == CELLSTOPROBE)
    {
        UpdateProbesVarPtrM2(name, cellCenterData);
    }
    else if (probeVariablesInterpolationMethod == NODESTOPROBE)
    {
        UpdateProbesVarPtrM3(name, cellCenterData);
    }
    else
    {
        TK_Exit::UnexpectedVarValue("probeVariablesInterpolationMethod", probeVariablesInterpolationMethod);
    }
}

LIB_EXPORT void Post_Probes::UpdateProbesVarPtrM1(const string &name, RDouble *cellCenterData)
{
    UnstructGrid *unstructGrid = UnstructGridCast(grid);
    int ZoneProbesNumber = unstructGrid->GetZoneProbesNumber();
    if (ZoneProbesNumber == 0) return;

    vector<vector<RDouble> > ProbesCoordinates = unstructGrid->GetProbesCoordinates();
    vector<int> ProbesCellID = unstructGrid->GetProbesCellID();

    RDouble *q_probes = new RDouble [ZoneProbesNumber];
    SetField(q_probes, 0.0, ZoneProbesNumber);

    for (int izoneProbe = 0; izoneProbe < ZoneProbesNumber; ++ izoneProbe)
    {
        int probeCell = ProbesCellID[izoneProbe];
        q_probes[izoneProbe] = cellCenterData[probeCell];
    }
    probeVariablesPtr->UpdateDataPtr(name, q_probes);
}

LIB_EXPORT void Post_Probes::UpdateProbesVarPtrM2(const string &name, RDouble *cellCenterData)
{
    UnstructGrid *unstructGrid = UnstructGridCast(grid);
    int ZoneProbesNumber = unstructGrid->GetZoneProbesNumber();
    if (ZoneProbesNumber == 0) return;

    vector<int> *cell2Cell = unstructGrid->GetCell2Cell();

    vector<int> ProbesCellID = unstructGrid->GetProbesCellID();

    vector<int> interpolationCellsID;
    vector<RDouble> numerator;

    RDouble *q_probes = new RDouble [ZoneProbesNumber];
    SetField(q_probes, 0.0, ZoneProbesNumber);

    for (int izoneProbe = 0; izoneProbe < ZoneProbesNumber; ++ izoneProbe)
    {
        int probeCell = ProbesCellID[izoneProbe];
        NecessaryVarForInterpolation(unstructGrid, izoneProbe, probeCell, numerator, interpolationCellsID, true);

        int nNeighbor = cell2Cell[probeCell].size();
        for (int iNeighbor = 0; iNeighbor < nNeighbor; ++ iNeighbor)
        {
            int neighborCell = cell2Cell[probeCell][iNeighbor];
            NecessaryVarForInterpolation(unstructGrid, izoneProbe, neighborCell, numerator, interpolationCellsID, true);
        }

        int nCells = interpolationCellsID.size();
        RDouble denominator = SUMMARY(numerator, nCells);

        for (int icell = 0; icell < nCells; ++ icell)
        {
            RDouble weightCoef = numerator[icell] / denominator;
            int cellID = interpolationCellsID[icell];
            q_probes[izoneProbe] += cellCenterData[cellID] * weightCoef;
        }
    }
    probeVariablesPtr->UpdateDataPtr(name, q_probes);
}

LIB_EXPORT void Post_Probes::UpdateProbesVarPtrM3(const string &name, RDouble *cellCenterData)
{
    UnstructGrid *unstructGrid = UnstructGridCast(grid);
    int ZoneProbesNumber = unstructGrid->GetZoneProbesNumber();
    if (ZoneProbesNumber == 0) return;

    int nTotalNode = unstructGrid->GetNTotalNode();
    int **cell2Node = unstructGrid->GetCell2NodeArray();
    int *nodeNumberOfEachCell = unstructGrid->GetNodeNumberOfEachCell();

    RDouble *qn = new RDouble [nTotalNode];
    if (name == "u" || name == "v" || name == "w")
    {
        CompNodeVarForVelocity(unstructGrid, qn, cellCenterData);
    }
    else
    {
        CompNodeVar(unstructGrid, qn, cellCenterData);
    }

    vector<int> ProbesCellID = unstructGrid->GetProbesCellID();

    vector<int> interpolationNodesID;
    vector<RDouble> numerator;

    RDouble *q_probes = new RDouble [ZoneProbesNumber];
    SetField(q_probes, 0.0, ZoneProbesNumber);

    int nodePosition = 0;
    for (int izoneProbe = 0; izoneProbe < ZoneProbesNumber; ++ izoneProbe)
    {
        int probeCell = ProbesCellID[izoneProbe];

        for (int iNode = 0; iNode < nodeNumberOfEachCell[probeCell]; ++ iNode)
        {
            int node = cell2Node[probeCell][iNode];
            NecessaryVarForInterpolation(unstructGrid, izoneProbe, node, numerator, interpolationNodesID, false);
        }

        int nNodes = interpolationNodesID.size();
        RDouble denominator = SUMMARY(numerator, nNodes);

        for (int inode = 0; inode < nNodes; ++ inode)
        {
            RDouble weightCoef = numerator[inode] / denominator;
            int nodeID = interpolationNodesID[inode];
            q_probes[izoneProbe] += qn[nodeID] * weightCoef;
        }
    }

    probeVariablesPtr->UpdateDataPtr(name, q_probes);

    delete [] qn;
}

LIB_EXPORT void Post_Probes::UpdateProbesVarPtr(const string &name, RDouble4D &cellCenterData, int index)
{
    if (!grid)
    {
        TK_Exit::ExceptionExit("Post_Probe::UpdateProbesVarPtr needs Geometry infomation.");
    }

    int probeVariablesInterpolationMethod = GlobalDataBase::GetIntParaFromDB("probeVariablesInterpolationMethod");

    if (probeVariablesInterpolationMethod == CELLVALUE)
    {
        UpdateProbesVarPtrM1(name, cellCenterData, index);
    }
    else if (probeVariablesInterpolationMethod == CELLSTOPROBE)
    {
        UpdateProbesVarPtrM2(name, cellCenterData, index);
    }
    else if (probeVariablesInterpolationMethod == NODESTOPROBE)
    {
        UpdateProbesVarPtrM3(name, cellCenterData, index);
    }
    else
    {
        TK_Exit::UnexpectedVarValue("probeVariablesInterpolationMethod", probeVariablesInterpolationMethod);
    }
}

LIB_EXPORT void Post_Probes::UpdateProbesVarPtrM1(const string &name, RDouble4D &cellCenterData, int index)
{
    StructGrid *structGrid = StructGridCast(grid);

    int ZoneProbesNumber = structGrid->GetZoneProbesNumber();
    if (ZoneProbesNumber == 0) return;

    vector<vector<RDouble> > ProbesCoordinates = structGrid->GetProbesCoordinates();
    vector<int> ProbesCellNI = structGrid->GetProbesCellNI();
    vector<int> ProbesCellNJ = structGrid->GetProbesCellNJ();
    vector<int> ProbesCellNK = structGrid->GetProbesCellNK();

    RDouble *q_probes = new RDouble [ZoneProbesNumber];
    SetField(q_probes, 0.0, ZoneProbesNumber);

    for (int izoneProbe = 0; izoneProbe < ZoneProbesNumber; ++ izoneProbe)
    {
        int probeNI = ProbesCellNI[izoneProbe];
        int probeNJ = ProbesCellNJ[izoneProbe];
        int probeNK = ProbesCellNK[izoneProbe];
        q_probes[izoneProbe] = cellCenterData(probeNI, probeNJ, probeNK, index);
    }

    probeVariablesPtr->UpdateDataPtr(name, q_probes);
}

LIB_EXPORT void Post_Probes::UpdateProbesVarPtrM2(const string &name, RDouble4D &cellCenterData, int index)
{
    StructGrid *structGrid = StructGridCast(grid);

    int ZoneProbesNumber = structGrid->GetZoneProbesNumber();
    if (ZoneProbesNumber == 0) return;

    vector<int> ProbesCellNI = structGrid->GetProbesCellNI();
    vector<int> ProbesCellNJ = structGrid->GetProbesCellNJ();
    vector<int> ProbesCellNK = structGrid->GetProbesCellNK();

    vector<vector<int> > interpolationCellsID;
    vector<RDouble> numerator;

    RDouble *q_probes = new RDouble [ZoneProbesNumber];
    SetField(q_probes, 0.0, ZoneProbesNumber);

    int nDim = GetDim();

    for (int izoneProbe = 0; izoneProbe < ZoneProbesNumber; ++ izoneProbe)
    {
        int probeNI = ProbesCellNI[izoneProbe];
        int probeNJ = ProbesCellNJ[izoneProbe];
        int probeNK = ProbesCellNK[izoneProbe];

        NecessaryVarForInterpolation(structGrid, izoneProbe, probeNI, probeNJ, probeNK, numerator, interpolationCellsID, true);

        for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
        {
            int i1, j1, k1;
            GetNsurfIndex(iSurface, i1, j1, k1);

            int il = probeNI - i1;
            int jl = probeNJ - j1;
            int kl = probeNK - k1;

            int ir = probeNI + i1;
            int jr = probeNJ + j1;
            int kr = probeNK + k1;

            NecessaryVarForInterpolation(structGrid, izoneProbe, il, jl, kl, numerator, interpolationCellsID, true);
            NecessaryVarForInterpolation(structGrid, izoneProbe, ir, jr, kr, numerator, interpolationCellsID, true);
        }

        int nCells = interpolationCellsID.size();
        RDouble denominator = SUMMARY(numerator, nCells);

        for (int icell = 0; icell < nCells; ++ icell)
        {
            RDouble weightCoef = numerator[icell] / denominator;
            int cellID_i = interpolationCellsID[icell][0];
            int cellID_j = interpolationCellsID[icell][1];
            int cellID_k = interpolationCellsID[icell][2];
            q_probes[izoneProbe] += cellCenterData(cellID_i, cellID_j, cellID_k, index) * weightCoef;
        }
    }

    probeVariablesPtr->UpdateDataPtr(name, q_probes);
}

LIB_EXPORT void Post_Probes::UpdateProbesVarPtrM3(const string &name, RDouble4D &cellCenterData, int index)
{
    StructGrid *structGrid = StructGridCast(grid);

    int ZoneProbesNumber = structGrid->GetZoneProbesNumber();
    if (ZoneProbesNumber == 0) return;

    int ni = structGrid->GetNI();
    int nj = structGrid->GetNJ();
    int nk = structGrid->GetNK();

    Range I(-1, ni+1);
    Range J(-1, nj+1);
    Range K(-1, nk+1);
    if (nk == 1) K.setRange(1, 1);
    Range M(0, 0);

    RDouble4D *qn = new RDouble4D(I, J, K, M, fortranArray);
    CompNodeVar(structGrid, *qn, 0, cellCenterData, index);

    vector<int> ProbesCellNI = structGrid->GetProbesCellNI();
    vector<int> ProbesCellNJ = structGrid->GetProbesCellNJ();
    vector<int> ProbesCellNK = structGrid->GetProbesCellNK();

    vector<vector<int> > interpolationNodesID;
    vector<RDouble> numerator;

    RDouble *q_probes = new RDouble [ZoneProbesNumber];
    SetField(q_probes, 0.0, ZoneProbesNumber);

    int nDim = GetDim();
    int i1 = 1;
    int j1 = 1;
    int k1 = 1;
    if (nk == 1) k1 = 0;

    for (int izoneProbe = 0; izoneProbe < ZoneProbesNumber; ++ izoneProbe)
    {
        int probeNI = ProbesCellNI[izoneProbe];
        int probeNJ = ProbesCellNJ[izoneProbe];
        int probeNK = ProbesCellNK[izoneProbe];

        for (int k = probeNK; k <= probeNK + k1; ++ k)
        {
            for (int j = probeNJ; j <= probeNJ + j1; ++ j)
            {
                for (int i = probeNI; i <= probeNI + i1; ++ i)
                {
                    NecessaryVarForInterpolation(structGrid, izoneProbe, i, j, k, numerator, interpolationNodesID, false);
                }
            }
        }

        int nNodes = interpolationNodesID.size();
        RDouble denominator = SUMMARY(numerator, nNodes);

        for (int inode = 0; inode < nNodes; ++ inode)
        {
            RDouble weightCoef = numerator[inode] / denominator;
            int nodeID_i = interpolationNodesID[inode][0];
            int nodeID_j = interpolationNodesID[inode][1];
            int nodeID_k = interpolationNodesID[inode][2];
            q_probes[izoneProbe] += (*qn)(nodeID_i, nodeID_j, nodeID_k, 0) * weightCoef;
        }
    }

    probeVariablesPtr->UpdateDataPtr(name, q_probes);
    delete qn;
}

LIB_EXPORT void Post_Probes::NecessaryVarForInterpolation(UnstructGrid *grid, int iprobe, int index, vector<RDouble> &numerator, vector<int> &interpolationIndexs, bool flag)
{
    vector<vector<RDouble> > ProbesCoordinates = grid->GetProbesCoordinates();
    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    RDouble *xn = grid->GetX();
    RDouble *yn = grid->GetY();
    RDouble *zn = grid->GetZ();

    RDouble distance = 0;
    RDouble odistance = 0;

    //! The range of r is from 0.5 to 3.0.
    const int r = 1;

    if (flag)
    {
        distance = DISTANCE(xcc[index] - ProbesCoordinates[iprobe][0], ycc[index] - ProbesCoordinates[iprobe][1], zcc[index] - ProbesCoordinates[iprobe][2]) + TINY;
    }
    else
    {
        distance = DISTANCE(xn[index] - ProbesCoordinates[iprobe][0], yn[index] - ProbesCoordinates[iprobe][1], zn[index] - ProbesCoordinates[iprobe][2]) + TINY;
    }

    odistance = pow(1.0 / distance, r);
    numerator.push_back(odistance);
    interpolationIndexs.push_back(index);
}

LIB_EXPORT void Post_Probes::NecessaryVarForInterpolation(StructGrid *grid, int iprobe, int i, int j, int k, vector<RDouble> &numerator, vector<vector<int> > &interpolationIndexs, bool flag)
{
    vector<vector<RDouble> > ProbesCoordinates = grid->GetProbesCoordinates();
    RDouble3D &xcc = *(grid->GetCellCenterX());
    RDouble3D &ycc = *(grid->GetCellCenterY());
    RDouble3D &zcc = *(grid->GetCellCenterZ());

    RDouble3D &xn = * grid->GetStructX();
    RDouble3D &yn = * grid->GetStructY();
    RDouble3D &zn = * grid->GetStructZ();

    RDouble distance = 0;
    RDouble odistance = 0;
    vector<int> index;

    //! The range of r is from 0.5 to 3.0.
    const int r = 1;

    if (flag)
    {
        distance = DISTANCE(xcc(i, j, k) - ProbesCoordinates[iprobe][0], ycc(i, j, k) - ProbesCoordinates[iprobe][1], zcc(i, j, k) - ProbesCoordinates[iprobe][2]) + TINY;
    }
    else
    {
        distance = DISTANCE(xn(i, j, k) - ProbesCoordinates[iprobe][0], yn(i, j, k) - ProbesCoordinates[iprobe][1], zn(i, j, k) - ProbesCoordinates[iprobe][2]) + TINY;
    }

    odistance = pow(1.0 / distance, r);
    numerator.push_back(odistance);
    index.push_back(i);
    index.push_back(j);
    index.push_back(k);
    interpolationIndexs.push_back(index);
    index.resize(0);
}

LIB_EXPORT void Post_Probes::GetAllProbesVarPtr(RDouble **q_probes)
{
    int varCount = 0;

    for (set<int>::iterator varIter = probeVariables.begin(); varIter != probeVariables.end(); ++ varIter)
    {
        int variableType = *varIter;
        string varName = GetProbesVariableName(variableType);
        q_probes[varCount ++] = (RDouble *)GetProbeVarPtr(varName);
    }
}

void Post_Probes::FreeAllProbeVarPtr()
{
    for (set<int>::iterator varIter = probeVariables.begin(); varIter != probeVariables.end(); ++ varIter)
    {
        int variableType = *varIter;
        string varName = GetProbesVariableName(variableType);
        RDouble *q_probes = (RDouble *)GetProbeVarPtr(varName);
        if (q_probes)
        {
            delete [] q_probes;
            //DeleteProbeVarPtr(varName);
        }
    }
}

DumpToFileASCII::DumpToFileASCII(vector<DataContainer *> datalist_in, string filename_in)
{
    this->dataList = datalist_in;
    this->fileName = filename_in;

    for (unsigned int izone = 0; izone < dataList.size(); ++ izone)
    {
        dataList[izone]->moveto_begin();
    }
    nVariables = 0;
}

DumpToFileASCII::~DumpToFileASCII()
{

}

void DumpToFileASCII::InitProbesVariables()
{
    int nProbeVariables = 6;
    int probeVariables[100] = {0, 1, 2, 3, 4, 5};
    nProbeVariables = GlobalDataBase::GetIntParaFromDB("nProbeVariables");
    GlobalDataBase::GetData("probeVariables", probeVariables, PHINT, nProbeVariables);

    Post_Probes *postProbes = new Post_Probes(nProbeVariables, probeVariables);

    int nVariablesToDump = postProbes->GetNumberofProbeVariables();

    SetnVariablesToDump(nVariablesToDump);

    vector<string> title_tecplot;
    postProbes->GetASCIIFileHeader(title_tecplot);
    SetVariables(title_tecplot.front());

    FreePointer(postProbes);
}

void DumpToFileASCII::Run()
{
    InitProbesVariables();
    WriteFlowFile();
}

void DumpToFileASCII::WriteFlowFile()
{
    int dataMonitorType = GlobalDataBase::GetIntParaFromDB("dataMonitorType");

    if (dataMonitorType == PROBESMONITOR)
    {
        WriteProbesFlowFile();
    }
    else if (dataMonitorType == LINESMONITOR)
    {
        WriteLinesFlowFile();
    }
    else if (dataMonitorType == SURFACESMONITOR)
    {
        WriteSurfacesFlowFile();
    }
    else
    {
        TK_Exit::UnexpectedVarValue("dataMonitorType", dataMonitorType);
    }
}

void DumpToFileASCII::WriteProbesFlowFile()
{
    uint_t nzones = dataList.size();

    for (int izone = 0; izone < nzones; ++ izone)
    {
        DataContainer *cdata = dataList[izone];

        int ZoneProbesNumber;
        PHRead(cdata, &ZoneProbesNumber, 1);

        if (ZoneProbesNumber == 0) continue;

        vector<vector<RDouble> > ProbesCoordinates;
        vector<int> ProbesLineID, ProbesSurfaceID;
        ProbesCoordinates.resize(ZoneProbesNumber);
        ProbesLineID.resize(ZoneProbesNumber);
        ProbesSurfaceID.resize(ZoneProbesNumber);
        int nDim;
        RDouble **qProbes = NewPointer2<RDouble> (nVariables, ZoneProbesNumber);

        for (int iProbe = 0; iProbe < ZoneProbesNumber; ++ iProbe)
        {
            if (GetFilename() == "")
            {
                return;
            }
            ostringstream ossFileName;
            string finalFileName = GetFilename();

            PHRead(cdata, &nDim, 1);
            PHRead(cdata, &ProbesLineID[iProbe], 1);
            PHRead(cdata, &ProbesSurfaceID[iProbe], 1);
            ProbesCoordinates[iProbe].resize(nDim);

            for (int iDim = 0; iDim < nDim; ++ iDim)
            {
                PHRead(cdata, &ProbesCoordinates[iProbe][iDim], 1);
                finalFileName = PHSPACE::AddSymbolToFileName(finalFileName, '_', ProbesCoordinates[iProbe][iDim]);
            }

            if (PHMPI::IsParallelRun())
            {
                int fileID = PHMPI::GetFileIndexofCurrentProcessor();
                finalFileName = PHSPACE::AddSymbolToFileName(finalFileName, '_', fileID);
            }

            ossFileName << finalFileName;

            string probesFlowFileName = ossFileName.str();
            fstream probesFlowFile;
            OpenFile(probesFlowFile, probesFlowFileName, ios_base::out|ios_base::app);
            probesFlowFile << setiosflags(ios::left);
            probesFlowFile << setiosflags(ios::scientific);
            probesFlowFile << setprecision(10);

            if (IfFileEmpty(probesFlowFile))
            {
                probesFlowFile << variables << endl;
            }

            int wordWidth = 8;

            int outIterStep = GlobalDataBase::GetIntParaFromDB("outnstep");

            probesFlowFile << setw(wordWidth) << outIterStep << "	";

            for (int iVariable = 0; iVariable < nVariables; ++ iVariable)
            {
                PHRead(cdata, qProbes[iVariable][iProbe]);
                probesFlowFile << setw(wordWidth) << qProbes[iVariable][iProbe] << "	";
            }

            RDouble step_time = PHSPACE::TIME_SPACE::GetWallTime();
            probesFlowFile << setw(wordWidth) << step_time << endl;

            CloseFile(probesFlowFile);
        }

        DelPointer2( qProbes );
    }
}

void DumpToFileASCII::WriteLinesFlowFile()
{
    uint_t nzones = dataList.size();

    for (int izone = 0; izone < nzones; ++ izone)
    {
        DataContainer *cdata = dataList[izone];

        int ZoneProbesNumber;
        PHRead(cdata, &ZoneProbesNumber, 1);

        if (ZoneProbesNumber == 0) continue;

        vector<vector<RDouble> > ProbesCoordinates;
        vector<int> ProbesLineID, ProbesSurfaceID;
        ProbesCoordinates.resize(ZoneProbesNumber);
        ProbesLineID.resize(ZoneProbesNumber);
        ProbesSurfaceID.resize(ZoneProbesNumber);
        int nDim;
        RDouble **qProbes = NewPointer2<RDouble> (nVariables, ZoneProbesNumber);

        for (int iProbe = 0; iProbe < ZoneProbesNumber; ++ iProbe)
        {
            if (GetFilename() == "")
            {
                return;
            }
            ostringstream ossFileName;
            string finalFileName = GetFilename();

            PHRead(cdata, &nDim, 1);
            PHRead(cdata, &ProbesLineID[iProbe], 1);
            PHRead(cdata, &ProbesSurfaceID[iProbe], 1);

            finalFileName = PHSPACE::AddSymbolToFileName(finalFileName, "_", "line", ProbesLineID[iProbe]);

            if (PHMPI::IsParallelRun())
            {
                int fileID = PHMPI::GetFileIndexofCurrentProcessor();
                finalFileName = PHSPACE::AddSymbolToFileName(finalFileName, '_', fileID);
            }

            ossFileName << finalFileName;

            string linesFlowFileName = ossFileName.str();
            fstream linesFlowFile;
            OpenFile(linesFlowFile, linesFlowFileName, ios_base::out|ios_base::app);
            linesFlowFile << setiosflags(ios::left);
            linesFlowFile << setiosflags(ios::scientific);
            linesFlowFile << setprecision(10);

            if (IfFileEmpty(linesFlowFile))
            {
                linesFlowFile << variables << endl;
            }

            int wordWidth = 8;

            int outIterStep = GlobalDataBase::GetIntParaFromDB("outnstep");
            linesFlowFile << setw(wordWidth) << outIterStep << "	";

            ProbesCoordinates[iProbe].resize(nDim);
            for (int iDim = 0; iDim < nDim; ++ iDim)
            {
                PHRead(cdata, &ProbesCoordinates[iProbe][iDim], 1);
                linesFlowFile << setw(wordWidth) << ProbesCoordinates[iProbe][iDim] << "	";
            }

            for (int iVariable = 0; iVariable < nVariables; ++ iVariable)
            {
                PHRead(cdata, qProbes[iVariable][iProbe]);
                linesFlowFile << setw(wordWidth) << qProbes[iVariable][iProbe] << "	";
            }

            RDouble step_time = PHSPACE::TIME_SPACE::GetWallTime();
            linesFlowFile << setw(wordWidth) << step_time << endl;

            CloseFile(linesFlowFile);
        }

        DelPointer2(qProbes);
    }
}

void DumpToFileASCII::WriteSurfacesFlowFile()
{
    uint_t nzones = dataList.size();

    for (int izone = 0; izone < nzones; ++ izone)
    {
        DataContainer *cdata = dataList[izone];

        int ZoneProbesNumber;
        PHRead(cdata, &ZoneProbesNumber, 1);

        if (ZoneProbesNumber == 0) continue;

        vector<vector<RDouble> > ProbesCoordinates;
        vector<int> ProbesLineID, ProbesSurfaceID;
        ProbesCoordinates.resize(ZoneProbesNumber);
        ProbesLineID.resize(ZoneProbesNumber);
        ProbesSurfaceID.resize(ZoneProbesNumber);
        int nDim;
        RDouble **qProbes = NewPointer2<RDouble> (nVariables, ZoneProbesNumber);

        for (int iProbe = 0; iProbe < ZoneProbesNumber; ++ iProbe)
        {
            if (GetFilename() == "")
            {
                return;
            }
            ostringstream ossFileName;
            string finalFileName = GetFilename();

            PHRead(cdata, &nDim, 1);
            PHRead(cdata, &ProbesLineID[iProbe], 1);
            PHRead(cdata, &ProbesSurfaceID[iProbe], 1);

            finalFileName = PHSPACE::AddSymbolToFileName(finalFileName, "_", "surface", ProbesSurfaceID[iProbe]);

            if (PHMPI::IsParallelRun())
            {
                int fileID = PHMPI::GetFileIndexofCurrentProcessor();
                finalFileName = PHSPACE::AddSymbolToFileName(finalFileName, '_', fileID);
            }

            ossFileName << finalFileName;

            string linesFlowFileName = ossFileName.str();
            fstream linesFlowFile;
            OpenFile(linesFlowFile, linesFlowFileName, ios_base::out|ios_base::app);
            linesFlowFile << setiosflags(ios::left);
            linesFlowFile << setiosflags(ios::scientific);
            linesFlowFile << setprecision(10);

            if (IfFileEmpty(linesFlowFile))
            {
                linesFlowFile << variables << endl;
            }

            int wordWidth = 8;

            int outIterStep = GlobalDataBase::GetIntParaFromDB("outnstep");
            linesFlowFile << setw(wordWidth) << outIterStep << "	";

            ProbesCoordinates[iProbe].resize(nDim);
            for (int iDim = 0; iDim < nDim; ++ iDim)
            {
                PHRead(cdata, &ProbesCoordinates[iProbe][iDim], 1);
                linesFlowFile << setw(wordWidth) << ProbesCoordinates[iProbe][iDim] << "	";
            }

            for (int iVariable = 0; iVariable < nVariables; ++ iVariable)
            {
                PHRead(cdata, qProbes[iVariable][iProbe]);
                linesFlowFile << setw(wordWidth) << qProbes[iVariable][iProbe] << "	";
            }

            RDouble step_time = PHSPACE::TIME_SPACE::GetWallTime();
            linesFlowFile << setw(wordWidth) << step_time << endl;

            CloseFile(linesFlowFile);
        }

        DelPointer2( qProbes );
    }
}

}