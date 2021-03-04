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
    probeVariablesMap.insert(pair<int, string>(PROBE_DENSITY    , "density"));
    probeVariablesMap.insert(pair<int, string>(PROBE_U          , "u"));
    probeVariablesMap.insert(pair<int, string>(PROBE_V          , "v"));
    probeVariablesMap.insert(pair<int, string>(PROBE_W          , "w"));
    probeVariablesMap.insert(pair<int, string>(PROBE_PRESSURE   , "pressure"));
    probeVariablesMap.insert(pair<int, string>(PROBE_TEMPERATURE, "temperature"));

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

    ostringstream oss;
    oss << "title = \"Probes Flow Fields Data of PHengLEI\"" << endl;
    oss << "Variables =" << endl;
    oss << "\"iter\"" << endl;
    set<int>::iterator varIter;
    for (varIter = probeVariables.begin(); varIter != probeVariables.end(); ++ varIter)
    {
        int variableType = * varIter;
        string varName = probeVariablesMap[variableType];
        string varNameTemp = "\"" + varName + "\"";
        oss << varNameTemp << endl;
    }
    oss << "\"WallTime\"" << endl;
    title_tecplot.push_back(oss.str());
    return;
}

LIB_EXPORT void Post_Probes::UpdateProbeVarPtr(const string &name, RDouble *NearestCellCenterData)
{
    if (!grid)
    {
        TK_Exit::ExceptionExit("Post_Probe::UpdateProbesVarPtr needs Geometry infomation.");
    }
    UnstructGrid * unstructGrid = UnstructGridCast(grid);
    int zoneProbesNumber = unstructGrid->GetZoneProbesNumber();
    if (zoneProbesNumber == 0) return;

    vector<vector<RDouble> > probesCoordinates = unstructGrid->GetProbesCoordinates();
    vector<int> probesNearestCellID = unstructGrid->GetProbesNearestCellID();

    RDouble *q_probes = new RDouble [zoneProbesNumber];

    for (int izoneProbe = 0; izoneProbe < zoneProbesNumber; ++ izoneProbe)
    {
        int probeCell = probesNearestCellID[izoneProbe];
        q_probes[izoneProbe] = NearestCellCenterData[probeCell];
    }
    probeVariablesPtr->UpdateDataPtr(name, q_probes);
}

LIB_EXPORT void Post_Probes::UpdateProbeVarPtr(const string &name, RDouble4D &NearestCellCenterData, int index)
{
    if (!grid)
    {
        TK_Exit::ExceptionExit("Post_Probe::UpdateProbesVarPtr needs Geometry infomation.");
    }

    StructGrid * structGrid = StructGridCast(grid);

    int zoneProbesNumber = structGrid->GetZoneProbesNumber();
    if (zoneProbesNumber == 0) return;

    vector<vector<RDouble> > probesCoordinates = structGrid->GetProbesCoordinates();
    vector<int> probesNearestNI = structGrid->GetProbesNearestNI();
    vector<int> probesNearestNJ = structGrid->GetProbesNearestNJ();
    vector<int> probesNearestNK = structGrid->GetProbesNearestNK();

    RDouble *q_probes = new RDouble [zoneProbesNumber];

    for (int izoneProbe = 0; izoneProbe < zoneProbesNumber; ++ izoneProbe)
    {
        int probeNI = probesNearestNI[izoneProbe];
        int probeNJ = probesNearestNJ[izoneProbe];
        int probeNK = probesNearestNK[izoneProbe];
        q_probes[izoneProbe] = NearestCellCenterData(probeNI, probeNJ, probeNK, index);
    }

    probeVariablesPtr->UpdateDataPtr(name, q_probes);
}

LIB_EXPORT void Post_Probes::GetAllProbeVarPtr(RDouble **q_probes)
{
    int varCount = 0;

    for (set<int>::iterator varIter = probeVariables.begin(); varIter != probeVariables.end(); ++ varIter)
    {
        int variableType = * varIter;
        string varName = GetProbeVariableName(variableType);
        q_probes[varCount ++] = (RDouble *)GetProbeVarPtr(varName);
    }
}

void Post_Probes::FreeAllProbeVarPtr()
{
    for (set<int>::iterator varIter = probeVariables.begin(); varIter != probeVariables.end(); ++ varIter)
    {
        int variableType = * varIter;
        string varName = GetProbeVariableName(variableType);
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
    nVariable = 0;
}

DumpToFileASCII::~DumpToFileASCII()
{

}

void DumpToFileASCII::InitProbeVariables()
{
    int nProbeVariables = 6;
    int probeVariables[100] = {0, 1, 2, 3, 4, 5};
    nProbeVariables = GlobalDataBase::GetIntParaFromDB("nProbeVariables");
    GlobalDataBase::GetData("probeVariables", probeVariables, PHINT, nProbeVariables);

    Post_Probes * postProbes = new Post_Probes(nProbeVariables, probeVariables);

    int nVariablesToDump = postProbes->GetNumberofProbeVariables();

    SetnVariablesToDump(nVariablesToDump);

    vector<string> title_tecplot;
    postProbes->GetASCIIFileHeader(title_tecplot);
    SetVariables(title_tecplot.front());

    FreePointer(postProbes);
}

void DumpToFileASCII::Run()
{
    InitProbeVariables();
    WriteProbeFlowFile();
}

void DumpToFileASCII::WriteProbeFlowFile()
{
    uint_t nzones = dataList.size();
    for (int izone = 0; izone < nzones; ++ izone)
    {
        DataContainer *cdata = dataList[izone];

        int zoneProbesNumber;
        PHRead(cdata, &zoneProbesNumber, 1);

        if (zoneProbesNumber == 0) continue;

        vector<vector<RDouble> > probesCoordinates;
        probesCoordinates.resize(zoneProbesNumber);
        int nDim;
        RDouble **qProbes = NewPointer2<RDouble> (nVariable, zoneProbesNumber);

        for (int iProbe = 0; iProbe < zoneProbesNumber; ++ iProbe)
        {
            if (GetFilename() == "")
            {
                return;
            }

            ostringstream ossFileName;

            string finalFileName = GetFilename();
            PHRead(cdata, &nDim, 1);

            probesCoordinates[iProbe].resize(nDim);

            for (int iDim = 0; iDim < nDim; ++ iDim)
            {
                PHRead(cdata, &probesCoordinates[iProbe][iDim], 1);
                finalFileName = PHSPACE::AddSymbolToFileName(finalFileName, '_', probesCoordinates[iProbe][iDim]);
            }

            if (PHMPI::IsParallelRun())
            {
                int fileID = PHMPI::GetFileIndexofCurrentProcessor();
                finalFileName = PHSPACE::AddSymbolToFileName(finalFileName, '_', fileID);
            }

            ossFileName << finalFileName;

            string probeFlowFileName = ossFileName.str();
            fstream probeFlowFile;
            OpenFile(probeFlowFile, probeFlowFileName, ios_base::out|ios_base::app);
            probeFlowFile << setiosflags(ios::left);
            probeFlowFile << setiosflags(ios::scientific);
            probeFlowFile << setprecision(10);

            if (IfFileEmpty(probeFlowFile))
            {
                probeFlowFile << Variables << endl;
            }

            int wordWidth = 8;

            int outIterStep = GlobalDataBase::GetIntParaFromDB("outnstep");

            probeFlowFile << setw(wordWidth) << outIterStep << "    ";

            for (int iVariable = 0; iVariable < nVariable; ++ iVariable)
            {
                PHRead(cdata, qProbes[iVariable][iProbe]);
                probeFlowFile << setw(wordWidth) << qProbes[iVariable][iProbe] << "    ";
            }

            RDouble step_time = PHSPACE::TIME_SPACE::GetWallTime();
            probeFlowFile << setw(wordWidth) << step_time << endl;

            CloseFile(probeFlowFile);
        }

        DelPointer2( qProbes );
    }
}

}
