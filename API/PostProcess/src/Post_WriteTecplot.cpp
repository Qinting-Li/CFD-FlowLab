#include "Post_WriteTecplot.h"
#include "PHIO.h"
#include "TK_Exit.h"
#include "TK_Parse.h"
#include "IO_FileName.h"
#include "TK_Log.h"
#include "Post_Visual.h"

namespace PHSPACE
{
Post_WriteTecplot::Post_WriteTecplot()
{
    this->visualDataList.resize(0);
}

Post_WriteTecplot::~Post_WriteTecplot()
{

}

void Post_WriteTecplot::Initialize()
{
    string visualfile = "tecflow.dat";
    GlobalDataBase::GetData("visualfile", &visualfile, PHSTRING, 1);

    string fname, fext;
    GetNameExt(visualfile, fname, fext, ".");

    VTKFileName = visualfile;
    visualFileName = visualfile;
    VTKvisual = false;

    if (flowType == AverageFlow)
    {
        visualFileName = AddSymbolToFileName(visualfile,"_Average");
    }

    if (flowType == AverageReynoldsStress)
    {
        visualFileName = AddSymbolToFileName(visualfile,"_ReynoldsStress");
    }

    string sentinelfilename = "results/sentinel";
    fstream sentinelfile;
    sentinelfile.open(sentinelfilename.c_str(), ios::in);

    if (!sentinelfile && !visualfile.empty() && (flowType == 0))
    {
        VTKFileName = visualfile + ".bak";
        VTKvisual = true;
    }

    int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
    if (isUnsteady)
    {
        int outnstep = GlobalDataBase::GetIntParaFromDB("outnstep");
        visualFileName = AddSymbolToFileName(visualFileName, '_', outnstep);
    }
}

void Post_WriteTecplot::StoreVisualizationData()
{
    using namespace PHMPI;
    int myid = GetCurrentProcessorID();
    int nZones = GetNumberofGlobalZones();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        DataContainer *cdata = new DataContainer();
        cdata->moveto_begin();

        int send_proc = GetZoneProcessorID(iZone);
        int recv_proc = GetZoneFileID(iZone);

        if (myid == send_proc)
        {
            StoreVisualizationData(iZone, cdata);
        }

        PH_Trade(cdata, send_proc, recv_proc, iZone);

        if (myid == recv_proc)
        {
            visualDataList.push_back(cdata);
        }
        else
        {
            delete cdata;
        }
    }
}

void Post_WriteTecplot::StoreVisualizationData(int zoneIndex, DataContainer *cdata)
{
    Grid *grid = flowFieldDataOnNode[zoneIndex]->GetGrid();
    if (grid->Type() == UNSTRUCTGRID)
    {
        int gridID = grid->GetZoneID();
        int dimension = grid->GetDim();
        int type = grid->Type();

        PHWrite(cdata, gridID);
        PHWrite(cdata, dimension);
        PHWrite(cdata, type);

        StoreUnsBoundaryVisualData(zoneIndex, cdata);

        if (WantVisualField(grid))
        {
            StoreUnsFieldVisualData(zoneIndex, cdata);
        }

        int TecioMission = Writecomplete;
        PHWrite(cdata, TecioMission);
    }
    else
    {
        int gridID = grid->GetZoneID();
        int dimension = grid->GetDim();
        int type = grid->Type();

        PHWrite(cdata, gridID);
        PHWrite(cdata, dimension);
        PHWrite(cdata, type);

        StoreStrBoundaryVisualData(zoneIndex, cdata);

        if (WantVisualField(grid))
        {
            StoreStrFieldVisualData(zoneIndex, cdata);
        }

        int TecioMission = Writecomplete;
        PHWrite(cdata, TecioMission);
    }
}

void Post_WriteTecplot::StoreStrBoundaryVisualData(int zoneIndex, DataContainer *cdata)
{
    Post_Visual *FieldData = flowFieldDataOnNode[zoneIndex];
    int nvarplot = FieldData->GetNumberofVisualVariables();
    int nVisualVariables = nvarplot;
    int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
    if (nchem)
    {
        int numberOfSpecies = GlobalDataBase::GetIntParaFromDB("numberOfSpecies");
        nVisualVariables = nVisualVariables + 2 * numberOfSpecies;
    }
    RDouble4D **qn = new RDouble4D *[nVisualVariables];
    FieldData->GetAllVisualNodeVarPtr(qn);

    int TecioMission = WriteBoundary;
    for (int iBoundary = 0; iBoundary < boundaryGrid.size(); ++ iBoundary)
    {
        if (boundaryGrid[iBoundary]->GetZoneID() != zoneIndex)
        {
            continue;
        }

        int bctype    = boundaryType[iBoundary];
        string bcName = boundaryName[iBoundary];

        PHWrite(cdata, TecioMission);
        PHWrite(cdata, bctype);
        cdata->write_string(bcName);

        int ni = StructGridCast(boundaryGrid[iBoundary])->GetNI();
        int nj = StructGridCast(boundaryGrid[iBoundary])->GetNJ();
        int nk = StructGridCast(boundaryGrid[iBoundary])->GetNK();

        PHWrite(cdata, ni);
        PHWrite(cdata, nj);
        PHWrite(cdata, nk);

        int nTotalNode = boundaryGrid[iBoundary]->GetNTotalNode();
        PHWrite(cdata, nTotalNode);

        RDouble *xBoundary = boundaryGrid[iBoundary]->GetX();
        RDouble *yBoundary = boundaryGrid[iBoundary]->GetY();
        RDouble *zBoundary = boundaryGrid[iBoundary]->GetZ();
        PHWrite(cdata, xBoundary, nTotalNode);
        PHWrite(cdata, yBoundary, nTotalNode);
        PHWrite(cdata, zBoundary, nTotalNode);

        for (int m = 0; m < nVisualVariables; ++ m)
        {
            RDouble *qBoundary = new RDouble[nTotalNode];

            for (int iNode = 0; iNode < nTotalNode; ++ iNode)
            {
                int originaI = nodeMap[iBoundary][0][iNode];
                int originaJ = nodeMap[iBoundary][1][iNode];
                int originaK = nodeMap[iBoundary][2][iNode];

                qBoundary[iNode] = (*qn[m])(originaI, originaJ, originaK, 0);
            }

            PHWrite(cdata, qBoundary, nTotalNode);
            delete []qBoundary;
        }
    }

    delete []qn;
}

void Post_WriteTecplot::StoreStrFieldVisualData(int zoneIndex, DataContainer *cdata)
{
    int nVisualVariables = flowFieldDataOnNode[zoneIndex]->GetNumberofVisualVariables();
    int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
    if (nchem)
    {
        int numberOfSpecies = GlobalDataBase::GetIntParaFromDB("numberOfSpecies");
        nVisualVariables = nVisualVariables + 2 * numberOfSpecies;
    }
    RDouble4D **qn = new RDouble4D *[nVisualVariables];
    flowFieldDataOnNode[zoneIndex]->GetAllVisualNodeVarPtr(qn);

    Grid *grid = flowFieldDataOnNode[zoneIndex]->GetGrid();
    StructGrid *strGrid = StructGridCast(grid);

    int TecioMission = WriteBlock;
    PHWrite(cdata, TecioMission);

    int ni = strGrid->GetNI();
    int nj = strGrid->GetNJ();
    int nk = strGrid->GetNK();

    PHWrite(cdata, ni);
    PHWrite(cdata, nj);
    PHWrite(cdata, nk);

    int nTotalNode = strGrid->GetNTotalNode();
    PHWrite(cdata, nTotalNode);

    RDouble *x = new RDouble[nTotalNode];
    RDouble *y = new RDouble[nTotalNode];
    RDouble *z = new RDouble[nTotalNode];

    RDouble3D &xx = *strGrid->GetStructX();
    RDouble3D &yy = *strGrid->GetStructY();
    RDouble3D &zz = *strGrid->GetStructZ();

    int count = 0;
    for (int k = 1; k <= nk; ++ k)
    {
        for (int j = 1; j <= nj; ++ j)
        {
            for (int i = 1; i <= ni; ++ i)
            {
                x[count] = xx(i, j, k);
                y[count] = yy(i, j, k);
                z[count] = zz(i, j, k);
                count ++;
            }
        }
    }

    PHWrite(cdata, x, nTotalNode);
    PHWrite(cdata, y, nTotalNode);
    PHWrite(cdata, z, nTotalNode);

    delete []x;
    delete []y;
    delete []z;

    count = 0;
    for (int m = 0; m < nVisualVariables; ++ m)
    {
        RDouble *qq = new RDouble[nTotalNode];
        for (int k = 1; k <= nk; ++ k)
        {
            for (int j = 1; j <= nj; ++ j)
            {
                for (int i = 1; i <= ni; ++ i)
                {
                    qq[count] = (*qn[m])(i, j, k, 0);
                    count ++;
                }
            }
        }
        PHWrite(cdata, qq, nTotalNode);
        delete []qq;
        count = 0;
    }

    delete []qn;
}

void Post_WriteTecplot::StoreUnsBoundaryVisualData(int zoneIndex, DataContainer *cdata)
{
    Post_Visual *FieldData = flowFieldDataOnNode[zoneIndex];
    int nvarplot         = FieldData->GetNumberofVisualVariables();
    int nVisualVariables = nvarplot;
    int nchem            = GlobalDataBase::GetIntParaFromDB("nchem");
    if (nchem)
    {
        int numberOfSpecies = GlobalDataBase::GetIntParaFromDB("numberOfSpecies");
        nVisualVariables = nVisualVariables + 2 * numberOfSpecies;
    }
    RDouble **qn = new RDouble *[nVisualVariables];
    FieldData->GetAllVisualNodeVarPtr(qn);

    UnstructGrid *grid = UnstructGridCast(FieldData->GetGrid());
    int *node_number_of_each_face = grid->GetNodeNumberOfEachFace();
    int *face2node = grid->GetFace2Node();
    int nBoundFace = grid->GetNBoundFace();
    int NumPts = grid->GetNTotalNode();

    RDouble **qFaceCenter = NewPointer2<RDouble>(nVisualVariables, nBoundFace);
    for (int m = 0; m < nVisualVariables; ++ m)
    {
        int nodeCount = 0;
        for (int iFace = 0; iFace < nBoundFace; ++ iFace)
        {
            RDouble value = 0.0;
            for (int iNode = 0; iNode < node_number_of_each_face[iFace]; ++ iNode)
            {
                int nodeIndex = face2node[nodeCount + iNode];
                value += qn[m][nodeIndex];
            }

            qFaceCenter[m][iFace] = value / node_number_of_each_face[iFace];
            nodeCount += node_number_of_each_face[iFace];
        }
    }

    int TecioMission = WriteBoundary;
    for (int iBoundary = 0; iBoundary < boundaryGrid.size(); ++ iBoundary)
    {
        if (boundaryGrid[iBoundary]->GetZoneID() != zoneIndex)
        {
            continue;
        }

        int bctype    = boundaryType[iBoundary];
        string bcName = boundaryName[iBoundary];

        PHWrite(cdata, TecioMission);
        PHWrite(cdata, bctype);
        cdata->write_string(bcName);

        int nTotalNode = boundaryGrid[iBoundary]->GetNTotalNode();
        int nTotalFace = boundaryGrid[iBoundary]->GetNTotalFace();
        int nTotalCell = boundaryGrid[iBoundary]->GetNTotalCell();
        PHWrite(cdata, nTotalNode);
        PHWrite(cdata, nTotalCell);
        PHWrite(cdata, nTotalFace);

        int TotalNumFaceNodes_Rect = nTotalCell * 4;
        PHWrite(cdata, &TotalNumFaceNodes_Rect, 1);

        int *node_number_of_each_cell = UnstructGridCast(boundaryGrid[iBoundary])->GetNodeNumberOfEachCell();
        PHWrite(cdata, node_number_of_each_cell, nTotalCell);

        PHWrite(cdata, nTotalNode);
        RDouble *xBoundary = boundaryGrid[iBoundary]->GetX();
        RDouble *yBoundary = boundaryGrid[iBoundary]->GetY();
        RDouble *zBoundary = boundaryGrid[iBoundary]->GetZ();
        PHWrite(cdata, xBoundary, nTotalNode);
        PHWrite(cdata, yBoundary, nTotalNode);
        PHWrite(cdata, zBoundary, nTotalNode);

        for (int m = 0; m < nVisualVariables; ++ m)
        {
            RDouble *qBoundary = new RDouble[nTotalNode];
            for (int iNode = 0; iNode < nTotalNode; ++ iNode)
            {
                int nodeIndex = nodeMap[iBoundary][0][iNode];

                if (nodeIndex < NumPts)
                {
                    qBoundary[iNode] = qn[m][nodeIndex];
                }
                else
                {
                    int faceIndex = nodeIndex - NumPts;
                    qBoundary[iNode] = qFaceCenter[m][faceIndex];
                }
            }

            PHWrite(cdata, qBoundary, nTotalNode);
            delete []qBoundary;
        }

        int *cell2node = UnstructGridCast(boundaryGrid[iBoundary])->GetCell2Node();
        PHWrite(cdata, cell2node, TotalNumFaceNodes_Rect);
    }

    delete []qn;
    DelPointer2(qFaceCenter);
}

void Post_WriteTecplot::StoreUnsFieldVisualData(int zoneIndex, DataContainer *cdata)
{
    int nVisualVariables = flowFieldDataOnNode[zoneIndex]->GetNumberofVisualVariables();
    int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
    if (nchem)
    {
        int numberOfSpecies = GlobalDataBase::GetIntParaFromDB("numberOfSpecies");
        nVisualVariables = nVisualVariables + 2 * numberOfSpecies;
    }
    RDouble **qn = new RDouble *[nVisualVariables];
    flowFieldDataOnNode[zoneIndex]->GetAllVisualNodeVarPtr(qn);

    UnstructGrid *grid = UnstructGridCast(flowFieldDataOnNode[zoneIndex]->GetGrid());
    int numberOfNodes = grid->GetNTotalNode();
    int numberOfCells = grid->GetNTotalCell();
    int numberOfFaces = grid->GetNTotalFace();

    int *faceNodeNumberContainer  = grid->GetNodeNumberOfEachFace();
    int *faceNodeIndexContainer   = grid->GetFace2Node();

    int *leftCellIndexContainer   = grid->GetLeftCellOfFace();
    int *rightCellIndexContainer  = grid->GetRightCellOfFace();

    int totalNumFaceNodes = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        totalNumFaceNodes += faceNodeNumberContainer[iFace];
    }

    int TecioMission = WriteBlock;
    PHWrite(cdata, TecioMission);
    PHWrite(cdata, numberOfNodes);
    PHWrite(cdata, numberOfCells);
    PHWrite(cdata, numberOfFaces);
    PHWrite(cdata, totalNumFaceNodes);
    PHWrite(cdata, faceNodeNumberContainer, numberOfFaces);

    RDouble *x = grid->GetX();
    RDouble *y = grid->GetY();
    RDouble *z = grid->GetZ();

    PHWrite(cdata, numberOfNodes);
    PHWrite(cdata, x, numberOfNodes);
    PHWrite(cdata, y, numberOfNodes);
    PHWrite(cdata, z, numberOfNodes);

    for (int iVarible = 0; iVarible < nVisualVariables; ++ iVarible)
    {
        RDouble *qq = qn[iVarible];
        PHWrite(cdata, qq, numberOfNodes);
    }

    PHWrite(cdata, faceNodeIndexContainer, totalNumFaceNodes);

    PHWrite(cdata, leftCellIndexContainer , numberOfFaces);
    PHWrite(cdata, rightCellIndexContainer, numberOfFaces);

    delete []qn;
}

void Post_WriteTecplot::WriteFile()
{
    int intervalStepPlot = GlobalDataBase::GetIntParaFromDB("intervalStepPlot");
    int outIterStep = GlobalDataBase::GetIntParaFromDB("outnstep");
    if (((outIterStep + 1) % intervalStepPlot == 0) && (outIterStep % intervalStepPlot != 0))
    {
        return;
    }

    using namespace PHMPI;
    int myid = GetCurrentProcessorID();

    if (VTKvisual == true && myid == server)
    {
        bool CharacteristicBoundary = false;
        DumpToVTK DumpToVTK(VTKDataList, VTKFileName, CharacteristicBoundary, this->flowType);
        DumpToVTK.Run();

        WriteSentinelFile();
    }

    if (visualDataList.size() > 0)
    {
        int visualfileType = GlobalDataBase::GetIntParaFromDB("visualfileType");
        if (visualfileType == Tecplot_Binary)
        {
#ifdef USE_TecplotLib
            DumpToTecio dumptotecio(visualDataList, visualFileName, this->flowType);
            dumptotecio.Run();
#else
            ostringstream WarningInfo;
            WarningInfo << "Warning : visualfileType = 0 is not right, dump flow data into ASCII file! \n";
            PrintToWindow(WarningInfo);

            DumpToTecplotASCII dumpToTecplotASCII(visualDataList, visualFileName, this->flowType);
            dumpToTecplotASCII.Run();
#endif
        }
        else if (visualfileType == Tecplot_ASCII)
        {
            DumpToTecplotASCII dumpToTecplotASCII(visualDataList, visualFileName, this->flowType);
            dumpToTecplotASCII.Run();
        }
        else
        {
            TK_Exit::UnexpectedVarValue("visualfileType", visualfileType);
        }
    }
}

void Post_WriteTecplot::ClearFieldData()
{
    for (int iData = 0; iData < flowFieldDataOnNode.size(); ++ iData)
    {
        if (flowFieldDataOnNode[iData])
        {
            FreePointer(flowFieldDataOnNode[iData]);
        }
    }

    for (unsigned int iZone = 0; iZone < VTKDataList.size(); ++ iZone)
    {
        delete VTKDataList[iZone];
    }
    VTKDataList.resize(0);

    for (unsigned int iZone = 0; iZone < visualDataList.size(); ++ iZone)
    {
        delete visualDataList[iZone];
    }
    visualDataList.resize(0);
}

}