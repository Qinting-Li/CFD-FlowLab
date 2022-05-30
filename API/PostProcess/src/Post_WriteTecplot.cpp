#include "Post_WriteTecplot.h"
#include "IO_FileName.h"
#include "TK_Parse.h"
#include "PHIO.h"
#include "HyList.h"
#include "TK_Exit.h"
#include "Glb_Dimension.h"
#include "Geo_UnstructBC.h"
#include "Geo_StructBC.h"
#include "Post_Visual.h"
#include "Pointer.h"
#include "Solver.h"
#include "Geo_UnstructGrid.h"
#include "TK_Log.h"

namespace PHSPACE
{

LIB_EXPORT Post_WriteTecplot::Post_WriteTecplot(int flowType)
{
    this->flowType = flowType;
}

LIB_EXPORT Post_WriteTecplot::~Post_WriteTecplot()
{

}

LIB_EXPORT void Post_WriteTecplot::Run(bool visualAllLevel)
{
    int nMGLevel = 1;
    if (visualAllLevel)
    {
        GlobalDataBase::GetData("nMGLevel", &nMGLevel, PHINT, 1);
    }

    for (int iLevel = 0; iLevel < nMGLevel; ++ iLevel)
    {
        string visualfile = "tecflow.dat";
        GlobalDataBase::GetData("visualfile", &visualfile, PHSTRING, 1);

        string fname, fext;
        GetNameExt(visualfile, fname, fext, ".");

        ActionKey *actkey = new ActionKey();
        actkey->level = iLevel;
        actkey->filename = AddDigitToFileName(fname, fext, ".", iLevel);
        actkey->tecfilename = AddDigitToFileName(fname, fext, ".", iLevel);
        actkey->VTKvisual = false;

        if (flowType == AverageFlow)
        {
            actkey->tecfilename = AddDigitToFileName(fname, fext, "_Average.", iLevel);
        }

        if (flowType == AverageReynoldsStress)
        {
            actkey->tecfilename = AddDigitToFileName(fname, fext, "_ReynoldsStress.", iLevel);
        }

        string sentinelfilename = "results/sentinel";
        fstream sentinelfile;
        sentinelfile.open(sentinelfilename.c_str(), ios::in);

        if (!sentinelfile && !visualfile.empty() && (flowType == 0))
        {
            actkey->filename = visualfile + ".bak";
            actkey->VTKvisual = true;
        }

        int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
        if (isUnsteady)
        {
            int outnstep = PHSPACE::GlobalDataBase::GetIntParaFromDB("outnstep");
            actkey->tecfilename = PHSPACE::AddSymbolToFileName(actkey->tecfilename, '_', outnstep);
        }
        ServerWrite(actkey);
        PostWrite(actkey);
        delete actkey;

        this->Clear();
    }
}

void Post_WriteTecplot::ServerWrite(ActionKey *actkey)
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    int recv_proc;
    int myid;

    for (int iZone =  0; iZone < nZones; ++ iZone)
    {
        DataContainer *cdata = actkey->GetData();
        cdata->moveto_begin();

        DataContainer *tecdata = actkey->GetTecData();
        tecdata->moveto_begin();

        int send_proc = GetZoneProcessorID(iZone);
        //recv_proc = GetServerProcessorID();
        recv_proc = GetZoneFileID(iZone);

        myid = GetCurrentProcessorID();
        int tag = GetSendRecvTag(actkey, iZone);

        if (myid == send_proc)
        {
            Grid *grid = GetGrid(iZone, actkey->level);
            int nVisualVariables = GlobalDataBase::GetIntParaFromDB("nVisualVariables");
            int visualVariables[100];
            GlobalDataBase::GetData("visualVariables", visualVariables, PHINT, nVisualVariables);
            Post_Visual *postVisualization = new Post_Visual(nVisualVariables, visualVariables, grid, flowType);

            int Nsolver = GlobalSolvers::GetNumberOfSolvers(iZone);
            for (int iSolver = 0; iSolver < Nsolver; ++ iSolver)
            {
                GlobalSolvers::GetSolver(iZone, iSolver)->ComputePostVisualVariables(postVisualization);
            }
            
            StoreVisualizationVar(actkey, postVisualization);
            FreePointer(postVisualization);
        }

        PH_Trade(actkey, send_proc, recv_proc, tag);
        PH_Trade(tecdata, send_proc, recv_proc, tag);

        if (myid == recv_proc)
        {
            vtkdatalist.push_back(cdata);
            DataContainer *cdata0 = new DataContainer();
            actkey->SetData(cdata0);

            datalist.push_back(tecdata);
            DataContainer *cdata1 = new DataContainer();
            actkey->SetTecData(cdata1);
        }
    }
}

void Post_WriteTecplot::StoreVisualizationVar(ActionKey *actkey, Post_Visual *postVisualization)
{
    Grid *grid_in = postVisualization->GetGrid();

    if (grid_in->Type() == UNSTRUCTGRID)
    {
        UnstructGrid *grid = UnstructGridCast(grid_in);

        int nvarplot = postVisualization->GetNumberofVisualVariables();

        int nVisualVariables = nvarplot;
        int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
        if (nchem)
        {
            int numberOfSpecies = GlobalDataBase::GetIntParaFromDB("numberOfSpecies");
            nVisualVariables = nVisualVariables + numberOfSpecies;
        }

        RDouble **qn = new RDouble *[nVisualVariables];
        postVisualization->GetAllVisualNodeVarPtr(qn);

        if (actkey->VTKvisual == true)
        {
            int GridID = grid->GetGridID()->GetIndex();
            PHWrite(actkey->GetData(), GridID);

            int dimension = grid->GetDim();
            PHWrite(actkey->GetData(), dimension);

            int type = grid->Type();
            PHWrite(actkey->GetData(), &type, 1);
            
            BoundaryVTKVisualization(grid, actkey, qn, nVisualVariables);

            if (WantVisualField(grid))
            {
                FieldVisualization(grid, actkey, qn, nVisualVariables);
            }

            int TecioMission = TEC_SPACE::Writecomplete;
            PHWrite(actkey->GetData(), &TecioMission, 1);
        }

        if (actkey->tecfilename != "")
        {
            int GridID = grid->GetGridID()->GetIndex();
            PHWrite(actkey->GetTecData(), GridID);

            int dimension = grid->GetDim();
            PHWrite(actkey->GetTecData(), dimension);

            int type = grid->Type();
            PHWrite(actkey->GetTecData(), &type, 1);

            BoundaryVisualization(grid, actkey, qn, nVisualVariables);

            if (WantVisualField(grid))
            {
                SaveDataForTecio(grid, actkey, qn, nVisualVariables);
            }

            int TecioMission = TEC_SPACE::Writecomplete;
            PHWrite(actkey->GetTecData(), &TecioMission, 1);
        }

        int visualSlice = 0;
        GlobalDataBase::GetData("visualSlice", &visualSlice, PHINT, 1);

        if (grid->GetDim() == THREE_D && visualSlice == 1)
        {
            SliceVisualization(grid, actkey, qn, nvarplot);
        }

        delete [] qn;
    }
    else if (grid_in->Type() == STRUCTGRID)
    {
        StructGrid *grid = StructGridCast(grid_in);

        int nvarplot = postVisualization->GetNumberofVisualVariables();

        int nVisualVariables = nvarplot;
        int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
        if (nchem)
        {
            int numberOfSpecies = GlobalDataBase::GetIntParaFromDB("numberOfSpecies");
            nVisualVariables = nVisualVariables + numberOfSpecies;
        }

        RDouble4D **qn = new RDouble4D *[nVisualVariables];
        postVisualization->GetAllVisualNodeVarPtr(qn);

        if (actkey->VTKvisual == true)
        {
            int GridID = grid->GetGridID()->GetIndex();
            PHWrite(actkey->GetData(), GridID);

            int dimension = grid->GetDim();
            PHWrite(actkey->GetData(), dimension);

            int type = grid->Type();
            PHWrite(actkey->GetData(), &type, 1);

            BoundaryVTKVisualization(grid, actkey, qn, nVisualVariables);

            if (WantVisualField(grid))
            {
                FieldVisualization(grid, actkey, qn, nVisualVariables);
            }

            int TecioMission = TEC_SPACE::Writecomplete;
            PHWrite(actkey->GetData(), &TecioMission, 1);
        }

        if (actkey->tecfilename != "")
        {
            int GridID = grid->GetGridID()->GetIndex();
            PHWrite(actkey->GetTecData(), GridID);

            int dimension = grid->GetDim();
            PHWrite(actkey->GetTecData(), dimension);

            int type = grid->Type();
            PHWrite(actkey->GetTecData(), &type, 1);

            BoundaryVisualization(grid, actkey, qn, nVisualVariables);

            if (WantVisualField(grid))
            {
                SaveDataForTecio(grid, actkey, qn, nVisualVariables);
            }

            int TecioMission = TEC_SPACE::Writecomplete;
            PHWrite(actkey->GetTecData(), &TecioMission, 1);
        }

        delete [] qn;
    }
    else
    {
        TK_Exit::ExceptionExit("Error!!! Task_GenerateVisualFile::Visualization grid type not right!");
    }
}

void Post_WriteTecplot::PostWrite(ActionKey *actkey)
{
    using namespace PHMPI;
    int myid = GetCurrentProcessorID();

    if (actkey->VTKvisual == true && myid == server)
    {
        bool CharacteristicBoundary = false;
        DumpToVTK DumpToVTK(vtkdatalist, actkey->filename, CharacteristicBoundary, this->flowType);
        DumpToVTK.Run();

        WriteSentinelFile(actkey);
    }

    if (datalist.size() > 0)
    {
        int visualfileType = GlobalDataBase::GetIntParaFromDB("visualfileType");
        if (visualfileType == 0)
        {
                ostringstream WarningInfo;
                WarningInfo << "Warning : visualfileType = 0 is not right, dump flow data into ASCII file! \n";
                PrintToWindow(WarningInfo);

                DumpToTecplotASCII dumpToTecplotASCII(datalist, actkey->tecfilename, this->flowType);
                dumpToTecplotASCII.Run();
        }
        else if (visualfileType == 1)
        {
            DumpToTecplotASCII dumpToTecplotASCII(datalist, actkey->tecfilename, this->flowType);
            dumpToTecplotASCII.Run();
        }
        else
        {
            TK_Exit::UnexpectedVarValue("visualfileType = ", visualfileType);
        }
    }

    for (unsigned int iZone = 0; iZone < vtkdatalist.size(); ++ iZone)
    {
        delete vtkdatalist[iZone];
    }

    for (unsigned int iZone = 0; iZone < datalist.size(); ++ iZone)
    {
        delete datalist[iZone];
    }
}

bool Post_WriteTecplot::WantVisualField(Grid *grid)
{
    int plotFieldType = 0;
    GlobalDataBase::GetData("plotFieldType", &plotFieldType, PHINT, 1);

    if (plotFieldType == TEC_SPACE::BlockVisual)
    {
        int dimension = GetDim();

        RDouble *pmin = grid->GetMinBox();
        RDouble *pmax = grid->GetMaxBox();

        RDouble lowerPlotFieldBox[3], upperPlotFieldBox[3];
        GlobalDataBase::GetData("lowerPlotFieldBox", &lowerPlotFieldBox, PHDOUBLE, 3);
        GlobalDataBase::GetData("upperPlotFieldBox", &upperPlotFieldBox, PHDOUBLE, 3);

        if ((pmin[0] > upperPlotFieldBox[0]) || (pmax[0] < lowerPlotFieldBox[0])) return false;

        if ((pmin[1] > upperPlotFieldBox[1]) || (pmax[1] < lowerPlotFieldBox[1])) return false;

        if (dimension == THREE_D)
        {
            if ((pmin[2] > upperPlotFieldBox[2]) || (pmax[2] < lowerPlotFieldBox[2])) return false;
        }
    }
    return grid->GetDim() == TWO_D || (grid->GetDim() == THREE_D && plotFieldType != TEC_SPACE::BoundaryVisual);
}

void Post_WriteTecplot::BoundaryVTKVisualization(Grid *grid_in, ActionKey *actkey, RDouble **qn, int nvarplot)
{
    UnstructGrid *grid = UnstructGridCast(grid_in);
    DataContainer *cdata = actkey->GetData();

    int nTotalNode = grid->GetNTotalNode();
    int nBoundFace = grid->GetNBoundFace();

    int *node_number_of_each_face = grid->GetNodeNumberOfEachFace();
    int *face2node = grid->GetFace2Node();

    RDouble *x = grid->GetX();
    RDouble *y = grid->GetY();
    RDouble *z = grid->GetZ();

    RDouble *xfc = grid->GetFaceCenterX();
    RDouble *yfc = grid->GetFaceCenterY();
    RDouble *zfc = grid->GetFaceCenterZ();

    RDouble **qFaceCenter = NewPointer2<RDouble>(nvarplot, nBoundFace);

    using namespace PHENGLEI;
    vector<vector<int> > face2nodelist(nBoundFace);
    vector<int> linkmap;

    set<pair<int, string> > bcnameMap;
    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();

    int nodeCount = 0;
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);

        string bcName = bcRegion->GetBCName();
        int bcType = bcRegion->GetBCType();

        bcnameMap.insert(pair<int, string>(bcType, bcName));

        for (int m = 0; m < nvarplot; ++ m)
        {
            RDouble value = 0.0;
            for (int iNode = 0; iNode < node_number_of_each_face[iFace]; ++ iNode)
            {
                int nodeIndex = face2node[nodeCount + iNode];
                value += qn[m][nodeIndex];
            }

            qFaceCenter[m][iFace] = value / node_number_of_each_face[iFace];
        }

        nodeCount += node_number_of_each_face[iFace];
    }

    int TecioMission;
    set<pair<int, string> >::iterator iter;
    for (iter = bcnameMap.begin(); iter != bcnameMap.end(); ++ iter)
    {
        if ((*iter).first == PHENGLEI::INTERFACE || (*iter).first < 0) continue;

        TecioMission = TEC_SPACE::WriteBoundary;
        PHWrite(cdata, &TecioMission, 1);

        int BCtype = (*iter).first;
        PHWrite(cdata, &BCtype, 1);

        string bcName = (*iter).second;
        cdata->write_string(bcName);

        face2nodelist.resize(0);
        linkmap.resize(0);

        GetFace2NodeList(grid, *iter, linkmap, face2nodelist);
        
        int NumPts = static_cast<int>(linkmap.size());
        int NumElements = static_cast<int>(face2nodelist.size());

        PHWrite(cdata, &NumPts, 1);
        PHWrite(cdata, &NumElements, 1);
        
        RDouble *xBoundary = new RDouble [NumPts];
        RDouble *yBoundary = new RDouble [NumPts];
        RDouble *zBoundary = new RDouble [NumPts];

        for (int iNode = 0; iNode < NumPts; ++ iNode)
        {
            int nodeIndex = linkmap[iNode];

            if (nodeIndex < nTotalNode)
            {
                xBoundary[iNode] = x[nodeIndex];
                yBoundary[iNode] = y[nodeIndex];
                zBoundary[iNode] = z[nodeIndex];
            }
            else
            {
                int faceIndex = nodeIndex - nTotalNode;

                xBoundary[iNode] = xfc[faceIndex];
                yBoundary[iNode] = yfc[faceIndex];
                zBoundary[iNode] = zfc[faceIndex];
            }
        }

        PHWrite(cdata, xBoundary, NumPts);
        PHWrite(cdata, yBoundary, NumPts);
        PHWrite(cdata, zBoundary, NumPts);

        delete [] xBoundary;
        delete [] yBoundary;
        delete [] zBoundary;

        for (int m = 0; m < nvarplot; ++ m)
        {
            RDouble *qBoundary = new RDouble [NumPts];
            
            for (int iNode = 0; iNode < NumPts; ++ iNode)
            {
                int nodeIndex = linkmap[iNode];
                
                if (nodeIndex < nTotalNode)
                {
                    qBoundary[iNode] = qn[m][nodeIndex];
                }
                else
                {
                    int faceIndex = nodeIndex - nTotalNode;
                    qBoundary[iNode] = qFaceCenter[m][faceIndex];
                }
            }

            PHWrite(cdata, qBoundary, NumPts);
            delete [] qBoundary;
        }

        int *cell2node = new int [NumElements * 4];
        int ncount = 0;
        for (int iCell = 0; iCell < NumElements; ++ iCell)
        {
            uint_t np = face2nodelist[iCell].size();
            int index0 = face2nodelist[iCell][0];
            for (int ip = 0; ip < np; ++ ip)
            {
                cell2node[ncount] = face2nodelist[iCell][ip] + 1;
                ncount ++;
            }

            for (uint_t ip = np; ip < 4; ++ ip)
            {
                cell2node[ncount] = index0 + 1;
                ncount ++;
            }
        }

        PHWrite(cdata, cell2node, NumElements * 4);
        delete [] cell2node;
    }

    DelPointer2(qFaceCenter);
}

void Post_WriteTecplot::BoundaryVisualization(Grid *grid_in, ActionKey *actkey, RDouble **qn, int nvarplot)
{
    UnstructGrid *grid = UnstructGridCast(grid_in);
    DataContainer *cdata = actkey->GetTecData();

    int nTotalNode = grid->GetNTotalNode();
    int nBoundFace = grid->GetNBoundFace();

    int *node_number_of_each_face = grid->GetNodeNumberOfEachFace();
    int *face2node = grid->GetFace2Node();

    RDouble *x = grid->GetX();
    RDouble *y = grid->GetY();
    RDouble *z = grid->GetZ();

    RDouble *xfc = grid->GetFaceCenterX();
    RDouble *yfc = grid->GetFaceCenterY();
    RDouble *zfc = grid->GetFaceCenterZ();

    RDouble **qFaceCenter = NewPointer2<RDouble>(nvarplot, nBoundFace);

    using namespace PHENGLEI;
    vector<vector<int> > face2nodelist(nBoundFace);
    vector<int> linkmap;

    set<pair<int, string> > bcnameMap;

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();

    int nodeCount = 0;
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);

        string bcName = bcRegion->GetBCName();
        int bcType = bcRegion->GetBCType();

        bcnameMap.insert(pair<int, string>(bcType, bcName));

        for (int m = 0; m < nvarplot; ++ m)
        {
            RDouble value = 0.0;
            for (int iNode = 0; iNode < node_number_of_each_face[iFace]; ++ iNode)
            {
                int nodeIndex = face2node[nodeCount + iNode];
                value += qn[m][nodeIndex];
            }

            qFaceCenter[m][iFace] = value / node_number_of_each_face[iFace];
        }

        nodeCount += node_number_of_each_face[iFace];
    }

    int TecioMission;
    set<pair<int, string> >::iterator iter;
    for (iter = bcnameMap.begin(); iter != bcnameMap.end(); ++ iter)
    {
        if ((*iter).first == PHENGLEI::INTERFACE || (*iter).first < 0) continue;

        face2nodelist.resize(0);
        linkmap.resize(0);

        GetFace2NodeList(grid, *iter, linkmap, face2nodelist);

        if (face2nodelist.size() == 0)
        {
            continue;
        }

        TecioMission = TEC_SPACE::WriteBoundary;
        PHWrite(cdata, &TecioMission, 1);

        int BCtype = (*iter).first;
        PHWrite(cdata, &BCtype, 1);

        string bcName = (*iter).second;
        cdata->write_string(bcName);

        /*std::ostringstream os_zone;
        os_zone << "\"" << grid->GetZoneID() << "BC=" << *iter << "\" ";
        string zone_title = os_zone.str();*/

        int NumPts      = static_cast<int>(linkmap.size());
        int NumElements = static_cast<int>(face2nodelist.size());
        int NumFaces    = NumElements * 4;
        PHWrite(cdata, &NumPts, 1);
        PHWrite(cdata, &NumElements, 1);
        PHWrite(cdata, &NumFaces, 1);

        int TotalNumFaceNodes_Rect = NumElements * 4;
        int *cell2node = new int [TotalNumFaceNodes_Rect];
        int *node_number_of_each_cell = new int [NumElements];

        int count = 0;
        for (std::size_t i = 0; i < face2nodelist.size(); ++ i)
        {
            node_number_of_each_cell[i] = 4;
            uint_t np = face2nodelist[i].size();
            int index0 = face2nodelist[i][0];
            for (int ip = 0; ip < np; ++ ip)
            {
                cell2node[count] = face2nodelist[i][ip];
                count ++;
            }

            for (uint_t ip = np; ip < 4; ++ ip)
            {
                cell2node[count] = index0;
                count ++;
            }
        }

        PHWrite(cdata, &TotalNumFaceNodes_Rect, 1);
        PHWrite(cdata, node_number_of_each_cell, NumElements);

        PHWrite(cdata, &NumPts, 1);

        RDouble *xx = new RDouble [NumPts];
        RDouble *yy = new RDouble [NumPts];
        RDouble *zz = new RDouble [NumPts];

        for (std::size_t iNode = 0; iNode < linkmap.size(); ++ iNode)
        {
            int nodeIndex = linkmap[iNode];

            if (nodeIndex < nTotalNode)
            {
                xx[iNode] = x[nodeIndex];
                yy[iNode] = y[nodeIndex];
                zz[iNode] = z[nodeIndex];
            }
            else
            {
                int faceIndex = nodeIndex - nTotalNode;

                xx[iNode] = xfc[faceIndex];
                yy[iNode] = yfc[faceIndex];
                zz[iNode] = zfc[faceIndex];
            }
        }

        PHWrite(cdata, xx, NumPts);
        PHWrite(cdata, yy, NumPts);
        PHWrite(cdata, zz, NumPts);

        for (int m = 0; m < nvarplot; ++ m)
        {
            RDouble *qBoundary = new RDouble [NumPts];

            for (int iNode = 0; iNode < NumPts; ++ iNode)
            {
                int nodeIndex = linkmap[iNode];
                
                if (nodeIndex < nTotalNode)
                {
                    qBoundary[iNode] = qn[m][nodeIndex];
                }
                else
                {
                    int faceIndex = nodeIndex - nTotalNode;
                    qBoundary[iNode] = qFaceCenter[m][faceIndex];
                }
            }

            PHWrite(cdata, qBoundary, NumPts);
            delete [] qBoundary;
        }

        PHWrite(cdata, cell2node, TotalNumFaceNodes_Rect);

        delete [] cell2node;
        delete [] node_number_of_each_cell;
        delete [] xx;
        delete [] yy;
        delete [] zz;
    }

    DelPointer2(qFaceCenter);
}

void Post_WriteTecplot::BoundaryVisualization(Grid *grid_in, ActionKey *actkey, RDouble4D **qn, int nvarplot)
{
    StructGrid *grid = StructGridCast(grid_in);
    DataContainer *cdata = actkey->GetTecData();
    
    int s_st[3] = {0}, s_ed[3] = {0};

    RDouble3D &xx = * grid->GetStructX();
    RDouble3D &yy = * grid->GetStructY();
    RDouble3D &zz = * grid->GetStructZ();

    StructBCSet *structBCSet = grid->GetStructBCSet();
    int nBCRegion = structBCSet->GetnBCRegion();

    int TecioMission;

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        StructBC *bcregion = structBCSet->GetBCRegion(iBCRegion);
        int bctype = bcregion->GetBCType();

        if (bctype == PHENGLEI::INTERFACE || bctype < 0) continue;

        TecioMission = TEC_SPACE::WriteBoundary;
        PHWrite(cdata, &TecioMission, 1);

        PHWrite(cdata, &bctype, 1);

        string bcName = bcregion->GetBCName();
        cdata->write_string(bcName);

        int *s_lr3d = bcregion->GetFaceDirectionIndex();
        int nDim = GetDim();
        
        for (int m = 0; m < nDim; ++ m)
        //for (int m = 0; m < 3; ++ m)
        {
            s_st[m] = bcregion->GetStartPoint(m);
            s_ed[m] = bcregion->GetEndPoint(m);
            if (bcregion->GetFaceDirection() == m)
            {
                if (s_st[m] > 1 || s_lr3d[m] == 1)
                {
                    s_st[m] += 1;
                    s_ed[m] += 1;
                }
            }
            else
            {
                s_ed[m] += 1;
            }
        }

        int ist = s_st[0];
        int ied = s_ed[0];
        int jst = s_st[1];
        int jed = s_ed[1];
        int kst = s_st[2];
        int ked = s_ed[2];

        Range I(ist, ied);
        Range J(jst, jed);
        Range K(kst, ked);

        if (nDim != THREE_D)
        {
            K.setRange(1, 1);
            kst = 1;
            ked = 1;
        }

        Int3D newindex(I, J, K, fortranArray);

        int count = 0;
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    newindex(i, j, k) = count;
                    ++ count;
                }
            }
        }

        int ii = ied - ist + 1;
        int jj = jed - jst + 1;
        int kk = ked - kst + 1;

        PHWrite(cdata, &ii, 1);
        PHWrite(cdata, &jj, 1);
        PHWrite(cdata, &kk, 1);

        int nTotalNode = count;
        PHWrite(cdata, &nTotalNode, 1);

        RDouble *x = new RDouble [nTotalNode];
        RDouble *y = new RDouble [nTotalNode];
        RDouble *z = new RDouble [nTotalNode];

        count = 0;
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
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

        delete [] x;
        delete [] y;
        delete [] z;

        count = 0;
        for (int n = 0; n < nvarplot; ++ n)
        {
            RDouble *qq = new RDouble [nTotalNode];
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    for (int i = ist; i <= ied; ++ i)
                    {
                        qq[count] = (*qn[n])(i, j, k, 0);
                        count ++;
                    }
                }
            }

            PHWrite(cdata, qq, nTotalNode);
            delete [] qq;
            count = 0;
        }
    }
}

void Post_WriteTecplot::BoundaryVTKVisualization(Grid *grid_in, ActionKey *actkey, RDouble4D **qn, int nvarplot)
{
    StructGrid *grid = StructGridCast(grid_in);
    DataContainer *cdata = actkey->GetData();

    int s_st[3] = {0}, s_ed[3] = {0};
    int nk = grid->GetNK();

    RDouble3D &xx = * grid->GetStructX();
    RDouble3D &yy = * grid->GetStructY();
    RDouble3D &zz = * grid->GetStructZ();
    
    StructBCSet *structBCSet = grid->GetStructBCSet();
    int nBCRegion = structBCSet->GetnBCRegion();

    int TecioMission;
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        StructBC *bcregion = structBCSet->GetBCRegion(iBCRegion);
        int bctype = bcregion->GetBCType();
        if (bctype == PHENGLEI::INTERFACE || bctype < 0) continue;

        TecioMission = TEC_SPACE::WriteBoundary;
        PHWrite(cdata, &TecioMission, 1);

        PHWrite(cdata, &bctype, 1);

        string bcName = bcregion->GetBCName();
        cdata->write_string(bcName);

        int *s_lr3d = bcregion->GetFaceDirectionIndex();
        int nDim = GetDim();
        
        for (int m = 0; m < nDim; ++ m)
        //for (int m = 0; m < 3; ++ m)
        {
            s_st[m] = bcregion->GetStartPoint(m);
            s_ed[m] = bcregion->GetEndPoint(m);
            if (bcregion->GetFaceDirection() == m)
            {
                if (s_st[m] > 1 || s_lr3d[m] == 1)
                {
                    s_st[m] += 1;
                    s_ed[m] += 1;
                }
            }
            else
            {
                s_ed[m] += 1;
            }
        }

        int ist = s_st[0];
        int ied = s_ed[0];
        int jst = s_st[1];
        int jed = s_ed[1];
        int kst = s_st[2];
        int ked = s_ed[2];

        Range I(ist, ied);
        Range J(jst, jed);
        Range K(kst, ked);

        if (nDim != THREE_D)
        {
            K.setRange(1, 1);
            kst = 1;
            ked = 1;
        }

        Int3D newindex(I, J, K, fortranArray);

        int nsurf = bcregion->GetFaceDirection() + 1;
        //int lr = bcregion->GetFaceLeftOrRightIndex();

        int count = 0;
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    newindex(i, j, k) = count;
                    ++ count;
                }
            }
        }

        int nTotalNode = count;
        int nTotalCell = 1;

        for (int m = 0; m < 3; ++ m)
        {
            nTotalCell *= (bcregion->GetEndPoint(m) - bcregion->GetStartPoint(m) + 1);
        }

        PHWrite(cdata, &nTotalNode, 1);
        PHWrite(cdata, &nTotalCell, 1);

        RDouble *xBoundary = new RDouble [nTotalNode];
        RDouble *yBoundary = new RDouble [nTotalNode];
        RDouble *zBoundary = new RDouble [nTotalNode];

        RDouble **qBoundary = new RDouble *[nvarplot];
        for (int m = 0; m < nvarplot; ++ m)
        {
            qBoundary[m] = new RDouble [nTotalNode];
        }

        int iNode = 0;
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    xBoundary[iNode] = xx(i, j, k);
                    yBoundary[iNode] = yy(i, j, k);
                    zBoundary[iNode] = zz(i, j, k);

                    for (int m = 0; m < nvarplot; ++ m)
                    {
                        qBoundary[m][iNode] = (*qn[m])(i, j, k, 0);
                    }

                    ++ iNode;
                }
            }
        }

        PHWrite(cdata, xBoundary, nTotalNode);
        PHWrite(cdata, yBoundary, nTotalNode);
        PHWrite(cdata, zBoundary, nTotalNode);

        delete [] xBoundary;
        delete [] yBoundary;
        delete [] zBoundary;

        for (int m = 0; m < nvarplot; ++ m)
        {
            PHWrite(cdata, qBoundary[m], nTotalNode);
            delete [] qBoundary[m];
        }

        delete [] qBoundary;

        int il1, jl1, kl1;
        GetNsurfIndex(nsurf, il1, jl1, kl1);

        il1 = 1 - il1;
        jl1 = 1 - jl1;
        kl1 = 1 - kl1;

        if (nk == 1) kl1 = 0;

        int *cell2node = new int [nTotalCell * 4];
        int ncount = 0;
        for (int k = kst; k <= ked - kl1; ++ k)
        {
            for (int j = jst; j <= jed - jl1; ++ j)
            {
                for (int i = ist; i <= ied - il1; ++ i)
                {
                    int p1, p2, p3, p4;

                    int il = i + il1;
                    int jl = j + jl1;
                    int kl = k + kl1;

                    if (nsurf == 1)
                    {
                        p1 = newindex(i, j , k) + 1;
                        p2 = newindex(i, jl, k) + 1;
                        p3 = newindex(i, jl, kl) + 1;
                        p4 = newindex(i, j , kl) + 1;
                    }
                    else if (nsurf == 2)
                    {
                        p1 = newindex(i , j, k) + 1;
                        p2 = newindex(il, j, k) + 1;
                        p3 = newindex(il, j, kl) + 1;
                        p4 = newindex(i , j, kl) + 1;
                    }
                    else
                    {
                        p1 = newindex(i , j , k) + 1;
                        p2 = newindex(il, j , k) + 1;
                        p3 = newindex(il, jl, k) + 1;
                        p4 = newindex(i , jl, k) + 1;
                    }

                    cell2node[ncount++] = p1;
                    cell2node[ncount++] = p2;
                    cell2node[ncount++] = p3;
                    cell2node[ncount++] = p4;
                }
            }
        }

        PHWrite(cdata, cell2node, nTotalCell * 4);
        delete [] cell2node;
    }
}

void Post_WriteTecplot::FieldVisualizationForVTK(Grid *grid_in, ActionKey *actkey, RDouble **qn, int nvarplot)
{
    int TecioMission = TEC_SPACE::WriteBlock;
    DataContainer *cdata = actkey->GetData();
    PHWrite(cdata, &TecioMission, 1);

    UnstructGrid *grid = UnstructGridCast(grid_in);

    int nTotalNode = grid->GetNTotalNode();
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotalFace = grid->GetNTotalFace();

    int *face2node = grid->GetFace2Node();
    int *left_cell_of_face = grid->GetLeftCellOfFace();
    int *right_cell_of_face = grid->GetRightCellOfFace();

    RDouble *x = grid->GetX();
    RDouble *y = grid->GetY();
    RDouble *z = grid->GetZ();

    //! Find out cellToFace.
    HyList<int> * cell2face = new HyList <int> [nTotalCell];
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        cell2face[iCell].SetAverNnode(4);
    }
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int le = left_cell_of_face[iFace];
        cell2face[le].insert(iFace);
    }
    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        int le = left_cell_of_face[iFace];
        int re = right_cell_of_face[iFace];
        cell2face[le].insert(iFace);
        cell2face[re].insert(iFace);
    }

    //! Find out cellToNode.
    HyList<int> *cell2node = new HyList <int> [nTotalCell];
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        cell2node[iCell].SetAverNnode(4);
    }

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        int face0 = cell2face[iCell][0];
        int ip1 = face2node[2 * face0];
        int ip2 = face2node[2 * face0 + 1];
        cell2node[iCell].insert(ip1);
        cell2node[iCell].insert(ip2);

        int targetNode = ip2;
        int checkNode  = ip1;
        
        int nFaceOfCell = cell2face[iCell].size();
        int isClockWise = false;
        while (! isClockWise)
        {
            for (int iFace = 1; iFace < nFaceOfCell; ++ iFace)
            {
                int face  = cell2face[iCell][iFace];
                int node1 = face2node[2 * face];
                int node2 = face2node[2 * face + 1];

                if (node1 == targetNode)
                {
                    if (node2 != checkNode)
                    {
                        cell2node[iCell].insert(node2);
                        targetNode = node2;
                    }
                    else
                    {
                        isClockWise = true;
                        break;
                    }
                }
                else if (node2 == targetNode)
                {
                    if (node1 != checkNode)
                    {
                        cell2node[iCell].insert(node1);
                        targetNode = node1;
                    }
                    else
                    {
                        isClockWise = true;
                        break;
                    }
                }
            }    //! Iterate over all faces of cell.
        }        //! Judge whether the nodes of cell have searched closed.
    }            //! iCell.

    PHWrite(cdata, &nTotalNode, 1);
    PHWrite(cdata, &nTotalCell, 1);

    PHWrite(cdata, x, nTotalNode);
    PHWrite(cdata, y, nTotalNode);
    PHWrite(cdata, z, nTotalNode);

    for (int m = 0; m < nvarplot; ++ m)
    {
        RDouble *qq = qn[m];
        PHWrite(cdata, qq, nTotalNode);
    }

    int *cell2node_vtk = new int [nTotalCell * 4];
    int ncount = 0;
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        if (cell2node[iCell].size() == 3)
        {
            for (int m = 0; m < 3; ++ m)
            {
                cell2node_vtk[ncount] = cell2node[iCell].GetData(m) + 1;
                ncount ++;
            }
            cell2node_vtk[ncount] = cell2node[iCell].GetData(2) + 1;
            ncount ++;
        }
        else if (cell2node[iCell].size() == 4)
        {
            for (int m = 0; m < 4; ++ m)
            {
                cell2node_vtk[ncount] = cell2node[iCell].GetData(m) + 1;
                ncount ++;
            }
        }
        else
        {
            cout << "Error: this function only support triangle and quad !\n";
            exit(0);
        }
    }

    PHWrite(cdata, cell2node_vtk, nTotalCell * 4);

    delete [] cell2node_vtk;
    delete [] cell2node;
    delete [] cell2face;
}

void Post_WriteTecplot::FieldVisualization(Grid *grid_in, ActionKey *actkey, RDouble **qn, int nvarplot)
{
    if (GetDim() == PHSPACE::TWO_D)
    {
        FieldVisualizationForVTK(grid_in, actkey, qn, nvarplot);
        return;
    }
}

void Post_WriteTecplot::FieldVisualization(Grid *grid_in, ActionKey *actkey, RDouble4D **qn, int nvarplot)
{
    DataContainer *cdata = actkey->GetData();

    StructGrid *grid = StructGridCast(grid_in);

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    RDouble3D &xx = * grid->GetStructX();
    RDouble3D &yy = * grid->GetStructY();
    RDouble3D &zz = * grid->GetStructZ();

    int nl = 5;
    GlobalDataBase::GetData("nl", &nl, PHINT, 1);

    int nchem = 0;
    GlobalDataBase::GetData("nchem", &nchem, PHINT, 1);

    RDouble refGama = 1.4;
    GlobalDataBase::GetData("refGama", &refGama, PHDOUBLE, 1);

    int nTotalCell = grid->GetNTotalCell();
    int nTotalNode = grid->GetNTotalNode();

    if (nk == 1)
    {
        int TecioMission = TEC_SPACE::WriteBlock;
        PHWrite(cdata, &TecioMission, 1);

        PHWrite(cdata, &nTotalNode, 1);
        PHWrite(cdata, &nTotalCell, 1);

        RDouble *xBoundary = new RDouble [nTotalNode];
        RDouble *yBoundary = new RDouble [nTotalNode];
        RDouble *zBoundary = new RDouble [nTotalNode];

        RDouble **qBoundary = new RDouble *[nvarplot];
        for (int m = 0; m < nvarplot; ++ m)
        {
            qBoundary[m] = new RDouble [nTotalNode];
        }

        int ist, ied, jst, jed, kst, ked;
        grid->GetNodeIterationIndex(ist, ied, jst, jed, kst, ked);

        int iNode = 0;
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    xBoundary[iNode] = xx(i, j, k);
                    yBoundary[iNode] = yy(i, j, k);
                    zBoundary[iNode] = zz(i, j, k);

                    for (int m = 0; m < nvarplot; ++ m)
                    {
                        qBoundary[m][iNode] = (*qn[m])(i, j, k, 0);
                    }

                    iNode ++ ;
                }
            }
        }

        PHWrite(cdata, xBoundary, nTotalNode);
        PHWrite(cdata, yBoundary, nTotalNode);
        PHWrite(cdata, zBoundary, nTotalNode);

        delete [] xBoundary;
        delete [] yBoundary;
        delete [] zBoundary;

        for (int m = 0; m < nvarplot; ++ m)
        {
            PHWrite(cdata, qBoundary[m], nTotalNode);
            delete [] qBoundary[m];
        }

        delete [] qBoundary;

        int *cell2node = new int [nTotalCell * 4];
        int ncount = 0;
        for (int j = 1; j <= nj - 1; ++ j)
        {
            for (int i = 1; i <= ni - 1; ++ i)
            {
                cell2node[ncount++] = i     + (j - 1) * ni;
                cell2node[ncount++] = i + 1 + (j - 1) * ni;
                cell2node[ncount++] = i + 1 + (j   ) * ni;
                cell2node[ncount++] = i     + (j   ) * ni;
            }
        }

        PHWrite(cdata, cell2node, nTotalCell * 4);
        delete [] cell2node;
    }
}

void Post_WriteTecplot::SaveDataForTecio(Grid *grid_in, ActionKey *actkey, RDouble **qn, int nvarplot)
{
    DataContainer *cdata = actkey->GetTecData();

    UnstructGrid *grid = UnstructGridCast(grid_in);

    //int * cellIBlank = grid->GetBlankIndex();
    int numberOfNodes = grid->GetNTotalNode();
    int numberOfCells = grid->GetNTotalCell();
    int numberOfFaces = grid->GetNTotalFace();

    int numberOfNodesInZone = 0;
    int numberOfFacesInZone = 0;
    int numberOfCellsInZone = 0;

    vector<int> keyNodes;
    vector<int> keyFaces;
    vector<int> nodeList;
    vector<int> cellList;

    keyNodes.resize(numberOfNodes);
    keyFaces.resize(numberOfFaces);
    nodeList.resize(numberOfNodes);
    cellList.resize(numberOfCells);

    SetField(keyNodes,  0, numberOfNodes);
    SetField(keyFaces,  0, numberOfFaces);
    SetField(nodeList, -1, numberOfNodes);
    SetField(cellList, -1, numberOfCells);

    int *keyActiveOfCells         = grid->GetBlankIndex();

    int **cellNodeIndexContainer  = grid->GetCell2NodeArray();
    int * cellNodeNumberContainer = grid->GetNodeNumberOfEachCell();

    int **cellFaceIndexContainer  = grid->GetCell2Face();
    int * cellFaceNumberContainer = grid->GetFaceNumberOfEachCell();

    int *faceNodeNumberContainer  = grid->GetNodeNumberOfEachFace();
    int *faceNodeIndexContainer   = grid->GetFace2Node();

    int *leftCellIndexContainer   = grid->GetLeftCellOfFace();
    int *rightCellIndexContainer  = grid->GetRightCellOfFace();

    for (int iCell = 0; iCell < numberOfCells; ++ iCell)
    {
        if (keyActiveOfCells[iCell] == 1)
        {
            cellList[iCell] = numberOfCellsInZone;

            for (int iNodes = 0; iNodes < cellNodeNumberContainer[iCell]; ++ iNodes)
            {
                keyNodes[cellNodeIndexContainer[iCell][iNodes]] = 1;
            }

            for (int iFaces = 0; iFaces < cellFaceNumberContainer[iCell]; ++ iFaces)
            {
                keyFaces[cellFaceIndexContainer[iCell][iFaces]] = 1;
            }
            numberOfCellsInZone ++;
        }
    }

    for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
    {
        if (keyNodes[iNode] == 1)
        {
            nodeList[iNode] = numberOfNodesInZone;
            numberOfNodesInZone ++;
        }
    }

    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        if (keyFaces[iFace] == 1)
        {
            numberOfFacesInZone ++;
        }
    }

    int totalNumFaceNodes = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        if (keyFaces[iFace] != 1) continue;
        totalNumFaceNodes += faceNodeNumberContainer[iFace];
    }

    if (numberOfNodesInZone == 0)
    {
        return;
    }

    int TecioMission = TEC_SPACE::WriteBlock;

    PHWrite(cdata, TecioMission);
    PHWrite(cdata, numberOfNodesInZone);
    PHWrite(cdata, numberOfCellsInZone);
    PHWrite(cdata, numberOfFacesInZone);
    PHWrite(cdata, totalNumFaceNodes);

    int count = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        if (keyFaces[iFace] != 1) continue;

        int faceNodeNumber = faceNodeNumberContainer[iFace];
        PHWrite(cdata, faceNodeNumber);
    }

    PHWrite(cdata, numberOfNodesInZone);

    RDouble *x = grid->GetX();
    RDouble *y = grid->GetY();
    RDouble *z = grid->GetZ();

    for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
    {
        if (keyNodes[iNode] != 1) continue;
        PHWrite(cdata, x[iNode]);
    }

    for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
    {
        if (keyNodes[iNode] != 1) continue;
        PHWrite(cdata, y[iNode]);
    }

    for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
    {
        if (keyNodes[iNode] != 1) continue;
        PHWrite(cdata, z[iNode]);
    }

    for (int iVarible = 0; iVarible < nvarplot; ++ iVarible)
    {
        RDouble *qq = qn[iVarible];
        for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
        {
            if (keyNodes[iNode] != 1) continue;
            PHWrite(cdata, qq[iNode]);
        }
    }

    count = 0;

    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        int faceNodeNumber = faceNodeNumberContainer[iFace];

        if (keyFaces[iFace] != 1)
        {
            count += faceNodeNumber;
            continue;
        }

        for (int iNode = 0; iNode < faceNodeNumberContainer[iFace]; ++ iNode)
        {
            int nodeIndex = nodeList[faceNodeIndexContainer[count ++]];
            PHWrite(cdata, nodeIndex);
        }
    }

    //! Left cell.
    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        if (keyFaces[iFace] != 1) continue;

        int cellIndex = cellList[leftCellIndexContainer[iFace]];
        PHWrite(cdata, cellIndex);
    }

    //! Right cell.
    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        int cellIndex;

        if (keyFaces[iFace] == 0) continue;

        int rightCell = rightCellIndexContainer[iFace];

        if (rightCell >= numberOfCells || rightCell < 0)
        {
            cellIndex = -1;
        }
        else
        {
            cellIndex = cellList[rightCell];
        }

        PHWrite(cdata, cellIndex);
    }
}

void Post_WriteTecplot::SaveDataForTecio(Grid *grid_in, ActionKey *actkey, RDouble4D **qn, int nl)
{
    DataContainer *cdata = actkey->GetTecData();

    StructGrid *grid = StructGridCast(grid_in);

    int TecioMission = TEC_SPACE::WriteBlock;
    PHWrite(cdata, &TecioMission, 1);

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    PHWrite(cdata, &ni, 1);
    PHWrite(cdata, &nj, 1);
    PHWrite(cdata, &nk, 1);

    int nTotalNode = grid->GetNTotalNode();
    PHWrite(cdata, &nTotalNode, 1);

    RDouble *x = new RDouble [nTotalNode];
    RDouble *y = new RDouble [nTotalNode];
    RDouble *z = new RDouble [nTotalNode];

    RDouble3D &xx = * grid->GetStructX();
    RDouble3D &yy = * grid->GetStructY();
    RDouble3D &zz = * grid->GetStructZ();

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

    delete [] x;
    delete [] y;
    delete [] z;

    count = 0;
    for (int m = 0; m < nl; ++ m)
    {
        RDouble *qq = new RDouble [nTotalNode];
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
        delete [] qq;
        count = 0;
    }
}

void Post_WriteTecplot::SliceVisualization(Grid *grid_in, ActionKey *actkey, RDouble **qn, int nvarplot)
{
    UnstructGrid *grid = UnstructGridCast(grid_in);
    DataContainer *cdata = actkey->GetTecData();

    int nBoundFace = grid->GetNBoundFace();

    int *node_number_of_each_face = grid->GetNodeNumberOfEachFace();
    int *face2node = grid->GetFace2Node();

    RDouble *x = grid->GetX();
    RDouble *y = grid->GetY();
    RDouble *z = grid->GetZ();

    vector<RDouble> xx;
    vector<RDouble> yy;
    vector<RDouble> zz;
    vector<vector<RDouble> > qq;
    qq.resize(nvarplot);

    int slice_dir;
    GlobalDataBase::GetData("slice_dir", &slice_dir, PHINT, 1);

    RDouble slice;
    GlobalDataBase::GetData("slice", &slice, PHDOUBLE, 1);

    RDouble *xyz = 0;
    if (slice_dir == X_DIR)
    {
        xyz = x;
    }
    else if (slice_dir == Y_DIR)
    {
        xyz = y;
    }
    else
    {
        xyz = z;
    }

    int ncount = 0;

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);
        int bcType = bcRegion->GetBCType();

        if (bcType != PHENGLEI::SOLID_SURFACE)
        {
            ncount += node_number_of_each_face[iFace];
            continue;
        }

        RDouble coor_min =   LARGE;
        RDouble coor_max = - LARGE;

        for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
        {
            int index = face2node[ncount + j];
            coor_min = MIN(xyz[index], coor_min);
            coor_max = MAX(xyz[index], coor_max);
        }

        if ((coor_min - slice) * (coor_max - slice) > 0.0)
        {
            ncount += node_number_of_each_face[iFace];
            continue;
        }

        for (int j = 0; j < node_number_of_each_face[iFace]; ++ j)
        {
            int Node_index0 = face2node[ncount + j % node_number_of_each_face[iFace]];
            int Node_index1 = face2node[ncount + (j + 1) % node_number_of_each_face[iFace]];

            if ((xyz[Node_index0] - slice) * (xyz[Node_index1] - slice) > 0.0) continue;

            if (slice_dir == X_DIR)
            {
                xx.push_back(slice);
            }
            else if (fabs(x[Node_index0] - x[Node_index1]) < SMALL)
            {
                xx.push_back(x[Node_index0]);
            }
            else
            {
                RDouble Node_x = (slice - xyz[Node_index0]) * (x[Node_index1] - x[Node_index0]) / (xyz[Node_index1] - xyz[Node_index0]) + x[Node_index0];
                xx.push_back(Node_x);
            }

            if (slice_dir == Y_DIR)
            {
                yy.push_back(slice);
            }
            else if (fabs(y[Node_index0] - y[Node_index1]) < SMALL)
            {
                yy.push_back(y[Node_index0]);
            }
            else
            {
                RDouble Node_y = (slice - xyz[Node_index0]) * (y[Node_index1] - y[Node_index0]) / (xyz[Node_index1] - xyz[Node_index0]) + y[Node_index0];
                yy.push_back(Node_y);
            }

            if (slice_dir == Z_DIR)
            {
                zz.push_back(slice);
            }
            else if (fabs(z[Node_index0] - z[Node_index1]) < SMALL)
            {
                zz.push_back(z[Node_index0]);
            }
            else
            {
                RDouble Node_z = (slice - xyz[Node_index0]) * (z[Node_index1] - z[Node_index0]) / (xyz[Node_index1] - xyz[Node_index0]) + z[Node_index0];
                zz.push_back(Node_z);
            }

            for (int m = 0; m < nvarplot; ++ m)
            {
                if (fabs(qn[m][Node_index0] - qn[m][Node_index1]) < SMALL)
                {
                    RDouble Node_q = qn[m][Node_index0];
                    qq[m].push_back(Node_q);
                }
                else
                {
                    RDouble Node_q = (slice - xyz[Node_index0]) * (qn[m][Node_index1] - qn[m][Node_index0]) / (xyz[Node_index1] - xyz[Node_index0]) + qn[m][Node_index0];
                    qq[m].push_back(Node_q);
                }
            }
        }

        ncount += node_number_of_each_face[iFace];
    }

    int NumPts = static_cast<int>(xx.size());
    PHWrite(cdata, &NumPts, 1);

    if (NumPts == 0) return;

    RDouble *xxx = &xx[0];
    PHWrite(cdata, xxx, NumPts);

    RDouble *yyy = &yy[0];
    PHWrite(cdata, yyy, NumPts);

    RDouble *zzz = &zz[0];
    PHWrite(cdata, zzz, NumPts);

    for (int m = 0; m < nvarplot; ++ m)
    {
        RDouble * q = &qq[m][0];
        PHWrite(cdata, q, NumPts);
    }
}

void Post_WriteTecplot::Clear()
{
    this->vtkdatalist.clear();
    this->datalist.clear();
}

}
