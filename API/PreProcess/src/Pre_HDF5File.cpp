#include "Pre_HDF5File.h"
#include "Math_BasisFunction.h"
#include "Glb_Dimension.h"
#include "IO_FileName.h"
#include "Zone.h"
#include "GridType.h"
#include "Geo_UnstructBC.h"
#include "PHIO.h"
#include "TK_Log.h"
#include "GlobalDataBase.h"
#include "TK_Parse.h"

namespace PHSPACE
{
using namespace PHMPI;

IO_HDF5Write::IO_HDF5Write(const string &gridFileName, Grid **grids, int nBlocks)
{
    this->gridFileName = gridFileName;
    this->grids = grids;
    this->nBlocks = nBlocks;

    gridsLocal = 0;
    nPartOfEachProcessor = 0;

    block_proc = 0;
    block_idx  = 0;
    block_type = 0;
    file_index = 0;
}

IO_HDF5Write::~IO_HDF5Write()
{
    if (gridsLocal)
    {
        delete [] gridsLocal;
        delete [] nPartOfEachProcessor;
    }

    if (block_proc)
    {
        delete [] block_proc;
        delete [] block_idx;
        delete [] block_type;
        delete [] file_index;
    }
}

hid_t IO_HDF5Write::CreateFile(const string &filename)
{
    hid_t file;
    file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    return file;
}

void IO_HDF5Write::Run()
{
    hid_t file;

    CreateTotalInfo();

    int myid = PHMPI::GetCurrentProcessorID();
    file = CreateFile(AddSymbolToFileName(gridFileName, "_", myid));

    WriteVersionInfo(file);
    WriteTotalInfo(file);

    int nBlocksLocal = nPartOfEachProcessor[myid];
    for (int izone = 0; izone < nBlocksLocal; ++ izone)
    {
        ostringstream oss;
        oss << "      Dumping zone " << izone << " ....\n";
        PrintToWindow(oss);

        if (block_type[izone] == STRUCTGRID)
        {
            WriteStructGrid(file, gridsLocal[izone]);
        }
        else
        {
            WriteUnstructGrid(file, gridsLocal[izone],izone);
        }
    }

    H5Fclose(file);
}

void IO_HDF5Write::WriteVersionInfo(hid_t loc_id)
{
    hid_t grpData;

    grpData = GetGroup(loc_id, "Version");

    RDouble Version = 2.0;
    CreateAndWriteData(grpData, "VersionData", 1, 1, PHDOUBLE, &Version);

    H5Gclose(grpData);
}

void IO_HDF5Write::CreateTotalInfo()
{
    int nBlocksLocal = 0;

    for (int iBlock = 0; iBlock < nBlocks; ++ iBlock)
    {
        if (grids[iBlock])
        {
            ++ nBlocksLocal;
        }
    }

    gridsLocal = new Grid *[nBlocksLocal]();

    nBlocksLocal = 0;
    for (int iBlock = 0; iBlock < nBlocks; ++ iBlock)
    {
        if (grids[iBlock])
        {
            gridsLocal[nBlocksLocal] = grids[iBlock];
            ++ nBlocksLocal;
        }
    }

    int nProc = PHMPI::GetNumberOfProcessor();
    int myid  = PHMPI::GetCurrentProcessorID();

    nPartOfEachProcessor = new int[nProc]();
    PH_Allgather(&nBlocksLocal, 1, nPartOfEachProcessor, 1);

    block_proc = new int[nBlocksLocal]();
    block_idx  = new int[nBlocksLocal]();
    block_type = new int[nBlocksLocal]();
    file_index = new int[nBlocksLocal]();

    int count = 0;
    for (int iProc = 0; iProc < myid; ++ iProc)
    {
        count += nPartOfEachProcessor[iProc];
    }

    for (int iBlock = 0; iBlock < nBlocksLocal; ++ iBlock)
    {
        block_proc[iBlock] = gridsLocal[iBlock]->GetGridID()->GetIndex();
        block_idx [iBlock] = iBlock + count;
        block_type[iBlock] = gridsLocal[iBlock]->Type();
        file_index[iBlock] = gridsLocal[iBlock]->GetFileIndex();
    }

    int *block_proc_inp = GetZoneProcessorID_INP();
    if (block_proc_inp)
    {
        for (int iBlock = 0; iBlock < nBlocksLocal; ++ iBlock)
        {
            block_proc[iBlock] = block_proc_inp[iBlock + count];
        }
    }

    /*if(parallelPartitionMethod == METIS_OPENMP_METHOD)
    {
        // For OpenMP case, reset the processor ID of sub-zone to the source grid's.
        int zoneID = gridSource->GetGridID()->GetIndex();
        int processorID = PHMPI::GetZoneProcessorID(zoneID);
        for (int izone = 0; izone < nBlocks; ++ izone)
        {
            block_proc[izone] = processorID;     
        }
    }*/
}

void IO_HDF5Write::WriteTotalInfo(hid_t loc_id)
{
    hid_t grpData;
    grpData = GetGroup(loc_id, "Information");

    int myid = PHMPI::GetCurrentProcessorID();
    int nBlocksLocal = nPartOfEachProcessor[myid];

    CreateAndWriteData(grpData, "nBlocks"   , 1, 1           , PHINT, &nBlocksLocal);
    CreateAndWriteData(grpData, "block_proc", 1, nBlocksLocal, PHINT, block_proc);
    CreateAndWriteData(grpData, "block_idx" , 1, nBlocksLocal, PHINT, block_idx);
    CreateAndWriteData(grpData, "block_type", 1, nBlocksLocal, PHINT, block_type);
    CreateAndWriteData(grpData, "file_index", 1, nBlocksLocal, PHINT, file_index);

    H5Gclose(grpData);
}

void IO_HDF5Write::WriteUnstructGrid(hid_t loc_id, Grid *gridIn, int iZone)
{
    hid_t grpData;
    string grpName;

    ostringstream oss;
    int gridID = block_idx[iZone];
    oss << "Grid-" << gridID;
    grpName = oss.str();
    grpData = GetGroup(loc_id, grpName);

    int iDimensions[3];
    iDimensions[0] = gridIn->GetNTotalNode();
    iDimensions[1] = gridIn->GetNTotalFace();
    iDimensions[2] = gridIn->GetNTotalCell();
    CreateAndWriteData(grpData, "iDimensions", 1, 3, PHINT, iDimensions);

    WriteGridCoordinates(grpData, gridIn);

    WriteFaceTopology(grpData, gridIn);

    WriteInterfaceInfo(grpData, gridIn);

    WriteCellTopology(grpData, gridIn);

    WriteUnstrBCName(grpData, gridIn);

    WriteInterpointInfo(grpData, gridIn);

    H5Gclose(grpData);
}

void IO_HDF5Write::WriteInterpointInfo(hid_t loc_id, Grid *gridIn)
{
    hid_t grpData;
    grpData = GetGroup(loc_id, "InterpointInfo");

    int nIPoint = 0;
    InterpointInformation *interpointInfo = gridIn->GetInterpointInfo();
    if (!interpointInfo)
    {
        CreateAndWriteData(grpData, "nIPoint", 1, 1, PHINT, &nIPoint);
        H5Gclose(grpData);
        return;
    }

    nIPoint = interpointInfo->GetNumberOfInterpoints();
    CreateAndWriteData(grpData, "nIPoint", 1, 1, PHINT, &nIPoint);

    int *interPoint2ZoneID       = interpointInfo->GetInterPoint2ZoneID();
    int *interPoint2InterPointID = interpointInfo->GetInterPoint2InterPointID();
    int *interPoint2GlobalPoint  = interpointInfo->GetInterPoint2GlobalPoint();
    int *cellNumberOfInterPoint  = interpointInfo->GetCellNumberOfInterPoint();
    int *totalZonesOfInterPoint  = interpointInfo->GetTotalZonesOfInterPoint();
    int *labelOfInterPoint       = interpointInfo->GetLabelOfInterPoint();

    CreateAndWriteData(grpData, "interPoint2ZoneID"      , 1, nIPoint, PHINT, interPoint2ZoneID);
    CreateAndWriteData(grpData, "interPoint2InterPointID", 1, nIPoint, PHINT, interPoint2InterPointID);
    CreateAndWriteData(grpData, "interPoint2GlobalPoint" , 1, nIPoint, PHINT, interPoint2GlobalPoint);
    CreateAndWriteData(grpData, "cellNumberOfInterPoint" , 1, nIPoint, PHINT, cellNumberOfInterPoint);
    CreateAndWriteData(grpData, "totalZonesOfInterPoint" , 1, nIPoint, PHINT, totalZonesOfInterPoint);
    CreateAndWriteData(grpData, "labelOfInterPoint"      , 1, nIPoint, PHINT, labelOfInterPoint);

    H5Gclose(grpData);
}

void IO_HDF5Write::WriteStructGrid(hid_t loc_id, Grid *gridIn)
{
    hid_t grpData;
    string grpName;

    ostringstream oss;
    int gridID = gridIn->GetZoneID();
    oss << "Grid-" << gridID;
    grpName = oss.str();
    grpData = GetGroup(loc_id, grpName);

    StructGrid *grid = StructGridCast(gridIn);

    int iDimensions[3];
    iDimensions[0] = grid->GetNI();
    iDimensions[1] = grid->GetNJ();
    iDimensions[2] = grid->GetNK();
    CreateAndWriteData(grpData, "iDimensions", 1, 3, PHINT, iDimensions);

    WriteGridCoordinates(grpData, gridIn);
    WriteBCRegion(grpData, gridIn);
    WriteInterfaceInfo(grpData, gridIn);
    WriteStrBCName(grpData, gridIn);
    WriteLnkInfo(grpData, gridIn, gridID);

    H5Gclose(grpData);
}

void IO_HDF5Write::WriteGridCoordinates(hid_t loc_id, Grid *gridIn)
{
    hid_t grpData;

    grpData = GetGroup(loc_id, "GridCoordinates");

    RDouble *x = gridIn->GetX();
    RDouble *y = gridIn->GetY();
    RDouble *z = gridIn->GetZ();
    int nTotalNode = gridIn->GetNTotalNode();

    int axisup = 1;
    GlobalDataBase::GetData("axisup", &axisup, PHINT, 1);
    if (axisup == 2)
    {
        gridIn->TurnZAxis2YAxis();
    }

    CreateAndWriteData(grpData, "CoordinateX", 1, nTotalNode, PHDOUBLE, x);
    CreateAndWriteData(grpData, "CoordinateY", 1, nTotalNode, PHDOUBLE, y);
    CreateAndWriteData(grpData, "CoordinateZ", 1, nTotalNode, PHDOUBLE, z);

    H5Gclose(grpData);
}

void IO_HDF5Write::WriteFaceTopology(hid_t loc_id, Grid *gridIn)
{
    hid_t grpData;

    grpData = GetGroup(loc_id, "FaceTopology");

    UnstructGrid *grid = UnstructGridCast(gridIn);

    int nTotalFace = grid->GetNTotalFace();
    int *nodeNumberOfEachFace = grid->GetNodeNumberOfEachFace();
    CreateAndWriteData(grpData, "nodeNumberOfEachFace", 1, nTotalFace, PHINT, nodeNumberOfEachFace);

    int faceCount = 0;
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        faceCount += nodeNumberOfEachFace[iFace];
    }

    int *face2Node = grid->GetFace2Node();
    CreateAndWriteData(grpData, "face2Node", 1, faceCount, PHINT, face2Node);

    int *leftCellOfFace  = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    CreateAndWriteData(grpData, "leftCellOfFace" , 1, nTotalFace, PHINT, leftCellOfFace);
    CreateAndWriteData(grpData, "rightCellOfFace", 1, nTotalFace, PHINT, rightCellOfFace);

    int nBoundFace = grid->GetNBoundFace();
    CreateAndWriteData(grpData, "nBoundFace", 1, 1, PHINT, &nBoundFace);

    UnstructBCSet *unstructBCRegion = grid->GetUnstructBCRegion();
    int nBCRegionUnstruct = unstructBCRegion->GetnBCRegion();
    int *tmpBcType = new int[nBoundFace];

    for (int iBCRegionUnstruct = 0; iBCRegionUnstruct < nBCRegionUnstruct; iBCRegionUnstruct++)
    {
        UnstructBC* bcRegion = unstructBCRegion->GetBCRegion(iBCRegionUnstruct);        
        vector<int>* faceIndex = bcRegion->GetFaceIndex();
        for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
        {
            // iFace is the face number in the set of faceIndex.
            int iFace = *iter;
            tmpBcType[iFace] = bcRegion->GetBCType();
        }     
    }

    /*int *tmpBcType = new int[nBoundFace];
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        UnstructBC *bcRegion = unstructBCRegion->GetBCRegion(bcRegionIDofBCFace[iFace]);
        tmpBcType[iFace] = bcRegion->GetBCType();
    }*/
    CreateAndWriteData(grpData, "BcType", 1, nBoundFace, PHINT, tmpBcType);

    delete [] tmpBcType;
    H5Gclose(grpData);
}

void IO_HDF5Write::WriteCellTopology(hid_t loc_id, Grid *gridIn)
{
    hid_t grpData;
    grpData = GetGroup(loc_id, "CellTopology");

    UnstructGrid *grid = UnstructGridCast(gridIn);

    int nTotalCell = grid->GetNTotalCell();
    int *nodeNumberOfEachCell = grid->GetNodeNumberOfEachCell();
    CreateAndWriteData(grpData, "nodeNumberOfEachCell", 1, nTotalCell, PHINT, nodeNumberOfEachCell);

    int *cell2Node = grid->GetCell2Node();
    int nodeCount = 0;
    for (int icell = 0; icell < nTotalCell; ++ icell)
    {
        nodeCount += nodeNumberOfEachCell[icell];
    }
    CreateAndWriteData(grpData, "cell2Node", 1, nodeCount, PHINT, cell2Node);

    H5Gclose(grpData);
}

void IO_HDF5Write::WriteUnstrBCName(hid_t loc_id, Grid *gridIn)
{
    hid_t grpData;

    grpData = GetGroup(loc_id, "BCName");

    UnstructGrid *grid = UnstructGridCast(gridIn);

    int nBoundFace = grid->GetNBoundFace();
    string *bcNameList = new string[nBoundFace];

    UnstructBCSet *unstructBCRegion = grid->GetUnstructBCRegion();
    int nBCRegionUnstruct = unstructBCRegion->GetnBCRegion();
    uint_long totalSize = 0;

    for (int iBCRegionUnstruct = 0; iBCRegionUnstruct < nBCRegionUnstruct; iBCRegionUnstruct++)
    {
        UnstructBC* bcRegion = unstructBCRegion->GetBCRegion(iBCRegionUnstruct);
        const string &bcName = bcRegion->GetBCName();

        vector<int>* faceIndex = bcRegion->GetFaceIndex();
        for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
        {
            // iFace is the face number in the set of faceIndex.
            int iFace = *iter;
            bcNameList[iFace] = bcName;

            totalSize += static_cast< uint_long > (bcName.size());
            totalSize += 1;
        }

    }

    char *bcNameChar = new char[totalSize];
    unsigned int count = 0;
    for (int iFace = 0; iFace < nBoundFace; iFace ++)
    {
        string &bcName = bcNameList[iFace];
        streamsize nameSize = bcName.size();
        for (unsigned int iChar = 0; iChar < nameSize; ++ iChar)
        {
            bcNameChar[count ++] = bcName[iChar];
        }
        bcNameChar[count ++] = '\0';
    }

    CreateAndWriteData(grpData, "BCName", 1, totalSize, PHSTRING, bcNameChar);

    delete [] bcNameList;
    delete [] bcNameChar;
    H5Gclose(grpData);
}

void IO_HDF5Write::WriteBCRegion(hid_t loc_id, Grid *gridIn)
{
    hid_t grpData;
    int iMin, iMax, jMin, jMax, kMin, kMax, bcType;
    int sLr, sNd;

    grpData = GetGroup(loc_id, "BCRegion");

    StructGrid *grid = StructGridCast(gridIn);
    StructBCSet *compositeBCRegion = grid->GetCompositeBCRegion();

    int nBCRegion = compositeBCRegion->GetnBCRegion();
    CreateAndWriteData(grpData, "nBCRegion", 1, 1, PHINT, &nBCRegion);

    int nBCTotal = nBCRegion;
    for (int ibcregion = 0; ibcregion < nBCRegion; ++ ibcregion)
    {
        StructBC *bcregion = compositeBCRegion->GetBCRegion(ibcregion);
        bcType = bcregion->GetBCType();
        if (bcType == -1)
        {
            nBCTotal ++;
        }
    }

    int **BCInfor = NewPointer2<int>(nBCTotal, 9);
    int iBCinfo = 0;
    for (int ibcregion = 0; ibcregion < nBCRegion; ++ ibcregion)
    {
        int index = 0;
        StructBC *bcregion = compositeBCRegion->GetBCRegion(ibcregion);
        bcregion->GetIJKRegion(iMin, iMax, jMin, jMax, kMin, kMax);
        bcType = bcregion->GetBCType();
        sLr = bcregion->GetFaceLeftOrRightIndex();
        sNd = bcregion->GetFaceDirection();

        BCInfor[iBCinfo][index++] = iMin;
        BCInfor[iBCinfo][index++] = iMax;
        BCInfor[iBCinfo][index++] = jMin;
        BCInfor[iBCinfo][index++] = jMax;
        BCInfor[iBCinfo][index++] = kMin;
        BCInfor[iBCinfo][index++] = kMax;
        BCInfor[iBCinfo][index++] = bcType;
        BCInfor[iBCinfo][index++] = sLr;
        BCInfor[iBCinfo][index++] = sNd;
        iBCinfo ++;

        if (bcType == -1)
        {
            index = 0;
            int *t_st = bcregion->GetTargetStart();
            int *t_ed = bcregion->GetTargetEnd();
            int neighborZoneIndex = bcregion->GetTargetRegionBlock();

            BCInfor[iBCinfo][index++] = t_st[0];
            BCInfor[iBCinfo][index++] = t_ed[0];
            BCInfor[iBCinfo][index++] = t_st[1];
            BCInfor[iBCinfo][index++] = t_ed[1];
            BCInfor[iBCinfo][index++] = t_st[2];
            BCInfor[iBCinfo][index++] = t_ed[2];
            BCInfor[iBCinfo][index++] = neighborZoneIndex;
            BCInfor[iBCinfo][index++] = 0;
            BCInfor[iBCinfo][index++] = 0;
            iBCinfo ++;
        }
    }

    CreateAndWriteData(grpData, "BCInfor", nBCTotal, 9, PHINT, BCInfor[0]);
    DelPointer2(BCInfor);

    H5Gclose(grpData);
}

void IO_HDF5Write::WriteStrBCName(hid_t loc_id, Grid *gridIn)
{
     hid_t grpData;

    grpData = GetGroup(loc_id, "BCName");

    string boundaryName[100] = {"NO_BOUNDARY_CONDITION", "EXTRAPOLATION", "SOLID_SURFACE", "SYMMETRY", 
                            "FARFIELD", "INFLOW", "OUTFLOW", "POLE", "Wall_8", "Wall_9", "Wall_10", "Wall_11", 
                            "Wall_12", "Wall_13", "Wall_14", "Wall_15", "Wall_16", "Wall_17", "Wall_18", "Wall_19",
                            "Wall_20", "Wall_21", "Wall_22", "Wall_23", "Wall_24", "Wall_25", "Wall_26", "Wall_27", 
                            "Wall_28", "Wall_29", "Wall_30", "Wall_31", "Wall_32", "Wall_33", "Wall_34", "Wall_35",
                            "Wall_36", "Wall_37", "Wall_38", "Wall_39", "Wall_40", "Wall_41", "Wall_42", "Wall_43",
                            "Wall_44", "Wall_45", "Wall_46", "Wall_47", "Wall_48", "Wall_49", "Wall_50"};
    boundaryName[52] = "PRESSURE_INLET";
    boundaryName[62] = "PRESSURE_OUTLET";
    boundaryName[71] = "POLE1";
    boundaryName[72] = "POLE2";
    boundaryName[73] = "POLE3";

    StructGrid *grid = StructGridCast(gridIn);
    StructBCSet *composite_bcregion = grid->GetCompositeBCRegion();

    int nBCRegion = static_cast<int> (composite_bcregion->GetnBCRegion());
    string *bcNameList = new string[nBCRegion];
    uint_long totalSize = 0;
    for (int ibcregion = 0; ibcregion < nBCRegion; ibcregion ++)
    {
        StructBC *bcRegion = composite_bcregion->GetBCRegion(ibcregion);

        string bcName = bcRegion->GetBCName();
        int bctype = bcRegion->GetBCType();
        if (bcName == "")
        {
            if (bctype == -1)
            {
                bcName = "Connection";
            }
            else
            {
                bcName = boundaryName[bctype];
                if ((bctype >7) && (bctype <51))
                {
                    bctype = 2;
                    bcRegion->SetBCType(bctype);
                }
            }
            bcRegion->SetBCName(bcName);
        }
        bcNameList[ibcregion] = bcName;

        totalSize += static_cast< uint_long > (bcNameList[ibcregion].size());
        totalSize += 1;
    }

    char *bcNameChar = new char[totalSize];
    unsigned int count = 0;
    for (int boco = 0; boco < nBCRegion; ++ boco)
    {
        string &bcName = bcNameList[boco];
        streamsize nameSize = bcName.size();
        for (unsigned int iChar = 0; iChar < nameSize; ++ iChar)
        {
            bcNameChar[count ++] = bcName[iChar];
        }
        bcNameChar[count ++] = '\0';
    }

    CreateAndWriteData(grpData, "BCName", 1, totalSize, PHSTRING, bcNameChar);

    delete [] bcNameChar;
    delete [] bcNameList;
    H5Gclose(grpData);
}

void IO_HDF5Write::WriteInterfaceInfo(hid_t loc_id, Grid *gridIn)
{
    hid_t grpData;

    grpData = GetGroup(loc_id, "InterfaceInfo");

    InterfaceInfo *interfaceInfo = gridIn->GetInterfaceInfo();

    int nIFace = 0;
    if (interfaceInfo)
    {
        nIFace = interfaceInfo->GetNIFace();
    }
    CreateAndWriteData(grpData, "nIFace", 1, 1, PHINT, &nIFace);

    if (nIFace)
    {
        int *interFace2ZoneID       = interfaceInfo->GetInterFace2ZoneID();
        int *interFace2InterFaceID  = interfaceInfo->GetInterFace2InterFaceID();
        int *interFace2BoundaryFace = interfaceInfo->GetInterFace2BoundaryFace();

        CreateAndWriteData(grpData, "interFace2ZoneID"      , 1, nIFace, PHINT, interFace2ZoneID);
        CreateAndWriteData(grpData, "interFace2InterFaceID" , 1, nIFace, PHINT, interFace2InterFaceID);
        CreateAndWriteData(grpData, "interFace2BoundaryFace", 1, nIFace, PHINT, interFace2BoundaryFace);
    }

    H5Gclose(grpData);
}

void IO_HDF5Write::WriteLnkInfo(hid_t loc_id, Grid *gridIn, int gridID)
{
    StructGrid *grid = StructGridCast(gridIn);
    StructBCSet *composite_bcregion = grid->GetCompositeBCRegion();

    int nBCLnkTotal = 0;
    int nBCRegion = composite_bcregion->GetnBCRegion();
    for (int ibcregion = 0; ibcregion < nBCRegion; ibcregion ++)
    {
        StructBC *structBC = composite_bcregion->GetBCRegion(ibcregion);

        int bctype = structBC->GetBCType();
        string bcName = structBC->GetBCName();
        if (bctype == PHENGLEI::INTERFACE && !(bcName == "Periodic_up" || bcName == "Periodic_down"))
        {
            nBCLnkTotal ++;
        }
    }

    if (!nBCLnkTotal)
    {
        return;
    }

    hid_t grpData;
    grpData = GetGroup(loc_id, "LnkInfo");

    int **LnkInfor = NewPointer2<int>(nBCLnkTotal * 2, 7);
    int iLnkinfo = 0;
    for (int ibcregion = 0; ibcregion < nBCRegion; ibcregion ++)
    {
        StructBC *structBC = composite_bcregion->GetBCRegion(ibcregion);

        int bctype = structBC->GetBCType();
        string bcName = structBC->GetBCName();
        if (bctype >= 0 || bcName == "Periodic_up" || bcName == "Periodic_down") continue;

        int index = 0;
        int icount = 0;
        int nbt = structBC->GetTargetRegionBlock() + 1;
        int *lnkInfo = structBC->GetInkInfo();

        LnkInfor[iLnkinfo][index++] = lnkInfo[icount++];
        LnkInfor[iLnkinfo][index++] = lnkInfo[icount++];
        LnkInfor[iLnkinfo][index++] = lnkInfo[icount++];
        LnkInfor[iLnkinfo][index++] = lnkInfo[icount++];
        LnkInfor[iLnkinfo][index++] = lnkInfo[icount++];
        LnkInfor[iLnkinfo][index++] = lnkInfo[icount++];
        LnkInfor[iLnkinfo][index++] = bctype;

        iLnkinfo ++;
        index = 0;

        LnkInfor[iLnkinfo][index++] = lnkInfo[icount++];
        LnkInfor[iLnkinfo][index++] = lnkInfo[icount++];
        LnkInfor[iLnkinfo][index++] = lnkInfo[icount++];
        LnkInfor[iLnkinfo][index++] = lnkInfo[icount++];
        LnkInfor[iLnkinfo][index++] = lnkInfo[icount++];
        LnkInfor[iLnkinfo][index++] = lnkInfo[icount++];
        LnkInfor[iLnkinfo][index++] = nbt;

        iLnkinfo ++;
    }

    CreateAndWriteData(grpData, "LnkInfor", nBCLnkTotal * 2, 7, PHINT, LnkInfor[0]);
    DelPointer2(LnkInfor);

    H5Gclose(grpData);
}

IO_HDF5Read::IO_HDF5Read(PHVectorString1D &gridNameListIn, Region *region_in)
{
    this->gridNameList = gridNameListIn;
    this->region = region_in;

    this->numberOfGridGroups = gridNameList.size();
    this->grids.resize(0);
    this->myProcID = PHMPI::GetCurrentProcessorID();
    this->serverProcID = 0;
    this->zoneStart = 0;
    this->presentGridGroup = 0;
    this->fileID = -1;
    this->isFileProcessor = false;
    this->taskType = GetTaskCode();
    this->isP2PMode = false;

    blockProcDump = 0;
    nZones        = 0;
    blockProc     = 0;
    blockType     = 0;
    blockIdx      = 0;
    blockProcGrid = 0;
    blockFileProc = 0;
    nZonesList.resize(numberOfGridGroups);
    localGridNameList.resize(numberOfGridGroups);
    version = 0.0;
    globalFileNumber = 0;
}

IO_HDF5Read::~IO_HDF5Read()
{
    delete [] blockProcDump;
    delete [] blockProc;
    delete [] blockType;
    delete [] blockIdx;
    delete [] blockProcGrid;
    delete [] blockFileProc;
}

void IO_HDF5Read::CheckNumberofFiles(int iGrid)
{
    string gridFileName = gridNameList[iGrid];

    string gridfileProc = gridFileName;
    gridfileProc = AddSymbolToFileName(gridfileProc, "_", myProcID);
    bool fileExist = FileExist(gridfileProc);

    int nFileProcessors = 0;
    if (!fileExist)
    {
        if (myProcID == serverProcID)
        {
            // The server processor must be file processor.
            TK_Exit::FileOpenErrorExit(gridfileProc);
        }
    }
    else
    {
        nFileProcessors = 1;
    }

    PH_AllReduce(&nFileProcessors, &globalFileNumber, 1, MPI_SUM);

    PHMPI::SetNumberofFiles(globalFileNumber);
    WriteLogFile("Global grid files number: ", globalFileNumber);

    nZonesList[iGrid].resize(globalFileNumber);

    localGridNameList[iGrid].resize(0);
    for (int iFile = 0; iFile < globalFileNumber; ++ iFile)
    {
        string gridfile = gridFileName;
        gridfile = AddSymbolToFileName(gridfile, "_", iFile);
        localGridNameList[iGrid].push_back(gridfile);
    }

    if (globalFileNumber > 1)
    {
        isP2PMode = true;
    }
}

void IO_HDF5Read::InitGridFileNameListAndFileProcessor()
{
    for (int iGrid = 0; iGrid < numberOfGridGroups; ++ iGrid)
    {
        CheckNumberofFiles(iGrid);
    }

    if (PHMPI::GetNumberOfProcessor() % PHMPI::GetNumberofFiles() != 0)
    {
        ostringstream oss;
        oss << "   ---- numberOfProcessor = " << PHMPI::GetNumberOfProcessor() << ", numberofFiles = " << PHMPI::GetNumberofFiles() << " ----\n" << "Abnormally Exit Program!!!\n";
        TK_Exit::ExceptionExit(oss.str());
    }

    int fileProcessorInterval = PHMPI::GetIntervalofProcessorFiles();
    isFileProcessor = (myProcID % fileProcessorInterval == 0);
    //isFileProcessor = (myProcID % fileProcessorInterval == 0 && myProcID < globalFileNumber * fileProcessorInterval);

    this->fileID = PHMPI::GetFileIndexofCurrentProcessor();
}

hid_t IO_HDF5Read::OpenH5File(const string &filename)
{
    hid_t file;

    file = OpenHDF5File(AddSymbolToFileName(filename, "_", this->fileID));
    return file;
}

bool IO_HDF5Read::Run()
{
    bool HDF5FileRead = false;

    InitGridFileNameListAndFileProcessor();
    ReadNumberOfZones();
    ReadTotalInfo();
    SetMixGridInfo();
    RedistributeZonesIfNeed();
    SetGlobalZoneLayout();
    ClassifyGridSystem();
    ReadEachGrid();
    CreateZones();
    ReadWallDist();

    HDF5FileRead = true;
    return HDF5FileRead;
}

void IO_HDF5Read::ReadWallDist()
{
    if (gridNameList.size() != 1)
    {
        return;
    }

    string WallDistFileName = ChangeExtensionOfFileName(gridNameList.front(), "wdt");

    bool fileExist = true;
    if (isFileProcessor)
    {
        fileExist = FileExist(AddSymbolToFileName(WallDistFileName, "_", this->fileID));
    }

    bool fileCheck = true;
    PH_AllReduce(&fileExist, &fileCheck, 1, MPI_MIN);
    if (!fileCheck)
    {
        return;
    }

    hid_t file;
    int zoneCount = 0;

    for (int iGrid = 0; iGrid < numberOfGridGroups; ++ iGrid)
    {
        for (int iFile = 0; iFile < localGridNameList[iGrid].size(); ++ iFile)
        {
            file = OpenHDF5File(ChangeExtensionOfFileName(localGridNameList[iGrid][iFile], "wdt"));
            if (file < 0)
            {
                return;
            }

            for (int iZone = 0; iZone < nZonesList[iGrid][iFile]; ++ iZone)
            {
                int gridIndex = blockIdx[zoneCount ++];
                if (myProcID != blockProc[gridIndex])
                {
                    continue;
                }

                ostringstream oss;
                oss << "Grid-" << gridIndex;
                string dataName = oss.str();

                Grid *grid = GetGrid(gridIndex, 0);
                grid->InitWallDist();

                int nTotalCell = grid->GetNTotalCell();
                RDouble *walldist = new RDouble[nTotalCell];
                ReadData(file, walldist, dataName);

                grid->ReadWalldist(walldist);
                delete [] walldist;
            }

            H5Fclose(file);
        }
    }
}

void IO_HDF5Read::ReadVersionInfo(hid_t loc_id)
{
    hid_t grpData;

    grpData = OpenGroup(loc_id, "Version");

    ReadData(grpData, &this->version, "VersionData");

    H5Gclose(grpData);
}

void IO_HDF5Read::ReadNumberOfZones()
{
    hid_t file;

    for (int iGrid = 0; iGrid < numberOfGridGroups; ++ iGrid)
    {
        for (int iFile = 0; iFile < globalFileNumber; ++ iFile)
        {
            file = OpenHDF5File(localGridNameList[iGrid][iFile]);
            if (file < 0)
            {
                TK_Exit::ExceptionExit("Grid file can not open !");
            }

            hid_t grpData;
            grpData = OpenGroup(file, "Information");

            int numberOfZones = 0;
            ReadData(grpData, &numberOfZones, "nBlocks");
            nZones += numberOfZones;
            nZonesList[iGrid][iFile] = numberOfZones;

            H5Gclose(grpData);
            H5Fclose(file);
        }
    }

    PHMPI::SetNumberofGlobalZones   (nZones);
    PHMPI::CreateZoneProcessorIDDump(nZones);
    PHMPI::CreateZoneProcessorID    (nZones);
    PHMPI::CreateZoneGridID         (nZones);
    PHMPI::CreateZoneGridType       (nZones);
    PHMPI::CreateZoneFileID         (nZones);

    int numberOfProcessor = PHMPI::GetNumberOfProcessor();
    if (nZones < numberOfProcessor)
    {
        ostringstream oss;
        oss << "   ---- numberOfProcessor = " << numberOfProcessor << ", numberOfTotalZones = " << nZones << " ----\n" << "Abnormally Exit Program!!!\n";
        TK_Exit::ExceptionExit(oss.str());
    }
}

void IO_HDF5Read::ReadNumberOfZonesOld(hid_t loc_id)
{
    hid_t grpData;
    grpData = OpenGroup(loc_id, "Information");

    int *nZonesListTemp = new int[globalFileNumber];
    SetField(nZonesListTemp, 0, globalFileNumber);

    int numberOfZones = 0;
    if (isFileProcessor)
    {
        ReadData(grpData, &numberOfZones, "nBlocks");
        nZonesListTemp[fileID] = numberOfZones;
        //MPI_Allgather(&numberOfZones, 1, MPI_INT, nZonesList, 1, MPI_INT, MPI_COMM_WORLD);
    }

    PH_AllReduce(&numberOfZones, &nZones, 1, MPI_SUM);
    //PH_AllReduce(nZonesListTemp, nZonesList, globalFileNumber, MPI_SUM);
    delete [] nZonesListTemp;

    PHMPI::SetNumberofGlobalZones   (nZones);
    PHMPI::CreateZoneProcessorIDDump(nZones);
    PHMPI::CreateZoneProcessorID    (nZones);
    PHMPI::CreateZoneGridID         (nZones);
    PHMPI::CreateZoneGridType       (nZones);
    PHMPI::CreateZoneFileID         (nZones);

    int numberOfProcessor = PHMPI::GetNumberOfProcessor();
    if (nZones < numberOfProcessor)
    {
        ostringstream oss;
        oss << "   ---- numberOfProcessor = " << numberOfProcessor << ", numberOfTotalZones = " << nZones << " ----\n" << "Abnormally Exit Program!!!\n";
        TK_Exit::ExceptionExit(oss.str());
    }

    H5Gclose(grpData);
}

void IO_HDF5Read::SetMixGridInfo()
{
    int nsimutask = GlobalDataBase::GetIntParaFromDB("nsimutask");
    if (!(nsimutask == CREATE_GRID || nsimutask == GRIDCHECK_SOLVER))
    {
        return;
    }

    if (nsimutask == CREATE_GRID)
    {
        int gridObj = GlobalDataBase::GetIntParaFromDB("gridobj");
        if (gridObj != PHSPACE::GRID_CONVERSION && gridObj != PHSPACE::GRID_STRTOUNS)
        {
            return;
        }

        int gridType = GlobalDataBase::GetIntParaFromDB("gridtype");
        if (gridType != PHSPACE::MIXGRID && gridObj != PHSPACE::GRID_STRTOUNS)
        {
            return;
        }
    }

    CreatetZoneProcessorID_INP(nZones);
    int *global_block_proc_inp = GetZoneProcessorID_INP();
    for (int izone = 0; izone < nZones; ++ izone)
    {
        global_block_proc_inp[izone] = blockProc[izone];
    }
}

void IO_HDF5Read::ReadTotalInfo()
{
    blockProc     = new int[nZones];
    blockType     = new int[nZones];
    blockIdx      = new int[nZones];
    blockFileProc = new int[nZones];

    SetField(blockProc    , 0, nZones);
    SetField(blockType    , 0, nZones);
    SetField(blockIdx     , 0, nZones);
    SetField(blockFileProc, 0, nZones);

    int processorNumInEachFile = PHMPI::GetNumberOfProcessor() / PHMPI::GetNumberofFiles();

    hid_t file;
    int count = 0;

    for (int iGrid = 0; iGrid < numberOfGridGroups; ++ iGrid)
    {
        for (int iFile = 0; iFile < globalFileNumber; ++ iFile)
        {
            file = OpenHDF5File(localGridNameList[iGrid][iFile]);

            hid_t grpData;
            grpData = OpenGroup(file, "Information");

            int *blockProc_local = ReadIntOneRow(grpData, "block_proc");
            int *blockIdx_local  = ReadIntOneRow(grpData, "block_idx");
            int *blockType_local = ReadIntOneRow(grpData, "block_type");

            for (int izone = 0; izone < nZonesList[iGrid][iFile]; ++ izone)
            {
                blockProc[count + izone] = blockProc_local[izone];
                blockType[count + izone] = blockType_local[izone];
                blockIdx [count + izone] = blockIdx_local [izone];

                blockFileProc[count + izone] = processorNumInEachFile * iFile;
            }
            delete [] blockProc_local;
            delete [] blockIdx_local ;
            delete [] blockType_local;

            count += nZonesList[iGrid][iFile];
            H5Gclose(grpData);
            H5Fclose(file);
        }
    }

    blockProcDump = new int[nZones];
    for (int izone = 0; izone < nZones; ++ izone)
    {
        blockProcDump[izone] = blockProc[izone];
    }
}

void IO_HDF5Read::ReadTotalInfoOld(hid_t loc_id)
{
    hid_t grpData;
    grpData = OpenGroup(loc_id, "Information");

    int *blockProcTemp     = new int[nZones];
    int *blockTypeTemp     = new int[nZones];
    int *blockIdxTemp      = new int[nZones];
    int *blockFileProcTemp = new int[nZones];

    blockProc     = new int[nZones];
    blockType     = new int[nZones];
    blockIdx      = new int[nZones];
    blockFileProc = new int[nZones];

    SetField(blockProcTemp    , 0, nZones);
    SetField(blockTypeTemp    , 0, nZones);
    SetField(blockIdxTemp     , 0, nZones);
    SetField(blockFileProcTemp, 0, nZones);

    if (isFileProcessor)
    {
        int *blockProc_local = ReadIntOneRow(grpData, "block_proc");
        int *blockIdx_local  = ReadIntOneRow(grpData, "block_idx");
        int *blockType_local = ReadIntOneRow(grpData, "block_type");

        int count = 0;
        for (int ifile = 0; ifile < fileID; ++ ifile)
        {
            count += nZonesList[0][ifile];
        }

        int myid = PHMPI::GetCurrentProcessorID();
        for (int izone = 0; izone < nZonesList[0][fileID]; ++ izone)
        {
            blockProcTemp[count + izone] = blockProc_local[izone];
            blockTypeTemp[count + izone] = blockType_local[izone];
            blockIdxTemp [count + izone] = blockIdx_local [izone];

            blockFileProcTemp[count + izone] = myid;
        }

        delete [] blockProc_local;
        delete [] blockIdx_local ;
        delete [] blockType_local;
    }

    PH_AllReduce(blockProcTemp    , blockProc    , nZones, MPI_SUM);
    PH_AllReduce(blockIdxTemp     , blockIdx     , nZones, MPI_SUM);
    PH_AllReduce(blockTypeTemp    , blockType    , nZones, MPI_SUM);
    PH_AllReduce(blockFileProcTemp, blockFileProc, nZones, MPI_SUM);

    delete [] blockProcTemp;
    delete [] blockIdxTemp ;
    delete [] blockTypeTemp;
    delete [] blockFileProcTemp;

    blockProcDump = new int[nZones];
    for (int izone = 0; izone < nZones; ++ izone)
    {
        blockProcDump[izone] = blockProc[izone];
    }

    H5Gclose(grpData);
}

void IO_HDF5Read::ReadEachGrid()
{
    int *blockProc = GetBlockProc();
    int *blockType = GetBlockType();

    hid_t file;
    zoneStart = 0;
    for (int iGrid = 0; iGrid < numberOfGridGroups; ++ iGrid)
    {
        presentGridGroup = iGrid;
        int nZonesOfGridGroup = 0;
        for (int iFile = 0; iFile < globalFileNumber; ++ iFile)
        {
            file = OpenHDF5File(localGridNameList[iGrid][iFile]);

            hid_t grpData;
            grpData = OpenGroup(file, "Information");
            int *blockIdx_local = ReadIntOneRow(grpData, "block_idx");
            H5Gclose(grpData);

            for (int iZone = 0; iZone < nZonesList[iGrid][iFile]; ++ iZone)
            {
                int gridIndex = blockIdx_local[iZone];
                if (myProcID != blockProc[gridIndex + zoneStart])
                {
                    AddGrid(0);
                    continue;
                }

                int gridType = blockType[gridIndex + zoneStart];
                if (gridType == PHSPACE::UNSTRUCTGRID)
                {
                    ReadUnstructGrid(file, gridIndex);
                }
                else
                {
                    ReadStructGrid(file, gridIndex);
                }
            }

            nZonesOfGridGroup += nZonesList[iGrid][iFile];
            delete [] blockIdx_local;
            H5Fclose(file);
        }
        zoneStart += nZonesOfGridGroup;
    }
}

void IO_HDF5Read::ReadEachGridOld(hid_t loc_id)
{
    int *blockProc = GetBlockProc();
    int *blockType = GetBlockType();

    int nZones = GetNZones();
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        if (myProcID != blockProc[iZone])
        {
            AddGrid(0);
            continue;
        }

        int gridType = blockType[iZone];
        if (gridType == PHSPACE::UNSTRUCTGRID)
        {
            ReadUnstructGrid(loc_id, iZone);
        }
        else
        {
            ReadStructGrid(loc_id, iZone);
        }
    }
}

void IO_HDF5Read::ReadGridCoordinates(hid_t loc_id, Grid *gridIn)
{
    PrintToWindow("    Reading coordinates ...\n");
    hid_t grpData;

    grpData = OpenGroup(loc_id, "GridCoordinates");

    RDouble *CoordinateX = ReadDoubleOneRow(grpData, "CoordinateX");
    RDouble *CoordinateY = ReadDoubleOneRow(grpData, "CoordinateY");
    RDouble *CoordinateZ = ReadDoubleOneRow(grpData, "CoordinateZ");

    gridIn->SetX(CoordinateX);
    gridIn->SetY(CoordinateY);
    gridIn->SetZ(CoordinateZ);

    int axisup = 1;
    PHSPACE::GlobalDataBase::GetData("axisup", &axisup, PHINT, 1);
    if (axisup == 2)
    {
        gridIn->TurnZAxis2YAxis();
    }

    gridIn->ComputeMinMaxBox();

    H5Gclose(grpData);
}

void IO_HDF5Read::ReadInterfaceInfo(hid_t loc_id, Grid *gridIn)
{
    PrintToWindow("    Reading interface ...\n");
    hid_t grpData;

    grpData = OpenGroup(loc_id, "InterfaceInfo");

    InterfaceInfo *interfaceInfo = 0;
    gridIn->SetInterfaceInfo(interfaceInfo);

    int nIFace = 0;
    ReadData(grpData, &nIFace, "nIFace");

    if (nIFace > 0)
    {
        interfaceInfo = new InterfaceInfo(nIFace, gridIn);
        gridIn->SetInterfaceInfo(interfaceInfo);

        int *interFace2ZoneID       = interfaceInfo->GetInterFace2ZoneID();
        int *interFace2InterFaceID  = interfaceInfo->GetInterFace2InterFaceID();
        int *interFace2BoundaryFace = interfaceInfo->GetInterFace2BoundaryFace();

        int *interFace2ZoneIDTemp       = ReadIntOneRow(grpData, "interFace2ZoneID");
        int *interFace2InterFaceIDTemp  = ReadIntOneRow(grpData, "interFace2InterFaceID");
        int *interFace2BoundaryFaceTemp = ReadIntOneRow(grpData, "interFace2BoundaryFace");

        for (int iface = 0; iface < nIFace; ++ iface)
        {
            interFace2ZoneID[iface]       = interFace2ZoneIDTemp[iface];
            interFace2InterFaceID[iface]  = interFace2InterFaceIDTemp[iface];
            interFace2BoundaryFace[iface] = interFace2BoundaryFaceTemp[iface];
        }

        delete [] interFace2ZoneIDTemp;
        delete [] interFace2InterFaceIDTemp;
        delete [] interFace2BoundaryFaceTemp;
    }

    H5Gclose(grpData);
}

void IO_HDF5Read::ReadStructGrid(hid_t loc_id, int iZone)
{
    PrintToWindow("Reading Structured Grid Of Zone ", iZone, "...\n");

    hid_t  grpData;
    string grpName;

    GridID *index = new GridID(iZone + zoneStart);
    Grid *grid = new StructGrid();

    grid->InitGrid(index, 0, PHSPACE::GetDim(), STRUCTGRID);
    grid->SetIBlock(presentGridGroup);

    ostringstream oss;
    oss << "Grid-" << iZone;
    grpName = oss.str();
    grpData = OpenGroup(loc_id, grpName);

    StructGrid *strGrid = StructGridCast(grid);
    int *iDimensions = ReadIntOneRow(grpData, "iDimensions");
    strGrid->SetNI(iDimensions[0]);
    strGrid->SetNJ(iDimensions[1]);
    strGrid->SetNK(iDimensions[2]);
    strGrid->SetBasicDimension();

    ReadGridCoordinates(grpData, grid);
    strGrid->SetArrayLayout();

    ReadStrBCName(grpData, grid);
    ReadBCRegion(grpData, grid);
    ReadInterfaceInfo(grpData, grid);

    grid->ReviseNeighborZoneIndex(zoneStart);
    AddGrid(grid);
    delete [] iDimensions;
    H5Gclose(grpData);
}

void IO_HDF5Read::ReadLnk(hid_t loc_id, int iZone)
{
    int numberOfOriginalBlocks = PHMPI::GetNumberofGlobalZones();
    Grid **OrdinaryGrid = new Grid *[numberOfOriginalBlocks];

    hid_t  grpData;
    string grpName;

    ostringstream oss;
    oss << "Grid-" << iZone;
    grpName = oss.str();
    grpData = OpenGroup(loc_id, grpName);

    OrdinaryGrid[iZone] = GetGrid(iZone, 0);

    ReadLnkInfo(grpData, OrdinaryGrid[iZone]);
    H5Gclose(grpData);
}

void IO_HDF5Read::ReadBCRegion(hid_t loc_id, Grid *gridIn)
{
    hid_t grpData;

    grpData = OpenGroup(loc_id, "BCRegion");

    StructGrid *grid = StructGridCast(gridIn);
    StructBCSet *compositeBCRegion = grid->GetCompositeBCRegion();
    int nBCRegion = compositeBCRegion->GetnBCRegion();

    int **BCInfor = ReadIntTwoRow(grpData, "BCInfor");

    int iMin, iMax, jMin, jMax, kMin, kMax, bcType;
    int sLr, sNd;
    int iBCinfo = 0;
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        int index = 0;
        StructBC *bcregion = compositeBCRegion->GetBCRegion(iBCRegion);

        iMin   = BCInfor[iBCinfo][index++];
        iMax   = BCInfor[iBCinfo][index++];
        jMin   = BCInfor[iBCinfo][index++];
        jMax   = BCInfor[iBCinfo][index++];
        kMin   = BCInfor[iBCinfo][index++];
        kMax   = BCInfor[iBCinfo][index++];
        bcType = BCInfor[iBCinfo][index++];
        sLr    = BCInfor[iBCinfo][index++];
        sNd    = BCInfor[iBCinfo][index++];
        iBCinfo ++;
        
        //! Restore periodic boundary bcType and correct read information only for PARTITION_GRID task.
        if (GetTaskCode() == PARTITION_GRID && PHSPACE::GlobalDataBase::GetIntParaFromDB("periodicType") != 0)
        {
            if (bcregion->GetBCName() == "Periodic_up"|| bcregion->GetBCName() == "Periodic_down")
            {
                bcType = PHENGLEI::PERIODIC;
                iBCinfo ++;
            }
        }

        bcregion->SetIJKRegion(iMin, iMax, jMin, jMax, kMin, kMax);
        bcregion->SetBCType(bcType);
        bcregion->SetFaceLeftOrRightIndex(sLr);
        bcregion->SetFaceDirection(sNd);
        bcregion->InitFaceDirectionIndex();

        if (bcType == -1)
        {
            index = 0;

            iMin = BCInfor[iBCinfo][index++];
            iMax = BCInfor[iBCinfo][index++];
            jMin = BCInfor[iBCinfo][index++];
            jMax = BCInfor[iBCinfo][index++];
            kMin = BCInfor[iBCinfo][index++];
            kMax = BCInfor[iBCinfo][index++];
            bcregion->SetTargetIJKRegion(iMin, iMax, jMin, jMax, kMin, kMax);

            int neighborZoneIndex = BCInfor[iBCinfo][index++];
            bcregion->SetTargetRegionBlock(neighborZoneIndex);

            iBCinfo ++;
        }
    }
    grid->SetBCFaceInfo();

    DelPointer2(BCInfor);
    H5Gclose(grpData);
}

void IO_HDF5Read::ReadLnkInfo(hid_t loc_id, Grid *gridIn)
{
    hid_t grpData;

    grpData = OpenGroup(loc_id, "LnkInfo");

    StructGrid *grid = StructGridCast(gridIn);
    StructBCSet *compositeBCRegion = grid->GetCompositeBCRegion();
    int nBCRegion = compositeBCRegion->GetnBCRegion();

    int **LnkInfor = ReadIntTwoRow(grpData, "LnkInfor");

    int iMin, iMax, jMin, jMax, kMin, kMax, nbt;
    int iLnkinfo = 0;
    int ndim = PHSPACE::GlobalDataBase::GetIntParaFromDB("ndim");

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        int index = 0;
        StructBC *bcregion = compositeBCRegion->GetBCRegion(iBCRegion);

        int bctype = bcregion->GetBCType();
        if (bctype >= 0) continue;

        iMin = LnkInfor[iLnkinfo][index++];
        iMax = LnkInfor[iLnkinfo][index++];
        jMin = LnkInfor[iLnkinfo][index++];
        jMax = LnkInfor[iLnkinfo][index++];

        if (ndim == THREE_D)
        {
            kMin = LnkInfor[iLnkinfo][index++];
            kMax = LnkInfor[iLnkinfo][index++];
        }
        else
        {
            kMin = 1;
            kMax = 1;
        }
        iLnkinfo++;

        bcregion->SetIJKRegion(iMin, iMax, jMin, jMax, kMin, kMax);

        index = 0;

        iMin = LnkInfor[iLnkinfo][index++];
        iMax = LnkInfor[iLnkinfo][index++];
        jMin = LnkInfor[iLnkinfo][index++];
        jMax = LnkInfor[iLnkinfo][index++];

        if (ndim == THREE_D)
        {
            kMin = LnkInfor[iLnkinfo][index++];
            kMax = LnkInfor[iLnkinfo][index++];
        }
        else
        {
            kMin = 1;
            kMax = 1;
        }

        nbt = LnkInfor[iLnkinfo][index++];
        bcregion->SetTargetIJKRegion(iMin, iMax, jMin, jMax, kMin, kMax);
        bcregion->SetTargetRegionBlock(nbt - 1);
        if (GetDim() == THREE_D)
        {
            bcregion->ComputeRelativeParameters();
        }

        iLnkinfo++;
    }

    DelPointer2(LnkInfor);
    H5Gclose(grpData);
}

void IO_HDF5Read::ReadStrBCName(hid_t loc_id, Grid *gridIn)
{
    hid_t grpData;
    
    grpData = OpenGroup(loc_id, "BCRegion");
    StructGrid *grid = StructGridCast(gridIn);

    int nBCRegion = 0;
    ReadData(grpData, &nBCRegion, "nBCRegion");

    grid->CreateCompositeBCRegion(nBCRegion);
    H5Gclose(grpData);

    grpData = OpenGroup(loc_id, "BCName");
    int gridIndex = grid->GetZoneID();
    StructBCSet *compositeBCRegion = grid->GetCompositeBCRegion();

    char *BCName = ReadCharData(grpData, "BCName");
    string *BCNameList = new string[nBCRegion];
    uint_long count = 0;

    for (int ibcregion = 0; ibcregion < nBCRegion; ++ ibcregion)
    {
        while (BCName[count] != '\0')
        {
            BCNameList[ibcregion].push_back(BCName[count]);
            ++ count;
        }
        ++ count;
    }

    for (int ibcregion = 0; ibcregion < nBCRegion; ++ ibcregion)
    {
        StructBC *bcregion = new StructBC(gridIndex, ibcregion);
        compositeBCRegion->SetBCRegion(ibcregion, bcregion);

        bcregion->SetBCName(BCNameList[ibcregion]);
    }

    delete [] BCNameList;
    delete [] BCName;
    H5Gclose(grpData);
}

void IO_HDF5Read::ReadUnstructGrid(hid_t loc_id, int iZone)
{
    PrintToWindow("Reading Unstructured Grid Of Zone ", iZone, "...\n");

    hid_t  grpData;
    string grpName;

    GridID *index = new GridID(iZone + zoneStart);
    Grid *grid = new UnstructGrid();

    grid->InitGrid(index, 0, PHSPACE::GetDim(), UNSTRUCTGRID);
    grid->SetIBlock(0);

    ostringstream oss;
    oss << "Grid-" << iZone;
    grpName = oss.str();
    grpData = OpenGroup(loc_id, grpName);

    int *iDimensions = ReadIntOneRow(grpData, "iDimensions");
    grid->SetNTotalNode(iDimensions[0]);
    grid->SetNTotalFace(iDimensions[1]);
    grid->SetNTotalCell(iDimensions[2]);

    PrintToWindow("  Grid Dimension  : ", grid->GetDim(), "\n");
    PrintToWindow("  Number of Points: ", iDimensions[0], "\n");
    PrintToWindow("  Number of Faces : ", iDimensions[1], "\n");
    PrintToWindow("  Number of Cells : ", iDimensions[2], "\n");

    ReadGridCoordinates(grpData, grid);
    ReadFaceTopology(grpData, grid);
    ReadInterfaceInfo(grpData, grid);
    ReadCellTopology(grpData, grid);
    ReadUnstrBCName(grpData, grid);
    ReadInterpointInfo(grpData, grid);

    grid->SetIBlock(presentGridGroup);
    grid->ReviseNeighborZoneIndex(zoneStart);

    AddGrid(grid);

    delete [] iDimensions;
    H5Gclose(grpData);
}

void IO_HDF5Read::ReadFaceTopology(hid_t loc_id, Grid *gridIn)
{
    PrintToWindow("    Reading face topology ...\n");
    hid_t grpData;

    grpData = OpenGroup(loc_id, "FaceTopology");

    UnstructGrid *grid = UnstructGridCast(gridIn);

    int *nodeNumberOfEachFace = ReadIntOneRow(grpData, "nodeNumberOfEachFace");
    grid->SetNodeNumberOfEachFace(nodeNumberOfEachFace);

    int *face2Node = ReadIntOneRow(grpData, "face2Node");
    grid->SetFace2Node(face2Node);

    int *leftCellOfFace  = ReadIntOneRow(grpData, "leftCellOfFace");
    int *rightCellOfFace = ReadIntOneRow(grpData, "rightCellOfFace");
    grid->SetLeftCellOfFace (leftCellOfFace);
    grid->SetRightCellOfFace(rightCellOfFace);

    int nTotalFace =grid->GetNTotalFace();
    int *nodePosi = new int[nTotalFace + 1];
    int nodeSum = 0;
    nodePosi[0] = 0;
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        nodeSum += nodeNumberOfEachFace[iFace];
        nodePosi[iFace+1] = nodeSum;
    }

    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        if (leftCellOfFace[iFace] < 0)
        {
            std::reverse(face2Node + nodePosi[iFace], face2Node + nodePosi[iFace+1]);
            SWAP(leftCellOfFace[iFace], rightCellOfFace[iFace]);
        }
    }
    delete [] nodePosi;

    int nBoundFace = 0;
    ReadData(grpData, &nBoundFace, "nBoundFace");
    grid->SetNBoundFace(nBoundFace);

    grid->SpecifyRightCellofBC();

    H5Gclose(grpData);
}

void IO_HDF5Read::ReadCellTopology(hid_t loc_id, Grid *gridIn)
{
    PrintToWindow("    Reading cell-to-node ...\n");
    hid_t grpData;

    grpData = OpenGroup(loc_id, "CellTopology");

    UnstructGrid *grid = UnstructGridCast(gridIn);

    int *cell2Node = ReadIntOneRow(grpData, "cell2Node");
    int *nodeNumberOfEachCell = ReadIntOneRow(grpData, "nodeNumberOfEachCell");

    grid->SetNodeNumberOfEachCell(nodeNumberOfEachCell);
    grid->SetCell2Node(cell2Node);

    H5Gclose(grpData);
}

void IO_HDF5Read::ReadUnstrBCName(hid_t loc_id, Grid *gridIn)
{
    PrintToWindow("    Reading boundary name ...\n");
    hid_t grpData;

    grpData = OpenGroup(loc_id, "BCName");

    UnstructGrid *grid = UnstructGridCast(gridIn);
    char *BCName = ReadCharData(grpData, "BCName");

    grpData = OpenGroup(loc_id, "FaceTopology");
    int *BcType = ReadIntOneRow(grpData, "BcType");

    int nBoundFace = grid->GetNBoundFace();
    string *BCNameList = new string[nBoundFace];
   
    uint_long count = 0;
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        while (BCName[count] != '\0')
        {
            BCNameList[iFace].push_back(BCName[count]);
            ++ count;
        }
        ++ count;
    }

    //! Create UnstructbcRegion.
    set<string> bcNameSet;
    for (int iFace = 0; iFace < nBoundFace; iFace++)
    {   
        bcNameSet.insert(BCNameList[iFace]);
    }

    int nBCRegionUnstruct = static_cast<int>(bcNameSet.size());
    
    grid->CreateUnstructBCRegion(nBCRegionUnstruct);

    UnstructBCSet *unstructBCRegion = grid->GetUnstructBCRegion();
    int* bcRegionIDofBCFace = new int[nBoundFace];
    set<string>::iterator iter;
    count = 0;
    for (iter = bcNameSet.begin(); iter != bcNameSet.end(); iter++)
    {
        UnstructBC *unstructBC = new UnstructBC(count);
        unstructBCRegion->SetBCRegion(count, unstructBC);
        unstructBC->SetBCName(*iter);

        for (int iFace = 0; iFace < nBoundFace; iFace++)
        {
            if (BCNameList[iFace] == unstructBC->GetBCName())
            {
                unstructBC->SetBCType(BcType[iFace]);
                bcRegionIDofBCFace[iFace] = count;
                vector<int> *faceIndex = unstructBC->GetFaceIndex();
                faceIndex->push_back(iFace);
            }
           
        }
        count++;
    }
    unstructBCRegion->SetBCFaceInfo(bcRegionIDofBCFace);
    delete [] BcType;
    delete [] BCNameList;
    delete [] BCName;
    H5Gclose(grpData);
}

void IO_HDF5Read::ReadInterpointInfo(hid_t loc_id, Grid *gridIn)
{
    if (isP2PMode)
    {
        return;
    }

    if (numberOfGridGroups != 1)
    {
        return;
    }

    InterpointInformation *ipointinfo = 0;
    gridIn->SetInterpointInfo(ipointinfo);

    hid_t grpData;
    grpData = OpenGroup(loc_id, "InterpointInfo");

    int nIPoint;
    ReadData(grpData, &nIPoint, "nIPoint");
    if (!nIPoint)
    {
        H5Gclose(grpData);
        return;
    }

    PrintToWindow("    Reading interpoint ...\n");

    ipointinfo = new InterpointInformation(nIPoint,gridIn);
    gridIn->SetInterpointInfo(ipointinfo);

    int *interPoint2ZoneID = ipointinfo->GetInterPoint2ZoneID();
    int *interPoint2InterPointID = ipointinfo->GetInterPoint2InterPointID();
    int *interPoint2GlobalPoint = ipointinfo->GetInterPoint2GlobalPoint();
    int *cellNumberOfInterPoint = ipointinfo->GetCellNumberOfInterPoint();
    int *totalZonesOfInterPoint = ipointinfo->GetTotalZonesOfInterPoint();
    int *labelOfInterPoint = ipointinfo->GetLabelOfInterPoint();

    ReadData(grpData, interPoint2ZoneID      , "interPoint2ZoneID");
    ReadData(grpData, interPoint2InterPointID, "interPoint2InterPointID");
    ReadData(grpData, interPoint2GlobalPoint , "interPoint2GlobalPoint");
    ReadData(grpData, cellNumberOfInterPoint , "cellNumberOfInterPoint");
    ReadData(grpData, totalZonesOfInterPoint , "totalZonesOfInterPoint");
    ReadData(grpData, labelOfInterPoint      , "labelOfInterPoint");

    H5Gclose(grpData);

    zoneConnectivityForPoint = new ZoneConnectivityForPoint();
}

bool IO_HDF5Read::IsZoneLayoutInitialized()
{
    if (this->blockIdx || this->blockProc || this->blockType)
    {
        return true;
    }
    else
    {
        return false;
    }
}

void IO_HDF5Read::RedistributeZonesIfNeed()
{
    int nZones = GetNZones();
    int numberOfProcessor = PHMPI::GetNumberOfProcessor();

    int *blockType     = GetBlockType();
    int *blockProc     = GetBlockProc();
    int *blockProcGrid = GetBlockProcGrid();
    int *blockProcDump = GetBlockProcDump();

    const int DETERMIN_BY_PARTITION = 0;
    const int REDISTRIBUTION        = 1;
    int zoneDistributionMethod = 0;
    int maxProc = blockProc[0];

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        maxProc = MAX(maxProc, blockProc[iZone]);
    }
    maxProc += 1;

    if (numberOfProcessor == maxProc)
    {
        zoneDistributionMethod = DETERMIN_BY_PARTITION;
    }
    else
    {
        zoneDistributionMethod = REDISTRIBUTION;
    }
    PHMPI::SetZoneDistributionMethod(zoneDistributionMethod);

    ostringstream oss;
    oss << "Max processor ID : " << maxProc << " of total zone " << nZones;
    if (zoneDistributionMethod == DETERMIN_BY_PARTITION)
    {
        oss << ", Zones' distribution is specified by partition ...";
    }
    else
    {
        oss << ", Zones' distribution is specified randomly ...";
    }
    oss << endl;
    WriteLogFile(oss);

    if (zoneDistributionMethod == REDISTRIBUTION)
    {
        if (nZones % numberOfProcessor == 0)
        {
            int nZonesPerProcessor = nZones / numberOfProcessor;
            int processorID = 0;
            for (int iZone = 0; iZone < nZones; ++ iZone)
            {
                blockProc[iZone] = processorID;
                if ((iZone + 1) % nZonesPerProcessor == 0)
                {
                    processorID ++;
                }
            }
        }
        else
        {
            PrintToWindow("Waining: nZones % numberOfProcessor != 0 !!");
            for (int izone = 0; izone < nZones; ++ izone)
            {
                blockProc[izone] = izone % numberOfProcessor;
            }
        }
    }
    else
    {
        if (numberOfProcessor == 1)
        {
            if (IsConvertGridToMixGrid())
            {
                if (blockType[0] == PHSPACE::STRUCTGRID)
                {
                    int maxproc = blockProc[0];
                    for (int izone = 0; izone < nZones; ++ izone)
                    {
                        if (blockProc[izone] > maxproc)
                        {
                            maxproc = blockProc[izone];
                        }
                    }
                    maxproc += 1;
                    PHMPI::SetNumberOfProcStructUsed(maxproc);
                }
                else
                {
                    for (int izone = 0; izone < nZones; ++ izone)
                    {
                        blockProcDump[izone] += PHMPI::GetNumberOfProcStructUsed();
                    }
                }
            }

            for (int izone = 0; izone < nZones; ++ izone)
            {
                blockProc[izone] = 0;
            }
        }
    }

    if (!blockProcGrid)
    {
        blockProcGrid = new int[nZones];
    }
    SetBlockProcGrid(blockProcGrid);

    int numberOfGridProcessor = PHMPI::GetNumberOfGridProcessor();
    if (numberOfGridProcessor != 0)
    {
        for (int iZone = 0; iZone < nZones; ++ iZone)
        {
            blockProcGrid[iZone] = iZone % numberOfGridProcessor + numberOfProcessor;
        }
    }
    else
    {
        for (int iZone = 0; iZone < nZones; ++ iZone)
        {
            blockProcGrid[iZone] = -1;
        }
    }
}

void IO_HDF5Read::SetGlobalZoneLayout()
{
    if (!PHMPI::GetZoneGridType())
    {
        PHMPI::SetNumberofGlobalZones(nZones);
        PHMPI::CreateZoneProcessorIDDump(nZones);
        PHMPI::CreateZoneProcessorID(nZones);
        PHMPI::CreateZoneGridID(nZones);
        PHMPI::CreateZoneGridType(nZones);
        PHMPI::CreateZoneFileID(nZones);
    }

    int *globalBlockProc     = PHMPI::GetZoneProcessorID();
    int *globalBlockProcGrid = PHMPI::GetZoneProcessorIDForGrid();
    int *globalBlockType     = PHMPI::GetZoneGridType();
    int *globalBlockIdx      = PHMPI::GetZoneGridID();
    int *globalBlockFile     = PHMPI::GetZoneFileID();
    int *globalBlockProcDump = PHMPI::GetZoneProcessorIDDump();

    set<int> zoneIDinCurrentProc;
    int myid = PHMPI::GetCurrentProcessorID();

    for (int izone = 0; izone < nZones; ++ izone)
    {
        globalBlockProc    [izone] = blockProc    [izone];
        globalBlockProcGrid[izone] = blockProcGrid[izone];
        globalBlockType    [izone] = blockType    [izone];
        globalBlockIdx     [izone] = blockIdx     [izone];
        globalBlockFile    [izone] = blockFileProc[izone];

        if (blockProc[izone] == myid || blockProcGrid[izone] == myid)
        {
            zoneIDinCurrentProc.insert(izone);
        }
    }

    ostringstream oss;
    oss << "Zone index in current processor : ";
    for (set<int>::iterator iter = zoneIDinCurrentProc.begin(); iter != zoneIDinCurrentProc.end(); ++ iter)
    {
        int partID = *iter;
        oss << partID << " ";
    }
    oss << endl;
    WriteLogFile(oss.str());

    for (int izone = 0; izone < nZones; ++ izone)
    {
        globalBlockProcDump[izone] = blockProcDump[izone];
    }
}

void IO_HDF5Read::ClassifyGridSystem()
{
    int nZones = PHMPI::GetNumberofGlobalZones();

    set<int> gridType;
    for (int izone = 0; izone < nZones; ++ izone)
    {
        gridType.insert(PHMPI::GetZoneGridType(izone));
    }

    int sysGridType = GetSystemGridType(gridType);
    GlobalDataBase::UpdateData("sys_gridtype", &sysGridType, PHINT, 1);

    vector<string> sysGridTypeName;
    sysGridTypeName.push_back("UNSTRUCTGRID");
    sysGridTypeName.push_back("STRUCTGRID");
    sysGridTypeName.push_back("MIXGRID");

    PrintToWindow("Grid Type: ", sysGridTypeName[sysGridType], "\n");
}

void IO_HDF5Read::AddGrid(Grid *grid)
{
    grids.push_back(grid);
}

void IO_HDF5Read::CreateZones()
{
    zoneStart = 0;
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int globalZoneIndex = zoneStart + iZone;

        int recvProc     = PHMPI::GetZoneProcessorID       (globalZoneIndex);
        int recvProcGrid = PHMPI::GetZoneProcessorIDForGrid(globalZoneIndex);

        if (myProcID == recvProc || myProcID == recvProcGrid)
        {
            int numberOfZoneInThisProcess = PHMPI::GetNumberofLocalZones();

            Grid *grid = grids[iZone];
            grid->GetGridID()->SetLocalIndex(numberOfZoneInThisProcess);

            Zone *zone = new Zone(globalZoneIndex);
            zone->AddGridAndCopyZonePara(grid);
            zone->AddGridToGlobal(globalZoneIndex);

            region->AddZone(zone);
        }
        else
        {
            region->AddZone(0);
            continue;
        }
    }
}

LIB_EXPORT void DumpHDF5Grid(const string &gridFileName, Grid **grids, int nBlocks)
{
    PrintToWindow("    Writing new HDF5 version grid to file ....\n");

    IO_HDF5Write writeGrid(gridFileName, grids, nBlocks);
    writeGrid.Run();
}

LIB_EXPORT bool ReadHDF5File(Region *region_in, PHVectorString1D &gridGroupNameList)
{
    bool HDF5FileRead = false;

    fstream file;
    string gridFileName = gridGroupNameList[0];
    OpenSerialFile(file, gridFileName, ios_base::in|ios_base::binary);
    int nZone;
    int sizeInt = sizeof(int);
    file.read((char *)&nZone, sizeInt);
    CloseFile(file);

    if (nZone <= static_cast <int> (1e8))
    {
        PrintToWindow("This grid file is old fts version!\n");
        WriteLogFile("This grid file is old fts version!\n");
        return HDF5FileRead;
    }
    else
    {
        PrintToWindow("This grid file is new HDF5 version!\n");
        WriteLogFile("This grid file is new HDF5 version!\n");
    }

    IO_HDF5Read readGrid(gridGroupNameList, region_in);
    HDF5FileRead = readGrid.Run();

    int isOldGrid = 0;
    GlobalDataBase::UpdateData("isOldGrid", &isOldGrid, PHINT, 1);

    return HDF5FileRead;
}

}