#include "Pre_UnstructGridPartition.h"
#include "IO_FileName.h"
#include "Glb_Dimension.h"
#ifdef USE_WINDOWS_X64
#else
#include "metis.h"
#endif
#include "Geo_UnstructBC.h"
#include "PHIO.h"
#include "Math_BasisFunction.h"
#include "Pre_GridConversion.h"
#include "Constants.h"
#include "IO_HDF5File.h"
#include "TK_Log.h"

using namespace std;

namespace PHSPACE
{
LIB_EXPORT Pre_UnstructGridPartition::Pre_UnstructGridPartition(UnstructGrid *gridSource, int nPartitionBlocks,
                                 int buildInterfaceMethod, string orgGridFileName, int parallelPartitionMethod)
{
    // if using ParMetis, the number of partitioned blocks is equal to the input.
    // else, it total nBlock = nProcessors * maxProcessor, where, maxProcessor in
    // this processor is the number of partitioned blocks.
    SetParallelPartitionMethod(parallelPartitionMethod);
    if (parallelPartitionMethod == METIS_METHOD || parallelPartitionMethod == METIS_OPENMP_METHOD)
    {
        this->maxProcessor = nPartitionBlocks;
        int nProcessors = PHMPI::GetNumberOfProcessor();
        if (nPartitionBlocks % nProcessors != 0)
        {
            ostringstream oss;
            oss << "Error: if using Metis to parallel partition, the number of zones wanted must INTER times "
                << "of the number of processors used !" << endl;
            PrintToWindow(oss);
        }
        this->nPartInEachProcIfMetis = nPartitionBlocks / nProcessors;
    }else
    {
        this->maxProcessor = nPartitionBlocks;
    }
    
    grids = 0;

    this->gridSource  = gridSource;
    this->npartmethod = buildInterfaceMethod;

    node2all = 0;
    face2all = 0;
    cell2all = 0;
    part = 0;
    cell_mark = 0;
    node_mark = 0;
    face_mark = 0;

    vtxdist = 0;
    interfaceIndexInBC_send    = 0;
    interfaceIndexInBC_recv    = 0;
    partIndexInNewPartion_recv = 0;
    part2GlobalZoneIndex       = 0;
    part2LocalZoneIndex        = 0;
    nIFaceInOrgBC              = new int[maxProcessor];

    this->orgGridFileName = orgGridFileName;

    parallelPartitionMethod = 1;

    partitionCellToNode = false;
    partitionFaceBC     = false;
    partitionWallDist   = false;
}

LIB_EXPORT Pre_UnstructGridPartition::~Pre_UnstructGridPartition()
{
    delete [] cell_mark;    cell_mark = NULL;
    delete [] face_mark;    face_mark = NULL;
    delete [] node_mark;    node_mark = NULL;
    delete [] part;    part = NULL;

    delete [] vtxdist;    vtxdist = NULL;

    delete [] interfaceIndexInBC_send;    interfaceIndexInBC_send = NULL;
    delete [] interfaceIndexInBC_recv;    interfaceIndexInBC_recv = NULL;
    delete [] partIndexInNewPartion_recv;    partIndexInNewPartion_recv = NULL;
    delete [] part2GlobalZoneIndex;    part2GlobalZoneIndex = NULL;
    delete [] part2LocalZoneIndex;    part2LocalZoneIndex = NULL;
    delete [] nIFaceInOrgBC;    nIFaceInOrgBC = NULL;

    delete [] pointNumber2zoneID;    pointNumber2zoneID = NULL;
    delete [] point2ZoneID;    point2ZoneID = NULL;
    delete [] pointPosition2zoneID;    pointPosition2zoneID = NULL;
    delete [] cellNumberForNodeValue;    cellNumberForNodeValue = NULL;

    FreeLocal2All();
}

void Pre_UnstructGridPartition::Init()
{
    CheckDoseNeedPartitionAdditionalData();
}

void Pre_UnstructGridPartition::CheckDoseNeedPartitionAdditionalData()
{
    if (this->orgGridFileName == "")
    {
        partitionCellToNode = false;
        partitionFaceBC     = false;
        return;
    }

    string originalCellToNodeFileName, originalFaceBCFileName, originalWallDistFileName;
    if (PHMPI::GetNumberOfProcessor() == 1)
    {
        originalCellToNodeFileName = ChangeExtensionOfFileName(orgGridFileName, "c2n");
        originalFaceBCFileName     = ChangeExtensionOfFileName(orgGridFileName, "bc");
        originalWallDistFileName   = ChangeExtensionOfFileName(orgGridFileName, "wdt");
    }
    else
    {
        if (parallelPartitionMethod == PARMETIS_METHOD)
        {
            partitionCellToNode = false;
            partitionFaceBC     = false;
            partitionWallDist   = false;
            return;
        }
        else
        {
            // Check if the first .c2n file exist or not.
            originalCellToNodeFileName = ChangeExtensionOfFileName(orgGridFileName, "c2n");
            originalFaceBCFileName     = ChangeExtensionOfFileName(orgGridFileName, "bc");
            originalWallDistFileName   = ChangeExtensionOfFileName(orgGridFileName, "wdt");
        }
    }
    
    //! Check if c2n file exist.
    partitionCellToNode = this->gridSource->ExistCell2Node();

    //! Check if face BC file exist.
    partitionFaceBC = SerialFileExist(originalFaceBCFileName);

    //! Check if walldist file exist.
    if (!this->gridSource->GetWallDist())
    {
        partitionWallDist = false;
    }
    else
    {
        partitionWallDist = true;
    }
}

Grid ** Pre_UnstructGridPartition::CreateAllPartGrid()
{
    this->CreatePartition();
    
    this->ComputeCellMark();

    this->ComputePoint2ZoneID();

    gridSource->computeBoundaryPointLabel();
    
    Grid **grids = new Grid * [maxProcessor];

    int nLocalZones = this->GetNumberofLocalZones();
    int interval = GetProgressInterval(nLocalZones, 10);
    int myid = PHMPI::GetCurrentProcessorID();
    int count = 0;

    ComputeCellNumberForNodeValue();

    for (int iproc = 0; iproc < maxProcessor; ++ iproc)
    {
        if (!DoseProcExist(iproc))
        {
            grids[iproc] = 0;
            continue;
        }

        ++ count;
        if (myid == PHMPI::GetServerProcessorID() && count > 0 && count % interval == 0)
        {
            ProgressMonitorToLogFile(count, nLocalZones, "grid partition reconstructing");
            ProgressMonitorToWindows(count, nLocalZones, "grid partition reconstructing");
        }

        GridID *index = new GridID(iproc);
        grids[iproc] = CreateGridGeneral(UNSTRUCTGRID, index, 0, GetDim());
        
        ComputeGridForPartition(static_cast<UnstructGrid*>(grids[iproc]), iproc);
        
        ComputeGlobalToLocalCellToNodeMapping(static_cast<UnstructGrid*>(grids[iproc]), iproc);
    }

    if (!IsParallelPartition())
    {
        ModifyInterpointInfoPart(grids);
    }
    else
    {
        PrintToWindow("Warning: Haven't consider inter-point for parallel partition!\n");
    }
    
    SetInterfaceInfoPartInOrgBC(grids);
    
    BuildInterfaceLinkForParallelPartition(grids);
    
    return grids;
}

void Pre_UnstructGridPartition::FreeLocal2All()
{
    delete [] node2all;
    delete [] face2all;
    delete [] cell2all;

    node2all = 0;
    face2all = 0;
    cell2all = 0;
}

LIB_EXPORT void Pre_UnstructGridPartition::PartitionGridSerial(const string &out_grid_file)
{
    Init();

    grids = this->CreateAllPartGrid();
    
    BuildInterfaceLink(grids, maxProcessor);
}

LIB_EXPORT void Pre_UnstructGridPartition::PartitionGridParallel(const string &out_grid_file)
{
    Init();

    GetVertexDist();

    grids = CreateAllPartGrid();

    // The following function is need only for serial partition.\n
    // This has been don in BuildInterfaceLinkForParallelPartition.
    //BuildInterfaceLink(grids, maxProcessor); 
}

void Pre_UnstructGridPartition::GetVertexDist()
{
    //if(parallelPartitionMethod == METIS_METHOD) return;

#if IDXTYPEWIDTH == 32
    PrintToWindow("  Using 32-bit METIS/PARMETIS lib ...\n");
#else
    PrintToWindow("  Using 64-bit METIS/PARMETIS lib ...\n");
#endif

    int nZones = PHMPI::GetNumberofGlobalZones();
    int nProc  = PHMPI::GetNumberOfProcessor();
    ASSERT(nProc == nZones);
    int *nTCellofEachZone = new int[nProc];

    FillnTCellofEachZone(nTCellofEachZone);

    vtxdist = new idx_t[nZones + 1];
    vtxdist[0] = 0;

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        vtxdist[iZone + 1] = vtxdist[iZone] + nTCellofEachZone[iZone];
    }

    if (PHMPI::GetCurrentProcessorID() == PHMPI::GetServerProcessorID())
    {
        int minNTCell = 10000000, maxNTCell = 0;
        for (int iZone = 0; iZone < nZones; ++ iZone)
        {
            minNTCell = MIN(minNTCell, nTCellofEachZone[iZone]);
            maxNTCell = MAX(maxNTCell, nTCellofEachZone[iZone]);
        }

        ostringstream oss;
        oss << "  Before partition, Min && Max number of cell : " << minNTCell << ", " << maxNTCell << endl;
        oss << "                    Min to Max nTCell ratio : " << (minNTCell * 1.0) / (maxNTCell * 1.0) << endl;
        PrintToWindow(oss);
        WriteLogFile(oss);
    }

    delete [] nTCellofEachZone;
}

void Pre_UnstructGridPartition::CreatePartition()
{
    int nTotalCell = gridSource->GetNTotalCell();
    int nTotalFace = gridSource->GetNTotalFace();
    int nTotalNode = gridSource->GetNTotalNode();

    cell_mark = new int[nTotalCell];    //all 2 local
    face_mark = new int[nTotalFace];    //all 2 local
    node_mark = new int[nTotalNode];    //all 2 local

    if (!IsParallelPartition())
    {
        part = CreatePartitionSerial(gridSource, maxProcessor);
    }
    else
    {
        part = CreatePartitionParallel(gridSource, maxProcessor, vtxdist);

        // The following section is used in parallel partition.
        // This pointers are used in communication.
        int nBoundFace = gridSource->GetNBoundFace();
        interfaceIndexInBC_send    = new int[nBoundFace];
        interfaceIndexInBC_recv    = new int[nBoundFace];
        partIndexInNewPartion_recv = new int[nBoundFace];
        for (int iface = 0; iface < nBoundFace; ++ iface)
        {
            interfaceIndexInBC_send[iface]    = -1;
            interfaceIndexInBC_recv[iface]    = -1;
            partIndexInNewPartion_recv[iface] = -1;
        }
    }

    SetProcessorList();
}

idx_t * Pre_UnstructGridPartition::CreatePartitionSerial(UnstructGrid *grid, int maxProcessor)
{
    if (maxProcessor < 2)
    {
        TK_Exit::ExceptionExit("The number of partitions should be greater than 1!\n");
    }

    int nTotalCell = grid->GetNTotalCell();
    int nTotalFace = grid->GetNTotalFace();
    int nBoundFace = grid->GetNBoundFace();

    idx_t *xadj   = new idx_t[nTotalCell + 1];
    idx_t *adjncy = new idx_t[2 * (nTotalFace - nBoundFace)];
    idx_t *part   = new idx_t[nTotalCell];

    for (int icell = 0; icell < nTotalCell; ++ icell)
    {
        part[icell] = -1;
    }

    Get_Xadj_Adjncy(grid, xadj, adjncy);
    MetisGrid(nTotalCell, xadj, adjncy, maxProcessor, part);

    delete [] xadj;
    delete [] adjncy;

    return part;
}

idx_t * Pre_UnstructGridPartition::CreatePartitionParallel(UnstructGrid *grid, int maxProcessor, idx_t *vtxdist)
{
    if (maxProcessor < 2)
    {
        TK_Exit::ExceptionExit("The number of partitions should be greater than 1!\n");
    }

    int nTotalCell = grid->GetNTotalCell();
    int nTotalFace = grid->GetNTotalFace();
    int nBoundFace = grid->GetNBoundFace();
    int nIFace = grid->GetNIFace();

    idx_t *xadj   = new idx_t[nTotalCell + 1];
    idx_t *adjncy = new idx_t[2 * (nTotalFace - nBoundFace + nIFace)];
    idx_t *part   = new idx_t[nTotalCell];
    for (int icell = 0; icell < nTotalCell; ++ icell)
    {
        part[icell] = -1;
    }

    int *neighborCell = new int[nIFace];
    CommunicateNeiborCell(neighborCell);

    Get_Xadj_Adjncy(grid, xadj, adjncy, neighborCell, vtxdist);
    
    if (parallelPartitionMethod == PARMETIS_METHOD)
    {
        PrintToWindow("Parallel grid partition by ParMETIS ...\n");

        ParMetisGrid(nTotalCell, xadj, adjncy, maxProcessor, part, vtxdist);
    }
    else if (parallelPartitionMethod == METIS_METHOD || parallelPartitionMethod == METIS_OPENMP_METHOD)
    {
        PrintToWindow("Parallel grid partition by METIS ...\n");

        MetisGrid(nTotalCell, xadj, adjncy, nPartInEachProcIfMetis, part);
        ConvertPartNumberToGlobal(part, nTotalCell);
    }

    delete [] xadj;
    delete [] adjncy;
    delete [] neighborCell;

    return part;
}

void Pre_UnstructGridPartition::MetisGrid(idx_t nTotalCell, idx_t *xadj, idx_t *adjncy, idx_t nparts, idx_t *part)
{
    WriteLogFile("Now beginning partition graph by METIS!");

    idx_t ncon     = 1;
    idx_t *vwgt    = 0;
    idx_t *vsize   = 0;
    idx_t *adjwgt  = 0;
    real_t *tpwgts = 0;
    real_t *ubvec  = 0;
    idx_t options[METIS_NOPTIONS];
    idx_t objval = 0;

    METIS_SetDefaultOptions(options);
  
    PrintToWindow("Now begining partition graph!\n");
  
    if (nparts > 8)
    {
        PrintToWindow("Using K-way Partitioning!\n");
        //METIS4.0
        //METIS_PartGraphKway(&nTotalCell,xadj,adjncy,vwgt,adjwgt,&wgtflag, &numflag, 
        //                         &nparts,options,&edgecut,part);

        //METIS5.0
       METIS_PartGraphKway(&nTotalCell, &ncon, xadj, adjncy, vwgt, vsize, adjwgt, 
                           &nparts, tpwgts, ubvec, options, &objval, part);

        //METIS PartGraphVKway (int *n, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *vsize, int *wgtflag,
        //int *numflag, int *nparts, int *options, int *volume, idxtype *part);

        //int METIS_PartGraphKway(idx_t *nvtxs, idx_t *ncon, idx_t *xadj, idx_t *adjncy, 
        //                        idx_t *vwgt, idx_t *vsize, idx_t *adjwgt, idx_t *nparts, 
        //                        real_t *tpwgts, real_t *ubvec, idx_t *options, idx_t *objval, 
        //                        idx_t *part)
    }
    else
    {
        PrintToWindow("Using Recursive Partitioning!\n");
        //METIS4.0
        //METIS_PartGraphRecursive(&nTotalCell,xadj,adjncy,vwgt,adjwgt,&wgtflag, &numflag, 
        //                         &nparts,options,&edgecut,part);

        METIS_PartGraphRecursive(&nTotalCell, &ncon, xadj, adjncy, vwgt, vsize, adjwgt, 
                                 &nparts, tpwgts, ubvec, options, &objval, part);

        //METIS5.0
        //int METIS_PartGraphRecursive(idx_t *nvtxs, idx_t *ncon, idx_t *xadj, 
        //                             idx_t *adjncy, idx_t *vwgt, idx_t *vsize, idx_t *adjwgt, 
        //                             idx_t *nparts, real_t *tpwgts, real_t *ubvec, idx_t *options, 
        //                             idx_t *objval, idx_t *part)
    }

    ostringstream oss;
    oss << "The interface number: " << objval << endl; 
    oss << "Partition is finished!\n";
    PrintToWindow(oss);
}

void Pre_UnstructGridPartition::ParMetisGrid(int nTotalCell, idx_t *xadj, idx_t *adjncy, idx_t nparts, idx_t *part, idx_t *vtxdist)
{
    WriteLogFile("Now beginning parallel partition graph by ParMETIS!");
    PrintToWindow("Now beginning parallel partition graph by ParMETIS!\n");

    idx_t ncon     = 1;
    idx_t *vwgt    = 0;
    idx_t *adjwgt  = 0;
    real_t *tpwgts = new real_t[nparts];
    for (int iPart = 0; iPart < nparts; ++ iPart)
    {
        tpwgts[iPart] = static_cast<real_t> (1.0 / nparts);
    }

    //! parmetisBalance: is used to specify the imbalance tolerance.
    //! 1     : perfect balance;
    //! nparts: perfect imbalance;
    //! 1.05  :recommended.
    RDouble parmetisBalance = GlobalDataBase::GetDoubleParaFromDB("parmetisBalance");

    real_t *ubvec = new real_t[ncon];
    for (int icon = 0; icon < ncon; ++ icon)
    {
        ubvec[icon] = static_cast<real_t> (parmetisBalance);    //modified by LIJIN,2019-2-19
    }
    idx_t options[METIS_NOPTIONS];
    options[0] = 0;
    idx_t wgtflag = 0;
    idx_t numflag = 0;
    idx_t edgecut = 0;

    MPI_Comm comm;
    MPI_Comm_dup(MPI_COMM_WORLD, &comm);
    
    if (PHMPI::GetNumberOfProcessor() == 1)
    {
        TK_Exit::ExceptionExit("Error: number of processors must great than one for parallel partition!\n");
    }
    else
    {
        ParMETIS_V3_PartKway(vtxdist, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag, &ncon,
            &nparts, tpwgts, ubvec, options, &edgecut, part, &comm);
    }

    delete [] ubvec;
    delete [] tpwgts;

    WriteLogFile("Partition is finished by Parallel Metis!");
    PrintToWindow("Partition is finished by Parallel Metis!\n");
}

void Pre_UnstructGridPartition::Get_Xadj_Adjncy(UnstructGrid *grid, idx_t *xadj, idx_t *adjncy)
{   
    int  nTotalCell = grid->GetNTotalCell();
    vector<int> *c2c = grid->GetCell2Cell();

    xadj[0]   = 0;
    int count = 0;
    for (int icell = 0; icell < nTotalCell; ++ icell)
    {
        xadj[icell+1] = xadj[icell] + c2c[icell].size();
        for (std::size_t j = 0; j < c2c[icell].size(); ++ j)
        {
            adjncy[count ++] = c2c[icell][j];
        }
    }
}

void Pre_UnstructGridPartition::Get_Xadj_Adjncy(UnstructGrid *grid, idx_t *xadj, idx_t *adjncy, int *neighborCell, idx_t *vtxdist)
{   
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotalFace = grid->GetNTotalFace();
    int nIFace = grid->GetNIFace();
    int *left_cell_of_face = grid->GetLeftCellOfFace();
    int *right_cell_of_face = grid->GetRightCellOfFace();

    InterfaceInfo *interfaceInfo = grid->GetInterfaceInfo();
    int *interFace2ZoneID = interfaceInfo->GetInterFace2ZoneID();
    int *interFace2BoundaryFace = interfaceInfo->GetInterFace2BoundaryFace();

    idx_t startGlobalIndex = 0;

    // First: compute the cell-to-cell connectivity relationship, include interface.
    vector<idx_t> *c2c = new vector<idx_t>[nTotalCell];
    idx_t le, re;

    // If boundary is an INTERFACE
    if (parallelPartitionMethod == PARMETIS_METHOD)
    {
        for (int iface = 0; iface < nIFace; ++ iface)
        {
            int neighborZone = interFace2ZoneID[iface];
            int bcFace = interFace2BoundaryFace[iface];
            le = left_cell_of_face[bcFace];
            re = neighborCell[iface];
            re += vtxdist[neighborZone];
            c2c[le].push_back(re);
        }

        int myid = PHMPI::GetCurrentProcessorID();
        startGlobalIndex = vtxdist[myid];
    }

    // Interior faces
    for (int iface = nBoundFace; iface < nTotalFace; ++ iface)
    {
        le = left_cell_of_face [iface];
        re = right_cell_of_face[iface];

        c2c[le].push_back(re + startGlobalIndex);
        c2c[re].push_back(le + startGlobalIndex);
    }

    // Second: compute the xadj and adjncy.
    xadj[0]   = 0;
    int count = 0;
    for (int icell = 0; icell < nTotalCell; ++ icell)
    {
        xadj[icell+1] = xadj[icell] + c2c[icell].size();
        for (std::size_t j = 0; j < c2c[icell].size(); ++ j)
        {
            adjncy[count ++] = c2c[icell][j];
        }
    }

    // Free memory.
    for (int icell = 0; icell < nTotalCell; ++ icell)
    {
        c2c[icell].clear();
    }
    delete [] c2c;
}

void Pre_UnstructGridPartition::ConvertPartNumberToGlobal(idx_t *part, int nTotalCell)
{
    // Each Processor owns one zone, and each one is partition to be same number of zones.
    // so, the start zone index in this processor(zone) is start from:
    // number of pre-zones * nPartInEachProcIfMetis in each original zones.
    int myid = PHMPI::GetCurrentProcessorID();
    int startZoneID = myid * this->nPartInEachProcIfMetis; 

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        part[iCell] += startZoneID;
    }
}

void Pre_UnstructGridPartition::CommunicateNeiborCell(int *neighborCell)
{
    WriteLogFile("Start to Communicate Neighbor Cells");

    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    // Each processor runs along the same order, which is from 0 to nZones.
    // Each simulation zone is looped over to judge if the communication is necessary.
    for (int izone = 0; izone < nZones; ++ izone)
    {
        ZoneNeighbor *zoneNeighbor = zoneConnectivity->GetZoneNeighbor(izone);
        std::size_t numberOfNeighbor = zoneNeighbor->GetNumberOfNeighbor();
        for (std::size_t ineighbor = 0; ineighbor < numberOfNeighbor; ++ ineighbor)
        {
            int neighborZoneIndex = zoneNeighbor->GetZoneIndexOfNeighbor(ineighbor);
            SwapInterface(neighborCell, izone, neighborZoneIndex);
        }
    }

    WriteLogFile("End to Communicate Neighbor Cells");
}

void Pre_UnstructGridPartition::SwapInterface(int *neighborCell, int izone, int jzone)
{
    using namespace PHMPI;

    //! izone is the current zone, and jzone is the target zone.
    int send_proc = GetZoneProcessorID(izone);
    int recv_proc = GetZoneProcessorID(jzone);

    int tag = izone;
    int myid = GetCurrentProcessorID();

    DataContainer *cdata = new DataContainer();
    cdata->moveto_begin();

    if (myid == send_proc)
    {
        UnstructGrid  *grid  = UnstructGridCast(GetGrid(izone, 0));
        InterfaceInfo *interfaceInfo = grid->GetInterfaceInfo();
        if (interfaceInfo)
        {
            int ineighbor = interfaceInfo->FindIthNeighbor(jzone);
            int nIFaceOfNeighbor = interfaceInfo->GetNIFaceOfNeighbor(ineighbor);
            int *faceIndexForSend = interfaceInfo->GetFaceIndexForSend(ineighbor);
            int iface,sourceCell;
            for (int ifaceLocal = 0; ifaceLocal < nIFaceOfNeighbor; ++ ifaceLocal)
            {
                iface = faceIndexForSend[ifaceLocal];
                grid->GetSourceIndex(iface, 1, sourceCell);
                cdata->write(&sourceCell, sizeof(int));
            }
        }
    }

    PH_Trade(cdata, send_proc, recv_proc, tag); 

    if (myid == recv_proc)
    {
        UnstructGrid *grid = UnstructGridCast(GetGrid(jzone, 0));
        InterfaceInfo *interfaceInfo = grid->GetInterfaceInfo();
        if (interfaceInfo)
        {
            int ineighbor = interfaceInfo->FindIthNeighbor(izone);
            int iface;
            int nIFaceOfNeighbor = interfaceInfo->GetNIFaceOfNeighbor(ineighbor);
            int *faceIndexForRecv = interfaceInfo->GetFaceIndexForRecv(ineighbor);
            cdata->moveto_begin();
            for (int ifaceLocal = 0; ifaceLocal < nIFaceOfNeighbor; ++ ifaceLocal)
            {
                iface = faceIndexForRecv[ifaceLocal];
                int cellIndex;
                cdata->read(&cellIndex, sizeof(int));
                neighborCell[iface] = cellIndex;
            }
        }
    }

    delete cdata;
}

void Pre_UnstructGridPartition::SetProcessorList()
{
    if (!IsParallelPartition()) return;

    int nTotalCell = this->gridSource->GetNTotalCell();
    
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        int partIndex = static_cast<int>(part[iCell]);

        processorList.insert(partIndex);
    }

    ostringstream oss;
    oss << "  Part index in processor " << PHMPI::GetCurrentProcessorID() << ": ";
    for (set<int>::iterator iter = processorList.begin(); iter != processorList.end(); ++ iter)
    {
        int partID = *iter;
        if (partID < 0)
        {
            TK_Exit::ExceptionExit("Error: partID < 0!");
        }
        oss << partID << " ";
    }
    oss << endl;
    WriteLogFile(oss.str());
    
    //! Compute the global zone index of each cell.
    int nLocalZones = this->GetNumberofLocalZones();

    int nZones = PHMPI::GetNumberofGlobalZones();
    int nProc  = PHMPI::GetNumberOfProcessor();
    ASSERT(nProc == nZones);
    int myid = PHMPI::GetCurrentProcessorID();
    int *nPartOfEachProcessor = new int[nZones];

    if (nProc == nZones)
    {
        //! Each processors has ONE zone.
        int *localSizeTemp = new int [nProc];
        SetField(localSizeTemp, 1, nProc);
        int **nPartListTemp = new int * [nProc];

        PHMPI::FillGobalDataCollectBcast(&nLocalZones, localSizeTemp, nPartListTemp);

        for (int iProc = 0; iProc < nProc; ++ iProc)
        {
            nPartOfEachProcessor[iProc] = nPartListTemp[iProc][0];
        }

        DelPointer2(nPartListTemp);
        DelPointer(localSizeTemp);
    }
    else
    {
        TK_Exit::ExceptionExit("number of processors != number of zones.");
        //PHMPI::FillGobalData(nLocalZones, nPartOfEachProcessor);
    }

    part2LocalZoneIndex  = new int[maxProcessor];
    part2GlobalZoneIndex = new int [maxProcessor];
    for (int iproc = 0; iproc < maxProcessor; ++ iproc)
    {
        part2LocalZoneIndex[iproc]  = -1;
        part2GlobalZoneIndex[iproc] = -1;
    }

    int localIndex = 0;
    for (int iproc = 0; iproc < maxProcessor; ++ iproc)
    {
        if (processorList.find(iproc) != processorList.end())
        {
            part2LocalZoneIndex[iproc] = localIndex ++;
        }
    }
    
    int count = 0;
    for (int iProc = 0; iProc < myid; ++ iProc)
    {
        count += nPartOfEachProcessor[iProc];
    }
    for (set<int>::iterator iter = processorList.begin(); iter != processorList.end(); ++ iter)
    {
        int partIndex = *iter;
        int localIndex = part2LocalZoneIndex[partIndex];
        part2GlobalZoneIndex[partIndex] = localIndex + count;
    }

    delete [] nPartOfEachProcessor;
}

int Pre_UnstructGridPartition::GetNumberofLocalZones()
{
    int nLocalZones = maxProcessor;
    if (parallelPartitionMethod == PARMETIS_METHOD)
    {
        nLocalZones = static_cast<int>(processorList.size());
    }
    else if (parallelPartitionMethod == METIS_METHOD || parallelPartitionMethod == METIS_OPENMP_METHOD)
    {
        nLocalZones = static_cast<int>(processorList.size());
    }

    return nLocalZones;
}

bool Pre_UnstructGridPartition::DoseProcExist(int iproc)
{
    if (!IsParallelPartition()) return true;

    if (processorList.find(iproc) == processorList.end())
    {
        return false;
    }
    else
    {
        return true;
    }
}

int Pre_UnstructGridPartition::GetMaxProcessorIndexInThisZone()
{
    int nTotalCell = this->gridSource->GetNTotalCell();
    int maxNumber = 0;
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        int partIndex = static_cast<int>(part[iCell]);
        maxNumber = MAX(maxNumber, partIndex);
    }
    return (maxNumber + 1);
}

void Pre_UnstructGridPartition::ComputeCellMark()
{
    int nTotalCell = gridSource->GetNTotalCell();

    int  maxProcessor = GetMaxProcessorIndexInThisZone();

    int *count_array = new int[maxProcessor];

    for (int i = 0; i < maxProcessor; ++ i)
    {
        count_array[i] = 0;
    }

    for (int icell = 0; icell < nTotalCell; ++ icell)
    {
        cell_mark[icell] = count_array[ part[icell] ] ++;
    }

    delete [] count_array;
}


//! Set the interface information of the interfaces that belong to the original BC faces.\n
//! The interfaces that belong to original interior faces have been set by SetInterfaceInfoPart().
void Pre_UnstructGridPartition::SetInterfaceInfoPartInOrgBC(Grid **grids)
{
    if (!IsParallelPartition()) return;

    if (npartmethod != 1) return;
    
    CommunicateInterfaceTopo();
    
    for (int iproc = 0; iproc < maxProcessor; ++ iproc)
    {
        if (!DoseProcExist(iproc))
        {
            grids[iproc] = 0;
            continue;
        }

        SetInterfaceInfoPartInOrgBC(iproc, UnstructGridCast(grids[iproc]));
    }
}

void Pre_UnstructGridPartition::SetInterfaceInfoPartInOrgBC(int izone, UnstructGrid *gridpart)
{
    if (npartmethod != 1) return;
    
    int nBoundFace = gridSource->GetNBoundFace();
    int nTotalCell = gridSource->GetNTotalCell();
    int *left_cell_of_face  = gridSource->GetLeftCellOfFace();
    int *right_cell_of_face = gridSource->GetRightCellOfFace();
    int le, re;
    CompositeUnstructBC **bcr = gridSource->GetBCRecord();

    //! Re-compute bc_face_mark for each sub-partition.
    //! The original face_mark have been overwritten.
    int *bc_face_mark = new int [nBoundFace];
    for (int iface = 0; iface < nBoundFace; ++ iface)
    {
        bc_face_mark[iface] = -2;
    }
    
    for (int iface = 0; iface < nBoundFace; ++ iface)
    {
        le = left_cell_of_face [iface];
        re = right_cell_of_face[iface];
        if (part[le] == izone)
        {
            bc_face_mark[iface] = -1;
        }
        else if (re < nTotalCell && part[re] == izone)
        {
            bc_face_mark[iface] = -1;
        }
    }
    
    //! Only deal with interfaces that belong to original BC faces.
    InterfaceInfo *local_iinfo = gridpart->GetInterfaceInfo();
    int *local_interFace2ZoneID       = local_iinfo->GetInterFace2ZoneID();
    int *local_interFace2InterFaceID  = local_iinfo->GetInterFace2InterFaceID();
    int *local_interFace2BoundaryFace = local_iinfo->GetInterFace2BoundaryFace();

    int i_nTFace = 0;
    int i_nIFace = 0;
    
    for (int iface = 0; iface < nBoundFace; ++ iface)
    {
        if (bc_face_mark[iface] == -1)
        {
            bc_face_mark[iface] = i_nTFace ++;

            int bcType = bcr[iface]->GetKey();
            if (bcType == PHENGLEI::INTERFACE)
            {
                local_interFace2ZoneID[i_nIFace]       = partIndexInNewPartion_recv[iface];
                local_interFace2InterFaceID[i_nIFace]  = interfaceIndexInBC_recv[iface];
                local_interFace2BoundaryFace[i_nIFace] = bc_face_mark[iface];

                i_nIFace ++;
            }
        }
    }
    delete [] bc_face_mark;
}

void Pre_UnstructGridPartition::BuildInterfaceLinkForParallelPartition(Grid **grids)
{
    if (!IsParallelPartition()) return;

    for (int izone = 0; izone < maxProcessor; ++ izone)
    {
        if (!grids[izone]) continue;

        InterfaceInfo *interfaceInfo = grids[izone]->GetInterfaceInfo();
        if (!interfaceInfo) continue;

        int *interFace2ZoneID      = interfaceInfo->GetInterFace2ZoneID();
        int *interFace2InterFaceID = interfaceInfo->GetInterFace2InterFaceID();

        UnstructGrid *grid = UnstructGridCast(grids[izone]);

        int nBoundFace = grid->GetNBoundFace();
        int nIFace = interfaceInfo->GetNIFace();
        int nphysi = nBoundFace - nIFace;    //! number of physical boundary = total number boundary face - number of interface.

        int *left_cell_of_face  = grid->GetLeftCellOfFace();
        int *right_cell_of_face = grid->GetRightCellOfFace();

        bool *isNeighborInterfaceFound = interfaceInfo->GetIsNeighborInterfaceFound();

        int nbz, le, re, ie;

        //! Interfaces count from 0.
        for (int iface = 0; iface < nIFace; ++ iface)
        {
            if (iface < nIFaceInOrgBC[izone])
            {
                //! The interfaces in the original BC \n
                //! have been built in function SetInterfaceInfoPartInOrgBC.
                continue;
            }

            if (isNeighborInterfaceFound[iface] == true)
            {
                //! The neighbor of this face has been found and paired, so skip it.
                continue;
            }
            
            //! Neighbor zone index.
            nbz = interFace2ZoneID[iface];
            le  = left_cell_of_face [iface + nphysi];
            re  = right_cell_of_face[iface + nphysi];
            ie  = MAX(le, re);
            if (nbz >= izone)
            {
                UnstructGrid *grid_t = UnstructGridCast(grids[nbz]);
                if (!FindMatchPart(grid_t, izone, iface, ie, interFace2InterFaceID))
                {
                    TK_Exit::ExceptionExit("!FindMatchPart\n");
                }
                else
                {
                    isNeighborInterfaceFound[iface] = true;
                }
            }
        }
    }

    //! Map the local part index to global zone index.
    for (int izone = 0; izone < maxProcessor; ++ izone)
    {
        if (!grids[izone]) continue;

        InterfaceInfo *interfaceInfo = grids[izone]->GetInterfaceInfo();
        if (!interfaceInfo) continue;

        int *interFace2ZoneID = interfaceInfo->GetInterFace2ZoneID();
        int nIFace = interfaceInfo->GetNIFace();
        for (int iface = 0; iface < nIFace; ++ iface)
        {
            if (iface < nIFaceInOrgBC[izone])
            {
                //! The interfaces in the original BC \n
                //! have been built in function SetInterfaceInfoPartInOrgBC.
                continue;
            }
            
            int localPartIndex = interFace2ZoneID[iface];
            interFace2ZoneID[iface] = part2GlobalZoneIndex[localPartIndex];
        }
    }
}

void Pre_UnstructGridPartition::CommunicateInterfaceTopo()
{
    if (!IsParallelPartition()) return;

    int nZones = PHMPI::GetNumberofGlobalZones();

    //! Each processor runs along the same order, which is from 0 to nZones.
    //! Each simulation zone is looped over to judge if the communication is necessary.
    for (int izone = 0; izone < nZones; ++ izone)
    {
        ZoneNeighbor *zoneNeighbor = zoneConnectivity->GetZoneNeighbor(izone);
        std::size_t numberOfNeighbor = zoneNeighbor->GetNumberOfNeighbor();

        for (std::size_t ineighbor = 0; ineighbor < numberOfNeighbor; ++ ineighbor)
        {
            int neighborZoneIndex = zoneNeighbor->GetZoneIndexOfNeighbor(ineighbor);

            SwapInterfaceTopo(izone, neighborZoneIndex);
        }
    }
}

void Pre_UnstructGridPartition::SwapInterfaceTopo(int izone, int jzone)
{
    using namespace PHMPI;

    //! izone is the current zone, and jzone is the target zone.
    int send_proc = GetZoneProcessorID(izone);
    int recv_proc = GetZoneProcessorID(jzone);

    int tag = izone;
    int myid = GetCurrentProcessorID();

    DataContainer *cdata = new DataContainer();

    if (myid == send_proc)
    {
        UnstructGrid  *grid = UnstructGridCast(GetGrid(izone, 0));
        InterfaceInfo *interfaceInfo = grid->GetInterfaceInfo();
        if (interfaceInfo)
        {
            int ineighbor = interfaceInfo->FindIthNeighbor(jzone);
            int nIFaceOfNeighbor = interfaceInfo->GetNIFaceOfNeighbor(ineighbor);
            int *faceIndexForSend = interfaceInfo->GetFaceIndexForSend(ineighbor);
            int *interFace2BoundaryFace = interfaceInfo->GetInterFace2BoundaryFace();
            int iface, le, partIndex, localInterfaceIndex, bcFaceIndex;
            cdata->moveto_begin();
            for (int ifaceLocal = 0; ifaceLocal < nIFaceOfNeighbor; ++ ifaceLocal)
            {
                iface = faceIndexForSend[ifaceLocal];
                bcFaceIndex = interFace2BoundaryFace[iface];    // BC face index in old grid.
                grid->GetSourceIndex(iface, 1, le);

                partIndex = part2GlobalZoneIndex[part[le]];    // global zone index in new partition.
                localInterfaceIndex = interfaceIndexInBC_send[bcFaceIndex]; // local interface index in new partition.

                cdata->write(&partIndex, sizeof(int));              // interFace2ZoneID
                cdata->write(&localInterfaceIndex, sizeof(int));    // interFace2InterFaceID
            }
        }
    }

    PH_Trade(cdata, send_proc, recv_proc, tag);   

    if (myid == recv_proc)
    {
        UnstructGrid *grid = UnstructGridCast(GetGrid(jzone, 0));
        InterfaceInfo *interfaceInfo = grid->GetInterfaceInfo();
        if (interfaceInfo)
        {
            int *interFace2BoundaryFace = interfaceInfo->GetInterFace2BoundaryFace();

            int ineighbor = interfaceInfo->FindIthNeighbor(izone);
            int nIFaceOfNeighbor = interfaceInfo->GetNIFaceOfNeighbor(ineighbor);
            int *faceIndexForRecv = interfaceInfo->GetFaceIndexForRecv(ineighbor);
            int iface, globalZoneIndexInNeighbor, localInterfaceIndexInNeighbor, bcFaceIndex;
            cdata->moveto_begin();
            for (int ifaceLocal = 0; ifaceLocal < nIFaceOfNeighbor; ++ ifaceLocal)
            {
                iface = faceIndexForRecv[ifaceLocal];
                bcFaceIndex = interFace2BoundaryFace[iface];
                
                cdata->read(&globalZoneIndexInNeighbor, sizeof(int));
                cdata->read(&localInterfaceIndexInNeighbor, sizeof(int));

                partIndexInNewPartion_recv[bcFaceIndex] = globalZoneIndexInNeighbor;
                interfaceIndexInBC_recv[bcFaceIndex] = localInterfaceIndexInNeighbor;
            }
        }
    }

    delete cdata;
}

bool Pre_UnstructGridPartition::FindMatchPart(UnstructGrid *grid, int izone_t, int iface_t, int ie_t, int *interface_interfaceid_t)
{
    bool found = false;

    InterfaceInfo *interfaceInfo = grid->GetInterfaceInfo();
    if (!interfaceInfo) return found;

    int nBoundFace = grid->GetNBoundFace();
    int nIFace = interfaceInfo->GetNIFace();
    int nphysi = nBoundFace - nIFace;    //! number of physical boundary = total number boundary face - number of interface.

    int *left_cell_of_face  = grid->GetLeftCellOfFace();
    int *right_cell_of_face = grid->GetRightCellOfFace();

    int *interFace2ZoneID      = interfaceInfo->GetInterFace2ZoneID();
    int *interFace2InterFaceID = interfaceInfo->GetInterFace2InterFaceID();
    
    bool *isNeighborInterfaceFound = interfaceInfo->GetIsNeighborInterfaceFound();

    int le, re, ie;

    //! Interfaces count from 0.
    for (int iface = 0; iface < nIFace; ++ iface)
    {
        le = left_cell_of_face [iface + nphysi];
        re = right_cell_of_face[iface + nphysi];
        ie = MAX(le, re);

        if (isNeighborInterfaceFound[iface] == true)
        {
            //! The neighbor of this face has been found and paired, so skip it.
            continue;
        }

        //! The interface iface target zone==izone.
        if ((interFace2ZoneID[iface]      == izone_t) &&
            (interFace2InterFaceID[iface] == ie_t) && 
            (ie  == interface_interfaceid_t[iface_t]))
        {
            //! The jface of current zone and the iface corresponding to target zone are all count from 0.
            interFace2InterFaceID  [iface  ] = iface_t;
            interface_interfaceid_t[iface_t] = iface;

            isNeighborInterfaceFound[iface] = true;

            found = true;
            break;
        } 
    }

    return found;
}

//! For each interface in the current zone,\n
//! Find out the interface-index in the neighbor zone.
void Pre_UnstructGridPartition::BuildInterfaceLink(Grid **grids, int maxProcessor)
{
    if (npartmethod == 1)
    {
        for (int iProc = 0; iProc < maxProcessor; ++ iProc)
        {
            ReconInterfaces(grids, iProc, maxProcessor);
        }
    }
    else
    {
        ConstructGlobalInterfaceLink(grids, maxProcessor);
    }
}

//! For each interface in the current zone,\n
//! Find out the interface-index in the neighbor zone.
void Pre_UnstructGridPartition::ReconInterfaces(Grid **grids, int izone, int maxProcessor)
{
    if (!grids[izone]) return;

    InterfaceInfo *interfaceInfo = grids[izone]->GetInterfaceInfo();
    if (!interfaceInfo) return;

    int *interFace2ZoneID      = interfaceInfo->GetInterFace2ZoneID();
    int *interFace2InterFaceID = interfaceInfo->GetInterFace2InterFaceID();

    UnstructGrid *grid = UnstructGridCast(grids[izone]);

    int nBoundFace = grid->GetNBoundFace();
    int nIFace = interfaceInfo->GetNIFace();
    int nphysi = nBoundFace - nIFace;    //! number of physical boundary = total number boundary face - number of interface.

    int *left_cell_of_face  = grid->GetLeftCellOfFace();
    int *right_cell_of_face = grid->GetRightCellOfFace();

    bool *isNeighborInterfaceFound = interfaceInfo->GetIsNeighborInterfaceFound();

    int nbz, le, re, ie;

    //! Interfaces count from 0.
    for (int iface = 0; iface < nIFace; ++ iface)
    {
        if (isNeighborInterfaceFound[iface] == true)
        {
            //! The neighbor of this face has been found and paired, so skip it.
            continue;
        }

        //! Neighbor zone index.
        nbz = interFace2ZoneID[iface];
        le  = left_cell_of_face [iface + nphysi];
        re  = right_cell_of_face[iface + nphysi];
        ie  = MAX(le, re);

        if (nbz >= izone)
        {
            UnstructGrid *grid_t = UnstructGridCast(grids[nbz]);
            if (!FindMatchPart(grid_t, izone, iface, ie, interFace2InterFaceID))
            {
                TK_Exit::ExceptionExit("!FindMatchPart\n");
            }
            else
            {
                isNeighborInterfaceFound[iface] = true;
            }
        }
    }
}

void Pre_UnstructGridPartition::SetCoorPart(UnstructGrid *gridpart)
{
    int nTotalNode = gridpart->GetNTotalNode();
    RDouble *x = new RDouble[nTotalNode];
    RDouble *y = new RDouble[nTotalNode];
    RDouble *z = new RDouble[nTotalNode];

    gridpart->SetX(x);
    gridpart->SetY(y);
    gridpart->SetZ(z);

    int nTNodeSrc = gridSource->GetNTotalNode();
    RDouble *x_src = gridSource->GetX();
    RDouble *y_src = gridSource->GetY();
    RDouble *z_src = gridSource->GetZ();

    int count = 0;
    for (int inode = 0; inode < nTNodeSrc; ++ inode)
    {
        if (node_mark[inode] > -1)
        {
            x[count] = x_src[inode];
            y[count] = y_src[inode];
            z[count] = z_src[inode];
            ++ count;
        }
    }

    if (count != nTotalNode)
    {
        cout << "error in SetCoorPart\n";
    }
}

void Pre_UnstructGridPartition::SetNodeNumberOfEachFaceAndFace2Node(UnstructGrid *gridpart)
{
    int *node_number_of_each_face_src = gridSource->GetNodeNumberOfEachFace();
    int *face2node_src = gridSource->GetFace2Node();

    int nTotalFace = gridpart->GetNTotalFace();
    int *node_number_of_each_face = new int[nTotalFace];
    gridpart->SetNodeNumberOfEachFace(node_number_of_each_face);

    //! Extract node_number_of_each_face from the original unstructured grid.
    int count = 0;
    for (int iface = 0; iface < nTotalFace; ++ iface)
    {
        int faceid = face2all[iface];
        node_number_of_each_face[iface] = node_number_of_each_face_src[faceid];
        count += node_number_of_each_face[iface];
    }

    int *face2node = new int[count];
    gridpart->SetFace2Node(face2node);

    int *face2NodeSubscript = gridSource->GetFace2NodeSubscript();

    count = 0;
    for (int iface = 0; iface < nTotalFace; ++ iface)
    {
        //! faceid is the index of global face.
        int faceid = face2all[iface];
        for (int jnode = face2NodeSubscript[faceid]; jnode < face2NodeSubscript[faceid + 1]; ++ jnode)
        {
            face2node[count ++] = node_mark[face2node_src[jnode]];
        }
    }
}

void Pre_UnstructGridPartition::SetFace2CellAndBC(int izone, UnstructGrid *gridpart)
{
    int nBFace_src = gridSource->GetNBoundFace();
    int *left_cell_of_face_src  = gridSource->GetLeftCellOfFace();
    int *right_cell_of_face_src = gridSource->GetRightCellOfFace();

    int nTotalFace = gridpart->GetNTotalFace();
    int nBoundFace = gridpart->GetNBoundFace();

    //! Face cell.
    int *left_cell_of_face  = new int[nTotalFace]();
    int *right_cell_of_face = new int[nTotalFace]();
    gridpart->SetLeftCellOfFace (left_cell_of_face);
    gridpart->SetRightCellOfFace(right_cell_of_face);

    CompositeUnstructBC **bcr_src = gridSource->GetBCRecord();

    CompositeUnstructBC **bcr = new CompositeUnstructBC * [nBoundFace];
    gridpart->SetBCRecord(bcr);

    int le_src = 0, re_src = 0, le = 0, re = 0, bctype = 0;
    for (int iface = 0; iface < nBoundFace; ++ iface)
    {
        int faceid = face2all[iface];

        le_src = left_cell_of_face_src [faceid];
        re_src = right_cell_of_face_src[faceid];

        string bcName = " ";
        if (faceid < nBFace_src)
        {       
            re = - 1;
            le = cell_mark[le_src];
            bctype = bcr_src[ faceid ]->GetKey();
            bcName = bcr_src[ faceid ]->GetBCName();
        }
        else
        {
            bctype = -1;
            bcName = "Interface";
            //! int face.
            if (part[le_src] == izone)
            {
                le = cell_mark[le_src];
                re = - 1;
            }
            else if (part[re_src] == izone)
            {
                re = cell_mark[re_src];
                le = - 1;
            }
            else
            {
                cout << "error in SetLinkPart\n";
            }
        }

        bcr[iface] = new CompositeUnstructBC();
        bcr[iface]->SetKey(bctype);
        bcr[iface]->SetBCName(bcName);

        left_cell_of_face [iface] = le;
        right_cell_of_face[iface] = re;
    }

    for (int iface = nBoundFace; iface < nTotalFace; ++ iface)
    {
        int faceid = face2all[iface];
        le_src = left_cell_of_face_src[faceid];
        re_src = right_cell_of_face_src[faceid];

        left_cell_of_face [iface] = cell_mark[le_src];
        right_cell_of_face[iface] = cell_mark[re_src];
    }
}

void Pre_UnstructGridPartition::SetLinkPart(int izone, UnstructGrid *gridpart)
{
    SetNodeNumberOfEachFaceAndFace2Node(gridpart);
    SetFace2CellAndBC(izone, gridpart);
}

void Pre_UnstructGridPartition::SetInterfaceInfoPart(int izone, UnstructGrid *gridpart)
{
    if (IsParallelPartition())
    {
        SetInterfaceInfoPartForParallelPartition(izone, gridpart);
        return;
    }

    if (npartmethod != 1) return;

    InterfaceInfo *interfaceInfo = gridpart->GetInterfaceInfo();

    int i_nIFace = interfaceInfo->GetNIFace();
    int i_nBFace = gridpart->GetNBoundFace();

    int nTotalFace = gridSource->GetNTotalFace();
    int nBoundFace = gridSource->GetNBoundFace();
    int *left_cell_of_face  = gridSource->GetLeftCellOfFace();
    int *right_cell_of_face = gridSource->GetRightCellOfFace();

    int le, re;

    int *interFace2ZoneID       = interfaceInfo->GetInterFace2ZoneID();
    int *interFace2CellID       = interfaceInfo->GetInterFace2CellID();
    int *interFace2InterFaceID  = interfaceInfo->GetInterFace2InterFaceID();
    int *interFace2BoundaryFace = interfaceInfo->GetInterFace2BoundaryFace();
    int *interFaceDirection     = interfaceInfo->GetInterFaceDirection();

    int nphysi = i_nBFace - i_nIFace;

    //! Here,the interFace2InterFaceID temporarily stores the cells index of neighbor zone.

    int count = 0;
    //! Interfaces in the original interior field.
    for (int iface = nBoundFace; iface < nTotalFace; ++ iface)
    {
        //! face_mark: the corresponding relationship from global face index to local face index.
        //! The following judgment condition means that this face is the interface of the zone.
        //! the index of a face is greater than 0 and less than i_nBFace in new partition zone, and it is the internal face of the original grid.
        //! So it must be interface.
        if (face_mark[iface] < i_nBFace && face_mark[iface] > -1)
        {
            ++ count;
            le = left_cell_of_face [iface];
            re = right_cell_of_face[iface];

            //! ii_local: the index of interface counts from 0.
            int ii_local = face_mark[iface] - nphysi;

            //! interface_interfaceid£ºtemporarily records local index of the neighbor interface cell.
            //! It will be modified to the interface index of neighbor zone.
            if (part[le] == izone)
            {
                interFaceDirection    [ii_local] = 1;
                interFace2ZoneID      [ii_local] = static_cast<int>(part[re]);
                interFace2CellID      [ii_local] = cell_mark[re];
                interFace2InterFaceID [ii_local] = cell_mark[re];
                interFace2BoundaryFace[ii_local] = face_mark[iface];
            }
            else if (part[re] == izone)
            {
                interFaceDirection    [ii_local] = - 1;
                interFace2ZoneID      [ii_local] = static_cast<int>(part[le]);
                interFace2CellID      [ii_local] = cell_mark[le];
                interFace2InterFaceID [ii_local] = cell_mark[le];
                interFace2BoundaryFace[ii_local] = face_mark[iface];
            }
            else
            {
                TK_Exit::ExceptionExit("Error: CAN NOT find neighbor interface cell!");
            }
        }
    }
}

void Pre_UnstructGridPartition::SetInterpointInfoPart(int izone, UnstructGrid *gridpart)
{
    if (npartmethod != 1) return;

    gridpart->computeBoundaryPointLabel();
    map<int, int> boundaryPointLabel_gridpart = gridpart->GetBoundaryPointLabel();
    map<int, int> boundaryPointLabel_gridSource = gridSource->GetBoundaryPointLabel();

    int i_nTNode = gridpart->GetNTotalNode();

    int i_nIPoint = 0;

    //! Compute the number of interpoints in current zone.
    for (int i = 0; i < i_nTNode; i ++)
    {
        int numberofnode2all = node2all[i];
        i_nIPoint += pointPosition2zoneID[numberofnode2all+1] - pointPosition2zoneID[numberofnode2all] -1;
    }

    if (i_nIPoint == 0)
    {
        return;
    }

    InterpointInformation *interpointInformation = new InterpointInformation(i_nIPoint);
    gridpart->SetInterpointInfo(interpointInformation);

    //! Only partition procedure need the data of isNeighborInterpointFound, so new it here.
    bool *isNeighborInterpointFound = new bool [i_nIPoint];
    PHSPACE::SetField(isNeighborInterpointFound, false, i_nIPoint);
    interpointInformation->SetIsNeighborInterpointFound(isNeighborInterpointFound);

    int *interPoint2ZoneID = interpointInformation->GetInterPoint2ZoneID();
    int *interPoint2InterPointID = interpointInformation->GetInterPoint2InterPointID();
    int *interPoint2GlobalPoint = interpointInformation->GetInterPoint2GlobalPoint();
    int *cellNumberOfInterPoint = interpointInformation->GetCellNumberOfInterPoint();
    int *totalZonesOfInterPoint = interpointInformation->GetTotalZonesOfInterPoint();
    int *labelOfInterPoint = interpointInformation->GetLabelOfInterPoint();

    int count = 0;
    for (int inode = 0; inode < i_nTNode; inode ++)
    {
        int numberofnode2all = node2all[inode];
        if (pointPosition2zoneID[numberofnode2all+1] - pointPosition2zoneID[numberofnode2all] == 1) continue;
        for (int i = pointPosition2zoneID[numberofnode2all]; i < pointPosition2zoneID[numberofnode2all+1]; i ++)
        {
            if (point2ZoneID[i] != izone)
            {
                interPoint2ZoneID[count] = point2ZoneID[i];
                interPoint2InterPointID[count] = node2all[inode];    //! Temporarily to stores the index of global points,it is modified to the index of interpoints later.
                interPoint2GlobalPoint[count] = inode;
                cellNumberOfInterPoint[count] = cellNumberForNodeValue[numberofnode2all];
                totalZonesOfInterPoint[count] = pointPosition2zoneID[numberofnode2all+1] - pointPosition2zoneID[numberofnode2all];
                if (boundaryPointLabel_gridpart[inode] != boundaryPointLabel_gridSource[numberofnode2all])
                {
                    labelOfInterPoint[count] = 0;
                }
                else
                {
                    labelOfInterPoint[count] = 1;
                }
                count ++;
            }
        }
    }
}

void Pre_UnstructGridPartition::ModifyInterpointInfoPart(Grid **grids)
{
    for (int iproc = 0; iproc < maxProcessor; ++ iproc)
    {
        UnstructGrid *gridpart = static_cast<UnstructGrid*>(grids[iproc]);
        InterpointInformation *interpointInformation = gridpart->GetInterpointInfo();

        int i_nIPoint = interpointInformation->GetNumberOfInterpoints();
        int *interPoint2ZoneID = interpointInformation->GetInterPoint2ZoneID();
        int *interPoint2InterPointID = interpointInformation->GetInterPoint2InterPointID();

        for (int iPoint = 0; iPoint < i_nIPoint; iPoint ++)
        {
            int neighbor_zone = interPoint2ZoneID[iPoint];
            if (iproc > neighbor_zone) continue;
            UnstructGrid *gridpart_neighbor = static_cast<UnstructGrid *>(grids[neighbor_zone]);
            int i_nIPoint_neightbor = gridpart_neighbor->GetInterpointInfo()->GetNumberOfInterpoints();
            int *interPoint2InterPointID_neightbor = gridpart_neighbor->GetInterpointInfo()->GetInterPoint2InterPointID();
            int *interPoint2ZoneID_neightbor = gridpart_neighbor->GetInterpointInfo()->GetInterPoint2ZoneID();
            
            for (int jPoint = 0; jPoint <i_nIPoint_neightbor; jPoint ++)
            {
                if (interPoint2InterPointID[iPoint] == interPoint2InterPointID_neightbor[jPoint] &&
                    interPoint2ZoneID_neightbor[jPoint] == iproc)
                {
                    interPoint2InterPointID[iPoint] = jPoint;
                    interPoint2InterPointID_neightbor[jPoint] = iPoint;
                    break;
                }
            }
        }
    }
}

void Pre_UnstructGridPartition::SetInterfaceInfoPartForParallelPartition(int izone, UnstructGrid *gridpart)
{
    if (npartmethod != 1) return;

    InterfaceInfo *interfaceInfo = gridpart->GetInterfaceInfo();

    int i_nIFace = interfaceInfo->GetNIFace();
    int i_nBFace = gridpart->GetNBoundFace();

    int nTotalFace = gridSource->GetNTotalFace();
    int nBoundFace = gridSource->GetNBoundFace();
    int *left_cell_of_face  = gridSource->GetLeftCellOfFace();
    int *right_cell_of_face = gridSource->GetRightCellOfFace();

    int le, re;

    int *interFace2ZoneID       = interfaceInfo->GetInterFace2ZoneID();
    int *interFace2InterFaceID  = interfaceInfo->GetInterFace2InterFaceID();
    int *interFace2BoundaryFace = interfaceInfo->GetInterFace2BoundaryFace();
    int *interFaceDirection     = interfaceInfo->GetInterFaceDirection();

    int nphysi = i_nBFace - i_nIFace;

    //! Set default value for parallel partition.
    for (int iface = 0; iface < i_nIFace; ++ iface)
    {
        interFace2BoundaryFace[iface] = -1;
        interFace2InterFaceID[iface]  = -1;
        interFace2ZoneID[iface]       = -1;
    }

    //! Here,the interFace2InterFaceID temporarily stores the cells index of neighbor zone.

    int count = 0;
    //! Interfaces in the original interior field.
    for (int iface = nBoundFace; iface < nTotalFace; ++ iface)
    {
        //! face_mark: the corresponding relationship from global face index to local face index.
        //! The following judgment condition means that this face is the interface of the zone.
        //! the index of a face is greater than 0 and less than i_nBFace in new partition zone, and it is the internal face of the original grid.
        //! So it must be interface.
        if (face_mark[iface] < i_nBFace && face_mark[iface] > -1)
        {
            ++ count;
            le = left_cell_of_face [iface];
            re = right_cell_of_face[iface];

            //! ii_local: the index of interface counts from 0.
            int ii_local = face_mark[iface] - nphysi;

            //! interface_interfaceid£ºtemporarily records local index of the neighbor interface cell.
            //! It will be modified to the interface index of neighbor zone.
            if (part[le] == izone)
            {
                interFaceDirection    [ii_local] = 1;
                interFace2ZoneID      [ii_local] = static_cast<int>(part[re]);
                interFace2InterFaceID [ii_local] = cell_mark[re];
                interFace2BoundaryFace[ii_local] = face_mark[iface];
            }
            else if (part[re] == izone)
            {
                interFaceDirection    [ii_local] = -1;
                interFace2ZoneID      [ii_local] = static_cast<int>(part[le]);
                interFace2InterFaceID [ii_local] = cell_mark[le];
                interFace2BoundaryFace[ii_local] = face_mark[iface];
            }
            else
            {
                TK_Exit::ExceptionExit("Error: CAN NOT find neighbor interface cell!");
            }
        }
    }

    //WriteLogFile(" nIFace in Interior = ", count);
}

void Pre_UnstructGridPartition::ComputeNodeMark(int izone, UnstructGrid *gridpart)
{
    int nTotalFace = gridSource->GetNTotalFace();
    int nTotalNode = gridSource->GetNTotalNode();
    int *node_number_of_each_face  = gridSource->GetNodeNumberOfEachFace();
    int *face2node = gridSource->GetFace2Node();

    for (int inode = 0; inode < nTotalNode; ++ inode)
    {
        node_mark[inode] = NOT_IN_THIS_ZONE;
    }

    //! Set all the node_mark belong to izone to -1.
    int count = 0;
    for (int iface = 0; iface < nTotalFace; ++ iface)
    {
        if (face_mark[iface] > -1)
        {
            for (int j = 0; j < node_number_of_each_face[iface]; ++ j)
            {
                node_mark[face2node[count ++]] = IN_THIS_ZONE;
            }
        }
        else
        {
            count += node_number_of_each_face[iface];
        }
    }

    int i_nTNode = 0;
    for (int inode = 0; inode < nTotalNode; ++ inode)
    {
        if (IsBelongtoThisZone(node_mark[inode]))
        {
            node_mark[inode] = i_nTNode ++;
        }
    }

    gridpart->SetNTotalNode(i_nTNode);
}

void Pre_UnstructGridPartition::ComputeFaceMark(int izone, UnstructGrid *gridpart)
{
    int nTotalCell = gridSource->GetNTotalCell();
    int nTotalFace = gridSource->GetNTotalFace();
    int nBoundFace = gridSource->GetNBoundFace();
    int *left_cell_of_face  = gridSource->GetLeftCellOfFace();
    int *right_cell_of_face = gridSource->GetRightCellOfFace();
    CompositeUnstructBC **bcr = gridSource->GetBCRecord();
    int le, re;

    for (int iface = 0; iface < nTotalFace; ++ iface)
    {
        face_mark[iface] = NOT_IN_THIS_ZONE;
    }

    //! The face belong to this zone is set to be 'IN_THIS_ZONE'.
    for (int iface = 0; iface < nTotalFace; ++ iface)
    {
        le = left_cell_of_face [iface];
        re = right_cell_of_face[iface];

        if (part[le] == izone)
        {
            face_mark[iface] = IN_THIS_ZONE;
        }
        else if (re < nTotalCell && part[re] == izone)
        {
            face_mark[iface] = IN_THIS_ZONE;
        }
    }

    int i_nTFace = 0;
    int i_nIFace = 0;

    //! BC faces.
    for (int iface = 0; iface < nBoundFace; ++ iface)
    {
        if (IsBelongtoThisZone(face_mark[iface]))
        {
            face_mark[iface] = i_nTFace ++;

            int bcType = bcr[iface]->GetKey();
            if (bcType == PHENGLEI::INTERFACE)
            {
                //! Parallel partition.
                interfaceIndexInBC_send[iface] = i_nIFace;
                i_nIFace ++;
            }
        }
    }

    nIFaceInOrgBC[izone] = i_nIFace;

    //! Interior interfaces.
    for (int iface = nBoundFace; iface < nTotalFace; ++ iface)
    {
        if (IsBelongtoThisZone(face_mark[iface]))
        {
            le = left_cell_of_face [iface];
            re = right_cell_of_face[iface];
            if (part[le] != part[re])
            {
                face_mark[iface] = i_nTFace ++;
                i_nIFace ++;
            }
        }
    }

    int i_nBFace = i_nTFace;

    //! Interior NON-interface.
    for (int iface = nBoundFace; iface < nTotalFace; ++ iface)
    {
        if (face_mark[iface] < 0)
        {
            //! Has not been set.
            if (IsBelongtoThisZone(face_mark[iface]))
            {
                le = left_cell_of_face [iface];
                re = right_cell_of_face[iface];

                if (part[le] == part[re])
                {
                    face_mark[iface] = i_nTFace ++;
                }
            }
        }
    }

    gridpart->SetNTotalFace(i_nTFace);
    gridpart->SetNBoundFace(i_nBFace);
    if (i_nIFace > 0)
    {
        InterfaceInfo *interfaceInfo = new InterfaceInfo(i_nIFace);
        gridpart->SetInterfaceInfo(interfaceInfo);

        //! Only partition procedure need the data of isNeighborInterfaceFound,
        //! so new it here.
        bool *isNeighborInterfaceFound = new bool [i_nIFace];
        PHSPACE::SetField(isNeighborInterfaceFound, false, i_nIFace);
        interfaceInfo->SetIsNeighborInterfaceFound(isNeighborInterfaceFound);
    }
}

int Pre_UnstructGridPartition::GetNTotalCell(int izone, idx_t *part)
{
    int nTotalCell = gridSource->GetNTotalCell();
    int count = 0;
    for (int icell = 0; icell < nTotalCell; ++ icell)
    {
        if (part[icell] == izone)
        {
            count ++;
        }
    }
    return count;
}

int * Pre_UnstructGridPartition::GetNode2All(int izone, UnstructGrid *gridpart)
{
    int i_nTNode = gridpart->GetNTotalNode();

    int *node2all = new int[i_nTNode];    //! node2all, record the connection relationship from nodes in partition zone to all nodes.

    int nTotalNode = gridSource->GetNTotalNode();

    for (int inode = 0; inode < nTotalNode; ++ inode)
    {
        if (node_mark[inode] > -1)
        {
            node2all[node_mark[inode]] = inode;
        }
    }

    return node2all;
}

int * Pre_UnstructGridPartition::GetFace2All(int izone, UnstructGrid *gridpart)
{
    int i_nTFace = gridpart->GetNTotalFace();

    int *face2all = new int[i_nTFace];    //! face2all,record the connection relationship from faces in partition zone to all faces(after ordering).

    int nTotalFace = gridSource->GetNTotalFace();

    for (int iface = 0; iface < nTotalFace; ++ iface)
    {
        if (face_mark[iface] > -1)
        {
            face2all[face_mark[iface]] = iface;
        }
    }

    return face2all;
}

int * Pre_UnstructGridPartition::GetCell2All(int izone, UnstructGrid *gridpart)
{
    int nTotalCell = gridpart->GetNTotalCell();

    int *cell2all = new int[nTotalCell];    //! cell2all,record the connection relationship from cells in partition zone to all cells.

    int nTotalCell_All = gridSource->GetNTotalCell();

    int count = 0;
    for (int icell = 0; icell < nTotalCell_All; ++ icell)
    {
        if (part[icell] == izone)
        {
            cell2all[count] = icell;
            ++ count;
        }
    }

    return cell2all;
}

void Pre_UnstructGridPartition::ComputePoint2ZoneID()
{
    int *n2c = gridSource->GetNode2Cell();
    int nTotalCell = gridSource->GetNTotalCell();
    int *nCPN = gridSource->GetCellNumberOfEachNode();
    int nTotalNode = gridSource->GetNTotalNode();

    set<int> *point2ZoneIDSet = new set<int>[nTotalNode];
    vector <int> pointsForPartition;
    vector <int> point2ZoneID_Vector;

    pointPosition2zoneID = new int [nTotalNode + 1];
    pointPosition2zoneID[0] = 0;

    int count = 0;
    for (int inode = 0; inode < nTotalNode; ++ inode)
    {
        //! node2ZoneSet is used to store and filter the repeated corresponding relationship from nodes to zone.
        for (int jcell=0; jcell < nCPN[inode]; ++ jcell)
        {
            if (!(n2c[count + jcell] < nTotalCell))
            {
                continue;
            }
                point2ZoneIDSet[inode].insert(static_cast<int>(part[n2c[count + jcell]]));
        }
        count += nCPN[inode];

        //! The size of the corresponding relationship from nodes to zone can not be determined,we can not creat static array to store.
        //! So it can only be stored with vector.
        for (set<int>::iterator iter = point2ZoneIDSet[inode].begin(); iter != point2ZoneIDSet[inode].end(); ++ iter)
        {
            pointsForPartition.push_back(inode);
            point2ZoneID_Vector.push_back((*iter));
        }
        pointPosition2zoneID[inode+1] = pointPosition2zoneID[inode] + static_cast<int>(point2ZoneIDSet[inode].size());
    }

    int pointsSize = static_cast<int>(pointsForPartition.size());
    pointNumber2zoneID = new int [pointsSize];
    point2ZoneID = new int [pointsSize];
    
    for (int iPoint =0; iPoint < pointsSize; iPoint ++)
    {
        pointNumber2zoneID[iPoint] = pointsForPartition[iPoint];
        point2ZoneID[iPoint] = point2ZoneID_Vector[iPoint];
    }

    delete [] point2ZoneIDSet;
}

void Pre_UnstructGridPartition::ComputeCellNumberForNodeValue()
{
    int nTotalNode = gridSource->GetNTotalNode();
    cellNumberForNodeValue = new int[nTotalNode];

    int nBoundFace = gridSource->GetNBoundFace();
    int nTotalFace = gridSource->GetNTotalFace();

    int *face2node = gridSource->GetFace2Node();
    int *node_number_of_each_face = gridSource->GetNodeNumberOfEachFace();

    for (int iNode = 0; iNode < nTotalNode; iNode ++)
    {
        cellNumberForNodeValue[iNode] = 0;
    }

    int pt;
    int nodepos = 0;
    for (int iface = 0; iface < nTotalFace; ++ iface)
    {
        for (int j = 0; j < node_number_of_each_face[iface]; ++ j)
        {
            pt = face2node[nodepos + j];
            //! From left.
            cellNumberForNodeValue[pt] += 1;
            //! From right.
            cellNumberForNodeValue[pt] += 1;
        }
        nodepos += node_number_of_each_face[iface];
    }

    nodepos = 0;
    int bctype;
    CompositeUnstructBC **bcr = gridSource->GetBCRecord();

    for (int iface = 0; iface < nBoundFace; ++ iface)
    {
        bctype = bcr[iface]->GetKey();
        if (bctype == PHENGLEI::SYMMETRY)
        {
            for (int j = 0; j < node_number_of_each_face[iface]; ++ j)
            {
                pt = face2node[nodepos + j];
                cellNumberForNodeValue[pt] = 0;
            }
        }
        nodepos += node_number_of_each_face[iface];
    }

    nodepos = 0;
    for (int iface = 0; iface < nBoundFace; ++ iface)
    {
        bctype = bcr[iface]->GetKey();
        if (bctype == PHENGLEI::SYMMETRY)
        {
            for (int j = 0; j < node_number_of_each_face[iface]; ++ j)
            {
                pt = face2node[nodepos + j];
                //! From left.
                cellNumberForNodeValue[pt] += 1;
                //! From right.
                cellNumberForNodeValue[pt] += 1;
            }
        }
        nodepos += node_number_of_each_face[iface];
    }

    nodepos = 0;
    for (int iface = 0; iface < nBoundFace; ++ iface)
    {
        bctype = bcr[iface]->GetKey();
        if (bctype == PHENGLEI::FARFIELD)
        {
            for (int j = 0; j < node_number_of_each_face[iface]; ++ j)
            {
                pt = face2node[nodepos + j];

                cellNumberForNodeValue[pt] = 0;
            }
        }
        nodepos += node_number_of_each_face[iface];
    }

    nodepos = 0;
    for (int iface = 0; iface < nBoundFace; ++ iface)
    {
        bctype = bcr[iface]->GetKey();
        if (bctype == PHENGLEI::FARFIELD)
        {
            for (int j = 0; j < node_number_of_each_face[iface]; ++ j)
            {
                pt = face2node[nodepos + j];
                //! From left.
                cellNumberForNodeValue[pt] += 1;
                //! From right.
                cellNumberForNodeValue[pt] += 1;
            }
        }
        nodepos += node_number_of_each_face[iface];
    }

    nodepos = 0;
    for (int iface = 0; iface < nBoundFace; ++ iface)
    {
        bctype = bcr[iface]->GetKey();
        if (bctype == PHENGLEI::SOLID_SURFACE)
        {
            for (int j = 0; j < node_number_of_each_face[iface]; ++ j)
            {
                pt = face2node[nodepos + j];

                cellNumberForNodeValue[pt] = 0;
            }
        }
        nodepos += node_number_of_each_face[iface];
    }

    nodepos = 0;
    for (int iface = 0; iface < nBoundFace; ++ iface)
    {
        bctype = bcr[iface]->GetKey();
        if (bctype == PHENGLEI::SOLID_SURFACE)
        {
            for (int j = 0; j < node_number_of_each_face[iface]; ++ j)
            {
                pt = face2node[nodepos + j];
                //! From left.
                cellNumberForNodeValue[pt] += 1;
                //! From right.
                cellNumberForNodeValue[pt] += 1;
            }
        }
        nodepos += node_number_of_each_face[iface];
    }
}

void Pre_UnstructGridPartition::SetNTotalCell(int izone, UnstructGrid *gridpart)
{
    int i_nTotalCell = GetNTotalCell(izone, part);
    gridpart->SetNTotalCell(i_nTotalCell);
}

void Pre_UnstructGridPartition::CreateLocal2All(int izone, UnstructGrid *gridpart)
{
    node2all = GetNode2All(izone, gridpart);
    face2all = GetFace2All(izone, gridpart);
    cell2all = GetCell2All(izone, gridpart);
}

void Pre_UnstructGridPartition::ComputeGridForPartition(UnstructGrid *gridpart, int izone)
{
    SetNTotalCell(izone, gridpart);
    ComputeFaceMark(izone, gridpart);
    ComputeNodeMark(izone, gridpart);
    
    CreateLocal2All(izone, gridpart);
    
    SetCoorPart(gridpart);
    SetLinkPart(izone, gridpart);
    SetInterfaceInfoPart(izone, gridpart);
    SetInterpointInfoPart(izone, gridpart);
    
    SetWallDist(gridpart);
    
    FreeLocal2All();
}

void Pre_UnstructGridPartition::SetWallDist(UnstructGrid *gridpart)
{
    if (!partitionWallDist) return;
    
    int nTotalCell = gridpart->GetNTotalCell();
    gridpart->InitWallDist();
    
    RDouble *walldist = gridpart->GetWallDist();
    RDouble *walldist_src = gridSource->GetWallDist();

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        int iCell2ALL = cell2all[iCell];
        walldist[iCell] = walldist_src[iCell2ALL];
    }
}

void Pre_UnstructGridPartition::WriteAdditionalData(const string &out_grid_file)
{
    WriteInterpointInfo(out_grid_file);

    //! Write cell to node connection.
    WriteCellToNode(out_grid_file);

    //! Write face boundary condition information.
    WriteFaceBC(out_grid_file);

    //! Write walldist information.
    WriteWalldist(out_grid_file);
}

void Pre_UnstructGridPartition::WriteWalldist(const string &out_grid_file)
{
    if (!partitionWallDist) return;
    string walldistfile = ChangeExtensionOfFileName(out_grid_file, "wdt");

    ActionKey *actkey = new ActionKey();
    actkey->filename = (walldistfile);

    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    string finalFileName = PHSPACE::AddSymbolToFileName(walldistfile, '_', currentProcessorIndex);
    actkey->filepos = H5Fcreate(finalFileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    using namespace PHMPI;
    for (int izone = 0; izone < maxProcessor; ++ izone)
    {
        if (!grids[izone]) continue;

        UnstructGrid *grid = UnstructGridCast(grids[izone]);
        WriteWalldist(actkey, grid);
    }
    
    H5Fclose(actkey->filepos);
    delete actkey;    actkey = NULL;
}

void Pre_UnstructGridPartition::WriteWalldist(ActionKey *actkey, UnstructGrid *gridPart)
{
    ostringstream oss;
    oss << "Grid-" << gridPart->GetZoneID();
    string dataName = oss.str();

    int nTotalCell = gridPart->GetNTotalCell();
    RDouble *walldist = gridPart->GetWallDist();

    CreateAndWriteData(actkey->filepos, dataName, 1, nTotalCell, PHDOUBLE, walldist);
}

void Pre_UnstructGridPartition::WriteCellToNode(const string &out_grid_file)
{
    if (!partitionCellToNode) return;
    string partitionedCellToNodeFileName = ChangeExtensionOfFileName(out_grid_file, "c2n");
    fstream file;
    PHSPACE::OpenFile(file, partitionedCellToNodeFileName, ios_base::out|ios_base::binary);

    PHWrite(file, maxProcessor);
    for (int iproc = 0; iproc < maxProcessor; ++ iproc)
    {
        if (!grids[iproc]) continue;

        UnstructGrid *gridPart = UnstructGridCast(grids[iproc]);
        WriteCellToNode(file, gridPart);
    }

    PHSPACE::CloseFile(file);
}

void Pre_UnstructGridPartition::WriteCellToNode(fstream &file, UnstructGrid *gridPart)
{
    VirtualFile *virtualFile = new VirtualFile(&file);
    virtualFile->BeginWriteWork();

    int *cellNodeIndexContainer  = gridPart->GetCell2Node();
    int *cellNodeNumberContainer = gridPart->GetNodeNumberOfEachCell();
    int iZone = gridPart->GetGridID()->GetIndex();

    int numberOfCells = gridPart->GetNTotalCell();

    PHWrite(virtualFile, iZone);
    PHWrite(virtualFile, numberOfCells);

    int iPosition = 0;
    for (int iCell = 0; iCell < numberOfCells; ++ iCell)
    {
        int nodeNumber = cellNodeNumberContainer[iCell];

        PHWrite(virtualFile, nodeNumber);

        for (int iNode = 0; iNode < nodeNumber; ++ iNode)
        {
            PHWrite(virtualFile, cellNodeIndexContainer[iPosition ++]);
        }
    }

    virtualFile->EndWriteWork();

    delete virtualFile;
}

void Pre_UnstructGridPartition::WriteInterpointInfo(const string &out_grid_file)
{
    if (IsParallelPartition() && (parallelPartitionMethod == PARMETIS_METHOD))
    {
        return;
    }

    string interpointInfoFileName = ChangeExtensionOfFileName(out_grid_file, "interpoint");
    fstream file;
    PHSPACE::OpenFile(file, interpointInfoFileName, ios_base::out|ios_base::binary);

    PHWrite(file, maxProcessor);
    for (int iproc = 0; iproc < maxProcessor; ++ iproc)
    {
        if (!grids[iproc])
        {
            continue;
        }
        UnstructGrid *gridPart = UnstructGridCast(grids[iproc]);
        WriteInterpointInfo(file, gridPart);
    }

    PHSPACE::CloseFile(file);
}

void Pre_UnstructGridPartition::WriteInterpointInfo(fstream &file, UnstructGrid *gridPart)
{
    VirtualFile *virtualFile = new VirtualFile(&file);
    virtualFile->BeginWriteWork();
    
    int iZone = gridPart->GetGridID()->GetIndex();
    
    InterpointInformation *interpointInformation = gridPart->GetInterpointInfo();
    
    int i_nIPoint = interpointInformation->GetNumberOfInterpoints();
    int *interPoint2ZoneID = interpointInformation->GetInterPoint2ZoneID();
    int *interPoint2InterPointID = interpointInformation->GetInterPoint2InterPointID();
    int *interPoint2GlobalPoint = interpointInformation->GetInterPoint2GlobalPoint();
    int *cellNumberOfInterPoint = interpointInformation->GetCellNumberOfInterPoint();
    int *totalZonesOfInterPoint = interpointInformation->GetTotalZonesOfInterPoint();
    int *labelOfInterPoint = interpointInformation->GetLabelOfInterPoint();

    PHWrite(virtualFile, iZone);
    PHWrite(virtualFile, i_nIPoint);
    
    for (int iPoint = 0; iPoint < i_nIPoint; iPoint ++)
    {
        PHWrite(virtualFile, interPoint2ZoneID[iPoint]);
        PHWrite(virtualFile, interPoint2InterPointID[iPoint]);
        PHWrite(virtualFile, interPoint2GlobalPoint[iPoint]);
        PHWrite(virtualFile, cellNumberOfInterPoint[iPoint]);
        PHWrite(virtualFile, totalZonesOfInterPoint[iPoint]);
        PHWrite(virtualFile, labelOfInterPoint[iPoint]);
    }

    virtualFile->EndWriteWork();

    delete virtualFile;
}

void Pre_UnstructGridPartition::WriteFaceBC(const string &out_grid_file)
{
    if (!partitionFaceBC) return;
    string partitionedFaceBCFileName = ChangeExtensionOfFileName(out_grid_file, "bc");
    fstream file;
    PHSPACE::OpenFile(file, partitionedFaceBCFileName, ios_base::out|ios_base::binary);

    PHWrite(file, maxProcessor);
    for (int iproc = 0; iproc < maxProcessor; ++ iproc)
    {
        if (!grids[iproc]) continue;

        UnstructGrid *gridPart = UnstructGridCast(grids[iproc]);
        WriteFaceBC(file, gridPart);
    }

    PHSPACE::CloseFile(file);
}

void Pre_UnstructGridPartition::WriteFaceBC(fstream &file, UnstructGrid *gridPart)
{
    gridPart->WriteFaceBC(file);
}

void Pre_UnstructGridPartition::ComputeGlobalToLocalCellToNodeMapping(UnstructGrid *gridpart, int izone)
{
    if (!partitionCellToNode) return;

    int *cellToNodeGlobal       = gridSource->GetCell2Node();
    int *nodeNumberOfCellGlobal = gridSource->GetNodeNumberOfEachCell();

    int numberOfCells = gridpart->GetNTotalCell();
    
    int *nodeNumberOfCellLocal = new int [numberOfCells];
    gridpart->SetNodeNumberOfEachCell(nodeNumberOfCellLocal);

    int numberOfGlobalCells = gridSource->GetNTotalCell();
    int *globalToLocalNodeIndexMapping = this->node_mark;
    idx_t *globalCellIndexToLocalZoneIndexMapping = this->part;

    int iCount = 0;
    int count = 0;
    for (int iCell = 0; iCell < numberOfGlobalCells; ++ iCell)
    {
        if (globalCellIndexToLocalZoneIndexMapping[iCell] == izone)
        {
            int nodeNumber = nodeNumberOfCellGlobal[iCell];
            nodeNumberOfCellLocal[iCount ++] = nodeNumber;
            count += nodeNumber;
        }
    }

    int *cellToNodeLocal = new int [count];
    gridpart->SetCell2Node(cellToNodeLocal);

    int iPositionGlobal = 0;
    count = 0;
    for (int iCell = 0; iCell < numberOfGlobalCells; ++ iCell)
    {
        int nodeNumber = nodeNumberOfCellGlobal[iCell];

        if (globalCellIndexToLocalZoneIndexMapping[iCell] == izone)
        {
            for (int iNode = 0; iNode < nodeNumber; ++ iNode)
            {
                int globalNodeIndex = cellToNodeGlobal[iPositionGlobal + iNode];
                int nodeIndex = globalToLocalNodeIndexMapping[globalNodeIndex];
                cellToNodeLocal[count ++] = nodeIndex;
            }
        }
        iPositionGlobal += nodeNumber;
    }
    gridpart->SetCell2Node(cellToNodeLocal);
    gridpart->SetNodeNumberOfEachCell(nodeNumberOfCellLocal);
}

bool Pre_UnstructGridPartition::IsParallelPartition()
{
    //if(PHMPI::GetNumberofGlobalZones() > 1)
    if (vtxdist)
    {
        return true;
    }
    else
    {
        return false;
    }
}

void Pre_UnstructGridPartition::FillnTCellofEachZone(int *&nTCellOfZones)
{
    WriteLogFile("Start to Fill nTotalCell of Each Zone");

    int nZones = PHMPI::GetNumberofGlobalZones();
    int myid   = PHMPI::GetCurrentProcessorID();

    int nTProcessor = PHMPI::GetNumberOfProcessor();

    if (nTProcessor == nZones)
    {
        //! Each processors has ONE zone.
        int nTCellInThisZone = gridSource->GetNTotalCell();

        int *localSizeTemp = new int [nTProcessor];
        SetField(localSizeTemp, 1, nTProcessor);
        int **nTCellListTemp = new int * [nTProcessor];

        PHMPI::FillGobalDataCollectBcast(&nTCellInThisZone, localSizeTemp, nTCellListTemp);

        for (int iProc = 0; iProc < nTProcessor; ++ iProc)
        {
            nTCellOfZones[iProc] = nTCellListTemp[iProc][0];
        }

        DelPointer2(nTCellListTemp);
        DelPointer(localSizeTemp);
    }
    else
    {
        //! Each processors has multiple zone.
        for (int iZone = 0; iZone < nZones; ++ iZone)
        {
            int procID = PHMPI::GetZoneProcessorID(iZone);
            int nTCellInThisZone;

            if (myid == procID)
            {
                nTCellInThisZone = GetGrid(iZone, 0)->GetNTotalCell();
            }

            PHMPI::PH_Bcast(&nTCellInThisZone, sizeof(int), procID);

            nTCellOfZones[iZone] = nTCellInThisZone;
        }
    }

    WriteLogFile("End to Fill nTCell of Each Zone");
}

LIB_EXPORT void Pre_UnstructGridPartition::DumpPartitionGrids(const string &out_grid_file)
{
    if (IsParallelPartition())
    {
        DumpParallelPartitionGrid(out_grid_file, grids, maxProcessor);
    }
    else
    {
        DumpSerialPartitionGrid(out_grid_file, grids, maxProcessor);
    }
}

void Pre_UnstructGridPartition::DumpSerialPartitionGrid(const string &gridFileName, Grid **grids, int nBlocks)
{
    DumpGrid(gridFileName, grids, maxProcessor);

    string partitionedGridFileName = AddSymbolToFileName(gridFileName, '_', PHMPI::GetCurrentProcessorID());

    int dumpOldGrid = GlobalDataBase::GetIntParaFromDB("dumpOldGrid");
    if (dumpOldGrid)
    {
        WriteAdditionalData(partitionedGridFileName);
    }

    WriteWalldist(gridFileName);
}

void Pre_UnstructGridPartition::DumpParallelPartitionGrid(const string &gridFileName, Grid **grids, int nBlocks)
{
    int dumpOldGrid = GlobalDataBase::GetIntParaFromDB("dumpOldGrid");
    if (!dumpOldGrid)
    {
        DumpGrid(gridFileName, grids, nBlocks);

        WriteWalldist(gridFileName);
        return;
    }

    //! Compute the actual number of parts in this processor.
    int nBlocksActual = 0;
    for (int iBlock = 0; iBlock < nBlocks; ++ iBlock)
    {
        if (grids[iBlock])
        {
            ++ nBlocksActual;
        }
    }
    
    Grid **gridsActual = new Grid * [nBlocksActual];
    nBlocksActual = 0;
    for (int iBlock = 0; iBlock < nBlocks; ++ iBlock)
    {
        if (grids[iBlock])
        {
            gridsActual[nBlocksActual] = grids[iBlock];
            ++ nBlocksActual;
        }
    }
    
    int nZones = PHMPI::GetNumberofGlobalZones();
    int nProc  = PHMPI::GetNumberOfProcessor();
    ASSERT(nProc == nZones);
    int myid = PHMPI::GetCurrentProcessorID();
    int *nPartOfEachProcessor = new int[nZones];

    if (nProc == nZones)
    {
        //! Each processors has ONE zone.
        int *localSizeTemp = new int [nProc];
        SetField(localSizeTemp, 1, nProc);
        int **nPartListTemp = new int * [nProc];

        PHMPI::FillGobalDataCollectBcast(&nBlocksActual, localSizeTemp, nPartListTemp);

        for (int iProc = 0; iProc < nProc; ++ iProc)
        {
            nPartOfEachProcessor[iProc] = nPartListTemp[iProc][0];
        }

        DelPointer2(nPartListTemp);
        DelPointer(localSizeTemp);
    }
    else
    {
        PHMPI::FillGobalData(nBlocksActual, nPartOfEachProcessor);
    }

    int *block_index = new int[nBlocksActual];
    int count = 0;
    for (int iProc = 0; iProc < myid; ++ iProc)
    {
        count += nPartOfEachProcessor[iProc];
    }
    for (int iBlock = 0; iBlock < nBlocksActual; ++ iBlock)
    {
        block_index[iBlock] = iBlock + count;
    }

    string partitionedGridFileName = AddSymbolToFileName(gridFileName, '_', PHMPI::GetCurrentProcessorID());
    DumpParallelPartitionGrid(partitionedGridFileName, gridsActual, nBlocksActual, block_index);
    WriteAdditionalData(partitionedGridFileName);

    delete [] gridsActual;
    delete [] nPartOfEachProcessor;
    delete [] block_index;
}

void Pre_UnstructGridPartition::DumpParallelPartitionGrid(const string &gridFileName, Grid **grids, int nBlocks, int *block_index_in)
{
    if (gridFileName == "") return;
    fstream file;
    PHSPACE::OpenFile(file, gridFileName, ios_base::out|ios_base::binary|ios_base::trunc);

    int *block_proc     = new int [nBlocks];
    int *block_idx      = new int [nBlocks];
    int *block_type     = new int [nBlocks];
    int *block_proc_inp = new int [nBlocks];

    for (int izone = 0; izone < nBlocks; ++ izone)
    {
        block_proc[izone] = grids[izone]->GetGridID()->GetIndex();  // Processor ID actually.
        block_idx [izone] = block_index_in[izone];
        block_type[izone] = grids[izone]->Type();
    }

    if (parallelPartitionMethod == METIS_OPENMP_METHOD)
    {
        //! For OpenMP case, reset the processor ID of sub-zone to the source grid's.
        int zoneID = gridSource->GetGridID()->GetIndex();
        int processorID = PHMPI::GetZoneProcessorID(zoneID);
        for (int izone = 0; izone < nBlocks; ++ izone)
        {
            block_proc[izone] = processorID;
        }
    }

    file.write(reinterpret_cast<char*>(&nBlocks), sizeof(int));
    file.write(reinterpret_cast<char*>(block_proc), nBlocks * sizeof(int));
    file.write(reinterpret_cast<char*>(block_idx), nBlocks * sizeof(int));
    file.write(reinterpret_cast<char*>(block_type), nBlocks * sizeof(int));
    
    for (int izone = 0; izone < nBlocks; ++ izone)
    {
        grids[izone]->WriteGrid(file);
    }

    file.close();
    file.clear();

    delete [] block_idx;
    delete [] block_proc;
    delete [] block_proc_inp;
    delete [] block_type;
}

}
