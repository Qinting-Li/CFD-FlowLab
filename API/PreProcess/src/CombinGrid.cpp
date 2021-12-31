#include "CombinGrid.h"
#include "Math_BasisFunction.h"
#include "Glb_Dimension.h"
#include "IO_FileName.h"
#include "Region.h"
#include "OverLappingGrid.h"
#include "Pre_GridConversion.h"
#include "PHIO.h"
#include "TK_Log.h"

using namespace std;

namespace PHSPACE
{

CombinGrid::CombinGrid(Region *region, Grid **grid_in, int nZones_in)
{
    this->region = region;
    this->grid_in = grid_in;
    nZones = nZones_in;

    nTotalCell = 0;
    nTotalFace = 0;
    nTotalNode = 0;
    nIFace = 0;
    nTotalNodeEachZone = new int [nZones];
    nTotalFaceEachZone = new int [nZones];
    nTotalCellEachZone = new int [nZones];
    nIFaces = new int [nZones];

    face2NewFaceafterReorder = 0;
    face2Interface = 0;
    face2InterfaceafterReorder = 0;
    interface2BoundaryFaceafterReorder = 0;
    isGlobalInterFace = 0;
}

CombinGrid::~CombinGrid()
{
    delete [] nTotalNodeEachZone;
    delete [] nTotalFaceEachZone;
    delete [] nTotalCellEachZone;
    delete [] nIFaces;
    delete [] left_cell_of_all_face;
    delete [] right_cell_of_all_face;
    delete [] node_number_of_all_face;
    delete [] isFaceMerged;
    delete [] face2Interface;
    delete [] isGlobalInterFace;
    delete [] face2InterfaceafterReorder;
    delete [] interface2BoundaryFaceafterReorder;
    delete [] face2NewFaceafterReorder;

    DelPointer2(node2all);
    DelPointer2(face2all);
    DelPointer2(cell2all);
}

void CombinGrid::CheckifIndependentZoneExist()
{
    int myid = PHMPI::GetCurrentProcessorID();

    int numberofGlobalZones = PHMPI::GetNumberofGlobalZones();
    int localIndex = 0;
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int gridIndex = grid_in[iZone]->GetZoneID();
        zoneIDinCurrentProc.insert(gridIndex);
        zoneIndexMap.insert(pair<int, int>(gridIndex ,localIndex));
        localIndex ++;
    }

    if (nZones != 1)
    {
        for (int iZone = 0; iZone < numberofGlobalZones; ++ iZone)
        {
            int processorID = PHMPI::GetZoneProcessorID(iZone);

            if (processorID == myid)
            {
                bool isNeighborinThisProc = false;

                ZoneNeighbor *globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
                std::size_t numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();
                
                for (std::size_t iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
                {
                    int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);

                    set<int>::iterator iter;
                    for (iter = zoneIDinCurrentProc.begin(); iter != zoneIDinCurrentProc.end(); ++ iter)
                    {
                        int zoneID = *iter;
                        if (neighborZone == zoneID)
                        {
                            isNeighborinThisProc = true;
                            break;
                        }
                    }
                }

                if (!isNeighborinThisProc)
                {
                    TK_Exit::ExceptionExit("Error: An single zone that couldn't merged to any others has been found!");
                }
            }
        }
    }

    WriteLogFile("All zones in this processor could be merged!");
    PrintToWindow("All zones in this processor could be merged!\n");
}

Grid * CombinGrid::CreateCombinationGrid()
{
    CheckifIndependentZoneExist();

    Init();

    //! Warning:
    //! All assumptions are that multi blocks in a process are merged into one block,
    //! If there is no common face between two blocks, a lot logic relationships need to be reconsidered.
    //! For example,MergeCell2All is the zone number of new grid.
    MergeCell2All();
    MergeNode2All();
    MergeFace2All();
    MergeInterface2All();

    int nTotalNode = GetNTotalNode();
    int nTotalFace = GetNTotalFace();
    int nTotalCell = GetNTotalCell();
    int nIFace = GetNIFace();

    ostringstream oss;
    oss << "After nodes faces cells merged ..." << "\n";
    oss << "  Number of Points     : " << nTotalNode << "\n";
    oss << "  Number of Faces      : " << nTotalFace << "\n";
    oss << "  Number of Cells      : " << nTotalCell << "\n";
    oss << "  Number of InterFaces : " << nIFace << "\n";
    PrintToWindow(oss);
    WriteLogFile(oss);

    int myid = PHMPI::GetCurrentProcessorID();
    Grid *grid_all = CreateGridGeneral(UNSTRUCTGRID, new GridID(myid), 0, GetDim());

    grid_all->SetNTotalNode(nTotalNode);
    grid_all->SetNTotalFace(nTotalFace);
    grid_all->SetNTotalCell(nTotalCell);

    return grid_all;
}

void CombinGrid::Init()
{
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        nTotalNodeEachZone[iZone] = grid_in[iZone]->GetNTotalNode();
        nTotalFaceEachZone[iZone] = grid_in[iZone]->GetNTotalFace();
        nTotalCellEachZone[iZone] = grid_in[iZone]->GetNTotalCell();
        nIFaces[iZone] = grid_in[iZone]->GetNIFace();

        nTotalCell += nTotalCellEachZone[iZone];
        nTotalFace += nTotalFaceEachZone[iZone];
        nIFace     += nIFaces[iZone];
        nTotalNode += nTotalNodeEachZone[iZone];
    }

    InitLocal2All();
}

void CombinGrid::InitLocal2All()
{
    WriteLogFile("Init local to all ...");
    InitNode2All();
    InitFace2All();
    InitCell2All();
}

void CombinGrid::RunCombinGrid()
{
    grid_all = UnstructGridCast(CreateCombinationGrid());

    ResetCoordinate();
    ResetNodeNumberOfEachFace();
    ResetFace2Node();
    ResetLeftRightCell();
    MergeCell2Node();
    MergeBCType();
    SwapInterfaceData();

    PrintToWindow("Combine grids finished ...\n");
}

void CombinGrid::MergeInterface2All()
{
    WriteLogFile("Merge interface2all ...");
    PrintToWindow("Merge interface2all ...\n");

    int nTotalFace = GetNTotalFace();
    isGlobalInterFace = new bool [nTotalFace];
    PHSPACE::SetField(isGlobalInterFace, false, nTotalFace);

    face2Interface = new int [nTotalFace];
    PHSPACE::SetField(face2Interface, -1, nTotalFace);

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        UnstructGrid *gridLocal = UnstructGridCast(grid_in[iZone]);
        CompositeUnstructBC **bcrLocal = gridLocal->GetBCRecord();
        int nBFaceLocal = gridLocal->GetNBoundFace();
        int *face2all = GetFace2All(iZone);
        
        int nMergedInterfaces = 0;
        for (int iFace = 0; iFace < nBFaceLocal; ++ iFace)
        {
            int globalFaceID = face2all[iFace];

            if (bcrLocal[iFace]->GetKey() == PHENGLEI::INTERFACE)
            {
                if (!isFaceMerged[globalFaceID])
                {
                    isGlobalInterFace[globalFaceID] = true;
                }else
                {
                    ++ nMergedInterfaces;
                }
            }
        }

        WriteLogFile("number of merged interfaces: ", nMergedInterfaces);
        if (nZones > 1 && nMergedInterfaces == 0)
        {
            TK_Exit::ExceptionExit("The sub-zone is independent of other zones, consider more!");
        }
    }

    int nIFaceGlobal = 0;
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        if(isGlobalInterFace[iFace])
        {
            face2Interface[iFace] = nIFaceGlobal;
            ++ nIFaceGlobal;
        }
    }

    this->nIFace = nIFaceGlobal;
}

void CombinGrid::MergeFace2All()
{
    WriteLogFile("Merge face2all ...");
    PrintToWindow("Merge face2all ...\n");

    left_cell_of_all_face = new int[nTotalFace];
    right_cell_of_all_face = new int[nTotalFace];
    node_number_of_all_face = new int[nTotalFace];

    isFaceMerged = new bool [nTotalFace];
    PHSPACE::SetField(isFaceMerged, false, nTotalFace);

    nTotalFace = 0;

    int globalFace_count = 0;
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        UnstructGrid *unstr_current = UnstructGridCast(grid_in[iZone]);

        int nBFace_current                    = unstr_current->GetNBoundFace();
        int *node_number_of_each_face_current = unstr_current->GetNodeNumberOfEachFace();
        int *face2node_current                = unstr_current->GetFace2Node();
        int *left_cell_of_face_current        = unstr_current->GetLeftCellOfFace();
        int *right_cell_of_face_current       = unstr_current->GetRightCellOfFace();
        CompositeUnstructBC **bcr_current     = unstr_current->GetBCRecord();
        InterfaceInfo *iinfo_current          = unstr_current->GetInterfaceInfo();
        int *interFace2ZoneID_current         = iinfo_current->GetInterFace2ZoneID();
        int *interFace2InterFaceID_current    = iinfo_current->GetInterFace2InterFaceID();

        int *node2all = GetNode2All(iZone);
        int *face2all = GetFace2All(iZone);
        int *cell2all = GetCell2All(iZone);

        int interface_current = 0;
        int node_sum = 0;
        for (int iface = 0; iface < nTotalFaceEachZone[iZone]; ++ iface)
        {
            int node_number_of_iface = node_number_of_each_face_current[iface];
            if (iface < nBFace_current)
            {
                int bctype = bcr_current[iface]->GetKey();
                if (bctype == PHENGLEI::INTERFACE)
                {
                    int GridID_Neighbor = interFace2ZoneID_current[interface_current];

                    bool isNeighborinThisProc = false;
                    set<int>::iterator iter;
                    for (iter = zoneIDinCurrentProc.begin(); iter != zoneIDinCurrentProc.end(); ++ iter)
                    {
                        int zoneID = *iter;
                        if (GridID_Neighbor == zoneID)
                        {
                            isNeighborinThisProc = true;
                            break;
                        }
                    }

                    int localIndexNeighbor = zoneIndexMap[GridID_Neighbor];
                    if ((isNeighborinThisProc) && (localIndexNeighbor < iZone))
                    {
                        int interface_Neighbor = interFace2InterFaceID_current[interface_current];
                        UnstructGrid *unstr_Neighbor = UnstructGridCast(grid_in[localIndexNeighbor]);
                        InterfaceInfo *iinfo_Neighbor = unstr_Neighbor->GetInterfaceInfo();
                        int *interFace2BoundaryFace_Neighbor = iinfo_Neighbor->GetInterFace2BoundaryFace();

                        int bcFace_Neighbor = interFace2BoundaryFace_Neighbor[interface_Neighbor];
                        int *face2all_Neighbor = GetFace2All(localIndexNeighbor);
                        face2all[iface] = face2all_Neighbor[bcFace_Neighbor];

                        right_cell_of_all_face[face2all_Neighbor[bcFace_Neighbor]] = cell2all[left_cell_of_face_current[iface]];

                        isFaceMerged[face2all_Neighbor[bcFace_Neighbor]] = true;
                        node_sum += node_number_of_iface;

                        ++ interface_current;
                        continue;
                    }

                    ++ interface_current;
                }
            }

            face2all[iface] = globalFace_count;

            node_number_of_all_face[globalFace_count] = node_number_of_iface;

            for (int iNode = 0; iNode < node_number_of_iface; ++ iNode)
            {
                int nodeindex = face2node_current[node_sum + iNode];
                allface2node.push_back(node2all[nodeindex]);
            }
            node_sum += node_number_of_iface;
            
            if (left_cell_of_face_current[iface] >= nTotalCellEachZone[iZone] || left_cell_of_face_current[iface] < 0)
            {
                left_cell_of_all_face[globalFace_count] = -1;
            }
            else
            {
                left_cell_of_all_face[globalFace_count] = cell2all[left_cell_of_face_current[iface]];
            }
            

            if (right_cell_of_face_current[iface] >= nTotalCellEachZone[iZone] || right_cell_of_face_current[iface] < 0)
            {
                right_cell_of_all_face[globalFace_count] = -1;
            }
            else
            {
                right_cell_of_all_face[globalFace_count] = cell2all[right_cell_of_face_current[iface]];
            }
            
            ++ globalFace_count;
        }

        nTotalFace = globalFace_count;
    }
}

void CombinGrid::MergeCell2All()
{
    WriteLogFile("Merge cell2all ...");
    PrintToWindow("Merge cell2all ...\n");

    int sum = 0;
    nTotalCell = 0;

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int *cell2all = GetCell2All(iZone);

        for (int i = 0; i < nTotalCellEachZone[iZone]; ++ i)
        {
            cell2all[i] = sum;
            sum ++;
        }

        nTotalCell = sum;
    }

    int minNTCell, maxNTCell;
    PH_Reduce(&nTotalCell, &minNTCell, 1, MPI_MIN);
    PH_Reduce(&nTotalCell, &maxNTCell, 1, MPI_MAX);

    if (PHMPI::GetCurrentProcessorID() == PHMPI::GetServerProcessorID())
    {
        ostringstream oss;
        oss << "Min && Max number of cell after merged : " << minNTCell << ", " << maxNTCell << endl;
        oss << "Min to Max nTCell ratio = " << (minNTCell * 1.0) / (maxNTCell * 1.0) << endl;
        PrintToWindow(oss);
        WriteLogFile(oss);
    }

    PH_Barrier();
}

void CombinGrid::MergeNode2All()
{
    WriteLogFile("Merge node2all ...");
    PrintToWindow("Merge node2all ...\n");

    int sum = 0;
    Globalx.reserve(nTotalNode);
    Globaly.reserve(nTotalNode);
    Globalz.reserve(nTotalNode);

    nTotalNode = 0;
    RDouble minDistance, Disx, Disy, Disz;
    minDistance = LARGE;
    bool discard;
    vector <int> bc_node_all;

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int *node2all = GetNode2All(iZone);

        UnstructGrid *unstr = UnstructGridCast(grid_in[iZone]);
        RDouble *x2 = grid_in[iZone]->GetX();
        RDouble *y2 = grid_in[iZone]->GetY();
        RDouble *z2 = grid_in[iZone]->GetZ();
        int *node_number_of_each_face = unstr->GetNodeNumberOfEachFace();
        int *face2node = unstr->GetFace2Node();
        minDistance = MIN(minDistance, grid_in[iZone]->CalMinEdgeLength());
        RDouble distanceTolerance = minDistance / 10;

        int nBoundFace = grid_in[iZone]->GetNBoundFace();

        int *isBCNode_iZone = new int[nTotalNodeEachZone[iZone]];
        PHSPACE::SetField(isBCNode_iZone, 0, nTotalNodeEachZone[iZone]);

        int nsum = 0;

        for (int iface = 0; iface < nBoundFace; ++ iface)
        {
            for (int iNode = 0; iNode < node_number_of_each_face[iface]; ++ iNode)
            {
                isBCNode_iZone[face2node[nsum+iNode]] = 1;
            }
            nsum += node_number_of_each_face[iface];
        }
        
        uint_t bc_node_number = bc_node_all.size();
        for (int iNode = 0; iNode < nTotalNodeEachZone[iZone]; ++ iNode)
        {
            discard = false;

            if (IsBcNode(isBCNode_iZone, iNode))
            {

                for (int jbcNode = 0; jbcNode < bc_node_number; jbcNode ++)
                {
                    Disx = fabs(x2[iNode] - Globalx[bc_node_all[jbcNode]]);
                    Disy = fabs(y2[iNode] - Globaly[bc_node_all[jbcNode]]);
                    Disz = fabs(z2[iNode] - Globalz[bc_node_all[jbcNode]]);

                    if (Disx < distanceTolerance && Disy < distanceTolerance && Disz < distanceTolerance)
                    {
                        node2all[iNode] = bc_node_all[jbcNode];
                        discard = true;
                        break;
                    }
                    else
                    {
                        continue;
                    }
                }
                if (discard != true)
                {
                    bc_node_all.push_back(sum);
                }
            }
            
            if (discard) continue;

            Globalx.push_back(x2[iNode]);
            Globaly.push_back(y2[iNode]);
            Globalz.push_back(z2[iNode]);
            node2all[iNode] = sum;
            sum ++ ;
        }
        nTotalNode = sum;
        delete [] isBCNode_iZone;
    }
}

void CombinGrid::MergeBCType()
{
    WriteLogFile("Merge boundary condition ...");
    PrintToWindow("Merge boundary condition ...\n");

    int nTFace_all = grid_all->GetNTotalFace();

    int nBFace_all = 0;
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        UnstructGrid *unstr = UnstructGridCast(grid_in[iZone]);
        int nBoundFace = unstr->GetNBoundFace();
        int *face2all = GetFace2All(iZone);

        CompositeUnstructBC **bcr = unstr->GetBCRecord();

        for (int iface = 0; iface < nBoundFace; ++ iface)
        {
            if (bcr[iface]->GetKey() == PHENGLEI::INTERFACE)
            {
                int globalFaceID = face2all[iface];
                if (isFaceMerged[globalFaceID])
                {
                    //! The interface has been merged, it no longer be boundary face now.
                    continue;
                }
            }
            ++ nBFace_all;
        }
    }

    grid_all->SetNBoundFace(nBFace_all);

    CompositeUnstructBC **bcr_all = new CompositeUnstructBC * [nTFace_all];

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        UnstructGrid *unstr = UnstructGridCast(grid_in[iZone]);
        int nBoundFace = unstr->GetNBoundFace();
        int *face2all = GetFace2All(iZone);

        CompositeUnstructBC **bcr = unstr->GetBCRecord();

        for (int iface = 0; iface < nBoundFace; ++ iface)
        {
            int globalFaceID = face2all[iface];
            if (bcr[iface]->GetKey() == PHENGLEI::INTERFACE)
            {
                if (isFaceMerged[globalFaceID])
                {
                    //! The interface has been merged, it no longer be boundary face now.
                    continue;
                }
            }

            bcr_all[globalFaceID] = new CompositeUnstructBC();
            bcr_all[globalFaceID]->SetKey(bcr[iface]->GetKey());
            bcr_all[globalFaceID]->SetBCName(bcr[iface]->GetBCName());
        }
    }

    vector< vector<int> > vface2node;
    vface2node.reserve(nTFace_all);
    int *node_number_of_each_face = grid_all->GetNodeNumberOfEachFace();
    int *face2node = grid_all->GetFace2Node();
    int *left_cell_of_face  = grid_all->GetLeftCellOfFace();
    int *right_cell_of_face = grid_all->GetRightCellOfFace();

    int nsum = 0;
    for (int iface = 0; iface < nTFace_all; ++ iface)
    {
        vector< int > eachface2node;
        eachface2node.reserve(node_number_of_each_face[iface]);
        for (int j = 0; j < node_number_of_each_face[iface]; ++ j)
        {
            eachface2node.push_back(face2node[nsum]);
            ++ nsum;
        }
        vface2node.push_back(eachface2node);
    }

    //! Swap faces make sure boundary faces are in the first section (iFace < nBoundFace),
    //! and the interior faces are in the second section (iFace >= nBoundFace).
    face2NewFaceafterReorder = new int [nTFace_all];
    for (int iFace = 0; iFace < nTFace_all; ++ iFace)
    {
        face2NewFaceafterReorder[iFace] = iFace;
    }

    int jstart = nBFace_all;
    for (int iface = 0; iface < nBFace_all; ++ iface)
    {
        if (right_cell_of_face[iface] != -1)
        {
            //! iFace is interior face, it must be reject out of nBoundFace.
            for (int jface = jstart; jface < nTFace_all; ++ jface)
            {
                if (right_cell_of_face[jface] == -1)
                {
                    SWAP(left_cell_of_face[iface], left_cell_of_face[jface]);
                    SWAP(right_cell_of_face[iface], right_cell_of_face[jface]);
                    SWAP(vface2node[iface], vface2node[jface]);
                    SWAP(node_number_of_each_face[iface], node_number_of_each_face[jface]);
                    SWAP(bcr_all[iface], bcr_all[jface]);

                    face2NewFaceafterReorder[iface] = jface;
                    face2NewFaceafterReorder[jface] = iface;

                    jstart = jface + 1;
                    break;
                }
            }
        }
    }

    nsum = 0;
    for (int iface = 0; iface < nTFace_all; ++ iface)
    {
        for (int inode = 0; inode < node_number_of_each_face[iface]; ++ inode)
        {
            face2node[nsum] = vface2node[iface][inode];
            ++ nsum;
        }
    }

    CompositeUnstructBC **bcr_all2 = new CompositeUnstructBC * [nBFace_all];
    grid_all->SetBCRecord(bcr_all2);

    WriteLogFile("nTFace_all = ", nTFace_all);
    WriteLogFile("nBFace_all = ", nBFace_all);

    face2InterfaceafterReorder = new int [nTFace_all];
    PHSPACE::SetField(face2InterfaceafterReorder, -1, nTFace_all);

    interface2BoundaryFaceafterReorder = new int [nBFace_all]; // nIFace of grid_all has not been set yet.
    PHSPACE::SetField(interface2BoundaryFaceafterReorder, -1, nBFace_all);

    int interfaceID = 0;
    for (int iface = 0; iface < nBFace_all; ++ iface)
    {
        bcr_all2[iface] = bcr_all[iface];

        if (bcr_all[iface]->GetKey() == PHENGLEI::INTERFACE)
        {
            interface2BoundaryFaceafterReorder[interfaceID] = iface;
            face2InterfaceafterReorder[iface] = interfaceID ++;
        }
    }
    delete [] bcr_all;

    nIFace = interfaceID;
    WriteLogFile("nIFace_all = ", nIFace);
    InterfaceInfo *iinfo = new InterfaceInfo(nIFace);
    grid_all->SetInterfaceInfo(iinfo);
}

bool CombinGrid::IsBcNode(int *p, int node_in)
{
    bool bc_node = false;

    if (p[node_in] == 1) bc_node = true;
    return bc_node;
}

void CombinGrid::InitNode2All()
{
    node2all = new int *[nZones];
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int nTotalNode = grid_in[iZone]->GetNTotalNode();
        node2all[iZone] = new int [nTotalNode];
    }
}

void CombinGrid::InitFace2All()
{
    face2all = new int *[nZones];
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int nTotalFace = grid_in[iZone]->GetNTotalFace();
        face2all[iZone] = new int [nTotalFace];
    }
}

void CombinGrid::InitIsBCFaceMerged()
{

}

void CombinGrid::InitCell2All()
{
    cell2all = new int *[nZones];
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int nTotalCell = grid_in[iZone]->GetNTotalCell();
        cell2all[iZone] = new int [nTotalCell];
    }
}

void CombinGrid::ResetCoordinate()
{
    WriteLogFile("Combine coordinates ...\n");
    PrintToWindow("Combine coordinates ...\n");

    int nTNode_all = grid_all->GetNTotalNode();

    RDouble *x_all = new RDouble[nTNode_all];
    RDouble *y_all = new RDouble[nTNode_all];
    RDouble *z_all = new RDouble[nTNode_all];

    vector<RDouble> Globalx = GetGlobalx();
    vector<RDouble> Globaly = GetGlobaly();
    vector<RDouble> Globalz = GetGlobalz();

    grid_all->SetX(x_all);
    grid_all->SetY(y_all);
    grid_all->SetZ(z_all);

    for (int inode = 0; inode < nTNode_all; ++ inode)
    {
        x_all[inode] = Globalx[inode];
        y_all[inode] = Globaly[inode];
        z_all[inode] = Globalz[inode];
    }

    grid_all->ComputeMinMaxBox();
}

void CombinGrid::ResetNodeNumberOfEachFace()
{
    WriteLogFile("Combine node_number_of_each_face ...\n");
    PrintToWindow("Combine node_number_of_each_face ...\n");

    int nTFace_all = grid_all->GetNTotalFace();

    int *node_number_of_each_face = new int[nTFace_all];

    grid_all->SetNodeNumberOfEachFace(node_number_of_each_face);

    int *node_number_of_all_face = GetNodeNumberOfAllFace();

    for (int inode = 0; inode < nTFace_all; ++ inode)
    {
        node_number_of_each_face[inode] = node_number_of_all_face[inode];
    }
}

void CombinGrid::ResetFace2Node()
{
    PrintToWindow("Combine face2node ...\n");
    WriteLogFile("Combine face2node ...\n");

    vector<int> allface2node = GetAllFace2Node();
    uint_t nsum = allface2node.size();

    int *face2node = new int [nsum];
    grid_all->SetFace2Node(face2node);
    
    for (int i = 0; i < nsum; ++ i)
    {
        face2node[i] = allface2node[i];
    }
}

void CombinGrid::ResetLeftRightCell()
{
    WriteLogFile("Reset Left Right Cell ...");
    PrintToWindow("Reset Left Right Cell ...\n");

    int nTFace_all = grid_all->GetNTotalFace();
    int *left_cell_of_all_face  = GetLeftCellOfAllFace();
    int *right_cell_of_all_face = GetRightCellOfAllFace();

    int *left_cell_of_face  = new int[nTFace_all];
    int *right_cell_of_face = new int[nTFace_all];
    grid_all->SetLeftCellOfFace (left_cell_of_face);
    grid_all->SetRightCellOfFace(right_cell_of_face);

    for (int iface = 0; iface < nTFace_all; ++ iface)
    {
        left_cell_of_face[iface]  = left_cell_of_all_face[iface];
        right_cell_of_face[iface] = right_cell_of_all_face[iface];
    }
}

void CombinGrid::MergeCell2Node()
{
    if (UnstructGridCast(grid_in[0])->IsCell2NodeExist() == false)
    {
        //! Judge if need merge cell2node.
        return;
    }
    WriteLogFile("Merge cell to node ...");
    PrintToWindow("Merge cell to node ...\n");

    vector< int > node_number_of_each_cell_all;
    node_number_of_each_cell_all.resize(0);
    node_number_of_each_cell_all.reserve(nTotalCell);

    vector< vector<cgsize_t> > cell2node_all;
    cell2node_all.resize(0);
    cell2node_all.reserve(nTotalCell);

    int cell2node_size = 0;

    for (int izone = 0; izone < nZones; ++ izone)
    {
        int *gnode_number_of_each_cell = UnstructGridCast(grid_in[izone])->GetNodeNumberOfEachCell();
        int *gcell2node = UnstructGridCast(grid_in[izone])->GetCell2Node();

        int *gnode2all = GetNode2All(izone);

        int count = 0;
        for (int icell = 0; icell < nTotalCellEachZone[izone]; ++ icell)
        {
            node_number_of_each_cell_all.push_back(gnode_number_of_each_cell[icell]);
            cell2node_size += gnode_number_of_each_cell[icell];

            vector<cgsize_t> single_cell;
            single_cell.resize(0);
            single_cell.reserve(gnode_number_of_each_cell[icell]);

            for (int inode = 0; inode < gnode_number_of_each_cell[icell]; ++ inode)
            {
                int node_index = gcell2node[count];
                ++ count;
                single_cell.push_back(gnode2all[node_index]);
            }
            cell2node_all.push_back(single_cell);
        }
    }

    int *node_number_of_each_cell = new int [nTotalCell];
    int *cell2node = new int [cell2node_size];

    int count = 0;
    for (cgsize_t icell = 0; icell < nTotalCell; ++ icell)
    {
        node_number_of_each_cell[icell] = node_number_of_each_cell_all[icell];
        for (int inode =0; inode < node_number_of_each_cell[icell]; ++ inode)
        {
            cell2node[count] = static_cast<int>(cell2node_all[icell][inode]);
            ++ count;
        }
    }
    grid_all->SetNodeNumberOfEachCell(node_number_of_each_cell);
    grid_all->SetCell2Node(cell2node);
}

Grid ** CombinGrid::GetGridAll() const
{
    Grid **grid_array = new Grid * [1];

    grid_array[0] = grid_all;

    return grid_array;
}

void CombinGrid::DumpCombinationGrid()
{
    string out_gfile;
    GlobalDataBase::GetData("out_gfile", &out_gfile, PHSTRING, 1);

    int codeOfOversetGrid = 0;
    GlobalDataBase::GetData("codeOfOversetGrid", &codeOfOversetGrid, PHINT, 1);

    using namespace PHMPI;

    int nZones = 1;

    Grid **grid_array = new Grid * [nZones];
    for (int izone = 0; izone < nZones; ++ izone)
    {
        grid_array[izone] = grid_all;
    }

    ConstructGlobalInterfaceLink(grid_array, nZones);
    DumpGrid(out_gfile, grid_array, nZones);

    delete [] grid_array;
}

void CombinGrid::CompressInterfaceData(int zoneID, DataContainer *sendBuffer, int neighborZone)
{
    int localZoneID = region->GetProcessGlobalZoneIndexToLocalZoneIndex(zoneID);
    Grid *grid = grid_in[localZoneID];
    InterfaceInfo *iinfo = grid->GetInterfaceInfo();
    if (!iinfo) return;

    int myid = PHMPI::GetCurrentProcessorID();
    int iNeighborZone                            = iinfo->FindIthNeighbor(neighborZone);
    int interfaceNumberBetweenTwoNeighboringZone = iinfo->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForSend          = iinfo->GetFaceIndexForSend(iNeighborZone);
    int *localInterfaceToBoundaryFaceIndexesMapping = iinfo->GetInterFace2BoundaryFace();
    sendBuffer->moveto_begin();
    
//  WriteLogFile("currentZone = ", zoneID);
//  WriteLogFile("neighborZone = ", neighborZone);

    int *face2all = GetFace2All(localZoneID);
    for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
    {
        int interFaceID = interfaceIndexContainerForSend[iLocalFace];
        int boundaryFaceID  = localInterfaceToBoundaryFaceIndexesMapping[interFaceID];
        int oldGlobalFaceID = face2all[boundaryFaceID];
        int newGlobalFaceID = face2NewFaceafterReorder[oldGlobalFaceID];

        //! interFace2ZoneID.
        int newZoneID = myid;
        PHWrite(sendBuffer, newZoneID);

        //! interFace2InterFaceID.
        int newInterfaceID = face2InterfaceafterReorder[newGlobalFaceID];
        PHWrite(sendBuffer, newInterfaceID);
    }
}

void CombinGrid::DecompressInterfaceData(int zoneID, DataContainer *recieveBuffer, int neighborZone)
{
    int localZoneID = region->GetProcessGlobalZoneIndexToLocalZoneIndex(zoneID);
    Grid *grid = grid_in[localZoneID];
    InterfaceInfo *iinfo = grid->GetInterfaceInfo();
    if (!iinfo) return;

    //! Interface information of local sub-zones.
    int iNeighborZone                            = iinfo->FindIthNeighbor(neighborZone);
    int interfaceNumberBetweenTwoNeighboringZone = iinfo->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForReceive       = iinfo->GetFaceIndexForRecv(iNeighborZone);
    int *localInterfaceToBoundaryFaceIndexesMapping = iinfo->GetInterFace2BoundaryFace();

    //! Interface information of global merged zone.
    InterfaceInfo *globalInterfaceInfor = grid_all->GetInterfaceInfo();
    int *interFace2ZoneID       = globalInterfaceInfor->GetInterFace2ZoneID();
    int *interFace2InterFaceID  = globalInterfaceInfor->GetInterFace2InterFaceID();
    int *interFace2BoundaryFace = globalInterfaceInfor->GetInterFace2BoundaryFace();

    recieveBuffer->moveto_begin();
    int *face2all = GetFace2All(localZoneID);
    for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
    {
        int interFaceID = interfaceIndexContainerForReceive[iLocalFace];
        int boundaryFaceID = localInterfaceToBoundaryFaceIndexesMapping[interFaceID];
        int oldGlobalFaceID = face2all[boundaryFaceID];
        int newGlobalFaceID = face2NewFaceafterReorder[oldGlobalFaceID];
        int newInterfaceID  = face2InterfaceafterReorder[newGlobalFaceID];

        //! interFace2ZoneID.
        int newNeighborZoneID;
        PHRead(recieveBuffer, newNeighborZoneID);

        //! interFace2InterFaceID
        int newNeighborInterfaceID;
        PHRead(recieveBuffer, newNeighborInterfaceID);

        if (newInterfaceID >= 0)
        {
            //! This interface has not been merged, so it also be interface in merged grid.
            interFace2ZoneID[newInterfaceID]       = newNeighborZoneID;
            interFace2InterFaceID[newInterfaceID]  = newNeighborInterfaceID;
            interFace2BoundaryFace[newInterfaceID] = newGlobalFaceID;
        }
    }
}

void CombinGrid::SwapInterfaceData()
{
    if (PHMPI::GetNumberOfProcessor() == 1) return;

    PrintToWindow("Swapping interface data ...\n");

    int currentProcessor = PHMPI::GetCurrentProcessorID();
    int numberOfGlobalZones = PHMPI::GetNumberofGlobalZones();
    vector< PH_Request > requestContainer;
    vector< vector< DataContainer * > > receivedDataBuffer;
    vector< vector< DataContainer * > > sendDataBuffer;
    receivedDataBuffer.resize(numberOfGlobalZones);
    sendDataBuffer.resize(numberOfGlobalZones);

    //! Step 0: Compressing data firstly.
    WriteLogFile("        Compressing data ...");
    for (int iZone = 0; iZone < numberOfGlobalZones; ++ iZone)
    {
        int sendProcessor = PHMPI::GetZoneProcessorID(iZone);

        ZoneNeighbor *globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
        std::size_t numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();
        if (currentProcessor == sendProcessor)
        {
            sendDataBuffer[iZone].resize(numberOfNeighbor);
        }

        for (std::size_t iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
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

                //! Compress the send information into the actkey.
                CompressInterfaceData(iZone, sendBuffer, neighborZone);
            }
        }
    }

    //! Step 1: Communication.
    WriteLogFile("        Communication ...");
    for (int iZone = 0; iZone < numberOfGlobalZones; ++ iZone)
    {
        ZoneNeighbor *globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
        std::size_t numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();

        //! Communicating.
        for (std::size_t iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
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
                receivedDataBuffer[ iZone ].push_back(receivedBuffer);
            }
            if (currentProcessor == sendProcessor)
            {
                sendBuffer = sendDataBuffer[iZone][iNeighbor];
            }

            int tag = iZone;

            //if (sendProcessor == receiveProcessor) continue;

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

                for (int iNeighbor = 0; iNeighbor < numberOfNeighborTemp; ++ iNeighbor)
                {
                    int neighborZoneTemp = globalNeighborZonesTemp->GetZoneIndexOfNeighbor(iNeighbor);
                    if (neighborZoneTemp == iZone)
                    {
                        neighborOrer = iNeighbor;
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
                    uint_t lastID = receivedDataBuffer[iZone].size();
                    SWAP(receivedDataBuffer[ iZone ][lastID - 1], sendDataBuffer[neighborZone][neighborOrer]);
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
    WriteLogFile("        Translating data ...");
    for (int iZone = 0; iZone < numberOfGlobalZones; ++ iZone)
    {
        ZoneNeighbor *globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
        std::size_t numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();

        int count = 0;
        for (std::size_t iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
        {
            int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);

            int receiveProcessor = PHMPI::GetZoneProcessorID(neighborZone);

            if (currentProcessor == receiveProcessor)
            {
                DataContainer *receiveData = receivedDataBuffer[iZone][count];

                //! Because get into here when currentProcessor==receiveProcessor,
                //! so, use zoneGridSolver of neighborZone, this zone is on the current processor!

                //! Decompress the interface data from datacontainer.
                DecompressInterfaceData(neighborZone, receiveData, iZone);

                ++ count;
            }
        }
    }

    WriteLogFile("        Free buffer data ...");
    //! Step4: Free the buffers.
    for (std::size_t iDim = 0; iDim < receivedDataBuffer.size(); ++ iDim)
    {
        for (std::size_t jDim = 0; jDim < receivedDataBuffer[iDim].size(); ++ jDim)
        {
            delete receivedDataBuffer[iDim][jDim];
        }
    }
    receivedDataBuffer.clear();
    for (std::size_t iDim = 0; iDim < sendDataBuffer.size(); ++ iDim)
    {
        for (std::size_t jDim = 0; jDim < sendDataBuffer[ iDim ].size(); ++ jDim)
        {
            delete sendDataBuffer[iDim][jDim];
        }
    }
    sendDataBuffer.clear();
    requestContainer.clear();
}

void CombinGrid::DumpMergedGrid(const string &gridFileName, Grid **grids, int nBlocks)
{
    PrintToWindow("    Writing grid file ....\n");

    if (!PHMPI::IsParallelRun())
    {
        TK_Exit::ExceptionExit("fileMode == PHSPACE::CS_MODE has not been considered!");
    }
    else
    {
        DumpMergedGridP2P(gridFileName, grids, nBlocks);
    }
}

void CombinGrid::DumpMergedGridP2P(const string &gridFileName, Grid **grids, int nBlocks)
{
    using namespace PHMPI;

    //! Notice: Directly writing to grid file will affect the number of grid blocks,so you can not dump it at will.
    SetNumberofGlobalZones(nBlocks);

    int nZones = GetNumberofGlobalZones();

    CreateZoneProcessorID(nZones);
    CreateZoneGridID(nZones);
    CreateZoneGridType(nZones);

    int *block_proc     = GetZoneProcessorID();
    int *block_idx      = GetZoneGridID();
    int *block_type     = GetZoneGridType();
    int *block_proc_inp = GetZoneProcessorID_INP();

    int myid = GetCurrentProcessorID();
    for (int izone = 0; izone < nZones; ++ izone)
    {
        block_proc[izone] = myid;
        block_idx [izone] = myid;
        block_type[izone] = grids[izone]->Type();
    }

    if (block_proc_inp)
    {
        for (int izone = 0; izone < nZones; ++ izone)
        {
            SetZoneProcessorID(izone, block_proc_inp[izone]);
        }
    }

    int m_block_proc = PHMPI::GetZoneDistributionMethod();

    for (int izone = 0; izone < nBlocks; ++ izone)
    {
        fstream file;
        string processorGridFileName = AddSymbolToFileName(gridFileName, '_', PHMPI::GetCurrentProcessorID());
        PHSPACE::OpenFile(file, processorGridFileName, ios_base::out|ios_base::binary|ios_base::trunc);

        int nBlocksLocal = 1;
        file.write(reinterpret_cast<char*>(&nBlocksLocal), sizeof(int));
        int blockIndex = block_idx[izone];
        int blockType  = block_type[izone];
        int blockProc  = block_proc[izone];
        if (m_block_proc == DETERMIN_BY_PARTITION && IsConvertGridToMixGrid())
        {
            int *block_proc_dump = GetZoneProcessorIDDump();
            int blockProcDump = block_proc_dump[izone];
            file.write(reinterpret_cast<char*>(&blockProcDump), nBlocksLocal * sizeof(int));
        }
        else
        {
            file.write(reinterpret_cast<char*>(&blockProc), nBlocksLocal * sizeof(int)); 
        }
        file.write(reinterpret_cast<char*>(&blockIndex), nBlocksLocal * sizeof(int));
        file.write(reinterpret_cast<char*>(&blockType), nBlocksLocal * sizeof(int));

        grids[izone]->WriteGrid(file);

        PHSPACE::CloseFile(file);
    }
}

}


