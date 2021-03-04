#include "WallDist.h"
#include "Math_BasisFunction.h"
#include "PHIO.h"
#include "Task_Simple.h"
#include "TK_Exit.h"

#ifdef USE_TecplotLib
#include "TECIO.h"
#endif

#include "Constants.h"
#include "Pre_WalldistCompute.h"
#include "IO_HDF5File.h"
#include "Geo_UnstructGrid.h"

using namespace std;

namespace PHSPACE
{
void GenerateWallDist()
{
    if (GetTaskCode() == DSMC_SOLVER || GetTaskCode() == UGKS_SOLVER || GetTaskCode() == HO_SOLVER) return;
    
    if (PHMPI::CurrentProcessorIsGridProcessor())
    {
        return;
    }

    int viscousType = 0;
    GlobalDataBase::GetData("viscousType", &viscousType, PHINT, 1);

    if (GetTaskCode() == CAL_WALL_DIST || (IsReadWallDist() == false))
    {
        AllocateWalldist();

        ComputeWallDist();

        WriteWallDist();

        //VisalWallDist();
    }

    if (viscousType > 1)
    {
        WallDistStatics();
    }
}

bool IsReadWallDist()
{
    using namespace PHMPI;

    int nZones = GetNumberofGlobalZones();
    int myid   = GetCurrentProcessorID();

    int isWDTExist = 0;
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int procID = GetZoneProcessorID(iZone);

        if (!CurrentProcessorIsGridProcessor())
        {
            if (myid == procID)
            {
                Grid * grid0 = GetGrid(iZone, 0);
                if (grid0->Type() == UNSTRUCTGRID)
                {
                    UnstructGrid *grid = UnstructGridCast(grid0);
                    if (grid->GetWallDist()) isWDTExist = 1;
                }
                else
                {
                    StructGrid * grid = StructGridCast(grid0);
                    if (grid->GetWallDist()) isWDTExist = 1;
                }
            }
        }
    }

    PH_CompareMaxMin(isWDTExist, 1);

    if (isWDTExist == 1)
    {
        return true;
    }

    return false;
}

void ComputeWallDist()
{
    int walldistMethod;
    GlobalDataBase::GetData("walldistMethod", &walldistMethod, PHINT, 1);

    int level = 0;

    Pre_WalldistCompute walldistCompute(walldistMethod, level);

    walldistCompute.Run();
}

void AllocateWalldist()
{
    ActionKey *actkey = new ActionKey();

    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    int myid = GetCurrentProcessorID();

    for (int izone = 0; izone < nZones; ++ izone)
    {
        int proc = GetZoneProcessorIDSepMode(izone);

        if (myid == proc)
        {
            Grid *grid = GetGrid(izone, actkey->level);
            grid->AllocateWalldist(actkey);
        }
    }
}

#ifdef USE_TecplotLib
void VisalWallDist()
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

    ostringstream oss;
    oss << "x y z dist";

    string Variable = oss.str();
    const char *Variables = Variable.c_str();

    string filename = "./in/WallDistVisual.plt";
    const char *pltname = filename.c_str();

    INTEGER4 Debug     = 1;
    INTEGER4 VIsDouble = 0;
    INTEGER4 FileType  = 0;
    INTEGER4 I;

    I = TECINI112((char*)"PHengLEI Grid Videotex",
                  (char*)Variables,
                  (char*)pltname,
                  (char*)".",
                  &FileType,
                  &Debug,
                  &VIsDouble);

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int izone = 0; izone < nLocalZones; ++ izone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(izone);
        UnstructGrid * grid = UnstructGridCast(GetGrid(zoneID, 0));

        INTEGER4 IMxOrNumPts, JMxOrNumElements, KMxOrNumFaces;
        INTEGER4 ZoneType;

        IMxOrNumPts = grid->GetNTotalNode();
        JMxOrNumElements = grid->GetNTotalCell();
        KMxOrNumFaces = grid->GetNTotalFace();

        int dimension = grid->GetDim();

        if (dimension == 2)
        {
            ZoneType = 6;
        }
        else
        {
            ZoneType = 7;
        }

        int GridID = grid->GetGridID()->GetIndex();

        RDouble *x = grid->GetX();
        RDouble *y = grid->GetY();
        RDouble *z = grid->GetZ();

        INTEGER4 ValueLocation[4] = {1, 1, 1, 0};

        INTEGER4 ICellMax           = 0;
        INTEGER4 JCellMax           = 0;
        INTEGER4 KCellMax           = 0;
        double   SolutionTime       = 0.0;
        INTEGER4 StrandID           = 0;
        INTEGER4 ParentZone         = 0;
        INTEGER4 IsBlock            = 1;
        INTEGER4 NumFaceConnections = 0;
        INTEGER4 FaceNeighborMode   = 0;

        INTEGER4 TotalNumFaceNodes_Rect = 0;

        //int * cell2node = grid->GetCell2Node();

        int * faceNodeNumberContainer = grid->GetNodeNumberOfEachFace();
        int * faceNodeIndexContainer  = grid->GetFace2Node();
    
        for (int iface = 0; iface < KMxOrNumFaces; ++ iface)
        {
            TotalNumFaceNodes_Rect += faceNodeNumberContainer[iface];
        }

        INTEGER4 *FaceNodes_Rect = new INTEGER4[TotalNumFaceNodes_Rect];
        for (INTEGER4 i = 0; i < TotalNumFaceNodes_Rect; ++ i)
        {
            FaceNodes_Rect[i] = faceNodeIndexContainer[i] + 1;
        }

        INTEGER4 NumConnBndryFaces_Rect  = 0;    /* interface */
        INTEGER4 TotalNumBndryConns_Rect = 0;
        INTEGER4 SharConn                = 0;

        ostringstream oss;

        //Zonetitle = Zone + zonenum;
        oss << "Zone"<< GridID;
    
        //const char * ZoneTitle = Zonetitle.c_str();
        string Zonetitle = oss.str();
        const char * ZoneTitle = Zonetitle.c_str();

        I = TECZNE112((char*)ZoneTitle,
                          &ZoneType,
                          &IMxOrNumPts,
                          &JMxOrNumElements,
                          &KMxOrNumFaces,
                          &ICellMax,
                          &JCellMax,
                          &KCellMax,
                          &SolutionTime,
                          &StrandID,
                          &ParentZone,
                          &IsBlock,
                          &NumFaceConnections,
                          &FaceNeighborMode,
                          &TotalNumFaceNodes_Rect,
                          &NumConnBndryFaces_Rect,
                          &TotalNumBndryConns_Rect,
                          NULL,
                          ValueLocation,
                          NULL,
                          &SharConn);

        INTEGER4 IsDouble = 1;
        I = TECDAT112(&IMxOrNumPts, x, &IsDouble);
        I = TECDAT112(&IMxOrNumPts, y, &IsDouble);
        I = TECDAT112(&IMxOrNumPts, z, &IsDouble);

        RDouble * walldist = grid->GetWallDist();
        I = TECDAT112(&JMxOrNumElements, walldist, &IsDouble);

        INTEGER4 *FaceLeftElems  = new INTEGER4[KMxOrNumFaces];
        INTEGER4 *FaceRightElems = new INTEGER4[KMxOrNumFaces];

        int * leftCellIndexContainer  = grid->GetLeftCellOfFace();
        int * rightCellIndexContainer = grid->GetRightCellOfFace();

        for (INTEGER4 iface = 0; iface < KMxOrNumFaces; ++ iface)
        {
            FaceLeftElems [iface] = leftCellIndexContainer [iface] + 1;
            FaceRightElems[iface] = rightCellIndexContainer[iface] + 1;
            if (FaceRightElems[iface] > JMxOrNumElements || FaceRightElems[iface] < 0) FaceRightElems[iface] = 0;
            if (FaceLeftElems[iface] > JMxOrNumElements || FaceLeftElems[iface] < 0) FaceLeftElems[iface] = 0;
        }

        INTEGER4 FaceBndryConnectionCounts = 0;
                 INTEGER4 FaceBndryConnectionElems = 0;
                 INTEGER4 FaceBndryConnectionZones = 0;

        I = TECPOLY112(faceNodeNumberContainer,
                       FaceNodes_Rect,
                       FaceLeftElems,
                       FaceRightElems,
                       &FaceBndryConnectionCounts,
                       &FaceBndryConnectionElems,
                       &FaceBndryConnectionZones);

        delete [] FaceLeftElems;
        delete [] FaceRightElems;
        delete [] FaceNodes_Rect;
    }

    I = TECEND112();
}
#endif

void WriteWallDist()
{
    string walldistfile = "walldist.dat";
    GlobalDataBase::GetData("walldistfile", &walldistfile, PHSTRING, 1);

    ActionKey *actkey = new ActionKey();
    actkey->filename = (walldistfile);

    ParallelCreateHDF5File(actkey);

    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    for (int izone = 0; izone < nZones; ++ izone)
    {
        int send_proc = GetZoneProcessorID(izone);
        int recv_proc = GetZoneFileID(izone);

        int myid = GetCurrentProcessorID();
        int tag = GetSendRecvTag(actkey, izone);

        if (myid == send_proc)
        {
            Grid *grid = GetGrid(izone, 0);
            grid->DumpWalldist(actkey);
        }

        PH_Trade(actkey, send_proc, recv_proc, tag);

        if (myid == recv_proc)
        {
            DumpWallDistToFile(actkey, izone);
        }
    }

    ParallelCloseHDF5File(actkey);
}

void DumpWallDistToFile(ActionKey *actkey, int iZone)
{
    DataContainer *cdata = actkey->GetData();
    cdata->moveto_begin();

    int nTotalCell = 0;
    cdata->read(&nTotalCell, sizeof(int));

    RDouble *walldist = new RDouble[nTotalCell];
    cdata->read(walldist, sizeof(RDouble) * nTotalCell);

    ostringstream oss;
    oss << "Grid-" << iZone;
    string dataName = oss.str();

    CreateAndWriteData(actkey->filepos, dataName, 1, nTotalCell, PHDOUBLE, walldist);

    /*test
    ostringstream oss2;
    oss2 << nTotalCell << endl;
    oss2.setf(ios::scientific);
    oss2.precision(10);
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        oss2 << iCell << " " << walldist[iCell] << endl;
    }
    WriteLogFile(oss2);
    //*/

    delete [] walldist;
}

void WallDistStatics()
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();
    int myid = GetCurrentProcessorID();

    RDouble minWalldist = PHSPACE::LARGE;
    RDouble maxWalldist = PHSPACE::SMALL;
    
    for (int izone = 0; izone < nZones; ++ izone)
    {
        int proc = GetZoneProcessorIDSepMode(izone);

        if (myid == proc)
        {
            Grid *grid = GetGrid(izone, 0);

            if (grid->Type() == UNSTRUCTGRID)
            {
                int nTotalCell = grid->GetNTotalCell();
                RDouble * walldist = UnstructGridCast(grid)->GetWallDist();
                for (int iCell = 0; iCell < nTotalCell; ++ iCell)
                {
                    minWalldist = MIN(minWalldist, walldist[iCell]);
                    maxWalldist = MAX(maxWalldist, walldist[iCell]);
                    if (maxWalldist > 1.0e8)
                    {
                        cout << iCell << endl;
                        ostringstream oss1;
                        oss1 << "iCell=" << iCell << endl;
                        TK_Exit::ExceptionExit(oss1);
                    }
                }
            }
            else
            {
                StructGrid * structGrid = StructGridCast(grid);

                RDouble3D & walldist = *structGrid->GetWallDist();

                int ist, ied, jst, jed, kst, ked;
                structGrid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

                for (int k = kst; k <= ked; ++ k)
                {
                    for (int j = jst; j <= jed; ++ j)
                    {
                        for (int i = ist; i <= ied; ++ i)
                        {
                            minWalldist = MIN(minWalldist, walldist(i, j, k));
                            maxWalldist = MAX(maxWalldist, walldist(i, j, k));
                        }
                    }
                }
            }
        }
    }

    PH_CompareMaxMin(minWalldist, 2);

    PH_CompareMaxMin(maxWalldist, 1);

    ostringstream oss;
    oss << "Min && Max wall distance : " << minWalldist << ", " << maxWalldist << endl;
    PrintToWindow( oss );
}

}