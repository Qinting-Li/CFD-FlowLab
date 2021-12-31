#include "WallDist.h"
#include "Math_BasisFunction.h"
#include "PHIO.h"
#include "Task_Simple.h"
#include "TK_Exit.h"
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
    int iLES = GlobalDataBase::GetIntParaFromDB("iLES");

    if (GetTaskCode() == CAL_WALL_DIST || (viscousType != INVISCID && IsReadWallDist() == false))
    {
        AllocateWalldist();

        ComputeWallDist();

        WriteWallDist();

        //VisualWallDist();
    }

    if (viscousType > LAMINAR || iLES == LES_SOLVER)
    {
        GhostWallDistance();

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
                Grid *grid0 = GetGrid(iZone, 0);
                if (grid0->Type() == UNSTRUCTGRID)
                {
                    UnstructGrid *grid = UnstructGridCast(grid0);
                    if (grid->GetWallDist()) isWDTExist = 1;
                }
                else
                {
                    StructGrid *grid = StructGridCast(grid0);
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
    delete actkey;
}

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
    delete actkey;
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
                RDouble *walldist = UnstructGridCast(grid)->GetWallDist();
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
                StructGrid *structGrid = StructGridCast(grid);

                RDouble3D &walldist = *structGrid->GetWallDist();

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

void GhostWallDistance()
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

    int level = 0;

    for (int izone = 0; izone < nLocalZones; ++ izone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(izone);

        Grid *grid = GetGrid(zoneID, level);

        if (grid->Type() == PHSPACE::STRUCTGRID)
        {
            StructGrid *structGrid = StructGridCast(grid);
            RDouble3D &walldist = *structGrid->GetWallDist();
        }
    }
}

}