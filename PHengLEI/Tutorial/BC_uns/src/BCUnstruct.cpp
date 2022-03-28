//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      BCUnstruct.cpp
//! @brief     It tests the class BCRecord of unstructured grid/solver Boundary Condition, 
//!			   check whether it can run as required or not.
//! @author    Meng Liyuan.

#include "TK_Log.h"
#include "Glb_Dimension.h"
#include "BC_Unstruct.h"
#include "Region.h"
#include "PHIO.h"
#include "GlobalDataBase.h"
#include "TK_Time.h"

using namespace PHSPACE;
using namespace PHMPI;
int main()
{
    //! Record time.
    setvbuf( stdout, NULL, _IONBF, 0 );
    double *ptr_time = NULL;
    PHSPACE::TIME_SPACE::GetWallTimeForResidual();
    PHSPACE::TIME_SPACE::GetWallTime();
    PHSPACE::TIME_SPACE::GetCPUTime(ptr_time);

    //! Perform the initialization operation of PHengLEI.
    ReadControlParameters();
    Initialization();

    //! Define output file which to store the results.
    fstream ofile;
    ParallelOpenFile(ofile, "./out/BC_UnsOutput.txt", ios::out);

    //! Read in the structured grid.
    Region * region = new Region();
    region->ReadGrid();

    //! Extract and dump out the boundary condition information of the unstructured grid.
    int iGrid = 0;
    Grid * grid = GetGrid(iGrid);
    UnstructGrid * UnstructGrid = UnstructGridCast(grid);
    int nBoundFace = UnstructGrid->GetNBoundFace();
    int bctype;
    string bcname;
    BCRecord ** bcr_src = UnstructGrid->GetBCRecord();
    BCRecord ** bcr     = new BCRecord * [nBoundFace];

    for ( int iFace = 0; iFace < nBoundFace; iFace ++ )
    {
        if ( iFace<=100 || iFace>=10000 )
        {
            bctype = -1;
            bcname = "Interface";
        }
        else
        {		
            bctype = bcr_src[iFace]->GetKey();
            bcname = bcr_src[iFace]->GetBCName();
        }
        bcr[iFace] = new BCRecord();
        bcr[iFace]->SetKey( bctype );
        bcr[iFace]->SetBCName(bcname);
        bctype = bcr[iFace]->GetKey();
        bcname = bcr[iFace]->GetBCName();
        ofile << "   The boundary condition type of  "<<iFace<< " is: "<< bctype << endl;
        ofile << "   The boundary condition name of  "<<iFace<< " is: "<< bcname << endl;

        BoundaryCondition * bc = 0;
        bcr_src[iFace]->SetBoundaryCondition(bc);
    }

    //! Close the output file and exit PHengLEI.
    ParallelCloseFile(ofile);

    PHMPI::FreeBlockData();
    TIME_SPACE::WriStatTime();
    PH_Finalize();
    return 0;
}