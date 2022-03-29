//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      GeoMultiGridInfoStruct.cpp
//! @brief     It calculates and extracts Multi-grid method information of structured grid.
//! @author    Meng Liyuan.

#include "TK_Log.h"
#include "Glb_Dimension.h"
#include "Geo_MultiGridInfo_Struct.h"
#include "Region.h"
#include "PHIO.h"
#include "TK_Time.h"
#include "GlobalDataBase.h"

using namespace PHSPACE;
using namespace PHMPI;
int main()
{ 
    //! Perform the initialization operation of PHengLEI.
    ReadControlParameters();
    Initialization();

    //! Define output file which to store the results.
    fstream ofile;
    ParallelOpenFile(ofile, "./out/StrMultiGridInfoOutput.txt", ios::out);	

    //! Read in the unstructured grid.
    //! Update all of parameters of each zone.
    //! Calculate the Multi-grid method information by InitGeometry function.
    Region *region = new Region();
    region->ReadGrid();
    region->UpdateAllZonesPara();
    region->InitGeometry();

    //! Extract and dump out the Multi-grid method information of the structured grid.
    int nZones = GetNumberofGlobalZones();
    for (int iGrid = 0; iGrid < nZones; ++ iGrid)	
    {
        ofile << "Structured MultiGridInfo of Zone "<< iGrid << " is: " << endl;

        Grid *grid = GetGrid(iGrid);
        int   istp = StructGridCast(grid)->GetMultigridStepI();
        int   jstp = StructGridCast(grid)->GetMultigridStepJ();
        int   kstp = StructGridCast(grid)->GetMultigridStepK();
        ofile << "   The multigrid method step I is:  " << istp << endl;
        ofile << "   The multigrid method step J is:  " << jstp << endl;
        ofile << "   The multigrid method step K is:  " << kstp << endl;
    }

    //! Close the output file and exit PHengLEI.
    ParallelCloseFile(ofile);

    PHMPI::FreeBlockData();
    TIME_SPACE::WriStatTime();
    PH_Finalize();
    return 0;
}