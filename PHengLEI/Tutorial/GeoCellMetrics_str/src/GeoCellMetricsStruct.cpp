//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      GeoCellMetricsStruct.cpp
//! @brief     It calculates the cell metrics of the structured grid,
//!			   and then extract the cell metrics information,
//!	           such as cell center and cell volume.
//! @author    Meng Liyuan.
#include "TK_Log.h"
#include "Glb_Dimension.h"
#include "Geo_CellMetrics_Struct.h"
#include "Region.h"
#include "PHIO.h"
#include "GlobalDataBase.h"
#include "TK_Time.h"

using namespace PHSPACE;
using namespace PHMPI;
int main()
{

    //! Perform the initialization operation of PHengLEI.
    ReadControlParameters();
    Initialization();

    //! Define output file which to store the results .
    fstream ofile;
    ParallelOpenFile(ofile, "./out/StrCellMetricsOutput.txt", ios::out);

    //! Read in the structured grid .
    //! Update all of parameters of each zone.
    //! Calculate the grid metrics by InitGeometry function.
    Region *region = new Region();
    region->ReadGrid();
    region->UpdateAllZonesPara();
    region->InitGeometry();

    //! Extract and dump out the cell metrics information of the structured grid.
    int nZones = GetNumberofGlobalZones();
    for (int iGrid = 0; iGrid < nZones; ++ iGrid)	
    {
        ofile << "Structured Grid CellMetrics of Zone "<< iGrid << " is: " << endl;

        Grid * grid = GetGrid(iGrid);
        int ni = StructGridCast(grid)->GetNI();
        int nj = StructGridCast(grid)->GetNJ();
        int nk = StructGridCast(grid)->GetNK();
        Range I(1, ni-1);
        Range J(1, nj-1);
        Range K(1, nk-1);
        if ( nk == 1 ) K.setRange(1, 1);
        RDouble3D & xcc           = *StructGridCast(grid)->GetCellCenterX();
        RDouble3D & ycc           = *StructGridCast(grid)->GetCellCenterY();
        RDouble3D & zcc           = *StructGridCast(grid)->GetCellCenterZ();
        RDouble3D & CellVolume    = *StructGridCast(grid)->GetCellVolume();
        RDouble3D & LargestLocalGridLength = *StructGridCast(grid)->GetLargestLocalGridLength();
        RDouble3D & SubgridLength = *StructGridCast(grid)->GetSubgridLength();

        int ist, ied, jst, jed, kst, ked;
        GetIJKRegion(I, J, K, ist, ied, jst, jed, kst, ked);
        for ( int k = kst; k <= ked; ++ k )       
        {
            for ( int j = jst; j <= jed; ++ j )         
            {
                for ( int i = ist; i <= ied; ++ i )
                {                                                                      
                    ofile << "   The cell center is: "     << xcc(i, j, k) << "  "<< ycc(i, j, k) << "  "<< zcc(i, j, k)<<endl; 
                    ofile << "   The cell volume is: "     << CellVolume(i, j, k)<< endl;
                    ofile << "   The largest local grid spacing is: " << LargestLocalGridLength(i, j, k)<< endl;	
                    ofile << "   The sub-grid length is: " << SubgridLength(i, j, k)<< endl;; 
                } 
            }
        }
    }

    //! Close the output file and exit PHengLEI.
    ParallelCloseFile(ofile);

    PHMPI::FreeBlockData();
    TIME_SPACE::WriStatTime();
    PH_Finalize();
    return 0;

}