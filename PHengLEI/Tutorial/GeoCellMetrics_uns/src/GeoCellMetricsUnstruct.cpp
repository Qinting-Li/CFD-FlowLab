//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      GeoCellMetricsUnstruct.cpp
//! @brief     It calculates the cell metrics of the unstructured grid,
//!			   and then extract the cell metrics information,
//!	           such as cell center and cell volume.
//! @author    Meng Liyuan.

#include "TK_Log.h"
#include "Glb_Dimension.h"
#include "Geo_CellMetrics_Unstruct.h"
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
    CompleteParamInfo();
    Initialization();

    //! Define output file which to store the results .
    fstream ofile;
    ParallelOpenFile(ofile, "./out/UnsCellMetricsOutput.txt", ios::out);

    //! Read in the unstructured grid .
    //! Update all of parameters of each zone 
    //! Calculate the grid metrics by InitGeometry function.
    Region *region = new Region();
    region->ReadGrid();
    region->UpdateAllZonesPara();
    region->InitGeometry();

    //! Extract and dump out the cell metrics information of the unstructured grid.
    int nZones = GetNumberofGlobalZones();

    for (int iGrid = 0; iGrid < nZones; ++ iGrid)	
    {
        ofile << "Unstructured Grid CellMetrics of Zone "<< iGrid << " is: " << endl;

        Grid * grid = GetGrid(iGrid);
        RDouble * xcc        = UnstructGridCast(grid)->GetCellCenterX();
        RDouble * ycc        = UnstructGridCast(grid)->GetCellCenterY();
        RDouble * zcc        = UnstructGridCast(grid)->GetCellCenterZ();
        RDouble * CellVolume = UnstructGridCast(grid)->GetCellVolume();

        int nTotalCell = grid->GetNTotalCell();
        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            ofile << "   The CellCenterX of "<< iCell << " is: " << xcc[iCell]        << endl;
            ofile << "   The CellCenterY of "<< iCell << " is: " << ycc[iCell]        << endl;
            ofile << "   The CellCenterZ of "<< iCell << " is: " << zcc[iCell]        << endl;
            ofile << "   The CellVolume of " << iCell << " is: " << CellVolume[iCell] << endl;
        }
    }

    //! Close the output file and exit PHengLEI.
    ParallelCloseFile(ofile);

    PHMPI::FreeBlockData();
    TIME_SPACE::WriStatTime();
    PH_Finalize();
    return 0;

}