//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      GeoLSQWeightUnstruct.cpp
//! @brief     It calculates the Least Square Weight of unstructured grid,
//!			   then extract the Weight information.
//! @author    Meng Liyuan.

#include "TK_Log.h"
#include "Glb_Dimension.h"
#include "Geo_LSQWeight_Unstruct.h"
#include "Region.h"
#include "PHIO.h"
#include "GlobalDataBase.h"
#include "TK_Time.h"
#include "WallDist.h"

using namespace PHSPACE;
using namespace PHMPI;
int main()
{
    //! Perform the initialization operation of PHengLEI.
    ReadControlParameters();
    CompleteParamInfo();
    Initialization();

    //! Define output file which to store the results.
    fstream ofile;
    ParallelOpenFile(ofile, "./out/UnsLSQWeightOutput.txt", ios::out);

    //! Read in the unstructured grid.
    //! Update all of parameters of each zone.
    //! Calculate the Least Square Weight by InitGeometry function.
    Region *region = new Region();
    region->ReadGrid();
    region->UpdateAllZonesPara();
    region->InitGeometry();
    region->SwapGeommetry();

    //! Extract and dump out the Least Square Weight information of the unstructured grid.
    int nZones = GetNumberofGlobalZones();
    for (int iGrid = 0; iGrid < nZones; ++ iGrid)	
    {

        ofile << "Unstructured Grid LSQWeight of Zones "<< iGrid << " is: " << endl;

        Grid * grid = GetGrid(iGrid);
        int nTotalCell = grid->GetNTotalCell();
        int nTotalFace = grid->GetNTotalFace();

        RDouble * iwt = UnstructGridCast(grid)->GetLeastSquareIWT();
        RDouble * ixx = UnstructGridCast(grid)->GetLeastSquareIXX();
        RDouble * iyy = UnstructGridCast(grid)->GetLeastSquareIYY();
        RDouble * izz = UnstructGridCast(grid)->GetLeastSquareIZZ();
        RDouble * ixy = UnstructGridCast(grid)->GetLeastSquareIXY();
        RDouble * ixz = UnstructGridCast(grid)->GetLeastSquareIXZ();
        RDouble * iyz = UnstructGridCast(grid)->GetLeastSquareIYZ();

        for (int iFace = 0; iFace < nTotalFace; ++ iFace)
        {
            ofile << "   The LeastSquareIWT of face "<< iFace << " is: " << iwt[iFace] << endl;
        }
        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            ofile << "   The LeastSquareIXX of cell "<< iCell << " is: " << ixx[iCell] << endl;
            ofile << "   The LeastSquareIYY of cell "<< iCell << " is: " << iyy[iCell] << endl;
            ofile << "   The LeastSquareIZZ of cell "<< iCell << " is: " << izz[iCell] << endl;
            ofile << "   The LeastSquareIXY of cell "<< iCell << " is: " << ixy[iCell] << endl;
            ofile << "   The LeastSquareIXZ of cell "<< iCell << " is: " << ixz[iCell] << endl;
            ofile << "   The LeastSquareIYZ of cell "<< iCell << " is: " << iyz[iCell] << endl;
        }
    }

    //! Close the output file and exit PHengLEI.
    ParallelCloseFile(ofile);

    PHMPI::FreeBlockData();
    TIME_SPACE::WriStatTime();
    PH_Finalize();
    return 0;
}