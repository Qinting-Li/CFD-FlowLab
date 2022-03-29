//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      GeoFaceMetricsStruct.cpp
//! @brief     It calculates the face metrics of the structured grid,
//!			   and then extract the face metrics information,
//!	           such as face center and face normal.
//! @author    Meng Liyuan.

#include "TK_Log.h"
#include "Glb_Dimension.h"
#include "Geo_FaceMetrics_Struct.h"
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
    ParallelOpenFile(ofile, "./out/StrFaceMetricsOutput.txt", ios::out);

    //! Read in the structured grid .
    //! Update all of parameters of each zone.
    //! Calculate the grid metrics by InitGeometry function.
    Region *region = new Region();
    region->ReadGrid();
    region->UpdateAllZonesPara();
    region->InitGeometry();


    //! Extract and dump out the face metrics information of the structured grid.
    int nZones = GetNumberofGlobalZones();
    for (int iGrid = 0; iGrid < nZones; ++ iGrid)	
    {
        ofile << "Structured Grid FaceMetrics of Zone "<< iGrid << " is: " << endl;

        Grid * grid = GetGrid(iGrid);
        int ni = StructGridCast(grid)->GetNI();
        int nj = StructGridCast(grid)->GetNJ();
        int nk = StructGridCast(grid)->GetNK();
        int nd = GetDim();
        /* Range I(-1,ni+1);
        Range J(-1,nj+1);
        Range K(-1,nk+1);
        if ( nk == 1 ) K.setRange(1,1);
        Range D(1,3);*/

        RDouble4D * xfn = StructGridCast(grid)->GetFaceNormalX();
        RDouble4D * yfn = StructGridCast(grid)->GetFaceNormalY();
        RDouble4D * zfn = StructGridCast(grid)->GetFaceNormalZ();
        RDouble4D * FaceArea = StructGridCast(grid)->GetFaceArea();
        for ( int d = 1; d <= nd; ++d )
        {
            for ( int k = 1; k <= nk - 1; ++ k )
            {
                for ( int j = 1; j <= nj - 1; ++ j )
                {
                    for ( int i = 1; i <= ni - 1; ++ i )
                    {
                        ofile << "   The X component of UNIT face normal is: " << (*xfn )(i, j, k, d)      << endl;
                        ofile << "   The Y component of UNIT face normal is: " << (*yfn )(i, j, k, d)      << endl;
                        ofile << "   The Z component of UNIT face normal is: " << (*zfn )(i, j, k, d)      << endl;
                        ofile << "   The face area of is:"                     << (*FaceArea )(i, j, k, d) << endl;
                    }
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