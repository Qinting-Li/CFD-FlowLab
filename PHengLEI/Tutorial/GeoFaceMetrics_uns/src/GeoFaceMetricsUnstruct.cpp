//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      GeoFaceMetricsUnstruct.cpp
//! @brief     It calculates the face metrics of the unstructured grid,
//!			   and then extract the face metrics information,
//!	           such as face center and face normal.
//! @author    Meng Liyuan.

#include "TK_Log.h"
#include "Glb_Dimension.h"
#include "Geo_FaceMetrics_Unstruct.h"
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

    //! Define output file which to store the results.
    fstream ofile;
    ParallelOpenFile(ofile, "./out/UnsFaceMetricsOutput.txt", ios::out);

    //! Read in the unstructured grid.
    //! Update all of parameters of each zone.
    //! Calculate the grid metrics by InitGeometry function.
    Region *region = new Region();
    region->ReadGrid();
    region->UpdateAllZonesPara();
    region->InitGeometry();

    //! Extract and dump out the face metrics information of the unstructured grid.
    int nZones = GetNumberofGlobalZones();
    for (int iGrid = 0; iGrid < nZones; ++ iGrid)	
    {
        ofile << "Unstructured Grid FaceMetrics of Zone "<< iGrid << " is: " << endl;

        Grid    * grid = GetGrid(iGrid);
        RDouble * xfc  = UnstructGridCast(grid)->GetFaceCenterX();
        RDouble * yfc  = UnstructGridCast(grid)->GetFaceCenterY();
        RDouble * zfc  = UnstructGridCast(grid)->GetFaceCenterZ();
        RDouble * xfn  = UnstructGridCast(grid)->GetFaceNormalX();
        RDouble * yfn  = UnstructGridCast(grid)->GetFaceNormalY();
        RDouble * zfn  = UnstructGridCast(grid)->GetFaceNormalZ();
        RDouble * FaceArea = UnstructGridCast(grid)->GetFaceArea();

        //int nTotalCell = grid->GetNTotalCell();
        int nTotalFace = grid->GetNTotalFace();
        for (int iFace = 0; iFace < nTotalFace; ++ iFace)
        {
            ofile << "   The FaceCenterX of "<< iFace << "is:" << xfc[iFace] << endl;
            ofile << "   The FaceCenterY of "<< iFace << "is:" << yfc[iFace] << endl;
            ofile << "   The FaceCenterZ of "<< iFace << "is:" << zfc[iFace] << endl;
            ofile << "   The FaceNormalX of "<< iFace << "is:" << xfn[iFace] << endl;
            ofile << "   The FaceNormalY of "<< iFace << "is:" << yfn[iFace] << endl;
            ofile << "   The FaceNormalZ of "<< iFace << "is:" << zfn[iFace] << endl;
            ofile << "   The FaceArea of "   << iFace << "is:" << FaceArea[iFace] << endl;
        }
    }

    //! Close the output file and exit PHengLEI.
    ParallelCloseFile(ofile);

    PHMPI::FreeBlockData();
    TIME_SPACE::WriStatTime();
    PH_Finalize();
    return 0;
}