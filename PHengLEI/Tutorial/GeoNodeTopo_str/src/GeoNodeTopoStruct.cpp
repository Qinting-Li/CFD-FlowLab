//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      GeoFaceTopoStruct.cpp
//! @brief     It calculates the node topology of the structured grid,
//!			   and then extract the node topology information,
//!	           such as node-cell and node-face connected relationship.
//! @author    Meng Liyuan.

#include "TK_Log.h"
#include "Glb_Dimension.h"
#include "Geo_NodeTopo_Struct.h"
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
    ParallelOpenFile(ofile, "./out/StrNodeTopologyOutput.txt", ios::out);

    //! Read in the structured grid.
    //! Update all of parameters of each zone.
    //! Calculate the grid topology by InitGeometry function.
    Region *region = new Region();
    region->ReadGrid();
    region->UpdateAllZonesPara();
    region->InitGeometry();

    //! Extract and dump out the node topology information of the structured grid.
    int nZones = GetNumberofGlobalZones();
    for (int iGrid = 0; iGrid < nZones; ++ iGrid)	
    {
        ofile << " Structured Grid Node Topology of Zone "<< iGrid << " is: " << endl;

        Grid *grid = GetGrid(iGrid);
        int ni = StructGridCast(grid)->GetNI();
        int nj = StructGridCast(grid)->GetNJ();
        int nk = StructGridCast(grid)->GetNK();

        ofile << "   The node dimension of I direction  is: " << ni << endl;
        ofile << "   The node dimension of J direction  is: " << nj << endl;
        ofile << "   The node dimension of K direction  is: " << nk << endl;

        RDouble3D * structx = StructGridCast(grid)->GetStructX();
        RDouble3D * structy = StructGridCast(grid)->GetStructY();
        RDouble3D * structz = StructGridCast(grid)->GetStructZ(); 

        for ( int k = 1; k <= nk; ++ k )
        {
            for ( int j = 1; j <= nj; ++ j )
            {
                for ( int i = 1; i <= ni; ++ i )
                {
                    ofile << "   The x coordinates array of the node is: " << (*structx)(i, j, k) << endl;
                    ofile << "   The y coordinates array of the node is: " << (*structy)(i, j, k) << endl;
                    ofile << "   The z coordinates array of the node is: " << (*structz)(i, j, k) << endl;
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