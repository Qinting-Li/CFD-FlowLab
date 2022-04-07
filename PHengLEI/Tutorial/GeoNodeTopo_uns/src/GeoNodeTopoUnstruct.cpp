//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      GeoFaceTopoUnstruct.cpp
//! @brief     It calculates the node topology of the unstructured grid,
//!			   and then extract the node topology information,
//!	           such as node-cell and node-face connected relationship.
//! @author    Meng Liyuan.

#include "TK_Log.h"
#include "TK_Time.h"
#include "Glb_Dimension.h"
#include "Geo_NodeTopo_Unstruct.h"
#include "Region.h"
#include "PHIO.h"
#include "GlobalDataBase.h"
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
    ParallelOpenFile(ofile, "./out/UnsNodeTopologyOutput.txt", ios::out);

    //! Read in the unstructured grid.
    //! Update all of parameters of each zone.
    //! Calculate the grid topology by InitGeometry function.
    Region *region = new Region();
    region->ReadGrid();
    region->UpdateAllZonesPara();
    region->InitGeometry();

    //! Extract and dump out the node topology information of the unstructured grid.
    int nZones = GetNumberofGlobalZones();
    for (int iGrid = 0; iGrid < nZones; ++ iGrid)	
    {
        ofile << "Unstructured Grid NodeTopology of Zones "<< iGrid << " is: " << endl;

        Grid * grid = GetGrid(iGrid);
        int  * cell_number_of_each_node = UnstructGridCast(grid)->GetCellNumberOfEachNode();
        int  * node2cell = UnstructGridCast(grid)->GetNode2Cell();

        int nTotalNode = grid->GetNTotalNode();
        for (int iNode = 0; iNode < nTotalNode; ++ iNode)
        {
            ofile << "   The Cells Number of Node "<< iNode << " is:" << cell_number_of_each_node[iNode] << endl;
            for (int icell = 0; icell < cell_number_of_each_node[iNode]; ++ icell)
            { 
                ofile << "   The Connected Cells Index of Node "<< iNode << " is:" << node2cell[icell] << endl;
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