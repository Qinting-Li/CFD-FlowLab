//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      GeoCellTopoUnstruct.cpp
//! @brief     It calculates the cell topology of the unstructured grid,
//!			   and then extract the cell topology information,
//!	           such as ell-face and cell-node connected relationship.
//! @author    Meng Liyuan.

#include "TK_Log.h"
#include "Glb_Dimension.h"
#include "TK_Time.h"
#include "Geo_CellTopo_Unstruct.h"
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
    ParallelOpenFile(ofile, "./out/UnsCellTopologyOutput.txt", ios::out);

    //! Read in the unstructured grid.
    //! Update all of parameters of each zone.
    //! Calculate the grid topology by InitGeometry function.
    Region *region = new Region();
    region->ReadGrid();
    region->UpdateAllZonesPara();
    region->InitGeometry();
    
    //! Extract and dump out the cell topology information of the unstructured grid.
    int nZones = GetNumberofGlobalZones();
    for (int iGrid = 0; iGrid < nZones; ++ iGrid)	
    {
        ofile << "Unstructured Grid Cell Topology of Zone "<< iGrid << " is: " << endl;

        Grid * grid = GetGrid(iGrid);
        int * face_number_of_each_cell = UnstructGridCast(grid)->GetFaceNumberOfEachCell();
        //int * node_number_of_each_cell = UnstructGridCast(grid)->GetNodeNumberOfEachCell();
        int * number_of_neighbor_cell = UnstructGridCast(grid)->GetNumberOfNeighborCell();  
        // int ** cell2face = UnstructGridCast(grid)->GetCell2Face();
        //int * cell2node = UnstructGridCast(grid)->GetCell2Node();
        vector<int> * cell2cell = UnstructGridCast(grid)->GetCell2Cell();

        int nTotalCell = grid->GetNTotalCell();

        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            ofile << "   The face number of cell "<< iCell << " is: " << face_number_of_each_cell[iCell] << endl;
            // ofile << "   The Node Number of "<< iCell << "is:" << node_number_of_each_cell[iCell] << endl;
            // ofile << "   The Node Index of "<< iCell << "is:" << cell2node[iCell] << endl;
            ofile << "   The neighbor cell number of "<< iCell << " is: " << number_of_neighbor_cell[iCell] << endl;
            /* for( int iface= 0; iface < face_number_of_each_cell[iCell]; ++ iface )
            {
            ofile << "   The face index of "<< iCell << " is: " << cell2face[iface] << endl;
            }

            for( int inode= 0; inode < node_number_of_each_cell[iCell]; ++ inode )
            {
            ofile << "   The  Nodes Index of "<< iCell << "is:" << cell2node[inode] << endl;
            }*/

            for( int icell= 0; icell < number_of_neighbor_cell[iCell]; ++ icell )
            {
                ofile << "   The neighbor cells index of "<< iCell << " is: " << cell2cell[iCell][icell] << endl;
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
