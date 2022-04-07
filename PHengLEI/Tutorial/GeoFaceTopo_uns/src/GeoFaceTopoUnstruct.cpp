//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      GeoFaceTopoUnstruct.cpp
//! @brief     It calculates the face topology of the unstructured grid,
//!			   and then extract the face topology information,
//!	           such as face-cell and face-node connected relationship.
//! @author    Meng Liyuan.

#include "TK_Log.h"
#include "TK_Exit.h"
#include "TK_Time.h"
#include "Glb_Dimension.h"
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
    ParallelOpenFile(ofile, "./out/UnsFaceTopologyOutput.txt", ios::out);

    //! Read in the unstructured grid.
    //! Update all of parameters of each zone.
    //! Calculate the grid topology by InitGeometry function.
    Region *region = new Region();
    region->ReadGrid();
    region->UpdateAllZonesPara();
    region->InitGeometry();

    //! Extract and dump out the face topology information of the unstructured grid.
    int nZones = GetNumberofGlobalZones();
    for (int iGrid = 0; iGrid < nZones; ++ iGrid)	
    {
        ofile << "Unstructured Grid Face Topology of Zone "<< iGrid << " is: " << endl;

        Grid * grid = GetGrid(iGrid);
        int *left_cell_of_face = UnstructGridCast(grid)->GetLeftCellOfFace();
        int *right_cell_of_face = UnstructGridCast(grid)->GetRightCellOfFace();;
        int *node_number_of_each_face = UnstructGridCast(grid)->GetNodeNumberOfEachFace();;
        int *face2node = UnstructGridCast(grid)->GetFace2Node();  

        int nTotalCell = grid->GetNTotalCell();
        int nTotalFace = grid->GetNTotalFace();

        for (int iFace = 0; iFace < nTotalFace; ++ iFace)
        {
            ofile << "   The left cell of face "<< iFace << " is: " << left_cell_of_face[iFace] << endl;
            ofile << "   The right cell of face "<< iFace << " is: " << right_cell_of_face[iFace] << endl;
            ofile << "   The node number of face "<< iFace << " is: " << node_number_of_each_face[iFace] << endl;
            for( int inode= 0; inode < node_number_of_each_face[iFace]; ++ inode )
            {
                ofile << "   The  node index of face "<< iFace << " is: " << face2node[inode] << endl;
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