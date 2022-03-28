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

void GetData(Grid *grid,fstream &ofile);
void SetData(Grid *grid,fstream &ofile);

int main()
{
    //! Perform the initialization operation of PHengLEI.
    ReadControlParameters();
    Initialization();

    //! Define output file which to store the results.
    fstream ofile;
    ParallelOpenFile(ofile, "./out/UT_GridMemberFunction.txt", ios::out);

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
        //int ni = StructGridCast(grid)->GetNI();
        //int nj = StructGridCast(grid)->GetNJ();
        //int nk = StructGridCast(grid)->GetNK();
        //int n1 = StructGridCast(grid)->GetDim();
        //int n2 = StructGridCast(grid)->GetLevel();
        //int n3 = StructGridCast(grid)->GetNTotalFace();
        //int n4 = StructGridCast(grid)->GetNTotalCell();
        //int n5 = StructGridCast(grid)->GetNBoundFace();
        //int n6 = StructGridCast(grid)->GetNIFace();
        //int n7 = StructGridCast(grid)->GetZoneID();
        //int n8 = StructGridCast(grid)->GetZoneLocalID();
        //int n = StructGridCast(grid)->Type();
        //ofile << "invoke function GetDim();        return value is: " << n1 << endl;
        //ofile << "invoke function GetLevel();      return value is: " << n2 << endl;
        //ofile << "invoke function GetNTotalFace()  return value is: " << n3 << endl;
        //ofile << "invoke function GetNTotalCell()  return value is: " << n4 << endl;
        //ofile << "invoke function GetNBoundFace()  return value is: " << n5 << endl;
        //ofile << "invoke function GetNIFace();     return value is: " << n6 << endl;
        //ofile << "invoke function GetZoneID();     return value is: " << n7 << endl;
        GetData(grid,ofile);
        SetData(grid,ofile);       
    }    
    
    //! Close the output file and exit PHengLEI.
    ParallelCloseFile(ofile);
    PHMPI::FreeBlockData();
    TIME_SPACE::WriStatTime();
    PH_Finalize();
    return 0;
}

void GetData(Grid *grid,fstream &ofile)
{
    ofile << "the first time invoke GetData() "<< endl;
    int ni = StructGridCast(grid)->GetNI();
    int nj = StructGridCast(grid)->GetNJ();
    int nk = StructGridCast(grid)->GetNK();

    int n1 = StructGridCast(grid)->GetDim();
    int n2 = StructGridCast(grid)->GetLevel();
    int n3 = StructGridCast(grid)->GetNTotalFace();
    int n4 = StructGridCast(grid)->GetNTotalCell();
    int n5 = StructGridCast(grid)->GetNBoundFace();
    int n6 = StructGridCast(grid)->GetNIFace();
    int n7 = StructGridCast(grid)->GetZoneID();
    int n8 = StructGridCast(grid)->GetZoneLocalID();

    int n  = StructGridCast(grid)->Type();

    ofile << "invoke GetData() "<< endl << endl;
    ofile << "invoke function GetDim();        return value is: " << n1 << endl;
    ofile << "invoke function GetLevel();      return value is: " << n2 << endl;
    ofile << "invoke function GetNTotalFace()  return value is: " << n3 << endl;
    ofile << "invoke function GetNTotalCell()  return value is: " << n4 << endl;
    ofile << "invoke function GetNBoundFace()  return value is: " << n5 << endl;
    ofile << "invoke function GetNIFace();     return value is: " << n6 << endl;
    ofile << "invoke function GetZoneID();     return value is: " << n7 << endl;

    //grid->Next();
    //grid->Prev();       
}

void SetData(Grid *grid,fstream &ofile)
{
    ofile << "invoke SetData() "<< endl << endl;
    StructGridCast(grid)->SetNI(50);
    StructGridCast(grid)->SetNJ(51);
    StructGridCast(grid)->SetNK(52);

    StructGridCast(grid)->SetDim(4);
    StructGridCast(grid)->SetLevel(5);
    StructGridCast(grid)->SetNTotalFace(6);
    StructGridCast(grid)->SetNTotalCell(7);
    StructGridCast(grid)->SetNBoundFace(8);
    StructGridCast(grid)->SetNIFace(9);
    StructGridCast(grid)->SetZoneID(10);
    StructGridCast(grid)->SetZoneLocalID(11);

    ofile << "invoke GetData() again"<< endl << endl;
    GetData(grid,ofile);
}