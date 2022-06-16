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
static int CountInvokeFun =0;
int main()
{
    //! Perform the initialization operation of PHengLEI.
    ReadControlParameters();
    Initialization();

    //! Define output file which to store the results.
    fstream ofile;
    ParallelOpenFile(ofile, "./out/UT_UnstrGridMemberFunction.txt", ios::out);

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
        GetData(grid,ofile);
        SetData(grid,ofile);
    }    
    
    delete region;

    //! Close the output file and exit PHengLEI.
    ParallelCloseFile(ofile);
    PHMPI::FreeBlockData();
    TIME_SPACE::WriStatTime();
    PH_Finalize();
    return 0;
}

void GetData(Grid *grid,fstream &ofile)
{
    if ( 1 == ++CountInvokeFun)
    {
        ofile << "the first time invoke GetData() "<< endl;
    }
    //ofile <<"CountInvokeFun = " << CountInvokeFun <<endl;
             
    int n1 = UnstructGridCast(grid)->GetDim();
    int n2 = UnstructGridCast(grid)->GetLevel();
    int n3 = UnstructGridCast(grid)->GetNTotalFace();
             
    int n4 = UnstructGridCast(grid)->GetNTotalCell();
    int n5 = UnstructGridCast(grid)->GetNBoundFace();
    int n6 = UnstructGridCast(grid)->GetNIFace();
    int n7 = UnstructGridCast(grid)->GetZoneID();
             
    int n  = UnstructGridCast(grid)->Type();

    ofile << "invoke GetData() : "<< endl ;
    ofile << "invoke function GetDim();        return value is: " << n1 << endl;
    ofile << "invoke function GetLevel();      return value is: " << n2 << endl;
    ofile << "invoke function GetNTotalFace()  return value is: " << n3 << endl;
    ofile << "invoke function GetNTotalCell()  return value is: " << n4 << endl;
    ofile << "invoke function GetNBoundFace()  return value is: " << n5 << endl;
    ofile << "invoke function GetNIFace();     return value is: " << n6 << endl;
    ofile << "invoke function GetZoneID();     return value is: " << n7 << endl<< endl;

    //grid->Next();
    //grid->Prev();       
}

void SetData(Grid *grid,fstream &ofile)
{
    ofile << endl;
    ofile << "invoke SetData() "<< endl ;
    
    UnstructGridCast(grid)->SetDim(4);
    UnstructGridCast(grid)->SetLevel(5);
    UnstructGridCast(grid)->SetNTotalFace(6);
    
    UnstructGridCast(grid)->SetNTotalCell(7);
    UnstructGridCast(grid)->SetNBoundFace(8);

    InterfaceInfo *iinfo = grid->GetInterfaceInfo();
    iinfo->SetNIFace(9);

    UnstructGridCast(grid)->SetZoneID(10);

    ofile << "       SetDim(4);"<<endl;
    ofile << "       SetLevel(5);"<<endl;
    ofile << "       SetNTotalFace(6);"<<endl;
    ofile << "       SetNTotalCell(7);"<<endl;
    ofile << "       SetNBoundFace(8);"<<endl;
    ofile << "       SetNIFace(9)"<<endl;
    ofile << "       SetZoneID(10);"<<endl<< endl;

    GetData(grid,ofile);
}