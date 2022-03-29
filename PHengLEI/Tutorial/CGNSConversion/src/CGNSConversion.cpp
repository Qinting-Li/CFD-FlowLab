//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      CGNSConversion.cpp
//! @brief     It tests the function of grid conversion from CGNS to HyperFlow,
//!			   then extract and dump out the grid information.
//! @author    Meng Liyuan.

#include "TK_Log.h"
#include "TK_Time.h"
#include "Glb_Dimension.h"
#include "Pre_CGNSToPHengLEI_Unstruct.h"
#include "GlobalDataBase.h"
#include "PHIO.h"
#include "Constants.h"

using namespace PHSPACE;
using namespace PHMPI;
int main()
{
    //! Perform the initialization operation of PHengLEI.
    ReadControlParameters();
    Initialization();
    SetDim(TWO_D);

    //! Define output file which to store the results.
    fstream file;
    ParallelOpenFile(file, "./out/circleOutput.txt", ios::out);

    //! Read in grid(CGNS) and convert .
    string cgnsGridFileName;
    GlobalDataBase::GetData("from_gfile", &cgnsGridFileName, PHSTRING, 1);
    string gridFileName;
    GlobalDataBase::GetData("out_gfile", &gridFileName, PHSTRING, 1);
    Pre_CGNSToPHengLEI_Unstruct cgnsGridConversion(cgnsGridFileName);
    cgnsGridConversion.Run();

    //! Extract and dump out the converted grid information  .
    int nBlocks = cgnsGridConversion.GetNumberofBlocks();
    file << "    nBlock = " << nBlocks << endl;	
    Grid ** grids = cgnsGridConversion.GetGrids();
    for ( int iBlock = 0; iBlock < nBlocks; ++ iBlock )
    {
        Grid * grid = &**grids;
        file << "    The grid type:  "             << grid->Type()           << endl;
        file << "    Dim = "                       << grid->GetDim()         << endl;
        file << "    nTNode = "                    << grid->GetNTotalNode()  << endl;
        file << "    nTCell = "                    << grid->GetNTotalCell()  << endl;
        file << "    nTFace = "                    << grid->GetNTotalFace()  << ",   nBFace = " << grid->GetNBoundFace() <<",   nIFace = " << grid->GetNIFace() << endl;	
        file << "    The Global Zone Grid index: " << grid->GetZoneID()      << endl;	
        file << "    The  Local Zone Grid index: " << grid->GetZoneLocalID() << endl;	

        UnstructGrid * unstrGrid = (UnstructGrid *) grid;
        double minDis, maxDis;
        unstrGrid->GetMinMaxDS(minDis, maxDis); 
        file<< "    For Grid " << iBlock << ",   minDis = " << minDis << ",   maxDis = " << maxDis << " " << endl;

        InterfaceInfo    *iinfo;
        InterfaceFields  *interfaceFields;
        OversetInfoProxy *GetOversetInfoProxy;
        GridID           *GetGridID;
        SimpleGrid       *GetOldGrid;
        iinfo               = grid->GetInterfaceInfo();
        interfaceFields     = grid->GetInterfaceFields();
        GetOversetInfoProxy = grid->GetOversetInfoProxy();
        GetGridID           = grid->GetGridID();
        GetOldGrid          = grid->GetOldGrid();

        file << "    Interface topology information : "               << iinfo                  << endl;
        file << "    Interface fields information:  "                 << interfaceFields        << endl;
        file << "    The Over-set or overlapping Grid information:  " << GetOversetInfoProxy    << endl;
        file << "    The Grid index : "                               << GetGridID              << endl;
        file << "    The Grid at the last time step : "               << GetOldGrid             << endl;
        file << "    The Coarser grid of the current grid: "          << grid->GetCoarseGrid()  << endl;
        file << "    Is this grid is the coarsest grid: "             << grid->IsCoarsestGrid() << endl;
        file << "    The Coarser level grid: "                        << grid->Next()           << endl;
        file << "    The Finer grid of the current grid:  "           << grid->GetFineGrid()    << endl;
        file << "    Is this grid is the finest grid: "               << grid->IsFinestGrid()   << endl;
        file << "    The Finer level grid: "                          << grid->Prev()           << endl;
    }	

    DumpOldGrid(gridFileName, grids, nBlocks);

    //! Clear the memory, and delete the converted computational grid. 
    cgnsGridConversion.Clear();
    FreeBlockData();

    //! Close the output file and exit PHengLEI.
    ParallelCloseFile(file);

    PHMPI::FreeBlockData();
    TIME_SPACE::WriStatTime();
    PH_Finalize();
    return 0;
}