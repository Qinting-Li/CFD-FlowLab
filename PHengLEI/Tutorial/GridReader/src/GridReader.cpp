//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      GridReader.cpp
//! @brief     It tests the functions of grid reader  
//!			   then extract and dump out the grid information.
//! @author    Meng Liyuan.

#include "IO_GridReader.h"
#include "PHMpi.h"
#include "TK_Log.h"
#include "TK_Time.h"
#include "GlobalDataBase.h"
#include "PHIO.h"
#include "Constants.h"
#include "TypeDefine.h"

using namespace PHSPACE;
using namespace PHMPI;
int main()
{
    //! Perform the initialization operation of PHengLEI.
    ReadControlParameters();
    Initialization();

    //! Define output file which to store the results.
    fstream ofile;
    ParallelOpenFile(ofile, "./out/ReadOut.txt", ios::out);

    //! Read the grid.
    string gridFileName;
    GlobalDataBase::GetData("from_gfile", &gridFileName, PHSTRING, 1);
    PHVectorString1D gridGroupNames;
    gridGroupNames.push_back( gridFileName );
    IO_GridReader gridReader(gridGroupNames, 2);
    gridReader.ReadGridsByCSMode();
    gridReader.ReadCellToNodeByCSMode();

    //! Dump out the grid information.
    int nZones = gridReader.GetNumberOfGrids();

    ofile << "  Number of Global Zones : "       << gridReader.GetNumberofGlobalZones() << endl;
    ofile << "  Number of grid files : "         << gridReader.GetNumberOfGridFile()    << endl;
    ofile << "  Number of Zones in processor "   << GetCurrentProcessorID()             << ": " << nZones << endl;

    for ( int iGrid = 0; iGrid < nZones; ++ iGrid )
    {
        Grid * grid = gridReader.GetGrid(iGrid);
        double minDis, maxDis;
        grid->GetMinMaxDS(minDis, maxDis);  

        ofile << "  The grid type:  " << grid->Type()           << endl;
        ofile << "  Dim = "           << grid->GetDim()         << endl;
        ofile << "  For Grid "        << iGrid << ", minDis = " << minDis        << ", " << "maxDis = "   << maxDis << " " << endl;
        ofile << "  nTNode = "        << grid->GetNTotalNode()  << ", nTCell = " << grid->GetNTotalCell() << endl;
        ofile << "  nTFace = "        << grid->GetNTotalFace()  << ", nBFace = " << grid->GetNBoundFace() << ", nIFace = " << grid->GetNIFace() << endl;    
    }

    //! Close the output file and exit PHengLEI.
    ParallelCloseFile(ofile);

    PHMPI::FreeBlockData();
    TIME_SPACE::WriStatTime();
    PH_Finalize();
    return 0;


}
