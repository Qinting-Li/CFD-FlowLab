//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      FluentToPHengLEI.cpp
//! @brief     It tests the function of grid conversion  
//!            from fluent format(.cas) to PHengLEI(.fts),
//!			   then extract and dump out the grid information.
//! @author    Meng Liyuan.

#include "TK_Log.h"
#include "Glb_Dimension.h"
#include "Pre_FluentToPHengLEI_Unstruct.h"
#include "GlobalDataBase.h"
#include "PHIO.h"
#include "TK_Time.h"
#include "Constants.h"

using namespace PHSPACE;
using namespace PHMPI;
int main()
{
    //! Record time.
    setvbuf( stdout, NULL, _IONBF, 0 );
    double *ptr_time = NULL;
    PHSPACE::TIME_SPACE::GetWallTimeForResidual();
    PHSPACE::TIME_SPACE::GetWallTime();
    PHSPACE::TIME_SPACE::GetCPUTime(ptr_time);

    //! Perform the initialization operation of PHengLEI.
    ReadControlParameters();
    Initialization();

    //! Define output file which to store the results.
    fstream ofile;
    ParallelOpenFile(ofile, "./out/FluentToPHengLEIOutput.txt", ios::out);

    //! Read in grid and convert.
    string fileName;
    GlobalDataBase::GetData("from_gfile",&fileName,PHSTRING,1);
    string gridFileName;
    GlobalDataBase::GetData("out_gfile",&gridFileName,PHSTRING,1);

    Pre_FluentToPHengLEI_Unstruct FluentToPHengLEI(fileName);
    FluentToPHengLEI.Run();

    //! Extract and dump out the grid information.
    int nBlocks = FluentToPHengLEI.GetNumberofBlocks();

    ofile << "    Number of Blocks: " << nBlocks << endl;

    Grid ** grids = FluentToPHengLEI.GetGrids();
    for ( int iBlock = 0; iBlock < nBlocks; ++ iBlock )
    {
        Grid * grid = grids[iBlock];

        double minDis, maxDis;
        grid->GetMinMaxDS(minDis, maxDis);

        ofile << "    For Block " << iBlock << ", minDis = " << minDis << ", " << "maxDis = " << maxDis << " " << endl;
        ofile << "    The grid type:  " << grid->Type() << endl;
        ofile << "    Dim = " << grid->GetDim() << endl;
        ofile << "    nTNode = " << grid->GetNTotalNode() << endl;
        ofile << "    nTCell = " << grid->GetNTotalCell() << endl;
        ofile << "    nTFace = " << grid->GetNTotalFace() << ",   nBFace = " << grid->GetNBoundFace() <<",   nIFace = " << grid->GetNIFace() << endl;	

    }
    
    DumpGrid(gridFileName, grids, nBlocks);
    FluentToPHengLEI.WriteAdditionalInformation(gridFileName);

    //! Clear the memory, and delete the converted computational grid. 
    FluentToPHengLEI.Clear(); 

    //! Close the output file and exit PHengLEI.
    ParallelCloseFile(ofile);

    PHMPI::FreeBlockData();
    TIME_SPACE::WriStatTime();
    PH_Finalize();
    return 0; 
}