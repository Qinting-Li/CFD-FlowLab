//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      FieldViewGridConversion.cpp
//! @brief     It tests the function of grid conversion  
//!            from field view format to HyperFlow,
//!			   then extract and dump out the grid information.
//! @author    Meng Liyuan.

#include "TK_Log.h"
#include "Glb_Dimension.h"
#include "Pre_FieldViewToPHengLEI_Unstruct.h"
#include "GlobalDataBase.h"
#include "PHIO.h"
#include "Constants.h"
#include "TK_Time.h"

using namespace PHSPACE;
using namespace PHMPI;
int main()
{
    //! Perform the initialization operation of PHengLEI.
    ReadControlParameters();
    Initialization();

    //! Define output file which to store the results.
    fstream ofile;
    ParallelOpenFile(ofile, "./out/FieldView_GridConversionOutput.txt", ios::out);

    //! Read in grid and convert.
    string fileName;
    GlobalDataBase::GetData("from_gfile", &fileName, PHSTRING, 1);
    string gridFileName;
    GlobalDataBase::GetData("out_gfile", &gridFileName, PHSTRING, 1);

    Pre_FieldViewToPHengLEI_Unstruct fieldViewGridConversion(fileName);
    fieldViewGridConversion.Run();

    //! Extract and dump out the grid information.
    int nBlocks = fieldViewGridConversion.GetNumberofBlocks();
    ofile << "    Number of Blocks: " << nBlocks << endl;

    Grid ** grids = fieldViewGridConversion.GetGrids();
    for ( int iBlock = 0; iBlock < nBlocks; ++ iBlock )
    {
        Grid * grid = grids[iBlock];
        double minDis, maxDis;
        grid->GetMinMaxDS(minDis, maxDis);

        ofile << "    For Block "       << iBlock << ", minDis = " << minDis << ", " << "maxDis = " << maxDis << " " << endl;
        ofile << "    The grid type:  " << grid->Type() << endl;
        ofile << "    Dim = "           << grid->GetDim() << endl;
        ofile << "    nTNode = "        << grid->GetNTotalNode()   << endl;
        ofile << "    nTCell = "        << grid->GetNTotalCell()   << endl;
        ofile << "    nTFace = "        << grid->GetNTotalFace()   << ",   nBFace = " << grid->GetNBoundFace() <<",   nIFace = " << grid->GetNIFace() << endl;	
    }

    DumpGrid(gridFileName, grids, nBlocks);
    fieldViewGridConversion.WriteAdditionalInformation(gridFileName, grids);

    //! Clear the memory, and delete the converted computational grid. 
    fieldViewGridConversion.Clear(); 

    //! Close the output file and exit PHengLEI.
    ParallelCloseFile(ofile);

    PHMPI::FreeBlockData();
    TIME_SPACE::WriStatTime();
    PH_Finalize();
    return 0; 
}