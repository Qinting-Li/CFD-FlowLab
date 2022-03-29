//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      CGNSConvertToPlot3D.cpp
//! @brief     It tests the function of grid conversion  
//!            from CGNS format to plot3D.
//! @author    Meng Liyuan.

#include "TK_Log.h"
#include "Pre_CGNSToPlot3D_Struct.h"
#include "PHIO.h"
#include "GlobalDataBase.h"
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
    ParallelOpenFile(ofile, "./out/CGNSToPlot3DOutput.txt", ios::out);

    //! Read in grid and convert.
    string gridFileName;
    GlobalDataBase::GetData("from_gfile", &gridFileName, PHSTRING, 1);
    Pre_CGNSToPlot3D_Struct cgnsGridConversionToPlot3D(gridFileName);
    cgnsGridConversionToPlot3D.Run();

    //! Close the output file and exit PHengLEI.
    ParallelCloseFile(ofile);

    PHMPI::FreeBlockData();
    TIME_SPACE::WriStatTime();
    PH_Finalize();
    return 0;
}