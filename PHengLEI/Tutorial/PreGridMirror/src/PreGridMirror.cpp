//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      PreGridMirror.cpp
//! @brief     It tests the function of grid mirror,
//!			   check whether it can run as required or not.
//! @author    Meng Liyuan.

#include "TK_Log.h"
#include "Glb_Dimension.h"
#include "Pre_GridMirror.h"
#include "PHIO.h"
#include "GlobalDataBase.h"
#include "TK_Time.h"
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
    ParallelOpenFile(ofile, "./out/GridMirrorOutput.txt", ios::out);

    //! Define parameters using in the process of grid mirror.
    //! Obtain the data from global database.
    string from_gfile;
    GlobalDataBase::GetData("from_gfile", &from_gfile, PHSTRING, 1);

    //! Mirror the Grid.
    Pre_GridMirror Pre_GridMirror(from_gfile);
    Pre_GridMirror.Run();

    //! Close the output file and exit PHengLEI.
    ParallelCloseFile(ofile);

    PHMPI::FreeBlockData();
    TIME_SPACE::WriStatTime();
    PH_Finalize();
    return 0;

}