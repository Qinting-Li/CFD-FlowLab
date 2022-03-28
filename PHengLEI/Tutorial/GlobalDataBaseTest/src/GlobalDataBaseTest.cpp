//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      GlobalDataBaseTest.cpp
//! @brief     It tests the functions of global database , 
//!            such as data acquisition and updating ,  
//!			   check whether they can run as required or not.
//! @author    Meng Liyuan.

#include "TK_Log.h"
#include "Region.h"

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
    ParallelOpenFile(ofile, "./out/GlobalDataBaseOutput.txt", ios::out);

    //! Obtain the data from the global database.
    string gridfile;
    GlobalDataBase::GetData("gridfile", &gridfile,PHSTRING, 1);
    string uns_gradient_name;
    GlobalDataBase::GetData("uns_gradient_name", &uns_gradient_name, PHSTRING,1);
    ofile << "The gridfile is : "         << gridfile << endl;
    ofile << "The uns_gradient_name is : "<< uns_gradient_name << endl;

    int mgrid ;
    GlobalDataBase::GetData("mgrid", &mgrid,PHINT,1);
    int iale ;
    GlobalDataBase::GetData("iale", &iale,PHINT, 1);
    int iviscous  ;
    GlobalDataBase::GetData("iviscous", &iviscous, PHINT, 1);
    ofile << "The mgrid is : "    << mgrid << endl;
    ofile << "The iale is : "     << iale << endl;
    ofile << "The iviscous is : " << iviscous << endl;

    string uns_limiter_name = GlobalDataBase::GetStrParaFromDB( "uns_limiter_name");
    int iovrlap = GlobalDataBase::GetIntParaFromDB( "iovrlap");
    double reference_mach_number = GlobalDataBase::GetDoubleParaFromDB("reference_mach_number");
    ofile << "The uns_limiter_name is : "     << uns_limiter_name << endl;
    ofile << "The iovrlap is : "              << iovrlap << endl;
    ofile << "The reference_mach_number is : "<< reference_mach_number << endl;

    //! Update the data of the global database.
    int iunsteady = 1;
    GlobalDataBase::UpdateData("iunsteady", &iunsteady,PHINT, 1);
    int read_wall_dist = 1;
    GlobalDataBase::UpdateDataPtr("read_wall_dist", &read_wall_dist);
    ofile << "The updated iunsteady is : "     << iunsteady << endl;
    ofile << "The updated read_wall_dist is : "<< read_wall_dist << endl;

    //! Close the output file and exit PHengLEI.
    ParallelCloseFile(ofile);

    PHMPI::FreeBlockData();
    TIME_SPACE::WriStatTime();
    PH_Finalize();
    return 0;
}