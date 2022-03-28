//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      WalldistCompute.cpp
//! @brief     It tests the function of wall distance computation for both structured and unstructured grid.
//! @author    Meng Liyuan.

#include "Pre_WalldistCompute.h"
#include "GlobalDataBase.h"
#include "Region.h"
#include "PHIO.h"
#include "TK_Time.h"
#include "Constants.h"
#include "WallDist.h"
#include "H5Fpublic.h"
#include "H5Ipublic.h"
#include "IO_HDF5File.h"
#include "iomanip"

using namespace PHSPACE;
using namespace PHMPI;

template<typename T>
string toString(T t)
{
    stringstream ss;
    ss.setf(ios::fixed);
    ss<<setprecision(6)<<t;
    return ss.str();
}

int main()
{
    //! Record time.
    setvbuf( stdout, NULL, _IONBF, 0 );
    double *ptr_time = NULL;
    PHSPACE::TIME_SPACE::GetWallTimeForResidual();
    PHSPACE::TIME_SPACE::GetWallTime();
    PHSPACE::TIME_SPACE::GetCPUTime(ptr_time);

    //! Read the parameters file and perform the initialization operation of PHengLEI.
    ReadControlParameters();
    CompleteParamInfo();
    Initialization();

    //! Read in the grid .
    //! Update all of parameters of each zone.
    Region *region = new Region();
    region->ReadGrid();
    region->UpdateAllZonesPara();	
    region->InitGrids();

    AllocateWalldist();

    //! Compute Wall distance for the grid.
    int walldistMethod;
    GlobalDataBase::GetData("walldistMethod", &walldistMethod, PHINT, 1);

    int level = 0;
    Pre_WalldistCompute WalldistCompute(walldistMethod, level);
    WalldistCompute.Run();
    WalldistCompute.ComputeNodeDistance();

    WriteWallDist();


    //! read sphere_0.wdt file ,and put it partly data into WalldistComputeOutput.txt
    hid_t file = H5Fopen("./out/sphere_0.wdt", H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0)
    {
        return 0;
    }
    ostringstream oss1;
    oss1<<"Grid-"<<0;
    string dataName = oss1.str();
    //! get walldist array.
    Grid *grid = GetGrid(0, 0);
    grid->InitWallDist();
    int nTotalCell = grid->GetNTotalCell();
    RDouble *walldist = new RDouble[nTotalCell];
    ReadData(file, walldist, dataName);
    grid->ReadWalldist(walldist);
    delete []walldist;
    
    fstream ofile;
    ParallelOpenFile(ofile, "./out/WalldistComputeOutput.txt", ios::out);
    double *walldist_sphere = UnstructGridCast( grid )->GetWallDist();
    std::ostringstream oss;
    int i=0;
    while (i<60)
    {
        oss <<setiosflags(ios::fixed);
        oss <<setprecision(10)<<  " the "<< i << "th  = " <<walldist_sphere[i] << endl;
        i++;
    }
    string str = oss.str();
    ofile << str << endl;; 
    ParallelCloseFile(ofile);

    //! Exit PHengLEI.
    PHMPI::FreeBlockData();
    TIME_SPACE::WriStatTime();
    PH_Finalize();
    return 0; 
}

