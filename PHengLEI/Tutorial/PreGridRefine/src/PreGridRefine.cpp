//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      PreGridRefine.cpp
//! @brief     It tests the function of grid Refining,
//!			   check whether it can run as required or not.
//! @author    Meng Liyuan.

#include "TK_Log.h"
#include "Glb_Dimension.h"
#include "Pre_GridRefine.h"
#include "PHIO.h"
#include "GlobalDataBase.h"
#include "TK_Time.h"
#include "H5Fpublic.h"
#include "IO_HDF5File.h"

using namespace PHSPACE;
using namespace PHMPI;
int main()
{
    //! Recording time.
    setvbuf( stdout, NULL, _IONBF, 0 );
    double *ptr_time = NULL;
    PHSPACE::TIME_SPACE::GetWallTimeForResidual();
    PHSPACE::TIME_SPACE::GetWallTime();
    PHSPACE::TIME_SPACE::GetCPUTime(ptr_time);

    //! Perform the initialization operation of PHengLEI.
    ReadControlParameters();
    Initialization();	

    fstream ofile;
    ParallelOpenFile(ofile, "./out/PreGridRefineOutput.txt", ios::out);

    //! Define parameters using in the process of grid refining.
    //! Obtain the data from global database.
    string from_gfile;
    GlobalDataBase::GetData("from_gfile", &from_gfile,PHSTRING,1);
    Pre_GridRefine GridRefine( from_gfile);
    GridRefine.Run();

    //! output fts file content。
    hid_t file, grpData;
    file = H5Fopen("./out/m6Refine_0.fts", H5F_ACC_RDONLY, H5P_DEFAULT);

    //! output ReadVersionInfo
    double version;
    grpData = OpenGroup(file, "Version");;
    ReadData(grpData, &version, "VersionData");
    ofile << "   VersionData =  " << version<< endl;
    hid_t  grpData1,grpData2;
    int iGrid = 0;

    //! output ReadBCRegion
    string grpName;
    ostringstream oss;
    oss<<"Grid-"<<iGrid;
    grpName = oss.str();
    ofile <<grpName<<" : "<<endl;
    grpData1 = OpenGroup(file, grpName);//grpName表示不同的grid

    int *iDimensions = 0 ;
    int iDimensions_total = 0;
    iDimensions = ReadIntOneRow(grpData1, "iDimensions");
    iDimensions_total = iDimensions[0]*iDimensions[1]*iDimensions[2];
    ofile <<" iDimensions["<<0<<"] = "<<iDimensions[0]<<endl;
    ofile <<" iDimensions["<<1<<"] = "<<iDimensions[1]<<endl;
    ofile <<" iDimensions["<<2<<"] = "<<iDimensions[2]<<endl<<endl;

    hid_t  grpData3;
    grpData3 = OpenGroup(grpData1, "GridCoordinates");
    RDouble *CoordinateX = ReadDoubleOneRow(grpData3, "CoordinateX");
    RDouble *CoordinateY = ReadDoubleOneRow(grpData3, "CoordinateY"); 
    RDouble *CoordinateZ = ReadDoubleOneRow(grpData3, "CoordinateZ");
    for (int i=0;i<=iDimensions_total;i+=10000)
    {
        ofile <<"  CoordinateX["<<i<<"] = "<<CoordinateX[i]<<endl;
        ofile <<"  CoordinateY["<<i<<"] = "<<CoordinateY[i]<<endl;
        ofile <<"  CoordinateZ["<<i<<"] = "<<CoordinateZ[i]<<endl;
    }

    int nIFace = 0;
    hid_t grpData4 = OpenGroup(grpData1, "InterfaceInfo");
    ReadData(grpData4, &nIFace, "nIFace");
    int *interFace2ZoneIDTemp       = ReadIntOneRow(grpData4, "interFace2ZoneID");
    int *interFace2InterFaceIDTemp  = ReadIntOneRow(grpData4, "interFace2InterFaceID");
    int *interFace2BoundaryFaceTemp = ReadIntOneRow(grpData4, "interFace2BoundaryFace");
    for (int i=0;i<nIFace;i+=500)
    {
        ofile <<"  interFace2ZoneIDTemp      ["<<i<<"] = "<<interFace2ZoneIDTemp[i]<<endl;
        ofile <<"  interFace2InterFaceIDTemp ["<<i<<"] = "<<interFace2InterFaceIDTemp[i]<<endl;
        ofile <<"  interFace2BoundaryFaceTemp["<<i<<"] = "<<interFace2BoundaryFaceTemp[i]<<endl;
    }

    hid_t grpData5 = OpenGroup(grpData1, "BCName");
    char *BCName   = ReadCharData(grpData5, "BCName");
    for (int i=0;i<64;i+=4)
    {
        ofile <<"  BCName["<<i<<"]       = "<<(int)BCName[i]<<endl;
    }

    //! output ReadTotalInfo
    hid_t  grpData0;
    grpData0 = OpenGroup(file, "Information");
    int *blockProc_local = ReadIntOneRow(grpData0, "block_proc");
    int *blockIdx_local  = ReadIntOneRow(grpData0, "block_idx");
    int *blockType_local = ReadIntOneRow(grpData0, "block_type");
    int *nBlocks_local   = ReadIntOneRow(grpData0, "nBlocks");
    for (int i=0;i<nBlocks_local[0];i++)
    {
        ofile <<"  blockProc_local["<<i<<"] = "<<  blockProc_local[i]<<endl;
        ofile <<"  blockIdx_local ["<<i<<"] = "<<  blockIdx_local [i]<<endl;
        ofile <<"  blockType_local["<<i<<"] = "<<  blockType_local[i]<<endl;
    }

    H5Gclose(grpData1);

    //! Close the output file and exit PHengLEI.
    ParallelCloseFile(ofile);

    PHMPI::FreeBlockData();
    TIME_SPACE::WriStatTime();
    PH_Finalize();
    return 0;

}