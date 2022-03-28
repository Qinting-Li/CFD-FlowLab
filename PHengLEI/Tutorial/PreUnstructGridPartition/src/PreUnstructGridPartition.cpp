//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      PreUnstructGridPartition.cpp
//! @brief     It tests the function of unstructured grid patition,
//!			   check whether it can run as required or not.
//! @author    Meng Liyuan.

#include "TK_Log.h"
#include "Glb_Dimension.h"
#include "Pre_UnstructGridPartition.h"
#include "Region.h"
#include "PHIO.h"
#include "GlobalDataBase.h"
#include "TK_Time.h"
#include "Constants.h"
#include "H5Fpublic.h"
#include "IO_HDF5File.h"

using namespace PHSPACE;
using namespace PHMPI;
int main()
{
    //! Perform the initialization operation of PHengLEI.
    ReadControlParameters();
    Initialization();

    //! Define output file which to store the results.
    fstream ofile;
    ParallelOpenFile(ofile, "./out/UnstructGridPartitionOutput.txt", ios::out);

    //! Define serial parameters using in the process of unstructured grid partition.
    //! Obtain the data from global database.
    string original_grid_file;
    GlobalDataBase::GetData("original_grid_file", &original_grid_file, PHSTRING, 1);
    string partition_grid_file;
    GlobalDataBase::GetData("partition_grid_file", &partition_grid_file, PHSTRING, 1);

    int maxproc;
    GlobalDataBase::GetData("maxproc", &maxproc, PHINT, 1);
    int npartmethod;
    GlobalDataBase::GetData("npartmethod", &npartmethod, PHINT, 1);
    int parallelPartitionMethod;
    GlobalDataBase::GetData("parallelPartitionMethod", &parallelPartitionMethod, PHINT, 1);

    //! Read in the original grid and partition the unstructured grid serially. 
    Region *region = new Region();
    region->ReadOrdinaryGrid(original_grid_file);

    Grid * grid = GetGrid(0, 0);
    Grid ** parititionGrids = 0; 

    //! Partition the unstructured Grid.
    Pre_UnstructGridPartition unstructGridParition(UnstructGridCast(grid), maxproc, npartmethod, original_grid_file, parallelPartitionMethod);
    unstructGridParition.PartitionGridSerial(partition_grid_file);	
    unstructGridParition.DumpPartitionGrids(partition_grid_file);
    parititionGrids = unstructGridParition.GetPartitionGrids();

    //! output fts file's content。
    hid_t file;
    hid_t grpData;
    file = H5Fopen("./out/rae2822_hybrid2d__4_0.fts", H5F_ACC_RDONLY, H5P_DEFAULT);

    //! output ReadVersionInfo
    double version;
    grpData = OpenGroup(file, "Version");;
    ReadData(grpData, &version, "VersionData");
    ofile << "   VersionData =  " << version<< endl;

    hid_t  grpData1,grpData2;
    for (int iGrid = 0; iGrid < 4; ++ iGrid)	
    {
        //! output ReadBCRegion
        string grpName;
        ostringstream oss;
        int nBCRegion = 0;
        oss<<"Grid-"<<iGrid;
        grpName = oss.str();
        ofile <<grpName<<" : "<<endl;
        grpData1 = OpenGroup(file, grpName);//grpName表示不同的grid

        //! output iDimensions
        int *iDimensions = 0 ;
        int iDimensions_total = 0;
        iDimensions = ReadIntOneRow(grpData1, "iDimensions");
        int sss1 = iDimensions[0];
        int sss2 = iDimensions[1];
        int sss3 = iDimensions[2];
        iDimensions_total = iDimensions[0]*iDimensions[1]*iDimensions[2];
        ofile <<" iDimensions["<<0<<"] = "<<iDimensions[0]<<endl;
        ofile <<" iDimensions["<<1<<"] = "<<iDimensions[1]<<endl;
        ofile <<" iDimensions["<<2<<"] = "<<iDimensions[2]<<endl<<endl;

        //! output CoordinateX.CoordinateY,CoordinateZ
        hid_t  grpData3;
        grpData3 = OpenGroup(grpData1, "GridCoordinates");
        RDouble *CoordinateX = ReadDoubleOneRow(grpData3, "CoordinateX");
        RDouble *CoordinateY = ReadDoubleOneRow(grpData3, "CoordinateY"); 
        RDouble *CoordinateZ = ReadDoubleOneRow(grpData3, "CoordinateZ");
        for (int i=0;i<iDimensions[0];i++)
        {
            ofile <<"  CoordinateX["<<i<<"] = "<<CoordinateX[i]<<endl;
            ofile <<"  CoordinateY["<<i<<"] = "<<CoordinateY[i]<<endl;
            ofile <<"  CoordinateZ["<<i<<"] = "<<CoordinateZ[i]<<endl;
        }

        //! output interFace2ZoneID,interFace2InterFaceID,interFace2BoundaryFace
        int nIFace = 0;
        hid_t grpData4 = OpenGroup(grpData1, "InterfaceInfo");
        ReadData(grpData4, &nIFace, "nIFace");
        int *interFace2ZoneIDTemp       = ReadIntOneRow(grpData4, "interFace2ZoneID");
        int *interFace2InterFaceIDTemp  = ReadIntOneRow(grpData4, "interFace2InterFaceID");
        int *interFace2BoundaryFaceTemp = ReadIntOneRow(grpData4, "interFace2BoundaryFace");
        for (int i=0;i<nIFace;i++)
        {
            ofile <<"  interFace2ZoneIDTemp      ["<<i<<"] = "<<interFace2ZoneIDTemp[i]<<endl;
            ofile <<"  interFace2InterFaceIDTemp ["<<i<<"] = "<<interFace2InterFaceIDTemp[i]<<endl;
            ofile <<"  interFace2BoundaryFaceTemp["<<i<<"] = "<<interFace2BoundaryFaceTemp[i]<<endl;
        }

        //! output BCName
        hid_t grpData5 = OpenGroup(grpData1, "BCName");
        char *BCName =   ReadCharData(grpData5, "BCName");
        for (int i=0;i<500;i++)
        {
            if(BCName[i]<0)  
            {
                break;    
            }
            ofile <<"  BCName["<<i<<"]       = "<<BCName[i]<<endl;
        }
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
        ofile << "  blockProc_local["<<i<<"] = " <<  blockProc_local[i] <<endl;
        ofile << "  blockIdx_local ["<<i<<"] = " <<  blockIdx_local [i] <<endl;
        ofile << "  blockType_local["<<i<<"] = " <<  blockType_local[i] <<endl;
    }

    H5Gclose(grpData1);

    //! Close the output file and exit PHengLEI.
    ParallelCloseFile(ofile);

    PHMPI::FreeBlockData();
    TIME_SPACE::WriStatTime();
    PH_Finalize();
    return 0;
}