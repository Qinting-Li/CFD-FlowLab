//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      PreStructGridPartition.cpp
//! @brief     It tests the function of structured grid patition,
//!			   check whether it can run as required or not.
//! @author    Meng Liyuan.

#include "TK_Log.h"
#include "Glb_Dimension.h"
#include "Pre_StrGridPartition.h"
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
    ParallelOpenFile(ofile, "./out/StructGridPartitionOutput.txt", ios::out);

    //! Define parameters using in the process of structured grid partition.
    //! Obtain the data from global database.
    int numberOfProcessors = PHSPACE::GlobalDataBase::GetIntParaFromDB( "maxproc" );

    //! Partition the Structured Grid.
    Pre_StrGridPartition StructGridPartition(numberOfProcessors);

    StructGridPartition.Run();

    //! output fts file's content¡£
    hid_t file;
    file = H5Fopen("./out/m6__4_0.fts", H5F_ACC_RDONLY, H5P_DEFAULT);

    hid_t  grpData;
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
        grpData1 = OpenGroup(file, grpName);

        int *iDimensions  = ReadIntOneRow(grpData1, "iDimensions");
        int iDimensions_total = iDimensions[0]*iDimensions[1]*iDimensions[2];
        ofile <<" iDimensions["<<0<<"] = "<<iDimensions[0]<<endl;
        ofile <<" iDimensions["<<1<<"] = "<<iDimensions[1]<<endl;
        ofile <<" iDimensions["<<2<<"] = "<<iDimensions[2]<<endl<<endl;

        hid_t  grpData3 = OpenGroup(grpData1, "GridCoordinates");
        RDouble *CoordinateX = ReadDoubleOneRow(grpData3, "CoordinateX");
        RDouble *CoordinateY = ReadDoubleOneRow(grpData3, "CoordinateY"); 
        RDouble *CoordinateZ = ReadDoubleOneRow(grpData3, "CoordinateZ");
        for (int i=0;i<=iDimensions_total;i+=100)
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
        for (int i=0;i<nIFace;i+=100)
        {
            ofile <<"  interFace2ZoneIDTemp      ["<<i<<"] = "<<interFace2ZoneIDTemp[i]<<endl;
            ofile <<"  interFace2InterFaceIDTemp ["<<i<<"] = "<<interFace2InterFaceIDTemp[i]<<endl;
            ofile <<"  interFace2BoundaryFaceTemp["<<i<<"] = "<<interFace2BoundaryFaceTemp[i]<<endl;
        }

        hid_t grpData5 = OpenGroup(grpData1, "BCName");
        char  *BCName  = ReadCharData(grpData5, "BCName");
        for (int i=0;i<61;i+=4)
        {
            if( (BCName[i]<0) )  
            {
                break;    
            }
            ofile <<"  BCName["<<i<<"]       = "<<(int)BCName[i]<<endl;
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