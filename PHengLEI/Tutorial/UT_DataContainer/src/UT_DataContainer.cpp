//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      GeoFaceTopoStruct.cpp
//! @brief     It calculates the node topology of the structured grid,
//!			   and then extract the node topology information,
//!	           such as node-cell and node-face connected relationship.
//! @author    Meng Liyuan.

#include "TK_Log.h"
#include "Glb_Dimension.h"
#include "Geo_NodeTopo_Struct.h"
#include "Region.h"
#include "PHIO.h"
#include "TK_Parse.h"
#include "GlobalDataBase.h"
#include "INCompSolver.h"
using namespace PHSPACE;
using namespace PHMPI;

int main()
{
    ReadControlParameters();
    Initialization();
    Region *region = new Region();
    region->ReadGrid();
    region->UpdateAllZonesPara();
    region->InitGeometry();

    fstream ifile,ofile;
    ParallelOpenFile(ofile, "./out/TK_DataContainer.txt", ios::out);

    string datafile = "./bin/Data.hypara";
    ifile.open(datafile.c_str(), ios_base::in);

    if ( !ifile )
    {
        cout << "could not open " << datafile << "\n";
        exit(0);
    }
    double outnstep1 = 0;
    double outnstep2 = 0;
    double outnstep3 = 0;
    double outnstep4 = 0;
    double outnstep5 = 0;
    Data_Param * testParam = new Data_Param();
    //! Read data from custom-defined parameter file. 	
    TK_Parse::ReadBasicData(ifile, testParam);
    testParam->GetData("outnstep1", &outnstep1, PHDOUBLE, 1);
    testParam->GetData("outnstep2", &outnstep2, PHDOUBLE, 1);
    testParam->GetData("outnstep3", &outnstep3, PHDOUBLE, 1);
    testParam->GetData("outnstep4", &outnstep4, PHDOUBLE, 1);

    ofile<<"Input data :" <<endl
        <<"outnstep1   "    <<outnstep1  <<endl
        <<"outnstep2   "    <<outnstep2  <<endl
        <<"outnstep3   "    <<outnstep3  <<endl
        <<"outnstep4   "    <<outnstep4  <<endl<<endl;      

    ActionKey *actkey =  new ActionKey();
    actkey->action = COMMCELLCENTERDATA;
    actkey->kind = GRID_BASED;
    actkey->level = 0;

    //fstream file;
    //OpenSerialFile( file, "./in/m6__4.fts", ios_base::in|ios_base::binary );
    //DataContainer * gridData = new DataContainer();
    //ReadFile( file, gridData );
    //DataContainer * cdata1  = gridData;
    //cdata1->moveto_begin();
    //move_pointer(cdata1,100);
    //int outnstep = 0;
    //cdata1->read(&outnstep,sizeof(int));

    Grid *grid = GetGrid(1);

    InterfaceInfo * iinfo = UnstructGridCast(grid)->GetInterfaceInfo();
    if ( iinfo == 0 ) return 0;

    int ineighbor = iinfo->FindIthNeighbor(actkey->ipos);
    int nIFaceOfNeighbor = iinfo->GetNIFaceOfNeighbor(ineighbor);
    int *faceIndexForSend = iinfo->GetFaceIndexForSend(ineighbor);

    DataContainer * cdata = actkey->GetData();
    //cdata->moveto_begin();
    //int outnstep = 0;
    //cdata->read(&outnstep,sizeof(int));

    RDouble * xcc = UnstructGridCast(grid)->GetCellCenterX();
    RDouble * ycc = UnstructGridCast(grid)->GetCellCenterY();
    RDouble * zcc = UnstructGridCast(grid)->GetCellCenterZ();
    RDouble * vol = UnstructGridCast(grid)->GetCellVolume();

    ofile <<"invoke write function!"<<endl<<endl;
    int iFace,s1;
    for ( int ifaceLocal = 0; ifaceLocal < nIFaceOfNeighbor; ++ ifaceLocal )
    {
        iFace = faceIndexForSend[ ifaceLocal ];
        UnstructGridCast(grid)->GetSourceIndex(iFace, 1, s1);

        xcc[s1]=outnstep1;
        ycc[s1]=outnstep2;
        zcc[s1]=outnstep3;
        vol[s1]=outnstep4;
       // xcc[s1]=rand();
       // ycc[s1]=rand();
       // zcc[s1]=rand();
       // vol[s1]=rand();

        cdata->write(&xcc[s1], sizeof(RDouble));
        cdata->write(&ycc[s1], sizeof(RDouble));
        cdata->write(&zcc[s1], sizeof(RDouble));
        cdata->write(&vol[s1], sizeof(RDouble));
    }

    cdata->moveto_begin();
    RDouble outnstep_read = 0;
    ofile <<"invoke read function!"<<endl<<endl;
    for (int i=0;i<nIFaceOfNeighbor;i=i+1)
    {
        cdata->read(&outnstep_read, sizeof(RDouble));
        ofile << "i = "<< i << "  outnstep_read  value is :" << outnstep_read <<endl;
    }
    ofile <<endl<<endl;
    cdata->moveto_begin();
    cdata->read(&outnstep_read, sizeof(RDouble));
    ofile << "invoke moveto_begin function ," << "  outnstep_read  value is :" << outnstep_read <<endl<<endl;
    //cdata->moveto_end();
    //cdata->read(&outnstep_read,sizeof(RDouble));
    //ofile << "moveto_end    "<< "  outnstep_read  value is :" << outnstep_read <<endl;
    cdata->moveto_begin();
    ofile <<"invoke movw_pointer function"<<endl<<endl;
    for (int i=0;i<10;i++)
    {
        move_pointer(cdata,i*sizeof(RDouble));
        cdata->read(&outnstep_read, sizeof(RDouble));
        ofile << "i = "<< i << "   outnstep_read  value is :" << outnstep_read <<endl;
    }

    ParallelCloseFile(ifile);
    ParallelCloseFile(ofile);

    delete region;
    delete testParam;
    delete actkey;

    return 0;
}

