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
    ParallelOpenFile(ofile, "./out/TK_Actionkey.txt", ios::out);

    string datafile = "./bin/Data.hypara";
    ifile.open(datafile.c_str(), ios_base::in);

    if ( !ifile )
    {
        cout << "could not open " << datafile << "\n";
        exit(0);
    }
    double outnstep1 = 0, outnstep2 = 0, outnstep3 = 0, outnstep4 = 0;
    double outnstep5 = 0, outnstep6 = 0, outnstep7 = 0, outnstep8 = 0;
    Data_Param * testParam = new Data_Param();
    //! Read data from custom-defined parameter file. 	
    TK_Parse::ReadBasicData(ifile, testParam);
    testParam->GetData("outnstep1", &outnstep1, PHDOUBLE, 1);
    testParam->GetData("outnstep2", &outnstep2, PHDOUBLE, 1);
    testParam->GetData("outnstep3", &outnstep3, PHDOUBLE, 1);
    testParam->GetData("outnstep4", &outnstep4, PHDOUBLE, 1);
    testParam->GetData("outnstep5", &outnstep5, PHDOUBLE, 1);
    testParam->GetData("outnstep6", &outnstep6, PHDOUBLE, 1);
    testParam->GetData("outnstep7", &outnstep7, PHDOUBLE, 1);
    testParam->GetData("outnstep8", &outnstep8, PHDOUBLE, 1);

    ActionKey *actkey =  new ActionKey();
    actkey->action = COMMCELLCENTERDATA;
    actkey->kind = GRID_BASED;
    actkey->level = 0;

    Grid *grid = GetGrid(1);

    InterfaceInfo * iinfo = UnstructGridCast(grid)->GetInterfaceInfo();
    if ( iinfo == 0 ) return 0;

    int ineighbor = iinfo->FindIthNeighbor(actkey->ipos);
    int nIFaceOfNeighbor = iinfo->GetNIFaceOfNeighbor(ineighbor);
    int *faceIndexForSend = iinfo->GetFaceIndexForSend(ineighbor);

    DataContainer * cdata = actkey->GetData();
    DataContainer * tecdata = actkey->GetTecData();

    RDouble * xcc = UnstructGridCast(grid)->GetCellCenterX();
    RDouble * ycc = UnstructGridCast(grid)->GetCellCenterY();
    RDouble * zcc = UnstructGridCast(grid)->GetCellCenterZ();
    RDouble * vol = UnstructGridCast(grid)->GetCellVolume();
    RDouble * tecdata1 = UnstructGridCast(grid)->GetCellCenterX();
    RDouble * tecdata2 = UnstructGridCast(grid)->GetCellCenterY();
    RDouble * tecdata3 = UnstructGridCast(grid)->GetCellCenterZ();
    RDouble * tecdata4 = UnstructGridCast(grid)->GetCellVolume();

    int iFace,s1;
    for ( int ifaceLocal = 0; ifaceLocal < nIFaceOfNeighbor; ++ ifaceLocal )
    {
        iFace = faceIndexForSend[ ifaceLocal ];
        UnstructGridCast(grid)->GetSourceIndex(iFace,1,s1);
        xcc[s1]=outnstep1;
        ycc[s1]=outnstep2;
        zcc[s1]=outnstep3;
        vol[s1]=outnstep4;

        cdata->write(&xcc[s1], sizeof(RDouble));
        cdata->write(&ycc[s1], sizeof(RDouble));
        cdata->write(&zcc[s1], sizeof(RDouble));
        cdata->write(&vol[s1], sizeof(RDouble));

    }

    RDouble outnstep_read11 = 0;
    cdata->moveto_begin();
    cdata->read(&outnstep_read11,sizeof(RDouble));
    //cout<<outnstep_read11<<endl;

    for ( int ifaceLocal = 0; ifaceLocal < nIFaceOfNeighbor; ++ ifaceLocal )
    {
        iFace = faceIndexForSend[ ifaceLocal ];
        UnstructGridCast(grid)->GetSourceIndex(iFace,1,s1);
        tecdata1[s1] = outnstep5;
        tecdata2[s1] = outnstep6;
        tecdata3[s1] = outnstep7;
        tecdata4[s1] = outnstep8;
        tecdata->write(&tecdata1[s1], sizeof(RDouble));
        tecdata->write(&tecdata2[s1], sizeof(RDouble));
        tecdata->write(&tecdata3[s1], sizeof(RDouble));
        tecdata->write(&tecdata4[s1], sizeof(RDouble));
    }

    ofile<<"Actkey invoke SetData function" <<endl;
    actkey->SetData(cdata);
    ofile<<"Actkey invoke SetTecData function" <<endl;
    actkey->SetTecData(tecdata);

    RDouble outnstep_read12 = 0;
    cdata->moveto_begin();
    cdata->read(&outnstep_read12,sizeof(RDouble));
    //cout<<outnstep_read12<<endl;

    DataContainer * cdata_data = new DataContainer();
    DataContainer * tecdata_data = new DataContainer();
    ofile<<"Actkey invoke GetData function" <<endl;
    cdata_data = actkey->GetData();
    ofile<<"Actkey invoke GetTecData function" <<endl;
    tecdata_data = actkey->GetTecData();

    int i=5;
    RDouble outnstep_read1 = 0;
    RDouble outnstep_read2 = 0;
    RDouble outnstep_read3 = 0;
    RDouble outnstep_read4 = 0;
    ofile<<"DataContainer object invoke read function to verify its value had been changed" <<endl;
    do 
    {
        if (i==5)
        {
            cdata->moveto_begin(); 
            tecdata->moveto_begin();
            cdata_data->moveto_begin();
            tecdata_data->moveto_begin();
        }
        i--;        
        ofile<<"i = "<< i <<endl;
        cdata->read(&outnstep_read1,sizeof(RDouble));
        ofile<<"source cdata ";
        ofile<<" outnstep_read1       = "<<outnstep_read1<<endl;
        tecdata->read(&outnstep_read2,sizeof(RDouble));
        ofile<<"source tecdata";
        ofile<<" outnstep_read2      = "<<outnstep_read2<<endl;
        cdata_data->read(&outnstep_read3,sizeof(RDouble));
        ofile<<"target cdata_data ";
        ofile<<" outnstep_read3  = "<<outnstep_read3<<endl;
        tecdata_data->read(&outnstep_read4,sizeof(RDouble));
        ofile<<"target tecdata_data";
        ofile<<" outnstep_read4 = "<<outnstep_read4<<endl;
    } while (i);


    ParallelCloseFile(ifile);
    ParallelCloseFile(ofile);

    delete region;
    delete testParam;
    delete actkey;
    delete cdata_data;
    delete tecdata_data;

    return 0;
}

