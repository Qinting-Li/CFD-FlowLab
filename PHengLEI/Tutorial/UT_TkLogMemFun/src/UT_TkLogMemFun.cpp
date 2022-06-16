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
#include "TK_Time.h"
#include "GlobalDataBase.h"
using namespace PHSPACE;
using namespace PHMPI;

int main()
{
    fstream fs;
    OpenLogFile(PHMPI::GetCurrentProcessorID(),fs);
    fs<<"OpenLogFile Test!"<<endl;

    CloseLogFile(fs);

    WriteLogFile("WriteLogFile String Test!");

    ostringstream ostr1;
    ostr1<<"WriteLogFile ostringstream Test!"<<endl;
    WriteLogFile(ostr1);

    WriteWarningLogFile("WriteWarningLogFile");//OpenWarningLogFile,OpenLogFile

    WriteLogFile("WriteLogFile by string para");//CloseLogFile
    ostringstream ostr2;
    ostr2<<"WriteLogFile by ostringstream para"<<endl;
    WriteLogFile(ostr2);

    string template_function="template_function";
    WriteLogFile(template_function,"_3_para");

    WriteLogFile(template_function,"_4","_para");

    WriteLogFile(template_function,"_5","_para","_XXX");

    //PHMakeDir("/test/LogTest");

    ProgressMonitorToWindows(1, 1000, "wall distance");

    ProgressMonitorToLogFile(1000, 1000, "wall distance");

    return 0;
}

