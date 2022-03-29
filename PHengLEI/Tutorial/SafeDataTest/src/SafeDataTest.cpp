//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      SafeDataTest.cpp
//! @brief     It tests functions of the basic class about data storage,
//!			   check whether they can run as required or not.
//! @author    Meng Liyuan.

#include <vector>
#include <iostream>
#include <fstream>
#include "TK_Log.h"
#include "Data_SafeData.h"
#include "Data_Param.h"
#include "TK_Parse.h"
#include "PHIO.h"
#include "TK_Time.h"
#include "Constants.h"
using namespace PHSPACE;
using namespace std;

int main()
{
    //! Record time.
    setvbuf( stdout, NULL, _IONBF, 0 );
    double *ptr_time = NULL;
    PHSPACE::TIME_SPACE::GetWallTimeForResidual();
    PHSPACE::TIME_SPACE::GetWallTime();
    PHSPACE::TIME_SPACE::GetCPUTime(ptr_time);

    //! Perform the initialization operation of PHengLEI.
    PHMPI::Initialization();

    //! Define the input and output files which to store parameters and results.
    fstream ifile,ofile;
    ParallelOpenFile(ofile, "./out/SafeDataTestOutput.txt", ios::out);
    string datafile = "./in/Data.hypara";
    ifile.open(datafile.c_str(), ios_base::in);
    if ( !ifile )
    {
        cout << "could not open " << datafile << "\n";
        exit(0);
    }

    //! Read data from custom-defined parameter file. 
    Data_Param *testParam = new Data_Param();
    TK_Parse::ReadBasicData(ifile, testParam);
    
    //! Define serial parameters using in the test process. 
    //! Obtain the data from the database.
    string name = "value";
    int type = 1,size = 1;
    testParam->GetData("name", &name, PHSTRING, 1);    
    testParam->GetData("type", &type, PHINT,    1);
    testParam->GetData("size", &size, PHINT,    1);

    //! Call the functions about data storage.
    Data_SafeData *datas =new Data_SafeData("value",&name,type,size);
    datas->ShowData();
    string dataname = datas->GetName();
    int    datatype = datas->GetType();
    int    datasize = datas->GetSize();    
    ofile << "The data's name is : " << dataname << endl;
    ofile << "The data's type is : " << datatype << endl;
    ofile << "The data's size is : " << datasize << endl;

    //! Close the output file and exit PHengLEI.
    ParallelCloseFile(ofile);

    PHMPI::FreeBlockData();
    TIME_SPACE::WriStatTime();
    PH_Finalize();
    return 0;
}