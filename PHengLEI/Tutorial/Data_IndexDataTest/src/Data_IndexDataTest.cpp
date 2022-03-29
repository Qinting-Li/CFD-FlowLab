//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      Data_IndexDataTest.cpp
//! @brief     It tests functions of the basic class about data storage that only stores the pointer to the data,
//!			   check whether they can run as required or not.
//! @author    Meng Liyuan.

#include <vector>
#include <iostream>
#include <fstream>
#include "TK_Log.h"
#include "Data_IndexData.h"
#include "Data_Param.h"
#include "TK_Parse.h"
#include "PHIO.h"
#include "TK_Time.h"
#include "Constants.h"

using namespace PHSPACE;
using namespace std;

int main()
{
    //! Define the input and output files which to store parameters and results.
    fstream ifile,ofile;
    ParallelOpenFile(ofile, "./out/Data_IndexDataTestOutput.txt", ios::out);
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
    testParam->GetData("name", &name, PHSTRING, 1);

    //! Call the functions about data storage.
    Data_IndexData *data = new Data_IndexData(name,&name);

    string dataname = data->GetName();

    data->GetData();
    
    ofile << "The data's name is : "<< dataname << endl;

    //! Close the output file.
    ParallelCloseFile(ofile);
    
    return 0;
}