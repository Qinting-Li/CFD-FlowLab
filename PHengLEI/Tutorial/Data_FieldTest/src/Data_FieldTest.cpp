//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      Data_FieldTest.cpp
//! @brief     It tests functions of the basic class that stores field data by name,
//!			   check whether they can run as required or not.
//! @author    Meng Liyuan.

#include <vector>
#include <iostream>
#include <fstream>
#include "TK_Log.h"
#include "Data_Field.h"
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
    ParallelOpenFile(ofile, "./out/Data_FieldTestOutput.txt", ios::out);
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
    string name;
    string ID;
    testParam->GetData("name", &name, PHSTRING, 1);
    testParam->GetData("ID", &ID, PHSTRING, 1);

    string nameNew = "valueNew";
    ofile << "The source data's name is : "<< name << endl;

    //! Call the functions about field data storage.
    Data_Field *Field = new Data_Field();

    Field->UpdateDataPtr("nameNew",&nameNew);
    Field->UpdateDataPtr(ID,&ID);
    
    ofile << "The new data's name is : "<< nameNew << endl;
    ofile << "The data's ID is : "<< ID << endl;

    Field->GetDataPtr(name);
    Field->DeleteDataPtr(name);

    //! Close the output file.
    ParallelCloseFile(ofile);
    
    return 0;
}