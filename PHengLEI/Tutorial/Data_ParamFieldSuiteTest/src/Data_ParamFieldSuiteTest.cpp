//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      Data_ParamFieldSuiteTest.cpp
//! @brief     It tests functions of the basic class that wraps Data_Param and Field's public interface,
//!			   check whether they can run as required or not.
//! @author    Meng Liyuan.

#include <vector>
#include <iostream>
#include <fstream>
#include "TK_Log.h"
#include "Data_ParamFieldSuite.h"
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
    ParallelOpenFile(ofile, "./out/Data_ParamFieldSuiteTestOutput.txt", ios::out);
    string datafile = "./in/Data.hypara";
    ifile.open(datafile.c_str(), ios_base::in);

    if ( !ifile )
    {
        cout << "could not open " << datafile << "\n";
        exit(0);
    }

    //! Read data from custom-defined parameter file. 
    Data_Param *Param = new Data_Param();
    TK_Parse::ReadBasicData(ifile, Param);
    
    //! Define serial parameters using in the test process. 
    //! Obtain the data from the database.
    int a = 2;
    double b = 3;
    Param->GetData("a", &a, PHINT, 1);
    Param->GetData("b", &b, PHDOUBLE, 1);
    
    //! Call the functions of the Data_ParamFieldSuite class.
    Data_ParamFieldSuite *datastore = new Data_ParamFieldSuite();

    a=3;
    datastore->UpdateData("a", &a, PHINT, 1);
    ofile<<"The updated value of a is: "<< a <<endl;

    b = 1.71828;

    datastore->UpdateData("b", &b, PHFLOAT, 1);
    datastore->GetData("b",&b, PHDOUBLE, 1);
    ofile<<"The value of b is: "<< b <<endl;

    datastore->UpdateDataPtr("b", &b);
    datastore->GetDataPtr("a");

    //! Close the output file.
    ParallelCloseFile(ofile);
    
    return 0;
}