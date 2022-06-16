//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      TK_parseTest.cpp
//! @brief     It tests the functions about database , 
//!            such as data acquisition and updating ,  
//!			   check whether they can run as required or not.
//! @author    Meng Liyuan.

#include <vector>
#include <iostream>
#include <fstream>
#include "TK_Log.h"
#include "PHIO.h"
#include "Data_Param.h"
#include "TK_Parse.h"
#include "TK_Time.h"
#include "Constants.h"

using namespace PHSPACE;
using namespace std;

int main()
{
    //! Define the input and output files which to store parameters and results.
    fstream ifile,ofile;
    ParallelOpenFile(ofile, "./out/TK_ParseTestOutput.txt", ios::out);
    string datafile = "./in/Data.hypara";
    ifile.open(datafile.c_str(), ios_base::in);

    if ( !ifile )
    {
        cout << "could not open " << datafile << "\n";
        exit(0);
    }

    Data_Param * testParam = new Data_Param();

    //! Read data from custom-defined parameter file. 	
    TK_Parse::ReadBasicData(ifile, testParam);
    double a, b, c, d;

    //! Define serial parameters. 
    //! Obtain the data from database.
    testParam->GetData("a", &a, PHDOUBLE, 1);
    testParam->GetData("b", &b, PHDOUBLE, 1);
    testParam->GetData("c", &c, PHDOUBLE, 1);
    testParam->GetData("d", &d, PHDOUBLE, 1);
    a = 4.5;
    testParam->UpdateData("a", &a, PHDOUBLE, 1);
    ofile<<"The updated value of a is: "<< a <<endl;
    ofile<<"The value of b is: "<< b <<endl;
    ofile<<"The value of c is: "<< c <<endl;

    //! Call the function which to check the existence of data.
    bool YorN = testParam->IsExist("d", PHDOUBLE, 1);
    if ( YorN == 0 )
        ofile<<"The value of d does not exist! "<<endl;
    else 
        ofile<<"The value of d exist! "<<endl;

    //! Close the input and output files.
    ParallelCloseFile(ifile);
    ParallelCloseFile(ofile);
    return 0;

}