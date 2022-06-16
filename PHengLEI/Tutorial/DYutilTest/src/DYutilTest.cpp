//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      DYutilTest.cpp
//! @brief     It tests the functions of the header file DYutil.h,
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
#include "Math_BasisFunction.h"

using namespace PHSPACE;
using namespace std;

int main()
{
    //! Define the input and output files which to store parameters and results .
    fstream ifile,ofile;
    ParallelOpenFile(ofile, "./out/DYutilOutput.txt", ios::out);
    string datafile = "./in/Data.hypara";
    ifile.open(datafile.c_str(), ios_base::in);

    if ( !ifile )
    {
        cout << "could not open " << datafile << "\n";
        exit(0);
    }

    Data_Param *testParam = new Data_Param();
    TK_Parse::ReadBasicData(ifile, testParam);

    //! Define serial parameters.
    double x1, y1, z1, x2, y2, z2;

    //! Obtain the data from database.
    testParam->GetData("x1", &x1, PHDOUBLE, 1);
    testParam->GetData("y1", &y1, PHDOUBLE, 1);
    testParam->GetData("z1", &z1, PHDOUBLE, 1);
    testParam->GetData("x2", &x2, PHDOUBLE, 1);
    testParam->GetData("y2", &y2, PHDOUBLE, 1);
    testParam->GetData("z2", &z2, PHDOUBLE, 1);

    //! Calculate the basis math functions.
    double Dist, Deter;
    Dist = ComputDist( x1, y1, z1, x2, y2, z2);
    ofile<<"The distance is : "<< Dist <<endl;

    double A[3] ={0,0,0};
    double B[3] ={0,0,0};
    double C[3] ={0,0,0};
    testParam->GetData("A", &A, PHDOUBLE, 3);
    testParam->GetData("B", &B, PHDOUBLE, 3);
    testParam->GetData("C", &C, PHDOUBLE, 3);
    Deter = Determinant( A[0], A[1], A[2], B[0], B[1], B[2], C[0], C[1], C[2]);
    ofile<<"The determinant is : "<< Deter <<endl;

    //! Close the input and output files.
    ParallelCloseFile(ifile);
    ParallelCloseFile(ofile);
    return 0;

}