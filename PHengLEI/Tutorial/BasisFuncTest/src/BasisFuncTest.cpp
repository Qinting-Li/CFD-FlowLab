//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      BasisFuncTest.cpp
//! @brief     It tests the basis math functions,
//!			   check whether they can run as required or not.
//! @author    Meng Liyuan.

#include <vector>
#include <iostream>
#include <fstream>
#include "TK_Log.h"
#include "Math_BasisFunction.h"
#include "PHIO.h"
#include "Data_Param.h"
#include "TK_Parse.h"
#include "Constants.h"

using namespace PHSPACE;
using namespace std;

int main()
{
    //! Define the input and output files which to store parameters and results.
    fstream ifile,ofile;
    ParallelOpenFile(ofile, "./out/BasisFuncOutput.txt", ios::out);
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
    double a,b,c;
    testParam->GetData("a", &a, PHDOUBLE, 1);
    testParam->GetData("b", &b, PHDOUBLE, 1);
    testParam->GetData("c", &c, PHDOUBLE, 1);

    //! Call the basis math functions.
    double mn, mx, abs, sqr, pow, sqrDis, Dis, sign;

    mn     = MIN(a, b);
    ofile<<"The Min is: "                   << mn <<endl;
    mx     = MAX(a, b);
    ofile<<"The Max is: "                   << mx <<endl; 
    abs    = ABS(a);
    ofile<<"The Absolute Value of a is: "   << abs <<endl;
    sqr    = SQR(a);
    ofile<<"The Square Value of a is: "     << sqr <<endl;
    pow    = POWER3(a);
    ofile<<"The three Power of a is: "      << pow <<endl;
    sqrDis = SQR(a, b, c);
    ofile<<"The Square distance of abc is: "<< sqrDis <<endl;
    Dis    = DISTANCE(a, b, c);
    ofile<<"The Distance of abc is: "       << Dis <<endl;
    SWAP(a, b);
    ofile<<"The Value of a is: "            << a <<endl;
    ofile<<"The Value of b is: "            << b <<endl;
    sign   = SIGN(a, b);
    ofile<<"The Sign Value of ab is: "      << sign <<endl;

    double v1[3] ={0,0,0};
    double v2[3] ={0,0,0};
    testParam->GetData("v1", &v1, PHDOUBLE, 3);
    testParam->GetData("v2", &v2, PHDOUBLE, 3);

    vector <double> crossVec(3,0);
    CrossProduct(v1, v2, &crossVec[0]);
    ofile<<"The cross product of two vectors is: "<<" ";

    vector<double>::iterator t;
    for(t=crossVec.begin();t!=crossVec.end();++t)
    {
        ofile<<*t<<" ";
    }

    //! Close the input and output files.
    ParallelCloseFile(ifile);
    ParallelCloseFile(ofile);
    return 0;
}