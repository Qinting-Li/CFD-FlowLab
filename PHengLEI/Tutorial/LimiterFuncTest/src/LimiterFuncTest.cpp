//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      LimiterFuncTest.cpp
//! @brief     It tests the limiter math functions,
//!			   check whether they can run as required or not.
//! @author    Meng Liyuan.

#include <vector>
#include <iostream>
#include <fstream>
#include "TK_Log.h"
#include "Math_Limiter.h"
#include "PHIO.h"
#include "GlobalDataBase.h"
#include "Data_Param.h"
#include "TK_Parse.h"
#include "TK_Time.h"
#include "Constants.h"
#include <bitset>
using namespace PHSPACE;
using namespace std;

int main()
{
    //! Define the input and output files which to store parameters and results.
    fstream ifile,ofile;
    ParallelOpenFile(ofile, "./out/LimiterFuncOutput.txt", ios::out);
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
    RFloat vmin,vmax,value,lx,ly,eps;
    testParam->GetData("value", &value, PHFLOAT, 1);
    testParam->GetData("vmin",  &vmin,  PHDOUBLE, 1);
    testParam->GetData("vmax",  &vmax,  PHDOUBLE, 1);

    testParam->GetData("lx",  &lx,  PHDOUBLE, 1);
    testParam->GetData("ly",  &ly,  PHDOUBLE, 1);
    testParam->GetData("eps", &eps, PHDOUBLE, 1);

    //! Call the basis math functions.
    double nolim, minmod, vanleer, mixMV, vanalbada, vanalbadaclz, smooth, venfun;

    ofile<<"The returned value of CheckLimit function is:  "<< CheckLimit( value, vmin, vmax)<<endl;

    nolim        = NOLIMITER( lx, ly );
    minmod       = MINMOD( lx, ly );
    vanleer      = VANLEER( lx, ly );
    mixMV        = MIX_MINMOD_VANLEER( lx, ly);
    vanalbada    = VANALBADA( lx, ly );
    vanalbadaclz = VANALBADACLZ( lx, ly );
    smooth       = SMOOTH( lx, ly );
    venfun       = VenFun(lx, ly, eps);
 
    //string str;
    //char a[50];
    //itoa(nolim,a,2);
    //cout<<a;
    bitset<32> t(nolim);
    cout<<"t = "<<t.to_ulong()<<"="<<t.to_string()<<""<<endl;

    ofile<<"The returned value of NOLIMITER function is:  "          << nolim <<endl;  
    ofile<<"The returned value of MINMOD function is:  "             << minmod <<endl;
    ofile<<"The returned value of VANLEER function is:  "            << vanleer <<endl;
    ofile<<"The returned value of MIX_MINMOD_VANLEER function is:  " << mixMV <<endl;
    ofile<<"The returned value of VANALBADA function is:  "          << vanalbada <<endl;
    ofile<<"The returned value of VANALBADACLZ function is:  "       << vanalbadaclz <<endl;
    ofile<<"There turned value of SMOOTH function is:  "             << vanalbadaclz <<endl;
    ofile<<"The returned value of VenFun function is:  "             << venfun <<endl;

    //! Close the input and output files.
    ParallelCloseFile(ifile);
    ParallelCloseFile(ofile);
    return 0;
}