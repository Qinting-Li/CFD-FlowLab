//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      TimeSpanTest.cpp
//! @brief     It tests interface class of time computing functions,
//!			   check whether they can be called as required or not.
//! @author    Meng Liyuan.

#include <iostream>
#include <string>
#include "TK_Time.h"
#include "PHIO.h"
#ifdef WIN32
#include "Windows.h"
#else
#include <sys/time.h>
#endif

using namespace PHSPACE;

int main()
{
    //! Define the input and output files which to store parameters and results. 
    fstream ifile,ofile;
    ParallelOpenFile(ofile, "./out/TK_TimeTestOutput.txt", ios::out);
    string datafile = "./in/Data.hypara";
    ifile.open(datafile.c_str(), ios_base::in);

    if ( !ifile )
    {
        cout << "could not open " << datafile << "\n";
        exit(0);
    }

#ifdef WIN32
    TimeSpan t;
    ofile << "Currently, testing time for windows ! " << endl;

    //! Call the time computing functions for windows.
    t.ResetTime();
    t.GetTimeSpan();
    t.GetTimeSpan();	
    Sleep(6000);
    t.ShowTimeSpan();
    t.ShowTimeSpanToLogfile();
#else

    TimeSpan ts;
    ofile << "Currently, testing time for linux ! " << endl;

    //! Call the time computing functions for linux.
    ts.ResetTime();
    ts.GetTimeSpan();
    ts.GetTimeSpan();
    ts.ShowTimeSpan();
    ts.ShowTimeSpanToLogfile();
#endif

    //! Close the input and output files.
    ParallelCloseFile(ifile);
    ParallelCloseFile(ofile);
    return 0;
}

