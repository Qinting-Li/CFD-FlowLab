//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      Data_UtilTest.cpp
//! @brief     It tests some util functions about copying data,
//!			   check whether they can run as required or not.
//! @author    Meng Liyuan.

#include <vector>
#include <iostream>
#include <fstream>
#include "TK_Log.h"
#include "Data_Util.h"
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
    ParallelOpenFile(ofile, "./out/Data_UtilTestOutput.txt", ios::out);
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

    int size = 1;
    double Atarget[20], Asource[20], Bsource[20];
    testParam->GetData("size", &size, PHINT, 1);
    testParam->GetData("Asource", &Asource, PHDOUBLE, size);
    testParam->GetData("Bsource", &Bsource, PHDOUBLE, size);
    void *Btarget = Bsource;

    //! Call the functions about data storage.
    Data_Util *Du;
    Du->copy_data_sub  (Btarget,  &Bsource, PHDOUBLE, size);
    Du->copy_data_exist(&Atarget, &Asource, PHDOUBLE, size);
    double *Btarget2 = (double *)Btarget;
    for ( int i = 0; i < size; ++ i )
    {
        ofile << "The Atarget["<< i<< "] is: " << Atarget[i]<<endl;
        ofile << "The Btarget["<< i<< "] is: " << Btarget2[i]<<endl;
    }

    //! Close the output file and exit PHengLEI.
    ParallelCloseFile(ifile);
    ParallelCloseFile(ofile);

    PHMPI::FreeBlockData();
    TIME_SPACE::WriStatTime();
    PH_Finalize();

    return 0;
}