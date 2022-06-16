//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      ExitTest.cpp
//! @brief     It tests several methods to terminate the program,
//!			   check whether they can quit the program as required or not.
//! @author    Meng Liyuan.

#include <iostream>
#include <string>
#include "TK_Exit.h"
#include "PHIO.h"
#include "TK_Log.h"
#include "Data_Param.h"
#include "TK_Parse.h"
#include "TK_Time.h"
#include "Constants.h"

using namespace PHSPACE;

//! Declare and define the functions of several methods that terminate the program .
void ExitPHengLEITest();
void ExceptionExitTest();
void UnexpectedVarValueTest();
void PrintDebugInfoExitTest(); 

void ExitPHengLEITest()
{
    PHMPI::Initialization();
    TK_Exit::ExitPHengLEI();
}
void ExceptionExitTest()
{
    PHMPI::Initialization();
    std::string Err = "Front edge is not exist!";
    TK_Exit::ExceptionExit(Err);
    PHMPI::FreeBlockData();
    TIME_SPACE::WriStatTime();
    PH_Finalize();
}

void UnexpectedVarValueTest()
{
    PHMPI::Initialization();
    TK_Exit::UnexpectedVarValue("force_component",0);
    PHMPI::FreeBlockData();
    TIME_SPACE::WriStatTime();
    PH_Finalize();
}

void PrintDebugInfoExitTest()
{
    PHMPI::Initialization();
    TK_Exit::PrintDebugInfoExit("Error: numberOfTotalCompositeFaces != 0");
    PHMPI::FreeBlockData();
    TIME_SPACE::WriStatTime();
    PH_Finalize();    
}
int main() 
{   
    //! Define the input and output files which to store parameters and results .
    fstream ifile,ofile;
    ParallelOpenFile(ofile, "./out/ExitTestOutput.txt", ios::out);
    string datafile = "./in/input.hypara";
    ifile.open(datafile.c_str(), ios_base::in);

    if ( !ifile )
    {
        cout << "could not open " << datafile << "\n";
        exit(0);
    }
    
    Data_Param *testParam = new Data_Param();
    TK_Parse::ReadBasicData(ifile, testParam);

    //! Define serial parameters. 
    //! Obtain the data from database.
    int nsituation = 0;
    testParam->GetData("nsituation",&nsituation,PHINT,1);

    //! Run function according to the different situation.
    if(nsituation == 1)
    {
        ofile << "The First Situation: Normal Exit PHengLEI !" << endl;
        ExitPHengLEITest();        
    }
    else if(nsituation == 2)
    {
        ofile << "The Second Situation: Abnormal Exit PHengLEI and Dumping Error Information to Log File !" << endl;
        ExceptionExitTest();
    }
    else if(nsituation == 3)
    {
        ofile << "The Third Situation: Normal Exit PHengLEI Because of A Variable's Vaule Is Not Right !" << endl;
        UnexpectedVarValueTest();        
    }
    else if (nsituation == 4)
    {
        ofile << "The Fourth Situation: Normal Exit PHengLEI and Print Debug Information to Log File or Screen !" << endl;
        PrintDebugInfoExitTest();        
    }
    else 
    {
        cout << "Please reset the nsituation value from 1 to 4 !" << endl;
    }

    //! Close the input and output files.
    ParallelCloseFile(ifile);
    ParallelCloseFile(ofile);
    return 0;

}

