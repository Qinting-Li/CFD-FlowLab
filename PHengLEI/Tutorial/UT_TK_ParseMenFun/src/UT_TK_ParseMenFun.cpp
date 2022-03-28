//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      GeoFaceTopoStruct.cpp
//! @brief     It calculates the node topology of the structured grid,
//!			   and then extract the node topology information,
//!	           such as node-cell and node-face connected relationship.
//! @author    Meng Liyuan.

#include "TK_Log.h"
#include "Glb_Dimension.h"
#include "Geo_NodeTopo_Struct.h"
#include "Region.h"
#include "PHIO.h"
#include "TK_Parse.h"
#include "GlobalDataBase.h"
using namespace PHSPACE;
using namespace PHMPI;

int main()
{
    fstream ifile,ofile;
    ParallelOpenFile(ofile, "./out/TK_ParseMenFunOutput.txt", ios::out);
    string datafile = "./bin/Data.hypara";
    ifile.open(datafile.c_str(), ios_base::in);

    if ( !ifile )
    {
        cout << "could not open " << datafile << "\n";
        exit(0);
    }

    Data_Param * testParam = new Data_Param();

    //! Read data from custom-defined parameter file. 	
    TK_Parse::ReadBasicData(ifile, testParam);
    ParallelCloseFile(ifile);

    string str_LowString ="";
    string str_UpString ="";
    string str_TrimBlanks="";
    string str_SkipLines="";
    string str_IsEmptyLine="";
    string str_ReadNextNonEmptyLine="";
    string str_FindNextWord="";
    string str_IsArrayParameter="";
    string str_ReadNextLine="";
    testParam->GetData("LowString", &str_LowString, PHDOUBLE, 1);
    testParam->GetData("UpString", &str_UpString, PHDOUBLE, 1);
    testParam->GetData("TrimBlanks", &str_TrimBlanks, PHDOUBLE, 1);
    testParam->GetData("SkipLines", &str_SkipLines, PHDOUBLE, 1);
    testParam->GetData("IsEmptyLine", &str_IsEmptyLine, PHDOUBLE, 1);
    testParam->GetData("ReadNextNonEmptyLine", &str_ReadNextNonEmptyLine, PHDOUBLE, 1);
    testParam->GetData("FindNextWord", &str_FindNextWord, PHDOUBLE, 1);
    testParam->GetData("IsArrayParameter", &str_IsArrayParameter, PHDOUBLE, 1);
    testParam->GetData("ReadNextLine", &str_ReadNextLine, PHDOUBLE, 1);

    ofile<<" LowString: "<<str_LowString<<endl;
    ofile<<" UpString: "<<str_UpString<<endl;
    ofile<<" TrimBlanks: "<<str_TrimBlanks<<endl;
    ofile<<" SkipLines: "<<str_SkipLines<<endl;
    ofile<<" IsEmptyLine: "<<str_IsEmptyLine<<endl;
    ofile<<" ReadNextNonEmptyLine: "<<str_ReadNextNonEmptyLine<<endl;
    ofile<<" FindNextWord: "<<str_FindNextWord<<endl;
    ofile<<" IsArrayParameter: "<<str_IsArrayParameter<<endl;
    ofile<<" ReadNextLine: "<<str_ReadNextLine<<endl<<endl;

    ofile << "invoke ToLowerCase"<<endl;
    ToLowerCase(str_LowString);
    ofile<<" str_LowString : "<<str_LowString<<endl<<endl;

    ofile << "invoke ToUpperCase"<<endl;
    ToUpperCase(str_UpString);
    ofile<<" str_UpString : "<<str_UpString<<endl<<endl;

    ofile << "invoke TrimBlanks"<<endl;
    string TrimBlanks_str = "     "+str_TrimBlanks+"    "+str_TrimBlanks+str_TrimBlanks+" ";
    TrimBlanks(TrimBlanks_str);
    ofile<<" str_TrimBlanks : "<<TrimBlanks_str<<endl<<endl;


    //ifile.open(datafile.c_str(), ios_base::in);
    //ofile << "invoke SkipLines"<<endl;
    //SkipLines(ifile,3);

    ofile << "invoke IsEmptyLine"<<endl;
    bool BoolParaIsEmptyLine = IsEmptyLine(str_IsEmptyLine);
    ofile<<" BoolParaIsEmptyLine: "<<BoolParaIsEmptyLine<<endl<<endl;
    //string word = "Find";
    //string separator =" \r\n\t";
    //string  FindNextWordStr = "";
    //FindNextWordStr = FindNextWord(str_IsEmptyLine,word,separator);
    //ofile<<" FindNextWordStr: "<<FindNextWordStr<<endl;
    ofile << "invoke IsArrayParameter"<<endl;
    bool BoolParaIsArrayParameter = IsArrayParameter(str_IsArrayParameter);
    ofile<<" BoolParaIsArrayParameter: "<<BoolParaIsArrayParameter<<endl<<endl;


    ParallelCloseFile(ifile);
    ParallelCloseFile(ofile);

    delete testParam;

    return 0;
}

