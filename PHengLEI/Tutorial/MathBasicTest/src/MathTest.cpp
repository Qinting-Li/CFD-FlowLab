//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      MathTest.cpp
//! @brief     It tests the math functions,
//!			   check whether they can run as required or not.
//! @author    Meng Liyuan.

#include <vector>
#include <iostream>
#include <fstream>
#include "TK_Log.h"
#include "Math_Limiter.h"
#include "Math_BasisFunction.h"
#include "PHIO.h"
#include "GlobalDataBase.h"
#include "Data_Param.h"
#include "TK_Parse.h"
#include "Constants.h"
using namespace PHSPACE;
using namespace std;

int main()
{
    //! Define the input and output files which to store parameters and results.
    fstream ifile,ofile;
    //ParallelOpenFile(ofile, "F://PHengLEI2020//API0//Math//test//out//output.txt", ios::out);
    ParallelOpenFile(ofile, "./out/output.txt", ios::out);
    string datafile = "./in/input.hypara";
    ifile.open(datafile.c_str(), ios_base::in);

    if ( !ifile )
    {
        ofile << "could not open " << datafile << "\n";
        exit(0);
    }

    //! Read data from custom-defined parameter file. 
    Data_Param *testParam = new Data_Param();
    TK_Parse::ReadBasicData(ifile, testParam);

    double a, b, c;               
    testParam->GetData("a", &a, PHDOUBLE, 1);
    testParam->GetData("b", &b, PHDOUBLE, 1);
    testParam->GetData("c", &c, PHDOUBLE, 1);
    ofile<<"The min value of a and b is: " << MIN(a,b)          <<endl;
    ofile<<"The Max value of a and b is: " << MAX(a,b)          <<endl;
    ofile<<"The Absolute Value of a is: "  << ABS(a)            <<endl;
    ofile<<"The square value of a is: "    << SQR(a)            <<endl;
    ofile<<"The three powers of a is: "    << POWER3(a)         <<endl;
    ofile<<"The square distance of a,b and c is: "<< SQR(a,b,c) <<endl;
    ofile<<"The Distance a,b and c is: "   << DISTANCE(a,b,c)   <<endl;  
    SWAP(a,b);  
    ofile<<" Swap the value of a and b: "  <<" a = "<< a <<" b = "<< b <<endl;
    ofile<<"The sign value of a and b is: "<< SIGN(a,b)         <<endl;
    
    int nTotal,nInterval;
    testParam->GetData("nTotal", &nTotal, PHINT, 1);
    testParam->GetData("nInterval", &nInterval, PHINT, 1);
    ofile<<"The interval value is: "<< GetProgressInterval(nTotal,nInterval) <<endl;

    double v1[3]  = {0.0,0.0,0.0};
    double v2[3]  = {0.0,0.0,0.0};
    double cp1[3] = {0.0,0.0,0.0};
    testParam->GetData("v1", &v1, PHDOUBLE, 3);
    testParam->GetData("v2", &v2, PHDOUBLE, 3);
    CrossProduct(v1, v2, cp1);
    ofile<<"The cross product of v1 and v2 is: "<<"("<<cp1[0]<<","<<cp1[1]<<","<<cp1[2]<<")" <<endl;

    double r1[3]  = {0.0,0.0,0.0};
    double r2[3]  = {0.0,0.0,0.0};
    double r3[3]  = {0.0,0.0,0.0};
    double r4[3]  = {0.0,0.0,0.0};
    double r1n[4] = {0.0,0.0,0.0,0.0};
    double r1X[3] = {0.0,0.0,0.0};
    double r1Y[3] = {0.0,0.0,0.0};
    double r1Z[3] = {0.0,0.0,0.0};
    double degree = 0.0;
    testParam->GetData("r1", &r1, PHDOUBLE, 3);
    testParam->GetData("r2", &r2, PHDOUBLE, 3);
    testParam->GetData("r3", &r3, PHDOUBLE, 3);
    testParam->GetData("r4", &r4, PHDOUBLE, 3);
    testParam->GetData("degree", &degree, PHDOUBLE, 1);
    NormalizeVector(r1,r1n,3);
    RotateAxisX(r1, r1X, degree);
    RotateAxisY(r1, r1Y, degree);
    RotateAxisZ(r1, r1Z, degree);
    ofile<<"The vector that has been normalized of r1: "    <<"("<<r1n[0]<<","<<r1n[1]<<","<<r1n[2]<<"," <<r1n[3]<<")" <<endl;
    ofile<<"The coordinate of vector rotation with X axis: "<<"("<<r1X[0]<<","<<r1X[1]<<","<<r1X[2]<<")" <<endl;
    ofile<<"The coordinate of vector rotation with Y axis: "<<"("<<r1Y[0]<<","<<r1Y[1]<<","<<r1Y[2]<<")" <<endl;
    ofile<<"The coordinate of vector rotation with Z axis: "<<"("<<r1Z[0]<<","<<r1Z[1]<<","<<r1Z[2]<<")" <<endl;

    double p1[2] = {0.0,0.0};
    double p2[2] = {0.0,0.0};
    double uN[3] = {0.0,0.0,0.0};
    testParam->GetData("p1", &p1, PHDOUBLE, 2);
    testParam->GetData("p2", &p2, PHDOUBLE, 2);
    UnitFaceNormal2D(p1, p1, uN);
    ofile<<"The unit normal direction is: "<<"("<<uN[0]<<","<<uN[1]<<","<<uN[2]<<")" <<endl;

    double fN[3] = {0.0,0.0,0.0};
    double area  = 0.0;
    double vol   = 0.0;
    Area2D         (r1 ,r2 ,r3, r4, vol);
    FaceNormal3D   (r1, r2, r3, r4, fN );
    Vol_Tetrahedron(r1, r2, r3, r4, vol);
    ofile<<"The area of quadrilateral is: "<< area <<endl;
    ofile<<"The face normal direction is: "<<"("<<fN[0]<<","<<fN[1]<<","<<fN[2]<<")" <<endl;
    ofile<<"The six-times of the real tetrahedron volume is: "<< vol <<endl;

    RFloat vmin,vmax,Value,lx,ly,eps;
    testParam->GetData("Value",&Value,3, 1);
    testParam->GetData("vmin", &vmin, 3, 1);
    testParam->GetData("vmax", &vmax, 3, 1);
    testParam->GetData("lx",   &lx,   3, 1);
    testParam->GetData("ly",   &ly,   3, 1);
    testParam->GetData("eps",  &eps,  3, 1);
    ofile<<"The returned value of CheckLimit function is:  "  << CheckLimit( Value,vmin,vmax) <<endl;
    ofile<<"The returned value of NOLIMITER function is:  "   << NOLIMITER(lx,ly) <<endl;
    ofile<<"The returned value of MINMOD function is:  "      << MINMOD(lx,ly) <<endl;
    ofile<<"The returned value of VANLEER function is:  "     << VANLEER(lx,ly) <<endl;
    ofile<<"The returned value of MIX_MINMOD_VANLEER function is:  "<< MIX_MINMOD_VANLEER(lx,ly) <<endl;
    ofile<<"The returned value of VANALBADA function is:  "   << VANALBADA( lx,ly ) <<endl;
    ofile<<"The returned value of VANALBADACLZ function is:  "<< VANALBADACLZ(lx,ly) <<endl;
    ofile<<"The returned value of SMOOTH function is:  "      << SMOOTH(lx,ly) <<endl;
    ofile<<"The returned value of SMOOTHTURB function is:  "  << SMOOTHTURB(lx,ly) <<endl;
    ofile<<"The returned value of THIRD_ORDER_SMOOTH function is:  "<< THIRD_ORDER_SMOOTH(lx,ly) <<endl;
    ofile<<"The returned value of THIRD_ORDER_MINMOD_SMOOTH function is:  "<< THIRD_ORDER_MINMOD_SMOOTH(lx,ly) <<endl;
    ofile<<"The returned value of VenFun function is:  "      << VenFun(lx, ly, eps) <<endl;

    //! Close the input and output files.
    ParallelCloseFile(ifile);
    ParallelCloseFile(ofile);
    return 0;
}