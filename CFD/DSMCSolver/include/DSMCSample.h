//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//          PPPPP  H   H  EEEEE  N    N  GGGGG  L      EEEEE  III         +
//          P   P  H   H  E      NN   N  G      L      E       I          +
//          PPPPP  HHHHH  EEEEE  N N  N  G  GG  L      EEEEE   I          +
//          P      H   H  E      N  N N  G   G  L      E       I          +
//          P      H   H  EEEEE  N    N  GGGGG  LLLLL  EEEEE  III         +
//------------------------------------------------------------------------+
//          Platform for Hybrid Engineering Simulation of Flows           +
//          China Aerodynamics Research and Development Center            +
//                     (C) Copyright, Since 2010                          +
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      DSMCSample.h
//! @brief     Sample operator for DSMC.
//! @author    Li Jin.

#pragma once
#include "Precision.h"

namespace PHSPACE
{

class Sample 
{
private:   
    RDouble ** CSR, ** SVIB, ** VAR, ** nVAR, ** VARS, ** nVARS;
    RDouble *** CS, *** CSV, *** CSVS, *** VARSP, ** nVARSP, ** DF;
    RDouble *** CSS1, *** CSS2;
    RDouble ** CSDV, ** CSDE, ** CSDVV;

public:
    Sample();
    ~Sample();
    void SetCSR(int size1,int size2);
    RDouble ** GetCSR() const {return this->CSR;}

    void SetCS(int size1, int size2, int size3);
    RDouble *** GetCS() const { return this->CS;}

    void SetCSV(int size1, int size2, int size3);
    RDouble *** GetCSV() const { return this->CSV;}

    void SetCSS1(int size1, int size2, int size3);
    RDouble *** GetCSS1() const { return this->CSS1;}

    void SetCSS2(int size1, int size2, int size3);
    RDouble *** GetCSS2() const { return this->CSS2;}

    RDouble ** GetVAR() const {return this->VAR;}
    RDouble *** GetVARSP() const {return this->VARSP;}
};

}



