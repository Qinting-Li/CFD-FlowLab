#include "DSMCSample.h"
#include "Pointer.h"
using namespace std;

namespace PHSPACE
{

Sample::Sample()
{
    //RDouble ** CSR, ** SVIB, ** VAR, ** nVAR, ** VARS, **nVARS;
    //RDouble ** CS, ** CSV, ** CSVS, ** VARSP, ** nVARSP, ** DF;
    //RDouble *** CSS;
    //RDouble ** CSDV, ** CSDE, ** CSDVV;
    CSR = NULL;
    CS = NULL;
    CSV = NULL;
    CSS1 = NULL;
    CSS2 = NULL;
    VAR = NULL;
    VARSP = NULL;
}


Sample::~Sample()
{
    DelPointer2 (CSR) ; CSR = NULL;
    DelPointer3 (CS) ; CS = NULL;
    DelPointer3 (CSV) ; CSV = NULL;
    DelPointer3 (CSS1) ; CSS1 = NULL;
    DelPointer3 (CSS2) ; CSS2= NULL;
    DelPointer2 (VAR) ; VAR = NULL;
    DelPointer3 (VARSP) ; VARSP = NULL;
}

void Sample::SetCSR(int nTotalCell,int nSpecies)
{
    CSR = NewPointer2<RDouble>(nTotalCell, nSpecies);
    CS = NewPointer3<RDouble>(8,nTotalCell, nSpecies);
    CSV = NewPointer3<RDouble>(1, nTotalCell, nSpecies);
    CSS1 = NewPointer3<RDouble>(8,nTotalCell, nSpecies);
    CSS2 = NewPointer3<RDouble>(8,nTotalCell, nSpecies);
    VAR = NewPointer2<RDouble>(22,nTotalCell);
    VARSP = NewPointer3<RDouble>(12,nTotalCell, nSpecies);
}

}



