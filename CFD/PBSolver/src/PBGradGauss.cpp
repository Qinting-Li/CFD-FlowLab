#include <iostream>
#include "PBGradGauss.h"
#include "Geo_UnstructGrid.h"

using namespace std;

namespace PHSPACE
{
PBGradGauss::PBGradGauss(const string &name):PBGrad(name)
{

}

PBGradGauss::~PBGradGauss()
{

}

void PBGradGauss::Calc(UnstructGrid *grid, string phiName, string dphidxName, string dphidyName, string dphidzName, bool isComm)
{
    RDouble *phi = reinterpret_cast<RDouble *>(grid->GetDataPtr(phiName));
    RDouble *wfL= reinterpret_cast<RDouble *>(grid->GetDataPtr("deltL"));
    RDouble *wfR = reinterpret_cast<RDouble *>(grid->GetDataPtr("deltR"));
    RDouble *dphidx = reinterpret_cast<RDouble *>(grid->GetDataPtr(dphidxName));
    RDouble *dphidy = reinterpret_cast<RDouble *>(grid->GetDataPtr(dphidyName));
    RDouble *dphidz = reinterpret_cast<RDouble *>(grid->GetDataPtr(dphidzName));
    grid->CompGradientGGCell2(phi, wfL, wfR, dphidx, dphidy, dphidz);
}

DEFINE_PBSOLVER_REGISTER_FUNCTION(PBGrad, PBGradGauss, GRAD_GAUSS);
}