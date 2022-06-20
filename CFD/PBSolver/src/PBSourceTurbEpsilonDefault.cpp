#include <iostream>
#include "PBSourceTurbEpsilonDefault.h"
#include "PBGrad.h"

using namespace std;

namespace PHSPACE
{
PBSourceTurbEpsilonDefault::PBSourceTurbEpsilonDefault(const string &name):PBSource(name), c1m_(1.44), c2m_(1.92)
{

}

PBSourceTurbEpsilonDefault::~PBSourceTurbEpsilonDefault()
{

}

void PBSourceTurbEpsilonDefault::InitMemory(UnstructGrid *grid,string phiName)
{

}

void PBSourceTurbEpsilonDefault::Discrt(UnstructGrid *grid,string phiName)
{
    int nTotalCell = grid->GetNTotalCell();
    RDouble *k = reinterpret_cast<RDouble *>(grid->GetDataPtr("kinetic"));
    RDouble *gen = reinterpret_cast<RDouble *>(grid->GetDataPtr("gen"));
    RDouble *rho = reinterpret_cast<RDouble *>(grid->GetDataPtr("rho"));
    RDouble *epsilon = reinterpret_cast<RDouble *>(grid->GetDataPtr("epsilon"));
    RDouble *vol = reinterpret_cast<RDouble *>(grid->GetCellVolume());
    RDouble *ap = reinterpret_cast<RDouble *>(grid->GetDataPtr("ap"));
    RDouble *sc = reinterpret_cast<RDouble *>(grid->GetDataPtr("sc"));

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        RDouble ke = vol[iCell] * epsilon[iCell] / (k[iCell] + 1e-20);
        sc[iCell] += this->c1m_ * gen[iCell] * ke;
        ap[iCell] += this->c2m_ * rho[iCell] * ke;
    }
}

DEFINE_PBSOLVER_REGISTER_FUNCTION(PBSource, PBSourceTurbEpsilonDefault, SOURCE_TURB_EPSILON_DEFAULT);
}