#include <iostream>
#include "PBSourceTurbSADefault.h"
#include "PBGrad.h"

using namespace std;

namespace PHSPACE
{
    PBSourceTurbSADefault::PBSourceTurbSADefault(const string &name):PBSource(name)
    {
    }

    PBSourceTurbSADefault::~PBSourceTurbSADefault()
    {

    }

    void PBSourceTurbSADefault::InitMemory(UnstructGrid *grid,string phiName)
    {

    }

    void PBSourceTurbSADefault::Discrt(UnstructGrid *grid,string phiName)
    {
        RDouble *k = reinterpret_cast<RDouble *>(grid->GetDataPtr(phiName));
        RDouble *gen = reinterpret_cast<RDouble *>(grid->GetDataPtr("gen"));
        RDouble *ap = reinterpret_cast<RDouble *>(grid->GetDataPtr("ap"));
        RDouble *sc = reinterpret_cast<RDouble *>(grid->GetDataPtr("sc"));
        RDouble *gra = reinterpret_cast<RDouble *>(grid->GetDataPtr("gra"));
        RDouble *des = reinterpret_cast<RDouble *>(grid->GetDataPtr("des"));
        RDouble *vol = reinterpret_cast<RDouble *>(grid->GetCellVolume());

        for (int iCell = 0; iCell < grid->GetNTotalCell(); ++ iCell)
        {
            sc[iCell] += (gen[iCell] + gra[iCell] ) * vol[iCell];
            ap[iCell] +=  des[iCell] / (k[iCell] + 1e-20) * vol[iCell];
        }
    }

    DEFINE_PBSOLVER_REGISTER_FUNCTION(PBSource, PBSourceTurbSADefault, SOURCE_TURB_SA_DEFAULT);
}