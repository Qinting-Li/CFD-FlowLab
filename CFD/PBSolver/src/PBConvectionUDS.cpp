#include <iostream>
#include <algorithm>
#include "PBConvectionUDS.h"
#include "Geo_UnstructGrid.h"
using namespace std;

namespace PHSPACE
{
    PBConvectionUDS::PBConvectionUDS(const string &name):PBConvection(name)
    {

    }

    PBConvectionUDS::~PBConvectionUDS()
    {

    }

    void PBConvectionUDS::Discrt(UnstructGrid *grid, string phiName, string apName, string scName, RDouble *v, RDouble *apV, RDouble *scV, RDouble *w, RDouble *apW, RDouble *scW)
    {
        int faceBegin = grid->GetNBoundFace();
        int faceEnd = grid->GetNTotalFace();
        RDouble *flux = reinterpret_cast<RDouble *>(grid->GetDataPtr("flux"));
        RDouble *ajUpper= reinterpret_cast<RDouble *>(grid->GetDataPtr("ajLeft"));
        RDouble *ajLower= reinterpret_cast<RDouble *>(grid->GetDataPtr("ajRight"));
        RDouble *phi = reinterpret_cast<RDouble *>(grid->GetDataPtr(phiName));
        RDouble *ap = reinterpret_cast<RDouble *>(grid->GetDataPtr(apName));
        RDouble *sc = reinterpret_cast<RDouble *>(grid->GetDataPtr(scName));

        for (int iFace = faceBegin; iFace < faceEnd; ++ iFace)
        {
            int iFacelocal = iFace - faceBegin;
            RDouble a12 = std::max(flux[iFace], 0.0);
            RDouble a21 = a12 - flux[iFace];
            ajUpper[iFacelocal] += a21;
            ajLower[iFacelocal] += a12;
        }
    }

    DEFINE_PBSOLVER_REGISTER_FUNCTION(PBConvection, PBConvectionUDS, CONVECTION_UPWIND);
}