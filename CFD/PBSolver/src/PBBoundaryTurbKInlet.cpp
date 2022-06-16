#include "PBBoundaryTurbKInlet.h"

using namespace std;

namespace PHSPACE
{
    PBBoundaryTurbKInlet::PBBoundaryTurbKInlet(const string &name):PBBoundaryTurbDefault(name)
    {

    }

    PBBoundaryTurbKInlet::~PBBoundaryTurbKInlet()
    {

    }

    void PBBoundaryTurbKInlet::Update(UnstructGrid *grid, string phiName, vector<int> *faceIndex, Data_Param *bcData)
    {

        RDouble kb = 0.0;
        RDouble *phi = reinterpret_cast<RDouble *>(grid->GetDataPtr("kinetic"));
        int nTotalCells = grid->GetNTotalCell();
        RDouble* phib = phi + nTotalCells;
        int *leftCellOfFace = grid->GetLeftCellOfFace();

        for (std::size_t iFace = 0; iFace < faceIndex->size(); ++ iFace)
        {
            int iFacelocal = (*faceIndex)[iFace];

            if (bcData)
            {
                if (bcData->IsExist(phiName, PHDOUBLE, 1))
                {
                    bcData->GetData(phiName, &kb, PHDOUBLE, 1);
                }
            }
            else
            {
                int le = leftCellOfFace[iFacelocal];
                kb = phi[le];
            }
            phib[iFacelocal] = kb;
        }
    }

    void PBBoundaryTurbKInlet::Discrt(UnstructGrid *grid, string phiName, vector<int> *faceIndex, Data_Param *bcData, RDouble *diffCoe)
    {
        int nTotalCells = grid->GetNTotalCell();
        RDouble *phi = reinterpret_cast<RDouble *>(grid->GetDataPtr(phiName));
        RDouble *phib = phi + nTotalCells;
        RDouble *diffCoeb = diffCoe + nTotalCells;
        RDouble *apPhi = reinterpret_cast<RDouble *>(grid->GetDataPtr("ap"));
        RDouble *scPhi = reinterpret_cast<RDouble *>(grid->GetDataPtr("sc"));
        RDouble *dphidx = reinterpret_cast<RDouble *>(grid->GetDataPtr("dphidx"));
        RDouble *dphidy = reinterpret_cast<RDouble *>(grid->GetDataPtr("dphidy"));
        RDouble *dphidz = reinterpret_cast<RDouble *>(grid->GetDataPtr("dphidz"));
        RDouble *fluxb = reinterpret_cast<RDouble *>(grid->GetDataPtr("flux"));

        for (std::size_t iFace = 0; iFace < faceIndex->size(); ++ iFace)
        {
            if (bcData)
            {

                int iFacelocal = (*faceIndex)[iFace]; //! The boundary element serial number is arranged in front of the element face.
                int le = grid->GetLeftCellOfFace()[iFacelocal];

                if (fluxb[iFacelocal] < 0)
                {
                    apPhi[le] -= fluxb[iFacelocal];
                    scPhi[le] -= fluxb[iFacelocal] * phib[iFacelocal];
                }
            }
        }
    }

    DEFINE_PBSOLVER_REGISTER_FUNCTION(PBBoundary, PBBoundaryTurbKInlet, TURB_K_INLET);
} // PB_SPACE