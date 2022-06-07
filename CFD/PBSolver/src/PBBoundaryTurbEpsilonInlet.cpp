#include "PBBoundaryTurbEpsilonInlet.h"

using namespace std;

namespace PHSPACE
{
    PBBoundaryTurbEpsilonInlet::PBBoundaryTurbEpsilonInlet(const string &name):PBBoundaryTurbDefault(name), cmu(0.09)
    {

    }

    PBBoundaryTurbEpsilonInlet::~PBBoundaryTurbEpsilonInlet()
    {

    }

    void PBBoundaryTurbEpsilonInlet::Update(UnstructGrid *grid, string phiName, vector<int> *faceIndex, Data_Param *bcData)
    {
        int nTotalCells = grid->GetNTotalCell();
        RDouble epsilonb = 0.0;
        RDouble *phi = reinterpret_cast<RDouble *>(grid->GetDataPtr("epsilon"));
        RDouble *phib = phi + nTotalCells;
        int *leftCellOfFace = grid->GetLeftCellOfFace();

        for (std::size_t iFace = 0; iFace < faceIndex->size(); ++ iFace)
        {
            const int iFacelocal = (*faceIndex)[iFace];

            if (bcData)
            {
                if (bcData->IsExist(phiName, PHDOUBLE, 1))
                {
                    bcData->GetData(phiName, &epsilonb, PHDOUBLE, 1);
                }
            }
            else
            {
                std::cout << "Please assign a value of epsilon on inlet" << std::endl;;
            }
            phib[iFacelocal] = epsilonb;
        }
    }

    void PBBoundaryTurbEpsilonInlet::Discrt(UnstructGrid *grid, string phiName, vector<int> *faceIndex, Data_Param *bcData, RDouble *diffCoe)
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

    DEFINE_PBSOLVER_REGISTER_FUNCTION(PBBoundary, PBBoundaryTurbEpsilonInlet, TURB_EPSILON_INLET);
} // PB_SPACE