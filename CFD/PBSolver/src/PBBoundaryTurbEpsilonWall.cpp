#include "PBBoundaryTurbEpsilonWall.h"

using namespace std;

namespace PHSPACE
{
PBBoundaryTurbEpsilonWall::PBBoundaryTurbEpsilonWall(const string &name)
    :PBBoundaryScalar(name), cmu75_(0.1643), cappa_(0.4187), cmu(0.09), cmu25_(0.5477), yptr_(11.225)
{

}

PBBoundaryTurbEpsilonWall::~PBBoundaryTurbEpsilonWall()
{

}

void PBBoundaryTurbEpsilonWall::Discrt(UnstructGrid *grid, string phiName, vector<int> *faceIndex, Data_Param *bcData, RDouble *diffCoe)
{
    int nTotalCells = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int *leftCellOfFace  = grid->GetLeftCellOfFace();
    int *faceNumberOfEachCell = grid->GetFaceNumberOfEachCell();
    RDouble *phi = reinterpret_cast<RDouble *>(grid->GetDataPtr(phiName));
    RDouble *phib = phi + nTotalCells;
    RDouble *diffCoeb = diffCoe + nTotalCells;
    RDouble *wDist = reinterpret_cast<RDouble *>(grid->GetWallDist());
    RDouble *k = reinterpret_cast<RDouble *>(grid->GetDataPtr("kinetic"));
    RDouble *ajLeft = reinterpret_cast<RDouble *>(grid->GetDataPtr("ajLeft"));
    RDouble *ajRight = reinterpret_cast<RDouble *>(grid->GetDataPtr("ajRight"));
    RDouble *apPhi = reinterpret_cast<RDouble *>(grid->GetDataPtr("ap"));
    RDouble *scPhi = reinterpret_cast<RDouble *>(grid->GetDataPtr("sc"));
    RDouble *dphidx = reinterpret_cast<RDouble *>(grid->GetDataPtr("dphidx"));
    RDouble *dphidy = reinterpret_cast<RDouble *>(grid->GetDataPtr("dphidy"));
    RDouble *dphidz = reinterpret_cast<RDouble *>(grid->GetDataPtr("dphidz"));

    RDouble *rho = reinterpret_cast<RDouble *>(grid->GetDataPtr("rho"));
    RDouble *visl = reinterpret_cast<RDouble *>(grid->GetDataPtr("visl"));

    for (std::size_t i = 0; i < faceIndex->size(); ++ i)
    {
        const int bface = (*faceIndex)[i];
        const int ic = leftCellOfFace[bface];

        RDouble cmute = this->cmu25_ * std::sqrt(k[ic]);
        RDouble yplus = rho[ic] * wDist[ic] * cmute / (visl[ic] + 1e-20);

        if (yplus > this->yptr_)
        {
            scPhi[ic] = this->cmu75_ * std::pow(k[ic], 1.5) / (this->cappa_ * wDist[ic] + 1e-20);

        }
        else
        {
            scPhi[ic] = (this->cmu * std::pow(k[ic], 2) * rho[ic]) / (visl[ic] + 1e-20);
        }
        apPhi[ic] = 1.0;

        for (int i = 0; i < faceNumberOfEachCell[ic]; ++ i)
        {
            int iFace = grid->GetCell2Face()[ic][i];
            int le    = leftCellOfFace[iFace];

            if (iFace < nBoundFace)
            {
                continue;
            }

            int iFacelocal = iFace - nBoundFace;

            //£¡Correct.
            if ((le == ic))
            {
                ajLeft[iFacelocal] = 0.0;
            }
            else
            {
                ajRight[iFacelocal] = 0.0;

            }
        }
    }
}

//£¡Correct:update epsilon.
void PBBoundaryTurbEpsilonWall::Update(UnstructGrid *grid, string phiName, vector<int> *faceIndex, Data_Param *bcData)
{
    int nTotalCells = grid->GetNTotalCell();
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    RDouble *phi = reinterpret_cast<RDouble *>(grid->GetDataPtr(phiName));
    RDouble *phib = phi + nTotalCells;
    for (size_t iFace = 0; iFace < faceIndex->size(); ++ iFace)
    {
        int iFacelocal = (*faceIndex)[iFace];
        int le = leftCellOfFace[iFacelocal];
        phib[iFacelocal] = phi[le];
    }
}

DEFINE_PBSOLVER_REGISTER_FUNCTION(PBBoundary, PBBoundaryTurbEpsilonWall, TURB_EPSILON_WALL);
} //! PB_SPACE