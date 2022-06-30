#include "PBBoundaryTurbEpsilonWall.h"

using namespace std;

namespace PHSPACE
{
PBBoundaryTurbEpsilonWall::PBBoundaryTurbEpsilonWall(const string &name):PBBoundaryScalar(name), cmu75_(0.1643),
    cappa_(0.4187),cmu(0.09), cmu25_(0.5477), yptr_(11.225), yvstar_(20.0), cl_(2.55)
{

}

PBBoundaryTurbEpsilonWall::~PBBoundaryTurbEpsilonWall()
{

}

void PBBoundaryTurbEpsilonWall::Update(UnstructGrid *grid, string phiName, vector<int> *faceIndex, Data_Param *bcData)
{
    int nTotalCell = grid->GetNTotalCell();
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    RDouble *xfc = grid->GetFaceCenterX();
    RDouble *yfc = grid->GetFaceCenterY();
    RDouble *zfc = grid->GetFaceCenterZ();
    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    RDouble *phi = reinterpret_cast<RDouble *>(grid->GetDataPtr(phiName));
    RDouble *phib = phi + nTotalCell;
    RDouble *wDist = reinterpret_cast<RDouble *>(grid->GetWallDist());
    RDouble *k = reinterpret_cast<RDouble *>(grid->GetDataPtr("kinetic"));
    RDouble *rho = reinterpret_cast<RDouble *>(grid->GetDataPtr("rho"));
    RDouble *visl = reinterpret_cast<RDouble *>(grid->GetDataPtr("visl"));
    RDouble *yplus = reinterpret_cast<RDouble*>(grid->GetDataPtr("yplus")) + nTotalCell;

    for (int i = 0; i < faceIndex->size(); ++i)
    {
        int bface = (*faceIndex)[i];
        int le = leftCellOfFace[bface];

        RDouble cmute = cmu25_ * sqrt(k[le]);

        RDouble dx = xfc[bface] - xcc[le];
        RDouble dy = yfc[bface] - xcc[le];
        RDouble dz = zfc[bface] - zcc[le];
        RDouble wd = wDist[le];
        RDouble yv = yptr_ * visl[le] / (rho[le] * pow(k[le], 0.5));

        phi[le] = cmu75_ * std::pow(k[le], 1.5) / (cappa_ * wDist[le]);
        phib[bface] =  phi[le];
    }
}


void PBBoundaryTurbEpsilonWall::Discrt(UnstructGrid *grid, string phiName, vector<int> *faceIndex, Data_Param *bcData, RDouble *diffCoe)
{
    int nTotalCells = grid->GetNTotalCell();
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int *leftCellOfFace  = grid->GetLeftCellOfFace();
    int *faceNumberOfEachCell = grid->GetFaceNumberOfEachCell();
    int **cell2face  = grid->GetCell2Face();

    RDouble *phi = reinterpret_cast<RDouble *>(grid->GetDataPtr(phiName));
    RDouble *phib = phi + nTotalCells;
    RDouble *diffCoeb = diffCoe + nTotalCells;
    RDouble *wDist = reinterpret_cast<RDouble *>(grid->GetWallDist());
    RDouble *k = reinterpret_cast<RDouble *>(grid->GetDataPtr("kinetic"));
    RDouble *ajLeft = reinterpret_cast<RDouble *>(grid->GetDataPtr("ajLeft"));
    RDouble *ajRight = reinterpret_cast<RDouble *>(grid->GetDataPtr("ajRight"));
    RDouble *apPhi = reinterpret_cast<RDouble *>(grid->GetDataPtr("ap"));
    RDouble *scPhi = reinterpret_cast<RDouble *>(grid->GetDataPtr("sc"));
    RDouble *yplus = reinterpret_cast<RDouble *>(grid->GetDataPtr("yplus")) + nTotalCells;
    RDouble *rho = reinterpret_cast<RDouble *>(grid->GetDataPtr("rho"));
    RDouble *visl = reinterpret_cast<RDouble *>(grid->GetDataPtr("visl"));

    for (int i = 0; i < faceIndex->size(); ++ i)
    {
        int bface = (*faceIndex)[i];
        int ic = leftCellOfFace[bface];

        RDouble cmute = cmu25_ * sqrt(k[ic]);
        RDouble yv = yptr_ * visl[ic] / (rho[ic] * cmute);

        scPhi[ic] = cmu75_ * pow(k[ic], 1.5) / (cappa_ * wDist[ic]);
        apPhi[ic] = 1.0;


        for (int i = 0; i < faceNumberOfEachCell[ic]; ++ i)
        {
            int iFace = cell2face[ic][i];

            if (iFace < nBoundFace)
            {
                continue;
            }

            int iFacelocal = iFace - nBoundFace;

            if (leftCellOfFace[iFace] == ic)
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

DEFINE_PBSOLVER_REGISTER_FUNCTION(PBBoundary, PBBoundaryTurbEpsilonWall, TURB_EPSILON_WALL);
} 