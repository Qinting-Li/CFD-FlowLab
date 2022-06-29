#include <cmath>
#include "PBBoundaryTurbKWall.h"

using namespace std;

namespace PHSPACE
{
PBBoundaryTurbKWall::PBBoundaryTurbKWall(const string &name):PBBoundaryScalarDiffusion(name),cmu_(0.09),
    cmu25_(0.5477), cmu75_(0.1643), yptr_(11.225), elog_(9.793), cappa_(0.4187)
{

}

PBBoundaryTurbKWall::~PBBoundaryTurbKWall()
{

}

void PBBoundaryTurbKWall::UpdateVistb(UnstructGrid* grid, string phiName, vector<int>* faceIndex, Data_Param* bcData)
{
    int nTotalCell = grid->GetNTotalCell();
    int * leftCellOfFace = grid->GetLeftCellOfFace();
    RDouble *wDist = reinterpret_cast<RDouble *>(grid->GetWallDist());
    RDouble *k = reinterpret_cast<RDouble *>(grid->GetDataPtr("kinetic"));
    RDouble *kb = reinterpret_cast<RDouble *>(grid->GetDataPtr("kinetic")) + nTotalCell;
    RDouble *rho = reinterpret_cast<RDouble *>(grid->GetDataPtr("rho"));
    RDouble *visl = reinterpret_cast<RDouble *>(grid->GetDataPtr("visl"));
    RDouble *vist = reinterpret_cast<RDouble *>(grid->GetDataPtr("vist"));
    RDouble *vistb = vist + nTotalCell;
    RDouble *mu = reinterpret_cast<RDouble *>(grid->GetDataPtr("mu"));
    RDouble *mub = mu + nTotalCell;
    RDouble *phi = reinterpret_cast<RDouble *>(grid->GetDataPtr(phiName));
    RDouble *phib = phi + nTotalCell;
    RDouble *yplus = reinterpret_cast<RDouble *>(grid->GetDataPtr("yplus")) + nTotalCell;

    for (int iFace = 0; iFace < faceIndex->size(); ++ iFace)
    {
        int iFacelocal = (*faceIndex)[iFace];
        int le = leftCellOfFace[iFacelocal];
        RDouble wd = wDist[le];
        RDouble cmute = cmu25_ * sqrt(k[le]);
        RDouble yplus = rho[le] * wd * cmute / (visl[le] + 1e-20);
        RDouble uplus = yplus;

        if (yplus > yptr_)
        {
            uplus = log(elog_ * yplus) / cappa_;
        }
        mub[iFacelocal] = visl[le] * yplus / uplus;

        phib[iFacelocal] = phi[le];
    }
}

void PBBoundaryTurbKWall::Discrt(UnstructGrid *grid, string phiName, vector<int> *faceIndex, Data_Param *bcData, RDouble *diffCoe)
{

}

void PBBoundaryTurbKWall::WallModification(UnstructGrid *grid, string phiName, vector<int> *faceIndex, Data_Param *bcData)
{
    int nTotalCell = grid->GetNTotalCell();
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    RDouble *xfn = grid->GetFaceNormalX();
    RDouble *yfn = grid->GetFaceNormalY();
    RDouble *zfn = grid->GetFaceNormalZ();
    RDouble *faceArea = grid->GetFaceArea();

    RDouble *tauWall = reinterpret_cast<RDouble *>(grid->GetDataPtr("tauWall"));
    RDouble *wDist = reinterpret_cast<RDouble *>(grid->GetWallDist());
    RDouble *gen = reinterpret_cast<RDouble *>(grid->GetDataPtr(phiName));
    RDouble *k = reinterpret_cast<RDouble *>(grid->GetDataPtr("kinetic"));
    RDouble *rho = reinterpret_cast<RDouble *>(grid->GetDataPtr("rho"));
    RDouble *visl = reinterpret_cast<RDouble *>(grid->GetDataPtr("visl"));
    RDouble *vist = reinterpret_cast<RDouble *>(grid->GetDataPtr("vist"));
    RDouble *vistb = vist+ nTotalCell;
    RDouble *u = reinterpret_cast<RDouble *>(grid->GetDataPtr("U"));
    RDouble *v = reinterpret_cast<RDouble *>(grid->GetDataPtr("V"));
    RDouble *w = reinterpret_cast<RDouble *>(grid->GetDataPtr("W"));
    RDouble *mu = reinterpret_cast<RDouble *>(grid->GetDataPtr("mu")) ;
    RDouble *mub = mu + nTotalCell;
    RDouble *yplus = reinterpret_cast<RDouble *>(grid->GetDataPtr("yplus")) + nTotalCell;
    RDouble *diffCoe = reinterpret_cast<RDouble *>(grid->GetDataPtr("diffCoe"));
    RDouble *diffCoeb = diffCoe + nTotalCell;
    RDouble *ub = u + nTotalCell;
    RDouble *vb = v + nTotalCell;
    RDouble *wb = w + nTotalCell;

    for (int iFace = 0; iFace < faceIndex->size(); ++ iFace)
    {
        int iFacelocal = (*faceIndex)[iFace];
        int le = leftCellOfFace[iFacelocal];
        RDouble wd = wDist[le];
        RDouble cmute = cmu25_ * sqrt(k[le]);
        yplus[iFacelocal] = rho[le] * wd * cmute / (visl[le]);

        RDouble uplus = yplus[iFacelocal];

        if (yplus[iFacelocal] > yptr_)
        {
            uplus = log(elog_ * yplus[iFacelocal]) / cappa_;
        }

        RDouble du = ub[iFacelocal] - u[le];
        RDouble dv = vb[iFacelocal] - v[le];
        RDouble dw = wb[iFacelocal] - w[le];

        RDouble area = faceArea[iFacelocal];

        RDouble ax = xfn[iFacelocal] * area;
        RDouble ay = yfn[iFacelocal] * area;
        RDouble az = zfn[iFacelocal] * area;

        RDouble area2 = ax * ax + ay * ay + az * az;
        RDouble dun = (du * ax + dv * ay + dw * az) / area2;
        RDouble taux = du - dun * ax;
        RDouble tauy = dv - dun * ay;
        RDouble tauz = dw - dun * az;
        RDouble tauwall = sqrt(taux * taux + tauy * tauy + tauz * tauz);

        gen[le] = 0;
        gen[le] += rho[le] * tauwall * cmute * cmute / (cappa_ * wd * uplus);
    }
}

DEFINE_PBSOLVER_REGISTER_FUNCTION(PBBoundary, PBBoundaryTurbKWall, TURB_K_WALL);
} // PB_SPACE