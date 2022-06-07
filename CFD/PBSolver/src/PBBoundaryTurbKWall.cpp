#include <cmath>
#include "PBBoundaryTurbKWall.h"

using namespace std;

namespace PHSPACE
{
    PBBoundaryTurbKWall::PBBoundaryTurbKWall(const string &name)
        :PBBoundaryScalarDiffusion(name), cmu_(0.09), cmu25_(0.5477), cmu75_(0.1643), yptr_(11.225), elog_(9.793), cappa_(0.4187)
    {

    }

    PBBoundaryTurbKWall::~PBBoundaryTurbKWall()
    {

    }

    //£¡Correct:add discrete item.
    void PBBoundaryTurbKWall::WallModification(UnstructGrid* grid, string phiName, vector<int>* faceIndex, Data_Param* bcData)
    {
        int nTotalCell = grid->GetNTotalCell();
        RDouble *wDist = reinterpret_cast<RDouble *>(grid->GetWallDist());
        RDouble *gen = reinterpret_cast<RDouble *>(grid->GetDataPtr(phiName));
        RDouble *k = reinterpret_cast<RDouble *>(grid->GetDataPtr("kinetic"));
        RDouble *kb = k+nTotalCell;
        RDouble *rho = reinterpret_cast<RDouble *>(grid->GetDataPtr("rho"));
        RDouble *visl = reinterpret_cast<RDouble *>(grid->GetDataPtr("visl"));
        RDouble *u = reinterpret_cast<RDouble *>(grid->GetDataPtr("U"));
        RDouble *v = reinterpret_cast<RDouble *>(grid->GetDataPtr("V"));
        RDouble *w = reinterpret_cast<RDouble *>(grid->GetDataPtr("W"));
        RDouble *mu = reinterpret_cast<RDouble *>(grid->GetDataPtr("mu"));
        RDouble *mub = reinterpret_cast<RDouble *>(grid->GetDataPtr("mu")) + nTotalCell;
        RDouble *ub = u + nTotalCell;
        RDouble *vb = v + nTotalCell;
        RDouble *wb = w + nTotalCell;
        RDouble *dkdx = reinterpret_cast<RDouble *>(grid->GetDataPtr("dkdx"));
        RDouble *dkdy = reinterpret_cast<RDouble *>(grid->GetDataPtr("dkdy"));
        RDouble *dkdz = reinterpret_cast<RDouble *>(grid->GetDataPtr("dkdz"));
        int *wallFlag = reinterpret_cast<int *>(ObjAllocateMemory("wallFlag", nTotalCell, 0.0, grid));
        int *tauWall = reinterpret_cast<int *>(ObjAllocateMemory("tauWall", nTotalCell, 0.0, grid));


        for (std::size_t iFace = 0; iFace < faceIndex->size(); ++ iFace)
        {
            int iFacelocal = (*faceIndex)[iFace];
            int le = grid->GetLeftCellOfFace()[iFacelocal];
            RDouble wd = wDist[le];
            RDouble cmute = this->cmu25_ * std::sqrt(k[le]);
            RDouble yplus = rho[le] * wd * cmute / (visl[le] + 1e-20);
            RDouble cc = visl[le] / wd;
            RDouble uplus = yplus;
            RDouble Fwall = rho[le] * cmute * cmute;
            mub[iFacelocal] = visl[le];
            RDouble du = ub[iFacelocal] - u[le];
            RDouble dv = vb[iFacelocal] - v[le];
            RDouble dw = wb[iFacelocal] - w[le];

            RDouble area = grid->GetFaceArea()[iFacelocal];

            RDouble ax = grid->GetFaceNormalX()[iFacelocal];
            RDouble ay = grid->GetFaceNormalY()[iFacelocal];
            RDouble az = grid->GetFaceNormalZ()[iFacelocal];


            RDouble taux = du - du * ax * ax - dv * ax * ay - dw * ax * az;
            RDouble tauy = dv - du * ay * ax - dv * ay * ay - dw * ay * az;
            RDouble tauz = dw - du * az * ax - dv * az * ay - dw * az * az;
            RDouble tauwall = std::sqrt(taux * taux + tauy * tauy + tauz * tauz);
            if (yplus > this->yptr_)
            {
                uplus = std::log(this->elog_ * yplus) / this->cappa_;
                mub[iFacelocal] = visl[le] * yplus / uplus;
                Fwall = mub[iFacelocal] * tauwall / wd;
                if (wallFlag[le] < 1)
                {
                    gen[le] += Fwall * cmute / (this->cappa_ * wd);
                    tauWall[le] = tauwall;
                    wallFlag[le] = 1;
                }
                else
                {
                    if (tauwall > tauWall[le])
                    {
                        tauWall[le] = tauwall;
                    }
                    gen[le] += Fwall * cmute / (this->cappa_ * wd);
                    wallFlag[le] += 1;
                    gen[le] /= wallFlag[le];
                }
            }
            else
            {
                gen[le] = 0;
            }
        }
    }

    //£¡Correct:Wall Function.
    void PBBoundaryTurbKWall::Discrt(UnstructGrid *grid, string phiName, vector<int> *faceIndex, Data_Param *bcData, RDouble *diffCoe)
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

        for (std::size_t iFace= 0; iFace < faceIndex->size(); ++ iFace)
        {
            if (bcData)
            {
                if (bcData->IsExist(phiName, PHDOUBLE, 1))
                {
                    const int iFacelocal = (*faceIndex)[iFace];  //! The boundary element serial number is arranged in front of the element face.
                    const int le = grid->GetLeftCellOfFace()[iFacelocal];
                    //! The boundary value is ranked behind the value on the cell. When it is passed in, the lower offset corresponds to the boundary sequence number.
                    RDouble phibValue = phib[iFacelocal]; 
                    RDouble diffCoebValue = diffCoeb[iFacelocal];
                    RDouble area = grid->GetFaceArea()[iFacelocal];

                    RDouble fanx = grid->GetFaceNormalX()[iFacelocal] * area;
                    RDouble fany = grid->GetFaceNormalY()[iFacelocal] * area;
                    RDouble fanz = grid->GetFaceNormalZ()[iFacelocal] * area;
                    RDouble dx = grid->GetFaceCenterX()[iFacelocal] - grid->GetCellCenterX()[le];
                    RDouble dy = grid->GetFaceCenterY()[iFacelocal] - grid->GetCellCenterY()[le];
                    RDouble dz = grid->GetFaceCenterZ()[iFacelocal] - grid->GetCellCenterZ()[le];

                    RDouble area2 = diffCoebValue * (fanx * fanx + fany * fany + fanz * fanz) / (fanx * dx + fany * dy + fanz * dz);
                    RDouble ax = diffCoebValue * fanx - dx * area2;
                    RDouble ay = diffCoebValue * fany - dy * area2;
                    RDouble az = diffCoebValue * fanz - dz * area2;

                    apPhi[le] += area2;
                    scPhi[le] += area2 * phib[iFacelocal] + dphidx[le] * ax + dphidy[le] * ay + dphidz[le] * az;
                }
            }
        }
    }

    DEFINE_PBSOLVER_REGISTER_FUNCTION(PBBoundary, PBBoundaryTurbKWall, TURB_K_WALL);
} // PB_SPACE