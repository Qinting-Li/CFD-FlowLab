#include <iostream>
#include "PBDiffusionNonorthogonal.h"

using namespace std;

namespace PHSPACE
{
PBDiffusionNonorthogonal::PBDiffusionNonorthogonal(const string &name):PBDiffusion(name)
{

}

PBDiffusionNonorthogonal::~PBDiffusionNonorthogonal()
{

}

void PBDiffusionNonorthogonal::Discrt(UnstructGrid *grid, string phiName, RDouble *diffCoe, RDouble *dphidx, RDouble *dphidy, RDouble *dphidz,
                                      RDouble *ajLeft, RDouble *ajRight,RDouble *scPhi, RDouble *dvdx, RDouble *dvdy, RDouble *dvdz,
                                      RDouble *scV, RDouble *dwdx, RDouble *dwdy, RDouble *dwdz, RDouble *scW)
{
    int faceBegin = grid->GetNBoundFace();
    int faceEnd = grid->GetNTotalFace();
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    RDouble *faceArea = grid->GetFaceArea();
    RDouble *xfn = grid->GetFaceNormalX();
    RDouble *yfn = grid->GetFaceNormalY();
    RDouble *zfn = grid->GetFaceNormalZ();
    RDouble *xfc = grid->GetFaceCenterX();
    RDouble *yfc = grid->GetFaceCenterY();
    RDouble *zfc = grid->GetFaceCenterZ();
    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    RDouble *deltaL = reinterpret_cast<RDouble *>(grid->GetDataPtr("deltL"));
    RDouble *deltaR = reinterpret_cast<RDouble *>(grid->GetDataPtr("deltR"));

    if (dvdx == NULL)
    {
        UnstructGrid *grid = UnstructGridCast(GetGrid(0));
        RDouble *cp;
        RDouble *sp;
        RDouble *u;
        RDouble *v;
        RDouble *w;
        if (phiName == "enthalpy")
        {
            cp = reinterpret_cast<RDouble *>(grid->GetDataPtr("cp"));
            sp = reinterpret_cast<RDouble *>(grid->GetDataPtr("ap"));
            u = reinterpret_cast<RDouble *>(grid->GetDataPtr("U"));
            v = reinterpret_cast<RDouble *>(grid->GetDataPtr("V"));
            w = reinterpret_cast<RDouble *>(grid->GetDataPtr("W"));
        }

        for (int iFace = faceBegin; iFace < faceEnd; ++ iFace)
        {
            int le = leftCellOfFace[iFace];
            int re = rightCellOfFace[iFace];

            RDouble ax = faceArea[iFace] * xfn[iFace];
            RDouble ay = faceArea[iFace] * yfn[iFace];
            RDouble az = faceArea[iFace] * zfn[iFace];

            RDouble ex = xcc[re] - xcc[le];
            RDouble ey = ycc[re] - ycc[le];
            RDouble ez = zcc[re] - zcc[le];

            RDouble wfL = deltaL[iFace];
            RDouble wfR = deltaR[iFace];

            RDouble diffCoeF = wfL * diffCoe[le] + wfR * diffCoe[re];

            RDouble aLR = diffCoeF * (ax * ax + ay * ay + az * az) / (ax * ex + ay * ey + az * ez);

            if (phiName == "enthalpy")
            {
                ajLeft[iFace - faceBegin] += aLR / cp[rightCellOfFace[iFace]];
                ajRight[iFace - faceBegin] += aLR / cp[leftCellOfFace[iFace]];

                sp[le] = sp[le] + ajRight[iFace - faceBegin] - ajLeft[iFace - faceBegin];
                sp[re] = sp[re] + ajLeft[iFace - faceBegin] - ajRight[iFace - faceBegin];

                RDouble usq = (u[le] * u[le] + v[le] * v[le] + w[le] * w[le]) / 2;
                scPhi[le] += aLR * usq / cp[le];
            }
            else
            {
                ajLeft[iFace - faceBegin] += aLR;
                ajRight[iFace - faceBegin] += aLR;
            }

            ax = diffCoeF * ax - ex * aLR;
            ay = diffCoeF * ay - ey * aLR;
            az = diffCoeF * az - ez * aLR;

            RDouble sLR = ax * (wfL * dphidx[le] + wfR * dphidx[re]) +
                ay * (wfL * dphidy[le] + wfR * dphidy[re]) +
                az * (wfL * dphidz[le] + wfR * dphidz[re]);

            scPhi[le] += sLR;
            scPhi[re] -= sLR;
        }
    }
    else
    {
        for (int iFace = faceBegin; iFace < faceEnd; ++ iFace)
        {
            int le = leftCellOfFace[iFace];
            int re = rightCellOfFace[iFace];

            RDouble ax = faceArea[iFace] * xfn[iFace];
            RDouble ay = faceArea[iFace] * yfn[iFace];
            RDouble az = faceArea[iFace] * zfn[iFace];

            RDouble ex = xcc[re] - xcc[le];
            RDouble ey = ycc[re] - ycc[le];
            RDouble ez = zcc[re] - zcc[le];

            RDouble wfL = deltaL[iFace];
            RDouble wfR = deltaR[iFace];

            RDouble diffCoeF = wfL * diffCoe[le] + wfR * diffCoe[re];

            RDouble aLR = diffCoeF * (ax * ax + ay * ay + az * az) / (ax * ex + ay * ey + az * ez);

            ajLeft[iFace - faceBegin] += aLR;
            ajRight[iFace - faceBegin] += aLR;

            ax = diffCoeF * ax - ex * aLR;
            ay = diffCoeF * ay - ey * aLR;
            az = diffCoeF * az - ez * aLR;

            RDouble sLR = ax * (wfL * dphidx[le] + wfR * dphidx[re]) +
                ay * (wfL * dphidy[le] + wfR * dphidy[re]) +
                az * (wfL * dphidz[le] + wfR * dphidz[re]);

            scPhi[le] += sLR;
            scPhi[re] -= sLR;

            sLR = ax * (wfL * dvdx[le] + wfR * dvdx[re]) +
                ay * (wfL * dvdy[le] + wfR * dvdy[re]) +
                az * (wfL * dvdz[le] + wfR * dvdz[re]);

            scV[le] += sLR;
            scV[re] -= sLR;

            sLR = ax * (wfL * dwdx[le] + wfR * dwdx[re]) +
                ay * (wfL * dwdy[le] + wfR * dwdy[re]) +
                az * (wfL * dwdz[le] + wfR * dwdz[re]);

            scW[le] += sLR;
            scW[re] -= sLR;
        }
    }
}

DEFINE_PBSOLVER_REGISTER_FUNCTION(PBDiffusion, PBDiffusionNonorthogonal, DIFFUSION_NON_ORTHOGONAL);
}