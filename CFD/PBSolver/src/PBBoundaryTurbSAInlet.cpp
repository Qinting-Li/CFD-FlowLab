#include "PBBoundaryTurbSAInlet.h"

using namespace std;

namespace PHSPACE
{
PBBoundaryTurbSAInlet::PBBoundaryTurbSAInlet(const string &name):PBBoundaryTurbDefault(name)
{

}

PBBoundaryTurbSAInlet::~PBBoundaryTurbSAInlet()
{

}

void PBBoundaryTurbSAInlet::Update(UnstructGrid *grid, string phiName, vector<int> *faceIndex, Data_Param *bcData)
{
    RDouble kb = 0.0;
    RDouble *phi = reinterpret_cast<RDouble *>(grid->GetDataPtr("kinetic"));
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int nTotalCell = grid->GetNTotalCell();
    RDouble *phib = phi + nTotalCell;
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

void PBBoundaryTurbSAInlet::Discrt(UnstructGrid *grid, string phiName, vector<int> *faceIndex, Data_Param *bcData, RDouble *diffCoe)
{
    int nTotalCell = grid->GetNTotalCell();
    int *leftCellOfFace = grid->GetLeftCellOfFace();

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

    RDouble *phi = reinterpret_cast<RDouble *>(grid->GetDataPtr(phiName));
    RDouble *phib = phi + nTotalCell;
    RDouble *diffCoeb = diffCoe + nTotalCell;
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
            if (bcData->IsExist(phiName, PHDOUBLE, 1))
            {
                int iFacelocal = (*faceIndex)[iFace];
                int le = leftCellOfFace[iFacelocal];
                RDouble phibValue = phib[iFacelocal];
                RDouble diffCoebValue = diffCoeb[iFacelocal];
                RDouble area = faceArea[iFacelocal];

                RDouble fanx = xfn[iFacelocal] * area;
                RDouble fany = yfn[iFacelocal] * area;
                RDouble fanz = zfn[iFacelocal] * area;
                RDouble dx = xfc[iFacelocal] - xcc[le];
                RDouble dy = yfc[iFacelocal] - ycc[le];
                RDouble dz = zfc[iFacelocal] - zcc[le];

                RDouble area2 = diffCoebValue * (fanx * fanx + fany * fany + fanz * fanz) / (fanx * dx + fany * dy + fanz * dz);
                RDouble ax = diffCoebValue * fanx - dx * area2;
                RDouble ay = diffCoebValue * fany - dy * area2;
                RDouble az = diffCoebValue * fanz - dz * area2;

                apPhi[le] += area2;
                scPhi[le] += area2 * phib[iFacelocal] + dphidx[le] * ax + dphidy[le] * ay + dphidz[le] * az;

                if (fluxb[iFacelocal] < 0)
                {
                    scPhi[le] -= fluxb[iFacelocal] * phib[iFacelocal];
                }
            }
            else
            {
                int iFacelocal = (*faceIndex)[iFace];
                int le = leftCellOfFace[iFacelocal];

                if (fluxb[iFacelocal] < 0)
                {
                    scPhi[le] -= fluxb[iFacelocal] * phib[iFacelocal];
                }
            }
        }
    }
}

DEFINE_PBSOLVER_REGISTER_FUNCTION(PBBoundary, PBBoundaryTurbSAInlet, TURB_SA_INLET);
} //! PB_SPACE