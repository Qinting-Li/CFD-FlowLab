#include "PBBoundaryTurbKPressureOutlet.h"

using namespace std;

namespace PHSPACE
{
PBBoundaryTurbKPressureOutlet::PBBoundaryTurbKPressureOutlet(const string &name):PBBoundaryTurbDefault(name)
{

}

PBBoundaryTurbKPressureOutlet::~PBBoundaryTurbKPressureOutlet()
{

}

void PBBoundaryTurbKPressureOutlet::Update(UnstructGrid *grid, string phiName, vector<int> *faceIndex, Data_Param *bcData)
{
    int nTotalCells = grid->GetNTotalCell();
    int *leftCellOfFace = grid->GetLeftCellOfFace();

    RDouble *xfc = grid->GetFaceCenterX();
    RDouble *yfc = grid->GetFaceCenterY();
    RDouble *zfc = grid->GetFaceCenterZ();
    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    RDouble kb = 0.0;
    RDouble *flux_ = reinterpret_cast<RDouble *>(grid->GetDataPtr("flux"));
    RDouble *phi = reinterpret_cast<RDouble *>(grid->GetDataPtr("kinetic"));
    RDouble *dphidx = reinterpret_cast<RDouble *>(grid->GetDataPtr("dphidx"));
    RDouble *dphidy = reinterpret_cast<RDouble *>(grid->GetDataPtr("dphidy"));
    RDouble *dphidz = reinterpret_cast<RDouble *>(grid->GetDataPtr("dphidz"));
    RDouble *phib = phi + nTotalCells;
    for (std::size_t iFace = 0; iFace < faceIndex->size(); ++ iFace)
    {
        int iFacelocal = (*faceIndex)[iFace];
        int le = leftCellOfFace[iFacelocal];
        RDouble faceCenterX = xfc[iFacelocal];
        RDouble faceCenterY = yfc[iFacelocal];
        RDouble faceCenterZ = zfc[iFacelocal];
        RDouble cellCenterXle = xcc[le];
        RDouble cellCenterYle = ycc[le];
        RDouble cellCenterZle = zcc[le];

        RDouble dx = faceCenterX - cellCenterXle;
        RDouble dy = faceCenterY - cellCenterYle;
        RDouble dz = faceCenterZ - cellCenterZle;
        if (flux_[iFacelocal] < 0.0)
        {
            if (bcData)
            {
                if (bcData->IsExist(phiName, PHDOUBLE, 1))
                {
                    bcData->GetData(phiName, &kb, PHDOUBLE, 1);
                    phib[iFacelocal] = RDouble(kb);
                }
                else
                {
                    std::cout << " Please assign a value " << std::endl;
                }
            }
        }
        else
        {
            phib[iFacelocal] = phi[le] + dphidx[le] * dx + dphidy[le] * dy + dphidz[le] * dz;
        }
    }
}

void PBBoundaryTurbKPressureOutlet::Discrt(UnstructGrid *grid, string phiName, vector<int> *faceIndex, Data_Param *bcData, RDouble *diffCoe)
{
    int nTotalCells = grid->GetNTotalCell();
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
            if (bcData->IsExist(phiName, PHDOUBLE, 1))
            {
                int iFacelocal = (*faceIndex)[iFace]; //! The boundary element serial number is arranged in front of the element face.
                int le = leftCellOfFace[iFacelocal];
                //! The boundary value is ranked behind the value on the cell. When it is passed in, the lower offset
                //! corresponds to the boundary sequence number.
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
                    apPhi[le] -= fluxb[iFacelocal];
                    scPhi[le] -= fluxb[iFacelocal] * phib[iFacelocal];
                }
            }
            else
            {
                int iFacelocal = (*faceIndex)[iFace]; //! The boundary element serial number is arranged in front of the element face.
                int le = leftCellOfFace[iFacelocal];

                if (fluxb[iFacelocal] < 0)
                {
                    apPhi[le] -= fluxb[iFacelocal];
                    scPhi[le] -= fluxb[iFacelocal] * phib[iFacelocal];
                }
            }
        }
    }
}

DEFINE_PBSOLVER_REGISTER_FUNCTION(PBBoundary, PBBoundaryTurbKPressureOutlet, TURB_K_PRESSUREOUTLET);
} //! PB_SPACE