#include <iostream>

#include "PBBoundaryScalarDiffusion.h"

namespace PHSPACE
{
    PBBoundaryScalarDiffusion::PBBoundaryScalarDiffusion(const string &name):PBBoundaryScalar(name)
    {

    }

    PBBoundaryScalarDiffusion::~PBBoundaryScalarDiffusion()
    {

    }

    void PBBoundaryScalarDiffusion::Discrt(UnstructGrid *grid, string phiName, vector<int> *faceIndex, Data_Param *bcData, RDouble *diffCoe)
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
        for (std::size_t iFace = 0; iFace < faceIndex->size(); ++ iFace)
        {
            int iFacelocal = (*faceIndex)[iFace]; //! The boundary element serial number is arranged in front of the element face.
            int le = grid->GetLeftCellOfFace()[iFacelocal];

            RDouble phibValue = phib[iFacelocal]; 
            RDouble diffCoebValue = diffCoeb[iFacelocal];
            RDouble area = grid->GetFaceArea()[iFacelocal];

            //! The boundary value is ranked behind the value on the cell. When it is passed in, the lower offset corresponds to the boundary sequence number.
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

    DEFINE_PBSOLVER_REGISTER_FUNCTION(PBBoundary, PBBoundaryScalarDiffusion, SCALAR_DIFFUSION);
}