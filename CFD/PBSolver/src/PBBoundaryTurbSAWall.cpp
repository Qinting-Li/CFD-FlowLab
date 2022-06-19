#include <cmath>
#include "PBBoundaryTurbSAWall.h"

using namespace std;

namespace PHSPACE
{
PBBoundaryTurbSAWall::PBBoundaryTurbSAWall(const string &name):PBBoundaryScalarDiffusion(name)
{

}

PBBoundaryTurbSAWall::~PBBoundaryTurbSAWall()
{

}

void PBBoundaryTurbSAWall::Discrt(UnstructGrid *grid, string phiName, vector<int> *faceIndex, Data_Param *bcData, RDouble *diffCoe)
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
    RDouble *visl = reinterpret_cast<RDouble *>(grid->GetDataPtr("visl"));
    RDouble *apPhi = reinterpret_cast<RDouble *>(grid->GetDataPtr("ap"));
    RDouble *scPhi = reinterpret_cast<RDouble *>(grid->GetDataPtr("sc"));
    RDouble *dphidx = reinterpret_cast<RDouble *>(grid->GetDataPtr("dphidx"));
    RDouble *dphidy = reinterpret_cast<RDouble *>(grid->GetDataPtr("dphidy"));
    RDouble *dphidz = reinterpret_cast<RDouble *>(grid->GetDataPtr("dphidz"));

    for (std::size_t iFace = 0; iFace < faceIndex->size(); ++ iFace)
    {
        if (bcData)
        {

            int iFacelocal = (*faceIndex)[iFace]; //! The boundary element serial number is arranged in front of the element face.
            int le = leftCellOfFace[iFacelocal];
            //! The boundary value is ranked behind the value on the cell. When it is passed in, the lower offset
            //! corresponds to the boundary sequence number.
            RDouble phibValue = 0.0; 
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
            scPhi[le] += area2 * phibValue + dphidx[le] * ax + dphidy[le] * ay + dphidz[le] * az;

        }
    }
}

DEFINE_PBSOLVER_REGISTER_FUNCTION(PBBoundary, PBBoundaryTurbSAWall, TURB_SA_WALL);
} //! PB_SPACE