#include "PBBoundaryEnergyOutlet.h"

using namespace std;

namespace PHSPACE
{
PBBoundaryEnergyOutlet::PBBoundaryEnergyOutlet(const string &name):PBBoundaryEnergyDefault(name)
{

}

PBBoundaryEnergyOutlet::~PBBoundaryEnergyOutlet()
{

}

void PBBoundaryEnergyOutlet::Discrt(UnstructGrid *grid, string phiName, vector<int> *faceIndex, Data_Param *bcData, RDouble *diffCoe)
{
    int nTotalCell = grid->GetNTotalCell();
    int *leftCellOfFace  = grid->GetLeftCellOfFace();

    RDouble *xfn  = grid->GetFaceNormalX();
    RDouble *yfn  = grid->GetFaceNormalY();
    RDouble *zfn  = grid->GetFaceNormalZ();
    RDouble *faceArea = grid->GetFaceArea();

    RDouble *xcc  = grid->GetCellCenterX();
    RDouble *ycc  = grid->GetCellCenterY();
    RDouble *zcc  = grid->GetCellCenterZ();

    RDouble *xfc  = grid->GetFaceCenterX();
    RDouble *yfc  = grid->GetFaceCenterY();
    RDouble *zfc  = grid->GetFaceCenterZ();

    RDouble *phi = reinterpret_cast<RDouble *>(grid->GetDataPtr("T"));
    RDouble *phib = phi + nTotalCell;
    RDouble *diffCoeb = diffCoe + nTotalCell;
    RDouble *cp = reinterpret_cast<RDouble *>(grid->GetDataPtr("cp"));
    RDouble *u = reinterpret_cast<RDouble *>(grid->GetDataPtr("U"));
    RDouble *v = reinterpret_cast<RDouble *>(grid->GetDataPtr("V"));
    RDouble *w = reinterpret_cast<RDouble *>(grid->GetDataPtr("W"));
    RDouble *apPhi = reinterpret_cast<RDouble *>(grid->GetDataPtr("ap"));
    RDouble *scPhi = reinterpret_cast<RDouble *>(grid->GetDataPtr("sc"));
    RDouble *dphidx = reinterpret_cast<RDouble *>(grid->GetDataPtr("dphidx"));
    RDouble *dphidy = reinterpret_cast<RDouble *>(grid->GetDataPtr("dphidy"));
    RDouble *dphidz = reinterpret_cast<RDouble *>(grid->GetDataPtr("dphidz"));
    RDouble *fluxb = reinterpret_cast<RDouble *>(grid->GetDataPtr("flux"));

    for (std::size_t iFace = 0; iFace < faceIndex->size(); ++ iFace)
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

        apPhi[le] += area2 / cp[le];
        scPhi[le] += area2 * phib[iFacelocal] + dphidx[le] * ax + dphidy[le] * ay + dphidz[le] * az;

        RDouble usq = (u[le] * u[le] + v[le] * v[le] + w[le] * w[le]) / 2;
        scPhi[le] += area2 * usq / cp[le];

        RDouble fluxbValue = -MIN(fluxb[iFacelocal], 0);

        scPhi[le] += fluxbValue * phib[iFacelocal];
    }
}

void PBBoundaryEnergyOutlet::Update(UnstructGrid *grid, std::string phiName,std::vector<int> *faceIndex, Data_Param *bcData)
{
    int nTotalCell = grid->GetNTotalCell();
    int *leftCellOfFace  = grid->GetLeftCellOfFace();
    RDouble *phi = reinterpret_cast<RDouble *>(grid->GetDataPtr("T"));
    RDouble *phib = phi + nTotalCell;
    RDouble *flux_ = reinterpret_cast<RDouble *>(grid->GetDataPtr("flux"));


    for (std::size_t iFace = 0; iFace < faceIndex->size(); ++ iFace)
    {
        RDouble spBound;

        int iFacelocal = (*faceIndex)[iFace];
        int le = leftCellOfFace[iFacelocal];

        if (flux_[iFacelocal] < 0.0)
        {
            if (bcData)
            {
                if (bcData->IsExist("T", PHDOUBLE, 1))
                {
                    bcData->GetData("T", &spBound, PHDOUBLE, 1);
                    phib[iFacelocal] = RDouble(spBound);
                }
                else
                {
                    std::cout << " Please assign a value " << std::endl;
                }
            }
        }
        else
        {
            phib[iFacelocal] = phi[le];
        }
    }
}


DEFINE_PBSOLVER_REGISTER_FUNCTION(PBBoundary, PBBoundaryEnergyOutlet, ENERGY_OUTLET);
}