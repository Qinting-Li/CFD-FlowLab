#include <iostream>

#include "PBBoundaryFlowSolidSurface.h"
#include "Geo_UnstructGrid.h"

using namespace std;

namespace PHSPACE
{
PBBoundaryFlowSolidSurface::PBBoundaryFlowSolidSurface(const string &name):PBBoundaryFlow(name)
{

}

PBBoundaryFlowSolidSurface::~PBBoundaryFlowSolidSurface()
{

}

void PBBoundaryFlowSolidSurface::CalcVelocityBoundaryCondition(UnstructGrid *grid, vector<int> *faceIndex)
{
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

    RDouble *mu   = reinterpret_cast<RDouble *>(grid->GetDataPtr("mu"));
    RDouble *rho  = reinterpret_cast<RDouble *>(grid->GetDataPtr("rho"));
    RDouble *u    = reinterpret_cast<RDouble *>(grid->GetDataPtr("U"));
    RDouble *v    = reinterpret_cast<RDouble *>(grid->GetDataPtr("V"));
    RDouble *w    = reinterpret_cast<RDouble *>(grid->GetDataPtr("W"));
    RDouble *dudx = reinterpret_cast<RDouble *>(grid->GetDataPtr("dudx"));
    RDouble *dudy = reinterpret_cast<RDouble *>(grid->GetDataPtr("dudy"));
    RDouble *dudz = reinterpret_cast<RDouble *>(grid->GetDataPtr("dudz"));
    RDouble *dvdx = reinterpret_cast<RDouble *>(grid->GetDataPtr("dvdx"));
    RDouble *dvdy = reinterpret_cast<RDouble *>(grid->GetDataPtr("dvdy"));
    RDouble *dvdz = reinterpret_cast<RDouble *>(grid->GetDataPtr("dvdz"));
    RDouble *dwdx = reinterpret_cast<RDouble *>(grid->GetDataPtr("dwdx"));
    RDouble *dwdy = reinterpret_cast<RDouble *>(grid->GetDataPtr("dwdy"));
    RDouble *dwdz = reinterpret_cast<RDouble *>(grid->GetDataPtr("dwdz"));

    RDouble *fluxb = reinterpret_cast<RDouble *>(grid->GetDataPtr("flux"));
    RDouble *apU = reinterpret_cast<RDouble *>(grid->GetDataPtr("apU"));
    RDouble *apV = reinterpret_cast<RDouble *>(grid->GetDataPtr("apV"));
    RDouble *apW = reinterpret_cast<RDouble *>(grid->GetDataPtr("apW"));
    RDouble *scU = reinterpret_cast<RDouble *>(grid->GetDataPtr("scU"));
    RDouble *scV = reinterpret_cast<RDouble *>(grid->GetDataPtr("scV"));
    RDouble *scW = reinterpret_cast<RDouble *>(grid->GetDataPtr("scW"));

    int nTotalCells = grid->GetNTotalCell();
    RDouble *mub = mu + nTotalCells;

    for (size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
    {
        int iFacelocal = (*faceIndex)[iFace];
        int le = leftCellOfFace[iFacelocal];
        int re = rightCellOfFace[iFacelocal];
        RDouble ub = u[re];
        RDouble vb = v[re];
        RDouble wb = w[re];
        RDouble muLocal = mub[iFacelocal];
        RDouble area = faceArea[iFacelocal];

        RDouble fanx = xfn[iFacelocal] * area;
        RDouble fany = yfn[iFacelocal] * area;
        RDouble fanz = zfn[iFacelocal] * area;
        RDouble dx = xfc[iFacelocal] - xcc[le];
        RDouble dy = yfc[iFacelocal] - ycc[le];
        RDouble dz = zfc[iFacelocal] - zcc[le];

        RDouble area2 = muLocal * (fanx * fanx + fany * fany + fanz * fanz) / (fanx * dx + fany * dy + fanz * dz);
        RDouble ax = muLocal * fanx - dx * area2;
        RDouble ay = muLocal * fany - dy * area2;
        RDouble az = muLocal * fanz - dz * area2;

        fluxb[iFacelocal] = 0.0;
        apU[le] = apU[le] + area2;
        apV[le] = apV[le] + area2;
        apW[le] = apW[le] + area2;
        scU[le] = scU[le] + area2 * ub + dudx[le] * ax + dudy[le] * ay + dudz[le] * az;
        scV[le] = scV[le] + area2 * vb + dvdx[le] * ax + dvdy[le] * ay + dvdz[le] * az;
        scW[le] = scW[le] + area2 * wb + dwdx[le] * ax + dwdy[le] * ay + dwdz[le] * az;
    }
}

void PBBoundaryFlowSolidSurface::UpdateVelocityBoundaryCondition(UnstructGrid *grid, vector<int> *faceIndex, Data_Param *bcData)
{
    RDouble ub, vb, wb;
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    RDouble *u = reinterpret_cast<RDouble *>(grid->GetDataPtr("U"));
    RDouble *v = reinterpret_cast<RDouble *>(grid->GetDataPtr("V"));
    RDouble *w = reinterpret_cast<RDouble *>(grid->GetDataPtr("W"));

    for (std::size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
    {
        int iFacelocal = (*faceIndex)[iFace];

        int re = rightCellOfFace[iFacelocal];

        if (bcData)
        {
            if (bcData->IsExist("flowU", PHDOUBLE, 1))
            {
                bcData->GetData("flowU", &ub, PHDOUBLE, 1);
            }
            else
            {
                ub = 0.0;
            }

            if (bcData->IsExist("flowV", PHDOUBLE, 1))
            {
                bcData->GetData("flowV", &vb, PHDOUBLE, 1);
            }
            else
            {
                vb = 0.0;
            }

            if (bcData->IsExist("flowW", PHDOUBLE, 1))
            {
                bcData->GetData("flowW", &wb, PHDOUBLE, 1);
            }
            else
            {
                wb = 0.0;
            }

            u[re] = RDouble(ub);
            v[re] = RDouble(vb);
            w[re] = RDouble(wb);
        }
    }
}

void PBBoundaryFlowSolidSurface::UpdatePressureBoundaryCondition(UnstructGrid *grid, vector<int> *faceIndex, Data_Param *bcData)
{
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    RDouble *xfc = grid->GetFaceCenterX();
    RDouble *yfc = grid->GetFaceCenterY();
    RDouble *zfc = grid->GetFaceCenterZ();

    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    RDouble *p = reinterpret_cast<RDouble *>(grid->GetDataPtr("P"));
    RDouble *dpdx = reinterpret_cast<RDouble *>(grid->GetDataPtr("dpdx"));
    RDouble *dpdy = reinterpret_cast<RDouble *>(grid->GetDataPtr("dpdy"));
    RDouble *dpdz = reinterpret_cast<RDouble *>(grid->GetDataPtr("dpdz"));
    if (GlobalDataBase::GetIntParaFromDB("bodyForceFlag"))
    {
        RDouble *bfU = reinterpret_cast<RDouble *>(grid->GetDataPtr("bfU"));
        RDouble *bfV = reinterpret_cast<RDouble *>(grid->GetDataPtr("bfV"));
        RDouble *bfW = reinterpret_cast<RDouble *>(grid->GetDataPtr("bfW"));
        RDouble *cellVolume = grid->GetCellVolume();

        for (size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
        {
            int iFacelocal = (*faceIndex)[iFace];
            int le = leftCellOfFace[iFacelocal];
            int re = rightCellOfFace[iFacelocal];
            RDouble faceCenterX = xfc[iFacelocal];
            RDouble faceCenterY = yfc[iFacelocal];
            RDouble faceCenterZ = zfc[iFacelocal];
            RDouble cellCenterXle = xcc[le];
            RDouble cellCenterYle = ycc[le];
            RDouble cellCenterZle = zcc[le];

            RDouble dx = faceCenterX - cellCenterXle;
            RDouble dy = faceCenterY - cellCenterYle;
            RDouble dz = faceCenterZ - cellCenterZle;
            p[re] = p[le] + bfU[le] * dx + bfV[le] * dy + bfW[le] * dz;
        }
    }
    else
    {
        for (size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
        {
            int iFacelocal = (*faceIndex)[iFace];
            int le = leftCellOfFace[iFacelocal];
            int re = rightCellOfFace[iFacelocal];
            RDouble faceCenterX = xfc[iFacelocal];
            RDouble faceCenterY = yfc[iFacelocal];
            RDouble faceCenterZ = zfc[iFacelocal];
            RDouble cellCenterXle = xcc[le];
            RDouble cellCenterYle = ycc[le];
            RDouble cellCenterZle = zcc[le];
            RDouble dx = faceCenterX - cellCenterXle;
            RDouble dy = faceCenterY - cellCenterYle;
            RDouble dz = faceCenterZ - cellCenterZle;
            p[re] = p[le]+ dpdx[le]*dx+ dpdy[le] * dy+ dpdz[le] * dz;
        }
    }
}

void PBBoundaryFlowSolidSurface::CalcPressureCorrectionBoundaryCondition(UnstructGrid *grid, vector<int> *faceIndex, Data_Param *bcData)
{

}

void PBBoundaryFlowSolidSurface::UpdatePressureCorrectionBoundaryCondition(UnstructGrid *grid, std::vector<int> *faceIndex, Data_Param *bcData)
{
    int* leftCellOfFace = grid->GetLeftCellOfFace();
    int* rightCellOfFace = grid->GetRightCellOfFace();

    RDouble *xfc = grid->GetFaceCenterX();
    RDouble *yfc = grid->GetFaceCenterY();
    RDouble *zfc = grid->GetFaceCenterZ();

    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    RDouble* pc= reinterpret_cast<RDouble*>((grid->GetDataPtr("pc")));
    RDouble* dpcdx = reinterpret_cast<RDouble*>(grid->GetDataPtr("dpcdx"));
    RDouble* dpcdy = reinterpret_cast<RDouble*>(grid->GetDataPtr("dpcdy"));
    RDouble* dpcdz = reinterpret_cast<RDouble*>(grid->GetDataPtr("dpcdz"));
    for (size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
    {
        int iFacelocal = (*faceIndex)[iFace];
        int le = leftCellOfFace[iFacelocal];
        int re = rightCellOfFace[iFacelocal];
        RDouble faceCenterX = xfc[iFacelocal];
        RDouble faceCenterY = yfc[iFacelocal];
        RDouble faceCenterZ = zfc[iFacelocal];
        RDouble cellCenterXle = xcc[le];
        RDouble cellCenterYle = ycc[le];
        RDouble cellCenterZle = zcc[le];
        RDouble dx = faceCenterX - cellCenterXle;
        RDouble dy = faceCenterY - cellCenterYle;
        RDouble dz = faceCenterZ - cellCenterZle;
        pc[re] = pc[le] + dpcdx[le] * dx + dpcdy[le] * dy + dpcdz[le] * dz;
    }
}

void PBBoundaryFlowSolidSurface::UpdateVelocityCorrectionBoundary(UnstructGrid *grid, vector<int> *faceIndex, Data_Param *bcData)
{

}

DEFINE_PBSOLVER_REGISTER_FUNCTION(PBBoundaryFlow, PBBoundaryFlowSolidSurface, FLOW_SOLID_SURFACE);
}
