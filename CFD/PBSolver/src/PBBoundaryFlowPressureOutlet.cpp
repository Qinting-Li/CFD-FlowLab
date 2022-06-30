#include <iostream>
#include <algorithm>

#include "PBBoundaryFlowPressureOutlet.h" 
#include "Geo_UnstructGrid.h"

using namespace std;

namespace PHSPACE
{
PBBoundaryFlowPressureOutlet::PBBoundaryFlowPressureOutlet(const string &name):PBBoundaryFlow(name)
{

}

PBBoundaryFlowPressureOutlet::~PBBoundaryFlowPressureOutlet()
{

}

void PBBoundaryFlowPressureOutlet::CalcVelocityBoundaryCondition(UnstructGrid *grid, vector<int> *faceIndex)
{
    int *leftCellOfFace  = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    RDouble *faceArea    = grid->GetFaceArea();
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
    RDouble *dudx = reinterpret_cast<RDouble *>(grid->GetDataPtr("dudx"));
    RDouble *dudy = reinterpret_cast<RDouble *>(grid->GetDataPtr("dudy"));
    RDouble *dudz = reinterpret_cast<RDouble *>(grid->GetDataPtr("dudz"));
    RDouble *dvdx = reinterpret_cast<RDouble *>(grid->GetDataPtr("dvdx"));
    RDouble *dvdy = reinterpret_cast<RDouble *>(grid->GetDataPtr("dvdy"));
    RDouble *dvdz = reinterpret_cast<RDouble *>(grid->GetDataPtr("dvdz"));
    RDouble *dwdx = reinterpret_cast<RDouble *>(grid->GetDataPtr("dwdx"));
    RDouble *dwdy = reinterpret_cast<RDouble *>(grid->GetDataPtr("dwdy"));
    RDouble *dwdz = reinterpret_cast<RDouble *>(grid->GetDataPtr("dwdz"));

    RDouble *u     = reinterpret_cast<RDouble *>(grid->GetDataPtr("U"));
    RDouble *v     = reinterpret_cast<RDouble *>(grid->GetDataPtr("V"));
    RDouble *w     = reinterpret_cast<RDouble *>(grid->GetDataPtr("W"));
    RDouble *fluxb = reinterpret_cast<RDouble *>(grid->GetDataPtr("flux"));
    RDouble *apU   = reinterpret_cast<RDouble *>(grid->GetDataPtr("apU"));
    RDouble *apV   = reinterpret_cast<RDouble *>(grid->GetDataPtr("apV"));
    RDouble *apW   = reinterpret_cast<RDouble *>(grid->GetDataPtr("apW"));
    RDouble *scU   = reinterpret_cast<RDouble *>(grid->GetDataPtr("scU"));
    RDouble *scV   = reinterpret_cast<RDouble *>(grid->GetDataPtr("scV"));
    RDouble *scW   = reinterpret_cast<RDouble *>(grid->GetDataPtr("scW"));

    for (std::size_t iFace = 0; iFace < faceIndex->size(); ++ iFace)
    {
        int iFacelocal = (*faceIndex)[iFace];
        int le = leftCellOfFace[iFacelocal];
        int re = rightCellOfFace[iFacelocal];

        RDouble ub = u[re];
        RDouble vb = v[re];
        RDouble wb = w[re];
        RDouble muLocal = mu[le];
        RDouble area = faceArea[iFacelocal];
        RDouble fanx = xfn[iFacelocal] * area;
        RDouble fany = yfn[iFacelocal] * area;
        RDouble fanz = zfn[iFacelocal] * area;
        RDouble dx = xfc[iFacelocal] - xcc[le];
        RDouble dy = yfc[iFacelocal] - ycc[le];
        RDouble dz = zfc[iFacelocal] - zcc[le];
        RDouble fluxbValue = fluxb[iFacelocal];
        RDouble area2 = muLocal * (fanx * fanx + fany * fany + fanz * fanz) /
            (fanx * dx + fany * dy + fanz * dz);
        RDouble ax = muLocal * fanx - dx * area2;
        RDouble ay = muLocal * fany - dy * area2;
        RDouble az = muLocal * fanz - dz * area2;

        apU[le] = apU[le] + area2;
        apV[le] = apV[le] + area2;
        apW[le] = apW[le] + area2;
        scU[le] = scU[le] + area2 * ub + dudx[le] * ax + dudy[le] * ay + dudz[le] * az;
        scV[le] = scV[le] + area2 * vb + dvdx[le] * ax + dvdy[le] * ay + dvdz[le] * az;
        scW[le] = scW[le] + area2 * wb + dwdx[le] * ax + dwdy[le] * ay + dwdz[le] * az;

        if (fluxbValue < 0)
        {
            RDouble a12 = MAX(fluxbValue, 0.0);
            RDouble a21 = a12 - fluxbValue;

            apU[le] = apU[le] + a21;
            apV[le] = apV[le] + a21;
            apW[le] = apW[le] + a21;
            scU[le] = scU[le] + a21 * ub;
            scV[le] = scV[le] + a21 * vb;
            scW[le] = scW[le] + a21 * wb;
        }
    }
}

void PBBoundaryFlowPressureOutlet::UpdatePressureBoundaryCondition(UnstructGrid *grid, vector<int> *faceIndex, Data_Param *bcData)
{
    RDouble pb;
    int *rightCellOfFace = grid->GetRightCellOfFace();
    RDouble *p = reinterpret_cast<RDouble *>(grid->GetDataPtr("P"));
    for (std::size_t iFace = 0; iFace < faceIndex->size(); ++ iFace)
    {
        const int iFacelocal = (*faceIndex)[iFace];
        int re = rightCellOfFace[iFacelocal];

        if (bcData)
        {
            if (bcData->IsExist("flowP", PHDOUBLE, 1))
            {
                bcData->GetData("flowP", &pb, PHDOUBLE, 1);
            }
            else
            {
                pb = 0.0;
            }

            p[re] = pb;
        }
    }
}

void PBBoundaryFlowPressureOutlet::CalcPressureCorrectionBoundaryCondition(UnstructGrid *grid, vector<int> *faceIndex, Data_Param *bcData)
{
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    RDouble *faceArea    = grid->GetFaceArea();
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
    RDouble *dunb = reinterpret_cast<RDouble *>(grid->GetDataPtr("dun"));
    RDouble *dpdx = reinterpret_cast<RDouble *>(grid->GetDataPtr("dpdx"));
    RDouble *dpdy = reinterpret_cast<RDouble *>(grid->GetDataPtr("dpdy"));
    RDouble *dpdz = reinterpret_cast<RDouble *>(grid->GetDataPtr("dpdz"));
    RDouble *hu   = reinterpret_cast<RDouble *>(grid->GetDataPtr("hu"));
    RDouble *hv   = reinterpret_cast<RDouble *>(grid->GetDataPtr("hv"));
    RDouble *hw   = reinterpret_cast<RDouble *>(grid->GetDataPtr("hw"));
    RDouble *apU  = reinterpret_cast<RDouble *>(grid->GetDataPtr("apU"));
    RDouble *apV  = reinterpret_cast<RDouble *>(grid->GetDataPtr("apV"));
    RDouble *apW  = reinterpret_cast<RDouble *>(grid->GetDataPtr("apW"));

    RDouble urfFlux = this->ObjGetDoubleParaFromDB("urfFlux", 0.8);

    RDouble *u = reinterpret_cast<RDouble *>(grid->GetDataPtr("U"));
    RDouble *v = reinterpret_cast<RDouble *>(grid->GetDataPtr("V"));
    RDouble *w = reinterpret_cast<RDouble *>(grid->GetDataPtr("W"));
    RDouble *p = reinterpret_cast<RDouble *>(grid->GetDataPtr("P"));
    RDouble *fluxb = reinterpret_cast<RDouble *>(grid->GetDataPtr("flux"));
    RDouble *ap = reinterpret_cast<RDouble *>(grid->GetDataPtr("ap"));
    RDouble *sc = reinterpret_cast<RDouble *>(grid->GetDataPtr("sc"));

    for (std::size_t iFace = 0; iFace < faceIndex->size(); ++ iFace)
    {
        int iFacelocal = (*faceIndex)[iFace];
        int le = leftCellOfFace[iFacelocal];
        int re = rightCellOfFace[iFacelocal];

        RDouble rurf = urfFlux / (1 + urfFlux);
        RDouble muLocal = mu[le];
        RDouble area = faceArea[iFacelocal];

        RDouble fanx = xfn[iFacelocal] * area;
        RDouble fany = yfn[iFacelocal] * area;
        RDouble fanz = zfn[iFacelocal] * area;
        RDouble dx = xfc[iFacelocal] - xcc[le];
        RDouble dy = yfc[iFacelocal] - ycc[le];
        RDouble dz = zfc[iFacelocal] - zcc[le];

        RDouble ad = fanx * dx + fany * dy + fanz * dz;

        RDouble ub = u[le] + apU[le] * fanx / ad *
            ( p[re] - p[le] - dpdx[le] * dx - dpdy[le] * dy - dpdz[le] * dz );
        RDouble vb = v[le] + apV[le] * fany / ad *
            ( p[re] - p[le] - dpdx[le] * dx - dpdy[le] * dy - dpdz[le] * dz );
        RDouble wb = w[le] + apW[le] * fanz / ad *
            ( p[re] - p[le] - dpdx[le] * dx - dpdy[le] * dy - dpdz[le] * dz );

        RDouble aj = rho[le] * ( hu[le] * fanx * fanx + hv[le] * fany * fany + hw[le] * fanz * fanz ) / ad;

        RDouble fluxbValue = rho[le] * ((ub * fanx + vb * fany + wb * fanz)+ rurf * dunb[iFacelocal]);
        fluxb[iFacelocal] = fluxbValue;

        u[re] = ub;
        v[re] = vb;
        w[re] = wb;

        ap[le] += aj;
        sc[le] -= fluxbValue;
    }
}

void PBBoundaryFlowPressureOutlet::UpdateVelocityCorrectionBoundary(UnstructGrid *grid, vector<int> *faceIndex, Data_Param *bcData)
{
    RDouble un, fluxbValue;

    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    RDouble *faceArea = grid->GetFaceArea();
    RDouble *xfn = grid->GetFaceNormalX();
    RDouble *yfn = grid->GetFaceNormalY();
    RDouble *zfn = grid->GetFaceNormalZ();

    RDouble *dpcdx = reinterpret_cast<RDouble *>(grid->GetDataPtr("dpcdx"));
    RDouble *dpcdy = reinterpret_cast<RDouble *>(grid->GetDataPtr("dpcdy"));
    RDouble *dpcdz = reinterpret_cast<RDouble *>(grid->GetDataPtr("dpcdz"));
    RDouble *hu = reinterpret_cast<RDouble *>(grid->GetDataPtr("hu"));
    RDouble *hv = reinterpret_cast<RDouble *>(grid->GetDataPtr("hv"));
    RDouble *hw = reinterpret_cast<RDouble *>(grid->GetDataPtr("hw"));

    RDouble *rho = reinterpret_cast<RDouble *>(grid->GetDataPtr("rho"));
    RDouble *u = reinterpret_cast<RDouble *>(grid->GetDataPtr("U"));
    RDouble *v = reinterpret_cast<RDouble *>(grid->GetDataPtr("V"));
    RDouble *w = reinterpret_cast<RDouble *>(grid->GetDataPtr("W"));
    RDouble *dunb = reinterpret_cast<RDouble *>(grid->GetDataPtr("dun"));
    RDouble *fluxb = reinterpret_cast<RDouble *>(grid->GetDataPtr("flux"));

    for (std::size_t iFace = 0; iFace < faceIndex->size(); ++ iFace)
    {
        int iFacelocal = (*faceIndex)[iFace];
        int le = leftCellOfFace[iFacelocal];
        int re = rightCellOfFace[iFacelocal];

        u[re] -= hu[le] * dpcdx[le];
        v[re] -= hv[le] * dpcdy[le];
        w[re] -= hw[le] * dpcdz[le];

        RDouble area = faceArea[iFacelocal];
        RDouble fanx = xfn[iFacelocal] * area;
        RDouble fany = yfn[iFacelocal] * area;
        RDouble fanz = zfn[iFacelocal] * area;

        fluxbValue = rho[le] * (u[re] * fanx + v[re] * fany + w[re] * fanz);
        fluxb[iFacelocal] = fluxbValue;

        un = u[re] * fanx + v[re] * fany + w[re] * fanz;
        dunb[iFacelocal] = fluxbValue / rho[le] - un;
    }
}

DEFINE_PBSOLVER_REGISTER_FUNCTION(PBBoundaryFlow, PBBoundaryFlowPressureOutlet, FLOW_PRESSURE_OUTLET);
}