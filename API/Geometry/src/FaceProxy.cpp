#include "FaceProxy.h"
#include "Geo_Interface.h"
#include "Geo_Interpoint.h"
#include "Geo_FaceMetrics_Unstruct.h"
#include "Geo_CellMetrics_Unstruct.h"
#include "Geo_DynamicGridMetrics_Unstruct.h"
#include "Geo_LSQWeight_Unstruct.h"
#include "Geo_UnstructBC.h"
#include "Pointer.h"
#include "GeomProxy.h"
#include "GeneralFieldProxy.h"
#include "Geo_UnstructGrid.h"
#include "Geo_FaceTopo_Unstruct.h"
#include "Math_Limiter.h"

namespace PHSPACE
{
FaceProxy::FaceProxy()
{
    ql   = 0;
    qr   = 0;
    flux = 0;

    tl   = 0;
    tr   = 0;

    gamal = 0;
    gamar = 0;

    deltl = 0;
    deltr = 0;

    lrtl  = 0;
    lrtr  = 0;

    ns_facevar   = 0;
    turb_facevar = 0;
    geom_proxy   = 0;

    next = 0;
    nlen = 0;
    nsize= 0;
    neqn = 0;
}

FaceProxy::~FaceProxy()
{
    if (next) delete next;

    DelPointer2(ql);
    DelPointer2(qr);
    DelPointer2(flux);

    DelPointer2(tl);
    DelPointer2(tr);

    delete [] gamal;    gamal = NULL;
    delete [] gamar;    gamar = NULL;

    delete [] deltl;    deltl = NULL;
    delete [] deltr;    deltr = NULL;

    delete [] lrtl;    lrtl = NULL;
    delete [] lrtr;    lrtr = NULL;

    delete ns_facevar;
    delete turb_facevar;
    delete geom_proxy;
}

void FaceProxy::Create(int nsize, int neqn)
{
    this->neqn  = neqn;
    this->nsize = nsize;

    ql   = NewPointer2<RDouble>(neqn, nsize);
    qr   = NewPointer2<RDouble>(neqn, nsize);
    flux = NewPointer2<RDouble>(neqn, nsize);

    gamal = new RDouble[nsize];
    gamar = new RDouble[nsize];

    deltl = new RDouble[nsize];
    deltr = new RDouble[nsize];

    lrtl  = new RDouble[nsize];
    lrtr  = new RDouble[nsize];

    tl = NewPointer2<RDouble>(3, nsize);
    tr = NewPointer2<RDouble>(3, nsize);
}

FaceProxy::value_type FaceProxy::GetQL()
{
    return ql;
}

FaceProxy::value_type FaceProxy::GetQR()
{
    return qr;
}

FaceProxy::value_type FaceProxy::GetFlux()
{
    return flux;
}

RDouble * FaceProxy::GetGamaL()
{
    return gamal;
}

RDouble * FaceProxy::GetGamaR()
{
    return gamar;
}

RDouble ** FaceProxy::GetLeftTemperature()
{
    return tl;
}

RDouble ** FaceProxy::GetRightTemperature()
{
    return tr;
}

RDouble * FaceProxy::GetPressureCoefficientL()
{
    return lrtl;
}

RDouble * FaceProxy::GetPressureCoefficientR()
{
    return lrtr;
}

RDouble * FaceProxy::GetWeightL()
{
    return deltl;
}

RDouble * FaceProxy::GetWeightR()
{
    return deltr;
}

NSFaceValue * FaceProxy::GetNSFaceValue()
{
    return ns_facevar;
}

void FaceProxy::SetNSFaceValue(NSFaceValue *facevar)
{
    if (ns_facevar) delete ns_facevar;
    ns_facevar = facevar;
}

TurbFaceValue * FaceProxy::GetTurbFaceValue()
{
    return turb_facevar;
}

void FaceProxy::SetTurbFaceValue(TurbFaceValue *facevar)
{
    if (turb_facevar) delete turb_facevar;
    turb_facevar = facevar;
}

GeomProxy * FaceProxy::GetGeomProxy()
{
    return geom_proxy;
}

void FaceProxy::SetGeomProxy(GeomProxy *geom_proxy)
{
    if (this->geom_proxy) delete this->geom_proxy;
    this->geom_proxy = geom_proxy;
}

NSFaceValue::NSFaceValue(int neqn, int ns, int nlen)
{
    rho_ds = new GeneralFieldProxy();
    hint_s = new GeneralFieldProxy();
    if (!(ns == 0 || nlen == 0))
    {
        rho_ds->Create(ns,nlen);
        hint_s->Create(ns,nlen);
    }

    kcp = new RDouble[nlen];
    mul = new RDouble[nlen];
    mut = new RDouble[nlen];

    prim = new GeneralFieldProxy();
    prim->Create(neqn,nlen);

    t = new GeneralFieldProxy();
    t->Create(1,nlen);
}

NSFaceValue::~NSFaceValue()
{
    delete hint_s;
    delete rho_ds;
    delete prim;
    delete t;

    delete [] kcp;    kcp = NULL;
    delete [] mul;    mul = NULL;
    delete [] mut;    mut = NULL;
}

TurbFaceValue::TurbFaceValue(int neqn, int nlen)
{
    mul = new RDouble[nlen];
    mut = new RDouble[nlen];

    mlt  = NewPointer2<RDouble>(neqn, nlen);
    prim = NewPointer2<RDouble>(neqn, nlen);
}

TurbFaceValue::~TurbFaceValue()
{
    delete [] mul;    mul = NULL;
    delete [] mut;    mut = NULL;

    DelPointer2(mlt);
    DelPointer2(prim);
}

void ReConstructQlQr_STR(Grid *grid_in, FaceProxy *face_proxy, int neqn, int nst, int ned)
{
    UnstructGrid *grid = UnstructGridCast(grid_in);

    int     *left_cell_of_face = grid->GetLeftCellOfFace();
    int     *right_cell_of_face = grid->GetRightCellOfFace();
    RDouble *xfc = grid->GetFaceCenterX();
    RDouble *yfc = grid->GetFaceCenterY();
    RDouble *zfc = grid->GetFaceCenterZ();
    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    RDouble **ql = face_proxy->GetQL();
    RDouble **qr = face_proxy->GetQR();

    int uns_limiter;
    GlobalDataBase::GetData("uns_limiter", &uns_limiter, PHINT, 1);

    RDouble *qltry = new RDouble[neqn];
    RDouble *qrtry = new RDouble[neqn];

    RDouble MUSCLCoefXk;
    GlobalDataBase::GetData("MUSCLCoefXk", &MUSCLCoefXk, PHDOUBLE, 1);
    RDouble MUSCLCoefXb;
    GlobalDataBase::GetData("MUSCLCoefXb", &MUSCLCoefXb, PHDOUBLE, 1);

    RDouble dx, dy, dz;
    RDouble c1, c2;

    c1 = fourth * (1.0 - MUSCLCoefXk);
    c2 = fourth * (1.0 + MUSCLCoefXk);

    RDouble (*LIMIT_FUN)(const RDouble &x, const RDouble &y);

    if (uns_limiter == ILMT_MINMOD)
    {
        LIMIT_FUN = MINMOD;
    }
    else if (uns_limiter == ILMT_VAN_ALBADA)
    {
        LIMIT_FUN = VANALBADA;
    }
    else if (uns_limiter == ILMT_SMOOTH)
    {
        LIMIT_FUN = SMOOTH;
    }
    else
    {
        LIMIT_FUN = VANALBADA;
    }

    RDouble **dqdx = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarX"));
    RDouble **dqdy = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarY"));
    RDouble **dqdz = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarZ"));

    int le,re,j;
    for (int iface = nst; iface < ned; ++ iface)
    {
        j  = iface - nst;
        le = left_cell_of_face [ iface ];
        re = right_cell_of_face[ iface ];

        for (int m = 0; m < neqn; ++ m)
        {
            qltry[m] = ql[m][j];
            qrtry[m] = qr[m][j];
        }

        dx = xfc[iface] - xcc[le];
        dy = yfc[iface] - ycc[le];
        dz = zfc[iface] - zcc[le];

        for (int m = 0; m < neqn; ++ m)
        {
            RDouble lx = LIMIT_FUN(dqdx[m][le],dqdx[m][re]);
            RDouble ly = LIMIT_FUN(dqdy[m][le],dqdy[m][re]);
            RDouble lz = LIMIT_FUN(dqdz[m][le],dqdz[m][re]);
            RDouble limdq = lx * dx + ly * dy + lz * dz;

            //limdq = MINMOD(limdq,0.5*(qr[m][j]-ql[m][j]));

            qltry[m] += limdq;
        }

        dx = xfc[iface] - xcc[re];
        dy = yfc[iface] - ycc[re];
        dz = zfc[iface] - zcc[re];

        for (int m = 0; m < neqn; ++ m)
        {
            RDouble lx = LIMIT_FUN(dqdx[m][le],dqdx[m][re]);
            RDouble ly = LIMIT_FUN(dqdy[m][le],dqdy[m][re]);
            RDouble lz = LIMIT_FUN(dqdz[m][le],dqdz[m][re]);
            RDouble limdq = lx * dx + ly * dy + lz * dz;

            //limdq = MINMOD(limdq,0.5*(ql[m][j]-qr[m][j]));

            qrtry[m] += limdq;
        }

        if (PositiveCheckForDensityAndPressure(qltry))
        {
            for (int m = 0; m < neqn; ++ m)
            {
                ql[m][j] = qltry[m];
            }
        }

        if (PositiveCheckForDensityAndPressure(qrtry))
        {
            for (int m = 0; m < neqn; ++ m)
            {
                qr[m][j] = qrtry[m];
            }
        }
    }

    delete [] qltry;    qltry = NULL;
    delete [] qrtry;    qrtry = NULL;
}

}

