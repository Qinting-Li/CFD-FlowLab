#include "ForceProcessStruct.h"
#include "Geo_StructBC.h"
#include "Glb_Dimension.h"
#include "GradientOperation.h"
#include "Post_ForceMoment.h"
#include <cmath>
#include "Geo_StructGrid.h"

namespace PHSPACE
{

ForceProcessStruct::ForceProcessStruct(PostProcessData *Data)
{
    this->ForceData = Data;
    CollectData = new DataContainer();
}

ForceProcessStruct::~ForceProcessStruct()
{
    delete CollectData;
}

void ForceProcessStruct::Run()
{
    InitValueField();
    CompValueField();
    InitForce();
    CollectionForce();
    CompForce();
}

void ForceProcessStruct::InitValueField()
{
    using namespace PHMPI;

    int nLocalZones = GetNumberofLocalZones();
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);
        StructGrid *grid = StructGridCast(GetGrid(zoneID));
        
        int ni = grid->GetNI();
        int nj = grid->GetNJ();
        int nk = grid->GetNK();
        Range I(-1, ni+1);
        Range J(-1, nj+1);
        Range K(-1, nk+1);
        if (nk == 1)
        {
            K.setRange(1, 1);
        }
        int nTemperatureModel = GlobalDataBase::GetIntParaFromDB("ntmodel");

        Range T(0, nTemperatureModel - 1);
        RDouble4D *t = new RDouble4D(I, J, K, T, fortranArray);
        grid->UpdateDataPtr("t", t);

        RDouble3D *visl = new RDouble3D(I, J, K, fortranArray);
        RDouble3D *vist = new RDouble3D(I, J, K, fortranArray);
        grid->UpdateDataPtr("visl", visl);
        grid->UpdateDataPtr("vist", vist);
    }
}

void ForceProcessStruct::InitForce()
{
    RDouble attack, sideslip;
    RDouble attackd = GlobalDataBase::GetDoubleParaFromDB("attackd");
    RDouble angleSlide = GlobalDataBase::GetDoubleParaFromDB("angleSlide");

    attack   = attackd * PI / 180.0;
    sideslip = angleSlide * PI / 180.0;

    RDouble forceReferenceArea   = GlobalDataBase::GetDoubleParaFromDB("forceReferenceArea");
    RDouble forceReferenceLength = GlobalDataBase::GetDoubleParaFromDB("forceReferenceLength");
    RDouble forceReferenceLengthSpanWise = GlobalDataBase::GetDoubleParaFromDB("forceReferenceLengthSpanWise");
    RDouble TorqueRefX = GlobalDataBase::GetDoubleParaFromDB("TorqueRefX");
    RDouble TorqueRefY = GlobalDataBase::GetDoubleParaFromDB("TorqueRefY");
    RDouble TorqueRefZ = GlobalDataBase::GetDoubleParaFromDB("TorqueRefZ");

    globalAerodynamicForce.Initialize(attack, sideslip, forceReferenceArea, forceReferenceLength, forceReferenceLengthSpanWise);
    globalAerodynamicForce.SetTorqueRef(TorqueRefX, TorqueRefY, TorqueRefZ);
}

void ForceProcessStruct::CollectionForce()
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int zone_proc = GetZoneProcessorID(iZone);
        int myid = GetCurrentProcessorID();

        if (myid == zone_proc)
        {
            AirForceCoef(iZone);
        }
    }
    globalAerodynamicForce.CollectOfAllProcessors();
}

void ForceProcessStruct::AirForceCoef(int iZone)
{
    const int ALL_PART = -1;
    AirForceCoefParts(iZone, ALL_PART);

    uint_t nTBC = GlobalBoundaryCondition::GetNumberofBC();
    for (int iPart = 0; iPart < nTBC; ++ iPart)
    {
        if (!GlobalBoundaryCondition::IsSolidWallBC(iPart)) continue;

        AirForceCoefParts(iZone, iPart);
    }
}

void ForceProcessStruct::AirForceCoefParts(int iZone, int partID)
{
    StructGrid *grid = StructGridCast(GetGrid(iZone));

    RDouble4D &q    = * reinterpret_cast<RDouble4D *>(grid->GetDataPtr("q"));
    RDouble3D &visl = * reinterpret_cast<RDouble3D *>(grid->GetDataPtr("visl"));

    RDouble refReNumber = GlobalDataBase::GetDoubleParaFromDB("refReNumber");
    RDouble oreynolds = 1.0 / refReNumber;

    int viscousType = GlobalDataBase::GetIntParaFromDB("viscousType");

    const int ALL_PART = -1;
    SimpleBC *boundaryCondition = 0;
    string partBCName = "";
    RDouble TorqueRefX, TorqueRefY, TorqueRefZ;

    //! used to compute hinge moment:
    int dumpHingeMoment = 0;
    RDouble localCoordAxis0[3];
    RDouble localCoordAxis1[3];
    RDouble hingeMoment = 0.0;

    if (partID == ALL_PART)
    {
        TorqueRefX = GlobalDataBase::GetDoubleParaFromDB("TorqueRefX");
        TorqueRefY = GlobalDataBase::GetDoubleParaFromDB("TorqueRefY");
        TorqueRefZ = GlobalDataBase::GetDoubleParaFromDB("TorqueRefZ");
    }
    else
    {
        boundaryCondition = GlobalBoundaryCondition::GetBC(partID);
        boundaryCondition->GetData("TorqueRefX", &TorqueRefX, PHDOUBLE, 1);
        boundaryCondition->GetData("TorqueRefY", &TorqueRefY, PHDOUBLE, 1);
        boundaryCondition->GetData("TorqueRefZ", &TorqueRefZ, PHDOUBLE, 1);
        partBCName = boundaryCondition->GetBCName();

        boundaryCondition->GetData("dumpHingeMoment", &dumpHingeMoment, PHINT, 1);
        boundaryCondition->GetData("localCoordAxis0", localCoordAxis0, PHDOUBLE, 3);
        boundaryCondition->GetData("localCoordAxis1", localCoordAxis1, PHDOUBLE, 3);
    }

    using namespace IDX;
    RDouble *prim_inf = reinterpret_cast<RDouble *>(GlobalDataBase::GetDataPtr("prim_inf"));
    RDouble poo = prim_inf[IP];

    RDouble4D &xfn  = * (grid->GetFaceNormalX());
    RDouble4D &yfn  = * (grid->GetFaceNormalY());
    RDouble4D &zfn  = * (grid->GetFaceNormalZ());
    RDouble4D &xfv  = * (grid->GetFaceVectorX());
    RDouble4D &yfv  = * (grid->GetFaceVectorY());
    RDouble4D &zfv  = * (grid->GetFaceVectorZ());
    RDouble4D &area = * (grid->GetFaceArea());
    RDouble3D &vol  = * (grid->GetCellVolume());

    RDouble cfx, cfy, cfz, cpx, cpy, cpz, cmx, cmy, cmz;
    cfx = 0.0;
    cfy = 0.0;
    cfz = 0.0;
    cpx = 0.0;
    cpy = 0.0;
    cpz = 0.0;
    cmx = 0.0;
    cmy = 0.0;
    cmz = 0.0;

    int ndim = GetDim();
    using namespace IDX;

    int ix = 0;

    CompositeStructBC *composite_bcregion = grid->GetCompositeBCRegion();
    int nBCRegion = composite_bcregion->GetnBCRegion();
    
    for (int ibcregion = 0; ibcregion < nBCRegion; ++ ibcregion)
    {
        StructBC *bcregion = composite_bcregion->GetBCRegion(ibcregion);
        int bctype = bcregion->GetBCType();
        
        int ist, ied, jst, jed, kst, ked;
        bcregion->GetIJKRegion(ist, ied, jst, jed, kst, ked);
        
        bool isNeedComputeAirForce = (bctype == PHENGLEI::SOLID_SURFACE);
        if (isNeedComputeAirForce && partID != ALL_PART)
        {
            const string &bcName = bcregion->GetBCName();
            isNeedComputeAirForce = isNeedComputeAirForce && (bcName == partBCName);
        }

        if (!isNeedComputeAirForce)
        {
            ix += (ied - ist + 1) * (jed - jst + 1) * (ked - kst + 1);
            continue;
        }

        int *s_lr3d = bcregion->GetFaceDirectionIndex();

        int nsurf = bcregion->GetFaceDirection() + 1;
        int lr    = bcregion->GetFaceLeftOrRightIndex();

        int id, jd, kd;
        GetBCFaceIDX(s_lr3d, id, jd, kd);

        int in, jn, kn, i1, j1, k1, il, jl, kl;
        RDouble pw, p1, p2, cp;
        RDouble dfx, dfy, dfz, dpx, dpy, dpz, dmx, dmy, dmz, fvsx, fvsy, fvsz;
        RDouble nxs, nys, nzs, ns, nx, ny, nz;
        RDouble dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz;
        RDouble xc, yc, zc;

        RDouble vsl, vst;
        RDouble vis, divv2p3;
        RDouble txx, tyy, tzz;
        RDouble txy, txz, tyz;

        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    in = i + id;
                    jn = j + jd;
                    kn = k + kd;

                    nxs = xfn (in, jn, kn, nsurf);
                    nys = yfn (in, jn, kn, nsurf);
                    nzs = zfn (in, jn, kn, nsurf);
                    ns  = area(in, jn, kn, nsurf);

                    nx = lr * nxs * ns;
                    ny = lr * nys * ns;
                    nz = lr * nzs * ns;

                    i1 = i - s_lr3d[0];
                    j1 = j - s_lr3d[1];
                    k1 = k - s_lr3d[2];

                    il = i + s_lr3d[0];
                    jl = j + s_lr3d[1];
                    kl = k + s_lr3d[2];

                    p2 = q(i1, j1, k1, IP);
                    p1 = q(i, j, k, IP);

                    //! pressure drag.
                    pw = p1;
                    cp = two * (pw - poo);

                    dpx = nx * cp;
                    dpy = ny * cp;
                    dpz = nz * cp;

                    dfx = dpx;
                    dfy = dpy;
                    dfz = dpz;

                    if (viscousType > INVISCID)
                    {
                        Geometry3D *geometry3D = new Geometry3D;

                        if (nsurf == 1)
                        {
                            //Get ith-facevector on the wall surface;
                            Get_GEO_I(i, j, k, in, in, ndim, geometry3D, xfv, yfv, zfv, vol);
                            //! Get derivative to the wall in i-direction of u, v, w.
                            DXDYDZ_I(q, geometry3D, i, j, k, ndim, IU, dudx, dudy, dudz, lr);
                            DXDYDZ_I(q, geometry3D, i, j, k, ndim, IV, dvdx, dvdy, dvdz, lr);
                            DXDYDZ_I(q, geometry3D, i, j, k, ndim, IW, dwdx, dwdy, dwdz, lr);
                            vsl = half * (visl(i, j, k) + visl(i+lr, j, k));
                        }
                        else if (nsurf == 2)
                        {
                            //! Get jth-facevector on the wall surface.
                            Get_GEO_J(i, j, k, jn, jn, ndim, geometry3D, xfv, yfv, zfv, vol);
                            //! Get derivative to the wall in j-direction of u, v, w.
                            DXDYDZ_J(q, geometry3D, i, j, k, ndim, IU, dudx, dudy, dudz, lr);
                            DXDYDZ_J(q, geometry3D, i, j, k, ndim, IV, dvdx, dvdy, dvdz, lr);
                            DXDYDZ_J(q, geometry3D, i, j, k, ndim, IW, dwdx, dwdy, dwdz, lr);
                            vsl = half * (visl(i, j, k) + visl(i, j+lr, k));
                        }
                        else
                        {
                            //! Get kth-facevector on the wall surface.
                            Get_GEO_K(i, j, k, kn, kn, ndim, geometry3D, xfv, yfv, zfv, vol);
                            //! Get derivative to the wall in k-direction of u, v, w.
                            DXDYDZ_K(q, geometry3D, i, j, k, ndim, IU, dudx, dudy, dudz, lr);
                            DXDYDZ_K(q, geometry3D, i, j, k, ndim, IV, dvdx, dvdy, dvdz, lr);
                            DXDYDZ_K(q, geometry3D, i, j, k, ndim, IW, dwdx, dwdy, dwdz, lr);
                            vsl = half * (visl(i, j, k) + visl(i, j, k+lr));
                        }

                        vsl = half * (visl(il, jl, kl) + visl(i, j, k));
                        vst = 0.0;
                        vis = vsl + vst;

                        divv2p3 = two3rd * (dudx + dvdy + dwdz);

                        //! stress components.
                        txx = vis * (two * dudx - divv2p3);
                        tyy = vis * (two * dvdy - divv2p3);
                        tzz = vis * (two * dwdz - divv2p3);
                        txy = vis * (dudy + dvdx);
                        txz = vis * (dudz + dwdx);
                        tyz = vis * (dvdz + dwdy);

                        fvsx = - two * (nx * txx + ny * txy + nz * txz) * oreynolds;
                        fvsy = - two * (nx * txy + ny * tyy + nz * tyz) * oreynolds;
                        fvsz = - two * (nx * txz + ny * tyz + nz * tzz) * oreynolds;

                        dfx += fvsx;
                        dfy += fvsy;
                        dfz += fvsz;

                        delete geometry3D;
                    }

                    grid->FaceCoor(in, jn, kn, nsurf, xc, yc, zc);

                    dmx = (yc - TorqueRefY) * dfz - (zc - TorqueRefZ) * dfy;
                    dmy = (zc - TorqueRefZ) * dfx - (xc - TorqueRefX) * dfz;
                    dmz = (xc - TorqueRefX) * dfy - (yc - TorqueRefY) * dfx;

                    cpx += dpx;
                    cpy += dpy;
                    cpz += dpz;

                    cfx += dfx;
                    cfy += dfy;
                    cfz += dfz;

                    cmx += dmx;
                    cmy += dmy;
                    cmz += dmz;

                    if (dumpHingeMoment)
                    {
                        Post_ForceMoment forceMoment;
                        RDouble point[3] = { xc, yc, zc };
                        RDouble faceForce[3] = { dfx, dfy, dfz };
                        hingeMoment += forceMoment.ComputeMoment(point, localCoordAxis0, localCoordAxis1, faceForce);
                    }

                    ++ ix;
                }
            }
        }
    }

    this->CollectData->write(&cfx, sizeof(RDouble));
    this->CollectData->write(&cfy, sizeof(RDouble));
    this->CollectData->write(&cfz, sizeof(RDouble));
    this->CollectData->write(&cpx, sizeof(RDouble));
    this->CollectData->write(&cpy, sizeof(RDouble));
    this->CollectData->write(&cpz, sizeof(RDouble));
    this->CollectData->write(&cmx, sizeof(RDouble));
    this->CollectData->write(&cmy, sizeof(RDouble));
    this->CollectData->write(&cmz, sizeof(RDouble));
    this->CollectData->write(&hingeMoment, sizeof(RDouble));
}

void ForceProcessStruct::CompForce()
{
    DataContainer *data = new DataContainer();

    globalAerodynamicForce.ComputeCoefficient();
    globalAerodynamicForce.ComputeCoefficient();

    DecodeData(data);

    delete data;
}

void ForceProcessStruct::DecodeData(DataContainer *data)
{
    data->moveto_begin();

    int numberofSolidWallPart = 0;
    PHRead(data, numberofSolidWallPart);
    ForceData->SetNWallPart(numberofSolidWallPart);

    int numberofTotalPart = 0;
    PHRead(data, numberofTotalPart);
    ForceData->SetNPart(numberofTotalPart);

    int *partBCID        = new int     [numberofTotalPart];
    RDouble *cfx         = new RDouble [numberofTotalPart];
    RDouble *cfy         = new RDouble [numberofTotalPart];
    RDouble *cfz         = new RDouble [numberofTotalPart];
    RDouble *cpx         = new RDouble [numberofTotalPart];
    RDouble *cpy         = new RDouble [numberofTotalPart];
    RDouble *cpz         = new RDouble [numberofTotalPart];
    RDouble *cmx         = new RDouble [numberofTotalPart];
    RDouble *cmy         = new RDouble [numberofTotalPart];
    RDouble *cmz         = new RDouble [numberofTotalPart];
    RDouble *cl_tot      = new RDouble [numberofTotalPart];
    RDouble *cd_tot      = new RDouble [numberofTotalPart];
    RDouble *cd_pr       = new RDouble [numberofTotalPart];
    RDouble *cd_sf       = new RDouble [numberofTotalPart];
    RDouble *cd_cl2pa    = new RDouble [numberofTotalPart];
    RDouble *xcp         = new RDouble [numberofTotalPart];
    RDouble *side        = new RDouble [numberofTotalPart];
    RDouble *hingeMoment = new RDouble [numberofTotalPart];

    PHRead(data, partBCID,    numberofTotalPart);
    PHRead(data, cfx,         numberofTotalPart);
    PHRead(data, cfy,         numberofTotalPart);
    PHRead(data, cfz,         numberofTotalPart);
    PHRead(data, cpx,         numberofTotalPart);
    PHRead(data, cpy,         numberofTotalPart);
    PHRead(data, cpz,         numberofTotalPart);
    PHRead(data, cmx,         numberofTotalPart);
    PHRead(data, cmy,         numberofTotalPart);
    PHRead(data, cmz,         numberofTotalPart);
    PHRead(data, cl_tot,      numberofTotalPart);
    PHRead(data, cd_tot,      numberofTotalPart);
    PHRead(data, cd_pr,       numberofTotalPart);
    PHRead(data, cd_sf,       numberofTotalPart);
    PHRead(data, cd_cl2pa,    numberofTotalPart);
    PHRead(data, xcp,         numberofTotalPart);
    PHRead(data, side,        numberofTotalPart);
    PHRead(data, hingeMoment, numberofTotalPart);

    ForceData->SetPartBCID(partBCID);

    ForceData->SetCfx(cfx);
    ForceData->SetCfy(cfy);
    ForceData->SetCfz(cfz);
    ForceData->SetCpx(cpx);
    ForceData->SetCpy(cpy);
    ForceData->SetCpz(cpz);
    ForceData->SetCmx(cmx);
    ForceData->SetCmy(cmy);
    ForceData->SetCmz(cmz);
    ForceData->SetClTot(cl_tot);
    ForceData->SetCdTot(cd_tot);
    ForceData->SetCdPr(cd_pr);
    ForceData->SetCdSf(cd_sf);
    ForceData->SetCdCl2Pa(cd_cl2pa);
    ForceData->SetXcp(xcp);
    ForceData->SetSide(side);
    ForceData->SetHingeMoment(hingeMoment);
}

void ForceProcessStruct::CompValueField()
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();

    #ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int izone = 0; izone < nLocalZones; ++ izone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(izone);

        Grid *grid = GetGrid(zoneID,0);

        ComputeGamaAndTemperature(grid);
        ComputeViscousLaminarCoefficient(grid);
    }
}

void ForceProcessStruct::ComputeGamaAndTemperature(Grid *grid_in)
{
    StructGrid *grid = StructGridCast(grid_in);

    RDouble4D &q = * reinterpret_cast<RDouble4D *>(grid->GetDataPtr("q"));
    RDouble4D &t = * reinterpret_cast<RDouble4D *>(grid->GetDataPtr("t"));

    int nm = 5;
    GlobalDataBase::GetData("nm", &nm, PHINT, 1);

    int nchem = 0;
    GlobalDataBase::GetData("nchem", &nchem, PHINT, 1);

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 2);

    using namespace IDX;

    RDouble omav = one;
    RDouble coefficientOfStateEquation = GlobalDataBase::GetDoubleParaFromDB("coefficientOfStateEquation");

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble &rm = q(i, j, k, IR);
                RDouble &pm = q(i, j, k, IP);
                t(i, j, k, ITT)  = pm / (coefficientOfStateEquation * rm * omav);
            }
        }
    }
}

void ForceProcessStruct::ComputeViscousLaminarCoefficient(Grid *grid_in)
{
    StructGrid *grid = StructGridCast(grid_in);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    RDouble3D &visl = * reinterpret_cast<RDouble3D *>(grid->GetDataPtr("visl"));
    RDouble4D &t    = * reinterpret_cast<RDouble4D *>(grid->GetDataPtr("t"));

    RDouble tsuth;
    GlobalDataBase::GetData("tsuth", &tsuth, PHDOUBLE, 1);

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    using namespace IDX;

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble tm = t(i, j, k, ITT);
                visl(i, j, k) = tm * sqrt(tm) * (1.0 + tsuth) / (tm + tsuth);
            }
        }
    }

    RDouble wallTemperature = GlobalDataBase::GetDoubleParaFromDB("wallTemperature");

    if (wallTemperature <= 0.0)
    {
        GhostCell3D(visl, ni, nj, nk);
    }
    else
    {
        int idx_bc[2];
        idx_bc[0] = 0;

        idx_bc[1] = ni;
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int m = 0; m <= 1; ++ m)
                {
                    int i = idx_bc[m];
                    RDouble tm = t(i, j, k, ITT);
                    visl(i, j, k) = tm * sqrt(tm) * (1.0 + tsuth) / (tm + tsuth);
                }
            }
        }

        idx_bc[1] = nj;
        for (int k = kst; k <= ked; ++ k)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                for (int m = 0; m <= 1; ++ m)
                {
                    int j = idx_bc[m];
                    RDouble tm = t(i, j, k, ITT);
                    visl(i, j, k) = tm * sqrt(tm) * (1.0 + tsuth) / (tm + tsuth);
                }
            }
        }

        if (GetDim() == THREE_D)
        {
            idx_bc[1] = nk;
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    for (int m = 0; m <= 1; ++ m)
                    {
                        int k = idx_bc[m];
                        RDouble tm = t(i, j, k, ITT);
                        visl(i, j, k) = tm * sqrt(tm) * (1.0 + tsuth) / (tm + tsuth);
                    }
                }
            }
        }
    }
}

}