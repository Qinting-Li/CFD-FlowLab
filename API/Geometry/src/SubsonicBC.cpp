#include "SubsonicBC.h"
#include "GridType.h"
#include "Geo_StructBC.h"
#include "Geo_UnstructBC.h"
#include "Constants.h"
#include "Geo_StructGrid.h"
#include "Geo_MultiGridInfo_Struct.h"
#include "Geo_UnstructGrid.h"
#include "TK_Exit.h"
#include "Geo_Interface.h"
#include "Geo_Interpoint.h"
#include "Geo_FaceTopo_Unstruct.h"
#include "Geo_NodeTopo_Struct.h"
#include "Geo_FaceMetrics_Struct.h"
#include "Geo_CellMetrics_Struct.h"
#include "Geo_DynamicGridMetrics_Struct.h"
#include "Geo_OversetGridTopo_Struct.h"
#include "Geo_FaceMetrics_Unstruct.h"
#include "Geo_CellMetrics_Unstruct.h"
#include "Geo_DynamicGridMetrics_Unstruct.h"
#include "Geo_LSQWeight_Unstruct.h"

using namespace std;
using namespace PHSPACE;
using namespace PHMPI;

namespace PHSPACE
{
    SubnicBCInfo outLetBCs;
    SubnicBCInfo inLetBCs;

void SubnicBCInfo::CheckSpecificBCExist(int specificBCID)
{

    int hasBCp,hasBC;

    hasBC  = 0;
    hasBCp = 0;

    int nZones = GetNumberofGlobalZones();

    if(  GetCurrentProcessorID() == server )
    {
        for ( int izone = 0; izone < nZones; ++ izone )
        {
            Grid *gridIn = GetGrid(izone,0);

            if ( !gridIn ) continue;

            int sys_gridtype = STRUCTGRID;
            sys_gridtype = GlobalDataBase::GetIntParaFromDB( "sys_gridtype" );

            if( sys_gridtype == PHSPACE::STRUCTGRID )
            {
                StructGrid *grid = StructGridCast(gridIn);

                CompositeStructBC * composite_bcregion = grid->GetCompositeBCRegion();
                int nBCRegion = composite_bcregion->GetnBCRegion();

                for ( int ibcregion = 0; ibcregion < nBCRegion; ++ ibcregion )
                {
                    StructBC * bcregion = composite_bcregion->GetBCRegion(ibcregion);
                    if ( bcregion->GetBCType() == specificBCID )
                    {
                        hasBCp = 1;
                        break;
                    }
                }
            }
            else if( sys_gridtype == PHSPACE::UNSTRUCTGRID )
            {
                UnstructGrid *grid = UnstructGridCast( gridIn );
                CompositeUnstructBC ** bcr = grid->GetBCRecord();
                for ( int iFace = 0; iFace < grid->GetNBoundFace(); ++ iFace )
                {
                    int type = bcr[ iFace ]->GetKey();

                    if ( type == specificBCID )
                    {
                        hasBCp = 1;
                        break;
                    }
                }
            }
        }
    }
    
#ifdef PH_PARALLEL
    MPI_Reduce( &hasBCp, &hasBC, 1, MPI_INT, MPI_MAX, server, MPI_COMM_WORLD );
    MPI_Bcast( &hasBC,1,MPI_INT,server,MPI_COMM_WORLD);
#endif
    if ( hasBC == 1 ) 
    {
        SetSimilarBCExistLabel();
    }
}

void SubnicBCInfo::Nondimensional()
{
    for (int i = 0; i < GetZoneNumberOfSameBCType(); ++ i) 
    {
        RDouble p0 = BCs[i].GetTotalPressure();
        RDouble t0 = BCs[i].GetTotalTemperature();

        RDouble refGama = 1.4;
        RDouble refMachNumber = 1.0;
        refGama = GlobalDataBase::GetDoubleParaFromDB("refGama");
        refMachNumber = GlobalDataBase::GetDoubleParaFromDB("refMachNumber");

        RDouble staticPressure    = GlobalDataBase::GetDoubleParaFromDB("refDimensionalPressure");
        RDouble staticTemperature = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");

        p0 /= (staticPressure * refGama * refMachNumber * refMachNumber);
        t0 /= staticTemperature;

        BCs[i].SetTotalPressure(p0);
        BCs[i].SetTotalTemperature(t0);
    }
}

void SubnicBCInfo::ReSetReferenceData()
{
    if ( hasSimilarBC == 1 ) return;

    RDouble p0Average = 0.0, t0Average = 0.0, mrAverage = 0.0;
    vector<RDouble> dirAverage(3,0.0); 
    int nZones = GetZoneNumberOfSameBCType();

    for ( int i = 0; i < nZones; ++ i ) 
    {
        p0Average += BCs[i].GetTotalPressure();
        t0Average += BCs[i].GetTotalTemperature();
        mrAverage += BCs[i].GetMassRate();

        for ( int j = 0; j < 3; ++ j )
        {
            dirAverage[j] += BCs[i].GetFlowDirection()[j];
        }
    }

    p0Average /= nZones;
    t0Average /= nZones;
    mrAverage /= nZones;

    for ( int j = 0; j < 3; ++ j )
    {
        dirAverage[j] /= nZones;
    }

    RDouble sum = sqrt( dirAverage[0]*dirAverage[0] + dirAverage[1]*dirAverage[1] + dirAverage[2]*dirAverage[2] );

    for ( int j = 0; j < 3; ++ j )
    {
        dirAverage[j] /= sum;
    }

    RDouble refGama = GlobalDataBase::GetDoubleParaFromDB("refGama");
    RDouble lref = GlobalDataBase::GetDoubleParaFromDB( "reynoldsReferenceLengthDimensional" );

    RDouble moo    = 1.0;
    RDouble mo2    = moo * moo;

    RDouble pref   = p0Average * pow((1.0 + (refGama - 1.0) / 2.0 * mo2), (- refGama / (refGama -1.0)));
    RDouble tref   = t0Average / (1.0 + (refGama - 1.0) / 2.0 * mo2);
    RDouble gasConstant = 287.053;
    RDouble rref   = pref / ( tref * gasConstant );
    RDouble ccoo   = sqrt(refGama * pref / rref);
    RDouble vref   = ccoo * moo;
    RDouble ccon   = ( 1.0 + 110.4 / 288.15 ) / ( 1.0 + 110.4 / tref );
    RDouble visloo = sqrt( tref / 288.15 ) * ccon * 1.7894e-5;
    RDouble refReNumber = rref * vref * lref / visloo;

    GlobalDataBase::UpdateData("refReNumber", &refReNumber, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("refMachNumber", &moo, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("refDimensionalTemperature", &tref, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("refDimensionalPressure", &pref, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("refDimensionalDensity", &rref, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("reference_velocity_dimensional", &vref, PHDOUBLE, 1);
}

void SubnicBCInfo::ReadBCData()
{
    using namespace PHMPI;
    if ( GetCurrentProcessorID() != server ) return;

    fstream file;
    string keyfile = GetBCDataFileName();

    file.open(keyfile.c_str(), ios_base::in);

    if ( !file )
    {
        TK_Exit::FileOpenErrorExit(keyfile);
    }

    int nZones,NO,bcID,bcType;
    RDouble p0 = 0, t0 = 0, massRate = 0;
    RDouble *fd = new RDouble [3]();

    file >> nZones;

    nbc = nZones;

    SubnicBC BC;

    for ( int i = 0; i < nZones; ++ i )
    {
        file >> NO >> bcID >> bcType;
        file >> p0 >> t0 >> massRate;
        file >> fd[0] >> fd[1] >> fd[2];

        BC.SetBCID( bcID );
        BC.SetBCType( bcType );
        BC.SetTotalPressure( p0 );
        BC.SetTotalTemperature( t0 );
        BC.SetFlowDirection( fd );

        BCs.push_back(BC);
    }
    delete [] fd;    fd = NULL;
    file.close();
    file.clear();
}

void SubnicBCInfo::BcastBCData()
{
#ifndef PH_PARALLEL
    return;
#endif
    int nZones = 0, bcID = -999, bcType = -999;
    RDouble p0,t0,massRate;
    RDouble *fd[3] ={0};
    nZones = GetZoneNumberOfSameBCType();

    using namespace PHMPI;

    MPI_Bcast(&nZones,1,MPI_INT,server,MPI_COMM_WORLD);

    if ( GetCurrentProcessorID() != server )
    {
        SubnicBC BC;
        for ( int i = 0; i < nZones; ++ i )
        {
            BCs.push_back( BC );
        }
    }
    
    for ( int i = 0; i < nZones; ++ i )
    {
        if ( GetCurrentProcessorID() == server )
        {
            bcID     = BCs[i].GetBCID();
            bcType   = BCs[i].GetBCType();
            p0       = BCs[i].GetTotalPressure();
            t0       = BCs[i].GetTotalTemperature();
            massRate = BCs[i].GetMassRate();
            *fd       = BCs[i].GetFlowDirection();
        }
        MPI_Bcast(&bcID    ,1,MPI_INT   ,server,MPI_COMM_WORLD);
        MPI_Bcast(&bcType  ,1,MPI_INT   ,server,MPI_COMM_WORLD);
        MPI_Bcast(&p0      ,1,MPI_DOUBLE,server,MPI_COMM_WORLD);
        MPI_Bcast(&t0      ,1,MPI_DOUBLE,server,MPI_COMM_WORLD);
        MPI_Bcast(&massRate,1,MPI_DOUBLE,server,MPI_COMM_WORLD);
        MPI_Bcast(fd       ,3,MPI_DOUBLE,server,MPI_COMM_WORLD);

        if ( GetCurrentProcessorID() != server )
        {
            BCs[i].SetBCID(bcID);
            BCs[i].SetBCType(bcType);
            BCs[i].SetTotalPressure(p0);
            BCs[i].SetTotalTemperature(t0);
            BCs[i].SetMassRate(massRate);
            BCs[i].SetFlowDirection(*fd);
        }
        cout << BCs[i].GetBCID() << "  " << GetCurrentProcessorID() << endl;
    }
    
}

void GetSubsonicBCData()
{
    int nsubsonicInlet  = 1;
    int nsubsonicOutlet = 0;

    nsubsonicInlet  = GlobalDataBase::GetIntParaFromDB( "nsubsonicInlet" );
    nsubsonicOutlet = GlobalDataBase::GetIntParaFromDB( "nsubsonicOutlet" );
    
    if ( nsubsonicInlet > 0 ) 
    {
        string inLetFileName = GlobalDataBase::GetStrParaFromDB( "inLetFileName" );
        inLetBCs.SetBCDataFileName(inLetFileName);

        inLetBCs.ReadBCData();
        inLetBCs.BcastBCData();
        
        inLetBCs.CheckSpecificBCExist( PHENGLEI::FARFIELD );
        inLetBCs.CheckSpecificBCExist( PHENGLEI::INFLOW );
        inLetBCs.ReSetReferenceData();
        inLetBCs.Nondimensional();
        inLetBCs.CalBCIDRange();

    }
    if ( nsubsonicOutlet > 0 ) 
    {
        string outLetFileName = GlobalDataBase::GetStrParaFromDB( "outLetFileName" );
        outLetBCs.SetBCDataFileName(outLetFileName);

        outLetBCs.ReadBCData();
        outLetBCs.BcastBCData();
        outLetBCs.Nondimensional();
        outLetBCs.CalBCIDRange();
        outLetBCs.CheckSpecificBCExist( PHENGLEI::OUTFLOW );
    }
}

bool SetSubsonicBC(Grid *gridIn, StructBC *bcr, int bctype)
{
    bool hasSubsonicBC = 1;

    if ( inLetBCs.GetBCIDRangLower() <= bctype &&  inLetBCs.GetBCIDRangUpper() >= bctype )
    {
        for ( int nZone = 0; nZone < inLetBCs.GetZoneNumberOfSameBCType(); ++ nZone )
        {
            SubnicBC BCS;
            BCS = inLetBCs.GetBCInfo()[nZone];

            if ( bctype == BCS.GetBCID() )
            {
                InFlowBC3D(gridIn, bcr, BCS);
                break;
            }
        }
    }
    else if( outLetBCs.GetBCIDRangLower() <= bctype && outLetBCs.GetBCIDRangUpper() >= bctype )
    {
        for ( int nZone = 0; nZone < outLetBCs.GetZoneNumberOfSameBCType(); ++ nZone )
        {
            SubnicBC BCS;
            BCS = outLetBCs.GetBCInfo()[nZone];

            if ( bctype == BCS.GetBCID() )
            {
                OutFlowBC3D(gridIn, bcr, BCS);
                break;
            }
        }
    }
    else
    {
        hasSubsonicBC = 0;
    }
    return hasSubsonicBC;
}

bool SetSubsonicBC(UnstructGrid * gridIn, CompositeUnstructBC ** bcr, int bctype, RDouble * primSource, RDouble * primTarget,
                   RDouble nxs, RDouble nys, RDouble nzs, RDouble t)
{
    bool hasSubsonicBC = 1;

    if (inLetBCs.GetBCIDRangLower() <= bctype &&  inLetBCs.GetBCIDRangUpper() >= bctype)
    {
        for (int nZone = 0; nZone < inLetBCs.GetZoneNumberOfSameBCType(); ++ nZone)
        {
            SubnicBC BCS;
            BCS = inLetBCs.GetBCInfo()[nZone];

            if (bctype == BCS.GetBCID())
            {
                InFlowBC3D(gridIn, bcr, BCS, primSource, primTarget, nxs, nys, nzs, t);
                break;
            }
        }
    }
    else if(outLetBCs.GetBCIDRangLower() <= bctype && outLetBCs.GetBCIDRangUpper() >= bctype)
    {
        for (int nZone = 0; nZone < outLetBCs.GetZoneNumberOfSameBCType(); ++ nZone)
        {
            SubnicBC BCS;
            BCS = outLetBCs.GetBCInfo()[nZone];

            if (bctype == BCS.GetBCID())
            {
                OutFlowBC3D(gridIn, bcr, BCS, primSource, primTarget, nxs, nys, nzs, t);
                break;
            }
        }
    }
    else
    {
        hasSubsonicBC = 0;
    }
    return hasSubsonicBC;
}

bool SetSubsonicBCTurb( int bctype, RDouble **q_turb, RDouble * turboo, int n_turb, int le, int re )
{
    bool hasSubsonicBC = 1;

    if ( inLetBCs.GetBCIDRangLower() <= bctype &&  inLetBCs.GetBCIDRangUpper() >= bctype )
    {
        for ( int nZone = 0; nZone < inLetBCs.GetZoneNumberOfSameBCType(); ++ nZone )
        {
            SubnicBC BCS;
            BCS = inLetBCs.GetBCInfo()[nZone];

            if ( bctype == BCS.GetBCID() )
            {
                InFlowBC3DTurb( bctype, q_turb, turboo, n_turb, le, re );
                break;
            }
        }
    }
    else if( outLetBCs.GetBCIDRangLower() <= bctype && outLetBCs.GetBCIDRangUpper() >= bctype )
    {
        for ( int nZone = 0; nZone < outLetBCs.GetZoneNumberOfSameBCType(); ++ nZone )
        {
            SubnicBC BCS;
            BCS = outLetBCs.GetBCInfo()[nZone];

            if ( bctype == BCS.GetBCID() )
            {
                OutFlowBC3DTurb( bctype, q_turb, turboo, n_turb, le, re );
                break;
            }
        }
    }
    else
    {
        hasSubsonicBC = 0;
    }
    return hasSubsonicBC;
}

void OutFlowBC3D(Grid *gridIn, StructBC *bcr, SubnicBC & outLet)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    int nl = 5;
    GlobalDataBase::GetData("nl", &nl, PHINT, 1);

    int nchem = 0;
    GlobalDataBase::GetData("nchem", &nchem, PHINT, 1);

    int neqn = nl + nchem;

    RDouble4D &q = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));

    int ist, ied, jst, jed, kst, ked;
    bcr->GetIJKRegion(ist, ied, jst, jed, kst, ked);

    int is, js, ks, it, jt, kt;

    RDouble pb = outLet.GetTotalPressure();

    for ( int k = kst; k <= ked; ++ k )
    {
        for ( int j = jst; j <= jed; ++ j )
        {
            for ( int i = ist; i <= ied; ++ i )
            {
                bcr->GetInsideCellIndex(i, j, k, is, js, ks, 1);

                RestrictIndex(is, js, ks, ni-1, nj-1, nk-1);

                for (int layer = 1; layer <= 2; ++ layer)
                {
                    bcr->GetGhostCellIndex(i, j, k, it, jt, kt, layer);

                    for ( int m = 0; m < neqn; ++ m )
                    {
                        q(it,jt,kt,m) = q(is,js,ks,m);
                    }
                    using namespace IDX;
                    q(it,jt,kt,IP) = pb;
                }
            }
        }
    }
}

void OutFlowBC3D(UnstructGrid * gridIn, CompositeUnstructBC ** bcr, SubnicBC & outLet, 
                 RDouble * primSource, RDouble * primTarget, 
                 RDouble nxs, RDouble nys, RDouble nzs, RDouble t)
{
    int nl    = GlobalDataBase::GetIntParaFromDB("nl");
    int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
    int neqn = nl + nchem;

    RDouble pb = outLet.GetTotalPressure();

    for (int m = 0; m < neqn; ++ m)
    {
        primTarget[m] = primSource[m];
    }
    using namespace IDX;
    primTarget[IP] = pb;
}

void OutFlowBC3DTurb( int bctype, RDouble **q_turb, RDouble * turboo, int n_turb, int le, int re )
{
    for ( int m = 0; m < n_turb; ++ m )
    {
        q_turb[m][re] = q_turb[m][le];
    }
}

void InFlowBC3D(Grid *gridIn, StructBC *bcr, SubnicBC & inLet)
{
    StructGrid *grid = StructGridCast(gridIn);

    RDouble4D &q = * reinterpret_cast< RDouble4D * > ( grid->GetDataPtr("q") );
    RDouble4D &t = * reinterpret_cast< RDouble4D * > ( grid->GetDataPtr("t") );

    int ist, ied, jst, jed, kst, ked;
    bcr->GetIJKRegion(ist, ied, jst, jed, kst, ked);

    RDouble4D &xfn = *(grid->GetFaceNormalX());
    RDouble4D &yfn = *(grid->GetFaceNormalY());
    RDouble4D &zfn = *(grid->GetFaceNormalZ());

    int it, jt, kt, in, jn, kn;
    int nsurf = bcr->GetFaceDirection() + 1;
    int lr    = bcr->GetFaceLeftOrRightIndex();

    RDouble nxs, nys, nzs;

    RDouble p0, t0;
    RDouble * fd;

    p0 = inLet.GetTotalPressure();
    t0 = inLet.GetTotalTemperature();
    fd = inLet.GetFlowDirection();

    RDouble refGama = 1.4, gama1, moo, gmogm1;
    GlobalDataBase::GetData("refGama", &refGama, PHDOUBLE, 1);
    GlobalDataBase::GetData("refMachNumber", &moo, PHDOUBLE, 1);
    gama1  = refGama - 1.0;
    gmogm1 = refGama / gama1;

    RDouble r1, u1, v1, w1, p1, t1;
    RDouble uvw, un1, cosq, qcta, a1, Rneg, ab, tb, pb, rb, uvwb, ub, vb, wb;

    RDouble a0 = sqrt(t0) / moo;
    RDouble a02 = a0 * a0;
    RDouble gmamo2 = refGama * moo * moo;

    for ( int k = kst; k <= ked; ++ k )
    {
        for ( int j = jst; j <= jed; ++ j )
        {
            for ( int i = ist; i <= ied; ++ i )
            {
                bcr->GetBoundaryFaceIndex(i, j, k, in, jn, kn);

                nxs = -lr * xfn(in, jn, kn, nsurf);
                nys = -lr * yfn(in, jn, kn, nsurf);
                nzs = -lr * zfn(in, jn, kn, nsurf);

                using namespace IDX;

                r1 = q(i,j,k,IR);
                u1 = q(i,j,k,IU);
                v1 = q(i,j,k,IV);
                w1 = q(i,j,k,IW);
                p1 = q(i,j,k,IP);
                t1 = t(i,j,k,ITT);

                uvw = sqrt( u1 * u1 + v1 * v1 + w1 * w1 );
                un1 = u1 * nxs + v1 * nys + w1 * nzs;

                cosq = - un1 / uvw;

                qcta = cosq * cosq + 2.0 / gama1;

                a1 = sqrt( t1 ) / moo;

                Rneg = un1 - 2.0 * a1 / gama1;

                ab = ( -Rneg - cosq * sqrt( qcta * a02 - gama1 / 2.0 * Rneg * Rneg ) ) / qcta;

                tb = t0 * pow( ( ab / a0 ), 2.0 );
                pb = p0 * pow( ( tb / t0 ), gmogm1 );
                rb = pb / tb * gmamo2;

                uvwb = sqrt( 2.0 * ( t0 - tb) / gama1 ) / moo;

                ub = uvwb * fd[0];
                vb = uvwb * fd[1];
                wb = uvwb * fd[2];

                for (int layer = 1; layer <= 2; ++ layer)
                {
                    bcr->GetGhostCellIndex(i, j, k, it, jt, kt, layer);

                    q(it,jt,kt,IR) = rb;
                    q(it,jt,kt,IU) = ub;
                    q(it,jt,kt,IV) = vb;
                    q(it,jt,kt,IW) = wb;
                    q(it,jt,kt,IP) = pb;
                }
            }
        }
    }
}

void InFlowBC3D(UnstructGrid * gridIn, CompositeUnstructBC ** bcr, SubnicBC & inLet, RDouble * primSource, RDouble * primTarget,
                RDouble nxs, RDouble nys, RDouble nzs, RDouble t)
{
    RDouble p0 = 0.0,t0 = 0.0;
    RDouble * fd;

    p0 = inLet.GetTotalPressure();
    t0 = inLet.GetTotalTemperature();
    fd = inLet.GetFlowDirection();

    RDouble refGama = 1.4;
    RDouble gama1, moo, gmogm1;
    GlobalDataBase::GetData("refGama", &refGama, PHDOUBLE, 1);
    GlobalDataBase::GetData("refMachNumber", &moo, PHDOUBLE, 1);
    gama1  = refGama - 1.0;
    gmogm1 = refGama / gama1;

    RDouble r1, u1, v1, w1, p1, t1;
    RDouble uvw, un1, cosq, qcta, a1, Rneg, ab, tb, pb, rb, uvwb, ub, vb, wb;

    RDouble a0 = sqrt(t0) / moo;
    RDouble a02 = a0 * a0;
    RDouble gmamo2 = refGama * moo * moo;

    using namespace IDX;
    r1 = primSource[IR];
    u1 = primSource[IU];
    v1 = primSource[IV];
    w1 = primSource[IW];
    p1 = primSource[IP];
    t1 = t;

    uvw = sqrt(u1 * u1 + v1 * v1 + w1 * w1);
    un1 = u1 * nxs + v1 * nys + w1 * nzs;

    cosq = - un1 / uvw;

    qcta = cosq * cosq + 2.0 / gama1;

    a1 = sqrt( t1 ) / moo;

    Rneg = un1 - 2.0 * a1 / gama1;

    ab = (-Rneg - cosq * sqrt(qcta * a02 - gama1 / 2.0 * Rneg * Rneg)) / qcta;

    tb = t0 * pow((ab / a0), 2.0);
    pb = p0 * pow((tb / t0), gmogm1);
    rb = pb / tb * gmamo2;

    uvwb = sqrt(2.0 * (t0 - tb) / gama1) / moo;

    ub = uvwb * fd[0];
    vb = uvwb * fd[1];
    wb = uvwb * fd[2];

    primTarget[IR] = rb;
    primTarget[IU] = ub;
    primTarget[IV] = vb;
    primTarget[IW] = wb;
    primTarget[IP] = pb;
}

void InFlowBC3DTurb( int bctype, RDouble **q_turb, RDouble * turboo, int n_turb, int le, int re )
{
    for ( int m = 0; m < n_turb; ++ m )
    {
        q_turb[m][re] = turboo[m];
    }
}


}
