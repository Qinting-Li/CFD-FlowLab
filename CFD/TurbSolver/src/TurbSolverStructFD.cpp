#include "TurbSolverStructFD.h"
#include "Param_TurbSolverStruct.h"
#include "JetBC.h"
#include "TK_Exit.h"
#include "Glb_Dimension.h"
#include "GradientOperation.h"

namespace PHSPACE
{
TurbSolverStrFD::TurbSolverStrFD()
{
}

TurbSolverStrFD::~TurbSolverStrFD()
{
}

void TurbSolverStrFD::Boundary(Grid *gridIn)
{
}

void TurbSolverStrFD::ComputeViscousCoeff(Grid *grid)
{
}

void TurbSolverStrFD::ComputeGradientCellCenter(Grid *gridIn)
{
}

void TurbSolverStrFD::GetDependentVariablesforStructHighOrder(Grid *gridIn)
{
    Param_TurbSolverStruct *parameters = GetControlParameters();
    int viscousType = parameters->GetViscousType();

    ObtainViscousCoefWithGhost(gridIn);
    ObtainqlqrLowOrder(gridIn);
    ObtainGradientCellCenter_Laminar(gridIn);

    ObtainViscosity(gridIn);
    ObtainBoundaryValue(gridIn);
    ObtainGradientCellCenter(gridIn);

    if(viscousType == 4)
    {
        ObtainCrossing(gridIn);
        ObtainBlending(gridIn);
    }
}

void TurbSolverStrFD::ObtainViscousCoefWithGhost(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    RDouble3D &visl = *reinterpret_cast <RDouble3D *> (grid->GetDataPtr("visl"));
    RDouble4D &q_laminar = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("q"));    

    Param_TurbSolverStruct *parameters = GetControlParameters();
    RDouble refMachNumber = parameters->GetRefMachNumber();
    const RDouble refGama = 1.4;
    RDouble nonDimensionalSutherlandTemperature = GlobalDataBase::GetDoubleParaFromDB("tsuth");

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    int ist, ied, jst, jed, kst, ked;
    GetIJKRegion(I, J, K, ist, ied, jst, jed, kst, ked);

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble rs = q_laminar(i, j, k, 0);
                RDouble ps = q_laminar(i, j, k, 4);
                RDouble tm = refGama * refMachNumber * refMachNumber * ps / rs;
                visl(i, j, k) = tm * sqrt(tm) * (1.0 + nonDimensionalSutherlandTemperature) / (tm + nonDimensionalSutherlandTemperature);
            }
        }
    }
}

void TurbSolverStrFD::ObtainqlqrLowOrder(Grid *gridIn)
{
    StructGrid * grid = StructGridCast(gridIn);

    FieldProxy *qlProxy = new FieldProxy();
    FieldProxy *qrProxy = new FieldProxy();

    int nDim = GetDim();
	for (int nsurf = 1; nsurf <= nDim; ++ nsurf)
	{
        if (nsurf == 1)
        {
            qlProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("ql1")));
            qrProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("qr1")));
        }
        else if (nsurf == 2)
        {
            qlProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("ql2")));
            qrProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("qr2")));
        }
        else
        {
            qlProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("ql3")));
            qrProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("qr3")));
        }
        
        GetFaceVariable(grid, qlProxy, qrProxy, nsurf);
        GetFaceVariableCorrectionAtPhysicalBoundary(grid, qlProxy, qrProxy, nsurf);
    }

    delete qlProxy;
    delete qrProxy;
}

void TurbSolverStrFD::GetFaceVariable(Grid *gridIn, FieldProxy *qlProxy, FieldProxy *qrProxy, int nsurf)
{
    StructGrid * grid = StructGridCast(gridIn);
	RDouble4D & q = * reinterpret_cast< RDouble4D * > (grid->GetDataPtr("q"));
	RDouble4D & ql = qlProxy->GetField_STR();
	RDouble4D & qr = qrProxy->GetField_STR();

	int ni = grid->GetNI();
	int nj = grid->GetNJ();
	int nk = grid->GetNK();

    int il1, jl1, kl1;
    GetNsurfIndex(nsurf, il1, jl1, kl1);
    
    int ist = 1;
    int jst = 1;
    int kst = 1;
	int ied = ni - 1 + il1;
	int jed = nj - 1 + jl1;
	int ked = nk - 1 + kl1;
	if (nk == 1)
    {
        kst = 1;
        ked = 1;
    }   

    RDouble PlateL[3];
    RDouble PlateR[3];
    
    for (int m = 1; m <= 3; ++m)
    {
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    for (int ip = 0; ip <= 2; ++ ip)
                    {
                        int iL = i + (ip - 2) * il1;
                        int jL = j + (ip - 2) * jl1;
                        int kL = k + (ip - 2) * kl1;

                        int iR = i + (1 - ip) * il1;
                        int jR = j + (1 - ip) * jl1;
                        int kR = k + (1 - ip) * kl1;
                        
                        PlateL[ip] = q(iL, jL, kL, m);
                        PlateR[ip] = q(iR, jR, kR, m);
                    }
                    ql(i - il1, j - jl1, k - kl1, m) = PlateL[1] + half * Limiter_ThirdSmooth(PlateL[2] - PlateL[1], PlateL[1] - PlateL[0]);
                    qr(i,       j,       k,       m) = PlateR[1] + half * Limiter_ThirdSmooth(PlateR[2] - PlateR[1], PlateR[1] - PlateR[0]);
                }
            }
        }
    }
}

RDouble TurbSolverStrFD::Limiter_ThirdSmooth(const RDouble &x, const RDouble &y)
{
   const RDouble eps2 = 1.0e-6;
	RDouble Smooth1 = x * ( y * y + 2.0 * eps2 ) +  y * ( 2.0 * x * x + eps2 );
	RDouble Smooth2 = 2.0 * x * x - x * y + 2.0 * y * y + 3.0 * eps2;
	RDouble Smooth = Smooth1 / Smooth2;
    return  Smooth;
}

void TurbSolverStrFD::GetFaceVariableCorrectionAtPhysicalBoundary(Grid *gridIn, FieldProxy *qlProxy, FieldProxy *qrProxy, int nsurf)
{
    StructGrid *grid = StructGridCast(gridIn);
    RDouble4D & ql = qlProxy->GetField_STR();
    RDouble4D & qr = qrProxy->GetField_STR();
    
    CompositeBCRegion *compositeBCRegion = grid->GetCompositeBCRegion();
    int nBCRegion = compositeBCRegion->GetnBCRegion();
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        BCRegion *bcregion = compositeBCRegion->GetBCRegion(iBCRegion);
        int BCType = bcregion->GetBCType();
        int nsurf_bc = bcregion->GetFaceDirection() + 1;
        if (!IsWall(BCType) || nsurf_bc != nsurf)
        {
            continue;
        }
        int ist, ied, jst, jed, kst, ked;
        bcregion->GetIJKRegion(ist, ied, jst, jed, kst, ked);
        
        int il1, jl1, kl1;
        GetNsurfIndex(nsurf, il1, jl1, kl1);

        int *s_lr3d = bcregion->GetFaceDirectionIndex();
        int id, jd, kd;
        GetBCFaceIDX(s_lr3d, id, jd, kd);

        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    int ibc = i + id;
                    int jbc = j + jd;
                    int kbc = k + kd;
                    
                    qr(ibc, jbc, kbc, 1) = 0.0;
                    qr(ibc, jbc, kbc, 2) = 0.0;
                    qr(ibc, jbc, kbc, 3) = 0.0;
                    
                    ql(ibc - il1, jbc - jl1, kbc - kl1, 1) = 0.0;
                    ql(ibc - il1, jbc - jl1, kbc - kl1, 2) = 0.0;
                    ql(ibc - il1, jbc - jl1, kbc - kl1, 3) = 0.0;
                }
            }
        }
    }
}

void TurbSolverStrFD::ObtainGradientCellCenter_Laminar(Grid *gridIn)
{
    StructGrid *strgrid = StructGridCast(gridIn);
    int ni = strgrid->GetNI();
    int nj = strgrid->GetNJ();
    int nk = strgrid->GetNK();

    RDouble4D &xfv = *(strgrid->GetFaceVectorX());
    RDouble4D &yfv = *(strgrid->GetFaceVectorY());
    RDouble4D &zfv = *(strgrid->GetFaceVectorZ());
    RDouble3D &volume      = *(strgrid->GetCellVolume());

    int ndim = GetDim();

    RDouble4D &q = *reinterpret_cast<RDouble4D *> (strgrid->GetDataPtr("q"));
    RDouble4D &dqdx = *reinterpret_cast<RDouble4D *> (strgrid->GetDataPtr("gradUVWTCellCenterX"));
    RDouble4D &dqdy = *reinterpret_cast<RDouble4D *> (strgrid->GetDataPtr("gradUVWTCellCenterY"));
    RDouble4D &dqdz = *reinterpret_cast<RDouble4D *> (strgrid->GetDataPtr("gradUVWTCellCenterZ"));

    dqdx = 0.0;
    dqdy = 0.0;
    dqdz = 0.0;

    for (int iSurface = 1; iSurface <= ndim; ++ iSurface)
    {
        int il1, jl1, kl1;
        GetNsurfIndex(iSurface, il1, jl1, kl1);

        int ist = 1;
        int ied = ni - 1 + il1;
        int jst = 1;
        int jed = nj - 1 + jl1;
        int kst = 1;
        int ked = nk - 1 + kl1;

        if (ndim == TWO_D)
        {
            ked = 1;
        }
        int il, jl, kl;

        //! The gradient of u,v,w.
        for (int m = 0; m < ndim; ++ m)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    for (int i = ist; i <= ied; ++ i)
                    {
                        il = i - il1;
                        jl = j - jl1;
                        kl = k - kl1;

                        RDouble phis = q(i,j,k,m+1) + q(il,jl,kl,m+1);

                        RDouble ddx = phis * xfv(i, j, k, iSurface);
                        RDouble ddy = phis * yfv(i, j, k, iSurface);
                        RDouble ddz = phis * zfv(i, j, k, iSurface);

                        dqdx(i ,j ,k ,m) -= ddx;
                        dqdy(i ,j ,k ,m) -= ddy;
                        dqdz(i ,j ,k ,m) -= ddz;

                        dqdx(il, jl, kl, m) += ddx;
                        dqdy(il, jl, kl, m) += ddy;
                        dqdz(il, jl, kl, m) += ddz;
                    }
                }
            }
        }
    }

    int ist = 1;
    int ied = ni-1;
    int jst = 1;
    int jed = nj-1;
    int kst = 1;
    int ked = nk-1;

    if (ndim == TWO_D) ked = 1;

    for (int m = 0; m < ndim; ++ m)
    {
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    RDouble oov = half / volume(i, j, k);
                    dqdx(i, j, k, m) *= oov;
                    dqdy(i, j, k, m) *= oov;
                    dqdz(i, j, k, m) *= oov;
                }
            }
        }
    }
}


void TurbSolverStrFD::ObtainBoundaryValue(Grid *gridIn)
{
    using namespace PHENGLEI;
    StructGrid *grid = StructGridCast(gridIn);

    Param_TurbSolverStruct *parameters = GetControlParameters();
    int wallFunctionType = parameters->GetWallFunctionType();

    CompositeBCRegion *compositeBCRregion = grid->GetCompositeBCRegion();
    int nBCRegion = compositeBCRregion->GetnBCRegion();

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        BCRegion *BCRegion = compositeBCRregion->GetBCRegion(iBCRegion);
        int BCType = BCRegion->GetBCType();

        if (BCType < 0)
        {
            continue;
        }
        else if (BCType == EXTRAPOLATION)
        {
            OutFlowBC(grid, BCRegion);
        }
        else if (IsWall(BCType))
        {
            if (wallFunctionType ==0)
            {
                VisWall(grid, BCRegion);
            }
            else if(wallFunctionType == 1)
            { 
                VisWallWithWallFunctionStandard(grid, BCRegion);
            }
            else if (wallFunctionType == 2)
            {
                VisWallWithWallFunctionPAB3D(grid, BCRegion);
            }
        }
        else if (BCType == SYMMETRY)
        {
            SymmetryBC(grid, BCRegion);
        }
        else if (BCType == FARFIELD)
        {
            FarFieldBC(grid, BCRegion);
        }
        else if (BCType == INFLOW)
        {
            InFlowBC(grid, BCRegion);
        }
        else if (BCType == OUTFLOW)
        {
            OutFlowBC(grid, BCRegion);
        }
        else if (BCType == POLE || BCType / 10 == POLE)
        {
            OutFlowBC(grid, BCRegion);
        }

        else if (BCType == EXTERNAL_BC)
        {
            OutFlowBC(grid, BCRegion);
        }
        //! Deal in symmetry wall method if there is no boundary ID specified.
        else if (ifIsJetBC(BCType))
        {
            SymmetryBC(grid, BCRegion);
        }
        else
        {
            ostringstream oss;
            oss << "Error : Illegal BCtype ID " << BCType << endl;
            TK_Exit::ExceptionExit(oss);
        }
    }
    CornerPoint(grid);
}

void TurbSolverStrFD::ObtainViscosity(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    RDouble3D &wallDistant = *grid->GetWallDist();

    RDouble4D &area   = *(grid->GetFaceArea()   );
    RDouble3D &volume = *(grid->GetCellVolume() );

    RDouble4D &qLaminar         = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q" ));
    RDouble3D &viscousLaminar   = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("visl"));
    RDouble3D &viscousTurbulent = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("vist"));
    RDouble4D &qTurbulent       = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb" ));
    RDouble3D &SST_F2           = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("SST_F2"));

    Param_TurbSolverStruct *parameters = GetControlParameters();
    RDouble eddyViscosityLimit = parameters->GetEddyViscosityLimit();

    int monitorVistmax = GlobalDataBase::GetIntParaFromDB("monitor_vistmax");

    string viscousName = parameters->GetViscousName();

    RDouble refReNumber = parameters->GetRefReNumber();
    RDouble reynoldsSquare = refReNumber * refReNumber;

    RDouble viscousTurbulentMaximum = 0.0;
    RDouble viscousTurbulentMinimum = 1.0e30;
    RDouble logsmall = log(SMALL);

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    int imax = 1;
    int jmax = 1;
    int kmax = 1;

    using namespace IDX;

    if (viscousName.substr(0,6) == "1eq-sa")
    {
        RDouble SA_cv1_cube = parameters->GetSA_cv1_cube();

        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    RDouble ld   = qTurbulent(i, j, k, ISA) * qLaminar(i, j, k, IR) / (viscousLaminar(i, j, k) + SMALL);
                    RDouble ld3  = ld * ld * ld;
                    RDouble fv1  = ld3 / (ld3 + SA_cv1_cube);
                    viscousTurbulent(i, j, k) = qLaminar(i, j, k, IR) * fv1 * qTurbulent(i, j, k, ISA);
                    viscousTurbulent(i, j, k) = MIN(eddyViscosityLimit, viscousTurbulent(i, j, k));
                    if (viscousTurbulentMaximum < viscousTurbulent(i, j, k))
                    {
                        viscousTurbulentMaximum =  viscousTurbulent(i, j, k);
                        imax = i;
                        jmax = j;
                        kmax = k;
                    }
                }
            }
        }
    }
    else if (viscousName.substr(0,6) == "2eq-kw")
    {
        //! Two equation models.
        if (viscousName.substr(0,17) == "2eq-kw-menter-sst")
        {
            //! Original Menter k-w-SST model.
            RDouble SST_a1    = parameters->GetSST_a1();
            RDouble SST_betaStar = parameters->GetSST_betaStar();

            Range II, JJ, KK;
            GetRange(ni, nj, nk, -2, 1, II, JJ, KK);

            RDouble4D &gradUVWTCellCenterX = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterX"));
            RDouble4D &gradUVWTCellCenterY = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterY"));
            RDouble4D &gradUVWTCellCenterZ = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterZ"));

            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    for (int i = ist; i <= ied; ++ i)
                    {
                        RDouble rho = ABS(qLaminar(i, j, k, IR)) + SMALL;
                        RDouble ke  = qTurbulent(i, j, k, IDX::IKE);
                        RDouble kw  = qTurbulent(i, j, k, IDX::IKW);
                        RDouble ds  = wallDistant(i, j, k);
                        RDouble ds2 = ds * ds;

                        RDouble dudx = gradUVWTCellCenterX(i, j, k, 0);
                        RDouble dudy = gradUVWTCellCenterY(i, j, k, 0);
                        RDouble dudz = gradUVWTCellCenterZ(i, j, k, 0);

                        RDouble dvdx = gradUVWTCellCenterX(i, j, k, 1);
                        RDouble dvdy = gradUVWTCellCenterY(i, j, k, 1);
                        RDouble dvdz = gradUVWTCellCenterZ(i, j, k, 1);

                        RDouble dwdx = gradUVWTCellCenterX(i, j, k, 2);
                        RDouble dwdy = gradUVWTCellCenterY(i, j, k, 2);
                        RDouble dwdz = gradUVWTCellCenterZ(i, j, k, 2);

                        //RDouble vort  = sqrt((dwdy - dvdz) * (dwdy - dvdz) +
                        //    (dudz - dwdx) * (dudz - dwdx) +
                        //    (dvdx - dudy) * (dvdx - dudy));

                        RDouble s11 = dudx; 
                        RDouble s22 = dvdy; 
                        RDouble s33 = dwdz;
                        RDouble s12 = half * (dudy + dvdx); 
                        RDouble s13 = half * (dudz + dwdx); 
                        RDouble s23 = half * (dvdz + dwdy);

                        RDouble sij2 = two * (s11*s11 + s22*s22 + s33*s33 + two*(s12*s12 + s13*s13 + s23*s23));  //! modulus of S.
                        RDouble Strain = sqrt(sij2);

                        RDouble part1 = two * sqrt(ke) / (SST_betaStar * kw * ds * refReNumber);
                        RDouble part2 = 500.0 * viscousLaminar(i, j, k) / (rho * kw * ds2 * reynoldsSquare);
                        RDouble arg2  = MAX(part1, part2);
                        RDouble f2    = tanh(static_cast<RDouble>(arg2 * arg2));
                        SST_F2(i, j, k) = f2;

                        RDouble &turbulentViscosity = viscousTurbulent(i, j, k);
                        turbulentViscosity = rho * ke / MAX(kw, Strain * f2 / (SST_a1 * refReNumber));
                        turbulentViscosity = MIN(eddyViscosityLimit, turbulentViscosity);
                        turbulentViscosity = MAX(1.E-5,turbulentViscosity);

                        if (viscousTurbulentMaximum < viscousTurbulent(i, j, k))
                        {
                            viscousTurbulentMaximum = turbulentViscosity;
                            imax = i;
                            jmax = j;
                            kmax = k;
                        }
                    }
                }
            }
        }
    }

    RDouble vist_max = viscousTurbulentMaximum;
    RDouble vist_min = viscousTurbulentMinimum;

    grid->UpdateData("vist_max", &vist_max, PHDOUBLE, 1);
    grid->UpdateData("vist_min", &vist_min, PHDOUBLE, 1);

#ifdef PH_PARALLEL
    using namespace PHMPI;
    int myid = GetCurrentProcessorID();
    if (myid == server)
    {
#endif
        if (monitorVistmax)
        {
            cout << "vistmax = " << viscousTurbulentMaximum << " " << imax << " " << jmax << " " << kmax << "\n";
        }
#ifdef PH_PARALLEL
    }
#endif
}

void TurbSolverStrFD::ObtainGradientCellCenter(Grid *gridIn)
{
    StructGrid *strgrid = StructGridCast(gridIn);
    int ni = strgrid->GetNI();
    int nj = strgrid->GetNJ();
    int nk = strgrid->GetNK();

    RDouble4D &faceVectorX = *(strgrid->GetFaceVectorX());
    RDouble4D &faceVectorY = *(strgrid->GetFaceVectorY());
    RDouble4D &faceVectorZ = *(strgrid->GetFaceVectorZ());
    RDouble3D &volume      = *(strgrid->GetCellVolume());

    int ndim = GetDim();

    Param_TurbSolverStruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    RDouble4D &qTurbulent    = *reinterpret_cast<RDouble4D *> (strgrid->GetDataPtr("q_turb"));
    RDouble4D &gradTurbulentCellCenterX = *reinterpret_cast<RDouble4D *> (strgrid->GetDataPtr("gradTurbulentCellCenterX"));
    RDouble4D &gradTurbulentCellCenterY = *reinterpret_cast<RDouble4D *> (strgrid->GetDataPtr("gradTurbulentCellCenterY"));
    RDouble4D &gradTurbulentCellCenterZ = *reinterpret_cast<RDouble4D *> (strgrid->GetDataPtr("gradTurbulentCellCenterZ"));

    gradTurbulentCellCenterX = 0.0;
    gradTurbulentCellCenterY = 0.0;
    gradTurbulentCellCenterZ = 0.0;

    for (int iSurface = 1; iSurface <= ndim; ++ iSurface)
    {
        int il1, jl1, kl1;
        GetNsurfIndex(iSurface, il1, jl1, kl1);

        int ist = 1;
        int ied = ni - 1 + il1;
        int jst = 1;
        int jed = nj - 1 + jl1;
        int kst = 1;
        int ked = nk - 1 + kl1;

        if (ndim == TWO_D)
        {
            ked = 1;
        }
        int il, jl, kl;

        //! The gradient of u,v,w.
        for (int m = 0; m < nTurbulentEquation; ++ m)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    for (int i = ist; i <= ied; ++ i)
                    {
                        il = i - il1;
                        jl = j - jl1;
                        kl = k - kl1;

                        RDouble phis = qTurbulent(i, j, k, m) + qTurbulent(il, jl, kl, m);

                        RDouble ddx = phis * faceVectorX(i, j, k, iSurface);
                        RDouble ddy = phis * faceVectorY(i, j, k, iSurface);
                        RDouble ddz = phis * faceVectorZ(i, j, k, iSurface);

                        RDouble &gradTurbccXLeft = gradTurbulentCellCenterX(i, j, k, m);
                        RDouble &gradTurbccYLeft = gradTurbulentCellCenterY(i, j, k, m);
                        RDouble &gradTurbccZLeft = gradTurbulentCellCenterZ(i, j, k, m);
                        gradTurbccXLeft -= ddx;
                        gradTurbccYLeft -= ddy;
                        gradTurbccZLeft -= ddz;

                        RDouble &gradTurbccXRight = gradTurbulentCellCenterX(il, jl, kl, m);
                        RDouble &gradTurbccYRight = gradTurbulentCellCenterY(il, jl, kl, m);
                        RDouble &gradTurbccZRight = gradTurbulentCellCenterZ(il, jl, kl, m);
                        gradTurbccXRight += ddx;
                        gradTurbccYRight += ddy;
                        gradTurbccZRight += ddz;
                    }
                }
            }
        }
    }

    int ist = 1;
    int ied = ni-1;
    int jst = 1;
    int jed = nj-1;
    int kst = 1;
    int ked = nk-1;

    if (ndim == TWO_D) ked = 1;

    for (int m = 0; m < nTurbulentEquation; ++ m)
    {
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    RDouble oov = half / volume(i, j, k);
                    RDouble &gradTurbccX = gradTurbulentCellCenterX(i, j, k, m);
                    RDouble &gradTurbccY = gradTurbulentCellCenterY(i, j, k, m);
                    RDouble &gradTurbccZ = gradTurbulentCellCenterZ(i, j, k, m);
                    gradTurbccX *= oov;
                    gradTurbccY *= oov;
                    gradTurbccZ *= oov;
                }
            }
        }
    }
}

void TurbSolverStrFD::ObtainCrossing(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);

    RDouble4D &qTurbulent = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb"));
    RDouble3D &cross = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("cross"));

    RDouble4D &faceNormalComponentX = *(grid->GetFaceNormalX());
    RDouble4D &faceNormalComponentY = *(grid->GetFaceNormalY());
    RDouble4D &faceNormalComponentZ = *(grid->GetFaceNormalZ());
    RDouble4D &area   = *(grid->GetFaceArea()   );
    RDouble3D &volume = *(grid->GetCellVolume() );

    GreenGeo3D greenGeo3D;

    int ist, ied, jst, jed, kst, ked;
	grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    int ndim = GetDim();
    int js = 1;
    int ks = 1;
    if (ndim == 2)
    {
        ks = 0;
    }
    if (ndim == 1)
    {
        js = 0;
    }

    RDouble dkedx,dkedy,dkedz,dkwdx,dkwdy,dkwdz;
    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                Get_SPM(i, j, k, js, ks, ndim, &greenGeo3D, faceNormalComponentX, faceNormalComponentY, faceNormalComponentZ, area, volume);
                DXDYDZ(qTurbulent, &greenGeo3D, i, j, k, js, ks, IDX::IKE, dkedx, dkedy, dkedz);
                DXDYDZ(qTurbulent, &greenGeo3D, i, j, k, js, ks, IDX::IKW, dkwdx, dkwdy, dkwdz);
                RDouble &crossCell = cross(i, j, k);
                crossCell = dkedx * dkwdx + dkedy * dkwdy + dkedz * dkwdz;
            }
        }
    }
}

void TurbSolverStrFD::ObtainBlending(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    RDouble4D &qLaminar       = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q" ));
    RDouble4D &qTurbulent     = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb" ));
    RDouble3D &viscousLaminar = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("visl"));

    RDouble3D &cross = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("cross"));
    RDouble3D &blend = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("blend"));

    RDouble3D &walldist = * StructGridCast(grid)->GetWallDist();

    Param_TurbSolverStruct *parameters = GetControlParameters();
    RDouble refReNumber  = parameters->GetRefReNumber();
    RDouble oRefReNumber = parameters->GetoRefReNumber();
    RDouble oRefReNumberSqr = oRefReNumber * oRefReNumber;

    RDouble betas = 0.09;
    RDouble sigw2 = 0.856;

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    using namespace IDX;

    //Compute maximum cross diffusion term across flowfield
    RDouble cdkwmax = 0.0;
    RDouble cdkwmin = 1.0e-10;

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble rho  = qLaminar(i, j, k,IR);
                RDouble ke   = qTurbulent(i, j, k, IKE);
                RDouble kw   = qTurbulent(i, j, k, IKW);
                RDouble crss = cross(i, j, k);
 
                RDouble miuLam = viscousLaminar(i, j, k);
                RDouble dist   = walldist(i, j, k);
                RDouble dist2  = dist * dist;

                //calculate cd_kw
                RDouble crossDiff = 2.0 * rho * sigw2 * crss / (kw + SMALL);

                //! Original Menter CD_kw calculation.
                RDouble cdkw = MAX(crossDiff, cdkwmin);

                //! Calculate arg1.
                RDouble term1 = sqrt(ABS(ke)) * oRefReNumber /(betas * dist * kw + SMALL);
                RDouble term2 = 500.0 * miuLam * oRefReNumberSqr / (rho * kw * dist2 + SMALL);
                RDouble term3 = MAX(term1, term2);
                RDouble term4 = 4.0 * rho * sigw2 * ke / (cdkw * dist2 + SMALL);
                RDouble arg1  = MIN(term3, term4);
     
                //! Compute the blending function fbsl.
                RDouble &blendCell = blend(i, j, k);
                blendCell = tanh(static_cast<RDouble>(arg1 * arg1 * arg1 * arg1));
            }
        }
    }
    GhostCell3D(blend, ni, nj, nk);
}


void TurbSolverStrFD::CornerPoint(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Param_TurbSolverStruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    RDouble4D &qTurbulent = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb"));
    FillCornerPoint3D(qTurbulent, ni, nj, nk, nTurbulentEquation);
}

void TurbSolverStrFD::InFlowBC(Grid *grid, BCRegion *BCRegion)
{
    Param_TurbSolverStruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();
    RDouble4D &qTurbulent = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb"));
    RDouble   *freeStreamTurbVar =  parameters->GetFreeStreamTurbVar();
    RDouble3D &viscousTurbulent = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("vist"));
    RDouble freeStreamViscosity = parameters->GetFreeStreamViscosity();

    int ist, ied, jst, jed, kst, ked;
    BCRegion->GetIJKRegion(ist, ied, jst, jed, kst, ked);

    int it, jt, kt;

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                for (int layer = 1; layer <= 2; ++ layer)
                {
					BCRegion->GetGhostCellIndex(i, j, k, it, jt, kt, layer);

                    for (int m = 0; m < nTurbulentEquation; ++ m)
                    {
                        qTurbulent(it, jt, kt, m) = freeStreamTurbVar[m];
                    }
                    viscousTurbulent(it, jt, kt) = freeStreamViscosity;
                }
            }
        }
    }
}

void TurbSolverStrFD::OutFlowBC(Grid *gridIn, BCRegion *BCRegion)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Param_TurbSolverStruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    RDouble4D &qTurbulent = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb"));
    RDouble3D &viscousTurbulent = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("vist"));

    int ist, ied, jst, jed, kst, ked;
    BCRegion->GetIJKRegion(ist, ied, jst, jed, kst, ked);

	int is1, js1, ks1, is2, js2, ks2;
	int it1, jt1, kt1, it2, jt2, kt2;

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
				BCRegion->GetInsideCellIndex(i, j, k, is1, js1, ks1, 1);
				BCRegion->GetInsideCellIndex(i, j, k, is2, js2, ks2, 2);

                RestrictIndex(is1, js1, ks1, ni-1, nj-1, nk-1);
                RestrictIndex(is2, js2, ks2, ni-1, nj-1, nk-1);

				BCRegion->GetGhostCellIndex(i, j, k, it1, jt1, kt1, 1);
				BCRegion->GetGhostCellIndex(i, j, k, it2, jt2, kt2, 2);

                    for (int m = 0; m < nTurbulentEquation; ++ m)
                    {
                        qTurbulent(it1, jt1, kt1, m) = two * qTurbulent(is1, js1, ks1, m) - qTurbulent(is2, js2, ks2, m);
                        qTurbulent(it2, jt2, kt2, m) = two * qTurbulent(it1, jt1, kt1, m) - qTurbulent(is1, js1, ks1, m);
                    }
                    viscousTurbulent(it1, jt1, kt1) = viscousTurbulent(is1, js1, ks1);
                    viscousTurbulent(it2, jt2, kt2) = viscousTurbulent(it1, jt1, kt1);
            }
        }
    }
}

void TurbSolverStrFD::VisWall(Grid *gridIn, BCRegion *BCRegion)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Param_TurbSolverStruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    RDouble4D &qLaminar           = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));
    RDouble4D &qTurbulent         = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb"));
    RDouble3D &viscousLaminars    = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("visl"));
    RDouble3D &viscousTurbulent   = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("vist"));

    RDouble3D &wallDistant = *StructGridCast(grid)->GetWallDist();

    RDouble refReNumber = parameters->GetRefReNumber();

    int ist, ied, jst, jed, kst, ked;
    BCRegion->GetIJKRegion(ist, ied, jst, jed, kst, ked);

    const RDouble beta1 = 0.075;

    using namespace IDX;

	int is1, js1, ks1, is2, js2, ks2;
	int it1, jt1, kt1, it2, jt2, kt2;

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
				BCRegion->GetInsideCellIndex(i, j, k, is1, js1, ks1, 1);
				BCRegion->GetInsideCellIndex(i, j, k, is2, js2, ks2, 2);

                RestrictIndex(is1, js1, ks1, ni-1, nj-1, nk-1);
                RestrictIndex(is2, js2, ks2, ni-1, nj-1, nk-1);

				BCRegion->GetGhostCellIndex(i, j, k, it1, jt1, kt1, 1);
				BCRegion->GetGhostCellIndex(i, j, k, it2, jt2, kt2, 2);

                if (nTurbulentEquation == 1)
                {
                    qTurbulent(it1, jt1, kt1, IDX::ISA) = 0.0;
                    qTurbulent(it2, jt2, kt2, IDX::ISA) = 0.0;
                }
                else if (nTurbulentEquation == 2)
                {
                    RDouble rho   = qLaminar(i, j, k, IR);

                    RDouble dist  = wallDistant(i, j, k);
                    RDouble dist2 = dist * dist;
                    RDouble muLaminar = viscousLaminars(i, j, k);

                    RDouble kw_wall = 60.0 * muLaminar / (rho * beta1 * dist2 * refReNumber * refReNumber);

                    qTurbulent(it1, jt1, kt1, IDX::IKE) = - qTurbulent(is1, js1, ks1, IDX::IKE);
                    qTurbulent(it2, jt2, kt2, IDX::IKE) =   qTurbulent(it1, jt1, kt1, IDX::IKE);
          
                    qTurbulent(it1, jt1, kt1, IDX::IKW) = 2.0 * kw_wall - qTurbulent(is1, js1, ks1, IDX::IKW);
                    qTurbulent(it2, jt2, kt2, IDX::IKW) = qTurbulent(it1, jt1, kt1, IDX::IKW);
                }
                viscousTurbulent(it1, jt1, kt1) = -viscousTurbulent(is1, js1, ks1);
                viscousTurbulent(it2, jt2, kt2) = -viscousTurbulent(is2, js2, ks2);
            }
        }
    }
}

void TurbSolverStrFD::VisWallWithWallFunctionStandard(Grid* gridIn, BCRegion* BCRegion)
{
    StructGrid* grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Param_TurbSolverStruct* parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    RDouble4D& qLaminar = *reinterpret_cast<RDouble4D*> (grid->GetDataPtr("q"));
    RDouble4D& qTurbulent = *reinterpret_cast<RDouble4D*> (grid->GetDataPtr("q_turb"));
    RDouble3D& viscousLaminars = *reinterpret_cast<RDouble3D*> (grid->GetDataPtr("visl"));
    RDouble3D& viscousTurbulent = *reinterpret_cast<RDouble3D*> (grid->GetDataPtr("vist"));

    RDouble4D& faceNormalX = *reinterpret_cast<RDouble4D*> (grid->GetFaceNormalX());
    RDouble4D& faceNormalY = *reinterpret_cast<RDouble4D*> (grid->GetFaceNormalY());
    RDouble4D& faceNormalZ = *reinterpret_cast<RDouble4D*> (grid->GetFaceNormalZ());

    RDouble3D& wallDistant = *StructGridCast(grid)->GetWallDist();

    RDouble refReNumber = parameters->GetRefReNumber();

    int ist, ied, jst, jed, kst, ked;
    BCRegion->GetIJKRegion(ist, ied, jst, jed, kst, ked);

    using namespace IDX;

    int iSurface = BCRegion->GetFaceDirection() + 1;

    const RDouble E_ = 9.793;
    const RDouble beta1 = 0.075;
    const RDouble kappa1 = 0.4178;
    const RDouble yPluslimit = 11.225;

    int is1, js1, ks1, is2, js2, ks2;
    int it1, jt1, kt1, it2, jt2, kt2;

    for (int k = kst; k <= ked; ++k)
    {
        for (int j = jst; j <= jed; ++j)
        {
            for (int i = ist; i <= ied; ++i)
            {

                BCRegion->GetInsideCellIndex(i, j, k, is1, js1, ks1, 1);
                BCRegion->GetInsideCellIndex(i, j, k, is2, js2, ks2, 2);

                RestrictIndex(is1, js1, ks1, ni - 1, nj - 1, nk - 1);
                RestrictIndex(is2, js2, ks2, ni - 1, nj - 1, nk - 1);

                BCRegion->GetGhostCellIndex(i, j, k, it1, jt1, kt1, 1);
                BCRegion->GetGhostCellIndex(i, j, k, it2, jt2, kt2, 2);

                RDouble rhoInside = qLaminar(is1, js1, ks1, IDX::IR);
                RDouble rhoGhost  = qLaminar(it1, jt1, kt1, IDX::IR);
                RDouble densityWall = 0.5 * (rhoInside + rhoGhost);
                RDouble up = qLaminar(is1, js1, ks1, IDX::IU);
                RDouble vp = qLaminar(is1, js1, ks1, IDX::IV);
                RDouble wp = qLaminar(is1, js1, ks1, IDX::IW);
                RDouble laminarViscosity = viscousLaminars(is1, js1, ks1);
                RDouble turbulentViscosity = viscousTurbulent(is1, js1, ks1);

                RDouble xfn = faceNormalX(is1, js1, ks1, iSurface);
                RDouble yfn = faceNormalY(is1, js1, ks1, iSurface);
                RDouble zfn = faceNormalZ(is1, js1, ks1, iSurface);

                RDouble vpn = up * xfn + vp * yfn + wp * zfn;
                RDouble vpt = sqrt(up * up + vp * vp + wp * wp - vpn * vpn);
                RDouble small = 1.0e-15;
                RDouble uvw = max(vpt, small);

                RDouble wallDistance = wallDistant(is1, js1, ks1);

                //! Obtain the reynold number.
                RDouble refReNumber = parameters->GetRefReNumber();

                RDouble taow = laminarViscosity * uvw / wallDistance / refReNumber;
                RDouble utao = sqrt(taow / densityWall);

                RDouble yPlus;
                RDouble epslon = 1.0;
                int iter = 0;
                while (epslon > 0.001 && ++iter < 20)
                {
                    yPlus = refReNumber * wallDistance * densityWall * utao / laminarViscosity;
                    yPlus = MAX(yPlus, yPluslimit);
                    RDouble utemp = utao;

                    utao = kappa1 * uvw / (log(E_ * yPlus));

                    epslon = abs(utao - utemp) / (abs(utemp) + 1.0e-20);
                }
                taow = densityWall * utao * utao;
                RDouble viscousTurbulenceWall = laminarViscosity * (taow / (uvw + small) / laminarViscosity * wallDistance * refReNumber - 1.0);
                viscousTurbulent(it1, jt1, kt1) = 2.0 * viscousTurbulenceWall - viscousTurbulent(is1, js1, ks1);
                viscousTurbulent(it2, jt2, kt2) = viscousTurbulent(it1, jt1, kt1);

                if (nTurbulentEquation == 1)
                {
                    qTurbulent(it1, jt1, kt1, IDX::ISA) = -qTurbulent(is1, js1, ks1, IDX::ISA);
                    qTurbulent(it2, jt2, kt2, IDX::ISA) =  qTurbulent(it1, jt1, kt1, IDX::ISA);
                }
                else if (nTurbulentEquation == 2)
                {
                    RDouble distanceSquare = wallDistance * wallDistance;

                    RDouble omgi = 60.0 * laminarViscosity / (rhoInside * beta1 * distanceSquare * refReNumber * refReNumber);
                    //RDouble omgo = utao / (0.126 * wallDistance * refReNumber);
                    RDouble omgo = 0.0;
                    RDouble kwWall = sqrt(omgi * omgi + omgo * omgo);

                    RDouble qTurbulentFirstCell = kwWall * viscousTurbulenceWall / densityWall;
                    //qTurbulent(is1, js1, ks1, IDX::IKE) = qTurbulentFirstCell;
                    qTurbulent(it1, jt1, kt1, IDX::IKE) = -qTurbulent(is1, js1, ks1, IDX::IKE);
                    qTurbulent(it2, jt2, kt2, IDX::IKE) = qTurbulent(it1, jt1, kt1, IDX::IKE);

                    //qTurbulent(is1, js1, ks1, IDX::IKW) = kwWall;
                    qTurbulent(it1, jt1, kt1, IDX::IKW) = 2.0 * kwWall - qTurbulent(is1, js1, ks1, IDX::IKW);
                    //qTurbulent(it1, jt1, kt1, IDX::IKW) = kwWall;
                    //qTurbulent(it1, jt1, kt1, IDX::IKW) = -qTurbulent(it1, jt1, kt1, IDX::IKE) * rhoGhost / viscousTurbulent(it1, jt1, kt1);
                    qTurbulent(it2, jt2, kt2, IDX::IKW) = qTurbulent(it1, jt1, kt1, IDX::IKW);
                    //qTurbulent(it1, jt1, kt1, IDX::IKE) = qTurbulent(it1, jt1, kt1, IDX::IKW) * viscousTurbulent(it1, jt1, kt1) / rhoGhost;
                    //qTurbulent(it2, jt2, kt2, IDX::IKE) = qTurbulent(it2, jt2, kt2, IDX::IKW) * viscousTurbulent(it2, jt2, kt2) / rhoGhost;
                }
            }
        }
    }
}

void TurbSolverStrFD::VisWallWithWallFunctionPAB3D(Grid* gridIn, BCRegion* BCRegion)
{
    StructGrid* grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    using namespace IDX;

    RDouble paba10[7] = { 2.354039, 0.1179840, -4.2899192e-04, 2.0404148e-06,-5.1775775e-09, 6.2687308e-12, -2.916958e-15 };
    RDouble paba11[5] = { 5.777191, 6.8756983e-02, -7.1582745e-06, 1.5594904e-09, -1.4865778e-13 };
    RDouble paba12[5] = { 31.08654, 5.0429072e-02, -2.0072314e-8 };
    RDouble Ylim = 10.0;
    RDouble beta1 = 0.075;

    Param_TurbSolverStruct* parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    RDouble4D& qLaminar = *reinterpret_cast<RDouble4D*> (grid->GetDataPtr("q"));
    RDouble4D& qTurbulent = *reinterpret_cast<RDouble4D*> (grid->GetDataPtr("q_turb"));
    RDouble4D& temperature = *reinterpret_cast<RDouble4D*> (grid->GetDataPtr("t"));
    RDouble3D& viscousLaminars = *reinterpret_cast<RDouble3D*> (grid->GetDataPtr("visl"));
    RDouble3D& viscousTurbulent = *reinterpret_cast<RDouble3D*> (grid->GetDataPtr("vist"));

    RDouble4D& faceNormalX = *reinterpret_cast<RDouble4D*> (grid->GetFaceNormalX());
    RDouble4D& faceNormalY = *reinterpret_cast<RDouble4D*> (grid->GetFaceNormalY());
    RDouble4D& faceNormalZ = *reinterpret_cast<RDouble4D*> (grid->GetFaceNormalZ());

    RDouble3D& wallDistant = *StructGridCast(grid)->GetWallDist();

    RDouble refReNumber = parameters->GetRefReNumber();
    RDouble refMachNumber = parameters->GetRefMachNumber();
    RDouble refGama = parameters->GetRefGama();

    int ist, ied, jst, jed, kst, ked;
    BCRegion->GetIJKRegion(ist, ied, jst, jed, kst, ked);

    using namespace IDX;

    int iSurface = BCRegion->GetFaceDirection() + 1;

    for (int k = kst; k <= ked; ++k)
    {
        for (int j = jst; j <= jed; ++j)
        {
            for (int i = ist; i <= ied; ++i)
            {
                int is1, js1, ks1, is2, js2, ks2;
                int it1, jt1, kt1, it2, jt2, kt2;

                BCRegion->GetInsideCellIndex(i, j, k, is1, js1, ks1, 1);
                BCRegion->GetInsideCellIndex(i, j, k, is2, js2, ks2, 2);

                RestrictIndex(is1, js1, ks1, ni - 1, nj - 1, nk - 1);
                RestrictIndex(is2, js2, ks2, ni - 1, nj - 1, nk - 1);

                BCRegion->GetGhostCellIndex(i, j, k, it1, jt1, kt1, 1);
                BCRegion->GetGhostCellIndex(i, j, k, it2, jt2, kt2, 2);

                RDouble rhoInside = qLaminar(is1, js1, ks1, IDX::IR);
                RDouble rhoGhost = qLaminar(it1, jt1, kt1, IDX::IR);
                RDouble up = qLaminar(is1, js1, ks1, IDX::IU);
                RDouble vp = qLaminar(is1, js1, ks1, IDX::IV);
                RDouble wp = qLaminar(is1, js1, ks1, IDX::IW);
                RDouble laminarViscosity = viscousLaminars(is1, js1, ks1);
                RDouble turbulentViscosity = viscousTurbulent(is1, js1, ks1);

                RDouble xfn = faceNormalX(is1, js1, ks1, iSurface);
                RDouble yfn = faceNormalY(is1, js1, ks1, iSurface);
                RDouble zfn = faceNormalZ(is1, js1, ks1, iSurface);

                RDouble vpn = up * xfn + vp * yfn + wp * zfn;
                RDouble vpt = sqrt(up * up + vp * vp + wp * wp - vpn * vpn);
                RDouble small = 1.0e-15;
                RDouble uvw = max(vpt, small);

                //RDouble uu = sqrt(up * up + vp * vp + wp * wp);
                RDouble& uu = uvw;

                RDouble wallDistance = wallDistant(is1, js1, ks1);

                RDouble refReNumber = parameters->GetRefReNumber(); //! Obtain the reynold number.

                RDouble rc = rhoInside * refReNumber * wallDistance * uu / laminarViscosity;
                RDouble xnplus = 0.0;
                //! The original one here is yPlus.
                if (rc <= 20.24)
                {
                    xnplus = sqrt(rc);
                }
                else if (rc <= 435.0)
                {
                    xnplus = paba10[0] + paba10[1] * rc + paba10[2] * rc * rc + paba10[3] * rc * rc * rc + paba10[4] * rc * rc * rc * rc + paba10[5] * rc * rc * rc * rc * rc + paba10[6] * rc * rc * rc * rc * rc * rc;
                }
                else if (rc <= 4000.0)
                {
                    xnplus = paba11[0] + paba11[1] * rc + paba11[2] * rc * rc + paba11[3] * rc * rc * rc + paba11[4] * rc * rc * rc * rc;
                }
                else if (rc > 4000.0)
                {
                    xnplus = paba12[0] + paba12[1] * rc + paba12[2] * rc * rc;
                }
                RDouble xnplussav = xnplus;

                //if (xnplus >= Ylim)
                {
                    //! Newton iteration to solve for nplus, assuming it is in log region:
                    for (int iter = 0; ; ++iter)
                    {
                        RDouble f = rc / xnplus - 2.44 * (log(xnplus)) - 5.2;
                        RDouble dfdn = -rc / (xnplus * xnplus) - 2.44 / xnplus;
                        RDouble delta = -f / dfdn;
                        xnplus = fabs(xnplus + delta);
                        if (iter > 10)
                        {
                            //! Revert back to approx series soln if Newton iteration fails.
                            xnplus = xnplussav;
                            break;
                        }
                        else if (abs(delta) < 1.e-3)
                        {
                            break;
                        }
                    }
                }

                RDouble viscousTurbulenceWall = laminarViscosity * (laminarViscosity * xnplus * xnplus / (wallDistance * rhoInside * uu * refReNumber) - 1.0);
                viscousTurbulenceWall = MAX(viscousTurbulenceWall, 0.0);
                viscousTurbulent(it1, jt1, kt1) = 2.0 * viscousTurbulenceWall - viscousTurbulent(is1, js1, ks1);
                viscousTurbulent(it2, jt2, kt2) = viscousTurbulent(it1, jt1, kt1);

                if (nTurbulentEquation == 1)
                {
                    qTurbulent(it1, jt1, kt1, IDX::ISA) = -qTurbulent(is1, js1, ks1, IDX::ISA);
                    qTurbulent(it2, jt2, kt2, IDX::ISA) =  qTurbulent(it1, jt1, kt1, IDX::ISA);
                }
                else if (nTurbulentEquation == 2)
                {
                    RDouble distanceSquare = wallDistance * wallDistance;

                    RDouble omgi = 60.0 * laminarViscosity / (rhoInside * beta1 * distanceSquare * refReNumber * refReNumber);
                    RDouble taow = laminarViscosity * uvw / wallDistance / refReNumber;
                    RDouble utao = sqrt(taow / rhoInside);
                    RDouble omgo = utao / (0.126 * wallDistance * refReNumber);
                    //RDouble omgo = 0.0;
                    RDouble kwWall = sqrt(omgi * omgi + omgo * omgo);

                    //qTurbulent(is1, js1, ks1, IDX::IKW) = kwWall;
                    qTurbulent(it1, jt1, kt1, IDX::IKW) = 2.0 * kwWall - qTurbulent(is1, js1, ks1, IDX::IKW);
                    qTurbulent(it2, jt2, kt2, IDX::IKW) = qTurbulent(it1, jt1, kt1, IDX::IKW);

                    qTurbulent(it1, jt1, kt1, IDX::IKE) = qTurbulent(it1, jt1, kt1, IDX::IKW) * viscousTurbulent(it1, jt1, kt1) / rhoGhost;
                    //qTurbulent(is1, js1, ks1, IDX::IKE) = qTurbulent(is1, js1, ks1, IDX::IKW) * viscousTurbulent(is1, js1, ks1) / rhoInside;
                    qTurbulent(it2, jt2, kt2, IDX::IKE) = qTurbulent(it1, jt1, kt1, IDX::IKE);
                }
            }
        }
    }
}

void TurbSolverStrFD::SymmetryBC(Grid *gridIn, BCRegion *BCRegion)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Param_TurbSolverStruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    RDouble4D &qTurbulent       = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb"));
    RDouble3D &viscousTurbulent = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("vist"));

    int ist, ied, jst, jed, kst, ked;
    BCRegion->GetIJKRegion(ist, ied, jst, jed, kst, ked);

	int is, js, ks, it, jt, kt;

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                for (int layer = 1; layer <= 2; ++ layer)
                {
					BCRegion->GetInsideCellIndex(i, j, k, is, js, ks, layer);
                    RestrictIndex(is, js, ks, ni-1, nj-1, nk-1);
					BCRegion->GetGhostCellIndex(i, j, k, it, jt, kt, layer);

                    for (int m = 0; m < nTurbulentEquation; ++ m)
                    {
                        qTurbulent(it, jt, kt, m) = qTurbulent(is, js, ks, m);
                    }
                    viscousTurbulent(it, jt, kt) = viscousTurbulent(is, js, ks);
                }
            }
        }
    }
}

void TurbSolverStrFD::FarFieldBC(Grid *gridIn, BCRegion *BCRegion)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    RDouble4D &faceNormalX  = *(grid->GetFaceNormalX());
    RDouble4D &faceNormalY  = *(grid->GetFaceNormalY());
    RDouble4D &faceNormalZ  = *(grid->GetFaceNormalZ());

    Param_TurbSolverStruct *parameters = GetControlParameters();
    RDouble *freeStreamTurbVar = parameters->GetFreeStreamTurbVar();

    RDouble4D &qTurbulent       = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb"));
    RDouble4D &qLaminar         = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));
    RDouble3D &gamma            = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("gama"));
    RDouble3D &viscousTurbulent = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("vist"));

    int nm = GlobalDataBase::GetIntParaFromDB("nm");
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    RDouble refGama = GlobalDataBase::GetDoubleParaFromDB("refGama");
    RDouble referenceMachNumber = parameters->GetRefMachNumber();
    RDouble freeStreamViscosity = parameters->GetFreeStreamViscosity();

    using namespace IDX;

    RDouble *primitiveVariablesInfinity = reinterpret_cast<RDouble *> (GlobalDataBase::GetDataPtr("prim_inf" ));
    RDouble rFarfield = primitiveVariablesInfinity[IR];
    RDouble uFarfield = primitiveVariablesInfinity[IU];
    RDouble vFarfield = primitiveVariablesInfinity[IV];
    RDouble wFarfield = primitiveVariablesInfinity[IW];
    RDouble pFarfield = primitiveVariablesInfinity[IP];

    int ist, ied, jst, jed, kst, ked;
    BCRegion->GetIJKRegion(ist, ied, jst, jed, kst, ked);

    int iSurface = BCRegion->GetFaceDirection() + 1;
    int leftOrRightIndex    = BCRegion->GetFaceLeftOrRightIndex();

    RDouble *prims1 = new RDouble [nm];

    int is1, js1, ks1, is2, js2, ks2;
	int it1, jt1, kt1, it2, jt2, kt2;

    using namespace IDX;

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
				int in, jn, kn;
				BCRegion->GetBoundaryFaceIndex(i, j, k, in, jn, kn);
                RDouble nxs = leftOrRightIndex * faceNormalX(in, jn, kn, iSurface);
                RDouble nys = leftOrRightIndex * faceNormalY(in, jn, kn, iSurface);
                RDouble nzs = leftOrRightIndex * faceNormalZ(in, jn, kn, iSurface);

				BCRegion->GetInsideCellIndex(i, j, k, is1, js1, ks1, 1);
				BCRegion->GetInsideCellIndex(i, j, k, is2, js2, ks2, 2);

                RestrictIndex(is1, js1, ks1, ni-1, nj-1, nk-1);
                RestrictIndex(is2, js2, ks2, ni-1, nj-1, nk-1);

				BCRegion->GetGhostCellIndex(i, j, k, it1, jt1, kt1, 1);
				BCRegion->GetGhostCellIndex(i, j, k, it2, jt2, kt2, 2);

                for (int m = 0; m < nm; ++ m)
                {
                    prims1[m] = qLaminar(is1, js1, ks1, m);
                }

                RDouble rin = prims1[IR];
                RDouble uin = prims1[IU];
                RDouble vin = prims1[IV];
                RDouble win = prims1[IW];
                RDouble pin = prims1[IP];

                RDouble vno = nxs * uFarfield + nys * vFarfield + nzs * wFarfield;
                RDouble vni = nxs * uin + nys * vin + nzs * win;
                RDouble vei = sqrt(uin * uin + vin * vin + win * win);

                RDouble gama = gamma(is1, js1, ks1);

                RDouble coo = sqrt(ABS(refGama * pFarfield / rFarfield));
                RDouble cin = sqrt(ABS(gama  * pin / rin));

                if (vei > cin)
                {
                    if (vni >= 0.0)
                    {
                        for (int m = 0; m < nTurbulentEquation; ++ m)
                        {
                            qTurbulent(it1, jt1, kt1, m) = qTurbulent(is1, js1, ks1, m);
                            qTurbulent(it2, jt2, kt2, m) = qTurbulent(is2, js2, ks2, m);
                        }
                        viscousTurbulent(it1, jt1, kt1) = viscousTurbulent(is1, js1, ks1);
                        viscousTurbulent(it2, jt2, kt2) = viscousTurbulent(is2, js2, ks2);
                    }
                    else
                    {
                        for (int m = 0; m < nTurbulentEquation; ++ m)
                        {
                            qTurbulent(it1, jt1, kt1,m) = freeStreamTurbVar[m];
                            qTurbulent(it2, jt2, kt2,m) = freeStreamTurbVar[m];
                        }
                        viscousTurbulent(it1, jt1, kt1) = freeStreamViscosity;
                        viscousTurbulent(it2, jt2, kt2) = freeStreamViscosity;
                    }
                    continue;
                }

                RDouble riemp = vni + 2.0 * cin / (gama  - 1.0);
                RDouble riemm = vno - 2.0 * coo / (refGama - 1.0);
                RDouble vnb = half * (riemp + riemm);
                RDouble cb  = fourth * (riemp - riemm) * (gama - 1.0);

                if (vnb  >= 0.0)
                {
                    //! Exit.
                    for (int m = 0; m < nTurbulentEquation; ++ m)
                    {
                        qTurbulent(it1, jt1, kt1, m) = qTurbulent(is1, js1, ks1, m);
                        qTurbulent(it2, jt2, kt2, m) = qTurbulent(is2, js2, ks2, m);
                    }
                    viscousTurbulent(it1, jt1, kt1) = viscousTurbulent(is1, js1, ks1);
                    viscousTurbulent(it2, jt2, kt2) = viscousTurbulent(it1, jt1, kt1);
                }
                else
                {
                    //! Inlet.
                    for (int m = 0; m < nTurbulentEquation; ++ m)
                    {
                        qTurbulent(it1, jt1, kt1, m) = freeStreamTurbVar[m];
                        qTurbulent(it2, jt2, kt2, m) = freeStreamTurbVar[m];
                    }
                    viscousTurbulent(it1, jt1, kt1) = freeStreamViscosity;
                    viscousTurbulent(it2, jt2, kt2) = freeStreamViscosity;
                }
            }
        }
    }
    delete []prims1;
}

void TurbSolverStrFD::UpdateFlowField(Grid *gridIn, FieldProxy *qProxy, FieldProxy *dqProxy)
{
    StructGrid *grid = StructGridCast(gridIn);

    RDouble4D &qTurbulent = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb"));
    RDouble4D &dqTurbulent = dqProxy->GetField_STR();

    Param_TurbSolverStruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    RDouble4D &qpmv = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q" ));

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    Range M(0, nTurbulentEquation - 1);

    using namespace IDX;
    if (nTurbulentEquation == 1)
    {
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    qTurbulent(i, j, k,ISA) = MAX(qTurbulent(i, j, k, ISA) + dqTurbulent(i, j, k, ISA), zero);
                }
            }
        }
    }
    else if (nTurbulentEquation == 2)
    {
        RDouble   *freeStreamTurbVar =  parameters->GetFreeStreamTurbVar();
        RDouble kelim = 1.0e-5 * freeStreamTurbVar[IKE];
        RDouble kwlim = 1.0e-5 * freeStreamTurbVar[IKW];

        RDouble ke,kw,dke,dkw,rho;

        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {                  
                    rho = qpmv(i, j, k,IR);
                    ke  = qTurbulent(i, j, k, IKE);
                    kw  = qTurbulent(i, j, k, IKW); 
                    dke = dqTurbulent(i, j, k, IKE);
                    dkw = dqTurbulent(i, j, k, IKW);
                    ke += dke / rho;
                    kw += dkw / rho;
                    ke = MAX(ke, kelim);
                    kw = MAX(kw, kwlim);
                    qTurbulent(i, j, k, IKE) = ke;
                    qTurbulent(i, j, k, IKW) = kw;
                }
            }
        }
    }
    ObtainViscosity(gridIn);
    ObtainBoundaryValue(gridIn);
}

void TurbSolverStrFD::Diagonal(Grid *gridIn)
{
    Param_TurbSolverStruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    if (nTurbulentEquation == 1)
    {
        SpectrumRadiusOfOneEquation(gridIn);
    }
    else if (nTurbulentEquation == 2)
    {
        SpectrumRadiusOfTwoEquation(gridIn);
    }
}

void TurbSolverStrFD::SpectrumRadiusOfOneEquation(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    RDouble4D &faceVectorX  = *(grid->GetFaceVectorX());
    RDouble4D &faceVectorY  = *(grid->GetFaceVectorY());
    RDouble4D &faceVectorZ  = *(grid->GetFaceVectorZ());
    RDouble4D &area         = *(grid->GetFaceArea());
    RDouble3D &volume       = *(grid->GetCellVolume());
    
    RDouble4D &qLaminar         = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"        ));
    RDouble3D &viscousLaminar   = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("visl"     ));
    RDouble4D &qTurbulent       = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb"   ));
    
    RDouble5D &matrixTurbLeft   = *reinterpret_cast<RDouble5D *> (grid->GetDataPtr("mat_turbl"));
    RDouble5D &matrixTurbRight  = *reinterpret_cast<RDouble5D *> (grid->GetDataPtr("mat_turbr"));
    Param_TurbSolverStruct *parameters = GetControlParameters();

    int nTurbulentEquation = parameters->GetNTurbulentEquation();
    RDouble turbCFLScale = parameters->GetTurbCFLScale();

    RDouble oRefReNumber = parameters->GetoRefReNumber();

    //! Initialize some constants.
    RDouble sigma  = 2.0 / 3.0;
    RDouble osigma = one / sigma;

    matrixTurbLeft  = 0.0;
    matrixTurbRight = 0.0;

    int nDim = GetDim();

    using namespace IDX;

    for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
    {
        int il1, jl1, kl1;
        GetNsurfIndex(iSurface, il1, jl1, kl1);
        
        int ist = 1 - il1;
        int jst = 1 - jl1;
        int kst = 1 - kl1;
        int ied = ni - 1 + il1;
        int jed = nj - 1 + jl1;
        int ked = nk - 1 + kl1;
        if (nk == 1)
        {
            kst = 1;
            ked = 1;
        }
        
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    int il = i + il1;
                    int jl = j + jl1;
                    int kl = k + kl1;

                    RDouble nx0 = faceVectorX(i, j, k, iSurface);
                    RDouble ny0 = faceVectorY(i, j, k, iSurface);
                    RDouble nz0 = faceVectorZ(i, j, k, iSurface);

                    RDouble nx1 = faceVectorX(il, jl, kl, iSurface);
                    RDouble ny1 = faceVectorY(il, jl, kl, iSurface);
                    RDouble nz1 = faceVectorZ(il, jl, kl, iSurface);

                    RDouble nx = 0.5 * (nx0 + nx1);
                    RDouble ny = 0.5 * (ny0 + ny1);
                    RDouble nz = 0.5 * (nz0 + nz1);
                    RDouble ns  = nx * nx +  ny * ny +  nz * nz;

                    RDouble ul  = qLaminar(i, j, k, IU);
                    RDouble vl  = qLaminar(i, j, k, IV);
                    RDouble wl  = qLaminar(i, j, k, IW);
                    RDouble Vn = ul * nx + vl * ny + wl * nz;
                    RDouble absVn = ABS(Vn);

                    RDouble rl = qLaminar(i, j, k, IR);
                    RDouble nul = qTurbulent(i, j, k, ISA);
                    RDouble nueff_l = ABS(viscousLaminar(i , j , k ) / rl + nul);

                    RDouble cellVolume  = volume(i, j, k);
                    RDouble oVolume = one / cellVolume;

                    matrixTurbLeft (i, j, k, ISA, iSurface) =   half * (Vn + absVn) + oRefReNumber * oVolume * (osigma * nueff_l * ns);
                    matrixTurbRight(i, j, k, ISA, iSurface) = - half * (Vn - absVn) + oRefReNumber * oVolume * (osigma * nueff_l * ns);
                }
            }
        }
    }

    RDouble4D &spectumTurbulent = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("spec_turb"));
    RDouble3D &dt               = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("dt"       ));

    RDouble dualTimeSpectrumC1 = zero;
    RDouble dualTimeSpectrumC2 = one;

    //! If flow is unsteady, need to count the contribution of unsteady.
    int isUnsteady = parameters->GetIsUnsteady();
    if (isUnsteady)
    {
        RDouble dualTimeCoefficient[7];
        const int methodOfDualTime = 3;
        ComputeDualTimeCoefficient(methodOfDualTime, dualTimeCoefficient);
        dualTimeSpectrumC1 = - dualTimeCoefficient[3];
        dualTimeSpectrumC2 =   dualTimeCoefficient[6];
    }

    int ist = 1;
    int ied = ni - 1;
    int jst = 1;
    int jed = nj - 1;
    int kst = 1;
    int ked = nk - 1;
    if (nk == 1)
    {
        kst = 1;
        ked = 1;
    }

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble &specTub = spectumTurbulent(i, j, k, ISA);
                specTub += dualTimeSpectrumC1 * volume(i, j, k) + dualTimeSpectrumC2 / (turbCFLScale * dt(i, j, k) + SMALL);
                
                for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
                {
                    specTub += matrixTurbLeft (i, j, k, ISA, iSurface) + matrixTurbRight(i, j, k, ISA, iSurface);
                }
            }
        }
    }
}

void TurbSolverStrFD::SpectrumRadiusOfTwoEquation(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    RDouble4D &faceVectorX  = *(grid->GetFaceVectorX());
    RDouble4D &faceVectorY  = *(grid->GetFaceVectorY());
    RDouble4D &faceVectorZ  = *(grid->GetFaceVectorZ());
    RDouble4D &area   = *(grid->GetFaceArea()  );
    RDouble3D &volume = *(grid->GetCellVolume());
    
    RDouble4D &qLaminar         = *reinterpret_cast< RDouble4D *> (grid->GetDataPtr("q"        ));
    RDouble3D &viscousLaminar   = *reinterpret_cast< RDouble3D *> (grid->GetDataPtr("visl"     ));
    RDouble3D &viscousTurbulent = *reinterpret_cast< RDouble3D *> (grid->GetDataPtr("vist"     ));
    RDouble3D &blend            = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("blend"));

    RDouble5D &matrixTurbLeft  = *reinterpret_cast< RDouble5D *> (grid->GetDataPtr("mat_turbl" ));
    RDouble5D &matrixTurbRight = *reinterpret_cast< RDouble5D *> (grid->GetDataPtr("mat_turbr" ));

    matrixTurbLeft  = 0.0;
    matrixTurbRight = 0.0;

    Param_TurbSolverStruct *parameters = GetControlParameters();
    RDouble turbCFLScale = parameters->GetTurbCFLScale();
    RDouble oRefReNumber = parameters->GetoRefReNumber();
    RDouble KW_sigmaK1   = parameters->GetKW_sigmaK1();
    RDouble KW_sigmaK2   = parameters->GetKW_sigmaK2();
    RDouble KW_sigmaW1   = parameters->GetKW_sigmaW1();
    RDouble KW_sigmaW2   = parameters->GetKW_sigmaW2();
    
    using namespace IDX;

    int mst = 0;
    int med = 1;

    RDouble work[2];

    int nDim = GetDim();
    for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
    {
        int il1, jl1, kl1;
        GetNsurfIndex(iSurface, il1, jl1, kl1);

        int ist = 1 - il1;
        int jst = 1 - jl1;
        int kst = 1 - kl1;
        int ied = ni - 1 + il1;
        int jed = nj - 1 + jl1;
        int ked = nk - 1 + kl1;
        if (nk == 1)
        {
            kst = 1;
            ked = 1;
        }

        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    int il = i + il1;
                    int jl = j + jl1;
                    int kl = k + kl1;

                    RDouble nx0 = faceVectorX(i, j, k, iSurface);
                    RDouble ny0 = faceVectorY(i, j, k, iSurface);
                    RDouble nz0 = faceVectorZ(i, j, k, iSurface);

                    RDouble nx1 = faceVectorX(il, jl, kl, iSurface);
                    RDouble ny1 = faceVectorY(il, jl, kl, iSurface);
                    RDouble nz1 = faceVectorZ(il, jl, kl, iSurface);

                    RDouble nx = 0.5 * (nx0 + nx1);
                    RDouble ny = 0.5 * (ny0 + ny1);
                    RDouble nz = 0.5 * (nz0 + nz1);
                    RDouble ns  = nx * nx + ny * ny + nz * nz;

                    RDouble ul  = qLaminar(i, j, k, IU);
                    RDouble vl  = qLaminar(i, j, k, IV);
                    RDouble wl  = qLaminar(i, j, k, IW);
                    RDouble Vn = ul * nx + vl * ny + wl * nz;
                    RDouble absVn = ABS(Vn);
                    
                    RDouble rhoAverage = qLaminar(i, j, k,IR);
                    RDouble muLaminarAverage = viscousLaminar(i ,j ,k);
                    RDouble muTurbulentAverage = viscousTurbulent(i ,j ,k);
                    
                    RDouble cblendAverage = blend(i, j, k);
                    RDouble tsigKTurbulent = cblendAverage * KW_sigmaK1 + (1.0 - cblendAverage) * KW_sigmaK2;
                    RDouble tsigWTurbulent = cblendAverage * KW_sigmaW1 + (1.0 - cblendAverage) * KW_sigmaW2;
                    
                    RDouble cellVolumeAverage  = volume(i, j, k);
                    RDouble oVolume = one / cellVolumeAverage;
                    
                    work[0] = (muLaminarAverage + muTurbulentAverage * tsigKTurbulent) / rhoAverage * ns * oVolume * oRefReNumber;
                    work[1] = (muLaminarAverage + muTurbulentAverage * tsigWTurbulent) / rhoAverage * ns * oVolume * oRefReNumber;
                    
                    for (int m = mst; m <= med; ++ m)
                    {
                        matrixTurbLeft (i, j, k, m, iSurface) =   half * (Vn + absVn) + work[m];
                        matrixTurbRight(i, j, k, m, iSurface) = - half * (Vn - absVn) + work[m];
                    }
                }
            }
        }
    }

    RDouble4D &spectumTurbulent = *reinterpret_cast< RDouble4D *> (grid->GetDataPtr("spec_turb"));
    RDouble3D &dt               = *reinterpret_cast< RDouble3D *> (grid->GetDataPtr("dt"       ));

    RDouble dualTimeSpectrumC1 = zero;
    RDouble dualTimeSpectrumC2 = one;

    //! If flow is unsteady, it needs to count the contribution of unsteady.
    int isUnsteady = parameters->GetIsUnsteady();
    if (isUnsteady)
    {
        RDouble dualTimeCoefficient[7];
        const int methodOfDualTime = 3;
        ComputeDualTimeCoefficient(methodOfDualTime, dualTimeCoefficient);
        dualTimeSpectrumC1 = - dualTimeCoefficient[3];
        dualTimeSpectrumC2 =   dualTimeCoefficient[6];
    }

    int ist = 1;
    int ied = ni - 1;
    int jst = 1;
    int jed = nj - 1;
    int kst = 1;
    int ked = nk - 1;
    if (nk == 1)
    {
        kst = 1;
        ked = 1;
    }

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble volumeCell = volume(i, j, k);
                RDouble dtCell = dt(i, j, k);
                RDouble specTurbDelta = dualTimeSpectrumC1 * volumeCell + dualTimeSpectrumC2 / (turbCFLScale * dtCell + SMALL);
                for (int m = mst; m <= med; ++ m)
                {
                    RDouble &specTur = spectumTurbulent(i, j, k, m);
                    specTur += specTurbDelta;

                    for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
                    {
                        specTur += matrixTurbLeft (i, j, k, m, iSurface) + matrixTurbRight(i, j, k, m, iSurface);
                    }
                }
            }
        }
    }
}

void TurbSolverStrFD::LUSGSInitializationStructHighOrder(Grid *gridIn, FieldProxy *dqProxy, FieldProxy *LUplusDQ)
{
    StructGrid *grid = StructGridCast(gridIn);   
    RDouble4D &spectumTurbulent = *reinterpret_cast< RDouble4D *> (grid->GetDataPtr("spec_turb"));
    RDouble4D &residualTurbulent = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_turb"));
    RDouble4D &dq = dqProxy->GetField_STR();
    RDouble4D &dRHS = LUplusDQ->GetField_STR();
    
    dRHS = 0.0;
    dq = 0.0;
    //! This statement is necessary.
    //! Otherwise debug failed to run for unknown reason.

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    int ist = 1;
    int ied = ni - 1;
    int jst = 1;
    int jed = nj - 1;
    int kst = 1;
    int ked = nk - 1;
    if (nk == 1)
    {
        kst = 1;
        ked = 1;
    }

    Param_TurbSolverStruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                for (int m = 0; m < nTurbulentEquation; ++ m)
                {
                    dq(i, j, k, m) = residualTurbulent(i, j, k, m) / spectumTurbulent(i, j, k, m);
                }
            }
        }
    }     
}

void TurbSolverStrFD::CompressDQ(DataContainer *&dataContainer, FieldProxy *fieldProxy, Grid *gridIn, const int &zoneIndex, const int &neighborZoneIndex, const int &nEquation)
{
    StructGrid *grid = StructGridCast(gridIn);
    int level = grid->GetLevel();
    InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
    if (!interfaceInformation)
    {
        return;
    }

    StructGrid *finestGrid = StructGridCast(grid->GetFinestGrid());
    InterfaceInfo *finestInterfaceInfo = finestGrid->GetInterfaceInfo();

    int iNeighborZone                            = finestInterfaceInfo->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = finestInterfaceInfo->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForSend          = finestInterfaceInfo->GetFaceIndexForSend(iNeighborZone);

    dataContainer->moveto_begin();

    RDouble4D &fieldSend = fieldProxy->GetField_STR();

    for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; -- iGhostLayer)
    {
        for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
        {
            int iFace = interfaceIndexContainerForSend[iLocalFace];
            int is, js, ks;
            finestGrid->GetSourceIndexIJK(iFace, iGhostLayer + 1, is, js, ks);
            finestGrid->RemapMultigridIJK(level, is, js, ks);
            for (int m = 0; m < nEquation; ++ m)
            {
                PHWrite(dataContainer, fieldSend(is, js, ks, m));
            }
        }
    }
}

void TurbSolverStrFD::DecompressDQ(DataContainer *&dataContainer, FieldProxy *fieldProxy, Grid *gridIn, const int &zoneIndex, const int &neighborZoneIndex, const int &nEquation)
{
    StructGrid *grid = StructGridCast(gridIn);
    int level = grid->GetLevel();
    InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
    if (!interfaceInformation)
    {
        return;
    }

    StructGrid *finestGrid = StructGridCast(grid->GetFinestGrid());
    InterfaceInfo *finestInterfaceInfo = finestGrid->GetInterfaceInfo();

    int iNeighborZone                            = finestInterfaceInfo->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = finestInterfaceInfo->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForReceive       = finestInterfaceInfo->GetFaceIndexForRecv(iNeighborZone);

    dataContainer->moveto_begin();

    RDouble4D &fieldRecv = fieldProxy->GetField_STR();

    for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; -- iGhostLayer)
    {
        for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
        {
            int iFace = interfaceIndexContainerForReceive[iLocalFace];
            int it, jt, kt;
            finestGrid->GetTargetIndexIJK(iFace, iGhostLayer + 1, it, jt, kt);
            finestGrid->RemapMultigridIJK(level, it, jt, kt);
            for (int m = 0; m < nEquation; ++ m)
            {
                PHRead(dataContainer, fieldRecv(it, jt, kt, m));
            }
        }
    }
}

void TurbSolverStrFD::SolveLUSGSForward(Grid *gridIn, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal)
{
    StructGrid *grid = StructGridCast(gridIn);

    RDouble4D &dqTurbulent = dqProxy->GetField_STR();
    RDouble4D &dRHS        = LUplusDQ->GetField_STR();

    Param_TurbSolverStruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    RDouble *dqNeighbor = new RDouble[nTurbulentEquation];
    RDouble *df_total   = new RDouble[nTurbulentEquation];
    RDouble *df         = new RDouble[nTurbulentEquation];
    RDouble *mat_abc    = new RDouble[nTurbulentEquation];

    RDouble4D &spectumTurbulent  = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("spec_turb"));
    RDouble4D &residualTurbulent = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_turb"));
    RDouble5D &matrixTurbulentLeft    = *reinterpret_cast<RDouble5D *> (grid->GetDataPtr("mat_turbl"));

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    int nDim = GetDim();

    int mst = 0;
    int med = nTurbulentEquation - 1;

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                for (int m = mst; m <= med; ++ m)
                {
                    df_total[m] = 0.0;
                }
                
                for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
                {
                    int il1, jl1, kl1;
                    GetNsurfIndex(iSurface, il1, jl1, kl1);

                    int il = i - il1;
                    int jl = j - jl1;
                    int kl = k - kl1;

                    for (int m = mst; m <= med; ++ m)
                    {
                        dqNeighbor[m] = dqTurbulent(il, jl, kl, m);
                        mat_abc[m] = matrixTurbulentLeft(il, jl, kl, m, iSurface);
                    }

                    Turb_MxDQ(mat_abc, dqNeighbor, nTurbulentEquation, df);

                    for (int m = 0; m < nTurbulentEquation; ++ m)
                    {
                        df_total[m] += df[m];
                    }
                }
                
                for (int m = mst; m <= med; ++ m)
                {
                    RDouble dqOld = dqTurbulent(i, j, k, m);

                    dqTurbulent(i, j, k, m) = (residualTurbulent(i, j, k, m) + dRHS(i, j, k, m) + df_total[m]) / spectumTurbulent(i, j, k, m);
                    dRHS(i, j, k, m) = df_total[m];

                    sweepNormal += SQR(dqTurbulent(i, j, k, m) - dqOld);
                }
            }
        }
    }

    delete []dqNeighbor;
    delete []df_total;
    delete []df;
    delete []mat_abc;
}

void TurbSolverStrFD::SolveLUSGSBackward(Grid *gridIn, FieldProxy * dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal)
{
    StructGrid *grid = StructGridCast(gridIn);

    RDouble4D &dqTurbulent = dqProxy->GetField_STR();
    RDouble4D &dRHS        = LUplusDQ->GetField_STR();

    Param_TurbSolverStruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    RDouble *dqNeighbor = new RDouble[nTurbulentEquation];
    RDouble *df_total   = new RDouble[nTurbulentEquation];
    RDouble *df         = new RDouble[nTurbulentEquation];
    RDouble *mat_abc    = new RDouble[nTurbulentEquation];

    RDouble4D &spectumTurbulent  = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("spec_turb"));
    RDouble4D &residualTurbulent = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_turb"));
    RDouble5D &matrixTurbulentRight   = *reinterpret_cast<RDouble5D *> (grid->GetDataPtr("mat_turbr"));

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    int nDim = GetDim();

    int mst = 0;
    int med = nTurbulentEquation - 1;

    for (int k = ked; k >= kst; -- k)
    {
        for (int j = jed; j >= jst; -- j)
        {
            for (int i = ied; i >= ist; -- i)
            {
                for (int m = mst; m <= med; ++ m)
                {
                    df_total[m] = 0.0;
                }
                
                for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
                {
                    int il1, jl1, kl1;
                    GetNsurfIndex(iSurface, il1, jl1, kl1);

                    int il = i + il1;
                    int jl = j + jl1;
                    int kl = k + kl1;

                    for (int m = mst; m <= med; ++ m)
                    {
                        dqNeighbor[m] = dqTurbulent(il, jl, kl, m);
                        mat_abc[m] = matrixTurbulentRight(il, jl, kl, m, iSurface);
                    }

                    Turb_MxDQ(mat_abc, dqNeighbor, nTurbulentEquation, df);

                    for (int m = mst; m <= med; ++ m)
                    {
                        df_total[m] += df[m];
                    }
                }
                
                for (int m = mst; m <= med; ++ m)
                {
                    RDouble dqOld  = dqTurbulent(i, j, k, m);

                    dqTurbulent(i, j, k, m) = (residualTurbulent(i, j, k, m) + dRHS(i, j, k, m) + df_total[m]) / spectumTurbulent(i, j, k, m);
                    dRHS(i, j, k, m) = df_total[m];

                    sweepNormal += SQR(dqTurbulent(i, j, k, m) - dqOld);
                }
            }
        }
    }

    delete []dqNeighbor;
    delete []df_total;
    delete []df;
    delete []mat_abc;
}



}