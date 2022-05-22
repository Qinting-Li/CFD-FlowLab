#include "Geo_UnstructGrid.h"
#include "Geo_StructGrid.h"
#include "BLUSGS.h"
#include "Flux_RoeEntropyCorrection.h"
#include "Gas.h"
#include "Glb_Dimension.h"
#include "Jacobian.h"

namespace PHSPACE
{
     using namespace GAS_SPACE; 
//! This is for blusgs according wzj reference in 2000 year.
//! This is blusgsforward.
void SolveBLUSGSForwardInHouse(UnstructGrid *gridIn, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal)
{
    int nTotalCell = gridIn->GetNTotalCell();
    int *leftCellofFace  = gridIn->GetLeftCellOfFace();
    int *rightCellofFace  = gridIn->GetRightCellOfFace();

    RDouble **dq  = dqProxy->GetField_UNS();
    RDouble **res = reinterpret_cast <RDouble **> (gridIn->GetDataPtr("res"));

    //! Get grid metrics.
    RDouble *xfn  = gridIn->GetFaceNormalX();
    RDouble *yfn  = gridIn->GetFaceNormalY();
    RDouble *zfn  = gridIn->GetFaceNormalZ();
    RDouble *area = gridIn->GetFaceArea();
    RDouble *vgn  = gridIn->GetFaceNormalVelocity();

    RDouble *xcc  = gridIn->GetCellCenterX();
    RDouble *ycc  = gridIn->GetCellCenterY();
    RDouble *zcc  = gridIn->GetCellCenterZ();

    //! Get flow variables.
    RDouble **q = reinterpret_cast< RDouble ** > (gridIn->GetDataPtr("q"));

    int nm = GlobalDataBase::GetIntParaFromDB("nm");
    RDouble refGama = GlobalDataBase::GetDoubleParaFromDB("refGama");
    RDouble refMachNumber = GlobalDataBase::GetDoubleParaFromDB("refMachNumber");
    RDouble prl = GlobalDataBase::GetDoubleParaFromDB("prl");
    RDouble prt = GlobalDataBase::GetDoubleParaFromDB("prt");

    int RoeEntropyFixMethod = 3;
    RDouble RoeEntropyFixCoef1 = 0.0001, RoeEntropyFixCoef2 = 0.0001;
    GetEntropyFixMethodAndCoefficients(RoeEntropyFixMethod, RoeEntropyFixCoef1, RoeEntropyFixCoef2);

    RDouble refReNumber = GlobalDataBase::GetDoubleParaFromDB("refReNumber");

    int viscousType = GlobalDataBase::GetIntParaFromDB("viscousType");

    RDouble *viscousLaminar, *viscousTurbulence;
    if (viscousType != INVISCID)
    {
        viscousLaminar = reinterpret_cast<RDouble *> (gridIn->GetDataPtr("visl"));
        viscousTurbulence = reinterpret_cast<RDouble *> (gridIn->GetDataPtr("vist"));
    }

    RDouble cp  = 1.0 / ((refGama - 1.0) * refMachNumber * refMachNumber);
    RDouble gm2 = refGama * refMachNumber * refMachNumber;

    RDouble oprl = 1.0 / prl;
    RDouble oprt = 1.0 / prt;

    //! Some intermediate variables.
    RDouble * fluxInviscid = new RDouble [ nm ]();
    RDouble * fluxRoe = new RDouble [ nm ]();
    RDouble * fluxViscous = new RDouble [ nm ]();
    RDouble * qL     = new RDouble [ nm ]();
    RDouble * qR     = new RDouble [ nm ]();
    RDouble * qm     = new RDouble [ nm ]();
    RDouble * dQNeighbor = new RDouble [ nm ]();
    RDouble * dqTry    = new RDouble [ nm ]();
    RDouble nxyz[5] = {0};

    RDouble *Ux = new RDouble [nm]();

    int **cell2face  = gridIn->GetCell2Face();
    int * faceNumberOfEachCell = gridIn->GetFaceNumberOfEachCell();

    RDouble *rhs0   = new RDouble[nm]();

    RDouble * dqOld  = new RDouble[ nm ]();
    RDouble **deltaFlux = LUplusDQ->GetField_UNS();

    using namespace IDX;

    MAT3D< RDouble > *impmat = reinterpret_cast< MAT3D< RDouble > * >(gridIn->GetDataPtr("impmat"));

    int le,re,face;
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        RDouble **A = impmat->GetT2D(iCell);

        for (int m = 0; m < nm; ++ m)
        {
            rhs0[m] = 0.0;

            //! Back up the old dq, to compute the convergence.
            //! it is initialized by zero or res in the first sweep,
            //! then it is updated by backward sweep, that is dq*.
            dqOld[m]       = dq [m][iCell];

            //! Here, the deltaFlux is computed in Upper backward sweep!
            //! res: b      deltaFlux: Ux, which is computed in backward sweep.
            //! the dq is not the real dq, now.
            //! Although the 'dq' changed can not the right one, it dosen't matter, since 
            //! the following only using the lower neighbor cells, whose 'dq' has been updated.
            //! It is convenient for precondition transform.
            dq[m][iCell] = res[m][iCell];

            Ux[m] = deltaFlux[m][iCell];

            //! Then reset it to zero to store the Lower forward sweep!
            deltaFlux[m][iCell] = 0.0;
        }

        for (int jFace = 0; jFace < faceNumberOfEachCell[iCell]; ++ jFace)
        {
            face  = cell2face[iCell][jFace];
            le    = leftCellofFace [ face ];
            re    = rightCellofFace[ face ];

            //! One of le and re must be cell itself.
            if (le > iCell || re > iCell) continue;
            //! Now its neighboring cell belongs to lower triangular.

            nxyz[0] = xfn [ face ];
            nxyz[1] = yfn [ face ];
            nxyz[2] = zfn [ face ];
            nxyz[3] = area[ face ];
            nxyz[4] = vgn [ face ];

            if (re == iCell)
            {
                //! if the re is the iCell itself, then swap with 'le'.
                //! Now its neighboring cell 're' belongs to lower triangular.
                SWAP(le, re);

                //! Swap the normal and the moving face velocity.
                nxyz[0] = - nxyz[0];
                nxyz[1] = - nxyz[1];
                nxyz[2] = - nxyz[2];

                nxyz[4] = - nxyz[4];
            }

            for (int m = 0; m < nm; ++ m)
            {
                qL[m] = q[m][le];
                qR[m] = q[m][re];
                dQNeighbor[m] = dq[m][re];
            }

            //! refGama is fixed here. It only suit for ideal gas.
            //! The inviscid Delta Flux.
            FluxDiffFace3D(qR, dQNeighbor, nxyz, refGama, fluxInviscid, nm);

            //! The Roe Delta Flux.
            RoeDisDiff3D(qL, qR, dQNeighbor, nxyz, refGama, RoeEntropyFixMethod, RoeEntropyFixCoef1, RoeEntropyFixCoef2, fluxRoe, nm);

            for (int m = 0; m < nm; ++ m)
            {
                rhs0[m] += half * nxyz[3] * (fluxInviscid[m] -  fluxRoe[m]);
            }

            if (viscousType != INVISCID)
            {
                //! The viscous Delta Flux.
                for (int m = 0; m < nm; ++ m)
                {
                    qm[m] = half * (qL[m] + qR[m]);
                }

                RDouble viscousLaminarFace = half * (viscousLaminar[le] + viscousLaminar[re]);
                RDouble viscousTurbulenceFace = half * (viscousTurbulence[le] + viscousTurbulence[re]);

                RDouble kcp = (viscousLaminarFace * oprl + viscousTurbulenceFace * oprt) * cp;
                RDouble viscosity = viscousLaminarFace + viscousTurbulenceFace;

                RDouble distance = ABS((xcc[le] - xcc[re]) * nxyz[0] + (ycc[le] - ycc[re])*nxyz[1] + (zcc[le] - zcc[re]) * nxyz[2]);

                RDouble viscousSpectrumRadius = viscosity * area[face] / (distance * refReNumber);

                RDouble tm = gm2 * qm[IP] / qm[IR];

                VisFluxDiffOneFace3D(qm, dQNeighbor, fluxViscous, viscosity, tm, kcp, refGama, nxyz);

                for (int m = 0; m < nm; ++ m)
                {
                    rhs0[m] -= viscousSpectrumRadius * fluxViscous[m];
                }
            }
        }

        //! dq = A^(-1)*{ (b - Ux) - Lx }.
        //! Note: the 'dq' has actually initialized by the residual at the beginning of this function.
        //! the 'rhs0' is the total Delta Flux of the neighbors.   
        //! it use LU decompose to solve the results A.
        //! rhs0: Lx.

        for (int m = 0; m < nm; ++ m)
        {
            dqTry[m] = dq[m][iCell] - Ux[m] - rhs0[m];
            deltaFlux[m][iCell] += rhs0[m];
        }

        ScalarLUSolve(A, dqTry, nm);

        for (int m = 0; m < nm; ++ m)
        {
            dq[m][iCell] = dqTry[m];
            sweepNormal += SQR(dq[m][iCell] - dqOld[m]);
        }
    }

    SetGhostDQLUSGSInHouse(gridIn, dq, q);

    delete [] Ux;    Ux = NULL;
    delete [] qL;    qL = NULL;
    delete [] qR;    qR = NULL;
    delete [] qm;    qm = NULL;
    delete [] dQNeighbor;    dQNeighbor = NULL;
    delete [] fluxInviscid;    fluxInviscid = NULL;
    delete [] fluxRoe;    fluxRoe = NULL;
    delete [] fluxViscous;    fluxViscous = NULL;
    delete [] dqTry;    dqTry = NULL;
    delete [] rhs0;    rhs0 = NULL;
    delete [] dqOld;    dqOld = NULL;
}

//! This is Blusgs Backward Sweep.
void SolveBLUSGSBackwardInHouse(UnstructGrid *gridIn, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace  = grid->GetRightCellOfFace();

    //! Here, the dq is the dq-star, which has been initialized in the forward sweep stage.
    //!       that it is solve by the lower-triangle equations.
    RDouble **dq  = dqProxy->GetField_UNS();
    RDouble **res = reinterpret_cast <RDouble **> (grid->GetDataPtr("res"));

    //! Get grid metrics.
    RDouble *xfn  = grid->GetFaceNormalX();
    RDouble *yfn  = grid->GetFaceNormalY();
    RDouble *zfn  = grid->GetFaceNormalZ();
    RDouble *area = grid->GetFaceArea();
    RDouble *vgn  = grid->GetFaceNormalVelocity();

    RDouble *xcc  = grid->GetCellCenterX();
    RDouble *ycc  = grid->GetCellCenterY();
    RDouble *zcc  = grid->GetCellCenterZ();

    //! Get flow variables.
    RDouble **q  = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q" ));

    int nm = GlobalDataBase::GetIntParaFromDB("nm");
    RDouble refGama = GlobalDataBase::GetDoubleParaFromDB("refGama");
    RDouble refMachNumber = GlobalDataBase::GetDoubleParaFromDB("refMachNumber");
    RDouble prl = GlobalDataBase::GetDoubleParaFromDB("prl");
    RDouble prt = GlobalDataBase::GetDoubleParaFromDB("prt");

    int RoeEntropyFixMethod = 3;
    RDouble RoeEntropyFixCoef1 = 0.0001, RoeEntropyFixCoef2 = 0.0001;
    GetEntropyFixMethodAndCoefficients(RoeEntropyFixMethod, RoeEntropyFixCoef1, RoeEntropyFixCoef2);

    RDouble refReNumber = GlobalDataBase::GetDoubleParaFromDB("refReNumber");

    int viscousType = GlobalDataBase::GetIntParaFromDB("viscousType");

    RDouble *viscousLaminar, *viscousTurbulence;
    if (viscousType != INVISCID)
    {
        viscousLaminar = reinterpret_cast<RDouble *> (grid->GetDataPtr("visl"));
        viscousTurbulence = reinterpret_cast<RDouble *> (grid->GetDataPtr("vist"));
    }

    //! Compute some intermediate variables to simplify calculations.
    RDouble cp  = 1.0 / ((refGama - 1.0) * refMachNumber * refMachNumber);
    RDouble gm2 = refGama * refMachNumber * refMachNumber;

    RDouble oprl = 1.0 / prl;
    RDouble oprt = 1.0 / prt;

    //! Some intermediate variables.
    RDouble *fluxInviscid = new RDouble [ nm ];
    RDouble *fluxRoe = new RDouble [ nm ];
    RDouble *fluxViscous = new RDouble [ nm ];
    RDouble *qL     = new RDouble [ nm ];
    RDouble *qR     = new RDouble [ nm ];
    RDouble *qm     = new RDouble [ nm ];
    RDouble *dQNeighbor = new RDouble [ nm ];
    RDouble *dqTry    = new RDouble [ nm ];
    RDouble nxyz[5];
    RDouble *Lx = new RDouble [nm];

    int **cell2face  = grid->GetCell2Face();
    int * faceNumberOfEachCell = grid->GetFaceNumberOfEachCell();

    RDouble *rhs0   = new RDouble[nm];

    RDouble * dqOld  = new RDouble[ nm ];
    RDouble **deltaFlux = LUplusDQ->GetField_UNS();

    using namespace IDX;
    int le,re,face;

    MAT3D< RDouble > *impmat = reinterpret_cast< MAT3D< RDouble > * >(grid->GetDataPtr("impmat"));


    for (int iCell = nTotalCell - 1; iCell >= 0; -- iCell)
    {
        RDouble **A = impmat->GetT2D(iCell);

        for (int m = 0; m < nm; ++ m)
        {
            rhs0[m] = 0.0;
        }

        for (int m = 0; m < nm; ++ m)
        {
            dqOld[m]       = dq [m][iCell];
            dq  [m][iCell] = res[m][iCell];
            Lx[m] = deltaFlux[m][iCell];
            deltaFlux[m][iCell] = 0.0;
        }

        for (int j = 0; j < faceNumberOfEachCell[iCell]; ++ j)
        {
            face  = cell2face[iCell][j];
            le    = leftCellofFace [ face ];
            re    = rightCellofFace[ face ];
            //! One of le and re must be cell itself.
            if (le < iCell || re < iCell) continue;

            //! Now its neighboring cell belongs to upper triangular.
            nxyz[0] = xfn [ face ];
            nxyz[1] = yfn [ face ];
            nxyz[2] = zfn [ face ];
            nxyz[3] = area [ face ];
            nxyz[4] = vgn[ face ];

            if (re == iCell)
            {
                SWAP(le, re);
                nxyz[0] = - nxyz[0];
                nxyz[1] = - nxyz[1];
                nxyz[2] = - nxyz[2];
                nxyz[4] = - nxyz[4];
            }

            for (int m = 0; m < nm; ++ m)
            {
                qL[m] = q[m][le];
                qR[m] = q[m][re];
                dQNeighbor[m] = dq[m][re];
            }


            //! The inviscid Delta Flux.
            FluxDiffFace3D(qR, dQNeighbor, nxyz, refGama, fluxInviscid, nm);

            //! The Roe Delta Flux.
            RoeDisDiff3D(qL, qR, dQNeighbor, nxyz, refGama, RoeEntropyFixMethod, RoeEntropyFixCoef1, RoeEntropyFixCoef2, fluxRoe, nm);

            for (int m = 0; m < nm; ++ m)
            {
                rhs0[m] += half * nxyz[3] * (fluxInviscid[m] -  fluxRoe[m]);
            }

            if (viscousType != INVISCID)
            {
                //! Viscous jacobian.
                for (int m = 0; m < nm; ++ m)
                {
                    qm[m] = half * (qL[m] + qR[m]);
                }

                RDouble viscousLaminarFace = half * (viscousLaminar[le] + viscousLaminar[re]);
                RDouble viscousTurbulenceFace = half * (viscousTurbulence[le] + viscousTurbulence[re]);

                RDouble kcp = (viscousLaminarFace * oprl + viscousTurbulenceFace * oprt) * cp;
                RDouble viscosity = viscousLaminarFace + viscousTurbulenceFace;

                RDouble Dis = ABS((xcc[le] - xcc[re]) * nxyz[0] + (ycc[le] - ycc[re])*nxyz[1] + (zcc[le] - zcc[re]) * nxyz[2]);

                RDouble delta = viscosity * area[face] / (Dis * refReNumber);

                RDouble tm = gm2 * qm[IP] / qm[IR];

                VisFluxDiffOneFace3D(qm, dQNeighbor, fluxViscous, viscosity, tm, kcp, refGama, nxyz);

                for (int m = 0; m < nm; ++ m)
                {
                    rhs0[m] -= delta * fluxViscous[m];
                }
            }
        }

        //! Note: the 'dq' has been updated by the forward sweep.
        //! the 'rhs0' is the total Delta Flux of the neighbors.
        //! x = A^(-1)*{(b - LX) - Ux}.
        //! rhs0: Ux.

        for (int m = 0; m < nm; ++ m)
        {
            dqTry[m] = dq[m][iCell] - Lx[m] - rhs0[m];
            deltaFlux[m][iCell] += rhs0[m];
        }

        ScalarLUSolve(A, dqTry, nm);

        for (int m = 0; m < nm; ++ m)
        {
            dq[m][iCell] = dqTry[m];
            sweepNormal += SQR(dq[m][iCell] - dqOld[m]);
        }
    }

    delete [] qL;    qL = NULL;
    delete [] qR;    qR = NULL;
    delete [] qm;    qm = NULL;
    delete [] dQNeighbor;    dQNeighbor = NULL;
    delete [] fluxInviscid;    fluxInviscid = NULL;
    delete [] fluxRoe;    fluxRoe = NULL;
    delete [] fluxViscous;    fluxViscous = NULL;
    delete [] dqTry;    dqTry = NULL;
    delete [] rhs0;    rhs0 = NULL;
    delete [] Lx;    Lx = NULL;
    delete [] dqOld;    dqOld = NULL;
}

void SolveBLUSGSForwardInHouse(StructGrid *gridIn, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal)
{
    FieldProxy * rightHandSideProxy = GetFieldProxy(gridIn, "res");
    FillField(dqProxy, rightHandSideProxy);
    delete rightHandSideProxy;
    //! To compute and save the diagonal matrix.
    ComputeDiagonalMatrice(gridIn);

    //! Obtain the basic variables.
    StructGrid *grid = StructGridCast(gridIn);
    RDouble4D &dq = dqProxy->GetField_STR();
    RDouble4D &gxfn = *(grid->GetFaceNormalX());
    RDouble4D &gyfn = *(grid->GetFaceNormalY());
    RDouble4D &gzfn = *(grid->GetFaceNormalZ());
    RDouble3D &gvol = *(grid->GetCellVolume());
    RDouble4D &gvgn = *(grid->GetFaceNormalVelocity());
    RDouble4D &garea = *(grid->GetFaceArea());

    //! Obtain the spectrum radius of inviscid flux and viscous flux.
    RDouble4D &inviscidSpectrumRadius = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("invSpectralRadius"));
    RDouble4D &viscousSpectrumRadius = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("visSpectralRadius"));
    RDouble4D &primitiveVariables = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));
    RDouble4D &primitiveTemperatures = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("t"));
    RDouble3D &specificHeatRatio = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("gama"));
    RDouble5D &diagonalMatrice = *reinterpret_cast<RDouble5D *> (grid->GetDataPtr("diagonalMatrice"));

    int viscousType = GlobalDataBase::GetIntParaFromDB("viscousType");
    RDouble3D &viscosityLaminar = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("visl"));
    RDouble3D &viscosityTurbulence = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("vist"));

    int nLaminar = GlobalDataBase::GetIntParaFromDB("nl");
    int nChemical = GlobalDataBase::GetIntParaFromDB("nchem");
    int nTemperatureModel = GlobalDataBase::GetIntParaFromDB("ntmodel");
    int iSimplifyViscousTerm = GlobalDataBase::GetIntParaFromDB("iSimplifyViscousTerm");

    int nElectronIndex = 0;
    if(nChemical)nElectronIndex = GlobalDataBase::GetIntParaFromDB("nElectronIndex");

    //! Compute the number of equations.
    int nEquation = nLaminar + nChemical + nTemperatureModel - 1;

    //! Allocate the memories.
    Range M(0, nEquation - 1);
    RDouble2D *tempMatrix = new RDouble2D(M, M, fortranArray);
    RDouble2D *diagonalMatrix = new RDouble2D(M, M, fortranArray);
    RDouble2D &tempMatrix2D = *reinterpret_cast<RDouble2D *> (tempMatrix);
    RDouble2D &diagonalMatrix2D = *reinterpret_cast<RDouble2D *> (diagonalMatrix);

    RDouble *primitiveVars = new RDouble[nEquation];
    RDouble *deltaQ = new RDouble[nEquation];
    RDouble *deltaFlux = new RDouble[nEquation];
    RDouble *temperature = new RDouble[nEquation];
    RDouble *rhsMatrix = new RDouble[nEquation];
    RDouble *residuals = new RDouble[nEquation];

    //! Obtain the indexes of grid nodes in matrix.
    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    //! Iteration: the forward sweep.
    int nDim = GetDim();
    int iL = 0, jL = 0, kL = 0;
    RDouble gama = 0.0, laminarViscosity = 0.0, turbulentViscosity = 0.0;
    RDouble xfn = 0.0, yfn = 0.0, zfn = 0.0, area = 0.0, vgn = 0.0, vol = 0.0;
    RDouble invSpectrumRadius = 0.0, visSpectrumRadius = 0.0;

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                //! Construct the diagonal matrix.
                for (int m = 0; m < nEquation; ++ m)
                {
                    for (int n = 0; n < nEquation; ++ n)
                    {
                        diagonalMatrix2D(m, n) = diagonalMatrice(i, j, k, m, n);
                    }
                }

                for (int m = 0; m < nEquation; ++ m)
                {
                    residuals[m] = dq(i, j, k, m);
                    rhsMatrix[m] = 0.0;
                }

                //! Construct the Lower-Block Jacobian matrix called L*dQ.
                for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
                {
                    int il1, jl1, kl1;
                    GetNsurfIndex(iSurface, il1, jl1, kl1);

                    //! Obtain the indexes of the neighbor cell in the lower section of the block.
                    iL = i - il1;
                    jL = j - jl1;
                    kL = k - kl1;
                    //! Obtain the primitive variables.
                    for (int m = 0; m < nEquation; ++ m)
                    {
                        primitiveVars[m] = primitiveVariables(iL, jL, kL, m);
                        deltaQ[m] = dq(iL, jL, kL, m);
                    }
                    gama = specificHeatRatio(iL, jL, kL);
                    if (viscousType > INVISCID)
                    {
                        laminarViscosity = viscosityLaminar(iL, jL, kL);
                        turbulentViscosity = viscosityTurbulence(iL, jL, kL);
                    }
                    else
                    {
                        laminarViscosity = 0.0;
                        turbulentViscosity = 0.0;
                    }
                    for (int m = 0; m < nTemperatureModel; ++ m)
                    {
                        temperature[m] = primitiveTemperatures(iL, jL, kL, m);
                    }
                    //! Obtain the normal vector of the face.
                    xfn = gxfn(i, j, k, iSurface);
                    yfn = gyfn(i, j, k, iSurface);
                    zfn = gzfn(i, j, k, iSurface);
                    vgn = gvgn(i, j, k, iSurface);
                    area = garea(i, j, k, iSurface);

                    //! Compute spectrum radius of inviscid flux and viscous flux.
                    if (iSimplifyViscousTerm == 1)
                    {
                        invSpectrumRadius = inviscidSpectrumRadius(iL, jL, kL, iSurface - 1);
                        visSpectrumRadius = viscousSpectrumRadius(iL, jL, kL, iSurface - 1);
                    }
                    else
                    {
                        invSpectrumRadius = ObtainInviscidSpectrumRadius(primitiveVars, xfn, yfn, zfn, area, vgn, gama);
                        vol = 0.5 * (gvol(i, j, k) + gvol(iL, jL, kL));
                        visSpectrumRadius = ObtainViscousSpectrumRadius(primitiveVars, area, vol, gama, laminarViscosity, turbulentViscosity);
                    }

                    //! To compute 0.5 * (M - R*I).
                    if (nTemperatureModel > 1)
                    {
                        ChemicalMXDQ(primitiveVars, temperature, xfn, yfn, zfn, area, vgn, deltaQ, deltaFlux, 5, nLaminar, nTemperatureModel, invSpectrumRadius, 1,nElectronIndex);
                    }
                    else
                    {
                        MXDQ_STD(primitiveVars, xfn, yfn, zfn, area, vgn, gama, deltaQ, deltaFlux, 5, nLaminar, invSpectrumRadius, 1);
                    }
                    for (int m = 0; m < nEquation; ++ m)
                    {
                        rhsMatrix[m] += -deltaFlux[m];
                    }

                    //! To compute 0.5 * (Mv - Rv*I).
                    if (viscousType > INVISCID)
                    {
                        if (iSimplifyViscousTerm == 1)      //! The Jacobian matrix Mv is set zero.
                        {
                            for (int m = 0; m < nEquation; ++ m)
                            {
                                rhsMatrix[m] += -0.5 * visSpectrumRadius * deltaQ[m];
                            }
                        }
                        else     //! The Jacobian matrix Mv is not set zero.
                        {
                            //MVXDQ_STD_LP(primitiveVars, xfn, yfn, zfn, area, vol, laminarViscosity, turbulentViscosity, gama, deltaQ, deltaFlux, 5, nLaminar, visSpectrumRadius, 1);
                            //for (int m = 0; m < nEquation; ++ m)
                            //{
                            //    rhsMatrix[m] += -deltaFlux[m];
                            //}

                            //! To obtain the Jacobian matrix of the viscous flux.
                            ObtainViscousMatrix(primitiveVars, temperature, xfn, yfn, zfn, area, vol, vgn, gama, laminarViscosity, turbulentViscosity, tempMatrix);
                            for (int m = 0; m < nEquation; ++ m)
                            {
                                RDouble fSum = 0.0;
                                for (int n = 0; n < nEquation; ++ n)
                                {
                                    fSum += -tempMatrix2D(m, n) * deltaQ[n];
                                }
                                rhsMatrix[m] += 0.5 * (fSum - visSpectrumRadius * deltaQ[m]);
                            }
                        }
                    }
                }

                //! Construct the right-hand side.
                for (int m = 0; m < nEquation; ++ m)
                {
                    rhsMatrix[m] = residuals[m] - rhsMatrix[m];
                }

                //! Iteration: compute the linear equation D*dQ = R - L*dQ using the SOR method.
                SolveLinearEquation(diagonalMatrix, rhsMatrix, nEquation, deltaQ);

                //! Update the values.
                for (int m = 0; m < nEquation; ++ m)
                {
                    dq(i, j, k, m) = deltaQ[m];
                }

            }
        }
    }

    //! Deallocate the memories.
    delete tempMatrix;
    delete diagonalMatrix;
    delete []primitiveVars;
    delete []deltaQ;
    delete []deltaFlux;
    delete []temperature;
    delete []rhsMatrix;
    delete []residuals;
}

void SolveBLUSGSBackwardInHouse(StructGrid *gridIn, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal)
{
    //! Obtain the basic variables.
    StructGrid *grid = StructGridCast(gridIn);
    RDouble4D &dq = dqProxy->GetField_STR();
    RDouble4D &gxfn = *(grid->GetFaceNormalX());
    RDouble4D &gyfn = *(grid->GetFaceNormalY());
    RDouble4D &gzfn = *(grid->GetFaceNormalZ());
    RDouble3D &gvol = *(grid->GetCellVolume());
    RDouble4D &gvgn = *(grid->GetFaceNormalVelocity());
    RDouble4D &garea = *(grid->GetFaceArea());

    //! Obtain the spectrum radius of inviscid flux and viscous flux.
    RDouble4D &inviscidSpectrumRadius = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("invSpectralRadius"));
    RDouble4D &viscousSpectrumRadius = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("visSpectralRadius"));
    RDouble4D &primitiveVariables = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));
    RDouble4D &primitiveTemperatures = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("t"));
    RDouble3D &specificHeatRatio = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("gama"));
    RDouble5D &diagonalMatrice = *reinterpret_cast<RDouble5D *> (grid->GetDataPtr("diagonalMatrice"));

    int viscousType = GlobalDataBase::GetIntParaFromDB("viscousType");
    RDouble3D &viscosityLaminar = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("visl"));
    RDouble3D &viscosityTurbulence = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("vist"));

    int nLaminar = GlobalDataBase::GetIntParaFromDB("nl");
    int nChemical = GlobalDataBase::GetIntParaFromDB("nchem");
    int nTemperatureModel = GlobalDataBase::GetIntParaFromDB("ntmodel");
    int iSimplifyViscousTerm = GlobalDataBase::GetIntParaFromDB("iSimplifyViscousTerm");

    int nElectronIndex = 0;
    if(nChemical)nElectronIndex = GlobalDataBase::GetIntParaFromDB("nElectronIndex");

    //! Compute the number of equations.
    int nEquation = nLaminar + nChemical + nTemperatureModel - 1; 

    //! Allocate the memories.
    Range M(0, nEquation - 1);
    RDouble2D *tempMatrix = new RDouble2D(M, M, fortranArray);
    RDouble2D *diagonalMatrix = new RDouble2D(M, M, fortranArray);
    RDouble2D &tempMatrix2D = *reinterpret_cast<RDouble2D *> (tempMatrix);
    RDouble2D &diagonalMatrix2D = *reinterpret_cast<RDouble2D *> (diagonalMatrix);

    RDouble *primitiveVars = new RDouble[nEquation];
    RDouble *deltaQ = new RDouble[nEquation];
    RDouble *deltaFlux = new RDouble[nEquation];
    RDouble *temperature = new RDouble[nEquation];
    RDouble *rhsMatrix = new RDouble[nEquation];

    //! Obtain the indexes of grid nodes in matrix.
    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    //! Iteration: the forward sweep.
    int nDim = GetDim();
    int iR = 0, jR = 0, kR = 0;
    RDouble gama = 0.0, laminarViscosity = 0.0, turbulentViscosity = 0.0;
    RDouble xfn = 0.0, yfn = 0.0, zfn = 0.0, area = 0.0, vgn = 0.0, vol = 0.0;
    RDouble invSpectrumRadius = 0.0, visSpectrumRadius = 0.0;

    for (int k = ked; k >= kst; -- k)
    {
        for (int j = jed; j >= jst; -- j)
        {
            for (int i = ied; i >= ist; -- i)
            {
                //! Construct the diagonal matrix.
                for (int m = 0; m < nEquation; ++ m)
                {
                    for (int n = 0; n < nEquation; ++ n)
                    {
                        diagonalMatrix2D(m, n) = diagonalMatrice(i, j, k, m, n);
                    }
                }

                for (int m = 0; m < nEquation; ++ m)
                {
                    rhsMatrix[m] = 0.0;
                }

                //! Construct the Upper-Block Jacobian matrix called U*dQ.
                for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
                {
                    int il1, jl1, kl1;
                    GetNsurfIndex(iSurface, il1, jl1, kl1);

                    //! Obtain the indexes of the neighbor cell in the upper section of the block.
                    iR = i + il1;
                    jR = j + jl1;
                    kR = k + kl1;
                    //! Obtain the primitive variables.
                    for (int m = 0; m < nEquation; ++ m)
                    {
                        primitiveVars[m] = primitiveVariables(iR, jR, kR, m);
                        deltaQ[m] = dq(iR, jR, kR, m);
                    }
                    gama = specificHeatRatio(iR, jR, kR);
                    if (viscousType > INVISCID)
                    {
                        laminarViscosity = viscosityLaminar(iR, jR, kR);
                        turbulentViscosity = viscosityTurbulence(iR, jR, kR);
                    }
                    else
                    {
                        laminarViscosity = 0.0;
                        turbulentViscosity = 0.0;
                    }
                    for (int m = 0; m < nTemperatureModel; ++ m)
                    {
                        temperature[m] = primitiveTemperatures(iR, jR, kR, m);
                    }
                    xfn = gxfn(iR, jR, kR, iSurface);
                    yfn = gyfn(iR, jR, kR, iSurface);
                    zfn = gzfn(iR, jR, kR, iSurface);
                    vgn = gvgn(iR, jR, kR, iSurface);
                    area = garea(iR, jR, kR, iSurface);

                    //! Compute spectrum radius of inviscid flux and viscous flux.
                    if (iSimplifyViscousTerm == 1)
                    {
                        invSpectrumRadius = inviscidSpectrumRadius(iR, jR, kR, iSurface - 1);
                        visSpectrumRadius = viscousSpectrumRadius(iR, jR, kR, iSurface - 1);
                    }
                    else
                    {
                        invSpectrumRadius = ObtainInviscidSpectrumRadius(primitiveVars, xfn, yfn, zfn, area, vgn, gama);
                        vol = 0.5 * (gvol(iR, jR, kR) + gvol(iR + il1, jR + jl1, kR + kl1));
                        visSpectrumRadius = ObtainViscousSpectrumRadius(primitiveVars, area, vol, gama, laminarViscosity, turbulentViscosity);
                    }

                    //! To compute 0.5 * (M - R*I).
                    if (nTemperatureModel > 1)
                    {
                        ChemicalMXDQ(primitiveVars, temperature, xfn, yfn, zfn, area, vgn, deltaQ, deltaFlux, 5, nLaminar, nTemperatureModel, invSpectrumRadius, -1,nElectronIndex);
                    }
                    else
                    {
                        MXDQ_STD(primitiveVars, xfn, yfn, zfn, area, vgn, gama, deltaQ, deltaFlux, 5, nLaminar, invSpectrumRadius, -1);
                    }
                    for (int m = 0; m < nEquation; ++ m)
                    {
                        rhsMatrix[m] += deltaFlux[m];
                    }

                    //! To compute 0.5 * (Mv - Rv*I).
                    if (viscousType > INVISCID)
                    {
                        if (iSimplifyViscousTerm == 1)
                        {
                            for (int m = 0; m < nEquation; ++ m) // The Jacobian matrix Mv is set zero.
                            {
                                rhsMatrix[m] += -0.5 * visSpectrumRadius * deltaQ[m];
                            }
                        }
                        else      //! The Jacobian matrix Mv is not set zero.
                        {
                            //MVXDQ_STD_LP(primitiveVars, xfn, yfn, zfn, area, vol, laminarViscosity, turbulentViscosity, gama, deltaQ, deltaFlux, 5, nLaminar, visSpectrumRadius, -1);
                            //for (int m = 0; m < nEquation; ++ m)
                            //{
                            //    rhsMatrix[m] += deltaFlux[m];
                            //}

                            //! To obtain the Jacobian matrix of the viscous flux.
                            ObtainViscousMatrix(primitiveVars, temperature, xfn, yfn, zfn, area, vol, vgn, gama, laminarViscosity, turbulentViscosity, tempMatrix);
                            for (int m = 0; m < nEquation; ++ m)
                            {
                                RDouble fSum = 0.0;
                                for (int n = 0; n < nEquation; ++ n)
                                {
                                    fSum += tempMatrix2D(m, n) * deltaQ[n];
                                }
                                rhsMatrix[m] += 0.5 * (fSum - visSpectrumRadius * deltaQ[m]);
                            }
                        }
                    }
                }

                //! Construct the right-hand side.
                for (int m = 0; m < nEquation; ++ m)
                {
                    RDouble fsum = 0.0;
                    for (int n = 0; n < nEquation; ++ n)
                    {
                        fsum += diagonalMatrix2D(m, n) * dq(i, j, k, n);
                    }
                    rhsMatrix[m] = fsum - rhsMatrix[m];
                }

                //! Iteration: compute the linear equation D*dQ = D*dQ* - U*dQ using the SOR method.
                SolveLinearEquation(diagonalMatrix, rhsMatrix, nEquation, deltaQ);

                //! Update the values.
                for (int m = 0; m < nEquation; ++ m)
                {
                    dq(i, j, k, m) = deltaQ[m];
                }

            }
        }
    }

    //! Deallocate the memories.
    delete tempMatrix;
    delete diagonalMatrix;
    delete []primitiveVars;
    delete []deltaQ;
    delete []deltaFlux;
    delete []temperature;
    delete []rhsMatrix;
}

void FluxDiffFace3D(RDouble *prim, RDouble *dQ, RDouble *nxyz, double gama, RDouble *df, int nm)
{
    using namespace IDX;

    double nx, ny, nz;
    nx = nxyz[0];
    ny = nxyz[1];
    nz = nxyz[2];

    double r,u,v,w,p,H,V2,Vn;
    r = prim[IR];
    u = prim[IU];
    v = prim[IV];
    w = prim[IW];
    p = prim[IP];

    V2 = u * u + v * v + w * w;
    Vn = nx * u + ny * v + nz * w;
    H = (gama / (gama-1)) * (p / r) + half * V2;

    RDouble b1 = -Vn * dQ[IR] + nx * dQ[IRU] + ny * dQ[IRV] + nz * dQ[IRW];
    RDouble b2 = (gama-1)*(half * V2 * dQ[IR] - u * dQ[IRU] - v * dQ[IRV] - w * dQ[IRW] + dQ[IRE]);

    df[IR ] = Vn * dQ[IR ] + b1;
    df[IRU] = Vn * dQ[IRU] + u * b1 + nx  * b2;
    df[IRV] = Vn * dQ[IRV] + v * b1 + ny  * b2;
    df[IRW] = Vn * dQ[IRW] + w * b1 + nz  * b2;
    df[IRE] = Vn * dQ[IRE] + H * b1 + Vn * b2;
}


void RoeDisDiff3D(RDouble *ql, RDouble *qr, RDouble *dQ, RDouble *nxyz, RDouble gama, int RoeEntropyFixMethod,
                  RDouble RoeEntropyFixCoef1, RDouble RoeEntropyFixCoef2, RDouble *flux, int nm)
{
    RDouble nx, ny, nz;

    RDouble r,u,v,w,h,c,Vn, absVel,halfV2,c2,ratio,coef;
    RDouble rl,ul,vl,wl,pl,hl,v2l;
    RDouble rr,ur,vr,wr,pr,hr,v2r;
    RDouble z0,z1, X1, X2, b1, b2;

    RDouble eigv1,eigv2,eigv3/*,eigv_max,tmp0,tmp1,tmp2*/;

    using namespace IDX;

    //! Get the face normal.
    nx = nxyz[0];
    ny = nxyz[1];
    nz = nxyz[2];   

    //! Get the the left and right cell value of the face for roe flux.
    rl = ql[IR];
    ul = ql[IU];
    vl = ql[IV];
    wl = ql[IW];
    pl = ql[IP];

    rr = qr[IR];
    ur = qr[IU];
    vr = qr[IV];
    wr = qr[IW];
    pr = qr[IP];

    if (rl <= 1.0e-5) rl = 1.0e-5;
    if (rr <= 1.0e-5) rr = 1.0e-5;

    v2l = ul * ul + vl * vl + wl * wl;
    v2r = ur * ur + vr * vr + wr * wr;

    RDouble gamm1  = gama - 1.0;
    RDouble gaminv = gama / gamm1;
    hl  = gaminv * pl / rl + half * v2l;
    hr  = gaminv * pr / rr + half * v2r;

    ratio = sqrt(rr / rl);
    coef  = 1.0 / (1.0 + ratio);

    r = sqrt(rl * rr);
    u = (ul + ur * ratio) * coef;
    v = (vl + vr * ratio) * coef;
    w = (wl + wr * ratio) * coef;
    h = (hl + hr * ratio) * coef;

    halfV2 = half * (u * u + v * v + w * w);
	absVel = sqrt(two * halfV2);
    Vn = nx * u + ny * v + nz * w;

    c2 = gamm1 * (h - halfV2);
    c = sqrt(ABS(c2));

    eigv1 = ABS(Vn);
    eigv2 = ABS(Vn + c);
    eigv3 = ABS(Vn - c);

    //! Limit the eigen values.
    Flux_RoeEntropyCorrection(eigv1, eigv2, eigv3, 0.0, 0.0, absVel, Vn, c, RoeEntropyFixMethod, RoeEntropyFixCoef1, RoeEntropyFixCoef2);

    z0      = u  * dQ[IRU] + v  * dQ[IRV] + w  * dQ[IRW];
    z1      = nx * dQ[IRU] + ny * dQ[IRV] + nz * dQ[IRW];

    X1 = (2*eigv1 - eigv2 - eigv3)/(2 * c2);
    X2 = (eigv2 - eigv3)/(2 * c);

    b1 = -halfV2 * dQ[IR] * gamm1 * X1 -Vn * dQ[IR] * X2 + z0 * gamm1 * X1 + z1 * X2 - dQ[IRE] * gamm1*X1;
    b2 = Vn * c2 * dQ[IR] * X1 +halfV2 * dQ[IR] * gamm1 * X2 - c2 * z1 * X1 - z0 * gamm1 * X2 + dQ[IRE] * gamm1 * X2;

    flux[IR ] = eigv1 * dQ[IR] + b1 * 1 + b2 * 0;
    flux[IRU] = eigv1 * dQ[IRU] + b1 * u + b2 * nx;
    flux[IRV] = eigv1 * dQ[IRV] + b1 * v + b2 * ny;
    flux[IRW] = eigv1 * dQ[IRW] + b1 * w + b2 * nz;
    flux[IRE] = (eigv2+eigv3) / 2 *dQ[IRE] + b1 * halfV2 + b2 * Vn + c2*z0*X1
        +c2*z1*X2/gamm1 - halfV2 * c2 * dQ[IR] * X1 - c2 * Vn * dQ[IR] * X2/gamm1;
}

void VisFluxDiffOneFace3D(RDouble *q, RDouble *dQ, RDouble *flux, RDouble vis, RDouble T, RDouble kcp, RDouble gama, RDouble *nxyz)
{
    using namespace IDX;
    RDouble r = q[IR];
    RDouble u = q[IU];
    RDouble v = q[IV];
    RDouble w = q[IW];
    RDouble p = q[IP];

    RDouble nx = nxyz[0];
    RDouble ny = nxyz[1];
    RDouble nz = nxyz[2];

    RDouble Vn = u * nx + v * ny + w * nz;
    RDouble V2 = u * u + v * v + w * w;

    RDouble t0 = nx * dQ[IRU] + ny * dQ[IRV] + nz * dQ[IRW];
    RDouble t1 = q[IU] * dQ[IRU] + q[IV] * dQ[IRV] + q[IW] * dQ[IRW];

    RDouble b1 = -dQ[IR] / r;
    RDouble b2 = 1 / r;
    RDouble b3 = (t0 - dQ[IR] * Vn) / 3 / r;

    RDouble y0 = kcp*T/vis/r;
    RDouble y1 = kcp*T*(gama-1)/vis/p;

    flux[IR ] = 0.0;
    flux[IRU] = b1 * u + b2 * dQ[IRU] + b3 * nx;
    flux[IRV] = b1 * v + b2 * dQ[IRV] + b3 * ny;
    flux[IRW] = b1 * w + b2 * dQ[IRW] + b3 * nz;
    flux[IRE] = u * flux[IRU] + v * flux[IRV] + w * flux[IRW] + (y1*half*V2-y0) * dQ[IR] - y1*t1 + y1 * dQ[IRE] ;
}

void GetEntropyFixMethodAndCoefficients(int &RoeEntropyFixMethod, RDouble &RoeEntropyFixCoef1, RDouble &RoeEntropyFixCoef2)
{
    RoeEntropyFixMethod = GlobalDataBase::GetIntParaFromDB("roeEntropyFixMethod");

    RDouble entropyAcoustic = 0.01, entropyConvective = 0.01;
    RDouble refMachNumber = GlobalDataBase::GetDoubleParaFromDB("refMachNumber");

    RDouble roeEntropyScale = GlobalDataBase::GetDoubleParaFromDB("roeEntropyScale");
    int ifLowSpeedPrecon = GlobalDataBase::GetIntParaFromDB("ifLowSpeedPrecon");
    if (ifLowSpeedPrecon)
    {
        //! Incompressible, precondition.
        entropyAcoustic = 1.0e-8;
    }
    else
    {
        //! Compressible.
        RDouble x0, x1, y0, y1;
        if (refMachNumber <= 0.6)
        {
            entropyAcoustic = 0.0001;
        }
        else if (refMachNumber <= 2.0)
        {
            x0 = 0.6;    x1 = 2.0;
            y0 = 0.0001; y1 = 0.01;
            entropyAcoustic = y0 * (refMachNumber - x1) / (x0 - x1) + y1 * (refMachNumber - x0) / (x1 - x0);
        }
        else if (refMachNumber <= 5.0)
        {
            x0 = 2.0;  x1 = 5.0;
            y0 = 0.01; y1 = 0.5;
            entropyAcoustic = y0 * (refMachNumber - x1) / (x0 - x1) + y1 * (refMachNumber - x0) / (x1 - x0);
        }
        else
        {
            entropyAcoustic = 1.0;
        }  

        //! Do not scale the entropy for incompressible precondition.
        entropyAcoustic *= roeEntropyScale;
    }

    entropyConvective = entropyAcoustic;
    //PrintToWindow("  Entropy fix coefficient for ROE scheme: ", entropyAcoustic, "\n");

    RoeEntropyFixCoef1 = entropyAcoustic;
    RoeEntropyFixCoef2 = entropyConvective;
}

void SetGhostDQLUSGSInHouse(UnstructGrid *gridIn, RDouble **dq, RDouble **q)
{
    int nBoundFace = gridIn->GetNBoundFace();
    int *leftCellofFace  = gridIn->GetLeftCellOfFace();
    int *rightCellofFace = gridIn->GetRightCellOfFace();

    //! Get grid metrics.
    RDouble *xfn  = gridIn->GetFaceNormalX();
    RDouble *yfn  = gridIn->GetFaceNormalY();
    RDouble *zfn  = gridIn->GetFaceNormalZ();

    UnstructBCSet *unstructBCRegion = gridIn->GetUnstructBCRegion();
    int *bcRegionIDofBCFace = unstructBCRegion->GetBCFaceInfo();

    RDouble *gama = reinterpret_cast<RDouble *> (gridIn->GetDataPtr("gama"));
    RDouble **t   = reinterpret_cast<RDouble **> (gridIn->GetDataPtr("t"));

    int nl    = GlobalDataBase::GetIntParaFromDB("nl");
    int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
    int ntmodel = GlobalDataBase::GetIntParaFromDB("ntmodel");
    int nEquation = nl + nchem + ntmodel - 1;

    int viscousType = GlobalDataBase::GetIntParaFromDB("viscousType");

    RDouble *primitiveVariable = new RDouble[nEquation];
    RDouble *Q = new RDouble[nEquation];

    RDouble gam1, temperature[3];
    RDouble rm, um, vm, wm, pm, ub, vb, wb;
    RDouble rvn, ru, rv, rw, p_new;

    using namespace IDX;

    int le, re, bcType;
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        UnstructBC *bcRegion = unstructBCRegion->GetBCRegion(bcRegionIDofBCFace[iFace]);

        le = leftCellofFace[iFace];
        re = rightCellofFace[iFace];

        //! Conservative Variables on level n+1.

        for (int m = 0; m < nEquation; ++ m)
        {
            primitiveVariable[m] = q[m][le];
        }

        temperature[ITT] = t[ITT][le];
        temperature[ITV] = temperature[ITT];
        temperature[ITE] = temperature[ITT];
        if (ntmodel > 1)     //! Consider Two-Temperature model and Three-Temperature model.
        {
            temperature[ITV] = t[ITV][le];
            temperature[ITE] = temperature[ITV];
            if (ntmodel == 3)
            {
                temperature[ITE] = t[ITE][le];
            }
        }
        gas->Primitive2Conservative(primitiveVariable, gama[le], temperature[ITV], temperature[ITE], Q);

        for (int m = 0; m < nl; ++ m)
        {
            Q[m] += dq[m][le];
        }
        for (int m = 0; m < ntmodel - 1; m++)
        {
            Q[nl + nchem + m] += dq[nl + nchem + m][le];
        }

        //! Non-conservative variables in re in level n.
        rm = q[IR][re];
        um = q[IU][re];
        vm = q[IV][re];
        wm = q[IW][re];
        pm = q[IP][re];

        int bcType = bcRegion->GetBCType();

        for (int m = 0; m < nl; ++ m)
        {
            dq[m][re] = dq[m][le];
        }
        for (int m = 0; m < ntmodel - 1; m++)
        {
            dq[nl + nchem + m][re] = dq[nl + nchem + m][le];
        }

        if (bcType == PHENGLEI::SYMMETRY || (IsWall(bcType) && !viscousType))
        {
            //! Inviscid wall or symmetry boundary.
            dq[IR][re] = dq[IR][le];
            rvn = 2.0 * (xfn[iFace]*dq[IRU][le] + yfn[iFace]*dq[IRV][le] + zfn[iFace]*dq[IRW][le]);
            dq[IRU][re] = dq[IRU][le] - xfn[iFace] * rvn;
            dq[IRV][re] = dq[IRV][le] - yfn[iFace] * rvn;
            dq[IRW][re] = dq[IRW][le] - zfn[iFace] * rvn;
            rvn   = 2.0 * (xfn[iFace]*Q[IRU] + yfn[iFace]*Q[IRV] + zfn[iFace]*Q[IRW]);
            ru    = Q[IRU] - xfn[iFace] * rvn;
            rv    = Q[IRV] - yfn[iFace] * rvn;
            rw    = Q[IRW] - zfn[iFace] * rvn;
            p_new = Q[IRE] - 0.5 * (Q[IRU]*Q[IRU] + Q[IRV]*Q[IRV] + Q[IRW]*Q[IRW]) / Q[IR];

            gam1 = gama[le] - 1.0;
            dq[IRE][re] = p_new - pm/gam1 + half*((ru*ru + rv*rv + rw*rw)/Q[IR] - rm*(um*um + vm*vm + wm*wm));
        }
        else if (IsWall(bcType) && viscousType)
        {
            ub = 0.0;
            vb = 0.0;
            wb = 0.0;
            dq[IR ][re] =   dq[IR ][le];
            dq[IRU][re] = - dq[IRU][le] + two * ub * dq[IR][le];
            dq[IRV][re] = - dq[IRV][le] + two * vb * dq[IR][le];
            dq[IRW][re] = - dq[IRW][le] + two * wb * dq[IR][le];
            dq[IRE][re] =   dq[IRE][le] + ub * (dq[IRU][re] - dq[IRU][le])
                + vb * (dq[IRV][re] - dq[IRV][le])
                + wb * (dq[IRW][re] - dq[IRW][le]);
        }
        else if (bcType == PHENGLEI::EXTRAPOLATION)
        {
            for (int m = 0; m < nl; ++ m)
            {
                dq[m][re] =  dq[m][le];
            }
            for (int m = 0; m < ntmodel - 1; m++)
            {
                dq[nl + nchem + m][re] = dq[nl + nchem + m][le];
            }
        }
        else
        {
            //! FIX_ALL set DQ = 0.
            //! Far_FIELD, PROFILE and FIX_V have not been implemented yet.
            for (int m = 0; m < nl; ++ m)
            {
                dq[m][re] = 0.0;
            }
            for (int m = 0; m < ntmodel - 1; m++)
            {
                dq[nl + nchem + m][re] = 0.0;
            }
        }
    }

    delete []primitiveVariable;
    delete []Q;
}

void ComputeDiagonalMatrice(StructGrid *gridIn)
{
    int iSimplifyViscousTerm = GlobalDataBase::GetIntParaFromDB("iSimplifyViscousTerm");
    //! Compute the number of equations.
    int nEquation = GlobalDataBase::GetIntParaFromDB("nl");

    RDouble4D &inviscidSpectrumRadius = *reinterpret_cast<RDouble4D *> (gridIn->GetDataPtr("invSpectralRadius"));
    RDouble4D &viscousSpectrumRadius = *reinterpret_cast<RDouble4D *> (gridIn->GetDataPtr("visSpectralRadius"));
    RDouble4D &primitiveVariables = *reinterpret_cast<RDouble4D *> (gridIn->GetDataPtr("q"));
    RDouble4D &primitiveTemperatures = *reinterpret_cast<RDouble4D *> (gridIn->GetDataPtr("t"));
    RDouble3D &deltaTime = *reinterpret_cast<RDouble3D *> (gridIn->GetDataPtr("dt"));
    RDouble3D &specificHeatRatio = *reinterpret_cast<RDouble3D *> (gridIn->GetDataPtr("gama"));

    int viscousType = GlobalDataBase::GetIntParaFromDB("viscousType");
    RDouble3D &viscosityLaminar = *reinterpret_cast<RDouble3D *> (gridIn->GetDataPtr("visl"));
    RDouble3D &viscosityTurbulence = *reinterpret_cast<RDouble3D *> (gridIn->GetDataPtr("vist"));

    RDouble4D &gxfn = *(gridIn->GetFaceNormalX());
    RDouble4D &gyfn = *(gridIn->GetFaceNormalY());
    RDouble4D &gzfn = *(gridIn->GetFaceNormalZ());
    RDouble3D &gvol = *(gridIn->GetCellVolume());
    RDouble4D &gvgn = *(gridIn->GetFaceNormalVelocity());
    RDouble4D &garea = *(gridIn->GetFaceArea());

    Range I, J, K, M(0, nEquation - 1);
    RDouble5D &diagonalMatrice = *reinterpret_cast<RDouble5D *> (gridIn->GetDataPtr("diagonalMatrice"));
    RDouble2D *tempMatrix = new RDouble2D(M, M, fortranArray);
    RDouble2D &tempMatrix2D = *reinterpret_cast<RDouble2D *> (tempMatrix);

    RDouble *primitiveVars = new RDouble[nEquation]();
    RDouble *temperature = new RDouble[nEquation]();
    int nTemperatureModel = GlobalDataBase::GetIntParaFromDB("ntmodel");

    int ist, ied, jst, jed, kst, ked;
    gridIn->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    int nDim = GetDim();
    int iR = 0, jR = 0, kR = 0;
    RDouble gama = 0.0, laminarViscosity = 0.0, turbulentViscosity = 0.0, vol1 = 0.0, vol2 = 0.0;
    RDouble xfn1 = 0.0, yfn1 = 0.0, zfn1 = 0.0, area1 = 0.0, vgn1 = 0.0;
    RDouble xfn2 = 0.0, yfn2 = 0.0, zfn2 = 0.0, area2 = 0.0, vgn2 = 0.0;
    RDouble invSpectrumRadius = 0.0, visSpectrumRadius = 0.0;

    //! To Compute diagonal matrix.
    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                //! Obtain the primitive variables.
                for (int m = 0; m < nEquation; ++ m)
                {
                    primitiveVars[m] = primitiveVariables(i, j, k, m);
                }
                gama = specificHeatRatio(i, j, k);
                if (viscousType > INVISCID)
                {
                    laminarViscosity = viscosityLaminar(i, j, k);
                    turbulentViscosity = viscosityTurbulence(i, j, k);
                }
                else
                {
                    laminarViscosity = 0.0;
                    turbulentViscosity = 0.0;
                }
                for (int m = 0; m < nTemperatureModel; ++ m)
                {
                    temperature[m] = primitiveTemperatures(i, j, k, m);
                }

                // initialize.
                for (int m = 0; m < nEquation; ++ m)
                {
                    for (int n = 0; n < nEquation; ++ n)
                    {
                        diagonalMatrice(i, j, k, m, n) = 0.0;
                    }
                }

                //! Compute the scalar coefficient of the diagonal matrix.
                for (int m = 0; m < nEquation; ++ m)
                {
                    diagonalMatrice(i, j, k, m, m) = 1.0 / deltaTime(i, j, k);
                }

                invSpectrumRadius = 0.0;
                visSpectrumRadius = 0.0;
                //! Construct the diagonal matrix called D.
                for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
                {
                    int il1, jl1, kl1;
                    GetNsurfIndex(iSurface, il1, jl1, kl1);
                    iR = i + il1;
                    jR = j + jl1;
                    kR = k + kl1;

                    xfn1 = gxfn(i, j, k, iSurface);
                    yfn1 = gyfn(i, j, k, iSurface);
                    zfn1 = gzfn(i, j, k, iSurface);
                    vgn1 = gvgn(i, j, k, iSurface);
                    area1 = garea(i, j, k, iSurface);

                    xfn2 = gxfn(iR, jR, kR, iSurface);
                    yfn2 = gyfn(iR, jR, kR, iSurface);
                    zfn2 = gzfn(iR, jR, kR, iSurface);
                    vgn2 = gvgn(iR, jR, kR, iSurface);
                    area2 = garea(iR, jR, kR, iSurface);

                    //! Compute spectrum radius of inviscid flux.
                    if (iSimplifyViscousTerm == 1)// Equivalent method to computing the spectrum radius.
                    {
                        invSpectrumRadius += inviscidSpectrumRadius(i, j, k, iSurface - 1);
                        visSpectrumRadius += viscousSpectrumRadius(i, j, k, iSurface - 1);
                    }
                    else     //! Employ the new method.
                    {
                        //! Compute spectrum radius of inviscid flux.
                        invSpectrumRadius += 0.5 * ObtainInviscidSpectrumRadius(primitiveVars, xfn1, yfn1, zfn1, area1, vgn1, gama);
                        invSpectrumRadius += 0.5 * ObtainInviscidSpectrumRadius(primitiveVars, xfn2, yfn2, zfn2, area2, vgn2, gama);

                        //! Compute spectrum radius of viscous flux.
                        vol1 = 0.5 * (gvol(i, j, k) + gvol(i - il1 , j - jl1, k - kl1));
                        vol2 = 0.5 * (gvol(i, j, k) + gvol(iR , jR, kR));
                        visSpectrumRadius += 0.5 * ObtainViscousSpectrumRadius(primitiveVars, area1, vol1, gama, laminarViscosity, turbulentViscosity);
                        visSpectrumRadius += 0.5 * ObtainViscousSpectrumRadius(primitiveVars, area2, vol2, gama, laminarViscosity, turbulentViscosity);
                    }

                    //! To obtain the Jacobian matrix of the inviscid flux using the side of face -1/2.
                    ObtainInviscidMatrix(primitiveVars, temperature, xfn1, yfn1, zfn1, area1, vgn1, gama, tempMatrix);
                    for (int m = 0; m < nEquation; ++ m)
                    {
                        for (int n = 0; n < nEquation; ++ n)
                        {
                            diagonalMatrice(i, j, k, m, n) += -0.5 * tempMatrix2D(m, n);
                        }
                    }

                    //! To obtain the Jacobian matrix of the inviscid flux using the side of face +1/2.
                    ObtainInviscidMatrix(primitiveVars, temperature, xfn2, yfn2, zfn2, area2, vgn2, gama, tempMatrix);
                    for (int m = 0; m < nEquation; ++ m)
                    {
                        for (int n = 0; n < nEquation; ++ n)
                        {
                            diagonalMatrice(i, j, k, m, n) += 0.5 * tempMatrix2D(m, n);
                        }
                    }

                    if (viscousType > INVISCID && iSimplifyViscousTerm != 1)
                    {
                        //! To obtain the Jacobian matrix of the viscous flux using the side of face -1/2.
                        ObtainViscousMatrix(primitiveVars, temperature, xfn1, yfn1, zfn1, area1, vol1, vgn1, gama, laminarViscosity, turbulentViscosity, tempMatrix);
                        for (int m = 0; m < nEquation; ++ m)
                        {
                            for (int n = 0; n < nEquation; ++ n)
                            {
                                diagonalMatrice(i, j, k, m, n) += 0.5 * tempMatrix2D(m, n);
                            }
                        }

                        //! To obtain the Jacobian matrix of the viscous flux using the side of face +1/2.
                        ObtainViscousMatrix(primitiveVars, temperature, xfn2, yfn2, zfn2, area2, vol2, vgn2, gama, laminarViscosity, turbulentViscosity, tempMatrix);
                        for (int m = 0; m < nEquation; ++ m)
                        {
                            for (int n = 0; n < nEquation; ++ n)
                            {
                                diagonalMatrice(i, j, k, m, n) += -0.5 * tempMatrix2D(m, n);
                            }
                        }
                    }
                }

                //! Compute the summation of spectrum radius.
                for (int m = 0; m < nEquation; ++ m)
                {
                    diagonalMatrice(i, j, k, m, m) += invSpectrumRadius + visSpectrumRadius;
                }

            }
        }
    }

    delete []primitiveVars;
    delete []temperature;
    delete tempMatrix;
}

RDouble ObtainInviscidSpectrumRadius(RDouble *primitiveVars, RDouble xfn, RDouble yfn, RDouble zfn, RDouble area, RDouble vgn, RDouble gama)
{
    RDouble invSpectrumRadius = 0.0;
    //! Obtain the primitive variables.
    using namespace IDX;
    RDouble uVelocity = primitiveVars[IU];
    RDouble vVelocity = primitiveVars[IV];
    RDouble wVelocity = primitiveVars[IW];
    //! Compute the relative velocity.
    RDouble rVelocity = uVelocity * xfn + vVelocity * yfn + wVelocity * zfn - vgn;

    //! Compute the sound velocity.
    RDouble density = primitiveVars[IR];
    RDouble pressure = primitiveVars[IP];
    RDouble soundVelocity = sqrt(gama * pressure / density);

    //! Compute the spectrum radius of inviscid flux Jacobian matrix.
    invSpectrumRadius = (fabs(rVelocity) + soundVelocity) * area;

    return invSpectrumRadius;
}

RDouble ObtainViscousSpectrumRadius(RDouble *primitiveVars, RDouble area, RDouble vol, RDouble gama, RDouble laminarViscosity, RDouble turbulentViscosity)
{
    RDouble visSpectrumRadius = 0.0;
    //! Obtain the primitive variables.
    using namespace IDX;
    RDouble density = primitiveVars[IR];
    RDouble refReNumber = GlobalDataBase::GetDoubleParaFromDB("refReNumber");
    RDouble prandtlLaminar = GlobalDataBase::GetDoubleParaFromDB("prl");
    RDouble prandtlTurbulence = GlobalDataBase::GetDoubleParaFromDB("prt");

    //! Compute the spectrum radius of viscous flux Jacobian matrix.
    visSpectrumRadius = MAX((4.0 / 3.0) * (laminarViscosity + turbulentViscosity), gama * (laminarViscosity / prandtlLaminar + turbulentViscosity / prandtlTurbulence));
    visSpectrumRadius *= (area * area) / (density * vol * refReNumber);

    return visSpectrumRadius;
}

void ObtainInviscidMatrix(RDouble *primitiveVars, RDouble *temperature, RDouble xfn, RDouble yfn, RDouble zfn, RDouble area, RDouble vgn, 
                          RDouble gama, RDouble2D *tempMatrix)
{
    int nLaminar = GlobalDataBase::GetIntParaFromDB("nl");
    int nChemical = GlobalDataBase::GetIntParaFromDB("nchem");
    int nTemperatureModel = GlobalDataBase::GetIntParaFromDB("ntmodel");
    //! Compute the number of equations.
    int neqn = nLaminar + nChemical + nTemperatureModel - 1;

    //! To store the Jacobian matrix of inviscid flux.
    RDouble2D &inviscidMatrix = *reinterpret_cast<RDouble2D *> (tempMatrix);

    using namespace IDX;
    //! Obtain the primitive variables.
    RDouble uVelocity, vVelocity, wVelocity, density, pressure, averageVelocity, rVelocity;
    density = primitiveVars[IR];
    pressure = primitiveVars[IP];
    uVelocity = primitiveVars[IU];
    vVelocity = primitiveVars[IV];
    wVelocity = primitiveVars[IW];
    averageVelocity = uVelocity * xfn + vVelocity * yfn + wVelocity * zfn;
    //! The relative velocity.
    rVelocity = averageVelocity - vgn;

    //! trTemperature is the translation-rotation temperature, vTemperature and eTemperature denote the temperatures \n
    //! of the vibration mode and electron mode, respectively.
    RDouble trTemperature, vTemperature, eTemperature;

    //! To obtain the value of temperature.
    trTemperature = temperature[ITT];
    if (nTemperatureModel > 1)
    {
        vTemperature = temperature[ITV];
        eTemperature = vTemperature;
        if (nTemperatureModel > 2)
        {
            eTemperature = temperature[ITE];
        }
    }
    else
    {
        vTemperature = trTemperature;
        eTemperature = trTemperature;
    }

    //! To compute the total enthalpy.
    RDouble squareVelocity, totalEnthalpy;
    squareVelocity = uVelocity * uVelocity + vVelocity * vVelocity + wVelocity * wVelocity;
    if (nChemical == 0)
    {
        totalEnthalpy = (gama / (gama - 1)) * pressure / density;
    }
    else
    {
        totalEnthalpy = gas->GetMixedGasEnthalpy(primitiveVars, trTemperature, vTemperature, eTemperature);
    }
    totalEnthalpy += 0.5 * squareVelocity;

    //! Construct the Jacobian matrix.
    RDouble alpha = 0.5 * (gama - 1.0) * squareVelocity;
    inviscidMatrix(0, 0) = -vgn;
    inviscidMatrix(1, 0) = alpha * xfn - uVelocity * averageVelocity;
    inviscidMatrix(2, 0) = alpha * yfn - vVelocity * averageVelocity;
    inviscidMatrix(3, 0) = alpha * zfn - wVelocity * averageVelocity;
    inviscidMatrix(4, 0) = (alpha - totalEnthalpy) * averageVelocity;

    //! The second column.
    inviscidMatrix(0, 1) = xfn;
    inviscidMatrix(1, 1) = (2.0 - gama) * uVelocity * xfn + rVelocity;
    inviscidMatrix(2, 1) = (1.0 - gama) * uVelocity * yfn + vVelocity * xfn;
    inviscidMatrix(3, 1) = (1.0 - gama) * uVelocity * zfn + wVelocity * xfn;
    inviscidMatrix(4, 1) = (1.0 - gama) * uVelocity * averageVelocity + totalEnthalpy * xfn;

    //! The third column.
    inviscidMatrix(0, 2) = yfn;
    inviscidMatrix(1, 2) = (1.0 - gama) * vVelocity * xfn + uVelocity * yfn;
    inviscidMatrix(2, 2) = (2.0 - gama) * vVelocity * yfn + rVelocity;
    inviscidMatrix(3, 2) = (1.0 - gama) * vVelocity * zfn + wVelocity * yfn;
    inviscidMatrix(4, 2) = (1.0 - gama) * vVelocity * averageVelocity + totalEnthalpy * yfn;

    //! The fourth column.
    inviscidMatrix(0, 3) = zfn;
    inviscidMatrix(1, 3) = (1.0 - gama) * wVelocity * xfn + uVelocity * zfn;
    inviscidMatrix(2, 3) = (1.0 - gama) * wVelocity * yfn + vVelocity * zfn;
    inviscidMatrix(3, 3) = (2.0 - gama) * wVelocity * zfn + rVelocity;
    inviscidMatrix(4, 3) = (1.0 - gama) * wVelocity * averageVelocity + totalEnthalpy * zfn;

    //! The fifth column.
    inviscidMatrix(0, 4) = 0.0;
    inviscidMatrix(1, 4) = (gama - 1.0) * xfn;
    inviscidMatrix(2, 4) = (gama - 1.0) * yfn;
    inviscidMatrix(3, 4) = (gama - 1.0) * zfn;
    inviscidMatrix(4, 4) = (gama - 1.0) * averageVelocity + rVelocity;

    //! To multiply the average area located in the center of the cell.
    for (int m = 0; m < neqn; ++ m)
    {
        for (int n = 0; n < neqn; ++ n)
        {
            inviscidMatrix(m, n) *= area;
        }
    }

}

void ObtainViscousMatrix(RDouble *primitiveVars, RDouble *temperature, RDouble xfn, RDouble yfn, RDouble zfn, RDouble area, RDouble vol, 
                         RDouble vgn, RDouble gama, RDouble laminarViscosity, RDouble turbulentViscosity, RDouble2D *tempMatrix)
{
    int nLaminar = GlobalDataBase::GetIntParaFromDB("nl");
    int nChemical = GlobalDataBase::GetIntParaFromDB("nchem");
    int nTemperatureModel = GlobalDataBase::GetIntParaFromDB("ntmodel");
    RDouble refReNumber = GlobalDataBase::GetDoubleParaFromDB("refReNumber");
    //! Compute the number of equations.
    int neqn = nLaminar + nChemical + nTemperatureModel - 1;

    //! To store the Jacobian matrix of inviscid flux.
    RDouble2D &viscousMatrix = *reinterpret_cast<RDouble2D *> (tempMatrix);

    using namespace IDX;
    //! Obtain the primary variables.
    RDouble uVelocity, vVelocity, wVelocity, density, pressure, averageVelocity;
    density = primitiveVars[IR];
    pressure = primitiveVars[IP];
    uVelocity = primitiveVars[IU];
    vVelocity = primitiveVars[IV];
    wVelocity = primitiveVars[IW];
    averageVelocity = uVelocity * xfn + vVelocity * yfn + wVelocity * zfn;

    //! trTemperature is the translation-rotation temperature, vTemperature and eTemperature denote the temperatures \n
    //! of the vibration mode and electron mode, respectively.
    RDouble trTemperature, vTemperature, eTemperature;

    //! To obtain the value of temperature.
    trTemperature = temperature[ITT];
    if (nTemperatureModel > 1)
    {
        vTemperature = temperature[ITV];
        eTemperature = vTemperature;
        if (nTemperatureModel > 2)
        {
            eTemperature = temperature[ITE];
        }
    }
    else
    {
        vTemperature = trTemperature;
        eTemperature = trTemperature;
    }

    //! ktr is the heat conductivity of translation and rotation mode, kv and ke are heat conductivities \n
    //! of vibration mode and electron mode, respectively. kcp denotes the total heat conductivity of mixed gas.
    RDouble ktr, kv, ke, kcp;
    RDouble oprl = GlobalDataBase::GetDoubleParaFromDB("oprl");
    RDouble oprt = GlobalDataBase::GetDoubleParaFromDB("oprt");

    //! Compute the heat conductivity.
    if (nChemical == 0) //The heat conductivity of perfect gas
    {
        RDouble R = gas->GetUniversalGasConstant();
        ktr = (laminarViscosity * oprl + turbulentViscosity * oprt) * gama / (gama - 1.0) * R;
        kv = 0.0;
        ke = 0.0;
    }
    else
    {
        gas->ComputeMixtureGasHeatConductivity(primitiveVars, trTemperature, vTemperature, eTemperature, ktr, kv, ke);
    }
    //! To compute the total heat conductivity and viscosity.
    kcp = ktr + kv + ke;
    RDouble viscosity = laminarViscosity + turbulentViscosity;

    //! To compute the square velocity.
    RDouble squareVelocity;
    squareVelocity = uVelocity * uVelocity + vVelocity * vVelocity + wVelocity * wVelocity;

    //! Construct the Jacobian matrix.
    RDouble alpha = viscosity / density;
    RDouble phi = kcp * (gama - 1.0) * pressure / trTemperature;
    viscousMatrix(0, 0) = 0.0;
    viscousMatrix(1, 0) = -alpha * (uVelocity + 1.0/3.0 * averageVelocity * xfn);
    viscousMatrix(2, 0) = -alpha * (vVelocity + 1.0/3.0 * averageVelocity * yfn);
    viscousMatrix(3, 0) = -alpha * (wVelocity + 1.0/3.0 * averageVelocity * zfn);
    viscousMatrix(4, 0) = -kcp * trTemperature / density + 0.5 * phi * squareVelocity
        -2.0 * alpha * (squareVelocity + 1.0/3.0 * averageVelocity * averageVelocity);

    //! The second column.
    viscousMatrix(0, 1) = 0.0;
    viscousMatrix(1, 1) = alpha * (1.0 + 1.0/3.0 * xfn * xfn);
    viscousMatrix(2, 1) = alpha * (0.0 + 1.0/3.0 * xfn * yfn);
    viscousMatrix(3, 1) = alpha * (0.0 + 1.0/3.0 * xfn * zfn);
    viscousMatrix(4, 1) = 2.0 * alpha * (uVelocity + 1.0/3.0 * averageVelocity * xfn) - phi * uVelocity;

    // The third column.
    viscousMatrix(0, 2) = 0.0;
    viscousMatrix(1, 2) = alpha * (0.0 + 1.0/3.0 * xfn * yfn);
    viscousMatrix(2, 2) = alpha * (1.0 + 1.0/3.0 * yfn * yfn);
    viscousMatrix(3, 2) = alpha * (0.0 + 1.0/3.0 * yfn * zfn);
    viscousMatrix(4, 2) = 2.0 * alpha * (vVelocity + 1.0/3.0 * averageVelocity * yfn) - phi * vVelocity;

    //! The fourth column.
    viscousMatrix(0, 3) = 0.0;
    viscousMatrix(1, 3) = alpha * (0.0 + 1.0/3.0 * xfn * zfn);
    viscousMatrix(2, 3) = alpha * (0.0 + 1.0/3.0 * yfn * zfn);
    viscousMatrix(3, 3) = alpha * (1.0 + 1.0/3.0 * zfn * zfn);
    viscousMatrix(4, 3) = 2.0 * alpha * (wVelocity + 1.0/3.0 * averageVelocity * zfn) - phi * wVelocity;

    //! The fifth column.
    viscousMatrix(0, 4) = 0.0;
    viscousMatrix(1, 4) = 0.0;
    viscousMatrix(2, 4) = 0.0;
    viscousMatrix(3, 4) = 0.0;
    viscousMatrix(4, 4) = phi;

    //! To multiply the average area located in the center of the cell.
    for (int m = 0; m < neqn; ++ m)
    {
        for (int n = 0; n < neqn; ++ n)
        {
            viscousMatrix(m, n) *= area * area / (vol * refReNumber);
        }
    }
}

void FillField(FieldProxy *fieldTargetProxy, FieldProxy *fieldSourceProxy)
{
    RDouble4D &fieldTarget = fieldTargetProxy->GetField_STR();
    RDouble4D &fieldSource = fieldSourceProxy->GetField_STR();

    fieldTarget = fieldSource;
}

FieldProxy *GetFieldProxy(StructGrid *gridIn, const string &fieldName)
{
    RDouble4D *field = reinterpret_cast <RDouble4D *> (gridIn->GetDataPtr(fieldName));

    FieldProxy *field_proxy = new FieldProxy();

    field_proxy->SetField_STR(field);

    return field_proxy;
}

void LHSMatrixsInHouse(UnstructGrid *gridIn)
{
    InitLHSMatrix(gridIn);

    LHSDiagMatrix(gridIn);

    InviscidMatrixs(gridIn);

    ViscousMatrixs(gridIn);

    MAT3D< RDouble > *impmat = reinterpret_cast< MAT3D< RDouble > * >(gridIn->GetDataPtr("impmat"));
    int nTotalCell = gridIn->GetNTotalCell();
    int nm = 5;
    GlobalDataBase::GetData("nm", &nm, PHINT, 1);

    LUDCMPBlockDiag(impmat, nTotalCell, nm);
}

void InitLHSMatrix(UnstructGrid *gridIn)
{
    int nTotalCell = gridIn->GetNTotalCell();

    MAT3D< RDouble > *impmat = reinterpret_cast< MAT3D< RDouble > * >(gridIn->GetDataPtr("impmat"));

    int nm = GlobalDataBase::GetIntParaFromDB("nm");

    if (!impmat)
    {
        impmat = new MAT3D< RDouble >(nTotalCell,nm,nm);
        gridIn->UpdateDataPtr("impmat", impmat);
    }

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        SetValue(impmat->GetT2D(iCell), nm, nm, static_cast<RDouble>(0.0));
    }
}

void LHSDiagMatrix(UnstructGrid *gridIn)
{
    int nTotalCell = gridIn->GetNTotalCell();

    MAT3D< RDouble > *impmat = reinterpret_cast< MAT3D< RDouble > * >(gridIn->GetDataPtr("impmat"));

    RDouble *dt  = reinterpret_cast< RDouble  * > (gridIn->GetDataPtr("dt"));

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        RDouble **A = impmat->GetT2D(iCell);

        RDouble vdt  = 1.0 / dt[iCell];

        A[0][0] = vdt;
        A[1][1] = vdt;
        A[2][2] = vdt;
        A[3][3] = vdt;
        A[4][4] = vdt;
    }
}


void InviscidMatrixs(UnstructGrid *gridIn)
{
    int nTotalFace = gridIn->GetNTotalFace();

    int nl = GlobalDataBase::GetIntParaFromDB("nl");

    FaceProxy   *faceProxy = CreateFaceProxyInHouse(gridIn);

    MAT3D< RDouble > *dfdq  = new MAT3D< RDouble >(SEGCTION_LENGTH, nl, nl);

    int localStart,localEnd;
    localStart = 0;
    do
    {
        localEnd = localStart + SEGCTION_LENGTH;
        if (localEnd > nTotalFace) localEnd = nTotalFace;

        GetQlQrInHouse(gridIn, faceProxy, localStart, localEnd);

        /* The Jacobian of the cell itself effect little.
        GetJacobianFull(gridIn, faceProxy->GetQL(), dfdq, localStart, localEnd);
        AddImpLHSLeft(gridIn, dfdq, localStart, localEnd);

        GetJacobianFull(gridIn, faceProxy->GetQR(), dfdq, localStart, localEnd);
        AddImpLHSRight(gridIn, dfdq, localStart, localEnd);//*/

        RoeDispMatrix(gridIn, faceProxy, dfdq, localStart, localEnd);

        AddImpLHSLeft(gridIn, dfdq, localStart, localEnd);
        dfdq->SwitchSign();
        AddImpLHSRight(gridIn, dfdq, localStart, localEnd);

        localStart = localEnd;
    } while (localStart < nTotalFace);

    delete faceProxy;
    delete dfdq;
}

//! ViscousMats computes Viscous jacobi matrix dFv/dQ for all cells.
void ViscousMatrixs(UnstructGrid *gridIn)
{
    int nTotalFace = gridIn->GetNTotalFace();

    int viscousType = GlobalDataBase::GetIntParaFromDB("viscousType");
    if (viscousType == INVISCID) return;

    int nl = 0;
    GlobalDataBase::GetData("nl", &nl, PHINT, 1);

    int nchem;
    GlobalDataBase::GetData("nchem", &nchem, PHINT, 1);

    FaceProxy   *faceProxy = CreateFaceProxyInHouse(gridIn);

    MAT3D< RDouble > *dfdq  = new MAT3D< RDouble >(SEGCTION_LENGTH, nl, nl);

    int localStart,localEnd;
    localStart = 0;
    do
    {
        localEnd = localStart + SEGCTION_LENGTH;
        if (localEnd > nTotalFace) localEnd = nTotalFace;

        GetQlQrInHouse(gridIn, faceProxy, localStart, localEnd);

        ViscousMatrix(gridIn, faceProxy, dfdq, localStart, localEnd);

        AddImpLHSLeft(gridIn, dfdq, localStart, localEnd);

        dfdq->SwitchSign();
        AddImpLHSRight(gridIn, dfdq, localStart, localEnd);

        localStart = localEnd;
    } while (localStart < nTotalFace);

    delete faceProxy;
    delete dfdq;
}

//! This is Viscous jacobi matrix dFv/dQ for conservative variables.
void ViscousMatrix(UnstructGrid *gridIn, FaceProxy *faceProxy, MAT3D< RDouble > *dfdq, int localStart, int localEnd)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    RDouble * ns   = grid->GetFaceArea();
    RDouble * xfn  = grid->GetFaceNormalX();
    RDouble * yfn  = grid->GetFaceNormalY();
    RDouble * zfn  = grid->GetFaceNormalZ();

    int nBoundFace = grid->GetNBoundFace();

    int * leftCellofFace = grid->GetLeftCellOfFace();
    int * rightCellofFace = grid->GetRightCellOfFace();

    RDouble **qL   = faceProxy->GetQL();
    RDouble **qR   = faceProxy->GetQR();

    RDouble refGama = GlobalDataBase::GetDoubleParaFromDB("refGama");

    RDouble * xcc = grid->GetCellCenterX();
    RDouble * ycc = grid->GetCellCenterY();
    RDouble * zcc = grid->GetCellCenterZ();

    RDouble * visl = reinterpret_cast< RDouble * > (grid->GetDataPtr("visl"));
    RDouble * vist = reinterpret_cast< RDouble * > (grid->GetDataPtr("vist"));

    RDouble refReNumber = GlobalDataBase::GetDoubleParaFromDB("refReNumber");
    RDouble prl = GlobalDataBase::GetDoubleParaFromDB("prl");
    RDouble prt = GlobalDataBase::GetDoubleParaFromDB("prt");
    RDouble refMachNumber = GlobalDataBase::GetDoubleParaFromDB("refMachNumber");
    RDouble csrvist = 1.0;

    RDouble vis,vis_l,vis_t,cp,kcp;
    RDouble r,u,v,w,p,t,Vn, halfV2;
    RDouble y0,t0,t1, Dis, delta;
    RDouble T;
    RDouble gm2 = refGama * refMachNumber * refMachNumber;

    RDouble oprl = 1.0 / prl;
    RDouble oprt = 1.0 / prt;

    //! Notice cp here, because chemical reaction is not considered.

    cp = 1.0 / ((refGama - 1.0) * refMachNumber * refMachNumber);

    using namespace IDX;

    for (int iFace = localStart; iFace < localEnd; ++ iFace)
    {
        int j  = iFace - localStart;

        RDouble nx = xfn[iFace];
        RDouble ny = yfn[iFace];
        RDouble nz = zfn[iFace];
        RDouble S = ns[iFace];

        int le  = leftCellofFace [ iFace ];
        int re  = rightCellofFace[ iFace ];

        r = half * (qL[IR][j] + qR[IR][j]);
        u = half * (qL[IU][j] + qR[IU][j]);
        v = half * (qL[IV][j] + qR[IV][j]);
        w = half * (qL[IW][j] + qR[IW][j]);
        p = half * (qL[IP][j] + qR[IP][j]);

        r = ABS(r);
        p = ABS(p);

        t = gm2 * p / r;
        Vn = nx * u + ny * v + nz * w;
        halfV2 = half * (u * u + u * v + u * w);

        vis_l = half * (visl[le] + visl[re]);
        vis_t = half * (vist[le] + vist[re]);
        vis_t *= csrvist;

        kcp = (vis_l * oprl + vis_t * oprt) * cp;
        vis = vis_l + vis_t;

        Dis = ABS((xcc[le] - xcc[re]) * nx + (ycc[le] - ycc[re])*ny + (zcc[le] - zcc[re]) * nz);
        y0 = S / Dis / refReNumber;

        //! Take the mu and density(r) out.
        delta = vis * y0 /r;

        //! Temperature.
        T = refGama * refMachNumber * refMachNumber * p / r;

        t0 = delta * kcp * T / vis;
        t1 = delta * kcp * T * r * (refGama - 1) / p / vis;

        RDouble **A = dfdq->GetT2D(j);

        if (iFace < nBoundFace)
        {
            //! Degenerated to LU-SGS.
            RDouble viscousSpectrumRadius = delta;

            SetValue(A, 5, 5, 0.0);
            A[0][0] = viscousSpectrumRadius;
            A[1][1] = viscousSpectrumRadius;
            A[2][2] = viscousSpectrumRadius;
            A[3][3] = viscousSpectrumRadius;
            A[4][4] = viscousSpectrumRadius;
        }
        else
        {
            //! First row.
            A[0][0] = 0.0;
            A[0][1] = 0.0;
            A[0][2] = 0.0;
            A[0][3] = 0.0;
            A[0][4] = 0.0;

            //! Second row.
            A[1][0] = -delta * (u + third * Vn * nx);
            A[1][1] = delta * (third * nx * nx + 1.0);
            A[1][2] = delta * third * nx * ny;
            A[1][3] = delta * third * nx * nz;
            A[1][4] = 0;

            //! Third row.
            A[2][0] = -delta * (v + third * Vn * ny);
            A[2][1] = delta * third * ny * nx;
            A[2][2] = delta * (third * ny * ny + 1.0);
            A[2][3] = delta * third * ny * nz;
            A[2][4] = 0;

            //! Fourth row.
            A[3][0] = -delta * (w + third * Vn * nz);
            A[3][1] = delta * third * nz * nx;
            A[3][2] = delta * third * nz * ny;
            A[3][3] = delta * (third * nz * nz + 1.0);
            A[3][4] = 0;

            //! Fifth row.
            //! Count the effect of area and heat transfer here.
            //! Only when n and fn point to the same side of the face.
            A[4][0] = u * A[1][0] + v * A[2][0] + w * A[3][0]-t0 + t1 * halfV2;
            A[4][1] = u * A[1][1] + v * A[2][1] + w * A[3][1]-t1 * u;
            A[4][2] = u * A[1][2] + v * A[2][2] + w * A[3][2]-t1 * v;
            A[4][3] = u * A[1][3] + v * A[2][3] + w * A[3][3]-t1 * w;
            A[4][4] = t1;
        }
    }
}

//! This is roe matrix for conservative variables.
void RoeDispMatrix(UnstructGrid *gridIn, FaceProxy *faceProxy, MAT3D< RDouble > *dfdq, int localStart, int localEnd)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nBoundFace = grid->GetNBoundFace();

    RDouble *ns  = grid->GetFaceArea();
    RDouble *nxn = grid->GetFaceNormalX();
    RDouble *nyn = grid->GetFaceNormalY();
    RDouble *nzn = grid->GetFaceNormalZ();

    RDouble **qL = faceProxy->GetQL();
    RDouble **qR = faceProxy->GetQR();

    RDouble refGama = GlobalDataBase::GetDoubleParaFromDB("refGama");

    int RoeEntropyFixMethod = 3;
    RDouble RoeEntropyFixCoef1 = 0.0001, RoeEntropyFixCoef2 = 0.0001;
    GetEntropyFixMethodAndCoefficients(RoeEntropyFixMethod, RoeEntropyFixCoef1, RoeEntropyFixCoef2);

    RDouble t0 = refGama - 1.0;

    RDouble r,u,v,w,h,Vn,c2,c,halfV2,absVel;  //! Roe average variables.
    RDouble rl,ul,vl,wl,pl,hl;
    RDouble rr,ur,vr,wr,pr,hr;
    RDouble eigv1,eigv2,eigv3;
    RDouble y0,y1,y2;
    RDouble ratio, coef;
    using namespace IDX;

    for (int iFace = localStart; iFace < localEnd; ++ iFace)
    {
        int j  = iFace - localStart;

        RDouble nx = nxn[iFace];
        RDouble ny = nyn[iFace];
        RDouble nz = nzn[iFace];
        RDouble S = ns [iFace];

        rl = qL[IR][j];
        ul = qL[IU][j];
        vl = qL[IV][j];
        wl = qL[IW][j];
        pl = qL[IP][j];

        rr = qR[IR][j];
        ur = qR[IU][j];
        vr = qR[IV][j];
        wr = qR[IW][j];
        pr = qR[IP][j];

        if (rl < 1.0e-5) rl = 1.0e-5;
        if (rr < 1.0e-5) rr = 1.0e-5;

        //! Note that not to miss the gama.
        hl = refGama * pl / rl / t0 + half * (ul * ul + vl * vl + wl * wl);
        hr = refGama * pr / rr / t0 + half * (ur * ur + vr * vr + wr * wr);

        ratio = sqrt(rr / rl);
        coef  = 1.0 / (1.0 + ratio);

        r = sqrt(rl * rr);
        u = (ul + ur * ratio) * coef;
        v = (vl + vr * ratio) * coef;
        w = (wl + wr * ratio) * coef;
        h = (hl + hr * ratio) * coef;
        halfV2 = half * (u * u + v * v + w * w);
		absVel = sqrt(two * halfV2);
        Vn = nx * u + ny * v + nz * w;

        c2 = t0 * (h - halfV2);
        c  = sqrt(ABS(c2));

        eigv1 = ABS(Vn     );
        eigv2 = ABS(Vn + c);
        eigv3 = ABS(Vn - c);

        if (iFace < nBoundFace)
        {
            RDouble **A = dfdq->GetT2D(j);
            SetValue(A, 5, 5, 0.0);

            RDouble maxEig = MAX(eigv1, eigv2);
            maxEig = MAX(maxEig, eigv3);

            //! Diagonal elements.
            maxEig *= half * S;

            A[0][0] = maxEig;
            A[1][1] = maxEig;
            A[2][2] = maxEig;
            A[3][3] = maxEig;
            A[4][4] = maxEig;
        }
        else
        {
            //! Limit the eigen values.
            //! Entropy correction, to limit the magnitude of the three eigenvalue.
            //! Warning: This function calling may affect the efficient!
            Flux_RoeEntropyCorrection(eigv1, eigv2, eigv3, 0.0, 0.0, absVel, Vn, c, RoeEntropyFixMethod, RoeEntropyFixCoef1, RoeEntropyFixCoef2);

            eigv1    *= half * S;
            eigv2    *= half * S;
            eigv3    *= half * S;
            //

            y0 = eigv1;
            y1 = half*(2.0*eigv1-eigv2-eigv3)/c2;
            y2 = half*(eigv2-eigv3)/c;

            RDouble **A = dfdq->GetT2D(j);

            //! First row.
            A[0][0]=y0- halfV2*t0*y1-Vn*y2;
            A[0][1]=u*t0*y1+ nx*y2;
            A[0][2]=v*t0*y1+ ny*y2;
            A[0][3]=w*t0*y1+ nz*y2;
            A[0][4]=-t0*y1;

            //! Second row.
            A[1][0]=(nx*Vn*c2-u*halfV2*t0)*y1+(nx*halfV2*t0-u*Vn)*y2;
            A[1][1]=(u*u*t0-c2*nx*nx)*y1+(nx*u-nx*u*t0)*y2+y0;
            A[1][2]=(v*u*t0-c2*ny*nx)*y1+(ny*u-nx*v*t0)*y2;
            A[1][3]=(w*u*t0-c2*nz*nx)*y1+(nz*u-nx*w*t0)*y2;
            A[1][4]=-u*t0*y1+nx*t0*y2;

            //! Third row.
            A[2][0]=(ny*Vn*c2-v*halfV2*t0)*y1+(ny*halfV2*t0-v*Vn)*y2;
            A[2][1]=(u*v*t0-c2*nx*ny)*y1+(nx*v-ny*u*t0)*y2;
            A[2][2]=(v*v*t0-c2*ny*ny)*y1+(ny*v-ny*v*t0)*y2+y0;
            A[2][3]=(w*v*t0-c2*nz*ny)*y1+(nz*v-ny*w*t0)*y2;
            A[2][4]=-v*t0*y1+ny*t0*y2;

            //! Fourth row.
            A[3][0]=(nz*Vn*c2-w*halfV2*t0)*y1+(nz*halfV2*t0-w*Vn)*y2;
            A[3][1]=(u*w*t0-c2*nx*nz)*y1+(nx*w-nz*u*t0)*y2;
            A[3][2]=(v*w*t0-c2*ny*nz)*y1+(ny*w-nz*v*t0)*y2;
            A[3][3]=(w*w*t0-c2*nz*nz)*y1+(nz*w-nz*w*t0)*y2+y0;
            A[3][4]=-w*t0*y1+nz*t0*y2;

            //! Fifth row.
            A[4][0]=(-halfV2 *halfV2 *t0+Vn*Vn*c2-halfV2*c2)*y1+(halfV2*Vn*(t0-1.0)-c2*Vn/t0)*y2;
            A[4][1]=(u*c2+u*halfV2*t0-nx*Vn*c2)*y1+(nx*halfV2+nx*c2/t0-u*Vn*t0)*y2;
            A[4][2]=(v*c2+v*halfV2*t0-ny*Vn*c2)*y1+(ny*halfV2+ny*c2/t0-v*Vn*t0)*y2;
            A[4][3]=(w*c2+w*halfV2*t0-nz*Vn*c2)*y1+(nz*halfV2+nz*c2/t0-w*Vn*t0)*y2;
            A[4][4]=-halfV2*t0*y1+Vn*t0*y2+half*(eigv2+eigv3);
        }
    }
}

//! This is inviscid jacobi matrix  dF/dQ for conservative variables.
void GetJacobianFull(UnstructGrid *gridIn, RDouble **q, MAT3D< RDouble > *dfdq, int localStart, int localEnd)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    RDouble * area = grid->GetFaceArea();
    RDouble * xfn  = grid->GetFaceNormalX();
    RDouble * yfn  = grid->GetFaceNormalY();
    RDouble * zfn  = grid->GetFaceNormalZ();

    RDouble refGama = GlobalDataBase::GetDoubleParaFromDB("refGama");

    RDouble t0 = refGama - 1.0 ;

    using namespace IDX;

    for (int iFace = localStart; iFace < localEnd; ++ iFace)
    {
        int j  = iFace - localStart;

        //! Count the effects of area and 0.5 in the fluxes here.
        RDouble nx = half * xfn[iFace] * area[iFace];
        RDouble ny = half * yfn[iFace] * area[iFace];
        RDouble nz = half * zfn[iFace] * area[iFace];

        RDouble r = q[IR][j];
        RDouble u = q[IU][j];
        RDouble v = q[IV][j];
        RDouble w = q[IW][j];
        RDouble p = q[IP][j];

        RDouble halfV2 = half * (u * u + v * v + w * w);
        RDouble Vn = nx * u + ny * v + nz * w;
        RDouble H = refGama * p / r / t0 + halfV2;

        RDouble **A = dfdq->GetT2D(j);

        A[0][0] = 0;
        A[0][1] = nx;
        A[0][2] = ny;
        A[0][3] = nz;
        A[0][4] = 0.0;

        A[1][0] = t0 * halfV2 * nx - u * Vn;
        A[1][1] = (1-t0)* u * nx + Vn;
        A[1][2] = -t0 * v * nx + u * ny;
        A[1][3] = -t0 * w * nx + u * nz;
        A[1][4] = t0 * nx;

        A[2][0] = t0 * halfV2 * ny - v * Vn;
        A[2][1] = -t0 * u * ny + v * nx;
        A[2][2] = (1-t0)* v * ny + Vn;
        A[2][3] = -t0 * w * ny + v * nz;
        A[2][4] = t0 * ny;

        A[3][0] = t0 * halfV2 * nz - w * Vn;
        A[3][1] = -t0 * u * nz + w * nx;
        A[3][2] = -t0 * v * nz + w * ny;
        A[3][3] = (1-t0) * w * nz + Vn;
        A[3][4] = t0 * nz;

        A[4][0] = (t0 * halfV2 - H) * Vn;
        A[4][1] = H * nx - t0 * u * Vn;
        A[4][2] = H * ny - t0 * v * Vn;
        A[4][3] = H * nz - t0 * w * Vn;
        A[4][4] = refGama * Vn;
    }
}

void AddImpLHSLeft(UnstructGrid *gridIn, MAT3D< RDouble > *dfdq, int localStart, int localEnd)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int * leftCellofFace = grid->GetLeftCellOfFace();

    int nMid  = localStart;

    MAT3D< RDouble > *impmat = reinterpret_cast< MAT3D< RDouble > * >(grid->GetDataPtr("impmat"));

    int nm = GlobalDataBase::GetIntParaFromDB("nm");

    int le,j;
    RDouble **A, **B;
    for (int i = nMid; i < localEnd; ++ i)
    {
        le = leftCellofFace[ i ];
        j  = i - localStart;

        A = impmat->GetT2D(le);
        B = dfdq->GetT2D(j);

        for (int ii = 0; ii < nm; ++ ii)
        {
            for (int jj = 0; jj < nm; ++ jj)
            {
                A[ii][jj] += B[ii][jj];
            }
        }
    }
}

void AddImpLHSRight(UnstructGrid *gridIn, MAT3D< RDouble > *dfdq, int localStart, int localEnd)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nBoundFace = grid->GetNBoundFace();

    int * rightCellofFace = grid->GetRightCellOfFace();

    int nMid  = localStart;
    if (localEnd  <= nBoundFace)
    {
        nMid = localEnd;
    }
    else if (localStart < nBoundFace)
    {
        nMid = nBoundFace;
    }

    MAT3D< RDouble > *impmat = reinterpret_cast< MAT3D< RDouble > * >(grid->GetDataPtr("impmat"));

    int nm = 5;
    GlobalDataBase::GetData("nm", &nm, PHINT, 1);

    int re,j;
    RDouble **A, **B;
    for (int i = nMid; i < localEnd; ++ i)
    {
        re = rightCellofFace[ i ];
        j  = i - localStart;

        A = impmat->GetT2D(re);
        B = dfdq->GetT2D(j);
        for (int ii = 0; ii < nm; ++ ii)
        {
            for (int jj = 0; jj < nm; ++ jj)
            {
                A[ii][jj] -= B[ii][jj];
            }
        }
    }
}

FaceProxy * CreateFaceProxyInHouse(UnstructGrid *gridIn)
{
    int nl    = GlobalDataBase::GetIntParaFromDB("nl");
    int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
    int nEquation = nl + nchem;

    FaceProxy *faceProxy = new FaceProxy();
    faceProxy->Create(SEGCTION_LENGTH, nEquation);
    return faceProxy;
}

void GetQlQrInHouse(UnstructGrid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd)
{
    int *leftCellofFace = gridIn->GetLeftCellOfFace();
    int *rightCellofFace = gridIn->GetRightCellOfFace();

    int nl    = GlobalDataBase::GetIntParaFromDB("nl");
    int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
    int nEquation  = nl + nchem;

    RDouble **q = reinterpret_cast<RDouble **> (gridIn->GetDataPtr("q"));

    RDouble **qL = faceProxy->GetQL();
    RDouble **qR = faceProxy->GetQR();

    //! iFace is the global ID, jFace is the local ID in this section.
    for (int iFace = localStart; iFace < localEnd; ++ iFace)
    {
        int le = leftCellofFace [iFace];
        int re = rightCellofFace[iFace];
        int jFace  = iFace - localStart;

        for (int m = 0; m < nEquation; ++ m)
        {
            qL[m][jFace] = q[m][le];
            qR[m][jFace] = q[m][re];
        }
    }
}


}