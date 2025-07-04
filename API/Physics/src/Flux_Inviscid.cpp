#include "Flux_Inviscid.h"
#include "Constants.h"
#include "TK_Exit.h"
#include "Flux_RoeEntropyCorrection.h"
#include "PHMatrix.h"
#include "Glb_Dimension.h" 
using namespace std;

#include "Gas.h"

namespace PHSPACE
{
using namespace GAS_SPACE;

int GetSchemeID(const string &scheme_name)
{
    int scheme_id = ISCHEME_ROE;
    if (scheme_name.substr(0,3) == "roe")
    {
        scheme_id = ISCHEME_ROE;
    }
    else if (scheme_name.substr(0,7) == "ausmpw+")
    {
        scheme_id = ISCHEME_AUSMPW_PLUS;
    }
    else if (scheme_name.substr(0,6) == "ausmpw")
    {
        scheme_id = ISCHEME_AUSMPW;
    }
    else if (scheme_name.substr(0,6) == "ausm+w")
    {
        scheme_id = ISCHEME_AUSM_W;
    }
    else if (scheme_name.substr(0,7) == "ausm+up")
    {
        scheme_id = ISCHEME_AUSMPUP;
    }
    else if (scheme_name.substr(0,5) == "ausm+")
    {
        scheme_id = ISCHEME_AUSMP;
    }
    else if (scheme_name.substr(0,6) == "ausmdv")
    {
        scheme_id = ISCHEME_AUSMDV;
    }
    else if (scheme_name.substr(0,14) == "vanleer_Rotate")
    {
        scheme_id = ISCHEME_Rotate;
    }
    else if (scheme_name.substr(0,7) == "vanleer")
    {
        scheme_id = ISCHEME_VANLEER;
    }
    else if (scheme_name.substr(0,6) == "steger")
    {
        scheme_id = ISCHEME_STEGER;
    }
    else if (scheme_name.substr(0,4) == "hlle")
    {
        scheme_id = ISCHEME_HLLE;
    }
    else if (scheme_name.substr(0,5) == "lax_f")
    {
        scheme_id = ISCHEME_LAX_FRIEDRICHS;
    }
    //! GMRES
    else if (scheme_name.substr(0,8) == "GMRESRoe")
    {
        scheme_id = ISCHEME_GMRES_ROE;
    }
    else if (scheme_name.substr(0,11) == "GMRESSteger")  // GMRES3D
    {
        scheme_id = ISCHEME_GMRES_Steger;
    }
    else
    {
        TK_Exit::ExceptionExit("Error: this inv-scheme is not exist !\n", true);
    }

    return scheme_id;
}

INVScheme GetINVScheme(int scheme_id)
{
    INVScheme inv_scheme = nullptr;
    switch (scheme_id)
    {
        case ISCHEME_ROE:
            inv_scheme = Roe_Scheme_ConservativeForm;
            break;
        case ISCHEME_AUSMDV:
            inv_scheme = Ausmdv_Scheme;
            break;
        case ISCHEME_AUSMPW:
            inv_scheme = Ausmpw_Scheme;
            break;
        case ISCHEME_AUSMPW_PLUS:
            inv_scheme = AusmpwPlus_Scheme;
            break;
        case ISCHEME_VANLEER:
            inv_scheme = Vanleer_Scheme;
            break;
        case ISCHEME_Rotate:
            inv_scheme = Rotate_Vanleer_Scheme;
            break;
        case ISCHEME_STEGER:
            inv_scheme = Steger_Scheme;
            break;
        case ISCHEME_HLLE:
            inv_scheme = Hlle_Scheme;
            break;
        #ifdef USE_GMRESSOLVER
        case ISCHEME_GMRES_ROE:    //! GMRES
            inv_scheme = GMRES_INVScheme;    //! GMRES3D
            break;
        case ISCHEME_GMRES_Steger:    //! GMRES3D
            inv_scheme = GMRES_INVScheme;
            break;
        #endif
        default:
            TK_Exit::ExceptionExit("this scheme doesn't exist\n");
            break;
    }

    return inv_scheme;
}

void Roe_Scheme_ConservativeForm(FaceProxy *face_proxy, InviscidSchemeParameter *invSchemePara)
{
    const RDouble LITTLE = 1.0E-10;
    int nTemperatureModel = invSchemePara->GetNumberOfTemperatureModel();
    int nEquation = invSchemePara->GetNumberOfTotalEquation();
    int nLaminar = invSchemePara->GetNumberOfLaminar();
    int nLength = invSchemePara->GetLength();
    int nm = invSchemePara->Getnm();
    int RoeEntropyFixMethod = invSchemePara->GetEntropyFixMethod();
    int nChemical = invSchemePara->GetNumberOfChemical();
    int nPrecondition = invSchemePara->GetIfPrecondition();
    int nIdealState = invSchemePara->GetIfIdealState();
    int isUnsteady = invSchemePara->GetIfIsUnsteady();
    RDouble entrFixCoefLeft = invSchemePara->GetLeftEntropyFixCoefficient();
    RDouble entrFixCoefRight = invSchemePara->GetRightEntropyFixCoefficient();
    RDouble **qL = invSchemePara->GetLeftQ();
    RDouble **qR = invSchemePara->GetRightQ();
    RDouble *gamaL = invSchemePara->GetLeftGama();
    RDouble *gamaR = invSchemePara->GetRightGama();
    RDouble *pressureCoeffL = face_proxy->GetPressureCoefficientL();
    RDouble *pressureCoeffR = face_proxy->GetPressureCoefficientR();
    RDouble *xfn   = invSchemePara->GetFaceNormalX();
    RDouble *yfn   = invSchemePara->GetFaceNormalY();
    RDouble *zfn   = invSchemePara->GetFaceNormalZ();
    RDouble *vgn   = invSchemePara->GetFaceVelocity();
    RDouble **flux = invSchemePara->GetFlux();
    RDouble **tl = invSchemePara->GetLeftTemperature();
    RDouble **tr = invSchemePara->GetRightTemperature();

    RDouble *timeCoeffL = nullptr;
    RDouble *timeCoeffR = nullptr;
    RDouble *preconCoeffL = nullptr;
    RDouble *preconCoeffR = nullptr;
    if (nPrecondition == 1)
    {
        preconCoeffL = face_proxy->GetPreconCoefficientL();
        preconCoeffR = face_proxy->GetPreconCoefficientR();
        if (isUnsteady)
        {
        timeCoeffL = face_proxy->GetTimeCoefficientL();
        timeCoeffR = face_proxy->GetTimeCoefficientR();
        }
    }

    using namespace IDX;
    using namespace GAS_SPACE;

    RDouble *prim       = new RDouble[nEquation]();
    RDouble *qConserveL = new RDouble[nEquation]();
    RDouble *qConserveR = new RDouble[nEquation]();
    RDouble *dq         = new RDouble[nEquation]();

    RDouble *primL = new RDouble[nEquation]();
    RDouble *primR = new RDouble[nEquation]();
    RDouble temperaturesL[3] = {0}, temperaturesR[3] = {0};

    for (int iFace = 0; iFace < nLength; ++ iFace)
    {
        for (int iEqn = 0; iEqn < nEquation; ++ iEqn)
        {
            primL[iEqn] = qL[iEqn][iFace];
            primR[iEqn] = qR[iEqn][iFace];
        }

        RDouble &rL = primL[IR];
        RDouble &uL = primL[IU];
        RDouble &vL = primL[IV];
        RDouble &wL = primL[IW];
        RDouble &pL = primL[IP];

        RDouble &rR = primR[IR];
        RDouble &uR = primR[IU];
        RDouble &vR = primR[IV];
        RDouble &wR = primR[IW];
        RDouble &pR = primR[IP];

        RDouble xsn = xfn[iFace];
        RDouble ysn = yfn[iFace];
        RDouble zsn = zfn[iFace];
        RDouble vb  = vgn[iFace];
        RDouble gmL = gamaL[iFace];
        RDouble gmR = gamaR[iFace];

        rL = MAX(rL, LITTLE);
        rR = MAX(rR, LITTLE);
        pL = MAX(pL, LITTLE);
        pR = MAX(pR, LITTLE);

        RDouble vSquareL = uL * uL + vL * vL + wL * wL;
        RDouble vSquareR = uR * uR + vR * vR + wR * wR;

        RDouble hint_L = (gmL / (gmL - one)) * (pL / rL);
        RDouble hint_R = (gmR / (gmR - one)) * (pR / rR);
        RDouble peL = 0.0, peR = 0.0;

        if (nChemical != 0 && nIdealState == 0)
        {
            for (int i = 0; i < 3; ++ i)
            {
                temperaturesL[i] = tl[i][iFace];
                temperaturesR[i] = tr[i][iFace];
            }
            gas->ComputeEnthalpyByPrimitive(primL, gmL, hint_L, temperaturesL);
            gas->ComputeEnthalpyByPrimitive(primR, gmR, hint_R, temperaturesR);

            if (nTemperatureModel > 1)
            {
                peL = gas->GetElectronPressure(primL, temperaturesL[ITE]);
                peR = gas->GetElectronPressure(primR, temperaturesR[ITE]);
            }
        }
        RDouble hL = hint_L + half * vSquareL;
        RDouble hR = hint_R + half * vSquareR;

        RDouble vnL  = xsn * uL  + ysn * vL  + zsn * wL  - vb;
        RDouble vnR  = xsn * uR  + ysn * vR  + zsn * wR  - vb;
        RDouble rvnL = rL * vnL;
        RDouble rvnR = rR * vnR;

        //! The first part of Roe scheme: F(QL, QR) = 0.5 * [ F(QL) + F(QR) ].
        flux[IR ][iFace] =  half * (rvnL                        + rvnR                     );
        flux[IRU][iFace] =  half * (rvnL * uL + xfn[iFace] * pL + rvnR * uR + xfn[iFace] * pR);
        flux[IRV][iFace] =  half * (rvnL * vL + yfn[iFace] * pL + rvnR * vR + yfn[iFace] * pR);
        flux[IRW][iFace] =  half * (rvnL * wL + zfn[iFace] * pL + rvnR * wR + zfn[iFace] * pR);
        flux[IRE][iFace] =  half * (rvnL * hL + vgn[iFace] * pL + rvnR * hR + vgn[iFace] * pR);

        //! The following loop is useless.
        for (int iLaminar = nm; iLaminar < nLaminar; ++ iLaminar)
        {
            flux[iLaminar][iFace] = half * (primL[iLaminar] * rvnL + primR[iLaminar] * rvnR);
        }

        for (int it = 1; it < nTemperatureModel; ++ it)
        {
            flux[nLaminar + it][iFace] = half * (primL[nLaminar + it] * rvnL + primR[nLaminar + it] * rvnR);
        }
        if (nChemical == 1)
        {
            flux[nLaminar + nTemperatureModel - 1][iFace] += half * (peL * (vnL + vgn[iFace]) + peR * (vnR + vgn[iFace]));
        }

        RDouble ratio = sqrt(rR / rL);
        RDouble roeAverageCoef  = one / (one + ratio);

        RDouble rRoeAverage = sqrt(rL * rR);
        RDouble uRoeAverage = (uL + uR * ratio) * roeAverageCoef;
        RDouble vRoeAverage = (vL + vR * ratio) * roeAverageCoef;
        RDouble wRoeAverage = (wL + wR * ratio) * roeAverageCoef;
        RDouble pRoeAverage = (pL + pR * ratio) * roeAverageCoef;

        prim[0] = rRoeAverage;
        prim[1] = uRoeAverage;
        prim[2] = vRoeAverage;
        prim[3] = wRoeAverage;
        prim[4] = pRoeAverage;

        RDouble gama = (gmL + gmR * ratio) * roeAverageCoef;
        RDouble gamm1 = gama - one;

        RDouble vSquareRoeAverage = uRoeAverage * uRoeAverage + vRoeAverage * vRoeAverage + wRoeAverage * wRoeAverage;
        RDouble absVel = sqrt(vSquareRoeAverage);
        RDouble hRoeAverage = gama / gamm1 * pRoeAverage / rRoeAverage + half * vSquareRoeAverage;

        RDouble vn = xsn * uRoeAverage + ysn * vRoeAverage + zsn * wRoeAverage - vb;
        RDouble theta = vn + vb;

        //! c2: square of sound speed.
        RDouble c2 = gamm1 * (hRoeAverage - half * vSquareRoeAverage);

        //! Sound speed: gama * pRoeAverage / rRoeAverage.
        RDouble cm = sqrt(ABS(c2));

        RDouble vntmp = 0.0, cnew2 = 0.0, cnew = 0.0, vnew = 0.0, timeCoeff = 0.0, preconCoeff = 0.0, alpha = 0.0;
        RDouble eigv1 = 0.0, eigv2 = 0.0, eigv3 = 0.0;
        if (nPrecondition == 0)
        {
            eigv1 = ABS(vn);
            eigv2 = ABS(vn + cm);
            eigv3 = ABS(vn - cm);
        }
        else
        {
            preconCoeff = (preconCoeffL[iFace] + preconCoeffR[iFace] * ratio) * roeAverageCoef;
            vntmp = vn;
            vnew  = half * vn * (one + preconCoeff);
            cnew  = half * sqrt(((preconCoeff - one) * vn) * ((preconCoeff - one) * vn) + four * preconCoeff * c2);

            if (isUnsteady)
            {
            timeCoeff = (timeCoeffL[iFace] + timeCoeffR[iFace] * ratio) * roeAverageCoef;
                vntmp *= timeCoeff;
                vnew  *= timeCoeff;
                cnew  *= timeCoeff;
            }

            eigv1 = ABS(vntmp);
            eigv2 = ABS(vnew + cnew);
            eigv3 = ABS(vnew - cnew);
        }

        RDouble EL = (one / (gmL - one)) * (pL / rL) + half * vSquareL;
        RDouble ER = (one / (gmR - one)) * (pR / rR) + half * vSquareR;

        //! The jump between left and right conservative variables.
        dq[IR]  = rR - rL;
        dq[IRU] = rR * uR - rL * uL;
        dq[IRV] = rR * vR - rL * vL;
        dq[IRW] = rR * wR - rL * wL;
        dq[IRE] = rR * ER - rL * EL; 

        if (nChemical != 0)
        {
            if (nTemperatureModel > 1)
            {
                gas->GetTemperature(primL, temperaturesL[ITT], temperaturesL[ITV], temperaturesL[ITE]);
                gas->GetTemperature(primR, temperaturesR[ITT], temperaturesR[ITV], temperaturesR[ITE]);
            }
            gas->Primitive2Conservative(primL, gmL, temperaturesL[ITV], temperaturesL[ITE], qConserveL);
            gas->Primitive2Conservative(primR, gmR, temperaturesR[ITV], temperaturesR[ITE], qConserveR);
        }

        for (int iLaminar = nm; iLaminar < nEquation; ++ iLaminar)
        {
            dq[iLaminar] = qConserveR[iLaminar] - qConserveL[iLaminar];
        }

        //! Entropy correction, to limit the magnitude of the three eigenvalue.
        //! Warning: This function calling may affect the efficient!
        if(nPrecondition == 0)
        {
            Flux_RoeEntropyCorrection(eigv1, eigv2, eigv3, pressureCoeffL[iFace], pressureCoeffR[iFace], absVel, vn, cm, RoeEntropyFixMethod, entrFixCoefLeft, entrFixCoefRight);
        }
        else
        {
            Flux_RoeEntropyCorrection(eigv1, eigv2, eigv3, pressureCoeffL[iFace], pressureCoeffR[iFace], absVel, vnew, cnew, RoeEntropyFixMethod, entrFixCoefLeft, entrFixCoefRight);
        }
        
        RDouble xi1 = 0.0, xi2 = 0.0, xi3 = 0.0, xi4 = 0.0;
        if (nPrecondition == 0) 
        {
            xi1 = (two * eigv1 - eigv2 - eigv3) / (two * c2);
            xi2 = (eigv2 - eigv3)/(two * cm);
        }
        else
        {
            xi1 = (two * eigv1 - eigv2 - eigv3) / (two * c2);
            xi2 = (eigv2 - eigv3)/(two * cnew);
            alpha = half * (one - preconCoeff);
            xi3 = two * alpha * eigv1 / c2;
            xi4 = alpha * vn * xi2 / c2;
        }

        RDouble dc = theta * dq[IR] - xsn * dq[IRU] - ysn * dq[IRV] - zsn * dq[IRW];

        RDouble c2dc = 0;

        c2dc = c2 * dc;

        RDouble dh = - gamm1 * (uRoeAverage * dq[IRU] + vRoeAverage * dq[IRV] + wRoeAverage * dq[IRW] - dq[IRE]) + half * gamm1 * vSquareRoeAverage * dq[IR];

        if (nChemical != 0)
        {
            gas->ComputeDHAndTotalEnthalpy(prim, gama, dq, dh, hRoeAverage);
        }

        if (nPrecondition == 0)
        {
            flux[IR ][iFace] -= half * (eigv1 * dq[IR ]                 -               dh  * xi1               -               dc  * xi2);
            flux[IRU][iFace] -= half * (eigv1 * dq[IRU] + (xsn   * c2dc - uRoeAverage * dh) * xi1 + (xsn   * dh - uRoeAverage * dc) * xi2);
            flux[IRV][iFace] -= half * (eigv1 * dq[IRV] + (ysn   * c2dc - vRoeAverage * dh) * xi1 + (ysn   * dh - vRoeAverage * dc) * xi2);
            flux[IRW][iFace] -= half * (eigv1 * dq[IRW] + (zsn   * c2dc - wRoeAverage * dh) * xi1 + (zsn   * dh - wRoeAverage * dc) * xi2);
            flux[IRE][iFace] -= half * (eigv1 * dq[IRE] + (theta * c2dc - hRoeAverage * dh) * xi1 + (theta * dh - hRoeAverage * dc) * xi2);

            for (int iLaminar = nm; iLaminar < nLaminar; ++ iLaminar)
            {
                prim[iLaminar] = half * (primL[iLaminar] + primR[iLaminar]);
                flux[iLaminar][iFace] -= half * (eigv1 * dq[iLaminar] - prim[iLaminar] * (dh * xi1 + dc * xi2));
            }
        }
        else
        {
            flux[IR ][iFace] -= half * (eigv1 * dq[IR ]               - dh / preconCoeff * (xi1 - xi3 + xi4)                                            -               dc  * xi2);
            flux[IRU][iFace] -= half * (eigv1 * dq[IRU] - uRoeAverage * dh / preconCoeff * (xi1 - xi3 + xi4) + xsn   * c2dc * (xi1 - xi4) + (xsn   * dh - uRoeAverage * dc) * xi2);
            flux[IRV][iFace] -= half * (eigv1 * dq[IRV] - vRoeAverage * dh / preconCoeff * (xi1 - xi3 + xi4) + ysn   * c2dc * (xi1 - xi4) + (ysn   * dh - vRoeAverage * dc) * xi2);
            flux[IRW][iFace] -= half * (eigv1 * dq[IRW] - wRoeAverage * dh / preconCoeff * (xi1 - xi3 + xi4) + zsn   * c2dc * (xi1 - xi4) + (zsn   * dh - wRoeAverage * dc) * xi2);
            flux[IRE][iFace] -= half * (eigv1 * dq[IRE] - hRoeAverage * dh / preconCoeff * (xi1 - xi3 + xi4) + theta * c2dc * (xi1 - xi4) + (theta * dh - hRoeAverage * dc) * xi2);

            for (int iLaminar = nm; iLaminar < nLaminar; ++ iLaminar)
            {
                prim[iLaminar] = half * (primL[iLaminar] + primR[iLaminar]);
                flux[iLaminar][iFace] -= half * (eigv1 * dq[iLaminar] - prim[iLaminar] * (dh / preconCoeff * (xi1 - xi3 + xi4) + dc * xi2));
            }
            //! Compute vibration-electron energy term.
            for (int it = 1; it < nTemperatureModel; ++ it)
            {
                prim[nLaminar + it] = half * (primL[nLaminar + it] + primR[nLaminar + it]);
                flux[nLaminar + it][iFace] -= half * (eigv1 * dq[nLaminar + it] - prim[nLaminar + it] * (dh / preconCoeff * (xi1 - xi3 + xi4) + dc * xi2));
            }
        }

    }

    delete [] prim;    prim = nullptr;
    delete [] primL;    primL = nullptr;
    delete [] primR;    primR = nullptr;
    delete [] qConserveL;    qConserveL = nullptr;
    delete [] qConserveR;    qConserveR = nullptr;
    delete [] dq;    dq = nullptr;
}
#ifdef USE_GMRESSOLVER
//! GMRESCSR
inline int GetIndexOfBlockJacobianMatrix(vector<int> &AI, vector<int> &AJ, int nEquation, int rCell, int qCell)
{
    int index;
    vector<int>::iterator result = find(AJ.begin()+AI[rCell], AJ.begin()+AI[rCell + 1], qCell);    //! find qcell
    if(result==AJ.begin()+AI[rCell+1])
    {
        printf("\nrcell %d cannot find qcell %d from:\n", rCell, qCell);
        for (int i = AI[rCell]; i < AI[rCell+1]; i++)
        {
            printf("%d  ", AJ[i]);
        }
        printf("\n");

        TK_Exit::ExceptionExit("Sparse matrix index is wrong");
    }
    index = distance(AJ.begin(), result);
    index *= nEquation;
    return index;
}
//! GMRES
void GMRES_Roe_Scheme_ConservativeForm_test(FaceProxy *face_proxy, InviscidSchemeParameter *invSchemePara)
{
    const RDouble LITTLE = 1.0E-10;
    int nTemperatureModel = invSchemePara->GetNumberOfTemperatureModel();
    int nEquation = invSchemePara->GetNumberOfTotalEquation();
    int nLaminar = invSchemePara->GetNumberOfLaminar();
    int nLength = invSchemePara->GetLength();
    int nm = invSchemePara->Getnm();
    int RoeEntropyFixMethod = invSchemePara->GetEntropyFixMethod();
    int nChemical = invSchemePara->GetNumberOfChemical();
    int nPrecondition = invSchemePara->GetIfPrecondition();
    int nIdealState = invSchemePara->GetIfIdealState();

    RDouble refMachNumber = invSchemePara->GetMachNumber();
    RDouble entrFixCoefLeft = invSchemePara->GetLeftEntropyFixCoefficient();
    RDouble entrFixCoefRight = invSchemePara->GetRightEntropyFixCoefficient();
    RDouble **qL = invSchemePara->GetLeftQ();
    RDouble **qR = invSchemePara->GetRightQ();
    RDouble *gamaL = invSchemePara->GetLeftGama();
    RDouble *gamaR = invSchemePara->GetRightGama();
    RDouble *pressureCoeffL = face_proxy->GetPressureCoefficientL();
    RDouble *pressureCoeffR = face_proxy->GetPressureCoefficientR();
    RDouble *xfn = invSchemePara->GetFaceNormalX();
    RDouble *yfn = invSchemePara->GetFaceNormalY();
    RDouble *zfn = invSchemePara->GetFaceNormalZ();
    RDouble *vgn = invSchemePara->GetFaceVelocity();
    RDouble **flux = invSchemePara->GetFlux();
    RDouble **tl = invSchemePara->GetLeftTemperature();
    RDouble **tr = invSchemePara->GetRightTemperature();
    RDouble mach2 = refMachNumber * refMachNumber;

    using namespace IDX;
    using namespace GAS_SPACE;

    RDouble *prim = new RDouble[nEquation]();
    RDouble *qConserveL = new RDouble[nEquation]();
    RDouble *qConserveR = new RDouble[nEquation]();
    RDouble *dq = new RDouble[nEquation]();

    RDouble *primL = new RDouble[nEquation]();
    RDouble *primR = new RDouble[nEquation]();
    RDouble temperaturesL[3] = { 0 }, temperaturesR[3] = { 0 };

    for (int iFace = 0; iFace < nLength; ++iFace)
    {
        for (int iEqn = 0; iEqn < nEquation; ++iEqn)
        {
            primL[iEqn] = qL[iEqn][iFace];
            primR[iEqn] = qR[iEqn][iFace];
        }

        RDouble &rL = primL[IR];
        RDouble &uL = primL[IU];
        RDouble &vL = primL[IV];
        RDouble &wL = primL[IW];
        RDouble &pL = primL[IP];

        RDouble &rR = primR[IR];
        RDouble &uR = primR[IU];
        RDouble &vR = primR[IV];
        RDouble &wR = primR[IW];
        RDouble &pR = primR[IP];

        RDouble xsn = xfn[iFace];
        RDouble ysn = yfn[iFace];
        RDouble zsn = zfn[iFace];
        RDouble vb = vgn[iFace];
        RDouble gmL = gamaL[iFace];
        RDouble gmR = gamaR[iFace];

        rL = MAX(rL, LITTLE);
        rR = MAX(rR, LITTLE);
        pL = MAX(pL, LITTLE);
        pR = MAX(pR, LITTLE);

        RDouble vSquareL = uL * uL + vL * vL + wL * wL;
        RDouble vSquareR = uR * uR + vR * vR + wR * wR;

        RDouble hint_L = (gmL / (gmL - one)) * (pL / rL);
        RDouble hint_R = (gmR / (gmR - one)) * (pR / rR);
        RDouble peL = 0.0, peR = 0.0;

        if (nChemical != 0 && nIdealState == 0)
        {
            for (int i = 0; i < 3; ++i)
            {
                temperaturesL[i] = tl[i][iFace];
                temperaturesR[i] = tr[i][iFace];
            }
            gas->ComputeEnthalpyByPrimitive(primL, gmL, hint_L, temperaturesL);
            gas->ComputeEnthalpyByPrimitive(primR, gmR, hint_R, temperaturesR);

            if (nTemperatureModel > 1)
            {
                peL = gas->GetElectronPressure(primL, temperaturesL[ITE]);
                peR = gas->GetElectronPressure(primR, temperaturesR[ITE]);
            }
        }
        RDouble hL = hint_L + half * vSquareL;
        RDouble hR = hint_R + half * vSquareR;

        RDouble vnL = xsn * uL + ysn * vL + zsn * wL - vb;
        RDouble vnR = xsn * uR + ysn * vR + zsn * wR - vb;
        RDouble rvnL = rL * vnL;
        RDouble rvnR = rR * vnR;

        //! The first part of Roe scheme: F(QL, QR) = 0.5 * [ F(QL) + F(QR) ].
        flux[IR][iFace] = half * (rvnL + rvnR);
        flux[IRU][iFace] = half * (rvnL * uL + xfn[iFace] * pL + rvnR * uR + xfn[iFace] * pR);
        flux[IRV][iFace] = half * (rvnL * vL + yfn[iFace] * pL + rvnR * vR + yfn[iFace] * pR);
        flux[IRW][iFace] = half * (rvnL * wL + zfn[iFace] * pL + rvnR * wR + zfn[iFace] * pR);
        flux[IRE][iFace] = half * (rvnL * hL + vgn[iFace] * pL + rvnR * hR + vgn[iFace] * pR);

        //! The following loop is useless.
        for (int iLaminar = nm; iLaminar < nLaminar; ++iLaminar)
        {
            flux[iLaminar][iFace] = half * (primL[iLaminar] * rvnL + primR[iLaminar] * rvnR);
        }

        for (int it = 1; it < nTemperatureModel; ++it)
        {
            flux[nLaminar + it][iFace] = half * (primL[nLaminar + it] * rvnL + primR[nLaminar + it] * rvnR);
        }
        if (nChemical == 1)
        {
            flux[nLaminar + nTemperatureModel - 1][iFace] += half * (peL * (vnL + vgn[iFace]) + peR * (vnR + vgn[iFace]));
        }

        RDouble ratio = sqrt(rR / rL);
        RDouble roeAverageCoef = one / (one + ratio);

        RDouble rRoeAverage = sqrt(rL * rR);
        RDouble uRoeAverage = (uL + uR * ratio) * roeAverageCoef;
        RDouble vRoeAverage = (vL + vR * ratio) * roeAverageCoef;
        RDouble wRoeAverage = (wL + wR * ratio) * roeAverageCoef;
        RDouble pRoeAverage = (pL + pR * ratio) * roeAverageCoef;

        RDouble gama = (gmL + gmR * ratio) * roeAverageCoef;
        RDouble gamm1 = gama - one;

        RDouble vSquareRoeAverage = uRoeAverage * uRoeAverage + vRoeAverage * vRoeAverage + wRoeAverage * wRoeAverage;
        RDouble absVel = sqrt(vSquareRoeAverage);
        RDouble hRoeAverage = gama / gamm1 * pRoeAverage / rRoeAverage + half * vSquareRoeAverage;

        RDouble vn = xsn * uRoeAverage + ysn * vRoeAverage + zsn * wRoeAverage - vb;
        RDouble theta = vn + vb;

        //! c2: square of sound speed.
        RDouble c2 = gamm1 * (hRoeAverage - half * vSquareRoeAverage);

        //! Sound speed: gama * pRoeAverage / rRoeAverage.
        RDouble cm = sqrt(ABS(c2));

        RDouble Ur2 = 0, cnew2 = 0, cnew = 0, vnew = 0, beta = 0, alpha = 0;
        RDouble eigv1 = 0, eigv2 = 0, eigv3 = 0;
        eigv1 = ABS(vn);
        if (nPrecondition == 0)
        {
            eigv2 = ABS(vn + cm);
            eigv3 = ABS(vn - cm);
        }
        else
        {
            beta = MIN(MAX(vSquareRoeAverage / c2, three * mach2), one);
            vnew = half * vn * (one + beta);
            cnew = half * sqrt(((beta - one) * vn) * ((beta - one) * vn) + four * beta * c2);
            cnew2 = cnew * cnew;
            eigv2 = ABS(vnew + cnew);
            eigv3 = ABS(vnew - cnew);
        }

        RDouble EL = (one / (gmL - one)) * (pL / rL) + half * vSquareL;
        RDouble ER = (one / (gmR - one)) * (pR / rR) + half * vSquareR;

        //! The jump between left and right conservative variables.
        dq[IR] = rR - rL;
        dq[IRU] = rR * uR - rL * uL;
        dq[IRV] = rR * vR - rL * vL;
        dq[IRW] = rR * wR - rL * wL;
        dq[IRE] = rR * ER - rL * EL;

        if (nChemical != 0)
        {
            if (nTemperatureModel > 1)
            {
                gas->GetTemperature(primL, temperaturesL[ITT], temperaturesL[ITV], temperaturesL[ITE]);
                gas->GetTemperature(primR, temperaturesR[ITT], temperaturesR[ITV], temperaturesR[ITE]);
            }
            gas->Primitive2Conservative(primL, gmL, temperaturesL[ITV], temperaturesL[ITE], qConserveL);
            gas->Primitive2Conservative(primR, gmR, temperaturesR[ITV], temperaturesR[ITE], qConserveR);
        }

        for (int iLaminar = nm; iLaminar < nEquation; ++iLaminar)
        {
            dq[iLaminar] = qConserveR[iLaminar] - qConserveL[iLaminar];
        }

        //! Entropy correction, to limit the magnitude of the three eigenvalue.
        //! Warning: This function calling may affect the efficient!

        Flux_RoeEntropyCorrection(eigv1, eigv2, eigv3, pressureCoeffL[iFace], pressureCoeffR[iFace], absVel, vn, cm, RoeEntropyFixMethod, entrFixCoefLeft, entrFixCoefRight);

        RDouble xi1 = 0.0, xi2 = 0.0, xi3 = 0.0, xi4 = 0.0;
        if (nPrecondition == 0)
        {
            xi1 = (two * eigv1 - eigv2 - eigv3) / (two * c2);
            xi2 = (eigv2 - eigv3) / (two * cm);
        }
        else
        {
            xi1 = (two * eigv1 - eigv2 - eigv3) / (two * cnew2);
            xi2 = (eigv2 - eigv3) / (two * cnew);
            alpha = half * (one - beta);
            xi3 = two * alpha * eigv1 / cnew2;
            xi4 = alpha * vn * xi2 / cnew2;
        }

        RDouble dc = theta * dq[IR] - xsn * dq[IRU] - ysn * dq[IRV] - zsn * dq[IRW];

        RDouble c2dc = 0;
        if (nPrecondition == 0)
        {
            c2dc = c2 * dc;
        }
        else
        {
            Ur2 = beta * c2;
            c2dc = cnew2 * dc;
        }

        RDouble dh = -gamm1 * (uRoeAverage * dq[IRU] + vRoeAverage * dq[IRV] + wRoeAverage * dq[IRW] - dq[IRE]) + half * gamm1 * vSquareRoeAverage * dq[IR];

        if (nChemical != 0)
        {
            gas->ComputeDHAndTotalEnthalpy(prim, gama, dq, dh, hRoeAverage);
        }

        if (nPrecondition == 0)
        {
            flux[IR][iFace] -= half * (eigv1 * dq[IR] - dh * xi1 - dc * xi2);
            flux[IRU][iFace] -= half * (eigv1 * dq[IRU] + (xsn * c2dc - uRoeAverage * dh) * xi1 + (xsn * dh - uRoeAverage * dc) * xi2);
            flux[IRV][iFace] -= half * (eigv1 * dq[IRV] + (ysn * c2dc - vRoeAverage * dh) * xi1 + (ysn * dh - vRoeAverage * dc) * xi2);
            flux[IRW][iFace] -= half * (eigv1 * dq[IRW] + (zsn * c2dc - wRoeAverage * dh) * xi1 + (zsn * dh - wRoeAverage * dc) * xi2);
            flux[IRE][iFace] -= half * (eigv1 * dq[IRE] + (theta * c2dc - hRoeAverage * dh) * xi1 + (theta * dh - hRoeAverage * dc) * xi2);

            for (int iLaminar = nm; iLaminar < nLaminar; ++iLaminar)
            {
                prim[iLaminar] = half * (primL[iLaminar] + primR[iLaminar]);
                flux[iLaminar][iFace] -= half * (eigv1 * dq[iLaminar] - prim[iLaminar] * (dh * xi1 + dc * xi2));
            }
        }
        else
        {
            flux[IR][iFace] -= half * (eigv1 * dq[IR] - dh * cnew2 / Ur2 * xi1 - dc * xi2 + dh * cnew2 / Ur2 * xi3 - dh * cnew2 / Ur2 * xi4);
            flux[IRU][iFace] -= half * (eigv1 * dq[IRU] + (xsn * c2dc - uRoeAverage * dh * cnew2 / Ur2) * xi1 + (xsn * dh - uRoeAverage * dc) * xi2 + uRoeAverage * dh * cnew2 / Ur2 * xi3 - (xsn * c2dc + uRoeAverage * dh * cnew2 / Ur2) * xi4);
            flux[IRV][iFace] -= half * (eigv1 * dq[IRV] + (ysn * c2dc - vRoeAverage * dh * cnew2 / Ur2) * xi1 + (ysn * dh - vRoeAverage * dc) * xi2 + vRoeAverage * dh * cnew2 / Ur2 * xi3 - (ysn * c2dc + vRoeAverage * dh * cnew2 / Ur2) * xi4);
            flux[IRW][iFace] -= half * (eigv1 * dq[IRW] + (zsn * c2dc - wRoeAverage * dh * cnew2 / Ur2) * xi1 + (zsn * dh - wRoeAverage * dc) * xi2 + wRoeAverage * dh * cnew2 / Ur2 * xi3 - (zsn * c2dc + wRoeAverage * dh * cnew2 / Ur2) * xi4);
            flux[IRE][iFace] -= half * (eigv1 * dq[IRE] + (theta * c2dc - hRoeAverage * dh * cnew2 / Ur2) * xi1 + (theta * dh - hRoeAverage * dc) * xi2 + hRoeAverage * dh * cnew2 / Ur2 * xi3 - (theta * c2dc + hRoeAverage * dh * cnew2 / Ur2) * xi4);

            for (int iLaminar = nm; iLaminar < nLaminar; ++iLaminar)
            {
                prim[iLaminar] = half * (primL[iLaminar] + primR[iLaminar]);
                flux[iLaminar][iFace] -= half * (eigv1 * dq[iLaminar] - prim[iLaminar] * (dh * cnew2 / Ur2 * (xi1 - xi3 + xi4) + dc * xi2));
            }
            //! Compute vibration-electron energy term.
            for (int it = 1; it < nTemperatureModel; ++it)
            {
                prim[nLaminar + it] = half * (primL[nLaminar + it] + primR[nLaminar + it]);
                flux[nLaminar + it][iFace] -= half * (eigv1 * dq[nLaminar + it] - prim[nLaminar + it] * (dh * cnew2 / Ur2 * (xi1 - xi3 + xi4) + dc * xi2));
            }
        }

    }

    delete [] prim;    prim = nullptr;
    delete [] primL;    primL = nullptr;
    delete [] primR;    primR = nullptr;
    delete [] qConserveL;    qConserveL = nullptr;
    delete [] qConserveR;    qConserveR = nullptr;
    delete [] dq;    dq = nullptr;
}


//! GMRES
void GMRES_Roe_Scheme_ConservativeForm_wholeJacobian(FaceProxy *face_proxy, InviscidSchemeParameter *invSchemePara)
{
    const RDouble LITTLE    = 1.0E-10;
    int nTemperatureModel   = invSchemePara->GetNumberOfTemperatureModel();
    int nEquation           = invSchemePara->GetNumberOfTotalEquation();
    int nLaminar            = invSchemePara->GetNumberOfLaminar();
    int nLength             = invSchemePara->GetLength();
    int nm                  = invSchemePara->Getnm();
    int RoeEntropyFixMethod = invSchemePara->GetEntropyFixMethod();
    int nChemical           = invSchemePara->GetNumberOfChemical();
    int nPrecondition       = invSchemePara->GetIfPrecondition();
    int nIdealState         = invSchemePara->GetIfIdealState();

    RDouble refMachNumber     = invSchemePara->GetMachNumber();
    RDouble entrFixCoefLeft   = invSchemePara->GetLeftEntropyFixCoefficient();
    RDouble entrFixCoefRight  = invSchemePara->GetRightEntropyFixCoefficient();
    RDouble **qL              = invSchemePara->GetLeftQ();
    RDouble **qR              = invSchemePara->GetRightQ();
    RDouble *gamaL            = invSchemePara->GetLeftGama();
    RDouble *gamaR            = invSchemePara->GetRightGama();
    RDouble *pressureCoeffL   = face_proxy->GetPressureCoefficientL();
    RDouble *pressureCoeffR   = face_proxy->GetPressureCoefficientR();
    RDouble *xfn              = invSchemePara->GetFaceNormalX();
    RDouble *yfn              = invSchemePara->GetFaceNormalY();
    RDouble *zfn              = invSchemePara->GetFaceNormalZ();
    RDouble *vgn              = invSchemePara->GetFaceVelocity();
    RDouble **flux            = invSchemePara->GetFlux();
    RDouble **tl              = invSchemePara->GetLeftTemperature();
    RDouble **tr              = invSchemePara->GetRightTemperature();
    int *leftcellindexofFace  = face_proxy->GetLeftCellIndexOfFace();
    int *rightcellindexofFace = face_proxy->GetRightCellIndexOfFace();
    GeomProxy *geomProxy      = face_proxy->GetGeomProxy();
    RDouble *areas            = geomProxy->GetFaceArea();
    RDouble **dRdq            = face_proxy->GetJacobianMatrix();
    RDouble mach2             = refMachNumber * refMachNumber;

    using namespace IDX;
    using namespace GAS_SPACE;

    RDouble *prim       = new RDouble[nEquation]();
    RDouble *qConserveL = new RDouble[nEquation]();
    RDouble *qConserveR = new RDouble[nEquation]();
    RDouble *qCvL       = new RDouble[nEquation]();
    RDouble *qCvR       = new RDouble[nEquation]();
    RDouble *qCvLpert   = new RDouble[nEquation]();
    RDouble *qCvRpert   = new RDouble[nEquation]();
    RDouble *dFluxdCv   = new RDouble[nEquation]();
    RDouble dCvpert;

    RDouble *primL = new RDouble[nEquation]();
    RDouble *primR = new RDouble[nEquation]();
    RDouble *fluxp = new RDouble[nEquation]();    //! flux proxy
    RDouble temperaturesL[3] = { 0 }, temperaturesR[3] = { 0 };
    RDouble perturbScale = 0.001;

    for (int iFace = 0; iFace < nLength; ++iFace)
    {
        int le          = leftcellindexofFace[iFace];
        int re          = rightcellindexofFace[iFace];
        RDouble area    = areas[iFace];
        for (int iEqn = 0; iEqn < nEquation; ++iEqn)
        {
            primL[iEqn] = qL[iEqn][iFace];
            primR[iEqn] = qR[iEqn][iFace];
        }

        RDouble gmL = gamaL[iFace];
        RDouble gmR = gamaR[iFace];

        gas->Primitive2Conservative(primL, gmL, temperaturesL[ITV], temperaturesL[ITE], qCvL);
        gas->Primitive2Conservative(primR, gmR, temperaturesR[ITV], temperaturesR[ITE], qCvR);

        gas->Conservative2Primitive(qCvL, gmL, primL, temperaturesL);
        gas->Conservative2Primitive(qCvR, gmR, primR, temperaturesR);

        Cal_GMRES_Roe_Scheme_Flux(primL, primR, qCvL, qCvR,iFace, fluxp,
                                  face_proxy,invSchemePara);

        flux[IR][iFace] = fluxp[IR];
        flux[IU][iFace] = fluxp[IU];
        flux[IV][iFace] = fluxp[IV];
        flux[IW][iFace] = fluxp[IW];
        flux[IP][iFace] = fluxp[IP];

    }

    delete [] prim;           prim = nullptr;
    delete [] primL;          primL = nullptr;
    delete [] primR;          primR = nullptr;
    delete [] qConserveL;     qConserveL = nullptr;
    delete [] qConserveR;     qConserveR = nullptr;
    delete [] qCvL;           qCvL = nullptr;
    delete [] dFluxdCv;       dFluxdCv = nullptr;
    delete [] qCvLpert;       qCvLpert = nullptr;
    delete [] qCvRpert;       qCvRpert = nullptr;
    delete [] qCvR;           qCvR = nullptr;
    delete [] fluxp;          fluxp = nullptr;
}
//! GMRESAD CSR format  GMRESCSR
//! GMRES for finite difference considering nBoundFace using primitive variables with auto difference
void GMRES_Roe_Scheme_ConservativeForm(FaceProxy *face_proxy, InviscidSchemeParameter *invSchemePara)
{
    const RDouble LITTLE    = 1.0E-10;
    int nTemperatureModel   = invSchemePara->GetNumberOfTemperatureModel();
    int nEquation           = invSchemePara->GetNumberOfTotalEquation();
    int nLaminar            = invSchemePara->GetNumberOfLaminar();
    int nLength             = invSchemePara->GetLength();
    int nm                  = invSchemePara->Getnm();
    int RoeEntropyFixMethod = invSchemePara->GetEntropyFixMethod();
    int nChemical           = invSchemePara->GetNumberOfChemical();
    int nPrecondition       = invSchemePara->GetIfPrecondition();
    int nIdealState         = invSchemePara->GetIfIdealState();

    //! GMRESBoundary
    int localStart  = face_proxy->GetlocalStart();
    int localEnd    = face_proxy->GetlocalEnd();
    int nBoundFace  = face_proxy->GetnBoundFace();
    int nTotalCell  = face_proxy->GetnTotalCell();
    int nMid;
    //! GMRESBCorrection
    vector<int>* wallFaceIndex = face_proxy->GetWallFaceIndex();

    RDouble refMachNumber     = invSchemePara->GetMachNumber();
    RDouble entrFixCoefLeft   = invSchemePara->GetLeftEntropyFixCoefficient();
    RDouble entrFixCoefRight  = invSchemePara->GetRightEntropyFixCoefficient();
    RDouble **qL              = invSchemePara->GetLeftQ();
    RDouble **qR              = invSchemePara->GetRightQ();
    RDouble **qLC             = invSchemePara->GetLeftQC();    //! GMRESPassQC
    RDouble **qRC             = invSchemePara->GetRightQC();   //! GMRESPassQC
    RDouble *gamaL            = invSchemePara->GetLeftGama();
    RDouble *gamaR            = invSchemePara->GetRightGama();
    RDouble *pressureCoeffL   = face_proxy->GetPressureCoefficientL();
    RDouble *pressureCoeffR   = face_proxy->GetPressureCoefficientR();
    RDouble *xfn              = invSchemePara->GetFaceNormalX();
    RDouble *yfn              = invSchemePara->GetFaceNormalY();
    RDouble *zfn              = invSchemePara->GetFaceNormalZ();
    RDouble *vgn              = invSchemePara->GetFaceVelocity();
    RDouble **flux            = invSchemePara->GetFlux();
    RDouble **tl              = invSchemePara->GetLeftTemperature();
    RDouble **tr              = invSchemePara->GetRightTemperature();
    int *leftcellindexofFace  = face_proxy->GetLeftCellIndexOfFace();
    int *rightcellindexofFace = face_proxy->GetRightCellIndexOfFace();
    GeomProxy *geomProxy      = face_proxy->GetGeomProxy();
    RDouble *areas            = geomProxy->GetFaceArea();
    RDouble **dRdq            = face_proxy->GetJacobianMatrix();
    //! GMRESBoundary
    RDouble **dDdP             = face_proxy->GetdDdPMatrix();
    //! GMRES CSR
    vector<int> AI = face_proxy->GetJacobianAI4GMRES();
    vector<int> AJ = face_proxy->GetJacobianAJ4GMRES();

    // RDouble   dRdqL[nEquation][nEquation];
    RDouble dRdqR[nEquation][nEquation];
    RDouble dDdPlocal[nEquation][nEquation];
    RDouble mach2 = refMachNumber * refMachNumber;

    using namespace IDX;
    using namespace GAS_SPACE;

    RDouble  dFluxdpL[nEquation][nEquation];
    RDouble  dFluxdpR[nEquation][nEquation];
    RDouble  dFluxdpwall[nEquation][nEquation];

    //! auto difference allocation
    ADReal *primL  = new ADReal[nEquation]();
    ADReal *primR  = new ADReal[nEquation]();
    ADReal *fluxp  = new ADReal[nEquation]();    //! flux proxy

    RDouble *primLC = new RDouble[nEquation]();    //! GMRESPassQC
    RDouble *primRC = new RDouble[nEquation]();    //! GMRESPassQC

    RDouble temperaturesL[3] = { 0 }, temperaturesR[3] = { 0 };
    RDouble perturbScale = 0.001;

    RDouble **dqdcvL= new RDouble*[nEquation];
    RDouble **dqdcvR= new RDouble*[nEquation];
    for(int indexI=0; indexI < nEquation; indexI++)
    {
        dqdcvL[indexI] = new RDouble[nEquation];
        dqdcvR[indexI] = new RDouble[nEquation];
        
        for(int indexJ =0; indexJ < nEquation; indexJ ++)
        {
            dqdcvL[indexI][indexJ] = 0.0;
            dqdcvR[indexI][indexJ] = 0.0;
        }
    }

    //! GMRESBoundary
     if (localStart >= nBoundFace)
     {
         nMid = localStart;    //! a bug 01.07
     }
     else if (localEnd <= nBoundFace)
     {
         //! If they are all boundary faces.
         nMid = localEnd;
     }
     else
     {
         //! Part of them are boundary faces.
         nMid = nBoundFace;
     }

     nMid = nMid - localStart;

     for(int iFace = 0;iFace < nMid; ++iFace)
     {
        int le          = leftcellindexofFace[iFace];
        int re          = rightcellindexofFace[iFace];
        int colidx;    //! colume index for GMRES CSR

        RDouble area    = areas[iFace];
        for (int iEqn = 0; iEqn < nEquation; ++iEqn)
        {
            primL[iEqn] = qL[iEqn][iFace];
            primR[iEqn] = qR[iEqn][iFace];

            primLC[iEqn] = qLC[iEqn][iFace];    //! GMRESPassQC
            primRC[iEqn] = qRC[iEqn][iFace];    //! GMRESPassQC
        }

        //! Define sequence order of the independent variables
        for (int iEqn = 0; iEqn < nEquation; iEqn++)
        {
            primL[iEqn].diff(iEqn, 10);
            primR[iEqn].diff(iEqn + nEquation, 10);
        }
        
        RDouble gmL = gamaL[iFace];
        RDouble gmR = gamaR[iFace];

        Cal_GMRES_Roe_Scheme_Flux_PV(primL, primR, iFace, fluxp,
                                  face_proxy,invSchemePara);

        flux[IR][iFace] = fluxp[IR].val();
        flux[IU][iFace] = fluxp[IU].val();
        flux[IV][iFace] = fluxp[IV].val();
        flux[IW][iFace] = fluxp[IW].val();
        flux[IP][iFace] = fluxp[IP].val();

        for (int m = 0; m < nEquation; m++)
        {
            for (int n = 0; n < nEquation; n++)
            {
                dFluxdpL[m][n] = fluxp[m].dx(n) * area;
                dFluxdpR[m][n] = fluxp[m].dx(n + nEquation) * area;
            }
        }

        //! Contributions from the bc
        int indexre = (re-nTotalCell)*nEquation;
        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                dDdPlocal[indexI][indexJ] = dDdP[indexI][indexre+indexJ];
            }
        }

        //! GMRESBCorrection
        if(find(wallFaceIndex->begin(),wallFaceIndex->end(),iFace) != wallFaceIndex->end())
        {
            //! GMRESPV , obtain the dqdcv for the left cell
            gas->dPrimitive2dConservative(primLC,gmL,dqdcvL);

            //! Average of the dFluxdpL and dFluxdpR
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ =0; indexJ < nEquation; indexJ ++)
                {
                    dFluxdpwall[indexI][indexJ] = 0.5*(dFluxdpL[indexI][indexJ]
                                                     + dFluxdpR[indexI][indexJ]);
                    dFluxdpL[indexI][indexJ]    =  dFluxdpwall[indexI][indexJ];
                    dFluxdpR[indexI][indexJ]    =  dFluxdpwall[indexI][indexJ];
                }
            }

            //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
            //! then minus to the corresponding location in dRdq.
            int indexf = re * nEquation;    //! first index
            int indexs = le * nEquation;    //! second index
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, le);
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        /* dRdq[indexf+indexI][indexs+indexJ] -= dFluxdpwall[indexI][indexK]*
                                                                   dqdcvL[indexK][indexJ]; */
                        dRdq[indexI][colidx + indexJ] -= dFluxdpwall[indexI][indexK] * dqdcvL[indexK][indexJ];
                    }
                }
            }

            //! contributions from the bc, 
            //! dFluxdpL = dFluxdpL + dFluxdPR*dDdPlocal
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        dFluxdpL[indexI][indexJ]  += dFluxdpwall[indexI][indexK]*
                                                       dDdPlocal[indexK][indexJ];
                    }
                }
            }

            //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
            //! then add to the corresponding location in dRdq.
            indexf = le * nEquation;    //! first index
            indexs = le * nEquation;    //! second index
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, le);
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        /* dRdq[indexf+indexI][indexs+indexJ] += dFluxdpL[indexI][indexK]*
                                                                 dqdcvL[indexK][indexJ]; */
                        dRdq[indexI][colidx + indexJ] += dFluxdpL[indexI][indexK] * dqdcvL[indexK][indexJ];
                    }
                }
            }

             //! GMRESPV , obtain the dqdcv for the right cell
            gas->dPrimitive2dConservative(primRC,gmR,dqdcvR);

            //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
            //! then minus to the corresponding location in dRdq.
            indexf = re * nEquation;    //! first index
            indexs = re * nEquation;    //! second index
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, re);
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        /*dRdq[indexf+indexI][indexs+indexJ] -= dFluxdpR[indexI][indexK]*
                                                                 dqdcvR[indexK][indexJ];*/
                        dRdq[indexI][colidx + indexJ] -= dFluxdpR[indexI][indexK] * dqdcvR[indexK][indexJ];
                    }
                }
            }

            //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
            //! then add to the corresponding location in dRdq.
            indexf = le * nEquation;    //! first index
            indexs = re * nEquation;    //! second index
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, re);
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        /* dRdq[indexf+indexI][indexs+indexJ] += dFluxdpR[indexI][indexK]*
                                                                 dqdcvR[indexK][indexJ]; */
                        dRdq[indexI][colidx + indexJ] += dFluxdpR[indexI][indexK] * dqdcvR[indexK][indexJ];
                    }
                }
            }
            
        }
        else
        {
            //! GMRESPV , obtain the dqdcv for the left cell
            gas->dPrimitive2dConservative(primLC,gmL,dqdcvL);

            //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
            //! then minus to the corresponding location in dRdq.
            int indexf = re * nEquation;    //! first index
            int indexs = le * nEquation;    //! second index
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, le);
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        /* dRdq[indexf+indexI][indexs+indexJ] -= dFluxdpL[indexI][indexK]*
                                                                 dqdcvL[indexK][indexJ]; */
                        dRdq[indexI][colidx + indexJ] -= dFluxdpL[indexI][indexK] * dqdcvL[indexK][indexJ];
                    }
                }
            }

            //! contributions from the bc, 
            //! dFluxdpL = dFluxdpL + dFluxdPR*dDdPlocal
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        dFluxdpL[indexI][indexJ]  += dFluxdpR[indexI][indexK]*
                                                    dDdPlocal[indexK][indexJ];
                    }
                }
            }

            //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
            //! then add to the corresponding location in dRdq.
            indexf = le * nEquation;    //! first index
            indexs = le * nEquation;    //! second index
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, le);
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        /* dRdq[indexf+indexI][indexs+indexJ] += dFluxdpL[indexI][indexK]*
                                                                 dqdcvL[indexK][indexJ]; */
                        dRdq[indexI][colidx + indexJ] += dFluxdpL[indexI][indexK] * dqdcvL[indexK][indexJ];
                    }
                }
            }

            //! GMRESPV , obtain the dqdcv for the right cell
            gas->dPrimitive2dConservative(primRC,gmR,dqdcvR);

            //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
            //! then minus to the corresponding location in dRdq.
            indexf = re * nEquation;    //! first index
            indexs = re * nEquation;    //! second index
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, re);
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        /* dRdq[indexf+indexI][indexs+indexJ] -= dFluxdpR[indexI][indexK]*
                                                                 dqdcvR[indexK][indexJ]; */
                        dRdq[indexI][colidx + indexJ] -= dFluxdpR[indexI][indexK] * dqdcvR[indexK][indexJ];
                    }
                }
            }

            //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
            //! then add to the corresponding location in dRdq.
            indexf = le * nEquation;    //! first index
            indexs = re * nEquation;    //! second index
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, re);
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        /* dRdq[indexf+indexI][indexs+indexJ] += dFluxdpR[indexI][indexK]*
                                                                 dqdcvR[indexK][indexJ]; */
                        dRdq[indexI][colidx + indexJ] += dFluxdpR[indexI][indexK] * dqdcvR[indexK][indexJ];
                    }
                }
            }
        }
     }

     for(int iFace = nMid;iFace < nLength; ++iFace)
     {
        int le          = leftcellindexofFace[iFace];
        int re          = rightcellindexofFace[iFace];
        int colidx;
        RDouble area    = areas[iFace];
        for (int iEqn = 0; iEqn < nEquation; ++iEqn)
        {
            primL[iEqn] = qL[iEqn][iFace];
            primR[iEqn] = qR[iEqn][iFace];

            primLC[iEqn] = qLC[iEqn][iFace];    //! GMRESPassQC
            primRC[iEqn] = qRC[iEqn][iFace];    //! GMRESPassQC
        }

        //! Define sequence order of the independent variables
        for (int iEqn = 0; iEqn < nEquation; iEqn++)
        {
            primL[iEqn].diff(iEqn, 10);
            primR[iEqn].diff(iEqn + nEquation, 10);
        }

        RDouble gmL = gamaL[iFace];
        RDouble gmR = gamaR[iFace];

        Cal_GMRES_Roe_Scheme_Flux_PV(primL, primR, iFace, fluxp,
                                  face_proxy,invSchemePara);

        flux[IR][iFace] = fluxp[IR].val();
        flux[IU][iFace] = fluxp[IU].val();
        flux[IV][iFace] = fluxp[IV].val();
        flux[IW][iFace] = fluxp[IW].val();
        flux[IP][iFace] = fluxp[IP].val();
       
        for (int m = 0; m < nEquation; m++)
         {
             for (int n = 0; n < nEquation; n++)
             {
                 dFluxdpL[m][n] = fluxp[m].dx(n) * area;
                 dFluxdpR[m][n] = fluxp[m].dx(n + nEquation) * area;
             }
         }
  
        //! GMRESPV , obtain the dqdcv for the left cell
        gas->dPrimitive2dConservative(primLC,gmL,dqdcvL);

        //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
        //! then minus to the corresponding location in dRdq.
        int indexf = re * nEquation;    //! first index
        int indexs = le * nEquation;    //! second index
        colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, le);
        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            { 
                for(int indexK=0; indexK < nEquation; indexK++)
                {
                    dRdq[indexI][colidx + indexJ] -= dFluxdpL[indexI][indexK] * dqdcvL[indexK][indexJ];
                }
            }
        }

        //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
        //! then add to the corresponding location in dRdq.
        indexf = le * nEquation;    //! first index
        indexs = le * nEquation;    //! second index
        colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, le);
        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                for(int indexK=0; indexK < nEquation; indexK++)
                {
                    dRdq[indexI][colidx + indexJ] += dFluxdpL[indexI][indexK] * dqdcvL[indexK][indexJ];
                }
            }
        }

        //! GMRESPV , obtain the dqdcv for the right cell
        gas->dPrimitive2dConservative(primRC,gmR,dqdcvR);

        //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
        //! then minus to the corresponding location in dRdq.
        indexf = re * nEquation;    //! first index
        indexs = re * nEquation;    //! second index
        colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, re);
        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                for(int indexK=0; indexK < nEquation; indexK++)
                {
                    dRdq[indexI][colidx + indexJ] -= dFluxdpR[indexI][indexK] * dqdcvR[indexK][indexJ];
                }
            }
        }

        //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
        //! then add to the corresponding location in dRdq.
        indexf = le * nEquation;    //! first index
        indexs = re * nEquation;    //! second index
        colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, re);
        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                for(int indexK=0; indexK < nEquation; indexK++)
                {
                    dRdq[indexI][colidx + indexJ] += dFluxdpR[indexI][indexK] * dqdcvR[indexK][indexJ];
                }
            }
        }
     }

     delete[] primL;     primL = nullptr;
    delete [] primR;     primR = nullptr;
    delete [] primLC;    primLC = nullptr;
    delete [] primRC;    primRC = nullptr;
    delete [] fluxp;     fluxp = nullptr;

    for(int index = 0; index < nEquation; index++)
    {
        delete [] dqdcvL[index];   dqdcvL[index] = nullptr;
        delete [] dqdcvR[index];   dqdcvR[index] = nullptr;
    }
    delete [] dqdcvL;   dqdcvL = nullptr;
    delete [] dqdcvR;   dqdcvR = nullptr;
}


//! GMRESJac1st-nolim GMRES3D
void GMRES_INVScheme(FaceProxy *face_proxy, InviscidSchemeParameter *invSchemePara)
{
    const RDouble LITTLE    = 1.0E-10;
    int nTemperatureModel   = invSchemePara->GetNumberOfTemperatureModel();
    int nEquation           = invSchemePara->GetNumberOfTotalEquation();
    int nLaminar            = invSchemePara->GetNumberOfLaminar();
    int nLength             = invSchemePara->GetLength();
    int nm                  = invSchemePara->Getnm();
    int RoeEntropyFixMethod = invSchemePara->GetEntropyFixMethod();
    int nChemical           = invSchemePara->GetNumberOfChemical();
    int nPrecondition       = invSchemePara->GetIfPrecondition();
    int nIdealState         = invSchemePara->GetIfIdealState();
    int limiterType         = invSchemePara->GetLimiterType();
    int unstructScheme      = GlobalDataBase::GetIntParaFromDB("uns_scheme");    //! GMRES3D
    UnstructBCSet *unstructBCSet       = invSchemePara->GetUnstructBCSet();    //! GMRES3D judge the type of the bc
    int *bcRegionIDofBCFace            = unstructBCSet->GetBCFaceInfo();    //! GMRES3D
    RDouble wallTemperature            = GlobalDataBase::GetDoubleParaFromDB("wallTemperature");
    RDouble refDimensionalTemperature  = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");
    RDouble coefficientOfStateEquation = gas->GetCoefficientOfStateEquation();

    //! GMRES1st2nd
    if( limiterType == ILMT_NOLIM || limiterType == ILMT_VENCAT )
    {
        int newIterStep = GlobalDataBase::GetIntParaFromDB("newnstep");
        int flowInitStep = GlobalDataBase::GetIntParaFromDB("flowInitStep");
        
        if(newIterStep <= flowInitStep)
        {
            limiterType = ILMT_FIRST;
        }
    }

    //! GMRESBoundary
    int localStart = face_proxy->GetlocalStart();
    int localEnd   = face_proxy->GetlocalEnd();
    int nBoundFace = face_proxy->GetnBoundFace();
    int nTotalCell = face_proxy->GetnTotalCell();
    int nMid;
    //! GMRESBCorrection
    vector<int>* wallFaceIndex = face_proxy->GetWallFaceIndex();

    RDouble refMachNumber     = invSchemePara->GetMachNumber();
    RDouble entrFixCoefLeft   = invSchemePara->GetLeftEntropyFixCoefficient();
    RDouble entrFixCoefRight  = invSchemePara->GetRightEntropyFixCoefficient();
    RDouble **qL              = invSchemePara->GetLeftQ();
    RDouble **qR              = invSchemePara->GetRightQ();
    RDouble **qLC             = invSchemePara->GetLeftQC();    //! GMRESPassQC
    RDouble **qRC             = invSchemePara->GetRightQC();   //! GMRESPassQC
    RDouble *gamaL            = invSchemePara->GetLeftGama();
    RDouble *gamaR            = invSchemePara->GetRightGama();
    RDouble *pressureCoeffL   = face_proxy->GetPressureCoefficientL();
    RDouble *pressureCoeffR   = face_proxy->GetPressureCoefficientR();
    RDouble *xfn              = invSchemePara->GetFaceNormalX();
    RDouble *yfn              = invSchemePara->GetFaceNormalY();
    RDouble *zfn              = invSchemePara->GetFaceNormalZ();
    RDouble *vgn              = invSchemePara->GetFaceVelocity();
    RDouble **flux            = invSchemePara->GetFlux();
    RDouble **tl              = invSchemePara->GetLeftTemperature();
    RDouble **tr              = invSchemePara->GetRightTemperature();
    int *leftcellindexofFace  = face_proxy->GetLeftCellIndexOfFace();
    int *rightcellindexofFace = face_proxy->GetRightCellIndexOfFace();
    GeomProxy* geomProxy      = face_proxy->GetGeomProxy();
    RDouble* areas            = geomProxy->GetFaceArea();

    //! GMRESnolim  xfc xcc should use the absolute face index
    RDouble *xfc               = invSchemePara->GetFaceCenterX();
    RDouble *yfc               = invSchemePara->GetFaceCenterY();
    RDouble *zfc               = invSchemePara->GetFaceCenterZ();
    RDouble *xcc               = invSchemePara->GetCellCenterX();
    RDouble *ycc               = invSchemePara->GetCellCenterY();
    RDouble *zcc               = invSchemePara->GetCellCenterZ();
    RDouble *nsx               = invSchemePara->GetFaceNormalXAbs();
    RDouble *nsy               = invSchemePara->GetFaceNormalYAbs();
    RDouble *nsz               = invSchemePara->GetFaceNormalZAbs();
    RDouble *ns                = invSchemePara->GetFaceAreaAbs();
    RDouble *vol               = invSchemePara->GetVolume();
    RDouble *gama              = invSchemePara->GetGama();
    RDouble **prim             = invSchemePara->GetPrimitive();
    vector<int> *neighborCells = invSchemePara->GetNeighborCells();
    vector<int> *neighborFaces = invSchemePara->GetNeighborFaces();
    vector<int> *neighborLR    = invSchemePara->GetNeighborLR();
    int *qlsign                = invSchemePara->GetLeftQSign();
    int *qrsign                = invSchemePara->GetRightQSign();

    //! GMRES2ndCorrection
    vector<int> BCLeftCells  = invSchemePara->GetBCLeftCells();
    vector<int> BCRightCells = invSchemePara->GetBCRightCells();
    vector<int> BCFaces      = invSchemePara->GetBCFaces();

    //! GMRESVenkat
    RDouble ** dqdx         = invSchemePara->GetGradientX();
    RDouble ** dqdy         = invSchemePara->GetGradientY();
    RDouble ** dqdz         = invSchemePara->GetGradientZ();
    int limitModel          = invSchemePara->GetLimitModel();
    int limitVector         = invSchemePara->GetLimitVector();
    bool usingVectorLimiter = (limitVector == 1);
    int limitnVariable      = 1;
    ADReal *limitL;
    ADReal *limitR;
    int     ivencat;
    RDouble venkatCoeff;
    RDouble averageVolume = invSchemePara->GetAverageVolume();
    RDouble nExp          = 3.0; 
    int dimension         = PHSPACE::GetDim();

    const int SCALAR_LIMITER    = 0;
    const int VECTOR_LIMITER    = 1;

    const int LIMIT_BY_RHO_P    = 0;
    const int LIMIT_BY_EACH_VAR = 1;

    int reconstructMethod       = GlobalDataBase::GetIntParaFromDB("reconmeth");
    const int RECONSTRUCT_USING_RESPECTIVE_LIMITER = 0;

    if ( limiterType == ILMT_VENCAT )
    {
        if( limitVector == 0  )
        {
            TK_Exit::ExceptionExit("SCALAR_LIMITER for Venkat is not yet implemented\n");
        }
        
        if (limitVector != 0 || limitModel != 0)
        {
            limitnVariable = nEquation;
        }

        limitL = new ADReal[limitnVariable];
        limitR = new ADReal[limitnVariable];

        ivencat     = GlobalDataBase::GetIntParaFromDB("ivencat");
        venkatCoeff = GlobalDataBase::GetDoubleParaFromDB("venkatCoeff");
    }

    RDouble dxL, dyL, dzL, dxR, dyR, dzR;
    int Lsign;
    int Rsign;

    RDouble **dRdq = face_proxy->GetJacobianMatrix();
    //! GMRESBoundary
    RDouble **dDdP = face_proxy->GetdDdPMatrix();
    //! GMRES CSR
    vector<int> AI = face_proxy->GetJacobianAI4GMRES();
    vector<int> AJ = face_proxy->GetJacobianAJ4GMRES();

    //! GMRESJac1st
    RDouble** dRdq1st = face_proxy->GetJacobianMatrix1st();
    vector<int> AI1st = face_proxy->GetJacobianAI1st4GMRES();
    vector<int> AJ1st = face_proxy->GetJacobianAJ1st4GMRES();
    int JacOrder = face_proxy->GetJacobianOrder();

    RDouble dRdqR[nEquation][nEquation];
    RDouble dDdPlocal[nEquation][nEquation];
    RDouble dGdqx, dGdqy, dGdqz;
    RDouble mach2 = refMachNumber * refMachNumber;

    using namespace IDX;
    using namespace GAS_SPACE;

    RDouble  dFluxdpL[nEquation][nEquation];
    RDouble  dFluxdpR[nEquation][nEquation];
    RDouble  dFluxdtmp[nEquation][nEquation];
    RDouble  dFluxdpwall[nEquation][nEquation];

    //! GMRESVenkat
    RDouble dFluxdgradxL[nEquation][nEquation];
    RDouble dFluxdgradyL[nEquation][nEquation];
    RDouble dFluxdgradzL[nEquation][nEquation];
    RDouble dFluxdgradxR[nEquation][nEquation];
    RDouble dFluxdgradyR[nEquation][nEquation];
    RDouble dFluxdgradzR[nEquation][nEquation];
    RDouble dFluxdminL[nEquation][nEquation];
    RDouble dFluxdmaxL[nEquation][nEquation];
    RDouble dFluxdminR[nEquation][nEquation];
    RDouble dFluxdmaxR[nEquation][nEquation];
    RDouble dGradtmp[nEquation][nEquation];
    RDouble dGradtmpBC[nEquation][nEquation];
    RDouble dFluxdMinMax[nEquation][nEquation];

    //! auto difference allocation
    ADReal *primL   = new ADReal[nEquation]();
    ADReal *primR   = new ADReal[nEquation]();
    ADReal *fluxp   = new ADReal[nEquation]();    //! flux proxy
    ADReal *primTry = new ADReal[nEquation]();

    //! GMRESVenkat
    ADReal *dqdxL = new ADReal[nEquation];
    ADReal *dqdyL = new ADReal[nEquation];
    ADReal *dqdzL = new ADReal[nEquation];

    ADReal *dqdxR = new ADReal[nEquation];
    ADReal *dqdyR = new ADReal[nEquation];
    ADReal *dqdzR = new ADReal[nEquation];

    ADReal dminL;
    ADReal dminR;
    ADReal dmaxL;
    ADReal dmaxR;

    // if ( limitVector != 0 || limitModel != 0 )
    ADReal *minL = new ADReal[nEquation]();
    ADReal *maxL = new ADReal[nEquation]();
    ADReal *minR = new ADReal[nEquation]();
    ADReal *maxR = new ADReal[nEquation]();

    int *IndexMinL = new int [nEquation];
    int *IndexMaxL = new int [nEquation];
    int *IndexMinR = new int [nEquation];
    int *IndexMaxR = new int [nEquation];
    int *LimitSignL = new int [nEquation];
    int *LimitSignR = new int [nEquation];

    int nIndependentVars;
    int halfIndependentVars ;
    if ( limitVector != 0 || limitModel != 0 )
    {
        nIndependentVars    = nEquation * 2 + nEquation  * 6 + nEquation * 4;
        halfIndependentVars = nIndependentVars / 2;
    }

    RDouble *primLC = new RDouble[nEquation]();    //! GMRESPassQC
    RDouble *primRC = new RDouble[nEquation]();    //! GMRESPassQC
    RDouble *primN  = new RDouble[nEquation]();

    RDouble temperaturesL[3] = { 0 }, temperaturesR[3] = { 0 };
    RDouble perturbScale = 0.001;

    RDouble **dqdcvL= new RDouble*[nEquation];
    RDouble **dqdcvR= new RDouble*[nEquation];
    RDouble **dqdcvN= new RDouble*[nEquation];    //! dqdcv for neighbors
    RDouble **dgraddqx = new RDouble*[nEquation];
    RDouble **dgraddqy = new RDouble*[nEquation];
    RDouble **dgraddqz = new RDouble*[nEquation];
    for(int indexI=0; indexI < nEquation; indexI++)
    {
        dqdcvL[indexI]   = new RDouble[nEquation];
        dqdcvR[indexI]   = new RDouble[nEquation];
        dqdcvN[indexI]   = new RDouble[nEquation];
        dgraddqx[indexI] = new RDouble[nEquation];
        dgraddqy[indexI] = new RDouble[nEquation];
        dgraddqz[indexI] = new RDouble[nEquation];

        for(int indexJ =0; indexJ < nEquation; indexJ ++)
        {
            dqdcvL[indexI][indexJ] = 0.0;
            dqdcvR[indexI][indexJ] = 0.0;
            dqdcvN[indexI][indexJ] = 0.0;
            dgraddqx[indexI][indexJ] = 0.0;
            dgraddqy[indexI][indexJ] = 0.0;
            dgraddqz[indexI][indexJ] = 0.0;
        }
    }

    //! GMRESBoundary
    if (localStart >= nBoundFace)
    {
        nMid = localStart;    //! a bug 01.07
    }
    else if (localEnd <= nBoundFace)
    {
        //! If they are all boundary faces.
        nMid = localEnd;
    }
    else
    {
        //! Part of them are boundary faces.
        nMid = nBoundFace;
    }

    nMid = nMid - localStart;

    for (int iFace = 0; iFace < nMid; ++iFace)
    {
        int le = leftcellindexofFace[iFace];
        int re = rightcellindexofFace[iFace];
        int colidx;    //! colume index for GMRES CSR
        int colidx1st; //! GMRESJac1st
        int kFace = iFace + localStart;    //! absolute face index

        UnstructBC* bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[kFace]);
        int bcType = bcRegion->GetBCType();

        RDouble area = areas[iFace];

        if ( limiterType == ILMT_VENCAT && 
             ( bcType == PHENGLEI::OUTFLOW || bcType == PHENGLEI::INTERFACE ) )
        {
            if (limitVector != 0 || limitModel != 0)
            {
                for (int nVar = 0; nVar < limitnVariable; ++ nVar)
                {
                    IndexMinL[nVar] = le;
                    IndexMaxL[nVar] = le;
                    RDouble pminL = prim[nVar][le];
                    RDouble pmaxL = prim[nVar][le];
                    for ( int index = 0; index < neighborCells[le].size(); index++ )
                    {
                        int neighborIndex       = neighborCells[le][index];
                        int neighborFaceIndex   = neighborFaces[le][index];
                        if( neighborFaceIndex >= nBoundFace )
                        {
                            RDouble  primNeighbor    = prim[nVar][neighborIndex];

                            if( pminL > primNeighbor)
                            {
                                pminL           = primNeighbor;
                                IndexMinL[nVar] = neighborIndex;
                            }

                            if( pmaxL < primNeighbor)
                            {
                                pmaxL = primNeighbor;
                                IndexMaxL[nVar] = neighborIndex;
                            }
                        }
                    }

                    IndexMinR[nVar] = re;
                    IndexMaxR[nVar] = re;
                    RDouble pminR = prim[nVar][re];
                    RDouble pmaxR = prim[nVar][re];
                    for ( int index = 0; index < neighborCells[re].size(); index++ )
                    {
                        int neighborIndex       = neighborCells[re][index];
                        int neighborFaceIndex   = neighborFaces[re][index];
                        if( neighborFaceIndex >= nBoundFace )
                        {
                            RDouble  primNeighbor    = prim[nVar][neighborIndex];

                            if( pminR > primNeighbor)
                            {
                                pminR           = primNeighbor;
                                IndexMinR[nVar] = neighborIndex;
                            }

                            if( pmaxR < primNeighbor)
                            {
                                pmaxR           = primNeighbor;
                                IndexMaxR[nVar] = neighborIndex;
                            }

                        }
                    }
                }
            }
        }
 
        //! GMRESVenkat pay attention to the interface boundary condition
        if ( limiterType == ILMT_VENCAT && 
             ( bcType == PHENGLEI::OUTFLOW || bcType == PHENGLEI::INTERFACE ) ) 
        {
            for (int iEqn = 0; iEqn < nEquation; ++iEqn)
            {
                primL[iEqn] = qL[iEqn][iFace];
                primR[iEqn] = qR[iEqn][iFace];

                primLC[iEqn] = qLC[iEqn][iFace];    //! GMRESPassQC
                primRC[iEqn] = qRC[iEqn][iFace];    //! GMRESPassQC

                dqdxL[iEqn] = dqdx[iEqn][le];
                dqdyL[iEqn] = dqdy[iEqn][le];
                dqdzL[iEqn] = dqdz[iEqn][le];
                dqdxR[iEqn] = dqdx[iEqn][re];
                dqdyR[iEqn] = dqdy[iEqn][re];
                dqdzR[iEqn] = dqdz[iEqn][re];

                minL[iEqn] = prim[iEqn][IndexMinL[iEqn]];
                maxL[iEqn] = prim[iEqn][IndexMaxL[iEqn]];
                minR[iEqn] = prim[iEqn][IndexMinR[iEqn]];
                maxR[iEqn] = prim[iEqn][IndexMaxR[iEqn]];
            }
        }
        else
        {
            for (int iEqn = 0; iEqn < nEquation; ++iEqn)
            {
                primL[iEqn] = qL[iEqn][iFace];
                primR[iEqn] = qR[iEqn][iFace];

                primLC[iEqn] = qLC[iEqn][iFace];    //! GMRESPassQC
                primRC[iEqn] = qRC[iEqn][iFace];    //! GMRESPassQC
            }
        }

        //! define sequence order of the independent variables
        //! GMRESVenkat
        if ( limiterType == ILMT_VENCAT && 
             (  bcType == PHENGLEI::OUTFLOW || bcType == PHENGLEI::INTERFACE ) )
        {
            //! define the sequence of the independent variables
            for (int m = 0; m < nEquation; m++)
            {
                primL[m].diff(m, nIndependentVars);
                primR[m].diff(halfIndependentVars + m, nIndependentVars);
            }

            for (int m = 0; m < nEquation; m++)
            {
                dqdxL[m].diff(nEquation + m, nIndependentVars);
                dqdyL[m].diff(nEquation + nEquation + m, nIndependentVars);
                dqdzL[m].diff(nEquation + 2 * nEquation + m, nIndependentVars);
                minL[m].diff(nEquation * 4 + m, nIndependentVars);
                maxL[m].diff(nEquation * 5 + m, nIndependentVars);
                dqdxR[m].diff(halfIndependentVars + nEquation + m, nIndependentVars);
                dqdyR[m].diff(halfIndependentVars + nEquation + nEquation + m, nIndependentVars);
                dqdzR[m].diff(halfIndependentVars + nEquation + 2 * nEquation + m, nIndependentVars);
                minR[m].diff(halfIndependentVars + nEquation * 4 + m, nIndependentVars);
                maxR[m].diff(halfIndependentVars + nEquation * 5 + m, nIndependentVars);
            }
        }
        else
        {
            for (int iEqn = 0; iEqn < nEquation; iEqn++)
            {
                primL[iEqn].diff(iEqn, 10);
                primR[iEqn].diff(iEqn + nEquation, 10);
            }
        }

        //! Boundaryqlqrfix
        if ( bcType == PHENGLEI::SYMMETRY )
        {
            ADReal vn = 2.0*( nsx[kFace]*primL[IU] + nsy[kFace]*primL[IV] + nsz[kFace]*primL[IW]);

            primR[IR] = primL[IR];
            primR[IU] = primL[IU] - vn * nsx[kFace];
            primR[IV] = primL[IV] - vn * nsy[kFace];
            primR[IW] = primL[IW] - vn * nsz[kFace];
            primR[IP] = primL[IP];
        }

        int viscousType = GlobalDataBase::GetIntParaFromDB("viscousType");

        if( bcType == PHENGLEI::SOLID_SURFACE && viscousType != INVISCID )
        {
            ADReal uWall = 0.0;
            ADReal vWall = 0.0;
            ADReal wWall = 0.0;
            Data_Param *bcData  = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[kFace])->GetBCParamDataBase();
            if(bcData)
            {
                if (bcData->IsExist("uWall", PHDOUBLE, 1))
                {
                    bcData->GetData("uWall", &uWall, PHDOUBLE, 1);
                    bcData->GetData("vWall", &vWall, PHDOUBLE, 1);
                    bcData->GetData("wWall", &wWall, PHDOUBLE, 1);
                }
            }

            ADReal velocityXWall = uWall;
            ADReal velocityYWall = vWall;
            ADReal velocityZWall = wWall;

            primL[IU] = velocityXWall;
            primL[IV] = velocityYWall;
            primL[IW] = velocityZWall;

            primR[IU] = velocityXWall;
            primR[IV] = velocityYWall;
            primR[IW] = velocityZWall;

            if (wallTemperature <= 0.0)
            {
                //! Viscous WALL, adiabatic.
            }
            else
            {
                //! Iso-thermal wall.
                ADReal tw = wallTemperature / refDimensionalTemperature;

                // for (int m = 0; m < nEquation; ++ m)
                // {
                //     primitiveVariable[m] = qL[m][jFace];
                // }

                ADReal omav = one;
                using namespace GAS_SPACE;
                // if (nChemical == 1)
                // {
                //     omav = gas->ComputeMolecularWeightReciprocal(primitiveVariable);
                // }

                ADReal pressureWall = primL[IP];
                ADReal rhoWall      = pressureWall / (coefficientOfStateEquation * tw * omav);

                primL[IR]  = rhoWall;
                primR[IR]  = rhoWall;
            }
        }

        //! GMRESVenkat 
        if ( limiterType == ILMT_VENCAT && 
             ( bcType == PHENGLEI::OUTFLOW || bcType == PHENGLEI::INTERFACE ) )
        {
            //! Initialize the limiter with one
            for ( int m = 0; m < limitnVariable; m++ )
            {
                limitL[m] = 1.0;
                limitR[m] = 1.0;
            }

            //! Limiter calculation
            if ( limitVector == SCALAR_LIMITER )
            {
                TK_Exit::ExceptionExit("SCALAR_LIMITER for Venkat is not yet implemented\n");
            }
            else if ( limitVector == VECTOR_LIMITER )
            {
                for (int nVar = 0; nVar < limitnVariable; ++ nVar)
                {
                    if( IndexMinL[nVar] == le)
                    {
                        dminL = 0;
                    }
                    else
                    {
                        dminL = minL[nVar] - primL[nVar];
                    }

                    if( IndexMaxL[nVar] == le)
                    {
                        dmaxL = 0;
                    }
                    else
                    {
                        dmaxL = maxL[nVar] - primL[nVar];
                    }

                    // !!!!! pay attention to the interface
                    // if( bcType == PHENGLEI::INTERFACE || bcType < 0 )
                    // {
                    //     dminL = MIN(dminL, primR[nVar]);
                    //     dmaxL = MIN(dmaxL, primR[nVar]);
                    // }

                    // for (int index = 0; index < neighborCells[le].size(); index++)
                    // {
                    //     int neighborIndex       = neighborCells[le][index];
                    //     if( neighborIndex != re)
                    //     {
                    //         ADReal  primNeighbor    = prim[nVar][neighborIndex];
                    //         dminL = MIN(dminL, primNeighbor);
                    //         dmaxL = MAX(dmaxL, primNeighbor);
                    //     }
                    // }

                    ADReal eps_tmp = venkatCoeff * venkatCoeff * venkatCoeff;
                    eps_tmp *= 1.0/pow(pow(averageVolume, 1.0/dimension), nExp);

                    ADReal epsCellL =  eps_tmp*primL[nVar]*primL[nVar];
                    epsCellL        *= pow(pow(vol[le],1.0/dimension),nExp);

                    LimitSignL[nVar] = 0;
                    for (int index = 0; index < neighborCells[le].size(); index++)
                    {
                        int FaceIndex = neighborFaces[le][index];
                        RDouble dx = xfc[FaceIndex] - xcc[le];
                        RDouble dy = yfc[FaceIndex] - ycc[le];
                        RDouble dz = zfc[FaceIndex] - zcc[le];

                        ADReal dqFace = dqdxL[nVar] * dx + dqdyL[nVar] * dy + dqdzL[nVar] * dz;

                        if( dqFace > SMALL )
                        {
                                ADReal x2   = dmaxL*dmaxL;
                                ADReal xy   = dmaxL*dqFace;
                                ADReal y2   = dqFace*dqFace; 
                                ADReal temp = (x2 + epsCellL + 2.0 * xy) / (x2 + 2.0 * y2 + xy + epsCellL);
                                
                                if( temp < limitL[nVar] )
                                {
                                     limitL[nVar]     = temp;
                                     LimitSignL[nVar] = 1;
                                }
                        }
                        else if( dqFace < -SMALL )
                        {
                                ADReal x2   = dminL*dminL;
                                ADReal xy   = dminL*dqFace;
                                ADReal y2   = dqFace*dqFace; 
                                ADReal temp = (x2 + epsCellL + 2.0 * xy) / (x2 + 2.0 * y2 + xy + epsCellL);
                                
                                if( temp < limitL[nVar] )
                                {
                                     limitL[nVar]       = temp;
                                     LimitSignL[nVar]   = -1;
                                }
                        }
                    }
                }
            }
            //! Reconstruct the left & right state using the limiter
            if (reconstructMethod == RECONSTRUCT_USING_RESPECTIVE_LIMITER)
            {
                //! reconstruct left state
                RDouble dxL = xfc[kFace] - xcc[le];
                RDouble dyL = yfc[kFace] - ycc[le];
                RDouble dzL = zfc[kFace] - zcc[le];

                for (int m = 0; m < nEquation; ++ m)
                {
                    primTry[m] = primL[m];
                }

                for (int m = 0; m < nEquation; ++ m)
                {
                    if(usingVectorLimiter)
                    {
                        primTry[m] += limitL[m] * (dqdxL[m] * dxL + dqdyL[m] * dyL + dqdzL[m] * dzL);
                    }
                    else
                    {
                        primTry[m] += limitL[0] * (dqdxL[m] * dxL + dqdyL[m] * dyL + dqdzL[m] * dzL);
                    }
                }

                if ( (primTry[IR] > 0.0)&&(primTry[IP] > 0.0) )
                {
                    for (int m = 0; m < nEquation; ++ m)
                    {
                        primL[m] = primTry[m];
                    }
                }

                //! reconstruct right state
                RDouble dxR = xfc[kFace] - xcc[re];
                RDouble dyR = yfc[kFace] - ycc[re];
                RDouble dzR = zfc[kFace] - zcc[re];

                for (int m = 0; m < nEquation; ++ m)
                {
                    primTry[m] = primR[m];
                }

                for (int m = 0; m < nEquation; ++ m)
                {
                    if(usingVectorLimiter)
                    {
                        primTry[m] += limitR[m] * (dqdxR[m] * dxR + dqdyR[m] * dyR + dqdzR[m] * dzR);
                    }
                    else
                    {
                        primTry[m] += limitR[0] * (dqdxR[m] * dxR + dqdyR[m] * dyR + dqdzR[m] * dzR);
                    }
                }

                if ( (primTry[IR] > 0.0)&&(primTry[IP] > 0.0) )
                {
                    for (int m = 0; m < nEquation; ++ m)
                    {
                        primR[m] = primTry[m];
                    }
                }
            }
            else
            {
                //! reconstruct left state
                RDouble dxL = xfc[kFace] - xcc[le];
                RDouble dyL = yfc[kFace] - ycc[le];
                RDouble dzL = zfc[kFace] - zcc[le];

                for (int m = 0; m < nEquation; ++ m)
                {
                    primTry[m] = primL[m];
                }

                for (int m = 0; m < nEquation; ++ m)
                {
                    if(usingVectorLimiter)
                    {
                        primTry[m] += MIN(limitL[m], limitR[m]) * (dqdxL[m] * dxL + dqdyL[m] * dyL + dqdzL[m] * dzL);
                    }
                    else
                    {
                        primTry[m] += MIN(limitL[0], limitR[0]) * (dqdxL[m] * dxL + dqdyL[m] * dyL + dqdzL[m] * dzL);
                    }
                }

                if ( (primTry[IR] > 0.0)&&(primTry[IP] > 0.0) )
                {
                    for (int m = 0; m < nEquation; ++ m)
                    {
                        primL[m] = primTry[m];
                    }
                }

                //! reconstruct right state
                RDouble dxR = xfc[kFace] - xcc[re];
                RDouble dyR = yfc[kFace] - ycc[re];
                RDouble dzR = zfc[kFace] - zcc[re];

                for (int m = 0; m < nEquation; ++ m)
                {
                    primTry[m] = primR[m];
                }

                for (int m = 0; m < nEquation; ++ m)
                {
                    if(usingVectorLimiter)
                    {
                        primTry[m] += MIN(limitR[m], limitL[m]) * (dqdxR[m] * dxR + dqdyR[m] * dyR + dqdzR[m] * dzR);
                    }
                    else
                    {
                        primTry[m] += MIN(limitR[0], limitL[0]) * (dqdxR[m] * dxR + dqdyR[m] * dyR + dqdzR[m] * dzR);
                    }
                }

                if ( (primTry[IR] > 0.0)&&(primTry[IP] > 0.0) )
                {
                    for (int m = 0; m < nEquation; ++ m)
                    {
                        primR[m] = primTry[m];
                    }
                }
            }
        }

        RDouble gmL = gamaL[iFace];
        RDouble gmR = gamaR[iFace];

        //! GMRES3D
        if( unstructScheme == ISCHEME_GMRES_ROE )
        {
            Cal_GMRES_Roe_Scheme_Flux_PV(primL, primR, iFace, fluxp,
                                        face_proxy, invSchemePara);
        }
        else if( unstructScheme == ISCHEME_GMRES_Steger )
        {
            Cal_GMRES_Steger_Scheme_Flux_PV(primL, primR, iFace, fluxp,
                                        face_proxy, invSchemePara);
        }
        else
        {
            TK_Exit::ExceptionExit("Error: this inv-scheme is not exist !\n", true);
        }
        

        flux[IR][iFace] = fluxp[IR].val();
        flux[IU][iFace] = fluxp[IU].val();
        flux[IV][iFace] = fluxp[IV].val();
        flux[IW][iFace] = fluxp[IW].val();
        flux[IP][iFace] = fluxp[IP].val();

        if ( limiterType == ILMT_VENCAT && 
             ( bcType == PHENGLEI::OUTFLOW || bcType == PHENGLEI::INTERFACE ) )
        {
            for (int m = 0; m < nEquation; m++)
            {
                for (int n = 0; n < nEquation; ++ n)
                {
                    dFluxdpL[m][n] = fluxp[m].dx(n) * area;
                    dFluxdpR[m][n] = fluxp[m].dx(halfIndependentVars + n) * area;
                }

                for (int n = 0; n < nEquation; n++)
                {
                    dFluxdgradxL[m][n] = fluxp[m].dx(nEquation + n) * area;
                    dFluxdgradyL[m][n] = fluxp[m].dx(nEquation + nEquation + n) * area;
                    dFluxdgradzL[m][n] = fluxp[m].dx(nEquation + 2 * nEquation  + n) * area;
                    dFluxdminL[m][n] = fluxp[m].dx(4 * nEquation  + n) * area;
                    dFluxdmaxL[m][n] = fluxp[m].dx(5 * nEquation  + n) * area;
                    dFluxdgradxR[m][n] = fluxp[m].dx(halfIndependentVars + nEquation + n) * area;
                    dFluxdgradyR[m][n] = fluxp[m].dx(halfIndependentVars + nEquation + nEquation + n) * area;
                    dFluxdgradzR[m][n] = fluxp[m].dx(halfIndependentVars + nEquation + 2 * nEquation + n) * area;
                    dFluxdminR[m][n] = fluxp[m].dx(halfIndependentVars + 4 * nEquation  + n) * area;
                    dFluxdmaxR[m][n] = fluxp[m].dx(halfIndependentVars + 5 * nEquation  + n) * area;
                }
            }
        }
        else
        {
            for (int m = 0; m < nEquation; m++)
            {
                for (int n = 0; n < nEquation; n++)
                {
                    dFluxdpL[m][n] = fluxp[m].dx(n) * area;
                    dFluxdpR[m][n] = fluxp[m].dx(n + nEquation) * area;
                }
            }
        }

        //! contributions from the bc
        int indexre = (re - nTotalCell) * nEquation;
        for (int indexI = 0; indexI < nEquation; indexI++)
        {
            for (int indexJ = 0; indexJ < nEquation; indexJ++)
            {
                dDdPlocal[indexI][indexJ] = dDdP[indexI][indexre + indexJ];
            }
        }

        //! GMRESBCorrection
        //if (find(wallFaceIndex->begin(), wallFaceIndex->end(), iFace) != wallFaceIndex->end())
        if ( bcType == PHENGLEI::SOLID_SURFACE || bcType == PHENGLEI::SYMMETRY )
        {
            //! GMRESPV , obtain the dqdcv for the left cell
            gas->dPrimitive2dConservative(primLC, gmL, dqdcvL);

            //! average of the dFluxdpL and dFluxdpR
            for (int indexI = 0; indexI < nEquation; indexI++)
            {
                for (int indexJ = 0; indexJ < nEquation; indexJ++)
                {
                    dFluxdpwall[indexI][indexJ] = 0.5 * (dFluxdpL[indexI][indexJ] + dFluxdpR[indexI][indexJ]);
                    dFluxdpL[indexI][indexJ] = dFluxdpwall[indexI][indexJ];
                    dFluxdpR[indexI][indexJ] = dFluxdpwall[indexI][indexJ];
                }
            }

            //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
            //! then minus to the corresponding location in dRdq.
            int indexf = re * nEquation;    //! first index
            int indexs = le * nEquation;    //! second index
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, le);
            for (int indexI = 0; indexI < nEquation; indexI++)
            {
                for (int indexJ = 0; indexJ < nEquation; indexJ++)
                {
                    for (int indexK = 0; indexK < nEquation; indexK++)
                    {
                        /* dRdq[indexf+indexI][indexs+indexJ] -= dFluxdpwall[indexI][indexK]*
                                                                dqdcvL[indexK][indexJ]; */
                        dRdq[indexI][colidx + indexJ] -= dFluxdpwall[indexI][indexK] * dqdcvL[indexK][indexJ];
                    }
                }
            }

            //! contributions from the bc,
            //! dFluxdpL = dFluxdpL + dFluxdPR*dDdPlocal
            for (int indexI = 0; indexI < nEquation; indexI++)
            {
                for (int indexJ = 0; indexJ < nEquation; indexJ++)
                {
                    for (int indexK = 0; indexK < nEquation; indexK++)
                    {
                        dFluxdpL[indexI][indexJ] += dFluxdpwall[indexI][indexK] *
                                                    dDdPlocal[indexK][indexJ];
                    }
                }
            }

            //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
            //! then add to the corresponding location in dRdq.
            indexf = le * nEquation;    //! first index
            indexs = le * nEquation;    //! second index
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, le);
            for (int indexI = 0; indexI < nEquation; indexI++)
            {
                for (int indexJ = 0; indexJ < nEquation; indexJ++)
                {
                    for (int indexK = 0; indexK < nEquation; indexK++)
                    {
                        /* dRdq[indexf+indexI][indexs+indexJ] += dFluxdpL[indexI][indexK]*
                                                                dqdcvL[indexK][indexJ]; */
                        dRdq[indexI][colidx + indexJ] += dFluxdpL[indexI][indexK] * dqdcvL[indexK][indexJ];
                    }
                }
            }

            //! GMRESPV , obtain the dqdcv for the right cell
            gas->dPrimitive2dConservative(primRC, gmR, dqdcvR);

            //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
            //! then minus to the corresponding location in dRdq.
            indexf = re * nEquation;    //! first index
            indexs = re * nEquation;    //! second index
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, re);
            for (int indexI = 0; indexI < nEquation; indexI++)
            {
                for (int indexJ = 0; indexJ < nEquation; indexJ++)
                {
                    for (int indexK = 0; indexK < nEquation; indexK++)
                    {
                        /*dRdq[indexf+indexI][indexs+indexJ] -= dFluxdpR[indexI][indexK]*
                                                                dqdcvR[indexK][indexJ];*/
                        dRdq[indexI][colidx + indexJ] -= dFluxdpR[indexI][indexK] * dqdcvR[indexK][indexJ];
                    }
                }
            }

            //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
            //! then add to the corresponding location in dRdq.
            indexf = le * nEquation;    //! first index
            indexs = re * nEquation;    //! second index
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, re);
            for (int indexI = 0; indexI < nEquation; indexI++)
            {
                for (int indexJ = 0; indexJ < nEquation; indexJ++)
                {
                    for (int indexK = 0; indexK < nEquation; indexK++)
                    {
                        /* dRdq[indexf+indexI][indexs+indexJ] += dFluxdpR[indexI][indexK]*
                                                                dqdcvR[indexK][indexJ]; */
                        dRdq[indexI][colidx + indexJ] += dFluxdpR[indexI][indexK] * dqdcvR[indexK][indexJ];
                    }
                }
            }

            if(JacOrder == 2)
            {
                //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
                //! then minus to the corresponding location in dRdq.
                colidx = GetIndexOfBlockJacobianMatrix(AI1st, AJ1st, nEquation, re, le);
                for (int indexI = 0; indexI < nEquation; indexI++)
                {
                    for (int indexJ = 0; indexJ < nEquation; indexJ++)
                    {
                        for (int indexK = 0; indexK < nEquation; indexK++)
                        {
                            dRdq1st[indexI][colidx + indexJ] -= dFluxdpwall[indexI][indexK] * dqdcvL[indexK][indexJ];
                        }
                    }
                }

                //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
                //! then add to the corresponding location in dRdq.
                indexf = le * nEquation;    //! first index
                indexs = le * nEquation;    //! second index
                colidx = GetIndexOfBlockJacobianMatrix(AI1st, AJ1st, nEquation, le, le);
                for (int indexI = 0; indexI < nEquation; indexI++)
                {
                    for (int indexJ = 0; indexJ < nEquation; indexJ++)
                    {
                        for (int indexK = 0; indexK < nEquation; indexK++)
                        {
                            dRdq1st[indexI][colidx + indexJ] += dFluxdpL[indexI][indexK] * dqdcvL[indexK][indexJ];
                        }
                    }
                }

                //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
                //! then minus to the corresponding location in dRdq.
                indexf = re * nEquation;    //! first index
                indexs = re * nEquation;    //! second index
                colidx = GetIndexOfBlockJacobianMatrix(AI1st, AJ1st, nEquation, re, re);
                for (int indexI = 0; indexI < nEquation; indexI++)
                {
                    for (int indexJ = 0; indexJ < nEquation; indexJ++)
                    {
                        for (int indexK = 0; indexK < nEquation; indexK++)
                        {
                            dRdq1st[indexI][colidx + indexJ] -= dFluxdpR[indexI][indexK] * dqdcvR[indexK][indexJ];
                        }
                    }
                }

                //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
                //! then add to the corresponding location in dRdq.
                indexf = le * nEquation;    //! first index
                indexs = re * nEquation;    //! second index
                colidx = GetIndexOfBlockJacobianMatrix(AI1st, AJ1st, nEquation, le, re);
                for (int indexI = 0; indexI < nEquation; indexI++)
                {
                    for (int indexJ = 0; indexJ < nEquation; indexJ++)
                    {
                        for (int indexK = 0; indexK < nEquation; indexK++)
                        {
                            dRdq1st[indexI][colidx + indexJ] += dFluxdpR[indexI][indexK] * dqdcvR[indexK][indexJ];
                        }
                    }
                }
            }
        }
        //! GMRESParallel 
        else if (bcType == PHENGLEI::INTERFACE)
        {
            if(JacOrder == 2)
            {
                 //! obtain the dqdcv for the left cell
                gas->dPrimitive2dConservative(primLC, gmL, dqdcvL);
                int indexf = le * nEquation;
                int dndexs = le * nEquation;
                colidx = GetIndexOfBlockJacobianMatrix(AI1st, AJ1st, nEquation, le, le);
                for (int indexI = 0; indexI < nEquation; indexI ++)
                {
                    for (int indexJ = 0; indexJ < nEquation; indexJ++)
                    {
                        for (int indexK = 0; indexK < nEquation; indexK++)
                        {
                            dRdq1st[indexI][colidx + indexJ] += dFluxdpL[indexI][indexK] * dqdcvL[indexK][indexJ];
                        }
                    }
                }
            }

            if(limiterType == ILMT_FIRST)
            {
                //! obtain the dqdcv for the left cell
                gas->dPrimitive2dConservative(primLC, gmL, dqdcvL);

                int indexf = le * nEquation;
                int indexs = le * nEquation;
                colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, le);
                for (int indexI = 0; indexI < nEquation; indexI++)
                {
                    for (int indexJ = 0; indexJ < nEquation; indexJ++)
                    {
                        for (int indexK = 0; indexK < nEquation; indexK++)
                        {
                            dRdq[indexI][colidx + indexJ] += dFluxdpL[indexI][indexK] * dqdcvL[indexK][indexJ];
                        }
                    }
                }

                indexf = re * nEquation;
                indexs = le * nEquation;
                colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, le);
                for (int indexI = 0; indexI < nEquation; indexI++)
                {
                    for (int indexJ = 0; indexJ < nEquation; indexJ++)
                    {
                        for (int indexK = 0; indexK < nEquation; indexK++)
                        {
                            dRdq[indexI][colidx + indexJ] -= dFluxdpL[indexI][indexK] * dqdcvL[indexK][indexJ];
                        }
                    }
                }
            }
            else if(limiterType ==ILMT_NOLIM)
            {
                                Lsign = qlsign[iFace];
                Rsign = qrsign[iFace];
                int jFace = iFace + localStart; // absolute face index

                if(Lsign > 0)
                {
                    dxL = xfc[jFace] - xcc[le];
                    dyL = yfc[jFace] - ycc[le];
                    dzL = zfc[jFace] - zcc[le];
                }
                else
                {
                    dxL = 0.0;
                    dyL = 0.0;
                    dzL = 0.0;
                }

                if(Rsign > 0)
                {
                    dxR = xfc[jFace] - xcc[re];
                    dyR = yfc[jFace] - ycc[re];
                    dzR = zfc[jFace] - zcc[re];
                }
                else
                {
                    dxR = 0.0;
                    dyR = 0.0;
                    dzR = 0.0;
                }

                //! obtain the dqdcv
                gas->dPrimitive2dConservative(primLC, gmL, dqdcvL);

                int indexf1 = le * nEquation;
                int indexf2 = re * nEquation;
                int indexs  = le * nEquation;
                int colidx1 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, le);
                int colidx2 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, le);
                for (int indexI = 0; indexI < nEquation; indexI++)
                {
                    for (int indexJ = 0; indexJ < nEquation; indexJ++)
                    {
                        for (int indexK = 0; indexK < nEquation; indexK++)
                        {
                            RDouble tmp = dFluxdpL[indexI][indexK] * dqdcvL[indexK][indexJ];
                            dRdq[indexI][colidx1 + indexJ] += tmp;
                            dRdq[indexI][colidx2 + indexJ] -= tmp;
                            // if(PHMPI::GetCurrentProcessorID() == 2 && le == 138)
                            // {
                            //     printf("dFluxdpR %lf, dqdcvR %lf\n", dFluxdpL[indexI][indexK], dqdcvL[indexK][indexJ]);
                            // }
                        }
                    }
                }

                //! considering the perturbation of the gradient in the limiter
                RDouble fnx = nsx[jFace];
                RDouble fny = nsy[jFace];
                RDouble fnz = nsz[jFace];
                RDouble farea = ns[jFace];

                RDouble volume = vol[re];
                dGdqx = -1.0 * 0.5 / volume * farea * fnx;
                dGdqy = -1.0 * 0.5 / volume * farea * fny;
                dGdqz = -1.0 * 0.5 / volume * farea * fnz;
                RDouble dGrad = dGdqx * dxR + dGdqy * dyR + dGdqz * dzR;

                for (int indexI = 0; indexI < nEquation; indexI ++)
                {
                    for (int indexJ = 0; indexJ < nEquation; indexJ++)
                    {
                        for (int indexK = 0; indexK < nEquation; indexK++)
                        {
                            RDouble tmp = dGrad * dFluxdpR[indexI][indexK] * dqdcvL[indexK][indexJ];
                            dRdq[indexI][colidx1 + indexJ] += tmp;
                            dRdq[indexI][colidx2 + indexJ] -= tmp;
                        }
                    }
                }

                for (int index = 0; index < neighborCells[le].size(); index++)
                {
                    int neighborIndex = neighborCells[le][index];
                    if(neighborIndex != re)
                    {
                        int faceIndex = neighborFaces[le][index]; // absolute face index
                        //! judge whether it is the boundary face as well and obtain its boundary type
                        bool isGhostCell = false;
                        UnstructBC *neighborBCRegion = nullptr;
                        int neighborBCType;
                        if (neighborIndex >= nTotalCell)
                        {
                            isGhostCell = true;
                            neighborBCRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[faceIndex]);
                            neighborBCType = neighborBCRegion->GetBCType();
                        }

                        int sign = neighborLR[le][index];
                        RDouble fnx = nsx[faceIndex];
                        RDouble fny = nsy[faceIndex];
                        RDouble fnz = nsz[faceIndex];
                        RDouble farea = ns[faceIndex];

                        RDouble volume = vol[le];
                        for (int m = 0; m < nEquation; m++)
                        {
                            primN[m] = prim[m][neighborIndex];
                        }

                        gas->dPrimitive2dConservative(primN, gama[neighborIndex], dqdcvN);

                        dGdqx = sign * 0.5 / volume * farea * fnx;
                        dGdqy = sign * 0.5 / volume * farea * fny;
                        dGdqz = sign * 0.5 / volume * farea * fnz;

                        dGrad = dGdqx * dxL + dGdqy * dyL + dGdqz * dzL;

                        int indexf1 = le * nEquation;
                        int indexf2 = re * nEquation;
                        int indexs  = neighborIndex * nEquation;
                        int colidx1 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, neighborIndex);
                        int colidx2 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, neighborIndex);
                        for (int indexI = 0; indexI < nEquation; indexI++)
                        {
                            for (int indexJ = 0; indexJ < nEquation; indexJ++)
                            {
                                for (int indexK = 0; indexK < nEquation; indexK++)
                                {
                                    RDouble tmp = dGrad * dFluxdpL[indexI][indexK] * dqdcvN[indexK][indexJ];
                                    dRdq[indexI][colidx1 + indexJ] += tmp;
                                    dRdq[indexI][colidx2 + indexJ] -= tmp;
                                }
                            }
                        }

                        if (isGhostCell && neighborBCType != PHENGLEI::INTERFACE)
                        {
                            int indexre = (neighborIndex - nTotalCell) * nEquation;
                            for (int indexI = 0; indexI < nEquation; indexI++)
                            {
                                for (int indexJ = 0; indexJ < nEquation; indexJ++)
                                {
                                    dDdPlocal[indexI][indexJ] = dDdP[indexI][indexre + indexJ];
                                }
                            }

                            //! HINT: this can be deleted
                            if(sign != 1.0) {
                                TK_Exit::ExceptionExit("sign is fault in Flux_Inviscid.cpp");
                            }

                            //! contribution from the bc
                            //! dFluxdpL = dFluxdpL + dFluxdpR * dDdPlocal
                            for (int indexI = 0; indexI < nEquation; indexI++)
                            {
                                for (int indexJ = 0; indexJ < nEquation; indexJ++)
                                {
                                    dFluxdtmp[indexI][indexJ] = 0.0;
                                    for (int indexK = 0; indexK < nEquation; indexK++)
                                    {
                                        dFluxdtmp[indexI][indexJ] += dGrad * dFluxdpL[indexI][indexK] * dDdPlocal[indexK][indexJ];
                                    }
                                }
                            }
                            
                            gas->dPrimitive2dConservative(primLC, gmL, dqdcvL);
                            colidx1 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, le);
                            colidx2 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, le);
                            for (int indexI = 0; indexI < nEquation; indexI++)
                            {
                                for (int indexJ = 0; indexJ < nEquation; indexJ++)
                                {
                                    for (int indexK = 0; indexK < nEquation; indexK++)
                                    {
                                        RDouble tmp = dFluxdtmp[indexI][indexK] * dqdcvL[indexK][indexJ];
                                        dRdq[indexI][colidx1 + indexJ] += tmp;
                                        dRdq[indexI][colidx2 + indexJ] -= tmp;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        else
        {
            if(JacOrder == 2)
            {
                //! GMRESPV , obtain the dqdcv for the left cell
                gas->dPrimitive2dConservative(primLC, gmL, dqdcvL);

                //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
                //! then minus to the corresponding location in dRdq.
                int indexf = re * nEquation;    //! first index
                int indexs = le * nEquation;    //! second index
                colidx = GetIndexOfBlockJacobianMatrix(AI1st, AJ1st, nEquation, re, le);
                for (int indexI = 0; indexI < nEquation; indexI++)
                {
                    for (int indexJ = 0; indexJ < nEquation; indexJ++)
                    {
                        for (int indexK = 0; indexK < nEquation; indexK++)
                        {
                            dRdq1st[indexI][colidx + indexJ] -= dFluxdpL[indexI][indexK] * dqdcvL[indexK][indexJ];
                        }
                    }
                }

                //! contributions from the bc,
                //! dFluxdpL = dFluxdpL + dFluxdPR*dDdPlocal
                for (int indexI = 0; indexI < nEquation; indexI++)
                {
                    for (int indexJ = 0; indexJ < nEquation; indexJ++)
                    {
                        for (int indexK = 0; indexK < nEquation; indexK++)
                        {
                            dFluxdpL[indexI][indexJ] += dFluxdpR[indexI][indexK] *
                                                        dDdPlocal[indexK][indexJ];
                        }
                    }
                }

                //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
                //! then add to the corresponding location in dRdq.
                indexf = le * nEquation;    //! first index
                indexs = le * nEquation;    //! second index
                colidx = GetIndexOfBlockJacobianMatrix(AI1st, AJ1st, nEquation, le, le);
                for (int indexI = 0; indexI < nEquation; indexI++)
                {
                    for (int indexJ = 0; indexJ < nEquation; indexJ++)
                    {
                        for (int indexK = 0; indexK < nEquation; indexK++)
                        {
                            dRdq1st[indexI][colidx + indexJ] += dFluxdpL[indexI][indexK] * dqdcvL[indexK][indexJ];
                        }
                    }
                }

                //! GMRESPV , obtain the dqdcv for the right cell
                gas->dPrimitive2dConservative(primRC, gmR, dqdcvR);

                //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
                //! then minus to the corresponding location in dRdq.
                indexf = re * nEquation;    //! first index
                indexs = re * nEquation;    //! second index
                colidx = GetIndexOfBlockJacobianMatrix(AI1st, AJ1st, nEquation, re, re);
                for (int indexI = 0; indexI < nEquation; indexI++)
                {
                    for (int indexJ = 0; indexJ < nEquation; indexJ++)
                    {
                        for (int indexK = 0; indexK < nEquation; indexK++)
                        {
                            dRdq1st[indexI][colidx + indexJ] -= dFluxdpR[indexI][indexK] * dqdcvR[indexK][indexJ];
                        }
                    }
                }

                //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
                //! then add to the corresponding location in dRdq.
                indexf = le * nEquation;    //! first index
                indexs = re * nEquation;    //! second index
                colidx = GetIndexOfBlockJacobianMatrix(AI1st, AJ1st, nEquation, le, re);
                for (int indexI = 0; indexI < nEquation; indexI++)
                {
                    for (int indexJ = 0; indexJ < nEquation; indexJ++)
                    {
                        for (int indexK = 0; indexK < nEquation; indexK++)
                        {
                            dRdq1st[indexI][colidx + indexJ] += dFluxdpR[indexI][indexK] * dqdcvR[indexK][indexJ];
                        }
                    }
                }

                if ( limiterType == ILMT_VENCAT && 
                    ( bcType == PHENGLEI::OUTFLOW || bcType == PHENGLEI::INTERFACE ) )
                {
                    for (int m = 0; m < nEquation; m++)
                    {
                        for (int n = 0; n < nEquation; ++ n)
                        {
                            dFluxdpL[m][n] = fluxp[m].dx(n) * area;
                            dFluxdpR[m][n] = fluxp[m].dx(halfIndependentVars + n) * area;
                        }

                    }
                }
                else
                {
                    for (int m = 0; m < nEquation; m++)
                    {
                        for (int n = 0; n < nEquation; n++)
                        {
                            dFluxdpL[m][n] = fluxp[m].dx(n) * area;
                            dFluxdpR[m][n] = fluxp[m].dx(n + nEquation) * area;
                        }
                    }
                }
            }

            if ( limiterType == ILMT_FIRST || 
                 (limiterType == ILMT_VENCAT && 
                    (bcType != PHENGLEI::OUTFLOW && bcType != PHENGLEI::INTERFACE )))
            {
                //! GMRESPV , obtain the dqdcv for the left cell
                gas->dPrimitive2dConservative(primLC, gmL, dqdcvL);

                //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
                //! then minus to the corresponding location in dRdq.
                int indexf = re * nEquation;    //! first index
                int indexs = le * nEquation;    //! second index
                colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, le);
                for (int indexI = 0; indexI < nEquation; indexI++)
                {
                    for (int indexJ = 0; indexJ < nEquation; indexJ++)
                    {
                        for (int indexK = 0; indexK < nEquation; indexK++)
                        {
                            dRdq[indexI][colidx + indexJ] -= dFluxdpL[indexI][indexK] * dqdcvL[indexK][indexJ];
                        }
                    }
                }

                //! contributions from the bc,
                //! dFluxdpL = dFluxdpL + dFluxdPR*dDdPlocal
                for (int indexI = 0; indexI < nEquation; indexI++)
                {
                    for (int indexJ = 0; indexJ < nEquation; indexJ++)
                    {
                        for (int indexK = 0; indexK < nEquation; indexK++)
                        {
                            dFluxdpL[indexI][indexJ] += dFluxdpR[indexI][indexK] *
                                                        dDdPlocal[indexK][indexJ];
                        }
                    }
                }

                //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
                //! then add to the corresponding location in dRdq.
                indexf = le * nEquation;    //! first index
                indexs = le * nEquation;    //! second index
                colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, le);
                for (int indexI = 0; indexI < nEquation; indexI++)
                {
                    for (int indexJ = 0; indexJ < nEquation; indexJ++)
                    {
                        for (int indexK = 0; indexK < nEquation; indexK++)
                        {

                            dRdq[indexI][colidx + indexJ] += dFluxdpL[indexI][indexK] * dqdcvL[indexK][indexJ];
                        }
                    }
                }

                //!z GMRESPV , obtain the dqdcv for the right cell
                gas->dPrimitive2dConservative(primRC, gmR, dqdcvR);

                //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
                //! then minus to the corresponding location in dRdq.
                indexf = re * nEquation;    //! first index
                indexs = re * nEquation;    //! second index
                colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, re);
                for (int indexI = 0; indexI < nEquation; indexI++)
                {
                    for (int indexJ = 0; indexJ < nEquation; indexJ++)
                    {
                        for (int indexK = 0; indexK < nEquation; indexK++)
                        {
                            dRdq[indexI][colidx + indexJ] -= dFluxdpR[indexI][indexK] * dqdcvR[indexK][indexJ];
                        }
                    }
                }

                //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
                //! then add to the corresponding location in dRdq.
                indexf = le * nEquation;    //! first index
                indexs = re * nEquation;    //! second index
                colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, re);
                for (int indexI = 0; indexI < nEquation; indexI++)
                {
                    for (int indexJ = 0; indexJ < nEquation; indexJ++)
                    {
                        for (int indexK = 0; indexK < nEquation; indexK++)
                        {
                            /* dRdq[indexf+indexI][indexs+indexJ] += dFluxdpR[indexI][indexK]*
                                                                    dqdcvR[indexK][indexJ]; */
                            dRdq[indexI][colidx + indexJ] += dFluxdpR[indexI][indexK] * dqdcvR[indexK][indexJ];
                        }
                    }
                }
            }
            else if (limiterType == ILMT_NOLIM)
            {
                Lsign = qlsign[iFace];
                Rsign = qrsign[iFace];
                int jFace = iFace + localStart;    //! absolute face index

                //! GMRESgrad
                if (Lsign > 0)
                {
                    dxL = xfc[jFace] - xcc[le];
                    dyL = yfc[jFace] - ycc[le];
                    dzL = zfc[jFace] - zcc[le];
                }
                else
                {
                    dxL = 0.0;
                    dyL = 0.0;
                    dzL = 0.0;
                }

                dxR = 0.0;
                dyR = 0.0;
                dzR = 0.0;
                //! GMRESPV , obtain the dqdcv for the left cell
                gas->dPrimitive2dConservative(primLC, gmL, dqdcvL);

                //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
                //! then minus to the corresponding location in dRdq.
                int indexf = re * nEquation;    //! first index
                int indexs = le * nEquation;    //! second index
                colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, le);
                for (int indexI = 0; indexI < nEquation; indexI++)
                {
                    for (int indexJ = 0; indexJ < nEquation; indexJ++)
                    {
                        for (int indexK = 0; indexK < nEquation; indexK++)
                        {
                            dRdq[indexI][colidx + indexJ] -= dFluxdpL[indexI][indexK] * dqdcvL[indexK][indexJ];
                        }
                    }
                }

                //! contributions from the bc,
                //! dFluxdpL = dFluxdpL + dFluxdPR*dDdPlocal
                for (int indexI = 0; indexI < nEquation; indexI++)
                {
                    for (int indexJ = 0; indexJ < nEquation; indexJ++)
                    {
                        dFluxdtmp[indexI][indexJ] = 0.0;
                        for (int indexK = 0; indexK < nEquation; indexK++)
                        {
                            dFluxdtmp[indexI][indexJ] += dFluxdpR[indexI][indexK] *
                                                        dDdPlocal[indexK][indexJ];
                        }
                    }
                }

                //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
                //! then add to the corresponding location in dRdq.
                indexf = le * nEquation;    //! first index
                indexs = le * nEquation;    //! second index
                colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, le);
                for (int indexI = 0; indexI < nEquation; indexI++)
                {
                    for (int indexJ = 0; indexJ < nEquation; indexJ++)
                    {
                        for (int indexK = 0; indexK < nEquation; indexK++)
                        {
                            dRdq[indexI][colidx + indexJ] += dFluxdtmp[indexI][indexK] * dqdcvL[indexK][indexJ];
                            dRdq[indexI][colidx + indexJ] += dFluxdpL[indexI][indexK] * dqdcvL[indexK][indexJ];
                        }
                    }
                }

                //! GMRESPV , obtain the dqdcv for the right cell
                gas->dPrimitive2dConservative(primRC, gmR, dqdcvR);

                //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
                //! then minus to the corresponding location in dRdq.
                indexf = re * nEquation;    //! first index
                indexs = re * nEquation;    //! second index
                colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, re);
                for (int indexI = 0; indexI < nEquation; indexI++)
                {
                    for (int indexJ = 0; indexJ < nEquation; indexJ++)
                    {
                        for (int indexK = 0; indexK < nEquation; indexK++)
                        {
                            dRdq[indexI][colidx + indexJ] -= dFluxdpR[indexI][indexK] * dqdcvR[indexK][indexJ];
                        }
                    }
                }

                //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
                //! then add to the corresponding location in dRdq.
                indexf = le * nEquation;    //! first index
                indexs = re * nEquation;    //! second index
                colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, re);
                for (int indexI = 0; indexI < nEquation; indexI++)
                {
                    for (int indexJ = 0; indexJ < nEquation; indexJ++)
                    {
                        for (int indexK = 0; indexK < nEquation; indexK++)
                        {
                            dRdq[indexI][colidx + indexJ] += dFluxdpR[indexI][indexK] * dqdcvR[indexK][indexJ];
                        }
                    }
                }

                //! considering the perturbation of the gradient in the limiter
                RDouble fnx = nsx[jFace];
                RDouble fny = nsy[jFace];
                RDouble fnz = nsz[jFace];
                RDouble farea = ns[jFace];
                RDouble volume = vol[re];

                //! dgradRdqL
                dGdqx = -1.0 * 0.5 / volume * farea * fnx;
                dGdqy = -1.0 * 0.5 / volume * farea * fny;
                dGdqz = -1.0 * 0.5 / volume * farea * fnz;

                //! dQre/dgradR*dgradRdqL
                RDouble dGrad = dGdqx * dxR + dGdqy * dyR + dGdqz * dzR;

                //! convert dGrad by right multiplying the dqdcvL,
                //! then add to the corresponding location in dRdq.
                indexf = le * nEquation; // first index
                indexs = le * nEquation; // second index
                colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, le);
                int colidx2 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, le);
                for (int indexI = 0; indexI < nEquation; indexI++)
                {
                    for (int indexJ = 0; indexJ < nEquation; indexJ++)
                    {
                        for (int indexK = 0; indexK < nEquation; indexK++)
                        {
                            dRdq[indexI][colidx + indexJ] += dGrad * dFluxdpR[indexI][indexK] * dqdcvL[indexK][indexJ];
                            dRdq[indexI][colidx2 + indexJ] -= dGrad * dFluxdpR[indexI][indexK] * dqdcvL[indexK][indexJ];
                        }
                    }
                }

                //! considering the perturbation of the gradient in the limiter
                volume = vol[le];

                //! dgradLdqR
                dGdqx = 1.0 * 0.5 / volume * farea * fnx;
                dGdqy = 1.0 * 0.5 / volume * farea * fny;
                dGdqz = 1.0 * 0.5 / volume * farea * fnz;

                dGrad = dGdqx * dxL + dGdqy * dyL + dGdqz * dzL;

                //! convert dGrad by right multiplying the dqdcvL,
                //! then add to the corresponding location in dRdq.

                for (int indexI = 0; indexI < nEquation; indexI++)
                {
                    for (int indexJ = 0; indexJ < nEquation; indexJ++)
                    {
                        dFluxdtmp[indexI][indexJ] = 0.0;
                        for (int indexK = 0; indexK < nEquation; indexK++)
                        {
                            dFluxdtmp[indexI][indexJ] += dGrad * dFluxdpL[indexI][indexK] *
                                                        dDdPlocal[indexK][indexJ];
                        }
                    }
                }

                colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, re);
                colidx2 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, re);
                int colidx3 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, le);
                for (int indexI = 0; indexI < nEquation; indexI++)
                {
                    for (int indexJ = 0; indexJ < nEquation; indexJ++)
                    {
                        for (int indexK = 0; indexK < nEquation; indexK++)
                        {
                            dRdq[indexI][colidx + indexJ] += dGrad * dFluxdpL[indexI][indexK] * dqdcvR[indexK][indexJ];
                            dRdq[indexI][colidx2 + indexJ] -= dGrad * dFluxdpL[indexI][indexK] * dqdcvR[indexK][indexJ];
                            dRdq[indexI][colidx3 + indexJ] += dFluxdtmp[indexI][indexK] * dqdcvL[indexK][indexJ];
                        }
                    }
                }

                //! the neighbors of the left cells
                for (int index = 0; index < neighborCells[le].size(); index++)
                {
                    int neighborIndex = neighborCells[le][index];
                    // std::cout << "Neighbors of the left cells: " << neighborIndex << " \n";

                    if (neighborIndex != re)
                    {
                        int Faceindex = neighborFaces[le][index];
                        //! judge whether it is the boundary face as well and obtain its boundary type
                        bool isGhostCell = false;
                        UnstructBC *neighborBCRegion = nullptr;
                        int neighborBCType;
                        if (neighborIndex >= nTotalCell)
                        {
                            isGhostCell = true;
                            neighborBCRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[Faceindex]);
                            neighborBCType = neighborBCRegion->GetBCType();
                        }

                        int sign = neighborLR[le][index];
                        RDouble fnx = nsx[Faceindex];
                        RDouble fny = nsy[Faceindex];
                        RDouble fnz = nsz[Faceindex];
                        RDouble farea = ns[Faceindex];
                        RDouble volume = vol[le];

                        for (int m = 0; m < nEquation; m++)
                        {
                            primN[m] = prim[m][neighborIndex];
                        }
                        gas->dPrimitive2dConservative(primN, gama[neighborIndex], dqdcvN);
                        //! dgradLdqN
                        dGdqx = sign * 0.5 / volume * farea * fnx;
                        dGdqy = sign * 0.5 / volume * farea * fny;
                        dGdqz = sign * 0.5 / volume * farea * fnz;

                        //! dQle/dgradL*dgradLdqN
                        dGrad = dGdqx * dxL + dGdqy * dyL + dGdqz * dzL;

                        int indexf = le * nEquation;            //! first index
                        int indexs = neighborIndex * nEquation; //! second index
                        int indexf2 = re * nEquation;
                        int colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, neighborIndex);
                        int colidx2 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, neighborIndex);

                        for (int indexI = 0; indexI < nEquation; indexI++)
                        {
                            for (int indexJ = 0; indexJ < nEquation; indexJ++)
                            {

                                for (int indexK = 0; indexK < nEquation; indexK++)
                                {
                                    //! dRL/dqLNeighbors
                                    dRdq[indexI][colidx + indexJ] += dGrad * dFluxdpL[indexI][indexK] *
                                                                    dqdcvN[indexK][indexJ];

                                    //! dRR/dqLNeighbors
                                    dRdq[indexI][colidx2 + indexJ] -= dGrad * dFluxdpL[indexI][indexK] *
                                                                    dqdcvN[indexK][indexJ];
                                }
                            }
                        }

                        if(isGhostCell && neighborBCType != PHENGLEI::INTERFACE)
                        {
                            int indexre = (neighborIndex - nTotalCell) * nEquation;
                            for (int indexI = 0; indexI < nEquation; indexI++)
                            {
                                for (int indexJ = 0; indexJ < nEquation; indexJ++)
                                {
                                    dDdPlocal[indexI][indexJ] = dDdP[indexI][indexre + indexJ];
                                }
                            }

                            //! HINT: this can be deleted
                            if(sign != 1.0)
                            {
                                TK_Exit::ExceptionExit("sign is falut in Flux_Inviscid.cpp 2");
                            }

                            //! contribution from the bc
                            //! dFluxdpL = dFluxdpL + dFluxdpR * dDdPlocal
                            for (int indexI = 0; indexI < nEquation; indexI++)
                            {
                                for (int indexJ = 0; indexJ < nEquation; indexJ++)
                                {
                                    dFluxdtmp[indexI][indexJ] = 0.0;
                                    for (int indexK = 0; indexK < nEquation; indexK++)
                                    {
                                        dFluxdtmp[indexI][indexJ] += dGrad * dFluxdpL[indexI][indexK] * dDdPlocal[indexK][indexJ];
                                    }
                                }
                            }

                            gas->dPrimitive2dConservative(primLC, gmL, dqdcvL);
                            colidx  = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, le);
                            colidx2 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, le);
                            for (int indexI = 0; indexI < nEquation; indexI++)
                            {
                                for (int indexJ = 0; indexJ < nEquation; indexJ++)
                                {
                                    for (int indexK = 0; indexK < nEquation; indexK++)
                                    {
                                        RDouble tmp = dFluxdtmp[indexI][indexK] * dqdcvL[indexK][indexJ];
                                        dRdq[indexI][colidx  + indexJ] += tmp;
                                        dRdq[indexI][colidx2 + indexJ] -= tmp;
                                    }
                                }
                            }
                    }
                }
              }
            }
            else if (limiterType == ILMT_VENCAT && 
                     ( bcType == PHENGLEI::OUTFLOW || bcType == PHENGLEI::INTERFACE ))
            {
                //! GMRESPV , obtain the dqdcv for the left cell
                gas->dPrimitive2dConservative(primLC, gmL, dqdcvL);

                //! GMRESPV , obtain the dqdcv for the right cell
                gas->dPrimitive2dConservative(primRC, gmR, dqdcvR);

                //! convert dFluxdp to dFluxdcv by right multiplying the dqdcvL/dqdcvR,
                //! then add/minus to the corresponding location in dRdq.

                int colidx1 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, le);
                // int colidx2 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, le);
                int colidx3 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, re);
                int colidx4 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, re);
                for (int indexI = 0; indexI < nEquation; indexI++)
                {
                    for (int indexJ = 0; indexJ < nEquation; indexJ++)
                    {
                        for (int indexK = 0; indexK < nEquation; indexK++)
                        {
                            dRdq[indexI][colidx1 + indexJ] -= dFluxdpL[indexI][indexK] * dqdcvL[indexK][indexJ];

                            dRdq[indexI][colidx3 + indexJ] -= dFluxdpR[indexI][indexK] * dqdcvR[indexK][indexJ];

                            dRdq[indexI][colidx4 + indexJ] += dFluxdpR[indexI][indexK] * dqdcvR[indexK][indexJ];
                        }
                    }
                }

                //! dFluxdpL = dFluxdpL + dFluxdPR*dDdPlocal
                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {
                        for(int indexK=0; indexK < nEquation; indexK++)
                        {
                            dFluxdpL[indexI][indexJ]  += dFluxdpR[indexI][indexK]*
                                                        dDdPlocal[indexK][indexJ];
                        }
                    }
                }

                int colidx2 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, le);
                for (int indexI = 0; indexI < nEquation; indexI++)
                {
                    for (int indexJ = 0; indexJ < nEquation; indexJ++)
                    {
                        for (int indexK = 0; indexK < nEquation; indexK++)
                        {
                            dRdq[indexI][colidx2 + indexJ] += dFluxdpL[indexI][indexK] * dqdcvL[indexK][indexJ];
                        }
                    }
                }

                //! considering the perturbation of the gradient in the flux calculation
                RDouble fnx    = nsx[kFace];
                RDouble fny    = nsy[kFace];
                RDouble fnz    = nsz[kFace];
                RDouble farea  = ns[kFace];
                RDouble volume = vol[re];

                #ifdef RightBCGradient
                //! dgradtmp = dfluxdgradxR*dgradRdqLx + dfluxdgradyR*dgradRdqLy + dfluxdgradzR*dgraddqLz
                //! dRdq[le][le] += dgradtmp*dqdcvL
                gas->dGradient2dPrimitive(primLC,-1,dgraddqx, 'x', fnx, fny, fnz, farea,volume, nEquation);
                gas->dGradient2dPrimitive(primLC,-1,dgraddqy, 'y', fnx, fny, fnz, farea,volume, nEquation);
                gas->dGradient2dPrimitive(primLC,-1,dgraddqz, 'z', fnx, fny, fnz, farea,volume, nEquation);

                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {
                        dGradtmp[indexI][indexJ] = 0;
                        for(int indexK=0; indexK < nEquation; indexK++)
                        {
                            dGradtmp[indexI][indexJ] += dFluxdgradxR[indexI][indexK] * dgraddqx[indexK][indexJ];
                            dGradtmp[indexI][indexJ] += dFluxdgradyR[indexI][indexK] * dgraddqy[indexK][indexJ];
                            dGradtmp[indexI][indexJ] += dFluxdgradzR[indexI][indexK] * dgraddqz[indexK][indexJ];
                        }
                    }
                }

                //! convert dGradtmp by right multiplying the dqdcvL,
                //! then add to the corresponding location in dRdq.
                // indexf = le * nEquation;    //! first index
                // indexs = le * nEquation;    //! second index
                colidx1 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, le);
                colidx2 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, le);

                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {
                        for(int indexK=0; indexK < nEquation; indexK++)
                        {
                            dRdq[indexI][colidx1+indexJ] += dGradtmp[indexI][indexK] * dqdcvL[indexK][indexJ];

                            dRdq[indexI][colidx2+indexJ] -= dGradtmp[indexI][indexK] * dqdcvL[indexK][indexJ];
                        }
                    }
                }
                #endif

                fnx    = nsx[kFace];
                fny    = nsy[kFace];
                fnz    = nsz[kFace];
                farea  = ns[kFace];
                volume = vol[le];
                // dgradtmp = dfluxdgradxL*dgradLdqRx + dfluxdgradyL*dgradLdqRy + dfluxdgradzL*dgradLdqRz
                // dRdq[le][re] += dgradtmp*dqdcvL

                gas->dGradient2dPrimitive(primRC,1,dgraddqx, 'x', fnx, fny, fnz, farea,volume, nEquation);
                gas->dGradient2dPrimitive(primRC,1,dgraddqy, 'y', fnx, fny, fnz, farea,volume, nEquation);
                gas->dGradient2dPrimitive(primRC,1,dgraddqz, 'z', fnx, fny, fnz, farea,volume, nEquation);

                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {
                        dGradtmp[indexI][indexJ] = 0;
                        
                        for(int indexK=0; indexK < nEquation; indexK++)
                        {
                            dGradtmp[indexI][indexJ] += dFluxdgradxL[indexI][indexK]*
                                                                    dgraddqx[indexK][indexJ];
                            dGradtmp[indexI][indexJ] += dFluxdgradyL[indexI][indexK]*
                                                                    dgraddqy[indexK][indexJ];
                            dGradtmp[indexI][indexJ] += dFluxdgradzL[indexI][indexK]*
                                                                    dgraddqz[indexK][indexJ];
                        }
                    }
                }

                //! convert dGradtmp by right multiplying the dqdcvL,
                //! then add to the corresponding location in dRdq.

                colidx1      = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, re);
                colidx2      = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, re);
                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {
                        for(int indexK=0; indexK < nEquation; indexK++)
                        {
                            dRdq[indexI][colidx1+indexJ] += dGradtmp[indexI][indexK]*
                                                                    dqdcvR[indexK][indexJ];

                            dRdq[indexI][colidx2+indexJ] -= dGradtmp[indexI][indexK]*
                                                                    dqdcvR[indexK][indexJ];
                        }
                    }
                }

                //! contributions from the bc, 
                //! dFluxdPR*dDdPlocal
                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {
                        dGradtmpBC[indexI][indexJ] = 0.0;
                        for(int indexK=0; indexK < nEquation; indexK++)
                        {
                            dGradtmpBC[indexI][indexJ]  += dGradtmp[indexI][indexK]*
                                                            dDdPlocal[indexK][indexJ];
                        }
                    }
                }

                colidx1 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, le);
                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {
                        for(int indexK=0; indexK < nEquation; indexK++)
                        {
                            dRdq[indexI][colidx1+indexJ] += dGradtmpBC[indexI][indexK]*
                                                                    dqdcvL[indexK][indexJ];
                        }
                    }
                }

                //! the neighbors of the left cells
                for(int index = 0; index < neighborCells[le].size(); index++)
                {
                    int neighborIndex = neighborCells[le][index];

                    if(neighborIndex != re)
                    {
                        int Faceindex  = neighborFaces[le][index];
                        int sign       = neighborLR[le][index];
                        RDouble fnx    = nsx[Faceindex];
                        RDouble fny    = nsy[Faceindex];
                        RDouble fnz    = nsz[Faceindex];
                        RDouble farea  = ns[Faceindex];
                        RDouble volume = vol[le];

                        for (int m = 0; m < nEquation; m++)
                        {
                            primN[m] = prim[m][neighborIndex];
                        }

                        gas->dPrimitive2dConservative(primN, gama[neighborIndex], dqdcvN);

                        gas->dGradient2dPrimitive(primN,sign,dgraddqx, 'x', fnx, fny, fnz, farea, volume, nEquation);
                        gas->dGradient2dPrimitive(primN,sign,dgraddqy, 'y', fnx, fny, fnz, farea, volume, nEquation);
                        gas->dGradient2dPrimitive(primN,sign,dgraddqz, 'z', fnx, fny, fnz, farea, volume, nEquation);

                        for(int indexI=0; indexI < nEquation; indexI++)
                        {
                            for(int indexJ=0; indexJ < nEquation; indexJ++)
                            {
                                dGradtmp[indexI][indexJ] = 0;

                                for(int indexK=0; indexK < nEquation; indexK++)
                                {
                                    dGradtmp[indexI][indexJ] += dFluxdgradxL[indexI][indexK]*
                                                                            dgraddqx[indexK][indexJ];
                                    dGradtmp[indexI][indexJ] += dFluxdgradyL[indexI][indexK]*
                                                                            dgraddqy[indexK][indexJ];
                                    dGradtmp[indexI][indexJ] += dFluxdgradzL[indexI][indexK]*
                                                                            dgraddqz[indexK][indexJ];
                                }
                            }
                        }

                        colidx1 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, neighborIndex);
                        colidx2 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, neighborIndex);
                        for(int indexI=0; indexI < nEquation; indexI++)
                        {
                            for(int indexJ=0; indexJ < nEquation; indexJ++)
                            {

                                for(int indexK=0; indexK < nEquation; indexK++)
                                {
                                    //! dRL/dqLNeighbors += dRL/dGradL*dGradL/dqLNeighbors ===> dFluxdgradL*dgraddq
                                    dRdq[indexI][colidx1+indexJ] += dGradtmp[indexI][indexK]*
                                                                            dqdcvN[indexK][indexJ];

                                    //! dRR/dqLNeighbors -= dRR/dGradL*dGradL/dqLNeighbors ===> dFluxdgradL*dgraddq
                                    dRdq[indexI][colidx2+indexJ] -= dGradtmp[indexI][indexK]*
                                                                            dqdcvN[indexK][indexJ];
                                }
                            }
                        }

                        //! GMRES2ndBCs
                        if( Faceindex < nBoundFace )
                        {
                            //! contributions from the bc
                            int indexre = (neighborIndex-nTotalCell)*nEquation;
                            for(int indexI=0; indexI < nEquation; indexI++)
                            {
                                for(int indexJ=0; indexJ < nEquation; indexJ++)
                                {
                                    dDdPlocal[indexI][indexJ] = dDdP[indexI][indexre+indexJ];
                                }
                            }

                            //! contributions from the bc, 
                            //! dFluxdPR*dDdPlocal
                            for(int indexI=0; indexI < nEquation; indexI++)
                            {
                                for(int indexJ=0; indexJ < nEquation; indexJ++)
                                {
                                    dGradtmpBC[indexI][indexJ] = 0.0;
                                    for(int indexK=0; indexK < nEquation; indexK++)
                                    {
                                        dGradtmpBC[indexI][indexJ]  += dGradtmp[indexI][indexK]*
                                                                        dDdPlocal[indexK][indexJ];
                                    }
                                }
                            }

                            gas->dPrimitive2dConservative(primLC, gmL, dqdcvL);
                            int colidx1 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, le);
                            int colidx2 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, le);
                            for(int indexI=0; indexI < nEquation; indexI++)
                            {
                                for(int indexJ=0; indexJ < nEquation; indexJ++)
                                {
                                    for(int indexK=0; indexK < nEquation; indexK++)
                                    {
                                        dRdq[indexI][colidx1+indexJ] += dGradtmpBC[indexI][indexK] * dqdcvL[indexK][indexJ];
                                        dRdq[indexI][colidx2+indexJ] -= dGradtmpBC[indexI][indexK] * dqdcvL[indexK][indexJ];
                                    }
                                }
                            }
                        }
                    }
                }

                #ifdef BCRightCells
                //! the neighbors of the right cells
                for(int index = 0; index < neighborCells[re].size(); index++)
                {
                    int neighborIndex = neighborCells[re][index];

                    if(neighborIndex != le)
                    {
                        int Faceindex  = neighborFaces[re][index];
                        int sign       = neighborLR[re][index];
                        RDouble fnx    = nsx[Faceindex];
                        RDouble fny    = nsy[Faceindex];
                        RDouble fnz    = nsz[Faceindex];
                        RDouble farea  = ns[Faceindex];
                        RDouble volume = vol[re];

                        for (int m = 0; m < nEquation; m++)
                        {
                            primN[m] = prim[m][neighborIndex];
                        }

                        gas->dPrimitive2dConservative(primN, gama[neighborIndex], dqdcvN);

                        gas->dGradient2dPrimitive(primN,sign,dgraddqx, 'x', fnx, fny, fnz, farea, volume, nEquation);
                        gas->dGradient2dPrimitive(primN,sign,dgraddqy, 'y', fnx, fny, fnz, farea, volume, nEquation);
                        gas->dGradient2dPrimitive(primN,sign,dgraddqz, 'z', fnx, fny, fnz, farea, volume, nEquation);

                        for(int indexI=0; indexI < nEquation; indexI++)
                        {
                            for(int indexJ=0; indexJ < nEquation; indexJ++)
                            {
                                dGradtmp[indexI][indexJ] = 0;
                                for(int indexK=0; indexK < nEquation; indexK++)
                                {
                                    dGradtmp[indexI][indexJ] += dFluxdgradxR[indexI][indexK]*
                                                                            dgraddqx[indexK][indexJ];
                                    dGradtmp[indexI][indexJ] += dFluxdgradyR[indexI][indexK]*
                                                                            dgraddqy[indexK][indexJ];
                                    dGradtmp[indexI][indexJ] += dFluxdgradzR[indexI][indexK]*
                                                                            dgraddqz[indexK][indexJ];
                                }
                            }
                        }

                        colidx1 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, neighborIndex);
                        colidx2 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, neighborIndex);
                        for(int indexI=0; indexI < nEquation; indexI++)
                        {
                            for(int indexJ=0; indexJ < nEquation; indexJ++)
                            {
                                for(int indexK=0; indexK < nEquation; indexK++)
                                {
                                    //! dRL/dqRNeighbors += dRL/dGradR*dGradR/dqRNeighbors ===> dFluxdgradR*dgraddq
                                    dRdq[indexI][colidx1+indexJ] += dGradtmp[indexI][indexK]*
                                                                            dqdcvN[indexK][indexJ];

                                    //! dRR/dqRNeighbors -= dRR/dGradR*dGradR/dqRNeighbors ===> dFluxdgradR*dgraddq
                                    dRdq[indexI][colidx2+indexJ] -= dGradtmp[indexI][indexK]*
                                                                            dqdcvN[indexK][indexJ];
                                }
                            }
                        }
                    }
                }
                #endif

                for (int nVar = 0; nVar < limitnVariable; ++ nVar)
                {
                    for( int n = 0; n < nEquation; n++ )
                    {
                        for( int m = 0; m < nEquation; m++ )
                        {
                            dFluxdMinMax[n][m] = 0.0;
                        }
                    }

                    if( LimitSignL[nVar] == 1 )
                    {
                        int indexmax = IndexMaxL[nVar];
                        if( indexmax != le )
                        {
                            for (int m = 0; m < nEquation; m++)
                            {
                                dFluxdMinMax[m][nVar]   = dFluxdmaxL[m][nVar];
                                primN[m]                = prim[m][indexmax];
                            }

                            gas->dPrimitive2dConservative(primN, gama[indexmax], dqdcvN);

                            int colidx1     = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, indexmax);
                            int colidx2     = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, indexmax);
                            for(int indexI=0; indexI < nEquation; indexI++)
                            {
                                for(int indexJ=0; indexJ < nEquation; indexJ++)
                                {
                                    for(int indexK=0; indexK < nEquation; indexK++)
                                    {
                                        dRdq[indexI][colidx1+indexJ] += dFluxdMinMax[indexI][indexK]*
                                                                                dqdcvN[indexK][indexJ];

                                        dRdq[indexI][colidx2+indexJ] -= dFluxdMinMax[indexI][indexK]*
                                                                                dqdcvN[indexK][indexJ];
                                    }
                                }
                            }

                            //! GMRES2ndBCs
                            if( indexmax >= nTotalCell )
                            {
                                //! contributions from the bc
                                int indexre = (indexmax-nTotalCell)*nEquation;
                                for(int indexI=0; indexI < nEquation; indexI++)
                                {
                                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                                    {
                                        dDdPlocal[indexI][indexJ] = dDdP[indexI][indexre+indexJ];
                                    }
                                }

                                for(int indexI=0; indexI < nEquation; indexI++)
                                {
                                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                                    {
                                        dGradtmpBC[indexI][indexJ] = 0.0;
                                        for(int indexK=0; indexK < nEquation; indexK++)
                                        {
                                            dGradtmpBC[indexI][indexJ]  += dFluxdMinMax[indexI][indexK]*
                                                                            dDdPlocal[indexK][indexJ];
                                        }
                                    }
                                }

                                gas->dPrimitive2dConservative(primLC, gmL, dqdcvL);
                                int colidx1     = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, le);
                                int colidx2     = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, le);
                                for(int indexI=0; indexI < nEquation; indexI++)
                                {
                                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                                    {
                                        for(int indexK=0; indexK < nEquation; indexK++)
                                        {
                                            dRdq[indexI][colidx1+indexJ] += dGradtmpBC[indexI][indexK] * dqdcvL[indexK][indexJ];
                                            dRdq[indexI][colidx2+indexJ] -= dGradtmpBC[indexI][indexK] * dqdcvL[indexK][indexJ];
                                        }
                                    }
                                }
                            }
                        }
                    }
                    else if ( LimitSignL[nVar] == -1 )
                    {
                        int indexmin = IndexMinL[nVar];
                        if( indexmin != le )
                        {
                            for (int m = 0; m < nEquation; m++)
                            {
                                dFluxdMinMax[m][nVar] = dFluxdminL[m][nVar];
                                primN[m]              = prim[m][indexmin];
                            }

                            gas->dPrimitive2dConservative(primN, gama[indexmin], dqdcvN);

                            int colidx1 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, indexmin);
                            int colidx2 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, indexmin);
                            for(int indexI=0; indexI < nEquation; indexI++)
                            {
                                for(int indexJ=0; indexJ < nEquation; indexJ++)
                                {
                                    for(int indexK=0; indexK < nEquation; indexK++)
                                    {
                                        dRdq[indexI][colidx1+indexJ] += dFluxdMinMax[indexI][indexK]*
                                                                                dqdcvN[indexK][indexJ];

                                        dRdq[indexI][colidx2+indexJ] -= dFluxdMinMax[indexI][indexK]*
                                                                                dqdcvN[indexK][indexJ];
                                    }
                                }
                            }

                            //! GMRES2ndBCs
                            if( indexmin >= nTotalCell )
                            {
                                //! contributions from the bc
                                int indexre = (indexmin-nTotalCell)*nEquation;
                                for(int indexI=0; indexI < nEquation; indexI++)
                                {
                                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                                    {
                                        dDdPlocal[indexI][indexJ] = dDdP[indexI][indexre+indexJ];
                                    }
                                }

                                //! contributions from the bc, 
                                //! dFluxdPR*dDdPlocal
                                for(int indexI=0; indexI < nEquation; indexI++)
                                {
                                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                                    {
                                        dGradtmpBC[indexI][indexJ] = 0.0;
                                        for(int indexK=0; indexK < nEquation; indexK++)
                                        {
                                            dGradtmpBC[indexI][indexJ]  += dFluxdMinMax[indexI][indexK]*
                                                                            dDdPlocal[indexK][indexJ];
                                        }
                                    }
                                }

                                gas->dPrimitive2dConservative(primLC, gmL, dqdcvL);
                                int colidx1     = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, le);
                                int colidx2     = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, le);
                                for(int indexI=0; indexI < nEquation; indexI++)
                                {
                                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                                    {
                                        for(int indexK=0; indexK < nEquation; indexK++)
                                        {
                                            dRdq[indexI][colidx1+indexJ] += dGradtmpBC[indexI][indexK] * dqdcvL[indexK][indexJ];
                                            dRdq[indexI][colidx2+indexJ] -= dGradtmpBC[indexI][indexK] * dqdcvL[indexK][indexJ];
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    for (int iFace = nMid; iFace < nLength; ++iFace)
    {
        int le = leftcellindexofFace[iFace];
        int re = rightcellindexofFace[iFace];
        int colidx;
        RDouble area = areas[iFace];
        int jFace = iFace + localStart;    //! absolute face index

        if ( limiterType == ILMT_VENCAT )
        {
            if (limitVector != 0 || limitModel != 0)
            {
                for (int nVar = 0; nVar < limitnVariable; ++ nVar)
                {
                    IndexMinL[nVar] = le;
                    IndexMaxL[nVar] = le;
                    RDouble pminL = prim[nVar][le];
                    RDouble pmaxL = prim[nVar][le];
                    for ( int index = 0; index < neighborCells[le].size(); index++ )
                    {
                        int neighborIndex     = neighborCells[le][index];
                        int neighborFaceIndex = neighborFaces[le][index];
                        if( neighborFaceIndex >= nBoundFace )
                        {
                            RDouble  primNeighbor    = prim[nVar][neighborIndex];

                            if( pminL > primNeighbor)
                            {
                                pminL           = primNeighbor;
                                IndexMinL[nVar] = neighborIndex;
                            }

                            if( pmaxL < primNeighbor)
                            {
                                pmaxL = primNeighbor;
                                IndexMaxL[nVar] = neighborIndex;
                            }
                        }
                    }

                    IndexMinR[nVar] = re;
                    IndexMaxR[nVar] = re;
                    RDouble pminR = prim[nVar][re];
                    RDouble pmaxR = prim[nVar][re];
                    for ( int index = 0; index < neighborCells[re].size(); index++ )
                    {
                        int neighborIndex       = neighborCells[re][index];
                        int neighborFaceIndex   = neighborFaces[re][index];
                        if( neighborFaceIndex >= nBoundFace )
                        {
                            RDouble  primNeighbor = prim[nVar][neighborIndex];

                            if( pminR > primNeighbor)
                            {
                                pminR           = primNeighbor;
                                IndexMinR[nVar] = neighborIndex;
                            }

                            if( pmaxR < primNeighbor)
                            {
                                pmaxR           = primNeighbor;
                                IndexMaxR[nVar] = neighborIndex;
                            }
                        }
                    }
                }
            }
        }

        //! GMRESVenkat
        if ( limiterType == ILMT_VENCAT )
        {
            for (int iEqn = 0; iEqn < nEquation; ++iEqn)
            {
                primL[iEqn]  = qL[iEqn][iFace];
                primR[iEqn]  = qR[iEqn][iFace];

                primLC[iEqn] = qLC[iEqn][iFace];    //! GMRESPassQC
                primRC[iEqn] = qRC[iEqn][iFace];    //! GMRESPassQC

                dqdxL[iEqn]  = dqdx[iEqn][le];
                dqdyL[iEqn]  = dqdy[iEqn][le];
                dqdzL[iEqn]  = dqdz[iEqn][le];
                dqdxR[iEqn]  = dqdx[iEqn][re];
                dqdyR[iEqn]  = dqdy[iEqn][re];
                dqdzR[iEqn]  = dqdz[iEqn][re];

                minL[iEqn]   = prim[iEqn][IndexMinL[iEqn]];
                maxL[iEqn]   = prim[iEqn][IndexMaxL[iEqn]];
                minR[iEqn]   = prim[iEqn][IndexMinR[iEqn]];
                maxR[iEqn]   = prim[iEqn][IndexMaxR[iEqn]];
            }
        }
        else
        {
            for (int iEqn = 0; iEqn < nEquation; ++iEqn)
            {
                primL[iEqn] = qL[iEqn][iFace];
                primR[iEqn] = qR[iEqn][iFace];

                primLC[iEqn] = qLC[iEqn][iFace];    //! GMRESPassQC
                primRC[iEqn] = qRC[iEqn][iFace];    //! GMRESPassQC
            }
        }

        //! define sequence order of the independent variables
        //! GMRESVenkat
        if ( limiterType == ILMT_VENCAT )
        {
            //! define the sequence of the independent variables
            for (int m = 0; m < nEquation; m++)
            {
                primL[m].diff(m, nIndependentVars);
                primR[m].diff(halfIndependentVars + m, nIndependentVars);
            }

            for (int m = 0; m < nEquation; m++)
            {
                dqdxL[m].diff(nEquation + m, nIndependentVars);
                dqdyL[m].diff(nEquation * 2 + m, nIndependentVars);
                dqdzL[m].diff(nEquation * 3 + m, nIndependentVars);
                minL[m].diff(nEquation * 4 + m, nIndependentVars);
                maxL[m].diff(nEquation * 5 + m, nIndependentVars);
                dqdxR[m].diff(halfIndependentVars + nEquation + m, nIndependentVars);
                dqdyR[m].diff(halfIndependentVars + nEquation * 2 + m, nIndependentVars);
                dqdzR[m].diff(halfIndependentVars + nEquation * 3 + m, nIndependentVars);
                minR[m].diff(halfIndependentVars + nEquation * 4 + m, nIndependentVars);
                maxR[m].diff(halfIndependentVars + nEquation * 5 + m, nIndependentVars);
            }

        }
        else
        {
            for (int iEqn = 0; iEqn < nEquation; iEqn++)
            {
                primL[iEqn].diff(iEqn, 10);
                primR[iEqn].diff(iEqn + nEquation, 10);
            }
        }

         //! GMRESVenkat 
        if ( limiterType == ILMT_VENCAT )
        {
            //! Initialize the limiter with one
            for ( int m = 0; m < limitnVariable; m++ )
            {
                limitL[m] = 1.0;
                limitR[m] = 1.0;
            }

            bool CalLimiterL = true;
            bool CalLimiterR = true;

            //! judge whether the le is adjacent to boundary face
            for (int index = 0; index < neighborCells[le].size(); index++)
            {
                int neighborFaceIndex   = neighborFaces[le][index];
                if( neighborFaceIndex < nBoundFace )
                {
                    UnstructBC* bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[neighborFaceIndex]);
                    int bcType = bcRegion->GetBCType();
                    if (bcType != PHENGLEI::INTERFACE && 
                        bcType != PHENGLEI::SYMMETRY  &&
                        bcType != PHENGLEI::SOLID_SURFACE &&
                        bcType != PHENGLEI::OUTFLOW)
                    {
                        for ( int m = 0; m < limitnVariable; m++ )
                        {
                            limitL[m]     = 0.0;
                            LimitSignL[m] = 0;
                        }

                        CalLimiterL = false;
                    }
                }
            }

            //! judge whether the re is adjacent to boundary face
            for (int index = 0; index < neighborCells[re].size(); index++)
            {
                int neighborFaceIndex   = neighborFaces[re][index];
                if( neighborFaceIndex < nBoundFace )
                {
                    UnstructBC* bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[neighborFaceIndex]);
                    int bcType = bcRegion->GetBCType();
                    if (bcType != PHENGLEI::INTERFACE && 
                        bcType != PHENGLEI::SYMMETRY  &&
                        bcType != PHENGLEI::SOLID_SURFACE &&
                        bcType != PHENGLEI::OUTFLOW)
                    {
                        for ( int m = 0; m < limitnVariable; m++ )
                        {
                            limitR[m]     = 0.0;
                            LimitSignR[m] = 0;
                        }

                        CalLimiterR = false;
                    }
                }
            }

            //! Limiter calculation
            if ( limitVector == SCALAR_LIMITER )
            {
                TK_Exit::ExceptionExit("SCALAR_LIMITER for Venkat is not yet implemented\n");
            }
            else if ( limitVector == VECTOR_LIMITER )
            {
                for (int nVar = 0; nVar < limitnVariable; ++ nVar)
                {
                    if(CalLimiterL)
                    {
                        if( IndexMinL[nVar] == le)
                        {
                            dminL = 0;
                        }
                        else
                        {
                            dminL   =   minL[nVar] - primL[nVar];
                        }

                        if( IndexMaxL[nVar] == le)
                        {
                            dmaxL = 0;
                        }
                        else
                        {
                            dmaxL = maxL[nVar] - primL[nVar];
                        }

                        // for (int index = 0; index < neighborCells[le].size(); index++)
                        // {
                        //     int neighborIndex       = neighborCells[le][index];
                        //     int neighborFaceIndex   = neighborFaces[le][index];
// 
                        //     if( neighborFaceIndex >= nBoundFace )
                        //     {
                        //         ADReal  primNeighbor    = prim[nVar][neighborIndex];
                        //         dminL = MIN(dminL, primNeighbor);
                        //         dmaxL = MAX(dmaxL, primNeighbor);
                        //     }
                        // }

                        ADReal eps_tmp = venkatCoeff * venkatCoeff * venkatCoeff;
                        eps_tmp *= 1.0/pow(pow(averageVolume, 1.0/dimension), nExp);

                        ADReal epsCellL =  eps_tmp*primL[nVar]*primL[nVar];
                        epsCellL        *= pow(pow(vol[le],1.0/dimension),nExp);

                        LimitSignL[nVar]   = 0;
                        for (int index = 0; index < neighborCells[le].size(); index++)
                        {
                            int FaceIndex = neighborFaces[le][index];
                            RDouble dx = xfc[FaceIndex] - xcc[le];
                            RDouble dy = yfc[FaceIndex] - ycc[le];
                            RDouble dz = zfc[FaceIndex] - zcc[le];

                            ADReal dqFace = dqdxL[nVar] * dx + dqdyL[nVar] * dy + dqdzL[nVar] * dz;

                            if( dqFace > SMALL )
                            {
                                    ADReal x2   = dmaxL*dmaxL;
                                    ADReal xy   = dmaxL*dqFace;
                                    ADReal y2   = dqFace*dqFace; 
                                    ADReal temp = (x2 + epsCellL + 2.0 * xy) / (x2 + 2.0 * y2 + xy + epsCellL);

                                    if( temp < limitL[nVar] )
                                    {
                                         limitL[nVar]     = temp;
                                         LimitSignL[nVar] = 1;
                                    }
                            }
                            else if( dqFace < -SMALL )
                            {
                                    ADReal x2   = dminL*dminL;
                                    ADReal xy   = dminL*dqFace;
                                    ADReal y2   = dqFace*dqFace; 
                                    ADReal temp = (x2 + epsCellL + 2.0 * xy) / (x2 + 2.0 * y2 + xy + epsCellL);

                                    if( temp < limitL[nVar] )
                                    {
                                         limitL[nVar]     = temp;
                                         LimitSignL[nVar] = -1;
                                    }
                            }
                        }
                    }

                    if(CalLimiterR)
                    {
                        if( IndexMinR[nVar] == re)
                        {
                            dminR = 0;
                        }
                        else
                        {
                            dminR   =  minR[nVar] - primR[nVar];
                        }

                        if( IndexMaxR[nVar] == re)
                        {
                            dmaxR = 0;
                        }
                        else
                        {
                            dmaxR = maxR[nVar] - primR[nVar];
                        }

                        // for (int index = 0; index < neighborCells[re].size(); index++)
                        // {
                        //     int neighborIndex       = neighborCells[re][index];
                        //     int neighborFaceIndex   = neighborFaces[re][index];
// 
                        //     if( neighborFaceIndex >= nBoundFace )
                        //     {
                        //         ADReal  primNeighbor    = prim[nVar][neighborIndex];
                        //         dminR = MIN(dminR, primNeighbor);
                        //         dmaxR = MAX(dmaxR, primNeighbor);
                        //     }
                        // }

                        ADReal eps_tmp = venkatCoeff * venkatCoeff * venkatCoeff;
                        eps_tmp *= 1.0/pow(pow(averageVolume, 1.0/dimension), nExp);

                        ADReal epsCellR =  eps_tmp*primR[nVar]*primR[nVar];
                        epsCellR        *= pow(pow(vol[re],1.0/dimension),nExp);

                        LimitSignR[nVar] = 0;
                        for (int index = 0; index < neighborCells[re].size(); index++)
                        {
                            int FaceIndex = neighborFaces[re][index];
                            RDouble dx = xfc[FaceIndex] - xcc[re];
                            RDouble dy = yfc[FaceIndex] - ycc[re];
                            RDouble dz = zfc[FaceIndex] - zcc[re];

                            ADReal dqFace = dqdxR[nVar] * dx + dqdyR[nVar] * dy + dqdzR[nVar] * dz;

                            if( dqFace > SMALL )
                            {
                                    ADReal x2   = dmaxR*dmaxR;
                                    ADReal xy   = dmaxR*dqFace;
                                    ADReal y2   = dqFace*dqFace; 
                                    ADReal temp = (x2 + epsCellR + 2.0 * xy) / (x2 + 2.0 * y2 + xy + epsCellR);

                                    if( temp < limitR[nVar] )
                                    {
                                        limitR[nVar]     = temp;
                                        LimitSignR[nVar] = 1;
                                    }
                            }
                            else if( dqFace < -SMALL )
                            {
                                    ADReal x2   = dminR*dminR;
                                    ADReal xy   = dminR*dqFace;
                                    ADReal y2   = dqFace*dqFace; 
                                    ADReal temp = (x2 + epsCellR + 2.0 * xy) / (x2 + 2.0 * y2 + xy + epsCellR);

                                    if( temp < limitR[nVar] )
                                    {
                                        limitR[nVar]     = temp;
                                        LimitSignR[nVar] = -1;
                                    }
                            }
                        }
                    }
                }
            }

            //! Reconstruct the left & right state using the limiter
            if (reconstructMethod == RECONSTRUCT_USING_RESPECTIVE_LIMITER)
            {
                //! reconstruct left state
                RDouble dxL = xfc[jFace] - xcc[le];
                RDouble dyL = yfc[jFace] - ycc[le];
                RDouble dzL = zfc[jFace] - zcc[le];

                for (int m = 0; m < nEquation; ++ m)
                {
                    primTry[m] = primL[m];
                }

                for (int m = 0; m < nEquation; ++ m)
                {
                    if(usingVectorLimiter)
                    {
                        primTry[m] += limitL[m] * (dqdxL[m] * dxL + dqdyL[m] * dyL + dqdzL[m] * dzL);
                    }
                    else
                    {
                        primTry[m] += limitL[0] * (dqdxL[m] * dxL + dqdyL[m] * dyL + dqdzL[m] * dzL);
                    }
                }

                if ( (primTry[IR] > 0.0)&&(primTry[IP] > 0.0) )
                {
                    for (int m = 0; m < nEquation; ++ m)
                    {
                        primL[m] = primTry[m];
                    }
                }

                //! reconstruct right state
                RDouble dxR = xfc[jFace] - xcc[re];
                RDouble dyR = yfc[jFace] - ycc[re];
                RDouble dzR = zfc[jFace] - zcc[re];

                for (int m = 0; m < nEquation; ++ m)
                {
                    primTry[m] = primR[m];
                }

                for (int m = 0; m < nEquation; ++ m)
                {
                    if(usingVectorLimiter)
                    {
                        primTry[m] += limitR[m] * (dqdxR[m] * dxR + dqdyR[m] * dyR + dqdzR[m] * dzR);
                    }
                    else
                    {
                        primTry[m] += limitR[0] * (dqdxR[m] * dxR + dqdyR[m] * dyR + dqdzR[m] * dzR);
                    }
                }

                if ( (primTry[IR] > 0.0)&&(primTry[IP] > 0.0) )
                {
                    for (int m = 0; m < nEquation; ++ m)
                    {
                        primR[m] = primTry[m];
                    }
                }
            }
            else
            {
                //! reconstruct left state
                RDouble dxL = xfc[jFace] - xcc[le];
                RDouble dyL = yfc[jFace] - ycc[le];
                RDouble dzL = zfc[jFace] - zcc[le];

                for (int m = 0; m < nEquation; ++ m)
                {
                    primTry[m] = primL[m];
                }

                for (int m = 0; m < nEquation; ++ m)
                {
                    if(usingVectorLimiter)
                    {
                        primTry[m] += MIN(limitL[m], limitR[m]) * (dqdxL[m] * dxL + dqdyL[m] * dyL + dqdzL[m] * dzL);
                    }
                    else
                    {
                        primTry[m] += MIN(limitL[0], limitR[0]) * (dqdxL[m] * dxL + dqdyL[m] * dyL + dqdzL[m] * dzL);
                    }
                }

                if ( (primTry[IR] > 0.0)&&(primTry[IP] > 0.0) )
                {
                    for (int m = 0; m < nEquation; ++ m)
                    {
                        primL[m] = primTry[m];
                    }
                }

                //! reconstruct right state
                RDouble dxR = xfc[jFace] - xcc[re];
                RDouble dyR = yfc[jFace] - ycc[re];
                RDouble dzR = zfc[jFace] - zcc[re];

                for (int m = 0; m < nEquation; ++ m)
                {
                    primTry[m] = primR[m];
                }

                for (int m = 0; m < nEquation; ++ m)
                {
                    if(usingVectorLimiter)
                    {
                        primTry[m] += MIN(limitR[m], limitL[m]) * (dqdxR[m] * dxR + dqdyR[m] * dyR + dqdzR[m] * dzR);
                    }
                    else
                    {
                        primTry[m] += MIN(limitR[0], limitL[0]) * (dqdxR[m] * dxR + dqdyR[m] * dyR + dqdzR[m] * dzR);
                    }
                }

                if ( (primTry[IR] > 0.0)&&(primTry[IP] > 0.0) )
                {
                    for (int m = 0; m < nEquation; ++ m)
                    {
                        primR[m] = primTry[m];
                    }
                }
            }
        }

        RDouble gmL = gamaL[iFace];
        RDouble gmR = gamaR[iFace];

        //! GMRES3D
        if ( unstructScheme == ISCHEME_GMRES_ROE )
        {
            Cal_GMRES_Roe_Scheme_Flux_PV(primL, primR, iFace, fluxp,
                                        face_proxy, invSchemePara);
        }
        else if ( unstructScheme == ISCHEME_GMRES_Steger )
        {
            Cal_GMRES_Steger_Scheme_Flux_PV(primL, primR, iFace, fluxp,
                                        face_proxy, invSchemePara);
        }
        else
        {
            TK_Exit::ExceptionExit("Error: this inv-scheme is not exist !\n", true);
        }

        flux[IR][iFace] = fluxp[IR].val();
        flux[IU][iFace] = fluxp[IU].val();
        flux[IV][iFace] = fluxp[IV].val();
        flux[IW][iFace] = fluxp[IW].val();
        flux[IP][iFace] = fluxp[IP].val();

        if ( limiterType == ILMT_VENCAT )
        {
            for (int m = 0; m < nEquation; m++)
            {
                for (int n = 0; n < nEquation; ++ n)
                {
                    dFluxdpL[m][n] = fluxp[m].dx(n) * area;
                    dFluxdpR[m][n] = fluxp[m].dx(halfIndependentVars + n) * area;
                }

                for (int n = 0; n < nEquation; n++)
                {
                    dFluxdgradxL[m][n] = fluxp[m].dx(nEquation + n) * area;
                    dFluxdgradyL[m][n] = fluxp[m].dx(nEquation + nEquation + n) * area;
                    dFluxdgradzL[m][n] = fluxp[m].dx(nEquation + 2 * nEquation  + n) * area;
                    dFluxdminL[m][n] = fluxp[m].dx(4 * nEquation  + n) * area;
                    dFluxdmaxL[m][n] = fluxp[m].dx(5 * nEquation  + n) * area;
                    dFluxdgradxR[m][n] = fluxp[m].dx(halfIndependentVars + nEquation + n) * area;
                    dFluxdgradyR[m][n] = fluxp[m].dx(halfIndependentVars + nEquation + nEquation + n) * area;
                    dFluxdgradzR[m][n] = fluxp[m].dx(halfIndependentVars + nEquation + 2 * nEquation + n) * area;
                    dFluxdminR[m][n] = fluxp[m].dx(halfIndependentVars + 4 * nEquation  + n) * area;
                    dFluxdmaxR[m][n] = fluxp[m].dx(halfIndependentVars + 5 * nEquation  + n) * area;
                }
            }
        }
        else
        {
            for (int m = 0; m < nEquation; m++)
            {
                for (int n = 0; n < nEquation; n++)
                {
                    dFluxdpL[m][n] = fluxp[m].dx(n) * area;
                    dFluxdpR[m][n] = fluxp[m].dx(n + nEquation) * area;
                }
            }
        }

//#define outputJac

#ifdef outputJac
        //! output the dFluxdminL dFluxdmaxL
        std::cout << "======================================================================================\n";
        std::cout << "le: " << le <<  " re: " << re << " \n";
        std::cout << "LimitSignL: \n";
        for( int n = 0; n < nEquation; n++ )
        {
            std::cout << LimitSignL[n] << "  ";
        }
        std::cout << "\n";

        std::cout << "IndexMinL: \n";
        for( int n = 0; n < nEquation; n++ )
        {
            std::cout << IndexMinL[n] << "  ";
        }
        std::cout << "\n";

        std::cout << "IndexMaxL: \n";
        for( int n = 0; n < nEquation; n++ )
        {
            std::cout << IndexMaxL[n] << "  ";
        }
        std::cout << "\n";

        std::cout << "dFluxdminL: \n";
        for( int n = 0; n < nEquation; n++ )
        {
            for( int m = 0; m < nEquation; m++ )
            {
                std::cout << dFluxdminL[n][m] << "  ";
            }
            std::cout << "\n";
        }
        std::cout << "\n";

        std::cout << "dFluxdmaxL: \n";
        for( int n = 0; n < nEquation; n++ )
        {
            for( int m = 0; m < nEquation; m++ )
            {
                std::cout << dFluxdmaxL[n][m] << "  ";
            }
            std::cout << "\n";
        }
        std::cout << "\n";

         std::cout << "LimitSignR: \n";
        for( int n = 0; n < nEquation; n++ )
        {
            std::cout << LimitSignR[n] << "  ";
        }
        std::cout << "\n";

        std::cout << "IndexMinR: \n";
        for( int n = 0; n < nEquation; n++ )
        {
            std::cout << IndexMinR[n] << "  ";
        }
        std::cout << "\n";

        std::cout << "IndexMaxR: \n";
        for( int n = 0; n < nEquation; n++ )
        {
            std::cout << IndexMaxR[n] << "  ";
        }
        std::cout << "\n";

        std::cout << "dFluxdminR: \n";
        for( int n = 0; n < nEquation; n++ )
        {
            for( int m = 0; m < nEquation; m++ )
            {
                std::cout << dFluxdminR[n][m] << "  ";
            }
            std::cout << "\n";
        }
        std::cout << "\n";

        std::cout << "dFluxdmaxR: \n";
        for( int n = 0; n < nEquation; n++ )
        {
            for( int m = 0; m < nEquation; m++ )
            {
                std::cout << dFluxdmaxR[n][m] << "  ";
            }
            std::cout << "\n";
        }
        std::cout << "\n";

#endif

        if (JacOrder == 2)
        {
            //! GMRESPV , obtain the dqdcv for the left cell
            gas->dPrimitive2dConservative(primLC, gmL, dqdcvL);

            //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
            //! then minus to the corresponding location in dRdq.
            int indexf = re * nEquation;    //! first index
            int indexs = le * nEquation;    //! second index
            colidx = GetIndexOfBlockJacobianMatrix(AI1st, AJ1st, nEquation, re, le);
            for (int indexI = 0; indexI < nEquation; indexI++)
            {
                for (int indexJ = 0; indexJ < nEquation; indexJ++)
                {
                    for (int indexK = 0; indexK < nEquation; indexK++)
                    {
                        dRdq1st[indexI][colidx + indexJ] -= dFluxdpL[indexI][indexK] * dqdcvL[indexK][indexJ];
                    }
                }
            }

            //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
            //! then add to the corresponding location in dRdq.
            indexf = le * nEquation;    //! first index
            indexs = le * nEquation;    //! second index
            colidx = GetIndexOfBlockJacobianMatrix(AI1st, AJ1st, nEquation, le, le);
            for (int indexI = 0; indexI < nEquation; indexI++)
            {
                for (int indexJ = 0; indexJ < nEquation; indexJ++)
                {
                    for (int indexK = 0; indexK < nEquation; indexK++)
                    {
                        dRdq1st[indexI][colidx + indexJ] += dFluxdpL[indexI][indexK] * dqdcvL[indexK][indexJ];
                    }
                }
            }

            //! GMRESPV , obtain the dqdcv for the right cell
            gas->dPrimitive2dConservative(primRC, gmR, dqdcvR);

            //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
            //! then minus to the corresponding location in dRdq.
            indexf = re * nEquation;    //! first index
            indexs = re * nEquation;    //! second index
            colidx = GetIndexOfBlockJacobianMatrix(AI1st, AJ1st, nEquation, re, re);
            for (int indexI = 0; indexI < nEquation; indexI++)
            {
                for (int indexJ = 0; indexJ < nEquation; indexJ++)
                {
                    for (int indexK = 0; indexK < nEquation; indexK++)
                    {
                        dRdq1st[indexI][colidx + indexJ] -= dFluxdpR[indexI][indexK] * dqdcvR[indexK][indexJ];
                    }
                }
            }

            //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
            //! then add to the corresponding location in dRdq.
            indexf = le * nEquation;    //! first index
            indexs = re * nEquation;    //! second index
            colidx = GetIndexOfBlockJacobianMatrix(AI1st, AJ1st, nEquation, le, re);
            for (int indexI = 0; indexI < nEquation; indexI++)
            {
                for (int indexJ = 0; indexJ < nEquation; indexJ++)
                {
                    for (int indexK = 0; indexK < nEquation; indexK++)
                    {
                        dRdq1st[indexI][colidx + indexJ] += dFluxdpR[indexI][indexK] * dqdcvR[indexK][indexJ];
                    }
                }
            }
        }

        if (limiterType == ILMT_FIRST)
        {
            //! GMRESPV , obtain the dqdcv for the left cell
            gas->dPrimitive2dConservative(primLC, gmL, dqdcvL);

            //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
            //! then minus to the corresponding location in dRdq.
            int indexf = re * nEquation;    //! first index
            int indexs = le * nEquation;    //! second index
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, le);
            for (int indexI = 0; indexI < nEquation; indexI++)
            {
                for (int indexJ = 0; indexJ < nEquation; indexJ++)
                {
                    for (int indexK = 0; indexK < nEquation; indexK++)
                    {
                        dRdq[indexI][colidx + indexJ] -= dFluxdpL[indexI][indexK] * dqdcvL[indexK][indexJ];
                    }
                }
            }

            //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
            //! then add to the corresponding location in dRdq.
            indexf = le * nEquation;    //! first index
            indexs = le * nEquation;    //! second index
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, le);
            for (int indexI = 0; indexI < nEquation; indexI++)
            {
                for (int indexJ = 0; indexJ < nEquation; indexJ++)
                {
                    for (int indexK = 0; indexK < nEquation; indexK++)
                    {
                        dRdq[indexI][colidx + indexJ] += dFluxdpL[indexI][indexK] * dqdcvL[indexK][indexJ];
                    }
                }
            }

            //! GMRESPV , obtain the dqdcv for the right cell
            gas->dPrimitive2dConservative(primRC, gmR, dqdcvR);

            //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
            //! then minus to the corresponding location in dRdq.
            indexf = re * nEquation;    //! first index
            indexs = re * nEquation;    //! second index
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, re);
            for (int indexI = 0; indexI < nEquation; indexI++)
            {
                for (int indexJ = 0; indexJ < nEquation; indexJ++)
                {
                    for (int indexK = 0; indexK < nEquation; indexK++)
                    {
                        dRdq[indexI][colidx + indexJ] -= dFluxdpR[indexI][indexK] * dqdcvR[indexK][indexJ];
                    }
                }
            }

            //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
            //! then add to the corresponding location in dRdq.
            indexf = le * nEquation;    //! first index
            indexs = re * nEquation;    //! second index
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, re);
            for (int indexI = 0; indexI < nEquation; indexI++)
            {
                for (int indexJ = 0; indexJ < nEquation; indexJ++)
                {
                    for (int indexK = 0; indexK < nEquation; indexK++)
                    {
                        dRdq[indexI][colidx + indexJ] += dFluxdpR[indexI][indexK] * dqdcvR[indexK][indexJ];
                    }
                }
            }
        }
        else if (limiterType == ILMT_NOLIM)
        {
            Lsign = qlsign[iFace];
            Rsign = qrsign[iFace];
            int jFace = iFace + localStart;    //! absolute face index

            if (Lsign > 0)
            {
                dxL = xfc[jFace] - xcc[le];
                dyL = yfc[jFace] - ycc[le];
                dzL = zfc[jFace] - zcc[le];
            }
            else
            {
                dxL = 0.0;
                dyL = 0.0;
                dzL = 0.0;
            }

            if (Rsign > 0)
            {
                dxR = xfc[jFace] - xcc[re];
                dyR = yfc[jFace] - ycc[re];
                dzR = zfc[jFace] - zcc[re];
            }
            else
            {
                dxR = 0.0;
                dyR = 0.0;
                dzR = 0.0;
            }

            //! GMRESPV , obtain the dqdcv for the left cell
            gas->dPrimitive2dConservative(primLC, gmL, dqdcvL);

            //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
            //! then minus to the corresponding location in dRdq.
            int indexf = re * nEquation;    //! first index
            int indexs = le * nEquation;    //! second index
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, le);
            for (int indexI = 0; indexI < nEquation; indexI++)
            {
                for (int indexJ = 0; indexJ < nEquation; indexJ++)
                {
                    for (int indexK = 0; indexK < nEquation; indexK++)
                    {
                        dRdq[indexI][colidx + indexJ] -= dFluxdpL[indexI][indexK] * dqdcvL[indexK][indexJ];
                    }
                }
            }

            //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
            //! then add to the corresponding location in dRdq.
            indexf = le * nEquation;    //! first index
            indexs = le * nEquation;    //! second index
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, le);
            for (int indexI = 0; indexI < nEquation; indexI++)
            {
                for (int indexJ = 0; indexJ < nEquation; indexJ++)
                {
                    for (int indexK = 0; indexK < nEquation; indexK++)
                    {
                        dRdq[indexI][colidx + indexJ] += dFluxdpL[indexI][indexK] * dqdcvL[indexK][indexJ];
                    }
                }
            }

            //! GMRESPV , obtain the dqdcv for the right cell
            gas->dPrimitive2dConservative(primRC, gmR, dqdcvR);

            //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
            //! then minus to the corresponding location in dRdq.
            indexf = re * nEquation;    //! first index
            indexs = re * nEquation;    //! second index
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, re);
            for (int indexI = 0; indexI < nEquation; indexI++)
            {
                for (int indexJ = 0; indexJ < nEquation; indexJ++)
                {
                    for (int indexK = 0; indexK < nEquation; indexK++)
                    {
                        dRdq[indexI][colidx + indexJ] -= dFluxdpR[indexI][indexK] * dqdcvR[indexK][indexJ];
                    }
                }
            }

            //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
            //! then add to the corresponding location in dRdq.
            indexf = le * nEquation;    //! first index
            indexs = re * nEquation;    //! second index
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, re);
            for (int indexI = 0; indexI < nEquation; indexI++)
            {
                for (int indexJ = 0; indexJ < nEquation; indexJ++)
                {
                    for (int indexK = 0; indexK < nEquation; indexK++)
                    {
                        dRdq[indexI][colidx + indexJ] += dFluxdpR[indexI][indexK] * dqdcvR[indexK][indexJ];
                    }
                }
            }

            //! considering the perturbation of the gradient in the limiter
            RDouble fnx = nsx[jFace];
            RDouble fny = nsy[jFace];
            RDouble fnz = nsz[jFace];
            RDouble farea = ns[jFace];
            RDouble volume = vol[re];

            //! dgradRdqL
            dGdqx = -1.0 * 0.5 / volume * farea * fnx;
            dGdqy = -1.0 * 0.5 / volume * farea * fny;
            dGdqz = -1.0 * 0.5 / volume * farea * fnz;

            //! dQre/dgradR*dgradRdqL
            RDouble dGrad = dGdqx * dxR + dGdqy * dyR + dGdqz * dzR;

            //! convert dGrad by right multiplying the dqdcvL,
            //! then add to the corresponding location in dRdq.
            indexf = le * nEquation;    //! first index
            indexs = le * nEquation;    //! second index
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, le);
            int colidx2 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, le);
            for (int indexI = 0; indexI < nEquation; indexI++)
            {
                for (int indexJ = 0; indexJ < nEquation; indexJ++)
                {
                    for (int indexK = 0; indexK < nEquation; indexK++)
                    {
                        dRdq[indexI][colidx + indexJ] += dGrad * dFluxdpR[indexI][indexK] * dqdcvL[indexK][indexJ];
                        dRdq[indexI][colidx2 + indexJ] -= dGrad * dFluxdpR[indexI][indexK] * dqdcvL[indexK][indexJ];
                    }
                }
            }

            //! considering the perturbation of the gradient in the limiter
            volume = vol[le];

            //! dgradLdqR
            dGdqx = 1.0 * 0.5 / volume * farea * fnx;
            dGdqy = 1.0 * 0.5 / volume * farea * fny;
            dGdqz = 1.0 * 0.5 / volume * farea * fnz;

            //! dQle/dgradL*dgradLdqR
            dGrad = dGdqx * dxL + dGdqy * dyL + dGdqz * dzL;

            //! convert dGrad by right multiplying the dqdcvL,
            //! then add to the corresponding location in dRdq.

            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, re);
            colidx2 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, re);
            for (int indexI = 0; indexI < nEquation; indexI++)
            {
                for (int indexJ = 0; indexJ < nEquation; indexJ++)
                {
                    for (int indexK = 0; indexK < nEquation; indexK++)
                    {
                        dRdq[indexI][colidx + indexJ] += dGrad * dFluxdpL[indexI][indexK] * dqdcvR[indexK][indexJ];
                        dRdq[indexI][colidx2 + indexJ] -= dGrad * dFluxdpL[indexI][indexK] * dqdcvR[indexK][indexJ];
                    }
                }
            }

            // std::cout << "iFace: " << iFace << " le, re " << le << " , " << re << "\n";
            // the neighbors of the left cells
            for (int index = 0; index < neighborCells[le].size(); index++)
            {
                int neighborIndex = neighborCells[le][index];
                // std::cout << "Neighbors of the left cells: " << neighborIndex << " \n";
                if (neighborIndex != re)
                {
                    int Faceindex = neighborFaces[le][index];
                    //! judge whether it is the boundary face and obtain its boundary type
                    bool isGhostCell = false;
                    UnstructBC *neighborBCRegion = nullptr;
                    int neighborBCType;
                    if(neighborIndex >= nTotalCell)
                    {
                        isGhostCell = true;
                        neighborBCRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[Faceindex]);
                        neighborBCType = neighborBCRegion->GetBCType();
                    }

                    int sign = neighborLR[le][index];
                    RDouble fnx = nsx[Faceindex];
                    RDouble fny = nsy[Faceindex];
                    RDouble fnz = nsz[Faceindex];
                    RDouble farea = ns[Faceindex];
                    RDouble volume = vol[le];

                    for (int m = 0; m < nEquation; m++)
                    {
                        primN[m] = prim[m][neighborIndex];
                    }

                    gas->dPrimitive2dConservative(primN, gama[neighborIndex], dqdcvN);

                    //! dgradLdqN
                    dGdqx = sign * 0.5 / volume * farea * fnx;
                    dGdqy = sign * 0.5 / volume * farea * fny;
                    dGdqz = sign * 0.5 / volume * farea * fnz;

                    //! dQle/dgradL*dgradLdqN
                    dGrad = dGdqx * dxL + dGdqy * dyL + dGdqz * dzL;

                    int indexf = le * nEquation;            // first index
                    int indexs = neighborIndex * nEquation; // second index
                    int indexf2 = re * nEquation;
                    int colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, neighborIndex);
                    int colidx2 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, neighborIndex);
                    // std::cout << "dRle/dqNeighbor: " << le << " " << neighborIndex <<  " " << colidx << " \n";
                    // std::cout << "dRre/dqNeighbor: " << re << " " << neighborIndex <<  " " << colidx2 << " \n";

                    for (int indexI = 0; indexI < nEquation; indexI++)
                    {
                        for (int indexJ = 0; indexJ < nEquation; indexJ++)
                        {
                            for (int indexK = 0; indexK < nEquation; indexK++)
                            {
                                //! dRL/dqLNeighbors
                                dRdq[indexI][colidx + indexJ] += dGrad * dFluxdpL[indexI][indexK] *
                                                                dqdcvN[indexK][indexJ];

                                //! dRR/dqLNeighbors
                                dRdq[indexI][colidx2 + indexJ] -= dGrad * dFluxdpL[indexI][indexK] *
                                                                dqdcvN[indexK][indexJ];
                            }
                        }
                    }

                    if (isGhostCell && neighborBCType != PHENGLEI::INTERFACE)
                    {
                        int indexre = (neighborIndex - nTotalCell) * nEquation;
                        for (int indexI = 0; indexI < nEquation; indexI++)
                        {
                            for (int indexJ = 0; indexJ < nEquation; indexJ++)
                            {
                                dDdPlocal[indexI][indexJ] = dDdP[indexI][indexre + indexJ];
                            }
                        }

                         if (sign != 1.0)
                        {
                            printf("sign %d\n", sign);
                            TK_Exit::ExceptionExit("sign is fault in Flux_Inviscid.cpp 3");
                        }

                        //! contribution from the bc
                        for (int indexI = 0; indexI < nEquation; indexI++)
                        {
                            for (int indexJ = 0; indexJ < nEquation; indexJ++)
                            {
                                dFluxdtmp[indexI][indexJ] = 0.0;
                                for (int indexK = 0; indexK < nEquation; indexK++)
                                {
                                    dFluxdtmp[indexI][indexJ] += dGrad * dFluxdpL[indexI][indexK] * dDdPlocal[indexK][indexJ];
                                }
                            }
                        }

                        gas->dPrimitive2dConservative(primLC, gmL, dqdcvL);
                        colidx  = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, le);
                        colidx2 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, le);
                        for (int indexI = 0; indexI < nEquation; indexI++)
                        {
                            for (int indexJ = 0; indexJ < nEquation; indexJ++)
                            {
                                for (int indexK = 0; indexK < nEquation; indexK++)
                                {
                                    RDouble tmp = dFluxdtmp[indexI][indexK] * dqdcvL[indexK][indexJ];
                                    dRdq[indexI][colidx  + indexJ] += tmp;
                                    dRdq[indexI][colidx2 + indexJ] -= tmp;
                                }
                            }
                        }
                    }
                    
                }
            }

            //! the neighbors of the right cells
            for (int index = 0; index < neighborCells[re].size(); index++)
            {
                int neighborIndex = neighborCells[re][index];
                // std::cout << "Neighbors of the right cells: " << neighborIndex << " \n";

                if (neighborIndex != le)
                {
                    int Faceindex = neighborFaces[re][index];
                    //! judge whether it is the boundary face and obtain its boundary type
                    bool isGhostCell = false;
                    UnstructBC *neighborBCRegion = nullptr;
                    int neighborBCType;
                    if (neighborIndex >= nTotalCell)
                    {
                        isGhostCell = true;
                        neighborBCRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[Faceindex]);
                        neighborBCType = neighborBCRegion->GetBCType();
                    }

                    int sign = neighborLR[re][index];
                    RDouble fnx = nsx[Faceindex];
                    RDouble fny = nsy[Faceindex];
                    RDouble fnz = nsz[Faceindex];
                    RDouble farea = ns[Faceindex];
                    RDouble volume = vol[re];

                    for (int m = 0; m < nEquation; m++)
                    {
                        primN[m] = prim[m][neighborIndex];
                    }

                    gas->dPrimitive2dConservative(primN, gama[neighborIndex], dqdcvN);

                    //! dgradRdqN
                    dGdqx = sign * 0.5 / volume * farea * fnx;
                    dGdqy = sign * 0.5 / volume * farea * fny;
                    dGdqz = sign * 0.5 / volume * farea * fnz;

                    // dQre/dgradR*dgradRdqN
                    dGrad = dGdqx * dxR + dGdqy * dyR + dGdqz * dzR;

                    int indexf = le * nEquation;            //! first index
                    int indexs = neighborIndex * nEquation; //! second index
                    int indexf2 = re * nEquation;
                    int colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, neighborIndex);
                    int colidx2 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, neighborIndex);
                    // std::cout << "dRle/dqNeighbor: " << le << " " << neighborIndex <<  " " << colidx << " \n";
                    // std::cout << "dRre/dqNeighbor: " << re << " " << neighborIndex <<  " " << colidx2 << " \n";
                    for (int indexI = 0; indexI < nEquation; indexI++)
                    {
                        for (int indexJ = 0; indexJ < nEquation; indexJ++)
                        {

                            for (int indexK = 0; indexK < nEquation; indexK++)
                            {
                                //! dRL/dqRNeighbors
                                dRdq[indexI][colidx + indexJ] += dGrad * dFluxdpR[indexI][indexK] *
                                                                dqdcvN[indexK][indexJ];

                                //! dRR/dqRNeighbors
                                dRdq[indexI][colidx2 + indexJ] -= dGrad * dFluxdpR[indexI][indexK] *
                                                                dqdcvN[indexK][indexJ];
                            }
                        }
                    }

                    if (isGhostCell && neighborBCType != PHENGLEI::INTERFACE)
                    {
                        int indexre = (neighborIndex - nTotalCell) * nEquation;
                        for (int indexI = 0; indexI < nEquation; indexI++)
                        {
                            for (int indexJ = 0; indexJ < nEquation; indexJ++)
                            {
                                dDdPlocal[indexI][indexJ] = dDdP[indexI][indexre + indexJ];
                            }
                        }

                        //! HINT: this can be deleted
                        if (sign != 1.0){
                            TK_Exit::ExceptionExit("sign is fault in Flux_Inviscid.cpp 4");
                        }

                        //! contribution from the bc
                        //! dFluxdpL = dFluxdpL + dFluxdpR * dDdPlocal
                        for (int indexI = 0; indexI < nEquation; indexI++)
                        {
                            for (int indexJ = 0; indexJ < nEquation; indexJ++)
                            {
                                dFluxdtmp[indexI][indexJ] = 0.0;
                                for (int indexK = 0; indexK < nEquation; indexK++)
                                {
                                    dFluxdtmp[indexI][indexJ] += dGrad * dFluxdpR[indexI][indexK] * dDdPlocal[indexK][indexJ];
                                }
                            }
                        }

                        gas->dPrimitive2dConservative(primRC, gmR, dqdcvR);
                        colidx  = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, re);
                        colidx2 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, re);
                        for (int indexI = 0; indexI < nEquation; indexI++)
                        {
                            for (int indexJ = 0; indexJ < nEquation; indexJ++)
                            {
                                for (int indexK = 0; indexK < nEquation; indexK++)
                                {
                                    RDouble tmp = dFluxdtmp[indexI][indexK] * dqdcvR[indexK][indexJ];
                                    dRdq[indexI][colidx  + indexJ] += tmp;
                                    dRdq[indexI][colidx2 + indexJ] -= tmp;
                                }
                            }
                        }
                    }
                }
            }

            #ifdef BC_effect_to_Interior
            //! judge whether le is in the left cell lists of the BC
            vector<int>::iterator result = find(BCLeftCells.begin(), BCLeftCells.end(), le);
            if (result != BCLeftCells.end())
            {
                int dist = distance(BCLeftCells.begin(), result);
                int rBCIndex = BCRightCells[dist];
                int FBCIndex = BCFaces[dist];

                //! contributions from the bc
                int indexre = (rBCIndex - nTotalCell) * nEquation;
                for (int indexI = 0; indexI < nEquation; indexI++)
                {
                    for (int indexJ = 0; indexJ < nEquation; indexJ++)
                    {
                        dDdPlocal[indexI][indexJ] = dDdP[indexI][indexre + indexJ];
                    }
                }

                RDouble fnx = nsx[FBCIndex];
                RDouble fny = nsy[FBCIndex];
                RDouble fnz = nsz[FBCIndex];
                RDouble farea = ns[FBCIndex];
                RDouble volume = vol[le];

                //! dgradLdqN
                dGdqx = 1.0 * 0.5 / volume * farea * fnx;
                dGdqy = 1.0 * 0.5 / volume * farea * fny;
                dGdqz = 1.0 * 0.5 / volume * farea * fnz;
                //! dQle/dgradL*dgradLdqN
                dGrad = dGdqx * dxL + dGdqy * dyL + dGdqz * dzL;

                //! contributions from the bc,
                //! dFluxdpL = dFluxdpL + dFluxdPR*dDdPlocal

                for (int indexI = 0; indexI < nEquation; indexI++)
                {
                    for (int indexJ = 0; indexJ < nEquation; indexJ++)
                    {
                        dFluxdtmp[indexI][indexJ] = 0.0;
                        for (int indexK = 0; indexK < nEquation; indexK++)
                        {
                            dFluxdtmp[indexI][indexJ] += dGrad * dFluxdpL[indexI][indexK] *
                                                        dDdPlocal[indexK][indexJ];
                        }
                    }
                }

                gas->dPrimitive2dConservative(primLC, gmL, dqdcvL);
                int colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, le);
                int colidx2 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, le);
                for (int indexI = 0; indexI < nEquation; indexI++)
                {
                    for (int indexJ = 0; indexJ < nEquation; indexJ++)
                    {
                        for (int indexK = 0; indexK < nEquation; indexK++)
                        {
                            dRdq[indexI][colidx + indexJ] += dFluxdtmp[indexI][indexK] * dqdcvL[indexK][indexJ];
                            dRdq[indexI][colidx2 + indexJ] -= dFluxdtmp[indexI][indexK] * dqdcvL[indexK][indexJ];
                        }
                    }
                }
            }

            //! judge whether re is in the left cell lists of the BC
            result = find(BCLeftCells.begin(), BCLeftCells.end(), re);
            if (result != BCLeftCells.end())
            {
                int dist = distance(BCLeftCells.begin(), result);
                int rBCIndex = BCRightCells[dist];
                int FBCIndex = BCFaces[dist];

                // contributions from the bc
                int indexre = (rBCIndex - nTotalCell) * nEquation;
                for (int indexI = 0; indexI < nEquation; indexI++)
                {
                    for (int indexJ = 0; indexJ < nEquation; indexJ++)
                    {
                        dDdPlocal[indexI][indexJ] = dDdP[indexI][indexre + indexJ];
                    }
                }

                RDouble fnx = nsx[FBCIndex];
                RDouble fny = nsy[FBCIndex];
                RDouble fnz = nsz[FBCIndex];
                RDouble farea = ns[FBCIndex];
                RDouble volume = vol[re];

                //! dgradLdqN
                dGdqx = 1.0 * 0.5 / volume * farea * fnx;
                dGdqy = 1.0 * 0.5 / volume * farea * fny;
                dGdqz = 1.0 * 0.5 / volume * farea * fnz;
                // dQle/dgradL*dgradLdqN
                dGrad = dGdqx * dxR + dGdqy * dyR + dGdqz * dzR;

                //! contributions from the bc,
                //! dFluxdpL = dFluxdpL + dFluxdPR*dDdPlocal

                for (int indexI = 0; indexI < nEquation; indexI++)
                {
                    for (int indexJ = 0; indexJ < nEquation; indexJ++)
                    {
                        dFluxdtmp[indexI][indexJ] = 0.0;
                        for (int indexK = 0; indexK < nEquation; indexK++)
                        {
                            dFluxdtmp[indexI][indexJ] += dGrad * dFluxdpR[indexI][indexK] *
                                                        dDdPlocal[indexK][indexJ];
                        }
                    }
                }

                gas->dPrimitive2dConservative(primRC, gmR, dqdcvR);
                int colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, re);
                int colidx2 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, re);
                for (int indexI = 0; indexI < nEquation; indexI++)
                {
                    for (int indexJ = 0; indexJ < nEquation; indexJ++)
                    {
                        for (int indexK = 0; indexK < nEquation; indexK++)
                        {
                            dRdq[indexI][colidx + indexJ] += dFluxdtmp[indexI][indexK] * dqdcvR[indexK][indexJ];
                            dRdq[indexI][colidx2 + indexJ] -= dFluxdtmp[indexI][indexK] * dqdcvR[indexK][indexJ];
                        }
                    }
                }
            }
            #endif
        }
        else if (limiterType == ILMT_VENCAT)
        {
            //! GMRESPV , obtain the dqdcv for the left cell
            gas->dPrimitive2dConservative(primLC, gmL, dqdcvL);

            //! GMRESPV , obtain the dqdcv for the right cell
            gas->dPrimitive2dConservative(primRC, gmR, dqdcvR);

            //! convert dFluxdp to dFluxdcv by right multiplying the dqdcvL/dqdcvR,
            //! then add/minus to the corresponding location in dRdq.

            int colidx1 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, le);
            int colidx2 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, le);
            int colidx3 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, re);
            int colidx4 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, re);
            for (int indexI = 0; indexI < nEquation; indexI++)
            {
                for (int indexJ = 0; indexJ < nEquation; indexJ++)
                {
                    for (int indexK = 0; indexK < nEquation; indexK++)
                    {
                        dRdq[indexI][colidx1 + indexJ] -= dFluxdpL[indexI][indexK] * dqdcvL[indexK][indexJ];

                        dRdq[indexI][colidx2 + indexJ] += dFluxdpL[indexI][indexK] * dqdcvL[indexK][indexJ];

                        dRdq[indexI][colidx3 + indexJ] -= dFluxdpR[indexI][indexK] * dqdcvR[indexK][indexJ];

                        dRdq[indexI][colidx4 + indexJ] += dFluxdpR[indexI][indexK] * dqdcvR[indexK][indexJ];
                    }
                }
            }

            //! considering the perturbation of the gradient in the flux calculation
            RDouble fnx    = nsx[jFace];
            RDouble fny    = nsy[jFace];
            RDouble fnz    = nsz[jFace];
            RDouble farea  = ns[jFace];
            RDouble volume = vol[re];
            //! dgradtmp = dfluxdgradxR*dgradRdqLx + dfluxdgradyR*dgradRdqLy + dfluxdgradzR*dgraddqLz
            //! dRdq[le][le] += dgradtmp*dqdcvL
            gas->dGradient2dPrimitive(primLC,-1,dgraddqx, 'x', fnx, fny, fnz, farea, volume, nEquation);
            gas->dGradient2dPrimitive(primLC,-1,dgraddqy, 'y', fnx, fny, fnz, farea, volume, nEquation);
            gas->dGradient2dPrimitive(primLC,-1,dgraddqz, 'z', fnx, fny, fnz, farea, volume, nEquation);

            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    dGradtmp[indexI][indexJ] = 0;
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        dGradtmp[indexI][indexJ] += dFluxdgradxR[indexI][indexK] * dgraddqx[indexK][indexJ];
                        dGradtmp[indexI][indexJ] += dFluxdgradyR[indexI][indexK] * dgraddqy[indexK][indexJ];
                        dGradtmp[indexI][indexJ] += dFluxdgradzR[indexI][indexK] * dgraddqz[indexK][indexJ];
                    }
                }
            }

            //! convert dGradtmp by right multiplying the dqdcvL,
            //! then add to the corresponding location in dRdq.
            //! indexf = le * nEquation; // first index
            //! indexs = le * nEquation; // second index
            colidx1 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, le);
            colidx2 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, le);

            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        dRdq[indexI][colidx1+indexJ] += dGradtmp[indexI][indexK] * dqdcvL[indexK][indexJ];

                        dRdq[indexI][colidx2+indexJ] -= dGradtmp[indexI][indexK] * dqdcvL[indexK][indexJ];
                    }
                }
            }

            fnx    = nsx[jFace];
            fny    = nsy[jFace];
            fnz    = nsz[jFace];
            farea  = ns[jFace];
            volume = vol[le];
            // dgradtmp = dfluxdgradxL*dgradLdqRx + dfluxdgradyL*dgradLdqRy + dfluxdgradzL*dgradLdqRz
            // dRdq[le][re] += dgradtmp*dqdcvL

            gas->dGradient2dPrimitive(primRC,1,dgraddqx, 'x', fnx, fny, fnz, farea,volume, nEquation);
            gas->dGradient2dPrimitive(primRC,1,dgraddqy, 'y', fnx, fny, fnz, farea,volume, nEquation);
            gas->dGradient2dPrimitive(primRC,1,dgraddqz, 'z', fnx, fny, fnz, farea,volume, nEquation);

            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    dGradtmp[indexI][indexJ] = 0;
                    
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        dGradtmp[indexI][indexJ] += dFluxdgradxL[indexI][indexK]*
                                                                dgraddqx[indexK][indexJ];
                        dGradtmp[indexI][indexJ] += dFluxdgradyL[indexI][indexK]*
                                                                dgraddqy[indexK][indexJ];
                        dGradtmp[indexI][indexJ] += dFluxdgradzL[indexI][indexK]*
                                                                dgraddqz[indexK][indexJ];
                    }
                }
            }

            //! convert dGradtmp by right multiplying the dqdcvL,
            //! then add to the corresponding location in dRdq.

            colidx1 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, re);
            colidx2 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, re);
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        
                        dRdq[indexI][colidx1+indexJ] += dGradtmp[indexI][indexK]*
                                                                dqdcvR[indexK][indexJ];

                        dRdq[indexI][colidx2+indexJ] -= dGradtmp[indexI][indexK]*
                                                                dqdcvR[indexK][indexJ];
                    }
                }
            }

            //! the neighbors of the left cells
            for(int index = 0; index < neighborCells[le].size(); index++)
            {
                int neighborIndex = neighborCells[le][index];

                if(neighborIndex != re)
                {
                    int Faceindex  = neighborFaces[le][index];
                    int sign       = neighborLR[le][index];
                    RDouble fnx    = nsx[Faceindex];
                    RDouble fny    = nsy[Faceindex];
                    RDouble fnz    = nsz[Faceindex];
                    RDouble farea  = ns[Faceindex];
                    RDouble volume = vol[le];

                    for (int m = 0; m < nEquation; m++)
                    {
                        primN[m] = prim[m][neighborIndex];
                    }

                    gas->dPrimitive2dConservative(primN, gama[neighborIndex], dqdcvN);

                    gas->dGradient2dPrimitive(primN,sign,dgraddqx, 'x', fnx, fny, fnz, farea, volume, nEquation);
                    gas->dGradient2dPrimitive(primN,sign,dgraddqy, 'y', fnx, fny, fnz, farea, volume, nEquation);
                    gas->dGradient2dPrimitive(primN,sign,dgraddqz, 'z', fnx, fny, fnz, farea, volume, nEquation);

                    for(int indexI=0; indexI < nEquation; indexI++)
                    {
                        for(int indexJ=0; indexJ < nEquation; indexJ++)
                        {
                            dGradtmp[indexI][indexJ] = 0;

                            for(int indexK=0; indexK < nEquation; indexK++)
                            {
                                dGradtmp[indexI][indexJ] += dFluxdgradxL[indexI][indexK]*
                                                                        dgraddqx[indexK][indexJ];
                                dGradtmp[indexI][indexJ] += dFluxdgradyL[indexI][indexK]*
                                                                        dgraddqy[indexK][indexJ];
                                dGradtmp[indexI][indexJ] += dFluxdgradzL[indexI][indexK]*
                                                                        dgraddqz[indexK][indexJ];
                            }
                        }
                    }

                    colidx1 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, neighborIndex);
                    colidx2 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, neighborIndex);
                    for(int indexI=0; indexI < nEquation; indexI++)
                    {
                        for(int indexJ=0; indexJ < nEquation; indexJ++)
                        {
                            for(int indexK=0; indexK < nEquation; indexK++)
                            {
                                //! dRL/dqLNeighbors += dRL/dGradL*dGradL/dqLNeighbors ===> dFluxdgradL*dgraddq
                                dRdq[indexI][colidx1+indexJ] += dGradtmp[indexI][indexK]*
                                                                        dqdcvN[indexK][indexJ];

                                //! dRR/dqLNeighbors -= dRR/dGradL*dGradL/dqLNeighbors ===> dFluxdgradL*dgraddq
                                dRdq[indexI][colidx2+indexJ] -= dGradtmp[indexI][indexK]*
                                                                        dqdcvN[indexK][indexJ];
                            }
                        }
                    }

                    //! GMRES2ndBCs
                    if( Faceindex < nBoundFace )
                    {
                        //! contributions from the bc
                        int indexre = (neighborIndex-nTotalCell)*nEquation;
                        for(int indexI=0; indexI < nEquation; indexI++)
                        {
                            for(int indexJ=0; indexJ < nEquation; indexJ++)
                            {
                                dDdPlocal[indexI][indexJ] = dDdP[indexI][indexre+indexJ];
                            }
                        }

                        //! contributions from the bc, 
                        //! dFluxdPR*dDdPlocal
                        for(int indexI=0; indexI < nEquation; indexI++)
                        {
                            for(int indexJ=0; indexJ < nEquation; indexJ++)
                            {
                                dGradtmpBC[indexI][indexJ] = 0.0;
                                for(int indexK=0; indexK < nEquation; indexK++)
                                {
                                    dGradtmpBC[indexI][indexJ]  += dGradtmp[indexI][indexK]*
                                                                    dDdPlocal[indexK][indexJ];
                                }
                            }
                        }

                        gas->dPrimitive2dConservative(primLC, gmL, dqdcvL);
                        int colidx1     = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, le);
                        int colidx2     = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, le);
                        for(int indexI=0; indexI < nEquation; indexI++)
                        {
                            for(int indexJ=0; indexJ < nEquation; indexJ++)
                            {
                                for(int indexK=0; indexK < nEquation; indexK++)
                                {
                                    dRdq[indexI][colidx1+indexJ] += dGradtmpBC[indexI][indexK] * dqdcvL[indexK][indexJ];
                                    dRdq[indexI][colidx2+indexJ] -= dGradtmpBC[indexI][indexK] * dqdcvL[indexK][indexJ];
                                }
                            }
                        }
                    }
                }
            }

            //! the neighbors of the right cells
            for(int index = 0; index < neighborCells[re].size(); index++)
            {
                int neighborIndex = neighborCells[re][index];

                if(neighborIndex != le)
                {
                    int Faceindex  = neighborFaces[re][index];
                    int sign       = neighborLR[re][index];
                    RDouble fnx    = nsx[Faceindex];
                    RDouble fny    = nsy[Faceindex];
                    RDouble fnz    = nsz[Faceindex];
                    RDouble farea  = ns[Faceindex];
                    RDouble volume = vol[re];

                    for (int m = 0; m < nEquation; m++)
                    {
                        primN[m] = prim[m][neighborIndex];
                    }

                    gas->dPrimitive2dConservative(primN, gama[neighborIndex], dqdcvN);

                    gas->dGradient2dPrimitive(primN,sign,dgraddqx, 'x', fnx, fny, fnz, farea, volume, nEquation);
                    gas->dGradient2dPrimitive(primN,sign,dgraddqy, 'y', fnx, fny, fnz, farea, volume, nEquation);
                    gas->dGradient2dPrimitive(primN,sign,dgraddqz, 'z', fnx, fny, fnz, farea, volume, nEquation);

                    for(int indexI=0; indexI < nEquation; indexI++)
                    {
                        for(int indexJ=0; indexJ < nEquation; indexJ++)
                        {
                            dGradtmp[indexI][indexJ] = 0;
                            for(int indexK=0; indexK < nEquation; indexK++)
                            {
                                dGradtmp[indexI][indexJ] += dFluxdgradxR[indexI][indexK]*
                                                                        dgraddqx[indexK][indexJ];
                                dGradtmp[indexI][indexJ] += dFluxdgradyR[indexI][indexK]*
                                                                        dgraddqy[indexK][indexJ];
                                dGradtmp[indexI][indexJ] += dFluxdgradzR[indexI][indexK]*
                                                                        dgraddqz[indexK][indexJ];
                            }
                        }
                    }

                    colidx1 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, neighborIndex);
                    colidx2 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, neighborIndex);
                    for(int indexI=0; indexI < nEquation; indexI++)
                    {
                        for(int indexJ=0; indexJ < nEquation; indexJ++)
                        {
                            for(int indexK=0; indexK < nEquation; indexK++)
                            {
                                //! dRL/dqRNeighbors += dRL/dGradR*dGradR/dqRNeighbors ===> dFluxdgradR*dgraddq
                                dRdq[indexI][colidx1+indexJ] += dGradtmp[indexI][indexK]*
                                                                        dqdcvN[indexK][indexJ];

                                //! dRR/dqRNeighbors -= dRR/dGradR*dGradR/dqRNeighbors ===> dFluxdgradR*dgraddq
                                dRdq[indexI][colidx2+indexJ] -= dGradtmp[indexI][indexK]*
                                                                        dqdcvN[indexK][indexJ];
                            }
                        }
                    }

                    //! GMRES2ndBCs
                    if ( Faceindex < nBoundFace )
                    {
                        //! contributions from the bc
                        int indexre = (neighborIndex-nTotalCell)*nEquation;
                        for(int indexI=0; indexI < nEquation; indexI++)
                        {
                            for(int indexJ=0; indexJ < nEquation; indexJ++)
                            {
                                dDdPlocal[indexI][indexJ] = dDdP[indexI][indexre+indexJ];
                            }
                        }

                        //! contributions from the bc, 
                        //! dFluxdPR*dDdPlocal
                        for(int indexI=0; indexI < nEquation; indexI++)
                        {
                            for(int indexJ=0; indexJ < nEquation; indexJ++)
                            {
                                dGradtmpBC[indexI][indexJ] = 0.0;
                                for(int indexK=0; indexK < nEquation; indexK++)
                                {
                                    dGradtmpBC[indexI][indexJ]  += dGradtmp[indexI][indexK]*
                                                                    dDdPlocal[indexK][indexJ];
                                }
                            }
                        }

                        gas->dPrimitive2dConservative(primRC, gmR, dqdcvR);
                        int colidx1     = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, re);
                        int colidx2     = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, re);
                        for(int indexI=0; indexI < nEquation; indexI++)
                        {
                            for(int indexJ=0; indexJ < nEquation; indexJ++)
                            {
                                for(int indexK=0; indexK < nEquation; indexK++)
                                {
                                    dRdq[indexI][colidx1+indexJ] += dGradtmpBC[indexI][indexK] * dqdcvR[indexK][indexJ];
                                    dRdq[indexI][colidx2+indexJ] -= dGradtmpBC[indexI][indexK] * dqdcvR[indexK][indexJ];
                                }
                            }
                        }
                    }
                }
            }

            #ifdef BC_effect_to_Interior
            //! judge whether le is in the left cell lists of the BC 
            vector<int>::iterator result = find(BCLeftCells.begin(), BCLeftCells.end(), le); 
            if(result != BCLeftCells.end())
            {
                int dist     = distance(BCLeftCells.begin(),result);
                int rBCIndex = BCRightCells[dist];
                int FBCIndex = BCFaces[dist];
                //! contributions from the bc
                int indexre = (rBCIndex-nTotalCell)*nEquation;
                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {
                        dDdPlocal[indexI][indexJ] = dDdP[indexI][indexre+indexJ];
                    }
                }

                RDouble fnx    = nsx[FBCIndex];
                RDouble fny    = nsy[FBCIndex];
                RDouble fnz    = nsz[FBCIndex];
                RDouble farea  = ns[FBCIndex];
                RDouble volume = vol[le];

                for (int m = 0; m < nEquation; m++)
                {
                    primN[m] = prim[m][rBCIndex];
                }

                gas->dGradient2dPrimitive(primN,1,dgraddqx, 'x', fnx, fny, fnz, farea, volume, nEquation);
                gas->dGradient2dPrimitive(primN,1,dgraddqy, 'y', fnx, fny, fnz, farea, volume, nEquation);
                gas->dGradient2dPrimitive(primN,1,dgraddqz, 'z', fnx, fny, fnz, farea, volume, nEquation);

                //! contributions from the bc, 
                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {
                        dGradtmp[indexI][indexJ] = 0;

                        for(int indexK=0; indexK < nEquation; indexK++)
                        {
                            dGradtmp[indexI][indexJ] += dFluxdgradxL[indexI][indexK]*
                                                                    dgraddqx[indexK][indexJ];
                            dGradtmp[indexI][indexJ] += dFluxdgradyL[indexI][indexK]*
                                                                    dgraddqy[indexK][indexJ];
                            dGradtmp[indexI][indexJ] += dFluxdgradzL[indexI][indexK]*
                                                                    dgraddqz[indexK][indexJ];
                        }
                    }
                }

                //! contributions from the bc, 
                //! dFluxdPR*dDdPlocal
                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {
                        dGradtmpBC[indexI][indexJ] = 0.0;
                        for(int indexK=0; indexK < nEquation; indexK++)
                        {
                            dGradtmpBC[indexI][indexJ]  += dGradtmp[indexI][indexK]*
                                                            dDdPlocal[indexK][indexJ];
                        }
                    }
                }

                gas->dPrimitive2dConservative(primLC, gmL, dqdcvL);
                int colidx1 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, le);
                int colidx2 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, le);
                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {
                        for(int indexK=0; indexK < nEquation; indexK++)
                        {
                            dRdq[indexI][colidx1+indexJ] += dGradtmpBC[indexI][indexK] * dqdcvL[indexK][indexJ];
                            dRdq[indexI][colidx2+indexJ] -= dGradtmpBC[indexI][indexK] * dqdcvL[indexK][indexJ];
                        }
                    }
                }
            }

            //! judge whether le is in the left cell lists of the BC 
            result = find(BCLeftCells.begin(), BCLeftCells.end(), re); 
            if(result != BCLeftCells.end())
            {
                int dist     = distance(BCLeftCells.begin(),result);
                int rBCIndex = BCRightCells[dist];
                int FBCIndex = BCFaces[dist];
                //! contributions from the bc
                int indexre = (rBCIndex-nTotalCell)*nEquation;
                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {
                        dDdPlocal[indexI][indexJ] = dDdP[indexI][indexre+indexJ];
                    }
                }

                RDouble fnx    = nsx[FBCIndex];
                RDouble fny    = nsy[FBCIndex];
                RDouble fnz    = nsz[FBCIndex];
                RDouble farea  = ns[FBCIndex];
                RDouble volume = vol[re];

                for (int m = 0; m < nEquation; m++)
                {
                    primN[m] = prim[m][rBCIndex];
                }

                gas->dGradient2dPrimitive(primN,1,dgraddqx, 'x', fnx, fny, fnz, farea, volume, nEquation);
                gas->dGradient2dPrimitive(primN,1,dgraddqy, 'y', fnx, fny, fnz, farea, volume, nEquation);
                gas->dGradient2dPrimitive(primN,1,dgraddqz, 'z', fnx, fny, fnz, farea, volume, nEquation);

                //! contributions from the bc, 
                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {
                        dGradtmp[indexI][indexJ] = 0;

                        for(int indexK=0; indexK < nEquation; indexK++)
                        {
                            dGradtmp[indexI][indexJ] += dFluxdgradxR[indexI][indexK]*
                                                                    dgraddqx[indexK][indexJ];
                            dGradtmp[indexI][indexJ] += dFluxdgradyR[indexI][indexK]*
                                                                    dgraddqy[indexK][indexJ];
                            dGradtmp[indexI][indexJ] += dFluxdgradzR[indexI][indexK]*
                                                                    dgraddqz[indexK][indexJ];
                        }
                    }
                }

                //! contributions from the bc, 
                //! dFluxdPR*dDdPlocal
                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {
                        dGradtmpBC[indexI][indexJ] = 0.0;
                        for(int indexK=0; indexK < nEquation; indexK++)
                        {
                            dGradtmpBC[indexI][indexJ]  += dGradtmp[indexI][indexK]*
                                                            dDdPlocal[indexK][indexJ];
                        }
                    }
                }

                gas->dPrimitive2dConservative(primRC, gmR, dqdcvR);
                int colidx1 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, re);
                int colidx2 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, re);
                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {
                        for(int indexK=0; indexK < nEquation; indexK++)
                        {
                            dRdq[indexI][colidx1+indexJ] += dGradtmpBC[indexI][indexK] * dqdcvR[indexK][indexJ];
                            dRdq[indexI][colidx2+indexJ] -= dGradtmpBC[indexI][indexK] * dqdcvR[indexK][indexJ];
                        }
                    }
                }
            }
            #endif

            for (int nVar = 0; nVar < limitnVariable; ++ nVar)
            {
                for( int n = 0; n < nEquation; n++ )
                {
                    for( int m = 0; m < nEquation; m++ )
                    {
                        dFluxdMinMax[n][m] = 0.0;
                    }
                }

                #ifdef outputJac
                    std::cout << "------ nVar: " << nVar << " ------\n";
                    std::cout << "LimitSignL["<<nVar<<"] = " << LimitSignL[nVar] << "  \n";
                    std::cout << "LimitSignR["<<nVar<<"] = " << LimitSignR[nVar] << "  \n";
                    std::cout << "dFluxdMinMax: \n";
                    for( int n = 0; n < nEquation; n++ )
                    {
                        for( int m = 0; m < nEquation; m++ )
                        {
                            std::cout << dFluxdMinMax[n][m] << "  ";
                        }
                        std::cout << "\n";
                    }
                    std::cout << "\n"; 
                #endif

                if( LimitSignL[nVar] == 1 )
                {
                    int indexmax = IndexMaxL[nVar];
                    if( indexmax != le )
                    {
                        for (int m = 0; m < nEquation; m++)
                        {
                            dFluxdMinMax[m][nVar] = dFluxdmaxL[m][nVar];
                            primN[m]              = prim[m][indexmax];
                        }

                        gas->dPrimitive2dConservative(primN, gama[indexmax], dqdcvN);

                        int colidx1 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, indexmax);
                        int colidx2 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, indexmax);
                        for(int indexI=0; indexI < nEquation; indexI++)
                        {
                            for(int indexJ=0; indexJ < nEquation; indexJ++)
                            {
                                for(int indexK=0; indexK < nEquation; indexK++)
                                {
                                    dRdq[indexI][colidx1+indexJ] += dFluxdMinMax[indexI][indexK]*
                                                                            dqdcvN[indexK][indexJ];

                                    dRdq[indexI][colidx2+indexJ] -= dFluxdMinMax[indexI][indexK]*
                                                                            dqdcvN[indexK][indexJ];
                                }
                            }
                        }
                        //! GMRES2ndBCs
                        if( indexmax >= nTotalCell )
                        {
                            //! contributions from the bc
                            int indexre = (indexmax-nTotalCell)*nEquation;
                            for(int indexI=0; indexI < nEquation; indexI++)
                            {
                                for(int indexJ=0; indexJ < nEquation; indexJ++)
                                {
                                    dDdPlocal[indexI][indexJ] = dDdP[indexI][indexre+indexJ];
                                }
                            }

                            for(int indexI=0; indexI < nEquation; indexI++)
                            {
                                for(int indexJ=0; indexJ < nEquation; indexJ++)
                                {
                                    dGradtmpBC[indexI][indexJ] = 0.0;
                                    for(int indexK=0; indexK < nEquation; indexK++)
                                    {
                                        dGradtmpBC[indexI][indexJ]  += dFluxdMinMax[indexI][indexK]*
                                                                        dDdPlocal[indexK][indexJ];
                                    }
                                }
                            }

                            gas->dPrimitive2dConservative(primLC, gmL, dqdcvL);
                            int colidx1     = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, le);
                            int colidx2     = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, le);
                            for(int indexI=0; indexI < nEquation; indexI++)
                            {
                                for(int indexJ=0; indexJ < nEquation; indexJ++)
                                {
                                    for(int indexK=0; indexK < nEquation; indexK++)
                                    {
                                        dRdq[indexI][colidx1+indexJ] += dGradtmpBC[indexI][indexK] * dqdcvL[indexK][indexJ];
                                        dRdq[indexI][colidx2+indexJ] -= dGradtmpBC[indexI][indexK] * dqdcvL[indexK][indexJ];
                                    }
                                }
                            }
                        }
                    }

                    #ifdef outputJac
                        if( indexmax == le )
                        {       
                                std::cout << "Left max: \n";
                                for( int n = 0; n < nEquation; n++ )
                                {
                                    for( int m = 0; m < nEquation; m++ )
                                    {
                                        std::cout << dFluxdMinMax[n][m] << "  ";
                                    }
                                    std::cout << "\n";
                                }
                                std::cout << "\n";
                        }
                    #endif
                }
                else if ( LimitSignL[nVar] == -1 )
                {
                    int indexmin = IndexMinL[nVar];
                    if( indexmin != le )
                    {
                        for (int m = 0; m < nEquation; m++)
                        {
                            dFluxdMinMax[m][nVar] = dFluxdminL[m][nVar];
                            primN[m]              = prim[m][indexmin];
                        }

                        gas->dPrimitive2dConservative(primN, gama[indexmin], dqdcvN);

                        int colidx1 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, indexmin);
                        int colidx2 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, indexmin);
                        for(int indexI=0; indexI < nEquation; indexI++)
                        {
                            for(int indexJ=0; indexJ < nEquation; indexJ++)
                            {
                                for(int indexK=0; indexK < nEquation; indexK++)
                                {
                                    
                                    dRdq[indexI][colidx1+indexJ] += dFluxdMinMax[indexI][indexK]*
                                                                            dqdcvN[indexK][indexJ];
                                    
                                    dRdq[indexI][colidx2+indexJ] -= dFluxdMinMax[indexI][indexK]*
                                                                            dqdcvN[indexK][indexJ];
                                }
                            }
                        }

                        //! GMRES2ndBCs
                        if( indexmin >= nTotalCell )
                        {
                            //! contributions from the bc
                            int indexre = (indexmin-nTotalCell)*nEquation;
                            for(int indexI=0; indexI < nEquation; indexI++)
                            {
                                for(int indexJ=0; indexJ < nEquation; indexJ++)
                                {
                                    dDdPlocal[indexI][indexJ] = dDdP[indexI][indexre+indexJ];
                                }
                            }

                            //! contributions from the bc, 
                            //! dFluxdPR*dDdPlocal
                            for(int indexI=0; indexI < nEquation; indexI++)
                            {
                                for(int indexJ=0; indexJ < nEquation; indexJ++)
                                {
                                    dGradtmpBC[indexI][indexJ] = 0.0;
                                    for(int indexK=0; indexK < nEquation; indexK++)
                                    {
                                        dGradtmpBC[indexI][indexJ]  += dFluxdMinMax[indexI][indexK]*
                                                                        dDdPlocal[indexK][indexJ];
                                    }
                                }
                            }

                            gas->dPrimitive2dConservative(primLC, gmL, dqdcvL);
                            int colidx1 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, le);
                            int colidx2 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, le);
                            for(int indexI=0; indexI < nEquation; indexI++)
                            {
                                for(int indexJ=0; indexJ < nEquation; indexJ++)
                                {
                                    for(int indexK=0; indexK < nEquation; indexK++)
                                    {
                                        dRdq[indexI][colidx1+indexJ] += dGradtmpBC[indexI][indexK] * dqdcvL[indexK][indexJ];
                                        dRdq[indexI][colidx2+indexJ] -= dGradtmpBC[indexI][indexK] * dqdcvL[indexK][indexJ];
                                    }
                                }
                            }
                        }
                    }

                    #ifdef outputJac
                        if( indexmin == le )
                        {
                            std::cout << "Left min: \n";
                            for( int n = 0; n < nEquation; n++ )
                            {
                                for( int m = 0; m < nEquation; m++ )
                                {
                                    std::cout << dFluxdMinMax[n][m] << "  ";
                                }
                                std::cout << "\n";
                            }
                            std::cout << "\n";
                        }
                    #endif
                }

                if( LimitSignR[nVar] == 1 )
                {
                    int indexmax = IndexMaxR[nVar];
                    if( indexmax != re )
                    {
                        for (int m = 0; m < nEquation; m++)
                        {
                            dFluxdMinMax[m][nVar] = dFluxdmaxR[m][nVar];
                            primN[m]              = prim[m][indexmax];
                        }

                        gas->dPrimitive2dConservative(primN, gama[indexmax], dqdcvN);

                        int colidx1     = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, indexmax);
                        int colidx2     = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, indexmax);
                        for(int indexI=0; indexI < nEquation; indexI++)
                        {
                            for(int indexJ=0; indexJ < nEquation; indexJ++)
                            {
                                for(int indexK=0; indexK < nEquation; indexK++)
                                {
                                    dRdq[indexI][colidx1+indexJ] += dFluxdMinMax[indexI][indexK]*
                                                                            dqdcvN[indexK][indexJ];

                                    dRdq[indexI][colidx2+indexJ] -= dFluxdMinMax[indexI][indexK]*
                                                                            dqdcvN[indexK][indexJ];
                                }
                            }
                        }

                        //! GMRES2ndBCs
                        if( indexmax >= nTotalCell )
                        {
                            //! contributions from the bc
                            int indexre = (indexmax-nTotalCell)*nEquation;
                            for(int indexI=0; indexI < nEquation; indexI++)
                            {
                                for(int indexJ=0; indexJ < nEquation; indexJ++)
                                {
                                    dDdPlocal[indexI][indexJ] = dDdP[indexI][indexre+indexJ];
                                }
                            }

                            for(int indexI=0; indexI < nEquation; indexI++)
                            {
                                for(int indexJ=0; indexJ < nEquation; indexJ++)
                                {
                                    dGradtmpBC[indexI][indexJ] = 0.0;
                                    for(int indexK=0; indexK < nEquation; indexK++)
                                    {
                                        dGradtmpBC[indexI][indexJ]  += dFluxdMinMax[indexI][indexK]*
                                                                        dDdPlocal[indexK][indexJ];
                                    }
                                }
                            }

                            gas->dPrimitive2dConservative(primRC, gmR, dqdcvR);
                            int colidx1     = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, re);
                            int colidx2     = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, re);
                            for(int indexI=0; indexI < nEquation; indexI++)
                            {
                                for(int indexJ=0; indexJ < nEquation; indexJ++)
                                {
                                    for(int indexK=0; indexK < nEquation; indexK++)
                                    {
                                        dRdq[indexI][colidx1+indexJ] += dGradtmpBC[indexI][indexK] * dqdcvR[indexK][indexJ];
                                        dRdq[indexI][colidx2+indexJ] -= dGradtmpBC[indexI][indexK] * dqdcvR[indexK][indexJ];
                                    }
                                }
                            }
                        }
                    }

                    #ifdef outputJac
                        if( indexmax == re )
                        {       
                                std::cout << "Right max: \n";
                                for( int n = 0; n < nEquation; n++ )
                                {
                                    for( int m = 0; m < nEquation; m++ )
                                    {
                                        std::cout << dFluxdMinMax[n][m] << "  ";
                                    }
                                    std::cout << "\n";
                                }
                                std::cout << "\n";
                        }
                    #endif
                }
                else if ( LimitSignR[nVar] == -1 )
                {
                    int indexmin = IndexMinR[nVar];
                    if( indexmin != re )
                    {
                        for (int m = 0; m < nEquation; m++)
                        {
                            dFluxdMinMax[m][nVar] = dFluxdminR[m][nVar];
                            primN[m]              = prim[m][indexmin];
                        }

                        gas->dPrimitive2dConservative(primN, gama[indexmin], dqdcvN);

                        int colidx1 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, indexmin);
                        int colidx2 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, indexmin);
                        for(int indexI=0; indexI < nEquation; indexI++)
                        {
                            for(int indexJ=0; indexJ < nEquation; indexJ++)
                            {
                                for(int indexK=0; indexK < nEquation; indexK++)
                                {
                                    dRdq[indexI][colidx1+indexJ] += dFluxdMinMax[indexI][indexK]*
                                                                            dqdcvN[indexK][indexJ];

                                    dRdq[indexI][colidx2+indexJ] -= dFluxdMinMax[indexI][indexK]*
                                                                            dqdcvN[indexK][indexJ];
                                }
                            }
                        }

                        //! GMRES2ndBCs
                        if( indexmin >= nTotalCell )
                        {
                            //! contributions from the bc
                            int indexre = (indexmin-nTotalCell)*nEquation;
                            for(int indexI=0; indexI < nEquation; indexI++)
                            {
                                for(int indexJ=0; indexJ < nEquation; indexJ++)
                                {
                                    dDdPlocal[indexI][indexJ] = dDdP[indexI][indexre+indexJ];
                                }
                            }

                            for(int indexI=0; indexI < nEquation; indexI++)
                            {
                                for(int indexJ=0; indexJ < nEquation; indexJ++)
                                {
                                    dGradtmpBC[indexI][indexJ] = 0.0;
                                    for(int indexK=0; indexK < nEquation; indexK++)
                                    {
                                        dGradtmpBC[indexI][indexJ] += dFluxdMinMax[indexI][indexK]*
                                                                        dDdPlocal[indexK][indexJ];
                                    }
                                }
                            }

                            gas->dPrimitive2dConservative(primRC, gmR, dqdcvR);
                            int colidx1 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, re);
                            int colidx2 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, re);
                            for(int indexI=0; indexI < nEquation; indexI++)
                            {
                                for(int indexJ=0; indexJ < nEquation; indexJ++)
                                {
                                    for(int indexK=0; indexK < nEquation; indexK++)
                                    {
                                        dRdq[indexI][colidx1+indexJ] += dGradtmpBC[indexI][indexK] * dqdcvR[indexK][indexJ];
                                        dRdq[indexI][colidx2+indexJ] -= dGradtmpBC[indexI][indexK] * dqdcvR[indexK][indexJ];
                                    }
                                }
                            }
                        }
                    }

                     #ifdef outputJac
                        if( indexmin == re )
                        {
                            std::cout << "Right min: \n";
                            for( int n = 0; n < nEquation; n++ )
                            {
                                for( int m = 0; m < nEquation; m++ )
                                {
                                    std::cout << dFluxdMinMax[n][m] << "  ";
                                }
                                std::cout << "\n";
                            }
                            std::cout << "\n";
                        }
                    #endif
                }
            }
        }
    }

    delete [] primL;    primL = nullptr;
    delete [] primR;    primR = nullptr;
    delete [] primN;    primN = nullptr;
    delete [] primLC;   primLC = nullptr;
    delete [] primRC;   primRC = nullptr;
    delete [] fluxp;    fluxp = nullptr;

    for(int index = 0; index < nEquation; index++)
    {
        delete [] dqdcvL[index];    dqdcvL[index] = nullptr;
        delete [] dqdcvR[index];    dqdcvR[index] = nullptr;
        delete [] dqdcvN[index];    dqdcvN[index] = nullptr;
        delete [] dgraddqx[index];    dgraddqx[index] = nullptr;
        delete [] dgraddqy[index];    dgraddqy[index] = nullptr;
        delete [] dgraddqz[index];    dgraddqz[index] = nullptr;
    }
    delete [] dqdcvL;    dqdcvL = nullptr;
    delete [] dqdcvR;    dqdcvR = nullptr;
    delete [] dqdcvN;    dqdcvN = nullptr;
    delete [] dgraddqx;  dgraddqx = nullptr;
    delete [] dgraddqy;  dgraddqy = nullptr;
    delete [] dgraddqz;  dgraddqz = nullptr;

    delete [] dqdxL;    dqdxL = nullptr;
    delete [] dqdyL;    dqdyL = nullptr;
    delete [] dqdzL;    dqdzL = nullptr;

    delete [] dqdxR;    dqdxR = nullptr;
    delete [] dqdyR;    dqdyR = nullptr;
    delete [] dqdzR;    dqdzR = nullptr;
    if ( limiterType == ILMT_VENCAT )
    {
        delete [] limitL;    limitL = nullptr;
        delete [] limitR;    limitR = nullptr;
    }
    delete [] primTry;    primTry = nullptr;

    delete [] IndexMinL;    IndexMinL = nullptr;
    delete [] IndexMaxL;    IndexMaxL = nullptr;
    delete [] IndexMinR;    IndexMinR = nullptr;
    delete [] IndexMaxR;    IndexMaxR = nullptr;

    delete [] minL;    minL = nullptr;
    delete [] maxL;    maxL = nullptr;
    delete [] minR;    minR = nullptr;
    delete [] maxR;    maxR = nullptr;

    delete [] LimitSignL;    LimitSignL = nullptr;
    delete [] LimitSignR;    LimitSignR = nullptr;
}

//! GMRESnolim
void GMRES_Roe_Scheme_ConservativeForm_nolim_2nd(FaceProxy *face_proxy, InviscidSchemeParameter *invSchemePara)
{
    const RDouble LITTLE    = 1.0E-10;
    int nTemperatureModel   = invSchemePara->GetNumberOfTemperatureModel();
    int nEquation           = invSchemePara->GetNumberOfTotalEquation();
    int nLaminar            = invSchemePara->GetNumberOfLaminar();
    int nLength             = invSchemePara->GetLength();
    int nm                  = invSchemePara->Getnm();
    int RoeEntropyFixMethod = invSchemePara->GetEntropyFixMethod();
    int nChemical           = invSchemePara->GetNumberOfChemical();
    int nPrecondition       = invSchemePara->GetIfPrecondition();
    int nIdealState         = invSchemePara->GetIfIdealState();
    int limiterType         = invSchemePara->GetLimiterType();
    // limiterType             = ILMT_FIRST;

    //! GMRESBoundary
    int localStart = face_proxy->GetlocalStart();
    int localEnd   = face_proxy->GetlocalEnd();
    int nBoundFace = face_proxy->GetnBoundFace();
    int nTotalCell = face_proxy->GetnTotalCell();
    int nMid;
    //! GMRESBCorrection
    vector<int> *wallFaceIndex = face_proxy->GetWallFaceIndex();

    RDouble refMachNumber     = invSchemePara->GetMachNumber();
    RDouble entrFixCoefLeft   = invSchemePara->GetLeftEntropyFixCoefficient();
    RDouble entrFixCoefRight  = invSchemePara->GetRightEntropyFixCoefficient();
    RDouble **qL              = invSchemePara->GetLeftQ();
    RDouble **qR              = invSchemePara->GetRightQ();
    RDouble **qLC             = invSchemePara->GetLeftQC();    //! GMRESPassQC
    RDouble **qRC             = invSchemePara->GetRightQC();   //! GMRESPassQC
    RDouble *gamaL            = invSchemePara->GetLeftGama();
    RDouble *gamaR            = invSchemePara->GetRightGama();
    RDouble *pressureCoeffL   = face_proxy->GetPressureCoefficientL();
    RDouble *pressureCoeffR   = face_proxy->GetPressureCoefficientR();
    RDouble *xfn              = invSchemePara->GetFaceNormalX();
    RDouble *yfn              = invSchemePara->GetFaceNormalY();
    RDouble *zfn              = invSchemePara->GetFaceNormalZ();
    RDouble *vgn              = invSchemePara->GetFaceVelocity();
    RDouble **flux            = invSchemePara->GetFlux();
    RDouble **tl              = invSchemePara->GetLeftTemperature();
    RDouble **tr              = invSchemePara->GetRightTemperature();
    int *leftcellindexofFace  = face_proxy->GetLeftCellIndexOfFace();
    int *rightcellindexofFace = face_proxy->GetRightCellIndexOfFace();
    GeomProxy *geomProxy      = face_proxy->GetGeomProxy();
    RDouble* areas            = geomProxy->GetFaceArea();

    // GMRESnolim  xfc xcc should use the absolute face index
    RDouble *xfc   = invSchemePara->GetFaceCenterX();
    RDouble *yfc   = invSchemePara->GetFaceCenterY();
    RDouble *zfc   = invSchemePara->GetFaceCenterZ();
    RDouble *xcc   = invSchemePara->GetCellCenterX();
    RDouble *ycc   = invSchemePara->GetCellCenterY();
    RDouble *zcc   = invSchemePara->GetCellCenterZ();
    RDouble *nsx   = invSchemePara->GetFaceNormalXAbs();
    RDouble *nsy   = invSchemePara->GetFaceNormalYAbs();
    RDouble *nsz   = invSchemePara->GetFaceNormalZAbs();
    RDouble *ns    = invSchemePara->GetFaceAreaAbs();
    RDouble *vol   = invSchemePara->GetVolume();
    RDouble *gama  = invSchemePara->GetGama();
    RDouble **prim = invSchemePara->GetPrimitive();
    vector<int> *neighborCells = invSchemePara->GetNeighborCells();
    vector<int> *neighborFaces = invSchemePara->GetNeighborFaces();
    vector<int> *neighborLR    = invSchemePara->GetNeighborLR();
    int *qlsign                = invSchemePara->GetLeftQSign();
    int *qrsign                = invSchemePara->GetRightQSign();

    //! GMRES2ndCorrection
    vector<int> BCLeftCells  = invSchemePara->GetBCLeftCells();
    vector<int> BCRightCells = invSchemePara->GetBCRightCells();
    vector<int> BCFaces      = invSchemePara->GetBCFaces();

    RDouble dxL, dyL, dzL, dxR, dyR, dzR;
    int Lsign;
    int Rsign;

    RDouble** dRdq = face_proxy->GetJacobianMatrix();
    //! GMRESBoundary
    RDouble** dDdP = face_proxy->GetdDdPMatrix();
    //! GMRES CSR
    vector<int> AI = face_proxy->GetJacobianAI4GMRES();
    vector<int> AJ = face_proxy->GetJacobianAJ4GMRES();

    //! RDouble   dRdqL[nEquation][nEquation];
    RDouble dRdqR[nEquation][nEquation];
    RDouble dDdPlocal[nEquation][nEquation];
    RDouble dGdqx, dGdqy, dGdqz;
    RDouble mach2 = refMachNumber * refMachNumber;

    using namespace IDX;
    using namespace GAS_SPACE;

    RDouble dFluxdpL[nEquation][nEquation];
    RDouble dFluxdpR[nEquation][nEquation];
    RDouble dFluxdtmp[nEquation][nEquation];
    RDouble dFluxdpwall[nEquation][nEquation];

    //! auto difference allocation
    ADReal *primL = new ADReal[nEquation]();
    ADReal *primR = new ADReal[nEquation]();
    ADReal *fluxp = new ADReal[nEquation]();    //! flux proxy

    RDouble *primLC = new RDouble[nEquation]();    //! GMRESPassQC
    RDouble *primRC = new RDouble[nEquation]();    //! GMRESPassQC
    RDouble *primN  = new RDouble[nEquation]();

    RDouble temperaturesL[3] = { 0 }, temperaturesR[3] = { 0 };
    RDouble perturbScale = 0.001;

    RDouble **dqdcvL= new RDouble*[nEquation];
    RDouble **dqdcvR= new RDouble*[nEquation];
    RDouble **dqdcvN= new RDouble*[nEquation];    //! dqdcv for neighbors
    RDouble **dgraddqx = new RDouble*[nEquation];
    RDouble **dgraddqy = new RDouble*[nEquation];
    RDouble **dgraddqz = new RDouble*[nEquation];
    for(int indexI=0; indexI < nEquation; indexI++)
    {
        dqdcvL[indexI]   = new RDouble[nEquation];
        dqdcvR[indexI]   = new RDouble[nEquation];
        dqdcvN[indexI]   = new RDouble[nEquation];
        dgraddqx[indexI] = new RDouble[nEquation];
        dgraddqy[indexI] = new RDouble[nEquation];
        dgraddqz[indexI] = new RDouble[nEquation];

        for(int indexJ =0; indexJ < nEquation; indexJ ++)
        {
            dqdcvL[indexI][indexJ] = 0.0;
            dqdcvR[indexI][indexJ] = 0.0;
            dqdcvN[indexI][indexJ] = 0.0;
            dgraddqx[indexI][indexJ] = 0.0;
            dgraddqy[indexI][indexJ] = 0.0;
            dgraddqz[indexI][indexJ] = 0.0;
        }
    }

    //! GMRESBoundary
     if (localStart >= nBoundFace)
     {
         nMid = localStart;    //! a bug 01.07
     }
     else if (localEnd <= nBoundFace)
     {
         //! If they are all boundary faces.
         nMid = localEnd;
     }
     else
     {
         //! Part of them are boundary faces.
         nMid = nBoundFace;
     }

     nMid = nMid - localStart;

     for(int iFace = 0;iFace < nMid; ++iFace)
     {
        int le = leftcellindexofFace[iFace];
        int re = rightcellindexofFace[iFace];
        int colidx;    //! colume index for GMRES CSR

        RDouble area    = areas[iFace];

        //! GMRESgrad
        for (int iEqn = 0; iEqn < nEquation; ++iEqn)
        {
            primL[iEqn] = qL[iEqn][iFace];
            primR[iEqn] = qR[iEqn][iFace];

            primLC[iEqn] = qLC[iEqn][iFace];    //! GMRESPassQC
            primRC[iEqn] = qRC[iEqn][iFace];    //! GMRESPassQC
        }

        //! define sequence order of the independent variables
        for (int iEqn = 0; iEqn < nEquation; iEqn++)
        {
            primL[iEqn].diff(iEqn, 10);
            primR[iEqn].diff(iEqn + nEquation, 10);
        }

        RDouble gmL = gamaL[iFace];
        RDouble gmR = gamaR[iFace];

        Cal_GMRES_Roe_Scheme_Flux_PV(primL, primR, iFace, fluxp,
                                  face_proxy,invSchemePara);

        flux[IR][iFace] = fluxp[IR].val();
        flux[IU][iFace] = fluxp[IU].val();
        flux[IV][iFace] = fluxp[IV].val();
        flux[IW][iFace] = fluxp[IW].val();
        flux[IP][iFace] = fluxp[IP].val();

        for (int m = 0; m < nEquation; m++)
        {
            for (int n = 0; n < nEquation; n++)
            {
                dFluxdpL[m][n] = fluxp[m].dx(n) * area;
                dFluxdpR[m][n] = fluxp[m].dx(n + nEquation) * area;
            }
        }

        //! contributions from the bc
        int indexre = (re-nTotalCell)*nEquation;
        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                dDdPlocal[indexI][indexJ] = dDdP[indexI][indexre+indexJ];
            }
        }

        //! GMRESBCorrection
        if(find(wallFaceIndex->begin(),wallFaceIndex->end(),iFace) != wallFaceIndex->end())
        {
            //! GMRESPV , obtain the dqdcv for the left cell
            gas->dPrimitive2dConservative(primLC,gmL,dqdcvL);

            //! average of the dFluxdpL and dFluxdpR
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ =0; indexJ < nEquation; indexJ ++)
                {
                    dFluxdpwall[indexI][indexJ] = 0.5*(dFluxdpL[indexI][indexJ]
                                                     + dFluxdpR[indexI][indexJ]);
                    dFluxdpL[indexI][indexJ]    =  dFluxdpwall[indexI][indexJ];
                    dFluxdpR[indexI][indexJ]    =  dFluxdpwall[indexI][indexJ];
                }
            }

            //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
            //! then minus to the corresponding location in dRdq.
            int indexf = re * nEquation;    //! first index
            int indexs = le * nEquation;    //! second index
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, le);
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        dRdq[indexI][colidx + indexJ] -= dFluxdpwall[indexI][indexK] * dqdcvL[indexK][indexJ];
                    }
                }
            }

            //! contributions from the bc, 
            //! dFluxdpL = dFluxdpL + dFluxdPR*dDdPlocal
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        dFluxdpL[indexI][indexJ]  += dFluxdpwall[indexI][indexK]*
                                                       dDdPlocal[indexK][indexJ];
                    }
                }
            }

            //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
            //! then add to the corresponding location in dRdq.
            indexf = le * nEquation;    //! first index
            indexs = le * nEquation;    //! second index
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, le);
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        dRdq[indexI][colidx + indexJ] += dFluxdpL[indexI][indexK] * dqdcvL[indexK][indexJ];
                    }
                }
            }

            //! GMRESPV , obtain the dqdcv for the right cell
            gas->dPrimitive2dConservative(primRC,gmR,dqdcvR);

            //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
            //! then minus to the corresponding location in dRdq.
            indexf = re * nEquation;    //! first index
            indexs = re * nEquation;    //! second index
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, re);
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        dRdq[indexI][colidx + indexJ] -= dFluxdpR[indexI][indexK] * dqdcvR[indexK][indexJ];
                    }
                }
            }

            //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
            //! then add to the corresponding location in dRdq.
            indexf = le * nEquation;    //! first index
            indexs = re * nEquation;    //! second index
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, re);
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        dRdq[indexI][colidx + indexJ] += dFluxdpR[indexI][indexK] * dqdcvR[indexK][indexJ];
                    }
                }
            }
            
        }
        else
        {
            if ( limiterType == ILMT_FIRST )
            {
                //! GMRESPV , obtain the dqdcv for the left cell
                gas->dPrimitive2dConservative(primLC,gmL,dqdcvL);

                //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
                //! then minus to the corresponding location in dRdq.
                int indexf = re * nEquation;    //! first index
                int indexs = le * nEquation;    //! second index
                colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, le);
                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {
                        for(int indexK=0; indexK < nEquation; indexK++)
                        {
                            dRdq[indexI][colidx + indexJ] -= dFluxdpL[indexI][indexK] * dqdcvL[indexK][indexJ];
                        }
                    }
                }

                //! contributions from the bc, 
                //! dFluxdpL = dFluxdpL + dFluxdPR*dDdPlocal
                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {
                        for(int indexK=0; indexK < nEquation; indexK++)
                        {
                            dFluxdpL[indexI][indexJ]  += dFluxdpR[indexI][indexK]*
                                                        dDdPlocal[indexK][indexJ];
                        }
                    }
                }

                //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
                //! then add to the corresponding location in dRdq.
                indexf = le * nEquation;    //! first index
                indexs = le * nEquation;    //! second index
                colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, le);
                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {
                        for(int indexK=0; indexK < nEquation; indexK++)
                        {
                            dRdq[indexI][colidx + indexJ] += dFluxdpL[indexI][indexK] * dqdcvL[indexK][indexJ];
                        }
                    }
                }

                 //! GMRESPV , obtain the dqdcv for the right cell
                gas->dPrimitive2dConservative(primRC,gmR,dqdcvR);

                //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
                //! then minus to the corresponding location in dRdq.
                indexf = re * nEquation;    //! first index
                indexs = re * nEquation;    //! second index
                colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, re);
                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {
                        for(int indexK=0; indexK < nEquation; indexK++)
                        {
                            dRdq[indexI][colidx + indexJ] -= dFluxdpR[indexI][indexK] * dqdcvR[indexK][indexJ];
                        }
                    }
                }

                //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
                //! then add to the corresponding location in dRdq.
                indexf = le * nEquation;    //! first index
                indexs = re * nEquation;    //! second index
                colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, re);
                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {
                        for(int indexK=0; indexK < nEquation; indexK++)
                        {
                            dRdq[indexI][colidx + indexJ] += dFluxdpR[indexI][indexK] * dqdcvR[indexK][indexJ];
                        }
                    }
                }
            }
            else if ( limiterType == ILMT_NOLIM )
            {
                Lsign = qlsign[iFace];
                Rsign = qrsign[iFace];
                int jFace =  iFace + localStart;    //! absolute face index

                //! GMRESgrad
                if( Lsign > 0)
                {
                    dxL = xfc[jFace] - xcc[le];
                    dyL = yfc[jFace] - ycc[le];
                    dzL = zfc[jFace] - zcc[le];
                }
                else
                {
                    dxL = 0.0;
                    dyL = 0.0;
                    dzL = 0.0;
                }

                 dxR = 0.0;
                 dyR = 0.0;
                 dzR = 0.0;

                //! GMRESPV , obtain the dqdcv for the left cell
                gas->dPrimitive2dConservative(primLC,gmL,dqdcvL);

                //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
                //! then minus to the corresponding location in dRdq.
                int indexf = re * nEquation;    //! first index
                int indexs = le * nEquation;    //! second index
                colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, le);
                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    { 
                        for(int indexK=0; indexK < nEquation; indexK++)
                        {
                            dRdq[indexI][colidx + indexJ] -= dFluxdpL[indexI][indexK] * dqdcvL[indexK][indexJ];
                        }
                    }
                }

                //! contributions from the bc, 
                //! dFluxdpL = dFluxdpL + dFluxdPR*dDdPlocal
                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {
                        dFluxdtmp[indexI][indexJ] = 0.0;
                        for(int indexK=0; indexK < nEquation; indexK++)
                        {
                            dFluxdtmp[indexI][indexJ]  +=   dFluxdpR[indexI][indexK]*
                                                            dDdPlocal[indexK][indexJ];
                        }
                    }
                }

                //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
                //! then add to the corresponding location in dRdq.
                indexf = le * nEquation;    //! first index
                indexs = le * nEquation;    //! second index
                colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, le);
                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {
                        for(int indexK=0; indexK < nEquation; indexK++)
                        {
                            dRdq[indexI][colidx + indexJ] += dFluxdtmp[indexI][indexK] * dqdcvL[indexK][indexJ];
                            dRdq[indexI][colidx + indexJ] += dFluxdpL[indexI][indexK] * dqdcvL[indexK][indexJ];
                        }
                    }
                }

                //! GMRESPV , obtain the dqdcv for the right cell
                gas->dPrimitive2dConservative(primRC,gmR,dqdcvR);

                //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
                //! then minus to the corresponding location in dRdq.
                indexf = re * nEquation;    //! first index
                indexs = re * nEquation;    //! second index
                colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, re);
                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {
                        for(int indexK=0; indexK < nEquation; indexK++)
                        {
                            dRdq[indexI][colidx + indexJ] -= dFluxdpR[indexI][indexK] * dqdcvR[indexK][indexJ];
                        }
                    }
                }

                //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
                //! then add to the corresponding location in dRdq.
                indexf = le * nEquation;    //! first index
                indexs = re * nEquation;    //! second index
                colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, re);
                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {
                        for(int indexK=0; indexK < nEquation; indexK++)
                        {
                            dRdq[indexI][colidx + indexJ] += dFluxdpR[indexI][indexK] * dqdcvR[indexK][indexJ];
                        }
                    }
                }

                //! considering the perturbation of the gradient in the limiter
                RDouble fnx    = nsx[jFace];
                RDouble fny    = nsy[jFace];
                RDouble fnz    = nsz[jFace];
                RDouble farea  = ns[jFace];
                RDouble volume = vol[re];

                //! dgradRdqL
                dGdqx = -1.0*0.5/volume*farea*fnx;
                dGdqy = -1.0*0.5/volume*farea*fny;
                dGdqz = -1.0*0.5/volume*farea*fnz;

                //! dQre/dgradR*dgradRdqL
                RDouble dGrad = dGdqx * dxR + dGdqy * dyR + dGdqz * dzR;

                //! convert dGrad by right multiplying the dqdcvL,
                //! then add to the corresponding location in dRdq.
                indexf = le * nEquation;    //! first index
                indexs = le * nEquation;    //! second index
                colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, le);
                int colidx2 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, le);
                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {
                        for(int indexK=0; indexK < nEquation; indexK++)
                        {
                            dRdq[indexI][colidx+indexJ] += dGrad* dFluxdpR[indexI][indexK] * dqdcvL[indexK][indexJ];
                            dRdq[indexI][colidx2+indexJ] -= dGrad* dFluxdpR[indexI][indexK] * dqdcvL[indexK][indexJ];
                        }
                    }
                }

                //! considering the perturbation of the gradient in the limiter
                volume = vol[le];

                //! dgradLdqR
                dGdqx = 1.0*0.5/volume*farea*fnx;
                dGdqy = 1.0*0.5/volume*farea*fny;
                dGdqz = 1.0*0.5/volume*farea*fnz;

                //! dQle/dgradL*dgradLdqR
                dGrad = dGdqx * dxL + dGdqy * dyL + dGdqz * dzL;

                //! convert dGrad by right multiplying the dqdcvL,
                //! then add to the corresponding location in dRdq.

                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {
                        dFluxdtmp[indexI][indexJ] = 0.0;
                        for(int indexK=0; indexK < nEquation; indexK++)
                        {
                            dFluxdtmp[indexI][indexJ]  += dGrad* dFluxdpL[indexI][indexK]*
                                                            dDdPlocal[indexK][indexJ];
                        }
                    }
                }

                colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, re);
                colidx2 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, re);
                int colidx3 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, le);
                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {
                        for(int indexK=0; indexK < nEquation; indexK++)
                        {
                            dRdq[indexI][colidx+indexJ] += dGrad* dFluxdpL[indexI][indexK] * dqdcvR[indexK][indexJ];
                            dRdq[indexI][colidx2+indexJ] -= dGrad* dFluxdpL[indexI][indexK] * dqdcvR[indexK][indexJ];
                            dRdq[indexI][colidx3+indexJ] += dFluxdtmp[indexI][indexK] * dqdcvL[indexK][indexJ];
                        }
                    }
                }

                //! the neighbors of the left cells
                for(int index = 0; index < neighborCells[le].size(); index++)
                {
                    int neighborIndex = neighborCells[le][index];
                    
                    if(neighborIndex != re)
                    {
                        int Faceindex  = neighborFaces[le][index];
                        int sign       = neighborLR[le][index];
                        RDouble fnx    = nsx[Faceindex];
                        RDouble fny    = nsy[Faceindex];
                        RDouble fnz    = nsz[Faceindex];
                        RDouble farea  = ns[Faceindex];
                        RDouble volume = vol[le];

                        for(int m =0; m < nEquation; m++)
                        {
                            primN[m] = prim[m][neighborIndex];
                        }
                        gas->dPrimitive2dConservative(primN,gama[neighborIndex],dqdcvN);
                        //! dgradLdqN
                        dGdqx = sign*0.5/volume*farea*fnx;
                        dGdqy = sign*0.5/volume*farea*fny;
                        dGdqz = sign*0.5/volume*farea*fnz;

                        //! dQle/dgradL*dgradLdqN
                        dGrad = dGdqx * dxL + dGdqy * dyL + dGdqz * dzL;

                        int indexf  = le * nEquation;    //! first index
                        int indexs  = neighborIndex * nEquation;    //! second index
                        int indexf2 = re * nEquation; 
                        int colidx  = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, neighborIndex);
                        int colidx2 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, neighborIndex);
                        
                        for(int indexI=0; indexI < nEquation; indexI++)
                        {
                            for(int indexJ=0; indexJ < nEquation; indexJ++)
                            {
                                for(int indexK=0; indexK < nEquation; indexK++)
                                {
                                    //! dRL/dqLNeighbors 
                                    dRdq[indexI][colidx+indexJ] += dGrad* dFluxdpL[indexI][indexK]*
                                                                             dqdcvN[indexK][indexJ];

                                    //! dRR/dqLNeighbors 
                                    dRdq[indexI][colidx2+indexJ] -= dGrad* dFluxdpL[indexI][indexK]*
                                                                             dqdcvN[indexK][indexJ];
                                }
                            }
                        }
                    }
                }

                //! the neighbors of the right cells
                for(int index = 0; index < neighborCells[re].size(); index++)
                {
                    int neighborIndex = neighborCells[re][index];
                    // std::cout << "Neighbors of the right cells: " << neighborIndex << " \n";
                    if(neighborIndex != le)
                    {
                        int Faceindex  = neighborFaces[re][index];
                        int sign       = neighborLR[re][index];
                        RDouble fnx    = nsx[Faceindex];
                        RDouble fny    = nsy[Faceindex];
                        RDouble fnz    = nsz[Faceindex];
                        RDouble farea  = ns[Faceindex];
                        RDouble volume = vol[re];

                        for(int m =0; m < nEquation; m++)
                        {
                            primN[m] = prim[m][neighborIndex];
                        }

                        gas->dPrimitive2dConservative(primN,gama[neighborIndex],dqdcvN);

                        //! dgradRdqN
                        dGdqx = sign*0.5/volume*farea*fnx;
                        dGdqy = sign*0.5/volume*farea*fny;
                        dGdqz = sign*0.5/volume*farea*fnz;

                        //! dQre/dgradR*dgradRdqN
                        dGrad = dGdqx * dxR + dGdqy * dyR + dGdqz * dzR;

                        int indexf  = le * nEquation; // first index
                        int indexs  = neighborIndex * nEquation; // second index
                        int indexf2 = re * nEquation; 
                        int colidx  = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, neighborIndex);
                        int colidx2 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, neighborIndex);
                        for(int indexI=0; indexI < nEquation; indexI++)
                        {
                            for(int indexJ=0; indexJ < nEquation; indexJ++)
                            {
                                for(int indexK=0; indexK < nEquation; indexK++)
                                {
                                    //! dRL/dqRNeighbors 
                                    dRdq[indexI][colidx+indexJ] += dGrad* dFluxdpR[indexI][indexK]*
                                                                             dqdcvN[indexK][indexJ];

                                    //! dRR/dqRNeighbors
                                    dRdq[indexI][colidx2+indexJ] -= dGrad* dFluxdpR[indexI][indexK]*
                                                                             dqdcvN[indexK][indexJ];
                                }
                            }
                        }
                }
            }
        }
     }

     for(int iFace = nMid;iFace < nLength; ++iFace)
     {
        int le = leftcellindexofFace[iFace];
        int re = rightcellindexofFace[iFace];
        int colidx;
        RDouble area  = areas[iFace];
        for (int iEqn = 0; iEqn < nEquation; ++iEqn)
        {
            primL[iEqn] = qL[iEqn][iFace];
            primR[iEqn] = qR[iEqn][iFace];

            primLC[iEqn] = qLC[iEqn][iFace];    //! GMRESPassQC
            primRC[iEqn] = qRC[iEqn][iFace];    //! GMRESPassQC
        }

        //! define sequence order of the independent variables
        for (int iEqn = 0; iEqn < nEquation; iEqn++)
        {
            primL[iEqn].diff(iEqn, 10);
            primR[iEqn].diff(iEqn + nEquation, 10);
        }

        RDouble gmL = gamaL[iFace];
        RDouble gmR = gamaR[iFace];

        Cal_GMRES_Roe_Scheme_Flux_PV(primL, primR, iFace, fluxp,
                                  face_proxy,invSchemePara);

        flux[IR][iFace] = fluxp[IR].val();
        flux[IU][iFace] = fluxp[IU].val();
        flux[IV][iFace] = fluxp[IV].val();
        flux[IW][iFace] = fluxp[IW].val();
        flux[IP][iFace] = fluxp[IP].val();

        for (int m = 0; m < nEquation; m++)
        {
             for (int n = 0; n < nEquation; n++)
             {
                 dFluxdpL[m][n] = fluxp[m].dx(n) * area;
                 dFluxdpR[m][n] = fluxp[m].dx(n + nEquation) * area;
             }
        }

        if (limiterType == ILMT_FIRST)
        {
            //! GMRESPV , obtain the dqdcv for the left cell
            gas->dPrimitive2dConservative(primLC,gmL,dqdcvL);

            //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
            //! then minus to the corresponding location in dRdq.
            int indexf = re * nEquation;    //! first index
            int indexs = le * nEquation;    //! second index
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, le);
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                { 
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        dRdq[indexI][colidx + indexJ] -= dFluxdpL[indexI][indexK] * dqdcvL[indexK][indexJ];
                    }
                }
            }

            //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
            //! then add to the corresponding location in dRdq.
            indexf = le * nEquation;    //! first index
            indexs = le * nEquation;    //! second index
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, le);
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        dRdq[indexI][colidx + indexJ] += dFluxdpL[indexI][indexK] * dqdcvL[indexK][indexJ];
                    }
                }
            }

            //! GMRESPV , obtain the dqdcv for the right cell
            gas->dPrimitive2dConservative(primRC,gmR,dqdcvR);

            //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
            //! then minus to the corresponding location in dRdq.
            indexf = re * nEquation;    //! first index
            indexs = re * nEquation;    //! second index
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, re);
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        dRdq[indexI][colidx + indexJ] -= dFluxdpR[indexI][indexK] * dqdcvR[indexK][indexJ];
                    }
                }
            }

            //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
            //! then add to the corresponding location in dRdq.
            indexf = le * nEquation;    //! first index
            indexs = re * nEquation;    //! second index
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, re);
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        dRdq[indexI][colidx + indexJ] += dFluxdpR[indexI][indexK] * dqdcvR[indexK][indexJ];
                    }
                }
            }
        }
        else if ( limiterType == ILMT_NOLIM )
        {
            Lsign = qlsign[iFace];
            Rsign = qrsign[iFace];
            int jFace =  iFace + localStart;    //! absolute face index

            if( Lsign > 0 )
            {
                dxL = xfc[jFace] - xcc[le];
                dyL = yfc[jFace] - ycc[le];
                dzL = zfc[jFace] - zcc[le];
            }
            else
            {
                dxL = 0.0;
                dyL = 0.0;
                dzL = 0.0;
            }

            if( Rsign > 0 )
            {
                dxR = xfc[jFace] - xcc[re];
                dyR = yfc[jFace] - ycc[re];
                dzR = zfc[jFace] - zcc[re];
            }
            else
            {
                dxR = 0.0;
                dyR = 0.0;
                dzR = 0.0;
            }

            //! GMRESPV , obtain the dqdcv for the left cell
            gas->dPrimitive2dConservative(primLC,gmL,dqdcvL);

            //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
            //! then minus to the corresponding location in dRdq.
            int indexf = re * nEquation;    //! first index
            int indexs = le * nEquation;    //! second index
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, le);
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                { 
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        dRdq[indexI][colidx + indexJ] -= dFluxdpL[indexI][indexK] * dqdcvL[indexK][indexJ];
                    }
                }
            }

            //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
            //! then add to the corresponding location in dRdq.
            indexf = le * nEquation;    //! first index
            indexs = le * nEquation;    //! second index
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, le);
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        dRdq[indexI][colidx + indexJ] += dFluxdpL[indexI][indexK] * dqdcvL[indexK][indexJ];
                    }
                }
            }

            //! GMRESPV , obtain the dqdcv for the right cell
            gas->dPrimitive2dConservative(primRC,gmR,dqdcvR);

            //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
            //! then minus to the corresponding location in dRdq.
            indexf = re * nEquation;    //! first index
            indexs = re * nEquation;    //! second index
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, re);
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        dRdq[indexI][colidx + indexJ] -= dFluxdpR[indexI][indexK] * dqdcvR[indexK][indexJ];
                    }
                }
            }

            //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
            //! then add to the corresponding location in dRdq.
            indexf = le * nEquation;    //! first index
            indexs = re * nEquation;    //! second index
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, re);
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        dRdq[indexI][colidx + indexJ] += dFluxdpR[indexI][indexK] * dqdcvR[indexK][indexJ];
                    }
                }
            }

            //! considering the perturbation of the gradient in the limiter
            RDouble fnx    = nsx[jFace];
            RDouble fny    = nsy[jFace];
            RDouble fnz    = nsz[jFace];
            RDouble farea  = ns[jFace];
            RDouble volume = vol[re];

            //! dgradRdqL
            dGdqx = -1.0*0.5/volume*farea*fnx;
            dGdqy = -1.0*0.5/volume*farea*fny;
            dGdqz = -1.0*0.5/volume*farea*fnz;

            //! dQre/dgradR*dgradRdqL
            RDouble dGrad = dGdqx * dxR + dGdqy * dyR + dGdqz * dzR;

            //! convert dGrad by right multiplying the dqdcvL,
            //! then add to the corresponding location in dRdq.
            indexf = le * nEquation;    //! first index
            indexs = le * nEquation;    //! second index
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, le);
            int colidx2 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, le);
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        dRdq[indexI][colidx+indexJ] += dGrad* dFluxdpR[indexI][indexK] * dqdcvL[indexK][indexJ];
                        dRdq[indexI][colidx2+indexJ] -= dGrad* dFluxdpR[indexI][indexK] * dqdcvL[indexK][indexJ];
                    }
                }
            }

            //! considering the perturbation of the gradient in the limiter
            volume      = vol[le];

            //! dgradLdqR
            dGdqx = 1.0*0.5/volume*farea*fnx;
            dGdqy = 1.0*0.5/volume*farea*fny;
            dGdqz = 1.0*0.5/volume*farea*fnz;

            //! dQle/dgradL*dgradLdqR
            dGrad = dGdqx * dxL + dGdqy * dyL + dGdqz * dzL;

            //! convert dGrad by right multiplying the dqdcvL,
            //! then add to the corresponding location in dRdq.

            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, re);
            colidx2 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, re);
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        dRdq[indexI][colidx+indexJ] += dGrad* dFluxdpL[indexI][indexK] * dqdcvR[indexK][indexJ];
                        dRdq[indexI][colidx2+indexJ] -= dGrad* dFluxdpL[indexI][indexK] * dqdcvR[indexK][indexJ];
                    }
                }
            }

            // std::cout << "iFace: " << iFace << " le, re " << le << " , " << re << "\n";
            // the neighbors of the left cells
            for(int index = 0; index < neighborCells[le].size(); index++)
            {
                int neighborIndex = neighborCells[le][index];
                // std::cout << "Neighbors of the left cells: " << neighborIndex << " \n";
                if(neighborIndex != re)
                {
                    int Faceindex  = neighborFaces[le][index];
                    int sign       = neighborLR[le][index];
                    RDouble fnx    = nsx[Faceindex];
                    RDouble fny    = nsy[Faceindex];
                    RDouble fnz    = nsz[Faceindex];
                    RDouble farea  = ns[Faceindex];
                    RDouble volume = vol[le];

                    for(int m =0; m < nEquation; m++)
                    {
                        primN[m] = prim[m][neighborIndex];
                    }

                    gas->dPrimitive2dConservative(primN,gama[neighborIndex],dqdcvN);

                    //! dgradLdqN
                    dGdqx = sign*0.5/volume*farea*fnx;
                    dGdqy = sign*0.5/volume*farea*fny;
                    dGdqz = sign*0.5/volume*farea*fnz;

                    //! dQle/dgradL*dgradLdqN
                    dGrad = dGdqx * dxL + dGdqy * dyL + dGdqz * dzL;

                    int indexf  = le * nEquation;    //! first index
                    int indexs  = neighborIndex * nEquation;    //! second index
                    int indexf2 = re * nEquation; 
                    int colidx  = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, neighborIndex);
                    int colidx2 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, neighborIndex);
                    // std::cout << "dRle/dqNeighbor: " << le << " " << neighborIndex <<  " " << colidx << " \n";
                    // std::cout << "dRre/dqNeighbor: " << re << " " << neighborIndex <<  " " << colidx2 << " \n";
                    
                    for(int indexI=0; indexI < nEquation; indexI++)
                    {
                        for(int indexJ=0; indexJ < nEquation; indexJ++)
                        {

                            for(int indexK=0; indexK < nEquation; indexK++)
                            {
                                //! dRL/dqLNeighbors 
                                dRdq[indexI][colidx+indexJ] += dGrad* dFluxdpL[indexI][indexK]*
                                                                         dqdcvN[indexK][indexJ];

                                //! dRR/dqLNeighbors 
                                dRdq[indexI][colidx2+indexJ] -= dGrad* dFluxdpL[indexI][indexK]*
                                                                         dqdcvN[indexK][indexJ];
                            }
                        }
                    }
                }
            }

            //! the neighbors of the right cells
            for(int index = 0; index < neighborCells[re].size(); index++)
            {
                int neighborIndex = neighborCells[re][index];
                // std::cout << "Neighbors of the right cells: " << neighborIndex << " \n";

                if(neighborIndex != le)
                {
                    int Faceindex  = neighborFaces[re][index];
                    int sign       = neighborLR[re][index];
                    RDouble fnx    = nsx[Faceindex];
                    RDouble fny    = nsy[Faceindex];
                    RDouble fnz    = nsz[Faceindex];
                    RDouble farea  = ns[Faceindex];
                    RDouble volume = vol[re];

                    for(int m =0; m < nEquation; m++)
                    {
                        primN[m] = prim[m][neighborIndex];
                    }

                    gas->dPrimitive2dConservative(primN,gama[neighborIndex],dqdcvN);

                    //! dgradRdqN
                    dGdqx = sign*0.5/volume*farea*fnx;
                    dGdqy = sign*0.5/volume*farea*fny;
                    dGdqz = sign*0.5/volume*farea*fnz;

                    //! dQre/dgradR*dgradRdqN
                    dGrad = dGdqx * dxR + dGdqy * dyR + dGdqz * dzR;

                    int indexf  = le * nEquation; // first index
                    int indexs  = neighborIndex * nEquation; // second index
                    int indexf2 = re * nEquation; 
                    int colidx  = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, neighborIndex);
                    int colidx2 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, neighborIndex);
                    // std::cout << "dRle/dqNeighbor: " << le << " " << neighborIndex <<  " " << colidx << " \n";
                    // std::cout << "dRre/dqNeighbor: " << re << " " << neighborIndex <<  " " << colidx2 << " \n";
                    for(int indexI=0; indexI < nEquation; indexI++)
                    {
                        for(int indexJ=0; indexJ < nEquation; indexJ++)
                        {

                            for(int indexK=0; indexK < nEquation; indexK++)
                            {
                                //! dRL/dqRNeighbors 
                                dRdq[indexI][colidx+indexJ] += dGrad* dFluxdpR[indexI][indexK]*
                                                                         dqdcvN[indexK][indexJ];

                                //! dRR/dqRNeighbors
                                dRdq[indexI][colidx2+indexJ] -= dGrad* dFluxdpR[indexI][indexK]*
                                                                         dqdcvN[indexK][indexJ];
                            }
                        }
                    }
                }
            }

            //! judge whether le is in the left cell lists of the BC 
            vector<int>::iterator result = find(BCLeftCells.begin(), BCLeftCells.end(), le); 
            if(result != BCLeftCells.end())
            {
                int dist     = distance(BCLeftCells.begin(),result);
                int rBCIndex = BCRightCells[dist];
                int FBCIndex = BCFaces[dist];

                //! contributions from the bc
                int indexre = (rBCIndex-nTotalCell)*nEquation;
                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {
                        dDdPlocal[indexI][indexJ] = dDdP[indexI][indexre+indexJ];
                    }
                }

                RDouble fnx    = nsx[FBCIndex];
                RDouble fny    = nsy[FBCIndex];
                RDouble fnz    = nsz[FBCIndex];
                RDouble farea  = ns[FBCIndex];
                RDouble volume = vol[le];

                //! dgradLdqN
                dGdqx = 1.0*0.5/volume*farea*fnx;
                dGdqy = 1.0*0.5/volume*farea*fny;
                dGdqz = 1.0*0.5/volume*farea*fnz;
                //! dQle/dgradL*dgradLdqN
                dGrad = dGdqx * dxL + dGdqy * dyL + dGdqz * dzL;

                //! contributions from the bc, 
                //! dFluxdpL = dFluxdpL + dFluxdPR*dDdPlocal

                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {
                        dFluxdtmp[indexI][indexJ] = 0.0;
                        for(int indexK=0; indexK < nEquation; indexK++)
                        {
                            dFluxdtmp[indexI][indexJ]  += dGrad* dFluxdpL[indexI][indexK]*
                                                            dDdPlocal[indexK][indexJ];
                        }
                    }
                }

                gas->dPrimitive2dConservative(primLC,gmL,dqdcvL);
                int colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, le);
                int colidx2 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, le);
                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {
                        for(int indexK=0; indexK < nEquation; indexK++)
                        {
                            dRdq[indexI][colidx+indexJ] += dFluxdtmp[indexI][indexK] * dqdcvL[indexK][indexJ];
                            dRdq[indexI][colidx2+indexJ] -= dFluxdtmp[indexI][indexK] * dqdcvL[indexK][indexJ];
                        }
                    }
                }
            }

            //! judge whether re is in the left cell lists of the BC 
            result = find(BCLeftCells.begin(), BCLeftCells.end(), re); 
            if(result != BCLeftCells.end())
            {
                int dist     = distance(BCLeftCells.begin(),result);
                int rBCIndex = BCRightCells[dist];
                int FBCIndex = BCFaces[dist];

                //! contributions from the bc
                int indexre = (rBCIndex-nTotalCell)*nEquation;
                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {
                        dDdPlocal[indexI][indexJ] = dDdP[indexI][indexre+indexJ];
                    }
                }

                RDouble fnx    = nsx[FBCIndex];
                RDouble fny    = nsy[FBCIndex];
                RDouble fnz    = nsz[FBCIndex];
                RDouble farea  = ns[FBCIndex];
                RDouble volume = vol[re];

                //! dgradLdqN
                dGdqx = 1.0*0.5/volume*farea*fnx;
                dGdqy = 1.0*0.5/volume*farea*fny;
                dGdqz = 1.0*0.5/volume*farea*fnz;
                //! dQle/dgradL*dgradLdqN
                dGrad = dGdqx * dxR + dGdqy * dyR + dGdqz * dzR;

                //! contributions from the bc, 
                //! dFluxdpL = dFluxdpL + dFluxdPR*dDdPlocal

                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {
                        dFluxdtmp[indexI][indexJ] = 0.0;
                        for(int indexK=0; indexK < nEquation; indexK++)
                        {
                            dFluxdtmp[indexI][indexJ]  += dGrad* dFluxdpR[indexI][indexK]*
                                                            dDdPlocal[indexK][indexJ];
                        }
                    }
                }

                gas->dPrimitive2dConservative(primRC,gmR,dqdcvR);
                int colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, re);
                int colidx2 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, re);
                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {
                        for(int indexK=0; indexK < nEquation; indexK++)
                        {
                            dRdq[indexI][colidx+indexJ] += dFluxdtmp[indexI][indexK] * dqdcvR[indexK][indexJ];
                            dRdq[indexI][colidx2+indexJ] -= dFluxdtmp[indexI][indexK] * dqdcvR[indexK][indexJ];
                        }
                    }
                }
            }
        }
     }

    delete [] primL;    primL = nullptr;
    delete [] primR;    primR = nullptr;
    delete [] primN;    primN = nullptr;
    delete [] primLC;   primLC = nullptr;
    delete [] primRC;   primRC = nullptr;
    delete [] fluxp;    fluxp = nullptr;

    for(int index = 0; index < nEquation; index++)
    {
        delete [] dqdcvL[index];    dqdcvL[index] = nullptr;
        delete [] dqdcvR[index];    dqdcvR[index] = nullptr;
        delete [] dqdcvN[index];    dqdcvN[index] = nullptr;
        delete [] dgraddqx[index];    dgraddqx[index] = nullptr;
        delete [] dgraddqy[index];    dgraddqy[index] = nullptr;
        delete [] dgraddqz[index];    dgraddqz[index] = nullptr;
    }
    delete [] dqdcvL;    dqdcvL = nullptr;
    delete [] dqdcvR;    dqdcvR = nullptr;
    delete [] dqdcvN;    dqdcvN = nullptr;
    delete [] dgraddqx;  dgraddqx = nullptr;
    delete [] dgraddqy;  dgraddqy = nullptr;
    delete [] dgraddqz;  dgraddqz = nullptr;
}
}

//! GMRESnolim
void GMRES_Roe_Scheme_ConservativeForm_nolim_Matrix(FaceProxy *face_proxy, InviscidSchemeParameter *invSchemePara)
{
    const RDouble LITTLE    = 1.0E-10;
    int nTemperatureModel   = invSchemePara->GetNumberOfTemperatureModel();
    int nEquation           = invSchemePara->GetNumberOfTotalEquation();
    int nLaminar            = invSchemePara->GetNumberOfLaminar();
    int nLength             = invSchemePara->GetLength();
    int nm                  = invSchemePara->Getnm();
    int RoeEntropyFixMethod = invSchemePara->GetEntropyFixMethod();
    int nChemical           = invSchemePara->GetNumberOfChemical();
    int nPrecondition       = invSchemePara->GetIfPrecondition();
    int nIdealState         = invSchemePara->GetIfIdealState();
    int limiterType         = invSchemePara->GetLimiterType();

    //! GMRESBoundary
    int localStart  = face_proxy->GetlocalStart();
    int localEnd    = face_proxy->GetlocalEnd();
    int nBoundFace  = face_proxy->GetnBoundFace();
    int nTotalCell  = face_proxy->GetnTotalCell();
    int nMid;
    //! GMRESBCorrection
    vector<int>* wallFaceIndex = face_proxy->GetWallFaceIndex();

    RDouble refMachNumber     = invSchemePara->GetMachNumber();
    RDouble entrFixCoefLeft   = invSchemePara->GetLeftEntropyFixCoefficient();
    RDouble entrFixCoefRight  = invSchemePara->GetRightEntropyFixCoefficient();
    RDouble **qL              = invSchemePara->GetLeftQ();
    RDouble **qR              = invSchemePara->GetRightQ();
    RDouble **qLC             = invSchemePara->GetLeftQC();    //! GMRESPassQC
    RDouble **qRC             = invSchemePara->GetRightQC();   //! GMRESPassQC
    RDouble *gamaL            = invSchemePara->GetLeftGama();
    RDouble *gamaR            = invSchemePara->GetRightGama();
    RDouble *pressureCoeffL   = face_proxy->GetPressureCoefficientL();
    RDouble *pressureCoeffR   = face_proxy->GetPressureCoefficientR();
    RDouble *xfn              = invSchemePara->GetFaceNormalX();
    RDouble *yfn              = invSchemePara->GetFaceNormalY();
    RDouble *zfn              = invSchemePara->GetFaceNormalZ();
    RDouble *vgn              = invSchemePara->GetFaceVelocity();
    RDouble **flux            = invSchemePara->GetFlux();
    RDouble **tl              = invSchemePara->GetLeftTemperature();
    RDouble **tr              = invSchemePara->GetRightTemperature();
    int *leftcellindexofFace  = face_proxy->GetLeftCellIndexOfFace();
    int *rightcellindexofFace = face_proxy->GetRightCellIndexOfFace();
    GeomProxy *geomProxy      = face_proxy->GetGeomProxy();
    RDouble *areas            = geomProxy->GetFaceArea();

    //! GMRESnolim  xfc xcc should use the absolute face index
    RDouble *xfc               = invSchemePara->GetFaceCenterX();
    RDouble *yfc               = invSchemePara->GetFaceCenterY();
    RDouble *zfc               = invSchemePara->GetFaceCenterZ();
    RDouble *xcc               = invSchemePara->GetCellCenterX();
    RDouble *ycc               = invSchemePara->GetCellCenterY();
    RDouble *zcc               = invSchemePara->GetCellCenterZ();
    RDouble *nsx               = invSchemePara->GetFaceNormalXAbs();
    RDouble *nsy               = invSchemePara->GetFaceNormalYAbs();
    RDouble *nsz               = invSchemePara->GetFaceNormalZAbs();
    RDouble *ns                = invSchemePara->GetFaceAreaAbs();
    RDouble *vol               = invSchemePara->GetVolume();
    RDouble *gama              = invSchemePara->GetGama();
    RDouble **prim             = invSchemePara->GetPrimitive();
    vector<int> *neighborCells = invSchemePara->GetNeighborCells();
    vector<int> *neighborFaces = invSchemePara->GetNeighborFaces();
    vector<int> *neighborLR    = invSchemePara->GetNeighborLR();
    int *qlsign                = invSchemePara->GetLeftQSign();
    int *qrsign                = invSchemePara->GetRightQSign();
    RDouble dxL, dyL, dzL, dxR, dyR, dzR;
    int Lsign;
    int Rsign;

    RDouble **dRdq = face_proxy->GetJacobianMatrix();
    //! GMRESBoundary
    RDouble **dDdP = face_proxy->GetdDdPMatrix();
    //! GMRES CSR
    vector<int> AI = face_proxy->GetJacobianAI4GMRES();
    vector<int> AJ = face_proxy->GetJacobianAJ4GMRES();

    RDouble dRdqR[nEquation][nEquation];
    RDouble dDdPlocal[nEquation][nEquation];
    RDouble dGdqx, dGdqy, dGdqz;
    RDouble mach2 = refMachNumber * refMachNumber;

    using namespace IDX;
    using namespace GAS_SPACE;

    RDouble dFluxdpL[nEquation][nEquation];
    RDouble dFluxdpR[nEquation][nEquation];
    RDouble dFluxdtmp[nEquation][nEquation];
    RDouble dFluxdpwall[nEquation][nEquation];

    //! auto difference allocation
    ADReal *primL  = new ADReal[nEquation]();
    ADReal *primR  = new ADReal[nEquation]();
    ADReal *fluxp  = new ADReal[nEquation]();    //! flux proxy

    RDouble *primLC = new RDouble[nEquation]();    //! GMRESPassQC
    RDouble *primRC = new RDouble[nEquation]();    //! GMRESPassQC
    RDouble *primN  = new RDouble[nEquation]();

    RDouble temperaturesL[3] = { 0 }, temperaturesR[3] = { 0 };
    RDouble perturbScale = 0.001;

    RDouble **dqdcvL= new RDouble*[nEquation];
    RDouble **dqdcvR= new RDouble*[nEquation];
    RDouble **dqdcvN= new RDouble*[nEquation];    //! dqdcv for neighbors
    RDouble **dgraddqx = new RDouble*[nEquation];
    RDouble **dgraddqy = new RDouble*[nEquation];
    RDouble **dgraddqz = new RDouble*[nEquation];
    for(int indexI=0; indexI < nEquation; indexI++)
    {
        dqdcvL[indexI]   = new RDouble[nEquation];
        dqdcvR[indexI]   = new RDouble[nEquation];
        dqdcvN[indexI]   = new RDouble[nEquation];
        dgraddqx[indexI] = new RDouble[nEquation];
        dgraddqy[indexI] = new RDouble[nEquation];
        dgraddqz[indexI] = new RDouble[nEquation];

        for(int indexJ =0; indexJ < nEquation; indexJ ++)
        {
            dqdcvL[indexI][indexJ] = 0.0;
            dqdcvR[indexI][indexJ] = 0.0;
            dqdcvN[indexI][indexJ] = 0.0;
            dgraddqx[indexI][indexJ] = 0.0;
            dgraddqy[indexI][indexJ] = 0.0;
            dgraddqz[indexI][indexJ] = 0.0;
        }
    }

    //! GMRESBoundary
     if (localStart >= nBoundFace)
     {
         nMid = localStart;    //1 a bug 01.07
     }
     else if (localEnd <= nBoundFace)
     {
         //! If they are all boundary faces.
         nMid = localEnd;
     }
     else
     {
         //! Part of them are boundary faces.
         nMid = nBoundFace;
     }

     nMid = nMid - localStart;

     for(int iFace = 0;iFace < nMid; ++iFace)
     {
        int le = leftcellindexofFace[iFace];
        int re = rightcellindexofFace[iFace];
        int colidx;    //! colume index for GMRES CSR

        RDouble area    = areas[iFace];
        for (int iEqn = 0; iEqn < nEquation; ++iEqn)
        {
            primL[iEqn] = qL[iEqn][iFace];
            primR[iEqn] = qR[iEqn][iFace];

            primLC[iEqn] = qLC[iEqn][iFace];    //! GMRESPassQC
            primRC[iEqn] = qRC[iEqn][iFace];    //! GMRESPassQC
        }

        //! define sequence order of the independent variables
        for (int iEqn = 0; iEqn < nEquation; iEqn++)
        {
            primL[iEqn].diff(iEqn, 10);
            primR[iEqn].diff(iEqn + nEquation, 10);
        }

        RDouble gmL = gamaL[iFace];
        RDouble gmR = gamaR[iFace];

        Cal_GMRES_Roe_Scheme_Flux_PV(primL, primR, iFace, fluxp,
                                  face_proxy,invSchemePara);

        flux[IR][iFace] = fluxp[IR].val();
        flux[IU][iFace] = fluxp[IU].val();
        flux[IV][iFace] = fluxp[IV].val();
        flux[IW][iFace] = fluxp[IW].val();
        flux[IP][iFace] = fluxp[IP].val();

        for (int m = 0; m < nEquation; m++)
        {
            for (int n = 0; n < nEquation; n++)
            {
                dFluxdpL[m][n] = fluxp[m].dx(n) * area;
                dFluxdpR[m][n] = fluxp[m].dx(n + nEquation) * area;
            }
        }
 
        //! contributions from the bc
        int indexre = (re-nTotalCell)*nEquation;
        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                dDdPlocal[indexI][indexJ] = dDdP[indexI][indexre+indexJ];
            }
        }

        //! GMRESBCorrection
        if(find(wallFaceIndex->begin(),wallFaceIndex->end(),iFace) != wallFaceIndex->end())
        {
            //! GMRESPV , obtain the dqdcv for the left cell
            gas->dPrimitive2dConservative(primLC,gmL,dqdcvL);

            //! average of the dFluxdpL and dFluxdpR
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ =0; indexJ < nEquation; indexJ ++)
                {
                    dFluxdpwall[indexI][indexJ] = 0.5*(dFluxdpL[indexI][indexJ]
                                                     + dFluxdpR[indexI][indexJ]);
                    dFluxdpL[indexI][indexJ]    =  dFluxdpwall[indexI][indexJ]; 
                    dFluxdpR[indexI][indexJ]    =  dFluxdpwall[indexI][indexJ];                               
                }
            }

            //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
            //! then minus to the corresponding location in dRdq.
            int indexf = re * nEquation;    //! first index
            int indexs = le * nEquation;    //! second index
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, le);
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        dRdq[indexI][colidx + indexJ] -= dFluxdpwall[indexI][indexK] * dqdcvL[indexK][indexJ];
                    }
                }
            }

            //! contributions from the bc, 
            //! dFluxdpL = dFluxdpL + dFluxdPR*dDdPlocal
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        dFluxdpL[indexI][indexJ]  += dFluxdpwall[indexI][indexK]*
                                                       dDdPlocal[indexK][indexJ];
                    }
                }
            }

            //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
            //! then add to the corresponding location in dRdq.
            indexf = le * nEquation;    //! first index
            indexs = le * nEquation;    //! second index
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, le);
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        dRdq[indexI][colidx + indexJ] += dFluxdpL[indexI][indexK] * dqdcvL[indexK][indexJ];
                    }
                }
            }

            //! GMRESPV , obtain the dqdcv for the right cell
            gas->dPrimitive2dConservative(primRC,gmR,dqdcvR);

            //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
            //! then minus to the corresponding location in dRdq.
            indexf = re * nEquation;    //! first index
            indexs = re * nEquation;    //! second index
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, re);
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        dRdq[indexI][colidx + indexJ] -= dFluxdpR[indexI][indexK] * dqdcvR[indexK][indexJ];
                    }
                }
            }

            //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
            //! then add to the corresponding location in dRdq.
            indexf = le * nEquation;    //! first index
            indexs = re * nEquation;    //! second index
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, re);
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        dRdq[indexI][colidx + indexJ] += dFluxdpR[indexI][indexK] * dqdcvR[indexK][indexJ];
                    }
                }
            }
        }
        else
        {
            if ( limiterType == ILMT_FIRST )
            {
                //! GMRESPV , obtain the dqdcv for the left cell
                gas->dPrimitive2dConservative(primLC,gmL,dqdcvL);

                //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
                //! then minus to the corresponding location in dRdq.
                int indexf = re * nEquation;    //! first index
                int indexs = le * nEquation;    //! second index
                colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, le);
                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {
                        for(int indexK=0; indexK < nEquation; indexK++)
                        {
                            dRdq[indexI][colidx + indexJ] -= dFluxdpL[indexI][indexK] * dqdcvL[indexK][indexJ];
                        }
                    }
                }

                //! contributions from the bc, 
                //! dFluxdpL = dFluxdpL + dFluxdPR*dDdPlocal
                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {
                        for(int indexK=0; indexK < nEquation; indexK++)
                        {
                            dFluxdpL[indexI][indexJ]  += dFluxdpR[indexI][indexK]*
                                                        dDdPlocal[indexK][indexJ];
                        }
                    }
                }

                //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
                //! then add to the corresponding location in dRdq.
                indexf = le * nEquation;    //! first index
                indexs = le * nEquation;    //! second index
                colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, le);
                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {
                        for(int indexK=0; indexK < nEquation; indexK++)
                        {
                            dRdq[indexI][colidx + indexJ] += dFluxdpL[indexI][indexK] * dqdcvL[indexK][indexJ];
                        }
                    }
                }

                //! GMRESPV , obtain the dqdcv for the right cell
                gas->dPrimitive2dConservative(primRC,gmR,dqdcvR);

                //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
                //! then minus to the corresponding location in dRdq.
                indexf = re * nEquation;    //! first index
                indexs = re * nEquation;    //! second index
                colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, re);
                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {
                        for(int indexK=0; indexK < nEquation; indexK++)
                        {
                            dRdq[indexI][colidx + indexJ] -= dFluxdpR[indexI][indexK] * dqdcvR[indexK][indexJ];
                        }
                    }
                }

                //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
                //! then add to the corresponding location in dRdq.
                indexf = le * nEquation;    //! first index
                indexs = re * nEquation;    //! second index
                colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, re);
                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {
                        for(int indexK=0; indexK < nEquation; indexK++)
                        {
                            dRdq[indexI][colidx + indexJ] += dFluxdpR[indexI][indexK] * dqdcvR[indexK][indexJ];
                        }
                    }
                }
            }
            else if ( limiterType == ILMT_NOLIM )
            {
                Lsign = qlsign[iFace];
                Rsign = qrsign[iFace];
                int jFace =  iFace + localStart;    //! absolute face index

                if( Lsign > 0)
                {
                    dxL = xfc[jFace] - xcc[le];
                    dyL = yfc[jFace] - ycc[le];
                    dzL = zfc[jFace] - zcc[le];
                }
                else
                {
                    dxL = 0.0;
                    dyL = 0.0;
                    dzL = 0.0;
                }

                if( Rsign > 0)
                {
                    dxR = xfc[jFace] - xcc[re];
                    dyR = yfc[jFace] - ycc[re];
                    dzR = zfc[jFace] - zcc[re];
                }
                else
                {
                    dxR = 0.0;
                    dyR = 0.0;
                    dzR = 0.0;
                }
                //! GMRESPV , obtain the dqdcv for the left cell
                gas->dPrimitive2dConservative(primLC,gmL,dqdcvL);

                //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
                //! then minus to the corresponding location in dRdq.
                int indexf = re * nEquation;    //! first index
                int indexs = le * nEquation;    //! second index
                colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, le);
                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    { 
                        for(int indexK=0; indexK < nEquation; indexK++)
                        {
                            dRdq[indexI][colidx + indexJ] -= dFluxdpL[indexI][indexK] * dqdcvL[indexK][indexJ];
                        }
                    }
                }

                //! contributions from the bc, 
                //! dFluxdpL = dFluxdpL + dFluxdPR*dDdPlocal
                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {
                        dFluxdtmp[indexI][indexJ] = 0.0;
                        for(int indexK=0; indexK < nEquation; indexK++)
                        {
                            dFluxdtmp[indexI][indexJ]  +=   dFluxdpR[indexI][indexK]*
                                                            dDdPlocal[indexK][indexJ];
                        }
                    }
                }

                //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
                //! then add to the corresponding location in dRdq.
                indexf = le * nEquation;    //! first index
                indexs = le * nEquation;    //! second index
                colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, le);
                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {
                        for(int indexK=0; indexK < nEquation; indexK++)
                        {
                            dRdq[indexI][colidx + indexJ] += dFluxdtmp[indexI][indexK] * dqdcvL[indexK][indexJ];
                        }
                    }
                }

                //! GMRESPV , obtain the dqdcv for the right cell
                gas->dPrimitive2dConservative(primRC,gmR,dqdcvR);

                //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
                //! then minus to the corresponding location in dRdq.
                indexf = re * nEquation;    //! first index
                indexs = re * nEquation;    //! second index
                colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, re);
                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {
                        for(int indexK=0; indexK < nEquation; indexK++)
                        {
                            dRdq[indexI][colidx + indexJ] -= dFluxdpR[indexI][indexK] * dqdcvR[indexK][indexJ];
                        }
                    }
                }

                //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
                //! then add to the corresponding location in dRdq.
                indexf = le * nEquation;    //! first index
                indexs = re * nEquation;    //! second index
                colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, re);
                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {
                        for(int indexK=0; indexK < nEquation; indexK++)
                        {
                            dRdq[indexI][colidx + indexJ] += dFluxdpR[indexI][indexK] * dqdcvR[indexK][indexJ];
                        }
                    }
                }

                //! considering the perturbation of the gradient in the limiter
                RDouble fnx    = nsx[jFace];
                RDouble fny    = nsy[jFace];
                RDouble fnz    = nsz[jFace];
                RDouble farea  = ns[jFace];
                RDouble volume = vol[re];

                //! dgradRdqL
                dGdqx = -1.0*0.5/volume*farea*fnx;
                dGdqy = -1.0*0.5/volume*farea*fny;
                dGdqz = -1.0*0.5/volume*farea*fnz;

                //! dQre/dgradR*dgradRdqL
                RDouble dGrad = dGdqx * dxR + dGdqy * dyR + dGdqz * dzR;

                //! convert dGrad by right multiplying the dqdcvL,
                //! then add to the corresponding location in dRdq.
                indexf = le * nEquation;    //! first index
                indexs = le * nEquation;    //! second index
                colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, le);
                int colidx2 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, le);
                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {
                        for(int indexK=0; indexK < nEquation; indexK++)
                        {
                            dRdq[indexI][colidx+indexJ] += dGrad* dFluxdpR[indexI][indexK] * dqdcvL[indexK][indexJ];
                            dRdq[indexI][colidx2+indexJ] -= dGrad* dFluxdpR[indexI][indexK] * dqdcvL[indexK][indexJ];
                        }
                    }
                }

                //! considering the perturbation of the gradient in the limiter
                volume = vol[le];

                //! dgradLdqR
                dGdqx = 1.0*0.5/volume*farea*fnx;
                dGdqy = 1.0*0.5/volume*farea*fny;
                dGdqz = 1.0*0.5/volume*farea*fnz;

                //! dQle/dgradL*dgradLdqR
                dGrad = dGdqx * dxL + dGdqy * dyL + dGdqz * dzL;

                //! convert dGrad by right multiplying the dqdcvL,
                //! then add to the corresponding location in dRdq.

                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {
                        dFluxdtmp[indexI][indexJ] = 0.0;
                        for(int indexK=0; indexK < nEquation; indexK++)
                        {
                            dFluxdtmp[indexI][indexJ]  += dGrad* dFluxdpL[indexI][indexK]*
                                                            dDdPlocal[indexK][indexJ];
                        }
                    }
                }

                colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, re);
                colidx2 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, re);
                int colidx3 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, le);
                for(int indexI=0; indexI < nEquation; indexI++)
                {
                    for(int indexJ=0; indexJ < nEquation; indexJ++)
                    {
                        for(int indexK=0; indexK < nEquation; indexK++)
                        {
                            dRdq[indexI][colidx+indexJ] += dGrad* dFluxdpL[indexI][indexK] * dqdcvR[indexK][indexJ];
                            dRdq[indexI][colidx2+indexJ] -= dGrad* dFluxdpL[indexI][indexK] * dqdcvR[indexK][indexJ];
                            dRdq[indexI][colidx3+indexJ] += dFluxdtmp[indexI][indexK] * dqdcvL[indexK][indexJ];
                        }
                    }
                }

                //! the neighbors of the left cells
                for(int index = 0; index < neighborCells[le].size(); index++)
                {
                    int neighborIndex = neighborCells[le][index];

                    if(neighborIndex != re)
                    {
                        int Faceindex  = neighborFaces[le][index];
                        int sign       = neighborLR[le][index];
                        RDouble fnx    = nsx[Faceindex];
                        RDouble fny    = nsy[Faceindex];
                        RDouble fnz    = nsz[Faceindex];
                        RDouble farea  = ns[Faceindex];
                        RDouble volume = vol[le];

                        for(int m =0; m < nEquation; m++)
                        {
                            primN[m] = prim[m][neighborIndex];
                        }
                        gas->dPrimitive2dConservative(primN,gama[neighborIndex],dqdcvN);
                        //! dgradLdqN
                        dGdqx = sign*0.5/volume*farea*fnx;
                        dGdqy = sign*0.5/volume*farea*fny;
                        dGdqz = sign*0.5/volume*farea*fnz;

                        //! dQle/dgradL*dgradLdqN
                        dGrad = dGdqx * dxL + dGdqy * dyL + dGdqz * dzL;

                        int indexf  = le * nEquation;    //! first index
                        int indexs  = neighborIndex * nEquation;    //! second index
                        int indexf2 = re * nEquation; 
                        int colidx  = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, neighborIndex);
                        int colidx2 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, neighborIndex);
                        for(int indexI=0; indexI < nEquation; indexI++)
                        {
                            for(int indexJ=0; indexJ < nEquation; indexJ++)
                            {
                                for(int indexK=0; indexK < nEquation; indexK++)
                                {
                                    //! dRL/dqLNeighbors 
                                    dRdq[indexI][colidx+indexJ] += dGrad* dFluxdpL[indexI][indexK]*
                                                                             dqdcvN[indexK][indexJ];

                                    //! dRR/dqLNeighbors 
                                    dRdq[indexI][colidx2+indexJ] -= dGrad* dFluxdpL[indexI][indexK]*
                                                                             dqdcvN[indexK][indexJ];
                                }
                            }
                        }
                    }
                }

                //! the neighbors of the right cells
                for(int index = 0; index < neighborCells[re].size(); index++)
                {
                    int neighborIndex = neighborCells[re][index];

                    if(neighborIndex != le)
                    {
                        int Faceindex  = neighborFaces[re][index];
                        int sign       = neighborLR[re][index];
                        RDouble fnx    = nsx[Faceindex];
                        RDouble fny    = nsy[Faceindex];
                        RDouble fnz    = nsz[Faceindex];
                        RDouble farea  = ns[Faceindex];
                        RDouble volume = vol[re];

                        for(int m =0; m < nEquation; m++)
                        {
                            primN[m] = prim[m][neighborIndex];
                        }

                        gas->dPrimitive2dConservative(primN,gama[neighborIndex],dqdcvN);

                        //! dgradRdqN
                        dGdqx = sign*0.5/volume*farea*fnx;
                        dGdqy = sign*0.5/volume*farea*fny;
                        dGdqz = sign*0.5/volume*farea*fnz;

                        //! dQre/dgradR*dgradRdqN
                        dGrad = dGdqx * dxR + dGdqy * dyR + dGdqz * dzR;

                        int indexf  = le * nEquation;    //! first index
                        int indexs  = neighborIndex * nEquation;    //! second index
                        int indexf2 = re * nEquation; 
                        int colidx  = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, neighborIndex);
                        int colidx2 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, neighborIndex);
                        for(int indexI=0; indexI < nEquation; indexI++)
                        {
                            for(int indexJ=0; indexJ < nEquation; indexJ++)
                            {

                                for(int indexK=0; indexK < nEquation; indexK++)
                                {
                                    //! dRL/dqRNeighbors 
                                    dRdq[indexI][colidx+indexJ] += dGrad* dFluxdpR[indexI][indexK]*
                                                                             dqdcvN[indexK][indexJ];

                                    //! dRR/dqRNeighbors
                                    dRdq[indexI][colidx2+indexJ] -= dGrad* dFluxdpR[indexI][indexK]*
                                                                             dqdcvN[indexK][indexJ];
                                }
                            }
                        }
                    }
                }
            }
        }
     }

     for(int iFace = nMid;iFace < nLength; ++iFace)
     {
        int le = leftcellindexofFace[iFace];
        int re = rightcellindexofFace[iFace];
        int colidx;
        RDouble area  = areas[iFace];
        for (int iEqn = 0; iEqn < nEquation; ++iEqn)
        {
            primL[iEqn] = qL[iEqn][iFace];
            primR[iEqn] = qR[iEqn][iFace];

            primLC[iEqn] = qLC[iEqn][iFace];    //! GMRESPassQC
            primRC[iEqn] = qRC[iEqn][iFace];    //! GMRESPassQC
        }

        //! define sequence order of the independent variables
        for (int iEqn = 0; iEqn < nEquation; iEqn++)
        {
            primL[iEqn].diff(iEqn, 10);
            primR[iEqn].diff(iEqn + nEquation, 10);
        }

        RDouble gmL = gamaL[iFace];
        RDouble gmR = gamaR[iFace];

        Cal_GMRES_Roe_Scheme_Flux_PV(primL, primR, iFace, fluxp,
                                  face_proxy,invSchemePara);

        flux[IR][iFace] = fluxp[IR].val();
        flux[IU][iFace] = fluxp[IU].val();
        flux[IV][iFace] = fluxp[IV].val();
        flux[IW][iFace] = fluxp[IW].val();
        flux[IP][iFace] = fluxp[IP].val();
       
        for (int m = 0; m < nEquation; m++)
        {
             for (int n = 0; n < nEquation; n++)
             {
                 dFluxdpL[m][n] = fluxp[m].dx(n) * area;
                 dFluxdpR[m][n] = fluxp[m].dx(n + nEquation) * area;
             }
        }
  
        if (limiterType == ILMT_FIRST)
        {
            //! GMRESPV , obtain the dqdcv for the left cell
            gas->dPrimitive2dConservative(primLC,gmL,dqdcvL);

            //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
            //! then minus to the corresponding location in dRdq.
            int indexf = re * nEquation;    //! first index
            int indexs = le * nEquation;    //! second index
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, le);
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                { 
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        dRdq[indexI][colidx + indexJ] -= dFluxdpL[indexI][indexK] * dqdcvL[indexK][indexJ];
                    }
                }
            }

            //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
            //! then add to the corresponding location in dRdq.
            indexf = le * nEquation;    //! first index
            indexs = le * nEquation;    //! second index
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, le);
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        dRdq[indexI][colidx + indexJ] += dFluxdpL[indexI][indexK] * dqdcvL[indexK][indexJ];
                    }
                }
            }

            //! GMRESPV , obtain the dqdcv for the right cell
            gas->dPrimitive2dConservative(primRC,gmR,dqdcvR);

            //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
            //! then minus to the corresponding location in dRdq.
            indexf = re * nEquation;    //! first index
            indexs = re * nEquation;    //! second index
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, re);
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        dRdq[indexI][colidx + indexJ] -= dFluxdpR[indexI][indexK] * dqdcvR[indexK][indexJ];
                    }
                }
            }

            //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
            //! then add to the corresponding location in dRdq.
            indexf = le * nEquation;    //! first index
            indexs = re * nEquation;    //! second index
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, re);
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        dRdq[indexI][colidx + indexJ] += dFluxdpR[indexI][indexK] * dqdcvR[indexK][indexJ];
                    }
                }
            }
        }
        else if ( limiterType == ILMT_NOLIM )
        {
            Lsign = qlsign[iFace];
            Rsign = qrsign[iFace];
            int jFace =  iFace + localStart;    //! absolute face index

            if( Lsign > 0 )
            {
                dxL = xfc[jFace] - xcc[le];
                dyL = yfc[jFace] - ycc[le];
                dzL = zfc[jFace] - zcc[le];
            }
            else
            {
                dxL = 0.0;
                dyL = 0.0;
                dzL = 0.0;
            }

            if( Rsign > 0 )
            {
                dxR = xfc[jFace] - xcc[re];
                dyR = yfc[jFace] - ycc[re];
                dzR = zfc[jFace] - zcc[re];
            }
            else
            {
                dxR = 0.0;
                dyR = 0.0;
                dzR = 0.0;
            }

            //! GMRESPV , obtain the dqdcv for the left cell
            gas->dPrimitive2dConservative(primLC,gmL,dqdcvL);

            //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
            //! then minus to the corresponding location in dRdq.
            int indexf = re * nEquation;    //! first index
            int indexs = le * nEquation;    //! second index
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, le);
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                { 
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        dRdq[indexI][colidx + indexJ] -= dFluxdpL[indexI][indexK] * dqdcvL[indexK][indexJ];
                    }
                }
            }

            //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
            //! then add to the corresponding location in dRdq.
            indexf = le * nEquation;    //! first index
            indexs = le * nEquation;    //! second index
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, le);
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        dRdq[indexI][colidx + indexJ] += dFluxdpL[indexI][indexK] * dqdcvL[indexK][indexJ];
                    }
                }
            }

            //! GMRESPV , obtain the dqdcv for the right cell
            gas->dPrimitive2dConservative(primRC,gmR,dqdcvR);

            //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
            //! then minus to the corresponding location in dRdq.
            indexf = re * nEquation; // first index
            indexs = re * nEquation; // second index
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, re);
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        dRdq[indexI][colidx + indexJ] -= dFluxdpR[indexI][indexK] * dqdcvR[indexK][indexJ];
                    }
                }
            }

            //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
            //! then add to the corresponding location in dRdq.
            indexf = le * nEquation;    //! first index
            indexs = re * nEquation;    //! second index
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, re);
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        dRdq[indexI][colidx + indexJ] += dFluxdpR[indexI][indexK] * dqdcvR[indexK][indexJ];
                    }
                }
            }

            //! considering the perturbation of the gradient in the limiter
            RDouble fnx    = nsx[jFace];
            RDouble fny    = nsy[jFace];
            RDouble fnz    = nsz[jFace];
            RDouble farea  = ns[jFace];
            RDouble volume = vol[re];

            //! dgradRdqL
            dGdqx = -1.0*0.5/volume*farea*fnx;
            dGdqy = -1.0*0.5/volume*farea*fny;
            dGdqz = -1.0*0.5/volume*farea*fnz;

            //! dQre/dgradR*dgradRdqL
            RDouble dGrad = dGdqx * dxR + dGdqy * dyR + dGdqz * dzR;

            //! convert dGrad by right multiplying the dqdcvL,
            //! then add to the corresponding location in dRdq.
            indexf = le * nEquation;    //! first index
            indexs = le * nEquation;    //! second index
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, le);
            int colidx2 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, le);
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        dRdq[indexI][colidx+indexJ] += dGrad* dFluxdpR[indexI][indexK] * dqdcvL[indexK][indexJ];
                        dRdq[indexI][colidx2+indexJ] -= dGrad* dFluxdpR[indexI][indexK] * dqdcvL[indexK][indexJ];
                    }
                }
            }

            //! considering the perturbation of the gradient in the limiter
            volume = vol[le];

            //! dgradLdqR
            dGdqx = 1.0*0.5/volume*farea*fnx;
            dGdqy = 1.0*0.5/volume*farea*fny;
            dGdqz = 1.0*0.5/volume*farea*fnz;

            //! dQle/dgradL*dgradLdqR
            dGrad = dGdqx * dxL + dGdqy * dyL + dGdqz * dzL;

            //! convert dGrad by right multiplying the dqdcvL,
            //! then add to the corresponding location in dRdq.

            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, re);
            colidx2 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, re);
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        dRdq[indexI][colidx+indexJ] += dGrad* dFluxdpL[indexI][indexK] * dqdcvR[indexK][indexJ];
                        dRdq[indexI][colidx2+indexJ] -= dGrad* dFluxdpL[indexI][indexK] * dqdcvR[indexK][indexJ];
                    }
                }
            }

            //! the neighbors of the left cells
            for(int index = 0; index < neighborCells[le].size(); index++)
            {
                int neighborIndex = neighborCells[le][index];

                if(neighborIndex != re)
                {
                    int Faceindex  = neighborFaces[le][index];
                    int sign       = neighborLR[le][index];
                    RDouble fnx    = nsx[Faceindex];
                    RDouble fny    = nsy[Faceindex];
                    RDouble fnz    = nsz[Faceindex];
                    RDouble farea  = ns[Faceindex];
                    RDouble volume = vol[le];

                    for(int m =0; m < nEquation; m++)
                    {
                        primN[m] = prim[m][neighborIndex];
                    }

                    gas->dPrimitive2dConservative(primN,gama[neighborIndex],dqdcvN);

                    //! dgradLdqN
                    dGdqx = sign*0.5/volume*farea*fnx;
                    dGdqy = sign*0.5/volume*farea*fny;
                    dGdqz = sign*0.5/volume*farea*fnz;

                    //! dQle/dgradL*dgradLdqN
                    dGrad = dGdqx * dxL + dGdqy * dyL + dGdqz * dzL;

                    int indexf  = le * nEquation; // first index
                    int indexs  = neighborIndex * nEquation; // second index
                    int indexf2 = re * nEquation; 
                    int colidx  = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, neighborIndex);
                    int colidx2 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, neighborIndex);
                    for(int indexI=0; indexI < nEquation; indexI++)
                    {
                        for(int indexJ=0; indexJ < nEquation; indexJ++)
                        {

                            for(int indexK=0; indexK < nEquation; indexK++)
                            {
                                //! dRL/dqLNeighbors 
                                dRdq[indexI][colidx+indexJ] += dGrad* dFluxdpL[indexI][indexK]*
                                                                         dqdcvN[indexK][indexJ];

                                //! dRR/dqLNeighbors 
                                dRdq[indexI][colidx2+indexJ] -= dGrad* dFluxdpL[indexI][indexK]*
                                                                         dqdcvN[indexK][indexJ];
                            }
                        }
                    }
                }
            }

            //! the neighbors of the right cells
            for(int index = 0; index < neighborCells[re].size(); index++)
            {
                int neighborIndex = neighborCells[re][index];

                if(neighborIndex != le)
                {
                    int Faceindex  = neighborFaces[re][index];
                    int sign       = neighborLR[re][index];
                    RDouble fnx    = nsx[Faceindex];
                    RDouble fny    = nsy[Faceindex];
                    RDouble fnz    = nsz[Faceindex];
                    RDouble farea  = ns[Faceindex];
                    RDouble volume = vol[re];

                    for(int m =0; m < nEquation; m++)
                    {
                        primN[m] = prim[m][neighborIndex];
                    }

                    gas->dPrimitive2dConservative(primN,gama[neighborIndex],dqdcvN);

                    //! dgradRdqN
                    dGdqx = sign*0.5/volume*farea*fnx;
                    dGdqy = sign*0.5/volume*farea*fny;
                    dGdqz = sign*0.5/volume*farea*fnz;

                    //! dQre/dgradR*dgradRdqN
                    dGrad = dGdqx * dxR + dGdqy * dyR + dGdqz * dzR;

                    int indexf  = le * nEquation;    //! first index
                    int indexs  = neighborIndex * nEquation;    //! second index
                    int indexf2 = re * nEquation; 
                    int colidx  = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, neighborIndex);
                    int colidx2 = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, neighborIndex);
                    for(int indexI=0; indexI < nEquation; indexI++)
                    {
                        for(int indexJ=0; indexJ < nEquation; indexJ++)
                        {

                            for(int indexK=0; indexK < nEquation; indexK++)
                            {
                                //! dRL/dqRNeighbors 
                                dRdq[indexI][colidx+indexJ] += dGrad* dFluxdpR[indexI][indexK]*
                                                                         dqdcvN[indexK][indexJ];

                                //! dRR/dqRNeighbors
                                dRdq[indexI][colidx2+indexJ] -= dGrad* dFluxdpR[indexI][indexK]*
                                                                         dqdcvN[indexK][indexJ];
                            }
                        }
                    }
                }
            }
        }
     }

    delete [] primL;    primL = nullptr;
    delete [] primR;    primR = nullptr;
    delete [] primN;    primN = nullptr;
    delete [] primLC;   primLC = nullptr;
    delete [] primRC;   primRC = nullptr;
    delete [] fluxp;    fluxp = nullptr;

    for(int index = 0; index < nEquation; index++)
    {
        delete [] dqdcvL[index];    dqdcvL[index] = nullptr;
        delete [] dqdcvR[index];    dqdcvR[index] = nullptr;
        delete [] dqdcvN[index];    dqdcvN[index] = nullptr;
        delete [] dgraddqx[index];    dgraddqx[index] = nullptr;
        delete [] dgraddqy[index];    dgraddqy[index] = nullptr;
        delete [] dgraddqz[index];    dgraddqz[index] = nullptr;
    }
    delete [] dqdcvL;    dqdcvL = nullptr;
    delete [] dqdcvR;    dqdcvR = nullptr;
    delete [] dqdcvN;    dqdcvN = nullptr;
    delete [] dgraddqx;    dgraddqx = nullptr;
    delete [] dgraddqy;    dgraddqy = nullptr;
    delete [] dgraddqz;    dgraddqz = nullptr;
}

//! GMRES FD CSR format  GMRESCSR
void GMRES_Roe_Scheme_ConservativeForm_CSR(FaceProxy *face_proxy, InviscidSchemeParameter *invSchemePara)
{
    const RDouble LITTLE    = 1.0E-10;
    int nTemperatureModel   = invSchemePara->GetNumberOfTemperatureModel();
    int nEquation           = invSchemePara->GetNumberOfTotalEquation();
    int nLaminar            = invSchemePara->GetNumberOfLaminar();
    int nLength             = invSchemePara->GetLength();
    int nm                  = invSchemePara->Getnm();
    int RoeEntropyFixMethod = invSchemePara->GetEntropyFixMethod();
    int nChemical           = invSchemePara->GetNumberOfChemical();
    int nPrecondition       = invSchemePara->GetIfPrecondition();
    int nIdealState         = invSchemePara->GetIfIdealState();

    //! GMRESBoundary
    int localStart = face_proxy->GetlocalStart();
    int localEnd   = face_proxy->GetlocalEnd();
    int nBoundFace = face_proxy->GetnBoundFace();
    int nTotalCell = face_proxy->GetnTotalCell();
    int nMid;
    //! GMRESBCorrection
    vector<int> *wallFaceIndex = face_proxy->GetWallFaceIndex();

    RDouble refMachNumber     = invSchemePara->GetMachNumber();
    RDouble entrFixCoefLeft   = invSchemePara->GetLeftEntropyFixCoefficient();
    RDouble entrFixCoefRight  = invSchemePara->GetRightEntropyFixCoefficient();
    RDouble **qL              = invSchemePara->GetLeftQ();
    RDouble **qR              = invSchemePara->GetRightQ();
    RDouble **qLC             = invSchemePara->GetLeftQC();    //! GMRESPassQC
    RDouble **qRC             = invSchemePara->GetRightQC();   //! GMRESPassQC
    RDouble *gamaL            = invSchemePara->GetLeftGama();
    RDouble *gamaR            = invSchemePara->GetRightGama();
    RDouble *pressureCoeffL   = face_proxy->GetPressureCoefficientL();
    RDouble *pressureCoeffR   = face_proxy->GetPressureCoefficientR();
    RDouble *xfn              = invSchemePara->GetFaceNormalX();
    RDouble *yfn              = invSchemePara->GetFaceNormalY();
    RDouble *zfn              = invSchemePara->GetFaceNormalZ();
    RDouble *vgn              = invSchemePara->GetFaceVelocity();
    RDouble **flux            = invSchemePara->GetFlux();
    RDouble **tl              = invSchemePara->GetLeftTemperature();
    RDouble **tr              = invSchemePara->GetRightTemperature();
    int *leftcellindexofFace  = face_proxy->GetLeftCellIndexOfFace();
    int *rightcellindexofFace = face_proxy->GetRightCellIndexOfFace();
    GeomProxy *geomProxy      = face_proxy->GetGeomProxy();
    RDouble *areas            = geomProxy->GetFaceArea();
    RDouble **dRdq            = face_proxy->GetJacobianMatrix();
    //! GMRESBoundary
    RDouble **dDdP = face_proxy->GetdDdPMatrix();
    //! GMRES CSR
    vector<int> AI = face_proxy->GetJacobianAI4GMRES();
    vector<int> AJ = face_proxy->GetJacobianAJ4GMRES();
    // RDouble   dRdqL[nEquation][nEquation];
    RDouble dRdqR[nEquation][nEquation];
    RDouble dDdPlocal[nEquation][nEquation];
    RDouble mach2 = refMachNumber * refMachNumber;

    using namespace IDX;
    using namespace GAS_SPACE;

    RDouble *primLpert = new RDouble[nEquation]();
    RDouble *primRpert = new RDouble[nEquation]();
    RDouble dFluxdpL[nEquation][nEquation];
    RDouble dFluxdpR[nEquation][nEquation];
    RDouble dFluxdpwall[nEquation][nEquation];
    RDouble dprimpert;
    RDouble *primL  = new RDouble[nEquation]();
    RDouble *primR  = new RDouble[nEquation]();
    RDouble *primLC = new RDouble[nEquation]();    //! GMRESPassQC
    RDouble *primRC = new RDouble[nEquation]();    //! GMRESPassQC
    RDouble *fluxp  = new RDouble[nEquation]();    //! flux proxy
    RDouble temperaturesL[3] = { 0 }, temperaturesR[3] = { 0 };
    RDouble perturbScale = 0.001;

    RDouble **dqdcvL= new RDouble*[nEquation];
    RDouble **dqdcvR= new RDouble*[nEquation];
    for(int indexI=0; indexI < nEquation; indexI++)
    {
        dqdcvL[indexI] = new RDouble[nEquation];
        dqdcvR[indexI] = new RDouble[nEquation];

        for(int indexJ =0; indexJ < nEquation; indexJ ++)
        {
            dqdcvL[indexI][indexJ] = 0.0;
            dqdcvR[indexI][indexJ] = 0.0;
        }
    }

    //! GMRESBoundary
     if (localStart >= nBoundFace)
     {
         nMid = localStart;    //! a bug 01.07
     }
     else if (localEnd <= nBoundFace)
     {
         //! If they are all boundary faces.
         nMid = localEnd;
     }
     else
     {
         //! Part of them are boundary faces.
         nMid = nBoundFace;
     }

     nMid = nMid - localStart;

     for(int iFace = 0;iFace < nMid; ++iFace)
     {
        int le = leftcellindexofFace[iFace];
        int re = rightcellindexofFace[iFace];
        int colidx;
        RDouble area  = areas[iFace];
        for (int iEqn = 0; iEqn < nEquation; ++iEqn)
        {
            primL[iEqn] = qL[iEqn][iFace];
            primR[iEqn] = qR[iEqn][iFace];

            primLC[iEqn] = qLC[iEqn][iFace];    //! GMRESPassQC
            primRC[iEqn] = qRC[iEqn][iFace];    //! GMRESPassQC
        }

        RDouble gmL = gamaL[iFace];
        RDouble gmR = gamaR[iFace];

        Cal_GMRES_Roe_Scheme_Flux_PV(primL, primR, iFace, fluxp,
                                  face_proxy,invSchemePara);

        flux[IR][iFace] = fluxp[IR];
        flux[IU][iFace] = fluxp[IU];
        flux[IV][iFace] = fluxp[IV];
        flux[IW][iFace] = fluxp[IW];
        flux[IP][iFace] = fluxp[IP];

        //! perturb the density of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            primLpert[index] = primL[index];
        }
        primLpert[IR] *= (1.0+perturbScale);
        dprimpert      = primLpert[IR] - primL[IR];

        Cal_GMRES_Roe_Scheme_Flux_PV(primLpert, primR, iFace, fluxp,
                                  face_proxy,invSchemePara);

        for(int index = 0; index < nEquation; index++)
        {
            dFluxdpL[index][IR] = (fluxp[index] - flux[index][iFace]) / dprimpert*area;
        }
        

        //! perturb the u of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            primLpert[index] = primL[index];
        }
        primLpert[IU] += perturbScale;
        dprimpert      = primLpert[IU] - primL[IU];

        Cal_GMRES_Roe_Scheme_Flux_PV(primLpert, primR, iFace, fluxp,
                                        face_proxy,invSchemePara);

        for(int index = 0; index < nEquation; index++)
        {
            dFluxdpL[index][IU] = (fluxp[index] - flux[index][iFace]) / dprimpert*area;
        }

        //! perturb the v of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            primLpert[index] = primL[index];
        }
        primLpert[IV] += perturbScale;
        dprimpert      = primLpert[IV] - primL[IV];

        Cal_GMRES_Roe_Scheme_Flux_PV(primLpert, primR, iFace, fluxp,
                                        face_proxy,invSchemePara);

        for(int index = 0; index < nEquation; index++)
        {
            dFluxdpL[index][IV] = (fluxp[index] - flux[index][iFace]) / dprimpert*area;
        }

        //! perturb the w of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            primLpert[index] = primL[index];
        }
        primLpert[IW] += perturbScale;
        dprimpert      = primLpert[IW] - primL[IW];

        Cal_GMRES_Roe_Scheme_Flux_PV(primLpert, primR, iFace, fluxp,
                                        face_proxy,invSchemePara);

        for(int index = 0; index < nEquation; index++)
        {
            dFluxdpL[index][IW] = (fluxp[index] - flux[index][iFace]) / dprimpert*area;
        }

        //! perturb the P of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            primLpert[index] = primL[index];
        }
        primLpert[IP] *= (1.0+perturbScale);
        dprimpert     = primLpert[IP] - primL[IP];

        Cal_GMRES_Roe_Scheme_Flux_PV(primLpert, primR, iFace, fluxp,
                                        face_proxy,invSchemePara);

        for(int index = 0; index < nEquation; index++)
        {
            dFluxdpL[index][IP] = (fluxp[index] - flux[index][iFace]) / dprimpert*area;
        }

        //! perturb the rho of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            primRpert[index] = primR[index];
        }
        primRpert[IR] *= (1.0+perturbScale);
        dprimpert      = primRpert[IR] - primR[IR];

        Cal_GMRES_Roe_Scheme_Flux_PV(primL, primRpert, iFace, fluxp,
                                        face_proxy,invSchemePara);

        for(int index = 0; index < nEquation; index++)
        {
            dFluxdpR[index][IR] = (fluxp[index] - flux[index][iFace]) / dprimpert*area;
        }

        //! perturb the u of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            primRpert[index] = primR[index];
        }
        primRpert[IU] += perturbScale;
        dprimpert      = primRpert[IU] - primR[IU];

        Cal_GMRES_Roe_Scheme_Flux_PV(primL, primRpert, iFace, fluxp,
                                        face_proxy,invSchemePara);

        for(int index = 0; index < nEquation; index++)
        {
            dFluxdpR[index][IU] = (fluxp[index] - flux[index][iFace]) / dprimpert*area;
        }

        //! perturb the v of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            primRpert[index] = primR[index];
        }
        primRpert[IV] += perturbScale;
        dprimpert      = primRpert[IV] - primR[IV];

        Cal_GMRES_Roe_Scheme_Flux_PV(primL, primRpert, iFace, fluxp,
                                        face_proxy,invSchemePara);

        for(int index = 0; index < nEquation; index++)
        {
            dFluxdpR[index][IV] = (fluxp[index] - flux[index][iFace]) / dprimpert*area;
        }

        //! perturb the w of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            primRpert[index] = primR[index];
        }
        primRpert[IW] += perturbScale;
        dprimpert      = primRpert[IW] - primR[IW];

        Cal_GMRES_Roe_Scheme_Flux_PV(primL, primRpert, iFace, fluxp,
                                        face_proxy,invSchemePara);

        for(int index = 0; index < nEquation; index++)
        {
            dFluxdpR[index][IW] = (fluxp[index] - flux[index][iFace]) / dprimpert*area;
        }

        //! perturb the P of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            primRpert[index] = primR[index];
        }
        primRpert[IP] *= (1.0+perturbScale);
        dprimpert      = primRpert[IP] - primR[IP];

        Cal_GMRES_Roe_Scheme_Flux_PV(primL, primRpert, iFace, fluxp,
                                        face_proxy,invSchemePara);

        for(int index = 0; index < nEquation; index++)
        {
            dFluxdpR[index][IP] = (fluxp[index] - flux[index][iFace]) / dprimpert*area;
        }

        //! contributions from the bc
        int indexre = (re-nTotalCell)*nEquation;
        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                dDdPlocal[indexI][indexJ] = dDdP[indexI][indexre+indexJ];
            }
        }

        //! GMRESBCorrection
        if(find(wallFaceIndex->begin(),wallFaceIndex->end(),iFace) != wallFaceIndex->end())
        {
            //! GMRESPV , obtain the dqdcv for the left cell
            gas->dPrimitive2dConservative(primLC,gmL,dqdcvL);

            //! average of the dFluxdpL and dFluxdpR
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                
                for(int indexJ =0; indexJ < nEquation; indexJ ++)
                {
                    dFluxdpwall[indexI][indexJ] = 0.5*(dFluxdpL[indexI][indexJ]
                                                     + dFluxdpR[indexI][indexJ]);
                    dFluxdpL[indexI][indexJ]    =  dFluxdpwall[indexI][indexJ];
                    dFluxdpR[indexI][indexJ]    =  dFluxdpwall[indexI][indexJ];
                }
            }

            //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
            //! then minus to the corresponding location in dRdq.
            int indexf = re * nEquation;    //! first index
            int indexs = le * nEquation;    //! second index
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, le);
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        dRdq[indexI][colidx+indexJ] -= dFluxdpwall[indexI][indexK]*
                                                            dqdcvL[indexK][indexJ];
                    }
                }
            }

            //! contributions from the bc, 
            //! dFluxdpL = dFluxdpL + dFluxdPR*dDdPlocal
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        dFluxdpL[indexI][indexJ]  += dFluxdpwall[indexI][indexK]*
                                                       dDdPlocal[indexK][indexJ];
                    }
                }
            }

            //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
            //! then add to the corresponding location in dRdq.
            indexf = le * nEquation;    //! first index
            indexs = le * nEquation;    //! second index
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, le);
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        dRdq[indexI][colidx+indexJ] += dFluxdpL[indexI][indexK]*
                                                         dqdcvL[indexK][indexJ];
                    }
                }
            }

            //! GMRESPV , obtain the dqdcv for the right cell
            gas->dPrimitive2dConservative(primRC,gmR,dqdcvR);

            //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
            //! then minus to the corresponding location in dRdq.
            indexf = re * nEquation;    //! first index
            indexs = re * nEquation;    //! second index
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, re);
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        dRdq[indexI][colidx+indexJ] -= dFluxdpR[indexI][indexK]*
                                                         dqdcvR[indexK][indexJ];
                    }
                }
            }

            //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
            //! then add to the corresponding location in dRdq.
            indexf = le * nEquation;    //! first index
            indexs = re * nEquation;    //! second index
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, re);
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        dRdq[indexI][colidx+indexJ] += dFluxdpR[indexI][indexK]*
                                                         dqdcvR[indexK][indexJ];
                    }
                }
            }
        }
        else
        {
            //! GMRESPV , obtain the dqdcv for the left cell
            gas->dPrimitive2dConservative(primLC,gmL,dqdcvL);

            //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
            //! then minus to the corresponding location in dRdq.
            int indexf = re * nEquation;    //! first index
            int indexs = le * nEquation;    //! second index
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, le);
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        dRdq[indexI][colidx+indexJ] -= dFluxdpL[indexI][indexK]*
                                                         dqdcvL[indexK][indexJ];
                    }
                }
            }

            //! contributions from the bc, 
            //! dFluxdpL = dFluxdpL + dFluxdPR*dDdPlocal
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        dFluxdpL[indexI][indexJ]  += dFluxdpR[indexI][indexK]*
                                                    dDdPlocal[indexK][indexJ];
                    }
                }
            }

            //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
            //! then add to the corresponding location in dRdq.
            indexf = le * nEquation;    //! first index
            indexs = le * nEquation;    //! second index
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, le);
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {

                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        dRdq[indexI][colidx+indexJ] += dFluxdpL[indexI][indexK]*
                                                         dqdcvL[indexK][indexJ];
                    }
                }
            }

            //! GMRESPV , obtain the dqdcv for the right cell
            gas->dPrimitive2dConservative(primRC,gmR,dqdcvR);

            //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
            //! then minus to the corresponding location in dRdq.
            indexf = re * nEquation;    //! first index
            indexs = re * nEquation;    //! second index
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, re);
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        dRdq[indexI][colidx+indexJ] -= dFluxdpR[indexI][indexK]*
                                                         dqdcvR[indexK][indexJ];
                    }
                }
            }

            //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
            //! then add to the corresponding location in dRdq.
            indexf = le * nEquation;    //! first index
            indexs = re * nEquation;    //! second index
            colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, re);
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        dRdq[indexI][colidx+indexJ] += dFluxdpR[indexI][indexK]*
                                                         dqdcvR[indexK][indexJ];
                    }
                }
            }
        }
     }

     for(int iFace = nMid;iFace < nLength; ++iFace)
     {
        int le = leftcellindexofFace[iFace];
        int re = rightcellindexofFace[iFace];
        int colidx;
        RDouble area  = areas[iFace];
        for (int iEqn = 0; iEqn < nEquation; ++iEqn)
        {
            primL[iEqn] = qL[iEqn][iFace];
            primR[iEqn] = qR[iEqn][iFace];

            primLC[iEqn] = qLC[iEqn][iFace];    //! GMRESPassQC
            primRC[iEqn] = qRC[iEqn][iFace];    //! GMRESPassQC
        }

        RDouble gmL = gamaL[iFace];
        RDouble gmR = gamaR[iFace];

        Cal_GMRES_Roe_Scheme_Flux_PV(primL, primR, iFace, fluxp,
                                  face_proxy,invSchemePara);

        flux[IR][iFace] = fluxp[IR];
        flux[IU][iFace] = fluxp[IU];
        flux[IV][iFace] = fluxp[IV];
        flux[IW][iFace] = fluxp[IW];
        flux[IP][iFace] = fluxp[IP];

        //! perturb the density of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            primLpert[index] = primL[index];
        }
        primLpert[IR] *= (1.0+perturbScale);
        dprimpert     = primLpert[IR] - primL[IR];

        Cal_GMRES_Roe_Scheme_Flux_PV(primLpert, primR, iFace, fluxp,
                                  face_proxy,invSchemePara);

        for(int index = 0; index < nEquation; index++)
        {
            dFluxdpL[index][IR] = (fluxp[index] - flux[index][iFace]) / dprimpert*area;
        }

        //! perturb the u of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            primLpert[index] = primL[index];
        }
        primLpert[IU] += perturbScale;
        dprimpert      = primLpert[IU] - primL[IU];

        Cal_GMRES_Roe_Scheme_Flux_PV(primLpert, primR, iFace, fluxp,
                                        face_proxy,invSchemePara);

        for(int index = 0; index < nEquation; index++)
        {
            dFluxdpL[index][IU] = (fluxp[index] - flux[index][iFace]) / dprimpert*area;
        }

        //! perturb the v of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            primLpert[index] = primL[index];
        }
        primLpert[IV] += perturbScale;
        dprimpert      = primLpert[IV] - primL[IV];

        Cal_GMRES_Roe_Scheme_Flux_PV(primLpert, primR, iFace, fluxp,
                                        face_proxy,invSchemePara);

        for(int index = 0; index < nEquation; index++)
        {
            dFluxdpL[index][IV] = (fluxp[index] - flux[index][iFace]) / dprimpert*area;
        }

        //! perturb the w of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            primLpert[index] = primL[index];
        }
        primLpert[IW] += perturbScale;
        dprimpert      = primLpert[IW] - primL[IW];

        Cal_GMRES_Roe_Scheme_Flux_PV(primLpert, primR, iFace, fluxp,
                                        face_proxy,invSchemePara);

        for(int index = 0; index < nEquation; index++)
        {
            dFluxdpL[index][IW] = (fluxp[index] - flux[index][iFace]) / dprimpert*area;
        }

        //! perturb the P of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            primLpert[index] = primL[index];
        }
        primLpert[IP] *= (1.0+perturbScale);
        dprimpert      = primLpert[IP] - primL[IP];

        Cal_GMRES_Roe_Scheme_Flux_PV(primLpert, primR, iFace, fluxp,
                                        face_proxy,invSchemePara);

        for(int index = 0; index < nEquation; index++)
        {
            dFluxdpL[index][IP] = (fluxp[index] - flux[index][iFace]) / dprimpert*area;
        }

        //! perturb the rho of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            primRpert[index] = primR[index];
        }
        primRpert[IR] *= (1.0+perturbScale);
        dprimpert      = primRpert[IR] - primR[IR];

        Cal_GMRES_Roe_Scheme_Flux_PV(primL, primRpert, iFace, fluxp,
                                        face_proxy,invSchemePara);

        for(int index = 0; index < nEquation; index++)
        {
            dFluxdpR[index][IR] = (fluxp[index] - flux[index][iFace])/dprimpert*area;
        }

        //! perturb the u of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            primRpert[index] = primR[index];
        }
        primRpert[IU] += perturbScale;
        dprimpert      = primRpert[IU] - primR[IU];

        Cal_GMRES_Roe_Scheme_Flux_PV(primL, primRpert, iFace, fluxp,
                                        face_proxy,invSchemePara);

        for(int index = 0; index < nEquation; index++)
        {
            dFluxdpR[index][IU] = (fluxp[index] - flux[index][iFace]) / dprimpert*area;
        }

        //! perturb the v of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            primRpert[index] = primR[index];
        }
        primRpert[IV] += perturbScale;
        dprimpert      = primRpert[IV] - primR[IV];

        Cal_GMRES_Roe_Scheme_Flux_PV(primL, primRpert, iFace, fluxp,
                                        face_proxy,invSchemePara);

        for(int index = 0; index < nEquation; index++)
        {
            dFluxdpR[index][IV] = (fluxp[index] - flux[index][iFace]) / dprimpert*area;
        }

        //! perturb the w of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            primRpert[index] = primR[index];
        }
        primRpert[IW] += perturbScale;
        dprimpert      = primRpert[IW] - primR[IW];

        Cal_GMRES_Roe_Scheme_Flux_PV(primL, primRpert, iFace, fluxp,
                                        face_proxy,invSchemePara);

        for(int index = 0; index < nEquation; index++)
        {
            dFluxdpR[index][IW] = (fluxp[index] - flux[index][iFace]) / dprimpert*area;
        }

        //! perturb the P of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            primRpert[index] = primR[index];
        }
        primRpert[IP] *= (1.0+perturbScale);
        dprimpert      = primRpert[IP] - primR[IP];

        Cal_GMRES_Roe_Scheme_Flux_PV(primL, primRpert, iFace, fluxp,
                                        face_proxy,invSchemePara);

        for(int index = 0; index < nEquation; index++)
        {
            dFluxdpR[index][IP] = (fluxp[index] - flux[index][iFace]) / dprimpert*area;
        }

        //! GMRESPV , obtain the dqdcv for the left cell
        gas->dPrimitive2dConservative(primLC,gmL,dqdcvL);

        //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
        //! then minus to the corresponding location in dRdq.
        int indexf = re * nEquation;    //! first index
        int indexs = le * nEquation;    //! second index
        colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, le);
        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                for(int indexK=0; indexK < nEquation; indexK++)
                {
                    dRdq[indexI][colidx+indexJ] -= dFluxdpL[indexI][indexK] *
                                                     dqdcvL[indexK][indexJ];
                }
            }
        }

        //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
        //! then add to the corresponding location in dRdq.
        indexf = le * nEquation;    //! first index
        indexs = le * nEquation;    //! second index
        colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, le);
        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                for(int indexK=0; indexK < nEquation; indexK++)
                {
                    dRdq[indexI][colidx+indexJ] += dFluxdpL[indexI][indexK] *
                                                     dqdcvL[indexK][indexJ];
                }
            }
        }

        //! GMRESPV , obtain the dqdcv for the right cell
        gas->dPrimitive2dConservative(primRC,gmR,dqdcvR);

        //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
        //! then minus to the corresponding location in dRdq.
        indexf = re * nEquation;    //! first index
        indexs = re * nEquation;    //! second index
        colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, re, re);
        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                
                for(int indexK=0; indexK < nEquation; indexK++)
                {
                    dRdq[indexI][colidx+indexJ] -= dFluxdpR[indexI][indexK]*
                                                             dqdcvR[indexK][indexJ];
                }
            }
        }

        //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
        //! then add to the corresponding location in dRdq.
        indexf = le * nEquation;    //! first index
        indexs = re * nEquation;    //! second index
        colidx = GetIndexOfBlockJacobianMatrix(AI, AJ, nEquation, le, re);
        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                for(int indexK=0; indexK < nEquation; indexK++)
                {
                    dRdq[indexI][colidx+indexJ] += dFluxdpR[indexI][indexK] *
                                                      dqdcvR[indexK][indexJ];
                }
            }
        }
     }

    delete [] primL;        primL = nullptr;
    delete [] primR;        primR = nullptr;
    delete [] primLC;       primLC = nullptr;
    delete [] primRC;       primRC = nullptr;
    delete [] primLpert;    primLpert = nullptr;
    delete [] primRpert;    primRpert = nullptr;
    delete [] fluxp;        fluxp = nullptr;

    for(int index = 0; index < nEquation; index++)
    {
        delete [] dqdcvL[index];    dqdcvL[index] = nullptr;
        delete [] dqdcvR[index];    dqdcvR[index] = nullptr;
    }
    delete [] dqdcvL;    dqdcvL = nullptr;
    delete [] dqdcvR;    dqdcvR = nullptr;
}

//! GMRESPV
//! GMRES for finite difference considering nBoundFace using primitive variables
void GMRES_Roe_Scheme_ConservativeForm_FD_PV(FaceProxy *face_proxy, InviscidSchemeParameter *invSchemePara)
{
    const RDouble LITTLE    = 1.0E-10;
    int nTemperatureModel   = invSchemePara->GetNumberOfTemperatureModel();
    int nEquation           = invSchemePara->GetNumberOfTotalEquation();
    int nLaminar            = invSchemePara->GetNumberOfLaminar();
    int nLength             = invSchemePara->GetLength();
    int nm                  = invSchemePara->Getnm();
    int RoeEntropyFixMethod = invSchemePara->GetEntropyFixMethod();
    int nChemical           = invSchemePara->GetNumberOfChemical();
    int nPrecondition       = invSchemePara->GetIfPrecondition();
    int nIdealState         = invSchemePara->GetIfIdealState();

    //! GMRESBoundary
    int localStart = face_proxy->GetlocalStart();
    int localEnd   = face_proxy->GetlocalEnd();
    int nBoundFace = face_proxy->GetnBoundFace();
    int nTotalCell = face_proxy->GetnTotalCell();
    int nMid;
    //! GMRESBCorrection
    vector<int> *wallFaceIndex = face_proxy->GetWallFaceIndex();

    RDouble refMachNumber     = invSchemePara->GetMachNumber();
    RDouble entrFixCoefLeft   = invSchemePara->GetLeftEntropyFixCoefficient();
    RDouble entrFixCoefRight  = invSchemePara->GetRightEntropyFixCoefficient();
    RDouble **qL              = invSchemePara->GetLeftQ();
    RDouble **qR              = invSchemePara->GetRightQ();
    RDouble **qLC             = invSchemePara->GetLeftQC();    //! GMRESPassQC
    RDouble **qRC             = invSchemePara->GetRightQC();   //! GMRESPassQC
    RDouble *gamaL            = invSchemePara->GetLeftGama();
    RDouble *gamaR            = invSchemePara->GetRightGama();
    RDouble *pressureCoeffL   = face_proxy->GetPressureCoefficientL();
    RDouble *pressureCoeffR   = face_proxy->GetPressureCoefficientR();
    RDouble *xfn              = invSchemePara->GetFaceNormalX();
    RDouble *yfn              = invSchemePara->GetFaceNormalY();
    RDouble *zfn              = invSchemePara->GetFaceNormalZ();
    RDouble *vgn              = invSchemePara->GetFaceVelocity();
    RDouble **flux            = invSchemePara->GetFlux();
    RDouble **tl              = invSchemePara->GetLeftTemperature();
    RDouble **tr              = invSchemePara->GetRightTemperature();
    int *leftcellindexofFace  = face_proxy->GetLeftCellIndexOfFace();
    int *rightcellindexofFace = face_proxy->GetRightCellIndexOfFace();
    GeomProxy *geomProxy      = face_proxy->GetGeomProxy();
    RDouble *areas            = geomProxy->GetFaceArea();
    RDouble **dRdq            = face_proxy->GetJacobianMatrix();
    //! GMRESBoundary
    RDouble **dDdP            = face_proxy->GetdDdPMatrix();
    //! RDouble   dRdqL[nEquation][nEquation];
    RDouble dRdqR[nEquation][nEquation];
    RDouble dDdPlocal[nEquation][nEquation];
    RDouble mach2 = refMachNumber * refMachNumber;

    using namespace IDX;
    using namespace GAS_SPACE;

    RDouble *primLpert = new RDouble[nEquation]();
    RDouble *primRpert = new RDouble[nEquation]();
    RDouble dFluxdpL[nEquation][nEquation];
    RDouble dFluxdpR[nEquation][nEquation];
    RDouble dFluxdpwall[nEquation][nEquation];
    RDouble dprimpert;
    RDouble *primL  = new RDouble[nEquation]();
    RDouble *primR  = new RDouble[nEquation]();
    RDouble *primLC = new RDouble[nEquation]();    //! GMRESPassQC
    RDouble *primRC = new RDouble[nEquation]();    //! GMRESPassQC
    RDouble *fluxp  = new RDouble[nEquation]();    //! flux proxy
    RDouble temperaturesL[3] = { 0 }, temperaturesR[3] = { 0 };
    RDouble perturbScale = 1e-6;

    RDouble **dqdcvL= new RDouble*[nEquation];
    RDouble **dqdcvR= new RDouble*[nEquation];
    for(int indexI=0; indexI < nEquation; indexI++)
    {
        dqdcvL[indexI] = new RDouble[nEquation];
        dqdcvR[indexI] = new RDouble[nEquation];

        for(int indexJ =0; indexJ < nEquation; indexJ ++)
        {
            dqdcvL[indexI][indexJ] = 0.0;
            dqdcvR[indexI][indexJ] = 0.0;
        }
    }

     //! GMRESBoundary
     if (localStart >= nBoundFace)
     {
         nMid = localStart;    //! a bug 01.07
     }
     else if (localEnd <= nBoundFace)    //! add else 
     {
         //! If they are all boundary faces.
         nMid = localEnd;
     }
     else
     {
         //! Part of them are boundary faces.
         nMid = nBoundFace;
     }

     nMid = nMid - localStart;

     for(int iFace = 0;iFace < nMid; ++iFace)
     {
        int le       = leftcellindexofFace[iFace];
        int re       = rightcellindexofFace[iFace];
        RDouble area = areas[iFace];
        for (int iEqn = 0; iEqn < nEquation; ++iEqn)
        {
            primL[iEqn] = qL[iEqn][iFace];
            primR[iEqn] = qR[iEqn][iFace];

            primLC[iEqn] = qLC[iEqn][iFace];    //! GMRESPassQC
            primRC[iEqn] = qRC[iEqn][iFace];    //! GMRESPassQC
        }

        RDouble gmL = gamaL[iFace];
        RDouble gmR = gamaR[iFace];

        Cal_GMRES_Roe_Scheme_Flux_PV(primL, primR, iFace, fluxp,
                                  face_proxy,invSchemePara);

        flux[IR][iFace] = fluxp[IR];
        flux[IU][iFace] = fluxp[IU];
        flux[IV][iFace] = fluxp[IV];
        flux[IW][iFace] = fluxp[IW];
        flux[IP][iFace] = fluxp[IP];

        //! perturb the density of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            primLpert[index] = primL[index];
        }
        primLpert[IR] *= (1.0+perturbScale);
        dprimpert      = primLpert[IR] - primL[IR];

        Cal_GMRES_Roe_Scheme_Flux_PV(primLpert, primR, iFace, fluxp,
                                  face_proxy,invSchemePara);

        for(int index = 0; index < nEquation; index++)
        {
            dFluxdpL[index][IR] = (fluxp[index] - flux[index][iFace]) / dprimpert*area;
        }

        //! perturb the u of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            primLpert[index] = primL[index];
        }
        primLpert[IU] += perturbScale;
        dprimpert      = primLpert[IU] - primL[IU];

        Cal_GMRES_Roe_Scheme_Flux_PV(primLpert, primR, iFace, fluxp,
                                        face_proxy,invSchemePara);

        for(int index = 0; index < nEquation; index++)
        {
            dFluxdpL[index][IU] = (fluxp[index] - flux[index][iFace]) / dprimpert*area;
        }

        //! perturb the v of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            primLpert[index] = primL[index];
        }
        primLpert[IV] += perturbScale;
        dprimpert      = primLpert[IV] - primL[IV];

        Cal_GMRES_Roe_Scheme_Flux_PV(primLpert, primR, iFace, fluxp,
                                        face_proxy,invSchemePara);

        for(int index = 0; index < nEquation; index++)
        {
            dFluxdpL[index][IV] = (fluxp[index] - flux[index][iFace]) / dprimpert*area;
        }

        //! perturb the w of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            primLpert[index] = primL[index];
        }
        primLpert[IW] += perturbScale;
        dprimpert      = primLpert[IW] - primL[IW];

        Cal_GMRES_Roe_Scheme_Flux_PV(primLpert, primR, iFace, fluxp,
                                        face_proxy,invSchemePara);

        for(int index = 0; index < nEquation; index++)
        {
            dFluxdpL[index][IW] = (fluxp[index] - flux[index][iFace]) / dprimpert*area;
        }

        //! perturb the P of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            primLpert[index] = primL[index];
        }
        primLpert[IP] *= (1.0+perturbScale);
        dprimpert      = primLpert[IP] - primL[IP];

        Cal_GMRES_Roe_Scheme_Flux_PV(primLpert, primR, iFace, fluxp,
                                        face_proxy,invSchemePara);

        for(int index = 0; index < nEquation; index++)
        {
            dFluxdpL[index][IP] = (fluxp[index] - flux[index][iFace]) / dprimpert*area;
        }

        //! perturb the rho of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            primRpert[index] = primR[index];
        }
        primRpert[IR] *= (1.0+perturbScale);
        dprimpert      = primRpert[IR] - primR[IR];

        Cal_GMRES_Roe_Scheme_Flux_PV(primL, primRpert, iFace, fluxp,
                                        face_proxy,invSchemePara);

        for(int index = 0; index < nEquation; index++)
        {
            dFluxdpR[index][IR] = (fluxp[index] - flux[index][iFace]) / dprimpert*area;
        }

        //! perturb the u of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            primRpert[index] = primR[index];
        }
        primRpert[IU] += perturbScale;
        dprimpert     = primRpert[IU] - primR[IU];

        Cal_GMRES_Roe_Scheme_Flux_PV(primL, primRpert, iFace, fluxp,
                                        face_proxy,invSchemePara);

        for(int index = 0; index < nEquation; index++)
        {
            dFluxdpR[index][IU] = (fluxp[index] - flux[index][iFace])/dprimpert*area;
        }

        //! perturb the v of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            primRpert[index] = primR[index];
        }
        primRpert[IV] += perturbScale;
        dprimpert      = primRpert[IV] - primR[IV];

        Cal_GMRES_Roe_Scheme_Flux_PV(primL, primRpert, iFace, fluxp,
                                        face_proxy,invSchemePara);

        for(int index = 0; index < nEquation; index++)
        {
            dFluxdpR[index][IV] = (fluxp[index] - flux[index][iFace]) / dprimpert*area;
        }

        //! perturb the w of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            primRpert[index] = primR[index];
        }
        primRpert[IW] += perturbScale;
        dprimpert     = primRpert[IW] - primR[IW];

        Cal_GMRES_Roe_Scheme_Flux_PV(primL, primRpert, iFace, fluxp,
                                        face_proxy,invSchemePara);

        for(int index = 0; index < nEquation; index++)
        {
            dFluxdpR[index][IW] = (fluxp[index] - flux[index][iFace]) / dprimpert*area;
        }

        //! perturb the P of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            primRpert[index] = primR[index];
        }
        primRpert[IP] *= (1.0+perturbScale);
        dprimpert      = primRpert[IP] - primR[IP];

        Cal_GMRES_Roe_Scheme_Flux_PV(primL, primRpert, iFace, fluxp,
                                        face_proxy,invSchemePara);

        for(int index = 0; index < nEquation; index++)
        {
            dFluxdpR[index][IP] = (fluxp[index] - flux[index][iFace]) / dprimpert*area;
        }

        //! contributions from the bc
        int indexre = (re-nTotalCell) * nEquation;
        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                dDdPlocal[indexI][indexJ] = dDdP[indexI][indexre+indexJ];
            }
        }

        //! GMRESBCorrection
        if(find(wallFaceIndex->begin(),wallFaceIndex->end(),iFace) != wallFaceIndex->end())
        {
            //! GMRESPV , obtain the dqdcv for the left cell
            gas->dPrimitive2dConservative(primLC,gmL,dqdcvL);

            //! average of the dFluxdpL and dFluxdpR
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                
                for(int indexJ =0; indexJ < nEquation; indexJ ++)
                {
                    dFluxdpwall[indexI][indexJ] = 0.5 * (dFluxdpL[indexI][indexJ]
                                                      + dFluxdpR[indexI][indexJ]);
                    dFluxdpL[indexI][indexJ]    =  dFluxdpwall[indexI][indexJ];
                    dFluxdpR[indexI][indexJ]    =  dFluxdpwall[indexI][indexJ];
                }
            }

            //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
            //! then minus to the corresponding location in dRdq.
            int indexf = re * nEquation;    //! first index
            int indexs = le * nEquation;    //! second index
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        dRdq[indexf+indexI][indexs+indexJ] -= dFluxdpwall[indexI][indexK]*
                                                                   dqdcvL[indexK][indexJ];
                    }
                }
            }

            //! contributions from the bc, 
            //! dFluxdpL = dFluxdpL + dFluxdPR*dDdPlocal
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        dFluxdpL[indexI][indexJ]  += dFluxdpwall[indexI][indexK]*
                                                       dDdPlocal[indexK][indexJ];
                    }
                }
            }

            //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
            //! then add to the corresponding location in dRdq.
            indexf = le * nEquation;    //! first index
            indexs = le * nEquation;    //! second index
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        dRdq[indexf+indexI][indexs+indexJ] += dFluxdpL[indexI][indexK]*
                                                                 dqdcvL[indexK][indexJ];
                    }
                }
            }

            //! GMRESPV , obtain the dqdcv for the right cell
            gas->dPrimitive2dConservative(primRC,gmR,dqdcvR);

            //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
            //! then minus to the corresponding location in dRdq.
            indexf = re * nEquation;    //! first index
            indexs = re * nEquation;    //! second index
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        dRdq[indexf+indexI][indexs+indexJ] -= dFluxdpR[indexI][indexK]*
                                                                dqdcvR[indexK][indexJ];
                    }
                }
            }

            //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
            //! then add to the corresponding location in dRdq.
            indexf = le * nEquation;    //! first index
            indexs = re * nEquation;    //! second index
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        dRdq[indexf+indexI][indexs+indexJ] += dFluxdpR[indexI][indexK]*
                                                                dqdcvR[indexK][indexJ];
                    }
                }
            }
        }
        else
        {
            //! GMRESPV , obtain the dqdcv for the left cell
            gas->dPrimitive2dConservative(primLC,gmL,dqdcvL);

            //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
            //! then minus to the corresponding location in dRdq.
            int indexf = re * nEquation;    //! first index
            int indexs = le * nEquation;    //! second index
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        dRdq[indexf+indexI][indexs+indexJ] -= dFluxdpL[indexI][indexK]*
                                                                dqdcvL[indexK][indexJ];
                    }
                }
            }

            //! contributions from the bc, 
            //! dFluxdpL = dFluxdpL + dFluxdPR*dDdPlocal
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        dFluxdpL[indexI][indexJ]  += dFluxdpR[indexI][indexK]*
                                                    dDdPlocal[indexK][indexJ];
                    }
                }
            }

            //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
            //! then add to the corresponding location in dRdq.
            indexf = le * nEquation;    //! first index
            indexs = le * nEquation;    //! second index
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        dRdq[indexf+indexI][indexs+indexJ] += dFluxdpL[indexI][indexK]*
                                                                dqdcvL[indexK][indexJ];
                    }
                }
            }

            //! GMRESPV , obtain the dqdcv for the right cell
            gas->dPrimitive2dConservative(primRC,gmR,dqdcvR);

            //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
            //! then minus to the corresponding location in dRdq.
            indexf = re * nEquation;    //! first index
            indexs = re * nEquation;    //! second index
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        dRdq[indexf+indexI][indexs+indexJ] -= dFluxdpR[indexI][indexK]*
                                                                dqdcvR[indexK][indexJ];
                    }
                }
            }

            //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
            //! then add to the corresponding location in dRdq.
            indexf = le * nEquation;    //! first index
            indexs = re * nEquation;    //! second index
            for(int indexI=0; indexI < nEquation; indexI++)
            {
                for(int indexJ=0; indexJ < nEquation; indexJ++)
                {
                    for(int indexK=0; indexK < nEquation; indexK++)
                    {
                        dRdq[indexf+indexI][indexs+indexJ] += dFluxdpR[indexI][indexK]*
                                                                dqdcvR[indexK][indexJ];
                    }
                }
            }
        }
     }

     for(int iFace = nMid;iFace < nLength; ++iFace)
     {
        int le       = leftcellindexofFace[iFace];
        int re       = rightcellindexofFace[iFace];
        RDouble area = areas[iFace];
        for (int iEqn = 0; iEqn < nEquation; ++iEqn)
        {
            primL[iEqn] = qL[iEqn][iFace];
            primR[iEqn] = qR[iEqn][iFace];

            primLC[iEqn] = qLC[iEqn][iFace];    //! GMRESPassQC
            primRC[iEqn] = qRC[iEqn][iFace];    //! GMRESPassQC
        }

        RDouble gmL = gamaL[iFace];
        RDouble gmR = gamaR[iFace];

        Cal_GMRES_Roe_Scheme_Flux_PV(primL, primR, iFace, fluxp,
                                  face_proxy,invSchemePara);

        flux[IR][iFace] = fluxp[IR];
        flux[IU][iFace] = fluxp[IU];
        flux[IV][iFace] = fluxp[IV];
        flux[IW][iFace] = fluxp[IW];
        flux[IP][iFace] = fluxp[IP];

        //! perturb the density of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            primLpert[index] = primL[index];
        }
        primLpert[IR] *= (1.0+perturbScale);
        dprimpert      = primLpert[IR] - primL[IR];

        Cal_GMRES_Roe_Scheme_Flux_PV(primLpert, primR, iFace, fluxp,
                                  face_proxy,invSchemePara);

        for(int index = 0; index < nEquation; index++)
        {
            dFluxdpL[index][IR] = (fluxp[index] - flux[index][iFace]) / dprimpert*area;
        }

        //! perturb the u of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            primLpert[index] = primL[index];
        }
        primLpert[IU] += perturbScale;
        dprimpert      = primLpert[IU] - primL[IU];

        Cal_GMRES_Roe_Scheme_Flux_PV(primLpert, primR, iFace, fluxp,
                                        face_proxy,invSchemePara);

        for(int index = 0; index < nEquation; index++)
        {
            dFluxdpL[index][IU] = (fluxp[index] - flux[index][iFace]) / dprimpert*area;
        }

        //! perturb the v of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            primLpert[index] = primL[index];
        }
        primLpert[IV] += perturbScale;
        dprimpert      = primLpert[IV] - primL[IV];

        Cal_GMRES_Roe_Scheme_Flux_PV(primLpert, primR, iFace, fluxp,
                                        face_proxy,invSchemePara);

        for(int index = 0; index < nEquation; index++)
        {
            dFluxdpL[index][IV] = (fluxp[index] - flux[index][iFace])/dprimpert*area;
        }

        //! perturb the w of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            primLpert[index] = primL[index];
        }
        primLpert[IW] += perturbScale;
        dprimpert      = primLpert[IW] - primL[IW];

        Cal_GMRES_Roe_Scheme_Flux_PV(primLpert, primR, iFace, fluxp,
                                        face_proxy,invSchemePara);

        for(int index = 0; index < nEquation; index++)
        {
            dFluxdpL[index][IW] = (fluxp[index] - flux[index][iFace]) / dprimpert*area;
        }

        //! perturb the P of the left cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            primLpert[index] = primL[index];
        }
        primLpert[IP] *= (1.0+perturbScale);
        dprimpert      = primLpert[IP] - primL[IP];

        Cal_GMRES_Roe_Scheme_Flux_PV(primLpert, primR, iFace, fluxp,
                                        face_proxy,invSchemePara);

        for(int index = 0; index < nEquation; index++)
        {
            dFluxdpL[index][IP] = (fluxp[index] - flux[index][iFace]) / dprimpert*area;
        }

        //! perturb the rho of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            primRpert[index] = primR[index];
        }
        primRpert[IR] *= (1.0+perturbScale);
        dprimpert      = primRpert[IR] - primR[IR];

        Cal_GMRES_Roe_Scheme_Flux_PV(primL, primRpert, iFace, fluxp,
                                        face_proxy,invSchemePara);

        for(int index = 0; index < nEquation; index++)
        {
            dFluxdpR[index][IR] = (fluxp[index] - flux[index][iFace]) / dprimpert*area;
        }

        //! perturb the u of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            primRpert[index] = primR[index];
        }
        primRpert[IU] += perturbScale;
        dprimpert     = primRpert[IU] - primR[IU];

        Cal_GMRES_Roe_Scheme_Flux_PV(primL, primRpert, iFace, fluxp,
                                        face_proxy,invSchemePara);

        for(int index = 0; index < nEquation; index++)
        {
            dFluxdpR[index][IU] = (fluxp[index] - flux[index][iFace]) / dprimpert*area;
        }

        //! perturb the v of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            primRpert[index] = primR[index];
        }
        primRpert[IV] += perturbScale;
        dprimpert      = primRpert[IV] - primR[IV];

        Cal_GMRES_Roe_Scheme_Flux_PV(primL, primRpert, iFace, fluxp,
                                        face_proxy,invSchemePara);

        for(int index = 0; index < nEquation; index++)
        {
            dFluxdpR[index][IV] = (fluxp[index] - flux[index][iFace]) / dprimpert*area;
        }

        //! perturb the w of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            primRpert[index] = primR[index];
        }
        primRpert[IW] += perturbScale;
        dprimpert      = primRpert[IW] - primR[IW];

        Cal_GMRES_Roe_Scheme_Flux_PV(primL, primRpert, iFace, fluxp,
                                        face_proxy,invSchemePara);

        for(int index = 0; index < nEquation; index++)
        {
            dFluxdpR[index][IW] = (fluxp[index] - flux[index][iFace]) / dprimpert*area;
        }

        //! perturb the P of the right cell of the face
        for(int index = 0; index < nEquation; index++)
        {
            primRpert[index] = primR[index];
        }
        primRpert[IP] *= (1.0+perturbScale);
        dprimpert      = primRpert[IP] - primR[IP];

        Cal_GMRES_Roe_Scheme_Flux_PV(primL, primRpert, iFace, fluxp,
                                        face_proxy,invSchemePara);

        for(int index = 0; index < nEquation; index++)
        {
            dFluxdpR[index][IP] = (fluxp[index] - flux[index][iFace]) / dprimpert*area;
        }

        //! GMRESPV , obtain the dqdcv for the left cell
        gas->dPrimitive2dConservative(primLC,gmL,dqdcvL);

        //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
        //! then minus to the corresponding location in dRdq.
        int indexf = re * nEquation;    //! first index
        int indexs = le * nEquation;    //! second index
        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                for(int indexK=0; indexK < nEquation; indexK++)
                {
                    dRdq[indexf+indexI][indexs+indexJ] -= dFluxdpL[indexI][indexK]*
                                                             dqdcvL[indexK][indexJ];
                }
            }
        }

        //! convert dFluxdpL to dFluxdcvL by right multiplying the dqdcvL,
        //! then add to the corresponding location in dRdq.
        indexf = le * nEquation;    //! first index
        indexs = le * nEquation;    //! second index
        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                for(int indexK=0; indexK < nEquation; indexK++)
                {
                    dRdq[indexf+indexI][indexs+indexJ] += dFluxdpL[indexI][indexK]*
                                                             dqdcvL[indexK][indexJ];
                }
            }
        }

        //! GMRESPV , obtain the dqdcv for the right cell
        gas->dPrimitive2dConservative(primRC,gmR,dqdcvR);

        //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
        //! then minus to the corresponding location in dRdq.
        indexf = re * nEquation;    //! first index
        indexs = re * nEquation;    //! second index
        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                for(int indexK=0; indexK < nEquation; indexK++)
                {
                    dRdq[indexf+indexI][indexs+indexJ] -= dFluxdpR[indexI][indexK]*
                                                            dqdcvR[indexK][indexJ];
                }
            }
        }

        //! convert dFluxdpR to dFluxdcvR by right multiplying the dqdcvR,
        //! then add to the corresponding location in dRdq.
        indexf = le * nEquation;    //! first index
        indexs = re * nEquation;    //! second index
        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                for(int indexK=0; indexK < nEquation; indexK++)
                {
                    dRdq[indexf+indexI][indexs+indexJ] += dFluxdpR[indexI][indexK]*
                                                            dqdcvR[indexK][indexJ];
                }
            }
        }
     }

    delete [] primL;        primL = nullptr;
    delete [] primR;        primR = nullptr;
    delete [] primLC;       primLC = nullptr;
    delete [] primRC;       primRC = nullptr;
    delete [] primLpert;    primLpert = nullptr;
    delete [] primRpert;    primRpert = nullptr;
    delete [] fluxp;        fluxp = nullptr;

    for(int index = 0; index < nEquation; index++)
    {
        delete [] dqdcvL[index];    dqdcvL[index] = nullptr;
        delete [] dqdcvR[index];    dqdcvR[index] = nullptr;
    }
    delete [] dqdcvL;    dqdcvL = nullptr;
    delete [] dqdcvR;    dqdcvR = nullptr;
}

//! GMRES for finite difference considering nBoundFace using conservative variables
void GMRES_Roe_Scheme_ConservativeForm_CV(FaceProxy *face_proxy, InviscidSchemeParameter *invSchemePara)
{
    const RDouble LITTLE    = 1.0E-10;
    int nTemperatureModel   = invSchemePara->GetNumberOfTemperatureModel();
    int nEquation           = invSchemePara->GetNumberOfTotalEquation();
    int nLaminar            = invSchemePara->GetNumberOfLaminar();
    int nLength             = invSchemePara->GetLength();
    int nm                  = invSchemePara->Getnm();
    int RoeEntropyFixMethod = invSchemePara->GetEntropyFixMethod();
    int nChemical           = invSchemePara->GetNumberOfChemical();
    int nPrecondition       = invSchemePara->GetIfPrecondition();
    int nIdealState         = invSchemePara->GetIfIdealState();

    //! MRESBoundary
    int localStart = face_proxy->GetlocalStart();
    int localEnd   = face_proxy->GetlocalEnd();
    int nBoundFace = face_proxy->GetnBoundFace();
    int nTotalCell = face_proxy->GetnTotalCell();
    int nMid;

    RDouble refMachNumber     = invSchemePara->GetMachNumber();
    RDouble entrFixCoefLeft   = invSchemePara->GetLeftEntropyFixCoefficient();
    RDouble entrFixCoefRight  = invSchemePara->GetRightEntropyFixCoefficient();
    RDouble **qL              = invSchemePara->GetLeftQ();
    RDouble **qR              = invSchemePara->GetRightQ();
    RDouble *gamaL            = invSchemePara->GetLeftGama();
    RDouble *gamaR            = invSchemePara->GetRightGama();
    RDouble *pressureCoeffL   = face_proxy->GetPressureCoefficientL();
    RDouble *pressureCoeffR   = face_proxy->GetPressureCoefficientR();
    RDouble *xfn              = invSchemePara->GetFaceNormalX();
    RDouble *yfn              = invSchemePara->GetFaceNormalY();
    RDouble *zfn              = invSchemePara->GetFaceNormalZ();
    RDouble *vgn              = invSchemePara->GetFaceVelocity();
    RDouble **flux            = invSchemePara->GetFlux();
    RDouble **tl              = invSchemePara->GetLeftTemperature();
    RDouble **tr              = invSchemePara->GetRightTemperature();
    int *leftcellindexofFace  = face_proxy->GetLeftCellIndexOfFace();
    int *rightcellindexofFace = face_proxy->GetRightCellIndexOfFace();
    GeomProxy *geomProxy      = face_proxy->GetGeomProxy();
    RDouble *areas            = geomProxy->GetFaceArea();
    RDouble **dRdq            = face_proxy->GetJacobianMatrix();
    //! GMRESBoundary
    RDouble** dDdP = face_proxy->GetdDdPMatrix();
    //! RDouble   dRdqL[nEquation][nEquation];
    RDouble dRdqR[nEquation][nEquation];
    RDouble dDdPlocal[nEquation][nEquation];
    RDouble mach2 = refMachNumber * refMachNumber;

    using namespace IDX;
    using namespace GAS_SPACE;

    RDouble *prim       = new RDouble[nEquation]();
    RDouble *qConserveL = new RDouble[nEquation]();
    RDouble *qConserveR = new RDouble[nEquation]();
    RDouble *qCvL       = new RDouble[nEquation]();
    RDouble *qCvR       = new RDouble[nEquation]();
    RDouble *qCvLpert   = new RDouble[nEquation]();
    RDouble *qCvRpert   = new RDouble[nEquation]();
    RDouble *dFluxdCv   = new RDouble[nEquation]();
    RDouble dCvpert;

    RDouble *primL = new RDouble[nEquation]();
    RDouble *primR = new RDouble[nEquation]();
    RDouble *fluxp = new RDouble[nEquation]();    //! flux proxy
    RDouble temperaturesL[3] = { 0 }, temperaturesR[3] = { 0 };
    RDouble perturbScale = 0.001;

     //! GMRESBoundary
     if (localStart >= nBoundFace)
     {
         nMid = 0;
     }
     if (localEnd <= nBoundFace)
     {
         //! If they are all boundary faces.
         nMid = localEnd;
     }
     else
     {
         //! Part of them are boundary faces.
         nMid = nBoundFace;
     }

     nMid = nMid - localStart;

     for(int iFace = 0;iFace < nMid; ++iFace)
     {
        int le       = leftcellindexofFace[iFace];
        int re       = rightcellindexofFace[iFace];
        RDouble area = areas[iFace];
        for (int iEqn = 0; iEqn < nEquation; ++iEqn)
        {
            primL[iEqn] = qL[iEqn][iFace];
            primR[iEqn] = qR[iEqn][iFace];
        }

        RDouble gmL = gamaL[iFace];
        RDouble gmR = gamaR[iFace];

        gas->Primitive2Conservative(primL, gmL, temperaturesL[ITV], temperaturesL[ITE], qCvL);
        gas->Primitive2Conservative(primR, gmR, temperaturesR[ITV], temperaturesR[ITE], qCvR);

        gas->Conservative2Primitive(qCvL, gmL, primL, temperaturesL);
        gas->Conservative2Primitive(qCvR, gmR, primR, temperaturesR);

        Cal_GMRES_Roe_Scheme_Flux(primL, primR, qCvL, qCvR,iFace, fluxp,
                                  face_proxy,invSchemePara);

        flux[IR][iFace] = fluxp[IR];
        flux[IU][iFace] = fluxp[IU];
        flux[IV][iFace] = fluxp[IV];
        flux[IW][iFace] = fluxp[IW];
        flux[IP][iFace] = fluxp[IP];

        //! perturb the density of the left cell of the face
        qCvLpert[IR] = qCvL[IR];
        qCvLpert[IU] = qCvL[IU];
        qCvLpert[IV] = qCvL[IV];
        qCvLpert[IW] = qCvL[IW];
        qCvLpert[IP] = qCvL[IP];
        qCvLpert[IR] *= (1.0+perturbScale);
        dCvpert      = qCvLpert[IR] - qCvL[IR];

        gas->Conservative2Primitive(qCvLpert, gmL, primL, temperaturesL);
        gas->Conservative2Primitive(qCvR, gmR, primR, temperaturesR);

        Cal_GMRES_Roe_Scheme_Flux(primL, primR, qCvLpert, qCvR, iFace, fluxp,
                                  face_proxy, invSchemePara);

        dFluxdCv[IR] = (fluxp[IR] - flux[IR][iFace]) / dCvpert;
        dFluxdCv[IU] = (fluxp[IU] - flux[IU][iFace]) / dCvpert;
        dFluxdCv[IV] = (fluxp[IV] - flux[IV][iFace]) / dCvpert;
        dFluxdCv[IW] = (fluxp[IW] - flux[IW][iFace]) / dCvpert;
        dFluxdCv[IP] = (fluxp[IP] - flux[IP][iFace]) / dCvpert;

        dRdq[le * nEquation + IR][le * nEquation + IR] += dFluxdCv[IR] * area;
        dRdq[le * nEquation + IU][le * nEquation + IR] += dFluxdCv[IU] * area;
        dRdq[le * nEquation + IV][le * nEquation + IR] += dFluxdCv[IV] * area;
        dRdq[le * nEquation + IW][le * nEquation + IR] += dFluxdCv[IW] * area;
        dRdq[le * nEquation + IP][le * nEquation + IR] += dFluxdCv[IP] * area;

        // for(int index=0; index < nEquation; index++)
        // {
        //     dRdqL[index][IR] = dFluxdCv[index] * area;
        // }

        dRdq[re * nEquation + IR][le * nEquation + IR] -= dFluxdCv[IR] * area;
        dRdq[re * nEquation + IU][le * nEquation + IR] -= dFluxdCv[IU] * area;
        dRdq[re * nEquation + IV][le * nEquation + IR] -= dFluxdCv[IV] * area;
        dRdq[re * nEquation + IW][le * nEquation + IR] -= dFluxdCv[IW] * area;
        dRdq[re * nEquation + IP][le * nEquation + IR] -= dFluxdCv[IP] * area;

        //! perturb the momentum rhou of the left cell of the face
        qCvLpert[IR] = qCvL[IR];
        qCvLpert[IU] = qCvL[IU];
        qCvLpert[IV] = qCvL[IV];
        qCvLpert[IW] = qCvL[IW];
        qCvLpert[IP] = qCvL[IP];
        qCvLpert[IU] += perturbScale;
        dCvpert = qCvLpert[IU] - qCvL[IU];

        gas->Conservative2Primitive(qCvLpert, gmL, primL, temperaturesL);
        gas->Conservative2Primitive(qCvR, gmR, primR, temperaturesR);

        Cal_GMRES_Roe_Scheme_Flux(primL, primR, qCvLpert, qCvR, iFace, fluxp,
                                  face_proxy, invSchemePara);

        dFluxdCv[IR] = (fluxp[IR] - flux[IR][iFace]) / dCvpert;
        dFluxdCv[IU] = (fluxp[IU] - flux[IU][iFace]) / dCvpert;
        dFluxdCv[IV] = (fluxp[IV] - flux[IV][iFace]) / dCvpert;
        dFluxdCv[IW] = (fluxp[IW] - flux[IW][iFace]) / dCvpert;
        dFluxdCv[IP] = (fluxp[IP] - flux[IP][iFace]) / dCvpert;

        dRdq[le * nEquation + IR][le * nEquation + IU] += dFluxdCv[IR] * area;
        dRdq[le * nEquation + IU][le * nEquation + IU] += dFluxdCv[IU] * area;
        dRdq[le * nEquation + IV][le * nEquation + IU] += dFluxdCv[IV] * area;
        dRdq[le * nEquation + IW][le * nEquation + IU] += dFluxdCv[IW] * area;
        dRdq[le * nEquation + IP][le * nEquation + IU] += dFluxdCv[IP] * area;

        // for(int index=0; index < nEquation; index++)
        // {
        //     dRdqL[index][IU] = dFluxdCv[index] * area;
        // }

        dRdq[re * nEquation + IR][le * nEquation + IU] -= dFluxdCv[IR] * area;
        dRdq[re * nEquation + IU][le * nEquation + IU] -= dFluxdCv[IU] * area;
        dRdq[re * nEquation + IV][le * nEquation + IU] -= dFluxdCv[IV] * area;
        dRdq[re * nEquation + IW][le * nEquation + IU] -= dFluxdCv[IW] * area;
        dRdq[re * nEquation + IP][le * nEquation + IU] -= dFluxdCv[IP] * area;

        //! perturb the momentum rhov of the left cell of the face
        qCvLpert[IR] = qCvL[IR];
        qCvLpert[IU] = qCvL[IU];
        qCvLpert[IV] = qCvL[IV];
        qCvLpert[IW] = qCvL[IW];
        qCvLpert[IP] = qCvL[IP];
        qCvLpert[IV] += perturbScale;
        dCvpert = qCvLpert[IV] - qCvL[IV];

        gas->Conservative2Primitive(qCvLpert, gmL, primL, temperaturesL);
        gas->Conservative2Primitive(qCvR, gmR, primR, temperaturesR);

        Cal_GMRES_Roe_Scheme_Flux(primL, primR, qCvLpert, qCvR, iFace, fluxp,
            face_proxy, invSchemePara);

        dFluxdCv[IR] = (fluxp[IR] - flux[IR][iFace]) / dCvpert;
        dFluxdCv[IU] = (fluxp[IU] - flux[IU][iFace]) / dCvpert;
        dFluxdCv[IV] = (fluxp[IV] - flux[IV][iFace]) / dCvpert;
        dFluxdCv[IW] = (fluxp[IW] - flux[IW][iFace]) / dCvpert;
        dFluxdCv[IP] = (fluxp[IP] - flux[IP][iFace]) / dCvpert;

        dRdq[le * nEquation + IR][le * nEquation + IV] += dFluxdCv[IR] * area;
        dRdq[le * nEquation + IU][le * nEquation + IV] += dFluxdCv[IU] * area;
        dRdq[le * nEquation + IV][le * nEquation + IV] += dFluxdCv[IV] * area;
        dRdq[le * nEquation + IW][le * nEquation + IV] += dFluxdCv[IW] * area;
        dRdq[le * nEquation + IP][le * nEquation + IV] += dFluxdCv[IP] * area;

        // for(int index=0; index < nEquation; index++)
        // {
        //     dRdqL[index][IV] = dFluxdCv[index] * area;
        // }

        dRdq[re * nEquation + IR][le * nEquation + IV] -= dFluxdCv[IR] * area;
        dRdq[re * nEquation + IU][le * nEquation + IV] -= dFluxdCv[IU] * area;
        dRdq[re * nEquation + IV][le * nEquation + IV] -= dFluxdCv[IV] * area;
        dRdq[re * nEquation + IW][le * nEquation + IV] -= dFluxdCv[IW] * area;
        dRdq[re * nEquation + IP][le * nEquation + IV] -= dFluxdCv[IP] * area;

        //! perturb the momentum rhow of the left cell of the face
        qCvLpert[IR] = qCvL[IR];
        qCvLpert[IU] = qCvL[IU];
        qCvLpert[IV] = qCvL[IV];
        qCvLpert[IW] = qCvL[IW];
        qCvLpert[IP] = qCvL[IP];
        qCvLpert[IW] += perturbScale;
        dCvpert = qCvLpert[IW] - qCvL[IW];

        gas->Conservative2Primitive(qCvLpert, gmL, primL, temperaturesL);
        gas->Conservative2Primitive(qCvR, gmR, primR, temperaturesR);

        Cal_GMRES_Roe_Scheme_Flux(primL, primR, qCvLpert, qCvR, iFace, fluxp,
            face_proxy, invSchemePara);

        dFluxdCv[IR] = (fluxp[IR] - flux[IR][iFace]) / dCvpert;
        dFluxdCv[IU] = (fluxp[IU] - flux[IU][iFace]) / dCvpert;
        dFluxdCv[IV] = (fluxp[IV] - flux[IV][iFace]) / dCvpert;
        dFluxdCv[IW] = (fluxp[IW] - flux[IW][iFace]) / dCvpert;
        dFluxdCv[IP] = (fluxp[IP] - flux[IP][iFace]) / dCvpert;

        dRdq[le * nEquation + IR][le * nEquation + IW] += dFluxdCv[IR] * area;
        dRdq[le * nEquation + IU][le * nEquation + IW] += dFluxdCv[IU] * area;
        dRdq[le * nEquation + IV][le * nEquation + IW] += dFluxdCv[IV] * area;
        dRdq[le * nEquation + IW][le * nEquation + IW] += dFluxdCv[IW] * area;
        dRdq[le * nEquation + IP][le * nEquation + IW] += dFluxdCv[IP] * area;

        // for(int index=0; index < nEquation; index++)
        // {
        //     dRdqL[index][IW] = dFluxdCv[index] * area;
        // }

        dRdq[re * nEquation + IR][le * nEquation + IW] -= dFluxdCv[IR] * area;
        dRdq[re * nEquation + IU][le * nEquation + IW] -= dFluxdCv[IU] * area;
        dRdq[re * nEquation + IV][le * nEquation + IW] -= dFluxdCv[IV] * area;
        dRdq[re * nEquation + IW][le * nEquation + IW] -= dFluxdCv[IW] * area;
        dRdq[re * nEquation + IP][le * nEquation + IW] -= dFluxdCv[IP] * area;

        //! perturb the rhoE of the left cell of the face
        qCvLpert[IR] = qCvL[IR];
        qCvLpert[IU] = qCvL[IU];
        qCvLpert[IV] = qCvL[IV];
        qCvLpert[IW] = qCvL[IW];
        qCvLpert[IP] = qCvL[IP];
        qCvLpert[IP] *= (1.0+perturbScale);
        dCvpert = qCvLpert[IP] - qCvL[IP];

        gas->Conservative2Primitive(qCvLpert, gmL, primL, temperaturesL);
        gas->Conservative2Primitive(qCvR, gmR, primR, temperaturesR);

        Cal_GMRES_Roe_Scheme_Flux(primL, primR, qCvLpert, qCvR, iFace, fluxp,
            face_proxy, invSchemePara);

        dFluxdCv[IR] = (fluxp[IR] - flux[IR][iFace]) / dCvpert;
        dFluxdCv[IU] = (fluxp[IU] - flux[IU][iFace]) / dCvpert;
        dFluxdCv[IV] = (fluxp[IV] - flux[IV][iFace]) / dCvpert;
        dFluxdCv[IW] = (fluxp[IW] - flux[IW][iFace]) / dCvpert;
        dFluxdCv[IP] = (fluxp[IP] - flux[IP][iFace]) / dCvpert;

        dRdq[le * nEquation + IR][le * nEquation + IP] += dFluxdCv[IR] * area;
        dRdq[le * nEquation + IU][le * nEquation + IP] += dFluxdCv[IU] * area;
        dRdq[le * nEquation + IV][le * nEquation + IP] += dFluxdCv[IV] * area;
        dRdq[le * nEquation + IW][le * nEquation + IP] += dFluxdCv[IW] * area;
        dRdq[le * nEquation + IP][le * nEquation + IP] += dFluxdCv[IP] * area;

        // for(int index=0; index < nEquation; index++)
        // {
        //     dRdqL[index][IP] = dFluxdCv[index] * area;
        // }

        dRdq[re * nEquation + IR][le * nEquation + IP] -= dFluxdCv[IR] * area;
        dRdq[re * nEquation + IU][le * nEquation + IP] -= dFluxdCv[IU] * area;
        dRdq[re * nEquation + IV][le * nEquation + IP] -= dFluxdCv[IV] * area;
        dRdq[re * nEquation + IW][le * nEquation + IP] -= dFluxdCv[IW] * area;
        dRdq[re * nEquation + IP][le * nEquation + IP] -= dFluxdCv[IP] * area;

        //! perturb the density of the right cell of the face
        qCvRpert[IR] = qCvR[IR];
        qCvRpert[IU] = qCvR[IU];
        qCvRpert[IV] = qCvR[IV];
        qCvRpert[IW] = qCvR[IW];
        qCvRpert[IP] = qCvR[IP];
        qCvRpert[IR] *= (1.0+perturbScale);
        dCvpert = qCvRpert[IR] - qCvR[IR];

        gas->Conservative2Primitive(qCvL, gmL, primL, temperaturesL);
        gas->Conservative2Primitive(qCvRpert, gmR, primR, temperaturesR);

        Cal_GMRES_Roe_Scheme_Flux(primL, primR, qCvL, qCvRpert, iFace, fluxp,
            face_proxy, invSchemePara);

        dFluxdCv[IR] = (fluxp[IR] - flux[IR][iFace]) / dCvpert;
        dFluxdCv[IU] = (fluxp[IU] - flux[IU][iFace]) / dCvpert;
        dFluxdCv[IV] = (fluxp[IV] - flux[IV][iFace]) / dCvpert;
        dFluxdCv[IW] = (fluxp[IW] - flux[IW][iFace]) / dCvpert;
        dFluxdCv[IP] = (fluxp[IP] - flux[IP][iFace]) / dCvpert;

        dRdq[le * nEquation + IR][re * nEquation + IR] += dFluxdCv[IR] * area;
        dRdq[le * nEquation + IU][re * nEquation + IR] += dFluxdCv[IU] * area;
        dRdq[le * nEquation + IV][re * nEquation + IR] += dFluxdCv[IV] * area;
        dRdq[le * nEquation + IW][re * nEquation + IR] += dFluxdCv[IW] * area;
        dRdq[le * nEquation + IP][re * nEquation + IR] += dFluxdCv[IP] * area;

        for(int index=0; index < nEquation; index++)
        {
            dRdqR[index][IR] = dFluxdCv[index] * area;
        }

        dRdq[re * nEquation + IR][re * nEquation + IR] -= dFluxdCv[IR] * area;
        dRdq[re * nEquation + IU][re * nEquation + IR] -= dFluxdCv[IU] * area;
        dRdq[re * nEquation + IV][re * nEquation + IR] -= dFluxdCv[IV] * area;
        dRdq[re * nEquation + IW][re * nEquation + IR] -= dFluxdCv[IW] * area;
        dRdq[re * nEquation + IP][re * nEquation + IR] -= dFluxdCv[IP] * area;

        //! perturb the momentum rhou of the right cell of the face
        qCvRpert[IR] = qCvR[IR];
        qCvRpert[IU] = qCvR[IU];
        qCvRpert[IV] = qCvR[IV];
        qCvRpert[IW] = qCvR[IW];
        qCvRpert[IP] = qCvR[IP];
        qCvRpert[IU] += perturbScale;
        dCvpert = qCvRpert[IU] - qCvR[IU];

        gas->Conservative2Primitive(qCvL, gmL, primL, temperaturesL);
        gas->Conservative2Primitive(qCvRpert, gmR, primR, temperaturesR);

        Cal_GMRES_Roe_Scheme_Flux(primL, primR, qCvL, qCvRpert, iFace, fluxp,
                                  face_proxy, invSchemePara);

        dFluxdCv[IR] = (fluxp[IR] - flux[IR][iFace]) / dCvpert;
        dFluxdCv[IU] = (fluxp[IU] - flux[IU][iFace]) / dCvpert;
        dFluxdCv[IV] = (fluxp[IV] - flux[IV][iFace]) / dCvpert;
        dFluxdCv[IW] = (fluxp[IW] - flux[IW][iFace]) / dCvpert;
        dFluxdCv[IP] = (fluxp[IP] - flux[IP][iFace]) / dCvpert;

        dRdq[le * nEquation + IR][re * nEquation + IU] += dFluxdCv[IR] * area;
        dRdq[le * nEquation + IU][re * nEquation + IU] += dFluxdCv[IU] * area;
        dRdq[le * nEquation + IV][re * nEquation + IU] += dFluxdCv[IV] * area;
        dRdq[le * nEquation + IW][re * nEquation + IU] += dFluxdCv[IW] * area;
        dRdq[le * nEquation + IP][re * nEquation + IU] += dFluxdCv[IP] * area;

        for(int index=0; index < nEquation; index++)
        {
            dRdqR[index][IU] = dFluxdCv[index] * area;
        }

        dRdq[re * nEquation + IR][re * nEquation + IU] -= dFluxdCv[IR] * area;
        dRdq[re * nEquation + IU][re * nEquation + IU] -= dFluxdCv[IU] * area;
        dRdq[re * nEquation + IV][re * nEquation + IU] -= dFluxdCv[IV] * area;
        dRdq[re * nEquation + IW][re * nEquation + IU] -= dFluxdCv[IW] * area;
        dRdq[re * nEquation + IP][re * nEquation + IU] -= dFluxdCv[IP] * area;

        //! perturb the momentum rhov of the right cell of the face
        qCvRpert[IR] = qCvR[IR];
        qCvRpert[IU] = qCvR[IU];
        qCvRpert[IV] = qCvR[IV];
        qCvRpert[IW] = qCvR[IW];
        qCvRpert[IP] = qCvR[IP];
        qCvRpert[IV] += perturbScale;
        dCvpert = qCvRpert[IV] - qCvR[IV];

        gas->Conservative2Primitive(qCvL, gmL, primL, temperaturesL);
        gas->Conservative2Primitive(qCvRpert, gmR, primR, temperaturesR);

        Cal_GMRES_Roe_Scheme_Flux(primL, primR, qCvL, qCvRpert, iFace, fluxp,
                                  face_proxy, invSchemePara);

        dFluxdCv[IR] = (fluxp[IR] - flux[IR][iFace]) / dCvpert;
        dFluxdCv[IU] = (fluxp[IU] - flux[IU][iFace]) / dCvpert;
        dFluxdCv[IV] = (fluxp[IV] - flux[IV][iFace]) / dCvpert;
        dFluxdCv[IW] = (fluxp[IW] - flux[IW][iFace]) / dCvpert;
        dFluxdCv[IP] = (fluxp[IP] - flux[IP][iFace]) / dCvpert;

        dRdq[le * nEquation + IR][re * nEquation + IV] += dFluxdCv[IR] * area;
        dRdq[le * nEquation + IU][re * nEquation + IV] += dFluxdCv[IU] * area;
        dRdq[le * nEquation + IV][re * nEquation + IV] += dFluxdCv[IV] * area;
        dRdq[le * nEquation + IW][re * nEquation + IV] += dFluxdCv[IW] * area;
        dRdq[le * nEquation + IP][re * nEquation + IV] += dFluxdCv[IP] * area;

        for(int index=0; index < nEquation; index++)
        {
            dRdqR[index][IV] = dFluxdCv[index] * area;
        }

        dRdq[re * nEquation + IR][re * nEquation + IV] -= dFluxdCv[IR] * area;
        dRdq[re * nEquation + IU][re * nEquation + IV] -= dFluxdCv[IU] * area;
        dRdq[re * nEquation + IV][re * nEquation + IV] -= dFluxdCv[IV] * area;
        dRdq[re * nEquation + IW][re * nEquation + IV] -= dFluxdCv[IW] * area;
        dRdq[re * nEquation + IP][re * nEquation + IV] -= dFluxdCv[IP] * area;

        //! perturb the momentum rhow of the right cell of the face
        qCvRpert[IR] = qCvR[IR];
        qCvRpert[IU] = qCvR[IU];
        qCvRpert[IV] = qCvR[IV];
        qCvRpert[IW] = qCvR[IW];
        qCvRpert[IP] = qCvR[IP];
        qCvRpert[IW] += perturbScale;
        dCvpert = qCvRpert[IW] - qCvR[IW];

        gas->Conservative2Primitive(qCvL, gmL, primL, temperaturesL);
        gas->Conservative2Primitive(qCvRpert, gmR, primR, temperaturesR);

        Cal_GMRES_Roe_Scheme_Flux(primL, primR, qCvL, qCvRpert, iFace, fluxp,
            face_proxy, invSchemePara);

        dFluxdCv[IR] = (fluxp[IR] - flux[IR][iFace]) / dCvpert;
        dFluxdCv[IU] = (fluxp[IU] - flux[IU][iFace]) / dCvpert;
        dFluxdCv[IV] = (fluxp[IV] - flux[IV][iFace]) / dCvpert;
        dFluxdCv[IW] = (fluxp[IW] - flux[IW][iFace]) / dCvpert;
        dFluxdCv[IP] = (fluxp[IP] - flux[IP][iFace]) / dCvpert;

        dRdq[le * nEquation + IR][re * nEquation + IW] += dFluxdCv[IR] * area;
        dRdq[le * nEquation + IU][re * nEquation + IW] += dFluxdCv[IU] * area;
        dRdq[le * nEquation + IV][re * nEquation + IW] += dFluxdCv[IV] * area;
        dRdq[le * nEquation + IW][re * nEquation + IW] += dFluxdCv[IW] * area;
        dRdq[le * nEquation + IP][re * nEquation + IW] += dFluxdCv[IP] * area;

        for(int index=0; index < nEquation; index++)
        {
            dRdqR[index][IW] = dFluxdCv[index] * area;
        }

        dRdq[re * nEquation + IR][re * nEquation + IW] -= dFluxdCv[IR] * area;
        dRdq[re * nEquation + IU][re * nEquation + IW] -= dFluxdCv[IU] * area;
        dRdq[re * nEquation + IV][re * nEquation + IW] -= dFluxdCv[IV] * area;
        dRdq[re * nEquation + IW][re * nEquation + IW] -= dFluxdCv[IW] * area;
        dRdq[re * nEquation + IP][re * nEquation + IW] -= dFluxdCv[IP] * area;

        //! perturb the momentum rhoE of the right cell of the face
        qCvRpert[IR] = qCvR[IR];
        qCvRpert[IU] = qCvR[IU];
        qCvRpert[IV] = qCvR[IV];
        qCvRpert[IW] = qCvR[IW];
        qCvRpert[IP] = qCvR[IP];
        qCvRpert[IP] *= (1.0+perturbScale);
        dCvpert = qCvRpert[IP] - qCvR[IP];

        gas->Conservative2Primitive(qCvL, gmL, primL, temperaturesL);
        gas->Conservative2Primitive(qCvRpert, gmR, primR, temperaturesR);

        Cal_GMRES_Roe_Scheme_Flux(primL, primR, qCvL, qCvRpert, iFace, fluxp,
            face_proxy, invSchemePara);

        dFluxdCv[IR] = (fluxp[IR] - flux[IR][iFace]) / dCvpert;
        dFluxdCv[IU] = (fluxp[IU] - flux[IU][iFace]) / dCvpert;
        dFluxdCv[IV] = (fluxp[IV] - flux[IV][iFace]) / dCvpert;
        dFluxdCv[IW] = (fluxp[IW] - flux[IW][iFace]) / dCvpert;
        dFluxdCv[IP] = (fluxp[IP] - flux[IP][iFace]) / dCvpert;

        dRdq[le * nEquation + IR][re * nEquation + IP] += dFluxdCv[IR] * area;
        dRdq[le * nEquation + IU][re * nEquation + IP] += dFluxdCv[IU] * area;
        dRdq[le * nEquation + IV][re * nEquation + IP] += dFluxdCv[IV] * area;
        dRdq[le * nEquation + IW][re * nEquation + IP] += dFluxdCv[IW] * area;
        dRdq[le * nEquation + IP][re * nEquation + IP] += dFluxdCv[IP] * area;

        for(int index=0; index < nEquation; index++)
        {
            dRdqR[index][IP] = dFluxdCv[index] * area;
        }

        dRdq[re * nEquation + IR][re * nEquation + IP] -= dFluxdCv[IR] * area;
        dRdq[re * nEquation + IU][re * nEquation + IP] -= dFluxdCv[IU] * area;
        dRdq[re * nEquation + IV][re * nEquation + IP] -= dFluxdCv[IV] * area;
        dRdq[re * nEquation + IW][re * nEquation + IP] -= dFluxdCv[IW] * area;
        dRdq[re * nEquation + IP][re * nEquation + IP] -= dFluxdCv[IP] * area;

        int indexre = (re-nTotalCell)*nEquation;
        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                dDdPlocal[indexI][indexJ] = dDdP[indexI][indexre+indexJ];
            }
        }

        int indexle = le * nEquation;
        for(int indexI=0; indexI < nEquation; indexI++)
        {
            for(int indexJ=0; indexJ < nEquation; indexJ++)
            {
                
                for(int indexK=0; indexK < nEquation; indexK++)
                {
                    dRdq[indexle+indexI][indexle+indexJ] += dRdqR[indexI][indexK]*
                                                             dDdPlocal[indexK][indexJ];
                }
            }
        }
     }

     for(int iFace = nMid;iFace < nLength; ++iFace)
     {
        int le       = leftcellindexofFace[iFace];
        int re       = rightcellindexofFace[iFace];
        RDouble area = areas[iFace];
        for (int iEqn = 0; iEqn < nEquation; ++iEqn)
        {
            primL[iEqn] = qL[iEqn][iFace];
            primR[iEqn] = qR[iEqn][iFace];
        }

        RDouble gmL = gamaL[iFace];
        RDouble gmR = gamaR[iFace];

        gas->Primitive2Conservative(primL, gmL, temperaturesL[ITV], temperaturesL[ITE], qCvL);
        gas->Primitive2Conservative(primR, gmR, temperaturesR[ITV], temperaturesR[ITE], qCvR);

        gas->Conservative2Primitive(qCvL, gmL, primL, temperaturesL);
        gas->Conservative2Primitive(qCvR, gmR, primR, temperaturesR);

        Cal_GMRES_Roe_Scheme_Flux(primL, primR, qCvL, qCvR,iFace, fluxp,
                                  face_proxy,invSchemePara);

        flux[IR][iFace] = fluxp[IR];
        flux[IU][iFace] = fluxp[IU];
        flux[IV][iFace] = fluxp[IV];
        flux[IW][iFace] = fluxp[IW];
        flux[IP][iFace] = fluxp[IP];

        //! perturb the density of the left cell of the face
        qCvLpert[IR] = qCvL[IR];
        qCvLpert[IU] = qCvL[IU];
        qCvLpert[IV] = qCvL[IV];
        qCvLpert[IW] = qCvL[IW];
        qCvLpert[IP] = qCvL[IP];
        qCvLpert[IR] *= (1.0+perturbScale);
        dCvpert     = qCvLpert[IR] - qCvL[IR];

        gas->Conservative2Primitive(qCvLpert, gmL, primL, temperaturesL);
        gas->Conservative2Primitive(qCvR, gmR, primR, temperaturesR);

        Cal_GMRES_Roe_Scheme_Flux(primL, primR, qCvLpert, qCvR, iFace, fluxp,
                                  face_proxy, invSchemePara);

        dFluxdCv[IR] = (fluxp[IR] - flux[IR][iFace])/dCvpert;
        dFluxdCv[IU] = (fluxp[IU] - flux[IU][iFace])/dCvpert;
        dFluxdCv[IV] = (fluxp[IV] - flux[IV][iFace])/dCvpert;
        dFluxdCv[IW] = (fluxp[IW] - flux[IW][iFace])/dCvpert;
        dFluxdCv[IP] = (fluxp[IP] - flux[IP][iFace])/dCvpert;

        dRdq[le * nEquation + IR][le * nEquation + IR] += dFluxdCv[IR] * area;
        dRdq[le * nEquation + IU][le * nEquation + IR] += dFluxdCv[IU] * area;
        dRdq[le * nEquation + IV][le * nEquation + IR] += dFluxdCv[IV] * area;
        dRdq[le * nEquation + IW][le * nEquation + IR] += dFluxdCv[IW] * area;
        dRdq[le * nEquation + IP][le * nEquation + IR] += dFluxdCv[IP] * area;

        dRdq[re * nEquation + IR][le * nEquation + IR] -= dFluxdCv[IR] * area;
        dRdq[re * nEquation + IU][le * nEquation + IR] -= dFluxdCv[IU] * area;
        dRdq[re * nEquation + IV][le * nEquation + IR] -= dFluxdCv[IV] * area;
        dRdq[re * nEquation + IW][le * nEquation + IR] -= dFluxdCv[IW] * area;
        dRdq[re * nEquation + IP][le * nEquation + IR] -= dFluxdCv[IP] * area;

        //! perturb the momentum rhou of the left cell of the face
        qCvLpert[IR] = qCvL[IR];
        qCvLpert[IU] = qCvL[IU];
        qCvLpert[IV] = qCvL[IV];
        qCvLpert[IW] = qCvL[IW];
        qCvLpert[IP] = qCvL[IP];
        qCvLpert[IU] += perturbScale;
        dCvpert = qCvLpert[IU] - qCvL[IU];

        gas->Conservative2Primitive(qCvLpert, gmL, primL, temperaturesL);
        gas->Conservative2Primitive(qCvR, gmR, primR, temperaturesR);

        Cal_GMRES_Roe_Scheme_Flux(primL, primR, qCvLpert, qCvR, iFace, fluxp,
                                  face_proxy, invSchemePara);

        dFluxdCv[IR] = (fluxp[IR] - flux[IR][iFace]) / dCvpert;
        dFluxdCv[IU] = (fluxp[IU] - flux[IU][iFace]) / dCvpert;
        dFluxdCv[IV] = (fluxp[IV] - flux[IV][iFace]) / dCvpert;
        dFluxdCv[IW] = (fluxp[IW] - flux[IW][iFace]) / dCvpert;
        dFluxdCv[IP] = (fluxp[IP] - flux[IP][iFace]) / dCvpert;

        dRdq[le * nEquation + IR][le * nEquation + IU] += dFluxdCv[IR] * area;
        dRdq[le * nEquation + IU][le * nEquation + IU] += dFluxdCv[IU] * area;
        dRdq[le * nEquation + IV][le * nEquation + IU] += dFluxdCv[IV] * area;
        dRdq[le * nEquation + IW][le * nEquation + IU] += dFluxdCv[IW] * area;
        dRdq[le * nEquation + IP][le * nEquation + IU] += dFluxdCv[IP] * area;

        dRdq[re * nEquation + IR][le * nEquation + IU] -= dFluxdCv[IR] * area;
        dRdq[re * nEquation + IU][le * nEquation + IU] -= dFluxdCv[IU] * area;
        dRdq[re * nEquation + IV][le * nEquation + IU] -= dFluxdCv[IV] * area;
        dRdq[re * nEquation + IW][le * nEquation + IU] -= dFluxdCv[IW] * area;
        dRdq[re * nEquation + IP][le * nEquation + IU] -= dFluxdCv[IP] * area;

        //! perturb the momentum rhov of the left cell of the face
        qCvLpert[IR] = qCvL[IR];
        qCvLpert[IU] = qCvL[IU];
        qCvLpert[IV] = qCvL[IV];
        qCvLpert[IW] = qCvL[IW];
        qCvLpert[IP] = qCvL[IP];
        qCvLpert[IV] += perturbScale;
        dCvpert = qCvLpert[IV] - qCvL[IV];

        gas->Conservative2Primitive(qCvLpert, gmL, primL, temperaturesL);
        gas->Conservative2Primitive(qCvR, gmR, primR, temperaturesR);

        Cal_GMRES_Roe_Scheme_Flux(primL, primR, qCvLpert, qCvR, iFace, fluxp,
            face_proxy, invSchemePara);

        dFluxdCv[IR] = (fluxp[IR] - flux[IR][iFace]) / dCvpert;
        dFluxdCv[IU] = (fluxp[IU] - flux[IU][iFace]) / dCvpert;
        dFluxdCv[IV] = (fluxp[IV] - flux[IV][iFace]) / dCvpert;
        dFluxdCv[IW] = (fluxp[IW] - flux[IW][iFace]) / dCvpert;
        dFluxdCv[IP] = (fluxp[IP] - flux[IP][iFace]) / dCvpert;

        dRdq[le * nEquation + IR][le * nEquation + IV] += dFluxdCv[IR] * area;
        dRdq[le * nEquation + IU][le * nEquation + IV] += dFluxdCv[IU] * area;
        dRdq[le * nEquation + IV][le * nEquation + IV] += dFluxdCv[IV] * area;
        dRdq[le * nEquation + IW][le * nEquation + IV] += dFluxdCv[IW] * area;
        dRdq[le * nEquation + IP][le * nEquation + IV] += dFluxdCv[IP] * area;

        dRdq[re * nEquation + IR][le * nEquation + IV] -= dFluxdCv[IR] * area;
        dRdq[re * nEquation + IU][le * nEquation + IV] -= dFluxdCv[IU] * area;
        dRdq[re * nEquation + IV][le * nEquation + IV] -= dFluxdCv[IV] * area;
        dRdq[re * nEquation + IW][le * nEquation + IV] -= dFluxdCv[IW] * area;
        dRdq[re * nEquation + IP][le * nEquation + IV] -= dFluxdCv[IP] * area;

        //! perturb the momentum rhow of the left cell of the face
        qCvLpert[IR] = qCvL[IR];
        qCvLpert[IU] = qCvL[IU];
        qCvLpert[IV] = qCvL[IV];
        qCvLpert[IW] = qCvL[IW];
        qCvLpert[IP] = qCvL[IP];
        qCvLpert[IW] += perturbScale;
        dCvpert = qCvLpert[IW] - qCvL[IW];

        gas->Conservative2Primitive(qCvLpert, gmL, primL, temperaturesL);
        gas->Conservative2Primitive(qCvR, gmR, primR, temperaturesR);

        Cal_GMRES_Roe_Scheme_Flux(primL, primR, qCvLpert, qCvR, iFace, fluxp,
            face_proxy, invSchemePara);

        dFluxdCv[IR] = (fluxp[IR] - flux[IR][iFace]) / dCvpert;
        dFluxdCv[IU] = (fluxp[IU] - flux[IU][iFace]) / dCvpert;
        dFluxdCv[IV] = (fluxp[IV] - flux[IV][iFace]) / dCvpert;
        dFluxdCv[IW] = (fluxp[IW] - flux[IW][iFace]) / dCvpert;
        dFluxdCv[IP] = (fluxp[IP] - flux[IP][iFace]) / dCvpert;

        dRdq[le * nEquation + IR][le * nEquation + IW] += dFluxdCv[IR] * area;
        dRdq[le * nEquation + IU][le * nEquation + IW] += dFluxdCv[IU] * area;
        dRdq[le * nEquation + IV][le * nEquation + IW] += dFluxdCv[IV] * area;
        dRdq[le * nEquation + IW][le * nEquation + IW] += dFluxdCv[IW] * area;
        dRdq[le * nEquation + IP][le * nEquation + IW] += dFluxdCv[IP] * area;

        dRdq[re * nEquation + IR][le * nEquation + IW] -= dFluxdCv[IR] * area;
        dRdq[re * nEquation + IU][le * nEquation + IW] -= dFluxdCv[IU] * area;
        dRdq[re * nEquation + IV][le * nEquation + IW] -= dFluxdCv[IV] * area;
        dRdq[re * nEquation + IW][le * nEquation + IW] -= dFluxdCv[IW] * area;
        dRdq[re * nEquation + IP][le * nEquation + IW] -= dFluxdCv[IP] * area;

        //! perturb the rhoE of the left cell of the face
        qCvLpert[IR] = qCvL[IR];
        qCvLpert[IU] = qCvL[IU];
        qCvLpert[IV] = qCvL[IV];
        qCvLpert[IW] = qCvL[IW];
        qCvLpert[IP] = qCvL[IP];
        qCvLpert[IP] *= (1.0+perturbScale);
        dCvpert = qCvLpert[IP] - qCvL[IP];

        gas->Conservative2Primitive(qCvLpert, gmL, primL, temperaturesL);
        gas->Conservative2Primitive(qCvR, gmR, primR, temperaturesR);

        Cal_GMRES_Roe_Scheme_Flux(primL, primR, qCvLpert, qCvR, iFace, fluxp,
            face_proxy, invSchemePara);

        dFluxdCv[IR] = (fluxp[IR] - flux[IR][iFace]) / dCvpert;
        dFluxdCv[IU] = (fluxp[IU] - flux[IU][iFace]) / dCvpert;
        dFluxdCv[IV] = (fluxp[IV] - flux[IV][iFace]) / dCvpert;
        dFluxdCv[IW] = (fluxp[IW] - flux[IW][iFace]) / dCvpert;
        dFluxdCv[IP] = (fluxp[IP] - flux[IP][iFace]) / dCvpert;

        dRdq[le * nEquation + IR][le * nEquation + IP] += dFluxdCv[IR] * area;
        dRdq[le * nEquation + IU][le * nEquation + IP] += dFluxdCv[IU] * area;
        dRdq[le * nEquation + IV][le * nEquation + IP] += dFluxdCv[IV] * area;
        dRdq[le * nEquation + IW][le * nEquation + IP] += dFluxdCv[IW] * area;
        dRdq[le * nEquation + IP][le * nEquation + IP] += dFluxdCv[IP] * area;

        dRdq[re * nEquation + IR][le * nEquation + IP] -= dFluxdCv[IR] * area;
        dRdq[re * nEquation + IU][le * nEquation + IP] -= dFluxdCv[IU] * area;
        dRdq[re * nEquation + IV][le * nEquation + IP] -= dFluxdCv[IV] * area;
        dRdq[re * nEquation + IW][le * nEquation + IP] -= dFluxdCv[IW] * area;
        dRdq[re * nEquation + IP][le * nEquation + IP] -= dFluxdCv[IP] * area;

        //! perturb the density of the right cell of the face
        qCvRpert[IR] = qCvR[IR];
        qCvRpert[IU] = qCvR[IU];
        qCvRpert[IV] = qCvR[IV];
        qCvRpert[IW] = qCvR[IW];
        qCvRpert[IP] = qCvR[IP];
        qCvRpert[IR] *= (1.0+perturbScale);
        dCvpert = qCvRpert[IR] - qCvR[IR];

        gas->Conservative2Primitive(qCvL, gmL, primL, temperaturesL);
        gas->Conservative2Primitive(qCvRpert, gmR, primR, temperaturesR);

        Cal_GMRES_Roe_Scheme_Flux(primL, primR, qCvL, qCvRpert, iFace, fluxp,
            face_proxy, invSchemePara);

        dFluxdCv[IR] = (fluxp[IR] - flux[IR][iFace]) / dCvpert;
        dFluxdCv[IU] = (fluxp[IU] - flux[IU][iFace]) / dCvpert;
        dFluxdCv[IV] = (fluxp[IV] - flux[IV][iFace]) / dCvpert;
        dFluxdCv[IW] = (fluxp[IW] - flux[IW][iFace]) / dCvpert;
        dFluxdCv[IP] = (fluxp[IP] - flux[IP][iFace]) / dCvpert;

        dRdq[le * nEquation + IR][re * nEquation + IR] += dFluxdCv[IR] * area;
        dRdq[le * nEquation + IU][re * nEquation + IR] += dFluxdCv[IU] * area;
        dRdq[le * nEquation + IV][re * nEquation + IR] += dFluxdCv[IV] * area;
        dRdq[le * nEquation + IW][re * nEquation + IR] += dFluxdCv[IW] * area;
        dRdq[le * nEquation + IP][re * nEquation + IR] += dFluxdCv[IP] * area;

        dRdq[re * nEquation + IR][re * nEquation + IR] -= dFluxdCv[IR] * area;
        dRdq[re * nEquation + IU][re * nEquation + IR] -= dFluxdCv[IU] * area;
        dRdq[re * nEquation + IV][re * nEquation + IR] -= dFluxdCv[IV] * area;
        dRdq[re * nEquation + IW][re * nEquation + IR] -= dFluxdCv[IW] * area;
        dRdq[re * nEquation + IP][re * nEquation + IR] -= dFluxdCv[IP] * area;

        //! perturb the momentum rhou of the right cell of the face
        qCvRpert[IR] = qCvR[IR];
        qCvRpert[IU] = qCvR[IU];
        qCvRpert[IV] = qCvR[IV];
        qCvRpert[IW] = qCvR[IW];
        qCvRpert[IP] = qCvR[IP];
        qCvRpert[IU] += perturbScale;
        dCvpert = qCvRpert[IU] - qCvR[IU];

        gas->Conservative2Primitive(qCvL, gmL, primL, temperaturesL);
        gas->Conservative2Primitive(qCvRpert, gmR, primR, temperaturesR);

        Cal_GMRES_Roe_Scheme_Flux(primL, primR, qCvL, qCvRpert, iFace, fluxp,
                                  face_proxy, invSchemePara);

        dFluxdCv[IR] = (fluxp[IR] - flux[IR][iFace]) / dCvpert;
        dFluxdCv[IU] = (fluxp[IU] - flux[IU][iFace]) / dCvpert;
        dFluxdCv[IV] = (fluxp[IV] - flux[IV][iFace]) / dCvpert;
        dFluxdCv[IW] = (fluxp[IW] - flux[IW][iFace]) / dCvpert;
        dFluxdCv[IP] = (fluxp[IP] - flux[IP][iFace]) / dCvpert;

        dRdq[le * nEquation + IR][re * nEquation + IU] += dFluxdCv[IR] * area;
        dRdq[le * nEquation + IU][re * nEquation + IU] += dFluxdCv[IU] * area;
        dRdq[le * nEquation + IV][re * nEquation + IU] += dFluxdCv[IV] * area;
        dRdq[le * nEquation + IW][re * nEquation + IU] += dFluxdCv[IW] * area;
        dRdq[le * nEquation + IP][re * nEquation + IU] += dFluxdCv[IP] * area;

        dRdq[re * nEquation + IR][re * nEquation + IU] -= dFluxdCv[IR] * area;
        dRdq[re * nEquation + IU][re * nEquation + IU] -= dFluxdCv[IU] * area;
        dRdq[re * nEquation + IV][re * nEquation + IU] -= dFluxdCv[IV] * area;
        dRdq[re * nEquation + IW][re * nEquation + IU] -= dFluxdCv[IW] * area;
        dRdq[re * nEquation + IP][re * nEquation + IU] -= dFluxdCv[IP] * area;

        //! perturb the momentum rhov of the right cell of the face
        qCvRpert[IR] = qCvR[IR];
        qCvRpert[IU] = qCvR[IU];
        qCvRpert[IV] = qCvR[IV];
        qCvRpert[IW] = qCvR[IW];
        qCvRpert[IP] = qCvR[IP];
        qCvRpert[IV] += perturbScale;
        dCvpert = qCvRpert[IV] - qCvR[IV];

        gas->Conservative2Primitive(qCvL, gmL, primL, temperaturesL);
        gas->Conservative2Primitive(qCvRpert, gmR, primR, temperaturesR);

        Cal_GMRES_Roe_Scheme_Flux(primL, primR, qCvL, qCvRpert, iFace, fluxp,
                                  face_proxy, invSchemePara);

        dFluxdCv[IR] = (fluxp[IR] - flux[IR][iFace]) / dCvpert;
        dFluxdCv[IU] = (fluxp[IU] - flux[IU][iFace]) / dCvpert;
        dFluxdCv[IV] = (fluxp[IV] - flux[IV][iFace]) / dCvpert;
        dFluxdCv[IW] = (fluxp[IW] - flux[IW][iFace]) / dCvpert;
        dFluxdCv[IP] = (fluxp[IP] - flux[IP][iFace]) / dCvpert;

        dRdq[le * nEquation + IR][re * nEquation + IV] += dFluxdCv[IR] * area;
        dRdq[le * nEquation + IU][re * nEquation + IV] += dFluxdCv[IU] * area;
        dRdq[le * nEquation + IV][re * nEquation + IV] += dFluxdCv[IV] * area;
        dRdq[le * nEquation + IW][re * nEquation + IV] += dFluxdCv[IW] * area;
        dRdq[le * nEquation + IP][re * nEquation + IV] += dFluxdCv[IP] * area;

        dRdq[re * nEquation + IR][re * nEquation + IV] -= dFluxdCv[IR] * area;
        dRdq[re * nEquation + IU][re * nEquation + IV] -= dFluxdCv[IU] * area;
        dRdq[re * nEquation + IV][re * nEquation + IV] -= dFluxdCv[IV] * area;
        dRdq[re * nEquation + IW][re * nEquation + IV] -= dFluxdCv[IW] * area;
        dRdq[re * nEquation + IP][re * nEquation + IV] -= dFluxdCv[IP] * area;

        //! perturb the momentum rhow of the right cell of the face
        qCvRpert[IR] = qCvR[IR];
        qCvRpert[IU] = qCvR[IU];
        qCvRpert[IV] = qCvR[IV];
        qCvRpert[IW] = qCvR[IW];
        qCvRpert[IP] = qCvR[IP];
        qCvRpert[IW] += perturbScale;
        dCvpert = qCvRpert[IW] - qCvR[IW];

        gas->Conservative2Primitive(qCvL, gmL, primL, temperaturesL);
        gas->Conservative2Primitive(qCvRpert, gmR, primR, temperaturesR);

        Cal_GMRES_Roe_Scheme_Flux(primL, primR, qCvL, qCvRpert, iFace, fluxp,
            face_proxy, invSchemePara);

        dFluxdCv[IR] = (fluxp[IR] - flux[IR][iFace]) / dCvpert;
        dFluxdCv[IU] = (fluxp[IU] - flux[IU][iFace]) / dCvpert;
        dFluxdCv[IV] = (fluxp[IV] - flux[IV][iFace]) / dCvpert;
        dFluxdCv[IW] = (fluxp[IW] - flux[IW][iFace]) / dCvpert;
        dFluxdCv[IP] = (fluxp[IP] - flux[IP][iFace]) / dCvpert;

        dRdq[le * nEquation + IR][re * nEquation + IW] += dFluxdCv[IR] * area;
        dRdq[le * nEquation + IU][re * nEquation + IW] += dFluxdCv[IU] * area;
        dRdq[le * nEquation + IV][re * nEquation + IW] += dFluxdCv[IV] * area;
        dRdq[le * nEquation + IW][re * nEquation + IW] += dFluxdCv[IW] * area;
        dRdq[le * nEquation + IP][re * nEquation + IW] += dFluxdCv[IP] * area;

        dRdq[re * nEquation + IR][re * nEquation + IW] -= dFluxdCv[IR] * area;
        dRdq[re * nEquation + IU][re * nEquation + IW] -= dFluxdCv[IU] * area;
        dRdq[re * nEquation + IV][re * nEquation + IW] -= dFluxdCv[IV] * area;
        dRdq[re * nEquation + IW][re * nEquation + IW] -= dFluxdCv[IW] * area;
        dRdq[re * nEquation + IP][re * nEquation + IW] -= dFluxdCv[IP] * area;

        //! perturb the momentum rhoE of the right cell of the face
        qCvRpert[IR] = qCvR[IR];
        qCvRpert[IU] = qCvR[IU];
        qCvRpert[IV] = qCvR[IV];
        qCvRpert[IW] = qCvR[IW];
        qCvRpert[IP] = qCvR[IP];
        qCvRpert[IP] *= (1.0+perturbScale);
        dCvpert = qCvRpert[IP] - qCvR[IP];

        gas->Conservative2Primitive(qCvL, gmL, primL, temperaturesL);
        gas->Conservative2Primitive(qCvRpert, gmR, primR, temperaturesR);

        Cal_GMRES_Roe_Scheme_Flux(primL, primR, qCvL, qCvRpert, iFace, fluxp,
            face_proxy, invSchemePara);

        dFluxdCv[IR] = (fluxp[IR] - flux[IR][iFace]) / dCvpert;
        dFluxdCv[IU] = (fluxp[IU] - flux[IU][iFace]) / dCvpert;
        dFluxdCv[IV] = (fluxp[IV] - flux[IV][iFace]) / dCvpert;
        dFluxdCv[IW] = (fluxp[IW] - flux[IW][iFace]) / dCvpert;
        dFluxdCv[IP] = (fluxp[IP] - flux[IP][iFace]) / dCvpert;

        dRdq[le * nEquation + IR][re * nEquation + IP] += dFluxdCv[IR] * area;
        dRdq[le * nEquation + IU][re * nEquation + IP] += dFluxdCv[IU] * area;
        dRdq[le * nEquation + IV][re * nEquation + IP] += dFluxdCv[IV] * area;
        dRdq[le * nEquation + IW][re * nEquation + IP] += dFluxdCv[IW] * area;
        dRdq[le * nEquation + IP][re * nEquation + IP] += dFluxdCv[IP] * area;

        dRdq[re * nEquation + IR][re * nEquation + IP] -= dFluxdCv[IR] * area;
        dRdq[re * nEquation + IU][re * nEquation + IP] -= dFluxdCv[IU] * area;
        dRdq[re * nEquation + IV][re * nEquation + IP] -= dFluxdCv[IV] * area;
        dRdq[re * nEquation + IW][re * nEquation + IP] -= dFluxdCv[IW] * area;
        dRdq[re * nEquation + IP][re * nEquation + IP] -= dFluxdCv[IP] * area;
     }

    delete [] prim;          prim = nullptr;
    delete [] primL;         primL = nullptr;
    delete [] primR;         primR = nullptr;
    delete [] qConserveL;    qConserveL = nullptr;
    delete [] qConserveR;    qConserveR = nullptr;
    delete [] qCvL;          qCvL = nullptr;
    delete [] dFluxdCv;      dFluxdCv = nullptr;
    delete [] qCvLpert;      qCvLpert = nullptr;
    delete [] qCvRpert;      qCvRpert = nullptr;
    delete [] qCvR;          qCvR = nullptr;
    delete [] fluxp;         fluxp = nullptr;
}

//! GMRES for finite difference
void GMRES_Roe_Scheme_ConservativeForm_FD(FaceProxy *face_proxy, InviscidSchemeParameter *invSchemePara)
{
    const RDouble LITTLE    = 1.0E-10;
    int nTemperatureModel   = invSchemePara->GetNumberOfTemperatureModel();
    int nEquation           = invSchemePara->GetNumberOfTotalEquation();
    int nLaminar            = invSchemePara->GetNumberOfLaminar();
    int nLength             = invSchemePara->GetLength();
    int nm                  = invSchemePara->Getnm();
    int RoeEntropyFixMethod = invSchemePara->GetEntropyFixMethod();
    int nChemical           = invSchemePara->GetNumberOfChemical();
    int nPrecondition       = invSchemePara->GetIfPrecondition();
    int nIdealState         = invSchemePara->GetIfIdealState();

    RDouble refMachNumber     = invSchemePara->GetMachNumber();
    RDouble entrFixCoefLeft   = invSchemePara->GetLeftEntropyFixCoefficient();
    RDouble entrFixCoefRight  = invSchemePara->GetRightEntropyFixCoefficient();
    RDouble **qL              = invSchemePara->GetLeftQ();
    RDouble **qR              = invSchemePara->GetRightQ();
    RDouble *gamaL            = invSchemePara->GetLeftGama();
    RDouble *gamaR            = invSchemePara->GetRightGama();
    RDouble *pressureCoeffL   = face_proxy->GetPressureCoefficientL();
    RDouble *pressureCoeffR   = face_proxy->GetPressureCoefficientR();
    RDouble *xfn              = invSchemePara->GetFaceNormalX();
    RDouble *yfn              = invSchemePara->GetFaceNormalY();
    RDouble *zfn              = invSchemePara->GetFaceNormalZ();
    RDouble *vgn              = invSchemePara->GetFaceVelocity();
    RDouble **flux            = invSchemePara->GetFlux();
    RDouble **tl              = invSchemePara->GetLeftTemperature();
    RDouble **tr              = invSchemePara->GetRightTemperature();
    int *leftcellindexofFace  = face_proxy->GetLeftCellIndexOfFace();
    int *rightcellindexofFace = face_proxy->GetRightCellIndexOfFace();
    GeomProxy *geomProxy      = face_proxy->GetGeomProxy();
    RDouble *areas            = geomProxy->GetFaceArea();
    RDouble **dRdq            = face_proxy->GetJacobianMatrix();
    RDouble mach2             = refMachNumber * refMachNumber;

    using namespace IDX;
    using namespace GAS_SPACE;

    RDouble *prim       = new RDouble[nEquation]();
    RDouble *qConserveL = new RDouble[nEquation]();
    RDouble *qConserveR = new RDouble[nEquation]();
    RDouble *qCvL       = new RDouble[nEquation]();
    RDouble *qCvR       = new RDouble[nEquation]();
    RDouble *qCvLpert   = new RDouble[nEquation]();
    RDouble *qCvRpert   = new RDouble[nEquation]();
    RDouble *dFluxdCv   = new RDouble[nEquation]();
    RDouble dCvpert;

    RDouble *primL = new RDouble[nEquation]();
    RDouble *primR = new RDouble[nEquation]();
    RDouble *fluxp = new RDouble[nEquation]();    //! flux proxy
    RDouble temperaturesL[3] = { 0 }, temperaturesR[3] = { 0 };
    RDouble perturbScale = 0.001;

    for (int iFace = 0; iFace < nLength; ++iFace)
    {
        int le       = leftcellindexofFace[iFace];
        int re       = rightcellindexofFace[iFace];
        RDouble area = areas[iFace];
        for (int iEqn = 0; iEqn < nEquation; ++iEqn)
        {
            primL[iEqn] = qL[iEqn][iFace];
            primR[iEqn] = qR[iEqn][iFace];
        }

        RDouble gmL = gamaL[iFace];
        RDouble gmR = gamaR[iFace];

        gas->Primitive2Conservative(primL, gmL, temperaturesL[ITV], temperaturesL[ITE], qCvL);
        gas->Primitive2Conservative(primR, gmR, temperaturesR[ITV], temperaturesR[ITE], qCvR);

        gas->Conservative2Primitive(qCvL, gmL, primL, temperaturesL);
        gas->Conservative2Primitive(qCvR, gmR, primR, temperaturesR);

        Cal_GMRES_Roe_Scheme_Flux(primL, primR, qCvL, qCvR,iFace, fluxp,
                                  face_proxy,invSchemePara);

        flux[IR][iFace] = fluxp[IR];
        flux[IU][iFace] = fluxp[IU];
        flux[IV][iFace] = fluxp[IV];
        flux[IW][iFace] = fluxp[IW];
        flux[IP][iFace] = fluxp[IP];

        //! perturb the density of the left cell of the face
        qCvLpert[IR] = qCvL[IR];
        qCvLpert[IU] = qCvL[IU];
        qCvLpert[IV] = qCvL[IV];
        qCvLpert[IW] = qCvL[IW];
        qCvLpert[IP] = qCvL[IP];
        qCvLpert[IR] *= (1.0+perturbScale);
        dCvpert     = qCvLpert[IR] - qCvL[IR];

        gas->Conservative2Primitive(qCvLpert, gmL, primL, temperaturesL);
        gas->Conservative2Primitive(qCvR, gmR, primR, temperaturesR);

        Cal_GMRES_Roe_Scheme_Flux(primL, primR, qCvLpert, qCvR, iFace, fluxp,
                                  face_proxy, invSchemePara);

        dFluxdCv[IR] = (fluxp[IR] - flux[IR][iFace])/dCvpert;
        dFluxdCv[IU] = (fluxp[IU] - flux[IU][iFace])/dCvpert;
        dFluxdCv[IV] = (fluxp[IV] - flux[IV][iFace])/dCvpert;
        dFluxdCv[IW] = (fluxp[IW] - flux[IW][iFace])/dCvpert;
        dFluxdCv[IP] = (fluxp[IP] - flux[IP][iFace])/dCvpert;

        dRdq[le * nEquation + IR][le * nEquation + IR] += dFluxdCv[IR] * area;
        dRdq[le * nEquation + IU][le * nEquation + IR] += dFluxdCv[IU] * area;
        dRdq[le * nEquation + IV][le * nEquation + IR] += dFluxdCv[IV] * area;
        dRdq[le * nEquation + IW][le * nEquation + IR] += dFluxdCv[IW] * area;
        dRdq[le * nEquation + IP][le * nEquation + IR] += dFluxdCv[IP] * area;

        dRdq[re * nEquation + IR][le * nEquation + IR] -= dFluxdCv[IR] * area;
        dRdq[re * nEquation + IU][le * nEquation + IR] -= dFluxdCv[IU] * area;
        dRdq[re * nEquation + IV][le * nEquation + IR] -= dFluxdCv[IV] * area;
        dRdq[re * nEquation + IW][le * nEquation + IR] -= dFluxdCv[IW] * area;
        dRdq[re * nEquation + IP][le * nEquation + IR] -= dFluxdCv[IP] * area;

        //! perturb the momentum rhou of the left cell of the face
        qCvLpert[IR] = qCvL[IR];
        qCvLpert[IU] = qCvL[IU];
        qCvLpert[IV] = qCvL[IV];
        qCvLpert[IW] = qCvL[IW];
        qCvLpert[IP] = qCvL[IP];
        qCvLpert[IU] += perturbScale;
        dCvpert = qCvLpert[IU] - qCvL[IU];

        gas->Conservative2Primitive(qCvLpert, gmL, primL, temperaturesL);
        gas->Conservative2Primitive(qCvR, gmR, primR, temperaturesR);

        Cal_GMRES_Roe_Scheme_Flux(primL, primR, qCvLpert, qCvR, iFace, fluxp,
                                  face_proxy, invSchemePara);

        dFluxdCv[IR] = (fluxp[IR] - flux[IR][iFace]) / dCvpert;
        dFluxdCv[IU] = (fluxp[IU] - flux[IU][iFace]) / dCvpert;
        dFluxdCv[IV] = (fluxp[IV] - flux[IV][iFace]) / dCvpert;
        dFluxdCv[IW] = (fluxp[IW] - flux[IW][iFace]) / dCvpert;
        dFluxdCv[IP] = (fluxp[IP] - flux[IP][iFace]) / dCvpert;

        dRdq[le * nEquation + IR][le * nEquation + IU] += dFluxdCv[IR] * area;
        dRdq[le * nEquation + IU][le * nEquation + IU] += dFluxdCv[IU] * area;
        dRdq[le * nEquation + IV][le * nEquation + IU] += dFluxdCv[IV] * area;
        dRdq[le * nEquation + IW][le * nEquation + IU] += dFluxdCv[IW] * area;
        dRdq[le * nEquation + IP][le * nEquation + IU] += dFluxdCv[IP] * area;

        dRdq[re * nEquation + IR][le * nEquation + IU] -= dFluxdCv[IR] * area;
        dRdq[re * nEquation + IU][le * nEquation + IU] -= dFluxdCv[IU] * area;
        dRdq[re * nEquation + IV][le * nEquation + IU] -= dFluxdCv[IV] * area;
        dRdq[re * nEquation + IW][le * nEquation + IU] -= dFluxdCv[IW] * area;
        dRdq[re * nEquation + IP][le * nEquation + IU] -= dFluxdCv[IP] * area;

        //! perturb the momentum rhov of the left cell of the face
        qCvLpert[IR] = qCvL[IR];
        qCvLpert[IU] = qCvL[IU];
        qCvLpert[IV] = qCvL[IV];
        qCvLpert[IW] = qCvL[IW];
        qCvLpert[IP] = qCvL[IP];
        qCvLpert[IV] += perturbScale;
        dCvpert = qCvLpert[IV] - qCvL[IV];

        gas->Conservative2Primitive(qCvLpert, gmL, primL, temperaturesL);
        gas->Conservative2Primitive(qCvR, gmR, primR, temperaturesR);

        Cal_GMRES_Roe_Scheme_Flux(primL, primR, qCvLpert, qCvR, iFace, fluxp,
            face_proxy, invSchemePara);

        dFluxdCv[IR] = (fluxp[IR] - flux[IR][iFace]) / dCvpert;
        dFluxdCv[IU] = (fluxp[IU] - flux[IU][iFace]) / dCvpert;
        dFluxdCv[IV] = (fluxp[IV] - flux[IV][iFace]) / dCvpert;
        dFluxdCv[IW] = (fluxp[IW] - flux[IW][iFace]) / dCvpert;
        dFluxdCv[IP] = (fluxp[IP] - flux[IP][iFace]) / dCvpert;

        dRdq[le * nEquation + IR][le * nEquation + IV] += dFluxdCv[IR] * area;
        dRdq[le * nEquation + IU][le * nEquation + IV] += dFluxdCv[IU] * area;
        dRdq[le * nEquation + IV][le * nEquation + IV] += dFluxdCv[IV] * area;
        dRdq[le * nEquation + IW][le * nEquation + IV] += dFluxdCv[IW] * area;
        dRdq[le * nEquation + IP][le * nEquation + IV] += dFluxdCv[IP] * area;

        dRdq[re * nEquation + IR][le * nEquation + IV] -= dFluxdCv[IR] * area;
        dRdq[re * nEquation + IU][le * nEquation + IV] -= dFluxdCv[IU] * area;
        dRdq[re * nEquation + IV][le * nEquation + IV] -= dFluxdCv[IV] * area;
        dRdq[re * nEquation + IW][le * nEquation + IV] -= dFluxdCv[IW] * area;
        dRdq[re * nEquation + IP][le * nEquation + IV] -= dFluxdCv[IP] * area;

        //! perturb the momentum rhow of the left cell of the face
        qCvLpert[IR] = qCvL[IR];
        qCvLpert[IU] = qCvL[IU];
        qCvLpert[IV] = qCvL[IV];
        qCvLpert[IW] = qCvL[IW];
        qCvLpert[IP] = qCvL[IP];
        qCvLpert[IW] += perturbScale;
        dCvpert = qCvLpert[IW] - qCvL[IW];

        gas->Conservative2Primitive(qCvLpert, gmL, primL, temperaturesL);
        gas->Conservative2Primitive(qCvR, gmR, primR, temperaturesR);

        Cal_GMRES_Roe_Scheme_Flux(primL, primR, qCvLpert, qCvR, iFace, fluxp,
            face_proxy, invSchemePara);

        dFluxdCv[IR] = (fluxp[IR] - flux[IR][iFace]) / dCvpert;
        dFluxdCv[IU] = (fluxp[IU] - flux[IU][iFace]) / dCvpert;
        dFluxdCv[IV] = (fluxp[IV] - flux[IV][iFace]) / dCvpert;
        dFluxdCv[IW] = (fluxp[IW] - flux[IW][iFace]) / dCvpert;
        dFluxdCv[IP] = (fluxp[IP] - flux[IP][iFace]) / dCvpert;

        dRdq[le * nEquation + IR][le * nEquation + IW] += dFluxdCv[IR] * area;
        dRdq[le * nEquation + IU][le * nEquation + IW] += dFluxdCv[IU] * area;
        dRdq[le * nEquation + IV][le * nEquation + IW] += dFluxdCv[IV] * area;
        dRdq[le * nEquation + IW][le * nEquation + IW] += dFluxdCv[IW] * area;
        dRdq[le * nEquation + IP][le * nEquation + IW] += dFluxdCv[IP] * area;

        dRdq[re * nEquation + IR][le * nEquation + IW] -= dFluxdCv[IR] * area;
        dRdq[re * nEquation + IU][le * nEquation + IW] -= dFluxdCv[IU] * area;
        dRdq[re * nEquation + IV][le * nEquation + IW] -= dFluxdCv[IV] * area;
        dRdq[re * nEquation + IW][le * nEquation + IW] -= dFluxdCv[IW] * area;
        dRdq[re * nEquation + IP][le * nEquation + IW] -= dFluxdCv[IP] * area;

        //! perturb the rhoE of the left cell of the face
        qCvLpert[IR] = qCvL[IR];
        qCvLpert[IU] = qCvL[IU];
        qCvLpert[IV] = qCvL[IV];
        qCvLpert[IW] = qCvL[IW];
        qCvLpert[IP] = qCvL[IP];
        qCvLpert[IP] *= (1.0+perturbScale);
        dCvpert = qCvLpert[IP] - qCvL[IP];

        gas->Conservative2Primitive(qCvLpert, gmL, primL, temperaturesL);
        gas->Conservative2Primitive(qCvR, gmR, primR, temperaturesR);

        Cal_GMRES_Roe_Scheme_Flux(primL, primR, qCvLpert, qCvR, iFace, fluxp,
            face_proxy, invSchemePara);

        dFluxdCv[IR] = (fluxp[IR] - flux[IR][iFace]) / dCvpert;
        dFluxdCv[IU] = (fluxp[IU] - flux[IU][iFace]) / dCvpert;
        dFluxdCv[IV] = (fluxp[IV] - flux[IV][iFace]) / dCvpert;
        dFluxdCv[IW] = (fluxp[IW] - flux[IW][iFace]) / dCvpert;
        dFluxdCv[IP] = (fluxp[IP] - flux[IP][iFace]) / dCvpert;

        dRdq[le * nEquation + IR][le * nEquation + IP] += dFluxdCv[IR] * area;
        dRdq[le * nEquation + IU][le * nEquation + IP] += dFluxdCv[IU] * area;
        dRdq[le * nEquation + IV][le * nEquation + IP] += dFluxdCv[IV] * area;
        dRdq[le * nEquation + IW][le * nEquation + IP] += dFluxdCv[IW] * area;
        dRdq[le * nEquation + IP][le * nEquation + IP] += dFluxdCv[IP] * area;

        dRdq[re * nEquation + IR][le * nEquation + IP] -= dFluxdCv[IR] * area;
        dRdq[re * nEquation + IU][le * nEquation + IP] -= dFluxdCv[IU] * area;
        dRdq[re * nEquation + IV][le * nEquation + IP] -= dFluxdCv[IV] * area;
        dRdq[re * nEquation + IW][le * nEquation + IP] -= dFluxdCv[IW] * area;
        dRdq[re * nEquation + IP][le * nEquation + IP] -= dFluxdCv[IP] * area;

        //! perturb the density of the right cell of the face
        qCvRpert[IR] = qCvR[IR];
        qCvRpert[IU] = qCvR[IU];
        qCvRpert[IV] = qCvR[IV];
        qCvRpert[IW] = qCvR[IW];
        qCvRpert[IP] = qCvR[IP];
        qCvRpert[IR] *= (1.0+perturbScale);
        dCvpert = qCvRpert[IR] - qCvR[IR];

        gas->Conservative2Primitive(qCvL, gmL, primL, temperaturesL);
        gas->Conservative2Primitive(qCvRpert, gmR, primR, temperaturesR);

        Cal_GMRES_Roe_Scheme_Flux(primL, primR, qCvL, qCvRpert, iFace, fluxp,
            face_proxy, invSchemePara);

        dFluxdCv[IR] = (fluxp[IR] - flux[IR][iFace]) / dCvpert;
        dFluxdCv[IU] = (fluxp[IU] - flux[IU][iFace]) / dCvpert;
        dFluxdCv[IV] = (fluxp[IV] - flux[IV][iFace]) / dCvpert;
        dFluxdCv[IW] = (fluxp[IW] - flux[IW][iFace]) / dCvpert;
        dFluxdCv[IP] = (fluxp[IP] - flux[IP][iFace]) / dCvpert;

        dRdq[le * nEquation + IR][re * nEquation + IR] += dFluxdCv[IR] * area;
        dRdq[le * nEquation + IU][re * nEquation + IR] += dFluxdCv[IU] * area;
        dRdq[le * nEquation + IV][re * nEquation + IR] += dFluxdCv[IV] * area;
        dRdq[le * nEquation + IW][re * nEquation + IR] += dFluxdCv[IW] * area;
        dRdq[le * nEquation + IP][re * nEquation + IR] += dFluxdCv[IP] * area;

        dRdq[re * nEquation + IR][re * nEquation + IR] -= dFluxdCv[IR] * area;
        dRdq[re * nEquation + IU][re * nEquation + IR] -= dFluxdCv[IU] * area;
        dRdq[re * nEquation + IV][re * nEquation + IR] -= dFluxdCv[IV] * area;
        dRdq[re * nEquation + IW][re * nEquation + IR] -= dFluxdCv[IW] * area;
        dRdq[re * nEquation + IP][re * nEquation + IR] -= dFluxdCv[IP] * area;

        //! perturb the momentum rhou of the right cell of the face
        qCvRpert[IR] = qCvR[IR];
        qCvRpert[IU] = qCvR[IU];
        qCvRpert[IV] = qCvR[IV];
        qCvRpert[IW] = qCvR[IW];
        qCvRpert[IP] = qCvR[IP];
        qCvRpert[IU] += perturbScale;
        dCvpert = qCvRpert[IU] - qCvR[IU];

        gas->Conservative2Primitive(qCvL, gmL, primL, temperaturesL);
        gas->Conservative2Primitive(qCvRpert, gmR, primR, temperaturesR);

        Cal_GMRES_Roe_Scheme_Flux(primL, primR, qCvL, qCvRpert, iFace, fluxp,
                                  face_proxy, invSchemePara);

        dFluxdCv[IR] = (fluxp[IR] - flux[IR][iFace]) / dCvpert;
        dFluxdCv[IU] = (fluxp[IU] - flux[IU][iFace]) / dCvpert;
        dFluxdCv[IV] = (fluxp[IV] - flux[IV][iFace]) / dCvpert;
        dFluxdCv[IW] = (fluxp[IW] - flux[IW][iFace]) / dCvpert;
        dFluxdCv[IP] = (fluxp[IP] - flux[IP][iFace]) / dCvpert;

        dRdq[le * nEquation + IR][re * nEquation + IU] += dFluxdCv[IR] * area;
        dRdq[le * nEquation + IU][re * nEquation + IU] += dFluxdCv[IU] * area;
        dRdq[le * nEquation + IV][re * nEquation + IU] += dFluxdCv[IV] * area;
        dRdq[le * nEquation + IW][re * nEquation + IU] += dFluxdCv[IW] * area;
        dRdq[le * nEquation + IP][re * nEquation + IU] += dFluxdCv[IP] * area;

        dRdq[re * nEquation + IR][re * nEquation + IU] -= dFluxdCv[IR] * area;
        dRdq[re * nEquation + IU][re * nEquation + IU] -= dFluxdCv[IU] * area;
        dRdq[re * nEquation + IV][re * nEquation + IU] -= dFluxdCv[IV] * area;
        dRdq[re * nEquation + IW][re * nEquation + IU] -= dFluxdCv[IW] * area;
        dRdq[re * nEquation + IP][re * nEquation + IU] -= dFluxdCv[IP] * area;

        //! perturb the momentum rhov of the right cell of the face
        qCvRpert[IR] = qCvR[IR];
        qCvRpert[IU] = qCvR[IU];
        qCvRpert[IV] = qCvR[IV];
        qCvRpert[IW] = qCvR[IW];
        qCvRpert[IP] = qCvR[IP];
        qCvRpert[IV] += perturbScale;
        dCvpert = qCvRpert[IV] - qCvR[IV];

        gas->Conservative2Primitive(qCvL, gmL, primL, temperaturesL);
        gas->Conservative2Primitive(qCvRpert, gmR, primR, temperaturesR);

        Cal_GMRES_Roe_Scheme_Flux(primL, primR, qCvL, qCvRpert, iFace, fluxp,
                                  face_proxy, invSchemePara);

        dFluxdCv[IR] = (fluxp[IR] - flux[IR][iFace]) / dCvpert;
        dFluxdCv[IU] = (fluxp[IU] - flux[IU][iFace]) / dCvpert;
        dFluxdCv[IV] = (fluxp[IV] - flux[IV][iFace]) / dCvpert;
        dFluxdCv[IW] = (fluxp[IW] - flux[IW][iFace]) / dCvpert;
        dFluxdCv[IP] = (fluxp[IP] - flux[IP][iFace]) / dCvpert;

        dRdq[le * nEquation + IR][re * nEquation + IV] += dFluxdCv[IR] * area;
        dRdq[le * nEquation + IU][re * nEquation + IV] += dFluxdCv[IU] * area;
        dRdq[le * nEquation + IV][re * nEquation + IV] += dFluxdCv[IV] * area;
        dRdq[le * nEquation + IW][re * nEquation + IV] += dFluxdCv[IW] * area;
        dRdq[le * nEquation + IP][re * nEquation + IV] += dFluxdCv[IP] * area;

        dRdq[re * nEquation + IR][re * nEquation + IV] -= dFluxdCv[IR] * area;
        dRdq[re * nEquation + IU][re * nEquation + IV] -= dFluxdCv[IU] * area;
        dRdq[re * nEquation + IV][re * nEquation + IV] -= dFluxdCv[IV] * area;
        dRdq[re * nEquation + IW][re * nEquation + IV] -= dFluxdCv[IW] * area;
        dRdq[re * nEquation + IP][re * nEquation + IV] -= dFluxdCv[IP] * area;

        //! perturb the momentum rhow of the right cell of the face
        qCvRpert[IR] = qCvR[IR];
        qCvRpert[IU] = qCvR[IU];
        qCvRpert[IV] = qCvR[IV];
        qCvRpert[IW] = qCvR[IW];
        qCvRpert[IP] = qCvR[IP];
        qCvRpert[IW] += perturbScale;
        dCvpert = qCvRpert[IW] - qCvR[IW];

        gas->Conservative2Primitive(qCvL, gmL, primL, temperaturesL);
        gas->Conservative2Primitive(qCvRpert, gmR, primR, temperaturesR);

        Cal_GMRES_Roe_Scheme_Flux(primL, primR, qCvL, qCvRpert, iFace, fluxp,
            face_proxy, invSchemePara);

        dFluxdCv[IR] = (fluxp[IR] - flux[IR][iFace]) / dCvpert;
        dFluxdCv[IU] = (fluxp[IU] - flux[IU][iFace]) / dCvpert;
        dFluxdCv[IV] = (fluxp[IV] - flux[IV][iFace]) / dCvpert;
        dFluxdCv[IW] = (fluxp[IW] - flux[IW][iFace]) / dCvpert;
        dFluxdCv[IP] = (fluxp[IP] - flux[IP][iFace]) / dCvpert;

        dRdq[le * nEquation + IR][re * nEquation + IW] += dFluxdCv[IR] * area;
        dRdq[le * nEquation + IU][re * nEquation + IW] += dFluxdCv[IU] * area;
        dRdq[le * nEquation + IV][re * nEquation + IW] += dFluxdCv[IV] * area;
        dRdq[le * nEquation + IW][re * nEquation + IW] += dFluxdCv[IW] * area;
        dRdq[le * nEquation + IP][re * nEquation + IW] += dFluxdCv[IP] * area;

        dRdq[re * nEquation + IR][re * nEquation + IW] -= dFluxdCv[IR] * area;
        dRdq[re * nEquation + IU][re * nEquation + IW] -= dFluxdCv[IU] * area;
        dRdq[re * nEquation + IV][re * nEquation + IW] -= dFluxdCv[IV] * area;
        dRdq[re * nEquation + IW][re * nEquation + IW] -= dFluxdCv[IW] * area;
        dRdq[re * nEquation + IP][re * nEquation + IW] -= dFluxdCv[IP] * area;

        //! perturb the momentum rhoE of the right cell of the face
        qCvRpert[IR] = qCvR[IR];
        qCvRpert[IU] = qCvR[IU];
        qCvRpert[IV] = qCvR[IV];
        qCvRpert[IW] = qCvR[IW];
        qCvRpert[IP] = qCvR[IP];
        qCvRpert[IP] *= (1.0+perturbScale);
        dCvpert = qCvRpert[IP] - qCvR[IP];

        gas->Conservative2Primitive(qCvL, gmL, primL, temperaturesL);
        gas->Conservative2Primitive(qCvRpert, gmR, primR, temperaturesR);

        Cal_GMRES_Roe_Scheme_Flux(primL, primR, qCvL, qCvRpert, iFace, fluxp,
            face_proxy, invSchemePara);

        dFluxdCv[IR] = (fluxp[IR] - flux[IR][iFace]) / dCvpert;
        dFluxdCv[IU] = (fluxp[IU] - flux[IU][iFace]) / dCvpert;
        dFluxdCv[IV] = (fluxp[IV] - flux[IV][iFace]) / dCvpert;
        dFluxdCv[IW] = (fluxp[IW] - flux[IW][iFace]) / dCvpert;
        dFluxdCv[IP] = (fluxp[IP] - flux[IP][iFace]) / dCvpert;

        dRdq[le * nEquation + IR][re * nEquation + IP] += dFluxdCv[IR] * area;
        dRdq[le * nEquation + IU][re * nEquation + IP] += dFluxdCv[IU] * area;
        dRdq[le * nEquation + IV][re * nEquation + IP] += dFluxdCv[IV] * area;
        dRdq[le * nEquation + IW][re * nEquation + IP] += dFluxdCv[IW] * area;
        dRdq[le * nEquation + IP][re * nEquation + IP] += dFluxdCv[IP] * area;

        dRdq[re * nEquation + IR][re * nEquation + IP] -= dFluxdCv[IR] * area;
        dRdq[re * nEquation + IU][re * nEquation + IP] -= dFluxdCv[IU] * area;
        dRdq[re * nEquation + IV][re * nEquation + IP] -= dFluxdCv[IV] * area;
        dRdq[re * nEquation + IW][re * nEquation + IP] -= dFluxdCv[IW] * area;
        dRdq[re * nEquation + IP][re * nEquation + IP] -= dFluxdCv[IP] * area;
}

    delete [] prim;          prim = nullptr;
    delete [] primL;         primL = nullptr;
    delete [] primR;         primR = nullptr;
    delete [] qConserveL;    qConserveL = nullptr;
    delete [] qConserveR;    qConserveR = nullptr;
    delete [] qCvL;          qCvL = nullptr;
    delete [] dFluxdCv;      dFluxdCv = nullptr;
    delete [] qCvLpert;      qCvLpert = nullptr;
    delete [] qCvRpert;      qCvRpert = nullptr;
    delete [] qCvR;          qCvR = nullptr;
    delete [] fluxp;         fluxp = nullptr;
}

//! GMRES for rhow
void GMRES_Roe_Scheme_ConservativeForm_rhow(FaceProxy *face_proxy, InviscidSchemeParameter *invSchemePara)
{
    const RDouble LITTLE    = 1.0E-10;
    int nTemperatureModel   = invSchemePara->GetNumberOfTemperatureModel();
    int nEquation           = invSchemePara->GetNumberOfTotalEquation();
    int nLaminar            = invSchemePara->GetNumberOfLaminar();
    int nLength             = invSchemePara->GetLength();
    int nm                  = invSchemePara->Getnm();
    int RoeEntropyFixMethod = invSchemePara->GetEntropyFixMethod();
    int nChemical           = invSchemePara->GetNumberOfChemical();
    int nPrecondition       = invSchemePara->GetIfPrecondition();
    int nIdealState         = invSchemePara->GetIfIdealState();

    RDouble refMachNumber     = invSchemePara->GetMachNumber();
    RDouble entrFixCoefLeft   = invSchemePara->GetLeftEntropyFixCoefficient();
    RDouble entrFixCoefRight  = invSchemePara->GetRightEntropyFixCoefficient();
    RDouble **qL              = invSchemePara->GetLeftQ();
    RDouble **qR              = invSchemePara->GetRightQ();
    RDouble *gamaL            = invSchemePara->GetLeftGama();
    RDouble *gamaR            = invSchemePara->GetRightGama();
    RDouble *pressureCoeffL   = face_proxy->GetPressureCoefficientL();
    RDouble *pressureCoeffR   = face_proxy->GetPressureCoefficientR();
    RDouble *xfn              = invSchemePara->GetFaceNormalX();
    RDouble *yfn              = invSchemePara->GetFaceNormalY();
    RDouble *zfn              = invSchemePara->GetFaceNormalZ();
    RDouble *vgn              = invSchemePara->GetFaceVelocity();
    RDouble **flux            = invSchemePara->GetFlux();
    RDouble **tl              = invSchemePara->GetLeftTemperature();
    RDouble **tr              = invSchemePara->GetRightTemperature();
    int *leftcellindexofFace  = face_proxy->GetLeftCellIndexOfFace();
    int *rightcellindexofFace = face_proxy->GetRightCellIndexOfFace();
    GeomProxy *geomProxy      = face_proxy->GetGeomProxy();
    RDouble *areas            = geomProxy->GetFaceArea();
    RDouble **dRdq            = face_proxy->GetJacobianMatrix();
    RDouble mach2             = refMachNumber * refMachNumber;

    using namespace IDX;
    using namespace GAS_SPACE;

    RDouble *prim       = new RDouble[nEquation]();
    RDouble *qConserveL = new RDouble[nEquation]();
    RDouble *qConserveR = new RDouble[nEquation]();
    RDouble *qCvL       = new RDouble[nEquation]();
    RDouble *qCvR       = new RDouble[nEquation]();
    RDouble *qCvLpert   = new RDouble[nEquation]();
    RDouble *qCvRpert   = new RDouble[nEquation]();
    RDouble *dFluxdCv   = new RDouble[nEquation]();
    RDouble dCvpert;

    RDouble *primL = new RDouble[nEquation]();
    RDouble *primR = new RDouble[nEquation]();
    RDouble *fluxp = new RDouble[nEquation]();    //! flux proxy
    RDouble temperaturesL[3] = { 0 }, temperaturesR[3] = { 0 };
    RDouble perturbScale = 0.001;

    for (int iFace = 0; iFace < nLength; ++iFace)
    {
        int le          = leftcellindexofFace[iFace];
        int re          = rightcellindexofFace[iFace];
        RDouble area    = areas[iFace];
        for (int iEqn = 0; iEqn < nEquation; ++iEqn)
        {
            primL[iEqn] = qL[iEqn][iFace];
            primR[iEqn] = qR[iEqn][iFace];
        }

        RDouble gmL = gamaL[iFace];
        RDouble gmR = gamaR[iFace];

        gas->Primitive2Conservative(primL, gmL, temperaturesL[ITV], temperaturesL[ITE], qCvL);
        gas->Primitive2Conservative(primR, gmR, temperaturesR[ITV], temperaturesR[ITE], qCvR);

        gas->Conservative2Primitive(qCvL, gmL, primL, temperaturesL);
        gas->Conservative2Primitive(qCvR, gmR, primR, temperaturesR);

        Cal_GMRES_Roe_Scheme_Flux(primL, primR, qCvL, qCvR,iFace, fluxp,
                                  face_proxy,invSchemePara);

        flux[IR][iFace] = fluxp[IR];
        flux[IU][iFace] = fluxp[IU];
        flux[IV][iFace] = fluxp[IV];
        flux[IW][iFace] = fluxp[IW];
        flux[IP][iFace] = fluxp[IP];

        //! perturb the density of the left cell of the face
        qCvLpert[IR] = qCvL[IR];
        qCvLpert[IU] = qCvL[IU];
        qCvLpert[IV] = qCvL[IV];
        qCvLpert[IW] = qCvL[IW];
        qCvLpert[IP] = qCvL[IP];
        qCvLpert[IR] *= (1.0+perturbScale);
        dCvpert     = qCvLpert[IR] - qCvL[IR];

        gas->Conservative2Primitive(qCvLpert, gmL, primL, temperaturesL);
        gas->Conservative2Primitive(qCvR, gmR, primR, temperaturesR);

        Cal_GMRES_Roe_Scheme_Flux(primL, primR, qCvLpert, qCvR, iFace, fluxp,
                                  face_proxy, invSchemePara);

        dFluxdCv[IR] = (fluxp[IR] - flux[IR][iFace])/dCvpert;
        dFluxdCv[IU] = (fluxp[IU] - flux[IU][iFace])/dCvpert;
        dFluxdCv[IV] = (fluxp[IV] - flux[IV][iFace])/dCvpert;
        dFluxdCv[IW] = (fluxp[IW] - flux[IW][iFace])/dCvpert;
        dFluxdCv[IP] = (fluxp[IP] - flux[IP][iFace])/dCvpert;

        dRdq[le * nEquation + IR][le * nEquation + IR] += dFluxdCv[IR] * area;
        dRdq[le * nEquation + IU][le * nEquation + IR] += dFluxdCv[IU] * area;
        dRdq[le * nEquation + IV][le * nEquation + IR] += dFluxdCv[IV] * area;
        dRdq[le * nEquation + IW][le * nEquation + IR] = 0;
        dRdq[le * nEquation + IP][le * nEquation + IR] += dFluxdCv[IP] * area;

        dRdq[re * nEquation + IR][le * nEquation + IR] -= dFluxdCv[IR] * area;
        dRdq[re * nEquation + IU][le * nEquation + IR] -= dFluxdCv[IU] * area;
        dRdq[re * nEquation + IV][le * nEquation + IR] -= dFluxdCv[IV] * area;
        dRdq[re * nEquation + IW][le * nEquation + IR] = 0;
        dRdq[re * nEquation + IP][le * nEquation + IR] -= dFluxdCv[IP] * area;

        //! perturb the momentum rhou of the left cell of the face
        qCvLpert[IR] = qCvL[IR];
        qCvLpert[IU] = qCvL[IU];
        qCvLpert[IV] = qCvL[IV];
        qCvLpert[IW] = qCvL[IW];
        qCvLpert[IP] = qCvL[IP];
        qCvLpert[IU] += perturbScale;
        dCvpert = qCvLpert[IU] - qCvL[IU];

        gas->Conservative2Primitive(qCvLpert, gmL, primL, temperaturesL);
        gas->Conservative2Primitive(qCvR, gmR, primR, temperaturesR);

        Cal_GMRES_Roe_Scheme_Flux(primL, primR, qCvLpert, qCvR, iFace, fluxp,
                                  face_proxy, invSchemePara);

        dFluxdCv[IR] = (fluxp[IR] - flux[IR][iFace]) / dCvpert;
        dFluxdCv[IU] = (fluxp[IU] - flux[IU][iFace]) / dCvpert;
        dFluxdCv[IV] = (fluxp[IV] - flux[IV][iFace]) / dCvpert;
        dFluxdCv[IW] = (fluxp[IW] - flux[IW][iFace]) / dCvpert;
        dFluxdCv[IP] = (fluxp[IP] - flux[IP][iFace]) / dCvpert;

        dRdq[le * nEquation + IR][le * nEquation + IU] += dFluxdCv[IR] * area;
        dRdq[le * nEquation + IU][le * nEquation + IU] += dFluxdCv[IU] * area;
        dRdq[le * nEquation + IV][le * nEquation + IU] += dFluxdCv[IV] * area;
        dRdq[le * nEquation + IW][le * nEquation + IU] = 0;
        dRdq[le * nEquation + IP][le * nEquation + IU] += dFluxdCv[IP] * area;

        dRdq[re * nEquation + IR][le * nEquation + IU] -= dFluxdCv[IR] * area;
        dRdq[re * nEquation + IU][le * nEquation + IU] -= dFluxdCv[IU] * area;
        dRdq[re * nEquation + IV][le * nEquation + IU] -= dFluxdCv[IV] * area;
        dRdq[re * nEquation + IW][le * nEquation + IU] = 0;
        dRdq[re * nEquation + IP][le * nEquation + IU] -= dFluxdCv[IP] * area;

        //! perturb the momentum rhov of the left cell of the face
        qCvLpert[IR] = qCvL[IR];
        qCvLpert[IU] = qCvL[IU];
        qCvLpert[IV] = qCvL[IV];
        qCvLpert[IW] = qCvL[IW];
        qCvLpert[IP] = qCvL[IP];
        qCvLpert[IV] += perturbScale;
        dCvpert = qCvLpert[IV] - qCvL[IV];

        gas->Conservative2Primitive(qCvLpert, gmL, primL, temperaturesL);
        gas->Conservative2Primitive(qCvR, gmR, primR, temperaturesR);

        Cal_GMRES_Roe_Scheme_Flux(primL, primR, qCvLpert, qCvR, iFace, fluxp,
            face_proxy, invSchemePara);

        dFluxdCv[IR] = (fluxp[IR] - flux[IR][iFace]) / dCvpert;
        dFluxdCv[IU] = (fluxp[IU] - flux[IU][iFace]) / dCvpert;
        dFluxdCv[IV] = (fluxp[IV] - flux[IV][iFace]) / dCvpert;
        dFluxdCv[IW] = (fluxp[IW] - flux[IW][iFace]) / dCvpert;
        dFluxdCv[IP] = (fluxp[IP] - flux[IP][iFace]) / dCvpert;

        dRdq[le * nEquation + IR][le * nEquation + IV] += dFluxdCv[IR] * area;
        dRdq[le * nEquation + IU][le * nEquation + IV] += dFluxdCv[IU] * area;
        dRdq[le * nEquation + IV][le * nEquation + IV] += dFluxdCv[IV] * area;
        dRdq[le * nEquation + IW][le * nEquation + IV] = 0;
        dRdq[le * nEquation + IP][le * nEquation + IV] += dFluxdCv[IP] * area;

        dRdq[re * nEquation + IR][le * nEquation + IV] -= dFluxdCv[IR] * area;
        dRdq[re * nEquation + IU][le * nEquation + IV] -= dFluxdCv[IU] * area;
        dRdq[re * nEquation + IV][le * nEquation + IV] -= dFluxdCv[IV] * area;
        dRdq[re * nEquation + IW][le * nEquation + IV] = 0;
        dRdq[re * nEquation + IP][le * nEquation + IV] -= dFluxdCv[IP] * area;

        //! perturb the momentum rhow of the left cell of the face
        qCvLpert[IR] = qCvL[IR];
        qCvLpert[IU] = qCvL[IU];
        qCvLpert[IV] = qCvL[IV];
        qCvLpert[IW] = qCvL[IW];
        qCvLpert[IP] = qCvL[IP];
        qCvLpert[IW] += perturbScale;
        dCvpert = qCvLpert[IW] - qCvL[IW];

        gas->Conservative2Primitive(qCvLpert, gmL, primL, temperaturesL);
        gas->Conservative2Primitive(qCvR, gmR, primR, temperaturesR);

        Cal_GMRES_Roe_Scheme_Flux(primL, primR, qCvLpert, qCvR, iFace, fluxp,
            face_proxy, invSchemePara);

        dFluxdCv[IR] = (fluxp[IR] - flux[IR][iFace]) / dCvpert;
        dFluxdCv[IU] = (fluxp[IU] - flux[IU][iFace]) / dCvpert;
        dFluxdCv[IV] = (fluxp[IV] - flux[IV][iFace]) / dCvpert;
        dFluxdCv[IW] = (fluxp[IW] - flux[IW][iFace]) / dCvpert;
        dFluxdCv[IP] = (fluxp[IP] - flux[IP][iFace]) / dCvpert;

        dRdq[le * nEquation + IR][le * nEquation + IW] = 0;
        dRdq[le * nEquation + IU][le * nEquation + IW] = 0;
        dRdq[le * nEquation + IV][le * nEquation + IW] = 0;
        dRdq[le * nEquation + IW][le * nEquation + IW] += 1;
        dRdq[le * nEquation + IP][le * nEquation + IW] = 0;

        dRdq[re * nEquation + IR][le * nEquation + IW] = 0;
        dRdq[re * nEquation + IU][le * nEquation + IW] = 0;
        dRdq[re * nEquation + IV][le * nEquation + IW] = 0;
        dRdq[re * nEquation + IW][le * nEquation + IW] -= 1;
        dRdq[re * nEquation + IP][le * nEquation + IW] = 0;

        //! perturb the rhoE of the left cell of the face
        qCvLpert[IR] = qCvL[IR];
        qCvLpert[IU] = qCvL[IU];
        qCvLpert[IV] = qCvL[IV];
        qCvLpert[IW] = qCvL[IW];
        qCvLpert[IP] = qCvL[IP];
        qCvLpert[IP] *= (1.0+perturbScale);
        dCvpert = qCvLpert[IP] - qCvL[IP];

        gas->Conservative2Primitive(qCvLpert, gmL, primL, temperaturesL);
        gas->Conservative2Primitive(qCvR, gmR, primR, temperaturesR);

        Cal_GMRES_Roe_Scheme_Flux(primL, primR, qCvLpert, qCvR, iFace, fluxp,
            face_proxy, invSchemePara);

        dFluxdCv[IR] = (fluxp[IR] - flux[IR][iFace]) / dCvpert;
        dFluxdCv[IU] = (fluxp[IU] - flux[IU][iFace]) / dCvpert;
        dFluxdCv[IV] = (fluxp[IV] - flux[IV][iFace]) / dCvpert;
        dFluxdCv[IW] = (fluxp[IW] - flux[IW][iFace]) / dCvpert;
        dFluxdCv[IP] = (fluxp[IP] - flux[IP][iFace]) / dCvpert;

        dRdq[le * nEquation + IR][le * nEquation + IP] += dFluxdCv[IR] * area;
        dRdq[le * nEquation + IU][le * nEquation + IP] += dFluxdCv[IU] * area;
        dRdq[le * nEquation + IV][le * nEquation + IP] += dFluxdCv[IV] * area;
        dRdq[le * nEquation + IW][le * nEquation + IP] = 0;
        dRdq[le * nEquation + IP][le * nEquation + IP] += dFluxdCv[IP] * area;

        dRdq[re * nEquation + IR][le * nEquation + IP] -= dFluxdCv[IR] * area;
        dRdq[re * nEquation + IU][le * nEquation + IP] -= dFluxdCv[IU] * area;
        dRdq[re * nEquation + IV][le * nEquation + IP] -= dFluxdCv[IV] * area;
        dRdq[re * nEquation + IW][le * nEquation + IP] = 0;
        dRdq[re * nEquation + IP][le * nEquation + IP] -= dFluxdCv[IP] * area;

        //! perturb the density of the right cell of the face
        qCvRpert[IR] = qCvR[IR];
        qCvRpert[IU] = qCvR[IU];
        qCvRpert[IV] = qCvR[IV];
        qCvRpert[IW] = qCvR[IW];
        qCvRpert[IP] = qCvR[IP];
        qCvRpert[IR] *= (1.0+perturbScale);
        dCvpert = qCvRpert[IR] - qCvR[IR];

        gas->Conservative2Primitive(qCvL, gmL, primL, temperaturesL);
        gas->Conservative2Primitive(qCvRpert, gmR, primR, temperaturesR);

        Cal_GMRES_Roe_Scheme_Flux(primL, primR, qCvL, qCvRpert, iFace, fluxp,
            face_proxy, invSchemePara);

        dFluxdCv[IR] = (fluxp[IR] - flux[IR][iFace]) / dCvpert;
        dFluxdCv[IU] = (fluxp[IU] - flux[IU][iFace]) / dCvpert;
        dFluxdCv[IV] = (fluxp[IV] - flux[IV][iFace]) / dCvpert;
        dFluxdCv[IW] = (fluxp[IW] - flux[IW][iFace]) / dCvpert;
        dFluxdCv[IP] = (fluxp[IP] - flux[IP][iFace]) / dCvpert;

        dRdq[le * nEquation + IR][re * nEquation + IR] += dFluxdCv[IR] * area;
        dRdq[le * nEquation + IU][re * nEquation + IR] += dFluxdCv[IU] * area;
        dRdq[le * nEquation + IV][re * nEquation + IR] += dFluxdCv[IV] * area;
        dRdq[le * nEquation + IW][re * nEquation + IR] = 0;
        dRdq[le * nEquation + IP][re * nEquation + IR] += dFluxdCv[IP] * area;

        dRdq[re * nEquation + IR][re * nEquation + IR] -= dFluxdCv[IR] * area;
        dRdq[re * nEquation + IU][re * nEquation + IR] -= dFluxdCv[IU] * area;
        dRdq[re * nEquation + IV][re * nEquation + IR] -= dFluxdCv[IV] * area;
        dRdq[re * nEquation + IW][re * nEquation + IR] = 0;
        dRdq[re * nEquation + IP][re * nEquation + IR] -= dFluxdCv[IP] * area;

        //! perturb the momentum rhou of the right cell of the face
        qCvRpert[IR] = qCvR[IR];
        qCvRpert[IU] = qCvR[IU];
        qCvRpert[IV] = qCvR[IV];
        qCvRpert[IW] = qCvR[IW];
        qCvRpert[IP] = qCvR[IP];
        qCvRpert[IU] += perturbScale;
        dCvpert = qCvRpert[IU] - qCvR[IU];

        gas->Conservative2Primitive(qCvL, gmL, primL, temperaturesL);
        gas->Conservative2Primitive(qCvRpert, gmR, primR, temperaturesR);

        Cal_GMRES_Roe_Scheme_Flux(primL, primR, qCvL, qCvRpert, iFace, fluxp,
                                  face_proxy, invSchemePara);

        dFluxdCv[IR] = (fluxp[IR] - flux[IR][iFace]) / dCvpert;
        dFluxdCv[IU] = (fluxp[IU] - flux[IU][iFace]) / dCvpert;
        dFluxdCv[IV] = (fluxp[IV] - flux[IV][iFace]) / dCvpert;
        dFluxdCv[IW] = (fluxp[IW] - flux[IW][iFace]) / dCvpert;
        dFluxdCv[IP] = (fluxp[IP] - flux[IP][iFace]) / dCvpert;

        dRdq[le * nEquation + IR][re * nEquation + IU] += dFluxdCv[IR] * area;
        dRdq[le * nEquation + IU][re * nEquation + IU] += dFluxdCv[IU] * area;
        dRdq[le * nEquation + IV][re * nEquation + IU] += dFluxdCv[IV] * area;
        dRdq[le * nEquation + IW][re * nEquation + IU] = 0;
        dRdq[le * nEquation + IP][re * nEquation + IU] += dFluxdCv[IP] * area;

        dRdq[re * nEquation + IR][re * nEquation + IU] -= dFluxdCv[IR] * area;
        dRdq[re * nEquation + IU][re * nEquation + IU] -= dFluxdCv[IU] * area;
        dRdq[re * nEquation + IV][re * nEquation + IU] -= dFluxdCv[IV] * area;
        dRdq[re * nEquation + IW][re * nEquation + IU] = 0;
        dRdq[re * nEquation + IP][re * nEquation + IU] -= dFluxdCv[IP] * area;

        //! perturb the momentum rhov of the right cell of the face
        qCvRpert[IR] = qCvR[IR];
        qCvRpert[IU] = qCvR[IU];
        qCvRpert[IV] = qCvR[IV];
        qCvRpert[IW] = qCvR[IW];
        qCvRpert[IP] = qCvR[IP];
        qCvRpert[IV] += perturbScale;
        dCvpert = qCvRpert[IV] - qCvR[IV];

        gas->Conservative2Primitive(qCvL, gmL, primL, temperaturesL);
        gas->Conservative2Primitive(qCvRpert, gmR, primR, temperaturesR);

        Cal_GMRES_Roe_Scheme_Flux(primL, primR, qCvL, qCvRpert, iFace, fluxp,
                                  face_proxy, invSchemePara);

        dFluxdCv[IR] = (fluxp[IR] - flux[IR][iFace]) / dCvpert;
        dFluxdCv[IU] = (fluxp[IU] - flux[IU][iFace]) / dCvpert;
        dFluxdCv[IV] = (fluxp[IV] - flux[IV][iFace]) / dCvpert;
        dFluxdCv[IW] = (fluxp[IW] - flux[IW][iFace]) / dCvpert;
        dFluxdCv[IP] = (fluxp[IP] - flux[IP][iFace]) / dCvpert;

        dRdq[le * nEquation + IR][re * nEquation + IV] += dFluxdCv[IR] * area;
        dRdq[le * nEquation + IU][re * nEquation + IV] += dFluxdCv[IU] * area;
        dRdq[le * nEquation + IV][re * nEquation + IV] += dFluxdCv[IV] * area;
        dRdq[le * nEquation + IW][re * nEquation + IV] = 0;
        dRdq[le * nEquation + IP][re * nEquation + IV] += dFluxdCv[IP] * area;

        dRdq[re * nEquation + IR][re * nEquation + IV] -= dFluxdCv[IR] * area;
        dRdq[re * nEquation + IU][re * nEquation + IV] -= dFluxdCv[IU] * area;
        dRdq[re * nEquation + IV][re * nEquation + IV] -= dFluxdCv[IV] * area;
        dRdq[re * nEquation + IW][re * nEquation + IV] = 0;
        dRdq[re * nEquation + IP][re * nEquation + IV] -= dFluxdCv[IP] * area;

        //! perturb the momentum rhow of the right cell of the face
        qCvRpert[IR] = qCvR[IR];
        qCvRpert[IU] = qCvR[IU];
        qCvRpert[IV] = qCvR[IV];
        qCvRpert[IW] = qCvR[IW];
        qCvRpert[IP] = qCvR[IP];
        qCvRpert[IW] += perturbScale;
        dCvpert = qCvRpert[IW] - qCvR[IW];

        gas->Conservative2Primitive(qCvL, gmL, primL, temperaturesL);
        gas->Conservative2Primitive(qCvRpert, gmR, primR, temperaturesR);

        Cal_GMRES_Roe_Scheme_Flux(primL, primR, qCvL, qCvRpert, iFace, fluxp,
            face_proxy, invSchemePara);

        dFluxdCv[IR] = (fluxp[IR] - flux[IR][iFace]) / dCvpert;
        dFluxdCv[IU] = (fluxp[IU] - flux[IU][iFace]) / dCvpert;
        dFluxdCv[IV] = (fluxp[IV] - flux[IV][iFace]) / dCvpert;
        dFluxdCv[IW] = (fluxp[IW] - flux[IW][iFace]) / dCvpert;
        dFluxdCv[IP] = (fluxp[IP] - flux[IP][iFace]) / dCvpert;

        dRdq[le * nEquation + IR][re * nEquation + IW] = 0;
        dRdq[le * nEquation + IU][re * nEquation + IW] = 0;
        dRdq[le * nEquation + IV][re * nEquation + IW] = 0;
        dRdq[le * nEquation + IW][re * nEquation + IW] += 1;
        dRdq[le * nEquation + IP][re * nEquation + IW] = 0;

        dRdq[re * nEquation + IR][re * nEquation + IW] = 0;
        dRdq[re * nEquation + IU][re * nEquation + IW] = 0;
        dRdq[re * nEquation + IV][re * nEquation + IW] = 0;
        dRdq[re * nEquation + IW][re * nEquation + IW] -= 1;
        dRdq[re * nEquation + IP][re * nEquation + IW] = 0;

        //! perturb the momentum rhoE of the right cell of the face
        qCvRpert[IR] = qCvR[IR];
        qCvRpert[IU] = qCvR[IU];
        qCvRpert[IV] = qCvR[IV];
        qCvRpert[IW] = qCvR[IW];
        qCvRpert[IP] = qCvR[IP];
        qCvRpert[IP] *= (1.0+perturbScale);
        dCvpert = qCvRpert[IP] - qCvR[IP];

        gas->Conservative2Primitive(qCvL, gmL, primL, temperaturesL);
        gas->Conservative2Primitive(qCvRpert, gmR, primR, temperaturesR);

        Cal_GMRES_Roe_Scheme_Flux(primL, primR, qCvL, qCvRpert, iFace, fluxp,
            face_proxy, invSchemePara);

        dFluxdCv[IR] = (fluxp[IR] - flux[IR][iFace]) / dCvpert;
        dFluxdCv[IU] = (fluxp[IU] - flux[IU][iFace]) / dCvpert;
        dFluxdCv[IV] = (fluxp[IV] - flux[IV][iFace]) / dCvpert;
        dFluxdCv[IW] = (fluxp[IW] - flux[IW][iFace]) / dCvpert;
        dFluxdCv[IP] = (fluxp[IP] - flux[IP][iFace]) / dCvpert;

        dRdq[le * nEquation + IR][re * nEquation + IP] += dFluxdCv[IR] * area;
        dRdq[le * nEquation + IU][re * nEquation + IP] += dFluxdCv[IU] * area;
        dRdq[le * nEquation + IV][re * nEquation + IP] += dFluxdCv[IV] * area;
        dRdq[le * nEquation + IW][re * nEquation + IP] = 0;
        dRdq[le * nEquation + IP][re * nEquation + IP] += dFluxdCv[IP] * area;

        dRdq[re * nEquation + IR][re * nEquation + IP] -= dFluxdCv[IR] * area;
        dRdq[re * nEquation + IU][re * nEquation + IP] -= dFluxdCv[IU] * area;
        dRdq[re * nEquation + IV][re * nEquation + IP] -= dFluxdCv[IV] * area;
        dRdq[re * nEquation + IW][re * nEquation + IP] = 0;
        dRdq[re * nEquation + IP][re * nEquation + IP] -= dFluxdCv[IP] * area;
    }

    delete [] prim;          prim = nullptr;
    delete [] primL;         primL = nullptr;
    delete [] primR;         primR = nullptr;
    delete [] qConserveL;    qConserveL = nullptr;
    delete [] qConserveR;    qConserveR = nullptr;
    delete [] qCvL;          qCvL = nullptr;
    delete [] dFluxdCv;      dFluxdCv = nullptr;
    delete [] qCvLpert;      qCvLpert = nullptr;
    delete [] qCvRpert;      qCvRpert = nullptr;
    delete [] qCvR;          qCvR = nullptr;
    delete [] fluxp;         fluxp = nullptr;
}

//! GMRESnolim
void cal_ReConstructFaceValue(RDouble* primcc, RDouble* dqdx, RDouble* dqdy, RDouble* dqdz,
                              RDouble dx, RDouble dy, RDouble dz, RDouble * primf,
                              InviscidSchemeParameter* invSchemePara)
{
    int nEquation = invSchemePara->GetNumberOfTotalEquation();

    for (int m = 0; m < nEquation; ++ m)
    {
        primf[m] = primcc[m];
        primf[m] += dqdx[m] * dx + dqdy[m] * dy + dqdz[m] * dz;
    }

    if(primf[IDX::IR] <= 0 && primf[IDX::IP] <= 0)
    {
        for (int m = 0; m < nEquation; ++ m)
        {
            primf[m] = primcc[m];
        }
    }
}

//! GMRESPV GMRESAD GMRES3D
template <typename T>
void Cal_GMRES_Steger_Scheme_Flux_PV(T *primL, T *primR, int iLen, T *flux,
                               FaceProxy *face_proxy, InviscidSchemeParameter *invSchemePara)
{
    int nTemperatureModel = invSchemePara->GetNumberOfTemperatureModel();
    int nEquation = invSchemePara->GetNumberOfTotalEquation();
    int nLaminar = invSchemePara->GetNumberOfLaminar();
    int nm = invSchemePara->Getnm();
    int EntropyFixMethod = invSchemePara->GetEntropyFixMethod();
    int nChemical = invSchemePara->GetNumberOfChemical();
    int nIdealState = invSchemePara->GetIfIdealState();

    int nElectronIndex = invSchemePara->GetIndexOfElectron();

    RDouble entrFixCoefLeft  = invSchemePara->GetLeftEntropyFixCoefficient();
    RDouble entrFixCoefRight = invSchemePara->GetRightEntropyFixCoefficient();
    RDouble *gamal = invSchemePara->GetLeftGama();
    RDouble *gamar = invSchemePara->GetRightGama();
    RDouble *xfn   = invSchemePara->GetFaceNormalX();
    RDouble *yfn   = invSchemePara->GetFaceNormalY();
    RDouble *zfn   = invSchemePara->GetFaceNormalZ();
    RDouble *vgn   = invSchemePara->GetFaceVelocity();
    RDouble **tl = invSchemePara->GetLeftTemperature();
    RDouble **tr = invSchemePara->GetRightTemperature();

    using namespace IDX;
    using namespace GAS_SPACE;

    T *fl = new T [nEquation];
    T *fr = new T [nEquation];
    T temperaturesL[3], temperaturesR[3];

    T rl = primL[IR];
    T ul = primL[IU];
    T vl = primL[IV];
    T wl = primL[IW];
    T pl = primL[IP];

    T rr = primR[IR];
    T ur = primR[IU];
    T vr = primR[IV];
    T wr = primR[IW];
    T pr = primR[IP];

    RDouble xsn = xfn[iLen];
    RDouble ysn = yfn[iLen];
    RDouble zsn = zfn[iLen];
    RDouble vb  = vgn[iLen];
    RDouble gmL = gamal[iLen];
    RDouble gmR = gamar[iLen];

    T v2l = ul * ul + vl * vl + wl * wl;
    T v2r = ur * ur + vr * vr + wr * wr;
    T hintLeft  = (gmL / (gmL - one)) * (pl / rl);
    T hintRight = (gmR / (gmR - one)) * (pr / rr);
    T pel = 0.0, per = 0.0;

    if (nChemical != 0 && nIdealState == 0)
    {
        /*
        for (int i = 0; i < 3; ++ i)
        {
            temperaturesL[i] = tl[i][iLen];
            temperaturesR[i] = tr[i][iLen];
        }
        gas->ComputeEnthalpyByPrimitive(priml, gmL, hintLeft, temperaturesL);
        gas->ComputeEnthalpyByPrimitive(primr, gmR, hintRight, temperaturesR);
        if (nTemperatureModel > 1)
        {
            pel = gas->GetElectronPressure(priml, temperaturesL[ITE]);
            per = gas->GetElectronPressure(primr, temperaturesR[ITE]);
        }*/
    }

    T hl = hintLeft  + half * v2l;
    T hr = hintRight + half * v2r;
    T el = hl - pl / (rl + SMALL);
    T er = hr - pr / (rr + SMALL);
    T vnl = xsn * ul + ysn * vl + zsn * wl - vb;
    T vnr = xsn * ur + ysn * vr + zsn * wr - vb;
    T c2l = gmL * pl / rl;
    T cl  = sqrt(ABS(c2l));
    T c2r = gmR * pr / rr;
    T cr  = sqrt(ABS(c2r));
    //! The eigenvalues of left and right.
    T eigvl1 = vnl;
    T eigvl2 = vnl + cl;
    T eigvl3 = vnl - cl;
    T eigvr1 = vnr;
    T eigvr2 = vnr + cr;
    T eigvr3 = vnr - cr;
    T absEigvl1 = ABS(eigvl1);
    T absEigvl2 = ABS(eigvl2);
    T absEigvl3 = ABS(eigvl3);
    T absEigvr1 = ABS(eigvr1);
    T absEigvr2 = ABS(eigvr2);
    T absEigvr3 = ABS(eigvr3);
    //! Entropy fix.
    //! Modified by clz : 2012-8-15.
    if (EntropyFixMethod == 1)
    {
        T maxEigv = MAX(absEigvl2, absEigvl3);
        T tmp0    = maxEigv * entrFixCoefLeft ;
        T tmp1    = maxEigv * entrFixCoefRight;
        absEigvl1 = MAX(tmp0, absEigvl1);
        absEigvl2 = MAX(tmp1, absEigvl2);
        absEigvl3 = MAX(tmp1, absEigvl3);
        maxEigv   = MAX(absEigvr2, absEigvr3);
        tmp0      = maxEigv * entrFixCoefLeft ;
        tmp1      = maxEigv * entrFixCoefRight;
        absEigvr1 = MAX(tmp0, absEigvr1);
        absEigvr2 = MAX(tmp1, absEigvr2);
        absEigvr3 = MAX(tmp1, absEigvr3);
    }
    else if (EntropyFixMethod == 2)
    {
        //RDouble tmp0 = MAX(entrFixCoefLeft , static_cast <RDouble> (0.02));
        //RDouble tmp1 = MAX(entrFixCoefRight, static_cast <RDouble> (0.02));
        //absEigvl1 = MAX(tmp0, absEigvl1);
        //absEigvl2 = MAX(tmp1, absEigvl2);
        //absEigvl3 = MAX(tmp1, absEigvl3);
               
        //absEigvr1 = MAX(tmp0, absEigvr1);
        //absEigvr2 = MAX(tmp1, absEigvr2);
        //absEigvr3 = MAX(tmp1, absEigvr3);
        absEigvl1 = sqrt(absEigvl1*absEigvl1 + 1.E-4);
        absEigvl2 = sqrt(absEigvl2*absEigvl2 + 1.E-4);
        absEigvl3 = sqrt(absEigvl3*absEigvl3 + 1.E-4);
        absEigvr1 = sqrt(absEigvr1*absEigvr1 + 1.E-4);
        absEigvr2 = sqrt(absEigvr2*absEigvr2 + 1.E-4);
        absEigvr3 = sqrt(absEigvr3*absEigvr3 + 1.E-4);
    }
    eigvl1 = half * (eigvl1 + absEigvl1);
    eigvl2 = half * (eigvl2 + absEigvl2);
    eigvl3 = half * (eigvl3 + absEigvl3);
    eigvr1 = half * (eigvr1 - absEigvr1);
    eigvr2 = half * (eigvr2 - absEigvr2);
    eigvr3 = half * (eigvr3 - absEigvr3);
    T c2gaml = c2l / gmL;
    T c2gamr = c2r / gmR;
    T xil1 = c2gaml * (eigvl1 + eigvl1 - eigvl2 - eigvl3) / (c2l + c2l);
    T xil2 = c2gaml * (eigvl2 - eigvl3) / (cl + cl);
    T xir1 = c2gamr * (eigvr1 + eigvr1 - eigvr2 - eigvr3) / (c2r + c2r);
    T xir2 = c2gamr * (eigvr2 - eigvr3) / (cr + cr);
    T eig_xi = eigvl1 - xil1;
    fl[IR ] = rl * (eig_xi                );
    fl[IRU] = rl * (eig_xi * ul + xil2 * xsn);
    fl[IRV] = rl * (eig_xi * vl + xil2 * ysn);
    fl[IRW] = rl * (eig_xi * wl + xil2 * zsn);
    fl[IRE] = rl * (eigvl1 * el - xil1 * hl + xil2 * (vnl + vb));
    T fel = - xil1 * pel;
    eig_xi = eigvr1 - xir1;
    fr[IR ] = rr * (eig_xi                );
    fr[IRU] = rr * (eig_xi * ur + xir2 * xsn);
    fr[IRV] = rr * (eig_xi * vr + xir2 * ysn);
    fr[IRW] = rr * (eig_xi * wr + xir2 * zsn);
    fr[IRE] = rr * (eigvr1 * er - xir1 * hr + xir2 * (vnr + vb));
    T fer = - xir1 * per;
    for (int iLaminar = nm; iLaminar < nEquation; ++ iLaminar)
    {
        fl[iLaminar] = primL[iLaminar] * fl[IR];
        fr[iLaminar] = primR[iLaminar] * fr[IR];
    }
    if (nChemical == 1)
    {
        fl[nLaminar + nTemperatureModel - 1] += fel;
        fr[nLaminar + nTemperatureModel - 1] += fer;
    }
    //! Compute vibration-electron energy term.
    for (int iLaminar = 0; iLaminar < nEquation; ++ iLaminar)
    {
        flux[iLaminar] = fl[iLaminar] + fr[iLaminar];
    }
    //! Modify the fluxes of the species whose transport equations need not to be computed.
    if (nChemical == 1)
    {
    
    }

    delete [] fl;    fl = nullptr;
    delete [] fr;    fr = nullptr;
}

//! GMRESPV GMRESAD
template <typename T>
void Cal_GMRES_Roe_Scheme_Flux_PV(T *primL, T *primR, int iFace, T *flux,
                               FaceProxy *face_proxy, InviscidSchemeParameter *invSchemePara)
{
    using namespace IDX;
    using namespace GAS_SPACE;
    const RDouble LITTLE     = 1.0E-10;
    int nTemperatureModel    = invSchemePara->GetNumberOfTemperatureModel();
    int nEquation            = invSchemePara->GetNumberOfTotalEquation();
    int nLaminar             = invSchemePara->GetNumberOfLaminar();
    int nm                   = invSchemePara->Getnm();
    RDouble *gamaL           = invSchemePara->GetLeftGama();
    RDouble *gamaR           = invSchemePara->GetRightGama();
    RDouble *pressureCoeffL  = face_proxy->GetPressureCoefficientL();
    RDouble *pressureCoeffR  = face_proxy->GetPressureCoefficientR();
    T entrFixCoefLeft        = invSchemePara->GetLeftEntropyFixCoefficient();
    T entrFixCoefRight       = invSchemePara->GetRightEntropyFixCoefficient();
    RDouble *xfn             = invSchemePara->GetFaceNormalX();
    RDouble *yfn             = invSchemePara->GetFaceNormalY();
    RDouble *zfn             = invSchemePara->GetFaceNormalZ();
    RDouble *vgn             = invSchemePara->GetFaceVelocity();
    RDouble refMachNumber    = invSchemePara->GetMachNumber();
    RDouble **tl             = invSchemePara->GetLeftTemperature();
    RDouble **tr             = invSchemePara->GetRightTemperature();
    int RoeEntropyFixMethod  = invSchemePara->GetEntropyFixMethod();
    int nChemical            = invSchemePara->GetNumberOfChemical();
    int nPrecondition        = invSchemePara->GetIfPrecondition();
    int nIdealState          = invSchemePara->GetIfIdealState();
    RDouble mach2            = refMachNumber * refMachNumber;
    T *prim                  = new T[nEquation]();
    RDouble *qConserveL      = new RDouble[nEquation]();
    RDouble *qConserveR      = new RDouble[nEquation]();
    T *dq                    = new T[nEquation]();
    RDouble temperaturesL[3] = {0}, temperaturesR[3] = {0};
    RDouble xsn              = xfn[iFace];
    RDouble ysn              = yfn[iFace];
    RDouble zsn              = zfn[iFace];
    RDouble vb               = vgn[iFace];
    RDouble gmL              = gamaL[iFace];
    RDouble gmR              = gamaR[iFace];
    T half = 0.5;
    T two = 2.0;
    T one = 1.0;
    T three = 3.0;

    T& rL     = primL[IR];
    T& uL     = primL[IU];
    T& vL     = primL[IV];
    T& wL     = primL[IW];
    T& pL     = primL[IP];

    T& rR     = primR[IR];
    T& uR     = primR[IU];
    T& vR     = primR[IV];
    T& wR     = primR[IW];
    T& pR     = primR[IP];

    rL = MAX(rL, LITTLE);
    rR = MAX(rR, LITTLE);
    pL = MAX(pL, LITTLE);
    pR = MAX(pR, LITTLE);

    T vSquareL = uL * uL + vL * vL + wL * wL;
    T vSquareR = uR * uR + vR * vR + wR * wR;

    T hint_L  = (gmL / (gmL - one)) * (pL / rL);
    T hint_R  = (gmR / (gmR - one)) * (pR / rR);
    T peL     = 0.0, peR = 0.0;

    if (nChemical != 0 && nIdealState == 0)
    {
       // for (int i = 0; i < 3; ++i)
       // {
       //     temperaturesL[i] = tl[i][iFace];
       //     temperaturesR[i] = tr[i][iFace];
       // }
       // gas->ComputeEnthalpyByPrimitive(primL, gmL, hint_L, temperaturesL);
       // gas->ComputeEnthalpyByPrimitive(primR, gmR, hint_R, temperaturesR);

       // if (nTemperatureModel > 1)
       // {
       //     peL = gas->GetElectronPressure(primL, temperaturesL[ITE]);
       //     peR = gas->GetElectronPressure(primR, temperaturesR[ITE]);
       // }
    }
    T hL = hint_L + half * vSquareL;
    T hR = hint_R + half * vSquareR;

    T vnL = xsn * uL + ysn * vL + zsn * wL - vb;
    T vnR = xsn * uR + ysn * vR + zsn * wR - vb;
    T rvnL = rL * vnL;
    T rvnR = rR * vnR;

    //! The first part of Roe scheme: F(QL, QR) = 0.5 * [ F(QL) + F(QR) ].
    flux[IR]  = half * (rvnL + rvnR);
    flux[IRU] = half * (rvnL * uL + xfn[iFace] * pL + rvnR * uR + xfn[iFace] * pR);
    flux[IRV] = half * (rvnL * vL + yfn[iFace] * pL + rvnR * vR + yfn[iFace] * pR);
    flux[IRW] = half * (rvnL * wL + zfn[iFace] * pL + rvnR * wR + zfn[iFace] * pR);
    flux[IRE] = half * (rvnL * hL + vgn[iFace] * pL + rvnR * hR + vgn[iFace] * pR);

    //! The following loop is useless.
    for (int iLaminar = nm; iLaminar < nLaminar; ++iLaminar)
    {
        flux[iLaminar] = half * (primL[iLaminar] * rvnL + primR[iLaminar] * rvnR);
    }

    for (int it = 1; it < nTemperatureModel; ++it)
    {
        flux[nLaminar + it] = half * (primL[nLaminar + it] * rvnL + primR[nLaminar + it] * rvnR);
    }
    if (nChemical == 1)
    {
        flux[nLaminar + nTemperatureModel - 1] += half * (peL * (vnL + vgn[iFace]) + peR * (vnR + vgn[iFace]));
    }

    T ratio = sqrt(rR / rL);
    T roeAverageCoef = one / (one + ratio);

    T rRoeAverage = sqrt(rL * rR);
    T uRoeAverage = (uL + uR * ratio) * roeAverageCoef;
    T vRoeAverage = (vL + vR * ratio) * roeAverageCoef;
    T wRoeAverage = (wL + wR * ratio) * roeAverageCoef;
    T pRoeAverage = (pL + pR * ratio) * roeAverageCoef;

    T gama = (gmL + gmR * ratio) * roeAverageCoef;
    T gamm1 = gama - one;

    T vSquareRoeAverage = uRoeAverage * uRoeAverage + vRoeAverage * vRoeAverage + wRoeAverage * wRoeAverage;
    T absVel = sqrt(vSquareRoeAverage);
    T hRoeAverage = gama / gamm1 * pRoeAverage / rRoeAverage + half * vSquareRoeAverage;

    T vn = xsn * uRoeAverage + ysn * vRoeAverage + zsn * wRoeAverage - vb;
    T theta = vn + vb;

    //! c2: square of sound speed.
    T c2 = gamm1 * (hRoeAverage - half * vSquareRoeAverage);

    //! Sound speed: gama * pRoeAverage / rRoeAverage.
    T cm = sqrt(abs(c2));    //! HINT: from ABS to abs

    T Ur2 = 0, cnew2 = 0, cnew = 0, vnew = 0, beta = 0, alpha = 0;
    T eigv1 = 0, eigv2 = 0, eigv3 = 0;
    eigv1 = abs(vn); //HINT: from ABS to abs
    if (nPrecondition == 0)
    {
        eigv2 = abs(vn + cm);    //! HINT: from ABS to abs
        eigv3 = abs(vn - cm);    //! HINT: from ABS to abs
    }
    else
    {
        beta = min(max(vSquareRoeAverage / c2, three * mach2), one);    //! HINT: MAX MIN
        vnew = half * vn * (one + beta);
        cnew = half * sqrt(((beta - one) * vn) * ((beta - one) * vn) + four * beta * c2);
        cnew2 = cnew * cnew;
        eigv2 = abs(vnew + cnew);    //! HINT: from ABS to abs
        eigv3 = abs(vnew - cnew);    //! HINT: from ABS to abs
    }

    T EL = (one / (gmL - one)) * (pL / rL) + half * vSquareL;
    T ER = (one / (gmR - one)) * (pR / rR) + half * vSquareR;

    //! The jump between left and right conservative variables.
    dq[IR]  = rR - rL;
    dq[IRU] = rR * uR - rL * uL;
    dq[IRV] = rR * vR - rL * vL;
    dq[IRW] = rR * wR - rL * wL;
    dq[IRE] = rR * ER - rL * EL; 

    if (nChemical != 0)
    {
        // if (nTemperatureModel > 1)
        // {
        //     gas->GetTemperature(primL, temperaturesL[ITT], temperaturesL[ITV], temperaturesL[ITE]);
        //     gas->GetTemperature(primR, temperaturesR[ITT], temperaturesR[ITV], temperaturesR[ITE]);
        // }
        // gas->Primitive2Conservative(primL, gmL, temperaturesL[ITV], temperaturesL[ITE], qConserveL);
        // gas->Primitive2Conservative(primR, gmR, temperaturesR[ITV], temperaturesR[ITE], qConserveR);
    }

    for (int iLaminar = nm; iLaminar < nEquation; ++iLaminar)
    {
        dq[iLaminar] = qConserveR[iLaminar] - qConserveL[iLaminar];
    }

    //! Entropy correction, to limit the magnitude of the three eigenvalue.
    //! Warning: This function calling may affect the efficient!
    T pressureCoeffadL = pressureCoeffL[iFace];
    T pressureCoeffadR = pressureCoeffR[iFace];
    Flux_RoeEntropyCorrection(eigv1, eigv2, eigv3, pressureCoeffadL, pressureCoeffadR, absVel, vn, cm, RoeEntropyFixMethod, entrFixCoefLeft, entrFixCoefRight);

    T xi1 = 0.0, xi2 = 0.0, xi3 = 0.0, xi4 = 0.0;
    if (nPrecondition == 0)
    {
        xi1 = (two * eigv1 - eigv2 - eigv3) / (two * c2);
        xi2 = (eigv2 - eigv3) / (two * cm);
    }
    else
    {
        xi1 = (two * eigv1 - eigv2 - eigv3) / (two * cnew2);
        xi2 = (eigv2 - eigv3) / (two * cnew);
        alpha = half * (one - beta);
        xi3 = two * alpha * eigv1 / cnew2;
        xi4 = alpha * vn * xi2 / cnew2;
    }

    T dc = theta * dq[IR] - xsn * dq[IRU] - ysn * dq[IRV] - zsn * dq[IRW];

    T c2dc = 0;
    if (nPrecondition == 0)
    {
        c2dc = c2 * dc;
    }
    else
    {
        Ur2 = beta * c2;
        c2dc = cnew2 * dc;
    }

    T dh = -gamm1 * (uRoeAverage * dq[IRU] + vRoeAverage * dq[IRV] + wRoeAverage * dq[IRW] - dq[IRE]) + half * gamm1 * vSquareRoeAverage * dq[IR];

    if (nChemical != 0)
    {
        // gas->ComputeDHAndTotalEnthalpy(prim, gama, dq, dh, hRoeAverage);
    }

    if (nPrecondition == 0)
    {
        flux[IR]     -= half * (eigv1 * dq[IR] - dh * xi1 - dc * xi2);
        flux[IRU]    -= half * (eigv1 * dq[IRU] + (xsn * c2dc - uRoeAverage * dh) * xi1 + (xsn * dh - uRoeAverage * dc) * xi2);
        flux[IRV]    -= half * (eigv1 * dq[IRV] + (ysn * c2dc - vRoeAverage * dh) * xi1 + (ysn * dh - vRoeAverage * dc) * xi2);
        flux[IRW]    -= half * (eigv1 * dq[IRW] + (zsn * c2dc - wRoeAverage * dh) * xi1 + (zsn * dh - wRoeAverage * dc) * xi2);
        flux[IRE]    -= half * (eigv1 * dq[IRE] + (theta * c2dc - hRoeAverage * dh) * xi1 + (theta * dh - hRoeAverage * dc) * xi2);

        for (int iLaminar = nm; iLaminar < nLaminar; ++iLaminar)
        {
            prim[iLaminar] = half * (primL[iLaminar] + primR[iLaminar]);
            flux[iLaminar] -= half * (eigv1 * dq[iLaminar] - prim[iLaminar] * (dh * xi1 + dc * xi2));
        }
    }
    else
    {
        flux[IR]     -= half * (eigv1 * dq[IR] - dh * cnew2 / Ur2 * xi1 - dc * xi2 + dh * cnew2 / Ur2 * xi3 - dh * cnew2 / Ur2 * xi4);
        flux[IRU]    -= half * (eigv1 * dq[IRU] + (xsn * c2dc - uRoeAverage * dh * cnew2 / Ur2) * xi1 + (xsn * dh - uRoeAverage * dc) * xi2 + uRoeAverage * dh * cnew2 / Ur2 * xi3 - (xsn * c2dc + uRoeAverage * dh * cnew2 / Ur2) * xi4);
        flux[IRV]    -= half * (eigv1 * dq[IRV] + (ysn * c2dc - vRoeAverage * dh * cnew2 / Ur2) * xi1 + (ysn * dh - vRoeAverage * dc) * xi2 + vRoeAverage * dh * cnew2 / Ur2 * xi3 - (ysn * c2dc + vRoeAverage * dh * cnew2 / Ur2) * xi4);
        flux[IRW]    -= half * (eigv1 * dq[IRW] + (zsn * c2dc - wRoeAverage * dh * cnew2 / Ur2) * xi1 + (zsn * dh - wRoeAverage * dc) * xi2 + wRoeAverage * dh * cnew2 / Ur2 * xi3 - (zsn * c2dc + wRoeAverage * dh * cnew2 / Ur2) * xi4);
        flux[IRE]    -= half * (eigv1 * dq[IRE] + (theta * c2dc - hRoeAverage * dh * cnew2 / Ur2) * xi1 + (theta * dh - hRoeAverage * dc) * xi2 + hRoeAverage * dh * cnew2 / Ur2 * xi3 - (theta * c2dc + hRoeAverage * dh * cnew2 / Ur2) * xi4);

        for (int iLaminar = nm; iLaminar < nLaminar; ++iLaminar)
        {
            prim[iLaminar] = half * (primL[iLaminar] + primR[iLaminar]);
            flux[iLaminar] -= half * (eigv1 * dq[iLaminar] - prim[iLaminar] * (dh * cnew2 / Ur2 * (xi1 - xi3 + xi4) + dc * xi2));
        }
        //! Compute vibration-electron energy term.
        for (int it = 1; it < nTemperatureModel; ++it)
        {
            prim[nLaminar + it] = half * (primL[nLaminar + it] + primR[nLaminar + it]);
            flux[nLaminar + it] -= half * (eigv1 * dq[nLaminar + it] - prim[nLaminar + it] * (dh * cnew2 / Ur2 * (xi1 - xi3 + xi4) + dc * xi2));
        }
    }

    delete [] dq;         dq = nullptr;
    delete [] qConserveL; qConserveL = nullptr;
    delete [] qConserveR; qConserveR = nullptr;
    delete [] prim;       prim = nullptr;
}

//! GMRES 
void Cal_GMRES_Roe_Scheme_Flux(RDouble *primL, RDouble *primR, RDouble *qCvL, RDouble *qCvR, int iFace, RDouble *flux,
                               FaceProxy *face_proxy, InviscidSchemeParameter *invSchemePara)
{
    using namespace IDX;
    using namespace GAS_SPACE;
    const RDouble LITTLE     = 1.0E-10;
    int nEquation            = invSchemePara->GetNumberOfTotalEquation();
    RDouble *gamaL           = invSchemePara->GetLeftGama();
    RDouble *gamaR           = invSchemePara->GetRightGama();
    RDouble *pressureCoeffL  = face_proxy->GetPressureCoefficientL();
    RDouble *pressureCoeffR  = face_proxy->GetPressureCoefficientR();
    RDouble entrFixCoefLeft  = invSchemePara->GetLeftEntropyFixCoefficient();
    RDouble entrFixCoefRight = invSchemePara->GetRightEntropyFixCoefficient();
    RDouble *xfn             = invSchemePara->GetFaceNormalX();
    RDouble *yfn             = invSchemePara->GetFaceNormalY();
    RDouble *zfn             = invSchemePara->GetFaceNormalZ();
    RDouble *vgn             = invSchemePara->GetFaceVelocity();
    RDouble refMachNumber    = invSchemePara->GetMachNumber();
    int RoeEntropyFixMethod  = invSchemePara->GetEntropyFixMethod();
    RDouble mach2            = refMachNumber * refMachNumber;
    RDouble *dq              = new RDouble[nEquation]();

    int nPrecondition = invSchemePara->GetIfPrecondition(); 
    RDouble xsn       = xfn[iFace];
    RDouble ysn       = yfn[iFace];
    RDouble zsn       = zfn[iFace];
    RDouble vb        = vgn[iFace];
    RDouble gmL       = gamaL[iFace];
    RDouble gmR       = gamaR[iFace];

    RDouble& rL  = qCvL[IR];
    RDouble& ruL = qCvL[IU];
    RDouble& rvL = qCvL[IV];
    RDouble& rwL = qCvL[IW];
    RDouble& rEL = qCvL[IP];

    RDouble& rR  = qCvR[IR];
    RDouble& ruR = qCvR[IU];
    RDouble& rvR = qCvR[IV];
    RDouble& rwR = qCvR[IW];
    RDouble& rER = qCvR[IP];

    RDouble& uL = primL[IU];
    RDouble& vL = primL[IV];
    RDouble& wL = primL[IW];
    RDouble& pL = primL[IP];

    RDouble& uR = primR[IU];
    RDouble& vR = primR[IV];
    RDouble& wR = primR[IW];
    RDouble& pR = primR[IP];

    rL = MAX(rL, LITTLE);
    rR = MAX(rR, LITTLE);
    pL = MAX(pL, LITTLE);
    pR = MAX(pR, LITTLE);

    // RDouble vSquareL = uL * uL + vL * vL + wL * wL;
    // RDouble vSquareR = uR * uR + vR * vR + wR * wR;

    RDouble hint_L = (gmL / (gmL - one)) * (pL / rL);
    RDouble hint_R = (gmR / (gmR - one)) * (pR / rR);
    RDouble peL    = 0.0, peR = 0.0;

    //if (nChemical != 0 && nIdealState == 0)
    //{
    //    for (int i = 0; i < 3; ++i)
    //    {
    //        temperaturesL[i] = tl[i][iFace];
    //        temperaturesR[i] = tr[i][iFace];
    //    }
    //    gas->ComputeEnthalpyByPrimitive(primL, gmL, hint_L, temperaturesL);
    //    gas->ComputeEnthalpyByPrimitive(primR, gmR, hint_R, temperaturesR);
    //
    //    if (nTemperatureModel > 1)
    //    {
    //        peL = gas->GetElectronPressure(primL, temperaturesL[ITE]);
    //        peR = gas->GetElectronPressure(primR, temperaturesR[ITE]);
    //    }
    //}
    //RDouble hL = hint_L + half * vSquareL;
    //RDouble hR = hint_R + half * vSquareR;

    RDouble vnL = xsn * uL + ysn * vL + zsn * wL - vb;
    RDouble vnR = xsn * uR + ysn * vR + zsn * wR - vb;
    //RDouble rvnL = rL * vnL;
    //RDouble rvnR = rR * vnR;

    //! The first part of Roe scheme: F(QL, QR) = 0.5 * [ F(QL) + F(QR) ].
    RDouble rvnL = ruL * xsn + rvL * ysn + rwL * zsn - rL * vb;
    RDouble rvnR = ruR * xsn + rvR * ysn + rwR * zsn - rR * vb;
    RDouble hL = (rEL + pL) / rL;
    RDouble hR = (rER + pR) / rR;

    flux[IR]  = half * (rvnL + rvnR);
    flux[IRU] = half * (rvnL * uL + xfn[iFace] * pL + rvnR * uR + xfn[iFace] * pR);
    flux[IRV] = half * (rvnL * vL + yfn[iFace] * pL + rvnR * vR + yfn[iFace] * pR);
    flux[IRW] = half * (rvnL * wL + zfn[iFace] * pL + rvnR * wR + zfn[iFace] * pR);
    flux[IRE] = half * (rvnL * hL + vgn[iFace] * pL + rvnR * hR + vgn[iFace] * pR);

    //! The following loop is useless.
    //for (int iLaminar = nm; iLaminar < nLaminar; ++iLaminar)
    //{
    //    flux[iLaminar][iFace] = half * (primL[iLaminar] * rvnL + primR[iLaminar] * rvnR);
    //}
    //
    //for (int it = 1; it < nTemperatureModel; ++it)
    //{
    //    flux[nLaminar + it][iFace] = half * (primL[nLaminar + it] * rvnL + primR[nLaminar + it] * rvnR);
    //}
    //if (nChemical == 1)
    //{
    //    flux[nLaminar + nTemperatureModel - 1][iFace] += half * (peL * (vnL + vgn[iFace]) + peR * (vnR + vgn[iFace]));
    //}

    RDouble ratio = sqrt(rR / rL);
    RDouble roeAverageCoef = one / (one + ratio);

    RDouble rRoeAverage = sqrt(rL * rR);
    RDouble uRoeAverage = (uL + uR * ratio) * roeAverageCoef;
    RDouble vRoeAverage = (vL + vR * ratio) * roeAverageCoef;
    RDouble wRoeAverage = (wL + wR * ratio) * roeAverageCoef;
    RDouble pRoeAverage = (pL + pR * ratio) * roeAverageCoef;

    RDouble gama = (gmL + gmR * ratio) * roeAverageCoef;
    RDouble gamm1 = gama - one;

    RDouble vSquareRoeAverage = uRoeAverage * uRoeAverage + vRoeAverage * vRoeAverage + wRoeAverage * wRoeAverage;
    RDouble absVel = sqrt(vSquareRoeAverage);
    RDouble hRoeAverage = gama / gamm1 * pRoeAverage / rRoeAverage + half * vSquareRoeAverage;

    RDouble vn = xsn * uRoeAverage + ysn * vRoeAverage + zsn * wRoeAverage - vb;
    RDouble theta = vn + vb;

    //! c2: square of sound speed.
    RDouble c2 = gamm1 * (hRoeAverage - half * vSquareRoeAverage);

    //! Sound speed: gama * pRoeAverage / rRoeAverage.
    RDouble cm = sqrt(ABS(c2));

    RDouble Ur2 = 0, cnew2 = 0, cnew = 0, vnew = 0, beta = 0, alpha = 0;
    RDouble eigv1 = 0, eigv2 = 0, eigv3 = 0;
    eigv1 = ABS(vn);
    if (nPrecondition == 0)
    {
        eigv2 = ABS(vn + cm);
        eigv3 = ABS(vn - cm);
    }
    else
    {
        beta = MIN(MAX(vSquareRoeAverage / c2, three * mach2), one);
        vnew = half * vn * (one + beta);
        cnew = half * sqrt(((beta - one) * vn) * ((beta - one) * vn) + four * beta * c2);
        cnew2 = cnew * cnew;
        eigv2 = ABS(vnew + cnew);
        eigv3 = ABS(vnew - cnew);
    }

    //RDouble EL = (one / (gmL - one)) * (pL / rL) + half * vSquareL;
    //RDouble ER = (one / (gmR - one)) * (pR / rR) + half * vSquareR;

    //! The jump between left and right conservative variables.
    dq[IR]  = rR - rL;
    dq[IRU] = ruR - ruL;
    dq[IRV] = rvR - rvL;
    dq[IRW] = rwR - rwL;
    dq[IRE] = rER - rEL;

    //if (nChemical != 0)
    //{
    //    if (nTemperatureModel > 1)
    //    {
    //        gas->GetTemperature(primL, temperaturesL[ITT], temperaturesL[ITV], temperaturesL[ITE]);
    //        gas->GetTemperature(primR, temperaturesR[ITT], temperaturesR[ITV], temperaturesR[ITE]);
    //    }
    //    gas->Primitive2Conservative(primL, gmL, temperaturesL[ITV], temperaturesL[ITE], qConserveL);
    //    gas->Primitive2Conservative(primR, gmR, temperaturesR[ITV], temperaturesR[ITE], qConserveR);
    //}
    //
    //for (int iLaminar = nm; iLaminar < nEquation; ++iLaminar)
    //{
    //    dq[iLaminar] = qConserveR[iLaminar] - qConserveL[iLaminar];
    //}

    //! Entropy correction, to limit the magnitude of the three eigenvalue.
    //! Warning: This function calling may affect the efficient!

    Flux_RoeEntropyCorrection(eigv1, eigv2, eigv3, pressureCoeffL[iFace], pressureCoeffR[iFace], absVel, vn, cm, RoeEntropyFixMethod, entrFixCoefLeft, entrFixCoefRight);

    RDouble xi1 = 0.0, xi2 = 0.0, xi3 = 0.0, xi4 = 0.0;
    if (nPrecondition == 0)
    {
        xi1 = (two * eigv1 - eigv2 - eigv3) / (two * c2);
        xi2 = (eigv2 - eigv3) / (two * cm);
    }
    else
    {
        xi1 = (two * eigv1 - eigv2 - eigv3) / (two * cnew2);
        xi2 = (eigv2 - eigv3) / (two * cnew);
        alpha = half * (one - beta);
        xi3 = two * alpha * eigv1 / cnew2;
        xi4 = alpha * vn * xi2 / cnew2;
    }

    RDouble dc = theta * dq[IR] - xsn * dq[IRU] - ysn * dq[IRV] - zsn * dq[IRW];

    RDouble c2dc = 0;
    if (nPrecondition == 0)
    {
        c2dc = c2 * dc;
    }
    else
    {
        Ur2 = beta * c2;
        c2dc = cnew2 * dc;
    }

    RDouble dh = -gamm1 * (uRoeAverage * dq[IRU] + vRoeAverage * dq[IRV] + wRoeAverage * dq[IRW] - dq[IRE]) + half * gamm1 * vSquareRoeAverage * dq[IR];

    //if (nChemical != 0)
    //{
    //    gas->ComputeDHAndTotalEnthalpy(prim, gama, dq, dh, hRoeAverage);
    //}

    if (nPrecondition == 0)
    {
        flux[IR]     -= half * (eigv1 * dq[IR] - dh * xi1 - dc * xi2);
        flux[IRU]    -= half * (eigv1 * dq[IRU] + (xsn * c2dc - uRoeAverage * dh) * xi1 + (xsn * dh - uRoeAverage * dc) * xi2);
        flux[IRV]    -= half * (eigv1 * dq[IRV] + (ysn * c2dc - vRoeAverage * dh) * xi1 + (ysn * dh - vRoeAverage * dc) * xi2);
        flux[IRW]    -= half * (eigv1 * dq[IRW] + (zsn * c2dc - wRoeAverage * dh) * xi1 + (zsn * dh - wRoeAverage * dc) * xi2);
        flux[IRE]    -= half * (eigv1 * dq[IRE] + (theta * c2dc - hRoeAverage * dh) * xi1 + (theta * dh - hRoeAverage * dc) * xi2);

        //for (int iLaminar = nm; iLaminar < nLaminar; ++iLaminar)
        //{
        //    prim[iLaminar] = half * (primL[iLaminar] + primR[iLaminar]);
        //    flux[iLaminar][iFace] -= half * (eigv1 * dq[iLaminar] - prim[iLaminar] * (dh * xi1 + dc * xi2));
        //}
    }
    else
    {
        flux[IR]     -= half * (eigv1 * dq[IR] - dh * cnew2 / Ur2 * xi1 - dc * xi2 + dh * cnew2 / Ur2 * xi3 - dh * cnew2 / Ur2 * xi4);
        flux[IRU]    -= half * (eigv1 * dq[IRU] + (xsn * c2dc - uRoeAverage * dh * cnew2 / Ur2) * xi1 + (xsn * dh - uRoeAverage * dc) * xi2 + uRoeAverage * dh * cnew2 / Ur2 * xi3 - (xsn * c2dc + uRoeAverage * dh * cnew2 / Ur2) * xi4);
        flux[IRV]    -= half * (eigv1 * dq[IRV] + (ysn * c2dc - vRoeAverage * dh * cnew2 / Ur2) * xi1 + (ysn * dh - vRoeAverage * dc) * xi2 + vRoeAverage * dh * cnew2 / Ur2 * xi3 - (ysn * c2dc + vRoeAverage * dh * cnew2 / Ur2) * xi4);
        flux[IRW]    -= half * (eigv1 * dq[IRW] + (zsn * c2dc - wRoeAverage * dh * cnew2 / Ur2) * xi1 + (zsn * dh - wRoeAverage * dc) * xi2 + wRoeAverage * dh * cnew2 / Ur2 * xi3 - (zsn * c2dc + wRoeAverage * dh * cnew2 / Ur2) * xi4);
        flux[IRE]    -= half * (eigv1 * dq[IRE] + (theta * c2dc - hRoeAverage * dh * cnew2 / Ur2) * xi1 + (theta * dh - hRoeAverage * dc) * xi2 + hRoeAverage * dh * cnew2 / Ur2 * xi3 - (theta * c2dc + hRoeAverage * dh * cnew2 / Ur2) * xi4);

        //for (int iLaminar = nm; iLaminar < nLaminar; ++iLaminar)
        //{
        //    prim[iLaminar] = half * (primL[iLaminar] + primR[iLaminar]);
        //    flux[iLaminar][iFace] -= half * (eigv1 * dq[iLaminar] - prim[iLaminar] * (dh * cnew2 / Ur2 * (xi1 - xi3 + xi4) + dc * xi2));
        //}
        ////! Compute vibration-electron energy term.
        //for (int it = 1; it < nTemperatureModel; ++it)
        //{
        //    prim[nLaminar + it] = half * (primL[nLaminar + it] + primR[nLaminar + it]);
        //    flux[nLaminar + it][iFace] -= half * (eigv1 * dq[nLaminar + it] - prim[nLaminar + it] * (dh * cnew2 / Ur2 * (xi1 - xi3 + xi4) + dc * xi2));
        //}
    }

    delete [] dq;    dq = nullptr;
}
#endif

void Roe_Scheme_PrimitiveForm(FaceProxy *faceProxy, InviscidSchemeParameter *inviscidSchemePara)
{
    const RDouble LITTLE = 1.0E-10;
    int nTemperatureModel = inviscidSchemePara->GetNumberOfTemperatureModel();
    int nEquation = inviscidSchemePara->GetNumberOfTotalEquation();
    int nLaminar = inviscidSchemePara->GetNumberOfLaminar();
    int nLength = inviscidSchemePara->GetLength();
    int nm = inviscidSchemePara->Getnm();
    int RoeEntropyFixMethod = inviscidSchemePara->GetEntropyFixMethod();
    int nChemical = inviscidSchemePara->GetNumberOfChemical();
    int nPrecondition = inviscidSchemePara->GetIfPrecondition();

    RDouble refMachNumber = inviscidSchemePara->GetMachNumber();
    RDouble entrFixCoefLeft = inviscidSchemePara->GetLeftEntropyFixCoefficient();
    RDouble entrFixCoefRight = inviscidSchemePara->GetRightEntropyFixCoefficient();
    RDouble **qL = inviscidSchemePara->GetLeftQ();
    RDouble **qR = inviscidSchemePara->GetRightQ();
    RDouble *gamaL = inviscidSchemePara->GetLeftGama();
    RDouble *gamaR = inviscidSchemePara->GetRightGama();

    RDouble *xfn   = inviscidSchemePara->GetFaceNormalX();
    RDouble *yfn   = inviscidSchemePara->GetFaceNormalY();
    RDouble *zfn   = inviscidSchemePara->GetFaceNormalZ();
    RDouble *vgn   = inviscidSchemePara->GetFaceVelocity();
    RDouble **flux = inviscidSchemePara->GetFlux();

    RDouble *pressureCoeffL = faceProxy->GetPressureCoefficientL();
    RDouble *pressureCoeffR = faceProxy->GetPressureCoefficientR();
    RDouble **tl = inviscidSchemePara->GetLeftTemperature();
    RDouble **tr = inviscidSchemePara->GetRightTemperature();
    RDouble mach2 = refMachNumber * refMachNumber;

    using namespace IDX;
    using namespace GAS_SPACE;

    RDouble *dq    = new RDouble[nEquation];

    RDouble *primL = new RDouble[nEquation];
    RDouble *primR = new RDouble[nEquation];
    RDouble temperaturesL[3], temperaturesR[3];

    for (int iFace = 0; iFace < nLength; ++ iFace)
    {
        for (int iEqn = 0; iEqn < nEquation; ++ iEqn)
        {
            primL[iEqn] = qL[iEqn][iFace];
            primR[iEqn] = qR[iEqn][iFace];
        }

        RDouble &rL = primL[IR];
        RDouble &uL = primL[IU];
        RDouble &vL = primL[IV];
        RDouble &wL = primL[IW];
        RDouble &pL = primL[IP];

        RDouble &rR = primR[IR];
        RDouble &uR = primR[IU];
        RDouble &vR = primR[IV];
        RDouble &wR = primR[IW];
        RDouble &pR = primR[IP];

        RDouble xsn = xfn[iFace];
        RDouble ysn = yfn[iFace];
        RDouble zsn = zfn[iFace];
        RDouble vb  = vgn[iFace];
        RDouble gmL = gamaL[iFace];
        RDouble gmR = gamaR[iFace];

        rL = MAX(rL, LITTLE);
        rR = MAX(rR, LITTLE);
        pL = MAX(pL, LITTLE);
        pR = MAX(pR, LITTLE);

        RDouble vSquareL = uL * uL + vL * vL + wL * wL;
        RDouble vSquareR = uR * uR + vR * vR + wR * wR;

        RDouble hint_L = (gmL / (gmL - one)) * (pL / rL);
        RDouble hint_R = (gmR / (gmR - one)) * (pR / rR);
        RDouble peL = 0.0, peR = 0.0;
        if (nChemical != 0)
        {
            for (int i = 0; i < 3; ++ i)
            {
                temperaturesL[i] = tl[i][iFace];
                temperaturesR[i] = tr[i][iFace];
            }
            gas->ComputeEnthalpyByPrimitive(primL, gmL, hint_L, temperaturesL);
            gas->ComputeEnthalpyByPrimitive(primR, gmR, hint_R, temperaturesR);

            if (nTemperatureModel > 1)
            {
                peL = gas->GetElectronPressure(primL, temperaturesL[ITE]);
                peR = gas->GetElectronPressure(primR, temperaturesR[ITE]);
            }
        }

        //! H = h + 0.5 * V^2.
        //! h = gam / (gam - 1) * p / rho.
        RDouble HL = hint_L + half * vSquareL;
        RDouble HR = hint_R + half * vSquareR;

        RDouble vnL  = xsn * uL  + ysn * vL  + zsn * wL  - vb;
        RDouble vnR  = xsn * uR  + ysn * vR  + zsn * wR  - vb;
        RDouble rvnL = rL * vnL;
        RDouble rvnR = rR * vnR;

        //! The first part of Roe scheme: F(QL, QR) = 0.5 * [ F(QL) + F(QR) ].
        flux[IR ][iFace] =  half * (rvnL                        + rvnR                     );
        flux[IRU][iFace] =  half * (rvnL * uL + xfn[iFace] * pL + rvnR * uR + xfn[iFace] * pR);
        flux[IRV][iFace] =  half * (rvnL * vL + yfn[iFace] * pL + rvnR * vR + yfn[iFace] * pR);
        flux[IRW][iFace] =  half * (rvnL * wL + zfn[iFace] * pL + rvnR * wR + zfn[iFace] * pR);
        flux[IRE][iFace] =  half * (rvnL * HL + vgn[iFace] * pL + rvnR * HR + vgn[iFace] * pR);

        //! The following loop is useless.
        for (int iLaminar = nm; iLaminar < nLaminar; ++ iLaminar)
        {
            flux[iLaminar][iFace] = half * (primL[iLaminar] * rvnL + primR[iLaminar] * rvnR);
        }
        for (int it = 1; it < nTemperatureModel; ++ it)
        {
            flux[nLaminar + it][iFace] = half * (primL[nLaminar + it] * rvnL + primR[nLaminar + it] * rvnR);
        }
        if (nChemical == 1)
        {
            flux[nLaminar + nTemperatureModel - 1][iFace] += half * (peL * (vnL + vgn[iFace]) + peR * (vnR + vgn[iFace]));
        }

        //! Compute diffusion terms following.

        //! Roe-average firstly. Using density, three velocity components and pressure
        RDouble ratio = sqrt(rR / rL);
        RDouble roeAverageCoef = one / (one + ratio);

        RDouble rRoeAverage = sqrt(rL * rR);
        RDouble uRoeAverage = (uL + uR * ratio) * roeAverageCoef;
        RDouble vRoeAverage = (vL + vR * ratio) * roeAverageCoef;
        RDouble wRoeAverage = (wL + wR * ratio) * roeAverageCoef;
        RDouble HRoeAverage = (HL + HR * ratio) * roeAverageCoef;

        RDouble gama = (gmL + gmR * ratio) * roeAverageCoef;
        RDouble gamm1 = gama - one;

        RDouble vSquareRoeAverage = uRoeAverage * uRoeAverage + vRoeAverage * vRoeAverage + wRoeAverage * wRoeAverage;
        RDouble absVel = sqrt(vSquareRoeAverage);
        RDouble vn = xsn * uRoeAverage + ysn * vRoeAverage + zsn * wRoeAverage - vb;

        //! c2: square of sound speed.
        RDouble c2 = gamm1 * (HRoeAverage - half * vSquareRoeAverage);

        //! Sound speed: gama * pRoeAverage / rRoeAverage.
        RDouble cm = sqrt(ABS(c2));

        RDouble prconCoeff;
        RDouble eigv1, eigv2, eigv3;
        eigv1 = ABS(vn);
        if (nPrecondition == 1)
        {
            prconCoeff = MIN(MAX(vSquareRoeAverage / c2, three * mach2), 1.0);
            vn = half * vn * (1.0 + prconCoeff);
            cm = half * sqrt(((1.0 - prconCoeff) * vn) * ((1.0 - prconCoeff) * vn) + 4.0 * prconCoeff * c2);
            c2 = cm * cm;
        }
        eigv2 = ABS(vn + cm);
        eigv3 = ABS(vn - cm);

        //! The jump between left and right primitive variables.
        dq[IR] = rR - rL;
        dq[IU] = uR - uL;
        dq[IV] = vR - vL;
        dq[IW] = wR - wL;
        dq[IP] = pR - pL;
        RDouble dVn = vnR - vnL;

        if (nChemical != 0)
        {
            TK_Exit::UnexpectedVarValue("Do not consider chemical = ", 1);
        }

        //! Entropy correction, to limit the magnitude of the three eigenvalue.
        //! Warning: This function calling may affect the efficient!

        Flux_RoeEntropyCorrection(eigv1, eigv2, eigv3, pressureCoeffL[iFace], pressureCoeffR[iFace], absVel, vn, cm, RoeEntropyFixMethod, entrFixCoefLeft, entrFixCoefRight);

        RDouble eigv1Term = 0.0, eigv2Term = 0.0, eigv3Term = 0.0;
        RDouble oSoundSpeed2 = 1.0 / (cm * cm);
        if (nPrecondition == 0) 
        {
            eigv1Term = dq[IR] - dq[IP] * oSoundSpeed2;
            eigv2Term = half * (dq[IP] + rRoeAverage * cm * dVn) * oSoundSpeed2;
            eigv3Term = half * (dq[IP] - rRoeAverage * cm * dVn) * oSoundSpeed2;
        }
        else
        {
            TK_Exit::UnexpectedVarValue("Do not consider nPrecondition = ", 1);
        }

        if (nChemical != 0)
        {
            TK_Exit::UnexpectedVarValue("Do not consider chemical = ", 1);
            //ComputeDHAndTotalEnthalpy(prim, gama, dq, dh, hRoeAverage);
        }

        RDouble fluxTerm111 = 1.0;
        RDouble fluxTerm112 = uRoeAverage;
        RDouble fluxTerm113 = vRoeAverage;
        RDouble fluxTerm114 = wRoeAverage;
        RDouble fluxTerm115 = half * vSquareRoeAverage;

        RDouble fluxTerm121 = 0.0;
        RDouble fluxTerm122 = dq[IU] - dVn * xsn;
        RDouble fluxTerm123 = dq[IV] - dVn * ysn;
        RDouble fluxTerm124 = dq[IW] - dVn * zsn;
        RDouble fluxTerm125 = uRoeAverage * dq[IU] + vRoeAverage * dq[IV] + wRoeAverage * dq[IW] - vn * dVn;

        RDouble fluxTerm21 = 1.0;
        RDouble fluxTerm22 = uRoeAverage + cm * xsn;
        RDouble fluxTerm23 = vRoeAverage + cm * ysn;
        RDouble fluxTerm24 = wRoeAverage + cm * zsn;
        RDouble fluxTerm25 = HRoeAverage + cm * vn;

        RDouble fluxTerm31 = 1.0;
        RDouble fluxTerm32 = uRoeAverage - cm * xsn;
        RDouble fluxTerm33 = vRoeAverage - cm * ysn;
        RDouble fluxTerm34 = wRoeAverage - cm * zsn;
        RDouble fluxTerm35 = HRoeAverage - cm * vn;

        if (nPrecondition == 0)
        {
            flux[IR ][iFace] -= half * (eigv1 * eigv1Term * fluxTerm111 + eigv1 * rRoeAverage * fluxTerm121 + eigv2 * eigv2Term * fluxTerm21 + eigv3 * eigv3Term * fluxTerm31);
            flux[IRU][iFace] -= half * (eigv1 * eigv1Term * fluxTerm112 + eigv1 * rRoeAverage * fluxTerm122 + eigv2 * eigv2Term * fluxTerm22 + eigv3 * eigv3Term * fluxTerm32);
            flux[IRV][iFace] -= half * (eigv1 * eigv1Term * fluxTerm113 + eigv1 * rRoeAverage * fluxTerm123 + eigv2 * eigv2Term * fluxTerm23 + eigv3 * eigv3Term * fluxTerm33);
            flux[IRW][iFace] -= half * (eigv1 * eigv1Term * fluxTerm114 + eigv1 * rRoeAverage * fluxTerm124 + eigv2 * eigv2Term * fluxTerm24 + eigv3 * eigv3Term * fluxTerm34);
            flux[IRE][iFace] -= half * (eigv1 * eigv1Term * fluxTerm115 + eigv1 * rRoeAverage * fluxTerm125 + eigv2 * eigv2Term * fluxTerm25 + eigv3 * eigv3Term * fluxTerm35);
        }
        else
        {
            TK_Exit::UnexpectedVarValue("Do not consider nPrecondition = ", 1);
        }

    }

    delete [] primL;    primL = nullptr;
    delete [] primR;    primR = nullptr;
    delete [] dq;    dq = nullptr;
}

//! clz begin
void Steger_Scheme(FaceProxy *face_proxy, InviscidSchemeParameter *invSchemePara)
{
    int nTemperatureModel = invSchemePara->GetNumberOfTemperatureModel();
    int nEquation = invSchemePara->GetNumberOfTotalEquation();
    int nLaminar = invSchemePara->GetNumberOfLaminar();
    int nLen = invSchemePara->GetLength();
    int nm = invSchemePara->Getnm();
    int EntropyFixMethod = invSchemePara->GetEntropyFixMethod();
    int nChemical = invSchemePara->GetNumberOfChemical();
    int nIdealState = invSchemePara->GetIfIdealState();

    int nElectronIndex = invSchemePara->GetIndexOfElectron();

    RDouble entrFixCoefLeft  = invSchemePara->GetLeftEntropyFixCoefficient();
    RDouble entrFixCoefRight = invSchemePara->GetRightEntropyFixCoefficient();
    RDouble **ql = invSchemePara->GetLeftQ();
    RDouble **qr = invSchemePara->GetRightQ();
    RDouble *gamal = invSchemePara->GetLeftGama();
    RDouble *gamar = invSchemePara->GetRightGama();
    RDouble *xfn   = invSchemePara->GetFaceNormalX();
    RDouble *yfn   = invSchemePara->GetFaceNormalY();
    RDouble *zfn   = invSchemePara->GetFaceNormalZ();
    RDouble *vgn   = invSchemePara->GetFaceVelocity();
    RDouble **flux = invSchemePara->GetFlux();
    RDouble **tl = invSchemePara->GetLeftTemperature();
    RDouble **tr = invSchemePara->GetRightTemperature();

    using namespace IDX;
    using namespace GAS_SPACE;

    RDouble *fl    = new RDouble[nEquation];
    RDouble *fr    = new RDouble[nEquation];
    RDouble *priml = new RDouble[nEquation];
    RDouble *primr = new RDouble[nEquation];
    RDouble temperaturesL[3], temperaturesR[3];

    for (int iLen = 0; iLen < nLen; ++ iLen)
    {
        for (int iEqn = 0; iEqn < nEquation; ++ iEqn)
        {
            priml[iEqn] = ql[iEqn][iLen];
            primr[iEqn] = qr[iEqn][iLen];
        }

        RDouble rl = priml[IR];
        RDouble ul = priml[IU];
        RDouble vl = priml[IV];
        RDouble wl = priml[IW];
        RDouble pl = priml[IP];

        RDouble rr = primr[IR];
        RDouble ur = primr[IU];
        RDouble vr = primr[IV];
        RDouble wr = primr[IW];
        RDouble pr = primr[IP];

        RDouble xsn = xfn[iLen];
        RDouble ysn = yfn[iLen];
        RDouble zsn = zfn[iLen];
        RDouble vb  = vgn[iLen];

        RDouble gmL = gamal[iLen];
        RDouble gmR = gamar[iLen];

        RDouble v2l = ul * ul + vl * vl + wl * wl;
        RDouble v2r = ur * ur + vr * vr + wr * wr;

        RDouble hintLeft  = (gmL / (gmL - one)) * (pl / rl);
        RDouble hintRight = (gmR / (gmR - one)) * (pr / rr);

        RDouble pel = 0.0, per = 0.0;

        if (nChemical != 0 && nIdealState == 0)
        {
            for (int i = 0; i < 3; ++ i)
            {
                temperaturesL[i] = tl[i][iLen];
                temperaturesR[i] = tr[i][iLen];
            }
            gas->ComputeEnthalpyByPrimitive(priml, gmL, hintLeft, temperaturesL);
            gas->ComputeEnthalpyByPrimitive(primr, gmR, hintRight, temperaturesR);

            if (nTemperatureModel > 1)
            {
                pel = gas->GetElectronPressure(priml, temperaturesL[ITE]);
                per = gas->GetElectronPressure(primr, temperaturesR[ITE]);
            }
        }

        RDouble hl = hintLeft  + half * v2l;
        RDouble hr = hintRight + half * v2r;
        RDouble el = hl - pl / (rl + SMALL);
        RDouble er = hr - pr / (rr + SMALL);

        RDouble vnl = xsn * ul + ysn * vl + zsn * wl - vb;
        RDouble vnr = xsn * ur + ysn * vr + zsn * wr - vb;

        RDouble c2l = gmL * pl / rl;
        RDouble cl  = sqrt(ABS(c2l));

        RDouble c2r = gmR * pr / rr;
        RDouble cr  = sqrt(ABS(c2r));

        //! The eigenvalues of left and right.
        RDouble eigvl1 = vnl;
        RDouble eigvl2 = vnl + cl;
        RDouble eigvl3 = vnl - cl;

        RDouble eigvr1 = vnr;
        RDouble eigvr2 = vnr + cr;
        RDouble eigvr3 = vnr - cr;

        RDouble absEigvl1 = ABS(eigvl1);
        RDouble absEigvl2 = ABS(eigvl2);
        RDouble absEigvl3 = ABS(eigvl3);

        RDouble absEigvr1 = ABS(eigvr1);
        RDouble absEigvr2 = ABS(eigvr2);
        RDouble absEigvr3 = ABS(eigvr3);

        //! Entropy fix.
        //! Modified by clz : 2012-8-15.
        if (EntropyFixMethod == 1)
        {
            RDouble maxEigv = MAX(absEigvl2, absEigvl3);
            RDouble tmp0    = maxEigv * entrFixCoefLeft ;
            RDouble tmp1    = maxEigv * entrFixCoefRight;
            absEigvl1       = MAX(tmp0, absEigvl1);
            absEigvl2       = MAX(tmp1, absEigvl2);
            absEigvl3       = MAX(tmp1, absEigvl3);

            maxEigv   = MAX(absEigvr2, absEigvr3);
            tmp0      = maxEigv * entrFixCoefLeft ;
            tmp1      = maxEigv * entrFixCoefRight;
            absEigvr1 = MAX(tmp0, absEigvr1);
            absEigvr2 = MAX(tmp1, absEigvr2);
            absEigvr3 = MAX(tmp1, absEigvr3);
        }
        else if (EntropyFixMethod == 2)
        {
            //RDouble tmp0 = MAX(entrFixCoefLeft , static_cast <RDouble> (0.02));
            //RDouble tmp1 = MAX(entrFixCoefRight, static_cast <RDouble> (0.02));

            //absEigvl1 = MAX(tmp0, absEigvl1);
            //absEigvl2 = MAX(tmp1, absEigvl2);
            //absEigvl3 = MAX(tmp1, absEigvl3);
                   
            //absEigvr1 = MAX(tmp0, absEigvr1);
            //absEigvr2 = MAX(tmp1, absEigvr2);
            //absEigvr3 = MAX(tmp1, absEigvr3);

            absEigvl1 = sqrt(absEigvl1*absEigvl1 + 1.E-4);
            absEigvl2 = sqrt(absEigvl2*absEigvl2 + 1.E-4);
            absEigvl3 = sqrt(absEigvl3*absEigvl3 + 1.E-4);

            absEigvr1 = sqrt(absEigvr1*absEigvr1 + 1.E-4);
            absEigvr2 = sqrt(absEigvr2*absEigvr2 + 1.E-4);
            absEigvr3 = sqrt(absEigvr3*absEigvr3 + 1.E-4);
        }
        eigvl1 = half * (eigvl1 + absEigvl1);
        eigvl2 = half * (eigvl2 + absEigvl2);
        eigvl3 = half * (eigvl3 + absEigvl3);

        eigvr1 = half * (eigvr1 - absEigvr1);
        eigvr2 = half * (eigvr2 - absEigvr2);
        eigvr3 = half * (eigvr3 - absEigvr3);

        RDouble c2gaml  = c2l / gmL;
        RDouble c2gamr  = c2r / gmR;

        RDouble xil1 = c2gaml * (eigvl1 + eigvl1 - eigvl2 - eigvl3) / (c2l + c2l);
        RDouble xil2 = c2gaml * (eigvl2 - eigvl3) / (cl + cl);

        RDouble xir1 = c2gamr * (eigvr1 + eigvr1 - eigvr2 - eigvr3) / (c2r + c2r);
        RDouble xir2 = c2gamr * (eigvr2 - eigvr3) / (cr + cr);

        RDouble eig_xi = eigvl1 - xil1;

        fl[IR ] = rl * (eig_xi                );
        fl[IRU] = rl * (eig_xi * ul + xil2 * xsn);
        fl[IRV] = rl * (eig_xi * vl + xil2 * ysn);
        fl[IRW] = rl * (eig_xi * wl + xil2 * zsn);
        fl[IRE] = rl * (eigvl1 * el - xil1 * hl + xil2 * (vnl + vb));
        RDouble fel = - xil1 * pel;

        eig_xi = eigvr1 - xir1;

        fr[IR ] = rr * (eig_xi                );
        fr[IRU] = rr * (eig_xi * ur + xir2 * xsn);
        fr[IRV] = rr * (eig_xi * vr + xir2 * ysn);
        fr[IRW] = rr * (eig_xi * wr + xir2 * zsn);
        fr[IRE] = rr * (eigvr1 * er - xir1 * hr + xir2 * (vnr + vb));
        RDouble fer = - xir1 * per;

        for (int iLaminar = nm; iLaminar < nEquation; ++ iLaminar)
        {
            fl[iLaminar] = priml[iLaminar] * fl[IR];
            fr[iLaminar] = primr[iLaminar] * fr[IR];
        }
        if (nChemical == 1)
        {
            fl[nLaminar + nTemperatureModel - 1] += fel;
            fr[nLaminar + nTemperatureModel - 1] += fer;
        }

        //! Compute vibration-electron energy term.
        for (int iLaminar = 0; iLaminar < nEquation; ++ iLaminar)
        {
            flux[iLaminar][iLen] = fl[iLaminar] + fr[iLaminar];
        }

        //! Modify the fluxes of the species whose transport equations need not to be computed.
        if (nChemical == 1)
        {
            flux[nLaminar][iLen] = 0.0;
            if (nElectronIndex >= 0)
            {
                flux[nLaminar - 1][iLen] = 0.0;
            }
        }
    }

    delete [] fl;    fl = nullptr;
    delete [] fr;    fr = nullptr;
    delete [] priml;    priml = nullptr;
    delete [] primr;    primr = nullptr;
}



void Vanleer_Scheme(FaceProxy *face_proxy, InviscidSchemeParameter *invSchemePara)
{
    int nTemperatureModel = invSchemePara->GetNumberOfTemperatureModel();
    int nEquation = invSchemePara->GetNumberOfTotalEquation();
    int nLaminar = invSchemePara->GetNumberOfLaminar();
    int nLen = invSchemePara->GetLength();
    int nm = invSchemePara->Getnm();
    int nChemical = invSchemePara->GetNumberOfChemical();
    int nPrecondition = invSchemePara->GetIfPrecondition();
    int nIdealState = invSchemePara->GetIfIdealState();
    int nElectronIndex = invSchemePara->GetIndexOfElectron();
    RDouble kPrec = invSchemePara->GetPreconCoefficient();

    RDouble **qL = invSchemePara->GetLeftQ();
    RDouble **qR = invSchemePara->GetRightQ();
    RDouble *gamal = invSchemePara->GetLeftGama();
    RDouble *gamar = invSchemePara->GetRightGama();
    RDouble *xfn   = invSchemePara->GetFaceNormalX();
    RDouble *yfn   = invSchemePara->GetFaceNormalY();
    RDouble *zfn   = invSchemePara->GetFaceNormalZ();
    RDouble *vgn   = invSchemePara->GetFaceVelocity();
    RDouble **flux = invSchemePara->GetFlux();
    RDouble refMachNumber = invSchemePara->GetMachNumber();
    RDouble mach2 = refMachNumber * refMachNumber;
    RDouble **tl = invSchemePara->GetLeftTemperature();
    RDouble **tr = invSchemePara->GetRightTemperature();
 
    using namespace IDX;
    using namespace GAS_SPACE;

    RDouble *fluxL = new RDouble[nEquation];
    RDouble *fluxR = new RDouble[nEquation];
    RDouble *primL = new RDouble[nEquation];
    RDouble *primR = new RDouble[nEquation];
    RDouble temperaturesL[3], temperaturesR[3];

    for (int iLen = 0; iLen < nLen; ++ iLen)
    {
        for (int iEqn = 0; iEqn < nEquation; ++ iEqn)
        {
            primL[iEqn] = qL[iEqn][iLen];
            primR[iEqn] = qR[iEqn][iLen];
        }

        RDouble &rL = primL[IR];
        RDouble &uL = primL[IU];
        RDouble &vL = primL[IV];
        RDouble &wL = primL[IW];
        RDouble &pL = primL[IP];

        RDouble &rR = primR[IR];
        RDouble &uR = primR[IU];
        RDouble &vR = primR[IV];
        RDouble &wR = primR[IW];
        RDouble &pR = primR[IP];

        //! It is convenient for later package flux solver,and adding variables.
        RDouble xsn = xfn[iLen];
        RDouble ysn = yfn[iLen];
        RDouble zsn = zfn[iLen];
        RDouble vb =  vgn[iLen];
        RDouble gmL = gamal[iLen];
        RDouble gmR = gamar[iLen];

        RDouble v2L = uL * uL + vL * vL + wL * wL;
        RDouble v2R = uR * uR + vR * vR + wR * wR;

        RDouble hintL = (gmL / (gmL - one)) * (pL / rL);
        RDouble hintR = (gmR / (gmR - one)) * (pR / rR);

        RDouble peL = 0.0, peR = 0.0;

        if (nChemical != 0 && nIdealState == 0)
        {
            for (int i = 0; i < 3; ++ i)
            {
                temperaturesL[i] = tl[i][iLen];
                temperaturesR[i] = tr[i][iLen];
            }
            gas->ComputeEnthalpyByPrimitive(primL, gmL, hintL, temperaturesL);
            gas->ComputeEnthalpyByPrimitive(primR, gmR, hintR, temperaturesR);

            if (nTemperatureModel > 1)
            {
                peL = gas->GetElectronPressure(primL, tl[2][iLen]);
                peR = gas->GetElectronPressure(primR, tr[2][iLen]);
            }
        }

        RDouble HL  = hintL + half * v2L;
        RDouble HR  = hintR + half * v2R;

        //! Full flux.
        RDouble vnL  = xsn * uL + ysn * vL + zsn * wL - vb;
        RDouble vnR  = xsn * uR + ysn * vR + zsn * wR - vb;
        RDouble rhoVnL = rL * vnL;
        RDouble rhoVnR = rR * vnR;

        RDouble c2l   = gmL * pL / rL;
        RDouble cL    = sqrt(ABS(c2l));
        RDouble machL = vnL / cL;
        RDouble mach2L = machL * machL;
        RDouble Mr2l =sqrt(MIN(MAX(mach2L, kPrec * mach2), 1.0));

        RDouble c2r   = gmR * pR / rR;
        RDouble cR    = sqrt(ABS(c2r));
        RDouble machR = vnR / cR;
        RDouble mach2R = machR * machR;
        RDouble Mr2r =sqrt(MIN(MAX(mach2R, kPrec * mach2), 1.0));

        RDouble fel = 0.0, fer = 0.0;
        if (machL > one)
        {
            //! Supersonic flow from left to right, upwind by left.
            fluxL[IR ] = rhoVnL;
            fluxL[IRU] = rhoVnL * uL + xsn * pL;
            fluxL[IRV] = rhoVnL * vL + ysn * pL;
            fluxL[IRW] = rhoVnL * wL + zsn * pL;
            fluxL[IRE] = rhoVnL * HL + vb  * pL;
            fel = peL * (vnL + vb);
        }
        else if (machL < - one)
        {
            //! Supersonic flow from right to left, upwind by right.
            fluxL[IR ] = zero;
            fluxL[IRU] = zero;
            fluxL[IRV] = zero;
            fluxL[IRW] = zero;
            fluxL[IRE] = zero;
            fel = zero;
        }
        else
        {
            //! Subsonic flow.
            if (nPrecondition == 0)
            {
                //! Without precondition.
                RDouble fmass = 0.25 * rL * cL * (machL + 1.0) * (machL + 1.0);
                RDouble tmp = (- vnL + 2.0 * cL) / gmL;
                fluxL[IR ] = fmass;
                fluxL[IRU] = fmass * (uL + xsn * tmp);
                fluxL[IRV] = fmass * (vL + ysn * tmp);
                fluxL[IRW] = fmass * (wL + zsn * tmp);
                fluxL[IRE] = fmass * (HL + vb  * tmp);
            }
            else
            {
                //! With precondition.
                RDouble fmass = 0.25 * rL * cL * (machL * machL + 1.0);
                RDouble tmp0  = - vnL / gmL;
                RDouble tmp1  = rL * cL / gmL;
                fluxL[IR ] = half * rhoVnL + Mr2l *  fmass;
                fluxL[IRU] = half * (rhoVnL * uL + xsn * pL) + Mr2l * (fmass * (uL + xsn * tmp0) + uL * tmp1);
                fluxL[IRV] = half * (rhoVnL * vL + ysn * pL) + Mr2l * (fmass * (vL + ysn * tmp0) + vL * tmp1);
                fluxL[IRW] = half * (rhoVnL * wL + zsn * pL) + Mr2l * (fmass * (wL + zsn * tmp0) + wL * tmp1);
                fluxL[IRE] = half * (rhoVnL * HL + vb  * pL) + Mr2l * (fmass * (HL + vb  * tmp0));
            }
        }

        //! Other equations besides N-S 5-equations, such as: multi-species, electronic.
        //! which are stored start from iEqua=5, 6, ...
        for (int iLaminar = nm; iLaminar < nLaminar; ++ iLaminar)
        {
            fluxL[iLaminar] = primL[iLaminar] * fluxL[IR];
        }

        //! Right contribution.
        if (machR < - one)
        {
            //! Supersonic flow from right to left, upwind by right.
            fluxR[IR ] = rhoVnR;
            fluxR[IRU] = rhoVnR * uR + xsn * pR;
            fluxR[IRV] = rhoVnR * vR + ysn * pR;
            fluxR[IRW] = rhoVnR * wR + zsn * pR;
            fluxR[IRE] = rhoVnR * HR + vb  * pR;
            fer = peR * (vnR + vb);
        }
        else if (machR > one)
        {
            //! Supersonic flow from left to right, upwind by left.
            fluxR[IR ] = zero;
            fluxR[IRU] = zero;
            fluxR[IRV] = zero;
            fluxR[IRW] = zero;
            fluxR[IRE] = zero;
            fer = zero;
        }
        else
        {
            //! Subsonic flow.
            if (nPrecondition == 0)
            {
                //! Without precondition.
                RDouble fmass = - 0.25 * rR * cR * (machR - 1.0) * (machR - 1.0);
                RDouble tmp = (- vnR - 2.0 * cR) / gmR;
                fluxR[IR ] = fmass;
                fluxR[IRU] = fmass * (uR + xsn * tmp);
                fluxR[IRV] = fmass * (vR + ysn * tmp);
                fluxR[IRW] = fmass * (wR + zsn * tmp);
                fluxR[IRE] = fmass * (HR + vb  * tmp);
            }
            else
            {
                //! With precondition.
                RDouble fmass = 0.25 * rR * cR * (machR * machR + 1.0);
                RDouble tmp0  = - vnR / gmR;
                RDouble tmp1  = rR * cR / gmR;
                fluxR[IR ] = half * rhoVnR  - Mr2r * fmass;
                fluxR[IRU] = half * (rhoVnR * uR + xsn * pR) - Mr2r * (fmass * (uR + xsn * tmp0) + uR * tmp1);
                fluxR[IRV] = half * (rhoVnR * vR + ysn * pR) - Mr2r * (fmass * (vR + ysn * tmp0) + vR * tmp1);
                fluxR[IRW] = half * (rhoVnR * wR + zsn * pR) - Mr2r * (fmass * (wR + zsn * tmp0) + wR * tmp1);
                fluxR[IRE] = half * (rhoVnR * HR + vb  * pR) - Mr2r * (fmass * (HR + vb  * tmp0));
            }
        }

        //! Other equations besides N-S 5-equations, such as: multi-species, electronic.
        //! which are stored start from iEqua=5, 6, ...
        for (int iLaminar = nm; iLaminar < nLaminar; ++ iLaminar)
        {
            fluxR[iLaminar] = primR[iLaminar] * fluxR[IR];
        }

        //! Compute vibration-electron energy term.
        for (int it = 1; it < nTemperatureModel; ++ it)
        {
            fluxL[nLaminar + it] = primL[nLaminar + it] * fluxL[IR];
            fluxR[nLaminar + it] = primR[nLaminar + it] * fluxR[IR];
        }
        if (nChemical == 1)
        {
            fluxL[nLaminar + nTemperatureModel - 1] += fel;
            fluxR[nLaminar + nTemperatureModel - 1] += fer;
        }

        for (int iLaminar = 0; iLaminar < nEquation; ++ iLaminar)
        {
            flux[iLaminar][iLen] = fluxL[iLaminar] + fluxR[iLaminar];
        }
        if (nChemical == 1)
        {
            flux[nLaminar][iLen] = 0.0;
            if (nElectronIndex >= 0)
            {
                flux[nLaminar - 1][iLen] = 0.0;
            }
        }
    }

    delete [] fluxL;    fluxL = nullptr;
    delete [] fluxR;    fluxR = nullptr;
    delete [] primL;    primL = nullptr;
    delete [] primR;    primR = nullptr;
}


void Rotate_Vanleer_Scheme(FaceProxy *face_proxy, InviscidSchemeParameter *invSchemePara)
{
    RDouble v2l,vnl,hint_l,hl,rvnl,c2l,cl,machl;
    RDouble v2r,vnr,hint_r,hr,rvnr,c2r,cr,machr;
    RDouble fmass, tmp;

    using namespace IDX;
    using namespace GAS_SPACE;
    int neqn     = invSchemePara->GetNumberOfTotalEquation();
    int nl       = invSchemePara->GetNumberOfLaminar();
    int nlen     = invSchemePara->GetLength();
    int nm       = invSchemePara->Getnm();
    RDouble **ql = invSchemePara->GetLeftQ();
    RDouble **qr = invSchemePara->GetRightQ();
    RDouble *gamal = invSchemePara->GetLeftGama();
    RDouble *gamar = invSchemePara->GetRightGama();
    RDouble *xfn   = invSchemePara->GetFaceNormalX();
    RDouble *yfn   = invSchemePara->GetFaceNormalY();
    RDouble *zfn   = invSchemePara->GetFaceNormalZ();
    RDouble *vgn   = invSchemePara->GetFaceVelocity();
    RDouble **flux = invSchemePara->GetFlux();

    RDouble *fl    = new RDouble[nl];
    RDouble *fr    = new RDouble[nl];
    RDouble *priml = new RDouble[neqn];
    RDouble *primr = new RDouble[neqn];

    for ( int i = 0; i < nlen; ++ i )
    {
        for ( int m = 0; m < nl; ++ m )
        {
            flux[m][i] = 0 ;
        }

        for ( int m = 0; m < neqn; ++ m )
        {
            priml[m] = ql[m][i];
            primr[m] = qr[m][i];
        }

        RDouble &rl = priml[IR];
        RDouble &ul = priml[IU];
        RDouble &vl = priml[IV];
        RDouble &wl = priml[IW];
        RDouble &pl = priml[IP];

        RDouble &rr = primr[IR];
        RDouble &ur = primr[IU];
        RDouble &vr = primr[IV];
        RDouble &wr = primr[IW];
        RDouble &pr = primr[IP];

        v2l = ul * ul + vl * vl + wl * wl;
        v2r = ur * ur + vr * vr + wr * wr;

        RDouble gmL = gamal[i];
        RDouble gmR = gamar[i];

        hint_l = (gmL / (gmL - one)) * (pl / rl);
        hint_r = (gmR / (gmR - one)) * (pr / rr);

        //gas->ComputeEnthalpyByPrimitive(priml, gamal[i], hint_l);
        //gas->ComputeEnthalpyByPrimitive(primr, gamar[i], hint_r);

        hl = hint_l + half * v2l;
        hr = hint_r + half * v2r;

        //! Face normal
        RDouble n[3];
        n[0] = xfn[i];
        n[1] = yfn[i];
        n[2] = zfn[i];

        RDouble n1[3];

        if((qr[IU][i] - ql[IU][i]) == 0 && (qr[IV][i] - ql[IV][i]) == 0 && (qr[IW][i] - ql[IW][i]) == 0)
        {
            n1[0] = 1;
            n1[1] = 0;
            n1[2] = 0;
        }
        else
        {
            RDouble sqdv = sqrt(pow((qr[IU][i] - ql[IU][i] ),2) + pow((qr[IV][i] - ql[IV][i] ),2) + pow((qr[IW][i] - ql[IW][i]  ),2));
            n1[0] = (qr[IU][i] - ql[IU][i])/sqdv;
            n1[1] = (qr[IV][i] - ql[IV][i])/sqdv;
            n1[2] = (qr[IW][i] - ql[IW][i])/sqdv;
        }

        RDouble angle1 = GetAngleOfTwoVector(n1, n, 3);

        if ( cos(angle1) < 0)
        {
            n1[0] = -n1[0];
            n1[1] = -n1[1];
            n1[2] = -n1[2];
            angle1 = GetAngleOfTwoVector(n1, n, 3);
        }
        RDouble n2[3];
        if (abs(n1[0] - xfn[i]) < 1.0e-8 && abs(n1[1] - yfn[i]) < 1.0e-8 && abs(n1[2] - zfn[i]) < 1.0e-8)
        {

            n2[0] = 0;
            n2[1] = 1;
            n2[2] = 0;
        }
        else
        {
            RDouble tmpV[3];
            CrossProduct(n1, n, tmpV);
            CrossProduct(tmpV, n1, n2);
        }
        RDouble normn2 = GetNorm(n2, 3);
        if ( cos(GetAngleOfTwoVector(n2, n, 3)) < 0)
        {
            n2[0] = -n2[0]/normn2;
            n2[1] = -n2[1]/normn2;
            n2[2] = -n2[2]/normn2;
        }
        else
        {
            n2[0] = n2[0]/normn2;
            n2[1] = n2[1]/normn2;
            n2[2] = n2[2]/normn2;
        }

        RDouble angle2 = GetAngleOfTwoVector(n2, n, 3);
        for(int dir= 1; dir <=2; ++dir)
        {
            if (dir == 1)
            {
                xfn[i] = n1[0];
                yfn[i] = n1[1];
                zfn[i] = n1[2];
            }
            else
            {
                xfn[i] = n2[0];
                yfn[i] = n2[1];
                zfn[i] = n2[2];
            }
            vnl  = xfn[i] * ul + yfn[i] * vl + zfn[i] * wl - vgn[i];
            vnr  = xfn[i] * ur + yfn[i] * vr + zfn[i] * wr - vgn[i];
            rvnl = rl * vnl;
            rvnr = rr * vnr;

            c2l   = gamal[i] * pl / rl;
            cl    = sqrt( ABS(c2l) );
            machl = vnl / cl;

            c2r   = gamar[i] * pr / rr;
            cr    = sqrt( ABS(c2r) );
            machr = vnr / cr;

            if ( machl > one )
            {
                fl[IR ] = rvnl;
                fl[IRU] = rvnl * ul + xfn[i] * pl;
                fl[IRV] = rvnl * vl + yfn[i] * pl;
                fl[IRW] = rvnl * wl + zfn[i] * pl;
                fl[IRE] = rvnl * hl + vgn[i] * pl;
            }
            else if ( machl < - one )
            {
                fl[IR ] = zero;
                fl[IRU] = zero;
                fl[IRV] = zero;
                fl[IRW] = zero;
                fl[IRE] = zero;
            }
            else
            {
                fmass = 0.25 * rl * cl * ( machl + 1.0 ) * ( machl + 1.0 );
                tmp   = ( - vnl + 2.0 * cl ) / gamal[i];
                fl[IR ] = fmass;
                fl[IRU] = fmass * ( ul + xfn[i] * tmp );
                fl[IRV] = fmass * ( vl + yfn[i] * tmp );
                fl[IRW] = fmass * ( wl + zfn[i] * tmp );
                fl[IRE] = fmass * ( hl + vgn[i] * tmp );
            }

            for ( int m = nm; m < nl; ++ m )
            {
                fl[m] = priml[m] * fl[IR ];
            }

            //! Right contribution
            if ( machr < - one )
            {
                fr[IR ] = rvnr                ;
                fr[IRU] = rvnr * ur + xfn[i] * pr;
                fr[IRV] = rvnr * vr + yfn[i] * pr;
                fr[IRW] = rvnr * wr + zfn[i] * pr;
                fr[IRE] = rvnr * hr + vgn[i] * pr;
            }
            else if ( machr > one )
            {
                fr[IR ] = zero;
                fr[IRU] = zero;
                fr[IRV] = zero;
                fr[IRW] = zero;
                fr[IRE] = zero;
            }
            else
            {
                fmass = - 0.25 * rr * cr * ( machr - 1.0 ) * ( machr - 1.0 );
                tmp   = ( - vnr - 2.0 * cr ) / gamar[i];
                fr[IR ] = fmass                        ;
                fr[IRU] = fmass * ( ur + xfn[i] * tmp );
                fr[IRV] = fmass * ( vr + yfn[i] * tmp );
                fr[IRW] = fmass * ( wr + zfn[i] * tmp );
                fr[IRE] = fmass * ( hr + vgn[i] * tmp );
            }

            for ( int m = nm; m < nl; ++ m )
            {
                fr[m] = primr[m] * fr[IR ];
            }

            for ( int m = 0; m < nl; ++ m )
            {
                if (dir ==1 )
                {
                    flux[m][i]  = flux[m][i] + cos(angle1) * (fl[m] + fr[m]) ;
                }
                else
                {
                    flux[m][i]  = flux[m][i] + cos(angle2) * (fl[m] + fr[m]) ;
                }
            }
        }
    }

    delete [] fl   ;    fl = nullptr;
    delete [] fr   ;    fr = nullptr;
    delete [] priml;    priml = nullptr;
    delete [] primr;    primr = nullptr;
}

void Ausmdv_Scheme(FaceProxy *face_proxy, InviscidSchemeParameter *invSchemePara)
{
    int nTemperatureModel = invSchemePara->GetNumberOfTemperatureModel();
    int nEquation = invSchemePara->GetNumberOfTotalEquation();
    int nLaminar = invSchemePara->GetNumberOfLaminar();
    int nLen = invSchemePara->GetLength();
    int nm = invSchemePara->Getnm();
    RDouble **ql = invSchemePara->GetLeftQ();
    RDouble **qr = invSchemePara->GetRightQ();
    RDouble **tl = invSchemePara->GetLeftTemperature();
    RDouble **tr = invSchemePara->GetRightTemperature();

    RDouble *gamal = invSchemePara->GetLeftGama();
    RDouble *gamar = invSchemePara->GetRightGama();
    RDouble *xfn = invSchemePara->GetFaceNormalX();
    RDouble *yfn = invSchemePara->GetFaceNormalY();
    RDouble *zfn = invSchemePara->GetFaceNormalZ();
    RDouble *vgn = invSchemePara->GetFaceVelocity();
    RDouble **flux = invSchemePara->GetFlux();
    //! Compute the number of species.
    int nChemical = invSchemePara->GetNumberOfChemical();
    int nElectronIndex = invSchemePara->GetIndexOfElectron();

    using namespace IDX;

    RDouble *priml = new RDouble[nEquation];
    RDouble *primr = new RDouble[nEquation];
    RDouble temperaturesL[3], temperaturesR[3];

    for (int iLen = 0; iLen < nLen; ++ iLen)
    {
        for (int iEqn = 0; iEqn < nEquation; ++ iEqn)
        {
            priml[iEqn] = ql[iEqn][iLen];
            primr[iEqn] = qr[iEqn][iLen];
        }

        RDouble &rl = priml[IR];
        RDouble &ul = priml[IU];
        RDouble &vl = priml[IV];
        RDouble &wl = priml[IW];
        RDouble &pl = priml[IP];

        RDouble &rr = primr[IR];
        RDouble &ur = primr[IU];
        RDouble &vr = primr[IV];
        RDouble &wr = primr[IW];
        RDouble &pr = primr[IP];

        RDouble gama = half * (gamal[iLen] + gamar[iLen]);

        RDouble hintLeft, hintRight;
        RDouble v2l = ul * ul + vl * vl + wl * wl;
        RDouble v2r = ur * ur + vr * vr + wr * wr;

        for (int i = 0; i < 3; ++ i)
        {
            temperaturesL[i] = tl[i][iLen];
            temperaturesR[i] = tr[i][iLen];
        }
        gas->ComputeEnthalpyByPrimitive(priml, gamal[iLen], hintLeft, temperaturesL);
        gas->ComputeEnthalpyByPrimitive(primr, gamar[iLen], hintRight, temperaturesR);

        RDouble peL = 0.0, peR = 0.0;
        if (nChemical > 0 && nTemperatureModel > 1)
        {
            peL = gas->GetElectronPressure(priml, temperaturesL[ITE]);
            peR = gas->GetElectronPressure(primr, temperaturesR[ITE]);
        }

        RDouble hl = hintLeft  + half * v2l;
        RDouble hr = hintRight + half * v2r;

        RDouble vnl = xfn[iLen] * ul + yfn[iLen] * vl + zfn[iLen] * wl - vgn[iLen];
        RDouble vnr = xfn[iLen] * ur + yfn[iLen] * vr + zfn[iLen] * wr - vgn[iLen];

        RDouble orl = 1.0 / (rl + SMALL);
        RDouble orr = 1.0 / (rr + SMALL);

        RDouble c2l =  gama * pl * orl;
        RDouble c2r =  gama * pr * orr;

        //! Here,the vn_l and vn_r are normal velocities in the inertial coordinate system.
        //! vnl and vnr are normal velocities relative to the moving interface.
        RDouble vn_l = vnl + vgn[iLen];
        RDouble vn_r = vnr + vgn[iLen];

        //! Fourth speed interface.
        RDouble cl = sqrt(c2l);
        RDouble cr = sqrt(c2r);
        RDouble cm = MAX(cl,cr);
        RDouble ocm = one / (cm + SMALL);

        RDouble ssw = zero;
        //! Determine shock switch.
        if ((vnl > cl && vnr < cr) || (vnl > -cl && vnr < -cr))
        {
            ssw = one;
        }
        RDouble ssw_a = one - ssw;

        //! Haenel / van Leer.
        RDouble ml  =  vnl / cl;
        RDouble mr  =  vnr / cr;

        RDouble mla =  ABS(ml);
        RDouble mra =  ABS(mr);

        RDouble mpl, ppl, mmr, pmr;

        RDouble mp1 = half * (ml + mla);
        RDouble mm1 = half * (mr - mra);

        if (mla >= one)
        {
            mpl = mp1;
            ppl = half * (one + ml / (mla + SMALL));
        }
        else
        {
            RDouble tmp = fourth * SQR(ml + one);
            mpl = tmp;
            ppl = tmp * (two - ml);
        }

        if (mra >= one)
        {
            mmr = mm1;
            pmr = half * (one - mr / (mra + SMALL));
        }
        else
        {
            RDouble tmp = fourth * SQR(mr - one);

            mmr = - tmp;
            pmr = tmp * (two + mr);
        }

        RDouble p12_v =  ppl * pl + pmr * pr;
        RDouble pe12 = ppl * peL + pmr * peR;
        RDouble rvnl_v = cl * rl * mpl;
        RDouble rvnr_v = cr * rr * mmr;

        //! AUSMDV.
        ml  =  vnl * ocm;
        mla =  ABS(ml);

        mr  =  vnr * ocm;
        mra =  ABS(mr);

        mp1 = half * (ml + mla);
        mm1 = half * (mr - mra);

        RDouble r1  = pl * orl;
        RDouble r2  = pr * orr;
        RDouble r3  = two / (r1 + r2);
        RDouble alphal = r1 * r3;
        RDouble alphar = r2 * r3;

        if (mla >= one)
        {
            mpl = mp1;
            ppl = half * (one + ml / (mla + SMALL));
        }
        else
        {
            RDouble tmp = fourth * SQR(ml + one);

            mpl = alphal * tmp + (one - alphal) * mp1;
            ppl = tmp * (two - ml);
        }

        if (mra >= one)
        {
            mmr = mm1;
            pmr = half * (one - mr / (mra + SMALL));
        }
        else
        {
            RDouble tmp = fourth * SQR(mr - one);
            mmr = - alphar * tmp + (one - alphar) * mm1;
            pmr = tmp * (two + mr);
        }
        RDouble rvnl_d = cm * rl * mpl;
        RDouble rvnr_d = cm * rr * mmr;

        RDouble rvn    = rvnl_d + rvnr_d;
        rvnl_d = half * (rvn + ABS(rvn));
        rvnr_d = half * (rvn - ABS(rvn));

        RDouble rvnl = rvnl_v;
        RDouble rvnr = rvnr_v;
        RDouble p12  = p12_v ;

        //! AUSMD/AUSMV switch.
        RDouble kSwitch = 10.0, cEfix = 0.125;
        RDouble s = half * MIN(static_cast <RDouble> (one), kSwitch * ABS(pr-pl) / MIN(pl,pr));

        RDouble rvvn = ssw * (rvnl_v * vn_l + rvnr_v * vn_r);
        rvvn += (half - s) * ssw_a *      (rvnl_d * vn_l +   rvnr_d * vn_r)
              + (half + s) * ssw_a * cm * (rl * mpl * vn_l + rr * mmr * vn_r);

        //! Entropy fix.
        bool atmp = ((vnl - cl) < 0.0) && ((vnr - cr) > 0.0);
        bool btmp = ((vnl + cl) < 0.0) && ((vnr + cr) > 0.0);

        if ((atmp && !btmp) || (!atmp && btmp))
        {
            if (atmp && !btmp)
            {
                RDouble tmp   = ssw_a * cEfix * ((vnr - cr) - (vnl - cl));
                rvnl += tmp * rl;
                rvnr -= tmp * rr;
            }
            else
            {
                RDouble tmp   = ssw_a * cEfix * ((vnr + cr) - (vnl + cl));
                rvnl += tmp * rl;
                rvnr -= tmp * rr;
            }
        }

        RDouble dul = ul - xfn[iLen] * vn_l;
        RDouble dvl = vl - yfn[iLen] * vn_l;
        RDouble dwl = wl - zfn[iLen] * vn_l;

        RDouble dur = ur - xfn[iLen] * vn_r;
        RDouble dvr = vr - yfn[iLen] * vn_r;
        RDouble dwr = wr - zfn[iLen] * vn_r;

        flux[IR ][iLen] = rvnl       + rvnr      ;
        flux[IRU][iLen] = rvnl * dul + rvnr * dur;
        flux[IRV][iLen] = rvnl * dvl + rvnr * dvr;
        flux[IRW][iLen] = rvnl * dwl + rvnr * dwr;
        flux[IRE][iLen] = rvnl * hl  + rvnr * hr ;

        for (int iLaminar = nm; iLaminar < nLaminar; ++ iLaminar)
        {
            flux[iLaminar][iLen] = rvnl * priml[iLaminar] + rvnr * primr[iLaminar];
        }
        //! Compute the vibration-electron energy term.
        for (int it = 1; it < nTemperatureModel; ++ it)
        {
            flux[nLaminar + it][iLen] = rvnl * priml[nLaminar + it] + rvnr * primr[nLaminar + it];
        }
        if (nChemical == 1)
        {
            flux[nLaminar + nTemperatureModel - 1][iLen] += peL * vnl + peR * vnr + pe12 * vgn[iLen];
        }

        rvvn += p12;

        flux[IRU][iLen] += rvvn * xfn[iLen];
        flux[IRV][iLen] += rvvn * yfn[iLen];
        flux[IRW][iLen] += rvvn * zfn[iLen];
        flux[IRE][iLen] += p12  * vgn[iLen];

        if (nChemical > 0)
        {
            flux[nLaminar][iLen] = 0.0;
            if (nElectronIndex >= 0)
            {
                flux[nLaminar - 1][iLen] = 0.0;
            }
        }
    }

    delete [] priml;    priml = nullptr;
    delete [] primr;    primr = nullptr;

}

void Ausmpw_Scheme(FaceProxy *face_proxy, InviscidSchemeParameter *invSchemePara)
{
    int nTemperatureModel = invSchemePara->GetNumberOfTemperatureModel();
    int nEquation = invSchemePara->GetNumberOfTotalEquation();
    int nLaminar = invSchemePara->GetNumberOfLaminar();
    int nLen = invSchemePara->GetLength();
    int nm = invSchemePara->Getnm();
    //! Compute the number of species.
    int nChemical = invSchemePara->GetNumberOfChemical();
    int nElectronIndex = invSchemePara->GetIndexOfElectron();

    RDouble **ql = invSchemePara->GetLeftQ();
    RDouble **qr = invSchemePara->GetRightQ();
    RDouble *gamal = invSchemePara->GetLeftGama();
    RDouble *gamar = invSchemePara->GetRightGama();
    RDouble *xfn   = invSchemePara->GetFaceNormalX();
    RDouble *yfn   = invSchemePara->GetFaceNormalY();
    RDouble *zfn   = invSchemePara->GetFaceNormalZ();
    RDouble *vgn   = invSchemePara->GetFaceVelocity();
    RDouble **flux = invSchemePara->GetFlux();
    RDouble **tl = invSchemePara->GetLeftTemperature();
    RDouble **tr = invSchemePara->GetRightTemperature();

    using namespace IDX;
    using namespace GAS_SPACE;

    RDouble *priml = new RDouble[nEquation];
    RDouble *primr = new RDouble[nEquation];
    RDouble temperaturesL[3], temperaturesR[3];

    for (int iLen = 0; iLen < nLen; ++ iLen)
    {
        for (int iEqn = 0; iEqn < nEquation; ++ iEqn)
        {
            priml[iEqn] = ql[iEqn][iLen];
            primr[iEqn] = qr[iEqn][iLen];
        }

        RDouble &rl = priml[IR];
        RDouble &ul = priml[IU];
        RDouble &vl = priml[IV];
        RDouble &wl = priml[IW];
        RDouble &pl = priml[IP];

        RDouble &rr = primr[IR];
        RDouble &ur = primr[IU];
        RDouble &vr = primr[IV];
        RDouble &wr = primr[IW];
        RDouble &pr = primr[IP];

        RDouble v2l = ul * ul + vl * vl + wl * wl;
        RDouble v2r = ur * ur + vr * vr + wr * wr;

        RDouble hint_l, hint_r;
        for (int i = 0; i < 3; ++ i)
        {
            temperaturesL[i] = tl[i][iLen];
            temperaturesR[i] = tr[i][iLen];
        }
        gas->ComputeEnthalpyByPrimitive(priml, gamal[iLen], hint_l, temperaturesL);
        gas->ComputeEnthalpyByPrimitive(primr, gamar[iLen], hint_r, temperaturesR);

        RDouble peL = 0.0, peR = 0.0;
        if (nChemical > 0 && nTemperatureModel > 1)
        {
            peL = gas->GetElectronPressure(priml, temperaturesL[ITE]);
            peR = gas->GetElectronPressure(primr, temperaturesR[ITE]);
        }

        RDouble hl  = hint_l + half * v2l;
        RDouble hr  = hint_r + half * v2r;

        RDouble vnl  = xfn[iLen] * ul + yfn[iLen] * vl + zfn[iLen] * wl - vgn[iLen];
        RDouble vnr  = xfn[iLen] * ur + yfn[iLen] * vr + zfn[iLen] * wr - vgn[iLen];

        RDouble orl = 1.0 / (rl + SMALL);
        RDouble orr = 1.0 / (rr + SMALL);

        RDouble c2l = gamal[iLen] * pl * orl;
        RDouble c2r = gamar[iLen] * pr * orr;
        RDouble cl  = sqrt(c2l);
        RDouble cr  = sqrt(c2r);
        RDouble cm  = half * (cl + cr);

        RDouble ocm = one / (cm + SMALL);

        RDouble ml  =  vnl * ocm;
        RDouble mr  =  vnr * ocm;

        RDouble fm4ml = FMSplit4(ml, zero,  one);
        RDouble fm4mr = FMSplit4(mr, zero, -one);

        RDouble fp5ml = FPSplit5(ml, zero,  one);
        RDouble fp5mr = FPSplit5(mr, zero, -one);

        RDouble m12  = fm4ml + fm4mr;
        RDouble p12  = fp5ml * pl + fp5mr * pr;
        RDouble pe12 = fp5ml * peL + fp5mr * peR;

        RDouble minplr = MIN(pl / (pr + SMALL), pr / (pl + SMALL));
        RDouble fw     = one - minplr * minplr * minplr;

        RDouble op12 = one / (p12 + SMALL);

        RDouble fl = zero;
        RDouble fr = zero;

        if (ABS(ml) < one)
        {
            fl = pl * op12 - one;
        }

        if (ABS(mr) < one)
        {
            fr = pr * op12 - one;
        }

        RDouble mpl, mmr;
        if (m12 >= zero)
        {
            RDouble one_fr = one + fr;
            mpl = fm4ml + fm4mr * ((one - fw) * one_fr - fl);
            mmr = fm4mr * fw * one_fr;
        }
        else
        {
            RDouble one_fl = one + fl;
            mpl = fm4ml * fw * one_fl;
            mmr = fm4mr + fm4ml * ((one - fw) * one_fl - fr);
        }

        flux[IR ][iLen] = cm * (mpl*rl    + mmr*rr );
        flux[IRU][iLen] = cm * (mpl*rl*ul + mmr*rr*ur) + p12 * xfn[iLen];
        flux[IRV][iLen] = cm * (mpl*rl*vl + mmr*rr*vr) + p12 * yfn[iLen];
        flux[IRW][iLen] = cm * (mpl*rl*wl + mmr*rr*wr) + p12 * zfn[iLen];
        flux[IRE][iLen] = cm * (mpl*rl*hl + mmr*rr*hr) + p12 * vgn[iLen];

        for (int iLaminar = nm; iLaminar < nLaminar; ++ iLaminar)
        {
            flux[iLaminar][iLen] = cm * (mpl*rl*priml[iLaminar] + mmr*rr*primr[iLaminar]);
        }

        //! Compute vibration-electron energy term.
        for (int it = 1; it < nTemperatureModel; ++ it)
        {
            flux[nLaminar + it][iLen] = cm * (mpl*rl*priml[nLaminar + it] + mmr*rr*primr[nLaminar + it]);
        }
        if (nChemical == 1)
        {
            flux[nLaminar + nTemperatureModel - 1][iLen] += cm * (mpl*peL + mmr*peR) + pe12 * vgn[iLen];
        }

        if (nChemical > 0)
        {
            flux[nLaminar][iLen] = 0.0;
            if (nElectronIndex >= 0)
            {
                flux[nLaminar - 1][iLen] = 0.0;
            }
        }
    }
    delete [] priml;    priml = nullptr;
    delete [] primr;    primr = nullptr;
}

void AusmpwPlus_Scheme(FaceProxy *face_proxy, InviscidSchemeParameter *invSchemePara)
{
    int nTemperatureModel = invSchemePara->GetNumberOfTemperatureModel();
    int nEquation = invSchemePara->GetNumberOfTotalEquation();
    int nLaminar = invSchemePara->GetNumberOfLaminar();
    int nLen = invSchemePara->GetLength();
    int nm = invSchemePara->Getnm();
    //! Compute the number of species.
    int nChemical = invSchemePara->GetNumberOfChemical();
    int nElectronIndex = invSchemePara->GetIndexOfElectron();
    RDouble AusmpwPlusLimiter = invSchemePara->GetAusmpwPlusLimiter();

    RDouble **ql = invSchemePara->GetLeftQ();
    RDouble **qr = invSchemePara->GetRightQ();
    RDouble *gamal = invSchemePara->GetLeftGama();
    RDouble *gamar = invSchemePara->GetRightGama();
    RDouble *xfn   = invSchemePara->GetFaceNormalX();
    RDouble *yfn   = invSchemePara->GetFaceNormalY();
    RDouble *zfn   = invSchemePara->GetFaceNormalZ();
    RDouble *vgn   = invSchemePara->GetFaceVelocity();
    RDouble **flux = invSchemePara->GetFlux();
    RDouble **tl = invSchemePara->GetLeftTemperature();
    RDouble **tr = invSchemePara->GetRightTemperature();

    using namespace IDX;
    using namespace GAS_SPACE;

    RDouble *priml = new RDouble[nEquation];
    RDouble *primr = new RDouble[nEquation];
    RDouble temperaturesL[3], temperaturesR[3];

    for (int iLen = 0; iLen < nLen; ++ iLen)
    {
        for (int iEqn = 0; iEqn < nEquation; ++ iEqn)
        {
            priml[iEqn] = ql[iEqn][iLen];
            primr[iEqn] = qr[iEqn][iLen];
        }

        RDouble &rl = priml[IR];
        RDouble &ul = priml[IU];
        RDouble &vl = priml[IV];
        RDouble &wl = priml[IW];
        RDouble &pl = priml[IP];

        RDouble &rr = primr[IR];
        RDouble &ur = primr[IU];
        RDouble &vr = primr[IV];
        RDouble &wr = primr[IW];
        RDouble &pr = primr[IP];

        RDouble v2l = ul * ul + vl * vl + wl * wl;
        RDouble v2r = ur * ur + vr * vr + wr * wr;

        RDouble hint_l, hint_r;
        for (int i = 0; i < 3; ++ i)
        {
            temperaturesL[i] = tl[i][iLen];
            temperaturesR[i] = tr[i][iLen];
        }
        gas->ComputeEnthalpyByPrimitive(priml, gamal[iLen], hint_l, temperaturesL);
        gas->ComputeEnthalpyByPrimitive(primr, gamar[iLen], hint_r, temperaturesR);

        RDouble peL = 0.0, peR = 0.0;
        if (nChemical > 0 && nTemperatureModel > 1)
        {
            peL = gas->GetElectronPressure(priml, temperaturesL[ITE]);
            peR = gas->GetElectronPressure(primr, temperaturesR[ITE]);
        }

        RDouble hl  = hint_l + half * v2l;
        RDouble hr  = hint_r + half * v2r;

        //RDouble gama2l = (gamal[iLen] - 1.0)/(gamal[iLen] + 1.0);
        //RDouble gama2r = (gamar[iLen] - 1.0)/(gamar[iLen] + 1.0);
        //RDouble cs  = half * ABS((gama2l + gama2r) * (hl + hr));             //!!!!
        //RDouble c12;
        //if (vnl + vnr > 0.0)
        //{
        //    c12 = cs / max(abs(vnl), sqrt(cs));
        //}
        //else
        //{
        //    c12 = cs / max(abs(vnr), sqrt(cs));
        //}

        RDouble al = sqrt(gamal[iLen] * (pl / (rl + SMALL)));
        RDouble ar = sqrt(gamar[iLen] * (pr / (rr + SMALL)));

        RDouble vnl= xfn[iLen] * ul + yfn[iLen] * vl + zfn[iLen] * wl - vgn[iLen];
        RDouble vnr= xfn[iLen] * ur + yfn[iLen] * vr + zfn[iLen] * wr - vgn[iLen];
        RDouble ml= abs(vnl / al);
        RDouble mr= abs(vnr / ar);

        RDouble minplr = MIN(pl / (pr + SMALL), pr / (pl + SMALL));

        RDouble c12, ocm;

        /*if( minplr < 0.5 && (( ml > 1.0 && mr < 1.0 ) || ( mr > 1.0 && ml < 1.0 )) )
        {
            c12 = max(sqrt(abs(vnl * vnr)), min(al, ar));

            if (vnl + vnr > 0.0)
            {
                c12 *= c12 / max(abs(vnl), c12);
            }
            else
            {
                c12 *= c12 / max(abs(vnr), c12);
            }
        }
        else
        {*/
            c12 = half * (al + ar);
        //}

        //if( minplr < 0.5)
        //{
        //    RDouble du = ur - ul;
        //    RDouble dv = vr - vl;
        //    RDouble dw = wr - wl;
        //    
        //    c12 = max(sqrt(du * du + dv * dv + dw * dw), SMALL);
        //    if(xfn[iLen] * du + yfn[iLen] * du + zfn[iLen] * du < 0.0)
        //    {
        //        c12 = -c12;
        //        du /= c12;
        //        dv /= c12;
        //        dw /= c12;
        //    }
        //    else
        //    {
        //        du /= c12;
        //        dv /= c12;
        //        dw /= c12;
        //    }

        //    vnl  = du * ul + dv * vl + dw * wl - vgn[iLen];
        //    vnr  = du * ur + dv * vr + dw * wr - vgn[iLen];

        //    ml = abs(vnl / al);
        //    mr = abs(vnr / ar);

        //    if(( ml > 1.0 && mr < 1.0 ) || ( mr > 1.0 && ml < 1.0 ))
        //    {
        //        c12 = max(sqrt(abs(vnl * vnr)), min(al, ar));
        //        if (vnl + vnr > 0.0)
        //        {
        //            c12 *= c12 / max(abs(vnl), c12);
        //        }
        //        else
        //        {
        //            c12 *= c12 / max(abs(vnr), c12);
        //        }
        //    }
        //    else
        //    {
        //        c12 = half * (al + ar);
        //    }
        //}
        //else
        //{
        //    c12 = half * (al + ar);
        //}

        //vnl  = xfn[iLen] * ul + yfn[iLen] * vl + zfn[iLen] * wl - vgn[iLen];
        //vnr  = xfn[iLen] * ur + yfn[iLen] * vr + zfn[iLen] * wr - vgn[iLen];

       // RDouble du = ur - ul;
       // RDouble dv = vr - vl;
       // RDouble dw = wr - wl;
       // c12 = max(sqrt(du * du + dv * dv + dw * dw), SMALL);
       // if(xfn[iLen] * du + yfn[iLen] * du + zfn[iLen] * du < 0.0)
       // {
       //     c12 = -c12;
       // }
       // du /= c12;
       // dv /= c12;
       // dw /= c12;
       // 
       // vnl  = du * ul + dv * vl + dw * wl - vgn[iLen];
       // vnr  = du * ur + dv * vr + dw * wr - vgn[iLen];

       // ml = abs(vnl / al);
       // mr = abs(vnr / ar);

       //if(( ml > 1.0 && mr < 1.0 ) || ( mr > 1.0 && ml < 1.0 ))
       // {
       //     c12 = max(sqrt(abs(vnl * vnr)), min(al, ar));
       //     if (vnl + vnr > 0.0)
       //     {
       //         c12 *= c12 / max(abs(vnl), c12);
       //     }
       //     else
       //     {
       //         c12 *= c12 / max(abs(vnr), c12);
       //     }
       // }
       // else
       // {
       //     c12 = half * (al + ar);
       // }

       // vnl  = xfn[iLen] * ul + yfn[iLen] * vl + zfn[iLen] * wl - vgn[iLen];
       // vnr  = xfn[iLen] * ur + yfn[iLen] * vr + zfn[iLen] * wr - vgn[iLen];

        ocm = 1.0 / c12;

        ml  =  vnl * ocm;
        mr  =  vnr * ocm;

        RDouble fm4ml = FMmpn(ml, one);
        RDouble fm4mr = FMmpn(mr, -one);

        RDouble fp5ml = FMppn(ml, one);
        RDouble fp5mr = FMppn(mr, -one);

        RDouble m12  = fm4ml + fm4mr;
        RDouble p12  = fp5ml * pl + fp5mr * pr;
        RDouble pe12 = fp5ml * peL + fp5mr * peR;

        RDouble fw     = one - minplr * minplr * minplr * AusmpwPlusLimiter;

        RDouble op12 = one / (p12 + SMALL);

        RDouble fl = zero;
        RDouble fr = zero;

        if (p12 != 0.0)
        {
            if (ABS(ml) <= one)
            {
                fl = pl * op12 - one;
            }

            if (ABS(mr) <= one)
            {
                fr = pr * op12 - one;
            }
        }

        RDouble mpl, mmr;
        if (m12 >= zero)
        {
            RDouble one_fr = one + fr;
            mpl = fm4ml + fm4mr * ((one - fw) * one_fr - fl);
            mmr = fm4mr * fw * one_fr;
        }
        else
        {
            RDouble one_fl = one + fl;
            mpl = fm4ml * fw * one_fl;
            mmr = fm4mr + fm4ml * ((one - fw) * one_fl - fr);
        }

        flux[IR ][iLen] = c12 * (mpl * rl      + mmr * rr);
        flux[IRU][iLen] = c12 * (mpl * rl * ul + mmr * rr * ur) + p12 * xfn[iLen];
        flux[IRV][iLen] = c12 * (mpl * rl * vl + mmr * rr * vr) + p12 * yfn[iLen];
        flux[IRW][iLen] = c12 * (mpl * rl * wl + mmr * rr * wr) + p12 * zfn[iLen];
        flux[IRE][iLen] = c12 * (mpl * rl * hl + mmr * rr * hr) + p12 * vgn[iLen];

        for (int iLaminar = nm; iLaminar < nLaminar; ++ iLaminar)
        {
            flux[iLaminar][iLen] = c12 * (mpl * rl * priml[iLaminar] + mmr * rr * primr[iLaminar]);
        }

        //! Compute vibration-electron energy term.
        for (int it = 1; it < nTemperatureModel; ++ it)
        {
            flux[nLaminar + it][iLen] = c12 * (mpl * rl * priml[nLaminar + it] + mmr * rr * primr[nLaminar + it]);
        }
        if (nChemical == 1)
        {
            flux[nLaminar + nTemperatureModel - 1][iLen] += c12 * (mpl * peL + mmr * peR) + pe12 * vgn[iLen];
        }

        if (nChemical > 0)
        {
            flux[nLaminar][iLen] = 0.0;
            if (nElectronIndex >= 0)
            {
                flux[nLaminar - 1][iLen] = 0.0;
            }
        }
    }
    delete [] priml;    priml = nullptr;
    delete [] primr;    primr = nullptr;
}

void Hlle_Scheme(FaceProxy *face_proxy, InviscidSchemeParameter *invSchemePara)
{
    int nTemperatureModel = invSchemePara->GetNumberOfTemperatureModel();
    int nEquation = invSchemePara->GetNumberOfTotalEquation();
    int nLaminar = invSchemePara->GetNumberOfLaminar();
    int nLen = invSchemePara->GetLength();
    int nm = invSchemePara->Getnm();
    //! Compute the number of species.
    int nChemical = invSchemePara->GetNumberOfChemical();
    int nElectronIndex = invSchemePara->GetIndexOfElectron();

    RDouble **ql = invSchemePara->GetLeftQ();
    RDouble **qr = invSchemePara->GetRightQ();
    RDouble *gamal = invSchemePara->GetLeftGama();
    RDouble *gamar = invSchemePara->GetRightGama();
    RDouble *xfn   = invSchemePara->GetFaceNormalX();
    RDouble *yfn   = invSchemePara->GetFaceNormalY();
    RDouble *zfn   = invSchemePara->GetFaceNormalZ();
    RDouble *vgn   = invSchemePara->GetFaceVelocity();
    RDouble **flux = invSchemePara->GetFlux();
    RDouble **tl = invSchemePara->GetLeftTemperature();
    RDouble **tr = invSchemePara->GetRightTemperature();

    using namespace IDX;
    using namespace GAS_SPACE;

    RDouble *qcsvl = new RDouble[nEquation];
    RDouble *qcsvr = new RDouble[nEquation];
    RDouble *fl    = new RDouble[nEquation];
    RDouble *fr    = new RDouble[nEquation];

    RDouble *priml = new RDouble[nEquation];
    RDouble *primr = new RDouble[nEquation];
    RDouble temperaturesL[3], temperaturesR[3];

    for (int iLen = 0; iLen < nLen; ++ iLen)
    {
        for (int iEqn = 0; iEqn < nEquation; ++ iEqn)
        {
            priml[iEqn] = ql[iEqn][iLen];
            primr[iEqn] = qr[iEqn][iLen];
        }

        RDouble &rl = priml[IR];
        RDouble &ul = priml[IU];
        RDouble &vl = priml[IV];
        RDouble &wl = priml[IW];
        RDouble &pl = priml[IP];

        RDouble &rr = primr[IR];
        RDouble &ur = primr[IU];
        RDouble &vr = primr[IV];
        RDouble &wr = primr[IW];
        RDouble &pr = primr[IP];

        RDouble v2l = ul * ul + vl * vl + wl * wl;
        RDouble v2r = ur * ur + vr * vr + wr * wr;

        RDouble hint_l, hint_r;
        for (int i = 0; i < 3; ++ i)
        {
            temperaturesL[i] = tl[i][iLen];
            temperaturesR[i] = tr[i][iLen];
        }
        gas->ComputeEnthalpyByPrimitive(priml, gamal[iLen], hint_l, temperaturesL);
        gas->ComputeEnthalpyByPrimitive(primr, gamar[iLen], hint_r, temperaturesR);

        RDouble hl  = hint_l + half * v2l;
        RDouble hr  = hint_r + half * v2r;

        RDouble ratio = sqrt(rr / rl);
        RDouble coef  = 1.0 / (1.0 + ratio);

        RDouble rm = sqrt(rl * rr);
        RDouble um = (ul + ur * ratio) * coef;
        RDouble vm = (vl + vr * ratio) * coef;
        RDouble wm = (wl + wr * ratio) * coef;
        RDouble pm = (pl + pr * ratio) * coef;
        
        RDouble gama  = (gamal[iLen] + gamar[iLen] * ratio) * coef;

        RDouble vnl  = xfn[iLen] * ul + yfn[iLen] * vl + zfn[iLen] * wl - vgn[iLen];
        RDouble vnr  = xfn[iLen] * ur + yfn[iLen] * vr + zfn[iLen] * wr - vgn[iLen];
        RDouble rvnl = rl * vnl;
        RDouble rvnr = rr * vnr;

        RDouble vn   = xfn[iLen] * um + yfn[iLen] * vm + zfn[iLen] * wm - vgn[iLen];
        RDouble cl   = sqrt(gamal[iLen] * pl / (rl + SMALL));
        RDouble cr   = sqrt(gamar[iLen] * pr / (rr + SMALL));
        RDouble cm   = sqrt(gama     * pm / (rm + SMALL));

        fl[IR ] = rvnl                      ;
        fl[IRU] = rvnl * ul + xfn[iLen] * pl;
        fl[IRV] = rvnl * vl + yfn[iLen] * pl;
        fl[IRW] = rvnl * wl + zfn[iLen] * pl;
        fl[IRE] = rvnl * hl + vgn[iLen] * pl;

        fr[IR ] = rvnr                      ;
        fr[IRU] = rvnr * ur + xfn[iLen] * pr;
        fr[IRV] = rvnr * vr + yfn[iLen] * pr;
        fr[IRW] = rvnr * wr + zfn[iLen] * pr;
        fr[IRE] = rvnr * hr + vgn[iLen] * pr;

        for (int iLaminar = nm; iLaminar < nLaminar; ++ iLaminar)
        {
            fl[iLaminar] = priml[iLaminar] * fl[IR];
            fr[iLaminar] = primr[iLaminar] * fr[IR];
        }

        if (nTemperatureModel > 1)
        {
            gas->GetTemperature(priml, temperaturesL[ITT], temperaturesL[ITV], temperaturesL[ITE]);
            gas->GetTemperature(primr, temperaturesR[ITT], temperaturesR[ITV], temperaturesR[ITE]);
        }
        gas->Primitive2Conservative(priml, gamal[iLen], temperaturesL[ITV], temperaturesL[ITE], qcsvl);
        gas->Primitive2Conservative(primr, gamar[iLen], temperaturesR[ITV], temperaturesR[ITE], qcsvr);

        //! Compute vibration-electron energy term.
        for (int it = 1; it < nTemperatureModel; ++ it)
        {
            fl[nLaminar + it] = priml[nLaminar + it] * fl[IR];
            fr[nLaminar + it] = primr[nLaminar + it] * fr[IR];
        }
        RDouble peL = 0.0, peR = 0.0;
        if (nChemical > 0 && nTemperatureModel > 1)
        {
            peL = gas->GetElectronPressure(priml, temperaturesL[ITE]);
            peR = gas->GetElectronPressure(primr, temperaturesR[ITE]);
            fl[nLaminar + nTemperatureModel - 1] += peL * (vnl + vgn[iLen]);
            fr[nLaminar + nTemperatureModel - 1] += peR * (vnr + vgn[iLen]);
        }

        RDouble bm = MIN(static_cast<RDouble>(zero), MIN(vn - cm, vnl - cl));
        RDouble bp = MAX(static_cast<RDouble>(zero), MAX(vn + cm, vnr + cr));
    
        RDouble c1 = bp * bm / (bp - bm + SMALL);
        RDouble c2 = - half * (bp + bm) / (bp - bm + SMALL);

        for (int iLaminar = 0; iLaminar < nEquation; ++ iLaminar)
        {
            flux[iLaminar][iLen] = half * (fl[iLaminar]    + fr[iLaminar]) 
                              + c1   * (qcsvr[iLaminar] - qcsvl[iLaminar]) 
                              + c2   * (fr[iLaminar]    - fl[iLaminar]);
        }
        if (nChemical == 1)
        {
            flux[nLaminar][iLen] = 0.0;
            if (nElectronIndex >= 0)
            {
                flux[nLaminar - 1][iLen] = 0.0;
            }
        }
    }

    delete [] qcsvl;    qcsvl = nullptr;
    delete [] qcsvr;    qcsvr = nullptr;
    delete [] fl;    fl = nullptr;
    delete [] fr;    fr = nullptr;
    delete [] priml;    priml = nullptr;
    delete [] primr;    primr = nullptr;
}

RDouble FMSplit1(const RDouble &mach, const RDouble &ipn)
{
    return half * (mach + ipn * ABS(mach));
}

RDouble FMSplit2(const RDouble &mach, const RDouble &ipn)
{
    return ipn * fourth * SQR(mach + ipn);
}

RDouble FMmpn(const RDouble &mach, const RDouble &ipn)
{
    if (ABS(mach) >= one)
    {
        return FMSplit1(mach,ipn);
    }
    return FMSplit2(mach,ipn);
}

RDouble FMppn(const RDouble &mach, const RDouble &ipn)
{
    RDouble macha = ABS(mach);
    if (macha >= one)
    {
        return  FMSplit1(mach,ipn) / mach ;
    }

    return FMSplit2(mach,ipn) * ((two * ipn - mach));
}

RDouble FMSplit4(const RDouble &mach, const RDouble &beta, const RDouble &ipn)
{
    if (ABS(mach) >= one)
    {
        return FMSplit1(mach,ipn);
    }
    return FMSplit2(mach,ipn) * (one - ipn * sixteen * beta * FMSplit2(mach,-ipn));
}

RDouble FPSplit5(const RDouble &mach, const RDouble &alpha, const RDouble &ipn)
{
    RDouble macha = ABS(mach);
    if (macha >= one)
    {
        return  FMSplit1(mach,ipn) / mach ;
    }

    return FMSplit2(mach,ipn) * ((two * ipn - mach) - ipn * sixteen * alpha * mach * FMSplit2(mach,-ipn));
}

void RoeAveragedVariables(RDouble rhoLeft, RDouble uLeft, RDouble vLeft, RDouble wLeft, RDouble pLeft, RDouble gamaLeft, 
                          RDouble rhoRight, RDouble uRight, RDouble vRight, RDouble wRight, RDouble pRight, RDouble gamaRight,
                          RDouble &rhoRoe, RDouble &uRoe, RDouble &vRoe, RDouble &wRoe, RDouble &pRoe)
{
    RDouble hLeft  = gamaLeft  / (gamaLeft  - 1.0) * pLeft / rhoLeft + half * (uLeft * uLeft + vLeft * vLeft + wLeft * wLeft);
    RDouble hRight = gamaRight / (gamaRight - 1.0) * pRight / rhoRight + half * (uRight * uRight + vRight * vRight + wRight * wRight);

    RDouble tmp0 = sqrt(rhoRight / rhoLeft);
    RDouble tmp1 = 1.0 / (1.0 + tmp0);

    RDouble rm = sqrt(rhoRight * rhoLeft);
    RDouble um = (uLeft + uRight * tmp0) * tmp1;
    RDouble vm = (vLeft + vRight * tmp0) * tmp1;
    RDouble wm = (wLeft + wRight * tmp0) * tmp1;
    RDouble hm = (hLeft + hRight * tmp0) * tmp1;
    RDouble pm = rm * (hm - half * SQR(um, vm, wm)) * (gamaLeft  - 1.0) / gamaLeft;

    rhoRoe = rm;
    uRoe = um;
    vRoe = vm;
    wRoe = wm;
    pRoe = pm;
}

InviscidSchemeParameter::InviscidSchemeParameter()
{
    gamal = 0;
    gamar = 0;

    lrtl = 0;
    lrtr = 0;

    tcl = 0;
    tcr = 0;

    pcl = 0;
    pcr = 0;

    tl   = 0;
    tr   = 0;

    ql   = 0;
    qr   = 0;
    flux = 0;

    qlc  = 0;    //! GMRESPassQC
    qrc  = 0;    //! GMRESPassQC

    xfn  = 0;
    yfn  = 0;
    zfn  = 0;
    area = 0;
    vgn  = 0;

    //! GMRESnolim
    xfc  = 0;
    yfc  = 0;
    zfc  = 0;
    xcc  = 0;
    ycc  = 0;
    zcc  = 0;
    dqdx = 0;
    dqdy = 0;
    dqdz = 0;
    nsx  = 0;
    nsy  = 0;
    nsz  = 0;
    ns   = 0;
    vol  = 0;
    prim = 0;
    gama = 0;
    neighborCells = 0;
    neighborFaces = 0;
    neighborLR    = 0;
    qlsign        = 0;
    qrsign        = 0;

    unstructBCSet = 0;    //! GMRES3D

}

InviscidSchemeParameter::~InviscidSchemeParameter()
{
    //! Do not free, since these pointers are not allocated by this class.
}


}

