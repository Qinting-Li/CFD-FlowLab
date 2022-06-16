#include "Flux_Inviscid.h"
#include "Constants.h"
#include "TK_Exit.h"
#include "Flux_RoeEntropyCorrection.h"
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
    //! Bell 20130418 add.
    else
    {
        TK_Exit::ExceptionExit("Error: this inv-scheme is not exist !\n", true);
    }
    //*/
    
    return scheme_id;
}

INVScheme GetINVScheme(int scheme_id)
{
    INVScheme inv_scheme = NULL;
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
    case ISCHEME_STEGER:
        inv_scheme = Steger_Scheme;
        break;
    case ISCHEME_HLLE:
        inv_scheme = Hlle_Scheme;
        break;
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

    RDouble refMachNumber = invSchemePara->GetMachNumber();
    RDouble entrFixCoefLeft = invSchemePara->GetLeftEntropyFixCoefficient();
    RDouble entrFixCoefRight = invSchemePara->GetRightEntropyFixCoefficient();
    RDouble **qL = invSchemePara->GetLeftQ();
    RDouble **qR = invSchemePara->GetRightQ();
    RDouble *gamaL = invSchemePara->GetLeftGama();
    RDouble *gamaR = invSchemePara->GetRightGama();
    RDouble *lrtL = face_proxy->GetPressureCoefficientL();
    RDouble *lrtR = face_proxy->GetPressureCoefficientR();
    RDouble *xfn   = invSchemePara->GetFaceNormalX();
    RDouble *yfn   = invSchemePara->GetFaceNormalY();
    RDouble *zfn   = invSchemePara->GetFaceNormalZ();
    RDouble *vgn   = invSchemePara->GetFaceVelocity();
    RDouble **flux = invSchemePara->GetFlux();
    RDouble **tl = invSchemePara->GetLeftTemperature();
    RDouble **tr = invSchemePara->GetRightTemperature();
    RDouble mach2 = refMachNumber * refMachNumber;

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
            vnew  = half * vn * (one + beta);
            cnew  = half * sqrt(((beta - one) * vn) * ((beta - one) * vn) + four * beta * c2);
            cnew2  = cnew * cnew;
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

        Flux_RoeEntropyCorrection(eigv1, eigv2, eigv3, lrtL[iFace], lrtR[iFace], absVel, vn, cm, RoeEntropyFixMethod, entrFixCoefLeft, entrFixCoefRight);

        RDouble xi1 = 0.0, xi2 = 0.0, xi3 = 0.0, xi4 = 0.0;
        if (nPrecondition == 0) 
        {
            xi1 = (two * eigv1 - eigv2 - eigv3) / (two * c2);
            xi2 = (eigv2 - eigv3)/(two * cm);
        }
        else
        {
            xi1 = (two * eigv1 - eigv2 - eigv3) / (two * cnew2);
            xi2 = (eigv2 - eigv3)/(two * cnew);
            alpha = half * (one - beta);
            xi3 = two * alpha * eigv1 / cnew2;
            xi4 =  alpha * vn * xi2 / cnew2;
        }

        RDouble dc = theta * dq[IR] - xsn * dq[IRU] - ysn * dq[IRV] - zsn * dq[IRW];

        RDouble c2dc = 0;
        if (nPrecondition == 0)
        {
            c2dc = c2 * dc;
        }
        else
        {
            Ur2  = beta * c2;
            c2dc = cnew2 * dc;
        }

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
            flux[IR ][iFace] -= half * (eigv1 * dq[IR ]                 -               dh * cnew2 / Ur2  * xi1               -               dc  * xi2 +               dh * cnew2 / Ur2 * xi3 -                                dh * cnew2 / Ur2  * xi4);
            flux[IRU][iFace] -= half * (eigv1 * dq[IRU] + (xsn   * c2dc - uRoeAverage * dh * cnew2 / Ur2) * xi1 + (xsn   * dh - uRoeAverage * dc) * xi2 + uRoeAverage * dh * cnew2 / Ur2 * xi3 - (xsn    * c2dc + uRoeAverage * dh * cnew2 / Ur2) * xi4);
            flux[IRV][iFace] -= half * (eigv1 * dq[IRV] + (ysn   * c2dc - vRoeAverage * dh * cnew2 / Ur2) * xi1 + (ysn   * dh - vRoeAverage * dc) * xi2 + vRoeAverage * dh * cnew2 / Ur2 * xi3 - (ysn    * c2dc + vRoeAverage * dh * cnew2 / Ur2) * xi4);
            flux[IRW][iFace] -= half * (eigv1 * dq[IRW] + (zsn   * c2dc - wRoeAverage * dh * cnew2 / Ur2) * xi1 + (zsn   * dh - wRoeAverage * dc) * xi2 + wRoeAverage * dh * cnew2 / Ur2 * xi3 - (zsn    * c2dc + wRoeAverage * dh * cnew2 / Ur2) * xi4);
            flux[IRE][iFace] -= half * (eigv1 * dq[IRE] + (theta * c2dc - hRoeAverage * dh * cnew2 / Ur2) * xi1 + (theta * dh - hRoeAverage * dc) * xi2 + hRoeAverage * dh * cnew2 / Ur2 * xi3 - (theta * c2dc + hRoeAverage * dh * cnew2 / Ur2) * xi4);

            for (int iLaminar = nm; iLaminar < nLaminar; ++ iLaminar)
            {
                prim[iLaminar] = half * (primL[iLaminar] + primR[iLaminar]);
                flux[iLaminar][iFace] -= half * (eigv1 * dq[iLaminar] - prim[iLaminar] * (dh * cnew2 / Ur2 * (xi1 - xi3 + xi4) + dc * xi2));
            }
            //! Compute vibration-electron energy term.
            for (int it = 1; it < nTemperatureModel; ++ it)
            {
                prim[nLaminar + it] = half * (primL[nLaminar + it] + primR[nLaminar + it]);
                flux[nLaminar + it][iFace] -= half * (eigv1 * dq[nLaminar + it] - prim[nLaminar + it] * (dh * cnew2 / Ur2 * (xi1 - xi3 + xi4) + dc * xi2));
            }
        }

    }

    delete []prim;    prim = NULL;
    delete []primL;    primL = NULL;
    delete []primR;    primR = NULL;
    delete []qConserveL;    qConserveL = NULL;
    delete []qConserveR;    qConserveR = NULL;
    delete []dq;    dq = NULL;
}

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

    RDouble *lrtL = faceProxy->GetPressureCoefficientL();
    RDouble *lrtR = faceProxy->GetPressureCoefficientR();
    RDouble **tl = inviscidSchemePara->GetLeftTemperature();
    RDouble **tr = inviscidSchemePara->GetRightTemperature();
    RDouble mach2 = refMachNumber * refMachNumber;

    using namespace IDX;
    using namespace GAS_SPACE;

    RDouble *dq         = new RDouble[nEquation];

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

        RDouble cnew2, cnew, vnew, beta;
        RDouble eigv1, eigv2, eigv3;
        eigv1 = ABS(vn);
        if (nPrecondition == 0)
        {
            eigv2 = ABS(vn + cm);
            eigv3 = ABS(vn - cm);
        }
        else
        {
            beta = MIN(MAX(vSquareRoeAverage / c2, three * mach2), one);
            vnew  = half * vn * (one + beta);
            cnew  = half * sqrt(((beta - one) * vn) * ((beta - one) * vn) + four * beta * c2);
            cnew2  = cnew * cnew;
            eigv2 = ABS(vnew + cnew);
            eigv3 = ABS(vnew - cnew);
        }

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



        Flux_RoeEntropyCorrection(eigv1, eigv2, eigv3, lrtL[iFace], lrtR[iFace], absVel, vn, cm, RoeEntropyFixMethod, entrFixCoefLeft, entrFixCoefRight);

        RDouble eigv1Term, eigv2Term, eigv3Term;
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

    delete []primL;
    delete []primR;
    delete []dq;
}

//clz begin
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

            maxEigv = MAX(absEigvr2, absEigvr3);
            tmp0     = maxEigv * entrFixCoefLeft ;
            tmp1     = maxEigv * entrFixCoefRight;
            absEigvr1   = MAX(tmp0, absEigvr1);
            absEigvr2   = MAX(tmp1, absEigvr2);
            absEigvr3   = MAX(tmp1, absEigvr3);
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

    delete []fl;
    delete []fr;
    delete []priml;
    delete []primr;
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

    RDouble **ql = invSchemePara->GetLeftQ();
    RDouble **qr = invSchemePara->GetRightQ();
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

        //! It is convenient for later package fluxsolver,and adding variables.
        RDouble xsn = xfn[iLen];
        RDouble ysn = yfn[iLen];
        RDouble zsn = zfn[iLen];
        RDouble vb =  vgn[iLen];
        RDouble gmL = gamal[iLen];
        RDouble gmR = gamar[iLen];

        RDouble v2l = ul * ul + vl * vl + wl * wl;
        RDouble v2r = ur * ur + vr * vr + wr * wr;

        RDouble hintLeft  = (gmL / (gmL - one)) * (pl / rl);
        RDouble hintRight = (gmR / (gmR - one)) * (pr / rr);

        RDouble peL = 0.0, peR = 0.0;

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
                peL = gas->GetElectronPressure(priml, tl[2][iLen]);
                peR = gas->GetElectronPressure(primr, tr[2][iLen]);
            }
        }

        RDouble hl  = hintLeft  + half * v2l;
        RDouble hr  = hintRight + half * v2r;

        //! Full flux.
        RDouble vnl  = xsn * ul + ysn * vl + zsn * wl - vb;
        RDouble vnr  = xsn * ur + ysn * vr + zsn * wr - vb;
        RDouble rvnl = rl * vnl;
        RDouble rvnr = rr * vnr;

        RDouble c2l   = gmL * pl / rl;
        RDouble cl    = sqrt(ABS(c2l));
        RDouble machl = vnl / cl;
        RDouble mach2l = machl * machl;

        RDouble c2r   = gmR * pr / rr;
        RDouble cr    = sqrt(ABS(c2r));
        RDouble machr = vnr / cr;
        RDouble mach2r = machr * machr;

        RDouble mach2m=MAX(mach2l,mach2r);
        RDouble coef =sqrt(MIN(MAX(mach2m, 3.0 * mach2), 1.0));

        RDouble fel = 0.0, fer = 0.0;
        if (machl > one)
        {
            fl[IR ] = rvnl;
            fl[IRU] = rvnl * ul + xsn * pl;
            fl[IRV] = rvnl * vl + ysn * pl;
            fl[IRW] = rvnl * wl + zsn * pl;
            fl[IRE] = rvnl * hl + vb  * pl;
            fel = peL * (vnl + vb);
        }
        else if (machl < - one)
        {
            fl[IR ] = zero;
            fl[IRU] = zero;
            fl[IRV] = zero;
            fl[IRW] = zero;
            fl[IRE] = zero;
            fel = zero;
        }
        else
        {
            if (nPrecondition == 0)
            {
                RDouble fmass = 0.25 * rl * cl * (machl + 1.0) * (machl + 1.0);
                RDouble tmp = (- vnl + 2.0 * cl) / gmL;
                fl[IR ] = fmass;
                fl[IRU] = fmass * (ul + xsn * tmp);
                fl[IRV] = fmass * (vl + ysn * tmp);
                fl[IRW] = fmass * (wl + zsn * tmp);
                fl[IRE] = fmass * (hl + vb  * tmp);
            }
            else
            {
                RDouble fmass = 0.25 * rl * cl * (machl * machl + 1.0);
                RDouble tmp0  = - vnl / gmL;
                RDouble tmp1  = rl * cl / gmL;
                fl[IR ] = half * rvnl + coef *  fmass;
                fl[IRU] = half * (rvnl * ul + xsn * pl) + coef * (fmass * (ul + xsn * tmp0) + ul * tmp1);
                fl[IRV] = half * (rvnl * vl + ysn * pl) + coef * (fmass * (vl + ysn * tmp0) + vl * tmp1);
                fl[IRW] = half * (rvnl * wl + zsn * pl) + coef * (fmass * (wl + zsn * tmp0) + wl * tmp1);
                fl[IRE] = half * (rvnl * hl + vb  * pl) + coef * (fmass * (hl + vb  * tmp0) + hl * tmp1);
                fel = half * peL * (vnl + vb) + coef * (fmass * (peL / rl + vb  * tmp0) + peL / rl * tmp1);
            }
        }

        for (int iLaminar = nm; iLaminar < nLaminar; ++ iLaminar)
        {
            fl[iLaminar] = priml[iLaminar] * fl[IR];
        }

        //! Right contribution.
        if (machr < - one)
        {
            fr[IR ] = rvnr                ;
            fr[IRU] = rvnr * ur + xsn * pr;
            fr[IRV] = rvnr * vr + ysn * pr;
            fr[IRW] = rvnr * wr + zsn * pr;
            fr[IRE] = rvnr * hr + vb  * pr;
            fer = peR * (vnr + vb);
        }
        else if (machr > one)
        {
            fr[IR ] = zero;
            fr[IRU] = zero;
            fr[IRV] = zero;
            fr[IRW] = zero;
            fr[IRE] = zero;
            fer = zero;
        }
        else
        {
            if (nPrecondition == 0)
            {
                RDouble fmass = - 0.25 * rr * cr * (machr - 1.0) * (machr - 1.0);
                RDouble tmp = (- vnr - 2.0 * cr) / gmR;
                fr[IR ] = fmass;
                fr[IRU] = fmass * (ur + xsn * tmp);
                fr[IRV] = fmass * (vr + ysn * tmp);
                fr[IRW] = fmass * (wr + zsn * tmp);
                fr[IRE] = fmass * (hr + vb  * tmp);
            }
            else
            {
                RDouble fmass = 0.25 * rr * cr * (machr * machr + 1.0);
                RDouble tmp0  = - vnr / gmR;
                RDouble tmp1  = rr * cr / gmR;
                fr[IR ] = half * rvnr  - coef *  fmass;
                fr[IRU] = half * (rvnr * ur + xsn * pr) - coef * (fmass * (ur + xsn * tmp0) + ur * tmp1);
                fr[IRV] = half * (rvnr * vr + ysn * pr) - coef * (fmass * (vr + ysn * tmp0) + vr * tmp1);
                fr[IRW] = half * (rvnr * wr + zsn * pr) - coef * (fmass * (wr + zsn * tmp0) + wr * tmp1);
                fr[IRE] = half * (rvnr * hr + vb  * pr) - coef * (fmass * (hr + vb  * tmp0) + hr * tmp1);
                fer = half * peR * (vnr + vb) - coef * (fmass * (peR / rr + vb  * tmp0) + peR / rr * tmp1);
            }
        }

        for (int iLaminar = nm; iLaminar < nLaminar; ++ iLaminar)
        {
            fr[iLaminar] = primr[iLaminar] * fr[IR];
        }

        //! Compute vibration-electron energy term.
        for (int it = 1; it < nTemperatureModel; ++ it)
        {
            fl[nLaminar + it] = priml[nLaminar + it] * fl[IR];
            fr[nLaminar + it] = primr[nLaminar + it] * fr[IR];
        }
        if (nChemical == 1)
        {
            fl[nLaminar + nTemperatureModel - 1] += fel;
            fr[nLaminar + nTemperatureModel - 1] += fer;
        }

        for (int iLaminar = 0; iLaminar < nEquation; ++ iLaminar)
        {
            flux[iLaminar][iLen] = fl[iLaminar] + fr[iLaminar];
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

    delete []fl;
    delete []fr;
    delete []priml;
    delete []primr;
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
        ml    =  vnl * ocm;
        mla   =  ABS(ml);

        mr    =  vnr * ocm;
        mra   =  ABS(mr);

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

    delete []priml;
    delete []primr;

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
    delete []priml;
    delete []primr;
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

        RDouble vnl  = xfn[iLen] * ul + yfn[iLen] * vl + zfn[iLen] * wl - vgn[iLen];
        RDouble vnr  = xfn[iLen] * ur + yfn[iLen] * vr + zfn[iLen] * wr - vgn[iLen];

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

        RDouble ml = abs(vnl / al);
        RDouble mr = abs(vnr / ar);

        RDouble c12, ocm;

        /*if( ( ml > 1.0 && mr < 1.0 ) || ( mr > 1.0 && ml < 1.0 ))
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
        /*}*/

        ocm = 1.0 / c12;

        ml  =  vnl * ocm;
        mr  =  vnr * ocm;

        RDouble fm4ml = FMmpn(ml, zero,  one);
        RDouble fm4mr = FMmpn(mr, zero, -one);

        RDouble fp5ml = FMppn(ml, zero,  one);
        RDouble fp5mr = FMppn(mr, zero, -one);

        RDouble m12  = fm4ml + fm4mr;
        RDouble p12  = fp5ml * pl + fp5mr * pr;
        RDouble pe12 = fp5ml * peL + fp5mr * peR;

        RDouble minplr = MIN(pl / (pr + SMALL), pr / (pl + SMALL));
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
    delete []priml;
    delete []primr;
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

    delete []qcsvl;
    delete []qcsvr;
    delete []fl;
    delete []fr;
    delete []priml;
    delete []primr;
}

RDouble FMSplit1(const RDouble &mach, const RDouble &ipn)
{
    return half * (mach + ipn * ABS(mach));
}

RDouble FMSplit2(const RDouble &mach, const RDouble &ipn)
{
    return ipn * fourth * SQR(mach + ipn);
}

RDouble FMmpn(const RDouble &mach, const RDouble &beta, const RDouble &ipn)
{
    if (ABS(mach) >= one)
    {
        return FMSplit1(mach,ipn);
    }
    return FMSplit2(mach,ipn);
}

RDouble FMppn(const RDouble &mach, const RDouble &beta, const RDouble &ipn)
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
    RDouble hRight  = gamaRight / (gamaRight - 1.0) * pRight / rhoRight + half * (uRight * uRight + vRight * vRight + wRight * wRight);

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

    tl   = 0;
    tr   = 0;

    ql   = 0;
    qr   = 0;
    flux = 0;

    xfn  = 0;
    yfn  = 0;
    zfn  = 0;
    area = 0;
    vgn  = 0;
}

InviscidSchemeParameter::~InviscidSchemeParameter()
{
    //! Do not free, since these pointers are not allocated by this class.
}




}

