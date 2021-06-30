#include "Jacobian.h"
#include "Constants.h"
#include "Math_BasisFunction.h"
#include "GlobalDataBase.h"
#include "Gas.h"

using namespace std;
namespace PHSPACE
{
//! The fuction computes the formula (M +- R) * dQ, M and R denote the Jacobian matrix and spectrum radius of inviscid flux,respectively,\n
//! dQ is the increment of conservative vector.
void MXDQ_STD(RDouble *prim, RDouble nxs, RDouble nys, RDouble nzs, RDouble area, RDouble vgn, RDouble gama, RDouble *dq, RDouble *f, 
              int nm, int nLaminar, RDouble radius, int ipn)
{
    //! The relevant variables related to vgn has modified by LiPeng on Feb 19, 2019
    using namespace IDX;
    using namespace GAS_SPACE;

    RDouble rm = prim[IDX::IR];
    RDouble um = prim[IDX::IU];
    RDouble vm = prim[IDX::IV];
    RDouble wm = prim[IDX::IW];
    RDouble pm = prim[IDX::IP];

    RDouble v2 = um * um + vm * vm + wm * wm;
    RDouble c2 = ABS(gama * pm / rm);
    RDouble cm = sqrt(c2);

    RDouble lmd1, lmd2, lmd3, x1, x2, cits;

    //! The variable cits equals to Ur=U-Ub, which can be referred to the formula (A.5) in appendix A of the PHengLEI Theory manual.
    cits = nxs * um + nys * vm + nzs * wm - vgn;

    RDouble ctgm = cits;
    RDouble cmgm = cm;

    lmd1 = ctgm;
    lmd2 = ctgm + cmgm;
    lmd3 = ctgm - cmgm;

    RDouble ipnrad = ipn * radius / area;

    lmd1 = half * (lmd1 + ipnrad);
    lmd2 = half * (lmd2 + ipnrad);
    lmd3 = half * (lmd3 + ipnrad);

    //!This term is actually equal to 0.0 
    x1 = (lmd1 + lmd1 - lmd2 - lmd3) / (c2 + c2);
    //!This term is actually equal to 0.5
    x2 = (lmd2 - lmd3) / (cm + cm);

    RDouble dc, dh, c2dc;
    //! dc=-b1,b1 is referred to the formula (A.7) in appendix A of the PHengLEI Theory manual.
    //! dc   = cits * dq[IR] - nxs * dq[IRU] - nys * dq[IRV] - nzs * dq[IRW];//!the previous form
    dc   = (cits + vgn) * dq[IDX::IR] - nxs * dq[IDX::IRU] - nys * dq[IDX::IRV] - nzs * dq[IDX::IRW];//!Modified by LiPeng on Feb 21, 2019.
    c2dc = c2 * dc;

    //! Compute total enthalpy hm and coefficient dh,dh=b2,see formula (A.7) and (A.8) in appendix A of the PHengLEI Theory manual.
    RDouble hm;
    gas->ComputeTotalEnthalpyAndDH(prim, gama, dq, hm, dh);

    //!This term is actually equal to (0.5 * dc)
    RDouble dhx1_dcx2 =  dh   * x1 + dc * x2;
    //!This term is actually equal to (0.5 * dh)
    RDouble c2dcx1_dhx2 = c2dc * x1 + dh * x2;

    //! The following expressions can be referred to formula (A.6) in appendix A of the PHengLEI Theory manual.
    f[IDX::IR ] = lmd1 * dq[IR ] - dhx1_dcx2;
    f[IDX::IRU] = lmd1 * dq[IRU] - um * dhx1_dcx2 + nxs  * c2dcx1_dhx2;
    f[IDX::IRV] = lmd1 * dq[IRV] - vm * dhx1_dcx2 + nys  * c2dcx1_dhx2;
    f[IDX::IRW] = lmd1 * dq[IRW] - wm * dhx1_dcx2 + nzs  * c2dcx1_dhx2;

    //! The previous form is not suitable for dynamic mesh.
    //! f[IRE] = lmd1 * dq[IRE] - hm * dhx1_dcx2 + cits * c2dcx1_dhx2;
    //! Modified by LiPeng on Feb 21, 2019.
    f[IDX::IRE] = lmd1 * dq[IDX::IRE] - hm * dhx1_dcx2 + (cits + vgn) * c2dcx1_dhx2;
    for (int iLaminar = nm; iLaminar < nLaminar; ++ iLaminar)
    {
        f[iLaminar] = lmd1 * dq[iLaminar] - prim[iLaminar] * dhx1_dcx2;
    }

    for (int iLaminar = 0; iLaminar < nLaminar; ++ iLaminar) //!Modified by LiPeng on Feb. 25, 2019
    {
        f[iLaminar] *= area;
    }
}

void MXDQ_STD(RDouble *prim, RDouble *gykb, RDouble *dq, RDouble *f, RDouble gama, int nm, int nLaminar, RDouble radius, int ipn)
{
    using namespace IDX;
    using namespace GAS_SPACE;

    RDouble rm  = prim[IR];
    RDouble um  = prim[IU];
    RDouble vm  = prim[IV];
    RDouble wm  = prim[IW];
    RDouble pm  = prim[IP];

    RDouble v2 = um * um + vm * vm + wm * wm;
    RDouble c2 = ABS(gama * pm / rm);
    RDouble cm = sqrt(c2);

    RDouble lmd1, lmd2, lmd3, x1, x2, cits;

    //! the variable cits equals to Ur=U-Ub, which can be referred to the formula (A.5) in appendix A of the PHengLEI Theory manual.
    cits = gykb[0] * um + gykb[1] * vm + gykb[2] * wm - gykb[4]; 

    RDouble ctgm = cits;
    RDouble cmgm = cm  ;

    lmd1 = ctgm;
    lmd2 = ctgm + cmgm;
    lmd3 = ctgm - cmgm;

    RDouble ipnrad = ipn * radius / gykb[3];
    //RDouble ipnrad = ipn * radius;  // org

    lmd1 = half * (lmd1 + ipnrad);
    lmd2 = half * (lmd2 + ipnrad);
    lmd3 = half * (lmd3 + ipnrad);

    x1 = (lmd1 + lmd1 - lmd2 - lmd3) / (c2 + c2);	//this term is actually equal to 0.0 
    x2 = (lmd2 - lmd3) / (cm + cm);					//this term is actually equal to 0.5

    RDouble dc,dh,c2dc;
    //! dc=-b1, b1 is referred to the formula (A.7) in appendix A of the PHengLEI Theory manual.
    //! the previous form
    //  dc   = cits * dq[IR] - nxs * dq[IRU] - nys * dq[IRV] - nzs * dq[IRW];
    //! modified by LiPeng on Feb 21, 2019.
    dc   = (cits + gykb[4]) * dq[IR] - gykb[0] * dq[IRU] - gykb[1] * dq[IRV] - gykb[2] * dq[IRW];
    c2dc = c2 * dc;

    //! compute total enthalpy hm and coefficient dh,dh=b2,see formula (A.7) and (A.8) in appendix A of the PHengLEI Theory manual.
    RDouble hm;
    gas->ComputeTotalEnthalpyAndDH(prim, gama, dq, hm, dh);

    RDouble dhx1_dcx2   =  dh   * x1 + dc * x2;//this term is actually equal to (0.5 * dc)
    RDouble c2dcx1_dhx2 =  c2dc * x1 + dh * x2;//this term is actually equal to (0.5 * dh)

    //the following expressions can be referred to formula (A.6) in appendix A of the PHengLEI Theory manual.
    f[IR ] = lmd1 * dq[IR ] - dhx1_dcx2;
    f[IRU] = lmd1 * dq[IRU] - um * dhx1_dcx2 + gykb[0]  * c2dcx1_dhx2;
    f[IRV] = lmd1 * dq[IRV] - vm * dhx1_dcx2 + gykb[1]  * c2dcx1_dhx2;
    f[IRW] = lmd1 * dq[IRW] - wm * dhx1_dcx2 + gykb[2]  * c2dcx1_dhx2;
    //! the previous form is not suitable for dynamic mesh.
    //! f[IRE] = lmd1 * dq[IRE] - hm * dhx1_dcx2 + cits * c2dcx1_dhx2;
    //! modified by LiPeng on Feb 21, 2019.
    f[IRE] = lmd1 * dq[IRE] - hm * dhx1_dcx2 + (cits + gykb[4]) * c2dcx1_dhx2;
    for (int iLaminar = nm; iLaminar < nLaminar; ++ iLaminar)
    {
        f[iLaminar] = lmd1 * dq[iLaminar] - prim[iLaminar] * dhx1_dcx2;
    }

    for (int iLaminar = 0; iLaminar < nLaminar; ++ iLaminar) //modified by LiPeng on Feb. 25, 2019
    {
        f[iLaminar] *= gykb[3];
    }
}

void ComputePreconMatrix(int sign, RDouble **GG, RDouble *prim, RDouble tm, RDouble gama)
{
    using namespace IDX;
    using namespace GAS_SPACE;

    RDouble rm, um, vm, wm, pm, Hm, mach2, v2, c2, rhot, Cp, beta, beta_new, Ep, theta, psi, gama1;    
    RDouble refMachNumber;
    GlobalDataBase::GetData("refMachNumber", &refMachNumber, PHDOUBLE, 1);

    switch(sign) 
    {
    case 0:  // precondition matrix 

        mach2 = refMachNumber * refMachNumber;
        gama1 = gama - 1.0;
        rm = prim[IR];
        um = prim[IU];
        vm = prim[IV];
        wm = prim[IW];
        pm = prim[IP];

        c2 = gama * pm / rm;
        v2 = um * um + vm * vm + wm * wm;

        beta = MIN(MAX(v2 / c2, 3.0 * mach2),1.0);
        Hm = gama / gama1 * pm / rm + 0.5 * v2;

        rhot = -rm / tm;
        Cp = 1.0/ (gama1 * mach2);
        beta_new = beta / (1 + gama1 * beta);
        Ep = 1/ (beta_new * c2);
        theta = mach2 / (Ep * tm);

        GG[IR][0] = theta;
        GG[IR][1] = 0.0;
        GG[IR][2] = 0.0;
        GG[IR][3] = 0.0;
        GG[IR][4] = rhot;

        GG[IU][0] = um * theta;
        GG[IU][1] = rm;
        GG[IU][2] = 0.0;
        GG[IU][3] = 0.0;
        GG[IU][4] = um * rhot;

        GG[IV][0] = vm * theta;
        GG[IV][1] = 0.0;
        GG[IV][2] = rm;
        GG[IV][3] = 0.0;
        GG[IV][4] = vm * rhot;

        GG[IW][0] = wm * theta;
        GG[IW][1] = 0.0;
        GG[IW][2] = 0.0;
        GG[IW][3] = rm;
        GG[IW][4] = wm * rhot;

        GG[IP][0] = Hm * theta - 1.0;
        GG[IP][1] = um * rm;
        GG[IP][2] = vm * rm;
        GG[IP][3] = wm * rm;
        GG[IP][4] = Hm * rhot + rm * Cp;

        break;

    case 1: // M * inverse precondition matrix

        mach2 = refMachNumber * refMachNumber;
        gama1 = gama - 1.0;
        rm = prim[IR];
        um = prim[IU];
        vm = prim[IV];
        wm = prim[IW];
        pm = prim[IP];

        c2 = gama * pm / rm;
        v2 = um * um + vm * vm + wm * wm;
        beta = MIN(MAX(v2 / c2, 3.0 * mach2), 1.0);
        Hm = gama / gama1 * pm / rm + 0.5 * v2;
        psi  = (beta - 1.0) * gama1 * mach2 / tm;

        GG[IR][0] = 1.0 + 0.5 * v2 * psi;
        GG[IR][1] = - um * psi;
        GG[IR][2] = - vm * psi;
        GG[IR][3] = - wm * psi;
        GG[IR][4] = psi;

        GG[IRU][0] = 0.5 * um * v2 * psi;
        GG[IRU][1] = 1.0 - um * um * psi;
        GG[IRU][2] = -um * vm * psi;
        GG[IRU][3] = -um * wm * psi;
        GG[IRU][4] = um * psi;

        GG[IRV][0] = 0.5 * vm * v2 * psi;
        GG[IRV][1] = -vm * um * psi;
        GG[IRV][2] = 1.0 - vm * vm * psi;
        GG[IRV][3] = -vm * wm * psi;
        GG[IRV][4] = vm * psi;

        GG[IRW][0] = 0.5 * wm * v2 * psi;
        GG[IRW][1] = -wm * um * psi;
        GG[IRW][2] = -wm * vm * psi;
        GG[IRW][3] = 1.0 - wm * wm * psi;
        GG[IRW][4] = wm * psi;

        GG[IRE][0] = 0.5 * Hm * v2 * psi;
        GG[IRE][1] = - um * Hm * psi;
        GG[IRE][2] = - vm * Hm * psi;
        GG[IRE][3] = - wm * Hm * psi;
        GG[IRE][4] = 1.0 + Hm * psi;

    }

}

void PreconMatrixDF(RDouble *df, RDouble **GG, int nm)
{
    using namespace IDX;

    RDouble * df_trans = new RDouble [nm];

    for (int m = 0; m < nm; ++ m)
    {
        df_trans[m]  = 0;
    }

    for (int m = 0; m < nm; ++ m)
    {
        df_trans[IR]  += df[m] * GG[IR][m];
        df_trans[IRU] += df[m] * GG[IRU][m];
        df_trans[IRV] += df[m] * GG[IRV][m];
        df_trans[IRW] += df[m] * GG[IRW][m];
        df_trans[IRE] += df[m] * GG[IRE][m];
    }

    for (int m = 0; m < nm; ++ m)
    {
        df[m] = df_trans[m];
    }

    delete []df_trans;
}

void PreconMatrixRes(RDouble **dq, RDouble **GG, int nm, int iCell)
{
    using namespace IDX;
    RDouble * dqTrans = new RDouble[nm];
    for (int m = 0; m < nm; ++ m)
    {
        dqTrans[m]  = 0;
    }

    for (int m = 0; m < nm; ++ m)
    {
        dqTrans[IR] += dq[m][iCell] * GG[IR ][m];
        dqTrans[IU] += dq[m][iCell] * GG[IRU][m];
        dqTrans[IV] += dq[m][iCell] * GG[IRV][m];
        dqTrans[IW] += dq[m][iCell] * GG[IRW][m];
        dqTrans[IP] += dq[m][iCell] * GG[IRE][m];
    }

    for (int m = 0; m < nm; ++ m)
    {
        dq[m][iCell] = dqTrans[m];
    }

    delete []dqTrans;
}

void ChemicalMXDQ(RDouble *primitiveVars, RDouble *temperature, RDouble nx, RDouble ny, RDouble nz, RDouble faceArea, RDouble moveVelocity,
                  RDouble *deltaQ, RDouble *flux, int nNSEquation, int nLaminar, int nTemperatureModel, RDouble radius, int iSign)
{
    using namespace IDX;
    using namespace GAS_SPACE;
    // The number of species.
    int nSpecies = nLaminar - nNSEquation + 1;
    int nElectronIndex = GlobalDataBase::GetIntParaFromDB("nElectronIndex");
    RDouble alpha = 0.0, beta = 0.0, gama = 0.0, *speciesBeta;
    speciesBeta = new RDouble[nSpecies];

    RDouble uVelocity = primitiveVars[IU];
    RDouble vVelocity = primitiveVars[IV];
    RDouble wVelocity = primitiveVars[IW];
    RDouble squareVelocity, absoluteVelocity, relativeVelocity;

    // The absolute velocity.
    absoluteVelocity = nx * uVelocity + ny * vVelocity + nz * wVelocity;
    // The relative velocity.
    relativeVelocity = absoluteVelocity - moveVelocity;
    // Compute squared velocity.
    squareVelocity = uVelocity * uVelocity + vVelocity * vVelocity + wVelocity * wVelocity;

    // Obtain the temperatures.
    RDouble transRotationTemperature = temperature[ITT];
    RDouble vibrationTemperature = transRotationTemperature;
    RDouble electronTemperature = transRotationTemperature;
    if (nTemperatureModel > 1)
    {
        vibrationTemperature = temperature[ITV];
        electronTemperature = vibrationTemperature;
        if (nTemperatureModel == 3)
        {
            electronTemperature = temperature[ITE];
        }
    }

    // Compute the total enthalpy.
    RDouble totalEnthalpy = gas->GetMixedGasEnthalpy(primitiveVars, transRotationTemperature, vibrationTemperature, electronTemperature);
    totalEnthalpy += 0.5 * squareVelocity;

    // Compute the coefficient b1 and b2 (see formula A.5 in PHengLEI Theory manual).
    RDouble b1 = nx * deltaQ[IRU] + ny * deltaQ[IRV] + nz * deltaQ[IRW] - absoluteVelocity * deltaQ[IR];
    RDouble dh = uVelocity * deltaQ[IRU] + vVelocity * deltaQ[IRV] + wVelocity * deltaQ[IRW] - deltaQ[IRE];
    RDouble b2 = 0.0, b3 = 0.0, electronPressure = 0.0;

    if (nTemperatureModel == 1)
    {
        // Compute the coefficients in the chemical reaction.
        gas->ComputeCoefficientInMXDQ(primitiveVars, alpha, beta, speciesBeta);

        b2 = beta * deltaQ[IR] - alpha * dh;
        for (int s = 0; s < nSpecies - 1; ++ s)
        {
            b2 += speciesBeta[s] * deltaQ[nNSEquation + s];
        }
    }
    else // Multi-temperature model.
    {
        // Compute the pressure of electron.
        electronPressure = gas->GetElectronPressure(primitiveVars, electronTemperature);
        RDouble *speciesGama = new RDouble[nSpecies];

        // Compute partial derivatives of pressure an electron pressure to conservative variables.
        gas->GetMultiTemperatureModelPartialDerivatives(primitiveVars, transRotationTemperature, electronTemperature,
                                                                   alpha, beta, gama, speciesBeta, speciesGama);

        if (nTemperatureModel == 2)
        {
            dh += deltaQ[nLaminar + 1];
            b2 += gama * deltaQ[nLaminar + 1];
            b3 += gama * deltaQ[nLaminar + 1];
        }
        else if (nTemperatureModel == 3)
        {
            b2 += gama * deltaQ[nLaminar + 2];
            b3 += gama * deltaQ[nLaminar + 2];
            dh += deltaQ[nLaminar + 1] + deltaQ[nLaminar + 2];
        }
        b2 += beta * deltaQ[IR] - alpha * dh;

        for (int iSpices = 0; iSpices < nSpecies - 1; ++ iSpices)
        {
            b2 += speciesBeta[iSpices] * deltaQ[nNSEquation + iSpices];
            b3 += speciesGama[iSpices] * deltaQ[nNSEquation + iSpices];
        }

        delete []speciesGama;
    }

    // Comoute the fluxes called M*dQ.
    flux[IR ] = relativeVelocity * deltaQ[IR ] + b1;
    flux[IRU] = relativeVelocity * deltaQ[IRU] + b1 * uVelocity     + b2 * nx;
    flux[IRV] = relativeVelocity * deltaQ[IRV] + b1 * vVelocity     + b2 * ny;
    flux[IRW] = relativeVelocity * deltaQ[IRW] + b1 * wVelocity     + b2 * nz;
    flux[IRE] = relativeVelocity * deltaQ[IRE] + b1 * totalEnthalpy + b2 * absoluteVelocity;
    for (int m = nNSEquation; m < nLaminar; ++ m)
    {
        flux[m] = relativeVelocity * deltaQ[m] + b1 * primitiveVars[m];
    }

    // Modify the fluxes of multi-temperature model.
    if (nTemperatureModel == 2) 
    {
        flux[nLaminar + 1] = relativeVelocity * deltaQ[nLaminar + 1] + b1 * (primitiveVars[nLaminar + 1] + electronPressure / primitiveVars[IR]) + b3 * absoluteVelocity;
    }
    else if (nTemperatureModel == 3)
    {
        flux[nLaminar + 1] = relativeVelocity * deltaQ[nLaminar + 1] + b1 * primitiveVars[nLaminar + 1];
        flux[nLaminar + 2] = relativeVelocity * deltaQ[nLaminar + 2] + b1 * (primitiveVars[nLaminar + 2] + electronPressure / primitiveVars[IR]) + b3 * absoluteVelocity;
    }

    // To add rM * dQ that is the product of the invscid spectrum radius rM and increment dQ.
    RDouble spectrumRadiusTerm = iSign * radius / faceArea;
    for (int m = 0; m < nLaminar; ++ m)
    {
        flux[m] += spectrumRadiusTerm * deltaQ[m];
        flux[m] *= 0.5 * faceArea;
    }

    flux[nLaminar] = 0.0;
    if (nElectronIndex >= 0)
    {
        flux[nLaminar - 1] = 0.0;
    }

    // Modify the fluxes of vibration and electron energy terms.
    for (int m = 0; m < nTemperatureModel - 1; ++ m)
    {
        flux[nLaminar + 1 + m] += spectrumRadiusTerm * deltaQ[nLaminar + 1 + m];
        flux[nLaminar + 1 + m] *= 0.5 * faceArea;
    }
    delete []speciesBeta;
}

//! The fuction added by LiPeng computes the formula (Mv +- Rv) * dQ, Mv and Rv denote the Jacobian matrix and spectrum radius of viscous flux,respectively,\n
//! dQ is the increment of conservative vector.
void MVXDQ_STD_LP(RDouble *prim, RDouble nxs, RDouble nys, RDouble nzs, RDouble area, RDouble vol, RDouble visl, RDouble vist, RDouble gama, 
                  RDouble *dq, RDouble *f, int nm, int nLaminar, RDouble radius, int ipn)
{
    using namespace IDX;
    using namespace GAS_SPACE;

    RDouble rm = prim[IR];//!the density
    RDouble um = prim[IU];//!the velocity u in x-direction
    RDouble vm = prim[IV];//!the velocity v in y-direction
    RDouble wm = prim[IW];//!the velocity w in z-direction
    RDouble pm = prim[IP];//!the pressure
    RDouble UU = um * um + vm * vm + wm * wm;
    RDouble Un = um * nxs + vm * nys + wm * nzs;
    RDouble Uq = um * dq[1] + vm * dq[2] + wm * dq[3];

    RDouble prl = GlobalDataBase::GetDoubleParaFromDB("prl");
    RDouble prt = GlobalDataBase::GetDoubleParaFromDB("prt");

    RDouble scl = GlobalDataBase::GetDoubleParaFromDB("sc_l");
    RDouble sct = GlobalDataBase::GetDoubleParaFromDB("sc_t");

    RDouble refReNumber = GlobalDataBase::GetDoubleParaFromDB("refReNumber"); //! Obtain the reynold number.
    RDouble coef = area * area / (vol * refReNumber);
    RDouble bb = (visl + vist)/rm;
    RDouble b1 = 1.0/3.0 * bb * (-Un * dq[0] + nxs * dq[1] + nys * dq[2] + nzs * dq[3]);
    RDouble phi = gama * (visl/prl + vist/prt) / rm;

    RDouble cp = 0.0, Tm = 0.0; //!the constant pressure specific heat and temperature.
    gas->GetTemperatureAndCP(prim, Tm, cp);

    int nChemical = GlobalDataBase::GetIntParaFromDB("nchem");
    int ns = nLaminar + nChemical - nm;//!the number of species
    RDouble Psi = -cp * (visl/prl + vist/prt) * Tm/rm;
    RDouble kn = 0.0;//!the kinetic energy of turbulence
    Psi += -2.0*bb * UU + phi * 0.5 * UU - bb * kn;

    RDouble *Ds = new RDouble[ns];
    RDouble *Phi_s = new RDouble[ns];
    RDouble *theta_s = new RDouble[ns];
    RDouble fsum = 0.0;

    if (nChemical == 1) //!the relevant terms of chemical reaction
    {
        RDouble K = visl/prl + vist/prt;
        RDouble rD = visl/scl + vist/sct;

        RDouble *hs = new RDouble[ns];
        gas->GetEverySpeciesEnthalpy(Tm, Tm, Tm, hs);//!Obtain the enthalpy of each species.

        //!Compute the coefficient theta_s, Phi_s and Ds.
        gas->ComputeCoefficientInMVXDQ(prim, K, rD, theta_s, Phi_s, Ds);

        Psi += Ds[ns - 1] * prim[nm + ns - 1] * hs[ns - 1];
        for (int i = 0; i < ns - 1; i++)
        {
            fsum += (Phi_s[i]/rm - phi * theta_s[i]) * dq[nm + i];
            Psi += Ds[i] * prim[nm + i] * hs[i];
            Psi += phi * prim[nm + i] * theta_s[i];
            Psi += -prim[nm + i] * Phi_s[i]/rm;
        }
        delete []hs;
    }

    //!To compute M*dQ
    //!f[IR ] = 0.0 * dq[IR ] + b1 * 0.0 - dq[0] * 0.0;
    f[IR ] = 0.0; //!Simplification
    f[IRU] = bb  * dq[IRU] + b1 * nxs - dq[0] * bb * um;
    f[IRV] = bb  * dq[IRV] + b1 * nys - dq[0] * bb * vm;
    f[IRW] = bb  * dq[IRW] + b1 * nzs - dq[0] * bb * wm;
    f[IRE] = phi * dq[IRE] + Uq * (2.0 * bb - phi) + b1 * (2.0*Un) + dq[0] * Psi + fsum; //!fsum denotes the computing result of species terms.
    for (int iLaminar = nm; iLaminar < nLaminar; ++ iLaminar)
    {
        //!f[i] = Ds[i - nm] * dq[i] + b1 * 0.0 + dq[0] * 0.0;
        f[iLaminar] = Ds[iLaminar - nm] * dq[iLaminar]; //!Simplification
    }

    for (int iLaminar = 0; iLaminar < nLaminar; ++ iLaminar)
    {
        f[iLaminar] *= 0.5 * coef;
    }

    //!To add rM * dQ that is the product of the viscous spectrum radius rM and increment dQ.
    RDouble ipnrad = 0.5 * ipn * radius;
    for (int iLaminar = 0; iLaminar < nLaminar; ++ iLaminar)
    {
        f[iLaminar] += ipnrad * dq[iLaminar];
    }

    //!Remove the dynamic memory space.
    delete []Ds;
    delete []Phi_s;
    delete []theta_s;
}

void PreconMatrixFreeDF(RDouble *q, RDouble *dq, RDouble *nxyz, RDouble **MG, RDouble gama, RDouble *df, int nm)
{
    using namespace IDX;

    RDouble gamm1 = gama - 1.0;
    RDouble nxs, nys, nzs, vgn, ns;
    RDouble refMachNumber = GlobalDataBase::GetDoubleParaFromDB("refMachNumber");

    nxs = nxyz[0];
    nys = nxyz[1];
    nzs = nxyz[2];
    ns  = nxyz[3];
    vgn = nxyz[4];

    RDouble rm,um,vm,wm,pm,em,cm,vn,v2,eig;
    RDouble mx,my,mz;

    RDouble mach2 = refMachNumber * refMachNumber;
    rm = q[IR];
    um = q[IU];
    vm = q[IV];
    wm = q[IW];
    pm = q[IP];
    cm = sqrt( gama * pm / rm );

    mx = rm * um;
    my = rm * vm;
    mz = rm * wm;
    em = pm / gamm1 + half * rm * ( um * um + vm * vm + wm * wm );

    //! Normal velocity and Eigenvalues.
    vn = nxs * um + nys * vm + nzs * wm - vgn;
    v2 = um * um + vm * vm + wm * wm;
    //eig= ABS( vn ) + sqrt( gama * pm / rm );

    RDouble c2  = cm * cm;
    //RDouble beta = MAX( MIN( v2 / c2,1.0 ),3.0 * mach2 );
    RDouble beta = MIN( MAX( v2 / c2,3.0 * mach2 ),1.0 );
    RDouble vnew  = half * vn * ( 1 + beta );
    RDouble cnew  = half * sqrt( ( ( beta - 1 ) * vn ) * ( ( beta - 1 ) * vn ) + 4 * beta * c2 );
    eig   = ( ABS( vnew ) +  cnew  );

    //! Need to find out the fluxes on level n.
    df[IR ] = - ( rm * vn );
    df[IRU] = - ( rm * um * vn + nxs * pm );
    df[IRV] = - ( rm * vm * vn + nys * pm );
    df[IRW] = - ( rm * wm * vn + nzs * pm );
    df[IRE] = - ( ( em + pm ) * vn + vgn * pm );

    //! Conservative variable on level n+1.
    rm += dq[IR ];
    mx += dq[IRU];
    my += dq[IRV];
    mz += dq[IRW];
    em += dq[IRE];

    um  = mx / rm;
    vm  = my / rm;
    wm  = mz / rm;

    vn  = nxs * um + nys * vm + nzs * wm - vgn;

    pm = gamm1 * ( em - half * rm * ( um * um + vm * vm + wm * wm ) );
    rm = ABS(rm);
    pm = ABS(pm);

    //! Now the flux difference due to dq.
    df[IR ] += ( rm * vn );
    df[IRU] += ( rm * um * vn + nxs * pm );
    df[IRV] += ( rm * vm * vn + nys * pm );
    df[IRW] += ( rm * wm * vn + nzs * pm );
    df[IRE] += ( ( em + pm ) * vn + vgn * pm );

    PreconMatrixDF(df, MG, nm);

    //! Subtract eigenvalue terms from the flux difference.

    for ( int m = 0; m < nm; ++ m )
    {
        df[m] -= eig * dq[m];
    }

    //! Divide by 2 out of necessity.
    for ( int m = 0; m < nm; ++ m )
    {
        df[m] *= half * ns;
    }
}

void PreconMatrixFreeDF(RDouble *q, RDouble *dq, RDouble *gykb,  RDouble **MG, RDouble tm, RDouble gama, RDouble *df, int nm, int ipn)
{
    using namespace IDX;

    RDouble gamm1 = gama - 1.0;
    RDouble refMachNumber = GlobalDataBase::GetDoubleParaFromDB("refMachNumber");

    RDouble rm,um,vm,wm,pm,em,cm,vn,v2,eig;
    RDouble mx,my,mz;

    RDouble mach2 = refMachNumber * refMachNumber;
    rm = q[IR];
    um = q[IU];
    vm = q[IV];
    wm = q[IW];
    pm = q[IP];
    //tm = gama * mach2 * pm / rm;
    cm = sqrt(gama * pm / rm);

    mx = rm * um;
    my = rm * vm;
    mz = rm * wm;
    em = pm / gamm1 + half * rm * (um * um + vm * vm + wm * wm);

    //! Normal velocity and Eigenvalues.
    vn = gykb[0] * um + gykb[1] * vm + gykb[2] * wm - gykb[4];
    v2 = um * um + vm * vm + wm * wm;

    RDouble c2  = cm * cm;
    RDouble beta = MIN(MAX(v2 / c2,3.0 * mach2),1.0);
    //RDouble beta = MIN(MAX(v2 / c2, mach2),1.0);
    RDouble vnew  = half * vn * (1 + beta);
    RDouble cnew  = half * sqrt(((beta - 1) * vn) * ((beta - 1) * vn) + 4 * beta * c2);
    eig   = (ABS(vnew) +  cnew);

    //! Need to find out the fluxes on leverl n.
    df[IR] = - (rm * vn);
    df[IRU] = - (rm * um * vn + gykb[0] * pm);
    df[IRV] = - (rm * vm * vn + gykb[1] * pm);
    df[IRW] = - (rm * wm * vn + gykb[2] * pm);
    df[IRE] = - ((em + pm) * vn + gykb[4] * pm);

    //! Conservative variable on level n+1.
    rm += dq[IR];
    mx += dq[IRU];
    my += dq[IRV];
    mz += dq[IRW];
    em += dq[IRE];


    um  = mx / rm;
    vm  = my / rm;
    wm  = mz / rm;

    vn  = gykb[0] * um + gykb[1] * vm + gykb[2] * wm - gykb[4];

    pm = gamm1 * (em - half * rm * (um * um + vm * vm + wm * wm));
    rm = ABS(rm);
    pm = ABS(pm);

    //! Now the flux difference due to dq.
    df[IR] += (rm * vn);
    df[IRU] += (rm * um * vn + gykb[0] * pm);
    df[IRV] += (rm * vm * vn + gykb[1] * pm);
    df[IRW] += (rm * wm * vn + gykb[2] * pm);
    df[IRE] += ((em + pm) * vn + gykb[4] * pm);

    ComputePreconMatrix(1, MG, q, tm, gama);
    PreconMatrixDF(df, MG, nm);

    //! Subtract eigenvalue terms from the flux difference.

    for (int m = 0; m < nm; ++ m)
    {
        df[m] += ipn * eig * dq[m];
    }

    //! Divide by 2 out of necessity.
    for (int m = 0; m < nm; ++ m)
    {
        df[m] *= half * gykb[3];
    }
}

void PreconMatrixFreeDF(RDouble *q, RDouble *dq, RDouble nxs, RDouble nys, RDouble nzs, RDouble area, RDouble vgn, RDouble **MG, RDouble gama, RDouble *df, int nm, int ipn)
{
    using namespace IDX;

    RDouble gamm1 = gama - 1.0;
    RDouble refMachNumber = GlobalDataBase::GetDoubleParaFromDB("refMachNumber");

    RDouble rm,um,vm,wm,pm,em,cm,vn,v2,eig;
    RDouble mx,my,mz;

    RDouble mach2 = refMachNumber * refMachNumber;
    rm = q[IR];
    um = q[IU];
    vm = q[IV];
    wm = q[IW];
    pm = q[IP];
    cm = sqrt(gama * pm / rm);

    mx = rm * um;
    my = rm * vm;
    mz = rm * wm;
    em = pm / gamm1 + half * rm * (um * um + vm * vm + wm * wm);

    //! Normal velocity and Eigenvalues.
    vn = nxs * um + nys * vm + nzs * wm - vgn;
    v2 = um * um + vm * vm + wm * wm;

    RDouble c2  = cm * cm;
    RDouble beta = MIN(MAX(v2 / c2,3.0 * mach2),1.0);
    RDouble vnew  = half * vn * (1 + beta);
    RDouble cnew  = half * sqrt(((beta - 1) * vn) * ((beta - 1) * vn) + 4 * beta * c2);
    eig   = (ABS(vnew) +  cnew);

    //! Need to find out the fluxes on level n.
    df[IR] = - (rm * vn);
    df[IRU] = - (rm * um * vn + nxs * pm);
    df[IRV] = - (rm * vm * vn + nys * pm);
    df[IRW] = - (rm * wm * vn + nzs * pm);
    df[IRE] = - ((em + pm) * vn + vgn * pm);

    //! Conservative variable on level n+1.
    rm += dq[IR];
    mx += dq[IRU];
    my += dq[IRV];
    mz += dq[IRW];
    em += dq[IRE];


    um  = mx / rm;
    vm  = my / rm;
    wm  = mz / rm;

    vn  = nxs * um + nys * vm + nzs * wm - vgn;

    pm = gamm1 * (em - half * rm * (um * um + vm * vm + wm * wm));
    rm = ABS(rm);
    pm = ABS(pm);

    //! Now the flux difference due to dq.
    df[IR] += (rm * vn);
    df[IRU] += (rm * um * vn + nxs * pm);
    df[IRV] += (rm * vm * vn + nys * pm);
    df[IRW] += (rm * wm * vn + nzs * pm);
    df[IRE] += ((em + pm) * vn + vgn * pm);

    PreconMatrixDF(df, MG, nm);

    //! Subtract eigenvalue terms from the flux difference.

    for (int m = 0; m < nm; ++ m)
    {
        df[m] += ipn * eig * dq[m];
    }

    //! Divide by 2 out of necessity.
    for (int m = 0; m < nm; ++ m)
    {
        df[m] *= half * area;
    }
}

}
