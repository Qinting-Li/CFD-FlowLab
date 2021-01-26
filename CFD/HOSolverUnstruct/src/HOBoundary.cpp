//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//          PPPPP  H   H  EEEEE  N    N  GGGGG  L      EEEEE  III         +  
//          P   P  H   H  E      NN   N  G      L      E       I          +
//          PPPPP  HHHHH  EEEEE  N N  N  G  GG  L      EEEEE   I          +
//          P      H   H  E      N  N N  G   G  L      E       I          +
//          P      H   H  EEEEE  N    N  GGGGG  LLLLL  EEEEE  III         +
//------------------------------------------------------------------------+
//          Platform for Hybrid Engineering Simulation of Flows           +
//          China Aerodynamics Research and Development Center            +
//                     (C) Copyright, Since 2010                          +
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file	HOBoundary.cpp
//! @brief	basis function, such as monomial, orthogonalized etc.\n
//!
//! @author		gongxq, liming.
//! @email       
//! @date		2019-11-09
//! @version	A001
//! @copyright

#include "HODefine.h"
#include "HOBoundary.h"
#include "Constants.h"
#include "Math_BasisFunction.h"

using PHSPACE::RDouble;
//using namespace PHSPACE;
namespace HOUnstruct
{
void SymmetryBoundary( int iGausspoint, int numberGaussPoint, RDouble xNormal, RDouble yNormal, RDouble zNormal, vector < vector < RDouble > > & conservationQ )
{
    RDouble rhom,um, vm, wm, pm, rhoem, vn;
    RDouble rhoR,uR, vR, wR;

    rhom  = conservationQ[0][iGausspoint];
    um    = conservationQ[1][iGausspoint]/rhom;
    vm    = conservationQ[2][iGausspoint]/rhom;
    wm    = conservationQ[3][iGausspoint]/rhom;
    rhoem = conservationQ[4][iGausspoint];
    pm    = HO_GAMAM1 * (rhoem - 0.5 * rhom * (um * um + vm * vm + wm * wm));

    vn = xNormal * um + yNormal * vm + zNormal * wm;

    int iGaussPointR = iGausspoint + numberGaussPoint;

    rhoR  = rhom;
    uR    = um - two * xNormal * vn;
    vR    = vm - two * yNormal * vn;
    wR    = wm - two * zNormal * vn;

    conservationQ[0][iGaussPointR] = rhoR;
    conservationQ[1][iGaussPointR] = rhoR*uR;
    conservationQ[2][iGaussPointR] = rhoR*vR;
    conservationQ[3][iGaussPointR] = rhoR*wR;
    conservationQ[4][iGaussPointR] = pm/HO_GAMAM1 + 0.5 * rhoR * (uR * uR + vR * vR + wR * wR);
}

void SymmetryBC(RDouble *prims, RDouble *primt, RDouble nxs, RDouble nys, RDouble nzs, RDouble vgn, int nl, int nchem)
{
    using namespace IDX;
    int neqn = nl + nchem;
    for (int m = 0; m < neqn; ++ m)
    {
        primt[m] = prims[m];
    }

    double um, vm, wm, vn;
    um = prims[IU];
    vm = prims[IV];
    wm = prims[IW];

    vn = nxs * um + nys * vm + nzs * wm - vgn;

    primt[IU] = um - two * nxs * vn;
    primt[IV] = vm - two * nys * vn;
    primt[IW] = wm - two * nzs * vn;
}

//! 3D.
//! Viscous adiabatic wall.
void ViscousAdiabaticWallBC(RDouble *prims, RDouble *primt, RDouble nxs, RDouble nys, RDouble nzs, RDouble xtn, RDouble ytn, RDouble ztn, RDouble refMachNumber, int nl, int nchem)
{
    using namespace IDX;
    int neqn = nl + nchem;
    for (int m = 0; m < neqn; ++ m)
    {
        primt[m] = prims[m];
    }

    primt[IU] = - prims[IU] + two * xtn;
    primt[IV] = - prims[IV] + two * ytn;
    primt[IW] = - prims[IW] + two * ztn;
}

void ViscousIsotropicWallBC(RDouble *prims, RDouble *primt, RDouble nxs, RDouble nys, RDouble nzs, RDouble xtn, RDouble ytn, RDouble ztn, RDouble refMachNumber, RDouble tw, int nl, int nchem)
{
    using namespace IDX;
    int neqn = nl + nchem;
    for (int m = 0; m < neqn; ++ m)
    {
        primt[m] = prims[m];
    }

    primt[IU] = - prims[IU] + two * xtn;
    primt[IV] = - prims[IV] + two * ytn;
    primt[IW] = - prims[IW] + two * ztn;

    RDouble omav = one;
    RDouble coefficientOfStateEquation = 1.0 / ( HO_GAMA * refMachNumber * refMachNumber );

    RDouble tm = (1.0 / coefficientOfStateEquation) * prims[IP] / prims[IR];

    RDouble rg, pg, tg; //g 代表ghost

    pg = prims[IP];
    tg = two * tw - tm;
    if (tg <= 0.0) tg = tw;

    rg = pg / (coefficientOfStateEquation * tg * omav);

    primt[IR] = rg;
    primt[IP] = pg;
}

void FarfieldBC(RDouble *prims, RDouble *primInflow, RDouble *primt, RDouble gama0, RDouble gama, RDouble nxs, RDouble nys, RDouble nzs, RDouble vgn, int nl, int nchem)
{
    using namespace IDX;
    int neqn = nl + nchem;

    RDouble rin,uin,vin,win,pin,cin,vni;
    RDouble roo,uoo,voo,woo,poo,coo,vno;
    RDouble lam1,lam2,lam3,gamm1,vei;
    RDouble vtx ,vty ,vtz ,entr ;
    RDouble riemm,riemp;
    RDouble rb,ub,vb,wb,pb,vnb,cb;

    //!内点
    rin = prims[IR];
    uin = prims[IU];
    vin = prims[IV];
    win = prims[IW];
    pin = prims[IP];

    //!远场
    roo = primInflow[IR];
    uoo = primInflow[IU];
    voo = primInflow[IV];
    woo = primInflow[IW];
    poo = primInflow[IP];

    vno = nxs * uoo + nys * voo + nzs * woo - vgn;
    vni = nxs * uin + nys * vin + nzs * win - vgn;

    coo = sqrt(ABS(gama0 * poo / roo));
    cin = sqrt(ABS(gama  * pin / rin));

    gamm1 = gama - one;

    lam1  = vni;
    lam2  = vni - cin;
    lam3  = vni + cin;

    vei = sqrt(uin * uin + vin * vin + win * win);

    //超声速
    if (vei > cin)
    {
        if (vni >= 0.0)
        {
            for (int m = 0; m < neqn; ++ m)
            {
                primt[m] = prims[m];
            }
        }
        else
        {
            for (int m = 0; m < neqn; ++ m)
            {
                primt[m] = primInflow[m];
            }
        }
        return;
    }

    //若程序执行到这里，表明一定是亚声速
    riemp = vni + 2.0 * cin / gamm1;
    riemm = vno - 2.0 * coo / gamm1;
    vnb =  half * ( riemp + riemm );
    cb  =  0.25 * ( riemp - riemm ) * gamm1;
    if (vnb >= 0.0)
    {
        // exit
        entr = pin / pow(rin, gama);
        vtx  = uin - nxs * vni;
        vty  = vin - nys * vni;
        vtz  = win - nzs * vni;

        /*       riemp = ABS( vni ) + 2.0 * cin / gamm1;
        riemm = ABS( vno ) - 2.0 * coo / gamm1;*/

    }
    else
    {
        //inlet
        entr = poo / pow(roo, gama);
        vtx  = uoo - nxs * vno;
        vty  = voo - nys * vno;
        vtz  = woo - nzs * vno;
        //cb  =   0.25 * ( riemp - riemm ) * gamm1;
    }

    rb  = pow((cb * cb / (entr * gama)), static_cast<RDouble>(one) / gamm1);
    ub  = vtx + nxs * vnb;
    vb  = vty + nys * vnb;
    wb  = vtz + nzs * vnb;
    pb  = cb * cb * rb / gama;

    primt[IR] = rb;
    primt[IU] = ub;
    primt[IV] = vb;
    primt[IW] = wb;
    primt[IP] = pb;
}

void FarFieldRiemannInvariants(RDouble *prims, RDouble *primInflow, RDouble *primt, RDouble gama0, RDouble gama, RDouble nxs, RDouble nys, RDouble nzs, RDouble vgn, int nl, int nchem)
{
	using namespace IDX;
    int neqn = nl + nchem;

	RDouble rin, uin, vin, win, pin, vnin, cinner, vein;
	RDouble roo, uoo, voo, woo, poo, vnoo, coo;
	RDouble rb , ub , vb , wb , pb , vnb , cb;

    //! inner point variables.
    rin = prims[ IR ];
    uin = prims[ IU ];
    vin = prims[ IV ];
    win = prims[ IW ];
    pin = prims[ IP ];

    //! infinite variables.
    roo = primInflow[ IR ];
    uoo = primInflow[ IU ];
    voo = primInflow[ IV ];
    woo = primInflow[ IW ];
    poo = primInflow[ IP ];

    vnin = nxs * uin + nys * vin + nzs * win - vgn;
    vnoo = nxs * uoo + nys * voo + nzs * woo - vgn;

    cinner = sqrt( ABS( gama  * pin / rin ) );
    coo    = sqrt( ABS( gama0 * poo / roo ) );

    vein   = sqrt( uin * uin + vin * vin + win * win );

    if ( vein > cinner )
    {
        //! supersonic.
        if ( vnin >= 0.0 )
        {
            //! supersonic outlet.
            for ( int m = 0; m < neqn; ++ m )
            {
                primt[ m ] = prims[ m ];
            }
        }
        else
        {
            //! supersonic inlet.
            for ( int m = 0; m < neqn; ++ m )
            {
                primt[ m ] = primInflow[ m ];
            }
        }
    }
    else
    {
        //! subsonic.
        RDouble gama1 = gama - one;
        RDouble riemp; //! Riemann invariant for inner.
        RDouble riemm; //! Riemann invariant for infinite.
        riemp = vnin + 2.0 * cinner / gama1;
        riemm = vnoo - 2.0 * coo    / gama1;

        vnb   = half   * ( riemp + riemm );
        cb    = fourth * ( riemp - riemm ) * gama1;

        double sb;                      //! entropy at the boundary.
        double uref, vref, wref, vnref; //! Outlet refers to inner point, and inlet refers to infinite.
        
        if ( vnb > 0.0 )
        {
            //! subsonic outlet.
            sb    = pin / pow( rin, gama );
            uref  = uin;
            vref  = vin;
            wref  = win;
            vnref = vnin;
        }
        else
        {
            //! subsonic inlet.
            sb    = poo / pow( roo, gama0 );
            uref  = uoo;
            vref  = voo;
            wref  = woo;
            vnref = vnoo;
        }

		rb  = pow( ( cb * cb / ( sb * gama ) ), 1.0 / gama1 );
        ub  = uref + nxs * ( vnb - vnref );
        vb  = vref + nys * ( vnb - vnref );
        wb  = wref + nzs * ( vnb - vnref );
		pb  = cb * cb * rb / gama0; //! The same to "pb = sb * pow( rb, gama0 );"

        primt[IR] = rb;
        primt[IU] = ub;
        primt[IV] = vb;
        primt[IW] = wb;
        primt[IP] = pb;
    }
}

void InflowBC(RDouble *prims, RDouble *primt, int nl, int nchem)
{
    int neqn = nl + nchem;
    for (int m = 0; m < neqn; ++ m)
    {
        primt[m] = prims[m];
    }
}

void OutflowBC(RDouble *prims, RDouble *primt, int nl, int nchem)
{
    int neqn = nl + nchem;
    for (int m = 0; m < neqn; ++ m)
    {
        primt[m] = prims[m];
    }
}

void PressureOutletBC(RDouble *prims, RDouble *primt, RDouble pres, int nl, int nchem)
{
    using namespace IDX;
    int neqn = nl + nchem;
    for (int m = 0; m < neqn; ++ m)
    {
        primt[m] = prims[m];
    }

    primt[IP] = pres;
}

void PressureOutletBC(RDouble *prims, RDouble *primt, int nl, int nchem)
{
    using namespace IDX;
    int neqn = nl + nchem;
    for ( int m = 0; m < neqn; ++ m )
    {
        primt[m] = prims[m];
    }

    primt[IP] = prims[IP];
    primt[IU] = fabs(prims[IU]);
}

void PressureInletBCRiemann(RDouble *prims, RDouble *primt, RDouble gama,  RDouble mach_inf,
                            RDouble nx, RDouble ny, RDouble nz, RDouble totalP, RDouble totalT, RDouble* angle, int nl, int nchem)
{
    using namespace IDX;
    int neqn = nl + nchem;
    for (int m = 0; m < neqn; ++ m)
    {
        primt[m] = prims[m];
    }

    RDouble alpha = angle[0], beita = angle[1], theta = angle[2];
    RDouble gama1 = gama - 1;
    RDouble p = prims[IP], u = prims[IU], v = prims[IV], w = prims[IW], r = prims[IR];

    RDouble cd = sqrt(fabs(gama * p / r));
    RDouble vn = u * nx + v * ny + w * nz;
    RDouble v2 = u * u + v * v + w * w;
    RDouble R = vn - 2. * cd / gama1;
    RDouble costheta = - vn / sqrt(v2);

    RDouble c02   = cd * cd + gama1 * v2 / 2.0;
    RDouble temp1 = (gama1 * costheta * costheta + 2) * c02 / (gama1 * R * R);
    RDouble temp2 = gama1 / 2.0;
    RDouble cb    = - R * gama1 * (1 + costheta * sqrt(fabs(temp1 - temp2))) / ( gama1 * costheta * costheta + 2);
    RDouble cb2   = cb * cb;

    RDouble Tb = MIN(totalT * (cb2 / c02), totalT * 0.999);
    RDouble Pb = MIN(totalP * pow(Tb / totalT, gama / gama1), totalP * 0.999);
    RDouble rhob = Pb * gama * mach_inf * mach_inf / Tb;
    RDouble temp3 = totalT / (gama * mach_inf * mach_inf) - Pb / rhob;
    RDouble Vb = sqrt( fabs(2 * (gama / gama1) * temp3));

    primt[IR] = rhob;
    primt[IU] = Vb * alpha;
    primt[IV] = Vb * beita;
    primt[IW] = Vb * theta;
    primt[IP] = Pb;
}

void PressureInletBCPlate (RDouble *prims, RDouble *primt, RDouble gama, RDouble mach_inf, RDouble totalP, RDouble totalT, RDouble* angle, int nl, int nchem)
{
    using namespace IDX;
    int neqn = nl + nchem;
    for (int m = 0; m < neqn; ++ m)
    {
        primt[m] = prims[m];
    }

    RDouble rm, um, vm, wm, V, pm, tm, am;
    RDouble alpha = angle[0], beita = angle[1], theta = angle[2];
    RDouble gama1 = gama - 1;

    pm = totalP;
    tm = totalT;
    rm = pm * gama * mach_inf * mach_inf / tm;

    am = sqrt(gama * pm / rm);
    V  = mach_inf * am;

    um = V * alpha;
    vm = V * beita;
    wm = V * theta;

    primt[IR] = rm;
    primt[IU] = um;
    primt[IV] = vm;
    primt[IW] = wm;
    primt[IP] = pm;
}
}

