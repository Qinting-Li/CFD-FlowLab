#include "Boundary.h"
#include "Constants.h"
#include "Math_BasisFunction.h"
#include "Gas.h"

using namespace std;

namespace PHSPACE
{
//! Subsonic inlet boundary: static pressure extrapolation method.
//! extrapolation static pressure --> local mach --> local static T --> density --> Sound velocity --> velocity.
void Pressure_InletBC (RDouble *prims, RDouble *primt, RDouble gama, RDouble mach_inf, RDouble totalP, RDouble totalT, RDouble* angle, int nl, int nchem)
{
	using namespace IDX;
	int neqn = nl + nchem;
	for ( int m = 0; m < neqn; ++ m )
	{
		primt[m] = prims[m];
	}

	 RDouble rm, um, vm, wm, V, pm, mach, tm, am;
	 RDouble alpha = angle[0], beita = angle[1], theta = angle[2];
	 RDouble gama1 = gama - 1;

	 pm = MIN(prims[IP], static_cast <RDouble> (totalP * 0.99));
	 mach = sqrt( (pow(totalP / pm, gama1 / gama ) - 1) * 2.0 / gama1);
	 tm = MIN(totalT / (1 + gama1 * mach * mach / 2.0), totalT * 0.99);
	 rm = pm * gama * mach_inf * mach_inf / tm;
	 am = sqrt(gama * pm / rm);
	 V  = mach * am;

	 um = V * alpha;
	 vm = V * beita;
	 wm = V * theta;

	 primt[IR] = rm;
	 primt[IU] = um;
	 primt[IV] = vm;
	 primt[IW] = wm;
	 primt[IP] = pm;
}

//! Subsonic inlet boundary: Riemann's no reflecting boundary method.
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

	RDouble Tb = MIN(totalT * (cb2 / c02), static_cast <RDouble> (totalT * 0.999));
	RDouble Pb = MIN(totalP * pow(Tb / totalT, gama / gama1), static_cast <RDouble> (totalP * 0.999));
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

void InflowBC(RDouble *prims, RDouble *primt, int nl, int nchem)
{
	int neqn = nl + nchem;
    for (int m = 0; m < neqn; ++ m)
    {
        primt[m] = prims[m];
    }
}

void Isentropic_Inflow_BC(RDouble *prims, RDouble *primt, RDouble gama, int nl, int nchem)
{
	using namespace IDX;
	int neqn = nl + nchem;
    RDouble rgo, pgo, rgn, pgn;
    rgo = primt[IR];
    pgo = primt[IP];

    pgn = prims[IP];
    rgn = rgo * pow( static_cast<RDouble>(pgn / pgo), static_cast<RDouble>(1.0 / gama) );

    primt[IR] = rgn;
    primt[IP] = pgn;
}

void OutflowBC(RDouble *prims, RDouble *primt, int nl, int nchem)
{
	int neqn = nl + nchem;
    for (int m = 0; m < neqn; ++ m)
    {
        primt[m] = prims[m];
    }
}

void MassOutflow_BC(RDouble *prims, RDouble *primt, RDouble mass_alpha, int nl, int nchem)
{
	using namespace IDX;
	int neqn = nl + nchem;
    for ( int m = 0; m < neqn; ++ m )
    {
        primt[m] = prims[m];
    }

    primt[IU] *= mass_alpha;
    primt[IV] *= mass_alpha;
    primt[IW] *= mass_alpha;
}

void SymmetryBC(RDouble *prims, RDouble *primt, RDouble nxs, RDouble nys, RDouble nzs, RDouble vgn, int nl, int nchem)
{
	using namespace IDX;
    int neqn = nl + nchem;
    for (int m = 0; m < neqn; ++ m)
    {
        primt[m] = prims[m];
    }

    RDouble um, vm, wm, vn;
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
	using namespace GAS_SPACE;
	RDouble coefficientofstateEquation = gas->GetCoefficientOfStateEquation();
    if (nchem == 1)
	{
        omav = gas->ComputeMolecularWeightReciprocal(prims);
	}
	
    RDouble tm = (1.0 / coefficientofstateEquation) * prims[IP] / prims[IR];

	RDouble rg, pg, tg;     //! g means ghost cell.

    pg = prims[IP];
    tg = two * tw - tm;
    if (tg <= 0.0) tg = tw;

	rg = pg / (coefficientofstateEquation * tg * omav);
    
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

    //! Inner point variables.
    rin = prims[IR];
    uin = prims[IU];
    vin = prims[IV];
    win = prims[IW];
    pin = prims[IP];

    //! Infinite variables.
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

    //! Supersonic.
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

    //! Subsonic.
    riemp = vni + 2.0 * cin / gamm1;
    riemm = vno - 2.0 * coo / gamm1;
    vnb =  half * ( riemp + riemm );
    cb  =  0.25 * ( riemp - riemm ) * gamm1;
    if (vnb >= 0.0)
    {
        entr = pin / pow(rin, gama);
        vtx  = uin - nxs * vni;
        vty  = vin - nys * vni;
        vtz  = win - nzs * vni;

 /*       riemp = ABS( vni ) + 2.0 * cin / gamm1;
        riemm = ABS( vno ) - 2.0 * coo / gamm1;*/

    }
    else
    {
        //! Inlet.
        entr = poo / pow(roo, gama);
        vtx  = uoo - nxs * vno;
        vty  = voo - nys * vno;
        vtz  = woo - nzs * vno;
//cb  =   0.25 * ( riemp - riemm ) * gamm1;
    }

    rb  = pow((cb * cb / (entr * gama)), one / gamm1);
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

void PreconFarfieldBC(RDouble *prims, RDouble *primInflow, RDouble *primt, RDouble gama0, RDouble gama, RDouble nxs, RDouble nys, RDouble nzs, RDouble vgn, RDouble refMachNumber, int nl, int nchem)
{
	using namespace IDX;
	int neqn = nl + nchem;

	RDouble rin,uin,vin,win,pin,tin,vni;
	RDouble roo,uoo,voo,woo,poo,too,vno;
	RDouble rb,ub,vb,wb,pb,tb;
	RDouble mach2 = refMachNumber * refMachNumber;

	//! Inner point variables.
	rin = prims[IR];
	uin = prims[IU];
	vin = prims[IV];
	win = prims[IW];
	pin = prims[IP];
	tin = gama * mach2 * pin / rin;

	//! Infinite variables.
	roo = primInflow[IR];
	uoo = primInflow[IU];
	voo = primInflow[IV];
	woo = primInflow[IW];
	poo = primInflow[IP];
	too = gama * mach2 * poo / roo;

	vno = nxs * uoo + nys * voo + nzs * woo - vgn;
	vni = nxs * uin + nys * vin + nzs * win - vgn;

	//! Subsonic.
	if (refMachNumber < 1.0)
	{
		if (vni < -1.e-10)
		{			   
		    ub  = uoo;
			vb  = voo;
			wb  = woo;
			pb  = pin;
			tb  = too;
			rb  = gama * mach2 * pb / tb;
		}
		else
		{
			ub  = uin;
			vb  = vin;
			wb  = win;
			pb  = poo;
			tb  = tin;
			rb  = gama * mach2 * pb / tb;
		}

	}
	else
	{
		if ( vni >= 0.0 )
		{
			ub  = uin;
			vb  = vin;
			wb  = win;
			pb  = pin;
			rb  = rin;
		}
		else
		{
			ub  = uoo;
			vb  = voo;
			wb  = woo;
			pb  = poo;
			rb  = roo;
		}
	}
	
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

    //! Inner point variables.
    rin = prims[ IR ];
    uin = prims[ IU ];
    vin = prims[ IV ];
    win = prims[ IW ];
    pin = prims[ IP ];

    //! Infinite variables.
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
        //! Supersonic.
        if ( vnin >= 0.0 )
        {
            //! Supersonic outlet.
            for ( int m = 0; m < neqn; ++ m )
            {
                primt[ m ] = prims[ m ];
            }
        }
        else
        {
            //! Supersonic inlet.
            for ( int m = 0; m < neqn; ++ m )
            {
                primt[ m ] = primInflow[ m ];
            }
        }
    }
    else
    {
        //! Subsonic.
        RDouble gama1 = gama - one;
        RDouble riemp;    //! Riemann invariant for inner.
        RDouble riemm;    //! Riemann invariant for infinite.
        riemp = vnin + 2.0 * cinner / gama1;
        riemm = vnoo - 2.0 * coo    / gama1;

        vnb   = half   * ( riemp + riemm );
        cb    = fourth * ( riemp - riemm ) * gama1;

        RDouble sb;                         //! Entropy at the boundary.
        RDouble uref, vref, wref, vnref;    //! Outlet refers to inner point, and inlet refers to infinite.
        
        if ( vnb > 0.0 )
        {
            //! Subsonic outlet.
            sb    = pin / pow( rin, gama );
            uref  = uin;
            vref  = vin;
            wref  = win;
            vnref = vnin;
        }
        else
        {
            //! Subsonic inlet.
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
		pb  = cb * cb * rb / gama0;    //! The same to "pb = sb * pow( rb, gama0 );".

        primt[IR] = rb;
        primt[IU] = ub;
        primt[IV] = vb;
        primt[IW] = wb;
        primt[IP] = pb;
    }
}

}
