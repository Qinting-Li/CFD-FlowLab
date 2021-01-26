#include "PerfectGas.h"
#include "Constants.h"
#include "BC.h"
#include "GlobalDataBase.h"
#include "TK_Exit.h"
#include "PHIO.h"
#include "Glb_Dimension.h"
#include "TK_Log.h"
#include "SubsonicBC.h"

namespace PHSPACE
{
namespace GAS_SPACE
{

PerfectGas::PerfectGas()
{
}

PerfectGas::~PerfectGas()
{
}

void PerfectGas::Primitive2Conservative(RDouble *prim, RDouble gama, RDouble Tv, RDouble Te, RDouble *q)
{
    using namespace IDX;
    RDouble em;

    //!obtain the total internal energy Em.
    ComputeInternalEnergy(prim, gama, Tv, Te, em);

    q[IR ] = prim[IR];
    q[IRU] = prim[IR] * prim[IU];
    q[IRV] = prim[IR] * prim[IV];
    q[IRW] = prim[IR] * prim[IW];
    q[IRE] = prim[IR] * em;
}

void PerfectGas::ComputeInternalEnergy(RDouble *primitiveVariables, RDouble gama, RDouble vibrationTemperature, RDouble electronTemperature, RDouble &totalEnergy)
{
    using namespace IDX;

    RDouble &um = primitiveVariables[IU];
    RDouble &vm = primitiveVariables[IV];
    RDouble &wm = primitiveVariables[IW];
    RDouble &density = primitiveVariables[IR];
    RDouble &pressure = primitiveVariables[IP];
    RDouble squareVelocity = um * um + vm * vm + wm * wm;

    totalEnergy = (pressure / density) / (gama - one) + half * squareVelocity;
}

void PerfectGas::Conservative2Primitive(RDouble *q, RDouble gama, RDouble *prim, RDouble temperature)
{
    using namespace IDX;
    RDouble density, oDensity, um, vm, wm, internalEnergy;
    density = q[IR];
    if (density <= zero)
    {
        prim[IR] = density;
        return;
    }

    oDensity  = 1.0 / density;
    um   = q[IU] * oDensity;
    vm   = q[IV] * oDensity;
    wm   = q[IW] * oDensity;
    internalEnergy  = q[IP] - half * density * (um * um + vm * vm + wm * wm);

    prim[IR] = density;
    prim[IU] = um;
    prim[IV] = vm;
    prim[IW] = wm;
    prim[IP] = (gama - one) * internalEnergy;

}

void PerfectGas::GetSpecificHeatRatio(RDouble *prim, RDouble &gama)
{
    gama = 1.4;
}

void PerfectGas::ComputeEnthalpyByPrimitive(RDouble *primitiveVariables, RDouble gama, RDouble &enthalpy, RDouble* temperatures)
{
    using namespace IDX;
    enthalpy = (gama / (gama - one)) * (primitiveVariables[IP] / primitiveVariables[IR]);
}

void PerfectGas::ComputeTotalEnthalpyAndDH(RDouble *primitiveVariables, RDouble &gama, RDouble *deltaQ, RDouble &totalEnthalpy, RDouble &deltaEnthalpy)
{
    // This function is used for perfect gas and single temperature model.
    // To compute the total enthalpy and the variable dh in the function MXDQ_STD(). dh=b2, b2 denotes the coefficient \n
    // of the vector M*dQ which can be referred to the forumla (A.7) and (A.8) in the appendix A of the PHengLEI Theory manual.
    using namespace IDX;

    // Obtain the primitive variables.
    RDouble &um = primitiveVariables[IU];
    RDouble &vm = primitiveVariables[IV];
    RDouble &wm = primitiveVariables[IW];
    RDouble &density = primitiveVariables[IR];
    RDouble &pressure = primitiveVariables[IP];

    // c1 is product coefficient of the variable dq[IR]£¬gama1 = gama-1¡£
    RDouble squareVelocity, c1, gama1;
    squareVelocity = um * um + vm * vm + wm * wm;
    gama1 = gama - 1.0;
    c1 = half * gama1 * squareVelocity;
    deltaEnthalpy = -gama1 * (um * deltaQ[IRU] + vm * deltaQ[IRV] + wm * deltaQ[IRW] - deltaQ[IRE]);

    // enthalpy is static enthalpy.
    RDouble enthalpy = (gama / gama1)* (pressure / density);

    deltaEnthalpy += c1 * deltaQ[IR];
    // h0 = h+v*v/2, total enthalpy equals to static enthalpy plus half of the squared velocity.
    totalEnthalpy = enthalpy + half * squareVelocity;
}

void PerfectGas::ComputeDHAndTotalEnthalpy(RDouble *primitiveVariables, RDouble &gama, RDouble *deltaQ, RDouble &deltaEnthalpy, RDouble &totalEnthalpy)
{
    using namespace IDX;

    RDouble squareVelocity, ae, af;
    RDouble enthalpy;
    // Obtain the primitive variables.
    RDouble &um = primitiveVariables[IU];
    RDouble &vm = primitiveVariables[IV];
    RDouble &wm = primitiveVariables[IW];
    RDouble &density = primitiveVariables[IR];
    RDouble &pressure = primitiveVariables[IP];

    squareVelocity = um * um + vm * vm + wm * wm;
    ae = gama - one;
    af = half * ae * squareVelocity;
    deltaEnthalpy = -ae * (um * deltaQ[IRU] + vm * deltaQ[IRV] + wm * deltaQ[IRW] - deltaQ[IRE]);
    enthalpy = (gama / (gama - one)) * (pressure / density);
  
    deltaEnthalpy += af * deltaQ[IR];
    // ht = h+v*v/2
    totalEnthalpy = enthalpy + half * squareVelocity;
}

RDouble PerfectGas::GetMixedGasEnthalpy(RDouble *primitiveVars, RDouble transRotationTemperature, RDouble vibrationTemperature, RDouble electronTemperature)
{
    RDouble  staticEnthaly = 3.5 * primitiveVars[IDX::IP] / primitiveVars[IDX::IR];
    return staticEnthaly;
}


void PerfectGas::ComputeViscosityByPrimitive(RDouble *primitiveVariables, RDouble &transRotationTemperature, RDouble &electronTemperature, RDouble &suthTemperature, RDouble &viscosity)
{
    viscosity = transRotationTemperature * sqrt(transRotationTemperature) * (1.0 + suthTemperature) / (transRotationTemperature + suthTemperature);
}

void PerfectGas::InitGlobalParameterOfNSEquation()
{
    //! Obtain the variable of "reynoldsReferenceLengthDimensional" according to "gridScaleFactor".
    InitCommonParameter();

    //! Read parameters from file of "boundary_condition.hypara". Notice that the five reference variables for the component aerodynamic coefficients integration should be updated with "reynoldsReferenceLengthDimensional".
    GlobalBoundaryCondition::ReadGlobalBoundaryCondition();

    GlobalBoundaryCondition::SetGlobalDataBaseByGlobalBC();

    GlobalBoundaryCondition::ChangeBCTypeByGlobalBC();

    InitReferenceParameter();

    ComputeReferenceParameter();

    PrintInflowConditionToWindow();

    ComputeProperReynoldsNumberForGrid();

    //! The five reference variables for the whole aerodynamic coefficients integration, which are from the file of "cfd_para.hypara", should be updated with "reynoldsReferenceLengthDimensional".
    ComputeOtherProperParameterForGrid();

    InitOtherParameter();
}

void PerfectGas::InitNumberOfSpecies()
{

    int numberOfSpecies = 0;
    GlobalDataBase::UpdateData( "numberOfSpecies", &numberOfSpecies, PHINT, 1 );

    int nm    = GlobalDataBase::GetIntParaFromDB( "nm"    );
    int nl = nm;
    GlobalDataBase::UpdateData( "nl", &nl, PHINT, 1 );
}


void PerfectGas::ComputeProperReynoldsNumberForGrid()
{
    RDouble refReNumber = referenceParameterDimensional.GetReynoldsNumber();
    RDouble reynoldsReferenceLengthDimensional = GlobalDataBase::GetDoubleParaFromDB("reynoldsReferenceLengthDimensional");

    refReNumber *= reynoldsReferenceLengthDimensional;
    RDouble oreynolds = 1.0 / refReNumber;

    referenceParameterDimensional.SetReynoldsNumber(refReNumber);
    GlobalDataBase::UpdateData("refReNumber", &refReNumber, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("oreynolds", &oreynolds, PHDOUBLE, 1);
}

void PerfectGas::InitCommonParameter()
{  
    RDouble reynoldsReferenceLengthDimensional = 1.0;
	    
	RDouble gridScaleFactor = GlobalDataBase::GetDoubleParaFromDB("gridScaleFactor");
			
	reynoldsReferenceLengthDimensional = gridScaleFactor;
		
	GlobalDataBase::UpdateData("reynoldsReferenceLengthDimensional", &reynoldsReferenceLengthDimensional, PHDOUBLE, 1);

}

void PerfectGas::InitOtherParameter()
{
    int nolstress, nrokplus;
    nolstress = -1;
    nrokplus  = -1;
    GlobalDataBase::UpdateData( "nolstress", &nolstress, PHINT, 1 );
    GlobalDataBase::UpdateData( "nrokplus" , &nrokplus , PHINT, 1 );

    int outnstep = 0;
    GlobalDataBase::UpdateData( "outnstep", &outnstep, PHINT, 1 );
}

void PerfectGas::InitReferenceParameter()
{
    InitNumberOfSpecies();

    int inflowParaType = GlobalDataBase::GetIntParaFromDB( "inflowParaType" );

    const int NONDIMENSIONCONDITION = 0;
    const int FLIGHTCONDITION       = 1;
    const int EXPERIMENTCONDITION   = 2;
    const int SUBSONICBCCONDITION   = 3;

    if ( inflowParaType == NONDIMENSIONCONDITION )
    {
        RDouble refMachNumber = GlobalDataBase::GetDoubleParaFromDB("refMachNumber");
        RDouble refReNumber   = GlobalDataBase::GetDoubleParaFromDB("refReNumber");
        RDouble refDimensionalTemperature = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");

        referenceParameterDimensional.SetMachNumber(refMachNumber);
        referenceParameterDimensional.SetReynoldsNumber(refReNumber);
        referenceParameterDimensional.SetTemperature(refDimensionalTemperature);
    }
    else if ( inflowParaType == FLIGHTCONDITION )
    {
        RDouble refMachNumber = GlobalDataBase::GetDoubleParaFromDB("refMachNumber");
        RDouble height        = GlobalDataBase::GetDoubleParaFromDB("height");

        referenceParameterDimensional.SetMachNumber(refMachNumber);

        SetAirInformationByDataBase();
    }
    else if ( inflowParaType == EXPERIMENTCONDITION )
    {
        RDouble refMachNumber = GlobalDataBase::GetDoubleParaFromDB("refMachNumber");
        RDouble refDimensionalTemperature = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");
        RDouble refDimensionalPressure = GlobalDataBase::GetDoubleParaFromDB("refDimensionalPressure");

        referenceParameterDimensional.SetMachNumber(refMachNumber);
        referenceParameterDimensional.SetTemperature(refDimensionalTemperature);
        referenceParameterDimensional.SetPressure(refDimensionalPressure);

        SetAirInformationByExpData();
    }
    else if ( inflowParaType == SUBSONICBCCONDITION )
    {
        RDouble refMachNumber = GlobalDataBase::GetDoubleParaFromDB("refMachNumber");
        RDouble refDimensionalTemperature = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");
        RDouble refDimensionalPressure = GlobalDataBase::GetDoubleParaFromDB("refDimensionalPressure");

        referenceParameterDimensional.SetMachNumber(refMachNumber);
        referenceParameterDimensional.SetTemperature(refDimensionalTemperature);
        referenceParameterDimensional.SetPressure(refDimensionalPressure);
    }
    else
    {
        TK_Exit::UnexpectedVarValue( "inflowParaType", inflowParaType );
    }

    RDouble prl = GlobalDataBase::GetDoubleParaFromDB("prl");
    RDouble prt = GlobalDataBase::GetDoubleParaFromDB("prt");

    RDouble oprl = 1.0 / prl;
    GlobalDataBase::UpdateData("oprl", &oprl, PHDOUBLE, 1 );
    RDouble oprt = 1.0 / prt;
    GlobalDataBase::UpdateData("oprt", &oprt, PHDOUBLE, 1 );
}

void PerfectGas::SetAirInformationByDataBase()
{
    RDouble height = GlobalDataBase::GetDoubleParaFromDB("height");

    RDouble refDimensionalTemperature;
    RDouble refDimensionalPressure;
    RDouble refDimensionalDensity;
    RDouble refDimensionalSonicSpeed;

    GetAirInfo(height, refDimensionalTemperature, refDimensionalPressure, refDimensionalDensity, refDimensionalSonicSpeed);

    referenceParameterDimensional.SetDensity(refDimensionalDensity);
    referenceParameterDimensional.SetPressure(refDimensionalPressure);
    referenceParameterDimensional.SetTemperature(refDimensionalTemperature);
    referenceParameterDimensional.SetSoundVelocity(refDimensionalSonicSpeed);

    GlobalDataBase::UpdateData("refDimensionalDensity", &refDimensionalDensity, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("refDimensionalPressure", &refDimensionalPressure, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("refDimensionalTemperature", &refDimensionalTemperature, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("refDimensionalSonicSpeed", &refDimensionalSonicSpeed, PHDOUBLE, 1);
}

void PerfectGas::SetAirInformationByExpData()
{
    RDouble refDimensionalTemperature;
    RDouble refDimensionalPressure;
    RDouble refDimensionalDensity;
    RDouble refDimensionalSonicSpeed;

    RDouble refMachNumber = GlobalDataBase::GetDoubleParaFromDB("refMachNumber");
    RDouble T0 = referenceParameterDimensional.GetTemperature();
    RDouble P0 = referenceParameterDimensional.GetPressure();

    RDouble refGama = GlobalDataBase::GetDoubleParaFromDB("refGama");
    RDouble R = 287.053;

    RDouble mi   = refGama / (refGama - 1.0);
    RDouble coef = 1.0 + (refGama - 1.0) * refMachNumber * refMachNumber / 2.0;

    RDouble T, P, rho, c;
    T = T0 / coef;
    P = P0 / pow(coef, mi);
    rho = P / R / T;
    c = sqrt(refGama * R * T);

    refDimensionalDensity = rho;
    refDimensionalPressure = P;
    refDimensionalTemperature = T;
    refDimensionalSonicSpeed = c;

    referenceParameterDimensional.SetDensity(refDimensionalDensity);
    referenceParameterDimensional.SetPressure(refDimensionalPressure);
    referenceParameterDimensional.SetTemperature(refDimensionalTemperature);
    referenceParameterDimensional.SetSoundVelocity(refDimensionalSonicSpeed);

    GlobalDataBase::UpdateData("refDimensionalDensity", &refDimensionalDensity, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("refDimensionalPressure", &refDimensionalPressure, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("refDimensionalTemperature", &refDimensionalTemperature, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("refDimensionalSonicSpeed", &refDimensionalSonicSpeed, PHDOUBLE, 1);
}

void PerfectGas::GetAirInfo(const RDouble &h1, RDouble &t, RDouble &p, RDouble &den, RDouble &a)
{//USSA 76 Model for atmosphere.
    RDouble h      = h1 * 1000.0;
    RDouble r      = 287.053;
    RDouble g0     = 9.80665;
    RDouble rp     = 6.37111e6;
    RDouble g      = (rp/(rp+h))*(rp/(rp+h))*g0;
    RDouble t0     = 288.15;
    RDouble p0     = 10.1325e2;
    RDouble rho0   = 1.225;
    RDouble t11    = 216.65;
    RDouble p11    = 2.2632e2;
    RDouble rho11  = 3.6392e-1;
    RDouble t20    = t11;
    RDouble p20    = 5.4747e1;
    RDouble rho20  = 8.8035e-2;
    RDouble t32    = 228.65;
    RDouble p32    = 8.6789;
    RDouble rho32  = 1.3225e-2;
    RDouble t47    = 270.65;
    RDouble p47    = 1.1090;
    RDouble rho47  = 1.4275e-3;
    RDouble t52    = t47;
    RDouble p52    = 5.8997e-1;
    RDouble rho52  = 7.5943e-4;
    RDouble t61    = 252.65;
    RDouble p61    = 1.8209e-1;
    RDouble rho61  = 2.5109e-4;
    RDouble t79    = 180.65;
    RDouble p79    = 1.0376e-2;
    RDouble rho79  = 2.0010e-5;
    RDouble t90    = t79;
    RDouble p90    = 1.6437e-3;
    RDouble rho90  = 3.4165e-6;
    RDouble t100   = 210.02;
    RDouble p100   = 3.0070e-4;
    RDouble rho100 = 5.6044e-7;
    RDouble t110   = 257.00;
    RDouble p110   = 7.3527e-5;
    RDouble rho110 = 9.7081e-8;
    RDouble t120   = 349.49;
    RDouble p120   = 2.5209e-5;
    RDouble rho120 = 2.2222e-8;
    RDouble t150   = 892.79;
    RDouble p150   = 5.0599e-6;
    RDouble rho150 = 2.0752e-9;
    RDouble t160   = 1022.20;
    RDouble p160   = 3.6929e-6;
    RDouble rho160 = 1.2336e-9;
    RDouble t170   = 1103.40;
    RDouble p170   = 2.7915e-6;
    RDouble rho170 = 7.8155e-10;
    RDouble t190   = 1205.40;
    RDouble p190   = 1.6845e-6;
    RDouble rho190 = 3.5807e-10;
    RDouble t230   = 1322.30;
    RDouble p230   = 6.7138e-7;
    RDouble rho230 = 1.5640e-10;
    RDouble t300   = 1432.10;
    RDouble p300   = 1.8828e-7;
    RDouble rho300 = 1.9159e-11;
    RDouble t400   = 1487.40;
    RDouble p400   = 4.0278e-8;
    RDouble rho400 = 2.8028e-12;
    RDouble t500   = 1499.20;
    RDouble p500   = 1.0949e-8;
    RDouble rho500 = 5.2148e-13;
    RDouble t600   = 1506.10;
    RDouble p600   = 3.4475e-9;
    RDouble rho600 = 1.1367e-13;
    RDouble t700   = 1507.60;
    RDouble p700   = 1.1908e-9;
    RDouble rho700 = 1.5270e-13;

    RDouble rho = 0.0;

    if (h <= 11019.0) 
    {
        RDouble al1 = ( t11 - t0 )/11019.0;
        t   = t0 + al1 * h;
        p   = p0 * pow( t/t0, -g/(r * al1 ) );
        rho = rho0 * pow( t/t0, -1.0 - g/(r*al1) );
    }
    else if (h <= 20063.0)
    {
        t   = t11;
        p   = p11   * exp( -g*(h-11019.0)/(r*t11) );
        rho = rho11 * exp( -g*(h-11019.0)/(r*t11) );
    }
    else if (h <= 32162.0)
    {
        RDouble al2 = (t32-t20)/(32162.0-20063.0);
        t   = t11+al2*(h-20063.0);
        p   = p20 * pow(t/t11, -g/(r*al2));
        rho = rho20 * pow(t/t11, -1.0-g/(r*al2));
    }
    else if (h <= 47350.0)
    {
        RDouble al3 = (t47-t32)/(47350.0-32162.0);
        t   = t32+al3*(h-32162.0);
        p   = p32 * pow(t/t32, -g/(r*al3));
        rho = rho32 * pow(t/t32, -1.0-g/(r*al3));
    }

    else if (h <= 52429.0)
    {
        t   = t47;
        p   = p47*exp(-g*(h-47350.0)/(r*t47));
        rho = rho47*exp(-g*(h-47350.0)/(r*t47));
    }
    else if(h <= 61591.0)
    {
        RDouble al4 =( t61-t52)/(61591.0-52429.0);
        t = t47+al4*(h-52429.0);
        p = p52 * pow(t/t47, -g/(r*al4));
        rho = rho52* pow(t/t47,-1.0-g/(r*al4));
    }
    else if(h <= 79994.0)
    {
        RDouble al5=(t79-t61)/(79994.0-61591.0);
        t=t61+al5*(h-61591.0);
        p=p61*pow(t/t61, -g/(r*al5));
        rho=rho61*pow(t/t61,-1.0-g/(r*al5));
    }
    else if(h <= 90000.0)
    {
        t=t79;
        p=p79*exp(-g*(h-79994.0)/(r*t79));
        rho=rho79*exp(-g*(h-79994.0)/(r*t79));
    }
    else if(h <= 100000.0)
    {
        RDouble  al6=(t100-t90)/10000.0;
        t=t79+al6*(h-90000.0);
        p=p90* pow(t/t79,-g/(r*al6));
        rho=rho90* pow(t/t79,-1.0-g/(r*al6));
    }
    else if(h <= 110000.0)
    {
        RDouble al7=(t110-t100)/10000.0;
        t=t100+al7*(h-100000.0);
        p=p100*pow(t/t100, -g/(r*al7));
        rho=rho100*(t/t100, -1.0-g/(r*al7));
    }
    else if(h <= 120000.0)
    {
        RDouble al8=(t120-t110)/10000.0;
        t=t110+al8*(h-110000.0);
        p=p110*pow(t/t110, -g/(r*al8));
        rho=rho110*pow(t/t110, -1.0-g/(r*al8));
    }
    else if(h <= 150000.0)
    {
        RDouble al9=(t150-t120)/30000.0;
        t=t120+al9*(h-120000.0);
        p=p120*pow(t/t120,-g/(r*al9));
        rho=rho120*pow(t/t120,-1.0-g/(r*al9));
    }
    else if(h <= 160000.0)
    {
        RDouble al10=(t160-t150)/10000.0;
        t=t150+al10*(h-150000.0);
        p=p150*pow(t/t150, -g/(r*al10));
        rho=rho150*pow(t/t150, -1.0-g/(r*al10));
    }
    else if(h <= 170000.0)
    {
        RDouble al11=(t170-t160)/10000.0;
        t=t160+al11*(h-160000.0);
        p=p160*pow(t/t160, -g/(r*al11));
        rho=rho160*pow(t/t160,-1.0-g/(r*al11));
    }
    else if(h <= 190000.0)
    {
        RDouble al12=(t190-t170)/20000.0;
        t=t170+al12*(h-170000.0);
        p=p170*pow(t/t170, -g/(r*al12));
        rho=rho170*pow(t/t170,-1.0-g/(r*al12));
    }
    else if(h <= 230000.0)
    {
        RDouble al13=(t230-t190)/40000.0;
        t=t190+al13*(h-190000.0);
        p=p190*pow(t/t190,-g/(r*al13));
        rho=rho190*pow(t/t190,-1.0-g/(r*al13));
    }
    else if(h <= 300000.0)
    {
        RDouble al14=(t300-t230)/70000.0;
        t=t230+al14*(h-230000.0);
        p=p230*pow(t/t230, -g/(r*al14));
        rho=rho230*pow(t/t230, -1.0-g/(r*al14));
    }
    else if(h <= 400000.0)
    {
        RDouble al15=(t400-t300)/100000.0;
        t=t300+al15*(h-300000.0);
        p=p300*pow(t/t300, -g/(r*al15));
        rho=rho300*pow(t/t300, -1.0-g/(r*al15));
    }
    else if(h <= 500000.0)
    {
        RDouble al16=(t500-t400)/100000.0;
        t=t400+al16*(h-400000.0);
        p=p400*pow(t/t400, -g/(r*al16));
        rho=rho400*pow(t/t400, -1.0-g/(r*al16));
    }
    else if(h <= 600000.0)
    {
        RDouble al17=(t600-t500)/100000.0;
        t=t500+al17*(h-500000.0);
        p=p500*pow(t/t500, -g/(r*al17));
        rho=rho500*pow(t/t500, -1.0-g/(r*al17));
    }
    else if(h <= 700000.0)
    {
        RDouble al18=(t700-t600)/100000.0;
        t=t600+al18*(h-600000.0);
        p=p600*pow(t/t600, -g/(r*al18));
        rho=rho600*pow(t/t600, -1.0-g/(r*al18));
    }
    else
    {
        TK_Exit::ExceptionExit("Error: height must be less than 700km!!!\n", true);
    }

    a   = sqrt(1.4*r*t);
    p   = p * 100.0;
    den = rho;
}

void PerfectGas::ComputeReferenceParameter()
{
    //!The first two variables is temperature independent, 
    //!depend on molecular weight and mass fraction of species if chemical reactions exist.
    ComputeReferenceMolecularInformation();


    ComputeReferenceGeneralGasConstant();

    //!The following five variables depend on temperature, 
    //!and mass fraction of species if chemical reactions exist.
    ComputeReferenceViscosityDimensional();

    ComputeReferenceSpecificHeatRatio();

    ComputeReferenceSoundVelocity();

    ComputeReferenceVelocity();

    ComputeCoefficientOfStateEquation();
    //!The above seven variables should be computed firstly.

    //!The variables in this function are computed by different methods according to "inflowParaType", 
    //!including density, pressure or reynolds.
    ComputeReferenceGasInformation();

    ComputeReferencePrimitive();
}

void PerfectGas::ComputeReferencePrimitive()
{
    RDouble refMachNumber = referenceParameterDimensional.GetMachNumber();
    RDouble reference_sound_velocity, reference_pressure;

    RDouble massoo_ref = 1.0;

    RDouble reference_density     = 1.0;
    RDouble reference_temperature = 1.0;

    reference_sound_velocity = 1.0 / refMachNumber;
    
    reference_pressure = coefficientOfStateEquation * 
        reference_density * reference_temperature / massoo_ref;

    RDouble attackd = GlobalDataBase::GetDoubleParaFromDB("attackd");
    RDouble angleSlide = GlobalDataBase::GetDoubleParaFromDB("angleSlide");

    RDouble attack   = attackd * PI / 180.0;
    RDouble sideslip = angleSlide * PI / 180.0;

    GlobalDataBase::UpdateData("attack"  ,&attack  ,PHDOUBLE,1);
    GlobalDataBase::UpdateData("sideslip",&sideslip,PHDOUBLE,1);

    using namespace IDX;
    int nLaminar = GlobalDataBase::GetIntParaFromDB("nl");
    RDouble *prim_inf = new RDouble[nLaminar];
    GlobalDataBase::UpdateDataPtr("prim_inf", prim_inf);

    prim_inf[IR] = reference_density;
    prim_inf[IU] = 1.0 * cos(attack) * cos(sideslip);
    prim_inf[IV] = 1.0 * sin(attack) * cos(sideslip);
    prim_inf[IW] = 1.0 *               sin(sideslip);
    prim_inf[IP] = reference_pressure;
}


void PerfectGas::ComputeReferenceGasInformation()
{
    int inflowParaType = GlobalDataBase::GetIntParaFromDB( "inflowParaType" );

    const int NONDIMENSIONCONDITION = 0;
    const int FLIGHTCONDITION       = 1;
    const int EXPERIMENTCONDITION   = 2;
    const int SUBSONICBCCONDITION   = 3;

    if ( inflowParaType == NONDIMENSIONCONDITION )    //!Obtain the referce density and referce pressure for nondimension condition.
    {
        ComputeDensityWithReynoldsNumber();
        ComputePressureInGasStateEquation();
    }
    else if ( inflowParaType == FLIGHTCONDITION || inflowParaType == EXPERIMENTCONDITION || inflowParaType == SUBSONICBCCONDITION )  //obtain the referce density for nondimension, and the Re.
    {
        NormalizeAirInformation();
        ComputeReferenceReynoldsNumber();

        if ( inflowParaType == SUBSONICBCCONDITION )
        {
            GetSubsonicBCData();
        }
    }
    else
    {
        TK_Exit::UnexpectedVarValue( "inflowParaType", inflowParaType );
    }
}

void PerfectGas::ComputeReferenceReynoldsNumber()
{
    RDouble refDimensionalDensity = referenceParameterDimensional.GetDensity();
    RDouble reference_velocity_dimensional  = referenceParameterDimensional.GetVelocity();
    RDouble reference_viscosity_dimensional = referenceParameterDimensional.GetViscosity();

    RDouble refReNumber = refDimensionalDensity * reference_velocity_dimensional / reference_viscosity_dimensional;
    RDouble oreynolds = 1.0 / refReNumber;

    referenceParameterDimensional.SetReynoldsNumber(refReNumber);
    GlobalDataBase::UpdateData("refReNumber", &refReNumber, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("oreynolds", &oreynolds, PHDOUBLE, 1);
}

void PerfectGas::NormalizeAirInformation()
{
    RDouble refDimensionalPressure = referenceParameterDimensional.GetPressure();
    RDouble refDimensionalTemperature = referenceParameterDimensional.GetTemperature();
    RDouble reference_average_generanl_gas_constant = referenceParameterDimensional.GetAverageGeneralGasConstant();

    RDouble refDimensionalDensity = refDimensionalPressure / reference_average_generanl_gas_constant / refDimensionalTemperature;

    referenceParameterDimensional.SetDensity(refDimensionalDensity);
    GlobalDataBase::UpdateData("refDimensionalDensity", &refDimensionalDensity, PHDOUBLE, 1);
}

void PerfectGas::ComputePressureInGasStateEquation()
{
    RDouble refDimensionalDensity = referenceParameterDimensional.GetDensity();
    RDouble reference_average_generanl_gas_constant = referenceParameterDimensional.GetAverageGeneralGasConstant();
    RDouble refDimensionalTemperature = referenceParameterDimensional.GetTemperature();

    RDouble refDimensionalPressure = refDimensionalDensity * reference_average_generanl_gas_constant * refDimensionalTemperature;

    referenceParameterDimensional.SetPressure(refDimensionalPressure);
    GlobalDataBase::UpdateData("refDimensionalPressure", &refDimensionalPressure, PHDOUBLE, 1);
}

void PerfectGas::ComputeDensityWithReynoldsNumber()
{
    RDouble refReNumber = referenceParameterDimensional.GetReynoldsNumber();
    RDouble reference_viscosity_dimensional = referenceParameterDimensional.GetViscosity();
    RDouble reference_velocity_dimensional = referenceParameterDimensional.GetVelocity();

    RDouble refDimensionalDensity = refReNumber * reference_viscosity_dimensional / reference_velocity_dimensional;

    referenceParameterDimensional.SetDensity(refDimensionalDensity);
    GlobalDataBase::UpdateData("refDimensionalDensity", &refDimensionalDensity, PHDOUBLE, 1);
}

void PerfectGas::ComputeCoefficientOfStateEquation()
{
    //! No matter whether chemical reactions exist or not.
    RDouble gama0 = referenceParameterDimensional.GetAverageSpecificHeatRatio();
    RDouble refMachNumber = referenceParameterDimensional.GetMachNumber();

    //! The variable is the non-dimensional value of 8.314/M£¬M denotes the non-dimensioanl molecular weight of the mixture gas.
    coefficientOfStateEquation = 1.0 / (gama0 * refMachNumber * refMachNumber); //! universal gas constant.
    GlobalDataBase::UpdateData("coefficientOfStateEquation", &coefficientOfStateEquation, PHDOUBLE, 1);
}

void PerfectGas::ComputeReferenceVelocity()
{
    RDouble refMachNumber = referenceParameterDimensional.GetMachNumber();
    RDouble refDimensionalSonicSpeed = referenceParameterDimensional.GetSoundVelocity();

    RDouble reference_velocity_dimensional = refDimensionalSonicSpeed * refMachNumber;
    referenceParameterDimensional.SetVelocity(reference_velocity_dimensional);
}

void PerfectGas::ComputeReferenceSoundVelocity()
{
    RDouble gama0 = referenceParameterDimensional.GetAverageSpecificHeatRatio();
    RDouble reference_general_gas_constant = referenceParameterDimensional.GetAverageGeneralGasConstant();
    RDouble refDimensionalTemperature = referenceParameterDimensional.GetTemperature();

    //! unit: m/s.
    RDouble refDimensionalSonicSpeed = sqrt(gama0 * reference_general_gas_constant * refDimensionalTemperature);

    referenceParameterDimensional.SetSoundVelocity(refDimensionalSonicSpeed);
    GlobalDataBase::UpdateData("refDimensionalSonicSpeed", &refDimensionalSonicSpeed, PHDOUBLE, 1);
}

void PerfectGas::ComputeReferenceSpecificHeatRatio()
{
    RDouble refGama = 1.4;
    refGama = GlobalDataBase::GetDoubleParaFromDB("refGama");
    referenceParameterDimensional.SetAverageSpecificHeatRatio(refGama);
}

void PerfectGas::ComputeReferenceViscosityDimensional()
{
    RDouble t0i  = 288.15;
    RDouble tsi  = 110.4;
    RDouble mu0i = 1.7894E-05;

    RDouble refDimensionalTemperature = referenceParameterDimensional.GetTemperature();
    RDouble reference_viscosity_dimensional = (t0i + tsi) / (refDimensionalTemperature + tsi) * pow(refDimensionalTemperature / t0i, 1.5) * mu0i;
    referenceParameterDimensional.SetViscosity(reference_viscosity_dimensional);

    RDouble tsuth = tsi / refDimensionalTemperature;
    GlobalDataBase::UpdateData("tsuth", &tsuth, PHDOUBLE, 1);
}

void PerfectGas::ComputeReferenceMolecularInformation()
{
    //!unit: kg/mol.
    RDouble reference_average_molecular_weight_dimensional = 28.9644 * 1.0e-3;
    referenceParameterDimensional.SetAverageMolecularWeight( reference_average_molecular_weight_dimensional );

}

void PerfectGas::ComputeReferenceGeneralGasConstant()
{
    RDouble reference_average_molecular_weight_dimensional   = referenceParameterDimensional.GetAverageMolecularWeight();
    RDouble reference_average_general_gas_constant          = rjmk / reference_average_molecular_weight_dimensional;

    referenceParameterDimensional.SetAverageGeneralGasConstant( reference_average_general_gas_constant );
}

void PerfectGas::ComputeOtherProperParameterForGrid()
{
    RDouble reynoldsReferenceLengthDimensional = GlobalDataBase::GetDoubleParaFromDB( "reynoldsReferenceLengthDimensional" );

    RDouble forceRefenenceLength = GlobalDataBase::GetDoubleParaFromDB("forceRefenenceLength");
    RDouble forceRefenenceLengthSpanWise = GlobalDataBase::GetDoubleParaFromDB("forceRefenenceLengthSpanWise");
    RDouble forceRefenenceArea = GlobalDataBase::GetDoubleParaFromDB("forceRefenenceArea");
    RDouble TorqueRefX = GlobalDataBase::GetDoubleParaFromDB("TorqueRefX");
    RDouble TorqueRefY = GlobalDataBase::GetDoubleParaFromDB("TorqueRefY");
    RDouble TorqueRefZ = GlobalDataBase::GetDoubleParaFromDB("TorqueRefZ");

    if ( GetDim() == THREE_D )
    {
        forceRefenenceArea /= reynoldsReferenceLengthDimensional*reynoldsReferenceLengthDimensional;
    }
    else if ( GetDim() == TWO_D )
    {
        forceRefenenceArea /= reynoldsReferenceLengthDimensional;
    }
    forceRefenenceLength /= reynoldsReferenceLengthDimensional;
    forceRefenenceLengthSpanWise /= reynoldsReferenceLengthDimensional;
    TorqueRefX /= reynoldsReferenceLengthDimensional;
    TorqueRefY /= reynoldsReferenceLengthDimensional;
    TorqueRefZ /= reynoldsReferenceLengthDimensional;

    GlobalDataBase::UpdateData("forceRefenenceLength", &forceRefenenceLength, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("forceRefenenceLengthSpanWise", &forceRefenenceLengthSpanWise, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("forceRefenenceArea", &forceRefenenceArea, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("TorqueRefX", &TorqueRefX, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("TorqueRefY", &TorqueRefY, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("TorqueRefZ", &TorqueRefZ, PHDOUBLE, 1);

    int isUnsteady = GlobalDataBase::GetIntParaFromDB( "iunsteady" );
    if ( isUnsteady == 1 )
    {
        PrintToWindow( "Unsteady Informations" , "\n" );

        RDouble physicalTimeStep               = GlobalDataBase::GetDoubleParaFromDB( "physicalTimeStep" );
        RDouble reference_velocity_dimensional = referenceParameterDimensional.GetVelocity();

        RDouble unit = reynoldsReferenceLengthDimensional / reference_velocity_dimensional;

        RDouble physicalTimeStep_dimensional = physicalTimeStep * unit;

        PrintToWindow( "  Physical Time Step       :    ", physicalTimeStep_dimensional , "s", "\n" );
    }
}

RDouble PerfectGas::GetUniversalGasConstant()
{
    //! Return the non-dimensional universal gas constant.
    return coefficientOfStateEquation;
}
}

}
