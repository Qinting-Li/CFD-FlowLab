#include "Param_TurbSolver.h"
#include "Transition.h"
#include "GlobalDataBase.h"
#include "Constants.h"

namespace PHSPACE
{

LIB_EXPORT Param_TurbSolver::Param_TurbSolver()
{

}

LIB_EXPORT Param_TurbSolver::~Param_TurbSolver()
{
    if (freeStreamTurbVar != NULL)
    {
        delete []freeStreamTurbVar;
        freeStreamTurbVar = NULL;
    }
}

LIB_EXPORT void Param_TurbSolver::Init()
{
    Param_CFDSolver::Init();

    turbCFLScale = GlobalDataBase::GetDoubleParaFromDB("turbCFLScale");

    DESType = GlobalDataBase::GetIntParaFromDB("DESType");

    transitionType = GlobalDataBase::GetIntParaFromDB("transitionType");

    SAProductType = GlobalDataBase::GetIntParaFromDB("SAProductType");

    SSTProductType = GlobalDataBase::GetIntParaFromDB("SSTProductType");

    nTurbulenceEquation = GlobalDataBase::GetIntParaFromDB("n_turb");

    eddyViscosityLimit = GlobalDataBase::GetDoubleParaFromDB("eddyViscosityLimit");

    wallFunctionType = GlobalDataBase::GetIntParaFromDB("wallFunctionType");

    SST_a1 = 0.31;

    freeStreamViscosity = GlobalDataBase::GetDoubleParaFromDB("freeStreamViscosity");

    using namespace IDX;
    RDouble * prim_inf = reinterpret_cast< RDouble * > (GlobalDataBase::GetDataPtr("prim_inf" ));
    RDouble reference_density_farfield = prim_inf[IR];

    RDouble refReNumber = 1.0e5;
    GlobalDataBase::GetData("refReNumber", &refReNumber, PHDOUBLE, 1);
   
    RDouble muoo = 1.0e-1;
    GlobalDataBase::GetData("muoo", &muoo, PHDOUBLE, 1);
    RDouble kwoo = 1.0e-8;
    GlobalDataBase::GetData("kwoo", &kwoo, PHDOUBLE, 1);

    int viscousType = 1;
    GlobalDataBase::GetData("viscousType", &viscousType, PHINT, 1);

    int neasm = -1;
    GlobalDataBase::GetData("neasm", &neasm, PHINT, 1);

    string viscousName;
    GlobalDataBase::GetData("viscousName", &viscousName, PHSTRING, 1);

    if (viscousName.substr(0, 12) == "easm-kw-2001" || viscousName.substr(0, 12) == "easm-kw-2003")
    {
        neasm = 2;
    }
    else if (viscousName.substr(0, 12) == "easm-kw-2005")
    {
        neasm = 3;
    }
    GlobalDataBase::UpdateData("neasm",&neasm, PHINT, 1);

    int nolstress, nrokplus;

    if (viscousType >= TWO_EQU && neasm > 0)
    {
        nolstress = 1;
        nrokplus  = 1;
    }
    else
    {
        nolstress = - 1;
        nrokplus  = - 1;
    }
    GlobalDataBase::UpdateData("nolstress",&nolstress, PHINT, 1);
    GlobalDataBase::UpdateData("nrokplus" ,&nrokplus , PHINT, 1);

    RDouble turb_fbeta, turb_fbetas, turb_cmu;

    RDouble SA_cv1;

    //! clz begin
    SSTProductLimit = 10.0;    //! Refer to the version of 2003.

    turb_fbeta  = 1.0;
    turb_fbetas = 1.0;
    SST_betaStar = 0.09;
    //! In fact£¬turb_cmu should be 0.09£¬but it is necessary because turb_cmu is generally absorbed.
    turb_cmu = 1.0;

    if (viscousName.substr(0,6) == "1eq-sa")
    {
        SA_cv1    = 7.1;
        KW_sigma  = 2.0/3.0;

        SA_cv1_cube = SA_cv1 * SA_cv1 * SA_cv1;
    }
    else if (viscousName.substr(0,6) == "2eq-kw")
    {
        if (viscousName.substr(0,18) == "2eq-kw-wilcox-1988")
        {
            KW_sigmaK   = 0.5;
            KW_sigmaW   = 0.5;
            SST_betaStar = 1.0;
            turb_cmu    = 0.09;

            if (neasm == 1)
            {
                KW_sigmaK   = 0.5;
                KW_sigmaW   = 0.5;
                SST_betaStar = 1.0;
                turb_cmu    = 0.09;
            }
        }
        else if (viscousName.substr(0,18) == "2eq-kw-wilcox-1993")
        {
            KW_sigmaK   = 1.0;
            KW_sigmaW   = 0.6;
        }
        else if (viscousName.substr(0,18) == "2eq-kw-wilcox-1998")
        {
            KW_sigmaK   = 0.5;
            KW_sigmaW   = 0.5;

            if (neasm == 1)
            {
                KW_sigmaK   = 0.5;
                KW_sigmaW   = 0.5;
                SST_betaStar = 1.0;
                turb_cmu    = 0.09;
            }
        }
        else if (viscousName.substr(0,18) == "2eq-kw-wilcox-2006")
        {
            KW_sigmaK   = 0.6;
            KW_sigmaW   = 0.5;
        }
        else if (viscousName.substr(0,14) == "2eq-kw-kok-tnt")
        {
            KW_sigmaK   = 0.667;
            KW_sigmaW   = 0.5;

            if (neasm == 1)
            {
                KW_sigmaK   = 1.0 / 1.4;
                KW_sigmaW   = 1.0 / 2.2;
                SST_betaStar = 1.0;
                turb_cmu    = 0.09;
            }
        }
        else if (viscousName.substr(0,17) == "2eq-kw-menter-sst")
        {
            KW_sigmaK1  = 0.85;
            KW_sigmaW1  = 0.5;
            SST_alphaw1 = 0.555;     //! = 5/9  refer to the version of SST-2003
            SST_beta1   = 0.075;

            KW_sigmaK2  = 1.0;
            KW_sigmaW2  = 0.856;
            SST_alphaw2 = 0.44;
            SST_beta2   = 0.0828;
            
            if (transitionType == IREGAMA)
            {
                transition_ce1 = 1.0;
                transition_ca1 = 2.0;
                transition_ce2 = 50.0;
                transition_ca2 = 0.06;
                transition_dct = 2.0;
                transition_df  = 1.0;
                transition_cct = 0.03;
                transition_s1  = 2.0;
                transition_ccf = 0.6;
            }
        }
        else if (viscousName.substr(0,17) == "2eq-kw-menter-bsl")
        {
            KW_sigmaK1  = 0.5;
            KW_sigmaW1  = 0.5;
            SST_alphaw1 = 0.5532;
            SST_beta1   = 0.075;

            KW_sigmaK2  = 1.0;
            KW_sigmaW2  = 0.856;
            SST_alphaw2 = 0.4404;
            SST_beta2   = 0.0828;
        }
        else if (viscousName.substr(0,10) == "2eq-kw-hyb")
        {
        }
        else if (viscousName.substr(0,10) == "2eq-kw-bdp")
        {
        }
    }
    else if (viscousName.substr(0,4) == "easm")
    {
        //! It is a specific easm, not a random collected one.
        if (viscousName.substr(0,12) == "easm-kw-2001")
        {
            KW_sigmaK   = 0.5;
            KW_sigmaW   = 0.5339;
            SST_betaStar = 1.0;
            turb_cmu    = 0.0895;
        }
        else if (viscousName.substr(0,12) == "easm-kw-2003")
        {
            KW_sigmaK   = 1.0;
            KW_sigmaW   = 0.5339;
            SST_betaStar = 1.0;
            turb_cmu    = 0.0895;
        }
        else if (viscousName.substr(0,12) == "easm-kw-2005")
        {
            KW_sigmaK1  = 1.1;
            KW_sigmaW1  = 0.53;
            SST_alphaw1 = 0.518;
            SST_beta1   = 0.0747;

            KW_sigmaK2  = 1.1;
            KW_sigmaW2  = 1.0;
            SST_alphaw2 = 0.44;
            SST_beta2   = 0.0828;

            SST_betaStar = 0.09;
        }
    }

    GlobalDataBase::UpdateData("turb_fbeta", &turb_fbeta, PHDOUBLE, 1);
    GlobalDataBase::UpdateData("turb_fbetas", &turb_fbetas, PHDOUBLE, 1);

    freeStreamTurbVar = new RDouble [nTurbulenceEquation];

    //! This code should be check for structured grid.
    if (viscousType == ONE_EQU)
    {
        RDouble muoo = GlobalDataBase::GetDoubleParaFromDB("muoo");
        freeStreamTurbVar[0] = muoo / (reference_density_farfield * turb_cmu + SMALL);
    }
    else if (viscousType == TWO_EQU)
    {
        RDouble kwoo = GlobalDataBase::GetDoubleParaFromDB("kwoo") / refReNumber;
        RDouble mutoo = 0.01;
        RDouble keoo = mutoo * kwoo / (turb_cmu * reference_density_farfield);
        freeStreamTurbVar[IKE] = keoo;
        freeStreamTurbVar[IKW] = kwoo;

        if (transitionType == IREGAMA)
        {
            RDouble TUoo = 0.18;
            RDouble Fcta = 1.0;
            freeStreamTurbVar[IGAMA] = 1.0;

            mutoo = GlobalDataBase::GetDoubleParaFromDB("freeStreamViscosity");
            turbulenceIntensity = GlobalDataBase::GetDoubleParaFromDB("turbIntensity");
            if (turbulenceIntensity >= SMALL)
            {
                keoo = 1.5 * (0.01 * turbulenceIntensity) * (0.01 * turbulenceIntensity);
                kwoo = keoo * (turb_cmu * reference_density_farfield) / mutoo;
                freeStreamTurbVar[IKE] = keoo;
                freeStreamTurbVar[IKW] = kwoo;
            }

            TUoo = TurbulenceIntensity(1.0, keoo);
            freeStreamTurbVar[IRECT] = EmpiricalCorrelationOfRectat(TUoo, Fcta);
        }
    }

    kindOfTurbSource = GlobalDataBase::GetIntParaFromDB("kindOfTurbSource");
    isWennScheme = GlobalDataBase::GetIntParaFromDB("isWennScheme");
}

}