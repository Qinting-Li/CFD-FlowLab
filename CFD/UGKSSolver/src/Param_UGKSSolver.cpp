#include "Param_UGKSSolver.h"
#include "GlobalDataBase.h"
#include "Constants.h"
#include "Pointer.h"
#include <cmath>

namespace PHSPACE
{
    Param_UGKSSolver::Param_UGKSSolver()
    {

    }

    Param_UGKSSolver::~Param_UGKSSolver()
    {
        DelPointer(primRef);
        DelPointer(us);
        DelPointer(vs);
        DelPointer(ws);
        DelPointer(wx);
        DelPointer(wy);
        DelPointer(wz);
    }

    void Param_UGKSSolver::CreateParaPtr()
    {
        us = new RDouble[uEnd + 1];
        vs = new RDouble[vEnd + 1];
        ws = new RDouble[wEnd + 1];
                             
        wx = new RDouble[uEnd + 1];
        wy = new RDouble[vEnd + 1];
        wz = new RDouble[wEnd + 1];
    }
    
    void Param_UGKSSolver::ReadParameter()
    {
        height = GlobalDataBase::GetDoubleParaFromDB("height");
        refGamma = GlobalDataBase::GetDoubleParaFromDB("refGama");
        refMachNumber = GlobalDataBase::GetDoubleParaFromDB("refMachNumber");
        refReNumber = GlobalDataBase::GetDoubleParaFromDB("refReNumber");
        attackd = GlobalDataBase::GetDoubleParaFromDB("attackd");
        angleSlide = GlobalDataBase::GetDoubleParaFromDB("angleSlide");
        refDimensionalTemperature = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");
        refDimensionalPressure = GlobalDataBase::GetDoubleParaFromDB("refDimensionalPressure");
        refDimensionalDensity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalDensity");
        refDimensionalVelocity  = GlobalDataBase::GetDoubleParaFromDB("reference_velocity_dimensional");
        wallTemperature = GlobalDataBase::GetDoubleParaFromDB("wallTemperature");
        prandtl  = GlobalDataBase::GetDoubleParaFromDB("prt");

        gridMagnificationInverse = GlobalDataBase::GetDoubleParaFromDB("forceReferenceLength");
        refLength = GlobalDataBase::GetDoubleParaFromDB("forceReferenceFLength");
        forceRefArea = GlobalDataBase::GetDoubleParaFromDB("forceReferenceArea");
        TorqueRefX = GlobalDataBase::GetDoubleParaFromDB("TorqueRefX");
        TorqueRefY = GlobalDataBase::GetDoubleParaFromDB("TorqueRefY");
        TorqueRefZ = GlobalDataBase::GetDoubleParaFromDB("TorqueRefZ");

        nStart = GlobalDataBase::GetIntParaFromDB("nstart");
        maxSimuStep = GlobalDataBase::GetIntParaFromDB("maxSimuStep");
        intervalStepFlow = GlobalDataBase::GetIntParaFromDB("intervalStepFlow");
        intervalStepPlot = GlobalDataBase::GetIntParaFromDB("intervalStepPlot");
        intervalStepForce = GlobalDataBase::GetIntParaFromDB("intervalStepForce");
        intervalStepRes = GlobalDataBase::GetIntParaFromDB("intervalStepRes");
        resSaveFile = GlobalDataBase::GetStrParaFromDB("resSaveFile");
        restartUGKSFile = GlobalDataBase::GetStrParaFromDB("restartNSFile");
        stopNow = GlobalDataBase::GetIntParaFromDB("stopnow");
        isDim = GlobalDataBase::GetIntParaFromDB("isdim");
        dataPos = GlobalDataBase::GetIntParaFromDB("data_pos");
        uEnd  = GlobalDataBase::GetIntParaFromDB("u_num");//
        vEnd  = GlobalDataBase::GetIntParaFromDB("v_num");
        wEnd  = GlobalDataBase::GetIntParaFromDB("w_num");
        neqn    = GlobalDataBase::GetIntParaFromDB("nm");
        uBegin  = 1;//
        vBegin  = 1;
        wBegin  = 1;
        nVelocity = (uEnd - uBegin + 1) *(vEnd - vBegin + 1)* (wEnd - wBegin + 1);
        GlobalDataBase::UpdateData("nVelocity",&nVelocity, PHDOUBLE, 1);
        cfl = GlobalDataBase::GetDoubleParaFromDB("cfl");
        nlhs = GlobalDataBase::GetIntParaFromDB("nlhs");
        cDom = GlobalDataBase::GetIntParaFromDB("cdom");

        RDouble timeMin = 1.e30;
        GlobalDataBase::UpdateData("timeMin",&timeMin, PHDOUBLE, 1);
    }
    
    void Param_UGKSSolver::SetPerfectgasParameter()
    {
        const int ns = 1;
        string *varName = new string[ns]();
        RDouble *cnInit = new RDouble[ns]();
        varName[0] = "Ar";
        cnInit[0] = 1;
        for(int i = 0; i < ns; i++)
        {
            cnInitMap.insert(pair<string,RDouble>(varName[i],cnInit[i]));
        }
        delete []  cnInit;    cnInit = NULL;
        delete []  varName;    varName = NULL;
    }
    
    void Param_UGKSSolver::SetDefaultParameter()
    {
        if(forceRefArea < 0)
        {
            cout << "The reference area(sref) can not be less than zero!" << endl;
        }
        if(gridMagnificationInverse < 0)
        {
            cout << "The reference length(lfref) can not be less than zero!" << endl;
        }
        if(refLength < 0)
        {
            cout << "The scale coefficient (lref) can not be less than zero!" << endl;
        }

        if (intervalStepFlow  < 0) intervalStepFlow  = 100;
        if (intervalStepPlot  < 0) intervalStepPlot  = 100;
        if (intervalStepForce < 0) intervalStepForce = 100;
        if (intervalStepRes   < 0) intervalStepRes   = 100;

        if(cfl < 1.0e-2)
            cfl = 0.01;
        else if(cfl > 1.0e+5)
            cfl = 5.0e3;

        if(nStart==0)
            nStart = 2;
        else if(nStart != 1 && nStart != 2 && nStart != 3)
            nStart = 0;

        int nOut = 0;
        RDouble wholeTime = 0.;
        if(nStart != 2)
        {
            nOut = 0;
            wholeTime = 0.0;
        }
        GlobalDataBase::UpdateData("nOut",&nOut, PHINT, 1);
        GlobalDataBase::UpdateData("wholeTime",&wholeTime, PHDOUBLE, 1);
    }

    void Param_UGKSSolver::SetScaleParameter()
    {
        RDouble unit = 1./gridMagnificationInverse; 
        forceRefArea *= unit * unit;
        refLength *= unit;
        TorqueRefX *= unit;
        TorqueRefY *= unit;
        TorqueRefZ *= unit;
    }
    
    void Param_UGKSSolver::ComputeInflowParameter()
    {
        //!initialize number of equation 
        int nm    = GlobalDataBase::GetIntParaFromDB("nm");
        int nl = nm;
        int neqn = nl;
        GlobalDataBase::UpdateData("nl",&nl, PHINT, 1);
        GlobalDataBase::UpdateData("neqn",&neqn, PHINT, 1);

        RDouble molarMassInverse = 1.0/0.03995;//! inverse of Ar's molar mass,unit is kg/mol.
        RDouble diameter = 4.17e-10;//! Ar's diameter,unit is m.
        RDouble mass = 6.63e-26;//! Ar's mass,unit is kg.
        RDouble R = 8.31434;//! universal gas constant,unit is J/(Mol*K).
        RDouble kBoltz = 1.380622e-23;//! Boltzmann constant£¬unit is J/(K).
        
        //RDouble attack   = attackd * PI / 180.0;
        RDouble soundVelocity = sqrt(refGamma * R * molarMassInverse * refDimensionalTemperature);
        if (refDimensionalVelocity < 0.0)
            refDimensionalVelocity = soundVelocity * refMachNumber;
        else
            refMachNumber = refDimensionalVelocity / soundVelocity;

        RDouble visloo = 5.0/16.0 * sqrt(mass * kBoltz * refDimensionalTemperature / PI) / (diameter * diameter);
        if (refReNumber < 0.0) 
            refReNumber = refDimensionalDensity * refDimensionalVelocity * gridMagnificationInverse * refLength / visloo;
        refReNumber = refReNumber / refLength;

        RDouble roo = 1.0;
        RDouble too = 1.0;
        RDouble uoo = 0.;//1.0 * cos(attack) * cos(sideslip);
        RDouble voo = 0;//1.0 * sin(attack) * cos(sideslip);
        RDouble woo = 0;//1.0 * sin(sideslip);
        RDouble poo = 1.0/(refGamma * refMachNumber * refMachNumber);
        GlobalDataBase::UpdateData("temperature_inflow",&too, PHDOUBLE, 1);

        using namespace IDX;
        primRef = new RDouble[neqn];
        primRef[IR] = roo;
        primRef[IU] = uoo;
        primRef[IV] = voo;
        primRef[IW] = woo;
        primRef[IP] = poo;
    }
    //GAUSS-HERMITE
     /*void Param_UGKSSolver::InitDiscreteVelocity(int uBegin,int uEnd,const RDouble uvwMin,const RDouble uvwMax,RDouble *uvw)
    {
        RDouble datar29[30] = {  0, -6.72869519860885, -5.99897128946382, -5.38964052196675, -4.84136365105916, -4.33147829381915, 
                                -3.84826679221362, -3.38464514109221, -2.93588250429013, -2.49858569101940, -2.07018107605343,
                                -1.64862291389232, -1.23221575508475, -0.81949868127091, -0.40916463639493,  0.00000000000000,
                                0.40916463639493,  0.81949868127091,  1.23221575508475,  1.64862291389232,  2.07018107605343,
                                2.49858569101940,  2.93588250429013,  3.38464514109221,  3.84826679221362,  4.33147829381915,
                                4.84136365105916,  5.38964052196675,  5.99897128946382,  6.72869519860885};
        for(int i = uBegin; i <= uEnd; ++ i)
        {
            RDouble coef = static_cast<RDouble>(uvwMax - uvwMin)/(2.0*fabs(datar29[1]));
            uvw[i] = datar29[i] * coef;
            //cout <<  uvw[i] <<endl;
        }
    }
    
    void Param_UGKSSolver::InitVelocityIntegrateCoef(int uBegin,int uEnd,const RDouble uvwMin,const RDouble uvwMax,RDouble *ws)
    {
         RDouble datar29[30] = {  0, -6.72869519860885, -5.99897128946382, -5.38964052196675, -4.84136365105916, -4.33147829381915, 
                                -3.84826679221362, -3.38464514109221, -2.93588250429013, -2.49858569101940, -2.07018107605343,
                                -1.64862291389232, -1.23221575508475, -0.81949868127091, -0.40916463639493,  0.00000000000000,
                                0.40916463639493,  0.81949868127091,  1.23221575508475,  1.64862291389232,  2.07018107605343,
                                2.49858569101940,  2.93588250429013,  3.38464514109221,  3.84826679221362,  4.33147829381915,
                                4.84136365105916,  5.38964052196675,  5.99897128946382,  6.72869519860885};
       RDouble dataw29[30] = {   0, 0.83939132420011,  0.65344293536984,  0.57349828373831,  0.52654651532625,  0.49509188270430,
                                0.47247844450338,  0.45553786063694,  0.44254485850123,  0.43247308063762,  0.42467426026128,
                                0.41872276848514,  0.41433372427760,  0.41131727936766,  0.40955227419762,  0.40897117463523,
                                0.40955227419762,  0.41131727936766,  0.41433372427760,  0.41872276848514,  0.42467426026128,
                                0.43247308063762,  0.44254485850123,  0.45553786063694,  0.47247844450338,  0.49509188270430,
                                0.52654651532625,  0.57349828373831,  0.65344293536984,  0.83939132420011};
       for(int i = uBegin; i <= uEnd; ++ i)
        {
            RDouble coef = static_cast<RDouble>(uvwMax - uvwMin)/(2.0*fabs(datar29[1]));
            ws[i] = dataw29[i] * coef;
            //cout <<  ws[i] <<endl;
        }
    }*/
    
    void Param_UGKSSolver::InitDiscreteVelocity(int uBegin,int uEnd,const RDouble uvwMin,const RDouble uvwMax,RDouble *uvw)
    {
        RDouble dv = (uvwMax - uvwMin)/(uEnd - uBegin);
        for(int i = uBegin; i <= uEnd; ++ i)
        {
            uvw[i] = uvwMin + (i - uBegin) * dv;
        }
    }
    
    void Param_UGKSSolver::InitVelocityIntegrateCoef(int uBegin,int uEnd,const RDouble uvwMin,const RDouble uvwMax,RDouble *ws)
    {
        RDouble dv = (uvwMax - uvwMin)/static_cast<RDouble>(uEnd - uBegin);
        ws[uBegin] = 14.0/45.0*dv;
        for(int j = uBegin + 1; j <= uEnd - 1;++ j)
        {
            ws[j] = 64.0/45.0*dv;
            if((j-1)%4 == 0)
            {
                ws[j] = 28.0/45.0 * dv;
            }
            else if((j-1)%4 == 2)
            {
                ws[j] = 24.0/45.0 * dv;
            }
        }
        ws[uEnd] = ws[uBegin];
    }
    
     void Param_UGKSSolver::InitDiscreteVelocityAndIntegrateCoef()
    {
        const RDouble velocityMax = 22.;
        const RDouble velocityMin = -22.;
        const RDouble uMin = velocityMin;
        const RDouble uMax = velocityMax;
        InitDiscreteVelocity(uBegin, uEnd, uMin, uMax, us);
        InitVelocityIntegrateCoef(uBegin, uEnd,uMin, uMax, wx);

        const RDouble vMin = velocityMin;
        const RDouble vMax = velocityMax;
        InitDiscreteVelocity(vBegin, vEnd, vMin, vMax, vs);
        InitVelocityIntegrateCoef(vBegin, vEnd, vMin, vMax, wy);

        const RDouble wMin = velocityMin;
        const RDouble wMax = velocityMax;
        InitDiscreteVelocity(wBegin, wEnd, wMin, wMax, ws);
        InitVelocityIntegrateCoef(wBegin, wEnd, wMin, wMax, wz);
    }

    
    void Param_UGKSSolver::Init()
    {
        ReadParameter();
        
        //SetPerfectgasParameter();

        SetDefaultParameter();

        SetScaleParameter();

        ComputeInflowParameter();
        
        CreateParaPtr();

        InitDiscreteVelocityAndIntegrateCoef();
    }
}
