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
//! @file      UGKSSimulation.h
//! @brief     It defines UGKS solver's parameters.
//! @author    

#pragma once
#include "LIB_Macro.h"
#include "TypeDefine.h"
namespace PHSPACE
{
    //! @brief Param_UGKSSolver defines UGKS solver's parameters.
    class Param_UGKSSolver
    {
    private:
        //! CFD Control Parameter.
        RDouble height;
        RDouble refGamma;
        RDouble refMachNumber;
        RDouble refReNumber;
        RDouble attackd;
        RDouble angleSlide;
        RDouble refDimensionalTemperature;
        RDouble refDimensionalPressure;
        RDouble refDimensionalDensity;
        RDouble refDimensionalVelocity;
        RDouble wallTemperature;
        RDouble prandtl;

        //! Force Parameter 
        RDouble gridMagnificationInverse;
        RDouble refLength;
        RDouble forceRefArea;
        RDouble TorqueRefX;
        RDouble TorqueRefY;
        RDouble TorqueRefZ;

        //! General Control Parameter.
        int nStart;
        int maxSimuStep;
        int intervalStepFlow;
        int intervalStepPlot;
        int intervalStepForce;
        int intervalStepRes;
        int stopNow;
        int isDim;
        int dataPos;
        int uEnd;
        int vEnd;
        int wEnd;
        int uBegin;
        int vBegin;
        int wBegin;
        int neqn;
        int nVelocity;

        //! Other Parameters
        RDouble cfl;
        RDouble globalMinTimeStep;
        int nlhs;
        int nLimiter;
        int cDom;
        //! PerfectGasMap
        RDoubleDBase cnInitMap;

        //! output
        string restartUGKSFile;
        string resSaveFile;
        string flowFieldTecFile;

        //! velocity space
        RDouble *us;
        RDouble *vs;
        RDouble *ws;
        RDouble *wx;
        RDouble *wy;
        RDouble *wz;
        RDouble *primRef;
    public:
        Param_UGKSSolver();
        ~Param_UGKSSolver();

    public:
        //! Init all parameters
        void Init();
        RDouble GetHeight() const;
        RDouble GetRefGamma() const;
        RDouble GetRefMachNumber() const;
        RDouble GetRefReNumber() const;
        RDouble GetAttackd() const;
        RDouble GetSideslipd() const;
        RDouble GetRefDimensionalTemperature() const;
        RDouble GetRefDimensionalPressure() const;
        RDouble GetRefDimensionalDensity() const;
        RDouble GetRefDimensionalVelocity() const;
        RDouble GetWallTemperature() const;
        RDouble GetPrandtl() const;

        RDouble GetForceRefLength() const;
        RDouble GetForceRefFLength() const;
        RDouble GetForceRefArea() const;
        RDouble GetTorqueRefX() const;
        RDouble GetTorqueRefY() const;
        RDouble GetTorqueRefZ() const;

        int GetNStart() const;
        int GetNumberOfEquations() const {return neqn;};
        int GetMaxSimuStep() const;

        //! Get interval steps of dumping restart
        int GetIntervalStepFlow() const;

        //! Get interval steps of dumping residuals
        int GetIntervalStepRes() const;

        //! Get interval steps of dumping flow field plot
        int GetIntervalStepPlot() const;

        //! Get interval steps of dumping force
        int GetIntervalStepForce() const;

        int GetStopNow() const;
        int GetIsDim() const;
        int GetDataPos() const;
        int GetUEnd() const;
        int GetVEnd() const;
        int GetWEnd() const;
        int GetUBegin() const;
        int GetVBegin() const;
        int GetWBegin() const;
        int GetNVelocity() const;

        RDouble GetCFL() const;
        RDouble GetGlobalMinTimeStep() const;
        int GetNlhs() const;
        int GetNLimiter() const;
        int GetCDom() const;

        const string & GetRestartUGKSFile() const;
        const string & GetResSaveFile() const;
        const string & GetFlowFieldTecFile() const;

        RDouble * GetDiscreteVelocityU() const;
        RDouble * GetDiscreteVelocityV() const;
        RDouble * GetDiscreteVelocityW() const;
        RDouble * GetVelocityIntegrateCoeX() const;
        RDouble * GetVelocityIntegrateCoeY() const;
        RDouble * GetVelocityIntegrateCoeZ() const;
        RDouble * GetPrimitiveVariable() const;
  private:
        //! read parameter.
        void ReadParameter();

        //! Set Perfect gas parameter.
        void SetPerfectgasParameter();

        //! Set defult parameter.
        void SetDefaultParameter();

        //! Set scale parameter.
        void SetScaleParameter();

         //! Compute inflow parameter.
        void ComputeInflowParameter();

        //! Init discrete velocity and integrate Coeffcient
        void InitDiscreteVelocityAndIntegrateCoef();

        //! Init discrete velocity
        //! @param[in] uBegin              minimum index of discrete velocity u
        //! @param[in] uEnd                maximum index of discrete velocity v
        //! @param[in] uvwMin              minimum discrete velocity u
        //! @param[in] uvwMax              maximum discrete velocity u
        //! @param[uvw] uvw                discrete velocity
        void InitDiscreteVelocity(int uBegin, int uEnd, const RDouble uvwMin, const RDouble uvwMax, RDouble *uvw);

        //! Init integrate Coeffcient
        //! @param[in] uBegin              minimum index of discrete velocity u
        //! @param[in] uEnd                maximum index of discrete velocity v
        //! @param[in] uvwMin              minimum discrete velocity u
        //! @param[in] uvwMax              maximum discrete velocity u
        //! @param[uvw] ws                 integrate Coeffcient
        void InitVelocityIntegrateCoef(int uBegin, int uEnd, const RDouble uvwMin, const RDouble uvwMax, RDouble *ws);

        //! Creat pointer of parameters
        void CreateParaPtr();
    };

    #include "Param_UGKSSolver.hxx"
}