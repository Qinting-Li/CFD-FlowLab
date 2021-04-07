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
//! @file      LES.h
//! @brief     LES Solver.
//! @author    ZhangZipei.

#pragma once
#include "TypeDefine.h"
#include "Complex.h"
#include "SpecDiffHybGrid.h"
#include "CompactDifferenceFirstDerivative.h"
using namespace std;

namespace PHSPACE
{

	class LES 
	{

	public:
		LES();
		~LES();

    public:
        Complex4D *complexDiffTau;
    
    private:
        RDouble1D *delta;
        RDouble1D *wallDistance, *wallDistancePlus;
        RDouble1D *vanDriestFunc;
        RDouble1D *dynamicCoef;
        RDouble1D *gkx, *gky;
        Complex4D *complexVarFilter, *testFilterGradUij, *complexTauij;
        RDouble4D *realVarFilter, *realTauij;        
        
    public:
        void InitLESData(SpecDiffHybGrid *);
        void Smagorinsky(SpecDiffHybGrid *, RDouble1D *, RDouble4D *);
        void DynamicSmagorinsky(SpecDiffHybGrid *, CompactDifferenceFirstDerivative *, RDouble1D *, RDouble3D *, RDouble4D *, Complex3D *, Complex3D *, Complex3D *);

	private:
        void AllocLESData(SpecDiffHybGrid *);
        void GetWallDistance(SpecDiffHybGrid *);
        void GetDelta(SpecDiffHybGrid *);
        void SetTestFilter(SpecDiffHybGrid *);  
        void SetFilterType(const double, double, double &);
        void FilterXYInSpectralSpace(PHComplex *, const int, SpecDiffHybGrid *);
        void GetDiffTau(SpecDiffHybGrid *, CompactDifferenceFirstDerivative *);

	};

}



