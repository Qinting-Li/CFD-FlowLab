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
//! @file      RescaleField.h
//! @brief     Rescale velocity every time step if necessary.
//! @author    ZhangZipei.

#pragma once
#include "TypeDefine.h"
#include "SpecDiffHybGrid.h"
using namespace std;

namespace PHSPACE
{

	class RescaleField 
	{

	public:
		RescaleField();
		~RescaleField();

	public:

	private: 
        double uavgDivideByUtau;
        RDouble1D *uAverage, *uFluctuation, *vFluctuation, *wFluctuation;
        RDouble1D *uAverageOfReference, *uFluctuationOfReference, *vFluctuationOfReference, *wFluctuationOfReference;
        RDouble2D *rescaleStatistic, *rescaleStatisticLocal;

    public:
        void InitRescaleData( SpecDiffHybGrid * );
        void GetUVWMagnitude(SpecDiffHybGrid *, Complex3D *, Complex3D *, Complex3D *);
        void RescaleUAverage(SpecDiffHybGrid *, Complex3D *);
        void RescaleFluctuation(SpecDiffHybGrid *, Complex3D *, Complex3D *, Complex3D *);
        void FixUAverage(SpecDiffHybGrid *, Complex3D *);

	private:
        void AllocRescaleData(SpecDiffHybGrid *);
		void GetUavg();
        void GetUVWMagnitudeOfReference(SpecDiffHybGrid *);
        void InputUVWMagnitudeOfReference( SpecDiffHybGrid * );
        void GetUVWLocal(SpecDiffHybGrid *, Complex3D *, Complex3D *, Complex3D *);
        void Interpolation( RDouble *, RDouble *, int, RDouble *, RDouble *, int );

	};

}



