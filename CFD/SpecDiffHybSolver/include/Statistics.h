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
//! @file      Statistics.h
//! @brief     statistics that can be computed from fourier planes.
//! @author    ZhangZipei.

#pragma once
#include "TypeDefine.h"
#include "SpecDiffHybGrid.h"

using namespace std;

namespace PHSPACE
{
	class Statistics 
	{

	public:
		Statistics();
		~Statistics();

	public:
        int inutavg, iz_half, naddnum;
        double tauWAverage;
        RDouble2D *statisticsAll, *statisticsAllLocal;
        RDouble2D *epsXYZ, *epsXYZLocal;
        RDouble3D *energyOfUUKx, *energyOfUUKy, *energyOfUUKxAll, *energyOfUUKyAll, *energyOfUUKxLocal, *energyOfUUKyLocal;

        // !<U>, <V>, <W>
        RDouble1D *uAverage, *vAverage, *wAverage;
        RDouble1D *uAverageTime, *vAverageTime, *wAverageTime;

        // !<uu>, <vv>, <ww>, <uv>, <uw>, <vw>
        RDouble1D *uSquare, *vSquare, *wSquare, *uMultiplyV, *uMultiplyW, *vMultiplyW;
        RDouble1D *uSquareTime, *vSquareTime, *wSquareTime, *uMultiplyVTime, *uMultiplyWTime, *vMultiplyWTime;

        // !dissipation rate
        RDouble1D *epsilon, *epsilonTime;
        RDouble2D *epsilonXYZ;
        // !turbulent viscosity
        RDouble1D *nutAverage;

        bool IsInitStatScratch;
    private: 
        string zGridPlot3dFile, zHalfGridPlot3dFile, statisticsNameFile, statisticsInformationFile, statisticsRestartFile, statisticsRestartTimeFile;
        int istats;       

    public:
        void InitStatisticsData(SpecDiffHybGrid *);
        void GetStatisticsLocal(SpecDiffHybGrid *, Complex3D *, Complex3D *, Complex3D *, Complex4D *);
        void OutputStatistics(SpecDiffHybGrid *);

    private:
        void InitStatisticsFileName();
        void OutputPlot3dGrid(SpecDiffHybGrid *);
        void AllocStatisticsdata(SpecDiffHybGrid *);
        void InitStatisticsDataValue(SpecDiffHybGrid *);
        void InputStatistics(SpecDiffHybGrid *);
        void WriteStatisticsName();
        void ReadStatisticsInformation();
        void ReadStatisticsRestart(SpecDiffHybGrid *);
        void ReadStatisticsSelected(SpecDiffHybGrid *);
        void WriteStatisticsInformation();
        void WriteStatisticsRestart(SpecDiffHybGrid *);
        void WriteStatisticsSelected(SpecDiffHybGrid *);
	};
}



