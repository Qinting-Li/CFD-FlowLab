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
//! @file      SpecDiffHybSolver.h
//! @brief     Solver of Spectral-Difference hybrid method.
//! @author    Zhang Zipei.

#pragma once

#include "TypeDefine.h"
#include "Complex.h"
using namespace std;

namespace PHSPACE
{
	class SpecDiffHybSolver 
	{
	public:
		SpecDiffHybSolver();
		~SpecDiffHybSolver();

	private:
        int istep;
		RDouble1D *nonlinearTermCoef, *pressureGradientTermCoef;
        Complex4D *convectiveTerm, *complexLaplacianVelocityTerm;
        Complex5D *complexPressureGradientTerm;

		Complex4D *complexViscousTerm, *complexNonlinearTerm, *complexVelocityGradient, *complexVelocitySquare;
        RDouble4D *realVelocityGradient, *realVelocitySquare;

        RDouble1D *realnut, *realnut_00;
        RDouble3D *realnut3d;
		Complex3D *complexU, *complexV, *complexW, *complexP, *complexDivergenceU;
		      
        string nonlinearTermFile, nonlinearTermBackupFile;
        string pressureGradientTermFile, pressureGradientTermBackupFile;
        string laplacianVelocityTermFile, laplacianVelocityTermBackupFile;
        string timeInformationRestartFile;
        string velocityRestartFile, velocityBackupFile;
        string pressureRestartFile, pressureBackupFile;
        string velocityNameFile, fieldNameFile;
        string screenOutFile;

        bool IsConvectiveTermFailed, IsPressureGradientFailed, IsLaplacianVelocityFailed, IsInitialization;

	public:
		void Run();
//		void GetFirstDerivative(RDouble1D *dvar);

	private:
		void Init();
		void InitTimeDependentData();
		void InitTimeDependentFileName();
		void AllocTimeDependentData();
        void InitGlobal();
        void InputOldTimeData();
        void ReadConvectiveTerm();
        void ReadTimeDependentInformation();
		void GetRHS();
		void Solve();
        void RescaleVelocityField(const int);
		void GetVelocityGradient();
		void GetNonlinearTerm();
		void GetViscousTerm();
        void InitFieldData();
		void AllocFieldData();
        void GetInitField();
        void InitFieldFileName();
        void WriteVelocityName();
        void InitProfile();
        void InitFromOldField();
        void OutputFieldRestartFile(const int);
        void GetAssumedSpectrum( const double, const double, double & );
        void GetAssumedSpectrumNearWall( const double, const double, double & );
        void ReadinVelocityRestartFile();
        void ReadinPressureRestartFile();
        void WriteVelocityRestartFile(const string);
        void WritePressureRestartFile(const string);
        void SetInitTimeDataCoef();
		void SetTimeDataCoef();
        void VanishingDIVproj();
        void OutputFieldName();
        void GetStatistics();
        void OutputResidual(const int);
        void OutputRestartData(const int);
        void OutputOldTimeData(const int);
        void WriteConvectiveTerm(const string);
        void WriteTimeDependentInformation();
        void mode00fix();

	private:
		class SpecDiffHybGrid * GridData;
		class ExplicitDifferenceBoundary * DiffBoundaryData;
		class CompactDifferenceFirstDerivative * CompactDifferenceFirstDerivativeData;
		class CompactDifferenceSecondDerivative * CompactDifferenceSecondDerivativeData;
		class PoissonSolver * PoissonSolverData;
		class CorrectWallDivergence * CorrectWallDivergenceData;
        class RescaleField * RescaleFieldData;
        class Statistics *StatisticsData;
        class LES *LESData;
	};

}