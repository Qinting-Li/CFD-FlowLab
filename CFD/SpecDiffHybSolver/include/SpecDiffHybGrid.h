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
//! @file      SpecDiffHybGrid.h
//! @brief     Grid data for SpecDiffHybSolver.
//! @author    ZhangZipei.

#pragma once
#include "TypeDefine.h"
#include "Complex.h"
//#include "CompactDifferenceFirstDerivative.h"
using namespace std;

namespace PHSPACE
{

	class SpecDiffHybGrid 
	{
	public:
		SpecDiffHybGrid();
		~SpecDiffHybGrid();
	private:   
		int nyquistX, nyquistY;		
		Int1D *physicalStartIndex, *physicalEndIndex, *physicalIndexSize;
//		Int1D *fourierStartIndexGlobal, *fourierEndIndexGlobal;
        Int1D *dimensionOfP3DFFT;		
        string gridFileName, zgridFileName;

	public:
		int nX, nY, nZ, nX2, nY2, nZ2, nXYAll, nXYZFourier, nXYZPhysical;
		int iStartFourierIndex , jStartFourierIndex , kStartFourierIndex , iEndFourierIndex , jEndFourierIndex , kEndFourierIndex ;
		int iStartPhysicalIndex, jStartPhysicalIndex, kStartPhysicalIndex, iEndPhysicalIndex, jEndPhysicalIndex, kEndPhysicalIndex;
		int x0Fall, y0Fall, z0Fall, xNFall, yNFall, zNFall;
		double realKNyquistX, realKNyquistY;
        int procIDWithKx0Ky0;
        Int1D *fourierStartIndex, *fourierEndIndex, *fourierIndexSize;
		Int2D *startAndEndFourierIndexOfAll, *startAndEndPhysicalIndexOfAll;
//		int **startAndEndAll, **startAndEndPhysical;
//				Int1D * physicalStart, ;
		RDouble1D *realKX, *realKY;
		Complex1D *complexKX, *complexKY;
		RDouble1D *realZ, *realDetaDz;

	private:
        void InitGridFileName();

	public:

		void InitGridData();
		void InitGridSize();
		void AllocGridData();
		void GetGrid();
        void SetZGrid();
		void Getiz_p2();
		void InitWavenumber();
//		void GetDetaDz(CompactDifferenceFirstDerivative *CompactDifferenceFirstDerivativeData);
		void InitP3DFFT();

	};


}



