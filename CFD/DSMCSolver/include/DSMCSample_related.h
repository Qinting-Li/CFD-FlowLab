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
//! @file      DSMCSample_related.h
//! @brief     Sample_related var definition for DSMC.
//! @author    Jia Chuan.
#pragma once

namespace PHSPACE
{

	using namespace std;
	class Sample_related 
	{
	private:

		double **CSR,**SVIB,**VAR,**nVAR,**VARS,**nVARS; 
		double ***CS,***CSV,***CSVS,***VARSP,***nVARSP,***DF;
		double ****CSS;
		double **CSDV,**CSDE,**CSDVV;



	public:

		Sample_related();
		~Sample_related();

	public:

		void SetCSR(double n, double value);
		double ** GetCSR(){return CSR;};

		void SetSVIB(double n, double value);
		double ** GetSVIB(){return SVIB;};

		void SetVAR(double n, double value);
		double ** GetVAR(){return VAR;};

		void SetnVAR(double n, double value);
		double ** GetnVAR(){return nVAR;};

		void SetVARS(double n, double value);
		double ** GetVARS(){return VARS;};

		void SetnVARS(double n, double value);
		double ** GetnVARS(){return nVARS;};

	public:

		void SetCS(double n, double value,double value_1);
		double *** GetCS(){return CS;};

		void SetCSV(double n, double value,double value_1);
		double *** GetCSV(){return CSV;};

		void SetCSVS(double n, double value,double value_1);
		double *** GetCSVS(){return CSVS;};

		void SetVARSP(double n, double value,double value_1);
		double *** GetVARSP(){return VARSP;};

		void SetnVARSP(double n, double value,double value_1);
		double *** GetnVARSP(){return nVARSP;};

		void SetDF(double n, double value,double value_1);
		double *** GetDF(){return DF;};

	public:

		void SetCSS(double n, double value,double value_1,double value_2);
		double **** GetCSS(){return CSS;};

	public:

		void SetCSDV(double n, double value);
		double ** GetCSDV(){return CSDV;};

		void SetCSDE(double n, double value);
		double ** GetCSDE(){return CSDE;};

		void SetCSDVV(double n, double value);
		double ** GetCSDVV(){return CSDVV;};

	};
}


   