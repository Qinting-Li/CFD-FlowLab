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
//! @file      DSMCBuffer_related.h
//! @brief     Buffer_related var definition for DSMC.
//! @author    Jia Chuan.

#pragma once

namespace PHSPACE
{

	using namespace std;
	class Buffer_related 
	{
	private:

		int MNB,NMB;
		int *IPS_B,*IPL_B,*IR_B;
		int **IPVIB_B;
		double *PRB;
		double **PPB,**PVB;

	public:
		Buffer_related();
		~Buffer_related();


	public:
		void setMNB(int MNB){this->MNB=MNB;}; 
		int getMNB(){return this->MNB;};

		void setNMB(int NMB){this->NMB=NMB;}; 
		int getNMB(){return this->NMB;};

	public:
		void setIPS_B(int *IPS_B){this->IPS_B=IPS_B;}; 
		int *getIPS_B(){return this->IPS_B;};

		void setIPL_B(int *IPL_B){this->IPL_B=IPL_B;}; 
		int *getIPL_B(){return this->IPL_B;};

		void setIR_B(int *IR_B){this->IR_B=IR_B;}; 
		int *getIR_B(){return this->IR_B;};

	public:
		void SetIPVIB_B(int n, int value);
		int ** GetIPVIB_B(){return IPVIB_B;};

	public:
		void setPRB(double *PRB){this->PRB=PRB;}; 
		double *getPRB(){return this->PRB;};

	public:
		void SetPPB(double n, double value);
		double ** GetPPB(){return PPB;};

		void SetPVB(double n, double value);
		double ** GetPVB(){return PVB;};
	};
}