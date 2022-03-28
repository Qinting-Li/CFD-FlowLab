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
//! @file      DSMCMolecule_related.h
//! @brief     Molecule_related var definition for DSMC.
//! @author    Jia Chuan.
#pragma once

namespace PHSPACE
{

	using namespace std;
	class Molecule_related 
	{
	private:

		int *IPL,*IPS,*IPCP;
		int **IPVIB;
		double *Erot,*Evib;
		double **PP, **PV;
    
	public:

		Molecule_related();
		~Molecule_related();

	public:
		
		void setIPL(int *IPL){this->IPL=IPL;}; 
		int *getIPL(){return this->IPL;};

		void setIPS(int *IPS){this->IPS=IPS;}; 
		int *getIPS(){return this->IPS;};

		void setIPCP(int *IPCP){this->IPCP=IPCP;}; 
		int *getIPCP(){return this->IPCP;};

	public:

		void SetIPVIB(int n, int value);
		int ** GetIPVIB(){return IPVIB;};

	public:

		void setErot(double *Erot){this->Erot=Erot;}; 
		double *getErot(){return this->Erot;};

	public:

		void SetPP(double n, double value);
		double ** GetPP(){return PP;};

		void SetPV(double n, double value);
		double ** GetPV(){return PV;};

	};

}