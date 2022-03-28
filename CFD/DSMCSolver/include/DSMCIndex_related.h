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
//! @file      DSMCIndex_related.h
//! @brief     Index_related var definition for DSMC.
//! @author    Jia Chuan.
#pragma once

namespace PHSPACE
{

	using namespace std;
	class Index_related 
	{
	private:

		int *IR;
		int *nTmpIdx;
		int **IC,**ISCG;
	public:
        Index_related ();
		~Index_related ();

	public:
		void setIR(int *IR){this->IR=IR;}; 
		int *getIR(){return this->IR;};

    public:
		void setnTmpIdx(int *nTmpIdx){this->nTmpIdx=nTmpIdx;}; 
		int *getnTmpIdx(){return this->nTmpIdx;};

	public:
		void SetIC(int n, int value);
		int ** GetIC(){return IC;};

		void SetISCG(int n, int value);
		int ** GetISCG(){return ISCG;};

	};
}

