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
//! @file      DSMCCollision_related.h
//! @brief     Collision_related var definition for DSMC.
//! @author    Jia Chuan.
#pragma once

namespace PHSPACE
{

	using namespace std;
	class Collision_related 
	{
	private:

		double VRC[3],VCM[3],VRR,VR,CVR;
		double **CCG;

	public:

		Collision_related();
		~Collision_related();

	public:
		void setVRC(double VRC[])
		{ for (int i=0;i<3,i++;)
		this->VRC[i]=VRC[i];
		}
		double getVRC(double VRC[])
		{ for (int i=0;i<3,i++;)
		return this->VRC[i];
		}
		void setVRR(double VRR){this->VRR=VRR;}; 
		double getVRR(){return this->VRR;};

		void setVR(double VR){this->VR=VR;}; 
		double getVR(){return this->VR;};

		void setCVR(double CVR){this->CVR=CVR;}; 
		double getCVR(){return this->CVR;};

	public:

		void SetCCG(double n, double value);
		double ** GetCCG(){return CCG;};

	};
}








