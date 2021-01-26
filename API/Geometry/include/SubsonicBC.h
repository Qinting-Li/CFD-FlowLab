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
//! @file      SubsonicBC.h
//! @brief     To set data structure for subsonic boundary conditions.
//! @author    Chen Liangzhong.

#pragma once
#include "Math_BasisFunction.h"
#include "Geo_Grid.h"
#include "Geo_UnstructGrid.h"
#include "BC_Unstruct.h"
#include "BC_Struct.h"
#include "Precision.h"

using namespace std;

namespace PHSPACE
{
	class SubnicBC
	{
	private:
        int bcID,bcType;
        RDouble pressureTotal;
        RDouble temperatureTotal;
        RDouble massRate;
		RDouble flowDirection[3];
        //vector<RDouble> flowDirection;
	public:
		SubnicBC()
		{
            //flowDirection = vector<RDouble>(3);
            bcID             = 0;
            bcType           = 0;
			flowDirection[0] = 1.0;
			flowDirection[1] = 0.0;
			flowDirection[2] = 0.0;
			pressureTotal    = 0.0;
			temperatureTotal = 0.0;
			massRate         = 0.0;
		}
		~SubnicBC()
		{			
		}

	public:
		int GetBCID()
		{
			return bcID;
		}
		void SetBCID(int id)
		{
			bcID = id;
		}
		int GetBCType()
		{
			return bcType;
		}
		void SetBCType(int type)
		{
			bcType = type;
		}
		RDouble GetTotalPressure()
		{
			return pressureTotal;
		}
		void SetTotalPressure(RDouble p0)
		{
			pressureTotal = p0;
		}
		RDouble GetTotalTemperature()
		{
			return temperatureTotal;
		}
		void SetTotalTemperature(RDouble t0)
		{
			temperatureTotal = t0;
		}
		RDouble GetMassRate()
		{
			return massRate;
		}
		void SetMassRate(RDouble mr)
		{
			massRate = mr;
		}
		/*vector<RDouble> GetFlowDirection()
		{
			return flowDirection;
		}*/
		RDouble * GetFlowDirection()
		{
			return flowDirection;
		}
		/*void SetFlowDirection(vector<RDouble> fd)
		{
			for ( int i = 0; i < 3; ++i )
			{
				flowDirection[i] = fd[i];
			}
		}*/
		void SetFlowDirection( RDouble * fd )
		{
			for ( int i = 0; i < 3; ++i )
			{
				flowDirection[i] = fd[i];
			}
		}
	};
	class SubnicBCInfo
	{
    private:
        string fileName;
        int bcIDMax,bcIDMin;
        int nbc,hasSimilarBC;
		vector<SubnicBC> BCs;

	public:
		SubnicBCInfo()
		{
			bcIDMax = -10000;
			bcIDMin = +10000;
			nbc     = 0;
			hasSimilarBC = 0;
		}
		~SubnicBCInfo()
		{
		}

	public:
		void CalBCIDRange()
		{
			for ( int i = 0; i < nbc; ++ i )
			{
				bcIDMax = MAX( bcIDMax, BCs[i].GetBCID() );
				bcIDMin = MIN( bcIDMin, BCs[i].GetBCID() );
			}
		}
		void SetSimilarBCExistLabel()
		{
			hasSimilarBC = 1;
		}
		int GetBCIDRangLower()
		{
			return bcIDMin;
		}
		int GetBCIDRangUpper()
		{
			return bcIDMax;
		}
		int GetZoneNumberOfSameBCType()
		{
			return nbc;
		}
		void SetBCDataFileName(string name)
		{
			fileName = name;
		}
		string GetBCDataFileName()
		{
			return fileName;
		}
		vector<SubnicBC> GetBCInfo()
		{
			return BCs;
		}
		void CheckSpecificBCExist(int specificBCID);
		void Nondimensional();
		void ReSetReferenceData();
		void ReadBCData();
		void BcastBCData();
	};

	void GetSubsonicBCData();
	bool SetSubsonicBC(Grid *grid_in, BCRegion *bcr, int bctype);
	bool SetSubsonicBC(UnstructGrid * grid_in, BCRecord ** bcr, int bctype, 
		                     RDouble * primSource, RDouble * primTarget,
							 RDouble nxs, RDouble nys, RDouble nzs,
							 RDouble t );
	bool SetSubsonicBCTurb( int bctype, RDouble **q_turb, RDouble * turboo, int n_turb, int le, int re );

	void InFlowBC3D(Grid *grid_in, BCRegion *bcr, SubnicBC & outLet);
	void OutFlowBC3D(Grid *grid_in, BCRegion *bcr, SubnicBC & outLet);

	void InFlowBC3D(UnstructGrid * grid_in, BCRecord ** bcr, SubnicBC & inLet, 
		                  RDouble * primSource, RDouble * primTarget, 
						  RDouble nxs, RDouble nys, RDouble nzs, RDouble t);
	void OutFlowBC3D(UnstructGrid * grid_in, BCRecord ** bcr, SubnicBC & outLet, 
						  RDouble * primSource, RDouble * primTarget, 
						  RDouble nxs, RDouble nys, RDouble nzs, RDouble t);
	void InFlowBC3DTurb( int bctype, RDouble **q_turb, RDouble * turboo, int n_turb, int le, int re );
	void OutFlowBC3DTurb( int bctype, RDouble **q_turb, RDouble * turboo, int n_turb, int le, int re );
}
