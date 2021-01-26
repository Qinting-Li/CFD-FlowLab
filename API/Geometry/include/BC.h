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
//! @file      BC.h
//! @brief     It defines the general Boundary Condition definition.
//! @author    Bell, He Xin.

#pragma once
#include "TypeDefine.h"
#include "Data_Param.h"
using namespace std;
namespace PHSPACE
{
namespace PHENGLEI
{
	const int PERIODIC              =   99;
	const int WAKE                  = - 2 ;
	const int INTERFACE             = - 1 ;
	const int NO_BOUNDARY_CONDITION =   0 ;
	const int EXTRAPOLATION         =   1 ;
	const int SOLID_SURFACE         =   2 ;
	const int SYMMETRY              =   3 ;
	const int FARFIELD              =   4 ;
	const int INFLOW                =   5 ;
	const int OUTFLOW               =   6 ;
	const int PRESSURE_INLET		=   52;
	const int PRESSURE_OUTLET		=   62;
	const int OUTFLOW_CONFINED      =   61;
	const int POLE                  =   7 ;
	const int POLE1                 =   71;
	const int POLE2                 =   72;
	const int POLE3                 =   73;
	const int GENERIC_1             =   8 ;
	const int GENERIC_2             =   9 ;
	const int GENERIC_3             =   10;
	const int INTERIOR              =   100;
	//! Structured patched grid by Wang Ziwei.
	//! If the boundary condition number of a boundary face between tSTART_PATCH_BC and END_PATCH_BC, the boundary face is considered as patching face.
	//! If the boundary condition numbers of two patching faces are same, there may be section patching between the two faces.
	//! Otherwise,if the boundary condition numbers of two patching faces are different,there is no patching relationship between the two faces.
	const int START_PATCH_BC             = 81 ;
	const int END_PATCH_BC             = 90 ;

	const int EXTERNAL_BC           =   9;    //! External overlap cells label.

	//! Boundary condition for special case.
	const int PRESSURE_INLET_PLATE	=   1001;
	const int PRESSURE_OUTLET_PLATE	=   1002;

	const int StartNumberOfBodyOfHyperFLOW = 11;
	const int ENDNumberOfBodyOfHyperFLOW   = 30;
	const int BODY1              =   11;
	const int BODY2              =   12;
	const int BODY3              =   13;
	const int BODY4              =   14;
	const int BODY5              =   15;
	const int BODY6              =   16;
	const int BODY7              =   17;
	const int BODY8              =   18;
	const int BODY9              =   19;
	const int BODY10             =   20;
	const int BODY11             =   21;
	const int BODY12             =   22;
	const int BODY13             =   23;
	const int BODY14             =   24;
	const int BODY15             =   25;
	const int BODY16             =   26;
	const int BODY17             =   27;
	const int BODY18             =   28;
	const int BODY19             =   29;
	const int BODY20             =   30;
	const int OVERSET            =   1000;
}

namespace GRIDGEN_SPACE
{
	const int INTERFACE             = - 1 ;
	const int NO_BOUNDARY_CONDITION =   0 ;
	const int SOLID_SURFACE         =   2 ;
	const int SYMMETRY              =   3 ;
	const int FARFIELD              =   4 ;
	const int INFLOW                =   5 ;
	const int OUTFLOW               =   6 ;
	const int POLE                  =   7 ;
	const int GENERIC_1             =   8 ;
	const int GENERIC_2             =   9 ;
	const int GENERIC_3             =   10;
	const int PERIODIC              =   99;   
}

namespace FLUENT_SPACE
{
    const int FLUENT_COMMENT    = 0;
    const int FLUENT_HEADER     = 1;
    const int FLUENT_DIMENSIONS = 2;
    const int FLUENT_NODES      = 10;
    const int FLUENT_CELLS      = 12;
    const int FLUENT_FACES      = 13;
    const int FLUENT_C39        = 39;
    const int FLUENT_C45        = 45;
    const int MaxFluentBCType    = 38;

    const int UNSPECIFIED        = 0;
    const int INTERIOR           = 2;
    const int WALL               = 3;
    const int PRESSURE_INLET     = 4;
    const int INLET_VENT         = 4;
    const int INTAKE_FAN         = 4;
    const int PRESSURE_OUTLET    = 5;
    const int EXHAUST_FAN        = 5;
    const int OUTLET_VENT        = 5;
    const int SYMMETRY           = 7;
    const int PERIODIC_SHADOW    = 8;
    const int PRESSURE_FAR_FIELD = 9;
    const int VELOCITY_INLET     = 10;
    const int PERIODIC           = 12;
    const int FAN                = 14;
    const int POROUS_JUMP        = 14;
    const int RADIATOR           = 14;
    const int MASS_FLOW_INLET    = 20;
    const int INTERFACE          = 24;
    const int PARENT             = 31;
    const int OUTFLOW            = 36;
    const int AXIS               = 37;
}

class BoundaryCondition
{
public:
	BoundaryCondition()
	{
		bcType  = 0;
		bcName  = "";
		bcParamDataBase = 0;
	}
	~BoundaryCondition()
	{
		delete bcParamDataBase;
	}

private:
	int bcType;
	string bcName;
	Data_Param * bcParamDataBase;

public:
	//! Get boundary conditions type.
	int GetBCType() const {return this->bcType;}
	Data_Param * GetBCParamDataBase() const { return this->bcParamDataBase; }
	const string & GetBCName() const {return this->bcName;}

	void GetData(const string &name,void *data, int type, int size)
	{
		bcParamDataBase->GetData(name, data, type, size);
	}

    int CheckData(const string &name)
    {
        return bcParamDataBase->CheckDataExist(name);
    }

	void SetBCType(int bcType) {this->bcType = bcType;}
	void SetBCName(const string & bcName) {this->bcName = bcName;}
	void SetBCParamDataBase(Data_Param * paraDB) {this->bcParamDataBase = paraDB; }
};

//! @brief Global boundary conditions, all of the global boundary conditions are stored here.
class GlobalBoundaryCondition
{
private:	
	//! Warning: it has not been deleted!.
	static vector<BoundaryCondition *> * globalBCList;

private:
	GlobalBoundaryCondition();
	~GlobalBoundaryCondition();

public:
	//! Get number of total BC.
	static uint_t GetNumberofBC()
	{
		if(!globalBCList)
			return 0;
		else
			return globalBCList->size(); 
	}

	//! Get number of BC with solid_wall BC type.
	static int GetNumberofSolidWallPart();

	//! Is solid wall.
	static bool IsSolidWallBC(const int & partID) 
	{return ( (*globalBCList)[partID]->GetBCType() == PHENGLEI::SOLID_SURFACE );}

	static BoundaryCondition * GetBC(const int & iBC) { return (*globalBCList)[iBC]; }
	
	//! Get global boundary conditions list.
	static vector<BoundaryCondition *> * GetGlobalBoundaryConditionList() {return globalBCList;}
    
	//! Read global boundary conditions.
	static void ReadGlobalBoundaryCondition();

	static void SetGlobalDataBaseByGlobalBC();

	static void ChangeBCTypeByGlobalBC();

private:
	static void ParseBCFromFile(fstream & file);
};

//! Judge whether a boundary condition type number belongs to patching boundary condition.
inline bool ISPatchInterfaceBC(const int& bctype)
{
	if ( bctype >= PHENGLEI::START_PATCH_BC && bctype <= PHENGLEI::END_PATCH_BC ) 
		return true;
	else 
		return false;
}
//! Judge whether a boundary condition type number belongs to interface.
bool IsInterface( const int& boundaryConditionType );

//! Judge whether a boundary condition type number belongs to solid surface.
inline bool IsWall(const int& bctype)
{
	return ( bctype == PHENGLEI::SOLID_SURFACE );
}

map<int,string> CreatFantasybctype();

}
