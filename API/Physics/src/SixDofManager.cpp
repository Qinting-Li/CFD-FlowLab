#include "PHMpi.h"
#include "AleManager.h"
#include "Force.h"
#include "SixDofManager.h"
#include "Pointer.h"
#include "GlobalDataBase.h"
#include "Constants.h"
#include "Force.h"
#include "PHIO.h"
#include "IO_FileName.h"
#include "IO_FileReader.h"
#include "TK_Exit.h"
#include "TK_Log.h"
#include "Math_BasisFunction.h"

namespace PHSPACE
{

RDouble ** MassCharacter::GetMomentOfInertiaTensor()
{
    RDouble ** momentOfInertiaTensor = PHSPACE::NewPointer2< RDouble > ( 3, 3 );

    momentOfInertiaTensor[ 0 ][ 0 ] =   momentOfInertiaTensorXX;
    momentOfInertiaTensor[ 0 ][ 1 ] = - momentOfInertiaTensorXY;
    momentOfInertiaTensor[ 0 ][ 2 ] = - momentOfInertiaTensorXZ;

    momentOfInertiaTensor[ 1 ][ 0 ] = - momentOfInertiaTensorXY;
    momentOfInertiaTensor[ 1 ][ 1 ] =   momentOfInertiaTensorYY;
    momentOfInertiaTensor[ 1 ][ 2 ] = - momentOfInertiaTensorYZ;

    momentOfInertiaTensor[ 2 ][ 0 ] = - momentOfInertiaTensorXZ;
    momentOfInertiaTensor[ 2 ][ 1 ] = - momentOfInertiaTensorYZ;
    momentOfInertiaTensor[ 2 ][ 2 ] =   momentOfInertiaTensorZZ;

    return momentOfInertiaTensor;
}

RDouble ** MassCharacter::GetInverseMomentOfInertiaTensor()
{
    RDouble ** inverseMomentOfInertiaTensor = PHSPACE::NewPointer2< RDouble > ( 3, 3 );

    RDouble d   =       momentOfInertiaTensorXX * momentOfInertiaTensorYY * momentOfInertiaTensorZZ -
                        momentOfInertiaTensorXX * momentOfInertiaTensorYZ * momentOfInertiaTensorYZ -
                        momentOfInertiaTensorYY * momentOfInertiaTensorXZ * momentOfInertiaTensorXZ -
                        momentOfInertiaTensorZZ * momentOfInertiaTensorXY * momentOfInertiaTensorXY -
                  2.0 * momentOfInertiaTensorXY * momentOfInertiaTensorXZ * momentOfInertiaTensorYZ;

    RDouble a11 = momentOfInertiaTensorYY * momentOfInertiaTensorZZ - momentOfInertiaTensorYZ * momentOfInertiaTensorYZ;
    RDouble a12 = momentOfInertiaTensorXY * momentOfInertiaTensorZZ + momentOfInertiaTensorYZ * momentOfInertiaTensorXZ;
    RDouble a13 = momentOfInertiaTensorXY * momentOfInertiaTensorYZ + momentOfInertiaTensorXZ * momentOfInertiaTensorYY;

    RDouble a21 = momentOfInertiaTensorXY * momentOfInertiaTensorZZ + momentOfInertiaTensorYZ * momentOfInertiaTensorXZ;
    RDouble a22 = momentOfInertiaTensorXX * momentOfInertiaTensorZZ - momentOfInertiaTensorXZ * momentOfInertiaTensorXZ;
    RDouble a23 = momentOfInertiaTensorXX * momentOfInertiaTensorYZ + momentOfInertiaTensorXZ * momentOfInertiaTensorXY;

    RDouble a31 = momentOfInertiaTensorXY * momentOfInertiaTensorYZ + momentOfInertiaTensorXZ * momentOfInertiaTensorYY;
    RDouble a32 = momentOfInertiaTensorXX * momentOfInertiaTensorYZ + momentOfInertiaTensorXZ * momentOfInertiaTensorXY;
    RDouble a33 = momentOfInertiaTensorXX * momentOfInertiaTensorYY - momentOfInertiaTensorXY * momentOfInertiaTensorXY;

    inverseMomentOfInertiaTensor[ 0 ][ 0 ] = a11 / d;
    inverseMomentOfInertiaTensor[ 0 ][ 1 ] = a12 / d;
    inverseMomentOfInertiaTensor[ 0 ][ 2 ] = a13 / d;

    inverseMomentOfInertiaTensor[ 1 ][ 0 ] = a21 / d;
    inverseMomentOfInertiaTensor[ 1 ][ 1 ] = a22 / d;
    inverseMomentOfInertiaTensor[ 1 ][ 2 ] = a23 / d;

    inverseMomentOfInertiaTensor[ 2 ][ 0 ] = a31 / d;
    inverseMomentOfInertiaTensor[ 2 ][ 1 ] = a32 / d;
    inverseMomentOfInertiaTensor[ 2 ][ 2 ] = a33 / d;

    return inverseMomentOfInertiaTensor;
}

MassCharacter::MassCharacter( RDouble mass, RDouble momentOfInertiaTensorXX, RDouble momentOfInertiaTensorYY, RDouble momentOfInertiaTensorZZ, RDouble momentOfInertiaTensorXY, RDouble momentOfInertiaTensorXZ, RDouble momentOfInertiaTensorYZ )
{
    this->mass = mass;
    this->momentOfInertiaTensorXX  = momentOfInertiaTensorXX;
    this->momentOfInertiaTensorYY  = momentOfInertiaTensorYY;
    this->momentOfInertiaTensorZZ  = momentOfInertiaTensorZZ;
    this->momentOfInertiaTensorXY  = momentOfInertiaTensorXY;
    this->momentOfInertiaTensorXZ  = momentOfInertiaTensorXZ;
    this->momentOfInertiaTensorYZ  = momentOfInertiaTensorYZ;
}

void MassCharacter::SetMassCharacter( RDouble mass, vector< RDouble > & massMatrix )
{
    this->mass = mass;
    this->momentOfInertiaTensorXX  = massMatrix[ 0 ];
    this->momentOfInertiaTensorYY  = massMatrix[ 1 ];
    this->momentOfInertiaTensorZZ  = massMatrix[ 2 ];
    this->momentOfInertiaTensorXY  = massMatrix[ 3 ];
    this->momentOfInertiaTensorXZ  = massMatrix[ 4 ];
    this->momentOfInertiaTensorYZ  = massMatrix[ 5 ];
}

void MassCharacter::SetMassCharacter( RDouble mass, RDouble momentOfInertiaTensorXX, RDouble momentOfInertiaTensorYY, RDouble momentOfInertiaTensorZZ, RDouble momentOfInertiaTensorXY, RDouble momentOfInertiaTensorXZ, RDouble momentOfInertiaTensorYZ )
{
    this->mass = mass;
    this->momentOfInertiaTensorXX  = momentOfInertiaTensorXX;
    this->momentOfInertiaTensorYY  = momentOfInertiaTensorYY;
    this->momentOfInertiaTensorZZ  = momentOfInertiaTensorZZ;
    this->momentOfInertiaTensorXY  = momentOfInertiaTensorXY;
    this->momentOfInertiaTensorXZ  = momentOfInertiaTensorXZ;
    this->momentOfInertiaTensorYZ  = momentOfInertiaTensorYZ;
}

void MassCharacter::GetMassCharacter( RDouble & mass, RDouble & momentOfInertiaTensorXX, RDouble & momentOfInertiaTensorYY, RDouble & momentOfInertiaTensorZZ,RDouble & momentOfInertiaTensorXY,RDouble & momentOfInertiaTensorXZ,RDouble & momentOfInertiaTensorYZ )
{
    mass = this->mass;
    momentOfInertiaTensorXX  = this->momentOfInertiaTensorXX;
    momentOfInertiaTensorYY  = this->momentOfInertiaTensorYY;
    momentOfInertiaTensorZZ  = this->momentOfInertiaTensorZZ;
    momentOfInertiaTensorXY  = this->momentOfInertiaTensorXY;
    momentOfInertiaTensorXZ  = this->momentOfInertiaTensorXZ;
    momentOfInertiaTensorYZ  = this->momentOfInertiaTensorYZ;
}

void MassCharacter::GetMomentInertia( vector< RDouble > & momentOfInertia )
{
    momentOfInertia[ 0 ]  = this->momentOfInertiaTensorXX;
    momentOfInertia[ 1 ]  = this->momentOfInertiaTensorYY;
    momentOfInertia[ 2 ]  = this->momentOfInertiaTensorZZ;
    momentOfInertia[ 3 ]  = this->momentOfInertiaTensorXY;
    momentOfInertia[ 4 ]  = this->momentOfInertiaTensorXZ;
    momentOfInertia[ 5 ]  = this->momentOfInertiaTensorYZ;
}

void MassCharacter::GetMomentInertia( RDouble & momentOfInertiaTensorXX, RDouble & momentOfInertiaTensorYY, RDouble & momentOfInertiaTensorZZ, RDouble & momentOfInertiaTensorXY, RDouble & momentOfInertiaTensorXZ, RDouble & momentOfInertiaTensorYZ )
{
    momentOfInertiaTensorXX  = this->momentOfInertiaTensorXX;
    momentOfInertiaTensorYY  = this->momentOfInertiaTensorYY;
    momentOfInertiaTensorZZ  = this->momentOfInertiaTensorZZ;
    momentOfInertiaTensorXY  = this->momentOfInertiaTensorXY;
    momentOfInertiaTensorXZ  = this->momentOfInertiaTensorXZ;
    momentOfInertiaTensorYZ  = this->momentOfInertiaTensorYZ;
}

SixDofParameter::SixDofParameter()
{
    massCenterX = 0.0;
    massCenterY = 0.0;
    massCenterZ = 0.0;

    eulerAngleAlpha = 0.0;
    eulerAngleBeta  = 0.0;
    eulerAngleGama  = 0.0; 

    massCenterVelocityX    = 0.0;
    massCenterVelocityY    = 0.0;
    massCenterVelocityZ    = 0.0;

    angularVelocityX = 0.0;
    angularVelocityY = 0.0;
    angularVelocityZ = 0.0; 

    referenceLength = GlobalDataBase::GetDoubleParaFromDB("referenceLength");
}

void SixDofParameter::GetNondimensionalMassCenter( RDouble & dimensionlessMassCenterX, RDouble & dimensionlessMassCenterY, RDouble & dimensionlessMassCenterZ )
{
    dimensionlessMassCenterX = this->massCenterX / referenceLength;
    dimensionlessMassCenterY = this->massCenterY / referenceLength;
    dimensionlessMassCenterZ = this->massCenterZ / referenceLength;
}

RDouble SixDofParameter::GetNondimensionalMassCenterX()
{
    return this->massCenterX / referenceLength;
}

RDouble SixDofParameter::GetNondimensionalMassCenterY()
{
    return this->massCenterY / referenceLength;
}

RDouble SixDofParameter::GetNondimensionalMassCenterZ()
{
    return this->massCenterZ / referenceLength;
}

void SixDofParameter::GetNondimensionalMassCenterVector( vector< RDouble > & nondimensionalMassCenterVector )
{
    nondimensionalMassCenterVector[ 0 ] = this->massCenterX / referenceLength;
    nondimensionalMassCenterVector[ 1 ] = this->massCenterY / referenceLength;
    nondimensionalMassCenterVector[ 2 ] = this->massCenterZ / referenceLength;
}

void SixDofParameter::GetGeneralVariableContainer( vector< RDouble > & generalVariableContainer )
{
    this->GetMassCenter        ( generalVariableContainer[ 0 ], generalVariableContainer[ 1  ], generalVariableContainer[ 2  ] );
    this->GetMassCenterVelocity( generalVariableContainer[ 3 ], generalVariableContainer[ 4  ], generalVariableContainer[ 5  ] );
    this->GetEulerAngle        ( generalVariableContainer[ 6 ], generalVariableContainer[ 7  ], generalVariableContainer[ 8  ] );
    this->GetAngularVelocity   ( generalVariableContainer[ 9 ], generalVariableContainer[ 10 ], generalVariableContainer[ 11 ] );
}

void SixDofParameter::SetGeneralVariableContainer( RDouble * generalVariableContainer )
{
    this->SetMassCenter        ( generalVariableContainer[ 0 ], generalVariableContainer[ 1  ], generalVariableContainer[ 2  ] );
    this->SetMassCenterVelocity( generalVariableContainer[ 3 ], generalVariableContainer[ 4  ], generalVariableContainer[ 5  ] );
    this->SetEulerAngle        ( generalVariableContainer[ 6 ], generalVariableContainer[ 7  ], generalVariableContainer[ 8  ] );
    this->SetAngularVelocity   ( generalVariableContainer[ 9 ], generalVariableContainer[ 10 ], generalVariableContainer[ 11 ] );
}

void SixDofParameter::SetGeneralVariableContainer( vector< RDouble > & massCenter, vector< RDouble >& attitudeAngle, 
    vector< RDouble >& massCenterVelocity, vector< RDouble >& angularVelocity)
{
    this->SetMassCenter        (massCenter        [0], massCenter        [1],         massCenter[2]);
    this->SetEulerAngle        (attitudeAngle     [0], attitudeAngle     [1],      attitudeAngle[2]);
    this->SetMassCenterVelocity(massCenterVelocity[0], massCenterVelocity[1], massCenterVelocity[2]);
    this->SetAngularVelocity   (angularVelocity   [0], angularVelocity   [1],    angularVelocity[2]);
}

void SixDofParameter::UpdateEulerAngle( RDouble * eulerAngleVariationContainer )
{
    eulerAngleGama  += eulerAngleVariationContainer[ 0 ];
    eulerAngleBeta  += eulerAngleVariationContainer[ 1 ];
    eulerAngleAlpha += eulerAngleVariationContainer[ 2 ];
}

void SixDofParameter::UpdateAngularVelocity( RDouble * angularVelocityVariationContainer )
{
    angularVelocityX += angularVelocityVariationContainer[ 0 ];
    angularVelocityY += angularVelocityVariationContainer[ 1 ];
    angularVelocityZ += angularVelocityVariationContainer[ 2 ];
}
//
//FYControlLaw::FYControlLaw( SixDofSolverManager * sixDofSolverManager ) : SixDofSolver ( sixDofSolverManager )
//{
//	integralValue = 0.0;
//
//	vector< RDouble > controlLawParameters;
//	GetRDoubleVectorFromDataBase( controlLawParameters, "controlLawParameters", 8 );
//	
//	k1                = controlLawParameters[ 0 ];
//	k2                = controlLawParameters[ 1 ];
//	k3                = controlLawParameters[ 2 ];
//	rudderLowerLimit  = controlLawParameters[ 3 ];
//	rudderUpperLimit  = controlLawParameters[ 4 ];
//	rudderVelocity    = controlLawParameters[ 5 ];
//	integralLoopLimit = controlLawParameters[ 6 ];
//	otherLimit        = controlLawParameters[ 7 ];
//
//	dimensionalAleTimeStep = PHSPACE::GetDimensionalAleTimeStep();
//
//	int numberOfControlStages = PHSPACE::GetIntegerParameterFromDataBase( "numberOfControlStages" );
//	GetRDoubleVectorFromDataBase( controlTargetAngleArray, "controlTargetAngleArray", numberOfControlStages );
//	GetRDoubleVectorFromDataBase( controlTargetTimeArray , "controlTargetTimeArray" , numberOfControlStages );
//}
//
//void FYControlLaw::ModifyTargetStateByStages()
//{
//	RDouble dimensionalAleTime = GetDimensionalAleSimulationTime();
//
//	int iStage = 0;
//	do
//	{
//		targetState = controlTargetAngleArray[ iStage ];
//	}while( dimensionalAleTime + 0.5 * dimensionalAleTimeStep > controlTargetTimeArray[ iStage ++ ] );
//}
//
//FY2DControlLaw::FY2DControlLaw( SixDofSolverManager * sixDofSolverManager ) : SixDofSolver ( sixDofSolverManager )
//{
//	integralValueA1  = 0.0;
//	integralValueA2  = 0.0;
//	integralValueB   = 0.0;
//	oldRudderAngleA1 = 0.0;
//	oldRudderAngleA2 = 0.0;
//	oldRudderAngleB  = 0.0;
//
//	vector< RDouble > controlLawParameters;
//	GetRDoubleVectorFromDataBase( controlLawParameters, "controlLawParameters", 16 );
//	
//	k1A                          = controlLawParameters[  0 ];
//	k2A                          = controlLawParameters[  1 ];
//	k3A                          = controlLawParameters[  2 ];
//	integralLoopLimitA           = controlLawParameters[  3 ];
//	largeAttackAngleLimitA       = controlLawParameters[  4 ];
//
//	k1B                          = controlLawParameters[  5 ];
//	k2B                          = controlLawParameters[  6 ];
//	k3B                          = controlLawParameters[  7 ];
//	deltaTargetStateLimitB       = controlLawParameters[  8 ];
//	integralIncreaseLimitB       = controlLawParameters[  9 ];
//	integralLoopLimitB           = controlLawParameters[ 10 ];
//	angularVelocityCommandLimitB = controlLawParameters[ 11 ];
//
//	rudderVelocityLimit          = controlLawParameters[ 12 ];
//	rudderAngleLimitA            = controlLawParameters[ 13 ];
//	rudderAngleLimitB            = controlLawParameters[ 14 ];
//	rudderAngleLimit             = controlLawParameters[ 15 ];
//
//	dimensionalAleTimeStep = PHSPACE::GetDimensionalAleTimeStep();
//
//	int numberOfControlStagesA = PHSPACE::GetIntegerParameterFromDataBase( "numberOfControlStagesA" );
//	GetRDoubleVectorFromDataBase( controlTargetAngleArrayA, "controlTargetAngleArrayA", numberOfControlStagesA );
//	GetRDoubleVectorFromDataBase( controlTargetTimeArrayA , "controlTargetTimeArrayA" , numberOfControlStagesA );
//
//	int numberOfControlStagesB = PHSPACE::GetIntegerParameterFromDataBase( "numberOfControlStagesB" );
//	GetRDoubleVectorFromDataBase( controlTargetAngleArrayB, "controlTargetAngleArrayB", numberOfControlStagesB );
//	GetRDoubleVectorFromDataBase( controlTargetTimeArrayB , "controlTargetTimeArrayB" , numberOfControlStagesB );
//}
//
//void FY2DControlLaw::ModifyTargetStateByStages()
//{
//	RDouble dimensionalAleTime = GetDimensionalAleSimulationTime();
//
//	int iStageA = 0;
//	do
//	{
//		targetStateA = controlTargetAngleArrayA[ iStageA ];
//	}while( dimensionalAleTime + 0.5 * dimensionalAleTimeStep > controlTargetTimeArrayA[ iStageA ++ ] );
//
//	int iStageB = 0;
//	do
//	{
//		targetStateB = controlTargetAngleArrayB[ iStageB ];
//	}while( dimensionalAleTime + 0.5 * dimensionalAleTimeStep > controlTargetTimeArrayB[ iStageB ++ ] );
//}
//
//void FYDIControlLaw::Initialize( int bodyIndexIn )
//{
//	SixDofSolver::Initialize( bodyIndexIn );
//
//	int fartherIndex = GetIntegerParameterFromDataBase( bodyIndex, "fartherIndex" );
//	MassCharacter * massCharacterFarther = sixDofSolverManager->GetMassCharacter( fartherIndex );
//	momentOfInertiaTensorZZ = massCharacterFarther->GetMomentOfInertiaTensorZZ();
//
//	RDouble referenceAreaOfAerodynamicForce    = PHSPACE::GetReferenceAreaOfAerodynamicForce();
//    RDouble referenceLengthOfAerodynamicMoment = PHSPACE::GetReferenceLengthOfAerodynamicMoment();
//	RDouble referenceDensity                   = PHSPACE::GetReferenceDensity();
//	RDouble referenceVelocity                  = PHSPACE::GetReferenceVelocity();
//	nodimensionalForceToDimensionalForce = half * referenceDensity * PHSPACE::SQR( referenceVelocity ) 
//		                                 * referenceAreaOfAerodynamicForce * referenceLengthOfAerodynamicMoment;
//
//	//本类采用飞行力学坐标系
//	k1                = k1;//k_a
//	k2                = k2;//k_q
//	k3                = k3;//k_Ia
//
//	rudderLowerLimit  = rudderLowerLimit  * DegToRad;
//	rudderUpperLimit  = rudderUpperLimit  * DegToRad;
//	rudderVelocity    = rudderVelocity    * DegToRad * dimensionalAleTimeStep;//从给定的每秒偏转值转换到每步
//
//	integralLoopLimit = integralLoopLimit / k3 / k2 / momentOfInertiaTensorZZ;//积分池大小限制
//	otherLimit        = otherLimit * DegToRad;                                //角速率指令限制
//
//	//气动数据库
//	fstream file;
//	string uDFSixDofFileName = GetStringParameterFromDataBase( bodyIndex, "uDFSixDofFileName" );
//	OpenFile( file, uDFSixDofFileName, ios_base::in );
//
//	file >> numberOfRudder;
//	file >> numberOfAlpha ;
//
//	rudderArray           .resize( numberOfRudder );
//	alphaArray            .resize( numberOfAlpha  );
//	rudderControlCoefArray.resize( numberOfAlpha  );
//	cmzArray              .resize( numberOfAlpha  );
//	for ( int i = 0; i < numberOfAlpha; ++ i )
//	{
//		cmzArray[ i ].resize( numberOfRudder );
//	}
//
//	for ( int i = 0; i < numberOfRudder; ++ i )
//	{
//		file >> rudderArray[ i ];
//	}
//	for ( int i = 0; i < numberOfAlpha; ++ i )
//	{
//		file >> alphaArray[ i ];
//	}
//	for ( int i = 0; i < numberOfAlpha; ++ i )
//	{
//		file >> rudderControlCoefArray[ i ];
//	}
//	for ( int i = 0; i < numberOfAlpha; ++ i )
//	{
//		for ( int j = 0; j < numberOfRudder; ++ j )
//		{
//			file >> cmzArray[ i ][ j ];
//		}
//	}
//}
//
//int FYDIControlLaw::FindInterpolationPosition( vector < RDouble > & targetArray, int lenthOfTargetArray, RDouble valueIn )
//{
//	int position = - 1;//如果只有一个数，没法插值，返回一个记号
//
//	if ( lenthOfTargetArray <= 2 )//如果最多只有两个数，也就不存在选择范围的问题
//	{
//		position = 0;
//	}
//	else
//	{
//		for ( int i = 0; i < lenthOfTargetArray - 1; ++ i )//需要i和i+1来插值，所以这里减1
//		{
//			if ( valueIn >= targetArray[ i ] && valueIn <= targetArray[ i + 1 ] 
//			  || valueIn <= targetArray[ i ] && valueIn >= targetArray[ i + 1 ] )
//			{
//				position = i;
//			}
//		}
//		//在起始值之左
//		if ( valueIn < targetArray[ 0 ] && targetArray[ 0 ] < targetArray[ 1 ]  //xIn小于起始值且递增
//		  || valueIn > targetArray[ 0 ] && targetArray[ 0 ] > targetArray[ 1 ] )//xIn大于起始值且递减
//		{
//			position = 0;
//		}
//		//在起始值之右
//		if ( valueIn < targetArray[ lenthOfTargetArray - 1 ] && targetArray[ lenthOfTargetArray - 1 ] < targetArray[ lenthOfTargetArray - 2 ]  //xIn大于末尾值且递增
//		  || valueIn > targetArray[ lenthOfTargetArray - 1 ] && targetArray[ lenthOfTargetArray - 1 ] > targetArray[ lenthOfTargetArray - 2 ] )//xIn小于末尾值且递减
//		{
//			position = lenthOfTargetArray - 2;
//		}
//	}
//
//	return position;
//}
//
//RDouble FYDIControlLaw::InterpolationCmde( RDouble alphaDeg )
//{
//	int xPosition = FindInterpolationPosition( alphaArray, numberOfAlpha, alphaDeg );
//
//	if ( xPosition >= 0 )
//	{
//		RDouble coefX = ( alphaDeg - alphaArray[ xPosition ] ) / ( alphaArray[ xPosition + 1 ] - alphaArray[ xPosition ] );
//		return rudderControlCoefArray[ xPosition ] + coefX * ( rudderControlCoefArray[ xPosition + 1 ] - rudderControlCoefArray[ xPosition ] );
//	}
//	else
//	{
//		return rudderControlCoefArray[ 0 ];
//	}
//}
//
//RDouble FYDIControlLaw::InterpolationCmz( RDouble alphaDeg, RDouble rudderDeg )
//{
//	int xPosition = FindInterpolationPosition( rudderArray, numberOfRudder, rudderDeg );
//	int yPosition = FindInterpolationPosition( alphaArray , numberOfAlpha , alphaDeg  );
//
//	int     iTempValue     =      0  ;
//	RDouble tempValue[ 2 ] = { 0, 0 };
//
//	for ( int tempPosition = yPosition; tempPosition <= yPosition + 1; ++ tempPosition )
//	{
//		RDouble coefX = ( rudderDeg - rudderArray[ xPosition ]) / ( rudderArray[ xPosition + 1 ] - rudderArray[ xPosition ] );
//		tempValue[ iTempValue ] = cmzArray[ tempPosition ][ xPosition     ]
//		              + coefX * ( cmzArray[ tempPosition ][ xPosition + 1 ]
//		                        - cmzArray[ tempPosition ][ xPosition     ] );
//		++ iTempValue;
//	}
//
//	RDouble coefY = ( alphaDeg - alphaArray[ yPosition ] ) / ( alphaArray[ yPosition + 1 ] - alphaArray[ yPosition ] );
//	return tempValue[ 0 ] + coefY * ( tempValue[ 1 ] - tempValue[ 0 ] );
//}
//
//void FYDIControlLaw::RudderRotationByDIControlLaw()
//{
//	//N-1时刻舵姿态
//	SixDofParameter * sixDofParameter1 = this->GetSixDofParameter( 1 );
//	RDouble oldMassCenterX     =   sixDofParameter1->GetMassCenterX    ();
//	RDouble oldMassCenterY     =   sixDofParameter1->GetMassCenterY    ();
//	RDouble oldMassCenterZ     =   sixDofParameter1->GetMassCenterZ    ();
//	RDouble oldEulerAngleAlpha =   sixDofParameter1->GetEulerAngleAlpha();
//	RDouble oldRudderAngle     = - sixDofParameter1->GetEulerAngleBeta ();//采用飞行力学坐标系，所以舵偏角应该反号
//	RDouble oldEulerAngleGama  =   sixDofParameter1->GetEulerAngleGama ();
//	//当前全弹姿态
//	int fartherIndex = GetIntegerParameterFromDataBase( bodyIndex, "fartherIndex" );
//	SixDofParameter * sixDofParameterFarther =   sixDofSolverManager   ->GetSixDofParameter ( fartherIndex );
//	RDouble eulerAngleAlphaFarther             = - sixDofParameterFarther->GetEulerAngleAlpha ();//采用飞行力学坐标系，俯仰角应该反号
//	RDouble angularVelocityZFarther            = - sixDofParameterFarther->GetAngularVelocityZ();//采用飞行力学坐标系，俯仰角速率应该反号
//
//	RDouble rudderAngle = 0.0;//n时刻舵偏角
//
//	ModifyTargetStateByStages(); targetState *= DegToRad;//输入的是攻角，与飞行力学坐标系同号
//
//	//当前攻角下舵零偏位置的俯仰力矩
//	RDouble momentCmz    = nodimensionalForceToDimensionalForce * InterpolationCmz ( eulerAngleAlphaFarther * RadToDeg, 0.0 );
//	//当前攻角下舵零偏位置的操纵力矩							 
//	RDouble momentRudder = nodimensionalForceToDimensionalForce * InterpolationCmde( eulerAngleAlphaFarther * RadToDeg );
//
//	//角速度指令积分项计算
//	integralValue = integralValue + dimensionalAleTimeStep * ( eulerAngleAlphaFarther - targetState );
//	//角速度积分池限制
//	integralValue = Saturation( integralValue, integralLoopLimit );
//
//	//慢回路给出的角速度指令
//	RDouble angularVelocityCommand = - k1 * ( eulerAngleAlphaFarther - targetState ) - k3 * integralValue;
//	//角速度指令限制
//	angularVelocityCommand         = Saturation( angularVelocityCommand, otherLimit );
//
//	//快回路给出的舵偏角指令
//	rudderAngle = ( - momentOfInertiaTensorZZ * k2 * ( angularVelocityZFarther - angularVelocityCommand ) - momentCmz ) / momentRudder * DegToRad;
//
//	//限制舵偏角角度限制
//	rudderAngle = Saturation( rudderAngle, rudderLowerLimit, rudderUpperLimit );
//	//限制舵偏角速率限制
//	rudderAngle = Saturation( rudderAngle - oldRudderAngle, rudderVelocity ) + oldRudderAngle;
//
//	cout << "eulerAngleAlphaFarther  = " << eulerAngleAlphaFarther  * RadToDeg << endl;
//	cout << "angularVelocityZFarther = " << angularVelocityZFarther * RadToDeg << endl;
//	cout << "rudderAngle             = " << rudderAngle             * RadToDeg << endl;
//	cout << "integralValue%          = " << integralValue / integralLoopLimit  << endl;
//
//	//将要重新填充的N时刻舵姿态
//	SixDofParameter * sixDofParameter0 = this->GetSixDofParameter( 0 );
//	sixDofParameter0->SetMassCenterX    ( oldMassCenterX     );
//	sixDofParameter0->SetMassCenterY    ( oldMassCenterY     );
//	sixDofParameter0->SetMassCenterZ    ( oldMassCenterZ     );
//	sixDofParameter0->SetEulerAngleAlpha( oldEulerAngleAlpha );
//	sixDofParameter0->SetEulerAngleBeta ( - rudderAngle      );//此时应从飞行力学坐标系转回CFD坐标系，所以应该反号
//	sixDofParameter0->SetEulerAngleGama ( oldEulerAngleGama  );
//}
//
//void FYDIControlLaw::RudderRotationByDI089ControlLaw()
//{
//	//N-1时刻舵姿态
//	SixDofParameter * sixDofParameter1 = this->GetSixDofParameter( 1 );
//	RDouble oldMassCenterX     =   sixDofParameter1->GetMassCenterX    ();
//	RDouble oldMassCenterY     =   sixDofParameter1->GetMassCenterY    ();
//	RDouble oldMassCenterZ     =   sixDofParameter1->GetMassCenterZ    ();
//	RDouble oldRudderAngle     = - sixDofParameter1->GetEulerAngleAlpha();//采用飞行力学坐标系，所以舵偏角应该反号
//	RDouble oldEulerAngleBeta  =   sixDofParameter1->GetEulerAngleBeta ();
//	RDouble oldEulerAngleGama  =   sixDofParameter1->GetEulerAngleGama ();
//	//当前全弹姿态
//	int fartherIndex = GetIntegerParameterFromDataBase( bodyIndex, "fartherIndex" );
//	SixDofParameter * sixDofParameterFarther =   sixDofSolverManager   ->GetSixDofParameter ( fartherIndex );
//	RDouble eulerAngleAlphaFarther             = - sixDofParameterFarther->GetEulerAngleAlpha ();//采用飞行力学坐标系，俯仰角应该反号
//	RDouble angularVelocityZFarther            = - sixDofParameterFarther->GetAngularVelocityZ();//采用飞行力学坐标系，俯仰角速率应该反号
//
//	RDouble rudderAngle = 0.0;//n时刻舵偏角
//
//	ModifyTargetStateByStages(); targetState *= DegToRad;//输入的是攻角，与飞行力学坐标系同号
//
//	//当前攻角下舵零偏位置的俯仰力矩
//	RDouble momentCmz    = nodimensionalForceToDimensionalForce * InterpolationCmz ( eulerAngleAlphaFarther * RadToDeg, 0.0 );
//	//当前攻角下舵零偏位置的操纵力矩							 
//	RDouble momentRudder = nodimensionalForceToDimensionalForce * InterpolationCmde( eulerAngleAlphaFarther * RadToDeg );
//
//	//角速度指令积分项计算
//	integralValue = integralValue + dimensionalAleTimeStep * ( eulerAngleAlphaFarther - targetState );
//	//角速度积分池限制
//	integralValue = Saturation( integralValue, integralLoopLimit );
//
//	//慢回路给出的角速度指令
//	RDouble angularVelocityCommand = - k1 * ( eulerAngleAlphaFarther - targetState ) - k3 * integralValue;
//	//角速度指令限制
//	angularVelocityCommand         = Saturation( angularVelocityCommand, otherLimit );
//
//	//快回路给出的舵偏角指令
//	rudderAngle = ( - momentOfInertiaTensorZZ * k2 * ( angularVelocityZFarther - angularVelocityCommand ) - momentCmz ) / momentRudder * DegToRad;
//
//	//限制舵偏角角度限制
//	rudderAngle = Saturation( rudderAngle, rudderLowerLimit, rudderUpperLimit );
//	//限制舵偏角速率限制
//	rudderAngle = Saturation( rudderAngle - oldRudderAngle, rudderVelocity ) + oldRudderAngle;
//
//	cout << "eulerAngleAlphaFarther  = " << eulerAngleAlphaFarther  * RadToDeg << endl;
//	cout << "angularVelocityZFarther = " << angularVelocityZFarther * RadToDeg << endl;
//	cout << "rudderAngle             = " << rudderAngle             * RadToDeg << endl;
//	cout << "integralValue%          = " << integralValue / integralLoopLimit  << endl;
//
//	//将要重新填充的N时刻舵姿态
//	SixDofParameter * sixDofParameter0 = this->GetSixDofParameter( 0 );
//	sixDofParameter0->SetMassCenterX    ( oldMassCenterX    );
//	sixDofParameter0->SetMassCenterY    ( oldMassCenterY    );
//	sixDofParameter0->SetMassCenterZ    ( oldMassCenterZ    );
//	sixDofParameter0->SetEulerAngleAlpha( - rudderAngle     );//此时应从飞行力学坐标系转回CFD坐标系，所以应该反号
//	sixDofParameter0->SetEulerAngleBeta ( oldEulerAngleBeta );
//	sixDofParameter0->SetEulerAngleGama ( oldEulerAngleGama );
//}
//
//void FYPIDControlLaw::Initialize( int bodyIndexIn )
//{
//	SixDofSolver::Initialize( bodyIndexIn );
//
//	int fartherIndex = GetIntegerParameterFromDataBase( bodyIndex, "fartherIndex" );
//	MassCharacter * massCharacterFarther = sixDofSolverManager->GetMassCharacter( fartherIndex );
//	massOfMissile   = massCharacterFarther->GetMass();
//
//	//本类采用飞行力学坐标系
//	k1                = k1;//KA
//	k2                = k2;//WI
//	k3                = k3;//KR
//
//	rudderLowerLimit  = rudderLowerLimit  * DegToRad;//上偏极限（上偏为负）
//	rudderUpperLimit  = rudderUpperLimit  * DegToRad;//下偏极限（下偏为正）
//	rudderVelocity    = rudderVelocity    * DegToRad * dimensionalAleTimeStep;//从给定的每秒偏转值转换到每步
//
//	integralLoopLimit = integralLoopLimit * DegToRad;//积分池大小限制
//	otherLimit        = otherLimit        * DegToRad;//大攻角舵效限制
//}
//
//void FYPIDControlLaw::RudderRotationByPIDAttitudeControlLaw()
//{
//	//N-1时刻舵姿态
//	SixDofParameter * sixDofParameter1 = this->GetSixDofParameter( 1 );
//	RDouble oldMassCenterX     =   sixDofParameter1->GetMassCenterX    ();
//	RDouble oldMassCenterY     =   sixDofParameter1->GetMassCenterY    ();
//	RDouble oldMassCenterZ     =   sixDofParameter1->GetMassCenterZ    ();
//	RDouble oldEulerAngleAlpha =   sixDofParameter1->GetEulerAngleAlpha();
//	RDouble oldRudderAngle     = - sixDofParameter1->GetEulerAngleBeta ();//采用飞行力学坐标系，所以舵偏角应该反号
//	RDouble oldEulerAngleGama  =   sixDofParameter1->GetEulerAngleGama ();
//	
//	//当前全弹姿态
//	int fartherIndex = GetIntegerParameterFromDataBase( bodyIndex, "fartherIndex" );
//	SixDofParameter * sixDofParameterFarther =   sixDofSolverManager   ->GetSixDofParameter ( fartherIndex );
//	RDouble eulerAngleAlphaFarther             = - sixDofParameterFarther->GetEulerAngleAlpha ();//采用飞行力学坐标系，俯仰角应该反号
//	RDouble angularVelocityZFarther            = - sixDofParameterFarther->GetAngularVelocityZ();//采用飞行力学坐标系，俯仰角速率应该反号
//	cout << "eulerAngleAlphaFarther  = " << eulerAngleAlphaFarther  * RadToDeg << endl;
//	cout << "angularVelocityZFarther = " << angularVelocityZFarther * RadToDeg << endl;
//	eulerAngleAlphaFarther  = atan( tan( eulerAngleAlphaFarther ) * sqrt( 2.0 ) / 2.0 );
//	angularVelocityZFarther = angularVelocityZFarther             * sqrt( 2.0 ) / 2.0  ;
//
//	ModifyTargetStateByStages();targetState *= DegToRad;//输入的是攻角，与飞行力学坐标系同号
//	targetState             = atan( tan( targetState            ) * sqrt( 2.0 ) / 2.0 );
//
//	RDouble rudderAngle = 0.0;//n时刻舵偏角
//
//	RDouble errorAttackAngle      = eulerAngleAlphaFarther - targetState;
//	RDouble integralIncreaseValue = ( errorAttackAngle * k1 + angularVelocityZFarther ) * dimensionalAleTimeStep;
//	integralValue                += integralIncreaseValue * k2;
//
//	//积分池限制
//	integralValue = Saturation( integralValue, integralLoopLimit );
//	//限制舵偏角角度限制
//	rudderAngle   = Saturation( ( integralValue + angularVelocityZFarther ) * k3, rudderLowerLimit, rudderUpperLimit );
//	//大攻角舵效影响限制
//	rudderAngle   = Saturation( eulerAngleAlphaFarther + rudderAngle, otherLimit ) - eulerAngleAlphaFarther;
//	//限制舵偏角速率限制
//	rudderAngle   = Saturation( rudderAngle - oldRudderAngle, rudderVelocity ) + oldRudderAngle;
//	//积分池限制
//	integralValue = Saturation( rudderAngle / k3 - angularVelocityZFarther, integralLoopLimit );
//
//	rudderAngle   = ( integralValue + angularVelocityZFarther ) * k3;
//
//	cout << "rudderAngle             = " << rudderAngle             * RadToDeg << endl;
//	cout << "integralValue%          = " << integralValue / integralLoopLimit  << endl;
//
//	//将要重新填充的N时刻舵姿态
//	SixDofParameter * sixDofParameter0 = this->GetSixDofParameter( 0 );
//	sixDofParameter0->SetMassCenterX    ( oldMassCenterX     );
//	sixDofParameter0->SetMassCenterY    ( oldMassCenterY     );
//	sixDofParameter0->SetMassCenterZ    ( oldMassCenterZ     );
//	sixDofParameter0->SetEulerAngleAlpha( oldEulerAngleAlpha );
//	sixDofParameter0->SetEulerAngleBeta ( - rudderAngle      );//此时应从飞行力学坐标系转回CFD坐标系，所以应该反号
//	sixDofParameter0->SetEulerAngleGama ( oldEulerAngleGama  );
//}
//
//void FYPIDControlLaw::RudderRotationByPIDAccelerationControlLaw()
//{
//	//N-1时刻舵姿态
//	SixDofParameter * sixDofParameter1 = this->GetSixDofParameter( 1 );
//	RDouble oldMassCenterX     =   sixDofParameter1->GetMassCenterX    ();
//	RDouble oldMassCenterY     =   sixDofParameter1->GetMassCenterY    ();
//	RDouble oldMassCenterZ     =   sixDofParameter1->GetMassCenterZ    ();
//	RDouble oldEulerAngleAlpha =   sixDofParameter1->GetEulerAngleAlpha();
//	RDouble oldRudderAngle     = - sixDofParameter1->GetEulerAngleBeta ();//采用飞行力学坐标系，所以舵偏角应该反号
//	RDouble oldEulerAngleGama  =   sixDofParameter1->GetEulerAngleGama ();
//	
//	//当前全弹姿态
//	int fartherIndex = GetIntegerParameterFromDataBase( bodyIndex, "fartherIndex" );
//	SixDofParameter * sixDofParameterFarther =   sixDofSolverManager   ->GetSixDofParameter ( fartherIndex );
//	RDouble eulerAngleAlphaFarther             = - sixDofParameterFarther->GetEulerAngleAlpha ();//采用飞行力学坐标系，俯仰角应该反号
//	RDouble angularVelocityZFarther            = - sixDofParameterFarther->GetAngularVelocityZ();//采用飞行力学坐标系，俯仰角速率应该反号
//	cout << "eulerAngleAlphaFarther  = " << eulerAngleAlphaFarther  * RadToDeg << endl;
//	cout << "angularVelocityZFarther = " << angularVelocityZFarther * RadToDeg << endl;
//	eulerAngleAlphaFarther  = atan( tan( eulerAngleAlphaFarther ) * sqrt( 2.0 ) / 2.0 );
//	angularVelocityZFarther = angularVelocityZFarther             * sqrt( 2.0 ) / 2.0  ;
//
//	vector< RDouble > aleForceContainerFather( 6 );
//	SixDofSolver * sixDofSolverFarther = sixDofSolverManager->GetSixDofSolver( fartherIndex );
//	sixDofSolverFarther->GetDimensionalAerodynamicForceAndMomentBodyFrame( aleForceContainerFather, 1 );
//	RDouble accelerationN1 = aleForceContainerFather[ 1 ]         * sqrt( 2.0 ) / 2.0 / massOfMissile;
//
//	ModifyTargetStateByStages(); targetState = targetState        * sqrt( 2.0 ) / 2.0  ;//指令过载
//
//	RDouble rudderAngle = 0.0;  //n时刻舵偏角
//
//	RDouble errorAttackAngle      = accelerationN1 - targetState;
//	RDouble integralIncreaseValue = ( errorAttackAngle * k1 + angularVelocityZFarther ) * dimensionalAleTimeStep;
//	integralValue                += integralIncreaseValue * k2;
//
//	//积分池限制
//	integralValue = Saturation( integralValue, integralLoopLimit );
//	//限制舵偏角角度限制
//	rudderAngle   = Saturation( ( integralValue + angularVelocityZFarther ) * k3, rudderLowerLimit, rudderUpperLimit );
//	//大攻角舵效影响限制
//	rudderAngle   = Saturation( eulerAngleAlphaFarther + rudderAngle, otherLimit ) - eulerAngleAlphaFarther;
//	//限制舵偏角速率限制
//	rudderAngle   = Saturation( rudderAngle - oldRudderAngle, rudderVelocity ) + oldRudderAngle;
//	//积分池限制
//	integralValue = Saturation( rudderAngle / k3 - angularVelocityZFarther, integralLoopLimit );
//
//	rudderAngle   = ( integralValue + angularVelocityZFarther ) * k3;
//
//	cout << "rudderAngle             = " << rudderAngle             * RadToDeg << endl;
//	cout << "integralValue%          = " << integralValue / integralLoopLimit  << endl;
//
//	//将要重新填充的N时刻舵姿态
//	SixDofParameter * sixDofParameter0 = this->GetSixDofParameter( 0 );
//	sixDofParameter0->SetMassCenterX    ( oldMassCenterX     );
//	sixDofParameter0->SetMassCenterY    ( oldMassCenterY     );
//	sixDofParameter0->SetMassCenterZ    ( oldMassCenterZ     );
//	sixDofParameter0->SetEulerAngleAlpha( oldEulerAngleAlpha );
//	sixDofParameter0->SetEulerAngleBeta ( - rudderAngle      );//此时应从飞行力学坐标系转回CFD坐标系，所以应该反号
//	sixDofParameter0->SetEulerAngleGama ( oldEulerAngleGama  );
//}
//
//void FYPIDControlLaw::RudderRotationByPID089ControlLaw()
//{
//	//N-1时刻舵姿态
//	SixDofParameter * sixDofParameter1 = this->GetSixDofParameter( 1 );
//	RDouble oldMassCenterX     =   sixDofParameter1->GetMassCenterX    ();
//	RDouble oldMassCenterY     =   sixDofParameter1->GetMassCenterY    ();
//	RDouble oldMassCenterZ     =   sixDofParameter1->GetMassCenterZ    ();
//	RDouble oldRudderAngle     = - sixDofParameter1->GetEulerAngleAlpha();//采用飞行力学坐标系，所以舵偏角应该反号
//	RDouble oldEulerAngleBeta  =   sixDofParameter1->GetEulerAngleBeta ();
//	RDouble oldEulerAngleGama  =   sixDofParameter1->GetEulerAngleGama ();
//	
//	//当前全弹姿态
//	int fartherIndex = GetIntegerParameterFromDataBase( bodyIndex, "fartherIndex" );
//	SixDofParameter * sixDofParameterFarther =   sixDofSolverManager   ->GetSixDofParameter ( fartherIndex );
//	RDouble eulerAngleAlphaFarther             = - sixDofParameterFarther->GetEulerAngleAlpha ();//采用飞行力学坐标系，俯仰角应该反号
//	RDouble angularVelocityZFarther            = - sixDofParameterFarther->GetAngularVelocityZ();//采用飞行力学坐标系，俯仰角速率应该反号
//	cout << "eulerAngleAlphaFarther  = " << eulerAngleAlphaFarther  * RadToDeg << endl;
//	cout << "angularVelocityZFarther = " << angularVelocityZFarther * RadToDeg << endl;
//
//	RDouble dimensionalAleTime = GetDimensionalAleSimulationTime();
//	if ( dimensionalAleTime < 1.0 )
//	{
//		targetState = 14.75;
//	}
//	else if ( dimensionalAleTime < 3.0 )
//	{
//		targetState = 13.395 + dimensionalAleTime * 1.355;
//	}
//	else if ( dimensionalAleTime >= 3.0 )
//	{
//		targetState = 17.46;
//	}
//	else
//	{
//		targetState = 0.0;
//	}
//	targetState *= DegToRad;//输入的是攻角，与飞行力学坐标系同号
//
//	RDouble rudderAngle = 0.0;//n时刻舵偏角
//
//	RDouble errorAttackAngle = eulerAngleAlphaFarther - targetState;
//	integralValue           += errorAttackAngle * dimensionalAleTimeStep * k2;
//
//	//积分池限制
//	integralValue = Saturation( integralValue, integralLoopLimit );
//	//限制舵偏角角度限制
//	rudderAngle   = Saturation( errorAttackAngle * k1 + angularVelocityZFarther * k3 + integralValue, rudderLowerLimit, rudderUpperLimit );
//	//限制舵偏角速率限制
//	rudderAngle   = Saturation( rudderAngle - oldRudderAngle, rudderVelocity ) + oldRudderAngle;
//
//	cout << "rudderAngle             = " << rudderAngle             * RadToDeg << endl;
//	cout << "integralValue%          = " << integralValue / integralLoopLimit  << endl;
//
//	//将要重新填充的N时刻舵姿态
//	SixDofParameter * sixDofParameter0 = this->GetSixDofParameter( 0 );
//	sixDofParameter0->SetMassCenterX    ( oldMassCenterX    );
//	sixDofParameter0->SetMassCenterY    ( oldMassCenterY    );
//	sixDofParameter0->SetMassCenterZ    ( oldMassCenterZ    );
//	sixDofParameter0->SetEulerAngleAlpha( - rudderAngle     );//此时应从飞行力学坐标系转回CFD坐标系，所以应该反号
//	sixDofParameter0->SetEulerAngleBeta ( oldEulerAngleBeta );
//	sixDofParameter0->SetEulerAngleGama ( oldEulerAngleGama );
//}
//
//void FYPID2DControlLaw::Initialize( int bodyIndexIn )
//{
//	SixDofSolver::Initialize( bodyIndexIn );
//
//	//本类采用飞行力学坐标系
//	k1A                          = k1A                                                             ;//俯仰偏航通道KA
//	k2A                          = k2A                                                             ;//俯仰偏航通道WI
//	k3A                          = k3A                                                             ;//俯仰偏航通道KR
//					       
//	integralLoopLimitA           = integralLoopLimitA           * DegToRad                         ;//俯仰偏航通道积分池限制
//	largeAttackAngleLimitA       = largeAttackAngleLimitA       * DegToRad                         ;//俯仰偏航通道大攻角舵效限制
//					       
//	k1B                          = k1B                                                             ;//滚转通道KGI
//	k2B                          = k2B                                                             ;//滚转通道KGA
//	k3B                          = k3B                                                             ;//滚转通道IX
//
//	deltaTargetStateLimitB       = deltaTargetStateLimitB       * DegToRad                         ;//滚转通道指令差量限制
//	integralIncreaseLimitB       = integralIncreaseLimitB       * DegToRad * dimensionalAleTimeStep;//滚转通道积分增量限制(与时间步长相关)
//	integralLoopLimitB           = integralLoopLimitB           * DegToRad                         ;//滚转通道积分池限制
//	angularVelocityCommandLimitB = angularVelocityCommandLimitB * DegToRad                         ;//滚转通道角速度指令限制
//
//	rudderVelocityLimit          = rudderVelocityLimit          * DegToRad * dimensionalAleTimeStep;//舵偏速率限制(与时间步长相关)
//	rudderAngleLimitA            = rudderAngleLimitA            * DegToRad                         ;//俯仰偏航通道舵偏限制
//	rudderAngleLimitB            = rudderAngleLimitB            * DegToRad                         ;//滚转通道舵偏限制
//	rudderAngleLimit             = rudderAngleLimit             * DegToRad                         ;//耦合舵偏限制
//}
//
//void FYPID2DControlLaw::RudderRotationByPID2DMissileControlLaw()
//{
//	//N-1时刻舵姿态
//	SixDofParameter * sixDofParameter1 = this->GetSixDofParameter( 1 );
//	RDouble oldMassCenterX     =   sixDofParameter1->GetMassCenterX    ();
//	RDouble oldMassCenterY     =   sixDofParameter1->GetMassCenterY    ();
//	RDouble oldMassCenterZ     =   sixDofParameter1->GetMassCenterZ    ();
//	RDouble oldEulerAngleAlpha =   sixDofParameter1->GetEulerAngleAlpha();
//	RDouble oldRudderAngle     = - sixDofParameter1->GetEulerAngleBeta ();//采用飞行力学坐标系，所以舵偏角应该反号
//	RDouble oldEulerAngleGama  =   sixDofParameter1->GetEulerAngleGama ();
//
//	//当前全机姿态
//	int fartherIndex = GetIntegerParameterFromDataBase( bodyIndex, "fartherIndex" );
//	SixDofParameter * sixDofParameterFarther =   sixDofSolverManager   ->GetSixDofParameter ( fartherIndex );
//	RDouble eulerAngleAlphaFarther             = - sixDofParameterFarther->GetEulerAngleAlpha ();//采用飞行力学坐标系，俯仰角应该反号
//	RDouble angularVelocityZFarther            = - sixDofParameterFarther->GetAngularVelocityZ();//采用飞行力学坐标系，俯仰角速率应该反号
//	RDouble eulerAngleGamaFarther              = - sixDofParameterFarther->GetEulerAngleGama  ();//采用飞行力学坐标系，滚转角应该反号
//	RDouble angularVelocityXFarther            = - sixDofParameterFarther->GetAngularVelocityX();//采用飞行力学坐标系，滚转角速率应该反号
//
//	//分解至俯仰偏航通道
//	RDouble eulerAngleAlphaFarther1  = atan( tan( eulerAngleAlphaFarther ) * cos( PI / 4.0 + eulerAngleGamaFarther ) );//俯仰
//	RDouble eulerAngleAlphaFarther2  = atan( tan( eulerAngleAlphaFarther ) * sin( PI / 4.0 + eulerAngleGamaFarther ) );//偏航
//	RDouble angularVelocityZFarther1 = angularVelocityZFarther             * cos( PI / 4.0 + eulerAngleGamaFarther )  ;//俯仰
//	RDouble angularVelocityZFarther2 = angularVelocityZFarther             * sin( PI / 4.0 + eulerAngleGamaFarther )  ;//偏航
//
//	//控制指令获取并分解至俯仰偏航通道
//	ModifyTargetStateByStages();
//	targetStateA *= DegToRad;//输入的是攻角，与飞行力学坐标系同号
//	targetStateB *= DegToRad;//输入的是飞行力学坐标系下的滚转角，左侧上抬为正
//	targetStateA1                    = atan( tan( targetStateA           ) * cos( PI / 4.0 + eulerAngleGamaFarther ) );//俯仰
//	targetStateA2                    = atan( tan( targetStateA           ) * sin( PI / 4.0 + eulerAngleGamaFarther ) );//偏航
//
//	////俯仰偏航通道
//	//俯仰偏航通道积分项
//	RDouble errorAttackAngle1       = eulerAngleAlphaFarther1 - targetStateA1;
//	RDouble errorAttackAngle2       = eulerAngleAlphaFarther2 - targetStateA2;
//	RDouble integralIncreaseValueA1 = ( errorAttackAngle1 * k1A + angularVelocityZFarther1 ) * dimensionalAleTimeStep;
//	RDouble integralIncreaseValueA2 = ( errorAttackAngle2 * k1A + angularVelocityZFarther2 ) * dimensionalAleTimeStep;
//	integralValueA1                += integralIncreaseValueA1 * k2A;
//	integralValueA2                += integralIncreaseValueA2 * k2A;
//	//俯仰偏航通道积分池限制
//	integralValueA1  = Saturation( integralValueA1, integralLoopLimitA );
//	integralValueA2  = Saturation( integralValueA2, integralLoopLimitA );
//	//俯仰偏航通道舵偏限制
//	rudderAngleA1    = Saturation( ( integralValueA1 + angularVelocityZFarther1 ) * k3A, rudderAngleLimitA );
//	rudderAngleA2    = Saturation( ( integralValueA2 + angularVelocityZFarther2 ) * k3A, rudderAngleLimitA );
//	//俯仰偏航通道大攻角舵效限制
//	rudderAngleA1    = Saturation( eulerAngleAlphaFarther1 + rudderAngleA1, largeAttackAngleLimitA ) - eulerAngleAlphaFarther1;
//	rudderAngleA2    = Saturation( eulerAngleAlphaFarther2 + rudderAngleA2, largeAttackAngleLimitA ) - eulerAngleAlphaFarther2;
//	//俯仰偏航舵偏速率限制
//	rudderAngleA1    = Saturation( rudderAngleA1 - oldRudderAngleA1, rudderVelocityLimit ) + oldRudderAngleA1;
//	rudderAngleA2    = Saturation( rudderAngleA2 - oldRudderAngleA2, rudderVelocityLimit ) + oldRudderAngleA2;
//	oldRudderAngleA1 = rudderAngleA1;
//	oldRudderAngleA2 = rudderAngleA2;
//	//俯仰偏航通道积分池限制
//	integralValueA1  = Saturation( rudderAngleA1 / k3A - angularVelocityZFarther1, integralLoopLimitA );
//	integralValueA2  = Saturation( rudderAngleA2 / k3A - angularVelocityZFarther2, integralLoopLimitA );
//	//俯仰偏航舵偏
//	rudderAngleA1    = ( integralValueA1 + angularVelocityZFarther1 ) * k3A;
//	rudderAngleA2    = ( integralValueA2 + angularVelocityZFarther2 ) * k3A;
//
//	////滚转通道
//	//滚转通道指令差量限制
//	RDouble errorRollingAngle       = Saturation( eulerAngleGamaFarther - targetStateB, deltaTargetStateLimitB );
//	//滚转通道积分增量限制
//	RDouble integralIncreaseValueB  = Saturation( errorRollingAngle * k1B * dimensionalAleTimeStep, integralIncreaseLimitB );
//	//滚转通道积分池限制
//	integralValueB                 += integralIncreaseValueB;
//	integralValueB                  = Saturation( integralValueB, integralLoopLimitB );
//	//滚转通道角速度指令限制
//	RDouble angularVelocityCommandB = Saturation( integralValueB + errorRollingAngle * k2B, angularVelocityCommandLimitB );
//	//滚转通道舵偏限制
//	rudderAngleB                    = Saturation( ( angularVelocityCommandB + angularVelocityXFarther ) * k3B, rudderAngleLimitB );
//	//滚转通道舵偏速率限制
//	rudderAngleB                    = Saturation( rudderAngleB - oldRudderAngleB, rudderVelocityLimit ) + oldRudderAngleB;
//	oldRudderAngleB                 = rudderAngleB;
//	////滚转通道积分池限制
//	//integralValueB                  = Saturation( rudderAngleB / k3B - angularVelocityXFarther - errorRollingAngle * k2B, integralLoopLimitB );
//	////滚转通道舵偏
//	//rudderAngleB                    = ( integralValueB + errorRollingAngle * k2B + angularVelocityXFarther ) * k3B;
//
//	cout << "eulerAngleAlphaFarther  = " << eulerAngleAlphaFarther    * RadToDeg << endl;
//	cout << "angularVelocityZFarther = " << angularVelocityZFarther   * RadToDeg << endl;
//	cout << "eulerAngleGamaFarther   = " << eulerAngleGamaFarther     * RadToDeg << endl;
//	cout << "angularVelocityXFarther = " << angularVelocityXFarther   * RadToDeg << endl;
//
//	cout << "rudderAngleA1           = " << rudderAngleA1             * RadToDeg << endl;
//	cout << "integralValueA1%        = " << integralValueA1 / integralLoopLimitA << endl;
//	cout << "rudderAngleA2           = " << rudderAngleA2             * RadToDeg << endl;
//	cout << "integralValueA2%        = " << integralValueA2 / integralLoopLimitA << endl;
//	cout << "rudderAngleB            = " << rudderAngleB              * RadToDeg << endl;
//	cout << "integralValueB%         = " << integralValueB  / integralLoopLimitB << endl;
//
//	//通道耦合
//	RDouble rudderAngle = 0.0;
//	if      ( bodyIndex == 2 )
//	{
//		rudderAngle = rudderAngleA1 - rudderAngleB;
//	}
//	else if ( bodyIndex == 3 )
//	{
//		rudderAngle = rudderAngleA2 - rudderAngleB;
//	}
//	else if ( bodyIndex == 4 )
//	{
//		rudderAngle = rudderAngleA2 + rudderAngleB;
//	}
//	else if ( bodyIndex == 5 )
//	{
//		rudderAngle = rudderAngleA1 + rudderAngleB;
//	}
//	else
//	{
//		StopProgram( "无法确定未知舵面等效舵偏到实际舵偏的组合关系" );
//	}
//
//	//耦合舵偏限制
//	rudderAngle = Saturation( rudderAngle, rudderAngleLimit );
//	//耦合舵偏速率限制
//	rudderAngle = Saturation( rudderAngle - oldRudderAngle, rudderVelocityLimit ) + oldRudderAngle;
//
//	cout << "rudderAngle             = " << rudderAngle               * RadToDeg << endl;
//
//	//将要重新填充的N时刻舵姿态
//	SixDofParameter * sixDofParameter0 = this->GetSixDofParameter( 0 );
//	sixDofParameter0->SetMassCenterX    ( oldMassCenterX     );
//	sixDofParameter0->SetMassCenterY    ( oldMassCenterY     );
//	sixDofParameter0->SetMassCenterZ    ( oldMassCenterZ     );
//	sixDofParameter0->SetEulerAngleAlpha( oldEulerAngleAlpha );
//	sixDofParameter0->SetEulerAngleBeta ( - rudderAngle      );//此时应从飞行力学坐标系转回CFD坐标系，所以应该反号
//	sixDofParameter0->SetEulerAngleGama ( oldEulerAngleGama  );
//}

SixDofSolver::SixDofSolver( SixDofSolverManager * sixDofSolverManager )
{
    this->sixDofSolverManager = sixDofSolverManager;

    bodyIndex = 0;

    AllocateAllVariables();

    DegToRad = PI / 180.0;
    RadToDeg = 180.0 / PI;
}

SixDofSolver::~SixDofSolver()
{
    DeallocateAllVariables();
}

void SixDofSolver::AllocateAllVariables()
{
    AllocateSixDofParameter();

    AllocateSixDofForceContainer();

    AllocateMassCharacter();
}

void SixDofSolver::DeallocateAllVariables()
{
    int numberOfTimeLayers = this->GetNumberOfTimeLayers();
    for ( int iTimeLayer = 0; iTimeLayer < numberOfTimeLayers; ++ iTimeLayer )
    {
        delete sixDofParameterContainer[ iTimeLayer ];
        sixDofParameterContainer[iTimeLayer] = NULL;
    }

    delete massCharacter;
    massCharacter = NULL;
}

void SixDofSolver::AllocateSixDofParameter()
{
    int numberOfTimeLayers = 3;
    SetNumberOfTimeLayers( numberOfTimeLayers );

    sixDofParameterContainer.resize( numberOfTimeLayers );

    for ( int iTimeLayer = 0; iTimeLayer < numberOfTimeLayers; ++ iTimeLayer )
    {
        sixDofParameterContainer[ iTimeLayer ] = new SixDofParameter();
    }
}

void SixDofSolver::AllocateMassCharacter()
{
    massCharacter = new MassCharacter();
}

void SixDofSolver::AllocateSixDofForceContainer()
{
    int numberOfTimeLayers = this->GetNumberOfTimeLayers();

    sixDofForceContainer.resize( numberOfTimeLayers );

    for ( int iTimeLayer = 0; iTimeLayer < numberOfTimeLayers; ++ iTimeLayer )
    {
        sixDofForceContainer[ iTimeLayer ].resize( 6 );
    }
}

void SixDofSolver::InitializeMassCharacter()
{
    MassCharacter * massCharacter = this->GetMassCharacter();

    RDouble mass = GetRDoubleParameterFromDataBase( bodyIndex, "mass" );

    vector< RDouble > massMatrix;
    GetRDoubleVectorFromDataBase( massMatrix, bodyIndex, "massMatrix", 6 );

    massCharacter->SetMassCharacter( mass, massMatrix );
}

void SixDofSolver::Initialize( int bodyIndexIn )
{
    this->bodyIndex = bodyIndexIn;

    this->InitializeMassCharacter();

    this->InitializeSixDofByParameter();
}

void SixDofSolver::Restart()
{
    fstream & file = PHSPACE::GetCurrentAleFile();
    this->InitializeSixDofByFiles( file );
}

void SixDofSolver::InitializeSixDofByParameter()
{
    vector< RDouble > massCenter;
    GetRDoubleVectorFromDataBase(massCenter, bodyIndex, "massCenter"  , 3 );

    vector< RDouble > attitudeAngle;
    GetRDoubleVectorFromDataBase(attitudeAngle, bodyIndex, "attitudeAngle", 3);

    vector< RDouble > massCenterVelocity;
    GetRDoubleVectorFromDataBase(massCenterVelocity, bodyIndex, "massCenterVelocity", 3 );

    vector< RDouble > angularVelocity;
    GetRDoubleVectorFromDataBase(angularVelocity, bodyIndex, "angularVelocity", 3);

    attitudeAngle[ 0 ] *= PI / 180;
    attitudeAngle[ 1 ] *= PI / 180;
    attitudeAngle[ 2 ] *= PI / 180;

    angularVelocity[ 0 ] *= PI / 180;
    angularVelocity[ 1 ] *= PI / 180;
    angularVelocity[ 2 ] *= PI / 180;

    for ( int iTimeLayer = 0; iTimeLayer < numberOfTimeLayers; ++ iTimeLayer )
    {
        SixDofParameter * sixDofParameter = this->GetSixDofParameter( iTimeLayer );
        sixDofParameter->SetGeneralVariableContainer(massCenter, attitudeAngle, massCenterVelocity, angularVelocity);
    }
}

void SixDofSolver::InitializeSixDofByFiles( fstream & file )
{
    const int numberOfEquations = 12;
    RDouble generalVariableContainer[ numberOfEquations ];

    for ( int iTimeLayer = 0; iTimeLayer < numberOfTimeLayers; ++ iTimeLayer )
    {
        PHSPACE::PHRead( file, generalVariableContainer, numberOfEquations );
        sixDofParameterContainer[ iTimeLayer ]->SetGeneralVariableContainer( generalVariableContainer );
    }
}

//void PHControlLaw::InitializeSixDofByFiles( fstream & file )
//{
//	SixDofSolver::InitializeSixDofByFiles( file );
//	PHSPACE::PHRead( file, integralValue );
//}
//
//void PH2DControlLaw::InitializeSixDofByFiles( fstream & file )
//{
//	SixDofSolver::InitializeSixDofByFiles( file );
//	PHSPACE::PHRead( file, integralValueA1  );
//	PHSPACE::PHRead( file, integralValueA2  );
//	PHSPACE::PHRead( file, integralValueB   );
//	PHSPACE::PHRead( file, oldRudderAngleA1 );
//	PHSPACE::PHRead( file, oldRudderAngleA2 );
//	PHSPACE::PHRead( file, oldRudderAngleB  );
//}

void SixDofSolver::Run()
{
    ComputeSixDofForceTimeSequenceContainer();

    int RBDMethod = PHSPACE::GetIntegerParameterFromDataBase( bodyIndex, "RBDMethod" );

    if      ( RBDMethod == - 1 )//给定运动规律
    {
        MovingByFile();
    }
    else if ( RBDMethod == 0   )//静止
    {
        Stationary();
    }
    else if ( RBDMethod == 1   )//六自由度运动
    {
        SolveSixDofEquation();
    }
    else if ( RBDMethod == 2   )//三自由度运动
    {
        ModifyForceAndSixDofInformationForXYSymmetryCase();
        SolveSixDofEquation();
    }
    else if ( RBDMethod == 11  )//X轴强迫运动
    {
        OscillatingInXDirection();
    }
    else if ( RBDMethod == 12  )//Y轴强迫运动
    {
        TK_Exit::ExceptionExit( "暂无Y轴强迫运动" );
        //OscillatingInYDirection();
    }
    else if ( RBDMethod == 13  )//Z轴强迫运动
    {
        TK_Exit::ExceptionExit( "暂无Z轴强迫运动" );
        //OscillatingInZDirection();
    }
    else if ( RBDMethod == 14  )//强迫俯仰运动
    {
        Pitching();
        //SolveSixDofEquation();
    }
    else if ( RBDMethod == 15  )//强迫偏航运动
    {
        Yawing();
    }
    else if ( RBDMethod == 16  )//强迫滚转运动
    {
        Rolling();
    }
    else if ( RBDMethod == 21  )//单自由度X轴运动
    {
        TK_Exit::ExceptionExit( "暂无单自由度X轴运动" );
        //ModifyForceAndSixDofInformationForXDirectionOnly();
        //SolveSixDofEquation();
    }
    else if ( RBDMethod == 22  )//单自由度Y轴运动
    {
        TK_Exit::ExceptionExit( "暂无单自由度Y轴运动" );
        //ModifyForceAndSixDofInformationForYDirectionOnly();
        //SolveSixDofEquation();
    }
    else if ( RBDMethod == 23  )//单自由度Z轴运动
    {
        TK_Exit::ExceptionExit( "暂无单自由度Z轴运动" );
        //ModifyForceAndSixDofInformationForZDirectionOnly();
        //SolveSixDofEquation();
    }
    else if ( RBDMethod == 24  )//单自由度俯仰运动
    {
        ModifyForceAndSixDofInformationForPitchingOnly();
        SolveSixDofEquation();
    }
    else if ( RBDMethod == 25  )//单自由度偏航运动
    {
        TK_Exit::ExceptionExit( "暂无单自由度偏航运动" );
        //ModifyForceAndSixDofInformationForYawingOnly();
        //SolveSixDofEquation();
    }
    else if ( RBDMethod == 26  )//单自由度滚转运动
    {
        ModifyForceAndSixDofInformationForRollingOnly();
        SolveSixDofEquation();
    }
    else if ( RBDMethod == 31  )//两自由度俯仰+滚转运动
    {
        ModifyForceAndSixDofInformationForPitchingAndRollingOnly();
        SolveSixDofEquation();
    }
    else if ( RBDMethod == 101 )//DI_Missile
    {
        RudderRotationByDIControlLaw();
    }
    else if ( RBDMethod == 102 )//DI_089
    {
        RudderRotationByDI089ControlLaw();
    }
    else if ( RBDMethod == 111 )//PID_Atd
    {
        RudderRotationByPIDAttitudeControlLaw();
    }
    else if ( RBDMethod == 112 )//PID_Acc
    {
        RudderRotationByPIDAccelerationControlLaw();
    }
    else if ( RBDMethod == 113 )//PID_089
    {
        RudderRotationByPID089ControlLaw();
    }
    else if ( RBDMethod == 211 )//PID2D_Missile
    {
        RudderRotationByPID2DMissileControlLaw();
    }
    else
    {
    TK_Exit::ExceptionExit( "无此运动方式" );
    }
}

void SixDofSolver::Post()
{
    UpdateUnsteadyTerm();
    DumpSixDofSolverInformation();
}

void SixDofSolver::UpdateUnsteadyTerm()
{
    SixDofParameter * sixDofParameter   = this->GetSixDofParameter( 0 );
    SixDofParameter * sixDofParameterN  = this->GetSixDofParameter( 1 );
    SixDofParameter * sixDofParameterN1 = this->GetSixDofParameter( 2 );

    ( * sixDofParameterN1 ) = ( * sixDofParameterN );
    ( * sixDofParameterN  ) = ( * sixDofParameter  );
}

void SixDofSolver::ComputeSixDofForceTimeSequenceContainer()
{
    int numberOfTimeLayers = this->GetNumberOfTimeLayers();

    for ( int iTimeLayer = 0; iTimeLayer < numberOfTimeLayers; ++ iTimeLayer )
    {
        vector< RDouble > & forceVector = this->GetForceVector( iTimeLayer );
        this->GetDimensionalTotalForceAndMomentOfBody( forceVector, iTimeLayer );
    }
}

void SixDofSolver::GetDimensionalTotalForceAndMomentOfBody( vector< RDouble > & dimensionalGeneralForceContainer, int iTimeLayer )
{
    this->GetDimensionalAerodynamicForceAndMomentOfBody( dimensionalGeneralForceContainer, iTimeLayer );
    this->GetDimensionalAdditionalForceAndMomentOfBody ( dimensionalGeneralForceContainer, iTimeLayer );
}

void SixDofSolver::GetMassCenterAndEulerAngleContainer( vector< RDouble > & massCenterContainer, vector< RDouble > & eulerAngleContainer, int iTimeLayer )
{
    SixDofParameter * sixDofParameter = this->GetSixDofParameter( iTimeLayer );

    sixDofParameter->GetMassCenter( massCenterContainer[ 0 ], massCenterContainer[ 1 ], massCenterContainer[ 2 ] );
    sixDofParameter->GetEulerAngleContainer( eulerAngleContainer );
}

void SixDofSolver::GetMomentAboutMassCenterInInertialFrame( vector< RDouble > & massCenterContainer, vector< RDouble > & totalForceVector, vector< RDouble > & momentVector, vector< RDouble > & momentAboutMassCenterInInertialFrame )
{
    RDouble momentChange [ 3 ];
    PHSPACE::CrossProduct( & massCenterContainer[ 0 ], & totalForceVector[ 0 ], momentChange );

    momentAboutMassCenterInInertialFrame[ 0 ] = momentVector[ 0 ] - momentChange[ 0 ];
    momentAboutMassCenterInInertialFrame[ 1 ] = momentVector[ 1 ] - momentChange[ 1 ];
    momentAboutMassCenterInInertialFrame[ 2 ] = momentVector[ 2 ] - momentChange[ 2 ];
}

//得到物体上的气动力（计算用，力在惯性系，力矩在体轴系）
void SixDofSolver::GetDimensionalAerodynamicForceAndMomentOfBody( vector< RDouble > & dimensionalGeneralAerodynamicForceContainer, int iTimeLayer )
{
    vector< RDouble > massCenterContainer( 3 );
    vector< RDouble > eulerAngleContainer( 3 );

    this->GetMassCenterAndEulerAngleContainer( massCenterContainer, eulerAngleContainer, iTimeLayer );

    vector< RDouble > dimensionalAerodynamicForce ( 3 );
    vector< RDouble > dimensionalAerodynamicMoment( 3 );
    
    GetDimensionalAerodynamicForceAndMoment( dimensionalAerodynamicForce, dimensionalAerodynamicMoment, bodyIndex, iTimeLayer );

    vector< RDouble > dimensionalAerodynamicMomentAboutMassCenter( 3 );
    this->GetMomentAboutMassCenterInInertialFrame( massCenterContainer, dimensionalAerodynamicForce, dimensionalAerodynamicMoment, dimensionalAerodynamicMomentAboutMassCenter );

    //将力矩从惯性系转换到体轴系
    vector< RDouble > dimensionalAerodynamicMomentAboutMassCenterInBodyFrame( 3 );
    PHSPACE::ConvertVectorFromInertialCoordinateSystemToBodyCoordinateSystem( dimensionalAerodynamicMomentAboutMassCenter, eulerAngleContainer, dimensionalAerodynamicMomentAboutMassCenterInBodyFrame );

    dimensionalGeneralAerodynamicForceContainer[ 0 ] = dimensionalAerodynamicForce[ 0 ];
    dimensionalGeneralAerodynamicForceContainer[ 1 ] = dimensionalAerodynamicForce[ 1 ];
    dimensionalGeneralAerodynamicForceContainer[ 2 ] = dimensionalAerodynamicForce[ 2 ];

    dimensionalGeneralAerodynamicForceContainer[ 3 ] = dimensionalAerodynamicMomentAboutMassCenterInBodyFrame[ 0 ];
    dimensionalGeneralAerodynamicForceContainer[ 4 ] = dimensionalAerodynamicMomentAboutMassCenterInBodyFrame[ 1 ];
    dimensionalGeneralAerodynamicForceContainer[ 5 ] = dimensionalAerodynamicMomentAboutMassCenterInBodyFrame[ 2 ];
}

//得到物体上的气动力（惯性系）
void SixDofSolver::GetDimensionalAerodynamicForceAndMomentInerFrame( vector< RDouble > & dimensionalGeneralAerodynamicForceContainer, int iTimeLayer )
{
    vector< RDouble > massCenterContainer( 3 );
    vector< RDouble > eulerAngleContainer( 3 );

    this->GetMassCenterAndEulerAngleContainer( massCenterContainer, eulerAngleContainer, iTimeLayer );

    vector< RDouble > dimensionalAerodynamicForce ( 3 );
    vector< RDouble > dimensionalAerodynamicMoment( 3 );
    
    GetDimensionalAerodynamicForceAndMoment( dimensionalAerodynamicForce, dimensionalAerodynamicMoment, bodyIndex, iTimeLayer );

    vector< RDouble > dimensionalAerodynamicMomentAboutMassCenter( 3 );
    this->GetMomentAboutMassCenterInInertialFrame( massCenterContainer, dimensionalAerodynamicForce, dimensionalAerodynamicMoment, dimensionalAerodynamicMomentAboutMassCenter );
    
    dimensionalGeneralAerodynamicForceContainer[ 0 ] = dimensionalAerodynamicForce[ 0 ];
    dimensionalGeneralAerodynamicForceContainer[ 1 ] = dimensionalAerodynamicForce[ 1 ];
    dimensionalGeneralAerodynamicForceContainer[ 2 ] = dimensionalAerodynamicForce[ 2 ];

    dimensionalGeneralAerodynamicForceContainer[ 3 ] = dimensionalAerodynamicMomentAboutMassCenter[ 0 ];
    dimensionalGeneralAerodynamicForceContainer[ 4 ] = dimensionalAerodynamicMomentAboutMassCenter[ 1 ];
    dimensionalGeneralAerodynamicForceContainer[ 5 ] = dimensionalAerodynamicMomentAboutMassCenter[ 2 ];
}

//得到物体上的气动力（体轴系）
void SixDofSolver::GetDimensionalAerodynamicForceAndMomentBodyFrame( vector< RDouble > & dimensionalGeneralAerodynamicForceContainer, int iTimeLayer )
{
    vector< RDouble > massCenterContainer( 3 );
    vector< RDouble > eulerAngleContainer( 3 );

    this->GetMassCenterAndEulerAngleContainer( massCenterContainer, eulerAngleContainer, iTimeLayer );

    vector< RDouble > dimensionalAerodynamicForce ( 3 );
    vector< RDouble > dimensionalAerodynamicMoment( 3 );
    
    GetDimensionalAerodynamicForceAndMoment( dimensionalAerodynamicForce, dimensionalAerodynamicMoment, bodyIndex, iTimeLayer );

    vector< RDouble > dimensionalAerodynamicMomentAboutMassCenter( 3 );
    this->GetMomentAboutMassCenterInInertialFrame( massCenterContainer, dimensionalAerodynamicForce, dimensionalAerodynamicMoment, dimensionalAerodynamicMomentAboutMassCenter );
    
    //将力、力矩从惯性系转换到体轴系
    vector< RDouble > dimensionalAerodynamicForceInBodyFrame( 3 );
    vector< RDouble > dimensionalAerodynamicMomentAboutMassCenterInBodyFrame( 3 );

    PHSPACE::ConvertVectorFromInertialCoordinateSystemToBodyCoordinateSystem( dimensionalAerodynamicForce, eulerAngleContainer, dimensionalAerodynamicForceInBodyFrame );
    PHSPACE::ConvertVectorFromInertialCoordinateSystemToBodyCoordinateSystem( dimensionalAerodynamicMomentAboutMassCenter, eulerAngleContainer, dimensionalAerodynamicMomentAboutMassCenterInBodyFrame );

    dimensionalGeneralAerodynamicForceContainer[ 0 ] = dimensionalAerodynamicForceInBodyFrame[ 0 ];
    dimensionalGeneralAerodynamicForceContainer[ 1 ] = dimensionalAerodynamicForceInBodyFrame[ 1 ];
    dimensionalGeneralAerodynamicForceContainer[ 2 ] = dimensionalAerodynamicForceInBodyFrame[ 2 ];

    dimensionalGeneralAerodynamicForceContainer[ 3 ] = dimensionalAerodynamicMomentAboutMassCenterInBodyFrame[ 0 ];
    dimensionalGeneralAerodynamicForceContainer[ 4 ] = dimensionalAerodynamicMomentAboutMassCenterInBodyFrame[ 1 ];
    dimensionalGeneralAerodynamicForceContainer[ 5 ] = dimensionalAerodynamicMomentAboutMassCenterInBodyFrame[ 2 ];
}

//附加气动力
void SixDofSolver::GetDimensionalAdditionalForceAndMomentOfBodyByParameter( vector< RDouble > & dimensionalGeneralForceContainer, int iTimeLayer )
{
    //将力从体轴系转换到惯性系；
    vector< RDouble > massCenterContainer( 3 );
    vector< RDouble > eulerAngleContainer( 3 );

    this->GetMassCenterAndEulerAngleContainer( massCenterContainer, eulerAngleContainer, iTimeLayer );

    vector< RDouble > dimensionalAerodynamicForce;
    GetRDoubleVectorFromDataBase( dimensionalAerodynamicForce, bodyIndex, "addedForce"  , 3 );

    vector< RDouble > dimensionalAerodynamicMonent;
    GetRDoubleVectorFromDataBase( dimensionalAerodynamicMonent, bodyIndex, "addedMoment" , 3 );

    vector< RDouble > dimensionalAerodynamicForceInInertialFrame( 3 );

    PHSPACE::ConvertVectorFromBodyCoordinateSystemToInertialCoordinateSystem( dimensionalAerodynamicForce, eulerAngleContainer, dimensionalAerodynamicForceInInertialFrame );

    dimensionalGeneralForceContainer[ 0 ] += dimensionalAerodynamicForceInInertialFrame[ 0 ];
    dimensionalGeneralForceContainer[ 1 ] += dimensionalAerodynamicForceInInertialFrame[ 1 ];
    dimensionalGeneralForceContainer[ 2 ] += dimensionalAerodynamicForceInInertialFrame[ 2 ];
    dimensionalGeneralForceContainer[ 3 ] += dimensionalAerodynamicMonent[ 0 ];
    dimensionalGeneralForceContainer[ 4 ] += dimensionalAerodynamicMonent[ 1 ];
    dimensionalGeneralForceContainer[ 5 ] += dimensionalAerodynamicMonent[ 2 ];
}

//外挂物分离的附加气动力
void SixDofSolver::GetDimensionalAdditionalForceAndMomentOfBodyWingStore( vector< RDouble > & dimensionalGeneralForceContainer, int iTimeLayer )
{
    //力矩在体轴系内，力在惯性系内；
    RDouble dimensionalAleSimulationTime = PHSPACE::GlobalDataBase::GetDoubleParaFromDB("dimensionalAleSimulationTime");

    RDouble mass = GetRDoubleParameterFromDataBase ( bodyIndex, "mass" );

    RDouble gravity = 9.8;

    //fy
    dimensionalGeneralForceContainer[ 1 ] -= gravity * mass;

    //当前的dimensionalAleSimulationTime是n+1时刻的物理时间，此处应当是n时刻的物理时间
    RDouble  timeOfNstep =  dimensionalAleSimulationTime - GetDimensionalAleTimeStep();

    if( timeOfNstep < 0.053 )
    {
        //fy
        dimensionalGeneralForceContainer[ 1 ] -= 53341;
        //mx
        dimensionalGeneralForceContainer[ 3 ] += 62;
        //mz
        dimensionalGeneralForceContainer[ 5 ] -= 12129;

    }
}

//外挂物分离的附加气动力
void SixDofSolver::GetDimensionalAdditionalForceAndMomentOfBody611C7( vector< RDouble > & dimensionalGeneralForceContainer, int iTimeLayer )
{
    //力矩在体轴系内，力在惯性系内；
    RDouble dimensionalAleSimulationTime = PHSPACE::GlobalDataBase::GetDoubleParaFromDB("dimensionalAleSimulationTime");

    RDouble mass = GetRDoubleParameterFromDataBase ( bodyIndex, "mass" );

    RDouble gravity = 9.8;

    //fy
    dimensionalGeneralForceContainer[ 1 ] -= gravity * mass;

    //当前的dimensionalAleSimulationTime是n+1时刻的物理时间，此处应当是n时刻的物理时间
    RDouble  timeOfNstep =  dimensionalAleSimulationTime - GetDimensionalAleTimeStep();

    if( timeOfNstep < 0.07 )
    {
        //fy
        //dimensionalGeneralForceContainer[ 1 ] -= 7478;
        dimensionalGeneralForceContainer[ 1 ] -= 8078;
        //mx
        //dimensionalGeneralForceContainer[ 3 ] += 62;
        //mz
        dimensionalGeneralForceContainer[ 5 ] -= 650.6;

    }
}

void SixDofSolver::GetDimensionalAdditionalForceAndMomentOfBody( vector< RDouble > & dimensionalGeneralForceContainer, int iTimeLayer )
{
    //GetDimensionalAdditionalForceAndMomentOfBodyByParameter( dimensionalGeneralForceContainer, iTimeLayer );

    //下面函数是针对机翼外挂物分离问题特定的
    //GetDimensionalAdditionalForceAndMomentOfBodyWingStore( dimensionalGeneralForceContainer, iTimeLayer );

    //611分离2017
    //GetDimensionalAdditionalForceAndMomentOfBody611C7( dimensionalGeneralForceContainer, iTimeLayer );
}

void SixDofSolver::ReadSixDofParameterFromFile()
{
    string uDFSixDofFileName = GetStringParameterFromDataBase( bodyIndex, "uDFSixDofFileName" );

    fstream infile( uDFSixDofFileName.c_str(), ios_base::in );
    
    string temp;
    int outerIterationSteps = GlobalDataBase::GetIntParaFromDB("outnstep");

    RDouble MassCenterX, MassCenterY, MassCenterZ;
    RDouble EulerAngleAlpha, EulerAngleBeta, EulerAngleGama;

    for ( int iLine = 0; iLine < outerIterationSteps - 1; ++ iLine )
    {
        PHSPACE::ReadNextNonEmptyLine( infile, temp );
    }

    infile >> MassCenterX;
    infile >> MassCenterY;
    infile >> MassCenterZ;

    for ( int i = 0; i < 3; ++ i )
    {
        infile >> temp;
    }

    infile >> EulerAngleAlpha;
    infile >> EulerAngleBeta;
    infile >> EulerAngleGama;

    for ( int i = 0; i < 3; ++ i )
    {
        infile >> temp;
    }

    SixDofParameter * sixDofParameter = this->GetSixDofParameter( 0 );

    sixDofParameter->SetMassCenterX( MassCenterX );
    sixDofParameter->SetMassCenterY( MassCenterY );
    sixDofParameter->SetMassCenterZ( MassCenterZ );
    
    sixDofParameter->SetEulerAngleAlpha( EulerAngleAlpha );
    sixDofParameter->SetEulerAngleBeta ( EulerAngleBeta );
    sixDofParameter->SetEulerAngleGama ( EulerAngleGama );
}

void SixDofSolver::Stationary()
{
    SixDofParameter * sixDofParameter0 = this->GetSixDofParameter( 0 );
    SixDofParameter * sixDofParameter1 = this->GetSixDofParameter( 1 );

    RDouble oldMassCenterX = sixDofParameter1->GetMassCenterX();
    RDouble oldMassCenterY = sixDofParameter1->GetMassCenterY();
    RDouble oldMassCenterZ = sixDofParameter1->GetMassCenterZ();

    RDouble oldEulerAngleAlpha = sixDofParameter1->GetEulerAngleAlpha();
    RDouble oldEulerAngleBeta  = sixDofParameter1->GetEulerAngleBeta();
    RDouble oldEulerAngleGama  = sixDofParameter1->GetEulerAngleGama();

    sixDofParameter0->SetMassCenterX( oldMassCenterX );
    sixDofParameter0->SetMassCenterY( oldMassCenterY );
    sixDofParameter0->SetMassCenterZ( oldMassCenterZ );

    sixDofParameter0->SetEulerAngleAlpha( oldEulerAngleAlpha );
    sixDofParameter0->SetEulerAngleBeta ( oldEulerAngleBeta );
    sixDofParameter0->SetEulerAngleGama ( oldEulerAngleGama );
}

void SixDofSolver::Rolling()
{
    SixDofParameter * sixDofParameter0 = this->GetSixDofParameter( 0 );
    SixDofParameter * sixDofParameter1 = this->GetSixDofParameter( 1 );

    RDouble oldMassCenterX     = sixDofParameter1->GetMassCenterX();
    RDouble oldMassCenterY     = sixDofParameter1->GetMassCenterY();
    RDouble oldMassCenterZ     = sixDofParameter1->GetMassCenterZ();

    RDouble oldEulerAngleAlpha = sixDofParameter1->GetEulerAngleAlpha();
    RDouble oldEulerAngleBeta  = sixDofParameter1->GetEulerAngleBeta();
    RDouble oldEulerAngleGama  = sixDofParameter1->GetEulerAngleGama();

    RDouble newMassCenterX     = sixDofParameter0->GetMassCenterX();
    RDouble newMassCenterY     = sixDofParameter0->GetMassCenterY();
    RDouble newMassCenterZ     = sixDofParameter0->GetMassCenterZ();

    RDouble newEulerAngleAlpha = sixDofParameter0->GetEulerAngleAlpha();
    RDouble newEulerAngleBeta  = sixDofParameter0->GetEulerAngleBeta();
    RDouble newEulerAngleGama  = sixDofParameter0->GetEulerAngleGama();

    RDouble aleAmplitude = PHSPACE::GetRDoubleParameterFromDataBase( bodyIndex, "amplitude" );
    RDouble alePeriod    = PHSPACE::GetRDoubleParameterFromDataBase( bodyIndex, "period"    );

    RDouble dimensionalAleSimulationTime = PHSPACE::GlobalDataBase::GetDoubleParaFromDB("dimensionalAleSimulationTime");
    RDouble dimensionalAleTimeStep       = PHSPACE::GetDimensionalAleTimeStep();

    RDouble deltAlpha = sin( 2.0 * PI * dimensionalAleSimulationTime / alePeriod ) - sin( 2.0 * PI * ( ( dimensionalAleSimulationTime - dimensionalAleTimeStep ) / alePeriod ) );
    
    deltAlpha = aleAmplitude * deltAlpha;

    newMassCenterX     = oldMassCenterX;
    newMassCenterY     = oldMassCenterY;
    newMassCenterZ     = oldMassCenterZ;

    newEulerAngleAlpha = oldEulerAngleAlpha + deltAlpha;
    newEulerAngleBeta  = oldEulerAngleBeta;
    newEulerAngleGama  = oldEulerAngleGama;

    sixDofParameter0->SetMassCenterX( newMassCenterX );
    sixDofParameter0->SetMassCenterY( newMassCenterY );
    sixDofParameter0->SetMassCenterZ( newMassCenterZ );

    sixDofParameter0->SetEulerAngleAlpha( newEulerAngleAlpha );
    sixDofParameter0->SetEulerAngleBeta ( newEulerAngleBeta );
    sixDofParameter0->SetEulerAngleGama ( newEulerAngleGama );
}

void SixDofSolver::Yawing()
{
    SixDofParameter * sixDofParameter0 = this->GetSixDofParameter( 0 );
    SixDofParameter * sixDofParameter1 = this->GetSixDofParameter( 1 );

    RDouble oldMassCenterX     = sixDofParameter1->GetMassCenterX();
    RDouble oldMassCenterY     = sixDofParameter1->GetMassCenterY();
    RDouble oldMassCenterZ     = sixDofParameter1->GetMassCenterZ();

    RDouble oldEulerAngleAlpha = sixDofParameter1->GetEulerAngleAlpha();
    RDouble oldEulerAngleBeta  = sixDofParameter1->GetEulerAngleBeta();
    RDouble oldEulerAngleGama  = sixDofParameter1->GetEulerAngleGama();

    RDouble newMassCenterX     = sixDofParameter0->GetMassCenterX();
    RDouble newMassCenterY     = sixDofParameter0->GetMassCenterY();
    RDouble newMassCenterZ     = sixDofParameter0->GetMassCenterZ();

    RDouble newEulerAngleAlpha = sixDofParameter0->GetEulerAngleAlpha();
    RDouble newEulerAngleBeta  = sixDofParameter0->GetEulerAngleBeta();
    RDouble newEulerAngleGama  = sixDofParameter0->GetEulerAngleGama();

    RDouble dimensionalAleSimulationTime = PHSPACE::GlobalDataBase::GetDoubleParaFromDB("dimensionalAleSimulationTime");
    RDouble dimensionalAleTimeStep       = PHSPACE::GetDimensionalAleTimeStep();

    RDouble aleAmplitude = PHSPACE::GetRDoubleParameterFromDataBase(bodyIndex, "amplitude");
    RDouble alePeriod = PHSPACE::GetRDoubleParameterFromDataBase(bodyIndex, "period");

    RDouble deltBeta = sin(2.0 * PI * dimensionalAleSimulationTime / alePeriod) - sin(2.0 * PI * ((dimensionalAleSimulationTime - dimensionalAleTimeStep) / alePeriod));

    deltBeta = aleAmplitude * deltBeta;

    newMassCenterX     = oldMassCenterX;
    newMassCenterY     = oldMassCenterY;
    newMassCenterZ     = oldMassCenterZ;

    newEulerAngleAlpha = oldEulerAngleAlpha;
    newEulerAngleBeta  = oldEulerAngleBeta + deltBeta;
    newEulerAngleGama  = oldEulerAngleGama;

    sixDofParameter0->SetMassCenterX( newMassCenterX );
    sixDofParameter0->SetMassCenterY( newMassCenterY );
    sixDofParameter0->SetMassCenterZ( newMassCenterZ );

    sixDofParameter0->SetEulerAngleAlpha( newEulerAngleAlpha );
    sixDofParameter0->SetEulerAngleBeta ( newEulerAngleBeta );
    sixDofParameter0->SetEulerAngleGama ( newEulerAngleGama );
}

void SixDofSolver::Pitching()
{
    SixDofParameter* sixDofParameter0 = this->GetSixDofParameter(0);
    SixDofParameter* sixDofParameter1 = this->GetSixDofParameter(1);

    RDouble oldMassCenterX = sixDofParameter1->GetMassCenterX();
    RDouble oldMassCenterY = sixDofParameter1->GetMassCenterY();
    RDouble oldMassCenterZ = sixDofParameter1->GetMassCenterZ();

    RDouble oldEulerAngleAlpha = sixDofParameter1->GetEulerAngleAlpha();
    RDouble oldEulerAngleBeta = sixDofParameter1->GetEulerAngleBeta();
    RDouble oldEulerAngleGama = sixDofParameter1->GetEulerAngleGama();

    RDouble oldAngularVelocityX = sixDofParameter1->GetAngularVelocityX();
    RDouble oldAngularVelocityY = sixDofParameter1->GetAngularVelocityY();
    RDouble oldAngularVelocityZ = sixDofParameter1->GetAngularVelocityZ();

    RDouble newMassCenterX = sixDofParameter0->GetMassCenterX();
    RDouble newMassCenterY = sixDofParameter0->GetMassCenterY();
    RDouble newMassCenterZ = sixDofParameter0->GetMassCenterZ();

    RDouble newEulerAngleAlpha = sixDofParameter0->GetEulerAngleAlpha();
    RDouble newEulerAngleBeta = sixDofParameter0->GetEulerAngleBeta();
    RDouble newEulerAngleGama = sixDofParameter0->GetEulerAngleGama();

    RDouble newAngularVelocityX = sixDofParameter0->GetAngularVelocityX();
    RDouble newAngularVelocityY = sixDofParameter0->GetAngularVelocityY();
    RDouble newAngularVelocityZ = sixDofParameter0->GetAngularVelocityZ();

    RDouble dimensionalAleSimulationTime = PHSPACE::GlobalDataBase::GetDoubleParaFromDB("dimensionalAleSimulationTime");
    RDouble dimensionalAleTimeStep = PHSPACE::GetDimensionalAleTimeStep();

    RDouble aleAmplitude = PHSPACE::GetRDoubleParameterFromDataBase(bodyIndex, "amplitude") * PI / 180.0;
    RDouble aleReduceFrequency = PHSPACE::GetRDoubleParameterFromDataBase(bodyIndex, "reduceFrequency");
    //RDouble alpha = sin(aleReduceFrequency * dimensionalAleSimulationTime);
    //RDouble angularVelocity = aleReduceFrequency * cos(aleReduceFrequency * dimensionalAleSimulationTime);
    //alpha = aleAmplitude * alpha;
    RDouble referenceLength = GlobalDataBase::GetDoubleParaFromDB("referenceLength");
    RDouble referenceVelocity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalVelocity");
    RDouble aleReduceFrequencyDimensional = aleReduceFrequency * referenceVelocity / referenceLength;

    //RDouble deltAlpha = sin(dimensionalAleSimulationTime * aleReduceFrequencyDimensional) - sin((dimensionalAleSimulationTime - dimensionalAleTimeStep) * aleReduceFrequencyDimensional);
    //deltAlpha = aleAmplitude * deltAlpha;
    RDouble alpha = sin(dimensionalAleSimulationTime * aleReduceFrequencyDimensional);
    alpha = aleAmplitude * alpha;
    RDouble angularVelocity = cos(dimensionalAleSimulationTime * aleReduceFrequencyDimensional);
    angularVelocity = aleAmplitude * aleReduceFrequencyDimensional * angularVelocity;

    newMassCenterX = oldMassCenterX;
    newMassCenterY = oldMassCenterY;
    newMassCenterZ = oldMassCenterZ;

    //newEulerAngleAlpha = oldEulerAngleAlpha + deltAlpha;
    newEulerAngleAlpha = alpha;
    newEulerAngleBeta  = oldEulerAngleBeta;
    newEulerAngleGama  = oldEulerAngleGama;

    newAngularVelocityX = oldAngularVelocityX;
    newAngularVelocityX = angularVelocity;
    newAngularVelocityY = oldAngularVelocityY;
    newAngularVelocityZ = oldAngularVelocityZ;

    sixDofParameter0->SetMassCenterX(newMassCenterX);
    sixDofParameter0->SetMassCenterY(newMassCenterY);
    sixDofParameter0->SetMassCenterZ(newMassCenterZ);

    sixDofParameter0->SetEulerAngleAlpha(newEulerAngleAlpha);
    sixDofParameter0->SetEulerAngleBeta(newEulerAngleBeta);
    sixDofParameter0->SetEulerAngleGama(newEulerAngleGama);

    sixDofParameter0->SetAngularVelocityX(newAngularVelocityX);
    sixDofParameter0->SetAngularVelocityY(newAngularVelocityY);
    sixDofParameter0->SetAngularVelocityZ(newAngularVelocityZ);
}

void SixDofSolver::Rotation_Rudder_Test1()
{
    SixDofParameter * sixDofParameter0 = this->GetSixDofParameter( 0 );
    SixDofParameter * sixDofParameter1 = this->GetSixDofParameter( 1 );

    RDouble oldMassCenterX     = sixDofParameter1->GetMassCenterX();
    RDouble oldMassCenterY     = sixDofParameter1->GetMassCenterY();
    RDouble oldMassCenterZ     = sixDofParameter1->GetMassCenterZ();

    RDouble oldEulerAngleAlpha = sixDofParameter1->GetEulerAngleAlpha();
    RDouble oldEulerAngleBeta  = sixDofParameter1->GetEulerAngleBeta ();
    RDouble oldEulerAngleGama  = sixDofParameter1->GetEulerAngleGama ();

    RDouble newMassCenterX     = sixDofParameter0->GetMassCenterX();
    RDouble newMassCenterY     = sixDofParameter0->GetMassCenterY();
    RDouble newMassCenterZ     = sixDofParameter0->GetMassCenterZ();

    RDouble newEulerAngleAlpha = sixDofParameter0->GetEulerAngleAlpha();
    RDouble newEulerAngleBeta  = sixDofParameter0->GetEulerAngleBeta ();
    RDouble newEulerAngleGama  = sixDofParameter0->GetEulerAngleGama ();

    RDouble dimensionalAleSimulationTime = PHSPACE::GlobalDataBase::GetDoubleParaFromDB("dimensionalAleSimulationTime");

    vector< RDouble >                   controlLawParameters;
    GetRDoubleVectorFromDataBase(        controlLawParameters, 1, "controlLawParameters", 2 );

    RDouble maximumRudderAngle        = controlLawParameters[ 0 ] * PI / 180.0;
    RDouble imposedRudderRotationTime = controlLawParameters[ 1 ];

    RDouble rudderAngle = maximumRudderAngle * 0.5 * 
                        ( 1.0 + cos( PI * ( dimensionalAleSimulationTime / imposedRudderRotationTime  - 1.0 ) ) );
    
    if ( dimensionalAleSimulationTime > imposedRudderRotationTime )
    {
        rudderAngle = maximumRudderAngle;
    }

    cout << "rudderAngle = " << rudderAngle * RadToDeg << endl;

    newMassCenterX     = oldMassCenterX;
    newMassCenterY     = oldMassCenterY;
    newMassCenterZ     = oldMassCenterZ;

    newEulerAngleAlpha = rudderAngle;
    newEulerAngleBeta  = oldEulerAngleBeta;
    newEulerAngleGama  = oldEulerAngleGama;

    sixDofParameter0->SetMassCenterX( newMassCenterX );
    sixDofParameter0->SetMassCenterY( newMassCenterY );
    sixDofParameter0->SetMassCenterZ( newMassCenterZ );

    sixDofParameter0->SetEulerAngleAlpha( newEulerAngleAlpha );
    sixDofParameter0->SetEulerAngleBeta ( newEulerAngleBeta  );
    sixDofParameter0->SetEulerAngleGama ( newEulerAngleGama  );
}

void SixDofSolver::OscillatingInXDirection()
{
    SixDofParameter * sixDofParameter0 = this->GetSixDofParameter( 0 );
    SixDofParameter * sixDofParameter1 = this->GetSixDofParameter( 1 );

    RDouble oldMassCenterX     = sixDofParameter1->GetMassCenterX();
    RDouble oldMassCenterY     = sixDofParameter1->GetMassCenterY();
    RDouble oldMassCenterZ     = sixDofParameter1->GetMassCenterZ();

    RDouble oldEulerAngleAlpha = sixDofParameter1->GetEulerAngleAlpha();
    RDouble oldEulerAngleBeta  = sixDofParameter1->GetEulerAngleBeta();
    RDouble oldEulerAngleGama  = sixDofParameter1->GetEulerAngleGama();

    RDouble newMassCenterX     = sixDofParameter0->GetMassCenterX();
    RDouble newMassCenterY     = sixDofParameter0->GetMassCenterY();
    RDouble newMassCenterZ     = sixDofParameter0->GetMassCenterZ();

    RDouble newEulerAngleAlpha = sixDofParameter0->GetEulerAngleAlpha();
    RDouble newEulerAngleBeta  = sixDofParameter0->GetEulerAngleBeta();
    RDouble newEulerAngleGama  = sixDofParameter0->GetEulerAngleGama();

    RDouble aleAmplitude = PHSPACE::GetRDoubleParameterFromDataBase( bodyIndex, "amplitude" );
    RDouble alePeriod    = PHSPACE::GetRDoubleParameterFromDataBase( bodyIndex, "period"    );

    RDouble dimensionalAleSimulationTime = PHSPACE::GlobalDataBase::GetDoubleParaFromDB("dimensionalAleSimulationTime");
    RDouble dimensionalAleTimeStep      = PHSPACE::GetDimensionalAleTimeStep();

    RDouble ds = aleAmplitude * ( sin( 2.0 * PI * dimensionalAleSimulationTime / alePeriod ) - sin( 2.0 * PI * ( ( dimensionalAleSimulationTime - dimensionalAleTimeStep ) / alePeriod ) ) );

    //ds = -1.0 * dimensionalAleTimeStep;

    RDouble dx = ds;
    RDouble dy = 0.0;
    RDouble dz = 0.0;

    newMassCenterX     = oldMassCenterX + dx;
    newMassCenterY     = oldMassCenterY + dy;
    newMassCenterZ     = oldMassCenterZ + dz;

    newEulerAngleAlpha = oldEulerAngleAlpha;
    newEulerAngleBeta  = oldEulerAngleBeta;
    newEulerAngleGama  = oldEulerAngleGama;

    sixDofParameter0->SetMassCenterX( newMassCenterX );
    sixDofParameter0->SetMassCenterY( newMassCenterY );
    sixDofParameter0->SetMassCenterZ( newMassCenterZ );

    sixDofParameter0->SetEulerAngleAlpha( newEulerAngleAlpha );
    sixDofParameter0->SetEulerAngleBeta ( newEulerAngleBeta  );
    sixDofParameter0->SetEulerAngleGama ( newEulerAngleGama  );
}

void SixDofSolver::ModifyForceAndSixDofInformationForXYSymmetryCase()
{
    int numberOfTimeLayers = this->GetNumberOfTimeLayers();

    for ( int iTimeLayer = 0; iTimeLayer < numberOfTimeLayers; ++ iTimeLayer )
    {
        vector< RDouble > & forceVector = this->GetForceVector( iTimeLayer );
        forceVector[ 2 ] = 0;
        forceVector[ 3 ] = 0;
        forceVector[ 4 ] = 0;
    }

    for ( int iTimeLayer = 0; iTimeLayer < numberOfTimeLayers; ++ iTimeLayer )
    {
        SixDofParameter * sixDofParameter = this->GetSixDofParameter( iTimeLayer );

        RDouble massCenterVelocityX, massCenterVelocityY, massCenterVelocityZ;
        RDouble angularVelocityX, angularVelocityY, angularVelocityZ;

        sixDofParameter->GetMassCenterVelocity( massCenterVelocityX, massCenterVelocityY, massCenterVelocityZ );
        sixDofParameter->GetAngularVelocity( angularVelocityX, angularVelocityY, angularVelocityZ );

        massCenterVelocityZ = 0.0;
        angularVelocityX    = 0.0;
        angularVelocityY    = 0.0;

        sixDofParameter->SetMassCenterVelocity( massCenterVelocityX, massCenterVelocityY, massCenterVelocityZ );
        sixDofParameter->SetAngularVelocity( angularVelocityX, angularVelocityY, angularVelocityZ );
    }
}

void SixDofSolver::ModifyForceAndSixDofInformationForPitchingOnly()
{
    int numberOfTimeLayers = this->GetNumberOfTimeLayers();

    for ( int iTimeLayer = 0; iTimeLayer < numberOfTimeLayers; ++ iTimeLayer )
    {
        vector< RDouble > & forceVector = this->GetForceVector( iTimeLayer );
        forceVector[ 0 ] = 0;
        forceVector[ 1 ] = 0;
        forceVector[ 2 ] = 0;
        forceVector[ 3 ] = 0;
        forceVector[ 4 ] = 0;
    }

    for ( int iTimeLayer = 0; iTimeLayer < numberOfTimeLayers; ++ iTimeLayer )
    {
        SixDofParameter * sixDofParameter = this->GetSixDofParameter( iTimeLayer );

        RDouble massCenterVelocityX, massCenterVelocityY, massCenterVelocityZ;
        RDouble angularVelocityX, angularVelocityY, angularVelocityZ;

        sixDofParameter->GetMassCenterVelocity( massCenterVelocityX, massCenterVelocityY, massCenterVelocityZ );
        sixDofParameter->GetAngularVelocity( angularVelocityX, angularVelocityY, angularVelocityZ );
        
        massCenterVelocityX = 0.0;
        massCenterVelocityY = 0.0;
        massCenterVelocityZ = 0.0;
        angularVelocityX    = 0.0;
        angularVelocityY    = 0.0;

        sixDofParameter->SetMassCenterVelocity( massCenterVelocityX, massCenterVelocityY, massCenterVelocityZ );
        sixDofParameter->SetAngularVelocity( angularVelocityX, angularVelocityY, angularVelocityZ );
    }
}

void SixDofSolver::ModifyForceAndSixDofInformationForRollingOnly()
{
    int numberOfTimeLayers = this->GetNumberOfTimeLayers();

    for ( int iTimeLayer = 0; iTimeLayer < numberOfTimeLayers; ++ iTimeLayer )
    {
        vector< RDouble > & forceVector = this->GetForceVector( iTimeLayer );
        forceVector[ 0 ] = 0;
        forceVector[ 1 ] = 0;
        forceVector[ 2 ] = 0;
        forceVector[ 4 ] = 0;
        forceVector[ 5 ] = 0;
    }

    for ( int iTimeLayer = 0; iTimeLayer < numberOfTimeLayers; ++ iTimeLayer )
    {
        SixDofParameter * sixDofParameter = this->GetSixDofParameter( iTimeLayer );

        RDouble massCenterVelocityX, massCenterVelocityY, massCenterVelocityZ;
        RDouble angularVelocityX, angularVelocityY, angularVelocityZ;

        sixDofParameter->GetMassCenterVelocity( massCenterVelocityX, massCenterVelocityY, massCenterVelocityZ );
        sixDofParameter->GetAngularVelocity( angularVelocityX, angularVelocityY, angularVelocityZ );
        
        massCenterVelocityX = 0.0;
        massCenterVelocityY = 0.0;
        massCenterVelocityZ = 0.0;
        angularVelocityY    = 0.0;
        angularVelocityZ    = 0.0;

        sixDofParameter->SetMassCenterVelocity( massCenterVelocityX, massCenterVelocityY, massCenterVelocityZ );
        sixDofParameter->SetAngularVelocity( angularVelocityX, angularVelocityY, angularVelocityZ );
    }
}

void SixDofSolver::ModifyForceAndSixDofInformationForPitchingAndRollingOnly()
{
    int numberOfTimeLayers = this->GetNumberOfTimeLayers();

    for ( int iTimeLayer = 0; iTimeLayer < numberOfTimeLayers; ++ iTimeLayer )
    {
        vector< RDouble > & forceVector = this->GetForceVector( iTimeLayer );
        forceVector[ 0 ] = 0;
        forceVector[ 1 ] = 0;
        forceVector[ 2 ] = 0;
        forceVector[ 4 ] = 0;
    }

    for ( int iTimeLayer = 0; iTimeLayer < numberOfTimeLayers; ++ iTimeLayer )
    {
        SixDofParameter * sixDofParameter = this->GetSixDofParameter( iTimeLayer );

        RDouble massCenterVelocityX, massCenterVelocityY, massCenterVelocityZ;
        RDouble angularVelocityX, angularVelocityY, angularVelocityZ;

        sixDofParameter->GetMassCenterVelocity( massCenterVelocityX, massCenterVelocityY, massCenterVelocityZ );
        sixDofParameter->GetAngularVelocity( angularVelocityX, angularVelocityY, angularVelocityZ );
        
        massCenterVelocityX = 0.0;
        massCenterVelocityY = 0.0;
        massCenterVelocityZ = 0.0;
        angularVelocityY    = 0.0;

        sixDofParameter->SetMassCenterVelocity( massCenterVelocityX, massCenterVelocityY, massCenterVelocityZ );
        sixDofParameter->SetAngularVelocity( angularVelocityX, angularVelocityY, angularVelocityZ );
    }
}

void SixDofSolver::ComputeCurrentSixDofVariables( vector< vector< RDouble > > & qContainer, vector< vector< RDouble > > & fluxContainer, int numberOfEquations )
{
    this->ComputeCurrentSixDofVariables( PHSPACE::HeaderPointer( qContainer   [ 0 ] ), 
                                         PHSPACE::HeaderPointer( qContainer   [ 1 ] ),
                                         PHSPACE::HeaderPointer( qContainer   [ 2 ] ),
                                         PHSPACE::HeaderPointer( fluxContainer[ 0 ] ),
                                         PHSPACE::HeaderPointer( fluxContainer[ 1 ] ), 
                                         PHSPACE::HeaderPointer( fluxContainer[ 2 ] ),
                                         numberOfEquations );
}

void SixDofSolver::ComputeCurrentSixDofVariables( RDouble * q, RDouble * qN, RDouble * qN1, RDouble * fluxM, RDouble * fluxN, RDouble * fluxN1, int numberOfEquations )
{
    RDouble a, b, c, d, e;
    RDouble relaxParameter = GlobalDataBase::GetDoubleParaFromDB( "relaxParameterOfKinetic" );
    PHSPACE::GetSixDofCoefficients( a, b, c, d, e, relaxParameter );

    RDouble coef1 = ( 1 - relaxParameter );
    RDouble coef2 = relaxParameter * a;
    RDouble coef3 = relaxParameter * b;
    RDouble coef4 = relaxParameter * c;
    RDouble coef5 = relaxParameter * d;
    RDouble coef6 = relaxParameter * e;

    for ( int iEquation = 0; iEquation < numberOfEquations; ++ iEquation )
    {
        q[ iEquation ] = coef1 * q     [ iEquation ] + 
                         coef2 * qN    [ iEquation ] + 
                         coef3 * qN1   [ iEquation ] + 
                         coef4 * fluxM [ iEquation ] + 
                         coef5 * fluxN [ iEquation ] + 
                         coef6 * fluxN1[ iEquation ];
    }
}

void SixDofSolver::SolveSixDofEquation()
{
    const int numberOfEquations = 12;
    int numberOfTimeLayers = this->GetNumberOfTimeLayers();

    vector< vector< RDouble > > sixDofVariableContainer;
    vector< vector< RDouble > > sixDofFluxContainer;

    PHSPACE::AllocateVector( sixDofVariableContainer, numberOfTimeLayers, numberOfEquations );
    PHSPACE::AllocateVector( sixDofFluxContainer    , numberOfTimeLayers, numberOfEquations );

    for ( int iTimeLayer = 0; iTimeLayer < numberOfTimeLayers; ++ iTimeLayer )
    {
        vector< RDouble > & forceVector     = this->GetForceVector    ( iTimeLayer );
        SixDofParameter * sixDofParameter = this->GetSixDofParameter( iTimeLayer );

        sixDofParameter->GetGeneralVariableContainer( sixDofVariableContainer[ iTimeLayer ] );

        GetSixDofFlux( sixDofParameter, massCharacter, forceVector, sixDofFluxContainer[ iTimeLayer ] );
    }

    ComputeCurrentSixDofVariables( sixDofVariableContainer, sixDofFluxContainer, numberOfEquations );

    SixDofParameter * sixDofParameter0 = this->GetSixDofParameter( 0 );

    sixDofParameter0->SetGeneralVariableContainer( & sixDofVariableContainer[ 0 ][ 0 ] );
}

//通过文件给定运动轨迹
void SixDofSolver::MovingByFile()
{
    SixDofParameter * sixDofParameter0 = this->GetSixDofParameter( 0 );

    RDouble dimensionalAleSimulationTime = PHSPACE::GlobalDataBase::GetDoubleParaFromDB("dimensionalAleSimulationTime");

    int     outIterationStep = GlobalDataBase::GetIntParaFromDB("outnstep");

    int aleStartStrategy = GlobalDataBase::GetIntParaFromDB("aleStartStrategy");
    if (aleStartStrategy == -1 )
    {
        outIterationStep = 1;
    }

    RDouble newMassCenterX    ;
    RDouble newMassCenterY    ;
    RDouble newMassCenterZ    ;

    RDouble newEulerAngleGama ;
    RDouble newEulerAngleBeta ;
    RDouble newEulerAngleAlpha;

    RDouble tmp;

    string uDFSixDofFileName = GetStringParameterFromDataBase( bodyIndex, "uDFSixDofFileName" );

    ios_base::openmode openMode = ios_base::in;
    fstream file;

    PHSPACE::OpenFile( file, uDFSixDofFileName, openMode );

    for( int iStep = 0; iStep < outIterationStep; ++ iStep )
    {
        file >> tmp;
        file >> newMassCenterX;
        file >> newMassCenterY;
        file >> newMassCenterZ;
        file >> newEulerAngleGama ;
        file >> newEulerAngleBeta ;
        file >> newEulerAngleAlpha;
    }

    newEulerAngleGama  *= PI / 180.0;
    newEulerAngleBeta  *= PI / 180.0;
    newEulerAngleAlpha *= PI / 180.0;

    sixDofParameter0->SetMassCenterX( newMassCenterX );
    sixDofParameter0->SetMassCenterY( newMassCenterY );
    sixDofParameter0->SetMassCenterZ( newMassCenterZ );

    sixDofParameter0->SetEulerAngleAlpha( newEulerAngleAlpha );
    sixDofParameter0->SetEulerAngleBeta ( newEulerAngleBeta  );
    sixDofParameter0->SetEulerAngleGama ( newEulerAngleGama  );
}

void SixDofSolver::GetDimensionalAerodynamicForceAndMoment( vector< RDouble > & dimensionalAerodynamicForce, vector< RDouble > & dimensionalAerodynamicMoment, int bodyIndex, int iTimeLayer )
{
    AleForceManager * aleForceManager = GetAleManager()->GetAleForceManager();

    BasicAleForce * basicAleForce = aleForceManager->GetBasicAleForce( bodyIndex, iTimeLayer );

    RDouble dimensionalForceCoefficient  = aleForceManager->GetDimensionalForceCoefficient ();
    RDouble dimensionalMomentCoefficient = aleForceManager->GetDimensionalMomentCoefficient();

    dimensionalAerodynamicForce [ 0 ] = basicAleForce->GetForceX()  * dimensionalForceCoefficient;
    dimensionalAerodynamicForce [ 1 ] = basicAleForce->GetForceY()  * dimensionalForceCoefficient;
    dimensionalAerodynamicForce [ 2 ] = basicAleForce->GetForceZ()  * dimensionalForceCoefficient;

    dimensionalAerodynamicMoment[ 0 ] = basicAleForce->GetMomentX() * dimensionalMomentCoefficient;
    dimensionalAerodynamicMoment[ 1 ] = basicAleForce->GetMomentY() * dimensionalMomentCoefficient;
    dimensionalAerodynamicMoment[ 2 ] = basicAleForce->GetMomentZ() * dimensionalMomentCoefficient;
}

void SixDofSolver::DumpSixDofSolverInformation()
{
    this->DumpSixDofParameterInformation();
    this->DumpAleForceInformation();
}

void SixDofSolver::DumpSixDofParameterInformation()
{
    string fileName = GlobalDataBase::GetStrParaFromDB("sixDofFileName");
    string sixDofFileName = PHSPACE::AddSymbolToFileName(fileName, bodyIndex + 1 );

    ios_base::openmode openMode = ios_base::out|ios_base::app;

    fstream file;
    PHSPACE::OpenFile( file, sixDofFileName, openMode );

    this->DumpSixDofParameterInformation( file );

    PHSPACE::CloseFile( file );
}

void SixDofSolver::DumpAleForceInformation()
{
    string fileName = GlobalDataBase::GetStrParaFromDB("sixDofFileName");
    string aleForceFileName = PHSPACE::AddSymbolToFileName(fileName, bodyIndex + 1 );

    ios_base::openmode openMode = ios_base::out|ios_base::app;
    fstream file;

    string fileName1 = PHSPACE::AddSymbolToFileName( aleForceFileName, "_body" );
    PHSPACE::OpenFile( file, fileName1, openMode );
    this->DumpAleForceInformationBodyFrame( file );
    PHSPACE::CloseFile( file );

    string fileName2 = PHSPACE::AddSymbolToFileName( aleForceFileName, "_iner" );
    PHSPACE::OpenFile( file, fileName2, openMode );
    this->DumpAleForceInformationInerFrame( file );
    PHSPACE::CloseFile( file );
}

//输出自由度信息
void SixDofSolver::DumpSixDofParameterInformation( fstream & file )
{
    RDouble massCenterX, massCenterY, massCenterZ;
    RDouble massCenterVelocityX, massCenterVelocityY, massCenterVelocityZ;
    RDouble eulerAngleAlpha, eulerAngleBeta, eulerAngleGama;
    RDouble angularVelocityX, angularVelocityY, angularVelocityZ;

    //string header = "# time		massCenterX, massCenterY, massCenterZ,  massCenterVelocityX		massCenterVelocityY      massCenterVelocityZ    eulerAngleGama       eulerAngleBeta       eulerAngleAlpha     angularVelocityX      angularVelocityY      angularVelocityZ";

    int iTimeLayer = 0;
    SixDofParameter * sixDofParameter =  this->GetSixDofParameter( iTimeLayer );

    sixDofParameter->GetMassCenter( massCenterX, massCenterY, massCenterZ );
    sixDofParameter->GetMassCenterVelocity( massCenterVelocityX, massCenterVelocityY, massCenterVelocityZ );
    sixDofParameter->GetEulerAngle( eulerAngleGama, eulerAngleBeta,eulerAngleAlpha );
    sixDofParameter->GetAngularVelocity( angularVelocityX, angularVelocityY, angularVelocityZ );

    int wordWidth = 16;
    file << setiosflags( ios::left );
    file << setiosflags( ios::scientific );
    file << setprecision( 6 );

    RDouble dimensionalAleSimulationTime = PHSPACE::GlobalDataBase::GetDoubleParaFromDB("dimensionalAleSimulationTime");

    file << dimensionalAleSimulationTime << "	";

    file << massCenterX << "	";
    file << massCenterY << "	";
    file << massCenterZ << "	";
    file << massCenterVelocityX << "	";
    file << massCenterVelocityY << "	";
    file << massCenterVelocityZ << "	";
    file << eulerAngleGama << "	";
    file << eulerAngleBeta << "	";
    file << eulerAngleAlpha << "	";
    file << angularVelocityX << "	";
    file << angularVelocityY << "	";
    file << angularVelocityZ << "	";

    file << "\n";
}

void SixDofSolver::DumpAleForceInformationBodyFrame( fstream & file )
{
    int numberOfAleForceComponents = 6;
    vector< RDouble > aleForceContainer( numberOfAleForceComponents );

    this->GetDimensionalAerodynamicForceAndMomentBodyFrame( aleForceContainer, 0 );

    int wordWidth = 16;
    file << setiosflags( ios::left );
    file << setiosflags( ios::scientific );
    file << setprecision( 6 );

    RDouble dimensionalAleSimulationTime = PHSPACE::GlobalDataBase::GetDoubleParaFromDB("dimensionalAleSimulationTime");

    file << dimensionalAleSimulationTime << "\t";

    for ( int iAleForceComponent = 0; iAleForceComponent < numberOfAleForceComponents; ++ iAleForceComponent )
    {
        file << aleForceContainer[ iAleForceComponent ] << "\t";
    }
    file << "\n";
}

void SixDofSolver::DumpAleForceInformationInerFrame( fstream & file )
{
    int numberOfAleForceComponents = 6;
    vector< RDouble > aleForceContainer( numberOfAleForceComponents );

    this->GetDimensionalAerodynamicForceAndMomentInerFrame( aleForceContainer, 0 );

    int wordWidth = 16;
    file << setiosflags( ios::left );
    file << setiosflags( ios::scientific );
    file << setprecision( 6 );

    RDouble dimensionalAleSimulationTime = PHSPACE::GlobalDataBase::GetDoubleParaFromDB("dimensionalAleSimulationTime");

    file << dimensionalAleSimulationTime << "\t";

    for ( int iAleForceComponent = 0; iAleForceComponent < numberOfAleForceComponents; ++ iAleForceComponent )
    {
        file << aleForceContainer[ iAleForceComponent ] << "\t";
    }
    file << "\n";
}

void SixDofSolver::DumpRestartSixDofSolver( fstream & file )
{
    const int numberOfEquations = 12;
    vector< RDouble >  generalVariableContainer( numberOfEquations );

    for ( int iTimeLayer = 0; iTimeLayer < numberOfTimeLayers; ++ iTimeLayer )
    {
        SixDofParameter * sixDofParameter = this->GetSixDofParameter( iTimeLayer );
        sixDofParameter->GetGeneralVariableContainer( generalVariableContainer );
        PHSPACE::PHWrite( file, generalVariableContainer, numberOfEquations );
    }
}
//
//void FYControlLaw::DumpRestartSixDofSolver( fstream & file )
//{
//	SixDofSolver::DumpRestartSixDofSolver( file );
//	PHSPACE::PHWrite( file, integralValue );
//}
//
//void FY2DControlLaw::DumpRestartSixDofSolver( fstream & file )
//{
//	SixDofSolver::DumpRestartSixDofSolver( file );
//	PHSPACE::PHWrite( file, integralValueA1  );
//	PHSPACE::PHWrite( file, integralValueA2  );
//	PHSPACE::PHWrite( file, integralValueB   );
//	PHSPACE::PHWrite( file, oldRudderAngleA1 );
//	PHSPACE::PHWrite( file, oldRudderAngleA2 );
//	PHSPACE::PHWrite( file, oldRudderAngleB  );
//}

void GetSixDofCoefficients( RDouble & a, RDouble & b, RDouble & c, RDouble & d, RDouble & e, RDouble & relaxParameter )
{
    int methodForKineticEquation = PHSPACE::GlobalDataBase::GetIntParaFromDB( "methodForKineticEquation" );
    relaxParameter = PHSPACE::GlobalDataBase::GetDoubleParaFromDB( "relaxParameterOfKinetic" );

    RDouble dimensionalAleTimeStep  = PHSPACE::GetDimensionalAleTimeStep();

    RDouble thet, xi, phi;

    if (      methodForKineticEquation == KINETIC_SCHEME::ADMAS_BASHFORTH_FIRST_ORDER  )
    {
        thet = 0.0;
        xi   = 0.0;
        phi  = 0.0;
        relaxParameter = 1.0;
    }
    else if ( methodForKineticEquation == KINETIC_SCHEME::ADMAS_BASHFORTH_SECOND_ORDER )
    {
        thet = 0.0;
        xi   = 0.0;
        phi  = 1.0 / 2.0;
        relaxParameter = 1.0;
    }
    else if ( methodForKineticEquation == KINETIC_SCHEME::IMPLICIT_EULER_FIRST_ORDER  )
    {
        thet = 1.0;
        xi   = 0.0;
        phi  = 0.0;
    }
    else if ( methodForKineticEquation == KINETIC_SCHEME::IMPLICIT_EULER_SECOND_ORDER  )
    {
        thet = 1.0;
        xi   = 1.0 / 2.0;
        phi  = 0.0;
    }
    else if ( methodForKineticEquation == KINETIC_SCHEME::ADMAS_MOULTON_SECOND_ORDER  )
    {
        thet = 1.0 / 2.0;
        xi   = 0.0;
        phi  = 0.0;
    }
    else if ( methodForKineticEquation == KINETIC_SCHEME::ADMAS_MOULTON_THIRD_ORDER )
    {
        thet = 5.0 / 12.0;
        xi   = 0.0;
        phi  = 1.0 / 12.0;
    }

    a = dimensionalAleTimeStep / ( 1.0 + xi ) * ( ( 1.0 + xi ) / dimensionalAleTimeStep + xi / dimensionalAleTimeStep );
    b = dimensionalAleTimeStep / ( 1.0 + xi ) * ( - xi / dimensionalAleTimeStep );
    c = dimensionalAleTimeStep / ( 1.0 + xi ) * ( thet );
    d = dimensionalAleTimeStep / ( 1.0 + xi ) * ( 1.0 - thet + phi );
    e = dimensionalAleTimeStep / ( 1.0 + xi ) * ( - phi );
}

//cxh:右端项（关键！）
void SixDofSolver::GetSixDofFlux( SixDofParameter * sixDofParameter, MassCharacter * massCharacter, vector< RDouble > & force, vector< RDouble > & generalSixDofFlux )
{
    this->ComputeRightHandSideOfRigidBodyTranslationalMotionEquation  ( sixDofParameter, massCharacter, force, generalSixDofFlux );
    this->ComputeRightHandSideOfRigidBodyTranslationalDynamicsEquation( sixDofParameter, massCharacter, force, generalSixDofFlux );
    this->ComputeRightHandSideOfRigidBodyRotationalMotionEquation     ( sixDofParameter, massCharacter, force, generalSixDofFlux );
    this->ComputeRightHandSideOfRigidBodyRotationalDynamicsEquation   ( sixDofParameter, massCharacter, force, generalSixDofFlux );
}

void SixDofSolver::ComputeRightHandSideOfRigidBodyTranslationalMotionEquation( SixDofParameter * sixDofParameter, MassCharacter * massCharacter, vector< RDouble > & force, vector< RDouble > & generalSixDofFlux )
{
    RDouble massCenterVelocityX,  massCenterVelocityY, massCenterVelocityZ;
    sixDofParameter->GetMassCenterVelocity( massCenterVelocityX, massCenterVelocityY, massCenterVelocityZ );

    //location
    generalSixDofFlux[ 0 ] = massCenterVelocityX;
    generalSixDofFlux[ 1 ] = massCenterVelocityY;
    generalSixDofFlux[ 2 ] = massCenterVelocityZ;
}

void SixDofSolver::ComputeRightHandSideOfRigidBodyTranslationalDynamicsEquation( SixDofParameter * sixDofParameter, MassCharacter * massCharacter, vector< RDouble > & force, vector< RDouble > & generalSixDofFlux )
{
    RDouble forceX = force[ 0 ];
    RDouble forceY = force[ 1 ];
    RDouble forcez = force[ 2 ];

    RDouble mass = massCharacter->GetMass();

    //velocity
    generalSixDofFlux[ 3 ] = forceX / mass;
    generalSixDofFlux[ 4 ] = forceY / mass;
    generalSixDofFlux[ 5 ] = forcez / mass;
}

void SixDofSolver::ComputeRightHandSideOfRigidBodyRotationalMotionEquation( SixDofParameter * sixDofParameter, MassCharacter * massCharacter, vector< RDouble > & force, vector< RDouble > & generalSixDofFlux )
{
    RDouble angularVelocityX, angularVelocityY, angularVelocityZ;
    RDouble eulerAngleAlpha, eulerAngleBeta, eulerAngleGama;

    sixDofParameter->GetEulerAngle( eulerAngleGama, eulerAngleBeta, eulerAngleAlpha );
    sixDofParameter->GetAngularVelocity( angularVelocityX, angularVelocityY, angularVelocityZ );

    //attitude
    generalSixDofFlux[ 8 ] = angularVelocityY * sin( eulerAngleGama  )                          + angularVelocityZ * cos( eulerAngleGama  );
    generalSixDofFlux[ 7 ] = angularVelocityY * cos( eulerAngleGama  ) / cos( eulerAngleAlpha ) - angularVelocityZ * sin( eulerAngleGama  ) / cos( eulerAngleAlpha );
    generalSixDofFlux[ 6 ] = angularVelocityX - generalSixDofFlux[ 7 ] * sin( eulerAngleAlpha );
}

void SixDofSolver::ComputeRightHandSideOfRigidBodyRotationalDynamicsEquation( SixDofParameter * sixDofParameter, MassCharacter * massCharacter, vector< RDouble > & force, vector< RDouble > & generalSixDofFlux )
{
    RDouble ** momentOfInertiaTensor        = massCharacter->GetMomentOfInertiaTensor();
    RDouble ** inverseMomentOfInertiaTensor = massCharacter->GetInverseMomentOfInertiaTensor();

    //rotation dRotation/dt

    vector< RDouble > angularVelocityVector( 3 );
    sixDofParameter->GetAngularVelocityContainer( angularVelocityVector );

    vector< RDouble > momentumMoment( 3 );
    PHSPACE::MatrixMultiply( momentOfInertiaTensor, angularVelocityVector, momentumMoment, 3, 3 );

    RDouble momentVector[ 3 ];
    momentVector[ 0 ] = force[ 3 ];
    momentVector[ 1 ] = force[ 4 ];
    momentVector[ 2 ] = force[ 5 ];

    vector< RDouble > temp1( 3 );
    PHSPACE::CrossProduct( angularVelocityVector, momentumMoment, temp1 );

    vector< RDouble > temp2( 3 );
    temp2[ 0 ] = momentVector[ 0 ] - temp1[ 0 ];
    temp2[ 1 ] = momentVector[ 1 ] - temp1[ 1 ];
    temp2[ 2 ] = momentVector[ 2 ] - temp1[ 2 ];

    vector< RDouble > dRotationDtVector( 3 );

    PHSPACE::MatrixMultiply( inverseMomentOfInertiaTensor, temp2, dRotationDtVector, 3, 3 );

    generalSixDofFlux[ 9  ] = dRotationDtVector[ 0 ];
    generalSixDofFlux[ 10 ] = dRotationDtVector[ 1 ];
    generalSixDofFlux[ 11 ] = dRotationDtVector[ 2 ];

    PHSPACE::DelPointer2( momentOfInertiaTensor        );
    PHSPACE::DelPointer2( inverseMomentOfInertiaTensor );
}

SixDofSolverManager::SixDofSolverManager()
{
    sixDofSolverContainer = 0;
    numberOfMovingBodies  = 0;
}

SixDofSolverManager::~SixDofSolverManager() 
{ 
    if ( sixDofSolverContainer )
    {
        for ( int iMovingBody = 0; iMovingBody < numberOfMovingBodies; ++ iMovingBody )
        {
            delete sixDofSolverContainer[ iMovingBody ];
            sixDofSolverContainer[iMovingBody] = NULL;
        }
    }
    delete [] sixDofSolverContainer;
    sixDofSolverContainer = NULL;
}

void SixDofSolverManager::AllocateSixDofSolverContainer()
{
    numberOfMovingBodies = PHSPACE::GlobalDataBase::GetIntParaFromDB( "numberOfMovingBodies" );

    sixDofSolverContainer = new SixDofSolver * [ numberOfMovingBodies ];
    for ( int iMovingBody = 0; iMovingBody < numberOfMovingBodies; ++ iMovingBody )
    {
        int RBDMethod = PHSPACE::GetIntegerParameterFromDataBase( iMovingBody, "RBDMethod" );
        
        if ( RBDMethod <= 100 )
        {
            sixDofSolverContainer[ iMovingBody ] = new SixDofSolver( this );
        }
        //else if ( RBDMethod == 101 || RBDMethod == 102 )
        //{
        //	sixDofSolverContainer[ iMovingBody ] = new FYDIControlLaw( this );
        //}
        //else if ( RBDMethod == 111 || RBDMethod == 112 || RBDMethod == 113 )
        //{
        //	sixDofSolverContainer[ iMovingBody ] = new FYPIDControlLaw( this );
        //}
        //else if ( RBDMethod == 211 )
        //{
        //	sixDofSolverContainer[ iMovingBody ] = new FYPID2DControlLaw( this );
        //}
        else
        {
            sixDofSolverContainer[ iMovingBody ] = new SixDofSolver( this );
        }
    }
}

void SixDofSolverManager::Initialize()
{
    AllocateSixDofSolverContainer();

    for ( int iMovingBody = 0; iMovingBody < numberOfMovingBodies; ++ iMovingBody )
    {
        sixDofSolverContainer[ iMovingBody ]->Initialize( iMovingBody );
    }
}

void SixDofSolverManager::Restart()
{
    for ( int iMovingBody = 0; iMovingBody < numberOfMovingBodies; ++ iMovingBody )
    {
        sixDofSolverContainer[ iMovingBody ]->Restart();
    }
}

void SixDofSolverManager::DumpRestart( fstream & file )
{
    for ( int iMovingBody = 0; iMovingBody < numberOfMovingBodies; ++ iMovingBody )
    {
        sixDofSolverContainer[ iMovingBody ]->DumpRestartSixDofSolver( file );
    }
}

void SixDofSolverManager::Run()
{
    using namespace PHMPI;
    int currentProcessor = GetCurrentProcessorID();
    int serverProcessor = GetServerProcessorID();

    vector< RDouble > sixDofVariableContainer;

    sixDofVariableContainer.resize( 12 );

    for ( int iMovingBody = 0; iMovingBody < numberOfMovingBodies; ++ iMovingBody )
    {
        SixDofParameter * sixDofParameter   = sixDofSolverContainer[ iMovingBody ]->GetSixDofParameter( 0 );

        if( currentProcessor == serverProcessor )
        {
            sixDofSolverContainer[ iMovingBody ]->Run();
            sixDofParameter->GetGeneralVariableContainer( sixDofVariableContainer );
        }

        PHMPI::PH_Bcast( &sixDofVariableContainer[ 0 ], 12*sizeof(RDouble), serverProcessor );

        sixDofParameter->SetGeneralVariableContainer( &sixDofVariableContainer[ 0 ] );
        //WriteLogFile("sixDofVariableContainer 8 = ", sixDofVariableContainer[8], true);
    }
}

void SixDofSolverManager::Post()
{
    for ( int iMovingBody = 0; iMovingBody < numberOfMovingBodies; ++ iMovingBody )
    {
        SixDofSolver * sixDofSolver = this->GetSixDofSolver( iMovingBody );
        sixDofSolver->UpdateUnsteadyTerm();
    }
    
    Dump();
}

void SixDofSolverManager::Dump()
{
    using namespace PHMPI;
    int currentProcessor = GetCurrentProcessorID();
    int serverProcessor = GetServerProcessorID();
    if(currentProcessor != serverProcessor)
    {
        return;
    }

    for ( int iMovingBody = 0; iMovingBody < numberOfMovingBodies; ++ iMovingBody )
    {
        SixDofSolver * sixDofSolver = this->GetSixDofSolver( iMovingBody );
        sixDofSolver->DumpSixDofSolverInformation();
    }
}

MassCharacter * SixDofSolverManager::GetMassCharacter( int bodyIndex )
{
    SixDofSolver * sixDofSolver = this->GetSixDofSolver( bodyIndex );
    return sixDofSolver->GetMassCharacter();
}

SixDofParameter * SixDofSolverManager::GetSixDofParameterN1( int bodyIndex )
{
    SixDofSolver * sixDofSolver = this->GetSixDofSolver( bodyIndex );
    return sixDofSolver->GetSixDofParameter( 2 );
}

SixDofParameter * SixDofSolverManager::GetSixDofParameterN ( int bodyIndex )
{
    SixDofSolver * sixDofSolver = this->GetSixDofSolver( bodyIndex );
    return sixDofSolver->GetSixDofParameter( 1 );
}

SixDofParameter * SixDofSolverManager::GetSixDofParameter  ( int bodyIndex )
{
    SixDofSolver * sixDofSolver = this->GetSixDofSolver( bodyIndex );
    return sixDofSolver->GetSixDofParameter( 0 );
}

SixDofParameter * GetSixDofParameter( int bodyIndex )
{
    SixDofSolverManager * sixDofSolverManager = GetAleManager()->GetSixDofSolverManager();

    return sixDofSolverManager->GetSixDofParameter( bodyIndex );
}

//计算系到体轴系
void ConvertVectorFromInertialCoordinateSystemToBodyCoordinateSystem( vector< RDouble > & vectorInInertialFrame, vector< RDouble > & angleVector, vector< RDouble > & massCenter, vector< RDouble > & vectorInBodyFrame )
{
    vectorInBodyFrame[ 0 ] = vectorInInertialFrame[ 0 ] - massCenter[ 0 ];
    vectorInBodyFrame[ 1 ] = vectorInInertialFrame[ 1 ] - massCenter[ 1 ];
    vectorInBodyFrame[ 2 ] = vectorInInertialFrame[ 2 ] - massCenter[ 2 ];

    PHSPACE::RotateAboutAxisY( vectorInBodyFrame, angleVector[ 1 ], vectorInBodyFrame );
    PHSPACE::RotateAboutAxisZ( vectorInBodyFrame, angleVector[ 2 ], vectorInBodyFrame );
    PHSPACE::RotateAboutAxisX( vectorInBodyFrame, angleVector[ 0 ], vectorInBodyFrame );
}

//计算系到体轴系
void ConvertVectorFromInertialCoordinateSystemToBodyCoordinateSystem( vector< RDouble > & vectorInInertialFrame, vector< RDouble > & angleVector, vector< RDouble > & vectorInBodyFrame )
{
    PHSPACE::RotateAboutAxisY( vectorInInertialFrame, angleVector[ 1 ], vectorInBodyFrame );
    PHSPACE::RotateAboutAxisZ( vectorInBodyFrame    , angleVector[ 2 ], vectorInBodyFrame );
    PHSPACE::RotateAboutAxisX( vectorInBodyFrame    , angleVector[ 0 ], vectorInBodyFrame );
}

//体轴系到计算系
void ConvertVectorFromBodyCoordinateSystemToInertialCoordinateSystem( vector< RDouble > & vectorInBodyFrame, vector< RDouble > & angleVector, vector< RDouble > & massCenter,  vector< RDouble > & vectorInInertialFrame )
{
    PHSPACE::RotateAboutAxisX( vectorInBodyFrame    , - angleVector[ 0 ], vectorInInertialFrame );
    PHSPACE::RotateAboutAxisZ( vectorInInertialFrame, - angleVector[ 2 ], vectorInInertialFrame );
    PHSPACE::RotateAboutAxisY( vectorInInertialFrame, - angleVector[ 1 ], vectorInInertialFrame );

    vectorInInertialFrame[ 0 ] += massCenter[ 0 ];
    vectorInInertialFrame[ 1 ] += massCenter[ 1 ];
    vectorInInertialFrame[ 2 ] += massCenter[ 2 ];
}

//体轴系到计算系
void ConvertVectorFromBodyCoordinateSystemToInertialCoordinateSystem( vector< RDouble > & vectorInBodyFrame, vector< RDouble > & angleVector, vector< RDouble > & vectorInInertialFrame )
{
    PHSPACE::RotateAboutAxisX( vectorInBodyFrame    , - angleVector[ 0 ], vectorInInertialFrame );
    PHSPACE::RotateAboutAxisZ( vectorInInertialFrame, - angleVector[ 2 ], vectorInInertialFrame );
    PHSPACE::RotateAboutAxisY( vectorInInertialFrame, - angleVector[ 1 ], vectorInInertialFrame );
}

//绕x轴旋转
void RotateAboutAxisX( vector< RDouble > & vectorIn, RDouble rotateAngle, vector< RDouble > & vectorOut )
{
    RDouble oldX = vectorIn[ 0 ];
    RDouble oldY = vectorIn[ 1 ];
    RDouble oldZ = vectorIn[ 2 ];
    vectorOut[ 0 ] = 1.0      * oldX + 0.0                * oldY + 0.0                * oldZ;
    vectorOut[ 1 ] = 0.0      * oldX + cos( rotateAngle ) * oldY + sin( rotateAngle ) * oldZ;
    vectorOut[ 2 ] = 0.0      * oldX - sin( rotateAngle ) * oldY + cos( rotateAngle ) * oldZ;
}

//绕y轴旋转
void RotateAboutAxisY( vector< RDouble > & vectorIn, RDouble rotateAngle, vector< RDouble > & vectorOut )
{
    RDouble oldX = vectorIn[ 0 ];
    RDouble oldY = vectorIn[ 1 ];
    RDouble oldZ = vectorIn[ 2 ];
    vectorOut[ 0 ] = cos( rotateAngle ) * oldX + 0.0      * oldY - sin( rotateAngle ) * oldZ;
    vectorOut[ 1 ] = 0.0                * oldX + 1.0      * oldY + 0.0                * oldZ;
    vectorOut[ 2 ] = sin( rotateAngle ) * oldX + 0.0      * oldY + cos( rotateAngle ) * oldZ;
}

//绕z轴旋转
void RotateAboutAxisZ( vector< RDouble > & vectorIn, RDouble rotateAngle, vector< RDouble > & vectorOut )
{
    RDouble oldX = vectorIn[ 0 ];
    RDouble oldY = vectorIn[ 1 ];
    RDouble oldZ = vectorIn[ 2 ];
    vectorOut[ 0 ] =   cos( rotateAngle ) * oldX + sin( rotateAngle ) * oldY + 0.0      * oldZ;
    vectorOut[ 1 ] = - sin( rotateAngle ) * oldX + cos( rotateAngle ) * oldY + 0.0      * oldZ;
    vectorOut[ 2 ] =   0.0                * oldX + 0.0                * oldY + 1.0      * oldZ;
}

}