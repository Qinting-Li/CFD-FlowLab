#pragma once
#include "Constants.h"
#include "Geo_StructGrid.h"
#include "PHSuperMatrix.h"
#include "Force.h"

namespace PHSPACE
{

class PHRigidBodyMovement
{
public:
	typedef PHMatrix< RDouble > PHSuperMatrix;
protected:
    //变量声明
	//飞行器属性，暂时没考虑变质量和变惯量问题；
	RDouble massOfAircraft;                           //  飞行器质量
	PHSuperMatrix * massCenterOfAircraft;              //  飞行器质心位置
    PHSuperMatrix * inertiaTensorOfAircraft;           //  飞行器转动惯量
	PHSuperMatrix * reverseOfInertiaTensor;            //  转动惯量逆矩阵

	// 与流场解算器进行交换的数据；
	// 输入：力和力矩系数
	PHSuperMatrix * aerodynamicForce;                  //  气动力
	PHSuperMatrix * aerodynamicMoment;                 //  气动力矩
	PHSuperMatrix * otherForceSum;                     //  其他外力合力（重力、推力）
	PHSuperMatrix * otherMomentSum;                    //  其他外力合力矩
	// 输入：流场解算器的物理推进无量纲时间步长
	RDouble currentMarchTimeStep;                     //  推进的时间步长
	// 输入：来流攻角、侧滑角
	RDouble initiativeAngleOfAttack;                  //  初始攻角
	RDouble initiativeAngleOfSlide;                   //  初始侧滑角
	// 输出：瞬时攻角、侧滑角
	RDouble temporaryAngleOfAttack;                   //  瞬时攻角
	RDouble temporaryAngleOfSlide;                    //  瞬时侧滑角
	// 输入/输出：网格信息
    StructGrid* grid;
	//	int ni, nj, nk;
	//坐标变换数据；
	PHSuperMatrix * eulerAngle;                        //  欧拉角
	PHSuperMatrix * eulerAngleRate;                    //  n-2, n-1, n时刻欧拉角速度
	PHSuperMatrix * eulerAngleAcceleration;            //  欧拉角加速度

	PHSuperMatrix * transformMatrixCurrent;            //  当前步计算系 -> 体轴系的坐标变换矩阵
	PHSuperMatrix * transformMatrixLastTime;               //  前一步计算系 -> 体轴系的坐标变换矩阵
	//运动量；
	PHSuperMatrix * momentOfInertia;                   //  n-2, n-1, n时刻体轴系下的角动量
	PHSuperMatrix * bodyAngleRate;                     //  体轴系（计算坐标系）角速度
	PHSuperMatrix * bodyAngleAcceleration;             //  体轴系（计算坐标系）角加速度

	//参考量；
	//输入：力和力矩系数
	RDouble densityRef;                               //  参考密度
	RDouble areaRef;                                  //  参考面积
	RDouble lengthRef;                                //  参考长度
public:
	PHRigidBodyMovement();
	~PHRigidBodyMovement();
protected:
	void UpdateCoordinate               ( Grid * gridIn );
	void UpdateGridVelocity             ( Grid * gridIn );
	void InitializeAircraftProperty     ();
	void InitializeEulerAngleInformation();
	void InitializeInertia      ();
	void CalculationOfMatrix    ();
	void NondimensionalOf6DOF   ();
	void GetReferenceValue      ();
	void AerodynamicsInCFDSolver();
	void UpdateDynamicInfo6DOF();
	void UpdateMomentOfInertia();
	void UpdateOfInertia();
	void CalculateAngleRateOfBody();
	void CalculateEulerAngleRate();
	void CalculateEulerAngle();
	void CalculateEulerAngleAcceleration();
	void CalculateAngleOfAttackAndSlide();
	PHSuperMatrix CalculateGridCellFaceVelocity(Grid* gridIn, int i, int j, int k, int nd );
	PHSuperMatrix GetCellFaceCoordinate(Grid* gridIn, int i, int j, int k, int nd );
	void UpdateCellArea(Grid* gridIn);
	void UpdateSingleZoneGridInformation( Grid * gridIn );
public:
	void InitializeSixDof();
	void MainCallOf6DOF();
};

void InitializeSixDof();
void MainCallOf6DOF();

}

