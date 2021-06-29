#pragma once

namespace PHSPACE
{

class MeshNodesSpring
{
public:
	MeshNodesSpring(){};
	~MeshNodesSpring(){};

public:
	int numberOfNodes;

public:
	void Initialize(){};
	void Run(){};
};

class MeshNodesMapping
{
public:
	MeshNodesMapping (){};
	~MeshNodesMapping(){};

public:
	int numberOfNodes;

public:
	void Initialize(){};
	void Run(){};
};

class MeshNodesRBF
{
public:
	MeshNodesRBF(){};
	~MeshNodesRBF(){};

public:
	int numberOfNodes;

public:
	void Initialize(){};
	void Run(){};
};

class SixDofParameter;
class DeformingParameter;
class Grid;

//最简单的运动网格类，实现的功能：变形+刚性运动。

class SimpleMovingMesh
{
public:
	SimpleMovingMesh();
	~SimpleMovingMesh();

public:
	int meshDoformingMethod;      //网格变形方法. 0-不变形；1-rbf；2-弹簧；3-mapping……

public:
    Grid* grid;

	int numberOfNodes;
	vector< RDouble > xInit;      //体坐标系内的初始坐标
	vector< RDouble > yInit;      //体坐标系内的初始坐标
	vector< RDouble > zInit;      //体坐标系内的初始坐标

	vector< RDouble > xPresent;   //当前在体坐标系内的坐标
	vector< RDouble > yPresent;   //当前在体坐标系内的坐标
	vector< RDouble > zPresent;   //当前在体坐标系内的坐标

public:
	void SetInitCoordinate( Grid* grid      );
	void CalInitCoordinate( SixDofParameter    * sixDof    );

	void Deforming        ( DeformingParameter * deforming );
	void RigidMoving      ( SixDofParameter    * sixDof    );
	void UpdateGrid       ();
};

//几个常见的“复杂”一点的运动网格
//从体坐标系，通过六自由度信息转换到计算坐标
//对于舵面控制问题：
//1-先在体轴系内完成旋转（舵偏）（将舵偏近似看成一种“变形”，）
//2-通过在弹体上的装配信息（装配六自由度）完成坐标转换
//3-通过弹体的六自由度信息转换到计算坐标

class OriginalMovingMesh
{
public:
	OriginalMovingMesh();
	~OriginalMovingMesh();

public:
	int bodyIndex;
	int fartherIndex;                        //“父”
	SixDofParameter    * location;         //“父”中的装配位置

	SixDofParameter    * presentSixdof;    //当前的刚体六自由的参数

	DeformingParameter * presentDeforming; //当前的变形参数

    Grid* grid;
	SimpleMovingMesh   * movingMesh;

public:
	void Initialize(Grid* grid, int bodyIndex );
	void Restart();
	void PreProcess();
	void PreProcessRigidMoving();
	void PreProcessDeforming();
	void Run();
};

class MovingMeshManager
{
public:
	MovingMeshManager();
	~MovingMeshManager();

public:
	int numberOfBodys;
	int numberOfZones;

	OriginalMovingMesh ** originalMovingMeshes;

	vector < OriginalMovingMesh * > innerMovingMesh;
	vector < OriginalMovingMesh * > outerMovingMesh;

public:
	void Initialize();
	void Restart();
	void Run();
	void Post();
};

}