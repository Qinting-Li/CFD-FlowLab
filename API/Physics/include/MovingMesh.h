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

// the father class of the moving mesh. the features includes the deforming and rigid motion.
class SimpleMovingMesh
{
public:
    SimpleMovingMesh();
    ~SimpleMovingMesh();

public:
    int meshDoformingMethod;      //! 0 - not doformed; 1 - rbf; 2 - spring; 3 - mapping¡­¡­

public:
    Grid *grid;

    int numberOfNodes;
    // the init coordinate of the body axis.
    vector< RDouble > xInit;
    vector< RDouble > yInit;
    vector< RDouble > zInit;

    // the current coordinate of the body axis.
    vector< RDouble > xPresent;
    vector< RDouble > yPresent;
    vector< RDouble > zPresent;

public:
    void SetInitCoordinate(Grid* grid     );
    void CalInitCoordinate(SixDofParameter    * sixDof   );

    void Deforming        (DeformingParameter * deforming);
    void RigidMoving      (SixDofParameter    * sixDof   );
    void UpdateGrid       ();
};


//! several common complex motion meshes.
//! transform to the calculating axis system from the body axis by the six DOF information.
//! for the rudder surface control.
//1 - first rotate in the body axis system.(rudder deviation)
//2 - implement the coordinate transform by the assembly information of the projectile.
//3 - transform to the calculating axis system by the six DOF information of the projectile.

class OriginalMovingMesh
{
public:
    OriginalMovingMesh();
    ~OriginalMovingMesh();

public:
    int bodyIndex;
    int fartherIndex;
    SixDofParameter    *location;         //! the assmbly positin in farther index.

    SixDofParameter    *presentSixdof;    //! the current six DOF parameters of the rigid body.

    DeformingParameter *presentDeforming; //! the current Deforming parameters.

    Grid* grid;
    SimpleMovingMesh   *movingMesh;

public:
    void Initialize(Grid* grid, int bodyIndex);
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