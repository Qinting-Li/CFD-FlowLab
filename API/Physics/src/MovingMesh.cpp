#include "PHMpi.h"
#include "AleManager.h"
#include "Force.h"
#include "SixDofManager.h"
#include "DeformingSolverManager.h"
#include "MovingMesh.h"
#include "Geo_Grid.h"
#include "Geo_UnstructGrid.h"
//#include "OversetConfigMulti.h"
#include "Region.h"
#include "IO_FileName.h"
#include "IO_FileReader.h"

namespace PHSPACE
{

SimpleMovingMesh::SimpleMovingMesh()
{
    meshDoformingMethod = 0;
    grid = NULL;
    numberOfNodes = 0;
}

SimpleMovingMesh::~SimpleMovingMesh()
{
    //if (grid != NULL) delete grid;
}

void SimpleMovingMesh::SetInitCoordinate(Grid* gridIn )
{
    this->grid = gridIn;

    RDouble *x = grid->GetX();
    RDouble *y = grid->GetY();
    RDouble *z = grid->GetZ();

    numberOfNodes = grid->GetNTotalNode();

    xInit.resize( numberOfNodes );
    yInit.resize( numberOfNodes );
    zInit.resize( numberOfNodes );

    xPresent.resize( numberOfNodes );
    yPresent.resize( numberOfNodes );
    zPresent.resize( numberOfNodes );

    for( int iNode = 0; iNode < numberOfNodes; ++ iNode )
    {
        xPresent[ iNode ] = xInit[ iNode ] = x[ iNode ];
        yPresent[ iNode ] = yInit[ iNode ] = y[ iNode ];
        zPresent[ iNode ] = zInit[ iNode ] = z[ iNode ];
    }
}

//根据sixdof，从inertialCoordinate转换到bodyCoordinate
void SimpleMovingMesh::CalInitCoordinate( SixDofParameter * sixDof )
{
    vector< RDouble > coordinateBody;
    vector< RDouble > coordinateInertial;
    vector< RDouble > angleVector;
    vector< RDouble > massCenter;

    coordinateBody.resize( 3 );
    coordinateInertial.resize( 3 );
    angleVector.resize( 3 );
    massCenter.resize( 3 );

    angleVector[ 0 ] = sixDof->GetEulerAngleGama();
    angleVector[ 1 ] = sixDof->GetEulerAngleBeta();
    angleVector[ 2 ] = sixDof->GetEulerAngleAlpha();
    massCenter [ 0 ] = sixDof->GetNondimensionalMassCenterX();
    massCenter [ 1 ] = sixDof->GetNondimensionalMassCenterY();
    massCenter [ 2 ] = sixDof->GetNondimensionalMassCenterZ();

    for( int iNode = 0; iNode < numberOfNodes; ++ iNode )
    {
        coordinateInertial[ 0 ] = xInit[ iNode ];
        coordinateInertial[ 1 ] = yInit[ iNode ];
        coordinateInertial[ 2 ] = zInit[ iNode ];

        ConvertVectorFromInertialCoordinateSystemToBodyCoordinateSystem( coordinateInertial, angleVector, massCenter, coordinateBody );

        xInit[ iNode ] = coordinateBody[ 0 ];
        yInit[ iNode ] = coordinateBody[ 1 ];
        zInit[ iNode ] = coordinateBody[ 2 ];
    }
}

void SimpleMovingMesh::Deforming( DeformingParameter * deforming )
{
    //暂时不考虑网格变形的问题
    for( int iNode = 0; iNode < numberOfNodes; ++ iNode )
    {
        xPresent[ iNode ] = xInit[ iNode ];
        yPresent[ iNode ] = yInit[ iNode ];
        zPresent[ iNode ] = zInit[ iNode ];
    }
}

//根据sixdof，从bodyCoordinate转换到presentCoordinate
void SimpleMovingMesh::RigidMoving( SixDofParameter * sixDof )
{
    vector< RDouble > coordinateBody;
    vector< RDouble > coordinateInertial;
    vector< RDouble > angleVector;
    vector< RDouble > massCenter;

    coordinateBody.resize( 3 );
    coordinateInertial.resize( 3 );
    angleVector.resize( 3 );
    massCenter.resize( 3 );

    angleVector[ 0 ] = sixDof->GetEulerAngleGama();
    angleVector[ 1 ] = sixDof->GetEulerAngleBeta();
    angleVector[ 2 ] = sixDof->GetEulerAngleAlpha();

    massCenter [ 0 ] = sixDof->GetNondimensionalMassCenterX();
    massCenter [ 1 ] = sixDof->GetNondimensionalMassCenterY();
    massCenter [ 2 ] = sixDof->GetNondimensionalMassCenterZ();

    for( int iNode = 0; iNode < numberOfNodes; ++ iNode )
    {
        coordinateBody[ 0 ] = xInit[ iNode ];
        coordinateBody[ 1 ] = yInit[ iNode ];
        coordinateBody[ 2 ] = zInit[ iNode ];

        ConvertVectorFromBodyCoordinateSystemToInertialCoordinateSystem( coordinateBody, angleVector, massCenter, coordinateInertial );

        xPresent[ iNode ] = coordinateInertial[ 0 ];
        yPresent[ iNode ] = coordinateInertial[ 1 ];
        zPresent[ iNode ] = coordinateInertial[ 2 ];
    }	
}

void SimpleMovingMesh::UpdateGrid()
{
    RDouble* x = grid->GetX();
    RDouble* y = grid->GetY();
    RDouble* z = grid->GetZ();

    int numberOfNodes = grid->GetNTotalNode();

    for( int iNode = 0; iNode < numberOfNodes; ++ iNode )
    {
        x[ iNode ] = xPresent[ iNode ];
        y[ iNode ] = yPresent[ iNode ];
        z[ iNode ] = zPresent[ iNode ];
    }
}

OriginalMovingMesh::OriginalMovingMesh()
{
    bodyIndex        = -1;
    fartherIndex     = -1;

    location         = NULL;
    movingMesh       = NULL;

    presentDeforming = NULL;
}

OriginalMovingMesh::~OriginalMovingMesh()
{
    if (location != NULL) { delete location; location = NULL; }
    if (movingMesh != NULL) { delete movingMesh; movingMesh = NULL; }
    //if (presentSixdof != NULL) { delete presentSixdof; presentSixdof = NULL; }
    //if (grid != NULL) { delete grid; grid = NULL; }
    if (presentDeforming != NULL) { delete presentDeforming; presentDeforming = NULL; }
}

void OriginalMovingMesh::Initialize( Grid * grid, int bodyIndex )
{
    this->grid      = grid;
    this->bodyIndex = bodyIndex;
    int numberOfMovingBodies = GlobalDataBase::GetIntParaFromDB("numberOfMovingBodies");
    if( bodyIndex >= numberOfMovingBodies) return;

    movingMesh = new SimpleMovingMesh;

    fartherIndex = GetIntegerParameterFromDataBase( bodyIndex, "fartherIndex" );
    if( fartherIndex != -1 )
    {
        vector< RDouble > configPamameter;
        GetRDoubleVectorFromDataBase( configPamameter, bodyIndex, "configPamameter", 6 );
        location  = new SixDofParameter;
        location->SetMassCenter( configPamameter[ 0 ], configPamameter[ 1 ], configPamameter[ 2 ] );
        location->SetEulerAngle( configPamameter[ 3 ], configPamameter[ 4 ], configPamameter[ 5 ] );
    }
}

void OriginalMovingMesh::PreProcessRigidMoving()
{
    presentSixdof = PHSPACE::GetSixDofParameter( bodyIndex );

    movingMesh->SetInitCoordinate( grid );

    //如果该网格块作为一个部件装配在一个物体上
    if( fartherIndex != -1 )
    {
        SixDofParameter * sixdofFarther = PHSPACE::GetSixDofParameter( fartherIndex );

        movingMesh->CalInitCoordinate( sixdofFarther );
        movingMesh->CalInitCoordinate( location );
    }

    movingMesh->CalInitCoordinate( presentSixdof );
}

void OriginalMovingMesh::PreProcessDeforming()
{

}

//动网格的预处理。得到体轴系内的坐标点；Mapping方法求面积坐标；rbf方法求矩阵
void OriginalMovingMesh::PreProcess()
{
    int numberOfMovingBodies = GlobalDataBase::GetIntParaFromDB("numberOfMovingBodies");
    if (bodyIndex >= numberOfMovingBodies) return;

    PreProcessRigidMoving();  //刚体运动，体坐标系内的坐标点

    PreProcessDeforming();
}	

void OriginalMovingMesh::Restart()
{
    int numberOfMovingBodies = GlobalDataBase::GetIntParaFromDB("numberOfMovingBodies");
    if (bodyIndex >= numberOfMovingBodies) return;

    presentDeforming = PHSPACE::GetDeformingParameter( bodyIndex );
    presentSixdof    = PHSPACE::GetSixDofParameter   ( bodyIndex );

    movingMesh->Deforming  ( presentDeforming );
    movingMesh->RigidMoving( presentSixdof    );

    if( fartherIndex != -1 )
    {
        SixDofParameter * sixdofFarther = PHSPACE::GetSixDofParameter( fartherIndex );

        movingMesh->RigidMoving( location      );
        movingMesh->RigidMoving( sixdofFarther );
    }

    movingMesh->UpdateGrid();
}

void OriginalMovingMesh::Run()
{
    int numberOfMovingBodies = GlobalDataBase::GetIntParaFromDB("numberOfMovingBodies");
    if (bodyIndex >= numberOfMovingBodies) return;

    presentDeforming = PHSPACE::GetDeformingParameter( bodyIndex );
    presentSixdof    = PHSPACE::GetSixDofParameter   ( bodyIndex );

    movingMesh->Deforming  ( presentDeforming );
    movingMesh->RigidMoving( presentSixdof    );

    if( fartherIndex != -1 )
    {
        SixDofParameter * sixdofFarther = PHSPACE::GetSixDofParameter( fartherIndex );

        movingMesh->RigidMoving( location      );
        movingMesh->RigidMoving( sixdofFarther );
    }

    movingMesh->UpdateGrid();
}

MovingMeshManager * movingMeshManager = NULL;

MovingMeshManager::MovingMeshManager()
{
    numberOfBodys        = 0;
    numberOfZones        = 0;

    originalMovingMeshes = NULL;
}

MovingMeshManager::~MovingMeshManager()
{
    for( int iZone = 0; iZone < numberOfZones; ++ iZone )
    {
        if (originalMovingMeshes[iZone])
        {
            delete originalMovingMeshes[iZone];
            originalMovingMeshes[iZone] = NULL;
        }
    }

    delete [] originalMovingMeshes;
    originalMovingMeshes = NULL;

    for( int num = 0; num < innerMovingMesh.size(); ++ num )
    {
        delete innerMovingMesh[ num ];
        innerMovingMesh[num] = NULL;
    }

    for( int num = 0; num < outerMovingMesh.size(); ++ num )
    {
        delete outerMovingMesh[ num ];
        outerMovingMesh[num] = NULL;
    }
}

void MovingMeshManager::Initialize()
{
    numberOfBodys = PHSPACE::GlobalDataBase::GetIntParaFromDB("numberOfMovingBodies");
    numberOfZones = PHMPI::GetNumberofGlobalZones();

    originalMovingMeshes = new OriginalMovingMesh * [ numberOfZones ];

    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    for( int iZone = 0; iZone < numberOfZones; ++ iZone )
    {
        int processorIndex = PHMPI::GetZoneProcessorID( iZone );
        if ( currentProcessorIndex == processorIndex )
        {
            originalMovingMeshes[ iZone ] = new OriginalMovingMesh;
        }
        else
        {
            originalMovingMeshes[ iZone ] = NULL;
        }
    }

    //区别哪一块网格属于哪一个物体，暂时根据读入的网格文件个数来判断。
    //进一步可以通过网格的体属性来判断，待完成。
    for( int iZone = 0; iZone < numberOfZones; ++ iZone )
    {
        int processorIndex = PHMPI::GetZoneProcessorID( iZone );
        if ( currentProcessorIndex == processorIndex )
        {
            Grid * grid =  GetGrid( iZone , 0);

            //就是这个地方，需要进一步完善交互策略！
            int bodyIndex = grid->GetIBlock();

            originalMovingMeshes[ iZone ]->Initialize( grid, bodyIndex );
        }
    }

    //如果是变形网格，就需要在这里进行一些整体上的处理，例如建立背景网格，建立物面点的RBF矩阵等
    //然后
    for( int iZone = 0; iZone < numberOfZones; ++ iZone )
    {
        int processorIndex = PHMPI::GetZoneProcessorID( iZone );
        if ( currentProcessorIndex == processorIndex )
        {
            originalMovingMeshes[ iZone ]->PreProcess();
        }
    }

    ////辅助网格的初始化: innerGrid
    //if( GetOversetConfigFactoryMulti() != 0 && GetReadInAuxiliaryInnerGrid() == 1 )
    //{	
    //	vector< UnstructGrid * > & innerGrids = GetOversetConfigFactoryMulti()->GetInnerGrid();

    //	for( int iBody = 0; iBody < numberOfBodys; ++ iBody )
    //	{
    //		if( innerGrids[ iBody ] != 0 )
    //		{
    //			UnstructGrid * innerGrid = innerGrids[ iBody ]; 

    //			OriginalMovingMesh * movingMesh = new OriginalMovingMesh();

    //			movingMesh->Initialize( innerGrid, iBody );
    //			movingMesh->PreProcess();

    //			innerMovingMesh.push_back( movingMesh );
    //		}
    //	}
    //}

    ////辅助网格的初始化: outerGrid
    //if( GetOversetConfigFactoryMulti() != 0 && GetReadInAuxiliaryOuterGrid() == 1 )
    //{
    //	vector< FYOversetAuxilaryGrid * > & outerAuxilaryGrids = GetOversetConfigFactoryMulti()->GetOuterAuxilaryGrid();

    //	for( int iBody = 0; iBody < numberOfBodys; ++ iBody )
    //	{
    //		UnstructGrid * outerGrid = outerAuxilaryGrids[ iBody ]->GetGrid();

    //		OriginalMovingMesh * movingMesh = new OriginalMovingMesh();

    //		movingMesh->Initialize( outerGrid, iBody );
    //		movingMesh->PreProcess();

    //		outerMovingMesh.push_back( movingMesh );
    //	}
    //}
}

void MovingMeshManager::Restart()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    for( int iZone = 0; iZone < numberOfZones; ++ iZone )
    {
        int processorIndex = PHMPI::GetZoneProcessorID( iZone );
        if ( currentProcessorIndex == processorIndex )
        {
            originalMovingMeshes[ iZone ]->Restart();
        }
    }

    //辅助网格的初始化
    for( int iInnerMesh = 0; iInnerMesh < innerMovingMesh.size(); ++ iInnerMesh )
    {
        innerMovingMesh[ iInnerMesh ]->Restart();
    }

    for( int iOuterMesh = 0; iOuterMesh < outerMovingMesh.size(); ++ iOuterMesh )
    {
        outerMovingMesh[ iOuterMesh ]->Restart();
    }
}

void MovingMeshManager::Run()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    for( int iZone = 0; iZone < numberOfZones; ++ iZone )
    {
        int processorIndex = PHMPI::GetZoneProcessorID( iZone );
        if ( currentProcessorIndex == processorIndex )
        {
            originalMovingMeshes[ iZone ]->Run();
        }
    }

    for( int iInnerMesh = 0; iInnerMesh < innerMovingMesh.size(); ++ iInnerMesh )
    {
        innerMovingMesh[ iInnerMesh ]->Run();
    }

    for( int iOuterMesh = 0; iOuterMesh < outerMovingMesh.size(); ++ iOuterMesh )
    {
        outerMovingMesh[ iOuterMesh ]->Run();
    }
}

void MovingMeshManager::Post()
{
    ;
}

}