//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//                                                                             +
//       Hybrid Platform for Engineering and Research of Flows (HyperFlow)     +
//                        面向流体工程与研究的混合CFD平台                      +
//                        （C）Zhang Laiping and He Xin                        +
//                    State Key Laboratory of Aerodynamics                     +
//                    Computational Aerodynamics Institute                     +
//            China Aerodynamics Research and Development Center               +
//                                 Since 2007                                  +
//                                                                             +
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#pragma once
#include "Precision.h"
#include "TypeDefine.h"
#include <fstream>
#include <vector>
using namespace std;

namespace PHSPACE
{
class UnstructGrid;
class ActionKey;
class DataContainer;
class FieldProxy;
class VirtualFile;

class GridManager
{
public:
    GridManager( UnstructGrid * unstructGrid );
    ~GridManager();
protected:
    UnstructGrid * unstructGrid;

    RDouble * cellCenterX;    // cell center data, including ghosts;
    RDouble * cellCenterY;    // cell center data, including ghosts;
    RDouble * cellCenterZ;    // cell center data, including ghosts;
    RDouble * faceNormalX;    // face unit normal;
    RDouble * faceNormalY;    // face unit normal;
    RDouble * faceNormalZ;    // face unit normal;
    RDouble * oldFaceNormalX; // n-1 grid face unit normal; 
    RDouble * oldFaceNormalY; // n-1 grid face unit normal; 
    RDouble * oldFaceNormalZ; // n-1 grid face unit normal; 
    RDouble * oldFaceArea;
    RDouble * cellVolume;     // face faceAreaContainer, cell cellVolume excluding ghosts;
    RDouble * faceCenterX;    // face center data;
    RDouble * faceCenterY;    // face center data;
    RDouble * faceCenterZ;    // face center data;
    RDouble * faceArea;       // face faceAreaContainer;

    RDouble * oldFaceCenterX; // old face center data;
    RDouble * oldFaceCenterY; // old face center data;
    RDouble * oldFaceCenterZ; // old face center data;

    RDouble * cellVolumeN1;   // cell cellVolume
    RDouble * cellVolumeN2;   // cell cellVolume
protected:
    RDouble * faceNormalVelocity;
    RDouble * faceVelocityX;
    RDouble * faceVelocityY;
    RDouble * faceVelocityZ;

    RDouble * oldFaceNormalVelocity;
    RDouble * oldFaceVelocityX;
    RDouble * oldFaceVelocityY;
    RDouble * oldFaceVelocityZ;

    RDouble * newFaceNormalVelocity;
    RDouble * newFaceVelocityX;
    RDouble * newFaceVelocityY;
    RDouble * newFaceVelocityZ;
protected:
    int * blankIndex;
protected:
    int * localToGlobalNodeIndexMapping;
    RDouble * globalCoordinateNew;
public:
    int * GetLocalToGlobalNodeIndexMapping() { return localToGlobalNodeIndexMapping; };
    RDouble * GetGlobalCoordinateNew() { return globalCoordinateNew; }

    void UpdateGridCoordinateByReadingFile                      ( int targetZone );
    void ComputeLocalToGlobalNodeIndexMappingDirectly           ( int targetZone );
    void ComputeLocalToGlobalNodeIndexMappingByReadingLayoutFile( int targetZone );
public:
    int *GetBlankIndex() { return blankIndex; };

    RDouble *GetCellCenterX() { return cellCenterX; };
    RDouble *GetCellCenterY() { return cellCenterY; };
    RDouble *GetCellCenterZ() { return cellCenterZ; };

    RDouble *GetFaceNormalX() { return faceNormalX;    };
    RDouble *GetFaceNormalY() { return faceNormalY;    };
    RDouble *GetFaceNormalZ() { return faceNormalZ;    };

    RDouble *GetOldFaceNormalX() { return oldFaceNormalX; };
    RDouble *GetOldFaceNormalY() { return oldFaceNormalY; };
    RDouble *GetOldFaceNormalZ() { return oldFaceNormalZ; };

    RDouble *GetFaceCenterX() { return faceCenterX;    };
    RDouble *GetFaceCenterY() { return faceCenterY;    };
    RDouble *GetFaceCenterZ() { return faceCenterZ;    };

    RDouble *GetOldFaceCenterX() { return oldFaceCenterX; };
    RDouble *GetOldFaceCenterY() { return oldFaceCenterY; };
    RDouble *GetOldFaceCenterZ() { return oldFaceCenterZ; };

    RDouble *GetFaceArea() { return faceArea;       };
    RDouble *GetOldFaceArea() { return oldFaceArea;    };

    RDouble *GetOldCellVolume() { return cellVolumeN1;   };
    RDouble *GetCellVolume( int iStep = 0 );
public:
    RDouble *GetFaceVelocityX() { return faceVelocityX;          };
    RDouble *GetFaceVelocityY() { return faceVelocityY;          };
    RDouble *GetFaceVelocityZ() { return faceVelocityZ;          };
    RDouble *GetFaceNormalVelocity() { return faceNormalVelocity;     };

    RDouble *GetOldFaceVelocityX() { return oldFaceVelocityX;       };
    RDouble *GetOldFaceVelocityY() { return oldFaceVelocityY;       };
    RDouble *GetOldFaceVelocityZ() { return oldFaceVelocityZ;       };
    RDouble *GetOldFaceNormalVelocity() { return oldFaceNormalVelocity;  };

    RDouble *GetNewFaceVelocityX() { return newFaceVelocityX;       };
    RDouble *GetNewFaceVelocityY() { return newFaceVelocityY;       };
    RDouble *GetNewFaceVelocityZ() { return newFaceVelocityZ;       };
    RDouble *GetNewFaceNormalVelocity() { return newFaceNormalVelocity;  };
public:
    void ComputeMassCharacterOfObject();
    void AllocateMetrics   ();
    void AllocateALEMetrics();
    void BackUpOldGrid();
    void ComputeMetrics(ActionKey* actkey);
    void ComputeGridFaceVelocity();
    void ComputeGridFaceVelocityDirectly();
    void ComputeGridFaceVelocitySatisfyingGCL();
    void ComputeGridFaceVelocity2DSatisfyingGCL();
    void ComputeGridFaceVelocity3DSatisfyingGCL();
    void DumpCellCenterFile();

protected:
    void ClosureCheck(RDouble *faceNormalX, RDouble *faceArea);
    void ComputeFaceNormal2D();
    void ComputeFaceCenter2D();
    void ComputeCellCenterVol2D(ActionKey* actkey);

    void ComputeFaceNormal3D();
    void ComputeFaceCenter3D();
    void ComputeCellCenterVol3D(ActionKey* actkey);

    void ComputeFaceNormal3DNew();
    void ComputeFaceCenter3DNew();
    void ComputeCellCenterVol3DNew();

public:
    void SkipReadLayoutFile( fstream & file, int targetZone );
    void ReadNewCoordinate();
    void UpdateGridNode();
public:
    void InitializeMovingGrids          ();
    void InitializeOldFaceNormalVelocity();
    void InitializeOldCellVolume        ();
public:
    void ReadMovingGrid();

    void UpdateOldCellVolume();
    void ConvergenceNorm( FieldProxy * q1Proxy, FieldProxy * q0Proxy, int numberOfEquations, RFloat & norm );
};

RDouble ComputeVolume4P( RDouble x1, RDouble y1, RDouble z1, RDouble x2, RDouble y2, RDouble z2, RDouble x3, RDouble y3, RDouble z3, RDouble x4, RDouble y4, RDouble z4 );

bool IsNotAllocated(void* pointer);

}
