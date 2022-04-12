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

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// 文件头部说明：                                                         +
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// FILE:       ParallelmeshGen.h                                          +
// AUTHOR(S):  Bell                                                       +
// PURPOSE:    Parallel unstrucutred mesh generation                      +
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#pragma once
#include "ZoneGridSolver.h"

namespace PHSPACE
{
class CGNSMultiZone;
class SuperElement;
class Zone;

typedef PHVector1D< Point3D > PointArray;
typedef vector< PointArray >  PointVector;

class ParallelMeshGen
{
public:
    ParallelMeshGen();
    ~ParallelMeshGen();
private:
    int nZones;
    SuperElement **super_elements;
    CGNSMultiZone *cgns_multizone;
    Region * region;
    ZoneGridSolver ** zoneGridSolverContainer;
    ZoneGridSolver * currentZoneGridSolver;
    Grid ** fantasyGrids;

    PHVector2D < int > ZERO_VECTOR2D;

    //ComputationalElementSolver * computationalElementSolver;

//  PHAdaptiveFlowFieldSolver * currentAdaptiveFlowFieldSolver;
//  vector< PHAdaptiveFlowFieldSolver * > adaptiveFlowFieldSolverVector;
//  PHAdaptiveFlowFieldSolver * currentAdaptiveFlowFieldSolver;
public:
    void ParallelGridGen();
    void ReadZoneGridFromFTS();
    void AllocateZoneGridSolverContainer();
    void ConstructZoneGridSolver();
    void VisulizationOriginalGrid();
    void RecoverZoneGridSolverFromFTS();
    void RecoverZoneGridSolverFromFTS(Zone * zone);
    void RecoverElements(Zone * zone);
    void RecoverPoints(Zone * zone);
    void RecoverFacess(Zone * zone);
    void RecoverMetrics(Zone * zone);
    void RecoverInterFaceInformation(Zone * zone);

    void BroadcastNeighborZonesToEachProcessor();
    void BroadcastNeighborZonesToEachProcessorNew();

    Zone * GetZone(int izone);
    bool ZoneBelongToCurrentProcessor(int zoneIndex);
    PHVector1D < int > GetNeighborZones(ZoneGridSolver * zoneGridSolver);
    void BroadcastToEachProcessor(const int zoneIndex, PHVector1D < int > & dataOfThisProcessor);
    void AddOneZoneNeighborToAllZonesNeighbor(PHVector1D < int > & dataOfThisZone, int zoneIndex, int targetZoneIndex);
    void InitialInterFaceInformation();
    void SwapNeighborsSendContent();
    PHVector1D < int > GetAllOtherZones(int zoneIndex);
    PHVector2D < int > & GetCurrentNeighborZones(int zoneIndex);
    ZoneGridSolver * GetZoneGridSolver(int iZone) { return zoneGridSolverContainer[ iZone ]; }

    void BuildElementAdaptProperty();
    void DoSomethingIfResetElementPropertyBySelfDefine();
    void ResetElementPropertyBySelfDefine();
    void ReComputeAdaptionProperty();

    void AdaptGrid();
    void AdaptGrid(ZoneGridSolver * zoneGridSolver);
    void GridAdaptation(ZoneGridSolver * zoneGridSolver);
    void SplitAllElements();
    void ConstructBasicData();
    void GenerateElementIndexMapping();

    void ShapePreserving();

    void ConstructGrid();
    void ReadGrid();
    void ResolveComputationGrid();
    void CreateGeneralMultiBlockGrid();
    int  GetGridType(int izone) const;
    Grid * CreateBasicGrid(int izone);

    int GetElementTypeOfFTS(const int & dimension, const int & nodeNumberOfCell);
    int GetFaceTypeOfFTS(const int & dimension, const int & nodeNumberOfCell);
    void SetCurrentZoneGridSolver(ZoneGridSolver * currentZoneGridSolver) { this->currentZoneGridSolver = currentZoneGridSolver; }

    void MeshQualityChecking();

    void GenerateAndDumpComputationalGrid();
    void AllocateComputationalGrid();
    void GenerateComputationalGrid();
    void DumpComputationalGrid();
    void DumpMultiZoneComputationalGrids(const string & fileName, Grid ** gridContainer, int numberOfZones);
    void DumpMultiZoneComputationalGridsCSMode(const string & fileName, Grid ** gridContainer, int numberOfZones);
    void DumpMultiZoneComputationalGridsP2PMode(const string & fileName, Grid ** gridContainer, int numberOfZones);
    void WriteFaceBC(const string &gridfileName, Grid ** gridContainer, int numberOfZones);
    void WriteCellToNode(const string & cellToNodeFileName);

private:
    void CreateAllSingleGrid(Grid ** grids);
    bool IfNeedProjection();
    void SetAdaptTropic();

    bool IsNeedInfectInterfaceIsotropicRefineType();
    bool IsNeedInfectInterfaceIsotropicRefineType(ZoneGridSolver * zoneGridSolver);

    void InitCellRefineType();
    void InitCellRefineType(ZoneGridSolver * zoneGridSolver);

    void FixCellRefineType();
    void FixCellRefineType(ZoneGridSolver * zoneGridSolver);

    void SwapInterfaceRefineType();
    void SwapInterfaceRefineType(int zoneIndex, int neighborZoneIndex, int sendProcessorIndex, int receiveProcessorIndex);

    void CompressVector(DataContainer *& dataContainer, int globalZoneIndex, int neighborZoneIndex);
    void DecompressVector(DataContainer *& dataContainer, int globalZoneIndex, int neighborZoneIndex);

    void InitGlobalBoundaryCondition();

    void PreparePointsCoordinates(PointVector *pointsCoord);
    void PrepareCheckIsSolidWallPoint(vector< vector<bool> > *isSolidWallPoint, vector< vector<bool> > *needDeform);
    void PrepareCheckIsSymmetryPoint(vector< vector<bool> > *isSolidWallPoint);
    void PrepareBoundaryDeformedPoint(vector< vector<bool> > *isSolidWallPoint, PointVector *zoneProjectVector,
        PointVector *boundaryDeformedPoints, PointVector *deformedVector);

    void UpdatePointsCoordinates(PointVector *pointsCoord, vector< vector<bool> > *needDeform);
};

void CompressFixCellType(ZoneGridSolver * zoneGridSolver, DataContainer *& dataContainer, int neighborZoneIndex);
void DecompressFixCellType(ZoneGridSolver * zoneGridSolver, DataContainer *& dataContainer, int neighborZoneIndex);

void CompressVectorIntoDataContainer(DataContainer *& dataContainer, PHVector1D < int > & data);
void AnasysVectorFromDataContainer(DataContainer *& dataContainer, PHVector1D < int > & data);

}
