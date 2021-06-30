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
//! @file      UnstructuredOversetConfig.cpp
//! @brief     Unstructured grid overset configuration.
//! @author    Chang Xinhua, He Kun.
#pragma once
#include "PHMpi.h"
using namespace std;

namespace PHSPACE
{
class OversetConfigFactoryMulti;
class OversetAuxilaryGrid;
class OversetInformationProxy;
class UnstructGrid;
class OversetConfig;
class InterCells;
class ZoneNodes;
class MRADTSearch;
class Pre_WalldistCompute;
class TimeSpan;

void RunOversetGridConfig();

void RunUnstructuredOversetGridConfig();    //! Perform an oversetConfig individually.

void InitializeOversetGridConfig();

void ReConfigUnstructOversetGrid();

void OutUnstructOversetGridConfig();

void FreeUnstructOversetGridConfig();

class OversetConfigFactoryMulti    //! Call OversetConfig as a whole, and complete the information communication between each zone.
{
public:
    OversetConfigFactoryMulti();
    ~OversetConfigFactoryMulti();

private:
    int numberOfBlocks;    //! Number of overlapping grids.
    int numberOfZones;     //! Total number of partitions in the grid.
    OversetConfig **oversetConfigMultiZone;    //! Corresponding to each grid partition.

private:
    ZoneNodes  **globalZoneNodes;     //! Set of all calculate points.
    InterCells **globalInterCells;    //! Set of all interpolation cells.
    TimeSpan *testTime;

private:
    vector< OversetAuxilaryGrid * > outerAuxilaryGrid;
    vector< UnstructGrid * >  innerGrid;

private:    //! This part of data is used to complete the attribute communication of the points on the interface.
    vector< int >           numberOfGlobalNodes;    //! Number of calculate points in the whole grid.
    vector< vector< int > > keyActiveOfGlobalInterNodes;    //! Property of points on the interface.

private:    //! Smallest box of each Block is used to improve the query efficiency.
    vector< vector< RDouble > > pMinOfBlocks;
    vector< vector< RDouble > > pMaxOfBlocks;

private:
    Pre_WalldistCompute *wallDistCompute;
    bool AdtMethodOrNod;

public:
    void Initialize();
    void OversetConfigMulti();

    void CalMinDistByGlobalNodes();
    void CalMinDistByAuxiliaryOuterGrid();

    void NodesOperator();
    void CellsOperator();
    void OutPutOversetFiles();
    void OutPutOversetInformation();

private:    //! Initialize.
    void InitializeGlobalValue();
    void SearchPhysicalAndInterfaceBoundary();
    void ReadInSklFiles();
    void ReadInInnerGrid();

    bool ReadInInnerGridH5(const string &fileName, UnstructGrid *gridIn);

    void ReadInOuterGrid();
    void GetGlobalInterfaceNodes();
    void GetGlobalInterfaceNodes(int iBlock);
    void ComputeNodeWalldist();

    void ReComputeWallDistance();
    void ComputeNodeAndCellTopology();

private:    //! Core process: points.
    void SetInitialValue();
    void BuildADTStructureAndMinBoxForSearch();
    void BuildADTStructure();
    void CalMinBoxForSearch();
    void SetBufferRegionOnOversetBoundary();
    void SetSearchRegionGlobal();
    void CalMinMaxBoxOfBlocks();
    void CreatGlobalNodesByFile();
    void WriteGlobalNodes();
    void ReadGlobalNodes();
    void CreatGlobalNodes();
    void CreatGlobalNodes_Old();
    void SpecifyTheDonorCellAndMinDistanceOfGlobalNodes();
    void SpecifyInitialInterNodesByDistance();
    void SpecifyNodesCharacterByDistanceAndInnerNodesCheck();
    void SpecifySpecicalCellsInteractWithWallBoundary();
    void SpecifyCellsInBodyByAuxiliaryInnerGrid();
    void FillTheActiveRegion();
    void SetNodeCharacterOnPhysicalBoundary();
    void SetNegativeForBufferRegion();
    void InterFaceBoundaryOptimize();
    void SetNegativeForInnerNodes();
    void SetNegativeForInnerAndInitialInterNodes();
private:    //! Core process: cells.
    void FindActiveCells();
    void FindInterCells();
    void CreatGlobalInterCells();
    void CreatGlobalInterCellsServerCollectAndBcast();
    void CreatGlobalInterCellsServerCollectAndBcastNonBlocking();
    void CreatGlobalInterCellsServerCollect();
    void InitializeInterCellByDataContainer(DataContainer *dataContainerRecv);
    void CreatGlobalInterCellsServerBcast();
    void CreatGlobalInterCells_Old();
    void FindDonorZoneAndDonorCell();
    void SelectDonorZoneAndDonorCell();

private:    //! Check.
    bool CheckDonorZoneAndDonorCell();
    bool CheckDonorZoneAndDonorCellServer();
    void PrintCellInformation();

private:    //! Communication of interface.
    void CharacterNodeToCell();
    void CharacterCellToNode();
    int  CharacterNodeToValidNodeByCell();

    void CommunicateNodeCharacterOnInterfaceBySKL(int key);
    void CommunicateNodeCharacterOnInterface();

private:    //! Reconstruction of interpolation information.
    void ResetOversetInformation();
    void GetInterAndDonorCellManagers();
    void GetInterAndDonorCellManagersNew();
    void GetNeighborZoneIndexContainerForTargetZone(int iZone, PHVectorInt1D &neighborZoneList);

private:    //! Output file.
    void OutKeyActive();
    void OutPutTecplot();
    void OutPutBoundaryCells();
    void OutInterCellInformation();
    void OutGlobalInterCells();
    void OutGlobalNodes();

public:
    int GetNumberOfBlocks(){ return numberOfBlocks; }

    vector< OversetAuxilaryGrid * > & GetOuterAuxilaryGrid() { return outerAuxilaryGrid; }

    vector< UnstructGrid * > & GetInnerGrid() { return innerGrid; }

    Pre_WalldistCompute * GetWallDistCompute() { return wallDistCompute; }
};

class OversetAuxilaryGrid
{
public:
    OversetAuxilaryGrid();
    ~OversetAuxilaryGrid();

public:
    int blockIndex;

    UnstructGrid *grid;
    vector < RDouble > nodeWallDistance;
    vector < RDouble > cellWallDistance;
    vector < RDouble > wallDistanceGradientX;
    vector < RDouble > wallDistanceGradientY;
    vector < RDouble > wallDistanceGradientZ;
    MRADTSearch  *adtSearch;
    Pre_WalldistCompute *walldistCompute;

public:
    void Initialize(UnstructGrid *gridIn, int iBlock, Pre_WalldistCompute *walldistComputeIn);
    void CalNodeWallDistance();
    void CalCellWallDistance();
    RDouble SearchAndCalMinDistByIDW(RDouble xc, RDouble yc, RDouble zc);
    RDouble SearchAndCalMinDistByVolumeCoordinate(RDouble xc, RDouble yc, RDouble zc);
    RDouble SearchAndCalMinDistByGradient(RDouble xc, RDouble yc, RDouble zc);
    void BuildADT();
    void CalMetrics();
    void CalGradientOfWallDistanceGGCell();
    void CalGradientOfWallDistanceGGNode();
    void CalGradientOfWallDistanceLSQ();

public:
    UnstructGrid * GetGrid() { return this->grid; }
    vector< RDouble > & GetNodeDistance() { return nodeWallDistance; }
    MRADTSearch * GetAdtSearch() { return adtSearch; }
};

class ZoneNodes 
{
public:
    ZoneNodes ();
    ~ZoneNodes ();

private:
    int numberOfNodes;
    int zoneIndex;
    int blockIndex;

    vector< RDouble > x;
    vector< RDouble > y;
    vector< RDouble > z;

    vector< RDouble > minDist;

    vector< RDouble > pMin;
    vector< RDouble > pMax;

public:
    void Initialize(int iZone, int iBlock, int numberOfNodes);

    void SetNumberOfNodes(int numberOfNodesIn) { this->numberOfNodes = numberOfNodesIn; }
    void SetZoneIndex    (int zoneIndex      ) { this->zoneIndex     = zoneIndex;       }
    void SetBlockIndex   (int blockIndex     ) { this->blockIndex    = blockIndex;      }
    
    int GetNumberOfNodes() { return numberOfNodes; }
    int GetZoneIndex    () { return zoneIndex;     }
    int GetBlockIndex   () { return blockIndex;    }

    vector< RDouble > & GetX() { return x; }
    vector< RDouble > & GetY() { return y; }
    vector< RDouble > & GetZ() { return z; }
    vector< RDouble > & GetMinDist() { return minDist; }

    vector< RDouble > & GetPMin() { return pMin; }
    vector< RDouble > & GetPMax() { return pMax; }
};

class InterCells
{
public:
    InterCells();
    ~InterCells();

private:
    int numberOfInterCells;
    int zoneIndex;
    int blockIndex;

    vector< int > donorZoneIndex;
    vector< int > donorCellIndex;
    vector< int > interCellIndex;

    vector< int > key;

    vector< RDouble > minDist;

    vector< RDouble > x;
    vector< RDouble > y;
    vector< RDouble > z;

public:
    void Initialize(int iZone, int iBlock, int numberOfInterCells);

    int GetZoneIndex()  { return zoneIndex;  }
    int GetBlockIndex() { return blockIndex; }

    int GetNumberOfInterCells() { return numberOfInterCells; }

    vector< int > & GetDonorZoneIndex() { return donorZoneIndex; }
    vector< int > & GetDonorCellIndex() { return donorCellIndex; }
    vector< int > & GetInterCellIndex() { return interCellIndex; }

    vector< RDouble > & GetX(){ return x; }
    vector< RDouble > & GetY(){ return y; }
    vector< RDouble > & GetZ(){ return z; }
    vector< RDouble > & GetMinDist() { return minDist; }
    vector< int > & GetKey() { return key; }

    void SetInterCellIndex(vector< int > interCellIndexIn) { this->interCellIndex = interCellIndexIn; }
    void SetCellCenterX   (vector< RDouble > xIn)          { this->x = xIn; }
    void SetCellCenterY   (vector< RDouble > yIn)          { this->y = yIn; }
    void SetCellCenterZ   (vector< RDouble > zIn)          { this->z = zIn; }
};

class OversetConfig
{
public:
    OversetConfig();
    ~OversetConfig();

public:    //! Necessary global information.
    int numberOfZones;
    int blockIndex;
    int zoneIndex;

    UnstructGrid *grid;
    ZoneNodes   **globalZoneNodes;
    InterCells  **globalInterCells;
    OversetInformationProxy   *oversetInformationProxy;    //! Information of interpolation calculated.
    OversetConfigFactoryMulti *oversetConfigFactoryMulti;

public:    //! Improve efficiency of query.
    MRADTSearch *aDTSearch;
    vector< RDouble > pMin;
    vector< RDouble > pMax;

    vector< RDouble > minDist;
public:
    int numberOfInterCells;

    vector< int > keyActiveOfNodes;
    vector< int > interCellIndex;

    vector< int > keyBufferOfNodes;
    vector< int > keyBufferOfCells;

    vector< int > keySearchOfNodes;    //! keySearch = 1£ºcells and points used for query.
    vector< int > keySearchOfCells;    //! keySearch = 0£ºpoints in the non overlapping areas (For example, point in the far field).

    bool *initialInterNodeOrNot;
    bool *innerCellsOrNot;
    bool *innerNodesOrNot;

public:
    vector< RDouble > nodeWalldistDistance;

public:    //! Boundary condition.
    vector< int > interfaceBoundaryNodesList;
    vector< int > physicalBoundaryNodesList;
    vector< int > innerBoundaryNodesList;

    vector< int > localToGlobalInterNodesIndex;    //! From local interNode to global interNode.
    vector< int > localToGlobalNodesIndex;         //! From local interNode to global Node.

public:    //! Initialize.
    void Initialize(OversetConfigFactoryMulti *oversetConfigFactoryMulti, UnstructGrid *gridIn, int nZone, int iZone, int iBlock, ZoneNodes  **globalZoneNodesIn, InterCells **globalInterCellsIn);
    void Allocate();
    void SetInitialValue();
    void BuildADTStructure();
    void CalculateMinBox();
    void SearchPhysicalAndInterfaceBoundary();
    void SearchBoundary();
    //void SearchInterfaceBoundary();
    void SetInitialValueForBufferRigion();
    void SetKeyBufferByKeyActive();

public:    //! Core operations: points.
    void ComputeNodeAndCellTopology();
    void CopyMinDist();
    void CalMinDistanceByAuxilaryGrid(vector< OversetAuxilaryGrid * > &auxilaryGrid);

    void SearchDonorCellAndComputeDistanceForGlobalNodes(int &numberOfSearchNodes, RDouble &time);
    bool TheTargetZoneOverLapWithMe(ZoneNodes *zoneNodes);
    RDouble SearchDonorCellAndComputeDistance(RDouble x, RDouble y, RDouble z);
    int  SearchDonorCell(RDouble xc, RDouble yc, RDouble zc);
    void SearchDonorCell(RDouble xc, RDouble yc, RDouble zc, int &donorCell);
    RDouble ComputeDistanceInDonorZone(int donorCell, RDouble xc, RDouble yc, RDouble RDouble);
    void SpecifyInitialInterNodesByDistance();
    void SpecifyNodesCharacterByDistanceAndInnerNodesCheck();
    void SpecifyCellsInBody(vector< MRADTSearch * > &adtSearchArray);
    void SpecifySpecicalCellsInteractWithWallBoundary();
    void SetNegativeForBufferRegion();
    void SetNegativeForInnerNodes();
    void SetNegativeForInnerAndInitialInterNodes();
    void SetNodeCharacterOnPhysicalBoundary();
    void FillTheRemainingNodes();
    
public:    //! Core operations: globalNodes.
    void CreatGlobalNodes(DataContainer *dataContainer);

public:    //! Core operations: points.
    void SetSearchRegionGlobal(vector< vector< RDouble > > &pMinOfBlocks, vector< vector< RDouble > > &pMaxOfBlocks);
    void SetKeySearchByKeyActive();

public:    //! Core operations: cells.
    void FindActiveCells();
    int  FindInterCells();
    void CalInterCellInformation(vector< int > &index,  vector< RDouble > &x,  vector< RDouble > &y,  vector< RDouble > &z);
    void FindInterCellsForTowOrderScheme();
    void FindDonorZoneAndDonorCell();
    int  FindDonorZoneAndDonorCellOfInterCell(int &donorCell, RDouble xc, RDouble yc, RDouble zc);

public:    //! Reconstruction of interpolation information.
    void InitializeOversetInformation();
    void PostProcessOversetInformation();

public:    //! Communication.
    void MarkGlobalInterfaceNodes(int iBlock, vector< int > &keyNodesGlobal);
    void CalLocalToGlobalInterfaceNodes(int iBlock, vector< int > &keyNodesGlobal);

    void CharacterNodeToCell();
    void CharacterCellToNode();
    void CharacterInterfaceToNode();
    int  CharacterNodeToValidNodeByCell();

    void UpdateInterfaceNodeCharacter(vector< vector < int > > &keyActiveOfGlobalInterNodes, int key);
    void DownInterfaceNodeCharacter(vector< vector< int > > &keyActiveOfGlobalInterNodes, int key);

public:    //! Check&Output.
    void ShowGrid();
    void ShowBoundaryCells();
    void PrintGridByOverSetConditionPoint(fstream &file, int key);
    void PrintGridByOverSetConditionCell(fstream &file, int key);
    void OutKeyBuffer();
    void OutKeyActive();
    void OutputLocalNodeWalldistance();

public:
    vector< int > & GetKeyActiveOfNodes() { return keyActiveOfNodes; }
    int           * GetKeyActiveOfCells();
    vector< int > & GetKeyBufferOfCells() { return keyBufferOfCells; }
    vector< int > & GetKeyBufferOfNodes() { return keyBufferOfNodes; }
    vector< int > & GetKeySearchOfCells() { return keySearchOfCells; }
    vector< int > & GetKeySearchOfNodes() { return keySearchOfNodes; }

    vector< int > & GetInterNodesList()               { return interfaceBoundaryNodesList; }
    vector< int > & GetLocalToGlobalNodesIndex()      { return localToGlobalNodesIndex; }
    vector< int > & GetLocalToGlobalInterNodesIndex() { return localToGlobalInterNodesIndex; }

    vector< RDouble > & GetPMin()                 { return pMin; }
    vector< RDouble > & GetPMax()                 { return pMax; }
    vector< RDouble > & GetNodeWalldistDistance() { return nodeWalldistDistance; }

    int GetBlockIndex()         { return blockIndex; }
    int GetNumberOfInterCells() { return numberOfInterCells; }
    vector< int > & GetInterCellIndex() { return interCellIndex; }

    OversetInformationProxy * GetOversetInformationProxy() { return oversetInformationProxy; }
    UnstructGrid * GetGrid() { return grid; }
};

RDouble GetToleranceForOversetSearch();
RDouble GetToleranceForOversetBox();
int GetTwoOrderInterpolationOrNot();
int GetKeyEnlargeOfNodes();
int GetReadInAuxiliaryInnerGrid();
int GetReadInAuxiliaryOuterGrid();
int GetReadInSklFileOrNot();
int GetSymetryOrNot();

template< typename T >
void PH_AllreduceSepModeForVector(vector< T > &sendingBuffer, vector< T > &receivingBuffer, int numberOfElements, MPI_Datatype mpiDataType, PH_Op op);

}