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
//! @file      UnstructuredOversetConfig.h
//! @brief     Unstructured grid overset configuration.
//! @author    Chang Xinhua, He Kun, Zhang Laiping.
#pragma once
#include "PHMpi.h"
#include "OversetKDTree.h"
using namespace std;

namespace PHSPACE
{
class OversetConfigFactoryMulti;
class OversetInformationProxy;
class UnstructGrid;
class OversetConfig;
class InterCells;
class ZoneSearchNodes;
class Pre_WalldistCompute;
class TimeSpan;

void RunOversetGridConfig();

void RunUnstructuredOversetGridConfig();    //! Perform an oversetConfig individually.

void InitializeOversetGridConfig();

void ReConfigUnstructOversetGrid();

void OutUnstructOversetGridConfig();

void FreeUnstructOversetGridConfig();

OversetConfigFactoryMulti *GetOversetConfigFactoryMulti();

class OversetConfigFactoryMulti    //! Call OversetConfig as a whole, and complete the information communication between each zone.
{
public:
    OversetConfigFactoryMulti();
    ~OversetConfigFactoryMulti();

private:
    int numberOfBlocks;    //! Number of overlapping grids.
    int numberOfZones;     //! Total number of partitions in the grid.
    OversetConfig **oversetConfigMultiZone;    //! Corresponding to each grid partition.
    int geometricDimension;
private:
    ZoneSearchNodes **zonesSearchNodes;     //! Set of all calculate points.
    InterCells **zonesInterpolateCells;    //! Set of all interpolation cells.
    TimeSpan *testTime;

private:
    vector< UnstructGrid * > innerGrid;

private:    //! This part of data is used to complete the attribute communication of the points on the interface.
    vector< int > numberOfGlobalNodes;    //! Number of calculate points in the whole grid.
    vector< vector< int > > keyActiveOfGlobalInterNodes;    //! Property of points on the interface.

    int *numberOfInterfaceNodesInBlocks;
    int **keyOfInterfaceNodesInBlocks;
private:    //! Smallest box of each Block is used to improve the query efficiency.
    RDouble **pMinOfBlocks;
    RDouble **pMaxOfBlocks;

public:
    void Initialize();
    void OversetConfigMulti();
    void FreeKeyBasicActiveNodes();
    void BuildNodesMinDistanceInOtherBlocks();
    void BuildNodesMinDistanceInOtherBlocks2();
    void BuildLocalToGlobalInterfaceNodesMaps();

    void NodesOperator();
    void CellsOperator();
    void OutPutOversetFiles();
    void OutPutOversetInformation();
    void OutPutOversetVisualizations(const string &name);
public:
    RDouble **GetMinBoxOfBlocks() {
        return this->pMinOfBlocks;
    }
    RDouble **GetMaxBoxOfBlocks() {
        return this->pMaxOfBlocks;
    }
    InterCells **GetZonesInterpolateCells() {
        return this->zonesInterpolateCells;
    }
    ZoneSearchNodes **GetZonesSearchNodes() {
        return this->zonesSearchNodes;
    }

    int **GetKeyOfInterfaceNodesInBlocks() {
        return this->keyOfInterfaceNodesInBlocks;
    }

private:    //! Initialize.
    void BuildNodesMinDistanceInThisBlocks(vector< RDouble > *&wallNodeCoordinateXIn,
        vector< RDouble > *&wallNodeCoordinateYIn,
        vector< RDouble > *&wallNodeCoordinateZIn);
    void ReadInInnerGrid();
    void SetNodesMinDistanceInThisBlocksByWalldist();
    void SetNodesMinDistanceInOtherBlocks();

    bool ReadInInnerGridH5(const string &fileName, UnstructGrid *gridIn);
    void BuildLocalToGlobalInterfaceNodesMap(int iBlockIn, int *numberOfNodesInBlocksIn);

    void GetGlobalInterfaceNodes();
    void GetGlobalInterfaceNodes(int iBlock);
    void ComputeNodeWalldist();
    void ComputeNodeWalldistInOtherBlock();
    void ReComputeWallDistance();
    void ComputeNodeAndCellTopology();

private:    //! Core process: points.
    void SetInitialValue();
    void CalMinBoxForSearch();
    void SetBufferRegionOnOversetBoundary();
    void BuildKDTreesBySearchRegions();
    void BuildNodesMinDistanceInOtherBlocksByKDTrees();
    void SetActivePropertyForSearchRegions();
    void SetSearchRegionByActivePropertys();
    void BuildMinMaxBoxOfZoneBySearchRegions();
    void BuildKDTreeByMinMaxBoxOfZones();
    void BuildInBlockCellsAndNodesByInnerAuxiliaryGridKDTrees(vector< GridKDTree * > &innerAuxiliaryGridKDTreesIn);

    void BuildMinMaxBoxOfBlocks();

    void BuildZonesSearchNodes();
    void SearchZonesSearchNodesInKDTrees();
    void AllReduceNodesMinDistanceInOtherBlock();
    void SetNodesMinDistanceInOtherBlockByZonesNodes();

    void SpecifyInitialInterNodesByDistance();
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
    void AllReduceZonesInterpolateCells();
private:    //! Core process: cells.
    void BuildInterpolateRelationships();
    void FindInterpolateRelationships();

private:    //! Check.
    void CheckInterpolateRelationships();
    bool CheckDonorZoneAndDonorCellServer();
    void PrintCellInformation();

private:    //! Communication of interface.
    void CharacterNodeToCell();
    void CharacterCellToNode();
    int  CharacterNodeToValidNodeByCell();

    void CommunicateNodePropertyByKeyOfInterfaceNodesInBlocks();
    void CommunicateNodeCharacterOnInterface();
    void CommunicateCellIBlank();

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
    int GetNumberOfBlocks() {
        return numberOfBlocks;
    }

    vector< UnstructGrid * > &GetInnerGrid() {
        return innerGrid;
    }
};

class ZoneSearchNodes
{
public:
    ZoneSearchNodes(int numberOfSearchNodesIn, int blockIndexIn);
    ~ZoneSearchNodes();

private:
    int numberOfNodes;
    int blockIndex;

    RDouble *x;
    RDouble *y;
    RDouble *z;

    RDouble *nodesMinDistanceToOtherBlock;

    RDouble *pMin;
    RDouble *pMax;

public:

    void SetNumberOfNodes(int numberOfNodesIn) {
        this->numberOfNodes = numberOfNodesIn;
    }
    void SetBlockIndex(int blockIndex) {
        this->blockIndex = blockIndex;
    }

    int GetNumberOfNodes() {
        return numberOfNodes;
    }
    int GetBlockIndex() {
        return blockIndex;
    }

    RDouble *GetX() {
        return x;
    }
    RDouble *GetY() {
        return y;
    }
    RDouble *GetZ() {
        return z;
    }
    RDouble *GetNodesMinDistanceToOtherBlock() {
        return nodesMinDistanceToOtherBlock;
    }

    RDouble *GetPMin() {
        return pMin;
    }
    RDouble *GetPMax() {
        return pMax;
    }
};

class InterCells
{
public:
    InterCells(int numberOfInterpolateCellsIn, int blockIndexIn);
    ~InterCells();

private:
    int numberOfInterCells;
    int blockIndex;

    int *donorZoneIndex;
    int *donorCellIndex;

    int *donorLevel;
    RDouble *cellMinDistance;

    RDouble *x;
    RDouble *y;
    RDouble *z;

public:
    void FreeUselessMemery();

    int GetBlockIndex() {
        return blockIndex;
    }

    int GetNumberOfInterCells() {
        return numberOfInterCells;
    }

    int *GetDonorZoneIndex() {
        return donorZoneIndex;
    }
    int *GetDonorCellIndex() {
        return donorCellIndex;
    }

    RDouble *GetX() {
        return x;
    }
    RDouble *GetY() {
        return y;
    }
    RDouble *GetZ() {
        return z;
    }
    RDouble *GetMinDist() {
        return cellMinDistance;
    }
    int *GetKey() {
        return donorLevel;
    }

    void SetCellCenterX(RDouble *xIn) {
        this->x = xIn;
    }
    void SetCellCenterY(RDouble *yIn) {
        this->y = yIn;
    }
    void SetCellCenterZ(RDouble *zIn) {
        this->z = zIn;
    }
};

class OversetConfig
{
    typedef void ( OversetConfig:: *BuildInBlockCellsAndNodesMethod )( vector< GridKDTree * > & );
public:
    OversetConfig();
    ~OversetConfig();

public:    //! Necessary global information.
    int numberOfZones;
    int blockIndex;
    int zoneIndex;
    int geometricDimension;

    UnstructGrid *grid;

    OversetInformationProxy *oversetInformationProxy;    //! Information of interpolation calculated.
    BuildInBlockCellsAndNodesMethod buildInBlockCellsAndNodesMethod;

public:    //! Improve efficiency of query.
    GridKDTree *gridKDTree;
    RDouble *minBoxOfZone;
    RDouble *maxBoxOfZone;

public:
    int numberOfInterCells;

    int *localToGlobalInterfaceNodesMap;
    int *keyActiveOfNodes;

    vector< int > interpolateCellList;

    int *keyBufferOfNodes;

    int *keyBufferOfCells;

    int *keySearchOfNodes;    //! keySearch = 1£ºcells and points used for query.

    int *keySearchOfCells;    //! keySearch = 0£ºpoints in the non overlapping areas (For example, point in the far field).

    int *keyBasicActiveNodes;

    //!the mark of the node which is inside of the body of other block
    //!0 - outside
    //!1 - inside
    int *keyInBlockNodes;

    RDouble *nodesMinDistanceInThisBlock;

    RDouble *nodesMinDistanceInOtherBlock;

public:
    vector< RDouble > nodeWalldistDistance;

public:    //! Boundary condition.
    vector< int > interfaceBoundaryNodeList;
    vector< int > physicalBoundaryNodeList;
    vector< int > innerBoundaryNodeList;

    vector< int > localToGlobalInterNodesIndex;    //! From local interNode to global interNode.
    vector< int > localToGlobalNodesIndex;         //! From local interNode to global Node.
public:
    void BuildBoundaryNodeList(vector< RDouble > &wallNodeCoordinateXIn,
        vector< RDouble > &wallNodeCoordinateYIn,
        vector< RDouble > &wallNodeCoordinateZIn);
    void BuildNodesMinDistanceInThisBlock(vector< RDouble > *wallNodeCoordinateXInBlocksIn,
        vector< RDouble > *wallNodeCoordinateYInBlocksIn,
        vector< RDouble > *wallNodeCoordinateZInBlocksIn);

    void BuildBoundaryNodeList();
    void SetNodesMinDistanceInThisBlocksByWalldist();
    void SetNodesMinDistanceInOtherBlocks();

public:    //! Initialize.
    void Initialize(UnstructGrid *gridIn, int nZone, int iZone, int iBlock);
    void SetInitialValue();
    void FreeKeyBasicActiveNode();
    void FreeKeyInBlockNodes();
    void SetInitialValueForBufferRigion();
    void SetKeyBufferByKeyActive();

public:    //! Core operations: points.
    void ComputeNodeAndCellTopology();
    void CopyMinDist(RDouble *nodesMinDistanceInOtherBlockIn);
    void BuildKDTreeByMinMaxBoxOfZone();
    void SearchZonesSearchNodesInKDTree();

    RDouble ComputeDistanceInDonorZone(RDouble *nodeCoordinateIn, int donorCell);
    void SpecifyInitialInterNodesByDistance();
    void SpecifySpecicalCellsInteractWithWallBoundary();
    void SetNegativeForBufferRegion();
    void SetNegativeForInnerNodes();
    void SetNegativeForInnerAndInitialInterNodes();
    void SetNodeCharacterOnPhysicalBoundary();
    void FillTheRemainingNodes();

public:    //! Core operations: globalNodes.
    void BuildMinMaxBoxOfZone();

public:    //! Core operations: points.
    void SetActivePropertyForSearchRegion();
    void SetSearchRegionByActiveProperty();
    void BuildMinMaxBoxOfZoneBySearchRegion();

public:    //! Core operations: cells.
    void FindActiveCells();
    void FindInterCells();
    int  FindDonorCell(RDouble *nodeCoordinateIn);
    void FindInterpolateRelationship();
    void CheckInterpolateRelationship(int iBlockIn,
        InterCells *zoneInterpolateCellsIn,
        int &numberOfInterpolateCellsInThisBlockIn,
        int &numberOfNegativeCellsInThisBlock_1In,
        int &numberOfNegativeCellsInThisBlock_2In,
        int &numberOfNegativeCellsInThisBlock_3In);

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
    void UpdateInterfaceNodeCharacter();
    void DownInterfaceNodeCharacter();

public:
    void BuildInBlockCellsAndNodesByInnerAuxiliaryGridKDTree(vector< GridKDTree * > &innerAuxiliaryGridKDTreesIn);

public:    //! Check&Output.
    void OutPutOversetVisualization(const string &name);
    void OutPutOversetVisualizationHeader(fstream &file);
    void OutPutOversetVisualizationCell(fstream &file);
    void OutPutOversetVisualizationPoint(fstream &file);
    void OutPutOversetVisualizationTail(fstream &file);
    void ShowGrid();
    void ShowBoundaryCells();
    void PrintGridByOverSetConditionPoint(fstream &file, int key);
    void PrintGridByOverSetConditionCell(fstream &file, int key);
    void OutKeyBuffer();
    void OutKeyActive();
    void OutputLocalNodeWalldistance();

    inline bool IfZoneSearchNodesIntersectionWithThisZone(RDouble *minBoxOfZoneSearchNodesIn, RDouble *maxBoxOfZoneSearchNodesIn)
    {
        for (int iDimension = 0; iDimension < geometricDimension; ++ iDimension)
        {
            if (maxBoxOfZoneSearchNodesIn[iDimension] < minBoxOfZone[iDimension]) return false;
            if (minBoxOfZoneSearchNodesIn[iDimension] > maxBoxOfZone[iDimension]) return false;
        }
        return true;
    }

public:
    int *GetKeyActiveOfNodes()
    {
        return keyActiveOfNodes;
    }
    int *GetKeyActiveOfCells();
    int *GetKeyBufferOfCells()
    {
        return keyBufferOfCells;
    }
    int *GetKeyBufferOfNodes()
    {
        return keyBufferOfNodes;
    }
    int *GetKeySearchOfCells()
    {
        return keySearchOfCells;
    }
    int *GetKeySearchOfNodes()
    {
        return keySearchOfNodes;
    }

    vector< int > &GetInterfaceBoundaryNodeList()
    {
        return interfaceBoundaryNodeList;
    }

    int *GetLocalToGlobalInterfaceNodesMap()
    {
        return localToGlobalInterfaceNodesMap;
    }

    RDouble *GetPMin()
    {
        return minBoxOfZone;
    }
    RDouble *GetPMax()
    {
        return maxBoxOfZone;
    }
    RDouble *GetNodeWalldistDistance()
    {
        return nodesMinDistanceInOtherBlock;
    }
    int GetBlockIndex() {
        return blockIndex;
    }
    int GetNumberOfInterCells()
    {
        return numberOfInterCells;
    }
    //vector< int > &GetInterCellIndex() { return interpolateCellList; }
    vector< int > &GetInterpolateCellList()
    {
        return this->interpolateCellList;
    }

    OversetInformationProxy *GetOversetInformationProxy()
    {
        return oversetInformationProxy;
    }
    UnstructGrid *GetGrid()
    {
        return grid;
    }
private:
    inline void OutPutValueName(fstream &file, const string &name)
    {
        file << "# " << name << endl;
    }

    inline void OutPutNullValue(fstream &file, int length)
    {
        int numberOfWordsInEachLine = 5;
        int iPrint = 0;
        for (int iValue = 0; iValue < length; ++ iValue)
        {
            file << 0 << " ";
            ++ iPrint;
            if (iPrint % numberOfWordsInEachLine == 0) file << endl;
        }
        if (iPrint % numberOfWordsInEachLine != 0) file << endl;
        file << endl;
    }

    inline void OutPutVectorValue(fstream &file, int length, vector< int > &value)
    {
        bool *valueFull = new bool[length];
        SetField(valueFull, false, length);
        int numberOfValue = value.size();
        for (int iValue = 0; iValue < numberOfValue; ++ iValue)
        {
            valueFull[value[iValue]] = true;
        }

        int numberOfWordsInEachLine = 5;
        int iPrint = 0;
        for (int iValue = 0; iValue < length; ++ iValue)
        {
            file << valueFull[iValue] << " ";
            ++ iPrint;
            if (iPrint % numberOfWordsInEachLine == 0) file << endl;
        }
        if (iPrint % numberOfWordsInEachLine != 0) file << endl;
        file << endl;

        delete[] valueFull;
    }

    template < typename T >
    inline void OutPutPointerValue(fstream &file, int length, T *value)
    {
        if (NULL == value)
        {
            OutPutNullValue(file, length);
        }
        else
        {
            int numberOfWordsInEachLine = 5;
            int iPrint = 0;
            for (int iValue = 0; iValue < length; ++ iValue)
            {
                file << value[iValue] << " ";
                ++ iPrint;
                if (iPrint % numberOfWordsInEachLine == 0) file << endl;
            }
            if (iPrint % numberOfWordsInEachLine != 0) file << endl;
            file << endl;
        }
    }
};

RDouble GetToleranceForOversetBox();
int GetTwoOrderInterpolationOrNot();
int GetKeyEnlargeOfNodes();
int GetReadInAuxiliaryInnerGrid();
int GetReadInAuxiliaryOuterGrid();
int GetReadInSklFileOrNot();
int GetSymetryOrNot();
int GetIfOutPutOversetVisualization();

template< typename T >
void PH_AllreduceSepModeForVector(vector< T > &sendingBuffer, vector< T > &receivingBuffer, int numberOfElements, MPI_Datatype mpiDataType, PH_Op op);

}