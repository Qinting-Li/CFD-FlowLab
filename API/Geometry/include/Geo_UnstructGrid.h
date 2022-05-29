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
//! @file      Geo_UnstructGrid.h
//! @brief     It defines the class 'UnstructGrid', which is type of Unstructured geometry grid.
//!            The inheriting order is: SimpleGrid->Grid->StructuredGrid/UnstructuredGrid
//! @author    Bell, He Xin.

#pragma once
#include "PHHeader.h"
#include "Geo_Grid.h"
#include "Geo_CellTopo_Unstruct.h"
#include "Geo_NodeTopo_Unstruct.h"
#include "Geo_FaceMetrics_Unstruct.h"
#include "Geo_FaceTopo_Unstruct.h"
#include "Geo_CellMetrics_Unstruct.h"
#include "Geo_DynamicGridMetrics_Unstruct.h"
#include "Geo_LSQWeight_Unstruct.h"
#include "Geo_UnstructBC.h"
#include "Geo_LineTopo_Unstruct.h"

#ifdef USE_IncompressibleSolver
#include "PHVector3.h"
#endif

namespace PHSPACE
{
class Connectivity
{
public:
    int *data;
    int *index;
    int n;
public:
    Connectivity();
    ~Connectivity();
};

class UnstructBCSet;
class UnstructBC;
class VirtualFile;
class OversetConfig;
class GridManager;

class UnstructGrid : public Grid
{
private:
    //! The boundary condition list.\n
    //! Each face has a Boundary Condition Record.
    UnstructBCSet ** bcr;

    UnstructBCSet *unstructBCSet;

    //! The volume condition list, without using at the moment.\n
    //UnstructBCSet **vcr;

    //! Unstructured node topology information:\n
    //! Defines the node-face, node-cell, node-node \n
    //! connected relationship of unstructured grid.
    //! 1: cells index list of each node.
    Geo_NodeTopo_Unstruct * nodeTopology;

    //! Unstructured face topology information:\n
    //! Defines the face-node connected relationship of unstructured grid.\n
    //! 1: right and left cell of each face.\n
    //! 2: nodes list of each face.
    Geo_FaceTopo_Unstruct * faceTopology;

    Geo_LineTopo_Unstruct * lineTopology;

    //! Unstructured cell topology information:\n
    //! Defines the cell-face, cell-node connected relationship of unstructured grid.\n
    //! 1: faces of each cell.\n
    //! 2: nodes of each cell.\n
    //! 3: neighbor cells of each cell.
    Geo_CellTopo_Unstruct * cellTopology;

    //! Unstructured face metrics information:\n
    //! Defines the face metrics of unstructured grid.\n
    //! 1: face center.\n
    //! 2: face normal.\n
    //! 3: face area.
    Geo_FaceMetrics_Unstruct * faceMetrics;

    //! Unstructured cell metrics information:\n
    //! Defines the cell metrics of unstructured grid.\n
    //! 1: cell center.\n
    //! 2: cell volume.
    Geo_CellMetrics_Unstruct * cellMetrics;

    //! Unstructured dynamic grid metrics information:\n
    //! Defines dynamic mesh metrics of unstructured grid.\n
    //! 1: face normal velocity.\n
    //! 2: face velocity.
    Geo_DynamicGridMetrics_Unstruct * dynamicGridMetrics;

    //! Least Square Weight of unstructured grid.\n
    Geo_LSQWeight_Unstruct * leastSquareWeights;

    //! Connection of cell index to coarse cell index.\n
    //! For multi-grid method, the original finest grid is coarsen,
    //! so, each cell belongs to a coarse cell on the coarse level.
    int * cell2coarsegridcell;

    //! There are two cells for each cell face.\n
    //! fc2cL: the order of the right cell on the neighbor cell list.
    //! fc2cR: the order of the left  cell on the neighbor cell list.
    int * fc2cL, * fc2cR;

    //! Wall distance: the nearest distance of a cell to the wall.
    RDouble  * walldist;

    //! Wall distance: the nearest distance of a node to the wall.
    RDouble * walldistNode;

    RDouble* normalDistanceC2C;

    //! the label for first layer cell of the wall boundary.
    int * wallCellLabel;

    //! Cell to cell connectivity, used in over-set unstructured grid ONLY.
    Connectivity *cell2cell_connectivity;

    //! Some special variables for Laplacian-ggnode reconstruction
    RDouble    * lamdax, * lamday, * lamdaz; 
    int        * knode;

    //! A lable to judge if the cell is in the over-set zone.\n
    //! For over-set grid method only.
    int *iBlank;

    //! A lable to judge if the point is boundary point .\n;
    //! if point i is boundary point, boundaryPointLabel[i] = 1; else boundaryPointLabel[i] = 0;
    //! if point j is interpoint but not boundary point, boundaryPointLabel[j] = 0;
    map <int, int> boundaryPointLabel;

    //! Boundary condition type array.
    vector< int > boundaryConditionType;

    vector< int > localToGlobalNodesMap;

    vector< RDouble > nodeWalldistDistance;

    OversetConfig * oversetConfig;

    //! Average volume of all mesh.
    RDouble averageVolume;

    GridManager *gridManager;
    //! Probes' cell index of current zone grid.
    vector <int> zoneProbesCellID;


public:

#ifdef USE_IncompressibleSolver
    //! This is values for SIMPLE algorithm.
    //! cellCenterVector container the three coordinates of the cells center. 
    Vector3D *cellCenterVector;

    //! faceNormalVector container the three units of the face normal directions. 
    Vector3D *faceNormalVector;
#endif

public:
    LIB_EXPORT UnstructGrid();
    LIB_EXPORT ~UnstructGrid();

public:
    GridManager * GetGridManager();
    //! Set the boundary condition records.
    void SetBCRecord(UnstructBCSet **bcr);

    //! Set the boundary condition region.
    void SetUnstructBCSet(UnstructBCSet* unstructBCSet);

    //! Set the left cell list of each face.
    void SetLeftCellOfFace(int *left_cell_of_face_in);

    //! Set the right cell list of each face.
    void SetRightCellOfFace(int *right_cell_of_face_in); 

    //! Set the number of nodes per face.
    void SetNodeNumberOfEachFace(int *node_number_of_each_face);

    //! Set the face to node, node index of each face.
    void SetFace2Node(int *face2node);

    //! Set each face's starting subscript in the face2node array.
    void SetFace2NodeSubscript(int *face2nodeSubscript);

    //! Set the Connection of cell index to coarse cell index.
    void SetCell2CoarseGridCell(int *cell2coarsegridcell);

    //! set the number of lines.
    void SetLI_nLine(int nLine);

    //! set the number of cells in all lines.
    void SetLI_nCellsInLine(int* nCellsInLines);

    //! set the cell No. of all lines.
    void SetLI_CellOfLine(int** CellOfLine);

    //! set the face No. of all lines.
    void SetLI_FaceOfLine(int** FaceOfLine);

    //! set the line No. of all cells.
    void SetLI_LineOfCell(int* LineOfCell);

    //! Set the order of the right cell on the neighbor cell list.
    void Setfc2cL(int *fc2cL);

    //! Set the order of the left  cell on the neighbor cell list.
    void Setfc2cR(int *fc2cR);

    void SetBlankIndex(int *iBlank);

    //! Set the probes cell index of each zone grid.
    void SetZoneProbesCellID(vector <int> zoneProbesCellID);

    //! Return the boundary condition records.
    UnstructBCSet **GetBCRecord() const;
    
    //! Create composite boundary condition region (allocate memory for structBCSet).
    LIB_EXPORT void CreateUnstructBCSet(int nBCRegion);

    //! Return composite boundary condition region pointer (structBCSet).
    UnstructBCSet * GetUnstructBCSet() const;
   
    RDouble * GetNormalDistanceOfC2C() const;


    vector <int> & GetBoundaryConditionType();

    vector <int> & GetLocalToGlobalNodesMap();

    //! Get the left cell list of each face.
    int * GetLeftCellOfFace() const;

    //! Get the right cell list of each face.
    int * GetRightCellOfFace() const;

    //! Get the number of nodes per face.
    int * GetNodeNumberOfEachFace() const;

    //! Get the face to node, node index of each face.
    int * GetFace2Node() const;

    //! Return the X pointer of face center.
    RDouble * GetFaceCenterX() const;

    //! Return the Y pointer of face center.
    RDouble * GetFaceCenterY() const;

    //! Return the Z pointer of face center.
    RDouble * GetFaceCenterZ() const;

    //! Return the X pointer of face UNIT normal.
    RDouble * GetFaceNormalX() const;

    //! Return the Y pointer of face UNIT normal.
    RDouble * GetFaceNormalY() const;

    //! Return the Z pointer of face UNIT normal.
    RDouble * GetFaceNormalZ() const;

    //! Return the face area.
    RDouble * GetFaceArea()    const;

    //! Return the X pointer of cell center.
    RDouble * GetCellCenterX() const;

    //! Return the Y pointer of cell center.
    RDouble * GetCellCenterY() const;

    //! Return the Z pointer of cell center.
    RDouble * GetCellCenterZ() const;

    //! Return the cell volume pointer.
    RDouble * GetCellVolume()  const;

    //! Return the X velocity pointer of face center.
    RDouble * GetFaceVelocityX() const;

    //! Return the Y velocity pointer of face center.
    RDouble * GetFaceVelocityY() const;

    //! Return the Z velocity pointer of face center.
    RDouble * GetFaceVelocityZ() const;

    //! Return the face normal velocity pointer.
    RDouble * GetFaceNormalVelocity() const;

    //! Return the cell volume of the last time step.
    RDouble * GetCellVolumeOld() const;

    //! Return the connection of cell index to coarse cell index.
    int * GetCell2CoarseGridCell() const;

    //! Return the order of the right cell on the neighbor cell list.
    int * Getfc2cL() const;

    //! Return the order of the left cell on the neighbor cell list.
    int * Getfc2cR() const;

    LIB_EXPORT void InitGrid(GridID *index, int level, int dim, int type);

    //! Get the Number of faces per cell.
    LIB_EXPORT int * GetFaceNumberOfEachCell();

    //! Get the cell to face: face index of each cell.
    LIB_EXPORT int ** GetCell2Face();

    //! Get Number of nodes per cell.
    LIB_EXPORT int * GetNodeNumberOfEachCell();

    LIB_EXPORT Geo_CellTopo_Unstruct * GetCellTopology() { return cellTopology; }

    LIB_EXPORT Geo_LineTopo_Unstruct * GetLineTopology() { return lineTopology; }

    LIB_EXPORT bool ExistCell2Node();

    //! Get the cell to node: node index of each cell.
    LIB_EXPORT int * GetCell2Node();
    //! Get the cell to node: node index of each cell.
    LIB_EXPORT int ** GetCell2NodeArray();

    //! Get the face to node, node index of each face.
    LIB_EXPORT int ** GetFace2NodeArray();

    //! Get each face's starting subscript in the face2node array.
    LIB_EXPORT int * GetFace2NodeSubscript();

    //! get the number of lines.
    int GetLI_nLine() const;

    //! get the number of cells in all lines.
    int* GetLI_nCellsInLine() const;

    //! get the cell No. of all lines.
    int** GetLI_CellOfLine() const;

    //! get the face No. of all lines.
    int** GetLI_FaceOfLine() const;

    //! get the line No. of all cells.
    int* GetLI_LineOfCell() const;

    //! Get the cell to node: node index of each cell.
    LIB_EXPORT bool IsCell2NodeExist();

    //! Get Number of neighbor cells per cell.
    LIB_EXPORT int * GetNumberOfNeighborCell();

    //! Get the Cell to cell: neighbor cell index of each cell.
    LIB_EXPORT vector<int> * GetCell2Cell();

    //! Get the cell to cell connectivity, used in over-set unstructured grid ONLY.
    LIB_EXPORT Connectivity * GetCell2CellConnectivity();

    //! Get the Number of connected cell per node.
    LIB_EXPORT int * GetCellNumberOfEachNode();

    //! Get Node to cell: connected cells index of per node.
    LIB_EXPORT int * GetNode2Cell();

    //! Get two-dimensional array of Node to cell connectivity.
    LIB_EXPORT int ** GetNode2CellArray();

    //! Compute the minimum distance of the grid edges.
    LIB_EXPORT RDouble CalMinEdgeLength();

    //! Compute the maximum distance of the grid edges.
    LIB_EXPORT RDouble CalMaxEdgeLength();

    //! Reset the left&&right cell index to -1 on BC faces.
    LIB_EXPORT void ReSpecifyBC();

    //! Get the largest local grid length, distance of neighbor cell centers.
    LIB_EXPORT RDouble * GetLargestLocalGridLength();

    //! Get the sub-grid length, used for DES simulation ONLY.
    LIB_EXPORT RDouble * GetSubgridLength();

    //! Get the boundary point label used for judge the point type.
    //! @param[out] boundaryPointLabel    the boundary point.\n.
    map<int,int> & GetBoundaryPointLabel();

    //! compute the boundary point label used for judge the point type.
    void computeBoundaryPointLabel();

    //! Set the cell to node: node index of each cell.
    void SetCell2Node(int * cell2node);

    //! Set the number of nodes in each cell.
    void SetNodeNumberOfEachCell(int * node_number_of_each_cell);

    //! Get the number of nodes in each cell.
    int * GetNodeNumberOfEachCell() const;

    //! Set the least square weights.
    void SetLeastSquareIWT(RDouble *iwt);

    //! Set the least square weights IXX.
    void SetLeastSquareIXX(RDouble *ixx);

    //! Set the least square weights IYY.
    void SetLeastSquareIYY(RDouble *iyy);

    //! Set the least square weights IZZ.
    void SetLeastSquareIZZ(RDouble *izz);

    //! Set the least square weights IXY.
    void SetLeastSquareIXY(RDouble *ixy);

    //! Set the least square weights IXZ.
    void SetLeastSquareIXZ(RDouble *ixz);

    //! Set the least square weights IYZ.
    void SetLeastSquareIYZ(RDouble *iyz);

    //! Set the face mark of the least square method.
    void SetFaceMark(char *fMark);

    //! Return the least square weights.
    RDouble * GetLeastSquareIWT() const;

    //! Return the least square weights IXX.
    RDouble * GetLeastSquareIXX() const;

    //! Return the least square weights IYY.
    RDouble * GetLeastSquareIYY() const;

    //! Return the least square weights IZZ.
    RDouble * GetLeastSquareIZZ() const;

    //! Return the least square weights IXY.
    RDouble * GetLeastSquareIXY() const;

    //! Return the least square weights IXZ.
    RDouble * GetLeastSquareIXZ() const;

    //! Return the least square weights IYZ.
    RDouble * GetLeastSquareIYZ() const;

    //! Return the face mark of the least square method.
    char * GetFaceMark() const;

    //! Return the wall distance: the nearest distance of a cell to the wall.
    RDouble * GetWallDist() const;

    //! Return the wall distance: the nearest distance of a node to the wall.
    RDouble* GetWallDistNode() const;

    //! Return the wall cell label for the first layer of wall boundary. 
    int * GetWallCellLabel();

    void InitWallDist();

    vector< RDouble > & GetNodeWallDist();

    //! Computes and returns the number of faces on the solid wall.
    int GetNumberOfWallCell();

    //! Get the probes cell index of each zone grid.
    vector <int> GetZoneProbesCellID() const;

    int * GetBlankIndex() const;

    RDouble * GetLamdax() const;

    RDouble * GetLamday() const;

    RDouble * GetLamdaz() const;

    int * Getknode() const; 

    RDouble *GetVolume(int istep = 0) const;

    //! Find out the points that lie on the solid wall.
    LIB_EXPORT void FindOutWallPoints(bool * isWallPoint, int & nWallPoints);

    //! Find out the faces that lie on the solid wall.
    LIB_EXPORT void FindOutWallFaces(bool * isWallFace, int & nWallFaces);

    //! Static the mesh skewness.\n
    //! Skewness, 0 represents bad mesh quality; 1 represents satisfactory quality.
    //! Generally, the minimum skewness must be greater than 0.001.
    LIB_EXPORT RDouble SkewnessSummary(ActionKey *actkey = 0);

    //! Write face boundary condition to file stream.
    LIB_EXPORT void WriteFaceBC(fstream & file);

    //! Allocate iBlank for unstructured grid.
    LIB_EXPORT void AllocateOversetGrid();

    //! Allocate overset storage for donor cells and inter cells.
    LIB_EXPORT void RegisterOversetField(const string & name, const int type, const int dimesion);

    //! Set the average volume.
    void SetAverageVolume(RDouble Rin);

    //! Get the average volume.
    LIB_EXPORT RDouble GetAverageVolume();

    void GetSourceIndex(int iFace, int ipos, int &s);
    void GetTargetIndex(int iFace, int ipos, int &t);

    //! Return the global point number map to the interpoint for send.
    //! @param[in] ipoint : the interpoint number.
    //! @param[in] ipos : the number of ghost cell layers.
    //! @param[out] s :the global point number
    void GetSourcePointIndex(int ipoint, int ipos, int &s);

    //! Return the global point number map to the interpoint for receive.
    //! @param[in] ipoint : the interpoint number.
    //! @param[in] ipos : the number of ghost cell layers.
    //! @param[out] s :the global point number
    void GetTargetPointIndex(int ipoint, int ipos, int &s);

    void Action(ActionKey *actkey);
    void TranslateAction(ActionKey *actkey);
    streamsize TranslateActionLength(ActionKey *actkey);

    void SetGhostCellExceptInterface(RDouble *f);
    void SetGhostCell(RDouble *f);
    void SetGhostCell(RDouble **f, const int &nVariables);
    void SetGhostCell(RDouble *f, RDouble value);

    void AllocateWalldist();

    void SpecifyRightCellofBC();

    //! Return the gradient computed by ggnode.
    //! @param[in] qnode : the flow value on the point.
    //! @param[out] dqdx/dqdy/dqdz :the gradient value in three direction.
    void CompGradientGGNode_NEW(RDouble *q, RDouble *qnode, RDouble *dqdx, RDouble *dqdy, RDouble *dqdz);
    void CompGradientGGNodeWeight(RDouble *q, RDouble *dqdx, RDouble *dqdy, RDouble *dqdz);
    void CompGradientGGCell(RDouble *q, RDouble *dqdx, RDouble *dqdy, RDouble *dqdz);
    void CompGradientGGModified2(RDouble *q, RDouble *dqdx, RDouble *dqdy, RDouble *dqdz);
    void CompGradientGGCellNew(RDouble *q, RDouble *dqdx, RDouble *dqdy, RDouble *dqdz);
    void CompGradientGGCellW(RDouble *q, RDouble *dqdx, RDouble *dqdy, RDouble *dqdz);
    void CompGradientLSQ   (RDouble *q, RDouble *dqdx, RDouble *dqdy, RDouble *dqdz);
    void CompGradientGGNode(RDouble *q, RDouble *dqdx, RDouble *dqdy, RDouble *dqdz);

    void CompGradient(RDouble *q, RDouble *dqdx, RDouble *dqdy, RDouble *dqdz);

    void CalcLaplacianWeitht();
    void FaceWeight(RDouble *deltl, RDouble *deltr, int nst, int ned);
    void CVGNorm(FieldProxy *q_1_proxy, FieldProxy *q_0_proxy, int neqn, RDouble &norm);

    void ComputeWeight();

    //! Compute node bctype in each zone.
    LIB_EXPORT void ComputeNodeBCType();

    void GetMinMaxDS(RDouble &dismin, RDouble &dismax);
    void GetBoundaryMinMax(RDouble *pmin, RDouble *pmax);

    // calculate the normal distance between two neighbor cells.
    LIB_EXPORT void CalNormalDistanceOfC2C();
    #ifdef USE_INHOUSE_VERSION
    LIB_EXPORT void ReorderCellForLineLUSGS();
    #endif
    void ReadGrid (fstream &file);
    void ReadGrid (VirtualFile *vfile);
    LIB_EXPORT void WriteGrid(fstream &file);
    void WriteGrid(VirtualFile *vfile);

    void Decode(DataContainer *cdata, int flag);

    //! Revise the zone index of neighbor zones.
    void ReviseNeighborZoneIndex(int zoneStart);

    void DecodeOversetGrid(DataContainer *cdata);
    void DecodeCellNode(DataContainer * cdata);

    //! Return the interpoint information from the DataContainer.
    void DecodeInterpointInfo(DataContainer *cdata);

    void DecodeWallDist(DataContainer * cdata);
    void DecodeBCFace(DataContainer * cdata);
    void Encode(DataContainer *cdata, int flag);
    void EncodeGrid(DataContainer *cdata);
    void EncodeOversetGrid(DataContainer *cdata);

    void UpdateCoordinate(DataContainer *cdata);

    void DXDYDZ(RDouble *f, int ie, RDouble &dfdx, RDouble &dfdy, RDouble &dfdz);
    void DXDYDZ(RDouble *f, int ie, int **cell2face, int *face_number_of_each_cell, RDouble &dfdx, RDouble &dfdy, RDouble &dfdz);
    void DXDYDZ_Face(RDouble *f, int iFace, int nodepos, RDouble fl, RDouble fr, RDouble &dfdx, RDouble &dfdy, RDouble &dfdz);
    void DXDYDZ_GG (RDouble *f, int ie, RDouble &dfdx, RDouble &dfdy, RDouble &dfdz);
    void DXDYDZ_LSQ(RDouble *f, int ie, RDouble &dfdx, RDouble &dfdy, RDouble &dfdz);

    void ComputeMetrics(ActionKey *actkey = 0);

    LIB_EXPORT void CoarseGrids(int maxLevel);
    void BuildGridLink(LinkStruct * link);
    void CalCNNCF();
    void ChangeBCType(int from_bctype, int to_bctype);
    void CompareMaxMinValue(RDouble &local_data, int flag);
    void ProcessBCInfo();
    void SetLnkInfo() {};
    void GridSurfaceVelocity(RDouble *xt, RDouble *yt, RDouble *zt);

    LIB_EXPORT void ComputeGridFaceVelocity();

    LIB_EXPORT void BackUpOldGrid();

    void ReSetBoundaryConditionByGlobalBC();

    void SetOversetConfig(OversetConfig *oversetConfig);

    void AllocateMetrics(ActionKey *actkey);
    void AllocateMetricsALE(ActionKey *actkey);

    void GetBcNamePair(set< pair<int, string> > &bcNamePairSet);
    void DumpWallFaceCenter(ActionKey *actkey);
    void InitVariableWallTemperature();
    void GetCellIBlank(ActionKey *actkey);
    void TranslateCellIBlank(ActionKey *actkey);
    streamsize TranslateCellIBlankLength(ActionKey *actkey);

private:
    int  CompNIFace();
    void ComputeCell2Face();
    void ComputeCell2Node();
    void ComputeCell2NodeArray();
    void ComputeFace2NodeArray();
    void ComputeFace2NodeSubscript();
    void ComputeNumberOfNeighborCell();
    void ComputeCell2Cell();
    void ComputeCell2CellConnectivity();
    void ComputeNode2Cell();
    void ComputeNode2CellArray();
    void ComputeCell2Node(set<int> *&cell2node_set);

    void ComputeLargestLocalGridLength();
    void ComputeSubgridLength();

    void DumpWalldist(ActionKey *actkey);
    void ReadWalldist(RDouble *wallDistIn);

    void GetCellCenter(ActionKey *actkey);
    void TranslateCellCenter(ActionKey *actkey);
    streamsize TranslateCellCenterLength(ActionKey *actkey);

    void SimpleAction(ActionKey *actkey);
    void UpdateVolold();
    void ClosureCheck(RDouble *xfn, RDouble *area);

    void CompFaceGradient(RDouble *q, RDouble *dqdx, RDouble *dqdy, RDouble *dqdz);
    void CompFaceGradientGauss2D(RDouble *q, RDouble *qn, RDouble *dqdx, RDouble *dqdy, RDouble *dqdz);
    void CompFaceGradientGauss  (RDouble *q, RDouble *qn, RDouble *dqdx, RDouble *dqdy, RDouble *dqdz);

    void ReadMovingGrid();
    void DecodeGrid(DataContainer *cdata);

    void ComputeFaceNormal2D();
    void ComputeFaceNormal3D();
    void ComputeFaceNormal3DNew();

    void ComputeFaceCenter2D(ActionKey *actkey);
    void ComputeFaceCenter3D(ActionKey *actkey);
    void ComputeFaceCenter3DNew(ActionKey *actkey);

    void ComputeCellCenterVol2D(ActionKey *actkey);
    void ComputeCellCenterVol3D(ActionKey *actkey);
    void ComputeCellCenterVol3DNew(ActionKey *actkey);
    
    void ComputeFaceCenterAreaNormal3DNew2(ActionKey *actkey);
    void ComputeCellCenterVol3DNew2(ActionKey *actkey);

public:
#ifdef USE_IncompressibleSolver
    //! This is function for SIMPLE algorithm.
    Vector3D * GetCellCenterVector() const;
    Vector3D * GetFaceNormalVector() const;
    void SetCellCenterVector(Vector3D *cellC);
    void SetFaceNormalVector(Vector3D *faceN);
#endif
};

#include "Geo_UnstructGrid.hxx"

void CoarseGridConnectUNMergeFace(UnstructGrid *fgrid, UnstructGrid *cgrid);

int PickOutBadSkewnessFaces(UnstructGrid * grid, bool * isBadFace = 0, bool * isBadCell = 0);

void FieldVisualization(Grid *grid_in, std::ostringstream &oss, vector<string> &title, RDouble **qn, int nl);
void FieldVisualizationForVTK(Grid *grid_in, std::ostringstream &oss, vector<string> &title, RDouble **qn, int nl);
void SaveDataForTecio_bac(Grid *grid_in, ActionKey *actkey, RDouble **qn, int nl);
void SaveDataForTecio(Grid *grid_in, ActionKey *actkey, RDouble **qn, int nl);

void CompNodeVarWeight(UnstructGrid *grid, RDouble *q_n, RDouble *q);    // Bell 20121220 add
void CompNodeVar(UnstructGrid *grid, RDouble *q_n, RDouble *q, bool isVelocityForPostVisual = false);
void CompNodeVar_new(UnstructGrid *grid, RDouble *q_n, RDouble *q);
void CompNodeVarForVisual(UnstructGrid *grid, RDouble *q_n, RDouble *q);
void CompNodeVarForVelocity(UnstructGrid *grid, RDouble *qNode, RDouble *q);
void ComputeNodeVariable(UnstructGrid *grid, RDouble **q_n, RDouble **q, int nEquation);
void FixBCNodeVar(UnstructGrid *grid, RDouble *q, RDouble *q_n, int *n_count, int bctype_in, bool twoside);
void FixBCNodeVar(UnstructGrid *grid, RDouble *q, RDouble *q_n, RDouble *n_count, int bctype_in, bool twoside);    // Bell 20121220 add
void CompNodeVarLimit(UnstructGrid *grid, RDouble *q_n, RDouble *q);

void GetFace2NodeList(UnstructGrid *grid, pair<int, string> iter, vector<int> &linkmap, vector < vector < int > > &face2nodelist);
void GetFace2NodeList(UnstructGrid *grid, int bcType, vector<int> &linkmap, vector < vector < int > > &face2nodelist);

void MinMaxDiffINF(Grid *grid_in, int ig, RDouble *dmin, RDouble *dmax, RDouble *q);

//! Communicate the interface variables at the boundary faces between different zones.
//! @param[in] grid              computational mesh.
//! @param[in] variable          variable pointer of which need to be communicated.
//! @param[in] variableName      variable name of which need to be communicated.
//! @param[in] variableDimension variable dimension of which need to be communicated.
//! -# 1 :  1-D array.
//! -# 2 :  2-D array.
//! -# n :  n-D array.
void CommunicateInterfaceValue(UnstructGrid * grid, RDouble ** variable, const string & variableName, int variableDimension);
void CommunicateInterfaceValue(UnstructGrid *grid, RDouble *variable, const string &variableName);

//! Communicate the interpoint variables at the interface between different zones.
//! @param[in] grid                computational mesh.
//! @param[in] variable            variable pointer of which need to be communicated.
//! @param[in] variableName        variable name of which need to be communicated.
//! @param[in] variableDimension   variable dimension of which need to be communicated.
//! -# 1 :  1-D array.
//! -# 2 :  2-D array.
//! -# n :  n-D array.
void CommunicateInterpointValue(UnstructGrid *grid, RDouble **variable, const string &variableName, int variableDimension);

void UploadInterfaceValue(UnstructGrid *grid, RDouble **f, const string &name, int neqn);
void UploadInterfaceValue(UnstructGrid *grid, RDouble *f, const string &name);
void DownloadInterfaceValue(UnstructGrid *grid, RDouble **f, const string &name, int neqn);
void DownloadInterfaceValue(UnstructGrid *grid, RDouble *f, const string &name);

//! Upload the interpoint variables for send.
//! @param[in] grid                computational mesh.
//! @param[in] f                   variable pointer of which need to be sent.
//! @param[in] name                variable name of which need to be sent.
//! @param[in] neqn                variable dimension of which need to be sent.
//! -# 1 :  1-D array.
//! -# 2 :  2-D array.
//! -# n :  n-D array.
void UploadInterpointValue(UnstructGrid *grid, RDouble **f, const string &name, int neqn);

//! down the interpoint variables for receive.
//! @param[in] grid    computational mesh.
//! @param[in] f       variable pointer of which need to be received.
//! @param[in] name    variable name of which need to be received.
//! @param[in] neqn    variable dimension of which need to be received.
//! -# 1 :  1-D array.
//! -# 2 :  2-D array.
//! -# n :  n-D array.
void DownloadInterpointValue(UnstructGrid *grid, RDouble **f, const string &name, int neqn);

void UploadOversetValue(UnstructGrid *grid, RDouble **field, const string &name, int neqn);
void UploadOversetValueDefault(UnstructGrid *grid, RDouble **field, const string &name, int neqn);
void UploadOversetValueByWeight(UnstructGrid *grid, RDouble **field, const string &name, int neqn);
void DownloadOversetValue(UnstructGrid *grid, RDouble **field, const string &name, int neqn);

void BarthLimiterINF (Grid *grid_in, int ig, RDouble *limit, RDouble *q, RDouble *dqdx, RDouble *dqdy, RDouble *dqdz);
void VencatLimiterINF(Grid *grid_in, int ig, RDouble *limit, RDouble *q, RDouble *dqdx, RDouble *dqdy, RDouble *dqdz);

void SmoothField(UnstructGrid *grid, RDouble *q);

void SetLimitBoundary(Grid *grid_in, RDouble *limit);
void addn2c(int ie, int ip, int *nCPN, struct linkbase **nc);
UnstructGrid * UnstructGridCast(Grid *grid_in);

LIB_EXPORT void ConstructFace2Node(int * nPointPerFace, int * nodesOfPerFace, uint_t nFaces, vector< vector<int> > & face2node);

void VisualizationMesh2D(UnstructGrid *grid, const string &fileName);
void VisualizationMesh3D(UnstructGrid *grid, const string &fileName);
void VisualizationEigenMeshofCoarseGrid3D(UnstructGrid *grid, const string &fileName);

void GetFaceToNodeList(UnstructGrid *grid, int boundaryConditionType, PHVectorInt1D &linkMap, PHVectorInt2D &faceToNodeList);

void ExtractSubGridFromGrid(UnstructGrid *grid, UnstructGrid *grid_Slice, int *isThisCellNeedExtract);

void DumpSliceGrid(UnstructGrid *unstructGrid, int sliceAxis, RDouble slicePosition, string fileName);

RDouble AspectRatio(RDouble dist, RDouble area, int nodeNumOfPerFace);
}
