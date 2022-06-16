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
#pragma once
#include "Pre_GridBase.h"
#include "Geo_PointFactory.h"

namespace PHSPACE
{
class ComputationalFaceSolver
{
public:
    ComputationalFaceSolver ();
    ~ComputationalFaceSolver();
public:
    PHVectorInt1D & GetLeftElementIndexContainer () { return *leftElementIndexContainer;  }
    PHVectorInt1D & GetRightElementIndexContainer() { return *rightElementIndexContainer; }

    PHVectorInt1D & GetOldLeftElementIndexContainer () { return *oldLeftElementIndexContainer;  }
    PHVectorInt1D & GetOldRightElementIndexContainer() { return *oldRightElementIndexContainer; }
public:
    PHVectorInt1D & GetFaceTypeContainer() { return *faceTypeContainer; }
    PHVectorInt1D & GetFaceBoundaryTypeContainer() { return *faceBoundaryTypeContainer; }
    vector<SimpleBC * > & GetFaceBoundaryCondition() { return *faceBoundaryCondition; }
public:
    int  GetNumberOfFaces() const;
    void SetNumberOfFaces(int numberOfFaces); //Green 20130703

    int  GetNumberOfBoundaryFaces() const; //Green
    void SetNumberOfBoundaryFaces(int numberOfBoundaryFaces);

    int  GetNumberOfElements() const { return numberOfElements; }
    void SetNumberOfElements(int numberOfElements) { this->numberOfElements = numberOfElements; }
public:
    void GenerateGhostField(RFloat **qField, RFloat **ghostQField, int numberOfEquations);
public:
    PHVectorInt2D & GetFaceNodeIndexesContainer   () { return *faceNodeIndexesContainer;  } //
    PHVectorInt2D & GetOldFaceNodeIndexesContainer() { return *oldFaceNodeIndexesContainer; } //

    //ComputationalElement中face编号与对应Element中Nodes编号，过渡
    PHVectorInt2D & GetIntermediateFaceNodeIndexesContainer() { return *intermediateFaceNodeIndexesContainer; }
public:
    PHVectorInt1D & GetComputationalToGlobalFaceIndexMapping   () { return computationalToGlobalFaceIndexMapping;    } //Green
    PHVectorInt1D & GetOldComputationalToGlobalFaceIndexMapping() { return oldComputationalToGlobalFaceIndexMapping; } //Green
protected:
    PHVectorInt1D *leftElementIndexContainer;     //实际计算所用关系
    PHVectorInt1D *rightElementIndexContainer;    //实际计算所用关系

    PHVectorInt1D *oldLeftElementIndexContainer;  //上一次计算信息，用于重起计算
    PHVectorInt1D *oldRightElementIndexContainer; //上一次计算信息，用于重起计算

    PHVectorInt1D *faceBoundaryTypeContainer; //
    PHVectorInt1D *faceTypeContainer;         //

    vector<SimpleBC * > *faceBoundaryCondition;

    PHVectorInt2D *intermediateFaceNodeIndexesContainer; //ComputationalElement中face编号与对应Element中Nodes编号，过渡
    PHVectorInt2D *faceNodeIndexesContainer;             //ComputationalElement中face编号与对应ComputationalElement中Nodes编号，完全转化
    PHVectorInt2D *oldFaceNodeIndexesContainer;          //上一次计算ComputationalElement中face编号与对应ComputationalElement中Nodes编号

    PHVectorInt1D computationalToGlobalFaceIndexMapping;
    PHVectorInt1D oldComputationalToGlobalFaceIndexMapping; //上一次计算信息，用于重起计算
protected:
    int numberOfElements;
    int numberOfFaces;         // Green
    int numberOfBoundaryFaces; // Green
};

class ElementPropertySolver;

const int ISOTROPIC_REFINE   = 0;
const int ANISOTROPIC_REFINE = 1;


class ElementFaceSolver
{
public:
    ElementFaceSolver(Geo_PointFactory *pointFactory);
    ~ElementFaceSolver();
protected:
    set< MultiIndex > *referenceFacesContainer;
    set< MultiIndex > *compositeFacesContainer;        //用于棱中点查找
    set< MultiIndex > *edgesWithMiddlePointsContainer; //计算单元的face中具有中点的棱

    PHVectorInt2D *faceNodeIndexesContainer;            //单元中face编号和Nodes编号对应关系
    PHVectorInt2D *nodeIndexesOfCompositeFaceContainer; //单元中复合face编号和Nodes编号对应关系  Green 20130630

    PHVectorInt1D *leftElementIndexContainer;           //指单元中face与左右单元对应关系
    PHVectorInt1D *rightElementIndexContainer;

    PHVectorInt1D *faceBoundaryTypeContainer;
    PHVectorInt1D *faceTypeContainer;
    vector<SimpleBC * > *faceBoundaryCondition;
    PHVectorInt1D *parentFaceIndexContainer;

    PHVectorInt1D *compositeFaceTypeContainer;      //
    PHVectorInt1D *middlePointOfEdgeIndexContainer; //计算单元的具有中点的棱的中点编号    Green20130629

    PHVectorInt2D *childFaceIndexContainer;

    PHVectorInt1D globalToComputationalFaceIndexMapping;

    PHVectorInt1D simpleToCompositeFaceIndexMapping; // Green20130629
    PHVectorInt1D compositeToSimpleFaceIndexMapping; // Green20130629

    PHVectorRDouble1D faceArea;
    PHVectorRDouble1D faceRBFWeight;

    PHVectorRDouble1D faceNormalX;
    PHVectorRDouble1D faceNormalY;
    PHVectorRDouble1D faceNormalZ;

    PHVectorRDouble1D faceCenterX;
    PHVectorRDouble1D faceCenterY;
    PHVectorRDouble1D faceCenterZ;
protected:
    ElementPropertySolver *elementPropertySolver;
    Geo_PointFactory *pointFactory;
    ComputationalFaceSolver *computationalFaceSolver;
public:
    ComputationalFaceSolver * GetComputationalFaceSolver() { return computationalFaceSolver; }
    ElementPropertySolver * GetElementPropertySolver() { return elementPropertySolver; }
    void SetElementPropertySolver(ElementPropertySolver *elementPropertySolver) { this->elementPropertySolver = elementPropertySolver; }

    bool IsCompositeFace(const int &faceIndex) { return GetChildFaceIndexContainer()[ faceIndex ].size() != 0; }

    int  FindFace(MultiIndex &face);

    int  FindCompositeFace      (MultiIndex &compositeFace      ); //Green20130629
    int  FindEdgeWithMiddlePoint(MultiIndex &edgeWithMiddlePoint); //Green20130629

    int  FindFace(PHVectorInt1D &faceNodeIndexContainer);

    void InsertFace(MultiIndex &face);

    void InsertCompositeFace(MultiIndex &compositeface);
    void InsertEdgeWithMiddlePoint(MultiIndex &edgeWithMiddlePoint);

    void ResizeFaceNumber(int nSize);
    void AddFaceHierachical(int faceIndex, int parentFaceIndex);
    void SetDefaultFaceInformation(int faceIndex, int faceType, PHVectorInt1D &faceNodeIndexes, int parentFaceIndex); //Green

    int  GetNumberOfFaces();
    int  GetNumberOfBoundaryFaces(); //暂时没用
    int  GetNumberOfTotalFaces() const { return static_cast<int>(faceNodeIndexesContainer->size()); }

    int  GetNumberOfTotalCompositeFaces       () const { return static_cast<int>(nodeIndexesOfCompositeFaceContainer->size()); }
    int  GetNumberOfTotalEdgesWithMiddlePoints() const { return static_cast<int>(middlePointOfEdgeIndexContainer->size());     }

    PHVectorRDouble1D & GetFaceArea() { return faceArea; }
    PHVectorRDouble1D & GetFaceRBFWeight() { return faceRBFWeight; }
    PHVectorRDouble1D & GetFaceNormalX() { return faceNormalX; }
    PHVectorRDouble1D & GetFaceNormalY() { return faceNormalY; }
    PHVectorRDouble1D & GetFaceNormalZ() { return faceNormalZ; }

    PHVectorRDouble1D & GetFaceCenterX() { return faceCenterX; }
    PHVectorRDouble1D & GetFaceCenterY() { return faceCenterY; }
    PHVectorRDouble1D & GetFaceCenterZ() { return faceCenterZ; }

    PHVectorInt2D & GetFaceNodeIndexesContainer         () { return *faceNodeIndexesContainer;            } //
    PHVectorInt2D & GetCompositeFaceNodeIndexesContainer() { return *nodeIndexesOfCompositeFaceContainer; } //

    set < MultiIndex > & GetReferenceFacesContainer       () { return *referenceFacesContainer;        } //Green
    set < MultiIndex > & GetCompositeFacesContainer       () { return *compositeFacesContainer;        } //Green 20130629
    set < MultiIndex > & GetEdgesWithMiddlePointsContainer() { return *edgesWithMiddlePointsContainer; } //Green 20130629

    PHVectorInt1D & GetLeftElementIndexContainer () { return *leftElementIndexContainer;  }
    PHVectorInt1D & GetRightElementIndexContainer() { return *rightElementIndexContainer; }

    PHVectorInt1D & GetFaceBoundaryTypeContainer   () { return *faceBoundaryTypeContainer; }
    vector<SimpleBC * > & GetFaceBoundaryCondition () { return *faceBoundaryCondition; }

    PHVectorInt1D & GetFaceTypeContainer       () { return *faceTypeContainer;        } //Green
    PHVectorInt1D & GetParentFaceIndexContainer() { return *parentFaceIndexContainer; }

    PHVectorInt1D & GetCompositeFaceTypeContainer     () { return *compositeFaceTypeContainer; }
    PHVectorInt1D & GetMiddlePointOfEdgeIndexContainer() { return *middlePointOfEdgeIndexContainer; }

    PHVectorInt2D & GetChildFaceIndexContainer() { return *childFaceIndexContainer; }

    void ScanNewElementFace        (PHVectorInt1D &elementNodeIndexes, int elementType, int elementIndex);
    void ScanNewElementBoundaryFace(set< int > &bcVertex, int boundaryConditionType, int numberOfElementOfThisBCType);
    int  AddElementFace(PHVectorInt1D &faceNodeIndexContainer, int faceTypeContainer, int parentFaceIndex);
    void MarkAllCandidateBoundaryFace();
    bool CheckBCFace(set< int > &bcVertex, PHVectorInt1D &faceNodeIndexContainer);
    int  GetRootFaceIndex(int iFace);
    
    void GenerateComputationalElementFaceInformation(); //Green
    void GenerateGlobalToComputationalElementFaceIndexMapping(); //Green
    void GenerateComputationalFaceNodeIndex(); //Green
    void GenerateEdgesWithMiddlePoints(int iCompositeFace); //Green20130629
protected:
    void ProcessBoundaryFace(int iFaceForComputation);
    void ProcessGereralFace (int iFaceForComputation);
    void ProcessHangingFace (int iFaceForComputation);
    void GenerateComputationalElementFaceInformationFirstPart ();
    void GenerateComputationalElementFaceInformationSecondPart();
    void GenerateGlobalToComputationalElementFaceIndexMapping(int iFace, int &iComputationalFaceCount, int &iBoundaryFace);

    void GenerateGlobalToComputationalElementFaceIndexMappingNoParentCase  (int iFace, int &iComputationalFaceCount, int &iBoundaryFace);
    void GenerateGlobalToComputationalElementFaceIndexMappingWithParentCase(int iFace, int &iComputationalFaceCount, int &iBoundaryFace);
public:
    PHVectorInt1D & GetGlobalToComputationalFaceIndexMapping() { return globalToComputationalFaceIndexMapping; } //Green

    PHVectorInt1D & GetSimpleToCompositeFaceIndexMapping() { return simpleToCompositeFaceIndexMapping; }
    PHVectorInt1D & GetCompositeToSimpleFaceIndexMapping() { return compositeToSimpleFaceIndexMapping; }

    void GenerateComputationalUnstructuredGrid(UnstructGrid *grid);
    void SetFaceBoundary(PHVectorInt1D &faceNodeIndexContainer, int faceType, int boundaryConditionType);
public:
    void DumpAdaptiveFaces                      (DataContainer *dataContainer);
    void WriteCompositeFacesContainer           (DataContainer *dataContainer);
    void WriteCompositeFaceToSimpleFaceContainer(DataContainer *dataContainer);

    void ReorderFaceToNodeIndexes(UnstructGrid *grid, PHVectorInt1D &leftElementIndexContainerIn, PHVectorInt1D &rightElementIndexContainerIn, PHVectorInt2D &faceNodeIndexesContainerIn);
    void ReorderFaceToCellIndexes(UnstructGrid *grid, PHVectorInt1D &leftElementIndexContainerIn, PHVectorInt1D &rightElementIndexContainerIn);
};

class ElementPropertySolver
{
public:
    ElementPropertySolver();
   ~ElementPropertySolver();
private:
    int numberOfBoundaryFaces;
    int oldNumberOfElements;
    
protected:
    PHVectorInt1D *elementTypeContainer;                //单元类型
    PHVectorInt1D *elementLevelContainer;               //单元级别
    PHVectorInt1D *elementModifiedStatusContainer;      //单元改动:NOCHANGE;SPLIT,分裂或不分裂
    PHVectorInt1D *parentElementIndexContainer;         //单元parent
    PHVectorInt1D *elementComputationalStatusContainer; //单元状态:ON;HIDDEN，计算或不计算
    PHVectorInt2D  elementNodeIndexContainer;           //单元index

    PHVectorInt2D *element2Face;                        // Cell to face.
    PHVectorInt1D *refineType;                          // Cell refine type, isotropic or anisotropic.
    PHVectorInt1D *BoundaryLevel;                       // boundary level.
    PHVectorInt1D *projectType;                         // Is the wall face project or not.
    PHVectorInt1D *wallBCID;                            // BC name ID(name order).

    PHVectorInt1D *computationalElementModifiedStatusContainer; //记录根据流场，上一次计算网格自适应状态
    PHVectorInt2D *childrenElementsIndexesContainer;

    PHVectorInt1D globalToComputationalElementIndexMapping;   //实际选用体单元编号与ComputationalElement编号对应
    PHVectorInt1D computationalToGlobalElementIndexMapping;   //反变换

    PHVectorInt1D oldGlobalToComputationalElementIndexMapping;   //上一次计算对应单元信息，用于重起计算
    PHVectorInt1D oldComputationalToGlobalElementIndexMapping;   

    PHVectorRDouble1D oldComputationalElementVolumeVector;

    vector<bool> isSolidWallCell;
protected:
    ElementFaceSolver *elementFaceSolver;
public:
    vector<bool> & GetIsSolidWallCell() { return this->isSolidWallCell; }
    void SetElement2Face(PHVectorInt2D *cell2face) { this->element2Face = cell2face; }
    PHVectorInt2D * GetElement2Face() { return this->element2Face; }

    PHVectorInt1D * GetRefineType() { return this->refineType; }
    PHVectorInt1D * GetBoundaryLevel() { return this->BoundaryLevel; }
    PHVectorInt1D * GetProjectType() { return this->projectType; }
    PHVectorInt1D * GetWallBCID() { return this->wallBCID; }

    ElementFaceSolver * GetElementFaceSolver() { return elementFaceSolver; }
    void SetElementFaceSolver(ElementFaceSolver *elementFaceSolver) { this->elementFaceSolver = elementFaceSolver; }
    void DumpElementToNode(VirtualFile *virtualFile);

    PHVectorInt1D & GetElementTypeContainer               (){ return *elementTypeContainer               ; } //返回引用类型函数
    PHVectorInt1D & GetElementLevelContainer              (){ return *elementLevelContainer              ; } //Green
    PHVectorInt1D & GetElementModifiedStatusContainer     (){ return *elementModifiedStatusContainer     ; } //Green
    PHVectorInt1D & GetElementComputationalStatusContainer(){ return *elementComputationalStatusContainer; } 
    PHVectorInt1D & GetParentElementIndexContainer        (){ return *parentElementIndexContainer        ; }  
    PHVectorInt2D & GetChildrenElementsIndexesContainer   (){ return *childrenElementsIndexesContainer   ; } 
    PHVectorInt2D & GetElementNodeIndexContainer          (){ return  elementNodeIndexContainer          ; }

    PHVectorInt1D & GetOldComputationalToGlobalElementIndexMapping() { return oldComputationalToGlobalElementIndexMapping; }
    PHVectorInt1D & GetOldGlobalToComputationalElementIndexMapping() { return oldGlobalToComputationalElementIndexMapping; }

    PHVectorInt1D & GetComputationalToGlobalElementIndexMapping() { return computationalToGlobalElementIndexMapping; }
    PHVectorInt1D & GetGlobalToComputationalElementIndexMapping() { return globalToComputationalElementIndexMapping; }

    PHVectorInt1D & GetChildrenElementsIndexes(int elementIndex);
    PHVectorInt1D & GetComputationalElementModifiedStatusContainer(){ return *computationalElementModifiedStatusContainer; } //Green

    PHVectorRDouble1D & GetOldComputationalElementVolumeVector() { return oldComputationalElementVolumeVector; }

    void SetElementType(int elementIndex, int elementType);
    int  GetElementType(int elementIndex) const;

    void SetSigleElementNodeIndexContainer(int elementIndex, PHVectorInt1D &elementNodeIndexes);
    PHVectorInt1D & GetElementNodeIndexes(int elementIndex)  { return  elementNodeIndexContainer[elementIndex]; }

    void SetNumberOfBoundaryFaces(const int numberIn) { this->numberOfBoundaryFaces = numberIn; }
    int  GetNumberOfBoundaryFaces() { return this->numberOfBoundaryFaces; }
    
    void SetOldNumberOfElements(const int numberIn) { this->oldNumberOfElements = numberIn; }
    int  GetOldNumberOfElements() { return this->oldNumberOfElements; }

    int  GetNumberOfTotalElements() const;
    void ResizeElementNumber(int numberOfTotalElements);

    void SetParentElementIndex(int elementIndex, int parentElementIndex);
    int  GetParentElementIndex(int elementIndex);

    void SetElementLevel(int elementIndex, int iLevel);
    int  GetElementLevel(int elementIndex);

    void SetElementModifiedStatus(int elementIndex, int elementModifiedStatus);
    int  GetElementModifiedStatus(int elementIndex);

    void SetComputationalStatusOfElement(int elementIndex, int computationalStatusOfElement);
    int  GetComputationalStatusOfElement(int elementIndex);

    void SetAllElementsModifiedStatusToNoChange();
    void SetAllElementsComputationalStatusToOnState();
    void SetAllElementsComputationalStatusToHiddenState();
    void SetAllElementsParentIndexToRootState();
    void SetAllElementsLevelToZero();

    void InitializeOldElementIndexMapping();
    void InitializeElementIndexAndState();

    void SplitAllElements(Geo_PointFactory *pointFactory);
    void SplitElement(int elementIndex, Geo_PointFactory *pointFactory);

    bool ElementNeedSplit(int elementIndex);

    void SplitElementFace(int elementIndex);
    void SplitCompositeFace(PHVectorInt1D &nodeIndexesOfCompositeFace, int compositeFaceType);

    bool ElementNeedMerge(int elementIndex);
    void MergeElements();

    void AddElement(PHVectorInt1D &elementNodeIndexes, int elementType, int levelOfElement, int parentElementIndex);
    void AddChildElement(int parentElementIndex, Geo_PointFactory *pointFactory);
    void ReactivateChildElement(int elementIndex);

    void HideElement(int elementIndex);

    void ScanAllElementsToReconstructFaceDataStructure(const bool hasBeenAdapted);

    void GenerateCellToNode(UnstructGrid *grid);
};

}
