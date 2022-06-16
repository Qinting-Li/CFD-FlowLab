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
// FILE:       XXX.cpp (XXX.h)                                            +
// AUTHOR(S):  Bell, Zhang Yang                                           +
// PURPOSE:    simulation object                                          +
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#pragma once
#include "TypeDefine.h"

namespace PHSPACE
{
class ZoneGridSolver;

class ComputationalElementSolver
{
public:
    ComputationalElementSolver ();
    ~ComputationalElementSolver();
protected:
    PHVectorInt1D leftElementIndexContainer;
    PHVectorInt1D rightElementIndexContainer;
    PHVectorInt2D faceNodeIndexesContainer;

    PHVector1D< set< int > > nodeToElement;
    PHVector1D< set< int > > elementToNode;
    PHVector1D< set< int > > elementToElement;

    int numberOfFaces;
    int numberOfNodes;
    int numberOfElements;
public:
    PHVectorInt1D & GetLeftElementIndexContainer () { return leftElementIndexContainer;  }
    PHVectorInt1D & GetRightElementIndexContainer() { return rightElementIndexContainer; }
    PHVectorInt2D & GetFaceNodeIndexesContainer  () { return faceNodeIndexesContainer;   }
public:
    int  GetNumberOfFaces() { return numberOfFaces; }
    void SetNumberOfFaces(int numberOfFaces) { this->numberOfFaces = numberOfFaces; }

    int  GetNumberOfNodes() { return numberOfNodes; }
    void SetNumberOfNodes(int numberOfNodes) { this->numberOfNodes = numberOfNodes; }

    int  GetNumberOfElements() { return numberOfElements; }
    void SetNumberOfElements(int numberOfElements) { this->numberOfElements = numberOfElements; }
public:
    PHVector1D< set< int > > & GetNodeToElement   () { return nodeToElement;    }
    PHVector1D< set< int > > & GetElementToNode   () { return elementToNode;    }
    PHVector1D< set< int > > & GetElementToElement() { return elementToElement; }
};

class AdaptivePropertySolver
{
public:
    AdaptivePropertySolver();
    ~AdaptivePropertySolver();
protected:
    ComputationalElementSolver *computationalElementSolver;
    ZoneGridSolver *zoneGridSolver;
public:
    ZoneGridSolver * GetZoneGridSolver() { return zoneGridSolver; }
    void SetZoneGridSolver(ZoneGridSolver *zoneGridSolver) { this->zoneGridSolver = zoneGridSolver; }
    ComputationalElementSolver * GetComputationalElementSolver() { return computationalElementSolver; }
public:
    void ReComputeAdaptionProperty();
protected:
    void Initialize();
    void InitializeNewStart();
    void InitializeReStart();
    void PrepareDataStructure();
    void ConstructNodeToElementStructure();
    void ConstructElementToElementStructure();
    void ReComputeAdaptionPropertyWarp();
    bool ReComputeAdaptionPropertyEveryIteration();
protected:
    void MergeElement(int computationalElementIndex, bool &iterationFlag);
    void SplitElement(int computationalElementIndex, bool &iterationFlag);

    void CheckParentElementForMerge  (int computationalElementIndex, bool &iterationFlag);
    void CheckChildElementForMerge   (int computationalElementIndex, bool &iterationFlag);
    void CheckNeighborElementForMerge(int computationalElementIndex, bool &iterationFlag);

    void SetComputationalElementModifiedStatus(int computationalElementIndex, int status);
    int  GetComputationalElementModifiedStatus(int computationalElementIndex);
protected:
    void SetCurrentElementsModifiedStatusByOldElementsModifiedStatus();
    void SetAllElementsModifiedStatusToNoChange();
    void SetAllElementsComputationalStatusToHiddenState();
    void SetCurrentElementsComputationalStatusToOnIfOldElementsModifiedStatusIsNoChange();

    void SetElementModifiedStatus(int elementIndex, int elementModifiedStatus);
};

}
