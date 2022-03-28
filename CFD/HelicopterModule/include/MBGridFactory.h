#pragma once
#include "Precision.h"
using namespace std;

namespace PHSPACE
{
class UnstructGrid;

class MBGridFactory
{
private:
    int numberOfNodes, numberOfCells, numberOfFaces;
    int numberOfPhysicalFaces, numberOfOuterFaces, numberOfInnerFaces;
    int * faceToCellRelation, * physicalType, * numberOfNodesOnFace, * nodeIndexListOnFace;
    int * geometricType, * faceIndexList;
    int * outerZoneList, * outerFaceList;
    int * iblank;
    RDouble * x, * y, * z;

    int numberOfLocalDonors;
    int * sourceRankList, * sourceCellList, * targetCellList;
    RDouble * lagrangeFactor;
public:
    MBGridFactory();
    ~MBGridFactory();
public:
    void UpdateDonorsSpace(int numberOfLocalDonors);
    void Register(UnstructGrid * grid);
private:
    bool IsInternalFace(int type);
    bool IsExternalFace(int type);
};

}