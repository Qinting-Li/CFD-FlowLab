#include "MBGridFactory.h"
#include "MBConstTable.h"
#include "Geo_UnstructBC.h"
#include "Geo_UnstructGrid.h"
using namespace std;

namespace PHSPACE
{
MBGridFactory::MBGridFactory()
{
    this->faceToCellRelation  = NULL;
    this->physicalType        = NULL;
    this->numberOfNodesOnFace = NULL;
    this->nodeIndexListOnFace = NULL;

    this->geometricType = NULL;
    this->faceIndexList = NULL;

    this->outerZoneList = NULL;
    this->outerFaceList = NULL;

    this->x = NULL;
    this->y = NULL;
    this->z = NULL;

    this->iblank = NULL;

    this->sourceRankList = NULL;
    this->sourceCellList = NULL;
    this->targetCellList = NULL;
    this->lagrangeFactor = NULL;
}

MBGridFactory::~MBGridFactory()
{
    delete [] faceToCellRelation ; faceToCellRelation  = NULL;
    delete [] physicalType       ; physicalType        = NULL;
    delete [] numberOfNodesOnFace; numberOfNodesOnFace = NULL;
    delete [] nodeIndexListOnFace; nodeIndexListOnFace = NULL;

    delete [] geometricType; geometricType = NULL;
    delete [] faceIndexList; faceIndexList = NULL;

    delete [] outerZoneList; outerZoneList = NULL;
    delete [] outerFaceList; outerFaceList = NULL;

    delete [] x; x = NULL;
    delete [] y; y = NULL;
    delete [] z; z = NULL;

    delete [] iblank; iblank = NULL;

    delete [] sourceRankList; sourceRankList = NULL;
    delete [] sourceCellList; sourceCellList = NULL;
    delete [] targetCellList; targetCellList = NULL;
    delete [] lagrangeFactor; lagrangeFactor = NULL;
}

void MBGridFactory::UpdateDonorsSpace(int numberOfLocalDonors)
{
    this->numberOfLocalDonors = numberOfLocalDonors;
    delete [] sourceRankList; sourceRankList = new int [numberOfLocalDonors];
    delete [] sourceCellList; sourceCellList = new int [numberOfLocalDonors];
    delete [] targetCellList; targetCellList = new int [numberOfLocalDonors];

    delete [] lagrangeFactor; lagrangeFactor = new RDouble [8 * numberOfLocalDonors];
}

void MBGridFactory::Register(UnstructGrid * grid)
{
    this->numberOfNodes = grid->GetNTotalNode();
    this->numberOfCells = grid->GetNTotalCell();
    this->numberOfFaces = grid->GetNTotalFace();

    this->numberOfOuterFaces    = grid->GetNIFace();
    this->numberOfPhysicalFaces = grid->GetNBoundFace() - numberOfOuterFaces;
    this->numberOfInnerFaces    = numberOfFaces - grid->GetNBoundFace();

    delete [] faceToCellRelation ; faceToCellRelation  = new int [2 * numberOfFaces];
    delete [] physicalType       ; physicalType        = new int [    numberOfFaces];
    delete [] numberOfNodesOnFace; numberOfNodesOnFace = new int [    numberOfFaces];
    delete [] nodeIndexListOnFace; nodeIndexListOnFace = new int [4 * numberOfFaces];

    int * leftCellOfFace  = grid->GetLeftCellOfFace();
    int * rightCellOfFace = grid->GetRightCellOfFace();
    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        int isl = 2 * iFace;
        faceToCellRelation[isl] = leftCellOfFace[iFace];

        int isr = isl + 1;
        if (rightCellOfFace[iFace] < numberOfCells)
        {
            faceToCellRelation[isr] = rightCellOfFace[iFace];
        }
        else
        {
            faceToCellRelation[isr] = MB_SPACE::GHOST_INDEX;
        }

        physicalType[iFace] = MB_SPACE::NEUTRAL_TYPE;
    }

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();

    int numberOfBcFaces = grid->GetNBoundFace();
    for (int iBcFace = 0; iBcFace < numberOfBcFaces; ++ iBcFace)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iBcFace]);

        int bcType = bcRegion->GetBCType();

        if (IsInternalFace(bcType))
        {
            physicalType[iBcFace] = MB_SPACE::INTERNAL_TYPE;
        }
        else if (IsExternalFace(bcType))
        {
            physicalType[iBcFace] = MB_SPACE::EXTERNAL_TYPE;
        }
        else
        {
            physicalType[iBcFace] = bcType;
        }
    }

    int * nodeListSize = grid->GetNodeNumberOfEachFace();
    memcpy(numberOfNodesOnFace, nodeListSize, numberOfFaces * sizeof(int));

    int * faceToNodeMapping = grid->GetFace2Node();
    int counter = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        int size  = numberOfNodesOnFace[iFace];
        int * dst = nodeIndexListOnFace + iFace * 4;
        int * src = faceToNodeMapping + counter;
        memcpy(dst, src, size * sizeof(int));
        counter += size;
    }

    delete [] iblank; iblank = new int [numberOfCells];

    delete [] geometricType; geometricType = new int [    numberOfCells];
    delete [] faceIndexList; faceIndexList = new int [6 * numberOfCells];

    int * nodeSizeOnCell = grid->GetNodeNumberOfEachCell();
    int ** cellToFaceMapping = grid->GetCell2Face();
    for (int iCell = 0; iCell < numberOfCells; ++ iCell)
    {
        int size = nodeSizeOnCell[iCell];
        if (size == 4)
        {
            geometricType[iCell] = MB_SPACE::TETRAHEDRON;
        }
        else if (size == 8)
        {
            geometricType[iCell] = MB_SPACE::HEXAHEDRON;
        }
        else if (size == 6)
        {
            geometricType[iCell] = MB_SPACE::THREE_PRISM;
        }
        else
        {
            geometricType[iCell] = MB_SPACE::PYRAMID;
        }

        memcpy(& faceIndexList[6 * iCell], cellToFaceMapping[iCell], size * sizeof(int));
    }

    delete [] outerZoneList; outerZoneList = new int [numberOfOuterFaces];
    delete [] outerFaceList; outerFaceList = new int [numberOfOuterFaces];

    InterfaceInfo * info = grid->GetInterfaceInfo();
    memcpy(outerZoneList, info->GetInterFace2ZoneID(), numberOfOuterFaces * sizeof(int));
    memcpy(outerFaceList, info->GetInterFace2InterFaceID(), numberOfOuterFaces * sizeof(int));

    delete [] x; x = new RDouble [numberOfNodes];
    delete [] y; y = new RDouble [numberOfNodes];
    delete [] z; z = new RDouble [numberOfNodes];

    memcpy(x, grid->GetX(), numberOfNodes * sizeof(RDouble));
    memcpy(y, grid->GetY(), numberOfNodes * sizeof(RDouble));
    memcpy(z, grid->GetZ(), numberOfNodes * sizeof(RDouble));

    return;
}

bool MBGridFactory::IsInternalFace(int type)
{
    return type == PHENGLEI::SOLID_SURFACE;
}

bool MBGridFactory::IsExternalFace(int type)
{
    return type == PHENGLEI::EXTERNAL_BC;
}

}