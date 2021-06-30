#include "UnstructuredGridQuickSearch.h"
#include "UnstructuredOversetConfig.h"
#include "Geo_UnstructGrid.h"
#include "Glb_Dimension.h"
#include "Pointer.h"

namespace PHSPACE
{

MRADTSearch::MRADTSearch(UnstructGrid *gridIn)
{
    this->grid = gridIn;

    oversetConfig = 0;

    int numberOfAllCells = grid->GetNTotalCell();

    numberOfSearchCells = 0;

    for (int iCell = 0; iCell < numberOfAllCells; ++ iCell)
    {
        searchCellsIndex.push_back(numberOfSearchCells);
        numberOfSearchCells ++;
    }

    this->numberOfDimension    = PHSPACE::GetDim();
    this->numberOfBoxDirection = 2 * numberOfDimension;

    this->globalBox = new RDouble [numberOfBoxDirection];
    this->cellBox   = NewPointer2< RDouble >(numberOfSearchCells, numberOfBoxDirection);
    
    this->leftSonIndexList  = new int [numberOfSearchCells];
    this->rightSonIndexList = new int [numberOfSearchCells];

    this->reorderCellIndexList = new int [numberOfSearchCells];

    SetField(leftSonIndexList , - 1, numberOfSearchCells);
    SetField(rightSonIndexList, - 1, numberOfSearchCells);

    for (int iCell = 0; iCell < numberOfSearchCells; ++ iCell)
    {
        reorderCellIndexList[iCell] = iCell;
    }
    rootCellIndex = 0;
}

MRADTSearch::MRADTSearch(UnstructGrid *gridIn, OversetConfig *oversetConfigIn)
{
    this->grid = gridIn;

    this->oversetConfig = oversetConfigIn;
    
    int numberOfAllCells = grid->GetNTotalCell();

    vector< int > &keySearchOfNodes = oversetConfig->GetKeySearchOfNodes();
    vector< int > &keySearchOfCells = oversetConfig->GetKeySearchOfCells();

    numberOfSearchCells = 0;

    for (int iCell = 0; iCell < numberOfAllCells; ++ iCell)
    {
        if (keySearchOfCells[iCell] == 1)
        {
            searchCellsIndex.push_back(iCell);
            numberOfSearchCells ++;
        }
    }

    //! numberOfSearchCells = 0，it shows that this grid is not in the overlapping area，so that there are no cells involved in the search.
    if (numberOfSearchCells == 0)
    {
        return;
    }

    this->numberOfDimension    = PHSPACE::GetDim();
    this->numberOfBoxDirection = 2 * numberOfDimension;

    this->globalBox = new RDouble [numberOfBoxDirection];
    this->cellBox   = NewPointer2< RDouble >(numberOfSearchCells, numberOfBoxDirection);
    
    this->leftSonIndexList  = new int [numberOfSearchCells];
    this->rightSonIndexList = new int [numberOfSearchCells];

    this->reorderCellIndexList = new int [numberOfSearchCells];

    SetField(leftSonIndexList , - 1, numberOfSearchCells);
    SetField(rightSonIndexList, - 1, numberOfSearchCells);

    for (int iCell = 0; iCell < numberOfSearchCells; ++ iCell)
    {
        reorderCellIndexList[iCell] = iCell;
    }

    rootCellIndex = 0;
}

MRADTSearch::~MRADTSearch()
{
    if (numberOfSearchCells == 0)
    {
        return;
    }

    delete [] globalBox;
    DelPointer2(cellBox);

    delete [] leftSonIndexList ;
    delete [] rightSonIndexList;

    delete [] reorderCellIndexList;
}

void MRADTSearch::BuildADTStructure()
{
    if (numberOfSearchCells == 0)
    {
        return;
    }

    ComputeGlobalAndCellBox();

    if (numberOfSearchCells == 1)
    {
        rootCellIndex = 0;
        return;
    }

    int iLay   = 0;
    int iStart = 0;
    int iEnd   = numberOfSearchCells - 1;

    RDouble *adtSortData = new RDouble [numberOfSearchCells];

    ADTSort(iLay, iStart, iEnd, adtSortData);

    rootCellIndex = reorderCellIndexList[iEnd / 2];

    delete [] adtSortData;
}

void MRADTSearch::ADTSort(int iLay, int iStart, int iEnd, RDouble *adtSortData)
{
    int iDirection = iLay % numberOfBoxDirection;

    int half = (iStart + iEnd) / 2;

    for (int i = iStart; i <= iEnd; ++ i)
    {
        adtSortData[i] = cellBox[reorderCellIndexList[i]][iDirection];
    }

    //! Quick sort.
    QuickSort(iStart, iEnd, adtSortData, half);

    //! Analyze the Left side.
    if (half - iStart == 0)    //! Nothing on the Left, left subtree is empty.
    {
        leftSonIndexList[reorderCellIndexList[half]] = - 1;
    }
    else if (half - iStart == 1)    //! Only one on the left, left subtree is clear.
    {
        leftSonIndexList[reorderCellIndexList[half]] = reorderCellIndexList[iStart];
    }
    else    //! Other cases, further analysis is needed on the left side, the root node of the left subtree is the left half position.
    {
        int newStart = iStart;
        int newEnd   = half - 1;
        ADTSort(iLay + 1, newStart, newEnd, adtSortData);
        leftSonIndexList[reorderCellIndexList[half]] = reorderCellIndexList[(newStart + newEnd) / 2];
    }

    //! Analyze the right side (Due to the rounding of half, the right subtree cannot be empty).
    if (half + 1 == iEnd)    //! Only one on the right, right subtree is clear.
    {
        rightSonIndexList[reorderCellIndexList[half]] = reorderCellIndexList[iEnd];
    }
    else    //! Other cases, further analysis is needed on the right side, the root node of the right subtree is the left half position.
    {
        int newStart = half + 1;
        int newEnd   = iEnd;
        ADTSort(iLay + 1, newStart, newEnd, adtSortData);
        rightSonIndexList[reorderCellIndexList[half]] = reorderCellIndexList[(newStart + newEnd) / 2];
    }
}

void MRADTSearch::QuickSort(int iStart, int iEnd, RDouble *adtSortData, int half)
{
    int last;

    if (iEnd - iStart < 1)
    {
        return;
    }

    SwapData(iStart, (iStart + iEnd) / 2, adtSortData);

    last = iStart;
    for (int i = iStart + 1; i <= iEnd; ++ i)
    {
        if (adtSortData[i] < adtSortData[iStart])
        {
            SwapData(++ last, i, adtSortData);
        }
    }
    SwapData(iStart, last, adtSortData);

    if (last == half)
    {
        return;
    }

    if (last > half)
    {
        QuickSort(iStart, last - 1, adtSortData, half);
    }

    if (last < half)
    {
        QuickSort(last + 1, iEnd, adtSortData, half);
    }
}

void MRADTSearch::SwapData(int iFirst, int iSecond, RDouble *adtSortData)
{
    RDouble tempAdtSortData;
    int     tempReorderCellIndex;

    tempAdtSortData      = adtSortData         [iFirst];
    tempReorderCellIndex = reorderCellIndexList[iFirst];

    adtSortData         [iFirst] = adtSortData         [iSecond];
    reorderCellIndexList[iFirst] = reorderCellIndexList[iSecond];

    adtSortData         [iSecond] = tempAdtSortData;
    reorderCellIndexList[iSecond] = tempReorderCellIndex;
}

void MRADTSearch::ComputeGlobalAndCellBox2D()
{
    RDouble *nodeCoordinateX = grid->GetX();
    RDouble *nodeCoordinateY = grid->GetY();
    RDouble *nodeCoordinateZ = grid->GetZ();

    int **cellNodeIndexContainer  = grid->GetCell2NodeArray();
    int * cellNodeNumberContainer = grid->GetNodeNumberOfEachCell();
    
    RDouble GREAT = 1.0E+20;

    ComputeGlobalBox2D();

    RDouble eps = GetToleranceForOversetBox();

    //! Calculate the cellBox.
    for (int iCell = 0; iCell < numberOfSearchCells; ++ iCell)
    {
        int cellIndex = searchCellsIndex[iCell];

        cellBox[iCell][0] = GREAT;cellBox[iCell][2] = - GREAT;
        cellBox[iCell][1] = GREAT;cellBox[iCell][3] = - GREAT;

        for (int iNodeInCell = 0; iNodeInCell < cellNodeNumberContainer[cellIndex]; ++ iNodeInCell)
        {
            int nodeIndex = cellNodeIndexContainer[cellIndex][iNodeInCell];
            cellBox[iCell][0] = min(cellBox[iCell][0], nodeCoordinateX[nodeIndex]); //! 本单元最小x.
            cellBox[iCell][1] = min(cellBox[iCell][1], nodeCoordinateY[nodeIndex]); //! 本单元最小y.
            cellBox[iCell][2] = max(cellBox[iCell][2], nodeCoordinateX[nodeIndex]); //! 本单元最大x.
            cellBox[iCell][3] = max(cellBox[iCell][3], nodeCoordinateY[nodeIndex]); //! 本单元最大y.
        }

        //! Keep out residuals.
        cellBox[iCell][0] -= eps;
        cellBox[iCell][1] -= eps;
        cellBox[iCell][2] += eps;
        cellBox[iCell][3] += eps;
    }
}

void MRADTSearch::ComputeGlobalBox2D()
{
    //! Calculate the globalBox. If there is oversetConfig, there is no need to solve it again.
    if (oversetConfig != 0)
    {
        globalBox[0] = oversetConfig->pMin[0];
        globalBox[1] = oversetConfig->pMin[1];
        globalBox[2] = oversetConfig->pMax[0];
        globalBox[3] = oversetConfig->pMax[1];

        return;
    }

    int numberOfNodes = grid->GetNTotalNode();
    int numberOfCells = grid->GetNTotalCell();

    int **cellNodeIndexContainer  = grid->GetCell2NodeArray();
    int * cellNodeNumberContainer = grid->GetNodeNumberOfEachCell();

    RDouble *x = grid->GetX();
    RDouble *y = grid->GetY();
    RDouble *z = grid->GetZ();

    RDouble xmin, ymin, zmin, xmax, ymax, zmax;
    xmin = ymin = zmin =   LARGE;
    xmax = ymax = zmax = - LARGE;

    for (int iCell = 0; iCell < numberOfCells; ++ iCell)
    {
        int cellNodeNumber = cellNodeNumberContainer[iCell];

        for (int iNode = 0; iNode < cellNodeNumber; ++ iNode)
        {
            int nodeIndex = cellNodeIndexContainer[iCell][iNode];

            xmin = PHSPACE::MIN(xmin, x[nodeIndex]);
            ymin = PHSPACE::MIN(ymin, y[nodeIndex]);
            zmin = PHSPACE::MIN(zmin, z[nodeIndex]);

            xmax = PHSPACE::MAX(xmax, x[nodeIndex]);
            ymax = PHSPACE::MAX(ymax, y[nodeIndex]);
            zmax = PHSPACE::MAX(zmax, z[nodeIndex]);
        }
    }

    //! minMaxBox keep out residuals.
    RDouble eps = GetToleranceForOversetBox();

    globalBox[0] = xmin - eps;
    globalBox[1] = ymin - eps;

    globalBox[2] = xmax + eps;
    globalBox[3] = ymax + eps;
}

void MRADTSearch::ComputeGlobalAndCellBox3D()
{
    RDouble *nodeCoordinateX = grid->GetX();
    RDouble *nodeCoordinateY = grid->GetY();
    RDouble *nodeCoordinateZ = grid->GetZ();

    int **cellNodeIndexContainer  = grid->GetCell2NodeArray();
    int * cellNodeNumberContainer = grid->GetNodeNumberOfEachCell();

    RDouble GREAT = 1.0E+20;

    RDouble eps = GetToleranceForOversetBox();

    //! Calculate globalBox.
    ComputeGlobalBox3D();

    //! Calculate cellBox.
    for (int iCell = 0; iCell < numberOfSearchCells; ++ iCell)
    {
        int cellIndex = searchCellsIndex[iCell];

        cellBox[iCell][0] = GREAT;cellBox[iCell][3] = - GREAT;
        cellBox[iCell][1] = GREAT;cellBox[iCell][4] = - GREAT;
        cellBox[iCell][2] = GREAT;cellBox[iCell][5] = - GREAT;

        for (int iNodeInCell = 0; iNodeInCell < cellNodeNumberContainer[cellIndex]; ++ iNodeInCell)
        {
            int nodeIndex = cellNodeIndexContainer[cellIndex][iNodeInCell];
            cellBox[iCell][0] = min(cellBox[iCell][0], nodeCoordinateX[nodeIndex]); //! 本单元最小x.
            cellBox[iCell][1] = min(cellBox[iCell][1], nodeCoordinateY[nodeIndex]); //! 本单元最小y.
            cellBox[iCell][2] = min(cellBox[iCell][2], nodeCoordinateZ[nodeIndex]); //! 本单元最小z.
            cellBox[iCell][3] = max(cellBox[iCell][3], nodeCoordinateX[nodeIndex]); //! 本单元最大x.
            cellBox[iCell][4] = max(cellBox[iCell][4], nodeCoordinateY[nodeIndex]); //! 本单元最大y.
            cellBox[iCell][5] = max(cellBox[iCell][5], nodeCoordinateZ[nodeIndex]); //! 本单元最大z.
        }

        //! Keep out residuals.
        cellBox[iCell][0] -= eps;
        cellBox[iCell][1] -= eps;
        cellBox[iCell][2] -= eps;
        cellBox[iCell][3] += eps;
        cellBox[iCell][4] += eps;
        cellBox[iCell][5] += eps;
    }
}

void MRADTSearch::ComputeGlobalBox3D()
{
    //! Calculate globalBox.
    if (oversetConfig != 0)
    {
        globalBox[0] = oversetConfig->pMin[0];
        globalBox[1] = oversetConfig->pMin[1];
        globalBox[2] = oversetConfig->pMin[2];
        globalBox[3] = oversetConfig->pMax[0];
        globalBox[4] = oversetConfig->pMax[1];
        globalBox[5] = oversetConfig->pMax[2];

        return;
    }

    int numberOfNodes = grid->GetNTotalNode();
    int numberOfCells = grid->GetNTotalCell();

    RDouble *x = grid->GetX();
    RDouble *y = grid->GetY();
    RDouble *z = grid->GetZ();

    int **cellNodeIndexContainer  = grid->GetCell2NodeArray();
    int * cellNodeNumberContainer = grid->GetNodeNumberOfEachCell();

    RDouble xmin, ymin, zmin, xmax, ymax, zmax;
    xmin = ymin = zmin =   LARGE;
    xmax = ymax = zmax = - LARGE;

    for (int iCell = 0; iCell < numberOfCells; ++ iCell)
    {
        int cellNodeNumber = cellNodeNumberContainer[iCell];

        for (int iNode = 0; iNode < cellNodeNumber; ++ iNode)
        {
            int nodeIndex = cellNodeIndexContainer[iCell][iNode];

            xmin = PHSPACE::MIN(xmin, x[nodeIndex]);
            ymin = PHSPACE::MIN(ymin, y[nodeIndex]);
            zmin = PHSPACE::MIN(zmin, z[nodeIndex]);

            xmax = PHSPACE::MAX(xmax, x[nodeIndex]);
            ymax = PHSPACE::MAX(ymax, y[nodeIndex]);
            zmax = PHSPACE::MAX(zmax, z[nodeIndex]);
        }
    }

    //! minMaxBox keep out residuals.
    RDouble eps = GetToleranceForOversetBox();

    globalBox[0] = xmin - eps;
    globalBox[1] = ymin - eps;
    globalBox[2] = zmin - eps;

    globalBox[3] = xmax + eps;
    globalBox[4] = ymax + eps;
    globalBox[5] = zmax + eps;
}

void MRADTSearch::ComputeGlobalAndCellBox()
{
    if (numberOfDimension == 2)
    {
        ComputeGlobalAndCellBox2D();
    }
    else if (numberOfDimension == 3)
    {
        ComputeGlobalAndCellBox3D();
    }
    else
    {
        cout << "暂不支持异次元重叠网格……\n" << endl;
        exit(0);
    }
}

void MRADTSearch::SearchInAdt(RDouble x, RDouble y, RDouble z, int &donorIndex)
{
    if (numberOfSearchCells == 0)
    {
        return;
    }

    RDouble *coor = new RDouble [numberOfDimension];

    int iLay = 0;
    int rootCellIndexNow = this->rootCellIndex;

    if (numberOfDimension == 2)
    {
        coor[0] = x; coor[1] = y;

        if (globalBox[0] <= coor[0] && globalBox[2] >= coor[0] &&
            globalBox[1] <= coor[1] && globalBox[3] >= coor[1])    //! Firstly, this point falls in the bounding box formed by the maximum and minimum values of the contribution grid.
        {
            SearchPoint2D(iLay, rootCellIndexNow, coor, donorIndex);
        }
    }
    else if (numberOfDimension == 3)
    {
        coor[0] = x; coor[1] = y; coor[2] = z;

        if (globalBox[0] <= coor[0] && globalBox[3] >= coor[0] &&
            globalBox[1] <= coor[1] && globalBox[4] >= coor[1] &&
            globalBox[2] <= coor[2] && globalBox[5] >= coor[2] )    //! Firstly, this point falls in the bounding box formed by the maximum and minimum values of the contribution grid.
        {
            SearchPoint3D(iLay, rootCellIndexNow, coor, donorIndex);
        }
    }
    else
    {
        cout << "暂不支持异次元重叠网格……\n" << endl;
        exit(0);
    }

    delete [] coor;
}

void MRADTSearch::InCellCheck(int rootCellIndexNow, RDouble *coor, int &donorIndex)
{
    int cellIndex = searchCellsIndex[rootCellIndexNow];

    int key;

    RDouble eps = GetToleranceForOversetSearch();

    RDouble xn, yn, zn;
    RDouble xf, yf, zf;
    RDouble xc, yc, zc;
    RDouble x1, y1, z1;
    RDouble x2, y2, z2;
    RDouble volume;

    int *leftCellIndexContainer = grid->GetLeftCellOfFace();

    int **cellFaceIndexContainer  = grid->GetCell2Face();
    int  *cellFaceNumberContainer = grid->GetFaceNumberOfEachCell();

    int  *faceNodeNumberContainer = grid->GetNodeNumberOfEachFace();
    int **faceNodeIndexVector     = grid->GetFace2NodeArray();

    RDouble *faceNormalXContainer = grid->GetFaceNormalX();
    RDouble *faceNormalYContainer = grid->GetFaceNormalY();
    RDouble *faceNormalZContainer = grid->GetFaceNormalZ();

    RDouble *faceCenterXContainer = grid->GetFaceCenterX();
    RDouble *faceCenterYContainer = grid->GetFaceCenterY();
    RDouble *faceCenterZContainer = grid->GetFaceCenterZ();

    RDouble *cellCenterXContainer = grid->GetCellCenterX();
    RDouble *cellCenterYContainer = grid->GetCellCenterY();
    RDouble *cellCenterZContainer = grid->GetCellCenterZ();

    RDouble *xContainer = grid->GetX();
    RDouble *yContainer = grid->GetY();
    RDouble *zContainer = grid->GetZ();

    int numberOfFacesInCell = cellFaceNumberContainer[cellIndex];
    
    if (numberOfDimension == 2)
    {
        key = 1;
        for (int iFace = 0; iFace < numberOfFacesInCell; ++ iFace)
        {
            int faceIndex = cellFaceIndexContainer[cellIndex][iFace];

            xn = faceNormalXContainer[faceIndex];
            yn = faceNormalYContainer[faceIndex];
            zn = faceNormalZContainer[faceIndex];

            xf = faceCenterXContainer[faceIndex];
            yf = faceCenterYContainer[faceIndex];
            zf = faceCenterZContainer[faceIndex];

            volume = (coor[0] - xf) * xn + (coor[1] - yf) * yn + (coor[2] - zf) * zn;
            if (leftCellIndexContainer[faceIndex] == cellIndex)
            {
                volume = - volume;
            }

            //! eps is small negative value, which is more likely to be judged as in the cell to avoid the situation that cannot be found;
            if (volume < - eps)
            {
                key = 0;
                break;    //! The point is located on the outside of a certain face, and certainly not in the interior of the cell, so the subsequent detection is skipped.
            }
        }
    }
    else
    {
        key = 0;

        xc = cellCenterXContainer[cellIndex];
        yc = cellCenterYContainer[cellIndex];
        zc = cellCenterZContainer[cellIndex];

        for (int iFace = 0; iFace < numberOfFacesInCell; ++ iFace)
        {
            int faceIndex = cellFaceIndexContainer[cellIndex][iFace];
            int faceNodeNumber = faceNodeNumberContainer[faceIndex];

            xf = faceCenterXContainer[faceIndex];
            yf = faceCenterYContainer[faceIndex];
            zf = faceCenterZContainer[faceIndex];

            for (int iNode = 0; iNode < faceNodeNumber; ++ iNode)
            {
                int nodeIndex1 = faceNodeIndexVector[faceIndex][iNode];
                int nodeIndex2 = faceNodeIndexVector[faceIndex][(iNode + 1) % faceNodeNumber];

                x1 = xContainer[nodeIndex1];
                y1 = yContainer[nodeIndex1];
                z1 = zContainer[nodeIndex1];
                x2 = xContainer[nodeIndex2];
                y2 = yContainer[nodeIndex2];
                z2 = zContainer[nodeIndex2];

                key = In4PointsCell(x1, y1, z1, x2, y2, z2, xf, yf, zf, xc, yc, zc, coor[0], coor[1], coor[2], eps);
                if (key == 1)    //! If the return value is 1, it indicates that it has been found and jump out of the loop.
                {
                    break;
                }
            }
            if (key == 1)    //! If the return value is 1, it indicates that it has been found and jump out of the loop.
            {
                break;
            }
        }
    }

    //! If the contribution cells are found, the distance will be calculated according to the area weight.
    if (key == 1)
    {
        donorIndex = cellIndex;
    }
}

int MRADTSearch::In4PointsCell(RDouble x1, RDouble y1, RDouble z1, RDouble x2, RDouble y2, RDouble z2,RDouble x3, RDouble y3, RDouble z3, RDouble x4, RDouble y4, RDouble z4, RDouble xt, RDouble yt, RDouble zt, RDouble eps)
{
    RDouble vol  = ABS(this->ComputeVolume4P(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4));
    
    RDouble vol1 = ABS(this->ComputeVolume4P(x1, y1, z1, x2, y2, z2, x3, y3, z3, xt, yt, zt)) / vol;
    RDouble vol2 = ABS(this->ComputeVolume4P(x2, y2, z2, x3, y3, z3, x4, y4, z4, xt, yt, zt)) / vol;
    RDouble vol3 = ABS(this->ComputeVolume4P(x3, y3, z3, x4, y4, z4, x1, y1, z1, xt, yt, zt)) / vol;
    RDouble vol4 = ABS(this->ComputeVolume4P(x4, y4, z4, x1, y1, z1, x2, y2, z2, xt, yt, zt)) / vol;

    if (ABS(vol1 + vol2 + vol3 + vol4 - 1.0) > eps)    //! eps is small positive value, which is more likely to be judged as in the cell to avoid the situation that cannot be found.
    {
        return 0;    //! Sum is larger than solving alone, indicating that it is in the outside, return 0.
    }
    else
    {
        return 1;    //! Sum is equal to solving alone, indicating that it is in the inside, return 1.
    }
}

void MRADTSearch::SearchPoint2D(int iLay, int rootCellIndexNow, RDouble *coor, int &donorIndex)
{
    int rootNew, iDirection;

    if (rootCellIndexNow == -1)
    {
        return;
    }

    iDirection = iLay % numberOfBoxDirection;

    if (cellBox[rootCellIndexNow][0] <= coor[0] && cellBox[rootCellIndexNow][2] >= coor[0] &&
        cellBox[rootCellIndexNow][1] <= coor[1] && cellBox[rootCellIndexNow][3] >= coor[1])    //! Check whether the point is in the root node bounding box.
    {
        //! Check if it is really in the cell.
        InCellCheck(rootCellIndexNow, coor, donorIndex);
        if (donorIndex != - 1)
        {
            return;
        }

        //! This point is in the bounding box of the root node, but not really in the cell, can only traverse the left and right subtree of the root node.
        //! But the process is fast, because there's only a simple coordinate comparison.
        rootNew = leftSonIndexList[rootCellIndexNow];
        SearchPoint2D(iLay + 1, rootNew, coor, donorIndex);
        if (donorIndex != - 1)
        {
            return;
        }
        rootNew = rightSonIndexList[rootCellIndexNow];
        SearchPoint2D(iLay + 1, rootNew, coor, donorIndex);
    }
    else
    {
        if (iDirection == 0 && cellBox[rootCellIndexNow][0] > coor[0] ||
            iDirection == 1 && cellBox[rootCellIndexNow][1] > coor[1])    //! The point is on the left/top/front side of the bounding box, and found it in the left subtree.
        {
            rootNew = leftSonIndexList[rootCellIndexNow];
            SearchPoint2D(iLay + 1, rootNew, coor, donorIndex);
        }
        else if (iDirection == 2 && cellBox[rootCellIndexNow][2] < coor[0] ||
                 iDirection == 3 && cellBox[rootCellIndexNow][3] < coor[1])    //! The point is on the right/bottom/back side of the bounding box, and found it in the right subtree.
        {
            rootNew = rightSonIndexList[rootCellIndexNow];
            SearchPoint2D(iLay + 1, rootNew, coor, donorIndex);
        }
        else    //! The point is outside the boundary of the iDirection of the bounding box, there is no way, can noly traverse the left and right subtree of the root node.
                //! But the process is fast, because there's only a simple coordinate comparison.
        {
            rootNew = leftSonIndexList[rootCellIndexNow];
            SearchPoint2D(iLay + 1, rootNew, coor, donorIndex);
            if (donorIndex != - 1)
            {
                return;
            }
            rootNew = rightSonIndexList[rootCellIndexNow];
            SearchPoint2D(iLay + 1, rootNew, coor, donorIndex);
        }
    }
}

void MRADTSearch::SearchPoint3D(int iLay, int rootCellIndexNow, RDouble *coor, int &donorIndex)
{
    int rootNew, iDirection;

    if (rootCellIndexNow == -1)
    {
        return;
    }

    iDirection = iLay % numberOfBoxDirection;

    RDouble x = coor[0];
    RDouble y = coor[1];
    RDouble z = coor[2];

    //RDouble dist = DISTANCE((x - 5.07062376), (y + 0.619199344), (z - 3.29163274));

    //if ((dist < 10e-4) && (searchCellsIndex[rootCellIndexNow] == 6317))
    //{
    // cout << "node coordinate = " << x << "  "  << y << "  "  << z << "  " << endl;

    // cout << "cellBox min x   = " << cellBox[rootCellIndexNow][0] << endl;
    // cout << "cellBox min y   = " << cellBox[rootCellIndexNow][1] << endl;
    // cout << "cellBox min z   = " << cellBox[rootCellIndexNow][2] << endl;
    //
    // cout << "cellBox max x   = " << cellBox[rootCellIndexNow][3] << endl;
    // cout << "cellBox max y   = " << cellBox[rootCellIndexNow][4] << endl;
    // cout << "cellBox max z   = " << cellBox[rootCellIndexNow][5] << endl;
    //}

    if (cellBox[rootCellIndexNow][0] <= coor[0] && cellBox[rootCellIndexNow][3] >= coor[0] &&
        cellBox[rootCellIndexNow][1] <= coor[1] && cellBox[rootCellIndexNow][4] >= coor[1] &&
        cellBox[rootCellIndexNow][2] <= coor[2] && cellBox[rootCellIndexNow][5] >= coor[2])    //! Check whether the point is in the root node bounding box.
    {
        //! Check if it is really in the cell.
        InCellCheck(rootCellIndexNow, coor, donorIndex);

        if (donorIndex != - 1)
        {
            return;
        }

        //! This point is in the bounding box of the root node, but not really in the cell, can only traverse the left and right subtree of the root node.
        //! But the process is fast, because there's only a simple coordinate comparison.
        rootNew = leftSonIndexList[rootCellIndexNow];
        SearchPoint3D(iLay + 1, rootNew, coor, donorIndex);
        if (donorIndex != - 1)
        {
            return;
        }
        rootNew = rightSonIndexList[rootCellIndexNow];
        SearchPoint3D(iLay + 1, rootNew, coor, donorIndex);
    }
    else
    {
        if (iDirection == 0 && cellBox[rootCellIndexNow][0] > coor[0] ||
            iDirection == 1 && cellBox[rootCellIndexNow][1] > coor[1] ||
            iDirection == 2 && cellBox[rootCellIndexNow][2] > coor[2])    //! The point is on the left/top/front side of the bounding box, and found it in the left subtree.
        {
            rootNew = leftSonIndexList[rootCellIndexNow];
            SearchPoint3D(iLay + 1, rootNew, coor, donorIndex);
        }
        else if (iDirection == 3 && cellBox[rootCellIndexNow][3] < coor[0] ||
                 iDirection == 4 && cellBox[rootCellIndexNow][4] < coor[1] ||
                 iDirection == 5 && cellBox[rootCellIndexNow][5] < coor[2])    //! The point is on the right/bottom/back side of the bounding box, and found it in the right subtree.
        {
            rootNew = rightSonIndexList[rootCellIndexNow];
            SearchPoint3D(iLay + 1, rootNew, coor, donorIndex);
        }
        else    //! The point is inside the boundary of the iDirection of the bounding box, there is no way, can noly traverse the left and right subtree of the root node.
                //! But the process is fast, because there's only a simple coordinate comparison.
        {
            rootNew = leftSonIndexList[rootCellIndexNow];
            SearchPoint3D(iLay + 1, rootNew, coor, donorIndex);
            if (donorIndex != - 1)
            {
                return;
            }
            rootNew = rightSonIndexList[rootCellIndexNow];
            SearchPoint3D(iLay + 1, rootNew, coor, donorIndex);
        }
    }
}

RDouble MRADTSearch::ComputeVolume4P(RDouble x1, RDouble y1, RDouble z1, RDouble x2, RDouble y2, RDouble z2, RDouble x3, RDouble y3, RDouble z3, RDouble x4, RDouble y4, RDouble z4)
{
    RDouble xx1, yy1, zz1, xx2, yy2, zz2, xx3, yy3, zz3, vol;

    xx1 = x2 - x1;
    yy1 = y2 - y1;
    zz1 = z2 - z1;
    xx2 = x3 - x1;
    yy2 = y3 - y1;
    zz2 = z3 - z1;
    xx3 = x4 - x1;
    yy3 = y4 - y1;
    zz3 = z4 - z1;

    vol = xx1 * yy2 * zz3 + xx2 * yy3 * zz1 + xx3 * yy1 * zz2 - xx3 * yy2 * zz1 - xx2 * yy1 * zz3 - xx1 * yy3 * zz2;
    vol /= 6.0;

    return vol;
}

}