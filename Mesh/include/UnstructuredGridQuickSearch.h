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
//! @file      UnstructuredGridQuickSearch.h
//! @brief     Quick search of host cells in unstructured grid.
//! @author    Ma Rong, Chang Xinghua.
#pragma once
#include "Precision.h"
#include <vector>
using namespace std;

namespace PHSPACE
{
class UnstructGrid;
class OversetConfig;

class MRADTSearch
{
public:
    MRADTSearch(UnstructGrid *gridIn);
    MRADTSearch(UnstructGrid *gridIn, OversetConfig *oversetConfig);
    ~MRADTSearch();
private:
    UnstructGrid *grid;    //! Grid pointer of this block.

    int numberOfDimension;
    int numberOfBoxDirection;

    int numberOfSearchCells;

    vector< int > searchCellsIndex;

    RDouble  *globalBox;
    RDouble **cellBox;
    
    int *leftSonIndexList;
    int *rightSonIndexList;
    int *reorderCellIndexList;

    int rootCellIndex;

    OversetConfig *oversetConfig;

public:
    UnstructGrid * GetGrid() const { return grid; }
public:
    void BuildADTStructure();

    void ComputeGlobalAndCellBox();
    void ComputeGlobalAndCellBox2D();
    void ComputeGlobalAndCellBox3D();
    void ComputeGlobalBox2D();
    void ComputeGlobalBox3D();

    void ADTSort(int iLay, int iStart, int iEnd, RDouble *adtSortData);

    void QuickSort(int iStart, int iEnd,    RDouble *adtSortData, int half);
    void SwapData (int iFirst, int iSecond, RDouble *adtSortData);

    void SearchInAdt(RDouble x, RDouble y, RDouble z, int &donorIndex);
    void SearchPoint2D(int iLay, int rootCellIndexNow, RDouble *coor, int &donorIndex);
    void SearchPoint3D(int iLay, int rootCellIndexNow, RDouble *coor, int &donorIndex);
    void InCellCheck(int rootCellIndexNow, RDouble *coor, int &donorIndex);

    int  In4PointsCell(RDouble x1, RDouble y1, RDouble z1,
                       RDouble x2, RDouble y2, RDouble z2,
                       RDouble x3, RDouble y3, RDouble z3,
                       RDouble x4, RDouble y4, RDouble z4,
                       RDouble xt, RDouble yt, RDouble zt, 
                       RDouble eps);

public:
    int GetNumberOfSearchCells() { return numberOfSearchCells; }
    RDouble ComputeVolume4P(RDouble x1, RDouble y1, RDouble z1, RDouble x2, RDouble y2, RDouble z2, RDouble x3, RDouble y3, RDouble z3, RDouble x4, RDouble y4, RDouble z4);
};

}