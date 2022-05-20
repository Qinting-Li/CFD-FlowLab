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
//! @file      Post_WriteTecplot.h
//! @brief     Write flow field into post visual file.
//! @author    Xu Gang.

#pragma once
#include "LIB_Macro.h"
#include "TypeDefine.h"

namespace PHSPACE
{
class ActionKey;
class Grid;
class DataContainer;
class Post_Visual;

class Post_WriteTecplot
{
public:
    LIB_EXPORT Post_WriteTecplot(int flowType = 0);
    LIB_EXPORT ~Post_WriteTecplot();

public:
    LIB_EXPORT void Run(bool visualAllLevel = false);

private:
    vector<DataContainer *> datalist;
    vector<DataContainer *> vtkdatalist;

    int flowType;

private:
    void ServerWrite(ActionKey *actkey);

    //! Store the visualization variables.
    void StoreVisualizationVar(ActionKey *actkey, Post_Visual *postVisulaization);

    void PostWrite(ActionKey *actkey);

    //!
    void BoundaryVisualization(Grid *grid_in, ActionKey *actkey, RDouble **qn, int nvarplot);

    //!
    void BoundaryVTKVisualization(Grid *grid_in, ActionKey *actkey, RDouble **qn, int nvarplot);

    //!
    void BoundaryVisualization(Grid *grid_in, ActionKey *actkey, RDouble4D **qn, int nvarplot);

    //!
    void BoundaryVTKVisualization(Grid *grid_in, ActionKey *actkey, RDouble4D **qn, int nvarplot);

    //!
    void SliceVisualization(Grid *grid_in, ActionKey *actkey, RDouble **qn, int nvarplot);

    //!
    void FieldVisualization(Grid *grid_in, ActionKey *actkey, RDouble **qn, int nvarplot);

    //!
    void FieldVisualizationForVTK(Grid *grid_in, ActionKey *actkey, RDouble **qn, int nvarplot);

    //!
    void FieldVisualization(Grid *grid_in, ActionKey *actkey, RDouble4D **qn, int nvarplot);

    //!
    void SaveDataForTecio(Grid *grid_in, ActionKey *actkey, RDouble **qn, int nvarplot);

    //!
    void SaveDataForTecio(Grid *grid_in, ActionKey *actkey, RDouble4D **qn, int nl);

    void Clear();
};

class Post_WriteTecplotByOriginalGrid
{
private:
    int nOriginalGrids;
    Grid **OriginalGrid;

public:
    Post_WriteTecplotByOriginalGrid();
    ~Post_WriteTecplotByOriginalGrid();

public:
    void Run();

private:
    void ServerWrite(ActionKey *actkey);
    void ComputerNodeDataOnPartitionGrid(int iZone, ActionKey *actkey);
    void CollectNodeData(ActionKey *actkey);
    void CollectStructedData(ActionKey *actkey, Post_Visual *postVisualization);

    void StoreVisualizationVar(ActionKey *actkey);
    void BoundaryVTKVisualization(Grid *grid_in, DataContainer *cdata, RDouble4D **qn, int nvarplot);
    void FieldVisualization(Grid *grid_in, DataContainer *cdata, RDouble4D **qn, int nvarplot);
    void BoundaryVisualization(Grid *grid_in, DataContainer *cdata, RDouble4D **qn, int nvarplot);
    void SaveDataForTecio(Grid *grid_in, DataContainer *cdata, RDouble4D **qn, int nl);

    void PostWrite(ActionKey *actkey);

    void UpdataNodeData(ActionKey *actkey, Grid *grid, RDouble4D **qn, int nVisualVariables);

    int GetOriginalGridProc(int iZone);

public:
    void Initialize();
    void SetOriginalGrid(Grid **gridsIn) { this->OriginalGrid = gridsIn; };
    void SetOriginalGridProc(int nOriginalGridsIn, int *originalGridProcOut);

private:
    int *originalGridProc;
    int *originalGridIndexOfEachGrid;
    vector <Post_Visual * > originalGridData;
    vector<DataContainer *> datalist;
    vector<DataContainer *> vtkdatalist;

};

Post_WriteTecplotByOriginalGrid *GetWriteTecplotByOriginalGrid();
bool WantVisualField(Grid *grid);


}