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

    bool WantVisualField(Grid *grid);

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

}