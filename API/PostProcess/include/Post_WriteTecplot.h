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
//! @brief     Write flow field into tecplot file.
//! @author    Xu Gang.

#pragma once
#include "Post_WriteVisualFile.h"

namespace PHSPACE
{

//! @brief Post_WriteTecplot class realize tecpolt visual file output function.\n
class Post_WriteTecplot : public Post_WriteVisualFile
{
public:
    Post_WriteTecplot();
    ~Post_WriteTecplot();

private:
    void Initialize();

    //! Collect visual data for tecplot visual file.
    void StoreVisualizationData();
    void StoreVisualizationData(int zoneIndex, DataContainer *cdata);
    void StoreStrBoundaryVisualData(int zoneIndex, DataContainer *cdata);
    void StoreStrFieldVisualData(int zoneIndex, DataContainer *cdata);
    void StoreUnsBoundaryVisualData(int zoneIndex, DataContainer *cdata);
    void StoreUnsFieldVisualData(int zoneIndex, DataContainer *cdata);

    //! Dump data into visual file.
    void WriteFile();

    //! Clear data.
    void ClearFieldData();

private:
    vector<DataContainer *> visualDataList;

};

}
