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
//! @file      Geo_NodeTopo_Unstruct.h
//! @brief     It defines the node topology of the unstructured grid,
//!            such as the node-face, node-cell, node-node connected relationship.
//! @author    Bell, He Xin.

#pragma once
#include "LIB_Macro.h"
namespace PHSPACE
{
//! @brief Geo_NodeTopo_Unstruct class defines the node-face, node-cell, node-node \n
//! connected relationship of unstructured grid.\n
//! 1: cells index list of each node.
class Geo_NodeTopo_Unstruct
{
private:
    //! Number of connected cells per node.
    int * cell_number_of_each_node;          

    //! Node to cell: connected cells index of per node.
    int * node2cell;

    //! 2D Node to cell: connected cells index of per node.
    int **node2CellArray;

public:
    LIB_EXPORT Geo_NodeTopo_Unstruct();
    LIB_EXPORT ~Geo_NodeTopo_Unstruct();

public:
    //! Set the number of connected cells per node.
    void SetCellNumberOfEachNode(int * nCPN);

    //! Set node to cell: connected cells index of per node.
    void SetNode2Cell(int * n2c);

    //! Get the number of connected cells per node.
    int * GetCellNumberOfEachNode() const;

    //! Get node to cell: connected cells index of per node.
    int * GetNode2Cell() const;

    //! Set 2-D node to cell: connected cells index of per node.
    void SetNode2CellArray(int ** node2Cell);

    //! Get 2-D node to cell: connected cells index of per node.
    int ** GetNode2CellArray() const;
};

#include "Geo_NodeTopo_Unstruct.hxx"





}