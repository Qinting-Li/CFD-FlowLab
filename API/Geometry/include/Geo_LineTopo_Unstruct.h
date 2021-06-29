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
//! @file      Geo_LineTopo_Unstruct.h
//! @brief     It defines the face topology of the unstructured grid,
//!              such as the face-cell and face-node connected relationship.
//! @author    WanYunBo.
#pragma once
#include "LIB_Macro.h"
#include "Pointer.h"
namespace PHSPACE
{
    //! @brief Geo_LineTopo_Unstruct class defines the line for line implicit LUSGS.\n
    class Geo_LineTopo_Unstruct
    {
    private:
        //! the number of lines.
        int nLine;

        //! the number of cells in all lines.
        int* nCellsInLines;

        //! the cell No. of all lines.
        int** CellOfLine;

        //! the face No. of all lines.
        int** FaceOfLine;

        //! the line No. of all cells.
        int* LineOfCell;

    public:
        LIB_EXPORT Geo_LineTopo_Unstruct();
        LIB_EXPORT ~Geo_LineTopo_Unstruct();

    public:
        //! set the number of lines.
        void SetLI_nLine(int nLine);

        //! set the number of cells in all lines.
        void SetLI_nCellsInLine(int* nCellsInLines);

        //! set the cell No. of all lines.
        void SetLI_CellOfLine(int** CellOfLine);

        //! set the face No. of all lines.
        void SetLI_FaceOfLine(int** FaceOfLine);

        //! set the line No. of all cells.
        void SetLI_LineOfCell(int* LineOfCell);

        //! get the number of lines.
        int GetLI_nLine() const;

        //! get the number of cells in all lines.
        int* GetLI_nCellsInLine() const;

        //! get the cell No. of all lines.
        int** GetLI_CellOfLine() const;

        //! get the face No. of all lines.
        int** GetLI_FaceOfLine() const;

        //! get the line No. of all cells.
        int* GetLI_LineOfCell() const;
    };

#include "Geo_LineTopo_Unstruct.hxx"
}
