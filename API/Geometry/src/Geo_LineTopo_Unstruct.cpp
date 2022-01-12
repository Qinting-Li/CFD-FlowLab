#include "Geo_LineTopo_Unstruct.h"

namespace PHSPACE
{
    LIB_EXPORT Geo_LineTopo_Unstruct::Geo_LineTopo_Unstruct()
    {
        nLine = 0;
        nCellsInLines = NULL;
        CellOfLine = NULL;
        FaceOfLine = NULL;
        LineOfCell = NULL;
    }

    LIB_EXPORT Geo_LineTopo_Unstruct::~Geo_LineTopo_Unstruct()
    {
        if (nCellsInLines != NULL)
        {
            DelPointer(nCellsInLines);
        }

        if (CellOfLine != NULL)
        {
            DelPointer2(CellOfLine);
        }
        
        if (FaceOfLine != NULL)
        {
            DelPointer2(FaceOfLine);
        }

        if (LineOfCell != NULL)
        {
            DelPointer(LineOfCell);
        }
    }

}