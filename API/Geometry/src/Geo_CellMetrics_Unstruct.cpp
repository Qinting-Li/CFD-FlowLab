#include "Geo_CellMetrics_Unstruct.h"

namespace PHSPACE
{
LIB_EXPORT Geo_CellMetrics_Unstruct::Geo_CellMetrics_Unstruct()
{
    xcc = nullptr;
    ycc = nullptr;
    zcc = nullptr;
    vol = nullptr;

    largestLocalGridLength = nullptr;
    subgridLength = nullptr;
}

LIB_EXPORT Geo_CellMetrics_Unstruct::~Geo_CellMetrics_Unstruct()
{
    delete [] xcc; xcc = nullptr;
    delete [] ycc; ycc = nullptr;
    delete [] zcc; zcc = nullptr;
    delete [] vol; vol = nullptr;

    delete [] largestLocalGridLength; largestLocalGridLength = nullptr;
    delete [] subgridLength; subgridLength = nullptr;
}

}