#include "Geo_NodeTopo_Unstruct.h"

namespace PHSPACE
{
LIB_EXPORT Geo_NodeTopo_Unstruct::Geo_NodeTopo_Unstruct()
{
    cell_number_of_each_node = 0;
    node2cell = 0;
    node2CellArray = 0;
}

LIB_EXPORT Geo_NodeTopo_Unstruct::~Geo_NodeTopo_Unstruct()
{
    if (cell_number_of_each_node != 0)
    {
        delete [] cell_number_of_each_node;
        cell_number_of_each_node = 0;
    }

    if (node2cell != 0)
    {
        delete [] node2cell;
        node2cell = 0;
    }

    if (node2CellArray != 0)
    {
        delete [] node2CellArray;
        node2CellArray = 0;
    }
}

}