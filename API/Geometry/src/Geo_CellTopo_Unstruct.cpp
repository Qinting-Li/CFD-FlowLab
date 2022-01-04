#include "Geo_CellTopo_Unstruct.h"
#include "Pointer.h"

namespace PHSPACE
{
LIB_EXPORT Geo_CellTopo_Unstruct::Geo_CellTopo_Unstruct()
{
    face_number_of_each_cell = nullptr;
    cell2face = nullptr;

    number_of_neighbor_cell = nullptr;
    cell2cell = nullptr;

    node_number_of_each_cell = nullptr;
    cell2node = nullptr;
    cell2nodeArray = nullptr;
}

LIB_EXPORT Geo_CellTopo_Unstruct::~Geo_CellTopo_Unstruct()
{
    DelPointer2(cell2face);
    DelPointer2(cell2nodeArray);

    delete [] face_number_of_each_cell;
    delete [] cell2node;

    delete [] number_of_neighbor_cell;
    delete [] cell2cell;

    if ( node_number_of_each_cell ) delete [] node_number_of_each_cell;
}

}