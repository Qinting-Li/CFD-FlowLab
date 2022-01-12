#include "Geo_FaceTopo_Unstruct.h"

namespace PHSPACE
{
LIB_EXPORT Geo_FaceTopo_Unstruct::Geo_FaceTopo_Unstruct()
{
    left_cell_of_face  = nullptr;
    right_cell_of_face = nullptr;
    node_number_of_each_face = nullptr;
    face2node = nullptr;

    face2nodeSubscript = nullptr;
    face2nodeArray = nullptr;
}

LIB_EXPORT Geo_FaceTopo_Unstruct::~Geo_FaceTopo_Unstruct()
{
    delete [] left_cell_of_face;    left_cell_of_face  = nullptr;
    delete [] right_cell_of_face;    right_cell_of_face = nullptr;
    delete [] node_number_of_each_face;    node_number_of_each_face = nullptr;
    delete [] face2node;    face2node = nullptr;

    if ( face2nodeSubscript != 0 )
    {
        delete [] face2nodeSubscript;    face2nodeSubscript = nullptr;
    }

    if ( face2nodeArray != 0 )
    {
        delete [] face2nodeArray;    face2nodeArray = nullptr;
    }
}

}