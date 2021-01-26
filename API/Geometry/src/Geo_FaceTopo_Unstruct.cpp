#include "Geo_FaceTopo_Unstruct.h"

namespace PHSPACE
{
LIB_EXPORT Geo_FaceTopo_Unstruct::Geo_FaceTopo_Unstruct()
{
	left_cell_of_face  = 0;
	right_cell_of_face = 0;
	node_number_of_each_face = 0;
	face2node		   = 0;

    face2nodeSubscript = 0;
	face2nodeArray     = 0;    
}

LIB_EXPORT Geo_FaceTopo_Unstruct::~Geo_FaceTopo_Unstruct()
{
	delete [] left_cell_of_face;
	delete [] right_cell_of_face;
	delete [] node_number_of_each_face;
	delete [] face2node;

    if( face2nodeSubscript != 0 )
    {
        delete [] face2nodeSubscript;
    }

	if( face2nodeArray != 0 )
	{
		delete [] face2nodeArray;
	}
}

}