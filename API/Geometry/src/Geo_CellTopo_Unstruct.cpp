#include "Geo_CellTopo_Unstruct.h"
#include "Pointer.h"

namespace PHSPACE
{
LIB_EXPORT Geo_CellTopo_Unstruct::Geo_CellTopo_Unstruct()
{
	face_number_of_each_cell	= 0;
	cell2face					= 0;

	number_of_neighbor_cell		= 0;
	cell2cell					= 0;

	node_number_of_each_cell    = 0;
	cell2node					= 0;
	cell2nodeArray              = 0;
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