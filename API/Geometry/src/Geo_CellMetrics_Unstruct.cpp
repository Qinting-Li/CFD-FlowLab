#include "Geo_CellMetrics_Unstruct.h"
namespace PHSPACE
{
LIB_EXPORT Geo_CellMetrics_Unstruct::Geo_CellMetrics_Unstruct()
{
	xcc = 0;
	ycc = 0;
	zcc = 0;
	vol = 0;

	largestLocalGridLength		= 0;
	subgridLength				= 0;
}

LIB_EXPORT Geo_CellMetrics_Unstruct::~Geo_CellMetrics_Unstruct()
{
	delete [] xcc; xcc = 0;
	delete [] ycc; ycc = 0;
	delete [] zcc; zcc = 0;
	delete [] vol; vol = 0;

	delete [] largestLocalGridLength; largestLocalGridLength = 0;
	delete [] subgridLength; subgridLength = 0;
}

}