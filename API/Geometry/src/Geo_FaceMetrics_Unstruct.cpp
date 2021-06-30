#include "Geo_FaceMetrics_Unstruct.h"

namespace PHSPACE
{
LIB_EXPORT Geo_FaceMetrics_Unstruct::Geo_FaceMetrics_Unstruct()
{
	xfc = 0;
	yfc = 0;
	zfc = 0;

	xfn = 0;
	yfn = 0;
	zfn = 0;

	area = 0;
}

LIB_EXPORT Geo_FaceMetrics_Unstruct::~Geo_FaceMetrics_Unstruct()
{
	delete [] xfc;
	delete [] yfc;
	delete [] zfc;

	delete [] xfn;
	delete [] yfn;
	delete [] zfn;

	delete [] area;
}

}