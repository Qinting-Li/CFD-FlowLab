#include "Geo_FaceMetrics_Unstruct.h"

namespace PHSPACE
{
LIB_EXPORT Geo_FaceMetrics_Unstruct::Geo_FaceMetrics_Unstruct()
{
    xfc = nullptr;
    yfc = nullptr;
    zfc = nullptr;

    xfn = nullptr;
    yfn = nullptr;
    zfn = nullptr;

    area = nullptr;
}

LIB_EXPORT Geo_FaceMetrics_Unstruct::~Geo_FaceMetrics_Unstruct()
{
    delete [] xfc;    xfc = nullptr;
    delete [] yfc;    yfc = nullptr;
    delete [] zfc;    zfc = nullptr;

    delete [] xfn;    xfn = nullptr;
    delete [] yfn;    yfn = nullptr;
    delete [] zfn;    zfn = nullptr;

    delete [] area;    area = nullptr;
}

}