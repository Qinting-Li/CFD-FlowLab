#include "Geo_DynamicGridMetrics_Unstruct.h"

namespace PHSPACE
{
LIB_EXPORT Geo_DynamicGridMetrics_Unstruct::Geo_DynamicGridMetrics_Unstruct()
{
    voln = NULL;
    vgn  = NULL;
    xfv  = NULL;
    yfv  = NULL;
    zfv  = NULL;
}

LIB_EXPORT Geo_DynamicGridMetrics_Unstruct::~Geo_DynamicGridMetrics_Unstruct()
{
    if (voln) delete [] voln; voln = NULL;
    if (vgn) delete [] vgn ; vgn  = NULL;
    if (xfv) delete [] xfv ; xfv  = NULL;
    if (yfv) delete [] yfv ; yfv  = NULL;
    if (zfv) delete [] zfv ; zfv  = NULL;
}

}