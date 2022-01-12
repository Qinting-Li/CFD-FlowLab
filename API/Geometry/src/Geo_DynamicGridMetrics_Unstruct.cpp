#include "Geo_DynamicGridMetrics_Unstruct.h"

namespace PHSPACE
{
LIB_EXPORT Geo_DynamicGridMetrics_Unstruct::Geo_DynamicGridMetrics_Unstruct()
{
    voln = nullptr;
    vgn  = nullptr;
    xfv  = nullptr;
    yfv  = nullptr;
    zfv  = nullptr;
}

LIB_EXPORT Geo_DynamicGridMetrics_Unstruct::~Geo_DynamicGridMetrics_Unstruct()
{
    if (voln) delete [] voln; voln = nullptr;
    if (vgn ) delete [] vgn ; vgn  = nullptr;
    if (xfv ) delete [] xfv ; xfv  = nullptr;
    if (yfv ) delete [] yfv ; yfv  = nullptr;
    if (zfv ) delete [] zfv ; zfv  = nullptr;
}

}