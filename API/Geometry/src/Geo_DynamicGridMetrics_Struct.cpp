#include "Geo_DynamicGridMetrics_Struct.h"

namespace PHSPACE
{

Geo_DynamicGridMetrics_Struct::Geo_DynamicGridMetrics_Struct()
{
    voln = 0;
    xtn = 0;
    ytn = 0;
    ztn = 0;
    vgn = 0;
}

Geo_DynamicGridMetrics_Struct::~Geo_DynamicGridMetrics_Struct()
{
    DeallocateAll();
}

}