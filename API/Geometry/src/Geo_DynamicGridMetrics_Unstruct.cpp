#include "Geo_DynamicGridMetrics_Unstruct.h"

namespace PHSPACE
{
LIB_EXPORT Geo_DynamicGridMetrics_Unstruct::Geo_DynamicGridMetrics_Unstruct()
{
	voln = 0;
	vgn  = 0;
	xtn  = 0;
	ytn  = 0;
	ztn  = 0;
}

LIB_EXPORT Geo_DynamicGridMetrics_Unstruct::~Geo_DynamicGridMetrics_Unstruct()
{
	if ( voln ) delete [] voln; voln = 0;
	if ( vgn  ) delete [] vgn ; vgn  = 0;
	if ( xtn  ) delete [] xtn ; xtn  = 0;
	if ( ytn  ) delete [] ytn ; ytn  = 0;
	if ( ztn  ) delete [] ztn ; ztn  = 0;
}

}