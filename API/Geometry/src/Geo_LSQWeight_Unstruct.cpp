#include "Geo_LSQWeight_Unstruct.h"

namespace PHSPACE
{
LIB_EXPORT Geo_LSQWeight_Unstruct::Geo_LSQWeight_Unstruct()
{
	iwt = 0;

	ixx = 0;
	iyy = 0;
	izz = 0;
	ixy = 0;
	ixz = 0;
	iyz = 0;

	fMark = 0;
}

LIB_EXPORT Geo_LSQWeight_Unstruct::~Geo_LSQWeight_Unstruct()
{
	delete [] iwt;

	delete [] ixx;
	delete [] iyy;
	delete [] izz;
	delete [] ixy;
	delete [] ixz;
	delete [] iyz;

	delete [] fMark;
}

}