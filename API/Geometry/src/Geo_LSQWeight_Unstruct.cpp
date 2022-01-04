#include "Geo_LSQWeight_Unstruct.h"

namespace PHSPACE
{
LIB_EXPORT Geo_LSQWeight_Unstruct::Geo_LSQWeight_Unstruct()
{
    iwt = nullptr;

    ixx = nullptr;
    iyy = nullptr;
    izz = nullptr;
    ixy = nullptr;
    ixz = nullptr;
    iyz = nullptr;

    fMark = nullptr;
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