#pragma once
using namespace std;

namespace PHSPACE
{
const int MBDIM = 3;

namespace MB_SPACE
{
const int GHOST_INDEX = -1;

const int NEUTRAL_TYPE  = 2002;
const int INTERNAL_TYPE = 2005;
const int EXTERNAL_TYPE = 2006;

const int TETRAHEDRON = 0;
const int HEXAHEDRON  = 1;
const int THREE_PRISM = 2;
const int PYRAMID     = 3;
}

}