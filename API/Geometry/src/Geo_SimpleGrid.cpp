#include "Geo_SimpleGrid.h"
#include "Math_BasisFunction.h"
#include "PHMpi.h"
using namespace std;

#pragma warning (disable:913)
namespace PHSPACE
{

SimpleGrid::SimpleGrid()
{
    x = NULL;
    y = NULL;
    z = NULL;

    pmin = new RDouble [3];
    pmax = new RDouble [3];

    pmin[0] = 0.0;
    pmin[1] = 0.0;
    pmin[2] = 0.0;

    pmax[0] = 1.0;
    pmax[1] = 1.0;
    pmax[2] = 1.0;

    this->has_cal_minmax = false;

    this->minEdgeLength = -1.0;
    this->maxEdgeLength = -1.0;

    nTotalNode = 0;
}

SimpleGrid::~SimpleGrid()
{
    FreeMemory();
}

SimpleGrid::SimpleGrid(const SimpleGrid &rhs)
{
    Copy(rhs);
    this->minEdgeLength = -1.0;
    this->maxEdgeLength = -1.0;
}

void SimpleGrid::FreeMemory()
{
    if (x != NULL) { delete [] x; x=NULL; }
    if (y != NULL) { delete [] y; y=NULL; }
    if (z != NULL) { delete [] z; z=NULL; }

    if (pmin != NULL) { delete [] pmin; pmin=NULL; }
    if (pmax != NULL) { delete [] pmax; pmax=NULL; }
}

void SimpleGrid::Copy(const SimpleGrid &rhs)
{
    nTotalNode = rhs.GetNTotalNode();

    if (nTotalNode <= 0)
    {
        PHMPI::FreeBlockData();
        cout << " FATAL Error for nTotalNode <= 0 in SimpleGrid Copy!\n" << endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    x = new RDouble [nTotalNode];
    y = new RDouble [nTotalNode];
    z = new RDouble [nTotalNode];

    std::copy(rhs.x, rhs.x+nTotalNode, x);
    std::copy(rhs.y, rhs.y+nTotalNode, y);
    std::copy(rhs.z, rhs.z+nTotalNode, z);

    pmin = new RDouble [3];
    pmax = new RDouble [3];

    std::copy(rhs.pmin, rhs.pmin+3, pmin);
    std::copy(rhs.pmax, rhs.pmax+3, pmax);
}

SimpleGrid & SimpleGrid::operator = (const SimpleGrid &rhs)
{
    if (this == &rhs) return *this;

    FreeMemory();

    Copy(rhs);

    return *this;
}

RDouble * SimpleGrid::GetMinBox()
{
    if (!has_cal_minmax)
    {
        ComputeMinMaxBox();
    }
    return pmin;
};

RDouble * SimpleGrid::GetMaxBox()
{
    if (!has_cal_minmax)
    {
        ComputeMinMaxBox();
    }
    return pmax;
};

LIB_EXPORT void SimpleGrid::ComputeMinMaxBox()
{
    this->has_cal_minmax = true;

    int nTotalNode = this->GetNTotalNode();

    RDouble *x = this->GetX();
    RDouble *y = this->GetY();
    RDouble *z = this->GetZ();

    RDouble xmin, ymin, zmin, xmax, ymax, zmax;
    xmin = ymin = zmin =   LARGE;
    xmax = ymax = zmax = - LARGE;

    for (int inode = 0; inode < nTotalNode; ++ inode)
    {
        xmin = MIN(xmin, x[inode]);
        ymin = MIN(ymin, y[inode]);
        zmin = MIN(zmin, z[inode]);

        xmax = MAX(xmax, x[inode]);
        ymax = MAX(ymax, y[inode]);
        zmax = MAX(zmax, z[inode]);
    }

    RDouble *pmin = this->GetMinBox();
    RDouble *pmax = this->GetMaxBox();

    pmin[0] = xmin;
    pmin[1] = ymin;
    pmin[2] = zmin;

    pmax[0] = xmax;
    pmax[1] = ymax;
    pmax[2] = zmax;
}

LIB_EXPORT RDouble SimpleGrid::CalMinEdgeLength()
{
    return 0.0;
}

LIB_EXPORT RDouble SimpleGrid::CalMaxEdgeLength()
{
    return 0.0;
}

}