#include "GlobalDataBase.h"
#include "Geo_CellMetrics_Struct.h"

namespace PHSPACE
{
Geo_CellMetrics_Struct::Geo_CellMetrics_Struct()
{
    xcc = nullptr;
    ycc = nullptr;
    zcc = nullptr;
    vol = nullptr;

    xlen = nullptr;
    ylen = nullptr;
    zlen = nullptr;

    jacobian = nullptr;

    largestLocalGridLength = nullptr;
    smallestLocalGridLength = nullptr;
    subgridLength = nullptr;
}

Geo_CellMetrics_Struct::~Geo_CellMetrics_Struct()
{
    DeallocateAll();

    delete largestLocalGridLength; largestLocalGridLength = nullptr;
    delete smallestLocalGridLength; smallestLocalGridLength = nullptr;
    delete subgridLength; subgridLength = nullptr;
}

void Geo_CellMetrics_Struct::AllocateMetrics(Range &I, Range &J, Range &K)
{
    DeallocateAll();
    xcc = new RDouble3D(I, J, K, fortranArray);
    ycc = new RDouble3D(I, J, K, fortranArray);
    zcc = new RDouble3D(I, J, K, fortranArray);
    vol = new RDouble3D(I, J, K, fortranArray);

    xlen = new RDouble3D(I, J, K, fortranArray);
    ylen = new RDouble3D(I, J, K, fortranArray);
    zlen = new RDouble3D(I, J, K, fortranArray);

    int SolverStructOrder = PHSPACE::GlobalDataBase::GetIntParaFromDB("SolverStructOrder");
    if (SolverStructOrder > 2)
    {
        jacobian = new RDouble3D(I, J, K, fortranArray);
    }
}
}