#include "PBBoundary.h"
#include "Data_Param.h"

namespace PHSPACE
{
    PBBoundary::PBBoundary(const string &name):PBObject(name)
    {

    }

    PBBoundary::~PBBoundary()
    {

    }

    DEFINE_PBSOLVER_REGISTER_FUNCTION(PBBoundary, PBBoundary, BOUNDARY_EMPTY);
}
