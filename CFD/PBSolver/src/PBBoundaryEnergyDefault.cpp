#include "PBBoundaryEnergyDefault.h"

namespace PHSPACE
{
    PBBoundaryEnergyDefault::PBBoundaryEnergyDefault(const string &name):PBBoundaryScalarDefault(name)
    {

    }

    PBBoundaryEnergyDefault::~PBBoundaryEnergyDefault()
    {

    }

    DEFINE_PBSOLVER_REGISTER_FUNCTION(PBBoundary, PBBoundaryEnergyDefault, ENERGY_DEFAULT);
}