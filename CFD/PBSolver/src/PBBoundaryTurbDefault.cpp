#include "PBBoundaryTurbDefault.h"

using namespace std;

namespace PHSPACE
{
PBBoundaryTurbDefault::PBBoundaryTurbDefault(const string &name):PBBoundaryScalarDefault(name)
{

}

PBBoundaryTurbDefault::~PBBoundaryTurbDefault()
{

}

DEFINE_PBSOLVER_REGISTER_FUNCTION(PBBoundary, PBBoundaryTurbDefault, TURB_DEFAULT);
} //! PB_SPACE