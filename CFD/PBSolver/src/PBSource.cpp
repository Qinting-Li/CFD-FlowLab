#include "PBSource.h"

namespace PHSPACE
{
    PBSource::PBSource(const string &name)
        :PBObject(name)
    {

    }

    PBSource::~PBSource()
    {

    }

    DEFINE_PBSOLVER_REGISTER_FUNCTION(PBSource, PBSource, SOURCE_EMPTY);
}