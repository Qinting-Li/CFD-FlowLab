#include "PBSource.h"

namespace PHSPACE
{
PBSource::PBSource(const string &name):PBObject(name)
{

}

PBSource::~PBSource()
{

}

void PBSource::InitMemory(UnstructGrid *grid, string phiName)
{

}

void PBSource::Discrt(UnstructGrid *grid, string phiName)
{

}

void PBSource::UpdateBodyForce(UnstructGrid *grid, string phiName)
{

}

DEFINE_PBSOLVER_REGISTER_FUNCTION(PBSource, PBSource, SOURCE_EMPTY);
}