#include "PBConvection.h"

namespace PHSPACE
{
PBConvection::PBConvection(const string &name):PBObject(name)
{

}

PBConvection::~PBConvection()
{

}

void PBConvection::Discrt(UnstructGrid *grid, string phiName, string apName, string scName, RDouble *v,RDouble *apV,
                          RDouble *scV, RDouble *w, RDouble *apW, RDouble *scW)
{

}

DEFINE_PBSOLVER_REGISTER_FUNCTION(PBConvection, PBConvection, CONVECTION_EMPTY);
} 