#include "PBDiffusion.h"

namespace PHSPACE
{
PBDiffusion::PBDiffusion(const string &name):PBObject(name)
{

}

PBDiffusion::~PBDiffusion()
{

}

void PBDiffusion::Discrt(UnstructGrid *grid, string phiName, RDouble *diffCoe, RDouble *dphidx, RDouble *dphidy, RDouble *dphidz,
                    RDouble *ajLeft, RDouble *ajRight, RDouble *scPhi, RDouble *dvdx, RDouble *dvdy, RDouble *dvdz,
                    RDouble *scV, RDouble *dwdx, RDouble *dwdy, RDouble *dwdz, RDouble *scW)
{

}

DEFINE_PBSOLVER_REGISTER_FUNCTION(PBDiffusion, PBDiffusion, DIFFUSION_EMPTY);
}