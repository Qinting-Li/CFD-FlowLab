#include "PBDiffusion.h"

namespace PHSPACE
{
    PBDiffusion::PBDiffusion(const string &name):PBObject(name)
    {

    }

    PBDiffusion::~PBDiffusion()
    {

    }

    DEFINE_PBSOLVER_REGISTER_FUNCTION(PBDiffusion, PBDiffusion, DIFFUSION_EMPTY);
}