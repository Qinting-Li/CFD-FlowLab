#include "Pointer.h"
#include "Precision.h"
#include "DSMCCollision.h"
using namespace std;

namespace PHSPACE
{

Collision::Collision(int nTotalCell)
{
    VRC.resize(3);
    VCM.resize(3);

    VRC[0] = VRC[1] = VRC[2] = 0.0;
    VCM[0] = VCM[1] = VCM[2] = 0.0;
    VRR = 0.0;VR = 0.0; CVR =0.0;
    CCG = NewPointer2<RDouble>(2,nTotalCell);

}

Collision::~Collision()
{
    DelPointer2 (CCG) ; CCG = NULL;

}

}



