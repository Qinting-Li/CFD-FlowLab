#include "PBBoundaryTurb.h"

using namespace std;

namespace PHSPACE
{
PBBoundaryTurb::PBBoundaryTurb(const string &name):PBBoundary(name)
{

}

PBBoundaryTurb::~PBBoundaryTurb()
{

}

void PBBoundaryTurb::Discrt(UnstructGrid *grid, string phiName, vector<int> *faceIndex, Data_Param *bcData, RDouble *diffCoe)
{

}

void PBBoundaryTurb::Update(UnstructGrid *grid, string phiName, vector<int> *faceIndex, Data_Param *bcData)
{

}

DEFINE_PBSOLVER_REGISTER_FUNCTION(PBBoundary, PBBoundaryTurb, TURB_EMPTY);
} // PB_SPACE