#include "PBBoundaryEnergy.h"

using namespace std;

namespace PHSPACE
{
PBBoundaryEnergy::PBBoundaryEnergy(const string &name):PBBoundaryScalar(name)
{

}

PBBoundaryEnergy::~PBBoundaryEnergy()
{

}

void PBBoundaryEnergy::Init(string phiName)
{

}

void PBBoundaryEnergy::Discrt(UnstructGrid *grid, string phiName, vector<int> *faceIndex, Data_Param *bcData, RDouble *diffCoe)
{

}

void PBBoundaryEnergy::Update(UnstructGrid *grid, string phiName, vector<int> *faceIndex, Data_Param *bcData)
{

}

DEFINE_PBSOLVER_REGISTER_FUNCTION(PBBoundary, PBBoundaryEnergy, ENERGY_EMPTY);
}