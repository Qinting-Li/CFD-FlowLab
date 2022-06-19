#include "PBBoundary.h"
#include "Data_Param.h"

namespace PHSPACE
{
PBBoundary::PBBoundary(const string &name):PBObject(name)
{

}

PBBoundary::~PBBoundary()
{

}

void PBBoundary::Init(string phiName)
{

}

void PBBoundary::Discrt(UnstructGrid *grid, string phiName, vector<int> *faceIndex, Data_Param *bcData, RDouble *diffCoe)
{

}

void PBBoundary::Update(UnstructGrid *grid, string phiName, vector<int> *faceIndex, Data_Param *bcData)
{

}

void PBBoundary::WallModification(UnstructGrid *grid, string phiName, vector<int> *faceIndex, Data_Param *bcData)
{

}

DEFINE_PBSOLVER_REGISTER_FUNCTION(PBBoundary, PBBoundary, BOUNDARY_EMPTY);
}
