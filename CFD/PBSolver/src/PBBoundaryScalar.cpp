#include "PBBoundaryScalar.h"

namespace PHSPACE
{
PBBoundaryScalar::PBBoundaryScalar(const string &name):PBBoundary(name)
{

}

PBBoundaryScalar::~PBBoundaryScalar()
{

}
void PBBoundaryScalar::Init(string phiName)
{

}

void PBBoundaryScalar::Discrt(UnstructGrid *grid, string phiName, vector<int> *faceIndex, Data_Param *bcData, RDouble *diffCoe)
{

}

void PBBoundaryScalar::Update(UnstructGrid *grid, string phiName, vector<int> *faceIndex, Data_Param *bcData)
{
    int nTotalCells = grid->GetNTotalCell();
    int *leftCellOfFace  = grid->GetLeftCellOfFace();
    RDouble *phi = reinterpret_cast<RDouble *>(grid->GetDataPtr(phiName));
    RDouble *phib = phi + nTotalCells;
    for (size_t iFace = 0; iFace < faceIndex->size(); ++ iFace)
    {
        int iFacelocal = (*faceIndex)[iFace];
        int le = leftCellOfFace[iFacelocal];
        phib[iFacelocal] = phi[le];
    }
}

DEFINE_PBSOLVER_REGISTER_FUNCTION(PBBoundary, PBBoundaryScalar, SCALAR_EMPTY);
} 