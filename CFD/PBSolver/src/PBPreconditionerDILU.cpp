#include "PBPreconditionerDILU.h"

namespace PHSPACE
{
PBPreconditionerDILU::PBPreconditionerDILU(const string &name):PBPreconditioner(name)
{

}

PBPreconditionerDILU::~PBPreconditionerDILU()
{

}

void PBPreconditionerDILU::Preconditioning(Region &region, vector<int> &zoneId, string &varStr, string &apStr, string &ajUpperStr,
                                           string &ajLowerStr, string &aiiStr, string &scStr, string &rsStr)
{

}

void PBPreconditionerDILU::PreconditioningProxy(int nFaces, int nCells, int *leftCellOfFace, int *rightCellOfFace,RDouble *ap,
                                                RDouble *ajLeft, RDouble *ajRight, RDouble *rs, string *phiName)
{
    for (int iCell = 0; iCell < nCells; ++ iCell)
    {
        rs[iCell] = ap[iCell];
    }

    for (int iFace = 0; iFace < nFaces; ++ iFace)
    {
        int le = leftCellOfFace[iFace];
        int re = rightCellOfFace[iFace];
        if (std::abs(rs[le]) > 1e-20)
        {
            rs[re] -= ajLeft[iFace] * ajRight[iFace] / rs[le];
        }
    }

    for (int iCell = 0; iCell < nCells; ++ iCell)
    {
        rs[iCell] = 1.0 / rs[iCell];
    }
}

void PBPreconditionerDILU::ForwardAndBackward(int nFaces, int nCells, int *leftCellOfFace, int *rightCellOfFace, RDouble *ap0,RDouble *ajLeft,
                                              RDouble *ajRight, RDouble *p, RDouble *rs, string *phiName)
{
    for (int iCell = 0; iCell < nCells; ++ iCell)
    {
        rs[iCell] = p[iCell];
    }

    for (int iFace = 0; iFace < nFaces; ++ iFace)
    {
        int le = leftCellOfFace[iFace];
        int re = rightCellOfFace[iFace];
        rs[re] += ajRight[iFace] * ap0[le] * rs[le];
    }

    for (int iCell = 0; iCell < nCells; ++ iCell)
    {
        rs[iCell] = rs[iCell] * ap0[iCell];
    }

    for (int iFace = nFaces - 1; iFace > -1; -- iFace)
    {
        int le = leftCellOfFace[iFace];
        int re = rightCellOfFace[iFace];
        rs[le] += ajLeft[iFace] * ap0[le] * rs[re];
    }

}

DEFINE_PBSOLVER_REGISTER_FUNCTION(PBPreconditioner, PBPreconditionerDILU, ITERSOLVPRECOND_DILU);
}