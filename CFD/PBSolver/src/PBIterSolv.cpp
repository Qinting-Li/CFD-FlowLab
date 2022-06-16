#include "PBIterSolv.h"
#include "PBPreconditioner.h"
#include "Task_ServerUpdateInterface.h"

namespace PHSPACE
{
    PBIterSolv::PBIterSolv(const string &name):PBObject(name)
    {
        preconditioner_ = PBPreconditionerFactory::CreatePBPreconditioner("DILU");
    }

    PBIterSolv::~PBIterSolv()
    {
        if (preconditioner_)
        {
            delete preconditioner_;
        }
    }

    void PBIterSolv::Matrix(int nFaces, int nCells, int *leftCellOfFace, int *rightCellOfFace, RDouble *ap, RDouble *ajLeft, RDouble *ajRight, RDouble *p, RDouble *rs, std::string *phiName)
    {
        for (int iCell = 0; iCell < nCells; ++ iCell)
        {
            rs[iCell] = ap[iCell] * p[iCell];
        }

        for (int iFace = 0; iFace < nFaces; ++ iFace)
        {
            int le = leftCellOfFace[iFace];
            int re = rightCellOfFace[iFace];
            rs[le] -= ajLeft[iFace] * p[re];
            rs[re] -= ajRight[iFace] * p[le];
        }

    }

    RDouble PBIterSolv::Dot(RDouble *b, RDouble *c, int rows)
    {
        RDouble ans = 0.0;
        for (int i = 0; i < rows; ++ i)
        {
            ans += b[i] * c[i];
        }
        return(ans);
    }

    void PBIterSolv::Rotg(RDouble *X, RDouble *Y, RDouble *c, RDouble *s)
    {
        RDouble temp;
        if (*Y == 0.0)
        {
            *c = 1.0;
            *s = 0.0;
        }
        else if (abs(*X) > abs(*Y))
        {
            temp = *Y / *X;
            *c = 1.0 / sqrt(1.0 + temp * temp);
            *s = temp * (*c);
        }
        else if (abs(*X) < abs(*Y))
        {
            temp = (*X) / (*Y);
            *s = 1.0 / sqrt(1.0 + temp * temp);
            *c = temp * (*s);
        }
        RDouble tmp = (*X) * (*c) + (*Y) * (*s);
        (*Y) = -(*s) * (*X) + (*c) * (*Y);
        (*X) = tmp;
    }

    void PBIterSolv::Rot(RDouble  *X, RDouble *Y, RDouble *c, RDouble *s)
    {
        RDouble tmp = (*X) * (*c) + (*Y) * (*s);
        (*Y) = -(*s) * (*X) + (*c) * (*Y);
        (*X) = tmp;
    }

    DEFINE_PBSOLVER_REGISTER_FUNCTION(PBIterSolv, PBIterSolv, ITERSOLV_EMPTY);
}