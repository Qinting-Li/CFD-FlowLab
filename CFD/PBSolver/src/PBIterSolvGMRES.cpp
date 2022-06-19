#include "PBIterSolvGMRES.h"
#include "PBPreconditioner.h"

namespace PHSPACE
{
PBIterSolvGMRES::PBIterSolvGMRES(const std::string &name):PBIterSolv(name)
{

}

PBIterSolvGMRES::~PBIterSolvGMRES()
{

}

void PBIterSolvGMRES::GMRESUpdate(RDouble **H, RDouble *x, RDouble *beta, std::vector<RDouble *> &Q, int iteration, int nCells)
{
    int lupe;
    for (lupe = iteration - 1; lupe >= 0; --lupe)
    {
        beta[lupe] = beta[lupe] / H[lupe][lupe];

        for (int innerlupe = lupe - 1; innerlupe >= 0; --innerlupe)
        {
            beta[innerlupe] -= beta[lupe] * H[lupe][innerlupe];
        }
    }
    for (lupe = 0; lupe < iteration; lupe++)
    {
        for (int iCell = 0; iCell < nCells; ++iCell)
        {
            Q[lupe][iCell] = Q[lupe][iCell] * beta[lupe];
        }

        for (int iCell = 0; iCell < nCells; ++iCell)
        {
            x[iCell] = Q[lupe][iCell] + x[iCell];
        }
    }
}

void PBIterSolvGMRES::Run(Region &region, vector<int> &zoneId, string varStr, string apStr, string ajUpperStr, string ajLowerStr,
                          string ajbUpperStr, string ajbLowerStr, string scStr, string rsStr, int maxSweep, RDouble tol)
{
    for (int iZone = 0; iZone < zoneId.size(); ++ iZone)
    {
        Zone *zone = region.GetZone(zoneId[iZone]);

        if (!zone)
        {
            continue;
        }

        UnstructGrid *grid = UnstructGridCast(zone->GetGeometry()->GetGrid());
        int nTotalCells = grid->GetNTotalCell();
        int nTotalFace = grid->GetNTotalFace();
        int nBoundFace = grid->GetNBoundFace();
        int nInteriorFace = nTotalFace - nBoundFace;

        int *leftCellOfFace = grid->GetLeftCellOfFace();
        int *rightCellOfFace = grid->GetRightCellOfFace();

        RDouble *apStr1 = reinterpret_cast<RDouble *>(grid->GetDataPtr(apStr));
        RDouble *ajUpperStr1 = reinterpret_cast<RDouble *>(grid->GetDataPtr(ajUpperStr));
        RDouble *ajLowerStr1 = reinterpret_cast<RDouble *>(grid->GetDataPtr(ajLowerStr));
        RDouble *scStr1 = reinterpret_cast<RDouble *>(grid->GetDataPtr(scStr));
        RDouble *rsStr1 = reinterpret_cast<RDouble *>(grid->GetDataPtr(rsStr));
        //! For the time being, it is only for a single zone.
        Calc(nInteriorFace, nTotalCells,
            grid->GetLeftCellOfFace() + grid->GetNBoundFace(),
            grid->GetRightCellOfFace() + grid->GetNBoundFace(),
            apStr1, ajUpperStr1, ajLowerStr1, scStr1, rsStr1, varStr, maxSweep, tol);
    }
}

void PBIterSolvGMRES::Calc(int nFaces, int nCells, int *leftCellOfFace, int *rightCellOfFace,RDouble *ap,RDouble *ajLeft,
                           RDouble *ajRight, RDouble *sc, RDouble *rs, string varStr, int maxSweep, RDouble tol)
{
    int iter = 0;
    RDouble res0, res;
    int krylovDemension = maxSweep;
    RDouble *apTmp = new RDouble[6 * nCells];
    RDouble *r = apTmp + nCells;
    RDouble *residual = r + nCells;
    RDouble *res_n = residual + nCells;
    RDouble *v = res_n + nCells;
    RDouble *q = v + nCells;
    RDouble *y = new RDouble[(krylovDemension + 1) * 2];
    RDouble *beta = y + krylovDemension + 1;
    RDouble **H = NewPointer2<RDouble>(krylovDemension, krylovDemension + 1);
    RDouble **givens = NewPointer2<RDouble>(krylovDemension, 2);

    for (int i = 0; i < krylovDemension; ++ i)
    {
        beta[i] = 0.0;
        for (int j = 0; j < 2; ++ j)
        {
            givens[i][j] = 0.0;
        }
    }
    beta[krylovDemension] = 0.0;
    std::vector<RDouble*> Q(krylovDemension + 1);

    preconditioner_->PreconditioningProxy(nFaces, nCells, leftCellOfFace, rightCellOfFace, ap, ajLeft, ajRight, apTmp);
    preconditioner_->ForwardAndBackward(nFaces, nCells, leftCellOfFace, rightCellOfFace, apTmp, ajLeft, ajRight, sc, r);

    RDouble norm_b = 0.0;
    for (int iCell = 0; iCell < nCells; ++ iCell)
    {
        norm_b += r[iCell] * r[iCell];
    }
    norm_b = sqrt(norm_b);

    res = norm_b;
    res0 = norm_b;
    if (norm_b < 1e-5)
    {
        norm_b = 1.0;
    }
    beta[0] = res0;

    RDouble one_norm = 1.0 / res0;

    for (int iCell = 0; iCell < nCells; ++iCell)
    {
        q[iCell] = r[iCell] * one_norm;
    }

    Q[0] = NewPointer<RDouble>(nCells);

    for (int iCell = 0; iCell < nCells; ++ iCell)
    {
        Q[0][iCell] = q[iCell];
    }

    while ((iter < krylovDemension) && (res > tol * res0))
    {
        Q[iter + 1] = NewPointer<RDouble>(nCells);

        Matrix(nFaces, nCells, leftCellOfFace, rightCellOfFace, ap, ajLeft, ajRight, q, v);
        preconditioner_->ForwardAndBackward(nFaces, nCells, leftCellOfFace, rightCellOfFace, apTmp, ajLeft, ajRight, v, v);

        for (int j = 0; j < iter + 1; j++)
        {
            RDouble ans = 0.0;
            for (int iCell = 0; iCell < nCells; ++ iCell)
            {
                ans += Q[j][iCell] * v[iCell];
            }
            H[iter][j] = ans;

            RDouble coe = -1.0 * H[iter][j];

            for (int iCell = 0; iCell < nCells; ++ iCell)
            {
                v[iCell] = v[iCell] + coe * Q[j][iCell];
            }
        }

        RDouble ans = 0.0;
        for (int iCell = 0; iCell < nCells; ++ iCell)
        {
            ans += v[iCell] * v[iCell];
        }
        H[iter][iter + 1] = sqrt(ans);
        RDouble coe = 1.0 / H[iter][iter + 1];

        for (int iCell = 0; iCell < nCells; ++ iCell)
        {
            q[iCell] = v[iCell] * coe;
        }

        for (int i = 0; i < nCells; ++ i)
        {
            Q[iter + 1][i] = q[i];
        }

        //! givens rotation for single rhs.
        for (int j = 0; j < iter; j++)
        {
            Rot(&H[iter][j], &H[iter][j + 1], &givens[j][0], &givens[j][1]);
        }

        Rotg(&H[iter][iter], &H[iter][iter + 1], &givens[iter][0], &givens[iter][1]);

        //! givens rotation for residual vector.
        Rot(&beta[iter], &beta[iter + 1], &givens[iter][0], &givens[iter][1]);

        res = abs(beta[iter + 1]);
        iter++;
    }
    GMRESUpdate(H, rs, beta, Q, iter, nCells);

    delete[] apTmp;
    delete[] y;
    DelPointer2(H);
    DelPointer2(givens);
    for (int i = 0; i <= iter; ++ i)
    {
        DelPointer(Q[i]);
    }
    std::vector<RDouble *>().swap(Q);

    std::cout << setiosflags(ios::left);
    std::cout << setprecision(5);
    std::cout << setiosflags(ios::scientific);
    std::cout << setw(12);
    std::cout << "GMRES: " << varStr << ",  SWEEP = " << iter << '\t' << "RES0 = " << res0 << '\t' << "RES = " << res << std::endl;
}

DEFINE_PBSOLVER_REGISTER_FUNCTION(PBIterSolv, PBIterSolvGMRES, ITERSOLV_GMRES);
}