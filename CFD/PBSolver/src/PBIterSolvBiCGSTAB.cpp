#include "PBIterSolvBiCGSTAB.h"
#include "PBPreconditioner.h"

namespace PHSPACE
{
PBIterSolvBiCGSTAB::PBIterSolvBiCGSTAB(const string &name):PBIterSolv(name)
{

}

PBIterSolvBiCGSTAB::~PBIterSolvBiCGSTAB()
{

}

void PBIterSolvBiCGSTAB::Run(Region &region, vector<int> &zoneId, string varStr, string apStr, string ajUpperStr, string ajLowerStr,
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
        int nTotalCell = grid->GetNTotalCell();
        int nTotalFace = grid->GetNTotalFace();
        int nBoundFace = grid->GetNBoundFace();
        int nInteriorFace = nTotalFace - nBoundFace;

        int *leftCellOfFace = grid->GetLeftCellOfFace();
        int *rightCellOfFace = grid->GetRightCellOfFace();

        int *leftCellOfFaceb = leftCellOfFace + nBoundFace;
        int *rightCellOfFaceb = rightCellOfFace + nBoundFace;

        RDouble *apStr1 = reinterpret_cast<RDouble *>(grid->GetDataPtr(apStr));
        RDouble *ajUpperStr1 = reinterpret_cast<RDouble *>(grid->GetDataPtr(ajUpperStr));
        RDouble *ajLowerStr1 = reinterpret_cast<RDouble *>(grid->GetDataPtr(ajLowerStr));
        RDouble *scStr1 = reinterpret_cast<RDouble *>(grid->GetDataPtr(scStr));
        RDouble *rsStr1 = reinterpret_cast<RDouble *>(grid->GetDataPtr(rsStr));

        //! For the time being, it is only for a single zone.
        Calc(nInteriorFace, nTotalCell, leftCellOfFaceb, rightCellOfFaceb, apStr1, ajUpperStr1, ajLowerStr1, scStr1, rsStr1, varStr, maxSweep, tol);
    }
}

void PBIterSolvBiCGSTAB::Calc(int nFaces, int nCells, int *leftCellOfFace, int *rightCellOfFace, RDouble *ap, RDouble *ajLeft,
                              RDouble *ajRight, RDouble *sc, RDouble *rs, string varStr, int maxSweep, RDouble tol)
{
    RDouble res0, res;
    int iter = 0;

    //! Initialize
    RDouble alpha, beta, omega;
    RDouble *r = new RDouble[8 * nCells];
    RDouble *p = r + nCells;
    RDouble *v = p + nCells;
    RDouble *t = v + nCells;
    RDouble *s = t + nCells;
    RDouble *rHat = s + nCells;
    RDouble *r_temp = rHat + nCells;
    RDouble *apTmp = r_temp + nCells;

    preconditioner_->PreconditioningProxy(nFaces, nCells, leftCellOfFace, rightCellOfFace, ap, ajLeft, ajRight, apTmp);
    preconditioner_->ForwardAndBackward(nFaces, nCells, leftCellOfFace, rightCellOfFace, apTmp, ajLeft, ajRight, sc, r);

    for (int iCell = 0; iCell < nCells; ++iCell)
    {
        p[iCell] = r[iCell];
        rHat[iCell] = r[iCell];
    }

    res0 = Dot(r, r, nCells);
    res = sqrt(res0);
    res0 = res;
    iter = 0;
    while (res > tol * res0 && iter < maxSweep)
    {
        iter++;
        Matrix(nFaces, nCells, leftCellOfFace, rightCellOfFace, ap, ajLeft, ajRight, p, v);
        preconditioner_->ForwardAndBackward(nFaces, nCells, leftCellOfFace, rightCellOfFace, apTmp, ajLeft, ajRight, v, v);

        alpha = Dot(rHat, r, nCells) / Dot(rHat, v, nCells);

        for (int iCell = 0; iCell < nCells; ++ iCell)
        {
            s[iCell] = r[iCell] - alpha * v[iCell];
        }
        Matrix(nFaces, nCells, leftCellOfFace, rightCellOfFace, ap, ajLeft, ajRight, s, t);
        preconditioner_->ForwardAndBackward(nFaces, nCells, leftCellOfFace, rightCellOfFace, apTmp, ajLeft, ajRight, t, t);

        omega = Dot(t, s, nCells) / Dot(t, t, nCells);

        for (int iCell = 0; iCell < nCells; ++ iCell)
        {
            rs[iCell] = rs[iCell] + alpha * p[iCell] + omega * s[iCell];
            r_temp[iCell] = r[iCell];
            r[iCell] = s[iCell] - omega * t[iCell];
        }

        res = Dot(r, r, nCells);
        res = sqrt(res);

        if (res > tol * res0)
        {
            beta = Dot(rHat, r, nCells);
            beta = beta / Dot(rHat, r_temp, nCells);
            beta = beta * (alpha / omega);
            for (int iCell = 0; iCell < nCells; ++iCell)
            {
                p[iCell] = r[iCell] + beta * (p[iCell] - omega * v[iCell]);
            }
        }
        else
        {
            continue;
        }
    }
    delete[] r;

    std::cout << setiosflags(ios::left);
    std::cout << setprecision(5);
    std::cout << setiosflags(ios::scientific);
    std::cout << setw(12);
    std::cout << "BiCGSTAB: " << varStr << ",  SWEEP = " << iter << '\t' << "RES0 = " << res0 << '\t' << "RES = " << res << std::endl;
}

DEFINE_PBSOLVER_REGISTER_FUNCTION(PBIterSolv, PBIterSolvBiCGSTAB, ITERSOLV_BiCGSTAB);
}