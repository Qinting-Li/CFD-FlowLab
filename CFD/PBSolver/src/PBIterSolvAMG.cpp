#include "PBIterSolvAMG.h"
#include "PBPreconditioner.h"

namespace PHSPACE
{
PBIterSolvAMG::PBIterSolvAMG(const string &name):PBIterSolv(name)
{
    this->groupsize = 2;
    this->entryNum = 0;
}

PBIterSolvAMG::~PBIterSolvAMG()
{

}

void PBIterSolvAMG::CSRsp(RDouble *ajLeft, RDouble *ajRight, int ncells, int nfaces, int *leftCellOfFace, int *rightCellOfFace, 
                               int *faceNumberOfEachCell, int **faceIndexOfEachCell, int boundFaces)
{
    entryNum = 0;
    int le, re, k;
    int *Counttmp = new int [ncells];
    int *faceNumberOfEachCell_ = new int [ncells];
    for (int iCell = 0; iCell < ncells; ++ iCell)
    {
        faceNumberOfEachCell_[iCell] = faceNumberOfEachCell[iCell];
        for (int iFace = 0; iFace < faceNumberOfEachCell_[iCell]; ++ iFace)
        {
            int fId = faceIndexOfEachCell[iCell][iFace];
            if (fId < boundFaces)
            {
                faceNumberOfEachCell_[iCell] -= 1;
            }
        }
    }

    for (int iCell = 0; iCell < ncells; ++ iCell)
    {
        entryNum += faceNumberOfEachCell_[iCell];
        Counttmp[iCell] = 0;
    }

    AMG_A = new RDouble [entryNum];
    IA = new int [ncells + 1];
    JA =new int [entryNum];

    for (int iCell = 0; iCell < ncells; ++ iCell)
    {
        IA[0] = 1;
        IA[iCell + 1] = IA[iCell] + faceNumberOfEachCell_[iCell];
    }

    for (int iFace = 0; iFace < nfaces; ++ iFace)
    {
        le = leftCellOfFace[iFace];
        re = rightCellOfFace[iFace];

        k = IA[le] + Counttmp[le];
        AMG_A[k - 1] = -ajLeft[iFace];
        JA[k - 1] = re + 1;
        Counttmp[le] += 1;

        k = IA[re] + Counttmp[re];
        AMG_A[k - 1] = -ajRight[iFace];
        JA[k - 1] = le + 1;
        Counttmp[re] += 1;
    }

    delete [] Counttmp;
    delete [] faceNumberOfEachCell_;
}

void PBIterSolvAMG::Run(Region &region, vector<int> &zoneId, string varStr, string apStr, string ajUpperStr, string ajLowerStr,
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

        faceNumberOfEachCell = grid->GetFaceNumberOfEachCell();
        faceIndexOfEachCell = grid->GetCell2Face();
        boundFaces = nBoundFace;

        //! For the time being, it is only for a single zone.
        Calc(nInteriorFace, nTotalCell, leftCellOfFaceb, rightCellOfFaceb, apStr1, ajUpperStr1, ajLowerStr1, scStr1, rsStr1, varStr, maxSweep, tol);
    }
}

void PBIterSolvAMG::Calc(int nfaces, int ncells, int *leftCellOfFace, int *rightCellOfFace, RDouble *ap, RDouble *ajLeft, RDouble *ajRight,
                         RDouble *sc, RDouble *rs,string varStr, int maxSweep, RDouble tol)
{
    int nfaces_ = nfaces;
    int ncells_ = ncells;
    int sweep = 0;
    RDouble res0 = 1.0;
    RDouble res1 = 1.0;
    RDouble *sp_ = const_cast<RDouble *>(ap);
    RDouble *sc_ = const_cast<RDouble *>(sc);
    CSRsp(ajLeft, ajRight, ncells, nfaces, leftCellOfFace, rightCellOfFace, faceNumberOfEachCell, faceIndexOfEachCell, boundFaces);
    PGAMG(&ncells_, &ncells_, &nfaces_, IA, JA, sp_, AMG_A, rs, sc_, &groupsize, &maxSweep, &tol, &sweep, &res0, &res1);
    delete [] AMG_A;
    delete [] IA;
    delete [] JA;

    std::cout << setiosflags(ios::left);
    std::cout << setprecision(5);
    std::cout << setiosflags(ios::scientific);
    std::cout << setw(12);
    std::cout << "AMG: " << varStr << ",  SWEEP = " << sweep << '\t' << "RES0 = " << res0 << '\t' << "RES = " << res1 << std::endl;
}

DEFINE_PBSOLVER_REGISTER_FUNCTION(PBIterSolv, PBIterSolvAMG, ITERSOLV_AMG);
}