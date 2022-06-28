#include <iostream>
#include <memory>
#include <cmath>
#include "PBSolverTurbK.h"
#include "PBBoundaryTurb.h"
#include "PBBoundaryScalar.h"
#include "PBGrad.h"

using namespace std;

namespace PHSPACE
{
PBSolverTurbK::PBSolverTurbK(const string &name):PBSolverScalar(name), prandtl(1.0), cmu(0.09), urfMu(1.0)
{
    urfMu = ObjGetDoubleParaFromDB("urfMu", 1.0);
}

PBSolverTurbK::~PBSolverTurbK()
{

}

void PBSolverTurbK::InitMemory()
{
    int nSolvers = GlobalDataBase::GetSizeOfData("TURB_K");
    SetNSolvers(nSolvers);
    SetControllerName("TURB_K");

    PBSolverScalar::InitMemory();

    AllocateMemory();
    InitField();
}

void PBSolverTurbK::AllocateMemory()
{
    PBSolverScalar::AllocateMemory();
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    visl = ObjAllocateMemory("visl", nTotal, 0.0, grid);
    vist = ObjAllocateMemory("vist", nTotal, 0.0, grid);
    vistb = vist + nTotalCell;
    gen = ObjAllocateMemory("gen", nTotal, 0.0, grid);
    tauWall = ObjAllocateMemory("tauWall", nTotalCell, 0.0, grid);
    yplus = ObjAllocateMemory("yplus", nTotal, 0.0, grid);
}

void PBSolverTurbK::IncompressibleInitial()
{
    int nTotalCell = grid->GetNTotalCell();
    int nTotalFace = grid->GetNTotalFace();
    int nBoundFace = grid->GetNBoundFace();
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    k = reinterpret_cast<RDouble *>(grid->GetDataPtr("kinetic"));
    epsilon = reinterpret_cast<RDouble *>(grid->GetDataPtr("epsilon"));

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        vist[iCell] = rho[iCell] * cmu * k[iCell] * k[iCell] / epsilon[iCell];
        visl[iCell] = mu[iCell];
        mu[iCell] = visl[iCell] + vist[iCell];
    }

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int le = leftCellOfFace[iFace];

        vistb[iFace] = vist[le];
        mub[iFace] = mu[le];

    }

    UpdateBoundary();
}

void PBSolverTurbK::Discrt()
{
    int nTotalCell = grid->GetNTotalCell();
    int nTotalFace = grid->GetNTotalFace();
    int nBoundFace = grid->GetNBoundFace();
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    CalcGenTerm();

    UpdateWallGen();

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        diffCoe[iCell] = visl[iCell] + vist[iCell] / prandtl;
    }
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int le = leftCellOfFace[iFace];
        mub[iFace] = mu[le];
        vistb[iFace] = vist[le];
        diffCoeb[iFace] = visl[le] + vistb[iFace] / prandtl;
    }
    PBSolverScalar::Discrt();
}

void PBSolverTurbK::CalcGenTerm()
{
    int nTotalCell = grid->GetNTotalCell();

    //! Calculate the gradient of U, V and W.
    RDouble *u = reinterpret_cast<RDouble*>(grid->GetDataPtr("U"));
    RDouble *v = reinterpret_cast<RDouble*>(grid->GetDataPtr("V"));
    RDouble *w = reinterpret_cast<RDouble*>(grid->GetDataPtr("W"));
    RDouble *dudx = reinterpret_cast<RDouble*>(grid->GetDataPtr("dudx"));
    RDouble *dudy = reinterpret_cast<RDouble*>(grid->GetDataPtr("dudy"));
    RDouble *dudz = reinterpret_cast<RDouble*>(grid->GetDataPtr("dudz"));
    RDouble *dvdx = reinterpret_cast<RDouble*>(grid->GetDataPtr("dvdx"));
    RDouble *dvdy = reinterpret_cast<RDouble*>(grid->GetDataPtr("dvdy"));
    RDouble *dvdz = reinterpret_cast<RDouble*>(grid->GetDataPtr("dvdz"));
    RDouble *dwdx = reinterpret_cast<RDouble*>(grid->GetDataPtr("dwdx"));
    RDouble *dwdy = reinterpret_cast<RDouble*>(grid->GetDataPtr("dwdy"));
    RDouble *dwdz = reinterpret_cast<RDouble*>(grid->GetDataPtr("dwdz"));
    RDouble *vist = reinterpret_cast<RDouble*>(grid->GetDataPtr("vist"));
    RDouble *gen = reinterpret_cast<RDouble*>(grid->GetDataPtr("gen"));
    RDouble* yplus = reinterpret_cast<RDouble*>(grid->GetDataPtr("yplus"));

    PBGrad *gradCalc = PBGradFactory::CreatePBGrad(ObjGetStrParaFromDB(gradCalcConfigName, "GAUSS"));
    gradCalc->Calc(grid, "U", "dudx", "dudy", "dudz");
    gradCalc->Calc(grid, "V", "dvdx", "dvdy", "dvdz");
    gradCalc->Calc(grid, "W", "dwdx", "dwdy", "dwdz");

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        RDouble _dudx = dudx[iCell];
        RDouble _dudy = dudy[iCell];
        RDouble _dudz = dudz[iCell];
        RDouble _dvdx = dvdx[iCell];
        RDouble _dvdy = dvdy[iCell];
        RDouble _dvdz = dvdz[iCell];
        RDouble _dwdx = dwdx[iCell];
        RDouble _dwdy = dwdy[iCell];
        RDouble _dwdz = dwdz[iCell];
        RDouble divv = _dudx + _dvdy + _dwdz;

        gen[iCell] = vist[iCell] * (2.0 * (pow(_dudx, 2) + pow(_dvdy, 2) + pow(_dwdz, 2)) +
            pow((_dvdx + _dudy), 2) + pow((_dwdx + _dudz), 2) + pow((_dwdy + _dvdz), 2));

        gen[iCell] = MAX(gen[iCell], 0.0);
    }

}

void PBSolverTurbK::UpdateAfterIterloop()
{
    int nTotalCell = grid->GetNTotalCell();
    int nTotalFace = grid->GetNTotalFace();
    int nBoundFace = grid->GetNBoundFace();
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    RDouble *mu = reinterpret_cast<RDouble *>(grid->GetDataPtr("mu"));
    RDouble *mub = mu + nTotalCell;
    RDouble *vistb = vist + nTotalCell;

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        RDouble vistold = vist[iCell];
        vist[iCell] = vistold + (rho[iCell] * cmu * k[iCell] * k[iCell] / epsilon[iCell] - vistold) * urfMu;
        mu[iCell] = visl[iCell] + vist[iCell];
    }

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegionUnstruct = unstructBCSet->GetnBCRegion();

    for (int iBCRegionUnstruct = 0; iBCRegionUnstruct < nBCRegionUnstruct; iBCRegionUnstruct++)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegionUnstruct);
        int bcType = bcRegion->GetBCType();
        vector<int> *faceIndex = bcRegion->GetFaceIndex();
        Data_Param *bcData = bcRegion->GetBCParamDataBase();
        static_cast<PBBoundary*>(boundCalc[iBCRegionUnstruct])->UpdateVistb(grid, ScalarName(),faceIndex, bcData);
    }
}

void PBSolverTurbK::ComputePostVisualVariables(Post_Visual *postVisualization)
{
    UnstructGrid *grid = UnstructGridCast(postVisualization->GetGrid());

    RDouble *k = reinterpret_cast<RDouble *> (grid->GetDataPtr("kinetic"));
    RDouble *e = reinterpret_cast<RDouble *> (grid->GetDataPtr("epsilon"));
    RDouble *visl = reinterpret_cast<RDouble *> (grid->GetDataPtr("visl"));
    RDouble *vist = reinterpret_cast<RDouble *> (grid->GetDataPtr("vist"));
    RDouble *gen = reinterpret_cast<RDouble *> (grid->GetDataPtr("gen"));
    RDouble *yplus = reinterpret_cast<RDouble *>(grid->GetDataPtr("yplus"));
    RDouble *wDist = reinterpret_cast<RDouble *>(grid->GetWallDist());

    postVisualization->UpdateVisualNodeVarPtr("VIS", vist);
    postVisualization->UpdateVisualNodeVarPtr("TE", k);
    postVisualization->UpdateVisualNodeVarPtr("ED", e);
    postVisualization->UpdateVisualNodeVarPtr("GENK", gen);
    postVisualization->UpdateVisualNodeVarPtr("YPLUS", yplus);
    postVisualization->UpdateVisualNodeVarPtr("WD", wDist);
}

void PBSolverTurbK::OutputResidual(ofstream &file, int flag)
{
    file << setiosflags(ios::left);
    file << setprecision(5);
    file << setiosflags(ios::scientific);

    if (flag == 1)
    {
        file << " Res TE" << '\t';
    }
    else
    {
        RDouble teResNow = 0.0;
        RDouble teRes = 0.0;

        grid->GetData(ScalarName() + "ResNow", &teRes, PHDOUBLE, 1);

        teResNow = teRes;

        file << teResNow << '\t';
    }
}
void PBSolverTurbK::UpdateWallGen()
{
    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegionUnstruct = unstructBCSet->GetnBCRegion();

    for (int iBCRegionUnstruct = 0; iBCRegionUnstruct < nBCRegionUnstruct; iBCRegionUnstruct++)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegionUnstruct);
        int bcType = bcRegion->GetBCType();
        vector<int> *faceIndex = bcRegion->GetFaceIndex();
        Data_Param *bcData = bcRegion->GetBCParamDataBase();

        static_cast<PBBoundaryScalar*>(boundCalc[iBCRegionUnstruct])->WallModification(grid, "gen", faceIndex, bcData);
    }
}

RDouble PBSolverTurbK::CacluteWallDistance(UnstructGrid *grid,int ic, int iFace)
{
    RDouble *xfn = grid->GetFaceNormalX();
    RDouble *yfn = grid->GetFaceNormalY();
    RDouble *zfn = grid->GetFaceNormalZ();
    RDouble *xfc = grid->GetFaceCenterX();
    RDouble *yfc = grid->GetFaceCenterY();
    RDouble *zfc = grid->GetFaceCenterZ();
    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    RDouble faceNormalX = xfn[iFace];
    RDouble faceNormalY = yfn[iFace];
    RDouble faceNormalZ = zfn[iFace];
    RDouble faceCenterX = xfc[iFace];
    RDouble faceCenterY = yfc[iFace];
    RDouble faceCenterZ = zfc[iFace];
    RDouble cellCenterX = xcc[ic];
    RDouble cellCenterY = ycc[ic];
    RDouble cellCenterZ = zcc[ic];
    RDouble dx = (faceCenterX - cellCenterX) * faceNormalX;
    RDouble dy = (faceCenterY - cellCenterY) * faceNormalY;
    RDouble dz = (faceCenterZ - cellCenterZ) * faceNormalZ;
    RDouble wdist = sqrt(pow(dx,2) + pow(dy, 2) + pow(dz, 2));
    return wdist;
}
DEFINE_PBSOLVER_REGISTER_FUNCTION(PBSolver, PBSolverTurbK, SOLVER_TURB_K);
}