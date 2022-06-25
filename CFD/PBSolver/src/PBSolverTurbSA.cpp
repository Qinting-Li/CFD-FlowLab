#include <iostream>
#include <memory>
#include <cmath>
#include "PBSolverTurbSA.h"
#include "PBBoundaryTurb.h"
#include "PBBoundaryScalar.h"
#include "PBGrad.h"

using namespace std;

namespace PHSPACE
{
PBSolverTurbSA::PBSolverTurbSA(const string &name):PBSolverScalar(name), delatv(0.6666),
    Cb2(0.622), urfMu(1.0),Cv1(7.1), Cb1(0.1355),kappa(0.4187),Cw2(0.3),Cw3(2)
{
    urfMu = ObjGetDoubleParaFromDB("urfMu", 1.0);
}

PBSolverTurbSA::~PBSolverTurbSA()
{

}

void PBSolverTurbSA::InitMemory()
{
    int nSolvers = GlobalDataBase::GetSizeOfData("TURB_SA");
    SetNSolvers(nSolvers);
    SetControllerName("TURB_SA");

    PBSolverScalar::InitMemory();

    AllocateMemory();
    InitField();

    diffCoe = mu;
}

//! Each solver applies for variable memory.
void PBSolverTurbSA::AllocateMemory()
{
    PBSolverScalar::AllocateMemory();

    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    visl = ObjAllocateMemory("visl", nTotal, 0.0, grid);
    vist = ObjAllocateMemory("vist", nTotal, 0.0, grid);
    gen = ObjAllocateMemory("gen", nTotal, 0.0, grid);
    des = ObjAllocateMemory("des", nTotal, 0.0, grid);
    gra = ObjAllocateMemory("gra", nTotal, 0.0, grid);
    fv1 = ObjAllocateMemory("fv1", nTotal, 0.0, grid);
    dkdx = ObjAllocateMemory("dkdx", nTotal, 0.0, grid);
    dkdy = ObjAllocateMemory("dkdy", nTotal, 0.0, grid);
    dkdz = ObjAllocateMemory("dkdz", nTotal, 0.0, grid);
    WD = ObjAllocateMemory("WD", nTotal, 0.0, grid);
}
//! The initialization of each solver depends on the associated variables of other solvers.
void PBSolverTurbSA::IncompressibleInitial()
{
    int nTotalCell = grid->GetNTotalCell();
    int nTotalFace = grid->GetNTotalFace();
    int nBoundFace = grid->GetNBoundFace();
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    k = reinterpret_cast<RDouble*>(grid->GetDataPtr("kinetic"));

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        RDouble chi = rho[iCell] * k[iCell] / diffCoe[iCell];
        fv1[iCell] = (pow(chi, 3)) / ((pow(chi, 3)) + (pow(Cv1, 3)));
        vist[iCell] = rho[iCell] * k[iCell] * fv1[iCell];
        visl[iCell] = diffCoe[iCell];
        diffCoe[iCell] = visl[iCell] + vist[iCell];//! Initial momentum equation diffusion coefficient
    }

    UpdateBoundary();
}

void PBSolverTurbSA::Discrt()
{
    int nTotalCell = grid->GetNTotalCell();
    int nTotalFace = grid->GetNTotalFace();
    int nBoundFace = grid->GetNBoundFace();
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    CalcGenAndDesTerm(grid);
    CalcGraTerm(grid);
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        diffCoe[iCell] = (visl[iCell] + rho[iCell] * k[iCell]) / delatv;
    }

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int le = leftCellOfFace[iFace];
        diffCoeb[iFace] = (visl[le] + rho[le] * k[le]) / delatv;
    }
    PBSolverScalar::Discrt();

}

void PBSolverTurbSA::CalcGenAndDesTerm(UnstructGrid* grid)
{
    int nTotalCell = grid->GetNTotalCell();
    //! Calculate the gradient of U, V and W.
    RDouble *u = reinterpret_cast<RDouble *>(grid->GetDataPtr("U"));
    RDouble *v = reinterpret_cast<RDouble *>(grid->GetDataPtr("V"));
    RDouble *w = reinterpret_cast<RDouble *>(grid->GetDataPtr("W"));
    RDouble *fv1 = reinterpret_cast<RDouble *>(grid->GetDataPtr("fv1"));
    RDouble *gen = reinterpret_cast<RDouble *>(grid->GetDataPtr("gen"));
    RDouble *des = reinterpret_cast<RDouble *>(grid->GetDataPtr("des"));
    RDouble *gra = reinterpret_cast<RDouble *>(grid->GetDataPtr("gra"));
    RDouble *WD = reinterpret_cast<RDouble *>(grid->GetDataPtr("WD"));
    RDouble *dudx = reinterpret_cast<RDouble *>(grid->GetDataPtr("dudx"));
    RDouble *dudy = reinterpret_cast<RDouble *>(grid->GetDataPtr("dudy"));
    RDouble *dudz = reinterpret_cast<RDouble *>(grid->GetDataPtr("dudz"));
    RDouble *dvdx = reinterpret_cast<RDouble *>(grid->GetDataPtr("dvdx"));
    RDouble *dvdy = reinterpret_cast<RDouble *>(grid->GetDataPtr("dvdy"));
    RDouble *dvdz = reinterpret_cast<RDouble *>(grid->GetDataPtr("dvdz"));
    RDouble *dwdx = reinterpret_cast<RDouble *>(grid->GetDataPtr("dwdx"));
    RDouble *dwdy = reinterpret_cast<RDouble *>(grid->GetDataPtr("dwdy"));
    RDouble *dwdz = reinterpret_cast<RDouble *>(grid->GetDataPtr("dwdz"));

    PBGrad *gradCalc = PBGradFactory::CreatePBGrad(ObjGetStrParaFromDB(gradCalcConfigName, "GAUSS"));
    gradCalc->Calc(grid, "U", "dudx", "dudy", "dudz");
    gradCalc->Calc(grid, "V", "dvdx", "dvdy", "dvdz");
    gradCalc->Calc(grid, "W", "dwdx", "dwdy", "dwdz");

    //! SA Wall treatment.
    RDouble* wDist = reinterpret_cast<RDouble *>(grid->GetWallDist());

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {

        RDouble wd = wDist[iCell];
        WD[iCell] = wd;

        RDouble _dudy = dudy[iCell];
        RDouble _dudz = dudz[iCell];
        RDouble _dvdx = dvdx[iCell];
        RDouble _dvdz = dvdz[iCell];
        RDouble _dwdx = dwdx[iCell];
        RDouble _dwdy = dwdy[iCell];

        //! Turbulent viscosity produces terms.
        RDouble chi = rho[iCell] * k[iCell] / visl[iCell];
        fv1[iCell] = (pow(chi, 3)) / ((pow(chi, 3)) + (pow(Cv1, 3)));

        RDouble fv2 = 1 - chi / (1 + chi * fv1[iCell]);
        RDouble S = sqrt(pow((_dwdy - _dvdz), 2) + pow((_dudz - _dwdx), 2) + pow((_dvdx - _dudy), 2));
        RDouble S_0 = S + k[iCell] * fv2 / (kappa * kappa * wd * wd);
        RDouble S_ = MAX(S_0, 0.3 * S);

        gen[iCell] = Cb1 * rho[iCell]*S_* k[iCell];
        gen[iCell] = MAX(gen[iCell], 0.0);

        //! Reduced turbulence viscosity.
        RDouble Cw1 = Cb1 / (kappa * kappa) + (1 + Cb2) / delatv;
        RDouble r = MIN((k[iCell]) / (S_ * kappa * kappa * wd * wd), 10.0);

        RDouble g = r + Cw2 * (pow(r, 6) - r);
        RDouble fw = g * pow(((1+Cw3)/(pow(g,6)+pow(Cw3,6))),1.0/6.0);

        des[iCell] = Cw1 * rho[iCell] * fw * (k[iCell] / wd) * (k[iCell] / wd);
    }
}

void PBSolverTurbSA::CalcGraTerm(UnstructGrid *grid)
{
    int nTotalCell = grid->GetNTotalCell();

    RDouble *dkdx = reinterpret_cast<RDouble *>(grid->GetDataPtr("dkdx"));
    RDouble *dkdy = reinterpret_cast<RDouble *>(grid->GetDataPtr("dkdy"));
    RDouble *dkdz = reinterpret_cast<RDouble *>(grid->GetDataPtr("dkdz"));
    RDouble *k = reinterpret_cast<RDouble *>(grid->GetDataPtr("kinetic"));
    RDouble *gra = reinterpret_cast<RDouble *>(grid->GetDataPtr("gra"));

    PBGrad *gradCalc = PBGradFactory::CreatePBGrad(ObjGetStrParaFromDB(gradCalcConfigName, "GAUSS"));
    gradCalc->Calc(grid,"kinetic","dkdx","dkdy","dkdz");

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        RDouble _dkdx = dkdx[iCell];
        RDouble _dkdy = dkdy[iCell];
        RDouble _dkdz = dkdz[iCell];
        gra[iCell] = Cb2 * rho[iCell] * (_dkdx * _dkdx + _dkdy * _dkdy + _dkdz * _dkdz) / delatv;
    }

}

void PBSolverTurbSA::UpdateAfterIterloop()
{
    int nTotalCell = grid->GetNTotalCell();
    int nTotalFace = grid->GetNTotalFace();
    int nBoundFace = grid->GetNBoundFace();
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    RDouble *vistb = vist + nTotalCell;

    for (int iCell = 0; iCell < nTotalCell; ++iCell)
    {
        RDouble vistold = vist[iCell];
        vist[iCell] = vistold + (rho[iCell] * k[iCell] * fv1[iCell] - vistold) * urfMu;
        diffCoe[iCell] = visl[iCell] + vist[iCell];
    }

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int le = leftCellOfFace[iFace];
        int re = rightCellOfFace[iFace];

        vistb[iFace] = vist[le];
        diffCoeb[iFace] = visl[le];//! for SA wall,vistb=0.
    }
}

void PBSolverTurbSA::ComputePostVisualVariables(Post_Visual *postVisualization)
{
    UnstructGrid *grid = UnstructGridCast(postVisualization->GetGrid());

    RDouble *gen = reinterpret_cast<RDouble *> (grid->GetDataPtr("mu"));

    postVisualization->UpdateVisualNodeVarPtr("VIS", vist);
}

void PBSolverTurbSA::OutputResidual(ofstream &file, int flag)
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

DEFINE_PBSOLVER_REGISTER_FUNCTION(PBSolver, PBSolverTurbSA, SOLVER_TURB_SA);
}