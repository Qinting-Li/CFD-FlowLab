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
    PBSolverTurbK::PBSolverTurbK(const string &name)
        :PBSolverScalar(name), prandtl(1.0), cmu(0.09), urfMu(1.0)
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

        int nTotal = grid->GetNTotalCell() + grid->GetNBoundFace();

        visl = ObjAllocateMemory("visl", nTotal, 0.0, grid);
        vist = ObjAllocateMemory("vist", nTotal, 0.0, grid);
        gen = ObjAllocateMemory("gen", nTotal, 0.0, grid);
    }

    void PBSolverTurbK::ReinitField()
    {
        k = reinterpret_cast<RDouble *>(grid->GetDataPtr("kinetic"));
        epsilon = reinterpret_cast<RDouble *>(grid->GetDataPtr("epsilon"));

        int *leftCellOfFace = grid->GetLeftCellOfFace();

        for (int iCell = 0; iCell < grid->GetNTotalCell(); ++ iCell)
        {
            vist[iCell] = rho[iCell] * cmu * k[iCell] * k[iCell] / epsilon[iCell];
            visl[iCell] = mu[iCell];
            mu[iCell] = visl[iCell] + vist[iCell];
        }

        for (int iFace = 0; iFace < grid->GetNBoundFace(); ++ iFace)
        {
            int le = leftCellOfFace[iFace];

            mu[le] = visl[le];
            mub[iFace] = mu[le];

        }

        UpdateBoundary();
    }

    void PBSolverTurbK::Discrt()
    {
        int *leftCellOfFace = grid->GetLeftCellOfFace();

        CalcGenTerm();

        for (int iCell = 0; iCell < grid->GetNTotalCell(); ++ iCell)
        {
            mu[iCell] = visl[iCell] + vist[iCell] / prandtl;
        }

        for (int iFace = 0; iFace < grid->GetNBoundFace(); ++ iFace)
        {
            int iCell = leftCellOfFace[iFace];
            mub[iFace] = visl[iCell] + vist[iFace] / prandtl;
        }

        diffCoe = mu;
        diffCoeb = mub;
        PBSolverScalar::Discrt();
    }

    void PBSolverTurbK::CalcGenTerm()
    {
        //! Calculate the gradient of U, V and W;
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

        PBGrad *gradCalc = PBGradFactory::CreatePBGrad(ObjGetStrParaFromDB(gradCalcConfigName, "GAUSS"));
        gradCalc->Calc(grid, "U", "dudx", "dudy", "dudz");
        gradCalc->Calc(grid, "V", "dvdx", "dvdy", "dvdz");
        gradCalc->Calc(grid, "W", "dwdx", "dwdy", "dwdz");

        for (int iCell = 0; iCell < grid->GetNTotalCell(); ++ iCell)
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

            gen[iCell] = vist[iCell] * (2.0 * (pow(_dudx, 2) + pow(_dvdy, 2) + pow(_dwdz, 2)) +
                pow((_dvdx + _dudy), 2) + pow((_dwdx + _dudz), 2) + pow((_dwdy + _dvdz), 2));

            gen[iCell] = MAX(gen[iCell], 0.0);
        }

        UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
        int nBCRegionUnstruct = unstructBCSet->GetnBCRegion();

        for (int iBCRegionUnstruct = 0; iBCRegionUnstruct < nBCRegionUnstruct; iBCRegionUnstruct++)
        {
            UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegionUnstruct);
            int bcType = bcRegion->GetBCType();

            if (bcType < 0)
            {
                continue;
            }

            static_cast<PBBoundaryScalar *>(boundCalc[iBCRegionUnstruct])->WallModification
            (grid, "gen", bcRegion->GetFaceIndex(), bcRegion->GetBCParamDataBase());
        }
    }

    void PBSolverTurbK::UpdateAfterIterloop()
    {
        int nTotalCells = grid->GetNTotalCell();

        RDouble *vistb = vist + nTotalCells;

        for (int iCell = 0; iCell < grid->GetNTotalCell(); ++ iCell)
        {
            RDouble vistold = vist[iCell];
            vist[iCell] = vistold + (rho[iCell] * cmu * k[iCell] * k[iCell] / epsilon[iCell] - vistold) * urfMu;
            diffCoe[iCell] = visl[iCell] + vist[iCell];
        }

        //£¡correct £ºy+ÅÐ¶Ï
        for (int iFace = 0; iFace < grid->GetNBoundFace(); ++ iFace)
        {
            int le = grid->GetLeftCellOfFace()[iFace];
            int re = grid->GetRightCellOfFace()[iFace];

            vistb[iFace] = vist[le];
            diffCoeb[iFace] = visl[le] + vistb[iFace];
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

        postVisualization->UpdateVisualNodeVarPtr("VIS", vist);
        postVisualization->UpdateVisualNodeVarPtr("TE", k);
        postVisualization->UpdateVisualNodeVarPtr("ED", e);
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

    DEFINE_PBSOLVER_REGISTER_FUNCTION(PBSolver, PBSolverTurbK, SOLVER_TURB_K);
}