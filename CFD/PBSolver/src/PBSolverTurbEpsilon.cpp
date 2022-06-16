#include <iostream>
#include "PBSolverTurbEpsilon.h"
#include "PBBoundaryScalar.h"

using namespace std;

namespace PHSPACE
{
    PBSolverTurbEpsilon::PBSolverTurbEpsilon(const string& name):PBSolverScalar(name), prandtl(1.3)
    {

    }

    PBSolverTurbEpsilon::~PBSolverTurbEpsilon()
    {

    }

    void PBSolverTurbEpsilon::InitMemory()
    {
        int nSolvers = GlobalDataBase::GetSizeOfData("TURB_EPSILON");
        SetNSolvers(nSolvers);
        SetControllerName("TURB_EPSILON");

        PBSolverScalar::InitMemory();

        AllocateMemory();
        InitField();
    }

    void PBSolverTurbEpsilon::Discrt()
    {
        RDouble *visl = reinterpret_cast<RDouble *>(grid->GetDataPtr("visl"));
        RDouble *vist = reinterpret_cast<RDouble *>(grid->GetDataPtr("vist"));
        int *leftCellOfFace = grid->GetLeftCellOfFace();

        for (int iCell = 0; iCell < grid->GetNTotalCell(); ++ iCell)
        {
            diffCoe[iCell] = visl[iCell] + vist[iCell] / prandtl;
        }

        for (int iFace = 0; iFace < grid->GetNBoundFace(); ++ iFace)
        {
            int le = leftCellOfFace[iFace];
            diffCoeb[iFace] = diffCoe[le];
        }

        //! Initialize the related variables of linear equations.
        InitMatrix();

        //! Calculate the diffusion term.
        CalcDiffusionTerm();

        //! Convection term calculation.
        CalcConvectionTerm();

        //! Source term calculation.
        CalcSourceTerm();

        //! Boundary calculation.
        CalcBoundaryConditions();

        //! Transient term calculation.
        CalcTransientTerm();

        // Form ap and sc.
        PrimeEquations();

        //! Calculate the residual.
        CalcRes();

        //! Relax ap.
        Relax();
    }

    void PBSolverTurbEpsilon::OutputResidual(ofstream &file, int flag)
    {
        file << setiosflags(ios::left);
        file << setprecision(5);
        file << setiosflags(ios::scientific);

        if (flag == 1)
        {
            file << " Res ED" << '\t';
        }
        else
        {
            RDouble edResNow = 0.0;
            RDouble edRes = 0.0;

            grid->GetData(ScalarName() + "ResNow", &edRes, PHDOUBLE, 1);

            edResNow = edRes;

            file << edResNow << '\t';
        }
    }

    DEFINE_PBSOLVER_REGISTER_FUNCTION(PBSolver, PBSolverTurbEpsilon, SOLVER_TURB_EPSILON);
}