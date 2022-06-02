#include "PBSolverEnergy.h"

namespace PHSPACE
{
    PBSolverEnergy::PBSolverEnergy(const string &name):PBSolverScalar(name)
    {

    }

    PBSolverEnergy::~PBSolverEnergy()
    {
    }

    void PBSolverEnergy::InitMemory()
    {
        string *phiName = new string[1];
        GlobalDataBase::GetData("ENERGY_SCALAR_NAME", phiName, PHSTRING, 1);
        SetScalarName(phiName[0]);
        SetGradientPhiName(phiName[0]);
        energyPrintName = new std::string[1];
        GlobalDataBase::GetData("energyPrintName", energyPrintName, PHDOUBLE, 1);
        SetPrintName(energyPrintName[0]);
        boundaryCalcConfigName = "energyBoundaryType";
        convCalcConfigName = "energyConvCalc";
        diffCalcConfigName = "energyDiffCalc";
        tranCalcConfigName = "energyTranCalc";
        sourceCalcConfigName = "energySourceCalc";
        gradCalcConfigName = "energyGradCalc";

        AllocateMemory();
        InitField();
    }

    void PBSolverEnergy::AllocateMemory()
    {
        PBSolverScalar::AllocateMemory();

        initk = ObjGetDoubleParaFromDB("initK", 0.0242);
        initCPg = ObjGetDoubleParaFromDB("initCPg", 1006.0);
        urfT = ObjGetDoubleParaFromDB("urfT", 1.0);
        urfH = ObjGetDoubleParaFromDB("urfH", 1.0);
        tMin = ObjGetDoubleParaFromDB("tMin", 1.0);
        tMax = ObjGetDoubleParaFromDB("tMax", 5000.0);

        int nTotalCells = grid->GetNTotalCell();
        int nBoundFaces = grid->GetNBoundFace();
        int nTotalFaces = grid->GetNTotalFace();
        int nTotalInFaces = nTotalFaces - nBoundFaces;
        unsigned int nTotal = nTotalCells + nBoundFaces;
        nTotalInFaces = nTotalInFaces;

        h = phi;
        hb = phib;

        initT = ObjGetDoubleParaFromDB("initT", 0.0);
        t = ObjAllocateMemory("T", nTotal, initT, grid);
        tb = t + nTotalCells;

        cpg = ObjAllocateMemory("cp", nTotal, initCPg, grid);
        cpgb = cpg + nTotalCells;
        grid->UpdateDataPtr("cpb", cpgb);

        k = ObjAllocateMemory("k", nTotal, initk, grid);
        kb = k + nTotalCells;
        grid->UpdateDataPtr("kb", kb);

        diffCoe = k;
        diffCoeb = k + nTotalCells;
    }

    void PBSolverEnergy::ReinitField()
    {
        UpdateBoundary();

        RDouble usq;
        int nTotalCell = grid->GetNTotalCell();
        int nTotalInFaces = grid->GetNTotalFace() - grid->GetNBoundFace();
        int *leftCellOfFace = grid->GetLeftCellOfFace();
        int *rightCellOfFace = grid->GetRightCellOfFace();
        int nTotal = nTotalCell + grid->GetNBoundFace();
        RDouble *u = reinterpret_cast<RDouble *>(grid->GetDataPtr("U"));
        RDouble *v = reinterpret_cast<RDouble *>(grid->GetDataPtr("V"));
        RDouble *w = reinterpret_cast<RDouble *>(grid->GetDataPtr("W"));

        for (int iCell = 0; iCell < nTotal; ++ iCell)
        {
            diffCoe[iCell] = k[iCell];

            usq = (u[iCell] * u[iCell] + v[iCell] * v[iCell] + w[iCell] * w[iCell]) / 2;
            h[iCell] = cpg[iCell] * t[iCell] + usq;
        }

        if (ObjGetIntParaFromDB("iunsteady", 0))
        {
            RDouble *hOld = reinterpret_cast<RDouble *>(grid->GetDataPtr(ScalarName() + "Old"));
            RDouble *hOldOld = reinterpret_cast<RDouble *>(grid->GetDataPtr(ScalarName() + "OldOld"));
            if (hOldOld == NULL)
            {
                for (int iCell = 0; iCell < nTotal; ++ iCell)
                {
                    hOld[iCell] = h[iCell];
                }
            }
            else
            {
                for (int iCell = 0; iCell < nTotal; ++ iCell)
                {
                    hOld[iCell] = h[iCell];
                    hOldOld[iCell] = h[iCell];
                }
            }
        }
    }

    void PBSolverEnergy::CalcDiffusionTerm()
    {
        SetGradientPhiName("T");
        PBSolverScalar::CalcDiffusionTerm();
    }

    void PBSolverEnergy::PrimeEquations()
    {
        phi = h;
        phib = hb;
        PBSolverScalar::PrimeEquations();
    }

    void PBSolverEnergy::Relax()
    {
        urf = urfH;
        PBSolverScalar::Relax();
    }

    void PBSolverEnergy::Update()
    {
        phi = h;
        phib = hb;
        PBSolverScalar::Update();
    }

    void PBSolverEnergy::UpdateProperties()
    {
        RDouble *u = reinterpret_cast<RDouble *>(grid->GetDataPtr("U"));
        RDouble *v = reinterpret_cast<RDouble *>(grid->GetDataPtr("V"));
        RDouble *w = reinterpret_cast<RDouble *>(grid->GetDataPtr("W"));
        RDouble *t = reinterpret_cast<RDouble *>(grid->GetDataPtr("T"));
        RDouble *h = reinterpret_cast<RDouble *>(grid->GetDataPtr("enthalpy"));
        RDouble *cp = reinterpret_cast<RDouble *>(grid->GetDataPtr("cp"));

        int nBoundFace = grid->GetNBoundFace();
        int nTotalCell = grid->GetNTotalCell();

        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            RDouble usq = (u[iCell] * u[iCell] + v[iCell] * v[iCell] + w[iCell] * w[iCell]) / 2;
            t[iCell] += urfT * ((h[iCell] - usq) / cp[iCell] - t[iCell]);

            if (t[iCell] < tMin)
            {
                t[iCell] = tMin;
                std::cout << "temperature is lower than the minimum !" << "ic = " << iCell << std::endl;
            }
            else if (t[iCell] > tMax)
            {
                t[iCell] = tMax;
                std::cout << "temperature is higher than the maximum !" << "ic = " << iCell << std::endl;
            }
        }

        UpdateBoundary();

        for (int iBFace = 0; iBFace < nBoundFace; ++ iBFace)
        {
            int iCell = nTotalCell + iBFace;
            RDouble usq = (u[iCell] * u[iCell] + v[iCell] * v[iCell] + w[iCell] * w[iCell]) / 2;
            h[iCell] = cpg[iCell] * t[iCell] + usq;
        }
    }

    void PBSolverEnergy::SetOutputs(vector<std::string> &phiNameList)
    {
        phiNameList.push_back("T");
        phiNameList.push_back("enthalpy");
    }

    bool PBSolverEnergy::SolutionIsConverged()
    {
        grid->GetData(ScalarName() + "ResNow", &resNow, PHDOUBLE, 1);
        grid->GetData(ScalarName() + "Res", &res, PHDOUBLE, 1);

        if ((resNow < res))
        {
            return false;
        }

        return true;
    }

    void PBSolverEnergy::ComputePostVisualVariables(Post_Visual *postVisualization)
    {
        UnstructGrid *grid = UnstructGridCast(postVisualization->GetGrid());

        RDouble *enthalpy = reinterpret_cast<RDouble*> (grid->GetDataPtr("enthalpy"));
        RDouble *T = reinterpret_cast<RDouble*> (grid->GetDataPtr("T"));

        postVisualization->UpdateVisualNodeVarPtr("enthalpy", enthalpy);
        postVisualization->UpdateVisualNodeVarPtr("T", T);
    }

    void PBSolverEnergy::OutputResidual(ofstream &file, int flag)
    {
        file << setiosflags(ios::left);
        file << setprecision(5);
        file << setiosflags(ios::scientific);

        if (flag == 1)
        {
            file << " Res H" << '\t';
        }
        else
        {
            RDouble hResNow = 0.0;
            RDouble hRes = 0.0;

            grid->GetData(ScalarName() + "ResNow", &hRes, PHDOUBLE, 1);

            hResNow = hRes;

            file << hResNow << '\t';
        }
    }

    DEFINE_PBSOLVER_REGISTER_FUNCTION(PBSolver, PBSolverEnergy, SOLVER_ENERGY);
}