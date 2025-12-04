#include <iostream>
#include "IncomKETurbEpsilonEqCalculator.h"




using namespace std;

namespace PHSPACE
{
IncomKETurbEpsilonEqCalculator::IncomKETurbEpsilonEqCalculator():IncomScalarEqCalculator()
{

}

IncomKETurbEpsilonEqCalculator::~IncomKETurbEpsilonEqCalculator()
{

}

void IncomKETurbEpsilonEqCalculator::SetDiffusionCoeff(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nTotalFace = grid->GetNTotalFace();
    int nBoundFace = grid->GetNBoundFace();
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    RDouble *visl = reinterpret_cast<RDouble *>(grid->GetDataPtr("visl"));
    RDouble *vist = reinterpret_cast<RDouble *>(grid->GetDataPtr("vist"));
    RDouble prandtl = 1.3;
    RDouble **DiffusionCoeff = reinterpret_cast <RDouble **> (GlobalDataBase::GetDataPtr("DiffusionCoeff"));
    int solverIndex = GetSolverIndex();

    for (int iCell = 0; iCell < nTotalCell; ++iCell)
    {
        DiffusionCoeff[solverIndex][iCell] = visl[iCell] + vist[iCell] / prandtl;
    }

    for (int iFace = 0; iFace < nBoundFace; ++iFace)
    {
        int le = leftCellOfFace[iFace];
        int re = rightCellOfFace[iFace];
        DiffusionCoeff[solverIndex][re] = DiffusionCoeff[solverIndex][le];
    }

    CommunicateAnInterfaceVar(DiffusionCoeff[solverIndex]);
}

void IncomKETurbEpsilonEqCalculator::InitFlowAsRestart(Grid *gridIn)
{
    IncomScalarEqCalculator::InitFlowAsRestart(gridIn);
    UnstructGrid *grid = UnstructGridCast(gridIn); 
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    RDouble initEpsilon = GlobalDataBase::GetDoubleParaFromDB("initEpsilon");

    string *varNameIncom = reinterpret_cast <string *> (GlobalDataBase::GetDataPtr("varNameIncom"));
    int solverIndex = GetSolverIndex();
    RDouble *epsilon = reinterpret_cast<RDouble *>(grid->GetDataPtr(varNameIncom[solverIndex]));
    PHSPACE::SetField(epsilon, initEpsilon, nTotal);
}

void IncomKETurbEpsilonEqCalculator::solveScalarEquation(Grid *gridIn, int iEquation)
{
    UnstructGrid *grid = UnstructGridCast(gridIn); 
    int solverIndex = GetSolverIndex();

    WallFunctionKE(grid);

    CalcGrad(grid);
    SetDiffusionCoeff(grid);

    constructMatrixACoeff(grid, solverIndex);
    constructBCoeff(grid, solverIndex);

    WallModification(grid);

    RDouble *diagMatrixCoeff = reinterpret_cast<RDouble *>(grid->GetDataPtr("diagMatrixCoeff"));
    RDouble *upperMatrixCoeff = reinterpret_cast<RDouble *>(grid->GetDataPtr("upperMatrixCoeff"));
    RDouble *lowerMatrixCoeff = reinterpret_cast<RDouble *>(grid->GetDataPtr("lowerMatrixCoeff"));
    RDouble *bCoeff = reinterpret_cast<RDouble *>(grid->GetDataPtr("bCoeff"));

    string *varNameIncom = reinterpret_cast <string *>(GlobalDataBase::GetDataPtr("varNameIncom"));
    string varName = varNameIncom[solverIndex];
    RDouble *q = reinterpret_cast<RDouble *>(grid->GetDataPtr(varName));
    calculateLinearEquation(grid, solverIndex, q, diagMatrixCoeff, bCoeff, upperMatrixCoeff , lowerMatrixCoeff);

    UpdateBCValue(grid);
    UpdateAfterIterloop(grid);
    UpdateProperties(grid);
}

void IncomKETurbEpsilonEqCalculator::calWallBC(Grid *gridIn, vector<int> *faceIndex, Data_Param *bcData, RDouble *diagMatrixCoeff, RDouble *bCoeff, RDouble *q,
    RDouble *dqdx, RDouble *dqdy, RDouble *dqdz, RDouble *mu, RDouble matrixCoeff, int iVariable)
{
    
}

void IncomKETurbEpsilonEqCalculator::WallModification(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    int nTotalCells = grid->GetNTotalCell();
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFaces = grid->GetNBoundFace() - grid->GetNIFace();
    int *leftCellOfFace  = grid->GetLeftCellOfFace();
    int *faceNumberOfEachCell = grid->GetFaceNumberOfEachCell();
    int **cell2face  = grid->GetCell2Face();

    RDouble *Epsilon = reinterpret_cast<RDouble *>(grid->GetDataPtr("Epsilon"));
    RDouble *upperMatrixCoeff = reinterpret_cast<RDouble *>(grid->GetDataPtr("upperMatrixCoeff"));
    RDouble *lowerMatrixCoeff = reinterpret_cast<RDouble *>(grid->GetDataPtr("lowerMatrixCoeff"));
    RDouble *diagMatrixCoeff = reinterpret_cast<RDouble *>(grid->GetDataPtr("diagMatrixCoeff"));
    RDouble *bCoeff = reinterpret_cast<RDouble *>(grid->GetDataPtr("bCoeff"));

    RDouble *soluRelaxCoeff = reinterpret_cast<RDouble *>(GlobalDataBase::GetDataPtr("soluRelaxCoeff"));
    RDouble RelaxCoeff = soluRelaxCoeff[IDX::S_IEPSILON];

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegionUnstruct = unstructBCSet->GetnBCRegion();
    for (int iBCRegionUnstruct = 0; iBCRegionUnstruct < nBCRegionUnstruct; iBCRegionUnstruct++)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegionUnstruct);
        int bcType = bcRegion->GetBCType();
        vector<int> *faceIndex = bcRegion->GetFaceIndex();

        if (bcType == PHENGLEI::SOLID_SURFACE)
        {
            for (int i = 0; i < faceIndex->size(); ++ i)
            {
                int bface = (*faceIndex)[i];
                int le = leftCellOfFace[bface];

                diagMatrixCoeff[le] = 1.0 / RelaxCoeff;
                bCoeff[le] = Epsilon[le] + (1.0-RelaxCoeff) * diagMatrixCoeff[le] * Epsilon[le];

                for (int j = 0; j < grid->GetFaceNumberOfEachCell()[le]; ++j)
                {
                    int iFace = grid->GetCell2Face()[le][j];

                    if (iFace < nBoundFaces)
                    {
                        continue;
                    }
                    else
                    {
                        if (grid->GetLeftCellOfFace()[iFace] == le)
                        {
                            upperMatrixCoeff[iFace] = 0.0;
                        }
                        else
                        {
                            lowerMatrixCoeff[iFace] = 0.0;
                        }
                    }
                }
            }
        }
    }
    
}

void IncomKETurbEpsilonEqCalculator::CalcOtherMatrixACoeff(Grid * gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    RDouble *k = reinterpret_cast<RDouble *>(grid->GetDataPtr("Kinetic"));
    RDouble *rho = reinterpret_cast<RDouble *>(grid->GetDataPtr("rho"));
    RDouble *Epsilon = reinterpret_cast<RDouble *>(grid->GetDataPtr("Epsilon"));
    RDouble *vol = reinterpret_cast<RDouble *>(grid->GetCellVolume());
    RDouble *diagMatrixCoeff = reinterpret_cast<RDouble *>(grid->GetDataPtr("diagMatrixCoeff"));
    RDouble c1m_ = 1.44;
    RDouble c2m_ = 1.92;
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        diagMatrixCoeff[iCell] += c2m_ * rho[iCell] * vol[iCell] * Epsilon[iCell] / (k[iCell] + 1e-20);
    }
}

void IncomKETurbEpsilonEqCalculator::CalcOtherbCoeff(Grid *gridIn, int iEquation)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    RDouble *k = reinterpret_cast<RDouble *>(grid->GetDataPtr("Kinetic"));
    RDouble *gen = reinterpret_cast<RDouble *>(grid->GetDataPtr("gen"));
    RDouble *rho = reinterpret_cast<RDouble *>(grid->GetDataPtr("rho"));
    RDouble *Epsilon = reinterpret_cast<RDouble *>(grid->GetDataPtr("Epsilon"));
    RDouble *vol = reinterpret_cast<RDouble *>(grid->GetCellVolume());
    RDouble *diagMatrixCoeff = reinterpret_cast<RDouble *>(grid->GetDataPtr("diagMatrixCoeff"));
    RDouble *bCoeff = reinterpret_cast<RDouble *>(grid->GetDataPtr("bCoeff"));
    RDouble c1m_ = 1.44;
    RDouble c2m_ = 1.92;
    int *leftCellOfFace = grid->GetLeftCellOfFace();

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        RDouble ke = vol[iCell] * Epsilon[iCell] / (k[iCell] + 1e-20);
        bCoeff[iCell] += c1m_ * gen[iCell] * ke;
    }
}

void IncomKETurbEpsilonEqCalculator::GetResidual(Grid *gridIn, vector<RDouble>& res)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    RDouble resNow = 0.0;

    grid->GetData("EpsilonResNow", &resNow, PHDOUBLE, 1);

    res.push_back(resNow);
}


void IncomKETurbEpsilonEqCalculator::InitialUnsteadyVar(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    PHString1D phiNameList;
    phiNameList.push_back("Epsilon");

    int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
    if (isUnsteady == 0)
    {
        return;
    }

    string TranCalcMethod = GlobalDataBase::GetStrParaFromDB("TranCalcMethod");
    if (TranCalcMethod == "IMPLICIT_EULER")
    {
        ImplicitEuler_ReInitTimeVar(grid, phiNameList);
    }
    else if (TranCalcMethod == "IMPLICIT_2ND_ORDER")
    {
        Implicit2ndOrder_ReInitTimeVar(grid, phiNameList);
    }

    PHString1D().swap(phiNameList);
}

void IncomKETurbEpsilonEqCalculator::UpdateUnsteadyVariable(Grid *grid)
{
    UnstructGrid *gridIn = UnstructGridCast(grid); 
    std::vector<std::string> phiNameList;
    phiNameList.push_back("Epsilon");

    int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
    if (isUnsteady == 0)
    {
        return;
    }

    string TranCalcMethod = GlobalDataBase::GetStrParaFromDB("TranCalcMethod");
    if (TranCalcMethod == "IMPLICIT_EULER")
    {
        ImplicitEuler_SaveOldTimeValue(gridIn, phiNameList);
    }
    else if (TranCalcMethod == "IMPLICIT_2ND_ORDER")
    {
        Implicit2ndOrder_SaveOldTimeValue(gridIn, phiNameList);
    }

}

void IncomKETurbEpsilonEqCalculator::UpdateBCValue(Grid *gridIn)
{
}

void IncomKETurbEpsilonEqCalculator::UpdateProperties(Grid *gridIn)
{

}

void IncomKETurbEpsilonEqCalculator::WallFunctionKE(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int *leftCellOfFace = grid->GetLeftCellOfFace();

    RDouble *u = reinterpret_cast<RDouble *>(grid->GetDataPtr("U"));
    RDouble *v = reinterpret_cast<RDouble *>(grid->GetDataPtr("V"));
    RDouble *w = reinterpret_cast<RDouble *>(grid->GetDataPtr("W"));
    RDouble *dudx = reinterpret_cast<RDouble *>(grid->GetDataPtr("dUdx"));
    RDouble *dudy = reinterpret_cast<RDouble *>(grid->GetDataPtr("dUdy"));
    RDouble *dudz = reinterpret_cast<RDouble *>(grid->GetDataPtr("dUdz"));
    RDouble *dvdx = reinterpret_cast<RDouble *>(grid->GetDataPtr("dVdx"));
    RDouble *dvdy = reinterpret_cast<RDouble *>(grid->GetDataPtr("dVdy"));
    RDouble *dvdz = reinterpret_cast<RDouble *>(grid->GetDataPtr("dVdz"));
    RDouble *dwdx = reinterpret_cast<RDouble *>(grid->GetDataPtr("dWdx"));
    RDouble *dwdy = reinterpret_cast<RDouble *>(grid->GetDataPtr("dWdy"));
    RDouble *dwdz = reinterpret_cast<RDouble *>(grid->GetDataPtr("dWdz"));
    RDouble *vist = reinterpret_cast<RDouble *>(grid->GetDataPtr("vist"));
    RDouble *visl = reinterpret_cast<RDouble *>(grid->GetDataPtr("visl"));
    RDouble *gen = reinterpret_cast<RDouble *>(grid->GetDataPtr("gen"));
    RDouble *Epsilon = reinterpret_cast<RDouble *>(grid->GetDataPtr("Epsilon"));
    RDouble *rho = reinterpret_cast<RDouble *>(grid->GetDataPtr("rho"));
    
    string GradCalcMethod = GlobalDataBase::GetStrParaFromDB("GradCalcMethod");
    
    GradientCalculation(grid, "U", "dUdx", "dUdy", "dUdz", GradCalcMethod);
    GradientCalculation(grid, "V", "dVdx", "dVdy", "dVdz", GradCalcMethod);
    GradientCalculation(grid, "W", "dWdx", "dWdy", "dWdz", GradCalcMethod);

    for (int iCell = 0; iCell < nTotalCell; ++iCell)
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
    }

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegionUnstruct = unstructBCSet->GetnBCRegion();

    for (int iBCRegionUnstruct = 0; iBCRegionUnstruct < nBCRegionUnstruct; iBCRegionUnstruct++)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegionUnstruct);
        int bcType = bcRegion->GetBCType();
        if (PHENGLEI::SOLID_SURFACE == bcType)
        {
            vector<int> *faceIndex = bcRegion->GetFaceIndex();
            for (int iFace = 0; iFace < faceIndex->size(); ++iFace)
            {
                int iFacelocal = (*faceIndex)[iFace];
                int le = leftCellOfFace[iFacelocal];

                gen[le] = 0.0;
                Epsilon[le] = 0.0;
            }
        }
    }

    for (int iBCRegionUnstruct = 0; iBCRegionUnstruct < nBCRegionUnstruct; iBCRegionUnstruct++)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegionUnstruct);
        int bcType = bcRegion->GetBCType();
        vector<int> *faceIndex = bcRegion->GetFaceIndex();
        Data_Param *bcData = bcRegion->GetBCParamDataBase();

        if (bcType == PHENGLEI::SOLID_SURFACE)
        {
            int nTotalCell = grid->GetNTotalCell();
            int *leftCellOfFace = grid->GetLeftCellOfFace();
            int *rightCellOfFace = grid->GetRightCellOfFace();
            int nBoundFace = grid->GetNBoundFace();
            RDouble *xfn = grid->GetFaceNormalX();
            RDouble *yfn = grid->GetFaceNormalY();
            RDouble *zfn = grid->GetFaceNormalZ();
            RDouble *WD = reinterpret_cast<RDouble *>(grid->GetDataPtr("wd"));
            RDouble *gen = reinterpret_cast<RDouble *>(grid->GetDataPtr("gen"));
            RDouble *k = reinterpret_cast<RDouble *>(grid->GetDataPtr("Kinetic"));
            RDouble *rho = reinterpret_cast<RDouble *>(grid->GetDataPtr("rho"));
            RDouble *visl = reinterpret_cast<RDouble *>(grid->GetDataPtr("visl"));
            RDouble *vist = reinterpret_cast<RDouble *>(grid->GetDataPtr("vist"));
            RDouble *u = reinterpret_cast<RDouble *>(grid->GetDataPtr("U"));
            RDouble *v = reinterpret_cast<RDouble *>(grid->GetDataPtr("V"));
            RDouble *w = reinterpret_cast<RDouble *>(grid->GetDataPtr("W"));
            RDouble *mu = reinterpret_cast<RDouble *>(grid->GetDataPtr("mu"));
            RDouble *Epsilon = reinterpret_cast<RDouble *>(grid->GetDataPtr("Epsilon"));
            RDouble *wallFlag = reinterpret_cast<RDouble *>(grid->GetDataPtr("wallFlag"));

            RDouble cmu = 0.09;
            RDouble cmu75_ = std::pow(cmu, 0.75);
            RDouble cmu25_ = std::pow(cmu, 0.25);
            RDouble cappa_ = 0.41;
            RDouble yptr_ = 11.06;

            for (int iFace = 0; iFace < faceIndex->size(); ++iFace)
            {
                int iFacelocal = (*faceIndex)[iFace];
                int le = leftCellOfFace[iFacelocal];
                int re = rightCellOfFace[iFacelocal];

                RDouble wd = WD[re];
                RDouble cmute = cmu25_ * sqrt(k[le]);
                RDouble yplus = rho[le] * wd * cmute / visl[le];
                
                RDouble du = u[le];
                RDouble dv = v[le];
                RDouble dw = w[le];
                RDouble ax = xfn[iFacelocal];
                RDouble ay = yfn[iFacelocal];
                RDouble az = zfn[iFacelocal];
                RDouble dun = (du * ax + dv * ay + dw * az);
                RDouble Ux = du - dun * ax;
                RDouble Uy = dv - dun * ay;
                RDouble Uz = dw - dun * az;
                RDouble Uwall = std::sqrt(Ux * Ux + Uy * Uy + Uz * Uz) / wd;

                if(yplus > yptr_)
                {
                    Epsilon[le] += (cmu75_ * std::pow(k[le], 1.5) / (cappa_ * wd)) / wallFlag[le];
                    gen[le] += (visl[le] + vist[re]) * Uwall * cmute / (cappa_ * wd) / wallFlag[le];
                }
                else
                {
                    Epsilon[le] += (cmu * rho[le] * std::pow(k[le], 2.0) / visl[le]) / wallFlag[le];
                    gen[le] += visl[le] * std::pow(Uwall, 2.0) / wallFlag[le];
                }
            }
        }
    }
}

void IncomKETurbEpsilonEqCalculator::CalcGenTerm(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();

    RDouble *u = reinterpret_cast<RDouble *>(grid->GetDataPtr("U"));
    RDouble *v = reinterpret_cast<RDouble *>(grid->GetDataPtr("V"));
    RDouble *w = reinterpret_cast<RDouble *>(grid->GetDataPtr("W"));
    RDouble *dudx = reinterpret_cast<RDouble *>(grid->GetDataPtr("dUdx"));
    RDouble *dudy = reinterpret_cast<RDouble *>(grid->GetDataPtr("dUdy"));
    RDouble *dudz = reinterpret_cast<RDouble *>(grid->GetDataPtr("dUdz"));
    RDouble *dvdx = reinterpret_cast<RDouble *>(grid->GetDataPtr("dVdx"));
    RDouble *dvdy = reinterpret_cast<RDouble *>(grid->GetDataPtr("dVdy"));
    RDouble *dvdz = reinterpret_cast<RDouble *>(grid->GetDataPtr("dVdz"));
    RDouble *dwdx = reinterpret_cast<RDouble *>(grid->GetDataPtr("dWdx"));
    RDouble *dwdy = reinterpret_cast<RDouble *>(grid->GetDataPtr("dWdy"));
    RDouble *dwdz = reinterpret_cast<RDouble *>(grid->GetDataPtr("dWdz"));
    RDouble *vist = reinterpret_cast<RDouble *>(grid->GetDataPtr("vist"));
    RDouble *visl = reinterpret_cast<RDouble *>(grid->GetDataPtr("visl"));
    RDouble *gen = reinterpret_cast<RDouble *>(grid->GetDataPtr("gen"));
    RDouble *yplus = reinterpret_cast<RDouble *>(grid->GetDataPtr("yplus"));
    RDouble *Epsilon = reinterpret_cast<RDouble *>(grid->GetDataPtr("Epsilon"));
    RDouble *rho = reinterpret_cast<RDouble *>(grid->GetDataPtr("rho"));
    
    string GradCalcMethod = GlobalDataBase::GetStrParaFromDB("GradCalcMethod");

    GradientCalculation(grid, "U", "dUdx", "dUdy", "dUdz", GradCalcMethod);
    GradientCalculation(grid, "V", "dVdx", "dVdy", "dVdz", GradCalcMethod);
    GradientCalculation(grid, "W", "dWdx", "dWdy", "dWdz", GradCalcMethod);

    for (int iCell = 0; iCell < nTotalCell; ++iCell)
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

        gen[iCell] = min(gen[iCell], 10.0 * rho[iCell] * Epsilon[iCell]);
    }

}

void IncomKETurbEpsilonEqCalculator::UpdateWallGen(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegionUnstruct = unstructBCSet->GetnBCRegion();

    for (int iBCRegionUnstruct = 0; iBCRegionUnstruct < nBCRegionUnstruct; iBCRegionUnstruct++)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegionUnstruct);
        int bcType = bcRegion->GetBCType();
        vector<int> *faceIndex = bcRegion->GetFaceIndex();
        Data_Param* bcData = bcRegion->GetBCParamDataBase();

        if (bcType == PHENGLEI::SOLID_SURFACE)

        {
            int nTotalCell = grid->GetNTotalCell();
            int *leftCellOfFace = grid->GetLeftCellOfFace();
            int *rightCellOfFace = grid->GetRightCellOfFace();
            int nBoundFace = grid->GetNBoundFace();
            RDouble *xfn = grid->GetFaceNormalX();
            RDouble *yfn = grid->GetFaceNormalY();
            RDouble *zfn = grid->GetFaceNormalZ();
            RDouble *WD = reinterpret_cast<RDouble *>(grid->GetDataPtr("wd"));
            RDouble *gen = reinterpret_cast<RDouble *>(grid->GetDataPtr("gen"));
            RDouble *k = reinterpret_cast<RDouble *>(grid->GetDataPtr("Kinetic"));
            RDouble *rho = reinterpret_cast<RDouble *>(grid->GetDataPtr("rho"));
            RDouble *visl = reinterpret_cast<RDouble *>(grid->GetDataPtr("visl"));
            RDouble *vist = reinterpret_cast<RDouble *>(grid->GetDataPtr("vist"));
            RDouble *vistb = vist + nTotalCell;
            RDouble *u = reinterpret_cast<RDouble *>(grid->GetDataPtr("U"));
            RDouble *v = reinterpret_cast<RDouble *>(grid->GetDataPtr("V"));
            RDouble *w = reinterpret_cast<RDouble *>(grid->GetDataPtr("W"));
            RDouble *mu = reinterpret_cast<RDouble *>(grid->GetDataPtr("mu"));
            RDouble *yplus = reinterpret_cast<RDouble *>(grid->GetDataPtr("yplus")) + nTotalCell;
            RDouble *Epsilon = reinterpret_cast<RDouble *>(grid->GetDataPtr("Epsilon"));
            RDouble *nu_eff_w = reinterpret_cast<RDouble *>(grid->GetDataPtr("nu_eff")) + nTotalCell;
            RDouble *wallFlag = reinterpret_cast<RDouble *>(grid->GetDataPtr("wallFlag"));

            RDouble cmu75_ = 0.1643;
            RDouble cappa_ = 0.4187;
            RDouble cmu = 0.09;
            RDouble cmu25_ = 0.5477;
            RDouble yptr_ = 11.225;
            RDouble yvstar_ = 20.0;
            RDouble cl_ = 2.55;
            RDouble elog_ = 9.793;

            for (int iFace = 0; iFace < faceIndex->size(); ++iFace)
            {
                int iFacelocal = (*faceIndex)[iFace];
                int le = leftCellOfFace[iFacelocal];
                int re = rightCellOfFace[iFacelocal];

                RDouble wd = WD[re];
                RDouble cmute = cmu25_ * sqrt(k[le]);
                yplus[iFacelocal] = rho[le] * wd * cmute / (visl[le]);

                RDouble du = u[re] - u[le];
                RDouble dv = v[re] - v[le];
                RDouble dw = w[re] - w[le];
                RDouble ax = xfn[iFacelocal];
                RDouble ay = yfn[iFacelocal];
                RDouble az = zfn[iFacelocal];
                RDouble dun = (du * ax + dv * ay + dw * az);
                RDouble Ux = du - dun * ax;
                RDouble Uy = dv - dun * ay;
                RDouble Uz = dw - dun * az;
                RDouble Uwall = std::sqrt(Ux * Ux + Uy * Uy + Uz * Uz);

                yplus[iFacelocal] = max(yplus[iFacelocal], yptr_);
                double tau_w = Uwall * cmute * rho[le] / (std::log(elog_ * yplus[iFacelocal]) / cappa_);
                double d = yplus[iFacelocal] * visl[le] / (rho[le] * cmute);
                nu_eff_w[iFacelocal] = tau_w * wd / Uwall;
                gen[le] += (tau_w * tau_w / (cappa_ * rho[le] * cmute * d)) / wallFlag[le];
                Epsilon[le] += (cmu75_ * std::pow(k[le], 1.5) / (cappa_ * d)) / wallFlag[le];
            }
        }
    }
}

}