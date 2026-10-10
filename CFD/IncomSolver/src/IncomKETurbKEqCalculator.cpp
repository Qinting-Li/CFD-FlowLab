#include <iostream>
#include <memory>
#include <cmath>
#include "IncomKETurbKEqCalculator.h"

using namespace std;

namespace PHSPACE
{
IncomKETurbKEqCalculator::IncomKETurbKEqCalculator():IncomScalarEqCalculator()
{
}

IncomKETurbKEqCalculator::~IncomKETurbKEqCalculator()
{

}

void IncomKETurbKEqCalculator::InitFlowAsRestart(Grid *gridIn)
{
    IncomScalarEqCalculator::InitFlowAsRestart(gridIn);
    UnstructGrid *grid = UnstructGridCast(gridIn); 
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    RDouble initKinetic = GlobalDataBase::GetDoubleParaFromDB("initKinetic");

    string *varNameIncom = reinterpret_cast <string *> (GlobalDataBase::GetDataPtr("varNameIncom"));
    int solverIndex = GetSolverIndex();
    RDouble *kinetic = reinterpret_cast<RDouble *>(grid->GetDataPtr(varNameIncom[solverIndex]));
    PHSPACE::SetField(kinetic, initKinetic, nTotal);

    RDouble *visl = reinterpret_cast<RDouble *> (grid->GetDataPtr("visl"));
    PHSPACE::SetField(visl, 0.0, nTotal);

    RDouble *vist = reinterpret_cast<RDouble *> (grid->GetDataPtr("vist"));
    PHSPACE::SetField(vist, 0.0, nTotal);

    RDouble *nu_eff = reinterpret_cast<RDouble *> (grid->GetDataPtr("nu_eff"));
    PHSPACE::SetField(nu_eff, 0.0, nTotal);

    RDouble *yplus = reinterpret_cast<RDouble *>(grid->GetDataPtr("yplus"));
    PHSPACE::SetField(yplus, 0.0, nTotal);

    RDouble *gen = reinterpret_cast<RDouble *> (grid->GetDataPtr("gen"));
    PHSPACE::SetField(gen, 0.0, nTotalCell);

    RDouble *wallFlag = reinterpret_cast<RDouble *> (grid->GetDataPtr("wallFlag"));
    PHSPACE::SetField(wallFlag, 0.0, nTotalCell);
}

void IncomKETurbKEqCalculator::AllocateGlobalVar(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn); 
    IncomScalarEqCalculator::AllocateGlobalVar(grid);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    RDouble *visl = NewPointer<RDouble>(nTotal);
    grid->UpdateDataPtr("visl", visl);
    
    RDouble *vist = NewPointer<RDouble>(nTotal);
    grid->UpdateDataPtr("vist", vist);

    RDouble *nu_eff = NewPointer<RDouble>(nTotal);
    grid->UpdateDataPtr("nu_eff", nu_eff);

    RDouble *yplus = NewPointer<RDouble>(nTotal);
    grid->UpdateDataPtr("yplus", yplus);

    RDouble *gen = NewPointer<RDouble>(nTotalCell);
    grid->UpdateDataPtr("gen", gen);

    RDouble *wallFlag = NewPointer<RDouble>(nTotalCell);
    grid->UpdateDataPtr("wallFlag", wallFlag);
}

void IncomKETurbKEqCalculator::IncompressibleInitial(Grid * gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn); 
    int nTotalCell = grid->GetNTotalCell();
    int nTotalFace = grid->GetNTotalFace();
    int nBoundFace = grid->GetNBoundFace();
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    RDouble *k = reinterpret_cast<RDouble *>(grid->GetDataPtr("Kinetic"));
    RDouble *Epsilon = reinterpret_cast<RDouble *>(grid->GetDataPtr("Epsilon"));
    RDouble *visl = reinterpret_cast<RDouble *>(grid->GetDataPtr("visl"));
    RDouble *vist = reinterpret_cast<RDouble *>(grid->GetDataPtr("vist"));
    RDouble *wallFlag = reinterpret_cast<RDouble *>(grid->GetDataPtr("wallFlag"));
    RDouble *nu_eff = reinterpret_cast<RDouble *>(grid->GetDataPtr("nu_eff"));
    RDouble *nu_eff_w = nu_eff + nTotalCell;
    RDouble cmu = 0.09; 
    RDouble *rho = reinterpret_cast<RDouble *>(grid->GetDataPtr("rho"));
    RDouble *mu = reinterpret_cast<RDouble *>(grid->GetDataPtr("mu"));
    RDouble *FaceFlux = reinterpret_cast<RDouble*>(grid->GetDataPtr("FaceFlux"));

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        visl[iCell] = mu[iCell];
        vist[iCell] = rho[iCell] * cmu * k[iCell] * k[iCell] / Epsilon[iCell];
    }

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int le = leftCellOfFace[iFace];
        int re = rightCellOfFace[iFace];
        vist[re] = vist[le];
        mu[re] = mu[le];
        nu_eff_w[iFace] = mu[le];
    }

    for (int iCell = 0; iCell < nTotalCell; iCell++)
    {
        wallFlag[iCell] = 0.0;
    }

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegionUnstruct = unstructBCSet->GetnBCRegion();

    for (int iBCRegionUnstruct = 0; iBCRegionUnstruct < nBCRegionUnstruct; iBCRegionUnstruct++)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegionUnstruct);
        int bcType = bcRegion->GetBCType();
        if (2 == bcType)
        {
            vector<int> *faceIndex = bcRegion->GetFaceIndex();

            for (int iFace = 0; iFace < faceIndex->size(); ++iFace)
            {
                int iFacelocal = (*faceIndex)[iFace];
                int leftCell = leftCellOfFace[iFacelocal];
                
                wallFlag[leftCell] += 1.0;
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
            for (size_t iFace = 0; iFace < faceIndex->size(); ++ iFace)
            {
                int iFacelocal = (*faceIndex)[iFace];
                int le = leftCellOfFace[iFacelocal];
                int re = rightCellOfFace[iFacelocal];
                k[re] = k[le];
                Epsilon[re] = Epsilon[le];
            }
        }
        else if (bcType == PHENGLEI::SYMMETRY)
        {
            for (size_t iFace = 0; iFace < faceIndex->size(); ++ iFace)
            {
                int iFacelocal = (*faceIndex)[iFace];
                int le = leftCellOfFace[iFacelocal];
                int re = rightCellOfFace[iFacelocal];
                k[re] = k[le];
                Epsilon[re] = Epsilon[le];
            }
        }
        else if (bcType == PHENGLEI::FARFIELD)
        {
            RDouble kb = 0.0;
            RDouble epsilonb = 0.0;
            for (std::size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
            {
                int iFacelocal = (*faceIndex)[iFace];

                if (bcData)
                {
                    if (bcData->IsExist("initKinetic", PHDOUBLE, 1))
                    {
                        bcData->GetData("initKinetic", &kb, PHDOUBLE, 1);
                    }
                }
                else
                {
                    std::cout << " No initial value has assigned to the boundary " << std::endl;
                }

                if (bcData)
                {
                    if (bcData->IsExist("initEpsilon", PHDOUBLE, 1))
                    {
                        bcData->GetData("initEpsilon", &epsilonb, PHDOUBLE, 1);
                    }
                }
                else
                {
                    std::cout << "Please assign a value of Epsilon on farield" << std::endl;;
                }
            }
            for (int iFace = 0; iFace < faceIndex->size(); ++iFace)
            {
                int iFacelocal = (*faceIndex)[iFace];
                int le = leftCellOfFace[iFacelocal];
                int re = rightCellOfFace[iFacelocal];
                if (FaceFlux[iFacelocal] < 0)
                {

                    k[re] = kb;
                    Epsilon[re] = epsilonb;
                }
                else
                {

                    k[re] = k[le];
                    Epsilon[re] = Epsilon[le];
                }
            }
        }
        else if (bcType == PHENGLEI::INFLOW)
        {
            RDouble initKinetic = 0.0;
            RDouble initEpsilon = 0.0;

            for (std::size_t iFace = 0; iFace < faceIndex->size(); ++ iFace)
            {
                int iFacelocal = (*faceIndex)[iFace];
                int le = leftCellOfFace[iFacelocal];
                int re = rightCellOfFace[iFacelocal];

                if (bcData)
                {
                    if (bcData->IsExist("initKinetic", PHDOUBLE, 1))
                    {
                        bcData->GetData("initKinetic", &initKinetic, PHDOUBLE, 1);
                    }
                }
                else
                {
                    std::cout << " No initKinetic initial value has assigned to the inlet boundary " << std::endl;
                }

                if (bcData)
                {
                    if (bcData->IsExist("initEpsilon", PHDOUBLE, 1))
                    {
                        bcData->GetData("initEpsilon", &initEpsilon, PHDOUBLE, 1);
                    }
                }
                else
                {
                    std::cout << "Please assign a value of Epsilon on inlet" << std::endl;;
                }

                if (FaceFlux[iFacelocal] < 0.0)
                {

                    k[re] = initKinetic;
                    Epsilon[re] = initEpsilon;
                }
                else
                {

                    k[re] = k[le];
                    Epsilon[re] = Epsilon[le];
                }
            }
        }
        else if (bcType == PHENGLEI::OUTFLOW)
        {
            RDouble initKinetic = 0.0;
            RDouble initEpsilon = 0.0;

            for (std::size_t iFace = 0; iFace < faceIndex->size(); ++ iFace)
            {
                int iFacelocal = (*faceIndex)[iFace];
                int le = leftCellOfFace[iFacelocal];
                int re = rightCellOfFace[iFacelocal];

                if (bcData)
                {
                    if (bcData->IsExist("initKinetic", PHDOUBLE, 1))
                    {
                        bcData->GetData("initKinetic", &initKinetic, PHDOUBLE, 1);
                    }
                }
                else
                {
                    std::cout << " No initKinetic initial value has assigned to the boundary " << std::endl;
                }

                if (bcData)
                {
                    if (bcData->IsExist("initEpsilon", PHDOUBLE, 1))
                    {
                        bcData->GetData("initEpsilon", &initEpsilon, PHDOUBLE, 1);
                    }
                }
                else
                {
                    std::cout << "Please assign a value of Epsilon on inlet" << std::endl;;
                }

                if (FaceFlux[iFacelocal] > 0.0)
                {
                    k[re] = initKinetic;
                    Epsilon[re] = initEpsilon;
                }
                else
                {
                    k[re] = k[le];
                    Epsilon[re] = Epsilon[le];
                }
            }
        }
        else if (bcType == PHENGLEI::PRESSURE_INLET)
        {
            RDouble kb = 0.0;
            for (std::size_t iFace = 0; iFace < faceIndex->size(); ++ iFace)
            {
                int iFacelocal = (*faceIndex)[iFace];
                int le = leftCellOfFace[iFacelocal];
                int re = rightCellOfFace[iFacelocal];
                if (bcData)
                {
                    if (bcData->IsExist("initKinetic", PHDOUBLE, 1))
                    {
                        bcData->GetData("initKinetic", &kb, PHDOUBLE, 1);
                    }
                }
                else
                {
                    std::cout << " No initial value has assigned to the boundary " << std::endl;
                }
                if (FaceFlux[iFacelocal] < 0)
                {
                    k[re] = kb;
                }
                else
                {
                    k[re] = k[le];
                }
            }
        }
        else if (bcType == PHENGLEI::PRESSURE_OUTLET)
        {
            RDouble kb = 0.0;
            for (std::size_t iFace = 0; iFace < faceIndex->size(); ++ iFace)
            {
                int iFacelocal = (*faceIndex)[iFace];
                int le = leftCellOfFace[iFacelocal];
                int re = rightCellOfFace[iFacelocal];

                if (FaceFlux[iFacelocal] < 0.0)
                {
                    if (bcData)
                    {
                        if (bcData->IsExist("initKinetic", PHDOUBLE, 1))
                        {
                            bcData->GetData("initKinetic", &kb, PHDOUBLE, 1);
                            k[re] = RDouble(kb);
                        }
                        else
                        {
                            std::cout << " Please assign a value " << std::endl;
                        }
                    }
                }
                else
                {
                    k[re] = k[le];
                }
            }
        }
        else 
        {
            for (size_t iFace = 0; iFace < faceIndex->size(); ++ iFace)
            {
                int iFacelocal = (*faceIndex)[iFace];
                int le = leftCellOfFace[iFacelocal];
                int re = rightCellOfFace[iFacelocal];
                k[re] = k[le];
                Epsilon[re] = Epsilon[le];
            }
        }
    }
    CommunicateAnInterfaceVar(k);
    CommunicateAnInterfaceVar(Epsilon);

    InitialUnsteadyVar(grid);
}

void IncomKETurbKEqCalculator::SetDiffusionCoeff(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nTotalFace = grid->GetNTotalFace();
    int nBoundFace = grid->GetNBoundFace();
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    RDouble *visl = reinterpret_cast<RDouble *>(grid->GetDataPtr("visl"));
    RDouble *vist = reinterpret_cast<RDouble *>(grid->GetDataPtr("vist"));
    RDouble prandtl = 1.0;
    RDouble **DiffusionCoeff = reinterpret_cast<RDouble **>(GlobalDataBase::GetDataPtr("DiffusionCoeff"));
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

void IncomKETurbKEqCalculator::UpdateProperties(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nTotalFace = grid->GetNTotalFace();
    int nBoundFace = grid->GetNBoundFace();
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    RDouble *k = reinterpret_cast<RDouble *>(grid->GetDataPtr("Kinetic"));
    RDouble *Epsilon = reinterpret_cast<RDouble *>(grid->GetDataPtr("Epsilon"));
    RDouble *visl = reinterpret_cast<RDouble *>(grid->GetDataPtr("visl"));
    RDouble *vist = reinterpret_cast<RDouble *>(grid->GetDataPtr("vist"));
    RDouble cmu = 0.09; 
    RDouble *rho = reinterpret_cast<RDouble *>(grid->GetDataPtr("rho"));
    RDouble *mu  = reinterpret_cast<RDouble *>(grid->GetDataPtr("mu"));

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

            RDouble *mu = reinterpret_cast<RDouble *>(grid->GetDataPtr("mu"));
            RDouble *nu_eff_w = reinterpret_cast<RDouble *>(grid->GetDataPtr("nu_eff")) + nTotalCell;
            RDouble *WD = reinterpret_cast<RDouble *>(grid->GetDataPtr("wd"));

            RDouble cmu = 0.09;
            RDouble cmu75_ = std::pow(cmu, 0.75);
            RDouble cmu25_ = std::pow(cmu, 0.25);
            RDouble cappa_ = 0.41;
            RDouble yptr_ = 11.06;

            for (int iFace = 0; iFace < faceIndex->size(); ++ iFace)
            {
                int iFacelocal = (*faceIndex)[iFace];
                int le = leftCellOfFace[iFacelocal];
                int re = rightCellOfFace[iFacelocal];

                RDouble wd = WD[re];
                RDouble cmute = cmu25_ * sqrt(k[le]);
                RDouble yplus = rho[le] * wd * cmute / visl[le];

                if(yplus > yptr_)
                {
                    vist[re] = (yplus / (log(yplus) / cappa_ + 5.25) - 1.0) * visl[le];
                    if(vist[re] < 0.0)
                    {
                        vist[re] = 0.0;
                    }
                    mu[re] = visl[le] + vist[re];
                }
                else
                {
                    vist[re] = 0.0;
                    mu[re] = visl[le];
                }
            }
        }
        else if(bcType == PHENGLEI::SYMMETRY || bcType == PHENGLEI::INFLOW || bcType == PHENGLEI::OUTFLOW)
        {
            for (int iFace = 0; iFace < faceIndex->size(); ++ iFace)
            {
                int iFacelocal = (*faceIndex)[iFace];
                int le = leftCellOfFace[iFacelocal];
                int re = rightCellOfFace[iFacelocal];
                vist[re] = vist[le];
                mu[re] = visl[le] + vist[re];
            }
        }
    }

    CommunicateAnInterfaceVar(mu);
    CommunicateAnInterfaceVar(vist);
}


void IncomKETurbKEqCalculator::CalcOtherMatrixACoeff(Grid * gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    RDouble *k = reinterpret_cast<RDouble *>(grid->GetDataPtr("Kinetic"));
    RDouble *gen = reinterpret_cast<RDouble *>(grid->GetDataPtr("gen"));
    RDouble *diagMatrixCoeff = reinterpret_cast<RDouble *>(grid->GetDataPtr("diagMatrixCoeff"));
    RDouble *bCoeff = reinterpret_cast<RDouble *>(grid->GetDataPtr("bCoeff"));
    RDouble *rho = reinterpret_cast<RDouble *>(grid->GetDataPtr("rho"));
    RDouble *Epsilon = reinterpret_cast<RDouble *>(grid->GetDataPtr("Epsilon"));
    RDouble *vol = reinterpret_cast<RDouble *>(grid->GetCellVolume());

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        diagMatrixCoeff[iCell] += rho[iCell] * Epsilon[iCell] / (k[iCell] + 1e-20) * vol[iCell];
    }
}

void IncomKETurbKEqCalculator::CalcOtherbCoeff(Grid *gridIn, int iEquation)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    int nTotalCell = grid->GetNTotalCell();
    RDouble *k = reinterpret_cast<RDouble *>(grid->GetDataPtr("Kinetic"));
    RDouble *gen = reinterpret_cast<RDouble *>(grid->GetDataPtr("gen"));
    RDouble *diagMatrixCoeff = reinterpret_cast<RDouble *>(grid->GetDataPtr("diagMatrixCoeff"));
    RDouble *bCoeff = reinterpret_cast<RDouble *>(grid->GetDataPtr("bCoeff"));
    RDouble *rho = reinterpret_cast<RDouble *>(grid->GetDataPtr("rho"));
    RDouble *Epsilon = reinterpret_cast<RDouble *>(grid->GetDataPtr("Epsilon"));
    RDouble *vol = reinterpret_cast<RDouble *>(grid->GetCellVolume());

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        bCoeff[iCell] += gen[iCell] * vol[iCell];
    }
}


void IncomKETurbKEqCalculator::GetResidual(Grid *gridIn, vector<RDouble>& res)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    RDouble resNow = 0.0;

    grid->GetData("KineticResNow", &resNow, PHDOUBLE, 1);

    res.push_back(resNow);
}

void IncomKETurbKEqCalculator::InitialUnsteadyVar(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    PHString1D phiNameList;
    phiNameList.push_back("Kinetic");

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

void IncomKETurbKEqCalculator::UpdateUnsteadyVariable(Grid *grid)
{
    UnstructGrid *gridIn = UnstructGridCast(grid); 
    std::vector<std::string> phiNameList;
    phiNameList.push_back("Kinetic");

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

void IncomKETurbKEqCalculator::UpdateBCValue(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn); 
    int nTotalCell = grid->GetNTotalCell();
    int nTotalFace = grid->GetNTotalFace();
    int nBoundFace = grid->GetNBoundFace();
    int *leftCellOfFace  = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    string varNameK = "Kinetic";
    string varNameE = "Epsilon";

    RDouble *kinetic = reinterpret_cast<RDouble *>(grid->GetDataPtr("Kinetic"));
    RDouble *epsilon = reinterpret_cast<RDouble *>(grid->GetDataPtr("Epsilon"));
    RDouble *visl = reinterpret_cast<RDouble *>(grid->GetDataPtr("visl"));
    RDouble *vist = reinterpret_cast<RDouble *>(grid->GetDataPtr("vist"));
    RDouble *rho = reinterpret_cast<RDouble *>(grid->GetDataPtr("rho"));
    RDouble *mu  = reinterpret_cast<RDouble *>(grid->GetDataPtr("mu"));
    RDouble *FaceFlux = reinterpret_cast<RDouble *>(grid->GetDataPtr("FaceFlux"));
    RDouble *wd = reinterpret_cast<RDouble *>(grid->GetDataPtr("wd"));
    RDouble cmu = 0.09;

    for (int iCell = 0; iCell < nTotalCell; iCell++)
    {
        if (kinetic[iCell] < 1e-14)
        {
            kinetic[iCell] = 1e-14;
        }

        if (epsilon[iCell] < 1e-20)
        {
            epsilon[iCell] = 1e-20;
        }

        vist[iCell] = rho[iCell] * cmu * kinetic[iCell] * kinetic[iCell] / epsilon[iCell];
        
        if(vist[iCell] > 1e5 * visl[iCell])
        {
            vist[iCell] = 1e5 * visl[iCell];
        }
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

        if (bcType == PHENGLEI::SOLID_SURFACE)
        {
            for (size_t iFace = 0; iFace < faceIndex->size(); ++ iFace)
            {
                int iFacelocal = (*faceIndex)[iFace];
                int le = leftCellOfFace[iFacelocal];
                int re = rightCellOfFace[iFacelocal];
                kinetic[re] = 0.0;
                epsilon[re] = 2.0 * visl[le] * kinetic[le] / rho[le] / pow(wd[re],2.0);
            }
        }
        else if (bcType == PHENGLEI::SYMMETRY)
        {
            for (size_t iFace = 0; iFace < faceIndex->size(); ++ iFace)
            {
                int iFacelocal = (*faceIndex)[iFace];
                int le = leftCellOfFace[iFacelocal];
                int re = rightCellOfFace[iFacelocal];
                kinetic[re] = kinetic[le];
                epsilon[re] = epsilon[le];
            }
        }
        else if (bcType == PHENGLEI::FARFIELD)
        {
            RDouble kb = 0.0;
            RDouble epsilonb = 0.0;
            for (std::size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
            {
                int iFacelocal = (*faceIndex)[iFace];

                if (bcData)
                {
                    if (bcData->IsExist("initKinetic", PHDOUBLE, 1))
                    {
                        bcData->GetData("initKinetic", &kb, PHDOUBLE, 1);
                    }
                }
                else
                {
                    std::cout << " No initial value has assigned to the boundary " << std::endl;
                }

                if (bcData)
                {
                    if (bcData->IsExist("initEpsilon", PHDOUBLE, 1))
                    {
                        bcData->GetData("initEpsilon", &epsilonb, PHDOUBLE, 1);
                    }
                }
                else
                {
                    std::cout << "Please assign a value of Epsilon on farield" << std::endl;;
                }
            }
            for (int iFace = 0; iFace < faceIndex->size(); ++iFace)
            {
                int iFacelocal = (*faceIndex)[iFace];
                int le = leftCellOfFace[iFacelocal];
                int re = rightCellOfFace[iFacelocal];
                if (FaceFlux[iFacelocal] < 0)
                {

                    kinetic[re] = kb;
                    epsilon[re] = epsilonb;
                }
                else
                {

                    kinetic[re] = kinetic[le];
                    epsilon[re] = epsilon[le];
                }
            }
        }
        else if (bcType == PHENGLEI::INFLOW)
        {
            RDouble initKinetic = 0.0;
            RDouble initEpsilon = 0.0;

            for (std::size_t iFace = 0; iFace < faceIndex->size(); ++ iFace)
            {
                int iFacelocal = (*faceIndex)[iFace];
                int le = leftCellOfFace[iFacelocal];
                int re = rightCellOfFace[iFacelocal];

                if (bcData)
                {
                    if (bcData->IsExist("initKinetic", PHDOUBLE, 1))
                    {
                        bcData->GetData("initKinetic", &initKinetic, PHDOUBLE, 1);
                    }
                }
                else
                {
                    std::cout << " No initKinetic initial value has assigned to the inlet boundary " << std::endl;
                }

                if (bcData)
                {
                    if (bcData->IsExist("initEpsilon", PHDOUBLE, 1))
                    {
                        bcData->GetData("initEpsilon", &initEpsilon, PHDOUBLE, 1);
                    }
                }
                else
                {
                    std::cout << "Please assign a value of Epsilon on inlet" << std::endl;;
                }

                if (FaceFlux[iFacelocal] < 0.0)
                {

                    kinetic[re] = initKinetic;
                    epsilon[re] = initEpsilon;
                }
                else
                {

                    kinetic[re] = kinetic[le];
                    epsilon[re] = epsilon[le];
                }
            }
        }
        else if (bcType == PHENGLEI::OUTFLOW)
        {
            RDouble initKinetic = 0.0;
            RDouble initEpsilon = 0.0;

            for (std::size_t iFace = 0; iFace < faceIndex->size(); ++ iFace)
            {
                int iFacelocal = (*faceIndex)[iFace];
                int le = leftCellOfFace[iFacelocal];
                int re = rightCellOfFace[iFacelocal];

                if (bcData)
                {
                    if (bcData->IsExist("initKinetic", PHDOUBLE, 1))
                    {
                        bcData->GetData("initKinetic", &initKinetic, PHDOUBLE, 1);
                    }
                }
                else
                {
                    std::cout << " No initKinetic initial value has assigned to the boundary " << std::endl;
                }

                if (bcData)
                {
                    if (bcData->IsExist("initEpsilon", PHDOUBLE, 1))
                    {
                        bcData->GetData("initEpsilon", &initEpsilon, PHDOUBLE, 1);
                    }
                }
                else
                {
                    std::cout << "Please assign a value of Epsilon on inlet" << std::endl;;
                }

                if (FaceFlux[iFacelocal] > 0.0)
                {
                    kinetic[re] = initKinetic;
                    epsilon[re] = initEpsilon;
                }
                else
                {
                    kinetic[re] = kinetic[le];
                    epsilon[re] = epsilon[le];
                }
            }
        }
        else if (bcType == PHENGLEI::PRESSURE_INLET)
        {
            RDouble kb = 0.0;
            for (std::size_t iFace = 0; iFace < faceIndex->size(); ++ iFace)
            {
                int iFacelocal = (*faceIndex)[iFace];
                int le = leftCellOfFace[iFacelocal];
                int re = rightCellOfFace[iFacelocal];
                if (bcData)
                {
                    if (bcData->IsExist("initKinetic", PHDOUBLE, 1))
                    {
                        bcData->GetData("initKinetic", &kb, PHDOUBLE, 1);
                    }
                }
                else
                {
                    std::cout << " No initial value has assigned to the boundary " << std::endl;
                }
                if (FaceFlux[iFacelocal] < 0)
                {
                    kinetic[re] = kb;
                }
                else
                {
                    kinetic[re] = kinetic[le];
                }
            }
        }
        else if (bcType == PHENGLEI::PRESSURE_OUTLET)
        {
            RDouble kb = 0.0;
            for (std::size_t iFace = 0; iFace < faceIndex->size(); ++ iFace)
            {
                int iFacelocal = (*faceIndex)[iFace];
                int le = leftCellOfFace[iFacelocal];
                int re = rightCellOfFace[iFacelocal];

                if (FaceFlux[iFacelocal] < 0.0)
                {
                    if (bcData)
                    {
                        if (bcData->IsExist("init" + varNameK, PHDOUBLE, 1))
                        {
                            bcData->GetData("init" + varNameK, &kb, PHDOUBLE, 1);
                            kinetic[re] = RDouble(kb);
                        }
                        else
                        {
                            std::cout << " Please assign a value " << std::endl;
                        }
                    }
                }
                else
                {
                    kinetic[re] = kinetic[le];
                }
            }
        }
        else 
        {
            for (size_t iFace = 0; iFace < faceIndex->size(); ++ iFace)
            {
                int iFacelocal = (*faceIndex)[iFace];
                int le = leftCellOfFace[iFacelocal];
                int re = rightCellOfFace[iFacelocal];
                kinetic[re] = kinetic[le];
                epsilon[re] = epsilon[le];
            }
        }
    }
    CommunicateAnInterfaceVar(kinetic);
    CommunicateAnInterfaceVar(epsilon);
}

void IncomKETurbKEqCalculator::calWallBC(Grid *gridIn, vector<int> *faceIndex, Data_Param *bcData, RDouble *diagMatrixCoeff, RDouble *bCoeff, RDouble *q,
                                        RDouble *dqdx, RDouble *dqdy, RDouble *dqdz, RDouble *mu, RDouble matrixCoeff, int iVariable)
{
    
}
}