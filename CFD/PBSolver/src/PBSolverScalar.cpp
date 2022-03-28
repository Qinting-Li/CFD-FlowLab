#include "PBSolverScalar.h"
#include "PBTransient.h"
#include "PBBoundary.h"
#include "PBGrad.h"
#include "PBConvection.h"
#include "PBDiffusion.h"
#include "PBSource.h"

namespace PHSPACE
{
PBSolverScalar::PBSolverScalar(const string &name):PBSolver(name),scalarName("phi")
{
    iterSolv = PBIterSolvFactory::CreatePBIterSolv(ObjGetStrParaFromDB(Name() + "_IterSolv", "BiCGSTAB"));
    maxSweep = ObjGetIntParaFromDB(Name() + "_MaxSweep", 30);
    iterSolvTol = ObjGetDoubleParaFromDB(Name() + "_IterSolvTol", 1e-3);
    maxIterStep = ObjGetIntParaFromDB(Name() + "_MaxIterStep", 1);
    minIterStep = ObjGetIntParaFromDB(Name() + "_MinIterStep", 1);
}

PBSolverScalar::~PBSolverScalar()
{
    if (iterSolv != NULL)
    {
        delete iterSolv;
    }
}

void PBSolverScalar::InitMemory()
{
    scalarNameList = new string[GetNSolvers()];
    GlobalDataBase::GetData(GetControllerName() + "_SCALAR_NAME", scalarNameList, PHSTRING, GetNSolvers());
    scalarName = scalarNameList[0];
    SetPrintName(scalarName); 
    SetGradientPhiName(scalarName);
    boundaryCalcConfigName = ScalarName() + "BoundaryType";
    convCalcConfigName = ScalarName() + "ConvCalc";
    diffCalcConfigName = ScalarName() + "DiffCalc";
    tranCalcConfigName = ScalarName() + "TranCalc";
    sourceCalcConfigName = ScalarName() + "SourceCalc";
    gradCalcConfigName = ScalarName() + "GradCalc";
    setFieldCalcConfigName;
}

void PBSolverScalar::AllocateMemory()
{
    grid = UnstructGridCast(GetGrid());

    //! Boundary calculation object creation.
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

        Data_Param *bcData = bcRegion->GetBCParamDataBase();
        std::string boundType;

        if (bcData)
        {
            if (bcData->IsExist(boundaryCalcConfigName, PHSTRING, 1))
            {
                bcData->GetData(boundaryCalcConfigName, &boundType, PHSTRING, 1);
            }
            PBBoundary *bound = PBBoundaryFactory::CreatePBBoundary(boundType);
            boundCalc.push_back(bound);
        }
    }

    //! Create a flow term calculation object.
    convCalc = PBConvectionFactory::CreatePBConvection(ObjGetStrParaFromDB(convCalcConfigName, STR(CONVECTION_EMPTY)));

    //! Diffusion term calculation object creation.
    diffCalc = PBDiffusionFactory::CreatePBDiffusion(ObjGetStrParaFromDB(diffCalcConfigName, STR(DIFFUSION_EMPTY)));

    //! Source term calculation object creation.
    int sourceCalcSize = GlobalDataBase::GetSizeOfData(sourceCalcConfigName);

    if (sourceCalcSize > 0)
    {
        std::string *sourceCalcs = new std::string[sourceCalcSize];
        GlobalDataBase::GetData(sourceCalcConfigName, sourceCalcs, PHSTRING, sourceCalcSize);

        for (int i = 0; i < sourceCalcSize; ++ i)
        {
            PBSource *source = PBSourceFactory::CreatePBSource(sourceCalcs[i]);
            sourceCalc.push_back(source);
        }

        delete[] sourceCalcs;
    }

    //! Transient term calculation object.
    if (!ObjGetIntParaFromDB("iunsteady", 0))
    {
        tranCalc = PBTransientFactory::CreatePBTransient(STR(TRANSIENT_EMPTY));
    }
    else
    {
        tranCalc = PBTransientFactory::CreatePBTransient(ObjGetStrParaFromDB(tranCalcConfigName, STR(TRANSIENT_EMPTY)));
    }

    //! Gradient calculation object.
    PBGrad *gradCalc = PBGradFactory::CreatePBGrad(ObjGetStrParaFromDB(gradCalcConfigName, "GAUSS"));
    this->gradCalc[ScalarName()] = gradCalc;

    //! Allocate variable memory.
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFaces = grid->GetNBoundFace();
    int nTotalFaces = grid->GetNTotalFace();
    int nTotalInFaces = nTotalFaces - nBoundFaces;
    this->nTotalInFaces = nTotalInFaces;
    unsigned int nTotal = nTotalCell + nBoundFaces;
    nTotalInFaces = nTotalInFaces;

    //! Variables registered on the grid.
    initPhi = ObjGetDoubleParaFromDB(ScalarName() + "InitValue", 0.0);
    phi = ObjAllocateMemory(ScalarName(), nTotal, initPhi, grid);
    phib = phi + nTotalCell;
    grid->UpdateDataPtr(ScalarName() + "b", phib);

    initRho = ObjGetDoubleParaFromDB("initRho", 1.0);
    rho = ObjAllocateMemory("rho", nTotal, initRho, grid);
    rhob = rho + nTotalCell;
    rhof = ObjAllocateMemory("rhof", nTotalFaces, initRho, grid);

    initMu = ObjGetDoubleParaFromDB("initMu", 1.0);
    mu = ObjAllocateMemory("mu", nTotal, initMu, grid);
    mub = mu + nTotalCell;
    diffCoe = ObjAllocateMemory("diffCoe", nTotal, initMu, grid);
    diffCoeb = diffCoe + nTotalCell;

    urf = ObjGetDoubleParaFromDB(ScalarName() + "Urf", 0.0);
    res = ObjGetDoubleParaFromDB(ScalarName() + "Res", 1e-2);

    phic = ObjAllocateMemory(ScalarName() + "c", nTotalCell, 0.0, grid);

    ap = ObjAllocateMemory("ap", nTotalCell, 0.0, grid);
    ajLeft = ObjAllocateMemory("ajLeft", nTotalInFaces, 0.0, grid);
    ajRight = ObjAllocateMemory("ajRight", nTotalInFaces, 0.0, grid);
    sc = ObjAllocateMemory("sc", nTotalCell, 0.0, grid);

    flux = ObjAllocateMemory("flux", nTotalFaces, 0.0, grid);
    faceFlux = flux + grid->GetNBoundFace();

    dphidx = ObjAllocateMemory("dphidx", nTotal, 0.0, grid);
    dphidy = ObjAllocateMemory("dphidy", nTotal, 0.0, grid);
    dphidz = ObjAllocateMemory("dphidz", nTotal, 0.0, grid);

    deltL = ObjAllocateMemory("deltL", nTotalFaces, 0.0, grid);
    deltR = ObjAllocateMemory("deltR", nTotalFaces, 0.0, grid);

    grid->FaceWeight(deltL, deltR, 0, nTotalFaces);
}

void PBSolverScalar::InitField()
{
    PHString1D phiNameList;
    phiNameList.push_back(ScalarName());
    PHDouble1D phiInitValueList;
    phiInitValueList.push_back(initPhi);
    tranCalc->AllocateMemory(grid, phiNameList, phiInitValueList);
}

void PBSolverScalar::IncompressibleInitial()
{
    UpdateBoundary();
}

void PBSolverScalar::UpdateUnsteadyVariable(Grid *grid)
{
    UnstructGrid *gridIn = UnstructGridCast(grid); 
    std::vector<std::string> phiNameList;
    phiNameList.push_back(ScalarName());
    tranCalc->SaveOldTimeValue(gridIn, phiNameList);
}

bool PBSolverScalar::SolutionIsConverged()
{
    return true;
}

void PBSolverScalar::OutputResidual(ofstream &file, int flag)
{

}

void PBSolverScalar::InitIncompressibleResultFile(Grid *grid, string fileName)
{
    string outputFileName = string("results/") + fileName;
    ofstream file(outputFileName.c_str(), ios::app);

    int innerIter = GlobalDataBase::GetIntParaFromDB("outnstep");
    if (innerIter == 1)
    {
        OutputResidual(file, 1);
    }
}

void PBSolverScalar::DumpIncompressibleResultFile(Grid *grid, string fileName)
{
    string outputFileName = string("results/") + fileName;
    ofstream file(outputFileName.c_str(), ios::app);

    OutputResidual(file, 2);
}

void PBSolverScalar::Discrt()
{
    //! Initialize the related variables of linear equations.
    InitMatrix();

    //! Calculate the diffusion term.
    CalcDiffusionTerm();

    //! Convection term calculation.
    CalcConvectionTerm();

    //! Boundary calculation.
    CalcBoundaryConditions();

    //! Source term calculation.
    CalcSourceTerm();

    //! Transient term calculation.
    CalcTransientTerm();

    //! Form ap and sc.
    PrimeEquations();

    //! Calculate the residual.
    CalcRes();

    //! Relax ap.
    Relax();
}

void PBSolverScalar::Update()
{
    int nTotalCell = grid->GetNTotalCell();
    UpdatePhi(ScalarName(), 0, nTotalCell, phic, phi);
}

void PBSolverScalar::UpdateBoundary()
{
    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegionUnstruct = unstructBCSet->GetnBCRegion();

    for (int iBCRegionUnstruct = 0; iBCRegionUnstruct < nBCRegionUnstruct; iBCRegionUnstruct++)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegionUnstruct);
        int bcType = bcRegion->GetBCType();
        vector<int> *faceIndex = bcRegion->GetFaceIndex();
        Data_Param *bcData = bcRegion->GetBCParamDataBase();

        if (bcType < 0)
        {
            continue;
        }

        //! Update boundary conditions.
        static_cast<PBBoundary *>(boundCalc[iBCRegionUnstruct])->Update(grid, ScalarName(), faceIndex, bcData);
    }
}

void PBSolverScalar::SetOutputs(std::vector<std::string> &phiNameList)
{
    phiNameList.push_back(ScalarName());
}

void PBSolverScalar::InitMatrix()
{
    int nTotalCell = grid->GetNTotalCell();
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        ap[iCell] = 0.0;
        sc[iCell] = 0.0;
        phic[iCell] = 0.0;
    }

    for (int iFace = 0; iFace < NTotalInFace(); ++ iFace)
    {
        ajLeft[iFace] = 0.0;
        ajRight[iFace] = 0.0;
    }
}

void PBSolverScalar::CalcGrad()
{
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;
    for (int iCell = 0; iCell < nTotal; ++ iCell)
    {
        dphidx[iCell] = 0.0;
        dphidy[iCell] = 0.0;
        dphidz[iCell] = 0.0;
    }
 
    gradCalc[ScalarName()]->Calc(grid, GradientPhiName(), "dphidx", "dphidy", "dphidz");
}

void PBSolverScalar::CalcBoundaryConditions()
{
    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegionUnstruct = unstructBCSet->GetnBCRegion();

    for (int iBCRegionUnstruct = 0; iBCRegionUnstruct < nBCRegionUnstruct; iBCRegionUnstruct++)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegionUnstruct);
        int bcType = bcRegion->GetBCType();
        vector<int> *faceIndex = bcRegion->GetFaceIndex();
        Data_Param *bcData = bcRegion->GetBCParamDataBase();

        if (bcType < 0)
        {
            continue;
        }

        static_cast<PBBoundary *>(boundCalc[iBCRegionUnstruct])->Discrt(grid, ScalarName(), faceIndex, bcData, diffCoe);
    }
}

void PBSolverScalar::CalcTransientTerm()
{
    
    tranCalc->Discrt(grid, ScalarName(), ScalarName() + "Old", ScalarName() + "OldOld", "rho", "rho");
}

void PBSolverScalar::CalcConvectionTerm()
{
    convCalc->Discrt(grid, ScalarName(), "ap", "sc");
}

void PBSolverScalar::CalcDiffusionTerm()
{
    if (diffCalc->Name() != STR(DIFFUSION_EMPTY))
    {
        CalcGrad();
    }

    diffCalc->Discrt(grid, ScalarName(), diffCoe, dphidx, dphidy, dphidz, ajLeft, ajRight, sc);
}

void PBSolverScalar::CalcSourceTerm()
{
    for (std::vector<PBSource *>::iterator iter = sourceCalc.begin(); iter != sourceCalc.end(); ++iter)
    {
        (*iter)->Discrt(grid, ScalarName());
    }
}

//! Correct;considering the influence of difference on coefficient.
void PBSolverScalar::PrimeEquations()
{
    int nTotalCell = grid->GetNTotalCell();
    int nTotalFace = grid->GetNTotalFace();
    int nBoundFace = grid->GetNBoundFace();
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        int iFacelocal = iFace - nBoundFace;
        int icL = leftCellOfFace[iFace];
        int icR = rightCellOfFace[iFace];
        ap[icL] += ajLeft[iFacelocal];
        ap[icR] += ajRight[iFacelocal];
        sc[icL] += ajLeft[iFacelocal] * phi[icR];
        sc[icR] += ajRight[iFacelocal] * phi[icL];
    }

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        sc[iCell] -= ap[iCell] * phi[iCell];
    }
}

void PBSolverScalar::CalcRes()
{
    int nTotalCell = grid->GetNTotalCell();
    resNow = 0.0;

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        resNow += std::abs(sc[iCell]);
    }

    grid->UpdateData(ScalarName() + "ResNow", &resNow, PHDOUBLE, 1);
}

void PBSolverScalar::Relax()
{
    int nTotalCell = grid->GetNTotalCell();
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        ap[iCell] = (1.0 + urf) * ap[iCell];
    }
}

void PBSolverScalar::SolveScalarEquation(Region &region)
{
    int nZones = PHMPI::GetNumberofGlobalZones();
    std::vector<int> zoneIdContainer(nZones);
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        zoneIdContainer[iZone] = iZone;
    }

    //! Solve the momentum equation and discrete equations on region.
    iterSolv->Run(region, zoneIdContainer, GetPrintName(), "ap", "ajLeft", "ajRight", "ajbLeft",
        "ajbRight", "sc", ScalarName() + "c",  maxSweep, iterSolvTol );
}

void PBSolverScalar::UpdatePhi(string phiName, int cellBegin, int cellEnd, RDouble *phic, RDouble *phi)
{
    for (int iCell = cellBegin; iCell < cellEnd; ++ iCell)
    {
        phi[iCell] += phic[iCell];
    }
}

void PBSolverScalar::SolutionIsConverged(int iZone, bool *flag)
{
    int innerIter = GlobalDataBase::GetIntParaFromDB("outnstep");
    if (innerIter > 0)
    {
        if ((resNow < res))
        {
            *flag = true;
        }
        else
        {
            *flag = false;
        }
    }
    else
    {
        *flag = false;
    }
}

DEFINE_PBSOLVER_REGISTER_FUNCTION(PBSolver, PBSolverScalar, SOLVER_SCALAR);
}
