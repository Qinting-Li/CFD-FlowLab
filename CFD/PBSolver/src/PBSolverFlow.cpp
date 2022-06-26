#include "PBSolverFlow.h"
#include "TK_Log.h"
#include "Post_Probes.h"
#include "Post_Visual.h"
#include "PBIterSolv.h"
#include "PBBoundary.h"
#include "PBTransient.h"
#include "PBConvection.h"
#include "PBDiffusion.h"
#include "PBSource.h"
#include "PBGrad.h"
#include "PBSolver.h"
#include "Geo_UnstructBC.h"
#include "Data_Param.h"
#include "PBBoundaryFlow.h"
#include "Post_WriteTecplot.h"

namespace PHSPACE
{

PBSolverFlow::PBSolverFlow(const string &name):PBSolver(name)
{
    phiNameList.push_back("U");
    phiNameList.push_back("V");
    phiNameList.push_back("W");
    phiNameList.push_back("P");

    //! Linear equations iterative solver configuration.
    iterSolvU = PBIterSolvFactory::CreatePBIterSolv(GlobalDataBase::GetStrParaFromDB("iterSolvU"));
    iterSolvV = PBIterSolvFactory::CreatePBIterSolv(GlobalDataBase::GetStrParaFromDB("iterSolvV"));
    iterSolvW = PBIterSolvFactory::CreatePBIterSolv(GlobalDataBase::GetStrParaFromDB("iterSolvW"));
    iterSolvP = PBIterSolvFactory::CreatePBIterSolv(GlobalDataBase::GetStrParaFromDB("iterSolvP"));
    maxSweepU = GlobalDataBase::GetIntParaFromDB("maxSweepU");
    maxSweepV = GlobalDataBase::GetIntParaFromDB("maxSweepV");
    maxSweepW = GlobalDataBase::GetIntParaFromDB("maxSweepW");
    maxSweepP = GlobalDataBase::GetIntParaFromDB("maxSweepP");
    iterSolvTolU = GlobalDataBase::GetDoubleParaFromDB("iterSolvTolU");
    iterSolvTolV = GlobalDataBase::GetDoubleParaFromDB("iterSolvTolV");
    iterSolvTolW = GlobalDataBase::GetDoubleParaFromDB("iterSolvTolW");
    iterSolvTolP = GlobalDataBase::GetDoubleParaFromDB("iterSolvTolP");
}

PBSolverFlow::~PBSolverFlow()
{
    delete iterSolvU;
    delete iterSolvV;
    delete iterSolvW;
    delete iterSolvP;
}

void PBSolverFlow::AllocateGlobalVar(Grid *grid)
{
    PrintToWindow("PBSolverFlow::AllocateGlobalVar\n");
}

void PBSolverFlow::DeAllocateGlobalVar(Grid *grid)
{
    PrintToWindow("PBSolverFlow::DeAllocateGlobalVar\n");
}

void PBSolverFlow::ZeroResiduals(Grid *grid)
{
    PrintToWindow("PBSolverFlow::ZeroResiduals\n");
}

void PBSolverFlow::TimeStep(Grid *grid)
{
    PrintToWindow("PBSolverFlow::TimeStep\n");
}

void PBSolverFlow::ComputeResidual(Grid *grid_in, FieldProxy *rhsProxy)
{
    PrintToWindow("PBSolverFlow::ComputeResidual\n");
}

void PBSolverFlow::SourceFlux(Grid *grid)
{
    PrintToWindow("PBSolverFlow::SourceFlux\n");
}

void PBSolverFlow::InviscidFlux(Grid *grid)
{
    PrintToWindow("PBSolverFlow::InviscidFlux\n");
}

void PBSolverFlow::ViscousFlux(Grid *grid)
{
    PrintToWindow("PBSolverFlow::ViscousFlux\n");
}

void PBSolverFlow::RungeKutta(Grid *grid, FieldProxy *rhsProxy)
{
    PrintToWindow("PBSolverFlow::RungeKutta\n");
}

void PBSolverFlow::Post()
{
    PrintToWindow("PBSolverFlow::Post\n");
}

void PBSolverFlow::PostSolve(Grid *grid, int stage, int level)
{
    int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");

    if (isUnsteady == 1)
    {
        RDouble currentTime = PHSPACE::GlobalDataBase::GetDoubleParaFromDB("currentTime");
        Post_WriteTecplot generateVisualFile;
        generateVisualFile.Run();
    }
}

void PBSolverFlow::InitIncompressibleResultFile(Grid *grid, string fileName)
{
    string outputFileName = string("results/") + fileName;
    ofstream file(outputFileName.c_str(), ios::app);

    int innerIter = GlobalDataBase::GetIntParaFromDB("outnstep");
    if (innerIter == 1)
    {

        if (ObjGetIntParaFromDB("isUnsteady", 0))
        {
            file << "Time" << "\t";
        }

        file << "IterStep" << "\t";

        OutputResidual(file, 1);
    }
}

void PBSolverFlow::DumpIncompressibleResultFile(Grid *grid, string fileName)
{
    string outputFileName = string("results/") + fileName;
    ofstream file(outputFileName.c_str(), ios::app);

    int innerIter = GlobalDataBase::GetIntParaFromDB("outnstep");
    if (this->ObjGetIntParaFromDB("isUnsteady", 0))
    {
        RDouble currentTime = GlobalDataBase::GetDoubleParaFromDB("currentTime");

        file << currentTime << "\t";
    }

    file << innerIter << "\t";

    OutputResidual(file, 2);
}

RDouble PBSolverFlow::UnsteadyConvergence(Grid *grid)
{
    if (GlobalDataBase::GetIntParaFromDB("iunsteady") == 0) { return 0; }
    RDouble dt = GlobalDataBase::GetDoubleParaFromDB("dt");
    return 0.0;
}

void PBSolverFlow::InitMemory()
{
    AllocateMemory();
    InitField();
}

void PBSolverFlow::AllocateMemory()
{
    initRho = ObjGetDoubleParaFromDB("initRho", 1.0);
    initMu = ObjGetDoubleParaFromDB("initMu", 1.0);
    initU = ObjGetDoubleParaFromDB("initU", 0.0);
    initV = ObjGetDoubleParaFromDB("initV", 0.0);
    initW = ObjGetDoubleParaFromDB("initW", 0.0);
    initP = ObjGetDoubleParaFromDB("initP", 0.0);
    urfU = ObjGetDoubleParaFromDB("urfU", 1.0);
    urfV = ObjGetDoubleParaFromDB("urfV", 1.0);
    urfW = ObjGetDoubleParaFromDB("urfW", 1.0);
    urfP = ObjGetDoubleParaFromDB("urfP", 1.0);
    urfFlux = ObjGetDoubleParaFromDB("urfFlux", 0.8);
    resU = ObjGetDoubleParaFromDB("resU", 1e-3);
    resV = ObjGetDoubleParaFromDB("resV", 1e-3);
    resW = ObjGetDoubleParaFromDB("resW", 1e-3);
    resP = ObjGetDoubleParaFromDB("resP", 1e-3);
    bodyForceFlag = ObjGetIntParaFromDB("bodyForceFlag", 0);
    dt = ObjGetDoubleParaFromDB("dt", 1.0);
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
        std::string flowType;

        if (bcData)
        {
            if (bcData->IsExist("flowType", PHSTRING, 1))
            {
                bcData->GetData("flowType", &flowType, PHSTRING, 1);
            }
            PBBoundary *bound = PBBoundaryFlowFactory::CreatePBBoundaryFlow(flowType);
            boundCalc.push_back(bound);
        }
    }

    //! Create a flow item calculation object.
    convCalc = PBConvectionFactory::CreatePBConvection(GlobalDataBase::GetStrParaFromDB("flowConvCalc"));

    //! Diffusion item calculation object creation.
    diffCalc = PBDiffusionFactory::CreatePBDiffusion(GlobalDataBase::GetStrParaFromDB("flowDiffCalc"));

    //! Source item calculation object creation.
    int sourceCalcSize = GlobalDataBase::GetSizeOfData("flowSourceCalc");
    std::string *sourceCalcs = new std::string[sourceCalcSize];
    GlobalDataBase::GetData("flowSourceCalc", sourceCalcs, PHSTRING, sourceCalcSize);

    for (int i = 0; i < sourceCalcSize; ++ i)
    {
        PBSource *source = PBSourceFactory::CreatePBSource(sourceCalcs[i]);
        sourceCalc.push_back(source);
    }

    delete[] sourceCalcs;

    //! Transient term calculation object.
    if (!ObjGetIntParaFromDB("iunsteady", 0))
    {
        tranCalc = PBTransientFactory::CreatePBTransient("EMPTY");
    }
    else
    {
        int timeStepNow = 0;
        GlobalDataBase::UpdateData("timeStepNow", &timeStepNow, PHINT, 1);
        RDouble startTime = GlobalDataBase::GetDoubleParaFromDB("startTime");
        GlobalDataBase::UpdateData("currentTime", &startTime, PHDOUBLE, 1);
        tranCalc = PBTransientFactory::CreatePBTransient(ObjGetStrParaFromDB("flowTranCalc", "EMPTY"));
    }

    //£¡Gradient calculation object.
    for (std::vector<string>::iterator iter = phiNameList.begin(); iter != phiNameList.end(); ++ iter)
    {
        std::string gradCalcConfigName = (*iter) + "GradCalc";
        PBGrad *gradCalc = PBGradFactory::CreatePBGrad(ObjGetStrParaFromDB(gradCalcConfigName, "GAUSS"));
        this->gradCalc[(*iter)] = gradCalc;
    }

    // Allocate variable memory
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFaces = grid->GetNBoundFace();
    unsigned int nTotal = nTotalCell + nBoundFaces;
    int nTotalInFaces = grid->GetNTotalFace() - nBoundFaces;
    this->nTotalInFaces = nTotalInFaces;

    // Variables registered on Grid
    u = ObjAllocateMemory("U", nTotal, initU, grid); // Boundary values come last
    ub = u + nTotalCell;

    v = ObjAllocateMemory("V", nTotal, initV, grid);
    vb = v + nTotalCell;

    w = ObjAllocateMemory("W", nTotal, initW, grid);
    wb = w + nTotalCell;

    p = ObjAllocateMemory("P", nTotal, initP, grid);
    pb = p + nTotalCell;

    uc = ObjAllocateMemory("uc", nTotalCell, 0.0, grid);

    vc = ObjAllocateMemory("vc", nTotalCell, 0.0, grid);

    wc = ObjAllocateMemory("wc", nTotalCell, 0.0, grid);

    pc = ObjAllocateMemory("pc", nTotal, 0.0, grid);
    pcb = pc + nTotalCell;

    uf = ObjAllocateMemory("Uf", grid->GetNTotalFace(), 0.0, grid);

    apU = ObjAllocateMemory("apU", nTotalCell, 0.0, grid);

    apV = ObjAllocateMemory("apV", nTotalCell, 0.0, grid);

    apW = ObjAllocateMemory("apW", nTotalCell, 0.0, grid);

    ap = ObjAllocateMemory("ap", nTotalCell, 0.0, grid);

    scU = ObjAllocateMemory("scU", nTotalCell, 0.0, grid);

    scV = ObjAllocateMemory("scV", nTotalCell, 0.0, grid);

    scW = ObjAllocateMemory("scW", nTotalCell, 0.0, grid);

    sc = ObjAllocateMemory("sc", nTotalCell, 0.0, grid);

    ajLeft = ObjAllocateMemory("ajLeft", nTotalInFaces, 0.0, grid);

    ajRight = ObjAllocateMemory("ajRight", nTotalInFaces, 0.0, grid);

    flux = ObjAllocateMemory("flux", grid->GetNTotalFace(), 0.0, grid);
    faceFlux = flux + grid->GetNBoundFace();

    if (bodyForceFlag == 1)
    {
        bfU = ObjAllocateMemory("bfU", nTotalCell, 0.0, grid);

        bfV = ObjAllocateMemory("bfV", nTotalCell, 0.0, grid);

        bfW = ObjAllocateMemory("bfW", nTotalCell, 0.0, grid);
    }

    dun = ObjAllocateMemory("dun", grid->GetNTotalFace(), 0.0, grid);

    rho = ObjAllocateMemory("rho", nTotal, initRho, grid);

    rhof = ObjAllocateMemory("rhof", grid->GetNTotalFace(), 0.0, grid);

    mu = ObjAllocateMemory("mu", nTotal, initMu, grid);

    hu = ObjAllocateMemory("hu", nTotalCell, 0.0, grid);

    hv = ObjAllocateMemory("hv", nTotalCell, 0.0, grid);

    hw = ObjAllocateMemory("hw", nTotalCell, 0.0, grid);

    dudx = ObjAllocateMemory("dudx", nTotal, 0.0, grid);

    dudy = ObjAllocateMemory("dudy", nTotal, 0.0, grid);

    dudz = ObjAllocateMemory("dudz", nTotal, 0.0, grid);

    dvdx = ObjAllocateMemory("dvdx", nTotal, 0.0, grid);

    dvdy = ObjAllocateMemory("dvdy", nTotal, 0.0, grid);

    dvdz = ObjAllocateMemory("dvdz", nTotal, 0.0, grid);

    dwdx = ObjAllocateMemory("dwdx", nTotal, 0.0, grid);

    dwdy = ObjAllocateMemory("dwdy", nTotal, 0.0, grid);

    dwdz = ObjAllocateMemory("dwdz", nTotal, 0.0, grid);

    dpdx = ObjAllocateMemory("dpdx", nTotal, 0.0, grid);

    dpdy = ObjAllocateMemory("dpdy", nTotal, 0.0, grid);

    dpdz = ObjAllocateMemory("dpdz", nTotal, 0.0, grid);

    dpcdx = ObjAllocateMemory("dpcdx", nTotal, 0.0, grid);

    dpcdy = ObjAllocateMemory("dpcdy", nTotal, 0.0, grid);

    dpcdz = ObjAllocateMemory("dpcdz", nTotal, 0.0, grid);

    if (bodyForceFlag)
    {
        RDouble *dbfdx = ObjAllocateMemory("dbfdx", grid->GetNTotalCell(), 0.0, grid);
        RDouble *dbfdy = ObjAllocateMemory("dbfdy", grid->GetNTotalCell(), 0.0, grid);
        RDouble *dbfdz = ObjAllocateMemory("dbfdz", grid->GetNTotalCell(), 0.0, grid);
    }
    deltL = ObjAllocateMemory("deltL", grid->GetNTotalFace(), 0.0, grid);
    deltR = ObjAllocateMemory("deltR", grid->GetNTotalFace(), 0.0, grid);
    grid->FaceWeight(deltL, deltR, 0, grid->GetNTotalFace());
}

void PBSolverFlow::InitField()
{
    PHString1D phiNameList;
    phiNameList.resize(4);
    PHDouble1D phiInitValueList;
    phiInitValueList.resize(4);

    phiNameList[0] = "U";
    phiNameList[1] = "V";
    phiNameList[2] = "W";
    phiNameList[3] = "P";

    phiInitValueList[0] = initU;
    phiInitValueList[1] = initV;
    phiInitValueList[2] = initW;
    phiInitValueList[3] = initP;

    tranCalc->AllocateMemory(grid, phiNameList, phiInitValueList);
}

void PBSolverFlow::IncompressibleInitial()
{
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFaces = grid->GetNBoundFace();
    int nTotalFaces = grid->GetNTotalFace();
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    RDouble *area    = grid->GetFaceArea();
    RDouble *xfn = grid->GetFaceNormalX();
    RDouble *yfn = grid->GetFaceNormalY();
    RDouble *zfn = grid->GetFaceNormalZ();

    //! Internal unit surface circulation
    for (int iFace = nBoundFaces; iFace < nTotalFaces; ++ iFace)
    {
        RDouble deltL = this->deltL[iFace];
        RDouble deltR = this->deltR[iFace];
        int le = leftCellOfFace[iFace];
        int re = rightCellOfFace[iFace];
        RDouble rhof = deltL * rho[le] + deltR * rho[re];
        RDouble uf = deltL * u[le] + deltR * u[re];
        RDouble vf = deltL * v[le] + deltR * v[re];
        RDouble wf = deltL * w[le] + deltR * w[re];
        RDouble faceArea = area[iFace];
        RDouble faceNormalX = xfn[iFace];
        RDouble faceNormalY = yfn[iFace];
        RDouble faceNormalZ = zfn[iFace];

        RDouble ax = faceArea * faceNormalX;
        RDouble ay = faceArea * faceNormalY;
        RDouble az = faceArea * faceNormalZ;
        this->rhof[iFace] = rhof;
        flux[iFace] = rhof * (uf * ax + vf * ay + wf * az);
        dun[iFace] = 0;
    }

    if (bodyForceFlag == 1)
    {
        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            bfU[iCell] = 0.0;
            bfV[iCell] = 0.0;
            bfW[iCell] = 0.0;
            dpdx[iCell] = 0.0;
            dpdy[iCell] = 0.0;
            dpdz[iCell] = 0.0;
        }
    }

    for (vector<PBSource *>::iterator iter = sourceCalc.begin(); iter != sourceCalc.end(); ++ iter)
    {
        (*iter)->InitMemory(grid, phiNameList[0]);
    }

    UpdateBoundaryPressureAndVelocity();

    //£¡ Initialize boundary flux
    for (int iFace = 0; iFace < nBoundFaces; ++ iFace)
    {
        int le = leftCellOfFace[iFace];
        int re = rightCellOfFace[iFace];
        RDouble uf = u[re];
        RDouble vf = v[re];
        RDouble wf = w[re];
        RDouble faceArea = area[iFace];
        RDouble faceNormalX = xfn[iFace];
        RDouble faceNormalY = yfn[iFace];
        RDouble faceNormalZ = zfn[iFace];

        RDouble ax = faceArea * faceNormalX;
        RDouble ay = faceArea * faceNormalY;
        RDouble az = faceArea * faceNormalZ;
        rhof[iFace] = rho[le];
        flux[iFace] = rho[le] * (uf * ax + vf * ay + wf * az);
        dun[iFace] = 0;
    }
}

void PBSolverFlow::UpdateUnsteadyVariable(Grid *grid)
{
    UnstructGrid *gridIn = UnstructGridCast(grid);
    tranCalc->SaveOldTimeValue(gridIn, phiNameList);
}

void PBSolverFlow::DiscrtMomentumEquations()
{
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFaces = grid->GetNBoundFace();
    int nTotalFaces = grid->GetNTotalFace();

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        apU[iCell] = 0.0;
        apV[iCell] = 0.0;
        apW[iCell] = 0.0;
        scU[iCell] = 0.0;
        scV[iCell] = 0.0;
        scW[iCell] = 0.0;
        uc[iCell] = 0.0;
        vc[iCell] = 0.0;
        wc[iCell] = 0.0;
    }

    for (int iFace = 0; iFace < NTotalInFace(); ++ iFace)
    {
        ajLeft[iFace] = 0.0;
        ajRight[iFace] = 0.0;
    }

    if (bodyForceFlag != 0)
    {
        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            bfU[iCell] = 0.0;
            bfV[iCell] = 0.0;
            bfW[iCell] = 0.0;
        }
    }

    for (vector<PBSource *>::iterator iter = sourceCalc.begin(); iter != sourceCalc.end(); ++ iter)
    {
        (*iter)->InitMemory(grid, phiNameList[0]);
    }

    //! Calculate the pressure gradient.
    CalcPressureGrad();

    //! Calculate the velocity gradient.
    CalcVelocityGrad();

    //! Calculate the diffusion term.
    CalcVelocityDiffusionTerm();

    //! Calculate the convection term.
    CalcVelocityConvectionTerm();

    //! Handle boundary conditions.
    CalcVelocityBoundaryConditions();

    //! Calculate the source term.
    CalcVelocitySourceTerm();

    //! Calculate the transient term.
    CalcVelocityTransientTerm();

    //! Form matrix
    PrimeVelocityEquations();

    //! Calculate the residual.
    CalcVelocityRes();

    //! Relax ap.
    RelaxVelocity();
}

void PBSolverFlow::SolveMomentumEquations(Region &region)
{
    int nZones = PHMPI::GetNumberofGlobalZones();
    std::vector<int> zoneIdContainer(nZones);
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        zoneIdContainer[iZone] = iZone;
    }

    //! Solve the momentum equation and discrete equations on region.
    iterSolvU->Run(region, zoneIdContainer, "U", "apU", "ajLeft", "ajRight", "ajbLeft", "ajbRight", "scU", "uc", maxSweepU, iterSolvTolU);

    iterSolvV->Run(region, zoneIdContainer, "V", "apV", "ajLeft", "ajRight", "ajbLeft", "ajbRight", "scV", "vc", maxSweepV, iterSolvTolV);

    iterSolvW->Run(region, zoneIdContainer, "W", "apW", "ajLeft", "ajRight", "ajbLeft", "ajbRight", "scW", "wc", maxSweepW, iterSolvTolW);

    std::vector<int>().swap(zoneIdContainer);
}

void PBSolverFlow::DiscrtPressureCorrectionEquation()
{
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFaces = grid->GetNBoundFace();
    int nTotalFaces = grid->GetNTotalFace();

    //! The correlation coefficient is stored for the pressure correction equation.
    CoefPrepForPres();

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        ap[iCell] = 0.0;
        sc[iCell] = 0.0;
        pc[iCell] = 0.0;
    }

    for (int iFace = 0; iFace < NTotalInFace(); ++ iFace)
    {
        ajLeft[iFace] = 0.0;
        ajRight[iFace] = 0.0;
    }

    //! Form ap,aj and sc.
    PrimePressureEquations();

    //! Calculate the residual.
    CalcPressureRes();
}

void PBSolverFlow::SolvePressureCorrectionEquation(Region &region)
{
    int nZones = PHMPI::GetNumberofGlobalZones();
    std::vector<int> zoneIdContainer(nZones);
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        zoneIdContainer[iZone] = iZone;
    }

    //! Solve the pressure correction equation and discrete equations on region.
    iterSolvP->Run(region, zoneIdContainer, "P", "ap", "ajLeft", "ajRight", "ajbLeft", "ajbRight", "sc", "pc", maxSweepP, iterSolvTolP);

    std::vector<int>().swap(zoneIdContainer);
}

void PBSolverFlow::UpdateVelocity()
{
    int nTotalCell = grid->GetNTotalCell();
    UpdatePhi("U", 0, nTotalCell, uc, u);
    UpdatePhi("V", 0, nTotalCell, vc, v);
    UpdatePhi("W", 0, nTotalCell, wc, w);
}

void PBSolverFlow::UpdatePressureAndVelocity()
{
    int nTotalCell = grid->GetNTotalCell();
    //! Solve the gradient of pressure correction.
    CalcPressureCorrectionGrad();

    //! Update pressure.
    if (ObjGetIntParaFromDB("refP", 0) == 1)
    {
        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            p[iCell] += urfP * (pc[iCell] - pc[0]);
        }
    }
    else
    {
        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            p[iCell] += urfP * pc[iCell];
        }
    }

    //! Update speed and flux.
    UpdateVelocityAndFlux();

    //! Update boundary pressure and velocity values.
    UpdateCorrectionBoundaryPressureAndVelocity();
}

void PBSolverFlow::UpdateBoundary()
{

}

void PBSolverFlow::ComputePostVisualVariables(Post_Visual *postVisualization)
{
    UnstructGrid *grid = UnstructGridCast(postVisualization->GetGrid());

    RDouble *u = reinterpret_cast<RDouble *> (grid->GetDataPtr("U"));
    RDouble *v = reinterpret_cast<RDouble *> (grid->GetDataPtr("V"));
    RDouble *w = reinterpret_cast<RDouble *> (grid->GetDataPtr("W"));
    RDouble *p = reinterpret_cast<RDouble *> (grid->GetDataPtr("P"));

    postVisualization->UpdateVisualNodeVarPtr("U", u);
    postVisualization->UpdateVisualNodeVarPtr("V", v);
    postVisualization->UpdateVisualNodeVarPtr("W", w);
    postVisualization->UpdateVisualNodeVarPtr("P", p);
}

void PBSolverFlow::ComputePostProbesVariables(Post_Probes* postProbesVar)
{

}

void PBSolverFlow::SolutionIsConverged(int iZone, bool *flag)
{
    int innerIter = GlobalDataBase::GetIntParaFromDB("outnstep");
    if (innerIter > 0)
    {
        if ((uResNow < resU) && (vResNow < resV) && (wResNow < resW) && (pResNow < resP))
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

void PBSolverFlow::OutputResidual(std::ofstream &file, int flag)
{
    file << setiosflags(ios::left);
    file << setprecision(5);
    file << setiosflags(ios::scientific);

    if (flag == 1)
    {
        file << " Res U" << '\t' << " Res V" << '\t' << " Res W" << '\t' << " Res P" << '\t';
    }
    else
    {
        RDouble uResNow = 0.0;
        RDouble vResNow = 0.0;
        RDouble wResNow = 0.0;
        RDouble pResNow = 0.0;

        RDouble uRes = 0.0;
        RDouble vRes = 0.0;
        RDouble wRes = 0.0;
        RDouble pRes = 0.0;

        grid->GetData("uResNow", &uRes, PHDOUBLE, 1);
        grid->GetData("vResNow", &vRes, PHDOUBLE, 1);
        grid->GetData("wResNow", &wRes, PHDOUBLE, 1);
        grid->GetData("pResNow", &pRes, PHDOUBLE, 1);

        uResNow = uRes;
        vResNow = vRes;
        wResNow = wRes;
        pResNow = pRes;

        file << uResNow << '\t' << vResNow << '\t' << wResNow << '\t' << pResNow << '\t';
    }
}

void PBSolverFlow::UpdatePhi(string phiName, int cellBegin, int cellEnd, RDouble *phic, RDouble *phi)
{
    for (int iCell = cellBegin; iCell < cellEnd; ++ iCell)
    {
        phi[iCell] += phic[iCell];
    }
}

void PBSolverFlow::UpdateProperties()
{

}

void PBSolverFlow::UpdateBoundaryPressureAndVelocity()
{
    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegionUnstruct = unstructBCSet->GetnBCRegion();
    int interfaceIndex = 0;

    for (int iBCRegionUnstruct = 0; iBCRegionUnstruct < nBCRegionUnstruct; iBCRegionUnstruct++)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegionUnstruct);
        int bcType = bcRegion->GetBCType();
        vector<int> *faceIndex = bcRegion->GetFaceIndex();
        Data_Param *bcData = bcRegion->GetBCParamDataBase();

        if (bcType < 0)
        {
            interfaceIndex += 1;
            continue;
        }

        //! Update the pressure boundary conditions.
        static_cast<PBBoundaryFlow *>(boundCalc[iBCRegionUnstruct - interfaceIndex])->UpdatePressureBoundaryCondition(grid, faceIndex, bcData);

        //! Update the velocity boundary conditions.
        static_cast<PBBoundaryFlow *>(boundCalc[iBCRegionUnstruct - interfaceIndex])->UpdateVelocityBoundaryCondition(grid, faceIndex, bcData);
    }
}

void PBSolverFlow::UpdateCorrectionBoundaryPressureAndVelocity()
{
    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegionUnstruct = unstructBCSet->GetnBCRegion();
    int interfaceIndex = 0;

    for (int iBCRegionUnstruct = 0; iBCRegionUnstruct < nBCRegionUnstruct; iBCRegionUnstruct++)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegionUnstruct);
        int bcType = bcRegion->GetBCType();
        vector<int> *faceIndex = bcRegion->GetFaceIndex();
        Data_Param *bcData = bcRegion->GetBCParamDataBase();

        if (bcType < 0)
        {
            interfaceIndex += 1;
            continue;
        }

        //! Update and correct the pressure boundary conditions.
        static_cast<PBBoundaryFlow *>(boundCalc[iBCRegionUnstruct - interfaceIndex])->UpdatePressureBoundaryCondition(grid, faceIndex, bcData);

        //! Update and correct the velocity boundary conditions.
        static_cast<PBBoundaryFlow *>(boundCalc[iBCRegionUnstruct - interfaceIndex])->UpdateVelocityCorrectionBoundary(grid, faceIndex, bcData);
    }
}

void PBSolverFlow::CalcPressureGrad()
{
    gradCalc["P"]->Calc(grid, "P","dpdx", "dpdy", "dpdz");

    if (bodyForceFlag)
    {
        int nTotalCell = grid->GetNTotalCell();
        int nBoundFace = grid->GetNBoundFace();
        int nTotalFace = grid->GetNTotalFace();
        int *leftCellOfFace = grid->GetLeftCellOfFace();
        int *rightCellOfFace = grid->GetRightCellOfFace();
        RDouble *faceArea    = grid->GetFaceArea();
        RDouble *xfn = grid->GetFaceNormalX();
        RDouble *yfn = grid->GetFaceNormalY();
        RDouble *zfn = grid->GetFaceNormalZ();
        RDouble *xfc = grid->GetFaceCenterX();
        RDouble *yfc = grid->GetFaceCenterY();
        RDouble *zfc = grid->GetFaceCenterZ();
        RDouble *xcc = grid->GetCellCenterX();
        RDouble *ycc = grid->GetCellCenterY();
        RDouble *zcc = grid->GetCellCenterZ();
        RDouble* dbfdx = reinterpret_cast<RDouble *>(grid->GetDataPtr("dbfdx"));
        RDouble* dbfdy = reinterpret_cast<RDouble *>(grid->GetDataPtr("dbfdy"));
        RDouble* dbfdz = reinterpret_cast<RDouble *>(grid->GetDataPtr("dbfdz"));

        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            dbfdx[iCell] = 0.0;
            dbfdy[iCell] = 0.0;
            dbfdz[iCell] = 0.0;
        }

        for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
        {
            int le = leftCellOfFace[iFace];
            int re = rightCellOfFace[iFace];
            RDouble area = faceArea[iFace];
            RDouble ax = area * xfn[iFace];
            RDouble ay = area * yfn[iFace];
            RDouble az = area * zfn[iFace];
            RDouble dxfL = xfc[iFace] - xcc[le];
            RDouble dyfL = yfc[iFace] - ycc[le];
            RDouble dzfL = zfc[iFace] - zcc[le];
            RDouble dxfR = xfc[iFace] - xcc[re];
            RDouble dyfR = yfc[iFace] - ycc[re];
            RDouble dzfR = zfc[iFace] - zcc[re];
            RDouble phif = deltL[iFace] * (dxfL * bfU[le] + dyfL * bfV[le] + dzfL * bfW[le]) +
                deltR[iFace] * (dxfR * bfU[re] + dyfR * bfV[re] + dzfR * bfW[re]);
            dbfdx[le] += ax * phif;
            dbfdy[le] += ay * phif;
            dbfdz[le] += az * phif;
            dbfdx[re] -= ax * phif;
            dbfdy[re] -= ay * phif;
            dbfdz[re] -= az * phif;
        }

        for (int ic = 0; ic < nTotalCell; ++ ic)
        {
            RDouble cv = grid->GetCellVolume()[ic];
            dpdx[ic] += dbfdx[ic] / cv;
            dpdy[ic] += dbfdy[ic] / cv;
            dpdz[ic] += dbfdz[ic] / cv;
        }
    }
}

void PBSolverFlow::CalcVelocityGrad()
{
    gradCalc["U"]->Calc(grid, "U", "dudx", "dudy", "dudz");
    gradCalc["V"]->Calc(grid, "V", "dvdx", "dvdy", "dvdz");
    gradCalc["W"]->Calc(grid, "W", "dwdx", "dwdy", "dwdz");
}

void PBSolverFlow::CalcVelocityBoundaryConditions()
{
    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegionUnstruct = unstructBCSet->GetnBCRegion();
    int interfaceIndex = 0;

    for (int iBCRegionUnstruct = 0; iBCRegionUnstruct < nBCRegionUnstruct; iBCRegionUnstruct++)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegionUnstruct);
        int bcType = bcRegion->GetBCType();
        vector<int> *faceIndex = bcRegion->GetFaceIndex();

        if (bcType < 0)
        {
            interfaceIndex += 1;
            continue;
        }

        static_cast<PBBoundaryFlow *>(boundCalc[iBCRegionUnstruct - interfaceIndex])->CalcVelocityBoundaryCondition(grid, faceIndex);
    }
}

void PBSolverFlow::CalcVelocityTransientTerm()
{
    tranCalc->Discrt(grid, "U", "UOld", "UOldOld", "rho", "rho");
}

void PBSolverFlow::CalcVelocityConvectionTerm()
{
    convCalc->Discrt(grid, "U", "apU", "scU", v, apV, scV, w, apW, scW);
}

void PBSolverFlow::CalcVelocityDiffusionTerm()
{
    diffCalc->Discrt(grid, "U", mu, dudx, dudy, dudz, ajLeft, ajRight, scU, dvdx, dvdy, dvdz, scV, dwdx, dwdy, dwdz, scW);
}

void PBSolverFlow::CalcVelocitySourceTerm()
{
    for (vector<PBSource *>::iterator iter = sourceCalc.begin(); iter != sourceCalc.end(); ++ iter)
    {
        (*iter)->Discrt(grid, "U");
    }
}

void PBSolverFlow::PrimeVelocityEquations()
{
    for (int iFace = grid->GetNBoundFace(); iFace < grid->GetNTotalFace(); ++ iFace)
    {
        int iFacelocal = iFace - grid->GetNBoundFace();
        int le = grid->GetLeftCellOfFace()[iFace];
        int re = grid->GetRightCellOfFace()[iFace];
        apU[le] += ajLeft[iFacelocal];
        apV[le] += ajLeft[iFacelocal];
        apW[le] += ajLeft[iFacelocal];
        apU[re] += ajRight[iFacelocal];
        apV[re] += ajRight[iFacelocal];
        apW[re] += ajRight[iFacelocal];
        scU[le] += ajLeft[iFacelocal] * u[re];
        scV[le] += ajLeft[iFacelocal] * v[re];
        scW[le] += ajLeft[iFacelocal] * w[re];
        scU[re] += ajRight[iFacelocal] * u[le];
        scV[re] += ajRight[iFacelocal] * v[le];
        scW[re] += ajRight[iFacelocal] * w[le];
    }

    for (int iCell = 0; iCell < grid->GetNTotalCell(); ++ iCell)
    {
        scU[iCell] -= apU[iCell] * u[iCell];
        scV[iCell] -= apV[iCell] * v[iCell];
        scW[iCell] -= apW[iCell] * w[iCell];
    }
}

void PBSolverFlow::CalcVelocityRes()
{
    uResNow = 0.0;
    vResNow = 0.0;
    wResNow = 0.0;

    for (int iCell = 0; iCell < grid->GetNTotalCell(); ++ iCell)
    {
        uResNow += std::abs(scU[iCell]);
        vResNow += std::abs(scV[iCell]);
        wResNow += std::abs(scW[iCell]);
    }

    grid->UpdateData("uResNow", &uResNow, PHDOUBLE, 1);
    grid->UpdateData("vResNow", &vResNow, PHDOUBLE, 1);
    grid->UpdateData("wResNow", &wResNow, PHDOUBLE, 1);
}

void PBSolverFlow::RelaxVelocity()
{
    for (int iCell = 0; iCell < grid->GetNTotalCell(); ++ iCell)
    {
        apU[iCell] = (1.0 + urfU) * apU[iCell];
        apV[iCell] = (1.0 + urfV) * apV[iCell];
        apW[iCell] = (1.0 + urfW) * apW[iCell];
    }
}

void PBSolverFlow::CoefPrepForPres()
{
    for (int iCell = 0; iCell < grid->GetNTotalCell(); ++ iCell)
    {
        hu[iCell] = apU[iCell];
        hv[iCell] = apV[iCell];
        hw[iCell] = apW[iCell];
    }

    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotalFace = grid->GetNTotalFace();
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        int iFacelocal = iFace - grid->GetNBoundFace();
        int le = leftCellOfFace[iFace];
        int re = rightCellOfFace[iFace];
        hu[le] -= ajLeft[iFacelocal];
        hu[re] -= ajRight[iFacelocal];
        hv[le] -= ajLeft[iFacelocal];
        hv[re] -= ajRight[iFacelocal];
        hw[le] -= ajLeft[iFacelocal];
        hw[re] -= ajRight[iFacelocal];
    }

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        RDouble vol = grid->GetCellVolume()[iCell];
        hu[iCell] = vol / hu[iCell];
        hv[iCell] = vol / hv[iCell];
        hw[iCell] = vol / hw[iCell];
        apU[iCell] = vol / apU[iCell];
        apV[iCell] = vol / apV[iCell];
        apW[iCell] = vol / apW[iCell];
    }
}

void PBSolverFlow::PrimePressureEquations()
{
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotalFace = grid->GetNTotalFace();
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    RDouble *area    = grid->GetFaceArea();
    RDouble *xfn = grid->GetFaceNormalX();
    RDouble *yfn = grid->GetFaceNormalY();
    RDouble *zfn = grid->GetFaceNormalZ();
    RDouble *xfc = grid->GetFaceCenterX();
    RDouble *yfc = grid->GetFaceCenterY();
    RDouble *zfc = grid->GetFaceCenterZ();
    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();
    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        int le = leftCellOfFace[iFace];
        int re = rightCellOfFace[iFace];

        RDouble rurf = urfFlux / (1.0 + urfFlux);

        RDouble faceArea = area[iFace];
        RDouble faceNormalX = xfn[iFace];
        RDouble faceNormalY = yfn[iFace];
        RDouble faceNormalZ = zfn[iFace];


        RDouble ax = faceArea * faceNormalX;
        RDouble ay = faceArea * faceNormalY;
        RDouble az = faceArea * faceNormalZ;
        RDouble cellCenterXre = xcc[re];
        RDouble cellCenterXle = xcc[le];
        RDouble cellCenterYre = ycc[re];
        RDouble cellCenterYle = ycc[le];
        RDouble cellCenterZre = zcc[re];
        RDouble cellCenterZle = zcc[le];

        RDouble dx = cellCenterXre - cellCenterXle;
        RDouble dy = cellCenterYre - cellCenterYle;
        RDouble dz = cellCenterZre - cellCenterZle;

        RDouble wfL = deltL[iFace];
        RDouble wfR = deltR[iFace];
        RDouble rhof = this->rhof[iFace];
        RDouble ad = ax * dx + ay * dy + az * dz;

        RDouble ubar_a = (wfL * u[le] + wfR * u[re]) * ax + (wfL * v[le] + wfR * v[re]) * ay + (wfL * w[le] + wfR * w[re]) * az;

        RDouble pDiff;
        RDouble gDiff;
        RDouble gddf;
        if (bodyForceFlag != 1)
        {
            pDiff = p[re] - p[le] - 0.5 * ((dpdx[le] + dpdx[re]) * dx + (dpdy[le] + dpdy[re]) * dy + (dpdz[le] + dpdz[re]) * dz);
        }
        else
        {
            RDouble cellCenterXre = xcc[re];
            RDouble cellCenterXle = xcc[le];
            RDouble cellCenterYre = ycc[re];
            RDouble cellCenterYle = ycc[le];
            RDouble cellCenterZre = zcc[re];
            RDouble cellCenterZle = zcc[le];
            RDouble faceCenterX = xfc[iFace];
            RDouble faceCenterY = yfc[iFace];
            RDouble faceCenterZ = zfc[iFace];

            RDouble dxfL = faceCenterX - cellCenterXle;
            RDouble dyfL = faceCenterY - cellCenterYle;
            RDouble dzfL = faceCenterZ - cellCenterZle;
            RDouble dxfR = faceCenterX - cellCenterXre;
            RDouble dyfR = faceCenterY - cellCenterYre;
            RDouble dzfR = faceCenterZ - cellCenterZre;

            pDiff = p[re] + (bfU[re] * dxfR + bfV[re] * dyfR + bfW[re] * dzfR) - (p[le] + (bfU[le] * dxfL + bfV[le] * dyfL + bfW[le] * dzfL))
                - 0.5 * (((dpdx[le] - bfU[le]) + (dpdx[re] - bfU[re])) * dx + ((dpdy[le] - bfV[le]) + (dpdy[re] - bfV[re])) * dy
                    + ((dpdz[le] - bfW[le]) + (dpdz[re] - bfW[re])) * dz);
        };

        pDiff /= ad;
        RDouble ddf = (apU[le] * wfL + apU[re] * wfR) * ax * ax + (apV[le] * wfL + apV[re] * wfR) * ay * ay + (apW[le] * wfL + apW[re] * wfR) * az * az;
        RDouble flux = rhof * (ubar_a - ddf * pDiff + rurf * dun[iFace]);
        RDouble aj = rhof * ((hu[le] * wfL + hu[re] * wfR) * ax * ax + (hv[le] * wfL + hv[re] * wfR) * ay * ay + (hw[le] * wfL + hw[re] * wfR) * az * az) / ad;
        int iFacelocal = iFace - nBoundFace;

        this->flux[iFace] = flux;
        ajLeft[iFacelocal] = aj;
        ajRight[iFacelocal] = aj;
        ap[le] += aj;
        ap[re] += aj;
        sc[le] -= flux;
        sc[re] += flux;
    }

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegionUnstruct = unstructBCSet->GetnBCRegion();
    int interfaceIndex = 0;

    for (int iBCRegionUnstruct = 0; iBCRegionUnstruct < nBCRegionUnstruct; iBCRegionUnstruct++)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegionUnstruct);
        int bcType = bcRegion->GetBCType();
        vector<int> *faceIndex = bcRegion->GetFaceIndex();
        Data_Param *bcData = bcRegion->GetBCParamDataBase();

        if (bcType < 0)
        {
            interfaceIndex += 1;
            continue;
        }

        static_cast<PBBoundaryFlow *>(boundCalc[iBCRegionUnstruct - interfaceIndex])->CalcPressureCorrectionBoundaryCondition(grid, faceIndex, bcData);
    }
}

void PBSolverFlow::CalcPressureRes()
{
    pResNow = 0.0;

    for (int iCell = 0; iCell < grid->GetNTotalCell(); ++ iCell)
    {
        pResNow += std::abs(sc[iCell]);
    }

    grid->UpdateData("pResNow", &pResNow, PHDOUBLE, 1);
}

void PBSolverFlow::CalcPressureCorrectionGrad()
{
    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegionUnstruct = unstructBCSet->GetnBCRegion();
    int interfaceIndex = 0;

    for (int iBCRegionUnstruct = 0; iBCRegionUnstruct < nBCRegionUnstruct; iBCRegionUnstruct++)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegionUnstruct);
        int bcType = bcRegion->GetBCType();
        vector<int> *faceIndex = bcRegion->GetFaceIndex();
        Data_Param *bcData = bcRegion->GetBCParamDataBase();

        if (bcType < 0)
        {
            interfaceIndex += 1;
            continue;
        }

        static_cast<PBBoundaryFlow *>(boundCalc[iBCRegionUnstruct - interfaceIndex])->UpdatePressureCorrectionBoundaryCondition(grid, faceIndex,bcData);

    }

    gradCalc["P"]->Calc(grid, "pc", "dpcdx", "dpcdy", "dpcdz");
}

void PBSolverFlow::UpdateVelocityAndFlux()
{
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotalFace = grid->GetNTotalFace();
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    RDouble *area    = grid->GetFaceArea();
    RDouble *xfn = grid->GetFaceNormalX();
    RDouble *yfn = grid->GetFaceNormalY();
    RDouble *zfn = grid->GetFaceNormalZ();
    RDouble *xfc = grid->GetFaceCenterX();
    RDouble *yfc = grid->GetFaceCenterY();
    RDouble *zfc = grid->GetFaceCenterZ();
    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        u[iCell] -= hu[iCell] * dpcdx[iCell];
        v[iCell] -= hv[iCell] * dpcdy[iCell];
        w[iCell] -= hw[iCell] * dpcdz[iCell];
    }

    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        int le = leftCellOfFace[iFace];
        int re = rightCellOfFace[iFace];
        RDouble wfL = deltL[iFace];
        RDouble wfR = deltR[iFace];

        RDouble faceArea = area[iFace];
        RDouble faceNormalX = xfn[iFace];
        RDouble faceNormalY = yfn[iFace];
        RDouble faceNormalZ = zfn[iFace];

        RDouble ax = faceArea * faceNormalX;
        RDouble ay = faceArea * faceNormalY;
        RDouble az = faceArea * faceNormalZ;
        RDouble uf = wfL * u[le] + wfR * u[re];
        RDouble vf = wfL * v[le] + wfR * v[re];
        RDouble wf = wfL * w[le] + wfR * w[re];
        RDouble un = uf * ax + vf * ay + wf * az;

        int iFacelocal = iFace - nBoundFace;
        flux[iFace] -= ajLeft[iFacelocal] * (pc[re] - pc[le]);
        dun[iFace] = flux[iFace] / rhof[iFace] - un;
    }
}

void PBSolverFlow::CheckFaceCells()
{
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotalFace = grid->GetNTotalFace();
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        int le = leftCellOfFace[iFace];
        int re = rightCellOfFace[iFace];

        if (le > re)
        {
            int tmp = le;
            le = re;
            re = tmp;
        }
    }

    cout << "CheckFaceCells" << endl;
}

DEFINE_PBSOLVER_REGISTER_FUNCTION(PBSolver, PBSolverFlow, SOLVER_FLOW);
}