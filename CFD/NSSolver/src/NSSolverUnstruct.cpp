#include "NSSolverUnstruct.h"
#include "MultiGridOperation.h"
#include "GradientOperation.h"
#include "Math_Limiter.h"
#include "AerodynamicsBasicRelations.h"
#include "Glb_Dimension.h"
#include "Constants.h"
#include <list>
#include "Param_NSSolverUnstruct.h"
#include "Post_ForceMoment.h"
#include "Residual.h"
#include "Jacobian.h"
#include "FieldProxy.h"
#include "TK_Exit.h"
#include "TK_Log.h"
#include "IO_FileName.h"
#include "PHIO.h"
#include "Force.h"
#include "AleForceManager.h"
#include <cmath>
#include "Flux_RoeEntropyCorrection.h"
#include "Gas.h"

using namespace std;

#pragma warning (disable:913)
namespace PHSPACE
{
using namespace GAS_SPACE;
using namespace IDX;

NSSolverUnstruct::NSSolverUnstruct()
{
    gradientNSField = 0;
    gradientTemperature = 0;
}

NSSolverUnstruct::~NSSolverUnstruct()
{
    DeAllocateGlobalVariables();
    FreeControlParameters();
}

void NSSolverUnstruct::ReadParameter()
{
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int viscousType = parameters->GetViscousType();

    string uns_scheme_name = parameters->GetSpaceSecheme();

    int uns_scheme = GetSchemeID(uns_scheme_name);
    GlobalDataBase::UpdateData("uns_scheme", &uns_scheme, PHINT, 1);

    if (uns_scheme == ISCHEME_ROE)
    {
        parameters->SetEntropyFixCoefficients();

        Grid *grid = GetGrid();
        int zoneID = grid->GetZoneID();
        RDouble entropyAcoustic = parameters->GetRoeEntropyFixCoef1();
        
        //! The output of Entropy fix coefficient for ROE scheme.
        if (zoneID == 0)  //! When grid->GetZoneID() = 0 and GetCurrentProcessorID() = 0.
        {
            PrintToWindow("Entropy fix coefficient for ROE scheme: ", entropyAcoustic, "\n");
        }
    }

    string uns_limiter_name;
    GlobalDataBase::GetData("uns_limiter_name", &uns_limiter_name, PHSTRING, 1);

    int limiterType = GetLimiterID(uns_limiter_name);
    parameters->SetLimiterType(limiterType);
    GlobalDataBase::UpdateData("uns_limiter", &limiterType, PHINT, 1);

    if (viscousType == INVISCID) return;

    string uns_vis_name = "default";
    GlobalDataBase::GetData("uns_vis_name", &uns_vis_name, PHSTRING, 1);

    int uns_vis;
    if (uns_vis_name.substr(0,4) == "test")
    {
        uns_vis = VIS_TEST;
    }
    else if (uns_vis_name.substr(0,4) == "aver")
    {
        uns_vis = VIS_AVER;
    }
    else if (uns_vis_name.substr(0,3) == "std")
    {
        uns_vis = VIS_STD;
    }
    else if (uns_vis_name.substr(0,4) == "new1")
    {
        uns_vis = VIS_NEW1;
    }
    else if (uns_vis_name.substr(0,4) == "new2")
    {
        uns_vis = VIS_NEW2;
    }
    else if (uns_vis_name.substr(0,4) == "eric")
    {
        uns_vis = VIS_ERIC;
    }
    else
    {
        uns_vis = VIS_STD;
    }

    GlobalDataBase::UpdateData("uns_vis", &uns_vis, PHINT, 1);
}

void NSSolverUnstruct::AllocateGlobalVar(Grid *gridIn)
{
    CFDSolver::AllocateGlobalVar(gridIn);

    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nNSEquation = parameters->GetNSEquationNumber();

    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    RDouble **q   = NewPointer2<RDouble>(nEquation, nTotal);
    RDouble **res = NewPointer2<RDouble>(nEquation, nTotal);
    RDouble *diagonal = new RDouble [nTotal];

    grid->UpdateDataPtr("q", q);    //! Scalar flow field variable(rho/u/v/w/p).
    grid->UpdateDataPtr("res", res);    //! Residual or right-hand side.
    grid->UpdateDataPtr("diagonal", diagonal);    //! Sacrificing space to improve efficiency.
    
    int nTotalNode = grid->GetNTotalNode();
    RDouble **qnode = NewPointer2<RDouble>(nEquation, nTotalNode);
    grid->UpdateDataPtr("qnode", qnode);    //! Scalar flow field variable(rho/u/v/w/p) at node.

    RDouble **tnode = NewPointer2<RDouble>(nTemperatureModel, nTotalNode);
    grid->UpdateDataPtr("tnode", tnode);    //! Static temperature at node.

    RDouble *nodeWeight = new RDouble[nTotalNode];
    grid->UpdateDataPtr("nodeWeight", nodeWeight);    //! Weight coefficient at node.

    RDouble *nodeWeightTrade = new RDouble[nTotalNode];
    PHSPACE::SetField(nodeWeightTrade, 0.0, nTotalNode);
    grid->UpdateDataPtr("nodeWeightTrade", nodeWeightTrade);

    if (grid->GetLevel() == 0)
    {
        int *nodeBCType = new int[nTotalNode];
        PHSPACE::SetField(nodeBCType, 0, nTotalNode);
        grid->UpdateDataPtr("nodeBCType", nodeBCType);
    }

    InterpointInformation *interpointInformation = grid->GetInterpointInfo();
    if (interpointInformation)
    {
        int numberOfInterpoints = interpointInformation->GetNumberOfInterpoints();
        RDouble **qInterPoint = NewPointer2<RDouble>(nEquation, numberOfInterpoints);
        grid->UpdateDataPtr("qInterPoint", qInterPoint);    //! Scalar flow field variable(rho/u/v/w/p) at interPoint.

        //RDouble **tInterPoint = NewPointer2<RDouble>(1, numberOfInterpoints);
        RDouble **tInterPoint = NewPointer2<RDouble>(nTemperatureModel, numberOfInterpoints);
        grid->UpdateDataPtr("tInterPoint", tInterPoint);    //! Static temperature at interPoint.
    }

    RDouble *rtem = new RDouble [nTotal];
    grid->UpdateDataPtr("rtem", rtem);    //! Pressure factor.

    RDouble *invSpectralRadius = new RDouble [nTotalCell];
    grid->UpdateDataPtr("invSpectralRadius", invSpectralRadius);    //! Inviscid spectral radius.

    RDouble *visSpectralRadius = new RDouble [nTotalCell];
    grid->UpdateDataPtr("visSpectralRadius", visSpectralRadius);    //! Viscous spectral radius.

    int nSpeciesNumber = parameters->GetNumberOfSpecies();
    int nChemicalRadius = parameters->GetNChemicalRadius();

    if (nChemical == 1 && nChemicalRadius == 1)
    {
        RDouble **chemSpectralRadius = NewPointer2<RDouble>(nSpeciesNumber,nTotalCell);
        grid->UpdateDataPtr("chemSpectralRadius"  , chemSpectralRadius);    //! Chemical reaction source item.
    }

    int isUnsteady = parameters->GetIsUnsteady();

    if (isUnsteady)
    {
        RDouble **q_unsteady_n1   = NewPointer2<RDouble>(nEquation, nTotal);
        RDouble **q_unsteady_n2   = NewPointer2<RDouble>(nEquation, nTotal);
        RDouble **res_unsteady_n1 = NewPointer2<RDouble>(nEquation, nTotal);
        RDouble **res_unsteady_n2 = NewPointer2<RDouble>(nEquation, nTotal);

        //! Store the historical sum of InviscidFlux, ViscousFlux, et al., exclusive of DualTimeSourceFlux, used in UpdateUnsteadyFlow for the Rn in the dualtime method.
        RDouble **res_unsteady_tmp = NewPointer2<RDouble>(nEquation, nTotal);

        grid->UpdateDataPtr("q_unsteady_n1"  , q_unsteady_n1 );
        grid->UpdateDataPtr("q_unsteady_n2"  , q_unsteady_n2 );
        grid->UpdateDataPtr("res_unsteady_n1", res_unsteady_n1);
        grid->UpdateDataPtr("res_unsteady_n2", res_unsteady_n2);
        grid->UpdateDataPtr("res_unsteady_tmp", res_unsteady_tmp);

        //! Statistical variables for unsteady simulation.
        RDouble **qAverage = NewPointer2<RDouble>(nEquation, nTotal); 
        grid->UpdateDataPtr("qAverage", qAverage);
        PHSPACE::SetField(qAverage, nEquation, nTotal, 0.0);

        RDouble **tauAverage = NewPointer2<RDouble>(6, nTotal); 
        grid->UpdateDataPtr("tauAverage", tauAverage);
        PHSPACE::SetField(tauAverage, 6, nTotal, 0.0);

        RDouble **q2Average = NewPointer2<RDouble>(6, nTotal); 
        grid->UpdateDataPtr("q2Average", q2Average);
        PHSPACE::SetField(q2Average, 6, nTotal, 0.0);
    }

    RDouble **t = NewPointer2<RDouble>(nTemperatureModel, nTotal);
    grid->UpdateDataPtr("t", t);    //! Static temperature.

    RDouble *dt  = new RDouble [nTotal];
    grid->UpdateDataPtr("dt", dt);    //! Time step.

    RDouble *gama  = new RDouble [nTotal];
    grid->UpdateDataPtr("gama", gama);    //! Ratio of specific heat coefficients at constant pressure and volume.

    int isViscous = parameters->IsViscous();

    if (isViscous)
    {
        RDouble *visl = new RDouble [nTotal];
        RDouble *vist = new RDouble [nTotal];

        PHSPACE::SetField(visl, 1.0, nTotal);
        PHSPACE::SetField(vist, 0.0, nTotal);

        grid->UpdateDataPtr("visl", visl);    //! Laminar viscous coefficient.
        grid->UpdateDataPtr("vist", vist);    //! Turbulence viscous coefficient.
    }

    int iLES = GlobalDataBase::GetIntParaFromDB("iLES");
    if (iLES == LES_SOLVER)
    {
        RDouble *subgridScaleEnergy = new RDouble[nTotal];
        PHSPACE::SetField(subgridScaleEnergy, 0.0, nTotal);
        grid->UpdateDataPtr("subgridScaleEnergy", subgridScaleEnergy);
    }

    RDouble **gradPrimtiveVarX = NewPointer2<RDouble>(nEquation, nTotal);
    RDouble **gradPrimtiveVarY = NewPointer2<RDouble>(nEquation, nTotal);
    RDouble **gradPrimtiveVarZ = NewPointer2<RDouble>(nEquation, nTotal);
    PHSPACE::SetField(gradPrimtiveVarX, nEquation, nTotal, 0.0);
    PHSPACE::SetField(gradPrimtiveVarY, nEquation, nTotal, 0.0);
    PHSPACE::SetField(gradPrimtiveVarZ, nEquation, nTotal, 0.0);

    grid->UpdateDataPtr("gradPrimtiveVarX", gradPrimtiveVarX);    //! Gradient of scalar flow field variable at cell, for x direction.
    grid->UpdateDataPtr("gradPrimtiveVarY", gradPrimtiveVarY);    //! Gradient of scalar flow field variable at cell, for y direction.
    grid->UpdateDataPtr("gradPrimtiveVarZ", gradPrimtiveVarZ);    //! Gradient of scalar flow field variable at cell, for z direction.

    RDouble **gradTemperatureX = NewPointer2<RDouble>(nTemperatureModel, nTotal);
    RDouble **gradTemperatureY = NewPointer2<RDouble>(nTemperatureModel, nTotal);
    RDouble **gradTemperatureZ = NewPointer2<RDouble>(nTemperatureModel, nTotal);
    PHSPACE::SetField(gradTemperatureX, nTemperatureModel, nTotal, 0.0);
    PHSPACE::SetField(gradTemperatureY, nTemperatureModel, nTotal, 0.0);
    PHSPACE::SetField(gradTemperatureZ, nTemperatureModel, nTotal, 0.0);
    grid->UpdateDataPtr("gradTemperatureX", gradTemperatureX);    //! Gradient of static temperature at cell, for x direction.
    grid->UpdateDataPtr("gradTemperatureY", gradTemperatureY);    //! Gradient of static temperature at cell, for y direction.
    grid->UpdateDataPtr("gradTemperatureZ", gradTemperatureZ);    //! Gradient of static temperature at cell, for z direction.

    systemgridtype = GlobalDataBase::GetIntParaFromDB("sys_gridtype");
    if (MIXGRID == systemgridtype)
    {
        InitiallimiterOnlyForMixGrid(grid);

        InitialInviscidfaceproxyOnlyForMixGrid(grid);

        if (isViscous)
        {
            InitialViscousfaceproxyOnlyForMixGrid(grid);
        }
    }
}

void NSSolverUnstruct::DeAllocateGlobalVar(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nChemical = parameters->GetChemicalFlag();

    int nChemicalRadius = parameters->GetNChemicalRadius();

    RDouble **q    = reinterpret_cast <RDouble **> (grid->GetDataPtr("q"));
    RDouble **res  = reinterpret_cast <RDouble **> (grid->GetDataPtr("res"));
    RDouble *diagonal = reinterpret_cast <RDouble *> (grid->GetDataPtr("diagonal"));

    DelPointer2(q);
    DelPointer2(res);
    delete [] diagonal;

    if (grid->GetLevel() == 0)
    {
        int *nodeBCType = reinterpret_cast <int *> (grid->GetDataPtr("nodeBCType"));
        delete [] nodeBCType;
    }

    RDouble *rtem = reinterpret_cast <RDouble *> (grid->GetDataPtr("rtem"));
    delete [] rtem;
    RDouble *invSpectralRadius = reinterpret_cast <RDouble *> (grid->GetDataPtr("invSpectralRadius"));
    delete [] invSpectralRadius;
    RDouble *visSpectralRadius = reinterpret_cast <RDouble *> (grid->GetDataPtr("visSpectralRadius"));
    delete [] visSpectralRadius;

    RDouble* DiagValue = reinterpret_cast <RDouble*> (grid->GetDataPtr("DiagValue"));
    if (DiagValue != NULL)
    {
        delete [] DiagValue;    DiagValue = NULL;
    }

    RDouble**** DMatrixStartByLine = reinterpret_cast<RDouble****>(grid->GetDataPtr("DMatrixStartByLine"));
    RDouble**** UMatrixStartByLine = reinterpret_cast<RDouble****>(grid->GetDataPtr("UMatrixStartByLine"));
    RDouble**** LMatrixStartByLine = reinterpret_cast<RDouble****>(grid->GetDataPtr("LMatrixStartByLine"));

    if (DMatrixStartByLine != NULL)
    {
        delete [] DMatrixStartByLine;    DMatrixStartByLine = NULL;
        delete [] UMatrixStartByLine;    UMatrixStartByLine = NULL;
        delete [] LMatrixStartByLine;    LMatrixStartByLine = NULL;
    }

    MAT3D< RDouble >* diagMatrix = reinterpret_cast<MAT3D< RDouble > *>(grid->GetDataPtr("diagMatrix"));
    MAT3D< RDouble >* lowerMatrix = reinterpret_cast<MAT3D< RDouble > *>(grid->GetDataPtr("lowerMatrix"));
    MAT3D< RDouble >* upperMatrix = reinterpret_cast<MAT3D< RDouble > *>(grid->GetDataPtr("upperMatrix"));

    if (diagMatrix != NULL)
    {
        delete diagMatrix;
        delete upperMatrix;
        delete lowerMatrix;
    }

    if (nChemical == 1 && nChemicalRadius == 1)
    {
        RDouble **chemSpectralRadius = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("chemSpectralRadius") );
        DelPointer2(chemSpectralRadius);
    }

    int isUnsteady = parameters->GetIsUnsteady();
    if (isUnsteady)
    {
        RDouble **q_unsteady_n1   = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q_unsteady_n1" ));
        RDouble **q_unsteady_n2   = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q_unsteady_n2" ));
        RDouble **res_unsteady_n1 = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res_unsteady_n1"));
        RDouble **res_unsteady_n2 = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res_unsteady_n2"));

        //! Store the historical sum of InviscidFlux, ViscousFlux, et al., exclusive of DualTimeSourceFlux, used in UpdateUnsteadyFlow for the Rn in the dualtime method.
        RDouble **res_unsteady_tmp = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res_unsteady_tmp"));

        DelPointer2(q_unsteady_n1 );
        DelPointer2(q_unsteady_n2 );
        DelPointer2(res_unsteady_n1);
        DelPointer2(res_unsteady_n2);
        DelPointer2(res_unsteady_tmp);

        //! Statistical variables for unsteady simulation.
        int ifStaticsFlowField = parameters->GetIfStaticsFlowField();
        RDouble **qAverage = reinterpret_cast< RDouble ** > (grid->GetDataPtr("qAverage"));
        DelPointer2(qAverage);

        //! Statistical variables for unsteady simulation.
        RDouble **tauAverage = reinterpret_cast< RDouble ** > (grid->GetDataPtr("tauAverage"));
        DelPointer2(tauAverage);

        RDouble **q2Average = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q2Average"));
        DelPointer2(q2Average);
    }

    RDouble **t   = reinterpret_cast< RDouble ** > (grid->GetDataPtr("t" ));
    DelPointer2(t );

    RDouble * dt  = reinterpret_cast< RDouble * > (grid->GetDataPtr("dt" ));
    delete [] dt;    dt = NULL;

    RDouble * gama  = reinterpret_cast< RDouble * > (grid->GetDataPtr("gama" ));
    delete [] gama;    gama = NULL;

    int viscousType = parameters->GetViscousType();

    if (viscousType > INVISCID)
    {
        RDouble * visl = reinterpret_cast< RDouble * > (grid->GetDataPtr("visl"));
        RDouble * vist = reinterpret_cast< RDouble * > (grid->GetDataPtr("vist"));
        delete [] visl;    visl = NULL;
        delete [] vist;    vist = NULL;
    }

    int iLES = GlobalDataBase::GetIntParaFromDB("iLES");
    if(iLES == LES_SOLVER)
    {
        RDouble * subgridScaleEnergy = reinterpret_cast< RDouble * > (grid->GetDataPtr("subgridScaleEnergy"));
        delete [] subgridScaleEnergy;
    }

    RDouble **gradPrimtiveVarX = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarX"));
    RDouble **gradPrimtiveVarY = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarY"));
    RDouble **gradPrimtiveVarZ = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarZ"));

    DelPointer2(gradPrimtiveVarX);
    DelPointer2(gradPrimtiveVarY);
    DelPointer2(gradPrimtiveVarZ);

    RDouble **gradTemperatureX = reinterpret_cast< RDouble ** > (grid->GetDataPtr("gradTemperatureX" ));
    RDouble **gradTemperatureY = reinterpret_cast< RDouble ** > (grid->GetDataPtr("gradTemperatureY" ));
    RDouble **gradTemperatureZ = reinterpret_cast< RDouble ** > (grid->GetDataPtr("gradTemperatureZ" ));

    DelPointer2(gradTemperatureX);
    DelPointer2(gradTemperatureY);
    DelPointer2(gradTemperatureZ);
}

bool NSSolverUnstruct::JudgeIfRestart()
{
    string restartNSFile = ".\results\flow.dat";
    GlobalDataBase::GetData("restartNSFile", &restartNSFile, PHSTRING, 1);

    if(PHMPI::IsParallelRun())
    {
        restartNSFile = PHSPACE::AddSymbolToFileName(restartNSFile, "_", 0);
    }

    bool restart_flag = false;

    ifstream infile(restartNSFile.c_str(), ios::in);
    if (infile)
    {
        restart_flag = true;
        infile.close();
        infile.clear();
    }

    return restart_flag;
}

bool NSSolverUnstruct::JudgeIfReadAverage()
{
    string restartNSFile = ".\results\flow.dat";
    GlobalDataBase::GetData("restartNSFile", &restartNSFile, PHSTRING, 1);
    string averageFlowFile = AddSymbolToFileName(restartNSFile, "_", "Average");

    if (PHMPI::IsParallelRun())
    {
        averageFlowFile = PHSPACE::AddSymbolToFileName(averageFlowFile, "_", 0);
    }

    bool restart_flag = false;

    ifstream infile(averageFlowFile.c_str(), ios::in);
    if (infile)
    {
        restart_flag = true;
        infile.close();
        infile.clear();
    }

    return restart_flag;
}

void NSSolverUnstruct::DumpRestart(ActionKey *actkey)
{
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int level = actkey->level;
    UnstructGrid * grid = UnstructGridCast(GetGrid(level));

    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    int outnstep = 0;
    GlobalDataBase::GetData("outnstep", &outnstep, PHINT, 1);

    int nBoundFace = grid->GetNBoundFace();
    int nTotalCell = grid->GetNTotalCell();
    int nTotal = nTotalCell + nBoundFace;

    DataContainer * cdata = actkey->GetData();
    cdata->moveto_begin();
    cdata->write(&outnstep,sizeof(int));

    for (int m = 0; m < nEquation; ++ m)
    {
        cdata->write(q[m],sizeof(RDouble) * nTotal);
    }

    int isUnsteady = parameters->GetIsUnsteady();

    if (!isUnsteady) return;

    RDouble physicalTime = 0;
    GlobalDataBase::GetData("physicalTime",&physicalTime, PHDOUBLE, 1);

    cdata->write(&physicalTime, sizeof(RDouble));

    RDouble **qn1 = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q_unsteady_n1"));
    RDouble **qn2 = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q_unsteady_n2"));

    for (int m = 0; m < nEquation; ++ m)
    {
        cdata->write(qn1[m], sizeof(RDouble) * nTotal);
    }

    for (int m = 0; m < nEquation; ++ m)
    {
        cdata->write(qn2[m],sizeof(RDouble) * nTotal);
    }

    RDouble **res   = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res"           ));
    RDouble **resn1 = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res_unsteady_n1"));
    RDouble **resn2 = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res_unsteady_n2"));

    for (int m = 0; m < nEquation; ++ m)
    {
        cdata->write(res[m], sizeof(RDouble) * nTotal);
    }

    for (int m = 0; m < nEquation; ++ m)
    {
        cdata->write(resn1[m], sizeof(RDouble) * nTotal);
    }

    for (int m = 0; m < nEquation; ++ m)
    {
        cdata->write(resn2[m],sizeof(RDouble) * nTotal);
    }
}

void NSSolverUnstruct::DumpRestartH5(ActionKey *actkey)
{
    using namespace PHMPI;
    int currentProcessorID = GetCurrentProcessorID();

    int version = 1;
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int isUnsteady  = parameters->GetIsUnsteady();
    int outputStepInterval = GlobalDataBase::GetIntParaFromDB("outnstep");
    if (currentProcessorID == GetServerProcessorID())
    {
        WriteData(actkey->filepos, &version, "Version");
        WriteData(actkey->filepos, &outputStepInterval, "outnstep");

        if (isUnsteady == 1)
        {
            RDouble physicalTime   = GlobalDataBase::GetDoubleParaFromDB("physicalTime");
            WriteData(actkey->filepos, &physicalTime, "physicalTime");

            if (IsNeedStatistics())
            {
                int nStatisticalStep = 0;
                GlobalDataBase::GetData("nStatisticalStep", &nStatisticalStep, PHINT, 1);
                WriteData(actkey->filepos, &nStatisticalStep, "nStatisticalStep");
            }
        }
    }

    Grid *grid = GetGrid(actkey->level);
    int gridID = grid->GetZoneID();

    int nEquation = GetNumberOfEquations();

    int nBoundFace = grid->GetNBoundFace();
    int nTotalCell = grid->GetNTotalCell();
    int nTotal = nTotalCell + nBoundFace;

    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    RDouble **primitiveVariables = NewPointer2<RDouble>(nEquation, nTotal);

    for (int iEquation = 0; iEquation < nEquation; ++ iEquation)
    {
        for (int iCell = 0; iCell < nTotal; ++ iCell)
        {
            primitiveVariables[iEquation][iCell] = q[iEquation][iCell];
        }
    }

    hid_t grploc;
    string grpName;

    ostringstream oss;
    oss<<"Group"<<gridID;
    grpName = oss.str();

    grploc = OpenGroup(actkey->filepos, grpName);

    WriteData(grploc, primitiveVariables[0], "q");
    WriteData(grploc, &nTotalCell, "nTotalCell");

    if (isUnsteady == 1)
    {
        RDouble **qn1 = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q_unsteady_n1"));
        for (int iEquation = 0; iEquation < nEquation; ++ iEquation)
        {
            for (int iCell = 0; iCell < nTotal; ++ iCell)
            {
                primitiveVariables[iEquation][iCell] = qn1[iEquation][iCell];
            }
        }
        WriteData(grploc, primitiveVariables[0], "q_unsteady_n1");

        RDouble **qn2 = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q_unsteady_n2"));
        for (int iEquation = 0; iEquation < nEquation; ++ iEquation)
        {
            for (int iCell = 0; iCell < nTotal; ++ iCell)
            {
                primitiveVariables[iEquation][iCell] = qn2[iEquation][iCell];
            }
        }
        WriteData(grploc, primitiveVariables[0], "q_unsteady_n2");

        RDouble **res   = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res"));
        for (int iEquation = 0; iEquation < nEquation; ++ iEquation)
        {
            for (int iCell = 0; iCell < nTotal; ++ iCell)
            {
                primitiveVariables[iEquation][iCell] = res[iEquation][iCell];
            }
        }
        WriteData(grploc, primitiveVariables[0], "res");

        RDouble **resn1 = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res_unsteady_n1"));
        for (int iEquation = 0; iEquation < nEquation; ++ iEquation)
        {
            for (int iCell = 0; iCell < nTotal; ++ iCell)
            {
                primitiveVariables[iEquation][iCell] = resn1[iEquation][iCell];
            }
        }
        WriteData(grploc, primitiveVariables[0], "res_unsteady_n1");

        RDouble **resn2 = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res_unsteady_n2"));
        for (int iEquation = 0; iEquation < nEquation; ++ iEquation)
        {
            for (int iCell = 0; iCell < nTotal; ++ iCell)
            {
                primitiveVariables[iEquation][iCell] = resn2[iEquation][iCell];
            }
        }
        WriteData(grploc, primitiveVariables[0], "res_unsteady_n2");

        RDouble **qAverage = reinterpret_cast< RDouble ** > (grid->GetDataPtr("qAverage"));
        for (int iEquation = 0; iEquation < nEquation; ++ iEquation)
        {
            for (int iCell = 0; iCell < nTotal; ++ iCell)
            {
                primitiveVariables[iEquation][iCell] = qAverage[iEquation][iCell];
            }
        }
        WriteData(grploc, primitiveVariables[0], "qAverage");

        RDouble **tauAverage = reinterpret_cast< RDouble ** > (grid->GetDataPtr("tauAverage"));
        WriteData(grploc, tauAverage[0], "tauAverage");

        RDouble **q2Average = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q2Average"));
        WriteData(grploc, q2Average[0], "q2Average");
    }

    H5Gclose(grploc);

    DelPointer2(primitiveVariables);
}

void NSSolverUnstruct::ReadRestartH5(ActionKey *actkey)
{
    Grid *grid = GetGrid(actkey->level);
    int gridID = grid->GetZoneID();

    Param_NSSolverUnstruct *parameters = GetControlParameters();

    int nEquation = GetNumberOfEquations();

    int outputStepInterval = 0;
    ReadData(actkey->filepos, &outputStepInterval, "outnstep");
    GlobalDataBase::UpdateData("outnstep", &outputStepInterval, PHINT, 1);

    int nBoundFace = grid->GetNBoundFace();
    int nTotalCell = grid->GetNTotalCell();
    int nTotal = nTotalCell + nBoundFace;

    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    RDouble **primitiveVariables = NewPointer2<RDouble>(nEquation, nTotal);

    hid_t grploc;
    string grpName;

    ostringstream oss;
    oss<<"Group"<<gridID;
    grpName = oss.str();

    grploc = OpenGroup(actkey->filepos, grpName);

    int nTotalCellRestart = 0;
    ReadData(grploc, &nTotalCellRestart, "nTotalCell");

    if(nTotalCellRestart != nTotalCell)
    {
        ostringstream erroeInfo;
        erroeInfo << " Error: the cell number in flow.dat is not equal to the cell number in grid file !" << endl;
        TK_Exit::ExceptionExit(erroeInfo.str());
    }

    ReadData(grploc, primitiveVariables[0], "q");
    for (int iEquation = 0; iEquation < nEquation; ++ iEquation)
    {
        for (int iCell = 0; iCell < nTotal; ++ iCell)
        {
            q[iEquation][iCell] = primitiveVariables[iEquation][iCell];
        }
    }

    int isUnsteady = parameters->GetIsUnsteady();
    if (!isUnsteady)
    {
        H5Gclose(grploc);
        DelPointer2(primitiveVariables);

        return;
    }

    RDouble **qn1 = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q_unsteady_n1"));
    RDouble **qn2 = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q_unsteady_n2"));
    RDouble **res   = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res"           ));
    RDouble **resn1 = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res_unsteady_n1"));
    RDouble **resn2 = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res_unsteady_n2"));
    RDouble **qAverage = reinterpret_cast< RDouble ** > (grid->GetDataPtr("qAverage"));
    RDouble **tauAverage = reinterpret_cast< RDouble ** > (grid->GetDataPtr("tauAverage"));
    RDouble **q2Average = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q2Average"));

    int ifStartFromSteadyResults = parameters->GetIfStartFromSteadyResults();
    RDouble physicalTime = 0.0;

    if (ifStartFromSteadyResults)
    {
        PrintToWindow("Restart from steady NS flow field, reset outer step to be zero!\n");
        //! Start from steady flow field.
        //! Reset the outer step when start from steady flow.
        outputStepInterval = 0;
        GlobalDataBase::UpdateData("outnstep", &outputStepInterval, PHINT, 1);

        for (int m = 0; m < nEquation; ++ m)
        {
            for (int iCell = 0; iCell < nTotalCell; ++ iCell)
            {
                qn1[m][iCell] = q[m][iCell];
                qn2[m][iCell] = q[m][iCell];
                res[m][iCell]   = 0.0;
                resn1[m][iCell] = 0.0;
                resn2[m][iCell] = 0.0;
            }
        }
    }
    else
    {
        ReadData(actkey->filepos, &physicalTime, "physicalTime");

        ReadData(grploc, primitiveVariables[0], "q_unsteady_n1");
        for (int iEquation = 0; iEquation < nEquation; ++ iEquation)
        {
            for (int iCell = 0; iCell < nTotal; ++ iCell)
            {
                qn1[iEquation][iCell] = primitiveVariables[iEquation][iCell];
            }
        }

        ReadData(grploc, primitiveVariables[0], "q_unsteady_n2");
        for (int iEquation = 0; iEquation < nEquation; ++ iEquation)
        {
            for (int iCell = 0; iCell < nTotal; ++ iCell)
            {
                qn2[iEquation][iCell] = primitiveVariables[iEquation][iCell];
            }
        }

        ReadData(grploc, primitiveVariables[0], "res");
        for (int iEquation = 0; iEquation < nEquation; ++ iEquation)
        {
            for (int iCell = 0; iCell < nTotal; ++ iCell)
            {
                res[iEquation][iCell] = primitiveVariables[iEquation][iCell];
            }
        }

        ReadData(grploc, primitiveVariables[0], "res_unsteady_n1");
        for (int iEquation = 0; iEquation < nEquation; ++ iEquation)
        {
            for (int iCell = 0; iCell < nTotal; ++ iCell)
            {
                resn1[iEquation][iCell] = primitiveVariables[iEquation][iCell];
            }
        }

        ReadData(grploc, primitiveVariables[0], "res_unsteady_n2");
        for (int iEquation = 0; iEquation < nEquation; ++ iEquation)
        {
            for (int iCell = 0; iCell < nTotal; ++ iCell)
            {
                resn2[iEquation][iCell] = primitiveVariables[iEquation][iCell];
            }
        }
    }

    int nStatisticalStep = 0;

    if(IsNeedStatistics())
    {
        if(ifStartFromSteadyResults)
        {
            for (int iEquation = 0; iEquation < nEquation; ++ iEquation)
            {
                for (int iCell = 0; iCell < nTotalCell; ++ iCell)
                {
                    qAverage[iEquation][iCell] = q[iEquation][iCell];
                }
            }
        }
        else
        {
            ReadData(grploc, qAverage[0], "qAverage");

            ReadData(actkey->filepos, &nStatisticalStep, "nStatisticalStep");
        }
    }

    if(IsNeedReynoldsStressStatistics())
    {
        if(ifStartFromSteadyResults)
        {
            for (int m = 0; m < 6; ++ m)
            {
                for (int iCell = 0; iCell < nTotalCell; ++ iCell)
                {
                    tauAverage[m][iCell] = 0.0;
                    q2Average[m][iCell] = 0.0;
                }
            }
        }
        else
        {
            ReadData(grploc, tauAverage[0], "tauAverage");

            ReadData(grploc, q2Average[0], "q2Average");
        }
    }

    GlobalDataBase::UpdateData("nStatisticalStep", &nStatisticalStep, PHINT, 1);

    GlobalDataBase::UpdateData("physicalTime",&physicalTime, PHDOUBLE, 1);

    H5Gclose(grploc);
    DelPointer2(primitiveVariables);
}

void NSSolverUnstruct::ReadRestart(ActionKey *actkey)
{
    int level = actkey->level;
    UnstructGrid * grid = UnstructGridCast(GetGrid(level));

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    DataContainer * cdata = actkey->GetData();
    cdata->moveto_begin();
    int outnstep = 0;
    cdata->read(&outnstep,sizeof(int));

    GlobalDataBase::UpdateData("outnstep", &outnstep, PHINT, 1);

    int nBoundFace = grid->GetNBoundFace();
    int nTotalCell = grid->GetNTotalCell();
    int nTotal = nTotalCell + nBoundFace;

    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));

    for (int m = 0; m < nEquation; ++ m)
    {
        cdata->read(q[m], sizeof(RDouble) * nTotal);
    }

    int isUnsteady = parameters->GetIsUnsteady();
    if (!isUnsteady) return;

    RDouble **qn1 = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q_unsteady_n1"));
    RDouble **qn2 = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q_unsteady_n2"));

    RDouble **res   = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res"           ));
    RDouble **resn1 = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res_unsteady_n1"));
    RDouble **resn2 = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res_unsteady_n2"));

    int ifStartFromSteadyResults = parameters->GetIfStartFromSteadyResults();

    RDouble physicalTime = 0.0;
    if (ifStartFromSteadyResults)
    {
        PrintToWindow("Restart from steady NS flow field, reset outer step to be zero!\n");
        //! Start from steady flow field.
        //! Reset the outer step when start from steady flow.
        outnstep = 0;
        GlobalDataBase::UpdateData("outnstep", &outnstep, PHINT, 1);

        for (int m = 0; m < nEquation; ++ m)
        {
            for (int iCell = 0; iCell < nTotalCell; ++ iCell)
            {
                qn1[m][iCell] = q[m][iCell];
                qn2[m][iCell] = q[m][iCell];

                res[m][iCell]   = 0.0;
                resn1[m][iCell] = 0.0;
                resn2[m][iCell] = 0.0;
            }
        }
    }
    else
    {
        //! Start from unsteady flow field.
        cdata->read(&physicalTime, sizeof(RDouble));

        for (int m = 0; m < nEquation; ++ m)
        {
            cdata->read(qn1[m], sizeof(RDouble) * nTotal);
        }

        for (int m = 0; m < nEquation; ++ m)
        {
            cdata->read(qn2[m], sizeof(RDouble) * nTotal);
        }

        for (int m = 0; m < nEquation; ++ m)
        {
            cdata->read(res[m],sizeof(RDouble) * nTotal);
        }

        for (int m = 0; m < nEquation; ++ m)
        {
            cdata->read(resn1[m], sizeof(RDouble) * nTotal);
        }

        for (int m = 0; m < nEquation; ++ m)
        {
            cdata->read(resn2[m], sizeof(RDouble) * nTotal);
        }
    }

    GlobalDataBase::UpdateData("physicalTime",&physicalTime, PHDOUBLE, 1);
}

void NSSolverUnstruct::InitFlowAsRestart()
{
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    UnstructGrid *grid = UnstructGridCast(GetGrid(0));
    int outnstep = 0;
    GlobalDataBase::UpdateData("outnstep",&outnstep, PHINT, 1);

    RDouble physicalTime = 0.0;
    GlobalDataBase::UpdateData("physicalTime",&physicalTime, PHDOUBLE, 1);

    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));

    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    RDouble *primitiveVarFarfield = parameters->GetPrimitiveVarFarfield();

    for (int m = 0; m < nEquation; ++ m)
    {
        for (int iCell = 0; iCell < nTotal; ++ iCell)
        {
            q[m][iCell] = primitiveVarFarfield[m];
        }
    }

    int isUnsteady = parameters->GetIsUnsteady();

    if (!isUnsteady)
    {
        RDouble **res = reinterpret_cast <RDouble **> (grid->GetDataPtr("res"));
        for (int m = 0; m < nEquation; ++ m)
        {
            for (int iCell = 0; iCell < nTotal; ++ iCell)
            {
                res[m][iCell] = 0.0;
            }
        }
        return;
    }

    RDouble **qn1 = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q_unsteady_n1"));
    RDouble **qn2 = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q_unsteady_n2"));

    for (int m = 0; m < nEquation; ++ m)
    {
        for (int iCell = 0; iCell < nTotal; ++ iCell)
        {
            qn1[m][iCell] = q[m][iCell];
        }
    }

    for (int m = 0; m < nEquation; ++ m)
    {
        for (int iCell = 0; iCell < nTotal; ++ iCell)
        {
            qn2[m][iCell] = q[m][iCell];
        }
    }

    RDouble **res   = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res"           ));
    RDouble **resn1 = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res_unsteady_n1"));
    RDouble **resn2 = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res_unsteady_n2"));

    //! Store the historical sum of InviscidFlux, ViscousFlux, et al., exclusive of DualTimeSourceFlux, used in UpdateUnsteadyFlow for the Rn in the dualtime method.
    RDouble **resTmp = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res_unsteady_tmp"));

    for (int m = 0; m < nEquation; ++ m)
    {
        for (int iCell = 0; iCell < nTotal; ++ iCell)
        {
            res  [m][iCell] = 0.0;
            resn1[m][iCell] = 0.0;
            resn2[m][iCell] = 0.0;
            resTmp[m][iCell] = 0.0;
        }
    }

    int nStatisticalStep = 0;
    GlobalDataBase::UpdateData("nStatisticalStep", &nStatisticalStep, PHINT, 1);

}

void NSSolverUnstruct::UploadInterfaceData(ActionKey *actkey)
{
    int level = actkey->level;
    UnstructGrid *grid = reinterpret_cast<UnstructGrid *>(GetGrid(level));
    InterfaceInfo *interfaceInfo = grid->GetInterfaceInfo();
    if (!interfaceInfo)
    {
        return;
    }

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    PHSPACE::UploadInterfaceValue(grid, q, "q",  nEquation);

    RDouble **t = reinterpret_cast< RDouble ** > (grid->GetDataPtr("t"));
    PHSPACE::UploadInterfaceValue(grid, t, "t", nTemperatureModel);
 }

void NSSolverUnstruct::DownloadInterfaceData(ActionKey *actkey)
{
    int level = actkey->level;
    UnstructGrid * grid = reinterpret_cast< UnstructGrid * >(GetGrid(level));
    InterfaceInfo * interfaceInfo = grid->GetInterfaceInfo();
    if (!interfaceInfo) return;

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    RDouble **q = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("q"  ) );
    PHSPACE::DownloadInterfaceValue(grid, q, "q", nEquation);

    RDouble **t = reinterpret_cast< RDouble ** > (grid->GetDataPtr("t" ));
    PHSPACE::DownloadInterfaceValue(grid, t, "t", nTemperatureModel);
}

void NSSolverUnstruct::UploadInterpointData(ActionKey *actkey)
{
    int level = actkey->level;
    UnstructGrid * grid = reinterpret_cast< UnstructGrid * >(GetGrid(level));
    InterpointInformation *interpointInformation = grid->GetInterpointInfo();
    if (interpointInformation == 0)
    {
        return;
    }
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    RDouble **qnode = reinterpret_cast<RDouble **>(grid->GetDataPtr("qnode"));
    PHSPACE::UploadInterpointValue(grid, qnode, "qnode",  nEquation);

    RDouble **tnode = reinterpret_cast<RDouble **>(grid->GetDataPtr("tnode"));
    PHSPACE::UploadInterpointValue(grid, tnode, "tnode", nTemperatureModel);
}

void NSSolverUnstruct::DownloadInterpointData(ActionKey *actkey)
{
    int level = actkey->level;
    UnstructGrid *grid = reinterpret_cast<UnstructGrid *>(GetGrid(level));
    InterpointInformation *interpointInformation = grid->GetInterpointInfo();
    if (interpointInformation == 0)
    {
        return;
    }
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    RDouble **qnode = reinterpret_cast<RDouble **>(grid->GetDataPtr("qnode"));
    PHSPACE::DownloadInterpointValue(grid, qnode, "qnode", nEquation);

    RDouble **tnode = reinterpret_cast<RDouble **> (grid->GetDataPtr("tnode"));
    PHSPACE::DownloadInterpointValue(grid, tnode, "tnode", nTemperatureModel);
}

void NSSolverUnstruct::CommunicationInterpointWeight(ActionKey *actkey)
{
    int level = actkey->level;

    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int myid = GetCurrentProcessorID();
        int tag = iZone + level * nZones;
        int send_proc = GetZoneProcessorID(iZone);

        ZoneNeighbor *zoneNeighborForPoint = zoneConnectivityForPoint->GetZoneNeighborForPoint(iZone);
        std::size_t numberOfNeighbor = zoneNeighborForPoint->GetNumberOfNeighbor();

        for (std::size_t iNeighborZone = 0; iNeighborZone < numberOfNeighbor; ++ iNeighborZone)
        {
            int neighborZoneIndex = zoneNeighborForPoint->GetZoneIndexOfNeighbor(iNeighborZone);
            int recv_proc = GetZoneProcessorID(neighborZoneIndex);

            DataContainer * cdata = new DataContainer();

            if (myid == send_proc)
            {
                UnstructGrid *grid = UnstructGridCast(PHSPACE::GetGrid(iZone, level));
                InterpointInformation *interpointInformation = grid->GetInterpointInfo();
                RDouble *nodeWeight = reinterpret_cast<RDouble *>(grid->GetDataPtr("nodeWeight"));

                cdata->moveto_begin();

                int iNeighbor = interpointInformation->FindIthNeighbor(neighborZoneIndex);
                int numberOfInterpointsForNeighbor = interpointInformation->GetNumberOfInterpointsForNeighbor(iNeighbor);
                int *pointIndexForSend = interpointInformation->GetPointIndexForSend(iNeighbor);
                int *interPoint2GlobalPoint = interpointInformation->GetInterPoint2GlobalPoint();

                int iPoint;
                int globalPoint;
                RDouble valueWeight;
                for (int ipointLocal = 0; ipointLocal < numberOfInterpointsForNeighbor; ++ ipointLocal)
                {
                    iPoint = pointIndexForSend[ipointLocal];
                    globalPoint = interPoint2GlobalPoint[iPoint];
                    valueWeight = nodeWeight[globalPoint];

                    cdata->write(&valueWeight, sizeof(RDouble));
                }
            }

            PH_Trade(cdata, send_proc, recv_proc, tag);

            if (myid == recv_proc)
            {
                UnstructGrid *gridNeighbor = UnstructGridCast(PHSPACE::GetGrid(neighborZoneIndex, level));
                InterpointInformation *interpointInformationNeighbor = gridNeighbor->GetInterpointInfo();
                RDouble *nodeWeightTrade = reinterpret_cast<RDouble *>(gridNeighbor->GetDataPtr("nodeWeightTrade"));

                cdata->moveto_begin();

                int iNeighbor = interpointInformationNeighbor->FindIthNeighbor(iZone);
                int numberOfInterpointsForNeighbor = interpointInformationNeighbor->GetNumberOfInterpointsForNeighbor(iNeighbor);
                int *pointIndexForReceive = interpointInformationNeighbor->GetPointIndexForReceive(iNeighbor);
                int *interPoint2GlobalPoint = interpointInformationNeighbor->GetInterPoint2GlobalPoint();

                int iPoint;
                int globalPoint;
                RDouble valueWeight;
                for (int j = 0; j < numberOfInterpointsForNeighbor; ++ j)
                {
                    iPoint = pointIndexForReceive[j];
                    globalPoint = interPoint2GlobalPoint[iPoint];
                    cdata->read(&valueWeight, sizeof(RDouble));

                    nodeWeightTrade[globalPoint] += valueWeight;
                }
            }

            delete cdata;
        }
    }
}

void NSSolverUnstruct::DownloadInterpointWeight(ActionKey *actkey)
{
    int level = actkey->level;

    using namespace PHMPI;

    int myid = GetCurrentProcessorID();
    int nZones = GetNumberofGlobalZones();
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int currentZoneProc = GetZoneProcessorID(iZone);

        if (myid != currentZoneProc)
        {
            continue;
        }

        UnstructGrid *grid = UnstructGridCast(PHSPACE::GetGrid(iZone, level));
        int nTotalNode = grid->GetNTotalNode();
        RDouble *nodeWeight = reinterpret_cast<RDouble *>(grid->GetDataPtr("nodeWeight"));
        RDouble *nodeWeightTrade = reinterpret_cast<RDouble *>(grid->GetDataPtr("nodeWeightTrade"));

        for (int iNode = 0; iNode < nTotalNode; ++ iNode)
        {
            nodeWeight[iNode] += nodeWeightTrade[iNode];
        }

        PHSPACE::SetField(nodeWeightTrade, 0.0, nTotalNode);
    }
}

void NSSolverUnstruct::UploadOversetData(ActionKey *actkey)
{
    int level = actkey->level;
    UnstructGrid * grid = reinterpret_cast< UnstructGrid * >(GetGrid(level));

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q" ));
    PHSPACE::UploadOversetValue(grid, q, "q",  nEquation);

    RDouble **t = reinterpret_cast< RDouble ** > (grid->GetDataPtr("t" ));
    PHSPACE::UploadOversetValue(grid, t, "t",  nTemperatureModel);
}

void NSSolverUnstruct::DownloadOversetData(ActionKey *actkey)
{
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int level = actkey->level;
    UnstructGrid * grid = reinterpret_cast< UnstructGrid * >(GetGrid(level));

    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q" ));
    PHSPACE::DownloadOversetValue(grid, q, "q", nEquation);

    RDouble **t = reinterpret_cast< RDouble ** > (grid->GetDataPtr("t" ));
    PHSPACE::DownloadOversetValue(grid, t, "t",  nTemperatureModel);
}

FieldProxy * NSSolverUnstruct::CreateFieldProxy(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    RDouble **field = NewPointer2<RDouble>(nEquation, nTotal);

    FieldProxy *fieldProxy = new FieldProxy();

    fieldProxy->SetField_UNS(field, true);

    return fieldProxy;
}

FieldProxy * NSSolverUnstruct::GetFieldProxy(Grid *gridIn, const string &field_name)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    RDouble ** field = reinterpret_cast< RDouble ** > (grid->GetDataPtr(field_name));

    FieldProxy * field_proxy = new FieldProxy();

    field_proxy->SetField_UNS(field);

    return field_proxy;
}

FieldProxy * NSSolverUnstruct::GetResidualProxy(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    RDouble **res = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res"));

    FieldProxy *res_proxy = new FieldProxy();

    res_proxy->SetField_UNS(res);

    return res_proxy;
}

void NSSolverUnstruct::RecoverResidual(Grid *gridIn, FieldProxy *rhsProxy)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    RDouble **res = reinterpret_cast<RDouble **> (grid->GetDataPtr("res"));

    if (grid->GetLevel() != 0) 
    {
        RDouble **rhs = rhsProxy->GetField_UNS();
        for (int m = 0; m < nEquation; ++ m)
        {
            for (int iCell = 0; iCell < nTotal; ++ iCell)
            {
                res[m][iCell] = rhs[m][iCell];
            }
        }
    }
}

void NSSolverUnstruct::StoreRhsByResidual(Grid *gridIn, FieldProxy *rhsProxy)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    RDouble **res = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res"));
    RDouble **rhs = rhsProxy->GetField_UNS();

    for (int m = 0; m < nEquation; ++ m)
    {
        for (int iCell = 0; iCell < nTotal; ++ iCell)
        {
            rhs[m][iCell] = res[m][iCell];
        }
    }
}

void NSSolverUnstruct::InitResidual(Grid *gridIn, FieldProxy *rhsProxy)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    RDouble **res = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res"));

    RDouble **rhs = rhsProxy->GetField_UNS();
    for (int m = 0; m < nEquation; ++ m)
    {
        for (int iCell = 0; iCell < nTotal; ++ iCell)
        {
            res[m][iCell] = - rhs[m][iCell];
        }
    }
}

void NSSolverUnstruct::FreeGradientProxy(Grid *grid)
{
}

//! Drive functions to calculate residuals for different orders and dimensions.
void NSSolverUnstruct::Turb_Sengy(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int viscousType = parameters->GetViscousType();

    if (viscousType >= TWO_EQU)
    {
        int nrokplus = 1;
        GlobalDataBase::GetData("nrokplus", &nrokplus, PHINT, 1);

        if (nrokplus > 0)
        {
            RDouble **res  = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res"));
            RDouble *sengy = reinterpret_cast< RDouble * > (grid->GetDataPtr("sengy"  ));
            for (int iCell = 0; iCell < nTotalCell; ++ iCell)
            {
                res[4][iCell] += sengy[iCell];
            }
        }
    }
}

//! Load flow variables stored in grid to q.
void NSSolverUnstruct::LoadQ(Grid *gridIn, FieldProxy *qProxy)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    RDouble **q = qProxy->GetField_UNS();

    RDouble **qold = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar = parameters->GetLaminarNumber();

    for (int m = 0; m < nLaminar; ++ m)
    {
        for (int iCell = 0; iCell < nTotal; ++ iCell)
        {
            q[m][iCell] = qold[m][iCell];
        }
    }
}

//! Restrict defect from fine to coarse grid.
void NSSolverUnstruct::RestrictDefect(Grid *fgrid_in, Grid *cgrid_in)
{
    UnstructGrid *fgrid = UnstructGridCast(fgrid_in);
    UnstructGrid *cgrid = UnstructGridCast(cgrid_in);

    int  fnTotalCell = fgrid->GetNTotalCell();
    int *cell2coarsegridcell    = fgrid->GetCell2CoarseGridCell();

    ZeroResiduals(cgrid);

    RDouble **cres = reinterpret_cast< RDouble ** > (cgrid->GetDataPtr("res"));
    RDouble **fres = reinterpret_cast< RDouble ** > (fgrid->GetDataPtr("res"));

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar = parameters->GetLaminarNumber();

    for (int m = 0; m < nLaminar; ++ m)
    {
        for (int iCell = 0; iCell < fnTotalCell; ++ iCell)
        {
            cres[m][ cell2coarsegridcell[iCell] ] -= fres[m][iCell];
        }
    }
}

//! Put correction back on the coarse grid.
void NSSolverUnstruct::PutCorrectionBack(Grid *gridIn, FieldProxy *qProxy)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();

    RDouble **q = qProxy->GetField_UNS();

    RDouble **qold = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar = parameters->GetLaminarNumber();

    for (int m = 0; m < nLaminar; ++ m)
    {
        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            qold[m][iCell] += q[m][iCell];
        }
    }
}

//! Correct variable in fine grid using correction in coarse grid.
void NSSolverUnstruct::CorrectFineGrid(Grid *fgrid_in, Grid *cgrid_in)
{
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int mgProlongationType = parameters->GetMgProlongationType();

    if(mgProlongationType == 2)
    {
        CorrectFineGridHighOrder(fgrid_in, cgrid_in);
    }
    else
    {
        CorrectFineGridZeroOrder(fgrid_in, cgrid_in);
    }
}

void NSSolverUnstruct::CorrectFineGridZeroOrder(Grid *fgrid_in, Grid *cgrid_in)
{
    UnstructGrid *fgrid = UnstructGridCast(fgrid_in);
    UnstructGrid *cgrid = UnstructGridCast(cgrid_in);

    int fnTCell = fgrid->GetNTotalCell();

    int * cell2coarsegridcell   = fgrid->GetCell2CoarseGridCell();

    RDouble **fq = reinterpret_cast< RDouble ** > (fgrid->GetDataPtr("q"));
    RDouble **cq = reinterpret_cast< RDouble ** > (cgrid->GetDataPtr("q"));

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();

    RDouble mgCorrectionLimit = parameters->GetMgCorrectionLimit();

    for (int m = 0; m < nLaminar; ++ m)
    {
        for (int iCell = 0; iCell < fnTCell; ++ iCell)
        {
            int coarseCell = cell2coarsegridcell[iCell];

            RDouble fgridValue = fq[m][iCell];
            RDouble correction = cq[m][coarseCell];
            RDouble DQ         = MIN(fabs(correction), fabs(fgridValue) * mgCorrectionLimit);
            fq[m][iCell]      += SIGN(DQ, correction);
        }
    }
}

void NSSolverUnstruct::CorrectFineGridHighOrder(Grid *fgrid_in, Grid *cgrid_in)
{
    UnstructGrid *fgrid = UnstructGridCast(fgrid_in);
    UnstructGrid *cgrid = UnstructGridCast(cgrid_in);

    int fnTNode = fgrid->GetNTotalNode();
    int fnTCell = fgrid->GetNTotalCell();
    int fnBFace = fgrid->GetNBoundFace();
    int fnTFace = fgrid->GetNTotalFace();

    int * cell2CoarseGridCell = fgrid->GetCell2CoarseGridCell();
    int * fleft_cell_of_face  = fgrid->GetLeftCellOfFace();
    int * fright_cell_of_face = fgrid->GetRightCellOfFace();

    int ** face2nodeArray = fgrid->GetFace2NodeArray();
    int * fnode_number_of_each_face = fgrid->GetNodeNumberOfEachFace();

    RDouble **fq = reinterpret_cast< RDouble ** > (fgrid->GetDataPtr("q"));
    RDouble **cq = reinterpret_cast< RDouble ** > (cgrid->GetDataPtr("q"));

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();

    RDouble minRatio[5], maxRatio[5];
    for (int m = 0; m < 5; ++ m)
    {
        minRatio[m] = LARGE;
        maxRatio[m] = -LARGE;
    }

    RDouble mgCorrectionLimit = parameters->GetMgCorrectionLimit();

    for (int m = 0; m < nLaminar; ++ m)
    {
        //! On coarse mesh, average the cell data onto nodes.
        RDouble * fqNode  = new RDouble [ fnTNode ];
        int * fnodeCount = new int [fnTNode];
        for (int iNode = 0; iNode < fnTNode; ++ iNode)
        {
            fqNode[iNode]     = 0.0;
            fnodeCount[iNode] = 0;
        }

        for (int iFace = 0; iFace < fnBFace; ++ iFace)
        {
            int le = fleft_cell_of_face[iFace];
            int leftCoarseCell = cell2CoarseGridCell[le];

            int nNode = fnode_number_of_each_face[iFace];
            for (int iNode = 0; iNode < nNode; ++ iNode)
            {
                int fnodeID = face2nodeArray[iFace][iNode];

                fqNode[fnodeID] += cq[m][leftCoarseCell];
                ++ fnodeCount[fnodeID];
            }
        }

        for (int iFace = fnBFace; iFace < fnTFace; ++ iFace)
        {
            int le = fleft_cell_of_face[iFace];
            int re = fright_cell_of_face[iFace];
            int leftCoarseCell  = cell2CoarseGridCell[le];
            int rightCoarseCell = cell2CoarseGridCell[re];

            int nNode = fnode_number_of_each_face[iFace];
            for (int iNode = 0; iNode < nNode; ++ iNode)
            {
                int fnodeID = face2nodeArray[iFace][iNode];

                fqNode[fnodeID] += cq[m][leftCoarseCell];
                fqNode[fnodeID] += cq[m][rightCoarseCell];

                fnodeCount[fnodeID] += 2;
            }
        }

        for (int iNode = 0; iNode < fnTNode; ++ iNode)
        {
            fqNode[iNode] /= fnodeCount[iNode];
        }
        delete [] fnodeCount;

        //! On fine mesh, average the nodes data to cells.
        int * fcellCount   = new int [fnTCell];
        RDouble * fcellData = new RDouble [fnTCell];
        for (int iCell = 0; iCell < fnTCell; ++ iCell)
        {
            fcellCount[iCell] = 0;
            fcellData[iCell]  = 0.0;
        }
        
        for (int iFace = 0; iFace < fnBFace; ++ iFace)
        {
            int le = fleft_cell_of_face[iFace];

            int nNode = fnode_number_of_each_face[iFace];
            for (int iNode = 0; iNode < nNode; ++ iNode)
            {
                int fnodeID = face2nodeArray[iFace][iNode];

                fcellData[le] += fqNode[fnodeID];
                ++ fcellCount[le];
            }
        }

        for (int iFace = fnBFace; iFace < fnTFace; ++ iFace)
        {
            int le = fleft_cell_of_face[iFace];
            int re = fright_cell_of_face[iFace];

            int nNode = fnode_number_of_each_face[iFace];
            for (int iNode = 0; iNode < nNode; ++ iNode)
            {
                int fnodeID = face2nodeArray[iFace][iNode];

                fcellData[le] += fqNode[fnodeID];
                ++ fcellCount[le];

                fcellData[re] += fqNode[fnodeID];
                ++ fcellCount[re];
            }
        }
        delete [] fqNode;

        //! Add the correction onto fine grid.
        for (int iCell = 0; iCell < fnTCell; ++ iCell)
        {
            RDouble fgridValue = fq[m][iCell];
            RDouble correction = fcellData[iCell] / fcellCount[iCell];
            RDouble DQ         = MIN(fabs(correction), fabs(fgridValue) * mgCorrectionLimit);
            fq[m][iCell]      += SIGN(DQ, correction);

            if(m == IDX::IP || m == IDX::IR)
            {
                if (fq[m][iCell] <= 1.0e-10)
                {
                    cout << " Negative found after correction, variable = " << m 
                         << ", find grid Value = " << fgridValue << ", correction = " << correction 
                         << ", fcellID = " << iCell << ", ccellID = " << cell2CoarseGridCell[iCell] 
                         << ", level = " << fgrid->GetLevel() << endl;
                }
            }
        }
        
        delete [] fcellCount;
        delete [] fcellData;
    }
}

void NSSolverUnstruct::InterpolatFineGrid(Grid *fineGrid, Grid *coarseGrid)
{
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    //! Interpolate vist only if for turbulent simulation is run on coarse grid.
    //!   - 0: do not run turbulent simulation on coarse grid.
    //!   - 1: run turbulent simulation on coarse grid.
    int isRunTurbOnCoarseGrid = 0;

    RDouble **fineQ   = reinterpret_cast<RDouble **> (fineGrid->GetDataPtr("q"));
    RDouble **coarseQ = reinterpret_cast<RDouble **> (coarseGrid->GetDataPtr("q"));

    for (int m = 0; m < nEquation; ++ m)
    {
        InterpolatQ(fineGrid, fineQ[m], coarseGrid, coarseQ[m]);
    }

    RDouble *fineGama   = reinterpret_cast<RDouble *> (fineGrid->GetDataPtr("gama"));
    RDouble *coarseGama = reinterpret_cast<RDouble *> (coarseGrid->GetDataPtr("gama"));

    InterpolatQ(fineGrid, fineGama, coarseGrid, coarseGama);

    int viscousType = parameters->GetViscousType();
    if (viscousType > INVISCID)
    {
        RDouble *fineViscousLaminar   = reinterpret_cast<RDouble *> (fineGrid->GetDataPtr("visl"));
        RDouble *coarseViscousLaminar = reinterpret_cast<RDouble *> (coarseGrid->GetDataPtr("visl"));

        InterpolatQ(fineGrid, fineViscousLaminar, coarseGrid, coarseViscousLaminar);

        if (viscousType > LAMINAR && isRunTurbOnCoarseGrid)
        {
            RDouble *fineViscousTurbulence   = reinterpret_cast<RDouble *> (fineGrid->GetDataPtr("vist"));
            RDouble *coarseViscousTurbulence = reinterpret_cast<RDouble *> (coarseGrid->GetDataPtr("vist"));

            InterpolatQ(fineGrid, fineViscousTurbulence, coarseGrid, coarseViscousTurbulence);
        }
    }
}

//! Load rhs to residual res stored in grid.
void NSSolverUnstruct::LoadResiduals(Grid *gridIn, FieldProxy *rhsProxy)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    RDouble **rhs = rhsProxy->GetField_UNS();

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    RDouble **res = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res"));

    for (int m = 0; m < nEquation; ++ m)
    {
        for (int iCell = 0; iCell < nTotal; ++ iCell)
        {
            res[m][iCell] = - rhs[m][iCell];
        }
    }
}

RDouble NSSolverUnstruct::UnsteadyConvergence(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int isUnsteady = parameters->GetIsUnsteady();

    if (!isUnsteady)
    {
        return zero;
    }

    RDouble **q   = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    RDouble **qn1 = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q_unsteady_n1"));
    RDouble **qn2 = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q_unsteady_n2"));
    RDouble **res = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res"));

    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();

    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    RDouble refGama = parameters->GetRefGama();

    RDouble *prim0 = new RDouble[nEquation];
    RDouble *prim1 = new RDouble[nEquation];
    RDouble *prim2 = new RDouble[nEquation];

    RDouble *qcsv0 = new RDouble[nEquation];
    RDouble *qcsv1 = new RDouble[nEquation];
    RDouble *qcsv2 = new RDouble[nEquation];

    RDouble Ttr[3] = {0.0}, Tv[3] = {0.0}, Te[3] = {0.0};

    using namespace GAS_SPACE;

    RDouble sum1 = zero;
    RDouble sum2 = zero;

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        for (int m = 0; m < nEquation; ++ m)
        {
            prim0[m] = q  [m][iCell];
            prim1[m] = qn1[m][iCell];
            prim2[m] = qn2[m][iCell];
        }

        if (nTemperatureModel > 1)
        {
            gas->GetTemperature(prim0, Ttr[0], Tv[0], Te[0]);
            gas->GetTemperature(prim1, Ttr[1], Tv[1], Te[1]);
            gas->GetTemperature(prim2, Ttr[2], Tv[2], Te[2]);
        }

        //! Primitive2Conservative function does not need to input gama in fact.
        //! It will be modified later.
        gas->Primitive2Conservative(prim0, refGama, Tv[0], Te[0], qcsv0);
        gas->Primitive2Conservative(prim1, refGama, Tv[1], Te[1], qcsv1);
        gas->Primitive2Conservative(prim2, refGama, Tv[2], Te[2], qcsv2);

        //! We need nLaminar here.
        for (int m = 0; m < nLaminar; ++ m)
        {
            //! Now res has been converted to dq. It means that residual is not right hand side term,but dq.
            RDouble dq_p =  res[m][iCell];       // qn+1,p+1 - qn+1,p
            RDouble dq_n =  qcsv0[m] - qcsv1[m]; // qn+1,p+1 - qn
            sum1 += dq_p * dq_p;
            sum2 += dq_n * dq_n;
        }
    }

    delete []prim0;
    delete []prim1;
    delete []prim2;
    delete []qcsv0;
    delete []qcsv1;
    delete []qcsv2;

    RDouble cvg = sqrt(ABS(sum1 / (sum2 + SMALL)));
    return cvg;
}

void NSSolverUnstruct::UpdateUnsteadyFlow(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int isUnsteady = parameters->GetIsUnsteady();

    if (!isUnsteady) return;

    RDouble **q   = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    RDouble **qn1 = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q_unsteady_n1"));
    RDouble **qn2 = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q_unsteady_n2"));

    RDouble **resn1 = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res_unsteady_n1"));
    RDouble **resn2 = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res_unsteady_n2"));
    RDouble **resTmp = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res_unsteady_tmp"));

    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    for (int m = 0; m < nEquation; ++ m)
    {
        for (int iCell = 0; iCell < nTotal; ++ iCell)
        {
            qn2[ m ][ iCell ] = qn1[ m ][ iCell ];
            qn1[ m ][ iCell ] = q  [ m ][ iCell ];

            resn2[ m ][ iCell ] = resn1 [ m ][ iCell] ;
            resn1[ m ][ iCell ] = resTmp[ m ][ iCell] ; //! Here the current outnstep is over, the value of the stored resTmp should be assigned to resn1 for the next outnstep. It should be noticed that resTmp only contain the inviscid and viscous flux.
        }
    }

    //! Statistical variables for unsteady simulation.
    if (IsNeedStatistics())
    {
        int nStatisticalStep   = GlobalDataBase::GetIntParaFromDB("nStatisticalStep");
        RDouble statisticalTimePeriod = GlobalDataBase::GetDoubleParaFromDB("statisticalTimePeriod");
        RDouble physicalTimeStep = GlobalDataBase::GetDoubleParaFromDB("physicalTimeStep");

        RDouble c1 = 1.0 / (nStatisticalStep * 1.0);
        if(statisticalTimePeriod > 0.0)
        {
            c1 = MAX(c1, physicalTimeStep / statisticalTimePeriod);
        }
        RDouble c2 = 1.0 - c1;

        RDouble **qAverage = reinterpret_cast< RDouble ** > (grid->GetDataPtr("qAverage"));
        for (int m = 0; m < nEquation; ++ m)
        {
            for (int iCell = 0; iCell < nTotal; ++ iCell)
            {
                qAverage[m][iCell] = c2 * qAverage[m][iCell] + c1 * q[m][iCell];
            }
        }
    }

    //! Added by zzp 202108, Statistical Reynolds stress for unsteady simulation.
    if (IsNeedReynoldsStressStatistics())
    {
        int nStatisticalStep   = GlobalDataBase::GetIntParaFromDB("nStatisticalStep");
        RDouble statisticalTimePeriod = GlobalDataBase::GetDoubleParaFromDB("statisticalTimePeriod");
        RDouble physicalTimeStep = GlobalDataBase::GetDoubleParaFromDB("physicalTimeStep");

        RDouble c1 = 1.0 / (nStatisticalStep * 1.0);
        if(statisticalTimePeriod > 0.0)
        {
            c1 = MAX(c1, physicalTimeStep / statisticalTimePeriod);
        }
        RDouble c2 = 1.0 - c1;

        RDouble **qAverage = reinterpret_cast< RDouble ** > (grid->GetDataPtr("qAverage"));
        RDouble **tauAverage = reinterpret_cast< RDouble ** > (grid->GetDataPtr("tauAverage"));
        RDouble **q2Average = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q2Average"));

        int statisticMethod = parameters->GetStatisticMethod();
        if(statisticMethod == 0)
        {
            for (int iCell = 0; iCell < nTotal; ++ iCell)
            {
                RDouble u = q[IU][iCell];
                RDouble v = q[IV][iCell];
                RDouble w = q[IW][iCell];

                q2Average[0][iCell] = c2 * q2Average[0][iCell] + c1 * u * u;
                q2Average[1][iCell] = c2 * q2Average[1][iCell] + c1 * v * v;
                q2Average[2][iCell] = c2 * q2Average[2][iCell] + c1 * w * w;
                q2Average[3][iCell] = c2 * q2Average[3][iCell] + c1 * u * v;
                q2Average[4][iCell] = c2 * q2Average[4][iCell] + c1 * u * w;
                q2Average[5][iCell] = c2 * q2Average[5][iCell] + c1 * v * w;
            }

            for (int iCell = 0; iCell < nTotal; ++ iCell)
            {
                RDouble uavg = qAverage[IU][iCell];
                RDouble vavg = qAverage[IV][iCell];
                RDouble wavg = qAverage[IW][iCell];

                tauAverage[0][iCell] = q2Average[0][iCell] - uavg * uavg;
                tauAverage[1][iCell] = q2Average[1][iCell] - vavg * vavg;
                tauAverage[2][iCell] = q2Average[2][iCell] - wavg * wavg;
                tauAverage[3][iCell] = q2Average[3][iCell] - uavg * vavg;
                tauAverage[4][iCell] = q2Average[4][iCell] - uavg * wavg;
                tauAverage[5][iCell] = q2Average[5][iCell] - vavg * wavg;
            }
        }

        else
        {
            for (int iCell = 0; iCell < nTotal; ++ iCell)
            {
                RDouble uprm = q[IU][iCell] - qAverage[IU][iCell];
                RDouble vprm = q[IV][iCell] - qAverage[IV][iCell];
                RDouble wprm = q[IW][iCell] - qAverage[IW][iCell];

                tauAverage[0][iCell] = c2 * tauAverage[0][iCell] + c1 * uprm * uprm;
                tauAverage[1][iCell] = c2 * tauAverage[1][iCell] + c1 * vprm * vprm;
                tauAverage[2][iCell] = c2 * tauAverage[2][iCell] + c1 * wprm * wprm;
                tauAverage[3][iCell] = c2 * tauAverage[3][iCell] + c1 * uprm * vprm;
                tauAverage[4][iCell] = c2 * tauAverage[4][iCell] + c1 * uprm * wprm;
                tauAverage[5][iCell] = c2 * tauAverage[5][iCell] + c1 * vprm * wprm;
            }
        }
    }
}

void NSSolverUnstruct::StoreHistoricalResidual(Grid *gridIn)
{
    ;
}

void NSSolverUnstruct::DualTimeSource(Grid *gridIn)
{
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int isUnsteady = parameters->GetIsUnsteady();
    if (!isUnsteady)
    {
        return;
    }

    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();

    RDouble **q     = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble **qn1   = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_unsteady_n1"));
    RDouble **qn2   = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_unsteady_n2"));

    RDouble **res   = reinterpret_cast<RDouble **> (grid->GetDataPtr("res"));
    RDouble **resn1 = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_unsteady_n1"));
    RDouble **resn2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_unsteady_n2"));
    RDouble **resTmp = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_unsteady_tmp"));

    RDouble *vol   = grid->GetCellVolume(); //! The volume at current timestep. If it is dualtime step method, it represents the volume at timestep of n+1.
    RDouble *voln1 = grid->GetCellVolume(); //! The volume at timestep of n  .
    RDouble *voln2 = grid->GetCellVolume(); //! the volume at timestep of n-1.

    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    RDouble refGama = parameters->GetRefGama();

    //! the primitive variables
    RDouble *prim0 = new RDouble[nEquation];
    RDouble *prim1 = new RDouble[nEquation];
    RDouble *prim2 = new RDouble[nEquation];

    //! the conservative variables
    RDouble *qcsv0 = new RDouble[nEquation];
    RDouble *qcsv1 = new RDouble[nEquation];
    RDouble *qcsv2 = new RDouble[nEquation];

    RDouble Ttr[3] = {0.0}, Tv[3] = {0.0}, Te[3] = {0.0};

    using namespace GAS_SPACE;

    //! Store the historical sum of InviscidFlux, ViscousFlux, et al., exclusive of DualTimeSourceFlux, used in UpdateUnsteadyFlow for the Rn in the dualtime method.
        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            for (int m = 0; m < nEquation; ++ m)
            {
                resTmp[m][iCell] = res [m][iCell];
            }
        }

    //! Computation of dualtime coefficients, including three coefficients for Residual of R(p), R(n) and R(n-1),
    //! and three coefficients for conservative variables of qcsv(p), qcsv(n), and qcsv(n-1),
    //! the fourth coefficient is also used for the real    time part of the diagonal data for the matrix of the left term when LUSGS is running,
    //! the last   coefficient is      used for the virtual time part of the diagonal data for the matrix of the left term when LUSGS is running.
    RDouble dualTimeCoefficient[7];
    const int methodOfDualTime = GlobalDataBase::GetIntParaFromDB("methodOfDualTime");
    ComputeDualTimeCoefficient(methodOfDualTime, dualTimeCoefficient);

    RDouble dualTimeResC1 = dualTimeCoefficient[0];
    RDouble dualTimeResC2 = dualTimeCoefficient[1];
    RDouble dualTimeResC3 = dualTimeCoefficient[2];
    RDouble dualTimeQC1   = dualTimeCoefficient[3];
    RDouble dualTimeQC2   = dualTimeCoefficient[4];
    RDouble dualTimeQC3   = dualTimeCoefficient[5];

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        for (int m = 0; m < nEquation; ++ m)
        {
            prim0[m] = q  [m][iCell];
            prim1[m] = qn1[m][iCell];
            prim2[m] = qn2[m][iCell];
        }
        if (nTemperatureModel > 1)
        {
            gas->GetTemperature(prim0, Ttr[0], Tv[0], Te[0]);
            gas->GetTemperature(prim1, Ttr[1], Tv[1], Te[1]);
            gas->GetTemperature(prim2, Ttr[2], Tv[2], Te[2]);
        }
        //! In fact, there is no need to input the variable of refGama for the function of Primitive2Conservative, it is to be improved.
        gas->Primitive2Conservative(prim0, refGama, Tv[0], Te[0], qcsv0);
        gas->Primitive2Conservative(prim1, refGama, Tv[1], Te[1], qcsv1);
        gas->Primitive2Conservative(prim2, refGama, Tv[2], Te[2], qcsv2);

        for (int m = 0; m < nEquation; ++ m)
        {
            RDouble dualSrcRes = dualTimeResC1 * res  [ m ][ iCell ] +
                                dualTimeResC2 * resn1[ m ][ iCell ] +
                                dualTimeResC3 * resn2[ m ][ iCell ];

            RDouble dualSrcQ   = dualTimeQC1 * qcsv0[ m ] * vol  [ iCell ] + 
                                dualTimeQC2 * qcsv1[ m ] * voln1[ iCell ] + 
                                dualTimeQC3 * qcsv2[ m ] * voln2[ iCell ];

            res[m][iCell] = dualSrcRes + dualSrcQ;
        }
    }

    delete []prim0;
    delete []prim1;
    delete []prim2;
    delete []qcsv0;
    delete []qcsv1;
    delete []qcsv2;
}

void NSSolverUnstruct::ChemicalSource(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    RDouble *vol  = grid->GetCellVolume();

    RDouble **q   = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    RDouble **t   = reinterpret_cast< RDouble ** > (grid->GetDataPtr("t"));
    RDouble **res = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res"));

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nm = parameters->GetNSEquationNumber();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    RDouble **src = NewPointer2<RDouble>(nLaminar, nTotalCell);

    using namespace GAS_SPACE;
    using namespace IDX;

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        gas->ChemicalSource(q, t, vol, src, iCell);
        for (int m = nm; m < nLaminar; ++ m)
        {
            res[m][iCell] += src[m - nm][iCell];
        }
    }

    DelPointer2(src);
}


void NSSolverUnstruct::ComputeFirstLayerGridHeight(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();

    RDouble *xfc = grid->GetFaceCenterX();
    RDouble *yfc = grid->GetFaceCenterY();
    RDouble *zfc = grid->GetFaceCenterZ();

    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    RDouble *xfn  = grid->GetFaceNormalX();
    RDouble *yfn  = grid->GetFaceNormalY();
    RDouble *zfn  = grid->GetFaceNormalZ();

    RDouble *vol  = grid->GetCellVolume();
    RDouble *area = grid->GetFaceArea();



    RDouble *firstLayerHeight = reinterpret_cast <RDouble *> (grid->GetDataPtr("firstLayerHeight"));

    using namespace IDX;

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegion = unstructBCSet->GetnBCRegion();

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();

        vector<int> *faceIndex = bcRegion->GetFaceIndex();
        for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
        {
            // iFace is the face number in the set of faceIndex.
            int iFace = *iter;
            if (!IsWall(bcType) && bcType != PHENGLEI::ABLATION_SURFACE)
            {
                continue;
            }

#ifdef USE_ALTERNATIVE_CODE
            RDouble cellVolume = gridCellVolume(iWall, jWall, kWall); //The volume of first cell.
            RDouble deltaHeight = half * cellVolume / surfaceArea;
#else
            //! The projection of the dy called the distance between the center of the cell and that of the surface on normal vector.
            RDouble deltaHeight = fabs((xcc - xfc) *xfn[iFace] + (ycc - yfc) * yfn[iFace] + (zcc - zfc) * zfn[iFace]);
#endif
        }
    }
}

void NSSolverUnstruct::ComputeGamaAndTemperature(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    Param_NSSolverUnstruct *parameters = GetControlParameters();

    RDouble **t   = reinterpret_cast< RDouble ** > (grid->GetDataPtr("t"));
    RDouble **q   = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    RDouble *gama = reinterpret_cast< RDouble *  > (grid->GetDataPtr("gama"));

    RDouble refGama = parameters->GetRefGama();

    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    int mTT = parameters->GetmTT();
    int mTV = parameters->GetmTV();
    int mTE = parameters->GetmTE();

    using namespace GAS_SPACE;
    using namespace IDX;

    if (nChemical == 1)
    {
        RDouble *prim = new RDouble[nEquation];
        for (int iCell = 0; iCell < nTotal; ++ iCell)
        {
            for (int m = 0; m < nEquation; ++ m)
            {
                prim[m] = q[m][iCell];
            }
            gas->GetSpecificHeatRatioAndTemperatute(prim, gama[iCell], t[mTT][iCell], t[mTV][iCell], t[mTE][iCell]);
            if (nTemperatureModel > 1)
            {
                t[ITV][iCell] = t[mTE][iCell];
                if (nTemperatureModel == 3)
                {
                    t[ITE][iCell] = t[mTE][iCell];
                }
            }
        }
        delete []prim;
    }
    else
    {
        RDouble coefficientofStateEquation = gas->GetCoefficientOfStateEquation();

        for (int iCell = 0; iCell < nTotal; ++ iCell)
        {
            RDouble &rm = q[IR][iCell];
            RDouble &pm = q[IP][iCell];
            gama[iCell] = refGama;
            t[ITT][iCell] = pm / (coefficientofStateEquation * rm);
        }
    }
}

void NSSolverUnstruct::ComputeNodeValue(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotalNode = grid->GetNTotalNode();

    RDouble *x = grid->GetX();
    RDouble *y = grid->GetY();
    RDouble *z = grid->GetZ();

    RDouble *xfc = grid->GetFaceCenterX();
    RDouble *yfc = grid->GetFaceCenterY();
    RDouble *zfc = grid->GetFaceCenterZ();

    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    RDouble **q = reinterpret_cast<RDouble **>(grid->GetDataPtr("q"));
    RDouble **t = reinterpret_cast<RDouble **>(grid->GetDataPtr("t"));

    RDouble     **qNode = reinterpret_cast<RDouble **>(grid->GetDataPtr("qnode"));
    RDouble     **tNode = reinterpret_cast<RDouble **>(grid->GetDataPtr("tnode"));
    RDouble *nodeWeight = reinterpret_cast<RDouble *>(grid->GetDataPtr("nodeWeight"));

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    int *nodeNumberOfEachFace = grid->GetNodeNumberOfEachFace();
    int **face2nodeArray      = grid->GetFace2NodeArray();

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegion = unstructBCSet->GetnBCRegion();

    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        nodeWeight[iNode] = 0.0;
        for(int m =0; m < nEquation; ++ m)
        {
            qNode[m][iNode] = 0; 
        }
        tNode[0][iNode] = 0;
    }

    //! Node BC type.
    int * nodeBC = NewPointer< int > (nTotalNode);
    PHSPACE::SetField(nodeBC, PHENGLEI::NO_BOUNDARY_CONDITION, nTotalNode);

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC* bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();
        if (bcType == PHENGLEI::SOLID_SURFACE)
        {
            vector<int>* faceIndex = bcRegion->GetFaceIndex();
            for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
            {
                //! iFace is the face number in the set of faceIndex.
                int iFace = *iter;
                for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
                {            
                    int point = face2nodeArray[iFace][jNode];
                    nodeBC[point] = bcType;
                } 
            }
        }
    }
   
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC* bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();
        if (bcType == PHENGLEI::FARFIELD)
        {
            vector<int>* faceIndex = bcRegion->GetFaceIndex();
            for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
            {
                //! iFace is the face number in the set of faceIndex.
                int iFace = *iter;

                //! Far field has the second highest priority.
                for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
                {            
                    int point = face2nodeArray[iFace][jNode];
                    if(nodeBC[point] == PHENGLEI::SOLID_SURFACE)
                    {
                        continue;
                    }

                    nodeBC[point] = bcType;
                } 
            }
        }
    }

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();

        vector<int> *faceIndex = bcRegion->GetFaceIndex();
        for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
        {
            //! iFace is the face number in the set of faceIndex.
            int iFace = *iter;
            if(bcType == PHENGLEI::SOLID_SURFACE || bcType == PHENGLEI::FARFIELD ||
                bcType == PHENGLEI::SYMMETRY || bcType == PHENGLEI::INTERFACE)
            {
                continue;
            }

            //! other BC, except symmetry and interface.
            for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
            {            
                int point = face2nodeArray[iFace][jNode];
                if(nodeBC[point] == PHENGLEI::SOLID_SURFACE || nodeBC[point] == PHENGLEI::FARFIELD)
                {
                    continue;
                }

                nodeBC[point] = bcType;
            } 
        }
    }

    //! Boundary faces.
    RDouble ** qOnBCFace = NewPointer2< RDouble >(nEquation + 1, nBoundFace);  // '1' is temperature.
    PHSPACE::SetField(qOnBCFace, nEquation + 1, nBoundFace, 0.0);

    //! Init q on BC faces.
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();

        vector<int> *faceIndex = bcRegion->GetFaceIndex();
        for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
        {
            // iFace is the face number in the set of faceIndex.
            int iFace = *iter;
            if (bcType == PHENGLEI::INTERFACE || bcType == PHENGLEI::SYMMETRY)
            {
                continue;
            }

            int le = leftCellOfFace [iFace];
            int re = rightCellOfFace[iFace];

            for (int m = 0; m < nEquation; ++ m)
            {
                qOnBCFace[m][iFace] = half * (q[m][le] + q[m][re]);
            }

            qOnBCFace[nEquation][iFace] = half * (t[0][le] + t[0][re]);
        }
    }

    //! Step1: solid surface.
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC* bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();
        if (bcType == PHENGLEI::SOLID_SURFACE)
        {
            vector<int>* faceIndex = bcRegion->GetFaceIndex();
            for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
            {
                //! iFace is the face number in the set of faceIndex.
                int iFace = *iter;
                for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
                {
                    int point = face2nodeArray[iFace][jNode];

                    RDouble dx = x[point] - xfc[iFace];
                    RDouble dy = y[point] - yfc[iFace];
                    RDouble dz = z[point] - zfc[iFace];
                    RDouble dist = DISTANCE(dx, dy, dz);
                    RDouble weightTemp = 1.0 / dist;                

                    for (int m = 0; m < nEquation; ++ m)
                    {
                        RDouble faceQtemp = qOnBCFace[m][iFace];
                        qNode[m][point] += faceQtemp * weightTemp;
                    }            
                    tNode[0][point] += qOnBCFace[nEquation][iFace] * weightTemp;

                    nodeWeight[point] += weightTemp;
                }
            }
        }
    }
    
    //! Step2: Far-field.
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC* bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();
        if (bcType == PHENGLEI::FARFIELD)
        {
            vector<int>* faceIndex = bcRegion->GetFaceIndex();
            for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
            {
                //! iFace is the face number in the set of faceIndex.
                int iFace = *iter;
                for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
                {
                    int point = face2nodeArray[iFace][jNode];
                    if(nodeBC[point] == PHENGLEI::SOLID_SURFACE)
                    {
                        continue;
                    }

                    RDouble dx = x[point] - xfc[iFace];
                    RDouble dy = y[point] - yfc[iFace];
                    RDouble dz = z[point] - zfc[iFace];
                    RDouble dist = DISTANCE(dx, dy, dz);
                    RDouble weightTemp = 1.0 / dist;                

                    for (int m = 0; m < nEquation; ++ m)
                    {
                        RDouble faceQtemp = qOnBCFace[m][iFace];
                        qNode[m][point] += faceQtemp * weightTemp;
                    }            
                    tNode[0][point] += qOnBCFace[nEquation][iFace] * weightTemp;

                    nodeWeight[point] += weightTemp;                
                }
            }
        }
    }

    //! Step3: other BC except solid surface/far field/symmetry/Interface.
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC* bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();
        if (bcType != PHENGLEI::SOLID_SURFACE && bcType != PHENGLEI::FARFIELD && 
            bcType != PHENGLEI::SYMMETRY      && bcType != PHENGLEI::INTERFACE)
        {
            vector<int>* faceIndex = bcRegion->GetFaceIndex();
            for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
            {
                //! iFace is the face number in the set of faceIndex.
                int iFace = *iter;
                for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
                {
                    int point = face2nodeArray[iFace][jNode];
                    if(nodeBC[point] == PHENGLEI::SOLID_SURFACE || nodeBC[point] == PHENGLEI::FARFIELD)
                    {
                        continue;
                    }

                    RDouble dx = x[point] - xfc[iFace];
                    RDouble dy = y[point] - yfc[iFace];
                    RDouble dz = z[point] - zfc[iFace];
                    RDouble dist = DISTANCE(dx, dy, dz);
                    RDouble weightTemp = 1.0 / dist;                

                    for (int m = 0; m < nEquation; ++ m)
                    {
                        RDouble faceQtemp = qOnBCFace[m][iFace];
                        qNode[m][point] += faceQtemp * weightTemp;
                    }            
                    tNode[0][point] += qOnBCFace[nEquation][iFace] * weightTemp;

                    nodeWeight[point] += weightTemp;                
                }
            }
        }
    }
   
    //! Step4: Now the interior points.
    //!        Importantly, the symmetry points are used as interior points.
    int ** cell2NodeArray       = grid->GetCell2NodeArray();
    int *  nodeNumberOfEachCell = grid->GetNodeNumberOfEachCell();
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        int nNode = nodeNumberOfEachCell[iCell];
        for (int jNode = 0; jNode < nNode; ++ jNode)
        {
            int point = cell2NodeArray[iCell][jNode];
            if(nodeBC[point] != PHENGLEI::NO_BOUNDARY_CONDITION)
            {
                continue;
            }

            //! NO_BC && Symmetry && Interface?
            RDouble dx = x[point] - xcc[iCell];
            RDouble dy = y[point] - ycc[iCell];
            RDouble dz = z[point] - zcc[iCell];
            RDouble dist = DISTANCE(dx, dy, dz);
            RDouble weightTemp = 1.0 / dist;                

            for (int m = 0; m < nEquation; ++ m)
            {
                RDouble cellQtemp = q[m][iCell];
                qNode[m][point] += cellQtemp * weightTemp;
            }            
            tNode[0][point] += t[0][iCell] * weightTemp;

            nodeWeight[point] += weightTemp; 
        }
    }

    //! The velocity modify has not been considered on the symmetry!

    //! Change the number of points divided by the interpoint to the number of points before partition,because it is a face loop,each point is calculated twice,so the number of points needs to be multiplied by 2.
    InterpointInformation *interPointInfor = grid->GetInterpointInfo();
    if (interPointInfor)
    {
        int numberOfInterpoints = interPointInfor->GetNumberOfInterpoints();
        int *interPoint2GlobalPoint = interPointInfor->GetInterPoint2GlobalPoint();
        int *labelOfInterPoint = interPointInfor->GetLabelOfInterPoint();

        for (int iPoint = 0; iPoint < numberOfInterpoints; ++ iPoint)
        {
            int globalPoint = interPoint2GlobalPoint[iPoint];
            if (labelOfInterPoint[iPoint] != 1)
            {
                for (int m = 0; m < nEquation; ++ m)
                {
                    qNode[m][globalPoint] = 0;
                }
                tNode[0][globalPoint] = 0;
            }
        }

        RDouble **qInterPoint = reinterpret_cast<RDouble **>(grid->GetDataPtr("qInterPoint"));
        RDouble **tInterPoint = reinterpret_cast<RDouble **>(grid->GetDataPtr("tInterPoint"));
        for (int iPoint = 0; iPoint < numberOfInterpoints; ++ iPoint)
        {
            for (int m = 0; m < nEquation; m++)
            {
                qInterPoint[m][iPoint] =0;
            }
            tInterPoint[0][iPoint] =0;
        }
    }

    DelPointer2(qOnBCFace);
    DelPointer(nodeBC);
}


void NSSolverUnstruct::ModifyNodeValue(Grid *grid_in)
{
    UnstructGrid *grid = UnstructGridCast(grid_in);
    int nTotalNode = grid->GetNTotalNode();

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    RDouble     **qnode = reinterpret_cast<RDouble **>(grid->GetDataPtr("qnode"));
    RDouble     **tnode = reinterpret_cast<RDouble **>(grid->GetDataPtr("tnode"));
    RDouble *nodeWeight = reinterpret_cast<RDouble  *>(grid->GetDataPtr("nodeWeight"));

    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        if (nodeWeight[iNode] > SMALL)
        {
            for (int m = 0; m < nEquation; m++)
            {
                qnode[m][iNode] /= nodeWeight[iNode];
            }
            tnode[0][iNode] /= nodeWeight[iNode];
        }
        else
        {
            for (int m = 0; m < nEquation; m++)
            {
                qnode[m][iNode] = 0.0;
            }
            tnode[0][iNode] = 0.0;
        }
    }
}

RDouble* NSSolverUnstruct::CompSoundSpeed(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    RDouble **q   = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    RDouble *gama = reinterpret_cast< RDouble * > (grid->GetDataPtr("gama"));
    RDouble *c = new RDouble[nTotal];

    using namespace IDX;

    for (int iCell = 0; iCell < nTotal; ++ iCell)
    {
        c[iCell] = sqrt(gama[iCell] * q[IP][iCell] / q[IR][iCell]);
    }

    return c;
}

RDouble** NSSolverUnstruct::ComputeDimensionalVariables(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    RDouble **q   = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    RDouble **t   = reinterpret_cast< RDouble ** > (grid->GetDataPtr("t"));

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    RDouble refDensity = parameters->GetRefDimensionalDensity();
    RDouble refTemperature = parameters->GetRefDimensionalTemperature();
    RDouble refVelocity = parameters->GetRefMachNumber() * parameters->GetRefDimensionalSonicSpeed();
    RDouble refPressure = refDensity * refVelocity * refVelocity;

    RDouble **dimensionalVariables   = NewPointer2<RDouble>(6, nTotal);

    for (int iCell = 0; iCell < nTotal; ++ iCell)
    {
        //! Compute the dimensional density.
        dimensionalVariables[0][iCell] = refDensity * q[IR][iCell];

        //! Compute the dimensional velocity.
        dimensionalVariables[1][iCell] = refVelocity * q[IU][iCell];
        dimensionalVariables[2][iCell] = refVelocity * q[IV][iCell];
        dimensionalVariables[3][iCell] = refVelocity * q[IW][iCell];

        //! Compute the dimensional pressure.
        dimensionalVariables[4][iCell] = refPressure * q[IP][iCell];

        //! Compute the dimensional temperature.
        dimensionalVariables[5][iCell] = refTemperature * t[ITT][iCell];
    }

    return dimensionalVariables;
}

RDouble* NSSolverUnstruct::CompMachNumber(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    RDouble **q   = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    RDouble *gama = reinterpret_cast< RDouble * > (grid->GetDataPtr("gama"));
    RDouble *mach = new RDouble[nTotal];

    using namespace IDX;

    RDouble rm,um,vm,wm,pm,c2,v2;
    for (int iCell = 0; iCell < nTotal; ++ iCell)
    {
        rm = q[IR][iCell];
        um = q[IU][iCell];
        vm = q[IV][iCell];
        wm = q[IW][iCell];
        pm = q[IP][iCell];
        c2 = abs(gama[iCell] * pm / rm);
        v2 = um * um + vm * vm + wm * wm;

        mach[iCell] = sqrt(v2 / c2);
    }
    return mach;
}

RDouble * NSSolverUnstruct::ComputeVorticitybyQCriteria(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;
    RDouble * vorticity = new RDouble [ nTotal ];

    if(!WantVisualField(grid))
    {
        //! Even if do not visual field vorticity, 
        //! it would be dumped on boundary, to compatible with the field variable.
        PHSPACE::SetField(vorticity, -1.0e10, nTotal);
        return vorticity;
    }

    RDouble **dqdx = reinterpret_cast <RDouble **> (gridIn->GetDataPtr("gradPrimtiveVarX"));
    RDouble **dqdy = reinterpret_cast <RDouble **> (gridIn->GetDataPtr("gradPrimtiveVarY"));
    RDouble **dqdz = reinterpret_cast <RDouble **> (gridIn->GetDataPtr("gradPrimtiveVarZ"));

    using namespace IDX;

    RDouble * dudxField = dqdx[ IU ];
    RDouble * dudyField = dqdy[ IU ];
    RDouble * dudzField = dqdz[ IU ];

    RDouble * dvdxField = dqdx[ IV ];
    RDouble * dvdyField = dqdy[ IV ];
    RDouble * dvdzField = dqdz[ IV ];

    RDouble * dwdxField = dqdx[ IW ];
    RDouble * dwdyField = dqdy[ IW ];
    RDouble * dwdzField = dqdz[ IW ];

    for (int iCell = 0; iCell < nTotal; ++ iCell) 
    {
        RDouble dudx = dudxField[ iCell ];
        RDouble dudy = dudyField[ iCell ];
        RDouble dudz = dudzField[ iCell ];

        RDouble dvdx = dvdxField[ iCell ];
        RDouble dvdy = dvdyField[ iCell ];
        RDouble dvdz = dvdzField[ iCell ];

        RDouble dwdx = dwdxField[ iCell ];
        RDouble dwdy = dwdyField[ iCell ];
        RDouble dwdz = dwdzField[ iCell ];

        RDouble q1 = dvdy * dwdz - dwdy * dvdz;
        RDouble q2 = dudx * dwdz - dwdx * dudz;
        RDouble q3 = dudx * dvdy - dvdx * dudy;

        vorticity[ iCell ] = q1 + q2 + q3; 
    }

    return vorticity;
}

void NSSolverUnstruct::ComputeVorticitybyQCriteria(Grid *gridIn, RDouble * vorticity_x, RDouble * vorticity_y, RDouble * vorticity_z, RDouble * vorticityMagnitude,  RDouble * strain_rate,  RDouble * Q_criteria)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    RDouble **dqdx = reinterpret_cast <RDouble **> (gridIn->GetDataPtr("gradPrimtiveVarX"));
    RDouble **dqdy = reinterpret_cast <RDouble **> (gridIn->GetDataPtr("gradPrimtiveVarY"));
    RDouble **dqdz = reinterpret_cast <RDouble **> (gridIn->GetDataPtr("gradPrimtiveVarZ"));
    using namespace IDX;

    RDouble * dudxField = dqdx[ IU ];
    RDouble * dudyField = dqdy[ IU ];
    RDouble * dudzField = dqdz[ IU ];

    RDouble * dvdxField = dqdx[ IV ];
    RDouble * dvdyField = dqdy[ IV ];
    RDouble * dvdzField = dqdz[ IV ];

    RDouble * dwdxField = dqdx[ IW ];
    RDouble * dwdyField = dqdy[ IW ];
    RDouble * dwdzField = dqdz[ IW ];

    RDouble Omx,Omy,Omz,Om;
    RDouble s11, s12, s13, s22, s23, s33, S;

    for (int iCell = 0; iCell < nTotal; ++ iCell) 
    {
        RDouble dudx = dudxField[ iCell ];
        RDouble dudy = dudyField[ iCell ];
        RDouble dudz = dudzField[ iCell ];

        RDouble dvdx = dvdxField[ iCell ];
        RDouble dvdy = dvdyField[ iCell ];
        RDouble dvdz = dvdzField[ iCell ];

        RDouble dwdx = dwdxField[ iCell ];
        RDouble dwdy = dwdyField[ iCell ];
        RDouble dwdz = dwdzField[ iCell ];


        s11 = dudx; s22 = dvdy; s33 = dwdz;
        s12 = half * (dudy + dvdx); s13 = half * (dudz + dwdx); s23 = half * (dvdz + dwdy);
        S = two * (s11*s11 + s22*s22 + s33*s33 + two*(s12*s12 + s13*s13 + s23*s23));
        S = sqrt (S);
        strain_rate[iCell] = S;

        Omx = dwdy - dvdz; 
        Omy = dudz - dwdx;
        Omz = dvdx - dudy;
        Om = DISTANCE(Omx, Omy, Omz);

        vorticity_x[iCell] = Omx; 
        vorticity_y[iCell] = Omy;
        vorticity_z[iCell] = Omz;
        vorticityMagnitude[iCell] = Om; 

        Q_criteria[iCell] = fourth * (Om * Om - S * S);
    }
}

RDouble * NSSolverUnstruct::ComputeCp(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;
    RDouble ** q   = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    RDouble *primitiveVarFarfield = parameters->GetPrimitiveVarFarfield();

    using namespace IDX;
    RDouble poo = primitiveVarFarfield[IP];

    //! Pressure cofficients: cp = [p - 1.0 / (gama * mach * mach)] / 2.0 = (p - poo) / 2.0.
    RDouble * cp = new RDouble [nTotal];
    for (int iCell = 0; iCell < nTotal; ++ iCell)
    {
        cp[iCell] = 2.0 * (q[IP][iCell] - poo);
    }

    return cp;
}

RDouble * NSSolverUnstruct::ComputeCpAverage(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;
    RDouble **qAverage = reinterpret_cast< RDouble ** > (grid->GetDataPtr("qAverage"));

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    RDouble *primitiveVarFarfield = parameters->GetPrimitiveVarFarfield();

    using namespace IDX;
    RDouble poo = primitiveVarFarfield[IP];

    //! Pressure cofficients: cp = [p - 1.0 / (gama * mach * mach)] / 2.0 = (p - poo) / 2.0.
    RDouble * cp = new RDouble [nTotal];
    for (int iCell = 0; iCell < nTotal; ++ iCell)
    {
        cp[iCell] = 2.0 * (qAverage[IP][iCell] - poo);
    }

    return cp;
}

void NSSolverUnstruct::GetQlQr(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *leftCellofFace = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    RDouble **q = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));

    RDouble **qL = faceProxy->GetQL();
    RDouble **qR = faceProxy->GetQR();

    //! iFace is the global ID, jFace is the local ID in this section.
    for (int iFace = localStart; iFace < localEnd; ++ iFace)
    {
        int le = leftCellofFace [iFace];
        int re = rightCellofFace[iFace];
        int jFace  = iFace - localStart;

        for (int m = 0; m < nEquation; ++ m)
        {
            qL[m][jFace] = q[m][le];
            qR[m][jFace] = q[m][re];
        }
    }
}

void NSSolverUnstruct::GetGamaLR(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *leftCellofFace = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();

    RDouble *gama = reinterpret_cast<RDouble *> (grid->GetDataPtr("gama"));

    RDouble *gamal = faceProxy->GetGamaL();
    RDouble *gamar = faceProxy->GetGamaR();

    //! iFace is the global ID, jFace is the local ID in this section.
    for (int iFace = localStart; iFace < localEnd; ++ iFace)
    {
        int le = leftCellofFace [iFace];
        int re = rightCellofFace[iFace];
        int jFace  = iFace - localStart;

        gamal[jFace] = gama[le];
        gamar[jFace] = gama[re];
    }
}

void NSSolverUnstruct::GetPressureFactorLR(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *leftCellofFace = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();

    RDouble *rtem = reinterpret_cast <RDouble *> (grid->GetDataPtr("rtem"));

    RDouble *pressureCoeffL = faceProxy->GetPressureCoefficientL();
    RDouble *pressureCoeffR = faceProxy->GetPressureCoefficientR();

    

    for (int iFace = localStart; iFace < localEnd; ++iFace)
    {
        int le = leftCellofFace[iFace];
        int re = rightCellofFace[iFace];
        int jFace = iFace - localStart;

        pressureCoeffL[jFace] = rtem[le];
        pressureCoeffR[jFace] = rtem[re];
    }

}

void NSSolverUnstruct::ZeroResiduals(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    RDouble **res = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res"));

    for (int m = 0; m < nEquation; ++ m)
    {
        PHSPACE::SetField(res[m], 0.0, nTotal);
    }
}

void NSSolverUnstruct::ComputeMinTimeStep(Grid *gridIn, RDouble &minDt, RDouble &maxDt)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    RDouble *dt = reinterpret_cast< RDouble *  > (grid->GetDataPtr("dt"));

    minDt = LARGE;
    maxDt = - LARGE;
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        minDt = MIN(minDt, dt[iCell]);
        maxDt = MAX(maxDt, dt[iCell]);
    }
}

void NSSolverUnstruct::ReduceMaxTimeStep(Grid *gridIn, RDouble globalMinDt)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    RDouble *volume = grid->GetCellVolume();
    RDouble *dt = reinterpret_cast< RDouble *  > (grid->GetDataPtr("dt"));

    RDouble ktmax = GlobalDataBase::GetDoubleParaFromDB("ktmax");

    if (ktmax > 0)
    {
        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            dt[iCell] = MIN(dt[iCell], ktmax * globalMinDt);
        }
    }

    //! The time step is divided by the volume, so it is need not divided other else.
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        dt[iCell] /= volume[iCell];
    }

    grid->SetGhostCell(dt);
}

//! Inviscid spectral radius = 0.5 * (|V.n| + c) * S.
//! Referenced in R. F. Chen && Z. J. Wang, AIAA Journal 2000, 38(12).
void NSSolverUnstruct::SpectrumRadiusInviscid(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();

    RDouble *xfn  = grid->GetFaceNormalX();
    RDouble *yfn  = grid->GetFaceNormalY();
    RDouble *zfn  = grid->GetFaceNormalZ();
    RDouble *area = grid->GetFaceArea();

    int nTotalCell = grid->GetNTotalCell();
    int nTotalFace = grid->GetNTotalFace();

    RDouble **q   = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble *gama = reinterpret_cast<RDouble *> (grid->GetDataPtr("gama"));
    RDouble *invSpectrum = reinterpret_cast<RDouble *> (grid->GetDataPtr("invSpectralRadius"));

    RDouble *rho, *u, *v, *w, *p;

    rho = q[IDX::IR];
    u   = q[IDX::IU];
    v   = q[IDX::IV];
    w   = q[IDX::IW];
    p   = q[IDX::IP];

    PHSPACE::SetField(invSpectrum, 0.0, nTotalCell);

    for(int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        int le = leftCellofFace [iFace];
        int re = rightCellofFace[iFace];

        RDouble rl  = rho[le];
        RDouble ul  = u[le];
        RDouble vl  = v[le];
        RDouble wl  = w[le];
        RDouble pl  = p[le];
        RDouble vnl = xfn[iFace] * ul + yfn[iFace] * vl + zfn[iFace] * wl;

        RDouble rr  = rho[re];
        RDouble ur  = u[re];
        RDouble vr  = v[re];
        RDouble wr  = w[re];
        RDouble pr  = p[re];
        RDouble vnr = xfn[iFace] * ur + yfn[iFace] * vr + zfn[iFace] * wr;

        RDouble specificHeatRatio = half * (gama[le] + gama[re]);
        RDouble pm = half * (pl + pr);
        RDouble rm = half * (rl + rr);
        RDouble vn = half * (vnl + vnr);
        RDouble cm = sqrt(specificHeatRatio * pm / rm);

        RDouble inviscidSpectrumRadius = half * area[iFace] * (ABS(vn) + cm);

        invSpectrum[le] += inviscidSpectrumRadius;
        if(re < nTotalCell)
        {
            //! The re cell is the interior cell.
            invSpectrum[re] += inviscidSpectrumRadius;
        }
    }
}

void NSSolverUnstruct::SpectrumRadiusViscous(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();

    RDouble *xfn  = grid->GetFaceNormalX();
    RDouble *yfn  = grid->GetFaceNormalY();
    RDouble *zfn  = grid->GetFaceNormalZ();
    RDouble *area = grid->GetFaceArea();
    RDouble *volume = grid->GetCellVolume();

    int nTotalCell = grid->GetNTotalCell();
    int nTotalFace = grid->GetNTotalFace();

    RDouble **q   = reinterpret_cast< RDouble **> (grid->GetDataPtr("q"));
    RDouble *gama = reinterpret_cast< RDouble * > (grid->GetDataPtr("gama"));
    RDouble *visSpectrum = reinterpret_cast< RDouble *> (grid->GetDataPtr("visSpectralRadius"));

    RDouble *rho, *u, *v, *w, *p;
    rho = q[IDX::IR];
    u   = q[IDX::IU];
    v   = q[IDX::IV];
    w   = q[IDX::IW];
    p   = q[IDX::IP];

    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    RDouble *viscousLaminar = reinterpret_cast< RDouble * > (grid->GetDataPtr("visl"));
    RDouble *viscousTurbulence = reinterpret_cast< RDouble * > (grid->GetDataPtr("vist"));

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    RDouble refReNumber = parameters->GetRefReNumber();

    RDouble prandtlLaminar    = parameters->GetPrandtlLaminar();    //! Default 0.72
    RDouble prandtlTurbulence = parameters->GetPrandtlTurbulence(); //! Default 0.9

    const RDouble fourthThird = 4.0 / 3.0;
    RDouble oprl = 1.0 / prandtlLaminar;
    RDouble oprt = 1.0 / prandtlTurbulence;

    PHSPACE::SetField(visSpectrum, 0.0, nTotalCell);

    int visSpetrumMethod = 2;
    if(visSpetrumMethod == 1)
    {
        for (int iFace = 0; iFace < nTotalFace; ++ iFace)
        {
            int le = leftCellofFace [iFace];
            int re = rightCellofFace[iFace];

            RDouble visLaminar = half * (viscousLaminar[le] + viscousLaminar[re]);
            RDouble visTurb    = half * (viscousTurbulence[le] + viscousTurbulence[re]);
            RDouble density    = half * (rho[le] + rho[re]);

            RDouble specificHeatRatio = half * (gama[le] + gama[re]);
            RDouble coefOfSpectrum1 = fourthThird * (visLaminar + visTurb);
            RDouble coefOfSpectrum2 = specificHeatRatio * (visLaminar * oprl + visTurb * oprt);

            //! Viscous spectral radius type 2.
            //! Referenced in Wang Gang, Chinese Journal of Aeronautics, 2012 (25): 33-41, eq. 13.
            RDouble coefOfSpectrum  = two * PHSPACE::MAX(coefOfSpectrum1, coefOfSpectrum2) / (refReNumber * density);

            RDouble faceArea = area[iFace];
            RDouble faceAera2 = PHSPACE::SQR(faceArea);

            RDouble viscousSpectrumRadius = faceAera2 * coefOfSpectrum;

            visSpectrum[le] += viscousSpectrumRadius / volume[le];
            if(re < nTotalCell)
            {
                visSpectrum[re] += viscousSpectrumRadius / volume[re];
            }
        }
    }
    else if (visSpetrumMethod == 2)
    {
        for (int iFace = 0; iFace < nTotalFace; ++ iFace)
        {
            int le = leftCellofFace [iFace];
            int re = rightCellofFace[iFace];

            RDouble dx  = xcc[re] - xcc[le];
            RDouble dy  = ycc[re] - ycc[le];
            RDouble dz  = zcc[re] - zcc[le];
            RDouble ds  = ABS(xfn[iFace] * dx + yfn[iFace] * dy + zfn[iFace] * dz);

            RDouble visLaminar = half * (viscousLaminar[le] + viscousLaminar[re]);
            RDouble visTurb    = half * (viscousTurbulence[le] + viscousTurbulence[re]);
            RDouble density    = half * (rho[le] + rho[re]);
            RDouble viscosity  = visLaminar + visTurb;
            RDouble faceArea  = area[iFace];

            //! Viscous spectral radius = 2 * miu/(rho.|n.dr|) * S.
            //! Referenced in R. F. Chen && Z. J. Wang, AIAA Journal 2000, 38(12).
            RDouble viscousSpectrumRadius  = 2.0 * viscosity / (density * ds * refReNumber + SMALL);
            viscousSpectrumRadius *= half * faceArea;

            visSpectrum[le] += viscousSpectrumRadius;
            if(re < nTotalCell)
            {
                visSpectrum[re] += viscousSpectrumRadius;
            }
        }
    }
}

void NSSolverUnstruct::PreconInvSpectrumRadius(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int      * leftCellofFace = grid->GetLeftCellOfFace();
    int      * rightCellofFace = grid->GetRightCellOfFace();

    RDouble * xfn  = grid->GetFaceNormalX();
    RDouble * yfn  = grid->GetFaceNormalY();
    RDouble * zfn  = grid->GetFaceNormalZ();
    RDouble * area = grid->GetFaceArea();

    int nTotalCell = grid->GetNTotalCell();
    int nTotalFace = grid->GetNTotalFace();

    RDouble **q          = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble *gama        = reinterpret_cast<RDouble  *> (grid->GetDataPtr("gama"));
    RDouble *invSpectrum = reinterpret_cast<RDouble  *> (grid->GetDataPtr("invSpectralRadius"));

    RDouble refMachNumber = GetControlParameters()->GetRefMachNumber();

    int le, re;

    RDouble *rho, *u, *v, *w, *p;

    using namespace IDX;

    rho = q[IR];
    u   = q[IU];
    v   = q[IV];
    w   = q[IW];
    p   = q[IP];

    PHSPACE::SetField(invSpectrum, 0.0, nTotalCell);

    RDouble half = 0.5;
    RDouble mach2 = refMachNumber * refMachNumber;
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        le = leftCellofFace [iFace];
        re = rightCellofFace[iFace];

        RDouble rl  = rho[le];
        RDouble ul  = u[le];
        RDouble vl  = v[le];
        RDouble wl  = w[le];
        RDouble pl  = p[le];
        RDouble vnl = xfn[iFace] * ul + yfn[iFace] * vl + zfn[iFace] * wl;
        RDouble v2l = ul * ul + vl * vl + wl * wl;

        RDouble rr  = rho[re];
        RDouble ur  = u[re];
        RDouble vr  = v[re];
        RDouble wr  = w[re];
        RDouble pr  = p[re];
        RDouble vnr = xfn[iFace] * ur + yfn[iFace] * vr + zfn[iFace] * wr;
        RDouble v2r = ur * ur + vr * vr + wr * wr;

        RDouble specificHeatRatio = half * (gama[le] + gama[re]);
        RDouble pm  = half * (pl + pr);
        RDouble rm  = half * (rl + rr);
        RDouble cm  = sqrt(specificHeatRatio * pm / rm);
        RDouble vn  = half * (vnl + vnr);
        RDouble v2  = half * (v2l + v2r);
        RDouble c2  = cm * cm;

        RDouble beta = MIN(MAX(v2 / c2,3.0 * mach2),1.0);
        RDouble vnew  = half * vn * (1 + beta);
        RDouble cnew  = half * sqrt(((beta - 1) * vn) * ((beta - 1) * vn) + 4 * beta * c2);
        RDouble inviscidSpectrumRadius = half * area[iFace] * (ABS(vnew) + cnew);

        invSpectrum[le] += inviscidSpectrumRadius;
        if(re < nTotalCell)
        {
            invSpectrum[re] += inviscidSpectrumRadius;
        }        
    }
}

void NSSolverUnstruct::SpectrumRadiusChemical(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    RDouble *vol = grid->GetCellVolume();

    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    RDouble **t = reinterpret_cast< RDouble ** > (grid->GetDataPtr("t"));
    RDouble **chemSpectrum = reinterpret_cast< RDouble ** > (grid->GetDataPtr("chemSpectralRadius"));

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nm = parameters->GetNSEquationNumber();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nSpeciesNumber = parameters->GetNumberOfSpecies();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    PHSPACE::SetField(chemSpectrum, nSpeciesNumber, nTotalCell, 0.0);

    RDouble **srs = NewPointer2<RDouble>(nSpeciesNumber, nTotalCell);

    using namespace IDX;

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        gas->ChemicalSpectrumRadius(q, t, vol, srs, iCell);

        for (int m = 0; m < nSpeciesNumber - 1; ++ m)
        {
            chemSpectrum[m][iCell] += srs[m][iCell];
        }
    }

    DelPointer2(srs);
}

void NSSolverUnstruct::LocalTimeStep(Grid *gridIn)
{

    UnstructGrid *grid = UnstructGridCast(gridIn);
    RDouble *volume  = grid->GetCellVolume();

    int nTotalCell = grid->GetNTotalCell();

    RDouble *dt   = reinterpret_cast< RDouble *  > (grid->GetDataPtr("dt"));
    RDouble *invSpectralRadius = reinterpret_cast<RDouble *> (grid->GetDataPtr("invSpectralRadius"));
    RDouble *visSpectralRadius = reinterpret_cast<RDouble *> (grid->GetDataPtr("visSpectralRadius"));

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    RDouble cfl = ComputeCFL();

    int isInMGIniting = GlobalDataBase::GetIntParaFromDB("isInMGIniting");
    if (!grid->IsFinestGrid() && !isInMGIniting)
    {
        RDouble mgCFLScale = parameters->GetMgCFLScale();
        cfl *= mgCFLScale;
    }

    //! Set dt to zero.
    PHSPACE::SetField(dt, 0.0, nTotalCell);

    //! Inviscid time step.
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        dt[iCell] = cfl * volume[iCell] / invSpectralRadius[iCell];
    }

    int isViscous = parameters->IsViscous();

    if (isViscous)
    {
        //! Viscous time step.
        RDouble *dtViscous = new RDouble [nTotalCell];
        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            dtViscous[iCell] = cfl * volume[iCell] / visSpectralRadius[iCell];
        }

        //! Consider the viscous time step.
        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            dt[iCell] = dt[iCell] * dtViscous[iCell] / (dt[iCell] + dtViscous[iCell]);
        }

        delete [] dtViscous;
    }
}

void NSSolverUnstruct::LocalTimeStep_Prince(Grid* gridIn)
{
    UnstructGrid* grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nTotalFace = grid->GetNTotalFace();
    int nBoundFace = grid->GetNBoundFace();
    int* leftCellofFace = grid->GetLeftCellOfFace();
    int* rightCellofFace = grid->GetRightCellOfFace();
    RDouble* xfc = grid->GetFaceCenterX();
    RDouble* yfc = grid->GetFaceCenterY();
    RDouble* zfc = grid->GetFaceCenterZ();
    RDouble* xcc = grid->GetCellCenterX();
    RDouble* ycc = grid->GetCellCenterY();
    RDouble* zcc = grid->GetCellCenterZ();
    RDouble* xfn = grid->GetFaceNormalX();
    RDouble* yfn = grid->GetFaceNormalY();
    RDouble* zfn = grid->GetFaceNormalZ();
    RDouble* area = grid->GetFaceArea();
    RDouble* vol = grid->GetCellVolume();

    RDouble** q = reinterpret_cast<RDouble**> (grid->GetDataPtr("q"));
    RDouble* visl = reinterpret_cast<RDouble*> (grid->GetDataPtr("visl"));
    RDouble* vist = reinterpret_cast<RDouble*> (grid->GetDataPtr("vist"));
    RDouble* dt = reinterpret_cast<RDouble*> (grid->GetDataPtr("dt"));
    //! Set dt to LARGE.
    PHSPACE::SetField(dt, LARGE, nTotalCell);

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    RDouble refGama = parameters->GetRefGama();
    RDouble isViscous = parameters->IsViscous();
    RDouble prandtlLaminar = parameters->GetPrandtlLaminar();
    RDouble prandtlTurbulence = parameters->GetPrandtlTurbulence();

    int le, re;
    RDouble c, Vn, eigv, dn;
    RDouble Coeff = 4.0;
    //! For boundary faces
    for (int iFace = 0; iFace < nBoundFace; iFace++)
    {
        le = leftCellofFace[iFace];

        c = sqrt(refGama * q[IP][le] / q[IR][le]);
        dn = fabs((xfc[iFace] - xcc[le]) * xfn[iFace] + (yfc[iFace] - ycc[le]) * yfn[iFace] + (zfc[iFace] - zcc[le]) * zfn[iFace]);
        Vn = ABS(q[IU][le] * xfn[iFace] + q[IV][le] * yfn[iFace] + q[IW][le] * zfn[iFace]);
        eigv = Vn + c;

        if (isViscous)
        {
            RDouble muoopr = visl[le] / prandtlLaminar + vist[le] / prandtlTurbulence;
            eigv += Coeff * refGama / q[IR][le] * muoopr * area[iFace] / vol[le];
        }
        dt[le] = MIN(dt[le], dn / eigv);
    }
    //! For interior faces
    for (int iFace = nBoundFace; iFace < nTotalFace; iFace++)
    {
        le = leftCellofFace[iFace];
        re = rightCellofFace[iFace];

        c = sqrt(refGama * q[IP][le] / q[IR][le]);
        dn = fabs((xfc[iFace] - xcc[le]) * xfn[iFace] + (yfc[iFace] - ycc[le]) * yfn[iFace] + (zfc[iFace] - zcc[le]) * zfn[iFace]);
        Vn = ABS(q[IU][le] * xfn[iFace] + q[IV][le] * yfn[iFace] + q[IW][le] * zfn[iFace]);
        eigv = Vn + c;
        if (isViscous)
        {
            RDouble muoopr = visl[le] / prandtlLaminar + vist[le] / prandtlTurbulence;
            eigv += Coeff * refGama / q[IR][le] * muoopr * area[iFace] / vol[le];
        }
        dt[le] = MIN(dt[le], dn / eigv);

        c = sqrt(refGama * q[IP][re] / q[IR][re]);
        dn = fabs((xfc[iFace] - xcc[re]) * xfn[iFace] + (yfc[iFace] - ycc[re]) * yfn[iFace] + (zfc[iFace] - zcc[re]) * zfn[iFace]);
        Vn = ABS(q[IU][re] * xfn[iFace] + q[IV][re] * yfn[iFace] + q[IW][re] * zfn[iFace]);
        eigv = Vn + c;
        if (isViscous)
        {
            RDouble muoopr = visl[re] / prandtlLaminar + vist[re] / prandtlTurbulence;
            eigv += Coeff * refGama / q[IR][re] * muoopr * area[iFace] / vol[re];
        }
        dt[re] = MIN(dt[re], dn / eigv);
    }

    RDouble cfl = ComputeCFL();
    for (int iCell = 0; iCell < nTotalCell; ++iCell)
    {
        dt[iCell] = dt[iCell] * cfl;
    }
}

void NSSolverUnstruct::GlobalTimeStep(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    RDouble * vol = grid->GetCellVolume();
    RDouble * dt  = reinterpret_cast< RDouble * > (grid->GetDataPtr("dt" ));

    int nTotalCell = grid->GetNTotalCell();

    RDouble dtau;
    GlobalDataBase::GetData("dtau", &dtau, PHDOUBLE, 1);

    //! Set dt to zero.
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        dt[iCell] = dtau / vol[iCell];
    }

    grid->SetGhostCell(dt);
    RDouble dtmin = dtau;
    GlobalDataBase::UpdateData("dtmin", &dtmin, PHDOUBLE, 1);
}

void NSSolverUnstruct::LocalGlobalTimeStep(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int * leftCellofFace = grid->GetLeftCellOfFace();
    int * rightCellofFace = grid->GetRightCellOfFace();

    RDouble * xfn  = grid->GetFaceNormalX();
    RDouble * yfn  = grid->GetFaceNormalY();
    RDouble * zfn  = grid->GetFaceNormalZ();
    RDouble * area = grid->GetFaceArea();
    RDouble * vgn  = grid->GetFaceNormalVelocity();
    RDouble * vol  = grid->GetCellVolume();

    int nTotalCell = grid->GetNTotalCell();
    int nTotalFace = grid->GetNTotalFace();
    int nBoundFace = grid->GetNBoundFace();

    RDouble **q   = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    RDouble *gama = reinterpret_cast< RDouble *  > (grid->GetDataPtr("gama"));
    RDouble *dt   = reinterpret_cast< RDouble *  > (grid->GetDataPtr("dt"));

    RDouble cfl = ComputeCFL();

    int le, re;
    RDouble um, vm, wm, cm;

    RDouble *rho, *u, *v, *w, *p;

    using namespace IDX;

    rho = q[IR];
    u   = q[IU];
    v   = q[IV];
    w   = q[IW];
    p   = q[IP];

    //! Set dt to zero.
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        dt[iCell] = 0.0;
    }

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        le = leftCellofFace[ iFace ];

        um      = u[le];
        vm      = v[le];
        wm      = w[le];
        cm      = sqrt(gama[le] * p[le] / rho[le]);
        dt[le] += (ABS(xfn[iFace] * um + yfn[iFace] * vm + zfn[iFace] * wm - vgn[iFace]) + cm) * area[iFace];
    }

    //! For interior faces.
    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        le = leftCellofFace[ iFace ];
        re = rightCellofFace[ iFace ];

        um      = u[le];
        vm      = v[le];
        wm      = w[le];
        cm      = sqrt(gama[le] * p[le] / rho[le]);
        dt[le] += (ABS(xfn[iFace] * um + yfn[iFace] * vm + zfn[iFace] * wm - vgn[iFace]) + cm) * area[iFace];

        um      = u[re];
        vm      = v[re];
        wm      = w[re];
        cm      = sqrt(gama[re] * p[re] / rho[re]);
        dt[re] += (ABS(xfn[iFace] * um + yfn[iFace] * vm + zfn[iFace] * wm - vgn[iFace]) + cm) * area[iFace];
    }

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        dt[iCell] = cfl * vol[iCell] / dt[iCell];
    }

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int viscousType = parameters->GetViscousType();

    if (viscousType > INVISCID)
    {
        RDouble * dtv;
        dtv = new RDouble [ nTotalCell ];

        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            dtv[iCell] = 0.0;
        }

        for (int iFace = 0; iFace < nBoundFace; ++ iFace)
        {
            le = leftCellofFace[ iFace ];
            re = rightCellofFace[ iFace ];

            dtv[le] += area[iFace] * area[iFace];
        }

        //! For interior faces.
        for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
        {
            le = leftCellofFace[ iFace ];
            re = rightCellofFace[ iFace ];

            dtv[le] += area[iFace] * area[iFace];
            dtv[re] += area[iFace] * area[iFace];
        }

        RDouble * visl = reinterpret_cast< RDouble * > (grid->GetDataPtr("visl"));
        RDouble * vist = reinterpret_cast< RDouble * > (grid->GetDataPtr("vist"));

        RDouble refReNumber  = parameters->GetRefReNumber();
        RDouble prandtlLaminar    = parameters->GetPrandtlLaminar();
        RDouble prandtlTurbulence = parameters->GetPrandtlTurbulence();

        RDouble coef1,coef2,coef,vis_l,vis_t;
        RDouble csrvist = 1.0;

        const RDouble foth = 4.0 / 3.0;

        RDouble oprl = 1.0 / prandtlLaminar;
        RDouble oprt = 1.0 / prandtlTurbulence;

        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            vis_l = visl[iCell];
            vis_t = csrvist * vist[iCell];
            coef1 = foth * (vis_l + vis_t);
            coef2 = gama[iCell] * (vis_l * oprl + vis_t * oprt);
            coef  = MAX(coef1, coef2);
            dtv[iCell] = cfl * rho[iCell] * refReNumber * vol[iCell] * vol[iCell] / (coef * dtv[iCell]);
        }

        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            dt[iCell] = dt[iCell] * dtv[iCell] / (dt[iCell] + dtv[iCell]);
        }

        delete [] dtv; 
    }

    RDouble dtmin = LARGE;
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        dtmin = MIN(dtmin, static_cast<RDouble>(dt[iCell]));
    }

    GlobalDataBase::UpdateData("dtmin", &dtmin, PHDOUBLE, 1);

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        dt[iCell] = dtmin / vol[iCell];
    }
}

void NSSolverUnstruct::FillField(Grid *gridIn, FieldProxy *field_proxy, RDouble value)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    RDouble **field = field_proxy->GetField_UNS();

    for (int m = 0; m < nEquation; ++ m)
    {
        for (int iCell = 0; iCell < nTotal; ++ iCell)
        {
            field[m][iCell] = value;
        }
    }
}

void NSSolverUnstruct::FillField(Grid *gridIn, FieldProxy *fieldTarget, FieldProxy *fieldSource)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    RDouble **fieldTar = fieldTarget->GetField_UNS();
    RDouble **fieldSrc = fieldSource->GetField_UNS();

    for (int m = 0; m < nEquation; ++ m)
    {
        for (int iCell = 0; iCell < nTotal; ++ iCell)
        {
            fieldTar[m][iCell] = fieldSrc[m][iCell];
        }
    }
}

void NSSolverUnstruct::RungeKuttaResidual(Grid *gridIn, FieldProxy *dqProxy, RDouble coef)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int      nTotalCell = grid->GetNTotalCell();

    RDouble **dq = dqProxy->GetField_UNS();
    RDouble **res = reinterpret_cast <RDouble **> (grid->GetDataPtr("res"));

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    RDouble * dt = reinterpret_cast <RDouble *> (grid->GetDataPtr("dt"));

    for (int iEquation = 0; iEquation < nEquation; ++ iEquation)
    {
        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            //! The volume has been divided in time step, so the dt is actually 'dt/vol'.
            dq[iEquation][iCell] = dt[iCell] * coef * res[iEquation][iCell];
        }
    }
}

//! Restrict from fine to coarse grid for all q.
void NSSolverUnstruct::RestrictAllQ(Grid *fgrid, Grid *cgrid)
{
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    RDouble **fq = reinterpret_cast< RDouble ** > (fgrid->GetDataPtr("q"));
    RDouble **cq = reinterpret_cast< RDouble ** > (cgrid->GetDataPtr("q"));

    for (int m = 0; m < nEquation; ++ m)
    {
        RestrictQ(fgrid, fq[m], cgrid, cq[m]);
    }

    RDouble *fgama = reinterpret_cast< RDouble * > (fgrid->GetDataPtr("gama"));
    RDouble *cgama = reinterpret_cast< RDouble * > (cgrid->GetDataPtr("gama"));

    RestrictQ(fgrid, fgama, cgrid, cgama);
    
    int viscousType = parameters->GetViscousType();
    if (viscousType > INVISCID)
    {
        RDouble * fvisl = reinterpret_cast< RDouble * > (fgrid->GetDataPtr("visl"));
        RDouble * cvisl = reinterpret_cast< RDouble * > (cgrid->GetDataPtr("visl"));
        RestrictQ(fgrid, fvisl, cgrid, cvisl);

        if (viscousType > LAMINAR)
        {
            RDouble * fvist = reinterpret_cast< RDouble * > (fgrid->GetDataPtr("vist"));
            RDouble * cvist = reinterpret_cast< RDouble * > (cgrid->GetDataPtr("vist"));
            RestrictQ(fgrid, fvist, cgrid, cvist);
        }
    }

    int isUnsteady = parameters->GetIsUnsteady();
    if (isUnsteady)
    {
        RDouble **fqn1 = reinterpret_cast< RDouble ** > (fgrid->GetDataPtr("q_unsteady_n1"));
        RDouble **fqn2 = reinterpret_cast< RDouble ** > (fgrid->GetDataPtr("q_unsteady_n2"));
        RDouble **cqn1 = reinterpret_cast< RDouble ** > (cgrid->GetDataPtr("q_unsteady_n1"));
        RDouble **cqn2 = reinterpret_cast< RDouble ** > (cgrid->GetDataPtr("q_unsteady_n2"));

        RDouble **fres   = reinterpret_cast< RDouble ** > (fgrid->GetDataPtr("res"           ));
        RDouble **fresn1 = reinterpret_cast< RDouble ** > (fgrid->GetDataPtr("res_unsteady_n1"));
        RDouble **fresn2 = reinterpret_cast< RDouble ** > (fgrid->GetDataPtr("res_unsteady_n2"));
        RDouble **cres   = reinterpret_cast< RDouble ** > (cgrid->GetDataPtr("res"           ));
        RDouble **cresn1 = reinterpret_cast< RDouble ** > (cgrid->GetDataPtr("res_unsteady_n1"));
        RDouble **cresn2 = reinterpret_cast< RDouble ** > (cgrid->GetDataPtr("res_unsteady_n2"));

        for (int m = 0; m < nEquation; ++ m)
        {
            RestrictQ(fgrid, fqn1[m], cgrid, cqn1[m]);
            RestrictQ(fgrid, fqn2[m], cgrid, cqn2[m]);

            RestrictQ(fgrid, fres[m],   cgrid, cres[m]);
            RestrictQ(fgrid, fresn1[m], cgrid, cresn1[m]);
            RestrictQ(fgrid, fresn2[m], cgrid, cresn2[m]);
        }
    }
}

void NSSolverUnstruct::LoadFlux(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nBoundFace = grid->GetNBoundFace();

    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar = parameters->GetLaminarNumber();

    RDouble **flux = faceProxy->GetFlux();

    //! Determine if there are boundary faces.
    int nMid  = localStart;
    if (localEnd <= nBoundFace)
    {
       //! If all boundary faces.
       nMid = localEnd;
    }
    else if (localStart < nBoundFace)
    {
       //! Part of them are boundary faces.
       nMid = nBoundFace;
    }

    //! 'rhs' := inviscid flux - viscid flux.
    //! 'res' is defined as '-rhs'.

    RDouble **res = reinterpret_cast<RDouble **> (grid->GetDataPtr("res"));

    //! For boundary faces, remember re is ghost cell.
    for (int iFace = localStart; iFace < nMid; ++ iFace)
    {
        int le = leftCellofFace[iFace];
        int jFace = iFace - localStart;
        for (int m = 0; m < nLaminar; ++ m)
        {
            res[m][le] -= flux[m][jFace];
        }
    }

    //! Interior faces.
    for (int iFace = nMid; iFace < localEnd; ++ iFace)
    {
        int le = leftCellofFace[iFace];
        int re = rightCellofFace[iFace];
        int jFace = iFace - localStart;

        for (int m = 0; m < nLaminar; ++ m)
        {
            res[m][le] -= flux[m][jFace];
            res[m][re] += flux[m][jFace];
        }
    }
}


void NSSolverUnstruct::InitialInviscidfaceproxyOnlyForMixGrid(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalFace = grid->GetNTotalFace();
    int nEquation = 5;

    GeomProxy *geomProxy = new GeomProxy();
    geomProxy->Create(nTotalFace);
    FillGeomProxy(grid, geomProxy, 0, nTotalFace);

    InviscidfaceProxyOnlyForMixGrid = new FaceProxy();
    InviscidfaceProxyOnlyForMixGrid->Create(nTotalFace, nEquation);
    InviscidfaceProxyOnlyForMixGrid->SetGeomProxy(geomProxy);

}

void NSSolverUnstruct::InitialViscousfaceproxyOnlyForMixGrid(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalFace = grid->GetNTotalFace();
    int nEquation = 5;

    ViscousfaceProxyOnlyForMixGrid = new FaceProxy();
    ViscousfaceProxyOnlyForMixGrid->Create(nTotalFace, nEquation);

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int numberOfSpecies = parameters->GetNumberOfSpecies();
    
    NSFaceValue *facevar = new NSFaceValue(nEquation, numberOfSpecies, nTotalFace);
    ViscousfaceProxyOnlyForMixGrid->SetNSFaceValue(facevar);

    ComputeFaceWeight(grid, ViscousfaceProxyOnlyForMixGrid, 0, nTotalFace);
}


void NSSolverUnstruct::UpdateQlQrOnlyforMixGrid(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalFace = grid->GetNTotalFace();

    int localStart = 0;
    int localEnd = nTotalFace;

    GetInviscidFaceValue(grid, InviscidfaceProxyOnlyForMixGrid, limiterOnlyForMixGrid, localStart, localEnd);
}

void NSSolverUnstruct::UpdateInviscidfluxOnlyforMixGrid(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalFace = grid->GetNTotalFace();

    int localStart = 0;
    int localEnd = nTotalFace;

    ComputeInviscidFlux(grid, InviscidfaceProxyOnlyForMixGrid, localStart, localEnd);
}

void NSSolverUnstruct::UpdateViscousfluxOnlyforMixGrid(Grid *gridIn)
{
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int viscousType = parameters->GetViscousType();
    if (viscousType == INVISCID) return;

    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalFace = grid->GetNTotalFace();

    int localStart = 0;
    int localEnd = nTotalFace;

    GetVisFaceValue(grid, ViscousfaceProxyOnlyForMixGrid, localStart, localEnd);

    ComputeVisflux(grid, ViscousfaceProxyOnlyForMixGrid, localStart, localEnd);
}

void NSSolverUnstruct::LoadGlobalfacefluxtoResOnlyforMixGrid(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nBoundFace = grid->GetNBoundFace();
    int nTotalFace = grid->GetNTotalFace();

    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar = parameters->GetLaminarNumber();

    RDouble **res = reinterpret_cast<RDouble **> (grid->GetDataPtr("res"));
    RDouble **Inviscidflux = InviscidfaceProxyOnlyForMixGrid->GetFlux();
    //! For boundary faces, remember re is ghost cell.
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int le = leftCellofFace[iFace];
        for (int m = 0; m < nLaminar; ++ m)
        {
            res[m][le] -= Inviscidflux[m][iFace];
        }
    }

    //! Interior faces.
    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        int le = leftCellofFace[iFace];
        int re = rightCellofFace[iFace];

        for (int m = 0; m < nLaminar; ++ m)
        {
            res[m][le] -= Inviscidflux[m][iFace];
            res[m][re] += Inviscidflux[m][iFace];
        }
    }

    int viscousType = parameters->GetViscousType();
    if (viscousType == INVISCID) return;

    RDouble **Viscousflux = ViscousfaceProxyOnlyForMixGrid->GetFlux();
    //! For boundary faces, remember re is ghost cell.
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int le = leftCellofFace[iFace];
        for (int m = 0; m < nLaminar; ++ m)
        {
            res[m][le] -= Viscousflux[m][iFace];
        }
    }

    //! Interior faces.
    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        int le = leftCellofFace[iFace];
        int re = rightCellofFace[iFace];

        for (int m = 0; m < nLaminar; ++ m)
        {
            res[m][le] -= Viscousflux[m][iFace];
            res[m][re] += Viscousflux[m][iFace];
        }
    }
}

void NSSolverUnstruct::ComputeViscousCoeff(Grid *gridIn)
{
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nChemical = parameters->GetChemicalFlag();
    if (nChemical == 0)
    {
        ComputeViscousCoefficientWithoutChemical(gridIn);
    }
    else
    {
        ComputeViscousCoefficient(gridIn);
    }

}

void NSSolverUnstruct::ComputeViscousCoefficient(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    RDouble nonDimensionalSutherlandTemperature = parameters->GetNonDimensionalSutherlandTemperature();
    RDouble wallTemperature = parameters->GetWallTemperature();
    int nLaminar = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nSpeciesNumber = parameters->GetNumberOfSpecies();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    RDouble *viscousLaminar = reinterpret_cast<RDouble *> (grid->GetDataPtr("visl"));
    RDouble *viscousTurbulence = reinterpret_cast <RDouble *> (grid->GetDataPtr("vist"));
    RDouble **q = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble **t = reinterpret_cast<RDouble **> (grid->GetDataPtr("t"));

    using namespace IDX;
    RDouble *primitiveVars  = new RDouble[nEquation];
    RDouble temperature, electronTemperature, viscosity;

    RDouble refTemperatureSutherland = GlobalDataBase::GetDoubleParaFromDB("tsuth");

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        temperature = t[ITT][iCell];
        for (int iEquation = 0; iEquation < nEquation; ++ iEquation)
        {
            primitiveVars[iEquation] = q[iEquation][iCell];
        }
        electronTemperature = temperature;
        if (nTemperatureModel == 2)
        {
            electronTemperature = t[ITV][iCell];
        }
        else if (nTemperatureModel == 3)
        {
            electronTemperature = t[ITE][iCell];
        }
        gas->ComputeViscosityByPrimitive(primitiveVars, temperature, electronTemperature, refTemperatureSutherland, viscosity);
        viscousLaminar[iCell] = viscosity;
    }
    if (wallTemperature <= 0.0)
    {
        grid->SetGhostCell(viscousLaminar);
    }
    else
    {
        for (int iCell = nTotalCell; iCell < nTotal; ++ iCell)
        {
            temperature = t[ITT][iCell];
            for (int iEquation = 0; iEquation < nEquation; ++ iEquation)
            {
                primitiveVars[iEquation] = q[iEquation][iCell];
            }
            electronTemperature = temperature;
            if (nTemperatureModel == 2)
            {
                electronTemperature = t[ITV][iCell];
            }
            else if (nTemperatureModel == 3)
            {
                electronTemperature = t[ITE][iCell];
            }
            gas->ComputeViscosityByPrimitive(primitiveVars, temperature, electronTemperature, refTemperatureSutherland, viscosity);
            viscousLaminar[iCell] = viscosity;
        }
    }

    delete []primitiveVars;
}

void NSSolverUnstruct::ComputeViscousCoefficientWithoutChemical(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    RDouble *viscousLaminar = reinterpret_cast<RDouble *> (grid->GetDataPtr("visl"));
    RDouble **t = reinterpret_cast<RDouble **> (grid->GetDataPtr("t"));

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    RDouble nonDimensionalSutherlandTemperature = parameters->GetNonDimensionalSutherlandTemperature();

    RDouble viscousLaminarMin;
    GlobalDataBase::GetData("visl_min", &viscousLaminarMin, PHDOUBLE, 1);

    using namespace IDX;

    for (int iCell = 0; iCell < nTotal; ++ iCell)
    {
        viscousLaminar[iCell] = t[ITT][iCell] * sqrt(t[ITT][iCell]) * (1.0 + nonDimensionalSutherlandTemperature) / (t[ITT][iCell] + nonDimensionalSutherlandTemperature);
        viscousLaminar[iCell] = MAX(viscousLaminarMin, viscousLaminar[iCell]);
    }
}

void NSSolverUnstruct::InitCGrid(Grid *fgrid_in, Grid *cgrid_in)
{
    UnstructGrid *fgrid = UnstructGridCast(fgrid_in);
    UnstructGrid *cgrid = UnstructGridCast(cgrid_in);
    int * cell2coarsegridcell  = fgrid->GetCell2CoarseGridCell();
    int fnTotalCell = fgrid->GetNTotalCell();
    int cnTCell = cgrid->GetNTotalCell();
    RDouble *cvol = cgrid->GetCellVolume();
    RDouble *fvol = fgrid->GetCellVolume();

    RDouble **fq = reinterpret_cast< RDouble ** > (fgrid->GetDataPtr("q"));
    RDouble **cq = reinterpret_cast< RDouble ** > (cgrid->GetDataPtr("q"));

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    int cn = cnTCell + cgrid->GetNBoundFace();

    for (int m = 0; m < nEquation; ++ m)
    {
        for (int iCell = 0; iCell < cn; ++ iCell)
        {
            cq[m][iCell] = 0.0;
        }
    }

    for (int m = 0; m < nEquation; ++ m)
    {
        for (int iCell = 0; iCell < fnTotalCell; ++ iCell)
        {
            cq[m][ cell2coarsegridcell[iCell] ] += fvol[iCell] * fq[m][iCell];
        }
    }
    for (int m = 0; m < nEquation; ++ m)
    {
        for (int iCell = 0; iCell < cnTCell; ++ iCell)
        {
            cq[m][iCell] /= cvol[iCell];
        }
    }
}

//! Drive individual functions to calculate inviscid fluxes.
void NSSolverUnstruct::InviscidFlux(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalFace = grid->GetNTotalFace();

    Limiter *limiter = 0;

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int limiterType = parameters->GetLimiterType();
    int isInIniting = GlobalDataBase::GetIntParaFromDB("isInIniting");
    bool isFirstOrder = (limiterType == ILMT_FIRST) || isInIniting;
    if (grid->IsFinestGrid() && !isFirstOrder)
    {
        limiter = CreateLimiter(grid);
    }

     //! Face proxy.
    FaceProxy *faceProxy = CreateFaceProxy(grid);
    faceProxy->SetGeomProxy(CreateGeomProxy(grid));

    int localStart, localEnd;
    localStart = 0;
    do 
    {
        localEnd = localStart + SEGCTION_LENGTH;
        if (localEnd > nTotalFace) 
        {
            localEnd = nTotalFace;
        }

        GetInviscidFaceValue(grid, faceProxy, limiter, localStart, localEnd);

        ComputeInviscidFlux(grid, faceProxy, localStart, localEnd);

        LoadFlux(grid, faceProxy, localStart, localEnd);

        localStart = localEnd;
    } while (localStart < nTotalFace);

    delete faceProxy;
    delete limiter;
}

void NSSolverUnstruct::CalPressureFactor(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    int nTotalCell = grid->GetNTotalCell();

    int **cell2face = grid->GetCell2Face();
    vector<int> * c2c = grid->GetCell2Cell();
    int *faceNumberOfEachCell = grid->GetFaceNumberOfEachCell();

    int *leftCellofFace = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();

    RDouble minPressureRatio = 1.3;
    RDouble maxPressureRatio = 2.6;
    RDouble pi = 3.141592653589;

    int RoeEntropyFixMethod = this->GetControlParameters()->GetRoeEntropyFixMethod();

    RDouble **q   = reinterpret_cast <RDouble **> (grid->GetDataPtr("q"));
    RDouble *rtem = reinterpret_cast <RDouble *> (grid->GetDataPtr("rtem"));

    if (RoeEntropyFixMethod == 2)
    {
        for (int iCell = 0; iCell < nTotalCell; ++iCell)
        {
            RDouble LeftPressureOfCell = zero;
            RDouble RightPressureOfCell = zero;
            rtem[iCell] = zero;
            for (int j = 0; j < faceNumberOfEachCell[iCell]; ++j)
            {
                int face = cell2face[iCell][j];
                int le = leftCellofFace[face];
                int re = rightCellofFace[face];

                LeftPressureOfCell += q[PHSPACE::IDX::IP][le];
                RightPressureOfCell += q[PHSPACE::IDX::IP][re];
            }
            rtem[iCell] += ABS(RightPressureOfCell - LeftPressureOfCell) / ABS(RightPressureOfCell + LeftPressureOfCell);
            rtem[iCell] /= static_cast <RDouble> (faceNumberOfEachCell[iCell]);
        }
    }

    if (RoeEntropyFixMethod == 6)
    {	
        for (int iCell = 0; iCell < nTotalCell; ++iCell)
        {
            int nNeighbor = c2c[iCell].size();

            RDouble minPressure = q[IP][iCell];
            RDouble maxPressure = q[IP][iCell];

            for (int j = 0; j < nNeighbor; ++j)
            {
                int neighborCell = c2c[iCell][j];
                if (q[IP][neighborCell] < minPressure)
                {
                    minPressure = q[IP][neighborCell];
                }
                if (q[IP][neighborCell] > maxPressure)
                {
                    maxPressure = q[IP][neighborCell];
                }
            }

            RDouble pressureRatio = maxPressure / minPressure;
            RDouble presFunc = min(one, max(zero, (maxPressureRatio - pressureRatio) / (maxPressureRatio - minPressureRatio)));
            RDouble pressCosFunc = cos(presFunc*pi);
            pressCosFunc = max(-1.0, min(pressCosFunc, one));

            rtem[iCell] = half * (one - pressCosFunc);
        }
    }

    grid->SetGhostCellExceptInterface(rtem);
}

void NSSolverUnstruct::GetHeatTransferCoeff(Grid *gridIn)
{
}

void NSSolverUnstruct::ViscousFlux(Grid *gridIn)
{
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int viscousType = parameters->GetViscousType();
    if (viscousType == INVISCID) return;

    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalFace = grid->GetNTotalFace();

    FaceProxy *faceProxy = CreateFaceProxy(grid);
    faceProxy->SetNSFaceValue(CreateNSFaceValue(grid));
    int localStart, localEnd;
    localStart = 0;
    do
    {
        localEnd = localStart + SEGCTION_LENGTH;
        if (localEnd > nTotalFace)
        {
            localEnd = nTotalFace;
        }

        ComputeFaceWeight(grid, faceProxy, localStart, localEnd);

        GetVisFaceValue(grid, faceProxy, localStart, localEnd);

        ComputeVisflux(grid, faceProxy, localStart, localEnd);

        LoadFlux(grid, faceProxy, localStart, localEnd);

        localStart = localEnd;
    } while (localStart < nTotalFace);

    delete faceProxy;
}

void NSSolverUnstruct::ComputeFaceWeight(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    RDouble *deltL = faceProxy->GetWeightL();
    RDouble *deltR = faceProxy->GetWeightR();

    grid->FaceWeight(deltL, deltR, localStart, localEnd);
}

void NSSolverUnstruct::GetVisFaceValue(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nBoundFace = grid->GetNBoundFace();
    int nTotalCell = grid->GetNTotalCell();

    NSFaceValue *faceVariable = faceProxy->GetNSFaceValue();

    RDouble *faceVelocityX = grid->GetFaceVelocityX();
    RDouble *faceVelocityY = grid->GetFaceVelocityY();
    RDouble *faceVelocityZ = grid->GetFaceVelocityZ();

    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    RDouble **q = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble **t = reinterpret_cast<RDouble **> (grid->GetDataPtr("t"));

    RDouble *viscousLaminar = reinterpret_cast<RDouble *> (grid->GetDataPtr("visl"));
    RDouble *viscousTurbulence = reinterpret_cast<RDouble *> (grid->GetDataPtr("vist"));

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    RDouble refGama = 1.4;
    GlobalDataBase::GetData("refGama", &refGama, PHDOUBLE, 1);

    RDouble oPrandtlLaminar    = parameters->GetoPrandtlLaminar();
    RDouble oPrandtlTurbulence = parameters->GetoPrandtlTurbulence();

    RDouble referenceMachNumber = parameters->GetRefMachNumber();

    RDouble tWall = parameters->GetWallTemperature();

    RDouble refDimensionalTemperature = parameters->GetRefDimensionalTemperature();

    int numberOfSpecies = parameters->GetNumberOfSpecies();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    RDouble **rhoDiffusion = faceVariable->GetRhoDS()->GetField();
    RDouble **hintSpecies = faceVariable->GetHintS()->GetField();
    RDouble **primitiveVariable = faceVariable->GetPrim()->GetField();
    RDouble **tMid   = faceVariable->GetT()->GetField();

    RDouble *kCp = faceVariable->GetKCP();   //! Heat conductivity coefficient, usually computed by Cp.
    RDouble *viscousLaminarFace = faceVariable->GetMUL();
    RDouble *viscousTurbulenceFace = faceVariable->GetMUT();

    RDouble referenceMachNumberSquare = referenceMachNumber * referenceMachNumber;

    using namespace IDX;

    RDouble *qL               = new RDouble[nEquation];
    RDouble *qR               = new RDouble[nEquation];
    RDouble *primface         = new RDouble[nEquation];
    RDouble *rhoDiffusionFace      = new RDouble[numberOfSpecies];
    RDouble *hintSpeciesOfFace = new RDouble[numberOfSpecies];

    int le,re,jFace;
    for (int iFace = localStart; iFace < localEnd; ++ iFace)
    {
        jFace = iFace - localStart;

        le = leftCellOfFace [iFace];
        re = rightCellOfFace[iFace];        

        viscousLaminarFace[jFace]   = half * (viscousLaminar[le] + viscousLaminar[re]);
        viscousTurbulenceFace[jFace] = half * (viscousTurbulence[le] + viscousTurbulence[re]);
    }

    int nMid = nBoundFace;
    if (localEnd < nBoundFace) nMid = localEnd;

    for (int iFace = localStart; iFace < localEnd; ++ iFace)
    {
        le = leftCellOfFace [iFace];
        re = rightCellOfFace[iFace];
        if (iFace < nBoundFace) re = iFace + nTotalCell;

        jFace = iFace - localStart;

        for (int m = 0; m < nEquation; ++ m)
        {
            primitiveVariable[m][jFace] = half * (q[m][le] + q[m][re]);
        }

        tMid[ITT][jFace] = half * (t[ITT][le] + t[ITT][re]);

        if (nChemical == 1)
        {
            using namespace GAS_SPACE;

            RDouble t1  = half * (t[ITT][le] + t[ITT][re]);

            for (int m = 0; m < nEquation; ++ m)
            {
                qL[m] = q[m][le];
                qR[m] = q[m][re];
                primface[m] = primitiveVariable[m][jFace];
            }

            gas->ComputeDensityDiffusionAndKCP(t1, primface, viscousLaminarFace[jFace], viscousTurbulenceFace[jFace], rhoDiffusionFace, hintSpeciesOfFace, kCp[jFace]);
            for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
            {
                rhoDiffusion[ispecies][jFace] = rhoDiffusionFace[ispecies];
            }

            for (int ispecies = 0; ispecies < numberOfSpecies; ++ ispecies)
            {
                hintSpecies[ispecies][jFace] = hintSpeciesOfFace[ispecies];
            }
        }
        else
        {
            RDouble cp = 1.0 / ((refGama - 1.0) * referenceMachNumberSquare);
            kCp[jFace] = (viscousLaminarFace[jFace] * oPrandtlLaminar + viscousTurbulenceFace[jFace] * oPrandtlTurbulence) * cp;
        }
    }

    int isUnsteady     = GlobalDataBase::GetIntParaFromDB("iunsteady");
    int isAle = GlobalDataBase::GetIntParaFromDB("codeOfAleModel");

    //! If there are some boundary faces, make modifications. 
    //! Check if there are boundary faces. If no, return.
    if (localStart < nBoundFace)
    {
        UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
        int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();
        
        RDouble tw = tWall / refDimensionalTemperature;
        for (int iFace = localStart; iFace < nMid; ++ iFace)
        {
            UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);
            int bcType = bcRegion->GetBCType();
            int jFace  = iFace - localStart;

            le = leftCellOfFace[iFace];
            re = rightCellOfFace[iFace];

            if (bcType == PHENGLEI::INTERFACE) continue;

            for (int m = 0; m < nEquation; ++ m)
            {
                RDouble tmp = half * (q[m][le] + q[m][re]);
                primitiveVariable[m][jFace] = tmp;
            }

            if (bcType == PHENGLEI::SOLID_SURFACE)
            {
                RDouble uWall = 0.0;
                RDouble vWall = 0.0;
                RDouble wWall = 0.0;
                Data_Param *bcData  = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace])->GetBCParamDataBase();
                if(bcData)
                {
                    if (bcData->IsExist("uWall", PHDOUBLE, 1))
                    {
                        bcData->GetData("uWall", &uWall, PHDOUBLE, 1);
                        bcData->GetData("vWall", &vWall, PHDOUBLE, 1);
                        bcData->GetData("wWall", &wWall, PHDOUBLE, 1);
                    }
                }
                
                RDouble velocityXWall = uWall;
                RDouble velocityYWall = vWall;
                RDouble velocityZWall = wWall;
                if (isUnsteady && isAle)
                {
                    velocityXWall = faceVelocityX[jFace] + uWall;
                    velocityYWall = faceVelocityY[jFace] + vWall;
                    velocityZWall = faceVelocityZ[jFace] + wWall;
                }

                primitiveVariable[IU][jFace] = velocityXWall;
                primitiveVariable[IV][jFace] = velocityYWall;
                primitiveVariable[IW][jFace] = velocityZWall;

                if (tw > 0.0)
                {
                    tMid[0][jFace] = tw;
                }
            }
        }
    }

    delete [] qL;
    delete [] qR;
    delete [] primface;
    delete [] rhoDiffusionFace;
    delete [] hintSpeciesOfFace;
}

void NSSolverUnstruct::ComputeVisflux(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nBoundFace = grid->GetNBoundFace();

    int *leftCellofFace = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();

    RDouble **dqdx = reinterpret_cast <RDouble **> (gridIn->GetDataPtr("gradPrimtiveVarX"));
    RDouble **dqdy = reinterpret_cast <RDouble **> (gridIn->GetDataPtr("gradPrimtiveVarY"));
    RDouble **dqdz = reinterpret_cast <RDouble **> (gridIn->GetDataPtr("gradPrimtiveVarZ"));

    RDouble **dtdx = reinterpret_cast <RDouble **> (gridIn->GetDataPtr("gradTemperatureX"));
    RDouble **dtdy = reinterpret_cast <RDouble **> (gridIn->GetDataPtr("gradTemperatureY"));
    RDouble **dtdz = reinterpret_cast <RDouble **> (gridIn->GetDataPtr("gradTemperatureZ"));

    RDouble  **flux = faceProxy->GetFlux();

    RDouble *deltaL = faceProxy->GetWeightL();
    RDouble *deltaR = faceProxy->GetWeightR();

    NSFaceValue *faceVariable = faceProxy->GetNSFaceValue();

    RDouble *xfn  = grid->GetFaceNormalX();
    RDouble *yfn  = grid->GetFaceNormalY();
    RDouble *zfn  = grid->GetFaceNormalZ();
    RDouble *area = grid->GetFaceArea();
    
    RDouble *xcc  = grid->GetCellCenterX();
    RDouble *ycc  = grid->GetCellCenterY();
    RDouble *zcc  = grid->GetCellCenterZ();

    RDouble *xfc  = grid->GetFaceCenterX();
    RDouble *yfc  = grid->GetFaceCenterY();
    RDouble *zfc  = grid->GetFaceCenterZ();

    RDouble **rhoDiffusion = faceVariable->GetRhoDS()->GetField();
    RDouble **hintSpecies = faceVariable->GetHintS()->GetField();
    RDouble **primitiveVariableFace = faceVariable->GetPrim()->GetField();
    RDouble **tm     = faceVariable->GetT()->GetField();

    RDouble *kCp     = faceVariable->GetKCP();
    RDouble *viscousLaminarFace   = faceVariable->GetMUL();
    RDouble *viscousTurbulenceFace = faceVariable->GetMUT();

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int numberOfSpecies = parameters->GetNumberOfSpecies();

    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    int nm = parameters->GetNSEquationNumber();

    RDouble refReNumber  = parameters->GetRefReNumber();
    RDouble oRefReNumber = parameters->GetoRefReNumber();

    int nolstress = GlobalDataBase::GetIntParaFromDB("nolstress");
    int nrokplus  = GlobalDataBase::GetIntParaFromDB("nrokplus");

    RDouble skewnessAngle = parameters->GetSkewnessAngle();

    RDouble **t    = reinterpret_cast<RDouble **> (grid->GetDataPtr("t"));
    RDouble **primitiveVariable = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));

    RDouble **qTurb = 0;
    RDouble ** aniss  = 0;

    if (nrokplus > 0)
    {
        qTurb = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
        aniss  = reinterpret_cast<RDouble **> (grid->GetDataPtr("aniss"));
    }

    RDouble *fvis = new RDouble[nLaminar];

    RDouble *dfdx = new RDouble[nEquation];
    RDouble *dfdy = new RDouble[nEquation];
    RDouble *dfdz = new RDouble[nEquation];

    RDouble *dfd1  = new RDouble[nEquation];
    RDouble *dfd2  = new RDouble[nEquation];
    RDouble *dfdn  = new RDouble[nEquation];
    RDouble *dfdt1 = new RDouble[nEquation];
    RDouble *dfdt2 = new RDouble[nEquation];

    RDouble *f1    = new RDouble[nEquation];
    RDouble *f2    = new RDouble[nEquation];
    RDouble *fMid  = new RDouble[nEquation];

    RDouble t1x, t1y, t1z, t2x, t2y, t2z;
    RDouble txx,tyy,tzz;
    RDouble txy,txz,tyz;
    
    bool isFineGrid = grid->IsFinestGrid();
    
    using namespace IDX;
    for (int iFace = localStart; iFace < localEnd; ++ iFace)
    {
        int le = leftCellofFace [iFace];
        int re = rightCellofFace[iFace];
        int jFace  = iFace - localStart;

        RDouble nxs = xfn[iFace];
        RDouble nys = yfn[iFace];
        RDouble nzs = zfn[iFace];

        //! Get first tangential vector on the face.
        if (abs(nxs) > SMALL)
        {
            t1x =   nys;
            t1y = - nxs;
            t1z =   0.0;
        }
        else if (abs(nys) > SMALL)
        {
            t1x = - nys;
            t1y =   nxs;
            t1z =   0.0;
        }
        else if (abs(nzs) > SMALL)
        {
            t1x =   0.0;
            t1y = - nzs;
            t1z =   nys;
        }
        else
        {
            for (int m = 0; m < nLaminar; ++ m)
            {
                flux[m][jFace] = 0.0;
            }
            continue;
        }

        //! Normalize the tangential vector.
        RDouble oNormal = 1.0 / DISTANCE(t1x, t1y, t1z);
        t1x *= oNormal;
        t1y *= oNormal;
        t1z *= oNormal;

        //! Get the second tangential vector by cross dot t1 to normal.
        t2x = nys * t1z - nzs * t1y;
        t2y = nzs * t1x - nxs * t1z;
        t2z = nxs * t1y - nys * t1x;

        RDouble dxL = xcc[le] - xfc[iFace];
        RDouble dyL = ycc[le] - yfc[iFace];
        RDouble dzL = zcc[le] - zfc[iFace];

        RDouble dxR = xcc[re] - xfc[iFace];
        RDouble dyR = ycc[re] - yfc[iFace];
        RDouble dzR = zcc[re] - zfc[iFace];

        RDouble dL = nxs * dxL + nys * dyL + nzs * dzL;
        RDouble dR = nxs * dxR + nys * dyR + nzs * dzR;

        RDouble distTemp = - dL / (DISTANCE(dxL, dyL, dzL) + SMALL);
        distTemp = MIN(distTemp,  1.0);
        distTemp = MAX(distTemp, -1.0);
        RDouble angle1 = asin(distTemp) * 180.0 / PI;
     
        distTemp = dR / (DISTANCE(dxR, dyR, dzR) + SMALL);
        distTemp = MIN(distTemp,  1.0);
        distTemp = MAX(distTemp, -1.0);
        RDouble angle2 = asin(distTemp) * 180.0 / PI;

        RDouble dxnL  = nxs * dL - dxL;
        RDouble dynL  = nys * dL - dyL;
        RDouble dznL  = nzs * dL - dzL;

        RDouble dxnR  = nxs * dR - dxR;
        RDouble dynR  = nys * dR - dyR;
        RDouble dznR  = nzs * dR - dzR;

        //! Quantities at points 1 and 2.
        for (int m = 0; m < nEquation; ++ m)
        {
            f1[m] = primitiveVariable[m][le];
            f2[m] = primitiveVariable[m][re];
        }

        RDouble tL = t[ITT][le];
        RDouble tR = t[ITT][re];

        for (int m = 0; m < nEquation; ++ m)
        {
            fMid[m] = half * (f1[m] + f2[m]);
        }
        RDouble tMid = half * (tL + tR);

        //! More accurate for skewness cells.
        //! If the cell degenerates too much, then interpolate the cell center value
        //! to the left and right points that vertical to the face normal.
        //! However, it will effect the robustness, so the skewnessAngle is suggested to be larger than 60 degree
        //! skewnessAngle, more large, more robust.
        //! For robust reason, this interpolation is NOT done on coarse grids.
        //! angle1/2 is defined as the complement of the angle between face normal and face-center/cell-center line.
        //! meaning that, (90 - angle1/2) is the angle between the face normal and the face-cell center line.
        //! angle1, angle2 = 0, represent completely orthogonal right angle, perfect quality.
        if (isFineGrid && angle1 > skewnessAngle && angle2 > skewnessAngle)
        {
            for (int m = 0; m < nEquation; ++ m)
            {
                f1[m] += dqdx[m][le] * dxnL + dqdy[m][le] * dynL + dqdz[m][le] * dznL;
                f2[m] += dqdx[m][re] * dxnR + dqdy[m][re] * dynR + dqdz[m][re] * dznR;
            }

            //! Only consider the rho limiter, may consider p limiter in future.
            tL += dtdx[ITT][le] * dxnL + dtdy[ITT][le] * dynL + dtdz[ITT][le] * dznL;
            tR += dtdx[ITT][re] * dxnR + dtdy[ITT][re] * dynR + dtdz[ITT][re] * dznR;

            if (tL < SMALL)
            {
                tL = t[ITT][le];
            }
            if (tR < SMALL)
            {
                tR = t[ITT][re];
            }

            //! Quantities at the face.
            for (int m = 0; m < nEquation; ++ m)
            {
                fMid[m] = primitiveVariableFace[m][jFace];
            }

            tMid = tm[ITT][jFace];
        }

        for (int m = 0; m < nEquation; ++ m)
        {
            dfdn[m] = 0.0;
        }
        RDouble dtdn = 0.0;

        if (angle1 > 0.0 && angle2 > 0.0 && ABS(dL) > TINY && ABS(dR) > TINY)
        {
            for (int m = 0; m < nEquation; ++ m)
            {
                dfd1[m] = (f1[m] - fMid[m]) / dL;
            }

            for (int m = 0; m < nEquation; ++ m)
            {
                dfd2[m] = (f2[m] - fMid[m]) / dR;
            }

            RDouble dtd1 = (tL - tMid) / dL;
            RDouble dtd2 = (tR - tMid) / dR;

            RDouble dtmp = dL * dL + dR * dR;
            RDouble weightL = dL * dL / dtmp;
            RDouble weightR = 1.0 - weightL;

            if (iFace < nBoundFace)
            {
                UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);
                int bcType = bcRegion->GetBCType();

                if (bcType != PHENGLEI::INTERFACE)
                {
                    weightL = 1.0;
                    weightR = 0.0;
                }
            }

            for (int m = 0; m < nEquation; ++ m)
            {
                dfdn[m] = weightL * dfd1[m] + weightR * dfd2[m];
            }

            dtdn = dtd1 * weightL + dtd2 * weightR;
        }

        RDouble weightL = deltaL[jFace];
        RDouble weightR = deltaR[jFace];

        for (int m = 0; m < nEquation; ++ m)
        {
            dfdx[m] = weightL * dqdx[m][le] + weightR * dqdx[m][re];
            dfdy[m] = weightL * dqdy[m][le] + weightR * dqdy[m][re];
            dfdz[m] = weightL * dqdz[m][le] + weightR * dqdz[m][re];
        }

        for (int m = IU; m <= IW; ++ m)
        {
            dfdt1[m] = t1x * dfdx[m] + t1y * dfdy[m] + t1z * dfdz[m];
            dfdt2[m] = t2x * dfdx[m] + t2y * dfdy[m] + t2z * dfdz[m];
        }

        //! Now true gradients. 
        for (int m = IU; m <= IW; ++ m)
        {
            dfdx[m] = nxs * dfdn[m] + t1x * dfdt1[m] + t2x * dfdt2[m];
            dfdy[m] = nys * dfdn[m] + t1y * dfdt1[m] + t2y * dfdt2[m];
            dfdz[m] = nzs * dfdn[m] + t1z * dfdt1[m] + t2z * dfdt2[m];
        }
        
        RDouble dudx = dfdx[IU];
        RDouble dudy = dfdy[IU];
        RDouble dudz = dfdz[IU];

        RDouble dvdx = dfdx[IV];
        RDouble dvdy = dfdy[IV];
        RDouble dvdz = dfdz[IV];

        RDouble dwdx = dfdx[IW];
        RDouble dwdy = dfdy[IW];
        RDouble dwdz = dfdz[IW];

        RDouble qNorm = 0.0;         
        if (nChemical == 1)
        {
            for (int iSpecies = 0; iSpecies < numberOfSpecies; ++ iSpecies)
            {
                qNorm += rhoDiffusion[iSpecies][jFace] * hintSpecies[iSpecies][jFace] * dfdn[iSpecies + nm];
            }

            for (int m = nm; m < nLaminar; ++ m)
            {
                fvis[m] = rhoDiffusion[m - nm][jFace] * dfdn[m];
            }
        }
        
        qNorm += kCp[jFace] * dtdn;
        RDouble divv2p3 = two3rd * (dudx + dvdy + dwdz);

        RDouble vis = viscousLaminarFace[jFace] + viscousTurbulenceFace[jFace];
        RDouble um  = fMid[1];
        RDouble vm  = fMid[2];
        RDouble wm  = fMid[3];

        //! Stress components.
        txx = vis * (two * dudx - divv2p3);
        tyy = vis * (two * dvdy - divv2p3);
        tzz = vis * (two * dwdz - divv2p3);
        txy = vis * (dudy + dvdx);
        txz = vis * (dudz + dwdx);
        tyz = vis * (dvdz + dwdy);

        if (nrokplus > 0)
        {
            RDouble rkl  = primitiveVariable[IR][le] * qTurb[0][le];
            RDouble rkr  = primitiveVariable[IR][re] * qTurb[0][re];
            RDouble rhok = half * (rkl + rkr) * refReNumber;
            if (nolstress > 0)
            {
                RDouble b11 = half * (aniss[0][le] + aniss[0][re]);
                RDouble b22 = half * (aniss[1][le] + aniss[1][re]);
                RDouble b12 = half * (aniss[2][le] + aniss[2][re]);
                RDouble b13 = half * (aniss[3][le] + aniss[3][re]);
                RDouble b23 = half * (aniss[4][le] + aniss[4][re]);
                RDouble b33 = - b11 - b22;
                txx = txx + b11 - two3rd * rhok;
                tyy = tyy + b22 - two3rd * rhok;
                tzz = tzz + b33 - two3rd * rhok;
                txy = txy + b12;
                txz = txz + b13;
                tyz = tyz + b23;
            }
            else
            {
                txx = txx - two3rd * rhok;
                tyy = tyy - two3rd * rhok;
                tzz = tzz - two3rd * rhok;
            }
        }

        fvis[IR ] = 0.0;
        fvis[IRU] = nxs * txx + nys * txy + nzs * txz;
        fvis[IRV] = nxs * txy + nys * tyy + nzs * tyz;
        fvis[IRW] = nxs * txz + nys * tyz + nzs * tzz;
        fvis[IRE] = um * fvis[IRU] + vm * fvis[IRV] + wm * fvis[IRW] + qNorm;

        for (int m = 0; m < nLaminar; ++ m)
        {
            flux[m][jFace] = - oRefReNumber * area[iFace] * fvis[m];
        }
    }

    delete []fvis;    fvis = NULL;

    delete []dfdx;    dfdx = NULL;
    delete []dfdy;    dfdy = NULL;
    delete []dfdz;    dfdz = NULL;

    delete []dfd1;    dfd1 = NULL;
    delete []dfd2;    dfd2 = NULL;
    delete []dfdn;    dfdn = NULL;
    delete []dfdt1;    dfdt1 = NULL;
    delete []dfdt2;    dfdt2 = NULL;

    delete []f1;    f1 = NULL;
    delete []f2;    f2 = NULL;
    delete []fMid;    fMid = NULL;
}

void NSSolverUnstruct::Boundary(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    RDouble wallTemperature = parameters->GetWallTemperature();
    int ifLowSpeedPrecon = parameters->GetIfLowSpeedPrecon();

    using namespace PHENGLEI;
    
    int isViscous = parameters->IsViscous();

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegion = unstructBCSet->GetnBCRegion();

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();
        if (bcType < 0)
        {
            continue;
        }
        else if (bcType == EXTRAPOLATION)
        {
            OutflowBCRegion(grid, bcRegion);
        }
        else if (IsWall(bcType))
        {
            if (!isViscous)
            {
                SymmetryBCRegion(grid, bcRegion);
            }
            else if (wallTemperature <= 0.0)
            {
                ViscousAdiabaticWallBCRegion(gridIn,bcRegion);
            }
            else
            {
                ViscousIsotropicWallBCRegion(grid, bcRegion);
            }
        }
        else if (bcType == SYMMETRY)
        {
            SymmetryBCRegion(grid, bcRegion);
        }
        else if (bcType == FARFIELD)
        {
            if (ifLowSpeedPrecon == 0)
            {
                FarfieldBCRegion(grid, bcRegion);
            }
            else
            {
                PreconFarfieldBCRegion(grid, bcRegion);
            }
        }
        else if (bcType == INFLOW)
        {
            InflowBCRegion(gridIn, bcRegion);
        }
        else if (bcType == OUTFLOW)
        {
            OutflowBCRegion(gridIn, bcRegion);
        }
        else if (bcType == PRESSURE_OUTLET)
        {
            PressureOutletBCRegion(gridIn, bcRegion);
        }
        else if (bcType == PRESSURE_INLET)
        {
            PressureInletBCRiemannRegion(gridIn, bcRegion);
        }
        else if (bcType == MASS_FLOW_INLET)
        {
            MassFlowInletBCRegion(gridIn, bcRegion);
        }
        else if (bcType == MASS_FLOW_OUTLET)
        {
            MassFlowOutletBCRegion(gridIn, bcRegion);
        }
        else if (bcType == OVERSET)
        {
            OutflowBCRegion(gridIn, bcRegion);
        }
        else
        {
            TK_Exit::ExceptionExit("Error: this boundary type does not exist!\n");
        }
    }
}

void NSSolverUnstruct::ComputeMassFlow(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    using namespace PHENGLEI;

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegion = unstructBCSet->GetnBCRegion();

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();
        
        if (bcType == INFLOW)
        {
            ComputeMassInFlow(gridIn, bcRegion);
        }
        else if (bcType == OUTFLOW)
        {
            ComputeMassOutFlow(gridIn, bcRegion);
        }
    }
}

void NSSolverUnstruct::ComputeMassInFlow(Grid *gridIn, UnstructBC *bcRegionUnstruct)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();
    RDouble *xfn = grid->GetFaceNormalX();
    RDouble *yfn = grid->GetFaceNormalY();
    RDouble *zfn = grid->GetFaceNormalZ();

    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));

    RDouble massIn  = 0.0;
    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();
    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        //! iFace is the face number in the set of faceIndex.
        int iFace = *iter;
        int le = leftCellofFace[iFace];
        int re = rightCellofFace[iFace];

        RDouble rmi, umi, vmi, wmi, vni;
        RDouble rmo, umo, vmo, wmo, vno;

        rmi = q[IR][le];
        umi = q[IU][le];
        vmi = q[IV][le];
        wmi = q[IW][le];
        vni = xfn[iFace] * umi + yfn[iFace] * vmi + zfn[iFace] * wmi;

        rmo = q[IR][re];
        umo = q[IU][re];
        vmo = q[IV][re];
        wmo = q[IW][re];
        vno = xfn[iFace] * umo + yfn[iFace] * vmo + zfn[iFace] * wmo;

        massIn -= half * (rmi * vni + rmo * vno);
    }
    GlobalDataBase::UpdateData("massIn", &massIn, PHDOUBLE, 1);
}

void NSSolverUnstruct::ComputeMassOutFlow(Grid *gridIn, UnstructBC *bcRegionUnstruct)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();
    RDouble *xfn = grid->GetFaceNormalX();
    RDouble *yfn = grid->GetFaceNormalY();
    RDouble *zfn = grid->GetFaceNormalZ();

    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));

    RDouble massOut = 0.0;
    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();
    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        //! iFace is the face number in the set of faceIndex.
        int iFace = *iter;
        int le = leftCellofFace[iFace];
        int re = rightCellofFace[iFace];

        RDouble rm, um, vm, wm, vn;
        rm = q[IR][le];
        um = q[IU][le];
        vm = q[IV][le];
        wm = q[IW][le];
        vn = xfn[iFace] * um + yfn[iFace] * vm + zfn[iFace] * wm;

        massOut += rm * vn;
    }
    GlobalDataBase::UpdateData("massOut", &massOut, PHDOUBLE, 1);
}

void NSSolverUnstruct::BoundaryQlQrFix(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nBoundFace = grid->GetNBoundFace();
    int nMid;

    //! Check if there are boundary faces. If no, return.
    if (localStart >= nBoundFace)
    {
        return;
    }

    if (localEnd <= nBoundFace)
    {
        //! If they are all boundary faces.
        nMid = localEnd;
    }
    else
    {
        //! Part of them are boundary faces.
        nMid = nBoundFace;
    }

    //! Get grid information.
    RDouble  *xfn  = grid->GetFaceNormalX();
    RDouble  *yfn  = grid->GetFaceNormalY();
    RDouble  *zfn  = grid->GetFaceNormalZ();
    RDouble *faceVelocityX = grid->GetFaceVelocityX();
    RDouble *faceVelocityY = grid->GetFaceVelocityY();
    RDouble *faceVelocityZ = grid->GetFaceVelocityZ();
    
    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    int isViscous = parameters->IsViscous();
    RDouble tWall = parameters->GetWallTemperature();
    RDouble refDimensionalTemperature = parameters->GetRefDimensionalTemperature();
    RDouble coefficientOfStateEquation = gas->GetCoefficientOfStateEquation();

    int isUnsteady     = GlobalDataBase::GetIntParaFromDB("iunsteady");
    int isAle = GlobalDataBase::GetIntParaFromDB("codeOfAleModel");

    using namespace PHSPACE;
    using namespace PHENGLEI;

    RDouble **qL = faceProxy->GetQL();
    RDouble **qR = faceProxy->GetQR();

    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();
    RDouble **q = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));

    RDouble *primitiveVariable = new RDouble[nEquation];

    using namespace IDX;

    for (int iFace = localStart; iFace < nMid; ++ iFace)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace]);
        int bcType = bcRegion->GetBCType();
        int jFace = iFace - localStart;

        if (bcType != PHENGLEI::SOLID_SURFACE && bcType != PHENGLEI::SYMMETRY)
        {
            continue;
        }

        int le = leftCellofFace[iFace];
        int re = rightCellofFace[iFace];

        for (int m = 0; m < nEquation; ++ m)
        {
            RFloat temp =  half * (q[m][le] + q[m][re]);
            qL[m][jFace] = temp;
            qR[m][jFace] = temp;
        }

        if (bcType == PHENGLEI::SYMMETRY)
        {
            RDouble vn = 2.0 * (xfn[iFace] * qL[IU][jFace] + yfn[iFace] * qL[IV][jFace] + zfn[iFace] * qL[IW][jFace]);

            qR[IR][jFace] = qL[IR][jFace];
            qR[IU][jFace] = qL[IU][jFace] - vn * xfn[iFace];
            qR[IV][jFace] = qL[IV][jFace] - vn * yfn[iFace];
            qR[IW][jFace] = qL[IW][jFace] - vn * zfn[iFace];
            qR[IP][jFace] = qL[IP][jFace];
        }

        if (bcType == PHENGLEI::SOLID_SURFACE && isViscous)
        {
            RDouble uWall = 0.0;
            RDouble vWall = 0.0;
            RDouble wWall = 0.0;                
            Data_Param *bcData  = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[iFace])->GetBCParamDataBase();
            if(bcData)
            {
                if (bcData->IsExist("uWall", PHDOUBLE, 1))
                {
                    bcData->GetData("uWall", &uWall, PHDOUBLE, 1);
                    bcData->GetData("vWall", &vWall, PHDOUBLE, 1);
                    bcData->GetData("wWall", &wWall, PHDOUBLE, 1);
                }
            }            

            RDouble velocityXWall = uWall;
            RDouble velocityYWall = vWall;
            RDouble velocityZWall = wWall;
            if (isUnsteady && isAle)
            {
                velocityXWall = faceVelocityX[iFace] + uWall;
                velocityYWall = faceVelocityY[iFace] + vWall;
                velocityZWall = faceVelocityZ[iFace] + wWall;
            }

            qL[IU][jFace] = velocityXWall;
            qL[IV][jFace] = velocityYWall;
            qL[IW][jFace] = velocityZWall;

            qR[IU][jFace] = velocityXWall;
            qR[IV][jFace] = velocityYWall;
            qR[IW][jFace] = velocityZWall;

            if (tWall <= 0.0)
            {
                //! Viscous WALL, adiabatic.
            }
            else
            {
                //! Iso-thermal wall.
                RDouble tw = tWall / refDimensionalTemperature;

                for (int m = 0; m < nEquation; ++ m)
                {
                    primitiveVariable[m] = qL[m][jFace];
                }

                RDouble omav = one;
                using namespace GAS_SPACE;
                if (nChemical == 1)
                {
                    omav = gas->ComputeMolecularWeightReciprocal(primitiveVariable);
                }

                RDouble pressureWall = primitiveVariable[IP];
                RDouble rhoWall  = pressureWall / (coefficientOfStateEquation * tw * omav);

                qL[IR][jFace] = rhoWall;
                qR[IR][jFace] = rhoWall;
            }
        }
    }

    delete []primitiveVariable;
}

void NSSolverUnstruct::GetInviscidFaceValue(Grid *gridIn, FaceProxy *faceProxy, Limiter * limiter, int localStart, int localEnd)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    GetQlQr(grid, faceProxy, localStart, localEnd);

    GetGamaLR(grid, faceProxy, localStart, localEnd);

	GetPressureFactorLR(grid, faceProxy, localStart, localEnd);

    ReConstructFaceValue(grid, faceProxy, limiter, localStart, localEnd);

    BoundaryQlQrFix(grid, faceProxy, localStart, localEnd);
}

//! Face values reconstruction by considering scalar limiters.
void NSSolverUnstruct::ReConstructFaceValue(Grid *gridIn, FaceProxy *faceProxy, Limiter * limiter, int localStart, int localEnd)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();

    RDouble *xfc  = grid->GetFaceCenterX();
    RDouble *yfc  = grid->GetFaceCenterY();
    RDouble *zfc  = grid->GetFaceCenterZ();

    RDouble *xcc  = grid->GetCellCenterX();
    RDouble *ycc  = grid->GetCellCenterY();
    RDouble *zcc  = grid->GetCellCenterZ();

    if(!grid->IsFinestGrid())
    {
        //! Reconstruction is not necessary for coarse grid.
        return;
    }

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int limiterType = parameters->GetLimiterType();
    int isInIniting = GlobalDataBase::GetIntParaFromDB("isInIniting");
    bool isFirstOrder = (limiterType == ILMT_FIRST) || isInIniting;
    if(isFirstOrder)
    {
        return;
    }

    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    bool isStructuredLimiter = (limiterType >= ILMT_STRUCT);
    if (isStructuredLimiter)
    { 
        //! Limiter for structured limiter.
        ReConstructQlQr_STR(gridIn, faceProxy, nEquation, localStart, localEnd);
        return;
    }

    RDouble *limit   = limiter->GetLimiter(0);
    RDouble **limit2D = new RDouble *[nEquation];

    int limitVector = parameters->GetLimitVector();
    bool usingVectorLimiter = (limitVector == 1);
    if (usingVectorLimiter)
    {
        for (int m = 0; m < nEquation; ++ m)
        {
            limit2D[m] = limiter->GetLimiter(m);
        }
    }

    RDouble **dqdx = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarX"));
    RDouble **dqdy = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarY"));
    RDouble **dqdz = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarZ"));

    RDouble **qL = faceProxy->GetQL();
    RDouble **qR = faceProxy->GetQR();

    RDouble *qTry = new RDouble[nEquation];

    int reconstructMethod = GlobalDataBase::GetIntParaFromDB("reconmeth");
    const int RECONSTRUCT_USING_RESPECTIVE_LIMITER = 0;
    if (reconstructMethod == RECONSTRUCT_USING_RESPECTIVE_LIMITER)
    {
        //! Q+, Q- reconstructing using their self limiters.
        for (int i = localStart; i < localEnd; ++ i)
        {
            int le,re,j;
            j  = i - localStart;
            le = leftCellofFace[i];
            re = rightCellofFace[i];

            RDouble dx = xfc[i] - xcc[le];
            RDouble dy = yfc[i] - ycc[le];
            RDouble dz = zfc[i] - zcc[le];

            for (int m = 0; m < nEquation; ++ m)
            {
                qTry[m] = qL[m][j];
            }

            for (int m = 0; m < nEquation; ++ m)
            {
                if(usingVectorLimiter)
                {
                    limit =  limit2D[m];
                }

                qTry[m] += limit[le] * (dqdx[m][le] * dx + dqdy[m][le] * dy + dqdz[m][le] * dz);
            }

            if (PositiveCheckForDensityAndPressure(qTry))
            {
                for (int m = 0; m < nEquation; ++ m)
                {
                    qL[m][j] = qTry[m];
                }
            }

            dx = xfc[i] - xcc[re];
            dy = yfc[i] - ycc[re];
            dz = zfc[i] - zcc[re];

            for (int m = 0; m < nEquation; ++ m)
            {
                qTry[m] = qR[m][j];
            }

            for (int m = 0; m < nEquation; ++ m)
            {
                if(usingVectorLimiter)
                {
                    limit =  limit2D[m];
                }

                qTry[m] += limit[re] * (dqdx[m][re] * dx + dqdy[m][re] * dy + dqdz[m][re] * dz);
            }

            if (PositiveCheckForDensityAndPressure(qTry))
            {
                for (int m = 0; m < nEquation; ++ m)
                {
                    qR[m][j] = qTry[m];
                }
            }
        }
    }
    else
    {
        //! Q+, Q- use the min limiter coefficients of left and right cell.
        for (int i = localStart; i < localEnd; ++ i)
        {
            int le,re,j;
            j  = i - localStart;
            le = leftCellofFace[i];
            re = rightCellofFace[i];

            RDouble dx = xfc[i] - xcc[le];
            RDouble dy = yfc[i] - ycc[le];
            RDouble dz = zfc[i] - zcc[le];

            for (int m = 0; m < nEquation; ++ m)
            {
                qTry[m] = qL[m][j];
            }

            for (int m = 0; m < nEquation; ++ m)
            {
                if(usingVectorLimiter)
                {
                    limit =  limit2D[m];
                }

                qTry[m] += MIN(limit[le], limit[re]) * (dqdx[m][le] * dx + dqdy[m][le] * dy + dqdz[m][le] * dz);
            }

            if (PositiveCheckForDensityAndPressure(qTry))
            {
                for (int m = 0; m < nEquation; ++ m)
                {
                    qL[m][j] = qTry[m];
                }
            }

            dx = xfc[i] - xcc[re];
            dy = yfc[i] - ycc[re];
            dz = zfc[i] - zcc[re];

            for (int m = 0; m < nEquation; ++ m)
            {
                qTry[m] = qR[m][j];
            }

            for (int m = 0; m < nEquation; ++ m)
            {
                if(usingVectorLimiter)
                {
                    limit =  limit2D[m];
                }

                qTry[m] += MIN(limit[le], limit[re]) * (dqdx[m][re] * dx + dqdy[m][re] * dy + dqdz[m][re] * dz);
            }

            if (PositiveCheckForDensityAndPressure(qTry))
            {
                for (int m = 0; m < nEquation; ++ m)
                {
                    qR[m][j] = qTry[m];
                }
            }
        }
    }
    delete []qTry;
    delete []limit2D;
}

void NSSolverUnstruct::ComputeInviscidFlux(Grid *grid, FaceProxy *faceProxy, int localStart, int localEnd)
{
    int unstructScheme = GlobalDataBase::GetIntParaFromDB("uns_scheme");
    int isInMGIniting = GlobalDataBase::GetIntParaFromDB("isInMGIniting");
    int isInIniting = GlobalDataBase::GetIntParaFromDB("isInIniting");
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nChemical = parameters->GetChemicalFlag();
    int precondition = parameters->GetIfLowSpeedPrecon();

    INVScheme inviscidScheme = GetINVScheme(unstructScheme);

    if ((!grid->IsFinestGrid() && isInMGIniting && !precondition) ||
        (grid->IsFinestGrid() && isInIniting   && !precondition))
    {
        //! The precondition has not been considered for vanleer method.
        //! so, only in the compressible multi-grid initialization state, run Vanleer.
        inviscidScheme = Vanleer_Scheme;
    }

    InviscidFluxWrap(grid, faceProxy, inviscidScheme, localStart, localEnd);
}

void NSSolverUnstruct::FillGeomProxy(Grid *gridIn, GeomProxy *geom_proxy, int localStart, int localEnd)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    RDouble *xfn  = grid->GetFaceNormalX();
    RDouble *yfn  = grid->GetFaceNormalY();
    RDouble *zfn  = grid->GetFaceNormalZ();
    RDouble *vgn  = grid->GetFaceNormalVelocity();
    RDouble *area = grid->GetFaceArea();

    std::copy(xfn  + localStart, xfn  + localEnd, geom_proxy->GetFaceNormalX());
    std::copy(yfn  + localStart, yfn  + localEnd, geom_proxy->GetFaceNormalY());
    std::copy(zfn  + localStart, zfn  + localEnd, geom_proxy->GetFaceNormalZ());
    std::copy(vgn  + localStart, vgn  + localEnd, geom_proxy->GetFaceVelocity());
    std::copy(area + localStart, area + localEnd, geom_proxy->GetFaceArea()   );
}

void NSSolverUnstruct::GetResidual(ActionKey *actkey)
{
    int level = actkey->level;
    UnstructGrid *grid = UnstructGridCast(GetGrid(level));

    RDouble **res = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res"));

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    PHSPACE::ComputeResidualonGrid(grid, actkey, res, nEquation);
}

void NSSolverUnstruct::InviscidFluxWrap(Grid *gridIn, FaceProxy *faceProxy, INVScheme inviscidScheme, int localStart, int localEnd)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar = parameters->GetLaminarNumber();
    int nFaceLength = localEnd - localStart;
    faceProxy->setsize(nFaceLength);

    RDouble **flux = faceProxy->GetFlux();
    GeomProxy *geomProxy = faceProxy->GetGeomProxy();

    if(systemgridtype == UNSTRUCTGRID)
    {
        FillGeomProxy(grid, geomProxy, localStart, localEnd);
    }
    else
    {
        //! mixgrid
        //! geomProxy has been set up once for ever.
    }
    

    //! Inviscid flux scheme parameters proxy.
    InviscidSchemeParameter inviscidSchemeParameterProxy;
    SetInviscidSchemeParameters(&inviscidSchemeParameterProxy, faceProxy, parameters);

    //! The flux is actually computed here!
    inviscidScheme(faceProxy, &inviscidSchemeParameterProxy);

    RDouble *area = geomProxy->GetFaceArea();

    //! iFace is the global ID, jFace is the local ID in this section.
    for (int iFace = localStart; iFace < localEnd; ++ iFace)
    {
        int jFace = iFace - localStart;

        RDouble &areas = area[jFace];

        for (int m = 0; m < nLaminar; ++ m)
        {
            flux[m][jFace] *= areas;
        }
    }
}

void NSSolverUnstruct::RegisterCFDSolverInterfaceField()
{
    NSSolver::RegisterCFDSolverInterfaceField();

    NSSolver::RegisterCFDSolverInterfaceFieldUnstruct();
}

void NSSolverUnstruct::RegisterCFDSolverInterpointField()
{
    NSSolver::RegisterCFDSolverInterpointField();
}

void NSSolverUnstruct::RegisterCFDSolverOversetField()
{
    NSSolver::RegisterCFDSolverOversetField();
}

void NSSolverUnstruct::SolutionFix(Grid *gridIn, RDouble *primitiveVariable, int iCell)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    RDouble **q = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    int **cell2face  = grid->GetCell2Face();
    int *faceNumberOfEachCell = grid->GetFaceNumberOfEachCell();

    int *leftCellofFace = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();

    using namespace GAS_SPACE;
    using namespace IDX;


    for (int m = 0; m < nEquation; ++ m)
    {
        primitiveVariable[m] = zero;
    }

    int iFace, le, re, neighborCell;
    RDouble volumeSum, ovolumeSum, vsur;
    volumeSum = 0.0;
    for (int icell2face = 0; icell2face < faceNumberOfEachCell[iCell]; ++ icell2face)
    {
        iFace = cell2face[iCell][icell2face];
        le    = leftCellofFace[iFace];
        re    = rightCellofFace[iFace];

        neighborCell = le;
        if (iCell == le) neighborCell = re;

        vsur = one;

        volumeSum += vsur;

        for (int m = 0; m < nEquation; ++ m)
        {
            RDouble f = q[m][neighborCell];
            if (m == IR || m == IP)
            {
                f = ABS(f);
            }
            primitiveVariable[m] += f * vsur;
        }
    }

    ovolumeSum = 1.0 / volumeSum;

    for (int m = 0; m < nEquation; ++ m)
    {
        primitiveVariable[m] *= ovolumeSum;
    }
}

void NSSolverUnstruct::UpdateFlowField(Grid *gridIn, FieldProxy *qProxy, FieldProxy *dqProxy)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();

    RDouble **q  = qProxy->GetField_UNS();
    RDouble **dq = dqProxy->GetField_UNS();

    RDouble **qnew = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble *gamma = reinterpret_cast<RDouble *> (grid->GetDataPtr("gama"));
    RDouble **t    = reinterpret_cast<RDouble **> (grid->GetDataPtr("t"));

    RDouble *xcc  = grid->GetCellCenterX();
    RDouble *ycc  = grid->GetCellCenterY();
    RDouble *zcc  = grid->GetCellCenterZ();

    Param_NSSolverUnstruct *parameters = GetControlParameters();

    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    //! Three modes for temperature, including the translational-rotational temperature, vibrational temperature and electron temperature.
    RDouble rm, pm, temperature[3];


    int *faceNumberOfEachCell = grid->GetFaceNumberOfEachCell();

    RDouble gama;
    using namespace GAS_SPACE;
    using namespace IDX;

    RDouble *primitiveVarFarfield = parameters->GetPrimitiveVarFarfield();
    RDouble roo = primitiveVarFarfield[IR];
    RDouble poo = primitiveVarFarfield[IP];
    RDouble densityLimit  = 1.0e-6 * roo;
    RDouble pressureLimit = 1.0e-6 * poo;

    RDouble *primitiveVariable = new RDouble[nEquation];
    RDouble *qtry = new RDouble[nEquation];

    int *negativeCellList = new int [nTotalCell];
    PHSPACE::SetField(negativeCellList, 0, nTotalCell);

    int nNegativeCell           = 0;
    int mostNegativeCell        = -1;
    RDouble mostNegativePressure = PHSPACE::LARGE;

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        for (int m = 0; m < nEquation; ++ m)
        {
            primitiveVariable[m] = q[m][iCell];
        }

        gama = gamma[iCell];
        //! Obtain the temperature.
        temperature[ITT] = t[ITT][iCell];
        if (nTemperatureModel == 1)
        {
            temperature[ITV] = temperature[ITT];
            temperature[ITE] = temperature[ITT];
        }
        else
        {
            temperature[ITV] = t[ITV][iCell];
            temperature[ITE] = temperature[ITV];
            if (nTemperatureModel == 3)
            {
                temperature[ITE] = t[ITE][iCell];
            }
        }
        gas->Primitive2Conservative(primitiveVariable, gama, temperature[ITV], temperature[ITE], qtry);

        for (int m = 0; m < nLaminar; ++ m)
        {
            RDouble dqTemp = dq[m][iCell];
            qtry[m] += dqTemp;
        }

        for (int m = 0; m < nTemperatureModel-1; ++ m)
        {
            qtry[nLaminar + nChemical + m] += dq[nLaminar + nChemical + m][iCell];
        }

        //! As input parameter, temperature is the initial temperature of the next iteration step.
        //! As returned parameter, temperature saves the temperatures of the next time-advancing step.
        gas->Conservative2Primitive(qtry, gama, primitiveVariable, temperature);

        rm = primitiveVariable[IR];
        pm = primitiveVariable[IP];

        if (rm < densityLimit || pm < pressureLimit)
        {
            if (pm < mostNegativePressure)
            {
                negativeCellList[iCell] = 1;
                mostNegativePressure = pm;
                mostNegativeCell     = iCell;
            }
            nNegativeCell += 1;
            SolutionFix(grid, primitiveVariable, iCell);
        }

        if (nChemical == 1)
        {
            gas->NormalizePrimitive(primitiveVariable);
        }

        for (int m = 0; m < nEquation; ++ m)
        {
            qnew[m][iCell] = primitiveVariable[m];
        }
    }

    if (nNegativeCell > 0)
    {
        UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
        int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();

        cout.setf(ios::scientific);
        cout.precision(4);
        cout << "      Warning: negative pressure appears in " << nNegativeCell << " cells";
        cout << ", iter = " << GlobalDataBase::GetIntParaFromDB("outnstep") ;
        cout << ", level = " << grid->GetLevel();
        cout << ", zone  = " << grid->GetZoneID();
        cout << ", cellID = " << mostNegativeCell;
        cout << ", the minimum pressure appears at: (" 
             << xcc[mostNegativeCell] << ", " << ycc[mostNegativeCell] << ", " << zcc[mostNegativeCell] << "). ";

        int nBoundFace = grid->GetNBoundFace();
        int **cell2face = grid->GetCell2Face();
        int isInteriorCell = 1;
        for (int jFace = 0; jFace < faceNumberOfEachCell[mostNegativeCell]; ++ jFace)
        {
            int faceID = cell2face[mostNegativeCell][jFace];
            if(faceID < nBoundFace)
            {
                if(isInteriorCell)
                {
                    cout << endl;
                }

                isInteriorCell = 0;
                UnstructBC *bcRegion = unstructBCSet->GetBCRegion(bcRegionIDofBCFace[faceID]);
                int bcType = bcRegion->GetBCType();
                cout << "            face bcType = " << bcType << endl;
            }
        }

        if(isInteriorCell)
        {
            cout << " is interior cell ..." << endl;
        }
    }

    DelPointer(negativeCellList);
    DelPointer(primitiveVariable);
    DelPointer(qtry);
}

void NSSolverUnstruct::SetGhostDQLUSGS(Grid *gridIn, RDouble **dq, RDouble **q)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nBoundFace = grid->GetNBoundFace();
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    int isViscous = parameters->IsViscous();

    using namespace PHSPACE;
    using namespace GAS_SPACE; 

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegion = unstructBCSet->GetnBCRegion();

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();

        if (bcType == PHENGLEI::SYMMETRY || (IsWall(bcType) && !isViscous))
        {
            SetSymmBCGhostDQLUSGS(grid, bcRegion, dq, q);
        }
        else if (IsWall(bcType) && isViscous)
        {
            SetWallBCGhostDQLUSGS(grid, bcRegion, dq, q);
        }
        else if (bcType == PHENGLEI::EXTRAPOLATION)
        {
            SetExtraBCGhostDQLUSGS(grid, bcRegion, dq, q);
        }else
        {
            vector<int> *faceIndex = bcRegion->GetFaceIndex();
            for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
            {
                // iFace is the face number in the set of faceIndex.
                int iFace = *iter;
                int le = leftCellofFace[iFace];
                int re = rightCellofFace[iFace];

                //! FIX_ALL set DQ = 0.
                //! Far_FIELD, PROFILE and FIX_V have not been implemented yet.
                for (int m = 0; m < nLaminar; ++ m)
                {
                    dq[m][re] = 0.0;
                }
                for (int m = 0; m < nTemperatureModel - 1; m++)
                {
                    dq[nLaminar + nChemical + m][re] = 0.0;
                }
            }
        }
    }
}

void NSSolverUnstruct::SetSymmBCGhostDQLUSGS(Grid *gridIn, UnstructBC *bcRegionUnstruct, RDouble **dq, RDouble **q)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();
    //! Get grid metrics.
    RDouble *xfn  = grid->GetFaceNormalX();
    RDouble *yfn  = grid->GetFaceNormalY();
    RDouble *zfn  = grid->GetFaceNormalZ();

    RDouble *gama = reinterpret_cast<RDouble *> (grid->GetDataPtr("gama"));
    RDouble **t   = reinterpret_cast<RDouble **> (grid->GetDataPtr("t"));
    Param_NSSolverUnstruct *parameters = GetControlParameters();    
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    RDouble *primitiveVariable = new RDouble[nEquation];
    RDouble *Q = new RDouble[nEquation];

    RDouble gam1,tv = 0.0,te = 0.0;
    RDouble rm, um, vm, wm, pm, ub, vb, wb;
    RDouble rvn, ru, rv, rw, p_new;
    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();
    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        //! iFace is the face number in the set of faceIndex.
        int iFace = *iter;
        int le = leftCellofFace[iFace];
        int re = rightCellofFace[iFace];

        for (int m = 0; m < nEquation; ++ m)
        {
            primitiveVariable[m] = q[m][le];
        }

        tv = t[ITT][le];
        te = tv;
        if (nTemperatureModel > 1)     //! Consider Two-Temperature model and Three-Temperature model.
        {
            tv = t[ITV][le];
            te = tv;
            if (nTemperatureModel == 3)
            {
                te = t[ITE][le];
            }
        }
        gas->Primitive2Conservative(primitiveVariable, gama[le], tv, te, Q);

        for (int m = 0; m < nLaminar; ++ m)
        {
            Q[m] += dq[m][le];
        }
        for (int m = 0; m < nTemperatureModel - 1; m++)
        {
            Q[nLaminar + nChemical + m] += dq[nLaminar + nChemical + m][le];
        }

        //! Non-conservative variables in re in level n.
        rm = q[IR][re];
        um = q[IU][re];
        vm = q[IV][re];
        wm = q[IW][re];
        pm = q[IP][re];

        for (int m = 0; m < nLaminar; ++ m)
        {
            dq[m][re] = dq[m][le];
        }
        for (int m = 0; m < nTemperatureModel - 1; m++)
        {
            dq[nLaminar + nChemical + m][re] = dq[nLaminar + nChemical + m][le];
        }

        //! Inviscid wall or symmetry boundary.
        dq[IR][re] = dq[IR][le];
        rvn = 2.0 * (xfn[iFace]*dq[IRU][le] + yfn[iFace]*dq[IRV][le] + zfn[iFace]*dq[IRW][le]);
        dq[IRU][re] = dq[IRU][le] - xfn[iFace] * rvn;
        dq[IRV][re] = dq[IRV][le] - yfn[iFace] * rvn;
        dq[IRW][re] = dq[IRW][le] - zfn[iFace] * rvn;
        rvn   = 2.0 * (xfn[iFace]*Q[IRU] + yfn[iFace]*Q[IRV] + zfn[iFace]*Q[IRW]);
        ru    = Q[IRU] - xfn[iFace] * rvn;
        rv    = Q[IRV] - yfn[iFace] * rvn;
        rw    = Q[IRW] - zfn[iFace] * rvn;
        p_new = Q[IRE] - 0.5 * (Q[IRU]*Q[IRU] + Q[IRV]*Q[IRV] + Q[IRW]*Q[IRW]) / Q[IR];

        gam1 = gama[le] - 1.0;
        dq[IRE][re] = p_new - pm/gam1 + half*((ru*ru + rv*rv + rw*rw)/Q[IR] - rm*(um*um + vm*vm + wm*wm));
    }
}

void NSSolverUnstruct::SetWallBCGhostDQLUSGS(Grid *gridIn, UnstructBC *bcRegionUnstruct, RDouble **dq, RDouble **q)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();
    //! Get grid metrics.
    RDouble *xfn  = grid->GetFaceNormalX();
    RDouble *yfn  = grid->GetFaceNormalY();
    RDouble *zfn  = grid->GetFaceNormalZ();

    RDouble *gama = reinterpret_cast<RDouble *> (grid->GetDataPtr("gama"));
    RDouble **t   = reinterpret_cast<RDouble **> (grid->GetDataPtr("t"));
    Param_NSSolverUnstruct *parameters = GetControlParameters();    
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    RDouble *primitiveVariable = new RDouble[nEquation];
    RDouble *Q = new RDouble[nEquation];

    RDouble gam1,tv = 0.0,te = 0.0;
    RDouble rm, um, vm, wm, pm, ub, vb, wb;
    RDouble rvn, ru, rv, rw, p_new;

    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();
    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        //! iFace is the face number in the set of faceIndex.
        int iFace = *iter;
        int le = leftCellofFace[iFace];
        int re = rightCellofFace[iFace];
    
        //! Conservative Variables on level n+1.
    
        for (int m = 0; m < nEquation; ++ m)
        {
            primitiveVariable[m] = q[m][le];
        }
    
        tv = t[ITT][le];
        te = tv;
        if (nTemperatureModel > 1)     //! Consider Two-Temperature model and Three-Temperature model.
        {
            tv = t[ITV][le];
            te = tv;
            if (nTemperatureModel == 3)
            {
                te = t[ITE][le];
            }
        }
        gas->Primitive2Conservative(primitiveVariable, gama[le], tv, te, Q);
    
        for (int m = 0; m < nLaminar; ++ m)
        {
            Q[m] += dq[m][le];
        }
        for (int m = 0; m < nTemperatureModel - 1; m++)
        {
            Q[nLaminar + nChemical + m] += dq[nLaminar + nChemical + m][le];
        }
    
        //! Non-conservative variables in re in level n.
        rm = q[IR][re];
        um = q[IU][re];
        vm = q[IV][re];
        wm = q[IW][re];
        pm = q[IP][re];

        for (int m = 0; m < nLaminar; ++ m)
        {
            dq[m][re] = dq[m][le];
        }
        for (int m = 0; m < nTemperatureModel - 1; m++)
        {
            dq[nLaminar + nChemical + m][re] = dq[nLaminar + nChemical + m][le];
        }
    
        ub = 0.0;
        vb = 0.0;
        wb = 0.0;
        dq[IR ][re] =   dq[IR ][le];
        dq[IRU][re] = - dq[IRU][le] + two * ub * dq[IR][le];
        dq[IRV][re] = - dq[IRV][le] + two * vb * dq[IR][le];
        dq[IRW][re] = - dq[IRW][le] + two * wb * dq[IR][le];
        dq[IRE][re] =   dq[IRE][le] + ub * (dq[IRU][re] - dq[IRU][le])
            + vb * (dq[IRV][re] - dq[IRV][le])
            + wb * (dq[IRW][re] - dq[IRW][le]);
    }
}

void NSSolverUnstruct::SetExtraBCGhostDQLUSGS(Grid *gridIn, UnstructBC *bcRegionUnstruct, RDouble **dq, RDouble **q)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();
    //! Get grid metrics.
    RDouble *xfn  = grid->GetFaceNormalX();
    RDouble *yfn  = grid->GetFaceNormalY();
    RDouble *zfn  = grid->GetFaceNormalZ();

    RDouble *gama = reinterpret_cast<RDouble *> (grid->GetDataPtr("gama"));
    RDouble **t   = reinterpret_cast<RDouble **> (grid->GetDataPtr("t"));
    Param_NSSolverUnstruct *parameters = GetControlParameters();    
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    RDouble *primitiveVariable = new RDouble[nEquation];
    RDouble *Q = new RDouble[nEquation];

    RDouble gam1,tv = 0.0,te = 0.0;
    RDouble rm, um, vm, wm, pm, ub, vb, wb;
    RDouble rvn, ru, rv, rw, p_new;

    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();
    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        // iFace is the face number in the set of faceIndex.
        int iFace = *iter;
        int le = leftCellofFace[iFace];
        int re = rightCellofFace[iFace];

        //! Conservative Variables on level n+1.

        for (int m = 0; m < nEquation; ++ m)
        {
            primitiveVariable[m] = q[m][le];
        }

        tv = t[ITT][le];
        te = tv;
        if (nTemperatureModel > 1)     //! Consider Two-Temperature model and Three-Temperature model.
        {
            tv = t[ITV][le];
            te = tv;
            if (nTemperatureModel == 3)
            {
                te = t[ITE][le];
            }
        }
        gas->Primitive2Conservative(primitiveVariable, gama[le], tv, te, Q);

        for (int m = 0; m < nLaminar; ++ m)
        {
            Q[m] += dq[m][le];
        }
        for (int m = 0; m < nTemperatureModel - 1; m++)
        {
            Q[nLaminar + nChemical + m] += dq[nLaminar + nChemical + m][le];
        }

        //! Non-conservative variables in re in level n.
        rm = q[IR][re];
        um = q[IU][re];
        vm = q[IV][re];
        wm = q[IW][re];
        pm = q[IP][re];

        for (int m = 0; m < nLaminar; ++ m)
        {
            dq[m][re] = dq[m][le];
        }
        for (int m = 0; m < nTemperatureModel - 1; m++)
        {
            dq[nLaminar + nChemical + m][re] = dq[nLaminar + nChemical + m][le];
        }

        for (int m = 0; m < nLaminar; ++ m)
        {
            dq[m][re] =  dq[m][le];
        }
        for (int m = 0; m < nTemperatureModel - 1; m++)
        {
            dq[nLaminar + nChemical + m][re] = dq[nLaminar + nChemical + m][le];
        }
    }
}

RDouble NSSolverUnstruct::GetRad(RDouble *prim, RDouble *nxyz, RDouble gama)
{
    using namespace IDX;

    RDouble &xfn  = nxyz[0];
    RDouble &yfn  = nxyz[1];
    RDouble &zfn  = nxyz[2];
    RDouble &area = nxyz[3];
    RDouble &vgn  = nxyz[4];

    RDouble &rm = prim[IR];
    RDouble &um = prim[IU];
    RDouble &vm = prim[IV];
    RDouble &wm = prim[IW];
    RDouble &pm = prim[IP];

    RDouble vn = xfn * um + yfn * vm + zfn * wm - vgn;
    RDouble eigen = ABS(vn) + sqrt(gama * pm / rm);
    return eigen * area;
}

void NSSolverUnstruct::SolveLUSGSForward(Grid *gridIn, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal, bool iAdvanceStep)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int *leftCellofFace = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();

    RDouble **dq = dqProxy->GetField_UNS();

    //! Get grid metrics.
    RDouble *xfn  = grid->GetFaceNormalX();
    RDouble *yfn  = grid->GetFaceNormalY();
    RDouble *zfn  = grid->GetFaceNormalZ();
    RDouble *area = grid->GetFaceArea();
    RDouble *vgn  = grid->GetFaceNormalVelocity();

    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    int *iBlank = grid->GetBlankIndex();

    //! Get flow variables.
    RDouble **q = reinterpret_cast <RDouble **> (grid->GetDataPtr("q"));
    RDouble **t = reinterpret_cast <RDouble **> (grid->GetDataPtr("t"));
    RDouble *gamma = reinterpret_cast <RDouble *> (grid->GetDataPtr("gama"));
    RDouble *diagonal = reinterpret_cast <RDouble *> (grid->GetDataPtr("diagonal"));
    RDouble **res = reinterpret_cast <RDouble **> (grid->GetDataPtr("res"));

    Param_NSSolverUnstruct *parameters = GetControlParameters();

    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    int nm = parameters->GetNSEquationNumber();

    int nElectronIndex = parameters->GetIndexOfElectron();

    RDouble refGama = parameters->GetRefGama();
    RDouble refReNumber  = parameters->GetRefReNumber();

    int ifLowSpeedPrecon = parameters->GetIfLowSpeedPrecon();

    bool isViscous = parameters->IsViscous();

    RDouble *viscousLaminar, *viscousTurbulence;
    if (isViscous)
    {
        viscousLaminar = reinterpret_cast<RDouble *> (grid->GetDataPtr("visl"));
        viscousTurbulence = reinterpret_cast<RDouble *> (grid->GetDataPtr("vist"));
    }

    //! Some intermediate variables.
    RDouble nxyz[5];
    RDouble *qNeighbor  = new RDouble[nEquation];
    RDouble *dqNeighbor = new RDouble[nEquation];
    RDouble *dqTry  = new RDouble[nEquation];
    RDouble *dqOld  = new RDouble [nEquation];
    RDouble *rhs0   = new RDouble[nEquation];
    RDouble *df     = new RDouble[nEquation];
    RDouble *facePrimitiveVariable = new RDouble [nEquation];
    RDouble *temperature = new RDouble[nTemperatureModel];
    RDouble *qPrecondition  = new RDouble [nEquation];
    RDouble **MG   = NewPointer2 <RDouble> (nEquation, nEquation);
    RDouble *Ux = new RDouble [nEquation];

    using namespace IDX;

    int **cell2face = grid->GetCell2Face();
    int *faceNumberOfEachCell = grid->GetFaceNumberOfEachCell();

    RDouble **deltaFlux = LUplusDQ->GetField_UNS();

    //! Now the Forward Sweep.
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        if (iBlank[iCell] <= 0)
        {
            for (int m = 0; m < nLaminar; ++ m)
            {
                dq[m][iCell] = 0.0;
            }
            continue;
        }

        for (int m = 0; m < nLaminar; ++ m)
        {
            rhs0[m] = 0.0;
            qPrecondition[m] = q[m][iCell];

            //! Back up the old dq, to compute the convergence.
            //! it is initialized by zero or res in the first sweep,
            //! then it is updated by backward sweep, that is dq*.
            dqOld[m] = dq[m][iCell];

            //! Here, the deltaFlux is computed in Upper backward sweep!
            //! res: b      deltaFlux: Ux, which is computed in backward sweep.
            //! the dq is not the real dq, now.
            //! Although the 'dq' changed can not the right one, it dosen't matter, since 
            //! the following only using the lower neighbor cells, whose 'dq' has been updated.
            //! It is convenient for precondition transform.
            dq[m][iCell] = res[m][iCell];

            Ux[m] = deltaFlux[m][iCell];

            //! Then reset it to zero to store the Lower forward sweep!
            deltaFlux[m][iCell] = 0.0;
        }

        if (ifLowSpeedPrecon != 0)
        {
            RDouble temperature = t[ITT][iCell];
            ComputePreconMatrix(1, MG, qPrecondition, temperature, refGama);
        }

        //! Compute the new Ux (deltaFlux).
        for (int jFace = 0; jFace < faceNumberOfEachCell[iCell]; ++ jFace)
        {
            int face = cell2face[iCell][jFace];
            int le   = leftCellofFace[face];
            int re   = rightCellofFace[face];

            //! One of le and re must be cell itself.
            if (le > iCell || re > iCell) continue;            

            nxyz[0] = xfn [face];
            nxyz[1] = yfn [face];
            nxyz[2] = zfn [face];
            nxyz[3] = area[face];
            nxyz[4] = vgn [face];

            if (re == iCell)
            {
                //! If the re is the iCell itself, then swap with 'le'.
                //! Now its neighboring cell 're' belongs to lower triangular.
                SWAP(le, re);

                //! Swap the normal and the moving face velocity.
                nxyz[0] = - nxyz[0];
                nxyz[1] = - nxyz[1];
                nxyz[2] = - nxyz[2];
                
                nxyz[4] = - nxyz[4];                
            }

            for (int m = 0; m < nEquation; ++ m)
            {
                qNeighbor[m] = q[m][re];
            }

            for (int m = 0; m < nLaminar; ++ m)
            {
                dqNeighbor[m] = dq[m][re];
            }

            //! Compute the ROE averaged face primitive variables.
            RDouble rhoRoe, uRoe, vRoe, wRoe, pRoe;
            RoeAveragedVariables(q[IDX::IR][le], q[IDX::IU][le], q[IDX::IV][le], q[IDX::IW][le], q[IDX::IP][le], gamma[le],
                                 q[IDX::IR][re], q[IDX::IU][re], q[IDX::IV][re], q[IDX::IW][re], q[IDX::IP][re], gamma[re],
                                 rhoRoe, uRoe, vRoe, wRoe, pRoe);

            facePrimitiveVariable[IDX::IR] = rhoRoe;
            facePrimitiveVariable[IDX::IU] = uRoe;
            facePrimitiveVariable[IDX::IV] = vRoe;
            facePrimitiveVariable[IDX::IW] = wRoe;
            facePrimitiveVariable[IDX::IP] = pRoe;

            //! The inviscid Delta Flux.
            if (ifLowSpeedPrecon == 0)
            {
                RDouble gama = gamma[re];

                RDouble rad = GetRad(facePrimitiveVariable, nxyz, gama);

                if (nChemical > 0)
                {
                    for (int m = 0; m < nTemperatureModel; ++ m)
                    {
                        temperature[m] = t[m][re];
                    }

                    ChemicalMXDQ(qNeighbor,temperature ,nxyz[0], nxyz[1], nxyz[2], nxyz[3], nxyz[4], dqNeighbor, df, nm, nLaminar, nTemperatureModel, rad, -1,nElectronIndex);
                }
                else
                {
                    MXDQ_STD(qNeighbor, nxyz[0], nxyz[1], nxyz[2], nxyz[3], nxyz[4], gama, dqNeighbor, df, nm, nLaminar, rad, -1);
                }
            }
            else
            {
                //! Compute the inviscid difference flux (df).
                PreconMatrixFreeDF(qNeighbor, dqNeighbor, nxyz, MG, refGama, df, nLaminar);
            }

            if (isViscous)
            {
                //! The viscous Delta Flux.
                RDouble distance;
                distance = ABS(xfn[face] * (xcc[re] - xcc[le])
                             + yfn[face] * (ycc[re] - ycc[le])
                             + zfn[face] * (zcc[re] - zcc[le]));

                RDouble viscousLaminarFace = half * (viscousLaminar[le] + viscousLaminar[re]);
                RDouble viscousTurbulenceFace = half * (viscousTurbulence[le] + viscousTurbulence[re]);
                RDouble density = half * (q[IDX::IR][le] + q[IDX::IR][re]);

                RDouble viscosity = viscousLaminarFace + viscousTurbulenceFace;
                RDouble faceArea  = area[face];
                RDouble viscousSpectrumRadius = 2.0 * viscosity / (density * distance * refReNumber + SMALL);
                viscousSpectrumRadius *= half * faceArea;

                for (int m = 0; m < nLaminar; ++ m)
                {
                    RDouble deltaVisFlux = viscousSpectrumRadius * dqNeighbor[m];
                    df[m] -= deltaVisFlux;
                }
            }

            //! Add Flux together.
            for (int m = 0; m < nLaminar; ++ m)
            {
                rhs0[m] += df[m];
            }
        }

        if(ifLowSpeedPrecon == 1)
        {
            //! Here, the 'dq' is 'res' actually.
            //! The reason why operate on 'dq' other than 'res', is that the 'dq' is a temporary variable.
            //! and the real 'res' is forbidden to be changed.
            PreconMatrixRes(dq, MG, nLaminar, iCell);
        }

        for (int m = 0; m < nLaminar; ++ m)
        {
            //! dq = { (b - Ux) - Lx } / D.
            //! Note: the 'dq' has actually initialized by the residual at the beginning of this function.
            //! the 'rhs0' is the total Delta Flux of the neighbors.            
            //! rhs0: Lx    diagonal: D.
            dqTry[m] = (dq[m][iCell] - Ux[m] - rhs0[m]) / diagonal[iCell];

            //! Store the lower forward sweep delta-flux, which will be used in the backward sweep.
            deltaFlux[m][iCell] += rhs0[m];
        }

        for (int m = 0; m < nLaminar; ++ m)
        {
            dq[m][iCell] = dqTry[m];

            sweepNormal += SQR(dq[m][iCell] - dqOld[m]);
        }
    } //! iCell loop.

    SetGhostDQLUSGS(grid, dq, q);
   
    delete []facePrimitiveVariable;
    delete []qNeighbor;
    delete []qPrecondition;
    delete []dqNeighbor;
    delete []dqTry;
    delete []dqOld;
    delete []rhs0;
    delete []df;
    DelPointer2(MG);
    delete []Ux;
}

void NSSolverUnstruct::SolveLUSGSBackward(Grid *gridIn, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal, bool iAdvanceStep)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();

    //! Here, the dq is the dq-star, which has been initialized in the forward sweep stage.
    //!       that it is solve by the lower-triangle equations.
    RDouble **dq = dqProxy->GetField_UNS();
    
    using namespace IDX;

    //! Get grid metrics.
    RDouble *xfn  = grid->GetFaceNormalX();
    RDouble *yfn  = grid->GetFaceNormalY();
    RDouble *zfn  = grid->GetFaceNormalZ();
    RDouble *area = grid->GetFaceArea();
    RDouble *vgn  = grid->GetFaceNormalVelocity();

    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    int *iBlank = grid->GetBlankIndex();

    //! Get flow variables.
    RDouble **q = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble **t = reinterpret_cast<RDouble **> (grid->GetDataPtr("t"));
    RDouble *gamma = reinterpret_cast<RDouble *> (grid->GetDataPtr("gama"));
    RDouble *diagonal = reinterpret_cast <RDouble *> (grid->GetDataPtr("diagonal"));
    RDouble **res = reinterpret_cast <RDouble **> (grid->GetDataPtr("res"));

    Param_NSSolverUnstruct *parameters = GetControlParameters();

    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    int nm    = parameters->GetNSEquationNumber();

    int nElectronIndex = parameters->GetIndexOfElectron();

    RDouble refGama = parameters->GetRefGama();
    RDouble refReNumber  = parameters->GetRefReNumber();

    int ifLowSpeedPrecon = parameters->GetIfLowSpeedPrecon();

    int isViscous = parameters->IsViscous();
    RDouble *viscousLaminar, *viscousTurbulence;
    if (isViscous)
    {
        viscousLaminar = reinterpret_cast<RDouble *> (grid->GetDataPtr("visl"));
        viscousTurbulence = reinterpret_cast<RDouble *> (grid->GetDataPtr("vist"));
    }

    //! Some intermediate variables.
    RDouble nxyz[5];
    RDouble *qNeighbor  = new RDouble [nEquation];
    RDouble *dqNeighbor = new RDouble [nEquation];
    RDouble *dqTry  = new RDouble [nEquation];
    RDouble *dqOld  = new RDouble [nEquation];
    RDouble *rhs0   = new RDouble [nEquation];
    RDouble *df     = new RDouble [nEquation];
    RDouble *facePrimitiveVariable = new RDouble [nEquation];
    RDouble *temperature = new RDouble[nTemperatureModel];
    RDouble *Lx = new RDouble [nEquation];
    RDouble *qPrecondition  = new RDouble [nEquation];

    RDouble **MG = NewPointer2 <RDouble> (nEquation, nEquation);
    
    int **cell2face = grid->GetCell2Face();
    int *faceNumberOfEachCell = grid->GetFaceNumberOfEachCell();

    RDouble **deltaFlux = LUplusDQ->GetField_UNS();

    //! Backward sweep.
    for (int iCell = nTotalCell - 1; iCell >= 0; -- iCell)
    {
        if (iBlank[iCell] <= 0)
        {
            for (int m = 0; m < nLaminar; ++ m)
            {
                dq[m][iCell] = 0.0;
            }
            continue;
        }

        for (int m = 0; m < nLaminar; ++ m)
        {
            rhs0[m] = 0.0;
            qPrecondition[m] = q[m][iCell];

            //! Back up the old dq, to compute the convergence.
            //! the 'dq' is dq*, which has been updated in forward.
            dqOld[m] = dq[m][iCell];

            //! the dq is not the real dq, now.
            //! it is convenient for precondition transform.
            dq[m][iCell] = res[m][iCell];

            Lx[m] = deltaFlux[m][iCell];

            deltaFlux[m][iCell] = 0.0;
        }

        if (ifLowSpeedPrecon != 0)
        {   
            RDouble temperature = t[ITT][iCell];
            ComputePreconMatrix(1, MG, qPrecondition, temperature, refGama);
        }

        for (int jFace = 0; jFace < faceNumberOfEachCell[iCell]; ++ jFace)
        {
            int face  = cell2face[iCell][jFace];
            int le    = leftCellofFace[face];
            int re    = rightCellofFace[face];

            //! One of le and re must be cell itself.
            if (le < iCell || re < iCell) continue;

            nxyz[0] = xfn [face];
            nxyz[1] = yfn [face];
            nxyz[2] = zfn [face];
            nxyz[3] = area[face];
            nxyz[4] = vgn [face];

            if (re == iCell)
            {
                //! If the re is the iCell itself, then swap with 'le'.
                //! Now its neighboring cell 're' belongs to upper triangular.
                SWAP(le, re);

                //! Swap the normal and the moving face velocity.
                nxyz[0] = - nxyz[0];
                nxyz[1] = - nxyz[1];
                nxyz[2] = - nxyz[2];

                nxyz[4] = - nxyz[4];
            }

            for (int m = 0; m < nEquation; ++ m)
            {
                qNeighbor [m] = q[m][re];
            }

            for (int m = 0; m < nLaminar; ++ m)
            {
                dqNeighbor[m] = dq[m][re];
            }

            //! Compute the ROE averaged face primitive variables.
            RDouble rhoRoe, uRoe, vRoe, wRoe, pRoe;
            RoeAveragedVariables(q[IDX::IR][le], q[IDX::IU][le], q[IDX::IV][le], q[IDX::IW][le], q[IDX::IP][le], gamma[le],
                q[IDX::IR][re], q[IDX::IU][re], q[IDX::IV][re], q[IDX::IW][re], q[IDX::IP][re], gamma[re],
                rhoRoe, uRoe, vRoe, wRoe, pRoe);

            facePrimitiveVariable[IDX::IR] = rhoRoe;
            facePrimitiveVariable[IDX::IU] = uRoe;
            facePrimitiveVariable[IDX::IV] = vRoe;
            facePrimitiveVariable[IDX::IW] = wRoe;
            facePrimitiveVariable[IDX::IP] = pRoe;

            if (ifLowSpeedPrecon  == 0 )
            {
                RDouble gama = gamma[re];

                RDouble rad = GetRad(facePrimitiveVariable,nxyz,gama);

                if (nChemical > 0)
                {
                    for (int m = 0; m < nTemperatureModel; ++ m)
                    {
                        temperature[m] = t[m][re];
                    }

                    ChemicalMXDQ(qNeighbor,temperature ,nxyz[0], nxyz[1], nxyz[2], nxyz[3], nxyz[4], dqNeighbor, df, nm, nLaminar, nTemperatureModel, rad, -1,nElectronIndex);
                }
                else
                {
                    MXDQ_STD(qNeighbor, nxyz[0], nxyz[1], nxyz[2], nxyz[3], nxyz[4], gama, dqNeighbor, df, nm, nLaminar, rad, -1);

                }
            }
            else
            {
                PreconMatrixFreeDF(qNeighbor, dqNeighbor, nxyz, MG, refGama, df, nLaminar);
            }

            if (isViscous)
            {
                RDouble dist;
                dist = ABS( xfn[face] * (xcc[re] - xcc[le])
                           + yfn[face] * (ycc[re] - ycc[le])
                           + zfn[face] * (zcc[re] - zcc[le]));

                RDouble vis_l = half * (viscousLaminar[le] + viscousLaminar[re]);
                RDouble vis_t = half * (viscousTurbulence[le] + viscousTurbulence[re]);
                RDouble density = half * (q[IDX::IR][le] + q[IDX::IR][re]);

                RDouble viscosity = vis_l + vis_t;
                RDouble faceArea  = area[face];
                RDouble viscousSpectrumRadius = 2.0 * viscosity / (density * dist * refReNumber + SMALL);
                viscousSpectrumRadius *= half * faceArea;

                for (int m = 0; m < nLaminar; ++ m)
                {
                    RDouble deltaVisFlux = viscousSpectrumRadius * dqNeighbor[m];
                    df[m] -= deltaVisFlux;
                }
                }

            //! Add Flux together.
            for (int m = 0; m < nLaminar; ++ m)
            {
                rhs0[m] += df[m];
            }
        }

        if(ifLowSpeedPrecon == 1)
        {
            //! Here, the 'dq' is 'res' actually.
            //! The reason why operate on 'dq' other than 'res', is that the 'dq' is a temporary variable.
            //! and the real 'res' is forbidden to be changed.
            PreconMatrixRes(dq, MG, nLaminar, iCell);
        }

        for (int m = 0; m < nLaminar; ++ m)
        {
            //! Note: the 'dq' has been updated by the forward sweep.
            //! the 'rhs0' is the total Delta Flux of the neighbors.
            //! x = {(b - LX) - Ux} / D.
            //! rhs0: Ux.    diagonal: D.
            dqTry[m] = (dq[m][iCell] - Lx[m] - rhs0[m]) / diagonal[iCell];

            //! Store the upper backward sweep delta-flux, which will be used in the forward sweep.
            deltaFlux[m][iCell] += rhs0[m];
        }

        //RestrictDQ(dqTry, nLaminar, dqlim);

        for (int m = 0; m < nLaminar; ++ m)
        {
            dq[m][iCell] = dqTry[m];

            sweepNormal += SQR(dq[m][iCell] - dqOld[m]);
        }
    }//! iCell loop.
   
    delete []facePrimitiveVariable;
    delete []qNeighbor;
    delete []dqNeighbor;
    delete []dqTry;
    delete []dqOld;
    delete []rhs0;
    delete []df;
    DelPointer2(MG);
    delete []Lx;
    delete []qPrecondition;
}


void NSSolverUnstruct::ZeroResidualOfSpecialCells(Grid *gridIn)
{
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int isOverset = parameters->GetIsOverLapping();
    if (!isOverset)
    {
        return;
    }

    UnstructGrid *grid = PHSPACE::UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();

    RDouble **residual = reinterpret_cast<RDouble **> (grid->GetDataPtr("res"));
    int nEquation = GetNumberOfEquations();
    int *iBlank = grid->GetBlankIndex();

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        if (iBlank[iCell] <= 0)
        {
            for (int m = 0; m < nEquation; ++ m)
            {
                residual[m][iCell] = 0.0;
            }
            continue;
        }
    }
    return;
}

void NSSolverUnstruct::SpectrumRadius(Grid *gridIn)
{
    int ifLowSpeedPrecon = GetControlParameters()->GetIfLowSpeedPrecon();
    if (ifLowSpeedPrecon == 0)
    {
        SpectrumRadiusInviscid(gridIn);
    }
    else
    {
        PreconInvSpectrumRadius(gridIn);
    }

    //! If the flow is viscous, it need to count the contribution from viscosity.
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int isViscous = parameters->GetViscousType();
    if (isViscous)
    {
        SpectrumRadiusViscous(gridIn);
    }

    int nChemical = parameters->GetChemicalFlag();
    int nChemicalRadius = parameters->GetNChemicalRadius();

    if (nChemical == 1 && nChemicalRadius == 1)
    {
        SpectrumRadiusChemical(UnstructGridCast(gridIn));
    }
}

void NSSolverUnstruct::Diagonal(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    RDouble *vol  = grid->GetCellVolume();

    RDouble *dt = reinterpret_cast<RDouble *> (grid->GetDataPtr("dt"));
    RDouble *invSpectralRadius = reinterpret_cast<RDouble *> (grid->GetDataPtr("invSpectralRadius"));
    RDouble *visSpectralRadius = reinterpret_cast<RDouble *> (grid->GetDataPtr("visSpectralRadius"));
    RDouble *diagonal = reinterpret_cast <RDouble *> (grid->GetDataPtr("diagonal"));

    using namespace IDX;

    //! If flow is unsteady, need to count the contribution of unsteady.
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int isUnsteady = parameters->GetIsUnsteady();
    if (!isUnsteady)
    {
        //! Part 1: Diagonal = volume / dt.
        //! notice: the dt[iCell] has divided vol, so the diagonal[iCell] has multiplied the vol by default.
        for (int iCell = 0; iCell < nTotal; ++ iCell)
        {
            diagonal[iCell] = 1.0 / dt[iCell];    //! notice: the dt[iCell] has divided vol, so the diagonal[iCell] has multiplied the vol by default.
        }
    }
    else
    {
        RDouble dualTimeCoefficient[7];
        const int methodOfDualTime = GlobalDataBase::GetIntParaFromDB("methodOfDualTime");
        ComputeDualTimeCoefficient(methodOfDualTime, dualTimeCoefficient);
        RDouble dualTimeSpectrumC1 = - dualTimeCoefficient[3];
        RDouble dualTimeSpectrumC2 =   dualTimeCoefficient[6];

        for (int iCell = 0; iCell < nTotal; ++ iCell)
        {
            diagonal[iCell] = dualTimeSpectrumC1 * vol[iCell] + (dualTimeSpectrumC2 / dt[iCell]);    //! notice: the dt[iCell] has divided vol, as same as above.
        }
    }

    //! Part 2: Diagonal += (|V.n| + c). (|V.n| + c) is inviscid spectral radius.
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        diagonal[iCell] += invSpectralRadius[iCell];
    }

    //! If the flow is viscous, it need to count the contribution from viscosity.
    int isViscous = parameters->IsViscous();
    if (isViscous)
    {
        //! Part 3: Diagonal += 2*miu/(rho.|n.dr|). 2*miu/(rho.|n.dr|) is viscous spectral radius.
        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            diagonal[iCell] += visSpectralRadius[iCell];
        }
    }

    int nm        = parameters->GetNSEquationNumber();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nChemicalRadius = parameters->GetNChemicalRadius();
    int nChemicalSource  = parameters->GetNChemicalSource();

    RDouble rad_chem;

    if (nChemical == 1 && nChemicalSource == 1 && nChemicalRadius == 1)
    {
        //! Part 4: Diagonal += chemical radius.
        RDouble **chemSpectralRadius = reinterpret_cast<RDouble **> (grid->GetDataPtr("chemSpectralRadius"));
        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            rad_chem = zero;
            for (int m = nm; m < nLaminar; ++ m)
            {
                rad_chem = MAX(nChemicalRadius * chemSpectralRadius[ m - nm ][ iCell ], rad_chem);
            }
            diagonal[iCell] += rad_chem;    //! Referring to the structure part, the dt[iCell] contained in the diagonal[iCell] can be improved for more stability here.
        }
    }
}

void NSSolverUnstruct::RestrictDQ(RDouble *dqtry, int nm, RDouble lim_max)
{
    for (int m = 0; m < nm; ++ m)
    {
        dqtry[m] = MAX(MIN(dqtry[m], lim_max),-lim_max);
    }
}

void NSSolverUnstruct::ModifyDensity(RDouble *dq, int nNSEquation,int nLaminar,int nEquation,RDouble &densityRelax)
{
    int nChemcalSourceModified = gas->GetnChemcalSourceModified();
    if(nChemcalSourceModified==0)
    {
        for (int m = 1; m < nNSEquation; ++ m)
        {
            dq[m] = 1.0;
        }
        for (int m = nNSEquation; m <= nLaminar; ++ m)
        {
            dq[m] = densityRelax;
        }
        for (int m = nLaminar+1; m < nEquation; ++ m)
        {
            dq[m] = 1.0;
        }
    }
    else
    {
        for (int m = 1; m < nNSEquation - 1; ++ m)
        {
            dq[m] = 1.0;
        }
        for (int m = nNSEquation - 1; m < nEquation; ++ m)
        {
            dq[m] = densityRelax;
        }
    }
}

void NSSolverUnstruct::AirForceCoef(ActionKey *actkey)
{
    const int ALL_PART = -1;
    AirForceCoefParts(actkey, ALL_PART);

    int Coordinate = LocalCoordinate;
    uint_t nTBC = GlobalBoundaryCondition::GetNumberOfBC();
    for (int iPart = 0; iPart < nTBC; ++ iPart)
    {
        if(!GlobalBoundaryCondition::IsSolidWallBC(iPart)) continue;

        //! Compute local part force.
        AirForceCoefParts(actkey, iPart, Coordinate);

        //! Compute global part force.
        AirForceCoefParts(actkey, iPart);
    }
}

void NSSolverUnstruct::AirForceCoefParts(ActionKey *actkey, int partID, int Coordinate)
{
    int level = actkey->level;

    using namespace PHMPI;

    UnstructGrid * grid = UnstructGridCast(GetGrid(level));

    RDouble **q   = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    RDouble *visl = reinterpret_cast< RDouble * > (grid->GetDataPtr("visl"));
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    RDouble oRefReNumber = parameters->GetoRefReNumber();

    int viscousType = parameters->GetViscousType();

    const int ALL_PART = -1;
    SimpleBC * boundaryCondition = NULL;
	string partBCName = "";
    RDouble TorqueRefX = 0.0, TorqueRefY = 0.0, TorqueRefZ = 0.0;

    //! Used to compute hinge moment:
    int dumpHingeMoment = 0;
    RDouble localCoordAxis0[3];
    RDouble localCoordAxis1[3];
    RDouble hingeMoment = 0.0;

    if(partID == ALL_PART)
    {
        TorqueRefX = GlobalDataBase::GetDoubleParaFromDB("TorqueRefX");
        TorqueRefY = GlobalDataBase::GetDoubleParaFromDB("TorqueRefY");
        TorqueRefZ = GlobalDataBase::GetDoubleParaFromDB("TorqueRefZ");
    }
    else if ((partID != ALL_PART) && (Coordinate == GlobalCoordinate))
    {
        boundaryCondition = GlobalBoundaryCondition::GetBC(partID);
        partBCName = boundaryCondition->GetBCName();

        TorqueRefX = GlobalDataBase::GetDoubleParaFromDB("TorqueRefX");
        TorqueRefY = GlobalDataBase::GetDoubleParaFromDB("TorqueRefY");
        TorqueRefZ = GlobalDataBase::GetDoubleParaFromDB("TorqueRefZ");
    }
    else if ((partID != ALL_PART) && (Coordinate == LocalCoordinate))
    {
        boundaryCondition = GlobalBoundaryCondition::GetBC(partID);
		partBCName = boundaryCondition->GetBCName();

        boundaryCondition->GetParamData("TorqueRefX", &TorqueRefX, PHDOUBLE, 1);
        boundaryCondition->GetParamData("TorqueRefY", &TorqueRefY, PHDOUBLE, 1);
        boundaryCondition->GetParamData("TorqueRefZ", &TorqueRefZ, PHDOUBLE, 1);
        
        
        boundaryCondition->GetParamData("dumpHingeMoment", &dumpHingeMoment, PHINT, 1);
        boundaryCondition->GetParamData("localCoordAxis0", localCoordAxis0, PHDOUBLE, 3);
        boundaryCondition->GetParamData("localCoordAxis1", localCoordAxis1, PHDOUBLE, 3);
    }

    int vis_run = 0;
    if (viscousType != INVISCID)
    {
        vis_run = 1;
    }

    using namespace IDX;
    RDouble *primitiveVarFarfield = parameters->GetPrimitiveVarFarfield();
    RDouble poo = primitiveVarFarfield[IP];

    RDouble *xfn  = grid->GetFaceNormalX();
    RDouble *yfn  = grid->GetFaceNormalY();
    RDouble *zfn  = grid->GetFaceNormalZ();
    RDouble *area = grid->GetFaceArea() ;

    RDouble *xfc = grid->GetFaceCenterX();
    RDouble *yfc = grid->GetFaceCenterY();
    RDouble *zfc = grid->GetFaceCenterZ();
    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    int * leftCellofFace = grid->GetLeftCellOfFace();
    int * rightCellofFace = grid->GetRightCellOfFace();

    RDouble cpx,cpy,cpz;
    cpx = 0.0;
    cpy = 0.0;
    cpz = 0.0;

    RDouble CA_f = 0.0;
    RDouble CA_p = 0.0;
    RDouble CN_f = 0.0;
    RDouble CN_p = 0.0;
    RDouble CZ_f = 0.0;
    RDouble CZ_p = 0.0;
    RDouble Cl_f = 0.0;
    RDouble Cl_p = 0.0;
    RDouble Cn_f = 0.0;
    RDouble Cn_p = 0.0;
    RDouble Cm_f = 0.0;
    RDouble Cm_p = 0.0;

    RDouble pw,cp;
    RDouble dfx,dfy,dfz,dpx,dpy,dpz,fvsx,fvsy,fvsz;
    RDouble nx,ny,nz;
    RDouble vis,txx,txy,txz,tyx,tyy,tyz,tzx,tzy,tzz;
    RDouble dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz;
    RDouble xc,yc,zc,dx,dy,dz,ods;

    RDouble **gradPrimtiveVarX = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarX"));
    RDouble **gradPrimtiveVarY = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarY"));
    RDouble **gradPrimtiveVarZ = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarZ"));

    using namespace IDX;

    //! Get the bcRegion information of unstruct grid.
    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegion = unstructBCSet->GetnBCRegion();

    //! Cycle all bcregions and find the solid surface type of boundary.
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();
        string bcName = bcRegion->GetBCName();

        //! If the bc type is wall, calculate the wall coefficient.
        if(bcType == PHENGLEI::SOLID_SURFACE)
        {
            if (partID != ALL_PART && bcName != partBCName)
            {
                continue;
            }
            vector<int> *faceIndex = bcRegion->GetFaceIndex();
            int iFace, le, re;
            for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
            {
                iFace = *iter;
                le = leftCellofFace [ iFace ];
                re = rightCellofFace[ iFace ];

                nx  = xfn[iFace] * area[iFace];
                ny  = yfn[iFace] * area[iFace];
                nz  = zfn[iFace] * area[iFace];

                dx  = xfc[iFace] - xcc[le];
                dy  = yfc[iFace] - ycc[le];
                dz  = zfc[iFace] - zcc[le];

                //! Pressure drag
                pw = q[IP][le];
                cp    = two * (pw - poo);

                dpx = nx * cp;
                dpy = ny * cp;
                dpz = nz * cp;

                CA_p += dpx;
                CN_p += dpy;
                CZ_p += dpz;

                dfx = dpx;
                dfy = dpy;
                dfz = dpz;

                xc = xfc[iFace];
                yc = yfc[iFace];
                zc = zfc[iFace];

                if (viscousType > INVISCID)
                {
                    dudx = gradPrimtiveVarX[IU][le];
                    dudy = gradPrimtiveVarY[IU][le];
                    dudz = gradPrimtiveVarZ[IU][le];

                    dvdx = gradPrimtiveVarX[IV][le];
                    dvdy = gradPrimtiveVarY[IV][le];
                    dvdz = gradPrimtiveVarZ[IV][le];

                    dwdx = gradPrimtiveVarX[IW][le];
                    dwdy = gradPrimtiveVarY[IW][le];
                    dwdz = gradPrimtiveVarZ[IW][le];

                    //! Gradient correction.
                    dx  = xcc[re] - xcc[le];
                    dy  = ycc[re] - ycc[le];
                    dz  = zcc[re] - zcc[le];
                    ods = 1.0 / sqrt(dx * dx + dy * dy + dz * dz);
                    dx  = dx * ods;
                    dy  = dy * ods;
                    dz  = dz * ods;

                    CorrectGradient(q[IU][le],q[IU][re],dudx,dudy,dudz,dx,dy,dz,ods);
                    CorrectGradient(q[IV][le],q[IV][re],dvdx,dvdy,dvdz,dx,dy,dz,ods);
                    CorrectGradient(q[IW][le],q[IW][re],dwdx,dwdy,dwdz,dx,dy,dz,ods);

                    //! It is unnecessary to consider because the turbulence viscosity coefficient is 0.
                    vis = visl[le];

                    const RDouble c23 = 2.0 / 3.0;
                    txx = vis * c23 * (2.0 * dudx - dvdy - dwdz);
                    tyy = vis * c23 * (2.0 * dvdy - dwdz - dudx);
                    tzz = vis * c23 * (2.0 * dwdz - dudx - dvdy);
                    txy = vis * (dudy + dvdx);
                    txz = vis * (dudz + dwdx);
                    tyz = vis * (dvdz + dwdy);
                    tyx = txy;
                    tzx = txz;
                    tzy = tyz;

                    fvsx = - two * (nx * txx + ny * tyx + nz * tzx) * oRefReNumber;
                    fvsy = - two * (nx * txy + ny * tyy + nz * tzy) * oRefReNumber;
                    fvsz = - two * (nx * txz + ny * tyz + nz * tzz) * oRefReNumber;

                    CA_f += fvsx;
                    CN_f += fvsy;
                    CZ_f += fvsz;

                    dfx += fvsx;
                    dfy += fvsy;
                    dfz += fvsz;

                    Cl_f += (yc - TorqueRefY) * fvsz - (zc - TorqueRefZ) * fvsy;
                    Cn_f += (zc - TorqueRefZ) * fvsx - (xc - TorqueRefX) * fvsz;
                    Cm_f += (xc - TorqueRefX) * fvsy - (yc - TorqueRefY) * fvsx;
                }

                Cl_p += (yc - TorqueRefY) * dpz - (zc - TorqueRefZ) * dpy;
                Cn_p += (zc - TorqueRefZ) * dpx - (xc - TorqueRefX) * dpz;
                Cm_p += (xc - TorqueRefX) * dpy - (yc - TorqueRefY) * dpx;

                cpx += dpx;
                cpy += dpy;
                cpz += dpz;

                if(dumpHingeMoment)
                {
                    Post_ForceMoment forceMoment;
                    RDouble point[3] = {xfc[iFace], yfc[iFace], zfc[iFace]};
                    RDouble faceForce[3] = {dfx, dfy, dfz};
                    hingeMoment += forceMoment.ComputeMoment(point, localCoordAxis0, localCoordAxis1, faceForce);
                }
            }
        }
    }

    DataContainer *cdata = actkey->GetData();

    cdata->write(&CA_f,sizeof(RDouble));
    cdata->write(&CA_p,sizeof(RDouble));
    cdata->write(&CN_f,sizeof(RDouble));
    cdata->write(&CN_p,sizeof(RDouble));
    cdata->write(&CZ_f,sizeof(RDouble));
    cdata->write(&CZ_p,sizeof(RDouble));

    cdata->write(&cpx,sizeof(RDouble));
    cdata->write(&cpy,sizeof(RDouble));
    cdata->write(&cpz,sizeof(RDouble));

    cdata->write(&Cl_f,sizeof(RDouble));
    cdata->write(&Cl_p,sizeof(RDouble));
    cdata->write(&Cn_f,sizeof(RDouble));
    cdata->write(&Cn_p,sizeof(RDouble));
    cdata->write(&Cm_f,sizeof(RDouble));
    cdata->write(&Cm_p,sizeof(RDouble));

    cdata->write(&hingeMoment, sizeof(RDouble));
}

void NSSolverUnstruct::AirForceCoefComponents(ActionKey *actkey, int partID, int Coordinate)
{
    int level = actkey->level;

    using namespace PHMPI;

    UnstructGrid *grid = UnstructGridCast(GetGrid(level));

    RDouble** q = reinterpret_cast<RDouble**> (grid->GetDataPtr("q"));
    RDouble* visl = reinterpret_cast<RDouble*> (grid->GetDataPtr("visl"));
    RDouble* vistl = reinterpret_cast<RDouble*> (grid->GetDataPtr("vist"));
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    RDouble oRefReNumber = parameters->GetoRefReNumber();

    int viscousType = parameters->GetViscousType();

    const int ALL_PART = -1;
    SimpleBC* boundaryCondition = 0;

    RDouble TorqueRefXGlobal = 0.0, TorqueRefYGlobal = 0.0, TorqueRefZGlobal = 0.0;
    RDouble TorqueRefXLocal = 0.0, TorqueRefYLocal = 0.0, TorqueRefZLocal = 0.0;

    TorqueRefXGlobal = GlobalDataBase::GetDoubleParaFromDB("TorqueRefX");
    TorqueRefYGlobal = GlobalDataBase::GetDoubleParaFromDB("TorqueRefY");
    TorqueRefZGlobal = GlobalDataBase::GetDoubleParaFromDB("TorqueRefZ");

    //! Used to compute hinge moment:
    int dumpHingeMoment = 0;
    RDouble localCoordAxis0[3];
    RDouble localCoordAxis1[3];
    RDouble hingeMoment = 0.0;

    int vis_run = 0;
    if (viscousType != INVISCID)
    {
        vis_run = 1;
    }

    using namespace IDX;
    RDouble* primitiveVarFarfield = parameters->GetPrimitiveVarFarfield();
    RDouble poo = primitiveVarFarfield[IP];

    RDouble* xfn = grid->GetFaceNormalX();
    RDouble* yfn = grid->GetFaceNormalY();
    RDouble* zfn = grid->GetFaceNormalZ();
    RDouble* area = grid->GetFaceArea();
    RDouble* vol = grid->GetCellVolume();

    RDouble* xfc = grid->GetFaceCenterX();
    RDouble* yfc = grid->GetFaceCenterY();
    RDouble* zfc = grid->GetFaceCenterZ();
    RDouble* xcc = grid->GetCellCenterX();
    RDouble* ycc = grid->GetCellCenterY();
    RDouble* zcc = grid->GetCellCenterZ();

    int* leftCellofFace = grid->GetLeftCellOfFace();
    int* rightCellofFace = grid->GetRightCellOfFace();

    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    RDouble cpx, cpy, cpz;
    cpx = 0.0;
    cpy = 0.0;
    cpz = 0.0;

    RDouble pw, cp;
    RDouble dfx, dfy, dfz, dpx, dpy, dpz, fvsx, fvsy, fvsz;
    RDouble nx, ny, nz;
    RDouble vis, txx, txy, txz, tyx, tyy, tyz, tzx, tzy, tzz;
    RDouble dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz;
    RDouble xc, yc, zc, dx, dy, dz, ods;

    RDouble **gradPrimtiveVarX = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarX"));
    RDouble **gradPrimtiveVarY = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarY"));
    RDouble **gradPrimtiveVarZ = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarZ"));

    using namespace IDX;

    int nsolid_surface = grid->GetNumberOfWallCell();

    //! Get the BCRegion information of unstruct grid.
    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegion = unstructBCSet->GetnBCRegion();

    int nTotalBC = GlobalBoundaryCondition::GetNumberOfBC();

    int numberOfAerodynamicForceComponents = GlobalDataBase::GetIntParaFromDB("numberOfAerodynamicForceComponents");

    vector< BasicAerodynamicForce* >* basicAerodynamicForceVector = CreateBasicAerodynamicForceVector();

    //! Cycle all bcregions and find the solid surface type of boundary.
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();
        string bcName = bcRegion->GetBCName();

        BasicAerodynamicForce tmpBasicAerodynamicForce;
        // if the bc type is wall, calculate the wall coefficient.
        if (bcType == PHENGLEI::SOLID_SURFACE)
        {
            for (int iPart = 0; iPart < nTotalBC; ++iPart)
            {
                boundaryCondition = GlobalBoundaryCondition::GetBC(iPart);
                string partName = boundaryCondition->GetBCName();
                if (partName != bcName)
                {
                    continue;
                }
                for (int iComponents = 0; iComponents < numberOfAerodynamicForceComponents; ++iComponents)
                {
                    if (FaceInAerodynamicForceComponents(bcName, iComponents))
                    {
                        tmpBasicAerodynamicForce.SetMomentReferencePoint(TorqueRefXGlobal, TorqueRefYGlobal, TorqueRefZGlobal);

                        vector<int>* faceIndex = bcRegion->GetFaceIndex();
                        int iFace, le, re;
                        for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
                        {
                            iFace = *iter;
                            le = leftCellofFace[iFace];
                            re = rightCellofFace[iFace];

                            nx = xfn[iFace] * area[iFace];
                            ny = yfn[iFace] * area[iFace];
                            nz = zfn[iFace] * area[iFace];

                            dx = xfc[iFace] - xcc[le];
                            dy = yfc[iFace] - ycc[le];
                            dz = zfc[iFace] - zcc[le];

                            //! Pressure drag
                            pw = q[IP][le];
                            cp = two * (pw - poo);

                            dpx = nx * cp;
                            dpy = ny * cp;
                            dpz = nz * cp;

                            tmpBasicAerodynamicForce.SetPressureAerodynamicForce(dpx, dpy, dpz);

                            dfx = dpx;
                            dfy = dpy;
                            dfz = dpz;

                            xc = xfc[iFace];
                            yc = yfc[iFace];
                            zc = zfc[iFace];

                            tmpBasicAerodynamicForce.SetPressureAerodynamicMoment(xc, yc, zc);

                            if (viscousType > INVISCID)
                            {
                                dudx = gradPrimtiveVarX[IU][le];
                                dudy = gradPrimtiveVarY[IU][le];
                                dudz = gradPrimtiveVarZ[IU][le];

                                dvdx = gradPrimtiveVarX[IV][le];
                                dvdy = gradPrimtiveVarY[IV][le];
                                dvdz = gradPrimtiveVarZ[IV][le];

                                dwdx = gradPrimtiveVarX[IW][le];
                                dwdy = gradPrimtiveVarY[IW][le];
                                dwdz = gradPrimtiveVarZ[IW][le];

                                //! Gradient correction.
                                dx = xcc[re] - xcc[le];
                                dy = ycc[re] - ycc[le];
                                dz = zcc[re] - zcc[le];
                                ods = 1.0 / sqrt(dx * dx + dy * dy + dz * dz);
                                dx = dx * ods;
                                dy = dy * ods;
                                dz = dz * ods;

                                CorrectGradient(q[IU][le], q[IU][re], dudx, dudy, dudz, dx, dy, dz, ods);
                                CorrectGradient(q[IV][le], q[IV][re], dvdx, dvdy, dvdz, dx, dy, dz, ods);
                                CorrectGradient(q[IW][le], q[IW][re], dwdx, dwdy, dwdz, dx, dy, dz, ods);

                                //! It is unnecessary to consider because the turbulence viscosity coefficient is 0.
                                vis = visl[le];

                                const RDouble c23 = 2.0 / 3.0;
                                txx = vis * c23 * (2.0 * dudx - dvdy - dwdz);
                                tyy = vis * c23 * (2.0 * dvdy - dwdz - dudx);
                                tzz = vis * c23 * (2.0 * dwdz - dudx - dvdy);
                                txy = vis * (dudy + dvdx);
                                txz = vis * (dudz + dwdx);
                                tyz = vis * (dvdz + dwdy);
                                tyx = txy;
                                tzx = txz;
                                tzy = tyz;

                                fvsx = -two * (nx * txx + ny * tyx + nz * tzx) * oRefReNumber;
                                fvsy = -two * (nx * txy + ny * tyy + nz * tzy) * oRefReNumber;
                                fvsz = -two * (nx * txz + ny * tyz + nz * tzz) * oRefReNumber;

                                tmpBasicAerodynamicForce.SetViscousAerodynamicForce(fvsx, fvsy, fvsz);

                                tmpBasicAerodynamicForce.SetViscousAerodynamicMoment(xc, yc, zc);
                            }
                            tmpBasicAerodynamicForce.ComputeResultantAerodynamicForce();

                            tmpBasicAerodynamicForce.ComputeResultantAerodynamicMoment();

                            (*basicAerodynamicForceVector)[iComponents]->AddAerodynamicForce(&tmpBasicAerodynamicForce);
                        }
                        (*basicAerodynamicForceVector)[numberOfAerodynamicForceComponents]->AddAerodynamicForce(&tmpBasicAerodynamicForce);

                    }
                }
            }
        }
    }

    PHSPACE::WriteAerodynamicForceToActionKey(actkey, basicAerodynamicForceVector);
    PHSPACE::FreeBasicAerodynamicForceVector(basicAerodynamicForceVector);
}

void NSSolverUnstruct::CpDistriCoef(ActionKey *actkey)
{
    int level = actkey->level;
    DataContainer *cdata = actkey->GetData();
    cdata->moveto_begin();

    UnstructGrid *grid = UnstructGridCast(GetGrid(level));
    int GridID = grid->GetGridID()->GetIndex();
    PHWrite(cdata, GridID);

    int nsolid_surface = grid->GetNumberOfWallCell();
    if ((nsolid_surface == 0) || (level != 0))
    {
        return;
    }

    using namespace PHMPI;
    using namespace PHENGLEI;
    using namespace IDX;

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    RDouble **q   = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    RDouble *visl = reinterpret_cast< RDouble * > (grid->GetDataPtr("visl"));

    RDouble *primitiveVarFarfield = parameters->GetPrimitiveVarFarfield();
    RDouble poo = primitiveVarFarfield[IP];

    RDouble oRefReNumber = parameters->GetoRefReNumber();
    RDouble refReNumber = parameters->GetRefReNumber();
    RDouble refMachNumber = parameters->GetRefMachNumber();

    int viscousType = parameters->GetViscousType();

    //! To obtain the reference values.
    RDouble refDimensionalDensity = parameters->GetRefDimensionalDensity();
    RDouble refDimensionalSonicSpeed = parameters->GetRefDimensionalSonicSpeed();
    RDouble refVelocity = refDimensionalSonicSpeed * refMachNumber;

    RDouble TorqueRefX, TorqueRefY, TorqueRefZ;
    TorqueRefX = GlobalDataBase::GetDoubleParaFromDB("TorqueRefX");
    TorqueRefY = GlobalDataBase::GetDoubleParaFromDB("TorqueRefY");
    TorqueRefZ = GlobalDataBase::GetDoubleParaFromDB("TorqueRefZ");

    RDouble *x = grid->GetX();
    RDouble *y = grid->GetY();
    RDouble *z = grid->GetZ();

    RDouble *xfn  = grid->GetFaceNormalX();
    RDouble *yfn  = grid->GetFaceNormalY();
    RDouble *zfn  = grid->GetFaceNormalZ();
    RDouble *area = grid->GetFaceArea() ;

    RDouble *xfc = grid->GetFaceCenterX();
    RDouble *yfc = grid->GetFaceCenterY();
    RDouble *zfc = grid->GetFaceCenterZ();
    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();
    int *node_number_of_each_face  = grid->GetNodeNumberOfEachFace();
    int *nodePosi  = grid->GetFace2NodeSubscript();
    int *face2node = grid->GetFace2Node();

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegion = unstructBCSet->GetnBCRegion();

    int nTotalCell = grid->GetNTotalCell();
    int nTotalNode = grid->GetNTotalNode();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    RDouble cfx, cfy, cfz, cpx, cpy, cpz, cmx, cmy, cmz;
    cfx = 0.0;
    cfy = 0.0;
    cfz = 0.0;
    cpx = 0.0;
    cpy = 0.0;
    cpz = 0.0;
    cmx = 0.0;
    cmy = 0.0;
    cmz = 0.0;

    RDouble AoA = 0, attack = 0, angleSlide = 0, angleSlided = 0;
    AoA = parameters->GetAoA();
    attack   = AoA * PI / 180.0;
    angleSlided = parameters->GetAngleOfSlide();
    angleSlide = angleSlided * PI / 180.0;

    RDouble sina = 0, cosa = 0, sinb = 0, cosb = 0;
    sina = sin(attack);
    cosa = cos(attack);
    sinb = sin(angleSlide);
    cosb = cos(angleSlide);

    RDouble **gradPrimtiveVarX = reinterpret_cast< RDouble ** > (grid->GetDataPtr("gradPrimtiveVarX"));
    RDouble **gradPrimtiveVarY = reinterpret_cast< RDouble ** > (grid->GetDataPtr("gradPrimtiveVarY"));
    RDouble **gradPrimtiveVarZ = reinterpret_cast< RDouble ** > (grid->GetDataPtr("gradPrimtiveVarZ"));

    RDouble pw = 0, cp = 0, pw_Dim = 0;
    RDouble dfx = 0, dfy = 0, dfz = 0, dpx = 0, dpy = 0, dpz = 0, dmx = 0, dmy = 0, dmz = 0, fvsx = 0, fvsy = 0, fvsz = 0;
    RDouble nx = 0, ny = 0, nz = 0;
    RDouble normalComponetX = 0, normalComponetY = 0, normalComponetZ = 0;
    RDouble vis = 0, txx = 0, txy = 0, txz = 0, tyx = 0, tyy = 0, tyz = 0, tzx = 0, tzy = 0, tzz = 0;
    RDouble dudx = 0, dudy = 0, dudz = 0, dvdx = 0, dvdy = 0, dvdz = 0, dwdx = 0, dwdy = 0, dwdz = 0;
    RDouble xc = 0, yc = 0, zc = 0, dx = 0, dy = 0, dz = 0, ods = 0;

    RDouble *wallDistance = grid->GetWallDist();

    RDouble *faceCP    = new RDouble[nBoundFace]();
    RDouble *faceCF    = new RDouble[nBoundFace]();
    RDouble *faceQ     = new RDouble[nBoundFace]();
    RDouble *faceYPlus = new RDouble[nBoundFace]();
    RDouble *facePW    = new RDouble[nBoundFace]();

    //! Compute cp and cf on wall face.
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();
        if (IsWall(bcType))
        {
            vector<int> *faceIndex = bcRegion->GetFaceIndex();
            for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
            {
                //! iFace is the face number in the set of faceIndex.
                int iFace = *iter;
                int le = leftCellofFace[iFace];
                int re = rightCellofFace[iFace];

                normalComponetX = xfn[iFace];
                normalComponetY = yfn[iFace];
                normalComponetZ = zfn[iFace];

                nx  = normalComponetX * area[iFace];
                ny  = normalComponetY * area[iFace];
                nz  = normalComponetZ * area[iFace];

                //! Pressure drag.
                dx  = xfc[iFace] - xcc[le];
                dy  = yfc[iFace] - ycc[le];
                dz  = zfc[iFace] - zcc[le];
                pw = q[IP][le];
                cp = two * (pw - poo);

                //! Change the nondimensional pressure into dimensional pressure(Pa).
                pw_Dim = pw * refDimensionalDensity * pow(refVelocity, 2);

                //! Note that the definition of cp is right here and now, but dpx and dfx have multiplied by area.
                dpx = nx * cp;
                dpy = ny * cp;
                dpz = nz * cp;

                dfx = dpx;
                dfy = dpy;
                dfz = dpz;

                fvsx = 0.0;
                fvsy = 0.0;
                fvsz = 0.0;

                if (viscousType > INVISCID)
                {
                    dudx = gradPrimtiveVarX[IU][le];
                    dudy = gradPrimtiveVarY[IU][le];
                    dudz = gradPrimtiveVarZ[IU][le];

                    dvdx = gradPrimtiveVarX[IV][le];
                    dvdy = gradPrimtiveVarY[IV][le];
                    dvdz = gradPrimtiveVarZ[IV][le];

                    dwdx = gradPrimtiveVarX[IW][le];
                    dwdy = gradPrimtiveVarY[IW][le];
                    dwdz = gradPrimtiveVarZ[IW][le];

                    //! Gradient correction.
                    dx  = xcc[re] - xcc[le];
                    dy  = ycc[re] - ycc[le];
                    dz  = zcc[re] - zcc[le];
                    ods = 1.0 / sqrt(dx * dx + dy * dy + dz * dz);
                    dx  = dx * ods;
                    dy  = dy * ods;
                    dz  = dz * ods;

                    CorrectGradient(q[IU][le],q[IU][re],dudx,dudy,dudz,dx,dy,dz,ods);
                    CorrectGradient(q[IV][le],q[IV][re],dvdx,dvdy,dvdz,dx,dy,dz,ods);
                    CorrectGradient(q[IW][le],q[IW][re],dwdx,dwdy,dwdz,dx,dy,dz,ods);

                    //! It is unnecessary to consider because the turbulence viscosity coefficient is 0.
                    vis = visl[le];

                    txx = vis * two3rd * (2.0 * dudx - dvdy - dwdz);
                    tyy = vis * two3rd * (2.0 * dvdy - dwdz - dudx);
                    tzz = vis * two3rd * (2.0 * dwdz - dudx - dvdy);
                    txy = vis * (dudy + dvdx);
                    txz = vis * (dudz + dwdx);
                    tyz = vis * (dvdz + dwdy);
                    tyx = txy;
                    tzx = txz;
                    tzy = tyz;

                    fvsx = - two * (nx * txx + ny * tyx + nz * tzx) * oRefReNumber;
                    fvsy = - two * (nx * txy + ny * tyy + nz * tzy) * oRefReNumber;
                    fvsz = - two * (nx * txz + ny * tyz + nz * tzz) * oRefReNumber;

                    dfx += fvsx;
                    dfy += fvsy;
                    dfz += fvsz;
                }
                xc = xfc[iFace];
                yc = yfc[iFace];
                zc = zfc[iFace];

                dmx = (yc - TorqueRefY) * dfz - (zc - TorqueRefZ) * dfy;
                dmy = (zc - TorqueRefZ) * dfx - (xc - TorqueRefX) * dfz;
                dmz = (xc - TorqueRefX) * dfy - (yc - TorqueRefY) * dfx;

                cpx += dpx;
                cpy += dpy;
                cpz += dpz;

                cfx += dfx;
                cfy += dfy;
                cfz += dfz;

                cmx += dmx;
                cmy += dmy;
                cmz += dmz;

                RDouble dir = cosa * fvsx + sina * fvsy;

                RDouble density = q[IR][le];

                //! Compute taoWall.
                RDouble cfn  = fvsx * normalComponetX + fvsy * normalComponetY + fvsz * normalComponetZ;
                RDouble cfxt = fvsx - cfn * normalComponetX;
                RDouble cfyt = fvsy - cfn * normalComponetY;
                RDouble cfzt = fvsz - cfn * normalComponetZ;

                RDouble cft = sqrt(cfxt * cfxt + cfyt * cfyt + cfzt * cfzt);

                RDouble cf = cft / area[iFace] * SIGN(1.0,dir);
                RDouble taoWall = half * cft / area[iFace];

                //! Compute yplus.
                RDouble yPlus = 0.0;
                if (viscousType > INVISCID)
                {
                    yPlus = refReNumber * sqrt(taoWall * density) * wallDistance[le] / vis;
                }

                faceCP   [iFace] = cp;
                faceCF   [iFace] = cf;
                faceYPlus[iFace] = yPlus;
                facePW   [iFace] = pw_Dim;
            }
        }
    }

    //! Compute Q on wall face.
    RDouble wallTemperature = parameters->GetWallTemperature();
    RDouble dimensionalQ = 0.0;
    if (wallTemperature >= 0)
    {
        RDouble refGama = parameters->GetRefGama();
        RDouble refMachNumber = parameters->GetRefMachNumber();
        RDouble refDimensionalTemperature = parameters->GetRefDimensionalTemperature();

        RDouble refDimensionalSonicSpeed = parameters->GetRefDimensionalSonicSpeed();
        RDouble refDimensionalDensity = parameters->GetRefDimensionalDensity();
        RDouble velocityInflow = refMachNumber * refDimensionalSonicSpeed;
        dimensionalQ = refDimensionalDensity * velocityInflow * velocityInflow * velocityInflow * 0.001;    //! 0.001: w --> kw

        RDouble prandtlLaminar = parameters->GetPrandtlLaminar();
        RDouble nonDimensionalSutherlandTemperature = parameters->GetNonDimensionalSutherlandTemperature();

        RDouble Tw  = wallTemperature / refDimensionalTemperature;
        RDouble Cv  = 1.0 / refGama / (refGama - 1.0) / refMachNumber / refMachNumber;
        RDouble tmu = pow(Tw, 1.5) * (1.0 + nonDimensionalSutherlandTemperature) / (Tw + nonDimensionalSutherlandTemperature);

        RDouble **t = reinterpret_cast< RDouble ** > (grid->GetDataPtr("t"));


        for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
        {
            UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
            int bcType = bcRegion->GetBCType();

            if (IsWall(bcType))
            {
                vector<int> *faceIndex = bcRegion->GetFaceIndex();
                for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
                {
                    //! iFace is the face number in the set of faceIndex.
                    int iFace = *iter;
                    int le = leftCellofFace[iFace];
                    int re = rightCellofFace[iFace];

                    RDouble deltT = t[ITT][re] - t[ITT][le];
                    RDouble dx = xcc[re] - xcc[le];
                    RDouble dy = ycc[re] - ycc[le];
                    RDouble dz = zcc[re] - zcc[le];
                    RDouble dd = dx * dx + dy * dy + dz * dz;

                    dd = - sqrt(dd);
                    deltT /= dd;

                    faceQ[iFace] = Cv * tmu * deltT * refGama * oRefReNumber / prandtlLaminar;
                }
            }
        }
    }

    //! Communicate data on interface.
    int *lableCell = new int[nTotalCell];
    std::fill_n(lableCell, nTotalCell, 0);

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();

        if (IsWall(bcType))
        {
            vector<int> *faceIndex = bcRegion->GetFaceIndex();
            for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
            {
                int iFace = *iter;
                int le = leftCellofFace[iFace];
                lableCell[le] = 1;
            }
        }
        if (bcType == INTERFACE)
        {
            vector<int> *faceIndex = bcRegion->GetFaceIndex();
            for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
            {
                int iFace = *iter;
                int cL = leftCellofFace[iFace];
                if (lableCell[cL] == 1)
                {
                    lableCell[cL] = 2;
                }
            }
        }
    }

    RDouble **cellQ2D = NewPointer2<RDouble>(1, nTotal);
    RDouble * cellQ   = cellQ2D[0];
    std::fill_n(cellQ, nTotal, 0);

    RDouble **cellCP2D = NewPointer2<RDouble>(1, nTotal);
    RDouble * cellCP   = cellCP2D[0];
    std::fill_n(cellCP, nTotal, 0);

    RDouble **cellCF2D = NewPointer2<RDouble>(1, nTotal);
    RDouble * cellCF   = cellCF2D[0];
    std::fill_n(cellCF, nTotal, 0);

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();

        if (IsWall(bcType))
        {
            vector<int> *faceIndex = bcRegion->GetFaceIndex();
            for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
            {
                int iFace = *iter;
                int le = leftCellofFace[iFace];
                if(lableCell[le] == 2)
                {
                    cellCF[le] = faceCF[iFace];
                    cellCP[le] = faceCP[iFace];
                    if (wallTemperature >= 0)
                    {
                        cellQ[le] = faceQ[iFace];
                    }
                }
            }
        }
    }

    delete []lableCell;    lableCell = NULL;

    PHSPACE::CommunicateInterfaceValue(grid, cellCP2D, "cellCP", 1);
    PHSPACE::CommunicateInterfaceValue(grid, cellCF2D, "cellCF", 1);
    if (wallTemperature >= 0)
    {
        PHSPACE::CommunicateInterfaceValue(grid, cellQ2D, "cellQ", 1);
    }

    //! Compute data on wall node.
    int *nCell = new int[nTotalNode]();
    std::fill_n(nCell, nTotalNode, 0);

    RDouble *nodeQ = new RDouble[nTotalNode]();
    std::fill_n(nodeQ, nTotalNode, 0.0);

    RDouble *nodeCP = new RDouble[nTotalNode]();
    std::fill_n(nodeCP, nTotalNode, 0.0);

    RDouble *nodeCF = new RDouble[nTotalNode]();
    std::fill_n(nodeCF, nTotalNode, 0.0);

    RDouble *nodeYPlus = new RDouble[nTotalNode]();
    std::fill_n(nodeYPlus, nTotalNode, 0.0);

    RDouble *nodePW = new RDouble[nTotalNode]();
    std::fill_n(nodePW, nTotalNode, 0.0);

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();

        if (IsWall(bcType))
        {
            vector<int> *faceIndex = bcRegion->GetFaceIndex();
            for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
            {
                int iFace = *iter;
                for (int iNode = nodePosi[iFace]; iNode < nodePosi[iFace+1]; ++ iNode)
                {
                    int nodeIndex = face2node[iNode];
                    nodeCP   [nodeIndex] += faceCP[iFace];
                    nodeCF   [nodeIndex] += faceCF[iFace];
                    nodeYPlus[nodeIndex] += faceYPlus[iFace];
                    nodePW   [nodeIndex] += facePW[iFace];
                    if (wallTemperature >= 0)
                    {
                        nodeQ[nodeIndex] += faceQ[iFace];
                    }
                    nCell[nodeIndex] ++;
                }
            }
        }
        if (bcType == INTERFACE)
        {
            vector<int> *faceIndex = bcRegion->GetFaceIndex();
            for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
            {
                int iFace = *iter;
                int cR = rightCellofFace[iFace];
                if (ABS(cellQ[cR]) > SMALL)
                {
                    for (int iNode = nodePosi[iFace]; iNode < nodePosi[iFace+1]; ++ iNode)
                    {
                        int nodeIndex = face2node[iNode];
                        if (nCell[nodeIndex] != 0)
                        {
                            nodeCP[nodeIndex] += cellCP[cR];
                            nodeCF[nodeIndex] += cellCF[cR];
                            if (wallTemperature >= 0)
                            {
                                nodeQ[nodeIndex] += cellQ[cR];
                            }
                            nCell[nodeIndex] ++;
                        }
                    }
                }
            }
        }
    }

    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        if (nCell[iNode])
        {
            nodeYPlus[iNode] /= nCell[iNode];
            nodePW   [iNode] /= nCell[iNode];
        }
    }

    //! For interface.
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();

        if (bcType == INTERFACE)
        {
            vector<int> *faceIndex = bcRegion->GetFaceIndex();
            for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
            {
                int iFace = *iter;
                int cR = rightCellofFace[iFace];
                if (ABS(cellQ[cR]) > SMALL)
                {
                    for (int iNode = nodePosi[iFace]; iNode < nodePosi[iFace+1]; ++ iNode)
                    {
                        int nodeIndex = face2node[iNode];
                        if (nCell[nodeIndex] != 0)
                        {
                            nodeCP[nodeIndex] += cellCP[cR];
                            nodeCF[nodeIndex] += cellCF[cR];
                            if (wallTemperature >= 0)
                            {
                                nodeQ[nodeIndex] += cellQ[cR];
                            }
                            nCell[nodeIndex] ++;
                        }
                    }
                }
            }
        }
    }

    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        if (nCell[iNode])
        {
            nodeCP[iNode] /= nCell[iNode];
            nodeCF[iNode] /= nCell[iNode];
            if (wallTemperature >= 0)
            {
                nodeQ[iNode] /= nCell[iNode];
            }
        }
    }

    delete []faceCP;    faceCP = NULL;
    delete []faceCF;    faceCF = NULL;
    delete []faceQ;    faceQ = NULL;
    delete []nCell;    nCell = NULL;
    delete []faceYPlus;    faceYPlus = NULL;
    delete []facePW;    facePW = NULL;
    DelPointer2(cellQ2D);
    DelPointer2(cellCP2D);
    DelPointer2(cellCF2D);

    //! Compute data on wall face center.
    RDouble *cpFaceCenter    = new RDouble[nBoundFace]();
    RDouble *cfFaceCenter    = new RDouble[nBoundFace]();
    RDouble *qFaceCenter     = new RDouble[nBoundFace]();
    RDouble *yPlusFaceCenter = new RDouble[nBoundFace]();
    RDouble *PWFaceCenter    = new RDouble[nBoundFace]();

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        if (face2node[iFace] > 4)
        {
            RDouble valueCP    = 0.0;
            RDouble valueCF    = 0.0;
            RDouble valueQ     = 0.0;
            RDouble valueYPlus = 0.0;
            RDouble valuePW    = 0.0;

            int nodeIndex;
            for (int iNode = nodePosi[iFace]; iNode < nodePosi[iFace+1]; ++ iNode)
            {
                nodeIndex   = face2node[iNode];
                valueCP    += nodeCP   [nodeIndex];
                valueCF    += nodeCF   [nodeIndex];
                valueYPlus += nodeYPlus[nodeIndex];
                valuePW    += nodePW   [nodeIndex];
                if (wallTemperature >= 0)
                {
                    valueQ += nodeQ[nodeIndex];
                }
            }
            cpFaceCenter[iFace]    = valueCP    / node_number_of_each_face[iFace];
            cfFaceCenter[iFace]    = valueCF    / node_number_of_each_face[iFace];
            qFaceCenter [iFace]    = valueQ     / node_number_of_each_face[iFace];
            yPlusFaceCenter[iFace] = valueYPlus / node_number_of_each_face[iFace];
            PWFaceCenter[iFace] = valuePW / node_number_of_each_face[iFace];
        }
    }

    //! Reconstruct grid mesh of wall.
    vector < vector < int > > face2nodelist(nBoundFace);
    vector < int > linkmap;
    GetFace2NodeList(grid, SOLID_SURFACE, linkmap, face2nodelist);

    int NumPts      = static_cast<int>(linkmap.size());
    int NumElements = static_cast<int>(face2nodelist.size());
    int NumFaces    = NumElements * 4;
    int TotalNumFaceNodes = NumElements * 4;
    int *cell2node        = new int[TotalNumFaceNodes]();

    int count = 0;
    for (std::size_t iFace = 0; iFace < face2nodelist.size(); ++ iFace)
    {
        uint_t nodes  = face2nodelist[iFace].size();
        int nodeIndex = face2nodelist[iFace][0];
        for (int iNode = 0; iNode < nodes; ++ iNode)
        {
            cell2node[count] = face2nodelist[iFace][iNode] + 1;
            count ++;
        }

        for (uint_t ip = nodes; ip < 4; ++ ip)
        {
            cell2node[count] = nodeIndex + 1;
            count ++;
        }
    }

    RDouble *coorX     = new RDouble[NumPts]();
    RDouble *coorY     = new RDouble[NumPts]();
    RDouble *coorZ     = new RDouble[NumPts]();
    RDouble *dumpCP    = new RDouble[NumPts]();
    RDouble *dumpCF    = new RDouble[NumPts]();
    RDouble *dumpQ     = new RDouble[NumPts]();
    RDouble *dumpQNon  = new RDouble[NumPts]();
    RDouble *dumpYPlus = new RDouble[NumPts]();
    RDouble *dumpPW    = new RDouble[NumPts]();

    for (std::size_t iNode = 0; iNode < linkmap.size(); ++ iNode)
    {
        int nodeIndex = linkmap[iNode];

        if (nodeIndex < nTotalNode)
        {
            coorX[iNode] = x[nodeIndex];
            coorY[iNode] = y[nodeIndex];
            coorZ[iNode] = z[nodeIndex];

            dumpCP   [iNode] = nodeCP   [nodeIndex];
            dumpCF   [iNode] = nodeCF   [nodeIndex];
            dumpYPlus[iNode] = nodeYPlus[nodeIndex];
            dumpPW   [iNode] = nodePW   [nodeIndex];
            dumpQ    [iNode] = nodeQ    [nodeIndex] * dimensionalQ;
            dumpQNon [iNode] = nodeQ    [nodeIndex];
        }
        else
        {
            int faceIndex = nodeIndex - nTotalNode;
            coorX[iNode] = xfc[faceIndex];
            coorY[iNode] = yfc[faceIndex];
            coorZ[iNode] = zfc[faceIndex];

            dumpCP   [iNode] = cpFaceCenter   [faceIndex];
            dumpCF   [iNode] = cfFaceCenter   [faceIndex];
            dumpYPlus[iNode] = yPlusFaceCenter[faceIndex];
            dumpPW   [iNode] = PWFaceCenter   [faceIndex];
            dumpQ    [iNode] = qFaceCenter    [faceIndex] * dimensionalQ;
            dumpQNon [iNode] = qFaceCenter    [faceIndex];
        }
    }

    delete []nodeQ;    nodeQ = NULL;
    delete []nodeCP;    nodeCP = NULL;
    delete []nodeCF;    nodeCF = NULL;
    delete []nodeYPlus;    nodeYPlus = NULL;
    delete []nodePW;    nodePW = NULL;
    delete []cpFaceCenter;    cpFaceCenter = NULL;
    delete []cfFaceCenter;    cfFaceCenter = NULL;
    delete []qFaceCenter;    qFaceCenter = NULL;
    delete []yPlusFaceCenter;    yPlusFaceCenter = NULL;
    delete []PWFaceCenter;    PWFaceCenter = NULL;

    int TecioMission = 1;
    PHWrite(cdata, &TecioMission, 1);

    PHWrite(cdata, &NumPts, 1);
    PHWrite(cdata, &NumElements, 1);
    PHWrite(cdata, &NumFaces, 1);

    PHWrite(cdata, NumPts);
    PHWrite(cdata, coorX, NumPts);
    PHWrite(cdata, coorY, NumPts);
    PHWrite(cdata, coorZ, NumPts);

    int ValueLocation = 1;

    PHWrite(cdata, NumPts);
    PHWrite(cdata, ValueLocation);
    PHWrite(cdata, dumpCP, NumPts);

    PHWrite(cdata, NumPts);
    PHWrite(cdata, ValueLocation);
    PHWrite(cdata, dumpCF, NumPts);

    PHWrite(cdata, NumPts);
    PHWrite(cdata, ValueLocation);
    PHWrite(cdata, dumpYPlus, NumPts);

    ValueLocation = 1;

    PHWrite(cdata, NumPts);
    PHWrite(cdata, ValueLocation);
    PHWrite(cdata, dumpQNon, NumPts);

    PHWrite(cdata, NumPts);
    PHWrite(cdata, ValueLocation);
    PHWrite(cdata, dumpQ, NumPts);

    PHWrite(cdata, NumPts);
    PHWrite(cdata, ValueLocation);
    PHWrite(cdata, dumpPW, NumPts);

    PHWrite(cdata, &TotalNumFaceNodes, 1);
    PHWrite(cdata, cell2node, TotalNumFaceNodes);

    delete []coorX;    coorX = NULL;
    delete []coorY;    coorY = NULL;
    delete []coorZ;    coorZ = NULL;
    delete []dumpCP;    dumpCP = NULL;
    delete []dumpCF;    dumpCF = NULL;
    delete []dumpQ;    dumpQ = NULL;
    delete []dumpQNon;    dumpQNon = NULL;
    delete []dumpYPlus;    dumpYPlus = NULL;
    delete []dumpPW;    dumpPW = NULL;
    delete []cell2node;    cell2node = NULL;
}

//! Dump out some flow field variables, for the case of turbulent flat plate of NASA Turbulence Center.
//! Inflow condition: Ma = 0.2, Re = 5.0e6.
void NSSolverUnstruct::Turbulence_Flat_Plate_Output_NASA(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalNode = grid->GetNTotalNode();
    RDouble ** qVist;

    //! 100.000 kPa = 1 bar ≈ 14.504 psi
    //! Condition for calculation
    //! Static Pressure (psia) 6.0 = 100000 / 14.504 * 6 Pa = 41367.8986 Pa
    //! Temperature(R) 700.0 = (5/9)*700 = 388.88K

    //! Dump out the Velocity.
    vector<string> title_tecplot;
    title_tecplot.push_back("title=\"Flow Fields of PHengLEI\"");
    title_tecplot.push_back("variables=");
    title_tecplot.push_back("\"y\"");
    title_tecplot.push_back("\"u\"");
    title_tecplot.push_back("\"y+\"");
    title_tecplot.push_back("\"u+\"");    
    int nvarplot;
    nvarplot = 2;

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    using namespace PHSPACE;
    RDouble refGama = parameters->GetRefGama();
    RDouble reference_mach_number = parameters->GetRefMachNumber();
    RDouble refReNumber = 1.0e5;
    refReNumber = parameters->GetRefReNumber();
    RDouble tsuth;
    GlobalDataBase::GetData("tsuth", &tsuth, PHDOUBLE, 1);

    fstream file;
    PHSPACE::OpenFile(file, "results/turblent_flat_plate.dat", ios_base::out);

    int wordwidth = 20;

    RDouble * vist = reinterpret_cast< RDouble * > (grid->GetDataPtr("vist"));

    int nline = 2;

    vector<RDouble> xPositionList(nline);
    xPositionList[0] = 0.97008;
    xPositionList[1] = 1.90334;    

    for (std::size_t i = 0; i < title_tecplot.size(); ++ i)
    {
        file << title_tecplot[i] << "\n";
    }

    RDouble ** nodevar = NewPointer2<RDouble>(nEquation + 1, nTotalNode);
    CreateNodeVariable(grid, nodevar);

    for (int m = 0; m < nline; ++ m)
    {
        PHCutPlane * cutpl = new PHCutPlane();
        CutPlane(grid, nodevar, nEquation, xPositionList[m], cutpl, X_DIR);

        cutpl->SortData(cutpl->y);

        int nT = cutpl->nTotalNode;
        file << "zone  T = \"x = " << xPositionList[m] << "\" \n";
        vector<RDouble> ylist(nT);
        for (int i = 0; i < nT; ++ i)
        {
            ylist[i] = cutpl->y[i];
        }
        sort(ylist.begin(),ylist.end());

        int yindex = 1;
        for (int i = 0; i < nT; ++ i)
        {
            if (ylist[i])
            {
                yindex = i;
                break;
            }
        }

        for (int i = 0; i < nT; ++ i)
        {
            if (ylist[yindex] == cutpl->y[i])
            {
                yindex = i;
                break;
            }
        }

        for (int i = 0; i < nT; ++ i)
        {
            RDouble rm,um,vm,wm,pm,tm,xm,ym,zm,mu,muw,dudy,tauw,utau,up,yp,rw;
            RDouble uw,vw,ww,pw,yw,tw;
            xm   = cutpl->x[i];
            ym   = cutpl->y[i];
            zm   = cutpl->z[i];

            rm   = cutpl->data[0][i];
            um   = cutpl->data[1][i];
            vm   = cutpl->data[2][i];
            wm   = cutpl->data[3][i];
            pm   = cutpl->data[4][i];
            tm   = refGama * reference_mach_number * reference_mach_number * pm / rm;
            mu   = tm * sqrt(tm) * (1.0 + tsuth) / (tm + tsuth);

            rw   = cutpl->data[0][yindex];
            uw   = cutpl->data[1][yindex];
            vw   = cutpl->data[2][yindex];
            ww   = cutpl->data[3][yindex];
            pw   = cutpl->data[4][yindex];
            tw   = refGama * reference_mach_number * reference_mach_number * pm / rm;

            muw  = tm * sqrt(tm) * (1.0 + tsuth) / (tm + tsuth);
            yw   = cutpl->y[yindex];

            dudy = uw / yw;
            tauw = muw * dudy;

            utau = sqrt(tauw/(rw * refReNumber));

            up   = um / utau;
            yp   = utau * ym * refReNumber / (mu / rm);

            file << setiosflags(ios::left);
            file << setiosflags(ios::scientific);
            file << setprecision(10);
            file << setw(wordwidth) << ym 
                 << setw(wordwidth) << um
                 << setw(wordwidth) << yp
                 << setw(wordwidth) << up;
            file << "\n";
        }
        delete cutpl;
    }
    file.close();
    file.clear();
    DelPointer2(nodevar);

    //! Dump out the friction.
    RDouble **q   = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q" ));
    PHCutPlane * cutpl = new PHCutPlane();
    CutPlaneCellCenter(grid, q, nEquation, 2.0e-6, cutpl, Y_DIR);
    cutpl->SortData(cutpl->x);

    int nT = cutpl->nTotalNode;
    if(nT > 0)
    {
        ostringstream frictionFileName;
        frictionFileName << "results/turblent_flat_plate_nasa_cf_" << grid->GetZoneID() + 1<< ".dat";
        PHSPACE::OpenFile(file, frictionFileName.str().c_str(), ios_base::out);

        file << "title=\"Flow Fields of PHengLEI\"\n";
        file << "variables = \"x\" \"cf\" \n";
        file << "zone  i = " << nT  << " \n";
        for (int i = 0; i < nT; ++ i)
        {
            RDouble rm,um,vm,wm,pm,tm,xm,ym,zm,mu,dudy;
            xm   = cutpl->x[i];
            ym   = cutpl->y[i];
            zm   = cutpl->z[i];

            rm   = cutpl->data[0][i];
            um   = cutpl->data[1][i];
            vm   = cutpl->data[2][i];
            wm   = cutpl->data[3][i];
            pm   = cutpl->data[4][i];
            tm   = refGama * reference_mach_number * reference_mach_number * pm / rm;
            mu   = tm * sqrt(tm) * (1.0 + tsuth) / (tm + tsuth);
            dudy = um / ym;

            file << setiosflags(ios::left);
            file << setiosflags(ios::scientific);
            file << setprecision(10);
            file << setw(wordwidth) << xm
                << setw(wordwidth) << 2 * mu * dudy / refReNumber
                << "\n";
        }

        frictionFileName.clear();
        
        PHSPACE::CloseFile(file);
    }

    delete cutpl;

    //! Dump out the turbulent viscosity of x = 0.97
    qVist = new RDouble *[1];    //! Tempal 2D array for vist, to dump both x and y value.
    qVist[0] = vist;
    cutpl = new PHCutPlane();
    CutPlaneCellCenter(grid,qVist,1,0.97,cutpl,X_DIR);
    cutpl->SortData(cutpl->y);
    nT = cutpl->nTotalNode;
    if(nT > 0)
    {
        ostringstream xSliceFileName;
        xSliceFileName << "results/turblent_flat_plate_nasa_vist_x0d97_" << grid->GetZoneID() + 1 << ".dat";
        PHSPACE::OpenFile(file, xSliceFileName.str().c_str(), ios_base::out);

        file << "title=\"turbulent viscosity at x = 0.97\"\n";
        file << "variables = \"vist\"    \"Y\"\n";
        file << "zone  i = " << nT  << " \n";
        for (int i = 0; i < nT; ++ i)
        {
            RDouble ym   = cutpl->y[i];
            RDouble vist = cutpl->data[0][i];

            file << setiosflags(ios::left);
            file << setiosflags(ios::scientific);
            file << setprecision(10);
            file << setw(wordwidth) << vist
                << setw(wordwidth) << ym
                << "\n";
        }

        xSliceFileName.clear();

        PHSPACE::CloseFile(file);
    }

    delete cutpl;
    delete [] qVist;

    //! Dump out the the MAXIMUM eddy viscosity as a function of x
    using namespace PHMPI;
    if (GetNumberofGlobalZones() > 1)
    {
        return;
    }
    string maxVistOfX = "results/turblent_flat_plate_nasa_vist_maximum.dat";
    PHSPACE::OpenFile(file, maxVistOfX.c_str(), ios_base::out);

    file << "title=\"the maximum turbulent viscosity\"\n";
    file << "variables = \"X\"    \"vist\"\n";

    qVist = new RDouble *[1];    //! Tempal 2D array for vist, to dump both x and y value.
    qVist[0] = vist;

    //! Get the number of x stand position.
    int nXStand = 0;
    int nBoundFace  = grid->GetNBoundFace();
    RDouble* xfc = grid->GetFaceCenterX();
   
    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegion = unstructBCSet->GetnBCRegion();
    
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC* bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();
        vector<int>* faceIndex = bcRegion->GetFaceIndex();
        if (bcType == PHENGLEI::SOLID_SURFACE || bcType == PHENGLEI::SYMMETRY)
        {
            for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
            {
                ++ nXStand;
            }

        }
    }

    //! Get the x stand position.
    RDouble* xOfStand = new RDouble [nXStand];
    nXStand = 0;

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC* bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();
        vector<int>* faceIndex = bcRegion->GetFaceIndex();
        if (bcType == PHENGLEI::SOLID_SURFACE || bcType == PHENGLEI::SYMMETRY)
        {
            for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
            {
                int iFace = *iter;
                xOfStand[nXStand++] = xfc[iFace];
            }

        }
    }

    vector< vector < RDouble > > var_list;
    vector < RDouble > x_list;
    vector < RDouble > var(1);
    cutpl = new PHCutPlane();    //! Store the maximum viscosity of each xcut stand.

    RDouble* maxVistOfStandX = new RDouble [nXStand];
    for (int iStand = 0; iStand < nXStand; ++ iStand)
    {
        RDouble xCut = xOfStand[iStand];
        PHCutPlane* cutpX = new PHCutPlane();    //! Store the viscosity of ONE x stand, x = xcut.

        RDouble maxVist = -LARGE;
        CutPlaneCellCenter(grid,qVist,1,xCut,cutpX,X_DIR);
        nT = cutpX->nTotalNode;
        for (int i = 0; i < nT; ++ i)
        {
            RDouble vist = cutpX->data[0][i];
            maxVist = MAX(maxVist, vist);
        }

        var[0] = maxVist;
        var_list.push_back(var);
        x_list.push_back(xCut);

        delete cutpX;
    }

    cutpl->neqn   = 1;
    cutpl->nTotalNode = static_cast<int>(var_list.size());
    cutpl->data   = new RDouble * [ 1 ];
    for (int m = 0; m < 1; ++ m)
    {
        cutpl->data[m] = new RDouble [ cutpl->nTotalNode ];
    }
    cutpl->x = new RDouble [ cutpl->nTotalNode ];
    cutpl->y = new RDouble [ cutpl->nTotalNode ];
    cutpl->z = new RDouble [ cutpl->nTotalNode ];
    for (int iNode = 0; iNode < cutpl->nTotalNode; ++ iNode)
    {
        cutpl->x[iNode] = x_list[iNode];
        cutpl->y[iNode] = 0;
        cutpl->z[iNode] = 0;
    }
    for (int iNode = 0; iNode < cutpl->nTotalNode; ++ iNode)
    {
        for (int m = 0; m < 1; ++ m)
        {
            cutpl->data[m][iNode] = var_list[iNode][m];
        }
    }

    cutpl->SortData(cutpl->x);

    nT = cutpl->nTotalNode;
    file << "zone  i = " << nT  << " \n";
    for (int i = 0; i < nT; ++ i)
    {
        RDouble xm   = cutpl->x[i];
        RDouble vist = cutpl->data[0][i];

        file << setiosflags(ios::left);
        file << setiosflags(ios::scientific);
        file << setprecision(10);
        file << setw(wordwidth) << xm
            << setw(wordwidth) << vist
            << "\n";
    }
    
    PHSPACE::CloseFile(file);
    delete cutpl;
    delete [] xOfStand;
    delete [] maxVistOfStandX;
    delete [] qVist;
}

void NSSolverUnstruct::CreateNodeVariable(Grid *gridIn, RDouble **nodevar)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    ComputeNodeVariable(grid, nodevar, q, nEquation);
}

void NSSolverUnstruct::ComputeAverageVisualNodeField(Grid *gridIn, RDouble **qn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    RDouble **qAverage = reinterpret_cast< RDouble ** > (grid->GetDataPtr("qAverage"));

    int varCount = 0;

    using namespace IDX;
    CompNodeVar(grid, qn[varCount ++], qAverage[IR]);
    CompNodeVarForVisual(grid, qn[varCount ++], qAverage[IU]);
    CompNodeVarForVisual(grid, qn[varCount ++], qAverage[IV]);
    CompNodeVarForVisual(grid, qn[varCount ++], qAverage[IW]);
    CompNodeVar(grid, qn[varCount ++], qAverage[IP]);

    RDouble * cp = ComputeCpAverage(gridIn);
    CompNodeVar(grid, qn[varCount ++], cp);
    delete [] cp;
}

void NSSolverUnstruct::SliceVisualization(Grid *gridIn, ActionKey *actkey, RDouble **qn, int nvarplot)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    DataContainer * cdata = actkey->GetTecData();

    int nBoundFace = grid->GetNBoundFace();
    int * node_number_of_each_face = grid->GetNodeNumberOfEachFace();
    int * face2node = grid->GetFace2Node();

    RDouble * x = grid->GetX();
    RDouble * y = grid->GetY();
    RDouble * z = grid->GetZ();

    vector < RDouble > xSlice;
    vector < RDouble > ySlice;
    vector < RDouble > zSlice;
    vector < vector <RDouble> > qnSlice;
    qnSlice.resize(nvarplot);

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int sliceAxis = parameters->GetSliceAxis();
    RDouble slicePostion = parameters->GetSlicePostion();

    RDouble *xyz = 0;
    if (sliceAxis == X_DIR)
    {
        xyz = x;
    }
    else if (sliceAxis == Y_DIR)
    {
        xyz = y;
    }
    else
    {
        xyz = z;
    }

    RDouble Small = 1.0e-10;
    int nCount = 0;

    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int nBCRegion = unstructBCSet->GetnBCRegion();

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
        int bcType = bcRegion->GetBCType();

        vector<int> *faceIndex = bcRegion->GetFaceIndex();
        for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
        {
            //! iFace is the face number in the set of faceIndex.
            int iFace = *iter;
            if (bcType != PHENGLEI::SOLID_SURFACE)
            {
                nCount += node_number_of_each_face[iFace];
                continue;
            }

            RDouble coor_min =   LARGE;
            RDouble coor_max = - LARGE;

            for (int iNode = 0; iNode < node_number_of_each_face[ iFace ]; ++ iNode)
            {
                int index = face2node[ nCount + iNode ];
                coor_min = MIN(xyz[ index ],coor_min);
                coor_max = MAX(xyz[ index ],coor_max);
            }

            if ((coor_min - slicePostion) * (coor_max - slicePostion) > 0.0)
            {
                nCount += node_number_of_each_face[ iFace ];
                continue;
            }

            for (int iNode = 0; iNode < node_number_of_each_face[ iFace ]; ++ iNode)
            {
                int Node_index0 = face2node[ nCount + iNode % node_number_of_each_face[ iFace ] ];
                int Node_index1 = face2node[ nCount + (iNode + 1) % node_number_of_each_face[ iFace ] ];

                if ((xyz[ Node_index0 ] - slicePostion) * (xyz[ Node_index1 ] - slicePostion) > 0.0) continue;

                if (fabs(xyz[ Node_index1 ] - xyz[ Node_index0 ]) < Small)
                {
                    xSlice.push_back(x[ Node_index0 ]);
                    xSlice.push_back(x[ Node_index1 ]);

                    ySlice.push_back(x[ Node_index0 ]);
                    ySlice.push_back(x[ Node_index1 ]);

                    zSlice.push_back(x[ Node_index0 ]);
                    zSlice.push_back(x[ Node_index1 ]);

                    continue;
                }

                if (sliceAxis == X_DIR)
                {
                    xSlice.push_back(slicePostion);
                }
                else if (fabs(x[ Node_index0 ] - x[ Node_index1 ]) < Small)
                {
                    xSlice.push_back(x[ Node_index0]);
                }
                else
                {
                    RDouble Node_x = (slicePostion - xyz[ Node_index0]) * (x[ Node_index1 ] - x[ Node_index0 ]) / (xyz[ Node_index1 ] - xyz[ Node_index0 ]) + x[ Node_index0 ];
                    xSlice.push_back(Node_x);
                }

                if (sliceAxis == Y_DIR)
                {
                    ySlice.push_back(slicePostion);
                }
                else if (fabs(y[ Node_index0 ] - y[ Node_index1 ]) < Small)
                {
                    ySlice.push_back(y[ Node_index0 ]);
                }
                else
                {
                    RDouble Node_y = (slicePostion - xyz[ Node_index0 ]) * (y[ Node_index1 ] - y[ Node_index0 ]) / (xyz[ Node_index1 ] - xyz[ Node_index0 ]) + y[ Node_index0 ];
                    ySlice.push_back(Node_y);
                }

                if (sliceAxis == Z_DIR)
                {
                    zSlice.push_back(slicePostion);
                }
                else if (fabs(z[ Node_index0 ] - z[ Node_index1 ]) < Small)
                {
                    zSlice.push_back(z[ Node_index0 ]);
                }
                else
                {
                    RDouble Node_z = (slicePostion - xyz[ Node_index0 ]) * (z[ Node_index1 ] - z[ Node_index0 ]) / (xyz[ Node_index1 ] - xyz[ Node_index0 ]) + z[ Node_index0 ];
                    zSlice.push_back(Node_z);
                }

                for (int iVarplot = 0; iVarplot < nvarplot; ++ iVarplot)
                {
                    if (fabs(qn[ iVarplot ][ Node_index0 ] - qn[ iVarplot ][ Node_index1 ]) < Small)
                    {
                        RDouble Node_q = qn[ iVarplot ][ Node_index0 ];
                        qnSlice[ iVarplot ].push_back(Node_q);
                    }
                    else
                    {
                        RDouble Node_q = (slicePostion - xyz[ Node_index0 ]) * (qn[ iVarplot ][ Node_index1 ] - qn[ iVarplot ][ Node_index0 ]) / (xyz[ Node_index1 ] - xyz[ Node_index0 ]) + qn[ iVarplot ][ Node_index0 ];
                        qnSlice[ iVarplot ].push_back(Node_q);
                    }
                }
            }

            nCount += node_number_of_each_face[ iFace ];
        }
    }

    size_t NumPts = xSlice.size();
    PHWrite(cdata, &NumPts, 1);

    if (NumPts == 0)  return;

    PHWrite(cdata, &xSlice[ 0 ], NumPts);
    PHWrite(cdata, &ySlice[ 0 ], NumPts);
    PHWrite(cdata, &zSlice[ 0 ], NumPts);
    for (int iVarplot = 0; iVarplot < nvarplot; ++ iVarplot)
    {
        PHWrite(cdata, &qnSlice[ iVarplot ][ 0 ], NumPts);
    }
}

void NSSolverUnstruct::BoundaryVisualization(Grid *gridIn, ActionKey *actkey, RDouble **qn, int nvarplot)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    DataContainer * cdata = actkey->GetTecData();

    int nBoundFace = grid->GetNBoundFace();

    RDouble * x = grid->GetX();
    RDouble * y = grid->GetY();
    RDouble * z = grid->GetZ();

    using namespace PHENGLEI;
    vector < vector < int > > face2nodelist(nBoundFace);
    vector < int > linkmap;
    set<int> bclist;

    set< pair<int, string> > bcnameMap;

    UnstructBCSet **bcr = grid->GetBCRecord();
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int bcType = bcr[iFace]->GetKey();
        bclist.insert(bcType);
        string bcname = bcr[iFace]->GetBCName();

        bcnameMap.insert(pair<int,string>(bcType, bcname));
    }
    
    int TecioMission;
    set< pair<int, string> >::iterator iter;
    for (iter = bcnameMap.begin(); iter != bcnameMap.end(); ++ iter)
    {
        if ((*iter).first == PHENGLEI::INTERFACE || (*iter).first < 0 ) continue;

        face2nodelist.resize(0);
        linkmap.resize(0);

        GetFace2NodeList(grid, *iter, linkmap, face2nodelist);

        if(linkmap.size() == 0)
        {
            continue;
        }

        TecioMission = TEC_SPACE::WriteBoundary;
        PHWrite(cdata, &TecioMission, 1);

        int BCtype = (*iter).first;
        PHWrite(cdata, &BCtype, 1);

        size_t NumPts = linkmap.size();
        size_t NumElements = face2nodelist.size();
        size_t NumFaces = NumElements * 4;
        PHWrite(cdata, &NumPts, 1);
        PHWrite(cdata, &NumElements, 1);
        PHWrite(cdata, &NumFaces, 1);

        size_t TotalNumFaceNodes_Rect = NumElements * 4;
        int *cell2node = new int [TotalNumFaceNodes_Rect];
        int *node_number_of_each_cell = new int [NumElements];

        int count = 0;
        for (std::size_t i = 0; i < face2nodelist.size(); ++ i)
        {
            node_number_of_each_cell[i] = 4;
            int np     = static_cast<int>(face2nodelist[i].size());
            int index0 = face2nodelist[i][0];
            for (int ip = 0; ip < np; ++ ip)
            {
                cell2node[count] = face2nodelist[i][ip];
                count ++;
            }

            for (int ip = np; ip < 4; ++ ip)
            {
                cell2node[count] = index0;
                count ++;
            }
        }

        PHWrite(cdata, &TotalNumFaceNodes_Rect, 1);
        PHWrite(cdata, node_number_of_each_cell, NumElements);

        PHWrite(cdata, &NumPts, 1);

        RDouble * xx = new RDouble [NumPts];
        RDouble * yy = new RDouble [NumPts];
        RDouble * zz = new RDouble [NumPts];

        for (std::size_t i = 0; i < linkmap.size(); ++ i)
        {
            int it = linkmap[i];
            xx[i] = x[it];
            yy[i] = y[it];
            zz[i] = z[it];
        }

        PHWrite(cdata, xx, NumPts);
        PHWrite(cdata, yy, NumPts);
        PHWrite(cdata, zz, NumPts);


        for (int m = 0; m < nvarplot; ++ m)
        {

            RDouble * qq = new RDouble [NumPts];

            for (std::size_t i = 0; i < linkmap.size(); ++ i)
            {
                int it = linkmap[i];
                qq[i] = qn[m][it];
            }

            PHWrite(cdata, qq, NumPts);
            delete [] qq;
        }

        PHWrite(cdata, cell2node, TotalNumFaceNodes_Rect);

        delete [] cell2node;
        delete [] node_number_of_each_cell;
        delete [] xx;
        delete [] yy;
        delete [] zz;
    }
}

void NSSolverUnstruct::VisualizationAverageFlow(ActionKey *actkey)
{
    int level = actkey->level;
    UnstructGrid * grid = UnstructGridCast(GetGrid(level));
    int nTotalNode = grid->GetNTotalNode();

    //! Average flow of DES statistics.
    int nvarplot = 6;

    RDouble ** qn = NewPointer2<RDouble>(nvarplot, nTotalNode);
    ComputeAverageVisualNodeField(grid, qn);

    int gridID = grid->GetGridID()->GetIndex();
    PHWrite(actkey->GetTecData(), gridID);

    int dimension = grid->GetDim();
    PHWrite(actkey->GetTecData(), dimension);

    int type = grid->Type();
    PHWrite(actkey->GetTecData(), &type, 1);

    BoundaryVisualization(grid, actkey, qn, nvarplot);

    if (WantVisualField(grid))
    {
        SaveDataForTecio(grid, actkey, qn, nvarplot);
    }

    int TecioMission = TEC_SPACE::Writecomplete;
    PHWrite(actkey->GetTecData(), &TecioMission, 1);

    int visualSlice = GlobalDataBase::GetIntParaFromDB("visualSlice");
    if (grid->GetDim() == THREE_D && visualSlice == 1)
    {
        SliceVisualization(grid, actkey, qn, nvarplot);
    }

    DelPointer2(qn);
}

Limiter *NSSolverUnstruct::CreateLimiter(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    Param_NSSolverUnstruct *parameters = GetControlParameters();

    int limitVariables = parameters->GetLimitVariables();
    int limitVector = parameters->GetLimitVector();
    int limiterType = parameters->GetLimiterType();

    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    int nVariables = 1;
    if (limitVector != 0 || limitVariables != 0)
    {
        nVariables = nEquation;
    }

    RDouble **q = reinterpret_cast<RDouble **>(grid->GetDataPtr("q"));

    Limiter *limiter = new Limiter(grid, q, nVariables, "limit");
    limiter->SetLimiterTypeID(limiterType);
    limiter->SetLimiterModel(limitVariables);
    limiter->SetLimiterVector(limitVector);
    limiter->InitLimitData();
    if (limiterType == ILMT_VENCAT)
    {
        int ivencat = GlobalDataBase::GetIntParaFromDB("ivencat");
        limiter->SetVencatLimiterType(ivencat);
    }

    limiter->Calculation(gridIn);

    return limiter;
}

FaceProxy * NSSolverUnstruct::CreateFaceProxy(Grid *gridIn)
{
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    FaceProxy *faceProxy = new FaceProxy();
    faceProxy->Create(SEGCTION_LENGTH, nEquation);
    return faceProxy;
}

GeomProxy * NSSolverUnstruct::CreateGeomProxy(Grid *gridIn)
{
    GeomProxy *geomProxy = new GeomProxy();
    geomProxy->Create(SEGCTION_LENGTH);
    return geomProxy;
}

NSFaceValue *NSSolverUnstruct::CreateNSFaceValue(Grid *gridIn)
{
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    int numberOfSpecies = parameters->GetNumberOfSpecies();

    NSFaceValue *facevar = new NSFaceValue(nEquation, numberOfSpecies, SEGCTION_LENGTH);
    return facevar;
}

RDouble NSSolverUnstruct::P0Wall(RDouble p1, RDouble gama, RDouble mach)
{
    RDouble p2, m2;
    if (mach > one)
    {
        p2 = p1 * ShockRelations::P2OverP1(gama, mach);
        m2 = ShockRelations::MachBehindShock(gama, mach);
        return p2  * IsentropicRelations::P0OverP(gama, m2);
    }
    else
    {
        return p1 * IsentropicRelations::P0OverP(gama, mach);
    }
}

RDouble NSSolverUnstruct::R0Wall(RDouble r1, RDouble gama, RDouble mach)
{
    RDouble r2, m2;
    if (mach > one)
    {
        r2 = r1 * ShockRelations::Rho2OverRho1(gama, mach);
        m2 = ShockRelations::MachBehindShock(gama, mach);
        return r2  * IsentropicRelations::Rho0OverRho(gama, m2);
    }
    else
    {
        return r1 * IsentropicRelations::Rho0OverRho(gama, mach);
    }
}

void NSSolverUnstruct::ComputeGradient(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int isViscous = parameters->IsViscous();

    if(grid->IsFinestGrid() || isViscous)
    {
        GetGradientField(gridIn);
    }
    else
    {
        return;
    }
}

void NSSolverUnstruct::GetGradientField(Grid *gridIn)
{
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    UnstructGrid *gridUnstruct = UnstructGridCast(gridIn);
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    RDouble **gradPrimtiveVarX = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("gradPrimtiveVarX"));
    RDouble **gradPrimtiveVarY = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("gradPrimtiveVarY"));
    RDouble **gradPrimtiveVarZ = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("gradPrimtiveVarZ"));
    FieldProxy *qProxy = GetFieldProxy(gridIn, "q");
    RDouble **q     = qProxy->GetField_UNS();
    RDouble **qnode = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("qnode"));
    string gradientName = parameters->GetGradientName();

    for (int m = 0; m < nEquation; ++ m)
    {
        if (gradientName == "ggnode" || gradientName == "ggnodelaplacian")
        {
            gridUnstruct->CompGradientGGNode_NEW(q[m], qnode[m], gradPrimtiveVarX[m], gradPrimtiveVarY[m], gradPrimtiveVarZ[m]);
        }
        else if (gradientName == "ggcell")
        {
            gridUnstruct->CompGradientGGCell(q[m], gradPrimtiveVarX[m], gradPrimtiveVarY[m], gradPrimtiveVarZ[m]);
        }
        else if (gradientName == "lsq")
        {
            gridUnstruct->CompGradientLSQ(q[m], gradPrimtiveVarX[m], gradPrimtiveVarY[m], gradPrimtiveVarZ[m]);
        }    
        else if (gradientName == "ggnode_weight")
        {
            gridUnstruct->CompGradientGGNodeWeight(q[m], gradPrimtiveVarX[m], gradPrimtiveVarY[m], gradPrimtiveVarZ[m]);
        }
        else if (gradientName == "gg_m2")
        {
            gridUnstruct->CompGradientGGModified2(q[m], gradPrimtiveVarX[m], gradPrimtiveVarY[m], gradPrimtiveVarZ[m]);
        }
        else if (gradientName == "ggcellnew")
        {
            gridUnstruct->CompGradientGGCellNew(q[m], gradPrimtiveVarX[m], gradPrimtiveVarY[m], gradPrimtiveVarZ[m]);
        }
        else if (gradientName == "ggcellw")
        {
            gridUnstruct->CompGradientGGCellW(q[m], gradPrimtiveVarX[m], gradPrimtiveVarY[m], gradPrimtiveVarZ[m]);
        }
        else
        {
            TK_Exit::ExceptionExit("No reconstruction method has been choosed ! /n");
        }
    }

    int isFVMOrFDM = GlobalDataBase::GetIntParaFromDB("isFVMOrFDM");
    if (systemgridtype == MIXGRID && isFVMOrFDM == FV_METHOD)
    {
        //! For mixgrid, gradient is communicated using another function.
    }
    else
    {
        PHSPACE::CommunicateInterfaceValue(gridUnstruct, gradPrimtiveVarX, "dqdx", nEquation);
        PHSPACE::CommunicateInterfaceValue(gridUnstruct, gradPrimtiveVarY, "dqdy", nEquation);
        PHSPACE::CommunicateInterfaceValue(gridUnstruct, gradPrimtiveVarZ, "dqdz", nEquation);
    }

    RDouble **gradTemperatureX = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("gradTemperatureX"));
    RDouble **gradTemperatureY = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("gradTemperatureY"));
    RDouble **gradTemperatureZ = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("gradTemperatureZ"));
    RDouble **t     = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("t"));
    RDouble **tnode = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("tnode"));

    for (int m = 0; m < nTemperatureModel; ++ m)
    {
        if (gradientName == "ggnode" || gradientName == "ggnodelaplacian")
        {
            gridUnstruct->CompGradientGGNode_NEW(t[m], tnode[m], gradTemperatureX[m], gradTemperatureY[m], gradTemperatureZ[m]);
        }
        else if (gradientName == "ggcell")
        {
            gridUnstruct->CompGradientGGCell(t[m], gradTemperatureX[m], gradTemperatureY[m], gradTemperatureZ[m]);
        }
        else if (gradientName == "lsq")
        {
            gridUnstruct->CompGradientLSQ(t[m], gradTemperatureX[m], gradTemperatureY[m], gradTemperatureZ[m]);
        }    
        else if (gradientName == "ggnode_weight")
        {
            gridUnstruct->CompGradientGGNodeWeight(t[m], gradTemperatureX[m], gradTemperatureY[m], gradTemperatureZ[m]);
        }
        else if (gradientName == "gg_m2")
        {
            gridUnstruct->CompGradientGGModified2(t[m], gradTemperatureX[m], gradTemperatureY[m], gradTemperatureZ[m]);
        }
        else if (gradientName == "ggcellnew")
        {
            gridUnstruct->CompGradientGGCellNew(t[m], gradTemperatureX[m], gradTemperatureY[m], gradTemperatureZ[m]);
        }
        else if (gradientName == "ggcellw")
        {
            gridUnstruct->CompGradientGGCellW(t[m], gradTemperatureX[m], gradTemperatureY[m], gradTemperatureZ[m]);
        }
        else
        {
            TK_Exit::ExceptionExit("No reconstruction method has been choosed ! /n");
        }
    }

    if (systemgridtype == MIXGRID && isFVMOrFDM == FV_METHOD)
    {
        //! For mixgrid, gradient is communicated using another function.
    }
    else
    {
        PHSPACE::CommunicateInterfaceValue(gridUnstruct, gradTemperatureX, "dtdx", nTemperatureModel);
        PHSPACE::CommunicateInterfaceValue(gridUnstruct, gradTemperatureY, "dtdy", nTemperatureModel);
        PHSPACE::CommunicateInterfaceValue(gridUnstruct, gradTemperatureZ, "dtdz", nTemperatureModel);
    }
    delete qProxy;
}

void NSSolverUnstruct::InitiallimiterOnlyForMixGrid(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    RDouble **q = reinterpret_cast<RDouble **>(grid->GetDataPtr("q"));
    Param_NSSolverUnstruct *parameters = GetControlParameters();

    int limitVector = parameters->GetLimitVector();
    int limiterType = parameters->GetLimiterType();

    int nVariables = 1;
    if (limitVector == 1)
    {
        nVariables = 5;
    }

    limiterOnlyForMixGrid = new Limiter(grid, q, nVariables, "limit");

    limiterOnlyForMixGrid->InitLimitData();

    if (limiterType == ILMT_VENCAT)
    {
        int ivencat = GlobalDataBase::GetIntParaFromDB("ivencat");
        limiterOnlyForMixGrid->SetVencatLimiterType(ivencat);
    }
}

void NSSolverUnstruct::ComputelimiterOnlyForMixGrid(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    RDouble **q = reinterpret_cast<RDouble **>(grid->GetDataPtr("q"));
    Param_NSSolverUnstruct *parameters = GetControlParameters();

    int limitVariables = parameters->GetLimitVariables();
    int limitVector = parameters->GetLimitVector();
    int limiterType = parameters->GetLimiterType();
    
    if (limiterType >= ILMT_STRUCT || limiterType == ILMT_NOLIM || limiterType == ILMT_FIRST)
    {
        return;
    }

    RDouble **dqdx = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarX"));
    RDouble **dqdy = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarY"));
    RDouble **dqdz = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarZ"));

    int nEquation = 5;

    const int SCALAR_LIMITER    = 0;
    const int VECTOR_LIMITER    = 1;
    
    const int LIMIT_BY_RHO_P    = 0;
    const int LIMIT_BY_EACH_VAR = 1;

    RDouble *limit_tmp = 0;
    if (limiterType == ILMT_BARTH)
    {
        if(limitVector == SCALAR_LIMITER)
        {
            //! Each variables has the same limiter value.
            limit_tmp = limiterOnlyForMixGrid->GetLimiter(0);
            if (limitVariables == LIMIT_BY_RHO_P)
            {
                limiterOnlyForMixGrid->BarthLimiter(limit_tmp, q[IR], dqdx[IR], dqdy[IR], dqdz[IR]);
                limiterOnlyForMixGrid->BarthLimiter(limit_tmp, q[IP], dqdx[IP], dqdy[IP], dqdz[IP]);
            }
            else if (limitVariables == LIMIT_BY_EACH_VAR)
            {
                for (int m = 0; m < nEquation; ++ m)
                {
                    limiterOnlyForMixGrid->BarthLimiter(limit_tmp, q[m], dqdx[m], dqdy[m], dqdz[m]);
                }
            }
        }
        else if(limitVector == VECTOR_LIMITER)
        {
            //! Limiter value of each variable is computed by itself.
            for (int m = 0; m < nEquation; ++ m)
            {
                limit_tmp = limiterOnlyForMixGrid->GetLimiter(m);
                limiterOnlyForMixGrid->BarthLimiter(limit_tmp, q[m], dqdx[m], dqdy[m], dqdz[m]);
            }
        }
    }
    else if (limiterType == ILMT_VENCAT)
    {
        if (limitVector == SCALAR_LIMITER)
        {
            //! Each variables has the same limiter value.
            limit_tmp = limiterOnlyForMixGrid->GetLimiter(0);
            if (limitVariables == LIMIT_BY_RHO_P)
            {
                limiterOnlyForMixGrid->VencatLimiter(limit_tmp, q[IR], dqdx[IR], dqdy[IR], dqdz[IR], IR);
                limiterOnlyForMixGrid->VencatLimiter(limit_tmp, q[IP], dqdx[IP], dqdy[IP], dqdz[IP], IP);
            }
            else if (limitVariables == LIMIT_BY_EACH_VAR)
            {
                for (int m = 0; m < nEquation; ++ m)
                {
                    limiterOnlyForMixGrid->VencatLimiter(limit_tmp, q[m], dqdx[m], dqdy[m], dqdz[m], m);
                }
            }
        }
        else if(limitVector == VECTOR_LIMITER)
        {
            //! Limiter value of each variable is computed by itself.
            for (int m = 0; m < nEquation; ++ m)
            {
                limit_tmp = limiterOnlyForMixGrid->GetLimiter(m);
                limiterOnlyForMixGrid->VencatLimiter(limit_tmp, q[m], dqdx[m], dqdy[m], dqdz[m], m);
            }
        }
    }
    else
    {
        TK_Exit::ExceptionExit("Error: this unstructured limiter is not exist !\n", true);
    }
    
    
    if (limitVector == SCALAR_LIMITER)
    {
        limit_tmp = limiterOnlyForMixGrid->GetLimiter(0);
        SetLimitBoundary(grid, limit_tmp);
    }
    else if(limitVector == VECTOR_LIMITER)
    {
        for (int m = 0; m < nEquation; ++ m)
        {
            limit_tmp = limiterOnlyForMixGrid->GetLimiter(m);
            SetLimitBoundary(grid, limit_tmp);
        }
    }
}



LIB_EXPORT void NSSolverUnstruct::ComputePostVisualVariables(Post_Visual * postVisualization)
{
    UnstructGrid *grid = UnstructGridCast(postVisualization->GetGrid());

    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    RDouble **t = reinterpret_cast< RDouble ** > (grid->GetDataPtr("t"));

    if (postVisualization->GetFlowType() == AverageFlow)
    {
        RDouble **qAverage = reinterpret_cast <RDouble **> (grid->GetDataPtr("qAverage"));

        string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_DENSITY);
        postVisualization->UpdateVisualNodeVarPtr(varName, qAverage[IDX::IR]);

        varName = postVisualization->GetVariableName(PHSPACE::VISUAL_U);
        postVisualization->UpdateVisualNodeVarPtr(varName, qAverage[IDX::IU]);

        varName = postVisualization->GetVariableName(PHSPACE::VISUAL_V);
        postVisualization->UpdateVisualNodeVarPtr(varName, qAverage[IDX::IV]);

        varName = postVisualization->GetVariableName(PHSPACE::VISUAL_W);
        postVisualization->UpdateVisualNodeVarPtr(varName, qAverage[IDX::IW]);

        varName = postVisualization->GetVariableName(PHSPACE::VISUAL_PRESSURE);
        postVisualization->UpdateVisualNodeVarPtr(varName, qAverage[IDX::IP]);

        return;
    }

    //! Added by zzp 202108, for ReynoldsStress output
    if (postVisualization->GetFlowType() == AverageReynoldsStress)
    {
        RDouble **tauAverage = reinterpret_cast <RDouble **> (grid->GetDataPtr("tauAverage"));

        string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_TAU_XX);
        postVisualization->UpdateVisualNodeVarPtr(varName, tauAverage[0]);

        varName = postVisualization->GetVariableName(PHSPACE::VISUAL_TAU_YY);
        postVisualization->UpdateVisualNodeVarPtr(varName, tauAverage[1]);

        varName = postVisualization->GetVariableName(PHSPACE::VISUAL_TAU_ZZ);
        postVisualization->UpdateVisualNodeVarPtr(varName, tauAverage[2]);

        varName = postVisualization->GetVariableName(PHSPACE::VISUAL_TAU_XY);
        postVisualization->UpdateVisualNodeVarPtr(varName, tauAverage[3]);

        varName = postVisualization->GetVariableName(PHSPACE::VISUAL_TAU_XZ);
        postVisualization->UpdateVisualNodeVarPtr(varName, tauAverage[4]);

        varName = postVisualization->GetVariableName(PHSPACE::VISUAL_TAU_YZ);
        postVisualization->UpdateVisualNodeVarPtr(varName, tauAverage[5]);

        return;
    }

    using namespace IDX;
    if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_DENSITY))
    {
        string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_DENSITY);
        postVisualization->UpdateVisualNodeVarPtr(varName, q[IDX::IR]);
    }
    if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_U))
    {
        string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_U);
        postVisualization->UpdateVisualNodeVarPtr(varName, q[IDX::IU]);
    }
    if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_V))
    {
        string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_V);
        postVisualization->UpdateVisualNodeVarPtr(varName, q[IDX::IV]);
    }
    if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_W))
    {
        string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_W);
        postVisualization->UpdateVisualNodeVarPtr(varName, q[IDX::IW]);
    }
    if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_PRESSURE))
    {
        string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_PRESSURE);
        postVisualization->UpdateVisualNodeVarPtr(varName, q[IDX::IP]);
    }

    if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_TEMPERATURE))
    {
        string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_TEMPERATURE);
        postVisualization->UpdateVisualNodeVarPtr(varName, t[IDX::ITT]);
    }    

    if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_MACH))
    {
        string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_MACH);
        RDouble * mach  = CompMachNumber(grid);
        postVisualization->UpdateVisualNodeVarPtr(varName, mach);
        delete [] mach;
    }

    if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_STREAMLINE_U))
    {
        string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_STREAMLINE_U);
        postVisualization->UpdateVisualNodeVarPtr(varName, q[IDX::IU]);
    }

    if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_STREAMLINE_V))
    {
        string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_STREAMLINE_V);
        postVisualization->UpdateVisualNodeVarPtr(varName, q[IDX::IV]);
    }

    if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_STREAMLINE_W))
    {
        string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_STREAMLINE_W);
        postVisualization->UpdateVisualNodeVarPtr(varName, q[IDX::IW]);
    }

    if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_TIME_STEP))
    {
        string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_TIME_STEP);
        RDouble * dt  = reinterpret_cast< RDouble * > (grid->GetDataPtr("dt" ));
        postVisualization->UpdateVisualNodeVarPtr(varName, dt);
    }
    if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_VOLUME))
    {
        string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_VOLUME);
        RDouble *vol  = grid->GetCellVolume();
        postVisualization->UpdateVisualNodeVarPtr(varName, vol);
    }
    if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_WALL_DIST))
    {
        string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_WALL_DIST);
        RDouble* wallDistance = grid->GetWallDist();
        postVisualization->UpdateVisualNodeVarPtr(varName, wallDistance);
    }
    if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_IBLANK))
    {
        string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_IBLANK);
        int *iBlank = reinterpret_cast<int *> (grid->GetBlankIndex());
        postVisualization->UpdateVisualNodeVarPtr(varName, iBlank);
    }

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int viscousType = parameters->GetViscousType();
    string viscousName = parameters->GetViscousName();
    int iLES = GlobalDataBase::GetIntParaFromDB("iLES");

    if (viscousType == LAMINAR && iLES == NOLES_SOLVER)
    {
        if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_VISCOSITY_LAMINAR))
        {
            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_VISCOSITY_LAMINAR);
            RDouble *visl = reinterpret_cast< RDouble * > (grid->GetDataPtr("visl"));
            postVisualization->UpdateVisualNodeVarPtr(varName, visl);
        }
    }
    else if (viscousType > LAMINAR || iLES == LES_SOLVER)
    {
        if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_VISCOSITY_LAMINAR))
        {
            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_VISCOSITY_LAMINAR);
            RDouble *visl = reinterpret_cast< RDouble * > (grid->GetDataPtr("visl"));
            postVisualization->UpdateVisualNodeVarPtr(varName, visl);
        }

        if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_VISCOSITY_TURBULENT))
        {
            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_VISCOSITY_TURBULENT);
            RDouble *vist = reinterpret_cast< RDouble * > (grid->GetDataPtr("vist"));
            postVisualization->UpdateVisualNodeVarPtr(varName, vist);
        }

        if (viscousName.substr(0,6) == "2eq-kw")
        {
            RDouble ** q_turb    = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q_turb" ));
            if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_MODELED_TKE))
            {
                string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_MODELED_TKE);
                postVisualization->UpdateVisualNodeVarPtr(varName, q_turb[IDX::IKE]);
            }
            if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_MODELED_DISSIPATION))
            {
                string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_MODELED_DISSIPATION);
                postVisualization->UpdateVisualNodeVarPtr(varName, q_turb[IDX::IKW]);
            }
            if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_INTERMITTENCY))
            {
                string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_INTERMITTENCY);
                postVisualization->UpdateVisualNodeVarPtr(varName, q_turb[IDX::IGAMA]);
            }
            if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_MOMENTUMTHICKREYNOLDS))
            {
                string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_MOMENTUMTHICKREYNOLDS);
                postVisualization->UpdateVisualNodeVarPtr(varName, q_turb[IDX::IRECT]);
            }
            if (viscousName.substr(0,17) == "2eq-kw-menter-sst")   //! For blending function F1 and F2
            {
                if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_SST_F1))
                {
                    RDouble * blend = reinterpret_cast< RDouble * > (grid->GetDataPtr("blend" ));   //! For blengding function F1
                    string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_SST_F1);
                    postVisualization->UpdateVisualNodeVarPtr(varName, blend);
                }
                if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_SST_F2))
                {
                    RDouble * SST_F2 = reinterpret_cast< RDouble * > (grid->GetDataPtr("SST_F2" ));   //! For blengding function F2
                    string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_SST_F2);
                    postVisualization->UpdateVisualNodeVarPtr(varName, SST_F2);
                }
            }
        }
    }

    int visualField = parameters->GetPlotFieldType();
    if(WantVisualField(grid) || visualField == TEC_SPACE::BlockVisual)
    {
        RDouble * vorticity_x, * vorticity_y, * vorticity_z, * vorticityMagnitude, * strain_rate, * Q_criteria;

        int nTotalCell = grid->GetNTotalCell();
        int nBoundFace = grid->GetNBoundFace();
        int nTotal = nTotalCell + nBoundFace;

        vorticity_x = new RDouble [nTotal];
        vorticity_y = new RDouble [nTotal];
        vorticity_z = new RDouble [nTotal];
        vorticityMagnitude = new RDouble [nTotal];
        strain_rate = new RDouble [nTotal];
        Q_criteria = new RDouble [nTotal];
        ComputeVorticitybyQCriteria(grid, vorticity_x, vorticity_y, vorticity_z, vorticityMagnitude, strain_rate, Q_criteria);

        if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_VORTICITY_X))
        {
            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_VORTICITY_X);
            postVisualization->UpdateVisualNodeVarPtr(varName, vorticity_x);
        }

        if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_VORTICITY_Y))
        {
            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_VORTICITY_Y);
            postVisualization->UpdateVisualNodeVarPtr(varName, vorticity_y);
        }

        if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_VORTICITY_Z))
        {
            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_VORTICITY_Z);
            postVisualization->UpdateVisualNodeVarPtr(varName, vorticity_z);
        }

        if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_VORTICITY_MAGNITUDE))
        {
            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_VORTICITY_MAGNITUDE);
            postVisualization->UpdateVisualNodeVarPtr(varName, vorticityMagnitude);
        }

        if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_STRAIN_RATE))
        {
            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_STRAIN_RATE);
            postVisualization->UpdateVisualNodeVarPtr(varName, strain_rate);
        }

        if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_Q_CRITERIA))
        {
            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_Q_CRITERIA);
            postVisualization->UpdateVisualNodeVarPtr(varName, Q_criteria);
        }

        delete []vorticity_x;
        delete []vorticity_y;
        delete []vorticity_z;
        delete []vorticityMagnitude;
        delete []strain_rate;
        delete []Q_criteria;
    }
    else if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_VORTICITY_X) 
         || postVisualization->IsNeedVisualization(PHSPACE::VISUAL_VORTICITY_Y)
         || postVisualization->IsNeedVisualization(PHSPACE::VISUAL_VORTICITY_Z)
         || postVisualization->IsNeedVisualization(PHSPACE::VISUAL_VORTICITY_MAGNITUDE)
         || postVisualization->IsNeedVisualization(PHSPACE::VISUAL_STRAIN_RATE)
         || postVisualization->IsNeedVisualization(PHSPACE::VISUAL_Q_CRITERIA))
    {
        TK_Exit::ExceptionExit("Error: can not visualization vorticity when plotFieldType = 0!");
    }

    if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_CP))
    {
        string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_CP);
        RDouble * cp = ComputeCp(grid);
        postVisualization->UpdateVisualNodeVarPtr(varName, cp);
        delete [] cp;
    }
    if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_GRADIENT_UX))
    {
        RDouble **dqdx = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarX"));        
        RDouble *dudx = dqdx[IU];
        postVisualization->UpdateVisualNodeVarPtr("gradientUx", dudx);       
    }
    if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_GRADIENT_UY))
    {
        RDouble **dqdy = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarY"));        
        RDouble *dudy = dqdy[IU];
        postVisualization->UpdateVisualNodeVarPtr("gradientUy", dudy);       
    }
    if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_GRADIENT_VX))
    {
        RDouble **dqdx = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarX"));        
        RDouble *dvdx = dqdx[IV];
        postVisualization->UpdateVisualNodeVarPtr("gradientVx", dvdx);       
    }
    if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_GRADIENT_VY))
    {
        RDouble **dqdy = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarY"));        
        RDouble *dvdy = dqdy[IV];
        postVisualization->UpdateVisualNodeVarPtr("gradientVy", dvdy);       
    }

    if (postVisualization->IsNeedVisualization(PHSPACE::VISUAL_DIMENSIONAL_DENSITY)
     || postVisualization->IsNeedVisualization(PHSPACE::VISUAL_DIMENSIONAL_U)
     || postVisualization->IsNeedVisualization(PHSPACE::VISUAL_DIMENSIONAL_V)
     || postVisualization->IsNeedVisualization(PHSPACE::VISUAL_DIMENSIONAL_W)
     || postVisualization->IsNeedVisualization(PHSPACE::VISUAL_DIMENSIONAL_PRESSURE)
     || postVisualization->IsNeedVisualization(PHSPACE::VISUAL_DIMENSIONAL_TEMPERATURE))
    {
        RDouble **dimensionalVariables = ComputeDimensionalVariables(grid);

        if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_DIMENSIONAL_DENSITY))
        {
            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_DIMENSIONAL_DENSITY);
            postVisualization->UpdateVisualNodeVarPtr(varName, dimensionalVariables[0]);
        }
        if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_DIMENSIONAL_U))
        {
            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_DIMENSIONAL_U);
            postVisualization->UpdateVisualNodeVarPtr(varName, dimensionalVariables[1]);
        }
        if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_DIMENSIONAL_V))
        {
            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_DIMENSIONAL_V);
            postVisualization->UpdateVisualNodeVarPtr(varName, dimensionalVariables[2]);
        }
        if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_DIMENSIONAL_W))
        {
            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_DIMENSIONAL_W);
            postVisualization->UpdateVisualNodeVarPtr(varName, dimensionalVariables[3]);
        }
        if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_DIMENSIONAL_PRESSURE))
        {
            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_DIMENSIONAL_PRESSURE);
            postVisualization->UpdateVisualNodeVarPtr(varName, dimensionalVariables[4]);
        }
        if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_DIMENSIONAL_TEMPERATURE))
        {
            string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_DIMENSIONAL_TEMPERATURE);
            postVisualization->UpdateVisualNodeVarPtr(varName, dimensionalVariables[5]);
        }

        DelPointer2(dimensionalVariables);
    }

    //! Chemical species must be the last.
    int nChemical = parameters->GetChemicalFlag();
    if(nChemical)
    {
        using namespace GAS_SPACE;

        int nm = parameters->GetNSEquationNumber();
        int nLaminar  = parameters->GetLaminarNumber();
        int nTemperatureModel = parameters->GetTemperatureModel();

        int nEquation = GetNumberOfEquations();

        string * varname = gas->GetNameOfSpecies();
        for (int m = nm; m < nEquation; ++ m)
        {
            string varName = varname[m-nm];
            postVisualization->UpdateVisualNodeVarPtr(varName, q[m]);
        }
    }
}

LIB_EXPORT void NSSolverUnstruct::ComputePostProbesVariables(Post_Probes * postProbesVar)
{
    UnstructGrid *grid = UnstructGridCast(postProbesVar->GetGrid());

    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    RDouble **t = reinterpret_cast< RDouble ** > (grid->GetDataPtr("t"));

    using namespace IDX;
    if(postProbesVar->IsNeedMonitoring(PHSPACE::PROBE_DENSITY))
    {
        string varName = postProbesVar->GetProbesVariableName(PHSPACE::PROBE_DENSITY);
        postProbesVar->UpdateProbesVarPtr(varName, q[IDX::IR]);
    }

    if(postProbesVar->IsNeedMonitoring(PHSPACE::PROBE_U))
    {
        string varName = postProbesVar->GetProbesVariableName(PHSPACE::PROBE_U);
        postProbesVar->UpdateProbesVarPtr(varName, q[IDX::IU]);
    }
    if(postProbesVar->IsNeedMonitoring(PHSPACE::PROBE_V))
    {
        string varName = postProbesVar->GetProbesVariableName(PHSPACE::PROBE_V);
        postProbesVar->UpdateProbesVarPtr(varName, q[IDX::IV]);
    }
    if(postProbesVar->IsNeedMonitoring(PHSPACE::PROBE_W))
    {
        string varName = postProbesVar->GetProbesVariableName(PHSPACE::PROBE_W);
        postProbesVar->UpdateProbesVarPtr(varName, q[IDX::IW]);
    }
    if(postProbesVar->IsNeedMonitoring(PHSPACE::PROBE_PRESSURE))
    {
        string varName = postProbesVar->GetProbesVariableName(PHSPACE::PROBE_PRESSURE);
        postProbesVar->UpdateProbesVarPtr(varName, q[IDX::IP]);
    }
    if(postProbesVar->IsNeedMonitoring(PHSPACE::PROBE_TEMPERATURE))
    {
        string varName = postProbesVar->GetProbesVariableName(PHSPACE::PROBE_TEMPERATURE);
        postProbesVar->UpdateProbesVarPtr(varName, t[IDX::ITT]);
    }
    if(postProbesVar->IsNeedMonitoring(PHSPACE::PROBE_MACH))
    {
        string varName = postProbesVar->GetProbesVariableName(PHSPACE::VISUAL_MACH);
        RDouble * mach  = CompMachNumber(grid);
        postProbesVar->UpdateProbesVarPtr(varName, mach);
        delete [] mach;
    }

}

LIB_EXPORT void NSSolverUnstruct::InitControlParameters()
{
    FreeControlParameters();
    controlParameters = new Param_NSSolverUnstruct();
    controlParameters->Init();
}

LIB_EXPORT Param_NSSolverUnstruct * NSSolverUnstruct::GetControlParameters() const
{
    return static_cast< Param_NSSolverUnstruct * > (controlParameters);
}

void NSSolverUnstruct::CompressDQ(DataContainer *&dataContainer, FieldProxy *fieldProxy, Grid *gridIn, const int &zoneIndex, const int &neighborZoneIndex, const int &nEquation)
{
    UnstructGrid * grid = UnstructGridCast(gridIn);

    InterfaceInfo * interfaceInformation = grid->GetInterfaceInfo();
    if (!interfaceInformation) return;

    int iNeighborZone                            = interfaceInformation->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = interfaceInformation->GetNIFaceOfNeighbor(iNeighborZone);
    int * interfaceIndexContainerForSend         = interfaceInformation->GetFaceIndexForSend(iNeighborZone);
    dataContainer->moveto_begin();

    RDouble **dq = fieldProxy->GetField_UNS();

    for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; -- iGhostLayer)
    {
    for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
    {
        int sourceCell;
        int iFace = interfaceIndexContainerForSend[ iLocalFace ];
            grid->GetSourceIndex(iFace, iGhostLayer + 1, sourceCell);
        for (int m = 0; m < nEquation; ++ m)
        {
            PHWrite(dataContainer, dq[m][sourceCell]);
        }
    }
    }
}

void NSSolverUnstruct::DecompressDQ(DataContainer *&dataContainer, FieldProxy *fieldProxy, Grid *gridIn, const int &zoneIndex, const int &neighborZoneIndex, const int &nEquation)
{
    UnstructGrid * grid = UnstructGridCast(gridIn);

    InterfaceInfo * interfaceInformation = grid->GetInterfaceInfo();
    if (!interfaceInformation) return;

    int iNeighborZone                            = interfaceInformation->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = interfaceInformation->GetNIFaceOfNeighbor(iNeighborZone);
    int * interfaceIndexContainerForReceive      = interfaceInformation->GetFaceIndexForRecv(iNeighborZone);

    dataContainer->moveto_begin();
    RDouble **dq = fieldProxy->GetField_UNS();

    for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; -- iGhostLayer)
    {
    for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
    {
        int iFace = interfaceIndexContainerForReceive[ iLocalFace ];
        int targetCell;
            grid->GetTargetIndex(iFace, iGhostLayer + 1, targetCell);

        for (int m = 0; m < nEquation; ++ m)
        {
            PHRead(dataContainer, dq[m][targetCell]);
        }
    }
    }
}

void NSSolverUnstruct::CompressSpecificArrayToInterface(DataContainer *&dataContainer, const string &fieldName, Grid *gridIn, const int &neighborZoneIndex, const int &nEquation)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
    if (!interfaceInformation)
    {
        return;
    }

    int iNeighborZone                            = interfaceInformation->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = interfaceInformation->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForSend          = interfaceInformation->GetFaceIndexForSend(iNeighborZone);
    dataContainer->moveto_begin();
    if (nEquation == 0)
    {
        RDouble *fieldSend = reinterpret_cast <RDouble *> (grid->GetDataPtr(fieldName));

        for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; -- iGhostLayer)
        {
            for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
            {
                int sourceCell;
                int iFace = interfaceIndexContainerForSend[iLocalFace];
                grid->GetSourceIndex(iFace, iGhostLayer + 1, sourceCell);
                PHWrite(dataContainer, fieldSend[sourceCell]);
            }
        }
    }
    else
    {
        RDouble **fieldSend = reinterpret_cast <RDouble **> (grid->GetDataPtr(fieldName));

        for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; -- iGhostLayer)
        {
            for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
            {
                int sourceCell;
                int iFace = interfaceIndexContainerForSend[iLocalFace];
                grid->GetSourceIndex(iFace, iGhostLayer + 1, sourceCell);
                for (int m = 0; m < nEquation; ++ m)
                {
                    PHWrite(dataContainer, fieldSend[m][sourceCell]);
                }
            }
        }
    }
}

void NSSolverUnstruct::DecompressArrayFromInterface(DataContainer *&dataContainer, const string &fieldName, Grid *gridIn, const int &neighborZoneIndex, const int &nEquation)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
    if (!interfaceInformation)
    {
        return;
    }

    int iNeighborZone                            = interfaceInformation->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = interfaceInformation->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForReceive       = interfaceInformation->GetFaceIndexForRecv(iNeighborZone);
    dataContainer->moveto_begin();
    if (nEquation == 0)
    {
        RDouble *fieldRecv = reinterpret_cast <RDouble *> (grid->GetDataPtr(fieldName));

        for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; -- iGhostLayer)
        {
            for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
            {
                int targetCell;
                int iFace = interfaceIndexContainerForReceive[iLocalFace];
                grid->GetTargetIndex(iFace, iGhostLayer + 1, targetCell);
                PHRead(dataContainer, fieldRecv[targetCell]);
            }
        }
    }
    else
    {
        RDouble **fieldRecv = reinterpret_cast <RDouble **> (grid->GetDataPtr(fieldName));

        for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; -- iGhostLayer)
        {
            for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
            {
                int iFace = interfaceIndexContainerForReceive[iLocalFace];
                int targetCell;
                grid->GetTargetIndex(iFace, iGhostLayer + 1, targetCell);
                for (int m = 0; m < nEquation; ++ m)
                {
                    PHRead(dataContainer, fieldRecv[m][targetCell]);
                }
            }
        }
    }
}

void NSSolverUnstruct::OutflowBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct)
{
    //! Get the basic variables for boundary calculating.
    UnstructGrid *grid   = UnstructGridCast(gridIn);
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();
    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));

    //! Get the number of equations.
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    int iFace, le, re;
    //! get the face indexes of this boundary condition region.
    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();

    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        //! iFace is the face number in the set of faceIndex.
        iFace = *iter;
        le = leftCellofFace[iFace];
        re = rightCellofFace[iFace];
        for (int m = 0; m < nEquation; ++ m)
        {
            q[m][re] = q[m][le];
        }
    }
}

void NSSolverUnstruct::SymmetryBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct)
{
    //! Get the basic variables for boundary calculating.
    UnstructGrid *grid   = UnstructGridCast(gridIn);
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();
    RDouble   *xfn       = grid->GetFaceNormalX();
    RDouble   *yfn       = grid->GetFaceNormalY();
    RDouble   *zfn       = grid->GetFaceNormalZ();
    RDouble   *vgn       = grid->GetFaceNormalVelocity();
    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));

    //! Get the number of equations.
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    int iFace, le, re;

    //! Get the face indexes of this boundary condition region.
    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();

    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        iFace = *iter;
        le = leftCellofFace[iFace];
        re = rightCellofFace[iFace];
        
        for (int m = 0; m < nEquation; ++ m)
        {
            q[m][re] = q[m][le];
        }

        RDouble uL = q[IU][le];
        RDouble vL = q[IV][le];
        RDouble wL = q[IW][le];
        RDouble vn = xfn[iFace] * uL + yfn[iFace] * vL + zfn[iFace] * wL - vgn[iFace];

        q[IU][re] = uL - two * xfn[iFace] * vn;
        q[IV][re] = vL - two * yfn[iFace] * vn;
        q[IW][re] = wL - two * zfn[iFace] * vn;
    }
}

void NSSolverUnstruct::ViscousIsotropicWallBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct)
{
    //! Get the basic variables for boundary calculating.
    UnstructGrid *grid   = UnstructGridCast(gridIn);
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();
    RDouble *faceVelocityX = grid->GetFaceVelocityX();
    RDouble *faceVelocityY = grid->GetFaceVelocityY();
    RDouble *faceVelocityZ = grid->GetFaceVelocityZ();
    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));

    //! Get the number of equations.
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    RDouble refDimensionalTemperature = parameters->GetRefDimensionalTemperature();

    //! Get the variables defined in boundarycondition.hypara file.
    Data_Param *bcParamDB = bcRegionUnstruct->GetBCParamDataBase();
    RDouble wallTemperature = 0.0;

    if (bcParamDB)
    {
        if (!bcParamDB->IsExist("wallTemperature", PHDOUBLE, 1))
        {
            wallTemperature = GlobalDataBase::GetDoubleParaFromDB("wallTemperature");
            bcParamDB->UpdateData("wallTemperature", &wallTemperature, PHDOUBLE, 1);
        }
        else
        {
            bcParamDB->GetData("wallTemperature", &wallTemperature, PHDOUBLE, 1);
        }
    }
    else
    {
        wallTemperature = GlobalDataBase::GetDoubleParaFromDB("wallTemperature");
    }
    wallTemperature /= refDimensionalTemperature;

    bool wallTempArrayExist = false;
    RDouble *wallTempArray;
    int countFace = 0;
    if (bcRegionUnstruct->CheckFieldData("wallTempArray"))
    {
        wallTempArrayExist = true;
        wallTempArray = reinterpret_cast <RDouble *> (bcRegionUnstruct->GetFieldDataPtr("wallTempArray"));
    }

    
    int iFace, le, re;
    //! Get the face indexes of this boundary condition region.
    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();

    RDouble molecularWeight = one;
    using namespace GAS_SPACE;
    RDouble coefficientofstateEquation = gas->GetCoefficientOfStateEquation();

    RDouble rGhost = 0.0, pGhost = 0.0, tGhost = 0.0;    //! rho/pressure/temprature on ghost cell.
    RDouble *qL = NULL;
    if (nChemical == 1)
    {
        qL = new RDouble [nEquation];
    }

    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        iFace = *iter;
        le = leftCellofFace[iFace];
        re = rightCellofFace[iFace];

        for (int m = 0; m < nEquation; ++ m)
        {
            q[m][re] = q[m][le];
        }

        q[IU][re] = - q[IU][le] + two * faceVelocityX[iFace];
        q[IV][re] = - q[IV][le] + two * faceVelocityY[iFace];
        q[IW][re] = - q[IW][le] + two * faceVelocityZ[iFace];

        if (nChemical == 1)
        {
            for (int m = 0; m < nEquation; ++ m)
            {
                qL[m] = q[m][le];
            }
            molecularWeight = gas->ComputeMolecularWeightReciprocal(qL);
        }

        RDouble tL = (1.0 / coefficientofstateEquation) * q[IP][le] / q[IR][le];

        //! Pressure keep the same on the wall normal direction.
        pGhost = q[IP][le];

        if (wallTempArrayExist)
        {
            wallTemperature = wallTempArray[countFace] / refDimensionalTemperature;
            countFace ++;
        }

        //! Extrapolate temperature from inside to outside.
        tGhost = two * wallTemperature - tL;
        if (tGhost <= 0.0) tGhost = wallTemperature;

        if (nChemical == 1)
        {
            tGhost = wallTemperature;
        }

        //! Then compute the density from the extrapolated temperature.
        rGhost = pGhost / (coefficientofstateEquation * tGhost * molecularWeight);

        q[IR][re] = rGhost;
        q[IP][re] = pGhost;
    }

    if (nChemical == 1)
    {
        delete [] qL;
    }
}

void NSSolverUnstruct::ViscousAdiabaticWallBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct)
{
    // Get the basic variables for boundary calculating.
    UnstructGrid *grid   = UnstructGridCast(gridIn);
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();
    RDouble *faceVelocityX = grid->GetFaceVelocityX();
    RDouble *faceVelocityY = grid->GetFaceVelocityY();
    RDouble *faceVelocityZ = grid->GetFaceVelocityZ();
    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    RDouble **t   = reinterpret_cast< RDouble ** > (grid->GetDataPtr("t"));

    //! Get the number of equations.
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    Data_Param *bcData = bcRegionUnstruct->GetBCParamDataBase();

    RDouble uWall = 0.0;
    RDouble vWall = 0.0;
    RDouble wWall = 0.0;

    if(bcData)
    {
        if (bcData->IsExist("uWall", PHDOUBLE, 1))
        {
            bcData->GetData("uWall", &uWall, PHDOUBLE, 1);
            bcData->GetData("vWall", &vWall, PHDOUBLE, 1);
            bcData->GetData("wWall", &wWall, PHDOUBLE, 1);
        }
    }    

    int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
    int isAle = GlobalDataBase::GetIntParaFromDB("codeOfAleModel");

    int iFace, le, re;
    //! Get the face indexes of this boundary condition region.
    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();

    RDouble molecularWeight = one;
    using namespace GAS_SPACE;
    RDouble coefficientofstateEquation = gas->GetCoefficientOfStateEquation();
    RDouble rGhost = 0.0, pGhost = 0.0, tGhost = 0.0;
    RDouble *qL = NULL;
    if (nChemical == 1)
    {
        qL = new RDouble [nEquation];
    }

    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        iFace = *iter;
        le = leftCellofFace[iFace];
        re = rightCellofFace[iFace];

        for (int m = 0; m < nEquation; ++ m)
        {
            q[m][re] = q[m][le];
        }

        RDouble velocityXWall = uWall;
        RDouble velocityYWall = vWall;
        RDouble velocityZWall = wWall;
        if (isUnsteady && isAle)
        {
            velocityXWall = faceVelocityX[iFace] + uWall;
            velocityYWall = faceVelocityY[iFace] + vWall;
            velocityZWall = faceVelocityZ[iFace] + wWall;
    }

        q[IU][re] = - q[IU][le] + two * velocityXWall;
        q[IV][re] = - q[IV][le] + two * velocityYWall;
        q[IW][re] = - q[IW][le] + two * velocityZWall;
    }
    if (nChemical == 1)
    {
        for (int m = 0; m < nEquation; ++ m)
        {
            qL[m] = q[m][le];
        }
        molecularWeight = gas->ComputeMolecularWeightReciprocal(qL);

        tGhost = t[ITT][le];
        pGhost = q[IP][le];
        //! Then compute the density from the extrapolated temperature.
        rGhost = pGhost / (coefficientofstateEquation * tGhost * molecularWeight);

        q[IR][re] = rGhost;
        q[IP][re] = pGhost;
    }
}

void NSSolverUnstruct::FarfieldBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct)
{
    //! Get the basic variables for boundary calculating.
    UnstructGrid *grid   = UnstructGridCast(gridIn);
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();
    RDouble   *xfn       = grid->GetFaceNormalX();
    RDouble   *yfn       = grid->GetFaceNormalY();
    RDouble   *zfn       = grid->GetFaceNormalZ();
    RDouble   *vgn       = grid->GetFaceNormalVelocity();
    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    RDouble *gama = reinterpret_cast< RDouble *  > (grid->GetDataPtr("gama"));

    //! Get the number of equations.
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    RDouble refGama = parameters->GetRefGama();

    int nEquation = GetNumberOfEquations();

    RDouble *primitiveVarFarfield = new RDouble [nEquation];
    Data_Param *bcData = bcRegionUnstruct->GetBCParamDataBase();
    bcData->GetData("primitiveVarFarfield", primitiveVarFarfield, PHDOUBLE, nEquation);

    int iFace, le, re;
    //! Get the face indexes of this boundary condition region.
    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();

    RDouble rInside, uInside, vInside, wInside, pInside, cInside, vNormalInside, vNormalOutside;
    RDouble gamm1, velocityInside;
    RDouble vTangentX, vTangentY, vTangentZ, entropyS;
    RDouble R_Minus, R_Plus;
    RDouble rBC, uBC, vBC, wBC, pBC, vNormalBC, cBC;

    //! Infinite variables.
    RDouble rOutside = primitiveVarFarfield[IR];
    RDouble uOutside = primitiveVarFarfield[IU];
    RDouble vOutside = primitiveVarFarfield[IV];
    RDouble wOutside = primitiveVarFarfield[IW];
    RDouble pOutside = primitiveVarFarfield[IP];
    RDouble cOutside = sqrt(ABS(refGama * pOutside / rOutside));

    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        iFace = *iter;
        le = leftCellofFace[iFace];
        re = rightCellofFace[iFace];

        for (int m = 0; m < nEquation; ++ m)
        {
            q[m][re] = q[m][le];
        }

        //! Inner point variables.
        rInside = q[IR][le];
        uInside = q[IU][le];
        vInside = q[IV][le];
        wInside = q[IW][le];
        pInside = q[IP][le];

        vNormalOutside = xfn[iFace] * uOutside + yfn[iFace] * vOutside + zfn[iFace] * wOutside - vgn[iFace];
        vNormalInside  = xfn[iFace] * uInside  + yfn[iFace] * vInside  + zfn[iFace] * wInside  - vgn[iFace];

        cInside = sqrt(ABS(gama[le] * pInside / rInside));
        gamm1   = gama[le] - one;

        velocityInside = sqrt(uInside * uInside + vInside * vInside + wInside * wInside);

        //! Supersonic.
        if (velocityInside > cInside)
        {
            if (vNormalInside >= 0.0)
            {
                for (int m = 0; m < nEquation; ++ m)
                {
                    q[m][re] = q[m][le];
                }
            }
            else
            {
                for (int m = 0; m < nEquation; ++ m)
                {
                    q[m][re] = primitiveVarFarfield[m];
                }
            }
            continue;
        }

         //! Subsonic.
        R_Plus  = vNormalInside + 2.0 * cInside / gamm1;
        R_Minus = vNormalOutside - 2.0 * cOutside / gamm1;
        vNormalBC =  half * (R_Plus + R_Minus);
        cBC =  0.25 * (R_Plus - R_Minus) * gamm1;
        if (vNormalBC >= 0.0)
        {
            entropyS = pInside / pow(rInside, gama[le]);
            vTangentX  = uInside - xfn[iFace] * vNormalInside;
            vTangentY  = vInside - yfn[iFace] * vNormalInside;
            vTangentZ  = wInside - zfn[iFace] * vNormalInside;
        }
        else
        {
            //! Inlet.
            entropyS = pOutside / pow(rOutside, gama[le]);
            vTangentX  = uOutside - xfn[iFace] * vNormalOutside;
            vTangentY  = vOutside - yfn[iFace] * vNormalOutside;
            vTangentZ  = wOutside - zfn[iFace] * vNormalOutside;
        }

        rBC = pow((cBC * cBC / (entropyS * gama[le])), one / gamm1);
        uBC = vTangentX + xfn[iFace] * vNormalBC;
        vBC = vTangentY + yfn[iFace] * vNormalBC;
        wBC = vTangentZ + zfn[iFace] * vNormalBC;
        pBC = cBC * cBC * rBC / gama[le];

        q[IR][re] = rBC;
        q[IU][re] = uBC;
        q[IV][re] = vBC;
        q[IW][re] = wBC;
        q[IP][re] = pBC;
    }
}

void NSSolverUnstruct::PreconFarfieldBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct)
{
    //! Get the basic variables for boundary calculating.
    UnstructGrid *grid   = UnstructGridCast(gridIn);
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();
    RDouble   *xfn       = grid->GetFaceNormalX();
    RDouble   *yfn       = grid->GetFaceNormalY();
    RDouble   *zfn       = grid->GetFaceNormalZ();
    RDouble   *vgn       = grid->GetFaceNormalVelocity();
    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    RDouble *gama = reinterpret_cast< RDouble *  > (grid->GetDataPtr("gama"));

    //! Get the number of equations.
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    RDouble refMachNumber = parameters->GetRefMachNumber();

    int nEquation = GetNumberOfEquations();

    RDouble *primitiveVarFarfield = new RDouble [nEquation];
    Data_Param *bcData = bcRegionUnstruct->GetBCParamDataBase();
    bcData->GetData("primitiveVarFarfield", primitiveVarFarfield, PHDOUBLE, nEquation);

    int iFace, le, re;
    //! Get the face indexes of this boundary condition region.
    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();

    RDouble rin,uin,vin,win,pin,tin,vni;
    RDouble too, vno;
    RDouble rb,ub,vb,wb,pb,tb;

    RDouble mach2 = refMachNumber * refMachNumber;

    //! Infinite variables.
    RDouble roo = primitiveVarFarfield[IR];
    RDouble uoo = primitiveVarFarfield[IU];
    RDouble voo = primitiveVarFarfield[IV];
    RDouble woo = primitiveVarFarfield[IW];
    RDouble poo = primitiveVarFarfield[IP];    

    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        iFace = *iter;
        le = leftCellofFace[iFace];
        re = rightCellofFace[iFace];

        for (int m = 0; m < nEquation; ++ m)
        {
            q[m][re] = q[m][le];
        }

        //! Inner point variables.
        rin = q[IR][le];
        uin = q[IU][le];
        vin = q[IV][le];
        win = q[IW][le];
        pin = q[IP][le];
        tin = gama[le] * mach2 * pin / rin;
        too = gama[le] * mach2 * poo / roo;

        vno = xfn[iFace] * uoo + yfn[iFace] * voo + zfn[iFace] * woo - vgn[iFace];
        vni = xfn[iFace] * uin + yfn[iFace] * vin + zfn[iFace] * win - vgn[iFace];

        //! Subsonic.
        if (refMachNumber < 1.0)
        {
            if (vni < -1.e-10)
            {
                ub  = uoo;
                vb  = voo;
                wb  = woo;
                pb  = pin;
                tb  = too;
                rb  = gama[le] * mach2 * pb / tb;
            }
            else
            {
                ub  = uin;
                vb  = vin;
                wb  = win;
                pb  = poo;
                tb  = tin;
                rb  = gama[le] * mach2 * pb / tb;
            }
        }
        else
        {
            if (vni >= 0.0)
            {
                ub  = uin;
                vb  = vin;
                wb  = win;
                pb  = pin;
                rb  = rin;
            }
            else
            {
                ub  = uoo;
                vb  = voo;
                wb  = woo;
                pb  = poo;
                rb  = roo;
            }
        }

        q[IR][re] = rb;
        q[IU][re] = ub;
        q[IV][re] = vb;
        q[IW][re] = wb;
        q[IP][re] = pb;
    }
}

void NSSolverUnstruct::InflowBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct)
{
    //! Get the basic variables for boundary calculating.
    UnstructGrid *grid   = UnstructGridCast(gridIn);
    int *rightCellofFace = grid->GetRightCellOfFace();
    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));

    //! Get the number of equations.
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    RDouble *primitiveVarFarfield = new RDouble [nEquation];
    Data_Param *bcData = bcRegionUnstruct->GetBCParamDataBase();
    bcData->GetData("primitiveVarFarfield", primitiveVarFarfield, PHDOUBLE, nEquation);

    int iFace, re;
    //! Get the face indexes of this boundary condition region.
    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();
   
    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        iFace = *iter;
        re = rightCellofFace[iFace];
        for (int m = 0; m < nEquation; ++ m)
        {
            q[m][re] = primitiveVarFarfield[m];
        }
    }
}

void NSSolverUnstruct::PressureOutletBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct)
{
    //! Get the basic variables for boundary calculating.
    UnstructGrid *grid   = UnstructGridCast(gridIn);
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();

    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));

    //! Get the number of equations.
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    Data_Param *bcData = bcRegionUnstruct->GetBCParamDataBase();

    RDouble staticPressure;
    bcData->GetData("staticPressure", &staticPressure, PHDOUBLE, 1);

    //! Dimensional.
    RDouble refDimensionalDensity = parameters->GetRefDimensionalDensity();
    RDouble refDimensionalSonicSpeed = parameters->GetRefDimensionalSonicSpeed();
    RDouble refMachNumber = parameters->GetRefMachNumber();
    RDouble velocity = refMachNumber * refDimensionalSonicSpeed;
    RDouble dynamic_pressure = refDimensionalDensity * velocity * velocity;

    staticPressure /= dynamic_pressure;

    int iFace, le, re;
    //! Get the face indexes of this boundary condition region.
    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();

    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        iFace = *iter;
        le = leftCellofFace[iFace];
        re = rightCellofFace[iFace];
        for (int m = 0; m < nEquation; ++ m)
        {
            q[m][re] = q[m][le];
        }
        q[IP][re] = staticPressure;
    }
}

void NSSolverUnstruct::MassFlowOutletBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct)
{
    calculateBCArea(gridIn, bcRegionUnstruct);
    calculateMassFluxRatio(gridIn, bcRegionUnstruct);

    using namespace IDX;
    //! Dimensional.
    Param_NSSolverUnstruct *parameters = GetControlParameters();  
    RDouble refMachNumber = parameters->GetRefMachNumber();
    RDouble refDimensionalDensity = parameters->GetRefDimensionalDensity();
    RDouble refGama = parameters->GetRefGama();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    Data_Param *bcData = bcRegionUnstruct->GetBCParamDataBase();
    RDouble massFluxRatio, massFlow, BCTotalArea;
    bcData->GetData("massFluxRatio", &massFluxRatio, PHDOUBLE, 1);
    bcData->GetData("massFlow", &massFlow, PHDOUBLE, 1);
    bcData->GetData("BCTotalArea", &BCTotalArea, PHDOUBLE, 1);

    //! Nondimensionalization.
    RDouble coefficientOfStateEquation = 1.0/(refGama * refMachNumber * refMachNumber);

    // fluxDensity : rho * u;
    RDouble fluxDensity = massFlow/(refDimensionalDensity * refDimensionalDensity * refMachNumber)/BCTotalArea;

    RDouble mfTbCoeff = 0.5 * (refGama - 1.0) / refGama / coefficientOfStateEquation * 
        (1.0 - massFluxRatio * massFluxRatio);

    //! Get the basic variables for boundary calculating.
    UnstructGrid *grid   = UnstructGridCast(gridIn);
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();
    RDouble *xfn = grid->GetFaceNormalX();
    RDouble *yfn = grid->GetFaceNormalY();
    RDouble *zfn = grid->GetFaceNormalZ();
    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    
    //! Get the face indexes of this boundary condition region.
    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();
    int iFace, le, re;

    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        iFace = *iter;
        le = leftCellofFace[iFace];
        re = rightCellofFace[iFace];

        for (int m = 0; m < nEquation; ++ m)
        {
            q[m][re] = q[m][le];
        }

        //! Calculate the temperature of the inside cell.
        RDouble staticTemperature = q[IP][le] / q[IR][le] / coefficientOfStateEquation;
        RDouble tb   = staticTemperature + mfTbCoeff * (q[IU][le] * q[IU][le] +
            q[IV][le] * q[IV][le] + q[IW][le] * q[IW][le]);

        if(tb < TINY) 
        {
            tb = staticTemperature;
        }

        RDouble pb   = q[IP][le] * pow(tb / staticTemperature, refGama / (refGama - 1.0));
        RDouble rb = pb / tb / coefficientOfStateEquation;

        RDouble ub   = q[IU][le] * massFluxRatio;
        RDouble vb   = q[IV][le] * massFluxRatio;
        RDouble wb   = q[IW][le] * massFluxRatio;

        //! Judge the unusual situation, Vn < 0.
        RDouble Vn  = ub * xfn[iFace] + vb * yfn[iFace] + wb * zfn[iFace];        
        if(Vn < 0.0)
        {
            Vn = fluxDensity / rb;
            ub = Vn * xfn[iFace];
            vb = Vn * yfn[iFace];
            wb = Vn * zfn[iFace];
        } 

        //! If the outlet is supersonic flow.
        RDouble cb = sqrt(refGama * coefficientOfStateEquation * tb);
        if(Vn >= cb)
        {
            rb = q[IR][le];
            ub = q[IU][le];
            vb = q[IV][le];
            wb = q[IW][le];
            pb = q[IP][le];
        }

        //! Allocate value to the ghost cell.
        q[IR][re] = rb;
        q[IU][re] = ub;
        q[IV][re] = vb;
        q[IW][re] = wb;
        q[IP][re] = pb;
    }
}

void NSSolverUnstruct::calculateBCArea(Grid *gridIn, UnstructBC *bcRegionUnstruct)
{
    Data_Param *bcData = bcRegionUnstruct->GetBCParamDataBase();

    //BC total area does not change with calculating, so it only calculates once.
    if (bcData->IsExist("BCTotalArea", PHDOUBLE, 1))
    {
        return;
    }

    using namespace PHMPI;
    string bcName = bcRegionUnstruct->GetBCName();
    RDouble BCTotalArea = 0.0;

    //! One processor may owe some zones, so calculate all BC area in one processor.
    int nLocalZones = GetNumberofLocalZones();
    for (int iZone = 0; iZone < nLocalZones; iZone++)
    {
        int iZoneID = GetLocalZoneIDToGlobalZoneID(iZone);
        Grid * iGrid = PHSPACE::GetGrid(iZoneID);
        UnstructGrid *grid = UnstructGridCast(iGrid);      
        
        UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
        int nBCRegion = unstructBCSet->GetnBCRegion();

        for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
        {
            UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
            string jBCName = bcRegion->GetBCName();

            //! Judge the same BC region according the BC name.
            if (jBCName == bcName)
            {
                int iFace;
                vector<int> *faceIndex = bcRegion->GetFaceIndex();
                RDouble   *area    = grid->GetFaceArea();
                for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
                {
                    iFace = *iter;
                    BCTotalArea += area[iFace];
                }
            }
        }
    }

    //! Get the total area in all processors with the same BC name.
    RDouble GLB_BCTotalArea;
    PH_AllReduce(&BCTotalArea, &GLB_BCTotalArea, 1, MPI_SUM);
    BCTotalArea = GLB_BCTotalArea;

    for (int iZone = 0; iZone < nLocalZones; iZone++)
    {
        int iZoneID = GetLocalZoneIDToGlobalZoneID(iZone);
        Grid * iGrid = PHSPACE::GetGrid(iZoneID);
        UnstructGrid *grid = UnstructGridCast(iGrid); 

        UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
        int nBCRegion = unstructBCSet->GetnBCRegion();

        for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
        {
            UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
            string jBCName = bcRegion->GetBCName();

            if (jBCName == bcName)
            {
                Data_Param *bcData = bcRegion->GetBCParamDataBase();
                bcData->UpdateData("BCTotalArea", &BCTotalArea, PHDOUBLE, 1);
            }
        }                       
    }   
}

void NSSolverUnstruct::calculateMassFluxRatio(Grid *gridIn, UnstructBC *bcRegionUnstruct)
{
    using namespace PHMPI;
    string bcName = bcRegionUnstruct->GetBCName();
    Data_Param * bcData = bcRegionUnstruct->GetBCParamDataBase();
    RDouble massFlow = 0.0;
    bcData->GetData("massFlow", &massFlow, PHDOUBLE, 1);
    RDouble massFlux = 0.0;

    //! One processor may owe some zones, so calculate all BC mass flux in one processor.
    //! In fact, massflux need calculte once in one step, but I dont find a good method to judge if it is
    //! Calculated or not, so it will calculate more times, may cost a little more time.
    int nLocalZones = GetNumberofLocalZones();
    for (int iZone = 0; iZone < nLocalZones; iZone++)
    {
        int iZoneID = GetLocalZoneIDToGlobalZoneID(iZone);
        Grid * iGrid = PHSPACE::GetGrid(iZoneID);
        UnstructGrid *grid = UnstructGridCast(iGrid);      

        UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
        int nBCRegion = unstructBCSet->GetnBCRegion();

        RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
        for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
        {
            UnstructBC *bcRegion = unstructBCSet->GetBCRegion(iBCRegion);
            string jBCName = bcRegion->GetBCName();
            
            //! Judge the same BC region according the BC name.
            if (jBCName == bcName)
            {
                int iFace;
                vector<int> *faceIndex = bcRegion->GetFaceIndex();
                RDouble *area = grid->GetFaceArea();
                int *leftCellofFace  = grid->GetLeftCellOfFace();
                RDouble *xfn = grid->GetFaceNormalX();
                RDouble *yfn = grid->GetFaceNormalY();
                RDouble *zfn = grid->GetFaceNormalZ();
                for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
                {
                    iFace = *iter;
                    int le = leftCellofFace[iFace];

                    RDouble vn_c1 = q[IU][le] * xfn[iFace] + q[IV][le] * yfn[iFace] + q[IW][le] * zfn[iFace];
                    massFlux += q[IR][le] * vn_c1 * area[iFace];
                }
            }
        }
    }
    
    //! Get the total area in all processors with the same BC name.
    RDouble GLB_MassFlux;
    PH_AllReduce(&massFlux, &GLB_MassFlux, 1, MPI_SUM);

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    RDouble refMachNumber = parameters->GetRefMachNumber();
    RDouble refDimensionalDensity = parameters->GetRefDimensionalDensity();
    RDouble refDimensionalSonicSpeed = parameters->GetRefDimensionalSonicSpeed();
    
    RDouble gridScaleFactor = GlobalDataBase::GetDoubleParaFromDB("gridScaleFactor");

    //! Calculte the dimensional mass flux in this BC name region.
    GLB_MassFlux *= refDimensionalDensity * refDimensionalSonicSpeed * refMachNumber * gridScaleFactor * gridScaleFactor;

    RDouble massFluxRatio =  massFlow / (GLB_MassFlux + TINY);
    massFluxRatio = MAX(0.8, MIN(massFluxRatio, 1.2));

    bcData->UpdateData("massFluxRatio", &massFluxRatio, PHDOUBLE, 1);
}

void NSSolverUnstruct::MassFlowInletBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct)
{
    calculateBCArea(gridIn, bcRegionUnstruct);

    //! Get the basic mesh and flow field information.
    UnstructGrid *grid   = UnstructGridCast(gridIn);
   
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();

    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));

    //! Get the number of equations.
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    RDouble *primitiveVarFarfield = parameters->GetPrimitiveVarFarfield();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    int iFace = 0, le = 0, re = 0;
    //! Get the face indexes of this boundary condition region.
    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();

    Data_Param *bcData = bcRegionUnstruct->GetBCParamDataBase();

    //! Get the basic variables for mass flow in boundary.
    RDouble direction_inlet[3] = {0.0}, massFlow = 0.0, totalTemperature = 0.0;
    bcData->GetData("massFlow", &massFlow, PHDOUBLE, 1);
    bcData->GetData("totalTemperature", &totalTemperature, PHDOUBLE, 1);
    bcData->GetData("direction_inlet", direction_inlet, PHDOUBLE, 3);

    RDouble BCTotalArea = 0.0;
    if (bcData->IsExist("BCTotalArea", PHDOUBLE, 1))
    {
        bcData->GetData("BCTotalArea", &BCTotalArea, PHDOUBLE, 1);
    }else
    {
        TK_Exit::ExceptionExit("BCTotalArea data does not exist!\n");
    }

    massFlow /= BCTotalArea;

    //! Dimensional.
    RDouble refMachNumber = parameters->GetRefMachNumber();
    RDouble refDimensionalDensity = parameters->GetRefDimensionalDensity();
    RDouble refDimensionalSonicSpeed = parameters->GetRefDimensionalSonicSpeed();
    RDouble refGama = parameters->GetRefGama();
    RDouble refDimensionalTemperature = parameters->GetRefDimensionalTemperature();

    //! Nondimensionalization.
    totalTemperature /= refDimensionalTemperature;
    RDouble coefficientOfStateEquation = 1.0/(refGama * refMachNumber * refMachNumber);

    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        iFace = *iter;
        le = leftCellofFace[iFace];
        re = rightCellofFace[iFace];

        for (int m = 0; m < nEquation; ++ m)
        {
            q[m][re] = primitiveVarFarfield[m];
        }

        q[IR][re]  = q[IR][le];
        RDouble Vnb   = massFlow/ (refDimensionalDensity*refDimensionalSonicSpeed*refMachNumber) / q[IR][le];  //注意：此处massFlow为单位面积流量

        q[IU][re] = Vnb * direction_inlet[0];
        q[IV][re] = Vnb * direction_inlet[1];
        q[IW][re] = Vnb * direction_inlet[2];

        //! Square of the stagnation sound speed.
        RDouble c2   = refGama * coefficientOfStateEquation * totalTemperature; 
        //! Square of the critical sound speed.
        RDouble c2_ct = 2.0 * c2 / (refGama+1.0); 
        //! Square of the Laval number.
        RDouble lamd2 = Vnb * Vnb / c2_ct;
        RDouble tb    = totalTemperature * (1.0 - (refGama - 1.0)/(refGama + 1.0) * lamd2);
        q[IP][re] = q[IR][le] * coefficientOfStateEquation * tb;
    }
}

void NSSolverUnstruct::PressureInletBCRiemannRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct)
{
    //! Get the basic variables for boundary calculating.
    UnstructGrid *grid   = UnstructGridCast(gridIn);
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();
    RDouble *xfn = grid->GetFaceNormalX();
    RDouble *yfn = grid->GetFaceNormalY();
    RDouble *zfn = grid->GetFaceNormalZ();
    RDouble **q  = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));

    //! Get the number of equations.
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = GetNumberOfEquations();

    int iFace = 0, le = 0, re = 0;
    //! Get the face indexes of this boundary condition region.
    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();

    Data_Param *bcData = bcRegionUnstruct->GetBCParamDataBase();

    RDouble direction_inlet[3], totalPressure, totalTemperature;
    bcData->GetData("totalPressure", &totalPressure, PHDOUBLE, 1);
    bcData->GetData("totalTemperature", &totalTemperature, PHDOUBLE, 1);
    bcData->GetData("direction_inlet", direction_inlet, PHDOUBLE, 3);

    RDouble *primitiveVarFarfield = new RDouble [nEquation];
    bcData->GetData("primitiveVarFarfield", primitiveVarFarfield, PHDOUBLE, nEquation);

    //! Dimensional.
    RDouble refDimensionalTemperature = parameters->GetRefDimensionalTemperature();
    RDouble refMachNumber = parameters->GetRefMachNumber();
    RDouble refDimensionalDensity = parameters->GetRefDimensionalDensity();
    RDouble refDimensionalSonicSpeed = parameters->GetRefDimensionalSonicSpeed();
    RDouble velocity = refMachNumber * refDimensionalSonicSpeed;
    RDouble dynamic_pressure = refDimensionalDensity * velocity * velocity;
    RDouble refGama = parameters->GetRefGama();
    RDouble gama1 = refGama - 1;

    //! Nondimensionalization.
    totalPressure /= dynamic_pressure;
    totalTemperature /= refDimensionalTemperature;
    RDouble coefficientOfStateEquation = 1.0/(refGama * refMachNumber * refMachNumber);

    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        iFace = *iter;
        le = leftCellofFace[iFace];
        re = rightCellofFace[iFace];

        for (int m = 0; m < nEquation; ++ m)
        {
            q[m][re] = primitiveVarFarfield[m];
        }

        //! Get inner cell flow field information.
        RDouble rInside = q[IR][le];
        RDouble V2Inside = q[IU][le] * q[IU][le] + q[IV][le] * q[IV][le]  + q[IW][le] * q[IW][le];
        RDouble pInside = q[IP][le];

        //! Calculate the total energy, total enthalpy and sound speed.
        RDouble totalEnergy = pInside / (rInside * gama1) + 0.5 * V2Inside;
        RDouble totalEnthalpy = (refGama * coefficientOfStateEquation / gama1) * totalTemperature;
        RDouble c2 = refGama * pInside / rInside;

        //! Calculate the riemann invariants R+.
        RDouble Riemann = 2.0 * sqrt(c2) / gama1;
        Riemann += q[IU][le] * xfn[iFace] + q[IV][le] * yfn[iFace] + q[IW][le] * zfn[iFace];

        RDouble theta = xfn[iFace] * direction_inlet[0] + yfn[iFace] * direction_inlet[1] + 
            zfn[iFace] * direction_inlet[2];

        RDouble totalc2 = gama1 * (totalEnthalpy - (totalEnergy + pInside / rInside)+0.5 * V2Inside) + c2;

        RDouble a =  1.0 + 0.5 * gama1 * theta * theta;
        RDouble b = -1.0 * gama1 * theta * Riemann;
        RDouble c =  0.5 * gama1 * Riemann * Riemann - 2.0 * totalc2 / gama1;

        RDouble d = b * b - 4.0 * a * c;
        d = sqrt(max(0.0, d));
        RDouble Vin   = (-b + d)/(2.0*a);
        Vin   = max(0.0, Vin);
        V2Inside = Vin * Vin;

        c2 = totalc2 - 0.5 * gama1 * V2Inside;

        //! Fix the mach number in case of mach > 1.0 .
        RDouble Mach2 = V2Inside / c2;
        Mach2 = min(1.0, Mach2);
        V2Inside = Mach2 * c2;
        Vin = sqrt(V2Inside);
        c2 = totalc2 - 0.5 * gama1 * V2Inside;

        RDouble tInside = c2 / (refGama * coefficientOfStateEquation);
        pInside = totalPressure * pow((tInside / totalTemperature), refGama / gama1);
        rInside = pInside/(coefficientOfStateEquation * tInside);

        q[IR][re] = rInside;
        q[IU][re] = Vin * direction_inlet[0];
        q[IV][re] = Vin * direction_inlet[1];
        q[IW][re] = Vin * direction_inlet[2];
        q[IP][re] = pInside;
    }
}

void NSSolverUnstruct::SwapNeighborsDQNonblocking(Grid *gridIn, FieldProxy *dqProxy)
{
    int currentProcessor = PHMPI::GetCurrentProcessorID();
    int nZones = PHMPI::GetNumberofGlobalZones();
    vector< PH_Request > requestContainer;
    vector< vector< DataContainer * > > receivedDataBuffer;
    vector< vector< DataContainer * > > sendDataBuffer;
    receivedDataBuffer.resize(nZones);
    sendDataBuffer.resize(nZones);
    int level = gridIn->GetLevel();

    //! Step 0: Compressing data firstly.
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int sendProcessor = PHMPI::GetZoneProcessorID(iZone);

        ZoneNeighbor * globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
        std::size_t numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();
        if(currentProcessor == sendProcessor)
        {
            sendDataBuffer[iZone].resize(numberOfNeighbor);
        }

        for (std::size_t iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
        {
            int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);

            int receiveProcessor = PHMPI::GetZoneProcessorID(neighborZone);

            if(currentProcessor != sendProcessor && currentProcessor != receiveProcessor)
            {
                continue;
            }

            //! Allocate and compress the buffers for sending.
            DataContainer * sendBuffer = 0;
            if(currentProcessor == sendProcessor)
            {
                sendBuffer = new DataContainer;
                sendDataBuffer[ iZone ][iNeighbor] = sendBuffer;

                //! Compress the send information into the actkey.
                UnstructGrid * grid = UnstructGridCast(PHSPACE::GetGrid(iZone, level));
                CompressDQ(sendBuffer, dqProxy, grid, iZone, neighborZone, 5);
            }
        }
    }

    //! Step 1: Communication.
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        ZoneNeighbor * globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
        std::size_t numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();

        //! Communicating.
        for (std::size_t iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
        {
            int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);

            int sendProcessor    = PHMPI::GetZoneProcessorID(iZone);
            int receiveProcessor = PHMPI::GetZoneProcessorID(neighborZone);

            if(currentProcessor != sendProcessor && currentProcessor != receiveProcessor)
            {
                continue;
            }

            //! Allocate the buffers for sending and receiving.
            DataContainer * receivedBuffer = 0;
            DataContainer * sendBuffer = 0;
            if(currentProcessor == receiveProcessor)
            {
                receivedBuffer = new DataContainer;
                receivedDataBuffer[ iZone ].push_back(receivedBuffer);
            }
            if(currentProcessor == sendProcessor)
            {
                sendBuffer = sendDataBuffer[iZone][iNeighbor];
            }

            int tag = iZone;

            //! Communication.
            if (currentProcessor == sendProcessor)
            {
                if (sendProcessor != receiveProcessor)
                {
                    streamsize nlen = sendBuffer->size();

                    //! Send the data to neighbors.
                    send(sendBuffer, receiveProcessor, nlen, requestContainer, tag);
                }
            }
            if (currentProcessor == receiveProcessor)
            {
                //! Get the neighbor order: which order that 'iZone' on neighbors of neighborZone.
                int neighborOrer = -1;
                ZoneNeighbor * globalNeighborZonesTemp = zoneConnectivity->GetZoneNeighbor(neighborZone);
                std::size_t numberOfNeighborTemp = globalNeighborZonesTemp->GetNumberOfNeighbor();

                for (std::size_t neighborID = 0; neighborID < numberOfNeighborTemp; ++ neighborID)
                {
                    int neighborZoneTemp = globalNeighborZonesTemp->GetZoneIndexOfNeighbor(neighborID);
                    if(neighborZoneTemp == iZone)
                    {
                        neighborOrer = static_cast<int>(neighborID);
                        break;
                    }
                }

                ASSERT(neighborOrer != -1);

                if (sendProcessor != receiveProcessor)
                {
                    //! The data length of the received data is same to the length that send to the neighbor.
                    CharVecSizeType nlen = sendDataBuffer[neighborZone][neighborOrer]->size();

                    //! Receive data from neighbors.
                    receive(receivedBuffer, sendProcessor, nlen, requestContainer, tag);
                }
                else
                {
                    SWAP(receivedDataBuffer[ iZone ].back(), sendDataBuffer[iZone][iNeighbor]);
                }
            }
        }
    }

    //! Step2: MPI waiting for Non-blocking communications.
    if(PHMPI::GetNumberOfProcessor() > 1)
    {
        PH_Wait(static_cast<int>(requestContainer.size()), &(requestContainer[0]));
    }

    //! Step3: Translate the data container.
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        ZoneNeighbor * globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
        std::size_t numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();

        int count = 0;
        for (std::size_t iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
        {
            int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);

            int receiveProcessor = PHMPI::GetZoneProcessorID(neighborZone);

            if (currentProcessor == receiveProcessor)
            {
                DataContainer * receiveData = receivedDataBuffer[ iZone ][ count ];

                //! Because get into here when currentProcessor==receiveProcessor,
                //! so, use zoneGridSolver of neighborZone, this zone is on the current processor!

                //! Decompress the interface data from data-container.
                UnstructGrid * grid = UnstructGridCast(PHSPACE::GetGrid(neighborZone, level));
                DecompressDQ(receiveData, dqProxy, grid, neighborZone, iZone, 5);

                ++ count;
            }
        }
    }

    //! Step4: Free the buffers.
    for (std::size_t iDim = 0; iDim < receivedDataBuffer.size(); ++ iDim)
    {
        for (std::size_t jDim = 0; jDim < receivedDataBuffer[ iDim ].size(); ++ jDim)
        {
            delete receivedDataBuffer[ iDim ][ jDim ];
        }
    }
    receivedDataBuffer.clear();
    for (std::size_t iDim = 0; iDim < sendDataBuffer.size(); ++ iDim)
    {
        for (std::size_t jDim = 0; jDim < sendDataBuffer[ iDim ].size(); ++ jDim)
        {
            delete sendDataBuffer[ iDim ][ jDim ];
        }
    }
    sendDataBuffer.clear();
    requestContainer.clear();
}

void NSSolverUnstruct::CompressGradientandLimitToInterfaceOnlyForMixGrid(DataContainer *&dataContainer, Grid *gridIn, const int &neighborZoneIndex)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
    if (!interfaceInformation)
    {
        return;
    }

    int iNeighborZone                            = interfaceInformation->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = interfaceInformation->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForSend          = interfaceInformation->GetFaceIndexForSend(iNeighborZone);
    dataContainer->moveto_begin();
    
    RDouble **dqdx = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarX"));
    RDouble **dqdy = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarY"));
    RDouble **dqdz = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarZ"));

    RDouble **dtdx  = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradTemperatureX"));
    RDouble **dtdy  = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradTemperatureY"));
    RDouble **dtdz  = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradTemperatureZ"));

    int limitVector = GlobalDataBase::GetIntParaFromDB("limitVector");

    RDouble *limit = 0;
    RDouble **limit2D = 0;
    if(limitVector == 0)
    {
        //! Scalar
        limit = limiterOnlyForMixGrid->GetLimiter(0);
    }
    else
    {
        limit2D = new RDouble *[5];

        for (int m = 0; m < 5; ++ m)
        {
            limit2D[m] = limiterOnlyForMixGrid->GetLimiter(m);
        }
    }
    
    for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
    {
        int sourceCell;
        int iFace = interfaceIndexContainerForSend[iLocalFace];
        grid->GetSourceIndex(iFace, 1, sourceCell);
        for (int m = 0; m <= 4; ++ m)
        {
            PHWrite(dataContainer, dqdx[m][sourceCell]);
            PHWrite(dataContainer, dqdy[m][sourceCell]);
            PHWrite(dataContainer, dqdz[m][sourceCell]);
        }

        PHWrite(dataContainer, dtdx[0][sourceCell]);
        PHWrite(dataContainer, dtdy[0][sourceCell]);
        PHWrite(dataContainer, dtdz[0][sourceCell]);

        if(limitVector == 0)
        {
            //! Scalar
            PHWrite(dataContainer, limit[sourceCell]);
        }
        else
        {
            PHWrite(dataContainer, limit2D[0][sourceCell]);
            PHWrite(dataContainer, limit2D[1][sourceCell]);
            PHWrite(dataContainer, limit2D[2][sourceCell]);
            PHWrite(dataContainer, limit2D[3][sourceCell]);
            PHWrite(dataContainer, limit2D[4][sourceCell]);
        }
    }

    delete [] limit2D;
}

void NSSolverUnstruct::DecompressGradientandLimitToInterfaceOnlyForMixGrid(DataContainer *&dataContainer, Grid *gridIn, const int &neighborZoneIndex)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
    if (!interfaceInformation)
    {
        return;
    }

    int iNeighborZone                            = interfaceInformation->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = interfaceInformation->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForReceive       = interfaceInformation->GetFaceIndexForRecv(iNeighborZone);
    dataContainer->moveto_begin();
    
    RDouble **dqdx = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarX"));
    RDouble **dqdy = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarY"));
    RDouble **dqdz = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarZ"));

    RDouble **dtdx  = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradTemperatureX"));
    RDouble **dtdy  = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradTemperatureY"));
    RDouble **dtdz  = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradTemperatureZ"));

    int limitVector = GlobalDataBase::GetIntParaFromDB("limitVector");

    RDouble *limit = 0;
    RDouble **limit2D = 0;
    if(limitVector == 0)
    {
        //! Scalar
        limit = limiterOnlyForMixGrid->GetLimiter(0);
    }
    else
    {
        limit2D = new RDouble *[5];

        for (int m = 0; m < 5; ++ m)
        {
            limit2D[m] = limiterOnlyForMixGrid->GetLimiter(m);
        }
    }

    for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
    {
        int iFace = interfaceIndexContainerForReceive[iLocalFace];
        int targetCell;
        grid->GetTargetIndex(iFace, 1, targetCell);
        for (int m = 0; m <= 4; ++ m)
        {
            PHRead(dataContainer, dqdx[m][targetCell]);
            PHRead(dataContainer, dqdy[m][targetCell]);
            PHRead(dataContainer, dqdz[m][targetCell]);
        }

        PHRead(dataContainer, dtdx[0][targetCell]);
        PHRead(dataContainer, dtdy[0][targetCell]);
        PHRead(dataContainer, dtdz[0][targetCell]);

        if(limitVector == 0)
        {
            //! Scalar
            PHRead(dataContainer, limit[targetCell]);
        }
        else
        {
            PHRead(dataContainer, limit2D[0][targetCell]);
            PHRead(dataContainer, limit2D[1][targetCell]);
            PHRead(dataContainer, limit2D[2][targetCell]);
            PHRead(dataContainer, limit2D[3][targetCell]);
            PHRead(dataContainer, limit2D[4][targetCell]);
        }
    }
    delete [] limit2D;
}

void NSSolverUnstruct::CompressQlQrToInterfaceOnlyForMixGrid(DataContainer *&dataContainer, Grid *gridIn, const int &neighborZoneIndex)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
    if (!interfaceInformation)
    {
        return;
    }

    int iNeighborZone                            = interfaceInformation->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = interfaceInformation->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForSend          = interfaceInformation->GetFaceIndexForSend(iNeighborZone);
    int * interFace2BoundaryFace = interfaceInformation->GetInterFace2BoundaryFace();
    dataContainer->moveto_begin();

    RDouble **qL = InviscidfaceProxyOnlyForMixGrid->GetQL();

    for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
    {
        int iFace = interfaceIndexContainerForSend[iLocalFace];
        int FaceID = interFace2BoundaryFace[iFace];
        
        for (int m = 0; m < 5; ++ m)
        {
            PHWrite(dataContainer, qL[m][FaceID]);
        }
    }
}

void NSSolverUnstruct::DecompressQlQrToInterfaceOnlyForMixGrid(DataContainer *&dataContainer, Grid *gridIn, const int &neighborZoneIndex)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
    if (!interfaceInformation)
    {
        return;
    }

    int iNeighborZone                            = interfaceInformation->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = interfaceInformation->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForReceive       = interfaceInformation->GetFaceIndexForRecv(iNeighborZone);
    int * interFace2BoundaryFace = interfaceInformation->GetInterFace2BoundaryFace();
    dataContainer->moveto_begin();

    RDouble **qR = InviscidfaceProxyOnlyForMixGrid->GetQR();

    for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
    {
        int iFace = interfaceIndexContainerForReceive[iLocalFace];
        int FaceID = interFace2BoundaryFace[iFace];
        
        for (int m = 0; m < 5; ++ m)
        {
            PHRead(dataContainer, qR[m][FaceID]);
        }
    }
}

void NSSolverUnstruct::CompressFacefluxToInterfaceOnlyForMixGrid(DataContainer *&dataContainer, Grid *gridIn, const int &neighborZoneIndex)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
    if (!interfaceInformation)
    {
        return;
    }

    int gridtypeofcurrentZone = grid->Type();
    int gridtypeofneighborZone = PHMPI::GetZoneGridType(neighborZoneIndex);
    if(gridtypeofcurrentZone == gridtypeofneighborZone)
    {
        RDouble uploadData = 0.0;
        dataContainer->moveto_begin();
        PHWrite(dataContainer, uploadData);
        return;
    }

    int iNeighborZone                            = interfaceInformation->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = interfaceInformation->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForSend          = interfaceInformation->GetFaceIndexForSend(iNeighborZone);
    int * interFace2BoundaryFace = interfaceInformation->GetInterFace2BoundaryFace();
    dataContainer->moveto_begin();

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int viscousType = parameters->GetViscousType();

    if(viscousType <= INVISCID)
    {
        //! Euler
        RDouble **Inviscidflux = InviscidfaceProxyOnlyForMixGrid->GetFlux();

        for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
        {
            int iFace = interfaceIndexContainerForSend[iLocalFace];
            int FaceID = interFace2BoundaryFace[iFace];
            
            for (int m = 0; m < 5; ++ m)
            {
                PHWrite(dataContainer, Inviscidflux[m][FaceID]);
            }
        }
    }
    else
    {
        //! Laminar
        RDouble **Inviscidflux = InviscidfaceProxyOnlyForMixGrid->GetFlux();
        RDouble **Viscousflux = ViscousfaceProxyOnlyForMixGrid->GetFlux();
        
        for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
        {
            int iFace = interfaceIndexContainerForSend[iLocalFace];
            int FaceID = interFace2BoundaryFace[iFace];
            
            for (int m = 0; m < 5; ++ m)
            {
                PHWrite(dataContainer, Inviscidflux[m][FaceID]);
                PHWrite(dataContainer, Viscousflux[m][FaceID]);
            }
        }
    }
}

void NSSolverUnstruct::DecompressFacefluxToInterfaceOnlyForMixGrid(DataContainer *&dataContainer, Grid *gridIn, const int &neighborZoneIndex)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
    if (!interfaceInformation)
    {
        return;
    }

    int gridtypeofcurrentZone = grid->Type();
    int gridtypeofneighborZone = PHMPI::GetZoneGridType(neighborZoneIndex);
    if(gridtypeofcurrentZone == gridtypeofneighborZone)
    {
        RDouble downdData = 0.0;
        dataContainer->moveto_begin();
        PHRead(dataContainer, downdData);
        return;
    }

    int iNeighborZone                            = interfaceInformation->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = interfaceInformation->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForReceive       = interfaceInformation->GetFaceIndexForRecv(iNeighborZone);
    int * interFace2BoundaryFace = interfaceInformation->GetInterFace2BoundaryFace();
    dataContainer->moveto_begin();

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int viscousType = parameters->GetViscousType();

    RDouble coefofstrflux = GlobalDataBase::GetDoubleParaFromDB("coefofstrflux");
    RDouble coefofunstrflux = 1.0 - coefofstrflux;

    if(viscousType <= INVISCID)
    {
        //! Euler
        RDouble **Inviscidflux = InviscidfaceProxyOnlyForMixGrid->GetFlux();

        RDouble Inviscidflux_str[5] = {0.0};
        for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
        {
            int iFace = interfaceIndexContainerForReceive[iLocalFace];
            int FaceID = interFace2BoundaryFace[iFace];
            
            for (int m = 0; m < 5; ++ m)
            {
                PHRead(dataContainer, Inviscidflux_str[m]);

                Inviscidflux[m][FaceID] = (coefofstrflux * Inviscidflux_str[m] + coefofunstrflux * Inviscidflux[m][FaceID]);
            }
        }
    }
    else
    {
        //! Laminar
        RDouble **Inviscidflux = InviscidfaceProxyOnlyForMixGrid->GetFlux();
        RDouble **Viscousflux = ViscousfaceProxyOnlyForMixGrid->GetFlux();

        RDouble Inviscidflux_str[5];
        RDouble Viscousflux_str[5];
        for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
        {
            int iFace = interfaceIndexContainerForReceive[iLocalFace];
            int FaceID = interFace2BoundaryFace[iFace];
            
            for (int m = 0; m < 5; ++ m)
            {
                PHRead(dataContainer, Inviscidflux_str[m]);
                PHRead(dataContainer, Viscousflux_str[m]);

                Inviscidflux[m][FaceID] = (coefofstrflux * Inviscidflux_str[m] + coefofunstrflux * Inviscidflux[m][FaceID]);
                Viscousflux[m][FaceID]  = (coefofstrflux * Viscousflux_str[m]  + coefofunstrflux * Viscousflux[m][FaceID]);
            }
        }
    }
}

}


