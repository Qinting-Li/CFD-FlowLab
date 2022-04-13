#include "TurbSolverUnstr.h"
#include "Geo_UnstructBC.h"
#include "IO_FileName.h"
#include "MultiGridOperation.h"
#include "PHMatrix.h"
#include "TK_Exit.h"
#include "Param_TurbSolverUnstruct.h"
#include "Residual.h"
#include "FieldProxy.h"
#include "Math_BasisFunction.h"
#include "Constants.h"
#include "TK_Log.h"
#include "IO_FileName.h"
#include "Transition.h"
using namespace std;

namespace PHSPACE
{

TurbSolverUnstr::TurbSolverUnstr()
{
    DESLength = 0;
    gradientTurbField = 0;
    gradientVelocity = 0;
}

TurbSolverUnstr::~TurbSolverUnstr()
{
    DeAllocateGlobalVariables();

    delete []DESLength; 
    DESLength = 0;

    FreeControlParameters();
}

void TurbSolverUnstr::AllocateGlobalVar(Grid *gridIn)
{
    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();
    int viscoueType = parameters->GetViscousType();
    if (nTurbulenceEquation == 0) return;

    int transitionType = parameters->GetTransitionType();

    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell    = grid->GetNTotalCell();
    int numberOfTotalFace    = grid->GetNTotalFace();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal        = numberOfTotalCell + numberOfBoundaryFace;

    RDouble **qTurbulence        = NewPointer2<RDouble>(nTurbulenceEquation, numberOfTotal);
    RDouble **residualTurbulence = NewPointer2<RDouble>(nTurbulenceEquation, numberOfTotal);
    RDouble **spectrumTurbulence = NewPointer2<RDouble>(nTurbulenceEquation, numberOfTotal);

    grid->UpdateDataPtr("q_turb"   , qTurbulence      );    //! Unknown variable of turbulent model.
    grid->UpdateDataPtr("res_turb" , residualTurbulence);    //! Residual or right-hand side.
    grid->UpdateDataPtr("spec_turb", spectrumTurbulence);    //! Spectral radius of turbulent model.

    RDouble **gradientTurbulenceX = NewPointer2<RDouble>(nTurbulenceEquation, numberOfTotal);
    RDouble **gradientTurbulenceY = NewPointer2<RDouble>(nTurbulenceEquation, numberOfTotal);
    RDouble **gradientTurbulenceZ = NewPointer2<RDouble>(nTurbulenceEquation, numberOfTotal);

    grid->UpdateDataPtr("gradTurbulenceX", gradientTurbulenceX);
    grid->UpdateDataPtr("gradTurbulenceY", gradientTurbulenceY);
    grid->UpdateDataPtr("gradTurbulenceZ", gradientTurbulenceZ);

    int numberOfTotalNode = grid->GetNTotalNode();
    RDouble **qTurbulenceNode = NewPointer2<RDouble>(nTurbulenceEquation, numberOfTotalNode);
    grid->UpdateDataPtr("qTurbNode", qTurbulenceNode);    //! Unknown variable of turbulent model at node.

    int *nodeValueSliceTurb = new int[numberOfTotalNode];
    grid->UpdateDataPtr("nodeValueSliceTurb", nodeValueSliceTurb);

    int *nodeValueSliceTurbTrade = new int[numberOfTotalNode];
    PHSPACE::SetField(nodeValueSliceTurbTrade, 0, numberOfTotalNode);
    grid->UpdateDataPtr("nodeValueSliceTurbTrade", nodeValueSliceTurbTrade);

    InterpointInformation *interpointInformation = grid->GetInterpointInfo();
    if (interpointInformation)
    {
        int numberOfInterPoint = interpointInformation->GetNumberOfInterpoints();
        RDouble **qTurbulenceInterPoint = NewPointer2<RDouble>(nTurbulenceEquation, numberOfInterPoint);
        grid->UpdateDataPtr("turb::qInterpoint", qTurbulenceInterPoint);    //! Unknown variable of turbulent model at interPoint.
    }

    int isUnsteady = parameters->GetIsUnsteady();
    if (isUnsteady)
    {
        RDouble **qTurbulenceUnsteadyN1 = NewPointer2<RDouble>(nTurbulenceEquation, numberOfTotal);
        RDouble **qTurbulenceUnsteadyN2 = NewPointer2<RDouble>(nTurbulenceEquation, numberOfTotal);
        RDouble **residualTurbulenceUnsteadyN1 = NewPointer2<RDouble>(nTurbulenceEquation, numberOfTotal);
        RDouble **residualTurbulenceUnsteadyN2 = NewPointer2<RDouble>(nTurbulenceEquation, numberOfTotal);

        //! Store the historical sum of InviscidFlux, ViscousFlux, et al., exclusive of DualTimeSourceFlux, used in UpdateUnsteadyFlow for the Rn in the dualtime method.
        RDouble **residualTurbulenceUnsteadyTemporary = NewPointer2<RDouble>(nTurbulenceEquation, numberOfTotal);

        grid->UpdateDataPtr("q_turb_unsteady_n1", qTurbulenceUnsteadyN1);
        grid->UpdateDataPtr("q_turb_unsteady_n2", qTurbulenceUnsteadyN2);
        grid->UpdateDataPtr("res_turb_unsteady_n1", residualTurbulenceUnsteadyN1);
        grid->UpdateDataPtr("res_turb_unsteady_n2", residualTurbulenceUnsteadyN2);
        grid->UpdateDataPtr("res_turb_unsteady_tmp", residualTurbulenceUnsteadyTemporary);

        //! Statistical variables for DES simulation.
        int ifStaticsFlowField = parameters->GetIfStaticsFlowField();
        if (ifStaticsFlowField > 0)
        {
            RDouble **qAverageTurbulence = NewPointer2<RDouble>(nTurbulenceEquation, numberOfTotal);
            grid->UpdateDataPtr("qAverageTurb", qAverageTurbulence);
        }
    }

    if (viscoueType == 4)
    {
        RDouble *cross = new RDouble[numberOfTotal];
        grid->UpdateDataPtr("cross", cross);    //! For the crossing term in two equation turbulence model.
        RDouble *blend = new RDouble[numberOfTotal];
        grid->UpdateDataPtr("blend", blend);    //! For blending function F1.
        RDouble *SST_F2 = new RDouble[numberOfTotal];
        grid->UpdateDataPtr("SST_F2", SST_F2);    //! For blending function F2.
        if (transitionType == IREGAMA)
        {
            RDouble *SpSdRatio = new RDouble[numberOfTotal];
            grid->UpdateDataPtr("SpSdRatio", SpSdRatio);
        }
    }

    int viscousType = parameters->GetViscousType();
    int neasm    = GlobalDataBase::GetIntParaFromDB("neasm");

    if (viscousType >= TWO_EQU && neasm > 0)
    {
        int naniss = 5;
        RDouble **aniss = NewPointer2<RDouble>(naniss, numberOfTotal);

        grid->UpdateDataPtr("aniss", aniss);

        RDouble *sengy = new RDouble[numberOfTotal];

        grid->UpdateDataPtr("sengy", sengy);
    }

    RDouble **matrixTurbulenceLeft  = NewPointer2<RDouble>(nTurbulenceEquation, numberOfTotalFace);
    RDouble **matrixTurbulenceRight = NewPointer2<RDouble>(nTurbulenceEquation, numberOfTotalFace);

    grid->UpdateDataPtr("mat_turbl", matrixTurbulenceLeft);
    grid->UpdateDataPtr("mat_turbr", matrixTurbulenceRight);
}

void TurbSolverUnstr::DeAllocateGlobalVar(Grid *gridIn)
{
    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();
    int viscoueType = parameters->GetViscousType();
    if (nTurbulenceEquation == 0) return;

    int transitionType = parameters->GetTransitionType();

    UnstructGrid *grid = UnstructGridCast(gridIn);

    RDouble **qTurbulence        = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"  ));
    RDouble **residualTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));
    RDouble **spectrumTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("spec_turb"));

    DelPointer2(qTurbulence      );
    DelPointer2(residualTurbulence);
    DelPointer2(spectrumTurbulence);

    RDouble **gradientTurbulenceX = reinterpret_cast<RDouble **> (grid->GetDataPtr("gradTurbulenceX"));
    RDouble **gradientTurbulenceY = reinterpret_cast<RDouble **> (grid->GetDataPtr("gradTurbulenceY"));
    RDouble **gradientTurbulenceZ = reinterpret_cast<RDouble **> (grid->GetDataPtr("gradTurbulenceZ"));

    DelPointer2(gradientTurbulenceX);
    DelPointer2(gradientTurbulenceY);
    DelPointer2(gradientTurbulenceZ);

    int isUnsteady = parameters->GetIsUnsteady();

    if (isUnsteady)
    {
        RDouble **qTurbulenceUnsteadyN1 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb_unsteady_n1"));
        RDouble **qTurbulenceUnsteadyN2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb_unsteady_n2"));
        RDouble **residualTurbulenceUnsteadyN1 = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb_unsteady_n1"));
        RDouble **residualTurbulenceUnsteadyN2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb_unsteady_n2"));
        RDouble **residualTurbulenceUnsteadyTemporary = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb_unsteady_tmp"));

        DelPointer2(qTurbulenceUnsteadyN1);
        DelPointer2(qTurbulenceUnsteadyN2);
        DelPointer2(residualTurbulenceUnsteadyN1);
        DelPointer2(residualTurbulenceUnsteadyN2);
        DelPointer2(residualTurbulenceUnsteadyTemporary);

        //! Statistical variables for DES simulation.
        int ifStaticsFlowField = parameters->GetIfStaticsFlowField();
        if (ifStaticsFlowField > 0)
        {
            RDouble **qAverageTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("qAverageTurb"));
            DelPointer2(qAverageTurbulence);
        }
    }

    if (viscoueType == 4)
    {
        RDouble *cross  = reinterpret_cast<RDouble *> (grid->GetDataPtr("cross"));
        delete []cross;
        RDouble *blend  = reinterpret_cast<RDouble *> (grid->GetDataPtr("blend"));
        delete []blend;
        RDouble *SST_F2 = reinterpret_cast<RDouble *> (grid->GetDataPtr("SST_F2"));
        delete []SST_F2;
        if (transitionType == IREGAMA)
        {
            RDouble *SpSdRatio = reinterpret_cast<RDouble *> (grid->GetDataPtr("SpSdRatio"));
            delete [] SpSdRatio;
        }
    }

    int viscousType = parameters->GetViscousType();
    int neasm    = GlobalDataBase::GetIntParaFromDB("neasm");

    if (viscousType >= TWO_EQU && neasm > 0)
    {
        RDouble **aniss = reinterpret_cast<RDouble **> (grid->GetDataPtr("aniss"));
        DelPointer2(aniss);

        RDouble *sengy = reinterpret_cast<RDouble * > (grid->GetDataPtr("sengy"));
        delete []sengy;
    }

    RDouble **matrixTurbulenceLeft  = reinterpret_cast<RDouble **> (grid->GetDataPtr("mat_turbl"));
    RDouble **matrixTurbulenceRight = reinterpret_cast<RDouble **> (grid->GetDataPtr("mat_turbr"));

    DelPointer2(matrixTurbulenceLeft);
    DelPointer2(matrixTurbulenceRight);
}

bool TurbSolverUnstr::JudgeIfRestart()
{
    string turbfile = ".\results\turb.dat";
    GlobalDataBase::GetData("turbfile", &turbfile, PHSTRING, 1);

    if (PHMPI::IsParallelRun())
    {
        turbfile = PHSPACE::AddSymbolToFileName(turbfile, "_", 0);
    }

    bool restart_flag = false;

    ifstream infile(turbfile.c_str(), ios::in);
    if (infile)
    {
        restart_flag = true;
        infile.close();
        infile.clear();
    }

    return restart_flag;
}

bool TurbSolverUnstr::JudgeIfReadAverage()
{
    string turbfile = ".\results\turb.dat";
    GlobalDataBase::GetData("turbfile", &turbfile, PHSTRING, 1);
    string averageTurbFile = AddSymbolToFileName(turbfile, "_", "Average");

    if (PHMPI::IsParallelRun())
    {
        averageTurbFile = PHSPACE::AddSymbolToFileName(averageTurbFile, "_", 0);
    }

    bool restart_flag = false;

    ifstream infile(averageTurbFile.c_str(), ios::in);
    if (infile)
    {
        restart_flag = true;
        infile.close();
        infile.clear();
    }

    return restart_flag;
}

void TurbSolverUnstr::DumpRestart(ActionKey *actkey)
{
    int level = actkey->level;
    UnstructGrid *grid = UnstructGridCast(GetGrid(level));
    int numberOfTotalCell    = grid->GetNTotalCell();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal        = numberOfTotalCell + numberOfBoundaryFace;

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();
    int outerIterationStep = GlobalDataBase::GetIntParaFromDB("outnstep");

    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));

    DataContainer *cdata = actkey->GetData();
    cdata->moveto_begin();
    cdata->write(&outerIterationStep, sizeof(int));

    for (int m = 0; m < nTurbulenceEquation; ++ m)
    {
        cdata->write(qTurbulence[m], sizeof(RDouble) * numberOfTotal);
    }

    int isUnsteady = parameters->GetIsUnsteady();
    if (!isUnsteady) return;

    RDouble **qTurbulenceUnsteadyN1 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb_unsteady_n1"));
    RDouble **qTurbulenceUnsteadyN2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb_unsteady_n2"));

    for (int m = 0; m < nTurbulenceEquation; ++ m)
    {
        cdata->write(qTurbulenceUnsteadyN1[m], sizeof(RDouble) * numberOfTotal);
    }

    for (int m = 0; m < nTurbulenceEquation; ++ m)
    {
        cdata->write(qTurbulenceUnsteadyN2[m], sizeof(RDouble) * numberOfTotal);
    }

    RDouble **residualTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));
    RDouble **residualTurbulenceUnsteadyN1 = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb_unsteady_n1"));
    RDouble **residualTurbulenceUnsteadyN2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb_unsteady_n2"));

    for (int m = 0; m < nTurbulenceEquation; ++ m)
    {
        cdata->write(residualTurbulence[m], sizeof(RDouble) * numberOfTotal);
    }

    for (int m = 0; m < nTurbulenceEquation; ++ m)
    {
        cdata->write(residualTurbulenceUnsteadyN1[m], sizeof(RDouble) * numberOfTotal);
    }

    for (int m = 0; m < nTurbulenceEquation; ++ m)
    {
        cdata->write(residualTurbulenceUnsteadyN2[m], sizeof(RDouble) * numberOfTotal);
    }
}

void TurbSolverUnstr::DumpRestartH5(ActionKey *actkey)
{
    using namespace PHMPI;
    int currentProcessorID = GetCurrentProcessorID();

    Grid *grid = GetGrid(actkey->level);
    int gridID = grid->GetZoneID();

    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotalCell    = grid->GetNTotalCell();
    int numberOfTotal        = numberOfTotalCell + numberOfBoundaryFace;

    int version = 1;

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();
    int outerIterationStep = GlobalDataBase::GetIntParaFromDB("outnstep");
    if (currentProcessorID == GetServerProcessorID())
    {
        WriteData(actkey->filepos, &version, "Version");
        WriteData(actkey->filepos, &outerIterationStep, "outnstep");

        if (IsNeedStatistics())
        {
            int nStatisticalStep = 0;
            GlobalDataBase::GetData("nStatisticalStep", &nStatisticalStep, PHINT, 1);
            WriteData(actkey->filepos, &nStatisticalStep, "nStatisticalStep");
        }
    }

    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
    RDouble **qTemporary  = NewPointer2<RDouble>(nTurbulenceEquation, numberOfTotal);

    for (int iEquation = 0; iEquation < nTurbulenceEquation; ++ iEquation)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            qTemporary[iEquation][iCell] = qTurbulence[iEquation][iCell];
        }
    }

    hid_t grploc;
    string grpName;

    ostringstream oss;
    oss<<"Group"<<gridID;
    grpName = oss.str();

    grploc = OpenGroup(actkey->filepos, grpName);

    WriteData(grploc, qTemporary[0], "q_turb");
    WriteData(grploc, &numberOfTotalCell, "nTotalCell");

    int isUnsteady = parameters->GetIsUnsteady();
    if (!isUnsteady)
    {
        H5Gclose(grploc);
        DelPointer2(qTemporary);

        return;
    }

    RDouble **qTurbulenceUnsteadyN1 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb_unsteady_n1"));
    for (int iEquation = 0; iEquation < nTurbulenceEquation; ++ iEquation)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            qTemporary[iEquation][iCell] = qTurbulenceUnsteadyN1[iEquation][iCell];
        }
    }

    WriteData(grploc, qTemporary[0], "q_turb_unsteady_n1");

    RDouble **qTurbulenceUnsteadyN2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb_unsteady_n2"));
    for (int iEquation = 0; iEquation < nTurbulenceEquation; ++ iEquation)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            qTemporary[iEquation][iCell] = qTurbulenceUnsteadyN2[iEquation][iCell];
        }
    }

    WriteData(grploc, qTemporary[0], "q_turb_unsteady_n2");

    RDouble **residualTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));
    for (int iEquation = 0; iEquation < nTurbulenceEquation; ++ iEquation)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            qTemporary[iEquation][iCell] = residualTurbulence[iEquation][iCell];
        }
    }

    WriteData(grploc, qTemporary[0], "res_turb");

    RDouble **residualTurbulenceUnsteadyN1 = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb_unsteady_n1"));
    for (int iEquation = 0; iEquation < nTurbulenceEquation; ++ iEquation)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            qTemporary[iEquation][iCell] = residualTurbulenceUnsteadyN1[iEquation][iCell];
        }
    }

    WriteData(grploc, qTemporary[0], "res_turb_unsteady_n1");

    RDouble **residualTurbulenceUnsteadyN2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb_unsteady_n2"));
    for (int iEquation = 0; iEquation < nTurbulenceEquation; ++ iEquation)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            qTemporary[iEquation][iCell] = residualTurbulenceUnsteadyN2[iEquation][iCell];
        }
    }

    WriteData(grploc, qTemporary[0], "res_turb_unsteady_n2");

    if (IsNeedStatistics())
    {
        RDouble **qAverageTurb = reinterpret_cast<RDouble **> (grid->GetDataPtr("qAverageTurb"));
        for (int iEquation = 0; iEquation < nTurbulenceEquation; ++ iEquation)
        {
            for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
            {
                qTemporary[iEquation][iCell] = qAverageTurb[iEquation][iCell];
            }
        }

        WriteData(grploc, qTemporary[0], "qAverageTurb");
    }

    H5Gclose(grploc);
    DelPointer2(qTemporary);
}

void TurbSolverUnstr::ReadRestartH5(ActionKey *actkey)
{
    Grid *grid = GetGrid(actkey->level);
    int gridID = grid->GetZoneID();

    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotalCell    = grid->GetNTotalCell();
    int numberOfTotal        = numberOfTotalCell + numberOfBoundaryFace;

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
    RDouble **qTemporary  = NewPointer2<RDouble>(nTurbulenceEquation, numberOfTotal);

    int outerIterationStepOfTurbulence = 0;
    ReadData(actkey->filepos, &outerIterationStepOfTurbulence, "outnstep");

    int outerIterationStepOfNS = GlobalDataBase::GetIntParaFromDB("outnstep");

    hid_t grploc;
    string grpName;

    ostringstream oss;
    oss<<"Group"<<gridID;
    grpName = oss.str();

    grploc = OpenGroup(actkey->filepos, grpName);

    int nTotalCellRestart = 0;
    ReadData(grploc, &nTotalCellRestart, "nTotalCell");
    if(nTotalCellRestart != numberOfTotalCell)
    {
        ostringstream erroeInfo;
        erroeInfo << " Error: the cell number in turb.dat is not equal to the cell number in grid file !" << endl;
        TK_Exit::ExceptionExit(erroeInfo.str());
    }

    ReadData(grploc, qTemporary[0], "q_turb");
    for (int iEquation = 0; iEquation < nTurbulenceEquation; ++ iEquation)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            qTurbulence[iEquation][iCell] = qTemporary[iEquation][iCell];
        }
    }

    int isUnsteady = parameters->GetIsUnsteady();
    if (isUnsteady)
    {
        RDouble **qTurbulenceUnsteadyN1 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb_unsteady_n1"));
        RDouble **qTurbulenceUnsteadyN2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb_unsteady_n2"));

        RDouble **residualTurbulence           = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));
        RDouble **residualTurbulenceUnsteadyN1 = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb_unsteady_n1"));
        RDouble **residualTurbulenceUnsteadyN2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb_unsteady_n2"));
        RDouble **qAverageTurb                 = reinterpret_cast<RDouble **> (grid->GetDataPtr("qAverageTurb"));

        //! Start from steady flow field and reset the outer step when start from steady flow.
        int ifStartFromSteadyResults = parameters->GetIfStartFromSteadyResults();
        if (ifStartFromSteadyResults)
        {
            PrintToWindow("Restart from steady Turbulent flow field, reset outer step to be zero!\n");

            outerIterationStepOfTurbulence = 0;
            GlobalDataBase::UpdateData("outnstep", &outerIterationStepOfTurbulence, PHINT, 1);

            for (int m = 0; m < nTurbulenceEquation; ++ m)
            {
                for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
                {
                    qTurbulenceUnsteadyN1[m][iCell] = qTurbulence[m][iCell];
                    qTurbulenceUnsteadyN2[m][iCell] = qTurbulence[m][iCell];

                    residualTurbulence[m][iCell]           = 0.0;
                    residualTurbulenceUnsteadyN1[m][iCell] = 0.0;
                    residualTurbulenceUnsteadyN2[m][iCell] = 0.0;
                }
            }
        }
        else
        {
            ReadData(grploc, qTemporary[0], "q_turb_unsteady_n1");
            for (int iEquation = 0; iEquation < nTurbulenceEquation; ++ iEquation)
            {
                for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
                {
                    qTurbulenceUnsteadyN1[iEquation][iCell] = qTemporary[iEquation][iCell];
                }
            }

            ReadData(grploc, qTemporary[0], "q_turb_unsteady_n2");
            for (int iEquation = 0; iEquation < nTurbulenceEquation; ++ iEquation)
            {
                for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
                {
                    qTurbulenceUnsteadyN2[iEquation][iCell] = qTemporary[iEquation][iCell];
                }
            }

            ReadData(grploc, qTemporary[0], "res_turb");
            for (int iEquation = 0; iEquation < nTurbulenceEquation; ++ iEquation)
            {
                for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
                {
                    residualTurbulence[iEquation][iCell] = qTemporary[iEquation][iCell];
                }
            }

            ReadData(grploc, qTemporary[0], "res_turb_unsteady_n1");
            for (int iEquation = 0; iEquation < nTurbulenceEquation; ++ iEquation)
            {
                for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
                {
                    residualTurbulenceUnsteadyN1[iEquation][iCell] = qTemporary[iEquation][iCell];
                }
            }

            ReadData(grploc, qTemporary[0], "res_turb_unsteady_n2");
            for (int iEquation = 0; iEquation < nTurbulenceEquation; ++ iEquation)
            {
                for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
                {
                    residualTurbulenceUnsteadyN2[iEquation][iCell] = qTemporary[iEquation][iCell];
                }
            }
        }

        int nStatisticalStep = 0;
        int ifStaticsFlowField = parameters->GetIfStaticsFlowField();
        int startStatisticStep = GlobalDataBase::GetIntParaFromDB("startStatisticStep");
        bool isReadAverageFlow = ifStaticsFlowField && (outerIterationStepOfTurbulence >= startStatisticStep);
        if(isReadAverageFlow && (!ifStartFromSteadyResults))
        {
            ReadData(grploc, qTemporary[0], "qAverageTurb");
            for (int iEquation = 0; iEquation < nTurbulenceEquation; ++ iEquation)
            {
                for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
                {
                    qAverageTurb[iEquation][iCell] = qTemporary[iEquation][iCell];
                }
            }

            ReadData(actkey->filepos, &nStatisticalStep, "nStatisticalStep");
            GlobalDataBase::UpdateData("nStatisticalStep", &nStatisticalStep, PHINT, 1);
        }
        else
        {
            for (int iEquation = 0; iEquation < nTurbulenceEquation; ++ iEquation)
            {
                for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
                {
                    qAverageTurb[iEquation][iCell] = qTurbulence[iEquation][iCell];
                }
            }

            GlobalDataBase::UpdateData("nStatisticalStep", &nStatisticalStep, PHINT, 1);
        }
    }

    H5Gclose(grploc);
    DelPointer2(qTemporary);

    CompareOutStepOfFlowfieldFile(outerIterationStepOfNS, outerIterationStepOfTurbulence);
}

void TurbSolverUnstr::ReadRestart(ActionKey *actkey)
{
    int level = actkey->level;
    UnstructGrid *grid = UnstructGridCast(GetGrid(level));

    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotalCell    = grid->GetNTotalCell();
    int numberOfTotal        = numberOfTotalCell + numberOfBoundaryFace;

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));

    DataContainer *cdata = actkey->GetData();
    cdata->moveto_begin();
    int outerIterationStepOfTurbulence = 0;
    cdata->read(&outerIterationStepOfTurbulence, sizeof(int));

    int outerIterationStepOfNS = GlobalDataBase::GetIntParaFromDB("outnstep");

    for (int m = 0; m < nTurbulenceEquation; ++ m)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            cdata->read(&qTurbulence[m][iCell], sizeof(RDouble));
        }
    }

    int isUnsteady = parameters->GetIsUnsteady();
    if (isUnsteady)
    {
        RDouble **qTurbulenceUnsteadyN1 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb_unsteady_n1"));
        RDouble **qTurbulenceUnsteadyN2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb_unsteady_n2"));

        RDouble **residualTurbulence           = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));
        RDouble **residualTurbulenceUnsteadyN1 = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb_unsteady_n1"));
        RDouble **residualTurbulenceUnsteadyN2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb_unsteady_n2"));

        int ifStartFromSteadyResults = parameters->GetIfStartFromSteadyResults();

        if (ifStartFromSteadyResults)
        {
            PrintToWindow("Restart from steady Turbulent flow field, reset outer step to be zero!\n");

            //! Start from steady flow field and reset the outer step when start from steady flow.
            outerIterationStepOfTurbulence = 0;
            GlobalDataBase::UpdateData("outnstep", &outerIterationStepOfTurbulence, PHINT, 1);

            for (int m = 0; m < nTurbulenceEquation; ++ m)
            {
                for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
                {
                    qTurbulenceUnsteadyN1[m][iCell] = qTurbulence[m][iCell];
                    qTurbulenceUnsteadyN2[m][iCell] = qTurbulence[m][iCell];

                    residualTurbulence          [m][iCell] = 0.0;
                    residualTurbulenceUnsteadyN1[m][iCell] = 0.0;
                    residualTurbulenceUnsteadyN2[m][iCell] = 0.0;
                }
            }
        }
        else
        {
            //! Start from unsteady flow field.
            for (int m = 0; m < nTurbulenceEquation; ++ m)
            {
                cdata->read(qTurbulenceUnsteadyN1[m], sizeof(RDouble) * numberOfTotal);
            }
            for (int m = 0; m < nTurbulenceEquation; ++ m)
            {
                cdata->read(qTurbulenceUnsteadyN2[m], sizeof(RDouble) * numberOfTotal);
            }

            for (int m = 0; m < nTurbulenceEquation; ++ m)
            {
                cdata->read(residualTurbulence[m], sizeof(RDouble) * numberOfTotal);
            }
            for (int m = 0; m < nTurbulenceEquation; ++ m)
            {
                cdata->read(residualTurbulenceUnsteadyN1[m], sizeof(RDouble) * numberOfTotal);
            }
            for (int m = 0; m < nTurbulenceEquation; ++ m)
            {
                cdata->read(residualTurbulenceUnsteadyN2[m], sizeof(RDouble) * numberOfTotal);
            }
        }
    }

    CompareOutStepOfFlowfieldFile(outerIterationStepOfNS, outerIterationStepOfTurbulence);
}

void TurbSolverUnstr::InitFlowAsRestart()
{
    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();
    if (nTurbulenceEquation == 0) return;

    UnstructGrid *grid = UnstructGridCast(GetGrid(0));

    //int outerIterationStep = 0;
    //GlobalDataBase::UpdateData("outnstep", &outerIterationStep, PHINT, 1);

    int numberOfTotalCell    = grid->GetNTotalCell();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal        = numberOfTotalCell + numberOfBoundaryFace;

    RDouble *freeStreamTurbVar =  parameters->GetFreeStreamTurbVar();
    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));

    for (int m = 0; m < nTurbulenceEquation; ++ m)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            qTurbulence[m][iCell] = freeStreamTurbVar[m];
        }
    }

    RDouble *SpSdRatio = reinterpret_cast<RDouble *> (grid->GetDataPtr("SpSdRatio"));
    int transitionType = parameters->GetTransitionType();
    if (transitionType == IREGAMA)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            SpSdRatio[iCell] = 1.0;
        }
    }

    int isUnsteady = parameters->GetIsUnsteady();
    if (!isUnsteady) return;

    RDouble **qTurbulenceUnsteadyN1 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb_unsteady_n1"));
    RDouble **qTurbulenceUnsteadyN2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb_unsteady_n2"));

    for (int m = 0; m < nTurbulenceEquation; ++ m)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            qTurbulenceUnsteadyN1[m][iCell] = qTurbulence[m][iCell];
        }
    }

    for (int m = 0; m < nTurbulenceEquation; ++ m)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            qTurbulenceUnsteadyN2[m][iCell] = qTurbulence[m][iCell];
        }
    }

    RDouble **residualTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));
    RDouble **residualTurbulenceUnsteadyN1 = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb_unsteady_n1"));
    RDouble **residualTurbulenceUnsteadyN2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb_unsteady_n2"));
    RDouble **residualTurbulenceUnsteadyTemporary = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb_unsteady_tmp"));

    for (int m = 0; m < nTurbulenceEquation; ++ m)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            residualTurbulence[m][iCell] = 0.0;
            residualTurbulenceUnsteadyN1[m][iCell] = 0.0;
            residualTurbulenceUnsteadyN2[m][iCell] = 0.0;
            residualTurbulenceUnsteadyTemporary[m][iCell] = 0.0;
        }
    }
}

void TurbSolverUnstr::ZeroResiduals(Grid *gridIn)
{
    UnstructGrid *grid = reinterpret_cast<UnstructGrid *>(gridIn);

    int numberOfGeneralCells = grid->GetNTotalCell() + grid->GetNBoundFace();

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    RDouble **residualTurbulence = reinterpret_cast<RDouble **>(grid->GetDataPtr("res_turb"));

    for (int iEquation = 0; iEquation < nTurbulenceEquation; ++ iEquation)
    {
        for (int iCell = 0; iCell < numberOfGeneralCells; ++ iCell)
        {
            residualTurbulence[iEquation][iCell] = 0.0;
        }
    }
    return;
}

void TurbSolverUnstr::InitSpectrum(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotalCell    = grid->GetNTotalCell();
    int numberOfTotal        = numberOfTotalCell + numberOfBoundaryFace;

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();
    if (nTurbulenceEquation == 0) return;

    RDouble turbCFLScale = parameters->GetTurbCFLScale();

    RDouble *vol  = grid->GetCellVolume();
    RDouble *dt = reinterpret_cast<RDouble *> (grid->GetDataPtr("dt"));
    RDouble **spectrumTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("spec_turb"));

    RDouble dualTimeSpectrumC1 = zero;
    RDouble dualTimeSpectrumC2 = one;

    //! If flow is unsteady, need to count the contribution of unsteady.
    int isUnsteady = parameters->GetIsUnsteady();
    if (isUnsteady)
    {
        RDouble dualTimeCoefficient[7];
        const int methodOfDualTime = GlobalDataBase::GetIntParaFromDB("methodOfDualTime");
        ComputeDualTimeCoefficient(methodOfDualTime, dualTimeCoefficient);
        dualTimeSpectrumC1 = - dualTimeCoefficient[3];
        dualTimeSpectrumC2 =   dualTimeCoefficient[6];
    }

    for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
    {
        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            spectrumTurbulence[m][iCell] = 0.0;
        }
    }

    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
    {
        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            spectrumTurbulence[m][iCell] = dualTimeSpectrumC1 * vol[iCell] + dualTimeSpectrumC2 / (turbCFLScale * dt[iCell] + SMALL);
        }
    }
}

void TurbSolverUnstr::CorrectFineGrid(Grid *fineGrid, Grid *coarseGrid)
{
    int numberOfTotalCellOnFineGrid = UnstructGridCast(fineGrid)->GetNTotalCell();
    int *cell2CoarseGridCell = UnstructGridCast(fineGrid)->GetCell2CoarseGridCell();

    RDouble **qTurbulenceFineGrid = reinterpret_cast<RDouble **> (fineGrid->GetDataPtr("q_turb"));
    RDouble **qTurbulenceCoarseGrid = reinterpret_cast<RDouble **> (coarseGrid->GetDataPtr("q_turb"));

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    for (int m = 0; m < nTurbulenceEquation; ++ m)
    {
        for (int iCell = 0; iCell <numberOfTotalCellOnFineGrid; ++ iCell)
        {
            qTurbulenceFineGrid[m][iCell] += qTurbulenceCoarseGrid[m][cell2CoarseGridCell[iCell]];
        }
    }
}

void TurbSolverUnstr::InterpolatFineGrid(Grid *fineGridIn, Grid *coarseGridIn)
{
    UnstructGrid *fineGrid = UnstructGridCast(fineGridIn);
    UnstructGrid *coarseGrid = UnstructGridCast(coarseGridIn);

    RDouble **qTurbulenceFineGrid = reinterpret_cast<RDouble **> (fineGrid->GetDataPtr("q_turb"));
    RDouble **qTurbulenceCoarseGrid = reinterpret_cast<RDouble **> (coarseGrid->GetDataPtr("q_turb"));

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    for (int m = 0; m < nTurbulenceEquation; ++ m)
    {
        InterpolatQ(fineGrid, qTurbulenceFineGrid[m], coarseGrid, qTurbulenceCoarseGrid[m]);
    }

    RDouble *viscousOfFineGrid = reinterpret_cast<RDouble *> (fineGrid->GetDataPtr("vist"));
    RDouble *viscousOfCoarseGrid = reinterpret_cast<RDouble *> (coarseGrid->GetDataPtr("vist"));
    InterpolatQ(fineGrid, viscousOfFineGrid, coarseGrid, viscousOfCoarseGrid);
}

void TurbSolverUnstr::PutCorrectionBack(Grid *gridIn, FieldProxy *qProxy)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell    = grid->GetNTotalCell();

    RDouble **q = qProxy->GetField_UNS();
    RDouble **qOld = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
    {
        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            qOld[m][iCell] += q[m][iCell];
        }
    }
}

void TurbSolverUnstr::RestrictDefect(Grid *fineGridIn, Grid *coarseGridIn)
{
    UnstructGrid *fineGrid = UnstructGridCast(fineGridIn);
    UnstructGrid *coarseGrid = UnstructGridCast(coarseGridIn);

    int numberOfTotalCellOnFineGrid = fineGrid->GetNTotalCell();
    int *cell2CoarseGridCell = fineGrid->GetCell2CoarseGridCell();

    ZeroResiduals(coarseGrid);

    RDouble **residualTuebulenceCoarseGrid = reinterpret_cast<RDouble **> (coarseGrid->GetDataPtr("res_turb"));
    RDouble **residualTuebulenceFineGrid = reinterpret_cast<RDouble **> (fineGrid->GetDataPtr("res_turb"));

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    for (int m = 0; m < nTurbulenceEquation; ++ m)
    {
        for (int iCell = 0; iCell <numberOfTotalCellOnFineGrid; ++ iCell)
        {
            residualTuebulenceCoarseGrid[m][cell2CoarseGridCell[iCell]] -= residualTuebulenceFineGrid[m][iCell];
        }
    }
}

void TurbSolverUnstr::LoadQ(Grid *gridIn, FieldProxy *qProxy)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell    = grid->GetNTotalCell();

    RDouble **qqTurbulence = qProxy->GetField_UNS();

    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
    {
        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            qqTurbulence[m][iCell] = qTurbulence[m][iCell];
        }
    }
}

void TurbSolverUnstr::UpdateFlowField(Grid *gridIn, FieldProxy *qProxy, FieldProxy *dqProxy)
{
    UnstructGrid *grid  = UnstructGridCast(gridIn);
    int numberOfTotalCell = grid->GetNTotalCell();

    RDouble **qPrimitiveVariable = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble **qTurbulence  = qProxy->GetField_UNS();
    RDouble **dqTurbulence = dqProxy->GetField_UNS();

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();
    int viscoueType = parameters->GetViscousType();

    using namespace IDX;

    RDouble *primitiveVariable = new RDouble [nTurbulenceEquation];

    if (viscoueType == 3)
    {
        //! Spalart-allmaras model.
        int nNegative = 0;

        RDouble nuoo  = 1.0;        //!Nondimensional inflow laminar viscosity / rho.
        RDouble nuMin = 1.0e-5 * nuoo;
        RDouble nuMax = 1.0e10 * nuoo;
        int newNstep = GlobalDataBase::GetIntParaFromDB("newnstep");
        if(newNstep < 1000)
        {
            nuMax = 1.0e5 * nuoo;
        }

        for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
        {
            if(fabs(dqTurbulence[ISA][iCell] / qTurbulence[ISA][iCell]) > half)
            {
                dqTurbulence[ISA][iCell] *= half * fabs(qTurbulence[ISA][iCell]) / fabs(dqTurbulence[ISA][iCell]);
            }

            qTurbulence[ISA][iCell] += dqTurbulence[ISA][iCell]; 

            if (qTurbulence[ISA][iCell] < 0.0)
            {
                RDouble *xcc = grid->GetCellCenterX();
                RDouble *ycc = grid->GetCellCenterY();
                RDouble *zcc = grid->GetCellCenterZ();
                vector<int> * c2c = grid->GetCell2Cell();
                RDouble *wallDist = this->GetLengthScale(gridIn);

                ++ nNegative;

                cout.setf(ios::scientific);
                cout.precision(8);

                int outerIterationStep = GlobalDataBase::GetIntParaFromDB("outnstep");
                cout << "    iter = " << outerIterationStep << endl;
                cout << "    Negative turbulence at " << " Zone = " << grid->GetZoneID()
                    << ", cell = " << iCell << ", " << " nTCell = " << numberOfTotalCell 
                    << ", qTurbulence = " << qTurbulence[ISA][iCell] << ", dq = " << dqTurbulence[ISA][iCell] << "\n";

                cout << "    Wall distance = " << wallDist[iCell] << ", cell center: " << xcc[iCell] << ", " << ycc[iCell] << ", " << zcc[iCell] << endl;
                cout << "        Neighbor cell: " << endl;

                int nNeighbor = c2c[iCell].size();
                for (int iNeighbor = 0; iNeighbor < nNeighbor; ++ iNeighbor)
                {
                    int neighborCell = c2c[iCell][iNeighbor];

                    cout << "        neighbor cell, " << neighborCell
                         << ", walldist = " << wallDist[neighborCell]
                         << ", center: " << xcc[neighborCell] << ", " << ycc[neighborCell] << ", " << zcc[neighborCell] << endl;
                }
            }

            qTurbulence[ISA][iCell] = MAX(qTurbulence[ISA][iCell], nuMin);
            qTurbulence[ISA][iCell] = MIN(qTurbulence[ISA][iCell], nuMax);
        }
    }
    else if (viscoueType == 4)
    {
        RDouble *freeStreamTurbVar = parameters->GetFreeStreamTurbVar();

        RDouble kelim = 1.0e-5 * freeStreamTurbVar[IKE];
        RDouble kwlim = 1.0e-5 * freeStreamTurbVar[IKW];

        int transitionType = parameters->GetTransitionType();

        for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
        {
            RDouble ke,kw,dke,dkw,rho;
            rho = qPrimitiveVariable[IR][iCell];
            ke  = qTurbulence[IKE][iCell];
            kw  = qTurbulence[IKW][iCell]; 
            dke = dqTurbulence[IKE][iCell];
            dkw = dqTurbulence[IKW][iCell];
            ke += dke/rho;
            kw += dkw/rho;
            ke = MAX(ke,kelim);
            kw = MAX(kw,kwlim);
            qTurbulence[IKE][iCell] = ke;
            qTurbulence[IKW][iCell] = kw;

            if (transitionType == IREGAMA)
            {
                RDouble &kgama = qTurbulence[IGAMA][iCell];
                RDouble &krect = qTurbulence[IRECT][iCell];
            
                RDouble dkgama = dqTurbulence[IGAMA][iCell];
                RDouble dkrect = dqTurbulence[IRECT][iCell];
            
                kgama += dkgama / rho;
                krect += dkrect / rho;
            
                kgama = MIN(kgama, 1.00);
                krect = MAX(krect, 20.0);
            }
        }
    }

    delete [] primitiveVariable;
}

void TurbSolverUnstr::SMoothTurbulencePoint(Grid *gridIn, RDouble *primitiveVariable, int i)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell    = grid->GetNTotalCell();

    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    int **cell2Face = grid->GetCell2Face();
    int *faceNumberOfEachCell = grid->GetFaceNumberOfEachCell();

    int le,re,ie,face;
    RDouble volumeSum = zero;
    for (int m = 0; m < nTurbulenceEquation; ++ m)
    {
        primitiveVariable[m] = zero;
    }

    for (int j = 0; j < faceNumberOfEachCell[i]; ++ j)
    {
        face = cell2Face[i][j];
        le   = leftCellOfFace[face];
        re   = rightCellOfFace[face];

        if (le != i)
        {
            ie = le;
        }
        else
        {
            ie = re;
        }

        if (ie >= numberOfTotalCell) continue;

        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            primitiveVariable[m] += qTurbulence[m][ie];
        }
        volumeSum += 1.0;
    }

    for (int m = 0; m < nTurbulenceEquation; ++ m)
    {
        primitiveVariable[m] /= volumeSum;
    }
}

void TurbSolverUnstr::SMoothTurbulence(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell    = grid->GetNTotalCell();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal        = numberOfTotalCell + numberOfBoundaryFace;

    int *leftCellOfFace  = grid->GetLeftCellOfFace();
    int *rightCellOfFace  = grid->GetRightCellOfFace();

    RDouble *volume = grid->GetCellVolume() ;

    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    RDouble **tmporary = NewPointer2<RDouble>(nTurbulenceEquation, numberOfTotal);

    int **cell2Face = grid->GetCell2Face();
    int *faceNumberOfEachCell = grid->GetFaceNumberOfEachCell();

    for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
    {
        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            qTurbulence[m][iCell] *= volume[iCell];
            tmporary[m][iCell]     = qTurbulence[m][iCell];
        }
    }

    RDouble volumeSum;

    int le,re,ie,face;
    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
    {
        volumeSum = volume[iCell];
        for (int j = 0; j < faceNumberOfEachCell[iCell]; ++ j)
        {
            face = cell2Face[iCell][j];
            le    = leftCellOfFace [face];
            re    = rightCellOfFace[face];
            
            if (le != iCell)
            {
                ie = le;
            }
            else
            {
                ie = re;
            }

            if (ie >= numberOfTotalCell) continue;

            for (int m = 0; m < nTurbulenceEquation; ++ m)
            {
                qTurbulence[m][iCell] += tmporary[m][ie];
            }
            volumeSum += volume[ie];
        }

        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            qTurbulence[m][iCell] /= volumeSum;
        }
    }

    DelPointer2(tmporary);
}

void TurbSolverUnstr::RungeKuttaResidual(Grid *gridIn, FieldProxy *dqProxy, RDouble coefficient)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell = grid->GetNTotalCell();

    RDouble **dqTurbulence = dqProxy->GetField_UNS();

    RDouble *dt = reinterpret_cast<RDouble *> (grid->GetDataPtr("dt"));

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
    {
        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            dqTurbulence[m][iCell] = dt[iCell] * coefficient * dqTurbulence[m][iCell];
        }
    }
}

void TurbSolverUnstr::SolveLUSGSForward(Grid *gridIn, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell = grid->GetNTotalCell();

    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    RDouble **dq = dqProxy->GetField_UNS();
    int *iBlank = grid->GetBlankIndex();

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    RDouble *df    = new RDouble [nTurbulenceEquation];
    RDouble *rhs0  = new RDouble [nTurbulenceEquation];
    RDouble *dqOld = new RDouble [nTurbulenceEquation];

    RDouble **spectrumTurbulence    = reinterpret_cast<RDouble **> (grid->GetDataPtr("spec_turb"));
    RDouble **matrixTurbulenceLeft  = reinterpret_cast<RDouble **> (grid->GetDataPtr("mat_turbl"));
    RDouble **matrixTurbulenceRight = reinterpret_cast<RDouble **> (grid->GetDataPtr("mat_turbr"));
    RDouble **residualTurbulence    = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));

    int **cell2Face = grid->GetCell2Face();
    int *faceNumberOfEachCell = grid->GetFaceNumberOfEachCell();

    RDouble **deltaFlux = LUplusDQ->GetField_UNS();
    RDouble *matrixUx = new RDouble [nTurbulenceEquation];

    //! Now the Forward Sweep
    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
    {
        if (iBlank[iCell] <= 0)
        {
            for (int m = 0; m < nTurbulenceEquation; ++ m)
            {
                dq[m][iCell] = 0.0;
            }
            continue;
        }

        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            rhs0[m] = 0.0;
        }

        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            //! Back up the old dq, to compute the convergency.
            dqOld[m] = dq[m][iCell];

            //! Here, the deltaFlux is computed in Upper backward sweep!
            //! res: b      deltaFlux: Ux, which is computed in backward sweep.
            //! the dq is not the real dq, now.
            //! Although the 'dq' changed can not the right one, it dosen't matter, since 
            //! the following only using the lower neighbor cells, whose 'dq' has been updated.
            //! It is convenient for precondition transform.
            dq[m][iCell] = residualTurbulence[m][iCell];
            matrixUx[m]  = deltaFlux[m][iCell];

            //! Then reset it to zero to store the Lower forward sweep!
            deltaFlux[m][iCell] = 0.0;
        }

        for (int jFace = 0; jFace < faceNumberOfEachCell[iCell]; ++ jFace)
        {
            int face = cell2Face[iCell][jFace];
            int le   = leftCellOfFace [face];
            int re   = rightCellOfFace [face];
            
            //! One of le and re must be cell itself.
            if (le > iCell || re > iCell) continue;            

            if (re == iCell)
            {
                //! Now its neighboring cell belongs to lower triangular
                SWAP(le, re);

                for (int m = 0; m < nTurbulenceEquation; ++ m)
                {
                    df[m] = matrixTurbulenceLeft[m][face] * dq[m][re];
                }
            }
            else
            {
                for (int m = 0; m < nTurbulenceEquation; ++ m)
                {
                    df[m] = matrixTurbulenceRight[m][face] * dq[m][re];
                }
            }

            //! Add Flux together
            for (int m = 0; m < nTurbulenceEquation; ++ m)
            {
                rhs0[m] += df[m];
            }
        }

        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            //! dq = { (b - Ux) - Lx } / D.
            //! Note: the 'dq' has actually initialized by the residual at the beginning of this function.
            //! the 'rhs0' is the total Delta Flux of the neighbors.            
            //! rhs0: Lx    diagonal: D.
            dq[m][iCell] = (dq[m][iCell] - matrixUx[m] - rhs0[m]) / spectrumTurbulence[m][iCell];

            //! Store the lower forward sweep delta-flux, which will be used in the backward sweep.
            deltaFlux[m][iCell] += rhs0[m];

            sweepNormal += SQR(dq[m][iCell] - dqOld[m]);
        }
    }

    delete []rhs0;
    delete []df;
    delete []dqOld;
    delete []matrixUx;
}

void TurbSolverUnstr::SolveLUSGSBackward(Grid *gridIn, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell = grid->GetNTotalCell();

    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    RDouble **dq = dqProxy->GetField_UNS();
    int *iBlank = grid->GetBlankIndex();

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    RDouble *df    = new RDouble [nTurbulenceEquation];
    RDouble *rhs0  = new RDouble [nTurbulenceEquation];
    RDouble *dqOld = new RDouble [nTurbulenceEquation];

    RDouble **spectrumTurbulence    = reinterpret_cast<RDouble **> (grid->GetDataPtr("spec_turb"));
    RDouble **matrixTurbulenceLeft  = reinterpret_cast<RDouble **> (grid->GetDataPtr("mat_turbl"));
    RDouble **matrixTurbulenceRight = reinterpret_cast<RDouble **> (grid->GetDataPtr("mat_turbr"));
    RDouble **residualTurbulence    = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));

    int **cell2Face = grid->GetCell2Face();
    int *faceNumberOfEachCell = grid->GetFaceNumberOfEachCell();

    RDouble **deltaFlux = LUplusDQ->GetField_UNS();
    RDouble *matrixLx = new RDouble [nTurbulenceEquation];

    // Backward Sweep
    for (int iCell = numberOfTotalCell - 1; iCell >= 0; -- iCell)
    {
        if (iBlank[iCell] <= 0)
        {
            for (int m = 0; m < nTurbulenceEquation; ++ m)
            {
                dq[m][iCell] = 0.0;
            }
            continue;
        }

        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            rhs0[m] = 0.0;
        }

        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            //! Back up the old dq, to compute the convergence.
            //! the 'dq' is dq*, which has been updated in forward.
            dqOld[m] = dq[m][iCell];

            //! the dq is not the real dq, now.
            //! it is convenient for precondition transform.
            dq[m][iCell] = residualTurbulence[m][iCell];

            matrixLx[m] = deltaFlux[m][iCell];

            deltaFlux[m][iCell] = 0.0;
        }

        for (int jFace = 0; jFace < faceNumberOfEachCell[iCell]; ++ jFace)
        {
            int face  = cell2Face[iCell][jFace];
            int le    = leftCellOfFace  [face];
            int re    = rightCellOfFace [face];

            //! One of le and re must be cell itself.
            if (le < iCell || re < iCell) continue;
            
            if (re == iCell)
            {
                //! Now its neighboring cell belongs to upper triangular
                SWAP(le, re);

                for (int m = 0; m < nTurbulenceEquation; ++ m)
                {
                    df[m] = matrixTurbulenceLeft[m][face] * dq[m][re];
                }
            }
            else
            {
                for (int m = 0; m < nTurbulenceEquation; ++ m)
                {
                    df[m] = matrixTurbulenceRight[m][face] * dq[m][re];
                }
            }

            for (int m = 0; m < nTurbulenceEquation; ++ m)
            {
                rhs0[m] += df[m];
            }
        }

        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            //! Note: the 'dq' has been updated by the forward sweep.
            //! the 'rhs0' is the total Delta Flux of the neighbors.
            //! x = {(b - LX) - Ux} / D.
            //! rhs0: Ux.    diagonal: D.
            dq[m][iCell] = (dq[m][iCell] - matrixLx[m] - rhs0[m]) / spectrumTurbulence[m][iCell];

            //! Store the upper backward sweep delta-flux, which will be used in the forward sweep.
            deltaFlux[m][iCell] += rhs0[m];

            sweepNormal += SQR(dq[m][iCell] - dqOld[m]);
        }
    }

    delete []rhs0;
    delete []df;
    delete []dqOld;
    delete []matrixLx;
}

// Bell 20130401 add
void TurbSolverUnstr::SetGhostDQLUSGS(Grid *gridIn, RDouble **dq)
{
    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();
    if (nTurbulenceEquation != 1)
    {
        cout <<"SetGhostDQLUSGS function is turned off !\n";
        return;
    }

    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *leftCellOfFace  = grid->GetLeftCellOfFace();
    int *rightCellOfFace  = grid->GetRightCellOfFace();

    CompositeUnstructBC *unstructBCRegion = grid->GetUnstructBCRegion();
    int nBCRegionUnstruct = unstructBCRegion->GetnBCRegion();

    for (int iBCRegionUnstruct = 0; iBCRegionUnstruct < nBCRegionUnstruct; iBCRegionUnstruct++)
    {
        UnstructBC *bcRegion = unstructBCRegion->GetBCRegion(iBCRegionUnstruct);
        int bcType = bcRegion->GetBCType();

        if (bcType < 0)
        {
            continue;
        }

        if (IsWall(bcType))
        {
            SetWallBCGhostDQLUSGS(grid, bcRegion, dq);
        }
        //FARFIELD, SYMMETRY, INFLOW and OUTFLOW can use the same module.
        else if (bcType == PHENGLEI::FARFIELD || bcType == PHENGLEI::SYMMETRY 
            || bcType == PHENGLEI::INFLOW   || bcType == PHENGLEI::OUTFLOW)
        {
            SetFarfieldBCGhostDQLUSGS(grid, bcRegion, dq);
        }
        else
        {
            vector<int> *faceIndex = bcRegion->GetFaceIndex();
            for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
            {
                int iFace = *iter;
                int re = rightCellOfFace[iFace];
                for (int m = 0; m < nTurbulenceEquation; ++ m)
                {
                    dq[m][re] = 0.0; 
                }
            }
        }
    }
}

void TurbSolverUnstr::SetWallBCGhostDQLUSGS(Grid *gridIn, UnstructBC *bcRegionUnstruct, RDouble **dq)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int *rightCellOfFace  = grid->GetRightCellOfFace();

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();
    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        int iFace = *iter;
        int re = rightCellOfFace[iFace];
        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            dq[m][re] = 0.0;
        }
    }
}

void TurbSolverUnstr::SetFarfieldBCGhostDQLUSGS(Grid *gridIn, UnstructBC *bcRegionUnstruct, RDouble **dq)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int *leftCellOfFace  = grid->GetLeftCellOfFace();
    int *rightCellOfFace  = grid->GetRightCellOfFace();

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();
    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();
    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        int iFace = *iter;
        int le = leftCellOfFace[iFace];
        int re = rightCellOfFace[iFace];
        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            dq[m][re] = dq[m][le]; 
        }
    }
}


void TurbSolverUnstr::ModifyResiduals(Grid *grid)
{
    ModifyResidualOnWallBoundary(grid);
}

void TurbSolverUnstr::LoadResiduals(Grid *gridIn, FieldProxy *rhsProxy)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell    = grid->GetNTotalCell();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal        = numberOfTotalCell + numberOfBoundaryFace;

    RDouble **rhs = rhsProxy->GetField_UNS();

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    RDouble **residualTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));

    for (int m = 0; m < nTurbulenceEquation; ++ m)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            residualTurbulence[m][iCell] = - rhs[m][iCell];
        }
    }
}

RDouble TurbSolverUnstr::UnsteadyConvergence(Grid *gridIn)
{
    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int isUnsteady = parameters->GetIsUnsteady();
    if (!isUnsteady) return zero;

    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell    = grid->GetNTotalCell();

    RDouble **qPrimitive  = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
    RDouble **qTurbulenceUnsteadyN1 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb_unsteady_n1"));
    RDouble **qTurbulenceUnsteadyN2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb_unsteady_n2"));
    RDouble **residualTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));

    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    RDouble *primitive0 = new RDouble[nTurbulenceEquation];
    RDouble *primitive1 = new RDouble[nTurbulenceEquation];
    RDouble *primitive2 = new RDouble[nTurbulenceEquation];

    RDouble *qConservative0 = new RDouble[nTurbulenceEquation];
    RDouble *qConservative1 = new RDouble[nTurbulenceEquation];
    RDouble *qConservative2 = new RDouble[nTurbulenceEquation];

    RDouble sum1 = zero;
    RDouble sum2 = zero;

    using namespace IDX;

    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
    {
        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            primitive0[m] = qTurbulence[m][iCell];
            primitive1[m] = qTurbulenceUnsteadyN1[m][iCell];
            primitive2[m] = qTurbulenceUnsteadyN2[m][iCell];
        }

        RDouble coefficient = 1.0;

        if (nTurbulenceEquation >= 2)
        {
            coefficient = qPrimitive[IR][iCell];
        }

        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            qConservative0[m] = coefficient * primitive0[m];
            qConservative1[m] = coefficient * primitive1[m];
            qConservative2[m] = coefficient * primitive2[m];
        }

        //! nl is needed.
        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            //! Here, res must be transformed into dq，that is, res is not equal to rhs but dq.
            RDouble dqp = residualTurbulence[m][iCell];        // qn+1,p+1 - qn+1,p
            RDouble dqn = qConservative0[m] - qConservative1[m];  // qn+1,p+1 - qn
            sum1 += dqp * dqp;
            sum2 += dqn * dqn;
        }
    }

    delete []primitive0;
    delete []primitive1;
    delete []primitive2;
    delete []qConservative0;
    delete []qConservative1;
    delete []qConservative2;

    RDouble cvg = sqrt(ABS(sum1 / (sum2 + SMALL)));
    return cvg;
}

void TurbSolverUnstr::UpdateUnsteadyFlow(Grid *gridIn)
{
    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int isUnsteady = parameters->GetIsUnsteady();
    if (!isUnsteady) return;

    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell    = grid->GetNTotalCell();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal        = numberOfTotalCell + numberOfBoundaryFace;

    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
    RDouble **qTurbulenceUnsteadyN1 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb_unsteady_n1"));
    RDouble **qTurbulenceUnsteadyN2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb_unsteady_n2"));

    RDouble **residualTurbulenceUnsteadyN1 = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb_unsteady_n1"));
    RDouble **residualTurbulenceUnsteadyN2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb_unsteady_n2"));
    RDouble **residualTurbulenceUnsteadyTemporary = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb_unsteady_tmp"));

    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    for (int m = 0; m < nTurbulenceEquation; ++ m)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            qTurbulenceUnsteadyN2[m][iCell] = qTurbulenceUnsteadyN1[m][iCell];
        }
    }

    for (int m = 0; m < nTurbulenceEquation; ++ m)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            qTurbulenceUnsteadyN1[m][iCell] = qTurbulence[m][iCell];
        }
    }

    for (int m = 0; m < nTurbulenceEquation; ++ m)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            residualTurbulenceUnsteadyN2[m][iCell] = residualTurbulenceUnsteadyN1[m][iCell];
        }
    }

    for (int m = 0; m < nTurbulenceEquation; ++ m)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            residualTurbulenceUnsteadyN1[m][iCell] = residualTurbulenceUnsteadyTemporary[m][iCell]; //! Here the current outerIterationStep is over, the value of the stored resTmp should be assigned to resn1 for the next outerIterationStep. It should be noticed that resTmp only contain the inviscid and viscous flux.
        }
    }

    //! Statistical variables for DES simulation.
    int ifStaticsFlowField = parameters->GetIfStaticsFlowField();
    if (ifStaticsFlowField > 0)
    {
        int outerIterationStep = GlobalDataBase::GetIntParaFromDB("outnstep");
        int startStatisticStep = GlobalDataBase::GetIntParaFromDB("startStatisticStep");

        if (outerIterationStep >= startStatisticStep)
        {
            int numberOfStatisticalStep = GlobalDataBase::GetIntParaFromDB("nStatisticalStep");
            RDouble statisticalTimePeriod = GlobalDataBase::GetDoubleParaFromDB("statisticalTimePeriod");
            RDouble physicalTimeStep = GlobalDataBase::GetDoubleParaFromDB("physicalTimeStep");

            RDouble c1 = 1.0 / (numberOfStatisticalStep * 1.0);
            if(statisticalTimePeriod > 0.0)
            {
                c1 = MAX(c1, physicalTimeStep / statisticalTimePeriod);
            }
            RDouble c2 = 1.0 - c1;

            RDouble **qAverageTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("qAverageTurb"));
            for (int m = 0; m < nTurbulenceEquation; ++ m)
            {
                for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
                {
                    //qAverageTurbulence[m][iCell] = c1 * qAverageTurbulence[m][iCell] + c2 * qTurbulence[m][iCell];
                    qAverageTurbulence[m][iCell] = c2 * qAverageTurbulence[m][iCell] + c1 * qTurbulence[m][iCell];
                }
            }
        }
    }
}

void TurbSolverUnstr::DualTimeSource(Grid *gridIn)
{
    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int isUnsteady = parameters->GetIsUnsteady();
    if (!isUnsteady) return;

    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell    = grid->GetNTotalCell();

    RDouble **qPrimitive  = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
    RDouble **qTurbulenceUnsteadyN1 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb_unsteady_n1"));
    RDouble **qTurbulenceUnsteadyN2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb_unsteady_n2"));

    RDouble **residualTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));
    RDouble **residualTurbulenceUnsteadyN1 = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb_unsteady_n1"));
    RDouble **residualTurbulenceUnsteadyN2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb_unsteady_n2"));
    RDouble **residualTurbulenceUnsteadyTemporary = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb_unsteady_tmp"));

    RDouble *volume   = grid->GetCellVolume();
    RDouble *volumeN1 = grid->GetCellVolume();
    RDouble *volumeN2 = grid->GetCellVolume();

    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    RDouble *primitive0 = new RDouble[nTurbulenceEquation];
    RDouble *primitive1 = new RDouble[nTurbulenceEquation];
    RDouble *primitive2 = new RDouble[nTurbulenceEquation];

    RDouble *qConservative0 = new RDouble[nTurbulenceEquation];
    RDouble *qConservative1 = new RDouble[nTurbulenceEquation];
    RDouble *qConservative2 = new RDouble[nTurbulenceEquation];

    //! Store the historical sum of InviscidFlux, ViscousFlux, et al., exclusive of DualTimeSourceFlux, used in UpdateUnsteadyFlow for the Rn in the dualtime method.
    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
    {
        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            residualTurbulenceUnsteadyTemporary[m][iCell] = residualTurbulence[m][iCell];
        }
    }

    //! Computation of dualtime coefficients, including three coefficients for Residual of R(p), R(n) and R(n-1);
    //! and three coefficients for conservative variables of qcsv(p), qcsv(n), and qcsv(n-1).
    RDouble dualTimeCoefficient[7];
    const int methodOfDualTime = GlobalDataBase::GetIntParaFromDB("methodOfDualTime");
    ComputeDualTimeCoefficient(methodOfDualTime, dualTimeCoefficient);

    RDouble dualTimeResC1 = dualTimeCoefficient[0];
    RDouble dualTimeResC2 = dualTimeCoefficient[1];
    RDouble dualTimeResC3 = dualTimeCoefficient[2];
    RDouble dualTimeQC1   = dualTimeCoefficient[3];
    RDouble dualTimeQC2   = dualTimeCoefficient[4];
    RDouble dualTimeQC3   = dualTimeCoefficient[5];

    using namespace IDX;

    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
    {
        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            primitive0[m] = qTurbulence[m][iCell];
            primitive1[m] = qTurbulenceUnsteadyN1[m][iCell];
            primitive2[m] = qTurbulenceUnsteadyN2[m][iCell];
        }

        RDouble coefficient = 1.0;

        if (nTurbulenceEquation >= 2)
        {
            coefficient  = qPrimitive[IR][iCell];
        }

        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            qConservative0[m] = coefficient * primitive0[m];
            qConservative1[m] = coefficient * primitive1[m];
            qConservative2[m] = coefficient * primitive2[m];
        }

        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            RDouble dualSrcRes = dualTimeResC1 * residualTurbulence[m][iCell] +
                                 dualTimeResC2 * residualTurbulenceUnsteadyN1[m][iCell] +
                                 dualTimeResC3 * residualTurbulenceUnsteadyN2[m][iCell];

            RDouble dualSrcQ   = dualTimeQC1 * qConservative0[m]* volume  [iCell] + 
                                 dualTimeQC2 * qConservative1[m]* volumeN1[iCell] + 
                                 dualTimeQC3 * qConservative2[m]* volumeN2[iCell];

            residualTurbulence[m][iCell] = dualSrcRes + dualSrcQ;
        }
    }

    delete []primitive0;
    delete []primitive1;
    delete []primitive2;
    delete []qConservative0;
    delete []qConservative1;
    delete []qConservative2;
}

void TurbSolverUnstr::FillField(Grid *gridIn, FieldProxy *fieldProxy, RDouble value)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell    = grid->GetNTotalCell();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal        = numberOfTotalCell + numberOfBoundaryFace;

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    RDouble **field = fieldProxy->GetField_UNS();

    for (int m = 0; m < nTurbulenceEquation; ++ m)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            field[m][iCell] = value;
        }
    }
}

void TurbSolverUnstr::FillField(Grid *gridIn, FieldProxy *field1Proxy, FieldProxy *field2Proxy)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell    = grid->GetNTotalCell();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal        = numberOfTotalCell + numberOfBoundaryFace;

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    RDouble **field1 = field1Proxy->GetField_UNS();
    RDouble **field2 = field2Proxy->GetField_UNS();

    for (int m = 0; m < nTurbulenceEquation; ++ m)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            field1[m][iCell] = field2[m][iCell];
        }
    }
}

FieldProxy * TurbSolverUnstr::CreateFieldProxy(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell    = grid->GetNTotalCell();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal        = numberOfTotalCell + numberOfBoundaryFace;

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    RDouble **field = NewPointer2<RDouble>(nTurbulenceEquation, numberOfTotal);

    FieldProxy *fieldProxy = new FieldProxy();

    fieldProxy->SetField_UNS(field, true);

    return fieldProxy;
}

FieldProxy * TurbSolverUnstr::GetFieldProxy(Grid *gridIn, const string &fieldName)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    RDouble **field = reinterpret_cast<RDouble **> (grid->GetDataPtr(fieldName));

    FieldProxy *fieldProxy = new FieldProxy();

    fieldProxy->SetField_UNS(field);

    return fieldProxy;
}

FieldProxy * TurbSolverUnstr::GetResidualProxy(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    RDouble **residualTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));

    FieldProxy *residualTurbulenceProxy = new FieldProxy();

    residualTurbulenceProxy->SetField_UNS(residualTurbulence);

    return residualTurbulenceProxy;
}

void TurbSolverUnstr::RecoverResidual(Grid *gridIn, FieldProxy *rhsProxy)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell    = grid->GetNTotalCell();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal        = numberOfTotalCell + numberOfBoundaryFace;

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    RDouble **residualTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));

    if (grid->GetLevel() != 0) 
    {
        RDouble **rhs = rhsProxy->GetField_UNS();
        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
            {
                residualTurbulence[m][iCell] = rhs[m][iCell];
            }
        }
    }
}

void TurbSolverUnstr::StoreRhsByResidual(Grid *gridIn, FieldProxy *rhsProxy)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell    = grid->GetNTotalCell();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal        = numberOfTotalCell + numberOfBoundaryFace;

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    RDouble **residualTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));
    RDouble **rhs = rhsProxy->GetField_UNS();

    for (int m = 0; m < nTurbulenceEquation; ++ m)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            rhs[m][iCell] = residualTurbulence[m][iCell];
        }
    }
}

void TurbSolverUnstr::InitResidual(Grid *gridIn, FieldProxy *rhsProxy)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell    = grid->GetNTotalCell();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal        = numberOfTotalCell + numberOfBoundaryFace;

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    RDouble **residualTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));
    RDouble **rhs = rhsProxy->GetField_UNS();

    for (int m = 0; m < nTurbulenceEquation; ++ m)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            residualTurbulence[m][iCell] = - rhs[m][iCell];
        }
    }
}

void TurbSolverUnstr::RestrictAllQ(Grid *fineGridIn, Grid *coarseGridIn)
{
    UnstructGrid *fineGrid   = UnstructGridCast(fineGridIn);
    UnstructGrid *coarseGrid = UnstructGridCast(coarseGridIn);

    RDouble **qTurbulenceFineGrid   = reinterpret_cast<RDouble **> (fineGrid->GetDataPtr("q_turb"));
    RDouble **qTurbulenceCoarseGrid = reinterpret_cast<RDouble **> (coarseGrid->GetDataPtr("q_turb"));

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    for (int m = 0; m < nTurbulenceEquation; ++ m)
    {
        RestrictQ(fineGrid, qTurbulenceFineGrid[m], coarseGrid, qTurbulenceCoarseGrid[m]);
    }

    RDouble *viscousOfFineGrid   = reinterpret_cast<RDouble *> (fineGrid->GetDataPtr("vist"));
    RDouble *viscousOfCoarseGrid = reinterpret_cast<RDouble *> (coarseGrid->GetDataPtr("vist"));
    RestrictQ(fineGrid, viscousOfFineGrid, coarseGrid, viscousOfCoarseGrid);
    //coarseGrid->SetGhostCellExceptInterface(cvist);
    coarseGrid->SetGhostCell(viscousOfCoarseGrid);
}

// Bell 20130513 mod
void TurbSolverUnstr::LoadFlux(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfBoundaryFace = grid->GetNBoundFace();

    int *leftCellOfFace  = grid->GetLeftCellOfFace();
    int *rightCellOfFace  = grid->GetRightCellOfFace();

    RDouble **flux = faceProxy->GetFlux();

    //! Determine if there are boundary faces.
    int nMid  = localStart;
    if (localEnd <= numberOfBoundaryFace)
    {
       //! If all boundary faces
       nMid = localEnd;
    }
    else if (localStart < numberOfBoundaryFace)
    {
       //! Part of them are boundary faces
       nMid = numberOfBoundaryFace;
    }

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    RDouble **residualTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));
    using namespace IDX;

    //! For boundary faces, remember re is ghost cell
    for (int i = localStart; i < nMid; ++ i)
    {
        int le,j;
        le = leftCellOfFace[i];
        j = i - localStart;
        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            residualTurbulence[m][le] -= flux[m][j];          
        }
    }

    //! Interior faces
    for (int i = nMid; i <localEnd; ++ i)
    {
        int le,re,j;
        le = leftCellOfFace[i];
        re = rightCellOfFace[i];
        j  = i - localStart;

        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            residualTurbulence[m][le] -= flux[m][j];
            residualTurbulence[m][re] += flux[m][j];
        }
    }
}

void TurbSolverUnstr::ModifyResidualOnWallBoundary(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    int modifyTurbulenceResidual = GlobalDataBase::GetIntParaFromDB("mod_turb_res");

    //! modify the residuals for the cells next to the wall
    if (modifyTurbulenceResidual != 1) return;

    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    RDouble **residualTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int viscoueType = parameters->GetViscousType();

    using namespace PHENGLEI;
    if (viscoueType == 3)
    {
        CompositeUnstructBC *unstructBCRegion = grid->GetUnstructBCRegion();
        int nBCRegionUnstruct = unstructBCRegion->GetnBCRegion();

        for (int iBCRegionUnstruct = 0; iBCRegionUnstruct < nBCRegionUnstruct; iBCRegionUnstruct ++)
        {
            UnstructBC *bcRegion = unstructBCRegion->GetBCRegion(iBCRegionUnstruct);
            int bcType = bcRegion->GetBCType();

            if (IsWall(bcType))
            {
                vector<int> *faceIndex = bcRegion->GetFaceIndex();
                for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter ++)
                {
                    // iFace is the face number in the set of faceIndex.
                    int iFace = *iter;
                    int le = leftCellOfFace[iFace];
                    int re = rightCellOfFace[iFace];
                    residualTurbulence[0][le] = 0.0;
                    residualTurbulence[0][re] = 0.0;    // Bell 20130508 add
                }
            }
        }
    }
}

void TurbSolverUnstr::ResetWallScalar(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *leftCellOfFace  = grid->GetLeftCellOfFace();
    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int viscoueType = parameters->GetViscousType();

    using namespace PHENGLEI;
    if (viscoueType == 3)
    {
        CompositeUnstructBC *unstructBCRegion = grid->GetUnstructBCRegion();
        int nBCRegionUnstruct = unstructBCRegion->GetnBCRegion();

        for (int iBCRegionUnstruct = 0; iBCRegionUnstruct < nBCRegionUnstruct; iBCRegionUnstruct++)
        {
            UnstructBC *bcRegion = unstructBCRegion->GetBCRegion(iBCRegionUnstruct);
            int bcType = bcRegion->GetBCType();

            if (IsWall(bcType))
            {
                vector<int> *faceIndex = bcRegion->GetFaceIndex();
                for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
                {
                    // iFace is the face number in the set of faceIndex.
                    int iFace = *iter;
                    int le = leftCellOfFace[iFace];
                    qTurbulence[0][le] = 0.0;
                }
            }
        }
    }
}

void TurbSolverUnstr::SourceFluxOneEquation(Grid *grid)
{
    ComputeLengthScaleofOneEquationModel(grid);

    //if (turbulentSourceKind == 0)
    //{
          if(TurbEquationFormation() == FORMATION_SA_CFL3D)
          {
              SourceFluxOneEquationOriginalCFL3D(grid);
          }
          else
          {
              SourceFluxOneEquationOriginal(grid);
          }        
}

// Bell 20130513 mod
void TurbSolverUnstr::SourceFluxOneEquationOriginal(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell    = grid->GetNTotalCell();

    RDouble **gradientPrimtiveVariableX = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarX"));
    RDouble **gradientPrimtiveVariableY = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarY"));
    RDouble **gradientPrimtiveVariableZ = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarZ"));

    RDouble *dudx = gradientPrimtiveVariableX[1];
    RDouble *dudy = gradientPrimtiveVariableY[1];
    RDouble *dudz = gradientPrimtiveVariableZ[1];

    RDouble *dvdx = gradientPrimtiveVariableX[2];
    RDouble *dvdy = gradientPrimtiveVariableY[2];
    RDouble *dvdz = gradientPrimtiveVariableZ[2];

    RDouble *dwdx = gradientPrimtiveVariableX[3];
    RDouble *dwdy = gradientPrimtiveVariableY[3];
    RDouble *dwdz = gradientPrimtiveVariableZ[3];

    RDouble **dqdx = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradTurbulenceX"));
    RDouble **dqdy = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradTurbulenceY"));
    RDouble **dqdz = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradTurbulenceZ"));

    RDouble **q         = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
    RDouble **residualTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));
    RDouble **spectrumTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("spec_turb"));
    RDouble *viscousLaminar = reinterpret_cast<RDouble *> (grid->GetDataPtr("visl"));

    RDouble *lengthScale = this->GetLengthScale(gridIn);
    RDouble *nxs,*nys,*nzs,*ns,*vol;

    nxs = grid->GetFaceNormalX();
    nys = grid->GetFaceNormalY();
    nzs = grid->GetFaceNormalZ();
    ns  = grid->GetFaceArea();
    vol = grid->GetCellVolume();

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    RDouble oRefReNumber = parameters->GetoRefReNumber();

    using namespace IDX;

    //!work    = temporary storage array, contains on input the
    //!          velocity gradients:
    //!          work(n,1) = ui,i
    //!          work(n,2) = u1,1**2 + u2,2**2 + u3,3**2
    //!          work(n,3) = (u1,2+u2,1)**2 + (u1,3+u3,1)**2 + (u2,3+u3,2)**2
    //!          work(n,4) = sqrt((u3,2-u2,3)**2 + (u1,3-u3,1)**2
    //!                                           + (u2,1-u1,2)**2)
    //!  calculation of dnu/dxi
    //!  dnu/dx in work(.,5), dnu/dy in work(.,6), dnu/dz in work(.,7)

    //von karman constant
    const RDouble karm   = 0.41;
    const RDouble cb1    = 0.1355;
    const RDouble cb2    = 0.622;
    const RDouble cv1    = 7.1;
    const RDouble cw2    = 0.3;
    const RDouble cw3    = 2.0;
    const RDouble sigma  = 2.0/3.0;
    const RDouble sdilim = 1.0e9;

    RDouble kap2,oKap2;
    kap2  = karm * karm;
    oKap2 = one / kap2;

    RDouble osigma = 1.0 / sigma;
    RDouble cw1    = cb1 * oKap2 + (one + cb2) / sigma;
    RDouble cv13   = cv1 * cv1 * cv1;
    RDouble cw36   = pow(cw3,6);
    RDouble or6    = one / six;

    RDouble production, diffsource, wallDest, source;
    RDouble cb2s, cw1Kap2, rs, std, ostd, nuetRs;
    RDouble prod, dest, prodp, destp, prde, prdep, part1, part2;
    RDouble distance2, oDistance2, xsi, xsi2, xsi3, fv1, fv2, ft2;
    RDouble r, r5, g, g6, fw, grd2;
    RDouble dsdnu, dfv2dk, dfwdg, dgdr, drdnu, dfwdnu;
    RDouble nueSA, volume, ftrans, olam, cmachfix;
    RDouble sBar, omega;
    RDouble saC2, saC3;

    saC2 = 0.7;
    saC3 = 0.9;

    cb2s = osigma * cb2;
    cw1Kap2 = cw1 * kap2;

    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
    {
        nueSA = qTurbulence[ISA][iCell];

        //! Omega: magnitude of the vorticity.
        omega = sqrt(SQR(dwdy[iCell] - dvdz[iCell]) +
                    SQR(dudz[iCell] - dwdx[iCell]) +
                    SQR(dvdx[iCell] - dudy[iCell]));   
       
        volume   = vol[iCell];
        cmachfix = nueSA * (dudx[iCell] + dvdy[iCell] + dwdz[iCell]);
        cmachfix = 0.0;
        olam     = q[IR][iCell] / (viscousLaminar[iCell] + SMALL);
        distance2= lengthScale[iCell] * lengthScale[iCell];
        oDistance2 = one / distance2;

        xsi      = nueSA * olam + SMALL;
        xsi3     = xsi * xsi * xsi;
        fv1      = xsi3 / (xsi3 + cv13);
        fv2      = one - xsi / (one + xsi * fv1);

        //! Here, absorb oRefReNumber into rs.
        rs     = nueSA * oKap2 * oDistance2 * oRefReNumber;
        nuetRs = nueSA * rs;

        //! sBar = fv2 * nueSA / (k2 * distance2).
        sBar    = rs * fv2;

        if (sBar >= - saC2 * omega)
        {
            std = omega + sBar;
        }
        else
        {
            std = omega + omega * (saC2 * saC2 * omega + saC3 * sBar) / ((saC3 - 2 *saC2) * omega - sBar);
        }

        //std     = MAX(str + rs * fv2, 0.3*str);
        ostd    = one / (std + SMALL);
        r       = rs * ostd;

        r    = MIN(r, static_cast<RDouble>(ten));
        r5   = pow(r, 5);
        g    = r + cw2 * (r5 * r - r);
        g6   = pow(g, 6);
        fw   = g * pow((static_cast<RDouble>(one) + cw36) / (g6 + cw36), or6);

        //! Calculate "negative" production term in order to assure a zero.
        //! Solution in the laminar region.

        ftrans = zero;
        xsi2   = xsi * xsi;
        ft2    = 0.0;       //! Do Not consider transition here!
        //ft2    = ct3 * exp(- ct4 * xsi2) * ftrans;
        grd2 = SQR(dqdx[ISA][iCell]) + SQR(dqdy[ISA][iCell]) + SQR(dqdz[ISA][iCell]);

        //! cb1  * omega * nueSA + cb1  * nueSA_rs * fv2.
        production  = cb1 * (one - ft2) * std * nueSA;       

        //! cb2s * rho * gradnue2.
        diffsource  = cb2s * grd2 * oRefReNumber;           

        //! ft2 = 0.0, cw1Kap2 = Cw1 * Kap2, so,
        //! cw1Kap2 * fw * nuetRs = (Cw1 * Kap2 * fw * nueSA) * (nueSA2 * oKap2 * oDistance2)
        //!                       = Cw1 * fw * (nueSA/d)^2
        wallDest   = (cw1Kap2 * fw - cb1 * ft2) * nuetRs; //! The Re has been included in nuetRs.

        diffsource = MIN(sdilim * production, diffsource);
        source     = production - wallDest;

        if(TurbEquationFormation() == FORMATION_SA_NSMB)
        {
            //source  += diffsource;                        //! Bell 20130513 delete for modsa.
        }
        else if(TurbEquationFormation() == FORMATION_ORGSA_SST)
        {
            source  += diffsource;
        }

        residualTurbulence[ISA][iCell] += (source) * volume;

        //! Spectrum.
        dfv2dk = (three * cv13 * pow(fv1/(xsi+SMALL), 2) - one) / pow(one + xsi * fv1, 2);
        dsdnu  = oRefReNumber * oKap2 * oDistance2 * (fv2 + xsi * dfv2dk);
        dfwdg   = fw / (g + SMALL) * (one - g6 / (g6 + cw36));
        dgdr    = one + cw2 * (six * r5 - one);
        drdnu   = oRefReNumber * oKap2 * oDistance2 * ostd * (one - nueSA * ostd * dsdnu);
        dfwdnu  = dfwdg * dgdr * drdnu;

        //dft2dnu = - two * ct3 * ct4 * xsi * exp (- ct4 * xsi2) * olam * ftrans;

        prod = cb1 * (one - ft2) * std;
        dest = (cw1*fw - cb1*ft2*oKap2) * nueSA * oDistance2 * oRefReNumber;

        prodp = cb1 * ((one - ft2) * dsdnu);
        destp = (cw1 * (fw + dfwdnu * nueSA) - cb1 * ft2 * oKap2) * oDistance2 * oRefReNumber;

        prde  = prod  - dest;
        prdep = prodp - destp;

        part1 = - half * (prde  - ABS(prde));
        part2 = - half * (prdep - ABS(prdep)) * nueSA;

        spectrumTurbulence[ISA][iCell] += (part1 + part2) * volume;        
    }
}

void TurbSolverUnstr::SourceFluxOneEquationOriginalCFL3D(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell    = grid->GetNTotalCell();

    RDouble **gradientPrimtiveVariableX = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarX"));
    RDouble **gradientPrimtiveVariableY = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarY"));
    RDouble **gradientPrimtiveVariableZ = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarZ"));

    RDouble *dudx = gradientPrimtiveVariableX[1];
    RDouble *dudy = gradientPrimtiveVariableY[1];
    RDouble *dudz = gradientPrimtiveVariableZ[1];

    RDouble *dvdx = gradientPrimtiveVariableX[2];
    RDouble *dvdy = gradientPrimtiveVariableY[2];
    RDouble *dvdz = gradientPrimtiveVariableZ[2];

    RDouble *dwdx = gradientPrimtiveVariableX[3];
    RDouble *dwdy = gradientPrimtiveVariableY[3];
    RDouble *dwdz = gradientPrimtiveVariableZ[3];

    RDouble **q         = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
    RDouble **residualTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));
    RDouble **spectrumTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("spec_turb"));
    RDouble *viscousLaminar = reinterpret_cast<RDouble *> (grid->GetDataPtr("visl"));

    RDouble *lengthScale = this->GetLengthScale(gridIn);
    RDouble *nxs,*nys,*nzs,*ns,*vol;

    nxs = grid->GetFaceNormalX();
    nys = grid->GetFaceNormalY();
    nzs = grid->GetFaceNormalZ();
    ns  = grid->GetFaceArea();
    vol = grid->GetCellVolume();

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    RDouble oRefReNumber = parameters->GetoRefReNumber();

    using namespace IDX;

    //von karman constant
    const RDouble karm   = 0.41;
    const RDouble cb1    = 0.1355;
    const RDouble cb2    = 0.622;
    const RDouble cv1    = 7.1;
    const RDouble cw2    = 0.3;
    const RDouble cw3    = 2.0;
    const RDouble sigma  = 2.0/3.0;

    RDouble kap2,oKap2;
    kap2  = karm * karm;
    oKap2 = one / kap2;

    RDouble osigma = 1.0 / sigma;
    RDouble cw1    = cb1 * oKap2 + (one + cb2) / sigma;
    RDouble cv13   = cv1 * cv1 * cv1;
    RDouble cw36   = pow(cw3, 6.0);
    RDouble or6    = one / six;

    RDouble ct3 = 1.2;
    RDouble ct4 = 0.5; 

    RDouble production, wallDest, source;
    RDouble cb2s, cw1Kap2, S, oS;
    RDouble distance2, oDistance2, xsi, xsi3, fv1, fv2, ft2;
    RDouble r, r5, g, g6, fw;
    RDouble nue, nueSA, volume, ftrans, rho;
    RDouble sBar, omega;
    RDouble saC2, saC3;

    saC2 = 0.7;
    saC3 = 0.9;

    cb2s = osigma * cb2;
    cw1Kap2 = cw1 * kap2;

    RDouble machinZero = pow(10.0, -machZeroExp + 1);

    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
    {
        nueSA = qTurbulence[ISA][iCell];
        rho   = q[IR][iCell];
        nue   = viscousLaminar[iCell] / (rho + SMALL);

        //! Omega: magnitude of the vorticity.
        omega = sqrt(SQR(dwdy[iCell] - dvdz[iCell]) +
                     SQR(dudz[iCell] - dwdx[iCell]) +
                     SQR(dvdx[iCell] - dudy[iCell]));   

        volume   = vol[iCell];
  
        distance2= lengthScale[iCell] * lengthScale[iCell];
        oDistance2 = one / distance2;

        xsi      = nueSA / nue;
        xsi3     = xsi * xsi * xsi;
        fv1      = xsi3 / (xsi3 + cv13);
        fv2      = one - xsi / (one + xsi * fv1);

        //! Important, this term must be non-dimensional with Reynolds.
        //! Here, absorb oRefReNumber into rs(term0): nueSA / (k2 * distance2).
        RDouble term0 = nueSA * oKap2 * oDistance2 * oRefReNumber;        

        //!sBar = fv2 * nueSA / (k2 * distance2).
        sBar    = term0 * fv2;

        if (sBar >= - saC2 * omega)
        {
            S = omega + sBar;
        }
        else
        {
            S = omega + omega * (saC2 * saC2 * omega + saC3 * sBar) / ((saC3 - 2 *saC2) * omega - sBar);
        }
        oS    = one / (S + SMALL);

        r     = term0 * oS;
        r    = MIN(r, static_cast<RDouble>(ten));

        r5   = pow(r, 5);
        g    = r + cw2 * (r5 * r - r);
        g    = MAX(g, machinZero);

        g6   = pow(g, 6);
        fw   = g * pow((static_cast<RDouble>(one) + cw36) / (g6 + cw36), or6);

        //! ftrans: transition or not -- 1/0.
        //! More work should be done to compare the influence of ft2.
        ftrans = zero;  // one.
        ct3   *= ftrans;
        ft2    = ct3 * exp(- ct4 * xsi * xsi);

        //! Using CFL3D terms.
        //! CFL3D user's manual v5.0, and, open source code.
        //! Here, consider omega only, since the sBar has been moved into 
        //! the wallDest term because they all have the Reynolds.
        RDouble term1_CFL3D = cb1 * (one - ft2) * omega;
        RDouble term2_CFL3D = cb1 * ((one - ft2) * fv2 + ft2) * oKap2 - cw1 * fw;

        production  = term1_CFL3D * nueSA;     
        wallDest    = term2_CFL3D * nueSA * nueSA * oDistance2 * oRefReNumber;
        source      = production + wallDest;

        residualTurbulence[ISA][iCell] += source * volume;
        
        //! Spectrum, reference from CFL3D codes.
        RDouble derivativePart1 = 2.0 * term2_CFL3D * nueSA * oDistance2;
        RDouble dfv1  = (fv1 - fv1 * fv1) * 3.0 / nueSA;
        RDouble dfv2  = (fv2 - 1.0) / nueSA + (1.0 - fv2) * (1.0 - fv2) * (fv1 / nueSA + dfv1);
        RDouble dft2  = -(2.0 * ct4 * nueSA / (nue * nue)) * ft2;
        RDouble drr   = r / nueSA - r * r * (fv2 / nueSA + dfv2); 
        RDouble dgg   = (1.0 - cw2 + 6.0 * cw2 * pow(r, 5.0)) * drr;

        g  = MAX(g, 10.0 * machinZero); // fix for single prec. via CFL3D via TLNS3D via FUN3D vi Tim Barth
        g6 = pow(g, 6);

        RDouble dfw   = pow((1.0 + cw36) / (g6 + cw36), or6) - 
                       (pow(1.0 + cw36, or6) / (pow(g6 + cw36,(7.0 / 6.0)))) * g6;      
        dfw *= dgg;  
        RDouble derivativePart2  = oDistance2 * nueSA * nueSA * 
            (cb1 * oKap2 * (dfv2 - ft2 * dfv2 - fv2 * dft2 + dft2)-cw1 * dfw);

        //! The source spectrum terms lie on left side of equation, so minus it.
        spectrumTurbulence[ISA][iCell] -= (derivativePart1 + derivativePart2) * oRefReNumber * volume;

        //! The following term is negatived. it may delay the converge, although may improve robust.
        //! Test it more in future.
        //spectrumTurbulence[ISA][iCell] -= term1_CFL3D * volume;
    }
}

void TurbSolverUnstr::SourceFluxTwoEquation(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    RDouble **primitiveVariables = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble *viscousLaminar = reinterpret_cast<RDouble * > (grid->GetDataPtr("visl"));
    RDouble *viscousTurbulence = reinterpret_cast<RDouble * > (grid->GetDataPtr("vist"));
    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
    RDouble **spectrumTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("spec_turb"));
    RDouble **residualTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));

    RDouble *wallDistance = grid->GetWallDist();
    RDouble *xc   = grid->GetCellCenterX();
    RDouble *xfn  = grid->GetFaceNormalX();
    RDouble *yfn  = grid->GetFaceNormalY();
    RDouble *zfn  = grid->GetFaceNormalZ();
    RDouble *area = grid->GetFaceArea();
    RDouble *vol  = grid->GetCellVolume();

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    RDouble refReNumber  = parameters->GetRefReNumber();
    RDouble oRefReNumber = parameters->GetoRefReNumber();

    string viscousName = parameters->GetViscousName();

    int neasm = -1;
    GlobalDataBase::GetData("neasm",&neasm, PHINT, 1);

    int SSTProductType = parameters->GetSSTProductType();

    int numberOfTotalCell    = grid->GetNTotalCell();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal        = numberOfTotalCell + numberOfBoundaryFace;

    RDouble cblend, density, ke, kw;
    RDouble sij2, divv, prodk, dissk, prodw, dissw, cdkww, crss, mut;
    RDouble srck, srcw, fskn, fswn, diak, diaw, oork, oorw, volume;
    /*RDouble turbXk, turbXw, turbXk2, turbFk2, turbFw;*/
    RDouble s11, s22, s33, s12, s13, s23, w12, w13, w23;
    RDouble vort2;

    RDouble turbFbetas, turbFbeta, turbBeta, turbAlphaw/*, turbSigd*/;

    RDouble *dudx = new RDouble[numberOfTotal]();
    RDouble *dudy = new RDouble[numberOfTotal]();
    RDouble *dudz = new RDouble[numberOfTotal]();
    
    RDouble *dvdx = new RDouble[numberOfTotal]();
    RDouble *dvdy = new RDouble[numberOfTotal]();
    RDouble *dvdz = new RDouble[numberOfTotal]();
    
    RDouble *dwdx = new RDouble[numberOfTotal]();
    RDouble *dwdy = new RDouble[numberOfTotal]();
    RDouble *dwdz = new RDouble[numberOfTotal]();

    using namespace IDX;
    grid->CompGradient(primitiveVariables[IU], dudx, dudy, dudz);
    grid->CompGradient(primitiveVariables[IV], dvdx, dvdy, dvdz);
    grid->CompGradient(primitiveVariables[IW], dwdx, dwdy, dwdz);

    RDouble **aniss = 0;
    RDouble *sengy  = 0;

    if (neasm > 0)
    {
        aniss = reinterpret_cast<RDouble **> (grid->GetDataPtr("aniss"));
        sengy = reinterpret_cast<RDouble * > (grid->GetDataPtr("sengy"));
    }

    int freeturbIntensitySRModify = 0;
    if (GlobalDataBase::IsExist("freeturbIntensitySRModify", PHINT, 1))
    {
        freeturbIntensitySRModify = GlobalDataBase::GetIntParaFromDB("freeturbIntensitySRModify");
    }
    double freeDecayXLocation = 0.0;
    if (GlobalDataBase::IsExist("freeDecayXLocation", PHDOUBLE, 1))
    {
        freeDecayXLocation = GlobalDataBase::GetDoubleParaFromDB("freeDecayXLocation");
    }

    using namespace IDX;
    RDouble refMaNumber = parameters->GetRefMachNumber();
    RDouble transition_ca1, transition_ca2, transition_ce1, transition_ce2, transition_cct, transition_s1, transition_ccf;
    RDouble SpSdlimit = MIN(1.0 + 0.00075 * refMaNumber * refMaNumber * refMaNumber * refMaNumber, 1.2);

    int isCC = 0;
    RDouble compressibleFix = 1.0 + 0.40 * pow(refMaNumber, 1.75);

    int transitionType = parameters->GetTransitionType();

    RDouble *primitiveVariablesInfinity = reinterpret_cast<RDouble*> (GlobalDataBase::GetDataPtr("prim_inf"));
    RDouble *freeStreamTurbVar = parameters->GetFreeStreamTurbVar();
    RDouble freeStreamViscosity = GlobalDataBase::GetDoubleParaFromDB("freeStreamViscosity");

    if (transitionType == IREGAMA)
    {
        transition_ca1 = parameters->GetTransition_ca1();
        transition_ca2 = parameters->GetTransition_ca2();
        transition_ce1 = parameters->GetTransition_ce1();
        transition_ce2 = parameters->GetTransition_ce2();
        transition_cct = parameters->GetTransition_cct();
        transition_s1  = parameters->GetTransition_s1();
        transition_ccf = parameters->GetTransition_ccf();
    }

    //! K-OMEGA MODELS
    if (viscousName.substr(0, 6) == "2eq-kw")
    {
        //! MENTER'S K-OMEGA SST MODEL
        if (viscousName.substr(0, 17) == "2eq-kw-menter-sst")
        {
            //turbulence model constants
            RDouble SSTProductLimit = parameters->GetSSTProductLimit();
            GlobalDataBase::GetData("turb_fbetas", &turbFbetas, PHDOUBLE, 1);
            GlobalDataBase::GetData("turb_fbeta" , &turbFbeta , PHDOUBLE, 1);
            RDouble SST_beta1 = parameters->GetSST_beta1();
            RDouble SST_beta2 = parameters->GetSST_beta2();
            RDouble SST_betaStar = parameters->GetSST_betaStar();
            RDouble KW_sigmaW2 = parameters->GetKW_sigmaW2();
            RDouble SST_alphaw1 = parameters->GetSST_alphaw1();
            RDouble SST_alphaw2 = parameters->GetSST_alphaw2();

            RDouble *cross = reinterpret_cast<RDouble  *> (grid->GetDataPtr("cross"));
            RDouble *blend = reinterpret_cast<RDouble  *> (grid->GetDataPtr("blend"));
            RDouble *SpSdRatio = NULL;
            if (transitionType == IREGAMA)
            {
                 SpSdRatio = reinterpret_cast<RDouble *> (grid->GetDataPtr("SpSdRatio"));
            }

            for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
            {
                cblend  = blend[iCell];
                crss    = cross[iCell];
                density = primitiveVariables[IR][iCell];
                ke      = qTurbulence[IKE][iCell];
                kw      = qTurbulence[IKW][iCell];
                mut     = viscousTurbulence[iCell];
                volume  = vol[iCell];

                turbBeta   = cblend * SST_beta1   + (1.0 - cblend) * SST_beta2;
                turbAlphaw = cblend * SST_alphaw1 + (1.0 - cblend) * SST_alphaw2;

                s11 = dudx[iCell];
                s22 = dvdy[iCell];
                s33 = dwdz[iCell];
                s12 = half * (dudy[iCell] + dvdx[iCell]);
                s13 = half * (dudz[iCell] + dwdx[iCell]);
                s23 = half * (dvdz[iCell] + dwdy[iCell]);
                w12 = half * (dudy[iCell] - dvdx[iCell]);
                w13 = half * (dudz[iCell] - dwdx[iCell]);
                w23 = half * (dvdz[iCell] - dwdy[iCell]);
                
                sij2 = two * (s11*s11 + s22*s22 + s33*s33 + two*(s12*s12 + s13*s13 + s23*s23));
                divv = s11 + s22 + s33;

                //  Re-gama transition model
                RDouble gmeff = 0.0;
                if (transitionType == IREGAMA)
                {
                    RDouble mul = viscousLaminar[iCell];
                    RDouble um = primitiveVariables[IU][iCell];
                    RDouble vm = primitiveVariables[IV][iCell];
                    RDouble wm = primitiveVariables[IW][iCell];
                    RDouble intermittency = qTurbulence[IGAMA][iCell];
                    RDouble Rectabar = qTurbulence[IRECT][iCell];

                    RDouble wallDist = wallDistance[iCell];

                    RDouble vorx = dvdz[iCell] - dwdy[iCell];
                    RDouble vory = dudz[iCell] - dwdx[iCell];
                    RDouble vorz = dudy[iCell] - dvdx[iCell];

                    RDouble vorticity = DISTANCE(vorx, vory, vorz);
                    RDouble strainRate = sqrt(sij2);

                    RDouble absU = MAX(DISTANCE(um, vm, wm), SMALL);

                    RDouble RT     = ViscosityRatio(density, mul, ke, kw, refReNumber);
                    RDouble Rev    = ReynoldsNumberBasedOnStrainRate(density, wallDist, mul, strainRate, refReNumber);
                    RDouble Rectac = TransitionOnsetMomentumThicknessReynoldsUnstrCaliBrated(Rectabar);

                    RDouble Rw     = ReynoldsNumberBasedOnDissipation(density, wallDist, mul, kw, refReNumber);
                    RDouble Flength = HighReynoldsCorrectionOfFlength(Rw, FlengthUnstrCaliBrated(Rectabar));

                    RDouble Fonset  = TransitionOnsetFunction(Rev, Rectac, RT);
                    RDouble Fturb   = ControlFunctionFturb(RT);

                    RDouble production1OfGama = transition_ca1 * Flength * density * strainRate * sqrt(intermittency * Fonset);
                    RDouble production2OfGama = production1OfGama * (-transition_ce1 * intermittency);

                    RDouble productionOfGama = production1OfGama + production2OfGama;

                    RDouble specGamap = MIN(0.0, 1.5 * production2OfGama + 0.5 * production1OfGama) / (density * intermittency);

                    RDouble destruction2OfGama = transition_ca2 * Fturb * density * vorticity * intermittency;
                    RDouble destruction1OfGama = transition_ce2 * intermittency * destruction2OfGama;

                    RDouble destructionOfGama = destruction1OfGama + destruction2OfGama;
                    RDouble specGamad = MAX(0.0, 2.0 * destruction1OfGama + destruction2OfGama) / (density * intermittency);
                    
                    residualTurbulence[IGAMA][iCell] += (productionOfGama - destructionOfGama) * volume;
                    spectrumTurbulence[IGAMA][iCell] += (specGamad - specGamap) * volume;

                    RDouble Fctat = 1.0, Fctatcf = 1.0;
                    BlendingFunctionOfFctat(intermittency, Rw, Rectabar, vorticity, mul, density, absU, wallDist, refReNumber, transition_ce2, Fctat, Fctatcf);
                    RDouble tscl = TimeScaleInSourceTerm(density, absU, mul, mut, refReNumber);
                    RDouble TU = TurbulenceIntensity(absU, ke);
                    gmeff = SeparationCorrectionOfIntermittency(intermittency, Rev, Rectac, RT, Fctat, transition_s1, refMaNumber);

                    RDouble dUds = AccelerationAlongStreamline(um, vm, wm, dudx[iCell], dudy[iCell], dudz[iCell], dvdx[iCell], dvdy[iCell], dvdz[iCell], dwdx[iCell], dwdy[iCell], dwdz[iCell]);

                    RDouble momentumThickness = 0.0;
                    RDouble Rectat;
                    const int subIter = 10;
                    RDouble lamdacta0, lamdacta = 1.0;
                    for (int iter = 0; iter < subIter; ++iter)
                    {
                        lamdacta0 = lamdacta;
                        lamdacta = PressureGradientFunction(density, momentumThickness, mul, dUds, refReNumber);
                        lamdacta = lamdacta * compressibleFix;
                        lamdacta = MAX(-0.1, MIN(0.1, lamdacta));
                        RDouble Flamdacta = EmpiricalCorrelationOfFlamdacta(TU, lamdacta);
                        Rectat = EmpiricalCorrelationOfRectat(TU, Flamdacta);
                        momentumThickness = MomentumThickness(density, absU, mul, Rectat, refReNumber);
                        if (abs(lamdacta0 - lamdacta) <= SMALL)
                        {
                            break;
                        }
                    }
                    RDouble omegaStreamWise = abs((dwdy[iCell] - dvdz[iCell]) * um + (dudz[iCell] - dwdx[iCell]) * vm + (dvdx[iCell] - dudy[iCell]) * wm)/absU;
                    RDouble heightCrossFlow = wallDist * omegaStreamWise / absU;
                    RDouble Rescf = ReynoldsNumberBasedOnScf(heightCrossFlow, absU, density, mul, mut, refReNumber);

                    RDouble coecommon = transition_cct * (1.0 - Fctat) / MAX(tscl, SMALL);

                    RDouble prodRe = coecommon * density * Rectat;
                    RDouble destRe = coecommon * density * Rectabar;

                    RDouble prodRecf = transition_cct * density * transition_ccf * MIN(Rescf - Rectabar, 0.0) * Fctatcf / tscl;
                    RDouble specRe = destRe / (density * Rectabar);

                    residualTurbulence[IRECT][iCell] += (prodRe + prodRecf - destRe) * volume;
                    spectrumTurbulence[IRECT][iCell] += specRe * volume;

                }

                if (neasm < 0)
                {
                    if (SSTProductType == 0)
                    {
                        //! Boussinesq approximation,full production.
                        prodk = mut * (sij2 - two3rd * divv * divv) * oRefReNumber - two3rd * density * ke * divv;
                    }
                    else
                    {
                        vort2 = four * (w12 * w12 + w13 * w13 + w23 * w23);
                        prodk = mut * (vort2) * oRefReNumber - two3rd * density * ke * divv;
                    }
                }
                dissk = turbFbetas * SST_betaStar * density * ke * kw * refReNumber;

                //!limit production
                prodk = MIN(prodk, SSTProductLimit * dissk);
                prodw = turbAlphaw * density / (mut + SMALL) * prodk;
                dissw = turbFbeta * turbBeta * density * kw * kw * refReNumber;
                cdkww = two * (one - cblend) * density * KW_sigmaW2 * crss / (kw + SMALL) * oRefReNumber;
                //  Re-gama transition model
                if (transitionType == IREGAMA && freeturbIntensitySRModify == 1 && xc[iCell] < freeDecayXLocation)
                {
                    //! S-R Modify
                    RDouble keamb = freeStreamTurbVar[IKE];
                    RDouble kwamb = freeStreamTurbVar[IKW];
                    prodk += turbFbetas * SST_betaStar * density * keamb * kwamb * refReNumber;
                    prodw += turbFbeta * turbBeta * density * kwamb * kwamb * refReNumber;
                }

                ////  Re-gama transition model
                if (transitionType == IREGAMA)
                {
                    prodk = CorrectionOfProductionInKEquation(gmeff, prodk);
                    dissk = CorrectionOfDestructionInKEquation(gmeff, dissk);
                    SpSdRatio[iCell] = MIN(MAX(prodk / (dissk + SMALL), 1.0), SpSdlimit);
                }

                srck  = prodk - dissk;
                srcw  = prodw - dissw + cdkww;

                residualTurbulence[IKE][iCell] += (srck) * volume;
                residualTurbulence[IKW][iCell] += (srcw) * volume;

                if (neasm > 0)
                {
                    sengy[iCell] = - (srck) * volume;
                }

                diak = - two * turbFbetas * SST_betaStar * kw;
                diaw = - two * turbFbeta  * turbBeta  * kw - ABS(cdkww) / (density * kw + SMALL);

                fskn = half * (srck - ABS(srck));
                fswn = half * (srcw - ABS(srcw));

                oork = 1.0 / (density * ke + SMALL);
                oorw = 1.0 / (density * kw + SMALL);

                spectrumTurbulence[IKE][iCell] += - (diak + fskn * oork) * volume;
                spectrumTurbulence[IKW][iCell] += - (diaw + fswn * oorw) * volume;
            }
        }
    }

    delete []dudx;
    delete []dudy;
    delete []dudz;

    delete []dvdx;
    delete []dvdy;
    delete []dvdz;

    delete []dwdx;
    delete []dwdy;
    delete []dwdz;
}

void TurbSolverUnstr::Crossing(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell    = grid->GetNTotalCell();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal        = numberOfTotalCell + numberOfBoundaryFace;

    RDouble *dkedx = new RDouble [numberOfTotal];
    RDouble *dkedy = new RDouble [numberOfTotal];
    RDouble *dkedz = new RDouble [numberOfTotal];
    RDouble *dkwdx = new RDouble [numberOfTotal];
    RDouble *dkwdy = new RDouble [numberOfTotal];
    RDouble *dkwdz = new RDouble [numberOfTotal];

    RDouble *cross = reinterpret_cast<RDouble  *> (grid->GetDataPtr("cross"));

    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
    grid->CompGradient(qTurbulence[0], dkedx, dkedy, dkedz);
    grid->CompGradient(qTurbulence[1], dkwdx, dkwdy, dkwdz);

    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
    {
        cross[iCell] = dkedx[iCell] * dkwdx[iCell] + dkedy[iCell] * dkwdy[iCell] + dkedz[iCell] * dkwdz[iCell];
    }
    delete []dkedx;
    delete []dkedy;
    delete []dkedz;

    delete []dkwdx;
    delete []dkwdy;
    delete []dkwdz;
}

void TurbSolverUnstr::Blending(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell    = grid->GetNTotalCell();

    RDouble **q      = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble *viscousLaminar   = reinterpret_cast<RDouble * > (grid->GetDataPtr("visl"));
    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
    RDouble *cross = reinterpret_cast<RDouble  *> (grid->GetDataPtr("cross"));
    RDouble *blend = reinterpret_cast<RDouble  *> (grid->GetDataPtr("blend"));

    RDouble *wallDistance = grid->GetWallDist();

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    RDouble refReNumber = parameters->GetRefReNumber();
    RDouble oReynolds = 1.0 / refReNumber;
    RDouble oReynoldsSqr = oReynolds * oReynolds;
    int transitionType = parameters->GetTransitionType();

    RDouble betas = 0.09;
    RDouble sigw2 = 0.856;

    using namespace IDX;

    //! Compute maximum cross diffusion term across flowfield
    RDouble cdkwmax, crossDiff;
    cdkwmax = 0.0;

    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
    {
        RDouble rho  = q[IR][iCell];
        RDouble kw   = qTurbulence[IKW][iCell];
        RDouble crss = cross[iCell];

        crossDiff = 2.0 * (rho) * crss * sigw2 / (kw + SMALL);
    }

    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
    {
        RDouble rho   = q[IR][iCell];
        RDouble ke    = qTurbulence[IKE][iCell];
        RDouble kw    = qTurbulence[IKW][iCell];
        RDouble crss  = cross[iCell];
        RDouble visL  = viscousLaminar[iCell];
        RDouble dist  = wallDistance[iCell];
        RDouble dist2 = dist * dist;
        RDouble cdkw, term1, term2, term3, term4, arg1;

        //! calculate cd_kw
        crossDiff = 2.0 * (rho) * sigw2 * crss / (kw + SMALL);

        //! Original Menter CD_kw calculation
        cdkw = MAX(crossDiff, 1.0e-10);

        //! calculate arg1
        term1 = sqrt(ABS(ke)) * oReynolds / (betas * dist * kw + SMALL);
        term2 = 500.0 * visL * oReynoldsSqr / (rho * kw * dist2 + SMALL);
        term3 = MAX(term1, term2);
        term4 = 4.0 * rho * sigw2 * ke / (cdkw * dist2 + SMALL);
        arg1  = MIN(term3, term4);

        //! Compute the blending function fbsl
        blend[iCell] = tanh(arg1 * arg1 * arg1 * arg1);

        if (transitionType == IREGAMA)
        {
            CorrectionOfBlendingFunctionInSST(rho, dist, visL, ke, refReNumber, blend[iCell]);
        }
    }

    grid->SetGhostCellExceptInterface(blend);
}

FaceProxy * TurbSolverUnstr::CreateFaceProxy(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast (gridIn);

    FaceProxy *faceProxy = CreateFaceProxyTurbulence(grid);
    faceProxy->SetNext(CreateFaceProxyNS(grid));

    return faceProxy;
}

FaceProxy * TurbSolverUnstr::CreateFaceProxyNS(Grid *gridIn)
{
    int nl;
    GlobalDataBase::GetData("nl", &nl, PHINT, 1);

    int nchem;
    GlobalDataBase::GetData("nchem", &nchem, PHINT, 1);

    int numberOfEquation = nl + nchem;

    FaceProxy *faceProxy = new FaceProxy();
    faceProxy->Create(SEGCTION_LENGTH, numberOfEquation);

    return faceProxy;
}

FaceProxy * TurbSolverUnstr::CreateFaceProxyTurbulence(Grid *gridIn)
{
    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    FaceProxy *faceProxy = new FaceProxy();
    faceProxy->Create(SEGCTION_LENGTH, nTurbulenceEquation);

    return faceProxy;
}

TurbFaceValue * TurbSolverUnstr::CreateTurbFaceValue(Grid *gridIn)
{
    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    TurbFaceValue *facevar = new TurbFaceValue(nTurbulenceEquation, SEGCTION_LENGTH);
    return facevar;
}

void TurbSolverUnstr::ComputeQTurbNodeValue(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    int numberOfTotalCell = grid->GetNTotalCell();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotalNode = grid->GetNTotalNode();
    int numberOfTotalFace = grid->GetNTotalFace();

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    RDouble **qTurbulenceNode = reinterpret_cast<RDouble **>(grid->GetDataPtr("qTurbNode"));
    RDouble     **qTurbulence = reinterpret_cast<RDouble **>(grid->GetDataPtr("q_turb"));
    int   *nodeValueSliceTurb = reinterpret_cast<int      *>(grid->GetDataPtr("nodeValueSliceTurb"));

    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    int *face2node = grid->GetFace2Node();
    int *nodeNumberOfEachFace = grid->GetNodeNumberOfEachFace();
 
    CompositeUnstructBC *unstructBCRegion = grid->GetUnstructBCRegion();
    int *bcRegionIDofBCFace = unstructBCRegion->GetBCFaceInfo();

    for (int iNode = 0; iNode < numberOfTotalNode; ++ iNode)
    {
        nodeValueSliceTurb[iNode] = 0;

        for(int m = 0; m < nTurbulenceEquation; ++ m)
        {
            qTurbulenceNode[m][iNode] = 0; 
        }
    }

    //! Boundary faces
    int nodePosition = 0;
    for (int iFace = 0; iFace < numberOfBoundaryFace; ++ iFace)
    {
        UnstructBC *bcRegion = unstructBCRegion->GetBCRegion(bcRegionIDofBCFace[iFace]);
        int bcType = bcRegion->GetBCType();

        int le = leftCellOfFace [iFace];
        int re = rightCellOfFace[iFace];
        if (bcType != PHENGLEI::INTERFACE)
        {
            for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
            {
                //! From left
                int point = face2node[nodePosition + jNode];

                for (int m = 0; m < nTurbulenceEquation; ++ m)
                {
                    qTurbulenceNode[m][point] += qTurbulence[m][le];
                }            
                nodeValueSliceTurb[point] += 1;

                //! From right
                for (int m = 0; m < nTurbulenceEquation; ++ m)
                {
                    qTurbulenceNode[m][point] += qTurbulence[m][re];
                }
                nodeValueSliceTurb[point] += 1;
            }
        }
        else
        {
            for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
            {
                int point = face2node[nodePosition + jNode];
                for (int m = 0; m < nTurbulenceEquation; ++ m)
                {
                    qTurbulenceNode[m][point] += qTurbulence[m][le];
                }
                nodeValueSliceTurb[point] += 1;
            }
        }        
        nodePosition += nodeNumberOfEachFace[iFace];
    }

    //! Interior faces
    for (int iFace = numberOfBoundaryFace; iFace < numberOfTotalFace; ++ iFace)
    {
        int le  = leftCellOfFace [iFace];
        int re  = rightCellOfFace[iFace];
        for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
        {
            //! From left
            int point = face2node[nodePosition + jNode];
            for (int m = 0; m < nTurbulenceEquation; ++ m)
            {
                qTurbulenceNode[m][point] += qTurbulence[m][le];
            }            
            ++ nodeValueSliceTurb[point];

            for (int m = 0; m < nTurbulenceEquation; ++ m)
            {
                qTurbulenceNode[m][point] += qTurbulence[m][re];
            }
            ++ nodeValueSliceTurb[point];
        }
        nodePosition += nodeNumberOfEachFace[iFace];
    }

    nodePosition = 0;
    for (int iFace = 0; iFace < numberOfBoundaryFace; ++ iFace)
    {
        UnstructBC *bcRegion = unstructBCRegion->GetBCRegion(bcRegionIDofBCFace[iFace]);
        int bcType = bcRegion->GetBCType();

        if (bcType == PHENGLEI::SYMMETRY)
        {
            for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
            {
                int point = face2node[nodePosition + jNode];
            
                for (int m = 0; m < nTurbulenceEquation; ++ m)
                {
                    qTurbulenceNode[m][point] = 0;
                }            
                nodeValueSliceTurb[point] = 0;
            }
        }
        nodePosition += nodeNumberOfEachFace[iFace];
    }

    nodePosition = 0;
    for (int iFace = 0; iFace < numberOfBoundaryFace; ++ iFace)
    {
        int le = leftCellOfFace[iFace];
        int re = iFace + numberOfTotalCell;
        
        UnstructBC *bcRegion = unstructBCRegion->GetBCRegion(bcRegionIDofBCFace[iFace]);
        int bcType = bcRegion->GetBCType();

        if (bcType == PHENGLEI::SYMMETRY)
        {
            for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
            {
                int point = face2node[nodePosition + jNode];
            
                for (int m = 0; m < nTurbulenceEquation; ++ m)
                {
                    qTurbulenceNode[m][point] += qTurbulence[m][le];
                }            
                ++ nodeValueSliceTurb[point];

                for (int m = 0; m < nTurbulenceEquation; ++ m)
                {
                    qTurbulenceNode[m][point] += qTurbulence[m][re];
                }
                ++ nodeValueSliceTurb[point];
            }
        }
        nodePosition += nodeNumberOfEachFace[iFace];
    }

    nodePosition = 0;
    for (int iFace = 0; iFace < numberOfBoundaryFace; ++ iFace)
    {
        UnstructBC *bcRegion = unstructBCRegion->GetBCRegion(bcRegionIDofBCFace[iFace]);
        int bcType = bcRegion->GetBCType();

        if (bcType == PHENGLEI::FARFIELD)
        {
            for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
            {
                int point = face2node[nodePosition + jNode];
            
                for (int m = 0; m < nTurbulenceEquation; ++ m)
                {
                    qTurbulenceNode[m][point] = 0;
                }            
                nodeValueSliceTurb[point] = 0;
            }
        }
        nodePosition += nodeNumberOfEachFace[iFace];
    }

    nodePosition = 0;
    for (int iFace = 0; iFace < numberOfBoundaryFace; ++ iFace)
    {
        int le = leftCellOfFace[iFace];
        int re = iFace + numberOfTotalCell;
       
        UnstructBC *bcRegion = unstructBCRegion->GetBCRegion(bcRegionIDofBCFace[iFace]);
        int bcType = bcRegion->GetBCType();

        if (bcType == PHENGLEI::FARFIELD)
        {
            for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
            {
                int point = face2node[nodePosition + jNode];
            
                for (int m = 0; m < nTurbulenceEquation; ++ m)
                {
                    qTurbulenceNode[m][point] += qTurbulence[m][le];
                }            
                ++ nodeValueSliceTurb[point];

                for (int m = 0; m < nTurbulenceEquation; ++ m)
                {
                    qTurbulenceNode[m][point] += qTurbulence[m][re];
                }
                ++ nodeValueSliceTurb[point];
            }
        }
        nodePosition += nodeNumberOfEachFace[iFace];
    }

    nodePosition = 0;
    for (int iFace = 0; iFace < numberOfBoundaryFace; ++ iFace)
    {
        UnstructBC *bcRegion = unstructBCRegion->GetBCRegion(bcRegionIDofBCFace[iFace]);
        int bcType = bcRegion->GetBCType();

        if (bcType == PHENGLEI::SOLID_SURFACE)
        {
            for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
            {
                int point = face2node[nodePosition + jNode];
            
                for (int m = 0; m < nTurbulenceEquation; ++ m)
                {
                    qTurbulenceNode[m][point] = 0;
                }            
                nodeValueSliceTurb[point] = 0;
            }
        }
        nodePosition += nodeNumberOfEachFace[iFace];
    }

    nodePosition = 0;
    for (int iFace = 0; iFace < numberOfBoundaryFace; ++ iFace)
    {
        int le = leftCellOfFace[iFace];
        int re = iFace + numberOfTotalCell;
        
        UnstructBC *bcRegion = unstructBCRegion->GetBCRegion(bcRegionIDofBCFace[iFace]);
        int bcType = bcRegion->GetBCType();

        if (bcType == PHENGLEI::SOLID_SURFACE)
        {
            for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
            {
                int point = face2node[nodePosition + jNode];
            
                for (int m = 0; m < nTurbulenceEquation; ++ m)
                {
                    qTurbulenceNode[m][point] += qTurbulence[m][le];
                }            
                ++ nodeValueSliceTurb[point];

                for (int m = 0; m < nTurbulenceEquation; ++ m)
                {
                    qTurbulenceNode[m][point] += qTurbulence[m][re];
                }
                ++ nodeValueSliceTurb[point];
            }
        }
        nodePosition += nodeNumberOfEachFace[iFace];
    }

    //! Change the number of points divided by the interpoint to the number of points before partition,because it is a face loop,each point is calculated twice,so the number of points needs to be multiplied by 2.
    InterpointInformation *interPointInfor = grid->GetInterpointInfo();
    if (interPointInfor)
    {
        int numberOfInterPoint = interPointInfor->GetNumberOfInterpoints();
        int *labelOfInterPoint = interPointInfor->GetLabelOfInterPoint();
        int *interPoint2GlobalPoint = interPointInfor->GetInterPoint2GlobalPoint();
        
        int globalPoint = 0;
        for (int iPoint = 0; iPoint < numberOfInterPoint; ++ iPoint)
        {
            globalPoint = interPoint2GlobalPoint[iPoint];
            if (labelOfInterPoint[iPoint] != 1)
            {
                for (int m = 0; m < nTurbulenceEquation;++ m)
                {    
                    globalPoint = interPoint2GlobalPoint[iPoint];
                    qTurbulenceNode[m][globalPoint] = 0;
                }
            }
        }
        RDouble **qTurbulenceInterPoint = reinterpret_cast<RDouble **> (grid->GetDataPtr("turb::qInterpoint"));
        for (int iPoint = 0; iPoint < numberOfInterPoint; ++ iPoint)
        {
            for (int m = 0; m < nTurbulenceEquation; ++ m)
            {
                qTurbulenceInterPoint[m][iPoint] = 0;
            }
        }
    }
}

void TurbSolverUnstr::ModifyQTurbNodeValue(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalNode = grid->GetNTotalNode();

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    RDouble **qTurbulenceNode = reinterpret_cast<RDouble **>(grid->GetDataPtr("qTurbNode"));
    int   *nodeValueSliceTurb = reinterpret_cast<int      *>(grid->GetDataPtr("nodeValueSliceTurb"));

    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        if (nodeValueSliceTurb[iNode] != 0)
        {
            for (int m = 0; m < nTurbulenceEquation; m++)
            {
                qTurbulenceNode[m][iNode] /= nodeValueSliceTurb[iNode];
            }
        }
        else
        {
            for (int m = 0; m < nTurbulenceEquation; m++)
            {
                qTurbulenceNode[m][iNode] = 0.0;
            }
        }
    }
}

void TurbSolverUnstr::ComputeFaceWeight(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    RDouble *deltL = faceProxy->GetWeightL();
    RDouble *deltR = faceProxy->GetWeightR();

    grid->FaceWeight(deltL, deltR, localStart, localEnd);
}

void TurbSolverUnstr::ViscousFlux(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalFace = grid->GetNTotalFace();

    FaceProxy *faceProxy = CreateFaceProxy(grid);
    faceProxy->SetTurbFaceValue(CreateTurbFaceValue(grid));

    int localStart,localEnd;
    localStart = 0;
    do
    {
        localEnd = localStart + SEGCTION_LENGTH;
        if (localEnd > numberOfTotalFace) localEnd = numberOfTotalFace;

        ComputeFaceWeight(grid, faceProxy, localStart, localEnd);

        GetVisFaceValue(grid, faceProxy, localStart, localEnd);

        ComputeVisflux(grid, faceProxy, localStart, localEnd);

        //LoadFlux(grid, faceProxy, localStart, localEnd);

        localStart = localEnd;
    } while (localStart < numberOfTotalFace);

    delete faceProxy;
}

void TurbSolverUnstr::GetVisFaceValue(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotalCell = grid->GetNTotalCell();

    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    TurbFaceValue *faceVariable = faceProxy->GetTurbFaceValue();

    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));

    RDouble *viscousLaminar = reinterpret_cast<RDouble *> (grid->GetDataPtr("visl"));
    RDouble *viscousTurbulence = reinterpret_cast<RDouble *> (grid->GetDataPtr("vist"));

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int viscoueType = parameters->GetViscousType();

    string viscousName = parameters->GetViscousName();

    RDouble *mul   = faceVariable->mul;
    RDouble *mut   = faceVariable->mut;
    RDouble **mlt  = faceVariable->mlt;

    int le,re,jFace;
    for (int iFace = localStart; iFace <localEnd; ++ iFace)
    {
        le = leftCellOfFace [iFace];
        re = rightCellOfFace[iFace];
        if (iFace < numberOfBoundaryFace) re = iFace + numberOfTotalCell;
        jFace  = iFace - localStart;

        mul[jFace] = half * (viscousLaminar[le] + viscousLaminar[re]);
        mut[jFace] = half * (viscousTurbulence[le] + viscousTurbulence[re]);
    }

    CompositeUnstructBC *unstructBCRegion = grid->GetUnstructBCRegion();
    int *bcRegionIDofBCFace = unstructBCRegion->GetBCFaceInfo();

    int nMid = numberOfBoundaryFace;
    if (localEnd < numberOfBoundaryFace) nMid = localEnd;

    if (localStart < numberOfBoundaryFace)
    {
        for (int iFace = localStart; iFace < nMid; ++ iFace)
        {
            UnstructBC *bcRegion = unstructBCRegion->GetBCRegion(bcRegionIDofBCFace[iFace]);
            int bcType = bcRegion->GetBCType();

            int jFace  = iFace - localStart;

            le = leftCellOfFace[iFace];

            if (IsInterface(bcType))
            {
                continue;
            }

            if (IsWall(bcType))
            {
                mul[jFace] = viscousLaminar[le];
                mut[jFace] = 0.0;
            }
            else
            {
                mul[jFace] = viscousLaminar[le];
                mut[jFace] = viscousTurbulence[le];
            }
        }
    }

    using namespace IDX;

    if (viscoueType == 3)
    {
        RDouble KW_sigma = parameters->GetKW_sigma();
        RDouble oSigma = 1.0 / KW_sigma;
        RDouble **qPrimitiveVariable = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));

        for (int iFace = localStart; iFace <localEnd; ++ iFace)
        {
            le = leftCellOfFace [iFace];
            re = rightCellOfFace[iFace];
            if (iFace < numberOfBoundaryFace) re = iFace + numberOfTotalCell;
            jFace  = iFace - localStart;

            RDouble oRho = two / (qPrimitiveVariable[IR][le] + qPrimitiveVariable[IR][re] + SMALL);
            mlt[ISA][jFace] = oSigma * (mul[jFace] * oRho + half * (qTurbulence[ISA][le] + qTurbulence[ISA][re]));
        }
    }
    else
    {
        if (viscousName.substr(0, 13) == "2eq-kw-menter" || viscousName.substr(0, 12) == "easm-kw-2005")
        {
            RDouble *blend = reinterpret_cast<RDouble  *> (grid->GetDataPtr("blend"));

            RDouble turbSigk,turbSigw;

            RDouble KW_sigmaK1 = parameters->GetKW_sigmaK1();
            RDouble KW_sigmaW1 = parameters->GetKW_sigmaW1();
            RDouble KW_sigmaK2 = parameters->GetKW_sigmaK2();
            RDouble KW_sigmaW2 = parameters->GetKW_sigmaW2();

            int transitionType = parameters->GetTransitionType();

            RDouble transition_dct = 0.0, transition_df = 0.0;
            if (transitionType == IREGAMA)
            {
                transition_dct = parameters->GetTransition_dct();
                transition_df = parameters->GetTransition_df();
            }

            for (int iFace = localStart; iFace <localEnd; ++ iFace)
            {
                le = leftCellOfFace [iFace];
                re = rightCellOfFace[iFace];
                if (iFace < numberOfBoundaryFace) re = iFace + numberOfTotalCell;
                jFace  = iFace - localStart;

                //RDouble cblend = half * (blend[le] + blend[re]);
                RDouble cblend = blend[le];
                turbSigk   = cblend * KW_sigmaK1 + (1.0 - cblend) * KW_sigmaK2;
                turbSigw   = cblend * KW_sigmaW1 + (1.0 - cblend) * KW_sigmaW2;

                mlt[IKE][jFace] = mul[jFace] + mut[jFace] * turbSigk;
                mlt[IKW][jFace] = mul[jFace] + mut[jFace] * turbSigw;

                //Re-gama transition model
                if (transitionType == IREGAMA)
                {
                    mlt[IGAMA][jFace] = mul[jFace] + mut[jFace] / transition_df;
                    mlt[IRECT][jFace] = (mul[jFace] + mut[jFace]) * transition_dct;
                }
            }
        }
        else
        {
            RDouble turbSigk,turbSigw;

            GlobalDataBase::GetData("turbSigk", &turbSigk, PHDOUBLE, 1);
            GlobalDataBase::GetData("turbSigw", &turbSigw, PHDOUBLE, 1);

            for (int iFace = localStart; iFace <localEnd; ++ iFace)
            {
                le = leftCellOfFace [iFace];
                re = rightCellOfFace[iFace];
                if (iFace < numberOfBoundaryFace) re = iFace + numberOfTotalCell;
                jFace  = iFace - localStart;

                mlt[IKE][jFace] = mul[jFace] + mut[jFace] * turbSigk;
                mlt[IKW][jFace] = mul[jFace] + mut[jFace] * turbSigw;
            }
        }
    }
}

//! Bell 20130513 mod.
void TurbSolverUnstr::ComputeVisflux(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    int numberOfBoundaryFace = grid->GetNBoundFace();

    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    CompositeUnstructBC *unstructBCRegion = grid->GetUnstructBCRegion();
    int *bcRegionIDofBCFace = unstructBCRegion->GetBCFaceInfo();

    RDouble **flux  = faceProxy->GetFlux();
    TurbFaceValue *faceVariable = faceProxy->GetTurbFaceValue();
    
    RDouble *xfn  = grid->GetFaceNormalX();
    RDouble *yfn  = grid->GetFaceNormalY();
    RDouble *zfn  = grid->GetFaceNormalZ();
    RDouble *area = grid->GetFaceArea();
    
    RDouble *cellCenterX = grid->GetCellCenterX();
    RDouble *cellCenterY = grid->GetCellCenterY();
    RDouble *cellCenterZ = grid->GetCellCenterZ();

    RDouble *xfc  = grid->GetFaceCenterX();
    RDouble *yfc  = grid->GetFaceCenterY();
    RDouble *zfc  = grid->GetFaceCenterZ();

    RDouble **mlt  = faceVariable->mlt;

    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));

    RDouble *viscousLaminar = reinterpret_cast<RDouble *> (grid->GetDataPtr("visl"));
    RDouble **qPrimitiveVariable = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));

    RDouble **residualTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));

    Param_TurbSolverUnstruct *parameters = GetControlParameters();

    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();
    RDouble oRefReNumber = parameters->GetoRefReNumber();

    const RDouble cb2    = 0.622;
    const RDouble sigma  = 2.0/3.0;
    const RDouble osigma = 1.0 / sigma;

    RDouble *dfd1  = new RDouble[nTurbulenceEquation];
    RDouble *dfd2  = new RDouble[nTurbulenceEquation];
    RDouble *dfdn  = new RDouble[nTurbulenceEquation];

    RDouble *f1    = new RDouble[nTurbulenceEquation];
    RDouble *f2    = new RDouble[nTurbulenceEquation];
    RDouble *fMid  = new RDouble[nTurbulenceEquation];

    int le, re, jFace;
    for (int iFace = localStart; iFace <localEnd; ++ iFace)
    {
        le = leftCellOfFace [iFace];
        re = rightCellOfFace[iFace];
        jFace  = iFace - localStart;

        RDouble nxs = xfn[iFace];
        RDouble nys = yfn[iFace];
        RDouble nzs = zfn[iFace];       

        RDouble dxL = cellCenterX[le] - xfc[iFace];
        RDouble dyL = cellCenterY[le] - yfc[iFace];
        RDouble dzL = cellCenterZ[le] - zfc[iFace];

        RDouble dxR = cellCenterX[re] - xfc[iFace];
        RDouble dyR = cellCenterY[re] - yfc[iFace];
        RDouble dzR = cellCenterZ[re] - zfc[iFace];

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

        //! quantities at points 1 and 2
        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            f1[m] = qTurbulence[m][le];
            f2[m] = qTurbulence[m][re];
        }

        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            fMid[m] = half * (f1[m] + f2[m]);
        }

        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            dfdn[m] = 0.0;
        }

        if (angle1 > 0.0 && angle2 > 0.0 && ABS(dL) > TINY && ABS(dR) > TINY)
        {
            for (int m = 0; m < nTurbulenceEquation; ++ m)
            {
                dfd1[m] = (f1[m] - fMid[m]) / dL;
            }

            for (int m = 0; m < nTurbulenceEquation; ++ m)
            {
                dfd2[m] = (f2[m] - fMid[m]) / dR;
            }

            RDouble dtmp = dL * dL + dR * dR;
            RDouble weightL = dL * dL / dtmp;
            RDouble weightR = 1.0 - weightL;

            for (int m = 0; m < nTurbulenceEquation; ++ m)
            {
                dfdn[m] = weightL * dfd1[m] + weightR * dfd2[m];
            }

            if (iFace < numberOfBoundaryFace)
            {
                UnstructBC *bcRegion = unstructBCRegion->GetBCRegion(bcRegionIDofBCFace[iFace]);
                int bcType = bcRegion->GetBCType();

                if (bcType != PHENGLEI::INTERFACE && bcType != PHENGLEI::SYMMETRY)
                {
                    RDouble normalDist = (cellCenterX[re] - cellCenterX[le]) * nxs +
                                         (cellCenterY[re] - cellCenterY[le]) * nys +
                                         (cellCenterZ[re] - cellCenterZ[le]) * nzs;

                    for (int m = 0; m < nTurbulenceEquation; ++ m)
                    {
                        dfdn[m] = (f2[m] - f1[m]) / normalDist;
                    }
                }
            }
        }
        
        if(TurbEquationFormation() == FORMATION_SA_CFL3D)
        {
            int m = IDX::ISA;

            RDouble viscousFace = half * (viscousLaminar[le] + viscousLaminar[re]);
            RDouble rhoFace     = half * (qPrimitiveVariable[IDX::IR][le] + qPrimitiveVariable[IDX::IR][re]);
            RDouble nueFace     = viscousFace / rhoFace;
            RDouble nueSAFace   = half * (qTurbulence[m][le] + qTurbulence[m][re]);

            //! Term1 && term2.           
            RDouble flux1 = - oRefReNumber * osigma * (nueFace + (1.0 + cb2) * nueSAFace) * dfdn[m] * area[iFace];
            RDouble flux2 =   oRefReNumber * osigma * cb2 * dfdn[m] * area[iFace];

            residualTurbulence[m][le] -= flux1;
            residualTurbulence[m][le] -= flux2 * qTurbulence[IDX::ISA][le];

            if(iFace >= numberOfBoundaryFace)
            {
                residualTurbulence[m][re] += flux1;
                residualTurbulence[m][re] += flux2 * qTurbulence[IDX::ISA][re];
            }
        }
        else if(TurbEquationFormation() == FORMATION_SA_NSMB)
        {     
            int m = IDX::ISA;

            RDouble flux1 = - oRefReNumber * (1.0 + cb2) * mlt[m][jFace] * dfdn[m] * area[iFace];
            RDouble flux2 =  oRefReNumber * cb2 * osigma * dfdn[m] * area[iFace];
            
            residualTurbulence[m][le] -= flux1;
            residualTurbulence[m][le] -= flux2 * (viscousLaminar[le] / qPrimitiveVariable[IDX::IR][le] + qTurbulence[IDX::ISA][le]);

            if(iFace >= numberOfBoundaryFace)
            {
                residualTurbulence[m][re] += flux1;
                residualTurbulence[m][re] += flux2 * (viscousLaminar[re] / qPrimitiveVariable[IDX::IR][re] + qTurbulence[IDX::ISA][re]);
            }
        }
        else if(TurbEquationFormation() == FORMATION_ORGSA_SST)
        {
            for (int m = 0; m < nTurbulenceEquation; ++ m)
            {
                flux[m][jFace] = - oRefReNumber * mlt[m][jFace] * dfdn[m] * area[iFace];
            }
            
            for (int m = 0; m < nTurbulenceEquation; ++ m)
            {
                residualTurbulence[m][le] -= flux[m][jFace];

                if(iFace >= numberOfBoundaryFace)
                {
                    residualTurbulence[m][re] += flux[m][jFace];
                }
            }
        }        
    }

    delete []dfd1;
    delete []dfd2;
    delete []dfdn;

    delete []f1;
    delete []f2;
    delete []fMid;
}

void TurbSolverUnstr::Diagonal(Grid *grid)
{
    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int viscoueType = parameters->GetViscousType();

    if (viscoueType == 3)
    {
        SpectrumRadiusOfOneEquation(grid);
    }
    else if (viscoueType == 4)
    {
        SpectrumRadiusOfTwoEquation(grid);
    }
}

void TurbSolverUnstr::SpectrumRadiusOfOneEquation(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalFace    = grid->GetNTotalFace();

    RDouble **q         = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
    RDouble *viscousLaminar = reinterpret_cast<RDouble *> (grid->GetDataPtr("visl"));

    Param_TurbSolverUnstruct *parameters = GetControlParameters();

    RDouble **spectrumTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("spec_turb"));
    RDouble **matrixTurbulenceLeft = reinterpret_cast<RDouble **> (grid->GetDataPtr("mat_turbl"));
    RDouble **matrixTurbulenceRight = reinterpret_cast<RDouble **> (grid->GetDataPtr("mat_turbr"));

    RDouble oRefReNumber = parameters->GetoRefReNumber();

    //!initialize some constants
    RDouble sigma  = 2.0 / 3.0;
    RDouble osigma = one / sigma;
    const RDouble cb2 = 0.622;

    int *leftCellOfFace  = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    RDouble *vol  = grid->GetCellVolume();
    RDouble *xfn  = grid->GetFaceNormalX();
    RDouble *yfn  = grid->GetFaceNormalY();
    RDouble *zfn  = grid->GetFaceNormalZ();
    RDouble *area = grid->GetFaceArea();

    using namespace IDX;

    RDouble invJacobianL, invJacobianR;
    for (int iFace = 0; iFace < numberOfTotalFace; ++ iFace)
    {
        int le = leftCellOfFace [iFace];
        int re = rightCellOfFace[iFace];

        //! Consider vgn here!

        //! Inviscid Jacobian.
        RDouble rhoL = q[IR][le];
        RDouble uL = q[IU][le];
        RDouble vL = q[IV][le];
        RDouble wL = q[IW][le];
        RDouble vnL = xfn[iFace] * uL + yfn[iFace] * vL + zfn[iFace] * wL;

        if (vnL > 0.0)
        {
            invJacobianL = rhoL * vnL * area[iFace];
        }
        else
        {
            invJacobianL = 0.0;
        }

        RDouble rhoR = q[IR][re];
        RDouble uR = q[IU][re];
        RDouble vR = q[IV][re];
        RDouble wR = q[IW][re];
        RDouble vnR = xfn[iFace] * uR + yfn[iFace] * vR + zfn[iFace] * wR;
        if (vnR < 0.0)
        {
            invJacobianR = rhoR * vnR * area[iFace];
        }
        else
        {
            invJacobianR = 0.0;
        }

        spectrumTurbulence[ISA][le] += invJacobianL;
        spectrumTurbulence[ISA][re] -= invJacobianR;

        matrixTurbulenceRight[ISA][iFace] = invJacobianR;
        matrixTurbulenceLeft[ISA][iFace]  = -invJacobianL;
        
        //! Viscous Jacobian.
        RDouble orL = 1.0 / (rhoL + SMALL);
        RDouble orR = 1.0 / (rhoR + SMALL);

        RDouble nueL = qTurbulence[ISA][le];
        RDouble nueR = qTurbulence[ISA][re];

        RDouble nueffL = viscousLaminar[le] * orL + ABS(nueL);
        RDouble nueffR = viscousLaminar[re] * orR + ABS(nueR);

        RDouble nueff   = half * (nueffL + nueffR);

        RDouble volume  = half * (vol[le] + vol[re]);
        RDouble oVolume = one / volume;

        RDouble ns2     = area[iFace] * area[iFace];
        RDouble lengthScale = ns2 * oVolume;

        if (TurbEquationFormation() == FORMATION_SA_CFL3D)
        {    
            RDouble viscousFace = half * (viscousLaminar[le] + viscousLaminar[re]);
            RDouble rhoFace     = half * (q[IDX::IR][le] + q[IDX::IR][re]);
            RDouble nueFace     = viscousFace / rhoFace;
            RDouble nueSAFace   = half * (qTurbulence[IDX::ISA][le] + qTurbulence[IDX::ISA][re]);         
            
            //! Term1 && term2.           
            RDouble term1 =   oRefReNumber * osigma * (nueFace + (1.0 + cb2) * nueSAFace) * lengthScale;
            RDouble term2 = - oRefReNumber * osigma * lengthScale;

            RDouble jocabianL = (term1 + term2 * nueL);
            RDouble jocabianR = (term1 + term2 * nueR);

            spectrumTurbulence[ISA][le] += jocabianL;
            spectrumTurbulence[ISA][re] += jocabianR;

            matrixTurbulenceRight[ISA][iFace] -= jocabianR;  // Right cell of iFace.
            matrixTurbulenceLeft[ISA][iFace]  -= jocabianL;  // Left cell of iFace.
        }
        else if (TurbEquationFormation() == FORMATION_SA_NSMB)
        {
            spectrumTurbulence[ISA][le] += oRefReNumber * oVolume * (1.0 + cb2) * (osigma * nueff * ns2);
            spectrumTurbulence[ISA][re] += oRefReNumber * oVolume * (1.0 + cb2) * (osigma * nueff * ns2);

            matrixTurbulenceLeft[ISA][iFace]  += - oRefReNumber * oVolume * (1.0 + cb2) * (osigma * nueff * ns2);
            matrixTurbulenceRight[ISA][iFace]  += - oRefReNumber * oVolume * (1.0 + cb2) * (osigma * nueff * ns2);

            //! 下面的nsmb第2项貌似也应该加到spectrumTurbulence中
            RDouble coeffLe, coeffRe;
            coeffLe = cb2 * osigma * nueffL;
            coeffRe = cb2 * osigma * nueffR;
            matrixTurbulenceLeft[ISA][iFace] += oRefReNumber * oVolume * ns2 * coeffLe;
            matrixTurbulenceRight[ISA][iFace] += oRefReNumber * oVolume * ns2 * coeffRe;
        } 
        else if (TurbEquationFormation() == FORMATION_ORGSA_SST)
        { 
            spectrumTurbulence[ISA][le] += oRefReNumber * oVolume * (osigma * nueff * ns2);
            spectrumTurbulence[ISA][re] += oRefReNumber * oVolume * (osigma * nueff * ns2);

            matrixTurbulenceLeft[ISA][iFace] += - oRefReNumber * oVolume * (osigma * nueff * ns2);
            matrixTurbulenceRight[ISA][iFace] += - oRefReNumber * oVolume * (osigma * nueff * ns2);
        }  
        //*/
    }
}

void TurbSolverUnstr::SpectrumRadiusOfTwoEquation(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalFace = grid->GetNTotalFace();

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    RDouble oRefReNumber = parameters->GetoRefReNumber();

    RDouble **q = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble *viscousLaminar = reinterpret_cast<RDouble * > (grid->GetDataPtr("visl"));
    RDouble *viscousTurbulence = reinterpret_cast<RDouble * > (grid->GetDataPtr("vist"));
    
    RDouble **spectrumTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("spec_turb"));
    RDouble **matrixTurbulenceLeft = reinterpret_cast<RDouble **> (grid->GetDataPtr("mat_turbl"));
    RDouble **matrixTurbulenceRight = reinterpret_cast<RDouble **> (grid->GetDataPtr("mat_turbr"));

    RDouble *vol  = grid->GetCellVolume();
    RDouble *xfn  = grid->GetFaceNormalX();
    RDouble *yfn  = grid->GetFaceNormalY();
    RDouble *zfn  = grid->GetFaceNormalZ();
    RDouble *area = grid->GetFaceArea();

    RDouble uL, vL, wL, vnL, absVnL;
    RDouble uR, vR, wR, vnR, absVnR;

    RDouble vn, absVn, ns2, rho, mul, mut, cblend, volume, oVolume;

    int *leftCellOfFace  = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    vector<RDouble> work(nTurbulenceEquation);

    using namespace IDX;

    for (int iFace = 0; iFace < numberOfTotalFace; ++ iFace)
    {
        int le = leftCellOfFace [iFace];
        int re = rightCellOfFace[iFace];

        //! Calculate velocity at the cell interface.
        uL = q[IU][le];
        vL = q[IV][le];
        wL = q[IW][le];
        vnL = xfn[iFace] * uL + yfn[iFace] * vL + zfn[iFace] * wL;

        absVnL = ABS(vnL);

        uR = q[IU][re];
        vR = q[IV][re];
        wR = q[IW][re];
        vnR = xfn[iFace] * uR + yfn[iFace] * vR + zfn[iFace] * wR;

        absVnR = ABS(vnR);

        vn = half * (vnL + vnR);
        absVn = ABS(vn);

        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            spectrumTurbulence[m][le] += area[iFace] * half * absVn;
            spectrumTurbulence[m][re] += area[iFace] * half * absVn;

            matrixTurbulenceRight[m][iFace]  = area[iFace] * half * ( vnR - absVn);
            matrixTurbulenceLeft[m][iFace]   = area[iFace] * half * (- vnL - absVn);
        }
    }

    RDouble *blend = reinterpret_cast<RDouble  *> (grid->GetDataPtr("blend"));

    RDouble turbSigk, turbSigw;

    RDouble KW_sigmaK1 = parameters->GetKW_sigmaK1();
    RDouble KW_sigmaW1 = parameters->GetKW_sigmaW1();
    RDouble KW_sigmaK2 = parameters->GetKW_sigmaK2();
    RDouble KW_sigmaW2 = parameters->GetKW_sigmaW2();

    int transitionType = parameters->GetTransitionType();

    RDouble transition_dct, transition_df;
    if (transitionType == IREGAMA)
    {
        transition_dct = parameters->GetTransition_dct();
        transition_df = parameters->GetTransition_df();
    }

    for (int iFace = 0; iFace < numberOfTotalFace; ++ iFace)
    {
        int le = leftCellOfFace [iFace];
        int re = rightCellOfFace[iFace];
        rho = half * (q[IR][le] + q[IR][re]);
        mul = half * (viscousLaminar[le] + viscousLaminar[re]);
        mut = half * (viscousTurbulence[le] + viscousTurbulence[re]);

        cblend    = half * (blend[le] + blend[re]);
        turbSigk = cblend * KW_sigmaK1 + (1.0 - cblend) * KW_sigmaK2;
        turbSigw = cblend * KW_sigmaW1 + (1.0 - cblend) * KW_sigmaW2;

        volume  = half * (vol[le] + vol[re]);
        oVolume = one / volume;

        ns2     = area[iFace] * area[iFace];
        work[IKE] = (mul + mut * turbSigk) / rho * ns2 * oVolume;
        work[IKW] = (mul + mut * turbSigw) / rho * ns2 * oVolume;

        if (transitionType == IREGAMA)
        {
              work[IGAMA] = (mul + mut / transition_df) / rho  * ns2 * oVolume;
              work[IRECT] = (mul + mut) * transition_dct / rho * ns2 * oVolume;
        }

        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            spectrumTurbulence[m][le] += oRefReNumber * work[m];
            spectrumTurbulence[m][re] += oRefReNumber * work[m];

            matrixTurbulenceLeft[m][iFace]  += - oRefReNumber * work[m];
            matrixTurbulenceRight[m][iFace] += - oRefReNumber * work[m];
        }
    }
}

void TurbSolverUnstr::GetQlQr(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    GetQlQrTurbulence(grid, faceProxy, localStart, localEnd);

    FaceProxy *faceProxyNS = faceProxy->GetNext();
    GetQlQrNS(grid, faceProxyNS, localStart, localEnd);

}

void TurbSolverUnstr::GetQlQrTurbulence(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    RDouble **qTurbulenceVarables = reinterpret_cast <RDouble **> (grid->GetDataPtr("q_turb"));

    RDouble **qL = faceProxy->GetQL();
    RDouble **qR = faceProxy->GetQR();

    for (int iFace = localStart; iFace <localEnd; ++ iFace)
    {
        int le, re, jFace;
        le = leftCellOfFace [iFace];
        re = rightCellOfFace[iFace];
        jFace  = iFace - localStart;

        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            qL[m][jFace] = qTurbulenceVarables[m][le];
            qR[m][jFace] = qTurbulenceVarables[m][re];
        }
    }

    //! Boundary fix.
    CompositeUnstructBC *unstructBCRegion = grid->GetUnstructBCRegion();
    int *bcRegionIDofBCFace = unstructBCRegion->GetBCFaceInfo();

    int nBoundFace = grid->GetNBoundFace();
    int nMid = 0;

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

    for (int iFace = localStart; iFace < nMid; ++ iFace)
    {
        UnstructBC *bcRegion = unstructBCRegion->GetBCRegion(bcRegionIDofBCFace[iFace]);
        int bcType = bcRegion->GetBCType();

        int jFace = iFace - localStart;

        if (bcType == PHENGLEI::INTERFACE)
        {
            continue;
        }

        if (bcType == PHENGLEI::SYMMETRY)
        {
            for (int m = 0; m < nTurbulenceEquation; ++ m)
            {
                qR[m][jFace] = qL[m][jFace];
            }
        }
        else
        {
            for (int m = 0; m < nTurbulenceEquation; ++ m)
            {
                RDouble temp = half * (qR[m][jFace] + qL[m][jFace]);
                qR[m][jFace] = temp;
                qL[m][jFace] = temp;
            }            
        }
    }
}

void TurbSolverUnstr::GetQlQrNS(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    int nl    = GlobalDataBase::GetIntParaFromDB("nl");
    int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
    int nEquation = nl + nchem;
   
    RDouble **q = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));

    RDouble **qL = faceProxy->GetQL();
    RDouble **qR = faceProxy->GetQR();

    for (int iFace = localStart; iFace <localEnd; ++ iFace)
    {
        int le, re, jFace;
        le = leftCellOfFace [iFace];
        re = rightCellOfFace[iFace];
        jFace  = iFace - localStart;

        for (int m = 0; m < nEquation; ++ m)
        {
            qL[m][jFace] = q[m][le];
            qR[m][jFace] = q[m][re];
        }
    }

    //! Boundary fix.
    CompositeUnstructBC *unstructBCRegion = grid->GetUnstructBCRegion();
    int *bcRegionIDofBCFace = unstructBCRegion->GetBCFaceInfo();
    
    int nBoundFace = grid->GetNBoundFace();
    int nMid = 0;

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

    for (int iFace = localStart; iFace < nMid; ++ iFace)
    {
        UnstructBC *bcRegion = unstructBCRegion->GetBCRegion(bcRegionIDofBCFace[iFace]);
        int bcType = bcRegion->GetBCType();
       
        int jFace = iFace - localStart;

        if (bcType == PHENGLEI::INTERFACE)
        {
            continue;
        }

        for (int m = 0; m < nEquation; ++ m)
        {
            RDouble temp = half * (qR[m][jFace] + qL[m][jFace]);
            qR[m][jFace] = temp;
            qL[m][jFace] = temp;
        }
    }
}

void TurbSolverUnstr::InviscidFlux(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalFace = grid->GetNTotalFace();

    FaceProxy *faceProxy = CreateFaceProxy(grid);

    int localStart, localEnd;

    localStart = 0;
    do
    {
        localEnd = localStart + SEGCTION_LENGTH;
        if (localEnd > numberOfTotalFace) localEnd = numberOfTotalFace;

        GetQlQr(grid, faceProxy, localStart, localEnd);

        ComputeInviscidFlux(grid, faceProxy, localStart, localEnd);

        LoadFlux(grid, faceProxy, localStart, localEnd);

        localStart = localEnd;

    } while (localStart < numberOfTotalFace);

    delete faceProxy;
}

void TurbSolverUnstr::ComputeInviscidFlux(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    NNDFlux(grid, faceProxy, localStart, localEnd);
}

void TurbSolverUnstr::NNDFlux(Grid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    RDouble *xfn  = grid->GetFaceNormalX();
    RDouble *yfn  = grid->GetFaceNormalY();
    RDouble *zfn  = grid->GetFaceNormalZ();
    RDouble *vgn  = grid->GetFaceNormalVelocity();
    RDouble *area = grid->GetFaceArea();

    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();
    int viscoueType = parameters->GetViscousType();

    RDouble **qL   = faceProxy->GetQL();
    RDouble **qR   = faceProxy->GetQR();
    RDouble **flux = faceProxy->GetFlux();

    FaceProxy *faceProxyNS = faceProxy->GetNext();

    RDouble **qPrimitiveVariableL = faceProxyNS->GetQL();
    RDouble **qPrimitiveVariableR = faceProxyNS->GetQR();

    RDouble *primitiveVariableL = new RDouble[nTurbulenceEquation];
    RDouble *primitiveVariableR = new RDouble[nTurbulenceEquation];
    RDouble *f = new RDouble[nTurbulenceEquation];

    RDouble rL,uL,vL,wL,vnL,vnLL,rR,uR,vR,wR,vnR,vnRR;
    RDouble rhoFlag = 0, rhoCoefficientL, rhoCoefficientR;

    using namespace IDX;

    if (viscoueType == 3)
    {
        //! rhoFlag, is used to compatible both original and conversative variables.
        //! For SA equation, using original variable; else, using conversative variables.
        rhoFlag = 1;
    }

    for (int iFace = localStart; iFace <localEnd; ++ iFace)
    {
        int le, re, jFace;
        jFace = iFace - localStart;
        le = leftCellOfFace [iFace];
        re = rightCellOfFace[iFace];

        rL = qPrimitiveVariableL[IR][jFace];
        uL = qPrimitiveVariableL[IU][jFace];
        vL = qPrimitiveVariableL[IV][jFace];
        wL = qPrimitiveVariableL[IW][jFace];

        rR = qPrimitiveVariableR[IR][jFace];
        uR = qPrimitiveVariableR[IU][jFace];
        vR = qPrimitiveVariableR[IV][jFace];
        wR = qPrimitiveVariableR[IW][jFace];

        vnL  = xfn[iFace] * uL + yfn[iFace] * vL + zfn[iFace] * wL - vgn[iFace];
        vnR  = xfn[iFace] * uR + yfn[iFace] * vR + zfn[iFace] * wR - vgn[iFace];

        vnLL = half * (vnL + ABS(vnL));
        vnRR = half * (vnR - ABS(vnR));

        rhoCoefficientL = rhoFlag * one + (one - rhoFlag) * rL;
        rhoCoefficientR = rhoFlag * one + (one - rhoFlag) * rR;

        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            primitiveVariableL[m] = rhoCoefficientL * qL[m][jFace];
            primitiveVariableR[m] = rhoCoefficientR * qR[m][jFace];
        }

        //! Consider vgn here!
        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            f[m] = (vnLL * primitiveVariableL[m] + vnRR * primitiveVariableR[m]);
        }

        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            flux[m][jFace] = f[m] * area[iFace];
        }
    }
   
    delete []primitiveVariableL;
    delete []primitiveVariableR;
    delete []f;
}

void TurbSolverUnstr::GetResidual(ActionKey *actkey)
{
    int level = actkey->level;
    UnstructGrid *grid = UnstructGridCast(GetGrid(level));

    RDouble **residualTurbulence  = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    PHSPACE::ComputeResidualonGrid(grid, actkey, residualTurbulence, nTurbulenceEquation);
}

void TurbSolverUnstr::CompressSpecificArrayToInterface(DataContainer *&dataContainer, const string &fieldName, Grid *gridIn, const int &neighborZoneIndex, const int &nEquation)
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
    if (0 == nEquation)
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

void TurbSolverUnstr::DecompressArrayFromInterface(DataContainer *&dataContainer, const string &fieldName, Grid *gridIn, const int &neighborZoneIndex, const int &nEquation)
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

    if (0 == nEquation)
    {
        RDouble *fieldReceive = reinterpret_cast <RDouble *> (grid->GetDataPtr(fieldName));

        for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; -- iGhostLayer)
        {
            for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
            {
                int targetCell;
                int iFace = interfaceIndexContainerForReceive[iLocalFace];
                grid->GetTargetIndex(iFace, iGhostLayer + 1, targetCell);
                PHRead(dataContainer, fieldReceive[targetCell]);
            }
        }
    }
    else
    {
        RDouble **fieldReceive = reinterpret_cast <RDouble **> (grid->GetDataPtr(fieldName));

        for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; -- iGhostLayer)
        {
            for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
            {
                int iFace = interfaceIndexContainerForReceive[iLocalFace];
                int targetCell;
                grid->GetTargetIndex(iFace, iGhostLayer + 1, targetCell);
                for (int m = 0; m < nEquation; ++ m)
                {
                    PHRead(dataContainer, fieldReceive[m][targetCell]);
                }
            }
        }
    }
}

void TurbSolverUnstr::Boundary(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    using namespace PHENGLEI;
    CompositeUnstructBC *unstructBCRegion = grid->GetUnstructBCRegion();
    int nBCRegionUnstruct = unstructBCRegion->GetnBCRegion();

    for (int iBCRegionUnstruct = 0; iBCRegionUnstruct < nBCRegionUnstruct; iBCRegionUnstruct++)
    {
        UnstructBC *bcRegion = unstructBCRegion->GetBCRegion(iBCRegionUnstruct);
        int bcType = bcRegion->GetBCType();

        if (bcType < 0)
        {
            continue;
        }
        if (bcType == EXTRAPOLATION || bcType == SYMMETRY || bcType == OUTFLOW || bcType == PRESSURE_OUTLET || bcType == MASS_FLOW_OUTLET)
        {
            OutflowBCRegion(grid, bcRegion);
        }
        else if(bcType == INFLOW || bcType == PRESSURE_INLET || bcType == MASS_FLOW_INLET)
        {
            InflowBCRegion(grid, bcRegion);
        }
        else if(IsWall(bcType))
        {
            WallBCRegion(grid, bcRegion);
        }
        else if (bcType == FARFIELD)
        {
            FarfieldBCRegion(grid, bcRegion);
        }else
        {
            OutflowBCRegion(grid, bcRegion);
        }
    }
}

void TurbSolverUnstr::OutflowBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *leftCellOfFace  = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();
    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
    RDouble *viscousTurbulence = reinterpret_cast<RDouble  *> (grid->GetDataPtr("vist"));
    
    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();
    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        // iFace is the face number in the set of faceIndex.
        int iFace = *iter;
        int le = leftCellOfFace[iFace];
        int re = rightCellOfFace[iFace];

        viscousTurbulence[re] = viscousTurbulence[le];
        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            qTurbulence[m][re] = qTurbulence[m][le];
        }
    }
}

void TurbSolverUnstr::InflowBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *leftCellOfFace  = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();
    RDouble *freeStreamTurbVar =  parameters->GetFreeStreamTurbVar();
    RDouble freeStreamViscosity = parameters->GetFreeStreamViscosity();
    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
    RDouble *viscousTurbulence = reinterpret_cast<RDouble  *> (grid->GetDataPtr("vist"));

    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();    
    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        // iFace is the face number in the set of faceIndex.
        int iFace = *iter;
        int le = leftCellOfFace[iFace];
        int re = rightCellOfFace[iFace];
    
        viscousTurbulence[re] = freeStreamViscosity;
        for (int m = 0; m < nTurbulenceEquation; ++ m)
        {
            qTurbulence[m][re] = freeStreamTurbVar[m];
        }
    }
}

void TurbSolverUnstr::WallBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *leftCellOfFace  = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    RDouble *wallDistance = grid->GetWallDist();

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int viscoueType = parameters->GetViscousType();
    RDouble refReNumber = parameters->GetRefReNumber();
    int transitionType = parameters->GetTransitionType();
    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
    RDouble **q      = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble *viscousLaminar = reinterpret_cast<RDouble  *> (grid->GetDataPtr("visl"));
    RDouble *viscousTurbulence = reinterpret_cast<RDouble  *> (grid->GetDataPtr("vist"));

    using namespace IDX;
    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();
    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        // iFace is the face number in the set of faceIndex.
        int iFace = *iter;
        int le = leftCellOfFace[iFace];
        int re = rightCellOfFace[iFace];

        viscousTurbulence[re] = -viscousTurbulence[le];
        if (viscoueType == 3)
        {
            qTurbulence[ISA][re] = - qTurbulence[ISA][le];
        }
        else if (viscoueType == 4)
        {
            RDouble beta1 = 0.075;
            RDouble rho   = q[IR][le];
            RDouble dist  = wallDistance[le];
            RDouble dist2 = dist * dist;

            RDouble kwWall = 60.0 * viscousLaminar[le] / (rho * beta1 * dist2 * refReNumber * refReNumber);

            qTurbulence[IKE][re] = -qTurbulence[IKE][le];
            //qTurbulence[IKW][re] = MAX(two * kw_wall - qTurbulence[IKW][le],SMALL);
            qTurbulence[IKW][re] = two * kwWall - qTurbulence[IKW][le];

            if (transitionType == IREGAMA)
            {
                qTurbulence[IGAMA][re] = qTurbulence[IGAMA][le];
                qTurbulence[IRECT][re] = qTurbulence[IRECT][le];
            }
        }
    }


}

void TurbSolverUnstr::FarfieldBCRegion(Grid *gridIn, UnstructBC *bcRegionUnstruct)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int *leftCellOfFace  = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    RDouble  *nxs  = grid->GetFaceNormalX();
    RDouble  *nys  = grid->GetFaceNormalY();
    RDouble  *nzs  = grid->GetFaceNormalZ();

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();
    RDouble refGama = parameters->GetRefGama();
    RDouble *freeStreamTurbVar =  parameters->GetFreeStreamTurbVar();
    RDouble freeStreamViscosity = parameters->GetFreeStreamViscosity();

    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
    RDouble **q      = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble *viscousTurbulence = reinterpret_cast<RDouble  *> (grid->GetDataPtr("vist"));
    RDouble *gamma = reinterpret_cast<RDouble *> (grid->GetDataPtr("gama"));

    using namespace IDX;
    RDouble *primitiveVariableInflow = reinterpret_cast<RDouble *> (GlobalDataBase::GetDataPtr("prim_inf"));
    RDouble roo = primitiveVariableInflow[IR];
    RDouble uoo = primitiveVariableInflow[IU];
    RDouble voo = primitiveVariableInflow[IV];
    RDouble woo = primitiveVariableInflow[IW];
    RDouble poo = primitiveVariableInflow[IP];
    
    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();
    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        // iFace is the face number in the set of faceIndex.
        int iFace = *iter;
        int le = leftCellOfFace[iFace];
        int re = rightCellOfFace[iFace];

        RDouble rin,uin,vin,win,pin;
        RDouble vno,vni,vei,coo,cin;
        RDouble riemp,riemm,cb,vnb;

        rin = q[IR][le];
        uin = q[IU][le];
        vin = q[IV][le];
        win = q[IW][le];
        pin = q[IP][le];

        vno = nxs[iFace] * uoo + nys[iFace] * voo + nzs[iFace] * woo;
        vni = nxs[iFace] * uin + nys[iFace] * vin + nzs[iFace] * win;
        vei = sqrt(uin * uin + vin * vin + win * win);

        RDouble gama = gamma[le];

        coo = sqrt(ABS(refGama * poo / roo));
        cin = sqrt(ABS(gama  * pin / rin));

        //! supersonic
        if (vei > cin)
        {
            if (vni >= 0.0)
            {
                //! exit
                viscousTurbulence[re] = viscousTurbulence[le];
                for (int m = 0; m < nTurbulenceEquation; ++ m)
                {
                    qTurbulence[m][re] = qTurbulence[m][le];
                }
            }
            else
            {
                //! inlet
                viscousTurbulence[re] = freeStreamViscosity;
                for (int m = 0; m < nTurbulenceEquation; ++ m)
                {
                    qTurbulence[m][re] = freeStreamTurbVar[m];
                }
            }
            continue;
        }

        //! subsonic
        riemp = vni + 2.0 * cin / (gama  - 1.0);
        riemm = vno - 2.0 * coo / (refGama - 1.0);
        vnb   = half   * (riemp + riemm);
        cb    = fourth * (riemp - riemm) * gama;

        if (vnb >= 0.0)
        {
            //! exit
            viscousTurbulence[re] = viscousTurbulence[le];
            for (int m = 0; m < nTurbulenceEquation; ++ m)
            {
                qTurbulence[m][re] = qTurbulence[m][le];
            }
        }
        else
        {
            //! inlet
            viscousTurbulence[re] = freeStreamViscosity;
            for (int m = 0; m < nTurbulenceEquation; ++ m)
            {
                qTurbulence[m][re] = freeStreamTurbVar[m];
            }
        }

    }
}

void TurbSolverUnstr::ComputeViscousCoeff(Grid *grid)
{
    Viscosity(grid);
}

void TurbSolverUnstr::Viscosity(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    RDouble **q = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));

    RDouble *viscousLaminar = reinterpret_cast<RDouble * > (grid->GetDataPtr("visl"));
    RDouble *viscousTurbulence = reinterpret_cast<RDouble * > (grid->GetDataPtr("vist"));

    RDouble *SST_F2 = reinterpret_cast<RDouble *> (grid->GetDataPtr("SST_F2"));

    RDouble *wallDistance = grid->GetWallDist();

    int numberOfTotalCell = grid->GetNTotalCell();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal = numberOfBoundaryFace + numberOfTotalCell;

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    RDouble eddyViscosityLimit = parameters->GetEddyViscosityLimit();
    int monitorViscousTurbulenceMax = GlobalDataBase::GetIntParaFromDB("monitor_vistmax");

    string viscousName  = parameters->GetViscousName();
    RDouble refReNumber = parameters->GetRefReNumber();
    RDouble reynoldsSquare = refReNumber * refReNumber;

    RDouble viscousTurbulenceMax = 0.0;
    RDouble viscousTurbulenceMin = 1.0e30;

    RDouble part1, part2, arg2, f2;
    RDouble ld, ld3, fv1;
    RDouble rho, ke, kw, ds, ds2;
    //RDouble s11, s22, s33, s12, s13, s23, sii, divv;

    int iMax = 0;
    int iMin = 0;

    using namespace IDX;

    if (viscousName.substr(0, 6) == "1eq-sa")
    {
        RDouble nuoo  = 1.0;        //!Nondimensional inflow laminar viscosity.
        RDouble nuMin = 1.0e-5 * nuoo;
        RDouble nuMax = 1.0e10 * nuoo;

        RDouble SA_cv1_cube = parameters->GetSA_cv1_cube();

        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            ld  = qTurbulence[0][iCell] * q[IR][iCell] / (viscousLaminar[iCell] + SMALL);
            ld3  = ld * ld * ld;
            fv1 = ld3 / (ld3 + SA_cv1_cube);
            viscousTurbulence[iCell] = q[IR][iCell] * fv1 * qTurbulence[0][iCell];
            //viscousTurbulence[iCell] = MIN(eddyViscosityLimit * viscousLaminar[iCell], viscousTurbulence[iCell]);  //viscousLaminar[iCell] * coef_kvist *
            
            viscousTurbulence[iCell] = MIN(nuMax, viscousTurbulence[iCell]);
            viscousTurbulence[iCell] = MAX(nuMin, viscousTurbulence[iCell]);

            if (iCell < numberOfTotalCell && viscousTurbulenceMax < viscousTurbulence[iCell])    // Bell 20130527 add iCell <numberOfTotalCell
            {
                viscousTurbulenceMax = viscousTurbulence[iCell];
                iMax = iCell;
            }
            if (iCell < numberOfTotalCell && viscousTurbulenceMin > viscousTurbulence[iCell])    // Bell 20130527 add iCell <numberOfTotalCell
            {
                viscousTurbulenceMin = viscousTurbulence[iCell];
                iMin = iCell;
            }
        }
    }
    else if (viscousName.substr(0, 6) == "2eq-kw")
    {
        if (viscousName.substr(0, 17) == "2eq-kw-menter-sst")
        {
            //! Original Menter k-w-SST model
            //! Important: no repalce! otherwise it can not repeat.
            RDouble SST_a1 = parameters->GetSST_a1();
            RDouble SST_betaStar = parameters->GetSST_betaStar();

            RDouble *dudx = new RDouble [numberOfTotal];
            RDouble *dudy = new RDouble [numberOfTotal];
            RDouble *dudz = new RDouble [numberOfTotal];

            RDouble *dvdx = new RDouble [numberOfTotal];
            RDouble *dvdy = new RDouble [numberOfTotal];
            RDouble *dvdz = new RDouble [numberOfTotal];

            RDouble *dwdx = new RDouble [numberOfTotal];
            RDouble *dwdy = new RDouble [numberOfTotal];
            RDouble *dwdz = new RDouble [numberOfTotal];

            grid->CompGradient(q[IU], dudx, dudy, dudz);
            grid->CompGradient(q[IV], dvdx, dvdy, dvdz);
            grid->CompGradient(q[IW], dwdx, dwdy, dwdz);

            RDouble *SpSdRatio = reinterpret_cast<RDouble *> (grid->GetDataPtr("SpSdRatio"));

            for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
            {
                rho = ABS(q[IR][iCell]) + SMALL;
                ke  = qTurbulence[IKE][iCell];
                kw  = qTurbulence[IKW][iCell];
                ds  = wallDistance[iCell];
                ds2 = ds * ds;

                RDouble s11 = dudx[iCell];
                RDouble s22 = dvdy[iCell];
                RDouble s33 = dwdz[iCell];
                RDouble s12 = half * (dudy[iCell] + dvdx[iCell]);
                RDouble s13 = half * (dudz[iCell] + dwdx[iCell]);
                RDouble s23 = half * (dvdz[iCell] + dwdy[iCell]);

                RDouble sij2 = two * (PHSPACE::SQR(s11, s22, s33) + two * PHSPACE::SQR(s12, s13, s23));

                part1 = 2.0 * sqrt(ke) / (SST_betaStar * kw * ds * refReNumber);
                part2 = 500.0 * viscousLaminar[iCell] / (rho * kw * ds2 * reynoldsSquare);
                arg2  = MAX(part1, part2);
                f2    = tanh(arg2 * arg2);
                SST_F2[iCell] = f2;

                //viscousTurbulence[iCell] = rho * ke / MAX(kw, vort * f2 / turbA1) * refReNumber;
                viscousTurbulence[iCell] = rho * ke / MAX(kw, sqrt(sij2) * f2 / (SST_a1 * refReNumber));

                int transitionType = parameters->GetTransitionType();
                if (transitionType == IREGAMA)
                {
                    viscousTurbulence[iCell] = viscousTurbulence[iCell] * sqrt(SpSdRatio[iCell]);
                }

                viscousTurbulence[iCell] = MIN(eddyViscosityLimit, viscousTurbulence[iCell]); //viscousLaminar[iCell] * coef_kvist * 

                if (viscousTurbulenceMax < viscousTurbulence[iCell]) 
                {
                    viscousTurbulenceMax = viscousTurbulence[iCell];
                    iMax = iCell;
                }
                if (viscousTurbulenceMin > viscousTurbulence[iCell])
                {
                    viscousTurbulenceMin = viscousTurbulence[iCell];
                    iMin = iCell;
                }
            }

            delete []dudx;
            delete []dudy;
            delete []dudz;

            delete []dvdx;
            delete []dvdy;
            delete []dvdz;

            delete []dwdx;
            delete []dwdy;
            delete []dwdz;

        }
    }

    RDouble vistMax = viscousTurbulenceMax;
    RDouble vistMin = viscousTurbulenceMin;

    grid->UpdateData("vist_max", &vistMax, PHDOUBLE, 1);
    grid->UpdateData("vist_min", &vistMin, PHDOUBLE, 1);

    if (monitorViscousTurbulenceMax)
    {
        grid->CompareMaxMinValue(viscousTurbulenceMax, 1);
        grid->CompareMaxMinValue(viscousTurbulenceMin, 2);

        ostringstream oss;
        cout << setprecision(6);
        cout << setiosflags(ios::scientific);
        oss  << "Max && min viscosity: " << viscousTurbulenceMax << " " << viscousTurbulenceMin << " " << "\n";
        PrintToWindow(oss);
    }
}

void TurbSolverUnstr::InitCGrid(Grid *fineGridIn, Grid *coarseGridIn)
{
    UnstructGrid *fineGrid = UnstructGridCast(fineGridIn);
    UnstructGrid *coarseGrid = UnstructGridCast(coarseGridIn);

    int *cell2CoarseGridCell = fineGrid->GetCell2CoarseGridCell();
    int numberOfTotalCellOnFineGrid = fineGrid->GetNTotalCell();
    int numberOfTotalCellOnCoarseGrid = coarseGrid->GetNTotalCell();

    RDouble *volumeOfCoarseGrid = coarseGrid->GetCellVolume();
    RDouble *volumeOfFineGrid = fineGrid->GetCellVolume();

    RDouble *walldistanceOfFineGrid = fineGrid->GetWallDist();

    if (!coarseGrid->GetWallDist())
    {
        coarseGrid->AllocateWalldist();
    }
    RDouble *walldistanceOfCoarseGrid = coarseGrid->GetWallDist();

    for (int iCell = 0; iCell <numberOfTotalCellOnCoarseGrid; ++ iCell)
    {
        walldistanceOfCoarseGrid[iCell] = 0.0;
    }

    for (int iCell = 0; iCell <numberOfTotalCellOnFineGrid; ++ iCell)
    {
        walldistanceOfCoarseGrid[cell2CoarseGridCell[iCell]] += volumeOfFineGrid[iCell] * walldistanceOfFineGrid[iCell];
    }

    for (int iCell = 0; iCell <numberOfTotalCellOnCoarseGrid; ++ iCell)
    {
        walldistanceOfCoarseGrid[iCell] /= volumeOfCoarseGrid[iCell];
    }
}

void TurbSolverUnstr::UploadInterfaceData(ActionKey *actkey)
{
    UploadInterfaceValue(actkey);
}

void TurbSolverUnstr::UploadInterfaceValue(ActionKey *actkey)
{
    int level = actkey->level;
    UnstructGrid *grid = reinterpret_cast<UnstructGrid *>(GetGrid(level));

    InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
    if (interfaceInformation == 0) return;

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();
    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
    PHSPACE::UploadInterfaceValue(grid, qTurbulence, "turb::q", nTurbulenceEquation);

    RDouble* viscousTurbulence = reinterpret_cast<RDouble*> (grid->GetDataPtr("vist"));
    PHSPACE::UploadInterfaceValue(grid, viscousTurbulence, "vist");
}

void TurbSolverUnstr::UploadInterpointData(ActionKey *actkey)
{
    int level = actkey->level;
    UnstructGrid *grid = reinterpret_cast<UnstructGrid *>(GetGrid(level));
    InterpointInformation *interPointInformation = grid->GetInterpointInfo();
    if (interPointInformation == 0)
    {
        return;
    }

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();
    RDouble **qTurbulenceNode = reinterpret_cast<RDouble **> (grid->GetDataPtr("qTurbNode"));
    PHSPACE::UploadInterpointValue(grid, qTurbulenceNode, "qTurbNode", nTurbulenceEquation);
}

void TurbSolverUnstr::DownloadInterpointData(ActionKey *actkey)
{
    int level = actkey->level;
    UnstructGrid *grid = reinterpret_cast<UnstructGrid *>(GetGrid(level));
    InterpointInformation *interPointInformation = grid->GetInterpointInfo();
    if (interPointInformation == 0)
    {
        return;
    }

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();
    RDouble **qTurbulenceNode = reinterpret_cast<RDouble **> (grid->GetDataPtr("qTurbNode"));
    PHSPACE::DownloadInterpointValue(grid, qTurbulenceNode, "qTurbNode", nTurbulenceEquation);
}

void TurbSolverUnstr::CommunicationInterpointWeight(ActionKey *actkey)
{
    int level = actkey->level;

    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();
    for (int izone = 0; izone < nZones; ++ izone)
    {
        int myid = GetCurrentProcessorID();
        int tag = izone + level * nZones;
        int send_proc = GetZoneProcessorID(izone);

        ZoneNeighbor *zoneNeighborForPoint = zoneConnectivityForPoint->GetZoneNeighborForPoint(izone);
        std::size_t numberOfNeighbor = zoneNeighborForPoint->GetNumberOfNeighbor();

        for (std::size_t iNeighborZone = 0; iNeighborZone < numberOfNeighbor; ++ iNeighborZone)
        {
            int neighborZoneIndex = zoneNeighborForPoint->GetZoneIndexOfNeighbor(iNeighborZone);
            int recv_proc = GetZoneProcessorID(neighborZoneIndex);

            DataContainer * cdata = new DataContainer();

            if (myid == send_proc)
            {
                UnstructGrid *grid = UnstructGridCast(PHSPACE::GetGrid(izone, level));
                InterpointInformation *interpointInformation = grid->GetInterpointInfo();
                int *nodeValueSliceTurb = reinterpret_cast<int *>(grid->GetDataPtr("nodeValueSliceTurb"));

                cdata->moveto_begin();

                int iNeighbor = interpointInformation->FindIthNeighbor(neighborZoneIndex);
                int numberOfInterpointsForNeighbor = interpointInformation->GetNumberOfInterpointsForNeighbor(iNeighbor);
                int *pointIndexForSend = interpointInformation->GetPointIndexForSend(iNeighbor);
                int *interPoint2GlobalPoint = interpointInformation->GetInterPoint2GlobalPoint();

                int iPoint;
                int globalPoint;
                int valueWeight;
                for (int ipointLocal = 0; ipointLocal < numberOfInterpointsForNeighbor; ++ ipointLocal)
                {
                    iPoint = pointIndexForSend[ipointLocal];
                    globalPoint = interPoint2GlobalPoint[iPoint];
                    valueWeight = nodeValueSliceTurb[globalPoint];

                    cdata->write(&valueWeight, sizeof(int));
                }
            }

            PH_Trade(cdata, send_proc, recv_proc, tag);

            if (myid == recv_proc)
            {
                UnstructGrid *gridNeighbor = UnstructGridCast(PHSPACE::GetGrid(neighborZoneIndex, level));
                InterpointInformation *interpointInformationNeighbor = gridNeighbor->GetInterpointInfo();
                int *nodeValueSliceTurbTrade = reinterpret_cast<int *>(gridNeighbor->GetDataPtr("nodeValueSliceTurbTrade"));

                cdata->moveto_begin();

                int iNeighbor = interpointInformationNeighbor->FindIthNeighbor(izone);
                int numberOfInterpointsForNeighbor = interpointInformationNeighbor->GetNumberOfInterpointsForNeighbor(iNeighbor);
                int *pointIndexForReceive = interpointInformationNeighbor->GetPointIndexForReceive(iNeighbor);
                int *interPoint2GlobalPoint = interpointInformationNeighbor->GetInterPoint2GlobalPoint();

                int iPoint;
                int globalPoint;
                int valueWeight;
                for (int j = 0; j < numberOfInterpointsForNeighbor; ++ j)
                {
                    iPoint = pointIndexForReceive[j];
                    globalPoint = interPoint2GlobalPoint[iPoint];
                    cdata->read(&valueWeight, sizeof(int));

                    nodeValueSliceTurbTrade[globalPoint] += valueWeight;
                }
            }

            delete cdata;
        }
    }
}

void TurbSolverUnstr::DownloadInterpointWeight(ActionKey *actkey)
{
    int level = actkey->level;

    using namespace PHMPI;

    int myid = GetCurrentProcessorID();
    int nZones = GetNumberofGlobalZones();
    for (int izone = 0; izone < nZones; ++ izone)
    {
        int currentZoneProc = GetZoneProcessorID(izone);

        if (myid != currentZoneProc)
        {
            continue;
        }

        UnstructGrid *grid = UnstructGridCast(PHSPACE::GetGrid(izone, level));
        int nTotalNode = grid->GetNTotalNode();
        int *nodeValueSliceTurb = reinterpret_cast<int *>(grid->GetDataPtr("nodeValueSliceTurb"));
        int *nodeValueSliceTurbTrade = reinterpret_cast<int *>(grid->GetDataPtr("nodeValueSliceTurbTrade"));

        for (int iNode = 0; iNode < nTotalNode; ++ iNode)
        {
            nodeValueSliceTurb[iNode] += nodeValueSliceTurbTrade[iNode];
        }

        PHSPACE::SetField(nodeValueSliceTurbTrade, 0, nTotalNode);
    }
}

// Bell 20120910 add
void TurbSolverUnstr::DownloadInterfaceData(ActionKey *actkey)
{
    DownloadInterfaceValue(actkey);
}
//

void TurbSolverUnstr::DownloadInterfaceValue(ActionKey *actkey)
{
    int level = actkey->level;
    UnstructGrid *grid = reinterpret_cast<UnstructGrid *>(GetGrid(level));
    InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
    if (interfaceInformation == 0) return;

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();
    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
    PHSPACE::DownloadInterfaceValue(grid, qTurbulence, "turb::q", nTurbulenceEquation); 

    RDouble *viscousTurbulence = reinterpret_cast<RDouble *> (grid->GetDataPtr("vist"));
    PHSPACE::DownloadInterfaceValue(grid, viscousTurbulence, "vist");
}

//! need further compllished!!!
void TurbSolverUnstr::UploadOversetData(ActionKey *actkey)
{
    int level = actkey->level;
    UnstructGrid *grid = reinterpret_cast<UnstructGrid *>(GetGrid(level));

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
    PHSPACE::UploadOversetValue(grid, qTurbulence, "turb::q", nTurbulenceEquation);
}

//! need further compllished!!!
void TurbSolverUnstr::DownloadOversetData(ActionKey *actkey)
{
    int level = actkey->level;
    UnstructGrid *grid = reinterpret_cast<UnstructGrid *>(GetGrid(level));

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));

    PHSPACE::DownloadOversetValue(grid, qTurbulence, "turb::q", nTurbulenceEquation);
}

//* Bell 20120910 add
void TurbSolverUnstr::DumpInterfaceInfoBell()
{
    int level = 0;
    UnstructGrid *grid = UnstructGridCast(GetGrid(level));
    
    InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
    if (!interfaceInformation) return;

    fstream file;
    std::ostringstream oss;
    oss << AbsoluteFilename("interface") << grid->GetZoneID() << ".dat";
    file.open(oss.str().c_str(), ios_base::out);
    file.setf(ios::scientific);
    file.precision(8);
    int seprate = 20;

    int numberOfInterFace = interfaceInformation->GetNIFace();
    if (numberOfInterFace == 0) return;

    int *interFace2ZoneID       = interfaceInformation->GetInterFace2ZoneID();
    int *interFace2InterFaceID  = interfaceInformation->GetInterFace2InterFaceID();
    int *interface2BoundaryFace = interfaceInformation->GetInterFace2BoundaryFace();

    int numberOfBoundaryFace = grid->GetNBoundFace();

    int *leftCellOfFace  = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    RDouble *viscousLaminar = reinterpret_cast<RDouble *> (grid->GetDataPtr("visl"));
    RDouble *viscousTurbulence = reinterpret_cast<RDouble *> (grid->GetDataPtr("vist"));

    RDouble **Q1;
    RDouble **Q2;

    RDouble **visLaminar = new RDouble *[1];
    RDouble **VisTurbulence = new RDouble *[1];
    visLaminar[0] = viscousLaminar;
    VisTurbulence[0] = viscousTurbulence;
    Q1 = visLaminar;
    Q2 = VisTurbulence;

    file << " numberOfInterFace = " << numberOfInterFace << " numberOfBoundaryFace = " << numberOfBoundaryFace << "numberOfTotalCell = " << grid->GetNTotalCell() << "\n";
    file << "i" << setw(seprate) << "nb_zone" << setw(seprate) << "le" << setw(seprate) << "re" << setw(seprate) << "tle" << setw(seprate)  << "tre" << setw(seprate)
         << "qTurbulence[L]" << setw(seprate) << "qTurbulence[R]"
         << endl;
    
    for (int iFace = 0; iFace < numberOfInterFace; ++ iFace)
    {
        int targetZone = interFace2ZoneID[iFace];
        int targetInterface = interFace2InterFaceID[iFace];
        UnstructGrid *targetGrid = UnstructGridCast(PHSPACE::GetGrid(targetZone,0));

        if(!targetGrid) return;

        InterfaceInfo *targetInterfaceInformation = targetGrid->GetInterfaceInfo();
       
        int *targetInterface2BoundaryFace = targetInterfaceInformation->GetInterFace2BoundaryFace();

        int *targetLeftCellOfFace = targetGrid->GetLeftCellOfFace();
        int *targetRightCellOfFace = targetGrid->GetRightCellOfFace();

        int le = leftCellOfFace [interface2BoundaryFace[iFace]];
        int re = rightCellOfFace[interface2BoundaryFace[iFace]];

        int tle = targetLeftCellOfFace [targetInterface2BoundaryFace[targetInterface]];
        int tre = targetRightCellOfFace[targetInterface2BoundaryFace[targetInterface]];

        file << iFace << setw(seprate) << targetZone << setw(seprate) << le << setw(seprate) << re << setw(seprate) << tle << setw(seprate) << tre << setw(seprate)
             << Q1[0][le]<< setw(seprate) << Q1[0][re] << setw(seprate)
             << Q2[0][le]<< setw(seprate) << Q2[0][re]
             << endl;

    }
    delete [] visLaminar[1];    visLaminar[1] = NULL;
    delete [] VisTurbulence[1];    VisTurbulence[1] = NULL;
    file.close();
    file.clear();
}
//*/

void TurbSolverUnstr::ComputeLengthScaleofOneEquationModel(Grid *gridIn)
{
    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int DESType = parameters->GetDESType();

    if (DESType == DES)
    {
        ComputeDESLength(gridIn);
    }
    else if (DESType == DDES)
    {
        ComputeDDESLength(gridIn);
    }
    else if (DESType == IDDES)
    {
        ComputeIDDESLength(gridIn);
    }
}

void TurbSolverUnstr::ComputeDESLength(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell = grid->GetNTotalCell();
    RDouble *wallDistance = grid->GetWallDist();

    RDouble *DESLength = this->GetDESLength();
    if(!DESLength)
    {
        DESLength = new RDouble [numberOfTotalCell];
    }
    
    RDouble *lengthOfLES = new RDouble [numberOfTotalCell];
    RDouble *lowReynoldsNumberCorrection = new RDouble [numberOfTotalCell];

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int DESType = parameters->GetDESType();
    ComputeLESLengthofSA(gridIn, lengthOfLES, lowReynoldsNumberCorrection, DESType);

    RDouble lengthRANS, lengthLES;
    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell) 
    {
        lengthRANS = wallDistance[iCell];
        lengthLES  = lengthOfLES[iCell];

        DESLength[iCell] = MIN(lengthRANS, lengthLES);
    }

    this->SetDESLength(DESLength);
    
    delete []lengthOfLES;
    delete []lowReynoldsNumberCorrection;
}

void TurbSolverUnstr::ComputeDDESLength(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell = grid->GetNTotalCell();
    RDouble *wallDist = grid->GetWallDist();

    RDouble *DDESLength = this->GetDESLength();
    if(!DDESLength)
    {
        DDESLength = new RDouble [numberOfTotalCell];
    }

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int DESType = parameters->GetDESType();

    RDouble *lengthOfLES = new RDouble [numberOfTotalCell];
    RDouble *lowReynoldsNumberCorrection = new RDouble [numberOfTotalCell];
    ComputeLESLengthofSA(gridIn, lengthOfLES, lowReynoldsNumberCorrection, DESType);

    RDouble **q     = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble *viscousLaminar  = reinterpret_cast<RDouble * > (grid->GetDataPtr("visl"));
    RDouble *viscousTurbulence = reinterpret_cast<RDouble * > (grid->GetDataPtr("vist"));

    RDouble **gradientPrimtiveVariableX = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarX"));
    RDouble **gradientPrimtiveVariableY = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarY"));
    RDouble **gradientPrimtiveVariableZ = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarZ"));

    RDouble *dudx = gradientPrimtiveVariableX[1];
    RDouble *dudy = gradientPrimtiveVariableY[1];
    RDouble *dudz = gradientPrimtiveVariableZ[1];

    RDouble *dvdx = gradientPrimtiveVariableX[2];
    RDouble *dvdy = gradientPrimtiveVariableY[2];
    RDouble *dvdz = gradientPrimtiveVariableZ[2];

    RDouble *dwdx = gradientPrimtiveVariableX[3];
    RDouble *dwdy = gradientPrimtiveVariableY[3];
    RDouble *dwdz = gradientPrimtiveVariableZ[3];

    RDouble oRefReNumber = parameters->GetoRefReNumber();
    const RDouble karm   = 0.41;
    const RDouble kap2   = karm * karm;

    using namespace IDX;
    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell) 
    {
        RDouble laminarViscosity = viscousLaminar[iCell];
        RDouble turbulentViscosity = viscousTurbulence[iCell];
        RDouble density = q[IR][iCell];

        RDouble wallDistance  = wallDist[iCell];
        RDouble wallDistance2 = SQR(wallDistance);
        RDouble lengthRans    = wallDistance;
        RDouble lengthLes     = lengthOfLES[iCell];

        RDouble sumOfGradient2 = SQR(dudx[iCell], dudy[iCell], dudz[iCell]) + 
            SQR(dvdx[iCell], dvdy[iCell], dvdz[iCell]) + 
            SQR(dwdx[iCell], dwdy[iCell], dwdz[iCell]);

        RDouble nueff = oRefReNumber * (laminarViscosity + turbulentViscosity) / density;
        RDouble rd    = nueff / (kap2 * wallDistance2 * MAX(sqrt(sumOfGradient2), 1.0e-20));
        RDouble fd    = 1.0 - tanh(POWER3(8.0 * rd));

        DDESLength[iCell] = lengthRans - fd * MAX(0.0, (lengthRans - lengthLes));
    }

    this->SetDESLength(DDESLength);

    delete [] lengthOfLES;
    delete [] lowReynoldsNumberCorrection;
}

void TurbSolverUnstr::ComputeIDDESLength(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell = grid->GetNTotalCell();
    RDouble *wallDist = grid->GetWallDist();
    RDouble *largestLocalGridLength = grid->GetLargestLocalGridLength();

    RDouble *IDDESLength = this->GetDESLength();
    if(!IDDESLength)
    {
        IDDESLength = new RDouble [numberOfTotalCell];
    }

    RDouble *lengthOfLES = new RDouble [numberOfTotalCell];
    RDouble *lowReynoldsNumberCorrection = new RDouble [numberOfTotalCell];

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int DESType = parameters->GetDESType();
    ComputeLESLengthofSA(gridIn, lengthOfLES, lowReynoldsNumberCorrection, DESType);

    RDouble **q     = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble *viscousLaminar  = reinterpret_cast<RDouble * > (grid->GetDataPtr("visl"));
    RDouble *viscousTurbulence = reinterpret_cast<RDouble * > (grid->GetDataPtr("vist"));

    RDouble **gradientPrimtiveVariableX = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarX"));
    RDouble **gradientPrimtiveVariableY = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarY"));
    RDouble **gradientPrimtiveVariableZ = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarZ"));

    RDouble *dudx = gradientPrimtiveVariableX[1];
    RDouble *dudy = gradientPrimtiveVariableY[1];
    RDouble *dudz = gradientPrimtiveVariableZ[1];

    RDouble *dvdx = gradientPrimtiveVariableX[2];
    RDouble *dvdy = gradientPrimtiveVariableY[2];
    RDouble *dvdz = gradientPrimtiveVariableZ[2];

    RDouble *dwdx = gradientPrimtiveVariableX[3];
    RDouble *dwdy = gradientPrimtiveVariableY[3];
    RDouble *dwdz = gradientPrimtiveVariableZ[3];

    RDouble oRefReNumber = parameters->GetoRefReNumber();
    const RDouble karm = 0.41;
    const RDouble kap2 = karm * karm;

    using namespace IDX;
    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell) 
    {
        RDouble laminarViscosity        = viscousLaminar[iCell];
        RDouble turbulentViscosity      = viscousTurbulence[iCell];
        RDouble density                 = q[IR][iCell];

        RDouble wallDistance            = wallDist[iCell];
        RDouble wallDistance2           = SQR(wallDistance);
        RDouble lengthRans              = wallDistance;
        RDouble lengthLes               = lengthOfLES[iCell];
        RDouble largestLocalGridSpacing = largestLocalGridLength[iCell];

        RDouble alf  = 0.25 - wallDistance / largestLocalGridSpacing;
        RDouble fb   = MIN(2.0 * exp(- 9.0 * alf * alf),  1.0);

        RDouble sumOfGradient2 = SQR(dudx[iCell], dudy[iCell], dudz[iCell]) + 
                                SQR(dvdx[iCell], dvdy[iCell], dvdz[iCell]) + 
                                SQR(dwdx[iCell], dwdy[iCell], dwdz[iCell]);

        RDouble nuet = oRefReNumber * turbulentViscosity / density;
        RDouble rdt  = nuet / (kap2 * wallDistance2 * MAX(sqrt(sumOfGradient2), 1.0e-10));

        RDouble fdt  = 1.0 - tanh(POWER3(8.0 * rdt));
        RDouble fdb  = MAX((1.0 - fdt), fb);

        RDouble ct   = 1.63;
        RDouble cL   = 3.55;
        RDouble nuel = oRefReNumber * laminarViscosity / density;
        RDouble rdl  = nuel / (kap2 * wallDistance2 * MAX(sqrt(sumOfGradient2), 1.0e-10));

        RDouble fl   = tanh(pow(cL * cL * rdl, 10.0));
        RDouble ft   = tanh(POWER3(ct * ct * rdt));

        RDouble fe2  = 1.0 - MAX(ft, fl);

        RDouble fe1 = 0.0;
        if (alf < 0.0)
        {
            fe1 = 2.0 * exp(- 9.0 * alf * alf);
        }
        else
        {
            fe1 = 2.0 * exp(- 11.09 * alf * alf);
        }

        RDouble fe = fe2 * MAX((fe1 - 1.0), 0.0) * lowReynoldsNumberCorrection[iCell];

        IDDESLength[iCell]  =  fdb * (1.0 + fe) * lengthRans  +  (1.0 - fdb) * lengthLes;
    }

    this->SetDESLength(IDDESLength);

    delete []lengthOfLES;
    delete []lowReynoldsNumberCorrection;
}

RDouble * TurbSolverUnstr::GetLengthScale(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int DESType = parameters->GetDESType();

    if (DESType != 0)
    {
        return this->GetDESLength();
    }
    else
    {
        return grid->GetWallDist();
    }
}

void TurbSolverUnstr::ReadStatisticalFlow(ActionKey *actkey)
{
    int level = actkey->level;
    UnstructGrid *grid = UnstructGridCast(GetGrid(level));

    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotalCell = grid->GetNTotalCell();
    int numberOfTotal = numberOfTotalCell + numberOfBoundaryFace;

    RDouble **qAverageTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("qAverageTurb"));

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulenceEquation = parameters->GetNTurbulenceEquation();

    DataContainer *cdata = actkey->GetData();
    cdata->moveto_begin();

    int numberOfStatisticalStep = 0;
    cdata->read(&numberOfStatisticalStep, sizeof(int));

    GlobalDataBase::UpdateData("nStatisticalStep", &numberOfStatisticalStep, PHINT, 1);

    for (int m = 0; m < nTurbulenceEquation; ++ m)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            cdata->read(&qAverageTurbulence[m][iCell], sizeof(RDouble));
        }
    }
}

void ComputeLESLengthofSA(Grid *gridIn, RDouble *lengthOfLES, RDouble *lowReynoldsNumberCorrection, int DESType)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell = grid->GetNTotalCell();

    //! Compute LES length.
    RDouble CoeffOfDES = 0.65;

    ComputeLowReynoldsNumberCorrection(gridIn, lowReynoldsNumberCorrection);

    if (DESType == IDDES)
    {
        RDouble *subgridLength = grid->GetSubgridLength();

        for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell) 
        {
            lengthOfLES[iCell] = CoeffOfDES * lowReynoldsNumberCorrection[iCell] * subgridLength[iCell];
        }
    }
    else
    {
        RDouble *largestLocalGridLength = grid->GetLargestLocalGridLength();

        for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell) 
        {
            lengthOfLES[iCell] = CoeffOfDES * lowReynoldsNumberCorrection[iCell] * largestLocalGridLength[iCell];
        }
    }
}

void ComputeLowReynoldsNumberCorrection(Grid *gridIn, RDouble *lowReynoldsNumberCorrection)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell = grid->GetNTotalCell();

    //! Compute low Reynolds number correction.
    RDouble **q      = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
    RDouble *viscousLaminar   = reinterpret_cast<RDouble * > (grid->GetDataPtr("visl"));

    const RDouble karm   = 0.41;
    const RDouble cb1    = 0.1355;
    const RDouble cb2    = 0.622;
    const RDouble cv1    = 7.1;
    const RDouble sigma  = 2.0/3.0;
    const RDouble cv13   =  cv1 * cv1 * cv1;
    const RDouble kap2   = karm * karm;
    const RDouble rkap2  = one / kap2;
    const RDouble cw1    = cb1 * rkap2 + (one + cb2) / sigma;
    const RDouble turbulentFwStar          = 0.424;

    using namespace IDX;
    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell) 
    {
        RDouble nuet = qTurbulence[ISA][iCell];
        RDouble olam = q[IR][iCell] / (viscousLaminar[iCell] + SMALL);
        RDouble xsi  = nuet * olam + SMALL;
        RDouble xsi3 = POWER3(xsi);

        RDouble turbulentFt2 = 0.0; 
        //int codeOfFt2 = 0; //PHSPACE::GetCodeOfFt2();
        //if (codeOfFt2)
        //{
        //    turbulentFt2 = ct3 * exp(- ct4 * xsi2);
        //}

        RDouble turbulentFv1 = xsi3 / (xsi3 + cv13);
        RDouble turbulentFv2 = one - xsi / (one + xsi * turbulentFv1);

        RDouble numerator   =  1.0  -  (cb1 / (cw1 * kap2 * turbulentFwStar)) * (turbulentFt2 + (1.0 - turbulentFt2) * turbulentFv2);
        RDouble denominator =  turbulentFv1 * MAX(1.0e-10, (1.0 - turbulentFt2));

        lowReynoldsNumberCorrection[iCell] = sqrt(MIN(1.0e2, numerator/denominator));
    }
}

void TurbSolverUnstr::ComputeGradient(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    if(grid->IsFinestGrid() && TurbEquationFormation() == FORMATION_ORGSA_SST)
    {
        GetGradientField(gridIn);
    }
}

void TurbSolverUnstr::GetGradientField(Grid *gridIn)
{
    UnstructGrid *gridUnstruct = UnstructGridCast(gridIn);
    int numberOfEquations = this->GetNumberOfEquations();

    RDouble **gradientTurbulenceX = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("gradTurbulenceX"));
    RDouble **gradientTurbulenceY = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("gradTurbulenceY"));
    RDouble **gradientTurbulenceZ = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("gradTurbulenceZ"));
    RDouble **qTurbulenceNode = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("qTurbNode"));
    RDouble **qTurbulence = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("q_turb"));

    for (int m = 0; m < numberOfEquations; ++ m)
    {
        string gradientName = GlobalDataBase::GetStrParaFromDB("gradientName");
        if (gradientName == "ggnode" || gradientName == "ggnodelaplacian")
        {
            gridUnstruct->CompGradientGGNode_NEW(qTurbulence[m], qTurbulenceNode[m], gradientTurbulenceX[m], gradientTurbulenceY[m], gradientTurbulenceZ[m]);
        }
        else if (gradientName == "ggcell")
        {
            gridUnstruct->CompGradientGGCell(qTurbulence[m], gradientTurbulenceX[m], gradientTurbulenceY[m], gradientTurbulenceZ[m]);
        }
        else if (gradientName == "lsq")
        {
            gridUnstruct->CompGradientLSQ(qTurbulence[m], gradientTurbulenceX[m], gradientTurbulenceY[m], gradientTurbulenceZ[m]);
        }    
        else if (gradientName == "ggnode_weight")
        {
            gridUnstruct->CompGradientGGNodeWeight(qTurbulence[m], gradientTurbulenceX[m], gradientTurbulenceY[m], gradientTurbulenceZ[m]);
        }
        else if (gradientName == "gg_m2")
        {
            gridUnstruct->CompGradientGGModified2(qTurbulence[m], gradientTurbulenceX[m], gradientTurbulenceY[m], gradientTurbulenceZ[m]);
        }
        else if (gradientName == "ggcellnew")
        {
            gridUnstruct->CompGradientGGCellNew(qTurbulence[m], gradientTurbulenceX[m], gradientTurbulenceY[m], gradientTurbulenceZ[m]);
        }
        else if (gradientName == "ggcellw")
        {
            gridUnstruct->CompGradientGGCellW(qTurbulence[m], gradientTurbulenceX[m], gradientTurbulenceY[m], gradientTurbulenceZ[m]);
        }
        else
        {
            TK_Exit::ExceptionExit("No reconstruction method has been choosed ! /n");
        }
    }

    PHSPACE::CommunicateInterfaceValue(gridUnstruct, gradientTurbulenceX, "turb::dqdx", numberOfEquations);
    PHSPACE::CommunicateInterfaceValue(gridUnstruct, gradientTurbulenceY, "turb::dqdy", numberOfEquations);
    PHSPACE::CommunicateInterfaceValue(gridUnstruct, gradientTurbulenceZ, "turb::dqdz", numberOfEquations);
}

LIB_EXPORT void TurbSolverUnstr::InitControlParameters()
{
    FreeControlParameters();
    controlParameters = new Param_TurbSolverUnstruct();
    controlParameters->Init();

    //! Define turbulent equation formation: Original SA / CFL3D SA / NSMB SA, or SST.
    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int viscoueType = parameters->GetViscousType();
    if (viscoueType == 3)
    {
        //Turb_Equation_Formation = FORMATION_SA_CFL3D;
        Turb_Equation_Formation = FORMATION_SA_NSMB;
    }
    else
    {
        Turb_Equation_Formation = FORMATION_ORGSA_SST;
    }

    //! Compute machin zero. reference from cfl3d.
    machZeroExp = 1;
    RDouble add, x11;
    RDouble compare = 1.0;
    for(int iCount = 0; iCount < 20;++ iCount)
    {
        add = 1.0;
        for(int n = 0; n < iCount + 1; ++ n)
        {
            add = add * 0.1;
        }

        x11 = compare +add;
        if(x11 == compare)
        {
            machZeroExp = iCount + 1;
            break;
        }
    }
}

LIB_EXPORT Param_TurbSolverUnstruct * TurbSolverUnstr::GetControlParameters()
{
    return static_cast<Param_TurbSolverUnstruct *>(controlParameters);
}

}