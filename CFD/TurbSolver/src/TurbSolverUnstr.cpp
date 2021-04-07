#include "TurbSolverUnstr.h"
#include "BC_Unstruct.h"
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
    int nTurbulentEquation = parameters->GetNTurbulentEquation();
	if (nTurbulentEquation == 0) return;

    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell    = grid->GetNTotalCell();
    int numberOfTotalFace    = grid->GetNTotalFace();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal        = numberOfTotalCell + numberOfBoundaryFace;

    RDouble **qTurbulence        = NewPointer2<RDouble>(nTurbulentEquation, numberOfTotal);
    RDouble **residualTurbulence = NewPointer2<RDouble>(nTurbulentEquation, numberOfTotal);
    RDouble **spectrumTurbulence = NewPointer2<RDouble>(nTurbulentEquation, numberOfTotal);

    grid->UpdateDataPtr("q_turb"   , qTurbulence       );
    grid->UpdateDataPtr("res_turb" , residualTurbulence);
    grid->UpdateDataPtr("spec_turb", spectrumTurbulence);

    RDouble **gradientTurbulenceX = NewPointer2<RDouble>(nTurbulentEquation, numberOfTotal);
    RDouble **gradientTurbulenceY = NewPointer2<RDouble>(nTurbulentEquation, numberOfTotal);
    RDouble **gradientTurbulenceZ = NewPointer2<RDouble>(nTurbulentEquation, numberOfTotal);

    grid->UpdateDataPtr("gradTurbulentX", gradientTurbulenceX);
    grid->UpdateDataPtr("gradTurbulentY", gradientTurbulenceY);
    grid->UpdateDataPtr("gradTurbulentZ", gradientTurbulenceZ);

    int numberOfTotalNode = grid->GetNTotalNode();
    RDouble **qTurbulenceNode = NewPointer2<RDouble>(nTurbulentEquation, numberOfTotalNode);
    grid->UpdateDataPtr("qTurbNode", qTurbulenceNode);	
    
    int *nodeValueSliceTurb = new int[numberOfTotalNode];
    grid->UpdateDataPtr("nodeValueSliceTurb", nodeValueSliceTurb);

    int *nodeValueSliceTurbTrade = new int[numberOfTotalNode];
    PHSPACE::SetField(nodeValueSliceTurbTrade, 0, numberOfTotalNode);
    grid->UpdateDataPtr("nodeValueSliceTurbTrade", nodeValueSliceTurbTrade);

    InterpointInformation *interpointInformation = grid->GetInterpointInfo();
    if (interpointInformation)
    {
        int numberOfInterPoint = interpointInformation->GetNumberOfInterpoints();
        RDouble **qTurbulenceInterPoint = NewPointer2<RDouble>(nTurbulentEquation, numberOfInterPoint);
        grid->UpdateDataPtr("turb::qInterpoint", qTurbulenceInterPoint);
    }

    int isUnsteady = parameters->GetIsUnsteady();

	if (isUnsteady != 0)
	{
		RDouble **qTurbulenceUnsteadyN1 = NewPointer2<RDouble>(nTurbulentEquation, numberOfTotal);
		RDouble **qTurbulenceUnsteadyN2 = NewPointer2<RDouble>(nTurbulentEquation, numberOfTotal);
		RDouble **residualTurbulenceUnsteadyN1 = NewPointer2<RDouble>(nTurbulentEquation, numberOfTotal);
		RDouble **residualTurbulenceUnsteadyN2 = NewPointer2<RDouble>(nTurbulentEquation, numberOfTotal);

        //! Store the historical sum of InviscidFlux, ViscousFlux, et al., exclusive of DualTimeSourceFlux, used in UpdateUnsteadyFlow for the Rn in the dualtime method.
		RDouble **residualTurbulenceUnsteadyTemporary = NewPointer2<RDouble>(nTurbulentEquation, numberOfTotal);

		grid->UpdateDataPtr("q_turb_unsteady_n1", qTurbulenceUnsteadyN1);
		grid->UpdateDataPtr("q_turb_unsteady_n2", qTurbulenceUnsteadyN2);
        grid->UpdateDataPtr("res_turb_unsteady_n1", residualTurbulenceUnsteadyN1);
        grid->UpdateDataPtr("res_turb_unsteady_n2", residualTurbulenceUnsteadyN2);
        grid->UpdateDataPtr("res_turb_unsteady_tmp", residualTurbulenceUnsteadyTemporary);

		//! Statistical variables for DES simulation.
		int ifStaticsFlowField = parameters->GetIfStaticsFlowField();
		if (ifStaticsFlowField > 0)
		{
			RDouble **qAverageTurbulence = NewPointer2<RDouble>(nTurbulentEquation, numberOfTotal);
			grid->UpdateDataPtr("qAverageTurb", qAverageTurbulence);
		}
	}

    if (nTurbulentEquation == 2)
    {
        RDouble *cross = new RDouble[numberOfTotal];
        grid->UpdateDataPtr("cross", cross);
        RDouble *blend = new RDouble[numberOfTotal];
        grid->UpdateDataPtr("blend", blend);
        RDouble *SST_F2 = new RDouble[numberOfTotal];
        grid->UpdateDataPtr("SST_F2", SST_F2);
    }

	int viscousType = parameters->GetViscousType();
	int neasm    = GlobalDataBase::GetIntParaFromDB("neasm");

    if (viscousType >= 4 && neasm > 0)
    {
        int naniss = 5;
        RDouble **aniss = NewPointer2<RDouble>(naniss, numberOfTotal);

        grid->UpdateDataPtr("aniss", aniss);

        RDouble *sengy = new RDouble[numberOfTotal];

        grid->UpdateDataPtr("sengy", sengy);
    }

    RDouble **matrixTurbulenceLeft  = NewPointer2<RDouble>(nTurbulentEquation, numberOfTotalFace);
    RDouble **matrixTurbulenceRight = NewPointer2<RDouble>(nTurbulentEquation, numberOfTotalFace);

    grid->UpdateDataPtr("mat_turbl", matrixTurbulenceLeft);
	grid->UpdateDataPtr("mat_turbr", matrixTurbulenceRight);
}

void TurbSolverUnstr::DeAllocateGlobalVar(Grid *gridIn)
{
    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();
	if (nTurbulentEquation == 0) return;

	UnstructGrid *grid = UnstructGridCast(gridIn);

    RDouble **qTurbulence        = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"   ));
    RDouble **residualTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb" ));
	RDouble **spectrumTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("spec_turb"));

    DelPointer2(qTurbulence       );
    DelPointer2(residualTurbulence);
	DelPointer2(spectrumTurbulence);

    RDouble **gradientTurbulenceX = reinterpret_cast<RDouble **> (grid->GetDataPtr("gradTurbulentX"));
    RDouble **gradientTurbulenceY = reinterpret_cast<RDouble **> (grid->GetDataPtr("gradTurbulentY"));
    RDouble **gradientTurbulenceZ = reinterpret_cast<RDouble **> (grid->GetDataPtr("gradTurbulentZ"));

    DelPointer2(gradientTurbulenceX);
    DelPointer2(gradientTurbulenceY);
    DelPointer2(gradientTurbulenceZ);

	int isUnsteady = parameters->GetIsUnsteady();

	if (isUnsteady != 0)
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

    if (nTurbulentEquation == 2)
    {
        RDouble *cross  = reinterpret_cast<RDouble *> (grid->GetDataPtr("cross"));
        delete []cross;
        RDouble *blend  = reinterpret_cast<RDouble *> (grid->GetDataPtr("blend"));
        delete []blend;
		RDouble *SST_F2 = reinterpret_cast<RDouble *> (grid->GetDataPtr("SST_F2"));
        delete []SST_F2;
    }

	int viscousType = parameters->GetViscousType();
	int neasm    = GlobalDataBase::GetIntParaFromDB("neasm");

    if (viscousType >= 4 && neasm > 0)
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
	int nTurbulentEquation = parameters->GetNTurbulentEquation();
	int outerIterationStep = GlobalDataBase::GetIntParaFromDB("outnstep");

	RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));

    DataContainer *cdata = actkey->GetData();
    cdata->moveto_begin();
    cdata->write(&outerIterationStep, sizeof(int));

	for (int m = 0; m < nTurbulentEquation; ++ m)
	{
		cdata->write(qTurbulence[m], sizeof(RDouble) * numberOfTotal);
	}

	int isUnsteady = parameters->GetIsUnsteady();
	if (isUnsteady == 0) return;

	RDouble **qTurbulenceUnsteadyN1 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb_unsteady_n1"));
	RDouble **qTurbulenceUnsteadyN2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb_unsteady_n2"));

	for (int m = 0; m < nTurbulentEquation; ++ m)
	{
		cdata->write(qTurbulenceUnsteadyN1[m], sizeof(RDouble) * numberOfTotal);
	}

	for (int m = 0; m < nTurbulentEquation; ++ m)
	{
		cdata->write(qTurbulenceUnsteadyN2[m], sizeof(RDouble) * numberOfTotal);
	}

	RDouble **residualTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));
	RDouble **residualTurbulenceUnsteadyN1 = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb_unsteady_n1"));
	RDouble **residualTurbulenceUnsteadyN2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb_unsteady_n2"));

	for (int m = 0; m < nTurbulentEquation; ++ m)
	{
		cdata->write(residualTurbulence[m], sizeof(RDouble) * numberOfTotal);
	}

	for (int m = 0; m < nTurbulentEquation; ++ m)
	{
		cdata->write(residualTurbulenceUnsteadyN1[m], sizeof(RDouble) * numberOfTotal);
	}

	for (int m = 0; m < nTurbulentEquation; ++ m)
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
    int nTurbulentEquation = parameters->GetNTurbulentEquation();
    int outerIterationStep = GlobalDataBase::GetIntParaFromDB("outnstep");
    if (currentProcessorID == GetServerProcessorID())
    {
        WriteData(actkey->filepos, &version, "Version");
        WriteData(actkey->filepos, &outerIterationStep, "outnstep");

        if (IsNeedDESStatistics())
        {
            int nStatisticalStep = 0;
            GlobalDataBase::GetData("nStatisticalStep", &nStatisticalStep, PHINT, 1);
            WriteData(actkey->filepos, &nStatisticalStep, "nStatisticalStep");
        }
    }

    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
    RDouble **qTemporary  = NewPointer2<RDouble>(nTurbulentEquation, numberOfTotal);

    for (int iEquation = 0; iEquation < nTurbulentEquation; ++ iEquation)
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
    if (isUnsteady == 0)
    {
        H5Gclose(grploc);
        DelPointer2(qTemporary);

        return;
    }

    RDouble **qTurbulenceUnsteadyN1 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb_unsteady_n1"));
    for (int iEquation = 0; iEquation < nTurbulentEquation; ++ iEquation)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            qTemporary[iEquation][iCell] = qTurbulenceUnsteadyN1[iEquation][iCell];
        }
    }

    WriteData(grploc, qTemporary[0], "q_turb_unsteady_n1");

    RDouble **qTurbulenceUnsteadyN2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb_unsteady_n2"));
    for (int iEquation = 0; iEquation < nTurbulentEquation; ++ iEquation)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            qTemporary[iEquation][iCell] = qTurbulenceUnsteadyN2[iEquation][iCell];
        }
    }

    WriteData(grploc, qTemporary[0], "q_turb_unsteady_n2");

    RDouble **residualTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));
    for (int iEquation = 0; iEquation < nTurbulentEquation; ++ iEquation)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            qTemporary[iEquation][iCell] = residualTurbulence[iEquation][iCell];
        }
    }

    WriteData(grploc, qTemporary[0], "res_turb");

    RDouble **residualTurbulenceUnsteadyN1 = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb_unsteady_n1"));
    for (int iEquation = 0; iEquation < nTurbulentEquation; ++ iEquation)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            qTemporary[iEquation][iCell] = residualTurbulenceUnsteadyN1[iEquation][iCell];
        }
    }

    WriteData(grploc, qTemporary[0], "res_turb_unsteady_n1");

    RDouble **residualTurbulenceUnsteadyN2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb_unsteady_n2"));
    for (int iEquation = 0; iEquation < nTurbulentEquation; ++ iEquation)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            qTemporary[iEquation][iCell] = residualTurbulenceUnsteadyN2[iEquation][iCell];
        }
    }

    WriteData(grploc, qTemporary[0], "res_turb_unsteady_n2");

    if (IsNeedDESStatistics())
    {
        RDouble **qAverageTurb = reinterpret_cast<RDouble **> (grid->GetDataPtr("qAverageTurb"));
        for (int iEquation = 0; iEquation < nTurbulentEquation; ++ iEquation)
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
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
    RDouble **qTemporary  = NewPointer2<RDouble>(nTurbulentEquation, numberOfTotal);

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
    for (int iEquation = 0; iEquation < nTurbulentEquation; ++ iEquation)
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

            for (int m = 0; m < nTurbulentEquation; ++ m)
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
            for (int iEquation = 0; iEquation < nTurbulentEquation; ++ iEquation)
            {
                for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
                {
                    qTurbulenceUnsteadyN1[iEquation][iCell] = qTemporary[iEquation][iCell];
                }
            }

            ReadData(grploc, qTemporary[0], "q_turb_unsteady_n2");
            for (int iEquation = 0; iEquation < nTurbulentEquation; ++ iEquation)
            {
                for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
                {
                    qTurbulenceUnsteadyN2[iEquation][iCell] = qTemporary[iEquation][iCell];
                }
            }

            ReadData(grploc, qTemporary[0], "res_turb");
            for (int iEquation = 0; iEquation < nTurbulentEquation; ++ iEquation)
            {
                for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
                {
                    residualTurbulence[iEquation][iCell] = qTemporary[iEquation][iCell];
                }
            }

            ReadData(grploc, qTemporary[0], "res_turb_unsteady_n1");
            for (int iEquation = 0; iEquation < nTurbulentEquation; ++ iEquation)
            {
                for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
                {
                    residualTurbulenceUnsteadyN1[iEquation][iCell] = qTemporary[iEquation][iCell];
                }
            }

            ReadData(grploc, qTemporary[0], "res_turb_unsteady_n2");
            for (int iEquation = 0; iEquation < nTurbulentEquation; ++ iEquation)
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
            for (int iEquation = 0; iEquation < nTurbulentEquation; ++ iEquation)
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
            for (int iEquation = 0; iEquation < nTurbulentEquation; ++ iEquation)
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
	int nTurbulentEquation = parameters->GetNTurbulentEquation();

	RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));

    DataContainer *cdata = actkey->GetData();
    cdata->moveto_begin();
    int outerIterationStepOfTurbulence = 0;
    cdata->read(&outerIterationStepOfTurbulence, sizeof(int));

	int outerIterationStepOfNS = GlobalDataBase::GetIntParaFromDB("outnstep");

	for (int m = 0; m < nTurbulentEquation; ++ m)
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

			for (int m = 0; m < nTurbulentEquation; ++ m)
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
			for (int m = 0; m < nTurbulentEquation; ++ m)
			{
				cdata->read(qTurbulenceUnsteadyN1[m], sizeof(RDouble) * numberOfTotal);
			}
			for (int m = 0; m < nTurbulentEquation; ++ m)
			{
				cdata->read(qTurbulenceUnsteadyN2[m], sizeof(RDouble) * numberOfTotal);
			}

			for (int m = 0; m < nTurbulentEquation; ++ m)
			{
				cdata->read(residualTurbulence[m], sizeof(RDouble) * numberOfTotal);
			}
			for (int m = 0; m < nTurbulentEquation; ++ m)
			{
				cdata->read(residualTurbulenceUnsteadyN1[m], sizeof(RDouble) * numberOfTotal);
			}
			for (int m = 0; m < nTurbulentEquation; ++ m)
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
    int nTurbulentEquation = parameters->GetNTurbulentEquation();
	if (nTurbulentEquation == 0) return;

	UnstructGrid *grid = UnstructGridCast(GetGrid(0));

	//int outerIterationStep = 0;
	//GlobalDataBase::UpdateData("outnstep", &outerIterationStep, PHINT, 1);

	int numberOfTotalCell    = grid->GetNTotalCell();
	int numberOfBoundaryFace = grid->GetNBoundFace();
	int numberOfTotal        = numberOfTotalCell + numberOfBoundaryFace;

	RDouble *freeStreamTurbVar =  parameters->GetFreeStreamTurbVar();
	RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));

	for (int m = 0; m < nTurbulentEquation; ++ m)
	{
		for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
		{
			qTurbulence[m][iCell] = freeStreamTurbVar[m];
		}
	}

	int isUnsteady = parameters->GetIsUnsteady();
    if (isUnsteady == 0) return;

	RDouble **qTurbulenceUnsteadyN1 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb_unsteady_n1"));
	RDouble **qTurbulenceUnsteadyN2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb_unsteady_n2"));

	for (int m = 0; m < nTurbulentEquation; ++ m)
	{
		for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
		{
			qTurbulenceUnsteadyN1[m][iCell] = qTurbulence[m][iCell];
		}
	}

	for (int m = 0; m < nTurbulentEquation; ++ m)
	{
		for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
		{
			qTurbulenceUnsteadyN2[m][iCell] = qTurbulence[m][iCell];
		}
	}

	RDouble **residualTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));
	RDouble **residualTurbulenceUnsteadyN1 = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb_unsteady_n1" ));
	RDouble **residualTurbulenceUnsteadyN2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb_unsteady_n2" ));
	RDouble **residualTurbulenceUnsteadyTemporary = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb_unsteady_tmp"));

	for (int m = 0; m < nTurbulentEquation; ++ m)
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
	int nTurbulentEquation = parameters->GetNTurbulentEquation();

	RDouble **residualTurbulence = reinterpret_cast<RDouble **>(grid->GetDataPtr("res_turb"));

	for (int iEquation = 0; iEquation < nTurbulentEquation; ++ iEquation)
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
	int nTurbulentEquation = parameters->GetNTurbulentEquation();
	if (nTurbulentEquation == 0) return;

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
        const int methodOfDualTime = 3;
        ComputeDualTimeCoefficient(methodOfDualTime, dualTimeCoefficient);
        dualTimeSpectrumC1 = - dualTimeCoefficient[3];
        dualTimeSpectrumC2 =   dualTimeCoefficient[6];
    }

    for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
    {
        for (int m = 0; m < nTurbulentEquation; ++ m)
        {
            spectrumTurbulence[m][iCell] = 0.0;
        }
    }

    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
    {
        for (int m = 0; m < nTurbulentEquation; ++ m)
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
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    for (int m = 0; m < nTurbulentEquation; ++ m)
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
	int nTurbulentEquation = parameters->GetNTurbulentEquation();

	for (int m = 0; m < nTurbulentEquation; ++ m)
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
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal        = numberOfTotalCell + numberOfBoundaryFace;

	RDouble **q = qProxy->GetField_UNS();
    RDouble **qOld = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
    {
		for (int m = 0; m < nTurbulentEquation; ++ m)
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
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    for (int m = 0; m < nTurbulentEquation; ++ m)
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
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal        = numberOfTotalCell + numberOfBoundaryFace;

	RDouble **qqTurbulence = qProxy->GetField_UNS();

    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
    {
		for (int m = 0; m < nTurbulentEquation; ++ m)
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
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    int iSmoothTurbulence = GlobalDataBase::GetIntParaFromDB("ismooth_turb");

    int nNegtiveMax = 0;

    RDouble turbulenceRelax = GlobalDataBase::GetDoubleParaFromDB("turb_relax");

    RDouble eddyViscosityLimit = parameters->GetEddyViscosityLimit();

	RDouble *volume = grid->GetCellVolume();

    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

	using namespace IDX;

	RDouble *primitiveVariable = new RDouble [nTurbulentEquation];

    if (nTurbulentEquation == 1)
    {
        //! Spalart-allmaras model.
		RDouble *viscousLaminar = reinterpret_cast<RDouble *> (grid->GetDataPtr("visl"));

		int nNegative = 0, nPositive = 0;
        int iMax = 0;

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
    else if (nTurbulentEquation == 2)
    {
		RDouble *freeStreamTurbVar = parameters->GetFreeStreamTurbVar();

        RDouble kelim = 1.0e-5 * freeStreamTurbVar[IKE];
        RDouble kwlim = 1.0e-5 * freeStreamTurbVar[IKW];

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
        }
    }

	delete []primitiveVariable;
}

void TurbSolverUnstr::SMoothTurbulencePoint(Grid *gridIn, RDouble *primitiveVariable, int i)
{
	UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell    = grid->GetNTotalCell();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal        = numberOfTotalCell + numberOfBoundaryFace;

    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
	RDouble *volume = grid->GetCellVolume() ;

    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
	int nTurbulentEquation = parameters->GetNTurbulentEquation();

    int **cell2Face = grid->GetCell2Face();
    int *faceNumberOfEachCell = grid->GetFaceNumberOfEachCell();

	int le,re,ie,face;
	RDouble volumeSum = zero;
	for (int m = 0; m < nTurbulentEquation; ++ m)
	{
		primitiveVariable[m] = zero;
	}

	for (int j = 0; j < faceNumberOfEachCell[i]; ++ j)
    {
        face = cell2Face[i][j];
        le    = leftCellOfFace[face];
        re    = rightCellOfFace[face];
			
        if (le != i)
		{
			ie = le;
		}
		else
		{
			ie = re;
		}

		if (ie >= numberOfTotalCell) continue;

		for (int m = 0; m < nTurbulentEquation; ++ m)
		{
			primitiveVariable[m] += qTurbulence[m][ie];
		}
		volumeSum += 1.0;
	}

	for (int m = 0; m < nTurbulentEquation; ++ m)
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
	int nTurbulentEquation = parameters->GetNTurbulentEquation();

    RDouble **tmporary = NewPointer2<RDouble>(nTurbulentEquation, numberOfTotal);

    int **cell2Face = grid->GetCell2Face();
    int *faceNumberOfEachCell = grid->GetFaceNumberOfEachCell();

	for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
	{
		for (int m = 0; m < nTurbulentEquation; ++ m)
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

			for (int m = 0; m < nTurbulentEquation; ++ m)
			{
				qTurbulence[m][iCell] += tmporary[m][ie];
			}
			volumeSum += volume[ie];
		}

		for (int m = 0; m < nTurbulentEquation; ++ m)
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

	RDouble *dt = reinterpret_cast<RDouble *> (grid->GetDataPtr("dt" ));

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
    {
		for (int m = 0; m < nTurbulentEquation; ++ m)
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

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    RDouble *df    = new RDouble [nTurbulentEquation];
    RDouble *rhs0  = new RDouble [nTurbulentEquation];
    RDouble *dqOld = new RDouble [nTurbulentEquation];

    RDouble **spectrumTurbulence    = reinterpret_cast<RDouble **> (grid->GetDataPtr("spec_turb"));
    RDouble **matrixTurbulenceLeft  = reinterpret_cast<RDouble **> (grid->GetDataPtr("mat_turbl"));
    RDouble **matrixTurbulenceRight = reinterpret_cast<RDouble **> (grid->GetDataPtr("mat_turbr"));
    RDouble **residualTurbulence    = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));

    int **cell2Face = grid->GetCell2Face();
    int *faceNumberOfEachCell = grid->GetFaceNumberOfEachCell();

    RDouble **deltaFlux = LUplusDQ->GetField_UNS();
    RDouble *matrixUx = new RDouble [nTurbulentEquation];

    //! Now the Forward Sweep
    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
    {
        for (int m = 0; m < nTurbulentEquation; ++ m)
        {
            rhs0[m] = 0.0;
        }

        for (int m = 0; m < nTurbulentEquation; ++ m)
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

                for (int m = 0; m < nTurbulentEquation; ++ m)
                {
                    df[m] = matrixTurbulenceLeft[m][face] * dq[m][re];
                }
            }
            else
            {
                for (int m = 0; m < nTurbulentEquation; ++ m)
                {
                    df[m] = matrixTurbulenceRight[m][face] * dq[m][re];
                }
            }

            //! Add Flux together
            for (int m = 0; m < nTurbulentEquation; ++ m)
            {
                rhs0[m] += df[m];
            }
        }

        for (int m = 0; m < nTurbulentEquation; ++ m)
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

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    RDouble *df    = new RDouble [nTurbulentEquation];
    RDouble *rhs0  = new RDouble [nTurbulentEquation];
    RDouble *dqOld = new RDouble [nTurbulentEquation];

    RDouble **spectrumTurbulence    = reinterpret_cast<RDouble **> (grid->GetDataPtr("spec_turb"));
    RDouble **matrixTurbulenceLeft  = reinterpret_cast<RDouble **> (grid->GetDataPtr("mat_turbl"));
    RDouble **matrixTurbulenceRight = reinterpret_cast<RDouble **> (grid->GetDataPtr("mat_turbr"));
    RDouble **residualTurbulence    = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));

    int **cell2Face = grid->GetCell2Face();
    int *faceNumberOfEachCell = grid->GetFaceNumberOfEachCell();

    RDouble **deltaFlux = LUplusDQ->GetField_UNS();
    RDouble *matrixLx = new RDouble [nTurbulentEquation];

    // Backward Sweep
    for (int iCell = numberOfTotalCell - 1; iCell >= 0; -- iCell)
    {
        for (int m = 0; m < nTurbulentEquation; ++ m)
        {
            rhs0[m] = 0.0;
        }

        for (int m = 0; m < nTurbulentEquation; ++ m)
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

                for (int m = 0; m < nTurbulentEquation; ++ m)
                {
                    df[m] = matrixTurbulenceLeft[m][face] * dq[m][re];
                }
            }
            else
            {
                for (int m = 0; m < nTurbulentEquation; ++ m)
                {
                    df[m] = matrixTurbulenceRight[m][face] * dq[m][re];
                }
            }

            for (int m = 0; m < nTurbulentEquation; ++ m)
            {
                rhs0[m] += df[m];
            }
        }

        for (int m = 0; m < nTurbulentEquation; ++ m)
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
	int nTurbulentEquation = parameters->GetNTurbulentEquation();
	if (nTurbulentEquation != 1)
	{
		cout <<"SetGhostDQLUSGS function is turned off !\n";
		return;
	}

	UnstructGrid *grid = UnstructGridCast(gridIn);
	int numberOfTotalCell    = grid->GetNTotalCell();
	int numberOfBoundaryFace = grid->GetNBoundFace();
	int numberOfTotal        = numberOfTotalCell + numberOfBoundaryFace;

	int      *leftCellOfFace  = grid->GetLeftCellOfFace();
	int      *rightCellOfFace  = grid->GetRightCellOfFace();

	//! Get grid metrices
	RDouble *xfn  = grid->GetFaceNormalX();
	RDouble *yfn  = grid->GetFaceNormalY();
	RDouble *zfn  = grid->GetFaceNormalZ();
	RDouble *area = grid->GetFaceArea();
	RDouble *vgn  = grid->GetFaceNormalVelocity();

	RDouble *cellCenterX = grid->GetCellCenterX();
	RDouble *cellCenterY = grid->GetCellCenterY();
	RDouble *cellCenterZ = grid->GetCellCenterZ();

	BCRecord **bcRecord = grid->GetBCRecord();	

	for (int i = 0; i < numberOfBoundaryFace; ++ i)
	{	
		int le     = leftCellOfFace[i];
		int re     = rightCellOfFace[i];

		int bcType = bcRecord[i]->GetKey();

		if (bcType < 0)
		{
			continue;
		}

		if (IsWall(bcType))
		{
			for (int m = 0; m < nTurbulentEquation; ++ m)
			{
				//dq[m][re] = - dq[m][le]; 
				dq[m][re] = 0.0;
			}
		} 
		else if (bcType == PHENGLEI::FARFIELD || bcType == PHENGLEI::SYMMETRY 
			  || bcType == PHENGLEI::INFLOW   || bcType == PHENGLEI::OUTFLOW)
		{
			for (int m = 0; m < nTurbulentEquation; ++ m)
			{
				dq[m][re] = dq[m][le]; 
			}
		}
		else
		{
			for (int m = 0; m < nTurbulentEquation; ++ m)
			{
				dq[m][re] = 0.0; 
			}
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
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    RDouble **residualTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));

    for (int m = 0; m < nTurbulentEquation; ++ m)
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
    if (isUnsteady == 0) return zero;

    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell    = grid->GetNTotalCell();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal        = numberOfTotalCell + numberOfBoundaryFace;

	RDouble **qPrimitive  = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
	RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
	RDouble **qTurbulenceUnsteadyN1 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb_unsteady_n1"));
	RDouble **qTurbulenceUnsteadyN2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb_unsteady_n2"));
	RDouble **residualTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));

    int nTurbulentEquation = parameters->GetNTurbulentEquation();

	RDouble *primitive0 = new RDouble[nTurbulentEquation];
	RDouble *primitive1 = new RDouble[nTurbulentEquation];
	RDouble *primitive2 = new RDouble[nTurbulentEquation];

	RDouble *qConservative0 = new RDouble[nTurbulentEquation];
	RDouble *qConservative1 = new RDouble[nTurbulentEquation];
	RDouble *qConservative2 = new RDouble[nTurbulentEquation];

	RDouble sum1 = zero;
	RDouble sum2 = zero;

	using namespace IDX;

	for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
	{
        for (int m = 0; m < nTurbulentEquation; ++ m)
		{
			primitive0[m] = qTurbulence[m][iCell];
			primitive1[m] = qTurbulenceUnsteadyN1[m][iCell];
			primitive2[m] = qTurbulenceUnsteadyN2[m][iCell];
		}

		RDouble coefficient = 1.0;

		if (nTurbulentEquation >= 2)
		{
			coefficient = qPrimitive[IR][iCell];
		}

		for (int m = 0; m < nTurbulentEquation; ++ m)
		{
			qConservative0[m] = coefficient * primitive0[m];
			qConservative1[m] = coefficient * primitive1[m];
			qConservative2[m] = coefficient * primitive2[m];
		}

		//! nl is needed.
        for (int m = 0; m < nTurbulentEquation; ++ m)
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
	if (isUnsteady == 0) return;

    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell    = grid->GetNTotalCell();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal        = numberOfTotalCell + numberOfBoundaryFace;

	RDouble *volume   = grid->GetCellVolume();
	RDouble *volumeN1 = grid->GetCellVolume();
	RDouble *volumeN2 = grid->GetCellVolume();

	RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
	RDouble **qTurbulenceUnsteadyN1 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb_unsteady_n1"));
	RDouble **qTurbulenceUnsteadyN2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb_unsteady_n2"));

	RDouble **residualTurbulenceUnsteadyN1 = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb_unsteady_n1" ));
	RDouble **residualTurbulenceUnsteadyN2 = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb_unsteady_n2" ));
	RDouble **residualTurbulenceUnsteadyTemporary = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb_unsteady_tmp"));

	int nTurbulentEquation = parameters->GetNTurbulentEquation();

    for (int m = 0; m < nTurbulentEquation; ++ m)
    {
		for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            qTurbulenceUnsteadyN2[m][iCell] = qTurbulenceUnsteadyN1[m][iCell];
        }
    }

    for (int m = 0; m < nTurbulentEquation; ++ m)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            qTurbulenceUnsteadyN1[m][iCell] = qTurbulence[m][iCell];
        }
    }

    for (int m = 0; m < nTurbulentEquation; ++ m)
    {
        for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
        {
            residualTurbulenceUnsteadyN2[m][iCell] = residualTurbulenceUnsteadyN1[m][iCell];
        }
    }

    for (int m = 0; m < nTurbulentEquation; ++ m)
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

			RDouble c1 = (numberOfStatisticalStep - 1.0) / (numberOfStatisticalStep * 1.0);
			RDouble c2 = 1.0 / (numberOfStatisticalStep * 1.0);

			RDouble **qAverageTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("qAverageTurb"));
			for (int m = 0; m < nTurbulentEquation; ++ m)
			{
				for (int iCell = 0; iCell < numberOfTotal; ++ iCell)
				{
					qAverageTurbulence[m][iCell] = c1 * qAverageTurbulence[m][iCell] + c2 * qTurbulence[m][iCell];
				}
			}
		}
	}
}

void TurbSolverUnstr::DualTimeSource(Grid *gridIn)
{
	Param_TurbSolverUnstruct *parameters = GetControlParameters();
	int isUnsteady = parameters->GetIsUnsteady();
    if (isUnsteady == 0) return;

    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell    = grid->GetNTotalCell();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal        = numberOfTotalCell + numberOfBoundaryFace;

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

    int nTurbulentEquation = parameters->GetNTurbulentEquation();

	RDouble *primitive0 = new RDouble[nTurbulentEquation];
	RDouble *primitive1 = new RDouble[nTurbulentEquation];
	RDouble *primitive2 = new RDouble[nTurbulentEquation];

	RDouble *qConservative0 = new RDouble[nTurbulentEquation];
	RDouble *qConservative1 = new RDouble[nTurbulentEquation];
	RDouble *qConservative2 = new RDouble[nTurbulentEquation];

    //! Store the historical sum of InviscidFlux, ViscousFlux, et al., exclusive of DualTimeSourceFlux, used in UpdateUnsteadyFlow for the Rn in the dualtime method.
    for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
	{
	    for (int m = 0; m < nTurbulentEquation; ++ m)
		{
			residualTurbulenceUnsteadyTemporary[m][iCell] = residualTurbulence[m][iCell];
		}
	}

    //! Computation of dualtime coefficients, including three coefficients for Residual of R(p), R(n) and R(n-1);
    //! and three coefficients for conservative variables of qcsv(p), qcsv(n), and qcsv(n-1).
    RDouble dualTimeCoefficient[7];
    const int methodOfDualTime = 3;
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
        for (int m = 0; m < nTurbulentEquation; ++ m)
		{
			primitive0[m] = qTurbulence[m][iCell];
			primitive1[m] = qTurbulenceUnsteadyN1[m][iCell];
			primitive2[m] = qTurbulenceUnsteadyN2[m][iCell];
		}

		RDouble coefficient = 1.0;

		if (nTurbulentEquation >= 2)
		{
			coefficient  = qPrimitive[IR][iCell];
		}

		for (int m = 0; m < nTurbulentEquation; ++ m)
		{
			qConservative0[m] = coefficient * primitive0[m];
			qConservative1[m] = coefficient * primitive1[m];
			qConservative2[m] = coefficient * primitive2[m];
		}

	    for (int m = 0; m < nTurbulentEquation; ++ m)
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
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

	RDouble **field = fieldProxy->GetField_UNS();

    for (int m = 0; m < nTurbulentEquation; ++ m)
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
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

	RDouble **field1 = field1Proxy->GetField_UNS();
	RDouble **field2 = field2Proxy->GetField_UNS();

    for (int m = 0; m < nTurbulentEquation; ++ m)
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
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

	RDouble **field = NewPointer2<RDouble>(nTurbulentEquation, numberOfTotal);

	FieldProxy *fieldProxy = new FieldProxy();

	fieldProxy->SetField_UNS(field, true);

	return fieldProxy;
}

FieldProxy * TurbSolverUnstr::GetFieldProxy(Grid *gridIn, const string &fieldName)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell    = grid->GetNTotalCell();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal        = numberOfTotalCell + numberOfBoundaryFace;

	RDouble **field = reinterpret_cast<RDouble **> (grid->GetDataPtr(fieldName));

	FieldProxy *fieldProxy = new FieldProxy();

	fieldProxy->SetField_UNS(field);

	return fieldProxy;
}

FieldProxy * TurbSolverUnstr::GetResidualProxy(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell    = grid->GetNTotalCell();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal        = numberOfTotalCell + numberOfBoundaryFace;

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
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

	RDouble **residualTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));

	if (grid->GetLevel() != 0) 
	{
		RDouble **rhs = rhsProxy->GetField_UNS();
		for (int m = 0; m < nTurbulentEquation; ++ m)
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
	int nTurbulentEquation = parameters->GetNTurbulentEquation();

	RDouble **residualTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));
	RDouble **rhs = rhsProxy->GetField_UNS();

	for (int m = 0; m < nTurbulentEquation; ++ m)
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
	int nTurbulentEquation = parameters->GetNTurbulentEquation();

	RDouble **residualTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));
	RDouble **rhs = rhsProxy->GetField_UNS();

	for (int m = 0; m < nTurbulentEquation; ++ m)
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
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    for (int m = 0; m < nTurbulentEquation; ++ m)
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
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    RDouble **residualTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));
	using namespace IDX;

    //! For boundary faces, remember re is ghost cell
    for (int i = localStart; i < nMid; ++ i)
    {
        int le,j;
        le = leftCellOfFace[i];
        j = i - localStart;
        for (int m = 0; m < nTurbulentEquation; ++ m)
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

        for (int m = 0; m < nTurbulentEquation; ++ m)
        {
            residualTurbulence[m][le] -= flux[m][j];
            residualTurbulence[m][re] += flux[m][j];
        }
    }
}

void TurbSolverUnstr::ModifyResidualOnWallBoundary(Grid *gridIn)
{
	UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell    = grid->GetNTotalCell();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal        = numberOfTotalCell + numberOfBoundaryFace;

    int modifyTurbulentResidual = GlobalDataBase::GetIntParaFromDB("mod_turb_res");

    //! modify the residuals for the cells next to the wall
    if (modifyTurbulentResidual != 1) return;

    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    BCRecord **bcRecord = grid->GetBCRecord();

    RDouble **residualTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
	int nTurbulentEquation = parameters->GetNTurbulentEquation();

    using namespace PHENGLEI;
    if (nTurbulentEquation == 1)
    {
        for (int iFace = 0; iFace < numberOfBoundaryFace; ++ iFace)
        {
            int le = leftCellOfFace [iFace];
            int re = rightCellOfFace[iFace];
            int bcType = bcRecord[iFace]->GetKey();

            if (IsWall(bcType))
            {
                residualTurbulence[0][le] = 0.0;	
				residualTurbulence[0][re] = 0.0;				// Bell 20130508 add
            }
        }
    }
}

void TurbSolverUnstr::ResetWallScalar(Grid *gridIn)
{
	UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell    = grid->GetNTotalCell();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal        = numberOfTotalCell + numberOfBoundaryFace;

    int      *leftCellOfFace  = grid->GetLeftCellOfFace();
    int      *rightCellOfFace  = grid->GetRightCellOfFace();
    BCRecord **bcRecord  = grid->GetBCRecord();

    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb" ));

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    using namespace PHENGLEI;

    if (nTurbulentEquation == 1)
    {
        for (int iFace = 0; iFace < numberOfBoundaryFace; ++ iFace)
        {
            int le = leftCellOfFace [iFace];
            int re = rightCellOfFace[iFace];
            int bcType = bcRecord[iFace]->GetKey();

            if (IsWall(bcType))
            {
                qTurbulence[0][le] = 0.0;
            }
        }
    }
}

void TurbSolverUnstr::SourceFluxOneEquation(Grid *grid)
{
    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    //int kindOfTurbSource = parameters->GetKindOfTurbSource();
   
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
    //}
    //else if (turbulentSourceKind == 1)
    //{
    //    SourceFluxOneEquationEdwards(grid);
    //}
    //else if (turbulentSourceKind == 2)
    //{
    //    SourceFluxOneEquationNew(grid);
    //}
}

// Bell 20130513 mod
void TurbSolverUnstr::SourceFluxOneEquationOriginal(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalCell    = grid->GetNTotalCell();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal        = numberOfTotalCell + numberOfBoundaryFace;

    RDouble **gradientPrimtiveVariableFaceX = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarFaceX"));
    RDouble **gradientPrimtiveVariableFaceY = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarFaceY"));
    RDouble **gradientPrimtiveVariableFaceZ = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarFaceZ"));

    RDouble *dudx = gradientPrimtiveVariableFaceX[1];
    RDouble *dudy = gradientPrimtiveVariableFaceY[1];
    RDouble *dudz = gradientPrimtiveVariableFaceZ[1];

    RDouble *dvdx = gradientPrimtiveVariableFaceX[2];
    RDouble *dvdy = gradientPrimtiveVariableFaceY[2];
    RDouble *dvdz = gradientPrimtiveVariableFaceZ[2];

    RDouble *dwdx = gradientPrimtiveVariableFaceX[3];
    RDouble *dwdy = gradientPrimtiveVariableFaceY[3];
    RDouble *dwdz = gradientPrimtiveVariableFaceZ[3];

    RDouble **dqdx = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradTurbulentX"));
    RDouble **dqdy = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradTurbulentY"));
    RDouble **dqdz = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradTurbulentZ"));

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

	RDouble *cellCenterX = grid->GetCellCenterX();
	RDouble *cellCenterY = grid->GetCellCenterY();
	RDouble *cellCenterZ = grid->GetCellCenterZ();

	Param_TurbSolverUnstruct *parameters = GetControlParameters();
    RDouble oRefReNumber = parameters->GetoRefReNumber();

    int SAProductType = parameters->GetSAProductType();
	
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
    const RDouble cdes   = - 1.0;
    const RDouble sdilim = 1.0e9;

    RDouble kap2,oKap2;
    kap2  = karm * karm;
    oKap2 = one / kap2;

    RDouble osigma = 1.0 / sigma;
    RDouble cw1    = cb1 * oKap2 + (one + cb2) / sigma;
    RDouble cv13   = cv1 * cv1 * cv1;
    RDouble cw36   = pow(cw3,6);
    RDouble or6    = one / six;

    RDouble ct3 = 1.2;
    RDouble ct4 = 0.5; 

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

	const RDouble EPS = 1.0e-10;

    cb2s = osigma * cb2;
    cw1Kap2 = cw1 * kap2;

	bool outStep = false;
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

        //! sBar = fv2 * nueSA / ( k2 * distance2).
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
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal        = numberOfTotalCell + numberOfBoundaryFace;

    RDouble **gradientPrimtiveVariableFaceX = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarFaceX"));
    RDouble **gradientPrimtiveVariableFaceY = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarFaceY"));
    RDouble **gradientPrimtiveVariableFaceZ = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarFaceZ"));

    RDouble *dudx = gradientPrimtiveVariableFaceX[1];
    RDouble *dudy = gradientPrimtiveVariableFaceY[1];
    RDouble *dudz = gradientPrimtiveVariableFaceZ[1];

    RDouble *dvdx = gradientPrimtiveVariableFaceX[2];
    RDouble *dvdy = gradientPrimtiveVariableFaceY[2];
    RDouble *dvdz = gradientPrimtiveVariableFaceZ[2];

    RDouble *dwdx = gradientPrimtiveVariableFaceX[3];
    RDouble *dwdy = gradientPrimtiveVariableFaceY[3];
    RDouble *dwdz = gradientPrimtiveVariableFaceZ[3];

    RDouble **dqdx = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradTurbulentX"));
    RDouble **dqdy = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradTurbulentY"));
    RDouble **dqdz = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradTurbulentZ"));

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

    RDouble *cellCenterX = grid->GetCellCenterX();
    RDouble *cellCenterY = grid->GetCellCenterY();
    RDouble *cellCenterZ = grid->GetCellCenterZ();

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    RDouble reynolds = parameters->GetRefReNumber();
    RDouble oRefReNumber = parameters->GetoRefReNumber();

    int SAProductType = parameters->GetSAProductType();

    using namespace IDX;

    //von karman constant
    const RDouble karm   = 0.41;
    const RDouble cb1    = 0.1355;
    const RDouble cb2    = 0.622;
    const RDouble cv1    = 7.1;
    const RDouble cw2    = 0.3;
    const RDouble cw3    = 2.0;
    const RDouble sigma  = 2.0/3.0;
    const RDouble cdes   = - 1.0;
    const RDouble sdilim = 1.0e9;

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

    const RDouble EPS = 1.0e-10;

    cb2s = osigma * cb2;
    cw1Kap2 = cw1 * kap2;

    RDouble machinZero = pow(10.0, -machZeroExp + 1);

    bool outStep = false;
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
        //! Here, absorb oRefReNumber into rs(term0): nueSA / ( k2 * distance2).
        RDouble term0 = nueSA * oKap2 * oDistance2 * oRefReNumber;        

        //!sBar = fv2 * nueSA / ( k2 * distance2).
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
        RDouble term1_CFL3D = cb1 * ( one - ft2) * omega;
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

    RDouble **primitiveVariables   = reinterpret_cast<RDouble **> (grid->GetDataPtr("q" ));
	RDouble *viscousLaminar        = reinterpret_cast<RDouble * > (grid->GetDataPtr("visl"));
    RDouble *viscousTurbulence = reinterpret_cast<RDouble * > (grid->GetDataPtr("vist"));
    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb" ));
	RDouble **spectrumTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("spec_turb"));
	RDouble **residualTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));

    RDouble *wallDistance = grid->GetWallDist();
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

    RDouble cblend, rho, ke, kw;
    RDouble sij2, divv, prodk, dissk, prodw, dissw, cdkww, crss, mut;
    RDouble srck, srcw, fskn, fswn, diak, diaw, oork, oorw, volume;
    /*RDouble turbXk, turbXw, turbXk2, turbFk2, turbFw;*/
    RDouble s11, s22, s33, s12, s13, s23, w12, w13, w23;
    RDouble vort2;

    RDouble turbFbetas, turbFbeta, turbBeta, turbAlphaw/*, turbSigd*/;
    //RDouble turbCmu;
    //RDouble** gradientPrimtiveVariableFaceX = reinterpret_cast <RDouble**> (grid->GetDataPtr("gradPrimtiveVarFaceX"));
    //RDouble** gradientPrimtiveVariableFaceY = reinterpret_cast <RDouble**> (grid->GetDataPtr("gradPrimtiveVarFaceY"));
    //RDouble** gradientPrimtiveVariableFaceZ = reinterpret_cast <RDouble**> (grid->GetDataPtr("gradPrimtiveVarFaceZ"));
    //
    //using namespace IDX;
    //RDouble* dudx = gradientPrimtiveVariableFaceX[IU];
    //RDouble* dudy = gradientPrimtiveVariableFaceY[IU];
    //RDouble* dudz = gradientPrimtiveVariableFaceZ[IU];
    //
    //RDouble* dvdx = gradientPrimtiveVariableFaceX[IV];
    //RDouble* dvdy = gradientPrimtiveVariableFaceY[IV];
    //RDouble* dvdz = gradientPrimtiveVariableFaceZ[IV];
    //
    //RDouble* dwdx = gradientPrimtiveVariableFaceX[IW];
    //RDouble* dwdy = gradientPrimtiveVariableFaceY[IW];
    //RDouble* dwdz = gradientPrimtiveVariableFaceZ[IW];

    RDouble *dudx = new RDouble[numberOfTotal];
    RDouble *dudy = new RDouble[numberOfTotal];
    RDouble *dudz = new RDouble[numberOfTotal];
    
    RDouble *dvdx = new RDouble[numberOfTotal];
    RDouble *dvdy = new RDouble[numberOfTotal];
    RDouble *dvdz = new RDouble[numberOfTotal];
    
    RDouble *dwdx = new RDouble[numberOfTotal];
    RDouble *dwdy = new RDouble[numberOfTotal];
    RDouble *dwdz = new RDouble[numberOfTotal];

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

	using namespace IDX;

    //====================
    //=== K-OMEGA MODELS
    //====================
    if (viscousName.substr(0, 6) == "2eq-kw")
    {
        //******************************
        //*** MENTER'S K-OMEGA SST MODEL
        //******************************
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
            //Crossing(grid, cross);

            RDouble *blend = reinterpret_cast<RDouble  *> (grid->GetDataPtr("blend"));
            //Blending(grid, cross, blend);

            for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell)
	        {
                cblend = blend[iCell];
                crss   = cross[iCell];
                rho    = primitiveVariables[IR][iCell];
                ke     = qTurbulence[0][iCell];
                kw     = qTurbulence[1][iCell];
                mut    = viscousTurbulence[iCell];
                volume = vol[iCell];

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
                
                sij2  = two * (s11*s11 + s22*s22 + s33*s33 + two*(s12*s12 + s13*s13 + s23*s23));
                divv  = s11 + s22 + s33;

                if (neasm < 0)
                {
                    if (SSTProductType == 0)
                    {
                        // Boussinesq approximation,full production
                        prodk = mut * (sij2 - two3rd * divv * divv) * oRefReNumber - two3rd * rho * ke * divv;
                    }
                    else
                    {
                        vort2 = four * (w12 * w12 + w13 * w13 + w23 * w23);
                        prodk = mut * (vort2) * oRefReNumber - two3rd * rho * ke * divv;
                    }
                }
                dissk = turbFbetas * SST_betaStar * rho * ke * kw * refReNumber;
                
                //!limit production
                prodk = MIN(prodk, SSTProductLimit * dissk);
                prodw = turbAlphaw * rho / (mut + SMALL) * prodk;
                dissw = turbFbeta * turbBeta * rho * kw * kw * refReNumber;
                cdkww = two * (one - cblend) * rho * KW_sigmaW2 * crss / (kw + SMALL) * oRefReNumber;

                srck  = prodk - dissk;
                srcw  = prodw - dissw + cdkww;

                residualTurbulence[0][iCell] += (srck) * volume;
                residualTurbulence[1][iCell] += (srcw) * volume;

                if (neasm > 0)
                {
                    sengy[iCell] = - (srck) * volume;
                }

                diak = - two * turbFbetas * SST_betaStar * kw;
                diaw = - two * turbFbeta  * turbBeta  * kw - ABS(cdkww) / (rho * kw + SMALL);

                fskn = half * (srck - ABS(srck));
                fswn = half * (srcw - ABS(srcw));

                oork = 1.0 / (rho * ke + SMALL);
                oorw = 1.0 / (rho * kw + SMALL);

                spectrumTurbulence[0][iCell] += - (diak + fskn * oork) * volume;
                spectrumTurbulence[1][iCell] += - (diaw + fswn * oorw) * volume;
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

    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb" ));
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
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal        = numberOfTotalCell + numberOfBoundaryFace;

    RDouble **q      = reinterpret_cast<RDouble **> (grid->GetDataPtr("q" ));
    RDouble *viscousLaminar   = reinterpret_cast<RDouble * > (grid->GetDataPtr("visl"));
    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb" ));
    RDouble *cross = reinterpret_cast<RDouble  *> (grid->GetDataPtr("cross"));
    RDouble *blend = reinterpret_cast<RDouble  *> (grid->GetDataPtr("blend"));

    RDouble *wallDistance = grid->GetWallDist();

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    RDouble refReNumber = parameters->GetRefReNumber();
    RDouble oReynolds = 1.0 / refReNumber;
    RDouble oReynoldsSqr = oReynolds * oReynolds;

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
        //cdkw      = MAX(crossDiff,cdkwmin);

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
    }

    grid->SetGhostCellExceptInterface(blend);
}

FaceProxy * TurbSolverUnstr::CreateFaceProxy(Grid *gridIn)
{
	UnstructGrid *grid = UnstructGridCast (gridIn);
    int numberOfTotalCell    = grid->GetNTotalCell();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal        = numberOfTotalCell + numberOfBoundaryFace;

	FaceProxy *faceProxy = CreateFaceProxyTurbulence(grid);
	faceProxy->SetNext(CreateFaceProxyNS(grid));

	return faceProxy;
}

FaceProxy * TurbSolverUnstr::CreateFaceProxyNS(Grid *gridIn)
{
	UnstructGrid *grid = UnstructGridCast (gridIn);
    int numberOfTotalCell    = grid->GetNTotalCell();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal        = numberOfTotalCell + numberOfBoundaryFace;

    int nl;
    GlobalDataBase::GetData("nl",&nl,PHINT,1);

    int nchem;
    GlobalDataBase::GetData("nchem",&nchem,PHINT,1);

	int numberOfEquation = nl + nchem;

	FaceProxy *faceProxy = new FaceProxy();
	faceProxy->Create(SEGCTION_LENGTH, numberOfEquation);

	return faceProxy;
}

FaceProxy * TurbSolverUnstr::CreateFaceProxyTurbulence(Grid *gridIn)
{
	UnstructGrid *grid = UnstructGridCast (gridIn);
    int numberOfTotalCell    = grid->GetNTotalCell();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal        = numberOfTotalCell + numberOfBoundaryFace;

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

	FaceProxy *faceProxy = new FaceProxy();
	faceProxy->Create(SEGCTION_LENGTH, nTurbulentEquation);

	return faceProxy;
}

TurbFaceValue * TurbSolverUnstr::CreateTurbFaceValue(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotalCell    = grid->GetNTotalCell();
    int numberOfTotal        = numberOfTotalCell + numberOfBoundaryFace;

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    TurbFaceValue *facevar = new TurbFaceValue(nTurbulentEquation, SEGCTION_LENGTH);
	return facevar;
}

void TurbSolverUnstr::ComputeQTurbNodeValue(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    int numberOfTotalCell = grid->GetNTotalCell();
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotal = numberOfTotalCell + numberOfBoundaryFace;
    int numberOfTotalNode = grid->GetNTotalNode();
    int numberOfTotalFace = grid->GetNTotalFace();

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    RDouble **qTurbulenceNode = reinterpret_cast<RDouble **>(grid->GetDataPtr("qTurbNode"));
    RDouble     **qTurbulence = reinterpret_cast<RDouble **>(grid->GetDataPtr("q_turb"));
    int   *nodeValueSliceTurb = reinterpret_cast<int      *>(grid->GetDataPtr("nodeValueSliceTurb"));

    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    int *face2node = grid->GetFace2Node();
    int *nodeNumberOfEachFace = grid->GetNodeNumberOfEachFace();
 
    BCRecord **bcRecord = grid->GetBCRecord();

    for (int iNode = 0; iNode < numberOfTotalNode; ++ iNode)
    {
        nodeValueSliceTurb[iNode] = 0;

        for(int m = 0; m < nTurbulentEquation; ++ m)
        {
            qTurbulenceNode[m][iNode] = 0; 
        }
    }

    //! Boundary faces
    int nodePosition = 0;
    for (int iFace = 0; iFace < numberOfBoundaryFace; ++ iFace)
    {
        int bcType = bcRecord[iFace]->GetKey();
        int le = leftCellOfFace [iFace];
        int re = rightCellOfFace[iFace];
        if (bcType != PHENGLEI::INTERFACE)
        {
            for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
            {
                //! From left
                int point = face2node[nodePosition + jNode];

                for (int m = 0; m < nTurbulentEquation; ++ m)
                {
                    qTurbulenceNode[m][point] += qTurbulence[m][le];
                }            
                nodeValueSliceTurb[point] += 1;

                //! From right
                for (int m = 0; m < nTurbulentEquation; ++ m)
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
                for (int m = 0; m < nTurbulentEquation; ++ m)
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
            for (int m = 0; m < nTurbulentEquation; ++ m)
            {
                qTurbulenceNode[m][point] += qTurbulence[m][le];
            }            
            ++ nodeValueSliceTurb[point];

            for (int m = 0; m < nTurbulentEquation; ++ m)
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
        int bcType = bcRecord[iFace]->GetKey();
        if (bcType == PHENGLEI::SYMMETRY)
        {
            for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
            {
                int point = face2node[nodePosition + jNode];
            
                for (int m = 0; m < nTurbulentEquation; ++ m)
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
        int bcType = bcRecord[iFace]->GetKey();
        if (bcType == PHENGLEI::SYMMETRY)
        {
            for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
            {
                int point = face2node[nodePosition + jNode];
            
                for (int m = 0; m < nTurbulentEquation; ++ m)
                {
                    qTurbulenceNode[m][point] += qTurbulence[m][le];
                }            
                ++ nodeValueSliceTurb[point];

                for (int m = 0; m < nTurbulentEquation; ++ m)
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
        int bcType = bcRecord[iFace]->GetKey();
        if (bcType == PHENGLEI::FARFIELD)
        {
            for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
            {
                int point = face2node[nodePosition + jNode];
            
                for (int m = 0; m < nTurbulentEquation; ++ m)
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
        int bcType = bcRecord[iFace]->GetKey();
        if (bcType == PHENGLEI::FARFIELD)
        {
            for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
            {
                int point = face2node[nodePosition + jNode];
            
                for (int m = 0; m < nTurbulentEquation; ++ m)
                {
                    qTurbulenceNode[m][point] += qTurbulence[m][le];
                }            
                ++ nodeValueSliceTurb[point];

                for (int m = 0; m < nTurbulentEquation; ++ m)
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
        int bcType = bcRecord[iFace]->GetKey();
        if (bcType == PHENGLEI::SOLID_SURFACE)
        {
            for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
            {
                int point = face2node[nodePosition + jNode];
            
                for (int m = 0; m < nTurbulentEquation; ++ m)
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
        int bcType = bcRecord[iFace]->GetKey();
        if (bcType == PHENGLEI::SOLID_SURFACE)
        {
            for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
            {
                int point = face2node[nodePosition + jNode];
            
                for (int m = 0; m < nTurbulentEquation; ++ m)
                {
                    qTurbulenceNode[m][point] += qTurbulence[m][le];
                }            
                ++ nodeValueSliceTurb[point];

                for (int m = 0; m < nTurbulentEquation; ++ m)
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
        int *cellNumberOfInterPoint = interPointInfor->GetCellNumberOfInterPoint();
        int *totalZonesOfInterPoint = interPointInfor->GetTotalZonesOfInterPoint();
        
        int globalPoint;
        for (int iPoint = 0; iPoint < numberOfInterPoint; ++ iPoint)
        {
            globalPoint = interPoint2GlobalPoint[iPoint];
            if (labelOfInterPoint[iPoint] != 1)
            {
                for (int m = 0; m < nTurbulentEquation;++ m)
                {    
                    globalPoint = interPoint2GlobalPoint[iPoint];
                    qTurbulenceNode[m][globalPoint] = 0;
                }
            }
        }
        RDouble **qTurbulenceInterPoint = reinterpret_cast<RDouble **> (grid->GetDataPtr("turb::qInterpoint"));
        for (int iPoint = 0; iPoint < numberOfInterPoint; ++ iPoint)
        {
            for (int m = 0; m < nTurbulentEquation; ++ m)
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
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    RDouble **qTurbulenceNode = reinterpret_cast<RDouble **>(grid->GetDataPtr("qTurbNode"));
    int   *nodeValueSliceTurb = reinterpret_cast<int      *>(grid->GetDataPtr("nodeValueSliceTurb"));

    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        if (nodeValueSliceTurb[iNode] != 0)
        {
            for (int m = 0; m < nTurbulentEquation; m++)
            {
                qTurbulenceNode[m][iNode] /= nodeValueSliceTurb[iNode];
            }
        }
        else
        {
            for (int m = 0; m < nTurbulentEquation; m++)
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
    int numberOfTotal = numberOfTotalCell + numberOfBoundaryFace;

    RDouble *xfc  = grid->GetFaceCenterX();
    RDouble *yfc  = grid->GetFaceCenterY();
    RDouble *zfc  = grid->GetFaceCenterZ();
    RDouble *cellCenterX  = grid->GetCellCenterX();
    RDouble *cellCenterY  = grid->GetCellCenterY();
    RDouble *cellCenterZ  = grid->GetCellCenterZ();
    RDouble *xfn  = grid->GetFaceNormalX();
    RDouble *yfn  = grid->GetFaceNormalY();
    RDouble *zfn  = grid->GetFaceNormalZ();
    RDouble *area = grid->GetFaceArea();
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

	TurbFaceValue *faceVariable = faceProxy->GetTurbFaceValue();

	RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));

    RDouble *viscousLaminar = reinterpret_cast<RDouble *> (grid->GetDataPtr("visl"));
	RDouble *viscousTurbulence = reinterpret_cast<RDouble *> (grid->GetDataPtr("vist"));

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

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

    BCRecord  **bcRecord = grid->GetBCRecord();

    int nMid = numberOfBoundaryFace;
    if (localEnd < numberOfBoundaryFace) nMid = localEnd;

    if (localStart < numberOfBoundaryFace)
    {
        for (int iFace = localStart; iFace < nMid; ++ iFace)
        {
            int bcType = bcRecord[iFace]->GetKey();
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

    if (nTurbulentEquation == 1)
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
    int numberOfTotalCell = grid->GetNTotalCell();
    int numberOfTotal = numberOfTotalCell + numberOfBoundaryFace;

	int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

	BCRecord **bcRecord  = grid->GetBCRecord();

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

	RDouble **prim = faceVariable->prim;
	RDouble **mlt  = faceVariable->mlt;

    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));

    RDouble *viscousLaminar = reinterpret_cast<RDouble *> (grid->GetDataPtr("visl"));
    RDouble **qPrimitiveVariable = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));

    RDouble **residualTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("res_turb"));
    RDouble **spectrumTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("spec_turb"));
    RDouble **matrixTurbulenceLeft = reinterpret_cast<RDouble **> (grid->GetDataPtr("mat_turbl"));
    RDouble **matrixTurbulenceRight = reinterpret_cast<RDouble **> (grid->GetDataPtr("mat_turbr"));

    Param_TurbSolverUnstruct *parameters = GetControlParameters();

	RDouble *freeStreamTurbVar =  parameters->GetFreeStreamTurbVar();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();
    RDouble oRefReNumber = parameters->GetoRefReNumber();

	const RDouble cb2    = 0.622;
	const RDouble sigma  = 2.0/3.0;
    const RDouble osigma = 1.0 / sigma;

    RDouble *dfd1  = new RDouble[nTurbulentEquation];
    RDouble *dfd2  = new RDouble[nTurbulentEquation];
    RDouble *dfdn  = new RDouble[nTurbulentEquation];

    RDouble *f1    = new RDouble[nTurbulentEquation];
    RDouble *f2    = new RDouble[nTurbulentEquation];
    RDouble *fMid  = new RDouble[nTurbulentEquation];
    bool isFineGrid = grid->IsFinestGrid();

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

        RDouble dxnL  = nxs * dL - dxL;
        RDouble dynL  = nys * dL - dyL;
        RDouble dznL  = nzs * dL - dzL;

        RDouble dxnR  = nxs * dR - dxR;
        RDouble dynR  = nys * dR - dyR;
        RDouble dznR  = nzs * dR - dzR;

        //! quantities at points 1 and 2
        for (int m = 0; m < nTurbulentEquation; ++ m)
        {
            f1[m] = qTurbulence[m][le];
            f2[m] = qTurbulence[m][re];
        }

        for (int m = 0; m < nTurbulentEquation; ++ m)
        {
            fMid[m] = half * (f1[m] + f2[m]);
        }

        for (int m = 0; m < nTurbulentEquation; ++ m)
        {
            dfdn[m] = 0.0;
        }

        if (angle1 > 0.0 && angle2 > 0.0 && ABS(dL) > TINY && ABS(dR) > TINY)
        {
            for (int m = 0; m < nTurbulentEquation; ++ m)
            {
                dfd1[m] = (f1[m] - fMid[m]) / dL;
            }

            for (int m = 0; m < nTurbulentEquation; ++ m)
            {
                dfd2[m] = (f2[m] - fMid[m]) / dR;
            }

            RDouble dtmp = dL * dL + dR * dR;
            RDouble weightL = dL * dL / dtmp;
            RDouble weightR = 1.0 - weightL;

            for (int m = 0; m < nTurbulentEquation; ++ m)
            {
                dfdn[m] = weightL * dfd1[m] + weightR * dfd2[m];
            }

            if (iFace < numberOfBoundaryFace)
            {
                int bcType = bcRecord[iFace]->GetKey();

                if (bcType != PHENGLEI::INTERFACE && bcType != PHENGLEI::SYMMETRY)
                {
                    RDouble normalDist = (cellCenterX[re] - cellCenterX[le]) * nxs +
                                         (cellCenterY[re] - cellCenterY[le]) * nys +
                                         (cellCenterZ[re] - cellCenterZ[le]) * nzs;

                    for (int m = 0; m < nTurbulentEquation; ++ m)
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
            for (int m = 0; m < nTurbulentEquation; ++ m)
            {
                flux[m][jFace] = - oRefReNumber * mlt[m][jFace] * dfdn[m] * area[iFace];
            }
            
            for (int m = 0; m < nTurbulentEquation; ++ m)
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
	int nTurbulentEquation = parameters->GetNTurbulentEquation();

	if (nTurbulentEquation == 1)
	{
		SpectrumRadiusOfOneEquation(grid);
	}
	else if (nTurbulentEquation == 2)
	{
		SpectrumRadiusOfTwoEquation(grid);
	}
}

void TurbSolverUnstr::SpectrumRadiusOfOneEquation(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfBoundaryFace = grid->GetNBoundFace();
    int numberOfTotalCell    = grid->GetNTotalCell();
    int numberOfTotalFace    = grid->GetNTotalFace();

    RDouble **q         = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
    RDouble *viscousLaminar = reinterpret_cast<RDouble *> (grid->GetDataPtr("visl"));

    RDouble *dt = reinterpret_cast<RDouble *> (grid->GetDataPtr("dt"));

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    RDouble turbCFLScale = parameters->GetTurbCFLScale();

    RDouble **spectrumTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("spec_turb"));
    RDouble **matrixTurbulenceLeft = reinterpret_cast<RDouble **> (grid->GetDataPtr("mat_turbl"));
    RDouble **matrixTurbulenceRight = reinterpret_cast<RDouble **> (grid->GetDataPtr("mat_turbr"));

    RDouble *gama = reinterpret_cast<RDouble *> (grid->GetDataPtr("gama"));

    RDouble oRefReNumber = parameters->GetoRefReNumber();

    //!initialize some constants
    RDouble sigma  = 2.0 / 3.0;
    RDouble osigma = one / sigma;
    const RDouble cb2 = 0.622;

    //RDouble xisens = one;
    //RDouble alpha  = one;
    //RDouble eps2   = one;

    //RDouble rhoL, uL, vL, wL, pL, vnL, orL;
    //RDouble rhoR, uR, vR, wR, pR, vnR, orR;
    //RDouble volume, oVolume, ns2;
    //RDouble nueL, nueR, nueffL, nueffR, nueff;

    int *leftCellOfFace  = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

    RDouble  *vol  = grid->GetCellVolume();
    RDouble  *xfn  = grid->GetFaceNormalX();
    RDouble  *yfn  = grid->GetFaceNormalY();
    RDouble  *zfn  = grid->GetFaceNormalZ();
    RDouble  *area = grid->GetFaceArea();
    BCRecord **bcRecord = grid->GetBCRecord();

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
        RDouble pL = q[IP][le];
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
        RDouble pR = q[IP][re];
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
    int numberOfTotalCell = grid->GetNTotalCell();
    int numberOfTotalFace = grid->GetNTotalFace();

	Param_TurbSolverUnstruct *parameters = GetControlParameters();
	int nTurbulentEquation = parameters->GetNTurbulentEquation();
	RDouble turbCFLScale = parameters->GetTurbCFLScale();

    RDouble oRefReNumber = parameters->GetoRefReNumber();

    RDouble **q = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
	RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));

	RDouble *viscousLaminar = reinterpret_cast<RDouble * > (grid->GetDataPtr("visl"));
	RDouble *viscousTurbulence = reinterpret_cast<RDouble * > (grid->GetDataPtr("vist"));
	
	RDouble *dt = reinterpret_cast<RDouble * > (grid->GetDataPtr("dt"));

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

    vector<RDouble> work(nTurbulentEquation);

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

        for (int m = 0; m < nTurbulentEquation; ++ m)
        {
		    spectrumTurbulence[m][le] += area[iFace] * half * absVn;
		    spectrumTurbulence[m][re] += area[iFace] * half * absVn;

            matrixTurbulenceRight[m][iFace]  = area[iFace] * half * (  vnR - absVn);
		    matrixTurbulenceLeft[m][iFace]   = area[iFace] * half * (- vnL - absVn);
        }
	}

    RDouble *blend = reinterpret_cast<RDouble  *> (grid->GetDataPtr("blend"));

    RDouble turbSigk, turbSigw;

    RDouble KW_sigmaK1 = parameters->GetKW_sigmaK1();
    RDouble KW_sigmaW1 = parameters->GetKW_sigmaW1();
    RDouble KW_sigmaK2 = parameters->GetKW_sigmaK2();
    RDouble KW_sigmaW2 = parameters->GetKW_sigmaW2();

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
	    work[0] = (mul + mut * turbSigk) / rho * ns2 * oVolume;
	    work[1] = (mul + mut * turbSigw) / rho * ns2 * oVolume;

        for (int m = 0; m < nTurbulentEquation; ++ m)
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
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    RDouble **qTurbulentVarables = reinterpret_cast <RDouble **> (grid->GetDataPtr("q_turb"));

    RDouble **qL = faceProxy->GetQL();
    RDouble **qR = faceProxy->GetQR();

    for (int iFace = localStart; iFace <localEnd; ++ iFace)
    {
        int le, re, jFace;
        le = leftCellOfFace [iFace];
        re = rightCellOfFace[iFace];
        jFace  = iFace - localStart;

        for (int m = 0; m < nTurbulentEquation; ++ m)
        {
            qL[m][jFace] = qTurbulentVarables[m][le];
            qR[m][jFace] = qTurbulentVarables[m][re];
        }
    }

    //! Boundary fix.
    BCRecord **bcRecord = grid->GetBCRecord();
    int nBoundFace = grid->GetNBoundFace();
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();
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

    for (int iFace = localStart; iFace < nMid; ++ iFace)
    {
        int BCType = bcRecord[iFace]->GetKey();
        int jFace = iFace - localStart;

        if (BCType == PHENGLEI::INTERFACE)
        {
            continue;
        }

        int le = leftCellofFace[iFace];
        int re = rightCellofFace[iFace];

        if (BCType == PHENGLEI::SYMMETRY)
        {
            for (int m = 0; m < nTurbulentEquation; ++ m)
            {
                qR[m][jFace] = qL[m][jFace];
            }
        }
        else
        {
            for (int m = 0; m < nTurbulentEquation; ++ m)
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
    BCRecord **bcRecord = grid->GetBCRecord();
    int nBoundFace = grid->GetNBoundFace();
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();
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

    for (int iFace = localStart; iFace < nMid; ++ iFace)
    {
        int BCType = bcRecord[iFace]->GetKey();
        int jFace = iFace - localStart;

        if (BCType == PHENGLEI::INTERFACE)
        {
            continue;
        }

        int le = leftCellofFace[iFace];
        int re = rightCellofFace[iFace];

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
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

	RDouble **qL   = faceProxy->GetQL();
	RDouble **qR   = faceProxy->GetQR();
	RDouble **flux = faceProxy->GetFlux();

	FaceProxy *faceProxyNS = faceProxy->GetNext();

	RDouble **qPrimitiveVariableL = faceProxyNS->GetQL();
	RDouble **qPrimitiveVariableR = faceProxyNS->GetQR();

    RDouble *primitiveVariableL = new RDouble[nTurbulentEquation];
    RDouble *primitiveVariableR = new RDouble[nTurbulentEquation];
    RDouble *f = new RDouble[nTurbulentEquation];

	RDouble rL,uL,vL,wL,vnL,vnLL,rR,uR,vR,wR,vnR,vnRR;
    RDouble rhoFlag = 0, rhoCoefficientL, rhoCoefficientR;

	using namespace IDX;

    if (nTurbulentEquation == 1)
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

		for (int m = 0; m < nTurbulentEquation; ++ m)
		{
			primitiveVariableL[m] = rhoCoefficientL * qL[m][jFace];
			primitiveVariableR[m] = rhoCoefficientR * qR[m][jFace];
		}

        //! Consider vgn here!
		for (int m = 0; m < nTurbulentEquation; ++ m)
		{
			f[m] = (vnLL * primitiveVariableL[m] + vnRR * primitiveVariableR[m]);
		}

        for (int m = 0; m < nTurbulentEquation; ++ m)
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
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

	PHSPACE::ComputeResidualonGrid(grid, actkey, residualTurbulence, nTurbulentEquation);
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

//void TurbSolverUnstr::CompressDQ(DataContainer *&dataContainer, FieldProxy *fieldProxy, Grid *gridIn, const int &zoneIndex, const int &neighborZoneIndex, const int &nEquation)
//{
//    UnstructGrid * grid = UnstructGridCast(gridIn);
//
//    InterfaceInfo * interfaceInformation = grid->GetInterfaceInfo();
//    if ( !interfaceInformation ) return;
//
//    int iNeighborZone                            = interfaceInformation->FindIthNeighbor(neighborZoneIndex);
//    int interfaceNumberBetweenTwoNeighboringZone = interfaceInformation->GetNIFaceOfNeighbor(iNeighborZone);
//    int * interfaceIndexContainerForSend         = interfaceInformation->GetFaceIndexForSend(iNeighborZone);
//    dataContainer->moveto_begin();
//
//    RDouble **dq = fieldProxy->GetField_UNS();
//
//    for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; -- iGhostLayer)
//    {
//        for ( int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace )
//        {
//            int sourceCell;
//            int iFace = interfaceIndexContainerForSend[ iLocalFace ];
//            grid->GetSourceIndex( iFace, iGhostLayer + 1, sourceCell );
//            for ( int m = 0; m < nEquation; ++ m )
//            {
//                PHWrite(dataContainer, dq[m][sourceCell]);
//            }
//        }
//    }
//}
//
//void TurbSolverUnstr::DecompressDQ(DataContainer *&dataContainer, FieldProxy *fieldProxy, Grid *gridIn, const int &zoneIndex, const int &neighborZoneIndex, const int &nEquation)
//{
//    UnstructGrid * grid = UnstructGridCast(gridIn);
//
//    InterfaceInfo * interfaceInformation = grid->GetInterfaceInfo();
//    if ( !interfaceInformation ) return;
//
//    int iNeighborZone                            = interfaceInformation->FindIthNeighbor( neighborZoneIndex );
//    int interfaceNumberBetweenTwoNeighboringZone = interfaceInformation->GetNIFaceOfNeighbor(iNeighborZone);
//    int * interfaceIndexContainerForReceive      = interfaceInformation->GetFaceIndexForRecv(iNeighborZone);
//
//    dataContainer->moveto_begin();
//    RDouble **dq = fieldProxy->GetField_UNS();
//
//    for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; -- iGhostLayer)
//    {
//        for ( int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace )
//        {
//            int iFace = interfaceIndexContainerForReceive[ iLocalFace ];
//            int targetCell;
//            grid->GetTargetIndex( iFace, iGhostLayer + 1, targetCell );
//
//            for ( int m = 0; m < nEquation; ++ m )
//            {
//                PHRead(dataContainer, dq[m][targetCell]);
//            }
//        }
//    }
//}

void TurbSolverUnstr::CreateNodeVariable(Grid *gridIn, RDouble **nodeVariable)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int numberOfTotalNode = grid->GetNTotalNode();

    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    int numberOfEquation = nTurbulentEquation;

    for (int m = 0; m < numberOfEquation; ++ m)
    {
        CompNodeVar(grid, nodeVariable[m], qTurbulence[m]);
    }
}

void TurbSolverUnstr::Boundary(Grid *gridIn)
{
	UnstructGrid *grid = UnstructGridCast(gridIn);

    int *leftCellOfFace  = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
	RDouble *wallDistance = grid->GetWallDist();

    RDouble  *nxs  = grid->GetFaceNormalX();
    RDouble  *nys  = grid->GetFaceNormalY();
    RDouble  *nzs  = grid->GetFaceNormalZ();

    BCRecord **bcRecord = grid->GetBCRecord();
    int numberOfBoundaryFace = grid->GetNBoundFace();

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    RDouble *freeStreamTurbVar =  parameters->GetFreeStreamTurbVar();
    RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
	RDouble **q      = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
	RDouble *gamma = reinterpret_cast<RDouble *> (grid->GetDataPtr("gama"));
    RDouble *viscousLaminar = reinterpret_cast<RDouble  *> (grid->GetDataPtr("visl"));
    RDouble *viscousTurbulence = reinterpret_cast<RDouble  *> (grid->GetDataPtr("vist"));

	int nTurbulentEquation = parameters->GetNTurbulentEquation();
	RDouble freeStreamViscosity = parameters->GetFreeStreamViscosity();

    RDouble refReNumber = parameters->GetRefReNumber();
    RDouble refGama = parameters->GetRefGama();

	using namespace IDX;
    using namespace PHENGLEI;

	RDouble *primitiveVariableInflow = reinterpret_cast<RDouble *> (GlobalDataBase::GetDataPtr("prim_inf"));
	RDouble roo = primitiveVariableInflow[IR];
	RDouble uoo = primitiveVariableInflow[IU];
	RDouble voo = primitiveVariableInflow[IV];
	RDouble woo = primitiveVariableInflow[IW];
	RDouble poo = primitiveVariableInflow[IP];

    RDouble gamm1 = refGama - 1.0;

    int count = 0;

    for (int iFace = 0; iFace < numberOfBoundaryFace; ++ iFace)
    {
        int le = leftCellOfFace [iFace];
        int re = rightCellOfFace[iFace];
        int bcType = bcRecord[iFace]->GetKey();

        if (bcType < 0)
        {
			//! This statement is very important, and related to the way this part is written.
			continue;
        }
        else if (bcType == EXTRAPOLATION)
        {
            viscousTurbulence[re] = viscousTurbulence[le];
            for (int m = 0; m < nTurbulentEquation; ++ m)
            {
                qTurbulence[m][re] = qTurbulence[m][le];
            }
        }
        else if (IsWall(bcType))
        {
            viscousTurbulence[re] = -viscousTurbulence[le];
            if (nTurbulentEquation == 1)
            {
                qTurbulence[0][re] = - qTurbulence[0][le];

				//qTurbulence[0][re] = 0.0;
            }
            else if (nTurbulentEquation == 2)
            {
                RDouble beta1 = 0.075;
                RDouble rho   = q[IR][le];
                RDouble dist  = wallDistance[le];
                RDouble dist2 = dist * dist;

				RDouble kwWall = 60.0 * viscousLaminar[le] / (rho * beta1 * dist2 * refReNumber * refReNumber);

			    qTurbulence[IKE][re] = -qTurbulence[IKE][le];
                //qTurbulence[IKW][re] = MAX(two * kw_wall - qTurbulence[IKW][le],SMALL);
				qTurbulence[IKW][re] = two * kwWall - qTurbulence[IKW][le];
            }
        }
        else if (bcType == SYMMETRY)
        {
            viscousTurbulence[re] = viscousTurbulence[le];
            for (int m = 0; m < nTurbulentEquation; ++ m)
            {
                qTurbulence[m][re] = qTurbulence[m][le];
            }
        }
        else if (bcType == FARFIELD)
        {
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
					for (int m = 0; m < nTurbulentEquation; ++ m)
					{
						qTurbulence[m][re] = qTurbulence[m][le];
					}
                }
                else
                {
					//! inlet
                    viscousTurbulence[re] = freeStreamViscosity;
					for (int m = 0; m < nTurbulentEquation; ++ m)
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
				for (int m = 0; m < nTurbulentEquation; ++ m)
				{
					qTurbulence[m][re] = qTurbulence[m][le];
				}
            }
            else
            {
				//! inlet
                viscousTurbulence[re] = freeStreamViscosity;
				for (int m = 0; m < nTurbulentEquation; ++ m)
				{
					qTurbulence[m][re] = freeStreamTurbVar[m];
				}
			}

        }
        else if (bcType == INFLOW)
        {
            viscousTurbulence[re] = freeStreamViscosity;
			for (int m = 0; m < nTurbulentEquation; ++ m)
			{
				qTurbulence[m][re] = freeStreamTurbVar[m];
			}
        }
        else if (bcType == PRESSURE_INLET)
        {
            viscousTurbulence[re] = freeStreamViscosity;
            for (int m = 0; m < nTurbulentEquation; ++ m)
            {
                qTurbulence[m][re] = freeStreamTurbVar[m];
            }
        }
        else if (bcType == MASS_FLOW_INLET)
        {
            viscousTurbulence[re] = freeStreamViscosity;
            for (int m = 0; m < nTurbulentEquation; ++ m)
            {
                qTurbulence[m][re] = freeStreamTurbVar[m];
            }
        }
        else if (bcType == OUTFLOW)
        {
            viscousTurbulence[re] = viscousTurbulence[le];
            for (int m = 0; m < nTurbulentEquation; ++ m)
            {
                qTurbulence[m][re] = qTurbulence[m][le];
            }
        }
        else if (bcType == PRESSURE_OUTLET)
        {
            viscousTurbulence[re] = viscousTurbulence[le];
            for (int m = 0; m < nTurbulentEquation; ++ m)
            {
                qTurbulence[m][re] = qTurbulence[m][le];
            }
        }
        else if (bcType == MASS_FLOW_OUTLET)
        {
            viscousTurbulence[re] = viscousTurbulence[le];
            for (int m = 0; m < nTurbulentEquation; ++ m)
            {
                qTurbulence[m][re] = qTurbulence[m][le];
            }
        }
        else
        {
            viscousTurbulence[re] = viscousTurbulence[le];
            for (int m = 0; m < nTurbulentEquation; ++ m)
            {
                qTurbulence[m][re] = qTurbulence[m][le];
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
    int monitorViscousTurbulentMax = GlobalDataBase::GetIntParaFromDB("monitor_vistmax");

    string viscousName  = parameters->GetViscousName();
	RDouble refReNumber = parameters->GetRefReNumber();
    RDouble reynoldsSquare = refReNumber * refReNumber;

    RDouble viscousTurbulenceMax = 0.0;
    RDouble viscousTurbulenceMin = 1.0e30;
	RDouble logSmall = log(SMALL);

    RDouble part1, part2, arg2, f2;
    RDouble ld, ld3, fv1;
    RDouble rho, ke, kw, ds, ds2;
    //RDouble s11, s22, s33, s12, s13, s23, sii, divv;

	int iMax = 0;
    int iMin = 0;

    RDouble *cellCenterX = grid->GetCellCenterX();
    RDouble *cellCenterY = grid->GetCellCenterY();
    RDouble *cellCenterZ = grid->GetCellCenterZ();

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

		    if (iCell < numberOfTotalCell && viscousTurbulenceMax < viscousTurbulence[iCell])		// Bell 20130527 add iCell <numberOfTotalCell
		    {
			    viscousTurbulenceMax = viscousTurbulence[iCell];
			    iMax = iCell;
		    }
		    if (iCell < numberOfTotalCell && viscousTurbulenceMin > viscousTurbulence[iCell])		// Bell 20130527 add iCell <numberOfTotalCell
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

	if (monitorViscousTurbulentMax)
	{
		grid->CompareMaxMinValue(viscousTurbulenceMax, 1);
		grid->CompareMaxMinValue(viscousTurbulenceMin, 2);

		ostringstream oss;
		cout << setprecision(6);
		cout << setiosflags(ios::scientific);
		oss  << "Max && min viscosity: " << viscousTurbulenceMax << " " << viscousTurbulenceMin << " " << "\n";
		PrintToWindow(oss);
	}

    //int nTurbulentEquation = parameters->GetNTurbulentEquation();
    //if (nTurbulentEquation == 2)
    //{
    //    grid->SetGhostCell(viscousTurbulence);
    //}
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
		coarseGrid->AllocateWalldist(0);
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
    int nTurbulentEquation = parameters->GetNTurbulentEquation();
	RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
    PHSPACE::UploadInterfaceValue(grid, qTurbulence, "turb::q", nTurbulentEquation);

    RDouble* viscousTurbulence = reinterpret_cast<RDouble*> (grid->GetDataPtr("vist"));
    PHSPACE::UploadInterfaceValue(grid, viscousTurbulence, "vist");

    //RDouble* gradientTurbulenceX = reinterpret_cast<RDouble*> (grid->GetDataPtr("gradTurbulentX"));
    //PHSPACE::UploadInterfaceValue(grid, viscousTurbulence, "turb::dqdx");
    //
    //RDouble* gradientTurbulenceY = reinterpret_cast<RDouble*> (grid->GetDataPtr("gradTurbulentY"));
    //PHSPACE::UploadInterfaceValue(grid, viscousTurbulence, "turb::dqdy");
    //
    //RDouble* gradientTurbulenceZ = reinterpret_cast<RDouble*> (grid->GetDataPtr("gradTurbulentZ"));
    //PHSPACE::UploadInterfaceValue(grid, viscousTurbulence, "turb::dqdz");

    //if (nTurbulentEquation == 2)
    //{
    //    RDouble *blend = reinterpret_cast<RDouble *> (grid->GetDataPtr("blend" ));
    //    PHSPACE::UploadInterfaceValue(grid, blend, "blend");
    //}
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
    int nTurbulentEquation = parameters->GetNTurbulentEquation();
    RDouble **qTurbulenceNode = reinterpret_cast<RDouble **> (grid->GetDataPtr("qTurbNode"));
    PHSPACE::UploadInterpointValue(grid, qTurbulenceNode, "qTurbNode", nTurbulentEquation);
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
    int nTurbulentEquation = parameters->GetNTurbulentEquation();
    RDouble **qTurbulenceNode = reinterpret_cast<RDouble **> (grid->GetDataPtr("qTurbNode"));
    PHSPACE::DownloadInterpointValue(grid, qTurbulenceNode, "qTurbNode", nTurbulentEquation);
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
                for (int ipointLocal = 0; ipointLocal < numberOfInterpointsForNeighbor; ++ ipointLocal )
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
    int nTurbulentEquation = parameters->GetNTurbulentEquation();
	RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
    PHSPACE::DownloadInterfaceValue(grid, qTurbulence, "turb::q", nTurbulentEquation); 

    RDouble *viscousTurbulence = reinterpret_cast<RDouble *> (grid->GetDataPtr("vist"));
    PHSPACE::DownloadInterfaceValue(grid, viscousTurbulence, "vist");
 
    //RDouble* gradientTurbulenceX = reinterpret_cast<RDouble*> (grid->GetDataPtr("gradTurbulentX"));
    //PHSPACE::DownloadInterfaceValue(grid, viscousTurbulence, "turb::dqdx");
    //
    //RDouble* gradientTurbulenceY = reinterpret_cast<RDouble*> (grid->GetDataPtr("gradTurbulentY"));
    //PHSPACE::DownloadInterfaceValue(grid, viscousTurbulence, "turb::dqdy");
    //
    //RDouble* gradientTurbulenceZ = reinterpret_cast<RDouble*> (grid->GetDataPtr("gradTurbulentZ"));
    //PHSPACE::DownloadInterfaceValue(grid, viscousTurbulence, "turb::dqdz");

    //if (nTurbulentEquation == 2)
    //{
    //    RDouble *blend = reinterpret_cast<RDouble *> (grid->GetDataPtr("blend" ));
    //    PHSPACE::DownloadInterfaceValue(grid, blend, "blend");
    //}
}

//! need further compllished!!!
void TurbSolverUnstr::UploadOversetData(ActionKey *actkey)
{
    int level = actkey->level;
    UnstructGrid *grid = reinterpret_cast<UnstructGrid *>(GetGrid(level));

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

	RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
    PHSPACE::UploadOversetValue(grid, qTurbulence, "turb::q", nTurbulentEquation);
}

//! need further compllished!!!
void TurbSolverUnstr::DownloadOversetData(ActionKey *actkey)
{
    int level = actkey->level;
    UnstructGrid *grid = reinterpret_cast<UnstructGrid *>(GetGrid(level));

    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

	RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));

    PHSPACE::DownloadOversetValue(grid, qTurbulence, "turb::q", nTurbulentEquation);
}

void TurbSolverUnstr::PrepareOversetInterfaceData(Data_ParamFieldSuite *datastore, InterfaceDataProxy *interfaceDataProxy)
{
    Param_TurbSolverUnstruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

	RDouble **qTurbulence = reinterpret_cast<RDouble **> (datastore->GetDataPtr("turb::q"));

	vector<RDouble **> &vectorData = interfaceDataProxy->GetVectorData();
	vector<int>   &vectorDimension = interfaceDataProxy->GetVectorDimension();

	vectorDimension.push_back(nTurbulentEquation);
	vectorData.push_back(qTurbulence);
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

	RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb" ));
	RDouble **qNS = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));

	RDouble *viscousLaminar = reinterpret_cast<RDouble *> (grid->GetDataPtr("visl"));
	RDouble *viscousTurbulence = reinterpret_cast<RDouble *> (grid->GetDataPtr("vist"));

	RDouble **Q1;
	RDouble **Q2;
	//Q1 = qTurbulence;
	//Q2 = qNS;

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
		
		//int targetNumberOfBoundaryFace = targetGrid->GetNBoundFace();

		//int *targetInterFace2ZoneID   = targetInterfaceInformation->GetInterFace2ZoneID();
		//int *interFace2InterFaceID_t  = targetInterfaceInformation->GetInterFace2InterFaceID();
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

	RDouble **q     = reinterpret_cast<RDouble **> (grid->GetDataPtr("q" ));
	RDouble *viscousLaminar  = reinterpret_cast<RDouble * > (grid->GetDataPtr("visl"));
	RDouble *viscousTurbulence = reinterpret_cast<RDouble * > (grid->GetDataPtr("vist"));

    RDouble **gradientPrimtiveVariableFaceX = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarFaceX"));
    RDouble **gradientPrimtiveVariableFaceY = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarFaceY"));
    RDouble **gradientPrimtiveVariableFaceZ = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarFaceZ"));

    RDouble *dudx = gradientPrimtiveVariableFaceX[1];
    RDouble *dudy = gradientPrimtiveVariableFaceY[1];
    RDouble *dudz = gradientPrimtiveVariableFaceZ[1];

    RDouble *dvdx = gradientPrimtiveVariableFaceX[2];
    RDouble *dvdy = gradientPrimtiveVariableFaceY[2];
    RDouble *dvdz = gradientPrimtiveVariableFaceZ[2];

    RDouble *dwdx = gradientPrimtiveVariableFaceX[3];
    RDouble *dwdy = gradientPrimtiveVariableFaceY[3];
    RDouble *dwdz = gradientPrimtiveVariableFaceZ[3];

	RDouble oRefReNumber = parameters->GetoRefReNumber();
	const RDouble karm   = 0.41;
	const RDouble kap2   = karm * karm;

	using namespace IDX;
	for (int iCell = 0; iCell < numberOfTotalCell; ++ iCell) 
	{
		RDouble laminarViscosity        = viscousLaminar[iCell];
		RDouble turbulentViscosity = viscousTurbulence[iCell];
		RDouble density                 = q[IR][iCell];

		RDouble wallDistance       = wallDist[iCell];
		RDouble wallDistance2           = SQR(wallDistance);
		RDouble lengthRans              = wallDistance;
		RDouble lengthLes          = lengthOfLES[iCell];

		RDouble sumOfGradient2     = SQR(dudx[iCell], dudy[iCell], dudz[iCell]) + 
			SQR(dvdx[iCell], dvdy[iCell], dvdz[iCell]) + 
			SQR(dwdx[iCell], dwdy[iCell], dwdz[iCell]);

		RDouble nueff = oRefReNumber * (laminarViscosity + turbulentViscosity) / density;
		RDouble rd    = nueff / (kap2 * wallDistance2 * MAX(sqrt(sumOfGradient2), 1.0e-20));
		RDouble fd    = 1.0 - tanh(POWER3(8.0 * rd));

		DDESLength[iCell]  =  lengthRans  -  fd * MAX(0.0, (lengthRans - lengthLes));
	}

	this->SetDESLength(DDESLength);

	delete []lengthOfLES;
	delete []lowReynoldsNumberCorrection;
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

	RDouble **q     = reinterpret_cast<RDouble **> (grid->GetDataPtr("q" ));
	RDouble *viscousLaminar  = reinterpret_cast<RDouble * > (grid->GetDataPtr("visl"));
	RDouble *viscousTurbulence = reinterpret_cast<RDouble * > (grid->GetDataPtr("vist"));

    RDouble **gradientPrimtiveVariableFaceX = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarFaceX"));
    RDouble **gradientPrimtiveVariableFaceY = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarFaceY"));
    RDouble **gradientPrimtiveVariableFaceZ = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarFaceZ"));

    RDouble *dudx = gradientPrimtiveVariableFaceX[1];
    RDouble *dudy = gradientPrimtiveVariableFaceY[1];
    RDouble *dudz = gradientPrimtiveVariableFaceZ[1];

    RDouble *dvdx = gradientPrimtiveVariableFaceX[2];
    RDouble *dvdy = gradientPrimtiveVariableFaceY[2];
    RDouble *dvdz = gradientPrimtiveVariableFaceZ[2];

    RDouble *dwdx = gradientPrimtiveVariableFaceX[3];
    RDouble *dwdy = gradientPrimtiveVariableFaceY[3];
    RDouble *dwdz = gradientPrimtiveVariableFaceZ[3];

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
		RDouble fb   = MIN( 2.0 * exp(- 9.0 * alf * alf),  1.0 );

		RDouble sumOfGradient2 = SQR(dudx[iCell], dudy[iCell], dudz[iCell]) + 
								SQR(dvdx[iCell], dvdy[iCell], dvdz[iCell]) + 
								SQR(dwdx[iCell], dwdy[iCell], dwdz[iCell]);

		RDouble nuet = oRefReNumber * turbulentViscosity / density;
		RDouble rdt  = nuet / (kap2 * wallDistance2 * MAX(sqrt(sumOfGradient2), 1.0e-10));

		RDouble fdt  = 1.0 - tanh( POWER3(8.0 * rdt) );
		RDouble fdb  = MAX((1.0 - fdt), fb);

		RDouble ct   = 1.63;
		RDouble cL   = 3.55;
		RDouble nuel = oRefReNumber * laminarViscosity / density;
		RDouble rdl  = nuel / (kap2 * wallDistance2 * MAX(sqrt(sumOfGradient2), 1.0e-10));

		RDouble fl   = tanh( pow(cL * cL * rdl, 10.0) );
		RDouble ft   = tanh( POWER3(ct * ct * rdt) );

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
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

	DataContainer *cdata = actkey->GetData();
	cdata->moveto_begin();

	int numberOfStatisticalStep = 0;
	cdata->read(&numberOfStatisticalStep, sizeof(int));

	GlobalDataBase::UpdateData("nStatisticalStep", &numberOfStatisticalStep,PHINT,1);

	for (int m = 0; m < nTurbulentEquation; ++ m)
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
	RDouble **q      = reinterpret_cast<RDouble **> (grid->GetDataPtr("q" ));
	RDouble **qTurbulence = reinterpret_cast<RDouble **> (grid->GetDataPtr("q_turb"));
	RDouble *viscousLaminar   = reinterpret_cast<RDouble * > (grid->GetDataPtr("visl"));
	RDouble *viscousTurbulence = reinterpret_cast<RDouble * > (grid->GetDataPtr("vist"));

	const RDouble ct3 = 1.2;
	const RDouble ct4 = 0.5;
	const RDouble karm   = 0.41;
	const RDouble cb1    = 0.1355;
	const RDouble cb2    = 0.622;
	const RDouble cv1    = 7.1;
	const RDouble cw2    = 0.3;
	const RDouble cw3    = 2.0;
	const RDouble sigma  = 2.0/3.0;
	const RDouble cdes   = - 1.0;
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
		RDouble xsi2 = SQR(xsi);

		RDouble turbulentFt2 = 0.0; 
		//int codeOfFt2 = 0; //PHSPACE::GetCodeOfFt2();
		//if (codeOfFt2)
		//{
		//	turbulentFt2 = ct3 * exp(- ct4 * xsi2);
		//}

		RDouble turbulentFv1 = xsi3 / (xsi3 + cv13);
		RDouble turbulentFv2 = one - xsi / (one + xsi * turbulentFv1);

		RDouble numerator   =  1.0  -  (cb1 / (cw1 * kap2 * turbulentFwStar)) * (turbulentFt2 + (1.0 - turbulentFt2) * turbulentFv2);
		RDouble denominator =  turbulentFv1 * MAX(1.0e-10, (1.0 - turbulentFt2));

		lowReynoldsNumberCorrection[iCell] = sqrt( MIN(1.0e2, numerator/denominator) );
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

    RDouble **gradientTurbulenceX = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("gradTurbulentX"));
    RDouble **gradientTurbulenceY = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("gradTurbulentY"));
    RDouble **gradientTurbulenceZ = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("gradTurbulentZ"));
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
    int nTurbulentEquation = parameters->GetNTurbulentEquation();
    if (nTurbulentEquation == 1)
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