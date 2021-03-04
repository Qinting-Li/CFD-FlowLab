#include "NSSolverUnstruct.h"
#include "MultiGridOperation.h"
#include "GradientOperation.h"
#include "Math_Limiter.h"
#include "JetBC.h"
#include "AerodynamicsBasicRelations.h"
#include "Glb_Dimension.h"
#include "Constants.h"
#include "SubsonicBC.h"
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
#include <cmath>
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
	Param_NSSolverUnstruct * parameters = GetControlParameters();
    int viscousType = parameters->GetViscousType();

    //string uns_scheme_name;
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

    if (viscousType == 0) return;

    string uns_vis_name = "default";
    GlobalDataBase::GetData("uns_vis_name", &uns_vis_name, PHSTRING, 1);

    int uns_vis;
    if ( uns_vis_name.substr(0,4) == "test" )
    {
        uns_vis = VIS_TEST;
    }
    else if ( uns_vis_name.substr(0,4) == "aver" )
    {
        uns_vis = VIS_AVER;
    }
    else if ( uns_vis_name.substr(0,3) == "std" )
    {
        uns_vis = VIS_STD;
    }
    else if ( uns_vis_name.substr(0,4) == "new1" )
    {
        uns_vis = VIS_NEW1;
    }
    else if ( uns_vis_name.substr(0,4) == "new2" )
    {
        uns_vis = VIS_NEW2;
    }
    else if ( uns_vis_name.substr(0,4) == "eric" )
    {
        uns_vis = VIS_ERIC;
    }
    else
    {
        uns_vis = VIS_STD;
    }

    GlobalDataBase::UpdateData("uns_vis",&uns_vis,PHINT,1);
}

void NSSolverUnstruct::AllocateGlobalVar(Grid *gridIn)
{
	CFDSolver::AllocateGlobalVar(gridIn);

    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

	Param_NSSolverUnstruct * parameters = GetControlParameters();
    int nNSEquation = parameters->GetNSEquationNumber();

    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();

    int nEquation = nLaminar + nChemical;

    RDouble **q   = NewPointer2<RDouble>(nEquation, nTotal);
    RDouble **res = NewPointer2<RDouble>(nEquation, nTotal);
    RDouble *spec = new RDouble [nTotal];

    grid->UpdateDataPtr("q", q);
    grid->UpdateDataPtr("res", res);
    grid->UpdateDataPtr("spec", spec);
    
    int nTotalNode = grid->GetNTotalNode();
    RDouble **qnode = NewPointer2<RDouble>(nEquation, nTotalNode);
    grid->UpdateDataPtr("qnode", qnode);

    RDouble **tnode = NewPointer2<RDouble>(1, nTotalNode);
    grid->UpdateDataPtr("tnode", tnode);

    RDouble *nodeWeight = new RDouble[nTotalNode];
    grid->UpdateDataPtr("nodeWeight", nodeWeight);

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
        grid->UpdateDataPtr("qInterPoint", qInterPoint);

        RDouble **tInterPoint = NewPointer2<RDouble>(1, numberOfInterpoints);
        grid->UpdateDataPtr("tInterPoint", tInterPoint);
    }

    RDouble *rtem = new RDouble [nTotal];
    grid->UpdateDataPtr("rtem"  , rtem);

    RDouble *srad = new RDouble [nTotalCell];
    grid->UpdateDataPtr("srad"  , srad);

    RDouble *srv = new RDouble [nTotalCell];
    grid->UpdateDataPtr("srv"  , srv);

    int ns = nEquation - nNSEquation;

    int nChemicalRadius = parameters->GetNChemicalRadius();

    if ( nChemical == 1 && nChemicalRadius == 1 )
    {
        RDouble **srs = NewPointer2<RDouble>(ns,nTotalCell);
		grid->UpdateDataPtr("srs"  , srs);
	}

    int isUnsteady = parameters->GetIsUnsteady();

	if (isUnsteady != 0)
	{
        RDouble **q_unsteady_n1   = NewPointer2<RDouble>(nEquation, nTotal);
        RDouble **q_unsteady_n2   = NewPointer2<RDouble>(nEquation, nTotal);
        RDouble **res_unsteady_n1 = NewPointer2<RDouble>(nEquation, nTotal);
        RDouble **res_unsteady_n2 = NewPointer2<RDouble>(nEquation, nTotal);

        //! Store the historical sum of InviscidFlux, ViscousFlux, et al., exclusive of DualTimeSourceFlux, used in UpdateUnsteadyFlow for the Rn in the dualtime method.
        RDouble **res_unsteady_tmp = NewPointer2<RDouble>(nEquation, nTotal);

		grid->UpdateDataPtr("q_unsteady_n1"  , q_unsteady_n1  );
		grid->UpdateDataPtr("q_unsteady_n2"  , q_unsteady_n2  );
		grid->UpdateDataPtr("res_unsteady_n1", res_unsteady_n1);
		grid->UpdateDataPtr("res_unsteady_n2", res_unsteady_n2);
        grid->UpdateDataPtr("res_unsteady_tmp", res_unsteady_tmp);

		//! Statistical variables for DES simulation.
		int ifStaticsFlowField = parameters->GetIfStaticsFlowField();
		if (ifStaticsFlowField > 0)
		{
            RDouble **qAverage = NewPointer2<RDouble>(nEquation, nTotal); 
			grid->UpdateDataPtr("qAverage", qAverage);
		}
	}

    int ntmodel = parameters->GetTemperatureModel();

    RDouble **t = NewPointer2<RDouble>(ntmodel, nTotal);
	grid->UpdateDataPtr("t"  , t);

    RDouble * dt  = new RDouble [nTotal];
	grid->UpdateDataPtr("dt" , dt );

    RDouble * gama  = new RDouble [nTotal];
	grid->UpdateDataPtr("gama", gama);

    int isViscous = parameters->IsViscous();

    if (isViscous)
    {
        RDouble *visl = new RDouble [nTotal];
        RDouble *vist = new RDouble [nTotal];

        PHSPACE::SetField(visl, 1.0, nTotal);
        PHSPACE::SetField(vist, 0.0, nTotal);

        grid->UpdateDataPtr("visl",visl);
		grid->UpdateDataPtr("vist",vist);
    }

    RDouble **gradPrimtiveVarFaceX = NewPointer2<RDouble>(nEquation, nTotal);
    RDouble **gradPrimtiveVarFaceY = NewPointer2<RDouble>(nEquation, nTotal);
    RDouble **gradPrimtiveVarFaceZ = NewPointer2<RDouble>(nEquation, nTotal);
    grid->UpdateDataPtr("gradPrimtiveVarFaceX", gradPrimtiveVarFaceX);
    grid->UpdateDataPtr("gradPrimtiveVarFaceY", gradPrimtiveVarFaceY);
    grid->UpdateDataPtr("gradPrimtiveVarFaceZ", gradPrimtiveVarFaceZ);

    RDouble **gradTemperatureFaceX = NewPointer2<RDouble>(ntmodel, nTotal);
    RDouble **gradTemperatureFaceY = NewPointer2<RDouble>(ntmodel, nTotal);
    RDouble **gradTemperatureFaceZ = NewPointer2<RDouble>(ntmodel, nTotal);
    grid->UpdateDataPtr("gradTemperatureFaceX", gradTemperatureFaceX);
    grid->UpdateDataPtr("gradTemperatureFaceY", gradTemperatureFaceY);
    grid->UpdateDataPtr("gradTemperatureFaceZ", gradTemperatureFaceZ);

    systemgridtype = GlobalDataBase::GetIntParaFromDB("sys_gridtype");
    if(MIXGRID == systemgridtype)
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
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

	Param_NSSolverUnstruct * parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();

    int nChemicalRadius = parameters->GetNChemicalRadius();

    RDouble **q    = reinterpret_cast <RDouble **> (grid->GetDataPtr("q"));
    RDouble **res  = reinterpret_cast <RDouble **> (grid->GetDataPtr("res"));
    RDouble *spec = reinterpret_cast <RDouble *> (grid->GetDataPtr("spec"));

    DelPointer2(q);
    DelPointer2(res);
    delete []spec;

    if (grid->GetLevel() == 0)
    {
        int *nodeBCType = reinterpret_cast <int *> (grid->GetDataPtr("nodeBCType"));
        delete []nodeBCType;
    }

    RDouble *rtem = reinterpret_cast <RDouble *> (grid->GetDataPtr("rtem"));
    delete []rtem;
    RDouble *srad = reinterpret_cast <RDouble *> (grid->GetDataPtr("srad"));
    delete []srad;
    RDouble *srv  = reinterpret_cast <RDouble *> (grid->GetDataPtr("srv"));
    delete []srv;

    if ( nChemical == 1 && nChemicalRadius == 1 )
    {
        RDouble **srs = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("srs") );
		DelPointer2(srs  );
	}

    int isUnsteady = parameters->GetIsUnsteady();
	if (isUnsteady != 0)
	{
        RDouble **q_unsteady_n1   = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("q_unsteady_n1"  ) );
        RDouble **q_unsteady_n2   = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("q_unsteady_n2"  ) );
        RDouble **res_unsteady_n1 = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("res_unsteady_n1") );
        RDouble **res_unsteady_n2 = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("res_unsteady_n2") );

        //! Store the historical sum of InviscidFlux, ViscousFlux, et al., exclusive of DualTimeSourceFlux, used in UpdateUnsteadyFlow for the Rn in the dualtime method.
        RDouble **res_unsteady_tmp = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("res_unsteady_tmp") );

		DelPointer2(q_unsteady_n1  );
		DelPointer2(q_unsteady_n2  );
		DelPointer2(res_unsteady_n1);
		DelPointer2(res_unsteady_n2);
		DelPointer2(res_unsteady_tmp);

		//! Statistical variables for DES simulation.
		int ifStaticsFlowField = parameters->GetIfStaticsFlowField();
		if (ifStaticsFlowField > 0)
		{
            RDouble **qAverage = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("qAverage") );
			DelPointer2(qAverage);
		}
	}

    RDouble **t   = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("t"  ) );
	DelPointer2(t  );

    RDouble * dt  = reinterpret_cast< RDouble * > ( grid->GetDataPtr("dt"  ) );
	delete [] dt;

    RDouble * gama  = reinterpret_cast< RDouble * > ( grid->GetDataPtr("gama"  ) );
	delete [] gama;

    int viscousType = parameters->GetViscousType();

    if ( viscousType > 0 )
    {
        RDouble * visl = reinterpret_cast< RDouble * > ( grid->GetDataPtr("visl") );
        RDouble * vist = reinterpret_cast< RDouble * > ( grid->GetDataPtr("vist") );
		delete [] visl;
		delete [] vist;
	}

    RDouble **gradPrimtiveVarFaceX = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarFaceX"));
    RDouble **gradPrimtiveVarFaceY = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarFaceY"));
    RDouble **gradPrimtiveVarFaceZ = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarFaceZ"));

    DelPointer2(gradPrimtiveVarFaceX);
    DelPointer2(gradPrimtiveVarFaceY);
    DelPointer2(gradPrimtiveVarFaceZ);

    RDouble **gradTemperatureFaceX = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("gradTemperatureFaceX"  ) );
    RDouble **gradTemperatureFaceY = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("gradTemperatureFaceY"  ) );
    RDouble **gradTemperatureFaceZ = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("gradTemperatureFaceZ"  ) );

    DelPointer2(gradTemperatureFaceX);
    DelPointer2(gradTemperatureFaceY);
    DelPointer2(gradTemperatureFaceZ);
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
    if ( infile )
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
        averageFlowFile = PHSPACE::AddSymbolToFileName( averageFlowFile, "_", 0 );
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
	Param_NSSolverUnstruct * parameters = GetControlParameters();
    int level = actkey->level;
	UnstructGrid * grid = UnstructGridCast(GetGrid(level));

    RDouble **q = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("q") );
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();

    int nEquation = nLaminar + nChemical;

    int outnstep = 0;
	GlobalDataBase::GetData("outnstep",&outnstep,PHINT,1);

    int nBoundFace = grid->GetNBoundFace();
    int nTotalCell = grid->GetNTotalCell();
	int nTotal = nTotalCell + nBoundFace;

    int nIFace = grid->GetNIFace();

    DataContainer * cdata = actkey->GetData();
    cdata->moveto_begin();
    cdata->write(&outnstep,sizeof(int));

	for ( int m = 0; m < nEquation; ++ m )
	{
        cdata->write(q[m],sizeof(RDouble) * nTotal);
	}

    int isUnsteady = parameters->GetIsUnsteady();

	if (isUnsteady == 0) return;

    RDouble physicalTime = 0;
    GlobalDataBase::GetData("physicalTime",&physicalTime, PHDOUBLE, 1);

    cdata->write(&physicalTime, sizeof(RDouble));

    RDouble **qn1 = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("q_unsteady_n1") );
    RDouble **qn2 = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("q_unsteady_n2") );

	for ( int m = 0; m < nEquation; ++ m )
	{
        cdata->write(qn1[m], sizeof(RDouble) * nTotal);
	}

	for ( int m = 0; m < nEquation; ++ m )
	{
        cdata->write(qn2[m],sizeof(RDouble) * nTotal);
	}

    RDouble **res   = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("res"            ) );
    RDouble **resn1 = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("res_unsteady_n1") );
    RDouble **resn2 = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("res_unsteady_n2") );

	for ( int m = 0; m < nEquation; ++ m )
	{
        cdata->write(res[m], sizeof(RDouble) * nTotal);
	}

	for ( int m = 0; m < nEquation; ++ m )
	{
        cdata->write(resn1[m], sizeof(RDouble) * nTotal);
	}

	for ( int m = 0; m < nEquation; ++ m )
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

            if (IsNeedDESStatistics())
            {
                int nStatisticalStep = 0;
                GlobalDataBase::GetData("nStatisticalStep", &nStatisticalStep, PHINT, 1);
                WriteData(actkey->filepos, &nStatisticalStep, "nStatisticalStep");
            }
        }
    }

    Grid *grid = GetGrid(actkey->level);
    int gridID = grid->GetZoneID();

    UnstructGrid *unstrGrid = UnstructGridCast(GetGrid(actkey->level));

    int nEquation = GetNumberOfEquations();

    int nBoundFace = grid->GetNBoundFace();
    int nTotalCell = grid->GetNTotalCell();
    int nTotal = nTotalCell + nBoundFace;

    RDouble **q = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("q") );
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
        RDouble **qn1 = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("q_unsteady_n1") );
        for (int iEquation = 0; iEquation < nEquation; ++ iEquation)
        {
            for (int iCell = 0; iCell < nTotal; ++ iCell)
            {
                primitiveVariables[iEquation][iCell] = qn1[iEquation][iCell];
            }
        }
        WriteData(grploc, primitiveVariables[0], "q_unsteady_n1");

        RDouble **qn2 = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("q_unsteady_n2") );
        for (int iEquation = 0; iEquation < nEquation; ++ iEquation)
        {
            for (int iCell = 0; iCell < nTotal; ++ iCell)
            {
                primitiveVariables[iEquation][iCell] = qn2[iEquation][iCell];
            }
        }
        WriteData(grploc, primitiveVariables[0], "q_unsteady_n2");

        RDouble **res   = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("res") );
        for (int iEquation = 0; iEquation < nEquation; ++ iEquation)
        {
            for (int iCell = 0; iCell < nTotal; ++ iCell)
            {
                primitiveVariables[iEquation][iCell] = res[iEquation][iCell];
            }
        }
        WriteData(grploc, primitiveVariables[0], "res");

        RDouble **resn1 = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("res_unsteady_n1") );
        for (int iEquation = 0; iEquation < nEquation; ++ iEquation)
        {
            for (int iCell = 0; iCell < nTotal; ++ iCell)
            {
                primitiveVariables[iEquation][iCell] = resn1[iEquation][iCell];
            }
        }
        WriteData(grploc, primitiveVariables[0], "res_unsteady_n1");

        RDouble **resn2 = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("res_unsteady_n2") );
        for (int iEquation = 0; iEquation < nEquation; ++ iEquation)
        {
            for (int iCell = 0; iCell < nTotal; ++ iCell)
            {
                primitiveVariables[iEquation][iCell] = resn2[iEquation][iCell];
            }
        }
        WriteData(grploc, primitiveVariables[0], "res_unsteady_n2");

        if (IsNeedDESStatistics())
        {
            RDouble **qAverage = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("qAverage") );
            for (int iEquation = 0; iEquation < nEquation; ++ iEquation)
            {
                for (int iCell = 0; iCell < nTotal; ++ iCell)
                {
                    primitiveVariables[iEquation][iCell] = qAverage[iEquation][iCell];
                }
            }
            WriteData(grploc, primitiveVariables[0], "qAverage");
        }
    }

    H5Gclose(grploc);

    DelPointer2(primitiveVariables);
}

void NSSolverUnstruct::ReadRestartH5(ActionKey *actkey)
{
    Grid *grid = GetGrid(actkey->level);
    int gridID = grid->GetZoneID();

    Param_NSSolverUnstruct * parameters = GetControlParameters();
    UnstructGrid *unstrGrid = UnstructGridCast(GetGrid(actkey->level));

    int nEquation = GetNumberOfEquations();

    int outputStepInterval = 0;
    ReadData(actkey->filepos, &outputStepInterval, "outnstep");
    GlobalDataBase::UpdateData("outnstep", &outputStepInterval, PHINT, 1);

    int nBoundFace = grid->GetNBoundFace();
    int nTotalCell = grid->GetNTotalCell();
    int nTotal = nTotalCell + nBoundFace;

    RDouble **q = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("q") );
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
    if (isUnsteady == 0)
    {
        H5Gclose(grploc);
        DelPointer2(primitiveVariables);

        return;
    }

    RDouble **qn1 = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("q_unsteady_n1") );
    RDouble **qn2 = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("q_unsteady_n2") );
    RDouble **res   = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("res"            ) );
    RDouble **resn1 = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("res_unsteady_n1") );
    RDouble **resn2 = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("res_unsteady_n2") );
    RDouble **qAverage = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("qAverage") );

    int ifStartFromSteadyResults = parameters->GetIfStartFromSteadyResults();
    RDouble physicalTime = 0.0;

    if (ifStartFromSteadyResults)
    {
        PrintToWindow("Restart from steady NS flow field, reset outer step to be zero!\n");
        //! Start from steady flow field.
        //! Reset the outer step when start from steady flow.
        outputStepInterval = 0;
        GlobalDataBase::UpdateData("outnstep",&outputStepInterval,PHINT,1);

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
    int ifStaticsFlowField = parameters->GetIfStaticsFlowField();
    int startStatisticStep = GlobalDataBase::GetIntParaFromDB("startStatisticStep");
    bool isReadAverageFlow = ifStaticsFlowField && (outputStepInterval >= startStatisticStep);
    if(isReadAverageFlow)
    {
        ReadData(grploc, primitiveVariables[0], "qAverage");
        for (int iEquation = 0; iEquation < nEquation; ++ iEquation)
        {
            for (int iCell = 0; iCell < nTotal; ++ iCell)
            {
                qAverage[iEquation][iCell] = primitiveVariables[iEquation][iCell];
            }
        }

        ReadData(actkey->filepos, &nStatisticalStep, "nStatisticalStep");
        GlobalDataBase::UpdateData("nStatisticalStep", &nStatisticalStep, PHINT, 1);
    }
    else
    {
        for (int iEquation = 0; iEquation < nEquation; ++ iEquation)
        {
            for (int iCell = 0; iCell < nTotalCell; ++ iCell)
            {
                qAverage[iEquation][iCell] = q[iEquation][iCell];
            }
        }

        GlobalDataBase::UpdateData("nStatisticalStep", &nStatisticalStep, PHINT, 1);
    }

    GlobalDataBase::UpdateData("physicalTime",&physicalTime, PHDOUBLE, 1);

    H5Gclose(grploc);
    DelPointer2(primitiveVariables);
}

void NSSolverUnstruct::ReadRestart(ActionKey *actkey)
{
    int level = actkey->level;
	UnstructGrid * grid = UnstructGridCast(GetGrid(level));

	Param_NSSolverUnstruct * parameters = GetControlParameters();
	int nLaminar  = parameters->GetLaminarNumber();
	int nChemical = parameters->GetChemicalFlag();
    int nEquation = nLaminar + nChemical;

    DataContainer * cdata = actkey->GetData();
    cdata->moveto_begin();
    int outnstep = 0;
    cdata->read(&outnstep,sizeof(int));

	GlobalDataBase::UpdateData("outnstep",&outnstep,PHINT,1);

    int nBoundFace = grid->GetNBoundFace();
    int nTotalCell = grid->GetNTotalCell();
	int nTotal = nTotalCell + nBoundFace;

    int nIFace = grid->GetNIFace();

    RDouble **q = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("q") );

    for (int m = 0; m < nEquation; ++ m)
	{
        cdata->read(q[m], sizeof(RDouble) * nTotal);
	}

    //ReadInterface(actkey);

    int isUnsteady = parameters->GetIsUnsteady();
	if (isUnsteady == 0) return;

    RDouble **qn1 = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("q_unsteady_n1") );
    RDouble **qn2 = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("q_unsteady_n2") );

    RDouble **res   = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("res"            ) );
    RDouble **resn1 = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("res_unsteady_n1") );
    RDouble **resn2 = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("res_unsteady_n2") );

	int ifStartFromSteadyResults = parameters->GetIfStartFromSteadyResults();

    RDouble physicalTime = 0.0;
	if (ifStartFromSteadyResults)
	{
		PrintToWindow("Restart from steady NS flow field, reset outer step to be zero!\n");
		//! Start from steady flow field.
		//! Reset the outer step when start from steady flow.
		outnstep = 0;
		GlobalDataBase::UpdateData("outnstep",&outnstep,PHINT,1);

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

		for ( int m = 0; m < nEquation; ++ m )
		{
            cdata->read(qn1[m], sizeof(RDouble) * nTotal);
		}

		for ( int m = 0; m < nEquation; ++ m )
		{
            cdata->read(qn2[m], sizeof(RDouble) * nTotal);
		}

		for ( int m = 0; m < nEquation; ++ m )
		{
            cdata->read(res[m],sizeof(RDouble) * nTotal);
		}

		for ( int m = 0; m < nEquation; ++ m )
		{
            cdata->read(resn1[m], sizeof(RDouble) * nTotal);
		}

		for ( int m = 0; m < nEquation; ++ m )
		{
            cdata->read(resn2[m], sizeof(RDouble) * nTotal);
		}
	}

    GlobalDataBase::UpdateData("physicalTime",&physicalTime, PHDOUBLE, 1);
}

void NSSolverUnstruct::InitFlowAsRestart()
{
	Param_NSSolverUnstruct * parameters = GetControlParameters();
	UnstructGrid *grid = UnstructGridCast( GetGrid(0) );
    int outnstep = 0;
    GlobalDataBase::UpdateData("outnstep",&outnstep, PHINT, 1);

	RDouble physicalTime = 0.0;
	GlobalDataBase::UpdateData("physicalTime",&physicalTime, PHDOUBLE, 1);

    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    RDouble **q = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("q") );

    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();

    int nEquation = nLaminar + nChemical;

	RDouble *primitiveVarFarfield = parameters->GetPrimitiveVarFarfield();

	for ( int m = 0; m < nEquation; ++ m )
	{
		for ( int iCell = 0; iCell < nTotal; ++ iCell )
		{
            q[m][iCell] = primitiveVarFarfield[m];
        }
    }

	int isUnsteady = parameters->GetIsUnsteady();

	if (isUnsteady == 0)
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

	RDouble **qn1 = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("q_unsteady_n1") );
	RDouble **qn2 = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("q_unsteady_n2") );

	for ( int m = 0; m < nEquation; ++ m )
	{
		for ( int iCell = 0; iCell < nTotal; ++ iCell )
		{
			qn1[m][iCell] = q[m][iCell];
		}
	}

	for ( int m = 0; m < nEquation; ++ m )
	{
		for ( int iCell = 0; iCell < nTotal; ++ iCell )
		{
			qn2[m][iCell] = q[m][iCell];
		}
	}

	RDouble **res   = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("res"            ) );
	RDouble **resn1 = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("res_unsteady_n1") );
	RDouble **resn2 = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("res_unsteady_n2") );

	//! Store the historical sum of InviscidFlux, ViscousFlux, et al., exclusive of DualTimeSourceFlux, used in UpdateUnsteadyFlow for the Rn in the dualtime method.
	RDouble **resTmp = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("res_unsteady_tmp") );

	for ( int m = 0; m < nEquation; ++ m )
	{
		for ( int iCell = 0; iCell < nTotal; ++ iCell )
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
    UnstructGrid *grid = reinterpret_cast<UnstructGrid *>( GetGrid(level));
    InterfaceInfo *iinfo = grid->GetInterfaceInfo();
    if (iinfo == 0)
    {
        return;
    }

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();
    int nEquation = nLaminar + nChemical + nTemperatureModel - 1;

    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    PHSPACE::UploadInterfaceValue(grid, q, "q",  nEquation);

    RDouble **t = reinterpret_cast< RDouble ** > (grid->GetDataPtr("t"));
    PHSPACE::UploadInterfaceValue(grid, t, "t", nTemperatureModel);

    //RDouble** gradPrimtiveVarFaceX = reinterpret_cast<RDouble**> (grid->GetDataPtr("gradPrimtiveVarFaceX"));
    //PHSPACE::UploadInterfaceValue(grid, gradPrimtiveVarFaceX, "dqdx", nEquation);
    //
    //RDouble** gradPrimtiveVarFaceY = reinterpret_cast<RDouble**> (grid->GetDataPtr("gradPrimtiveVarFaceY"));
    //PHSPACE::UploadInterfaceValue(grid, gradPrimtiveVarFaceY, "dqdy", nEquation);
    //
    //RDouble** gradPrimtiveVarFaceZ = reinterpret_cast<RDouble**> (grid->GetDataPtr("gradPrimtiveVarFaceZ"));
    //PHSPACE::UploadInterfaceValue(grid, gradPrimtiveVarFaceZ, "dqdz", nEquation);
    //
    //RDouble** gradTemperatureFaceX = reinterpret_cast<RDouble**> (grid->GetDataPtr("gradTemperatureFaceX"));
    //PHSPACE::UploadInterfaceValue(grid, gradTemperatureFaceX, "dtdx", nTemperatureModel);
    //
    //RDouble** gradTemperatureFaceY = reinterpret_cast<RDouble**> (grid->GetDataPtr("gradTemperatureFaceY"));
    //PHSPACE::UploadInterfaceValue(grid, gradTemperatureFaceY, "dtdy", nTemperatureModel);
    //
    //RDouble** gradTemperatureFaceZ = reinterpret_cast<RDouble**> (grid->GetDataPtr("gradTemperatureFaceZ"));
    //PHSPACE::UploadInterfaceValue(grid, gradTemperatureFaceZ, "dtdz", nTemperatureModel);
    //RDouble *viscousLaminar = reinterpret_cast <RDouble *> (grid->GetDataPtr("visl"));
    //PHSPACE::UploadInterfaceValue(grid, viscousLaminar, "visl");
}

void NSSolverUnstruct::DownloadInterfaceData(ActionKey *actkey)
{
    int level = actkey->level;
    UnstructGrid * grid = reinterpret_cast< UnstructGrid * >( GetGrid(level) );
	InterfaceInfo * iinfo = grid->GetInterfaceInfo();
	if ( iinfo == 0 ) return;

	Param_NSSolverUnstruct * parameters = GetControlParameters();
	int nLaminar  = parameters->GetLaminarNumber();
	int nChemical = parameters->GetChemicalFlag();
	int nEquation = nLaminar + nChemical;

	RDouble **q = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("q"  ) );
    PHSPACE::DownloadInterfaceValue(grid, q, "q", nEquation);

	int nTemperatureModel = parameters->GetTemperatureModel();

    RDouble **t = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("t"  ) );
    PHSPACE::DownloadInterfaceValue(grid, t, "t", nTemperatureModel);

    //RDouble** gradPrimtiveVarFaceX = reinterpret_cast<RDouble**> (grid->GetDataPtr("gradPrimtiveVarFaceX"));
    //PHSPACE::DownloadInterfaceValue(grid, gradPrimtiveVarFaceX, "dqdx", nEquation);
    //
    //RDouble** gradPrimtiveVarFaceY = reinterpret_cast<RDouble**> (grid->GetDataPtr("gradPrimtiveVarFaceY"));
    //PHSPACE::DownloadInterfaceValue(grid, gradPrimtiveVarFaceY, "dqdy", nEquation);
    //
    //RDouble** gradPrimtiveVarFaceZ = reinterpret_cast<RDouble**> (grid->GetDataPtr("gradPrimtiveVarFaceZ"));
    //PHSPACE::DownloadInterfaceValue(grid, gradPrimtiveVarFaceZ, "dqdz", nEquation);
    //
    //RDouble** gradTemperatureFaceX = reinterpret_cast<RDouble**> (grid->GetDataPtr("gradTemperatureFaceX"));
    //PHSPACE::DownloadInterfaceValue(grid, gradTemperatureFaceX, "dtdx", nTemperatureModel);
    //
    //RDouble** gradTemperatureFaceY = reinterpret_cast<RDouble**> (grid->GetDataPtr("gradTemperatureFaceY"));
    //PHSPACE::DownloadInterfaceValue(grid, gradTemperatureFaceY, "dtdy", nTemperatureModel);
    //
    //RDouble** gradTemperatureFaceZ = reinterpret_cast<RDouble**> (grid->GetDataPtr("gradTemperatureFaceZ"));
    //PHSPACE::DownloadInterfaceValue(grid, gradTemperatureFaceZ, "dtdz", nTemperatureModel);
    //RDouble *viscousLaminar = reinterpret_cast <RDouble *> (grid->GetDataPtr("visl"));
    //PHSPACE::DownloadInterfaceValue(grid, viscousLaminar, "visl");
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
    Param_NSSolverUnstruct * parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nEquation = nLaminar + nChemical;

    RDouble **qnode = reinterpret_cast<RDouble **>(grid->GetDataPtr("qnode"));
    PHSPACE::UploadInterpointValue(grid, qnode, "qnode",  nEquation);

    int ntmodel = parameters->GetTemperatureModel();

    RDouble **tnode = reinterpret_cast<RDouble **>(grid->GetDataPtr("tnode"));
    PHSPACE::UploadInterpointValue(grid, tnode, "tnode", ntmodel);
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
    Param_NSSolverUnstruct * parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nEquation = nLaminar + nChemical;

    RDouble **qnode = reinterpret_cast<RDouble **>(grid->GetDataPtr("qnode"));
    PHSPACE::DownloadInterpointValue(grid, qnode, "qnode", nEquation);

    int ntmodel = parameters->GetTemperatureModel();

    RDouble **tnode = reinterpret_cast<RDouble **> (grid->GetDataPtr("tnode"));
    PHSPACE::DownloadInterpointValue(grid, tnode, "tnode", ntmodel);
}

void NSSolverUnstruct::CommunicationInterpointWeight(ActionKey *actkey)
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
                RDouble *nodeWeight = reinterpret_cast<RDouble *>(grid->GetDataPtr("nodeWeight"));

                cdata->moveto_begin();

                int iNeighbor = interpointInformation->FindIthNeighbor(neighborZoneIndex);
                int numberOfInterpointsForNeighbor = interpointInformation->GetNumberOfInterpointsForNeighbor(iNeighbor);
                int *pointIndexForSend = interpointInformation->GetPointIndexForSend(iNeighbor);
                int *interPoint2GlobalPoint = interpointInformation->GetInterPoint2GlobalPoint();

                int iPoint;
                int globalPoint;
                RDouble valueWeight;
                for (int ipointLocal = 0; ipointLocal < numberOfInterpointsForNeighbor; ++ ipointLocal )
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

                int iNeighbor = interpointInformationNeighbor->FindIthNeighbor(izone);
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
    for (int izone = 0; izone < nZones; ++ izone)
    {
        int currentZoneProc = GetZoneProcessorID(izone);

        if (myid != currentZoneProc)
        {
            continue;
        }

        UnstructGrid *grid = UnstructGridCast(PHSPACE::GetGrid(izone, level));
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
    UnstructGrid * grid = reinterpret_cast< UnstructGrid * >( GetGrid(level) );

	Param_NSSolverUnstruct * parameters = GetControlParameters();
	int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();

    int nEquation = nLaminar + nChemical;

	RDouble **q = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("q"  ) );
    PHSPACE::UploadOversetValue(grid, q, "q",  nEquation);

	int ntmodel = parameters->GetTemperatureModel();

	RDouble **t = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("t"  ) );
    PHSPACE::UploadOversetValue(grid, t, "t",  ntmodel);
}

void NSSolverUnstruct::DownloadOversetData(ActionKey *actkey)
{
	Param_NSSolverUnstruct * parameters = GetControlParameters();
    int level = actkey->level;
    UnstructGrid * grid = reinterpret_cast< UnstructGrid * >( GetGrid(level) );

	int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();

    int nEquation = nLaminar + nChemical;

	RDouble **q = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("q"  ) );
    PHSPACE::DownloadOversetValue(grid, q, "q", nEquation);

    int ntmodel = parameters->GetTemperatureModel();

	RDouble **t = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("t"  ) );
    PHSPACE::DownloadOversetValue(grid, t, "t",  ntmodel);
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
    int nEquation = nLaminar + nChemical;

	RDouble **field = NewPointer2<RDouble>(nEquation, nTotal);

	FieldProxy *fieldProxy = new FieldProxy();

	fieldProxy->SetField_UNS(field, true);

	return fieldProxy;
}

FieldProxy * NSSolverUnstruct::GetFieldProxy(Grid *gridIn, const string &field_name)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

	RDouble ** field = reinterpret_cast< RDouble ** > ( grid->GetDataPtr(field_name) );

	FieldProxy * field_proxy = new FieldProxy();

	field_proxy->SetField_UNS(field);

	return field_proxy;
}

FieldProxy * NSSolverUnstruct::GetResidualProxy(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

	RDouble **res = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("res") );

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
    int nEquation = nLaminar + nChemical;

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
    int nEquation = nLaminar + nChemical;

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

    int nEquation = nLaminar + nChemical;

	RDouble **res = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("res") );

	RDouble **rhs = rhsProxy->GetField_UNS();
	for ( int m = 0; m < nEquation; ++ m )
	{
		for ( int iCell = 0; iCell < nTotal; ++ iCell )
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
    int nBoundFace = grid->GetNBoundFace();
    int nTotalCell = grid->GetNTotalCell();
    int nTotal = nTotalCell + nBoundFace;

	Param_NSSolverUnstruct * parameters = GetControlParameters();
	int viscousType = parameters->GetViscousType();

    if ( viscousType >= 4 )
    {
	    int nrokplus = 1;
	    GlobalDataBase::GetData("nrokplus",&nrokplus,PHINT,1);

        if ( nrokplus > 0 )
        {
            RDouble **res  = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("res") );
            RDouble *sengy = reinterpret_cast< RDouble * > ( grid->GetDataPtr("sengy"   ) );
            for ( int iCell = 0; iCell < nTotalCell; ++ iCell )
            {
                res[4][iCell] += sengy[iCell];
            }
        }
    }
}

//! Load flow variables stored in grid to q.
void NSSolverUnstruct::LoadQ(Grid *gridIn, FieldProxy *q_proxy)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

	RDouble **q = q_proxy->GetField_UNS();

    RDouble **qold = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("q") );

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar = parameters->GetLaminarNumber();

    for ( int m = 0; m < nLaminar; ++ m )
    {
        for ( int iCell = 0; iCell < nTotal; ++ iCell )
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

    RDouble **cres = reinterpret_cast< RDouble ** > ( cgrid->GetDataPtr("res") );
    RDouble **fres = reinterpret_cast< RDouble ** > ( fgrid->GetDataPtr("res") );

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar = parameters->GetLaminarNumber();

    for ( int m = 0; m < nLaminar; ++ m )
    {
        for ( int iCell = 0; iCell < fnTotalCell; ++ iCell )
        {
            cres[m][ cell2coarsegridcell[iCell] ] -= fres[m][iCell];
        }
    }
}

//! Put correction back on the coarse grid.
void NSSolverUnstruct::PutCorrectionBack(Grid *gridIn, FieldProxy *q_proxy)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

	RDouble **q = q_proxy->GetField_UNS();

    RDouble **qold = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("q") );

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar = parameters->GetLaminarNumber();

    for ( int m = 0; m < nLaminar; ++ m )
    {
        for ( int iCell = 0; iCell < nTotalCell; ++ iCell )
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

	int cnTCell = cgrid->GetNTotalCell();
	int cnBFace = cgrid->GetNBoundFace();
	int cnTotal = cnTCell + cnBFace;

	RDouble **fq = reinterpret_cast< RDouble ** > ( fgrid->GetDataPtr("q") );
	RDouble **cq = reinterpret_cast< RDouble ** > ( cgrid->GetDataPtr("q") );

    Param_NSSolverUnstruct *parameters = GetControlParameters();
	int nLaminar  = parameters->GetLaminarNumber();
	int nChemical = parameters->GetChemicalFlag();
	int nEquation = nLaminar + nChemical;

	int fnBFace = fgrid->GetNBoundFace();
	int fnTotal = fnBFace + fnTCell;	

	RDouble mgCorrectionLimit = parameters->GetMgCorrectionLimit();

	for ( int m = 0; m < nLaminar; ++ m )
	{
		for ( int iCell = 0; iCell < fnTCell; ++ iCell )
		{
			int coarseCell = cell2coarsegridcell[iCell];

			RDouble fgridValue = fq[m][iCell];
			RDouble correction = cq[m][coarseCell];
			RDouble DQ         = MIN(fabs(correction), fabs(fgridValue) * mgCorrectionLimit);
			fq[m][iCell]      += SIGN(DQ, correction);
		}
	}

	//FreeGradientProxy(cgrid);
}

void NSSolverUnstruct::CorrectFineGridHighOrder(Grid *fgrid_in, Grid *cgrid_in)
{
	UnstructGrid *fgrid = UnstructGridCast(fgrid_in);
	UnstructGrid *cgrid = UnstructGridCast(cgrid_in);

	int fnTNode = fgrid->GetNTotalNode();
	int fnTCell = fgrid->GetNTotalCell();
	int fnBFace = fgrid->GetNBoundFace();
	int fnTFace = fgrid->GetNTotalFace();
	int fnTotal = fnBFace + fnTCell;

	int * cell2CoarseGridCell = fgrid->GetCell2CoarseGridCell();
	int * fleft_cell_of_face  = fgrid->GetLeftCellOfFace();
	int * fright_cell_of_face = fgrid->GetRightCellOfFace();

	int ** face2nodeArray = fgrid->GetFace2NodeArray();
	int * fnode_number_of_each_face = fgrid->GetNodeNumberOfEachFace();

	RDouble **fq = reinterpret_cast< RDouble ** > ( fgrid->GetDataPtr("q") );
	RDouble **cq = reinterpret_cast< RDouble ** > ( cgrid->GetDataPtr("q") );

    Param_NSSolverUnstruct *parameters = GetControlParameters();
	int nLaminar  = parameters->GetLaminarNumber();
	int nChemical = parameters->GetChemicalFlag();
	int nEquation = nLaminar + nChemical;

	RDouble minRatio[5], maxRatio[5];
	for (int m = 0; m < 5; ++ m)
	{
		minRatio[m] = LARGE;
		maxRatio[m] = -LARGE;
	}

	RDouble mgCorrectionLimit = parameters->GetMgCorrectionLimit();

	for ( int m = 0; m < nLaminar; ++ m )
	{
		//if( m == 1 || m == 2 || m == 3) continue;
		//! On coarse mesh, average the cell data onto nodes.
		RDouble * fqNode  = new RDouble [ fnTNode ];
		int * fnodeCount = new int [fnTNode];
		for (int iNode = 0; iNode < fnTNode; ++ iNode )
		{
			fqNode[iNode]     = 0.0;
			fnodeCount[iNode] = 0;
		}

		for (int iFace = 0; iFace < fnBFace; ++ iFace)
		{
			int le = fleft_cell_of_face[iFace];
			int leftCoarseCell = cell2CoarseGridCell[le];

			int nNode = fnode_number_of_each_face[iFace];
			for (int iNode = 0; iNode < nNode; ++ iNode )
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
			for (int iNode = 0; iNode < nNode; ++ iNode )
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
			for (int iNode = 0; iNode < nNode; ++ iNode )
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
			for (int iNode = 0; iNode < nNode; ++ iNode )
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
		for ( int iCell = 0; iCell < fnTCell; ++ iCell )
		{		
			RDouble fgridValue = fq[m][iCell];
			RDouble correction = fcellData[iCell] / fcellCount[iCell];
			RDouble DQ         = MIN(fabs(correction), fabs(fgridValue) * mgCorrectionLimit);
			fq[m][iCell]      += SIGN(DQ, correction);

			if(m == IDX::IP || m == IDX::IR)
			{
				if ( fq[m][iCell] <= 1.0e-10 )
				{
					cout << " Negative found after correction, variable = " << m 
						 << ", find grid Value = " << fgridValue	<< ", correction = " << correction						
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
    int nEquation = nLaminar + nChemical;

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
    if (viscousType > 0)
    {
        RDouble *fineViscousLaminar   = reinterpret_cast<RDouble *> (fineGrid->GetDataPtr("visl"));
        RDouble *coarseViscousLaminar = reinterpret_cast<RDouble *> (coarseGrid->GetDataPtr("visl"));

        InterpolatQ(fineGrid, fineViscousLaminar, coarseGrid, coarseViscousLaminar);

        if (viscousType > 1 && isRunTurbOnCoarseGrid)
        {
            RDouble *fineViscousTurbulent   = reinterpret_cast<RDouble *> (fineGrid->GetDataPtr("vist"));
            RDouble *coarseViscousTurbulent = reinterpret_cast<RDouble *> (coarseGrid->GetDataPtr("vist"));

            InterpolatQ(fineGrid, fineViscousTurbulent, coarseGrid, coarseViscousTurbulent);
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
    int nEquation = nLaminar + nChemical;

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
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int isUnsteady = parameters->GetIsUnsteady();

    if (isUnsteady == 0)
    {
        return zero;
    }

    RDouble **q   = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("q") );
    RDouble **qn1 = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("q_unsteady_n1") );
    RDouble **qn2 = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("q_unsteady_n2") );
    RDouble **res = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("res") );

    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();

    int ntmodel = parameters->GetTemperatureModel();

    int nEquation = nLaminar + nChemical + ntmodel - 1;

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

        if (ntmodel > 1)
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

    RDouble cvg = sqrt( ABS( sum1 / ( sum2 + SMALL ) ) );
    return cvg;
}

void NSSolverUnstruct::UpdateUnsteadyFlow(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

	Param_NSSolverUnstruct * parameters = GetControlParameters();
    int isUnsteady = parameters->GetIsUnsteady();

	if (isUnsteady == 0) return;

	RDouble *vol   = grid->GetCellVolume();
	RDouble *voln1 = grid->GetCellVolume();
	RDouble *voln2 = grid->GetCellVolume();

	RDouble **q   = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("q") );
	RDouble **qn1 = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("q_unsteady_n1") );
	RDouble **qn2 = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("q_unsteady_n2") );

	RDouble **resn1 = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("res_unsteady_n1") );
	RDouble **resn2 = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("res_unsteady_n2") );
	RDouble **resTmp = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("res_unsteady_tmp") );

    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();

	int nEquation = nLaminar + nChemical;

    for ( int m = 0; m < nEquation; ++ m )
    {
        for ( int iCell = 0; iCell < nTotal; ++ iCell )
        {
            qn2[ m ][ iCell ] = qn1[ m ][ iCell ];
            qn1[ m ][ iCell ] = q  [ m ][ iCell ];

            resn2[ m ][ iCell ] = resn1 [ m ][ iCell] ;
            resn1[ m ][ iCell ] = resTmp[ m ][ iCell] ; //! Here the current outnstep is over, the value of the stored resTmp should be assigned to resn1 for the next outnstep. It should be noticed that resTmp only contain the inviscid and viscous flux.
        }
    }

	//! Statistical variables for DES simulation.
	int ifStaticsFlowField = parameters->GetIfStaticsFlowField();
	if (ifStaticsFlowField != 0)
	{
		int outnstep = GlobalDataBase::GetIntParaFromDB("outnstep");	
		int startStatisticStep = GlobalDataBase::GetIntParaFromDB("startStatisticStep");

		if( outnstep >= startStatisticStep )
		{
			int nStatisticalStep   = GlobalDataBase::GetIntParaFromDB("nStatisticalStep");

			RDouble c1 = (nStatisticalStep - 1.0) / (nStatisticalStep * 1.0);
			RDouble c2 = 1.0 / (nStatisticalStep * 1.0);

			RDouble **qAverage = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("qAverage") );
			for ( int m = 0; m < nEquation; ++ m )
			{
				for ( int iCell = 0; iCell < nTotal; ++ iCell )
				{
					qAverage[m][iCell] = c1 * qAverage[m][iCell] + c2 * q[m][iCell];
				}
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
    if (isUnsteady == 0)
    {
        return;
    }

    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

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
    int ntmodel = parameters->GetTemperatureModel();
    int nEquation = nLaminar + nChemical + ntmodel - 1;

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
    const int methodOfDualTime = 3;
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
        if (ntmodel > 1)
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
}


void NSSolverUnstruct::ComputeGamaAndTemperature(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    Param_NSSolverUnstruct *parameters = GetControlParameters();

    RDouble **t   = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("t") );
    RDouble **q   = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    RDouble *gama = reinterpret_cast< RDouble *  > (grid->GetDataPtr("gama"));

    RDouble refGama = parameters->GetRefGama();

    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int ntmodel = parameters->GetTemperatureModel();
    int nEquation = nLaminar + nChemical + ntmodel - 1;

    using namespace GAS_SPACE;
    using namespace IDX;

    if (nChemical == 1)
    {
        RDouble Tv = 0.0, Te = 0.0;
        RDouble *prim = new RDouble[nEquation];
        for (int iCell = 0; iCell < nTotal; ++ iCell)
        {
            for (int m = 0; m < nEquation; ++ m)
            {
                prim[m] = q[m][iCell];
            }
            gas->GetSpecificHeatRatioAndTemperatute(prim, gama[iCell], t[ITT][iCell], Tv, Te);
            if (ntmodel > 1)
            {
                t[ITV][iCell] = Tv;
                if (ntmodel == 3)
                {
                    t[ITE][iCell] = Te;
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
            t[ITT][iCell] = pm / ( coefficientofStateEquation * rm);
        }
    }
}

// void NSSolverUnstruct::ComputeNodeValue(Grid *gridIn)
// {
//     UnstructGrid *grid = UnstructGridCast(gridIn);
// 
//     int nTotalCell = grid->GetNTotalCell();
//     int nBoundFace = grid->GetNBoundFace();
//     int nTotal = nTotalCell + nBoundFace;
//     int nTotalNode = grid->GetNTotalNode();
//     int nTotalFace = grid->GetNTotalFace();
// 
//     RDouble **q = reinterpret_cast<RDouble **>(grid->GetDataPtr("q"));
//     RDouble **t = reinterpret_cast<RDouble **>(grid->GetDataPtr("t"));
// 
//     RDouble     **qNode = reinterpret_cast<RDouble **>(grid->GetDataPtr("qnode"));
//     RDouble     **tNode = reinterpret_cast<RDouble **>(grid->GetDataPtr("tnode"));
//     RDouble *nodeWeight = reinterpret_cast<RDouble *>(grid->GetDataPtr("nodeWeight"));
// 
//     Param_NSSolverUnstruct *parameters = GetControlParameters();
//     int nLaminar  = parameters->GetLaminarNumber();
//     int nChemical = parameters->GetChemicalFlag();
//     int nEquation = nLaminar + nChemical;
// 
//     int *leftCellOfFace = grid->GetLeftCellOfFace();
//     int *rightCellOfFace = grid->GetRightCellOfFace();
//     int *face2node = grid->GetFace2Node();
//     int *nodeNumberOfEachFace = grid->GetNodeNumberOfEachFace();
// 
//     BCRecord **bcRecord = grid->GetBCRecord();
// 
//     for (int iNode = 0; iNode < nTotalNode; ++ iNode)
//     {
//         nodeWeight[iNode] = 0.0;
//         for(int m =0; m < nEquation; ++ m)
//         {
//             qNode[m][iNode] = 0; 
//         }
//         tNode[0][iNode] = 0;
//     }
// 
//     //! Boundary faces.
//     int nodePosition = 0;
//     for (int iFace = 0; iFace < nBoundFace; ++ iFace)
//     {
//         int bcType = bcRecord[iFace]->GetKey();
//         int le = leftCellOfFace [iFace];
//         int re = rightCellOfFace[iFace];
//         if (bcType != PHENGLEI::INTERFACE)
//         {
//             for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
//             {
//                 //! From left.
//                 int point = face2node[nodePosition + jNode];
// 
//                 for (int m = 0; m < nEquation; ++ m)
//                 {
//                     qNode[m][point] += q[m][le];
//                 }            
//                 tNode[0][point] += t[0][le];
//                 nodeWeight[point] += 1.0;
// 
//                 //! From right.
//                 for (int m = 0; m < nEquation; ++ m)
//                 {
//                     qNode[m][point] += q[m][re];
//                 }
//                 tNode[0][point] += t[0][re];
//                 nodeWeight[point] += 1.0;
//             }
//         }
//         else
//         {
//             for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
//             {
//                 int point = face2node[nodePosition + jNode];
//                 for (int m = 0; m < nEquation; ++ m)
//                 {
//                     qNode[m][point] += q[m][le];
//                 }
//                 tNode[0][point] += t[0][le];
//                 nodeWeight[point] += 1.0;
//             }
//         }        
//         nodePosition += nodeNumberOfEachFace[iFace];
//     }
// 
//     //! Interior faces.
//     for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
//     {
//         int le = leftCellOfFace [iFace];
//         int re = rightCellOfFace[iFace];
//         for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
//         {
//             //! From left.
//             int point = face2node[nodePosition + jNode];
// 
//             for (int m = 0; m < nEquation; m++)
//             {
//                 qNode[m][point] += q[m][le];
//             }            
//             tNode[0][point] += t[0][le];
//             nodeWeight[point] += 1.0;
// 
//             for (int m = 0; m < nEquation; ++ m)
//             {
//                 qNode[m][point] += q[m][re];
//             }
//             tNode[0][point] += t[0][re];
//             nodeWeight[point] += 1.0;
//         }
//         nodePosition += nodeNumberOfEachFace[iFace];
//     }
// 
//     nodePosition = 0;
//     for (int iFace = 0; iFace < nBoundFace; ++ iFace)
//     {
//         int bcType = bcRecord[iFace]->GetKey();
//         if (bcType == PHENGLEI::SYMMETRY)
//         {
//             for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
//             {
//                 int point = face2node[nodePosition + jNode];
// 
//                 for (int m = 0; m < nEquation; ++ m)
//                 {
//                     qNode[m][point] = 0;
//                 }            
//                 tNode[0][point] = 0;
//                 nodeWeight[point] = 0.0;
//             }
//         }
//         nodePosition += nodeNumberOfEachFace[iFace];
//     }
// 
//     nodePosition = 0;
//     for (int iFace = 0; iFace < nBoundFace; ++ iFace)
//     {
//         int le = leftCellOfFace[iFace];
//         int re = iFace + nTotalCell;
//         int bcType = bcRecord[iFace]->GetKey();
//         if (bcType == PHENGLEI::SYMMETRY)
//         {
//             for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
//             {
//                 int point = face2node[nodePosition + jNode];    
//                 for (int m = 0; m < nEquation; ++ m)
//                 {
//                     qNode[m][point] += q[m][le];
//                 }            
//                 tNode[0][point] += t[0][le];
//                 nodeWeight[point] += 1.0;
// 
//                 for (int m = 0; m < nEquation; ++ m)
//                 {
//                     qNode[m][point] += q[m][re];
//                 }
//                 tNode[0][point] += t[0][re];
//                 nodeWeight[point] += 1.0;
//             }
//         }
//         nodePosition += nodeNumberOfEachFace[iFace];
//     }
// 
//     nodePosition = 0;
//     for (int iFace = 0; iFace < nBoundFace; ++ iFace)
//     {
//         int bcType = bcRecord[iFace]->GetKey();
//         if (bcType == PHENGLEI::FARFIELD)
//         {
//             for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
//             {
//                 int point = face2node[nodePosition + jNode];
// 
//                 for (int m = 0; m < nEquation; ++ m)
//                 {
//                     qNode[m][point] = 0;
//                 }            
//                 tNode[0][point] = 0;
//                 nodeWeight[point] = 0.0;
//             }
//         }
//         nodePosition += nodeNumberOfEachFace[iFace];
//     }
// 
//     nodePosition = 0;
//     for (int iFace = 0; iFace < nBoundFace; ++ iFace)
//     {
//         int le = leftCellOfFace[iFace];
//         int re = iFace + nTotalCell;
//         int bcType = bcRecord[iFace]->GetKey();
//         if (bcType == PHENGLEI::FARFIELD )
//         {
//             for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
//             {
//                 int point = face2node[nodePosition + jNode];
//                 for (int m = 0; m < nEquation; ++ m)
//                 {
//                     qNode[m][point] += q[m][le];
//                 }
//                 tNode[0][point] += t[0][le];
//                 nodeWeight[point] += 1.0;
// 
//                 for (int m = 0; m < nEquation; ++ m)
//                 {
//                     qNode[m][point] += q[m][re];
//                 }
//                 tNode[0][point] += t[0][re];
//                 nodeWeight[point] += 1.0;
//             }
//         }
//         nodePosition += nodeNumberOfEachFace[iFace];
//     }
// 
//     nodePosition = 0;
//     for (int iFace = 0; iFace < nBoundFace; ++ iFace)
//     {
//         int bcType = bcRecord[iFace]->GetKey();
//         if (bcType == PHENGLEI::SOLID_SURFACE)
//         {
//             for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
//             {
//                 int point = face2node[nodePosition + jNode];
//                 for (int m = 0; m < nEquation; m++)
//                 {
//                     qNode[m][point] = 0;
//                 }
//                 tNode[0][point] = 0;
//                 nodeWeight[point] = 0.0;
//             }
//         }
//         nodePosition += nodeNumberOfEachFace[iFace];
//     }
// 
//     nodePosition = 0;
//     for (int iFace = 0; iFace < nBoundFace; ++ iFace)
//     {
//         int le = leftCellOfFace[iFace];
//         int re = iFace + nTotalCell;
//         int bcType = bcRecord[iFace]->GetKey();
//         if (bcType == PHENGLEI::SOLID_SURFACE)
//         {
//             for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
//             {
//                 int point = face2node[nodePosition + jNode];
//                 for (int m = 0; m < nEquation; m++)
//                 {
//                     qNode[m][point] += q[m][le];
//                 }            
//                 tNode[0][point] += t[0][le];
//                 nodeWeight[point] += 1.0;
// 
//                 for (int m = 0; m < nEquation; m++)
//                 {
//                     qNode[m][point] += q[m][re];
//                 }
//                 tNode[0][point] += t[0][re];
//                 nodeWeight[point] += 1.0;
//             }
//         }
//         nodePosition += nodeNumberOfEachFace[iFace];
//     }
// 
//     //! Change the number of points divided by the interpoint to the number of points before partition,because it is a face loop,each point is calculated twice,so the number of points needs to be multiplied by 2.
//     InterpointInformation *interPointInfor = grid->GetInterpointInfo();
//     if (interPointInfor)
//     {
//         int numberOfInterpoints = interPointInfor->GetNumberOfInterpoints();
//         int *interPoint2GlobalPoint = interPointInfor->GetInterPoint2GlobalPoint();
//         int *cellNumberOfInterPoint = interPointInfor->GetCellNumberOfInterPoint();
//         int *labelOfInterPoint = interPointInfor->GetLabelOfInterPoint();
//         int *totalZonesOfInterPoint = interPointInfor->GetTotalZonesOfInterPoint();
// 
//         int globalPoint;
//         /*if (grid->GetLevel() == 0)
//         {
//             for (int iPoint = 0; iPoint < numberOfInterpoints; ++ iPoint)
//             {
//                 globalPoint = interPoint2GlobalPoint[iPoint];
//                 nCount[globalPoint] = cellNumberOfInterPoint[iPoint];
//             }
//         }
//         else
//         {
//             for (int iPoint = 0; iPoint < numberOfInterpoints - 1; ++ iPoint)
//             {
//                 globalPoint = interPoint2GlobalPoint[iPoint];
//                 if (globalPoint == interPoint2GlobalPoint[iPoint + 1])
//                 {
//                     continue;
//                 }
//                 nCount[globalPoint] *= totalZonesOfInterPoint[iPoint];
//             }
//             globalPoint = interPoint2GlobalPoint[numberOfInterpoints - 1];
//             nCount[globalPoint] *= totalZonesOfInterPoint[numberOfInterpoints - 1];
//         }*/
// 
//         for (int iPoint = 0; iPoint < numberOfInterpoints; ++ iPoint)
//         {
//             globalPoint = interPoint2GlobalPoint[iPoint];
//             if (labelOfInterPoint[iPoint] != 1)
//             {
//                 for (int m = 0; m < nEquation; ++ m)
//                 {
//                     qNode[m][globalPoint] = 0;
//                 }
//                 tNode[0][globalPoint] = 0;
//             }
//         }
// 
//         RDouble **qInterPoint = reinterpret_cast<RDouble **>(grid->GetDataPtr("qInterPoint"));
//         RDouble **tInterPoint = reinterpret_cast<RDouble **>(grid->GetDataPtr("tInterPoint"));
//         for (int iPoint = 0; iPoint < numberOfInterpoints; ++ iPoint)
//         {
//             for (int m = 0; m < nEquation; m++)
//             {
//                 qInterPoint[m][iPoint] =0;
//             }
//             tInterPoint[0][iPoint] =0;
//         }
//     }
// 
//     /*for (int iNode = 0; iNode < nTotalNode; ++ iNode)
//     {
//         if (nodeWeight[iNode] != 0)
//         {
//             for (int m = 0; m < nEquation; m++)
//             {
//                 qNode[m][iNode] /= nodeWeight[iNode];
//             }
//             tNode[0][iNode] /= nodeWeight[iNode];
//         }
//         else
//         {
//             for (int m = 0; m < nEquation; m++)
//             {
//                 qNode[m][iNode] = 0.0;
//             }
//             tNode[0][iNode] = 0.0;
//         }
//     }*/
// }

void NSSolverUnstruct::ComputeNodeValue(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;
    int nTotalNode = grid->GetNTotalNode();
    int nTotalFace = grid->GetNTotalFace();

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
    int nEquation = nLaminar + nChemical;

    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();
    int *nodeNumberOfEachFace = grid->GetNodeNumberOfEachFace();
    int **face2nodeArray      = grid->GetFace2NodeArray();

    BCRecord **bcRecord = grid->GetBCRecord();

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

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int bcType = bcRecord[iFace]->GetKey();
        if(bcType == PHENGLEI::SOLID_SURFACE)
        {
            //! Solid surface has the highest priority.
            for (int jNode = 0; jNode < nodeNumberOfEachFace[iFace]; ++ jNode)
            {            
                int point = face2nodeArray[iFace][jNode];
                nodeBC[point] = bcType;
            } 
        }       
    }

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int bcType = bcRecord[iFace]->GetKey();
        if(bcType == PHENGLEI::FARFIELD)
        {
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

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int bcType = bcRecord[iFace]->GetKey();
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


    //! Boundary faces.
    RDouble ** qOnBCFace = NewPointer2< RDouble >(nEquation + 1, nBoundFace);  // '1' is temperature.
    PHSPACE::SetField(qOnBCFace, nEquation + 1, nBoundFace, 0.0);

    //! Init q on BC faces.
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int bcType = bcRecord[iFace]->GetKey();
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

    //! Step1: solid surface.
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int bcType = bcRecord[iFace]->GetKey();

        if (bcType == PHENGLEI::SOLID_SURFACE)
        {
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

    //! Step2: Far-field.
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int bcType = bcRecord[iFace]->GetKey();

        if (bcType == PHENGLEI::FARFIELD)
        {
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

    //! Step3: other BC except solid surface/far field/symmetry/Interface.
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int bcType = bcRecord[iFace]->GetKey();

        if (bcType != PHENGLEI::SOLID_SURFACE && bcType != PHENGLEI::FARFIELD && 
            bcType != PHENGLEI::SYMMETRY      && bcType != PHENGLEI::INTERFACE)
        {
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

    //////////////////////////////////////////////////////////////////////////
    //! The velocity modify has not been considered on the symmetry!!!!!!!!!!
    //////////////////////////////////////////////////////////////////////////


    //! Change the number of points divided by the interpoint to the number of points before partition,because it is a face loop,each point is calculated twice,so the number of points needs to be multiplied by 2.
    InterpointInformation *interPointInfor = grid->GetInterpointInfo();
    if (interPointInfor)
    {
        int numberOfInterpoints = interPointInfor->GetNumberOfInterpoints();
        int *interPoint2GlobalPoint = interPointInfor->GetInterPoint2GlobalPoint();
        int *cellNumberOfInterPoint = interPointInfor->GetCellNumberOfInterPoint();
        int *labelOfInterPoint = interPointInfor->GetLabelOfInterPoint();
        int *totalZonesOfInterPoint = interPointInfor->GetTotalZonesOfInterPoint();

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

    int nEquation = nLaminar + nChemical;

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

    RDouble **q   = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("q") );
    RDouble *gama = reinterpret_cast< RDouble * > ( grid->GetDataPtr("gama") );
    RDouble *c = new RDouble[nTotal];

	using namespace IDX;

    for ( int iCell = 0; iCell < nTotal; ++ iCell )
	{
        c[iCell] = sqrt( gama[iCell] * q[IP][iCell] / q[IR][iCell] );
	}

    return c;
}

RDouble* NSSolverUnstruct::CompMachNumber(Grid *gridIn)
{
	UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    RDouble **q   = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("q") );
    RDouble *gama = reinterpret_cast< RDouble * > ( grid->GetDataPtr("gama") );
    RDouble *mach = new RDouble[nTotal];

	using namespace IDX;

    RDouble rm,um,vm,wm,pm,c2,v2;
    for ( int iCell = 0; iCell < nTotal; ++ iCell )
	{
        rm = q[IR][iCell];
        um = q[IU][iCell];
        vm = q[IV][iCell];
        wm = q[IW][iCell];
        pm = q[IP][iCell];
        c2 = abs(gama[iCell] * pm / rm);
        v2 = um * um + vm * vm + wm * wm;

        mach[iCell] = sqrt( v2 / c2 );
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

    RDouble **dqdx = reinterpret_cast <RDouble **> (gridIn->GetDataPtr("gradPrimtiveVarFaceX"));
    RDouble **dqdy = reinterpret_cast <RDouble **> (gridIn->GetDataPtr("gradPrimtiveVarFaceY"));
    RDouble **dqdz = reinterpret_cast <RDouble **> (gridIn->GetDataPtr("gradPrimtiveVarFaceZ"));

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

	for ( int iCell = 0; iCell < nTotal; ++ iCell ) 
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

/*
	RDouble * vorticityMagnitude = new RDouble [ nTotal ];

	if(!WantVisualField(grid))
	{
		//! Even if do not visual field vorticity, 
		//! it would be dumped on boundary, to compatible with the field variable.
		PHSPACE::SetField(vorticityMagnitude, -1.0e10, nTotal);
	}
*/
    RDouble **dqdx = reinterpret_cast <RDouble **> (gridIn->GetDataPtr("gradPrimtiveVarFaceX"));
    RDouble **dqdy = reinterpret_cast <RDouble **> (gridIn->GetDataPtr("gradPrimtiveVarFaceY"));
    RDouble **dqdz = reinterpret_cast <RDouble **> (gridIn->GetDataPtr("gradPrimtiveVarFaceZ"));
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

	for ( int iCell = 0; iCell < nTotal; ++ iCell ) 
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
		s12 = half * ( dudy + dvdx ); s13 = half * ( dudz + dwdx ); s23 = half * ( dvdz + dwdy );
		S = two * ( s11*s11 + s22*s22 + s33*s33 + two*(s12*s12 + s13*s13 + s23*s23) );
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

		Q_criteria[iCell] = fourth * ( Om * Om - S * S);
	}
}

RDouble * NSSolverUnstruct::ComputeCp(Grid *gridIn)
{
	UnstructGrid *grid = UnstructGridCast(gridIn);
	int nTotalCell = grid->GetNTotalCell();
	int nBoundFace = grid->GetNBoundFace();
	int nTotal = nTotalCell + nBoundFace;
	RDouble ** q   = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("q") );

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
	RDouble **qAverage = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("qAverage") );

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
    int nEquation = nLaminar + nChemical;

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

void NSSolverUnstruct::ZeroResiduals(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar = parameters->GetLaminarNumber();

    RDouble **res = reinterpret_cast< RDouble ** > (grid->GetDataPtr("res"));

    for (int m = 0; m < nLaminar; ++ m)
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

	if ( ktmax > 0 )
	{
		for ( int iCell = 0; iCell < nTotalCell; ++ iCell )
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
    int nBoundFace = grid->GetNBoundFace();

    RDouble **q   = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
    RDouble *gama = reinterpret_cast<RDouble *> (grid->GetDataPtr("gama"));
    RDouble *invSpectrum = reinterpret_cast<RDouble *> (grid->GetDataPtr("srad"));

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
		RDouble cl  = sqrt(gama[le] * pl / rl);
		RDouble vnl = xfn[iFace] * ul + yfn[iFace] * vl + zfn[iFace] * wl;

		RDouble rr  = rho[re];
		RDouble ur  = u[re];
		RDouble vr  = v[re];
		RDouble wr  = w[re];
		RDouble pr  = p[re];
		RDouble cr  = sqrt(gama[re] * pr / rr);
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

//! All in one.
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
	int nBoundFace = grid->GetNBoundFace();

	RDouble **q   = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
	RDouble *gama = reinterpret_cast< RDouble *  > (grid->GetDataPtr("gama"));
    RDouble *visSpectrum = reinterpret_cast< RDouble *  > (grid->GetDataPtr("srv"));

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
	RDouble *viscousTurbulent = reinterpret_cast< RDouble * > (grid->GetDataPtr("vist"));

    Param_NSSolverUnstruct *parameters = GetControlParameters();
	RDouble refReNumber = parameters->GetRefReNumber();

	RDouble prandtlLaminar    = parameters->GetPrandtlLaminar();    //! default 0.72
	RDouble prandtlTurbulence = parameters->GetPrandtlTurbulence(); //! default 0.9

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
			RDouble visTurb    = half * (viscousTurbulent[le] + viscousTurbulent[re]);
			RDouble density    = half * (rho[le] + rho[re]);
			RDouble viscosity  = visLaminar + visTurb;	

			RDouble specificHeatRatio = half * (gama[le] + gama[re]);
			RDouble coefOfSpectrum1 = fourthThird * (visLaminar + visTurb);
			RDouble coefOfSpectrum2 = specificHeatRatio * (visLaminar * oprl + visTurb * oprt);
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
	else if(visSpetrumMethod == 2)
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
			RDouble visTurb    = half * (viscousTurbulent[le] + viscousTurbulent[re]);
			RDouble density    = half * (rho[le] + rho[re]);
			RDouble viscosity  = visLaminar + visTurb;		
            //WriteLogFile("visTurb = ", viscousTurbulent[le], viscousTurbulent[re], true);

			RDouble faceArea  = area[iFace];
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
	RDouble * vol  = grid->GetCellVolume();

	int nTotalCell = grid->GetNTotalCell();
	int nTotalFace = grid->GetNTotalFace();
	int nBoundFace = grid->GetNBoundFace();

	RDouble **q          = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
	RDouble *gama        = reinterpret_cast<RDouble  *> (grid->GetDataPtr("gama"));
	RDouble **t          = reinterpret_cast<RDouble **> (grid->GetDataPtr("t"));
    RDouble *invSpectrum = reinterpret_cast<RDouble  *> (grid->GetDataPtr("srad"));

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
		RDouble pl	= p[le];
		RDouble cl  = sqrt(gama[le] * pl / rl);
		RDouble vnl = xfn[iFace] * ul + yfn[iFace] * vl + zfn[iFace] * wl;
		RDouble v2l	= ul * ul + vl * vl + wl * wl;

		RDouble rr  = rho[re];
		RDouble ur  = u[re];
		RDouble vr  = v[re];
		RDouble wr  = w[re];
		RDouble pr  = p[re];
		RDouble cr  = sqrt(gama[re] * pr / rr);
		RDouble vnr	= xfn[iFace] * ur + yfn[iFace] * vr + zfn[iFace] * wr;
		RDouble v2r = ur * ur + vr * vr + wr * wr;

		RDouble specificHeatRatio = half * ( gama[le] + gama[re] );
		RDouble pm  = half * ( pl + pr );
		RDouble rm  = half * ( rl + rr );
		RDouble cm  = sqrt( specificHeatRatio * pm / rm );
		RDouble vn  = half * ( vnl + vnr );	
		RDouble v2  = half * ( v2l + v2r );	
		RDouble c2  = cm * cm;	

		//RDouble beta = MAX( MIN( v2 / c2,1.0 ),3.0 * mach2 );
		RDouble beta = MIN( MAX( v2 / c2,3.0 * mach2 ),1.0 );
		RDouble vnew  = half * vn * ( 1 + beta );
		RDouble cnew  = half * sqrt( ( ( beta - 1 ) * vn ) * ( ( beta - 1 ) * vn ) + 4 * beta * c2 );
		RDouble inviscidSpectrumRadius = half * area[iFace] * ( ABS( vnew ) + cnew );

		invSpectrum[le] += inviscidSpectrumRadius;
		if(re < nTotalCell)
		{
			invSpectrum[re] += inviscidSpectrumRadius;
		}        
	}
    //grid->SetGhostCell(invSpectrum);
}

void NSSolverUnstruct::LocalTimeStep(Grid *gridIn)
{
	UnstructGrid *grid = UnstructGridCast(gridIn);

	RDouble *volume  = grid->GetCellVolume();

	int nTotalCell = grid->GetNTotalCell();

	RDouble **q   = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
	RDouble *gama = reinterpret_cast< RDouble *  > (grid->GetDataPtr("gama"));
	RDouble *dt   = reinterpret_cast< RDouble *  > (grid->GetDataPtr("dt"));
    RDouble *srad = reinterpret_cast<RDouble *> (grid->GetDataPtr("srad"));
    RDouble *srv  = reinterpret_cast<RDouble *> (grid->GetDataPtr("srv"));

    Param_NSSolverUnstruct *parameters = GetControlParameters();
	RDouble cfl = ComputeCFL();

	int isInMGIniting = GlobalDataBase::GetIntParaFromDB("isInMGIniting");
	if(!grid->IsFinestGrid() && !isInMGIniting)
	{
		RDouble mgCFLScale = parameters->GetMgCFLScale();
		cfl *= mgCFLScale;
	}

	//! Set dt to zero.
    PHSPACE::SetField(dt, 0.0, nTotalCell);

	//! Inviscid time step.
	for (int iCell = 0; iCell < nTotalCell; ++ iCell)
	{
		dt[iCell] = cfl * volume[iCell] / srad[iCell];
	}

    int isViscous = parameters->IsViscous();

	if (isViscous)
	{
		//! Viscous time step.
        RDouble *dtViscous = new RDouble [nTotalCell];
		for (int iCell = 0; iCell < nTotalCell; ++ iCell)
		{
			dtViscous[iCell] = cfl * volume[iCell] / srv[iCell];
		}

		//! Consider the viscous time step.
		for (int iCell = 0; iCell < nTotalCell; ++ iCell)
		{
			dt[iCell] = dt[iCell] * dtViscous[iCell] / (dt[iCell] + dtViscous[iCell]);
		}

		delete []dtViscous;
	}
}

void NSSolverUnstruct::GlobalTimeStep(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

    RDouble * vol = grid->GetCellVolume();
	RDouble * dt  = reinterpret_cast< RDouble * > ( grid->GetDataPtr("dt"  ) );

    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();

    RDouble dtau;
    GlobalDataBase::GetData("dtau",&dtau,PHDOUBLE,1);

    //! Set dt to zero.
    for ( int iCell = 0; iCell < nTotalCell; ++ iCell )
    {
        dt[iCell] = dtau / vol[iCell];
    }

	grid->SetGhostCell(dt);
    RDouble dtmin = dtau;
    GlobalDataBase::UpdateData("dtmin",&dtmin,PHDOUBLE,1);
}

void NSSolverUnstruct::LocalGlobalTimeStep(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int      * leftCellofFace = grid->GetLeftCellOfFace();
    int      * rightCellofFace = grid->GetRightCellOfFace();

    RDouble * xfn  = grid->GetFaceNormalX();
    RDouble * yfn  = grid->GetFaceNormalY();
    RDouble * zfn  = grid->GetFaceNormalZ();
    RDouble * area = grid->GetFaceArea();
	RDouble * vgn  = grid->GetFaceNormalVelocity();
    RDouble * vol  = grid->GetCellVolume();

    int nTotalCell = grid->GetNTotalCell();
    int nTotalFace = grid->GetNTotalFace();
    int nBoundFace = grid->GetNBoundFace();

    RDouble **q   = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("q") );
    RDouble *gama = reinterpret_cast< RDouble *  > ( grid->GetDataPtr("gama") );
    RDouble *dt   = reinterpret_cast< RDouble *  > ( grid->GetDataPtr("dt") );

    RDouble cfl = ComputeCFL();

    //RDouble refGama = parameters->GetRefGama();

    int le,re;
    RDouble um, vm, wm, cm;

    RDouble *rho, *u, *v, *w, *p;

	using namespace IDX;

    rho = q[IR];
    u   = q[IU];
    v   = q[IV];
    w   = q[IW];
    p   = q[IP];

    //! Set dt to zero.
    for ( int iCell = 0; iCell < nTotalCell; ++ iCell )
    {
        dt[iCell] = 0.0;
    }

    BCRecord ** bcr = grid->GetBCRecord();

    for ( int iFace = 0; iFace < nBoundFace; ++ iFace )
    {
        le = leftCellofFace[ iFace ];

        um      = u[le];
        vm      = v[le];
        wm      = w[le];
        cm      = sqrt( gama[le] * p[le] / rho[le] );
        dt[le] += ( ABS( xfn[iFace] * um + yfn[iFace] * vm + zfn[iFace] * wm - vgn[iFace] ) + cm ) * area[iFace];
    }

    //! For interior faces.
    for ( int iFace = nBoundFace; iFace < nTotalFace; ++ iFace )
    {
        le = leftCellofFace[ iFace ];
        re = rightCellofFace[ iFace ];

        um      = u[le];
        vm      = v[le];
        wm      = w[le];
        cm      = sqrt( gama[le] * p[le] / rho[le] );
        dt[le] += ( ABS( xfn[iFace] * um + yfn[iFace] * vm + zfn[iFace] * wm - vgn[iFace] ) + cm ) * area[iFace];

        um      = u[re];
        vm      = v[re];
        wm      = w[re];
        cm      = sqrt( gama[re] * p[re] / rho[re] );
        dt[re] += ( ABS( xfn[iFace] * um + yfn[iFace] * vm + zfn[iFace] * wm - vgn[iFace] ) + cm ) * area[iFace];
    }

    for ( int iCell = 0; iCell < nTotalCell; ++ iCell )
    {
        dt[iCell] = cfl * vol[iCell] / dt[iCell];
    }

	Param_NSSolverUnstruct * parameters = GetControlParameters();
    int viscousType = parameters->GetViscousType();

    if (viscousType > 0)
    {
        //cout << " viscousType > 0 \n";
        RDouble * dtv;
        dtv = new RDouble [ nTotalCell ];

        for ( int iCell = 0; iCell < nTotalCell; ++ iCell )
        {
            dtv[iCell] = 0.0;
        }

        for ( int iFace = 0; iFace < nBoundFace; ++ iFace )
        {
            le = leftCellofFace[ iFace ];
            re = rightCellofFace[ iFace ];

            dtv[le] += area[iFace] * area[iFace];
        }

        //! For interior faces.
        for ( int iFace = nBoundFace; iFace < nTotalFace; ++ iFace )
        {
            le = leftCellofFace[ iFace ];
            re = rightCellofFace[ iFace ];

            dtv[le] += area[iFace] * area[iFace];
            dtv[re] += area[iFace] * area[iFace];
        }

        RDouble * visl = reinterpret_cast< RDouble * > ( grid->GetDataPtr("visl") );
        RDouble * vist = reinterpret_cast< RDouble * > ( grid->GetDataPtr("vist") );

        RDouble refReNumber  = parameters->GetRefReNumber();
        RDouble oRefReNumber = parameters->GetoRefReNumber();
		RDouble prandtlLaminar    = parameters->GetPrandtlLaminar();
		RDouble prandtlTurbulence = parameters->GetPrandtlTurbulence();

		RDouble coef1,coef2,coef,vis_l,vis_t;
		RDouble csrv = 1.0, csrvist = 1.0;

        const RDouble foth = 4.0 / 3.0;

	    RDouble oprl = 1.0 / prandtlLaminar;
	    RDouble oprt = 1.0 / prandtlTurbulence;

        for ( int iCell = 0; iCell < nTotalCell; ++ iCell )
        {
			vis_l = visl[iCell];
			vis_t = csrvist * vist[iCell];
			coef1 = foth * ( vis_l + vis_t );
			coef2 = gama[iCell] * ( vis_l * oprl + vis_t * oprt );
            coef  = MAX( coef1, coef2 );
            dtv[iCell] = cfl * rho[iCell] * refReNumber * vol[iCell] * vol[iCell] / ( coef * dtv[iCell] );
        }

        for ( int iCell = 0; iCell < nTotalCell; ++ iCell )
        {
            dt[iCell] = dt[iCell] * dtv[iCell] / ( dt[iCell] + dtv[iCell] );
        }

        delete [] dtv; 
    }

    RDouble dtmin = LARGE;
    for ( int iCell = 0; iCell < nTotalCell; ++ iCell )
    {
        dtmin = MIN( dtmin, static_cast<RDouble>( dt[iCell] ) );
    }

    GlobalDataBase::UpdateData("dtmin",&dtmin,PHDOUBLE,1);

    for ( int iCell = 0; iCell < nTotalCell; ++ iCell )
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

    int nEquation = nLaminar + nChemical;

	RDouble **field = field_proxy->GetField_UNS();

    for ( int m = 0; m < nEquation; ++ m )
    {
        for ( int iCell = 0; iCell < nTotal; ++ iCell )
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
    int nEquation = nLaminar + nChemical;

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
    int nEquation = nLaminar + nChemical;

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

    int nEquation = nLaminar + nChemical;

    RDouble **fq = reinterpret_cast< RDouble ** > ( fgrid->GetDataPtr("q") );
    RDouble **cq = reinterpret_cast< RDouble ** > ( cgrid->GetDataPtr("q") );

    for ( int m = 0; m < nEquation; ++ m )
    {
        RestrictQ(fgrid, fq[m], cgrid, cq[m]);
    }

	RDouble *fgama = reinterpret_cast< RDouble * > ( fgrid->GetDataPtr("gama") );
	RDouble *cgama = reinterpret_cast< RDouble * > ( cgrid->GetDataPtr("gama") );

	RestrictQ(fgrid, fgama, cgrid, cgama);
	
	int viscousType = parameters->GetViscousType();
    if ( viscousType > 0 )
    {
		RDouble * fvisl = reinterpret_cast< RDouble * > ( fgrid->GetDataPtr("visl") );
        RDouble * cvisl = reinterpret_cast< RDouble * > ( cgrid->GetDataPtr("visl") );
        RestrictQ(fgrid, fvisl, cgrid, cvisl);

		if ( viscousType > 1 )
		{
			RDouble * fvist = reinterpret_cast< RDouble * > ( fgrid->GetDataPtr("vist") );
			RDouble * cvist = reinterpret_cast< RDouble * > ( cgrid->GetDataPtr("vist") );
			RestrictQ(fgrid, fvist, cgrid, cvist);
		}
	}

	int isUnsteady = parameters->GetIsUnsteady();
	if (isUnsteady)
	{
		RDouble **fqn1 = reinterpret_cast< RDouble ** > ( fgrid->GetDataPtr("q_unsteady_n1") );
		RDouble **fqn2 = reinterpret_cast< RDouble ** > ( fgrid->GetDataPtr("q_unsteady_n2") );
		RDouble **cqn1 = reinterpret_cast< RDouble ** > ( cgrid->GetDataPtr("q_unsteady_n1") );
		RDouble **cqn2 = reinterpret_cast< RDouble ** > ( cgrid->GetDataPtr("q_unsteady_n2") );

		RDouble **fres   = reinterpret_cast< RDouble ** > ( fgrid->GetDataPtr("res"            ) );
		RDouble **fresn1 = reinterpret_cast< RDouble ** > ( fgrid->GetDataPtr("res_unsteady_n1") );
		RDouble **fresn2 = reinterpret_cast< RDouble ** > ( fgrid->GetDataPtr("res_unsteady_n2") );
		RDouble **cres   = reinterpret_cast< RDouble ** > ( cgrid->GetDataPtr("res"            ) );
		RDouble **cresn1 = reinterpret_cast< RDouble ** > ( cgrid->GetDataPtr("res_unsteady_n1") );
		RDouble **cresn2 = reinterpret_cast< RDouble ** > ( cgrid->GetDataPtr("res_unsteady_n2") );

		for ( int m = 0; m < nEquation; ++ m )
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
    int nTotalCell = grid->GetNTotalCell();
    int nTotal = nTotalCell + nBoundFace;

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
        for ( int m = 0; m < nLaminar; ++ m )
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
    if (viscousType == 0) return;

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
        for ( int m = 0; m < nLaminar; ++ m )
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
    if (viscousType == 0) return;

    RDouble **Viscousflux = ViscousfaceProxyOnlyForMixGrid->GetFlux();
    //! For boundary faces, remember re is ghost cell.
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int le = leftCellofFace[iFace];
        for ( int m = 0; m < nLaminar; ++ m )
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
    ComputeViscousCoefficientWithoutChemical(gridIn);
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
    GlobalDataBase::GetData("visl_min",&viscousLaminarMin,PHDOUBLE,1);

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

    RDouble **fq = reinterpret_cast< RDouble ** > ( fgrid->GetDataPtr("q") );
    RDouble **cq = reinterpret_cast< RDouble ** > ( cgrid->GetDataPtr("q") );

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();

    int nEquation = nLaminar + nChemical;
    int cn = cnTCell + cgrid->GetNBoundFace();

    for ( int m = 0; m < nEquation; ++ m )
    {
        for ( int iCell = 0; iCell < cn; ++ iCell )
        {
            cq[m][iCell] = 0.0;
        }
    }

    for ( int m = 0; m < nEquation; ++ m )
    {
        for ( int iCell = 0; iCell < fnTotalCell; ++ iCell )
        {
            cq[m][ cell2coarsegridcell[iCell] ] += fvol[iCell] * fq[m][iCell];
        }
    }
    for ( int m = 0; m < nEquation; ++ m )
    {
        for ( int iCell = 0; iCell < cnTCell; ++ iCell )
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

    int nBoundFace = grid->GetNBoundFace();
    int nTotalCell = grid->GetNTotalCell();
    int nTotalFace = grid->GetNTotalFace();
    int nTotal = nTotalCell + nBoundFace;

    int **cell2face = grid->GetCell2Face();
    int *faceNumberOfEachCell = grid->GetFaceNumberOfEachCell();

    int *leftCellofFace = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();

    RDouble **q   = reinterpret_cast <RDouble **> (grid->GetDataPtr("q"));
    RDouble *rtem = reinterpret_cast <RDouble *> (grid->GetDataPtr("rtem"));

    RDouble pressureSensorCoefficient = third;  // for 3-D
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

            LeftPressureOfCell  += q[PHSPACE::IDX::IP][le];
            RightPressureOfCell += q[PHSPACE::IDX::IP][re];
        }
        rtem[iCell] = ABS(RightPressureOfCell - LeftPressureOfCell) / ABS(RightPressureOfCell + LeftPressureOfCell);
        rtem[iCell] /= static_cast <RDouble> (faceNumberOfEachCell[iCell]);
    }
    grid->SetGhostCellExceptInterface(rtem);
}

void NSSolverUnstruct::ViscousFlux(Grid *gridIn)
{
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int viscousType = parameters->GetViscousType();
    if (viscousType == 0) return;

    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nBoundFace = grid->GetNBoundFace();
    int nTotalCell = grid->GetNTotalCell();
    int nTotalFace = grid->GetNTotalFace();
    int nTotal = nTotalCell + nBoundFace;

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
    int nTotal = nTotalCell + nBoundFace;

    NSFaceValue *faceVariable = faceProxy->GetNSFaceValue();

    RDouble * xfc  = grid->GetFaceCenterX();
    RDouble * yfc  = grid->GetFaceCenterY();
    RDouble * zfc  = grid->GetFaceCenterZ();
    RDouble * xcc  = grid->GetCellCenterX();
    RDouble * ycc  = grid->GetCellCenterY();
    RDouble * zcc  = grid->GetCellCenterZ();

	RDouble * xtn  = grid->GetNodeVelocityX();
	RDouble * ytn  = grid->GetNodeVelocityY();
	RDouble * ztn  = grid->GetNodeVelocityZ();

    RDouble *area = grid->GetFaceArea();
    int *leftCellOfFace = grid->GetLeftCellOfFace();
    int *rightCellOfFace = grid->GetRightCellOfFace();

	RDouble **q   = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));
	RDouble **t   = reinterpret_cast<RDouble **> (grid->GetDataPtr("t"));

    RDouble *viscousLaminar = reinterpret_cast<RDouble *> (grid->GetDataPtr("visl"));
	RDouble *viscousTurbulent = reinterpret_cast<RDouble *> (grid->GetDataPtr("vist"));

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    RDouble refGama = 1.4;
    GlobalDataBase::GetData("refGama", &refGama, PHDOUBLE, 1);
    //RDouble refGama = parameters->GetRefGama();

    RDouble oPrandtlLaminar    = parameters->GetoPrandtlLaminar();
    RDouble oPrandtlTurbulence = parameters->GetoPrandtlTurbulence();

    RDouble referenceMachNumber = parameters->GetRefMachNumber();

    RDouble tWall = parameters->GetWallTemperature();

    RDouble refDimensionalTemperature = parameters->GetRefDimensionalTemperature();

    int numberOfSpecies = parameters->GetNumberOfSpecies();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nEquation = nLaminar + nChemical;

    int nm = parameters->GetNSEquationNumber();

    RDouble **rhoDiffusion = faceVariable->GetRhoDS()->GetField();
    RDouble **hintSpecies = faceVariable->GetHintS()->GetField();
    RDouble **primitiveVariable = faceVariable->GetPrim()->GetField();
    RDouble **tMid   = faceVariable->GetT()->GetField();

    RDouble *kCp = faceVariable->GetKCP();   //! Heat conductivity coefficient, usually computed by Cp.
    RDouble *viscousLaminarFace = faceVariable->GetMUL();
    RDouble *viscousTurbulentFace = faceVariable->GetMUT();

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
        viscousTurbulentFace[jFace] = half * (viscousTurbulent[le] + viscousTurbulent[re]);
    }

    int nMid = nBoundFace;
    if ( localEnd < nBoundFace ) nMid = localEnd;

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

            gas->ComputeDensityDiffusionAndKCP(t1, primface, viscousLaminarFace[jFace], viscousTurbulentFace[jFace], rhoDiffusionFace, hintSpeciesOfFace, kCp[jFace]);
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
            //RDouble cp = 1.0 / ((refGama - 1.0) * referenceMachNumberSquare);
            kCp[jFace] = (viscousLaminarFace[jFace] * oPrandtlLaminar + viscousTurbulentFace[jFace] * oPrandtlTurbulence) * cp;
        }
    }

    //! If there are some boundary faces, make modifications. 
    //! Check if there are boundary faces. If no, return.
    if (localStart < nBoundFace)
    {
        BCRecord **bcRecord = grid->GetBCRecord();
        RDouble tw = tWall / refDimensionalTemperature;
        for (int iFace = localStart; iFace < nMid; ++ iFace)
        {
            int bcType = bcRecord[iFace]->GetKey();
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
                Data_Param *bcData  = bcRecord[iFace]->GetBCRegion()->GetBCDataParam();
                if(bcData)
                {
                    if (bcData->IsExist("uWall", PHDOUBLE, 1))
                    {
                        bcData->GetData("uWall", &uWall, PHDOUBLE, 1);
                        bcData->GetData("vWall", &vWall, PHDOUBLE, 1);
                        bcData->GetData("wWall", &wWall, PHDOUBLE, 1);
                    }
                }
                

                primitiveVariable[IU][jFace] = xtn[iFace] + uWall;
                primitiveVariable[IV][jFace] = ytn[iFace] + vWall;
                primitiveVariable[IW][jFace] = ztn[iFace] + wWall;

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
    int nTotalCell = grid->GetNTotalCell();
    int nTotal = nTotalCell + nBoundFace;

	int *leftCellofFace = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();

	BCRecord **bcRecord = grid->GetBCRecord();

    RDouble **dqdx = reinterpret_cast <RDouble **> (gridIn->GetDataPtr("gradPrimtiveVarFaceX"));
    RDouble **dqdy = reinterpret_cast <RDouble **> (gridIn->GetDataPtr("gradPrimtiveVarFaceY"));
    RDouble **dqdz = reinterpret_cast <RDouble **> (gridIn->GetDataPtr("gradPrimtiveVarFaceZ"));

    RDouble **dtdx = reinterpret_cast <RDouble **> (gridIn->GetDataPtr("gradTemperatureFaceX"));
    RDouble **dtdy = reinterpret_cast <RDouble **> (gridIn->GetDataPtr("gradTemperatureFaceY"));
    RDouble **dtdz = reinterpret_cast <RDouble **> (gridIn->GetDataPtr("gradTemperatureFaceZ"));

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
	RDouble *viscousTurbulentFace = faceVariable->GetMUT();

	Param_NSSolverUnstruct *parameters = GetControlParameters();
    int numberOfSpecies = parameters->GetNumberOfSpecies();

    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nEquation = nLaminar + nChemical;

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

    if ( nrokplus > 0 )
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
                int bcType = bcRecord[iFace]->GetKey();

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

        RDouble vis = viscousLaminarFace[jFace] + viscousTurbulentFace[jFace];
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

        for ( int m = 0; m < nLaminar; ++ m )
        {
            flux[m][jFace] = - oRefReNumber * area[iFace] * fvis[m];
        }
    }

    delete []fvis;

    delete []dfdx;
    delete []dfdy;
    delete []dfdz;

    delete []dfd1;
    delete []dfd2;
    delete []dfdn;
    delete []dfdt1;
    delete []dfdt2;

    delete []f1;
    delete []f2;
    delete []fMid;
}

void NSSolverUnstruct::Boundary(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);

	Param_NSSolverUnstruct *parameters = GetControlParameters();
    RDouble wallTemperature = parameters->GetWallTemperature();
	int ifLowSpeedPrecon = parameters->GetIfLowSpeedPrecon();

    using namespace PHENGLEI;
	
    int isViscous = parameters->IsViscous();

    BCRegionUnstruct **bcRegionUnstruct = grid->GetBCRegionUnstruct();
    int nBCRegionUnstruct = grid->GetNBCRegionUnstruct();
    for (int iBCRegionUnstruct = 0; iBCRegionUnstruct < nBCRegionUnstruct; iBCRegionUnstruct++)
    {
        int bcType = bcRegionUnstruct[iBCRegionUnstruct]->GetBCType();
        if (bcType < 0)
        {
            continue;
        }
        else if (bcType == EXTRAPOLATION)
        {
            OutflowBCRegion(grid, bcRegionUnstruct[iBCRegionUnstruct]);
        }
        else if (IsWall(bcType))
        {
            if (!isViscous)
            {
                SymmetryBCRegion(grid, bcRegionUnstruct[iBCRegionUnstruct]);
            }
            else if (wallTemperature <= 0.0)
            {
                ViscousAdiabaticWallBCRegion(gridIn,bcRegionUnstruct[iBCRegionUnstruct]);
            }
            else
            {
                ViscousIsotropicWallBCRegion(grid, bcRegionUnstruct[iBCRegionUnstruct]);
            }
        }
        else if (bcType == SYMMETRY)
        {
            SymmetryBCRegion(grid, bcRegionUnstruct[iBCRegionUnstruct]);
        }
        else if (bcType == FARFIELD)
        {
            if (ifLowSpeedPrecon == 0)
            {
                FarfieldBCRegion(grid, bcRegionUnstruct[iBCRegionUnstruct]);
            }
            else
            {
                PreconFarfieldBCRegion(grid, bcRegionUnstruct[iBCRegionUnstruct]);
            }
        }
        else if (bcType == INFLOW)
        {
            InflowBCRegion(gridIn, bcRegionUnstruct[iBCRegionUnstruct]);
        }
        else if (bcType == OUTFLOW)
        {
            OutflowBCRegion(gridIn, bcRegionUnstruct[iBCRegionUnstruct]);
        }
        else if (bcType == PRESSURE_OUTLET)
        {
            PressureOutletBCRegion(gridIn, bcRegionUnstruct[iBCRegionUnstruct]);
        }
        else if (bcType == PRESSURE_INLET)
        {
            PressureInletBCRiemannRegion(gridIn, bcRegionUnstruct[iBCRegionUnstruct]);
        }
        else if (bcType == MASS_FLOW_INLET)
        {
            MassFlowInletBCRegion(gridIn, bcRegionUnstruct[iBCRegionUnstruct]);
        }
        else if (bcType == MASS_FLOW_OUTLET)
        {
            MassFlowOutletBCRegion(gridIn, bcRegionUnstruct[iBCRegionUnstruct]);
        }
        else if ( bcType == OVERSET )
        {
            FarfieldBCRegion(gridIn, bcRegionUnstruct[iBCRegionUnstruct]);
        }
        else if (ifIsJetBC(bcType))
        {
            JetBC3DRegion(gridIn, bcRegionUnstruct[iBCRegionUnstruct]);
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
    int       *leftCellofFace   = grid->GetLeftCellOfFace();
    int       *rightCellofFace   = grid->GetRightCellOfFace();
    RDouble   *xfn    = grid->GetFaceNormalX();
    RDouble   *yfn    = grid->GetFaceNormalY();
    RDouble   *zfn    = grid->GetFaceNormalZ();
	RDouble   *vgn    = grid->GetFaceNormalVelocity();
    RDouble   *xtn    = grid->GetNodeVelocityX();
    RDouble   *ytn    = grid->GetNodeVelocityY();
    RDouble   *ztn    = grid->GetNodeVelocityZ();
    BCRecord  **bcr   = grid->GetBCRecord();

    int nBoundFace    = grid->GetNBoundFace();
    int nTotalCell    = grid->GetNTotalCell();

    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
	RDouble *gama = reinterpret_cast< RDouble * > (grid->GetDataPtr("gama"));

	Param_NSSolverUnstruct * parameters = GetControlParameters();
    int nm    = parameters->GetNSEquationNumber();

    int nLaminar  = parameters->GetLaminarNumber();
	int nChemical = parameters->GetChemicalFlag();
    
    RDouble refGama = parameters->GetRefGama();

    RDouble wallTemperature = parameters->GetWallTemperature();
    RDouble refDimensionalTemperature = parameters->GetRefDimensionalTemperature();
	RDouble tw = wallTemperature / refDimensionalTemperature;

    RDouble *prims = new RDouble[nLaminar+nChemical];
    RDouble *primt = new RDouble[nLaminar+nChemical];

	using namespace IDX;
    using namespace PHENGLEI;

	RDouble massIn  = 0.0;
	RDouble massOut = 0.0;

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int le = leftCellofFace [iFace];
        int re = rightCellofFace[iFace];
        int bcType = bcr[iFace]->GetKey();

        for (int m = 0; m < nLaminar + nChemical; ++ m)
        {
            prims[m] = q[m][le];
            primt[m] = q[m][re];
        }

        if (bcType < 0)
        {
            continue;
        }
        else if (bcType == INFLOW)
        {
			RDouble rmi, umi, vmi, wmi, vni;
            RDouble rmo, umo, vmo, wmo, vno;

			rmi = prims[IR];
			umi = prims[IU];
			vmi = prims[IV];
			wmi = prims[IW];
			vni = xfn[iFace] * umi + yfn[iFace] * vmi + zfn[iFace] * wmi;

			rmo = primt[IR];
			umo = primt[IU];
			vmo = primt[IV];
			wmo = primt[IW];
			vno = xfn[iFace] * umo + yfn[iFace] * vmo + zfn[iFace] * wmo;

			massIn -= half * (rmi * vni + rmo * vno);
        }
        else if (bcType == OUTFLOW)
        {
			RDouble rm, um, vm, wm, vn;
			rm = prims[IR];
			um = prims[IU];
			vm = prims[IV];
			wm = prims[IW];
			vn = xfn[iFace] * um + yfn[iFace] * vm + zfn[iFace] * wm;

			massOut += rm * vn;
        }
    }

    delete [] prims;
    delete [] primt;

    GlobalDataBase::UpdateData("massIn" ,&massIn ,PHDOUBLE,1);
	GlobalDataBase::UpdateData("massOut",&massOut,PHDOUBLE,1);
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
    RDouble  *vgn  = grid->GetFaceNormalVelocity();
    RDouble  *xtn  = grid->GetNodeVelocityX();
    RDouble  *ytn  = grid->GetNodeVelocityY();
    RDouble  *ztn  = grid->GetNodeVelocityZ();
    BCRecord **bcRecord = grid->GetBCRecord();

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nEquation = nLaminar + nChemical;

    int isViscous = parameters->IsViscous();
    RDouble tWall = parameters->GetWallTemperature();
    RDouble refDimensionalTemperature = parameters->GetRefDimensionalTemperature();
    RDouble refGama = parameters->GetRefGama();
    RDouble coefficientOfStateEquation = gas->GetCoefficientOfStateEquation();

    using namespace PHSPACE;
    using namespace PHENGLEI;

    RDouble **qL = faceProxy->GetQL();
    RDouble **qR = faceProxy->GetQR();

    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();
    RDouble **q = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));

    RDouble *primitiveVariable = new RDouble[nLaminar + nChemical];

    using namespace IDX;

    for (int iFace = localStart; iFace < nMid; ++ iFace)
    {
        int BCType = bcRecord[iFace]->GetKey();
        int jFace = iFace - localStart;

        if (BCType != PHENGLEI::SOLID_SURFACE && BCType != PHENGLEI::SYMMETRY)
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

        if (BCType == PHENGLEI::SYMMETRY)
        {
            RDouble vn = 2.0 * (xfn[iFace] * qL[IU][jFace] + yfn[iFace] * qL[IV][jFace] + zfn[iFace] * qL[IW][jFace]);

            qR[IR][jFace] = qL[IR][jFace];
            qR[IU][jFace] = qL[IU][jFace] - vn * xfn[iFace];
            qR[IV][jFace] = qL[IV][jFace] - vn * yfn[iFace];
            qR[IW][jFace] = qL[IW][jFace] - vn * zfn[iFace];
            qR[IP][jFace] = qL[IP][jFace];
        }

        if (BCType == PHENGLEI::SOLID_SURFACE && isViscous)
        {
            RDouble uWall = 0.0;
            RDouble vWall = 0.0;
            RDouble wWall = 0.0;                
            Data_Param *bcData  = bcRecord[iFace]->GetBCRegion()->GetBCDataParam();
            if(bcData)
            {
                if (bcData->IsExist("uWall", PHDOUBLE, 1))
                {
                    bcData->GetData("uWall", &uWall, PHDOUBLE, 1);
                    bcData->GetData("vWall", &vWall, PHDOUBLE, 1);
                    bcData->GetData("wWall", &wWall, PHDOUBLE, 1);
                }
            }            

            qL[IU][jFace] = xtn[iFace] + uWall;
            qL[IV][jFace] = ytn[iFace] + vWall;
            qL[IW][jFace] = ztn[iFace] + wWall;

            qR[IU][jFace] = xtn[iFace] + uWall;
            qR[IV][jFace] = ytn[iFace] + vWall;
            qR[IW][jFace] = ztn[iFace] + wWall;

            if (tWall <= 0.0)
            {
                //! Viscous WALL, adiabatic.
            }
            else
            {
                //! iso-thermal wall.
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

    ReConstructFaceValue(grid, faceProxy, limiter, localStart, localEnd);

    BoundaryQlQrFix(grid, faceProxy, localStart, localEnd);
}

/************************************************************************/
/*        Face values reconstruction by considering scalar limiters     */
/*							 Bell 20121227 modify						*/
/************************************************************************/
void NSSolverUnstruct::ReConstructFaceValue(Grid *gridIn, FaceProxy *faceProxy, Limiter * limiter, int localStart, int localEnd)
{
	UnstructGrid *grid = UnstructGridCast(gridIn);

	int nBoundFace = grid->GetNBoundFace();
	int nTotalCell = grid->GetNTotalCell();

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
    int nEquation = nLaminar + nChemical;

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

    RDouble **dqdx = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarFaceX"));
    RDouble **dqdy = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarFaceY"));
    RDouble **dqdz = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarFaceZ"));

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
    int precondition = parameters->GetIfLowSpeedPrecon();

    INVScheme inviscidScheme = GetINVScheme(unstructScheme);

    if ((!grid->IsFinestGrid() && isInMGIniting && !precondition) ||
        ( grid->IsFinestGrid() && isInIniting   && !precondition) )
    {
        //! The precondition has not been considered for vanleer method.
        //! so, only in the compressible multi-grid initialization state, run Vanleer.
        inviscidScheme = Vanleer_Scheme;

        int timeIntegration = parameters->GetTimeIntegration();
        if (timeIntegration == Block_LU_SGS)
        {
            //! Block LUSGS is only suitable for transonic flow. so if we takes BLUSGS, 
            //! the inviscid need to take Roe scheme.
            inviscidScheme = Roe_Scheme_ConservativeForm;
        }
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

	std::copy(xfn  + localStart, xfn  + localEnd, geom_proxy->GetFaceNormalX() );
	std::copy(yfn  + localStart, yfn  + localEnd, geom_proxy->GetFaceNormalY() );
	std::copy(zfn  + localStart, zfn  + localEnd, geom_proxy->GetFaceNormalZ() );
	std::copy(vgn  + localStart, vgn  + localEnd, geom_proxy->GetFaceVelocity());
	std::copy(area + localStart, area + localEnd, geom_proxy->GetFaceArea()    );
}

void NSSolverUnstruct::GetResidual(ActionKey *actkey)
{
    int level = actkey->level;
    UnstructGrid *grid = UnstructGridCast(GetGrid(level));

    RDouble **res = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("res") );

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int ntmodel = parameters->GetTemperatureModel();
    int nEquation = nLaminar + nChemical + ntmodel - 1;

    PHSPACE::ComputeResidualonGrid(grid, actkey, res, nEquation);
}

void NSSolverUnstruct::InviscidFluxWrap(Grid *gridIn, FaceProxy *faceProxy, INVScheme inviscidScheme, int localStart, int localEnd)
{
	UnstructGrid *grid = UnstructGridCast(gridIn);

    Param_NSSolverUnstruct *parameters = GetControlParameters();
	int nLaminar = parameters->GetLaminarNumber();
	int nFaceLength = localEnd - localStart;
	faceProxy->setsize(nFaceLength);
    RDouble **res = reinterpret_cast<RDouble**> (grid->GetDataPtr("res"));

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

		for ( int m = 0; m < nLaminar; ++ m )
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

	Param_NSSolverUnstruct * parameters = GetControlParameters();
	int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();

	int nEquation = nLaminar + nChemical;

	CFDSolver::RegisterOversetField( "dqdx", PHDOUBLE, nEquation );
	CFDSolver::RegisterOversetField( "dqdy", PHDOUBLE, nEquation );
	CFDSolver::RegisterOversetField( "dqdz", PHDOUBLE, nEquation );

	CFDSolver::RegisterOversetField( "limit", PHDOUBLE, 1 );

	int ntmodel = parameters->GetTemperatureModel();
	CFDSolver::RegisterOversetField( "dtdx", PHDOUBLE, ntmodel );
	CFDSolver::RegisterOversetField( "dtdy", PHDOUBLE, ntmodel );
	CFDSolver::RegisterOversetField( "dtdz", PHDOUBLE, ntmodel );

	CFDSolver::RegisterOversetField( "Qcell", PHDOUBLE, 1 );

	CFDSolver::RegisterOversetField( "dq", PHDOUBLE, nEquation );
}

void NSSolverUnstruct::SolutionFix(Grid *gridIn, RDouble *primitiveVariable, int iCell)
{
	UnstructGrid *grid = UnstructGridCast(gridIn);

	RDouble **q = reinterpret_cast<RDouble **> (grid->GetDataPtr("q"));

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nm = parameters->GetNSEquationNumber();
	int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nEquation = nLaminar + nChemical;

    int **cell2face  = grid->GetCell2Face();
    int *faceNumberOfEachCell = grid->GetFaceNumberOfEachCell();

    int *leftCellofFace = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();

    using namespace GAS_SPACE;
    using namespace IDX;

    int np = 0;

	for (int m = 0; m < nEquation; ++ m)
	{
		primitiveVariable[m] = zero;
	}

	int iFace, le, re, isur;
	RDouble volumeSum, ovolumeSum, vsur;
	volumeSum = 0.0;
	for (int icell2face = 0; icell2face < faceNumberOfEachCell[iCell]; ++ icell2face)
	{
		iFace  = cell2face[iCell][icell2face];
		le    = leftCellofFace[iFace];
		re    = rightCellofFace[iFace];

        isur = le;
        if (iCell == le) isur = re;

		//vsur = vol[isur];
		vsur = one;

		volumeSum += vsur;

		for (int m = 0; m < nEquation; ++ m)
		{
			RDouble f = q[m][isur];
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
    RDouble *vol  = grid->GetCellVolume();

    int *leftCellofFace = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();

    Param_NSSolverUnstruct *parameters = GetControlParameters();

    int nm      = parameters->GetNSEquationNumber();

    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int ntmodel = parameters->GetTemperatureModel();
    int nEquation = nLaminar + nChemical + ntmodel - 1;

    RDouble refGama = parameters->GetRefGama();
    RDouble refMachNumber = parameters->GetRefMachNumber();

    RDouble rm, pm, tm, tv, te;

    int nmodcell = 0;
    int nnegtive = 0;

    int **cell2face  = grid->GetCell2Face();
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

    RDouble p0 = P0Wall(poo, refGama, refMachNumber);
    RDouble r0 = R0Wall(roo, refGama, refMachNumber);

    int nNegativeCell           = 0;
    int mostNegativeCell        = -1;
    RDouble mostNegativePressure = PHSPACE::LARGE;

    int zIndex = gridIn->GetZoneID() + 1;
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        for (int m = 0; m < nEquation; ++ m)
        {
            primitiveVariable[m] = q[m][iCell];
        }

        gama = gamma[iCell];
        tm   = t[ITT][iCell];
        tv = tm;
        te = tm;
        if (ntmodel > 1)
        {
            tv = t[ITV][iCell];
            te = tv;
            if (ntmodel == 3)
            {
                te = t[ITE][iCell];
            }
        }
        gas->Primitive2Conservative(primitiveVariable, gama, tv, te, qtry);

        for (int m = 0; m < nLaminar; ++ m)
        {
            qtry[m] += dq[m][iCell];
        }
        for (int m = 0; m < ntmodel-1; ++ m)
        {
            qtry[nLaminar + nChemical + m] += dq[nLaminar + nChemical + m][iCell];
        }

        gas->Conservative2Primitive(qtry, gama, primitiveVariable, tm);

        rm = primitiveVariable[IR];
        pm = primitiveVariable[IP];

        if (rm < densityLimit || pm < pressureLimit)
        {
            if(pm < mostNegativePressure)
            {
                mostNegativePressure = pm;
                mostNegativeCell     = iCell;
            }

            //cout << "Negative rho or pressure found, level = " << grid->GetLevel() << ", cell = " << iCell << endl;
            //TK_Exit::ExceptionExit("Negative!");

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

    if ( nNegativeCell > 0 )
    {
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
                BCRecord **bcr = grid->GetBCRecord();
                int bcType = bcr[faceID]->GetKey();
                cout << "            face bcType = " << bcType << endl;
            }
        }

        if(isInteriorCell)
        {
            cout << " is interior cell ..." << endl;
        }
    }

    delete []primitiveVariable;
    delete []qtry;
}

void NSSolverUnstruct::SetGhostDQLUSGS(Grid *gridIn, RDouble **dq, RDouble **q)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();

    //! Get grid metrics.
    RDouble *xfn  = grid->GetFaceNormalX();
    RDouble *yfn  = grid->GetFaceNormalY();
    RDouble *zfn  = grid->GetFaceNormalZ();
    RDouble *area = grid->GetFaceArea();
    RDouble *vgn  = grid->GetFaceNormalVelocity();

    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    BCRecord **bcr = grid->GetBCRecord();

    RDouble *gama = reinterpret_cast<RDouble *> (grid->GetDataPtr("gama"));
    RDouble **t   = reinterpret_cast<RDouble **> (grid->GetDataPtr("t"));

    Param_NSSolverUnstruct *parameters = GetControlParameters();

    int nLaminar  = parameters->GetLaminarNumber();
	int nChemical = parameters->GetChemicalFlag();
    int ntmodel = parameters->GetTemperatureModel();
    int nEquation = nLaminar + nChemical + ntmodel - 1;

    int isViscous = parameters->IsViscous();

    RDouble *primitiveVariable = new RDouble[nEquation];
    RDouble *Q = new RDouble[nEquation];

    using namespace PHSPACE;
    using namespace GAS_SPACE; 

    RDouble gam1,tv = 0.0,te = 0.0;
    RDouble rm, um, vm, wm, pm, ub, vb, wb;
    RDouble rvn, ru, rv, rw, p_new;

    using namespace IDX;

    int le, re, bcType;
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        le = leftCellofFace[iFace];
        re = rightCellofFace[iFace];

        //! Conservative Variables on level n+1.

        for (int m = 0; m < nEquation; ++ m)
        {
            primitiveVariable[m] = q[m][le];
        }

        tv = t[ITT][le];
        te = tv;
        if (ntmodel > 1)     //! Consider Two-Temperature model and Three-Temperature model.
        {
            tv = t[ITV][le];
            te = tv;
            if (ntmodel == 3)
            {
                te = t[ITE][le];
            }
        }
        gas->Primitive2Conservative(primitiveVariable, gama[le], tv, te, Q);

        for (int m = 0; m < nLaminar; ++ m)
        {
            Q[m] += dq[m][le];
        }
        for (int m = 0; m < ntmodel - 1; m++)
        {
            Q[nLaminar + nChemical + m] += dq[nLaminar + nChemical + m][le];
        }

        //! Non-conservative variables in re in level n.
        rm = q[IR][re];
        um = q[IU][re];
        vm = q[IV][re];
        wm = q[IW][re];
        pm = q[IP][re];

        bcType = bcr[iFace]->GetKey();

        for (int m = 0; m < nLaminar; ++ m)
        {
            dq[m][re] = dq[m][le];
        }
        for (int m = 0; m < ntmodel - 1; m++)
        {
            dq[nLaminar + nChemical + m][re] = dq[nLaminar + nChemical + m][le];
        }

        if (bcType == PHENGLEI::SYMMETRY || (IsWall(bcType) && !isViscous))
        {
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
        else if (IsWall(bcType) && isViscous)
        {
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
        else if (bcType == PHENGLEI::EXTRAPOLATION)
        {
            for (int m = 0; m < nLaminar; ++ m)
            {
                dq[m][re] =  dq[m][le];
            }
            for (int m = 0; m < ntmodel - 1; m++)
            {
                dq[nLaminar + nChemical + m][re] = dq[nLaminar + nChemical + m][le];
            }
        }
        else
        {
            //! FIX_ALL set DQ = 0.
            //! Far_FIELD, PROFILE and FIX_V have not been implemented yet.
            for (int m = 0; m < nLaminar; ++ m)
            {
                dq[m][re] = 0.0;
            }
            for (int m = 0; m < ntmodel - 1; m++)
            {
                dq[nLaminar + nChemical + m][re] = 0.0;
            }
        }
    }

    delete []primitiveVariable;
    delete []Q;
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

void NSSolverUnstruct::SolveLUSGSForward(Grid *gridIn, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal)
{
	UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;
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
    RDouble *spec = reinterpret_cast <RDouble *> (grid->GetDataPtr("spec"));
    RDouble **res = reinterpret_cast <RDouble **> (grid->GetDataPtr("res"));

	Param_NSSolverUnstruct *parameters = GetControlParameters();

    int nLaminar  = parameters->GetLaminarNumber();
	int nChemical = parameters->GetChemicalFlag();
    int nEquation = nLaminar + nChemical;

    int nm = parameters->GetNSEquationNumber();
   
    RDouble refGama = parameters->GetRefGama();
    RDouble refReNumber  = parameters->GetRefReNumber();
    RDouble oRefReNumber = parameters->GetoRefReNumber();

    RDouble csrv = 1.0;   

	int ifLowSpeedPrecon = parameters->GetIfLowSpeedPrecon();

    bool isViscous = parameters->IsViscous();

    RDouble *viscousLaminar, *viscousTurbulent;
    if (isViscous)
    {
        viscousLaminar = reinterpret_cast<RDouble *> (grid->GetDataPtr("visl"));
        viscousTurbulent = reinterpret_cast<RDouble *> (grid->GetDataPtr("vist"));
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
	RDouble *qPrecondition  = new RDouble [nEquation];
	RDouble **MG   = NewPointer2 <RDouble> (nEquation, nEquation);
    RDouble *Ux = new RDouble [nEquation];

	using namespace IDX;

 	//SpectrumRadius(grid,spec);

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

				MXDQ_STD(qNeighbor, nxyz[0], nxyz[1], nxyz[2], nxyz[3], nxyz[4], gama, dqNeighbor, df, nm, nLaminar, rad, -1);
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
				RDouble viscousTurbulentFace = half * (viscousTurbulent[le] + viscousTurbulent[re]);
				RDouble density = half * (q[IDX::IR][le] + q[IDX::IR][re]);

				RDouble viscosity = viscousLaminarFace + viscousTurbulentFace;
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
            //! rhs0: Lx        spec: D.
            dqTry[m] = (dq[m][iCell] - Ux[m] - rhs0[m]) / spec[iCell];

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

void NSSolverUnstruct::SolveLUSGSBackward(Grid *gridIn, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal)
{
	UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;
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
    RDouble *spec = reinterpret_cast <RDouble *> (grid->GetDataPtr("spec"));
    RDouble **res = reinterpret_cast <RDouble **> (grid->GetDataPtr("res"));

	Param_NSSolverUnstruct *parameters = GetControlParameters();

    int nLaminar  = parameters->GetLaminarNumber();
	int nChemical = parameters->GetChemicalFlag();
    int nEquation = nLaminar + nChemical;

    int nm    = parameters->GetNSEquationNumber();
   
    RDouble refGama = parameters->GetRefGama();
    RDouble refReNumber  = parameters->GetRefReNumber();
    RDouble oRefReNumber = parameters->GetoRefReNumber();

    RDouble csrv = 1.0;

	int ifLowSpeedPrecon = parameters->GetIfLowSpeedPrecon();

    int isViscous = parameters->IsViscous();
    RDouble *viscousLaminar, *viscousTurbulent;
    if (isViscous)
    {
        viscousLaminar = reinterpret_cast<RDouble *> (grid->GetDataPtr("visl"));
        viscousTurbulent = reinterpret_cast<RDouble *> (grid->GetDataPtr("vist"));
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

			if ( ifLowSpeedPrecon  == 0  )
		    {
				RDouble gama = gamma[re];

				RDouble rad = GetRad(facePrimitiveVariable,nxyz,gama);

				MXDQ_STD(qNeighbor, nxyz[0], nxyz[1], nxyz[2], nxyz[3], nxyz[4], gama, dqNeighbor, df, nm, nLaminar, rad, -1);
		    }
		    else
		    {
				PreconMatrixFreeDF(qNeighbor, dqNeighbor, nxyz, MG, refGama, df, nLaminar);
		    }	

            if (isViscous)
            {
                RDouble dist;
                dist = ABS(  xfn[face] * ( xcc[re] - xcc[le] )
                           + yfn[face] * ( ycc[re] - ycc[le] )
                           + zfn[face] * ( zcc[re] - zcc[le] ) );

				RDouble vis_l = half * (viscousLaminar[le] + viscousLaminar[re]);
				RDouble vis_t = half * (viscousTurbulent[le] + viscousTurbulent[re]);
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
            //! rhs0: Ux.               spec: D.
            dqTry[m] = (dq[m][iCell] - Lx[m] - rhs0[m]) / spec[iCell];

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

    // If the flow is viscous, it need to count the contribution from viscosity.
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
    }
}

void NSSolverUnstruct::Spectrum(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast( gridIn );
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    int *leftCellofFace = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();

    RDouble *xfn  = grid->GetFaceNormalX();
    RDouble *yfn  = grid->GetFaceNormalY();
    RDouble *zfn  = grid->GetFaceNormalZ();
    RDouble *area = grid->GetFaceArea();
	RDouble *vgn  = grid->GetFaceNormalVelocity();
	RDouble *vol  = grid->GetCellVolume();

    RDouble *dt   = reinterpret_cast<RDouble *> (grid->GetDataPtr("dt"));
    RDouble *srad = reinterpret_cast<RDouble *> (grid->GetDataPtr("srad"));
    RDouble *srv  = reinterpret_cast<RDouble *> (grid->GetDataPtr("srv"));
    RDouble *spec = reinterpret_cast <RDouble *> (grid->GetDataPtr("spec"));

	using namespace IDX;

	// If flow is unsteady, need to count the contribution of unsteady.
	Param_NSSolverUnstruct *parameters = GetControlParameters();
	int isUnsteady = parameters->GetIsUnsteady();
	if (!isUnsteady)
	{
		for (int iCell = 0; iCell < nTotal; ++ iCell)
		{
			spec[iCell] = 1.0 / dt[iCell]; //! notice: the dt[iCell] has divided vol, so the spec[iCell] has multiplied the vol by default.
		}
	}
	else
	{
        RDouble dualTimeCoefficient[7];
        const int methodOfDualTime = 3;
        ComputeDualTimeCoefficient(methodOfDualTime, dualTimeCoefficient);
        RDouble dualTimeSpectrumC1 = - dualTimeCoefficient[3];
        RDouble dualTimeSpectrumC2 =   dualTimeCoefficient[6];

		for (int iCell = 0; iCell < nTotal; ++ iCell)
		{
            spec[iCell] = dualTimeSpectrumC1 * vol[iCell] + (dualTimeSpectrumC2 / dt[iCell]); //! notice: the dt[iCell] has divided vol, as same as above.
		}
	}

	for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        spec[iCell] += srad[iCell];
    }

    // If the flow is viscous, it need to count the contribution from viscosity.
    int isViscous = parameters->IsViscous();

    if (isViscous)
    {
        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            spec[iCell] += srv[iCell];
        }
    }

	int nm        = parameters->GetNSEquationNumber();
    int nLaminar  = parameters->GetLaminarNumber();
	int nChemical = parameters->GetChemicalFlag();
    int nChemicalRadius = parameters->GetNChemicalRadius();

    RDouble rad_chem;

    if (nChemical == 1 && nChemicalRadius == 1)
	{
        RDouble **srs = reinterpret_cast<RDouble **> (grid->GetDataPtr("srs"));
        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
			rad_chem = zero;
            for (int m = nm; m < nLaminar; ++ m)
			{
                rad_chem = MAX(nChemicalRadius * srs[ m ][ iCell ], rad_chem);
			}
            spec[iCell] += rad_chem; //! Referring to the structure part, the dt[iCell] contained in the spec[iCell] can be improved for more stability here.

        }
	}
    //grid->SetGhostCell(spec, nLaminar);
}

void NSSolverUnstruct::RestrictDQ(RDouble *dqtry, int nm, RDouble lim_max)
{
	for ( int m = 0; m < nm; ++ m )
	{
		dqtry[m] = MAX(MIN(dqtry[m], lim_max),-lim_max);
	}
}

void NSSolverUnstruct::AirForceCoef(ActionKey *actkey)
{
    const int ALL_PART = -1;
    AirForceCoefParts(actkey, ALL_PART);

    int Coordinate = LocalCoordinate;
    uint_t nTBC = GlobalBoundaryCondition::GetNumberofBC();
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
    //fstream file;

    UnstructGrid * grid = UnstructGridCast(GetGrid(level));

    RDouble **q   = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("q") );
    RDouble *visl = reinterpret_cast< RDouble * > ( grid->GetDataPtr("visl") );
    Param_NSSolverUnstruct * parameters = GetControlParameters();
    RDouble oRefReNumber = parameters->GetoRefReNumber();

    int viscousType = parameters->GetViscousType();

    const int ALL_PART = -1;
    BoundaryCondition * boundaryCondition = 0;
    string partBCName = "";

    RDouble TorqueRefX, TorqueRefY, TorqueRefZ;

    //! Used to compute hinge moment:
    int dumpHingeMoment = 0;
    RDouble localCoordAxis0[3];
    RDouble localCoordAxis1[3];
    RDouble hingeMoment = 0.0;

    if(partID == ALL_PART)
    {
        TorqueRefX = parameters->GetTorqueRefX();
        TorqueRefY = parameters->GetTorqueRefY();
        TorqueRefZ = parameters->GetTorqueRefZ();
    }
    else if ((partID != ALL_PART) && (Coordinate == GlobalCoordinate))
    {
        boundaryCondition = GlobalBoundaryCondition::GetBC(partID);
        partBCName = boundaryCondition->GetBCName();

        TorqueRefX = parameters->GetTorqueRefX();
        TorqueRefY = parameters->GetTorqueRefY();
        TorqueRefZ = parameters->GetTorqueRefZ();
    }
    else if ((partID != ALL_PART) && (Coordinate == LocalCoordinate))
    {
        boundaryCondition = GlobalBoundaryCondition::GetBC(partID);
        boundaryCondition->GetData("TorqueRefX", &TorqueRefX, PHDOUBLE, 1);
        boundaryCondition->GetData("TorqueRefY", &TorqueRefY, PHDOUBLE, 1);
        boundaryCondition->GetData("TorqueRefZ", &TorqueRefZ, PHDOUBLE, 1);
        partBCName = boundaryCondition->GetBCName();
        
        boundaryCondition->GetData("dumpHingeMoment", &dumpHingeMoment, PHINT, 1);
        boundaryCondition->GetData("localCoordAxis0", localCoordAxis0, PHDOUBLE, 3);
        boundaryCondition->GetData("localCoordAxis1", localCoordAxis1, PHDOUBLE, 3);
    }

    int vis_run = 0;
    if (viscousType != INVISCID)
    {
        vis_run = 1;
        //vis_run = 0;
    }

    using namespace IDX;
    RDouble *primitiveVarFarfield = parameters->GetPrimitiveVarFarfield();
    RDouble poo = primitiveVarFarfield[IP];

    RDouble *xfn  = grid->GetFaceNormalX();
    RDouble *yfn  = grid->GetFaceNormalY();
    RDouble *zfn  = grid->GetFaceNormalZ();
    RDouble *area = grid->GetFaceArea() ;
    RDouble *vol  = grid->GetCellVolume();

    RDouble *xfc = grid->GetFaceCenterX();
    RDouble *yfc = grid->GetFaceCenterY();
    RDouble *zfc = grid->GetFaceCenterZ();
    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    int * leftCellofFace = grid->GetLeftCellOfFace();
    int * rightCellofFace = grid->GetRightCellOfFace();

    BCRecord **bcr = grid->GetBCRecord();

    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

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

    RDouble **gradPrimtiveVarFaceX = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarFaceX"));
    RDouble **gradPrimtiveVarFaceY = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarFaceY"));
    RDouble **gradPrimtiveVarFaceZ = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarFaceZ"));

    using namespace IDX;

    int nsolid_surface = grid->GetNumberOfWallCell();

    // get the BCRegion information of unstruct grid.
    BCRegionUnstruct **bcRegionUnstruct = grid->GetBCRegionUnstruct();
    int nBCRegionUnstruct = grid->GetNBCRegionUnstruct();

    // cycle all bcregions and find the solid surface type of boundary.
    for (int iBCRegionUnstruct = 0; iBCRegionUnstruct < nBCRegionUnstruct; iBCRegionUnstruct++)
    {
        int bcType = bcRegionUnstruct[iBCRegionUnstruct]->GetBCType();
        string bcName = bcRegionUnstruct[iBCRegionUnstruct]->GetBCName();

        // if the bc type is wall, calculate the wall coefficient.
        if(bcType == PHENGLEI::SOLID_SURFACE)
        {
            if (partID != ALL_PART && bcName != partBCName)
            {
                continue;
            }
            vector<int> *faceIndex = bcRegionUnstruct[iBCRegionUnstruct]->GetFaceIndex();
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

                //pressure drag
                //pw    = q[IP][le] + dpdx[iFace] * dx + dpdy[iFace] * dy + dpdz[iFace] * dz;
                //if ( pw < 0.0 ) pw = q[IP][le];
                pw = q[IP][le];
                cp    = two * ( pw - poo );

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

                if ( viscousType > 0 )
                {
                    dudx = gradPrimtiveVarFaceX[IU][le];
                    dudy = gradPrimtiveVarFaceY[IU][le];
                    dudz = gradPrimtiveVarFaceZ[IU][le];

                    dvdx = gradPrimtiveVarFaceX[IV][le];
                    dvdy = gradPrimtiveVarFaceY[IV][le];
                    dvdz = gradPrimtiveVarFaceZ[IV][le];

                    dwdx = gradPrimtiveVarFaceX[IW][le];
                    dwdy = gradPrimtiveVarFaceY[IW][le];
                    dwdz = gradPrimtiveVarFaceZ[IW][le];

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
                    txx = vis * c23 * ( 2.0 * dudx - dvdy - dwdz );
                    tyy = vis * c23 * ( 2.0 * dvdy - dwdz - dudx );
                    tzz = vis * c23 * ( 2.0 * dwdz - dudx - dvdy );
                    txy = vis * ( dudy + dvdx );
                    txz = vis * ( dudz + dwdx );
                    tyz = vis * ( dvdz + dwdy );
                    tyx = txy;
                    tzx = txz;
                    tzy = tyz;

                    fvsx = - two * ( nx * txx + ny * tyx + nz * tzx ) * oRefReNumber;
                    fvsy = - two * ( nx * txy + ny * tyy + nz * tzy ) * oRefReNumber;
                    fvsz = - two * ( nx * txz + ny * tyz + nz * tzz ) * oRefReNumber;

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

                if( dumpHingeMoment )
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
    RDouble **q   = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("q") );
    RDouble *visl = reinterpret_cast< RDouble * > ( grid->GetDataPtr("visl") );

    RDouble *primitiveVarFarfield = parameters->GetPrimitiveVarFarfield();
    RDouble roo = primitiveVarFarfield[IR];
    RDouble uoo = primitiveVarFarfield[IU];
    RDouble voo = primitiveVarFarfield[IV];
    RDouble woo = primitiveVarFarfield[IW];
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
    TorqueRefX = parameters->GetTorqueRefX();
    TorqueRefY = parameters->GetTorqueRefY();
    TorqueRefZ = parameters->GetTorqueRefZ();

    RDouble *x = grid->GetX();
    RDouble *y = grid->GetY();
    RDouble *z = grid->GetZ();

    RDouble *xfn  = grid->GetFaceNormalX();
    RDouble *yfn  = grid->GetFaceNormalY();
    RDouble *zfn  = grid->GetFaceNormalZ();
    RDouble *area = grid->GetFaceArea() ;
    RDouble *vol  = grid->GetCellVolume();

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
    BCRecord **bcr = grid->GetBCRecord();

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

    RDouble AoA, attack, angleSlide, angleSlided;
    AoA = parameters->GetAoA();
    attack   = AoA * PI / 180.0;
    angleSlided = parameters->GetAngleOfSlide();
    angleSlide = angleSlided * PI / 180.0;

    RDouble sina, cosa, sinb, cosb;
    sina = sin(attack);
    cosa = cos(attack);
    sinb = sin(angleSlide);
    cosb = cos(angleSlide);

    RDouble **gradPrimtiveVarFaceX = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("gradPrimtiveVarFaceX"  ) );
    RDouble **gradPrimtiveVarFaceY = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("gradPrimtiveVarFaceY"  ) );
    RDouble **gradPrimtiveVarFaceZ = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("gradPrimtiveVarFaceZ"  ) );

    RDouble pw, cp, pw_Dim;
    RDouble dfx, dfy, dfz, dpx, dpy, dpz, dmx, dmy, dmz, fvsx, fvsy, fvsz;
    RDouble nx, ny, nz;
    RDouble normalComponetX, normalComponetY, normalComponetZ;
    RDouble vis, txx, txy, txz, tyx, tyy, tyz, tzx, tzy, tzz;
    RDouble dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz;
    RDouble xc, yc, zc, dx, dy, dz, ods;

    RDouble *wallDistance = grid->GetWallDist();

    RDouble *faceCP    = new RDouble[nBoundFace];
    RDouble *faceCF    = new RDouble[nBoundFace];
    RDouble *faceQ     = new RDouble[nBoundFace];
    RDouble *faceYPlus = new RDouble[nBoundFace];
    RDouble *facePW    = new RDouble[nBoundFace];

    //! Compute cp and cf on wall face.
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int bcType = bcr[iFace]->GetKey();
        if (IsWall(bcType))
        {
            int le,re;
            le = leftCellofFace [ iFace ];
            re = rightCellofFace[ iFace ];

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

            if (viscousType > 0)
            {
                dudx = gradPrimtiveVarFaceX[IU][le];
                dudy = gradPrimtiveVarFaceY[IU][le];
                dudz = gradPrimtiveVarFaceZ[IU][le];

                dvdx = gradPrimtiveVarFaceX[IV][le];
                dvdy = gradPrimtiveVarFaceY[IV][le];
                dvdz = gradPrimtiveVarFaceZ[IV][le];

                dwdx = gradPrimtiveVarFaceX[IW][le];
                dwdy = gradPrimtiveVarFaceY[IW][le];
                dwdz = gradPrimtiveVarFaceZ[IW][le];

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

                txx = vis * two3rd * ( 2.0 * dudx - dvdy - dwdz );
                tyy = vis * two3rd * ( 2.0 * dvdy - dwdz - dudx );
                tzz = vis * two3rd * ( 2.0 * dwdz - dudx - dvdy );
                txy = vis * ( dudy + dvdx );
                txz = vis * ( dudz + dwdx );
                tyz = vis * ( dvdz + dwdy );
                tyx = txy;
                tzx = txz;
                tzy = tyz;

                fvsx = - two * ( nx * txx + ny * tyx + nz * tzx ) * oRefReNumber;
                fvsy = - two * ( nx * txy + ny * tyy + nz * tzy ) * oRefReNumber;
                fvsz = - two * ( nx * txz + ny * tyz + nz * tzz ) * oRefReNumber;

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
            if (viscousType > 0)
            {
                yPlus = refReNumber * sqrt(taoWall * density) * wallDistance[le] / vis;
            }

            faceCP   [iFace] = cp;
            faceCF   [iFace] = cf;
            faceYPlus[iFace] = yPlus;
            facePW   [iFace] = pw_Dim;
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
        dimensionalQ = refDimensionalDensity * velocityInflow * velocityInflow * velocityInflow * 0.001;		//! 0.001: w --> kw

        RDouble prandtlLaminar = parameters->GetPrandtlLaminar();
        RDouble nonDimensionalSutherlandTemperature = parameters->GetNonDimensionalSutherlandTemperature();

        RDouble Tw  = wallTemperature / refDimensionalTemperature;
        RDouble Cv  = 1.0 / refGama / (refGama - 1.0) / refMachNumber / refMachNumber;
        RDouble tmu = pow(Tw, 1.5) * ( 1.0 + nonDimensionalSutherlandTemperature ) / ( Tw + nonDimensionalSutherlandTemperature );

        RDouble **t = reinterpret_cast< RDouble ** > (grid->GetDataPtr("t"));

        for (int iFace = 0; iFace < nBoundFace; ++ iFace)
        {
            int bcType = bcr[iFace]->GetKey();
            if (IsWall( bcType ))
            {
                int le = leftCellofFace [iFace];
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

    //! Communicate data on interface.
    int *lableCell = new int[nTotalCell];
    std::fill_n(lableCell, nTotalCell, 0);

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int bcType = bcr[iFace]->GetKey();
        if (IsWall(bcType))
        {
            int cL = leftCellofFace[iFace];
            lableCell[cL] = 1;
        }
    }

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int bcType = bcr[iFace]->GetKey();
        if (bcType == INTERFACE)
        {
            int cL = leftCellofFace[iFace];
            if (lableCell[cL] == 1)
            {
                lableCell[cL] = 2;
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

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int bcType = bcr[iFace]->GetKey();
        if (IsWall(bcType))
        {
            int cL = leftCellofFace[iFace];
            if(lableCell[cL] == 2)
            {
                cellCF[cL] = faceCF[iFace];
                cellCP[cL] = faceCP[iFace];
                if (wallTemperature >= 0)
                {
                    cellQ[cL] = faceQ[iFace];
                }
            }
        }
    }
    delete []lableCell;

    PHSPACE::CommunicateInterfaceValue(grid, cellCP2D, "cellCP", 1);
    PHSPACE::CommunicateInterfaceValue(grid, cellCF2D, "cellCF", 1);
    if (wallTemperature >= 0)
    {
        PHSPACE::CommunicateInterfaceValue(grid, cellQ2D, "cellQ", 1);
    }

    //! Compute data on wall node.
    int *nCell = new int[nTotalNode];
    std::fill_n(nCell, nTotalNode, 0);

    RDouble *nodeQ = new RDouble[nTotalNode];
    std::fill_n(nodeQ, nTotalNode, 0.0);

    RDouble *nodeCP = new RDouble[nTotalNode];
    std::fill_n(nodeCP, nTotalNode, 0.0);

    RDouble *nodeCF = new RDouble[nTotalNode];
    std::fill_n(nodeCF, nTotalNode, 0.0);

    RDouble *nodeYPlus = new RDouble[nTotalNode];
    std::fill_n(nodeYPlus, nTotalNode, 0.0);

    RDouble *nodePW = new RDouble[nTotalNode];
    std::fill_n(nodePW, nTotalNode, 0.0);

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int bcType = bcr[iFace]->GetKey();
        if (IsWall(bcType))
        {
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

    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        if (nCell[iNode])
        {
            nodeYPlus[iNode] /= nCell[iNode];
            nodePW   [iNode] /= nCell[iNode];
        }
    }

    //! For interface.
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        int bcType = bcr[iFace]->GetKey();
        if (bcType == INTERFACE)
        {
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

    delete []faceCP;
    delete []faceCF;
    delete []faceQ;
    delete []nCell;
    delete []faceYPlus;
    delete []facePW;
    DelPointer2(cellQ2D);
    DelPointer2(cellCP2D);
    DelPointer2(cellCF2D);

    //! Compute data on wall face center.
    RDouble *cpFaceCenter    = new RDouble[nBoundFace];
    RDouble *cfFaceCenter    = new RDouble[nBoundFace];
    RDouble *qFaceCenter     = new RDouble[nBoundFace];
    RDouble *yPlusFaceCenter = new RDouble[nBoundFace];
    RDouble *PWFaceCenter    = new RDouble[nBoundFace];

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
    int *cell2node        = new int[TotalNumFaceNodes];

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

    RDouble *coorX     = new RDouble[NumPts];
    RDouble *coorY     = new RDouble[NumPts];
    RDouble *coorZ     = new RDouble[NumPts];
    RDouble *dumpCP    = new RDouble[NumPts];
    RDouble *dumpCF    = new RDouble[NumPts];
    RDouble *dumpQ     = new RDouble[NumPts];
    RDouble *dumpQNon  = new RDouble[NumPts];
    RDouble *dumpYPlus = new RDouble[NumPts];
    RDouble *dumpPW    = new RDouble[NumPts];

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

            dumpCP   [iNode] = cpFaceCenter   [nodeIndex];
            dumpCF   [iNode] = cfFaceCenter   [nodeIndex];
            dumpYPlus[iNode] = yPlusFaceCenter[nodeIndex];
            dumpPW   [iNode] = PWFaceCenter   [nodeIndex];
            dumpQ    [iNode] = qFaceCenter    [nodeIndex] * dimensionalQ;
            dumpQNon [iNode] = qFaceCenter    [nodeIndex];
        }
    }

    delete []nodeQ;
    delete []nodeCP;
    delete []nodeCF;
    delete []nodeYPlus;
    delete []nodePW;
    delete []cpFaceCenter;
    delete []cfFaceCenter;
    delete []qFaceCenter;
    delete []yPlusFaceCenter;
    delete []PWFaceCenter;

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

    delete []coorX;
    delete []coorY;
    delete []coorZ;
    delete []dumpCP;
    delete []dumpCF;
    delete []dumpQ;
    delete []dumpQNon;
    delete []dumpYPlus;
    delete []dumpPW;
    delete []cell2node;
}

// Dump out some flow field variables, for the case of turbulent flat plate of NASA Turbulent Center.
// Inflow condition: Ma = 0.2, Re = 5.0e6.
void NSSolverUnstruct::Turbulent_Flat_Plate_Output_NASA(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalNode = grid->GetNTotalNode();
    RDouble ** qVist;

    //100.000 kPa = 1 bar ≈ 14.504 psi
    //condition for calculation
    //Static Pressure (psia) 6.0 = 100000 / 14.504 * 6 Pa = 41367.8986 Pa
    //Temperature(R) 700.0 = (5/9)*700 = 388.88K

    //////////////////////////////////////////////////////////////////////////
    // Dump out the Velocity
    //////////////////////////////////////////////////////////////////////////
    vector<string> title_tecplot;
    title_tecplot.push_back("title=\"Flow Fields of PHengLEI\"");
    title_tecplot.push_back("variables=");
    title_tecplot.push_back("\"y\"");
    title_tecplot.push_back("\"u\"");
    title_tecplot.push_back("\"y+\"");
    title_tecplot.push_back("\"u+\"");    
    int nvarplot;
    nvarplot = 2;

    Param_NSSolverUnstruct * parameters = GetControlParameters();
    int nm = parameters->GetNSEquationNumber();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nEquation = nLaminar + nChemical;

    using namespace PHSPACE;
    RDouble refGama = parameters->GetRefGama();
    RDouble reference_mach_number = parameters->GetRefMachNumber();
    RDouble twall = parameters->GetWallTemperature();
    RDouble refDimensionalTemperature = parameters->GetRefDimensionalTemperature();
    RDouble refReNumber = 1.0e5;
    refReNumber = parameters->GetRefReNumber();
    RDouble tsuth;
    GlobalDataBase::GetData("tsuth",&tsuth,PHDOUBLE,1);

    fstream file;
    PHSPACE::OpenFile(file, "results/turblent_flat_plate.dat", ios_base::out);

    int wordwidth = 20;
    RDouble prandtlLaminar    = parameters->GetPrandtlLaminar();

    RDouble * visl = reinterpret_cast< RDouble * > ( grid->GetDataPtr("visl") );
    RDouble * vist = reinterpret_cast< RDouble * > ( grid->GetDataPtr("vist") );

    int nline = 2;

    vector<RDouble> xPositionList(nline);
    xPositionList[0] = 0.97008;
    xPositionList[1] = 1.90334;    

    for ( std::size_t i = 0; i < title_tecplot.size(); ++ i )
    {
        file << title_tecplot[i] << "\n";
    }

    RDouble ** nodevar = NewPointer2<RDouble>(nEquation + 1, nTotalNode);
    CreateNodeVariable(grid, nodevar);

    for ( int m = 0; m < nline; ++ m )
    {
        PHCutPlane * cutpl = new PHCutPlane();
        CutPlane(grid, nodevar, nEquation, xPositionList[m], cutpl, X_DIR);

        cutpl->SortData(cutpl->y);

        int nT = cutpl->nTotalNode;
        file << "zone  T = \"x = " << xPositionList[m] << "\" \n";
        //找到壁面第一层的点
        vector<RDouble> ylist(nT);
        //RDouble tmp[10000],tmp1[10000];
        for ( int i = 0; i < nT; ++ i )
        {
            ylist[i] = cutpl->y[i];
            //tmp[i] = ylist[i];
        }
        sort(ylist.begin(),ylist.end());

        for ( int i = 0; i < nT; ++ i )
        {
            //tmp1[i] = ylist[i];
        }

        int yindex = 1;
        for ( int i = 0; i < nT; ++ i )
        {
            if ( ylist[i] )
            {
                yindex = i;
                break;
            }
        }

        for ( int i = 0; i < nT; ++ i )
        {
            if ( ylist[yindex] == cutpl->y[i] )
            {
                yindex = i;
                break;
            }
        }

        for ( int i = 0; i < nT; ++ i )
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
            mu   = tm * sqrt( tm ) * ( 1.0 + tsuth ) / ( tm + tsuth );

            rw   = cutpl->data[0][yindex];
            uw   = cutpl->data[1][yindex];
            vw   = cutpl->data[2][yindex];
            ww   = cutpl->data[3][yindex];
            pw   = cutpl->data[4][yindex];
            tw   = refGama * reference_mach_number * reference_mach_number * pm / rm;

            muw  = tm * sqrt( tm ) * ( 1.0 + tsuth ) / ( tm + tsuth );
            yw   = cutpl->y[yindex];

            dudy = uw / yw;
            tauw = muw * dudy;

            utau = sqrt(tauw/(rw * refReNumber));

            //注意这里的定义
            up   = um / utau;
            yp   = utau * ym * refReNumber / ( mu / rm );

            //up   = um / sqrt( tauw / ( rm * reynolds) );
            //yp   = sqrt(tauw * reynolds * rm ) * ym / mu;

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

    //////////////////////////////////////////////////////////////////////////
    // Dump out the friction
    //////////////////////////////////////////////////////////////////////////
    RDouble **q   = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("q"  ) );
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
        for ( int i = 0; i < nT; ++ i )
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
            mu   = tm * sqrt( tm ) * ( 1.0 + tsuth ) / ( tm + tsuth );
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

    //////////////////////////////////////////////////////////////////////////
    // Dump out the turbulent viscosity of x = 0.97
    //////////////////////////////////////////////////////////////////////////
    qVist = new RDouble *[1];				// Tempal 2D array for vist, to dump both x and y value.
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
        for ( int i = 0; i < nT; ++ i )
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

    //////////////////////////////////////////////////////////////////////////
    // Dump out the the MAXIMUM eddy viscosity as a function of x
    //////////////////////////////////////////////////////////////////////////
    using namespace PHMPI;
    if(GetNumberofGlobalZones() > 1)
    {
        return;
    }
    string maxVistOfX = "results/turblent_flat_plate_nasa_vist_maximum.dat";
    PHSPACE::OpenFile(file, maxVistOfX.c_str(), ios_base::out);

    file << "title=\"the maximum turbulent viscosity\"\n";
    file << "variables = \"X\"    \"vist\"\n";

    qVist = new RDouble *[1];				// Tempal 2D array for vist, to dump both x and y value.
    qVist[0] = vist;

    // Get the number of x stand position.
    int nXStand = 0;
    int nBoundFace  = grid->GetNBoundFace();
    RDouble* xfc = grid->GetFaceCenterX();
    BCRecord ** bcr = grid->GetBCRecord();	
    for ( int iFace = 0; iFace < nBoundFace; ++ iFace )
    {
        int type = bcr[ iFace ]->GetKey();
        if ( type == PHENGLEI::SOLID_SURFACE || type == PHENGLEI::SYMMETRY)
        {
            ++ nXStand;
        }
    }

    // Get the x stand position.
    RDouble* xOfStand = new RDouble [nXStand];
    nXStand = 0;
    for ( int iFace = 0; iFace < nBoundFace; ++ iFace )
    {
        int type = bcr[ iFace ]->GetKey();
        if ( type == PHENGLEI::SOLID_SURFACE || type == PHENGLEI::SYMMETRY)
        {
            xOfStand[nXStand++] = xfc[iFace];
        }
    }

    vector< vector < RDouble > > var_list;
    vector < RDouble > x_list;
    vector < RDouble > var(1);
    cutpl = new PHCutPlane();				// Store the maximum viscosity of each xcut stand.

    RDouble* maxVistOfStandX = new RDouble [nXStand];
    for (int iStand = 0; iStand < nXStand; ++ iStand)
    {
        RDouble xCut = xOfStand[iStand];
        PHCutPlane* cutpX = new PHCutPlane();	// Store the viscosity of ONE x stand, x = xcut.

        RDouble maxVist = -LARGE;
        CutPlaneCellCenter(grid,qVist,1,xCut,cutpX,X_DIR);
        nT = cutpX->nTotalNode;
        for ( int i = 0; i < nT; ++ i )
        {			
            RDouble vist = cutpX->data[0][i];
            maxVist = MAX(maxVist, vist);
        }

        var[0] = maxVist;
        var_list.push_back(var);
        x_list.push_back( xCut );

        delete cutpX;
    }

    cutpl->neqn   = 1;
    cutpl->nTotalNode = static_cast<int>(var_list.size());
    cutpl->data   = new RDouble * [ 1 ];
    for ( int m = 0; m < 1; ++ m )
    {
        cutpl->data[m] = new RDouble [ cutpl->nTotalNode ];
    }
    cutpl->x = new RDouble [ cutpl->nTotalNode ];
    cutpl->y = new RDouble [ cutpl->nTotalNode ];
    cutpl->z = new RDouble [ cutpl->nTotalNode ];
    for ( int inode = 0; inode < cutpl->nTotalNode; ++ inode )
    {
        cutpl->x[inode] = x_list[inode];
        cutpl->y[inode] = 0;
        cutpl->z[inode] = 0;
    }
    for ( int inode = 0; inode < cutpl->nTotalNode; ++ inode )
    {
        for ( int m = 0; m < 1; ++ m )
        {
            cutpl->data[m][inode] = var_list[inode][m];
        }
    }

    cutpl->SortData(cutpl->x);

    nT = cutpl->nTotalNode;
    file << "zone  i = " << nT  << " \n";
    for ( int i = 0; i < nT; ++ i )
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

    RDouble **q = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("q") );

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nEquation = nLaminar + nChemical;

    ComputeNodeVariable(grid, nodevar, q, nEquation);
}

void NSSolverUnstruct::ComputeAverageVisualNodeField(Grid *gridIn, RDouble **qn)
{
	UnstructGrid *grid = UnstructGridCast(gridIn);

    Param_NSSolverUnstruct *parameters = GetControlParameters();
	int nLaminar = parameters->GetLaminarNumber();

    RDouble **qAverage = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("qAverage") );

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
    if ( sliceAxis == X_DIR )
    {
        xyz = x;
    }
    else if ( sliceAxis == Y_DIR )
    {
        xyz = y;
    }
    else
    {
        xyz = z;
    }

    RDouble Small	= 1.0e-10;
    int nCount		= 0;
    BCRecord **bcr	= grid->GetBCRecord();
    for ( int iFace = 0; iFace < nBoundFace; ++ iFace )
    {
        int bcType = bcr[ iFace ]->GetKey();

        if ( bcType != PHENGLEI::SOLID_SURFACE )
        {
            nCount += node_number_of_each_face[ iFace ];
            continue;
        }

        RDouble coor_min =   LARGE;
        RDouble coor_max = - LARGE;

        for ( int iNode = 0; iNode < node_number_of_each_face[ iFace ]; ++ iNode )
        {
            int index = face2node[ nCount + iNode ];
            coor_min = MIN( xyz[ index ],coor_min );
            coor_max = MAX( xyz[ index ],coor_max );
        }

        if ( (coor_min - slicePostion) * (coor_max - slicePostion) > 0.0 )
        {
            nCount += node_number_of_each_face[ iFace ];
            continue;
        }

        for ( int iNode = 0; iNode < node_number_of_each_face[ iFace ]; ++ iNode )
        {
            int Node_index0 = face2node[ nCount + iNode % node_number_of_each_face[ iFace ] ];
            int Node_index1 = face2node[ nCount + (iNode + 1) % node_number_of_each_face[ iFace ] ];

            if ( (xyz[ Node_index0 ] - slicePostion) * (xyz[ Node_index1 ] - slicePostion) > 0.0 ) continue;

            if ( fabs( xyz[ Node_index1 ] - xyz[ Node_index0 ] ) < Small )
            {
                xSlice.push_back(x[ Node_index0 ]);
                xSlice.push_back(x[ Node_index1 ]);

                ySlice.push_back(x[ Node_index0 ]);
                ySlice.push_back(x[ Node_index1 ]);

                zSlice.push_back(x[ Node_index0 ]);
                zSlice.push_back(x[ Node_index1 ]);

                continue;
            }

            if ( sliceAxis == X_DIR )
            {
                xSlice.push_back(slicePostion);
            }
            else if ( fabs( x[ Node_index0 ] - x[ Node_index1 ] ) < Small)
            {
                xSlice.push_back(x[ Node_index0] );
            }
            else
            {
                RDouble Node_x = (slicePostion - xyz[ Node_index0] ) * ( x[ Node_index1 ] - x[ Node_index0 ] ) / ( xyz[ Node_index1 ] - xyz[ Node_index0 ] ) + x[ Node_index0 ];
                xSlice.push_back(Node_x);
            }

            if ( sliceAxis == Y_DIR )
            {
                ySlice.push_back(slicePostion);
            }
            else if ( fabs( y[ Node_index0 ] - y[ Node_index1 ] ) < Small)
            {
                ySlice.push_back(y[ Node_index0 ]);
            }
            else
            {
                RDouble Node_y = (slicePostion - xyz[ Node_index0 ]) * ( y[ Node_index1 ] - y[ Node_index0 ] ) / ( xyz[ Node_index1 ] - xyz[ Node_index0 ] ) + y[ Node_index0 ];
                ySlice.push_back(Node_y);
            }

            if ( sliceAxis == Z_DIR )
            {
                zSlice.push_back(slicePostion);
            }
            else if ( fabs( z[ Node_index0 ] - z[ Node_index1 ] ) < Small)
            {
                zSlice.push_back(z[ Node_index0 ]);
            }
            else
            {
                RDouble Node_z = (slicePostion - xyz[ Node_index0 ]) * ( z[ Node_index1 ] - z[ Node_index0 ] ) / ( xyz[ Node_index1 ] - xyz[ Node_index0 ] ) + z[ Node_index0 ];
                zSlice.push_back(Node_z);
            }

            for ( int iVarplot = 0; iVarplot < nvarplot; ++ iVarplot )
            {
                if ( fabs( qn[ iVarplot ][ Node_index0 ] - qn[ iVarplot ][ Node_index1 ] ) < Small)
                {
                    RDouble Node_q = qn[ iVarplot ][ Node_index0 ];
                    qnSlice[ iVarplot ].push_back(Node_q);
                }
                else
                {
                    RDouble Node_q = (slicePostion - xyz[ Node_index0 ]) * ( qn[ iVarplot ][ Node_index1 ] - qn[ iVarplot ][ Node_index0 ] ) / ( xyz[ Node_index1 ] - xyz[ Node_index0 ] ) + qn[ iVarplot ][ Node_index0 ];
                    qnSlice[ iVarplot ].push_back(Node_q);
                }
            }
        }

        nCount += node_number_of_each_face[ iFace ];
    }

    size_t NumPts = xSlice.size();
    PHWrite( cdata, &NumPts, 1 );

    if ( NumPts == 0 )  return;

    PHWrite( cdata, &xSlice[ 0 ], NumPts );
    PHWrite( cdata, &ySlice[ 0 ], NumPts );
    PHWrite( cdata, &zSlice[ 0 ], NumPts );
    for ( int iVarplot = 0; iVarplot < nvarplot; ++ iVarplot )
    {
        PHWrite( cdata, &qnSlice[ iVarplot ][ 0 ], NumPts );
    }
}

void NSSolverUnstruct::BoundaryVisualization(Grid *gridIn, ActionKey *actkey, RDouble **qn, int nvarplot)
{
	UnstructGrid *grid = UnstructGridCast(gridIn);
	DataContainer * cdata = actkey->GetTecData();

    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();

    int * node_number_of_each_face = grid->GetNodeNumberOfEachFace();
    int * leftCellofFace = grid->GetLeftCellOfFace();
    int * rightCellofFace = grid->GetRightCellOfFace();
    int * face2node = grid->GetFace2Node();

    RDouble * x = grid->GetX();
    RDouble * y = grid->GetY();
    RDouble * z = grid->GetZ();

    using namespace PHENGLEI;
    vector < vector < int > > face2nodelist(nBoundFace);
    vector < int > linkmap;
    set<int> bclist;

	set< pair<int, string> > bcnameMap;

	BCRecord **bcr = grid->GetBCRecord();
	for ( int iFace = 0; iFace < nBoundFace; ++ iFace )
	{
	    int bcType = bcr[iFace]->GetKey();
	    bclist.insert(bcType);
		string bcname = bcr[iFace]->GetBCName();

		bcnameMap.insert(pair<int,string>(bcType, bcname));
	}
	
	int TecioMission;
	set< pair<int, string> >::iterator iter;
    for ( iter = bcnameMap.begin(); iter != bcnameMap.end(); ++ iter )
	{
		if ( (*iter).first == PHENGLEI::INTERFACE || (*iter).first < 0  ) continue;

        face2nodelist.resize(0);
        linkmap.resize(0);

        GetFace2NodeList(grid, *iter, linkmap, face2nodelist);

		if( linkmap.size() == 0 )
		{
			continue;
		}

		TecioMission = TEC_SPACE::WriteBoundary;
		PHWrite(cdata, &TecioMission, 1);

		int BCtype = (*iter).first;
		PHWrite(cdata, &BCtype, 1);

        /*std::ostringstream os_zone;
        os_zone << "\"" << grid->GetZoneID() << "BC=" << *iter << "\" ";
        string zone_title = os_zone.str();*/

		size_t NumPts		= linkmap.size();
		size_t NumElements = face2nodelist.size();
		size_t NumFaces	= NumElements * 4;
		PHWrite(cdata, &NumPts, 1);
		PHWrite(cdata, &NumElements, 1);
		PHWrite(cdata, &NumFaces, 1);

        size_t TotalNumFaceNodes_Rect     = NumElements * 4;
		int * cell2node				   = new int [TotalNumFaceNodes_Rect];
		int * node_number_of_each_cell = new int [NumElements];

		int count = 0;
		for ( std::size_t i = 0; i < face2nodelist.size(); ++ i )
        {
			node_number_of_each_cell[i] = 4;
            int np     = static_cast<int>(face2nodelist[i].size());
            int index0 = face2nodelist[i][0];
            for ( int ip = 0; ip < np; ++ ip )
            {
                cell2node[count] = face2nodelist[i][ip];
				count ++;
            }

            for ( int ip = np; ip < 4; ++ ip )
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

        for ( std::size_t i = 0; i < linkmap.size(); ++ i )
        {
            int it = linkmap[i];
            xx[i] = x[it];
            yy[i] = y[it];
            zz[i] = z[it];
        }

		PHWrite(cdata, xx, NumPts);
		PHWrite(cdata, yy, NumPts);
		PHWrite(cdata, zz, NumPts);


		for ( int m = 0; m < nvarplot; ++ m )
		{

			RDouble * qq = new RDouble [NumPts];

			for ( std::size_t i = 0; i < linkmap.size(); ++ i )
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
	//int nVisualVariables = 6;
	//int visualVariables = rho, u, v, w, p, cp;
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

	if ( WantVisualField(grid) )
	{
		SaveDataForTecio(grid, actkey, qn, nvarplot);
	}

	int TecioMission = TEC_SPACE::Writecomplete;
	PHWrite(actkey->GetTecData(), &TecioMission, 1);

	int visualSlice = GlobalDataBase::GetIntParaFromDB("visualSlice");
	if ( grid->GetDim() == THREE_D && visualSlice == 1 )
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
    int nEquation = nLaminar + nChemical;

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
	UnstructGrid *grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nEquation = nLaminar + nChemical;

	FaceProxy *faceProxy = new FaceProxy();
	faceProxy->Create(SEGCTION_LENGTH, nEquation);
	return faceProxy;
}

GeomProxy * NSSolverUnstruct::CreateGeomProxy(Grid *gridIn)
{
	UnstructGrid * grid = UnstructGridCast(gridIn);
    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

	GeomProxy *geomProxy = new GeomProxy();
	geomProxy->Create(SEGCTION_LENGTH);
	return geomProxy;
}

NSFaceValue *NSSolverUnstruct::CreateNSFaceValue(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    int nBoundFace = grid->GetNBoundFace();
    int nTotalCell = grid->GetNTotalCell();
    int nTotal = nTotalCell + nBoundFace;

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();

    int nEquation = nLaminar + nChemical;

    int numberOfSpecies = parameters->GetNumberOfSpecies();

	NSFaceValue *facevar = new NSFaceValue(nEquation, numberOfSpecies, SEGCTION_LENGTH);
	return facevar;
}

RDouble NSSolverUnstruct::P0Wall(RDouble p1, RDouble gama, RDouble mach)
{
	RDouble p2, m2;
	if ( mach > one )
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
	if ( mach > one )
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

    //Gradient *gradientNS = 0;
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
    int nEquation = GetNumberOfEquations();

    RDouble **gradPrimtiveVarFaceX  = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("gradPrimtiveVarFaceX"));
    RDouble **gradPrimtiveVarFaceY  = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("gradPrimtiveVarFaceY"));
    RDouble **gradPrimtiveVarFaceZ  = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("gradPrimtiveVarFaceZ"));
    FieldProxy *qProxy = GetFieldProxy(gridIn, "q");
    RDouble **q     = qProxy->GetField_UNS();
    RDouble **qnode = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("qnode"));

    for (int m = 0; m < nEquation; ++ m)
    {
        string gradientName = parameters->GetGradientName();
        if (gradientName == "ggnode" || gradientName == "ggnodelaplacian")
        {
            gridUnstruct->CompGradientGGNode_NEW(q[m], qnode[m], gradPrimtiveVarFaceX[m], gradPrimtiveVarFaceY[m], gradPrimtiveVarFaceZ[m]);
        }
        else if (gradientName == "ggcell")
        {
            gridUnstruct->CompGradientGGCell(q[m], gradPrimtiveVarFaceX[m], gradPrimtiveVarFaceY[m], gradPrimtiveVarFaceZ[m]);
        }
        else if (gradientName == "lsq")
        {
            gridUnstruct->CompGradientLSQ(q[m], gradPrimtiveVarFaceX[m], gradPrimtiveVarFaceY[m], gradPrimtiveVarFaceZ[m]);
        }    
        else if (gradientName == "ggnode_weight")
        {
            gridUnstruct->CompGradientGGNodeWeight(q[m], gradPrimtiveVarFaceX[m], gradPrimtiveVarFaceY[m], gradPrimtiveVarFaceZ[m]);
        }
        else if (gradientName == "gg_m2")
        {
            gridUnstruct->CompGradientGGModified2(q[m], gradPrimtiveVarFaceX[m], gradPrimtiveVarFaceY[m], gradPrimtiveVarFaceZ[m]);
        }
        else if (gradientName == "ggcellnew")
        {
            gridUnstruct->CompGradientGGCellNew(q[m], gradPrimtiveVarFaceX[m], gradPrimtiveVarFaceY[m], gradPrimtiveVarFaceZ[m]);
        }
        else if (gradientName == "ggcellw")
        {
            gridUnstruct->CompGradientGGCellW(q[m], gradPrimtiveVarFaceX[m], gradPrimtiveVarFaceY[m], gradPrimtiveVarFaceZ[m]);
        }
        else
        {
            TK_Exit::ExceptionExit("No reconstruction method has been choosed ! /n");
        }
    }

    int ifvfd = GlobalDataBase::GetIntParaFromDB("ifvfd");
    if(systemgridtype == MIXGRID && ifvfd == 0)
    {
        //! for mixgrid, gradient is communicated using another function.
    }
    else
    {
        PHSPACE::CommunicateInterfaceValue(gridUnstruct, gradPrimtiveVarFaceX, "dqdx", nEquation);
        PHSPACE::CommunicateInterfaceValue(gridUnstruct, gradPrimtiveVarFaceY, "dqdy", nEquation);
        PHSPACE::CommunicateInterfaceValue(gridUnstruct, gradPrimtiveVarFaceZ, "dqdz", nEquation);
    }

    int nTemperatureModel = parameters->GetTemperatureModel();

    RDouble **gradTemperatureFaceX  = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("gradTemperatureFaceX"));
    RDouble **gradTemperatureFaceY  = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("gradTemperatureFaceY"));
    RDouble **gradTemperatureFaceZ  = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("gradTemperatureFaceZ"));
    RDouble **t     = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("t"));
    RDouble **tnode = reinterpret_cast <RDouble **> (gridUnstruct->GetDataPtr("tnode"));

    for (int m = 0; m < nTemperatureModel; ++ m)
    {
        string gradientName = parameters->GetGradientName();
        if (gradientName == "ggnode" || gradientName == "ggnodelaplacian")
        {
            gridUnstruct->CompGradientGGNode_NEW(t[m], tnode[m], gradTemperatureFaceX[m], gradTemperatureFaceY[m], gradTemperatureFaceZ[m]);
        }
        else if (gradientName == "ggcell")
        {
            gridUnstruct->CompGradientGGCell(t[m], gradTemperatureFaceX[m], gradTemperatureFaceY[m], gradTemperatureFaceZ[m]);
        }
        else if (gradientName == "lsq")
        {
            gridUnstruct->CompGradientLSQ(t[m], gradTemperatureFaceX[m], gradTemperatureFaceY[m], gradTemperatureFaceZ[m]);
        }    
        else if (gradientName == "ggnode_weight")
        {
            gridUnstruct->CompGradientGGNodeWeight(t[m], gradTemperatureFaceX[m], gradTemperatureFaceY[m], gradTemperatureFaceZ[m]);
        }
        else if (gradientName == "gg_m2")
        {
            gridUnstruct->CompGradientGGModified2(t[m], gradTemperatureFaceX[m], gradTemperatureFaceY[m], gradTemperatureFaceZ[m]);
        }
        else if (gradientName == "ggcellnew")
        {
            gridUnstruct->CompGradientGGCellNew(t[m], gradTemperatureFaceX[m], gradTemperatureFaceY[m], gradTemperatureFaceZ[m]);
        }
        else if (gradientName == "ggcellw")
        {
            gridUnstruct->CompGradientGGCellW(t[m], gradTemperatureFaceX[m], gradTemperatureFaceY[m], gradTemperatureFaceZ[m]);
        }
        else
        {
            TK_Exit::ExceptionExit("No reconstruction method has been choosed ! /n");
        }
    }

    if(systemgridtype == MIXGRID && ifvfd == 0)
    {
        //! for mixgrid, gradient is communicated using another function.
    }
    else
    {
        PHSPACE::CommunicateInterfaceValue(gridUnstruct, gradTemperatureFaceX, "dtdx", nTemperatureModel);
        PHSPACE::CommunicateInterfaceValue(gridUnstruct, gradTemperatureFaceY, "dtdy", nTemperatureModel);
        PHSPACE::CommunicateInterfaceValue(gridUnstruct, gradTemperatureFaceZ, "dtdz", nTemperatureModel);
    }
}

void NSSolverUnstruct::InitiallimiterOnlyForMixGrid(Grid *gridIn)
{
    UnstructGrid *grid = UnstructGridCast(gridIn);
    RDouble **q = reinterpret_cast<RDouble **>(grid->GetDataPtr("q"));
    Param_NSSolverUnstruct *parameters = GetControlParameters();

    int limitVariables = parameters->GetLimitVariables();
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

    RDouble **dqdx = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarFaceX"));
    RDouble **dqdy = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarFaceY"));
    RDouble **dqdz = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarFaceZ"));

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
    else if( limitVector == VECTOR_LIMITER )
    {
        for ( int m = 0; m < nEquation; ++ m )
        {
            limit_tmp = limiterOnlyForMixGrid->GetLimiter(m);
            SetLimitBoundary(grid, limit_tmp);
        }
    }
}



LIB_EXPORT void NSSolverUnstruct::ComputePostVisualVariables(Post_Visual * postVisualization)
{
	UnstructGrid *grid = UnstructGridCast(postVisualization->GetGrid());

	RDouble **q = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("q") );
	RDouble **t = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("t") );

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
	if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_TIME_STEP))
	{
		string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_TIME_STEP);
		RDouble * dt  = reinterpret_cast< RDouble * > ( grid->GetDataPtr("dt"  ) );
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

	//! Bell 20130319 add.
	Param_NSSolverUnstruct * parameters = GetControlParameters();
	int viscousType = parameters->GetViscousType();
	string viscousName = parameters->GetViscousName();

	if ( viscousType == 1 )
	{
		if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_VISCOSITY_LAMINAR))
		{
			string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_VISCOSITY_LAMINAR);
			RDouble *visl = reinterpret_cast< RDouble * > ( grid->GetDataPtr("visl") );
			postVisualization->UpdateVisualNodeVarPtr(varName, visl);
		}
	}
	else if ( viscousType > 1 )
	{
		if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_VISCOSITY_LAMINAR))
		{
			string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_VISCOSITY_LAMINAR);
			RDouble *visl = reinterpret_cast< RDouble * > ( grid->GetDataPtr("visl") );
			postVisualization->UpdateVisualNodeVarPtr(varName, visl);
		}

		if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_VISCOSITY_TURBULENT))
		{
			string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_VISCOSITY_TURBULENT);
			RDouble *vist = reinterpret_cast< RDouble * > ( grid->GetDataPtr("vist") );
			postVisualization->UpdateVisualNodeVarPtr(varName, vist);
		}

		if (viscousName.substr(0,6) == "2eq-kw")
		{
			RDouble ** q_turb    = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("q_turb"  ) );
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
			if (viscousName.substr(0,17) == "2eq-kw-menter-sst")   // for blending function F1 and F2
			{
				if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_SST_F1))
				{
					RDouble * blend = reinterpret_cast< RDouble * > ( grid->GetDataPtr("blend"  ) );   // for blengding function F1
					string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_SST_F1);
					postVisualization->UpdateVisualNodeVarPtr(varName, blend);
				}
				if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_SST_F2))
				{
					RDouble * SST_F2 = reinterpret_cast< RDouble * > ( grid->GetDataPtr("SST_F2"  ) );   // for blengding function F2
					string varName = postVisualization->GetVariableName(PHSPACE::VISUAL_SST_F2);
					postVisualization->UpdateVisualNodeVarPtr(varName, SST_F2);
				}
			}
		}
	}

    int visualField = parameters->GetIsPlotVolumeField();
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
		TK_Exit::ExceptionExit("Error: can not visualization vorticity when isPlotVolumeField = 0!");
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
        RDouble **dqdx = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarFaceX"));        
        RDouble *dudx = dqdx[IU];
        postVisualization->UpdateVisualNodeVarPtr("gradientUx", dudx);       
    }
    if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_GRADIENT_UY))
    {
        RDouble **dqdy = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarFaceY"));        
        RDouble *dudy = dqdy[IU];
        postVisualization->UpdateVisualNodeVarPtr("gradientUy", dudy);       
    }
    if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_GRADIENT_VX))
    {
        RDouble **dqdx = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarFaceX"));        
        RDouble *dvdx = dqdx[IV];
        postVisualization->UpdateVisualNodeVarPtr("gradientVx", dvdx);       
    }
    if(postVisualization->IsNeedVisualization(PHSPACE::VISUAL_GRADIENT_VY))
    {
        RDouble **dqdy = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarFaceY"));        
        RDouble *dvdy = dqdy[IV];
        postVisualization->UpdateVisualNodeVarPtr("gradientVy", dvdy);       
    }

	//! Chemical species must be the last.
	int nChemical = parameters->GetChemicalFlag();
	if(nChemical)
	{
		using namespace GAS_SPACE;

		int nm = parameters->GetNSEquationNumber();
		int nLaminar  = parameters->GetLaminarNumber();
		int nEquation = nLaminar + nChemical;

		string * varname = gas->GetNameOfSpecies();
		for ( int m = nm; m < nEquation; ++ m )
		{
			string varName = varname[m-nm];
			postVisualization->UpdateVisualNodeVarPtr(varName, q[m]);
		}
	}
}

LIB_EXPORT void NSSolverUnstruct::ComputePostProbesVariables(Post_Probes * postProbesVar)
{
    UnstructGrid *grid = UnstructGridCast(postProbesVar->GetGrid());

    RDouble **q = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("q") );
    RDouble **t = reinterpret_cast< RDouble ** > ( grid->GetDataPtr("t") );

    using namespace IDX;
    if(postProbesVar->IsNeedMonitoring(PHSPACE::PROBE_DENSITY))
    {
        string varName = postProbesVar->GetProbeVariableName(PHSPACE::PROBE_DENSITY);
        postProbesVar->UpdateProbeVarPtr(varName, q[IDX::IR]);
    }

    if(postProbesVar->IsNeedMonitoring(PHSPACE::PROBE_U))
    {
        string varName = postProbesVar->GetProbeVariableName(PHSPACE::PROBE_U);
        postProbesVar->UpdateProbeVarPtr(varName, q[IDX::IU]);
    }
    if(postProbesVar->IsNeedMonitoring(PHSPACE::PROBE_V))
    {
        string varName = postProbesVar->GetProbeVariableName(PHSPACE::PROBE_V);
        postProbesVar->UpdateProbeVarPtr(varName, q[IDX::IV]);
    }
    if(postProbesVar->IsNeedMonitoring(PHSPACE::PROBE_W))
    {
        string varName = postProbesVar->GetProbeVariableName(PHSPACE::PROBE_W);
        postProbesVar->UpdateProbeVarPtr(varName, q[IDX::IW]);
    }
    if(postProbesVar->IsNeedMonitoring(PHSPACE::PROBE_PRESSURE))
    {
        string varName = postProbesVar->GetProbeVariableName(PHSPACE::PROBE_PRESSURE);
        postProbesVar->UpdateProbeVarPtr(varName, q[IDX::IP]);
    }
    if(postProbesVar->IsNeedMonitoring(PHSPACE::PROBE_TEMPERATURE))
    {
        string varName = postProbesVar->GetProbeVariableName(PHSPACE::PROBE_TEMPERATURE);
        postProbesVar->UpdateProbeVarPtr(varName, t[IDX::ITT]);
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
	if ( !interfaceInformation ) return;

	int iNeighborZone                            = interfaceInformation->FindIthNeighbor(neighborZoneIndex);
	int interfaceNumberBetweenTwoNeighboringZone = interfaceInformation->GetNIFaceOfNeighbor(iNeighborZone);
	int * interfaceIndexContainerForSend         = interfaceInformation->GetFaceIndexForSend(iNeighborZone);
	dataContainer->moveto_begin();

    RDouble **dq = fieldProxy->GetField_UNS();

    for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; -- iGhostLayer)
    {
	for ( int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace )
	{
		int sourceCell;
		int iFace = interfaceIndexContainerForSend[ iLocalFace ];
	    	grid->GetSourceIndex( iFace, iGhostLayer + 1, sourceCell );
		for ( int m = 0; m < nEquation; ++ m )
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
	if ( !interfaceInformation ) return;

	int iNeighborZone                            = interfaceInformation->FindIthNeighbor( neighborZoneIndex );
	int interfaceNumberBetweenTwoNeighboringZone = interfaceInformation->GetNIFaceOfNeighbor(iNeighborZone);
	int * interfaceIndexContainerForReceive      = interfaceInformation->GetFaceIndexForRecv(iNeighborZone);

	dataContainer->moveto_begin();
    RDouble **dq = fieldProxy->GetField_UNS();

    for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; -- iGhostLayer)
    {
	for ( int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace )
	{
		int iFace = interfaceIndexContainerForReceive[ iLocalFace ];
		int targetCell;
	    	grid->GetTargetIndex( iFace, iGhostLayer + 1, targetCell );

		for ( int m = 0; m < nEquation; ++ m )
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
    if (0 == nEquation)
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

void NSSolverUnstruct::OutflowBCRegion(Grid *gridIn, BCRegionUnstruct *bcRegionUnstruct)
{
    // get the basic variables for boundary calculating.
    UnstructGrid *grid   = UnstructGridCast(gridIn);
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();
    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));

    // get the number of equations.
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int neqn = nLaminar + nChemical;

    int iFace, le, re;	
    //! get the face indexes of this boundary condition region.
    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();

    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        // iFace is the face number in the set of faceIndex.
        iFace = *iter;
        le = leftCellofFace[iFace];
        re = rightCellofFace[iFace];		
        for (int m = 0; m < neqn; ++ m)
        {
            q[m][re] = q[m][le];
        }
    }
}

void NSSolverUnstruct::SymmetryBCRegion(Grid *gridIn, BCRegionUnstruct *bcRegionUnstruct)
{
    // get the basic variables for boundary calculating.
    UnstructGrid *grid   = UnstructGridCast(gridIn);
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();
    RDouble   *xfn       = grid->GetFaceNormalX();
    RDouble   *yfn       = grid->GetFaceNormalY();
    RDouble   *zfn       = grid->GetFaceNormalZ();
    RDouble   *vgn       = grid->GetFaceNormalVelocity();
    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));

    // get the number of equations.
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int neqn = nLaminar + nChemical;
    
    int iFace, le, re;	
    RDouble um, vm, wm, vn;
    //! get the face indexes of this boundary condition region.
    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();
    
    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        iFace = *iter;
        le = leftCellofFace[iFace];
        re = rightCellofFace[iFace];		
        
        for (int m = 0; m < neqn; ++ m)
        {
            q[m][re] = q[m][le];
        }

        um = q[IU][le];
        vm = q[IV][le];
        wm = q[IW][le];
        vn = xfn[iFace] * um + yfn[iFace] * vm + zfn[iFace] * wm - vgn[iFace];
               
        q[IU][re] = um - two * xfn[iFace] * vn;
        q[IV][re] = vm - two * yfn[iFace] * vn;
        q[IW][re] = wm - two * zfn[iFace] * vn;
    }
}

void NSSolverUnstruct::ViscousIsotropicWallBCRegion(Grid *gridIn, BCRegionUnstruct *bcRegionUnstruct)
{
    // get the basic variables for boundary calculating.
    UnstructGrid *grid   = UnstructGridCast(gridIn);
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();
    RDouble   *xtn    = grid->GetNodeVelocityX();
    RDouble   *ytn    = grid->GetNodeVelocityY();
    RDouble   *ztn    = grid->GetNodeVelocityZ();
    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));

    // get the number of equations.
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int neqn = nLaminar + nChemical;
    RDouble refDimensionalTemperature = parameters->GetRefDimensionalTemperature();

    //get the variables defined in boundarycondition.hypara file.
    Data_Param * bcData = bcRegionUnstruct->GetBCDataParam();
    RDouble wallTemperature;
    if(bcData)
    {
        if (!bcData->IsExist("wallTemperature", PHDOUBLE, 1))
        {
            wallTemperature = GlobalDataBase::GetDoubleParaFromDB("wallTemperature");
            bcData->UpdateData("wallTemperature", &wallTemperature, PHDOUBLE, 1);
        }
        else
        {
            bcData->GetData("wallTemperature", &wallTemperature, PHDOUBLE, 1);
        }
    }
    else
    {
        wallTemperature = GlobalDataBase::GetDoubleParaFromDB("wallTemperature");
    }
    wallTemperature /= refDimensionalTemperature;

    int iFace, le, re;	
    //! get the face indexes of this boundary condition region.
    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();

    RDouble omav = one;
    using namespace GAS_SPACE;
    RDouble coefficientofstateEquation = gas->GetCoefficientOfStateEquation();

    RDouble rg, pg, tg;     //! g means ghost cell.
    RDouble * prims;
    if (nChemical == 1)
    {
        prims = new RDouble [neqn];
    }

    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        iFace = *iter;
        le = leftCellofFace[iFace];
        re = rightCellofFace[iFace];	

        for (int m = 0; m < neqn; ++ m)
        {
            q[m][re] = q[m][le];
        }

        q[IU][re] = - q[IU][le] + two * xtn[iFace];
        q[IV][re] = - q[IV][le] + two * ytn[iFace];
        q[IW][re] = - q[IW][le] + two * ztn[iFace];

        if (nChemical == 1)
        {
            for (int m = 0; m < neqn; ++ m)
            {
                prims[m] = q[m][le];
            }
            omav = gas->ComputeMolecularWeightReciprocal(prims);
        }

        RDouble tm = (1.0 / coefficientofstateEquation) * q[IP][le] / q[IR][le];

        pg = q[IP][le];
        tg = two * wallTemperature - tm;
        if (tg <= 0.0) tg = wallTemperature;

        rg = pg / (coefficientofstateEquation * tg * omav);

        q[IR][re] = rg;
        q[IP][re] = pg;
    }

    if (nChemical == 1)
    {
        delete [] prims;
    }
}

void NSSolverUnstruct::ViscousAdiabaticWallBCRegion(Grid *gridIn, BCRegionUnstruct *bcRegionUnstruct)
{
    // get the basic variables for boundary calculating.
    UnstructGrid *grid   = UnstructGridCast(gridIn);
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();
    RDouble   *xtn    = grid->GetNodeVelocityX();
    RDouble   *ytn    = grid->GetNodeVelocityY();
    RDouble   *ztn    = grid->GetNodeVelocityZ();
    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));

    // get the number of equations.
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int neqn = nLaminar + nChemical;

    Data_Param *bcData = bcRegionUnstruct->GetBCDataParam();

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

    int iFace, le, re;	
    //! get the face indexes of this boundary condition region.
    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();
    
    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        iFace = *iter;
        le = leftCellofFace[iFace];
        re = rightCellofFace[iFace];	

        for (int m = 0; m < neqn; ++ m)
        {
            q[m][re] = q[m][le];
        }

        q[IU][re] = - q[IU][le] + two * (xtn[iFace] + uWall);
        q[IV][re] = - q[IV][le] + two * (ytn[iFace] + vWall);
        q[IW][re] = - q[IW][le] + two * (ztn[iFace] + wWall);
    }
}

void NSSolverUnstruct::FarfieldBCRegion(Grid *gridIn, BCRegionUnstruct *bcRegionUnstruct)
{
    // get the basic variables for boundary calculating.
    UnstructGrid *grid   = UnstructGridCast(gridIn);
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();
    RDouble   *xfn       = grid->GetFaceNormalX();
    RDouble   *yfn       = grid->GetFaceNormalY();
    RDouble   *zfn       = grid->GetFaceNormalZ();
    RDouble   *vgn       = grid->GetFaceNormalVelocity();
    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    RDouble *gama = reinterpret_cast< RDouble *  > (grid->GetDataPtr("gama"));

    // get the number of equations.
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    RDouble refGama = parameters->GetRefGama();
    RDouble *primitiveVarFarfield = parameters->GetPrimitiveVarFarfield();
    
    int neqn = nLaminar + nChemical;

    int iFace, le, re;
    //! get the face indexes of this boundary condition region.
    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();

    RDouble rin,uin,vin,win,pin,cin,vni;
    RDouble vno;
    RDouble lam1,lam2,lam3,gamm1,vei;
    RDouble vtx ,vty ,vtz ,entr ;
    RDouble riemm,riemp;
    RDouble rb,ub,vb,wb,pb,vnb,cb;   
    
    //! Infinite variables.
    RDouble roo = primitiveVarFarfield[IR];
    RDouble uoo = primitiveVarFarfield[IU];
    RDouble voo = primitiveVarFarfield[IV];
    RDouble woo = primitiveVarFarfield[IW];
    RDouble poo = primitiveVarFarfield[IP];
    RDouble coo = sqrt(ABS(refGama * poo / roo));

    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        iFace = *iter;
        le = leftCellofFace[iFace];
        re = rightCellofFace[iFace];	

        //! Inner point variables.
        rin = q[IR][le];
        uin = q[IU][le];
        vin = q[IV][le];
        win = q[IW][le];
        pin = q[IP][le];

        vno = xfn[iFace] * uoo + yfn[iFace] * voo + zfn[iFace] * woo - vgn[iFace];
        vni = xfn[iFace] * uin + yfn[iFace] * vin + zfn[iFace] * win - vgn[iFace];

        cin = sqrt(ABS(gama[le]  * pin / rin));
        gamm1 = gama[le] - one;

        lam1  = vni;
        lam2  = vni - cin;
        lam3  = vni + cin;
        vei = sqrt(uin * uin + vin * vin + win * win);

        //! Supersonic.
        if (vei > cin)
        {
            if (vni >= 0.0)
            {
                for (int m = 0; m < neqn; ++ m)
                {
                    q[m][re] = q[m][le];
                }
            }
            else
            {
                for (int m = 0; m < neqn; ++ m)
                {
                    q[m][re] = primitiveVarFarfield[m];
                }
            }
            continue;
        }

         //! Subsonic.
        riemp = vni + 2.0 * cin / gamm1;
        riemm = vno - 2.0 * coo / gamm1;
        vnb =  half * ( riemp + riemm );
        cb  =  0.25 * ( riemp - riemm ) * gamm1;
        if (vnb >= 0.0)
        {
            entr = pin / pow(rin, gama[le]);
            vtx  = uin - xfn[iFace] * vni;
            vty  = vin - yfn[iFace] * vni;
            vtz  = win - zfn[iFace] * vni;
        }
        else
        {
            //! Inlet.
            entr = poo / pow(roo, gama[le]);
            vtx  = uoo - xfn[iFace] * vno;
            vty  = voo - yfn[iFace] * vno;
            vtz  = woo - zfn[iFace] * vno;
        }

        rb  = pow((cb * cb / (entr * gama[le])), one / gamm1);
        ub  = vtx + xfn[iFace] * vnb;
        vb  = vty + yfn[iFace] * vnb;
        wb  = vtz + zfn[iFace] * vnb;
        pb  = cb * cb * rb / gama[le];

        q[IR][re] = rb;
        q[IU][re] = ub;
        q[IV][re] = vb;
        q[IW][re] = wb;
        q[IP][re] = pb;
    }
}

void NSSolverUnstruct::PreconFarfieldBCRegion(Grid *gridIn, BCRegionUnstruct *bcRegionUnstruct)
{
    // get the basic variables for boundary calculating.
    UnstructGrid *grid   = UnstructGridCast(gridIn);
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();
    RDouble   *xfn       = grid->GetFaceNormalX();
    RDouble   *yfn       = grid->GetFaceNormalY();
    RDouble   *zfn       = grid->GetFaceNormalZ();
    RDouble   *vgn       = grid->GetFaceNormalVelocity();
    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
    RDouble *gama = reinterpret_cast< RDouble *  > (grid->GetDataPtr("gama"));

    // get the number of equations.
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    RDouble refGama = parameters->GetRefGama();
    RDouble *primitiveVarFarfield = parameters->GetPrimitiveVarFarfield();
    RDouble refMachNumber = parameters->GetRefMachNumber();

    int iFace, le, re;
    //! get the face indexes of this boundary condition region.
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
            if ( vni >= 0.0 )
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

void NSSolverUnstruct::InflowBCRegion(Grid *gridIn, BCRegionUnstruct *bcRegionUnstruct)
{
    // get the basic variables for boundary calculating.
    UnstructGrid *grid   = UnstructGridCast(gridIn);
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();
    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));

    // get the number of equations.
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    RDouble *primitiveVarFarfield = parameters->GetPrimitiveVarFarfield();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int neqn = nLaminar + nChemical;

    int iFace, re;	
    //! get the face indexes of this boundary condition region.
    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();
   
    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        iFace = *iter;
        re = rightCellofFace[iFace];		
        for (int m = 0; m < neqn; ++ m)
        {
            q[m][re] = primitiveVarFarfield[m];
        }
    }
}

void NSSolverUnstruct::PressureOutletBCRegion(Grid *gridIn, BCRegionUnstruct *bcRegionUnstruct)
{
    // get the basic variables for boundary calculating.
    UnstructGrid *grid   = UnstructGridCast(gridIn);
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();
    RDouble   *xfn    = grid->GetFaceNormalX();
    RDouble   *yfn    = grid->GetFaceNormalY();
    RDouble   *zfn    = grid->GetFaceNormalZ();
    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));

    // get the number of equations.
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int neqn = nLaminar + nChemical;

    Data_Param *bcData = bcRegionUnstruct->GetBCDataParam();

    RDouble staticPressure;
    bcData->GetData("staticPressure", &staticPressure, PHDOUBLE, 1);

    //! Dimensional.
    RDouble refDimensionalDensity = parameters->GetRefDimensionalDensity();
    RDouble refDimensionalSonicSpeed = parameters->GetRefDimensionalSonicSpeed();
    RDouble refMachNumber = parameters->GetRefMachNumber();
    RDouble refGama = parameters->GetRefGama();
    RDouble velocity = refMachNumber * refDimensionalSonicSpeed;
    RDouble dynamic_pressure = refDimensionalDensity * velocity * velocity;

    staticPressure /= dynamic_pressure;

    int iFace, le, re;	
    //! get the face indexes of this boundary condition region.
    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();

    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        iFace = *iter;
        le = leftCellofFace[iFace];
        re = rightCellofFace[iFace];		
        for (int m = 0; m < neqn; ++ m)
        {
            q[m][re] = q[m][le];
        }

        //RDouble cc = sqrt( refGama * q[IP][le] / q[IR][le]);
        //RDouble vn = q[IU][le] * xfn[iFace] + q[IV][le] * yfn[iFace] + q[IW][le] * zfn[iFace];
        //if(vn < cc)
        //{
            q[IP][re] = staticPressure;
        //}
        //else
        //{
        //    q[IP][re] = q[IP][le];
        //}
        
    }
}

void NSSolverUnstruct::MassFlowOutletBCRegion(Grid *gridIn, BCRegionUnstruct *bcRegionUnstruct)
{
    calculateBCArea(gridIn, bcRegionUnstruct);
    calculateMassFluxRatio(gridIn, bcRegionUnstruct);

    using namespace IDX;
    //! Dimensional.
    Param_NSSolverUnstruct *parameters = GetControlParameters();  
    RDouble refDimensionalTemperature = parameters->GetRefDimensionalTemperature();
    RDouble refMachNumber = parameters->GetRefMachNumber();
    RDouble refDimensionalDensity = parameters->GetRefDimensionalDensity();
    RDouble refDimensionalSonicSpeed = parameters->GetRefDimensionalSonicSpeed();			
    RDouble refGama = parameters->GetRefGama();

    Data_Param *bcData = bcRegionUnstruct->GetBCDataParam();
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

    // get the basic variables for boundary calculating.
    UnstructGrid *grid   = UnstructGridCast(gridIn);
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();
    RDouble *xfn = grid->GetFaceNormalX();
    RDouble *yfn = grid->GetFaceNormalY();
    RDouble *zfn = grid->GetFaceNormalZ();
    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
  	
    //! get the face indexes of this boundary condition region.
    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();
    int iFace, le, re;

    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        iFace = *iter;
        le = leftCellofFace[iFace];
        re = rightCellofFace[iFace];	

        //! calculate the temperature of the inside cell.
        RDouble staticTemperature = q[IP][le] / q[IR][le] / coefficientOfStateEquation;
        RDouble tb   = staticTemperature + mfTbCoeff * (q[IU][le] * q[IU][le] +
            q[IV][le] * q[IV][le] + q[IW][le] * q[IW][le]);

        if(tb < TINY) 
        {
            tb = staticTemperature;
        }

        RDouble pb   = q[IP][le] * pow( tb / staticTemperature, refGama / (refGama - 1.0));
        RDouble rb = pb / tb / coefficientOfStateEquation;

        RDouble ub   = q[IU][le] * massFluxRatio;
        RDouble vb   = q[IV][le] * massFluxRatio;
        RDouble wb   = q[IW][le] * massFluxRatio;

        //! judge the unusual situation, Vn < 0.
        RDouble Vn  = ub * xfn[iFace] + vb * yfn[iFace] + wb * zfn[iFace];        
        if(Vn < 0.0)
        {
            Vn = fluxDensity / rb;
            ub = Vn * xfn[iFace];
            vb = Vn * yfn[iFace];
            wb = Vn * zfn[iFace];
        } 

        //! if the outlet is supersonic flow.
        RDouble cb = sqrt(refGama * coefficientOfStateEquation * tb);
        if(Vn >= cb)
        {
            rb = q[IR][le];
            ub = q[IU][le];
            vb = q[IV][le];
            wb = q[IW][le];
            pb = q[IP][le];
        }

        //! allocate value to the ghost cell.
        q[IR][re] = rb;
        q[IU][re] = ub;
        q[IV][re] = vb;
        q[IW][re] = wb;
        q[IP][re] = pb;
    }
}

void NSSolverUnstruct::calculateBCArea(Grid *gridIn, BCRegionUnstruct *bcRegionUnstruct)
{
    Data_Param *bcData = bcRegionUnstruct->GetBCDataParam();

    //BC total area does not change with calculating, so it only calculates once.
    if (bcData->IsExist("BCTotalArea", PHDOUBLE, 1))
    {
        return;
    }

    using namespace PHMPI;
    int zoneID = gridIn->GetZoneID();
    string bcName = bcRegionUnstruct->GetBCName();
    RDouble BCTotalArea = 0.0;

    //one processor may owe some zones, so calculate all BC area in one processor.
    int nLocalZones = GetNumberofLocalZones();
    for (int iZone = 0; iZone < nLocalZones; iZone++)
    {
        int iZoneID = GetLocalZoneIDToGlobalZoneID(iZone);
        Grid * iGrid = PHSPACE::GetGrid(iZoneID);
        UnstructGrid *grid = UnstructGridCast(iGrid);      
        
        BCRegionUnstruct **bcRegionUnstruct = grid->GetBCRegionUnstruct();
        int nBCRegionUnstruct = grid->GetNBCRegionUnstruct();
        for (int iBCRegionUnstruct = 0; iBCRegionUnstruct < nBCRegionUnstruct; iBCRegionUnstruct++)
        {
            string jBCName = bcRegionUnstruct[iBCRegionUnstruct]->GetBCName();

            //judge the same BC region according the BC name.
            if (jBCName == bcName)
            {
                int iFace;
                vector<int> *faceIndex = bcRegionUnstruct[iBCRegionUnstruct]->GetFaceIndex();
                RDouble   *area    = grid->GetFaceArea();
                for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
                {
                    iFace = *iter;
                    BCTotalArea += area[iFace];
                }
            }
        }
    }

    //get the total area in all processors with the same BC name.
    RDouble GLB_BCTotalArea;
    PH_AllReduce(&BCTotalArea, &GLB_BCTotalArea, 1, MPI_SUM);
    BCTotalArea = GLB_BCTotalArea;

    for (int iZone = 0; iZone < nLocalZones; iZone++)
    {
        int iZoneID = GetLocalZoneIDToGlobalZoneID(iZone);
        Grid * iGrid = PHSPACE::GetGrid(iZoneID);
        UnstructGrid *grid = UnstructGridCast(iGrid); 

        BCRegionUnstruct **bcRegionUnstruct = grid->GetBCRegionUnstruct();
        int nBCRegionUnstruct = grid->GetNBCRegionUnstruct();
        for (int iBCRegionUnstruct = 0; iBCRegionUnstruct < nBCRegionUnstruct; iBCRegionUnstruct++)
        {
            string jBCName = bcRegionUnstruct[iBCRegionUnstruct]->GetBCName();

            if (jBCName == bcName)
            {
                Data_Param *bcData = bcRegionUnstruct[iBCRegionUnstruct]->GetBCDataParam();
                bcData->UpdateData("BCTotalArea", &BCTotalArea, PHDOUBLE, 1);
            }
        }                       
    }   
}

void NSSolverUnstruct::calculateMassFluxRatio(Grid *gridIn, BCRegionUnstruct *bcRegionUnstruct)
{
    using namespace PHMPI;
    int zoneID = gridIn->GetZoneID();
    string bcName = bcRegionUnstruct->GetBCName();
    Data_Param * bcData = bcRegionUnstruct->GetBCDataParam();
    RDouble massFlow = 0.0;
    bcData->GetData("massFlow", &massFlow, PHDOUBLE, 1);
    RDouble massFlux = 0.0;

    //one processor may owe some zones, so calculate all BC mass flux in one processor.
    //in fact, massflux need calculte once in one step, but I dont find a good method to judge if it is
    //calculated or not, so it will calculate more times, may cost a little more time.
    int nLocalZones = GetNumberofLocalZones();
    for (int iZone = 0; iZone < nLocalZones; iZone++)
    {
        int iZoneID = GetLocalZoneIDToGlobalZoneID(iZone);
        Grid * iGrid = PHSPACE::GetGrid(iZoneID);
        UnstructGrid *grid = UnstructGridCast(iGrid);      

        BCRegionUnstruct **bcRegionUnstructLocal = grid->GetBCRegionUnstruct();
        int nBCRegionUnstruct = grid->GetNBCRegionUnstruct();
        RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));
        for (int iBCRegionUnstruct = 0; iBCRegionUnstruct < nBCRegionUnstruct; iBCRegionUnstruct++)
        {
            string jBCName = bcRegionUnstructLocal[iBCRegionUnstruct]->GetBCName();
            
            //judge the same BC region according the BC name.
            if (jBCName == bcName)
            {
                int iFace;
                vector<int> *faceIndex = bcRegionUnstructLocal[iBCRegionUnstruct]->GetFaceIndex();
                RDouble *area = grid->GetFaceArea();
                int *leftCellofFace  = grid->GetLeftCellOfFace();
                RDouble *xfn = grid->GetFaceNormalX();
                RDouble *yfn = grid->GetFaceNormalY();
                RDouble *zfn = grid->GetFaceNormalZ();
                for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
                {
                    iFace = *iter;
                    int le = leftCellofFace[iFace];

                    RDouble vn_c1 = q[IU][le] * xfn[iFace] + q[IU][le] * yfn[iFace] + q[IW][le] * zfn[iFace];
                    massFlux += q[IR][le] * vn_c1 * area[iFace];
                }
            }
        }
    }
    
    //get the total area in all processors with the same BC name.
    RDouble GLB_MassFlux;
    PH_AllReduce(&massFlux, &GLB_MassFlux, 1, MPI_SUM);

    Param_NSSolverUnstruct *parameters = GetControlParameters();
    RDouble refMachNumber = parameters->GetRefMachNumber();
    RDouble refDimensionalDensity = parameters->GetRefDimensionalDensity();
    RDouble refDimensionalSonicSpeed = parameters->GetRefDimensionalSonicSpeed();
    
    //calculte the dimensional mass flux in this BC name region.
    GLB_MassFlux *= (refDimensionalDensity * refDimensionalSonicSpeed * refMachNumber);
    RDouble massFluxRatio =  massFlow / (GLB_MassFlux + TINY);
    massFluxRatio = MAX(0.8, MIN(massFluxRatio, 1.2));

    bcData->UpdateData("massFluxRatio", &massFluxRatio, PHDOUBLE, 1);
}

void NSSolverUnstruct::MassFlowInletBCRegion(Grid *gridIn, BCRegionUnstruct *bcRegionUnstruct)
{
    calculateBCArea(gridIn, bcRegionUnstruct);

    //!  get the basic mesh and flow field information.
    UnstructGrid *grid   = UnstructGridCast(gridIn);
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();
    RDouble *xfn  = grid->GetFaceNormalX();
    RDouble *yfn  = grid->GetFaceNormalY();
    RDouble *zfn  = grid->GetFaceNormalZ();
    RDouble *area = grid->GetFaceArea();
    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));

    // get the number of equations.
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int neqn = nLaminar + nChemical;

    int iFace, le, re;	
    //!  get the face indexes of this boundary condition region.
    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();

    Data_Param *bcData = bcRegionUnstruct->GetBCDataParam();

    // get the basic variables for mass flow in boundary.
    RDouble direction_inlet[3], massFlow, totalTemperature;
    bcData->GetData("massFlow", &massFlow, PHDOUBLE, 1);
    bcData->GetData("totalTemperature", &totalTemperature, PHDOUBLE, 1);
    bcData->GetData("direction_inlet", direction_inlet, PHDOUBLE, 3);

    RDouble BCTotalArea;
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
    RDouble refDimensionalPressure = parameters->GetRefDimensionalPressure();
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

        q[IR][re]  = q[IR][le];
        RDouble Vnb   = massFlow/ (refDimensionalDensity*refDimensionalSonicSpeed*refMachNumber) / q[IR][le];  //注意：此处massFlow为单位面积流量

        q[IU][re] = Vnb * direction_inlet[0];
        q[IV][re] = Vnb * direction_inlet[1];
        q[IW][re] = Vnb * direction_inlet[2];

        //! square of the stagnation sound speed.
        RDouble c2   = refGama * coefficientOfStateEquation * totalTemperature; 
        //! square of the critical sound speed.
        RDouble c2_ct = 2.0 * c2 / (refGama+1.0); 
        //! square of the Laval number.
        RDouble lamd2 = Vnb * Vnb / c2_ct;
        RDouble tb    = totalTemperature * (1.0 - (refGama - 1.0)/(refGama + 1.0) * lamd2);
        q[IP][re] = q[IR][le] * coefficientOfStateEquation * tb;
    }
}

void NSSolverUnstruct::PressureInletBCRiemannRegion(Grid *gridIn, BCRegionUnstruct *bcRegionUnstruct)
{
    // get the basic variables for boundary calculating.
    UnstructGrid *grid   = UnstructGridCast(gridIn);
    int *leftCellofFace  = grid->GetLeftCellOfFace();
    int *rightCellofFace = grid->GetRightCellOfFace();
    RDouble *xfn = grid->GetFaceNormalX();
    RDouble *yfn = grid->GetFaceNormalY();
    RDouble *zfn = grid->GetFaceNormalZ();
    RDouble **q  = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));

    // get the number of equations.
    Param_NSSolverUnstruct *parameters = GetControlParameters();
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int neqn = nLaminar + nChemical;

    int iFace, le, re;	
    //! get the face indexes of this boundary condition region.
    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();

    Data_Param *bcData = bcRegionUnstruct->GetBCDataParam();

    RDouble direction_inlet[3], totalPressure, totalTemperature;
    bcData->GetData("totalPressure", &totalPressure, PHDOUBLE, 1);
    bcData->GetData("totalTemperature", &totalTemperature, PHDOUBLE, 1);
    bcData->GetData("direction_inlet", direction_inlet, PHDOUBLE, 3);
    
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

        //get inner cell flow field information.
        RDouble rin = q[IR][le];
        RDouble V2in = q[IU][le] * q[IU][le] + q[IV][le] * q[IV][le]  + q[IW][le] * q[IW][le];
        RDouble pin = q[IP][le];

        //calculate the total energy, total enthalpy and sound speed.
        RDouble totalEnergy = pin / (rin * gama1) + 0.5 * V2in;
        RDouble totalEnthalpy = (refGama * coefficientOfStateEquation / gama1) * totalTemperature;
        RDouble c2 = refGama * pin / rin;

        //calculate the riemann invariants R+.
        RDouble Riemann = 2.0 * sqrt(c2) / gama1;
        Riemann += q[IU][le] * xfn[iFace] + q[IV][le] * yfn[iFace] + q[IW][le] * zfn[iFace];

        RDouble theta = xfn[iFace] * direction_inlet[0] + yfn[iFace] * direction_inlet[1] + 
            zfn[iFace] * direction_inlet[2];

        RDouble totalc2 = gama1 * (totalEnthalpy - (totalEnergy + pin / rin)+0.5 * V2in) + c2;

        RDouble a =  1.0 + 0.5 * gama1 * theta * theta;
        RDouble b = -1.0 * gama1 * theta * Riemann;
        RDouble c =  0.5 * gama1 * Riemann * Riemann - 2.0 * totalc2 / gama1;

        RDouble d = b * b - 4.0 * a * c;
        d = sqrt(max(0.0, d));
        RDouble Vin   = (-b + d)/(2.0*a);
        Vin   = max(0.0, Vin);
        V2in = Vin * Vin;

        c2 = totalc2 - 0.5 * gama1 * V2in;

        //fix the mach number in case of mach > 1.0 .
        RDouble Mach2 = V2in / c2;
        Mach2 = min(1.0, Mach2);
        V2in = Mach2 * c2;
        Vin = sqrt(V2in);
        c2 = totalc2 - 0.5 * gama1 * V2in;

        RDouble tin = c2 / (refGama * coefficientOfStateEquation);
        pin = totalPressure * pow((tin / totalTemperature), refGama / gama1);
        rin = pin/(coefficientOfStateEquation * tin);

        q[IR][re] = rin;
        q[IU][re] = Vin * direction_inlet[0];
        q[IV][re] = Vin * direction_inlet[1];
        q[IW][re] = Vin * direction_inlet[2];
        q[IP][re] = pin;
    }
}

void NSSolverUnstruct::SwapNeighborsDQNonblocking(Grid *gridIn, FieldProxy *dqProxy)
{
	int currentProcessor = PHMPI::GetCurrentProcessorID();
	int numberOfSimulationZones = PHMPI::GetNumberofGlobalZones();
	vector< PH_Request > requestContainer;
	vector< vector< DataContainer * > > receivedDataBuffer;
	vector< vector< DataContainer * > > sendDataBuffer;
	receivedDataBuffer.resize( numberOfSimulationZones );
	sendDataBuffer.resize( numberOfSimulationZones );
	int level = gridIn->GetLevel();

	//! Step 0: Compressing data firstly.
	for ( int iZone = 0; iZone < numberOfSimulationZones; ++ iZone )
	{
		int sendProcessor = PHMPI::GetZoneProcessorID( iZone );		

		ZoneNeighbor * globalNeighborZones = zoneConnectivity->GetZoneNeighbor( iZone );
		std::size_t numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();
		if( currentProcessor == sendProcessor )
		{
			sendDataBuffer[iZone].resize(numberOfNeighbor);
		}

		for ( std::size_t iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor )
		{
			int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);

			int receiveProcessor = PHMPI::GetZoneProcessorID( neighborZone );

			if( currentProcessor != sendProcessor && currentProcessor != receiveProcessor )
			{
				continue;
			}

			//! Allocate and compress the buffers for sending.
			DataContainer * sendBuffer = 0;
			if( currentProcessor == sendProcessor )
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
	for ( int iZone = 0; iZone < numberOfSimulationZones; ++ iZone )
	{
		ZoneNeighbor * globalNeighborZones = zoneConnectivity->GetZoneNeighbor( iZone );
		std::size_t numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();

		//! Communicating.
		for ( std::size_t iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor )
		{
			int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);

			int sendProcessor    = PHMPI::GetZoneProcessorID( iZone );		
			int receiveProcessor = PHMPI::GetZoneProcessorID( neighborZone );

			if( currentProcessor != sendProcessor && currentProcessor != receiveProcessor )
			{
				continue;
			}

			//! Allocate the buffers for sending and receiving.
			DataContainer * receivedBuffer = 0;
			DataContainer * sendBuffer = 0;
			if( currentProcessor == receiveProcessor )
			{
				receivedBuffer = new DataContainer;
				receivedDataBuffer[ iZone ].push_back( receivedBuffer );
			}
			if( currentProcessor == sendProcessor )
			{
				sendBuffer = sendDataBuffer[iZone][iNeighbor];				
			}

			int tag = iZone;

			//if ( sendProcessor == receiveProcessor ) continue;

			//! Communication.
			if ( currentProcessor == sendProcessor )
			{
				if ( sendProcessor != receiveProcessor )
				{
					streamsize nlen = sendBuffer->size();

					//! Send the data to neighbors.
					send(sendBuffer, receiveProcessor, nlen, requestContainer, tag);
				}
			}
			if ( currentProcessor == receiveProcessor )
			{
				//! Get the neighbor order: which order that 'iZone' on neighbors of neighborZone.
				int neighborOrer = -1;
				ZoneNeighbor * globalNeighborZonesTemp = zoneConnectivity->GetZoneNeighbor( neighborZone );
				std::size_t numberOfNeighborTemp = globalNeighborZonesTemp->GetNumberOfNeighbor();

				for ( std::size_t neighborID = 0; neighborID < numberOfNeighborTemp; ++ neighborID )
				{
					int neighborZoneTemp = globalNeighborZonesTemp->GetZoneIndexOfNeighbor(neighborID);
					if(neighborZoneTemp == iZone)
					{
						neighborOrer = static_cast<int>(neighborID);
						break;
					}
				}

				ASSERT(neighborOrer != -1);

				if ( sendProcessor != receiveProcessor )
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
	if( PHMPI::GetNumberOfProcessor() > 1 )
	{
		PH_Wait( static_cast<int>(requestContainer.size()), &(requestContainer[0]) );
	}

	//! Step3: Translate the data container.
	for ( int iZone = 0; iZone < numberOfSimulationZones; ++ iZone )
	{
		ZoneNeighbor * globalNeighborZones = zoneConnectivity->GetZoneNeighbor( iZone );
		std::size_t numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();

		int count = 0;
		for ( std::size_t iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor )
		{
			int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);

			int receiveProcessor = PHMPI::GetZoneProcessorID( neighborZone );

			if ( currentProcessor == receiveProcessor )
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
	for ( std::size_t iDim = 0; iDim < receivedDataBuffer.size(); ++ iDim )
	{
		for ( std::size_t jDim = 0; jDim < receivedDataBuffer[ iDim ].size(); ++ jDim )
		{
			delete receivedDataBuffer[ iDim ][ jDim ];
		}
	}
	receivedDataBuffer.clear();
	for ( std::size_t iDim = 0; iDim < sendDataBuffer.size(); ++ iDim )
	{
		for ( std::size_t jDim = 0; jDim < sendDataBuffer[ iDim ].size(); ++ jDim )
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
    
    RDouble **dqdx = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarFaceX"));
    RDouble **dqdy = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarFaceY"));
    RDouble **dqdz = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarFaceZ"));

    RDouble **dtdx  = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradTemperatureFaceX"));
    RDouble **dtdy  = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradTemperatureFaceY"));
    RDouble **dtdz  = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradTemperatureFaceZ"));

    int limitVector = GlobalDataBase::GetIntParaFromDB("limitVector");

    RDouble *limit = 0;
    RDouble **limit2D = 0;
    if(limitVector == 0)
    {
        //! scalar
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
            //! scalar
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
    
    RDouble **dqdx = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarFaceX"));
    RDouble **dqdy = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarFaceY"));
    RDouble **dqdz = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradPrimtiveVarFaceZ"));

    RDouble **dtdx  = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradTemperatureFaceX"));
    RDouble **dtdy  = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradTemperatureFaceY"));
    RDouble **dtdz  = reinterpret_cast <RDouble **> (grid->GetDataPtr("gradTemperatureFaceZ"));

    int limitVector = GlobalDataBase::GetIntParaFromDB("limitVector");

    RDouble *limit = 0;
    RDouble **limit2D = 0;
    if(limitVector == 0)
    {
        //! scalar
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
            //! scalar
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

    if(viscousType <= 0)
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
        //! laminar
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

    if(viscousType <= 0)
    {
        //! Euler
        RDouble **Inviscidflux = InviscidfaceProxyOnlyForMixGrid->GetFlux();

        RDouble Inviscidflux_str[5];
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
        //! laminar
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

