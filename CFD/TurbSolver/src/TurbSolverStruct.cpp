#include "TurbSolverStruct.h"
#include "OversetGridFactory.h"
#include "BC_Struct.h"
#include "JetBC.h"
#include "GradientOperation.h"
#include "Math_Limiter.h"
#include "TK_Exit.h"
#include "Param_TurbSolverStruct.h"
#include "Residual.h"
#include "FieldProxy.h"
#include "Glb_Dimension.h"
#include "TK_Log.h"
#include "IO_FileName.h"

#ifdef USE_ARRAY_OPT
#include "PHOpt.h"
#endif

using namespace std;

namespace PHSPACE
{
TurbSolverStr::TurbSolverStr()
{
    DESLength = 0;
}

TurbSolverStr::~TurbSolverStr()
{
    DeAllocateGlobalVariables();

    delete DESLength;
    DESLength = NULL;

    FreeControlParameters();
}

void TurbSolverStr::AllocateGlobalVar(Grid *gridIn)
{
    Param_TurbSolverStruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    if (nTurbulentEquation == 0)
    {
        return;
    }

    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    Range IFACE, JFACE, KFACE;
    GetRange(ni, nj, nk, -2, 2, IFACE, JFACE, KFACE);

    Range M(0, nTurbulentEquation - 1);

    RDouble4D *qTurbulent = new RDouble4D(I, J, K, M, fortranArray);
    RDouble4D *residualTurbulent  = new RDouble4D(I, J, K, M, fortranArray);
    RDouble4D *spectumTurbulent = new RDouble4D(I, J, K, M, fortranArray);

    grid->UpdateDataPtr("q_turb"   , qTurbulent);
    grid->UpdateDataPtr("res_turb" , residualTurbulent );
    grid->UpdateDataPtr("spec_turb", spectumTurbulent);

    RDouble4D *gradTurbulentX = new RDouble4D(IFACE, JFACE, KFACE, M, fortranArray);
    RDouble4D *gradTurbulentY = new RDouble4D(IFACE, JFACE, KFACE, M, fortranArray);
    RDouble4D *gradTurbulentZ = new RDouble4D(IFACE, JFACE, KFACE, M, fortranArray);
    grid->UpdateDataPtr("gradTurbulentX", gradTurbulentX);
    grid->UpdateDataPtr("gradTurbulentY", gradTurbulentY);
    grid->UpdateDataPtr("gradTurbulentZ", gradTurbulentZ);

    RDouble4D *gradTurbulentCellCenterX = new RDouble4D(I, J, K, M, fortranArray);
    RDouble4D *gradTurbulentCellCenterY = new RDouble4D(I, J, K, M, fortranArray);
    RDouble4D *gradTurbulentCellCenterZ = new RDouble4D(I, J, K, M, fortranArray);
    grid->UpdateDataPtr("gradTurbulentCellCenterX", gradTurbulentCellCenterX);
    grid->UpdateDataPtr("gradTurbulentCellCenterY", gradTurbulentCellCenterY);
    grid->UpdateDataPtr("gradTurbulentCellCenterZ", gradTurbulentCellCenterZ);

    RDouble4D *rotTurbgradValueX = new RDouble4D(I, J, K, M, fortranArray);
    RDouble4D *rotTurbgradValueY = new RDouble4D(I, J, K, M, fortranArray);
    RDouble4D *rotTurbgradValueZ = new RDouble4D(I, J, K, M, fortranArray);
    grid->UpdateDataPtr("rotTurbgradValueX", rotTurbgradValueX);
    grid->UpdateDataPtr("rotTurbgradValueY", rotTurbgradValueY);
    grid->UpdateDataPtr("rotTurbgradValueZ", rotTurbgradValueZ);

    int isUnsteady = parameters->GetIsUnsteady();
    if (isUnsteady)
    {
        RDouble4D *qTurbUnsteadyn1 = new RDouble4D(I, J, K, M, fortranArray);
        RDouble4D *qTurbUnsteadyn2 = new RDouble4D(I, J, K, M, fortranArray);

        RDouble4D *residualTurbUnsteadyn1 = new RDouble4D(I, J, K, M, fortranArray);
        RDouble4D *residualTurbUnsteadyn2 = new RDouble4D(I, J, K, M, fortranArray);

        //! Store the historical sum of InviscidFlux, ViscousFlux, et al., exclusive of DualTimeSourceFlux, used in UpdateUnsteadyFlow for the Rn in the dualtime method.
        RDouble4D *residualTurbUnsteadyTemporary = new RDouble4D(I, J, K, M, fortranArray);

        grid->UpdateDataPtr("q_turb_unsteady_n1"   , qTurbUnsteadyn1);
        grid->UpdateDataPtr("q_turb_unsteady_n2"   , qTurbUnsteadyn2);
        grid->UpdateDataPtr("res_turb_unsteady_n1" , residualTurbUnsteadyn1);
        grid->UpdateDataPtr("res_turb_unsteady_n2" , residualTurbUnsteadyn2);
        grid->UpdateDataPtr("res_turb_unsteady_tmp", residualTurbUnsteadyTemporary);

        //! Statistical variables for DES simulation.
        int ifStaticsFlowField = parameters->GetIfStaticsFlowField();
        if (ifStaticsFlowField > 0)
        {
            RDouble4D *qAverageTurb = new RDouble4D(I, J, K, M, fortranArray);
            grid->UpdateDataPtr("qAverage_turb", qAverageTurb);
        }
    }

    if (nTurbulentEquation == 2)
    {
        RDouble3D *cross = new RDouble3D(I, J, K, fortranArray);
        grid->UpdateDataPtr("cross", cross);
        RDouble3D *blend = new RDouble3D(I, J, K, fortranArray);  //SST_F1
        grid->UpdateDataPtr("blend", blend);
        RDouble3D *SST_F2 = new RDouble3D(I, J, K, fortranArray); //SST_F2
        grid->UpdateDataPtr("SST_F2", SST_F2);
    }

    int numberOfDimensions = GetDim();
    Range D(1, numberOfDimensions);

    RDouble5D *matrixTurbLeft, *matrixTurbRight;

    matrixTurbLeft = new RDouble5D(I, J, K, M, D, fortranArray);
    matrixTurbRight = new RDouble5D(I, J, K, M, D, fortranArray);

    grid->UpdateDataPtr("mat_turbl", matrixTurbLeft);
    grid->UpdateDataPtr("mat_turbr", matrixTurbRight);
}

void TurbSolverStr::DeAllocateGlobalVar(Grid *gridIn)
{
    Param_TurbSolverStruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();
    if (nTurbulentEquation == 0)
    {
        return;
    }

    StructGrid *grid = StructGridCast(gridIn);

    RDouble4D *qTurbulent        = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb"   ));
    RDouble4D *residualTurbulent = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_turb" ));
    RDouble4D *spectumTurbulent  = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("spec_turb"));

    delete qTurbulent;
    delete residualTurbulent;
    delete spectumTurbulent;

    RDouble4D *gradTurbulentX = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradTurbulentX"));
    RDouble4D *gradTurbulentY = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradTurbulentY"));
    RDouble4D *gradTurbulentZ = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradTurbulentZ"));

    delete gradTurbulentX;
    delete gradTurbulentY;
    delete gradTurbulentZ;

    RDouble4D *gradTurbulentCellCenterX = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradTurbulentCellCenterX"));
    RDouble4D *gradTurbulentCellCenterY = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradTurbulentCellCenterY"));
    RDouble4D *gradTurbulentCellCenterZ = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradTurbulentCellCenterZ"));

    delete gradTurbulentCellCenterX;
    delete gradTurbulentCellCenterY;
    delete gradTurbulentCellCenterZ;

    RDouble4D *rotTurbgradValueX = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("rotTurbgradValueX"));
    RDouble4D *rotTurbgradValueY = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("rotTurbgradValueY"));
    RDouble4D *rotTurbgradValueZ = reinterpret_cast <RDouble4D *> (grid->GetDataPtr("rotTurbgradValueZ"));
    delete rotTurbgradValueX;
    delete rotTurbgradValueY;
    delete rotTurbgradValueZ;

    int isUnsteady = parameters->GetIsUnsteady();
    if (isUnsteady)
    {
        RDouble4D *qTurbUnsteadyn1        = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb_unsteady_n1"));
        RDouble4D *qTurbUnsteadyn2        = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb_unsteady_n2"));
        RDouble4D *residualTurbUnsteadyn1 = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_turb_unsteady_n1"));
        RDouble4D *residualTurbUnsteadyn2 = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_turb_unsteady_n2"));

        //! Store the historical sum of InviscidFlux, ViscousFlux, et al., exclusive of DualTimeSourceFlux, used in UpdateUnsteadyFlow for the Rn in the dualtime method.
        RDouble4D *residualTurbUnsteadyTemporary = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_turb_unsteady_tmp"));

        delete qTurbUnsteadyn1;
        delete qTurbUnsteadyn2;
        delete residualTurbUnsteadyn1;
        delete residualTurbUnsteadyn2;
        delete residualTurbUnsteadyTemporary;

        //! Statistical variables for DES simulation.
        int ifStaticsFlowField = parameters->GetIfStaticsFlowField();
        if (ifStaticsFlowField > 0)
        {
            RDouble4D *qAverageTurb = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("qAverage_turb"));
            delete qAverageTurb;
        }
    }

    if (nTurbulentEquation == 2)
    {
        RDouble3D *cross  = reinterpret_cast<RDouble3D *> (grid->GetDataPtr("cross" ));
        delete cross;
        RDouble3D *blend  = reinterpret_cast<RDouble3D *> (grid->GetDataPtr("blend" ));
        delete blend;
        RDouble3D *SST_F2 = reinterpret_cast<RDouble3D *> (grid->GetDataPtr("SST_F2" ));
        delete SST_F2;
    }

    RDouble5D *matrixTurbLeft, *matrixTurbRight;
    matrixTurbLeft  = reinterpret_cast<RDouble5D *> (grid->GetDataPtr("mat_turbl"));
    matrixTurbRight = reinterpret_cast<RDouble5D *> (grid->GetDataPtr("mat_turbr"));

    delete matrixTurbLeft;
    delete matrixTurbRight;
}

bool TurbSolverStr::JudgeIfRestart()
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

bool TurbSolverStr::JudgeIfReadAverage()
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

void TurbSolverStr::DumpRestart(ActionKey *actkey)
{
    int level = actkey->level;
    StructGrid *grid = StructGridCast(GetGrid(level));

    int ist, ied, jst, jed, kst, ked;
	grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 2);

    Param_TurbSolverStruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    Range M(0, nTurbulentEquation - 1);
    int mst = M.first();
    int med = M.last();

    int outerIterStep = GlobalDataBase::GetIntParaFromDB("outnstep");

    DataContainer *dataContainer = actkey->GetData();
    dataContainer->moveto_begin();
    dataContainer->write(&outerIterStep, sizeof(int));

    RDouble4D &qTurbulent = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb"));

    for (int m = mst; m <= med; ++ m)
    {
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    dataContainer->write(&qTurbulent(i, j, k, m), sizeof(RDouble));
                }
            }
        }
    }

    int isUnsteady = parameters->GetIsUnsteady();
    if (isUnsteady)
    {
        RDouble4D &qTurbUnsteadyn1 = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb_unsteady_n1"));
        RDouble4D &qTurbUnsteadyn2 = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb_unsteady_n2"));

        for (int m = mst; m <= med; ++ m)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    for (int i = ist; i <= ied; ++ i)
                    {
                        dataContainer->write(&qTurbUnsteadyn1(i, j, k, m), sizeof(RDouble));
                    }
                }
            }
        }

        for (int m = mst; m <= med; ++ m)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    for (int i = ist; i <= ied; ++ i)
                    {
                        dataContainer->write(&qTurbUnsteadyn2(i, j, k, m), sizeof(RDouble));
                    }
                }
            }
        }

        RDouble4D &residualTurbulent      = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_turb"));
        RDouble4D &residualTurbUnsteadyn1 = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_turb_unsteady_n1"));
        RDouble4D &residualTurbUnsteadyn2 = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_turb_unsteady_n2"));

        for (int m = mst; m <= med; ++ m)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    for (int i = ist; i <= ied; ++ i)
                    {
                        dataContainer->write(&residualTurbulent(i, j, k, m), sizeof(RDouble));
                    }
                }
            }
        }

        for (int m = mst; m <= med; ++ m)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    for (int i = ist; i <= ied; ++ i)
                    {
                        dataContainer->write(&residualTurbUnsteadyn1(i, j, k, m), sizeof(RDouble));
                    }
                }
            }
        }

        for (int m = mst; m <= med; ++ m)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    for (int i = ist; i <= ied; ++ i)
                    {
                        dataContainer->write(&residualTurbUnsteadyn2(i, j, k, m), sizeof(RDouble));
                    }
                }
            }
        }
    }
}

void TurbSolverStr::DumpRestartH5(ActionKey *actkey)
{
    using namespace PHMPI;
    int currentProcessorID = GetCurrentProcessorID();

    int version = 1;
    int outIterStep = GlobalDataBase::GetIntParaFromDB("outnstep");
    if (currentProcessorID == GetServerProcessorID())
    {
        WriteData(actkey->filepos, &version, "Version");
        WriteData(actkey->filepos, &outIterStep, "outnstep");

        if (IsNeedDESStatistics())
        {
            int nStatisticalStep = 0;
            GlobalDataBase::GetData("nStatisticalStep", &nStatisticalStep, PHINT, 1);
            WriteData(actkey->filepos, &nStatisticalStep, "nStatisticalStep");
        }
    }

    Grid *grid = GetGrid(actkey->level);
    int gridID = grid->GetZoneID();

    StructGrid *strGrid = StructGridCast(GetGrid(actkey->level));
    int nTotalCell = strGrid->GetNTotalCell();

    int ist, ied, jst, jed, kst, ked;
    strGrid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 2);

    Param_TurbSolverStruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    Range M(0, nTurbulentEquation - 1);
    int mst = M.first();
    int med = M.last();

    RDouble4D &qTurbulent = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb"));
    RDouble3D &viscousTurbulent = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("vist"));

    hid_t grploc;
    string grpName;

    ostringstream oss;
    oss<<"Group"<<gridID;
    grpName = oss.str();

    grploc = OpenGroup(actkey->filepos, grpName);

    WriteData(grploc, &nTotalCell, "nTotalCell");
    WriteData(grploc, &qTurbulent(ist, jst, kst, mst), "q_turb");
    WriteData(grploc, &viscousTurbulent(ist, jst, kst), "visturb");

    int isUnsteady = parameters->GetIsUnsteady();
    if (isUnsteady == 1)
    {
        RDouble4D &qTurbUnsteadyn1 = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb_unsteady_n1"));
        WriteData(grploc, &qTurbUnsteadyn1(ist, jst, kst, mst), "q_turb_unsteady_n1");

        RDouble4D &qTurbUnsteadyn2 = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb_unsteady_n2"));
        WriteData(grploc, &qTurbUnsteadyn2(ist, jst, kst, mst), "q_turb_unsteady_n2");

        RDouble4D &residualTurbulent   = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_turb"));
        WriteData(grploc, &residualTurbulent(ist, jst, kst, mst), "res_turb");

        RDouble4D &residualTurbUnsteadyn1 = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_turb_unsteady_n1"));
        WriteData(grploc, &residualTurbUnsteadyn1(ist, jst, kst, mst), "res_turb_unsteady_n1");

        RDouble4D &residualTurbUnsteadyn2 = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_turb_unsteady_n2"));
        WriteData(grploc, &residualTurbUnsteadyn2(ist, jst, kst, mst), "res_turb_unsteady_n2");

        if (IsNeedDESStatistics())
        {
            RDouble4D &qAverageTurb = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("qAverage_turb"));
            WriteData(grploc, &qAverageTurb(ist, jst, kst, mst), "qAverageTurb");
        }
    }

    H5Gclose(grploc);
}

void TurbSolverStr::ReadRestartH5(ActionKey *actkey)
{
    Grid *grid = GetGrid(actkey->level);
    int gridID = grid->GetZoneID();

    StructGrid *strGrid = StructGridCast(GetGrid(actkey->level));
    int nTotalCell = strGrid->GetNTotalCell();

    int ist, ied, jst, jed, kst, ked;
    strGrid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 2);

    Param_TurbSolverStruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    Range M(0, nTurbulentEquation - 1);
    int mst = M.first();
    int med = M.last();

    RDouble4D &qTurbulent = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb"));
    RDouble3D &viscousTurbulent = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("vist"));

    int outIterstepofNS = GlobalDataBase::GetIntParaFromDB("outnstep");
    int outIterStepofTurb = 0;
    ReadData(actkey->filepos, &outIterStepofTurb, "outnstep");

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
        erroeInfo << " Error: the cell number in turb.dat is not equal to the cell number in grid file !" << endl;
        TK_Exit::ExceptionExit(erroeInfo.str());
    }

    ReadData(grploc, &qTurbulent(ist, jst, kst, mst), "q_turb");
    ReadData(grploc, &viscousTurbulent(ist, jst, kst), "visturb");

    int isUnsteady = parameters->GetIsUnsteady();
    if (isUnsteady)
    {
        RDouble4D &qTurbUnsteadyn1        = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb_unsteady_n1"));
        RDouble4D &qTurbUnsteadyn2        = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb_unsteady_n2"));
        RDouble4D &residualTurbulent      = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_turb"));
        RDouble4D &residualTurbUnsteadyn1 = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_turb_unsteady_n1"));
        RDouble4D &residualTurbUnsteadyn2 = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_turb_unsteady_n2"));
        RDouble4D &qAverageTurb = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("qAverage_turb"));

        int ifStartFromSteadyResults = parameters->GetIfStartFromSteadyResults();

        if (ifStartFromSteadyResults)
        {
            if (grid->GetGridID()->GetIndex() == 0)
            {
                PrintToWindow("Restart from steady Turbulent flow field, reset outer step to be zero!\n");
            }

            outIterStepofTurb = 0;
            GlobalDataBase::UpdateData("outnstep", &outIterStepofTurb, PHINT, 1);

            qTurbUnsteadyn1 = qTurbulent;
            qTurbUnsteadyn2 = qTurbulent;

            residualTurbulent      = 0.0;
            residualTurbUnsteadyn1 = 0.0;
            residualTurbUnsteadyn2 = 0.0;
        }
        else
        {
            ReadData(grploc, &qTurbUnsteadyn1(ist, jst, kst, mst), "q_turb_unsteady_n1");
            ReadData(grploc, &qTurbUnsteadyn2(ist, jst, kst, mst), "q_turb_unsteady_n2");
            ReadData(grploc, &residualTurbulent(ist, jst, kst, mst), "res_turb");
            ReadData(grploc, &residualTurbUnsteadyn1(ist, jst, kst, mst), "res_turb_unsteady_n1");
            ReadData(grploc, &residualTurbUnsteadyn2(ist, jst, kst, mst), "res_turb_unsteady_n2");
        }

        int nStatisticalStep = 0;
        int ifStaticsFlowField = parameters->GetIfStaticsFlowField();
        int startStatisticStep = GlobalDataBase::GetIntParaFromDB("startStatisticStep");
        bool isReadAverageFlow = ifStaticsFlowField && (outIterStepofTurb >= startStatisticStep);
        if(isReadAverageFlow && (!ifStartFromSteadyResults))
        {
            ReadData(grploc, &qAverageTurb(ist, jst, kst, mst), "qAverageTurb");

            ReadData(actkey->filepos, &nStatisticalStep, "nStatisticalStep");
            GlobalDataBase::UpdateData("nStatisticalStep", &nStatisticalStep, PHINT, 1);
        }
        else
        {
            qAverageTurb = qTurbulent;
            GlobalDataBase::UpdateData("nStatisticalStep", &nStatisticalStep, PHINT, 1);
        }
    }

    H5Gclose(grploc);

    CompareOutStepOfFlowfieldFile(outIterstepofNS, outIterStepofTurb);
}

void TurbSolverStr::ReadRestart(ActionKey *actkey)
{
    int level = actkey->level;
    StructGrid * grid = StructGridCast(GetGrid(level));

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 2);

    Param_TurbSolverStruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    RDouble4D &qTurbulent = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb"));

    Range M(0, nTurbulentEquation - 1);
    int mst = M.first();
    int med = M.last();

    DataContainer *dataContainer = actkey->GetData();
    dataContainer->moveto_begin();
    int outnstepofTurb = 0;
    dataContainer->read(&outnstepofTurb, sizeof(int));

    int outIterstepofNS = GlobalDataBase::GetIntParaFromDB("outnstep");

    for (int m = mst; m <= med; ++ m)
    {
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    dataContainer->read(&qTurbulent(i, j, k, m), sizeof(RDouble));
                }
            }
        }
    }

    int isUnsteady = parameters->GetIsUnsteady();
    if (isUnsteady)
    {
        RDouble4D &qTurbulentn1 = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb_unsteady_n1"));
        RDouble4D &qTurbulentn2 = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb_unsteady_n2"));

        RDouble4D &residualTurbulent   = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_turb"           ));
        RDouble4D &residualTurbulentn1 = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_turb_unsteady_n1"));
        RDouble4D &residualTurbulentn2 = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_turb_unsteady_n2"));

        int ifStartFromSteadyResults = parameters->GetIfStartFromSteadyResults();

        if (ifStartFromSteadyResults)
        {
            if (grid->GetGridID()->GetIndex() == 0)
            {
                PrintToWindow("Restart from steady Turbulent flow field, reset outer step to be zero!\n");
            }

            //! Start from steady flow field.
            //! Reset the outer step when start from steady flow.
            outnstepofTurb = 0;
            GlobalDataBase::UpdateData("outnstep", &outnstepofTurb, PHINT, 1);

            qTurbulentn1 = qTurbulent;
            qTurbulentn2 = qTurbulent;

            residualTurbulent   = 0.0;
            residualTurbulentn1 = 0.0;
            residualTurbulentn2 = 0.0;
        }
        else
        {
            //! Start from unsteady flow field.
            for (int m = mst; m <= med; ++ m)
            {
                for (int k = kst; k <= ked; ++ k)
                {
                    for (int j = jst; j <= jed; ++ j)
                    {
                        for (int i = ist; i <= ied; ++ i)
                        {
                            dataContainer->read(&qTurbulentn1(i, j, k, m), sizeof(RDouble));
                        }
                    }
                }
            }

            for (int m = mst; m <= med; ++ m)
            {
                for (int k = kst; k <= ked; ++ k)
                {
                    for (int j = jst; j <= jed; ++ j)
                    {
                        for (int i = ist; i <= ied; ++ i)
                        {
                            dataContainer->read(&qTurbulentn2(i, j, k, m), sizeof(RDouble));
                        }
                    }
                }
            }

            for (int m = mst; m <= med; ++ m)
            {
                for (int k = kst; k <= ked; ++ k)
                {
                    for (int j = jst; j <= jed; ++ j)
                    {
                        for (int i = ist; i <= ied; ++ i)
                        {
                            dataContainer->read(&residualTurbulent(i, j, k, m), sizeof(RDouble));
                        }
                    }
                }
            }
         
            for (int m = mst; m <= med; ++ m)
            {
                for (int k = kst; k <= ked; ++ k)
                {
                    for (int j = jst; j <= jed; ++ j)
                    {
                        for (int i = ist; i <= ied; ++ i)
                        {
                            dataContainer->read(&residualTurbulentn1(i, j, k ,m), sizeof(RDouble));
                        }
                    }
                }
            }
         
            for (int m = mst; m <= med; ++ m)
            {
                for (int k = kst; k <= ked; ++ k)
                {
                    for (int j = jst; j <= jed; ++ j)
                    {
                        for (int i = ist; i <= ied; ++ i)
                        {
                            dataContainer->read(&residualTurbulentn2(i, j, k ,m), sizeof(RDouble));
                        }
                    }
                }
            }
        }
    }

    CompareOutStepOfFlowfieldFile(outIterstepofNS, outnstepofTurb);
}

void TurbSolverStr::GetResidual(ActionKey *actkey)
{
    int level = actkey->level;
    StructGrid *grid = StructGridCast(GetGrid(level));

    RDouble4D &residualTurbulent = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_turb"));
    
    Param_TurbSolverStruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    PHSPACE::ComputeResidualonGrid(grid, actkey, residualTurbulent, nTurbulentEquation);
}

void TurbSolverStr::ZeroResiduals(Grid *gridIn)
{
    Param_TurbSolverStruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();
    if (nTurbulentEquation == 0)
    {
        return;
    }

    RDouble4D &residualTurbulent = *reinterpret_cast<RDouble4D *> (gridIn->GetDataPtr("res_turb"));
    residualTurbulent = 0;
}

void TurbSolverStr::InitSpectrum(Grid *gridIn)
{
    RDouble4D &spectumTurbulent = *reinterpret_cast<RDouble4D *> (gridIn->GetDataPtr("spec_turb"));
    spectumTurbulent = 0;
}

void TurbSolverStr::InitFlowAsRestart()
{
    Param_TurbSolverStruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();
    if (nTurbulentEquation == 0)
    {
        return;
    }

    StructGrid *grid = StructGridCast(GetGrid(0));

    //int outIterStep = 0;
    //GlobalDataBase::UpdateData("outnstep", &outIterStep, PHINT, 1);

    RDouble   *freeStreamTurbVar =  parameters->GetFreeStreamTurbVar();
    RDouble4D &qTurbulent = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb"));

    RDouble3D &viscousTurbulent = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("vist"));

    RDouble freeStreamViscosity = parameters->GetFreeStreamViscosity();

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 2);

    Range M(0, nTurbulentEquation - 1);
    int mst = M.first();
    int med = M.last();

    for (int m = mst; m <= med; ++ m)
    {
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    qTurbulent(i, j, k, m) = freeStreamTurbVar[m];
                }
            }
        }
    }

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                viscousTurbulent(i, j, k) = freeStreamViscosity;
            }
        }
    }

    int isUnsteady = parameters->GetIsUnsteady();
    if (isUnsteady)
    {
        RDouble4D &qTurbUnsteadyn1 = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb_unsteady_n1"));
        RDouble4D &qTurbUnsteadyn2 = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb_unsteady_n2"));

        qTurbUnsteadyn1 = qTurbulent;
        qTurbUnsteadyn2 = qTurbUnsteadyn1;

        RDouble4D &residualTurbulent      = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_turb"));
        RDouble4D &residualTurbUnsteadyn1 = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_turb_unsteady_n1"));
        RDouble4D &residualTurbUnsteadyn2 = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_turb_unsteady_n2"));

        residualTurbulent      = 0.0;
        residualTurbUnsteadyn1 = 0.0;
        residualTurbUnsteadyn2 = 0.0;

        //! Store the historical sum of InviscidFlux, ViscousFlux, et al., exclusive of DualTimeSourceFlux, used in UpdateUnsteadyFlow for the Rn in the dualtime method.
        RDouble4D &residualTurbUnsteadyTemporary = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_turb_unsteady_tmp"));

        residualTurbUnsteadyTemporary = 0.0;
    }
}

void TurbSolverStr::FillField(Grid *gridIn, FieldProxy *fieldProxy, RDouble value)
{
    RDouble4D &field = fieldProxy->GetField_STR();

    field = value;
}

void TurbSolverStr::FillField(Grid *gridIn, FieldProxy *fieldProxyTarget, FieldProxy *fieldProxySource)
{
    RDouble4D &fieldTarget = fieldProxyTarget->GetField_STR();
    RDouble4D &fieldSource = fieldProxySource->GetField_STR();

    fieldTarget = fieldSource;
}

FieldProxy *TurbSolverStr::CreateFieldProxy(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Param_TurbSolverStruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);
    Range M(0, nTurbulentEquation - 1);

    RDouble4D *field = new RDouble4D(I, J ,K, M, fortranArray);

    FieldProxy *fieldProxy = new FieldProxy();
    fieldProxy->SetField_STR(field, true);
    return fieldProxy;
}

FieldProxy *TurbSolverStr::CreateFlowFieldProxy(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    int numberOfLaminar = GlobalDataBase::GetIntParaFromDB("nl");
    int numberOfChemical = GlobalDataBase::GetIntParaFromDB("nchem");
    int numberOfEquation = numberOfLaminar + numberOfChemical;

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);
    Range M(0, numberOfEquation - 1);

    RDouble4D *field = new RDouble4D(I, J, K, M, fortranArray);

    FieldProxy *fieldProxy = new FieldProxy();
    fieldProxy->SetField_STR(field, true);
    return fieldProxy;
}

FieldProxy *TurbSolverStr::GetFieldProxy(Grid *gridIn, const string &fieldName)
{
    Grid *grid = gridIn;

    RDouble4D *field = reinterpret_cast<RDouble4D *> (grid->GetDataPtr(fieldName));

    FieldProxy *fieldProxy = new FieldProxy();

    fieldProxy->SetField_STR(field);

    return fieldProxy;
}

FieldProxy *TurbSolverStr::GetResidualProxy(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);

    RDouble4D &residualTurbulent = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_turb"));

    FieldProxy *residualProxy = new FieldProxy();

    residualProxy->SetField_STR(&residualTurbulent);

    return residualProxy;
}

void TurbSolverStr::StoreRhsByResidual(Grid *gridIn, FieldProxy *rightHandSideProxy)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Param_TurbSolverStruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    RDouble4D &residualTurbulent = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_turb"));

    RDouble4D &rightHandSide = rightHandSideProxy->GetField_STR();

    rightHandSide = residualTurbulent;
}

void TurbSolverStr::InitResidual(Grid *gridIn, FieldProxy *rightHandSideProxy)
{
    StructGrid *grid = StructGridCast(gridIn);

    Param_TurbSolverStruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 2);

	Range M(0, nTurbulentEquation - 1);
    int mst = M.first();
    int med = M.last();

    RDouble4D &residualTurbulent = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_turb"));
    RDouble4D &rightHandSide = rightHandSideProxy->GetField_STR();

    for (int m = mst; m <= med; ++ m)
    {
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    RDouble &resTurb = residualTurbulent(i, j, k, m);
                    resTurb = -rightHandSide(i, j, k, m);
                }
            }
        }
    }
}

void TurbSolverStr::RecoverResidual(Grid *gridIn, FieldProxy *rightHandSideProxy)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    RDouble4D &residualTurbulent  = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_turb"));
    RDouble4D &rightHandSide = rightHandSideProxy->GetField_STR();

    if (grid->GetLevel() != 0) 
    {
        residualTurbulent = rightHandSide;
    }
}

LIB_EXPORT void TurbSolverStr::InitControlParameters()
{
    FreeControlParameters();
    controlParameters = new Param_TurbSolverStruct();
    controlParameters->Init();
}

LIB_EXPORT Param_TurbSolverStruct *TurbSolverStr::GetControlParameters()
{
    return static_cast <Param_TurbSolverStruct *> (controlParameters);
}

void TurbSolverStr::LoadQ(Grid *gridIn,  FieldProxy *qTrubProxy)
{
    RDouble4D &qTurbulent = *reinterpret_cast<RDouble4D *> (gridIn->GetDataPtr("q_turb"));
    RDouble4D &qTurbInProxy = qTrubProxy->GetField_STR();

    qTurbInProxy = qTurbulent;
}

void TurbSolverStr::RungeKuttaResidual(Grid *gridIn, FieldProxy *dqTurbProxy, RDouble RKCoeff)
{
    StructGrid *grid = StructGridCast(gridIn);

    RDouble4D &dqTurb = dqTurbProxy->GetField_STR();

    Param_TurbSolverStruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    RDouble3D &dt = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("dt"));

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    Range M(0, nTurbulentEquation - 1);
    int mst = M.first();
    int med = M.last();

    for ( int m = mst; m <= med; ++ m )
    {
        for ( int k = kst; k <= ked; ++ k )
        {
            for ( int j = jst; j <= jed; ++ j )
            {
                for ( int i = ist; i <= ied; ++ i )
                {
                    RDouble &dq = dqTurb(i, j, k, m);
                    dq *= dt(i, j, k) * RKCoeff;
                }
            }
        }
    }
}

//! Load rhs to residual res stored in grid.
void TurbSolverStr::LoadResiduals(Grid *gridIn, FieldProxy *rightHandSideProxy)
{
    StructGrid *grid = StructGridCast(gridIn);

    RDouble4D &rightHandSide = rightHandSideProxy->GetField_STR();

    RDouble4D &residualTurbulent = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_turb"));
    
    Param_TurbSolverStruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 2);

    Range M(0, nTurbulentEquation - 1);
    int mst = M.first();
    int med = M.last();

    for (int m = mst; m <= med; ++ m)
    {
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    RDouble &resTurb = residualTurbulent(i, j, k, m);
                    resTurb = - rightHandSide(i, j, k, m);
                }
            }
        }
    }
}

RDouble TurbSolverStr::UnsteadyConvergence(Grid *gridIn)
{
    Param_TurbSolverStruct *parameters = GetControlParameters();
    int isUnsteady = parameters->GetIsUnsteady();
    if (isUnsteady == 0)
    {
        return zero;
    }

    StructGrid *grid = StructGridCast(gridIn);

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    RDouble4D &qLaminar        = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));
    RDouble4D &qTurbulent      = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb"));
    RDouble4D &qTurbUnsteadyn1 = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb_unsteady_n1"));
    RDouble4D &qTurbUnsteadyn2 = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb_unsteady_n2"));
    RDouble4D &residualTurbulent = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_turb"));

    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    RDouble *qTurb0 = new RDouble[nTurbulentEquation];
    RDouble *qTurb1 = new RDouble[nTurbulentEquation];
    RDouble *qTurb2 = new RDouble[nTurbulentEquation];

    RDouble *qTurbConserve0 = new RDouble[nTurbulentEquation];
    RDouble *qTurbConserve1 = new RDouble[nTurbulentEquation];
    RDouble *qTurbConserve2 = new RDouble[nTurbulentEquation];

    using namespace IDX;

    RDouble summedPrimitive = zero;
    RDouble summedConserve  = zero;

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                for (int m = 0; m < nTurbulentEquation; ++ m)
                {
                    qTurb0[m] = qTurbulent(i, j, k, m);
                    qTurb1[m] = qTurbUnsteadyn1(i, j, k, m);
                    qTurb2[m] = qTurbUnsteadyn2(i, j, k, m);
                }

                RDouble rho = 1.0;

                if (nTurbulentEquation >= 2)
                {
                    rho = qLaminar(i, j, k, IR);
                }

                for (int m = 0; m < nTurbulentEquation; ++ m)
                {
                    qTurbConserve0[m] = rho * qTurb0[m];
                    qTurbConserve1[m] = rho * qTurb1[m];
                    qTurbConserve2[m] = rho * qTurb2[m];
                }

                for (int m = 0; m < nTurbulentEquation; ++ m)
                {
                    RDouble dqTurbPseudo = residualTurbulent(i, j, k, m);                   // qn+1,p+1 - qn+1,p
                    RDouble dqTurbPhysical = qTurbConserve0[m] - qTurbConserve1[m];  // qn+1,p+1 - qn
                    summedPrimitive += dqTurbPseudo * dqTurbPseudo;
                    summedConserve  += dqTurbPhysical * dqTurbPhysical;
                }
            }
        }
    }

    delete []qTurb0;
    delete []qTurb1;
    delete []qTurb2;
    delete []qTurbConserve0;
    delete []qTurbConserve1;
    delete []qTurbConserve2;

    RDouble unsteadyConvergence = sqrt(ABS(summedPrimitive / (summedConserve + SMALL)));
    return unsteadyConvergence;
}

void TurbSolverStr::UpdateUnsteadyFlow(Grid *gridIn)
{
    Param_TurbSolverStruct *parameters = GetControlParameters();
    
    int isUnsteady = parameters->GetIsUnsteady();
    if (isUnsteady == 0)
    {
        return;
    }

    StructGrid *grid = StructGridCast(gridIn);

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 2);

    //! The volume at current timestep. If it is dualtime step method, it represents the volume at timestep of n+1.
    RDouble3D &volume           = grid->GetVolume(0);
    //! The volume at timestep of n  .
    RDouble3D &volumeUnsteadyn1 = grid->GetVolume(1);
    //! The volume at timestep of n-1.
    RDouble3D &volumeUnsteadyn2 = grid->GetVolume(2);

    RDouble4D &qTurbulent = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb"));
    RDouble4D &qTurbUnsteadyn1 = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb_unsteady_n1"));
    RDouble4D &qTurbUnsteadyn2 = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb_unsteady_n2"));

    RDouble4D &residualTurbUnsteadyn1        = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_turb_unsteady_n1" ));
    RDouble4D &residualTurbUnsteadyn2        = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_turb_unsteady_n2" ));
    RDouble4D &residualTurbUnsteadyTemporary = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_turb_unsteady_tmp"));

    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    qTurbUnsteadyn2 = qTurbUnsteadyn1;
    qTurbUnsteadyn1 = qTurbulent  ;

    volumeUnsteadyn2 = volumeUnsteadyn1;
    volumeUnsteadyn1 = volume  ;

    //! Here the current outnstep is over, the value of the stored resTmp should be assigned to resn1 for the next outnstep.
    //! It should be noticed that residualUnsteadyTemporary only contain the inviscid and viscous flux.
    residualTurbUnsteadyn2 = residualTurbUnsteadyn1;
    residualTurbUnsteadyn1 = residualTurbUnsteadyTemporary;

    //! Statistical variables for DES simulation.
    int ifStaticsFlowField = parameters->GetIfStaticsFlowField();
    if (ifStaticsFlowField > 0)
    {
        int outIterStep        = GlobalDataBase::GetIntParaFromDB("outnstep");
        int startStatisticStep = GlobalDataBase::GetIntParaFromDB("startStatisticStep");

        if (outIterStep >= startStatisticStep)
        {
            int nStatisticalStep = GlobalDataBase::GetIntParaFromDB("nStatisticalStep");

            RDouble coeff1 = (nStatisticalStep - 1.0) / (nStatisticalStep * 1.0);
            RDouble coeff2 =                     1.0  / (nStatisticalStep * 1.0);

            RDouble4D &qAverageTurb = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("qAverage_turb"));

            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    for (int i = ist; i <= ied; ++ i)
                    {
                        for (int m = 0; m < nTurbulentEquation; ++ m)
                        {
                            qAverageTurb(i, j, k, m) = coeff1 * qAverageTurb(i, j, k, m) + coeff2 * qTurbulent(i, j, k, m);
                        }
                    }
                }
            }
        }
    }
}

void TurbSolverStr::DualTimeSource(Grid *gridIn)
{
    Param_TurbSolverStruct *parameters = GetControlParameters();
    int isUnsteady = parameters->GetIsUnsteady();
    if (isUnsteady == 0)
    {
        return;
    }

    StructGrid *grid = StructGridCast(gridIn);

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    RDouble4D &qLaminar = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));

    RDouble4D &qTurbulent      = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb"));
    RDouble4D &qTurbUnsteadyn1 = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb_unsteady_n1"));
    RDouble4D &qTurbUnsteadyn2 = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb_unsteady_n2"));

    RDouble4D &residualTurbulent             = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_turb"));
    RDouble4D &residualTurbUnsteadyn1        = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_turb_unsteady_n1" ));
    RDouble4D &residualTurbUnsteadyn2        = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_turb_unsteady_n2" ));
    RDouble4D &residualTurbUnsteadyTemporary = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_turb_unsteady_tmp"));

    //! The volume at current timestep. If it is dualtime step method, it represents the volume at timestep of n+1.
    RDouble3D &volume           = grid->GetVolume(0);
    //! The volume at timestep of n  .
    RDouble3D &volumeUnsteadyn1 = grid->GetVolume(1);
    //! the volume at timestep of n-1.
    RDouble3D &volumeUnsteadyn2 = grid->GetVolume(2);

    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    //! The primitive variables.
    RDouble *qTurbPrimtive0 = new RDouble[nTurbulentEquation];
    RDouble *qTurbPrimtive1 = new RDouble[nTurbulentEquation];
    RDouble *qTurbPrimtive2 = new RDouble[nTurbulentEquation];

    //! The conservative variables.
    RDouble *qTurbConserve0 = new RDouble[nTurbulentEquation];
    RDouble *qTurbConserve1 = new RDouble[nTurbulentEquation];
    RDouble *qTurbConserve2 = new RDouble[nTurbulentEquation];

    //! Store the historical sum of InviscidFlux, ViscousFlux, et al., exclusive of DualTimeSourceFlux,
    //! used in UpdateUnsteadyFlow for the Rn in the dualtime method.
    residualTurbUnsteadyTemporary = residualTurbulent;

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

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                for (int m = 0; m < nTurbulentEquation; ++ m)
                {
                    qTurbPrimtive0[m] = qTurbulent          (i, j, k, m);
                    qTurbPrimtive1[m] = qTurbUnsteadyn1(i, j, k, m);
                    qTurbPrimtive2[m] = qTurbUnsteadyn2(i, j, k, m);
                }

                RDouble roe = 1.0;

                if (nTurbulentEquation >= 2)
                {
                    roe  = qLaminar(i, j, k, IR);
                }

                for (int m = 0; m < nTurbulentEquation; ++ m)
                {
                    qTurbConserve0[m] = roe * qTurbPrimtive0[m];
                    qTurbConserve1[m] = roe * qTurbPrimtive1[m];
                    qTurbConserve2[m] = roe * qTurbPrimtive2[m];
                }

                for (int m = 0; m < nTurbulentEquation; ++ m)
                {
                    RDouble dualSrcRes = dualTimeResC1 * residualTurbulent(i, j, k, m)           +
                                         dualTimeResC2 * residualTurbUnsteadyn1(i, j, k, m) +
                                         dualTimeResC3 * residualTurbUnsteadyn2(i, j, k, m);

                    RDouble dualSrcQ  = dualTimeQC1   * qTurbConserve0[m] * volume(i, j, k)           +
                                        dualTimeQC2   * qTurbConserve1[m] * volumeUnsteadyn1(i, j, k) +
                                        dualTimeQC3   * qTurbConserve2[m] * volumeUnsteadyn2(i, j, k);

                    residualTurbulent(i, j ,k, m) = dualSrcRes + dualSrcQ;
                }
            }
        }
    }

    delete []qTurbPrimtive0;
    delete []qTurbPrimtive1;
    delete []qTurbPrimtive2;
    delete []qTurbConserve0;
    delete []qTurbConserve1;
    delete []qTurbConserve2;
}

RDouble TurbSolverStr::GetDistance(Grid *gridIn, int i, int j, int k)
{
    StructGrid *grid = StructGridCast(gridIn);

    RDouble3D &xx = *grid->GetStructX();
    RDouble3D &yy = *grid->GetStructY();
    RDouble3D &zz = *grid->GetStructZ();

    RDouble dx = xx(i + 1, j, k) - xx(i, j, k);
    RDouble dy = yy(i + 1, j, k) - yy(i, j, k);
    RDouble dz = zz(i + 1, j, k) - zz(i, j, k);

    RDouble dsDirI = (dx * dx + dy * dy + dz * dz);

    dx = xx(i, j + 1, k) - xx(i, j, k);
    dy = yy(i, j + 1, k) - yy(i, j, k);
    dz = zz(i, j + 1, k) - zz(i, j, k);

    RDouble dsDirJ = (dx * dx + dy * dy + dz * dz);

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    RDouble dsDirK;
    if (nk == 1)
    {
        dsDirK = 0;
    }
    else
    {
        dx = xx(i, j, k + 1) - xx(i, j, k);
        dy = yy(i, j, k + 1) - yy(i, j, k);
        dz = zz(i, j, k + 1) - zz(i, j, k);
        dsDirK = (dx * dx + dy * dy + dz * dz);
    }

    return sqrt(MAX(MAX(dsDirI, dsDirJ), dsDirK));

}

void TurbSolverStr::SourceFluxOneEquation(Grid *grid)
{
    Param_TurbSolverStruct *parameters = GetControlParameters();
    //int kindOfTurbSource = parameters->GetKindOfTurbSource();

    ComputeLengthScaleofOneEquationModel(grid);

    //if (turbSourceKind == 0)
    //{
        SourceFluxOneEquationOriginal(grid);
    //}
    //else if (turbSourceKind == 1)
    //{
        //SourceFluxOneEquationNew(grid);
    //}
    //else
    //{
    //    TK_Exit::ExceptionExit("No this turbulence Source Kind! \n");
    //}
}

void TurbSolverStr::SourceFluxOneEquationOriginal(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    RDouble4D &qLaminar          = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));
    RDouble4D &qTurbulent        = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb"  ));
    RDouble4D &spectumTurbulent  = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("spec_turb"));
    RDouble4D &residualTurbulent = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_turb"));
    RDouble3D &viscousLaminar    = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("visl"));

    RDouble3D *lengthScale = this->GetLengthScale(gridIn);

    RDouble3D &walldist = *grid->GetWallDist();
    RDouble3D &volume   = *(grid->GetCellVolume() );

    Param_TurbSolverStruct *parameters = GetControlParameters();
    RDouble oRefReNumber = parameters->GetoRefReNumber();

    //! von karman constant.
    const RDouble karm   = 0.41;
    const RDouble cb1    = 0.1355;
    const RDouble cb2    = 0.622;
    const RDouble cv1    = 7.1;
    const RDouble cw2    = 0.3;
    const RDouble cw3    = 2.0;
    const RDouble sigma  = 2.0/3.0;
    const RDouble cdes   = - 1.0;
    const RDouble sdilim = 1.0e9;

    RDouble oSigma = one / sigma;
    RDouble oKarmQuadratic = one / (karm * karm);
    RDouble cw1    = cb1 * oKarmQuadratic + (one + cb2) / sigma;
    RDouble cv13   = cv1 * cv1 * cv1;
    RDouble cw32   = cw3 * cw3;
    RDouble cw36   = cw32 * cw32 * cw32;
    RDouble or6    = one / six;

    RDouble ct3 = 1.2;
    RDouble ct4 = 0.5;

    int SAProductType = parameters->GetSAProductType();

    using namespace IDX;

    Range II, JJ, KK;
    GetRange(ni, nj, nk, -2, 1, II, JJ, KK);

    RDouble4D &gradUVWTCellCenterX = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterX"));
    RDouble4D &gradUVWTCellCenterY = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterY"));
    RDouble4D &gradUVWTCellCenterZ = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterZ"));

    RDouble4D &gradTurbulentCellCenterX = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradTurbulentCellCenterX"));
    RDouble4D &gradTurbulentCellCenterY = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradTurbulentCellCenterY"));
    RDouble4D &gradTurbulentCellCenterZ = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradTurbulentCellCenterZ"));

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble dudx = gradUVWTCellCenterX(i, j, k, 0);
                RDouble dudy = gradUVWTCellCenterY(i, j, k, 0);
                RDouble dudz = gradUVWTCellCenterZ(i, j, k, 0);

                RDouble dvdx = gradUVWTCellCenterX(i, j, k, 1);
                RDouble dvdy = gradUVWTCellCenterY(i, j, k, 1);
                RDouble dvdz = gradUVWTCellCenterZ(i, j, k, 1);

                RDouble dwdx = gradUVWTCellCenterX(i, j, k, 2);
                RDouble dwdy = gradUVWTCellCenterY(i, j, k, 2);
                RDouble dwdz = gradUVWTCellCenterZ(i, j, k, 2);

                RDouble nu      = qTurbulent(i, j, k, ISA);
                RDouble vol     = volume(i, j, k);
                RDouble density = qLaminar(i, j, k, IR);
                RDouble vislam  = viscousLaminar(i, j, k);
                RDouble dn      = (*lengthScale)(i, j, k);

                RDouble work1 = dudx + dvdy + dwdz;
                RDouble work2 = dudx * dudx + dvdy * dvdy + dwdz * dwdz;
                RDouble work3 = (dudy + dvdx) * (dudy + dvdx) + 
                                (dvdz + dwdy) * (dvdz + dwdy) +
                                (dwdx + dudz) * (dwdx + dudz);
                RDouble work4 = sqrt((dwdy - dvdz) * (dwdy - dvdz) +
                                     (dudz - dwdx) * (dudz - dwdx) +
                                     (dvdx - dudy) * (dvdx - dudy));

                RDouble stp, str;
                if (SAProductType == 1)
                {
                    //! Classical production (instead of rotational).
                    stp   = two * work2 + work3 - two3rd * work1 * work1;
                    work4 = sqrt(ABS(stp));
                }
                else if (SAProductType == 4)
                {
                    //! Correction on the rotational production term for laminar vortex.
                    stp   = two * work2 + work3 - two3rd * work1 * work1;
                    str   = work4;
                    work4 = str - two * MAX(zero,str - sqrt(ABS(stp)));
                }
                else if (SAProductType == 5)
                {
                    //! Higher correction on the rotational production term for laminar vortex.
                    stp   = two * work2 + work3 - two3rd * work1 * work1;
                    str   = work4;
                    stp   = ABS(stp);
                    work4 = str * (one - MAX(zero, (str * str - stp) / (str * str + stp + SMALL)));
                }

                str = work4;

                RDouble odnQuadratic = one / (dn * dn);
                RDouble kyQuadratic  = oKarmQuadratic * odnQuadratic;
                RDouble oLam = one / (vislam / density);
                RDouble ld   = MAX (nu * oLam, 1E-20);
                RDouble ld3  = POWER3(ld);
                RDouble fv1  = ld3 / (ld3 + cv13);
                RDouble fv2  = one - ld / (one + ld * fv1);
                RDouble std  = MAX(str + nu * fv2 * kyQuadratic * oRefReNumber, 1.E-20);
                RDouble product = cb1 * std * nu;
                RDouble dproduct_dnu = cb1 * str;

                RDouble oostd = one / std;
                RDouble nuQuadratic = nu * nu;
                RDouble r = MIN(nu * oostd * kyQuadratic * oRefReNumber, 10.0);
                RDouble r3 = POWER3(r);
                RDouble r6 = r3 * r3;
                RDouble g = r + cw2 * (r6 - r);
                RDouble g3 = g * g * g;
                RDouble g6   = g3 * g3;
                g6   = MAX (g6, 1E-10);
                RDouble fw   = g * pow(static_cast<RDouble>((one + cw36) / (g6 + cw36)), static_cast<RDouble>(or6));
                RDouble destruc = cw1 * fw * nuQuadratic * odnQuadratic * oRefReNumber;

                RDouble ddestruc_dnu = two * (cw1 * fw - cb1 * oKarmQuadratic * fv2) * odnQuadratic * nu * oRefReNumber;
                ddestruc_dnu = abs (ddestruc_dnu);

                RDouble dkedx = gradTurbulentCellCenterX(i, j, k, ISA);
                RDouble dkedy = gradTurbulentCellCenterY(i, j, k, ISA);
                RDouble dkedz = gradTurbulentCellCenterZ(i, j, k, ISA);
                RDouble gradientQuadratic  = dkedx * dkedx + dkedy * dkedy + dkedz * dkedz;
                RDouble sdi = oSigma * cb2 * gradientQuadratic * oRefReNumber;
                RDouble dsdi_dnu = 0.0;

                RDouble &resTurb = residualTurbulent(i, j, k, ISA);
                RDouble &spectTurb = spectumTurbulent(i, j, k, ISA);
                resTurb   +=  (product - destruc + sdi) * vol;
                spectTurb +=  (dproduct_dnu + ddestruc_dnu + dsdi_dnu) * vol;
            }
        }
    }
}

//! Hekun test by using chant SA  source model.
void TurbSolverStr::SourceFluxOneEquationNew(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    RDouble4D &qLaminar            = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));
    RDouble4D &qTurbulent          = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb"  ));
    RDouble4D &spectumTurbulent    = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("spec_turb"));
    RDouble4D &residualTurbulent   = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_turb"));
    RDouble3D &viscousLaminar      = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("visl"));

    RDouble3D *lengthScale = this->GetLengthScale(gridIn);

    RDouble3D &walldist = * grid->GetWallDist();
    RDouble3D &volume   = *(grid->GetCellVolume() );

    Param_TurbSolverStruct *parameters = GetControlParameters();
    RDouble refReNumber = parameters->GetRefReNumber();
    RDouble oRefReNumber = parameters->GetoRefReNumber();

    //! von karman constant.
    const RDouble karm   = 0.41;
    const RDouble cb1    = 0.1355;
    const RDouble cb2    = 0.622;
    const RDouble cv1    = 7.1;
    const RDouble cw2    = 0.3;
    const RDouble cw3    = 2.0;
    const RDouble sigma  = 2.0/3.0;
    const RDouble sdilim = 1.0e9;

    RDouble osigma = 1.0/sigma;
    RDouble okarm2 = one / (karm * karm);
    RDouble cw1    = cb1 * okarm2 + (one + cb2) * osigma;
    RDouble cv13   = cv1 * cv1 * cv1;
    RDouble cw32   = cw3 * cw3;
    RDouble cw36   = cw32 * cw32 * cw32;
    RDouble or6    = one / six;
    RDouble cv2    = 5.0;

    RDouble ct3 = 1.2;
    RDouble ct4 = 0.5;

    using namespace IDX;

    Range II, JJ, KK;
    GetRange(ni, nj, nk, -2, 1, II, JJ, KK);

    RDouble4D &gradUVWTCellCenterX = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterX"));
    RDouble4D &gradUVWTCellCenterY = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterY"));
    RDouble4D &gradUVWTCellCenterZ = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterZ"));

    RDouble4D &gradTurbulentCellCenterX = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradTurbulentCellCenterX"));
    RDouble4D &gradTurbulentCellCenterY = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradTurbulentCellCenterY"));
    RDouble4D &gradTurbulentCellCenterZ = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradTurbulentCellCenterZ"));

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble dudx = gradUVWTCellCenterX(i, j, k, 0);
                RDouble dudy = gradUVWTCellCenterY(i, j, k, 0);
                RDouble dudz = gradUVWTCellCenterZ(i, j, k, 0);

                RDouble dvdx = gradUVWTCellCenterX(i, j, k, 1);
                RDouble dvdy = gradUVWTCellCenterY(i, j, k, 1);
                RDouble dvdz = gradUVWTCellCenterZ(i, j, k, 1);

                RDouble dwdx = gradUVWTCellCenterX(i, j, k, 2);
                RDouble dwdy = gradUVWTCellCenterY(i, j, k, 2);
                RDouble dwdz = gradUVWTCellCenterZ(i, j, k, 2);

                RDouble work1 = dudx + dvdy + dwdz;
                RDouble work2 = dudx * dudx + dvdy * dvdy + dwdz * dwdz;
                RDouble work3 = (dudy + dvdx) * (dudy + dvdx) +
                                (dvdz + dwdy) * (dvdz + dwdy) +
                                (dwdx + dudz) * (dwdx + dudz);
                RDouble work4 = sqrt((dwdy - dvdz) * (dwdy - dvdz) +
                                     (dudz - dwdx) * (dudz - dwdx) +
                                     (dvdx - dudy) * (dvdx - dudy));

                RDouble stp   = two * work2 + work3 - two3rd * work1 * work1;
                RDouble str   = work4;
                work4 = str - two * MAX(zero,str - sqrt(ABS(stp)));
                str    = work4;
                RDouble vtbar  = qTurbulent(i, j, k, ISA);
                RDouble vol    = volume(i, j, k);

                RDouble vsl    = viscousLaminar(i, j, k);
                RDouble rho    = qLaminar(i, j, k,IR);
                RDouble romu   = rho * vtbar;
                RDouble olam   = rho / (vsl + SMALL);

                RDouble dn = (*lengthScale)(i, j, k);

                RDouble dn2 = dn * dn;
                RDouble od2 = one / dn2;

                //! Spalart mod. Use limited version of ld.
                //! All the functions are well-behaved as ld -> 0.
                RDouble ld  = MAX(olam * vtbar, SMALL);
                RDouble ld3  = POWER3(ld);
                //! Spalart mod. New damping functions so Stilde is always positive.
                RDouble fv1 = ld3 / (ld3 + cv13);
                RDouble fvt = 1.0 / (1.0 + ld / cv2);
                RDouble fv2 = POWER3(fvt);
                RDouble fv3 = (1 + ld * fv1) / ld * (1 - fv2);
                RDouble omg = str;

                RDouble ki2      = karm * karm;
                RDouble reki2dn2 = ki2 * dn2 * refReNumber;
                RDouble sbar     = fv3 * omg + vtbar / reki2dn2 * fv2;
                sbar             = MAX(sbar, SMALL);
                RDouble std      = sbar;
                RDouble oostd    = SIGN(one,std) / (ABS(std) + SMALL);
                RDouble rbar     = vtbar * oostd / reki2dn2;
                rbar             = MIN(rbar, 10.0);
                RDouble rbar2    = rbar * rbar;
                RDouble rbar5    = rbar2 * rbar2 * rbar;
                RDouble gbar     = MAX(rbar * (1.0 + cw2 * (rbar5 - 1.0)), SMALL);
                RDouble g3       = POWER3(gbar);
                RDouble g6       = g3 * g3;
                g6       = MAX (g6, 1E-20);
                RDouble fw       = pow(static_cast<RDouble>((1.0 + cw36) / (1.0 + pow(cw3 / gbar, 6.0))), 1.0 / 6.0);
                RDouble ft2      = 0.0;
                RDouble cpmu     = cb1 * (1.0 - ft2) * str;
                RDouble cdmu     = (cw1 * fw - cb1 * fv2 / ki2) * vtbar * od2 * oRefReNumber;

                //! Calculate "negative" production term in order to assure a zero.
                //! Solution in the laminar region.

                RDouble ld2  = ld * ld;
                RDouble srct = cpmu - cdmu;

                RDouble dkedx = gradTurbulentCellCenterX(i, j, k, ISA);
                RDouble dkedy = gradTurbulentCellCenterY(i, j, k, ISA);
                RDouble dkedz = gradTurbulentCellCenterZ(i, j, k, ISA);

                //! Non linear diffusion term.
                RDouble grd2 = dkedx * dkedx + dkedy * dkedy + dkedz * dkedz;
                RDouble sdi  = osigma * cb2 * grd2 * oRefReNumber;

                RDouble &resTurb = residualTurbulent(i, j, k, ISA);
                RDouble &spectTurb = spectumTurbulent(i, j, k, ISA);
                resTurb +=  (srct * vtbar + sdi) * vol;
                spectTurb  +=  (abs(cpmu) + 2.0 * abs(cdmu)) * vol;
            }
        }
    }
}

void TurbSolverStr::Crossing(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);

    RDouble4D &qTurbulent = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb"));
    RDouble3D &cross = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("cross"));

    RDouble4D &faceNormalComponentX = *(grid->GetFaceNormalX());
    RDouble4D &faceNormalComponentY = *(grid->GetFaceNormalY());
    RDouble4D &faceNormalComponentZ = *(grid->GetFaceNormalZ());
    RDouble4D &area   = *(grid->GetFaceArea()   );
    RDouble3D &volume = *(grid->GetCellVolume() );

    GreenGeo3D greenGeo3D;

    int ist, ied, jst, jed, kst, ked;
	grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    int ndim = GetDim();
    int js = 1;
    int ks = 1;
    if (ndim == 2)
    {
        ks = 0;
    }
    if (ndim == 1)
    {
        js = 0;
    }

    RDouble dkedx,dkedy,dkedz,dkwdx,dkwdy,dkwdz;
    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                Get_SPM(i, j, k, js, ks, ndim, &greenGeo3D, faceNormalComponentX, faceNormalComponentY, faceNormalComponentZ, area, volume);
                DXDYDZ(qTurbulent, &greenGeo3D, i, j, k, js, ks, IDX::IKE, dkedx, dkedy, dkedz);
                DXDYDZ(qTurbulent, &greenGeo3D, i, j, k, js, ks, IDX::IKW, dkwdx, dkwdy, dkwdz);
                RDouble &crossCell = cross(i, j, k);
                crossCell = dkedx * dkwdx + dkedy * dkwdy + dkedz * dkwdz;
            }
        }
    }
}

void TurbSolverStr::Blending(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    RDouble4D &qLaminar       = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q" ));
    RDouble4D &qTurbulent     = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb" ));
    RDouble3D &viscousLaminar = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("visl"));

    RDouble3D &cross = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("cross"));
    RDouble3D &blend = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("blend"));

    RDouble3D &walldist = * StructGridCast(grid)->GetWallDist();

    Param_TurbSolverStruct *parameters = GetControlParameters();
    RDouble refReNumber  = parameters->GetRefReNumber();
    RDouble oRefReNumber = parameters->GetoRefReNumber();
    RDouble oRefReNumberSqr = oRefReNumber * oRefReNumber;

    RDouble betas = 0.09;
    RDouble sigw2 = 0.856;

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    using namespace IDX;

    //Compute maximum cross diffusion term across flowfield
    RDouble cdkwmax = 0.0;
    RDouble cdkwmin = 1.0e-10;

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble rho  = qLaminar(i, j, k,IR);
                RDouble ke   = qTurbulent(i, j, k, IKE);
                RDouble kw   = qTurbulent(i, j, k, IKW);
                RDouble crss = cross(i, j, k);
 
                RDouble miuLam = viscousLaminar(i, j, k);
                RDouble dist   = walldist(i, j, k);
                RDouble dist2  = dist * dist;

                //calculate cd_kw
                RDouble crossDiff = 2.0 * rho * sigw2 * crss / (kw + SMALL);

                //! Original Menter CD_kw calculation.
                RDouble cdkw = MAX(crossDiff, cdkwmin);

                //! Calculate arg1.
                RDouble term1 = sqrt(ABS(ke)) * oRefReNumber /(betas * dist * kw + SMALL);
                RDouble term2 = 500.0 * miuLam * oRefReNumberSqr / (rho * kw * dist2 + SMALL);
                RDouble term3 = MAX(term1, term2);
                RDouble term4 = 4.0 * rho * sigw2 * ke / (cdkw * dist2 + SMALL);
                RDouble arg1  = MIN(term3, term4);
     
                //! Compute the blending function fbsl.
                RDouble &blendCell = blend(i, j, k);
                blendCell = tanh(static_cast<RDouble>(arg1 * arg1 * arg1 * arg1));
            }
        }
    }
    grid->GhostCell3DExceptInterface(blend);
    //GhostCell3D(blend, ni, nj, nk);
    //PHSPACE::CommunicateInterfaceValue(grid, &blend, "blend");
}

void TurbSolverStr::SourceFluxTwoEquation(Grid *gridIn)
{
    //! This routine calculates the source terms relative to k and tau,
    //! epsilon or omega equations for (k-tau), (k-epsilon) or (k-omega)
    //! two-equation models.

    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    RDouble4D &qPrimitive        = *reinterpret_cast< RDouble4D * > (grid->GetDataPtr("q" ));
    RDouble3D &viscousLaminar    = *reinterpret_cast< RDouble3D * > (grid->GetDataPtr("visl"));
    RDouble3D &viscousTurbulent  = *reinterpret_cast< RDouble3D * > (grid->GetDataPtr("vist"));
    RDouble4D &qTurbulent        = *reinterpret_cast< RDouble4D * > (grid->GetDataPtr("q_turb" ));
    RDouble4D &spectumTurbulent  = *reinterpret_cast< RDouble4D * > (grid->GetDataPtr("spec_turb"));
    RDouble4D &residualTurbulent = *reinterpret_cast< RDouble4D * > (grid->GetDataPtr("res_turb"));

    RDouble3D &walldist = * grid->GetWallDist();

    RDouble4D &xfn  = *(grid->GetFaceNormalX());
    RDouble4D &yfn  = *(grid->GetFaceNormalY());
    RDouble4D &zfn  = *(grid->GetFaceNormalZ());
    RDouble4D &area = *(grid->GetFaceArea()   );
    RDouble3D &vol  = *(grid->GetCellVolume() );

    Param_TurbSolverStruct *parameters = GetControlParameters();
    RDouble refReNumber = parameters->GetRefReNumber();
    RDouble oRefReNumber = parameters->GetoRefReNumber();

    string viscousName = parameters->GetViscousName();

    int neasm        = GlobalDataBase::GetIntParaFromDB("neasm");
    int SSTProductType = parameters->GetSSTProductType();

    //! To determine if DES used or not.
    int isdes = 0;

    const RDouble cdes_ke = 0.61, cdes_kw = 0.78;
    RDouble cdes,dxyz,desdxyz,cdissk,cdissw,tmp,lkw;

    //! To determine if compressible correction used or not.
    int isCC = 0;

    const RDouble a1 = 1.0, a2 = 0.4, a3 = 0.2, mt02 = 0.25*0.25;
    RDouble mt2, c2, f2mt2, a1f2mt2, pdbar, prodk_cc, prodw_cc, dissk_cc, dissw_cc;
    RDouble3D &gama = *reinterpret_cast< RDouble3D * > (grid->GetDataPtr("gama"));

    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();
    int nTotal = nTotalCell + nBoundFace;

    RDouble cblend, rho, ke, kw;
    RDouble sij2, divv, prodk, dissk, prodw, dissw, cdkww, crss, mut;
    RDouble srck, srcw, volume;
    RDouble dudx0, dudy0, dudz0, dvdx0, dvdy0, dvdz0, dwdx0, dwdy0, dwdz0;  
    RDouble s11, s22, s33, s12, s13, s23, w12, w13, w23;  
    RDouble vort2;
    RDouble cprodk, cprodw, ccdkww, source_dk_neg, source_dw_neg;

    RDouble turb_fbetas,turb_fbeta,turb_beta,turb_alphaw;

    Range II, JJ, KK;
    GetRange(ni, nj, nk, -2, 1, II, JJ, KK);

    RDouble4D &gradUVWTCellCenterX = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterX"));
    RDouble4D &gradUVWTCellCenterY = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterY"));
    RDouble4D &gradUVWTCellCenterZ = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterZ"));

    //RDouble3D dudx(II, JJ, KK, fortranArray);
    //RDouble3D dudy(II, JJ, KK, fortranArray);
    //RDouble3D dudz(II, JJ, KK, fortranArray);

    
    //RDouble3D dvdx(II, JJ, KK, fortranArray);
    //RDouble3D dvdy(II, JJ, KK, fortranArray);
    //RDouble3D dvdz(II, JJ, KK, fortranArray);


    //RDouble3D dwdx(II, JJ, KK, fortranArray);
    //RDouble3D dwdy(II, JJ, KK, fortranArray);
    //RDouble3D dwdz(II, JJ, KK, fortranArray);
    
    using namespace IDX;

    //dudx(II, JJ, KK) = gradUVWTCellCenterX(II, JJ, KK, 0);
    //dudy(II, JJ, KK) = gradUVWTCellCenterY(II, JJ, KK, 0);
    //dudz(II, JJ, KK) = gradUVWTCellCenterZ(II, JJ, KK, 0);


    //dvdx(II, JJ, KK) = gradUVWTCellCenterX(II, JJ, KK, 1);
    //dvdy(II, JJ, KK) = gradUVWTCellCenterY(II, JJ, KK, 1);
    //dvdz(II, JJ, KK) = gradUVWTCellCenterZ(II, JJ, KK, 1);


    //dwdx(II, JJ, KK) = gradUVWTCellCenterX(II, JJ, KK, 2);
    //dwdy(II, JJ, KK) = gradUVWTCellCenterY(II, JJ, KK, 2);
    //dwdz(II, JJ, KK) = gradUVWTCellCenterZ(II, JJ, KK, 2);

    RDouble3D &sourceEnergy = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("sengy"));

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    //====================
    //=== K-OMEGA MODELS
    //====================
    if (viscousName.substr(0,6) == "2eq-kw")
    {
        //******************************
        //*** MENTER'S K-OMEGA SST MODEL
        //******************************
        if (viscousName.substr(0,13) == "2eq-kw-menter")
        {
            //turbulence model constants
            RDouble SSTProductLimit = parameters->GetSSTProductLimit();
            GlobalDataBase::GetData("turb_fbetas",&turb_fbetas,PHDOUBLE,1);
            GlobalDataBase::GetData("turb_fbeta" ,&turb_fbeta ,PHDOUBLE,1);
            

            RDouble SST_beta1 = parameters->GetSST_beta1();
            RDouble SST_beta2 = parameters->GetSST_beta2();
            RDouble SST_betaStar = parameters->GetSST_betaStar();
            RDouble KW_sigmaW2 = parameters->GetKW_sigmaW2();

            RDouble SST_alphaw1 = parameters->GetSST_alphaw1();
            RDouble SST_alphaw2 = parameters->GetSST_alphaw2();

            RDouble3D &cross = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("cross"));
            RDouble3D &blend = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("blend"));

            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    for (int i = ist; i <= ied; ++ i)
                    {
                        cblend = blend(i,j,k);
                        crss   = cross(i,j,k);
                        rho    = qPrimitive(i,j,k,IR);
                        ke     = qTurbulent(i,j,k,IKE);
                        kw     = qTurbulent(i,j,k,IKW);
                        mut    = viscousTurbulent(i,j,k);
                        volume = vol(i,j,k);

                        // to determine if compressible correction used or not                        
                        a1f2mt2 = 0.0;
                        
                        if (isCC == 1)
                        {
                            c2      = gama(i, j, k) * qPrimitive(i, j, k, IP) / rho;
                            mt2     = (ke + ke) / c2;
                            f2mt2   = (1.0 - cblend) * MAX(mt2 - mt02, 0.0);
                            a1f2mt2 = a1 * f2mt2;
                        }

                        turb_beta   = cblend * SST_beta1   + (1.0 - cblend) * SST_beta2;
                        turb_alphaw = cblend * SST_alphaw1 + (1.0 - cblend) * SST_alphaw2;

                        dudx0 = gradUVWTCellCenterX(i,j,k,0); dudy0 = gradUVWTCellCenterY(i,j,k,0); dudz0 = gradUVWTCellCenterZ(i,j,k,0);
                        dvdx0 = gradUVWTCellCenterX(i,j,k,1); dvdy0 = gradUVWTCellCenterY(i,j,k,1); dvdz0 = gradUVWTCellCenterZ(i,j,k,1);
                        dwdx0 = gradUVWTCellCenterX(i,j,k,2); dwdy0 = gradUVWTCellCenterY(i,j,k,2); dwdz0 = gradUVWTCellCenterZ(i,j,k,2);	
                        
                        s11 = dudx0; s22 = dvdy0; s33 = dwdz0;
                        s12 = half * (dudy0 + dvdx0); s13 = half * (dudz0 + dwdx0); s23 = half * (dvdz0 + dwdy0);
                        w12 = half * (dudy0 - dvdx0); w13 = half * (dudz0 - dwdx0); w23 = half * (dvdz0 - dwdy0);

                        sij2 = two * (s11*s11 + s22*s22 + s33*s33 + two*(s12*s12 + s13*s13 + s23*s23));  //! Modulus of S.
                        divv  = s11 + s22 + s33;                                   //! Divergence of velocity.

                        prodk = mut * (sij2 - two3rd * divv * divv) * oRefReNumber - two3rd * rho * ke * divv;  //! SST-2003 standard version.

                        if (SSTProductType != 0) //! Produck term based on vorticity.
                        {
                            vort2 = four * (w12 * w12 + w13 * w13 + w23 * w23);           //! Vorticity.
                            prodk = mut * vort2 * oRefReNumber - two3rd * rho * ke * divv;   //! SST-V.
                        }

                        dissk = turb_fbetas * SST_betaStar * rho * ke * kw * refReNumber;  //! Distructive term of TKE.
                        dissw = turb_fbeta * turb_beta * rho * kw * kw * refReNumber;	   //! Distructive term of Omega.

                        cdissk = dissk / (rho * ke);
                        cdissw = two * dissw / (rho * kw);

                        //! To determine if DES used or not.
                        if (isdes == 1)    //! DES-97.
                        {
                            cdes = cblend * cdes_kw + (1.0 - cblend) * cdes_ke;
                            dxyz = GetDistance(grid,i,j,k);
                            desdxyz = cdes * dxyz;
                            tmp = sqrt(ke) / (SST_betaStar * kw);
                            lkw = tmp / (1.0 + a1f2mt2);
                            cdes   = MIN(lkw, desdxyz);
                            cdissk *= tmp / cdes;
                        }

                        prodk = MIN(prodk, SSTProductLimit * dissk);
                        prodw = turb_alphaw * rho / MAX(mut,1.E-20) * prodk;
                        const RDouble zero = 0.0;
                        cprodk = two3rd * (max(divv, zero) + divv * divv * oRefReNumber / kw);

                        //! Produck term based on vorticity.
                        if (SSTProductType != 0)
                        {
                            cprodk = two3rd * max(divv, zero);
                        }
                        cprodw = turb_alphaw * two3rd * max(divv, zero);

                        cdkww = (one - cblend) * two * rho * KW_sigmaW2 * crss * oRefReNumber / kw;
                        ccdkww = ABS(cdkww / rho / kw);

                        srck  = prodk - dissk;
                        srcw  = prodw  - dissw + cdkww;

                        residualTurbulent(i, j, k, IKE) += srck * volume;
                        residualTurbulent(i, j, k, IKW) += srcw * volume;

                        if (neasm > 0)
                        {
                            sourceEnergy(i, j, k) = - srck * volume;
                        }

                        source_dk_neg = ABS(cprodk + cdissk);
                        source_dw_neg = ABS(cprodw + cdissw + ccdkww);

                        spectumTurbulent(i, j, k, IKE) += source_dk_neg * volume;
                        spectumTurbulent(i, j, k, IKW) += source_dw_neg * volume;

                        //! To determine if compressible correction used or not.
                        if (isCC == 1)
                        {
                            if (isdes == 1)
                            {
                                dissk_cc = 0.0;
                            }
                            else
                            {
                                dissk_cc = a1f2mt2 * dissk;
                            }
                            dissw_cc = -a1f2mt2 * SST_betaStar /turb_beta * dissw;

                            pdbar = f2mt2 * (a3 * dissk - a2 * prodk);

                            prodk_cc = pdbar;

                            prodw_cc = -pdbar * rho / mut * refReNumber;

                            residualTurbulent(i, j, k, IKE) += (prodk_cc - dissk_cc) * volume;
                            residualTurbulent(i, j, k, IKW) += (prodw_cc - dissw_cc) * volume;
                        }
                    }
                }
            }
        }
    }
}

void TurbSolverStr::Diagonal(Grid *gridIn)
{
    Param_TurbSolverStruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    if (nTurbulentEquation == 1)
    {
        SpectrumRadiusOfOneEquation(gridIn);
    }
    else if (nTurbulentEquation == 2)
    {
        SpectrumRadiusOfTwoEquation(gridIn);
    }
}

void TurbSolverStr::SpectrumRadiusOfOneEquation(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    RDouble4D &faceVectorX  = *(grid->GetFaceVectorX());
    RDouble4D &faceVectorY  = *(grid->GetFaceVectorY());
    RDouble4D &faceVectorZ  = *(grid->GetFaceVectorZ());
    RDouble4D &area         = *(grid->GetFaceArea());
    RDouble3D &volume       = *(grid->GetCellVolume());

    RDouble4D &spectumTurbulent = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("spec_turb"));
    RDouble4D &qLaminar         = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"        ));
    RDouble3D &viscousLaminar   = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("visl"     ));
    RDouble4D &qTurbulent       = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb"   ));
    RDouble3D &dt               = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("dt"       ));
                               
    RDouble5D &matrixTurbLeft   = *reinterpret_cast<RDouble5D *> (grid->GetDataPtr("mat_turbl"));
    RDouble5D &matrixTurbRight  = *reinterpret_cast<RDouble5D *> (grid->GetDataPtr("mat_turbr"));
    Param_TurbSolverStruct *parameters = GetControlParameters();

    int nTurbulentEquation = parameters->GetNTurbulentEquation();
    RDouble turbCFLScale = parameters->GetTurbCFLScale();

    RDouble oRefReNumber = parameters->GetoRefReNumber();

    RDouble csrv = 1.0;

    //! Initialize some constants.
    RDouble sigma  = 2.0 / 3.0;
    RDouble osigma = one / sigma;

    matrixTurbLeft  = 0;
    matrixTurbRight = 0;

    int nDim = GetDim();

    using namespace IDX;

    int il, jl, kl;
    for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
    {
            int il1, jl1, kl1;
            GetNsurfIndex(iSurface, il1, jl1, kl1);

            int ist, ied, jst, jed, kst, ked;
        grid->GetFaceIterationIndex(ist, ied, jst, jed, kst, ked, iSurface);

            for (int k = kst; k <= ked; ++ k)
            {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    il = i + il1;
                    jl = j + jl1;
                    kl = k + kl1;

                    //RDouble nx = faceVectorX(il, jl, kl, iSurface);
                    //RDouble ny = faceVectorY(il, jl, kl, iSurface);
                    //RDouble nz = faceVectorZ(il, jl, kl, iSurface);

                    //RDouble ul  = qLaminar(i, j, k, IU);
                    //RDouble vl  = qLaminar(i, j, k, IV);
                    //RDouble wl  = qLaminar(i, j, k, IW);
                    //RDouble vnl = nx * ul + ny * vl + nz * wl;
                    //RDouble absVnLeft = ABS(vnl);
			    
                    //RDouble ur  = qLaminar(il, jl, kl, IU);
                    //RDouble vr  = qLaminar(il, jl, kl, IV);
                    //RDouble wr  = qLaminar(il, jl, kl, IW);
                    //RDouble vnr = nx * ur + ny * vr + nz * wr;
                    //RDouble absVnRight = ABS(vnr);
			    
                    //RDouble absVn = half * (absVnLeft + absVnRight);


                    RDouble nx0 = faceVectorX(i, j, k, iSurface);
                    RDouble ny0 = faceVectorY(i, j, k, iSurface);
                    RDouble nz0 = faceVectorZ(i, j, k, iSurface);

                    RDouble nx1 = faceVectorX(il, jl, kl, iSurface);
                    RDouble ny1 = faceVectorY(il, jl, kl, iSurface);
                    RDouble nz1 = faceVectorZ(il, jl, kl, iSurface);

                    RDouble nx = 0.5 * (nx0 + nx1);
                    RDouble ny = 0.5 * (ny0 + ny1);
                    RDouble nz = 0.5 * (nz0 + nz1);

                    RDouble ul  = qLaminar(i, j, k, IU);
                    RDouble vl  = qLaminar(i, j, k, IV);
                    RDouble wl  = qLaminar(i, j, k, IW);
                    RDouble Vn = ul * nx + vl * ny + wl * nz;
                    RDouble absVn = ABS(Vn);

                    //spectumTurbulent(i , j , k , ISA) += half * absVn;
                    //spectumTurbulent(il, jl, kl, ISA) += half * absVn;

                    spectumTurbulent(i , j , k , ISA) += absVn;
                    			    
                    matrixTurbLeft(il, jl, kl, ISA, iSurface)  += half * (- Vn - absVn);
                    matrixTurbRight(i, j, k, ISA, iSurface)    += half * (  Vn - absVn);

                    //matrixTurbLeft(i, j, k, ISA, iSurface)  += half * (- Vn - absVn);
                    //matrixTurbRight(i, j, k, ISA, iSurface) += half * (  Vn - absVn);

                    //! Calculate jacobians of the diffusion operator at the cell face i+1/2.
                    //RDouble rl = qLaminar(i , j , k , IR);
                    //RDouble rr = qLaminar(il, jl, kl, IR);

                    //RDouble orl = 1.0 / (rl + SMALL);
                    //RDouble orr = 1.0 / (rr + SMALL);

                    //RDouble nuLeft    = qTurbulent(i , j , k , ISA);
                    //RDouble nuRight   = qTurbulent(il, jl, kl, ISA);

                    //RDouble nueff_l = ABS(viscousLaminar(i , j , k ) * orl + nuLeft);
                    //RDouble nueff_r = ABS(viscousLaminar(il, jl, kl) * orr + nuRight);
                    //RDouble nueff   = half * (nueff_l + nueff_r);

                    //RDouble cellVolume  = half * (volume(i, j, k) + volume(il, jl, kl));
                    //RDouble oVolume = one / cellVolume;

                    //RDouble ns  = area(il, jl, kl, iSurface);
                    //RDouble ns2 = ns * ns;
                    //ns2 = 2.0 * ns2;

                    //spectumTurbulent(i , j , k , ISA) += oRefReNumber * oVolume * (osigma * nueff * ns2);
                    //spectumTurbulent(il, jl, kl, ISA) += oRefReNumber * oVolume * (osigma * nueff * ns2);


                    RDouble rl = qLaminar(i, j, k, IR);
                    RDouble nul = qTurbulent(i, j, k, ISA);
                    RDouble nueff_l = ABS(viscousLaminar(i , j , k ) / rl + nul);

                    RDouble cellVolume  = volume(i, j, k);
                    RDouble oVolume = one / cellVolume;
                    RDouble ns  = nx * nx +  ny * ny +  nz * nz;
                    RDouble ns2 = 2.0 * ns;

                    spectumTurbulent(i , j , k , ISA) += oRefReNumber * oVolume * (osigma * nueff_l * ns2);

                    matrixTurbLeft(il, jl, kl, ISA, iSurface)  += - half * oRefReNumber * oVolume * (osigma * nueff_l * ns2);
                    matrixTurbRight(i, j, k, ISA, iSurface)    += - half * oRefReNumber * oVolume * (osigma * nueff_l * ns2);

                }
            }
        }
    }

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

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

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble &specTub = spectumTurbulent(i, j, k, ISA);
                specTub += dualTimeSpectrumC1 * volume(i, j, k) + dualTimeSpectrumC2 / (turbCFLScale * dt(i, j, k) + SMALL);
            }
        }
    }
    GhostCell3D(spectumTurbulent, ni, nj, nk, nTurbulentEquation);
}

void TurbSolverStr::SpectrumRadiusOfTwoEquation(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    RDouble4D &faceVectorX  = *(grid->GetFaceVectorX());
    RDouble4D &faceVectorY  = *(grid->GetFaceVectorY());
    RDouble4D &faceVectorZ  = *(grid->GetFaceVectorZ());
    RDouble4D &area   = *(grid->GetFaceArea()  );
    RDouble3D &volume = *(grid->GetCellVolume());

    RDouble4D &spectumTurbulent = *reinterpret_cast< RDouble4D *> (grid->GetDataPtr("spec_turb"));
    RDouble4D &qLaminar         = *reinterpret_cast< RDouble4D *> (grid->GetDataPtr("q"        ));
    RDouble3D &viscousLaminar   = *reinterpret_cast< RDouble3D *> (grid->GetDataPtr("visl"     ));
    RDouble3D &viscousTurbulent = *reinterpret_cast< RDouble3D *> (grid->GetDataPtr("vist"     ));
    RDouble4D &qTurbulent       = *reinterpret_cast< RDouble4D *> (grid->GetDataPtr("q_turb"   ));
    RDouble3D &dt               = *reinterpret_cast< RDouble3D *> (grid->GetDataPtr("dt"       ));

    RDouble5D &matrixTurbLeft  = *reinterpret_cast< RDouble5D *> (grid->GetDataPtr("mat_turbl" ));
    RDouble5D &matrixTurbRight = *reinterpret_cast< RDouble5D *> (grid->GetDataPtr("mat_turbr" ));

    Param_TurbSolverStruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();
    RDouble turbCFLScale = parameters->GetTurbCFLScale();
    string viscousName = parameters->GetViscousName();

    Range M(0, nTurbulentEquation - 1);
    int mst = M.first();
    int med = M.last();

    RDouble oRefReNumber = parameters->GetoRefReNumber();

    RDouble csrv = 1.0;

    matrixTurbLeft  = 0.0;
    matrixTurbRight = 0.0;

    using namespace IDX;

    RDouble dualTimeSpectrumC1 = zero;
    RDouble dualTimeSpectrumC2 = one;

    //! If flow is unsteady, it needs to count the contribution of unsteady.
    int isUnsteady = parameters->GetIsUnsteady();
    if (isUnsteady)
    {
        RDouble dualTimeCoefficient[7];
        const int methodOfDualTime = 3;
        ComputeDualTimeCoefficient(methodOfDualTime, dualTimeCoefficient);
        dualTimeSpectrumC1 = - dualTimeCoefficient[3];
        dualTimeSpectrumC2 =   dualTimeCoefficient[6];
    }

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble volumeCell = volume(i, j, k);
                RDouble dtCell = dt(i, j, k);
                RDouble specTurbDelta = dualTimeSpectrumC1 * volumeCell + dualTimeSpectrumC2 / (turbCFLScale * dtCell + SMALL);
                for (int m = mst; m <= med; ++ m)
                {
                    RDouble &specTur = spectumTurbulent(i, j, k, m);
                    specTur += specTurbDelta;
                }
            }
        }
    }

    int nDim = GetDim();

    vector<RDouble> work(nTurbulentEquation);

    //int il, jl, kl;
    for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
    {
        int il1, jl1, kl1;
        GetNsurfIndex(iSurface, il1, jl1, kl1);

        int ist, ied, jst, jed, kst, ked;
        grid->GetFaceIterationIndex(ist, ied, jst, jed, kst, ked, iSurface);

        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    int il = i + il1;
                    int jl = j + jl1;
                    int kl = k + kl1;

                    RDouble nx0 = faceVectorX(i, j, k, iSurface);
                    RDouble ny0 = faceVectorY(i, j, k, iSurface);
                    RDouble nz0 = faceVectorZ(i, j, k, iSurface);

                    RDouble nx1 = faceVectorX(il, jl, kl, iSurface);
                    RDouble ny1 = faceVectorY(il, jl, kl, iSurface);
                    RDouble nz1 = faceVectorZ(il, jl, kl, iSurface);

                    RDouble nx = 0.5 * (nx0 + nx1);
                    RDouble ny = 0.5 * (ny0 + ny1);
                    RDouble nz = 0.5 * (nz0 + nz1);

                    RDouble ul  = qLaminar(i, j, k, IU);
                    RDouble vl  = qLaminar(i, j, k, IV);
                    RDouble wl  = qLaminar(i, j, k, IW);
                    RDouble Vn = ul * nx + vl * ny + wl * nz;
                    RDouble absVn = ABS(Vn);

                    //RDouble ul  = qLaminar(i, j, k, IU);
                    //RDouble vl  = qLaminar(i, j, k, IV);
                    //RDouble wl  = qLaminar(i, j, k, IW);
                    //RDouble vnl = nx * ul + ny * vl + nz * wl;
                    //RDouble absVnleft = ABS(vnl);

                    //RDouble ur  = qLaminar(il, jl, kl, IU);
                    //RDouble vr  = qLaminar(il, jl, kl, IV);
                    //RDouble wr  = qLaminar(il, jl, kl, IW);
                    //RDouble vnr = nx * ur + ny * vr + nz * wr;
                    //RDouble absVnRight = ABS(vnr);

                    ////RDouble absVn = half * (absVnleft + absVnRight);
                    //RDouble absVn = MAX(absVnleft, absVnRight);

                    for (int m = M.first(); m <= M.last(); ++ m)
                    {
                        //RDouble &specTurbLeft  = spectumTurbulent(i, j, k, m);
                        //RDouble &specTurbRight = spectumTurbulent(il, jl, kl, m);
                        //specTurbLeft += half * absVn;
                        //specTurbRight += half * absVn;

                        RDouble &specTurb  = spectumTurbulent(i, j, k, m);
                        specTurb += absVn; 

                        RDouble &matTurbLeft  = matrixTurbLeft(il, jl, kl, m, iSurface);
                        RDouble &matTurbRight = matrixTurbRight(i, j, k, m, iSurface);
                        matTurbLeft  += half * (- Vn - absVn);
                        matTurbRight += half * (  Vn - absVn);
                    }
                }
            }
        }
    }

    if (viscousName.substr(0,13) == "2eq-kw-menter" || viscousName.substr(0,12) == "easm-kw-2005")
    {
        RDouble3D &blend = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("blend"));

        RDouble KW_sigmaK1 = parameters->GetKW_sigmaK1();
        RDouble KW_sigmaK2 = parameters->GetKW_sigmaK2();
        RDouble KW_sigmaW1 = parameters->GetKW_sigmaW1();
        RDouble KW_sigmaW2 = parameters->GetKW_sigmaW2();

        //int il, jl, kl;
        for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
        {
            int il1, jl1, kl1;
            GetNsurfIndex(iSurface, il1, jl1, kl1);

            int ist, ied, jst, jed, kst, ked;
            grid->GetFaceIterationIndex(ist, ied, jst, jed, kst, ked, iSurface);

			Range M(0, nTurbulentEquation - 1);
            int mst = M.first();
            int med = M.last();

            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    for (int i = ist; i <= ied; ++ i)
                    {
                        int il = i + il1;
                        int jl = j + jl1;
                        int kl = k + kl1;

                        RDouble rhoAverage = qLaminar(i, j, k,IR);
                        RDouble muLaminarAverage = viscousLaminar(i ,j ,k);
                        RDouble muTurbulentAverage = viscousTurbulent(i ,j ,k);

                        RDouble cblendAverage = blend(i, j, k);
                        RDouble tsigKTurbulent = cblendAverage * KW_sigmaK1 + (1.0 - cblendAverage) * KW_sigmaK2;
                        RDouble tsigWTurbulent = cblendAverage * KW_sigmaW1 + (1.0 - cblendAverage) * KW_sigmaW2;

                        RDouble cellVolumeAverage  = volume(i, j, k);
                        RDouble oVolume = one / cellVolumeAverage;

                        RDouble nx0 = faceVectorX(i, j, k, iSurface);
                        RDouble ny0 = faceVectorY(i, j, k, iSurface);
                        RDouble nz0 = faceVectorZ(i, j, k, iSurface);

                        RDouble nx1 = faceVectorX(il, jl, kl, iSurface);
                        RDouble ny1 = faceVectorY(il, jl, kl, iSurface);
                        RDouble nz1 = faceVectorZ(il, jl, kl, iSurface);

                        RDouble nx = 0.5 * (nx0 + nx1);
                        RDouble ny = 0.5 * (ny0 + ny1);
                        RDouble nz = 0.5 * (nz0 + nz1);

                        RDouble ns  = nx * nx + ny * ny + nz * nz;
                        RDouble ns2 = 2.0 * ns;

                        //RDouble rhoAverage = half * (qLaminar(i, j, k,IR) + qLaminar(il, jl, kl, IR));
                        //RDouble muLaminarAverage = half * (viscousLaminar(i ,j ,k) + viscousLaminar(il, jl, kl));
                        //RDouble muTurbulentAverage = half * (viscousTurbulent(i ,j ,k) + viscousTurbulent(il, jl, kl));

                        //RDouble cblendAverage    = half * (blend(i, j, k) + blend(il, jl, kl));
                        //RDouble tsigKTurbulent = cblendAverage * KW_sigmaK1 + (1.0 - cblendAverage) * KW_sigmaK2;
                        //RDouble tsigWTurbulent = cblendAverage * KW_sigmaW1 + (1.0 - cblendAverage) * KW_sigmaW2;

                        //RDouble cellVolumeAverage  = half * (volume(i, j, k) + volume(il, jl, kl));
                        //RDouble oVolume = one / cellVolumeAverage;

                        //RDouble ns  = area(il, jl, kl, iSurface);
                        //RDouble ns2 = ns * ns;
                        //ns2    = 2.0 * ns2;

                        work[0] = (muLaminarAverage + muTurbulentAverage * tsigKTurbulent) / rhoAverage * ns2 * oVolume * oRefReNumber;
                        work[1] = (muLaminarAverage + muTurbulentAverage * tsigWTurbulent) / rhoAverage * ns2 * oVolume * oRefReNumber;

                        for (int m = mst; m <= med; ++ m)
                        {
                            
                            RDouble &specTurb  = spectumTurbulent(i, j, k, m);
                            specTurb += work[m];    

                            RDouble &matTurbLeft  = matrixTurbLeft(il, jl, kl, m, iSurface);
                            RDouble &matTurbRight = matrixTurbRight(i, j, k, m, iSurface);
                            matTurbLeft  += - half * work[m];
                            matTurbRight += - half * work[m];

                            ////RDouble &specTurbLeft  = spectumTurbulent(i, j, k, m);
                            ////RDouble &specTurbRight = spectumTurbulent(il, jl, kl, m);

                            //specTurbLeft += csrv * oRefReNumber * work[m];
                            //specTurbRight += csrv * oRefReNumber * work[m];

                            //RDouble &matTurbLeft  = matrixTurbLeft(il, jl, kl, m, iSurface);
                            //RDouble &matTurbRight = matrixTurbRight(il, jl, kl, m, iSurface);
                            //matTurbLeft  += - half * oRefReNumber * work[m];
                            //matTurbRight += - half * oRefReNumber * work[m];
                        }
                    }
                }
            }
        }
    }
}

void TurbSolverStr::ChangeTurbQ(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);

    Param_TurbSolverStruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();
    if (nTurbulentEquation == 1)
    {
        return;
    }

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 2);

    RDouble4D &qTurbulent = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb"));
    RDouble4D &qLaminar   = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));

    RDouble rho;

    Range M(0, nTurbulentEquation - 1);
    int mst = M.first();
    int med = M.last();

    using namespace IDX;

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                rho = qLaminar(i, j, k, IR) + SMALL;
                for (int m = mst; m <= med; ++ m)
                {
                    qTurbulent(i, j, k, m) *= rho;
                }
            }
        }
    }
}

void TurbSolverStr::RecoverTurbQ(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);

    Param_TurbSolverStruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();
    if (nTurbulentEquation == 1)
    {
        return;
    }

    RDouble4D &qTurbulent = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb"));
    RDouble4D &qLaminar   = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 2);

    Range M(0, nTurbulentEquation - 1);

    int mst = M.first();
    int med = M.last();

    using namespace IDX;

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble rho = qLaminar(i, j, k,IR) + SMALL;
                for (int m = mst; m <= med; ++ m)
                {
                    RDouble &qTurb  = qTurbulent(i, j, k, m);
                    qTurb /= rho;
                }
            }
        }
    }
}

void TurbSolverStr::InviscidFlux(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);

    FieldProxy *qTurbulentLeftPproxy  = CreateFieldProxy(grid);
    FieldProxy *qTurbulentRightPproxy = CreateFieldProxy(grid);
    FieldProxy *fluxProxy = CreateFieldProxy(grid);

    FieldProxy *qLaminarLeftProxy = new FieldProxy();
    FieldProxy *qLaminarRightProxy = new FieldProxy();

    ChangeTurbQ(grid);

    int nDim = GetDim();
    for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
    {

        if (iSurface == 1)
		{
			qLaminarLeftProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("ql1")));
			qLaminarRightProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("qr1")));
		}
		else if (iSurface == 2)
		{
			qLaminarLeftProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("ql2")));
			qLaminarRightProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("qr2")));
		}
		else
		{
			qLaminarLeftProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("ql3")));
			qLaminarRightProxy->SetField_STR(reinterpret_cast <RDouble4D *> (grid->GetDataPtr("qr3")));
		}	

        //! Interpolate den*k den*w.
        GetInviscidFaceValue(grid, qTurbulentLeftPproxy, qTurbulentRightPproxy, iSurface);
        
        //! Interpolate u ,v, w, for turb equation.
        //GetInvFaceQValueforTurb(grid, qLaminarLeftProxy, qLaminarRightProxy, iSurface);

        //! Correct the variables at the face where the BC condition is set to be wall.
        CorrectFaceVar(grid, qTurbulentLeftPproxy, qTurbulentRightPproxy, qLaminarLeftProxy, qLaminarRightProxy, iSurface);

        ComputeInviscidFlux(grid, qTurbulentLeftPproxy, qTurbulentRightPproxy, qLaminarLeftProxy, qLaminarRightProxy, fluxProxy, iSurface);
        
        LoadFlux(grid, fluxProxy, iSurface);
    }

    RecoverTurbQ(grid);

    delete qTurbulentLeftPproxy;
    delete qTurbulentRightPproxy;

    delete qLaminarLeftProxy;
    delete qLaminarRightProxy;

    delete fluxProxy;
}

void TurbSolverStr::GetInviscidFaceValue(Grid *gridIn, FieldProxy *qLeftProxy, FieldProxy *qRightProxy, int iSurface)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    Param_TurbSolverStruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    Range M(0, nTurbulentEquation - 1);

    RDouble4D &qTurbulent = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb" ));

    RDouble4D &qLeft = qLeftProxy->GetField_STR();
    RDouble4D &qRight = qRightProxy->GetField_STR();

    RDouble MUSCLCoefXk = GlobalDataBase::GetDoubleParaFromDB("MUSCLCoefXk");

    int il1, jl1, kl1;
    GetNsurfIndex(iSurface, il1, jl1, kl1);

    I = Range(1 - il1, ni - 1 + il1);
    J = Range(1 - jl1, nj - 1 + jl1);
    K = Range(1 - kl1, nk - 1 + kl1);
    if (nk == 1)
    {
        K.setRange(1,1);
    }

    RDouble c1 = 0.25 * (1.0 - MUSCLCoefXk);
    RDouble c2 = 0.25 * (1.0 + MUSCLCoefXk);

    RDouble4D dql(I, J, K, M, fortranArray);
    RDouble4D dqr(I, J, K, M, fortranArray);
    RDouble4D dqlr(I, J, K, M, fortranArray);

    dql(I, J, K, M) = qTurbulent(I, J, K, M) - qTurbulent(I-il1, J-jl1, K-kl1, M);
    dqr(I, J, K, M) = qTurbulent(I+il1, J+jl1, K+kl1, M) - qTurbulent(I, J, K, M);
    dqlr(I, J, K, M) = Smoothturb(dql(I, J, K, M), dqr(I, J, K, M));

    qLeft (I, J, K, M) = qTurbulent(I, J, K, M) + 0.25 * dqlr(I, J, K, M) * ((1.0 - MUSCLCoefXk * dqlr(I, J, K, M)) * dql(I, J, K, M) + (1.0 + MUSCLCoefXk * dqlr(I, J, K, M)) * dqr(I, J, K, M));
    qRight(I, J, K, M) = qTurbulent(I, J, K, M) - 0.25 * dqlr(I, J, K, M) * ((1.0 - MUSCLCoefXk * dqlr(I, J, K, M)) * dqr(I, J, K, M) + (1.0 + MUSCLCoefXk * dqlr(I, J, K, M)) * dql(I, J, K, M));
}

void TurbSolverStr:: GetInvFaceQValueforTurb(Grid *gridIn, FieldProxy *qLeftProxy, FieldProxy *qRightProxy, int iSurface)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    int neqn = GetNumberOfEquations();

    Range M(1, 3);

    FieldProxy *qProxy = GetFieldProxy(grid, "q");

    RDouble4D &q      = qProxy->GetField_STR();
    RDouble4D &qLeft  = qLeftProxy->GetField_STR();
    RDouble4D &qRight = qRightProxy->GetField_STR();

    int il1, jl1, kl1;
    GetNsurfIndex(iSurface, il1, jl1, kl1);

    I = Range(1 - il1, ni - 1 + il1);
    J = Range(1 - jl1, nj - 1 + jl1);
    K = Range(1 - kl1, nk - 1 + kl1);
    if (nk == 1)
    {
        K.setRange(1, 1);
    }

    RDouble4D dql(I, J, K, M, fortranArray);
    RDouble4D dqr(I, J, K, M, fortranArray);
    RDouble4D dqlr(I, J, K, M, fortranArray);

    RDouble MUSCLCoefXk = GlobalDataBase::GetDoubleParaFromDB("MUSCLCoefXk");

    RDouble c1 = fourth * (1.0 - MUSCLCoefXk);
    RDouble c2 = fourth * (1.0 + MUSCLCoefXk);

    dql(I, J, K, M) = q(I, J, K, M) - q(I-il1, J-jl1, K-kl1, M);
    dqr(I, J, K, M) = q(I+il1, J+jl1, K+kl1, M) - q(I, J, K, M);
    dqlr(I, J, K, M) = Smooth(dql(I, J, K, M), dqr(I, J, K, M));

    qLeft (I, J, K, M) = q(I, J, K, M) + fourth * dqlr(I, J, K, M) * ((1.0 - MUSCLCoefXk * dqlr(I, J, K, M)) * dql(I, J, K, M) + (1.0 + MUSCLCoefXk * dqlr(I, J, K, M)) * dqr(I, J, K, M));
    qRight(I, J, K, M) = q(I, J, K, M) - fourth * dqlr(I, J, K, M) * ((1.0 - MUSCLCoefXk * dqlr(I, J, K, M)) * dqr(I, J, K, M) + (1.0 + MUSCLCoefXk * dqlr(I, J, K, M)) * dql(I, J, K, M));

    delete qProxy;
}

void TurbSolverStr:: CorrectFaceVar(Grid *gridIn, FieldProxy *qLeftProxy, FieldProxy *qRightProxy, FieldProxy *qLaminarLeftProxy, FieldProxy *qLaminarRightProxy, int iSurface)
{
    StructGrid *grid = StructGridCast(gridIn);

    Param_TurbSolverStruct *parameters = GetControlParameters();
    int viscousType = parameters->GetViscousType();

    using namespace IDX;
    CompositeBCRegion *compositBCRegion = grid->GetCompositeBCRegion();
    int nBCRegion = compositBCRegion->GetnBCRegion();

    int nTurbulentEquation = parameters->GetNTurbulentEquation();
    Range M(0, nTurbulentEquation - 1);

    RDouble4D &qLaminarLeft = qLaminarLeftProxy->GetField_STR();   // u,v,w.
    RDouble4D &qLaminarRight = qLaminarRightProxy->GetField_STR();

    RDouble4D &qLeftTurbulent  = qLeftProxy->GetField_STR();   // variables in turb eqautions.
    RDouble4D &qRightTurbulent = qRightProxy->GetField_STR();

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        BCRegion *BCRegion = compositBCRegion->GetBCRegion(iBCRegion);
        int BCType = BCRegion->GetBCType();

        if (!(IsWall(BCType) && viscousType > 0))
        {
            continue;
        }

        int *faceDirectionIndex  = BCRegion->GetFaceDirectionIndex();
        int iSurfOfBC = BCRegion->GetFaceDirection() + 1;
        int leftOrRightIndex       = BCRegion->GetFaceLeftOrRightIndex();

        int ist, ied, jst, jed, kst, ked;
        BCRegion->GetIJKRegion(ist, ied, jst, jed, kst, ked);

        if (iSurfOfBC != iSurface)
        {
            continue;
        }

        int il, jl, kl;
        int id, jd, kd;

        GetBCFaceIDX(faceDirectionIndex, id, jd, kd);

        kst = kst + kd;
        ked = ked + kd;
        jst = jst + jd;
        jed = jed + jd;
        ist = ist + id;
        ied = ied + id;

        int di,dj,dk;

        di = ABS(faceDirectionIndex[0]);
        dj = ABS(faceDirectionIndex[1]);
        dk = ABS(faceDirectionIndex[2]);

        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    il = i - di;
                    jl = j - dj;
                    kl = k - dk; 
                    
                    //qLaminarLeft(il, jl, kl, IU) = zero;
                    //qLaminarLeft(il, jl, kl, IV) = zero;
                    //qLaminarLeft(il, jl, kl, IW) = zero;  //for velocity only 
                    //qLaminarRight(i, j, k, IU) = zero;
                    //qLaminarRight(i, j, k, IV) = zero;
                    //qLaminarRight(i, j, k, IW) = zero; 
                    
                    //! Nest is for turb equations.
                    if (nTurbulentEquation == 1)  //! For SA model.
                    {	
                        qLeftTurbulent(il, jl, kl, ISA) = zero; 
                        qRightTurbulent(i, j, k,ISA) = zero; 
                    }
                    
                    if (nTurbulentEquation == 2)  //! For k-w-SST model.
                    {	
                        qLeftTurbulent(il, jl, kl, IKE) = zero; 
                        qRightTurbulent(i, j, k,IKE) = zero;
                    }
                }
            }
        }
    }
}

void TurbSolverStr::ComputeInviscidFlux(Grid *gridIn, FieldProxy *qLeftProxy, FieldProxy *qRightProxy, FieldProxy *qLaminarLeftProxy, FieldProxy *qLaminarRightProxy, FieldProxy *fluxProxy, int iSurface)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Param_TurbSolverStruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    int str_scheme = GlobalDataBase::GetIntParaFromDB("str_scheme");

    NNDFlux(grid, qLeftProxy, qRightProxy, qLaminarLeftProxy, qLaminarRightProxy, fluxProxy, iSurface);
}

void TurbSolverStr::NNDFlux(Grid *gridIn, FieldProxy *qLeftProxy, FieldProxy *qRightProxy, FieldProxy *qLaminarLeftProxy, FieldProxy *qLaminarRightProxy, FieldProxy *fluxProxy, int iSurface)
{
    StructGrid *grid = StructGridCast(gridIn);

    RDouble4D &qLeft = qLeftProxy->GetField_STR();
    RDouble4D &qRight = qRightProxy->GetField_STR();

    RDouble4D &qLaminarLeft = qLaminarLeftProxy->GetField_STR();
    RDouble4D &qLaminarRight = qLaminarRightProxy->GetField_STR();

    RDouble4D &flux  = fluxProxy->GetField_STR();

    RDouble4D &faceNormalX  = *(grid->GetFaceNormalX());
    RDouble4D &faceNormalY  = *(grid->GetFaceNormalY());
    RDouble4D &faceNormalZ  = *(grid->GetFaceNormalZ());
    RDouble4D &sarea = *(grid->GetFaceArea()   );
    RDouble3D &svol  = *(grid->GetCellVolume() );
    RDouble4D &svgn  = *(grid->GetFaceNormalVelocity());

    Param_TurbSolverStruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    int il1, jl1, kl1;
    GetNsurfIndex(iSurface, il1, jl1, kl1);

    int ist, ied, jst, jed, kst, ked;
    grid->GetFaceIterationIndex(ist, ied, jst, jed, kst, ked, iSurface);

    Range M(0, nTurbulentEquation - 1);
    int mst = M.first();
    int med = M.last();

    using namespace IDX;

    RDouble *priml = new RDouble[nTurbulentEquation];
    RDouble *primr = new RDouble[nTurbulentEquation];

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
           for (int i = ist; i <= ied; ++ i)
            {
                int kl = k + kl1;
                int jl = j + jl1;
                int il = i + il1;

                RDouble xfn  = faceNormalX (il, jl, kl, iSurface);
                RDouble yfn  = faceNormalY (il, jl, kl, iSurface);
                RDouble zfn  = faceNormalZ (il, jl, kl, iSurface);
                RDouble area = sarea(il, jl, kl, iSurface);
                RDouble vgn  = svgn (il, jl, kl, iSurface);

                for (int m = mst; m <= med; ++ m)
                {
                    priml[m] = qLeft(i, j ,k , m);
                    primr[m] = qRight(il, jl, kl, m);
                }

                RDouble ul = qLaminarLeft(i, j, k, IU);
                RDouble vl = qLaminarLeft(i, j, k, IV);
                RDouble wl = qLaminarLeft(i, j, k, IW);

                RDouble ur = qLaminarRight(il, jl, kl, IU);
                RDouble vr = qLaminarRight(il, jl, kl, IV);
                RDouble wr = qLaminarRight(il, jl, kl, IW);

                RDouble vnl  = xfn * ul + yfn * vl + zfn * wl - vgn;
                RDouble vnr  = xfn * ur + yfn * vr + zfn * wr - vgn;

                RDouble vnla = half * (vnl + ABS(vnl));
                RDouble vnra = half * (vnr - ABS(vnr));

                for (int m = mst; m <= med; ++ m)
                {
                    flux(il, jl, kl, m) = (vnla * priml[m] + vnra * primr[m]) * area;
                }
            }
        }
    }
    delete []priml;
    delete []primr;
}

void TurbSolverStr::LoadFlux(Grid *gridIn,  FieldProxy *fluxProxy, int iSurface)
{
    StructGrid *grid = StructGridCast(gridIn);

    RDouble4D &flux = fluxProxy->GetField_STR();

    Param_TurbSolverStruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    RDouble4D &residualTurbulent = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_turb"));

    int il1, jl1, kl1;
    GetNsurfIndex(iSurface, il1, jl1, kl1);

    int il, jl, kl;

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    Range M(0, nTurbulentEquation - 1);
    int mst = M.first();
    int med = M.last();

    for (int m = mst; m <= med; ++ m)
    {
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    il = i + il1;
                    jl = j + jl1;
                    kl = k + kl1;
                    RDouble &resTurb = residualTurbulent(i, j, k,m);
                    RDouble &fluxLeft = flux(il, jl, kl,m);
                    RDouble &fluxRight = flux(i, j, k,m);
                    resTurb -= (fluxLeft - fluxRight);
                }
            }
        }
    }
}

void TurbSolverStr::ComputeGradientCellCenter(Grid *gridIn)
{
    StructGrid *strgrid = StructGridCast(gridIn);
    int ni = strgrid->GetNI();
    int nj = strgrid->GetNJ();
    int nk = strgrid->GetNK();

    RDouble4D &faceVectorX = *(strgrid->GetFaceVectorX());
    RDouble4D &faceVectorY = *(strgrid->GetFaceVectorY());
    RDouble4D &faceVectorZ = *(strgrid->GetFaceVectorZ());
    RDouble3D &volume      = *(strgrid->GetCellVolume());

    int ndim = GetDim();

    Param_TurbSolverStruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    RDouble4D &qTurbulent    = *reinterpret_cast<RDouble4D *> (strgrid->GetDataPtr("q_turb"));
    RDouble4D &gradTurbulentCellCenterX = *reinterpret_cast<RDouble4D *> (strgrid->GetDataPtr("gradTurbulentCellCenterX"));
    RDouble4D &gradTurbulentCellCenterY = *reinterpret_cast<RDouble4D *> (strgrid->GetDataPtr("gradTurbulentCellCenterY"));
    RDouble4D &gradTurbulentCellCenterZ = *reinterpret_cast<RDouble4D *> (strgrid->GetDataPtr("gradTurbulentCellCenterZ"));

    gradTurbulentCellCenterX = 0.0;
    gradTurbulentCellCenterY = 0.0;
    gradTurbulentCellCenterZ = 0.0;

    for (int iSurface = 1; iSurface <= ndim; ++ iSurface)
    {
        int il1, jl1, kl1;
        GetNsurfIndex(iSurface, il1, jl1, kl1);

        int ist = 1;
        int ied = ni - 1 + il1;
        int jst = 1;
        int jed = nj - 1 + jl1;
        int kst = 1;
        int ked = nk - 1 + kl1;

        if (ndim == TWO_D)
        {
            ked = 1;
        }
        int il, jl, kl;

        //! The gradient of u,v,w.
        for (int m = 0; m < nTurbulentEquation; ++ m)
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    for (int i = ist; i <= ied; ++ i)
                    {
                        il = i - il1;
                        jl = j - jl1;
                        kl = k - kl1;

                        RDouble phis = qTurbulent(i, j, k, m) + qTurbulent(il, jl, kl, m);

                        RDouble ddx = phis * faceVectorX(i, j, k, iSurface);
                        RDouble ddy = phis * faceVectorY(i, j, k, iSurface);
                        RDouble ddz = phis * faceVectorZ(i, j, k, iSurface);

                        RDouble &gradTurbccXLeft = gradTurbulentCellCenterX(i, j, k, m);
                        RDouble &gradTurbccYLeft = gradTurbulentCellCenterY(i, j, k, m);
                        RDouble &gradTurbccZLeft = gradTurbulentCellCenterZ(i, j, k, m);
                        gradTurbccXLeft -= ddx;
                        gradTurbccYLeft -= ddy;
                        gradTurbccZLeft -= ddz;

                        RDouble &gradTurbccXRight = gradTurbulentCellCenterX(il, jl, kl, m);
                        RDouble &gradTurbccYRight = gradTurbulentCellCenterY(il, jl, kl, m);
                        RDouble &gradTurbccZRight = gradTurbulentCellCenterZ(il, jl, kl, m);
                        gradTurbccXRight += ddx;
                        gradTurbccYRight += ddy;
                        gradTurbccZRight += ddz;
                    }
                }
            }
        }
    }

    int ist = 1;
    int ied = ni-1;
    int jst = 1;
    int jed = nj-1;
    int kst = 1;
    int ked = nk-1;

    if (ndim == TWO_D) ked = 1;

    for (int m = 0; m < nTurbulentEquation; ++ m)
    {
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    RDouble oov = half / volume(i, j, k);
                    RDouble &gradTurbccX = gradTurbulentCellCenterX(i, j, k, m);
                    RDouble &gradTurbccY = gradTurbulentCellCenterY(i, j, k, m);
                    RDouble &gradTurbccZ = gradTurbulentCellCenterZ(i, j, k, m);
                    gradTurbccX *= oov;
                    gradTurbccY *= oov;
                    gradTurbccZ *= oov;
                }
            }
        }
    }
}

void TurbSolverStr::CompressSpecificArrayToInterface(DataContainer *&dataContainer, const string &fieldName, Grid *gridIn, const int &neighborZoneIndex, const int &nEquation)
{
    StructGrid *grid = StructGridCast(gridIn);
    int level = grid->GetLevel();
    InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
    if (!interfaceInformation)
    {
        return;
    }

    StructGrid *finestGrid = StructGridCast(grid->GetFinestGrid());
    InterfaceInfo *finestInterfaceInfo = finestGrid->GetInterfaceInfo();

    int iNeighborZone                            = finestInterfaceInfo->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = finestInterfaceInfo->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForSend          = finestInterfaceInfo->GetFaceIndexForSend(iNeighborZone);

    dataContainer->moveto_begin();

    if (0 == nEquation)
    {
        RDouble3D &fieldSend = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr(fieldName));

        for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; -- iGhostLayer)
        {
            for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
            {
                int iFace = interfaceIndexContainerForSend[iLocalFace];
                int is, js, ks;
                finestGrid->GetSourceIndexIJK(iFace, iGhostLayer + 1, is, js, ks);
                finestGrid->RemapMultigridIJK(level, is, js, ks);
                PHWrite(dataContainer, fieldSend(is, js, ks));
            }
        }
    }
    else
    {
        RDouble4D &fieldSend = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr(fieldName));

        for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; -- iGhostLayer)
        {
            for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
            {
                int iFace = interfaceIndexContainerForSend[iLocalFace];
                int is, js, ks;
                finestGrid->GetSourceIndexIJK(iFace, iGhostLayer + 1, is, js, ks);
                finestGrid->RemapMultigridIJK(level, is, js, ks);
                for (int m = 0; m < nEquation; ++ m)
                {
                    PHWrite(dataContainer, fieldSend(is, js, ks, m));
                }
            }
        }
    }
}

void TurbSolverStr::DecompressArrayFromInterface(DataContainer *&dataContainer, const string &fieldName, Grid *gridIn, const int &neighborZoneIndex, const int &nEquation)
{
    StructGrid *grid = StructGridCast(gridIn);
    int level = grid->GetLevel();
    InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
    if (!interfaceInformation)
    {
        return;
    }

    StructGrid *finestGrid = StructGridCast(grid->GetFinestGrid());
    InterfaceInfo *finestInterfaceInfo = finestGrid->GetInterfaceInfo();

    int iNeighborZone                            = finestInterfaceInfo->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = finestInterfaceInfo->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForReceive       = finestInterfaceInfo->GetFaceIndexForRecv(iNeighborZone);

    dataContainer->moveto_begin();

    if (0 == nEquation)
    {
        RDouble3D &fieldRecv = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr(fieldName));

        for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; -- iGhostLayer)
        {
            for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
            {
                int iFace = interfaceIndexContainerForReceive[iLocalFace];
                int it, jt, kt;
                finestGrid->GetTargetIndexIJK(iFace, iGhostLayer + 1, it, jt, kt);
                finestGrid->RemapMultigridIJK(level, it, jt, kt);
                PHRead(dataContainer, fieldRecv(it, jt, kt));
            }
        }
    }
    else
    {
        RDouble4D &fieldRecv = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr(fieldName));

        for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; -- iGhostLayer)
        {
            for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
            {
                int iFace = interfaceIndexContainerForReceive[iLocalFace];
                int it, jt, kt;
                finestGrid->GetTargetIndexIJK(iFace, iGhostLayer + 1, it, jt, kt);
                finestGrid->RemapMultigridIJK(level, it, jt, kt);
                for (int m = 0; m < nEquation; ++ m)
                {
                    PHRead(dataContainer, fieldRecv(it, jt, kt, m));
                }
            }
        }
    }
}

void TurbSolverStr::RotateVectorFromInterface(Grid *gridIn, const int &neighborZoneIndex, const int &nEquation)
{
    StructGrid *grid = StructGridCast(gridIn);
    int level = grid->GetLevel();
    InterfaceInfo *interfaceInformation = grid->GetInterfaceInfo();
    if (!interfaceInformation)
    {
        return;
    }

    StructGrid *finestGrid = StructGridCast(grid->GetFinestGrid());
    InterfaceInfo *finestInterfaceInfo = finestGrid->GetInterfaceInfo();
    CompositeBCRegion * composite_bcregion = finestGrid->GetCompositeBCRegion();
    RDouble4D &rotTurbgradValueX = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("rotTurbgradValueX"));
    RDouble4D &rotTurbgradValueY = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("rotTurbgradValueY"));
    RDouble4D &rotTurbgradValueZ = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("rotTurbgradValueZ"));

    rotTurbgradValueX = 0.0;
    rotTurbgradValueY = 0.0;
    rotTurbgradValueZ = 0.0;

    int iNeighborZone                            = finestInterfaceInfo->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = finestInterfaceInfo->GetNIFaceOfNeighbor(iNeighborZone);
    int *interfaceIndexContainerForReceive       = finestInterfaceInfo->GetFaceIndexForRecv(iNeighborZone);

    double rotationAngle        = PHSPACE::GlobalDataBase::GetDoubleParaFromDB("rotationAngle");
    rotationAngle = rotationAngle*PI/180;
    if(nEquation > 0)
    {
        RDouble4D &fieldRecvTurbX = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradTurbulentCellCenterX"));
        RDouble4D &fieldRecvTurbY = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradTurbulentCellCenterY"));
        RDouble4D &fieldRecvTurbZ = *reinterpret_cast <RDouble4D *> (grid->GetDataPtr("gradTurbulentCellCenterZ"));

        for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; -- iGhostLayer)
        {
            for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
            {
                int iFace = interfaceIndexContainerForReceive[iLocalFace];
                int it, jt, kt;
                finestGrid->GetTargetIndexIJK(iFace, iGhostLayer + 1, it, jt, kt);
                finestGrid->RemapMultigridIJK(level, it, jt, kt);
                int *ibcregions = composite_bcregion->GetIFaceInfo();
                int ibcregion = ibcregions[iFace];
                BCRegion * bcregion = composite_bcregion->GetBCRegion(ibcregion);
                string bcName = bcregion->GetBCName();
                if ( bcName == "Periodic_up" )
                {
                    for (int m = 0; m < nEquation; ++ m)
                    {
                        rotTurbgradValueY(it, jt, kt, m) = fieldRecvTurbY(it, jt, kt, m) * cos(2*PI-rotationAngle) - fieldRecvTurbZ(it, jt, kt, m) * sin(2*PI-rotationAngle);
                        rotTurbgradValueZ(it, jt, kt, m) = fieldRecvTurbY(it, jt, kt, m) * sin(2*PI-rotationAngle) + fieldRecvTurbZ(it, jt, kt, m) * cos(2*PI-rotationAngle);
                    }							
                }
                else if( bcName == "Periodic_down" )
                {
                    for (int m = 0; m < nEquation; ++ m)
                    {
                        rotTurbgradValueY(it, jt, kt, m) = fieldRecvTurbY(it, jt, kt, m) * cos(rotationAngle) - fieldRecvTurbZ(it, jt, kt, m) * sin(rotationAngle);
                        rotTurbgradValueZ(it, jt, kt, m) = fieldRecvTurbY(it, jt, kt, m) * sin(rotationAngle) + fieldRecvTurbZ(it, jt, kt, m) * cos(rotationAngle);
                    }
                }

                for (int m = 0; m < nEquation; ++ m)
                {
                    fieldRecvTurbY(it, jt, kt, m) = rotTurbgradValueY(it, jt, kt, m);
                    fieldRecvTurbZ(it, jt, kt, m) = rotTurbgradValueZ(it, jt, kt, m);
                }
            }
        }
    }
}

void TurbSolverStr::ReconGradVolWeightedwithCorrection(Grid *gridIn, int iSurface)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    int nEquation = GetNumberOfEquations();

    RDouble4D &qTurbulent = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb"));

    RDouble4D &gradTurbulentX = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradTurbulentX"));
    RDouble4D &gradTurbulentY = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradTurbulentY"));
    RDouble4D &gradTurbulentZ = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradTurbulentZ"));

    gradTurbulentX = 0.0;
    gradTurbulentY = 0.0;
    gradTurbulentZ = 0.0;

    RDouble4D &gradTurbulentCellCenterX = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradTurbulentCellCenterX"));
    RDouble4D &gradTurbulentCellCenterY = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradTurbulentCellCenterY"));
    RDouble4D &gradTurbulentCellCenterZ = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradTurbulentCellCenterZ"));

    RDouble3D &xcc = *(grid->GetCellCenterX());
    RDouble3D &ycc = *(grid->GetCellCenterY());
    RDouble3D &zcc = *(grid->GetCellCenterZ());

    int il1, jl1, kl1;
    GetNsurfIndex(iSurface, il1, jl1, kl1);

    Range I(1, ni - 1 + il1);
    Range J(1, nj - 1 + jl1);
    Range K(1, nk - 1 + kl1);
    if (nk == 1)
    {
        K.setRange(1, 1);
    }

    int ist, ied, jst, jed, kst, ked;
    GetIJKRegion(I, J, K, ist, ied, jst, jed, kst, ked);

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble dx = xcc(i, j, k) - xcc(i-il1,j-jl1,k-kl1);
                RDouble dy = ycc(i, j, k) - ycc(i-il1,j-jl1,k-kl1);
                RDouble dz = zcc(i, j, k) - zcc(i-il1,j-jl1,k-kl1);
                RDouble dr = sqrt(dx * dx + dy * dy + dz * dz);

                dx = dx / dr;
                dy = dy / dr;
                dz = dz / dr;

                for (int m = 0; m < nEquation; ++ m)
                {
                    RDouble a = half * (gradTurbulentCellCenterX(i, j, k, m) + gradTurbulentCellCenterX(i - il1, j - jl1, k - kl1, m));
                    RDouble b = half * (gradTurbulentCellCenterY(i, j, k, m) + gradTurbulentCellCenterY(i - il1, j - jl1, k - kl1, m));
                    RDouble c = half * (gradTurbulentCellCenterZ(i, j, k, m) + gradTurbulentCellCenterZ(i - il1, j - jl1, k - kl1, m));
                    RDouble dqdr = (qTurbulent(i, j, k, m) - qTurbulent(i - il1, j - jl1, k - kl1, m)) / dr;
                    RDouble d = a * dx + b * dy + c * dz - dqdr;

                    RDouble &gradTurbX = gradTurbulentX(i, j, k, m);
                    RDouble &gradTurbY = gradTurbulentY(i, j, k, m);
                    RDouble &gradTurbZ = gradTurbulentZ(i, j, k, m);
                    gradTurbX = a - d * dx;
                    gradTurbY = b - d * dy;
                    gradTurbZ = c - d * dz;
                }
            }
        }
    }

    GetGradientAtFace_CorrectionAtPhysicalBoundary(grid, iSurface);
}

void TurbSolverStr:: GetGradientAtFace_CorrectionAtPhysicalBoundary(Grid *gridIn, int iSurface)
{
    StructGrid *grid = StructGridCast(gridIn);

    RDouble4D &xfn  = *(grid->GetFaceNormalX());
    RDouble4D &yfn  = *(grid->GetFaceNormalY());
    RDouble4D &zfn  = *(grid->GetFaceNormalZ());
    RDouble3D &vol  = *(grid->GetCellVolume());

    int nEquation = GetNumberOfEquations();

    RDouble4D &gradTurbulentFaceX = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradTurbulentX"));
    RDouble4D &gradTurbulentFaceY = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradTurbulentY"));
    RDouble4D &gradTurbulentFaceZ = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradTurbulentZ"));

    RDouble4D &gradTurbulentCellCenterX = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradTurbulentCellCenterX"));
    RDouble4D &gradTurbulentCellCenterY = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradTurbulentCellCenterY"));
    RDouble4D &gradTurbulentCellCenterZ = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradTurbulentCellCenterZ"));

    CompositeBCRegion *compositeBCRegion = grid->GetCompositeBCRegion();
    int nBCRegion = compositeBCRegion->GetnBCRegion();

    using namespace PHENGLEI;

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        BCRegion *bcregion = compositeBCRegion->GetBCRegion(iBCRegion);
        int BCType = bcregion->GetBCType();
        int nsurf_bc = bcregion->GetFaceDirection() + 1;

        int ist, ied, jst, jed, kst, ked;
        bcregion->GetIJKRegion(ist, ied, jst, jed, kst, ked);

        if (BCType < 0 || nsurf_bc != iSurface) // inner P TO P BC in the current direction
        {
            continue;
        }
        else if (BCType == SYMMETRY) // symmetry BC
        {
            RDouble xfnSign; 
            RDouble yfnSign; 
            RDouble zfnSign; 

            xfnSign = ABS(xfn(ist,jst,kst,iSurface)) - 0.1;
            yfnSign = ABS(yfn(ist,jst,kst,iSurface)) - 0.1;
            zfnSign = ABS(zfn(ist,jst,kst,iSurface)) - 0.1;

            if (xfnSign > 0.0) // for Y-Z symmetry plane
            {
                for (int k = kst; k <= ked; ++ k)
                {
                    for (int j = jst; j <= jed; ++ j)
                    {
                        for (int i = ist; i <= ied; ++ i)
                        {
							int ibc1, jbc1, kbc1;
							bcregion->GetBoundaryFaceIndex(i, j, k, ibc1, jbc1, kbc1);

                            for (int m = 0; m < nEquation; ++ m)
                            {
                                gradTurbulentFaceX(ibc1, jbc1, kbc1, m) = 0.0                                 ;
                                gradTurbulentFaceY(ibc1, jbc1, kbc1, m) = gradTurbulentCellCenterY(i, j, k, m);
                                gradTurbulentFaceZ(ibc1, jbc1, kbc1, m) = gradTurbulentCellCenterZ(i, j, k, m);
                            }
                        }
                    }
                }
            }

            if (yfnSign > 0.0) // for X-Z symmetry plane
            {
                for (int k = kst; k <= ked; ++ k)
                {
                    for (int j = jst; j <= jed; ++ j)
                    {
                        for (int i = ist; i <= ied; ++ i)
                        {
							int ibc1, jbc1, kbc1;
							bcregion->GetBoundaryFaceIndex(i, j, k, ibc1, jbc1, kbc1);

                            for (int m = 0; m < nEquation; ++ m)
                            {
                                gradTurbulentFaceX(ibc1, jbc1, kbc1, m) = gradTurbulentCellCenterX(i, j, k, m);
                                gradTurbulentFaceY(ibc1, jbc1, kbc1, m) = 0.0                                 ;
                                gradTurbulentFaceZ(ibc1, jbc1, kbc1, m) = gradTurbulentCellCenterZ(i, j, k, m);
                            }
                        }
                    }
                }
            }

            if (zfnSign > 0.0) // for X-Y symmetry plane
            {
                for (int k = kst; k <= ked; ++ k)
                {
                    for (int j = jst; j <= jed; ++ j)
                    {
                        for (int i = ist; i <= ied; ++ i)
                        {
							int ibc1, jbc1, kbc1;
							bcregion->GetBoundaryFaceIndex(i, j, k, ibc1, jbc1, kbc1);

                            for (int m = 0; m < nEquation; ++ m)
                            {
                                gradTurbulentFaceX(ibc1, jbc1, kbc1, m) = gradTurbulentCellCenterX(i, j, k, m);
                                gradTurbulentFaceY(ibc1, jbc1, kbc1, m) = gradTurbulentCellCenterY(i, j, k, m);
                                gradTurbulentFaceZ(ibc1, jbc1, kbc1, m) = 0.0                                 ;
                            }
                        }
                    }
                }
            }
        }
        else
        {
            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    for (int i = ist; i <= ied; ++ i)
                    {       
						int ibc1, jbc1, kbc1;
						bcregion->GetBoundaryFaceIndex(i, j, k, ibc1, jbc1, kbc1);
                    
                        for (int m = 0; m < nEquation; ++ m)
                        {
                            gradTurbulentFaceX(ibc1, jbc1, kbc1, m) = gradTurbulentCellCenterX(i, j, k, m);
                            gradTurbulentFaceY(ibc1, jbc1, kbc1, m) = gradTurbulentCellCenterY(i, j, k, m);
                            gradTurbulentFaceZ(ibc1, jbc1, kbc1, m) = gradTurbulentCellCenterZ(i, j, k, m);
                        }                        
                    }
                }
            }
        }
    }
}

void TurbSolverStr::ViscousFlux(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);

    FieldProxy *fluxProxy = CreateFieldProxy(grid);

    int nDim = GetDim();
    for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
    {
        ReconGradVolWeightedwithCorrection(gridIn, iSurface);

        CompVisFlux(grid, fluxProxy, iSurface);

        LoadFlux(grid, fluxProxy, iSurface);
    }

    delete fluxProxy;
}

void TurbSolverStr::CompVisFlux(Grid *gridIn, FieldProxy *fluxProxy, int iSurface)
{
    Param_TurbSolverStruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();
    if (nTurbulentEquation == 1)
    {
        //! One equation turbulence model.
        CompVisFluxOneEquation(gridIn, fluxProxy, iSurface);
    }
    else if (nTurbulentEquation == 2)
    {
        //! Two equation turbulence model.
        CompVisFluxTwoEquation(gridIn, fluxProxy, iSurface);
    }
}

void TurbSolverStr::CompVisFluxOneEquation(Grid *gridIn, FieldProxy *fluxProxy, int iSurface)
{
    StructGrid *grid = StructGridCast(gridIn);

    RDouble4D &gradientTurbulentX = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradTurbulentX"));
    RDouble4D &gradientTurbulentY = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradTurbulentY"));
    RDouble4D &gradientTurbulentZ = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradTurbulentZ"));

    RDouble4D &flux = fluxProxy->GetField_STR();

    RDouble4D &faceNormalX  = *(grid->GetFaceNormalX());
    RDouble4D &faceNormalY  = *(grid->GetFaceNormalY());
    RDouble4D &faceNormalZ  = *(grid->GetFaceNormalZ());
    RDouble4D &area   = *(grid->GetFaceArea()   );

    RDouble4D &qLaminar = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q" ));
    RDouble4D &qTurbulent = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb"));
    RDouble3D &viscousLaminar = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("visl"));

    Param_TurbSolverStruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    RDouble oRefReNumber = parameters->GetoRefReNumber();

    RDouble KW_sigma = parameters->GetKW_sigma();

    RDouble oSigmaTurbulent = 1.0 / KW_sigma;

    RDouble *fvis = new RDouble[nTurbulentEquation];

    int il1, jl1, kl1;
    GetNsurfIndex(iSurface, il1, jl1, kl1);

    int ist, ied, jst, jed, kst, ked;
    grid->GetFaceIterationIndex(ist, ied, jst, jed, kst, ked, iSurface);

    using namespace IDX;

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                int il = i + il1;
                int jl = j + jl1;
                int kl = k + kl1;

                RDouble nxs = faceNormalX(il, jl, kl, iSurface);
                RDouble nys = faceNormalY(il, jl, kl, iSurface);
                RDouble nzs = faceNormalZ(il, jl, kl, iSurface);
                RDouble nss =         area(il, jl, kl, iSurface);

                RDouble dkedx  = gradientTurbulentX(il, jl, kl, ISA);
                RDouble dkedy  = gradientTurbulentY(il, jl, kl, ISA);
                RDouble dkedz  = gradientTurbulentZ(il, jl, kl, ISA);

                RDouble dkedn = nxs * dkedx + nys * dkedy + nzs * dkedz;

                RDouble orl = 1.0 / (qLaminar(i , j , k , IR) + SMALL);
                RDouble orr = 1.0 / (qLaminar(il, jl, kl, IR) + SMALL);

                //! flux.
                RDouble keViscous   = half * (viscousLaminar(i, j, k) * orl + viscousLaminar(il, jl, kl) * orr + qTurbulent(i, j, k, ISA) + qTurbulent(il, jl, kl, ISA));
                fvis[0]  = - oRefReNumber * (oSigmaTurbulent * keViscous * dkedn * nss);

                for (int m = 0; m < nTurbulentEquation; ++ m)
                {
                    RDouble &fluxface = flux(il, jl, kl, m);
                    fluxface = fvis[m];
                }
            }
        }
    }
    delete []fvis;
}

void TurbSolverStr::CompVisFluxTwoEquation(Grid *gridIn, FieldProxy *fluxProxy, int iSurface)
{
    StructGrid *grid = StructGridCast(gridIn);

    RDouble4D &gradientTurbulentX = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradTurbulentX"));
    RDouble4D &gradientTurbulentY = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradTurbulentY"));
    RDouble4D &gradientTurbulentZ = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradTurbulentZ"));

    RDouble4D &flux   = fluxProxy->GetField_STR();

    RDouble4D &faceNormalX = *(grid->GetFaceNormalX());
    RDouble4D &faceNormalY = *(grid->GetFaceNormalY());
    RDouble4D &faceNormalZ = *(grid->GetFaceNormalZ());
    RDouble4D &area        = *(grid->GetFaceArea()   );
    RDouble3D &vol         = *(grid->GetCellVolume() );

    RDouble4D &qLaminar         = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q" ));
    RDouble4D &qTurbulent       = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb"));
    RDouble3D &viscousLaminar   = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("visl"));
    RDouble3D &viscousTurbulent = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("vist"));

    Param_TurbSolverStruct *parameters = GetControlParameters();
    string viscousName = parameters->GetViscousName();

    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    RDouble oRefReNumber = parameters->GetoRefReNumber();

    RDouble *fvis = new RDouble[nTurbulentEquation];

    int il1, jl1, kl1;
    GetNsurfIndex(iSurface, il1, jl1, kl1);

    int ist, ied, jst, jed, kst, ked;
    grid->GetFaceIterationIndex(ist, ied, jst, jed, kst, ked, iSurface);

    using namespace IDX;

    if (viscousName.substr(0,13) == "2eq-kw-menter" || viscousName.substr(0,12) == "easm-kw-2005")
    {
        RDouble3D &blend = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("blend"));

        RDouble KW_sigmaK1 = parameters->GetKW_sigmaK1();
        RDouble KW_sigmaW1 = parameters->GetKW_sigmaW1();
        RDouble KW_sigmaW2 = parameters->GetKW_sigmaW2();
        RDouble KW_sigmaK2 = parameters->GetKW_sigmaK2();

        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    int il = i + il1;
                    int jl = j + jl1;
                    int kl = k + kl1;

                    RDouble nxs = faceNormalX(il, jl, kl, iSurface);
                    RDouble nys = faceNormalY(il, jl, kl, iSurface);
                    RDouble nzs = faceNormalZ(il, jl, kl, iSurface);
                    RDouble nss =        area(il, jl, kl, iSurface);

                    RDouble dkedx = gradientTurbulentX(il, jl, kl, IKE);
                    RDouble dkedy = gradientTurbulentY(il, jl, kl, IKE);
                    RDouble dkedz = gradientTurbulentZ(il, jl, kl, IKE);

                    RDouble dkwdx = gradientTurbulentX(il, jl, kl, IKW);
                    RDouble dkwdy = gradientTurbulentY(il, jl, kl, IKW);
                    RDouble dkwdz = gradientTurbulentZ(il, jl, kl, IKW);

                    RDouble cblend = half * (blend(i, j, k) + blend(il, jl, kl));
                    RDouble sigkTurbulent = cblend * KW_sigmaK1 + (1.0 - cblend) * KW_sigmaK2;
                    RDouble sigwTurbulent = cblend * KW_sigmaW1 + (1.0 - cblend) * KW_sigmaW2;
                    RDouble muLaminar = half * (viscousLaminar(i, j, k) + viscousLaminar(il, jl, kl));
                    RDouble muTurbulent = half * (viscousTurbulent(i, j, k) + viscousTurbulent(il, jl, kl));

                    RDouble mlt1 = muLaminar + muTurbulent * sigkTurbulent;
                    RDouble mlt2 = muLaminar + muTurbulent * sigwTurbulent;
                    
                    RDouble dkedn = nxs * dkedx + nys * dkedy + nzs * dkedz;
                    RDouble dkwdn = nxs * dkwdx + nys * dkwdy + nzs * dkwdz;

                    fvis[0] = - oRefReNumber * mlt1 * dkedn * nss;
                    fvis[1] = - oRefReNumber * mlt2 * dkwdn * nss;

                    for (int m = 0; m < nTurbulentEquation; ++ m)
                    {
                        RDouble &fluxface = flux(il, jl, kl, m);
                        fluxface = fvis[m];
                    }
                }
            }
        }
    }
    else
    {
        RDouble KW_sigmaK = parameters->GetKW_sigmaK();
        RDouble KW_sigmaW = parameters->GetKW_sigmaW();

        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    int il = i + il1;
                    int jl = j + jl1;
                    int kl = k + kl1;

                    RDouble nxs = faceNormalX(il, jl, kl, iSurface);
                    RDouble nys = faceNormalY(il, jl, kl, iSurface);
                    RDouble nzs = faceNormalZ(il, jl, kl, iSurface);
                    RDouble nss =        area(il, jl, kl, iSurface);

                    RDouble dkedx = gradientTurbulentX(il, jl, kl, IKE);
                    RDouble dkedy = gradientTurbulentY(il, jl, kl, IKE);
                    RDouble dkedz = gradientTurbulentZ(il, jl, kl, IKE);

                    RDouble dkwdx = gradientTurbulentX(il, jl, kl, IKW);
                    RDouble dkwdy = gradientTurbulentY(il, jl, kl, IKW);
                    RDouble dkwdz = gradientTurbulentZ(il, jl, kl, IKW);
           
                    RDouble muLaminar = half * (viscousLaminar(i, j, k) + viscousLaminar(il, jl, kl));
                    RDouble muTurbulent = half * (viscousTurbulent(i, j, k) + viscousTurbulent(il, jl, kl));

                    RDouble mlt1 = muLaminar + muTurbulent * KW_sigmaK;
                    RDouble mlt2 = muLaminar + muTurbulent * KW_sigmaW;
                    
                    RDouble dkedn = nxs * dkedx + nys * dkedy + nzs * dkedz;
                    RDouble dkwdn = nxs * dkwdx + nys * dkwdy + nzs * dkwdz;

                    fvis[IKE] = - oRefReNumber * mlt1 * dkedn * nss;
                    fvis[IKW] = - oRefReNumber * mlt2 * dkwdn * nss;

                    for (int m = 0; m < nTurbulentEquation; ++ m)
                    {
                        RDouble &fluxface = flux(il, jl, kl, m);
                        fluxface = fvis[m];
                    }
                }
            }
        }
    }

    delete []fvis;
}

void TurbSolverStr::ZeroResidualOfSpecialCells(Grid * gridIn)
{
    Param_TurbSolverStruct *parameters = GetControlParameters();
    int isOverLapping = parameters->GetIsOverLapping();
    if (isOverLapping == 0)
    {
        return;
    }

    RDouble4D &residualTurbulent = *reinterpret_cast<RDouble4D *> (gridIn->GetDataPtr("res_turb"));
    StructGrid *grid = PHSPACE::StructGridCast(gridIn);
    int numberOfEquations = GetNumberOfEquations();

    int *cellTypeContainer = grid->GetCellTypeContainer();

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    int cellLabel = 0;
    for (int k = kst; k <= ked; ++k)
    {
        for (int j = jst; j <= jed; ++j)
        {
            for (int i = ist; i <= ied; ++i)
            {
                if (cellTypeContainer[cellLabel] != GENERIC_COLOR)
                {
                    for (int m = 0; m < numberOfEquations; ++ m)
                    {
                        residualTurbulent(i, j, k, m) = 0.0;
                    }
                }
                cellLabel += 1;
            }
        }
    }
}

void TurbSolverStr::SolveLUSGSForward(Grid *gridIn, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal)
{
    StructGrid *grid = StructGridCast(gridIn);

    int *cellTypeContainer = grid->GetCellTypeContainer();

    RDouble4D &dqTurbulent = dqProxy->GetField_STR();
    RDouble4D &deltaFlux = LUplusDQ->GetField_STR();

    Param_TurbSolverStruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    RDouble *dqNeighbor = new RDouble[nTurbulentEquation + 1];
    RDouble *prim       = new RDouble[nTurbulentEquation + 1];
    RDouble *rhs0       = new RDouble[nTurbulentEquation + 1];
    RDouble *dflux      = new RDouble[nTurbulentEquation + 1];
    RDouble *mat_abc    = new RDouble[nTurbulentEquation + 1];
    RDouble *dqOld      = new RDouble[nTurbulentEquation + 1];
    RDouble *Ux         = new RDouble[nTurbulentEquation + 1];

    RDouble4D &spectumTurbulent  = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("spec_turb"));
    RDouble4D &residualTurbulent = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_turb"));
    RDouble5D &matrixTurbulentLeft    = *reinterpret_cast<RDouble5D *> (grid->GetDataPtr("mat_turbl"));

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    int nDim = GetDim();

	Range M(0, nTurbulentEquation - 1);
    int mst = M.first();
    int med = M.last();

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                for (int m = mst; m <= med; ++ m)
                {
                    rhs0[m] = 0.0;
                }

                for (int m = mst; m <= med; ++ m)
                {
                    //! Back up the old dq, to compute the convergence.
                    //! it is initialized by zero or res in the first sweep,
                    //! then it is updated by backward sweep, that is dq*.
                    dqOld[m] = dqTurbulent(i, j, k, m);

                    //! Here, the deltaFlux is computed in Upper backward sweep!
                    //! res: b      deltaFlux: Ux, which is computed in backward sweep.
                    //! the dq is not the real dq, now.
                    //! Although the 'dq' changed can not the right one, it dosen't matter, since 
                    //! the following only using the lower neighbor cells, whose 'dq' has been updated.
                    //! It is convenient for precondition transform.
                    RDouble &dq = dqTurbulent(i, j, k, m);
                    dq = residualTurbulent(i, j, k, m);

                    Ux[m] = deltaFlux(i, j, k, m);

                    //! Then reset it to zero to store the Lower forward sweep!
                    RDouble &df = deltaFlux(i, j, k, m);
                    df = 0.0;
                }

                for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
                {
                    int il1, jl1, kl1;
                    GetNsurfIndex(iSurface, il1, jl1, kl1);

                    int il = i - il1;
                    int jl = j - jl1;
                    int kl = k - kl1;

                    for (int m = mst; m <= med; ++ m)
                    {
                        dqNeighbor[m] = dqTurbulent(il, jl, kl, m);
                    }

                    for (int m = mst; m <= med; ++ m)
                    {
                        mat_abc[m] = matrixTurbulentLeft(il, jl, kl, m, iSurface);
                    }

                    Turb_MxDQ(mat_abc, dqNeighbor, nTurbulentEquation, dflux);

                    for (int m = 0; m < nTurbulentEquation; ++ m)
                    {
                        rhs0[m] += dflux[m];
                    }
                }
                
                for (int m = mst; m <= med; ++ m)
                {
                    //! dq = { (b - Ux) - Lx } / D.
                    //! Note: the 'dq' has actually initialized by the residual at the beginning of this function.
                    //! the 'rhs0' is the total Delta Flux of the neighbors.            
                    //! rhs0: Lx    diagonal: D.
                    dqTurbulent(i, j, k, m) = (dqTurbulent(i, j, k, m) - Ux[m] - rhs0[m]) / spectumTurbulent(i, j, k, m);

                    //! Store the lower forward sweep delta-flux, which will be used in the backward sweep.
                    deltaFlux(i, j, k, m) += rhs0[m];

                    sweepNormal += SQR(dqTurbulent(i, j, k, m) - dqOld[m]);
                }
            }
        }
    }

    delete []dqNeighbor;
    delete []prim;
    delete []rhs0;
    delete []dflux;
    delete []mat_abc;
    delete []dqOld;
    delete []Ux;
}

void TurbSolverStr::SolveLUSGSBackward(Grid *gridIn, FieldProxy * dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal)
{
    StructGrid *grid = StructGridCast(gridIn);

    int *cellTypeContainer = grid->GetCellTypeContainer();

    RDouble4D &dqTurbulent = dqProxy->GetField_STR();
    RDouble4D &deltaFlux = LUplusDQ->GetField_STR();

    Param_TurbSolverStruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    RDouble *dqNeighbor = new RDouble[nTurbulentEquation + 1];
    RDouble *dqOld      = new RDouble[nTurbulentEquation + 1];
    RDouble *prim       = new RDouble[nTurbulentEquation + 1];
    RDouble *rhs0       = new RDouble[nTurbulentEquation + 1];
    RDouble *dflux      = new RDouble[nTurbulentEquation + 1];
    RDouble *mat_abc    = new RDouble[nTurbulentEquation + 1];
    RDouble *Lx         = new RDouble[nTurbulentEquation + 1];

    RDouble4D &spectumTurbulent  = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("spec_turb"));
    RDouble4D &residualTurbulent = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("res_turb"));
    RDouble5D &matrixTurbulentRight   = *reinterpret_cast<RDouble5D *> (grid->GetDataPtr("mat_turbr"));

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    Range M(0, nTurbulentEquation - 1);

    int nDim = GetDim();
    int il, jl, kl;

    int mst = M.first();
    int med = M.last();

    for (int k = ked; k >= kst; -- k)
    {
        for (int j = jed; j >= jst; -- j)
        {
            for (int i = ied; i >= ist; -- i)
            {
                for (int m = mst; m <= med; ++ m)
                {
                    rhs0[m] = 0.0;
                }

                for (int m = mst; m <= med; ++ m)
                {
                    //! Back up the old dq, to compute the convergence.
                    //! The 'dq' is dq*, which has been updated in forward.
                    dqOld[m] = dqTurbulent(i, j, k, m);

                    //! The dq is not the real dq, now.
                    //! it is convenient for precondition transform.
                    dqTurbulent(i, j, k, m) = residualTurbulent(i, j, k, m);

                    Lx[m] = deltaFlux(i, j, k, m);

                    deltaFlux(i, j, k, m) = 0.0;
                }

                for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
                {
                    int il1, jl1, kl1;
                    GetNsurfIndex(iSurface, il1, jl1, kl1);

                    il = i + il1;
                    jl = j + jl1;
                    kl = k + kl1;

                    for (int m = mst; m <= med; ++ m)
                    {
                        dqNeighbor[m] = dqTurbulent(il, jl, kl, m);
                    }

                    for (int m = mst; m <= med; ++ m)
                    {
                        mat_abc[m] = matrixTurbulentRight(il, jl, kl, m, iSurface);
                    }

                    Turb_MxDQ(mat_abc, dqNeighbor, nTurbulentEquation, dflux);

                    for (int m = mst; m <= med; ++ m)
                    {
                        rhs0[m] += dflux[m];
                    }
                }
                
                for (int m = mst; m <= med; ++ m)
                {
                    //! Note: the 'dq' has been updated by the forward sweep.
                    //! The 'rhs0' is the total Delta Flux of the neighbors.
                    //! x = {(b - LX) - Ux} / D.
                    //! rhs0: Ux.    diagonal: D.
                    dqTurbulent(i, j, k, m) = (dqTurbulent(i, j, k, m) - Lx[m] - rhs0[m]) / spectumTurbulent(i, j, k, m);

                    //! Store the upper backward sweep delta-flux, which will be used in the forward sweep.
                    deltaFlux(i, j, k, m) += rhs0[m];

                    sweepNormal += SQR(dqTurbulent(i, j, k, m) - dqOld[m]);
                }
            }
        }
    }

    delete []dqNeighbor;
    delete []dqOld;
    delete []prim;
    delete []rhs0;
    delete []dflux;
    delete []mat_abc;
    delete []Lx;
}

void TurbSolverStr::UpdateFlowField(Grid *gridIn, FieldProxy *qProxy, FieldProxy *dqProxy)
{
    StructGrid *grid = StructGridCast(gridIn);

    RDouble4D &qTurbulent  = qProxy->GetField_STR();
    RDouble4D &dqTurbulent = dqProxy->GetField_STR();

    Param_TurbSolverStruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    RDouble4D &qpmv = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q" ));

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    Range M(0, nTurbulentEquation - 1);

    using namespace IDX;
    if (nTurbulentEquation == 1)
    {
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    qTurbulent(i, j, k,ISA) = MAX(qTurbulent(i, j, k, ISA) + dqTurbulent(i, j, k, ISA), zero);
                }
            }
        }
    }
    else if (nTurbulentEquation == 2)
    {
        RDouble   *freeStreamTurbVar =  parameters->GetFreeStreamTurbVar();
        RDouble kelim = 1.0e-5 * freeStreamTurbVar[IKE];
        RDouble kwlim = 1.0e-5 * freeStreamTurbVar[IKW];

        RDouble ke,kw,dke,dkw,rho;

        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {                  
                    rho = qpmv(i, j, k,IR);
                    ke  = qTurbulent(i, j, k, IKE);
                    kw  = qTurbulent(i, j, k, IKW); 
                    dke = dqTurbulent(i, j, k, IKE);
                    dkw = dqTurbulent(i, j, k, IKW);
                    ke += dke / rho;
                    kw += dkw / rho;
                    ke = MAX(ke, kelim);
                    kw = MAX(kw, kwlim);
                    qTurbulent(i, j, k, IKE) = ke;
                    qTurbulent(i, j, k, IKW) = kw;
                }
            }
        }
    }
}

void TurbSolverStr::ComputeViscousCoeff(Grid *grid)
{
        Viscosity(grid);
}

void TurbSolverStr::Viscosity(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    RDouble3D &wallDistant = *grid->GetWallDist();

    RDouble4D &area   = *(grid->GetFaceArea()   );
    RDouble3D &volume = *(grid->GetCellVolume() );

    RDouble4D &qLaminar         = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q" ));
    RDouble3D &viscousLaminar   = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("visl"));
    RDouble3D &viscousTurbulent = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("vist"));
    RDouble4D &qTurbulent       = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb" ));
    RDouble3D &SST_F2           = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("SST_F2"));

    Param_TurbSolverStruct *parameters = GetControlParameters();
    RDouble eddyViscosityLimit = parameters->GetEddyViscosityLimit();

    int monitorVistmax = GlobalDataBase::GetIntParaFromDB("monitor_vistmax");

    string viscousName = parameters->GetViscousName();

    RDouble refReNumber = parameters->GetRefReNumber();
    RDouble reynoldsSquare = refReNumber * refReNumber;

    RDouble viscousTurbulentMaximum = 0.0;
    RDouble viscousTurbulentMinimum = 1.0e30;
    RDouble logsmall = log(SMALL);

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    int imax = 1;
    int jmax = 1;
    int kmax = 1;

    using namespace IDX;

    if (viscousName.substr(0,6) == "1eq-sa")
    {
        RDouble SA_cv1_cube = parameters->GetSA_cv1_cube();

        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    RDouble ld   = qTurbulent(i, j, k, ISA) * qLaminar(i, j, k, IR) / (viscousLaminar(i, j, k) + SMALL);
                    RDouble ld3  = ld * ld * ld;
                    RDouble fv1  = ld3 / (ld3 + SA_cv1_cube);
                    viscousTurbulent(i, j, k) = qLaminar(i, j, k, IR) * fv1 * qTurbulent(i, j, k, ISA);
                    viscousTurbulent(i, j, k) = MIN(eddyViscosityLimit, viscousTurbulent(i, j, k));
                    if (viscousTurbulentMaximum < viscousTurbulent(i, j, k))
                    {
                        viscousTurbulentMaximum =  viscousTurbulent(i, j, k);
                        imax = i;
                        jmax = j;
                        kmax = k;
                    }
                }
            }
        }
    }
    else if (viscousName.substr(0,6) == "2eq-kw")
    {
        //! Two equation models.
        if (viscousName.substr(0,17) == "2eq-kw-menter-sst")
        {
            //! Original Menter k-w-SST model.
            RDouble SST_a1    = parameters->GetSST_a1();
            RDouble SST_betaStar = parameters->GetSST_betaStar();

            Range II, JJ, KK;
            GetRange(ni, nj, nk, -2, 1, II, JJ, KK);

            RDouble4D &gradUVWTCellCenterX = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterX"));
            RDouble4D &gradUVWTCellCenterY = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterY"));
            RDouble4D &gradUVWTCellCenterZ = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterZ"));

            for (int k = kst; k <= ked; ++ k)
            {
                for (int j = jst; j <= jed; ++ j)
                {
                    for (int i = ist; i <= ied; ++ i)
                    {
                        RDouble rho = ABS(qLaminar(i, j, k, IR)) + SMALL;
                        RDouble ke  = qTurbulent(i, j, k, IDX::IKE);
                        RDouble kw  = qTurbulent(i, j, k, IDX::IKW);
                        RDouble ds  = wallDistant(i, j, k);
                        RDouble ds2 = ds * ds;

                        RDouble dudx = gradUVWTCellCenterX(i, j, k, 0);
                        RDouble dudy = gradUVWTCellCenterY(i, j, k, 0);
                        RDouble dudz = gradUVWTCellCenterZ(i, j, k, 0);

                        RDouble dvdx = gradUVWTCellCenterX(i, j, k, 1);
                        RDouble dvdy = gradUVWTCellCenterY(i, j, k, 1);
                        RDouble dvdz = gradUVWTCellCenterZ(i, j, k, 1);

                        RDouble dwdx = gradUVWTCellCenterX(i, j, k, 2);
                        RDouble dwdy = gradUVWTCellCenterY(i, j, k, 2);
                        RDouble dwdz = gradUVWTCellCenterZ(i, j, k, 2);

                        //RDouble vort  = sqrt((dwdy - dvdz) * (dwdy - dvdz) +
                        //    (dudz - dwdx) * (dudz - dwdx) +
                        //    (dvdx - dudy) * (dvdx - dudy));

                        RDouble s11 = dudx; 
                        RDouble s22 = dvdy; 
                        RDouble s33 = dwdz;
                        RDouble s12 = half * (dudy + dvdx); 
                        RDouble s13 = half * (dudz + dwdx); 
                        RDouble s23 = half * (dvdz + dwdy);

                        RDouble sij2 = two * (s11*s11 + s22*s22 + s33*s33 + two*(s12*s12 + s13*s13 + s23*s23));  //! modulus of S.
                        RDouble Strain = sqrt(sij2);

                        RDouble part1 = two * sqrt(ke) / (SST_betaStar * kw * ds * refReNumber);
                        RDouble part2 = 500.0 * viscousLaminar(i, j, k) / (rho * kw * ds2 * reynoldsSquare);
                        RDouble arg2  = MAX(part1, part2);
                        RDouble f2    = tanh(static_cast<RDouble>(arg2 * arg2));
                        SST_F2(i, j, k) = f2;

                        RDouble &turbulentViscosity = viscousTurbulent(i, j, k);
                        turbulentViscosity = rho * ke / MAX(kw, Strain * f2 / (SST_a1 * refReNumber));
                        turbulentViscosity = MIN(eddyViscosityLimit, turbulentViscosity);
                        turbulentViscosity = MAX(1.E-5,turbulentViscosity);

                        if (viscousTurbulentMaximum < viscousTurbulent(i, j, k))
                        {
                            viscousTurbulentMaximum = turbulentViscosity;
                            imax = i;
                            jmax = j;
                            kmax = k;
                        }
                    }
                }
            }
        }
    }

    //GhostCell3D(viscousTurbulent, ni, nj, nk);
    //grid->GhostCell3DExceptInterface(viscousTurbulent);

    RDouble vist_max = viscousTurbulentMaximum;
    RDouble vist_min = viscousTurbulentMinimum;

    grid->UpdateData("vist_max", &vist_max, PHDOUBLE, 1);
    grid->UpdateData("vist_min", &vist_min, PHDOUBLE, 1);

#ifdef PH_PARALLEL
    using namespace PHMPI;
    int myid = GetCurrentProcessorID();
    if (myid == server)
    {
#endif
        if (monitorVistmax)
        {
            cout << "vistmax = " << viscousTurbulentMaximum << " " << imax << " " << jmax << " " << kmax << "\n";
        }
#ifdef PH_PARALLEL
    }
#endif
}

void TurbSolverStr::Boundary(Grid *gridIn)
{
    using namespace PHENGLEI;
    StructGrid *grid = StructGridCast(gridIn);

    Param_TurbSolverStruct *parameters = GetControlParameters();
    int wallFunctionType = parameters->GetWallFunctionType();

    CompositeBCRegion *compositeBCRregion = grid->GetCompositeBCRegion();
    int nBCRegion = compositeBCRregion->GetnBCRegion();

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        BCRegion *BCRegion = compositeBCRregion->GetBCRegion(iBCRegion);
        int BCType = BCRegion->GetBCType();

        if (BCType < 0)
        {
            continue;
        }
        else if (BCType == EXTRAPOLATION)
        {
            OutFlowBC(grid, BCRegion);
        }
        else if (IsWall(BCType))
        {
            if (wallFunctionType ==0)
            {
                VisWall(grid, BCRegion);
            }
            else if(wallFunctionType == 1)
            { 
                VisWallWithWallFunctionStandard(grid, BCRegion);
            }
            else if (wallFunctionType == 2)
            {
                VisWallWithWallFunctionPAB3D(grid, BCRegion);
            }
        }
        else if (BCType == SYMMETRY)
        {
            SymmetryBC(grid, BCRegion);
        }
        else if (BCType == FARFIELD)
        {
            FarFieldBC(grid, BCRegion);
        }
        else if (BCType == INFLOW)
        {
            InFlowBC(grid, BCRegion);
        }
        else if (BCType == PRESSURE_INLET)
        {
            InFlowBC(grid, BCRegion);
        }
        else if (BCType == MASS_FLOW_INLET)
        {
            InFlowBC(grid, BCRegion);
        }
        else if (BCType == OUTFLOW)
        {
            OutFlowBC(grid, BCRegion);
        }
        else if (BCType == PRESSURE_OUTLET)
        {
            OutFlowBC(grid, BCRegion);
        }
        else if (BCType == MASS_FLOW_OUTLET)
        {
            OutFlowBC(grid, BCRegion);
        }
        else if (BCType == POLE || BCType / 10 == POLE)
        {
            OutFlowBC(grid, BCRegion);
        }
        else if (BCType == EXTERNAL_BC)
        {
            OutFlowBC(grid, BCRegion);
        }
        //! Deal in symmetry wall method if there is no boundary ID specified.
        else if (ifIsJetBC(BCType))
        {
            SymmetryBC(grid, BCRegion);
        }
        else
        {
            ostringstream oss;
            oss << "Error : Illegal BCtype ID " << BCType << endl;
            TK_Exit::ExceptionExit(oss);
        }
    }
    CornerPoint(grid);
}

void TurbSolverStr::CornerPoint(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Param_TurbSolverStruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    RDouble4D &qTurbulent = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb"));
    FillCornerPoint3D(qTurbulent, ni, nj, nk, nTurbulentEquation);
}

void TurbSolverStr::InFlowBC(Grid *grid, BCRegion *BCRegion)
{
    Param_TurbSolverStruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();
    RDouble4D &qTurbulent = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb"));
    RDouble   *freeStreamTurbVar =  parameters->GetFreeStreamTurbVar();
    RDouble3D &viscousTurbulent = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("vist"));
    RDouble freeStreamViscosity = parameters->GetFreeStreamViscosity();

    int ist, ied, jst, jed, kst, ked;
    BCRegion->GetIJKRegion(ist, ied, jst, jed, kst, ked);

    int it, jt, kt;

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                for (int layer = 1; layer <= 2; ++ layer)
                {
					BCRegion->GetGhostCellIndex(i, j, k, it, jt, kt, layer);

                    for (int m = 0; m < nTurbulentEquation; ++ m)
                    {
                        qTurbulent(it, jt, kt, m) = freeStreamTurbVar[m];
                    }
                    viscousTurbulent(it, jt, kt) = freeStreamViscosity;
                }
            }
        }
    }
}

void TurbSolverStr::OutFlowBC(Grid *gridIn, BCRegion *BCRegion)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Param_TurbSolverStruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    RDouble4D &qTurbulent = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb"));
    RDouble3D &viscousTurbulent = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("vist"));

    int ist, ied, jst, jed, kst, ked;
    BCRegion->GetIJKRegion(ist, ied, jst, jed, kst, ked);

	int is1, js1, ks1, is2, js2, ks2;
	int it1, jt1, kt1, it2, jt2, kt2;

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
				BCRegion->GetInsideCellIndex(i, j, k, is1, js1, ks1, 1);
				BCRegion->GetInsideCellIndex(i, j, k, is2, js2, ks2, 2);

                RestrictIndex(is1, js1, ks1, ni-1, nj-1, nk-1);
                RestrictIndex(is2, js2, ks2, ni-1, nj-1, nk-1);

				BCRegion->GetGhostCellIndex(i, j, k, it1, jt1, kt1, 1);
				BCRegion->GetGhostCellIndex(i, j, k, it2, jt2, kt2, 2);

                    for (int m = 0; m < nTurbulentEquation; ++ m)
                    {
                        qTurbulent(it1, jt1, kt1, m) = two * qTurbulent(is1, js1, ks1, m) - qTurbulent(is2, js2, ks2, m);
                        qTurbulent(it2, jt2, kt2, m) = two * qTurbulent(it1, jt1, kt1, m) - qTurbulent(is1, js1, ks1, m);
                    }
                    viscousTurbulent(it1, jt1, kt1) = viscousTurbulent(is1, js1, ks1);
                    viscousTurbulent(it2, jt2, kt2) = viscousTurbulent(it1, jt1, kt1);
            }
        }
    }
}

void TurbSolverStr::SymmetryBC(Grid *gridIn, BCRegion *BCRegion)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Param_TurbSolverStruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    RDouble4D &qTurbulent       = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb"));
    RDouble3D &viscousTurbulent = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("vist"));

    int ist, ied, jst, jed, kst, ked;
    BCRegion->GetIJKRegion(ist, ied, jst, jed, kst, ked);

	int is, js, ks, it, jt, kt;

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                for (int layer = 1; layer <= 2; ++ layer)
                {
					BCRegion->GetInsideCellIndex(i, j, k, is, js, ks, layer);
                    RestrictIndex(is, js, ks, ni-1, nj-1, nk-1);
					BCRegion->GetGhostCellIndex(i, j, k, it, jt, kt, layer);

                    for (int m = 0; m < nTurbulentEquation; ++ m)
                    {
                        qTurbulent(it, jt, kt, m) = qTurbulent(is, js, ks, m);
                    }
                    viscousTurbulent(it, jt, kt) = viscousTurbulent(is, js, ks);
                }
            }
        }
    }
}

void TurbSolverStr::VisWall(Grid *gridIn, BCRegion *BCRegion)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Param_TurbSolverStruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    RDouble4D &qLaminar           = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));
    RDouble4D &qTurbulent         = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb"));
    RDouble3D &viscousLaminars    = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("visl"));
    RDouble3D &viscousTurbulent   = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("vist"));

    RDouble3D &wallDistant = *StructGridCast(grid)->GetWallDist();

    RDouble refReNumber = parameters->GetRefReNumber();

    int ist, ied, jst, jed, kst, ked;
    BCRegion->GetIJKRegion(ist, ied, jst, jed, kst, ked);

    const RDouble beta1 = 0.075;

    using namespace IDX;

	int is1, js1, ks1, is2, js2, ks2;
	int it1, jt1, kt1, it2, jt2, kt2;

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
				BCRegion->GetInsideCellIndex(i, j, k, is1, js1, ks1, 1);
				BCRegion->GetInsideCellIndex(i, j, k, is2, js2, ks2, 2);

                RestrictIndex(is1, js1, ks1, ni-1, nj-1, nk-1);
                RestrictIndex(is2, js2, ks2, ni-1, nj-1, nk-1);

				BCRegion->GetGhostCellIndex(i, j, k, it1, jt1, kt1, 1);
				BCRegion->GetGhostCellIndex(i, j, k, it2, jt2, kt2, 2);

                if (nTurbulentEquation == 1)
                {
                    qTurbulent(it1, jt1, kt1, IDX::ISA) = 0.0;
                    qTurbulent(it2, jt2, kt2, IDX::ISA) = 0.0;
                }
                else if (nTurbulentEquation == 2)
                {
                    RDouble rho   = qLaminar(i, j, k, IR);

                    RDouble dist  = wallDistant(i, j, k);
                    RDouble dist2 = dist * dist;
                    RDouble muLaminar = viscousLaminars(i, j, k);

                    RDouble kw_wall = 60.0 * muLaminar / (rho * beta1 * dist2 * refReNumber * refReNumber);

                    qTurbulent(it1, jt1, kt1, IDX::IKE) = - qTurbulent(is1, js1, ks1, IDX::IKE);
                    qTurbulent(it2, jt2, kt2, IDX::IKE) =   qTurbulent(it1, jt1, kt1, IDX::IKE);
          
                    qTurbulent(it1, jt1, kt1, IDX::IKW) = 2.0 * kw_wall - qTurbulent(is1, js1, ks1, IDX::IKW);
                    qTurbulent(it2, jt2, kt2, IDX::IKW) = qTurbulent(it1, jt1, kt1, IDX::IKW);
                }
                viscousTurbulent(it1, jt1, kt1) = -viscousTurbulent(is1, js1, ks1);
                viscousTurbulent(it2, jt2, kt2) = -viscousTurbulent(is2, js2, ks2);
            }
        }
    }
}

void TurbSolverStr::VisWallWithWallFunctionStandard(Grid* gridIn, BCRegion* BCRegion)
{
    StructGrid* grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Param_TurbSolverStruct* parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    RDouble4D& qLaminar = *reinterpret_cast<RDouble4D*> (grid->GetDataPtr("q"));
    RDouble4D& qTurbulent = *reinterpret_cast<RDouble4D*> (grid->GetDataPtr("q_turb"));
    RDouble3D& viscousLaminars = *reinterpret_cast<RDouble3D*> (grid->GetDataPtr("visl"));
    RDouble3D& viscousTurbulent = *reinterpret_cast<RDouble3D*> (grid->GetDataPtr("vist"));

    RDouble4D& faceNormalX = *reinterpret_cast<RDouble4D*> (grid->GetFaceNormalX());
    RDouble4D& faceNormalY = *reinterpret_cast<RDouble4D*> (grid->GetFaceNormalY());
    RDouble4D& faceNormalZ = *reinterpret_cast<RDouble4D*> (grid->GetFaceNormalZ());

    RDouble3D& wallDistant = *StructGridCast(grid)->GetWallDist();

    RDouble refReNumber = parameters->GetRefReNumber();

    int ist, ied, jst, jed, kst, ked;
    BCRegion->GetIJKRegion(ist, ied, jst, jed, kst, ked);

    using namespace IDX;

    int iSurface = BCRegion->GetFaceDirection() + 1;

    const RDouble E_ = 9.793;
    const RDouble beta1 = 0.075;
    const RDouble kappa1 = 0.4178;
    const RDouble yPluslimit = 11.225;

    int is1, js1, ks1, is2, js2, ks2;
    int it1, jt1, kt1, it2, jt2, kt2;

    for (int k = kst; k <= ked; ++k)
    {
        for (int j = jst; j <= jed; ++j)
        {
            for (int i = ist; i <= ied; ++i)
            {

                BCRegion->GetInsideCellIndex(i, j, k, is1, js1, ks1, 1);
                BCRegion->GetInsideCellIndex(i, j, k, is2, js2, ks2, 2);

                RestrictIndex(is1, js1, ks1, ni - 1, nj - 1, nk - 1);
                RestrictIndex(is2, js2, ks2, ni - 1, nj - 1, nk - 1);

                BCRegion->GetGhostCellIndex(i, j, k, it1, jt1, kt1, 1);
                BCRegion->GetGhostCellIndex(i, j, k, it2, jt2, kt2, 2);

                RDouble rhoInside = qLaminar(is1, js1, ks1, IDX::IR);
                RDouble rhoGhost  = qLaminar(it1, jt1, kt1, IDX::IR);
                RDouble densityWall = 0.5 * (rhoInside + rhoGhost);
                RDouble up = qLaminar(is1, js1, ks1, IDX::IU);
                RDouble vp = qLaminar(is1, js1, ks1, IDX::IV);
                RDouble wp = qLaminar(is1, js1, ks1, IDX::IW);
                RDouble laminarViscosity = viscousLaminars(is1, js1, ks1);
                RDouble turbulentViscosity = viscousTurbulent(is1, js1, ks1);

                RDouble xfn = faceNormalX(is1, js1, ks1, iSurface);
                RDouble yfn = faceNormalY(is1, js1, ks1, iSurface);
                RDouble zfn = faceNormalZ(is1, js1, ks1, iSurface);

                RDouble vpn = up * xfn + vp * yfn + wp * zfn;
                RDouble vpt = sqrt(up * up + vp * vp + wp * wp - vpn * vpn);
                RDouble small = 1.0e-15;
                RDouble uvw = max(vpt, small);

                RDouble wallDistance = wallDistant(is1, js1, ks1);

                //! Obtain the reynold number.
                RDouble refReNumber = parameters->GetRefReNumber();

                RDouble taow = laminarViscosity * uvw / wallDistance / refReNumber;
                RDouble utao = sqrt(taow / densityWall);

                RDouble yPlus;
                RDouble epslon = 1.0;
                int iter = 0;
                while (epslon > 0.001 && ++iter < 20)
                {
                    yPlus = refReNumber * wallDistance * densityWall * utao / laminarViscosity;
                    yPlus = MAX(yPlus, yPluslimit);
                    RDouble utemp = utao;

                    utao = kappa1 * uvw / (log(E_ * yPlus));

                    epslon = abs(utao - utemp) / (abs(utemp) + 1.0e-20);
                }
                taow = densityWall * utao * utao;
                RDouble viscousTurbulenceWall = laminarViscosity * (taow / (uvw + small) / laminarViscosity * wallDistance * refReNumber - 1.0);
                viscousTurbulent(it1, jt1, kt1) = 2.0 * viscousTurbulenceWall - viscousTurbulent(is1, js1, ks1);
                viscousTurbulent(it2, jt2, kt2) = viscousTurbulent(it1, jt1, kt1);

                if (nTurbulentEquation == 1)
                {
                    qTurbulent(it1, jt1, kt1, IDX::ISA) = -qTurbulent(is1, js1, ks1, IDX::ISA);
                    qTurbulent(it2, jt2, kt2, IDX::ISA) =  qTurbulent(it1, jt1, kt1, IDX::ISA);
                }
                else if (nTurbulentEquation == 2)
                {
                    RDouble distanceSquare = wallDistance * wallDistance;

                    RDouble omgi = 60.0 * laminarViscosity / (rhoInside * beta1 * distanceSquare * refReNumber * refReNumber);
                    //RDouble omgo = utao / (0.126 * wallDistance * refReNumber);
                    RDouble omgo = 0.0;
                    RDouble kwWall = sqrt(omgi * omgi + omgo * omgo);

                    RDouble qTurbulentFirstCell = kwWall * viscousTurbulenceWall / densityWall;
                    //qTurbulent(is1, js1, ks1, IDX::IKE) = qTurbulentFirstCell;
                    qTurbulent(it1, jt1, kt1, IDX::IKE) = -qTurbulent(is1, js1, ks1, IDX::IKE);
                    qTurbulent(it2, jt2, kt2, IDX::IKE) = qTurbulent(it1, jt1, kt1, IDX::IKE);

                    //qTurbulent(is1, js1, ks1, IDX::IKW) = kwWall;
                    qTurbulent(it1, jt1, kt1, IDX::IKW) = 2.0 * kwWall - qTurbulent(is1, js1, ks1, IDX::IKW);
                    //qTurbulent(it1, jt1, kt1, IDX::IKW) = kwWall;
                    //qTurbulent(it1, jt1, kt1, IDX::IKW) = -qTurbulent(it1, jt1, kt1, IDX::IKE) * rhoGhost / viscousTurbulent(it1, jt1, kt1);
                    qTurbulent(it2, jt2, kt2, IDX::IKW) = qTurbulent(it1, jt1, kt1, IDX::IKW);
                    //qTurbulent(it1, jt1, kt1, IDX::IKE) = qTurbulent(it1, jt1, kt1, IDX::IKW) * viscousTurbulent(it1, jt1, kt1) / rhoGhost;
                    //qTurbulent(it2, jt2, kt2, IDX::IKE) = qTurbulent(it2, jt2, kt2, IDX::IKW) * viscousTurbulent(it2, jt2, kt2) / rhoGhost;
                }
            }
        }
    }
}

void TurbSolverStr::VisWallWithWallFunctionPAB3D(Grid* gridIn, BCRegion* BCRegion)
{
    StructGrid* grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    using namespace IDX;

    RDouble paba10[7] = { 2.354039, 0.1179840, -4.2899192e-04, 2.0404148e-06,-5.1775775e-09, 6.2687308e-12, -2.916958e-15 };
    RDouble paba11[5] = { 5.777191, 6.8756983e-02, -7.1582745e-06, 1.5594904e-09, -1.4865778e-13 };
    RDouble paba12[5] = { 31.08654, 5.0429072e-02, -2.0072314e-8 };
    RDouble Ylim = 10.0;
    RDouble beta1 = 0.075;

    Param_TurbSolverStruct* parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    RDouble4D& qLaminar = *reinterpret_cast<RDouble4D*> (grid->GetDataPtr("q"));
    RDouble4D& qTurbulent = *reinterpret_cast<RDouble4D*> (grid->GetDataPtr("q_turb"));
    RDouble4D& temperature = *reinterpret_cast<RDouble4D*> (grid->GetDataPtr("t"));
    RDouble3D& viscousLaminars = *reinterpret_cast<RDouble3D*> (grid->GetDataPtr("visl"));
    RDouble3D& viscousTurbulent = *reinterpret_cast<RDouble3D*> (grid->GetDataPtr("vist"));

    RDouble4D& faceNormalX = *reinterpret_cast<RDouble4D*> (grid->GetFaceNormalX());
    RDouble4D& faceNormalY = *reinterpret_cast<RDouble4D*> (grid->GetFaceNormalY());
    RDouble4D& faceNormalZ = *reinterpret_cast<RDouble4D*> (grid->GetFaceNormalZ());

    RDouble3D& wallDistant = *StructGridCast(grid)->GetWallDist();

    RDouble refReNumber = parameters->GetRefReNumber();
    RDouble refMachNumber = parameters->GetRefMachNumber();
    RDouble refGama = parameters->GetRefGama();

    int ist, ied, jst, jed, kst, ked;
    BCRegion->GetIJKRegion(ist, ied, jst, jed, kst, ked);

    using namespace IDX;

    int iSurface = BCRegion->GetFaceDirection() + 1;

    for (int k = kst; k <= ked; ++k)
    {
        for (int j = jst; j <= jed; ++j)
        {
            for (int i = ist; i <= ied; ++i)
            {
                int is1, js1, ks1, is2, js2, ks2;
                int it1, jt1, kt1, it2, jt2, kt2;

                BCRegion->GetInsideCellIndex(i, j, k, is1, js1, ks1, 1);
                BCRegion->GetInsideCellIndex(i, j, k, is2, js2, ks2, 2);

                RestrictIndex(is1, js1, ks1, ni - 1, nj - 1, nk - 1);
                RestrictIndex(is2, js2, ks2, ni - 1, nj - 1, nk - 1);

                BCRegion->GetGhostCellIndex(i, j, k, it1, jt1, kt1, 1);
                BCRegion->GetGhostCellIndex(i, j, k, it2, jt2, kt2, 2);

                RDouble rhoInside = qLaminar(is1, js1, ks1, IDX::IR);
                RDouble rhoGhost = qLaminar(it1, jt1, kt1, IDX::IR);
                RDouble up = qLaminar(is1, js1, ks1, IDX::IU);
                RDouble vp = qLaminar(is1, js1, ks1, IDX::IV);
                RDouble wp = qLaminar(is1, js1, ks1, IDX::IW);
                RDouble laminarViscosity = viscousLaminars(is1, js1, ks1);
                RDouble turbulentViscosity = viscousTurbulent(is1, js1, ks1);

                RDouble xfn = faceNormalX(is1, js1, ks1, iSurface);
                RDouble yfn = faceNormalY(is1, js1, ks1, iSurface);
                RDouble zfn = faceNormalZ(is1, js1, ks1, iSurface);

                RDouble vpn = up * xfn + vp * yfn + wp * zfn;
                RDouble vpt = sqrt(up * up + vp * vp + wp * wp - vpn * vpn);
                RDouble small = 1.0e-15;
                RDouble uvw = max(vpt, small);

                //RDouble uu = sqrt(up * up + vp * vp + wp * wp);
                RDouble& uu = uvw;

                RDouble wallDistance = wallDistant(is1, js1, ks1);

                RDouble refReNumber = parameters->GetRefReNumber(); //! Obtain the reynold number.

                RDouble rc = rhoInside * refReNumber * wallDistance * uu / laminarViscosity;
                RDouble xnplus = 0.0;
                //! The original one here is yPlus.
                if (rc <= 20.24)
                {
                    xnplus = sqrt(rc);
                }
                else if (rc <= 435.0)
                {
                    xnplus = paba10[0] + paba10[1] * rc + paba10[2] * rc * rc + paba10[3] * rc * rc * rc + paba10[4] * rc * rc * rc * rc + paba10[5] * rc * rc * rc * rc * rc + paba10[6] * rc * rc * rc * rc * rc * rc;
                }
                else if (rc <= 4000.0)
                {
                    xnplus = paba11[0] + paba11[1] * rc + paba11[2] * rc * rc + paba11[3] * rc * rc * rc + paba11[4] * rc * rc * rc * rc;
                }
                else if (rc > 4000.0)
                {
                    xnplus = paba12[0] + paba12[1] * rc + paba12[2] * rc * rc;
                }
                RDouble xnplussav = xnplus;

                //if (xnplus >= Ylim)
                {
                    //! Newton iteration to solve for nplus, assuming it is in log region:
                    for (int iter = 0; ; ++iter)
                    {
                        RDouble f = rc / xnplus - 2.44 * (log(xnplus)) - 5.2;
                        RDouble dfdn = -rc / (xnplus * xnplus) - 2.44 / xnplus;
                        RDouble delta = -f / dfdn;
                        xnplus = fabs(xnplus + delta);
                        if (iter > 10)
                        {
                            //! Revert back to approx series soln if Newton iteration fails.
                            xnplus = xnplussav;
                            break;
                        }
                        else if (abs(delta) < 1.e-3)
                        {
                            break;
                        }
                    }
                }

                RDouble viscousTurbulenceWall = laminarViscosity * (laminarViscosity * xnplus * xnplus / (wallDistance * rhoInside * uu * refReNumber) - 1.0);
                viscousTurbulenceWall = MAX(viscousTurbulenceWall, 0.0);
                viscousTurbulent(it1, jt1, kt1) = 2.0 * viscousTurbulenceWall - viscousTurbulent(is1, js1, ks1);
                viscousTurbulent(it2, jt2, kt2) = viscousTurbulent(it1, jt1, kt1);

                if (nTurbulentEquation == 1)
                {
                    qTurbulent(it1, jt1, kt1, IDX::ISA) = -qTurbulent(is1, js1, ks1, IDX::ISA);
                    qTurbulent(it2, jt2, kt2, IDX::ISA) =  qTurbulent(it1, jt1, kt1, IDX::ISA);
                }
                else if (nTurbulentEquation == 2)
                {
                    RDouble distanceSquare = wallDistance * wallDistance;

                    RDouble omgi = 60.0 * laminarViscosity / (rhoInside * beta1 * distanceSquare * refReNumber * refReNumber);
                    RDouble taow = laminarViscosity * uvw / wallDistance / refReNumber;
                    RDouble utao = sqrt(taow / rhoInside);
                    RDouble omgo = utao / (0.126 * wallDistance * refReNumber);
                    //RDouble omgo = 0.0;
                    RDouble kwWall = sqrt(omgi * omgi + omgo * omgo);

                    //qTurbulent(is1, js1, ks1, IDX::IKW) = kwWall;
                    qTurbulent(it1, jt1, kt1, IDX::IKW) = 2.0 * kwWall - qTurbulent(is1, js1, ks1, IDX::IKW);
                    qTurbulent(it2, jt2, kt2, IDX::IKW) = qTurbulent(it1, jt1, kt1, IDX::IKW);

                    qTurbulent(it1, jt1, kt1, IDX::IKE) = qTurbulent(it1, jt1, kt1, IDX::IKW) * viscousTurbulent(it1, jt1, kt1) / rhoGhost;
                    //qTurbulent(is1, js1, ks1, IDX::IKE) = qTurbulent(is1, js1, ks1, IDX::IKW) * viscousTurbulent(is1, js1, ks1) / rhoInside;
                    qTurbulent(it2, jt2, kt2, IDX::IKE) = qTurbulent(it1, jt1, kt1, IDX::IKE);
                }
            }
        }
    }
}

void TurbSolverStr::FarFieldBC(Grid *gridIn, BCRegion *BCRegion)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    RDouble4D &faceNormalX  = *(grid->GetFaceNormalX());
    RDouble4D &faceNormalY  = *(grid->GetFaceNormalY());
    RDouble4D &faceNormalZ  = *(grid->GetFaceNormalZ());

    Param_TurbSolverStruct *parameters = GetControlParameters();
    RDouble *freeStreamTurbVar = parameters->GetFreeStreamTurbVar();

    RDouble4D &qTurbulent       = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb"));
    RDouble4D &qLaminar         = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));
    RDouble3D &gamma            = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("gama"));
    RDouble3D &viscousTurbulent = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("vist"));

    int nm = GlobalDataBase::GetIntParaFromDB("nm");
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    RDouble refGama = GlobalDataBase::GetDoubleParaFromDB("refGama");
    RDouble referenceMachNumber = parameters->GetRefMachNumber();
    RDouble freeStreamViscosity = parameters->GetFreeStreamViscosity();

    using namespace IDX;

    RDouble *primitiveVariablesInfinity = reinterpret_cast<RDouble *> (GlobalDataBase::GetDataPtr("prim_inf" ));
    RDouble rFarfield = primitiveVariablesInfinity[IR];
    RDouble uFarfield = primitiveVariablesInfinity[IU];
    RDouble vFarfield = primitiveVariablesInfinity[IV];
    RDouble wFarfield = primitiveVariablesInfinity[IW];
    RDouble pFarfield = primitiveVariablesInfinity[IP];

    int ist, ied, jst, jed, kst, ked;
    BCRegion->GetIJKRegion(ist, ied, jst, jed, kst, ked);

    int iSurface = BCRegion->GetFaceDirection() + 1;
    int leftOrRightIndex    = BCRegion->GetFaceLeftOrRightIndex();

    RDouble *prims1 = new RDouble [nm];

    int is1, js1, ks1, is2, js2, ks2;
	int it1, jt1, kt1, it2, jt2, kt2;

    using namespace IDX;

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
				int in, jn, kn;
				BCRegion->GetBoundaryFaceIndex(i, j, k, in, jn, kn);
                RDouble nxs = leftOrRightIndex * faceNormalX(in, jn, kn, iSurface);
                RDouble nys = leftOrRightIndex * faceNormalY(in, jn, kn, iSurface);
                RDouble nzs = leftOrRightIndex * faceNormalZ(in, jn, kn, iSurface);

				BCRegion->GetInsideCellIndex(i, j, k, is1, js1, ks1, 1);
				BCRegion->GetInsideCellIndex(i, j, k, is2, js2, ks2, 2);

                RestrictIndex(is1, js1, ks1, ni-1, nj-1, nk-1);
                RestrictIndex(is2, js2, ks2, ni-1, nj-1, nk-1);

				BCRegion->GetGhostCellIndex(i, j, k, it1, jt1, kt1, 1);
				BCRegion->GetGhostCellIndex(i, j, k, it2, jt2, kt2, 2);

                for (int m = 0; m < nm; ++ m)
                {
                    prims1[m] = qLaminar(is1, js1, ks1, m);
                }

                RDouble rin = prims1[IR];
                RDouble uin = prims1[IU];
                RDouble vin = prims1[IV];
                RDouble win = prims1[IW];
                RDouble pin = prims1[IP];

                RDouble vno = nxs * uFarfield + nys * vFarfield + nzs * wFarfield;
                RDouble vni = nxs * uin + nys * vin + nzs * win;
                RDouble vei = sqrt(uin * uin + vin * vin + win * win);

                RDouble gama = gamma(is1, js1, ks1);

                RDouble coo = sqrt(ABS(refGama * pFarfield / rFarfield));
                RDouble cin = sqrt(ABS(gama  * pin / rin));

                if (vei > cin)
                {
                    if (vni >= 0.0)
                    {
                        for (int m = 0; m < nTurbulentEquation; ++ m)
                        {
                            qTurbulent(it1, jt1, kt1, m) = qTurbulent(is1, js1, ks1, m);
                            qTurbulent(it2, jt2, kt2, m) = qTurbulent(is2, js2, ks2, m);
                        }
                        viscousTurbulent(it1, jt1, kt1) = viscousTurbulent(is1, js1, ks1);
                        viscousTurbulent(it2, jt2, kt2) = viscousTurbulent(is2, js2, ks2);
                    }
                    else
                    {
                        for (int m = 0; m < nTurbulentEquation; ++ m)
                        {
                            qTurbulent(it1, jt1, kt1,m) = freeStreamTurbVar[m];
                            qTurbulent(it2, jt2, kt2,m) = freeStreamTurbVar[m];
                        }
                        viscousTurbulent(it1, jt1, kt1) = freeStreamViscosity;
                        viscousTurbulent(it2, jt2, kt2) = freeStreamViscosity;
                    }
                    continue;
                }

                RDouble riemp = vni + 2.0 * cin / (gama  - 1.0);
                RDouble riemm = vno - 2.0 * coo / (refGama - 1.0);
                RDouble vnb = half * (riemp + riemm);
                RDouble cb  = fourth * (riemp - riemm) * (gama - 1.0);

                if (vnb  >= 0.0)
                {
                    //! Exit.
                    for (int m = 0; m < nTurbulentEquation; ++ m)
                    {
                        qTurbulent(it1, jt1, kt1, m) = qTurbulent(is1, js1, ks1, m);
                        qTurbulent(it2, jt2, kt2, m) = qTurbulent(is2, js2, ks2, m);
                    }
                    viscousTurbulent(it1, jt1, kt1) = viscousTurbulent(is1, js1, ks1);
                    viscousTurbulent(it2, jt2, kt2) = viscousTurbulent(it1, jt1, kt1);
                }
                else
                {
                    //! Inlet.
                    for (int m = 0; m < nTurbulentEquation; ++ m)
                    {
                        qTurbulent(it1, jt1, kt1, m) = freeStreamTurbVar[m];
                        qTurbulent(it2, jt2, kt2, m) = freeStreamTurbVar[m];
                    }
                    viscousTurbulent(it1, jt1, kt1) = freeStreamViscosity;
                    viscousTurbulent(it2, jt2, kt2) = freeStreamViscosity;
                }
            }
        }
    }
    delete []prims1;
}

void TurbSolverStr::WallFunctionStandard(Grid *gridIn, BCRegion *BCRegion, int is, int js, int ks, int it, int jt, int kt, int i, int j, int k)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    RDouble4D &faceNormalX = *(grid->GetFaceNormalX());
    RDouble4D &faceNormalY = *(grid->GetFaceNormalY());
    RDouble4D &faceNormalZ = *(grid->GetFaceNormalZ());

    RDouble4D &qTurbulent       = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb"));
    RDouble4D &qLaminar         = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));
    RDouble3D &viscousLaminars  = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("visl"));
    RDouble3D &viscousTurbulent = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("vist"));
    RDouble3D &walldist = *StructGridCast(grid)->GetWallDist();

    Param_TurbSolverStruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    int iSurface = BCRegion->GetFaceDirection() + 1;

    using namespace IDX;

    RDouble rhoInside = qLaminar(is, js, ks, IDX::IR);
    RDouble rhoGhost = abs(qLaminar(it, jt, kt, IDX::IR));
    RDouble densityWall = 0.5*(rhoInside + rhoGhost);
    RDouble up = qLaminar(is, js, ks, IDX::IU);
    RDouble vp = qLaminar(is, js, ks, IDX::IV);
    RDouble wp = qLaminar(is, js, ks, IDX::IW);
    RDouble laminarViscosity   = viscousLaminars(is, js, ks);
    RDouble turbulentViscosity = viscousTurbulent(is, js, ks);
    
    RDouble xfn = faceNormalX(is, js, ks, iSurface);
    RDouble yfn = faceNormalY(is, js, ks, iSurface);
    RDouble zfn = faceNormalZ(is, js, ks, iSurface);

    RDouble vpn    = up * xfn + vp * yfn + wp * zfn;
    RDouble vpt    = sqrt(up * up + vp * vp + wp * wp - vpn * vpn);
    RDouble small  = 1.0e-15;
    RDouble uvw    = max(vpt, small);

    RDouble wallDistance   = walldist(is, js, ks);

    //! Obtain the reynolds number.
    RDouble refReNumber = parameters->GetRefReNumber();

    RDouble taow  = laminarViscosity * uvw / wallDistance / refReNumber;
    RDouble utao  = sqrt(taow / densityWall);

    RDouble E_     = 9.793;
    RDouble beta1  = 0.075;
    RDouble kappa1 = 0.4178;

    const RDouble yPluslimit = 11.225;

    RDouble yPlus;
    RDouble epslon = 1.0;
    int iter = 0;
    while (epslon > 0.001)
    {
        yPlus = refReNumber * wallDistance * densityWall * utao / laminarViscosity;
        yPlus = MAX(yPlus, yPluslimit);
        RDouble utemp = utao;

        utao = kappa1 * uvw / (log(E_ * yPlus));

        epslon = abs(utao-utemp)/(abs(utemp)+1.0e-20);
        if (++ iter > 20)
        {
            break;
        }
    }
    taow = densityWall * utao * utao;
    RDouble viscousTurbulenceWall = laminarViscosity * (taow / (uvw + small) / laminarViscosity * wallDistance * refReNumber - 1.0);
    viscousTurbulent(it, jt, kt) = 2.0 * viscousTurbulenceWall - viscousTurbulent(is, js, ks);

    if (nTurbulentEquation == 1)
    {
        qTurbulent(it, jt, kt, 0) = -qTurbulent(is, js, ks, 0);
    }
    else if (nTurbulentEquation == 2)
    {
        RDouble omgi = 60.0 * laminarViscosity / (beta1 * wallDistance * wallDistance * refReNumber * refReNumber);
        RDouble omgo = utao / (0.126 * wallDistance * refReNumber);
        //qTurbulent(is, js, ks, 1) = sqrt(omgi * omgi + omgo * omgo);
        qTurbulent(it, jt, kt, 1) = 2.0 * sqrt(omgi * omgi + omgo * omgo) - qTurbulent(is, js, ks, 1);
        qTurbulent(is, js, ks, 0) = qTurbulent(is, js, ks, 1) * viscousTurbulent(is, js, ks) / rhoInside;
        qTurbulent(it, jt, kt, 0) = qTurbulent(it, jt, kt, 1) * viscousTurbulent(it, jt, kt) / rhoGhost;
    }
}

void TurbSolverStr::WallFunctionPab3D(Grid *gridIn, BCRegion *BCRegion, int is, int js, int ks, int it, int jt, int kt, int i, int j, int k)
{
    RDouble paba10[7] = {2.354039, 0.1179840, -4.2899192e-04, 2.0404148e-06,-5.1775775e-09, 6.2687308e-12, -2.916958e-15};
    RDouble paba11[5] = {5.777191, 6.8756983e-02, -7.1582745e-06, 1.5594904e-09, -1.4865778e-13};
    RDouble paba12[5] = {31.08654, 5.0429072e-02, -2.0072314e-8};
    RDouble Ylim = 10.0;
    RDouble beta1 = 0.075;

    StructGrid *grid = StructGridCast(gridIn);

    RDouble4D &qTurbulent        = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb"));
    RDouble4D &qLaminar          = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q"));
    RDouble3D &viscousLaminar    = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("visl"));
    RDouble3D &viscousTurbulence = *reinterpret_cast<RDouble3D *> (grid->GetDataPtr("vist"));
    RDouble4D &temperature       = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("t"));
    RDouble3D &wallDistance      = *StructGridCast(grid)->GetWallDist();

    Param_TurbSolverStruct *parameters = GetControlParameters();

    int nm = GlobalDataBase::GetIntParaFromDB("nm");
    int nTurbulentEquation = parameters->GetNTurbulentEquation();

    int iSurface = BCRegion->GetFaceDirection() + 1;
    using namespace IDX;

    RDouble up = qLaminar(is, js, ks, IDX::IU);
    RDouble vp = qLaminar(is, js, ks, IDX::IV);
    RDouble wp = qLaminar(is, js, ks, IDX::IW);
    RDouble laminarViscosity  = viscousLaminar(is, js, ks);
    RDouble turbulentViscosity = viscousTurbulence(is, js, ks);

    RDouble temperatureWall = GlobalDataBase::GetDoubleParaFromDB("wallTemperature");
    RDouble refDimensionalTemperature = parameters->GetRefDimensionalTemperature();
    temperatureWall = temperatureWall / refDimensionalTemperature;

    RDouble PrandtlNumber = GlobalDataBase::GetDoubleParaFromDB("prl");
    RDouble rWF = pow(static_cast<RDouble>(PrandtlNumber), 1.0/3.0);

    RDouble refMachNumber = parameters->GetRefMachNumber();
    RDouble refGama = parameters->GetRefGama();

    RDouble uu = sqrt(up * up + vp * vp + wp * wp);

    RDouble wallTemperature  = 1.0;
    if (temperatureWall < 0.0)
    {
        RDouble moo2 = refMachNumber * refMachNumber;
        RDouble cp  = 1.0 / (refGama - 1.0) / moo2;
        RDouble tm  = temperature(is, js, ks, IDX::ITT);
        wallTemperature = tm + 0.5 * rWF * uu * uu / cp;
    }
    else
    {
        wallTemperature = temperatureWall;
    }

    RDouble pressure = qLaminar(is, js, ks, IDX::IP);
    //RDouble density = refGama * refMachNumber * refMachNumber * pressure / wallTemperature;
    RDouble temperatureSuthland = GlobalDataBase::GetDoubleParaFromDB("tsuth");
    RDouble rhoInside = qLaminar(is, js, ks, IDX::IR);
    RDouble rhoGhost = abs(qLaminar(it, jt, kt, IDX::IR));
    RDouble densityWall = 0.5*(rhoInside + rhoGhost);

    RDouble wallDist   = wallDistance(is, js, ks);

    RDouble refReNumber = parameters->GetRefReNumber(); //! Obtain the reynolds number.
    RDouble taow   = laminarViscosity * uu / wallDist / refReNumber;
    RDouble utao   = sqrt(taow / densityWall);
    RDouble yPlus  = refReNumber * wallDist * densityWall * utao / laminarViscosity;
    
    RDouble rc = densityWall * refReNumber * wallDist * uu / laminarViscosity;
    RDouble xnplus = 0.0;
    //! The original one here is yPlus.
    if (rc <= 20.24)
    {
        xnplus = sqrt(rc);
    }
    else if (rc <= 435.0)
    {
        xnplus = paba10[0] + paba10[1]*rc + paba10[2]*rc*rc + paba10[3]*rc*rc*rc + paba10[4]*rc*rc*rc*rc + paba10[5]*rc*rc*rc*rc*rc + paba10[6]*rc*rc*rc*rc*rc*rc;
    }
    else if (rc <= 4000.0)
    {
        xnplus = paba11[0] + paba11[1]*rc + paba11[2]*rc*rc + paba11[3]*rc*rc*rc + paba11[4]*rc*rc*rc*rc;
    }
    else if (rc > 4000.0)
    {
        xnplus = paba12[0] + paba12[1]*rc + paba12[2]*rc*rc;
    }
    RDouble xnplussav = xnplus;
    if (xnplus >= Ylim) 
    {
        //! Newton iteration to solve for nplus, assuming it is in log region:
        for (int iter = 0; ; ++ iter)
        {
            RDouble f = rc / xnplus - 2.44 * (log(xnplus)) - 5.2;
            RDouble dfdn = -rc / (xnplus * xnplus) - 2.44 / xnplus;
            RDouble delta = -f / dfdn;
            xnplus = fabs(xnplus + delta);
            if (iter > 20)
            {
                //! Revert back to approx series soln if Newton iteration fails.
                xnplus = xnplussav;
                break;
            }
            else if (abs(delta) < 1.e-3)
            {
                break;
            }
        }
    }

    RDouble xmut = laminarViscosity * (laminarViscosity * xnplus * xnplus / (wallDist * densityWall * uu * refReNumber) - 1.0);
    viscousTurbulence(it, jt, kt) = 2.0 * xmut - viscousTurbulence(is, js, ks);

    if (nTurbulentEquation == 1)
    {
        qTurbulent(it, jt, kt,0) = -qTurbulent(is, js, ks, 0);
    }
    else if (nTurbulentEquation == 2)
    {
        RDouble omgi = 60.0 * laminarViscosity / (beta1 * wallDist * wallDist * refReNumber * refReNumber);
        taow  = laminarViscosity * uu / wallDist / refReNumber;
        utao   = sqrt(taow / densityWall);
        RDouble omgo = utao / (0.126 * wallDist * refReNumber);
        //qTurbulent(is, js, ks, 1) = sqrt(omgi * omgi + omgo * omgo);
        qTurbulent(it, jt, kt, 1) = 2.0 * sqrt(omgi * omgi + omgo * omgo) - qTurbulent(is, js, ks, 1);
        qTurbulent(is, js, ks, 0) = qTurbulent(is, js, ks, 1) * viscousTurbulence(is, js, ks) / rhoInside;
        qTurbulent(it, jt, kt, 0) = qTurbulent(it, jt, kt, 1) * viscousTurbulence(it, jt, kt) / rhoGhost;
    }
}

void TurbSolverStr::UploadInterfaceData(ActionKey *actkey)
{
    UploadInterfaceValue(actkey);
}

void TurbSolverStr::UploadInterfaceValue(ActionKey *actkey)
{
    int level = actkey->level;
    StructGrid *grid = reinterpret_cast<StructGrid *> (GetGrid(level));
    InterfaceInfo *interfaceInfo = grid->GetInterfaceInfo();
    if (interfaceInfo == 0)
    {
        return;
    }

    Param_TurbSolverStruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();
    RDouble4D *qTurbulent = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb"));
    PHSPACE::UploadInterfaceValue(grid, qTurbulent, "turb::q", nTurbulentEquation);

    RDouble3D *viscousTurbulence = reinterpret_cast<RDouble3D *> (grid->GetDataPtr("vist"));
    PHSPACE::UploadInterfaceValue(grid, viscousTurbulence, "vist");

    //if (nTurbulentEquation == 2)
    //{
    //    RDouble3D* blend = reinterpret_cast<RDouble3D*> (grid->GetDataPtr("blend"));
    //    PHSPACE::UploadInterfaceValue(grid, blend, "blend");
    //}
}

void TurbSolverStr::DownloadInterfaceData(ActionKey *actkey)
{
    DownloadInterfaceValue(actkey);
}

void TurbSolverStr::DownloadInterfaceValue(ActionKey *actkey)
{
    int level = actkey->level;
    StructGrid *grid = reinterpret_cast<StructGrid *> (GetGrid(level));
    InterfaceInfo *interfaceInfo = grid->GetInterfaceInfo();
    if (interfaceInfo == 0)
    {
        return;
    }

    Param_TurbSolverStruct *parameters = GetControlParameters();
    int nTurbulentEquation = parameters->GetNTurbulentEquation();
    RDouble4D *qTurbulent = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb"));
    PHSPACE::DownloadInterfaceValue(grid, qTurbulent, "turb::q",  nTurbulentEquation);

    RDouble3D *viscousTurbulence = reinterpret_cast<RDouble3D *> (grid->GetDataPtr("vist"));
    PHSPACE::DownloadInterfaceValue(grid, viscousTurbulence, "vist");

    //if (nTurbulentEquation == 2)
    //{
    //    RDouble3D* blend = reinterpret_cast<RDouble3D*> (grid->GetDataPtr("blend"));
    //    PHSPACE::DownloadInterfaceValue(grid, blend, "blend");
    //}
}

//void TurbSolverStr::CompressDQ(DataContainer*& dataContainer, FieldProxy* fieldProxy, Grid* gridIn, const int& zoneIndex, const int& neighborZoneIndex, const int& nEquation)
//{
//    StructGrid* grid = StructGridCast(gridIn);
//    int level = grid->GetLevel();
//    InterfaceInfo* interfaceInformation = grid->GetInterfaceInfo();
//    if (!interfaceInformation)
//    {
//        return;
//    }
//
//    StructGrid* finestGrid = StructGridCast(grid->GetFinestGrid());
//    InterfaceInfo* finestInterfaceInfo = finestGrid->GetInterfaceInfo();
//
//    int iNeighborZone = finestInterfaceInfo->FindIthNeighbor(neighborZoneIndex);
//    int interfaceNumberBetweenTwoNeighboringZone = finestInterfaceInfo->GetNIFaceOfNeighbor(iNeighborZone);
//    int* interfaceIndexContainerForSend = finestInterfaceInfo->GetFaceIndexForSend(iNeighborZone);
//
//    dataContainer->moveto_begin();
//
//    RDouble4D& fieldSend = fieldProxy->GetField_STR();
//
//    for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; --iGhostLayer)
//    {
//        for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++iLocalFace)
//        {
//            int iFace = interfaceIndexContainerForSend[iLocalFace];
//            int is, js, ks;
//            finestGrid->GetSourceIndexIJK(iFace, iGhostLayer + 1, is, js, ks);
//            finestGrid->RemapMultigridIJK(level, is, js, ks);
//            for (int m = 0; m < nEquation; ++m)
//            {
//                PHWrite(dataContainer, fieldSend(is, js, ks, m));
//            }
//        }
//    }
//}
//
//void TurbSolverStr::DecompressDQ(DataContainer*& dataContainer, FieldProxy* fieldProxy, Grid* gridIn, const int& zoneIndex, const int& neighborZoneIndex, const int& nEquation)
//{
//    StructGrid* grid = StructGridCast(gridIn);
//    int level = grid->GetLevel();
//    InterfaceInfo* interfaceInformation = grid->GetInterfaceInfo();
//    if (!interfaceInformation)
//    {
//        return;
//    }
//
//    StructGrid* finestGrid = StructGridCast(grid->GetFinestGrid());
//    InterfaceInfo* finestInterfaceInfo = finestGrid->GetInterfaceInfo();
//
//    int iNeighborZone = finestInterfaceInfo->FindIthNeighbor(neighborZoneIndex);
//    int interfaceNumberBetweenTwoNeighboringZone = finestInterfaceInfo->GetNIFaceOfNeighbor(iNeighborZone);
//    int* interfaceIndexContainerForReceive = finestInterfaceInfo->GetFaceIndexForRecv(iNeighborZone);
//
//    dataContainer->moveto_begin();
//
//    RDouble4D& fieldRecv = fieldProxy->GetField_STR();
//
//    for (int iGhostLayer = GetNumberOfGhostCellLayers() - 1; iGhostLayer >= 0; --iGhostLayer)
//    {
//        for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++iLocalFace)
//        {
//            int iFace = interfaceIndexContainerForReceive[iLocalFace];
//            int it, jt, kt;
//            finestGrid->GetTargetIndexIJK(iFace, iGhostLayer + 1, it, jt, kt);
//            finestGrid->RemapMultigridIJK(level, it, jt, kt);
//            for (int m = 0; m < nEquation; ++m)
//            {
//                PHRead(dataContainer, fieldRecv(it, jt, kt, m));
//            }
//        }
//    }
//}

void TurbSolverStr::ComputeLengthScaleofOneEquationModel(Grid *gridIn)
{
    Param_TurbSolverStruct *parameters = GetControlParameters();
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

void TurbSolverStr::ComputeDESLength(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, 0, -1, I, J, K);

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    RDouble3D *wallDistant = grid->GetWallDist();

    RDouble3D *DESLength = this->GetDESLength();

    Param_TurbSolverStruct *parameters = GetControlParameters();
    int DESType = parameters->GetDESType();

    if (!DESLength)
    {
        DESLength = new RDouble3D (I, J, K, fortranArray);
    }

    RDouble3D *LESLength = new RDouble3D (I, J, K, fortranArray);
    RDouble3D *lowReynoldsNumberCorrection = new RDouble3D (I, J, K, fortranArray);
    ComputeLESLengthofSA(gridIn, LESLength, lowReynoldsNumberCorrection, DESType);

    RDouble lengthRANS, lengthLES;
    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                lengthRANS = (*wallDistant)(i, j, k);
                lengthLES  = (*LESLength  )(i, j, k);
                (*DESLength)(i, j, k) = MIN(lengthRANS, lengthLES);
            }
        }
    }

    this->SetDESLength(DESLength);

    delete LESLength;
    delete lowReynoldsNumberCorrection;
}

void TurbSolverStr::ComputeDDESLength(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, 0, -1, I, J, K);

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    RDouble3D *wallDistant = grid->GetWallDist();

    RDouble3D *DDESLength = this->GetDESLength();

    Param_TurbSolverStruct *parameters = GetControlParameters();
    int DESType = parameters->GetDESType();

    if (!DDESLength)
    {
        DDESLength = new RDouble3D (I, J, K, fortranArray);
    }

    RDouble3D *LESLength = new RDouble3D (I, J, K, fortranArray);
    RDouble3D *lowReynoldsNumberCorrection = new RDouble3D(I, J, K, fortranArray);
    ComputeLESLengthofSA(gridIn, LESLength, lowReynoldsNumberCorrection, DESType);

    RDouble4D *qLaminar         = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q" ));
    RDouble3D *viscousLaminar   = reinterpret_cast<RDouble3D *> (grid->GetDataPtr("visl"));
    RDouble3D *viscousTurbulent = reinterpret_cast<RDouble3D *> (grid->GetDataPtr("vist"));

    Range II, JJ, KK;
    GetRange(ni, nj, nk, -2, 1, II, JJ, KK);

    RDouble4D &gradUVWTCellCenterX = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterX"));
    RDouble4D &gradUVWTCellCenterY = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterY"));
    RDouble4D &gradUVWTCellCenterZ = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterZ"));

    using namespace IDX;

    RDouble oRefReNumber = parameters->GetoRefReNumber();

    const RDouble karm   = 0.41;
    const RDouble kap2   = karm * karm;

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble laminarViscosity   = (*viscousLaminar)(i, j, k);
                RDouble turbulentViscosity = (*viscousTurbulent)(i, j, k);
                RDouble density            = (*qLaminar)(i, j, k, IR);

                RDouble wallDistance  = (*wallDistant)(i, j, k);
                RDouble wallDistance2 = SQR(wallDistance);
                RDouble lengthRans    = wallDistance;
                RDouble lengthLes     = (*LESLength)(i, j, k);

                RDouble dudx = gradUVWTCellCenterX(i, j, k, 0);
                RDouble dudy = gradUVWTCellCenterY(i, j, k, 0);
                RDouble dudz = gradUVWTCellCenterZ(i, j, k, 0);

                RDouble dvdx = gradUVWTCellCenterX(i, j, k, 1);
                RDouble dvdy = gradUVWTCellCenterY(i, j, k, 1);
                RDouble dvdz = gradUVWTCellCenterZ(i, j, k, 1);

                RDouble dwdx = gradUVWTCellCenterX(i, j, k, 2);
                RDouble dwdy = gradUVWTCellCenterY(i, j, k, 2);
                RDouble dwdz = gradUVWTCellCenterZ(i, j, k, 2);

                RDouble sumOfgradient2 = SQR(dudx, dudy, dudz) + 
                                         SQR(dvdx, dvdy, dvdz) + 
                                         SQR(dwdx, dwdy, dwdz);

                RDouble nueff = oRefReNumber * (laminarViscosity + turbulentViscosity) / density;
                RDouble rd    = nueff / (kap2 * wallDistance2 * MAX(sqrt(sumOfgradient2), 1.0e-20));
                RDouble fd    = 1.0 - tanh(static_cast<RDouble>(POWER3(8.0 * rd)));

                (*DDESLength)(i, j, k) = lengthRans - fd * MAX(0.0, (lengthRans - lengthLes));

            }
        }
    }

    this->SetDESLength(DDESLength);

    delete LESLength;
    delete lowReynoldsNumberCorrection;
}

void TurbSolverStr::ComputeIDDESLength(Grid *gridIn)
{
    StructGrid * grid = StructGridCast(gridIn);
    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, 0, -1, I, J, K);

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    RDouble3D *wallDistant = grid->GetWallDist();

    RDouble3D *largestLocalGridLength = grid->GetLargestLocalGridLength();

    RDouble3D *IDDESLength = this->GetDESLength();

    Param_TurbSolverStruct *parameters = GetControlParameters();
    int DESType = parameters->GetDESType();

    if (!IDDESLength)
    {
        IDDESLength = new RDouble3D (I, J, K, fortranArray);
    }

    RDouble3D *LESLength = new RDouble3D (I, J, K, fortranArray);
    RDouble3D *lowReynoldsNumberCorrection = new RDouble3D (I, J, K, fortranArray);
    ComputeLESLengthofSA(gridIn, LESLength, lowReynoldsNumberCorrection, DESType);

    RDouble4D *qLaminar         = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q" ));
    RDouble3D *viscousLaminar   = reinterpret_cast<RDouble3D *> (grid->GetDataPtr("visl"));
    RDouble3D *viscousTurbulent = reinterpret_cast<RDouble3D *> (grid->GetDataPtr("vist"));

    Range II, JJ, KK;
    GetRange(ni, nj, nk, -2, 1, II, JJ, KK);

    RDouble4D &gradUVWTCellCenterX = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterX"));
    RDouble4D &gradUVWTCellCenterY = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterY"));
    RDouble4D &gradUVWTCellCenterZ = *reinterpret_cast<RDouble4D *> (grid->GetDataPtr("gradUVWTCellCenterZ"));

    using namespace IDX;

    RDouble oRefReNumber = parameters->GetoRefReNumber();

    const RDouble karm   = 0.41;
    const RDouble kap2   = karm * karm;

    RDouble ct = 1.63;
    RDouble cl = 3.55;

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble laminarViscosity   = (*viscousLaminar)(i, j, k);
                RDouble turbulentViscosity = (*viscousTurbulent)(i, j, k);
                RDouble density            = (*qLaminar)(i, j, k,IR);
                RDouble wallDistance       = (*wallDistant)(i, j, k);
                RDouble wallDistance2      = SQR(wallDistance);
                RDouble lengthRans         = wallDistance;
                RDouble lengthLes          = (*LESLength)(i, j, k);
                RDouble largestLocalGridSpacing = (*largestLocalGridLength)(i, j, k);

                RDouble alf = 0.25 - wallDistance / largestLocalGridSpacing;
                RDouble fb  = MIN( 2.0 * exp(- 9.0 * alf * alf),  1.0 );

                RDouble dudx = gradUVWTCellCenterX(i, j, k, 0);
                RDouble dudy = gradUVWTCellCenterY(i, j, k, 0);
                RDouble dudz = gradUVWTCellCenterZ(i, j, k, 0);

                RDouble dvdx = gradUVWTCellCenterX(i, j, k, 1);
                RDouble dvdy = gradUVWTCellCenterY(i, j, k, 1);
                RDouble dvdz = gradUVWTCellCenterZ(i, j, k, 1);

                RDouble dwdx = gradUVWTCellCenterX(i, j, k, 2);
                RDouble dwdy = gradUVWTCellCenterY(i, j, k, 2);
                RDouble dwdz = gradUVWTCellCenterZ(i, j, k, 2);

                RDouble sumOfgradient2 = SQR(dudx, dudy, dudz) + 
                                         SQR(dvdx, dvdy, dvdz) + 
                                         SQR(dwdx, dwdy, dwdz);

                RDouble nuet = oRefReNumber * turbulentViscosity / density;
                RDouble rdt  = nuet / (kap2 * wallDistance2 * MAX(sqrt(sumOfgradient2), 1.0e-10));
                RDouble ft   = tanh(static_cast<RDouble>(POWER3(ct * ct * rdt)));

                RDouble fdt  = 1.0 - tanh(static_cast<RDouble>(POWER3(8.0 * rdt)));
                RDouble fdb  = MAX((1.0 - fdt), fb);

                RDouble nuel = oRefReNumber * laminarViscosity / density;
                RDouble rdl  = nuel / (kap2 * wallDistance2 * MAX(sqrt(sumOfgradient2), 1.0e-10));
                RDouble fl   = tanh(static_cast<RDouble>(pow(static_cast<RDouble>(cl * cl * rdl), 10.0)));

                RDouble fe2  = 1.0 - MAX(ft, fl);

                RDouble fe1  = 0.0;
                if (alf < 0.0)
                {
                    fe1 = 2.0 * exp(- 9.0 * alf * alf);
                }
                else
                {
                    fe1 = 2.0 * exp(- 11.09 * alf * alf);
                }

                RDouble fe                   = fe2 * MAX((fe1 - 1.0), 0.0) * (*lowReynoldsNumberCorrection)(i, j, k);

                (*IDDESLength)(i, j, k)       = fdb * (1.0 + fe) * lengthRans  +  (1.0 - fdb) * lengthLes;

            }
        }
    }

    this->SetDESLength(IDDESLength);

    delete LESLength;
    delete lowReynoldsNumberCorrection;
}

RDouble3D *TurbSolverStr::GetLengthScale(Grid *gridIn)
{
    StructGrid *grid = StructGridCast(gridIn);

    Param_TurbSolverStruct *parameters = GetControlParameters();
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

void ComputeLESLengthofSA(Grid *gridIn, RDouble3D *LESLength, RDouble3D *lowReynoldsNumberCorrection, int DESType)
{
    StructGrid *grid = StructGridCast(gridIn);

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    //! Compute LES length.
    RDouble CDES = 0.65;

    ComputeLowReynoldsNumberCorrection(gridIn, lowReynoldsNumberCorrection);

    if (DESType == IDDES)
    {
        RDouble3D *subGridLength = grid->GetSubgridLength();
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    (*LESLength)(i, j, k) = CDES * (*lowReynoldsNumberCorrection)(i, j, k) * (*subGridLength)(i, j, k);
                }
            }
        }
    }
    else
    {
        RDouble3D *largestLocalGridLength = grid->GetLargestLocalGridLength();
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    (*LESLength)(i, j, k) = CDES * (*lowReynoldsNumberCorrection)(i, j, k) * (*largestLocalGridLength)(i, j, k);
                }
            }
        }
    }
}

void ComputeLowReynoldsNumberCorrection(Grid *gridIn, RDouble3D *lowReynoldsNumberCorrection)
{
    StructGrid * grid = StructGridCast(gridIn);

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    //! Compute low Reynolds number correction.
    RDouble4D *qLaminar       = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q" ));
    RDouble3D *viscousLaminar = reinterpret_cast<RDouble3D *> (grid->GetDataPtr("visl"));
    RDouble4D *qTurbulent     = reinterpret_cast<RDouble4D *> (grid->GetDataPtr("q_turb"));

    const RDouble ct3   = 1.2;
    const RDouble ct4   = 0.5;
    const RDouble karm  = 0.41;
    const RDouble cb1   = 0.1355;
    const RDouble cb2   = 0.622;
    const RDouble cv1   = 7.1;
    const RDouble cw2   = 0.3;
    const RDouble cw3   = 2.0;
    const RDouble sigma = 2.0/3.0;
    const RDouble cdes  = - 1.0;
    const RDouble cv13  = POWER3(cv1 * cv1 * cv1);
    const RDouble kap2  = karm * karm;
    const RDouble rkap2 = one / kap2;
    const RDouble cw1   = cb1 * rkap2 + (one + cb2) / sigma;
    const RDouble turbulentFwStar = 0.424;

    using namespace IDX;
    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble nuet = (*qTurbulent)(i, j, k, ISA);
                RDouble olam = (*qLaminar )(i, j, k, IR) / ((* viscousLaminar)(i, j, k) + SMALL);
                RDouble xsi  = nuet * olam + SMALL;
                RDouble xsi3 = POWER3(xsi);
                RDouble xsi2 = SQR(xsi);
        
                RDouble turbulentFt2 = 0.0; 
                int codeOfFt2 = 0;
                if (codeOfFt2)
                {
                    turbulentFt2 = ct3 * exp(- ct4 * xsi2);
                }
        
                RDouble turbulentFv1 = xsi3 / (xsi3 + cv13);
                RDouble turbulentFv2 = one - xsi / (one + xsi * turbulentFv1);
        
                RDouble numerator   =  1.0  -  (cb1 / (cw1 * kap2 * turbulentFwStar)) * (turbulentFt2 + (1.0 - turbulentFt2) * turbulentFv2);
                RDouble denominator =  turbulentFv1 * MAX(1.0e-10, (1.0 - turbulentFt2));
        
                (*lowReynoldsNumberCorrection)(i, j, k) = sqrt( MIN(1.0e2, numerator / denominator) );
            }
        }
    }
}

}