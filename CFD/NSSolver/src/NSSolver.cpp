#include "NSSolver.h"
#include "Constants.h"
#include "Glb_Dimension.h"
#include "PHMpi.h"
#include "Param_NSSolver.h"
#include "Math_BasisFunction.h"
#include "PHIO.h"

#ifdef USE_HYFLOW_CODE
#include "RealGas.h"
#else
#include "Gas.h"
#endif
#include "TK_Log.h"

using namespace std;

namespace PHSPACE
{

using namespace GAS_SPACE;
NSSolver::NSSolver()
{
}

NSSolver::~NSSolver()
{
}

LIB_EXPORT Param_NSSolver *NSSolver::GetControlParameters() const
{
    return static_cast <Param_NSSolver *> (controlParameters);
}

int NSSolver::GetNumberOfEquations()
{
    Param_NSSolver *parameters = GetControlParameters();
    int nLaminar = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();

    int nEquation = nLaminar;
    if (nChemical > 0)
    {
        nEquation += nChemical + nTemperatureModel - 1;
    }

    return nEquation;
}

void NSSolver::TimeStep(Grid *grid)
{
    Param_NSSolver *parameters = GetControlParameters();

    int ifLocalTimeStep = parameters->GetIfLocalTimeStep();

    if (ifLocalTimeStep == 0)
    {
        LocalTimeStep(grid);
    }
    else if (ifLocalTimeStep == 1)
    {
        GlobalTimeStep(grid);
    }
    else
    {
        LocalGlobalTimeStep(grid);
    }
}

void NSSolver::UpdateResiduals(Grid *grid)
{
    InviscidFlux(grid);

    ViscousFlux(grid);

    SourceFlux(grid);

    ZeroResidualOfSpecialCells(grid);

    FreeGradientProxy(grid);    //! It is null for structured grid.
}

//! Get the file name for Residual dumping.
const string NSSolver::GetResidualFileName()
{
    string resSaveFile = "res.dat";
    GlobalDataBase::GetData("resSaveFile", &resSaveFile, PHSTRING, 1);
    return resSaveFile;
}

void NSSolver::Post()
{
    using namespace PHMPI;

    //! To dump restart data for continual simulation.
    ActionKey *actkeyDumpRestartData = new ActionKey();
    FillActionKey(actkeyDumpRestartData, DUMP_RESTART, 0);
    DumpRestartData(actkeyDumpRestartData);
    delete actkeyDumpRestartData;

    //! To dump residual data.
    ActionKey *actkeyDumpResidual = new ActionKey();
    FillActionKey(actkeyDumpResidual, DUMP_RESIDUAL, 0);
    DumpResidual(actkeyDumpResidual);
    delete actkeyDumpResidual;

    //! To dump airforce coefficient data.
    ActionKey *actkeyDumpAirForceCoef = new ActionKey();
    FillActionKey(actkeyDumpAirForceCoef, DUMP_AIR_FORCE_COEF, 0);
    DumpAirForceCoef(actkeyDumpAirForceCoef);
    delete actkeyDumpAirForceCoef;

    //! To dump wall airforce coefficient data.
    ActionKey *actkeyDumpWallAircoef = new ActionKey();
    FillActionKey(actkeyDumpWallAircoef, DUMP_CP_DISTRI, 0);
    DumpCpDistri(actkeyDumpWallAircoef);
    delete actkeyDumpWallAircoef;

    if (IsNeedStatistics())
    {
        //! To output Average data by tecplot format.
        ActionKey *actkeyTecOutAverageFlow = new ActionKey();
        FillActionKey(actkeyTecOutAverageFlow, VISUALIZATION_AVERAGE_FLOW, 0);
        TecOutAverageFlow(actkeyTecOutAverageFlow);
        delete actkeyTecOutAverageFlow;
    }

    if (IsNeedReynoldsStressStatistics())
    {
        //! To output Average Reynolds stress by tecplot format.
        ActionKey *actkeyTecOutAverageReynoldsStress = new ActionKey();
        FillActionKey(actkeyTecOutAverageReynoldsStress, VISUALIZATION_AVERAGE_ReynoldsStress, 0);
        TecOutAverageReynoldsStress(actkeyTecOutAverageReynoldsStress);
        delete actkeyTecOutAverageReynoldsStress;
    }
}

void NSSolver::PostSolve(Grid *grid, int stage, int level)
{
    CommunicationInterfaceData();
}

void NSSolver::DumpResultFile(Grid *grid, int level)
{
    Param_NSSolver *parameters = GetControlParameters();
    int intervalStepFlow = parameters->GetIntervalStepFlow();

    int intervalStepPlot = parameters->GetIntervalStepPlot();
    int outnstep = GlobalDataBase::GetIntParaFromDB("outnstep");

    int sys_gridtype = 0;
    GlobalDataBase::GetData("sys_gridtype", &sys_gridtype, PHINT, 1);

    int nProtectData = 0;
    if(GlobalDataBase::IsExist("nProtectData", PHINT, 1))
    {
        nProtectData = GlobalDataBase::GetIntParaFromDB("nProtectData");
    }

    bool isSubIterationDump = true;
    int isUnsteady = parameters->GetIsUnsteady();
    if (isUnsteady)
    {
        GlobalDataBase::GetData("isSubIterationDump", &isSubIterationDump, PHBOOL, 1);
    }

    if (outnstep % intervalStepFlow == 0 && isSubIterationDump)
    {
        int dumpStandardModel = 0;
        GlobalDataBase::GetIntParaFromDB("dumpStandardModel");
        if(dumpStandardModel == 1)
        {
            Turbulence_Flat_Plate_Output_NASA(grid);
        }
    }

    if (outnstep % intervalStepFlow == 0 && level == 0 && isSubIterationDump)
    {
        if(nProtectData ==0)
        {
            //! To dump restart data for continual simulation.
            ActionKey *actkeyDumpRestartData = new ActionKey();
            FillActionKey(actkeyDumpRestartData, DUMP_RESTART, level);
            DumpRestartData(actkeyDumpRestartData);
            delete actkeyDumpRestartData;
        }

        else
        {
            int RestartFile = 0;
            ActionKey *actkeyDumpRestartData = new ActionKey();
            FillActionKey(actkeyDumpRestartData, DUMP_RESTART, level);
            DumpProtectData(actkeyDumpRestartData, RestartFile);
            delete actkeyDumpRestartData;

            if(outnstep > 1)
            {
                if(RestartFile == 0)
                {
                    string path = "./results/flow1.dat";
                    remove(path.c_str());

                    if(PHMPI::IsParallelRun())
                    {
                        string path = "./results/flow1_0.dat";
                        remove(path.c_str());
                    }
                }

                else if(RestartFile == 1)
                {
                    string path = "./results/flow0.dat";
                    remove(path.c_str());

                    if(PHMPI::IsParallelRun())
                    {
                        string path = "./results/flow0_0.dat";
                        remove(path.c_str());
                    }
                }
            }
        }
    }

    if (outnstep % intervalStepPlot == 0 && isSubIterationDump)
    {
        if (IsNeedStatistics())
        {
            //! To output Average data by tecplot format.
            ActionKey *actkeyTecOutAverageFlow = new ActionKey();
            FillActionKey(actkeyTecOutAverageFlow, VISUALIZATION_AVERAGE_FLOW, level);
            TecOutAverageFlow(actkeyTecOutAverageFlow);
            delete actkeyTecOutAverageFlow;
        }

        if (IsNeedReynoldsStressStatistics())
        {
            //! To output Average Reynolds stress by tecplot format.
            ActionKey *actkeyTecOutAverageReynoldsStress = new ActionKey();
            FillActionKey(actkeyTecOutAverageReynoldsStress, VISUALIZATION_AVERAGE_ReynoldsStress, level);
            TecOutAverageReynoldsStress(actkeyTecOutAverageReynoldsStress);
            delete actkeyTecOutAverageReynoldsStress;
        }
    }

    int intervalStepRes = parameters->GetIntervalStepRes();
    if ((outnstep == 1 || outnstep % intervalStepRes == 0) && (!isUnsteady || !isSubIterationDump))
    {
        //! To dump residual data.
        ActionKey *actkeyDumpResidual = new ActionKey();
        FillActionKey(actkeyDumpResidual, DUMP_RESIDUAL, level);
        DumpResidual(actkeyDumpResidual);
        delete actkeyDumpResidual;
    }

    int intervalStepForce = parameters->GetIntervalStepForce();
    if (intervalStepForce > 0)
    {
        if (outnstep % intervalStepForce == 0 && isSubIterationDump)
        {
            if (level == 0)
            {
                //! To dump airforce coefficient data.
                ActionKey *actkeyDumpAirForceCoef = new ActionKey();
                FillActionKey(actkeyDumpAirForceCoef, DUMP_AIR_FORCE_COEF, level);
                DumpAirForceCoef(actkeyDumpAirForceCoef);
                delete actkeyDumpAirForceCoef;

                //! To dump wall airforce coefficient data.
                ActionKey *actkeyDumpWallAircoef = new ActionKey();
                FillActionKey(actkeyDumpWallAircoef, DUMP_CP_DISTRI, 0);
                DumpCpDistri(actkeyDumpWallAircoef);
                delete actkeyDumpWallAircoef;

                //! To dump surface information data.
                ActionKey *actkeyDumpSurfaceInfo = new ActionKey();
                FillActionKey(actkeyDumpSurfaceInfo, DUMP_SURFACE_INFO, 0);
                DumpSurfaceInfo(actkeyDumpSurfaceInfo);
                delete actkeyDumpSurfaceInfo;

                if (sys_gridtype == STRUCTGRID)
                {
                    //! To dump surface max and average heat flux.
                    ActionKey *actkeyDumpHeatFlux = new ActionKey();
                    FillActionKey(actkeyDumpHeatFlux, DUMP_HEATFLUX, 0);
                    DumpHeatFlux(actkeyDumpHeatFlux);
                    delete actkeyDumpHeatFlux;
                }
            }
        }
    }
}

void NSSolver::CheckResult(Grid *grid, int level)
{
    //! Compute the residual variation.
    ActionKey *actkeyCheckResidual = new ActionKey();
    FillActionKey(actkeyCheckResidual, DUMP_RESIDUAL, level);
    CheckResidual(actkeyCheckResidual);
    delete actkeyCheckResidual;

    Param_NSSolver *parameters = GetControlParameters();
    int nChemical = parameters->GetChemicalFlag();
    int nTModel = parameters->GetTemperatureModel();
    int iBranchCode = parameters->GetSelfAdaptionSolveFlag();
    int interResStep = parameters->GetIntervalStepRes();
    int outnstep = GlobalDataBase::GetIntParaFromDB("outnstep");
    int nCurrentStep = GlobalDataBase::GetIntParaFromDB("nCurrentStep");
    int nMaxSimuStep = GlobalDataBase::GetIntParaFromDB("maxSimuStep");
    RDouble monitorResidual = GlobalDataBase::GetDoubleParaFromDB("monitorResidual");

    int nPauseKey = 0;
    bool isConvergent = false;
    if (iBranchCode == 1)
    {
        int nIterFirstStep = 1000, nIterSecondStep = 2000, nIterThirdStep = 2000;
        if (GlobalDataBase::IsExist("nIterFirstStep", PHINT, 1))
        {
            nIterFirstStep = GlobalDataBase::GetIntParaFromDB("nIterFirstStep");
        }
        if (GlobalDataBase::IsExist("nIterSecondStep", PHINT, 1))
        {
            nIterSecondStep = GlobalDataBase::GetIntParaFromDB("nIterSecondStep");
        }
        if (GlobalDataBase::IsExist("nIterThirdStep", PHINT, 1))
        {
            nIterThirdStep = GlobalDataBase::GetIntParaFromDB("nIterThirdStep");
        }

        int nIteration = nMaxSimuStep;
        if (nCurrentStep == 1)
        {
            nIteration = nIterFirstStep;
        }
        else if (nCurrentStep == 2)
        {
            nIteration = nIterFirstStep + nIterSecondStep;
        }
        else if (nCurrentStep == 3)
        {
            nIteration = nIterFirstStep + nIterSecondStep + nIterThirdStep;
        }

        if (outnstep >= nIteration)
        {
            isConvergent = true;
        }
    }
    else if (iBranchCode == 2)
    {
        RDouble firstStepError = 0.01, secondStepError = 0.001, thirdStepError = 0.001;
        if (GlobalDataBase::IsExist("firstStepError", PHDOUBLE, 1))
        {
            firstStepError = GlobalDataBase::GetDoubleParaFromDB("firstStepError");
        }
        if (GlobalDataBase::IsExist("secondStepError", PHDOUBLE, 1))
        {
            secondStepError = GlobalDataBase::GetDoubleParaFromDB("secondStepError");
        }
        if (GlobalDataBase::IsExist("thirdStepError", PHDOUBLE, 1))
        {
            thirdStepError = GlobalDataBase::GetDoubleParaFromDB("thirdStepError");
        }

        if (nCurrentStep == 1)
        {
            if (monitorResidual <= firstStepError)
            {
                isConvergent = true;
            }
        }
        else if (nCurrentStep == 2)
        {
            if (monitorResidual <= secondStepError)
            {
                isConvergent = true;
            }
        }
        else if (nCurrentStep == 3)
        {
            if (monitorResidual <= thirdStepError)
            {
                isConvergent = true;
            }
        }
    }

    if (isConvergent)
    {
        nPauseKey = 1;
        PHSPACE::GlobalDataBase::UpdateData("nPauseKey", &nPauseKey, PHINT, 1);
    }

    int nCurrentProcessorID = PHMPI::GetCurrentProcessorID();
    if (nCurrentProcessorID == 0 && (outnstep % interResStep == 0))
    {
        string keyVariable = GlobalDataBase::GetStrParaFromDB("keyVariable");
        ostringstream oss;
        oss<<"The residual variation of "<<keyVariable<<" is: ";
        oss << setiosflags(ios::left);
        oss << setprecision(5);
        oss << setiosflags(ios::scientific);
        oss << setiosflags(ios::showpoint);
        oss<<monitorResidual<<endl;
        cout<<oss.str();
        //cout<<"The residual variation of "<<keyVariable<<" is: "<<monitorResidual<<endl;
    }
}

void NSSolver::DumpCFLNumber(Grid *grid, int level)
{
    //! Obtain the CFL number.
    ActionKey *actkeyObtainCFLNumber = new ActionKey();
    FillActionKey(actkeyObtainCFLNumber, DUMP_RESIDUAL, level);
    //! Obtain the minimum CFL number and maximum CFL number.
    ObtainCFLNumber(actkeyObtainCFLNumber);

    Param_NSSolver *parameters = GetControlParameters();
    int interResStep = parameters->GetIntervalStepRes();
    int outnstep = GlobalDataBase::GetIntParaFromDB("outnstep");

    int nCurrentProcessorID = PHMPI::GetCurrentProcessorID();
    if (nCurrentProcessorID == 0 && (outnstep % interResStep == 0))
    {
        RDouble globalMinCFL = GlobalDataBase::GetDoubleParaFromDB("globalMinCFL");
        RDouble globalMaxCFL = GlobalDataBase::GetDoubleParaFromDB("globalMaxCFL");

        //! Dump to file.
        string saveHeatRes = "./results/dumpCFL.plt";
        if (GlobalDataBase::IsExist("saveHeatRes", PHSTRING, 1))
        {
            saveHeatRes = GlobalDataBase::GetStrParaFromDB("saveHeatRes");
        }
        actkeyObtainCFLNumber->filename = saveHeatRes;
        if (actkeyObtainCFLNumber->filename == "")
        {
            return;
        }
        ios_base::openmode openmode = ios_base::out|ios_base::app;
        actkeyObtainCFLNumber->openmode = openmode;
        //Open file.
        ParallelOpenFile(actkeyObtainCFLNumber);
        fstream &file = *actkeyObtainCFLNumber->file;
        std::ostringstream oss;
        if (IfFileEmpty(file))
        {
            vector<string> title_tecplot;
            title_tecplot.push_back("Title=\"The change of CFL number\"");
            title_tecplot.push_back("Variables=");
            title_tecplot.push_back("\"Iterations\"");
            title_tecplot.push_back("\"globalMinCFL\"");
            title_tecplot.push_back("\"globalMaxCFL\"");
            for (std::size_t i = 0; i < title_tecplot.size(); ++ i)
            {
                oss << title_tecplot[i] << "\n";
            }
        }
        oss << setiosflags(ios::left);
        oss << setprecision(5);
        //oss << setiosflags(ios::scientific);
        oss << setiosflags(ios::showpoint);

        oss << setw(7) << outnstep << "    ";
        oss << globalMinCFL << "    "<<globalMaxCFL<<endl;

        WriteASCIIFile(file, oss.str());

        //Close file.
        ParallelCloseFile(actkeyObtainCFLNumber);
    }

    delete actkeyObtainCFLNumber;
}

void NSSolver::CheckSurfaceHeatingChange(Grid *grid, int level)
{
    //! Compute the surface heating change.
    ActionKey *actkeyCheckHeatChange = new ActionKey();
    FillActionKey(actkeyCheckHeatChange, DUMP_RESIDUAL, level);
    ComputeSurfaceHeatingChange(actkeyCheckHeatChange);

    Param_NSSolver *parameters = GetControlParameters();
    int interResStep = parameters->GetIntervalStepRes();
    int outnstep = GlobalDataBase::GetIntParaFromDB("outnstep");

    int nCurrentProcessorID = PHMPI::GetCurrentProcessorID();
    if (nCurrentProcessorID == 0 && (outnstep % interResStep == 0))
    {
        RDouble monitorHeatingChange = GlobalDataBase::GetDoubleParaFromDB("monitorHeatingChange");
        ostringstream oss1;
        oss1 << setiosflags(ios::left);
        oss1 << setprecision(5);
        oss1 << setiosflags(ios::scientific);
        oss1 << setiosflags(ios::showpoint);
        oss1 << "The maximum change of surface heating ratio is "<<monitorHeatingChange<<endl;
        cout<<oss1.str();
        //cout<<"The maximum change of surface heating ratio is "<<monitorHeatingChange<<endl;

        //! Dump to file.
        string saveHeatRes = "./results/heatResidual.plt";
        if (GlobalDataBase::IsExist("saveHeatRes", PHSTRING, 1))
        {
            saveHeatRes = GlobalDataBase::GetStrParaFromDB("saveHeatRes");
        }
        actkeyCheckHeatChange->filename = saveHeatRes;
        if (actkeyCheckHeatChange->filename == "")
        {
            return;
        }
        ios_base::openmode openmode = ios_base::out|ios_base::app;
        actkeyCheckHeatChange->openmode = openmode;
        //Open file.
        ParallelOpenFile(actkeyCheckHeatChange);
        fstream &file = *actkeyCheckHeatChange->file;
        std::ostringstream oss;
        if (IfFileEmpty(file))
        {
            vector<string> title_tecplot;
            title_tecplot.push_back("Title=\"The Surface Heating Change\"");
            title_tecplot.push_back("Variables=");
            title_tecplot.push_back("\"Iterations\"");
            title_tecplot.push_back("\"Heating Change\"");
            for (std::size_t i = 0; i < title_tecplot.size(); ++ i)
            {
                oss << title_tecplot[i] << "\n";
            }
        }
        oss << setiosflags(ios::left);
        oss << setprecision(5);
        oss << setiosflags(ios::scientific);
        oss << setiosflags(ios::showpoint);

        oss << setw(7) << outnstep << "    ";
        oss << monitorHeatingChange <<endl;

        WriteASCIIFile(file, oss.str());

        //Close file.
        ParallelCloseFile(actkeyCheckHeatChange);
    }

    delete actkeyCheckHeatChange;
}

void NSSolver::FillActionKey(ActionKey *actkey, int action, int level)
{
    actkey->action     = action;
    actkey->solver     = NS_EQUATION;
    actkey->solverID   = this->GetIndex();
    actkey->kind       = SOLVER_BASED;
    actkey->level      = level;
}

RDouble NSSolver::ComputeCFL()
{
    int iterationStep;

    Param_NSSolver *parameters = GetControlParameters();
    //! Steady.
    int nMGLevel = parameters->GetNMGLevel();
    if (nMGLevel <= 1)
    {
        //! Finest grid.
        GlobalDataBase::GetData("outnstep", &iterationStep, PHINT, 1);
    }
    else
    {
        //! Coarse grid.
        GlobalDataBase::GetData("newnstep", &iterationStep, PHINT, 1);
    }

    int cflNstep    = parameters->GetCFLVaryStep();
    RDouble cflStart = parameters->GetCFLStart();
    RDouble cflEnd   = parameters->GetCFLEnd();

    RDouble cfl = cflStart + (cflEnd - cflStart) * MIN (iterationStep, cflNstep) / cflNstep;

    return cfl;
}

void NSSolver::SourceFlux(Grid *grid)
{
    Param_NSSolver *parameters = GetControlParameters();
    int nChemical = parameters->GetChemicalFlag();
    int nChemicalSource = parameters->GetNChemicalSource();

    if (nChemical == 1 && nChemicalSource == 1)
    {
        ChemicalSource(grid);
    }

    //! It is null for Structure grid calculation.
    Turb_Sengy(grid);

    //! Dual time step source.
    DualTimeSource(grid);
}

bool NSSolver::WantVisualField(Grid *grid)
{
    Param_NSSolver *parameters = GetControlParameters();
    int visualField = parameters->GetPlotFieldType();

    if (visualField == TEC_SPACE::BlockVisual)
    {
        int dimension = PHSPACE::GetDim();

        RDouble *boxMin = grid->GetMinBox();
        RDouble *boxMax = grid->GetMaxBox();

        const RDouble *visualBlockMin = parameters->GetLowerPlotFieldBox();
        const RDouble *visualBlockMax = parameters->GetUpperPlotFieldBox();

        if ((boxMin[0] > visualBlockMax[0]) || (boxMax[0] < visualBlockMin[0]))
        {
            return false;
        }

        if ((boxMin[1] > visualBlockMax[1]) || (boxMax[1] < visualBlockMin[1]))
        {
            return false;
        }

        if (dimension == THREE_D)
        {
            if ((boxMin[2] > visualBlockMax[2]) || (boxMax[2] < visualBlockMin[2]))
            {
                return false;
            }
        }
    }

    return grid->GetDim() == TWO_D || (grid->GetDim() == THREE_D && visualField != TEC_SPACE::BoundaryVisual);
}

void NSSolver::InitCoarseGridsFlow()
{
    Grid *grid = GetGrid();
    Grid *coarseGrid = grid->GetCoarseGrid();

    if(!coarseGrid)
    {
        return;
    }

    ComputeGamaAndTemperature(coarseGrid);

    Boundary(coarseGrid);

    ComputeGamaAndTemperature(coarseGrid);

    ComputeNodeValue(coarseGrid);

    //! dual time step source.
    DualTimeSource(grid);
}

void NSSolver::InitDependentVariables()
{
    Grid *fineGrid = GetGrid();

    //! These variables are need for boundary condition.
    ComputeGamaAndTemperature(fineGrid);

    Boundary(fineGrid);

    //! After boundary condition, the gama and temperature would be changed.
    ComputeGamaAndTemperature(fineGrid);

    Param_NSSolver *parameters = GetControlParameters();
    int isViscous = parameters->IsViscous();
    if (isViscous)
    {
        int nChemical = parameters->GetChemicalFlag();
        if (nChemical == 0)
        {
            ComputeViscousCoefficientWithoutChemical(fineGrid);
        }
        else
        {
            ComputeViscousCoefficient(fineGrid);
            #ifdef USE_HYFLOW_CODE
            GAS_SPACE::ComputeViscousLaminarCoefficient(StructGridCast(fineGrid));
            #endif
        }
    }

    ComputeNodeValue(fineGrid);
}

void NSSolver::RegisterCFDSolverInterfaceField()
{
    Param_NSSolver *parameters = GetControlParameters();
    int nEquation = GetNumberOfEquations();
    CFDSolver::RegisterInterfaceField("q", PHDOUBLE, nEquation);

    int nTemperatureModel = parameters->GetTemperatureModel();
    CFDSolver::RegisterInterfaceField("t", PHDOUBLE, nTemperatureModel);
    
    int viscousType = parameters->GetViscousType();
    int numberOfSpecies = parameters->GetNumberOfSpecies();
    int *visualVariablesType = postVisualWall->GetVisualVariablesType();
    int nWallVariables = postVisualWall->GetVisualVariablesNumber();
    int dataNumber = nWallVariables;
    for (int m = 0; m < nWallVariables; ++ m)
    {
        if (visualVariablesType[m] == VISUAL_WALL_NS)
        {
            dataNumber += numberOfSpecies - 1;
        }
    }
    CFDSolver::RegisterInterfaceField("wallCell", PHDOUBLE, 1);
    CFDSolver::RegisterInterfaceField("wallAircoefData", PHDOUBLE, dataNumber);
}

void NSSolver::ReleaseCFDSolverInterfaceField()
{
    CFDSolver::ReleaseInterfaceField();
}

void NSSolver::RegisterCFDSolverInterfaceFieldUnstruct()
{
    Param_NSSolver *parameters = GetControlParameters();
    
    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();

    int nEquation = GetNumberOfEquations();

    CFDSolver::RegisterInterfaceField("dqdx", PHDOUBLE, nEquation);
    CFDSolver::RegisterInterfaceField("dqdy", PHDOUBLE, nEquation);
    CFDSolver::RegisterInterfaceField("dqdz", PHDOUBLE, nEquation);

    CFDSolver::RegisterInterfaceField("limit", PHDOUBLE, nEquation);

    int ntmodel = parameters->GetTemperatureModel();
    CFDSolver::RegisterInterfaceField("dtdx", PHDOUBLE, ntmodel);
    CFDSolver::RegisterInterfaceField("dtdy", PHDOUBLE, ntmodel);
    CFDSolver::RegisterInterfaceField("dtdz", PHDOUBLE, ntmodel);

    CFDSolver::RegisterInterfaceField("cellQ" , PHDOUBLE, 1);
    CFDSolver::RegisterInterfaceField("cellCP", PHDOUBLE, 1);
    CFDSolver::RegisterInterfaceField("cellCF", PHDOUBLE, 1);

    CFDSolver::RegisterInterfaceField("dq", PHDOUBLE, nEquation);
}

void NSSolver::RegisterCFDSolverInterpointField()
{
    Param_NSSolver *parameters = GetControlParameters();

    int nLaminar  = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    
    int nEquation = GetNumberOfEquations();

    CFDSolver::RegisterInterpointField("qnode", PHDOUBLE, nEquation);
    
    int ntmodel = parameters->GetTemperatureModel();
    CFDSolver::RegisterInterpointField("tnode", PHDOUBLE, ntmodel);
}

void NSSolver::RegisterCFDSolverOversetField()
{
    Param_NSSolver *parameters = GetControlParameters();
    int nEquation = GetNumberOfEquations();
    CFDSolver::RegisterOversetField("q", PHDOUBLE, nEquation);
    int nTemperatureModel = parameters->GetTemperatureModel();
    CFDSolver::RegisterOversetField("t", PHDOUBLE, nTemperatureModel);
}

void NSSolver::AllocateOversetInterfaceVar(int nIOversetFace, Data_ParamFieldSuite *dataStore)
{
    if (nIOversetFace <= 0)
    {
        return;
    }

    Param_NSSolver *parameters = GetControlParameters();
    int nTemperatureModel = parameters->GetTemperatureModel();
    int nEquation = GetNumberOfEquations();

    RDouble **q = NewPointer2 <RDouble> (nEquation, nIOversetFace);
    dataStore->UpdateDataPtr("q", q);

    RDouble **t = NewPointer2 <RDouble> (nTemperatureModel, nIOversetFace);
    dataStore->UpdateDataPtr("t", t);
}

void NSSolver::DeAllocateOversetInterfaceVar(Data_ParamFieldSuite *dataStore)
{
    RDouble **q = reinterpret_cast <RDouble **> (dataStore->GetDataPtr("q"));
    DelPointer2(q);

    RDouble **t = reinterpret_cast <RDouble **> (dataStore->GetDataPtr("t"));
    DelPointer2(t);
}

void NSSolver::PrepareOversetInterfaceData(Data_ParamFieldSuite *dataStore, InterfaceDataProxy *interfaceDataProxy)
{
    Param_NSSolver *parameters = GetControlParameters();
    int nTemperatureModel = parameters->GetTemperatureModel();
    int nEquation = GetNumberOfEquations();

    RDouble **q = reinterpret_cast< RDouble ** > (dataStore->GetDataPtr("q"));

    vector <RDouble **> &vectorData = interfaceDataProxy->GetVectorData();
    vector <int> & vectorDimension = interfaceDataProxy->GetVectorDimension();

    vectorDimension.push_back(nEquation);
    vectorData.push_back(q);
    
    RDouble **t = reinterpret_cast <RDouble **> (dataStore->GetDataPtr("t"));
    
    vectorDimension.push_back(nTemperatureModel);
    vectorData.push_back(t);
}

void NSSolver::GetResidual(ActionKey *actkey)
{

}

void SetInviscidSchemeParameters(InviscidSchemeParameter *invSchemePara, FaceProxy *faceProxy, Param_NSSolver *parameters)
{
    int nm = parameters->GetNSEquationNumber();
    int nLaminar = parameters->GetLaminarNumber();
    int nChemical = parameters->GetChemicalFlag();
    int nTemperatureModel = parameters->GetTemperatureModel();
    int nEquation = nLaminar + nChemical + nTemperatureModel - 1;
    int nPrecondition = parameters->GetIfLowSpeedPrecon();
    int nIdealState = parameters->GetIfIdealGasState();
    int nNitrogenIndex = parameters->GetIndexOfNitrogen();
    int nElectronIndex = parameters->GetIndexOfElectron();

    int RoeEntropyFixMethod  = parameters->GetRoeEntropyFixMethod();
    RDouble entrFixCoefLeft  = parameters->GetRoeEntropyFixCoef1();
    RDouble entrFixCoefRight = parameters->GetRoeEntropyFixCoef2();
    RDouble refMachNumber = parameters->GetRefMachNumber();

    //! Scalar parameters.
    invSchemePara->Setnm(nm);
    invSchemePara->SetNumberOfLaminar(nLaminar);
    invSchemePara->SetNumberOfChemical(nChemical);
    invSchemePara->SetNumberOfTemperatureModel(nTemperatureModel);
    invSchemePara->SetNumberOfTotalEquation(nEquation);
    invSchemePara->SetIfPrecondition(nPrecondition);
    invSchemePara->SetIfIdealState(nIdealState);
    invSchemePara->SetEntropyFixMethod(RoeEntropyFixMethod);
    invSchemePara->SetMachNumber(refMachNumber);
    invSchemePara->SetEntropyFixCoefficients(entrFixCoefLeft, entrFixCoefRight);

    invSchemePara->SetLength(faceProxy->size());
    invSchemePara->SetIndexOfNitrogen(nNitrogenIndex);
    invSchemePara->SetIndexOfElectron(nElectronIndex);

    //! Pointer parameters.
    invSchemePara->SetLeftAndRightQ(faceProxy->GetQL(), faceProxy->GetQR());
    invSchemePara->SetLeftAndRightTemperature(faceProxy->GetLeftTemperature(), faceProxy->GetRightTemperature());
    invSchemePara->SetFlux(faceProxy->GetFlux());
    invSchemePara->SetLeftAndRightPressureCoefficient(faceProxy->GetPressureCoefficientL(), faceProxy->GetPressureCoefficientR());
    invSchemePara->SetLeftAndRightGama(faceProxy->GetGamaL(), faceProxy->GetGamaR());

    GeomProxy *geomProxy = faceProxy->GetGeomProxy();

    invSchemePara->SetFaceNormal(geomProxy->GetFaceNormalX(), geomProxy->GetFaceNormalY(), geomProxy->GetFaceNormalZ());
    invSchemePara->SetFaceArea(geomProxy->GetFaceArea());
    invSchemePara->SetFaceVelocity(geomProxy->GetFaceVelocity());
}

}
