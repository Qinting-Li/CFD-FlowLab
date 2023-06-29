#include "GlobalDataBase.h"
#include "CFDSimulation.h"
#include "Simulation.h"
#include "TK_Exit.h"
#include "Constants.h"
#include "PHIO.h"
#include "TK_Log.h"
#include "Residual.h"
#include "Post_Visual.h"

using namespace std;

namespace PHSPACE
{
void RunPHengLEI()
{
    //! Initialize global variables and MPI.

    PHMPI::Initialization();

    //! Read control parameters and do self-checking.
    ReadControlParameters();

    ShowWelcomeTitle();

    CompleteParamInfo();
       
    ControlParametersPreliminarySelfCheck();

    //! Start CFD simulation.
    int dimension = GlobalDataBase::GetIntParaFromDB("ndim");
    Simulation *simulation = new Simulation(dimension);

    simulation->Start();

    delete simulation;
}

void ShowWelcomeTitle()
{
    using namespace PHMPI;
    PH_Barrier();

    int myid = GetCurrentProcessorID();
    if (myid == PHMPI::GetServerProcessorID())
    {
        ostringstream oss;
        oss << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
        oss << "$#                                                     #$" << endl;
        oss << "$#  PPPPP H   H  EEEE N     N  GGGGG  L     EEEEE III  #$" << endl;
        oss << "$#  P   P H   H  E    N N   N  G      L     E      I   #$" << endl;
        oss << "$#  PPPPP HHHHH  EEEE N  N  N  G  GG  L     EEEEE  I   #$" << endl;
        oss << "$#  P     H   H  E    N   N N  G   G  L     E      I   #$" << endl;
        oss << "$#  P     H   H  EEEE N     N  GGGGG  LLLLL EEEEE III  #$" << endl;
        oss << "$#                                                     #$" << endl;
        oss << "$#=====================================================#$" << endl;
        oss << "$#                                                     #$" << endl;
        oss << "$# Platform for Hybrid Engineering Simulation of Flows #$" << endl;
        oss << "$#                                                     #$" << endl;
        oss << "$#       Computational Aerodynamics Institute          #$" << endl;
        oss << "$#  China Aerodynamics Research and Development Center #$" << endl;
        oss << "$#            (C) Copyright, Since 2010                #$" << endl;
        oss << "$#                  PHengLEI 2306                      #$" << endl;
        oss << "$#                                                     #$" << endl;
        oss << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;

        cout << oss.str();
    }
}

void ParametersSelfCheckForCFD()
{
    //! Is in Multi-Grid initialization state.
    int isInMGIniting = 0;
    GlobalDataBase::UpdateData("isInMGIniting", &isInMGIniting, PHINT, 1);

    //! Is in fine-Grid initialization state.
    int isInIniting = 1;
    GlobalDataBase::UpdateData("isInIniting", &isInIniting, PHINT, 1);

    int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");

    int flowInitStep = PHSPACE::GlobalDataBase::GetIntParaFromDB("flowInitStep");
    int maxSimuStep = PHSPACE::GlobalDataBase::GetIntParaFromDB("maxSimuStep");
    if (flowInitStep > maxSimuStep || isUnsteady)
    {
        //! If the maximum simulation step is small than the initialization step,
        //! do not init then.
        flowInitStep = 0;
        GlobalDataBase::UpdateData("flowInitStep", &flowInitStep, PHINT, 1);
    }

    int visualVariablesSize = GlobalDataBase::GetSizeOfData("visualVariables");
    int nVisualVariables = visualVariablesSize;
    GlobalDataBase::UpdateData("nVisualVariables", &nVisualVariables, PHINT, 1);

    int *visualVariables = new int[nVisualVariables];
    GlobalDataBase::GetData("visualVariables", visualVariables, PHINT, nVisualVariables);
    int viscousType = GlobalDataBase::GetIntParaFromDB("viscousType");
    int iLES = GlobalDataBase::GetIntParaFromDB("iLES");
    for (int i = 0; i < nVisualVariables; ++ i)
    {
        if (visualVariables[i] >= VISUAL_VISCOSITY_LAMINAR && visualVariables[i] < VISUAL_CP)
        {
            if (viscousType == INVISCID)
            {
                ostringstream oss;
                oss << "The array \"visualVariables\" is error. The euler solver cannot export the data related to viscosity." << endl;
                TK_Exit::ExceptionExit(oss);
            }
        }

        if (visualVariables[i] == VISUAL_VISCOSITY_TURBULENT)
        {
            if (viscousType == LAMINAR && iLES == NOLES_SOLVER)
            {
                ostringstream oss;
                oss << "The array \"visualVariables\" is error. The laminar solver cannot export the data related to turbulence." << endl;
                TK_Exit::ExceptionExit(oss);
            }
        }

        if (visualVariables[i] == VISUAL_INTERMITTENCY || visualVariables[i] == VISUAL_MOMENTUMTHICKREYNOLDS)
        {
            if (viscousType == ONE_EQU)
            {
                ostringstream oss;
                oss << "The array \"visualVariables\" is error. The solver cannot export the transition data with SA model." << endl;
                TK_Exit::ExceptionExit(oss);
            }
            else if (viscousType == TWO_EQU)
            {
                int transitionType = GlobalDataBase::GetIntParaFromDB("transitionType");
                if (transitionType != IREGAMA)
                {
                    ostringstream oss;
                    oss << "transitionType = " << transitionType << endl;
                    oss << "The array \"visualVariables\" is error. The solver cannot export the transition data with transition type Not Equal 2." << endl;
                    TK_Exit::ExceptionExit(oss);
                }
            }
        }
    }

    delete[] visualVariables;

    if (isUnsteady)
    {
        int minSubIter = GlobalDataBase::GetIntParaFromDB("min_sub_iter");
        if (minSubIter <= 0)
        {
            ostringstream oss;
            oss << "The variable \"min_sub_iter\" must be greater than 1." << endl;
            TK_Exit::ExceptionExit(oss);
        }

        int maxSubIter = GlobalDataBase::GetIntParaFromDB("max_sub_iter");
        if (minSubIter > maxSubIter)
        {
            ostringstream oss;
            oss << "Warning: The variable \"min_sub_iter\" cannot be greater than \"max_sub_iter\"! Now let them be equal: max_sub_iter = min_sub_iter!" << endl;
            PrintToWindow(oss.str());
            WriteLogFile(oss.str());

            maxSubIter = minSubIter;
            GlobalDataBase::UpdateData("max_sub_iter", &maxSubIter, PHINT, 1);
        }
    }

    //! If use Block-LUSGS, ROE scheme must be selected.
    int tscheme = GlobalDataBase::GetIntParaFromDB("tscheme");
    if (tscheme == Block_LU_SGS)
    {
        string uns_scheme_name = GlobalDataBase::GetStrParaFromDB("uns_scheme_name");
        if (uns_scheme_name.substr(0, 3) != "roe")
        {
            ostringstream oss;
            oss << "Error: if BLU is used, ROE scheme must be selected!" << endl;
            PrintToWindow(oss.str());
            WriteLogFile(oss.str());
            TK_Exit::ExceptionExit(oss);
        }
    }
    //! Judge if use wenn scheme, provided by Li Qin.
    int isWennScheme = 0;
    string str_limiter_name = GlobalDataBase::GetStrParaFromDB("str_limiter_name");
    if (str_limiter_name == "weno3_js"
        || str_limiter_name == "wenn3_zm"
        || str_limiter_name == "wenn3_zes2"
        || str_limiter_name == "wenn3_zes3"
        || str_limiter_name == "wenn3_prm211")
    {
        isWennScheme = 1;
    }
    GlobalDataBase::UpdateData("isWennScheme", &isWennScheme, PHINT, 1);
}

void ControlParametersPreliminarySelfCheck()
{
    int taskType = GetTaskCode();

    switch (taskType)
    {
        case SOLVE_FIELD:
            ParametersSelfCheckForCFD();
            break;
        case CREATE_GRID:
            break;
        case CAL_WALL_DIST:
            break;
        case PARTITION_GRID:
            break;
        case OVERSETGRID_VIEW:
            break;
        case OVERSET_CONFIG:
            break;
        case HO_SOLVER:
            break;
        case INTEGRATIVE_SOLVER:
            ParametersSelfCheckForCFD();
            break;
        case POST_PROCESSING:
            break;
        case SPECDIFFHYB_SOLVER:
            break;
        case SPEC_SOLVER:
            break;
 #ifdef USE_LBMSolverMPI
        case LBM_SOLVER_MPI:
            break;
        #endif
        #ifdef USE_LBMSolverOMP
        case LBM_SOLVER_OMP:
            break;
        #endif
        default:
            TK_Exit::UnexpectedVarValue("taskType = ", taskType);
            break;
    }
}

}



