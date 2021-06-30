#include "GlobalDataBase.h"
#include "CFDSimulation.h"
#include "Simulation.h"
#include "TK_Exit.h"
#include "Constants.h"
#include "PHIO.h"
#include "TK_Log.h"

using namespace std;

namespace PHSPACE
{
void RunPHengLEI()
{
    //! Initialize global variables and MPI.
    PHMPI::Initialization();

    ShowWelcomeTitle();

    //! Read control parameters and do self-checking.
    ReadControlParameters();

    CompleteParamInfo();

    ControlParametersPreliminarySelfCheck();

    InitFantasyMPIForSepMode(); //! This function is useless and will be abandoned later.

    //! Start CFD simulations.
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

    int flowInitStep = PHSPACE::GlobalDataBase::GetIntParaFromDB("flowInitStep");
    int maxSimuStep  = PHSPACE::GlobalDataBase::GetIntParaFromDB("maxSimuStep");
    if (flowInitStep > maxSimuStep)
    {
        //! If the maximum simulation step is small than the initialization step,
        //! do not init then.
        flowInitStep = 0;
        GlobalDataBase::UpdateData("flowInitStep", &flowInitStep, PHINT, 1);
    }

    int visualVariablesSize = GlobalDataBase::GetSizeOfData("visualVariables");
    int nVisualVariables    = visualVariablesSize;
    GlobalDataBase::UpdateData("nVisualVariables", &nVisualVariables, PHINT, 1);

    int *visualVariables = new int [nVisualVariables];
    GlobalDataBase::GetData("visualVariables", visualVariables, PHINT, nVisualVariables);
    int viscousType = GlobalDataBase::GetIntParaFromDB("viscousType");
    int iLES = GlobalDataBase::GetIntParaFromDB("iLES");
    if (viscousType == 0)
    {
        for (int i = 0; i < nVisualVariables; ++ i)
        {
            if (visualVariables[i] == 7 || visualVariables[i] == 8 || visualVariables[i] == 9 || visualVariables[i] == 10)
            {
                ostringstream oss;
                oss << "The array \"visualVariables\" is error. The euler solver cannot export the data related to viscosity." << endl;
                TK_Exit::ExceptionExit(oss);
            }
        }
    }
    else if (viscousType == 1 && iLES == 0)
    {
        for (int i = 0; i < nVisualVariables; ++ i)
        {
            if (visualVariables[i] == 8)
            {
                ostringstream oss;
                oss << "The array \"visualVariables\" is error. The laminar solver cannot export the data related to turbulence." << endl;
                TK_Exit::ExceptionExit(oss);
            }
        }
    }
    else if (viscousType == 3)
    {
        for (int i = 0; i < nVisualVariables; ++ i)
        {
            if (visualVariables[i] == 35 || visualVariables[i] == 36)
            {
                ostringstream oss;
                oss << "The array \"visualVariables\" is error. The solver cannot export the transition data with SA model." << endl;
                TK_Exit::ExceptionExit(oss);
            }
        }
    }
    else if (viscousType == 4)
    {
        int transitionType = GlobalDataBase::GetIntParaFromDB("transitionType");
        for (int i = 0; i < nVisualVariables; ++ i)
        {
            if (transitionType /= 2)
            {
                if (visualVariables[i] == 35 || visualVariables[i] == 36)
                {
                    ostringstream oss;
                    oss << "The array \"visualVariables\" is error. The solver cannot export the transition data without transition type." << endl;
                    TK_Exit::ExceptionExit(oss);
                }
            }
        }
    }

    delete [] visualVariables;

    int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
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
            WriteLogFile (oss.str());

            maxSubIter = minSubIter;
            GlobalDataBase::UpdateData("max_sub_iter", &maxSubIter, PHINT, 1);
        }
    }
#ifdef USE_OVERSET
    int codeOfOversetGrid = GlobalDataBase::GetIntParaFromDB("codeOfOversetGrid");
    if( codeOfOversetGrid == 1 )
    {
        int dimension = GlobalDataBase::GetIntParaFromDB("ndim");
        if ( dimension < 3 )
        {
            ostringstream oss;
            oss << "The variable \"ndim\" must be equal to three with overset grid open ." << endl;
            TK_Exit::ExceptionExit(oss);
        }
    }
#endif
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
        case REPOSITORY:
            break;
        case OVERSETGRID_VIEW:
            break;
        case OVERSET_CONFIG:
            break;
        case CONVERT_PARTITION:
            break;
        case HO_SOLVER:
            break;
        case INTEGRATIVE_SOLVER:
            ParametersSelfCheckForCFD();
            break;
        case SPLIT_OVERSET_GRID:
            break;
        case POST_PROCESSING:
            break;
        case SPECDIFFHYB_Solver:
            break;
        default:
            TK_Exit::UnexpectedVarValue("taskType = ", taskType);
            break;
    }
}

}


