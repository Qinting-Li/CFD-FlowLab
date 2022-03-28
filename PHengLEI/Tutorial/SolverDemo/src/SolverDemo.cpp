//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      SolverDemo.cpp
//! @brief     It tests whether the whole solve process from grid conversion to calculation can work well,
//!			   furthermore, gives users some guidance.
//! @author    Meng Liyuan.

#include "PHMpi.h"
#include "TK_Log.h"
#include "TK_Time.h"
#include "Region.h"
#include "CFDSolver.h"
#include "Controller.h"
#include "Glb_Dimension.h"
#include "Pre_CGNSToPHengLEI_Unstruct.h"
#include "Task.h"
#include "PHIO.h"
#include "Constants.h"
#include "Gas.h"
#include "RealGas.h"
#include "Zone.h"
#include "OverLappingGrid.h"
#include "MBReactor.h"

using namespace PHSPACE;
using namespace PHMPI;
using namespace GAS_SPACE;
int main()
{
    //! Perform the initialization operation of PHengLEI.
    ReadControlParameters();
    Initialization();

    InitGlobalValuesSolver();

    SetDim(TWO_D);	

    //! Define output file which to store the results.
    fstream file;
    ParallelOpenFile(file, "./out/SolverDemoOutput.txt", ios::out);

    string cgnsGridFileName = "./in/cylinder.cgns";
    string gridFileName = "./out/cylinder.fts";
    GlobalDataBase::UpdateData("gridfile", &gridFileName, PHSTRING, 1);

    //! Read in grid(CGNS) and convert.
    Pre_CGNSToPHengLEI_Unstruct cgnsGridConversion(cgnsGridFileName);
    cgnsGridConversion.Run();

    int nBlocks = cgnsGridConversion.GetNumberofBlocks();
    Grid ** grids = cgnsGridConversion.GetGrids();
    ConstructOverlappingInfo(grids,nBlocks,0);
    DumpGrid(gridFileName, grids, nBlocks);
    cgnsGridConversion.Clear();
    FreeBlockData();

    //! Read in the original grid
    Region * region = new Region();
    //region->ReadOrdinaryGrid(gridFileName);
    region->ReadGrid();

    region->UpdateAllZonesPara();
    region->InitGeometry();
    region->SwapGeommetry();

    // Dump out information.
    int nZones = PHMPI::GetNumberofGlobalZones();
    ostringstream oss;
    oss << "  Number of Zones in processor " << GetCurrentProcessorID() << ": " << nZones << endl;
    for ( int izone = 0; izone < nZones; ++ izone )
    {
        Zone * zone = region->GetZone(izone); 
        if ( !zone ) continue;
        Grid * grid = zone->GetGeometry()->GetGrid(0);

        double minDis, maxDis;
        grid->GetMinMaxDS(minDis, maxDis);
        oss << "    For Grid " << izone << ", minDis = " << minDis << ", " << "maxDis = " << maxDis << " " << endl;
        oss << "    nTNode = " << grid->GetNTotalNode() << ", nTCell = " << grid->GetNTotalCell() << endl;
    }
    PHSPACE::WriteLogFile(oss.str());

    // Set parameters,such as Euler/Laminar/Turbulence,discretization scheme,limiter and so on.	
    int iviscous= 0; 
    string visname = "Euler";
    string uns_scheme_name = "vanleer";
    string uns_limiter_name="nolim";
    GlobalDataBase::UpdateData("iviscous",&iviscous,PHINT,1);
    GlobalDataBase::UpdateData("visname",&visname,PHSTRING,1);	
    GlobalDataBase::UpdateData("uns_scheme_name",&uns_scheme_name,PHSTRING,1);
    GlobalDataBase::UpdateData("uns_limiter_name",&uns_limiter_name,PHSTRING,1);

    
    // Create solver,and load to solution work flow.

    for ( int izone = 0; izone < nZones; ++ izone )
    {
        Zone * zone = region->GetZone(izone); 
        if ( !zone ) continue;


        MYReactor *s= new MYReactor( NULL, NULL );
        CFDSolver * cfdSolver = s->GetCFDSolver(izone);

        //CFDSolver * cfdSolver = new CFDSolver();

        cfdSolver->SetGeometry(zone->GetGeometry());
        cfdSolver->SetIndex(0);

        zone->AddSolver(cfdSolver);
        zone->AddGridToGlobal(izone);
    }
    GlobalSolve();


    // Solver iteration .
    int maxInteration = 10;
    int iteration     = 0 ;
    while ( iteration < maxInteration )
    {
        iteration += 1;

        int nSolver = Zone::GetNumberOfSolvers();

        for ( int isolver = 0; isolver < nSolver; ++ isolver )
        {
            Controller * controller = new Controller(isolver, region);

            controller->SolveSteadyField();

            controller->PostSolve();

            delete controller;
        }
    }
    delete region;

    file << " The whole solve process works well !" << endl;

    //! Close the output file and exit PHengLEI.
    ParallelCloseFile(file);

    PHMPI::FreeBlockData();
    TIME_SPACE::WriStatTime();
    PH_Finalize();

    return 0;
}