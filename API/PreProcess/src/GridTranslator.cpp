#include <cmath>
#include "GridTranslator.h"
#include "Pre_GridConversion.h"
#include "MixGrid.h"
#include "Task_WriteBinaryTecplot.h"
#include "BC_Struct.h"
#include "MultigridFactory.h"
#include "PHIO.h"

using namespace std;

namespace PHSPACE
{

GridTranslator::GridTranslator()
{
    GlobalDataBase::GetData("codeOfOversetGrid",&codeOfOversetGrid,PHINT,1);

    gridGroup = new vector< StructGrid * >();

    structFaceGroup = new vector< StructFace * >();
}

GridTranslator::GridTranslator(Grid **& strgrids, int &nBlocks)
{
    numberOfZones = nBlocks;

    GlobalDataBase::GetData("codeOfOversetGrid",&codeOfOversetGrid,PHINT,1);

    gridGroup = new vector< StructGrid * >();
    gridGroup->resize( nBlocks, NULL );

    for ( int iZone = 0; iZone < numberOfZones; ++ iZone )
    {
        //GridID *index = new GridID(iZone);
        //StructGrid * grid;
        //grid = StructGridCast(CreateGridGeneral(STRUCTGRID, index, 0, GetDim()));
        ( * gridGroup )[iZone] = StructGridCast(strgrids[iZone]);
        //grid->CopyGrid(StructGridCast(strgrids[iZone]));
    }

    structFaceGroup = new vector< StructFace * >();
}

GridTranslator::~GridTranslator()
{
    for ( int iZone = 0; iZone < numberOfZones; ++ iZone )
    {
        delete ( * structFaceGroup )[iZone];
    }
    delete structFaceGroup;	structFaceGroup = NULL;

    /*for ( int iZone = 0; iZone < numberOfZones; ++ iZone )
    {
        delete ( * gridGroup )[iZone];
    }*/	
    delete gridGroup;	gridGroup = NULL;
}

void GridTranslator::Run()
{
    ProcessBoundaryCondition();

    GenerateStructFaceGroup();

    ColorateStructFaceGroup();

    CollectInterfaceInformation();

    //CheckGrid();

    //OutputMultiZoneComputationalGrids();

    OutputGrid();

    ProcessOversetGrid();
}

void GridTranslator::ReadGridgenFiles()
{
    //ReadGridgenCoordinateFile();

    //ReadGridgenBoundaryFile();

    //ProcessBoundaryCondition();

    //return;
}

void GridTranslator::ReadGridgenCoordinateFile()
{
}

void GridTranslator::ReadGridgenBoundaryFile()
{
}

void GridTranslator::ProcessOversetGrid()
{
    return;
}

void GridTranslator::ProcessBoundaryCondition()
{
    for ( int iZone = 0; iZone < numberOfZones; ++ iZone )
    {
        CompositeBCRegion * compositeBCRegion = ( * gridGroup )[iZone]->GetCompositeBCRegion();
        
        compositeBCRegion->ProcessBCInfo();
    }
}

void GridTranslator::GenerateStructFaceGroup()
{
    structFaceGroup->resize( numberOfZones, NULL );

    for ( int iZone = 0; iZone < numberOfZones; ++ iZone )
    {
        ( * structFaceGroup )[iZone] = new StructFace( ( * gridGroup )[iZone] );
    }
}

void GridTranslator::ColorateStructFaceGroup()
{
    for ( int iZone = 0; iZone < numberOfZones; ++ iZone )
    {
        
        StructFace * structFace = ( * structFaceGroup )[iZone];
        structFace->SetStructFaceGroup( this->structFaceGroup );
        structFace->PaintSelf();
    }

    for ( int iZone = 0; iZone < numberOfZones; ++ iZone )
    {
        StructFace * structFace = ( * structFaceGroup )[iZone];
        structFace->PaintNeighbor();
    }
}

void GridTranslator::CollectInterfaceInformation()
{
    for ( int iZone = 0; iZone < numberOfZones; ++ iZone )
    {
        cout << "iZone = " << iZone << endl;

        int numberOfInterfaces = ( * structFaceGroup )[iZone]->ComputeNumberOfInterfaces();

        cout << " numberOfInterfaces = " << numberOfInterfaces << endl;

        StructGrid * structuredGrid = ( * gridGroup )[iZone];

        structuredGrid->SetNIFace( numberOfInterfaces );

        if ( numberOfInterfaces == 0 )
        {
            structuredGrid->SetInterfaceInfo( NULL );
        }
        else
        {
            structuredGrid->SetInterfaceInfo( new InterfaceInfo( numberOfInterfaces ) );
        }

        InterfaceInfo * interfaceInfo = structuredGrid->GetInterfaceInfo();

        ( * structFaceGroup )[iZone]->CollectInterfaceInformation( interfaceInfo );
    }
}

void GridTranslator::OutputGrid()
{
    string out_gfile = "grid.fts";
    int tasktype = GetTaskCode();
    if (tasktype == PHSPACE::PARTITION_GRID)
    {
        GlobalDataBase::GetData("partition_grid_file",&out_gfile,PHSTRING,1);
    }
    else if(tasktype == PHSPACE::CREATE_GRID || tasktype == SPLIT_OVERSET_GRID)
    {
        GlobalDataBase::GetData("out_gfile",&out_gfile,PHSTRING,1);
    }

    Grid **grids = new Grid *[numberOfZones];
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        grids[iZone] = (*gridGroup)[iZone];
    }

    DumpGrid(out_gfile, grids, numberOfZones);
}

void GridTranslator::OutputMultiZoneComputationalGrids()
{
    fstream file;
    int tasktype = GetTaskCode();
    if (tasktype == PHSPACE::PARTITION_GRID)
    {
        string partition_grid_file = "grid.fts";
        GlobalDataBase::GetData("partition_grid_file",&partition_grid_file,PHSTRING,1);
        OpenSerialFile(file, partition_grid_file,ios_base::out|ios_base::binary|ios_base::trunc);
    }
    else if(tasktype == PHSPACE::CREATE_GRID || tasktype == SPLIT_OVERSET_GRID)
    {
        string out_gfile = "grid.fts";
        GlobalDataBase::GetData("out_gfile",&out_gfile,PHSTRING,1);
        OpenSerialFile(file, out_gfile.c_str(),ios_base::out|ios_base::binary|ios_base::trunc);
    }

    using namespace PHMPI;

    CreateZoneProcessorID(numberOfZones);
    CreateZoneGridID(numberOfZones);
    CreateZoneGridType(numberOfZones);

    int * zoneProcessorIndexContainer = GetZoneProcessorID();
    int * zoneIndexContainer	= GetZoneGridID();
    int * zoneTypeContainer		= GetZoneGridType();
    int * block_proc_inp		= GetZoneProcessorID_INP();

    for ( int iZone = 0; iZone < numberOfZones; ++ iZone )
    {
        zoneProcessorIndexContainer[iZone] = iZone;
        zoneIndexContainer[iZone] = iZone;
        zoneTypeContainer[iZone] = ( * gridGroup )[iZone]->Type();
    }

    if ( block_proc_inp )
    {
        for ( int izone = 0; izone < numberOfZones; ++ izone )
        {
            SetZoneProcessorID(izone, block_proc_inp[izone]);
        }
    }

    PHWrite( file, numberOfZones );
    PHWrite( file, & zoneProcessorIndexContainer[0], numberOfZones );
    PHWrite( file, & zoneIndexContainer[0], numberOfZones );
    PHWrite( file, & zoneTypeContainer[0], numberOfZones );

    for ( int iZone = 0; iZone < numberOfZones; ++ iZone )
    {
        ( * gridGroup )[iZone]->WriteGrid( file );
    }

    CloseFile( file );
}

void GridTranslator::CheckGrid()
{
    for ( int izone = 0; izone < numberOfZones; ++ izone )
    {
        StructGrid *grid = StructGridCast(( * gridGroup )[izone]);
        InterfaceInfo *iinfo = grid->GetInterfaceInfo();

        int nIFace = 0;
        if ( iinfo )
        {
            nIFace = iinfo->GetNIFace();
        }

        if (nIFace == 0) continue;

        int *interFace2ZoneID       = iinfo->GetInterFace2ZoneID();
        int *interFace2InterFaceID  = iinfo->GetInterFace2InterFaceID();
        int *interFace2BoundaryFace = iinfo->GetInterFace2BoundaryFace();

        int target_zone, itface;
        int is, js, ks, nsurf, s_lr, it, jt, kt, t_lr;
        RDouble s_xfc, s_yfc, s_zfc, t_xfc, t_yfc, t_zfc;
        int il1, jl1, kl1;
  
        for ( int iface = 0; iface < nIFace; ++ iface )
        {
            target_zone = interFace2ZoneID[iface];
            if (target_zone < izone) continue;
            
            grid->GetSourceIndexIJK_Nsurf_LR(iface,1,is,js,ks,nsurf,s_lr);
            if (s_lr == 1)
            {
                GetNsurfIndex(nsurf, il1, jl1, kl1);
                is = is + il1;
                js = js + jl1;
                ks = ks + kl1;
            }
            grid->FaceCoor(is,js,ks,nsurf,s_xfc,s_yfc,s_zfc);

            itface = interFace2InterFaceID[iface];

            if (target_zone == izone)
            {
                grid->GetSourceIndexIJK_Nsurf_LR(itface,1,it,jt,kt,nsurf,t_lr);
                if (t_lr == 1)
                {
                    GetNsurfIndex(nsurf, il1, jl1, kl1);
                    it = it + il1;
                    jt = jt + jl1;
                    kt = kt + kl1;
                }
                grid->FaceCoor(it,jt,kt,nsurf,t_xfc,t_yfc,t_zfc);	
            }
            else
            {
                StructGrid *tgrid = StructGridCast(( * gridGroup )[target_zone]);
                tgrid->GetSourceIndexIJK_Nsurf_LR(itface,1,it,jt,kt,nsurf,t_lr);
                if (t_lr == 1)
                {
                    GetNsurfIndex(nsurf, il1, jl1, kl1);
                    it = it + il1;
                    jt = jt + jl1;
                    kt = kt + kl1;
                }
                tgrid->FaceCoor(it,jt,kt,nsurf,t_xfc,t_yfc,t_zfc);
            }
            
            if ( (fabs( s_xfc - t_xfc ) > 1.0e-8) || (fabs( s_yfc - t_yfc ) > 1.0e-8) || (fabs( s_zfc - t_zfc ) > 1.0e-8) )
            {
                TK_Exit::ExceptionExit(" Grid Link Error! ");
            }
        }		
    }
}

}


