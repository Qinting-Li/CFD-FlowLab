//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      BCStruct.cpp
//! @brief     It tests the functions of structured grid boundary condition , 
//!			   check whether they can run as required or not.
//! @author    Meng Liyuan.

#include "TK_Log.h"
#include "Glb_Dimension.h"
#include "BC_Struct.h"
#include "Region.h"
#include "PHIO.h"
#include "GlobalDataBase.h"
#include "TK_Time.h"

using namespace PHSPACE;
using namespace PHMPI;
int main()
{
    //! Record time.
    setvbuf( stdout, NULL, _IONBF, 0 );
    double *ptr_time = NULL;
    PHSPACE::TIME_SPACE::GetWallTimeForResidual();
    PHSPACE::TIME_SPACE::GetWallTime();
    PHSPACE::TIME_SPACE::GetCPUTime(ptr_time);

    //! Perform the initialization operation of PHengLEI.
    ReadControlParameters();
    Initialization();

    //! Define output file which to store the results.
    fstream ofile;
    ParallelOpenFile(ofile, "./out/BC_StrOutput.txt", ios::out);

    //! Read in the structured grid.
    Region *region = new Region();
    region->ReadGrid();

    //! Extract and dump out the boundary condition information of the structured grid.
    int nZones = GetNumberofGlobalZones();
    for (int iGrid = 0; iGrid < nZones; ++ iGrid)	
    {
        ofile << "Structured Grid CellMetrics of Zone "<< iGrid << " is: " << endl;
        Grid * grid = GetGrid(iGrid);

        CompositeBCRegion * composite_bcregion = StructGridCast(grid)->GetCompositeBCRegion();
        int  nBCRegion  = composite_bcregion->GetnBCRegion();
        int  nBoundFace = composite_bcregion->GetNBoundFace();
        int  nIFace     = composite_bcregion->GetNIFace();

        for ( int ibcregion = 0; ibcregion < nBCRegion; ibcregion ++ )
        {
            BCRegion * BCRegion = composite_bcregion->GetBCRegion(ibcregion);
            int bctype = BCRegion->GetBCType();
            ofile << "   The boundary condition type of region "<<ibcregion<< " is: "<< bctype << endl;
            if ( bctype >= 0 ) continue;

            int dimension = GetDim();
            int m = 2;

            int * s_st = BCRegion->GetStartPoint();
            int * s_ed = BCRegion->GetEndPoint();
            for ( int i = 0; i < dimension; i ++ )
            {
                ofile << "   The start point coordinate for " << i << " dimension of source region "<<ibcregion<<" is: " << s_st[i]<< endl;
                ofile << "   The end point coordinate for "   << i << " dimension of source region "<<ibcregion<<" is: "<< s_ed[i]<< endl;	
            }

            int s_st2 = BCRegion->GetStartPoint(m);
            int s_ed2 = BCRegion->GetEndPoint(m);
            ofile << "   The start point coordinate for 2_th dimension of source region "<<ibcregion<<" is: " << s_st2<< endl;
            ofile << "   The end point coordinate for 2_th dimension of source region "  <<ibcregion<<" is: "<< s_ed2<< endl;	

            int * t_st = BCRegion->GetTargetStart();
            int * t_ed = BCRegion->GetTargetEnd();
            for ( int i = 0; i < dimension; i ++ )
            {
                ofile << "   The start point coordinate for " << i << " dimension of target region " <<ibcregion<<" is: " << t_st[i]<< endl;
                ofile << "   The end point coordinate for "   << i << " dimension of target region " <<ibcregion<<" is: "<< t_ed[i]<< endl;	
            }

            int t_st2 = BCRegion->GetTargetStart(m);
            int t_ed2 = BCRegion->GetTargetEnd(m);
            ofile << "   The start point coordinate for 2_th dimension of target region "<<ibcregion<<" is: " << t_st2<< endl;
            ofile << "   The end point coordinate for 2_th dimension of target region "  <<ibcregion<<" is: "<< t_ed2<< endl;	

            int nbt  = BCRegion->GetTargetRegionBlock();
            int nbs  = BCRegion->GetRegionBlock();
            int s_lr = BCRegion->GetFaceLeftOrRightIndex();
            int t_lr = BCRegion->GetFaceLeftOrRightIndexOfTargetBlock();
            int s_nd = BCRegion->GetFaceDirection();	
            int t_nd = BCRegion->GetFaceDirectionOfTargetBlock();
            ofile << "   The source region block is: "<< nbs << endl;	
            ofile << "   The boundary face left or right index of source region block is: "<< s_lr << endl;
            ofile << "   The boundary face left or right index of target region block is: "<< t_lr << endl;
            ofile << "   The boundary face direction of source region block is: "<< s_nd << endl;
            ofile << "   The boundary face direction of target region block is: "<< t_nd << endl;

        }
    }

    //! Close the output file and exit PHengLEI.
    ParallelCloseFile(ofile);

    PHMPI::FreeBlockData();
    TIME_SPACE::WriStatTime();
    PH_Finalize();
    return 0;
}