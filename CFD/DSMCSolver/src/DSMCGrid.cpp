#include "DSMCGrid.h"
#include "Geo_UnstructBC.h"
#include "TK_Exit.h"
#include "TK_Log.h"
#include <cmath>
using namespace std;

namespace PHSPACE
{
DSMCUnstructGrid::DSMCUnstructGrid()
{
    nWallFace = 0; 
    nInletFace = 0; 
    nOutletFace = 0; 
    nSymmetryFace = 0; 
    nInnerFace = 0;
    nMolecule = 0;

    wallFace2BoundaryFace = NULL;
    inletFace2BoundaryFace = NULL;
    outletFace2BoundaryFace = NULL;
    symmetryFace2BoundaryFace = NULL;
    interFace2BoundaryFace = NULL;
    cellWeight = NULL;
}

DSMCUnstructGrid::~DSMCUnstructGrid()
{
    delete [] cellWeight;   cellWeight = NULL;
    delete [] wallFace2BoundaryFace; wallFace2BoundaryFace = NULL;
    delete [] inletFace2BoundaryFace;
    delete [] outletFace2BoundaryFace;
    delete [] symmetryFace2BoundaryFace;
    delete [] interFace2BoundaryFace;
}

void DSMCUnstructGrid::BoundaryCount(UnstructGrid * grid)
{
    using namespace PHMPI;
    using namespace GRIDGEN_SPACE;

    int nBoundFace = grid->GetNBoundFace();
    int nIFace = grid->GetNIFace();
    //cout<<nTotalFace<<" "<<nBoundFace<<" "<<nIFace<<endl;

    UnstructBCSet **bcr = grid->GetBCRecord();

    for(int iface = 0; iface < nBoundFace; ++ iface)
    {
        int bctype = bcr[iface]->GetKey();

        switch (bctype)
        {
            case SOLID_SURFACE :
                nWallFace += 1;
                break;
            case INFLOW :
                nInletFace += 1;
                break;
            case OUTFLOW :
                nOutletFace += 1;
                break;
            case SYMMETRY :
                nSymmetryFace += 1;
                break;
            case INTERFACE :
                nInnerFace += 1;
                break;
            default:
                TK_Exit::UnexpectedVarValue("bctype", bctype);
                break;
        }
    }
    //cout<<nWallFace<<endl;
    //cout<<nInletFace<<endl;
    //cout<<nOutletFace<<endl;
    //cout<<nSymmetryFace<<endl;
    //cout<<nInnerFace<<endl;

    SetnWallFace(nWallFace);
    SetnInletFace(nInletFace);
    SetnOutletFace(nOutletFace);
    SetnSymmetryFace(nSymmetryFace);
    SetnInnerFace(nInnerFace);

    if(nInnerFace != nIFace) 
    {
        TK_Exit::ExceptionExit("nInnerFace != nIFace ");
    }

    /*int * wallFace2BoundaryFace = new int [nWallFace];
    int * inletFace2BoundaryFace = new int [nInletFace];
    int * outletFace2BoundaryFace = new int [nOutletFace];
    int * symmetryFace2BoundaryFace = new int [nSymmetryFace];
    int * interFace2BoundaryFace = new int [nInnerFace];*/

    //SetWallFace2BoundaryFace(wallFace2BoundaryFace);
    //SetInletFace2BoundaryFace(inletFace2BoundaryFace);
    //SetOutletFace2BoundaryFace(outletFace2BoundaryFace);
    //SetSymmetryFace2BoundaryFace(symmetryFace2BoundaryFace);
    //SetInterFace2BoundaryFace(interFace2BoundaryFace);

    //delete [] wallFace2BoundaryFace;
    //delete [] inletFace2BoundaryFace;
    //delete [] outletFace2BoundaryFace;
    //delete [] symmetryFace2BoundaryFace;
    //delete [] interFace2BoundaryFace;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// 函数头部说明：                                                            
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:   例如：CalCellWeight()                   
// AUTHOR(S):  Li Jin            
// PURPOSE:    计算网格单元权重因子
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void DSMCUnstructGrid::CalCellWeight(UnstructGrid * grid,RDouble referenceCellVolume)
{
    int nTotalCell = grid->GetNTotalCell();
    RDouble * cellVolume = grid->GetCellVolume();

    cellWeight = new RDouble [nTotalCell];

    for(int icell = 0; icell < nTotalCell; ++icell)
    {
        //if(cellVolume[icell] > referenceCellVolume)
        //{
            cellWeight[icell] = cellVolume[icell]/referenceCellVolume;
        //}
        //else
        //{
        //    cellWeight[icell] = 1.0;
        //}
    }
}

RDouble CalMaxCellVolume()
{
    using namespace PHMPI;
    RDouble maxCellVolume = 0.0;
    int nZones = GetNumberofGlobalZones();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int zoneID = PHMPI::GetLocalZoneIDToGlobalZoneID(iZone);

        UnstructGrid *grid = UnstructGridCast(GetGrid(zoneID));

        int nTotalCell = grid->GetNTotalCell();

        RDouble *vol = grid->GetCellVolume();

        for(int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            if(vol[iCell] > maxCellVolume) maxCellVolume = vol[iCell];
        }
    }
    
    ostringstream oss;
    oss << "maxCellVolume:" << maxCellVolume;
    WriteLogFile(oss.str());
    return maxCellVolume;
}

RDouble CalMinCellVolume()
{
    using namespace PHMPI;
    RDouble minCellVolume = 1.0e30;
    int nZones = GetNumberofGlobalZones();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int zoneID = PHMPI::GetLocalZoneIDToGlobalZoneID(iZone);

        UnstructGrid * grid = UnstructGridCast(GetGrid(zoneID));

        int nTotalCell = grid->GetNTotalCell();

        RDouble * vol = grid->GetCellVolume();

        for(int icell = 0; icell < nTotalCell; ++ icell)
        {
            if(vol[icell] < minCellVolume) minCellVolume = vol[icell];
        }
    }

    ostringstream oss;
    oss<<"minCellVolume:"<<minCellVolume;
    WriteLogFile(oss.str());
    return minCellVolume;
}

RDouble referenceCellVolume = 0.0;

RDouble GetReferenceCellVolume()
{
    return referenceCellVolume;
}

void CalReferenceCellVolume()
{
    RDouble maxCellVolume = CalMaxCellVolume();
    RDouble minCellVolume = CalMinCellVolume();
    referenceCellVolume = sqrt(maxCellVolume/minCellVolume)*minCellVolume;

    ostringstream oss;
    oss<<"referenceCellVolume:"<<referenceCellVolume;
    WriteLogFile(oss.str());

}

}

