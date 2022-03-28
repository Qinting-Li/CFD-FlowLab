#include "Mesh_Repair.h"
#include "GridType.h"
#include "TK_Log.h"
#include "Geometry.h"

namespace PHSPACE
{
LIB_EXPORT Mesh_Repair::Mesh_Repair(UnstructGrid * originalGridIn)
{
    this->originalGrid = originalGridIn;

    this->repairedGrid = 0;
}

LIB_EXPORT Mesh_Repair::~Mesh_Repair()
{

}

LIB_EXPORT void Mesh_Repair::Repair()
{
    int nTotalFace = originalGrid->GetNTotalFace();
    int nBoundFace = originalGrid->GetNBoundFace();
    int nTotalCell = originalGrid->GetNTotalCell();

    int * left_cell_of_face  = originalGrid->GetLeftCellOfFace();
    int * right_cell_of_face = originalGrid->GetRightCellOfFace();

    int level = originalGrid->GetLevel();

    bool * isBadFace = new bool [nTotalFace];
    bool * isBadCell = new bool [nTotalCell];

    PickOutBadSkewnessFaces(originalGrid, isBadFace, isBadCell);

    bool * unMerged = new bool [nTotalCell];
    for (int icell = 0; icell < nTotalCell; ++ icell)
    {
        unMerged[icell] = true;
    }

    // now merging
    int *cell2coarsegridcell = originalGrid->GetCell2CoarseGridCell();
    if(!cell2coarsegridcell)
    {
        cell2coarsegridcell = new int[ nTotalCell ];
        originalGrid->SetCell2CoarseGridCell(cell2coarsegridcell);
    }

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        cell2coarsegridcell[iCell] = - 1;
    }

    int cnTCell = 0;
    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        if (isBadFace[iFace])
        {
            int le = left_cell_of_face[ iFace ];
            int re = right_cell_of_face[ iFace ];
            if (unMerged[le] && unMerged[re])
            {
                unMerged[le] = false;
                unMerged[re] = false;

                cell2coarsegridcell[le] = cnTCell;
                cell2coarsegridcell[re] = cnTCell ++;
            }
        }
    }

    // some cells may not have been merged, merge them now
    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        int le = left_cell_of_face [ iFace ];
        int re = right_cell_of_face[ iFace ];

        if (cell2coarsegridcell[ le ] == - 1) 
        {
            cell2coarsegridcell[ le ] = cnTCell ++;
        }

        if (cell2coarsegridcell[ re ] == - 1) 
        {
            cell2coarsegridcell[ re ] = cnTCell ++;
        }
    }

    delete [] unMerged;

    repairedGrid = new UnstructGrid();
    repairedGrid->InitGrid(new GridID(*originalGrid->GetGridID()), level + 1, originalGrid->GetDim(), UNSTRUCTGRID);

    repairedGrid->SetNTotalCell(cnTCell);
    int nTotalNode = originalGrid->GetNTotalNode();
    repairedGrid->SetNTotalNode(nTotalNode);
    RDouble *x, *y, *z;
    RDouble * xx = originalGrid->GetX();
    RDouble * yy = originalGrid->GetY();
    RDouble * zz = originalGrid->GetZ();

    x = new RDouble[ nTotalNode ];
    y = new RDouble[ nTotalNode ];
    z = new RDouble[ nTotalNode ];

    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        x[iNode] = xx[iNode];
        y[iNode] = yy[iNode];
        z[iNode] = zz[iNode];
    }
    repairedGrid->SetX(x);
    repairedGrid->SetY(y);
    repairedGrid->SetZ(z);

    delete [] isBadFace;
    delete [] isBadCell;
}

LIB_EXPORT Grid * RepairGrid(Grid * originalGridIn)
{

    int gridType = originalGridIn->Type();
    if (gridType == STRUCTGRID)
    {
        return 0;
    }
    else
    {
        //VisualizationMesh2D(this, "./grid/original.plt");
        UnstructGrid * grid = (static_cast<UnstructGrid*>(originalGridIn));
        int nTotalCell = grid->GetNTotalCell();
        int *cell2coarsegridcell = grid->GetCell2CoarseGridCell();
        if(!cell2coarsegridcell)
        {
            cell2coarsegridcell = new int[ nTotalCell ];
            grid->SetCell2CoarseGridCell(cell2coarsegridcell);
        }

        UnstructGrid * originalGrid = grid;
        UnstructGrid * repairedGrid = originalGrid;
        int repairTimes = 0;
        while(PickOutBadSkewnessFaces(originalGrid) && repairTimes <= 1000)
        {
            ++ repairTimes;

            PrintToWindow(" repair time: ", repairTimes, "...\n");

            Mesh_Repair * meshRepaire = new Mesh_Repair(originalGrid);
            meshRepaire->Repair();
            UnstructGrid * repairedGridTemp = meshRepaire->GetRepairedGrid();
            delete meshRepaire;

            CoarseGridConnectUNMergeFace(originalGrid, repairedGridTemp);
            repairedGridTemp->ComputeMetrics();

            originalGrid->SetCoarseGrid(repairedGridTemp);
            repairedGridTemp->SetFineGrid(originalGrid);
            repairedGridTemp->SetLevel(originalGrid->GetLevel()+1);
            repairedGridTemp->SetBCRecord(originalGrid->GetBCRecord());

            //! re-connect the cell to coarse cell information.
            int *org_cell2coarsegridcell = originalGrid->GetCell2CoarseGridCell();
            if(originalGrid != grid)
            {
                for (int iCell = 0; iCell < nTotalCell; ++ iCell)
                {
                    int orgCellID    = cell2coarsegridcell[iCell];
                    int cgridCellID  = org_cell2coarsegridcell[orgCellID];
                    cell2coarsegridcell[iCell] = cgridCellID;
                }
            }

            UnstructGrid * tempGrid = originalGrid;

            originalGrid = repairedGridTemp;
            repairedGrid = repairedGridTemp;

            if(repairTimes > 1)
            {
                delete tempGrid;
            }
        }

        //VisualizationMesh2D(repairedGrid, "./grid/repairedGrid.plt");

        return repairedGrid;
    }

}

LIB_EXPORT void RepairGrid(Zone * zone, ActionKey * actkey)
{

    PHGeometry *geom = zone->GetGeometry();
    Grid * originalGrid = geom->GetGrid(0);
    originalGrid->ComputeMetrics();
    Grid * repairGrid = RepairGrid(originalGrid);

    if(repairGrid != originalGrid)
    {
        zone->SetInterface(originalGrid, repairGrid);
    }

    DataContainer * cdata = actkey->GetData();

    cdata->moveto_begin();
    repairGrid->Encode(cdata, 0);
}

int PickOutBadSkewnessFaces(UnstructGrid * grid, bool * isBadFace, bool * isBadCell)
{
    const RDouble PI = 3.141592653589793;

    RDouble angTorelance  = 0.001;
    RDouble dotp1, x1, y1, z1, dis1, dotp2, x2, y2, z2, dis2, angle1, angle2, angle;
    RDouble area_scale = TINY;
    RDouble minAng = PHSPACE::LARGE;

    int le, re;
    int * left_cell_of_face = grid->GetLeftCellOfFace();
    int * right_cell_of_face = grid->GetRightCellOfFace();
    RDouble * xcc = grid->GetCellCenterX();
    RDouble * ycc = grid->GetCellCenterY();
    RDouble * zcc = grid->GetCellCenterZ();

    RDouble * xfc = grid->GetFaceCenterX();
    RDouble * yfc = grid->GetFaceCenterY();
    RDouble * zfc = grid->GetFaceCenterZ();

    RDouble * nxs = grid->GetFaceNormalX();
    RDouble * nys = grid->GetFaceNormalY();
    RDouble * nzs = grid->GetFaceNormalZ();
    RDouble * ns  = grid->GetFaceArea();

    int nTotalFace = grid->GetNTotalFace();
    int nBoundFace = grid->GetNBoundFace();

    if(isBadFace)
    {
        for (int iFace = 0; iFace < nTotalFace; ++ iFace)
        {
            isBadFace[iFace] = false;
        }
    }

    if(isBadCell)
    {
        int nTotalCell = grid->GetNTotalCell();
        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            isBadCell[iCell] = false;
        }
    }

    int nBadFace = 0;
    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        if (ns[iFace] <= area_scale) continue;
        le = left_cell_of_face [iFace];
        re = right_cell_of_face[iFace];

        x1 = xfc[iFace] - xcc[le];
        y1 = yfc[iFace] - ycc[le];
        z1 = zfc[iFace] - zcc[le];
        dis1  = sqrt(x1 * x1 + y1 * y1 + z1 * z1);
        dotp1 = (nxs[iFace] * x1 + nys[iFace] * y1 + nzs[iFace] * z1) / (dis1 + TINY);
        if (dotp1 >  1.0) dotp1 =   1.0;
        if (dotp1 < -1.0) dotp1 = - 1.0;
        angle1 = asin(dotp1) * 180.0 / PI;
        if(angle1 <= angTorelance)
        {
            if(isBadCell)
            {
                isBadCell[le] = true;
            }
        }

        x2 = xcc[re] - xfc[iFace];
        y2 = ycc[re] - yfc[iFace];
        z2 = zcc[re] - zfc[iFace];
        dis2  = sqrt(x2 * x2 + y2 * y2 + z2 * z2);
        dotp2 = (nxs[iFace] * x2 + nys[iFace] * y2 + nzs[iFace] * z2) / (dis2 + TINY);      
        if (dotp2 >  1.0) dotp2 =   1.0;
        if (dotp2 < -1.0) dotp2 = - 1.0;
        angle2 = asin(dotp2) * 180.0 / PI;
        if(angle2 <= angTorelance)
        {
            if(isBadCell)
            {
                isBadCell[re] = true;
            }
        }

        angle  = MIN(asin(dotp1), asin(dotp2)) * 180.0 / PI;
        minAng = MIN(minAng, angle);
        if (angle <= angTorelance)
        {
            ++ nBadFace;
            if(isBadFace)
            {
                isBadFace[iFace] = true;
            }
        }
    }

    if(nBadFace > 0)
    {
        WriteLogFile("    Local minimum skewness angle is ", minAng, true);
    }

    return nBadFace;
}



}