#include "Pre_GridRefine.h"
#include "Glb_Dimension.h"
#include "Pre_GridConversion.h"
#include "Zone.h"

using namespace std;

namespace PHSPACE
{
LIB_EXPORT Pre_GridRefine::Pre_GridRefine(const string &from_gfile)
{
    this->gridFileName = from_gfile;
}

LIB_EXPORT Pre_GridRefine::~Pre_GridRefine()
{
    for (int izone = 0; izone < numberOfBlocks; ++ izone)
    {
        delete NewGrids[izone];
    }
    delete [] NewGrids;
}

LIB_EXPORT void Pre_GridRefine::Run()
{
    using namespace PHMPI;

    string out_gfile = "grid.fts";
    GlobalDataBase::GetData("out_gfile", &out_gfile, PHSTRING, 1);
    GlobalDataBase::UpdateData("gridfile", &gridFileName, PHSTRING, 1);

    cout << "  Task type: " << "RefineGrid" << endl;
    cout << " From_gfile: " << gridFileName << endl;
    cout << "  Out_gfile: " << out_gfile    << endl << endl;

    Region *region = new Region();

    region->ReadGrid();

    numberOfBlocks = PHMPI::GetNumberofGlobalZones();
    for (int iZone = 0; iZone < numberOfBlocks; ++ iZone)
    {
        Zone *zone = region->GetZone(iZone);
        if (!zone) continue;

        Grid *grid = zone->GetGeometry()->GetGrid();
        OrdinaryGrid.push_back(grid);
    }

    cout << endl << "End reading grid ... " << endl << endl;

    RefineGrid();
    RefineBCRegion();

    cout << endl << "Refine grid finished ... " << endl << endl;

    ConstructGlobalInterfaceLink(NewGrids, numberOfBlocks);
    DumpGrid(out_gfile, NewGrids, numberOfBlocks);

    delete region;
}

void Pre_GridRefine::RefineGrid()
{
    NewGrids = new Grid * [numberOfBlocks];

    for (int izone = 0; izone < numberOfBlocks; izone ++)
    {
        StructGrid *grid = StructGridCast(OrdinaryGrid[izone]);
        RefineGridStr(izone, grid);
    }
}

void Pre_GridRefine::RefineGridStr(int izone, StructGrid *grid_in)
{
    int ii, jj, kk;

    GridID *index_new = new GridID(izone);
    Grid *grid_new = CreateGridGeneral(STRUCTGRID, index_new, 0, GetDim());
    NewGrids[izone] = grid_new;
    StructGrid *strgrid_new = StructGridCast(grid_new);

    int ni = grid_in->GetNI();
    int nj = grid_in->GetNJ();
    int nk = grid_in->GetNK();

    strgrid_new->SetNI((ni-1)*2+1);
    strgrid_new->SetNJ((nj-1)*2+1);
    strgrid_new->SetNK((nk-1)*2+1);
    strgrid_new->SetBasicDimension();

    int nTotalNode = strgrid_new->GetNTotalNode();
    RDouble *nx = new RDouble [nTotalNode];
    RDouble *ny = new RDouble [nTotalNode];
    RDouble *nz = new RDouble [nTotalNode];
    strgrid_new->SetX(nx);
    strgrid_new->SetY(ny);
    strgrid_new->SetZ(nz);
    strgrid_new->SetArrayLayout();

    int ist, ied, jst, jed, kst, ked;
    grid_in->GetNodeIterationIndex(ist, ied, jst, jed, kst, ked);

    RDouble3D &xx = * grid_in->GetStructX();
    RDouble3D &yy = * grid_in->GetStructY();
    RDouble3D &zz = * grid_in->GetStructZ();

    RDouble3D &nxx = * strgrid_new->GetStructX();
    RDouble3D &nyy = * strgrid_new->GetStructY();
    RDouble3D &nzz = * strgrid_new->GetStructZ();

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                ii = (i-1) * 2 + 1;
                jj = (j-1) * 2 + 1;
                kk = (k-1) * 2 + 1;
                nxx(ii, jj, kk) = xx(i, j, k);
                nyy(ii, jj, kk) = yy(i, j, k);
                nzz(ii, jj, kk) = zz(i, j, k);
            }
        }
    }

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied-1; ++ i)
            {
                ii = (i-1) * 2 + 2;
                jj = (j-1) * 2 + 1;
                kk = (k-1) * 2 + 1;
                nxx(ii, jj, kk) = 0.5 * (nxx(ii-1, jj, kk) + nxx(ii+1, jj, kk));
                nyy(ii, jj, kk) = 0.5 * (nyy(ii-1, jj, kk) + nyy(ii+1, jj, kk));
                nzz(ii, jj, kk) = 0.5 * (nzz(ii-1, jj, kk) + nzz(ii+1, jj, kk));
            }
        }
    }

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed-1; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                ii = (i-1) * 2 + 1;
                jj = (j-1) * 2 + 2;
                kk = (k-1) * 2 + 1;
                nxx(ii, jj, kk) = 0.5 * (nxx(ii, jj-1, kk) + nxx(ii, jj+1, kk));
                nyy(ii, jj, kk) = 0.5 * (nyy(ii, jj-1, kk) + nyy(ii, jj+1, kk));
                nzz(ii, jj, kk) = 0.5 * (nzz(ii, jj-1, kk) + nzz(ii, jj+1, kk));
            }
        }
    }

    for (int k = kst; k <= ked-1; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                ii = (i-1) * 2 + 1;
                jj = (j-1) * 2 + 1;
                kk = (k-1) * 2 + 2;
                nxx(ii, jj, kk) = 0.5 * (nxx(ii, jj, kk-1) + nxx(ii, jj, kk+1));
                nyy(ii, jj, kk) = 0.5 * (nyy(ii, jj, kk-1) + nyy(ii, jj, kk+1));
                nzz(ii, jj, kk) = 0.5 * (nzz(ii, jj, kk-1) + nzz(ii, jj, kk+1));
            }
        }
    }

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed-1; ++ j)
        {
            for (int i = ist; i <= ied-1; ++ i)
            {
                ii = (i-1) * 2 + 2;
                jj = (j-1) * 2 + 2;
                kk = (k-1) * 2 + 1;
                nxx(ii, jj, kk) = 0.25 * (nxx(ii-1, jj, kk) + nxx(ii+1, jj, kk) + nxx(ii, jj-1, kk) + nxx(ii, jj+1, kk));
                nyy(ii, jj, kk) = 0.25 * (nyy(ii-1, jj, kk) + nyy(ii+1, jj, kk) + nyy(ii, jj-1, kk) + nyy(ii, jj+1, kk));
                nzz(ii, jj, kk) = 0.25 * (nzz(ii-1, jj, kk) + nzz(ii+1, jj, kk) + nzz(ii, jj-1, kk) + nzz(ii, jj+1, kk));
            }
        }
    }

    for (int k = kst; k <= ked-1; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied-1; ++ i)
            {
                ii = (i-1) * 2 + 2;
                jj = (j-1) * 2 + 1;
                kk = (k-1) * 2 + 2;
                nxx(ii, jj, kk) = 0.25 * (nxx(ii-1, jj, kk) + nxx(ii+1, jj, kk) + nxx(ii, jj, kk-1) + nxx(ii, jj, kk+1));
                nyy(ii, jj, kk) = 0.25 * (nyy(ii-1, jj, kk) + nyy(ii+1, jj, kk) + nyy(ii, jj, kk-1) + nyy(ii, jj, kk+1));
                nzz(ii, jj, kk) = 0.25 * (nzz(ii-1, jj, kk) + nzz(ii+1, jj, kk) + nzz(ii, jj, kk-1) + nzz(ii, jj, kk+1));
            }
        }
    }

    for (int k = kst; k <= ked-1; ++ k)
    {
        for (int j = jst; j <= jed-1; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                ii = (i-1) * 2 + 1;
                jj = (j-1) * 2 + 2;
                kk = (k-1) * 2 + 2;
                nxx(ii, jj, kk) = 0.25 * (nxx(ii, jj-1, kk) + nxx(ii, jj+1, kk) + nxx(ii, jj, kk-1) + nxx(ii, jj, kk+1));
                nyy(ii, jj, kk) = 0.25 * (nyy(ii, jj-1, kk) + nyy(ii, jj+1, kk) + nyy(ii, jj, kk-1) + nyy(ii, jj, kk+1));
                nzz(ii, jj, kk) = 0.25 * (nzz(ii, jj-1, kk) + nzz(ii, jj+1, kk) + nzz(ii, jj, kk-1) + nzz(ii, jj, kk+1));
            }
        }
    }

    for (int k = kst; k <= ked-1; ++ k)
    {
        for (int j = jst; j <= jed-1; ++ j)
        {
            for (int i = ist; i <= ied-1; ++ i)
            {
                ii = (i-1) * 2 + 2;
                jj = (j-1) * 2 + 2;
                kk = (k-1) * 2 + 2;
                nxx(ii, jj, kk) = 0.125 * (nxx(ii-1, jj-1, kk-1) + nxx(ii+1, jj-1, kk-1) + nxx(ii-1, jj+1, kk-1) + nxx(ii+1, jj+1, kk-1)
                                         + nxx(ii-1, jj-1, kk+1) + nxx(ii+1, jj-1, kk+1) + nxx(ii-1, jj+1, kk+1) + nxx(ii+1, jj+1, kk+1));
                nyy(ii, jj, kk) = 0.125 * (nyy(ii-1, jj-1, kk-1) + nyy(ii+1, jj-1, kk-1) + nyy(ii-1, jj+1, kk-1) + nyy(ii+1, jj+1, kk-1)
                                         + nyy(ii-1, jj-1, kk+1) + nyy(ii+1, jj-1, kk+1) + nyy(ii-1, jj+1, kk+1) + nyy(ii+1, jj+1, kk+1));
                nzz(ii, jj, kk) = 0.125 * (nzz(ii-1, jj-1, kk-1) + nzz(ii+1, jj-1, kk-1) + nzz(ii-1, jj+1, kk-1) + nzz(ii+1, jj+1, kk-1)
                                         + nzz(ii-1, jj-1, kk+1) + nzz(ii+1, jj-1, kk+1) + nzz(ii-1, jj+1, kk+1) + nzz(ii+1, jj+1, kk+1));
            }
        }
    }

    strgrid_new->ComputeMinMaxBox();
}

void Pre_GridRefine::RefineBCRegion()
{
    for (int izone = 0; izone < numberOfBlocks; izone ++)
    {
        StructGrid *grid = StructGridCast(OrdinaryGrid[izone]);
        StructGrid *strgrid_new = StructGridCast(NewGrids[izone]);

        strgrid_new->CopyCompositeBCRegion(grid->GetCompositeBCRegion());
        
        CompositeStructBC *rhs = strgrid_new->GetCompositeBCRegion();

        int numberOfBoundaryFaces = rhs->GetnBCRegion();

        for (int ibcregion = 0; ibcregion < numberOfBoundaryFaces; ibcregion ++)
        {
            StructBC *bcregion = rhs->GetBCRegion(ibcregion);

            int bctype = bcregion->GetBCType();

            int *s_st = bcregion->GetStartPoint();
            int *s_ed = bcregion->GetEndPoint();
            int *t_st = bcregion->GetTargetStart();
            int *t_ed = bcregion->GetTargetEnd();
            int *s_lr3d = bcregion->GetFaceDirectionIndex();
            int *t_lr3d = bcregion->GetFaceDirectionIndexOfTargetBlock();

            for (int m = 0; m < GetDim(); ++ m)
            {
                if (s_st[m] == s_ed[m])
                {
                    if (s_st[m] > 1 || s_lr3d[m] == 1)
                    //if ( s_st[m] != 1 )
                    {
                        s_st[m] += 1;
                        s_ed[m] += 1;
                    }
                }
                else
                {
                    s_ed[m] += 1;
                }
            }

            for (int m = 0; m < GetDim(); ++ m)
            {
                s_st[m] = (s_st[m]-1) * 2 + 1;
                s_ed[m] = (s_ed[m]-1) * 2 + 1;
            }

            if (bctype == PHENGLEI::INTERFACE)
            {
                for (int m = 0; m < GetDim(); ++ m)
                {
                    if (t_st[m] == t_ed[m])
                    {
                        if (t_st[m] > 1 || t_lr3d[m] == 1)
                        //if ( t_st[m] != 1 )
                        {
                            t_st[m] += 1;
                            t_ed[m] += 1;
                        }
                    }
                    else
                    {
                        t_ed[m] += 1;
                    }
                }

                for (int m = 0; m < GetDim(); ++ m)
                {
                    t_st[m] = (t_st[m]-1) * 2 + 1;
                    t_ed[m] = (t_ed[m]-1) * 2 + 1;
                }
            }
        }

        strgrid_new->ProcessBCInfo();
    }
}

}
