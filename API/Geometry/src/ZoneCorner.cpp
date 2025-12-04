#include "ZoneCorner.h"
#include "GridType.h"
#include "Geo_Grid.h"
#include "Geo_StructGrid.h"
#include "Geo_UnstructGrid.h"
#include "TK_Log.h"
#include "Task.h"
#include "Glb_Dimension.h"
using namespace std;

namespace PHSPACE
{

CornerGridIndex::CornerGridIndex()
{

}

CornerGridIndex::~CornerGridIndex()
{

}

void CornerGridIndex::SetCornerGirdIndex(Grid *grid, int &zoneCornerIndex, set<int> &zoneCornerNeighbor)
{
    int gridType = grid->Type();
    if (gridType == PHSPACE::UNSTRUCTGRID)
    {
        UnstructGrid *unstructGrid = UnstructGridCast(grid);
        this->SetCornerGirdIndex(unstructGrid, zoneCornerIndex, zoneCornerNeighbor);
    }
    else if (gridType == PHSPACE::STRUCTGRID)
    {
        StructGrid *structGrid = StructGridCast(grid);
        this->SetCornerGirdIndex(structGrid, zoneCornerIndex, zoneCornerNeighbor);
    }
    else if (gridType == MIXGRID)
    {
        //! Mix grid.
        ostringstream oss;
        oss << " Error: gridType MIXGRID " << endl;
        TK_Exit::ExceptionExit(oss);
    }
}

void CornerGridIndex::SetCornerGirdIndex(StructGrid *structGrid, int &zoneCornerIndex, set<int> &zoneCornerNeighbor)
{
    //! Try not to use StructBC search here, 
    //! because for zones with no neighbors, 
    //! the default critical is 0, causing repetitive errors.
    //! For example, zone 0 is surrounded by a boundary of 1, 
    //! but its StructBCSet has 4 structBCs, 
    //! so the three additional structBCs will default to neighbor 0.
    //! Get boundary condition.
    //! Here we put the old identification method as follows:

    //! Store the neighbor of corner zone information.
    map<int, vector<int> > repeatPart;
    ostringstream ossCorner;
    ossCorner << " --- Check SetCornerGirdIndex On Struct grid CornerIndex --- " << "\n";
    ossCorner << " zone Center id : " << structGrid->GetZoneID() << "\n";
    ossCorner << " zone Corner id : " << zoneCornerIndex << "\n";
    //! Loop for the neighbor of corner.
    for (set<int>::iterator iNeighbor = zoneCornerNeighbor.begin(); iNeighbor != zoneCornerNeighbor.end(); ++iNeighbor)
    {
        ossCorner << " neighbor of Corner " << *iNeighbor << "\n";
    }
    StructBCSet *structBCSet = structGrid->GetStructBCSet();
    int nBC = structBCSet->GetnBCRegion();
    for (uint_t iBC = 0; iBC < nBC; ++iBC)
    {
        //! The bc of center zone.
        StructBC *structBC = structBCSet->GetBCRegion(iBC);
        int bcType = structBC->GetBCType();
        if (bcType == PHENGLEI::INTERFACE)
        {
            //! The neighbor zone index of zone Center.
            int zoneNeighborIndex = structBC->GetTargetRegionBlock();
            ossCorner << "  ibc " << iBC  << " bcType " << bcType  << "\n";
            ossCorner << "   zoneNeighborIndex of center " << zoneNeighborIndex << "\n";
            int *s_lr3d = structBC->GetFaceDirectionIndex();
            ossCorner << "   FaceDirectionIndex : " << s_lr3d[0] << " " << s_lr3d[1] << " " << s_lr3d[2] << "\n";

            set<int>::iterator iterZoneCornerNeighbor = zoneCornerNeighbor.find(zoneNeighborIndex);
            if (iterZoneCornerNeighbor != zoneCornerNeighbor.end())
            {
                int *s_lr3d = structBC->GetFaceDirectionIndex();
                vector<int> s_lr3dVector;
                s_lr3dVector.push_back(s_lr3d[0]);
                s_lr3dVector.push_back(s_lr3d[1]);
                s_lr3dVector.push_back(s_lr3d[2]);

                bool check = repeatPart.insert(make_pair(zoneNeighborIndex, s_lr3dVector)).second;
                if (!check)
                {
                    ostringstream ossError;
                    ossError << "Error: insert repeatPart: " << zoneNeighborIndex << endl;
                    TK_Exit::ExceptionExit(ossError);
                }
            }
            else
            {
                continue;
            }
        }
    }

    //! The index of i,j,k
    //! These same to StructBC::s_lr3d[3].
    //! These is the corner id, eg [-1,-1,0] and [1,1,1]
    int s_lr3di, s_lr3dj, s_lr3dk;
    s_lr3di = 0;
    s_lr3dj = 0;
    s_lr3dk = 0;

    for (map<int, vector<int> >::iterator iterRepeat = repeatPart.begin(); iterRepeat != repeatPart.end(); ++iterRepeat)
    {
        int gridGlobalIndex = iterRepeat->first;
        //! s_lr3d && t_lr3d stores the composite information of both s_nd/t_nd && s_lr/t_lr.
        //! It is used to uniform the loop format only. The values are equal to the following:
        vector<int> s_lr3dVector = iterRepeat->second;
        s_lr3di += s_lr3dVector[0];
        s_lr3dj += s_lr3dVector[1];
        s_lr3dk += s_lr3dVector[2];
    }

    this->indexStr[0] = s_lr3di;
    this->indexStr[1] = s_lr3dj;
    this->indexStr[2] = s_lr3dk;

    ossCorner << "  indexStr : " << indexStr[0] << " " << indexStr[1] << " " << indexStr[2] << "\n" ;
    ossCorner << endl;

}

void CornerGridIndex::SetCornerGirdIndex(UnstructGrid *unstructGrid, int &zoneCornerIndex, set<int> &zoneCornerNeighbor)
{

}

int *CornerGridIndex::GetIndex(Grid *grid)
{
    int gridType = grid->Type();
    int *index = this->GetIndex(gridType);
    return index;
}

int *CornerGridIndex::GetIndex(int gridType)
{
    int *index;
    index = 0;
    if (gridType == PHSPACE::UNSTRUCTGRID)
    {
        index = this->indexUnstr;
    }
    else if (gridType == PHSPACE::STRUCTGRID)
    {
        index = this->indexStr;
    }
    else if (gridType == MIXGRID)
    {
        //! Mix grid.
        ostringstream oss;
        oss << " Error: gridType MIXGRID " << endl;
        TK_Exit::ExceptionExit(oss);
    }
    return index;
}

ZoneCornerIndex::ZoneCornerIndex()
{

}

ZoneCornerIndex::~ZoneCornerIndex()
{

}

void ZoneCornerIndex::SetZoneCornerIndex(int iCornerZoneIndex)
{
    this->zoneCornerIndex = iCornerZoneIndex;
}

void ZoneCornerIndex::AddZoneCornerNeighborIndex(set<int> indexOfCornerNeighbor)
{
    this->zoneCornerNeighbor = indexOfCornerNeighbor;
}

void ZoneCornerIndex::SetCornerGirdIndex(Grid *grid)
{
    this->cornerGridIndex.SetCornerGirdIndex(grid,this->zoneCornerIndex,this->zoneCornerNeighbor);
}

set<int> &ZoneCornerIndex::GetZoneCornerNeighbor()
{
    return this->zoneCornerNeighbor;
}

int *ZoneCornerIndex::GetCornerGridIndex(int gridType)
{
    return this->cornerGridIndex.GetIndex(gridType);
}

int *ZoneCornerIndex::GetCornerGridIndex(Grid *grid)
{
    return this->cornerGridIndex.GetIndex(grid);
}

int ZoneCornerIndex::GetCornerZoneIndex()
{
    return this->zoneCornerIndex;
}

ZoneCorner::ZoneCorner()
{
    this->globalZoneIndex = -1;
    this->gridType = -1;
}

ZoneCorner::~ZoneCorner()
{

}

void ZoneCorner::AddNeighborOfNeighbor(int iNeighborZoneIndex, set<int> indexOfNeighbor)
{
    bool check = neighborOfNeighbor.insert(make_pair(iNeighborZoneIndex, indexOfNeighbor)).second;
    if (!check)
    {
        ostringstream oss;
        oss << "Error: insert ZoneCorner map: " << "iNeighbor = " << iNeighborZoneIndex << endl;
        TK_Exit::ExceptionExit(oss);
    }
}

void ZoneCorner::AddCornerIndexForPoint(ZoneNeighbor *zoneNeighbor)
{
    //! To do later.
}

void ZoneCorner::CalcCornerIndex()
{
    int iCorner = 0;

    //! Loop for the 1st index of Neighbor zone.
    for (map<int, set<int> >::iterator iterSetNeighbor = neighborOfNeighbor.begin();
        iterSetNeighbor != neighborOfNeighbor.end();
        ++iterSetNeighbor)
    {
        //! Loop for the 2ed index of Neighbor zone.
        for (map<int, set<int> >::iterator iterSetNeighborFind = neighborOfNeighbor.begin();
            iterSetNeighborFind != neighborOfNeighbor.end();
            ++iterSetNeighborFind)
        {
            //! Loop for the 3rd index of Neighbor zone.
            for (map<int, set<int> >::iterator iterSetNeighborLast = neighborOfNeighbor.begin();
                iterSetNeighborLast != neighborOfNeighbor.end();
                ++iterSetNeighborLast)
            {
                //! The global 1st zone index of neighbor.
                int iNeighborZoneIndex = iterSetNeighbor->first;

                //! The global 2ed zone index of anthor neighbor.
                int iNeighborZoneIndexFind = iterSetNeighborFind->first;

                //! The global 3rd zone index of anthor neighbor.
                int iNeighborZoneIndexLast = iterSetNeighborLast->first;

                if (iNeighborZoneIndex == iNeighborZoneIndexFind
                    && iNeighborZoneIndex == iNeighborZoneIndexLast)
                {
                    continue;
                }

                for (set<int>::iterator iterNeighborNeighbor = iterSetNeighbor->second.begin();
                    iterNeighborNeighbor != iterSetNeighbor->second.end();
                    ++iterNeighborNeighbor)
                {
                    int iNeighborNeighborZoneIndex = *iterNeighborNeighbor;

                    if (iNeighborNeighborZoneIndex == globalZoneIndex)
                    {
                        continue;
                    }

                    for (set<int>::iterator iterNeighborNeighborFind = iterSetNeighborFind->second.begin();
                        iterNeighborNeighborFind != iterSetNeighborFind->second.end();
                        ++iterNeighborNeighborFind)
                    {
                        int iNeighborNeighborZoneIndexFind = *iterNeighborNeighborFind;

                        if (iNeighborNeighborZoneIndexFind == globalZoneIndex)
                        {
                            continue;
                        }

                        for (set<int>::iterator iterNeighborNeighborLast = iterSetNeighborLast->second.begin();
                            iterNeighborNeighborLast != iterSetNeighborLast->second.end();
                            ++iterNeighborNeighborLast)
                        {
                            int iNeighborNeighborZoneIndexLast = *iterNeighborNeighborLast;

                            if (iNeighborNeighborZoneIndexLast == globalZoneIndex)
                            {
                                continue;
                            }

                            if (iNeighborNeighborZoneIndex == iNeighborNeighborZoneIndexFind &&
                                iNeighborNeighborZoneIndex == iNeighborNeighborZoneIndexLast)
                            {
                                //! So first, we need to check to see if there are duplicate id.
                                bool ifExist = this->isExistCornerIndex(iNeighborNeighborZoneIndex);

                                if (!ifExist)
                                {

                                    ZoneCornerIndex zoneCornerIndex;
                                    zoneCornerIndex.SetZoneCornerIndex(iNeighborNeighborZoneIndex);

                                    set<int> zoneCornerNeighbor;
                                    zoneCornerNeighbor.insert(iNeighborZoneIndex);
                                    zoneCornerNeighbor.insert(iNeighborZoneIndexFind);
                                    zoneCornerNeighbor.insert(iNeighborZoneIndexLast);

                                    zoneCornerIndex.AddZoneCornerNeighborIndex(zoneCornerNeighbor);

                                    bool check = cornerInfo.insert(make_pair(iNeighborNeighborZoneIndex, zoneCornerIndex)).second;
                                    iCorner += 1;
                                    if (!check)
                                    {
                                        ostringstream oss;
                                        oss << "Error: insert ZoneCorner Index: " << "iCorner = " << iCorner << endl;
                                        TK_Exit::ExceptionExit(oss);
                                    }
                                }
                                else
                                {
                                    continue;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void ZoneCorner::SetZoneCornerIndexGridType(int type)
{
    this->gridType = type;
}

void ZoneCorner::SetZoneCornerIndexByGrid(Grid *grid)
{
    //! Loop for each corner on current zone.
    for (map<int, ZoneCornerIndex>::iterator iterCorner = cornerInfo.begin(); iterCorner != cornerInfo.end();++iterCorner)
    {
        iterCorner->second.SetCornerGirdIndex(grid);
    }
}

void ZoneCorner::SetGlobalZoneIndex(int iZone)
{
    this->globalZoneIndex = iZone;
}

map<int, set<int> > &ZoneCorner::GetSetIndexNeighborOfNeighbor()
{
    return neighborOfNeighbor;
}

map<int, ZoneCornerIndex> &ZoneCorner::GetCornerInfo()
{
    return cornerInfo;
}

bool ZoneCorner::isExistCornerIndex(int iCornerIndex)
{
    bool  ifExist = false;

    map<int, ZoneCornerIndex>::iterator iterCornerIndex = this->cornerInfo.find(iCornerIndex);

    if (iterCornerIndex != this->cornerInfo.end())
    {
        ifExist = true;
    }
    else
    {
        ifExist = false;
    }

    return ifExist;
}

int ZoneCorner::GetNumOfCorner()
{
    return this->cornerInfo.size();
}

int ZoneCorner::GetZoneIndexOfCorner(int iCorner)
{
    int zoneIndex = 0;
    int countCorner = -1;
    for (map<int, ZoneCornerIndex>::iterator iter = this->cornerInfo.begin(); iter != cornerInfo.end();
        ++iter)
    {
        countCorner += 1;
        if (countCorner == iCorner)
        {
            zoneIndex = iter->first;
        }
    }

    return zoneIndex;
}

map<int, ZoneCornerIndex >::iterator ZoneCorner::GetCornerIterBegin()
{
    map<int, ZoneCornerIndex >::iterator iter = this->cornerInfo.begin();
    return iter;
}

map<int, ZoneCornerIndex >::iterator ZoneCorner::GetCornerIterEnd()
{
    map<int, ZoneCornerIndex >::iterator iter = this->cornerInfo.end();
    return iter;
}

void MainTaskBlockingCorner(ActionKey *actkey)
{
#ifdef USE_LagrangianParticle
    bool useParSolver = GlobalDataBase::IsExist("iParticleModel", 1, 1);
    if (useParSolver)
    {
        int iParticleModel = GlobalDataBase::GetIntParaFromDB("iParticleModel");
        if (iParticleModel == 1)
        {
            using namespace PHMPI;

            int nZones = GetNumberofGlobalZones();

            for (int iZone = 0; iZone < nZones; ++iZone)
            {
                //! One thing to watch out for here is 
                //! when there is no Corner.
                //! So the ghostand corner layers of 
                //! StructX must be assigned.
                //! XCC and XFN, etc., are Ghost with only the first layer.
                //! For the case without Corner, 
                //! the geometry values of the second level Ghostand corner 
                //! are not needed
                //! (but the variable values of the second level Ghost are needed).
                ZoneCorner *zoneCorner = zoneConnectivity->GetZoneCorner(iZone);
                int nCorner = zoneCorner->GetNumOfCorner();

                using namespace FILLCORNER;
                FillGridCorner(GetGrid(iZone, 0));

                for (map<int, ZoneCornerIndex >::iterator iterCorner = zoneCorner->GetCornerIterBegin();
                    iterCorner != zoneCorner->GetCornerIterEnd();
                    ++iterCorner)
                {
                    PH_Interface(actkey, iZone, iterCorner->second.GetCornerZoneIndex());
                }
            }
        }
    }
#endif
}

namespace FILLCORNER
{


void FillGridCorner(Grid *grid)
{
    int gridType = grid->Type();

    if (gridType == PHSPACE::STRUCTGRID)
    {
        StructGrid *structgrid = StructGridCast(grid);

        FillGridCorner(structgrid);
    }
    else if (gridType == PHSPACE::UNSTRUCTGRID)
    {
        UnstructGrid *unstructgrid = UnstructGridCast(grid);
        FillGridCorner(unstructgrid);
    }
    else if (gridType == PHSPACE::MIXGRID)
    {
        //! Mixing
        TK_Exit::UnexpectedVarValue(" gridType = ", PHSPACE::MIXGRID);
    }
}

void FillGridCorner(StructGrid *grid)
{
    //! variable       ||  ghost  || corner 
    //! vol(-1:ni+1)         2         2
    //! xcc(-1:ni+1)         1         0
    //! xx(-1:ni+1)          2         2
    //! xfn(-1:ni+1)         1         0
    //! area,xfv seem to xfn

    RDouble4D &xfn = *(grid->GetFaceNormalX());
    RDouble4D &yfn = *(grid->GetFaceNormalY());
    RDouble4D &zfn = *(grid->GetFaceNormalZ());
    RDouble4D &area = *(grid->GetFaceArea());
    RDouble3D &vol = *(grid->GetCellVolume());
    RDouble4D &xfv = *(grid->GetFaceVectorX());
    RDouble4D &yfv = *(grid->GetFaceVectorY());
    RDouble4D &zfv = *(grid->GetFaceVectorZ());

    //! =================================================================
    //! ===                    Fill Cell center.                      ===
    //! =================================================================
    RDouble3D &xcc = *(grid->GetCellCenterX());
    RDouble3D &ycc = *(grid->GetCellCenterY());
    RDouble3D &zcc = *(grid->GetCellCenterZ());

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 0);

    //! for corner.
    vector<int> xCornerIndex;
    vector<int> yCornerIndex;
    vector<int> zCornerIndex;

    xCornerIndex.push_back(ist - 1);
    xCornerIndex.push_back(ied + 1);
    xCornerIndex.push_back(ist - 2);
    xCornerIndex.push_back(ied + 2);

    yCornerIndex.push_back(jst - 1);
    yCornerIndex.push_back(jed + 1);
    yCornerIndex.push_back(jst - 2);
    yCornerIndex.push_back(jed + 2);

    if (GetDim() == two)
    {
        zCornerIndex.push_back(1);
    }
    else
    {
        zCornerIndex.push_back(kst - 2);
        zCornerIndex.push_back(kst - 1);
        zCornerIndex.push_back(ked + 1);
        zCornerIndex.push_back(ked + 2);
    }

    for (vector<int>::iterator iIter = xCornerIndex.begin(); iIter != xCornerIndex.end(); ++iIter)
    {
        for (vector<int>::iterator jIter = yCornerIndex.begin(); jIter != yCornerIndex.end(); ++jIter)
        {
            for (vector<int>::iterator kIter = zCornerIndex.begin(); kIter != zCornerIndex.end(); ++kIter)
            {
                int i, j, k;
                i = *iIter;
                j = *jIter;
                k = *kIter;
                grid->CenterCoor(i, j, k, xcc(i, j, k), ycc(i, j, k), zcc(i, j, k));
            }
        }
    }
    xCornerIndex.clear();
    yCornerIndex.clear();
    zCornerIndex.clear();

    //! for ghost.
    for (int k = kst; k <= ked; ++k)
    {
        for (int j = jst; j <= jed; ++j)
        {
            int i;

            i = ist - 1;
            grid->CenterCoor(i, j, k, xcc(i, j, k), ycc(i, j, k), zcc(i, j, k));

            i = ist - 2;
            grid->CenterCoor(i, j, k, xcc(i, j, k), ycc(i, j, k), zcc(i, j, k));

            i = ied + 1;
            grid->CenterCoor(i, j, k, xcc(i, j, k), ycc(i, j, k), zcc(i, j, k));

            i = ied + 2;
            grid->CenterCoor(i, j, k, xcc(i, j, k), ycc(i, j, k), zcc(i, j, k));
        }
    }

    for (int k = kst; k <= ked; ++k)
    {
        for (int i = ist; i <= ied; ++i)
        {
            int j;

            j = jst - 1;
            grid->CenterCoor(i, j, k, xcc(i, j, k), ycc(i, j, k), zcc(i, j, k));

            j = jst - 2;
            grid->CenterCoor(i, j, k, xcc(i, j, k), ycc(i, j, k), zcc(i, j, k));

            j = jed + 1;
            grid->CenterCoor(i, j, k, xcc(i, j, k), ycc(i, j, k), zcc(i, j, k));

            j = jed + 2;
            grid->CenterCoor(i, j, k, xcc(i, j, k), ycc(i, j, k), zcc(i, j, k));
        }
    }

    if (GetDim() == THREE_D)
    {
        for (int j = jst; j <= jed; ++j)
        {
            for (int i = ist; i <= ied; ++i)
            {
                int k;

                k = kst - 1;
                grid->CenterCoor(i, j, k, xcc(i, j, k), ycc(i, j, k), zcc(i, j, k));

                k = kst - 2;
                grid->CenterCoor(i, j, k, xcc(i, j, k), ycc(i, j, k), zcc(i, j, k));

                k = ked + 1;
                grid->CenterCoor(i, j, k, xcc(i, j, k), ycc(i, j, k), zcc(i, j, k));

                k = ked + 2;
                grid->CenterCoor(i, j, k, xcc(i, j, k), ycc(i, j, k), zcc(i, j, k));
            }
        }
    }

    //! =================================================================
    //! ===                    Fill face vector.                      ===
    //! =================================================================
    //! Fill xfn,xfv and area,corner 1
    int isurf;
    int i1, j1, k1;
    int i2, j2, k2;

    if (GetDim() == TWO_D)
    {
        zfn = 0.0;
        zfv = 0.0;

        //! For corner 1.
        isurf = 1;
        i1 = 0; i2 = 0;
        j1 = 0; j2 = 1;
        k1 = 1; k2 = 1;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 1;
        i1 = 0; i2 = 0;
        j1 = nj; j2 = nj - 1;
        k1 = 1; k2 = 1;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 1;
        i1 = ni; i2 = ni;
        j1 = 0; j2 = 1;
        k1 = 1; k2 = 1;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 1;
        i1 = ni; i2 = ni;
        j1 = nj; j2 = nj - 1;
        k1 = 1; k2 = 1;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 2;
        i1 = 0; i2 = 1;
        j1 = 0; j2 = 0;
        k1 = 1; k2 = 1;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 2;
        i1 = 0; i2 = 1;
        j1 = nj; j2 = nj;
        k1 = 1; k2 = 1;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 2;
        i1 = ni; i2 = ni - 1;
        j1 = 0; j2 = 0;
        k1 = 1; k2 = 1;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 2;
        i1 = ni; i2 = ni - 1;
        j1 = nj; j2 = nj;
        k1 = 1; k2 = 1;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        //! for ghost 2.
        for (int j = jst-1; j <= jed + 1; ++j)
        {
            isurf = 2;
            i1 = -1; i2 = 0;
            k1 = 1; k2 = 1;

            xfn(i1, j, k1, isurf) = xfn(i2, j, k2, isurf);
            yfn(i1, j, k1, isurf) = yfn(i2, j, k2, isurf);
            zfn(i1, j, k1, isurf) = zfn(i2, j, k2, isurf);
            area(i1, j, k1, isurf) = area(i2, j, k2, isurf);
            xfv(i1, j, k1, isurf) = xfv(i2, j, k2, isurf);
            yfv(i1, j, k1, isurf) = yfv(i2, j, k2, isurf);
            zfv(i1, j, k1, isurf) = zfv(i2, j, k2, isurf);

            isurf = 1;

            xfn(i1, j, k1, isurf) = xfn(i2, j, k2, isurf);
            yfn(i1, j, k1, isurf) = yfn(i2, j, k2, isurf);
            zfn(i1, j, k1, isurf) = zfn(i2, j, k2, isurf);
            area(i1, j, k1, isurf) = area(i2, j, k2, isurf);
            xfv(i1, j, k1, isurf) = xfv(i2, j, k2, isurf);
            yfv(i1, j, k1, isurf) = yfv(i2, j, k2, isurf);
            zfv(i1, j, k1, isurf) = zfv(i2, j, k2, isurf);

        }

        for (int i = ist - 1; i <= ied + 1; ++i)
        {
            isurf = 1;
            j1 = -1; j2 = 0;
            k1 = 1; k2 = 1;

            xfn(i, j1, k1, isurf) = xfn(i, j2, k2, isurf);
            yfn(i, j1, k1, isurf) = yfn(i, j2, k2, isurf);
            zfn(i, j1, k1, isurf) = zfn(i, j2, k2, isurf);
            area(i, j1, k1, isurf) = area(i, j2, k2, isurf);
            xfv(i, j1, k1, isurf) = xfv(i, j2, k2, isurf);
            yfv(i, j1, k1, isurf) = yfv(i, j2, k2, isurf);
            zfv(i, j1, k1, isurf) = zfv(i, j2, k2, isurf);

            isurf = 2;

            xfn(i, j1, k1, isurf) = xfn(i, j2, k2, isurf);
            yfn(i, j1, k1, isurf) = yfn(i, j2, k2, isurf);
            zfn(i, j1, k1, isurf) = zfn(i, j2, k2, isurf);
            area(i, j1, k1, isurf) = area(i, j2, k2, isurf);
            xfv(i, j1, k1, isurf) = xfv(i, j2, k2, isurf);
            yfv(i, j1, k1, isurf) = yfv(i, j2, k2, isurf);
            zfv(i, j1, k1, isurf) = zfv(i, j2, k2, isurf);
        }

        //! Second corner
        isurf = 1;
        i1 = -1; i2 = -1;
        j1 = -1; j2 = 0;
        k1 = 1; k2 = 1;
        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 1;
        i1 = ni + 1; i2 = ni + 1;
        j1 = -1; j2 = 0;
        k1 = 1; k2 = 1;
        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 1;
        i1 = -1; i2 = -1;
        j1 = nj + 1; j2 = nj;
        k1 = 1; k2 = 1;
        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 1;
        i1 = ni + 1; i2 = ni + 1;
        j1 = nj + 1; j2 = nj;
        k1 = 1; k2 = 1;
        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 2;
        i1 = -1; i2 = 0;
        j1 = -1; j2 = -1;
        k1 = 1; k2 = 1;
        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 2;
        i1 = -1; i2 = 0;
        j1 = nj + 1; j2 = nj + 1;
        k1 = 1; k2 = 1;
        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 2;
        i1 = ni + 1; i2 = ni;
        j1 = -1; j2 = -1;
        k1 = 1; k2 = 1;
        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 2;
        i1 = ni + 1; i2 = ni;
        j1 = nj + 1; j2 = nj + 1;
        k1 = 1; k2 = 1;
        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

    }
    else if (GetDim() == THREE_D)
    {
        isurf = 1;
        i1 = 0; i2 = 0;
        j1 = 0; j2 = 1;
        k1 = 0; k2 = 0;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 1;
        i1 = 0; i2 = 0;
        j1 = 0; j2 = 1;
        k1 = nk; k2 = nk;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 1;
        i1 = ni; i2 = ni;
        j1 = 0; j2 = 1;
        k1 = 0; k2 = 0;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 1;
        i1 = ni; i2 = ni;
        j1 = 0; j2 = 1;
        k1 = nk; k2 = nk;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 1;
        i1 = 0; i2 = 0;
        j1 = nj; j2 = nj - 1;
        k1 = 0; k2 = 0;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 1;
        i1 = 0; i2 = 0;
        j1 = nj; j2 = nj - 1;
        k1 = nk; k2 = nk;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 1;
        i1 = ni; i2 = ni;
        j1 = nj; j2 = nj - 1;
        k1 = 0; k2 = 0;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 1;
        i1 = ni; i2 = ni;
        j1 = nj; j2 = nj - 1;
        k1 = nk; k2 = nk;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 1;
        i1 = 0; i2 = 0;
        j1 = 0; j2 = 0;
        k1 = 0; k2 = 1;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 1;
        i1 = 0; i2 = 0;
        j1 = nj; j2 = nj;
        k1 = 0; k2 = 1;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 1;
        i1 = ni; i2 = ni;
        j1 = 0; j2 = 0;
        k1 = 0; k2 = 1;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 1;
        i1 = ni; i2 = ni;
        j1 = nj; j2 = nj;
        k1 = 0; k2 = 1;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 1;
        i1 = 0; i2 = 0;
        j1 = 0; j2 = 0;
        k1 = nk; k2 = nk - 1;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 1;
        i1 = 0; i2 = 0;
        j1 = nj; j2 = nj;
        k1 = nk; k2 = nk - 1;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 1;
        i1 = ni; i2 = ni;
        j1 = 0; j2 = 0;
        k1 = nk; k2 = nk - 1;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 1;
        i1 = ni; i2 = ni;
        j1 = nj; j2 = nj;
        k1 = nk; k2 = nk - 1;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 2;
        i1 = 0; i2 = 1;
        j1 = 0; j2 = 0;
        k1 = 0; k2 = 0;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 2;
        i1 = 0; i2 = 1;
        j1 = 0; j2 = 0;
        k1 = nk; k2 = nk;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 2;
        i1 = 0; i2 = 1;
        j1 = nj; j2 = nj;
        k1 = 0; k2 = 0;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 2;
        i1 = 0; i2 = 1;
        j1 = nj; j2 = nj;
        k1 = nk; k2 = nk;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 2;
        i1 = ni; i2 = ni - 1;
        j1 = 0; j2 = 0;
        k1 = 0; k2 = 0;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 2;
        i1 = ni; i2 = ni - 1;
        j1 = 0; j2 = 0;
        k1 = nk; k2 = nk;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 2;
        i1 = ni; i2 = ni - 1;
        j1 = nj; j2 = nj;
        k1 = nk; k2 = nk;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 2;
        i1 = 0; i2 = 0;
        j1 = 0; j2 = 0;
        k1 = 0; k2 = 1;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 2;
        i1 = 0; i2 = 0;
        j1 = nj; j2 = nj;
        k1 = 0; k2 = 1;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 2;
        i1 = ni; i2 = ni;
        j1 = 0; j2 = 0;
        k1 = 0; k2 = 1;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 2;
        i1 = ni; i2 = ni;
        j1 = nj; j2 = nj;
        k1 = 0; k2 = 1;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 2;
        i1 = 0; i2 = 0;
        j1 = 0; j2 = 0;
        k1 = nk; k2 = nk - 1;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 2;
        i1 = 0; i2 = 0;
        j1 = nj; j2 = nj;
        k1 = nk; k2 = nk - 1;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 2;
        i1 = ni; i2 = ni;
        j1 = nj; j2 = nj;
        k1 = nk; k2 = nk - 1;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 3;
        i1 = 0; i2 = 1;
        j1 = 0; j2 = 0;
        k1 = 0; k2 = 0;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 3;
        i1 = 0; i2 = 1;
        j1 = 0; j2 = 0;
        k1 = nk; k2 = nk;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 3;
        i1 = 0; i2 = 1;
        j1 = nj; j2 = nj;
        k1 = 0; k2 = 0;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 3;
        i1 = 0; i2 = 1;
        j1 = nj; j2 = nj;
        k1 = nk; k2 = nk;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 3;
        i1 = ni; i2 = ni - 1;
        j1 = 0; j2 = 0;
        k1 = 0; k2 = 0;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 3;
        i1 = ni; i2 = ni - 1;
        j1 = 0; j2 = 0;
        k1 = nk; k2 = nk;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 3;
        i1 = ni; i2 = ni - 1;
        j1 = nj; j2 = nj;
        k1 = 0; k2 = 0;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 3;
        i1 = ni; i2 = ni - 1;
        j1 = nj; j2 = nj;
        k1 = nk; k2 = nk;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 3;
        i1 = 0; i2 = 0;
        j1 = 0; j2 = 1;
        k1 = 0; k2 = 0;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 3;
        i1 = 0; i2 = 0;
        j1 = 0; j2 = 1;
        k1 = nk; k2 = nk;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 3;
        i1 = ni; i2 = ni;
        j1 = 0; j2 = 1;
        k1 = nk; k2 = nk;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 3;
        i1 = 0; i2 = 0;
        j1 = nj; j2 = nj - 1;
        k1 = 0; k2 = 0;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 3;
        i1 = 0; i2 = 0;
        j1 = nj; j2 = nj - 1;
        k1 = nk; k2 = nk;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 3;
        i1 = ni; i2 = ni;
        j1 = nj; j2 = nj - 1;
        k1 = nk; k2 = nk;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        //! Second Ghost
        //! for -1 0
        int i, j, k;
        isurf = 1;
        i1 = -1; i2 = 0;
        j1 = 0; j2 = 0;
        k1 = 1; k2 = 1;
        for (k = 0; k <= nk; ++k)
        {
            for (j = 0; j <= nj; ++j)
            {
                isurf = 1;
                i1 = -1; i2 = 0;
                k1 = 1; k2 = 1;

                xfn(i1, j, k, isurf) = xfn(i2, j, k, isurf);
                yfn(i1, j, k, isurf) = yfn(i2, j, k, isurf);
                zfn(i1, j, k, isurf) = zfn(i2, j, k, isurf);
                area(i1, j, k, isurf) = area(i2, j, k, isurf);
                xfv(i1, j, k, isurf) = xfv(i2, j, k, isurf);
                yfv(i1, j, k, isurf) = yfv(i2, j, k, isurf);
                zfv(i1, j, k, isurf) = zfv(i2, j, k, isurf);

                isurf = 2;

                xfn(i1, j, k, isurf) = xfn(i2, j, k, isurf);
                yfn(i1, j, k, isurf) = yfn(i2, j, k, isurf);
                zfn(i1, j, k, isurf) = zfn(i2, j, k, isurf);
                area(i1, j, k, isurf) = area(i2, j, k, isurf);
                xfv(i1, j, k, isurf) = xfv(i2, j, k, isurf);
                yfv(i1, j, k, isurf) = yfv(i2, j, k, isurf);
                zfv(i1, j, k, isurf) = zfv(i2, j, k, isurf);

                isurf = 3;

                xfn(i1, j, k, isurf) = xfn(i2, j, k, isurf);
                yfn(i1, j, k, isurf) = yfn(i2, j, k, isurf);
                zfn(i1, j, k, isurf) = zfn(i2, j, k, isurf);
                area(i1, j, k, isurf) = area(i2, j, k, isurf);
                xfv(i1, j, k, isurf) = xfv(i2, j, k, isurf);
                yfv(i1, j, k, isurf) = yfv(i2, j, k, isurf);
                zfv(i1, j, k, isurf) = zfv(i2, j, k, isurf);

            }

            for (i = 0; i <= ni; ++i)
            {
                j1 = -1; j2 = 0;
                isurf = 1;

                xfn(i, j1, k, isurf) = xfn(i, j2, k, isurf);
                yfn(i, j1, k, isurf) = yfn(i, j2, k, isurf);
                zfn(i, j1, k, isurf) = zfn(i, j2, k, isurf);
                area(i, j1, k, isurf) = area(i, j2, k, isurf);
                xfv(i, j1, k, isurf) = xfv(i, j2, k, isurf);
                yfv(i, j1, k, isurf) = yfv(i, j2, k, isurf);
                zfv(i, j1, k, isurf) = zfv(i, j2, k, isurf);

                isurf = 2;

                xfn(i, j1, k, isurf) = xfn(i, j2, k, isurf);
                yfn(i, j1, k, isurf) = yfn(i, j2, k, isurf);
                zfn(i, j1, k, isurf) = zfn(i, j2, k, isurf);
                area(i, j1, k, isurf) = area(i, j2, k, isurf);
                xfv(i, j1, k, isurf) = xfv(i, j2, k, isurf);
                yfv(i, j1, k, isurf) = yfv(i, j2, k, isurf);
                zfv(i, j1, k, isurf) = zfv(i, j2, k, isurf);

                isurf = 3;

                xfn(i, j1, k, isurf) = xfn(i, j2, k, isurf);
                yfn(i, j1, k, isurf) = yfn(i, j2, k, isurf);
                zfn(i, j1, k, isurf) = zfn(i, j2, k, isurf);
                area(i, j1, k, isurf) = area(i, j2, k, isurf);
                xfv(i, j1, k, isurf) = xfv(i, j2, k, isurf);
                yfv(i, j1, k, isurf) = yfv(i, j2, k, isurf);
                zfv(i, j1, k, isurf) = zfv(i, j2, k, isurf);
            }
        }

        for (i = 0; i <= ni; ++i)
        {
            for (j = 0; j <= nj; ++j)
            {
                k1 = -1; k2 = 0;
                isurf = 1;

                xfn(i, j, k1, isurf) = xfn(i, j, k2, isurf);
                yfn(i, j, k1, isurf) = yfn(i, j, k2, isurf);
                zfn(i, j, k1, isurf) = zfn(i, j, k2, isurf);
                area(i, j, k1, isurf) = area(i, j, k2, isurf);
                xfv(i, j, k1, isurf) = xfv(i, j, k2, isurf);
                yfv(i, j, k1, isurf) = yfv(i, j, k2, isurf);
                zfv(i, j, k1, isurf) = zfv(i, j, k2, isurf);

                isurf = 2;

                xfn(i, j, k1, isurf) = xfn(i, j, k2, isurf);
                yfn(i, j, k1, isurf) = yfn(i, j, k2, isurf);
                zfn(i, j, k1, isurf) = zfn(i, j, k2, isurf);
                area(i, j, k1, isurf) = area(i, j, k2, isurf);
                xfv(i, j, k1, isurf) = xfv(i, j, k2, isurf);
                yfv(i, j, k1, isurf) = yfv(i, j, k2, isurf);
                zfv(i, j, k1, isurf) = zfv(i, j, k2, isurf);

                isurf = 3;

                xfn(i, j, k1, isurf) = xfn(i, j, k2, isurf);
                yfn(i, j, k1, isurf) = yfn(i, j, k2, isurf);
                zfn(i, j, k1, isurf) = zfn(i, j, k2, isurf);
                area(i, j, k1, isurf) = area(i, j, k2, isurf);
                xfv(i, j, k1, isurf) = xfv(i, j, k2, isurf);
                yfv(i, j, k1, isurf) = yfv(i, j, k2, isurf);
                zfv(i, j, k1, isurf) = zfv(i, j, k2, isurf);
            }
        }

        isurf = 1;
        i1 = -1; i2 = -1;
        j1 = -1; j2 = 0;
        k1 = -1; k2 = -1;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 1;
        i1 = -1; i2 = -1;
        j1 = -1; j2 = 0;
        k1 = nk + 1; k2 = nk + 1;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 1;
        i1 = ni + 1; i2 = ni + 1;
        j1 = -1; j2 = 0;
        k1 = -1; k2 = -1;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);


        isurf = 1;
        i1 = ni + 1; i2 = ni + 1;
        j1 = -1; j2 = 0;
        k1 = nk + 1; k2 = nk + 1;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 1;
        i1 = -1; i2 = -1;
        j1 = nj + 1; j2 = nj;
        k1 = -1; k2 = -1;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 1;
        i1 = -1; i2 = -1;
        j1 = nj + 1; j2 = nj;
        k1 = nk + 1; k2 = nk + 1;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 1;
        i1 = ni + 1; i2 = ni + 1;
        j1 = nj + 1; j2 = nj;
        k1 = -1; k2 = -1;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 1;
        i1 = ni + 1; i2 = ni + 1;
        j1 = nj + 1; j2 = nj;
        k1 = nk + 1; k2 = nk + 1;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 1;
        i1 = -1; i2 = -1;
        j1 = -1; j2 = -1;
        k1 = -1; k2 = 0;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 1;
        i1 = -1; i2 = -1;
        j1 = nj + 1; j2 = nj + 1;
        k1 = -1; k2 = 0;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 1;
        i1 = ni + 1; i2 = ni + 1;
        j1 = -1; j2 = -1;
        k1 = -1; k2 = 0;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 1;
        i1 = ni + 1; i2 = ni + 1;
        j1 = nj + 1; j2 = nj + 1;
        k1 = -1; k2 = 0;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 1;
        i1 = -1; i2 = -1;
        j1 = -1; j2 = -1;
        k1 = nk + 1; k2 = nk;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 1;
        i1 = -1; i2 = -1;
        j1 = nj + 1; j2 = nj + 1;
        k1 = nk + 1; k2 = nk;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 1;
        i1 = ni + 1; i2 = ni + 1;
        j1 = -1; j2 = -1;
        k1 = nk + 1; k2 = nk;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 1;
        i1 = ni + 1; i2 = ni + 1;
        j1 = nj + 1; j2 = nj + 1;
        k1 = nk + 1; k2 = nk;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 2;
        i1 = -1; i2 = 0;
        j1 = -1; j2 = -1;
        k1 = -1; k2 = -1;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 2;
        i1 = -1; i2 = 0;
        j1 = -1; j2 = -1;
        k1 = nk + 1; k2 = nk + 1;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 2;
        i1 = -1; i2 = 0;
        j1 = nj + 1; j2 = nj + 1;
        k1 = -1; k2 = -1;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 2;
        i1 = -1; i2 = 0;
        j1 = nj + 1; j2 = nj + 1;
        k1 = nk + 1; k2 = nk + 1;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 2;
        i1 = ni + 1; i2 = ni;
        j1 = -1; j2 = -1;
        k1 = -1; k2 = -1;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 2;
        i1 = ni + 1; i2 = ni;
        j1 = -1; j2 = -1;
        k1 = nk + 1; k2 = nk + 1;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 2;
        i1 = ni + 1; i2 = ni;
        j1 = nj + 1; j2 = nj + 1;
        k1 = nk + 1; k2 = nk + 1;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 2;
        i1 = -1; i2 = -1;
        j1 = -1; j2 = -1;
        k1 = -1; k2 = 0;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 2;
        i1 = -1; i2 = -1;
        j1 = nj + 1; j2 = nj + 1;
        k1 = -1; k2 = 0;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 2;
        i1 = ni + 1; i2 = ni + 1;
        j1 = -1; j2 = -1;
        k1 = -1; k2 = 0;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 2;
        i1 = ni + 1; i2 = ni + 1;
        j1 = nj + 1; j2 = nj + 1;
        k1 = -1; k2 = 0;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 2;
        i1 = -1; i2 = -1;
        j1 = -1; j2 = -1;
        k1 = nk + 1; k2 = nk;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 2;
        i1 = -1; i2 = -1;
        j1 = nj + 1; j2 = nj + 1;
        k1 = nk + 1; k2 = nk;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 2;
        i1 = ni + 1; i2 = ni + 1;
        j1 = nj + 1; j2 = nj + 1;
        k1 = nk + 1; k2 = nk;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 3;
        i1 = -1; i2 = 0;
        j1 = -1; j2 = -1;
        k1 = -1; k2 = -1;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 3;
        i1 = -1; i2 = 0;
        j1 = -1; j2 = -1;
        k1 = nk + 1; k2 = nk + 1;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 3;
        i1 = -1; i2 = 0;
        j1 = nj + 1; j2 = nj + 1;
        k1 = -1; k2 = -1;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 3;
        i1 = -1; i2 = 0;
        j1 = nj + 1; j2 = nj + 1;
        k1 = nk + 1; k2 = nk + 1;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 3;
        i1 = ni + 1; i2 = ni;
        j1 = -1; j2 = -1;
        k1 = -1; k2 = -1;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 3;
        i1 = ni + 1; i2 = ni;
        j1 = -1; j2 = -1;
        k1 = nk + 1; k2 = nk + 1;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 3;
        i1 = ni + 1; i2 = ni;
        j1 = nj + 1; j2 = nj + 1;
        k1 = -1; k2 = -1;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 3;
        i1 = ni + 1; i2 = ni;
        j1 = nj + 1; j2 = nj + 1;
        k1 = nk + 1; k2 = nk + 1;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 3;
        i1 = -1; i2 = -1;
        j1 = -1; j2 = 0;
        k1 = -1; k2 = -1;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 3;
        i1 = -1; i2 = -1;
        j1 = -1; j2 = 0;
        k1 = nk + 1; k2 = nk + 1;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 3;
        i1 = ni + 1; i2 = ni + 1;
        j1 = -1; j2 = 0;
        k1 = nk + 1; k2 = nk + 1;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 3;
        i1 = -1; i2 = -1;
        j1 = nj + 1; j2 = nj;
        k1 = -1; k2 = -1;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 3;
        i1 = -1; i2 = -1;
        j1 = nj + 1; j2 = nj;
        k1 = nk + 1; k2 = nk + 1;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);

        isurf = 3;
        i1 = ni + 1; i2 = ni + 1;
        j1 = nj + 1; j2 = nj;
        k1 = nk + 1; k2 = nk + 1;

        xfn(i1, j1, k1, isurf) = xfn(i2, j2, k2, isurf);
        yfn(i1, j1, k1, isurf) = yfn(i2, j2, k2, isurf);
        zfn(i1, j1, k1, isurf) = zfn(i2, j2, k2, isurf);
        area(i1, j1, k1, isurf) = area(i2, j2, k2, isurf);
        xfv(i1, j1, k1, isurf) = xfv(i2, j2, k2, isurf);
        yfv(i1, j1, k1, isurf) = yfv(i2, j2, k2, isurf);
        zfv(i1, j1, k1, isurf) = zfv(i2, j2, k2, isurf);
    }

}

void FillGridCorner(UnstructGrid *grid)
{

}

}

}