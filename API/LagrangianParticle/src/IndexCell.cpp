#include "IndexCell.h"
#include "GridType.h"
#include "Glb_Dimension.h"
#include "TK_Exit.h"
#include "Geo_Grid.h"
#include "Geo_StructGrid.h"
#include "Geo_UnstructGrid.h"
#include "TK_Log.h"

using namespace std;

namespace PHSPACE
{
IndexCell::IndexCell()
{
    this->cellParticleID = new map<int, IndexParticleOfCell*>;
    this->nCell = -1;
}

IndexCell::~IndexCell()
{
    //! Destructors can cause many problems of deletion by mistake.
    //! We must be careful.

    if (0 == this)
    {
        return;
    }

    if (0 == cellParticleID)
    {
        return;
    }

    map<int, IndexParticleOfCell*>::iterator iter = cellParticleID->begin();
    while (iter != cellParticleID->end())
    {
        delete [] iter->second;
        iter->second = nullptr;
        cellParticleID->erase(iter++);
    }
    delete cellParticleID;
    cellParticleID = nullptr;
}

void IndexCell::InitNumOfCell(Grid *gridIn)
{
    int gridType = gridIn->Type();

    if (gridType == PHSPACE::STRUCTGRID)
    {
        StructGrid *grid = StructGridCast(gridIn);
        int nI, nJ, nK;
        grid->GetND(nI, nJ, nK);

        //! Note here, nCell include ghost cell.
        //int nLayers = GetNumberOfGhostCellLayers();
        this->nLayersBC = 1;

        int nDim = GetDim();
        if (TWO_D == nDim)
        {
            //! The cell all.
            this->nCell = (nI - 1 + 2* nLayersBC) * (nJ - 1 + 2* nLayersBC);
            //! The cell near the boundary.
            this->nCell -= (nI - 1 - 2 * nLayersBC) * (nJ - 1 - 2 * nLayersBC);
        }
        else
        {
            //! The cell all.
            this->nCell = (nI - 1 + 2 * nLayersBC) * (nJ - 1 + 2 * nLayersBC) * (nK - 1 + 2 * nLayersBC);
            //! The cell near the boundary.
            this->nCell -= (nI - 1 - 2 * nLayersBC) * (nJ - 1 - 2 * nLayersBC) * (nK - 1 - 2 * nLayersBC);
        }
    }
    else if (PHSPACE::UNSTRUCTGRID == gridType)
    {
        UnstructGrid *grid = UnstructGridCast(gridIn);
        this->nCell = grid->GetNTotalCell();
    }
    else if (PHSPACE::MIXGRID == gridType)
    {
        //! Mixing
        TK_Exit::UnexpectedVarValue(" gridType = ", PHSPACE::MIXGRID);
    }
}

int IndexCell::GetNumOfCell()
{
    return this->nCell;
}

void IndexCell::InitCellID(Grid *grid)
{
    int gridType = grid->Type();
    if (PHSPACE::STRUCTGRID == gridType)
    {
        StructGrid *structGrid = StructGridCast(grid);
        InitCellID(structGrid);
    }
    else if (PHSPACE::UNSTRUCTGRID == gridType)
    {
        UnstructGrid *unStructGrid = UnstructGridCast(grid);
        InitCellID(unStructGrid);
    }
    else if (PHSPACE::MIXGRID == gridType)
    {
        //! Mixing
        TK_Exit::UnexpectedVarValue(" gridType = ", PHSPACE::MIXGRID);
    }
}

void IndexCell::InitCellID(StructGrid *grid)
{
    bool isEmpty = this->cellParticleID->empty();

    if (isEmpty)
    {
        //! The map cellParticleID has not been init.
        //! Loop for each cell to init this map.
        //! Index for cell center.
        //! Index start and end of cell in current zone, layer is 2.
        int ist, ied;
        int jst, jed;
        int kst, ked;
        grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, this->nLayersBC);
        int nLayers = GetNumberOfGhostCellLayers();
        //! Get index IJK of particle cell position in current zone.
        int nI, nJ, nK;
        grid->GetND(nI, nJ, nK);
        nI -= 1;
        nJ -= 1;
        nK -= 1;
        if (TWO_D == GetDim())
        {
            nK = 1;
        }

        int count = 0;
        ostringstream ossLog;
        for (int k = kst; k <= ked; ++k)
        {
            for (int j = jst; j <= jed; ++j)
            {
                for (int i = ist; i <= ied; ++i)
                {
                    bool inI = i >= 1 + this->nLayersBC && i <= nI - this->nLayersBC;
                    bool inJ = j >= 1 + this->nLayersBC && j <= nJ - this->nLayersBC;
                    bool inK = k >= 1 + this->nLayersBC && k <= nK - this->nLayersBC;

                    bool inZone = (inI && inJ && inK);
                    if (TWO_D == GetDim())
                    {
                        inZone = inI && inJ;
                    }
                    if (!inZone)
                    {
                        IndexParticleOfCell *indexParticleOfCell = new IndexParticleOfCell();
                        int iCell;
                        GetIDCellOfParticleByIJKSturct(nI, nJ, nK, i, j,  k, nLayers, iCell);
                        //! Insert indexParticleOfCell and check.
                        bool check = this->cellParticleID->insert(make_pair(iCell, indexParticleOfCell)).second;

                        if (!check)
                        {
                            ostringstream ossError;
                            ossError << "Error: insert init cell: " << "Cell ID = " << iCell << endl;
                            TK_Exit::ExceptionExit(ossError);
                        }
                        ossLog << " iCell " << iCell << " i " << i << " j " << j << " k " << k << "\n";
                        count++;
                    }
                }
            }
        }

        ossLog << " sum count " << count << " nCell " << this->nCell;
        ossLog << endl;
        //PrintToWindow(ossLog); 

        if (count != this->nCell)
        {
            ostringstream ossError;
            ossError << "Error: the number of cell: " << " count " << count;
            ossError << " nCell " << this->nCell << endl;
            TK_Exit::ExceptionExit(ossError);
        }
    }
    else
    {
        //! The map cellParticleID had been init.
        TK_Exit::UnexpectedVarValue(" Un expect map size = ", this->cellParticleID->size());
    }
}

void IndexCell::InitCellID(UnstructGrid *grid)
{
    bool isEmpty = this->cellParticleID->empty();

    if (isEmpty)
    {
        //! The map cellParticleID has not been init.
        //! Loop for each cell to init this map.
    }
    else
    {
        //! The map cellParticleID had been init.
        TK_Exit::UnexpectedVarValue(" Un expect map size = ", this->cellParticleID->size());
    }
}

int IndexCell::InitIndexOfParticleInCell(Grid *grid, int idParticle,OnePointVariable *onePointVariable, bool &inZone )
{
    //! Get coordinate for each particle.
    int IDCellOfParticle = SearchParticleCell(grid, onePointVariable->GetOnePointCoordinate(),inZone);
    if (inZone)
    {
        this->AddParticleID(IDCellOfParticle, idParticle);
        return IDCellOfParticle;
    }

    return IDCellOfParticle;
}

void IndexCell::InitIndexOfParticleInCell(Grid *grid, map<int, OnePointVariable*> *particlePointVariable)
{
    bool isEmpty = particlePointVariable->empty();
    if (isEmpty)
    {
        cout << "The Init particle point Group is empty" << endl;
        return;
    }

    int nParticle = static_cast<int>(particlePointVariable->size());
    map<int, OnePointVariable*>::iterator iter;
    for (iter = particlePointVariable->begin(); iter != particlePointVariable->end(); ++iter)
    {
        //! Get coordinate for each particle.
        SimplePointer<RDouble> *onePointCoordinate = iter->second->GetOnePointCoordinate();
        //! Seach current particle.
        bool inZone;
        int IDCellOfParticle = SearchParticleCell(grid, onePointCoordinate, inZone);

        if (inZone)
        {
            //! Add particle ID.
            this->AddParticleID(IDCellOfParticle, iter->first);
        }
        else
        {
            //! This particle is not in current zone
            continue;
        }
    }
}

int IndexCell::SearchParticleCell(Grid *grid, SimplePointer<RDouble> *onePointCoordinate,bool &inZone)
{
    //! This Part is refer to ProbeData::SearchRealCell.
    int level = 0;
    int gridType = grid->Type();
    int IDCellOfParticle = -1;
    if (PHSPACE::STRUCTGRID == gridType)
    {
        StructGrid *structGrid = StructGridCast(grid);
        IDCellOfParticle = SearchParticleCell(structGrid, onePointCoordinate,inZone);
    }
    else if (PHSPACE::UNSTRUCTGRID == gridType)
    {
        UnstructGrid *unStructGrid = UnstructGridCast(grid);
        IDCellOfParticle = SearchParticleCell(unStructGrid, onePointCoordinate,inZone);
    }
    else if (PHSPACE::MIXGRID == gridType)
    {
        //! Mixing
        TK_Exit::UnexpectedVarValue(" gridType = ", PHSPACE::MIXGRID);
    }
    return IDCellOfParticle;
}

int IndexCell::SearchParticleCell(StructGrid *grid, SimplePointer<RDouble> *onePointCoordinate, bool &inZone)
{
    //! Index start and end of cell in current zone.
    int ist, ied;
    int jst, jed;
    int kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    //! Index start and end of surface in .
    //! From 1 to nDim.
    int nDim = GetDim();
    int iSurfaceStart = 1;
    int iSurfaceEnd = nDim;

    //! Get cell center.
    //! Note here, GetCellCenterX() return RDouble3D*,
    //! So RDouble3D &xcc is a refer to RDouble3D*.
    RDouble3D &xcc = *(grid->GetCellCenterX());
    RDouble3D &ycc = *(grid->GetCellCenterY());
    RDouble3D &zcc = *(grid->GetCellCenterZ());

    //! Get face normal.
    RDouble4D &xfn = *(grid->GetFaceNormalX());
    RDouble4D &yfn = *(grid->GetFaceNormalY());
    RDouble4D &zfn = *(grid->GetFaceNormalZ());

    //! The local index of surface in one cell.
    //! The index start from 0 to nDim-1;
    int iSurLocal, jSurLocal, kSurLocal;
    //! The local index for parallel surface of local cell.
    int iSurLocalParallel, jSurLocalParallel, kSurLocalParallel;

    //! The bool for particle in current cell.
    //! 1 -- the particle is in current cell.
    //! 0 -- else.
    bool isParticleInCell = false;

    //! To judge if particle is on surface.
    RDouble small = 1.0e-8;

    int indexForCell = -1;
    int indexForParticleInCell = -1;

    int iP, jP, kP;

    //! First, loop for index of cell in struct grid.
    for (int k = kst; k <= ked; ++k)
    {
        for (int j = jst; j <= jed; ++j)
        {
            for (int i = ist; i <= ied; ++i)
            {
                //! Add cell index.
                ++indexForCell;

                //! Second. loop for index of surface.
                for (int iSurface = iSurfaceStart; iSurface <= iSurfaceEnd; ++iSurface)
                {
                    //! Get the normal vector of surface
                    RDouble xfn1 = xfn(i, j, k, iSurface);
                    RDouble yfn1 = yfn(i, j, k, iSurface);
                    RDouble zfn1 = zfn(i, j, k, iSurface);

                    //! The center coordinate of surface.
                    RDouble xfc1, yfc1, zfc1;
                    //! Get the center coodinate of surface.
                    grid->FaceCoor(i, j, k, iSurface, xfc1, yfc1, zfc1);

                    //! Calculate the dot of Vector from particle to surface center..
                    RDouble particleToFaceX = (*onePointCoordinate)[0] - xfc1;
                    RDouble particleToFaceY = (*onePointCoordinate)[1] - yfc1;
                    RDouble particleToFaceZ = (*onePointCoordinate)[2] - zfc1;
                    RDouble dotParticleToSurface = particleToFaceX * xfn1
                                                 + particleToFaceY * yfn1
                                                 + particleToFaceZ * zfn1;
                    RDouble distanceParticleToFace = DISTANCE(particleToFaceX, particleToFaceY, particleToFaceZ);

                    //! Calculate the dote of Vector from cell center to surface center.
                    RDouble centerToFaceX = xcc(i, j, k) - xfc1;
                    RDouble centerToFaceY = ycc(i, j, k) - yfc1;
                    RDouble centerToFaceZ = zcc(i, j, k) - zfc1;
                    RDouble dotCenterToSurface = centerToFaceX * xfn1
                                               + centerToFaceY * yfn1
                                               + centerToFaceZ * zfn1;

                    //! Get the local index of parallel surface.
                    //! Get the local index of surface in current one cell.
                    GetNsurfIndex(iSurface, iSurLocal, jSurLocal, kSurLocal);
                    //! iSurLocalParallel is parallel surface.
                    iSurLocalParallel = i + iSurLocal;
                    jSurLocalParallel = j + jSurLocal;
                    kSurLocalParallel = k + kSurLocal;

                    RDouble xfn2 = xfn(iSurLocalParallel, jSurLocalParallel, kSurLocalParallel, iSurface);
                    RDouble yfn2 = yfn(iSurLocalParallel, jSurLocalParallel, kSurLocalParallel, iSurface);
                    RDouble zfn2 = zfn(iSurLocalParallel, jSurLocalParallel, kSurLocalParallel, iSurface);

                    //! The center coordinate of parallel surface.
                    RDouble xfc2, yfc2, zfc2;
                    grid->FaceCoor(iSurLocalParallel, jSurLocalParallel, kSurLocalParallel, iSurface, xfc2, yfc2, zfc2);

                    //! Calculate the dote of Vector from cell center to parallel surface center.
                    RDouble particleToParallelFaceX = (*onePointCoordinate)[0] - xfc2;
                    RDouble particleToParallelFaceY = (*onePointCoordinate)[1] - yfc2;
                    RDouble particleToParallelFaceZ = (*onePointCoordinate)[2] - zfc2;
                    RDouble dotParticleToParallelSurface = particleToParallelFaceX * xfn2
                                                         + particleToParallelFaceY * yfn2
                                                         + particleToParallelFaceZ * zfn2;
                    RDouble distanceParticleToParallelFace = DISTANCE(particleToParallelFaceX, particleToParallelFaceY, particleToParallelFaceZ);

                    //! Calculate the dote of Vector from cell center to surface center.
                    RDouble centerToParallelFaceX = xcc(i, j, k) - xfc2;
                    RDouble centerToParallelFaceY = ycc(i, j, k) - yfc2;
                    RDouble centerToParallelFaceZ = zcc(i, j, k) - zfc2;
                    RDouble dotCenterToParallelSurface = centerToParallelFaceX * xfn2
                        + centerToParallelFaceY * yfn2
                        + centerToParallelFaceZ * zfn2;

                    //! Judge if the particle in current cell.
                    //! Note here, break will end for-loop of iSurface.
                    //! So only isParticleInCell == true for all iSurface, the indexForCell can be add.
                    if (dotParticleToSurface * dotCenterToSurface > zero)
                    {
                        isParticleInCell = true;
                    }
                    else if (ABS(dotParticleToSurface / distanceParticleToFace) < small)
                    {
                        //! Particle is close to surface.
                        isParticleInCell = true;
                    }
                    else
                    {
                        isParticleInCell = false;
                        break;
                    }

                    //! For parallel surface.
                    if (dotParticleToParallelSurface * dotCenterToParallelSurface > zero)
                    {
                        isParticleInCell = true;
                    }
                    else if (ABS(dotParticleToParallelSurface / distanceParticleToParallelFace) < small)
                    {
                        //! Particle is close to parallel surface.
                        isParticleInCell = true;
                    }
                    else
                    {
                        isParticleInCell = false;
                        break;
                    }
                }
                
                //! After loop for surface.
                if (isParticleInCell)
                {
                    iP = i;
                    jP = j;
                    kP = k;
                    indexForParticleInCell = indexForCell;
                    break;
                }        
            }
            if (isParticleInCell) break;
        }
        if (isParticleInCell) break;
    }

    //! Get index IJK of particle position in current zone.
    int nI, nJ, nK;
    grid->GetND(nI, nJ, nK);
    nI -= 1;
    nJ -= 1;
    nK -= 1;
    if (TWO_D == GetDim())
    {
        nK = 1;
    }

    GetIDCellOfParticleByIJKSturct(nI, nJ, nK, iP, jP, kP, GetNumberOfGhostCellLayers(), indexForCell);

    indexForParticleInCell = indexForCell;

    if (isParticleInCell)
    {
        inZone = true;
    }
    else
    {
        inZone = false;
    }
    return indexForParticleInCell;
}

int IndexCell::SearchParticleCell(UnstructGrid *grid, SimplePointer<RDouble> *onePointCoordinate,bool &inZone)
{
    using namespace PHMPI;
    int gridtype = grid->Type();
    int zoneID = grid->GetGridID()->GetIndex();

    int nTotalCell = grid->GetNTotalCell();

    //! Get cell center.
    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    //! Get face normal.
    RDouble *xfn = grid->GetFaceNormalX();
    RDouble *yfn = grid->GetFaceNormalY();
    RDouble *zfn = grid->GetFaceNormalZ();

    //! Get the center coordinate of surface.
    //! Note here is different from struct grid.
    RDouble *xfc = grid->GetFaceCenterX();
    RDouble *yfc = grid->GetFaceCenterY();
    RDouble *zfc = grid->GetFaceCenterZ();

    //! The number of face in each cell.
    //! Note here , is a 2 dim array.
    int *faceNumberOfEachCell = grid->GetFaceNumberOfEachCell();

    //! The global ID of face in each cell
    int **cell2Face = grid->GetCell2Face();

    //! The bool for particle in current cell.
    //! 1 -- the particle is in current cell.
    //! 0 -- else.
    bool isParticleInCell = false;

    //! To judge if particle is on surface.
    RDouble small = 1.0e-8;

    int indexForCell = -1;
    int indexForParticleInCell = -1;
    for (int iCell = 0; iCell < nTotalCell; ++iCell)
    {
        ++indexForCell;
        for (int iFace = 0; iFace < faceNumberOfEachCell[iCell]; ++iFace)
        {
            //! The global ID of face.
            int faceID = cell2Face[iCell][iFace];

            //! Calculate the dot of Vector from particle to surface center..
            RDouble particleToFaceX = (*onePointCoordinate)[0] - xfc[faceID];
            RDouble particleToFaceY = (*onePointCoordinate)[1] - yfc[faceID];
            RDouble particleToFaceZ = (*onePointCoordinate)[2] - zfc[faceID];
            RDouble dotParticleToSurface = particleToFaceX * xfn[faceID]
                + particleToFaceY * yfn[faceID]
                + particleToFaceZ * zfn[faceID];
            RDouble distanceParticleToFace = DISTANCE(particleToFaceX, particleToFaceY, particleToFaceZ);

            //! Calculate the dote of Vector from cell center to surface center.
            RDouble centerToFaceX = xcc[iCell] - xfc[iCell];
            RDouble centerToFaceY = ycc[iCell] - yfc[iCell];
            RDouble centerToFaceZ = zcc[iCell] - zfc[iCell];
            RDouble dotCenterToSurface = centerToFaceX * xfn[iCell]
                + centerToFaceY * yfn[iCell]
                + centerToFaceZ * zfn[iCell];

            //! Judge if the particle in current cell.
            if (dotParticleToSurface * dotCenterToSurface > zero)
            {
                isParticleInCell = true;
            }
            else if (ABS(dotParticleToSurface / distanceParticleToFace) < small)
            {
                //! Particle is close to surface.
                isParticleInCell = true;
            }
            else
            {
                isParticleInCell = false;
                break;
            }

            indexForParticleInCell++;
        }

        if (isParticleInCell)
        {
            indexForParticleInCell = indexForCell;
            break;
        }
    }

    if (isParticleInCell)
    {
        inZone = true;
    }
    else
    {
        inZone = false;
    }

    return indexForParticleInCell;
}

void IndexCell::SetBoundaryCellIndex(Grid *grid)
{
    int gridType = grid->Type();

    if (PHSPACE::STRUCTGRID == gridType)
    {
        StructGrid *structgrid = StructGridCast(grid);
        SetBoundaryCellIndex(structgrid);
    }
    else if (PHSPACE::UNSTRUCTGRID == gridType)
    {
        UnstructGrid *unstructgrid = UnstructGridCast(grid);
        SetBoundaryCellIndex(unstructgrid);
    }
    else if (PHSPACE::MIXGRID == gridType)
    {
        //! Mixing
        TK_Exit::UnexpectedVarValue(" gridType = ", PHSPACE::MIXGRID);
    }
}

void IndexCell::SetBoundaryCellIndex(StructGrid *grid)
{
    const int nDimVar = 3;

    //! Get bc face index and bc type.
    int nBCFace = grid->GetNBoundFace();
    StructBCSet *structBCSet = grid->GetStructBCSet();
    //! Loop for each bc on current zone.
    //int iLocalFace = 0;
    int nBCRegion = structBCSet->GetnBCRegion();
    int bcSet[nDimVar * 2];
    for (int iBC = 0; iBC < nDimVar * 2; ++iBC)
    {
        bcSet[iBC] = 0;
    }

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++iBCRegion)
    {
        StructBC *structBC = structBCSet->GetBCRegion(iBCRegion);
        int bcType = structBC->GetBCType();
        int particleBCType = structBC->GetParticleBCType();

        int snd = structBC->GetFaceDirection();
        int slr = structBC->GetFaceLeftOrRightIndex();
        if (0 == snd)
        {
            if (-1 == slr)
            {
                bcSet[0] = particleBCType;
            }
            else if (1 == slr)
            {
                bcSet[1] = particleBCType;
            }
        }
        else if (1 == snd)
        {
            if (-1 == slr)
            {
                bcSet[2] = particleBCType;
            }
            else if (1 == slr)
            {
                bcSet[3] = particleBCType;
            }
        }
        else if (2 == snd)
        {
            if (-1 == slr)
            {
                bcSet[4] = particleBCType;
            }
            else if (1 == slr)
            {
                bcSet[5] = particleBCType;
            }
        }
    }
    //! bcRegionIDofIFace only in interface.

    //! Index for cell center.
    //! Index start and end of cell in current zone, layer is 2.
    int ist, ied;
    int jst, jed;
    int kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, this->nLayersBC);
    int nLayers = GetNumberOfGhostCellLayers();
    //! Get index IJK of particle cell position in current zone.
    int nI, nJ, nK;
    int iP, jP, kP;
    grid->GetND(nI, nJ, nK);
    nI -= 1;
    nJ -= 1;
    nK -= 1;
    if (TWO_D == GetDim())
    {
        nK = 1;
    }
    int nCellDir[nDimVar];
    nCellDir[0] = nI;
    nCellDir[1] = nJ;
    nCellDir[2] = nK;

    int nFaceOfCell = GetDim() * 2;

    //! Set cell type and bc type for ghost with corner.
    for (iterCell iter = this->cellParticleID->begin(); iter != this->cellParticleID->end(); ++iter)
    {
        int cellID = iter->first;
        GetIndexCellIJKStruct(nI, nJ, nK, iP, jP, kP, nLayers, cellID);
        int indexCellPoint[nDimVar];
        indexCellPoint[0] = iP;
        indexCellPoint[1] = jP;
        indexCellPoint[2] = kP;

        IndexParticleOfCell *indexParticleOfCell = iter->second;
        //! init as NO_PARTICLE_BC on each face.
        indexParticleOfCell->InitBCIndex(nFaceOfCell);

        using namespace INDEX_CELLTYPE;

        //! i:0,1, j:2,3, k:4,5
        for (int iDim = 0; iDim < GetDim(); ++iDim)
        {
            //! Ghost whitout corner.
            int iBCFace = 0;
            if (0 == indexCellPoint[iDim])
            {
                iBCFace = iDim * 2 + 1;
                indexParticleOfCell->SetBCType(iBCFace, bcSet[iDim * 2]);
            }

            if (1 == indexCellPoint[iDim])
            {
                iBCFace = iDim * 2;
                indexParticleOfCell->SetBCType(iBCFace, bcSet[iDim * 2]);
            }

            if (nCellDir[iDim] == indexCellPoint[iDim])
            {
                iBCFace = iDim * 2 + 1;
                indexParticleOfCell->SetBCType(iBCFace, bcSet[iDim * 2 + 1]);
            }

            if (nCellDir[iDim] + 1 == indexCellPoint[iDim] )
            {
                iBCFace = iDim * 2;
                indexParticleOfCell->SetBCType(iBCFace, bcSet[iDim * 2 + 1]);
            }
        }

        //! for ghost in.
        for (int iDim = 0; iDim < GetDim(); ++iDim)
        {
            if (1 == indexCellPoint[iDim]  || nCellDir[iDim] == indexCellPoint[iDim])
            {
                indexParticleOfCell->SetCellType(GHOST_IN);
            }
        }

        //! for ghost out.
        for (int iDim = 0; iDim < GetDim(); ++iDim)
        {
            if (0 == indexCellPoint[iDim] || nCellDir[iDim] + 1 == indexCellPoint[iDim])
            {
                indexParticleOfCell->SetCellType(GHOST_OUT);
            }
        }

        //! For corner.
        if (TWO_D == GetDim())
        {
            bool corner1, corner2, corner3, corner4;

            //! corner in.
            corner1 = (iP == 1 && jP == 1);
            corner2 = (iP == 1 && jP == nJ);
            corner3 = (iP == nI && jP == 1);
            corner4 = (iP == nI && jP == nJ);

            if (corner1 || corner2 || corner3 || corner4)
            {
                indexParticleOfCell->SetCellType(CORNER_IN);
            }

            //! corner out.
            corner1 = (iP == 0 && jP == 0);
            corner2 = (iP == 0 && jP == nJ + 1);
            corner3 = (iP == nI + 1 && jP == 0);
            corner4 = (iP == nI + 1 && jP == nJ + 1);

            if (corner1 || corner2 || corner3 || corner4)
            {
                indexParticleOfCell->SetCellType(CORNER_OUT);
                //! judge if this is corner cell on mpi.
            }
        }
        else if (THREE_D == GetDim())
        {
            bool corner1, corner2, corner3, corner4;
            bool corner5, corner6, corner7, corner8;

            corner1 = (iP == 1 && jP == 1 && kP == 1);
            corner2 = (iP == 1 && jP == 1 && kP == nK);
            corner3 = (iP == 1 && jP == nJ && kP == 1);
            corner4 = (iP == 1 && jP == nJ && kP == nK);

            corner5 = (iP == nI && jP == 1 && kP == 1);
            corner6 = (iP == nI && jP == 1 && kP == nK);
            corner7 = (iP == nI && jP == nJ && kP == 1);
            corner8 = (iP == nI && jP == nJ && kP == nK);

            if (corner1 || corner2 || corner3 || corner4 ||
                corner5 || corner6 || corner7 || corner8)
            {
                indexParticleOfCell->SetCellType(CORNER_IN);
            }

            corner1 = (iP == 0 && jP == 0 && kP == 0);
            corner2 = (iP == 0 && jP == 0 && kP == nK + 1);
            corner3 = (iP == 0 && jP == nJ + 1 && kP == 0);
            corner4 = (iP == 0 && jP == nJ + 1 && kP == nK + 1);

            corner5 = (iP == nI + 1 && jP == 0 && kP == 0);
            corner6 = (iP == nI + 1 && jP == 0 && kP == nK + 1);
            corner7 = (iP == nI + 1 && jP == nJ + 1 && kP == 0);
            corner8 = (iP == nI + 1 && jP == nJ + 1 && kP == nK + 1);

            if (corner1 || corner2 || corner3 || corner4 ||
                corner5 || corner6 || corner7 || corner8)
            {
                indexParticleOfCell->SetCellType(CORNER_OUT);
            }
        }
    }

    ostringstream ossLog;
    ossLog << " --- check cellType and bcType of particle --- " << "\n";
    ossLog << " nI " << nI << " nJ " << nJ << " nK " << nK << "\n";
    for (iterCell iter = this->cellParticleID->begin(); iter != this->cellParticleID->end(); ++iter)
    {
        int cellID = iter->first;
        GetIndexCellIJKStruct(nI, nJ, nK, iP, jP, kP, nLayers, cellID);
        int indexCell[nDimVar];
        indexCell[0] = iP;
        indexCell[1] = jP;
        indexCell[2] = kP;

        IndexParticleOfCell *indexParticleOfCell = iter->second;
        map<int, int>* bcIndex = indexParticleOfCell->GetBCType();

        using namespace INDEX_CELLTYPE;
        int cellType = indexParticleOfCell->GetCellType();

        bool condition1 = cellType == CORNER_IN;
        bool condition2 = cellType == CORNER_OUT;
        bool condition3 = cellType == GHOST_IN;
        bool condition4 = cellType == GHOST_OUT;

        if (!condition2)
        {
            continue;
        }

        ossLog << " iP jP kP : " << iP << " " << jP << " " << kP << "\n";
        ossLog << "  cellID " << cellID << "\n";
        ossLog << "  cellType : " << indexParticleOfCell->GetCellType() << "\n";
        for (map<int, int>::iterator iterFace = bcIndex->begin(); iterFace != bcIndex->end(); ++iterFace)
        {
            ossLog << "   iBCFace " << iterFace->first << " bcType " << iterFace->second << "\n";
        }

        ossLog << " " << "\n";
    }

    ossLog << endl;
}

void IndexCell::SetBoundaryCellIndex(UnstructGrid *grid)
{

}

void IndexCell::AddParticleID(int idCell, int idParticle)
{
    iterCell iter;
    iter = this->cellParticleID->find(idCell);
    if (iter != this->cellParticleID->end())
    {
        //! If cellParticleID had been init.
        //! The cell ID shoule exist in map.
        iter->second->AddParticleID(idParticle);
    }
    else
    {
        //! The cell ID is not exist in map.
        TK_Exit::UnexpectedVarValue(" Un expect idCell = ", idCell);
    }
}

void IndexCell::RemoveParticleID(int idCell, int idParticle)
{
    iterCell iter;
    iter = this->cellParticleID->find(idCell);
    if (iter != this->cellParticleID->end())
    {
        //! If cellParticleID had been init.
        //! The cell ID shoule exist in map.
        iter->second->RemoveParticleID(idParticle);
    }
    else
    {
        //! The cell ID is not exist in map.
        TK_Exit::UnexpectedVarValue(" Un expect idCell = ", idCell);
    }
}

bool IndexCell::IsCellEmpty(int idCell)
{
    iterCell iter;
    iter = this->cellParticleID->find(idCell);
    if (iter != this->cellParticleID->end())
    {
        //! If cellParticleID had been init.
        //! The cell ID shoule exist in map.
        return iter->second->IsEmpty();
    }
    else
    {
        //! The cell ID is not exist in map.
        TK_Exit::UnexpectedVarValue(" Un expect idCell = ", idCell);
    }
    return false;
}

bool IndexCell::IsCellInMap(int idCell)
{
    iterCell iter;
    iter = this->cellParticleID->find(idCell);
    if (iter != this->cellParticleID->end())
    {
        //! If cellParticleID had been init.
        //! The cell ID shoule exist in map.
        return true;
    }
    else
    {
        return false;
    }
}

map<int, IndexParticleOfCell*> *IndexCell::GetCellParticleID()
{
    return this->cellParticleID;
}

bool IndexCell::IsParticleOnCell(int idCell, int idParticle)
{
    if (this->IsCellEmpty(idCell))
    {
        return false;
    }
    else
    {
        iterCell iter;
        iter = this->cellParticleID->find(idCell);
        if (iter != this->cellParticleID->end())
        {
            //! If cellParticleID had been init.
            //! The cell ID shoule exist in map.
            iter->second->IsExistParticleInCell(idParticle);
        }
        else
        {
            //! The cell ID is not exist in map.
            TK_Exit::UnexpectedVarValue(" Un expect idCell = ", idCell);
        }
    }
    return false;
}

int IndexCell::GetNumOfParticleOnCell(int idCell)
{
    if (this->IsCellEmpty(idCell))
    {
        //TK_Exit::UnexpectedVarValue(" Empty idCell = ", idCell);
        return 0;
    }
    else
    {
        iterCell iter;
        iter = this->cellParticleID->find(idCell);
        if (iter != this->cellParticleID->end())
        {
            //! If cellParticleID had been init.
            //! The cell ID shoule exist in map.
            return iter->second->GetNumOfParticleOfCell();
        }
        else
        {
            //! The cell ID is not exist in map.
            TK_Exit::UnexpectedVarValue(" Un expect idCell = ", idCell);
        }
    }
    return 0;
}

int *IndexCell::GetIDOfParticleForCell(int idCell)
{
    if (this->IsCellEmpty(idCell))
    {
        //!TK_Exit::UnexpectedVarValue(" Empty idCell = ", idCell);
        int *particleID = new int;
        *particleID = -1;
        return particleID;
    }
    else
    {
        iterCell iter;
        iter = this->cellParticleID->find(idCell);
        if (iter != this->cellParticleID->end())
        {
            //! If cellParticleID had been init.
            //! The cell ID shoule exist in map.
            return iter->second->GetParticleIDOfCell();
        }
        else
        {
            //! The cell ID is not exist in map.
            TK_Exit::UnexpectedVarValue(" Un expect idCell = ", idCell);
        }
    }
    return NULL;
}

void IndexCell::DeleteMap()
{
    iterCell iter = cellParticleID->begin();
    while(iter != cellParticleID->end())
    {
        //! Because the IndexParticle here is init by new.
        //! And the IndexParticle 's set is init by new.
        //! So first ,we delete IndexPaticleOfCell.
        //! This delete here will call IndexPaticleOfCell's desturctor function.
        delete iter->second;
        //! Second, we make the pointer is NULL.
        iter->second = NULL;
        //! At last, we delete the node for map.
        cellParticleID->erase(iter++);
    }

    delete cellParticleID;
    cellParticleID = nullptr;
}

bool IndexCell::SetParticleInCellInfo(Grid *grid, OnePointVariable *onePointVariable)
{
    bool ifParticleInIndexCell = false;
    int gridType = grid->Type();
    if (PHSPACE::STRUCTGRID == gridType)
    {
        StructGrid *structGrid = StructGridCast(grid);
        ifParticleInIndexCell = SetParticleInCellInfo(structGrid, onePointVariable);
    }
    else if (PHSPACE::UNSTRUCTGRID == gridType)
    {
        UnstructGrid *unStructGrid = UnstructGridCast(grid);
        ifParticleInIndexCell = SetParticleInCellInfo(unStructGrid, onePointVariable);
    }
    else if (PHSPACE::MIXGRID == gridType)
    {
        //! Mixing
        TK_Exit::UnexpectedVarValue(" gridType = ", PHSPACE::MIXGRID);
    }
    return ifParticleInIndexCell;
}

bool IndexCell::SetParticleInCellInfo(StructGrid *grid, OnePointVariable *onePointVariable)
{
    bool ifParticleInIndexCell = false;

    int nI, nJ, nK;
    int iP, jP, kP;
    grid->GetND(nI, nJ, nK);
    nI -= 1;
    nJ -= 1;
    nK -= 1;
    if (TWO_D == GetDim())
    {
        nK = 1;
    }
    int cellID = onePointVariable->GetCellID();
    GetIndexCellIJKStruct(nI, nJ, nK, iP, jP, kP, GetNumberOfGhostCellLayers(), cellID);

    RDouble particleDiameter = onePointVariable->GetOnePointDiameter();
    SPDouble &particleCoordinate = *(onePointVariable->GetOnePointCoordinate());
    //! Get face normal.
    RDouble4D &xfn = *(grid->GetFaceNormalX());
    RDouble4D &yfn = *(grid->GetFaceNormalY());
    RDouble4D &zfn = *(grid->GetFaceNormalZ());

    int iSurLocal, jSurLocal, kSurLocal;

    map<int, IndexParticleOfCell*>::iterator iterCell = this->cellParticleID->find(cellID);
    if (iterCell != this->cellParticleID->end())
    {
        ifParticleInIndexCell = true;

        using namespace INDEX_CELLTYPE;
        int cellType = iterCell->second->GetCellType();

        map<int, int> *bcIndex = iterCell->second->GetBCType();

        //! Loop for iBCFace.
        for (map<int, int>::iterator iterFace = bcIndex->begin(); iterFace != bcIndex->end(); ++iterFace)
        {
            int iBCFace = iterFace->first;
            int particleBCType = iterFace->second;

            //! iBCFace in 2D : i:0,1, j:2,3.
            //! iBCFace in 3D : i:0,1, j:2,3, k:4,5
            //! iDir : thr direction of face in i,j,k.
            //! ilr : Left or right of boundary condition face (0 or 1).
            //! ilrDir : seem to ilr, 1 left, -1 right, 
            //!          which is to make the face normal always point into the cell.
            int iDir = 0;
            int ilr = 0;
            int ilrDir = 0;
            for (int iDim = 1; iDim <= GetDim(); ++iDim)
            {
                //! The direction of boundary condition face.
                if (iBCFace <= 2 * iDim - 1 && iBCFace > 2 * (iDim - 1) - 1)
                {
                    iDir = iDim;
                    ilr = iBCFace - (2 * (iDim - 1) - 1);
                    //! ilr : 0,left, 1,right.
                    ilr = ilr - 1;
                    if (ilr < 0 || ilr > 1)
                    {
                        ostringstream ossError;
                        ossError << " error on ilr = " << ilr << endl;
                        TK_Exit::ExceptionExit(ossError);
                    }
                    else
                    {
                        if (0 == ilr)
                        {
                            ilrDir = 1;
                        }
                        else
                        {
                            ilrDir = -1;
                        }
                    }
                }
            }

            //! Get the surface (or parallel) index.
            GetNsurfIndex(iDir, iSurLocal, jSurLocal, kSurLocal);
            iSurLocal = iP + ilr * iSurLocal;
            jSurLocal = jP + ilr * jSurLocal;
            kSurLocal = kP + ilr * kSurLocal;

            const int nDimVar = 3;
            RDouble fn[nDimVar];
            fn[0] = xfn(iP, jP, kP, iDir) * ilrDir;
            fn[1] = yfn(iP, jP, kP, iDir) * ilrDir;
            fn[2] = zfn(iP, jP, kP, iDir) * ilrDir;

            //! The center coordinate of surface.
            RDouble fc[nDimVar];
            //! Get the center coodinate of surface.
            grid->FaceCoor(iSurLocal, jSurLocal, kSurLocal, iDir, fc[0], fc[1], fc[2]);

            RDouble particleDiameter = onePointVariable->GetOnePointDiameter();
            SPDouble &particleCoordinate = *(onePointVariable->GetOnePointCoordinate());
            SPDouble &particleVelocity = *(onePointVariable->GetOnePointVelocity());

            RDouble distanceParticleToFace = 0.0;
            RDouble surEqD = -fn[0] * fc[0] - fn[1] * fc[1] - fn[2] * fc[2];
            for (int iDim = 0; iDim < nDimVar; ++iDim)
            {
                distanceParticleToFace += fn[iDim] * particleCoordinate[iDim];
            }
            distanceParticleToFace += surEqD;
            distanceParticleToFace = ABS(distanceParticleToFace);
            distanceParticleToFace /= DISTANCE(fn[0], fn[1], fn[2]);

            if (PARTICLEBCSPACE::FULLY_ELASTIC_COLLISION_WALL == particleBCType)
            {
                if (GHOST_IN == cellType || CORNER_IN == cellType)
                {
                    RDouble checkDistance = particleDiameter * 0.5;
                    if (distanceParticleToFace < checkDistance)
                    {
                        onePointVariable->SetIfOneRK(1);

                        ostringstream ossLog;
                        ossLog << " --- check indexCell Wall in --- " << "\n";
                        ossLog << " zoneGlobalID : " << grid->GetZoneID() << "\n";
                        ossLog << " particle ID : " << onePointVariable->GetParticleID() << "\n";
                        ossLog << " particleCoordinate : ";
                        ossLog << particleCoordinate[0] << " ";
                        ossLog << particleCoordinate[1] << " ";
                        ossLog << particleCoordinate[2] << " ";
                        ossLog << "\n";
                        ossLog << " distanceParticleToFace : " << distanceParticleToFace << "\n";
                        ossLog << " fn " << fn[0] << " " << fn[1] << " " << fn[2] << "\n";
                        ossLog << " surEqD " << surEqD << "\n";
                        ossLog << " iDir " << iDir << "\n";
                        ossLog << " ilr " << ilr << "\n";
                        ossLog << " ip jp kp " << iP << " " << jP << " " << kP << "\n";
                        ossLog << " cellID : " << cellID << "\n";
                        ossLog << " celltYPE : " << cellType << "\n";
                        ossLog << " iBCFace : " << iBCFace << "\n";
                        ossLog << endl;
                        bool writeEachProcessor = true;
                    }
                }
                else if (GHOST_OUT == cellType || CORNER_OUT == cellType)
                {
                    onePointVariable->SetIfOneRK(1);

                    ostringstream ossLog;
                    ossLog << " --- check indexCell wall Out --- " << "\n";
                    ossLog << " zoneGlobalID : " << grid->GetZoneID() << "\n";
                    ossLog << " particle ID : " << onePointVariable->GetParticleID() << "\n";
                    ossLog << " particleCoordinate : ";
                    ossLog << particleCoordinate[0] << " ";
                    ossLog << particleCoordinate[1] << " ";
                    ossLog << particleCoordinate[2] << " ";
                    ossLog << "\n";
                    ossLog << " distanceParticleToFace : " << distanceParticleToFace << "\n";
                    ossLog << " fn " << fn[0] << " " << fn[1] << " " << fn[2] << "\n";
                    ossLog << " surEqD " << surEqD << "\n";
                    ossLog << " iDir " << iDir << "\n";
                    ossLog << " ilr " << ilr << "\n";
                    ossLog << " ip jp kp " << iP << " " << jP << " " << kP << "\n";
                    ossLog << " cellID : " << cellID << "\n";
                    ossLog << " celltYPE : " << cellType << "\n";
                    ossLog << " iBCFace : " << iBCFace << "\n";
                    ossLog << endl;
                    bool writeEachProcessor = true;
                }
            }
            else if (PARTICLEBCSPACE::REMOVE_PARTICLE == particleBCType)
            {
                if (GHOST_OUT == cellType || CORNER_OUT == cellType)
                {
                    bool inZone = false;
                    onePointVariable->SetIsParticleInZone(inZone);
                    //! Note that for areas outside the zone, 
                    //! if it is a wall, the RK loop should be exited immediately.
                    onePointVariable->SetIfOneRK(1);

                    ostringstream ossLog;
                    ossLog << " --- check indexCell REMOVE_PARTICLE OUT--- " << "\n";
                    ossLog << " zoneGlobalID : " << grid->GetZoneID() << "\n";
                    ossLog << " particle ID : " << onePointVariable->GetParticleID() << "\n";
                    ossLog << " particleCoordinate : ";
                    ossLog << particleCoordinate[0] << " ";
                    ossLog << particleCoordinate[1] << " ";
                    ossLog << particleCoordinate[2] << " ";
                    ossLog << "\n";
                    ossLog << " distanceParticleToFace : " << distanceParticleToFace << "\n";
                    ossLog << " fn " << fn[0] << " " << fn[1] << " " << fn[2] << "\n";
                    ossLog << " surEqD " << surEqD << "\n";
                    ossLog << " iDir " << iDir << "\n";
                    ossLog << " ilr " << ilr << "\n";
                    ossLog << " ip jp kp " << iP << " " << jP << " " << kP << "\n";
                    ossLog << " cellID : " << cellID << "\n";
                    ossLog << " celltYPE : " << cellType << "\n";
                    ossLog << " iBCFace : " << iBCFace << "\n";
                    ossLog << endl;
                    bool writeEachProcessor = true;
                }
            }
            else if (PARTICLEBCSPACE::REMOVE_PARTICLE_AS_INIT == particleBCType ||
                PARTICLEBCSPACE::REMOVE_PARTICLE_AS_INIT_ONLY_COORD == particleBCType)
            {
                bool inZone = false;
                onePointVariable->SetIsParticleInZone(inZone);
                onePointVariable->SetIfOneRK(1);
            }
            else if (PARTICLEBCSPACE::INTERFACE == particleBCType)
            {
                if (GHOST_OUT == cellType || CORNER_OUT == cellType)
                {
                    this->AddParticleID(cellID, onePointVariable->GetParticleID());

                    bool inZone = false;
                    onePointVariable->SetIsParticleInZone(inZone);
                }
            }
        }
    }

    return ifParticleInIndexCell;
}

bool IndexCell::SetParticleInCellInfo(UnstructGrid *grid, OnePointVariable *onePointVariable)
{
    bool ifParticleInIndexCell = false;
    return ifParticleInIndexCell;
}


void GetIndexCellIJKStruct(int &nI, int &nJ, int &nK, int &iP, int &jP, int &kP, int nGhost, int &cellID)
{
    //! nI,nJ,nK is cell number on IJK, is not node number (not include ghost and corner).
    //! Note here, cellID is start from 0.
    int nGhostX, nGhostY, nGhostZ, nGhostC;

    nGhostX = 2;
    nGhostY = 2;

    if (1 == nK)
    {
        nGhostZ = 0;
        nGhostC = 4;
    }
    else
    {
        nGhostZ = 2;
        nGhostC = 8;
    }

    int nIJ = (nI + nGhost * nGhostX) * (nJ + nGhost * nGhostY);

    kP = floor((cellID + 1.0) / nIJ) + 1;
    jP = floor((cellID + 1.0 - (kP - 1.0) * nIJ) / (nI + nGhostX * nGhost)) + 1;
    iP = cellID + 1 - (kP - 1) * nIJ - (jP - 1) * (nI + nGhost * nGhostX);

    iP -= nGhostX;
    jP -= nGhostY;
    kP -= nGhostZ;
}

void GetIDCellOfParticleByIJKSturct(int &nI, int &nJ, int &nK, int &iP, int &jP, int &kP, int nGhost, int &cellID)
{
    //! nI,nJ,nK is cell number on IJK, is not node number.
    //! Note here, cellID is start from 0.
    int nGhostX, nGhostY, nGhostZ, nGhostC;

    nGhostX = 2;
    nGhostY = 2;

    if (1 == nK)
    {
        nGhostZ = 0;
        nGhostC = 4;
    }
    else
    {
        nGhostZ = 2;
        nGhostC = 8;
    }

    int nIJ = (nI + nGhost * nGhostX) * (nJ + nGhost * nGhostY);

    cellID = iP + nGhost * nGhostX / 2 + (jP + nGhost * nGhostY / 2 - 1) * (nI + nGhost * nGhostX) + (kP + nGhost * nGhostZ / 2 - 1) * nIJ - 1;
}


}