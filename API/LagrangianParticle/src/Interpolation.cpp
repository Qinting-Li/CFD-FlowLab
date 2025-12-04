#include "Interpolation.h"
#include "TK_Exit.h"
#include "GridType.h"
#include "Glb_Dimension.h"

namespace PHSPACE
{

vector< ParticleInterpolation* > *LocalParticleInterpolation::localInterpolationList = new vector< ParticleInterpolation* >;
vector<int> *LocalParticleInterpolation::initInterpolation = new vector<int>;

ParticleInterpolation *LocalParticleInterpolation::GetInterpolation(int iZone)
{
    return (*localInterpolationList)[iZone];
}

bool LocalParticleInterpolation::JudgeIfInitInterpolation()
{
    using namespace PHMPI;
    //! The number of zones on current region.
    int nLocalZones = PHMPI::GetNumberofLocalZones();

    int checkInit = 0;

    for (int iZone = 0; iZone < nLocalZones; ++iZone)
    {
        checkInit += (*initInterpolation)[iZone];
    }

    if (checkInit != nLocalZones)
    {
        return false;
    }
    else
    {
        return true;
    }
}

void LocalParticleInterpolation::InitParticleInterpolation(Grid *grid)
{
    string particleInterpolation_struct = GlobalDataBase::GetStrParaFromDB("particleInterpolation_struct");

    //! The plan of Interpolation
    int gridType = grid->Type();
    if (PHSPACE::STRUCTGRID == gridType)
    {
        if ("TrilinearInterpolation" == particleInterpolation_struct)
        {
            ParticleInterpolation *particleInterpolation = new ParticleTrilinearInterpolation(grid);
            localInterpolationList->push_back(particleInterpolation);
            initInterpolation->push_back(1);
        }
        else if ("UserDefine" == particleInterpolation_struct)
        {
            ParticleInterpolation *particleInterpolation = new ParticleInterpolation(grid);
            initInterpolation->push_back(0);
        }
    }
    else if (PHSPACE::UNSTRUCTGRID == gridType)
    {
        //! This part has not been completed, and will be added later.
        ParticleInterpolation *particleInterpolation = new ParticleInterpolation(grid);
        initInterpolation->push_back(0);
    }
    else
    {
        ParticleInterpolation *particleInterpolation = new ParticleInterpolation(grid);
        initInterpolation->push_back(0);
    }
}

ParticleInterpolation::ParticleInterpolation()
{
    this->grid = 0;
}

ParticleInterpolation::ParticleInterpolation(Grid *gridIn)
{
    this->grid = gridIn;
}

Grid *ParticleInterpolation::GetGrid()
{
    return this->grid;
}

void ParticleInterpolation::SetGrid(Grid *gridIn)
{
    this->grid = gridIn;
}

void ParticleInterpolation::ComputeNodeVar()
{

}

void ParticleInterpolation::ComputePointVar()
{

}

ParticleTrilinearInterpolation::ParticleTrilinearInterpolation()
{
    nPoint = 8;
    pointToInterpolation = new vector<RDouble*>;
}

ParticleTrilinearInterpolation::ParticleTrilinearInterpolation(Grid *gridIn)
{
    this->grid = gridIn;
    nPoint = 8;
    pointToInterpolation = new vector<RDouble*>;
}

ParticleTrilinearInterpolation::~ParticleTrilinearInterpolation()
{
    DeletePointToInterpolation();
}

void ParticleTrilinearInterpolation::DeletePointToInterpolation()
{
    delete pointToInterpolation;
}

void ParticleTrilinearInterpolation::InitPointToInterpolation()
{
    int nDim = 3;
    for (int iPoint = 0; iPoint < nPoint; ++iPoint)
    {
        RDouble *point = new RDouble(nDim);
        pointToInterpolation->push_back(point);
    }
}

void ParticleTrilinearInterpolation::SetPointToInterpolation(vector<RDouble*> *pointToInterpolationIn)
{
    this->pointToInterpolation = pointToInterpolationIn;
}

void ParticleTrilinearInterpolation::ComputeNodeVar(RDouble4D &q, RDouble4D &qNode, int &nVar)
{
    StructGrid *structGrid = StructGridCast(this->grid);

    int ist, ied, jst, jed, kst, ked;
    //! Cell id start from 1, 0 is ghost id
    structGrid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 1);

    RDouble3D &structX = *(structGrid->GetStructX());
    RDouble3D &structY = *(structGrid->GetStructY());
    RDouble3D &structZ = *(structGrid->GetStructZ());

    RDouble3D &xcc = *(structGrid->GetCellCenterX());
    RDouble3D &ycc = *(structGrid->GetCellCenterY());
    RDouble3D &zcc = *(structGrid->GetCellCenterZ());

    //! Temporary variable of linear interpolation.
    RDouble f0, f1, f2, f3;

    //! Note here,ied is cell index.
    ied = ied + 1;
    jed = jed + 1;

    int dI, dJ, dK;
    dI = 1;
    dJ = 1;
    dK = 1;
    if (TWO_D == GetDim())
    {
        dK = 0;
    }
    else
    {
        ked = ked + 1;
    }

    for (int k = kst; k <= ked; ++k)
    {
        for (int j = jst; j <= jed; ++j)
        {
            for (int i = ist; i <= ied; ++i)
            {
                for (int m = 0; m < nVar; ++m)
                {
                    f0 = LinearInterpolation(xcc(i - dI, j - dJ, k - dK), q(i - dI, j - dJ, k - dK, m), xcc(i, j - dJ, k - dK), q(i, j - dJ, k - dK, m), structX(i, j, k));
                    f1 = LinearInterpolation(xcc(i - dI, j     , k - dK), q(i - dI, j     , k - dK, m), xcc(i, j     , k - dK), q(i, j     , k - dK, m), structX(i, j, k));
                    f2 = LinearInterpolation(xcc(i - dI, j - dJ, k     ), q(i - dI, j - dJ, k     , m), xcc(i, j - dJ, k     ), q(i, j - dJ, k     , m), structX(i, j, k));
                    f3 = LinearInterpolation(xcc(i - dI, j     , k     ), q(i - dI, j     , k     , m), xcc(i, j     , k     ), q(i, j     , k     , m), structX(i, j, k));
                    
                    //! For y direction.
                    f0 = LinearInterpolation(ycc(i, j - dJ, k - dK), f0, ycc(i, j, k - dK), f1, structY(i, j, k));
                    f1 = LinearInterpolation(ycc(i, j - dJ, k     ), f2, ycc(i, j, k     ), f3, structY(i, j, k));

                    //! For k direction.
                    if (TWO_D == GetDim())
                    {
                        qNode(i, j, k, m) = f0;
                    }
                    else
                    {
                        f0 = LinearInterpolation(zcc(i, j, k - dK), f0, zcc(i, j, k), f1, structZ(i, j, k));
                        qNode(i, j, k, m) = f0;
                    }
                }
            }
        }
    }
}

void ParticleTrilinearInterpolation::ComputePointVar(RDouble4D &qNode, ParticlePointGroup *particlePointGroup, LDouble *lagrangianVariable, int *varIndex, int *inIndex, int &nVar)
{
    StructGrid *structGrid = StructGridCast(this->grid);
    int ist, ied, jst, jed, kst, ked;
    structGrid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 1);

    RDouble3D &structX = *(structGrid->GetStructX());
    RDouble3D &structY = *(structGrid->GetStructY());
    RDouble3D &structZ = *(structGrid->GetStructZ());

    //! Temporary variable of linear interpolation.
    RDouble f0, f1, f2, f3;

    int dI, dJ, dK;
    dI = 1;
    dJ = 1;
    dK = 1;
    if (TWO_D == GetDim())
    {
        dK = 0;
    }

    //! Get index IJK of particle position in current zone.
    int nI, nJ, nK;
    // int iP, jP, kP;
    structGrid->GetND(nI, nJ, nK);
    nI -= 1;
    nJ -= 1;
    nK -= 1;
    if (TWO_D == GetDim())
    {
        nK = 1;
    }
    
    int cellID;
    //! Loop for grid cells on current zone.
    for (int k = kst; k <= ked; ++k)
    {
        for (int j = jst; j <= jed; ++j)
        {
            for (int i = ist; i <= ied; ++i)
            {
                //! Get the num of particle on current cell.

                GetIDCellOfParticleByIJKSturct(nI, nJ, nK, i, j, k, 2, cellID);
                int nParticleOfCell = particlePointGroup->GetNumOfParticleOfCell(cellID);

                //! Loop for particle on current cell.
                //! If there are no particle on current cell.
                //! The loop will over auto,so don's need breal.
                for (int iLocalParticle = 0; iLocalParticle < nParticleOfCell; ++iLocalParticle)
                {
                    //! The list of global particleID on current grid cell.
                    int *particleIDList = particlePointGroup->GetIDOfParticleOfCell(cellID);
                    //! Get the global particleID.
                    int particleID = particleIDList[iLocalParticle];

                    //! Get the coordinate of iParticle-th particle.
                    //! If here we use:
                    //! SimplePointer<RDouble> coordinate = *(particlePointGroup->GetParticleCoordinate(particleID));
                    //! It will call destructor after for loop.
                    //! So we need to use &,which will not call destructor function.
                    SimplePointer<RDouble> &coordinate = *(particlePointGroup->GetParticleCoordinate(particleID));

                    for (int m = 0; m < nVar; ++m)
                    {
                        //! For x direction.
                        // cout << structX(i, j, k) << " " << structX(i + dI, j, k) << " " << coordinate[0] << endl;
                        f0 = LinearInterpolation(structX(i, j, k), qNode(i, j, k, varIndex[m]), structX(i + dI, j, k), qNode(i + dI, j, k, varIndex[m]), coordinate[0]);
                        f1 = LinearInterpolation(structX(i, j + dJ, k), qNode(i, j + dJ, k, varIndex[m]), structX(i + dI, j + dJ, k), qNode(i + dI, j + dJ, k, varIndex[m]), coordinate[0]);
                        f2 = LinearInterpolation(structX(i, j, k + dK), qNode(i, j, k + dK, varIndex[m]), structX(i + dI, j, k + dK), qNode(i + dI, j, k + dK, varIndex[m]), coordinate[0]);
                        f3 = LinearInterpolation(structX(i, j + dJ, k + dK), qNode(i, j + dJ, k + dK, varIndex[m]), structX(i + dI, j + dJ, k + dK), qNode(i + dI, j + dJ, k + dK, varIndex[m]), coordinate[0]);

                        //! For y direction.
                        f0 = LinearInterpolation(structY(i, j, k), f0, structY(i, j + dJ, k), f1, coordinate[1]);
                        f1 = LinearInterpolation(structY(i, j, k + dK), f2, structY(i, j + dJ, k + dK), f3, coordinate[1]);

                        //! For k direction.
                        if (TWO_D == GetDim())
                        {
                            lagrangianVariable->CopyValue(particleID, f0, inIndex[m]);
                        }
                        else
                        {
                            f0 = LinearInterpolation(structZ(i, j, k), f0, structZ(i, j, k + dK), f1, coordinate[2]);
                            lagrangianVariable->CopyValue(particleID, f0, inIndex[m]);
                        }
                    }
                    
                    if (nParticleOfCell > 0)
                    {

                    }
                }
            }
        }
    }
}

RDouble ParticleTrilinearInterpolation::LinearInterpolation(RDouble &x0, RDouble &f0, RDouble &x1, RDouble &f1, RDouble &x)
{
    if ((x - x0) * (x - x1) > 0)
    {
        TK_Exit::UnexpectedVarValue(" x = ", x);
    }
    return f0 * (x1 - x) / (x1 - x0) + f1 * (x - x0) / (x1 - x0);
}

}