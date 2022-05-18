#include "Geo_MultiGridInfo_Struct.h"
#include "Geo_SimpleBC.h"
#include "GridType.h"
#include "Geo_NodeTopo_Struct.h"
#include "Geo_FaceMetrics_Struct.h"
#include "Geo_CellMetrics_Struct.h"
#include "Geo_DynamicGridMetrics_Struct.h"
#include "Geo_OversetGridTopo_Struct.h"
#include "PHIO.h"
#include "LinkStruct.h"
#include "Glb_Dimension.h"
#include "Constants.h"
#include "Math_Limiter.h"
#include "PHHeader.h"
#include "OversetGridFactory.h"
#include "Geo_StructBC.h"
#include "Geo_Interface.h"
#include "Geo_Interpoint.h"
#include "TK_Log.h"
#include "TK_Parse.h"
#include "FieldProxy.h"
using namespace std;

#pragma warning (disable:913)
namespace PHSPACE
{
int mide[3][3];

class InitializingClass
{
public:
    InitializingClass();
    ~InitializingClass();
public:
    void InitConstants();
};

InitializingClass::InitializingClass()
{
    InitConstants();
}

InitializingClass::~InitializingClass()
{
}

void InitializingClass::InitConstants()
{
    for (int i = 0; i < 3; ++ i)
    {
        for (int j = 0; j < 3; ++ j)
        {
            mide[i][j] = 0;
        }
    }
    mide[0][0] = 1;
    mide[1][1] = 1;
    mide[2][2] = 1;
}

InitializingClass AnInstance;

StructGrid * StructGridCast(Grid *gridIn)
{
    return static_cast <StructGrid *> (gridIn);
}

LIB_EXPORT StructGrid::StructGrid()
{
    nodeTopology = NULL;
    faceMetrics = NULL;
    cellMetrics = NULL;
    multigridInfo = NULL;
    dynamicGridMetrics = NULL;
    walldist = NULL;
    iBlank3D = NULL;
    gridUnstr = NULL;
    wallCellLabel = NULL;
    composite_bcregion = NULL;
    oversetGridTopology = NULL;

    zoneProbesCellNI.resize(0);
    zoneProbesCellNJ.resize(0);
    zoneProbesCellNK.resize(0);

    ordinaryGridIndex = -1;
    for (int iDim = 0; iDim < THREE_D; ++ iDim)
    {
        ordinaryDimStartIndex[iDim] = 0;
        ordinaryDimEndIndex  [iDim] = 0;
    }
}

LIB_EXPORT void StructGrid::InitGrid(GridID *index, int level, int dim, int type)
{
    Grid::InitGrid(index, level, dim, STRUCTGRID);

    nodeTopology  = new Geo_NodeTopo_Struct();
    faceMetrics   = new Geo_FaceMetrics_Struct();
    cellMetrics   = new Geo_CellMetrics_Struct();
    multigridInfo = new Geo_MultiGridInfo_Struct();
    dynamicGridMetrics = new Geo_DynamicGridMetrics_Struct();

    int isOverset = GlobalDataBase::GetIntParaFromDB("codeOfOversetGrid");
    if (isOverset)
    {
        gridUnstr = new UnstructGrid();
    }

    walldist = 0;
    composite_bcregion = 0;
    oversetGridTopology = new Geo_OversetGridTopo_Struct();
}

LIB_EXPORT StructGrid::~StructGrid()
{
    if (nodeTopology) { delete nodeTopology; nodeTopology = NULL; }
    if (faceMetrics) { delete faceMetrics; faceMetrics = NULL; }
    if (cellMetrics) { delete cellMetrics; cellMetrics = NULL; }
    if (multigridInfo) { delete multigridInfo; multigridInfo = NULL; }
    if (dynamicGridMetrics) { delete dynamicGridMetrics; dynamicGridMetrics = NULL; }
    if (gridUnstr) { delete gridUnstr; gridUnstr = NULL; }
    if (oversetGridTopology) { delete oversetGridTopology; oversetGridTopology = NULL; }
    delete walldist;
    delete composite_bcregion;

    zoneProbesCellNI.clear();
    zoneProbesCellNJ.clear();
    zoneProbesCellNK.clear();

    //! Added by Guo Yongheng 20160818.
    //DeleteLinkPointStructure();
    //DeleteCellContainer();
    //DeleteCoreParameterContainer();
}

LIB_EXPORT void StructGrid::CreateCompositeBCRegion(int nBCRegion)
{
    if (!composite_bcregion)
    {
        composite_bcregion = new CompositeStructBC(this->GetZoneID());
    }

    composite_bcregion->CreateBCRegion(nBCRegion);
}

LIB_EXPORT void StructGrid::SetBCFaceInfo()
{ 
    this->composite_bcregion->SetIFaceInfo();
    int nBoundFace = composite_bcregion->GetNBoundFace();

    int ni = this->GetNI();
    int nj = this->GetNJ();
    int nk = this->GetNK();
    int ndim = GetDim();

    if (ndim == 2)
    {
        if (nBoundFace != 2 * ((ni-1) + (nj-1)))
        {
            PHMPI::FreeBlockData();
            cout << "根据region计算的边界面数nBoundFace与网格块本身边界面数不一致!\n" << endl;
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
    }
    else if (ndim == 3)
    {
        if (nBoundFace != 2 * ((ni-1)*(nj-1) + (ni-1)*(nk-1) + (nj-1)*(nk-1)))
        {
            PHMPI::FreeBlockData();
            cout << "根据region计算的边界面数nBoundFace与网格块本身边界面数不一致!\n" << endl;
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
    }

    this->SetNBoundFace(nBoundFace);
}

LIB_EXPORT void StructGrid::CopyCompositeBCRegion(CompositeStructBC *composite_bcregion)
{
    if (!this->composite_bcregion)
    {
        this->composite_bcregion = new CompositeStructBC(this->GetZoneID());
    }
    this->composite_bcregion->CopyCompositeBCRegion(composite_bcregion);
}

LIB_EXPORT void StructGrid::CopyGrid(StructGrid *grid)
{
    CopyXYZ(grid->GetNI(), grid->GetNJ(), grid->GetNK(), grid->GetX(), grid->GetY(), grid->GetZ());
    CopyCompositeBCRegion(grid->GetCompositeBCRegion());
}

LIB_EXPORT void StructGrid::GetSourceIndexIJK(int iface, int ipos, int &i, int &j, int &k)
{
    this->composite_bcregion->GetSourceIndexIJK(iface, ipos, i, j, k);
}

LIB_EXPORT void StructGrid::GetSourceIndexIJK_Nsurf_LR(int iface, int ipos, int &i, int &j, int &k, int &nsurf, int &s_lr)
{
    this->composite_bcregion->GetSourceIndexIJK_Nsurf_LR(iface, ipos, i, j, k, nsurf, s_lr);
}

LIB_EXPORT void StructGrid::GetTargetIndexIJK(int iface, int ipos, int &i, int &j, int &k)
{
    this->composite_bcregion->GetTargetIndexIJK(iface, ipos, i, j, k);
}

LIB_EXPORT void StructGrid::CopyXYZ(int ni, int nj, int nk, RDouble *xx, RDouble *yy, RDouble *zz)
{
    this->SetNI(ni);
    this->SetNJ(nj);
    this->SetNK(nk);

    this->SetBasicDimension();

    int nTotalNode = this->GetNTotalNode();

    RDouble *x = new RDouble [nTotalNode];
    RDouble *y = new RDouble [nTotalNode];
    RDouble *z = new RDouble [nTotalNode];

    this->SetX(x);
    this->SetY(y);
    this->SetZ(z);
    this->SetArrayLayout();

    for (int inode = 0; inode < nTotalNode; ++ inode)
    {
        x[inode] = xx[inode];
        y[inode] = yy[inode];
        z[inode] = zz[inode];
    }
    cout << "ni, nj, nk = " << ni << " " << nj  << " " << nk << "\n";
    ComputeMinMaxBox();
}

LIB_EXPORT RDouble3D & StructGrid::GetVolume(int istep) const
{
    return * GetCellVolume();
}

LIB_EXPORT void StructGrid::Action(ActionKey *actkey)
{
    switch (actkey->action)
    {
    case WRITE_WALL_DIST:
        DumpWalldist(actkey);
        break;
    case COMMCELLCENTERDATA:
        GetCellCenter(actkey);
        break;
    case COMMCELLIBLANK:
        GetCellIBlank(actkey);
        break;
    case SIMPLE_ACTION:
        SimpleAction(actkey);
        break;
    case ALLOCATE_WALL_DIST:
        AllocateWalldist();
        break;
    default:
        break;
    }
}

void StructGrid::SimpleAction(ActionKey *actkey)
{
    RDouble value;
    GetData(actkey->taskname, &value, PHDOUBLE, 1);

    DataContainer *cdata = actkey->GetData();
    cdata->moveto_begin();
    cdata->write(&value, sizeof(RDouble));
}

LIB_EXPORT void StructGrid::ActionReflect(ActionTag *acttag)
{
}

LIB_EXPORT void StructGrid::TranslateAction(ActionKey *actkey)
{
    if (actkey->action == COMMCELLCENTERDATA)
    {
        TranslateCellCenter(actkey);
    }
    else if (actkey->action == COMMCELLIBLANK)
    {
        TranslateCellIBlank(actkey);
    }
}

LIB_EXPORT streamsize StructGrid::TranslateActionLength(ActionKey *actkey)
{
    if (actkey->action == COMMCELLCENTERDATA)
    {
        return TranslateCellCenterLength(actkey);
    }
    return 0;
}

LIB_EXPORT void StructGrid::AllocateWalldist()
{
    int ni = this->GetNI();
    int nj = this->GetNJ();
    int nk = this->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);

    walldist = new RDouble3D(I, J, K, fortranArray);
}

void StructGrid::ReadWalldist(RDouble *wallDistIn)
{
    RDouble3D &walldist = * this->GetWallDist();

    int ist, ied, jst, jed, kst, ked;
    this->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    int count = 0;
    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                walldist(i, j, k) = wallDistIn[count];
                count ++;
            }
        }
    }
}

LIB_EXPORT void StructGrid::WriteFaceBC(fstream &file)
{
    VirtualFile *virtualFile = new VirtualFile(&file);
    virtualFile->BeginWriteWork();

    int iZone = this->GetGridID()->GetIndex();
    PHWrite(virtualFile, iZone);

    CompositeStructBC *composite_bcregion = this->GetCompositeBCRegion();

    int nBCRegion = static_cast<int> (composite_bcregion->GetnBCRegion());

    PHWrite(virtualFile, nBCRegion);

    string *bcNameList = new string [nBCRegion];

    //StringSizeType totalSize = 0;
    uint_long totalSize = 0;
    for (int ibcregion = 0; ibcregion < nBCRegion; ibcregion ++)
    {
        StructBC *structBC = composite_bcregion->GetBCRegion(ibcregion);

        string bcName = structBC->GetBCName();
        int bctype = structBC->GetBCType();
        if (bcName == "")
        {
            if (bctype == -1)
            {
                bcName = "Connection";
            }
            else
            {
                bcName = composite_bcregion->GetBCName(bctype);
                if ((bctype > 7) && (bctype < 51))
                {
                    bctype = 2;
                    structBC->SetBCType(bctype);
                }
            }
            structBC->SetBCName(bcName);
        }

        bcNameList[ibcregion] = bcName;

        //totalSize += bcNameList[ibcregion].size();
        totalSize += static_cast< uint_long > (bcNameList[ibcregion].size());
        totalSize += 1;
    }

    char *bcNameChar = new char [totalSize];
    unsigned int count = 0;
    for (int boco = 0; boco < nBCRegion; ++ boco)
    {
        string &bcName = bcNameList[boco];
        streamsize nameSize = bcName.size();
        for (unsigned int iChar = 0; iChar < nameSize; ++ iChar)
        {
            bcNameChar[count ++] = bcName[iChar];
        }
        bcNameChar[count ++] = '\0';
    }

    PHWrite(virtualFile, totalSize);
    PHWrite(virtualFile, bcNameChar, static_cast<int>(totalSize));

    delete [] bcNameList;
    delete [] bcNameChar;

    virtualFile->EndWriteWork();

    delete virtualFile;
}

void StructGrid::DumpWalldist(ActionKey *actkey)
{
    int nTotalCell = this->GetNTotalCell();

    DataContainer *cdata = actkey->GetData();
    cdata->moveto_begin();

    RDouble3D &walldist = * this->GetWallDist();

    int ist, ied, jst, jed, kst, ked;
    this->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    cdata->write(&nTotalCell, sizeof(int));

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                cdata->write(&walldist(i, j, k), sizeof(RDouble));
            }
        }
    }
}

LIB_EXPORT void StructGrid::CompNodeVar(RDouble4D &q_n, RDouble4D &q)
{
    int nk = this->GetNK();

    int nl = 5;
    GetData("nl", &nl, PHINT, 1);

    int il1 = 1;
    int jl1 = 1;
    int kl1 = 1;
    if (nk == 1) kl1 = 0;

    int il, jl, kl;

    int ist, ied, jst, jed, kst, ked;
    this->GetNodeIterationIndex(ist, ied, jst, jed, kst, ked);

    Range M(0, nl-1);
    int mst = M.first();
    int med = M.last();

    for (int m = mst; m <= med; ++ m)
    {
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    il = i - il1;
                    jl = j - jl1;
                    kl = k - kl1;
                    q_n(i, j, k, m) = eighth * (q(il, jl, kl, m) + q(il, j, kl, m) + q(i, jl, kl, m) + q(i, j, kl, m) +
                                                q(il, jl, k , m) + q(il, j, k , m) + q(i, jl, k , m) + q(i, j, k , m));
                }
            }
        }
    }
}

LIB_EXPORT void StructGrid::FaceCoor(int i, int j, int k, int nsurf, RDouble &xc, RDouble &yc, RDouble &zc)
{
    int nk = GetNK();
    int il1 = 1;
    int jl1 = 1;
    int kl1 = 1;
    if (nk == 1) kl1 = 0;

    int il, jl, kl;

    il = i + il1;
    jl = j + jl1;
    kl = k + kl1;

    RDouble3D &structx = * this->GetStructX();
    RDouble3D &structy = * this->GetStructY();
    RDouble3D &structz = * this->GetStructZ();

    if (nsurf == 1)
    {
        xc = fourth * (structx(i, j, k) + structx(i, jl, k) + structx(i, jl, kl) + structx(i, j, kl));
        yc = fourth * (structy(i, j, k) + structy(i, jl, k) + structy(i, jl, kl) + structy(i, j, kl));
        zc = fourth * (structz(i, j, k) + structz(i, jl, k) + structz(i, jl, kl) + structz(i, j, kl));
    }
    else if (nsurf == 2)
    {
        xc = fourth * (structx(i, j, k) + structx(il, j, k) + structx(il, j, kl) + structx(i, j, kl));
        yc = fourth * (structy(i, j, k) + structy(il, j, k) + structy(il, j, kl) + structy(i, j, kl));
        zc = fourth * (structz(i, j, k) + structz(il, j, k) + structz(il, j, kl) + structz(i, j, kl));
    }
    else if (nsurf == 3)
    {
        xc = fourth * (structx(i, j, k) + structx(il, j, k) + structx(il, jl, k) + structx(i, jl, k));
        yc = fourth * (structy(i, j, k) + structy(il, j, k) + structy(il, jl, k) + structy(i, jl, k));
        zc = fourth * (structz(i, j, k) + structz(il, j, k) + structz(il, jl, k) + structz(i, jl, k));
    }
}

LIB_EXPORT void StructGrid::CenterCoor(int i, int j, int k, RDouble &xc, RDouble &yc, RDouble &zc)
{
    RDouble3D &xx = * this->GetStructX();
    RDouble3D &yy = * this->GetStructY();
    RDouble3D &zz = * this->GetStructZ();
    int nk = GetNK();

    int il1 = 1;
    int jl1 = 1;
    int kl1 = 1;
    if (nk == 1) kl1 = 0;

    int il, jl, kl;

    il = i + il1;
    jl = j + jl1;
    kl = k + kl1;

    xc = eighth * (xx(i, j, k) + xx(il, j, k) + xx(il, jl, k) + xx(i, jl, k)
                  + xx(i, j, kl) + xx(il, j, kl) + xx(il, jl, kl) + xx(i, jl, kl));
    yc = eighth * (yy(i, j, k) + yy(il, j, k) + yy(il, jl, k) + yy(i, jl, k)
                  + yy(i, j, kl) + yy(il, j, kl) + yy(il, jl, kl) + yy(i, jl, kl));
    zc = eighth * (zz(i, j, k) + zz(il, j, k) + zz(il, jl, k) + zz(i, jl, k)
                  + zz(i, j, kl) + zz(il, j, kl) + zz(il, jl, kl) + zz(i, jl, kl));
}

LIB_EXPORT void StructGrid::SetArrayLayout()
{
    RDouble *x = GetX();
    RDouble *y = GetY();
    RDouble *z = GetZ();
    if (x && y && z)
    {
        nodeTopology->SetArrayLayout(x, y ,z);
    }
}

LIB_EXPORT void StructGrid::AllocateMetrics(ActionKey *actkey)
{
    int ni = this->GetNI();
    int nj = this->GetNJ();
    int nk = this->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);
    Range D(1, 3);

    faceMetrics->AllocateMetrics(I, J, K, D);
    cellMetrics->AllocateMetrics(I, J, K);

    RDouble4D *vgn = dynamicGridMetrics->GetFaceNormalVelocity();
    if (vgn == NULL) vgn = new RDouble4D(I, J, K, Range(1, 3), fortranArray);

    AllocateMetricsALE(actkey);

    (*vgn) = zero;
    dynamicGridMetrics->SetFaceNormalVelocity(vgn);
}

void StructGrid::AllocateMetricsALE(ActionKey *actkey)
{
    int ni = this->GetNI();
    int nj = this->GetNJ();
    int nk = this->GetNK();

    Range IFACE, JFACE, KFACE;
    GetRange(ni, nj, nk, -2, 2, IFACE, JFACE, KFACE);
    Range I, J, K;
    GetRange(ni, nj, nk, -2, 1, I, J, K);
    Range D(1, 3);

    //! codeOfAleModel = 0: no ALE method;
    //!        1: ALE method for non-moving grids;
    //!        2: ALE method for moving grids;
    //!        3: ALE method for deforming grids.

    int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");

    int isAle = GlobalDataBase::GetIntParaFromDB("codeOfAleModel");

    RDouble4D *faceVelocityX = dynamicGridMetrics->GetFaceVelocityX();
    RDouble4D *faceVelocityY = dynamicGridMetrics->GetFaceVelocityY();
    RDouble4D *faceVelocityZ = dynamicGridMetrics->GetFaceVelocityZ();
    RDouble3D *voln = dynamicGridMetrics->GetCellVolumeOld();

    if (isUnsteady && isAle)
    {
        if (faceVelocityX == NULL)
        {
            faceVelocityX = new RDouble4D(I, J, K, Range(1, 3), fortranArray);
            (*faceVelocityX) = zero;
        }
        if (faceVelocityY == NULL)
        {
            faceVelocityY = new RDouble4D(I, J, K, Range(1, 3), fortranArray);
            (*faceVelocityY) = zero;
        }
        if (faceVelocityZ == NULL)
        {
            faceVelocityZ = new RDouble4D(I, J, K, Range(1, 3), fortranArray);
            (*faceVelocityZ) = zero;
        }
    }

    if (isUnsteady && isAle)
    {
        if (voln == NULL) voln = new RDouble3D(I, J, K, fortranArray);
    }
    else
    {
        if (!isAle)
        {
            //voln = vol;
            voln = 0;    //! This is not sure.
        }
    }

    dynamicGridMetrics->SetFaceVelocityX(faceVelocityX);
    dynamicGridMetrics->SetFaceVelocityY(faceVelocityY);
    dynamicGridMetrics->SetFaceVelocityZ(faceVelocityZ);
    dynamicGridMetrics->SetCellVolumeOld(voln);
}

LIB_EXPORT void StructGrid::ComputeMetrics(ActionKey *actkey)
{
    int isFVMOrFDM = PHSPACE::GlobalDataBase::GetIntParaFromDB("isFVMOrFDM");

    if (isFVMOrFDM == FD_METHOD)
    {
        AllocateMetricsStructHighOrder(actkey);
        if (GetDim() == TWO_D)
        {
            ComputeMetricsStructHighOrder2D(actkey);
            GhostMetricsStructHighOrder2D(actkey);
        }
        else
        {
            ComputeMetricsStructHighOrder3D(actkey);
            GhostMetricsStructHighOrder3D(actkey);
        }
    }
    else
    {
        AllocateMetrics(actkey);
        if (GetDim() == TWO_D)
        {
            ComputeMetrics2D(actkey);
        }
        else
        {
            ComputeMetrics3D(actkey);
        }
    }
}

void StructGrid::ComputeMetrics3D(ActionKey *actkey)
{
    RDouble r1[3], r2[3], r3[3], r4[3], r5[3], r6[3], r7[3], r8[3];
    RDouble vkc[3], vet[3], vct[3], vkcn[4], vetn[4], vctn[4];
    RDouble v1, v2, v3, v4, v5, v6, vv;
    int i, j, k;
    int ni = GetNI();
    int nj = GetNJ();
    int nk = GetNK();

    RDouble volmin = LARGE;
    RDouble volmax = - LARGE;
    int imin = 1;
    int jmin = 1;
    int kmin = 1;
    int imax = 1;
    int jmax = 1;
    int kmax = 1;

    //! Set surfaces and volume to zero.
    RDouble4D *xfn  = GetFaceNormalX();
    RDouble4D *yfn  = GetFaceNormalY();
    RDouble4D *zfn  = GetFaceNormalZ();
    RDouble4D *area = GetFaceArea();
    RDouble3D *vol  = GetCellVolume();
    RDouble4D *xfv  = GetFaceVectorX();
    RDouble4D *yfv  = GetFaceVectorY();
    RDouble4D *zfv  = GetFaceVectorZ();
    (*xfn) = 0.0;
    (*yfn) = 0.0;
    (*zfn) = 0.0;
    (*area) = 0.0;
    (*vol) = 0.0;
    (*xfv) = 0.0;
    (*yfv) = 0.0;
    (*zfv) = 0.0;

    std::ostringstream oss;

    RDouble3D &xx = * this->GetStructX();
    RDouble3D &yy = * this->GetStructY();
    RDouble3D &zz = * this->GetStructZ();

    int zIndex = this->GetZoneID();
    for (k = 1; k <= nk - 1; ++ k)
    {
        for (j = 1; j <= nj - 1; ++ j)
        {
            for (i = 1; i <= ni - 1; ++ i)
            {
                r1[0] = xx(i  , j  , k );
                r1[1] = yy(i  , j  , k );
                r1[2] = zz(i  , j  , k );

                r2[0] = xx(i+1, j  , k );
                r2[1] = yy(i+1, j  , k );
                r2[2] = zz(i+1, j  , k );

                r3[0] = xx(i  , j+1, k );
                r3[1] = yy(i  , j+1, k );
                r3[2] = zz(i  , j+1, k );

                r4[0] = xx(i  , j  , k+1);
                r4[1] = yy(i  , j  , k+1);
                r4[2] = zz(i  , j  , k+1);

                r5[0] = xx(i+1, j+1, k );
                r5[1] = yy(i+1, j+1, k );
                r5[2] = zz(i+1, j+1, k );

                r6[0] = xx(i  , j+1, k+1);
                r6[1] = yy(i  , j+1, k+1);
                r6[2] = zz(i  , j+1, k+1);

                r7[0] = xx(i+1, j  , k+1);
                r7[1] = yy(i+1, j  , k+1);
                r7[2] = zz(i+1, j  , k+1);

                r8[0] = xx(i+1, j+1, k+1);
                r8[1] = yy(i+1, j+1, k+1);
                r8[2] = zz(i+1, j+1, k+1);

                Vol_Tetrahedron(r1, r2, r5, r8, v1);
                Vol_Tetrahedron(r1, r2, r8, r7, v2);
                Vol_Tetrahedron(r1, r3, r8, r5, v3);
                Vol_Tetrahedron(r1, r3, r6, r8, v4);
                Vol_Tetrahedron(r1, r4, r8, r6, v5);
                Vol_Tetrahedron(r1, r4, r7, r8, v6);

                vv = sixth * (ABS(v1) + ABS(v2) + ABS(v3) + ABS(v4) + ABS(v5) + ABS(v6));

                (*vol)(i, j, k) = vv;
                if (vv <= 0.0)
                {
                    oss << " vol < 0 " << zIndex << " " << i << " " << j << " " << k << vv << "\n";
                }

                if (vv < volmin)
                {
                    volmin = vv;
                    imin = i;
                    jmin = j;
                    kmin = k;
                }

                if (vv > volmax)
                {
                    volmax = vv;
                    imax = i;
                    jmax = j;
                    kmax = k;
                }

                FaceNormal3D(r1, r3, r6, r4, vkc);
                FaceNormal3D(r1, r4, r7, r2, vet);
                FaceNormal3D(r1, r2, r5, r3, vct);

                NormalizeVector(vkc, vkcn, 3);
                NormalizeVector(vet, vetn, 3);
                NormalizeVector(vct, vctn, 3);

                (*xfn)(i, j, k, 1) = vkcn[0];
                (*yfn)(i, j, k, 1) = vkcn[1];
                (*zfn)(i, j, k, 1) = vkcn[2];
                (*area)(i, j, k, 1) = vkcn[3];
                (*xfv)(i, j, k, 1) = vkc[0];
                (*yfv)(i, j, k, 1) = vkc[1];
                (*zfv)(i, j, k, 1) = vkc[2];

                (*xfn)(i, j, k, 2) = vetn[0];
                (*yfn)(i, j, k, 2) = vetn[1];
                (*zfn)(i, j, k, 2) = vetn[2];
                (*area)(i, j, k, 2) = vetn[3];
                (*xfv)(i, j, k, 2) = vet[0];
                (*yfv)(i, j, k, 2) = vet[1];
                (*zfv)(i, j, k, 2) = vet[2];

                (*xfn)(i, j, k, 3) = vctn[0];
                (*yfn)(i, j, k, 3) = vctn[1];
                (*zfn)(i, j, k, 3) = vctn[2];
                (*area)(i, j, k, 3) = vctn[3];
                (*xfv)(i, j, k, 3) = vct[0];
                (*yfv)(i, j, k, 3) = vct[1];
                (*zfv)(i, j, k, 3) = vct[2];
            }
        }
    }

    //! Compute the grid derivatives on the face of i + 1 = ni.
    i = ni - 1;
    for (k = 1; k <= nk - 1; ++ k)
    {
        for (j = 1; j <= nj - 1; ++ j)
        {
            r2[0] = xx(i + 1, j    , k   );
            r2[1] = yy(i + 1, j    , k   );
            r2[2] = zz(i + 1, j    , k   );

            r5[0] = xx(i + 1, j + 1, k   );
            r5[1] = yy(i + 1, j + 1, k   );
            r5[2] = zz(i + 1, j + 1, k   );

            r7[0] = xx(i + 1, j    , k + 1);
            r7[1] = yy(i + 1, j    , k + 1);
            r7[2] = zz(i + 1, j    , k + 1);

            r8[0] = xx(i + 1, j + 1, k + 1);
            r8[1] = yy(i + 1, j + 1, k + 1);
            r8[2] = zz(i + 1, j + 1, k + 1);

            FaceNormal3D(r2, r5, r8, r7, vkc);
            NormalizeVector(vkc, vkcn, 3);

            (*xfn)(i+1, j, k, 1) = vkcn[0];
            (*yfn)(i+1, j, k, 1) = vkcn[1];
            (*zfn)(i+1, j, k, 1) = vkcn[2];
            (*area)(i+1, j, k, 1) = vkcn[3];
            (*xfv)(i+1, j, k, 1) = vkc[0];
            (*yfv)(i+1, j, k, 1) = vkc[1];
            (*zfv)(i+1, j, k, 1) = vkc[2];
        }
    }

    //! Compute the grid derivatives on the face of j + 1 = nj.
    j = nj - 1;
    for (k = 1; k <= nk - 1; ++ k)
    {
        for (i = 1; i <= ni - 1; ++ i)
        {
            r3[0] = xx(i    , j + 1, k   );
            r3[1] = yy(i    , j + 1, k   );
            r3[2] = zz(i    , j + 1, k   );

            r5[0] = xx(i + 1, j + 1, k   );
            r5[1] = yy(i + 1, j + 1, k   );
            r5[2] = zz(i + 1, j + 1, k   );

            r6[0] = xx(i    , j + 1, k + 1);
            r6[1] = yy(i    , j + 1, k + 1);
            r6[2] = zz(i    , j + 1, k + 1);

            r8[0] = xx(i + 1, j + 1, k + 1);
            r8[1] = yy(i + 1, j + 1, k + 1);
            r8[2] = zz(i + 1, j + 1, k + 1);

            FaceNormal3D(r3, r6, r8, r5, vet);
            NormalizeVector(vet, vetn, 3);

            (*xfn)(i, j+1, k, 2) = vetn[0];
            (*yfn)(i, j+1, k, 2) = vetn[1];
            (*zfn)(i, j+1, k, 2) = vetn[2];
            (*area)(i, j+1, k, 2) = vetn[3];
            (*xfv)(i, j+1, k, 2) = vet[0];
            (*yfv)(i, j+1, k, 2) = vet[1];
            (*zfv)(i, j+1, k, 2) = vet[2];
        }
    }

    //! Compute the grid derivatives on the face of k + 1 = nk.
    k = nk - 1;
    for (j = 1; j <= nj - 1; ++ j)
    {
        for (i = 1; i <= ni - 1; ++ i)
        {
            r4[0] = xx(i    , j    , k + 1);
            r4[1] = yy(i    , j    , k + 1);
            r4[2] = zz(i    , j    , k + 1);

            r6[0] = xx(i    , j + 1, k + 1);
            r6[1] = yy(i    , j + 1, k + 1);
            r6[2] = zz(i    , j + 1, k + 1);

            r7[0] = xx(i + 1, j    , k + 1);
            r7[1] = yy(i + 1, j    , k + 1);
            r7[2] = zz(i + 1, j    , k + 1);

            r8[0] = xx(i + 1, j + 1, k + 1);
            r8[1] = yy(i + 1, j + 1, k + 1);
            r8[2] = zz(i + 1, j + 1, k + 1);

            FaceNormal3D(r4, r7, r8, r6, vct);
            NormalizeVector(vct, vctn, 3);

            (*xfn)(i, j, k+1, 3) = vctn[0];
            (*yfn)(i, j, k+1, 3) = vctn[1];
            (*zfn)(i, j, k+1, 3) = vctn[2];
            (*area)(i, j, k+1, 3) = vctn[3];
            (*xfv)(i, j, k+1, 3) = vct[0];
            (*yfv)(i, j, k+1, 3) = vct[1];
            (*zfv)(i, j, k+1, 3) = vct[2];
        }
    }

    GhostCell3D(*vol, ni, nj, nk);
    GhostMetrics3D();

    ComputeCellCenter();

    StreamToActionKey(actkey, oss);
}

void StructGrid::ComputeCellCenter()
{
    RDouble3D &xcc = * this->GetCellCenterX();
    RDouble3D &ycc = * this->GetCellCenterY();
    RDouble3D &zcc = * this->GetCellCenterZ();

    int ist, ied, jst, jed, kst, ked;
    this->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                CenterCoor(i, j, k, xcc(i, j, k), ycc(i, j, k), zcc(i, j, k));
            }
        }
    }

    RDouble xfc, yfc, zfc;

    int nsurf;
    nsurf = 1;
    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            int i;
            i = ist;

            FaceCoor(i, j, k, nsurf, xfc, yfc, zfc);
            
            xcc(i-1, j, k) = two * xfc - xcc(i, j, k);
            ycc(i-1, j, k) = two * yfc - ycc(i, j, k);
            zcc(i-1, j, k) = two * zfc - zcc(i, j, k);

            i = ied;

            FaceCoor(i+1, j, k, nsurf, xfc, yfc, zfc);
            
            xcc(i+1, j, k) = two * xfc - xcc(i, j, k);
            ycc(i+1, j, k) = two * yfc - ycc(i, j, k);
            zcc(i+1, j, k) = two * zfc - zcc(i, j, k);
        }
    }

    nsurf = 2;
    for (int k = kst; k <= ked; ++ k)
    {
        for (int i = ist; i <= ied; ++ i)
        {
            int j;
            j = jst;

            FaceCoor(i, j, k, nsurf, xfc, yfc, zfc);
            
            xcc(i, j-1, k) = two * xfc - xcc(i, j, k);
            ycc(i, j-1, k) = two * yfc - ycc(i, j, k);
            zcc(i, j-1, k) = two * zfc - zcc(i, j, k);

            j = jed;

            FaceCoor(i, j+1, k, nsurf, xfc, yfc, zfc);
            
            xcc(i, j+1, k) = two * xfc - xcc(i, j, k);
            ycc(i, j+1, k) = two * yfc - ycc(i, j, k);
            zcc(i, j+1, k) = two * zfc - zcc(i, j, k);
        }
    }

    if (this->GetDim() == THREE_D)
    {
        nsurf = 3;
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                int k;
                k = kst;

                FaceCoor(i, j, k, nsurf, xfc, yfc, zfc);
                
                xcc(i, j, k-1) = two * xfc - xcc(i, j, k);
                ycc(i, j, k-1) = two * yfc - ycc(i, j, k);
                zcc(i, j, k-1) = two * zfc - zcc(i, j, k);

                k = ked;

                FaceCoor(i, j, k+1, nsurf, xfc, yfc, zfc);
                
                xcc(i, j, k+1) = two * xfc - xcc(i, j, k);
                ycc(i, j, k+1) = two * yfc - ycc(i, j, k);
                zcc(i, j, k+1) = two * zfc - zcc(i, j, k);
            }
        }
    }
}

void StructGrid::ComputeMetrics2D(ActionKey *actkey)
{
    RDouble4D &xfn  = * (this->GetFaceNormalX());
    RDouble4D &yfn  = * (this->GetFaceNormalY());
    RDouble4D &zfn  = * (this->GetFaceNormalZ());
    RDouble4D &area = * (this->GetFaceArea());
    RDouble3D &vol  = * (this->GetCellVolume());
    RDouble4D &xfv  = * (this->GetFaceVectorX());
    RDouble4D &yfv  = * (this->GetFaceVectorY());
    RDouble4D &zfv  = * (this->GetFaceVectorZ());

    RDouble3D &xx = * this->GetStructX();
    RDouble3D &yy = * this->GetStructY();

    zfn = 0.0;
    zfv = 0.0;

    int ni = this->GetNI();
    int nj = this->GetNJ();

    RDouble r1[3], r2[3], r3[3], r4[3];
    RDouble vkc[3], vet[3];
    RDouble vvol;
    int i, j;

    RDouble volmin, volmax;
    volmin = LARGE;
    volmax = - LARGE;
    int imin, jmin, imax, jmax;
    imin = 1;
    jmin = 1;
    imax = 1;
    jmax = 1;

    std::ostringstream oss;

    int k = 1;
    int negativecell = 0;
    for (j = 1; j <= nj - 1; ++ j)
    {
        for (i = 1; i <= ni - 1; ++ i)
        {
            r1[0] = xx(i  , j, k);
            r1[1] = yy(i  , j, k);
            r1[2] = 0.0;

            r2[0] = xx(i+1, j, k);
            r2[1] = yy(i+1, j, k);
            r2[2] = 0.0;

            r3[0] = xx(i+1, j+1, k);
            r3[1] = yy(i+1, j+1, k);
            r3[2] = 0.0;

            r4[0] = xx(i, j+1, k);
            r4[1] = yy(i, j+1, k);
            r4[2] = 0.0;

            UnitFaceNormal2D(r4, r1, vkc);
            xfn (i, j, k, 1) = - vkc[0];
            yfn (i, j, k, 1) = - vkc[1];
            area(i, j, k, 1) =   vkc[2];
            xfv (i, j, k, 1) = - vkc[0] * vkc[2];
            yfv (i, j, k, 1) = - vkc[1] * vkc[2];

            UnitFaceNormal2D(r1, r2, vet);
            xfn (i, j, k, 2) = - vet[0];
            yfn (i, j, k, 2) = - vet[1];
            area(i, j, k, 2) =   vet[2];
            xfv (i, j, k, 2) = - vet[0] * vet[2];
            yfv (i, j, k, 2) = - vet[1] * vet[2];

            Area2D(r1, r2, r3, r4, vvol);
            vol(i, j, k) = vvol;

            if (vvol <= 0.0)
            {
                vvol = - vvol;
                vol(i, j, k) = vvol;
                ++ negativecell;
                oss << " Warning: negative cell occur : " << "(i, j, k) = " << i << " " << j << " " << k << endl;
                if (PHMPI::GetCurrentProcessorID() == PHMPI::GetServerProcessorID())
                {
                    PrintToWindow(oss);
                }
            }

            if (vvol < volmin)
            {
                volmin = vvol;
                imin = i;
                jmin = j;
            }

            if (vvol > volmax)
            {
                volmax = vvol;
                imax = i;
                jmax = j;
            }
        }
    }

    int numberOfTotalCell = this->GetNTotalCell();
    if (negativecell == numberOfTotalCell)
    {
        //oss << cell << " cells have negative vols \n";
        PHMPI::FreeBlockData();
        cout << "ERROR: cells have negative vols, maybe the i j k normal direction has set error! " << endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    //GhostCell2D(vol, ni, nj, 1);
    GhostCell3D(vol, ni, nj, 1);

    //! Compute the grid derivatives on the face of i = ni (value range of i is 1 : ni).
    i = ni;
    for (j = 1; j <= nj - 1; ++ j)
    {
        r1[0] = xx(i, j, k);
        r1[1] = yy(i, j, k);
        r1[2] = 0.0;

        r4[0] = xx(i, j+1, k);
        r4[1] = yy(i, j+1, k);
        r4[2] = 0.0;

        UnitFaceNormal2D(r4, r1, vkc);

        xfn (i, j, k, 1) = - vkc[0];
        yfn (i, j, k, 1) = - vkc[1];
        area(i, j, k, 1) =   vkc[2];
        xfv (i, j, k, 1) = - vkc[0] * vkc[2];
        yfv (i, j, k, 1) = - vkc[1] * vkc[2];
    }

    //! Compute the grid derivatives on the face of j = nj (value range of j is 1 : nj).
    j = nj;
    for (i = 1; i <= ni - 1; ++ i)
    {
        r1[0] = xx(i, j, k);
        r1[1] = yy(i, j, k);
        r1[2] = 0.0;

        r2[0] = xx(i+1, j, k);
        r2[1] = yy(i+1, j, k);
        r2[2] = 0.0;

        UnitFaceNormal2D(r1, r2, vet);

        xfn (i, j, k, 2) = - vet[0];
        yfn (i, j, k, 2) = - vet[1];
        area(i, j, k, 2) =   vet[2];
        xfv (i, j, k, 2) = - vet[0] * vet[2];
        yfv (i, j, k, 2) = - vet[1] * vet[2];
    }
    GhostMetrics2D();
    ComputeCellCenter();

//    oss << "Block " << zIndex << "\n";
//    oss << "min volume is " << imin << " " << jmin << " " << ": " << volmin << "\n";
//    oss << "max volume is " << imax << " " << jmax << " " << ": " << volmax << "\n";

    //StreamToActionKey(actkey, oss);
}

void StructGrid::GhostMetrics2D()
{
    int ni = this->GetNI();
    int nj = this->GetNJ();

    RDouble4D &xfn  = * (this->GetFaceNormalX());
    RDouble4D &yfn  = * (this->GetFaceNormalY());
    RDouble4D &area = * (this->GetFaceArea());
    RDouble4D &xfv  = * (this->GetFaceVectorX());
    RDouble4D &yfv  = * (this->GetFaceVectorY());

    int k = 1;
    int i, j;
    RDouble nx, ny, ns, ons;
    //! Assign initial value to the virtual points, and introduce the Green formula to modify.
    for (j = 1; j <= nj; ++ j)
    {
        xfn (0, j, k, 2) = xfn (1, j, k, 2);
        yfn (0, j, k, 2) = yfn (1, j, k, 2);
        area(0, j, k, 2) = area(1, j, k, 2);
        xfv (0, j, k, 2) = xfv (1, j, k, 2);
        yfv (0, j, k, 2) = yfv (1, j, k, 2);

        xfn (ni, j, k, 2) = xfn (ni-1, j, k, 2);
        yfn (ni, j, k, 2) = yfn (ni-1, j, k, 2);
        area(ni, j, k, 2) = area(ni-1, j, k, 2);
        xfv (ni, j, k, 2) = xfv (ni-1, j, k, 2);
        yfv (ni, j, k, 2) = yfv (ni-1, j, k, 2);
    }

    for (i = 1; i <= ni; ++ i)
    {
        xfn (i, 0, k, 1) = xfn (i, 1, k, 1);
        yfn (i, 0, k, 1) = yfn (i, 1, k, 1);
        area(i, 0, k, 1) = area(i, 1, k, 1);
        xfv (i, 0, k, 1) = xfv (i, 1, k, 1);
        yfv (i, 0, k, 1) = yfv (i, 1, k, 1);

        xfn (i, nj, k, 1) = xfn (i, nj-1, k, 1);
        yfn (i, nj, k, 1) = yfn (i, nj-1, k, 1);
        area(i, nj, k, 1) = area(i, nj-1, k, 1);
        xfv (i, nj, k, 1) = xfv (i, nj-1, k, 1);
        yfv (i, nj, k, 1) = yfv (i, nj-1, k, 1);
    }

    //! Introduce the Green formula to modify.
    for (j = 1; j <= nj - 1; ++ j)
    {
        nx = xfn(1, j,   k, 1) * area(1, j,   k, 1)
           - xfn(0, j,   k, 2) * area(0, j,   k, 2)
           + xfn(0, j+1, k, 2) * area(0, j+1, k, 2);
        ny = yfn(1, j,   k, 1) * area(1, j,   k, 1)
           - yfn(0, j,   k, 2) * area(0, j,   k, 2)
           + yfn(0, j+1, k, 2) * area(0, j+1, k, 2);

        ns = sqrt(nx * nx + ny * ny);
        ons = 1.0 / (ns + TINY);

        xfn (0, j, k, 1) = nx * ons;
        yfn (0, j, k, 1) = ny * ons;
        area(0, j, k, 1) = ns;
        xfv (0, j, k, 1) = nx;
        yfv (0, j, k, 1) = ny;

        nx = xfn(ni, j,   k, 1) * area(ni, j,   k, 1)
           + xfn(ni, j,   k, 2) * area(ni, j,   k, 2)
           - xfn(ni, j+1, k, 2) * area(ni, j+1, k, 2);
        ny = yfn(ni, j,   k, 1) * area(ni, j,   k, 1)
           + yfn(ni, j,   k, 2) * area(ni, j,   k, 2)
           - yfn(ni, j+1, k, 2) * area(ni, j+1, k, 2);
        ns = sqrt(nx * nx + ny * ny);
        ons = 1.0 / (ns + TINY);

        xfn (ni+1, j, k, 1) = nx * ons;
        yfn (ni+1, j, k, 1) = ny * ons;
        area(ni+1, j, k, 1) = ns;
        xfv (ni+1, j, k, 1) = nx;
        yfv (ni+1, j, k, 1) = ny;
    }

    for (i = 1; i <= ni - 1; ++ i)
    {
        nx = xfn(i,   1, k, 2) * area(i,   1, k, 2)
           - xfn(i,   0, k, 1) * area(i,   0, k, 1)
           + xfn(i+1, 0, k, 1) * area(i+1, 0, k, 1);
        ny = yfn(i,   1, k, 2) * area(i,   1, k, 2)
           - yfn(i,   0, k, 1) * area(i,   0, k, 1)
           + yfn(i+1, 0, k, 1) * area(i+1, 0, k, 1);
        ns = sqrt(nx * nx + ny * ny);
        ons = 1.0 / (ns + TINY);

        xfn (i, 0, k, 2) = nx * ons;
        yfn (i, 0, k, 2) = ny * ons;
        area(i, 0, k, 2) = ns;
        xfv (i, 0, k, 2) = nx;
        yfv (i, 0, k, 2) = ny;

        nx = xfn(i,   nj, k, 2) * area(i,   nj, k, 2)
           + xfn(i,   nj, k, 1) * area(i,   nj, k, 1)
           - xfn(i+1, nj, k, 1) * area(i+1, nj, k, 1);
        ny = yfn(i,   nj, k, 2) * area(i,   nj, k, 2)
           + yfn(i,   nj, k, 1) * area(i,   nj, k, 1)
           - yfn(i+1, nj, k, 1) * area(i+1, nj, k, 1);
        ns = sqrt(nx * nx + ny * ny);
        ons = 1.0 / (ns + TINY);

        xfn (i, nj+1, k, 2) = nx * ons;
        yfn (i, nj+1, k, 2) = ny * ons;
        area(i, nj+1, k, 2) = ns;
        xfv (i, nj+1, k, 2) = nx;
        yfv (i, nj+1, k, 2) = ny;
    }
}

void StructGrid::ReadMovingGrid()
{
    fstream file;
    std::ostringstream oss;
    oss << "move_grid" << this->GetZoneID() << ".dat";
    string filename = oss.str();

    file.open(filename.c_str(), ios_base::in|ios_base::binary);

    RDouble *x = this->GetX();
    RDouble *y = this->GetY();
    RDouble *z = this->GetZ();
    RDouble3D &voln = * this->GetCellVolumeOld();

    int nTotalNode = this->GetNTotalNode();
    file.read(reinterpret_cast<char*>(x), nTotalNode * sizeof(RDouble));
    file.read(reinterpret_cast<char*>(y), nTotalNode * sizeof(RDouble));
    file.read(reinterpret_cast<char*>(z), nTotalNode * sizeof(RDouble));

    int ist, ied, jst, jed, kst, ked;
    this->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, 2);

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                file.read(reinterpret_cast<char*>(&voln(i, j, k)), sizeof(RDouble));
            }
        }
    }
    file.close();
}

LIB_EXPORT void StructGrid::CVGNorm(FieldProxy *q_1_proxy, FieldProxy *q_0_proxy, int neqn, RDouble &norm)
{
    int nTotalCell = this->GetNTotalCell();

    RDouble diff;

    RDouble4D &q_1 = q_1_proxy->GetField_STR();
    RDouble4D &q_0 = q_0_proxy->GetField_STR();

    int ist, ied, jst, jed, kst, ked;
    this->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    norm = zero;
    for (int m = 0; m < neqn; ++ m)
    {
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    diff = q_1(i, j, k, m) - q_0(i, j, k, m);
                    norm = norm + diff * diff;
                }
            }
        }
    }

    norm = sqrt(norm / (nTotalCell * neqn));
}

void StructGrid::UpdateVolold()
{
    int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");

    int isAle = GlobalDataBase::GetIntParaFromDB("codeOfAleModel");

    if (!isUnsteady || !isAle) return;
    
    RDouble3D &vol  = * this->GetCellVolume();
    RDouble3D &voln = * this->GetCellVolumeOld();

    voln = vol;
}

void StructGrid::GhostMetrics3D()
{
    //! Assign initial value to the virtual points.
    RDouble nx, ny, nz, ns, ons;

    int i, j, k;
    int ni = GetNI();
    int nj = GetNJ();
    int nk = GetNK();
    RDouble4D *xfn = GetFaceNormalX();
    RDouble4D *yfn = GetFaceNormalY();
    RDouble4D *zfn = GetFaceNormalZ();
    RDouble4D *area = GetFaceArea();
    RDouble4D *xfv = GetFaceVectorX();
    RDouble4D *yfv = GetFaceVectorY();
    RDouble4D *zfv = GetFaceVectorZ();
    for (k = 1; k <= nk; ++ k)
    {
        for (j = 1; j <= nj; ++ j)
        {
            (*xfn)(0 , j, k, 2) = (*xfn)(1   , j, k, 2);
            (*yfn)(0 , j, k, 2) = (*yfn)(1   , j, k, 2);
            (*zfn)(0 , j, k, 2) = (*zfn)(1   , j, k, 2);
            (*area)(0 , j, k, 2) = (*area)(1   , j, k, 2);
            (*xfv)(0 , j, k, 2) = (*xfv)(1   , j, k, 2);
            (*yfv)(0 , j, k, 2) = (*yfv)(1   , j, k, 2);
            (*zfv)(0 , j, k, 2) = (*zfv)(1   , j, k, 2);

            (*xfn)(ni, j, k, 2) = (*xfn)(ni-1, j, k, 2);
            (*yfn)(ni, j, k, 2) = (*yfn)(ni-1, j, k, 2);
            (*zfn)(ni, j, k, 2) = (*zfn)(ni-1, j, k, 2);
            (*area)(ni, j, k, 2) = (*area)(ni-1, j, k, 2);
            (*xfv)(ni, j, k, 2) = (*xfv)(ni-1, j, k, 2);
            (*yfv)(ni, j, k, 2) = (*yfv)(ni-1, j, k, 2);
            (*zfv)(ni, j, k, 2) = (*zfv)(ni-1, j, k, 2);

            (*xfn)(0 , j, k, 3) = (*xfn)(1   , j, k, 3);
            (*yfn)(0 , j, k, 3) = (*yfn)(1   , j, k, 3);
            (*zfn)(0 , j, k, 3) = (*zfn)(1   , j, k, 3);
            (*area)(0 , j, k, 3) = (*area)(1   , j, k, 3);
            (*xfv)(0 , j, k, 3) = (*xfv)(1   , j, k, 3);
            (*yfv)(0 , j, k, 3) = (*yfv)(1   , j, k, 3);
            (*zfv)(0 , j, k, 3) = (*zfv)(1   , j, k, 3);

            (*xfn)(ni, j, k, 3) = (*xfn)(ni-1, j, k, 3);
            (*yfn)(ni, j, k, 3) = (*yfn)(ni-1, j, k, 3);
            (*zfn)(ni, j, k, 3) = (*zfn)(ni-1, j, k, 3);
            (*area)(ni, j, k, 3) = (*area)(ni-1, j, k, 3);
            (*xfv)(ni, j, k, 3) = (*xfv)(ni-1, j, k, 3);
            (*yfv)(ni, j, k, 3) = (*yfv)(ni-1, j, k, 3);
            (*zfv)(ni, j, k, 3) = (*zfv)(ni-1, j, k, 3);
        }
    }

    for (i = 1; i <= ni; ++ i)
    {
        for (k = 1; k <= nk; ++ k)
        {
            (*xfn)(i, 0, k, 1) = (*xfn)(i, 1, k, 1);
            (*yfn)(i, 0, k, 1) = (*yfn)(i, 1, k, 1);
            (*zfn)(i, 0, k, 1) = (*zfn)(i, 1, k, 1);
            (*area)(i, 0, k, 1) = (*area)(i, 1, k, 1);
            (*xfv)(i, 0, k, 1) = (*xfv)(i, 1, k, 1);
            (*yfv)(i, 0, k, 1) = (*yfv)(i, 1, k, 1);
            (*zfv)(i, 0, k, 1) = (*zfv)(i, 1, k, 1);

            (*xfn)(i, nj, k, 1) = (*xfn)(i, nj-1, k, 1);
            (*yfn)(i, nj, k, 1) = (*yfn)(i, nj-1, k, 1);
            (*zfn)(i, nj, k, 1) = (*zfn)(i, nj-1, k, 1);
            (*area)(i, nj, k, 1) = (*area)(i, nj-1, k, 1);
            (*xfv)(i, nj, k, 1) = (*xfv)(i, nj-1, k, 1);
            (*yfv)(i, nj, k, 1) = (*yfv)(i, nj-1, k, 1);
            (*zfv)(i, nj, k, 1) = (*zfv)(i, nj-1, k, 1);

            (*xfn)(i, 0, k, 3) = (*xfn)(i, 1, k, 3);
            (*yfn)(i, 0, k, 3) = (*yfn)(i, 1, k, 3);
            (*zfn)(i, 0, k, 3) = (*zfn)(i, 1, k, 3);
            (*area)(i, 0, k, 3) = (*area)(i, 1, k, 3);
            (*xfv)(i, 0, k, 3) = (*xfv)(i, 1, k, 3);
            (*yfv)(i, 0, k, 3) = (*yfv)(i, 1, k, 3);
            (*zfv)(i, 0, k, 3) = (*zfv)(i, 1, k, 3);

            (*xfn)(i, nj, k, 3) = (*xfn)(i, nj-1, k, 3);
            (*yfn)(i, nj, k, 3) = (*yfn)(i, nj-1, k, 3);
            (*zfn)(i, nj, k, 3) = (*zfn)(i, nj-1, k, 3);
            (*area)(i, nj, k, 3) = (*area)(i, nj-1, k, 3);
            (*xfv)(i, nj, k, 3) = (*xfv)(i, nj-1, k, 3);
            (*yfv)(i, nj, k, 3) = (*yfv)(i, nj-1, k, 3);
            (*zfv)(i, nj, k, 3) = (*zfv)(i, nj-1, k, 3);
        }
    }

    for (j = 1; j <= nj; ++ j)
    {
        for (i = 1; i <= ni; ++ i)
        {
            (*xfn)(i, j, 0, 1) = (*xfn)(i, j, 1, 1);
            (*yfn)(i, j, 0, 1) = (*yfn)(i, j, 1, 1);
            (*zfn)(i, j, 0, 1) = (*zfn)(i, j, 1, 1);
            (*area)(i, j, 0, 1) = (*area)(i, j, 1, 1);
            (*xfv)(i, j, 0, 1) = (*xfv)(i, j, 1, 1);
            (*yfv)(i, j, 0, 1) = (*yfv)(i, j, 1, 1);
            (*zfv)(i, j, 0, 1) = (*zfv)(i, j, 1, 1);

            (*xfn)(i, j, nk, 1) = (*xfn)(i, j, nk-1, 1);
            (*yfn)(i, j, nk, 1) = (*yfn)(i, j, nk-1, 1);
            (*zfn)(i, j, nk, 1) = (*zfn)(i, j, nk-1, 1);
            (*area)(i, j, nk, 1) = (*area)(i, j, nk-1, 1);
            (*xfv)(i, j, nk, 1) = (*xfv)(i, j, nk-1, 1);
            (*yfv)(i, j, nk, 1) = (*yfv)(i, j, nk-1, 1);
            (*zfv)(i, j, nk, 1) = (*zfv)(i, j, nk-1, 1);

            (*xfn)(i, j, 0, 2) = (*xfn)(i, j, 1, 2);
            (*yfn)(i, j, 0, 2) = (*yfn)(i, j, 1, 2);
            (*zfn)(i, j, 0, 2) = (*zfn)(i, j, 1, 2);
            (*area)(i, j, 0, 2) = (*area)(i, j, 1, 2);
            (*xfv)(i, j, 0, 2) = (*xfv)(i, j, 1, 2);
            (*yfv)(i, j, 0, 2) = (*yfv)(i, j, 1, 2);
            (*zfv)(i, j, 0, 2) = (*zfv)(i, j, 1, 2);

            (*xfn)(i, j, nk, 2) = (*xfn)(i, j, nk-1, 2);
            (*yfn)(i, j, nk, 2) = (*yfn)(i, j, nk-1, 2);
            (*zfn)(i, j, nk, 2) = (*zfn)(i, j, nk-1, 2);
            (*area)(i, j, nk, 2) = (*area)(i, j, nk-1, 2);
            (*xfv)(i, j, nk, 2) = (*xfv)(i, j, nk-1, 2);
            (*yfv)(i, j, nk, 2) = (*yfv)(i, j, nk-1, 2);
            (*zfv)(i, j, nk, 2) = (*zfv)(i, j, nk-1, 2);
        }
    }

    //! Introduce the Green formula to modify.
    for (k = 1; k < nk; ++ k)
    {
        for (j = 1; j < nj; ++ j)
        {
            nx = (*xfn)(1, j  , k  , 1) * (*area)(1, j  , k  , 1)
               + (*xfn)(0, j+1, k  , 2) * (*area)(0, j+1, k  , 2)
               - (*xfn)(0, j  , k  , 2) * (*area)(0, j  , k  , 2)
               + (*xfn)(0, j  , k+1, 3) * (*area)(0, j  , k+1, 3)
               - (*xfn)(0, j  , k  , 3) * (*area)(0, j  , k  , 3);

            ny = (*yfn)(1, j  , k  , 1) * (*area)(1, j  , k  , 1)
               + (*yfn)(0, j+1, k  , 2) * (*area)(0, j+1, k  , 2)
               - (*yfn)(0, j  , k  , 2) * (*area)(0, j  , k  , 2)
               + (*yfn)(0, j  , k+1, 3) * (*area)(0, j  , k+1, 3)
               - (*yfn)(0, j  , k  , 3) * (*area)(0, j  , k  , 3);

            nz = (*zfn)(1, j  , k  , 1) * (*area)(1, j  , k  , 1)
               + (*zfn)(0, j+1, k  , 2) * (*area)(0, j+1, k  , 2)
               - (*zfn)(0, j  , k  , 2) * (*area)(0, j  , k  , 2)
               + (*zfn)(0, j  , k+1, 3) * (*area)(0, j  , k+1, 3)
               - (*zfn)(0, j  , k  , 3) * (*area)(0, j  , k  , 3);

            ns  = sqrt(nx * nx + ny * ny + nz * nz);
            ons = 1.0 / (ns + TINY);

            (*xfn)(0, j, k, 1) = nx * ons;
            (*yfn)(0, j, k, 1) = ny * ons;
            (*zfn)(0, j, k, 1) = nz * ons;
            (*area)(0, j, k, 1) = ns;
            (*xfv)(0, j, k, 1) = nx;
            (*yfv)(0, j, k, 1) = ny;
            (*zfv)(0, j, k, 1) = nz;

            nx = (*xfn)(ni, j  , k  , 1) * (*area)(ni, j  , k  , 1)
               - (*xfn)(ni, j+1, k  , 2) * (*area)(ni, j+1, k  , 2)
               + (*xfn)(ni, j  , k  , 2) * (*area)(ni, j  , k  , 2)
               - (*xfn)(ni, j  , k+1, 3) * (*area)(ni, j  , k+1, 3)
               + (*xfn)(ni, j  , k  , 3) * (*area)(ni, j  , k  , 3);

            ny = (*yfn)(ni, j  , k  , 1) * (*area)(ni, j  , k  , 1)
               - (*yfn)(ni, j+1, k  , 2) * (*area)(ni, j+1, k  , 2)
               + (*yfn)(ni, j  , k  , 2) * (*area)(ni, j  , k  , 2)
               - (*yfn)(ni, j  , k+1, 3) * (*area)(ni, j  , k+1, 3)
               + (*yfn)(ni, j  , k  , 3) * (*area)(ni, j  , k  , 3);

            nz = (*zfn)(ni, j  , k  , 1) * (*area)(ni, j  , k  , 1)
               - (*zfn)(ni, j+1, k  , 2) * (*area)(ni, j+1, k  , 2)
               + (*zfn)(ni, j  , k  , 2) * (*area)(ni, j  , k  , 2)
               - (*zfn)(ni, j  , k+1, 3) * (*area)(ni, j  , k+1, 3)
               + (*zfn)(ni, j  , k  , 3) * (*area)(ni, j  , k  , 3);

            ns  = sqrt(nx * nx + ny * ny + nz * nz);
            ons = 1.0 / (ns + TINY);

            (*xfn)(ni+1, j, k, 1) = nx * ons;
            (*yfn)(ni+1, j, k, 1) = ny * ons;
            (*zfn)(ni+1, j, k, 1) = nz * ons;
            (*area)(ni+1, j, k, 1) = ns;
            (*xfv)(ni+1, j, k, 1) = nx;
            (*yfv)(ni+1, j, k, 1) = ny;
            (*zfv)(ni+1, j, k, 1) = nz;
        }
    }

    for (i = 1; i <= ni; ++ i)
    {
        for (k = 1; k <= nk; ++ k)
        {
            nx = (*xfn)(i  , 1, k  , 2) * (*area)(i  , 1, k  , 2)
               + (*xfn)(i  , 0, k+1, 3) * (*area)(i  , 0, k+1, 3)
               - (*xfn)(i  , 0, k  , 3) * (*area)(i  , 0, k  , 3)
               + (*xfn)(i+1, 0, k  , 1) * (*area)(i+1, 0, k  , 1)
               - (*xfn)(i  , 0, k  , 1) * (*area)(i  , 0, k  , 1);

            ny = (*yfn)(i  , 1, k  , 2) * (*area)(i  , 1, k  , 2)
               + (*yfn)(i  , 0, k+1, 3) * (*area)(i  , 0, k+1, 3)
               - (*yfn)(i  , 0, k  , 3) * (*area)(i  , 0, k  , 3)
               + (*yfn)(i+1, 0, k  , 1) * (*area)(i+1, 0, k  , 1)
               - (*yfn)(i  , 0, k  , 1) * (*area)(i  , 0, k  , 1);

            nz = (*zfn)(i  , 1, k  , 2) * (*area)(i  , 1, k  , 2)
               + (*zfn)(i  , 0, k+1, 3) * (*area)(i  , 0, k+1, 3)
               - (*zfn)(i  , 0, k  , 3) * (*area)(i  , 0, k  , 3)
               + (*zfn)(i+1, 0, k  , 1) * (*area)(i+1, 0, k  , 1)
               - (*zfn)(i  , 0, k  , 1) * (*area)(i  , 0, k  , 1);

            ns  = sqrt(nx * nx + ny * ny + nz * nz);
            ons = 1.0 / (ns + TINY);

            (*xfn)(i, 0, k, 2) = nx * ons;
            (*yfn)(i, 0, k, 2) = ny * ons;
            (*zfn)(i, 0, k, 2) = nz * ons;
            (*area)(i, 0, k, 2) = ns;
            (*xfv)(i, 0, k, 2) = nx;
            (*yfv)(i, 0, k, 2) = ny;
            (*zfv)(i, 0, k, 2) = nz;

            nx = (*xfn)(i  , nj, k  , 2) * (*area)(i  , nj, k  , 2)
               - (*xfn)(i  , nj, k+1, 3) * (*area)(i  , nj, k+1, 3)
               + (*xfn)(i  , nj, k  , 3) * (*area)(i  , nj, k  , 3)
               - (*xfn)(i+1, nj, k  , 1) * (*area)(i+1, nj, k  , 1)
               + (*xfn)(i  , nj, k  , 1) * (*area)(i  , nj, k  , 1);

            ny = (*yfn)(i  , nj, k  , 2) * (*area)(i  , nj, k  , 2)
               - (*yfn)(i  , nj, k+1, 3) * (*area)(i  , nj, k+1, 3)
               + (*yfn)(i  , nj, k  , 3) * (*area)(i  , nj, k  , 3)
               - (*yfn)(i+1, nj, k  , 1) * (*area)(i+1, nj, k  , 1)
               + (*yfn)(i  , nj, k  , 1) * (*area)(i  , nj, k  , 1);

            nz = (*zfn)(i  , nj, k  , 2) * (*area)(i  , nj, k  , 2)
               - (*zfn)(i  , nj, k+1, 3) * (*area)(i  , nj, k+1, 3)
               + (*zfn)(i  , nj, k  , 3) * (*area)(i  , nj, k  , 3)
               - (*zfn)(i+1, nj, k  , 1) * (*area)(i+1, nj, k  , 1)
               + (*zfn)(i  , nj, k  , 1) * (*area)(i  , nj, k  , 1);

            ns  = sqrt(nx * nx + ny * ny + nz * nz);
            ons = 1.0 / (ns + TINY);

            (*xfn)(i, nj+1, k, 2) = nx * ons;
            (*yfn)(i, nj+1, k, 2) = ny * ons;
            (*zfn)(i, nj+1, k, 2) = nz * ons;
            (*area)(i, nj+1, k, 2) = ns;
            (*xfv)(i, nj+1, k, 2) = nx;
            (*yfv)(i, nj+1, k, 2) = ny;
            (*zfv)(i, nj+1, k, 2) = nz;
        }
    }

    for (j = 1; j <= nj; ++ j)
    {
        for (i = 1; i <= ni; ++ i)
        {
            nx = (*xfn)(i  , j  , 1, 3) * (*area)(i  , j  , 1, 3)
               + (*xfn)(i+1, j  , 0, 1) * (*area)(i+1, j  , 0, 1)
               - (*xfn)(i  , j  , 0, 1) * (*area)(i  , j  , 0, 1)
               + (*xfn)(i  , j+1, 0, 2) * (*area)(i  , j+1, 0, 2)
               - (*xfn)(i  , j  , 0, 2) * (*area)(i  , j  , 0, 2);

            ny = (*yfn)(i  , j  , 1, 3) * (*area)(i  , j  , 1, 3)
               + (*yfn)(i+1, j  , 0, 1) * (*area)(i+1, j  , 0, 1)
               - (*yfn)(i  , j  , 0, 1) * (*area)(i  , j  , 0, 1)
               + (*yfn)(i  , j+1, 0, 2) * (*area)(i  , j+1, 0, 2)
               - (*yfn)(i  , j  , 0, 2) * (*area)(i  , j  , 0, 2);

            nz = (*zfn)(i  , j  , 1, 3) * (*area)(i  , j  , 1, 3)
               + (*zfn)(i+1, j  , 0, 1) * (*area)(i+1, j  , 0, 1)
               - (*zfn)(i  , j  , 0, 1) * (*area)(i  , j  , 0, 1)
               + (*zfn)(i  , j+1, 0, 2) * (*area)(i  , j+1, 0, 2)
               - (*zfn)(i  , j  , 0, 2) * (*area)(i  , j  , 0, 2);

            ns = sqrt(nx * nx + ny * ny + nz * nz);
            ons = 1.0 / (ns + TINY);

            (*xfn)(i, j, 0, 3) = nx * ons;
            (*yfn)(i, j, 0, 3) = ny * ons;
            (*zfn)(i, j, 0, 3) = nz * ons;
            (*area)(i, j, 0, 3) = ns;
            (*xfv)(i, j, 0, 3) = nx;
            (*yfv)(i, j, 0, 3) = ny;
            (*zfv)(i, j, 0, 3) = nz;

            nx = (*xfn)(i  , j  , nk, 3) * (*area)(i  , j  , nk, 3)
               - (*xfn)(i+1, j  , nk, 1) * (*area)(i+1, j  , nk, 1)
               + (*xfn)(i  , j  , nk, 1) * (*area)(i  , j  , nk, 1)
               - (*xfn)(i  , j+1, nk, 2) * (*area)(i  , j+1, nk, 2)
               + (*xfn)(i  , j  , nk, 2) * (*area)(i  , j  , nk, 2);

            ny = (*yfn)(i  , j  , nk, 3) * (*area)(i  , j  , nk, 3)
               - (*yfn)(i+1, j  , nk, 1) * (*area)(i+1, j  , nk, 1)
               + (*yfn)(i  , j  , nk, 1) * (*area)(i  , j  , nk, 1)
               - (*yfn)(i  , j+1, nk, 2) * (*area)(i  , j+1, nk, 2)
               + (*yfn)(i  , j  , nk, 2) * (*area)(i  , j  , nk, 2);

            nz = (*zfn)(i  , j  , nk, 3) * (*area)(i  , j  , nk, 3)
               - (*zfn)(i+1, j  , nk, 1) * (*area)(i+1, j  , nk, 1)
               + (*zfn)(i  , j  , nk, 1) * (*area)(i  , j  , nk, 1)
               - (*zfn)(i  , j+1, nk, 2) * (*area)(i  , j+1, nk, 2)
               + (*zfn)(i  , j  , nk, 2) * (*area)(i  , j  , nk, 2);

            ns  = sqrt(nx * nx + ny * ny + nz * nz);
            ons = 1.0 / (ns + TINY);

            (*xfn)(i, j, nk+1, 3) = nx * ons;
            (*yfn)(i, j, nk+1, 3) = ny * ons;
            (*zfn)(i, j, nk+1, 3) = nz * ons;
            (*area)(i, j, nk+1, 3) = ns;
            (*xfv)(i, j, nk+1, 3) = nx;
            (*yfv)(i, j, nk+1, 3) = ny;
            (*zfv)(i, j, nk+1, 3) = nz;
        }
    }
}

LIB_EXPORT void StructGrid::SetLnkInfo()
{
    this->composite_bcregion->SetLnkInfo();
}

LIB_EXPORT void StructGrid::ProcessBCInfo()
{
    this->composite_bcregion->ProcessBCInfo();
    this->SetBCFaceInfo();
}

LIB_EXPORT void StructGrid::CoarseGrids(int nlmax)
{
    //++++++++++++++++++++++++++++++++++++++++++++++
    //+ The finest grid level = 0
    //+ The coarse grid level = 1, 2, 3, ...
    //++++++++++++++++++++++++++++++++++++++++++++++

    int iLevel = this->GetLevel();
    if (iLevel > nlmax)
    {
        this->SetCoarseGrid(0);
        return;
    }

    int ni = this->GetNI();
    int nj = this->GetNJ();
    int nk = this->GetNK();

    //++++++++++++++++++++++++++++++++++++++++++++++
    //+ Take a look at the rules: ni = 3
    //+                           ni/2 + 1 = 2
    //+                           ni/2 + 1 = 2
    //+                           ni = 5
    //+                           ni/2 + 1 = 3
    //+                           ni/2 + 1 = 2
    //+                           ni/2 + 1 = 2
    //+                           ni = 11
    //+                           ni/2 + 1 = 6
    //+                           ni/2 + 1 = 4
    //++++++++++++++++++++++++++++++++++++++++++++++

    //! mgconst(6) = 1 standard multigrid - window and grid sizes need
    //                 to be a multiple of 2.
    //               2 heterogeneous multigrid.
    //               3 anisotropic multigrid.

    int mgtype = 2;

    int istp, jstp, kstp;
    istp = 2;
    jstp = 2;
    kstp = 2;

    if (MAX((ni - 1), 1) % 2 != 0)
    {
        if (mgtype <= 1) return;
        istp = 1;
    }

    if (MAX((nj - 1), 1) % 2 != 0)
    {
        if (mgtype <= 1) return;
        jstp = 1;
    }

    if (MAX((nk - 1), 1) % 2 != 0)
    {
        if (mgtype <= 1) return;
        kstp = 1;
    }

    int imin, jmin, kmin, imax, jmax, kmax;

    CompositeStructBC *compositeBCRegion = this->GetCompositeBCRegion();
    int nBCRegion = compositeBCRegion->GetnBCRegion();

    for (int ibcregion = 0; ibcregion < nBCRegion; ++ ibcregion)
    {
        StructBC *bcregion = compositeBCRegion->GetBCRegion(ibcregion);

        bcregion->GetIJKRegion(imin, imax, jmin, jmax, kmin, kmax);

        int nsurf = bcregion->GetFaceDirection() + 1;

        if (nsurf == 1)
        {
            if ((jmin - 1) % 2 != 0)
            {
                if (mgtype <= 1) return;
                jstp = 1;
            }
            if (jmax % 2 != 0)
            {
                if (mgtype <= 1) return;
                jstp = 1;
            }

            if ((kmin - 1) % 2 != 0)
            {
                if (mgtype <= 1) return;
                kstp = 1;
            }
            if (kmax % 2 != 0)
            {
                if (mgtype <= 1) return;
                kstp = 1;
            }
        }
        else if (nsurf == 2)
        {
            if ((kmin - 1) % 2 != 0)
            {
                if (mgtype <= 1) return;
                kstp = 1;
            }
            if (kmax % 2 != 0)
            {
                if (mgtype <= 1) return;
                kstp = 1;
            }

            if ((imin - 1) % 2 != 0)
            {
                if (mgtype <= 1) return;
                istp = 1;
            }
            if (imax % 2 != 0)
            {
                if (mgtype <= 1) return;
                istp = 1;
            }
        }
        else if (nsurf == 3)
        {
            if ((imin - 1) % 2 != 0)
            {
                if (mgtype <= 1) return;
                istp = 1;
            }
            if (imax % 2 != 0)
            {
                if (mgtype <= 1) return;
                istp = 1;
            }

            if ((jmin - 1) % 2 != 0)
            {
                if (mgtype <= 1) return;
                jstp = 1;
            }
            if (jmax % 2 != 0)
            {
                if (mgtype <= 1) return;
                jstp = 1;
            }
        }
        else
        {
            PHMPI::FreeBlockData();
            cout << "FATAL error in void StructGrid::CoarseGrids(int nlmax)!\n" << endl;
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
    }

    this->SetMultigridStepI(istp);
    this->SetMultigridStepJ(jstp);
    this->SetMultigridStepK(kstp);

    //! After the previous process, the program runs here, it must meet the needs.
    int cni, cnj, cnk;
    cni = (ni - 1) / istp + 1;
    cnj = (nj - 1) / jstp + 1;
    cnk = (nk - 1) / kstp + 1;

    StructGrid *cgrid = new StructGrid();
    cgrid->InitGrid(new GridID(* this->GetGridID()), level + 1, this->GetDim(), STRUCTGRID);

    //! coor1 = Coordinates of the corner points of the control
    //!         volumes of the fine grid.
    //! n11   = Number of cells in I-direction of the fine grid.
    //! n21   = Number of cells in J-direction of the fine grid.
    //! n31   = Number of cells in K-direction of the fine grid.
    //! n10   = Number of cells in I-direction of the coarse grid.
    //! n20   = Number of cells in J-direction of the coarse grid.
    //! n30   = Number of cells in K-direction of the coarse grid.

    int cnTNode = cni * cnj * cnk;
    Range cIC(1, cni-1);
    Range cJC(1, cnj-1);
    Range cKC(1, cnk-1);

    if (cgrid->GetDim() == TWO_D) cKC.setRange(1, 1);
    int cnTCell = cIC.length() * cJC.length() * cKC.length();
    int cnTFace = (cnj - 1) * (cnk - 1) * cni + (cnk - 1) * (cni - 1) * cnj + (cni - 1) * (cnj - 1) * cnk;

    cgrid->SetNTotalNode(cnTNode);
    cgrid->SetNTotalFace(cnTFace);
    cgrid->SetNTotalCell(cnTCell);

    RDouble *cx = new RDouble [cnTNode];
    RDouble *cy = new RDouble [cnTNode];
    RDouble *cz = new RDouble [cnTNode];

    cgrid->SetX(cx);
    cgrid->SetY(cy);
    cgrid->SetZ(cz);

    cgrid->SetNI(cni);
    cgrid->SetNJ(cnj);
    cgrid->SetNK(cnk);

    cgrid->SetArrayLayout();
    RDouble3D &xx = * (this->GetStructX());
    RDouble3D &yy = * (this->GetStructY());
    RDouble3D &zz = * (this->GetStructZ());

    RDouble3D &cxx = * (cgrid->GetStructX());
    RDouble3D &cyy = * (cgrid->GetStructY());
    RDouble3D &czz = * (cgrid->GetStructZ());

    int i0, j0, k0;

    k0 = 1;
    for (int k = 1; k <= nk; k += kstp)
    {
        j0 = 1;
        for (int j = 1; j <= nj; j += jstp)
        {
            i0 = 1;
            for (int i = 1; i <= ni; i += istp)
            {
                cxx(i0, j0, k0) = xx(i, j, k);
                cyy(i0, j0, k0) = yy(i, j, k);
                czz(i0, j0, k0) = zz(i, j, k);
                i0 += 1;
            }
            j0 += 1;
        }
        k0 += 1;
    }

    int nb = this->GetZoneID();

    cgrid->CreateCompositeBCRegion(nBCRegion);
    CompositeStructBC *c_composite_bcregion = cgrid->GetCompositeBCRegion();

    int bctype;
    string bcname;
    int s_lr, s_nd;

    for (int ibcregion = 0; ibcregion < nBCRegion; ++ ibcregion)
    {
        StructBC *bcregion = composite_bcregion->GetBCRegion(ibcregion);

        bcregion->GetIJKRegion(imin, imax, jmin, jmax, kmin, kmax);
        bctype = bcregion->GetBCType();
        bcname = bcregion->GetBCName();
        s_lr   = bcregion->GetFaceLeftOrRightIndex();
        s_nd   = bcregion->GetFaceDirection();

        int nsurf = bcregion->GetFaceDirection() + 1;

        StructBC *cbcregion = new StructBC(nb, ibcregion);
        c_composite_bcregion->SetBCRegion(ibcregion, cbcregion);

        cbcregion->ComputeMultiGridIJKRegion(imin, imax, jmin, jmax, kmin, kmax, istp, jstp, kstp, nsurf);
        cbcregion->SetBCType(bctype);
        cbcregion->SetBCName(bcname);
        cbcregion->SetFaceLeftOrRightIndex(s_lr);
        cbcregion->SetFaceDirection(s_nd);
        cbcregion->InitFaceDirectionIndex();
    }

    cgrid->SetBCFaceInfo();

    this->SetCoarseGrid(cgrid);
    cgrid->SetFineGrid(this);
    cgrid->SetLevel(this->GetLevel() + 1);

    cgrid->ComputeMetrics();
}

LIB_EXPORT void StructGrid::ReadXYZ(fstream &file)
{
    RDouble3D &xx = * this->GetStructX();
    RDouble3D &yy = * this->GetStructY();
    RDouble3D &zz = * this->GetStructZ();

    int ist, ied, jst, jed, kst, ked;
    this->GetNodeIterationIndex(ist, ied, jst, jed, kst, ked);

    int fileFormat = GlobalDataBase::GetIntParaFromDB("fileformat");
    //! fileformat : Grid file format.
    //!              0 -- BINARY
    //!              1 -- ASCII
    if (fileFormat == 1)
    {
        ReadXYZASCII(file, xx, ist, ied, jst, jed, kst, ked);
        ReadXYZASCII(file, yy, ist, ied, jst, jed, kst, ked);
        ReadXYZASCII(file, zz, ist, ied, jst, jed, kst, ked);
    }
    else if (fileFormat == 0)
    {
        ReadXYZ(file, xx, ist, ied, jst, jed, kst, ked);
        ReadXYZ(file, yy, ist, ied, jst, jed, kst, ked);
        ReadXYZ(file, zz, ist, ied, jst, jed, kst, ked);
    }

    ComputeMinMaxBox();
}

void StructGrid::ReadXYZ(fstream &file, RDouble3D &coor, int ist, int ied, int jst, int jed, int kst, int ked)
{
    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                file.read(reinterpret_cast<char *>(&coor(i, j, k)), sizeof(RDouble));
            }
        }
    }
}

void StructGrid::ReadXYZASCII(fstream &file, RDouble3D &coor, int ist, int ied, int jst, int jed, int kst, int ked)
{
    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                file >> coor(i, j, k);
            }
        }
    }
}

void StructGrid::WriteXYZ(fstream &file, RDouble3D &coor, int ist, int ied, int jst, int jed, int kst, int ked)
{
    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                file.write(reinterpret_cast<char *>(&coor(i, j, k)), sizeof(RDouble));
            }
        }
    }
}

LIB_EXPORT void StructGrid::WriteXYZ(fstream &file, int ist, int ied, int jst, int jed, int kst, int ked)
{
    RDouble3D &xx = * this->GetStructX();
    RDouble3D &yy = * this->GetStructY();
    RDouble3D &zz = * this->GetStructZ();

    WriteXYZ(file, xx, ist, ied, jst, jed, kst, ked);
    WriteXYZ(file, yy, ist, ied, jst, jed, kst, ked);
    WriteXYZ(file, zz, ist, ied, jst, jed, kst, ked);
}

LIB_EXPORT void StructGrid::SetBasicDimension()
{
    int ni = this->GetNI();
    int nj = this->GetNJ();
    int nk = this->GetNK();

    Range IC, JC, KC;
    GetRange(ni, nj, nk, 0, -1, IC, JC, KC);

    int nTotalNode, nTotalFace, nTotalCell;

    nTotalNode = ni * nj * nk;
    nTotalCell = IC.length() * JC.length() * KC.length();

    nTotalFace = (nj - 1) * (nk - 1) * ni + (nk - 1) * (ni - 1) * nj + (ni - 1) * (nj - 1) * nk;

    if (this->GetDim() == TWO_D)
    {
        nTotalFace = (nj - 1) * ni + (ni - 1) * nj;
    }

    this->SetNTotalNode(nTotalNode);
    this->SetNTotalFace(nTotalFace);
    this->SetNTotalCell(nTotalCell);

    if (PHMPI::GetCurrentProcessorID() == 0)
    {
        ostringstream oss;
        oss << "  Grid Dimension  : " << " " << this->GetDim() << " " << "\n";
        oss << "  Number of Points: " << " " << nTotalNode << " " << "\n";
        oss << "  Number of Faces : " << " " << nTotalFace << " " << "\n";
        oss << "  Number of Cells : " << " " << nTotalCell << " " << "\n";
        cout << oss.str();
    }
}

LIB_EXPORT void StructGrid::SetBasicDimension(int ni, int nj, int nk)
{
    this->SetNI(ni);
    this->SetNJ(nj);
    this->SetNK(nk);
    Range I(1, ni - 1), J(1, nj - 1), K(1, nk - 1);

    nTotalNode = ni * nj * nk;
    nTotalCell = I.length() * J.length() * K.length();
    nBoundFace = (ni - 1) * (nj - 1) * nk + (ni - 1) * nj * (nk - 1) + ni * (nj - 1) * (nk - 1);

    return;
}

LIB_EXPORT void StructGrid::ReadGrid(fstream &file)
{
    VirtualFile *vfile = new VirtualFile(&file);

    vfile->BeginReadWork();

    ReadGrid(vfile);

    vfile->EndReadWork();

    delete vfile;
}

LIB_EXPORT void StructGrid::WriteGrid(fstream &file)
{
    VirtualFile *vfile = new VirtualFile(&file);

    vfile->BeginWriteWork();

    WriteGrid(vfile);

    vfile->EndWriteWork();

    delete vfile;
}

void StructGrid::DecodeGrid(DataContainer *cdata)
{
    VirtualFile *vfile = new VirtualFile(cdata);

    vfile->BeginReadWork();

    ReadGrid(vfile);

    vfile->EndReadWork();

    delete vfile;
}

void StructGrid::EncodeGrid(DataContainer *cdata)
{
    VirtualFile *vfile = new VirtualFile(cdata);

    vfile->BeginWriteWork();

    WriteGrid(vfile);

    vfile->EndWriteWork();

    delete vfile;
}

void StructGrid::ReadGrid(VirtualFile *vfile)
{
    if (PHMPI::GetCurrentProcessorID() == 0)
    {
        ostringstream oss;
        oss << "Reading Structured Grid Of Zone " << " " << GetZoneID() << " " << "...\n";
        cout << oss.str();
    }
    
    int ni = 0;
    int nj = 0;
    int nk = 0;

    vfile->read(reinterpret_cast<char*>(&ni), sizeof(int));
    vfile->read(reinterpret_cast<char*>(&nj), sizeof(int));
    vfile->read(reinterpret_cast<char*>(&nk), sizeof(int));
    SetNI(ni);
    SetNJ(nj);
    SetNK(nk);
    SetBasicDimension();

    int nTotalNode = this->GetNTotalNode();

    RDouble *x, *y, *z;
    x = new RDouble [nTotalNode];
    y = new RDouble [nTotalNode];
    z = new RDouble [nTotalNode];

    this->SetX(x);
    this->SetY(y);
    this->SetZ(z);

    vfile->read(reinterpret_cast<char*>(x), nTotalNode * sizeof(RDouble));
    vfile->read(reinterpret_cast<char*>(y), nTotalNode * sizeof(RDouble));
    vfile->read(reinterpret_cast<char*>(z), nTotalNode * sizeof(RDouble));

    int axisup = 1;
    PHSPACE::GlobalDataBase::GetData("axisup", &axisup, PHINT, 1);
    if (axisup == 2)
    {
        TurnZAxis2YAxis();
    }

    ComputeMinMaxBox();
    SetArrayLayout();

    int sizeint = sizeof(int);
    int nBCRegion;
    vfile->read(reinterpret_cast<char*>(&nBCRegion), sizeint);

    this->CreateCompositeBCRegion(nBCRegion);
    CompositeStructBC *composite_bcregion = this->GetCompositeBCRegion();

    int nb = this->GetZoneID();
    int imin, imax, jmin, jmax, kmin, kmax, bctype;
    int s_lr, s_nd;
    //! Read the source region imin, imax, jmin, jmax, kmin, kmax, bctype.
    for (int ibcregion = 0; ibcregion < nBCRegion; ++ ibcregion)
    {
        StructBC *bcregion = new StructBC(nb, ibcregion);
        composite_bcregion->SetBCRegion(ibcregion, bcregion);
        vfile->read(reinterpret_cast<char*>(&imin), sizeint);
        vfile->read(reinterpret_cast<char*>(&imax), sizeint);
        vfile->read(reinterpret_cast<char*>(&jmin), sizeint);
        vfile->read(reinterpret_cast<char*>(&jmax), sizeint);
        vfile->read(reinterpret_cast<char*>(&kmin), sizeint);
        vfile->read(reinterpret_cast<char*>(&kmax), sizeint);

        vfile->read(reinterpret_cast<char*>(&bctype), sizeint);
        vfile->read(reinterpret_cast<char*>(&s_lr), sizeint);
        vfile->read(reinterpret_cast<char*>(&s_nd), sizeint);

        bcregion->SetIJKRegion(imin, imax, jmin, jmax, kmin, kmax);
        bcregion->SetBCType(bctype);
        bcregion->SetFaceLeftOrRightIndex(s_lr);
        bcregion->SetFaceDirection(s_nd);
        bcregion->InitFaceDirectionIndex();
    }
    this->SetBCFaceInfo();

    int nIFace;
    vfile->read(reinterpret_cast<char*>(&nIFace), sizeint);
    InterfaceInfo *interfaceInfo = 0;
    this->SetInterfaceInfo(interfaceInfo);

    if (nIFace > 0)
    {
        interfaceInfo = new InterfaceInfo(nIFace, this);
        this->SetInterfaceInfo(interfaceInfo);
        int *interFace2ZoneID       = interfaceInfo->GetInterFace2ZoneID();
        int *interFace2InterFaceID  = interfaceInfo->GetInterFace2InterFaceID();
        int *interFace2BoundaryFace = interfaceInfo->GetInterFace2BoundaryFace();

        vfile->read(reinterpret_cast<char*>(interFace2ZoneID     ), nIFace * sizeof(int));
        vfile->read(reinterpret_cast<char*>(interFace2InterFaceID), nIFace * sizeof(int));
        vfile->read(reinterpret_cast<char*>(interFace2BoundaryFace), nIFace * sizeof(int));
    }
}

void StructGrid::WriteGrid(VirtualFile *vfile)
{
    int ni = GetNI();
    int nj = GetNJ();
    int nk = GetNK();
    vfile->write(reinterpret_cast<char*>(&ni), sizeof(int));
    vfile->write(reinterpret_cast<char*>(&nj), sizeof(int));
    vfile->write(reinterpret_cast<char*>(&nk), sizeof(int));

    RDouble *x = this->GetX();
    RDouble *y = this->GetY();
    RDouble *z = this->GetZ();

    int nTotalNode = this->GetNTotalNode();

    int axisup = 1;
    GlobalDataBase::GetData("axisup", &axisup, PHINT, 1);
    if (axisup == 2)
    {
        TurnZAxis2YAxis();
    }

    vfile->write(reinterpret_cast<char*>(x), nTotalNode * sizeof(RDouble));
    vfile->write(reinterpret_cast<char*>(y), nTotalNode * sizeof(RDouble));
    vfile->write(reinterpret_cast<char*>(z), nTotalNode * sizeof(RDouble));

    int imin, imax, jmin, jmax, kmin, kmax, bctype;
    int s_lr, s_nd;

    CompositeStructBC *composite_bcregion = this->GetCompositeBCRegion();
    int nBCRegion = composite_bcregion->GetnBCRegion();

    vfile->write(reinterpret_cast<char*>(&nBCRegion), sizeof(int));

    //! Read the source region imin,imax,jmin,jmax,kmin,kmax,bctype.
    for (int ibcregion = 0; ibcregion < nBCRegion; ++ ibcregion)
    {
        StructBC *bcregion = composite_bcregion->GetBCRegion(ibcregion);
        bcregion->GetIJKRegion(imin, imax, jmin, jmax, kmin, kmax);
        bctype = bcregion->GetBCType();
        s_lr = bcregion->GetFaceLeftOrRightIndex();
        s_nd = bcregion->GetFaceDirection();

        vfile->write(reinterpret_cast<char*>(&imin), sizeof(int));
        vfile->write(reinterpret_cast<char*>(&imax), sizeof(int));
        vfile->write(reinterpret_cast<char*>(&jmin), sizeof(int));
        vfile->write(reinterpret_cast<char*>(&jmax), sizeof(int));
        vfile->write(reinterpret_cast<char*>(&kmin), sizeof(int));
        vfile->write(reinterpret_cast<char*>(&kmax), sizeof(int));
        vfile->write(reinterpret_cast<char*>(&bctype), sizeof(int));
        vfile->write(reinterpret_cast<char*>(&s_lr), sizeof(int));
        vfile->write(reinterpret_cast<char*>(&s_nd), sizeof(int));
    }

    InterfaceInfo *interfaceInfo = this->GetInterfaceInfo();

    int nIFace = 0;
    if (interfaceInfo)
    {
        nIFace = interfaceInfo->GetNIFace();
    }

    vfile->write(reinterpret_cast<char*>(&nIFace), sizeof(int));
    if (interfaceInfo)
    {
        int *interFace2ZoneID       = interfaceInfo->GetInterFace2ZoneID();
        int *interFace2InterFaceID  = interfaceInfo->GetInterFace2InterFaceID();
        int *interFace2BoundaryFace = interfaceInfo->GetInterFace2BoundaryFace();

        vfile->write(reinterpret_cast<char*>(interFace2ZoneID     ), nIFace * sizeof(int));
        vfile->write(reinterpret_cast<char*>(interFace2InterFaceID), nIFace * sizeof(int));
        vfile->write(reinterpret_cast<char*>(interFace2BoundaryFace), nIFace * sizeof(int));
    }
}

LIB_EXPORT void StructGrid::AllocateOversetGrid()
{
    oversetGridTopology->CreateCellContainer(nTotalCell);
    UnstructGrid *gridUnstr = this->GetUnstrGrid();
    if (!gridUnstr)
    {
        return;
    }
    gridUnstr->AllocateOversetGrid();
    int *iBlank = gridUnstr->GetBlankIndex();
    oversetGridTopology->SetCellContainer(iBlank);
}

LIB_EXPORT void StructGrid::CopyOversetFromUnstructGrid()
{
    UnstructGrid *gridUnstr = this->GetUnstrGrid();
    if (!gridUnstr) return;

    int *iBlank = gridUnstr->GetBlankIndex();
    CopyIBlank3D(iBlank);
}

LIB_EXPORT int *StructGrid::GetCellTypeContainer()
{
    return oversetGridTopology->GetCellTypeContainer();
}

LIB_EXPORT void StructGrid::GetBoundaryMinMax(RDouble *pmin, RDouble *pmax)
{
    using namespace PHSPACE;

    RDouble3D &xx = * this->GetStructX();
    RDouble3D &yy = * this->GetStructY();
    RDouble3D &zz = * this->GetStructZ();

    pmin[0] = LARGE;
    pmin[1] = LARGE;
    pmin[2] = LARGE;

    pmax[0] = - LARGE;
    pmax[1] = - LARGE;
    pmax[2] = - LARGE;

    int imin, imax, jmin, jmax, kmin, kmax, bctype;

    CompositeStructBC *composite_bcregion = this->GetCompositeBCRegion();
    int nBCRegion = composite_bcregion->GetnBCRegion();

    for (int ibcregion = 0; ibcregion < nBCRegion; ++ ibcregion)
    {
        StructBC *bcregion = composite_bcregion->GetBCRegion(ibcregion);
        bcregion->GetIJKRegion(imin, imax, jmin, jmax, kmin, kmax);
        bctype = bcregion->GetBCType();

        int nsurf = bcregion->GetFaceDirection() + 1;

        int il1, jl1, kl1;
        this->GetNsurfIndex( il1, jl1, kl1, nsurf);

        Range I(imin, imax + il1);
        Range J(jmin, jmax + jl1);
        Range K(kmin, kmax + kl1);

        int ist, ied, jst, jed, kst, ked;
        GetIJKRegion(I, J, K, ist, ied, jst, jed, kst, ked);

        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    RDouble xm, ym, zm;
                    xm = xx(i, j, k);
                    ym = yy(i, j, k);
                    zm = zz(i, j, k);
                    pmin[0] = MIN(static_cast<RDouble>(pmin[0]), xm);
                    pmin[1] = MIN(static_cast<RDouble>(pmin[1]), ym);
                    pmin[2] = MIN(static_cast<RDouble>(pmin[2]), zm);

                    pmax[0] = MAX(static_cast<RDouble>(pmax[0]), xm);
                    pmax[1] = MAX(static_cast<RDouble>(pmax[1]), ym);
                    pmax[2] = MAX(static_cast<RDouble>(pmax[2]), zm);
                }
            }
        }
    }
}

LIB_EXPORT void StructGrid::GetMinMaxDS(RDouble &dismin, RDouble &dismax)
{
    using namespace PHSPACE;
    //! Check the minimum grid distance.
    dismin =   LARGE;
    dismax = - LARGE;

    RDouble3D &xs = * this->GetStructX();
    RDouble3D &ys = * this->GetStructY();
    RDouble3D &zs = * this->GetStructZ();

    int nDim = GetDim();

    int ist, ied, jst, jed, kst, ked;
    this->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    RDouble dx[12], dy[12], dz[12];
    if (nDim == THREE_D)
    {
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    dx[0]  = xs(i+1, j  , k ) - xs(i  , j  , k );
                    dy[0]  = ys(i+1, j  , k ) - ys(i  , j  , k );
                    dz[0]  = zs(i+1, j  , k ) - zs(i  , j  , k );

                    dx[1]  = xs(i+1, j+1, k ) - xs(i+1, j  , k );
                    dy[1]  = ys(i+1, j+1, k ) - ys(i+1, j  , k );
                    dz[1]  = zs(i+1, j+1, k ) - zs(i+1, j  , k );

                    dx[2]  = xs(i  , j+1, k ) - xs(i+1, j+1, k );
                    dy[2]  = ys(i  , j+1, k ) - ys(i+1, j+1, k );
                    dz[2]  = zs(i  , j+1, k ) - zs(i+1, j+1, k );

                    dx[3]  = xs(i  , j  , k ) - xs(i  , j+1, k );
                    dy[3]  = ys(i  , j  , k ) - ys(i  , j+1, k );
                    dz[3]  = zs(i  , j  , k ) - zs(i  , j+1, k );

                    dx[4]  = xs(i+1, j  , k+1) - xs(i  , j  , k+1);
                    dy[4]  = ys(i+1, j  , k+1) - ys(i  , j  , k+1);
                    dz[4]  = zs(i+1, j  , k+1) - zs(i  , j  , k+1);

                    dx[5]  = xs(i+1, j+1, k+1) - xs(i+1, j  , k+1);
                    dy[5]  = ys(i+1, j+1, k+1) - ys(i+1, j  , k+1);
                    dz[5]  = zs(i+1, j+1, k+1) - zs(i+1, j  , k+1);

                    dx[6]  = xs(i  , j+1, k+1) - xs(i+1, j+1, k+1);
                    dy[6]  = ys(i  , j+1, k+1) - ys(i+1, j+1, k+1);
                    dz[6]  = zs(i  , j+1, k+1) - zs(i+1, j+1, k+1);

                    dx[7]  = xs(i  , j  , k+1) - xs(i  , j+1, k+1);
                    dy[7]  = ys(i  , j  , k+1) - ys(i  , j+1, k+1);
                    dz[7]  = zs(i  , j  , k+1) - zs(i  , j+1, k+1);

                    dx[8]  = xs(i  , j  , k+1) - xs(i  , j  , k );
                    dy[8]  = ys(i  , j  , k+1) - ys(i  , j  , k );
                    dz[8]  = zs(i  , j  , k+1) - zs(i  , j  , k );

                    dx[9]  = xs(i+1, j  , k+1) - xs(i+1, j  , k );
                    dy[9]  = ys(i+1, j  , k+1) - ys(i+1, j  , k );
                    dz[9]  = zs(i+1, j  , k+1) - zs(i+1, j  , k );

                    dx[10] = xs(i+1, j+1, k+1) - xs(i+1, j+1, k );
                    dy[10] = ys(i+1, j+1, k+1) - ys(i+1, j+1, k );
                    dz[10] = zs(i+1, j+1, k+1) - zs(i+1, j+1, k );

                    dx[11] = xs(i  , j+1, k+1) - xs(i  , j+1, k );
                    dy[11] = ys(i  , j+1, k+1) - ys(i  , j+1, k );
                    dz[11] = zs(i  , j+1, k+1) - zs(i  , j+1, k );

                    for (int m = 0; m < 12; ++ m)
                    {
                        RDouble ds = sqrt(dx[m] * dx[m] + dy[m] * dy[m] + dz[m] * dz[m]);
                        if (ds <= SAME_POINT_TOL) continue;
                        dismin = MIN(dismin, ds);
                        dismax = MAX(dismax, ds);
                    }
                }
            }
        }
    }
    else
    {
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    dx[0] = xs(i+1, j  , k) - xs(i  , j  , k);
                    dy[0] = ys(i+1, j  , k) - ys(i  , j  , k);

                    dx[1] = xs(i+1, j+1, k) - xs(i+1, j  , k);
                    dy[1] = ys(i+1, j+1, k) - ys(i+1, j  , k);

                    dx[2] = xs(i  , j+1, k) - xs(i+1, j+1, k);
                    dy[2] = ys(i  , j+1, k) - ys(i+1, j+1, k);

                    dx[3] = xs(i  , j  , k) - xs(i  , j+1, k);
                    dy[3] = ys(i  , j  , k) - ys(i  , j+1, k);
                    for (int m = 0; m < 4; ++ m)
                    {
                        RDouble ds = sqrt(dx[m] * dx[m] + dy[m] * dy[m]);
                        if (ds <= SAME_POINT_TOL) continue;
                        dismin = MIN(dismin, ds);
                        dismax = MAX(dismax, ds);
                    }
                }
            }
        }
    }
}

LIB_EXPORT void StructGrid::ToTecplot()
{
    fstream file;

    file.open("teststructgrid.plt", ios_base::out);

    int ni = this->GetNI();
    int nj = this->GetNJ();
    int nk = this->GetNK();

    file << "variables = \"x\" \"y\" \"z\" \n";

    file << "zone i = " << ni << " j = " << nj << " k = " << nk << " f = block \n";

    RDouble3D &xx = * GetStructX();
    RDouble3D &yy = * GetStructY();
    RDouble3D &zz = * GetStructZ();
    int count = 0;
    for (int k = 1; k <= nk; ++ k)
    {
        for (int j = 1; j <= nj; ++ j)
        {
            for (int i = 1; i <= ni; ++ i)
            {
                file << xx(i, j, k) << " ";
                ++ count;
                if (count % 5 == 0) file << "\n";
            }
        }
    }

    for (int k = 1; k <= nk; ++ k)
    {
        for (int j = 1; j <= nj; ++ j)
        {
            for (int i = 1; i <= ni; ++ i)
            {
                file << yy(i, j, k) << " ";
                ++ count;
                if (count % 5 == 0) file << "\n";
            }
        }
    }

    for (int k = 1; k <= nk; ++ k)
    {
        for (int j = 1; j <= nj; ++ j)
        {
            for (int i = 1; i <= ni; ++ i)
            {
                file << zz(i, j, k) << " ";
                ++ count;
                if (count % 5 == 0) file << "\n";
            }
        }
    }

    file.close();
    file.clear();

    PHMPI::FreeBlockData();
    cout << " " << endl;
    MPI_Abort(MPI_COMM_WORLD, -1);
}

LIB_EXPORT void StructGrid::GridSurfaceVelocity(RDouble *xt, RDouble *yt, RDouble *zt)
{
    int ni = this->GetNI();
    int nj = this->GetNJ();
    int nk = this->GetNK();

    int ndim = GetDim();

    Range I, J, K;
    GetRange(ni, nj, nk, 0, 0, I, J, K);

    Range M(1, ndim);

    RDouble3D xxt(xt, I, J, K, neverDeleteData, fortranArray);
    RDouble3D yyt(yt, I, J, K, neverDeleteData, fortranArray);
    RDouble3D zzt(zt, I, J, K, neverDeleteData, fortranArray);

    RDouble4D &faceNormalX = *(this->GetFaceNormalX());
    RDouble4D &faceNormalY = *(this->GetFaceNormalY());
    RDouble4D &faceNormalZ = *(this->GetFaceNormalZ());

    RDouble4D &faceVelocityX = *(this->GetFaceVelocityX());
    RDouble4D &faceVelocityY = *(this->GetFaceVelocityY());
    RDouble4D &faceVelocityZ = *(this->GetFaceVelocityZ());
    RDouble4D &faceNormalVelocity = *(this->GetFaceNormalVelocity());

    int ist, ied, jst, jed, kst, ked;
    ist = 1;
    ied = ni;
    jst = 1;
    jed = nj;
    kst = 1;
    ked = nk;

    if (nk == 1) ked = 1;

    int i1 = 1;
    int j1 = 1;
    int k1 = 1;

    if (nk == 1) k1 = 0;
    if (ndim == THREE_D)
    {
        int nd = 1;
        for (int k = kst; k < ked; ++ k)
        {
            for (int j = jst; j < jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    //! Face 1: (i,j,k), (i,j,k+1), (i,j+1,k), (i,j+1,k+1).
                    RDouble nx = faceNormalX(i, j, k, nd);
                    RDouble ny = faceNormalY(i, j, k, nd);
                    RDouble nz = faceNormalZ(i, j, k, nd);

                    RDouble a  = fourth * (xxt(i, j, k) + xxt(i, j, k + k1) + xxt(i, j + j1, k) + xxt(i, j + j1, k + k1));
                    RDouble b  = fourth * (yyt(i, j, k) + yyt(i, j, k + k1) + yyt(i, j + j1, k) + yyt(i, j + j1, k + k1));
                    RDouble c  = fourth * (zzt(i, j, k) + zzt(i, j, k + k1) + zzt(i, j + j1, k) + zzt(i, j + j1, k + k1));
                    RDouble d  =  (a * nx + b * ny + c * nz);

                    faceVelocityX(i, j, k, nd) = a;
                    faceVelocityY(i, j, k, nd) = b;
                    faceVelocityZ(i, j, k, nd) = c;
                    faceNormalVelocity(i, j, k, nd) = d;
                }
            }
        }

        //! Face 2: (i,j,k), (i,j,k+1), (i+1,j,k), (i+1,j,k+1).
        nd = 2;
        for (int k = kst; k < ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i < ied; ++ i)
                {
                    RDouble nx = faceNormalX(i, j, k, nd);
                    RDouble ny = faceNormalY(i, j, k, nd);
                    RDouble nz = faceNormalZ(i, j, k, nd);

                    RDouble a = fourth * (xxt(i, j, k) + xxt(i, j, k + k1) + xxt(i + i1, j, k) + xxt(i + i1, j, k + k1));
                    RDouble b = fourth * (yyt(i, j, k) + yyt(i, j, k + k1) + yyt(i + i1, j, k) + yyt(i + i1, j, k + k1));
                    RDouble c = fourth * (zzt(i, j, k) + zzt(i, j, k + k1) + zzt(i + i1, j, k) + zzt(i + i1, j, k + k1));
                    RDouble d  =  (a * nx + b * ny + c * nz);

                    faceVelocityX(i, j, k, nd) = a;
                    faceVelocityY(i, j, k, nd) = b;
                    faceVelocityZ(i, j, k, nd) = c;
                    faceNormalVelocity(i, j, k, nd) = d;
                }
            }
        }

        //! Face 3: (i,j,k), (i+1,j,k), (i,j+1,k), (i+1,j+1,k).
        nd = 3;
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j < jed; ++ j)
            {
                for (int i = ist; i < ied; ++ i)
                {
                    RDouble nx = faceNormalX(i, j, k, nd);
                    RDouble ny = faceNormalY(i, j, k, nd);
                    RDouble nz = faceNormalZ(i, j, k, nd);

                    RDouble a = fourth * (xxt(i, j, k) + xxt(i + i1, j, k) + xxt(i, j + j1, k) + xxt(i + i1, j + j1, k));
                    RDouble b = fourth * (yyt(i, j, k) + yyt(i + i1, j, k) + yyt(i, j + j1, k) + yyt(i + i1, j + j1, k));
                    RDouble c = fourth * (zzt(i, j, k) + zzt(i + i1, j, k) + zzt(i, j + j1, k) + zzt(i + i1, j + j1, k));
                    RDouble d  =  (a * nx + b * ny + c * nz);

                    faceVelocityX(i, j, k, nd) = a;
                    faceVelocityY(i, j, k, nd) = b;
                    faceVelocityZ(i, j, k, nd) = c;
                    faceNormalVelocity(i, j, k, nd) = d;
                }
            }
        }
    }
    else if (ndim == TWO_D)
    {
        int nd = 1;
        int k = 1;
        {
            for (int j = jst; j < jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    //! Face 1: (i,j,k), (i,j,k+1), (i,j+1,k), (i,j+1,k+1).
                    RDouble nx = faceNormalX(i, j, k, nd);
                    RDouble ny = faceNormalY(i, j, k, nd);
                    RDouble nz = faceNormalZ(i, j, k, nd);

                    RDouble a = fourth * (xxt(i, j, k) + xxt(i, j, k + k1) + xxt(i, j + j1, k) + xxt(i, j + j1, k + k1));
                    RDouble b = fourth * (yyt(i, j, k) + yyt(i, j, k + k1) + yyt(i, j + j1, k) + yyt(i, j + j1, k + k1));
                    RDouble c = fourth * (zzt(i, j, k) + zzt(i, j, k + k1) + zzt(i, j + j1, k) + zzt(i, j + j1, k + k1));
                    RDouble d  =  (a * nx + b * ny + c * nz);

                    faceVelocityX(i, j, k, nd) = a;
                    faceVelocityY(i, j, k, nd) = b;
                    faceVelocityZ(i, j, k, nd) = c;
                    faceNormalVelocity(i, j, k, nd) = d;
                }
            }
        }
        //! Face 2: (i,j,k), (i,j,k+1), (i+1,j,k), (i+1,j,k+1).
        nd = 2;
        k = 1;
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i < ied; ++ i)
                {
                    RDouble nx = faceNormalX(i, j, k, nd);
                    RDouble ny = faceNormalY(i, j, k, nd);
                    RDouble nz = faceNormalZ(i, j, k, nd);

                    RDouble a = fourth * (xxt(i, j, k) + xxt(i, j, k + k1) + xxt(i + i1, j, k) + xxt(i + i1, j, k + k1));
                    RDouble b = fourth * (yyt(i, j, k) + yyt(i, j, k + k1) + yyt(i + i1, j, k) + yyt(i + i1, j, k + k1));
                    RDouble c = fourth * (zzt(i, j, k) + zzt(i, j, k + k1) + zzt(i + i1, j, k) + zzt(i + i1, j, k + k1));
                    RDouble d  =  (a * nx + b * ny + c * nz);

                    faceVelocityX(i, j, k, nd) = a;
                    faceVelocityY(i, j, k, nd) = b;
                    faceVelocityZ(i, j, k, nd) = c;
                    faceNormalVelocity(i, j, k, nd) = d;
                }
            }
        }
    }
}

LIB_EXPORT void StructGrid::ComputeGridFaceVelocity()
{
    int nTotalNode = this->GetNTotalNode();
    RDouble *xt = new RDouble[nTotalNode];
    RDouble *yt = new RDouble[nTotalNode];
    RDouble *zt = new RDouble[nTotalNode];

    GridVerticeVelocity(xt, yt, zt);
    GridSurfaceVelocity(xt, yt, zt);

    delete [] xt;
    delete [] yt;
    delete [] zt;
}

LIB_EXPORT void StructGrid::EncodeOversetGrid(DataContainer *cdata)
{

}

LIB_EXPORT void StructGrid::Encode(DataContainer *cdata, int flag)
{
    if (flag == 0)
    {
        EncodeGrid(cdata);
    }
    else if (flag == 1)
    {
        EncodeOversetGrid(cdata);
    }
}

LIB_EXPORT void StructGrid::Decode(DataContainer *cdata, int flag)
{
    if (flag == 0)
    {
        DecodeGrid(cdata);
    }
    else if (flag == 1)
    {
        DecodeOversetGrid(cdata);
    }
}

LIB_EXPORT void StructGrid::ReviseNeighborZoneIndex(int zoneStart)
{
    InterfaceInfo *interfaceInfo = this->GetInterfaceInfo();

    if (!interfaceInfo) return;

    int nIFace = interfaceInfo->GetNIFace();

    int *interFace2ZoneID = interfaceInfo->GetInterFace2ZoneID();

    for (int iFace = 0; iFace < nIFace; ++ iFace)
    {
        interFace2ZoneID[iFace] += zoneStart;
    }

    InterpointInformation *interpointInformation = this->GetInterpointInfo();

    if (!interpointInformation) return;

    int nIPoint = interpointInformation->GetNumberOfInterpoints();

    int *interPoint2ZoneID = interpointInformation->GetInterPoint2ZoneID();

    for (int iPoint = 0; iPoint < nIPoint; ++ iPoint)
    {
        interPoint2ZoneID[iPoint] += zoneStart;
    }
}

LIB_EXPORT void StructGrid::DecodeOversetGrid(DataContainer *cdata)
{

}

LIB_EXPORT void StructGrid::UpdateCoordinate(DataContainer *cdata)
{
    cdata->moveto_begin();
    cout << "读入结构网格动网格数据文件......\n";

    int nTotalNode = this->GetNTotalNode();

    RDouble *x = this->GetX();
    RDouble *y = this->GetY();
    RDouble *z = this->GetZ();

    cdata->read(reinterpret_cast<char *>(x), nTotalNode * sizeof(RDouble));
    cdata->read(reinterpret_cast<char *>(y), nTotalNode * sizeof(RDouble));
    cdata->read(reinterpret_cast<char *>(z), nTotalNode * sizeof(RDouble));

    ComputeMinMaxBox();
}

void StructGrid::GetCellCenter(ActionKey *actkey)
{
    InterfaceInfo *interfaceInfo = this->GetInterfaceInfo();
    if (!interfaceInfo) return;

    StructGrid *fgrid = StructGridCast(this->GetFinestGrid());
    InterfaceInfo *fiinfo = fgrid->GetInterfaceInfo();

    int level = this->GetLevel();

    RDouble3D &vol = * this->GetCellVolume();
    RDouble3D &xcc = * this->GetCellCenterX();
    RDouble3D &ycc = * this->GetCellCenterY();
    RDouble3D &zcc = * this->GetCellCenterZ();

    RDouble3D &xx = * this->GetStructX();
    RDouble3D &yy = * this->GetStructY();
    RDouble3D &zz = * this->GetStructZ();

    int ineighbor = fiinfo->FindIthNeighbor(actkey->ipos);
    int *faceIndexForSend = fiinfo->GetFaceIndexForSend(ineighbor);
    int nIFaceOfNeighbor = fiinfo->GetNIFaceOfNeighbor(ineighbor);

    DataContainer *cdata = actkey->GetData();
    cdata->moveto_begin();

    int iface, is, js, ks;
    bool ifCommunicateNodeLocation = CommunicateNodeLocation();
    for (int ifaceLocal = 0; ifaceLocal < nIFaceOfNeighbor; ++ ifaceLocal)
    {
        iface = faceIndexForSend[ifaceLocal];
        fgrid->GetSourceIndexIJK(iface, 1 , is, js, ks);
        fgrid->RemapMultigridIJK(level, is, js, ks);

        cdata->write(&xcc(is, js, ks), sizeof(RDouble));
        cdata->write(&ycc(is, js, ks), sizeof(RDouble));
        cdata->write(&zcc(is, js, ks), sizeof(RDouble));
        cdata->write(&vol(is, js, ks), sizeof(RDouble));

        if (!ifCommunicateNodeLocation)
        {
            continue;
        }

        int IJKofNodes_1st_Layer[4][3];
        int IJKofNodes_2nd_Layer[4][3];

        GetNodeIJKfromSourceIndex(fgrid, iface, 1, IJKofNodes_1st_Layer);
        GetNodeIJKfromSourceIndex(fgrid, iface, 2, IJKofNodes_2nd_Layer);
        
        for (int nodeID = 0; nodeID < 4; ++ nodeID)
        {
            cdata->write(&xx(IJKofNodes_1st_Layer[nodeID][0], IJKofNodes_1st_Layer[nodeID][1], IJKofNodes_1st_Layer[nodeID][2]), sizeof(RDouble));
            cdata->write(&yy(IJKofNodes_1st_Layer[nodeID][0], IJKofNodes_1st_Layer[nodeID][1], IJKofNodes_1st_Layer[nodeID][2]), sizeof(RDouble));
            cdata->write(&zz(IJKofNodes_1st_Layer[nodeID][0], IJKofNodes_1st_Layer[nodeID][1], IJKofNodes_1st_Layer[nodeID][2]), sizeof(RDouble));

            cdata->write(&xx(IJKofNodes_2nd_Layer[nodeID][0], IJKofNodes_2nd_Layer[nodeID][1], IJKofNodes_2nd_Layer[nodeID][2]), sizeof(RDouble));
            cdata->write(&yy(IJKofNodes_2nd_Layer[nodeID][0], IJKofNodes_2nd_Layer[nodeID][1], IJKofNodes_2nd_Layer[nodeID][2]), sizeof(RDouble));
            cdata->write(&zz(IJKofNodes_2nd_Layer[nodeID][0], IJKofNodes_2nd_Layer[nodeID][1], IJKofNodes_2nd_Layer[nodeID][2]), sizeof(RDouble));
        }
    }
}

void StructGrid::TranslateCellCenter(ActionKey *actkey)
{
    InterfaceInfo *interfaceInfo = this->GetInterfaceInfo();
    if (!interfaceInfo) return;

    StructGrid *fgrid = StructGridCast(this->GetFinestGrid());
    int level = this->GetLevel();
    CompositeStructBC *composite_bcregion = fgrid->GetCompositeBCRegion();
    int periodicType = PHSPACE::GlobalDataBase::GetIntParaFromDB("periodicType");
    double rotationAngle = PHSPACE::GlobalDataBase::GetDoubleParaFromDB("rotationAngle");
    rotationAngle = rotationAngle * PI / 180;

    RDouble3D &vol = * this->GetCellVolume();
    RDouble3D &xcc = * this->GetCellCenterX();
    RDouble3D &ycc = * this->GetCellCenterY();
    RDouble3D &zcc = * this->GetCellCenterZ();

    RDouble3D &xx = * this->GetStructX();
    RDouble3D &yy = * this->GetStructY();
    RDouble3D &zz = * this->GetStructZ();

    InterfaceInfo *fiinfo = fgrid->GetInterfaceInfo();
    int ineighbor = fiinfo->FindIthNeighbor(actkey->ipos);
    int nIFaceOfNeighbor = fiinfo->GetNIFaceOfNeighbor(ineighbor);
    int *faceIndexForRecv = fiinfo->GetFaceIndexForRecv(ineighbor);

    int iface, it, jt, kt;

    DataContainer *cdata = actkey->GetData();
    cdata->moveto_begin();

    bool ifCommunicateNodeLocation = CommunicateNodeLocation();
    for (int ifaceLocal = 0; ifaceLocal < nIFaceOfNeighbor; ++ ifaceLocal)
    {
        iface = faceIndexForRecv[ifaceLocal];
        fgrid->GetTargetIndexIJK(iface, 1 , it, jt, kt);
        fgrid->RemapMultigridIJK(level, it, jt, kt);

        cdata->read(&xcc(it, jt, kt), sizeof(RDouble));
        cdata->read(&ycc(it, jt, kt), sizeof(RDouble));
        cdata->read(&zcc(it, jt, kt), sizeof(RDouble));
        cdata->read(&vol(it, jt, kt), sizeof(RDouble));

        int *ibcregions = composite_bcregion->GetIFaceInfo();
        int ibcregion = ibcregions[iface];
        StructBC *bcregion = composite_bcregion->GetBCRegion(ibcregion);
        string bcName = bcregion->GetBCName();
        if (periodicType == 2)
        {
            RDouble rotcell[3] = {0.0, 0.0, 0.0};
            if (bcName == "Periodic_up")
            {
                rotcell[1] = ycc(it, jt, kt) * cos(2 * PI - rotationAngle) - zcc(it, jt, kt) * sin(2 * PI - rotationAngle);
                rotcell[2] = ycc(it, jt, kt) * sin(2 * PI - rotationAngle) + zcc(it, jt, kt) * cos(2 * PI - rotationAngle);
            }
            else if (bcName == "Periodic_down")
            {
                rotcell[1] = ycc(it, jt, kt) * cos(rotationAngle) - zcc(it, jt, kt) * sin(rotationAngle);
                rotcell[2] = ycc(it, jt, kt) * sin(rotationAngle) + zcc(it, jt, kt) * cos(rotationAngle);
            }
            ycc(it, jt, kt) = rotcell[1];
            zcc(it, jt, kt) = rotcell[2];
        }

        if (!ifCommunicateNodeLocation)
        {
            continue;
        }
        RDouble Data_1st_Layer[4][3] = {0.0};
        RDouble Data_2nd_Layer[4][3] = {0.0};

        for (int nodeID = 0; nodeID < 4; ++ nodeID)
        {
            cdata->read(&Data_1st_Layer[nodeID][0], sizeof(RDouble));
            cdata->read(&Data_1st_Layer[nodeID][1], sizeof(RDouble));
            cdata->read(&Data_1st_Layer[nodeID][2], sizeof(RDouble));

            cdata->read(&Data_2nd_Layer[nodeID][0], sizeof(RDouble));
            cdata->read(&Data_2nd_Layer[nodeID][1], sizeof(RDouble));
            cdata->read(&Data_2nd_Layer[nodeID][2], sizeof(RDouble));
            if (periodicType == 2)
            {
                RDouble rotData_1st[3] = {0.0};
                RDouble rotData_2nd[3] = {0.0};

                if (bcName == "Periodic_up")
                {
                    rotData_1st[1] = Data_1st_Layer[nodeID][1] * cos(2 * PI - rotationAngle) - Data_1st_Layer[nodeID][2] * sin(2 * PI - rotationAngle);
                    rotData_1st[2] = Data_1st_Layer[nodeID][1] * sin(2 * PI - rotationAngle) + Data_1st_Layer[nodeID][2] * cos(2 * PI - rotationAngle);

                    rotData_2nd[1] = Data_2nd_Layer[nodeID][1] * cos(2 * PI - rotationAngle) - Data_2nd_Layer[nodeID][2] * sin(2 * PI - rotationAngle);
                    rotData_2nd[2] = Data_2nd_Layer[nodeID][1] * sin(2 * PI - rotationAngle) + Data_2nd_Layer[nodeID][2] * cos(2 * PI - rotationAngle);
                }
                else if (bcName == "Periodic_down")
                {
                    rotData_1st[1] = Data_1st_Layer[nodeID][1] * cos(rotationAngle) - Data_1st_Layer[nodeID][2] * sin(rotationAngle);
                    rotData_1st[2] = Data_1st_Layer[nodeID][1] * sin(rotationAngle) + Data_1st_Layer[nodeID][2] * cos(rotationAngle);

                    rotData_2nd[1] = Data_2nd_Layer[nodeID][1] * cos(rotationAngle) - Data_2nd_Layer[nodeID][2] * sin(rotationAngle);
                    rotData_2nd[2] = Data_2nd_Layer[nodeID][1] * sin(rotationAngle) + Data_2nd_Layer[nodeID][2] * cos(rotationAngle);

                }
                Data_1st_Layer[nodeID][1] = rotData_1st[1];
                Data_1st_Layer[nodeID][2] = rotData_1st[2];

                Data_2nd_Layer[nodeID][1] = rotData_2nd[1];
                Data_2nd_Layer[nodeID][2] = rotData_2nd[2];
            }
        }  

        int IJKofNodes_1st_Layer[4][3];
        int IJKofNodes_2nd_Layer[4][3];
        GetNodeIJKfromTargetIndex(fgrid, iface, 1, IJKofNodes_1st_Layer);
        GetNodeIJKfromTargetIndex(fgrid, iface, 2, IJKofNodes_2nd_Layer);

        RDouble XYZofNodes_1st_Layer[4][3];
        for (int nodeID = 0; nodeID < 4; ++ nodeID)
        {
            XYZofNodes_1st_Layer[nodeID][0] = xx(IJKofNodes_1st_Layer[nodeID][0], IJKofNodes_1st_Layer[nodeID][1], IJKofNodes_1st_Layer[nodeID][2]);
            XYZofNodes_1st_Layer[nodeID][1] = yy(IJKofNodes_1st_Layer[nodeID][0], IJKofNodes_1st_Layer[nodeID][1], IJKofNodes_1st_Layer[nodeID][2]);
            XYZofNodes_1st_Layer[nodeID][2] = zz(IJKofNodes_1st_Layer[nodeID][0], IJKofNodes_1st_Layer[nodeID][1], IJKofNodes_1st_Layer[nodeID][2]);
        }
        
        int RelationBetweenTwoPairsofNodes[4];
        GetRelationBetweenTwoPairsofNodes(RelationBetweenTwoPairsofNodes, XYZofNodes_1st_Layer, Data_1st_Layer);
        
        for (int nodeID = 0; nodeID < 4; ++ nodeID)
        {
            int abcdID = RelationBetweenTwoPairsofNodes[nodeID];
            xx(IJKofNodes_2nd_Layer[nodeID][0], IJKofNodes_2nd_Layer[nodeID][1], IJKofNodes_2nd_Layer[nodeID][2]) = Data_2nd_Layer[abcdID][0];
            yy(IJKofNodes_2nd_Layer[nodeID][0], IJKofNodes_2nd_Layer[nodeID][1], IJKofNodes_2nd_Layer[nodeID][2]) = Data_2nd_Layer[abcdID][1];
            zz(IJKofNodes_2nd_Layer[nodeID][0], IJKofNodes_2nd_Layer[nodeID][1], IJKofNodes_2nd_Layer[nodeID][2]) = Data_2nd_Layer[abcdID][2];
        }
    }
}

void StructGrid::GetCellIBlank(ActionKey *actkey)
{
    GetUnstrGrid()->GetCellIBlank(actkey);
    return;
}

void StructGrid::TranslateCellIBlank(ActionKey *actkey)
{
    GetUnstrGrid()->TranslateCellIBlank(actkey);
    return;

    //InterfaceInfo * interfaceInfo = this->GetInterfaceInfo();
    //if (!interfaceInfo) return;

    //int ineighbor = interfaceInfo->FindIthNeighbor(actkey->ipos);

    //int nIFaceOfNeighbor = interfaceInfo->GetNIFaceOfNeighbor(ineighbor) / 2;
    //int *faceIndexForRecv = interfaceInfo->GetFaceIndexForRecv(ineighbor);

    //int *iBlank = this->GetUnstrGrid()->GetBlankIndex();
    //Int3D &cellMapStrToUns = *reinterpret_cast<Int3D *> (GetDataPtr("cellMapStrToUns"));

    //DataContainer *cdata = actkey->GetData();
    //cdata->moveto_begin();

    //int iFace, it, jt, kt;
    //for (int ifaceLocal = 0; ifaceLocal < nIFaceOfNeighbor; ++ ifaceLocal)
    //{
    //    iFace = faceIndexForRecv[ifaceLocal];
    //    this->GetTargetIndexIJK(iFace, 1 , it, jt, kt);
    //    int iCell = cellMapStrToUns(it, jt, kt);

    //    cdata->read(&iBlank[iCell],sizeof(int));
    //}
}

bool StructGrid::CommunicateNodeLocation()
{
    bool Rwithghost = true;

    int tscheme = GlobalDataBase::GetIntParaFromDB("tscheme");
    if (tscheme != LU_SGS)
    {
        Rwithghost = false;
        return Rwithghost;
    }

    int isFVMOrFDM = GlobalDataBase::GetIntParaFromDB("isFVMOrFDM");
    if (isFVMOrFDM != FV_METHOD)
    {
        Rwithghost = false;
        return Rwithghost;
    }

    int systemGridType = PHSPACE::GlobalDataBase::GetIntParaFromDB("sys_gridtype");
    if (systemGridType != STRUCTGRID)
    {
        Rwithghost = false;
        return Rwithghost;
    }

    int isOverset = PHSPACE::GlobalDataBase::GetIntParaFromDB("codeOfOversetGrid");
    if (isOverset)
    {
        Rwithghost = false;
        return Rwithghost;
    }

    int nMGLevel = PHSPACE::GlobalDataBase::GetIntParaFromDB("nMGLevel");
    if (nMGLevel != 1)
    {
        Rwithghost = false;
        return Rwithghost;
    }

    int compressible = GlobalDataBase::GetIntParaFromDB("compressible");
    if (compressible != 1)
    {
        Rwithghost = false;
        return Rwithghost;
    }

    int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");
    if (isUnsteady)
    {
        Rwithghost = false;
        return Rwithghost;
    }

    int nChemical = GlobalDataBase::GetIntParaFromDB("nchem");
    if (nChemical != 0)
    {
        Rwithghost = false;
        return Rwithghost;
    }

    int ifLowSpeedPrecon = GlobalDataBase::GetIntParaFromDB("ifLowSpeedPrecon");
    if (ifLowSpeedPrecon != 0)
    {
        Rwithghost = false;
        return Rwithghost;
    }
    int wallFunctionType = GlobalDataBase::GetIntParaFromDB("wallFunctionType");
    if (wallFunctionType != WALLFUNCTION::NONE)
    {
        Rwithghost = false;
        return Rwithghost;
    }
    return Rwithghost;
}

void StructGrid::GetNodeIJKfromSourceIndex(StructGrid *fgrid, int iface, int layer, int IJKofNodes[4][3])
{
    int is, js, ks;
    fgrid->GetSourceIndexIJK(iface, layer, is, js, ks);

    int id, jd, kd;
    int i_lr, j_lr, k_lr;
    fgrid->GetSourceIndexIJK_fornodetransfer(iface, id, jd, kd, i_lr, j_lr, k_lr);

    int dijk1[3];
    int dijk2[3];
    int dijk3[3];

    if (i_lr != 0)
    {
        //! I_face.
        dijk1[0] = 0;
        dijk1[1] = 1;
        dijk1[2] = 0;

        dijk2[0] = 0;
        dijk2[1] = 0;
        dijk2[2] = 1;
    }
    else if (j_lr != 0)
    {
        //! J_face.
        dijk1[0] = 1;
        dijk1[1] = 0;
        dijk1[2] = 0;

        dijk2[0] = 0;
        dijk2[1] = 0;
        dijk2[2] = 1;
    }
    else
    {
        //! K_face.
        dijk1[0] = 1;
        dijk1[1] = 0;
        dijk1[2] = 0;

        dijk2[0] = 0;
        dijk2[1] = 1;
        dijk2[2] = 0;
    }

    if (GetDim() == 2)
    {
        dijk1[2] = 0;
        dijk2[2] = 0;
    }

    dijk3[0] = dijk1[0] + dijk2[0];
    dijk3[1] = dijk1[1] + dijk2[1];
    dijk3[2] = dijk1[2] + dijk2[2];

    IJKofNodes[0][0] = is + id;
    IJKofNodes[0][1] = js + jd;
    IJKofNodes[0][2] = ks + kd;

    IJKofNodes[1][0] = IJKofNodes[0][0] + dijk1[0];
    IJKofNodes[1][1] = IJKofNodes[0][1] + dijk1[1];
    IJKofNodes[1][2] = IJKofNodes[0][2] + dijk1[2];

    IJKofNodes[2][0] = IJKofNodes[0][0] + dijk2[0];
    IJKofNodes[2][1] = IJKofNodes[0][1] + dijk2[1];
    IJKofNodes[2][2] = IJKofNodes[0][2] + dijk2[2];

    IJKofNodes[3][0] = IJKofNodes[0][0] + dijk3[0];
    IJKofNodes[3][1] = IJKofNodes[0][1] + dijk3[1];
    IJKofNodes[3][2] = IJKofNodes[0][2] + dijk3[2];
}

void StructGrid::GetNodeIJKfromTargetIndex(StructGrid *fgrid, int iface, int layer, int IJKofNodes[4][3])
{
    int it, jt, kt;
    fgrid->GetTargetIndexIJK(iface, layer, it, jt, kt);

    int id, jd, kd;
    int i_lr, j_lr, k_lr;
    fgrid->GetTargetIndexIJK_fornodetransfer(iface, id, jd, kd, i_lr, j_lr, k_lr);

    int dijk1[3];
    int dijk2[3];
    int dijk3[3];

    if (i_lr != 0)
    {
        //! I_face.
        dijk1[0] = 0;
        dijk1[1] = 1;
        dijk1[2] = 0;

        dijk2[0] = 0;
        dijk2[1] = 0;
        dijk2[2] = 1;
    }
    else if (j_lr != 0)
    {
        //! J_face.
        dijk1[0] = 1;
        dijk1[1] = 0;
        dijk1[2] = 0;

        dijk2[0] = 0;
        dijk2[1] = 0;
        dijk2[2] = 1;
    }
    else
    {
        //! K_face.
        dijk1[0] = 1;
        dijk1[1] = 0;
        dijk1[2] = 0;

        dijk2[0] = 0;
        dijk2[1] = 1;
        dijk2[2] = 0;
    }

    if (GetDim() == 2)
    {
        dijk1[2] = 0;
        dijk2[2] = 0;
    }

    dijk3[0] = dijk1[0] + dijk2[0];
    dijk3[1] = dijk1[1] + dijk2[1];
    dijk3[2] = dijk1[2] + dijk2[2];

    IJKofNodes[0][0] = it + id;
    IJKofNodes[0][1] = jt + jd;
    IJKofNodes[0][2] = kt + kd;

    IJKofNodes[1][0] = IJKofNodes[0][0] + dijk1[0];
    IJKofNodes[1][1] = IJKofNodes[0][1] + dijk1[1];
    IJKofNodes[1][2] = IJKofNodes[0][2] + dijk1[2];

    IJKofNodes[2][0] = IJKofNodes[0][0] + dijk2[0];
    IJKofNodes[2][1] = IJKofNodes[0][1] + dijk2[1];
    IJKofNodes[2][2] = IJKofNodes[0][2] + dijk2[2];

    IJKofNodes[3][0] = IJKofNodes[0][0] + dijk3[0];
    IJKofNodes[3][1] = IJKofNodes[0][1] + dijk3[1];
    IJKofNodes[3][2] = IJKofNodes[0][2] + dijk3[2];
}

void StructGrid::GetRelationBetweenTwoPairsofNodes(int RelationBetweenTwoPairsofNodes[4], RDouble Nodes_ABCD[4][3], RDouble Nodes_abcd[4][3])
{
    for (int nodeID = 0; nodeID < 4; ++ nodeID)
    {
        RDouble length2_Aa = 0.0;
        RDouble length2_Ab = 0.0;
        RDouble length2_Ac = 0.0;
        RDouble length2_Ad = 0.0;

        for (int m = 0; m <= 2; ++ m)
        {
            length2_Aa += (Nodes_ABCD[nodeID][m] - Nodes_abcd[0][m]) * (Nodes_ABCD[nodeID][m] - Nodes_abcd[0][m]);
            length2_Ab += (Nodes_ABCD[nodeID][m] - Nodes_abcd[1][m]) * (Nodes_ABCD[nodeID][m] - Nodes_abcd[1][m]);
            length2_Ac += (Nodes_ABCD[nodeID][m] - Nodes_abcd[2][m]) * (Nodes_ABCD[nodeID][m] - Nodes_abcd[2][m]);
            length2_Ad += (Nodes_ABCD[nodeID][m] - Nodes_abcd[3][m]) * (Nodes_ABCD[nodeID][m] - Nodes_abcd[3][m]);
        }

        RDouble length2_Min = length2_Aa;
        RelationBetweenTwoPairsofNodes[nodeID] = 0;
        
        if (length2_Ab < length2_Min)
        {
            length2_Min = length2_Ab;
            RelationBetweenTwoPairsofNodes[nodeID] = 1;
        }
        if (length2_Ac < length2_Min)
        {
            length2_Min = length2_Ac;
            RelationBetweenTwoPairsofNodes[nodeID] = 2;
        }
        if (length2_Ad < length2_Min)
        {
            length2_Min = length2_Ad;
            RelationBetweenTwoPairsofNodes[nodeID] = 3;
        }
    }
}

LIB_EXPORT void StructGrid::GetSourceIndexIJK_fornodetransfer(int iface, int &id, int &jd, int &kd, int &i_lr, int &j_lr, int &k_lr)
{
    this->composite_bcregion->GetSourceIndexIJK_fornodetransfer(iface, id, jd, kd, i_lr, j_lr, k_lr);
}

LIB_EXPORT void StructGrid::GetTargetIndexIJK_fornodetransfer(int iface, int &id, int &jd, int &kd, int &i_lr, int &j_lr, int &k_lr)
{
    this->composite_bcregion->GetTargetIndexIJK_fornodetransfer(iface, id, jd, kd, i_lr, j_lr, k_lr);
}


void StructGrid::ComputeMetrics_BC()
{
    if (GetDim() == 2)
    {
        ComputeMetrics_BC_2D();
    }
    else
    {
        ComputeMetrics_BC_3D();
    }
}

void StructGrid::ComputeMetrics_BC_2D()
{
    RDouble4D &xfn  = * (this->GetFaceNormalX());
    RDouble4D &yfn  = * (this->GetFaceNormalY());
    RDouble4D &zfn  = * (this->GetFaceNormalZ());
    RDouble4D &area = * (this->GetFaceArea());
    RDouble4D &xfv  = * (this->GetFaceVectorX());
    RDouble4D &yfv  = * (this->GetFaceVectorY());
    RDouble4D &zfv  = * (this->GetFaceVectorZ());

    RDouble3D &xx = * this->GetStructX();
    RDouble3D &yy = * this->GetStructY();

    zfn = 0.0;
    zfv = 0.0;

    int ni = this->GetNI();
    int nj = this->GetNJ();

    RDouble r1[3], r2[3], r4[3];
    RDouble vkc[3], vet[3];
    int i, j;
    const int k = 1;

    i = 0;
    for (j = 1; j <= nj - 1; ++ j)
    {
        r1[0] = xx(i, j, k);
        r1[1] = yy(i, j, k);
        r1[2] = 0.0;

        r4[0] = xx(i, j+1, k);
        r4[1] = yy(i, j+1, k);
        r4[2] = 0.0;

        UnitFaceNormal2D(r4, r1, vkc);

        xfn (i, j, k, 1) = - vkc[0];
        yfn (i, j, k, 1) = - vkc[1];
        area(i, j, k, 1) =   vkc[2];
        xfv (i, j, k, 1) = - vkc[0] * vkc[2];
        yfv (i, j, k, 1) = - vkc[1] * vkc[2];
    }

    i = ni+1;
    for (j = 1; j <= nj - 1; ++ j)
    {
        r1[0] = xx(i, j, k);
        r1[1] = yy(i, j, k);
        r1[2] = 0.0;

        r4[0] = xx(i, j+1, k);
        r4[1] = yy(i, j+1, k);
        r4[2] = 0.0;

        UnitFaceNormal2D(r4, r1, vkc);

        xfn (i, j, k, 1) = - vkc[0];
        yfn (i, j, k, 1) = - vkc[1];
        area(i, j, k, 1) =   vkc[2];
        xfv (i, j, k, 1) = - vkc[0] * vkc[2];
        yfv (i, j, k, 1) = - vkc[1] * vkc[2];
    }

    i = 0;
    for (j = 1; j <= nj; ++ j)
    {
        r1[0] = xx(i, j, k);
        r1[1] = yy(i, j, k);
        r1[2] = 0.0;

        r2[0] = xx(i+1, j, k);
        r2[1] = yy(i+1, j, k);
        r2[2] = 0.0;

        UnitFaceNormal2D(r1, r2, vet);

        xfn (i, j, k, 2) = - vet[0];
        yfn (i, j, k, 2) = - vet[1];
        area(i, j, k, 2) =   vet[2];
        xfv (i, j, k, 2) = - vet[0] * vet[2];
        yfv (i, j, k, 2) = - vet[1] * vet[2];
    }

    i = ni;
    for (j = 1; j <= nj; ++ j)
    {
        r1[0] = xx(i, j, k);
        r1[1] = yy(i, j, k);
        r1[2] = 0.0;

        r2[0] = xx(i+1, j, k);
        r2[1] = yy(i+1, j, k);
        r2[2] = 0.0;

        UnitFaceNormal2D(r1, r2, vet);

        xfn (i, j, k, 2) = - vet[0];
        yfn (i, j, k, 2) = - vet[1];
        area(i, j, k, 2) =   vet[2];
        xfv (i, j, k, 2) = - vet[0] * vet[2];
        yfv (i, j, k, 2) = - vet[1] * vet[2];
    }

    j = 0;
    for (i = 1; i <= ni - 1; ++ i)
    {
        r1[0] = xx(i, j, k);
        r1[1] = yy(i, j, k);
        r1[2] = 0.0;

        r2[0] = xx(i+1, j, k);
        r2[1] = yy(i+1, j, k);
        r2[2] = 0.0;

        UnitFaceNormal2D(r1, r2, vet);

        xfn (i, j, k, 2) = - vet[0];
        yfn (i, j, k, 2) = - vet[1];
        area(i, j, k, 2) =   vet[2];
        xfv (i, j, k, 2) = - vet[0] * vet[2];
        yfv (i, j, k, 2) = - vet[1] * vet[2];
    }

    j = nj+1;
    for (i = 1; i <= ni - 1; ++ i)
    {
        r1[0] = xx(i, j, k);
        r1[1] = yy(i, j, k);
        r1[2] = 0.0;

        r2[0] = xx(i+1, j, k);
        r2[1] = yy(i+1, j, k);
        r2[2] = 0.0;

        UnitFaceNormal2D(r1, r2, vet);

        xfn (i, j, k, 2) = - vet[0];
        yfn (i, j, k, 2) = - vet[1];
        area(i, j, k, 2) =   vet[2];
        xfv (i, j, k, 2) = - vet[0] * vet[2];
        yfv (i, j, k, 2) = - vet[1] * vet[2];
    }

    j = 0;
    for (i = 1; i <= ni; ++ i)
    {
        r1[0] = xx(i, j, k);
        r1[1] = yy(i, j, k);
        r1[2] = 0.0;

        r4[0] = xx(i, j+1, k);
        r4[1] = yy(i, j+1, k);
        r4[2] = 0.0;

        UnitFaceNormal2D(r4, r1, vkc);

        xfn (i, j, k, 1) = - vkc[0];
        yfn (i, j, k, 1) = - vkc[1];
        area(i, j, k, 1) =   vkc[2];
        xfv (i, j, k, 1) = - vkc[0] * vkc[2];
        yfv (i, j, k, 1) = - vkc[1] * vkc[2];
    }

    j = nj;
    for (i = 1; i <= ni; ++ i)
    {
        r1[0] = xx(i, j, k);
        r1[1] = yy(i, j, k);
        r1[2] = 0.0;

        r4[0] = xx(i, j+1, k);
        r4[1] = yy(i, j+1, k);
        r4[2] = 0.0;

        UnitFaceNormal2D(r4, r1, vkc);

        xfn (i, j, k, 1) = - vkc[0];
        yfn (i, j, k, 1) = - vkc[1];
        area(i, j, k, 1) =   vkc[2];
        xfv (i, j, k, 1) = - vkc[0] * vkc[2];
        yfv (i, j, k, 1) = - vkc[1] * vkc[2];
    }
}

void StructGrid::ComputeMetrics_BC_3D()
{
    RDouble4D &xfn  = * (this->GetFaceNormalX());
    RDouble4D &yfn  = * (this->GetFaceNormalY());
    RDouble4D &zfn  = * (this->GetFaceNormalZ());
    RDouble4D &area = * (this->GetFaceArea());
    RDouble4D &xfv  = * (this->GetFaceVectorX());
    RDouble4D &yfv  = * (this->GetFaceVectorY());
    RDouble4D &zfv  = * (this->GetFaceVectorZ());

    RDouble3D &xx = * this->GetStructX();
    RDouble3D &yy = * this->GetStructY();
    RDouble3D &zz = * this->GetStructZ();

    int ni = this->GetNI();
    int nj = this->GetNJ();
    int nk = this->GetNK();

    RDouble r2[3], r3[3], r4[3], r5[3], r6[3], r7[3], r8[3];
    RDouble vkc[3], vet[3], vct[3], vkcn[4], vetn[4], vctn[4];

    int ibc[2];
    int jbc[2];
    int kbc[2];

    ibc[0] = 0;
    ibc[1] = ni+1;
    
    jbc[0] = 0;
    jbc[1] = nj+1;
    
    kbc[0] = 0;
    kbc[1] = nk+1;

    for (int bcindex = 0; bcindex <= 1; ++ bcindex)
    {
        int i = ibc[bcindex];
        for (int k = 1; k <= nk - 1; ++ k)
        {
            for (int j = 1; j <= nj - 1; ++ j)
            {
                r2[0] = xx(i, j    , k   );
                r2[1] = yy(i, j    , k   );
                r2[2] = zz(i, j    , k   );
                                     
                r5[0] = xx(i, j + 1, k   );
                r5[1] = yy(i, j + 1, k   );
                r5[2] = zz(i, j + 1, k   );
                                     
                r7[0] = xx(i, j    , k + 1);
                r7[1] = yy(i, j    , k + 1);
                r7[2] = zz(i, j    , k + 1);
                                     
                r8[0] = xx(i, j + 1, k + 1);
                r8[1] = yy(i, j + 1, k + 1);
                r8[2] = zz(i, j + 1, k + 1);
                
                FaceNormal3D(r2, r5, r8, r7, vkc);
                NormalizeVector(vkc, vkcn, 3);
                
                xfn (i, j, k, 1) = vkcn[0];
                yfn (i, j, k, 1) = vkcn[1];
                zfn (i, j, k, 1) = vkcn[2];
                area(i, j, k, 1) = vkcn[3];
                xfv (i, j, k, 1) = vkc[0];
                yfv (i, j, k, 1) = vkc[1];
                zfv (i, j, k, 1) = vkc[2];
            }
        } 
    }

    for (int bcindex = 0; bcindex <= 1; ++ bcindex)
    {
        int j = jbc[bcindex];        
        for (int k = 1; k <= nk - 1; ++ k)
        {
            for (int i = 1; i <= ni - 1; ++ i)
            {
                r3[0] = xx(i    , j, k   );
                r3[1] = yy(i    , j, k   );
                r3[2] = zz(i    , j, k   );
                                     
                r5[0] = xx(i + 1, j, k   );
                r5[1] = yy(i + 1, j, k   );
                r5[2] = zz(i + 1, j, k   );
                                     
                r6[0] = xx(i    , j, k + 1);
                r6[1] = yy(i    , j, k + 1);
                r6[2] = zz(i    , j, k + 1);
                                     
                r8[0] = xx(i + 1, j, k + 1);
                r8[1] = yy(i + 1, j, k + 1);
                r8[2] = zz(i + 1, j, k + 1);
                
                FaceNormal3D(r3, r6, r8, r5, vet);
                NormalizeVector(vet, vetn, 3);
                
                xfn (i, j, k, 2) = vetn[0];
                yfn (i, j, k, 2) = vetn[1];
                zfn (i, j, k, 2) = vetn[2];
                area(i, j, k, 2) = vetn[3];
                xfv (i, j, k, 2) = vet[0];
                yfv (i, j, k, 2) = vet[1];
                zfv (i, j, k, 2) = vet[2];
            }
        }
    }

    for (int bcindex = 0; bcindex <= 1; ++ bcindex)
    {
        int k = kbc[bcindex];
        for (int j = 1; j <= nj - 1; ++ j)
        {
            for (int i = 1; i <= ni - 1; ++ i)
            {
                r4[0] = xx(i    , j    , k);
                r4[1] = yy(i    , j    , k);
                r4[2] = zz(i    , j    , k);
                                         
                r6[0] = xx(i    , j + 1, k);
                r6[1] = yy(i    , j + 1, k);
                r6[2] = zz(i    , j + 1, k);
                                         
                r7[0] = xx(i + 1, j    , k);
                r7[1] = yy(i + 1, j    , k);
                r7[2] = zz(i + 1, j    , k);
                                         
                r8[0] = xx(i + 1, j + 1, k);
                r8[1] = yy(i + 1, j + 1, k);
                r8[2] = zz(i + 1, j + 1, k);
                
                FaceNormal3D(r4, r7, r8, r6, vct);
                NormalizeVector(vct, vctn, 3);
                
                xfn (i, j, k, 3) = vctn[0];
                yfn (i, j, k, 3) = vctn[1];
                zfn (i, j, k, 3) = vctn[2];
                area(i, j, k, 3) = vctn[3];
                xfv (i, j, k, 3) = vct[0];
                yfv (i, j, k, 3) = vct[1];
                zfv (i, j, k, 3) = vct[2];
            }
        }
    }

    ibc[0] = 0;
    ibc[1] = ni;
    
    jbc[0] = 0;
    jbc[1] = nj;
    
    kbc[0] = 0;
    kbc[1] = nk;

    for (int bcindex = 0; bcindex <= 1; ++ bcindex)
    {
        int i = ibc[bcindex];       
        
        for (int k = 1; k <= nk - 1; ++ k)
        {
            for (int j = 1; j <= nj; ++ j)
            {
                r3[0] = xx(i    , j, k   );
                r3[1] = yy(i    , j, k   );
                r3[2] = zz(i    , j, k   );
                                     
                r5[0] = xx(i + 1, j, k   );
                r5[1] = yy(i + 1, j, k   );
                r5[2] = zz(i + 1, j, k   );
                                     
                r6[0] = xx(i    , j, k + 1);
                r6[1] = yy(i    , j, k + 1);
                r6[2] = zz(i    , j, k + 1);
                                     
                r8[0] = xx(i + 1, j, k + 1);
                r8[1] = yy(i + 1, j, k + 1);
                r8[2] = zz(i + 1, j, k + 1);
                
                FaceNormal3D(r3, r6, r8, r5, vet);
                NormalizeVector(vet, vetn, 3);
                
                xfn (i, j, k, 2) = vetn[0];
                yfn (i, j, k, 2) = vetn[1];
                zfn (i, j, k, 2) = vetn[2];
                area(i, j, k, 2) = vetn[3];
                xfv (i, j, k, 2) = vet[0];
                yfv (i, j, k, 2) = vet[1];
                zfv (i, j, k, 2) = vet[2];
            }
        }
        
        for (int j = 1; j <= nj - 1; ++ j)
        {
            for (int k = 1; k <= nk; ++ k)
            {
                r4[0] = xx(i    , j    , k);
                r4[1] = yy(i    , j    , k);
                r4[2] = zz(i    , j    , k);
                                         
                r6[0] = xx(i    , j + 1, k);
                r6[1] = yy(i    , j + 1, k);
                r6[2] = zz(i    , j + 1, k);
                                         
                r7[0] = xx(i + 1, j    , k);
                r7[1] = yy(i + 1, j    , k);
                r7[2] = zz(i + 1, j    , k);
                                         
                r8[0] = xx(i + 1, j + 1, k);
                r8[1] = yy(i + 1, j + 1, k);
                r8[2] = zz(i + 1, j + 1, k);
                
                FaceNormal3D(r4, r7, r8, r6, vct);
                NormalizeVector(vct, vctn, 3);
                
                xfn (i, j, k, 3) = vctn[0];
                yfn (i, j, k, 3) = vctn[1];
                zfn (i, j, k, 3) = vctn[2];
                area(i, j, k, 3) = vctn[3];
                xfv (i, j, k, 3) = vct[0];
                yfv (i, j, k, 3) = vct[1];
                zfv (i, j, k, 3) = vct[2];
            }
        }
    }

    for (int bcindex = 0; bcindex <= 1; ++ bcindex)
    {
        int j = jbc[bcindex];
        for (int k = 1; k <= nk - 1; ++ k)
        {
            for (int i = 1; i <= ni; ++ i)
            {
                r2[0] = xx(i, j    , k   );
                r2[1] = yy(i, j    , k   );
                r2[2] = zz(i, j    , k   );
                                     
                r5[0] = xx(i, j + 1, k   );
                r5[1] = yy(i, j + 1, k   );
                r5[2] = zz(i, j + 1, k   );
                                     
                r7[0] = xx(i, j    , k + 1);
                r7[1] = yy(i, j    , k + 1);
                r7[2] = zz(i, j    , k + 1);
                                     
                r8[0] = xx(i, j + 1, k + 1);
                r8[1] = yy(i, j + 1, k + 1);
                r8[2] = zz(i, j + 1, k + 1);
                
                FaceNormal3D(r2, r5, r8, r7, vkc);
                NormalizeVector(vkc, vkcn, 3);
                
                xfn (i, j, k, 1) = vkcn[0];
                yfn (i, j, k, 1) = vkcn[1];
                zfn (i, j, k, 1) = vkcn[2];
                area(i, j, k, 1) = vkcn[3];
                xfv (i, j, k, 1) = vkc[0];
                yfv (i, j, k, 1) = vkc[1];
                zfv (i, j, k, 1) = vkc[2];
            }
        }    
        
        for (int i = 1; i <= ni - 1; ++ i)
        {
            for (int k = 1; k <= nk; ++ k)
            {
                r4[0] = xx(i    , j    , k);
                r4[1] = yy(i    , j    , k);
                r4[2] = zz(i    , j    , k);
                                         
                r6[0] = xx(i    , j + 1, k);
                r6[1] = yy(i    , j + 1, k);
                r6[2] = zz(i    , j + 1, k);
                                         
                r7[0] = xx(i + 1, j    , k);
                r7[1] = yy(i + 1, j    , k);
                r7[2] = zz(i + 1, j    , k);
                                         
                r8[0] = xx(i + 1, j + 1, k);
                r8[1] = yy(i + 1, j + 1, k);
                r8[2] = zz(i + 1, j + 1, k);
                
                FaceNormal3D(r4, r7, r8, r6, vct);
                NormalizeVector(vct, vctn, 3);
                
                xfn (i, j, k, 3) = vctn[0];
                yfn (i, j, k, 3) = vctn[1];
                zfn (i, j, k, 3) = vctn[2];
                area(i, j, k, 3) = vctn[3];
                xfv (i, j, k, 3) = vct[0];
                yfv (i, j, k, 3) = vct[1];
                zfv (i, j, k, 3) = vct[2];
            }
        }
    }

    for (int bcindex = 0; bcindex <= 1; ++ bcindex)
    {
        int k = kbc[bcindex];
        for (int j = 1; j <= nj - 1; ++ j)
        {
            for (int i = 1; i <= ni; ++ i)
            {
                r2[0] = xx(i, j    , k   );
                r2[1] = yy(i, j    , k   );
                r2[2] = zz(i, j    , k   );
                                     
                r5[0] = xx(i, j + 1, k   );
                r5[1] = yy(i, j + 1, k   );
                r5[2] = zz(i, j + 1, k   );
                                     
                r7[0] = xx(i, j    , k + 1);
                r7[1] = yy(i, j    , k + 1);
                r7[2] = zz(i, j    , k + 1);
                                     
                r8[0] = xx(i, j + 1, k + 1);
                r8[1] = yy(i, j + 1, k + 1);
                r8[2] = zz(i, j + 1, k + 1);
                
                FaceNormal3D(r2, r5, r8, r7, vkc);
                NormalizeVector(vkc, vkcn, 3);
                
                xfn (i, j, k, 1) = vkcn[0];
                yfn (i, j, k, 1) = vkcn[1];
                zfn (i, j, k, 1) = vkcn[2];
                area(i, j, k, 1) = vkcn[3];
                xfv (i, j, k, 1) = vkc[0];
                yfv (i, j, k, 1) = vkc[1];
                zfv (i, j, k, 1) = vkc[2];
            }
        }
        
        for (int j = 1; j <= nj; ++ j)
        {
            for (int i = 1; i <= ni - 1; ++ i)
            {
                r3[0] = xx(i    , j, k   );
                r3[1] = yy(i    , j, k   );
                r3[2] = zz(i    , j, k   );
                                     
                r5[0] = xx(i + 1, j, k   );
                r5[1] = yy(i + 1, j, k   );
                r5[2] = zz(i + 1, j, k   );
                                     
                r6[0] = xx(i    , j, k + 1);
                r6[1] = yy(i    , j, k + 1);
                r6[2] = zz(i    , j, k + 1);
                                     
                r8[0] = xx(i + 1, j, k + 1);
                r8[1] = yy(i + 1, j, k + 1);
                r8[2] = zz(i + 1, j, k + 1);
                
                FaceNormal3D(r3, r6, r8, r5, vet);
                NormalizeVector(vet, vetn, 3);
                
                xfn (i, j, k, 2) = vetn[0];
                yfn (i, j, k, 2) = vetn[1];
                zfn (i, j, k, 2) = vetn[2];
                area(i, j, k, 2) = vetn[3];
                xfv (i, j, k, 2) = vet[0];
                yfv (i, j, k, 2) = vet[1];
                zfv (i, j, k, 2) = vet[2];
            }
        }
    }
}

streamsize StructGrid::TranslateCellCenterLength(ActionKey *actkey)
{
    InterfaceInfo *interfaceInfo = this->GetInterfaceInfo();
    if (!interfaceInfo) return 0;

    StructGrid *fgrid = StructGridCast(this->GetFinestGrid());

    InterfaceInfo *fiinfo = fgrid->GetInterfaceInfo();
    int ineighbor = fiinfo->FindIthNeighbor(actkey->ipos);
    int nIFaceOfNeighbor = fiinfo->GetNIFaceOfNeighbor(ineighbor);

    streamsize nlen = nIFaceOfNeighbor * sizeof(RDouble) * 4;
    return nlen;

// for (int ifaceLocal = 0; ifaceLocal < nIFaceOfNeighbor; ++ ifaceLocal)
// {
//     iface = faceIndexForRecv[ifaceLocal];
//     fgrid->GetTargetIndexIJK(iface,1,it,jt,kt);
//     fgrid->RemapMultigridIJK(level,it,jt,kt);
// 
//     cdata->read(&xcc(it,jt,kt),sizeof(RDouble));
//     cdata->read(&ycc(it,jt,kt),sizeof(RDouble));
//     cdata->read(&zcc(it,jt,kt),sizeof(RDouble));
//     RDouble tmp;
//     cdata->read(&tmp,sizeof(RDouble));
// }
}

LIB_EXPORT void StructGrid::ChangeBCType(int from_bctype, int to_bctype)
{
    CompositeStructBC *composite_bcregion = this->GetCompositeBCRegion();
    int nBCRegion = composite_bcregion->GetnBCRegion();

    for (int ibcregion = 0; ibcregion < nBCRegion; ++ ibcregion)
    {
        StructBC *bcregion = composite_bcregion->GetBCRegion(ibcregion);

        if (bcregion->GetBCType() == from_bctype)
        {
            bcregion->SetBCType(to_bctype);
        }
    }
}

int StructGrid::CompNIFace()
{
    int nIFace = 0;
    int imin, imax, jmin, jmax, kmin, kmax, bctype;

    CompositeStructBC *composite_bcregion = this->GetCompositeBCRegion();
    int nBCRegion = composite_bcregion->GetnBCRegion();

    for (int ibcregion = 0; ibcregion < nBCRegion; ++ ibcregion)
    {
        StructBC *bcregion = composite_bcregion->GetBCRegion(ibcregion);
        bcregion->GetIJKRegion(imin, imax, jmin, jmax, kmin, kmax);
        bctype = bcregion->GetBCType();

        if (IsInterface(bctype))
        {
            nIFace += (imax - imin + 1) * (jmax - jmin + 1) * (kmax - kmin + 1);
        }
    }

    return nIFace;
}

void StructGrid::DecodeWallDist(DataContainer *cdata)
{
    cdata->moveto_begin();

    RDouble3D &walldist = * this->GetWallDist();

    int ist, ied, jst, jed, kst, ked;
    this->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                cdata->read(&walldist(i, j, k), sizeof(RDouble));
            }
        }
    }
}

LIB_EXPORT void StructGrid::RemapMultigridIJK(int level, int &i, int &j, int &k)
{
    StructGrid *grid = this;
    while (grid->GetCoarseGrid() && grid->GetLevel() < level)
    {
        int istp = grid->GetMultigridStepI();
        int jstp = grid->GetMultigridStepJ();
        int kstp = grid->GetMultigridStepK();

        i = MIN((i + 1) / istp, i);
        j = MIN((j + 1) / jstp, j);
        k = MIN((k + 1) / kstp, k);

        grid = StructGridCast(grid->GetCoarseGrid());
    };
}

Int3D * StructGrid::GetWallCellLabel()
{
    if (wallCellLabel == NULL)
    {
        int ni = this->GetNI();
        int nj = this->GetNJ();
        int nk = this->GetNK();
        Range I, J, K;
        GetRange(ni, nj, nk, 0, -1, I, J, K);
        wallCellLabel = new Int3D(I, J, K, fortranArray);
        Int3D &wallLabel = *this->wallCellLabel;

        wallLabel(I, J, K) = 0;

        CompositeStructBC *compositeBCRegion = this->GetCompositeBCRegion();
        int nBCRegion = compositeBCRegion->GetnBCRegion();
        for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
        {
            StructBC *structBC = compositeBCRegion->GetBCRegion(iBCRegion);

            int is, js, ks;
            int iStart, iEnd, jStart, jEnd, kStart, kEnd;
            structBC->GetIJKRegion(iStart, iEnd, jStart, jEnd, kStart, kEnd);

            for (int k = kStart; k <= kEnd; ++ k)
            {
                for (int j = jStart; j <= jEnd; ++ j)
                {
                    for (int i = iStart; i <= iEnd; ++ i)
                    {
                        //! First cell index inside flowfield.
                        structBC->GetInsideCellIndex(i, j, k, is, js, ks, 1);
                        RestrictIndex(is, js, ks, ni - 1, nj - 1, nk - 1);

                        wallLabel(is, js, ks) = 1;
                    }
                }
            }
        }
    }
    return wallCellLabel;
}

void StructGrid::InitWallDist()
{
    int ni = this->GetNI();
    int nj = this->GetNJ();
    int nk = this->GetNK();

    Range I(1, ni-1);
    Range J(1, nj-1);
    Range K(1, nk-1);

    if (nk == 1) K.setRange(1,1);

    walldist = new RDouble3D(I, J, K, fortranArray);
}

void StructGrid::InitIBlank3D()
{
    if (iBlank3D) delete iBlank3D;
    int ni = this->GetNI();
    int nj = this->GetNJ();
    int nk = this->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, -1, 0, I, J, K);

    iBlank3D = new Int3D(I, J, K, fortranArray);

    int ist, ied, jst, jed, kst, ked;
    this->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, ONE_GHOST_LAYER);

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                (*iBlank3D)(i, j, k) = GENERIC_COLOR;
            }
        }
    }
}

void StructGrid::CopyIBlank3D(int *iBlank)
{
    int ni = GetNI();
    int nj = GetNJ();
    int nk = GetNK();
    Range I, J, K;
    GetRange(ni, nj, nk, -1, 0, I, J, K);

    Int3D &cellMapStrToUns = *reinterpret_cast<Int3D *> (GetDataPtr("cellMapStrToUns"));
    int ist, ied, jst, jed, kst, ked;
    this->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, ONE_GHOST_LAYER);

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                int iCell = cellMapStrToUns(i, j, k);
                (*iBlank3D)(i, j, k) = iBlank[iCell];

            }
        }
    }
}
void StructGrid::SetIBlank(int *iBlank)
{
    int ni = GetNI();
    int nj = GetNJ();
    int nk = GetNK();
    Int3D &cellMapStrToUns = *reinterpret_cast<Int3D *> (GetDataPtr("cellMapStrToUns"));

    Range I, J, K;
    GetRange(ni, nj, nk, -1, 0, I, J, K);

    int ist, ied, jst, jed, kst, ked;
    this->GetCellIterationIndex(ist, ied, jst, jed, kst, ked, ONE_GHOST_LAYER);

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                int iCell = cellMapStrToUns(i, j, k);
                iBlank[iCell] = (*iBlank3D)(i, j, k);

            }
        }
    }
}

Int3D *StructGrid::GetIBlank3D()
{
    return this->iBlank3D;
}

void StructGrid::FillCoreIndex(vector< int > &core_index, int iCore, int jCore, int kCore)
{
    int ni = GetNI();
    int nj = GetNJ();

    int ic = ni + 1;
    int jc = nj + 1;
    core_index[0] = ic * jc * kCore + ic * jCore + iCore;
    core_index[1] = ic * jc * kCore + ic * jCore + (iCore + 1);
    core_index[2] = ic * jc * kCore + ic * (jCore + 1) + iCore;
    core_index[3] = ic * jc * (kCore + 1) + ic * jCore + iCore;
    core_index[4] = ic * jc * kCore + ic * (jCore + 1) + (iCore + 1);
    core_index[5] = ic * jc * (kCore + 1) + ic * jCore + (iCore + 1);
    core_index[6] = ic * jc * (kCore + 1) + ic * (jCore + 1) + iCore;
    core_index[7] = ic * jc * (kCore + 1) + ic * (jCore + 1) + (iCore + 1);

    return;
}

void StructGrid::ComputePhysicalCoordinate(vector< vector< RDouble > > &amer, vector< int > &core_index, vector< RDouble > &xx, vector< RDouble > &yy,  vector< RDouble > &zz)
{
    RDouble *xCoreContainer = GetXCoreContainer();
    RDouble *yCoreContainer = GetYCoreContainer();
    RDouble *zCoreContainer = GetZCoreContainer();
    for (int m = 0; m < 8; ++ m)
    {
        RDouble xg = xCoreContainer[core_index[m]];
        RDouble yg = yCoreContainer[core_index[m]];
        RDouble zg = zCoreContainer[core_index[m]];

        xx[m] = localFrame[0][0] + xg * localFrame[0][1] + yg * localFrame[0][2] + zg * localFrame[0][3];
        yy[m] = localFrame[1][0] + xg * localFrame[1][1] + yg * localFrame[1][2] + zg * localFrame[1][3];
        zz[m] = localFrame[2][0] + xg * localFrame[2][1] + yg * localFrame[2][2] + zg * localFrame[2][3];
    }

    return;
}

void StructGrid::ComputePhysicalCoordinate(vector< vector< RDouble > > &localFrame, RDouble xx, RDouble yy, RDouble zz, vector< RDouble > &globalCoordinate)
{
    for (int i = 0; i < 3; ++ i)
    {
        globalCoordinate[i] = localFrame[i][0] + xx * localFrame[i][1] + yy * localFrame[i][2] + zz * localFrame[i][3];
    }
    return;
}

void StructGrid::GaussElimination(vector< vector< RDouble > > &a, vector< RDouble > &b, int n, vector< RDouble > &x, int &flag)
{
    int is = 0;
    flag = 1;
    vector< int > js(n);

    for (int k = 0; k < n-1; ++ k)
    {
        RDouble d = 0.0;
        for (int i = k; i < n; ++ i)
        {
            for (int j = k; j < n; ++ j)
            {
                if (abs(a[i][j]) > d)
                {
                    d = abs(a[i][j]);
                    js[k] = j;
                    is = i;
                }
            }
        }

        if (d < SMALL)
        {
            flag = 0;
        }
        else
        {
            if (js[k] != k)
            {
                for (int i = 0; i < n; ++ i)
                {
                    RDouble t = a[i][k];
                    a[i][k] = a[i][js[k]];
                    a[i][js[k]] = t;
                }
            }

            if (is != k)
            {
                for (int j = k; j < n; ++ j)
                {
                    RDouble t = a[k][j];
                    a[k][j] = a[is][j];
                    a[is][j] = t;
                }

                RDouble t = b[k];
                b[k] = b[is];
                b[is] = t;
            }
        }

        if (flag == 0) return;

        RDouble sfn = 1.0 / a[k][k];
        for (int j = k+1; j < n; ++ j)
        {
            a[k][j] *= sfn;
        }

        b[k] *= sfn;

        for (int i = k+1; i < n; ++ i)
        {
            for (int j = k+1; j < n; ++ j)
            {
                a[i][j] -= a[i][k] * a[k][j];
            }
            b[i] -= a[i][k] * b[k];
        }
    }

    if (abs(a[n-1][n-1]) < SMALL)
    {
        flag = 0;
        return;
    }

    x[n-1] = b[n-1] / a[n-1][n-1];
    for (int i = n-2; i >= 0; -- i)
    {
        RDouble t = 0;
        for (int j = i+1; j < n; ++ j)
        {
            t += a[i][j] * x[j];
        }
        x[i] = b[i] - t;
    }

    js[n-1] = n - 1;
    for (int k = n-1; k >= 0; -- k)
    {
        if (js[k] != k)
        {
            RDouble t = x[k];
            x[k] = x[js[k]];
            x[js[k]] = t;
        }
    }

    return;
}

void StructGrid::ComputeJacobiMatrix(vector< RDouble > &a, vector< RDouble > &b, vector< RDouble > &c, vector< RDouble > &ksai, vector< vector< RDouble > > &amer)
{
    amer[0][0] = a[1] + a[4] * ksai[1] + a[5] * ksai[2] + a[7] * ksai[1] * ksai[2];
    amer[0][1] = a[2] + a[4] * ksai[0] + a[6] * ksai[2] + a[7] * ksai[0] * ksai[2];
    amer[0][2] = a[3] + a[5] * ksai[0] + a[6] * ksai[1] + a[7] * ksai[0] * ksai[1];

    amer[1][0] = b[1] + b[4] * ksai[1] + b[5] * ksai[2] + b[7] * ksai[1] * ksai[2];
    amer[1][1] = b[2] + b[4] * ksai[0] + b[6] * ksai[2] + b[7] * ksai[0] * ksai[2];
    amer[1][2] = b[3] + b[5] * ksai[0] + b[6] * ksai[1] + b[7] * ksai[0] * ksai[1];

    amer[2][0] = c[1] + c[4] * ksai[1] + c[5] * ksai[2] + c[7] * ksai[1] * ksai[2];
    amer[2][1] = c[2] + c[4] * ksai[0] + c[6] * ksai[2] + c[7] * ksai[0] * ksai[2];
    amer[2][2] = c[3] + c[5] * ksai[0] + c[6] * ksai[1] + c[7] * ksai[0] * ksai[1];

    return;
}

void StructGrid::ComputeTempContainer(vector< RDouble > &ksai, vector< RDouble > &temp)
{
    temp[0] = 1.0;
    temp[1] = ksai[0];
    temp[2] = ksai[1];
    temp[3] = ksai[2];
    temp[4] = ksai[0] * ksai[1];
    temp[5] = ksai[0] * ksai[2];
    temp[6] = ksai[1] * ksai[2];
    temp[7] = ksai[0] * ksai[1] * ksai[2];

    return;
}

void StructGrid::FillCoefficientContainer(vector< vector< RDouble > > &source, vector< RDouble > &target, int i)
{
    target[0] =  source[0][i];
    target[1] = -source[0][i] + source[1][i];
    target[2] = -source[0][i] + source[2][i];
    target[3] = -source[0][i] + source[3][i];
    target[4] =  source[0][i] - source[1][i] + source[4][i] - source[2][i];
    target[5] =  source[0][i] - source[1][i] - source[3][i] + source[5][i];
    target[6] =  source[0][i] - source[2][i] - source[3][i] + source[6][i];
    target[7] = -source[0][i] + source[1][i] - source[4][i] + source[2][i] + source[3][i] - source[5][i] + source[7][i] - source[6][i];

    return;
}

/*LIB_EXPORT void StructGrid::InitLinkPointLabel()
{
    iLinkPointLabel = 0;
    jLinkPointLabel = 0;
    kLinkPointLabel = 0;
}*/

LIB_EXPORT void StructGrid::CreateLinkPointStructure()
{
    oversetGridTopology->CreateLinkPointStructure(nTotalNode);
}

LIB_EXPORT void StructGrid::SetLinkPointLabel(int *iLabel, int iDirection)
{
    int *iLinkPointLabel = GetILinkPointLabel();
    int *jLinkPointLabel = GetJLinkPointLabel();
    int *kLinkPointLabel = GetKLinkPointLabel();
    if (iDirection == 0)
    {
        for (int ip = 0; ip < nTotalNode; ++ ip)
        {
            iLinkPointLabel[ip] = iLabel[ip];
        }
    }
    else if (iDirection == 1)
    {
        for (int ip = 0; ip < nTotalNode; ++ ip)
        {
            jLinkPointLabel[ip] = iLabel[ip];
        }
    }
    else
    {
        for (int ip = 0; ip < nTotalNode; ++ ip)
        {
            kLinkPointLabel[ip] = iLabel[ip];
        }
    }
}

LIB_EXPORT void StructGrid::SetHingedPointContainer(int *iLabel)
{
    int numberOfCores = GetNumberOfCores();
    int *hingedPointContainer = GetHingedPointContainer();
    for (int iCore = 0; iCore < numberOfCores; ++ iCore)
    {
        hingedPointContainer[iCore] = iLabel[iCore];
    }
}

LIB_EXPORT void StructGrid::SetCoreCoordinateContainer(RDouble *sLabel, int iDirection)
{
    int numberOfCores = GetNumberOfCores();
    RDouble *xCoreContainer = GetXCoreContainer();
    RDouble *yCoreContainer = GetYCoreContainer();
    RDouble *zCoreContainer = GetZCoreContainer();

    if (iDirection == 0)
    {
        for (int iCore = 0; iCore < numberOfCores; ++ iCore)
        {
            xCoreContainer[iCore] = sLabel[iCore];
        }
    }
    else if (iDirection == 1)
    {
        for (int iCore = 0; iCore < numberOfCores; ++ iCore)
        {
            yCoreContainer[iCore] = sLabel[iCore];
        }
    }
    else
    {
        for (int iCore = 0; iCore < numberOfCores; ++ iCore)
        {
            zCoreContainer[iCore] = sLabel[iCore];
        }
    }
}

LIB_EXPORT void StructGrid::CreateCoreParameterContainer()
{
    int ni = GetNI();
    int nj = GetNJ();
    int nk = GetNK();
    int numberOfCores = (ni+1) * (nj+1) * (nk+1);
    oversetGridTopology->SetNumberOfCores(numberOfCores);
    oversetGridTopology->CreateCoreParameterContainer(numberOfCores);
}

LIB_EXPORT void StructGrid::GatherCellsByColor(int &pCell, int color)
{
    int *iBlank = GetCellTypeContainer();
    for (int i = 0; i < nTotalCell; ++ i)
    {
        if (iBlank[i] == color)
        {
            pCell += 1;
        }
    }

    return;
}

LIB_EXPORT void StructGrid::ResetCellType()
{
    int *iBlank = GetCellTypeContainer();
    for (int pCounter = 0; pCounter < nTotalCell; ++ pCounter)
    {
        if (iBlank[pCounter] == SECOND_COLOR)
        {
            iBlank[pCounter] = OVERSET_COLOR;
        }
    }

    return;
}

LIB_EXPORT void StructGrid::ProbeExternalOverlapCells(int &counter)
{
    using namespace PHENGLEI;
    int ni = GetNI();
    int nj = GetNJ();
    int nk = GetNK();
    int *iBlank = GetCellTypeContainer();
    int nBcRegions = composite_bcregion->GetnBCRegion();
    for (int i = 0; i < nBcRegions; ++ i)
    {
        StructBC *bcRegion = composite_bcregion->GetBCRegion(i);

        int bcType = bcRegion->GetBCType();

        if (bcType != EXTERNAL_BC) continue;

        bcRegion->ProbeExternalOverlapCells(iBlank, ni - 1, nj - 1, nk - 1, counter);
    }

    return;
}

LIB_EXPORT void StructGrid::ProbeHoleCell(HoleGrid *holeGrid)
{
    RDouble3D *xcc = GetCellCenterX();
    RDouble3D *ycc = GetCellCenterY();
    RDouble3D *zcc = GetCellCenterZ();
    int *iBlank = GetCellTypeContainer();

    RDouble *holeBox = holeGrid->GetBox();
    for (int i = 0; i < 3; ++ i)
    {
        if (holeBox[2*i] > pmax[i] || holeBox[2*i+1] < pmin[i]) return;
    }

    int ist, ied, jst, jed, kst, ked;
    this->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    RDouble y_far = 2.0 * pmax[1] - pmin[1];

    int ic = 0;
    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble xx = (* xcc)(i, j, k);
                RDouble yy = (* ycc)(i, j, k);
                RDouble zz = (* zcc)(i, j, k);

                iBlank[ic] = holeGrid->GetCellColor(xx, yy, zz, y_far);

                ic += 1;
            }
        }
    }
}

LIB_EXPORT void StructGrid::ProbeBorderCell(int iType, int jType, vector< int > *friendPoint, vector< int > *friendZone)
{
    using namespace PHMPI;
    int ni = GetNI();
    int nj = GetNJ();
    int nk = GetNK();

    int *iLinkPointLabel = GetILinkPointLabel();
    int *jLinkPointLabel = GetJLinkPointLabel();
    int *kLinkPointLabel = GetKLinkPointLabel();
    int *iBlank = GetCellTypeContainer();

    int ist, ied, jst, jed, kst, ked;
    this->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    int ip, jp, kp;
    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                int cellLabel = (ni-1) * (nj-1) * (k-1) + (ni-1) * (j-1) + (i-1);
                int pointLabel = ni * nj * k + ni * j + i;
                if (iBlank[cellLabel] == iType)
                {
                    ip = i-1;
                    jp = j;
                    kp = k;

                    if (i == 1)
                    {
                        if (iLinkPointLabel[pointLabel] >= 0)
                        {
                            int globalPointLabel = iLinkPointLabel[pointLabel];
                            int zoneIndex = GetZoneIndexAccordingToGlobalPointLabel(globalPointLabel);
                            StructGrid *grid = PHSPACE::GetStructGrid(zoneIndex);
                                
                            if (grid)
                            {
                                grid->ColorOverlapCell(globalPointLabel, jType);
                            }
                            else
                            {
                                int processorIndex = GetZoneProcessorID(zoneIndex);
                                friendPoint[processorIndex].push_back(globalPointLabel);
                                friendZone[processorIndex].push_back(zoneIndex);
                            }
                        }
                    }
                    else
                    {
                        int twinCell = (ni-1) * (nj-1) * (kp-1) + (ni-1) * (jp-1) + (ip-1);
                        if (iBlank[twinCell] == GENERIC_COLOR)
                        {
                            iBlank[twinCell] = jType;
                        }
                    }

                    ip = i+1;
                    jp = j;
                    kp = k;

                    if (i == ied)
                    {
                        if (iLinkPointLabel[pointLabel] >= 0)
                        {
                            int globalPointLabel = iLinkPointLabel[pointLabel];
                            int zoneIndex = GetZoneIndexAccordingToGlobalPointLabel(globalPointLabel);
                                
                            StructGrid *grid = PHSPACE::GetStructGrid(zoneIndex);

                            if (grid)
                            {
                                grid->ColorOverlapCell(globalPointLabel, jType);
                            }
                            else
                            {
                                int processorIndex = GetZoneProcessorID(zoneIndex);
                                friendPoint[processorIndex].push_back(globalPointLabel);
                                friendZone[processorIndex].push_back(zoneIndex);
                            }
                        }
                    }
                    else
                    {
                        int twinCell = (ni-1) * (nj-1) * (kp-1) + (ni-1) * (jp-1) + (ip-1);
                        if (iBlank[twinCell] == GENERIC_COLOR)
                        {
                            iBlank[twinCell] = jType;
                        }
                    }

                    ip = i;
                    jp = j-1;
                    kp = k;

                    if (j == 1)
                    {
                        if (jLinkPointLabel[pointLabel] >= 0)
                        {
                            int globalPointLabel = jLinkPointLabel[pointLabel];
                            int zoneIndex = GetZoneIndexAccordingToGlobalPointLabel(globalPointLabel);
                                
                            StructGrid *grid = PHSPACE::GetStructGrid(zoneIndex);

                            if (grid)
                            {
                                grid->ColorOverlapCell(globalPointLabel, jType);
                            }
                            else
                            {
                                int processorIndex = GetZoneProcessorID(zoneIndex);
                                friendPoint[processorIndex].push_back(globalPointLabel);
                                friendZone[processorIndex].push_back(zoneIndex);
                            }
                        }
                    }
                    else
                    {
                        int twinCell = (ni-1) * (nj-1) * (kp-1) + (ni-1) * (jp-1) + (ip-1);
                        if (iBlank[twinCell] == GENERIC_COLOR)
                        {
                            iBlank[twinCell] = jType;
                        }
                    }

                    ip = i;
                    jp = j+1;
                    kp = k;

                    if (j == jed)
                    {
                        if (jLinkPointLabel[pointLabel] >= 0)
                        {
                            int globalPointLabel = jLinkPointLabel[pointLabel];
                            int zoneIndex = GetZoneIndexAccordingToGlobalPointLabel(globalPointLabel);
                                
                            StructGrid *grid = PHSPACE::GetStructGrid(zoneIndex);

                            if (grid)
                            {
                                grid->ColorOverlapCell(globalPointLabel, jType);
                            }
                            else
                            {
                                int processorIndex = GetZoneProcessorID(zoneIndex);
                                friendPoint[processorIndex].push_back(globalPointLabel);
                                friendZone[processorIndex].push_back(zoneIndex);
                            }
                        }
                    }
                    else
                    {
                        int twinCell = (ni-1) * (nj-1) * (kp-1) + (ni-1) * (jp-1) + (ip-1);
                        if (iBlank[twinCell] == GENERIC_COLOR)
                        {
                            iBlank[twinCell] = jType;
                        }
                    }

                    ip = i;
                    jp = j;
                    kp = k-1;

                    if (k == 1)
                    {
                        if (kLinkPointLabel[pointLabel] >= 0)
                        {
                            int globalPointLabel = kLinkPointLabel[pointLabel];
                            int zoneIndex = GetZoneIndexAccordingToGlobalPointLabel(globalPointLabel);
                                
                            StructGrid *grid = PHSPACE::GetStructGrid(zoneIndex);

                            if (grid)
                            {
                                grid->ColorOverlapCell(globalPointLabel, jType);
                            }
                            else
                            {
                                int processorIndex = GetZoneProcessorID(zoneIndex);
                                friendPoint[processorIndex].push_back(globalPointLabel);
                                friendZone[processorIndex].push_back(zoneIndex);
                            }
                        }
                    }
                    else
                    {
                        int twinCell = (ni-1) * (nj-1) * (kp-1) + (ni-1) * (jp-1) + (ip-1);
                        if (iBlank[twinCell] == GENERIC_COLOR)
                        {
                            iBlank[twinCell] = jType;
                        }
                    }

                    ip = i;
                    jp = j;
                    kp = k+1;

                    if (k == ked)
                    {
                        if (kLinkPointLabel[pointLabel] >= 0)
                        {
                            int globalPointLabel = kLinkPointLabel[pointLabel];
                            int zoneIndex = GetZoneIndexAccordingToGlobalPointLabel(globalPointLabel);
                                
                            StructGrid *grid = PHSPACE::GetStructGrid(zoneIndex);

                            if (grid)
                            {
                                grid->ColorOverlapCell(globalPointLabel, jType);
                            }
                            else
                            {
                                int processorIndex = GetZoneProcessorID(zoneIndex);
                                friendPoint[processorIndex].push_back(globalPointLabel);
                                friendZone[processorIndex].push_back(zoneIndex);
                            }
                        }
                    }
                    else
                    {
                        int twinCell = (ni-1) * (nj-1) * (kp-1) + (ni-1) * (jp-1) + (ip-1);
                        if (iBlank[twinCell] == GENERIC_COLOR)
                        {
                            iBlank[twinCell] = jType;
                        }
                    }
                }
            }
        }
    }
}

LIB_EXPORT void StructGrid::ColorOverlapCell(int gp, int type)
{
    int zoneStartPointLabel = GetZoneStartPointLabel();
    int *iBlank = GetCellTypeContainer();
    int ni = GetNI();
    int nj = GetNJ();
    int nk = GetNK();
    int mm = ni * nj;
    int ms = gp - zoneStartPointLabel;
    int kPoint = ms / mm;
    int other = ms % mm;
    int jPoint = other / ni;
    int iPoint = other % ni;
    int cellLabel = (ni-1) * (nj-1) * (kPoint-1) + (ni-1) * (jPoint-1) + (iPoint-1);
    if (iBlank[cellLabel] == GENERIC_COLOR)
    {
        iBlank[cellLabel] = type;
    }
}

LIB_EXPORT PHDouble2D & StructGrid::GetLocalFrame()
{
    int ni = GetNI();
    int nj = GetNJ();
    int nk = GetNK();
    localFrame.resize(3);
    for (int i = 0; i < 3; ++ i)
    {
        localFrame[i].resize(4);
    }

    int iHarf = ni / 2;
    int jHarf = nj / 2;
    int kHarf = nk / 2;

    int originalLabel = 0;
    int iDirectLabel = iHarf;
    int jDirectLabel = ni * jHarf;
    int kDirectLabel = ni * nj * kHarf;

    localFrame[0][0] = x[originalLabel];
    localFrame[0][1] = x[iDirectLabel] - x[originalLabel];
    localFrame[0][2] = x[jDirectLabel] - x[originalLabel];
    localFrame[0][3] = x[kDirectLabel] - x[originalLabel];

    localFrame[1][0] = y[originalLabel];
    localFrame[1][1] = y[iDirectLabel] - y[originalLabel];
    localFrame[1][2] = y[jDirectLabel] - y[originalLabel];
    localFrame[1][3] = y[kDirectLabel] - y[originalLabel];

    localFrame[2][0] = z[originalLabel];
    localFrame[2][1] = z[iDirectLabel] - z[originalLabel];
    localFrame[2][2] = z[jDirectLabel] - z[originalLabel];
    localFrame[2][3] = z[kDirectLabel] - z[originalLabel];

    return localFrame;
}

LIB_EXPORT void StructGrid::BuildGridLink(LinkStruct *link)
{
    using namespace PHSPACE;
    int zoneIndex = this->GetZoneID();

    cout << " zoneIndex = " << zoneIndex << "\n";
    int nIFace = this->CompNIFace();

    InterfaceInfo *interfaceInfo = 0;
    if (nIFace == 0)
    {
        this->SetInterfaceInfo(interfaceInfo);
        return;
    }

    interfaceInfo = new InterfaceInfo(nIFace);
    this->SetInterfaceInfo(interfaceInfo);
    int *interFace2BoundaryFace = interfaceInfo->GetInterFace2BoundaryFace();

    RDouble3D &xx = * this->GetStructX();
    RDouble3D &yy = * this->GetStructY();
    RDouble3D &zz = * this->GetStructZ();

    VVInt &facemap = link->GetFaceMap();

    facemap[zoneIndex].resize(nIFace);

    int ijk[4][3] = {0};
    int irot[3] = {0}, jrot[3] = {0}, krot[3] = {0};

    GetIJKROT(ijk, irot, jrot, krot);
    RDouble xlist[4] = {0}, ylist[4] = {0}, zlist[4] = {0};

    int il, jl, kl;
    int ir, jr, kr;

    int nlist = 4;
    if (this->GetDim() == 2) nlist = 2;

    vector<int> pindex(nlist);

    int imin, imax, jmin, jmax, kmin, kmax, bctype;

    int nsurf;

    LinkStruct::AdtTree &coor_tree = link->GetCoordinateTree();
    set < DataStruct_Sort< VInt > > &facelist = link->GetFaceList();
    VVInt &zoneid = link->GetZoneID();
    VVInt &faceid = link->GetFaceID();
    RDouble tol = link->GetTolerance();

    int iface = 0;
    uint_t fcount = zoneid.size();
    int pcount = coor_tree.GetNodeNum();

    int ibface = 0;

    int periodicType = PHSPACE::GlobalDataBase::GetIntParaFromDB("periodicType");
    RDouble translationLength[3] = {0.0};
    GlobalDataBase::GetData("translationLength", &translationLength, PHDOUBLE, 3);
    double rotationAngle = PHSPACE::GlobalDataBase::GetDoubleParaFromDB("rotationAngle");
    rotationAngle = rotationAngle * PI / 180;

    CompositeStructBC *composite_bcregion = this->GetCompositeBCRegion();
    int nBCRegion = composite_bcregion->GetnBCRegion();

    for (int ibcregion = 0; ibcregion < nBCRegion; ++ ibcregion)
    {
        StructBC *bcregion = composite_bcregion->GetBCRegion(ibcregion);
        bcregion->GetIJKRegion(imin, imax, jmin, jmax, kmin, kmax);
        bctype = bcregion->GetBCType();
        nsurf  = bcregion->GetFaceDirection() + 1;

        if (!IsInterface(bctype))
        {
            ibface += (imax - imin + 1) * (jmax - jmin + 1) * (kmax - kmin + 1);
            continue;
        }

        for (int k = kmin; k <= kmax; ++ k)
        {
            for (int j = jmin; j <= jmax; ++ j)
            {
                for (int i = imin; i <= imax; ++ i)
                {
                    bcregion->GetBoundaryFaceIndex(i, j, k, il, jl, kl);
                
                    ir = irot[nsurf - 1];
                    jr = jrot[nsurf - 1];
                    kr = krot[nsurf - 1];

                    GetFaceCoorList(il, jl, kl, ir, jr, kr, ijk, nlist, xlist, ylist, zlist, xx, yy, zz);
                    
                    if (bcregion->GetBCName() == "Periodic_up")
                    {
                        if (periodicType == 1)
                        {
                            for (int m = 0; m < nlist; ++ m)
                            {
                                xlist[m] += translationLength[0];
                                ylist[m] += translationLength[1];
                                zlist[m] += translationLength[2];
                            }
                        }
                        else if (periodicType == 2)
                        {
                            RDouble rotPoint[3]={0.0, 0.0, 0.0};
                            for (int m = 0; m < nlist; ++ m)
                            {
                                rotPoint[0] = xlist[m];
                                rotPoint[1] = ylist[m] * cos(rotationAngle) - zlist[m] * sin(rotationAngle);
                                rotPoint[2] = ylist[m] * sin(rotationAngle) + zlist[m] * cos(rotationAngle);
                                xlist[m] = rotPoint[0];
                                ylist[m] = rotPoint[1];
                                zlist[m] = rotPoint[2];
                            }
                        }
                    }

                    GetCoorIndexList(&coor_tree, tol, pcount, xlist, ylist, zlist, nlist, pindex);
                    Create_Link_Info(pindex, zoneIndex, iface, fcount, facelist, zoneid, faceid, link);

                    interFace2BoundaryFace[iface] = ibface;

                    ++ ibface;
                    ++ iface;
                }
            }
        }
    }

    uint_t nface = facelist.size();
    cout << "iface = " << iface << "\n";
}

LIB_EXPORT int StructGrid::GetNumberOfWallCell()
{
    int nsolid_surface = 0;

    CompositeStructBC *composite_bcregion = this->GetCompositeBCRegion();
    int nBCRegion = composite_bcregion->GetnBCRegion();

    for (int ibcregion = 0; ibcregion < nBCRegion; ++ ibcregion)
    {
        StructBC *bcregion = composite_bcregion->GetBCRegion(ibcregion);
        int bctype = bcregion->GetBCType();

        if (IsWall(bctype) || bctype == PHENGLEI::ABLATION_SURFACE)
        {
            int ist, ied, jst, jed, kst, ked;
            bcregion->GetIJKRegion(ist, ied, jst, jed, kst, ked);

            nsolid_surface += (ied - ist + 1) * (jed - jst + 1) * (ked - kst + 1);
        }
    }

    return nsolid_surface;
}

LIB_EXPORT void StructGrid::GetWallBCInfo(int &wallBCNum, int &maxCellNum)
{
    int nCell = 0;
    wallBCNum = 0;
    maxCellNum = 0;
    CompositeStructBC *composite_bcregion = this->GetCompositeBCRegion();
    int nBCRegion = composite_bcregion->GetnBCRegion();

    for (int ibcregion = 0; ibcregion < nBCRegion; ++ ibcregion)
    {
        StructBC *bcregion = composite_bcregion->GetBCRegion(ibcregion);
        int bctype = bcregion->GetBCType();

        if (IsWall(bctype) || bctype == PHENGLEI::ABLATION_SURFACE)
        {
            int ist, ied, jst, jed, kst, ked;
            bcregion->GetIJKRegion(ist, ied, jst, jed, kst, ked);

            nCell = (ied - ist + 1) * (jed - jst + 1) * (ked - kst + 1);
            wallBCNum ++;
            maxCellNum = MAX(maxCellNum, nCell);
        }
    }
}

LIB_EXPORT void StructGrid::ComputeWeight()
{

}

LIB_EXPORT void StructGrid::ComputeUnitCellCapsule(BackgroundTree *backgroundTree, int &iUnitCell)
{
    vector< int > core_index(8);
    vector< RDouble > xx(8), yy(8), zz(8);

    vector< vector< RDouble > > &amer = this->GetLocalFrame();

    int iZone = GetZoneID();
    int ni = GetNI();
    int nj = GetNJ();
    int nk = GetNK();
    for (int kCore = 0; kCore < nk; ++ kCore)
    {
        for (int jCore = 0; jCore < nj; ++ jCore)
        {
            for (int iCore = 0; iCore < ni; ++ iCore)
            {
                backgroundTree->SetUnitCellLocation(iUnitCell, iZone, iCore, jCore, kCore);

                FillCoreIndex(core_index, iCore, jCore, kCore);

                ComputePhysicalCoordinate(amer, core_index, xx, yy, zz);

                backgroundTree->SetUnitCellParameter(iUnitCell, xx, yy, zz);

                iUnitCell += 1;
            }
        }
    }
    
    return;
}

LIB_EXPORT void StructGrid::CollectOversetCells(OversetCellCollector *oversetCellCollector, int myid)
{
    int *iBlank = GetCellTypeContainer();
    RDouble3D &xx = * this->GetCellCenterX();
    RDouble3D &yy = * this->GetCellCenterY();
    RDouble3D &zz = * this->GetCellCenterZ();

    int ist, ied, jst, jed, kst, ked;
    this->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    vector< int > localCellIndex;
    vector< RDouble > xCellContainer, yCellContainer, zCellContainer;

    int counter = 0;
    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                if (iBlank[counter] == OVERSET_COLOR)
                {
                    localCellIndex.push_back(counter);

                    xCellContainer.push_back(xx(i, j, k));
                    yCellContainer.push_back(yy(i, j, k));
                    zCellContainer.push_back(zz(i, j, k));
                }
                counter += 1;
            }
        }
    }

    int zoneIndex = GetZoneID();
    oversetCellCollector->CollectOversetCells(zoneIndex, localCellIndex, xCellContainer, yCellContainer, zCellContainer);
}

LIB_EXPORT void StructGrid::ProbeNearestPoint(RDouble xx, RDouble yy, RDouble zz, vector< RDouble > &ss, vector< int > &gp)
{
    int zoneStartPointLabel = GetZoneStartPointLabel();
    RDouble *xCoreContainer = GetXCoreContainer();
    RDouble *yCoreContainer = GetYCoreContainer();
    RDouble *zCoreContainer = GetZCoreContainer();
    vector< RDouble > globalCoordinate(3);
    vector< vector< RDouble > > &localFrame = this->GetLocalFrame();
    int ni = GetNI();
    int nj = GetNJ();
    int nk = GetNK();
    int delta_k = max(1, nk-2);
    for (int kCore = 1; kCore < nk; kCore += delta_k)
    {
        for (int jCore = 1; jCore < nj; ++ jCore)
        {
            for (int iCore = 1; iCore < ni; ++ iCore)
            {
                int coreLabel = (ni+1) * (nj+1) * kCore + (ni+1) * jCore + iCore;
                int nodeLabel = zoneStartPointLabel + ni * nj * kCore + ni * jCore + iCore;

                gp.push_back(nodeLabel);

                RDouble xCore = xCoreContainer[coreLabel];
                RDouble yCore = yCoreContainer[coreLabel];
                RDouble zCore = zCoreContainer[coreLabel];

                ComputePhysicalCoordinate(localFrame, xCore, yCore, zCore, globalCoordinate);

                RDouble dX = globalCoordinate[0] - xx;
                RDouble dY = globalCoordinate[1] - yy;
                RDouble dZ = globalCoordinate[2] - zz;
                RDouble distance = sqrt(dX * dX + dY * dY + dZ * dZ);

                ss.push_back(distance);
            }
        }
    }

    int delta_j = max(1, nj-2);
    for (int jCore = 1; jCore < nj; jCore += delta_j)
    {
        for (int kCore = 1; kCore < nk; ++ kCore)
        {
            for (int iCore = 1; iCore < ni; ++ iCore)
            {
                int coreLabel = (ni+1) * (nj+1) * kCore + (ni+1) * jCore + iCore;
                int nodeLabel = zoneStartPointLabel + ni * nj * kCore + ni * jCore + iCore;

                gp.push_back(nodeLabel);

                RDouble xCore = xCoreContainer[coreLabel];
                RDouble yCore = yCoreContainer[coreLabel];
                RDouble zCore = zCoreContainer[coreLabel];

                ComputePhysicalCoordinate(localFrame, xCore, yCore, zCore, globalCoordinate);

                RDouble dX = globalCoordinate[0] - xx;
                RDouble dY = globalCoordinate[1] - yy;
                RDouble dZ = globalCoordinate[2] - zz;
                RDouble distance = sqrt(dX * dX + dY * dY + dZ * dZ);

                ss.push_back(distance);
            }
        }
    }

    int delta_i = max(1, ni-2);
    for (int iCore = 1; iCore < ni; iCore += delta_i)
    {
        for (int kCore = 1; kCore < nk; ++ kCore)
        {
            for (int jCore = 1; jCore < nj; ++ jCore)
            {
                int coreLabel = (ni+1) * (nj+1) * kCore + (ni+1) * jCore + iCore;
                int nodeLabel = zoneStartPointLabel + ni * nj * kCore + ni * jCore + iCore;

                gp.push_back(nodeLabel);

                RDouble xCore = xCoreContainer[coreLabel];
                RDouble yCore = yCoreContainer[coreLabel];
                RDouble zCore = zCoreContainer[coreLabel];

                ComputePhysicalCoordinate(localFrame, xCore, yCore, zCore, globalCoordinate);

                RDouble dX = globalCoordinate[0] - xx;
                RDouble dY = globalCoordinate[1] - yy;
                RDouble dZ = globalCoordinate[2] - zz;
                RDouble distance = sqrt(dX * dX + dY * dY + dZ * dZ);

                ss.push_back(distance);
            }
        }
    }

    return;
}

LIB_EXPORT void StructGrid::ComputeLocalCoordinateByNewtonMethod(int iCore, int jCore, int kCore, RDouble xx, RDouble yy, RDouble zz, vector< RDouble > &ss, vector< int > &node)
{
    RDouble *xCoreContainer = GetXCoreContainer();
    RDouble *yCoreContainer = GetYCoreContainer();
    RDouble *zCoreContainer = GetZCoreContainer();
    int *hingedPointContainer = GetHingedPointContainer();
    const RDouble LOWER = -0.001, UPPER = 1.001, EPS = 1.e-4;
    const int ITERATION_NUM = 20, NODE_NUM = 8;

    vector< vector< RDouble > > globalCoordinate(8, vector< RDouble >(3));
    vector< RDouble > a(8), b(8), c(8), temp(8);

    vector< int > core_index(NODE_NUM);
    FillCoreIndex(core_index, iCore, jCore, kCore);

    vector< vector< RDouble > > &localFrame = this->GetLocalFrame();
    for (int i = 0; i < NODE_NUM; ++ i)
    {
        int coreLabel = core_index[i];
        ComputePhysicalCoordinate(localFrame, xCoreContainer[coreLabel], yCoreContainer[coreLabel], zCoreContainer[coreLabel], globalCoordinate[i]);
    }

    FillCoefficientContainer(globalCoordinate, a, 0);
    FillCoefficientContainer(globalCoordinate, b, 1);
    FillCoefficientContainer(globalCoordinate, c, 2);

    vector< RDouble > ksai(3, 0.5), delta(3), rhs(3);
    vector< vector< RDouble > > amer(3, vector<RDouble>(3));

    int flag = 0;
    for (int iter = 0; iter < ITERATION_NUM; ++ iter)
    {
        ComputeTempContainer(ksai, temp);

        rhs[0] = -inner_product(a.begin(), a.end(), temp.begin(), -xx);
        rhs[1] = -inner_product(b.begin(), b.end(), temp.begin(), -yy);
        rhs[2] = -inner_product(c.begin(), c.end(), temp.begin(), -zz);

        ComputeJacobiMatrix(a, b, c, ksai, amer);

        GaussElimination(amer, rhs, 3, delta, flag);

        for (int i = 0; i < 3; ++ i)
        {
            ksai[i] += delta[i];
        }

        if (abs(delta[0]) + abs(delta[1]) + abs(delta[2]) < EPS) break;
    }

    if (flag == 1 && LOWER < ksai[0] && ksai[0] < UPPER && LOWER < ksai[1] && ksai[1] < UPPER && LOWER < ksai[2] && ksai[2] < UPPER)
    {
        ss = ksai;
        node.resize(NODE_NUM);
        for (int iCell = 0; iCell < NODE_NUM; ++ iCell)
        {
            node[iCell] = hingedPointContainer[core_index[iCell]];
        }
    }

    return;
}

LIB_EXPORT RDouble StructGrid::CalMinEdgeLength()
{
    if (this->minEdgeLength > 0.0)
    {
        return this->minEdgeLength;
    }

    PHMPI::FreeBlockData();
    cout << "    Error: have not calculate the edge length for structured grid!\n" << endl;
    MPI_Abort(MPI_COMM_WORLD, -1);

    return 0;
}

LIB_EXPORT RDouble StructGrid::CalMaxEdgeLength()
{
    if (this->maxEdgeLength > 0.0)
    {
        return this->maxEdgeLength;
    }

    maxEdgeLength = PHSPACE::SMALL;

    int ni = this->GetNI();
    int nj = this->GetNJ();
    int nk = this->GetNK();

    RDouble3D &xx = *this->GetStructX();
    RDouble3D &yy = *this->GetStructY();
    RDouble3D &zz = *this->GetStructZ();

    Range I, J, K;
    GetRange(ni, nj, nk, 0, -1, I, J, K);
    int ist, ied, jst, jed, kst, ked;
    this->GetNodeIterationIndex(ist, ied, jst, jed, kst, ked);

    int il1 = 1;
    int jl1 = 1;
    int kl1 = 1;

    RDouble dx, dy, dz, ds;
    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                int il = i + il1;
                if (il <= ied)
                {
                    dx = xx(i, j, k) - xx(il, j, k);
                    dy = yy(i, j, k) - yy(il, j, k);
                    dz = zz(i, j, k) - zz(il, j, k);
                    ds = DISTANCE(dx, dy, dz);
                    maxEdgeLength = MAX(maxEdgeLength, ds);
                }

                int jl = j + jl1;
                if (jl <= jed)
                {
                    dx = xx(i, j, k) - xx(i, jl, k);
                    dy = yy(i, j, k) - yy(i, jl, k);
                    dz = zz(i, j, k) - zz(i, jl, k);
                    ds = DISTANCE(dx, dy, dz);
                    maxEdgeLength = MAX(maxEdgeLength, ds);

                }

                if (ked == 1) continue;
                int kl = k + kl1;
                if (kl <= ked)
                {
                    dx = xx(i, j, k) - xx(i, j, kl);
                    dy = yy(i, j, k) - yy(i, j, kl);
                    dz = zz(i, j, k) - zz(i, j, kl);
                    ds = DISTANCE(dx, dy, dz);
                    maxEdgeLength = MAX(maxEdgeLength, ds);
                }
            }
        }
    }

    return maxEdgeLength;
}

void StructGrid::DecodeBCFace(DataContainer *cdata)
{
    VirtualFile *virtualFile = new VirtualFile(cdata);

    virtualFile->BeginReadWork();

    int zoneIndex;
    PHRead(virtualFile, zoneIndex);

    int numberOfBCRegions;
    PHRead(virtualFile, numberOfBCRegions);
    ASSERT(numberOfBCRegions == this->GetCompositeBCRegion()->GetnBCRegion());

    uint_long totalSize;
    PHRead(virtualFile, totalSize);

    char *bcNameChar = new char [totalSize];
    PHRead(virtualFile, bcNameChar, static_cast<int>(totalSize));
    
    string *bcNameList = new string [numberOfBCRegions];
    uint_long count = 0;
    for (int ibcregion = 0; ibcregion < numberOfBCRegions; ++ ibcregion)
    {
        if (count < totalSize)
        {
            while (bcNameChar[count] != '\0')
            {
                bcNameList[ibcregion].push_back(bcNameChar[count]);
                ++ count;
            }
            ++ count;
        }
    }

    //PHRead(virtualFile, numberOfBCRegions);
    //int *bcTypeList = new int [numberOfBCRegions];
    //PHRead(virtualFile, bcTypeList, numberOfBCRegions);

    virtualFile->EndReadWork();

    CompositeStructBC *composite_bcregion = this->GetCompositeBCRegion();
    int nBCRegion = composite_bcregion->GetnBCRegion();
    for (int ibcregion = 0; ibcregion < nBCRegion; ++ ibcregion)
    {
        StructBC *bcregion = composite_bcregion->GetBCRegion(ibcregion);
        bcregion->SetBCName(bcNameList[ibcregion]);

        int bctype = bcregion->GetBCType();
        if ((bctype > 7) && (bctype < 51))
        {
            bctype = 2;
            bcregion->SetBCType(bctype);
        }
    }

    //delete [] bcTypeList;
    delete [] bcNameList;
    delete [] bcNameChar;
    delete virtualFile;
}

void StructGrid::DecodeBcDirFace(DataContainer *cdata)
{
    VirtualFile *virtualFile = new VirtualFile(cdata);

    virtualFile->BeginReadWork();

    int zoneIndex;
    PHRead(virtualFile, zoneIndex);

    int numberOfBCRegions;
    PHRead(virtualFile, numberOfBCRegions);

    CompositeStructBC *composite_bcregion = this->GetCompositeBCRegion();
    ASSERT(numberOfBCRegions == composite_bcregion->GetnBCRegion());

    for (int ibcregion = 0; ibcregion < numberOfBCRegions; ++ ibcregion)
    {
        StructBC *structBC = composite_bcregion->GetBCRegion(ibcregion);
        int *s_dir3d = structBC->GetFaceMatchingTargetDirIndex();
        virtualFile->read(reinterpret_cast<char*>(s_dir3d), 3 * sizeof(int));
    }

    virtualFile->EndReadWork();

    delete virtualFile;
}

void StructGrid::GhostCell3DExceptInterface(RDouble3D &TargetArray3D)
{
    int ni = this->GetNI();
    int nj = this->GetNJ();
    int nk = this->GetNK();

    CompositeStructBC *compositeBCRegion = this->GetCompositeBCRegion();
    int nBCRegion = compositeBCRegion->GetnBCRegion();
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        StructBC *structBC = compositeBCRegion->GetBCRegion(iBCRegion);
        //int BCType = structBC->GetBCType();
        //if (BCType == PHENGLEI::INTERFACE || BCType < 0)
        //{
        //    continue;
        //}

        int iStart, iEnd, jStart, jEnd, kStart, kEnd;
        structBC->GetIJKRegion(iStart, iEnd, jStart, jEnd, kStart, kEnd);

        int *faceDirectionIndex3d = structBC->GetFaceDirectionIndex();
        for (int k = kStart; k <= kEnd; ++ k)
        {
            for (int j = jStart; j <= jEnd; ++ j)
            {
                for (int i = iStart; i <= iEnd; ++ i)
                {
                    //! First cell index inside flowfield.
                    int firstInsideLayer = 0;
                    int iStartInsideLayer1 = i - faceDirectionIndex3d[0] * firstInsideLayer;
                    int jStartInsideLayer1 = j - faceDirectionIndex3d[1] * firstInsideLayer;
                    int kStartInsideLayer1 = k - faceDirectionIndex3d[2] * firstInsideLayer;

                    RestrictIndex(iStartInsideLayer1, jStartInsideLayer1, kStartInsideLayer1, ni-1, nj-1, nk-1);
                    for (int nGhostLayer = 1; nGhostLayer <= 2; ++ nGhostLayer)
                    {
                        int iGhostLayer = i + faceDirectionIndex3d[0] * nGhostLayer;
                        int jGhostLayer = j + faceDirectionIndex3d[1] * nGhostLayer;
                        int kGhostLayer = k + faceDirectionIndex3d[2] * nGhostLayer;

                        TargetArray3D(iGhostLayer, jGhostLayer, kGhostLayer) = TargetArray3D(iStartInsideLayer1, jStartInsideLayer1, kStartInsideLayer1);
                    }
                }
            }
        }
    }
}

void StructGrid::GhostCell3DExceptInterface(RDouble4D &TargetArray4D, const int &nVariables)
{
    int ni = this->GetNI();
    int nj = this->GetNJ();
    int nk = this->GetNK();

    CompositeStructBC *compositeBCRegion = this->GetCompositeBCRegion();
    int nBCRegion = compositeBCRegion->GetnBCRegion();
    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        StructBC *structBC = compositeBCRegion->GetBCRegion(iBCRegion);
        //int BCType = structBC->GetBCType();
        //if (BCType == PHENGLEI::INTERFACE || BCType < 0)
        //{
        //    continue;
        //}

        int iStart, iEnd, jStart, jEnd, kStart, kEnd;
        structBC->GetIJKRegion(iStart, iEnd, jStart, jEnd, kStart, kEnd);

        int *faceDirectionIndex3d = structBC->GetFaceDirectionIndex();
        for (int k = kStart; k <= kEnd; ++ k)
        {
            for (int j = jStart; j <= jEnd; ++ j)
            {
                for (int i = iStart; i <= iEnd; ++ i)
                {
                    //! First cell index inside flowfield.
                    int firstInsideLayer = 0;
                    int iStartInsideLayer1 = i - faceDirectionIndex3d[0] * firstInsideLayer;
                    int jStartInsideLayer1 = j - faceDirectionIndex3d[1] * firstInsideLayer;
                    int kStartInsideLayer1 = k - faceDirectionIndex3d[2] * firstInsideLayer;

                    RestrictIndex(iStartInsideLayer1, jStartInsideLayer1, kStartInsideLayer1, ni-1, nj-1, nk-1);

                    for (int nGhostLayer = 1; nGhostLayer <= 2; ++ nGhostLayer)
                    {
                        int iGhostLayer = i + faceDirectionIndex3d[0] * nGhostLayer;
                        int jGhostLayer = j + faceDirectionIndex3d[1] * nGhostLayer;
                        int kGhostLayer = k + faceDirectionIndex3d[2] * nGhostLayer;
                        for (int iVar = 0; iVar < nVariables; ++ iVar)
                        {
                            TargetArray4D(iGhostLayer, jGhostLayer, kGhostLayer, iVar) = 
                            TargetArray4D(iStartInsideLayer1, jStartInsideLayer1, kStartInsideLayer1, iVar);
                        }
                    }
                }
            }
        }
    }
}

void StructGrid::GetBcNamePair(set< pair<int, string> > &bcNamePairSet)
{
    CompositeStructBC *composite_bcregion = this->GetCompositeBCRegion();
    int  nBCRegion = composite_bcregion->GetnBCRegion();

    for (int ibcregion = 0; ibcregion < nBCRegion; ibcregion ++)
    {
        StructBC *structBC = composite_bcregion->GetBCRegion(ibcregion);

        string bcName = structBC->GetBCName();
        int bctype = structBC->GetBCType();

        if (bctype == PHENGLEI::INTERFACE) continue;
        pair<int, string> bcNamePair(bctype, bcName);
        bcNamePairSet.insert(bcNamePair);
    }
}

void StructGrid::DumpWallFaceCenter(ActionKey *actkey)
{
    DataContainer *cdata = actkey->GetData();
    cdata->moveto_begin();

    std::ostringstream oss;

    int girdIndex = this->GetZoneID();
    if (girdIndex == 0)
    {
        int globalTotalWallFace = GlobalDataBase::GetIntParaFromDB("GlobalTotalWallFace");
        oss << globalTotalWallFace << "\n";
    }

    int nsolid_surface = this->GetNumberOfWallCell();
    if ((nsolid_surface == 0) || (level != 0))
    {
        string str = oss.str();
        cdata->write(const_cast <char *> (str.data()), str.size() * sizeof(char));
        return;
    }

    oss << setiosflags(ios::left);
    oss << setprecision(5);
    oss << setiosflags(ios::scientific);
    oss << setiosflags(ios::showpoint);

    CompositeStructBC *compositeBCRegion = this->GetCompositeBCRegion();
    int nBCRegion = compositeBCRegion->GetnBCRegion();

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        StructBC *structBC = compositeBCRegion->GetBCRegion(iBCRegion);
        int BCType = structBC->GetBCType();

        if (!IsWall(BCType))
        {
            continue;
        }

        int ist, ied, jst, jed, kst, ked;
        structBC->GetIJKRegion(ist, ied, jst, jed, kst, ked);
        int nSurface = structBC->GetFaceDirection() + 1;

        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    int iWall, jWall, kWall;
                    structBC->GetBoundaryFaceIndex(i, j, k, iWall, jWall, kWall);

                    RDouble xCenterWall, yCenterWall, zCenterWall;
                    this->FaceCoor(iWall, jWall, kWall, nSurface, xCenterWall, yCenterWall, zCenterWall);

                    oss << xCenterWall << "	";
                    oss << yCenterWall << "	";
                    oss << zCenterWall << "	";
                    oss <<"\n";
                }
            }
        }
    }

    string str = oss.str();
    cdata->write(const_cast <char *> (str.data()), str.size() * sizeof(char));
}

void CommunicateInterfaceValue(StructGrid *grid, RDouble4D *variable, const string &variableName, int variableDimension)
{
    //! Upload the interface data to the buffer.
    PHSPACE::UploadInterfaceValue(grid, variable, variableName, variableDimension);

    //! Download the interface data from the buffer.
    PHSPACE::DownloadInterfaceValue(grid, variable, variableName, variableDimension);
}

void CommunicateInterfaceValue(StructGrid *grid, RDouble3D *variable, const string &variableName)
{
    //! Upload the interface data to the buffer.
    PHSPACE::UploadInterfaceValue(grid, variable, variableName);

    //! Download the interface data from the buffer.
    PHSPACE::DownloadInterfaceValue(grid, variable, variableName);
}

void UploadInterfaceValue(StructGrid *grid, RDouble4D *f, const string &name, int neqn)
{
    InterfaceInfo *interfaceInfo = grid->GetInterfaceInfo();
    if (!interfaceInfo) return;

    StructGrid *fgrid = StructGridCast(grid->GetFinestGrid());
    InterfaceInfo *fiinfo = fgrid->GetInterfaceInfo();
    int level = grid->GetLevel();

    for (int ngst = GetNumberOfGhostCellLayers() - 1; ngst >= 0; -- ngst)
    {
        Data_ParamFieldSuite *dsend = interfaceInfo->GetSendDataStorage(ngst);
        RDouble **fg = reinterpret_cast<RDouble **> (dsend->GetDataPtr(name));
        if (fg == NULL)
        {
            cout<<name<<" is null"<<endl;
        }

        int nIFace = fiinfo->GetNIFace();
        
        for (int iface = 0; iface < nIFace; ++ iface)
        {
            int is, js, ks;
            fgrid->GetSourceIndexIJK(iface, ngst + 1, is, js, ks);
            fgrid->RemapMultigridIJK(level, is, js, ks);

            for (int m = 0; m < neqn; ++ m)
            {
                fg[m][iface] = (*f)(is, js, ks, m);
            }
        }
    }
}

void UploadInterfaceValue(StructGrid *grid, RDouble3D *f, const string &name)
{
    InterfaceInfo *interfaceInfo = grid->GetInterfaceInfo();
    if (!interfaceInfo) return;

    StructGrid *fgrid = StructGridCast(grid->GetFinestGrid());
    InterfaceInfo *fiinfo = fgrid->GetInterfaceInfo();
    int level = grid->GetLevel();

    int m = 0;
    for (int ngst = GetNumberOfGhostCellLayers() - 1; ngst >= 0; -- ngst)
    {
        Data_ParamFieldSuite *dsend = interfaceInfo->GetSendDataStorage(ngst);
        RDouble **fg = reinterpret_cast<RDouble **> (dsend->GetDataPtr(name));

        int nIFace = fiinfo->GetNIFace();
        
        for (int iface = 0; iface < nIFace; ++ iface)
        {
            int is, js, ks;
            fgrid->GetSourceIndexIJK(iface, ngst + 1, is, js, ks);
            fgrid->RemapMultigridIJK(level, is, js, ks);

            fg[m][iface] = (*f)(is, js, ks);
        }
    }
}

void DownloadInterfaceValue(StructGrid *grid, RDouble4D *f, const string &name, int neqn)
{
    InterfaceInfo *interfaceInfo = grid->GetInterfaceInfo();
    if (!interfaceInfo) return;

    StructGrid *fgrid = StructGridCast(grid->GetFinestGrid());
    InterfaceInfo *fiinfo = fgrid->GetInterfaceInfo();
    int level = grid->GetLevel();

    CompositeStructBC *composite_bcregion = fgrid->GetCompositeBCRegion();
    int periodicType = PHSPACE::GlobalDataBase::GetIntParaFromDB("periodicType");
    double rotationAngle = PHSPACE::GlobalDataBase::GetDoubleParaFromDB("rotationAngle");
    rotationAngle = rotationAngle * PI / 180;

    for (int ngst = GetNumberOfGhostCellLayers() - 1; ngst >= 0; -- ngst)
    {
        Data_ParamFieldSuite *drecv = interfaceInfo->GetRecvDataStorage(ngst);

        RDouble **fg = reinterpret_cast<RDouble **> (drecv->GetDataPtr(name));

        int nIFace = fiinfo->GetNIFace();
        for (int iface = 0; iface < nIFace; ++ iface)
        {
            int it, jt, kt;
            fgrid->GetTargetIndexIJK(iface, ngst+1, it, jt, kt);
            fgrid->RemapMultigridIJK(level, it, jt, kt);

            //! Rotate vector for periodic boundary at the end of communication for InterfaceValue.
            int *ibcregions = composite_bcregion->GetIFaceInfo();
            int ibcregion = ibcregions[iface];
            StructBC *bcregion = composite_bcregion->GetBCRegion(ibcregion);
            string bcName = bcregion->GetBCName();
            //! Rotate process only for vector in rotational periodicity.
            if (neqn >= 4 && periodicType == 2)
            {
                if (bcName == "Periodic_up" || bcName == "Periodic_down")
                {
                    RDouble rotfg[3] = {0.0, 0.0, 0.0};
                    if (bcName == "Periodic_up")
                    {
                        rotfg[0] = fg[1][iface];
                        rotfg[1] = fg[2][iface] * cos(2 * PI - rotationAngle) - fg[3][iface] * sin(2 * PI - rotationAngle);
                        rotfg[2] = fg[2][iface] * sin(2 * PI - rotationAngle) + fg[3][iface] * cos(2 * PI - rotationAngle);
                    }
                    else if (bcName == "Periodic_down")
                    {
                        rotfg[0] = fg[1][iface];
                        rotfg[1] = fg[2][iface] * cos(rotationAngle) - fg[3][iface] * sin(rotationAngle);
                        rotfg[2] = fg[2][iface] * sin(rotationAngle) + fg[3][iface] * cos(rotationAngle);
                    }

                    for (int m = 1; m <= 3; ++ m)
                    {
                        fg[m][iface] = rotfg[m-1];
                    }
                }
            }
            for (int m = 0; m < neqn; ++ m)
            {
                (*f)(it, jt, kt, m) = fg[m][iface];
            }
        }
    }
}

void DownloadInterfaceValue(StructGrid *grid, RDouble3D *f, const string &name)
{
    InterfaceInfo *interfaceInfo = grid->GetInterfaceInfo();
    if (!interfaceInfo) return;

    StructGrid *fgrid = StructGridCast(grid->GetFinestGrid());
    InterfaceInfo *fiinfo = fgrid->GetInterfaceInfo();
    int level = grid->GetLevel();

    int m = 0;
    for (int ngst = GetNumberOfGhostCellLayers() - 1; ngst >= 0; -- ngst)
    {
        Data_ParamFieldSuite *drecv = interfaceInfo->GetRecvDataStorage(ngst);

        RDouble **fg = reinterpret_cast<RDouble **> (drecv->GetDataPtr(name));
        int nIFace = fiinfo->GetNIFace();
        for (int iface = 0; iface < nIFace; ++ iface)
        {
            int it, jt, kt;
            fgrid->GetTargetIndexIJK(iface, ngst + 1, it, jt, kt);
            fgrid->RemapMultigridIJK(level, it, jt, kt);
            (*f)(it, jt, kt) = fg[m][iface];
        }
    }
}

void UploadOversetValue(StructGrid *grid, RDouble4D *field, const string &name, int neqn)
{
    StructGrid *fgrid = StructGridCast(grid->GetFinestGrid());
    OversetInformationProxy * oversetInformationProxy = grid->GetUnstrGrid()->GetOversetInformationProxy();
    OversetDataProxy * oversetDataProxy = oversetInformationProxy->GetOversetDataProxyForSend();

    int ni = fgrid->GetNI();
    int nj = fgrid->GetNJ();
    int nk = fgrid->GetNK();

    Data_ParamFieldSuite * dataStorage = oversetDataProxy->GetDataStorage();

    RDouble **fieldStorage = reinterpret_cast<RDouble **> (dataStorage->GetDataPtr(name));
    int numberOfNeighbors = oversetDataProxy->GetNumberOfNeighbors();

    for (int neighborIndex = 0; neighborIndex < numberOfNeighbors; ++ neighborIndex)
    {
        int   numberOfCells           = oversetDataProxy->GetNumberOfCells          (neighborIndex);
        int * cellIndexContainer      = oversetDataProxy->GetCellIndexContainer     (neighborIndex);
        int * cellStorageIndexMapping = oversetDataProxy->GetCellStorageIndexMapping(neighborIndex);

        for (int iCell = 0; iCell < numberOfCells; ++ iCell)
        {
            int storageAddress = cellStorageIndexMapping[iCell];
            int cellIndex = cellIndexContainer[iCell];
            int i, j, k;

            DecodeIJK(cellIndex, i, j, k, ni-1, nj-1, nk-1);

            for (int iEquation = 0; iEquation < neqn; ++ iEquation)
            {
                fieldStorage[ iEquation ][ storageAddress ] = (*field)(i, j, k, iEquation);
            }
        }
    }
}

void UploadOversetValue(StructGrid *grid, RDouble3D *field, const string &name)
{
    StructGrid *fgrid = StructGridCast(grid->GetFinestGrid());
    OversetInformationProxy * oversetInformationProxy = grid->GetUnstrGrid()->GetOversetInformationProxy();
    OversetDataProxy * oversetDataProxy = oversetInformationProxy->GetOversetDataProxyForSend();

    int ni = fgrid->GetNI();
    int nj = fgrid->GetNJ();
    int nk = fgrid->GetNK();

    Data_ParamFieldSuite * dataStorage = oversetDataProxy->GetDataStorage();

    RDouble **fieldStorage = reinterpret_cast<RDouble **> (dataStorage->GetDataPtr(name));
    int numberOfNeighbors = oversetDataProxy->GetNumberOfNeighbors();

    int iEquation = 0;
    for (int neighborIndex = 0; neighborIndex < numberOfNeighbors; ++ neighborIndex)
    {
        int   numberOfCells           = oversetDataProxy->GetNumberOfCells          (neighborIndex);
        int * cellIndexContainer      = oversetDataProxy->GetCellIndexContainer     (neighborIndex);
        int * cellStorageIndexMapping = oversetDataProxy->GetCellStorageIndexMapping(neighborIndex);

        for (int iCell = 0; iCell < numberOfCells; ++ iCell)
        {
            int storageAddress = cellStorageIndexMapping[iCell];
            int cellIndex = cellIndexContainer[iCell];
            int i, j, k;

            DecodeIJK(cellIndex, i, j, k, ni-1, nj-1, nk-1);

            fieldStorage[ iEquation ][ storageAddress ] = (*field)(i, j, k);
        }
    }
}

void DownloadOversetValue(StructGrid *grid, RDouble4D *field, const string &name, int neqn)
{
    if (field == 0) return;
    StructGrid *fgrid = StructGridCast(grid->GetFinestGrid());
    OversetInformationProxy * oversetInformationProxy = grid->GetUnstrGrid()->GetOversetInformationProxy();

    int gridLevelIndex = grid->GetLevel();

    int ni = fgrid->GetNI();
    int nj = fgrid->GetNJ();
    int nk = fgrid->GetNK();

    OversetDataProxy * oversetDataProxy = oversetInformationProxy->GetOversetDataProxyForReceive();

    Data_ParamFieldSuite * dataStorage = oversetDataProxy->GetDataStorage();

    RDouble **fieldStorage = reinterpret_cast<RDouble **> (dataStorage->GetDataPtr(name));

    int numberOfNeighborZones = oversetDataProxy->GetNumberOfNeighbors();

    for (int neighborIndex = 0; neighborIndex < numberOfNeighborZones; ++ neighborIndex)
    {
        int   numberOfCells           = oversetDataProxy->GetNumberOfCells          (neighborIndex);
        int * cellIndexContainer      = oversetDataProxy->GetCellIndexContainer     (neighborIndex);
        int * cellStorageIndexMapping = oversetDataProxy->GetCellStorageIndexMapping(neighborIndex);

        for (int iCell = 0; iCell < numberOfCells; ++ iCell)
        {
            int storageAddress = cellStorageIndexMapping[ iCell ];
            int cellIndex      = cellIndexContainer     [ iCell ];
            int i, j, k;

            DecodeIJK(cellIndex, i, j, k, ni-1, nj-1, nk-1);

            for (int iEquation = 0; iEquation < neqn; ++ iEquation)
            {
                (*field)(i, j, k, iEquation) = fieldStorage[ iEquation ][ storageAddress ];
            }
        }
    }
}

void DownloadOversetValue(StructGrid *grid, RDouble3D *field, const string &name)
{
    if (field == 0) return;
    StructGrid *fgrid = StructGridCast(grid->GetFinestGrid());
    OversetInformationProxy * oversetInformationProxy = grid->GetUnstrGrid()->GetOversetInformationProxy();

    int gridLevelIndex = grid->GetLevel();

    int ni = fgrid->GetNI();
    int nj = fgrid->GetNJ();
    int nk = fgrid->GetNK();

    OversetDataProxy * oversetDataProxy = oversetInformationProxy->GetOversetDataProxyForReceive();

    Data_ParamFieldSuite * dataStorage = oversetDataProxy->GetDataStorage();

    RDouble **fieldStorage = reinterpret_cast<RDouble **> (dataStorage->GetDataPtr(name));

    int numberOfNeighborZones = oversetDataProxy->GetNumberOfNeighbors();
    int iEquation = 0;
    for (int neighborIndex = 0; neighborIndex < numberOfNeighborZones; ++ neighborIndex)
    {
        int   numberOfCells           = oversetDataProxy->GetNumberOfCells          (neighborIndex);
        int * cellIndexContainer      = oversetDataProxy->GetCellIndexContainer     (neighborIndex);
        int * cellStorageIndexMapping = oversetDataProxy->GetCellStorageIndexMapping(neighborIndex);

        for (int iCell = 0; iCell < numberOfCells; ++ iCell)
        {
            int storageAddress = cellStorageIndexMapping[ iCell ];
            int cellIndex      = cellIndexContainer     [ iCell ];
            int i, j, k;

            DecodeIJK(cellIndex, i, j, k, ni-1, nj-1, nk-1);

            (*field)(i, j, k) = fieldStorage[ iEquation ][ storageAddress ];
        }
    }
}

void MinMaxDiffINF(StructGrid *grid, int ig, RDouble *dmin, RDouble *dmax, RDouble4D *q, int m)
{
    InterfaceInfo *interfaceInfo = grid->GetInterfaceInfo();
    if (!interfaceInfo) return;

    StructGrid *fgrid = StructGridCast(grid->GetFinestGrid());
    InterfaceInfo *fiinfo = fgrid->GetInterfaceInfo();

    int level = grid->GetLevel();
    int nIFace = fiinfo->GetNIFace();

    int ndim = GetDim();
    int id = 1;
    int jd = 1;
    int kd = 1;
    if (ndim == TWO_D) kd = 0;

    int idx[6], jdx[6], kdx[6];

    //! Find the maximum and minimum in the neighbor of each cell.
    for (int iface = 0; iface < nIFace; ++ iface)
    {
        int is, js, ks;
        fgrid->GetSourceIndexIJK(iface, ig+1, is, js, ks);
        fgrid->RemapMultigridIJK(level, is, js, ks);

        dmin[iface] = (*q)(is, js, ks, m);
        dmax[iface] = (*q)(is, js, ks, m);

        idx[0] = is + id;
        jdx[0] = js;
        kdx[0] = ks;

        idx[1] = is - id;
        jdx[1] = js;
        kdx[1] = ks;

        idx[2] = is;
        jdx[2] = js + jd;
        kdx[2] = ks;

        idx[3] = is;
        jdx[3] = js - jd;
        kdx[3] = ks;

        idx[4] = is;
        jdx[4] = js;
        kdx[4] = ks + kd;

        idx[5] = is;
        jdx[5] = js;
        kdx[5] = ks - kd;

        for (int n = 0; n < 6; ++ n)
        {
            dmin[iface] = MIN(dmin[iface], (*q)(idx[n], jdx[n], kdx[n], m));
            dmax[iface] = MAX(dmax[iface], (*q)(idx[n], jdx[n], kdx[n], m));
        }

        dmin[iface] -= (*q)(is, js, ks, m);
        dmax[iface] -= (*q)(is, js, ks, m);
    }
}

void BarthLimiterINF(StructGrid *grid, int ig, RDouble *limit, RDouble4D *q, RDouble *dqdx, RDouble *dqdy, RDouble *dqdz, int m)
{
    InterfaceInfo *interfaceInfo = grid->GetInterfaceInfo();
    if (!interfaceInfo) return;

    StructGrid *fgrid = StructGridCast(grid->GetFinestGrid());
    InterfaceInfo *fiinfo = fgrid->GetInterfaceInfo();

    int level = grid->GetLevel();
    int nIFace = fiinfo->GetNIFace();

    //! Need temporary arrays for the differences between the value in every cell and
    //! maximum/minimum in the neighboring cells.
    RDouble *dmin = new RDouble [nIFace];
    RDouble *dmax = new RDouble [nIFace];

    //! Find the the differences for q.
    MinMaxDiffINF(grid, ig, dmin, dmax, q, m);

    RDouble4D &gxfn = * (grid->GetFaceNormalX());
    RDouble4D &gyfn = * (grid->GetFaceNormalY());
    RDouble4D &gzfn = * (grid->GetFaceNormalZ());

    RDouble dx, dy, dz, ds, dot;
    RDouble dq_face;

    RDouble xcc, ycc, zcc, xfc, yfc, zfc, xfn, yfn, zfn;

    int ndim = GetDim();
    int id = 1;
    int jd = 1;
    int kd = 1;
    if (ndim == TWO_D) kd = 0;

    int ndmax = 6;
    if (ndim == TWO_D) ndmax = 4;

    int idx[6], jdx[6], kdx[6], ndx[6], ipn[6];

    //! Find the maximum and minimum in the neighbor of each cell.
    for (int iface = 0; iface < nIFace; ++ iface)
    {
        int is, js, ks;
        fgrid->GetSourceIndexIJK(iface, ig + 1, is, js, ks);
        fgrid->RemapMultigridIJK(level, is, js, ks);

        grid->CenterCoor(is, js, ks, xcc, ycc, zcc);

        idx[0] = is;
        jdx[0] = js;
        kdx[0] = ks;
        ndx[0] = 1;
        ipn[0] = -1;

        idx[1] = is + id;
        jdx[1] = js;
        kdx[1] = ks;
        ndx[1] = 1;
        ipn[1] = 1;

        idx[2] = is;
        jdx[2] = js;
        kdx[2] = ks;
        ndx[2] = 2;
        ipn[2] = -1;

        idx[3] = is;
        jdx[3] = js + jd;
        kdx[3] = ks;
        ndx[3] = 2;
        ipn[3] = 1;

        if (ndim == THREE_D)
        {
            idx[4] = is;
            jdx[4] = js;
            kdx[4] = ks;
            ndx[4] = 3;
            ipn[4] = -1;

            idx[5] = is;
            jdx[5] = js;
            kdx[5] = ks + kd;
            ndx[5] = 3;
            ipn[5] = 1;
        }

        for (int n = 0; n < ndmax; ++ n)
        {
            grid->FaceCoor(idx[n], jdx[n], kdx[n], ndx[n], xfc, yfc, zfc);

            xfn = gxfn(idx[n], jdx[n], kdx[n], ndx[n]);
            yfn = gyfn(idx[n], jdx[n], kdx[n], ndx[n]);
            zfn = gzfn(idx[n], jdx[n], kdx[n], ndx[n]);

            dx = xfc - xcc;
            dy = yfc - ycc;
            dz = zfc - zcc;
            dq_face = dqdx[iface] * dx + dqdy[iface] * dy + dqdz[iface] * dz;

            ds  = sqrt(dx * dx + dy * dy + dz * dz);
            dot = ipn[n] * (xfn * dx + yfn * dy + zfn * dz) / (ds + SMALL);

            if (dot < 0.0) limit[iface] = 0.0;

            if (dq_face > dmax[iface])
            {
                limit[iface] = MIN(limit[iface], dmax[iface] / dq_face);
            }
            else if (dq_face < dmin[iface])
            {
                limit[iface] = MIN(limit[iface], dmin[iface] / dq_face);
            }
        }

    }
    delete [] dmin;
    delete [] dmax;
}

void VencatLimiterINF(StructGrid *grid, int ig, RDouble *limit, RDouble4D *q, RDouble *dqdx, RDouble *dqdy, RDouble *dqdz, int m)
{
    InterfaceInfo *interfaceInfo = grid->GetInterfaceInfo();
    if (!interfaceInfo) return;

    StructGrid *fgrid = StructGridCast(grid->GetFinestGrid());
    InterfaceInfo *fiinfo = fgrid->GetInterfaceInfo();

    int level = grid->GetLevel();
    int nIFace = fiinfo->GetNIFace();

    //! Need temporary arrays for the differences between the value in every cell and
    //! maximum/minimum in the neighboring cells.
    RDouble *dmin = new RDouble [nIFace];
    RDouble *dmax = new RDouble [nIFace];

    // Find the the differences for q
    MinMaxDiffINF(grid, ig, dmin, dmax, q, m);

    RDouble4D &gxfn = * (grid->GetFaceNormalX());
    RDouble4D &gyfn = * (grid->GetFaceNormalY());
    RDouble4D &gzfn = * (grid->GetFaceNormalZ());

    RDouble dx, dy, dz, ds, dot;
    RDouble dq_face, eps, tmp;

    RDouble xcc, ycc, zcc, xfc, yfc, zfc, xfn, yfn, zfn;

    RDouble venkatCoeff = 1.0e-5;
    GlobalDataBase::GetData("venkatCoeff", &venkatCoeff, PHDOUBLE, 1);

    int ndim = GetDim();
    int id = 1;
    int jd = 1;
    int kd = 1;
    if (ndim == TWO_D) kd = 0;

    int ndmax = 6;
    if (ndim == TWO_D) ndmax = 4;

    int idx[6], jdx[6], kdx[6], ndx[6], ipn[6];

    //! Find the maximum and minimum in the neighbor of each cell.
    for (int iface = 0; iface < nIFace; ++ iface)
    {
        int is, js, ks;
        fgrid->GetSourceIndexIJK(iface, ig + 1, is, js, ks);
        fgrid->RemapMultigridIJK(level, is, js, ks);

        grid->CenterCoor(is, js, ks, xcc, ycc, zcc);

        idx[0] = is;
        jdx[0] = js;
        kdx[0] = ks;
        ndx[0] = 1;
        ipn[0] = -1;

        idx[1] = is + id;
        jdx[1] = js;
        kdx[1] = ks;
        ndx[1] = 1;
        ipn[1] = 1;

        idx[2] = is;
        jdx[2] = js;
        kdx[2] = ks;
        ndx[2] = 2;
        ipn[2] = -1;

        idx[3] = is;
        jdx[3] = js + jd;
        kdx[3] = ks;
        ndx[3] = 2;
        ipn[3] = 1;

        if (ndim == THREE_D)
        {
            idx[4] = is;
            jdx[4] = js;
            kdx[4] = ks;
            ndx[4] = 3;
            ipn[4] = -1;

            idx[5] = is;
            jdx[5] = js;
            kdx[5] = ks + kd;
            ndx[5] = 3;
            ipn[5] = 1;
        }

        for (int n = 0; n < ndmax; ++ n)
        {
            grid->FaceCoor(idx[n], jdx[n], kdx[n], ndx[n], xfc, yfc, zfc);

            xfn = gxfn(idx[n], jdx[n], kdx[n], ndx[n]);
            yfn = gyfn(idx[n], jdx[n], kdx[n], ndx[n]);
            zfn = gzfn(idx[n], jdx[n], kdx[n], ndx[n]);

            dx = xfc - xcc;
            dy = yfc - ycc;
            dz = zfc - zcc;
            dq_face = dqdx[iface] * dx + dqdy[iface] * dy + dqdz[iface] * dz;

            ds = sqrt(dx * dx + dy * dy + dz * dz);
            dot = ipn[n] * (xfn * dx + yfn * dy + zfn * dz) / (ds + SMALL);

            if (dot < 0.0) limit[iface] = 0.0;

            eps = venkatCoeff * ds;
            eps = eps * eps * eps;

            if (dq_face > SMALL)
            {
                tmp = VenFun(dmax[iface], dq_face, eps);
                limit[iface] = MIN(limit[iface], tmp);
            }
            else if (dq_face < - SMALL)
            {
                tmp = VenFun(dmin[iface], dq_face, eps);
                limit[iface] = MIN(limit[iface], tmp);
            }
        }
    }
    delete [] dmin;
    delete [] dmax;
}

void EncodeIJK(int &index, int i, int j, int k, int ni, int nj, int nk)
{
    index = (i - 1) + (j - 1) * ni + (k - 1) * (ni * nj);
}

void DecodeIJK(int index, int &i, int &j, int &k, int ni, int nj, int nk)
{
    k = index / (ni * nj) + 1;
    index -= (k - 1) * ni * nj;
    j = index / ni + 1;
    i = index - (j - 1) * ni + 1;
}

void GetIJKRegion(const Range &I, const Range &J, const Range &K, int &ist, int &ied, int &jst, int &jed, int &kst, int &ked)
{
    ist = I.first();
    ied = I.last();
    jst = J.first();
    jed = J.last();
    kst = K.first();
    ked = K.last();
}

void GetNsurfIndex(int nsurf, int &i, int &j, int &k)
{
    i = mide[0][nsurf-1];
    j = mide[1][nsurf-1];
    k = mide[2][nsurf-1];
}

void GetRange(int ni, int nj, int nk, int startShift, int endShift, Range &I, Range &J, Range &K)
{
    I.setRange(1 + startShift, ni + endShift);
    J.setRange(1 + startShift, nj + endShift);
    K.setRange(1 + startShift, nk + endShift);

    if (nk == 1)
    {
        K.setRange(1, 1);
    }

    return;
}

void GetRange(Range &I, Range &J, Range &K, int ist, int ied, int jst, int jed, int kst, int ked)
{
    I.setRange(ist, ied);
    J.setRange(jst, jed);
    K.setRange(kst, ked);

    return;
}

void StructGrid::SetOrdinaryGridIndex(int ordinaryGridIndexIn)
{
    this->ordinaryGridIndex = ordinaryGridIndexIn;
}

void StructGrid::SetOrdinaryDimStartIndex(int *ordinaryDimStartIndexIn)
{
    for (int iDim = 0; iDim < THREE_D; ++ iDim)
    {
        this->ordinaryDimStartIndex[iDim] = ordinaryDimStartIndexIn[iDim];
    }
}

void StructGrid::SetOrdinaryDimEndIndex(int *ordinaryDimEndIndexIn)
{
    for (int iDim = 0; iDim < THREE_D; ++ iDim)
    {
        this->ordinaryDimEndIndex[iDim] = ordinaryDimEndIndexIn[iDim];
    }
}

int StructGrid::GetOrdinaryGridIndex()
{
    return this->ordinaryGridIndex;
}

int *StructGrid::GetOrdinaryDimStartIndex()
{
    return ordinaryDimStartIndex;
}

int *StructGrid::GetOrdinaryDimEndIndex()
{
    return ordinaryDimEndIndex;
}

void StructGrid::GetLocalIndex(int gp, int &ip, int &jp, int &kp) const
{
    int zoneStartPointLabel = GetZoneStartPointLabel();
    int ni = nodeTopology->GetNI();
    int nj = nodeTopology->GetNJ();
    int mm = ni * nj;
    int ms = gp - zoneStartPointLabel;
    kp = ms / mm;
    int other = ms % mm;
    jp = other / ni;
    ip = other % ni;
}

LIB_EXPORT RDouble3D * StructGrid::GetLargestLocalGridLength()
{
    if (!cellMetrics->GetLargestLocalGridLength())
    {
        ComputeLargestLocalGridLength();
    }
    return cellMetrics->GetLargestLocalGridLength();
}

LIB_EXPORT RDouble3D * StructGrid::GetSmallestLocalGridLength()
{
    if (!cellMetrics->GetSmallestLocalGridLength())
    {
        ComputeSmallestLocalGridLength();
    }
    return cellMetrics->GetSmallestLocalGridLength();
}

LIB_EXPORT RDouble3D * StructGrid::GetSubgridLength()
{
    if (!cellMetrics->GetSubgridLength())
    {
        ComputeSubgridLength();
    }
    return cellMetrics->GetSubgridLength();
}

void StructGrid::ComputeLargestLocalGridLength()
{
    int ni = this->GetNI();
    int nj = this->GetNJ();
    int nk = this->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, 0, -1, I, J, K);

    int ist, ied, jst, jed, kst, ked;
    this->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    RDouble3D *xcc = this->GetCellCenterX();
    RDouble3D *ycc = this->GetCellCenterY();
    RDouble3D *zcc = this->GetCellCenterZ();

    RDouble3D *largestLocalGridLength = new RDouble3D(I, J, K, fortranArray);
    PHSPACE::SetField(largestLocalGridLength, -LARGE, ist, ied, jst, jed, kst, ked);

    RDouble dist = 0.0;
    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                if (i < ied)
                {
                    dist = DISTANCE((*xcc)(i, j, k) - (*xcc)(i+1, j, k), (*ycc)(i, j, k) - (*ycc)(i+1, j, k), (*zcc)(i, j, k) - (*zcc)(i+1, j, k));
                    (*largestLocalGridLength)(i,   j, k) = MAX((*largestLocalGridLength)(i,   j, k), dist);
                    (*largestLocalGridLength)(i+1, j, k) = MAX((*largestLocalGridLength)(i+1, j, k), dist);
                }
                if (j < jed)
                {
                    dist = DISTANCE((*xcc)(i, j, k) - (*xcc)(i, j+1, k), (*ycc)(i, j, k) - (*ycc)(i, j+1, k), (*zcc)(i, j, k) - (*zcc)(i, j+1, k));
                    (*largestLocalGridLength)(i, j,   k) = MAX((*largestLocalGridLength)(i, j,   k), dist);
                    (*largestLocalGridLength)(i, j+1, k) = MAX((*largestLocalGridLength)(i, j+1, k), dist);
                }
                if (k < ked)
                {
                    dist = DISTANCE((*xcc)(i, j, k) - (*xcc)(i, j, k+1), (*ycc)(i, j, k) - (*ycc)(i, j, k+1), (*zcc)(i, j, k) - (*zcc)(i, j, k+1));
                    (*largestLocalGridLength)(i, j,   k) = MAX((*largestLocalGridLength)(i, j,   k), dist);
                    (*largestLocalGridLength)(i, j, k+1) = MAX((*largestLocalGridLength)(i, j, k+1), dist);
                }
            }
        }
    }
    this->cellMetrics->SetLargestLocalGridLength(largestLocalGridLength);
}

void StructGrid::ComputeSmallestLocalGridLength()
{
    int ni = this->GetNI();
    int nj = this->GetNJ();
    int nk = this->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, 0, -1, I, J, K);

    int ist, ied, jst, jed, kst, ked;
    this->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    RDouble3D *xcc = this->GetCellCenterX();
    RDouble3D *ycc = this->GetCellCenterY();
    RDouble3D *zcc = this->GetCellCenterZ();

    RDouble3D *smallestLocalGridLength = new RDouble3D(I, J, K, fortranArray);
    PHSPACE::SetField(smallestLocalGridLength, LARGE, ist, ied, jst, jed, kst, ked);

    RDouble dist = 0.0;
    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                if (i < ied)
                {
                    dist = DISTANCE((*xcc)(i, j, k) - (*xcc)(i+1, j, k), (*ycc)(i, j, k) - (*ycc)(i+1, j, k), (*zcc)(i, j, k) - (*zcc)(i+1, j, k));
                    (*smallestLocalGridLength)(i,   j, k) = MIN((*smallestLocalGridLength)(i,   j, k), dist);
                    (*smallestLocalGridLength)(i+1, j, k) = MIN((*smallestLocalGridLength)(i+1, j, k), dist);
                }
                if (j < jed)
                {
                    dist = DISTANCE((*xcc)(i, j, k) - (*xcc)(i, j+1, k), (*ycc)(i, j, k) - (*ycc)(i, j+1, k), (*zcc)(i, j, k) - (*zcc)(i, j+1, k));
                    (*smallestLocalGridLength)(i, j,   k) = MIN((*smallestLocalGridLength)(i, j,   k), dist);
                    (*smallestLocalGridLength)(i, j+1, k) = MIN((*smallestLocalGridLength)(i, j+1, k), dist);
                }
                if (k < ked)
                {
                    dist = DISTANCE((*xcc)(i, j, k) - (*xcc)(i, j, k+1), (*ycc)(i, j, k) - (*ycc)(i, j, k+1), (*zcc)(i, j, k) - (*zcc)(i, j, k+1));
                    (*smallestLocalGridLength)(i, j,   k) = MIN((*smallestLocalGridLength)(i, j,   k), dist);
                    (*smallestLocalGridLength)(i, j, k+1) = MIN((*smallestLocalGridLength)(i, j, k+1), dist);
                }
            }
        }
    }
    this->cellMetrics->SetSmallestLocalGridLength(smallestLocalGridLength);
}

void StructGrid::ComputeSubgridLength()
{
    int ni = this->GetNI();
    int nj = this->GetNJ();
    int nk = this->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, 0, -1, I, J, K);

    int ist, ied, jst, jed, kst, ked;
    this->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    RDouble3D *walldist = this->GetWallDist();
    RDouble3D *largestLocalGridLength = this->GetLargestLocalGridLength();
    RDouble3D *smallestLocalGridLength = this->GetSmallestLocalGridLength();

    RDouble3D *subgridLength = new RDouble3D(I, J, K, fortranArray);

    const RDouble cw = 0.15;
    RDouble wallDistance, largestLocalGridSpacing, smallestLocalGridSpacing, subgridLengthScale;
    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                wallDistance = (*walldist)(i, j, k);
                largestLocalGridSpacing = (*largestLocalGridLength)(i, j, k);
                smallestLocalGridSpacing = (*smallestLocalGridLength)(i, j, k);

                //! The change ratio of the normal distance is ignored bellow!!! Re-consider it future!
                subgridLengthScale = MIN(MAX(cw * MAX(wallDistance, largestLocalGridSpacing), smallestLocalGridSpacing), largestLocalGridSpacing);

                (*subgridLength)(i, j, k) = subgridLengthScale;
            }
        }
    }

    this->cellMetrics->SetSubgridLength(subgridLength);
}

void StructGrid::InitVariableWallTemperature()
{
    fstream file;
    string wallTemperatureFile = "";
    GlobalDataBase::GetData("wallTemperaturefile", &wallTemperatureFile, PHSTRING, 1);
    OpenFile(file, wallTemperatureFile, ios_base::in|ios_base::binary);

    int dataNumber = 0;
    PHRead(file, dataNumber);

    int varNumber = 0;
    PHRead(file, varNumber);

    RDouble *pmin = this->GetMinBox();
    RDouble *pmax = this->GetMaxBox();
    vector < vector <RDouble> > dataInCurrentGrid;
    dataInCurrentGrid.resize(0);

    for (int iData = 0; iData < dataNumber; ++ iData)
    {
        float xCoor, yCoor, zCoor;
        float tWall;

        PHRead(file, xCoor);
        PHRead(file, yCoor);
        PHRead(file, zCoor);
        PHRead(file, tWall);

        RDouble coordX = (double)xCoor;
        RDouble coordY = (double)yCoor;
        RDouble coordZ = (double)zCoor;

        if (coordX < pmin[0] || coordX > pmax[0]) continue;
        if (coordY < pmin[1] || coordY > pmax[1]) continue;
        if (coordZ < pmin[2] || coordZ > pmax[2]) continue;

        vector <RDouble> currentData;
        currentData.push_back(coordX);
        currentData.push_back(coordY);
        currentData.push_back(coordZ);
        currentData.push_back((double)tWall);

        dataInCurrentGrid.push_back(currentData);
    }

    CloseFile(file);

    int dataNumberInCurrentGrid = dataInCurrentGrid.size();
    if (dataNumberInCurrentGrid == 0)
    {
        return;
    }

    int *dataUsed = new int[dataNumberInCurrentGrid];
    std::fill_n(dataUsed, dataNumberInCurrentGrid, 0);

    RDouble3D &structx = * this->GetStructX();
    RDouble3D &structy = * this->GetStructY();
    RDouble3D &structz = * this->GetStructZ();
    RDouble4D &faceArea= *(this->GetFaceArea());

    CompositeStructBC *compositeBCRegion = this->GetCompositeBCRegion();
    int nBCRegion = compositeBCRegion->GetnBCRegion();

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        StructBC *structBC = compositeBCRegion->GetBCRegion(iBCRegion);
        int BCType = structBC->GetBCType();

        if (!IsWall(BCType) && BCType != PHENGLEI::ABLATION_SURFACE)
        {
            continue;
        }

        int ist, ied, jst, jed, kst, ked;
        structBC->GetIJKRegion(ist, ied, jst, jed, kst, ked);
        Range I(ist, ied);
        Range J(jst, jed);
        Range K(kst, ked);
        RDouble3D *wallTempArray = new RDouble3D(I, J, K, fortranArray);

        int nk = GetNK();
        int il1 = 1;
        int jl1 = 1;
        int kl1 = 1;
        if (nk == 1) kl1 = 0;
        int faceDirection = structBC->GetFaceDirection() + 1;

        int il, jl, kl;

        int findData = 0;
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    bool flag = false;

                    int iWall, jWall, kWall;
                    structBC->GetBoundaryFaceIndex(i, j, k, iWall, jWall, kWall);
                    RDouble surfaceArea = faceArea(iWall, jWall, kWall, faceDirection);

                    il = i + il1;
                    jl = j + jl1;
                    kl = k + kl1;
                    vector <vector <RDouble> > nodeCoor;
                    if (faceDirection == 1)
                    {
                        vector <RDouble> oneNode;
                        oneNode.push_back(structx(i, j, k));
                        oneNode.push_back(structy(i, j, k));
                        oneNode.push_back(structz(i, j, k));
                        nodeCoor.push_back(oneNode);

                        oneNode.resize(0);
                        oneNode.push_back(structx(i, jl, k));
                        oneNode.push_back(structy(i, jl, k));
                        oneNode.push_back(structz(i, jl, k));
                        nodeCoor.push_back(oneNode);

                        oneNode.resize(0);
                        oneNode.push_back(structx(i, jl, kl));
                        oneNode.push_back(structy(i, jl, kl));
                        oneNode.push_back(structz(i, jl, kl));
                        nodeCoor.push_back(oneNode);

                        oneNode.resize(0);
                        oneNode.push_back(structx(i, j, kl));
                        oneNode.push_back(structy(i, j, kl));
                        oneNode.push_back(structz(i, j, kl));
                        nodeCoor.push_back(oneNode);
                    }
                    else if (faceDirection == 2)
                    {
                        vector <RDouble> oneNode;
                        oneNode.push_back(structx(i, j, k));
                        oneNode.push_back(structy(i, j, k));
                        oneNode.push_back(structz(i, j, k));
                        nodeCoor.push_back(oneNode);

                        oneNode.resize(0);
                        oneNode.push_back(structx(il, j, k));
                        oneNode.push_back(structy(il, j, k));
                        oneNode.push_back(structz(il, j, k));
                        nodeCoor.push_back(oneNode);

                        oneNode.resize(0);
                        oneNode.push_back(structx(il, j, kl));
                        oneNode.push_back(structy(il, j, kl));
                        oneNode.push_back(structz(il, j, kl));
                        nodeCoor.push_back(oneNode);

                        oneNode.resize(0);
                        oneNode.push_back(structx(i, j, kl));
                        oneNode.push_back(structy(i, j, kl));
                        oneNode.push_back(structz(i, j, kl));
                        nodeCoor.push_back(oneNode);
                    }
                    else if (faceDirection == 3)
                    {
                        vector <RDouble> oneNode;
                        oneNode.push_back(structx(i, j, k));
                        oneNode.push_back(structy(i, j, k));
                        oneNode.push_back(structz(i, j, k));
                        nodeCoor.push_back(oneNode);

                        oneNode.resize(0);
                        oneNode.push_back(structx(il, j, k));
                        oneNode.push_back(structy(il, j, k));
                        oneNode.push_back(structz(il, j, k));
                        nodeCoor.push_back(oneNode);

                        oneNode.resize(0);
                        oneNode.push_back(structx(il, jl, k));
                        oneNode.push_back(structy(il, jl, k));
                        oneNode.push_back(structz(il, jl, k));
                        nodeCoor.push_back(oneNode);

                        oneNode.resize(0);
                        oneNode.push_back(structx(i, jl, k));
                        oneNode.push_back(structy(i, jl, k));
                        oneNode.push_back(structz(i, jl, k));
                        nodeCoor.push_back(oneNode);
                    }

                    for (int iData = 0; iData < dataNumberInCurrentGrid; ++ iData)
                    {
                        if (dataUsed[iData])
                        {
                            continue;
                        }

                        RDouble coordX = dataInCurrentGrid[iData][0];
                        RDouble coordY = dataInCurrentGrid[iData][1];
                        RDouble coordZ = dataInCurrentGrid[iData][2];

                        vector <vector <RDouble> > lineVec;
                        lineVec.resize(4);

                        for (int iNode = 0; iNode < 4; ++ iNode)
                        {
                            RDouble dx = nodeCoor[iNode][0] - coordX;
                            RDouble dy = nodeCoor[iNode][1] - coordY;
                            RDouble dz = nodeCoor[iNode][2] - coordZ;

                            lineVec[iNode].push_back(dx);
                            lineVec[iNode].push_back(dy);
                            lineVec[iNode].push_back(dz);
                        }

                        RDouble areaSum = 0.0;
                        for (int iTri = 0; iTri < 4; ++ iTri)
                        {
                            int lineIndex1 =  iTri;
                            int lineIndex2 = (iTri + 1) % 4;

                            RDouble p[3];
                            CrossProduct(&lineVec[lineIndex1][0], &lineVec[lineIndex2][0], p);
                            RDouble areaCurrent = DISTANCE(p[0], p[1], p[2]);

                            areaSum += ABS(areaCurrent);
                        }

                        areaSum = half * areaSum;

                        RDouble eps = (areaSum - surfaceArea) / surfaceArea;
                        if (ABS(eps) < 1e-4)
                        {
                            flag = true;
                        }

                        if (flag)
                        {
                            dataUsed[iData] = 1;
                            RDouble tWall = dataInCurrentGrid[iData][3];
                            (*wallTempArray)(i, j, k) = tWall;
                            findData ++;
                            break;
                        }
                    }
                }
            }
        }

        /*int dataSize = (ied - ist + 1)
                     * (jed - jst + 1)
                     * (ked - kst + 1);

        if (dataSize != findData)
        {
            WriteLogFile("Error !!!!");
        }*/

        structBC->UpdateFieldDataPtr("wallTempArray" , wallTempArray);
    }

    delete []dataUsed;
}

void CompNodeVar(Grid *grid_in, RDouble4D &qn, int m, RDouble4D &q, int n, bool isVelocityForPostVisual)
{
    StructGrid *grid = StructGridCast(grid_in);
    int nk = grid->GetNK();

    int ist, ied, jst, jed, kst, ked;
    grid->GetNodeIterationIndex(ist, ied, jst, jed, kst, ked);

    int il1 = 1;
    int jl1 = 1;
    int kl1 = 1;
    if (nk == 1) kl1 = 0;

    int il, jl, kl;

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                il = i - il1;
                jl = j - jl1;
                kl = k - kl1;
                qn(i, j, k, m) = eighth * (q(il, jl, kl, n) + q(il, j, kl, n) + q(i, jl, kl, n) + q(i, j, kl, n) +
                                           q(il, jl, k,  n) + q(il, j, k , n) + q(i, jl, k,  n) + q(i, j, k,  n));
            }
        }
    }

    if (!isVelocityForPostVisual)
    {
        return;
    }

    CompositeStructBC *composite_bcregion = grid->GetCompositeBCRegion();
    int nBCRegion = composite_bcregion->GetnBCRegion();
    for (int ibcregion = 0; ibcregion < nBCRegion; ++ ibcregion)
    {
        StructBC *bcregion = composite_bcregion->GetBCRegion(ibcregion);
        int bctype = bcregion->GetBCType();

        if (!IsWall(bctype) && bctype != PHENGLEI::ABLATION_SURFACE)
        {
            continue;
        }

        bcregion->GetIJKRegion(ist, ied, jst, jed, kst, ked);

        int nsurf = bcregion->GetFaceDirection() + 1;

        int inode, jnode, knode;
        for (int k = kst; k <= ked; ++ k)
        {
            for (int j = jst; j <= jed; ++ j)
            {
                for (int i = ist; i <= ied; ++ i)
                {
                    bcregion->GetBoundaryFaceIndex(i, j, k, inode, jnode, knode);

                    int il, jl, kl;
                    il = inode + il1;
                    jl = jnode + jl1;
                    kl = knode + kl1;

                    if (nsurf == 1)
                    {
                        qn(inode, jnode, knode, m) = q(i, j, k, n);
                        qn(inode, jl   , knode, m) = q(i, j, k, n);
                        qn(inode, jl   , kl   , m) = q(i, j, k, n);
                        qn(inode, jnode, kl   , m) = q(i, j, k, n);
                    }
                    else if (nsurf == 2)
                    {
                        qn(inode, jnode, knode, m) = q(i, j, k, n);
                        qn(il   , jnode, knode, m) = q(i, j, k, n);
                        qn(il   , jnode, kl   , m) = q(i, j, k, n);
                        qn(inode, jnode, kl   , m) = q(i, j, k, n);
                    }
                    else if (nsurf == 3)
                    {
                        qn(inode, jnode, knode, m) = q(i, j, k, n);
                        qn(il   , jnode, knode, m) = q(i, j, k, n);
                        qn(il   , jl   , knode, m) = q(i, j, k, n);
                        qn(inode, jl   , knode, m) = q(i, j, k, n);
                    }
                }
            }
        }
    }
}

void CompNodeVar(Grid *grid_in, RDouble4D &qn, int m, RDouble3D &q)
{
    StructGrid *grid = StructGridCast(grid_in);
    int nk = grid->GetNK();

    int ist, ied, jst, jed, kst, ked;
    grid->GetNodeIterationIndex(ist, ied, jst, jed, kst, ked);

    int il1 = 1;
    int jl1 = 1;
    int kl1 = 1;
    if (nk == 1) kl1 = 0;

    int il, jl, kl;

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                il = i - il1;
                jl = j - jl1;
                kl = k - kl1;
                qn(i, j, k, m) = eighth * (q(il, jl, kl) + q(il, j, kl) + q(i, jl, kl) + q(i, j, kl) +
                                           q(il, jl, k) + q(il, j, k) + q(i, jl, k) + q(i, j, k));
            }
        }
    }
}

void ModifyWallTemperature(Grid *grid_in, RDouble4D &qn)
{
    if (GetDim() == TWO_D)
    {
        return;
    }

    StructGrid *grid = StructGridCast(grid_in);

    CompositeStructBC *compositeBCRegion = grid->GetCompositeBCRegion();
    int nBCRegion = compositeBCRegion->GetnBCRegion();

    for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
    {
        StructBC *structBC = compositeBCRegion->GetBCRegion(iBCRegion);
        int BCType = structBC->GetBCType();

        if (!(IsWall(BCType)))
        {
            continue;
        }

        RDouble wallTemperature;
        structBC->GetParamData("wallTemperature", &wallTemperature, PHDOUBLE, 1);

        if (wallTemperature < 0.0)
        {
            continue;
        }

        int *s_st = structBC->GetStartPoint();
        int *s_ed = structBC->GetEndPoint();
        int *s_lr3d = structBC->GetFaceDirectionIndex();

        int nodeIndexStart[3];
        int nodeIndexEnd[3];

        for (int iDim = 0; iDim < 3; ++ iDim)
        {
            nodeIndexStart[iDim] = s_st[iDim];
            nodeIndexEnd  [iDim] = s_ed[iDim];
        }

        for (int m = 0; m < GetDim(); ++ m)
        {
            if (structBC->GetFaceDirection() == m)
            {
                if (s_st[m] > 1 || s_lr3d[m] == 1)
                {
                    nodeIndexStart[m] += 1;
                    nodeIndexEnd  [m] += 1;
                }
            }
            else
            {
                nodeIndexEnd[m] += 1;
            }
        }

        int ist = nodeIndexStart[0];
        int ied = nodeIndexEnd  [0];

        int jst = nodeIndexStart[1];
        int jed = nodeIndexEnd  [1];

        int kst = nodeIndexStart[2];
        int ked = nodeIndexEnd  [2];

        int nsurf = structBC->GetFaceDirection();
        if (nsurf == 0)
        {
            for (int iNode = jst; iNode <= jed; ++ iNode)
            {
                qn(ist, iNode, kst, 0) = qn(ist, iNode, kst + 1, 0);
                qn(ist, iNode, ked, 0) = qn(ist, iNode, ked - 1, 0);
            }

            for (int iNode = kst; iNode <= ked; ++ iNode)
            {
                qn(ist, jst, iNode, 0) = qn(ist, jst + 1, iNode, 0);
                qn(ist, jed, iNode, 0) = qn(ist, jed - 1, iNode, 0);
            }

            qn(ist, jst, kst, 0) = qn(ist, jst + 1, kst + 1, 0);
            qn(ist, jed, kst, 0) = qn(ist, jed - 1, kst + 1, 0);
            qn(ist, jst, ked, 0) = qn(ist, jst + 1, ked - 1, 0);
            qn(ist, jed, ked, 0) = qn(ist, jed - 1, ked - 1, 0);
        }

        if (nsurf == 1)
        {
            for (int iNode = ist; iNode <= ied; ++ iNode)
            {
                qn(iNode, jst, kst, 0) = qn(iNode, jst, kst + 1, 0);
                qn(iNode, jst, ked, 0) = qn(iNode, jst, ked - 1, 0);
            }

            for (int iNode = kst; iNode <= ked; ++ iNode)
            {
                qn(ist, jst, iNode, 0) = qn(ist + 1, jst, iNode, 0);
                qn(ied, jst, iNode, 0) = qn(ied - 1, jst, iNode, 0);
            }

            qn(ist, jst, kst, 0) = qn(ist + 1, jst, kst + 1, 0);
            qn(ied, jst, kst, 0) = qn(ied - 1, jst, kst + 1, 0);
            qn(ist, jst, ked, 0) = qn(ist + 1, jst, ked - 1, 0);
            qn(ied, jst, ked, 0) = qn(ied - 1, jst, ked - 1, 0);
        }

        if (nsurf == 2)
        {
            for (int iNode = ist; iNode <= ied; ++ iNode)
            {
                qn(iNode, jst, kst, 0) = qn(iNode, jst + 1, kst, 0);
                qn(iNode, jed, kst, 0) = qn(iNode, jed - 1, kst, 0);
            }

            for (int iNode = jst; iNode <= jed; ++ iNode)
            {
                qn(ist, iNode, kst, 0) = qn(ist + 1, iNode, kst, 0);
                qn(ied, iNode, kst, 0) = qn(ied - 1, iNode, kst, 0);
            }

            qn(ist, jst, kst, 0) = qn(ist + 1, jst + 1, kst, 0);
            qn(ied, jst, kst, 0) = qn(ied - 1, jst + 1, kst, 0);
            qn(ist, jed, kst, 0) = qn(ist + 1, jed - 1, kst, 0);
            qn(ied, jed, kst, 0) = qn(ied - 1, jed - 1, kst, 0);
        }
    }
}

void GetIJKROT(int ijk[4][3], int *irot, int *jrot, int *krot)
{
    ijk[0][0] = 0;
    ijk[0][1] = 0;
    ijk[0][2] = 0;

    ijk[1][0] = 0;
    ijk[1][1] = 1;
    ijk[1][2] = 0;

    ijk[2][0] = 0;
    ijk[2][1] = 1;
    ijk[2][2] = 1;

    ijk[3][0] = 0;
    ijk[3][1] = 0;
    ijk[3][2] = 1;

    irot[0] = 0;
    jrot[0] = 1;
    krot[0] = 2;

    irot[1] = 1;
    jrot[1] = 0;
    krot[1] = 2;

    irot[2] = 1;
    jrot[2] = 2;
    krot[2] = 0;
}

}