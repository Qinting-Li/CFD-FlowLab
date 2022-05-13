#include "Geo_StructGrid.h"
#include "IO_FileName.h"
#include "Constants.h"
#include "GridType.h"
#include "Geo_StructBC.h"
#include "Geo_Interface.h"
#include "Geo_Interpoint.h"
#include "LinkStruct.h"
#include "TK_Exit.h"
#include "OversetInformation.h"
#include "PHIO.h"

using namespace std;

namespace PHSPACE
{

Grid::Grid()
{
    zoneProbesNumber = 0;
    zoneProbesGlobalID.resize(0);
    zoneProbesLineID.resize(0);
    zoneProbesSurfaceID.resize(0);
    zoneProbesCoordinates.resize(0);
}

LIB_EXPORT void Grid::InitGrid(GridID *index, int level, int dim, int type)
{
    this->dimension = dim;
    this->type  = type;
    this->level = level;
    this->index = index;

    cGrid = NULL;
    fGrid = NULL;

    interfaceInfo = NULL;
    interpointInformation = NULL;
    interfaceFields = NULL;
    interpointFields = NULL;
    oversetInfoProxy = NULL;
    oversetInformationProxy = NULL;

    own_database = 0;
    gField = new Data_Field();
    gPara = NULL;

    oldGrid = NULL;
    fileIndexCurrentGridBelong = 0;
}

Grid::~Grid()
{
    delete index;    index = NULL;
    delete oversetInfoProxy;    oversetInfoProxy = NULL;
    delete gField;    gField = NULL;

    if (own_database)
    {
        delete gPara;    gPara = NULL;
    }

    delete oldGrid;    oldGrid = NULL;

    //if (IsFinestGrid())
    //{
    //    FreeInterfaceField();
    //}

    if (!interfaceFields){delete interfaceFields;    interfaceFields = NULL;};
    delete interpointFields;    interpointFields = NULL;

    delete oversetInformationProxy;    oversetInformationProxy = NULL;

    if (!interfaceInfo){delete interfaceInfo;    interfaceInfo = NULL;};
    delete interpointInformation;    interpointInformation = NULL;

    zoneProbesGlobalID.clear();
    zoneProbesLineID.clear();
    zoneProbesSurfaceID.clear();
    zoneProbesCoordinates.clear();
}

LIB_EXPORT Grid * Grid::GetFinestGrid()
{
    if (!this->GetFineGrid())
    {
        return this;
    }
    else
    {
        return this->GetFineGrid()->GetFinestGrid();
    }
}

LIB_EXPORT void Grid::BackUpOldGrid()
{
    if (!oldGrid) oldGrid = new SimpleGrid();
    *oldGrid = *this;
}

void Grid::GridVerticeVelocity(RDouble *xt, RDouble *yt, RDouble *zt)
{
    RDouble *x = this->GetX();
    RDouble *y = this->GetY();
    RDouble *z = this->GetZ();

    SimpleGrid *oldGrid = this->GetOldGrid();

    RDouble *oldX = oldGrid->GetX();
    RDouble *oldY = oldGrid->GetY();
    RDouble *oldZ = oldGrid->GetZ();

    int numberOfNodes = this->GetNTotalNode();

    RDouble physicalTimeStep = GlobalDataBase::GetDoubleParaFromDB("physicalTimeStep");

    for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
    {
        xt[iNode] = (x[iNode] - oldX[iNode]) / physicalTimeStep;
        yt[iNode] = (y[iNode] - oldY[iNode]) / physicalTimeStep;
        zt[iNode] = (z[iNode] - oldZ[iNode]) / physicalTimeStep;
    }
}

LIB_EXPORT void Grid::RotTransVel(RDouble *xyzref, RDouble *dangle, RDouble *tcoord, RDouble **angle0, RDouble **coord0)
{
    //! This subroutine calculates the new set of coordinates after a rotation and a translation.

    RDouble psi, theta, phi, cp, sp, ct, st, cf, sf;
    RDouble xm, ym, zm, xp, yp, zp, dx, dy, dz;
    RDouble m11, m12, m13, m21, m22, m23, m31, m32, m33;

    int nTotalNode = this->GetNTotalNode();

    RDouble *x = this->GetX();
    RDouble *y = this->GetY();
    RDouble *z = this->GetZ();

    psi   = dangle[0];
    theta = dangle[1];
    phi   = dangle[2];

    dx = tcoord[0];
    dy = tcoord[1];
    dz = tcoord[2];
    
    xp = coord0[0][0] - dx;
    yp = coord0[1][0] - dy;
    zp = coord0[2][0] - dz;
    
    cp = cos(psi);
    sp = sin(psi);
    ct = cos(theta);
    st = sin(theta);
    cf = cos(phi);
    sf = sin(phi);

    m11 =   ct * cp;
    m12 = - ct * sp;
    m13 =   st;
    m21 =   cf * sp + sf * st * cp;
    m22 =   cf * cp - sf * st * sp;
    m23 = - sf * ct;
    m31 =   sf * sp - cf * st * cp;
    m32 =   sf * cp + cf * st * sp;
    m33 =   cf * ct;

    //! Update the grid.
    for (int inode = 0; inode < nTotalNode; ++ inode)
    {
        xm = x[inode] - xp;
        ym = y[inode] - yp;
        zm = z[inode] - zp;

        x[inode] = m11 * xm + m12 * ym + m13 * zm + xp + dx;
        y[inode] = m21 * xm + m22 * ym + m23 * zm + yp + dy;
        z[inode] = m31 * xm + m32 * ym + m33 * zm + zp + dz;
    }

    //! Update xyzref array.
    xm = xyzref[0] - xp;
    ym = xyzref[1] - yp;
    zm = xyzref[2] - zp;

    xyzref[0] = m11 * xm + m12 * ym + m13 * zm + xp + dx;
    xyzref[1] = m21 * xm + m22 * ym + m23 * zm + yp + dy;
    xyzref[2] = m31 * xm + m32 * ym + m33 * zm + zp + dz;
}

LIB_EXPORT void Grid::ComputeGDVel(RDouble *xt, RDouble *yt, RDouble *zt, RDouble **angle0, RDouble **coord0)
{
    //! This subroutine calculates the grid velocity at the cell vertices.
    
    //! angle0 = Euler's angles and their time derivatives.
    //! coord0 = position of the reference point abd their time derivatives.

    RDouble psip, theta, thetap, phi, phip, omega1, omega2, omega3;

    int nTotalNode = this->GetNTotalNode();

    RDouble *x = this->GetX();
    RDouble *y = this->GetY();
    RDouble *z = this->GetZ();

    //! Compute the components of the rotation vector.
    //! Omega (psip,thetap,phip) in the frame (k0,jh,i) where
    //! k0 = -sin theta * i + cos theta sin phi * j + cos theta cos phi * k
    //! jh = cos phi * j - sin phi * k
    //! Omega(p,q,r) in the frame (i,j,k) angle0(1,.) = (psi,theta,phi)
    //! psi is not set since it is not needed.

    psip   = angle0[0][1];
    theta  = angle0[1][0];
    thetap = angle0[1][1];
    phi    = angle0[2][0];
    phip   = angle0[2][1];

    omega1 = - psip * sin(theta)            + phip;
    omega2 =   psip * cos(theta) * sin(phi) + thetap * cos(phi);
    omega3 =   psip * cos(theta) * cos(phi) - thetap * sin(phi);

    //! Calculate the grid velocity with the formula dOM =dOP+Omega*PM.
    for (int inode = 0; inode < nTotalNode; ++ inode)
    {
        xt[inode] = coord0[0][1] + omega2 * (z[inode] - coord0[2][0]) - omega3 * (y[inode] - coord0[1][0]);
        yt[inode] = coord0[1][1] + omega3 * (x[inode] - coord0[0][0]) - omega1 * (z[inode] - coord0[2][0]);
        zt[inode] = coord0[2][1] + omega1 * (y[inode] - coord0[1][0]) - omega2 * (x[inode] - coord0[0][0]);
    }
}

LIB_EXPORT void Grid::TurnZAxis2YAxis()
{
    //! Turn the Y axis to be up, if the original Z axis is up.
    int nTotalNode = this->GetNTotalNode();

    RDouble *x, *y, *z;
    x = this->GetX();
    y = this->GetY();
    z = this->GetZ();

    RDouble tmp;
    for (int inode = 0; inode < nTotalNode; ++ inode)
    {
        tmp      = y[inode];
        y[inode] = z[inode];
        z[inode] = - tmp;
    }

    RDouble xmin, ymin, zmin, xmax, ymax, zmax;
    xmin = ymin = zmin =   LARGE;
    xmax = ymax = zmax = - LARGE;

    for (int inode = 0; inode < nTotalNode; ++ inode)
    {
        xmin = MIN(xmin, x[inode]);
        ymin = MIN(ymin, y[inode]);
        zmin = MIN(zmin, z[inode]);

        xmax = MAX(xmax, x[inode]);
        ymax = MAX(ymax, y[inode]);
        zmax = MAX(zmax, z[inode]);
    }

    ostringstream oss;
    oss << "After turning Z axis up coordinate to Y axis up coordinate\n";

    oss << setiosflags(ios::right);
    oss << setprecision(8);
    oss << setiosflags(ios::scientific);
    oss << setiosflags(ios::showpoint);
    int wordwidth = 16;
    oss << " xmin = " << setw(wordwidth) << xmin << " xmax = " << setw(wordwidth) << xmax << "\n";
    oss << " ymin = " << setw(wordwidth) << ymin << " ymax = " << setw(wordwidth) << ymax << "\n";
    oss << " zmin = " << setw(wordwidth) << zmin << " zmax = " << setw(wordwidth) << zmax << "\n";

    if (PHMPI::GetCurrentProcessorID() == PHMPI::GetServerProcessorID())
    {
        cout << oss.str();
    }
}

LIB_EXPORT void Grid::RegisterInterfaceField(const string &name, const int type, const int dimesion, const int solverID)
{
    if (!interfaceFields) return;

    interfaceFields->RegisterField(name, type, dimesion, solverID);
}

LIB_EXPORT void Grid::ReleaseInterfaceField()
{
    if (interfaceFields)
    {
        interfaceFields->RemoveAllField();
    }
}

LIB_EXPORT void Grid::RegisterInterpointField(const string &name, const int type, const int dimesion, const int solverID)
{
    if (!interpointFields)
    {
        return;
    }

    interpointFields->RegisterField(name, type, dimesion, solverID);
}

LIB_EXPORT void Grid::RegisterOversetField(const string &name, const int type, const int dimesion)
{
}

void Grid::FreeInterfaceField()
{
    if (!interfaceFields) return;
    if (interfaceFields->Size() == 0) return;

    interfaceFields->FreeInterfaceVar();
}

LIB_EXPORT void Grid::InitMovingGrids()
{
    int isAle = GlobalDataBase::GetIntParaFromDB("codeOfAleModel");
    if (isAle)
    {
        //ostringstream oss;
        //oss << "Warning: in Grid::InitMovingGrids, nstart parameter has been deleted from cfd_para" << " !!!" << endl;
        //cout << oss.str();

        //int nstart = GlobalDataBase::GetIntParaFromDB("nstart");
        int nstart = 0;

        if (nstart == 0 || nstart == 3)
        {
            //! nstart==0: start all over again.
            //! nstart==3: Unsteady-->Dynamics mesh.
            UpdateVolold();
        }
    }
}

//! Check grid validity.
LIB_EXPORT void Grid::GridValidityCheck()
{
    //! Interface check.
    int nIFace = this->GetNIFace();
    if (PHMPI::GetNumberOfProcessor() > 1 && nIFace == 0)
    {
        //! Parallel computing, but without any interface faces.
        ostringstream oss;
        oss << "Error: myid = " << PHMPI::GetCurrentProcessorID() 
            << ", zone index = " << this->GetZoneID() 
            << ", nIFace = 0 for parallel computing!!"
            << endl;

        PHMPI::FreeBlockData();
        cout << oss.str() << endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
}

void Grid::RotateTranslation()
{
}

void CheckGrid(Grid **grids, int nBlocks)
{
    for (int iZone = 0; iZone < nBlocks; ++ iZone)
    {
        StructGrid *grid = StructGridCast(grids[iZone]);
        InterfaceInfo *interfaceInfo = grid->GetInterfaceInfo();

        int nIFace = 0;
        if (interfaceInfo)
        {
            nIFace = interfaceInfo->GetNIFace();
        }

        if (nIFace == 0) continue;

        int *interFace2ZoneID       = interfaceInfo->GetInterFace2ZoneID();
        int *interFace2InterFaceID  = interfaceInfo->GetInterFace2InterFaceID();

        int target_zone, itface;
        int is, js, ks, nsurf, s_lr, it, jt, kt, t_lr;
        RDouble s_xfc, s_yfc, s_zfc, t_xfc, t_yfc, t_zfc;
        int il1, jl1, kl1;
  
        for (int iface = 0; iface < nIFace; ++ iface)
        {
            target_zone = interFace2ZoneID[iface];
            if (target_zone < iZone) continue;
            
            grid->GetSourceIndexIJK_Nsurf_LR(iface, 1, is, js, ks, nsurf, s_lr);
            if (s_lr == 1)
            {
                GetNsurfIndex(nsurf, il1, jl1, kl1);
                is = is + il1;
                js = js + jl1;
                ks = ks + kl1;
            }
            grid->FaceCoor(is, js, ks, nsurf, s_xfc, s_yfc, s_zfc);

            itface = interFace2InterFaceID[iface];

            if (target_zone == iZone)
            {
                grid->GetSourceIndexIJK_Nsurf_LR(itface, 1, it, jt, kt, nsurf, t_lr);
                if (t_lr == 1)
                {
                    GetNsurfIndex(nsurf, il1, jl1, kl1);
                    it = it + il1;
                    jt = jt + jl1;
                    kt = kt + kl1;
                }
                grid->FaceCoor(it, jt, kt, nsurf, t_xfc, t_yfc, t_zfc);
            }
            else
            {
                StructGrid *tgrid = StructGridCast(grids[target_zone]);
                tgrid->GetSourceIndexIJK_Nsurf_LR(itface, 1, it, jt, kt, nsurf, t_lr);
                if (t_lr == 1)
                {
                    GetNsurfIndex(nsurf, il1, jl1, kl1);
                    it = it + il1;
                    jt = jt + jl1;
                    kt = kt + kl1;
                }
                tgrid->FaceCoor(it, jt, kt, nsurf, t_xfc, t_yfc, t_zfc);
            }
            
            if ((fabs(s_xfc - t_xfc) > 1.0e-8) || (fabs(s_yfc - t_yfc) > 1.0e-8) || (fabs(s_zfc - t_zfc) > 1.0e-8))
            {
                PHMPI::FreeBlockData();
                cout << " Grid Link Error! " << endl;
                MPI_Abort(MPI_COMM_WORLD, -1);
            }
        }
    }
}

bool IsConvertGridToMixGrid()
{
    bool flag = false;
    if (GetTaskCode() == CREATE_GRID)
    {
        int gridobj = 0;
        GlobalDataBase::GetData("gridobj", &gridobj, PHINT, 1);

        if (gridobj == 1)
        {
            int gridtype;
            GlobalDataBase::GetData("gridtype", &gridtype, PHINT, 1);

            if (gridtype == MIXGRID)
            {
                flag = true;
            }
        }
    }
    return flag;
}

LIB_EXPORT void Create_Link_Info(VInt &pindex, int iZone, int iface, uint_t &fcount, set< DataStruct_Sort< VInt > > &face_list, VVInt &zoneid, VVInt &faceid, LinkStruct *link)
{
    sort(pindex.begin(), pindex.end());
    DataStruct_Sort< VInt > face(pindex, static_cast<int>(fcount));
    set< DataStruct_Sort< VInt > >::iterator iter = face_list.find(face);

    VVInt &facemap = link->GetFaceMap();

    if (iter == face_list.end())
    {
        face_list.insert(face);
        facemap[iZone][iface] = static_cast<int>(fcount);
        VInt l_zoneid;
        VInt l_faceid;
        l_zoneid.push_back(iZone);
        l_faceid.push_back(iface);

        zoneid.push_back(l_zoneid);
        faceid.push_back(l_faceid);

        ++ fcount;
    }
    else
    {
        facemap[iZone][iface] = iter->index;
        zoneid[iter->index].push_back(iZone);
        faceid[iter->index].push_back(iface);
    }
}

LIB_EXPORT void GetFaceCoorList(VInt &index, int nlist, RDouble *xlist, RDouble *ylist, RDouble *zlist, RDouble *x, RDouble *y, RDouble *z)
{
    for (int m = 0; m < nlist; ++ m)
    {
        xlist[m] = x[index[m]];
        ylist[m] = y[index[m]];
        zlist[m] = z[index[m]];
    }
}

LIB_EXPORT void GetFaceCoorList(int il, int jl, int kl, int ir, int jr, int kr, int ijk[4][3], int nlist,
                                RDouble *xlist, RDouble *ylist, RDouble *zlist, RDouble3D &x, RDouble3D &y, RDouble3D &z)
{
    for (int m = 0; m < nlist; ++ m)
    {
        int i = il + ijk[m][ir];
        int j = jl + ijk[m][jr];
        int k = kl + ijk[m][kr];

        xlist[m] = x(i, j, k);
        ylist[m] = y(i, j, k);
        zlist[m] = z(i, j, k);
    }
}

LIB_EXPORT void GetCoorIndexList(DataStruct_AdtTree<int, RDouble> *adt_tree, RDouble &diff_tolerance, int &pcount, RDouble *xlist, RDouble *ylist, RDouble *zlist, int nlist, vector<int> &pindex)
{
    RDouble coor[3], minwin[3], maxwin[3];

    for (int m = 0; m < nlist; ++ m)
    {
        coor[0] = xlist[m];
        coor[1] = ylist[m];
        coor[2] = zlist[m];

        minwin[0] = coor[0] - diff_tolerance;
        minwin[1] = coor[1] - diff_tolerance;
        minwin[2] = coor[2] - diff_tolerance;

        maxwin[0] = coor[0] + diff_tolerance;
        maxwin[1] = coor[1] + diff_tolerance;
        maxwin[2] = coor[2] + diff_tolerance;

        GetCoorIndex(adt_tree, coor, minwin, maxwin, pcount, pindex[m]);
    }
}

LIB_EXPORT int GetCoorIndex(DataStruct_AdtTree<int, RDouble> *adt_tree, RDouble *coor, RDouble *minwin, RDouble *maxwin, int &count, int &index)
{
    typedef DataStruct_AdtTree<int, RDouble> AdtTree;
    typedef DataStruct_AdtNode<int, RDouble> AdtNode;

    AdtNode *node;
    AdtTree::AdtNodeList node_list;

    node_list.resize(0);
    adt_tree->FindNodesInRegion(minwin, maxwin, node_list);

    if (node_list.size() == 0)
    {
        node = new AdtNode(3, coor, count);
        adt_tree->AddNode(node);
        index = count;
        ++ count;
        return 0;
    }
    else
    {
        if (node_list.size() > 1)
        {
            ostringstream oss;
            oss << " impossible node_list.size() = " << node_list.size() << "\n";
            TK_Exit::ExceptionExit(oss.str());
        }
        node  = node_list[0];
        index = node->GetData();
        return 1;
    }
}

LIB_EXPORT void ShiftMinMaxBox(RDouble *pmin, RDouble *pmax, RDouble tol)
{
    pmin[0] -= tol;
    pmin[1] -= tol;
    pmin[2] -= tol;

    pmax[0] += tol;
    pmax[1] += tol;
    pmax[2] += tol;
}

}
