#include "Pre_GridBase.h"
#include "Math_BasisFunction.h"
#include "IO_FileName.h"
#include "HyList.h"
#include "TK_Exit.h"
#include "Geo_UnstructBC.h"
#include "Constants.h"
#include "TK_Log.h"
#include "Pre_Reorder.h"


using namespace std;

namespace PHSPACE
{

CGNSRawCoor::CGNSRawCoor()
{
    type = new DataType_t[3];
    coor[0] = 0;
    coor[1] = 0;
    coor[2] = 0;
    nTotalNode = 0;
}

CGNSRawCoor::~CGNSRawCoor()
{
    DeAllocateData();
    delete [] type;
}

void CGNSRawCoor::AllocateData(int icoor, cgsize_t nTotalNode, DataType_t type)
{
    this->type[icoor] = type;
    this->nTotalNode = nTotalNode;
    if (type == RealSingle)
    {
        coor[icoor] = new RDouble[nTotalNode];
    }
    else
    {
        coor[icoor] = new RDouble[nTotalNode];
    }
}

void CGNSRawCoor::SetAllData(RDouble *x, RDouble *y, RDouble *z)
{
    RDouble *xyz[3];

    xyz[0] = x;
    xyz[1] = y;
    xyz[2] = z;

    for (int icoor = 0; icoor < 3; ++ icoor)
    {
        DataType_t type = this->type[icoor];
        SetData(icoor, type, xyz[icoor]);
    }

    //! Compute the max && min box.
    RDouble xMin, yMin, zMin, xMax, yMax, zMax;
    xMin = yMin = zMin =   LARGE;
    xMax = yMax = zMax = - LARGE;

    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        xMin = MIN(xMin, x[iNode]);
        yMin = MIN(yMin, y[iNode]);
        zMin = MIN(zMin, z[iNode]);

        xMax = MAX(xMax, x[iNode]);
        yMax = MAX(yMax, y[iNode]);
        zMax = MAX(zMax, z[iNode]);
    }

    ostringstream oss;
    oss << "Min && Max box of the CGNS grid: \n";
    oss << setiosflags(ios::right);
    oss << setprecision(8);
    oss << setiosflags(ios::scientific);
    oss << setiosflags(ios::showpoint);
    int wordwidth = 16;
    oss << "  xMin = " << setw(wordwidth) << xMin << ", xMax = " << setw(wordwidth) << xMax << "\n";
    oss << "  yMin = " << setw(wordwidth) << yMin << ", yMax = " << setw(wordwidth) << yMax << "\n";
    oss << "  zMin = " << setw(wordwidth) << zMin << ", zMax = " << setw(wordwidth) << zMax << "\n";
    oss << "\n";
    PrintToWindow(oss);
}

void CGNSRawCoor::SetData(int icoor, DataType_t type, RDouble *x)
{
    if (type == RealSingle)
    {
        RDouble *data = static_cast<RDouble *>(coor[icoor]);
        for (int iNode = 0; iNode < nTotalNode; ++ iNode)
        {
            x[iNode] = data[iNode];
        }
    }
    else
    {
        RDouble *data = static_cast<RDouble *>(coor[icoor]);
        for (int iNode = 0; iNode < nTotalNode; ++ iNode)
        {
            x[iNode] = data[iNode];
        }
    }
}

void CGNSRawCoor::DeAllocateData()
{
    for (int icoor = 0; icoor < 3; ++ icoor)
    {
        int type = this->type[icoor];
        if (type == RealSingle)
        {
            RDouble *data = static_cast<RDouble *>(coor[icoor]);
            delete [] data;
        }
        else
        {
            RDouble *data = static_cast<RDouble *>(coor[icoor]);
            delete [] data;
        }
    }
}

MultiIndex::MultiIndex()
{
    this->size = 0;
    this->id   = -1;
    data       = 0;
}

MultiIndex::MultiIndex(int size, int id)
{
    this->size = size;
    this->id   = id;
    data       = new int [size];
}

MultiIndex::MultiIndex(const MultiIndex &rhs)
{
    this->size = rhs.size;
    this->id   = rhs.id;
    this->data = new int [size];

    for (int i = 0; i < size; ++ i)
    {
        this->data[i] = rhs.data[i];
    }
}

MultiIndex & MultiIndex::operator = (const MultiIndex &rhs)
{
    if (this == &rhs) return *this;

    if (this->size != rhs.size)
    {
        delete [] this->data;
        this->size = rhs.size;
        this->data = new int [this->size];
    }

    this->id = rhs.id;
    for (int i = 0; i < size; ++ i)
    {
        this->data[i] = rhs.data[i];
    }

    return *this;
}

MultiIndex::~MultiIndex()
{
    delete [] data;
}

bool MultiIndex::operator < (const MultiIndex &rhs) const
{
    if (this->size != rhs.size) return this->size < rhs.size;
    for (int i = 0; i < this->size; ++ i)
    {
        if (this->data[i] != rhs.data[i])
        {
            return this->data[i] < rhs.data[i];
        }
    }
    return false;
}

void MultiIndex::SetData(int *data)
{
    for (int i = 0; i < size; ++ i)
    {
        this->data[i] = data[i];
    }
}

BaseElement::BaseElement()
{
    size         = NofValidElementTypes;
    elementface  = new CElemFace[size];
    elementpoint = new cgsize_t[size];
    InitElementPoint();
    InitElementFace();
};

BaseElement::~BaseElement()
{
    delete [] elementface;
    delete [] elementpoint;
};

void BaseElement::InitElementPoint()
{
    elementpoint[NODE ]    = 1;
    elementpoint[BAR_2]    = 2;
    elementpoint[TRI_3]    = 3;
    elementpoint[TRI_6]    = 6;
    elementpoint[QUAD_4]   = 4;
    elementpoint[QUAD_9]   = 9;
    elementpoint[TETRA_4]  = 4;
    elementpoint[TETRA_10] = 10;
    elementpoint[PYRA_5 ]  = 5;
    elementpoint[PYRA_14]  = 14;
    elementpoint[PENTA_6]  = 6;
    elementpoint[PENTA_18] = 18;
    elementpoint[HEXA_8 ]  = 8;
    elementpoint[HEXA_27]  = 27;
};

void BaseElement::InitElementFace()
{
    elementface[NODE].Init(NODE);
    elementface[BAR_2].Init(BAR_2);
    elementface[TRI_3].Init(TRI_3);
    elementface[QUAD_4].Init(QUAD_4);
    elementface[TETRA_4].Init(TETRA_4);
    elementface[TETRA_10].Init(TETRA_10);
    elementface[PYRA_5].Init(PYRA_5);
    elementface[PYRA_14].Init(PYRA_14);
    elementface[PENTA_6].Init(PENTA_6);
    elementface[PENTA_18].Init(PENTA_18);
    elementface[HEXA_8].Init(HEXA_8);
    elementface[HEXA_27].Init(HEXA_27);
};

CElemFace::CElemFace()
{
    face_num     = 0;
    face_clkwise = 0;
    facept_num   = 0;
    face_id      = 0;
}

CElemFace::~CElemFace()
{
    if (face_id)
    {
        for (int i = 0; i < face_num; ++ i)
        {
            delete [] face_id[i];
        }
    }
    delete [] face_id;
    delete [] face_clkwise;
    delete [] facept_num;
}

void CElemFace::Init(int elem_type)
{
    if (elem_type == NODE)
    {
        face_num = 1;
        face_id    = new int * [face_num];
        facept_num = new int [face_num];

        facept_num[0] = 1;

        for (int i = 0; i < face_num; ++ i)
        {
            face_id[i] = new int [facept_num[i]];
        }

        face_id[0][0] = 0;
    }
    else if (elem_type == BAR_2)
    {
        face_num = 1;
        face_id    = new int * [face_num];
        facept_num = new int [face_num];

        facept_num[0] = 2;

        for (int i = 0; i < face_num; ++ i)
        {
            face_id[i] = new int [facept_num[i]];
        }

        face_id[0][0] = 0;
        face_id[0][1] = 1;
    }
    else if (elem_type == TRI_3)
    {
        face_num = 3;
        face_id    = new int * [face_num];
        facept_num = new int [face_num];

        facept_num[0] = 2;
        facept_num[1] = 2;
        facept_num[2] = 2;

        for (int i = 0; i < face_num; ++ i)
        {
            face_id[i] = new int [facept_num[i]];
        }

        face_id[0][0] = 0;
        face_id[0][1] = 1;

        face_id[1][0] = 1;
        face_id[1][1] = 2;

        face_id[2][0] = 2;
        face_id[2][1] = 0;
    }
    else if (elem_type == QUAD_4)
    {
        face_num = 4;
        face_id    = new int * [face_num];
        facept_num = new int [face_num];

        facept_num[0] = 2;
        facept_num[1] = 2;
        facept_num[2] = 2;
        facept_num[3] = 2;

        for (int i = 0; i < face_num; ++ i)
        {
            face_id[i] = new int [facept_num[i]];
        }

        face_id[0][0] = 0;
        face_id[0][1] = 1;

        face_id[1][0] = 1;
        face_id[1][1] = 2;

        face_id[2][0] = 2;
        face_id[2][1] = 3;

        face_id[3][0] = 3;
        face_id[3][1] = 0;
    }
    else if (elem_type == TETRA_4)
    {
        face_num = 4;
        face_id    = new int * [face_num];
        facept_num = new int [face_num];

        facept_num[0] = 3;
        facept_num[1] = 3;
        facept_num[2] = 3;
        facept_num[3] = 3;

        for (int i = 0; i < face_num; ++ i)
        {
            face_id[i] = new int [facept_num[i]];
        }

        face_id[0][0] = 0;
        face_id[0][1] = 2;
        face_id[0][2] = 1;

        face_id[1][0] = 0;
        face_id[1][1] = 1;
        face_id[1][2] = 3;

        face_id[2][0] = 1;
        face_id[2][1] = 2;
        face_id[2][2] = 3;

        face_id[3][0] = 2;
        face_id[3][1] = 0;
        face_id[3][2] = 3;
    }
    else if (elem_type == TETRA_10)
    {
        face_num = 4;
        face_id = new int * [face_num];
        facept_num = new int[face_num];

        facept_num[0] = 6;
        facept_num[1] = 6;
        facept_num[2] = 6;
        facept_num[3] = 6;

        for (int i = 0; i < face_num; ++ i)
        {
            face_id[i] = new int[facept_num[i]];
        }

        face_id[0][0] = 0;
        face_id[0][1] = 1;
        face_id[0][2] = 2;
        face_id[0][3] = 4;
        face_id[0][4] = 5;
        face_id[0][5] = 6;

        face_id[1][0] = 0;
        face_id[1][1] = 1;
        face_id[1][2] = 3;
        face_id[1][3] = 4;
        face_id[1][4] = 8;
        face_id[1][5] = 7;

        face_id[2][0] = 1;
        face_id[2][1] = 2;
        face_id[2][2] = 3;
        face_id[2][3] = 5;
        face_id[2][4] = 9;
        face_id[2][5] = 8;

        face_id[3][0] = 2;
        face_id[3][1] = 3;
        face_id[3][2] = 0;
        face_id[3][3] = 9;
        face_id[3][4] = 7;
        face_id[3][5] = 6;
    }
    else if (elem_type == PYRA_5)
    {
        face_num = 5;
        face_id    = new int * [face_num];
        facept_num = new int [face_num];

        facept_num[0] = 4;
        facept_num[1] = 3;
        facept_num[2] = 3;
        facept_num[3] = 3;
        facept_num[4] = 3;

        for (int i = 0; i < face_num; ++ i)
        {
            face_id[i] = new int [facept_num[i]];
        }

        face_id[0][0] = 0;
        face_id[0][1] = 3;
        face_id[0][2] = 2;
        face_id[0][3] = 1;

        face_id[1][0] = 0;
        face_id[1][1] = 1;
        face_id[1][2] = 4;

        face_id[2][0] = 1;
        face_id[2][1] = 2;
        face_id[2][2] = 4;

        face_id[3][0] = 2;
        face_id[3][1] = 3;
        face_id[3][2] = 4;

        face_id[4][0] = 3;
        face_id[4][1] = 0;
        face_id[4][2] = 4;
    }
    else if (elem_type == PYRA_14)
    {
        face_num = 5;
        face_id = new int * [face_num];
        facept_num = new int[face_num];

        facept_num[0] = 9;
        facept_num[1] = 6;
        facept_num[2] = 6;
        facept_num[3] = 6;
        facept_num[4] = 6;

        for (int i = 0; i < face_num; ++ i)
        {
            face_id[i] = new int[facept_num[i]];
        }

        face_id[0][0] = 0;
        face_id[0][1] = 3;
        face_id[0][2] = 2;
        face_id[0][3] = 1;
        face_id[0][4] = 8;
        face_id[0][5] = 7;
        face_id[0][6] = 6;
        face_id[0][7] = 5;
        face_id[0][8] = 13;

        face_id[1][0] = 0;
        face_id[1][1] = 1;
        face_id[1][2] = 4;
        face_id[1][3] = 5;
        face_id[1][4] = 10;
        face_id[1][5] = 9;

        face_id[2][0] = 1;
        face_id[2][1] = 2;
        face_id[2][2] = 4;
        face_id[2][3] = 6;
        face_id[2][4] = 11;
        face_id[2][5] = 10;

        face_id[3][0] = 2;
        face_id[3][1] = 3;
        face_id[3][2] = 4;
        face_id[3][3] = 7;
        face_id[3][4] = 12;
        face_id[3][5] = 11;

        face_id[4][0] = 3;
        face_id[4][1] = 0;
        face_id[4][2] = 4;
        face_id[4][3] = 8;
        face_id[4][4] = 9;
        face_id[4][5] = 12;
    }
    else if (elem_type == PENTA_6)
    {
        face_num = 5;
        face_id    = new int * [face_num];
        facept_num = new int [face_num];

        facept_num[0] = 3;
        facept_num[1] = 3;
        facept_num[2] = 4;
        facept_num[3] = 4;
        facept_num[4] = 4;

        for (int i = 0; i < face_num; ++ i)
        {
            face_id[i] = new int [facept_num[i]];
        }

        face_id[0][0] = 0;
        face_id[0][1] = 2;
        face_id[0][2] = 1;

        face_id[1][0] = 3;
        face_id[1][1] = 4;
        face_id[1][2] = 5;

        face_id[2][0] = 0;
        face_id[2][1] = 1;
        face_id[2][2] = 4;
        face_id[2][3] = 3;

        face_id[3][0] = 1;
        face_id[3][1] = 2;
        face_id[3][2] = 5;
        face_id[3][3] = 4;

        face_id[4][0] = 2;
        face_id[4][1] = 0;
        face_id[4][2] = 3;
        face_id[4][3] = 5;
    }
    else if (elem_type == PENTA_18)
    {
        face_num = 5;
        face_id = new int * [face_num];
        facept_num = new int[face_num];

        facept_num[0] = 9;
        facept_num[1] = 9;
        facept_num[2] = 9;
        facept_num[3] = 6;
        facept_num[4] = 6;

        for (int i = 0; i < face_num; ++ i)
        {
            face_id[i] = new int[facept_num[i]];
        }

        face_id[0][0] = 0;
        face_id[0][1] = 1;
        face_id[0][2] = 4;
        face_id[0][3] = 3;
        face_id[0][4] = 6;
        face_id[0][5] = 10;
        face_id[0][6] = 12;
        face_id[0][7] = 9;
        face_id[0][8] = 15;

        face_id[1][0] = 1;
        face_id[1][1] = 2;
        face_id[1][2] = 5;
        face_id[1][3] = 4;
        face_id[1][4] = 7;
        face_id[1][5] = 11;
        face_id[1][6] = 13;
        face_id[1][7] = 10;
        face_id[1][8] = 16;

        face_id[2][0] = 2;
        face_id[2][1] = 0;
        face_id[2][2] = 3;
        face_id[2][3] = 5;
        face_id[2][4] = 8;
        face_id[2][5] = 9;
        face_id[2][6] = 14;
        face_id[2][7] = 11;
        face_id[2][8] = 17;

        face_id[3][0] = 0;
        face_id[3][1] = 2;
        face_id[3][2] = 1;
        face_id[3][3] = 8;
        face_id[3][4] = 7;
        face_id[3][5] = 6;

        face_id[4][0] = 3;
        face_id[4][1] = 4;
        face_id[4][2] = 5;
        face_id[4][3] = 12;
        face_id[4][4] = 13;
        face_id[4][5] = 14;
    }
    else if (elem_type == HEXA_8)
    {
        face_num = 6;
        face_id    = new int * [face_num];
        facept_num = new int [face_num];

        facept_num[0] = 4;
        facept_num[1] = 4;
        facept_num[2] = 4;
        facept_num[3] = 4;
        facept_num[4] = 4;
        facept_num[5] = 4;

        for (int i = 0; i < face_num; ++ i)
        {
            face_id[i] = new int [facept_num[i]];
        }

        face_id[0][0] = 0;
        face_id[0][1] = 3;
        face_id[0][2] = 2;
        face_id[0][3] = 1;

        face_id[1][0] = 4;
        face_id[1][1] = 5;
        face_id[1][2] = 6;
        face_id[1][3] = 7;

        face_id[2][0] = 3;
        face_id[2][1] = 0;
        face_id[2][2] = 4;
        face_id[2][3] = 7;

        face_id[3][0] = 0;
        face_id[3][1] = 1;
        face_id[3][2] = 5;
        face_id[3][3] = 4;

        face_id[4][0] = 1;
        face_id[4][1] = 2;
        face_id[4][2] = 6;
        face_id[4][3] = 5;

        face_id[5][0] = 2;
        face_id[5][1] = 3;
        face_id[5][2] = 7;
        face_id[5][3] = 6;
    }
    else if (elem_type == HEXA_27)
    {
        face_num = 6;
        face_id = new int * [face_num];
        facept_num = new int[face_num];

        facept_num[0] = 9;
        facept_num[1] = 9;
        facept_num[2] = 9;
        facept_num[3] = 9;
        facept_num[4] = 9;
        facept_num[5] = 9;

        for (int i = 0; i < face_num; ++ i)
        {
            face_id[i] = new int[facept_num[i]];
        }

        face_id[0][0] = 0;
        face_id[0][1] = 3;
        face_id[0][2] = 2;
        face_id[0][3] = 1;
        face_id[0][4] = 11;
        face_id[0][5] = 10;
        face_id[0][6] = 9;
        face_id[0][7] = 8;
        face_id[0][8] = 20;

        face_id[1][0] = 0;
        face_id[1][1] = 1;
        face_id[1][2] = 5;
        face_id[1][3] = 4;
        face_id[1][4] = 8;
        face_id[1][5] = 13;
        face_id[1][6] = 16;
        face_id[1][7] = 12;
        face_id[1][8] = 21;

        face_id[2][0] = 1;
        face_id[2][1] = 2;
        face_id[2][2] = 6;
        face_id[2][3] = 5;
        face_id[2][4] = 9;
        face_id[2][5] = 14;
        face_id[2][6] = 17;
        face_id[2][7] = 13;
        face_id[2][8] = 22;

        face_id[3][0] = 2;
        face_id[3][1] = 3;
        face_id[3][2] = 7;
        face_id[3][3] = 6;
        face_id[3][4] = 10;
        face_id[3][5] = 15;
        face_id[3][6] = 18;
        face_id[3][7] = 14;
        face_id[3][8] = 23;

        face_id[4][0] = 0;
        face_id[4][1] = 4;
        face_id[4][2] = 7;
        face_id[4][3] = 3;
        face_id[4][4] = 12;
        face_id[4][5] = 19;
        face_id[4][6] = 15;
        face_id[4][7] = 11;
        face_id[4][8] = 24;

        face_id[5][0] = 4;
        face_id[5][1] = 5;
        face_id[5][2] = 6;
        face_id[5][3] = 7;
        face_id[5][4] = 16;
        face_id[5][5] = 17;
        face_id[5][6] = 18;
        face_id[5][7] = 19;
        face_id[5][8] = 25;
    }
    else
    {
        TK_Exit::ExceptionExit(" CElemFace::Init, Unknown type\n");
    }
}

BaseBCType::BaseBCType()
{
    ibcmax = NofValidBCTypes;
    Init();
};

void BaseBCType::Init()
{
    using namespace PHSPACE;

    bctypemap.insert(pair<int, int>(BCTypeUserDefined, PHENGLEI::OVERSET             ));
    bctypemap.insert(pair<int, int>(BCSymmetryPlane  , PHENGLEI::SYMMETRY            ));
    bctypemap.insert(pair<int, int>(BCInflow         , PHENGLEI::INFLOW              ));
    bctypemap.insert(pair<int, int>(BCOutflow        , PHENGLEI::OUTFLOW             ));
    bctypemap.insert(pair<int, int>(BCTunnelOutflow  , PHENGLEI::OUTFLOW             ));
    bctypemap.insert(pair<int, int>(BCWall           , PHENGLEI::SOLID_SURFACE       ));
    bctypemap.insert(pair<int, int>(BCWallInviscid   , PHENGLEI::SOLID_SURFACE       ));
    bctypemap.insert(pair<int, int>(BCWallViscous    , PHENGLEI::SOLID_SURFACE       ));
    bctypemap.insert(pair<int, int>(BCDegenerateLine , PHENGLEI::POLE                ));
    bctypemap.insert(pair<int, int>(BCTypeNull       , PHENGLEI::NO_BOUNDARY_CONDITION));
    bctypemap.insert(pair<int, int>(BCFarfield       , PHENGLEI::FARFIELD            ));
    bctypemap.insert(pair<int, int>(BCGeneral        , PHENGLEI::GENERIC_1           ));

    fts2cgns.insert(pair<int, int>(PHENGLEI::OVERSET              , BCTypeUserDefined));
    fts2cgns.insert(pair<int, int>(PHENGLEI::SYMMETRY             , BCSymmetryPlane ));
    fts2cgns.insert(pair<int, int>(PHENGLEI::INFLOW               , BCInflow        ));
    fts2cgns.insert(pair<int, int>(PHENGLEI::OUTFLOW              , BCOutflow       ));
    fts2cgns.insert(pair<int, int>(PHENGLEI::SOLID_SURFACE        , BCWall          ));
    fts2cgns.insert(pair<int, int>(PHENGLEI::POLE                 , BCDegenerateLine));
    fts2cgns.insert(pair<int, int>(PHENGLEI::NO_BOUNDARY_CONDITION, BCTypeNull      ));
    fts2cgns.insert(pair<int, int>(PHENGLEI::FARFIELD             , BCFarfield      ));

    int maxNumberOfBody = 20;
    for (int iBody = 0; iBody < maxNumberOfBody; iBody ++)
    {
        bctypemap.insert(pair<int, int>(NofValidBCTypes + iBody, PHENGLEI::StartNumberOfBodyOfHyperFLOW + iBody));
        fts2cgns. insert(pair<int, int>(PHENGLEI::StartNumberOfBodyOfHyperFLOW + iBody, NofValidBCTypes + iBody));
    }
}

BaseBCType::~BaseBCType()
{
}

RawGrid::RawGrid()
{
    x          = 0;
    y          = 0;
    z          = 0;
    nTotalNode = 0;
}

RawGrid::~RawGrid()
{
    //delete [] x;
    //delete [] y;
    //delete [] z;
}

Base_Grid_Conn::Base_Grid_Conn()
{
    left_cell_of_face        = new vector<cgsize_t>;
    right_cell_of_face       = new vector<cgsize_t>;
    face_location            = new vector<int>;
    bcFace                   = new vector<CGNS_BCFace *>;
    node_number_of_each_face = new vector<int>;
    face2node                = new vector< vector<cgsize_t> >;
    faceset                  = new set<MultiIndex>;
    node_number_of_each_cell = new vector<int>;
    cell2node                = new vector< vector<cgsize_t> >;
    nTotalNode               = 0;
    nTotalCell               = 0;
    nBoundFace               = 0;
    nTotalFace               = 0;
}

Base_Grid_Conn::~Base_Grid_Conn()
{
    delete left_cell_of_face;
    delete right_cell_of_face;
    delete face_location;
    delete node_number_of_each_face;
    delete face2node ;
    delete faceset;
    delete node_number_of_each_cell;
    delete cell2node;

    vector<CGNS_BCFace *>::iterator iter;
    for (iter = bcFace->begin(); iter != bcFace->end(); ++ iter)
    {
        delete *iter;
    }
    delete bcFace;
}

void Base_Grid_Conn::FreeFaceSet()
{
    delete faceset;
    faceset = 0;
}

void Base_Grid_Conn::FreeFace2Node()
{
    delete face2node;
    face2node = 0;
}

void Base_Grid_Conn::FreeCell2Node()
{
    delete cell2node;
    cell2node = 0;
}

void Base_Grid_Conn::FreeCellOfFace()
{
    delete left_cell_of_face;
    left_cell_of_face = 0;

    delete right_cell_of_face;
    right_cell_of_face = 0;
}

void Base_Grid_Conn::FreeNodeNumberOfEachFace()
{
    delete node_number_of_each_face;
    node_number_of_each_face = 0;
}

void Base_Grid_Conn::FreeNodeNumberOfEachCell()
{
    delete node_number_of_each_cell;
    node_number_of_each_cell = 0;
}

CGNSBase::CGNSBase()
{
    elem_type        = 0;
    istart           = 0;
    iend             = 0;
    conn_list        = 0;
    nBCElem          = 0;
    bc_type          = 0;
    bc_conn_list     = 0;
    boundaryName     = 0;
    ElementDataSize  = 0;
    mixed_conn_index = 0;
}

CGNSBase::~CGNSBase()
{
    delete [] boundaryName;
    delete [] elem_type;
    delete [] istart;
    delete [] iend;
    delete [] nBCElem;
    delete [] bc_type;
    delete [] ElementDataSize;

    if (conn_list)
    {
        for (int i = 0; i < nSections; ++ i)
        {
            delete [] conn_list[i];
        }
    }
    delete [] conn_list;

    if (mixed_conn_index)
    {
        for (int i = 0; i < nSections; ++ i)
        {
            delete [] mixed_conn_index[i];
        }
    }
    delete [] mixed_conn_index;

    if (bc_conn_list)
    {
        for (int iBCRegion = 0; iBCRegion < nBCRegions; ++ iBCRegion)
        {
            delete [] bc_conn_list[iBCRegion];
        }
    }
    delete [] bc_conn_list;
}

void CGNSBase::CreateElement(int nSections)
{
    this->nSections = nSections;

    conn_list = new cgsize_t * [nSections];
    istart    = new cgsize_t [nSections];
    iend      = new cgsize_t [nSections];
    elem_type = new int [nSections];
    ElementDataSize = new cgsize_t [nSections];
    mixed_conn_index= new cgsize_t * [nSections];
    cout << " nSections = " << nSections << "\n";
}

void CGNSBase::CreateBC(int nBCRegions)
{
    this->nBCRegions = nBCRegions;
    bc_conn_list     = new cgsize_t * [nBCRegions];
    nBCElem          = new int [nBCRegions];
    bc_type          = new int [nBCRegions];
    bc_elem_set_type = new int [nBCRegions];
    bc_grid_location = new int [nBCRegions];
    boundaryName     = new string[nBCRegions];
}

void CGNSBase::PreProcess(Base_Grid_Conn *gconn, BaseElement *base_elem)
{
    cout << "PreProcess ...\n";
    gconn->SetNTotalCell(nTotalCell);
    gconn->SetNTotalNode(nTotalNode);
    gconn->SetNBoundFace(nBoundFace);

    for (int i = 0; i < nSections; ++ i)
    {
        istart[i] -= 1;
        iend[i]   -= 1;
    }

    cgsize_t *base_elem_elem_pt = base_elem->GetElementPoint();

    cgsize_t size;
    for (int i = 0; i < nSections; ++ i)
    {
        if (elem_type[i] == MIXED)
        {
            size = ElementDataSize[i];
        }
        else
        {
            cgsize_t nelem = iend[i] - istart[i] + 1;
            cgsize_t pt    = base_elem_elem_pt[ elem_type[i] ];
            size  = nelem * pt;
        }

        for (cgsize_t j = 0; j < size; ++ j)
        {
            conn_list[i][j] -= 1;
        }
    }

    for (int iBCRegion = 0; iBCRegion < nBCRegions; ++ iBCRegion)
    {
        //! It is true for all situation.
        for (int j = 0; j < nBCElem[iBCRegion]; ++ j)
        {
            bc_conn_list[iBCRegion][j] -= 1;
        }
    }

    //! For boundary BC's, it should be checked index base.
    int ibase = 1;

    for (int iBCRegion = 0; iBCRegion < nBCRegions; ++ iBCRegion)
    {
        for (int j = 0; j < nBCElem[iBCRegion]; ++ j)
        {
            if (bc_conn_list[iBCRegion][j] < nTotalCell)
            {
                ibase = 0;
                break;
            }
        }
        if (ibase == 0) break;
    }

    if (ibase == 0)
    {
        for (int iBCRegion = 0; iBCRegion < nBCRegions; ++ iBCRegion)
        {
            if ((this->bc_grid_location[iBCRegion] != Vertex))
            {
                for (int j = 0; j < nBCElem[iBCRegion]; ++ j)
                {
                    bc_conn_list[iBCRegion][j] += nTotalCell;
                }
            }
        }
    }
}

void CGNSBase::GetElem(cgsize_t *elem_source, cgsize_t iCell, int ptnum, cgsize_t *&elem)
{
    elem = &elem_source[ptnum * iCell];
}

void CGNSBase::GetFace(Base_Grid_Conn *gconn, cgsize_t *elem, cgsize_t iCell, CElemFace *elemface)
{
    int face_list[100], face_now[100];
    int maxface = elemface->GetFaceNumber();

    vector<cgsize_t> single_face(100);
    single_face.resize(0);

    vector< vector<cgsize_t> > *face2node = gconn->GetFace2Node();
    vector<cgsize_t> *left_cell_of_face   = gconn->GetLeftCellOfFace();
    vector<cgsize_t> *right_cell_of_face  = gconn->GetRightCellOfFace();
    vector<int> *node_number_of_each_face = gconn->GetNodeNumberOfEachFace();
    vector<CGNS_BCFace *> *bcFace         = gconn->GetBCFace();
    vector<int> *face_location            = gconn->GetFaceLocation();
    set<MultiIndex> *faceSet              = gconn->GetFaceSet();

    int *elemface_facept_num = elemface->GetFacePointNumber();
    int **elemface_face_id   = elemface->GetFaceID();

    for (int iFace = 0; iFace < maxface; ++ iFace)
    {
        int facept_num = elemface_facept_num[iFace];
        int *face_id   = elemface_face_id[iFace];

        MultiIndex faceIndex(facept_num, static_cast<int>(left_cell_of_face->size()));
        for (int j = 0; j < facept_num; ++ j)
        {
            face_list[j] = static_cast<int>(elem[face_id[j]]);
            face_now[j]  = face_list[j];
        }
        sort(face_list, face_list + facept_num);
        faceIndex.SetData(face_list);

        set<MultiIndex>::iterator iter = faceSet->find(faceIndex);

        if (iter == faceSet->end())
        {
            faceSet->insert(faceIndex);
            left_cell_of_face->push_back(iCell);
            right_cell_of_face->push_back(-1);

            CGNS_BCFace *boundaryFace = new CGNS_BCFace(static_cast<int>(bcFace->size()), 0, "");
            bcFace->push_back(boundaryFace);

            face_location->push_back(iFace);

            single_face.resize(facept_num);
            for (int j = 0; j < facept_num; ++ j)
            {
                single_face[j] = face_now[j];
            }

            face2node->push_back(single_face);
            node_number_of_each_face->push_back(facept_num);
        }
        else
        {
            (*right_cell_of_face)[iter->id] = iCell;
        }
    }
}

void CGNSBase::PostProcess(Base_Grid_Conn *gconn)
{
    cout << "  Boundary Faces Reordering ..." << endl;
    gconn->FreeFaceSet();

    vector<cgsize_t> *left_cell_of_face  = gconn->GetLeftCellOfFace();
    vector<cgsize_t> *right_cell_of_face = gconn->GetRightCellOfFace();

    cgsize_t nTotalFace = right_cell_of_face->size();

    gconn->SetNTotalFace(nTotalFace);

    nBoundFace = 0;
    for (cgsize_t iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        if ((*right_cell_of_face)[iFace] == -1)
        {
            ++ nBoundFace;
        }
    }

    gconn->SetNBoundFace(nBoundFace);

    vector< vector<cgsize_t> > *face2node = gconn->GetFace2Node();
    vector<int> *node_number_of_each_face = gconn->GetNodeNumberOfEachFace();
    vector<CGNS_BCFace *> *bcFace = gconn->GetBCFace();

    //! Reorder faces, make sure the boundary faces lie in the first part.
    cgsize_t jstart = nBoundFace;
    for (cgsize_t iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        if ((*right_cell_of_face)[iFace] != -1)
        {
            //! It needs to swap.
            for (cgsize_t jFace = jstart; jFace < nTotalFace; ++ jFace)
            {
                if ((*right_cell_of_face)[jFace] == -1)
                {
                    //! Swap i,j.
                    SWAP((*left_cell_of_face)[iFace], (*left_cell_of_face)[jFace]);
                    SWAP((*right_cell_of_face)[iFace], (*right_cell_of_face)[jFace]);

                    //! Swap bctype && bcname only, do not swap faceID.
                    //! So, did not use SWAP((*bcFace)[iFace], (*bcFace)[jFace]) here.
                    ASSERT((*bcFace)[iFace]->GetFaceID() == iFace);
                    int iFaceBCType = (*bcFace)[iFace]->GetBCType();
                    int jFaceBCType = (*bcFace)[jFace]->GetBCType();
                    SWAP(iFaceBCType, jFaceBCType);
                    (*bcFace)[iFace]->SetBCType(iFaceBCType);
                    (*bcFace)[jFace]->SetBCType(jFaceBCType);
                    string iFaceBCName = (*bcFace)[iFace]->GetBCName();
                    string jFaceBCName = (*bcFace)[jFace]->GetBCName();
                    SWAP(iFaceBCName, jFaceBCName);
                    (*bcFace)[iFace]->SetBCName(iFaceBCName);
                    (*bcFace)[jFace]->SetBCName(jFaceBCName);    

                    SWAP((*face2node)[iFace], (*face2node)[jFace]);
                    SWAP((*node_number_of_each_face)[iFace], (*node_number_of_each_face)[jFace]);

                    jstart = jFace + 1;
                    break;
                }
            }
        }
    }
}

void CGNSBase::GetISection(cgsize_t iCell, int &isection)
{
    cgsize_t iabscell = iCell;
    isection = 0;
    for (int i = 0; i < nSections; ++ i)
    {
        cgsize_t ist = istart[i];
        cgsize_t ied = iend  [i];
        if (ist <= iabscell && iabscell <= ied)
        {
            isection = i;
            break;
        }
    }
}

void CGNSBase::GetBCFace(BaseElement *base_elem, cgsize_t iCell, int &facept_num, cgsize_t *&bcface)
{
    cgsize_t iabscell = iCell;
    int isection = 0;

    GetISection(iCell, isection);
    
    ElementType_t etype = static_cast<ElementType_t>(elem_type[isection]);
    cgsize_t *base_elem_elem_pt = base_elem->GetElementPoint();

    if (etype != MIXED)
    {
        facept_num = base_elem_elem_pt[etype];
        GetElem(conn_list[isection], iabscell - istart[isection], facept_num, bcface);
    }
    else
    {
        int iNode;
        iNode = static_cast<int>(mixed_conn_index[isection][iabscell - istart[isection]]);
        etype = static_cast<ElementType_t>(conn_list[isection][iNode] + 1);
        facept_num = base_elem_elem_pt[etype];
        GetBCFaceMixed(conn_list[isection], iabscell - istart[isection], isection, bcface);
    }
}

void CGNSBase::GetBCFaceMixed(cgsize_t *elem_source, cgsize_t iCell, int isection, cgsize_t *&elem)
{
    elem = &elem_source[mixed_conn_index[isection][iCell] + 1];
}

void CGNSBase::GetElem(BaseElement *base_elem, cgsize_t iCell, int &etype, cgsize_t *&elem)
{
    int isection = 0;
    int elempt_num;

    GetISection(iCell, isection);
    etype = elem_type[isection];

    cgsize_t *base_elem_elem_pt = base_elem->GetElementPoint();

    if (etype == MIXED)
    {
        etype = static_cast<ElementType_t>(conn_list[isection][mixed_conn_index[isection][iCell]] + 1);
        elem = &conn_list[isection][mixed_conn_index[isection][iCell] + 1];
    }
    else
    {
        elempt_num = base_elem_elem_pt[etype];
        GetElem(conn_list[isection], iCell - istart[isection], elempt_num, elem);
    }
}

void CGNSBase::ProcessBCFace(Base_Grid_Conn *gconn, BaseElement *base_elem, 
                             cgsize_t iElement, BCType_t bctype, const string &BCName)
{
    cgsize_t *bcface = 0;
    int facept_num;

    GetBCFace(base_elem, iElement, facept_num, bcface);
    
    if (!MarkBCFace(gconn, bcface, facept_num, bctype, BCName))
    {
        ostringstream oss;
        oss << " ERROR in MarkBCFace : iElement = " << iElement << "\n";
        TK_Exit::ExceptionExit(oss.str());
    }
}

bool CGNSBase::MarkBCFace(Base_Grid_Conn *gconn, cgsize_t *bcface, int facept_num, BCType_t bctype, const string &BCName)
{
    int face_list[100];
    MultiIndex idx(facept_num, 0);

    for (int i = 0; i < facept_num; ++ i)
    {
        face_list[i] = static_cast<int>(bcface[i]);
    }

    sort(face_list, face_list + facept_num);

    idx.SetData(face_list);

    set<MultiIndex> *faceset = gconn->GetFaceSet();

    set<MultiIndex>::iterator iter = faceset->find(idx);
    vector<CGNS_BCFace *> *bcFace = gconn->GetBCFace();

    if (iter == faceset->end())
    {
        cout << "impossible, boundary face should be in the total face\n";
        for (int i = 0; i < facept_num; ++ i)
        {
            cout << bcface[ i ] << " ";
        }
        cout << "\n";
        return false;
    }
    else
    {
        (*bcFace)[iter->id]->SetBCType(bctype);
        (*bcFace)[iter->id]->SetBCName(BCName);
        return true;
    }
}

void CGNSBase::ComputeGridConn(Base_Grid_Conn *gconn, BaseElement *base_elem)
{
    PreProcess(gconn, base_elem);
    ComputeElemFace(gconn, base_elem);
    ComputeBCFace(gconn, base_elem);
    ComputeCellToNode(gconn, base_elem);
    PostProcess(gconn);
	
	int gridReorder = PHSPACE::GlobalDataBase::GetIntParaFromDB("gridReorder");
	
	if (gridReorder == 1)
	{
		using namespace Reorder;
		//compute bandwidth
		printf("Before ReOrder...\n");
		CallBandwidthEvaluate(gconn);
		CallMetricLeftRightCellOfFace(gconn);
        //! CELLREORDER
		CallReorderCellLabel(gconn);
		//! FACEREORDER
		CallReorderFaceLabel(gconn);
		//! metrics for left_cell_of_face and right_cell_of_face
		printf("After ReOrder...\n");
		CallBandwidthEvaluate(gconn);
		CallMetricLeftRightCellOfFace(gconn);
	}
}

void CGNSBase::ComputeCellToNode(Base_Grid_Conn *gconn, BaseElement *base_elem)
{
    int nSections = this->GetNumberOfSections();
    int *base_cgns_elem_type = this->GetElementType();
    int nTotalCell = static_cast<int>(this->GetNTotalCell());
    cgsize_t **conn_list = this->GetElementConnectionList();

    vector<int> *node_number_of_each_cell = gconn->GetNodeNumberOfEachCell();
    vector< vector<cgsize_t> > *cell2node = gconn->GetCell2Node();
    //cell2node->reserve(nTotalCell);
    cell2node->resize(nTotalCell);
    //node_number_of_each_cell->reserve(nTotalCell);
    node_number_of_each_cell->resize(nTotalCell);

    cgsize_t *base_elem_elem_pt = base_elem->GetElementPoint();

    bool flag = false;
    for (int sec = 0; sec < nSections; ++ sec)
    {
        cgsize_t *base_cgns_istart = this->GetIndexOfStart();
        cgsize_t *base_cgns_iend   = this->GetIndexOfEnd();
        ElementType_t etype = static_cast<ElementType_t>(base_cgns_elem_type[sec]);

        if (etype != MIXED)
        {
            int nPoints = base_elem_elem_pt[etype];

            int count = 0;
            for (cgsize_t iCell = base_cgns_istart[sec]; iCell <= base_cgns_iend[sec]; ++ iCell)
            {
                if (iCell >= nTotalCell)
                {
                    flag = true;
                    break;
                }

                //node_number_of_each_cell->push_back(nPoints);
                (*node_number_of_each_cell)[iCell] = nPoints;

                vector<cgsize_t> single_cell;
                single_cell.resize(0);
                single_cell.reserve(nPoints);

                for (cgsize_t jPoint = 0; jPoint < nPoints; ++ jPoint)
                {
                    single_cell.push_back(conn_list[sec][count]);
                    count ++;
                }
                //cell2node->push_back(single_cell);
                (*cell2node)[iCell] = single_cell;
            }
            if (flag) break;
        }
        else
        {
            int count = 0;
            for (cgsize_t iCell = base_cgns_istart[sec]; iCell <= base_cgns_iend[sec]; ++ iCell)
            {
                if (iCell >= nTotalCell)
                {
                    flag = true;
                    break;
                }

                int nPoints = base_elem_elem_pt[conn_list[sec][count] + 1];
                count ++;

                //node_number_of_each_cell->push_back(nPoints);
                (*node_number_of_each_cell)[iCell] = nPoints;

                vector<cgsize_t> single_cell;
                single_cell.resize(0);
                single_cell.reserve(nPoints);

                for (cgsize_t jPoint = 0; jPoint < nPoints; ++ jPoint)
                {
                    single_cell.push_back(conn_list[sec][count]);
                    count ++;
                }
                //cell2node->push_back(single_cell);
                (*cell2node)[iCell] = single_cell;
            }
            if (flag) break;
        }
    }

    if (conn_list)
    {
        for (int i = 0; i < nSections; ++ i)
        {
            delete [] conn_list[i];
        }
    }
    delete [] conn_list;

    delete [] elem_type;
    delete [] istart;
    delete [] iend;
    delete [] nBCElem;
    delete [] bc_type;
    delete [] ElementDataSize;
}

void CGNSBase::ComputeElemFace(Base_Grid_Conn *gconn, BaseElement *base_elem)
{
    std::cout << "Reconstruction faces ..." << endl;
    bool intoBCSection = false;
    cgsize_t *elem = 0;
    int count = 0;

    CElemFace *face_list = base_elem->GetElementFace();
    cgsize_t *base_elem_elem_pt = base_elem->GetElementPoint();
    for (int iSection = 0; iSection < nSections; ++ iSection)
    {
        if (istart[iSection] < nTotalCell)
        {
            //! Only dump information for interior computational cells.
            //! Do not dump for BC sections.
            std::cout << "  Process_Grid: isection = " << iSection << " nSections = " << nSections << "\n";
            std::cout << "    iStart = " << istart[iSection] << " " << ", iEnd = " << iend[iSection] << "\n";
        }
        else
        {
            std::cout << "  Get into BC sections from section " << iSection << endl;
        }

        //! Reconstruct NON-BC faces ONLY.
        for (cgsize_t iCell = istart[iSection]; iCell <= iend[iSection]; ++ iCell)
        {
            ElementType_t etype = static_cast<ElementType_t>(elem_type[iSection]);

            if (iCell >= nTotalCell)
            {
                intoBCSection = true;
                break;
            }

            ++ count;
            if (count % 200000 == 0) std::cout << "  iCell = " << count << "\n";
            
            if (etype == MIXED)
            {
                etype = static_cast<ElementType_t>(conn_list[iSection][mixed_conn_index[iSection][iCell - istart[iSection]]] + 1);
                elem = &conn_list[iSection][mixed_conn_index[iSection][iCell - istart[iSection]] + 1];
            }
            else
            {
                int ptnum = base_elem_elem_pt[etype];
                GetElem(conn_list[iSection], iCell - istart[iSection], ptnum, elem);
            }

            GetFace(gconn, elem, iCell, &face_list[etype]);
        }
        if (intoBCSection) break;
    }
}

void CGNSBase::ComputeBCFace(Base_Grid_Conn *gconn, BaseElement *base_elem)
{
    std::cout << "Boundary Faces Reconstruction ..." << endl;

    vector<cgsize_t> &right_cell_of_face = *gconn->GetRightCellOfFace();
    cgsize_t nTotalFace = right_cell_of_face.size();
    cgsize_t nBoundFace = 0;
    for (cgsize_t iFace = 0; iFace < static_cast<cgsize_t>(nTotalFace); ++ iFace)
    {
        if (right_cell_of_face[iFace] == -1)
        {
            nBoundFace ++;
        }
    }
    cgsize_t MaxCell = this->nTotalCell + nBoundFace;

    for (int iBCRegion = 0; iBCRegion < nBCRegions; ++ iBCRegion)
    {
        BCType_t bctype = static_cast<BCType_t> (bc_type[iBCRegion]);
        GridLocation_t igr = static_cast<GridLocation_t>(bc_grid_location[iBCRegion]);

        string BCName = "";
        if (boundaryName)
        {
            BCName = boundaryName[iBCRegion];
            if (bctype < NofValidBCTypes)
            {
                std::cout << "  Boundary Condition : the " << iBCRegion + 1 
                    << " BCName = "   << BCName 
                    << ", BCType = "  << BCTypeName[bctype]
                    << "\n";
            }
            else
            {
                std::cout << "  Boundary Condition : the " << iBCRegion + 1 
                    << " BCName = " << BCName << "(SelfDefined)\n";
            }
        }

        if (bc_elem_set_type[iBCRegion] == PointRange)
        {
            if (igr == EdgeCenter || igr == FaceCenter)
            {
                bc_elem_set_type[iBCRegion] = ElementRange;
            }
        }
        else if (bc_elem_set_type[iBCRegion] == PointList)
        {
            if (igr == EdgeCenter || igr == FaceCenter)
            {
                bc_elem_set_type[iBCRegion] = ElementList;
            }
        }

        if (bc_elem_set_type[iBCRegion] == ElementRange)
        {
            if (bc_conn_list[iBCRegion][0] >= MaxCell)
            {
                continue;
            }

            cgsize_t startElement = bc_conn_list[iBCRegion][0];
            cgsize_t endElement   = bc_conn_list[iBCRegion][1];
            for (cgsize_t iElement = startElement; iElement <= endElement; ++ iElement)
            {
                ProcessBCFace(gconn, base_elem, iElement, bctype, BCName);
            }
        }
        else if (bc_elem_set_type[iBCRegion] == ElementList)
        {
            if (bc_conn_list[iBCRegion][0] >= MaxCell)
            {
                continue;
            }

            for (int i = 0; i < nBCElem[iBCRegion]; ++ i)
            {
                ProcessBCFace(gconn, base_elem, bc_conn_list[iBCRegion][i], bctype, BCName);
            }
        }
        //! This situation has not been encountered, we can not ensure the code is right! 
        else if (bc_elem_set_type[iBCRegion] == PointRange ||
                 bc_elem_set_type[iBCRegion] == PointList)
        {
            ProcessBCVertex(gconn, base_elem, bc_elem_set_type[iBCRegion], 
                nBCElem[iBCRegion], bc_conn_list[iBCRegion], bctype, BCName);
        }
        else
        {
            ostringstream oss;
            oss << "Error: this situation has not been considered, for points set type is " << PointSetTypeName[bc_elem_set_type[iBCRegion]];
            oss << " in BCRegions " << iBCRegion << endl;
            TK_Exit::ExceptionExit(oss);
        }
    }

    MarkBCInterface(gconn);

    if (bc_conn_list)
    {
        for (int iBCRegion = 0; iBCRegion < nBCRegions; ++ iBCRegion)
        {
            delete [] bc_conn_list[iBCRegion];
        }
    }
    delete [] bc_conn_list;
}

void CGNSBase::MarkBCInterface(Base_Grid_Conn *gconn)
{
    vector<cgsize_t> &right_cell_of_face = *gconn->GetRightCellOfFace();

    uint_t nTotalFace = right_cell_of_face.size();
    int nBoundFace = 0;
    for (cgsize_t iFace = 0; iFace <static_cast<cgsize_t>(nTotalFace); ++ iFace)
    {
        if (right_cell_of_face[iFace] == -1)
        {
            nBoundFace ++;
        }
    }
}

void CGNSBase::ProcessBCVertex(Base_Grid_Conn *gconn, BaseElement *base_elem, 
                               int bc_elem_set_type, int nBCElem, cgsize_t *bc_conn_list, int bctype, const string &BCName)
{
    set<cgsize_t> bc_vertex;
    if (bc_elem_set_type == PointRange)
    {
        //ProcessBCVertex(gconn, base_elem, bc_elem_set_type, nBCElem, bc_conn_list, bctype);
        for (int i = static_cast<int>(bc_conn_list[0]); i <= static_cast<int>(bc_conn_list[1]); ++ i)
        {
            bc_vertex.insert(i);
        }
    }
    else
    {
        for (int i = 0; i < nBCElem; ++ i)
        {
            bc_vertex.insert(bc_conn_list[i]);
        }
    }

    ProcessBCVertex(gconn, base_elem, bc_vertex, bctype, BCName);
}

void CGNSBase::GetFace(BaseElement *base_elem, cgsize_t iCell, int &facept_num, cgsize_t *&bcface)
{
    cgsize_t iabscell = iCell;
    int isection = 0;

    GetISection(iCell, isection);
    ElementType_t etype = static_cast<ElementType_t>(elem_type[isection]);
    cgsize_t *base_elem_elem_pt = base_elem->GetElementPoint();
    facept_num = base_elem_elem_pt[etype];
    GetElem(conn_list[isection], iabscell - istart[isection], facept_num, bcface);
}

void CGNSBase::GetFace(cgsize_t *elem, int loc, int *bc_face, int &nfacept, CElemFace *elemface)
{
    int * elemface_facept_num = elemface->GetFacePointNumber();
    int **elemface_face_id    = elemface->GetFaceID();

    nfacept      = elemface_facept_num[loc];
    int *face_id = elemface_face_id[loc];

    for (int i = 0; i < nfacept; ++ i)
    {
        bc_face[i] = static_cast<int>(elem[face_id[i]]);
    }
}

bool CGNSBase::CheckBCFace(int *bc_face, int nfacept, set<cgsize_t> &bc_vertex)
{
    set<cgsize_t>::iterator iter;
    for (int i = 0; i < nfacept; ++ i)
    {
        iter = bc_vertex.find(bc_face[i]);
        if (iter == bc_vertex.end())
        {
            return false;
        }
    }
    return true;
}

void CGNSBase::ProcessBCVertex(Base_Grid_Conn *gconn, BaseElement *base_elem, 
                               set<cgsize_t> &bc_vertex, int bctype, const string &BCName)
{
    vector<cgsize_t> &left_cell_of_face  = *gconn->GetLeftCellOfFace();
    vector<cgsize_t> &right_cell_of_face = *gconn->GetRightCellOfFace();

    uint_t nTotalFace = right_cell_of_face.size();
    int nBoundFace = 0;
    for (cgsize_t iFace = 0; iFace < static_cast<cgsize_t>(nTotalFace); ++ iFace)
    {
        if (right_cell_of_face[iFace] == -1)
        {
            nBoundFace ++;
        }
    }

    cout << "ProcessBCVertex BCtype = " << bctype << " ........................................\n";
    cout << "All Candidate bc num is : " << nBoundFace << "\n";

    cgsize_t *elem;
    cgsize_t iCell;
    int etype, loc, nfacept;
    const int maxfacept = 100;
    int bc_face[maxfacept];

    vector<CGNS_BCFace *> *bcFace = gconn->GetBCFace();
    vector<int> *face_location = gconn->GetFaceLocation();

    CElemFace *face_list = base_elem->GetElementFace();

    for (cgsize_t iFace = 0; iFace < static_cast<cgsize_t>(nTotalFace); ++ iFace)
    {
        int old_bc_type = (*bcFace)[iFace]->GetBCType();
        if (right_cell_of_face[iFace] == - 1 && old_bc_type == 0)
        {
            iCell = left_cell_of_face[iFace];
            loc   = (*face_location)[iFace];

            GetElem(base_elem, iCell, etype, elem);
            GetFace(elem, loc, bc_face, nfacept, &face_list[etype]);

            if (CheckBCFace(bc_face, nfacept, bc_vertex))
            {
                (*bcFace)[iFace]->SetBCType(bctype);
                (*bcFace)[iFace]->SetBCName(BCName);
                if (bctype == 0)
                {
                    cout << "impossible bctype == 0 in ProcessBCVertex\n";
                }
            }
        }
    }
}

CGNSFactory::CGNSFactory(int nBlocks)
{
    this->nBlocks = nBlocks;
    base_conn = new Base_Grid_Conn *[nBlocks];
    for (int iZone = 0; iZone < nBlocks; ++ iZone)
    {
        base_conn[iZone] = new Base_Grid_Conn();
    }

    rawgrid = new RawGrid *[nBlocks];
    for (int iZone = 0; iZone < nBlocks; ++ iZone)
    {
        rawgrid[iZone] = new RawGrid();
    }

    base_cgns = new CGNSBase *[nBlocks];
    for (int iZone = 0; iZone < nBlocks; ++ iZone)
    {
        base_cgns[iZone] = new CGNSBase();
    }

    base_elem = new BaseElement();
    fantasy_bc = new BaseBCType();
}

CGNSFactory::~CGNSFactory()
{
    delete base_elem;
    delete fantasy_bc;

    for (int iZone = 0; iZone < nBlocks; ++ iZone)
    {
        delete base_conn[iZone];
    }
    delete [] base_conn;

    for (int iZone = 0; iZone < nBlocks; ++ iZone)
    {
        delete rawgrid[iZone];
    }
    delete [] rawgrid;

    //for (int iZone = 0; iZone < nBlocks; ++ iZone)
    //{
    //    delete base_cgns[iZone];
    //}
    delete [] base_cgns;
}

void CGNSFactory::ConvertGrid2Fantasy()
{
    //ReLableSelfDefineBC();

    for (int iZone = 0; iZone < nBlocks; ++ iZone)
    {
        std::cout << "Converting Grid to PHengLEI for Zone " << iZone << " of Total Blocks " << nBlocks << "\n";
        ComputeGridConn(this->GetCGNSBase(iZone), base_conn[iZone]);
    }

    DumpSelfDefinePartInformation();
}

void CGNSFactory::ComputeGridConn(CGNSBase *cgnsgrid, Base_Grid_Conn *gconn)
{
    cgnsgrid->ComputeGridConn(gconn, this->base_elem);
}

void CGNSFactory::ReLableSelfDefineBC()
{
    int countOfSelfDefine = 0;
    for (int iBlock = 0; iBlock < nBlocks; ++ iBlock)
    {
        CGNSBase *cgnsGrid = this->GetCGNSBase(iBlock);

        int *bc_type = cgnsGrid->GetBCType();
        int nBCRegions = cgnsGrid->GetNumberOfBCRegions();

        for (int iBCRegion = 0; iBCRegion < nBCRegions; ++ iBCRegion)
        {
            BCType_t bctype = static_cast<BCType_t>(bc_type[iBCRegion]);
            //! Since the number of custom boundary condition is 1, all of them are forced to add the maximum number NofValidBCTypes.
            if (bctype == BCTypeUserDefined)
            {
                bctype = static_cast<BCType_t> (NofValidBCTypes + countOfSelfDefine);
                ++ countOfSelfDefine;
            }

            if (bctype >= NofValidBCTypes)
            {
                bc_type[iBCRegion] = bctype;
            }
        }
    }
}

void CGNSFactory::DumpSelfDefinePartInformation()
{
    //! Collect all of the boundary condition.
    set<string> nameContainer;
    vector<int> bcType;
    vector<string> boundaryName;

    for (int iZone = 0; iZone < nBlocks; ++ iZone)
    {
        CGNSBase *cgnsgrid = this->GetCGNSBase(iZone);

        int *bc_type = cgnsgrid->GetBCType();
        string *bc_name = cgnsgrid->GetBCName();
        int nBCRegions = cgnsgrid->GetNumberOfBCRegions();

        for (int iBC = 0; iBC < nBCRegions; ++ iBC)
        {
            int cgnsBCType = bc_type[iBC];
            int fantasyBCType = this->GetBaseBCType()->GetBCType(cgnsBCType);

            string name = bc_name[iBC];

            set<string> :: iterator nameIter;

            nameIter = nameContainer.find(name);
            //if(nameIter == nameContainer.end())
            if (nameIter == nameContainer.end() && fantasyBCType >= PHENGLEI::StartNumberOfBodyOfHyperFLOW && fantasyBCType <= PHENGLEI::ENDNumberOfBodyOfHyperFLOW)
            {
                boundaryName.push_back(name);
                bcType.push_back(fantasyBCType);
            }
        }
    }

    uint_t nPart = bcType.size();
    if (!nPart)
    {
        return;
    }

    //! Dump out.
    string fileName;
    GlobalDataBase::GetData("from_gfile", &fileName, PHSTRING, 1);
    fileName = ChangeExtensionOfFileName(fileName, "part");

    fstream outfile(fileName.c_str(), ios::out);
    if (!outfile)
    {
        ostringstream oss;
        oss << "Can not open file " << fileName << endl;
        TK_Exit::ExceptionExit(oss.str());
    }

    RDouble dm = 0.0, dl = 0.0, dn3 = 0.0;
    RDouble Lref = 1.0, Sref = 1.0, xref = 0.0, yref = 0.0, zref = 0.0;

    outfile << nPart << endl;

    outfile << dm << " " << dl << " " << dn3 << endl;

    vector< int > :: iterator typeIter = bcType.begin();
    vector< string > :: iterator nameIter = boundaryName.begin();

    for (int iPart = 0; iPart < nPart; ++ iPart)
    {
        int type = *typeIter;
        string name = *nameIter;

        outfile << name << "	" << type << endl;
        outfile << "		" << xref << "	" << yref << "	" << zref << "	" << Sref << "	" << Lref << endl;

        ++ typeIter;
        ++ nameIter;
    }

    outfile.close();
    outfile.clear();
}

void CGNS2UnsGrid(CGNSFactory *factory_cgns, int iZone, UnstructGrid *grid)
{
    Base_Grid_Conn *gconn   = factory_cgns->GetBaseGridConnection(iZone);
    RawGrid        *rawgrid = factory_cgns->GetRawGrid(iZone);
    CGNSBase *cgnsGrid = factory_cgns->GetCGNSBase(iZone);

    int nTotalCell = static_cast<int>(gconn->GetNTotalCell());
    int nTotalNode = static_cast<int>(gconn->GetNTotalNode());
    int nBoundFace = static_cast<int>(gconn->GetNBoundFace());

    grid->SetNTotalCell(nTotalCell);
    grid->SetNTotalNode(nTotalNode);
    grid->SetNBoundFace(nBoundFace);

    RDouble *x = new RDouble[nTotalNode];
    RDouble *y = new RDouble[nTotalNode];
    RDouble *z = new RDouble[nTotalNode];

    grid->SetX(x);
    grid->SetY(y);
    grid->SetZ(z);

    RDouble *xx = rawgrid->GetX();
    RDouble *yy = rawgrid->GetY();
    RDouble *zz = rawgrid->GetZ();

    for (int iNode = 0; iNode < nTotalNode; ++ iNode)
    {
        x[iNode] = xx[iNode];
        y[iNode] = yy[iNode];
        z[iNode] = zz[iNode];
    }

    delete [] xx;
    delete [] yy;
    delete [] zz;

    cout << "  Boundary condition setting ......\n";
    //! Set the boundary conditions.
    HyList<int> cgnsBCList;
    HyList<int> fantasyBCList;

    map<int, string> Fantasybctype = CreatFantasybctype();

    vector<CGNS_BCFace *> *bcFace = gconn->GetBCFace();

    int nBCRegionUnstruct = cgnsGrid->GetNumberOfBCRegions();
    string *bcRegionName = cgnsGrid->GetBCName();

    for (cgsize_t iFace = 0; iFace < nBoundFace; ++iFace)
    {
        const string &bcName = (*bcFace)[iFace]->GetBCName();
        int cgnsBCType = (*bcFace)[iFace]->GetBCType();
        int fantasyBCType = factory_cgns->GetBaseBCType()->GetBCType(cgnsBCType);

        cgnsBCList.insert(cgnsBCType);
        fantasyBCList.insert(fantasyBCType);

        delete (*bcFace)[iFace];
    }
    string *boundaryName = cgnsGrid->GetBCName();

    grid->CreateUnstructBCSet(nBCRegionUnstruct);
    UnstructBCSet *unstructBCSet = grid->GetUnstructBCSet();
    int* bcRegionIDofBCFace = new int[nBoundFace];
    for (int iBCRegionUnstruct = 0; iBCRegionUnstruct < nBCRegionUnstruct; iBCRegionUnstruct++)
    {
        UnstructBC *unstructBC = new UnstructBC(iBCRegionUnstruct);
        unstructBCSet->SetBCRegion(iBCRegionUnstruct, unstructBC);
        unstructBC->SetBCName(boundaryName[iBCRegionUnstruct]);
        
        for (cgsize_t iface = 0; iface < nBoundFace; ++iface)
        {
            const string &bcName = (*bcFace)[iface]->GetBCName();
            int cgnsBCType = (*bcFace)[iface]->GetBCType();
            int fantasyBCType = factory_cgns->GetBaseBCType()->GetBCType(cgnsBCType);
  
            cgnsBCList.insert(cgnsBCType);
            fantasyBCList.insert(fantasyBCType);

            cgnsBCList.insert(cgnsBCType);
            fantasyBCList.insert(fantasyBCType);
            if (bcName == unstructBC->GetBCName())
            {
                unstructBC->SetBCType(fantasyBCType);
                bcRegionIDofBCFace[iface] = iBCRegionUnstruct;
                vector<int> *faceIndex = unstructBC->GetFaceIndex();
                faceIndex->push_back(iface);
            }
        }
    }
    unstructBCSet->SetBCFaceInfo(bcRegionIDofBCFace);
    cout << "After converting CGNS grid,\n";
    cout << "Number of Original Boundary Condition: " << cgnsBCList.size() << "\n";
    for (int iBC = 0; iBC < cgnsBCList.size(); ++ iBC)
    {
        cout<< "  The "<< iBC + 1 << " CGNS BC Name is: " << BCTypeName[cgnsBCList[ iBC ]] << endl;
    }

    //map<int, string> Fantasybctype = CreatFantasybctype();

    cout << "Number of Final Boundary Condition: " << fantasyBCList.size() << "\n";
    for (int iBC = 0; iBC < fantasyBCList.size(); ++ iBC)
    {
        cout << "  The " << iBC + 1 << " PHengLEI BC Name is: ";
        cout << Fantasybctype[fantasyBCList[iBC]] << "\n";
    }
    //     for (set<int>::iterator iter = s2.begin(); iter != s2.end(); ++ iter)
    //     {
    //         cout << *iter << " ";
    //     }
    //cout << "\n";

    vector<cgsize_t> *gconn_left_cell_of_face  = gconn->GetLeftCellOfFace();
    vector<cgsize_t> *gconn_right_cell_of_face = gconn->GetRightCellOfFace();

    vector<int> *gconn_node_number_of_each_face = gconn->GetNodeNumberOfEachFace();

    int nTotalFace = static_cast<int>(gconn_left_cell_of_face->size());
    grid->SetNTotalFace(nTotalFace);
    cout << "Number of total faces: " << nTotalFace << "\n";

    int *node_number_of_each_face = new int [nTotalFace];
    for (cgsize_t iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        node_number_of_each_face[iFace] = (*gconn_node_number_of_each_face)[iFace];
    }
    grid->SetNodeNumberOfEachFace(node_number_of_each_face);

    gconn->FreeNodeNumberOfEachFace();

    long long int face2node_size = 0;
    for (cgsize_t iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        face2node_size += node_number_of_each_face[iFace];
    }
    int *face2node = new int [face2node_size];
    grid->SetFace2Node(face2node);

    vector< vector<cgsize_t> > *gface2node = gconn->GetFace2Node();

    long long int count = 0;
    for (cgsize_t iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        for (std::size_t j = 0; j < (*gface2node)[iFace].size(); ++ j)
        {
            face2node[count] = static_cast<int>((*gface2node)[iFace][j]);
            ++ count;
        }
    }
    gconn->FreeFace2Node();
    //cout << " Set face2node finished\n";

    int *left_cell_of_face  = new int[nTotalFace];
    int *right_cell_of_face = new int[nTotalFace];
    grid->SetLeftCellOfFace (left_cell_of_face);
    grid->SetRightCellOfFace(right_cell_of_face);
    for (cgsize_t iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        left_cell_of_face [iFace] = static_cast<int>((*gconn_left_cell_of_face)[iFace]);
        //! Set the label of ghostcell.
        right_cell_of_face[iFace] = static_cast<int>((*gconn_right_cell_of_face)[iFace]);
    }

    gconn->FreeCellOfFace();

    int *node_number_of_each_cell = new int [nTotalCell];

    vector<int> *gconn_node_number_of_each_cell = gconn->GetNodeNumberOfEachCell();

    for (cgsize_t iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        node_number_of_each_cell[iCell] = (*gconn_node_number_of_each_cell)[iCell];
    }
    grid->SetNodeNumberOfEachCell(node_number_of_each_cell);
    //cout << " Set node_number_of_each_cell finished\n";

    gconn->FreeNodeNumberOfEachCell();

    long long int cell2node_size =0;

    for (cgsize_t iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        cell2node_size += node_number_of_each_cell[iCell];
    }
    int *cell2node = new int [cell2node_size];
    grid->SetCell2Node(cell2node);

    vector< vector<cgsize_t> > *gcell2node = gconn->GetCell2Node();

    count = 0;
    for (cgsize_t iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        for (int iNode =0; iNode < node_number_of_each_cell[iCell]; ++ iNode)
        {
            cell2node[count] = static_cast<int>((*gcell2node)[iCell][iNode]);
            ++ count;
        }
    }
    gconn->FreeCell2Node();
}

}
