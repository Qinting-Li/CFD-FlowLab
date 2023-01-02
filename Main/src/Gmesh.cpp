#include "Gmesh.h"
#include "TK_Parse.h"
#include "Math_BasisFunction.h"
#include "Pre_GridBase.h"
using namespace std;

#pragma warning (disable:913)
namespace PHSPACE
{
void ElemReorderGmsh2CGNS(int etype, int *node, int *node_cgns, int nvertex)
{
    for (int i = 0; i < nvertex; ++ i)
    {
        node_cgns[i] = node[i];
    }

    if (etype == PYRA_14)
    {
        //! Pyramid
        node_cgns[6 ] = node[8 ];
        node_cgns[7 ] = node[10]; 
        node_cgns[8 ] = node[6 ];
        node_cgns[9 ] = node[7 ];
        node_cgns[10] = node[9 ];
    }
    else if (etype == TETRA_10)
    {
        //! Tetrahedron
        node_cgns[8 ] = node[9 ];
        node_cgns[9 ] = node[8 ]; 
    }
    else if (etype == PENTA_15)
    {
        //! TriangularPrism
        node_cgns[7 ] = node[9 ];
        node_cgns[8 ] = node[7 ];
        node_cgns[9 ] = node[8 ];
        node_cgns[13] = node[14];
        node_cgns[14] = node[13];
    }
    else if (etype == PENTA_18)
    {
        //! TriangularPrism
        node_cgns[7 ] = node[9 ];
        node_cgns[8 ] = node[7 ];
        node_cgns[9 ] = node[8 ];
        node_cgns[13] = node[14];
        node_cgns[14] = node[13];
        node_cgns[16] = node[17];
        node_cgns[17] = node[16];
    }
    else if (etype == HEXA_20)
    {
        //! Hexahedron
        node_cgns[9 ] = node[11];
        node_cgns[10] = node[13];
        node_cgns[11] = node[9 ];
        node_cgns[12] = node[10];
        node_cgns[13] = node[12];
        node_cgns[17] = node[18];
        node_cgns[18] = node[19];
        node_cgns[19] = node[17];
    }
    else if (etype == HEXA_27)
    {
        //! Hexahedron
        node_cgns[9 ] = node[11];
        node_cgns[10] = node[13];
        node_cgns[11] = node[9 ];
        node_cgns[12] = node[10];
        node_cgns[13] = node[12];
        node_cgns[17] = node[18];
        node_cgns[18] = node[19];
        node_cgns[19] = node[17];

        node_cgns[22] = node[23];
        node_cgns[23] = node[24];
        node_cgns[24] = node[22];
    }
    else
    {
        ;
    }
}

GmeshElem::GmeshElem()
{
    ntype      = 0;
    conn_list  = 0;
    type       = 0;
    istart     = 0;
    iend       = 0;
    number     = 0;
    nTotalCell = 0;
}

GmeshElem::~GmeshElem()
{
    //conn_list = new int * [ntype];
    for (int i = 0; i < ntype; ++ i)
    {
        delete [] conn_list[i];
    }
    delete [] conn_list;
    delete [] type;
    delete [] istart;
    delete [] iend;
    delete [] number;
}

void GmeshElem::Create(int ntype)
{
    this->ntype = ntype;
    conn_list = new int * [ntype];
    type      = new int [ntype];
    istart    = new int [ntype];
    iend      = new int [ntype];
    number    = new int [ntype];
}

void CreateValueSet(vector<int> &value_container, set<int> &value_set)
{
    for (std::size_t i = 0; i < value_container.size(); ++ i)
    { 
        value_set.insert(value_container[i]);
    }
}

void GmeshElem::TreatElements(CGNSFactory *factory_cgns, vector<int> &type_container, vector<int> &number_container, vector<int> &index_container, vector<int> &bctype_container, int *e_nvertex)
{
    //找出所有不同的单元类型
    set<int> elem_type_set;
    for (std::size_t i = 0; i < type_container.size(); ++ i)
    { 
        int etype = type_container[i];
        elem_type_set.insert(etype);
    }

    int nSections = static_cast<int>(elem_type_set.size());

    this->Create(nSections);

    Set2Array(elem_type_set, type);

    for (int m = 0; m < nSections; ++ m)
    {
        int nelem_this_type = 0;
        for (std::size_t i = 0; i < type_container.size(); ++ i)
        {
            int etype     = type_container[i];
            int nelem_num = number_container[i];
            if (etype == type[m])
            {
                nelem_this_type += nelem_num;
            }
        }
        int nvertex  = e_nvertex[ type[m] ];
        conn_list[m] = new int[ nelem_this_type * nvertex ];
        number[m]    = nelem_this_type;
    }

    nTotalCell = 0;
    for (int m = 0; m < nSections; ++ m)
    {
        nTotalCell += number[m];
    }

    istart[0] = 1;
    for (int m = 0; m < nSections; ++ m)
    {
        iend[m] = istart[m] + number[m] - 1;
        if (m < nSections - 1)
        {
            istart[m+1] = iend[m] + 1;
        }
    }

    for (int m = 0; m < nSections; ++ m)
    {
        int icount = 0;
        int jcount = 0;
        for (std::size_t i = 0; i < type_container.size(); ++ i)
        {
            int etype     = type_container[i];
            int nelem_num = number_container[i];
            int nvertex   = e_nvertex[ etype ];
            int ptmax     = nelem_num * nvertex;

            if (etype == type[m])
            {
                for (int j = 0; j < ptmax; ++ j)
                {
                    conn_list[m][icount+j] = index_container[jcount+j];
                }
                icount += ptmax;
            }
            jcount += ptmax;
        }
    }

    //找出所有不同的边界类型
    set<int> bc_type_set;
    CreateValueSet(bctype_container, bc_type_set);

    //删除内边界
    bc_type_set.erase(0);

    //! Determine the number of boundary conditions for this zone.
    int nBCRegions = static_cast<int>(bc_type_set.size());
    int *bctype_array = new int [nBCRegions];

    Set2Array(bc_type_set, bctype_array);

    CGNSBase *base_cgns = factory_cgns->GetCGNSBase(0);
    base_cgns->CreateBC(nBCRegions);

    int * base_cgns_nBCElem          = base_cgns->GetNumberOfBCElements();
    int * base_cgns_bc_elem_set_type = base_cgns->GetBoundaryElementType();
    int * base_cgns_bc_grid_location = base_cgns->GetBoundaryGridLocation();
    int **base_cgns_bc_conn_list     = (int**) base_cgns->GetBoundaryElementConnectionList();
    cout << "    Warning: the cgsize_t type of conn_list has been convert to int type, when deal with GMSH grid!!!" 
         << endl;
    PointSetType_t ptsetType = PointList;

    //! Loop over the number of boundary conditions.
    for (int iBCRegion = 0; iBCRegion < nBCRegions; ++ iBCRegion)
    {
        GridLocation_t igr = FaceCenter;
        cout << " Grid Location Name = " << GridLocationName[igr] << "\n";

        if (igr == FaceCenter)
        {
            cout << "\nGridLocation = FaceCenter means BC data refers to elements, not nodes\n";
        }
        else if (igr == Vertex)
        {
            cout << "\nGridLocation = Vertex means BC data refers to nodes, not elements\n";
        }

        int bctype = bctype_array[iBCRegion];
        int nBCElem = 0;
        for (std::size_t i = 0; i < bctype_container.size(); ++ i)
        {
            if (bctype == bctype_container[i])
            {
                ++ nBCElem;
            }
        }

        base_cgns_nBCElem         [iBCRegion] = nBCElem;
        base_cgns_bc_conn_list    [iBCRegion] = new int [nBCElem];
        //base_cgns->bc_type        [iBCRegion] = bocoType;
        base_cgns_bc_elem_set_type[iBCRegion] = ptsetType;
        base_cgns_bc_grid_location[iBCRegion] = igr;

        cout << "PointSetType_t = " << ptsetType << "\n";

        if (nBCElem == 2)
        {
            cout << base_cgns_bc_conn_list[iBCRegion][0] << " " << base_cgns_bc_conn_list[iBCRegion][1] << "\n";
        }
        else
        {
            for (int j = 0; j < MIN(20, nBCElem); ++ j)
            {
                cout << base_cgns_bc_conn_list[iBCRegion][j] << " ";
            }
            cout << "\n";
        }
    }

    delete [] bctype_array;
}

Gmesh::Gmesh()
{
    //vector<string> keyword;
    keyword.push_back("MeshFormat");
    keyword.push_back("Nodes");
    keyword.push_back("$Elements");
    keyword.push_back("PhysicalNames");
    keyword.push_back("NodeData");
    keyword.push_back("ElementData");
    keyword.push_back("ElementNodeData");

    elem_pt = new int [NofValidElementTypes];

    elem_pt[NODE   ] = 1;
    elem_pt[BAR_2  ] = 2;
    elem_pt[TRI_3  ] = 3;
    elem_pt[QUAD_4 ] = 4;
    elem_pt[TETRA_4] = 4;
    elem_pt[PYRA_5 ] = 5;
    elem_pt[PENTA_6] = 6;
    elem_pt[HEXA_8 ] = 8;

    NofValidElementTypesGmsh = MSH_NUM_TYPE;
    TypeGmsh2CGNS = new int [NofValidElementTypesGmsh+1];

    InitGmsh2CGNSMap();

    gmesh_nvertex = new int [NofValidElementTypesGmsh];

    InitGmshNVertex();

    cgns_nvertex = new int [NofValidElementTypesGmsh];

    InitCGNSNVertex();
}

Gmesh::~Gmesh()
{
    delete [] elem_pt;
    delete [] TypeGmsh2CGNS;
    delete [] gmesh_nvertex;
    delete [] cgns_nvertex;
}

void Gmesh::InitGmshNVertex()
{
    for (int i = 0; i < NofValidElementTypesGmsh; ++ i)
    {
        gmesh_nvertex[i] = 1;
    }

    gmesh_nvertex[ PHSPACE::PH_MSH_LIN_2   ] = 2;
    gmesh_nvertex[ PHSPACE::PH_MSH_TRI_3   ] = 3;
    gmesh_nvertex[ PHSPACE::PH_MSH_QUA_4   ] = 4;
    gmesh_nvertex[ PHSPACE::PH_MSH_TET_4   ] = 4;
    gmesh_nvertex[ PHSPACE::PH_MSH_HEX_8   ] = 8;
    gmesh_nvertex[ PHSPACE::PH_MSH_PRI_6   ] = 6;
    gmesh_nvertex[ PHSPACE::PH_MSH_PYR_5   ] = 5;
    gmesh_nvertex[ PHSPACE::PH_MSH_LIN_3   ] = 3;
    gmesh_nvertex[ PHSPACE::PH_MSH_TRI_6   ] = 6;
    gmesh_nvertex[ PHSPACE::PH_MSH_QUA_9   ] = 9;
    gmesh_nvertex[ PHSPACE::PH_MSH_TET_10  ] = 10;
    gmesh_nvertex[ PHSPACE::PH_MSH_HEX_27  ] = 27;
    gmesh_nvertex[ PHSPACE::PH_MSH_PRI_18  ] = 18;
    gmesh_nvertex[ PHSPACE::PH_MSH_PYR_14  ] = 14;
    gmesh_nvertex[ PHSPACE::PH_MSH_PNT     ] = 1;
    gmesh_nvertex[ PHSPACE::PH_MSH_QUA_8   ] = 8;
    gmesh_nvertex[ PHSPACE::PH_MSH_HEX_20  ] = 20;
    gmesh_nvertex[ PHSPACE::PH_MSH_PRI_15  ] = 15;
    gmesh_nvertex[ PHSPACE::PH_MSH_PYR_13  ] = 13;
}

void Gmesh::InitGmsh2CGNSMap()
{
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_LIN_2   ] = BAR_2;
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TRI_3   ] = TRI_3;
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_QUA_4   ] = QUAD_4;
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TET_4   ] = TETRA_4;
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_HEX_8   ] = HEXA_8;
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_PRI_6   ] = PENTA_6;
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_PYR_5   ] = PYRA_5;
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_LIN_3   ] = BAR_3;
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TRI_6   ] = TRI_6;
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_QUA_9   ] = QUAD_9;
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TET_10  ] = TETRA_10;
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_HEX_27  ] = HEXA_27;
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_PRI_18  ] = PENTA_18;
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_PYR_14  ] = PYRA_14;
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_PNT     ] = NODE;
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_QUA_8   ] = QUAD_8;
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_HEX_20  ] = HEXA_20;
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_PRI_15  ] = PENTA_15;
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_PYR_13  ] = PYRA_13;
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TRI_9   ] = TRI_6;       //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TRI_10  ] = TRI_6;       //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TRI_12  ] = TRI_6;       //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TRI_15  ] = TRI_6;       //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TRI_15I ] = TRI_6;       //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TRI_21  ] = TRI_6;       //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_LIN_4   ] = BAR_3;       //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_LIN_5   ] = BAR_3;       //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_LIN_6   ] = BAR_3;       //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TET_20  ] = TETRA_10;    //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TET_35  ] = TETRA_10;    //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TET_56  ] = TETRA_10;    //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TET_34  ] = TETRA_10;    //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TET_52  ] = TETRA_10;    //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_POLYG_  ] = NGON_n;      //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_POLYH_  ] = NGON_n;      //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_QUA_16  ] = QUAD_9;      //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_QUA_25  ] = QUAD_9;      //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_QUA_36  ] = QUAD_9;      //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_QUA_12  ] = QUAD_9;      //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_QUA_16I ] = QUAD_9;      //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_QUA_20  ] = QUAD_9;      //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TRI_28  ] = TRI_6;       //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TRI_36  ] = TRI_6;       //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TRI_45  ] = TRI_6;       //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TRI_55  ] = TRI_6;       //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TRI_66  ] = TRI_6;       //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_QUA_49  ] = QUAD_9;      //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_QUA_64  ] = QUAD_9;      //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_QUA_81  ] = QUAD_9;      //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_QUA_100 ] = QUAD_9;      //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_QUA_121 ] = QUAD_9;      //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TRI_18  ] = TRI_6;       //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TRI_21I ] = TRI_6;       //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TRI_24  ] = TRI_6;       //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TRI_27  ] = TRI_6;       //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TRI_30  ] = TRI_6;       //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_QUA_24  ] = QUAD_9;      //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_QUA_28  ] = QUAD_9;      //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_QUA_32  ] = QUAD_9;      //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_QUA_36I ] = QUAD_9;      //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_QUA_40  ] = QUAD_9;      //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_LIN_7   ] = BAR_3;       //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_LIN_8   ] = BAR_3;       //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_LIN_9   ] = BAR_3;       //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_LIN_10  ] = BAR_3;       //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_LIN_11  ] = BAR_3;       //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_LIN_B   ] = BAR_3;       //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TRI_B   ] = TRI_6;       //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_POLYG_B ] = NGON_n;      //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_LIN_C   ] = BAR_3;       //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TET_84  ] = TETRA_10;    //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TET_120 ] = TETRA_10;    //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TET_165 ] = TETRA_10;    //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TET_220 ] = TETRA_10;    //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TET_286 ] = TETRA_10;    //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_HEX_64  ] = HEXA_20;     //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_HEX_125 ] = HEXA_20;     //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_HEX_196 ] = HEXA_20;     //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TET_74  ] = TETRA_10;    //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TET_100 ] = TETRA_10;    //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TET_130 ] = TETRA_10;    //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TET_164 ] = TETRA_10;    //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_TET_202 ] = TETRA_10;    //CGNS没有对应类型
    TypeGmsh2CGNS[ PHSPACE::PH_MSH_NUM_TYPE] = NGON_n;      //CGNS没有对应类型
}

void Gmesh::InitCGNSNVertex()
{
    for (int i = 0; i < NofValidElementTypes; ++ i)
    {
        cgns_nvertex[i] = 1;
    }

    cgns_nvertex[ ElementTypeNull         ] = 0;
    cgns_nvertex[ ElementTypeUserDefined  ] = 0;
    cgns_nvertex[ NODE                    ] = NPE_NODE;
    cgns_nvertex[ BAR_2                   ] = NPE_BAR_2;
    cgns_nvertex[ BAR_3                   ] = NPE_BAR_3;
    cgns_nvertex[ TRI_3                   ] = NPE_TRI_3;
    cgns_nvertex[ TRI_6                   ] = NPE_TRI_6;
    cgns_nvertex[ QUAD_4                  ] = NPE_QUAD_4;
    cgns_nvertex[ QUAD_8                  ] = NPE_QUAD_8;
    cgns_nvertex[ QUAD_9                  ] = NPE_QUAD_9;
    cgns_nvertex[ TETRA_4                 ] = NPE_TETRA_4;
    cgns_nvertex[ TETRA_10                ] = NPE_TETRA_10;
    cgns_nvertex[ PYRA_5                  ] = NPE_PYRA_5;
    cgns_nvertex[ PYRA_14                 ] = NPE_PYRA_13;
    cgns_nvertex[ PENTA_6                 ] = NPE_PYRA_14;
    cgns_nvertex[ PENTA_15                ] = NPE_PENTA_6;
    cgns_nvertex[ PENTA_18                ] = NPE_PENTA_15; 
    cgns_nvertex[ HEXA_8                  ] = NPE_PENTA_18;
    cgns_nvertex[ HEXA_20                 ] = NPE_HEXA_8;
    cgns_nvertex[ HEXA_27                 ] = NPE_HEXA_20;
    cgns_nvertex[ MIXED                   ] = NPE_HEXA_27;
    cgns_nvertex[ PYRA_13                 ] = NPE_MIXED;
    cgns_nvertex[ NGON_n                  ] = NPE_NGON_n;
    cgns_nvertex[ NFACE_n                 ] = NPE_NFACE_n;
}

void Gmesh::Read(fstream &file, CGNSFactory *factory_cgns)
{
    ReadMeshFormat(file);
    ReadNodes(file, factory_cgns);
    ReadElements(file, factory_cgns);
}

void Gmesh::ReadMeshFormat(fstream &file)
{
    string line,word;
    string separator = " =\t\r\n#$,;";
    do
    {
        getline(file,line);
        line = FindNextWord(line,word,separator);
        transform(word.begin(), word.end(), word.begin(), ToLower());
    } while (word != "meshformat");

    getline(file,line);

    line = FindNextWord(line,word,separator);
    from_string< RDouble >(version_number,word,std::dec);

    line = FindNextWord(line,word,separator);
    from_string< int >(file_type,word,std::dec);

    line = FindNextWord(line,word,separator);
    from_string< int >(data_size,word,std::dec);


    cout << "version-number = " << version_number << "\n";
    cout << "file-type      = " << file_type      << "\n";
    cout << "data-size      = " << data_size      << "\n";
}

void Gmesh::ReadNodes(fstream &file, CGNSFactory *factory_cgns)
{
    string line,word;
    string separator = " =\t\r\n#$,;";
    do
    {
        getline(file,line);
        line = FindNextWord(line,word,separator);
        transform(word.begin(), word.end(), word.begin(), ToLower());
    } while (word != "nodes");

    getline(file,line);

    line = FindNextWord(line,word,separator);
    from_string< int >(nTotalNode,word,std::dec);

    int nBlocks = 1;

    for (int iZone = 0; iZone < nBlocks; ++ iZone)
    {
        cout << " iZone = " << iZone << " nBlocks = " << nBlocks << "\n";

        RDouble *x = new RDouble[ nTotalNode ];
        RDouble *y = new RDouble[ nTotalNode ];
        RDouble *z = new RDouble[ nTotalNode ];

        RawGrid *rawgrid = factory_cgns->GetRawGrid(iZone);

        rawgrid->SetX(x);
        rawgrid->SetY(y);
        rawgrid->SetZ(z);

        RDouble xmin,ymin,zmin,xmax,ymax,zmax;
        xmin = ymin = zmin =   LARGE;
        xmax = ymax = zmax = - LARGE;

        cout << " Read grid coordinates \n";

        int iskip = 500000;
        int index;

        for (int iNode = 0; iNode < nTotalNode; ++ iNode)
        {
            if (iNode % iskip == 0)
            {
                cout << " iNode = " << iNode << " nTotalNode = " << nTotalNode << "\n";
            }
            getline(file,line);
            line = FindNextWord(line,word,separator);
            from_string< int >(index,word,std::dec);

            line = FindNextWord(line,word,separator);
            from_string< RDouble >(x[iNode],word,std::dec);

            line = FindNextWord(line,word,separator);
            from_string< RDouble >(y[iNode],word,std::dec);

            line = FindNextWord(line,word,separator);
            from_string< RDouble >(z[iNode],word,std::dec);

            xmin = MIN(xmin, x[iNode]);
            ymin = MIN(ymin, y[iNode]);
            zmin = MIN(zmin, z[iNode]);

            xmax = MAX(xmax, x[iNode]);
            ymax = MAX(ymax, y[iNode]);
            zmax = MAX(zmax, z[iNode]);
        }

        cout << setiosflags(ios::right);
        cout << setprecision(8);
        cout << setiosflags(ios::scientific);
        cout << setiosflags(ios::showpoint);
        int wordwidth = 16;
        cout << " xmin = " << setw(wordwidth) << xmin << " xmax = " << setw(wordwidth) << xmax << "\n";
        cout << " ymin = " << setw(wordwidth) << ymin << " ymax = " << setw(wordwidth) << ymax << "\n";
        cout << " zmin = " << setw(wordwidth) << zmin << " zmax = " << setw(wordwidth) << zmax << "\n";
    }
}

//$MeshFormat
//2.2 0 8
//$EndMeshFormat
//$Nodes
//6 six mesh nodes:
//1 0.0 0.0 0.0 node #1: coordinates (0.0, 0.0, 0.0)
//2 1.0 0.0 0.0 node #2: coordinates (1.0, 0.0, 0.0)
//3 1.0 1.0 0.0 etc.
//4 0.0 1.0 0.0
//5 2.0 0.0 0.0
//6 2.0 1.0 0.0
//$EndNodes
//$Elements
//2 two elements:
//1 3 2 99 2 1 2 3 4 quad #1: type 3, physical 99, elementary 2, nodes 1 2 3 4
//2 3 2 99 2 2 5 6 3 quad #2: type 3, physical 99, elementary 2, nodes 2 5 6 3
//$EndElements
//$NodeData
//1 one string tag:
//"A scalar view" the name of the view ("A scalar view")
//1 one real tag:
//0.0 the time value (0.0)
//3 three integer tags:
//0 the time step (0; time steps always start at 0)
//1 1-component (scalar) field
//6 six associated nodal values
//1 0.0 value associated with node #1 (0.0)
//2 0.1 value associated with node #2 (0.1)
//3 0.2 etc.
//4 0.0
//5 0.2
//6 0.4
//$EndNodeData

//2 two elements:
//1 3 2 99 2 1 2 3 4 quad #1: type 3, physical 99, elementary 2, nodes 1 2 3 4
//2 3 2 99 2 2 5 6 3 quad #2: type 3, physical 99, elementary 2, nodes 2 5 6 3

//elm-number
//is the number (index) of the n-th element in the mesh; elm-number must be a
//postive (non-zero) integer. Note that the elm-numbers do not necessarily have
//to form a dense nor an ordered sequence.
//elm-type
//number-of-tags
//gives the number of integer tags that follow for the n-th element. By default,
//the first tag is the number of the physical entity to which the element belongs;
//the second is the number of the elementary geometrical entity to which the
//element belongs; the third is the number of mesh partitions to which the element
//belongs, followed by the partition ids (negative partition ids indicate ghost
//cells). A zero tag is equivalent to no tag.

//number-of-elements
//elm-number elm-type number-of-tags < tag > ... node-number-list

//quad #1: type 3, physical 99, elementary 2, nodes 1 2 3 4
//#1 : 3-node triangle.type 
//elm-number elm-type number-of-tags physical    elementary 
//1           2          2            0          2
//1 2 2 0 2 1 2 3 
//2 2 2 0 2 2 4 3 
//3 2 2 0 2 5 6 7  
//
//35 27 2 1 0 25 17 26 27 28
//36 27 2 1 0 17 18 29 30 31
//37 27 2 1 0 18 19 32 33 34
//38 27 2 1 0 19 11 35 36 37
//39 27 2 1 0 11 12 38 39 40
//40 27 2 1 0 12 23 41 42 43 

void Gmesh::ReadElements(fstream &file, CGNSFactory *factory_cgns)
{
    string line,word;
    string separator = " =\t\r\n#$,;";
    do
    {
        getline(file,line);
        line = FindNextWord(line,word,separator);
        transform(word.begin(), word.end(), word.begin(), ToLower());
    } while (word != "elements");

    getline(file,line);

    //所有单元，包括边界单元
    int elements = 0;

    line = FindNextWord(line,word,separator);
    from_string< int >(elements,word,std::dec);

    cout << " Read grid coordinates \n";

    int elem_index,gmsh_elem_type;
    int number_of_tags;
    int elem_tag[3];

    vector<int> elem_type_container;
    vector<int> elem_number_container;
    vector<int> elem_index_container;
    vector<int> elem_bctype_container;

    int *node_gmsh = new int [NofValidElementTypesGmsh]();
    int *node_cgns = new int [NofValidElementTypesGmsh]();

    GmeshElem *gmshelem = new GmeshElem();

    for (int i = 0; i < elements; ++ i)
    {
        getline(file,line);
        line = FindNextWord(line,word,separator);
        from_string< int >(elem_index,word,std::dec);

        line = FindNextWord(line,word,separator);
        from_string< int >(gmsh_elem_type,word,std::dec);

        line = FindNextWord(line,word,separator);
        from_string< int >(number_of_tags,word,std::dec);

        for (int j = 0; j < number_of_tags; ++ j)
        {
            line = FindNextWord(line,word,separator);
            from_string< int >(elem_tag[j],word,std::dec);
        }

        elem_bctype_container.push_back(elem_tag[0]);

        int cgns_elem_type = TypeGmsh2CGNS[ gmsh_elem_type ];
        int cgns_vnum = cgns_nvertex[ cgns_elem_type ];
        int gesh_vnum = gmesh_nvertex[ gmsh_elem_type ];

        elem_type_container.push_back(cgns_elem_type);
        elem_number_container.push_back(1);

        int nodeindex;
        for (int j = 0; j < gesh_vnum; ++ j)
        {
            line = FindNextWord(line,word,separator);
            from_string< int >(nodeindex,word,std::dec);
            node_gmsh[j] = nodeindex;
        }

        ElemReorderGmsh2CGNS(cgns_elem_type, node_gmsh, node_cgns, gesh_vnum);
        for (int j = 0; j < cgns_vnum; ++ j)
        {
            elem_index_container.push_back(node_cgns[j]);
        }
    }

    gmshelem->TreatElements(factory_cgns, elem_type_container, elem_number_container, elem_index_container, elem_bctype_container, cgns_nvertex);

    delete gmshelem;    gmshelem = NULL;
    delete [] node_gmsh;    node_gmsh = NULL;
    delete [] node_cgns;    node_cgns = NULL;
}

}
