#include "GridAdaptationCGNS.h"
#include "Glb_Dimension.h"
#include "HyList.h"
#include "TK_Exit.h"
#include "Pre_GridConversion.h"
#include "Geo_StructGrid.h"
#include "Geo_Element.h"
using namespace std;

#pragma warning (disable:913)
namespace PHSPACE
{

int del(int x, int y);

CGNSection::CGNSection(CGNSZone *zone)
{
    this->zone = zone;
    connection_list = 0;
};

CGNSection::~CGNSection()
{
    delete [] connection_list;
};

void CGNSection::CreateConnectionList(int element_node_number)
{
    cgsize_t element_number = iend - istart + 1;
    connection_list = new cgsize_t [ element_number * element_node_number ];
}

cgsize_t * CGNSection::GetElementIndex(cgsize_t ielement) const
{
    int element_node_number  = GetElementNodeNumbers(element_type);
    return &connection_list[(ielement-istart) * element_node_number];
}

void CGNSection::GetElementIndex(int ielement, vector<cgsize_t> &element_index)
{
    int element_node_number  = GetElementNodeNumbers(element_type);
    cgsize_t *element_node        = GetElementIndex(ielement);
    element_index.resize(0);
    for (int i = 0; i < element_node_number; ++ i)
    {
        element_index.push_back(element_node[i]);
    }
}

void CGNSection::SetElementTypeAndIndex(vector<int> &element_type, vector< vector<int> > &element_index)
{
    cgsize_t istart = GetStart();
    cgsize_t iend   = GetEnd();

    int section_element_type = GetElementType();
    int element_node_number  = GetElementNodeNumbers(section_element_type);

    for (cgsize_t i = istart; i <= iend; ++ i)
    {
        element_type.push_back(section_element_type);
        cgsize_t * elem = GetElementIndex(i);

        vector<int> elem_index(elem, elem + element_node_number);
        element_index.push_back(elem_index);
    }
}

void CGNSection::Convert2InerDataStandard()
{
    istart -= 1;
    iend   -= 1;

    cgsize_t element_number = iend - istart + 1;
    int element_node_number  = GetElementNodeNumbers(element_type);

    cgsize_t connection_size = element_number * element_node_number;
        
    for (int i = 0; i < connection_size; ++ i)
    {
        connection_list[i] -= 1;
    }
}

void CGNSection::ReadSection()
{
    int fileid = zone->GetFileIndex();
    int baseid = zone->GetBaseIndex();
    int zoneid = zone->GetZoneID();
    int sectionid = GetSectionIndex();

    cgsize_t istart, iend;
    int nbndry, iparent_flag;

    ElementType_t etype;

    cg_section_read(fileid, baseid, zoneid, sectionid, section_name, &etype, &istart, &iend, &nbndry, &iparent_flag);

    cout << "\nReading section data...\n";
    cout << "   section name = " << section_name << "\n";
    cout << "   section type = " << ElementTypeName[etype] << "\n";
    cout << "   istart,iend =  " << istart << " " << iend << "\n";

    SetStart(istart);
    SetEnd  (iend);
    SetElementType(etype);

    int element_node_number = GetElementNodeNumbers(etype);

    CreateConnectionList(element_node_number);
    cgsize_t *conn_list = GetConnectionList();

    // Read the connectivity. Again, the node numbering of the 
    // connectivities start at 1. If internally a starting index 
    // of 0 is used (typical for C-codes) 1 must be substracted 
    // from the connectivities read. 
    cg_elements_read(fileid, baseid, zoneid, sectionid, conn_list, NULL);
}

CGNSectionProxy::CGNSectionProxy(CGNSZone *zone)
{
    this->zone = zone;
    nsections = 0;
    sections  = 0;
};

CGNSectionProxy::~CGNSectionProxy()
{
    for (int i = 0; i < nsections; ++ i)
    {
        delete sections[i];
    }
    delete [] sections;
};

void CGNSectionProxy::Create(int nsections)
{
    this->nsections = nsections;
    sections = new CGNSection *[nsections];

    for (int isection = 0; isection < nsections; ++ isection)
    {
        sections[isection] = new CGNSection(zone);
        sections[isection]->SetSectionIndex(isection+1);
    }
};

CGNSection * CGNSectionProxy::GetSectionByElementIndex(int ielement) const
{
    for (int i = 0; i < nsections; ++ i)
    {
        CGNSection *section = sections[i];
        if (section->GetStart() <= ielement && ielement <= section->GetEnd())
        {
            return section;
        }
    }
    return 0;
}

void CGNSectionProxy::Convert2InerDataStandard()
{
    int nsection = GetNSection();
    for (int isection = 0; isection < nsection; ++ isection)
    {
        CGNSection * section = GetSection(isection);
        section->Convert2InerDataStandard();
    }
}

int del(int x, int y)
{
    if (ABS(x) == ABS(y))
    {
        return 1;
    }
    return 0;
}

CGNSBCRegion::CGNSBCRegion(CGNSZone *zone)
{
    this->zone = zone;
    bc_connection_list = 0;
};

CGNSBCRegion::~CGNSBCRegion()
{
    delete [] bc_connection_list;
};

void CGNSBCRegion::CreateConnectionList()
{
    if (zone->GetZoneType() == Unstructured)
    {
        bc_connection_list = new cgsize_t [ bc_element_number ];
    }
    else
    {
        if (bc_element_number == 2)
        {
            bc_element_number = 6;
        }
        bc_connection_list = new cgsize_t [ bc_element_number ];
    }
}

void CGNSBCRegion::PrintConnectionList()
{
    if (zone->GetZoneType() == Unstructured)
    {
        if (bc_element_number == 2)
        {
            cout << " " << bc_connection_list[0] << " " << bc_connection_list[1] << "\n";
        }
        else
        {
            for (int i = 0; i < MIN(20,bc_element_number); ++ i)
            {
                cout << " " << bc_connection_list[i];
            }
            cout << "\n";
        }

    }
    else
    {
        int size = 6;
        int dim = zone->GetCellDimension();
        if (dim == 2) size = 4;

        for (int i = 0; i < size; ++ i)
        {
            cout << " " << bc_connection_list[i];
        }
        cout << "\n";
    }
}

void CGNSBCRegion::SetStructuredGridBoundary(StructBC * bcregion)
{
    cgsize_t imin,imax,jmin,jmax,kmin,kmax;
    int dim = zone->GetCellDimension();
    if (dim == 2)
    {
        imin = bc_connection_list[0];
        jmin = bc_connection_list[1];
        kmin = 1;
        imax = bc_connection_list[2];
        jmax = bc_connection_list[3];
        kmax = 1;
    }
    else
    {
        imin = bc_connection_list[0];
        jmin = bc_connection_list[1];
        kmin = bc_connection_list[2];
        imax = bc_connection_list[3];
        jmax = bc_connection_list[4];
        kmax = bc_connection_list[5];
    }


    BaseBCType *bc_convert = new BaseBCType();

    bcregion->SetIJKRegion(static_cast<int>(imin),static_cast<int>(imax),static_cast<int>(jmin),static_cast<int>(jmax),static_cast<int>(kmin),static_cast<int>(kmax));
    bcregion->SetBCType(bc_convert->GetBCType(bc_type));

    delete bc_convert;
}

void CGNSBCRegion::ReadOrdinaryBoundaryRegion(int ibcregion)
{
    cgsize_t npts;
    int normal_index[3], normal_list, ndataset;
    cgsize_t norm_list_flag;
    DataType_t     norm_data_type;
    BCType_t       bctype;
    PointSetType_t iptset;
    GridLocation_t igrid_location;

    int fileid = zone->GetFileIndex();
    int baseid = zone->GetBaseIndex();
    int zoneid = zone->GetZoneID();

    cg_goto(fileid, baseid, "Zone_t", zoneid, "ZoneBC_t", 1, "BC_t", ibcregion, "end");
    cg_gridlocation_read(&igrid_location);

    cout << " Grid Location Name = " << GridLocationName[igrid_location] << "\n";
    if (igrid_location == FaceCenter)
    {
        cout << "\nGridLocation = FaceCenter means BC data refers to elements, not nodes\n";
    }
    else if (igrid_location == Vertex)
    {
        cout << "\nGridLocation = Vertex means BC data refers to nodes, not elements\n";

        //cout << " Now Change Vertex to FaceCenter for incompatible reasons\n";
        //igrid_location = FaceCenter;
    }
    // Read the info for this boundary condition.
    cg_boco_info(fileid, baseid, zoneid, ibcregion, bcregion_name, &bctype, &iptset, &npts, normal_index,
                &norm_list_flag, &norm_data_type, &ndataset);


    cout << " Boundary Name = " << bcregion_name << "\n";
    cout << " Boundary Condition Name = " << BCTypeName[bctype] << "\n";
    cout << " Boundary Element Number or Boundary Point Number = " << npts << "\n";

    SetBCElementNumber(static_cast<int>(npts));
    SetBCType         (bctype       );
    SetBCPointSetType (iptset       );
    SetBCGridLocation (igrid_location);
    CreateConnectionList();
    SetNTotalCell(zone->GetNTotalCell());

    cgsize_t *bc_conn_list = GetConnectionList();

    cout << " PointSetType_t = " << iptset << "\n";

    // Read the element ID’s.
    cg_boco_read(fileid, baseid, zoneid, ibcregion, bc_conn_list, &normal_list);
    PrintConnectionList();
}

void CGNSBCRegion::Read1to1BoundaryRegion(int i1to1)
{
    cgsize_t ipntsdonor[6];
    int itranfrm[3];

    CGNSTraits::char_33 donorname;

    int fileid = zone->GetFileIndex();
    int baseid = zone->GetBaseIndex();
    int zoneid = zone->GetZoneID();

    //Zone Connectivity
    cg_goto(fileid, baseid, "Zone_t", zoneid, "ZoneGridConnectivity_t", 1, "GridConnectivity1to1_t", 1, "end");

    SetBCElementNumber(6               );
    SetBCType         (BCTypeUserDefined);
    SetBCPointSetType (PointRange      );
    SetBCGridLocation (FaceCenter      );
    CreateConnectionList();
    SetNTotalCell(zone->GetNTotalCell());

    cgsize_t *ipnts = GetConnectionList();
    cg_1to1_read(fileid, baseid, zoneid, i1to1, bcregion_name, donorname, ipnts, ipntsdonor, itranfrm);

    int transform[3][3];

    //For 3-D, the
    //transformation matrix T is constructed from Transform = [±a, ±b, ±c] as follows:
    //T =
    //[sgn(a) del(a ? 1) sgn(b) del(b ? 1) sgn(c) del(c ? 1)]
    //[sgn(a) del(a ? 2) sgn(b) del(b ? 2) sgn(c) del(c ? 2)]
    //[sgn(a) del(a ? 3) sgn(b) del(b ? 3) sgn(c) del(c ? 3)]

    for (int i = 0; i < zone->GetCellDimension(); ++ i)
    {
        for (int j = 0; j < zone->GetCellDimension(); ++ j)
        {
            transform[i][j] = SIGN(1,itranfrm[j]) * del(itranfrm[j], i+1);
        }
    }
}

int CGNSBCRegion::DetermineBase()
{
    for (int i = 0; i < bc_element_number; ++ i)
    {
        if (bc_connection_list[i] < nTotalCell)
        {
            return 0;
        };
    }
    return 1;
}

void CGNSBCRegion::Convert2InerDataStandard()
{
    //对于各种情况均成立
    for (int i = 0; i < bc_element_number; ++ i)
    {
        bc_connection_list[i] -= 1;
    }
}

void CGNSBCRegion::ShiftBCRegion()
{
    if (bc_grid_location != Vertex)
    {
        for (int i = 0; i < bc_element_number; ++ i)
        {
            bc_connection_list[i] += nTotalCell;
        }
    }
}

void CGNSBCRegion::ResolveElementBoundaryFace(ElementFaceProxy *element_face_proxy)
{
    set<cgsize_t> bc_vertex;
    if (bc_grid_location == Vertex)
    {
        ProcessBCVertex(bc_vertex);
    }
    else
    {
        ProcessBCFace (bc_vertex);
    }
    element_face_proxy->ResolveElementBoundaryFace(bc_vertex, bc_type);
}

void CGNSBCRegion::ProcessBCVertex(set<cgsize_t> &bc_vertex)
{
    if (bc_point_set_type == ElementRange || bc_point_set_type == PointRange)
    {
        for (cgsize_t ibc = bc_connection_list[0]; ibc <= bc_connection_list[1]; ++ ibc)
        {
            bc_vertex.insert(ibc);
        }
    }
    else
    {
        for (cgsize_t ibc = 0; ibc < bc_element_number; ++ ibc)
        {
            bc_vertex.insert(bc_connection_list[ibc]);
        }
    }
}

void CGNSBCRegion::ProcessBCFace(set<cgsize_t> &bc_vertex)
{
    vector<cgsize_t> face_node;
    if (bc_point_set_type == ElementRange || bc_point_set_type == PointRange)
    {
        for (cgsize_t ibc = bc_connection_list[0]; ibc <= bc_connection_list[1]; ++ ibc)
        {
            zone->GetElementIndex(static_cast<int>(ibc), face_node);
            for (std::size_t inode = 0; inode < face_node.size(); ++ inode)
            {
                bc_vertex.insert(face_node[inode]);
            }
        }
    }
    else
    {
        for (cgsize_t ibc = 0; ibc < bc_element_number; ++ ibc)
        {
            zone->GetElementIndex(static_cast<int>(bc_connection_list[ibc]), face_node);
            for (std::size_t inode = 0; inode < face_node.size(); ++ inode)
            {
                bc_vertex.insert(face_node[inode]);
            }
        }
    }
}

CGNSBCRegionProxy::CGNSBCRegionProxy(CGNSZone *zone)
{
    this->zone           = zone;
    nbcregions           = 0;
    bcregions            = 0;
    n1to1                = 0;
    n_ordinary_bcregions = 0;
};

CGNSBCRegionProxy::~CGNSBCRegionProxy()
{
    for (int ibcregion = 0; ibcregion < nbcregions; ++ ibcregion)
    {
        delete bcregions[ibcregion];
    }
    delete [] bcregions;
};

void CGNSBCRegionProxy::Create(int n_ordinary_bcregions, int n1to1)
{
    this->n_ordinary_bcregions = n_ordinary_bcregions;
    this->n1to1                = n1to1;
    this->nbcregions = n_ordinary_bcregions + n1to1;
    bcregions = new CGNSBCRegion *[nbcregions];

    for (int ibcregion = 0; ibcregion < nbcregions; ++ ibcregion)
    {
        bcregions[ibcregion] = new CGNSBCRegion(zone);
    }
};

void CGNSBCRegionProxy::Convert2InerDataStandard()
{
    int nbcregions = GetNBCRegion();
    for (int ibcregion = 0; ibcregion < nbcregions; ++ ibcregion)
    {
        CGNSBCRegion * bcregion = GetBCRegion(ibcregion);
        bcregion->Convert2InerDataStandard();
    }

    int bcregion_base = 1;
    for (int ibcregion = 0; ibcregion < nbcregions; ++ ibcregion)
    {
        CGNSBCRegion * bcregion = GetBCRegion(ibcregion);
        if (!bcregion->DetermineBase())
        {
            bcregion_base = 0;
            break;
        }
    }

    if (bcregion_base == 0)
    {
        for (int ibcregion = 0; ibcregion < nbcregions; ++ ibcregion)
        {
            CGNSBCRegion * bcregion = GetBCRegion(ibcregion);
            bcregion->ShiftBCRegion();
        }
    }
}

CGNSZone::CGNSZone()
{
    grid           = 0;
    section_proxy  = 0;
    bcregion_proxy = 0;
};

CGNSZone::~CGNSZone()
{
    delete grid;
    delete section_proxy;
    delete bcregion_proxy;
};

void CGNSZone::Create()
{
    grid           = new SimpleGrid();
    section_proxy  = new CGNSectionProxy(this);
    bcregion_proxy = new CGNSBCRegionProxy(this);
};

void CGNSZone::CreateSections(int nsections)
{
    section_proxy->Create(nsections);
};

void CGNSZone::CreateBCRegions(int n_ordinary_bcregions, int n1to1)
{
    bcregion_proxy->Create(n_ordinary_bcregions, n1to1); 
};

void CGNSZone::SetDimension()
{
    if (zone_type == Structured)
    {
        //结构网格赋值

        // lower range index
        irmin[0] = 1;
        irmin[1] = 1;
        irmin[2] = 1;

        // upper range index of vertices
        irmax[0] = 1;
        irmax[1] = 1;
        irmax[2] = 1;

        cell_size[0] = 1;
        cell_size[1] = 1;
        cell_size[2] = 1;

        // upper range index of vertices
        // vertex size
        if (cell_dimension == 2)
        {
            irmax[0] = isize[0][0];
            irmax[1] = isize[0][1];

            // cell size
            cell_size[0] = static_cast<int>(isize[0][2]);
            cell_size[1] = static_cast<int>(isize[1][0]);
        }
        else
        {
            irmax[0] = isize[0][0];
            irmax[1] = isize[0][1];
            irmax[2] = isize[0][2];

            // cell size
            cell_size[0] = static_cast<int>(isize[1][0]);
            cell_size[1] = static_cast<int>(isize[1][1]);
            cell_size[2] = static_cast<int>(isize[1][2]);
        }

        nTotalNode = static_cast<int>(irmax[0] * irmax[1] * irmax[2]);
        nTotalCell = cell_size[0] * cell_size[1] * cell_size[2];
    }
    else
    {
        irmin[0] = 1;
        irmin[1] = 0;
        irmin[2] = 0;

        irmax[0] = isize[0][0];
        irmax[1] = 0;
        irmax[2] = 0;

        //注意，下面和一般的isize[1][0] 有所区别，因为这里我们定义isize[3][3]
        //通常对于非结构定义isize[3][1]，我们这样定义是为了统一处理

        cell_size[0] = static_cast<int>(isize[0][1]);

        nTotalNode = static_cast<int>(irmax[0]);
        nTotalCell = cell_size[0];
    }

    cout << " nTotalNode = " << nTotalNode << " nTotalCell = " << nTotalCell << "\n";
}

void CGNSZone::CopyStructuredGrid(Grid *grid_in)
{
    StructGrid *grid = StructGridCast (grid_in);

    SimpleGrid * simple_grid = GetSimpleGrid();

    cgsize_t * nijk = GetIMax();
    grid->CopyXYZ(static_cast<int>(nijk[0]), static_cast<int>(nijk[1]), static_cast<int>(nijk[2]), simple_grid->GetX(), simple_grid->GetY(), simple_grid->GetZ());

    int nBCRegion = bcregion_proxy->GetNBCRegion();

    grid->CreateCompositeBCRegion(nBCRegion);
    CompositeStructBC * composite_bcregion = grid->GetCompositeBCRegion();

    for (int ibcregion = 0; ibcregion < nBCRegion; ++ ibcregion)
    {
        StructBC * bcregion = new StructBC(GetZoneID()-1, ibcregion);
        composite_bcregion->SetBCRegion(ibcregion, bcregion);

        CGNSBCRegion *cgns_bcregion = GetBCRegion(ibcregion);
        cgns_bcregion->SetStructuredGridBoundary(bcregion);
    }

    grid->ProcessBCInfo();
}

void CGNSZone::ReadGrid()
{
    //读取网格坐标
    ReadGridCoordinates();

    //读取网格单元连接关系
    ReadGridConnectivities();

    //读取边界条件
    ReadGridBoundaryConditions();
}

void CGNSZone::ReadGridCoordinates()
{
    //Check the zone type
    cg_zone_type(fileid, baseid, zoneid, &zone_type);

    //Determine the number of vertices and volume elements in this zone
    cg_zone_read(fileid, baseid, zoneid, zonename, isize[0]);

    SetDimension();

    ReadCoor();

    //ShowStructuredGrid(simple_grid, irmax[0], irmax[1], irmax[2]);
}

void CGNSZone::ReadGridBoundaryConditions()
{
    // find out how many general interfaces there are in this zone
    int n1to1, nOrdinaryBCRegions;
    cg_n1to1(fileid, baseid, zoneid, &n1to1);

    // Determine the number of boundary conditions for this zone.
    cg_nbocos(fileid, baseid, zoneid, &nOrdinaryBCRegions); 
    CreateBCRegions(nOrdinaryBCRegions, n1to1);

    // Loop over the number of boundary conditions.
    for (int ibcregion = 1; ibcregion <= nOrdinaryBCRegions; ++ ibcregion)
    {
        CGNSBCRegion *cgns_bcregion = GetBCRegion(ibcregion-1);
        cgns_bcregion->ReadOrdinaryBoundaryRegion(ibcregion);
    }

    for (int i1to1 = 1; i1to1 <= n1to1; ++ i1to1)
    {
        CGNSBCRegion *cgns_bcregion = Get1to1BCRegion(i1to1-1);
        cgns_bcregion->Read1to1BoundaryRegion(i1to1);
    }
}

void CGNSZone::ReadGridConnectivities()
{
    if (GetZoneType() == Structured) return;

    // Determine the number of sections for this zone. Note that
    // surface elements can be stored in a volume zone, but they
    // are NOT taken into account in the number obtained from 
    // cg_zone_read.
    int nSections;
    cg_nsections(fileid, baseid, zoneid, &nSections);

    cout << " nSections = " << nSections << "\n";

    CreateSections(nSections);

    // Loop over the number of sections and read the element
    // connectivities. As CGNS starts the numbering at 1 the
    // for-loop starts at 1 as well.

    for (int isection = 1; isection <= nSections; ++ isection)
    {
        CGNSection *cgns_section = GetSection(isection-1);
        cgns_section->ReadSection();
    }
}

void CGNSZone::ShowStructuredGrid(SimpleGrid *simple_grid, int ni, int nj, int nk)
{
    GridID *index = new GridID(0);
    StructGrid * grid = StructGridCast(CreateGridGeneral(STRUCTGRID, index, 0, GetDim()));
    grid->CopyXYZ(ni, nj, nk, simple_grid->GetX(), simple_grid->GetY(), simple_grid->GetZ());

    grid->ToTecplot();
    delete grid;
}

void CGNSZone::ReadCoor()
{
    int nCoords;
    //Determine the number and names of the coordinates.
    cg_ncoords(fileid, baseid, zoneid, &nCoords);

    int nTotalNode = GetNTotalNode();

    //暂时假定grid为非结构，以后再说
    RDouble *x = new RDouble[ nTotalNode ];
    RDouble *y = new RDouble[ nTotalNode ];
    RDouble *z = new RDouble[ nTotalNode ];

    SimpleGrid * simple_grid = GetSimpleGrid();

    simple_grid->SetNTotalNode(nTotalNode);
    simple_grid->SetX(x);
    simple_grid->SetY(y);
    simple_grid->SetZ(z);

    CGNSRawCoor * raw_coor = new CGNSRawCoor();

    for (int icoor = 0; icoor < nCoords; ++ icoor)
    {
        DataType_t dataType;
        cg_coord_info(fileid, baseid, zoneid, icoor + 1, &dataType, coordname);
        raw_coor->AllocateData(icoor, nTotalNode, dataType);

        //Read the x-, y-, z-coordinates.
        cg_coord_read(fileid, baseid, zoneid, coordname, dataType, GetIMin(), GetIMax(), raw_coor->GetCoor(icoor));
    }

    raw_coor->SetAllData(x, y, z);
    delete raw_coor;
}

void CGNSZone::SetElementTypeAndIndex(vector<int> &element_type, vector< vector<int> > &element_index)
{
    int nsection = GetNSection();
    for (int isection = 0; isection < nsection; ++ isection)
    {
        CGNSection * section = GetSection(isection);
        section->SetElementTypeAndIndex(element_type, element_index);
    }
}

void CGNSZone::Convert2InerDataStandard()
{
    this->section_proxy->Convert2InerDataStandard();
    this->bcregion_proxy->Convert2InerDataStandard();
}

void CGNSZone::BuildPoints(PointOperation *point_operation)
{
    int nTotalNode = grid->GetNTotalNode();
    RDouble *x = grid->GetX();
    RDouble *y = grid->GetY();
    RDouble *z = grid->GetZ();

    for (int inode = 0; inode < nTotalNode; ++ inode)
    {
        PointOperation::point_type point(x[inode],y[inode],z[inode],inode);
        point_operation->AddPoint(point);
    }
}

void CGNSZone::ResolveElementBoundaryFace(ElementFaceProxy *element_face_proxy)
{
    int nbcregions = GetNBCRegion();
    for (int ibcregion = 0; ibcregion < nbcregions; ++ ibcregion)
    {
        CGNSBCRegion *bcregion = GetBCRegion(ibcregion);

        bcregion->ResolveElementBoundaryFace(element_face_proxy);
    }
}

void CGNSZone::GetElementIndex(int ielement, vector<cgsize_t> &element_index)
{
    CGNSection * cgns_section = section_proxy->GetSectionByElementIndex(ielement);
    cgns_section->GetElementIndex(ielement, element_index);
}

CGNSMultiZone::CGNSMultiZone()
{
    nzones = 0;
    zones  = 0;
}

CGNSMultiZone::~CGNSMultiZone()
{
    for (int i = 0; i < nzones; ++ i)
    {
        delete zones[i];
    }
    delete [] zones;
}

void CGNSMultiZone::InitAllZones()
{
    this->nzones = PHMPI::GetNumberofGlobalZones();
    zones = new CGNSZone *[nzones];
    for (int izone = 0; izone < nzones; ++ izone)
    {
        zones[izone] = new CGNSZone();
        zones[izone]->Create();
        zones[izone]->SetCellDimension(cellDim);
        zones[izone]->SetPhysDimension(physDim);
        zones[izone]->SetFileIndex(fileid);
        zones[izone]->SetBaseIndex(baseid);
        zones[izone]->SetZoneID(izone+1); //zoneindex start from 1
    }
}

void CGNSMultiZone::ReadGrid(const string &filename)
{
    //Open the CGNS for reading and check if the file was found.
    if (cg_open(filename.c_str(), CG_MODE_READ, &fileid) != CG_OK)
    {
        ostringstream oss;
        oss << cg_get_error() << "\n";
        TK_Exit::ExceptionExit(oss.str());
    }

    //Determine the of bases in the grid
    cg_nbases(fileid, &baseid);

    //Read the number of zones in the grid.
    cg_nzones(fileid, baseid, &nzones);

    //Check the cell and physical dimensions of the bases.
    cg_base_read(fileid, baseid, basename, &cellDim, &physDim);

    InitAllZones();

    for (int izone = 1; izone <= nzones; ++ izone)
    {
        CGNSZone *cgns_zone = GetZone(izone-1);
        cgns_zone->ReadGrid();
    }

    //close CGNS file
    cg_close(fileid);
}

PointOperation::PointOperation()
{
}

PointOperation::~PointOperation()
{
}

ElementFaceProxy::ElementFaceProxy(PointOperation *point_operation)
{
    this->point_operation = point_operation;
    reference_faces    = new set< MultiIndex >;
    face_to_nodes      = new vector< vector<int> >;
    face_to_cell_left  = new vector< int >;
    face_to_cell_right = new vector< int >;
    face_boundary_type = new vector< int >;
    face_element_type  = new vector< int >;
    face_parent        = new vector< int >;
    face_child         = new vector< vector<int> >;

    boundary_face_bctype       = new vector< int >;
    boundary_face_element_type = new vector< int >;
    boundary_face_to_nodes     = new vector< vector<int> >;
    nTotalCell = 0;
}

ElementFaceProxy::~ElementFaceProxy()
{
    delete reference_faces;
    delete face_to_nodes;
    delete face_to_cell_left;
    delete face_to_cell_right;
    delete face_boundary_type;
    delete face_element_type;
    delete face_parent;
    delete face_child;

    delete boundary_face_bctype;
    delete boundary_face_element_type;
    delete boundary_face_to_nodes;
}

void ElementFaceProxy::FreeMemory()
{
    delete reference_faces;
    delete face_to_nodes;
    delete face_to_cell_left;
    delete face_to_cell_right;
    delete face_boundary_type;
    delete face_element_type;
    delete face_parent;
    delete face_child;

    reference_faces    = new set< MultiIndex >;
    face_to_nodes      = new vector< vector<int> >;
    face_to_cell_left  = new vector< int >;
    face_to_cell_right = new vector< int >;
    face_boundary_type = new vector< int >;
    face_element_type  = new vector< int >;
}

int ElementFaceProxy::FindFace(MultiIndex &face)
{
    set< MultiIndex >::iterator iter = reference_faces->find(face);
    if (iter == reference_faces->end())
    {
        return -1;
    }
    return iter->id;
}

int ElementFaceProxy::FindFace(vector<int> &facenodes)
{
    int face_point_number = static_cast<int>(facenodes.size());

    vector<int> face_of_sort = facenodes;
    sort(face_of_sort.begin(), face_of_sort.end());

    MultiIndex face_of_multi_index(face_point_number, GetNTotalFace());
    face_of_multi_index.SetData(&face_of_sort[0]);

    return FindFace(face_of_multi_index);
}

void ElementFaceProxy::InsertFace(MultiIndex &face)
{
    reference_faces->insert(face);
}

void ElementFaceProxy::ResizeFaceNumber(int n)
{
    face_to_nodes->resize(n);
    face_to_cell_left->resize(n);
    face_to_cell_right->resize(n);
    face_boundary_type->resize(n);
    face_element_type->resize(n);

    face_parent->resize(n);
    face_child->resize(n);
}

void ElementFaceProxy::ResolveElementFace(vector<int> &element_index, int element_type, int ielement)
{
    BasicElement *basic_element = GetBasicElement(element_type);
    uint_t face_number = basic_element->GetElementFaceNumber();

    vector<int> element_face;
    //vector<int> sort_face;

    for (int iface = 0; iface < face_number; ++ iface)
    {
        vector<int> & face_index = basic_element->GetElementFace(iface);
        int face_type = basic_element->GetFaceElementType(iface);
         
        int face_point_number = static_cast<int>(face_index.size());

        element_face.resize(0);
        for (int i = 0; i < face_point_number; ++ i)
        {
            element_face.push_back(element_index[ face_index[i] ]);
        }

        vector<int> face_of_sort = element_face;
        sort(face_of_sort.begin(), face_of_sort.end());

        MultiIndex face_of_multi_index(face_point_number, GetNTotalFace());
        face_of_multi_index.SetData(&face_of_sort[0]);

        int id = FindFace(face_of_multi_index);
        if (id == -1)
        {
            InsertFace(face_of_multi_index);
            face_to_cell_left ->push_back(ielement);
            face_to_cell_right->push_back(-1);
            face_boundary_type->push_back(0);
            face_element_type->push_back(face_type);
            face_to_nodes->push_back(element_face);

            face_parent->push_back(-1);
            face_child->resize(face_child->size()+1);
        }
        else
        {
            //face已经存在
            //如果此face has child
            int f2cr = (*face_to_cell_right)[id];

            bool isright = (f2cr == ielement);
            if ((*face_child)[id].size() != 0)
            {
                //遍历所有child face
                for (std::size_t ichild = 0; ichild < (*face_child)[id].size(); ++ ichild)
                {
                    int child_face_id = (*face_child)[id][ichild];
                    if (isright)
                    {
                        (*face_to_cell_right)[child_face_id] = (*face_to_cell_right)[id];
                    }
                    else
                    {
                        (*face_to_cell_right)[child_face_id] = ielement;
                    }
                }
            }
            else
            {
                if ((*face_to_cell_left)[id] == -1)
                {
                    //说明此面虽然存在，但是由于种种原因没有处理
                    (*face_to_cell_left)[id] = ielement;
                }
                else
                {
                    if ((*face_to_cell_right)[id] == -1)
                    {
                        //此面无parent
                        if ((*face_parent)[id] == -1 && (*face_to_cell_left)[id] != ielement)
                        {
                            (*face_to_cell_right)[id] = ielement;
                        }
                    }
                }
            }
        }
    }
}

void ElementFaceProxy::AddFaceHierachical(int parent_face_id, int child_face_id)
{
    (*face_parent)[child_face_id] = parent_face_id;
    (*face_child)[parent_face_id].push_back(child_face_id);
}

void ElementFaceProxy::SetDefaultFaceInfo(int index, vector<int> &facenodes)
{
    (*face_to_cell_left)[index] = -1;
    (*face_to_cell_right)[index] = -1;
    (*face_boundary_type)[index] =  0;
    (*face_element_type)[index] = -1;
    (*face_to_nodes    )[index] = facenodes;
}

int ElementFaceProxy::AddElementFace(vector<int> &facenodes, int parent_face_id)
{
    int face_point_number = static_cast<int>(facenodes.size());

    vector<int> face_of_sort = facenodes;
    sort(face_of_sort.begin(), face_of_sort.end());

    MultiIndex face_of_multi_index(face_point_number, GetNTotalFace());
    face_of_multi_index.SetData(&face_of_sort[0]);

    int child_face_id = FindFace(face_of_multi_index);
    if (child_face_id == -1) 
    {
        InsertFace(face_of_multi_index);
        child_face_id = GetNTotalFace();
        ResizeFaceNumber(child_face_id + 1);
        SetDefaultFaceInfo(child_face_id, facenodes);
        AddFaceHierachical(parent_face_id, child_face_id);
    }
    return child_face_id;
}

bool ElementFaceProxy::CheckBCFace(set<cgsize_t> &bc_vertex, vector<int> &face_node)
{
    uint_t size = face_node.size();
    for (int i = 0; i < size; ++ i)
    {
        set<cgsize_t>::iterator iter = bc_vertex.find(face_node[i]);
        if (iter == bc_vertex.end())
        {
            return false;
        }
    }
    return true;
}

void ElementFaceProxy::ResolveElementBoundaryFace(set<cgsize_t> &bc_vertex, int bctype)
{
    int nTotalFace = this->GetNTotalFace();
    for (int iface = 0; iface < nTotalFace; ++ iface)
    {
        int original_bc_type = (*face_boundary_type)[iface];
        int right_cell_of_face = (*face_to_cell_right)[iface];

        if (right_cell_of_face == - 1 && original_bc_type == 0)
        {
            if (CheckBCFace(bc_vertex, (*face_to_nodes)[iface]))
            { 
                (*face_boundary_type)[iface] = bctype;

                boundary_face_bctype->push_back(bctype);
                boundary_face_element_type->push_back((*face_element_type)[iface]);
                boundary_face_to_nodes->push_back((*face_to_nodes)[iface]);
            }
        }
    }

    point_operation->SetBoundaryPoint(bc_vertex, bctype);

}

void ElementFaceProxy::ResolveElementBoundaryFace(int iadapt)
{
    //for (int i = 0; i < iadapt; ++ i)
    //{
    //    SplitElementBoundaryFace();
    //}
    MarkElementBoundaryFace();
}

void ElementFaceProxy::MarkElementBoundaryFace()
{
    uint_t size = boundary_face_bctype->size();
    for (int ibcface = 0; ibcface < size; ++ ibcface)
    {
        int bctype    = (*boundary_face_bctype)[ibcface];
        int face_type = (*boundary_face_element_type)[ibcface];

        //检测单元是否是基本单元
        if (!ISBasicElement(face_type)) continue;

        vector<int> face_node = (*boundary_face_to_nodes)[ibcface];

        SetFaceBoundary(face_node, face_type, bctype);
    }
}

int ElementFaceProxy::GetRootFaceID(int iface)
{
    int parent_face_id = (*face_parent)[iface];
    int next_parent_face_id = -1;
    if (parent_face_id != -1)
    {
        next_parent_face_id = (*face_parent)[parent_face_id];
        if (next_parent_face_id == -1)
        {
            return parent_face_id;
        }
        else
        {
            return GetRootFaceID(next_parent_face_id);
        }
    }
    else
    {
        return iface;
    }
}

void ElementFaceProxy::MarkAllCandidateBoundaryFace()
{
    int nTotalFace = this->GetNTotalFace();

    for (int iface = 0; iface < nTotalFace; ++ iface)
    {
        int re = (*face_to_cell_right)[iface];

        if (re == -1)
        {
            //说明具有候选资格
            if ((*face_child)[iface].size() == 0)
            {
                //说明无子面
                int parent_face_id = (*face_parent)[iface];
                //说明此单元面是简单物理面，应该已经mark
                if (parent_face_id == -1) continue;

                int bctype = (*face_boundary_type)[parent_face_id];
                if (bctype == 0)
                {
                    //说明为内部边界
                    (*face_to_cell_right)[iface] = (*face_to_cell_right)[parent_face_id];
                }
                else
                {
                    (*face_boundary_type)[iface] = (*face_boundary_type)[parent_face_id];
                }
            }
        }
    }
}

void ElementFaceProxy::SetFaceBoundary(vector<int> &face_node, int face_type, int bctype)
{
    sort(face_node.begin(), face_node.end());

    MultiIndex face_of_multi_index(static_cast<int>(face_node.size()));
    face_of_multi_index.SetData(&face_node[0]);

    int id = FindFace(face_of_multi_index);
    if (id == -1)
    {
        ostringstream oss;
        oss << "impossible ! can't FindFace(face_of_multi_index, iter) \n";
        oss << "bctype = " << bctype << "\n";
        oss << " face_node size = " << face_node.size() << "\n";
        for (std::size_t i = 0; i < face_node.size(); ++ i)
        {
            oss << face_node[i] << " ";
        }
        oss << "\n";
        TK_Exit::ExceptionExit(oss.str());
    }
    else
    {
        //face已经存在
        (*face_boundary_type)[id] = bctype;
    }
}

void ElementFaceProxy::SplitElementBoundaryFace()
{
    /*
    //分裂边界的过程不加入整体的单元面，只对边界单元调整
    int nTotalFace = this->GetNTotalFace();
    int size = boundary_face_bctype->size();
    for (int ibcface = 0; ibcface < size; ++ ibcface)
    {
        int bctype    = (*boundary_face_bctype)[ibcface];
        int face_type = (*boundary_face_element_type)[ibcface];

        //检测单元是否是基本单元
        if (!ISBasicElement(face_type)) continue;

        //element boundary face refine
        //tri_3 -> tri_6

        int face_type_composite   = GetRefineElemType   (face_type         );
        int node_number           = GetElementNodeNumbers(face_type         );
        int node_number_composite = GetElementNodeNumbers(face_type_composite);

        //将单元的类型改为复合单元的类型
        (*boundary_face_element_type)[ibcface] = face_type_composite;
        //将单元的节点数增加到符合单元的节点数
        //这里必须是引用！！！否则你想改谁的值呢？老大！
        vector<int> &face_node = (*boundary_face_to_nodes)[ibcface];
        face_node.resize(node_number_composite);

        for (int i = node_number; i < node_number_composite; ++ i)
        {
            //找出所有相关点，如果为棱中点，则相关点为棱端点
            //如果为面，则相关点为四个角点甚至更多
            vector<int> &mplist = GetRelatedPointListForMiddlePointComputation(face_type_composite, i);

            PointOperation::point_type centerpoint =  point_operation->GetCenterPoint(face_node, mplist);
            //如果中心点已经存在，即无论是棱中点还是面心，都只需找出此点在整体点列里面的序号
            //如果中心点不存在，依然要将此新点加入，并返回序号index
            int index = point_operation->GetPointIndex(centerpoint);
            face_node[i] = index;

            point_operation->CheckBCPoint(face_node, mplist, index);
        }
        AddChildBoundaryFace(face_node, face_type_composite, bctype);
    }
    */
}

void ElementFaceProxy::AddChildBoundaryFace(vector<int> &face_node, int face_type, int bctype)
{
    //AddChildBoundaryFace会导致容器失效，所以face_node不可用，需要重新赋值
    vector<int> face_node_fixed = face_node;
    BasicElement *basic_element = GetBasicElement(face_type);

    //子单元数目
    uint_t subcell_number = GetChildElementNumber(face_type);

    vector<int> ptlist; 
    for (int isubcell = 0; isubcell < subcell_number; ++ isubcell)
    {
        //子单元的类型
        int subcell_type        = basic_element->GetChildElementType(isubcell);
        //子单元相对于当前单元的下标
        vector<int> &child_index = basic_element->GetChildIndex(isubcell); 
        //子单元节点数
        int node_number = GetElementNodeNumbers(subcell_type);

        ptlist.resize(0);
        for (int inode = 0; inode < node_number; ++ inode)
        {
            ptlist.push_back(face_node_fixed[ child_index[inode] ]);
        }

        boundary_face_bctype->push_back(bctype);
        boundary_face_element_type->push_back(subcell_type);
        boundary_face_to_nodes->push_back(ptlist);
    }
}

int ElementFaceProxy::GetNTFace() const
{
    int nTotalFace = this->GetNTotalFace();
    int nTFace = 0;
    for (int iface = 0; iface < nTotalFace; ++ iface)
    {
        if ((*face_child)[iface].size() == 0)
        {
            ++ nTFace;
        }
    }
    return nTFace;
}

int ElementFaceProxy::GetNBoundFace()  const
{
    int nTotalFace = this->GetNTotalFace();
    int nBoundFace = 0;
    for (int iface = 0; iface < nTotalFace; ++ iface)
    {
        int f2cr = (*face_to_cell_right)[iface];
        if ((f2cr == - 1) && ((*face_child)[iface].size() == 0))
        {
            ++ nBoundFace;
        }
    }
    return nBoundFace;
}

void ElementFaceProxy::Normalize(vector< int > &element_all_to_local_map)
{
    int nTotalFace = this->GetNTotalFace();

    for (int iface = 0; iface < nTotalFace; ++ iface)
    {
        int le = (*face_to_cell_left)[iface];
        int re = (*face_to_cell_right)[iface];

        (*face_to_cell_left)[iface] = element_all_to_local_map[le];
        if (re == -1)
        {
            continue;
        }
        (*face_to_cell_right)[iface] = element_all_to_local_map[re];
    }
}

void ElementFaceProxy::SwapElementFace(int iface, int jface)
{
    //Swap Face iface, jface
    SWAP((*face_to_cell_left)[iface], (*face_to_cell_left)[jface]);
    SWAP((*face_to_cell_right)[iface], (*face_to_cell_right)[jface]);
    SWAP((*face_boundary_type)[iface], (*face_boundary_type)[jface]);
    SWAP((*face_to_nodes)[iface]     , (*face_to_nodes)[jface]     );
}

void ElementFaceProxy::ReorderElementFace()
{
    int nTFace = this->GetNTFace();
    int nBoundFace = this->GetNBoundFace();

    int isort_start = nBoundFace;
    for (int iface = 0; iface < nBoundFace; ++ iface)
    {
        if ((*face_to_cell_right)[iface] == -1) continue;
        //满足条件的属于内边界，应该移到后面去（nBoundFace之后）

        //需要交换
        for (int jface = isort_start; jface < nTFace; ++ jface)
        {
            //第j个面属于边界面
            if ((*face_to_cell_right)[jface] == -1)
            {
                SwapElementFace(iface, jface);
                isort_start = jface + 1;
                break;
            }
        }
    }
}

void ElementFaceProxy::PrintBoundaryFace()
{
    int nTFace = this->GetNTFace();

    fstream file;
    file.open("boundaryface.plt",ios_base::out);

    for (int iface = 0; iface < nTFace; ++ iface)
    {
        int re = (*face_to_cell_right)[iface];

        if (re == -1)
        {
            file << "iface = " << iface << " nTFace = " << nTFace;
            for (std::size_t i = 0; i < (*face_to_nodes)[iface].size(); ++ i)
            {
                file << " " << (*face_to_nodes)[iface][i];
            }
            file << "\n";
        }
    }
    file.close();
    file.clear();
}

map<int,string> ElementFaceProxy::CreatFantasybctype()
{
    map<int,string> Fantasybctype;

    Fantasybctype[-2]  = "WAKE";
    Fantasybctype[-1]  = "INTERFACE";
    Fantasybctype[0]   = "NO_BOUNDARY_CONDITION";
    Fantasybctype[1]   = "EXTRAPOLATION";
    Fantasybctype[2]   = "SOLID_SURFACE";
    Fantasybctype[3]   = "SYMMETRY";
    Fantasybctype[4]   = "FARFIELD";
                      
    Fantasybctype[5]   = "INFLOW";
    Fantasybctype[52]  = "PRESSURE_INLET";

    Fantasybctype[6]   = "OUTFLOW";
    Fantasybctype[61]  = "OUTFLOW_CONFINED";
    Fantasybctype[62]  = "PRESSURE_OUTLET";

    Fantasybctype[7]   = "POLE";
    Fantasybctype[71]  = "POLE1";
    Fantasybctype[72]  = "POLE2";
    Fantasybctype[73]  = "POLE3";

    Fantasybctype[8]   = "GENERIC_1";
    Fantasybctype[9]   = "GENERIC_2";
    Fantasybctype[10]  = "GENERIC_3";

    Fantasybctype[11]  = "StartNumberOfBodyOfHyperFLOW";

    Fantasybctype[12]  = "BODY";
    Fantasybctype[13]  = "BODY";
    Fantasybctype[14]  = "BODY";
    Fantasybctype[15]  = "BODY";
    Fantasybctype[16]  = "BODY";
    Fantasybctype[17]  = "BODY";
    Fantasybctype[18]  = "BODY";
    Fantasybctype[19]  = "BODY";
    Fantasybctype[20]  = "BODY";
    Fantasybctype[21]  = "BODY";
    Fantasybctype[22]  = "BODY";
    Fantasybctype[23]  = "BODY";
    Fantasybctype[24]  = "BODY";
    Fantasybctype[25]  = "BODY";
    Fantasybctype[26]  = "BODY";
    Fantasybctype[27]  = "BODY";
    Fantasybctype[28]  = "BODY";
    Fantasybctype[29]  = "BODY";
    Fantasybctype[30]  = "ENDNumberOfBodyOfHyperFLOW";

    Fantasybctype[99]  = "PERIODIC";
    Fantasybctype[100] = "INTERIOR";

    Fantasybctype[1001] = "PRESSURE_INLET_PLATE";
    Fantasybctype[1002] = "PRESSURE_OUTLET_PLATE";

    //const int EXTERNAL_BC           =   9;
    ////Added by Guo Yongheng 20160815: external overlap cells label

    return Fantasybctype;
}

void ElementFaceProxy::CreateUnstructGrid(UnstructGrid *grid)
{
    cout << "设置边界条件......\n";
    //设置边界条件
    int nBoundFace     = this->GetNBoundFace();
    int nTotalFace = this->GetNTotalFace();
    CompositeUnstructBC **bcrs = new CompositeUnstructBC * [ nBoundFace ];

    BaseBCType *bc_convert = new BaseBCType();

    set<int> s1,s2;
    HyList<int> bcList;    // Bell 20131204 add
    int icount = 0;
    for (int iface = 0; iface < nTotalFace; ++ iface)
    {
        int f2cr = (*face_to_cell_right)[iface];
        if ((f2cr == - 1) && ((*face_child)[iface].size() == 0))
        {
            bcrs[icount] = new CompositeUnstructBC();
            int bctype = bc_convert->GetBCType((*face_boundary_type)[iface]);
            bcrs[icount]->SetKey(bctype);
            s1.insert((*face_boundary_type)[iface]);
            s2.insert(bctype);
            bcList.insert(bctype);
            ++ icount;
        } 
    }
    grid->SetBCRecord(bcrs);

    delete bc_convert;

    cout << "LOOK AT THIS .............................\n";
    cout << "original bc num are " << s1.size() << "\n";
    for (set<int>::iterator iter = s1.begin(); iter != s1.end(); ++ iter)
    {
        cout << *iter << " ";
    }
    cout << "\n";

    map<int,string> Fantasybctype = CreatFantasybctype();

    cout << "final bc num are " << s2.size() << "\n";
    for (int iBC = 0; iBC < bcList.size(); ++ iBC)
    {
        cout<< " the "<< iBC + 1<< " BCName is ";
        cout << Fantasybctype[bcList[ iBC ]] << "\n";
    }
//     for (set<int>::iterator iter = s2.begin(); iter != s2.end(); ++ iter)
//     {
//         cout << *iter << " ";
//     }
    //cout << "\n";

    int nTFace = GetNTFace();
    grid->SetNTotalFace(nTFace);
    cout << "nTFace = " << nTFace <<"\n";

    int * node_number_of_each_face = new int [nTFace];
    icount = 0;
    //先处理边界
    for (int iface = 0; iface < nTotalFace; ++ iface)
    {
        int f2cr = (*face_to_cell_right)[iface];
        if ((f2cr == - 1) && ((*face_child)[iface].size() == 0))
        {
            node_number_of_each_face[icount] = GetFacePointNumber(iface);
            ++ icount;
        }
    }

    //再处理内部面
    for (int iface = 0; iface < nTotalFace; ++ iface)
    {
        int f2cr = (*face_to_cell_right)[iface];
        if ((f2cr != - 1) && ((*face_child)[iface].size() == 0))
        {
            node_number_of_each_face[icount] = GetFacePointNumber(iface);
            ++ icount;
        }
    }
    grid->SetNodeNumberOfEachFace(node_number_of_each_face);
    cout << " Set node_number_of_each_face finished\n";

    int face2node_size = 0;
    icount = 0;
    //先处理边界
    for (int iface = 0; iface < nTotalFace; ++ iface)
    {
        int f2cr = (*face_to_cell_right)[iface];
        if ((f2cr == - 1) && ((*face_child)[iface].size() == 0))
        {
            face2node_size += node_number_of_each_face[icount];
            ++ icount;
        }
    }

    //再处理内部面
    for (int iface = 0; iface < nTotalFace; ++ iface)
    {
        int f2cr = (*face_to_cell_right)[iface];
        if ((f2cr != - 1) && ((*face_child)[iface].size() == 0))
        {
            face2node_size += node_number_of_each_face[icount];
            ++ icount;
        }
    }

    int * face2node  = new int [ face2node_size ];
    grid->SetFace2Node(face2node);

    int count = 0;
    icount = 0;

    //先处理边界
    for (int iface = 0; iface < nTotalFace; ++ iface)
    {
        int f2cr = (*face_to_cell_right)[iface];
        if ((f2cr == - 1) && ((*face_child)[iface].size() == 0))
        {
            uint_t size = (*face_to_nodes)[iface].size();
            for (int n = 0; n < size; ++ n)
            {
                face2node[count] = (*face_to_nodes)[iface][n];
                ++ count;
            }
            ++ icount;
        }
    }

    //再处理内部面
    for (int iface = 0; iface < nTotalFace; ++ iface)
    {
        int f2cr = (*face_to_cell_right)[iface];
        if ((f2cr != - 1) && ((*face_child)[iface].size() == 0))
        {
            uint_t size = (*face_to_nodes)[iface].size();
            for (int n = 0; n < size; ++ n)
            {
                face2node[count] = (*face_to_nodes)[iface][n];
                ++ count;
            }
            ++ icount;
        }
    }

    cout << " Set face2node finished\n";

    int * left_cell_of_face = new int[ nTFace ];
    int * right_cell_of_face = new int[ nTFace ];
    grid->SetLeftCellOfFace(left_cell_of_face);
    grid->SetRightCellOfFace(right_cell_of_face);

    icount = 0;
    for (int iface = 0; iface < nTotalFace; ++ iface)
    {
        int f2cr = (*face_to_cell_right)[iface];
        if ((f2cr == - 1) && ((*face_child)[iface].size() == 0))
        {
            //先考虑边界面
            left_cell_of_face[ icount ] = (*face_to_cell_left)[iface];
            right_cell_of_face[ icount ] = (*face_to_cell_right)[iface];
            ++ icount;
        }
    }

    icount = nBoundFace;
    for (int iface = 0; iface < nTotalFace; ++ iface)
    {
        int f2cr = (*face_to_cell_right)[iface];
        if ((f2cr != - 1) && ((*face_child)[iface].size() == 0))
        {
            //先考虑边界面
            left_cell_of_face[ icount ] = (*face_to_cell_left)[iface];
            right_cell_of_face[ icount ] = (*face_to_cell_right)[iface];
            ++ icount;
        }
    }

}

SuperElement::SuperElement(CGNSZone *cgns_zone, int iadapt)
{
    this->cgns_zone  = cgns_zone;
    this->iadapt     = iadapt;
    point_operation    = new PointOperation();
    element_face_proxy = new ElementFaceProxy(point_operation);
    elem_type   = new vector<int>;
    elem_level  = new vector<int>;
    elem_modify = new vector<int>;
    elem_parent = new vector<int>;
    elem_mat    = new vector<int>;
}

SuperElement::~SuperElement()
{
    delete point_operation;
    delete element_face_proxy;
    delete elem_type;
    delete elem_level;
    delete elem_modify;
    delete elem_parent;
    delete elem_mat;
}

void SuperElement::CreateFantasyGrid(Grid *grid, int grid_type)
{
    if (grid_type == Unstructured)
    {
        CreateFantasyUnstructGrid(grid);
    }
    else
    {
        CreateFantasyStructGrid(grid);
    }
}

void SuperElement::CreateFantasyStructGrid(Grid *grid_in)
{
    cgns_zone->CopyStructuredGrid(grid_in);
}

void SuperElement::CreateFantasyUnstructGrid(Grid *grid_in)
{
    UnstructGrid *grid = UnstructGridCast (grid_in);
    int nTotalCell = element_face_proxy->GetNTotalCell();
    int nTotalNode = point_operation->GetPointNumber();
    int nBoundFace = element_face_proxy->GetNBoundFace();

    grid->SetNTotalCell(nTotalCell);
    grid->SetNTotalNode(nTotalNode);
    grid->SetNBoundFace(nBoundFace);

    RDouble *x = new RDouble[nTotalNode];
    RDouble *y = new RDouble[nTotalNode];
    RDouble *z = new RDouble[nTotalNode];

    grid->SetX(x);
    grid->SetY(y);
    grid->SetZ(z);

    for (int inode = 0; inode < nTotalNode; ++ inode)
    {
        x[inode] = point_operation->GetPoint(inode).X();
        y[inode] = point_operation->GetPoint(inode).Y();
        z[inode] = point_operation->GetPoint(inode).Z();
    }

    element_face_proxy->CreateUnstructGrid(grid);

    cout << "所有的计算信息准备完毕\n";
}

void SuperElement::InitElement()
{
    cgns_zone->Convert2InerDataStandard();
    cgns_zone->SetElementTypeAndIndex(*elem_type,elem_index);
    cgns_zone->BuildPoints(point_operation);
    SetElementSplitInfo();
}

void SuperElement::SetElementSplitInfo()
{
    //当前单元数量,在加密或稀疏化的时候，网格单元数量可能会变化
    uint_t element_number = this->GetElementNumber();
    this->elem_modify->resize(element_number);

    for (int ie = 0; ie < element_number; ++ ie)
    {
        (*elem_modify)[ie] = PHSPACE::NOCHANGE;
    }
    (*elem_modify)[0] = PHSPACE::SPLIT;
    //(*elem_modify)[3] = PHSPACE::SPLIT;
}

bool SuperElement::ElementNeedSplit(int ielement)
{
    if ((*elem_modify)[ielement] == PHSPACE::SPLIT)
    {
        return true;
    }
    return false;
}

void SuperElement::ResolveElementFace(vector<int> &element_index, int element_type, int ielement)
{
    element_face_proxy->ResolveElementFace(element_index, element_type, ielement);
}

void SuperElement::ResolveElementFace()
{
    //当前单元数量,在加密或稀疏化的时候，网格单元数量可能会变化
    uint_t element_number = this->GetElementNumber();

    int count = 0;
    for (int ie = 0; ie < element_number; ++ ie)
    {
        int elem_type = GetElementType(ie);

        ////检测单元是否是基本单元
        if (!ISBasicElement(elem_type)) continue;

        ++ count;
        if (count % 200000 == 0) cout << " ie = " << ie << " total element number = " << element_number << "\n";

        ResolveElementFace(elem_index[ie], elem_type, ie);
    }

    cout << "ResolveElementFace face number = " << element_face_proxy->GetNTFace() << "\n";
}

void SuperElement::ResolveElementBoundaryFace()
{
    cgns_zone->ResolveElementBoundaryFace(element_face_proxy);
    element_face_proxy->MarkElementBoundaryFace();
}

void SuperElement::ResolveElementBoundaryFaceFinal()
{
    element_face_proxy->ResolveElementBoundaryFace(iadapt);
}

int SuperElement::GetElementType(int ie) const
{
    vector<int> &elem_type = *this->elem_type;
    return elem_type[ie];
}

void SuperElement::SetElementType(uint_t ie, int new_elem_type)
{
    vector<int> &elem_type = *this->elem_type;
    elem_type[ie] = new_elem_type;
}

void SuperElement::SetElementIndex(uint_t ie, vector<int> &ptlist)
{
    elem_index[ie] = ptlist;
}

void SuperElement::SetElementParent(uint_t ie, int parent)
{
    vector<int> &elem_parent = *this->elem_parent;
    elem_parent[ie] = parent;
}

void SuperElement::SetElementLevel(uint_t ie, int level)
{
    vector<int> &elem_level = *this->elem_level;
    elem_level[ie] = level;
}

void SuperElement::SetElementSplitProperty(uint_t ie, int split_property)
{
    (*elem_modify)[ie] = split_property;
}

int SuperElement::GetLevel(int ie) const
{
    vector<int> &elem_level = *this->elem_level;
    return elem_level[ie];
}

uint_t SuperElement::GetElementNumber() const
{
    return elem_type->size();
}

void SuperElement::ResizeElementNumber(uint_t n)
{
    //所有与单元相关的都要同步变动容量
    this->elem_index.resize(n);
    this->elem_type->resize(n);
    this->elem_level->resize(n);
    this->elem_modify->resize(n);
    this->elem_mat->resize(n);
    this->elem_parent->resize(n);
}

void SuperElement::AddElement(vector<int> &ptlist, int type, int level, int parent)
{
    uint_t index = GetElementNumber();
    ResizeElementNumber(index + 1);

    SetElementIndex(index,ptlist);
    SetElementType(index,type);
    SetElementParent(index,parent);
    SetElementLevel(index,level); 
    SetElementSplitProperty(index, PHSPACE::SPLIT);
}

void SuperElement::HideElement(int ie)
{
/*  1. Flag the cell with <HIDDEN>. */
/*  2. Delete from <ptoe> of all points of it. */
/*  3. Delete all middle point structure of THE cell. */
/*  4. Add the array of <son[]>. */

    vector<int> &elem_modify = *this->elem_modify;
    elem_modify[ie] = PHSPACE::HIDDEN;
}

void SuperElement::AddChildElement(int ie)
{
    int type  = GetElementType(ie);
    //int level = GetLevel(ie);
    int level = 0;

    BasicElement *basic_element = GetBasicElement(type);

    //子单元数目
    uint_t subcell_number = GetChildElementNumber(type);

    vector<int> ptlist; 
    for (int i = 0; i < subcell_number; ++ i)
    {
        //子单元的类型
        int subcell_type        = basic_element->GetChildElementType(i);
        //子单元相对于当前单元的下标
        vector<int> child_index = basic_element->GetChildIndex(i); 
        //子单元节点数
        int node_number = GetElementNodeNumbers(subcell_type);

        //当前单元的节点标号
        vector<int> & parent_index = GetElementIndex(ie);

        ptlist.resize(0);
        for (int j = 0; j < node_number; ++ j)
        {
            ptlist.push_back(parent_index[ child_index[j] ]);
        }

        //这里有一个很严重的问题，resize后可能会导致引用失效，切记切记!

        AddElement(ptlist, subcell_type, level + 1, ie);
    }
}

void SuperElement::SplitElement()
{
    //当前单元数量,在加密或稀疏化的时候，网格单元数量可能会变化
    uint_t element_number = this->GetElementNumber();

    for (int ie = 0; ie < element_number; ++ ie)
    {
        SplitElement(ie);
    }
}

void SuperElement::CreateElementGlobal2LocalMap()
{
    //当前单元数量,在加密或稀疏化的时候，网格单元数量可能会变化
    uint_t element_number = this->GetElementNumber();

    element_all_to_local_map.resize(element_number);

    int count = 0;
    for (int ie = 0; ie < element_number; ++ ie)
    {
        int elem_type = GetElementType(ie);

        //检测单元是否是基本单元
        if (!ISBasicElement(elem_type))
        {
            element_all_to_local_map[ie] = -1;
        }
        else
        {
            element_all_to_local_map[ie] = count;
            ++ count;
        }
    }

    element_face_proxy->SetNTotalCell(count);
}

void SuperElement::SplitElement(int ie)
{
    /*
    //实际上既然加密cell的策略就是所有棱中点处和面心中点处都要加点
    //那么不如索性一次加入所有点，如果有变化，可以形成一定的策略
    int elem_type = GetElementType(ie);

    //如果本单元不需要加密，则直接返回
    vector<int> &elem_modify = *this->elem_modify;
    //if (elem_modify[ie] != AdaptEnum::SPLIT) return;

    //检测单元是否是基本单元
    if (!ISBasicElement(elem_type)) return;
    if (!ElementNeedSplit(ie)) return;

    //element refine
    //tri_3 -> tri_6

    int elem_type_composite   = GetRefineElemType   (elem_type         );
    int node_number           = GetElementNodeNumbers(elem_type         );
    int node_number_composite = GetElementNodeNumbers(elem_type_composite);

    //将单元的类型改为复合单元的类型
    SetElementType(ie, elem_type_composite);
    //将单元的节点数增加到符合单元的节点数
    elem_index[ie].resize(node_number_composite);

    for (int i = node_number; i < node_number_composite; ++ i)
    {
        //找出所有相关点，如果为棱中点，则相关点为棱端点
        //如果为面，则相关点为四个角点甚至更多
        vector<int> &mplist = GetRelatedPointListForMiddlePointComputation(elem_type_composite, i);

        PointOperation::point_type centerpoint = point_operation->GetCenterPoint(elem_index[ie], mplist);
        //如果中心点已经存在，即无论是棱中点还是面心，都只需找出此点在整体点列里面的序号
        //如果中心点不存在，依然要将此新点加入，并返回序号index
        int index = point_operation->GetPointIndex(centerpoint);
        elem_index[ie][i] = index;
        point_operation->CheckBCPoint(elem_index[ie], mplist, index);
    }

    AddChildElement(ie);

    SplitElementFace(ie); 
    */
}

int SuperElement::AddElementFace(vector<int> &facenodes, int parent_face_id)
{
    return element_face_proxy->AddElementFace(facenodes, parent_face_id);
}

void SuperElement::SplitCompositeFace(vector<int> &facenodes, int composite_face_type)
{
    BasicElement *basic_element = GetBasicElement(composite_face_type);
    //子单元数目，在这里实际上代表着子面的数目，因为基本单元实际上是面
    uint_t child_element_number = GetChildElementNumber(composite_face_type);

    vector<int> simple_face_nodes;
    int simple_face_type        = GetSimpleElementType(composite_face_type);
    int simple_face_node_number = GetElementNodeNumbers(simple_face_type);
    for (int inode = 0; inode < simple_face_node_number; ++ inode)
    {
        simple_face_nodes.push_back(facenodes[inode]);
    }

    int parent_face_id = element_face_proxy->FindFace(simple_face_nodes);

    //vector<int> element_face;
    //vector<int> sort_face;

    vector<int> ptlist; 
 
    //这里的iface实际上是单元，只不过表现形式为面的形式
    for (int ichildelem = 0; ichildelem < child_element_number; ++ ichildelem)
    {
        //子单元的类型
        int child_element_type = basic_element->GetChildElementType(ichildelem);

        //子单元相对于当前单元的下标
        vector<int> child_index = basic_element->GetChildIndex(ichildelem); 

        //子单元节点数
        int child_node_number = GetElementNodeNumbers(child_element_type);

        //当前单元的节点编号，这里为facenodes，相当于父单元的节点编号

        ptlist.resize(0);
        for (int j = 0; j < child_node_number; ++ j)
        {
            ptlist.push_back(facenodes[ child_index[j] ]);
        }

        //每一个物理面包含child_face_number个子面，应该
        int child_face_id = AddElementFace(ptlist, parent_face_id);
    }
}

void SuperElement::SplitElementFace(int ielement)
{
    int type  = GetElementType(ielement);
    //int level = GetLevel(ielement);

    //获取当前单元类型
    BasicElement *basic_element = GetBasicElement(type);

    //获取本单元物理面的个数
    uint_t element_physical_face_number = basic_element->GetElementPhysicsFaceNumber();

    //对每个物理面进行子面的分裂
    for (int icompositeface = 0; icompositeface < element_physical_face_number; ++ icompositeface)
    {
        int composite_face_type    = basic_element->GetElementPhysicsFaceType(icompositeface);
        vector<int> facenodeslocal = basic_element->GetElementPhysicsFace(icompositeface);
        vector<int> facenodes;
        vector<int> & parent_index = GetElementIndex(ielement);
        for (std::size_t inode = 0; inode < facenodeslocal.size(); ++ inode)
        {
            facenodes.push_back(parent_index[ facenodeslocal[inode] ]);
        }

        //对单个物理面进行子面的分裂
        SplitCompositeFace(facenodes, composite_face_type);
    }
}

void SuperElement::FreeMemory()
{
    element_face_proxy->FreeMemory();
}

void SuperElement::PrintBoundaryFace()
{
    element_face_proxy->PrintBoundaryFace();
}

void SuperElement::MarkAllCandidateBoundaryFace()
{
    element_face_proxy->MarkAllCandidateBoundaryFace();
}

int SuperElement::GetGridType() const
{
    return cgns_zone->GetZoneType();
}

void SuperElement::ResolveComputationStructuredGrid()
{
}

void SuperElement::ResolveComputationUnstructuredGrid()
{
    InitElement();
    //解析所有单元面
    //我们可以发现单元面的解析实际上并不依赖于CGNS数据结构
    ResolveElementFace();
    //解析单元边界面
    //我们可以发现边界面的解析实际上依赖于CGNS数据结构
    ResolveElementBoundaryFace();

    GridAdaption();

    ConstructBasicData();
}

void SuperElement::GridAdaption()
{
    for (int i = 0; i < iadapt; ++ i)
    {
        SplitElement();
    }
}

void SuperElement::ConstructBasicData()
{
    //FreeMemory();
    //再次解析所有单元面
    ResolveElementFace();
    MarkAllCandidateBoundaryFace();
    //PrintBoundaryFace
     
    //解析单元边界面
    //ResolveElementBoundaryFaceFinal();

    CreateElementGlobal2LocalMap();

    Normalize();
}

void SuperElement::ResolveComputationGrid()
{
    if (GetGridType() == Unstructured)
    {
        ResolveComputationUnstructuredGrid();
    }
    else
    {
        ResolveComputationStructuredGrid();
    }
}

void SuperElement::Normalize()
{
    element_face_proxy->Normalize(element_all_to_local_map);
}

CGNSGridFactory::CGNSGridFactory()
{
    nzones         = 0;
    super_elements = 0;
    cgns_multizone = NULL;
}

CGNSGridFactory::~CGNSGridFactory()
{
    for (int izone = 0; izone < nzones; ++ izone)
    {
        delete super_elements[izone];
    }
    delete [] super_elements;
}

void CGNSGridFactory::ConstructGrid()
{
    ReadGrid();

    //从读入的标准网格信息（CGNS）解析出实际计算需要的网格信息
    ResolveComputationGrid();

    CreateGeneralMultiBlockGrid();
}

void CGNSGridFactory::ReadGrid()
{
    cgns_multizone = new CGNSMultiZone();

    string from_gfile = "grid.fts";
    GlobalDataBase::GetData("from_gfile", &from_gfile, PHSTRING, 1);

    cgns_multizone->ReadGrid(from_gfile);
}

void CGNSGridFactory::CreateGeneralMultiBlockGrid()
{
    Grid ** grids  = new Grid * [nzones];

    CreateAllSingleGrid(grids);

    ConstructGlobalInterfaceLink(grids, nzones);

    string out_gfile = "grid.fts";
    GlobalDataBase::GetData("out_gfile", &out_gfile, PHSTRING, 1);

    DumpGrid(out_gfile, grids, nzones);

    for (int izone = 0; izone < nzones; ++ izone)
    {
        delete grids[izone];
    }
    delete [] grids;
}

int CGNSGridFactory::GetGridType(int izone) const
{
    return cgns_multizone->GetZone(izone)->GetZoneType();
}

Grid * CGNSGridFactory::CreateBasicGrid(int izone)
{
    int grid_type = GetGridType(izone);
    if (grid_type == Unstructured)
    {
        return CreateGridGeneral(UNSTRUCTGRID, new GridID(izone), 0, GetDim());
    }
    else
    {
        return CreateGridGeneral(STRUCTGRID, new GridID(izone), 0, GetDim());
    }
}

void CGNSGridFactory::CreateAllSingleGrid(Grid ** grids)
{
    for (int izone = 0; izone < nzones; ++ izone)
    {
        Grid *grid  = CreateBasicGrid(izone);
        grids[izone] = grid;
        super_elements[izone]->CreateFantasyGrid(grid, GetGridType(izone));
    }
}

void CGNSGridFactory::ResolveComputationGrid()
{
    nzones = cgns_multizone->GetNZones();

    int iadapt = 0;
    GlobalDataBase::GetData("iadapt", &iadapt, PHINT, 1);

    super_elements = new SuperElement * [nzones];
    for (int izone = 0; izone < nzones; ++ izone)
    {
        super_elements[izone] = new SuperElement(cgns_multizone->GetZone(izone), iadapt);
    }

    for (int izone = 0; izone < nzones; ++ izone)
    {
        super_elements[izone]->ResolveComputationGrid();
    }
}



void PointOperation::CheckBCPoint(vector<int> &elem, vector<int> &point_list, int index)
{
    uint_t node_number = point_list.size();

    int bcp = 0;

    for (int i = 0; i < node_number; ++ i)
    {
        if (PointWithBoundary(elem[ point_list[i] ], BCWall))
        {
            ++ bcp;
        }
    }

    if (bcp >= 2)
    {
        InsertPointBCType(index, BCWall);
    }
}

PointOperation::point_type PointOperation::CirclePoint(vector<int> &elem, vector<int> &point_list)
{
    RDouble xc = zero;
    RDouble yc = zero;
    RDouble zc = zero;

    RDouble radius = one;

    uint_t node_number = point_list.size();

    for (int i = 0; i < node_number; ++ i)
    {
        PointOperation::point_type pt = GetPoint(elem[ point_list[i] ]);
        xc += pt.X();
        yc += pt.Y();
        zc += pt.Z();
    }

    RDouble o_node_number = one / (node_number + SMALL);

    xc *= o_node_number;
    yc *= o_node_number;
    zc *= o_node_number;

    RDouble osc = one / (sqrt(xc * xc + yc * yc + zc * zc) + SMALL);

    RDouble x0 = xc * osc * radius;
    RDouble y0 = yc * osc * radius;
    RDouble z0 = zc * osc * radius;

    PointOperation::point_type center_point(x0,y0,z0);
    return center_point;
}

void PointOperation::SetBoundaryPoint(set<cgsize_t> &bc_vertex, int bctype)
{
    for(set<cgsize_t>::iterator iter = bc_vertex.begin(); iter != bc_vertex.end(); ++ iter)
    {
        point_bc_type[*iter].insert(bctype);
    }
}

bool PointOperation::PointWithBoundary(int index, int bctype)
{
    set<int>::iterator iter = point_bc_type[index].find(bctype);
    if (iter != point_bc_type[index].end())
    {
        return true;
    }
    return false;
}

void PointOperation::InsertPointBCType(int index, int bctype)
{
    point_bc_type[index].insert(bctype);
}


PointOperation::point_type PointOperation::GetCenterPoint(vector<int> &elem, vector<int> &point_list)
{
    RDouble xc = zero;
    RDouble yc = zero;
    RDouble zc = zero;

    uint_t node_number = point_list.size();

    int bcp = 0;

    for (int i = 0; i < node_number; ++ i)
    {
        if (PointWithBoundary(elem[ point_list[i] ], BCWall))
        {
            ++ bcp;
        }
    }

    if (bcp >= 2)
    {
        //return CirclePoint(elem, point_list);
    }

    for (int i = 0; i < node_number; ++ i)
    {
        PointOperation::point_type pt = GetPoint(elem[ point_list[i] ]);
        xc += pt.X();
        yc += pt.Y();
        zc += pt.Z();
    }

    RDouble o_node_number = one / (node_number + SMALL);

    xc *= o_node_number;
    yc *= o_node_number;
    zc *= o_node_number;

    PointOperation::point_type center_point(xc,yc,zc);
    return center_point;
}



}

