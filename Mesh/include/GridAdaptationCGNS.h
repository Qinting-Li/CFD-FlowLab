//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//          PPPPP  H   H  EEEEE  N    N  GGGGG  L      EEEEE  III         +
//          P   P  H   H  E      NN   N  G      L      E       I          +
//          PPPPP  HHHHH  EEEEE  N N  N  G  GG  L      EEEEE   I          +
//          P      H   H  E      N  N N  G   G  L      E       I          +
//          P      H   H  EEEEE  N    N  GGGGG  LLLLL  EEEEE  III         +
//------------------------------------------------------------------------+
//          Platform for Hybrid Engineering Simulation of Flows           +
//          China Aerodynamics Research and Development Center            +
//                     (C) Copyright, Since 2010                          +
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      GridAdaptationCGNS.h
//! @brief     Explain this file briefly.
//! @author    He Xin.

#pragma once
#include "Pre_GridBase.h"

namespace PHSPACE
{
typedef enum
{
    SPLIT,
    NOCHANGE,
    MERGE_ELEMENT,
    DELETE,
    HIDDEN,
    ON, //Green
    BOUNDARYFACE, //
    HANGINGFACE, //
    GENERALFACE
} AdaptEnum;


typedef enum
{
    DENSITYGRADIENT,
    U_VELOCITYGRADIENT,
    V_VELOCITYGRADIENT,
    W_VELOCITYGRADIENT,
    MACHNUMBERGRADIENT,
    PRESSUREGRADIENT,
    VORTICITYMAGNITUDE,
    VORTICITYQCRITERIA,
    ADAPTIVEDISSIPATIVECOEFFICIENT,
} AdaptApproachEnum;

class Grid;
class SimpleGrid;
class StructBC;
class PointOperation;
class ElementFaceProxy;
class CGNSZone;

class CGNSTraits
{
public:
    typedef char char_33[33];
};

class CGNSection
{
public:
    CGNSection(CGNSZone *zone);
    ~CGNSection();
private:
    int element_type;
    cgsize_t istart, iend;
    cgsize_t * connection_list;
    CGNSTraits::char_33 section_name;
    CGNSZone *zone;
    int sectionid;
    int ncoords;
public:
    void SetSectionIndex(int sectionid) { this->sectionid = sectionid; };
    int  GetSectionIndex() const { return sectionid; }

    void SetStart(cgsize_t istart) { this->istart = istart; };
    void SetEnd  (cgsize_t iend ) { this->iend   = iend;   };
    cgsize_t GetStart() const { return istart; };
    cgsize_t GetEnd  () const { return iend;   };

    void SetElementType(int element_type) { this->element_type = element_type; };
    int  GetElementType() const { return element_type; };

    void CreateConnectionList(int element_node_number);
    cgsize_t * GetConnectionList() const { return connection_list; }

    void ReadSection();
public:
    cgsize_t *GetElementIndex(cgsize_t ielem) const;
    void GetElementIndex(int ielement, vector<cgsize_t> &element_index);
    void SetElementTypeAndIndex(vector<int> &element_type, vector< vector<int> > &element_index);
    void Convert2InerDataStandard();
};

class CGNSectionProxy
{
public:
    CGNSectionProxy(CGNSZone *zone);
    ~CGNSectionProxy();
private:
    int nsections;
    CGNSection **sections;
    CGNSZone *zone;
public:
    void Create(int nsections);
    int GetNSection() const { return nsections; };
    CGNSection *GetSection(int isection) const {  return sections[isection]; }
    CGNSection *GetSectionByElementIndex(int ielement) const;
public:
    void Convert2InerDataStandard();
};

class CGNSBCRegion
{
public:
    CGNSBCRegion(CGNSZone *zone);
    ~CGNSBCRegion();
private:
    int bc_type;
    int bc_element_number;
    int bc_point_set_type;
    int bc_grid_location;
    cgsize_t *bc_connection_list;
    int nTotalCell;

    CGNSTraits::char_33 bcregion_name;
    CGNSZone *zone;
public:
    void SetBCType(int bc_type) { this->bc_type = bc_type; };
    int  GetBCType() const { return bc_type; };

    void SetBCPointSetType(int bc_point_set_type) { this->bc_point_set_type = bc_point_set_type; };
    int  GetBCPointSetType() const { return bc_point_set_type; };

    void SetBCElementNumber(int bc_element_number) { this->bc_element_number = bc_element_number; };
    int  GetBCElementNumber() const { return bc_element_number; };

    void SetBCGridLocation(int bc_grid_location) { this->bc_grid_location = bc_grid_location; };
    int  GetBCGridLocation() const { return bc_grid_location; };

    void CreateConnectionList();
    cgsize_t * GetConnectionList() const { return bc_connection_list; }

    void SetNTotalCell(int nTotalCell) { this->nTotalCell = nTotalCell; }
    int  GetNTotalCell() const { return nTotalCell; }
public:
    void Convert2InerDataStandard();
    int  DetermineBase();
    void ShiftBCRegion();
    void ResolveElementBoundaryFace(ElementFaceProxy *element_face_proxy);
    void ProcessBCVertex(set<cgsize_t> &bc_vertex);
    void ProcessBCFace  (set<cgsize_t> &bc_vertex);
    void PrintConnectionList();
    void SetStructuredGridBoundary(StructBC * bcregion);
public:
    void ReadOrdinaryBoundaryRegion(int ibcregion);
    void Read1to1BoundaryRegion(int i1to1);
};

class CGNSBCRegionProxy
{
public:
    CGNSBCRegionProxy(CGNSZone *zone);
    ~CGNSBCRegionProxy();
private:
    int nbcregions, n1to1, n_ordinary_bcregions;
    CGNSBCRegion **bcregions;
    CGNSZone *zone;
public:
    void Create(int n_ordinary_bcregions, int n1to1);
    int GetNBCRegion() const { return nbcregions; };

    int GetN1to1BCRegion() const { return n1to1; }
    void SetN1to1BCRegion(int n1to1) { this->n1to1 = n1to1; }

    int GetNOrdinaryBCRegion() const { return n_ordinary_bcregions; }
    void SetNOrdinaryBCRegion(int n_ordinary_bcregions) { this->n_ordinary_bcregions = n_ordinary_bcregions; }

    CGNSBCRegion *GetBCRegion(int ibcregion) const {  return bcregions[ibcregion]; }
    CGNSBCRegion *Get1to1BCRegion(int i1to1) const { return bcregions[n_ordinary_bcregions + i1to1]; }
public:
    void Convert2InerDataStandard();

};

class CGNSZone
{
public:
    CGNSZone();
    ~CGNSZone();
private:
    SimpleGrid        *grid;
    CGNSectionProxy   *section_proxy;
    CGNSBCRegionProxy *bcregion_proxy;
    int nTotalNode,nTotalCell;
    ZoneType_t zone_type;
    int fileid,baseid,zoneid;
    int cell_dimension, phys_dimension;
    cgsize_t irmin[3],irmax[3];
    int cell_size[3];
    cgsize_t isize[3][3];
    CGNSTraits::char_33 zonename;
    CGNSTraits::char_33 coordname;
    int n1to1, n_ordinary_bcregions;
public:
    cgsize_t  * GetISize() { return isize[0]; };
    cgsize_t  * GetIMin() { return irmin; };
    cgsize_t  * GetIMax() { return irmax; };
    void SetDimension();

    void SetFileIndex(int fileid) { this->fileid = fileid; };
    void SetBaseIndex(int baseid) { this->baseid = baseid; };
    void SetZoneID(int zoneid) { this->zoneid = zoneid; };

    int GetZoneID()  const { return zoneid; }
    int GetFileIndex()  const { return fileid; }
    int GetBaseIndex()  const { return baseid; }

    void CopyStructuredGrid(Grid *grid_in);

    void ReadGrid();

    void ReadGridCoordinates();
    void ReadGridConnectivities();
    void ReadGridBoundaryConditions();

    void ReadCoor();

    void ShowStructuredGrid(SimpleGrid *simple_grid, int ni, int nj, int nk);

    void SetNTotalNode(int nTotalNode) { this->nTotalNode = nTotalNode; }
    void SetNTotalCell(int nTotalCell) { this->nTotalCell = nTotalCell; }

    void SetZoneType(ZoneType_t zone_type) { this->zone_type = zone_type; }
    ZoneType_t GetZoneType() const { return zone_type; }

    int  GetNTotalNode() const { return nTotalNode; }
    int  GetNTotalCell() const { return nTotalCell; }

    void SetCellDimension(int cell_dimension) { this->cell_dimension = cell_dimension; }
    int  GetCellDimension() const { return cell_dimension; }

    void SetPhysDimension(int phys_dimension) { this->phys_dimension = phys_dimension; }
    int  GetPhysDimension() const { return phys_dimension; }
public:
    void Create();
    void CreateSections(int nsections);
    void CreateBCRegions(int n_ordinary_bcregions, int n1to1);
    int GetNSection () const { return section_proxy->GetNSection();   };
    int GetNBCRegion() const { return bcregion_proxy->GetNBCRegion(); };
    int GetN1to1BCRegion() const { return n1to1; }
    int GetNOrdinaryBCRegion() const { return n_ordinary_bcregions; }
    CGNSection   *GetSection (int isection) const {  return section_proxy->GetSection(isection); }
    CGNSBCRegion *GetBCRegion(int ibcregion) const {  return bcregion_proxy->GetBCRegion(ibcregion); }

    CGNSBCRegion *Get1to1BCRegion(int i1to1) const {  return bcregion_proxy->Get1to1BCRegion(i1to1); }
public:
    SimpleGrid *GetSimpleGrid() const { return grid; };
    CGNSectionProxy *GetSectionProxy() const { return section_proxy; }
    void SetElementTypeAndIndex(vector<int> &element_type, vector< vector<int> > &element_index);
    void Convert2InerDataStandard();
    void BuildPoints(PointOperation *point_operation);

    void ResolveElementBoundaryFace(ElementFaceProxy *element_face_proxy);
    void GetElementIndex(int ielement, vector<cgsize_t> &element_index);
};

class UnstructGrid;

class CGNSMultiZone
{
private:
    int fileid, baseid;
    int istart, iend, nbndry, iparent_flag;
    int cellDim, physDim;
    int nTotalNode, nTotalCell, nCoords, nSections;    int normalIndex, normListFlag, nDataSet;
    CGNSTraits::char_33 zonename, sectionname;
    CGNSTraits::char_33 basename, coordname;
private:
    int nzones;
    CGNSZone **zones;
public:
    CGNSMultiZone();
    ~CGNSMultiZone();
public:
    int GetNZones() const { return nzones; }
    CGNSZone *GetZone(int iZone) const {  return zones[iZone]; };
public:
    void InitAllZones();
    void ReadGrid(const string &filename);
public:
    void Convert2InerDataStandard();
    void Convert2CGNSDataStandard(){};
};

class BasicFacet
{
private:
    int *elems;
public:
    BasicFacet(int n);
    ~BasicFacet();
public:
    // operator[]
    int& operator[](int i)
    {
        return elems[i];
    }

    const int& operator[](int i) const
    {
        return elems[i];
    }
}; 

class PointOperation
{
public:
    typedef RDouble value_type;
    typedef Geo_Point< value_type > point_type;
    typedef set< point_type >::iterator iterator;
private:
    set   < point_type > point_list;
    vector< point_type > point_array;
    vector< set<int> > point_bc_type;
public:
    point_type & GetPoint(int id) { return point_array[id]; }
    void SetPoint(int id, point_type &point) { point_array[id] = point; }
    bool FindPoint(point_type &point, iterator &iter)
    {
        iter = point_list.find(point);
        if (iter == point_list.end())
        {
            return false;
        }
        return true;
    }

    void AddPoint(point_type &point, int bctype = 0)
    {
        point_list.insert(point);
        point_array.push_back(point);

        set<int> tmp;
        tmp.insert(bctype);
        point_bc_type.push_back(tmp);
    }

    int GetPointNumber() { return static_cast<int>(point_array.size()); }

    int GetPointIndex(PointOperation::point_type &point)
    {
        PointOperation::iterator iter;
        if (FindPoint(point, iter))
        {
            return iter->ID();
        }
        else
        {
            //此点不存在，将此新点加入点集，并返回index
            int index = GetPointNumber();
            point.SetID(index);
            AddPoint(point);
            return index;
        }
    }
    PointOperation::point_type GetCenterPoint(vector<int> &elem, vector<int> &point_list);
    PointOperation::point_type CirclePoint   (vector<int> &elem, vector<int> &point_list);
    void CheckBCPoint(vector<int> &elem, vector<int> &point_list, int index);
    void SetBoundaryPoint(set<cgsize_t> &bc_vertex, int bctype);
    bool PointWithBoundary(int index, int bctype);
    void InsertPointBCType(int index, int bctype);
public:
    PointOperation();
    ~PointOperation();
public:
};

class ElementFaceProxy
{
public:
    ElementFaceProxy(PointOperation *point_operation);
    ~ElementFaceProxy();
private:
    set< MultiIndex > *reference_faces;
    vector< vector<int> > *face_to_nodes;
    vector< int > *face_to_cell_left;
    vector< int > *face_to_cell_right;
    vector< int > *face_boundary_type;
    vector< int > *face_element_type;
    vector< int > *face_parent;
    vector< vector<int> > *face_child;

    vector< int > *boundary_face_bctype;
    vector< int > *boundary_face_element_type;
    vector< vector<int> > *boundary_face_to_nodes;

    PointOperation *point_operation;
    int nTotalCell;
    map<int,string> CreatFantasybctype();
public:
    int  FindFace(MultiIndex &face);
    int  FindFace(vector<int> &facenodes);
    void InsertFace(MultiIndex &face);
    void ResizeFaceNumber(int n);
    void AddFaceHierachical(int parent_face_id, int child_face_id);
    void SetDefaultFaceInfo(int index, vector<int> &facenodes);
public:
    int  GetNTFace() const;
    int  GetNTotalFace() const { return static_cast<int>(face_to_nodes->size()); }
    void ResolveElementFace(vector<int> &element_index, int element_type, int ielement);
    void ResolveElementBoundaryFace(set<cgsize_t> &bc_vertex, int bctype);
    void ResolveElementBoundaryFace(int iadapt);
    int  AddElementFace(vector<int> &facenodes, int parent_face_id);
    void SplitElementBoundaryFace();
    void MarkElementBoundaryFace();
    void MarkAllCandidateBoundaryFace();
    bool CheckBCFace(set<cgsize_t> &bc_vertex, vector<int> &face_node);
    int  GetRootFaceID(int iface);
    int  GetNBoundFace() const;
    int  GetFacePointNumber(int iface) { return static_cast<int>((*face_to_nodes)[iface].size()); }
    int  GetNTotalCell() const { return nTotalCell; }
    void SetNTotalCell(int nTotalCell) { this->nTotalCell = nTotalCell; }
public: 
    void Normalize(vector< int > &element_all_to_local_map);
    void CreateUnstructGrid(UnstructGrid *grid);
    void ReorderElementFace();
    void SwapElementFace(int iface, int jface);
    void FreeMemory();
    void PrintBoundaryFace();
    void AddChildBoundaryFace(vector<int> &face_node, int face_type, int bctype);
    void SetFaceBoundary     (vector<int> &face_node, int face_type, int bctype);
};

class SuperElement
{
public:
    typedef vector<int> value_type;
private: 
    value_type *elem_type;   //单元类型
    value_type *elem_level;  //单元级别
    value_type *elem_modify; //单元改动
    value_type *elem_parent; //单元parent
    value_type *elem_mat;
    vector< vector<int> > elem_index; //单元index
    vector< int > element_all_to_local_map;
private:
    PointOperation *point_operation;
    ElementFaceProxy *element_face_proxy;
    CGNSZone *cgns_zone;
    int iadapt;
public:
    SuperElement(CGNSZone *cgns_zone, int iadapt = 0);
    ~SuperElement();
public:
    int  GetElementType(int ie) const;
    void SetElementType(uint_t ie, int new_elem_type);
    int  GetLevel(int ie) const;

    vector<int> & GetElementIndex(int ie) { return elem_index[ie]; };

    void SetElementIndex(uint_t ie, vector<int> &ptlist);
    void SetElementParent(uint_t ie, int parent);
    void SetElementLevel(uint_t ie, int level);
    void SetElementSplitProperty(uint_t ie, int split_property);

    uint_t  GetElementNumber() const;
    void ResizeElementNumber(uint_t n);
public:
    void CreateFantasyGrid(Grid *grid, int grid_type);
    void CreateFantasyUnstructGrid(Grid *grid_in);
    void CreateFantasyStructGrid  (Grid *grid_in);
public:
    void InitElement();
    void SetElementSplitInfo();
    void SplitCompositeFace(vector<int> &facenodes, int composite_face_type);
    void ResolveElementFace();
    void ResolveElementFace(vector<int> &element_index, int element_type, int ielement);
    void ResolveElementBoundaryFace();
    void ResolveElementBoundaryFaceFinal();
    void Normalize();
    void FreeMemory();
    void PrintBoundaryFace();

    bool ElementNeedSplit(int ielement);
    void SplitElementFace(int ielement);

    void GridAdaption();
    void ConstructBasicData();
    void MarkAllCandidateBoundaryFace();

    void SplitElement();
    void SplitElement(int ie);
    void CreateElementGlobal2LocalMap();

    void AddElement(vector<int> &ptlist, int type, int level, int parent);
    void AddChildElement(int ie);
    void HideElement(int ie);

    int  AddElementFace(vector<int> &facenodes, int parent_face_id);

    void ResolveComputationGrid();
private:
    void ResolveComputationStructuredGrid();
    void ResolveComputationUnstructuredGrid();
    int  GetGridType() const; 
};

class CGNSGridFactory
{
public:
    CGNSGridFactory();
    ~CGNSGridFactory();
private:
    int nzones;
    SuperElement **super_elements;
    CGNSMultiZone *cgns_multizone;
public:
    void ConstructGrid();
    void ReadGrid();
    void ResolveComputationGrid();
    void CreateGeneralMultiBlockGrid();
    int  GetGridType(int iZone) const;
    Grid * CreateBasicGrid(int iZone);
private:
    void CreateAllSingleGrid(Grid ** grids);
};

/*PointOperation::point_type GetCenterPoint(vector<int> &elem, vector<int> &point_list);*/

}
