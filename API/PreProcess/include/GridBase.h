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
//! @file      GridBase.h
//! @brief     Explain this file briefly.
//! @author    He Xin.

#pragma once
#include <vector>
#include <set>
#include <map>
#include <string>
#include "Precision.h"
#include "cgnslib.h"
#include "Geo_UnstructGrid.h"

using namespace std;

namespace PHSPACE
{

class CGNSRawCoor
{
public:
	void *GetCoor(int icoor) { return coor[icoor]; };
private:
    cgsize_t nTotalNode;
	DataType_t *type;
	void *coor[3];
public:
	CGNSRawCoor();
	~CGNSRawCoor();
public:
	void SetAllData(RDouble *x, RDouble *y, RDouble *z);
	void AllocateData(int icoor, cgsize_t nTotalNode, DataType_t type);
private:
	void SetData(int icoor, DataType_t type, RDouble *x);
	void DeAllocateData();
};

class RawGrid
{
public:
	RawGrid();
	~RawGrid();
private:
    cgsize_t nTotalNode;
	RDouble *x, *y, *z;
public:
	cgsize_t GetNTotalNode() const { return nTotalNode; }
	void SetNTotalNode(cgsize_t nTotalNode) { this->nTotalNode = nTotalNode; }
	RDouble *GetX() const { return x; }
	RDouble *GetY() const { return y; }
	RDouble *GetZ() const { return z; }
	void SetX(RDouble *x) { this->x = x; };
	void SetY(RDouble *y) { this->y = y; };
	void SetZ(RDouble *z) { this->z = z; };

};

class MultiIndex
{
public:
    int size, id;
    int *data;
public:
    MultiIndex();
    MultiIndex(int size, int id = 0);
    MultiIndex(const MultiIndex & rhs);
    MultiIndex & operator = (const MultiIndex & rhs);
    ~MultiIndex();
public:
    bool operator < (const MultiIndex &rhs) const;
	void SetData(int *data);
};

class CElemFace
{
public:
    CElemFace();
    ~CElemFace();
private:
    int face_num;
    int *face_clkwise;
    int *facept_num;
    int **face_id;
public:
	int **GetFaceID() const { return face_id; }
	int * GetFacePointNumber() const { return facept_num; }
	int   GetFaceNumber() const { return face_num; }
public:
    void Init(int elem_type);
};

class BaseElement
{
public:
    BaseElement();
    ~BaseElement();
private:
    int size;
	CElemFace *elementface;
    cgsize_t *elementpoint;
public:
    cgsize_t *GetElementPoint() const { return elementpoint; }
	CElemFace *GetElementFase() const { return elementface; }
private:
    void InitElementFace();
    void InitElementPoint();
};

class BaseBCType
{
public:
    BaseBCType();
    ~BaseBCType();
	void Init();
private:
	int ibcmax;
    map<int,int> bctypemap;
	map<int,int> fts2cgns;
public:
	int GetBCType(int bctype) { return bctypemap[bctype]; }
	int Fantasy2CGNS(int bctype) { return fts2cgns[bctype]; }
};

class Base_Grid_Conn;

class CGNSBase
{
public:
	CGNSBase();
	~CGNSBase();
public:
	void SetNTotalCell(cgsize_t nTotalCell) { this->nTotalCell = nTotalCell; }
	cgsize_t  GetNTotalCell() const { return nTotalCell; }

	void SetNTotalNode(cgsize_t nTotalNode) { this->nTotalNode = nTotalNode; }
	cgsize_t  GetNTotalNode() const { return nTotalNode; }

	void SetNBoundFace(cgsize_t nBoundFace) { this->nBoundFace = nBoundFace; }
	cgsize_t  GetNBoundFace() const { return nBoundFace; }
private:
	cgsize_t nTotalNode, nTotalCell, nBoundFace;
	int nSections, nBCRegions;
	string* boundaryName;
	cgsize_t *istart, *iend;
	cgsize_t *ElementDataSize;
	int *elem_type;
	cgsize_t **conn_list;
	int *nBCElem;
	cgsize_t **mixed_conn_index;

	int *bc_type;
	int *bc_elem_set_type, *bc_grid_location;
	cgsize_t **bc_conn_list;
public:
	string* GetBCName() const { return boundaryName; }
	void    SetBCName(string* nameIn) { boundaryName = nameIn; } 

	cgsize_t ** GetBoundaryElementConnectionList() const { return bc_conn_list; }
	void   SetBoundaryElementConnectionList(cgsize_t **bc_conn_list) { this->bc_conn_list = bc_conn_list; }

	int * GetBCType() const { return bc_type; };
	void  SetBCType(int *bc_type) { this->bc_type = bc_type; };

	int * GetBoundaryGridLocation() const { return bc_grid_location; };
	void  SetBoundaryGridLocation(int *bc_grid_location) { this->bc_grid_location = bc_grid_location; };

	int * GetBoundaryElementType() const { return bc_elem_set_type; };
	void  SetBoundaryElementType(int *bc_elem_set_type) { this->bc_elem_set_type = bc_elem_set_type; };

	int * GetNumberOfBCElements() const { return nBCElem; };
	void  SetNumberOfBCElements(int *nBCElem) { this->nBCElem = nBCElem; };

	cgsize_t ** GetElementConnectionList() const { return conn_list; };
	cgsize_t ** GetMixedConnectionIndex() const { return mixed_conn_index; };
	int  GetNumberOfSections() const { return nSections; };
	void SetNumberOfSections(int nSections) { this->nSections = nSections; };

	int  GetNumberOfBCRegions() const { return nBCRegions; };
	void SetNumberOfBCRegions(int nBCRegions) { this->nBCRegions = nBCRegions; };

	cgsize_t  * GetIndexOfStart() const { return istart; };
	cgsize_t  * GetIndexOfEnd() const { return iend; };
	cgsize_t  * GetElementDataSize() const { return ElementDataSize; };

	int  * GetElementType() const { return elem_type; };
public:
	void CreateElement(int nSections);
	void CreateBC(int nBCRegions);
	void PreProcess(Base_Grid_Conn *gconn, BaseElement *base_elem);
	void PostProcess(Base_Grid_Conn *gconn);
    void GetISection(cgsize_t icell, int &isection);

    bool CheckBCFace(int *bc_face, int nfacept, set<cgsize_t> &bc_vertex);

    void ComputeGridConn(Base_Grid_Conn *gconn, BaseElement *base_elem);
    void GetElem(cgsize_t *elem_source, cgsize_t icell, int ptnum, cgsize_t *&elem);
	void GetBCFaceMixed(cgsize_t *elem_source, cgsize_t icell, int isection, cgsize_t *&elem);
    void GetElem(BaseElement *base_elem, cgsize_t icell, int &etype, cgsize_t *&elem);
    void GetFace(BaseElement *base_elem, cgsize_t icell, int &facept_num, cgsize_t *&bcface);
	void GetFace(Base_Grid_Conn *gconn, cgsize_t *elem, cgsize_t icell, CElemFace *elemface);
    void GetFace(cgsize_t *elem, int loc, int *bc_face, int &nfacept, CElemFace *elemface);

	void ComputeElemFace(Base_Grid_Conn *gconn, BaseElement *base_elem);
    void ComputeBCFace  (Base_Grid_Conn *gconn, BaseElement *base_elem);
	void ComputeCellToNode(Base_Grid_Conn *gconn, BaseElement *base_elem);
    bool MarkBCFace(Base_Grid_Conn *gconn, cgsize_t *bcface, int facept_num, BCType_t bctype, const string & BCName);
    void MarkBCInterface(Base_Grid_Conn *gconn);

    void GetBCFace(BaseElement *base_elem, cgsize_t icell, int &facept_num, cgsize_t *&bcface);
    void ProcessBCFace  (Base_Grid_Conn *gconn, BaseElement *base_elem, 
		cgsize_t icell, BCType_t bctype, const string & BCName);
    void ProcessBCVertex(Base_Grid_Conn *gconn, BaseElement *base_elem, 
		int bc_elem_set_type, int nBCElem, cgsize_t *bc_conn_list, int bctype, const string & BCName);
    void ProcessBCVertex(Base_Grid_Conn *gconn, BaseElement *base_elem, 
		set<cgsize_t> &bc_vertex, int bctype, const string & BCName);
};

class CGNS_BCFace
{
public:
	CGNS_BCFace(int faceID, int bcType, const string & bcName)
	{
		this->faceID = faceID;
		this->bcType = bcType;
		this->bcName = bcName;
	}

private:
	int faceID;
	int bcType;
	string bcName;
	//int bcID;

public:
	int GetFaceID() const {return this->faceID;}
	int GetBCType() const {return this->bcType;}
	const string & GetBCName() const {return this->bcName;}

	void SetFaceID(int faceID) {this->faceID = faceID;}
	void SetBCType(int bcType) {this->bcType = bcType;}
	void SetBCName(const string & bcName) {this->bcName = bcName;}

	CGNS_BCFace & operator = (const CGNS_BCFace & rhs)
	{
		if ( this == &rhs ) return *this;

		this->faceID = rhs.GetFaceID();
		this->bcType = rhs.GetBCType();
		this->bcName = rhs.GetBCName();

		return *this;
	}
};

class Base_Grid_Conn
{
public:
	Base_Grid_Conn();
	~Base_Grid_Conn();
public:
	void SetNTotalFace(cgsize_t nTotalFace) { this->nTotalFace = nTotalFace; }
	cgsize_t  GetNTotalFace() const { return nTotalFace; }

	void SetNBoundFace(cgsize_t nBoundFace) { this->nBoundFace = nBoundFace; }
	cgsize_t  GetNBoundFace() const { return nBoundFace; }

	void SetNTotalCell(cgsize_t nTotalCell) { this->nTotalCell = nTotalCell; }
	cgsize_t  GetNTotalCell() const { return nTotalCell; }

	void SetNTotalNode(cgsize_t nTotalNode) { this->nTotalNode = nTotalNode; }
	cgsize_t  GetNTotalNode() const { return nTotalNode; }

	vector< vector<cgsize_t> > * GetFace2Node() const { return face2node; }
	vector<CGNS_BCFace*> * GetBCFace() const { return this->bcFace; }
	vector<cgsize_t> * GetLeftCellOfFace() const { return left_cell_of_face; };
	vector<cgsize_t> * GetRightCellOfFace() const { return right_cell_of_face; };
	vector<int> * GetFaceLocation() const { return face_location; };
	vector<int> * GetNodeNumberOfEachFace() const { return node_number_of_each_face; }
	set<MultiIndex> * GetFaceSet() const { return faceset; } ;
	vector<int> * GetNodeNumberOfEachCell() const { return node_number_of_each_cell; }
	vector< vector<cgsize_t> > * GetCell2Node() const { return cell2node; }

	void FreeFaceSet();
	void FreeFace2Node();
	void FreeCell2Node();
	void FreeCellOfFace();
	void FreeNodeNumberOfEachFace();
	void FreeNodeNumberOfEachCell();

private:
	cgsize_t nTotalNode, nTotalCell, nBoundFace, nTotalFace;
	set<MultiIndex> * faceset;
	vector<int> *face_location;
	vector<CGNS_BCFace *> * bcFace;
	vector<cgsize_t> *left_cell_of_face, *right_cell_of_face;
	vector<int> *node_number_of_each_face;
	vector< vector<cgsize_t> > *face2node;
	vector<int> *node_number_of_each_cell;
	vector< vector<cgsize_t> > *cell2node;
};

class CGNSFactory
{
private:
    int nBlocks;
    CGNSBase       **base_cgns;
    Base_Grid_Conn **base_conn;
	BaseElement    * base_elem;
	BaseBCType     * fantasy_bc;
	RawGrid        **rawgrid;
public:
	CGNSFactory(int nBlocks);
    ~CGNSFactory();
public:
	Base_Grid_Conn **GetBaseGridConnection() const { return base_conn; };
	Base_Grid_Conn *GetBaseGridConnection(int izone) const { return base_conn[izone]; };
	BaseElement    *GetBaseElement() const { return base_elem; };
	BaseBCType     *GetBaseBCType() const { return fantasy_bc; };
	RawGrid  * GetRawGrid(int izone) { return rawgrid[izone]; } 
	CGNSBase * GetCGNSBase(int izone = 0) { return base_cgns[izone]; };
	int GetNumberofBlocks() const { return nBlocks; }
    void ConvertGrid2Fantasy();
	void DumpSelfDefinePartInformation();
	void SetNumberofBlocks(int nBlocks_in) { this->nBlocks = nBlocks_in; };
private:
	void ComputeGridConn(CGNSBase *cgnsgrid, Base_Grid_Conn *gconn);
	void ReLableSelfDefineBC();
};

void CGNS2UnsGrid  (CGNSFactory *factory_cgns, int izone, UnstructGrid *grid);

}
