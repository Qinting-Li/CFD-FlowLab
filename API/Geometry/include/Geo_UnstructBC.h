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
//! @file      Geo_UnstructBC.h
//! @brief     It is the base class of geometry operation for boundary condition.
//!            The inheriting order is: SimpleBC -> StructuredBC/UnstructuredBC.
//! @author    Zhang Yong, Bell, He Xin.

#pragma once
#include "Geo_SimpleBC.h"
using namespace std;

namespace PHSPACE
{
class SimpleBC;

//class BCRegionUnstruct is used to store a single boundary Region.
//i.e. including some boundary faces with a same name.
class UnstructBC : public SimpleBC
{
private:
    //! The faceIndex set store the all face numbers in this boundary region.
    vector <int> * faceIndex;
    int zoneID;
    int bcRegionID;

public:
    vector <int> * GetFaceIndex() const { return this->faceIndex; }
    void SetFaceIndex(vector<int> *faceIndex) { this->faceIndex = faceIndex; }

public:
    //! copy constructor function of this class.
    UnstructBC(UnstructBC &bcRegionUnstruct)
    {
        SetBCParamDataBase(bcRegionUnstruct.GetBCParamDataBase());
        SetBCType(bcRegionUnstruct.GetBCType());
        SetBCName(bcRegionUnstruct.GetBCName());
    }
    //! constructor function of this class.
    UnstructBC()
    {
        faceIndex = new vector<int>;
    }
    //! destructor function of this class.

    UnstructBC(int zoneID, int bcRegionID)
    {
        this->zoneID = zoneID;
        this->bcRegionID = bcRegionID;
        faceIndex = new vector<int>;
    }

    ~UnstructBC()
    {
        faceIndex->clear();
        delete faceIndex;
    }
};

class CompositeUnstructBC
{
private:
    int key;    //! The bc key, used to store an integer type.
    string bcName;
    SimpleBC * boundaryCondition;
    UnstructBC * bcRegion;
    
    //! Number of zone (block) index this BC in.
    int zoneID;
    
    //! The BCRegions list in the current CompositeUnstructBC.
    vector<UnstructBC *> * bcRegions;

    //! The bcRegions id of boundary condition faces.
    int *bcRegionIDofBCFace; 

public:
    void SetKey(int key) { this->key = key; }
    int GetKey() const { return key; }    //! Return the key.
    const string & GetBCName() const { return this->bcName; }
    void SetBCName(const string &bcName) { this->bcName = bcName; }
    void SetBCRegion(uint_t iBCRegion, UnstructBC *unstructBC) { (*bcRegions)[iBCRegion] = unstructBC; }
    
    UnstructBC * GetBCRegion(uint_t iBCRegion) const{ return (*bcRegions)[iBCRegion];}

    SimpleBC * GetBoundaryCondition() const;
    void SetBoundaryCondition(SimpleBC *boundaryCondition) { this->boundaryCondition = boundaryCondition; }

    void CreatenBCRegion(uint_t nBCRegion) {bcRegions = new vector<UnstructBC *>(nBCRegion);}
    int GetnBCRegion() const {  return static_cast<int>(bcRegions->size()); }

    void SetBCFaceInfo(int* bcRegionIDofBCFace) {this->bcRegionIDofBCFace = bcRegionIDofBCFace;}

    //! Get the bcregion id by interfaceIndexContainerForReceive value.
    int * GetBCFaceInfo() const {  return this->bcRegionIDofBCFace; }
public:
    CompositeUnstructBC(const CompositeUnstructBC &rhs)
    {
        key = rhs.key;
        bcName = rhs.GetBCName();
        this->boundaryCondition = rhs.GetBoundaryCondition();
    }
    //CompositeUnstructBC() { key = 0; boundaryCondition = 0; bcRegion = 0; }
    CompositeUnstructBC(int zoneID = 0) { this->zoneID = zoneID; bcRegionIDofBCFace = NULL; key = 0; boundaryCondition = 0; bcRegion = 0; }
  
    ~CompositeUnstructBC()
   {
       if (bcRegions)
	   {
		   for (std::size_t i = 0; i < bcRegions->size(); ++i)
		   {
			   delete (*bcRegions)[i];
		   }
	   }
       delete bcRegions;    bcRegions = NULL;

       if (bcRegionIDofBCFace != NULL)
       {
           delete bcRegionIDofBCFace;
           bcRegionIDofBCFace = NULL;
       } 
   }
};

}