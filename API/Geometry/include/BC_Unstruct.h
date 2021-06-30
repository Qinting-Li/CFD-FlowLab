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
//! @file      BC_Unstruct.h
//! @brief     It defines unstructured grid/solver Boundary Condition.
//! @author    Bell, He Xin.

#pragma once
#include "BC.h"
using namespace std;

namespace PHSPACE
{
class BoundaryCondition;

//class BCRegionUnstruct is used to store a single boundary Region.
//i.e. including some boundary faces with a same name.
class BCRegionUnstruct
{
private:
    //! The bc type, used to store an integer type.
    int bcType;

    //! The bcName store the name of this boundary region.
    string bcName;

    //! The bcData store some variables used in this boundary region.
    Data_Param * bcData;

    //! The faceIndex set store the all face numbers in this boundary region.
    vector <int> * faceIndex;
public:
    void SetBCType(int bcType) { this->bcType = bcType; }
    int GetBCType() const { return bcType; };    //! Return the bcType.

    const string & GetBCName() const { return this->bcName; }
    void SetBCName(const string &bcName) { this->bcName = bcName; }

    Data_Param * GetBCDataParam() const { return this->bcData; }
    void SetBCDataParam(Data_Param *bcData) { this->bcData = bcData; }

    vector <int> * GetFaceIndex() const { return this->faceIndex; }
    void SetFaceIndex(vector<int> *faceIndex) { this->faceIndex = faceIndex; }

public:
    //! copy constructor function of this class.
    BCRegionUnstruct(const BCRegionUnstruct &bcRegionUnstruct)
    {
        bcData = bcRegionUnstruct.GetBCDataParam();
        bcType = bcRegionUnstruct.GetBCType();
        bcName = bcRegionUnstruct.GetBCName();
    }
    //! constructor function of this class.
    BCRegionUnstruct()
    {
        bcData = 0;    //! new Data_Param();
        faceIndex = new vector<int>;
        bcType = 0;
    }
    //! destructor function of this class.
    ~BCRegionUnstruct()
    {
        //delete bcData;

        faceIndex->clear();
        delete faceIndex;
    }
};

class BCRecord
{
private:
    int key;    //! The bc key, used to store an integer type.
    string bcName;
    BoundaryCondition * boundaryCondition;
    BCRegionUnstruct * bcRegion;
public:
    void SetKey(int key) { this->key = key; }
    int GetKey() const { return key; }    //! Return the key.
    const string & GetBCName() const { return this->bcName; }
    void SetBCName(const string &bcName) { this->bcName = bcName; }
    void SetBCRegion(BCRegionUnstruct *bcRegionIn) { this->bcRegion = bcRegionIn; }
    BCRegionUnstruct * GetBCRegion() const { return this->bcRegion; }
    BoundaryCondition * GetBoundaryCondition() const;
    void SetBoundaryCondition(BoundaryCondition *boundaryCondition) { this->boundaryCondition = boundaryCondition; }
public:
    BCRecord(const BCRecord &rhs)
    {
        key = rhs.key;
        bcName = rhs.GetBCName();
        this->boundaryCondition = rhs.GetBoundaryCondition();
    }
    BCRecord() { key = 0; boundaryCondition = 0; bcRegion = 0; }
   ~BCRecord() {}
};

}