//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//          PPPPP  H   H  EEEEE  N    N  GGGGG  L      EEEEE  III         +
//          P   P  H   H  E      NN   N  G      L      E       I          +
//          PPPPP  HHHHH  EEEEE  N N  N  G  GG  L      EEEEE   I          +
//          P      H   H  E      N  N N  G   G  L      E       I          +
//          P      H   H  EEEEE  N    N  GGGGG  LLLLL  EEEEE  III         +
//------------------------------------------------------------------------+
//          Platform for Hybrid Engineering Simulation of Flows           +
//          China Aerodynamics Research and Development Center             +
//                     (C) Copyright, Since 2010                          +
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      LagrangianVariable.h
//! @brief     It is the class of lagrangian varibale.
//! @author    Lei Yinghaonan (Lanzhou University).

#pragma once
#include "IndexCell.h"
#include "TK_Exit.h"

using namespace std;

namespace PHSPACE
{

template <class T>
class LagrangianVariable
{
private:
    //! The lagrangian variable map.
    //! The first varibal is the global id for particle for all zones.
    //! The second varibal is each lagrangian point variable on this zones.
    map<int, T* > *lagrangianVariable;

    //! The mapping relationship of cell local id and lagrangian global id .
    IndexCell *cellLagrangianID;

    //! The size of T*.
    int nDimVariable;

public:
    //! Constructor function, new map and but don't IndexCell.
    LagrangianVariable();

    //! Constructor function, new map and copy IndexCell's address.
    LagrangianVariable(IndexCell *rhs);

    //! Destructor function , delete the map.
    ~LagrangianVariable();
    
public:
    //! New the elements of map and add to map.
    //! Note this function 's map key is started from 0 t0 nVar-1;
    //! So this function does not save the local varID,
    //! which is very simple when only serial computing.
    void InitLagrangianVariable(int &nVar, const int &nDim);
    //! New the elements of map and add to map.
    //! Note this function's map key is started by int *varIDList;
    void InitLagrangianVariable(int *varIDList,int &nVar, const int &nDim);
    //! New the elements of map and add to map.
    //! copy the value of map by a other LagrangianVariabl.
    void InitLagrangianVariable(LagrangianVariable<T> *dataIn);

    //! Copy the value by only one value for all variableID on all iDim.
    void CopyValue(const T &value);
    //! Copy the value by only one value for all variableID on iDim.
    void CopyValue(const T &value, const int &iDim);
    //! Copy tthe value by only one value for variableID all iDim.
    void CopyValue(int &variableID, const T &value);
    //! Copy tthe value by only one value for variableID on iDim.
    void CopyValue(int &variableID, const T &value, const int &iDim);
    //! Copy the value of var on all variableID on all iDim.
    void CopyValue(T *var);
    //! Copy the iDim-th value of var[iDim] on all variableID on iDim.
    void CopyValue(T *var, const int &iDim);
    //! Copy the value of var by variableID on all iDim.
    void CopyValue(int &variableID, T *var);
    //! Copy the iDim-th value of var[iDim] by variableID on iDim.
    void CopyValue(int &variableID, T *var, const int &iDim);

    //! Set the address of current cellLagrangianID.
    void SetIndexCell(IndexCell *rhs);
    //! SetNumOfDimOfVariable.
    void SetNumOfDim(int &nDim);

    //! Get the init.
    int GetInit();
    //! Get the number of lagrangian point;
    int GetNumOfLagrangianPoint();
    //! Get the dimemsion of lagrangian variable.
    int GetDimOfVariable();
    //! Get the value of data by particle id.
    int GetValue(int &varID,int &iDim);
    //! Get the value of data by local index.
    int GetIndexValue(const int iVar, int &iDim);
    //! Get the map of lagrangianVariable.
    map<int, T* > *GetMap();
    //! Get the pointer of current cellLagrangianID.
    IndexCell *GetIndexCell();

    void RemoveParticle(int &particleID);
    void AddParticle(int &particleID);

    //! print the value of varibale by varID.
    void PrintData(int varID);

    T *operator ()(int particleID)
    {
        return (lagrangianVariable->find(particleID))->second;
    }

    T &operator ()(int particleID,int iDim)
    {
        return (lagrangianVariable->find(particleID))->second[iDim];
    }

private:
    //! Add the variable to map.
    //! Note that we should ensure that 
    //! only particlepointgroup has permission to modify indexcell.
    void AddVariable2Map(int &variableID, T *lagrangianVariable);

    //! Delete the map of lagrangianVariable.
    void DeleteMap();

    //! Delete data.
    void DeleteData();

};

typedef LagrangianVariable<int> Lint;
typedef LagrangianVariable<RDouble> LDouble;

template <class T>
LagrangianVariable<T>::LagrangianVariable()
{
    //! Note that most of the indexcells are indexcells of particlepointgroup, 
    //! so there is no room for celllagrangianid here. 
    //! If it is an independent indexcell, you need a separate new to open it.
    this->lagrangianVariable = new map<int, T* >;
    this->cellLagrangianID = NULL;
    this->nDimVariable = 0;
}

template <class T>
LagrangianVariable<T>::LagrangianVariable(IndexCell *rhs)
{
    this->lagrangianVariable = new map<int, T* >;
    this->SetIndexCell(rhs);
    this->nDimVariable = 0;
}

template <class T>
LagrangianVariable<T>::~LagrangianVariable()
{
    this->DeleteData();
}

template <class T>
void LagrangianVariable<T>::InitLagrangianVariable(int &nVar, const int &nDim)
{
    //! Our purpose here is to assign the value of another map to a new map.
    if (!this->lagrangianVariable->empty())
    {
        ostringstream oss;
        oss << "The map variable is not empty" << endl;
        TK_Exit::ExceptionExit(oss);
    }

    int varID;
    for (int iVar = 0; iVar < nVar; ++iVar)
    {
        varID = iVar;
        //! Parentheses indicate initialization to 0
        T *lagrangianVariable = new T[nDim]();
        this->AddVariable2Map(varID, lagrangianVariable);
    }
    this->nDimVariable = nDim;
}

template <class T>
void LagrangianVariable<T>::InitLagrangianVariable(int *varIDList, int &nVar, const int &nDim)
{
    //! Our purpose here is to assign the value of another map to a new map.
    //if (!this->lagrangianVariable->empty())
    //{
    //    ostringstream oss;
    //    oss << "The map variable is not empty" << endl;
    //    TK_Exit::ExceptionExit(oss);
    //}

    int varID;
    for (int iVar = 0; iVar < nVar; ++iVar)
    {
        varID = varIDList[iVar];
        T *lagrangianVariable = new T[nDim]();
        this->AddVariable2Map(varID, lagrangianVariable);
    }
    this->nDimVariable = nDim;
}

template <class T>
void LagrangianVariable<T>::InitLagrangianVariable(LagrangianVariable<T> *dataIn)
{
    //! Our purpose here is to assign the value of another map to a new map.
    if (!this->lagrangianVariable->empty())
    {
        ostringstream oss;
        oss << "The map variable is not empty" << endl;
        TK_Exit::ExceptionExit(oss);
    }

    this->nDimVariable = dataIn->GetDimOfVariable();

    map<int, T* > *dataInMap = dataIn->GetMap();
    typename map<int, T* >::iterator iter;
    for (iter = dataInMap->begin(); iter != dataInMap->end(); ++iter)
    {
        //! Parentheses indicate initialization to 0
        T *lagrangianVariable = new T[nDimVariable]();
        //! Note here, only copy value.
        for (int iDim = 0; iDim < nDimVariable; ++iDim)
        {
            lagrangianVariable[iDim] = iter->second[iDim];
        }
        int particleID = iter->first;
        this->AddVariable2Map(particleID, lagrangianVariable);
    }
}

template <class T>
void LagrangianVariable<T>::CopyValue(const T &value)
{
    typename map<int, T* >::iterator iter;
    for (iter = lagrangianVariable->begin(); iter != lagrangianVariable->end(); ++iter)
    {
        for (int iDim = 0; iDim < nDimVariable; ++iDim)
        {
            iter->second[iDim] = value;
        }
    }
}

template <class T>
void LagrangianVariable<T>::CopyValue(const T &value, const int &iDim)
{
    if (iDim > this->nDimVariable)
    {
        ostringstream oss;
        oss << "Error: insert variable dimension: " << iDim << endl;
        oss << "current variable dimension: " << this->nDimVariable << endl;
        TK_Exit::ExceptionExit(oss);
    }

    typename map<int, T* >::iterator iter;
    for (iter = lagrangianVariable->begin(); iter != lagrangianVariable->end(); ++iter)
    {
        iter->second[iDim] = value;
    }
}

template <class T>
void LagrangianVariable<T>::CopyValue(int &variableID, const T &value)
{
    typename map<int, T* >::iterator iter;
    iter = this->lagrangianVariable->find(variableID);
    if (iter != lagrangianVariable->end())
    {
        for (int iDim = 0; iDim < nDimVariable; ++iDim)
        {
            iter->second[iDim] = value;
        }
    }
    else
    {
        ostringstream oss;
        oss << "Error: insert variable: " << "variableID = " << variableID << endl;
        TK_Exit::ExceptionExit(oss);
    }
}

template <class T>
void LagrangianVariable<T>::CopyValue(int &variableID, const T &value, const int &iDim)
{
    if (iDim > this->nDimVariable)
    {
        ostringstream oss;
        oss << "Error: insert variable dimension: " << iDim << endl;
        oss << "current variable dimension: " << this->nDimVariable << endl;
        TK_Exit::ExceptionExit(oss);
    }

    typename map<int, T* >::iterator iter;
    iter = this->lagrangianVariable->find(variableID);
    if (iter != lagrangianVariable->end())
    {
        iter->second[iDim] = value;
    }
    else
    {
        ostringstream oss;
        oss << "Error: insert variable: " << "variableID = " << variableID << endl;
        TK_Exit::ExceptionExit(oss);
    }
}

template <class T>
void LagrangianVariable<T>::CopyValue(T *var)
{
    typename map<int, T* >::iterator iter;
    for (iter = lagrangianVariable->begin(); iter != lagrangianVariable->end(); ++iter)
    {
        for (int iDim = 0; iDim < nDimVariable; ++iDim)
        {
            iter->second[iDim] = var[iDim];
        }
    }
}

template <class T>
void LagrangianVariable<T>::CopyValue(T *var, const int &iDim)
{
    if (iDim > this->nDimVariable)
    {
        ostringstream oss;
        oss << "Error: insert variable dimension: " << iDim << endl;
        oss << "current variable dimension: " << this->nDimVariable << endl;
        TK_Exit::ExceptionExit(oss);
    }

    typename map<int, T* >::iterator iter;
    for (iter = lagrangianVariable->begin(); iter != lagrangianVariable->end(); ++iter)
    {
        iter->second[iDim] = var[iDim];
    }
}

template <class T>
void LagrangianVariable<T>::CopyValue(int &variableID, T *var)
{
    typename map<int, T* >::iterator iter;
    iter = this->lagrangianVariable->find(variableID);
    if (iter != lagrangianVariable->end())
    {
        for (int iDim = 0; iDim < nDimVariable; ++iDim)
        {
            iter->second[iDim] = var[iDim];
        }
    }
    else
    {
        ostringstream oss;
        oss << "Error: insert variable: " << "variableID = " << variableID << endl;
        TK_Exit::ExceptionExit(oss);
    }
}

template <class T>
void LagrangianVariable<T>::CopyValue(int &variableID, T *var, const int &iDim)
{
    if (iDim > this->nDimVariable)
    {
        ostringstream oss;
        oss << "Error: insert variable dimension: " << iDim << endl;
        oss << "current variable dimension: " << this->nDimVariable << endl;
        TK_Exit::ExceptionExit(oss);
    }

    typename map<int, T* >::iterator iter;
    iter = this->lagrangianVariable->find(variableID);
    if (iter != lagrangianVariable->end())
    {
        iter->second[iDim] = var[iDim];
    }
    else
    {
        ostringstream oss;
        oss << "Error: insert variable: " << "variableID = " << variableID << endl;
        TK_Exit::ExceptionExit(oss);
    }
}

template <class T>
void LagrangianVariable<T>::SetIndexCell(IndexCell *rhs)
{
    this->cellLagrangianID = rhs;
}

template <class T>
IndexCell *LagrangianVariable<T>::GetIndexCell()
{
    return this->cellLagrangianID;
}

template <class T>
void LagrangianVariable<T>::RemoveParticle(int &particleID)
{
    typename map<int, T* >::iterator iter;
    iter = this->lagrangianVariable->find(particleID);

    //cout << this->lagrangianVariable->size() << endl;
    if (iter != lagrangianVariable->end())
    {
        //cout << iter->second << endl;
        //cout << *(iter->second) << endl;
        if (iter->second != NULL)
        {
            delete iter->second;
            iter->second = NULL;
        }

        lagrangianVariable->erase(iter);
    }
    else
    {
        cout << "Un exist particleID  = " << particleID << endl;
    }
}

template <class T>
void LagrangianVariable<T>::AddParticle(int &particleID)
{
    //! Parentheses indicate initialization to 0
    T *data = new T[nDimVariable]();
    this->AddVariable2Map(particleID, data);
}

template <class T>
void LagrangianVariable<T>::DeleteData()
{
    if (this == 0)
    {
        return;
    }

    //! Destructors can cause many problems of deletion by mistake.
    //! We must be careful.
    this->DeleteMap();

    //! If IndexCell is opened separately by new, 
    //! 'this->DeleteIndexCell()' is required here. 
    //! On the contrary, do not delete indexcell here.
    this->cellLagrangianID = NULL;
    this->nDimVariable = 0;
}

template <class T>
void LagrangianVariable<T>::SetNumOfDim(int &nDim)
{
    this->nDimVariable = nDim;
}

template <class T>
int LagrangianVariable<T>::GetInit()
{
    return this->init;
}

template <class T>
int LagrangianVariable<T>::GetNumOfLagrangianPoint()
{
    return this->lagrangianVariable->size();
}

template <class T>
int LagrangianVariable<T>::GetDimOfVariable()
{
    return this->nDimVariable;
}

template <class T>
int LagrangianVariable<T>::GetValue(int &varID, int &iDim)
{
    typename map<int, T* >::iterator iter;
    iter = this->lagrangianVariable->find(varID);
    return iter->second[iDim];
}

template <class T>
int LagrangianVariable<T>::GetIndexValue(const int iVar, int &iDim)
{
    T value;
    int count = 0;
    typename map<int, T* >::iterator iter;
    for (iter = lagrangianVariable->begin(); iter != lagrangianVariable->end(); ++iter)
    {
        if (count == iVar)
        {
            value = iter->second[iDim];
            break;
        }
        count++;
    }
    return value;
}

template <class T>
map<int, T* > *LagrangianVariable<T>::GetMap()
{
    return this->lagrangianVariable;
}

template <class T>
void LagrangianVariable<T>::PrintData(int varID)
{
    typename map<int, T* >::iterator iter;
    iter = this->lagrangianVariable->find(varID);
    if (iter != lagrangianVariable->end())
    {
        cout << "check variableID  = " << varID << endl;
        for (int iDim = 0; iDim < nDimVariable; ++iDim)
        {
            cout << iter->second[iDim] << endl;
        }
    }
    else
    {
        cout << "Un exist particleID  = " << varID << endl;
    }
}

template <class T>
void LagrangianVariable<T>::AddVariable2Map(int &variableID, T *lagrangianVariable)
{
    bool check = this->lagrangianVariable->insert(make_pair(variableID, lagrangianVariable)).second;
    if (!check)
    {
        ostringstream oss;
        oss << "Error: insert variable: " << "variableID = " << variableID << endl;
        TK_Exit::ExceptionExit(oss);
    }
}

template <class T>
void LagrangianVariable<T>::DeleteMap()
{
    if (this == 0)
    {
        return;
    }

    typename map<int, T *>::iterator iter;
    iter = lagrangianVariable->begin();
    while (iter != lagrangianVariable->end())
    {
        delete [] iter->second;
        iter->second = NULL;
        lagrangianVariable->erase(iter++);
    }

    delete lagrangianVariable;
    lagrangianVariable = NULL;
}

}