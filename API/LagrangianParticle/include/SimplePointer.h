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
//! @file      SimplePointer.h
//! @brief     It is the basic pointer structure,
//!            witch is similar to Pointer.h and SimpleGrid.h
//! @author    Lei Yinghaonan (Lanzhou University).

#pragma once
#include "Precision.h"
#include "TypeDefine.h"
#include "Data_Util.h"

namespace PHSPACE
{

template <class T >
class SimplePointer
{
private:
    T *data;
public:
    //! Default Construct function.
    //! The pointer is initialized to 0.
    SimplePointer();

    //! Construct function.
    //! New a size T *data.
    SimplePointer(const int &size);

    //! Destructor function.
    //! If new , delete data and set null.
    ~SimplePointer();

public:
    //! The function for set data size.
    void SetNull();
    void SetSize(const int &size);
    void SetSize(const int &size, T &data);
    void SetResize(const int &size);

    //! Make the addresses of two data the same directly.
    //! If dataIn is an array, the first address of pointer is the same.
    void SetData(const SimplePointer &dataIn);

    //! Copy the value of dataIn from 0 to size.
    void CopyData(SimplePointer &dataIn,const int &size);
    //! Copy the value of dataIn from 0 to size.
    void CopyData(const T dataIn, const int &size);
    //! Copy the value of dataIn by index-th.
    void CopyData(const T *dataIn, int &indexSt , int &indexEd);
    //! Copy only one value of dataIn by indexDataIn-th on data's indexData-th.
    void CopyDataOnIndex(const T *dataIn, int &indexDataIn,int &indexData);

    //! The data value is determined by a single value.
    void InitData(const T &value,const int &size);
    //! The data value is determined by a single value.
    //! Set this data value from indexSt-th to indexEd-th.
    void InitData(const T &value, int &indexSt, int &indexEd);

    void PrintData2Window(const int &size);
    void PrintData2Window(int indexSt, int indexEd);

    void DotValue(const T &value, const int &size);

    T *GetPointer();
    T *GetPointer(int &index);

    T &GetData();
    T &GetData(int &index);

    T &operator[](int index)
    {
        return (data[index]);
    }

    T *operator()(int index)
    {
        return *(data[index]);
    }

private:
    void DeleteData();

};

typedef SimplePointer<int> SPint;
typedef SimplePointer<RDouble> SPDouble;

template <class T>
SimplePointer<T>::SimplePointer()
{
    SetNull();
}

template <class T>
SimplePointer<T>::SimplePointer(const int &size)
{
    SetSize(size);
}

template <class T>
SimplePointer<T>::~SimplePointer()
{
    this->DeleteData();
}

template <class T>
void SimplePointer< T >::SetNull()
{
    //! When C++ does not explicitly initialize non static local variables or member variables, 
    //! their values are random (garbage values)
    //! For static or global variables, 
    //! the initial value will be 0 (built-in value type or pointer) or false (bool). 
    //! For the class object here, the internal is still the same as above
    //! So here, please note that,
    //! if we doesn't specify the pointer to an address, 
    //! the value of pointer is not zeros.
    //! by Lei Yinghaonan, 2021/7/27
    this->data = NULL;
}

template <class T>
void SimplePointer< T >::SetSize(const int &size)
{
    if (this->data != 0)
    {
        SetNull();
    }
    this->data = new T[size];
}

template <class T>
void  SimplePointer< T >::SetSize(const int &size, T &data)
{
    if (this->data != 0)
    {
        SetNull();
    }
    this->data = new T[size];
    for (int iDim = 0; iDim < size; ++iDim)
    {
        this->data[iDim] = data;
    }
}

template <class T>
void SimplePointer< T >::SetResize(const int  &size)
{
    if (this->data != 0)
    {
        DeleteData();
    }
    this->data = new T[size];
}

template <class T>
void SimplePointer< T >::SetData(const SimplePointer &dataIn)
{
    this->data = dataIn->data;
}

template <class T>
void SimplePointer< T >::CopyData(SimplePointer &dataIn,const int &size)
{
    for (int iter = 0; iter < size; ++iter)
    {
        this->data[iter] = dataIn.GetData(iter);
    }
}

template <class T>
void SimplePointer< T >::CopyData(const T dataIn, const int &size)
{
    for (int iter = 0; iter < size; ++iter)
    {
        data[iter] = dataIn;
    }
}

template <class T>
void SimplePointer< T >::CopyData(const T *dataIn,  int &indexSt , int &indexEd)
{
    for (int iter = indexSt; iter < indexEd; ++iter)
    {
        data[iter] = dataIn[iter];
    }
}

template <class T>
void SimplePointer< T >::CopyDataOnIndex(const T *dataIn, int &indexDataIn, int &indexData)
{
    this->data[indexData] = dataIn[indexDataIn];
}

template <class T>
void SimplePointer< T >::InitData(const T &value,const int &size)
{
    for (int iter = 0; iter < size; ++iter)
    {
        data[iter] = value;
    }
}

template <class T>
void SimplePointer< T >::InitData(const T &value, int &indexSt , int &indexEd)
{
    for (int iter = indexSt; iter < indexEd; ++iter)
    {
        this->data[iter] = value;
    }
}

template <class T>
void SimplePointer<T>::PrintData2Window(const int &size)
{
    for (int iter = 0; iter < size; ++iter)
    {
        cout << " iter = " << iter << " data = " << data[iter] << endl;
    }
}

template <class T>
void SimplePointer<T>::PrintData2Window(int indexSt , int indexEd)
{
    for (int iter = indexSt; iter <= indexEd; ++iter)
    {
        cout << " iter = " << iter << " data = " << data[iter] << endl;
    }
}

template <class T>
T *SimplePointer< T >::GetPointer()
{
    return data;
}

template <class T>
T *SimplePointer< T >::GetPointer(int &index)
{
    return &(data[index]);
}

template <class T>
T &SimplePointer< T >::GetData()
{
    return data[0];
}

template <class T>
T &SimplePointer< T >::GetData(int &index)
{
    return data[index];
}

template <class T>
void SimplePointer<T>::DeleteData()
{
    if (this->data == 0)
    {
        return;
    }
    else
    {
        Data_Util::delete_void<T>(this->data);
        SetNull();
    }
}

template <class T>
void SimplePointer<T>::DotValue(const T &value, const int &size)
{
    for (int iter = 0; iter < size; ++iter)
    {
        data[iter] = data[iter]  *value;
    }
}

}