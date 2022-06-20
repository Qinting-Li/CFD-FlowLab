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
//! @file      PBObject.h
//! @brief     The base class for all PBSolver objects.
//! @author    Wan Yunbo, Ming Pingjian, Meng Liyuan, Sun Huawen, Bell.

#pragma once

#include <map>
#include <string>
#include <sstream>
#include <iostream>
#include <vector>

#include "GlobalDataBase.h"
#include "Constants.h"
#include "Precision.h"
#include "Geo_Grid.h"


#ifndef PBOBJECT
#define PBOBJECT

namespace PHSPACE
{
#define _STR(x) #x
#define STR(x) _STR(x)

//! @brief The base class for all PBSolver objects.
class PBObject
{
public:
    PBObject(const std::string &name) : name(name) {}
    virtual ~PBObject() {}

public:
    //! @brief Returns the object name¡£
    //! @return name¡£
    inline const std::string& Name() { return name; }

protected:
    const std::string name_; //!object name¡£

    //! Get one string type data from database by name.
    std::string ObjGetStrParaFromDB(const std::string &name, std::string defaultValue)
    {
        if (GlobalDataBase::GetDataPara()->IsExist(name, PHSTRING, 1))
        {
            return(GlobalDataBase::GetStrParaFromDB(name));
        }
        else
        {
            return(defaultValue);
        }
    }

    //! Get one integer type data from database by name.
    int ObjGetIntParaFromDB(const std::string &name, int defaultValue)
    {
        if (GlobalDataBase::GetDataPara()->IsExist(name, PHINT, 1))
        {
            return(GlobalDataBase::GetIntParaFromDB(name));
        }
        else
        {
            return(defaultValue);
        }
    }

    //! Get one double type data from database by name.
    double ObjGetDoubleParaFromDB(const std::string &name, double defaultValue)
    {
        if (GlobalDataBase::GetDataPara()->IsExist(name, PHDOUBLE, 1))
        {
            return(GlobalDataBase::GetDoubleParaFromDB(name));
        }
        else
        {
            return(defaultValue);
        }
    }

    //! @brief Allocate variable memory.
    //! @param name 
    //! @param size 
    //! @param grid 
    //! @return param pointers
    inline RDouble *ObjAllocateMemory(const std::string &name, const int size, const RDouble init = 0.0, Grid *grid = NULL)
    {
        if (grid != NULL)
        {
            if (grid->GetDataPtr(name))
            {
                return((RDouble*)grid->GetDataPtr(name));
            }
        }

        RDouble *phi = new RDouble[size];

        if (grid != NULL)
        {
            grid->UpdateDataPtr(name, phi);
        }

        for (int i = 0; i < size; ++ i)
        {
            phi[i] = init;
        }

        return(phi);
    }

private:
    std::string name; //! Object name.
};

#define DEFINE_PBOBJECT_FACTORY(ObjectTypeName) \
typedef ObjectTypeName *(* Creator##ObjectTypeName)(const std::string &name); \
typedef std::map<std::string, Creator##ObjectTypeName> CreatorRegistry##ObjectTypeName; \
class ObjectTypeName##Factory \
{ \
public: \
    static CreatorRegistry##ObjectTypeName &Registry() \
    { \
        static CreatorRegistry##ObjectTypeName *registry = new CreatorRegistry##ObjectTypeName(); \
        return *registry; \
    } \
    static void AddCreator(const std::string &type, Creator##ObjectTypeName creator) \
    { \
        CreatorRegistry##ObjectTypeName &registry = Registry(); \
        if(registry.count(type) == 0) \
        { \
            registry[type] = creator; \
        } \
        else \
        { \
            std::cerr << type << " is already registered." << std::endl; \
        } \
    } \
    static ObjectTypeName *Create##ObjectTypeName(const std::string &type) \
    { \
        CreatorRegistry##ObjectTypeName &registry = Registry(); \
        if(registry.count(type) == 0) \
        { \
            std::cerr << "Unknown " << std::string(#ObjectTypeName) << " type: " << type << std::endl; \
            return(NULL); \
        } \
        else \
        { \
            return(registry[type](type)); \
        } \
    } \
private: \
    ObjectTypeName##Factory(){} \
}; \
class ObjectTypeName##Register \
{ \
public: \
    ObjectTypeName##Register(const std::string &type, Creator##ObjectTypeName creator) \
    { \
        ObjectTypeName##Factory::AddCreator(type, creator); \
    } \
};

#define _DEFINE_PBSOLVER_REGISTER_OBJECT(ObjectTypeName, DerivedTypeName, type) \
static ObjectTypeName##Register register_##ObjectTypeName##type(#type, Creator_##ObjectTypeName##type);

#define _DECLARE_PBSOLVER_REGISTER_FUNCTION(ObjectTypeName, DerivedTypeName, type) \
ObjectTypeName *Creator_##ObjectTypeName##type(const std::string &name);

#define _DEFINE_PBSOLVER_REGISTER_FUNCTION(ObjectTypeName, DerivedTypeName, type) \
ObjectTypeName *Creator_##ObjectTypeName##type(const std::string &name) \
{ \
    return(new DerivedTypeName(name));\
}

#define DEFINE_PBSOLVER_REGISTER_OBJECT(ObjectTypeName, DerivedTypeName, type) \
_DEFINE_PBSOLVER_REGISTER_OBJECT(ObjectTypeName, DerivedTypeName, type)

#define DECLARE_PBSOLVER_REGISTER_FUNCTION(ObjectTypeName, DerivedTypeName, type) \
_DECLARE_PBSOLVER_REGISTER_FUNCTION(ObjectTypeName, DerivedTypeName, type)

#define DEFINE_PBSOLVER_REGISTER_FUNCTION(ObjectTypeName, DerivedTypeName, type) \
_DEFINE_PBSOLVER_REGISTER_FUNCTION(ObjectTypeName, DerivedTypeName, type)
}

#endif 
