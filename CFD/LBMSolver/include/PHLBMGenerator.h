/////////////////////////////////////////////////////////////////////
// to add header.
/////////////////////////////////////////////////////////////////////
// author: Chen, Songze. email: jacksongze@hust.edu.cn
/////////////////////////////////////////////////////////////////////
// participant 1: 
/////////////////////////////////////////////////////////////////////

#pragma once

#ifndef PHLBMGENERATOR_H
#define PHLBMGENERATOR_H
#include <CommonAPI.h>
#include <PHLBM.h>
#include <memory>

namespace phlbm 
{
  /** 利用风雷的文本处理工具读取几何类型数据，生成nnwlbm::Geometry的派生对象
  */
  template<typename Base>
  class Input {
  public:
    /** 返回类型名称. 通过name函数内部定义字符串实现一种类型一个名称，
    * 该名称须于读入的几何类型的name一致
    */
    virtual ::std::string name() const = 0;
  private:
    /** 虚函数实现某种数据类型的二进制数据读入功能 
    * \param ndb [in,out]经风雷导入的数据 
    */
    virtual std::shared_ptr<Base> InputPrivate(
      const phlbm::NamedDataBlock &ndb, // [in,out]生成几何的数据
      const std::map<std::string, phlbm::NamedDataBlock> &ndbMap // [in,out]所有经风雷导入的数据 
    ) const = 0;
  public:
    template<typename Base1>
    friend class Generator;
  };

  template<typename Base, typename Derived>
  class In : public Input<Base> {
  public:
    virtual ::std::string name() const;
  private:
    virtual std::shared_ptr<Base> InputPrivate(
      const phlbm::NamedDataBlock &ndb, // [in,out]生成几何的数据
      const std::map<std::string, phlbm::NamedDataBlock> &ndbMap // [in,out]所有经风雷导入的数据 
    ) const;
  };

  template<typename Base>
  class InputMap {
  public:
    int identity() { return 1; }
    /** 智能指针类型别名 */
    typedef std::shared_ptr<Input<Base>> inptr;


    inline bool emplace(
      const ::std::string &id, /**< [in]读入工具名称 */
      inptr ptr /**< [in]某种数据的读入工具指针 */
    ) 
    {
      return readerMap.emplace(id, ptr).second;
    }

    inline bool replace(
      const ::std::string &id, /**< [in]读入工具名称 */
      inptr ptr /**< [in]某种数据的读入工具指针 */
    ) 
    { 
      bool exist = !(readerMap.emplace(id, ptr).second);
      if(exist) readerMap[id] = ptr;
      return exist;
    }

    inline inptr get(
      const ::std::string &id /**< [in]读入工具名称 */
    ) const 
    {
      auto it = readerMap.find(id);
      if(it != readerMap.end()) return it->second;
      else return nullptr;
    }

  private:
    ::std::map<std::string, inptr> readerMap;
  };

  template<typename Base, typename T0>
  typename std::enable_if<::std::is_base_of<Base, T0>::value, bool>::type
    EnrollInput(
      InputMap<Base> &reader
    ) {
    auto ptr = std::make_shared<In<Base, T0>>();
    if (!(reader.emplace(ptr->name(), ptr))) {
      std::cout << "GeoInput : " << ptr->name() << " already exists in the reader" << std::endl;
      return false;
    }
    return true;
  }

  template<typename Base, typename T0>
  typename std::enable_if<::std::is_base_of<::phlbm::Input<Base>, T0>::value, bool>::type
    EnrollGeoInput(
      InputMap<Base> &reader
    ) {
    auto ptr = std::make_shared<T0>();
    if (!(reader.emplace(ptr->name(), ptr))) {
      std::cout << "GeoInput : " << ptr->name() << " already exists in the reader" << std::endl;
      return false;
    }
    return true;
  }

  template<
    typename Base,
    typename T0, 
    typename ... Ts,
    class = typename ::std::enable_if< sizeof...(Ts)>=1 >::type
  >
    void EnrollInput(
      InputMap<Base> &reader
    ) {
    EnrollInput<Base, T0>(reader);
    EnrollInput<Base, Ts ... >(reader);
  }

  template<typename Base>
  class Generator
  {
    InputMap<Base> inputMap;
  public:
    Generator(); 
    virtual ::std::shared_ptr<Base> operator() (
      const phlbm::NamedDataBlock &ndb, // [in,out]生成几何的数据
      const std::map<std::string, phlbm::NamedDataBlock> &ndbMap // [in,out]所有经风雷导入的数据 
    ) const final;
    const InputMap<Base> &GetInputMap() const { return inputMap; }
    InputMap<Base> &GetInputMap() { return inputMap; }
  };

  //typedef void (*EnrollReader)(InputMap<>& p);

  //extern COMMONAPI InputMap geoGenerator;

  //void COMMONAPI prepare_phlbm_geo_generator();
}

#endif

