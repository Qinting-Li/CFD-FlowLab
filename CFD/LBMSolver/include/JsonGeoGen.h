/////////////////////////////////////////////////////////////////////
// to add header.
/////////////////////////////////////////////////////////////////////
// author: Chen, Songze. email: jacksongze@hust.edu.cn
/////////////////////////////////////////////////////////////////////
// participant 1: 
/////////////////////////////////////////////////////////////////////

#pragma once

#ifndef PHGEOGEN_H
#define PHGEOGEN_H
#include <CommonAPI.h>
#include <LBMGeometry.h>
#include "JsonReader.h"
#include <string>
#include <ostream>
#include <memory>

namespace nnwlbm 
{
  /** 利用风雷的文本处理工具读取几何类型数据，生成nnwlbm::Geometry的派生对象
   */
  class COMMONAPI GeoInput {
  public:
    /** 返回类型名称. 通过name函数内部定义字符串实现一种类型一个名称，
    * 该名称须于读入的几何类型的name一致
    */
    virtual ::std::string name() const = 0;
  private:
    /** 虚函数实现某种数据类型的二进制数据读入功能 
    * \param ndb [in,out]经风雷导入的数据 
    */
    virtual std::shared_ptr<nnwlbm::Geometry> InputGeoPrivate(
      const phlbm::NamedDataBlock &ndb, // [in,out]生成几何的数据
      const std::vector<phlbm::NamedDataBlock> &ndbVec // [in,out]所有经风雷导入的数据 
    ) const = 0;

  public:
    friend class GeoGenerator;
  };

  /** 该模板函数用于简化几何类生成工具的编写
  * 具体用法：
  * 1.已有nnwlbm::Geometry的派生类A，可通过该模板类快速建立GeoInput的派生类GeoIn<A>
  * 2.写出函数体
  * ~~~~
  * std::shared_ptr<nnwlbm::Geometry> GeoIn<A>::InputGeoPrivate(
  *   phlbm::NamedDataBlock& ndb // [in,out]经风雷导入的数据 
  * ) const { // 加入读入文件的代码  }
  * ~~~~
  * 3.系统通过GeoGenerator类调用InputGeoPrivate()函数
  */
  template<typename T>
  class GeoIn : public GeoInput{
  public:
    /** 返回类型名称. 通过name函数内部定义字符串实现一种类型一个名称，
    * 该名称须于读入的几何类的name一致
    */
    virtual ::std::string name() const;

  private:
    /** 虚函数实现某种数据类型的几何生成功能 
    * \param ndb [in,out]经风雷导入的数据 
    */
    virtual std::shared_ptr<nnwlbm::Geometry> InputGeoPrivate(
      const phlbm::NamedDataBlock& ndb, // [in,out]生成几何的数据
      const std::vector<phlbm::NamedDataBlock> &ndbVec // [in,out]所有经风雷导入的数据 
    ) const;
  };

  /** 生成几何的入口类。负责维护所有几何生成工具。
   */
  class COMMONAPI GeoGenerator {
  public:
    int identity() { return 1; }
    /** 智能指针类型别名 */
    typedef std::shared_ptr<GeoInput> inptr;

    /** 生成几何类的入口程序
    * \param ndb [in]读入的文件名
    * \param path [in]文件路径，默认路径为当前目录
    * \return 返回由几何生成数据产生的DataClass对象的智能指针
    */
    virtual ::std::shared_ptr<nnwlbm::Geometry> gen(
      const phlbm::NamedDataBlock &ndb, // [in,out]生成几何的数据
      const std::vector<phlbm::NamedDataBlock> &ndbVec // [in,out]所有经风雷导入的数据 
    ) const final;

    /** 插入名字为id的某种数据类型的几何生成工具的指针ptr 
    * \param id [in]读入工具名称
    * \param ptr [in]某种几何类的生成工具指针
    * \return 返回bool变量，true表示不存在名为id的对象，成功插入；
    * false表示已存在名为id的对象，没有正确插入
    */
    inline bool emplace(
      const ::std::string &id, /**< [in]读入工具名称 */
      inptr ptr /**< [in]某种数据的读入工具指针 */
    ) 
    {
      return readerMap.emplace(id, ptr).second;
    }

    /** 替换名字为id的几何生成工具的指针ptr 
    * \param id [in]读入工具名称
    * \param ptr [in]某种几何类的生成工具指针
    * \return 返回bool变量，true表示已存在名为id的对象，替换成功；
    * false表示不存在名为id的对象，直接插入。
    */
    inline bool replace(
      const ::std::string &id, /**< [in]读入工具名称 */
      inptr ptr /**< [in]某种数据的读入工具指针 */
    ) 
    { 
      bool exist = !(readerMap.emplace(id, ptr).second);
      if(exist) readerMap[id] = ptr;
      return exist;
    }

    /** 取出名字为id的几何生成工具 
    * \param id [in]读入工具名称
    * \return 返回读入工具的智能指针
    */
    inline inptr get(
      const ::std::string &id /**< [in]读入工具名称 */
    ) const 
    {
      auto it = readerMap.find(id);
      if(it != readerMap.end()) return it->second;
      else return nullptr;
    }

  private:
    /** 用于存放几何生成工具，建立工具与名字之间的映射 */
    ::std::map<std::string, inptr> readerMap;
  };

  /** 注册各种几何生成的模板函数的入口函数
  */
  template<typename T0>
  typename std::enable_if<::std::is_base_of<::nnwlbm::Geometry, T0>::value, bool>::type
    EnrollGeoInput(
      GeoGenerator& reader
    ) {
    auto ptr = std::make_shared<GeoIn<T0>>();
    if (!(reader.emplace(ptr->name(), ptr))) {
      std::cout << "GeoInput : " << ptr->name() << " already exists in the reader" << std::endl;
      return false;
    }
    return true;
  }

  template<typename T0>
  typename std::enable_if<::std::is_base_of<::phlbm::GeoInput, T0>::value, bool>::type
    EnrollGeoInput(
      GeoGenerator& reader
    ) {
    auto ptr = std::make_shared<T0>();
    if (!(reader.emplace(ptr->name(), ptr))) {
      std::cout << "GeoInput : " << ptr->name() << " already exists in the reader" << std::endl;
      return false;
    }
    return true;
  }

  template<
    typename T0, 
    typename ... Ts,
    class = typename ::std::enable_if< sizeof...(Ts)>=1 >::type
  >
  void EnrollGeoInput(
    GeoGenerator &reader
  ) {
    EnrollGeoInput<T0>(reader);
    EnrollGeoInput< Ts ... >(reader);
  }
  /** 注册函数的公共接口 
  * \param p 测试类列表的指针，用于存储各种测试类的对象，由外部提供。
  */
  typedef void (*EnrollGeoReader)(GeoGenerator& p);

  extern COMMONAPI GeoGenerator geoGenerator;

  /** 注册phlbm标准二进制读入工具 
  * 在使用phlbm::GenGenerator之前，需调用本函数。
  */
  void COMMONAPI prepare_phlbm_geo_generator();
};

#endif

