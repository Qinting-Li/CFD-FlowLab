///////////////////////////////////////////////////////////////
// author : Chen, Songze. email: jacksongze@hust.edu.cn
///////////////////////////////////////////////////////////////
// participant 1 : He, Jiaming. email: 
///////////////////////////////////////////////////////////////
// tester 1 : 
///////////////////////////////////////////////////////////////

// BinBase.h
#include "DevEnv.h"

#ifndef BINBASE_H
#define BINBASE_H

#include "CommonAPI.h"
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <memory>
#include <iomanip>

extern "C" int COMMONAPI create_directory(const std::string &directoryPath);

namespace nnwlbm {
  /** \brief nnwlbm中所有数据类型的基类，用于维护二进制输出接口。
  * DataClass的派生类A需要实现二进制输出的功能。具体如下：
  * 1.重载name()函数，定义派生类的数据名称
  * 2.重载 函数，
  *   OutputBinPrivate( 
  *     std::ofstream& ostr,  
  *     const std::string& filename, 
  *     const std::string& path 
  *   ) const 
  * 3.在派生类A外部，建立一个BinInput的派生类B用于读入A的二进制数据
  */
  class COMMONAPI DataClass {
    public:
      virtual ~DataClass(){}
      /** \brief 静态成员函数实现二进制输出功能。
      * 在指定路径下将数据ref输出为二进制文件
      * \param ref [in]需要输出的数据，类型为DataClass的对象，
      * \param filename [in]输出的文件名
      * \param path [in]输出的路径，默认为当前目录
      * \return 返回bool变量表示输出状态，正常为true，异常为false
      */
      static bool output_bin(
        const DataClass &ref, /**< [in]需要输出的数据 */
        const std::string &filename,  /**< [in]输出的文件名 */
        const std::string &path = ".\\" /**< [in]输出的路径，默认为当前目录 */
      );

    private:
      /** \brief 输出二进制文件头。
      * 文件头的前4字节储存一个整型数据，表示此后的字符串长度，该字符串表示类型的名称（见name()函数）
      * \param ostr [in,out]二进制文件输出流 
      * \param s [in] 文件头的内容，一般为输出类型的名称 
      * \return 返回bool变量表示输出状态，正常为true，异常为false
      */
      static bool OutputHead(
        std::ofstream &ostr, /**< [in,out]二进制文件输出流 */
        const std::string &s  /**< [in] 文件头的内容，一般为输出类型的名称 */
      );

      /** \brief 虚函数输出当前对象的二进制文件数据。
      * DataClass的派生类实现自己的二进制输出
      * \param ostr [in,out]二进制文件输出流 
      * \param filename [in]输出的文件名
      * \param path [in]输出的路径，默认为当前目录
      * \return 返回bool变量表示输出状态，正常为true，异常为false
      */

      virtual bool OutputBinPrivate(
        std::ofstream &ostr,  /**< [in,out]二进制文件输出流 */
        const std::string &filename, /**< [in]输出的文件名 */
        const std::string &path /**< [in]输出的路径，默认为当前目录 */
      ) const ;

      /** 返回类型名称. 通过name函数内部定义字符串实现一种类型一个名称 */
      virtual ::std::string name() const = 0;
  };

  /** nnwlbm某种数据类型的二进制文件的读入工具的基类。
   */
  class COMMONAPI BinInput {
  public:
    /** 返回类型名称. 通过name函数内部定义字符串实现一种类型一个名称，
    * 该名称须于读入的Data种类的name一致
    */
    virtual ::std::string name() const = 0;
  private:
    /** 虚函数实现某种数据类型的二进制数据读入功能 
    * \param istr [in,out]二进制文件输入流 
    * \param filename [in]读入的文件名
    * \param path [in]文件路径，默认路径为当前目录
    */
    virtual std::shared_ptr<DataClass> InputBinPrivate(
      std::ifstream &istr, /**< [in,out]二进制文件输入流 */
      const std::string &filename, /**< [in]读入的文件名 */
      const std::string &path=".\\" /**< [in]文件路径，默认路径为当前目录 */
    ) const = 0;

  public:
    friend class BinReader;
  };

  /** 该模板函数用于简化二进制文件读入类的书写
  * 具体用法：
  * 1.已有DataClass的派生类A，可通过该模板类快速建立BinInput的派生类BinIn<A>，实现读入A的二进制文件
  * 2.写出函数体
  * std::shared_ptr<DataClass> BinIn<A>::InputBinPrivate(
  *   std::ifstream& istr,
  *   const std::string& filename,
  *   const std::string& path
  * ) const { // 加入读入文件的代码  }
  * 3.系统通过BinReader类调用InputBinPrivate()函数
  */
  template<typename T>
  class BinIn : public BinInput{
  public:
    /** 返回类型名称. 通过name函数内部定义字符串实现一种类型一个名称，
    * 该名称须于读入的Data种类的name一致
    */
    virtual ::std::string name() const ;

  private:
    /** 虚函数实现某种数据类型的二进制数据读入功能
    * \param istr [in,out]二进制文件输入流
    * \param filename [in]读入的文件名
    * \param path [in]文件路径，默认路径为当前目录
    */
    virtual std::shared_ptr<DataClass> InputBinPrivate(
      std::ifstream &istr, /**< [in,out]二进制文件输入流 */
      const std::string &filename, /**< [in]读入的文件名 */
      const std::string &path=".\\" /**< [in]文件路径，默认路径为当前目录 */
    ) const;
  };

  /** 读入nnwlbm的二进制文件的入口类。负责维护所有数据类型的二进制读入工具。
   */
  class COMMONAPI BinReader {
  public:
    /** 智能指针类型别名 */
    typedef std::shared_ptr<BinInput> inptr;

    /** 读入DataClass为基类的二进制文件。
    * \param filename [in]读入的文件名
    * \param path [in]文件路径，默认路径为当前目录
    * \return 返回由二进制文件数据产生的DataClass对象的智能指针
    */
    virtual ::std::shared_ptr<DataClass> read(
      const std::string &filename, /**< [in]文件名 */
      const std::string &path = ".\\" /**< [in]文件路径 */
    ) const final;

    /** 插入名字为id的某种数据类型的二进制读入工具的指针ptr 
    * \param id [in]读入工具名称
    * \param ptr [in]某种数据的读入工具指针
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

    /** 替换名字为id的二进制读入工具的指针ptr 
    * \param id [in]读入工具名称
    * \param ptr [in]某种数据的读入工具指针
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

    /** 取出名字为id的二进制读入工具 
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

  protected:
    /** 读入DataClass的二进制文件头 
    * \param istr [in,out] 二进制文件输入流
    * \return 返回字符串，表示二进制文件头
    */
    inline ::std::string InputHead(
      std::ifstream &istr /**< [in,out] 二进制文件输入流 */
    ) const
    {
      //istr.setf( std::ios::binary );
      char tempMem[256];
      istr.read(tempMem, 4);
      unsigned nameLength = *((unsigned*)tempMem);
      if(nameLength > 256) return std::string();
      istr.read(tempMem, nameLength);
      return std::string(tempMem, nameLength);
    }

  private:
    /** 用于存放读入工具，建立工具与名字之间的映射 */
    ::std::map<std::string, inptr> readerMap;
  };

  /** 注册各种二进制读入类的模板函数的出口函数
  */
  template<typename T0>
  typename std::enable_if<::std::is_base_of<::nnwlbm::DataClass, T0>::value, bool>::type
    enroll_binInput(
      BinReader &reader
    ) {
    auto ptr = std::make_shared<BinIn<T0>>();
    if (!(reader.emplace(ptr->name(), ptr))) {
      std::cout << "BinInput : " << ptr->name() << " already exists in the reader" << std::endl;
      return false;
    }
    return true;
  }

  template<typename T0>
  typename std::enable_if<::std::is_base_of<::nnwlbm::BinInput, T0>::value, bool>::type
    enroll_binInput(
      BinReader &reader
    ) {
    auto ptr = std::make_shared<T0>();
    if (!(reader.emplace(ptr->name(), ptr))) {
      std::cout << "BinInput : " << ptr->name() << " already exists in the reader" << std::endl;
      return false;
    }
    return true;
  }

  template<
    typename T0, 
    typename ... Ts,
    class = typename ::std::enable_if< sizeof...(Ts)>=1 >::type
  >
  void enroll_binInput(
    BinReader &reader
  ) {
    enroll_binInput<T0>(reader);
    enroll_binInput< Ts ... >(reader);
  }
  /** 注册函数的公共接口 
  * \param p 测试类列表的指针，用于存储各种测试类的对象，由外部提供。
  */
  typedef void (*EnrollBinReader)(BinReader &p);

  extern COMMONAPI BinReader binReader;
  typedef ::std::shared_ptr<DataClass> Dptr;
  typedef ::std::map<::std::string, Dptr> DataMap;

  /** 注册nnwlbm标准二进制读入工具 
  * 在使用nnwlbm::binReader之前，需调用本函数。
  */
  void COMMONAPI prepare_nnwlbm_bin_reader();

  /** 注册用户自定义的二进制读入工具 
  * \param r 类型为nnwlbm::BinReader，用于存放二进制读入工具。
  * \param dllName 类型为字符串，表示二进制读入工具所在动态链接库的名字。
  * \return 返回bool型，若动态链接库读入失败，返回false，否则返回true。
  */
  bool COMMONAPI enroll_bin_reader_from(BinReader &r, const std::string &dllName);
}

#endif
// end of BinBase.h
