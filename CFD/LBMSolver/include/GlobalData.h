///////////////////////////////////////////////////////////////
// author : He, Jiaming. email: 2457837217@qq.com
///////////////////////////////////////////////////////////////
// participant 1 : Chen, Songze. email: jacksongze@hust.edu.cn
///////////////////////////////////////////////////////////////

//GlobalData.h
#ifndef GLOBALDATA_H
#define GLOBALDATA_H
#include "CommonAPI.h"
#include <memory>
#include "BinBase.h"

/**nnwlbm命名空间*/
namespace nnwlbm {

  /**存储数据类,在map容器中存储指向某类型数据类对象的智能指针*/
  class COMMONAPI GlobalData : public DataClass
  {
  public:
    /**构造函数*/
    GlobalData() {}
    /**析构函数*/
    virtual ~GlobalData() {}
    /**复制函数*/
    GlobalData(
      const GlobalData& inObj
    ) : dataMap(inObj.dataMap) {}

  public:
    /** \brief 根据要取出的指针所对应的键取出map容器中对应的指针。
    * \param key [in] 要取出的指针所对应的键
    * \return 成功取出则返回输入的键所对应的指针，失败则返回空指针
    */
    inline Dptr GetMember(
      const ::std::string &key  /**< [in] 要取出的指针所对应的键*/
    ) const
    {
      auto it = dataMap.find(key);
      if (it != dataMap.end()) return it->second;
      else return nullptr;
    }

    /** \brief 将指针及其对应的键存入map容器中。
    * \param key [in] 要存入的指针对应的键
    * \param ptr [in] 要存入map的指针
    * \return 成功存入则返回true,未成功存入则返回false
    */
    inline bool emplace(
      const ::std::string &key,    /**< [in] 要存入的指针对应的键*/
      Dptr ptr                     /**< [in] 要存入map的指针*/
    )
    {
      return dataMap.emplace(key, ptr).second;
    }

    /** \brief 将map容器中指定的键所对应的指针用新指针进行替换。
    * 若没有与指定的键对应的指针，则将新指针与之对应
    * 若有与指定的键对应的指针，则将该指针替换为新指针
    * \param key [in] 指定的键
    * \param ptr [in] 新指针
    * \return 无论上述哪种情况都返回true
    */
    inline bool replace(
      const ::std::string &key,    /**< [in] 指定的键*/
      Dptr ptr                     /**< [in] 新指针*/
    )
    {
      bool isBool = dataMap.emplace(key, ptr).second;
      if (isBool)
      {
        return true;
      }
      else
      {
        dataMap[key] = ptr;
        return true;
      }
    }

    /**取出存储数据类对象的map容器*/
    inline const DataMap &GetMap() const { return dataMap; }

    /** \brief 模板函数。
    * 将输入的基类指针转换为T类型的指针
    * \param ptr [in,out] 输入的基类指针
    * \return 成功则返回T类型的指针，失败则返回空指针
    */
    template<class T>
    static T *GetPtr(
      Dptr ptr       /**< [in,out] 输入的基类指针*/
    )
    {
      if (ptr == nullptr) return nullptr;
      else return dynamic_cast<T*>(ptr.get());
    }
    /** \brief 模板函数。
    * 将输入的键在map容器中对应的基类指针转换成T类型指针
    * \param key [in] 进行类型转换的指针所对应的键
    * \return 成功则返回T类型的指针，失败则返回空指针
    */
    template<class T>
    T *GetPtr(
      const ::std::string &key      /**< [in] 进行类型转换的指针所对应的键*/
    )
    {
      auto ptr = GetMember(key);
      return GetPtr<T>(ptr);
    }
    /** \brief 模板函数。
    * 根据输入的键，将map容器中其对应的指针转换成T类型后返回其所指对象的引用
    * \param key [in] 要转换类型并输出所指对象引用的指针所对应的键
    * \return 成功则返回T类型指针所指对象的引用，失败则系统报错
    */
    template<class T>
    T &get_ref(
      const ::std::string &key     /**< [in] 要转换类型并输出所指对象引用的指针所对应的键*/
    ) const
    {
      T *ans = nullptr;
      auto ptr = GetMember(key);
      if (ptr == nullptr) {
        ::std::cout << "Cannot find " << key << ::std::endl;
#ifdef _WIN32
        system("pause");
#endif // _WIN32
        exit(1);
      }
      ans = dynamic_cast<T*>(ptr.get());
      if (ans == nullptr)
      {
        ::std::cout << "Cannot convert " << key << " to " << "target " << "class " << ::std::endl;
#ifdef _WIN32
        system("pause");
#endif // _WIN32
        exit(1);
      }
      else
        return *ans;
    }
  private:
    /** \brief 虚函数读出当前对象的二进制文件数据。
    * 读出方式取决与需要输出的数据的存储方式
    * \param ostr [in,out] 二进制文件输出流
    * \param filename [in] 输出的文件名
    * \param path [in] 输出的路径，默认为当前目录
    * \return 输出状态正常返回true，异常返回false
    */
    virtual bool OutputBinPrivate(
      std::ofstream &ostr,          /**< [in,out] 二进制文件输出流 */
      const std::string &filename,  /**< [in] 输出的文件名 */
      const std::string &path       /**< [in] 输出的路径，默认为当前目录 */
    )const;

    /**虚函数返回类型名，实现一种类型一个名称*/
    virtual ::std::string name() const;
  private:
    /**map容器，存储指向数据的指针*/
    DataMap dataMap;
  };
}
#endif
  // end of GlobalData.h
