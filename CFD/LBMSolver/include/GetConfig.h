///////////////////////////////////////////////////////////////
// author : He, Jiaming. email: 2457837217@qq.com
///////////////////////////////////////////////////////////////
// participant 1 : Wang,Junquan. email: as2697920806@gmail.com
///////////////////////////////////////////////////////////////
//GetConfig.h
#ifndef GETCONFIG_H
#define GETCONFIG_H
#include "CommonAPI.h"
#include<iostream>
#include<fstream>
#include<map>
#include<string>
/**nnwlbm命名空间*/
namespace nnwlbm {
  /**读入参数类,在map容器中存储从文件中读入的参数*/
  class COMMONAPI GetConfig
  {
  public:
    /**构造函数*/
    GetConfig() {}
    /**析构函数*/
    virtual ~GetConfig() {}
  public:
    /** 判断字符是否为空格
    * \param c [in] 需要判断的字符
    * \return  返回bool变量，true表示该字符是空格符；
    * false表示该字符不是空格符
    */
    static bool is_space
    (
      const char &c  /**< [in] 需要判断的字符 */
    );
    /** 判断字符是否为注释符"#"
    * \param c [in] 需要判断的字符
    * \return  返回bool变量，true表示该字符是注释符；
    * false表示该字符不是注释符
    */
    bool is_commentChar
    (
      const char &c  /**< [in] 需要判断的字符 */
    )const;
    /** 把字符串左边和右边的空格去掉
    * \param str [in,out] 需要处理的字符串
    */
    void trim
    (
      std::string &str  /**< [in,out] 需要处理的字符串 */
    )const;
    /** 分析一行字符串，如果是注释行，则不处理，如果是k-v行，则提取出key-value
    * \param line [in] 需要处理的字符串
    * \param key [in] key值
    * \param value [in] value值
    * \return  返回bool变量，true表示提取成功；
    * false表示无法提取
    */
    bool AnalyseLine
    (
      const std::string &line,  /**< [in] 需要处理的字符串*/
      std::string &key,  /**< [in] key值*/
      std::string &value  /**< [in] value值*/
    )const;
    /** 读取config文件并保存到map中
    * \param filename [in] 读入的config文件名称
    * \return  返回bool变量，true表示config文件读取成功；
    * false表示config读取失败，无法将内容存储到容器中
    */
    bool ReadConfig
    (
      const std::string &filename /**< [in] 读入的config文件名称 */
    );
    /** 将map中的数据输出到屏幕
    */
    void PrintConfig
    (
      void
    )const;
    /** 在容器m中找到key值对应的value值，并将其打印到屏幕
    * \param key [in] 关联容器中的key值
    * \return  返回对应的value值
    */
    double FindInConfig
    (
      const std::string &key  /**< [in] 关联容器中的key值 */
    );
  private:
    /**map容器，存储config文件内容*/
    std::map<std::string, double>ConfigMap;
  public:
    /**取出存储config文件内容的map容器*/
    std::map<std::string, double> &GetMap() { return ConfigMap; }
  };
}
#endif
// end of GetConfig.h
