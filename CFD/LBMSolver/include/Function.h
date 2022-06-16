//////////////////////////////////////////////////////////
//Function.h
//////////////////////////////////////////////////////////
//author:Guan, Qingshan. email:qingsg@hust.edu.cn
//////////////////////////////////////////////////////////
//participant 1:He, Shuangwu. email:shuangwuh@gmail.com.
//participant 2:Wang, Junquan. email:1935405719@qq.com.
//////////////////////////////////////////////////////////
//tester 1:He, Shuangwu. email:shuangwuh@gmail.com.
//tester 2:Wang, Junquan. email:1935405719@qq.com.

#ifndef FUNCTION_H
#define FUNCTION_H
#include "LBMGeometry.h"
namespace nnwlbm {
  /**Function基类，继承自DataClass
  * 构造一个函数，具有定义域和定义法则；
  * 私有成员包括维度和定义域。
  */
  class  COMMONAPI Function : public DataClass
  {
  public:
    /** Function的构造函数
    * 通过输入维度和定义域，构造出一个函数；
    * \param oDim 输出的维度，size_t类型；
    * \param geo 定义域，Geometry类型；
    */
    explicit Function(
      const size_t &oDim, 
      std::shared_ptr<Geometry>geo = std::shared_ptr<EmptyGeo>()
    ) : geometry(geo), outputDim(oDim) { }

    Function(const Function &fun) : geometry(fun.geometry), outputDim(fun.outputDim) { }

    virtual ~Function() {}

    /** 设置定义域
    * \param geoIn [in] 几何体内部为函数的定义域
    */
    void SetDomainofFunction(std::shared_ptr<const Geometry> geoIn) { geometry = geoIn; }
    std::shared_ptr<const Geometry> GetDomainofFunction() const { return geometry; }

    /**
    * 定义域检测函数；
    * 通过输入所计算的点坐标数组，计算是否在定义域内，如果不在，则输出”不在定义域内“；
    * \param *p [in] 所计算的点；
    * \return 返回bool变量表示输出状态，正常为true，异常为false。
    */
    inline bool is_applied(const double *p)const {
      if (geometry->CalRelativeLocation(p, 1e-9) >= geo::boundary) return true;
      std::cout << "不在定义域内" << std::endl;
      return false;
    }

    /**
    *返回类型名称. 通过name函数内部定义字符串实现一种类型一个名称。
    */
    virtual ::std::string name() const ;

    /**
    * 通过调用CalVarPrivate函数获得因变量；
    * \param *var [in,out]因变量所要改变的值；
    * \param *pos [in]自变量；
    * \return 返回bool变量表示输出状态，正常为true，异常为false。
    */
    inline bool cal_var(double *var, const double *pos) const
    {
      if (is_applied(pos)) return CalVarPrivate(var, pos);
      return false;
    }

    /**
    * 输出维度内联函数，维度从Geometry中获得；
    * 维度从Geometry中的GetDim函数中获得。
    */
    inline size_t get_input_dim() const { return (geometry)->GetDim(); }

    /**
    * 输出维度内联函数；
    * 将从Geometry获得的维度输出；
    */
    inline const size_t &GetOutputDim() const { return outputDim; }

  private:
    /** \brief 虚函数输出当前对象的二进制文件数据。
     * \param ostr [in,out]二进制文件输出流
     * \param filename [in]输出的文件名
     * \param path [in]输出的路径，默认为当前目录
     * \return 返回bool变量表示输出状态，正常为true，异常为false
     */
    virtual bool OutputBinPrivate(
      std::ofstream &ostr,  /**< [in,out]二进制文件输出流 */
      const std::string &filename, /**< [in]输出的文件名 */
      const std::string &path /**< [in]输出的路径，默认为当前目录 */
    )const; 

  protected:
    /**
    * 主函数，输入自变量获得因变量；
    * \param *var [in,out] 需要给因变量赋的值；
    * \param *pos [in]自变量；
    * \return 返回bool变量表示输出状态，正常为true，异常为false。
    */
    virtual bool CalVarPrivate(double *var, const double *pos) const = 0;

    /**
    * outputDim 维度；
    * Geometry 定义域。
    */
    size_t outputDim;
    std::shared_ptr<const Geometry> geometry;
  };

  /** 常函数
  */
  class COMMONAPI ConstFun : public Function
  {
  public:
    /**
    *返回二进制输出名称. 通过name函数内部的常字符串定义名称。
    */
    virtual ::std::string name()const;

    /**
    * 常量值。
    */
    std::vector<double> var;

    /**
    * ConstFun构造函数
    * 输入一个double值，一个Goemetry值,和一个vector实现构造一个ConstFun；
    * \param dim 输出维度；
    * \param geo 定义域：
    * \param forVar 固定的因变量值。
    */
    explicit ConstFun(const std::shared_ptr<Geometry> &geo, const std::vector<double> &forVar) :
      Function(forVar.size(), geo), var(forVar) {}
  private:
    /** \brief 虚函数输出当前对象的二进制文件数据。
     * \param ostr [in,out]二进制文件输出流
     * \param filename [in]输出的文件名
     * \param path [in]输出的路径，默认为当前目录
     * \return 返回bool变量表示输出状态，正常为true，异常为false
     */
    virtual bool OutputBinPrivate(
      std::ofstream &ostr,  /**< [in,out]二进制文件输出流 */
      const std::string &filename, /**< [in]输出的文件名 */
      const std::string &path /**< [in]输出的路径，默认为当前目录 */
    )const;
  protected:
    /**
    * 主函数，输入两个double数组；
    * \param *va [in,out]因变量坐标数组；
    * \param *pos [in]自变量坐标数组；
    * 通过因变量得到自变量的计算过程；
    * \return 返回bool变量表示输出状态，正常为true，异常为false。
    */
    virtual bool CalVarPrivate(double *va, const double *pos) const;
  };

  /**线性函数，继承自Function
  * 私有数组包括斜率和截距。
  */
  class COMMONAPI LinearFun : public Function
  {
  public:
    /** linearFun的构造函数
    * 要求输入的几何维度等于slopeIn和offsetIn的大小
    * \param geo 几何体，表示函数的定义域
    * \param slopeIn 输入斜率或梯度
    * \param offsetIn 输入截距
    */
    explicit LinearFun(
      const std::shared_ptr<Geometry> &geo, 
      const std::vector<double> &slopeIn, 
      const double &offsetIn 
    );
  private:
    std::vector<double> slope; /**< 斜率 */
    const double offset; /**< 截距 */
  public:
    /**
    *返回类型名称. 通过name函数内部定义字符串实现一种类型一个名称。
    */
    virtual ::std::string name()const;

    /**
    * 返回斜率
    */
    inline const std::vector<double> get_slope() const { return slope; }

    /**
    * 返回截距
    */
    inline const double get_offset()const { return offset; }

  private:
    /** \brief 虚函数输出当前对象的二进制文件数据。
     * \param ostr [in,out]二进制文件输出流
     * \param filename [in]输出的文件名
     * \param path [in]输出的路径，默认为当前目录
     * \return 返回bool变量表示输出状态，正常为true，异常为false
     */
    virtual bool OutputBinPrivate(
      std::ofstream &ostr,  /**< [in,out]二进制文件输出流 */
      const std::string &filename, /**< [in]输出的文件名 */
      const std::string &path /**< [in]输出的路径，默认为当前目录 */
    )const;

  protected:
    /**
    * 主函数，输入两个double数组；
    * \param *va [in,out]因变量数组；
    * \param *pos [in]自变量数组；
    * 通过自变量得到因变量的计算过程；
    * \return 返回bool变量表示输出状态，正常为true，异常为false。
    */
    virtual bool CalVarPrivate(double *va, const double *pos) const ;
  };

  class COMMONAPI AssemblyFun : virtual public Function
  {
  public:
    explicit AssemblyFun(const std::vector<std::shared_ptr<Function>> &funVecIn);
  protected:
    virtual bool CalVarPrivate(
      double *var, 
      const double *pos
    ) const;

    size_t CalOutputDim() const;
  private:
    std::vector<std::shared_ptr<Function>> funVec;
  };
}


#endif
