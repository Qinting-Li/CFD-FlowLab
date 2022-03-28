///////////////////////////////////////////////////////////////
// author : Wang,Junquan. email: as2697920806@gmail.com
///////////////////////////////////////////////////////////////
// participant 1 : Chen,Songze. email: 
///////////////////////////////////////////////////////////////
// tester 1 : 
///////////////////////////////////////////////////////////////

#ifndef BOUNDARYCONDITION_H
#define BOUNDARYCONDITION_H

#include "CommonAPI.h"
#include "Lattice.h"
#include "Collision.h"

namespace nnwlbm
{

  /** \brief
    * 单点边界条件(SingleNodeBc)类：作为单点边界条件类的父类 \n
    * operator(): 调用以执行具体的边界条件格式 \n
    * 示例，生成一个名为bounceCondition的对象：
    ~~~~~{.cpp}
    SingleNodeBc<LC> bounceCondition;
    ~~~~~
    * 示例，调用operator()执行具体的边界条件格式：
    ~~~~~{.cpp}
    bounceCondition();
    ~~~~~
  */
  template<typename LC>
  class SingleNodeBc 
  {
  public:
    explicit SingleNodeBc() {}
    virtual ~SingleNodeBc() {}
    virtual bool operator ()(
      double *pdfOut,  // [out]分布函数数组
      const double *pdfIn,  // [in]碰撞后的分布函数
      const size_t &i,  // [in]离散速度方向
      const double *macroVar // [in] 与边界相关的宏观量 
    ) const { return false; }
  };

  /** \brief
    * BounceBack类：用于执行简单的标准碰撞边界条件 \n
    * operator()函数: 调用以执行简单的标准碰撞边界条件 \n
    * 示例，生成一个名为bounceBack的对象：
    ~~~~~{.cpp}
    BounceBack<nnwlbm::D2Q9> bounceBack;
    ~~~~~
    * 示例，调用operator()执行标准碰撞边界条件：
    ~~~~~{.cpp}
    bounceBack(pdfOut,pdfIn,i,macroVar);
    ~~~~~
  */
  template<typename LC>
  class BounceBack : virtual public SingleNodeBc<LC>
  {
  public:
    /** \brief BounceBack 类的构造函数:
     */
    explicit BounceBack() {}
    ~BounceBack() {}


    /** \brief 实现离散速度方向上的反弹格式
     *\param pdfOut [out] 分布函数数组
     *\param pdfIn [in] 碰撞后的分布函数
     *\param i [in]离散速度方向
     *\param macroVar 格点上的宏观量
     *\return 返回bool变量true
    */
    virtual bool operator ()(
      double *pdfOut,  // [out]分布函数数组
      const double *pdfIn,  // [in]碰撞后的分布函数
      const size_t &i,  // [in]离散速度方向
      const double *macroVar // [in] 宏观量 
    ) const 
    {
      pdfOut[i] = pdfIn[LC::InvertedDirection(i)];
      return true;
    }
  };

  /** \brief
    * MovingBounceBack类：用于执行移动边界的反弹边界条件 \n
    * 数据包括：<格子速度、移动边界上的宏观速度> \n
    * operator(): 调用以执行移动边界的反弹边界条件 \n
    * 示例，生成一个名为movingBounceBack的对象：
    ~~~~~{.cpp}
    double wallMacro[2]={1.0,0.0};
    MovingBounceBack<nnwlbm::D2Q9> movingBounceBack(LC::velocity(),wallMacro);
    ~~~~~
    * 示例，调用operator()执行边界条件：
    ~~~~~{.cpp}
    movingBounceBack(pdfOut,pdfIn,i,macroVar);
    ~~~~~
  */
  template<typename LC>
  class MovingBounceBack : virtual public SingleNodeBc<LC>
  {
  public:
    /** \brief MovingBounceBack 类的构造函数:
     *\param latticeVel [in] 格子速度
     *\param wVel [in] 移动边界的宏观速度
    */
    explicit MovingBounceBack(const double &latticeVel, const double *wVel) 
      : latticeVelocity(latticeVel) 
    { 
      for(size_t i = 0; i < LC::GetDim(); i++ ) wallVel[i] = wVel[i];
    }

    ~MovingBounceBack() {}

    /** \brief 实现离散速度方向上的移动边界反弹格式
     *\param pdfOut [out] 分布函数数组
     *\param pdfIn [in] 碰撞后的分布函数
     *\param i [in]离散速度方向
     *\param macroVar 格点处的宏观量
     *\return 返回bool变量true
     */
    virtual bool operator ()(
      double *pdfOut,  // [out]分布函数数组
      const double *pdfIn,  // [in]碰撞后的分布函数
      const size_t &i,  // [in]离散速度方向
      const double *macroVar // [in] 宏观量
    ) const 
    {
      int invI = LC::InvertedDirection(i);
      double ev = 0.0;
      for(size_t iDim = 0; iDim < LC::GetDim(); iDim++) {
        ev += LC::direction(invI)[iDim] * wallVel[iDim];
      }
      pdfOut[i] = pdfIn[LC::InvertedDirection(i)] 
        - 6.0 * LC::weight(i) * macroVar[0]  * ev/latticeVelocity;

      return true;
    }
  protected:
    double latticeVelocity; //格子速度
    double wallVel[LC::GetDim()]; //移动边界的宏观速度
  };

//==================================================================
  class SimpleBC
  {
  public:
    explicit SimpleBC() {}
    virtual ~SimpleBC() {}
    virtual bool operator ()(
      double *pdfOut,  // [out]分布函数数组
      const double *pdfIn,  // [in]碰撞后的分布函数
      const size_t &i
      ) const {
      return false;
    }
  };

  /** \brief
    * BounceBackBc类：用于执行简单的标准碰撞边界条件 \n
    * 数据包括：<离散速度方向编号> \n
    * operator(): 调用以执行简单的标准碰撞边界条件 \n
    * 示例，生成一个名为bounceBack的对象：
    ~~~~~{.cpp}
    nnwlbm::BounceBackBc<nnwlbm::D2Q9> bounceBack();
    ~~~~~
    * 示例，调用operator()执行标准碰撞边界条件：
    ~~~~~{.cpp}
    bounceBack(pdfOut,pdfIn,i);
    ~~~~~
  */
  template<typename LC>
  class BounceBackBc : virtual public SimpleBC
  {
  public:
    /** \brief BounceBackBc 类的构造函数:
     *\param iIn [in] 离散速度方向编号
    */
    explicit BounceBackBc() {}
    ~BounceBackBc() {}

    /** \brief 实现离散速度方向上的标准反弹格式
     *\param pdfOut [out] 分布函数数组
     *\param pdfIn [in] 碰撞后的分布函数
     *\param macroVar [in] 该处的宏观量
     *\return 返回bool型变量true
     */
    virtual bool operator ()(
      double *pdfOut,  // [out]分布函数数组
      const double *pdfIn,  // [in]碰撞后的分布函数
      const size_t &i // [in] 与边界相关的宏观量 
      ) const
    {
      pdfOut[i] = pdfIn[LC::InvertedDirection(i)];
      return true;
    }
  };

  /** \brief
    * MovingBounceBackBc类：用于执行移动边界对应的碰撞边界条件 \n
    * 数据包括：<格子速度、移动边界的宏观速度、离散速度方向编号> \n
    * operator(): 调用以执行移动边界碰撞边界条件 \n
    * 示例，生成一个名为movingBounceBack的对象：
    ~~~~~{.cpp}
    double macroWall[nnwlbm::D2Q9::GetDim()] = { 6.0,0.0 };
    nnwlbm::MovingBounceBackBc<nnwlbm::D2Q9> movingBounceBack(i, lc.velocity(), macroWall);
    ~~~~~
    * 示例，调用operator()执行标准碰撞边界条件：
    ~~~~~{.cpp}
    movingBounceBack(pdf, pdfPost, macro);
    ~~~~~
  */

  template<typename LC>
  class MovingBounceBackBc : virtual public SimpleBC
  {
  public:
    /** \brief MovingBounceBackBc 类的构造函数:
     *\param iIn [in] 离散速度方向编号
     *\param latticeVel [in] 格子速度
     *\param wVel [in] 移动边界对应的宏观速度
    */
    explicit MovingBounceBackBc(const double& macroRhoIn, const double& latticeVel, const double* wVel)
      : latticeVelocity(latticeVel), macroRho(macroRhoIn)
    {
      for (size_t i = 0; i < LC::GetDim(); i++) wallVel[i] = wVel[i];
    }

    ~MovingBounceBackBc() {}

    /** \brief 实现离散速度方向上的移动边界反弹格式
     *\param pdfOut [out] 分布函数数组
     *\param pdfIn [in] 碰撞后的分布函数
     *\param macroVar [in] 边界格点的宏观量
     *\return 返回bool型变量true
     */
    virtual bool operator ()(
      double *pdfOut,  // [out]分布函数数组
      const double *pdfIn,  // [in]碰撞后的分布函数
      const size_t &i  //[in]离散速度方向
      ) const
    {
      int invI = LC::InvertedDirection(i);
      double ev = 0.0;

      for (size_t iDim = 0; iDim < LC::GetDim(); iDim++) {
        ev += LC::direction(invI)[iDim] * wallVel[iDim];
      }
      pdfOut[i] = pdfIn[LC::InvertedDirection(i)] - 6.0 * LC::weight(i)*macroRho * ev / latticeVelocity;

      return true;
    }
  protected:
    double latticeVelocity;  //格子速度
    double wallVel[LC::GetDim()];  //移动边界上的宏观速度
    double macroRho;
  };

  /** \brief
    * NonEquilibriumExtrapolation类：用于执行非平衡外推边界条件 \n
    * 数据包括：<相邻格点的分布函数指针、相邻格点的宏观量指针、给定的边界宏观量、离散速度方向编号> \n
    * operator(): 调用以执行非平衡外推边界条件 \n
    * 示例，生成一个名为nonEquilibriumExtrapolation的对象：
    ~~~~~{.cpp}
    nnwlbm::NonEquilibriumExtrapolation<nnwlbm::D2Q9> nonEquilibriumExtrapolation(i, nearNodePdf, nearNodeMacro, nodeMacro);
    ~~~~~
    * 示例，调用operator()执行非平衡外推边界条件：
    ~~~~~{.cpp}
     nonEquilibriumExtrapolation(nodePdf, nodePdfPost, i);
    ~~~~~
  */
  template<typename LC>
  class NonEquilibriumExtrapolation : virtual public SimpleBC
  {
  public:
    /** \brief NonEquilibriumExtrapolation 类的构造函数:
     *\param iIn [in] 离散速度方向编号
     *\param nearNodePdfIn [in] 相邻格点的分布函数
     *\param nearNodeMacroIn [in] 相邻格点的宏观量
     *\param wVel [in] 边界给定的宏观速度
    */
    explicit NonEquilibriumExtrapolation(
      const double *nearNodePdfIn,
      const double *wVel,
      const double latticeVelIn) :latticeVel(latticeVelIn), nearNodePdf(nearNodePdfIn) {
      for (size_t i = 0; i < LC::GetDim(); i++) boundaryVel[i] = wVel[i];
    }


    ~NonEquilibriumExtrapolation() {}

    /** \brief 实现离散速度方向上的非平衡外推格式
     *\param pdfOut [out] 分布函数数组
     *\param pdfIn [in] 碰撞后的分布函数
     *\param macroVar [in] 宏观量
     *\return 返回bool型变量true
    */
    virtual bool operator ()(
      double *pdfOut,  // [out]分布函数数组
      const double *pdfIn,  // [in]碰撞后的分布函数
      const size_t &i
      ) const
    {
      nnwlbm::Equilibrium<LC> equilibrium;

      LC lc(latticeVel);
      double temp[LC::GetDim() + 1];
      nnwlbm::cal_macro_var<LC>(temp, nearNodePdf, lc);

      double boundaryMacroVar[LC::GetDim() + 1];
      boundaryMacroVar[0] = temp[0];
      for (size_t i = 0; i < LC::GetDim(); i++)
        boundaryMacroVar[i + 1] = boundaryVel[i];

      pdfOut[i] = equilibrium(boundaryMacroVar, i, latticeVel) + nearNodePdf[i] - equilibrium(temp, i, latticeVel);

      return true;
    }
  protected:
    const double *nearNodePdf;   //存储相邻格点上的分布函数 
    double boundaryVel[LC::GetDim()];  //边界给定的宏观量
    const double latticeVel;
  };
  
  /** \brief
    * NonEquilibriumOfPressure类：用于执行非平衡外推边界条件 \n
    * 给定密度（即压力）\n
    * 数据包括：<相邻格点的分布函数指针、相邻格点的宏观量指针、给定的边界宏观量、离散速度方向编号> \n
  */
  template<typename LC>
  class NonEquilibriumOfPressure : virtual public SimpleBC
  {
  public:
    /** \brief NonEquilibriumExtrapolation 类的构造函数:
     *\param iIn [in] 离散速度方向编号
     *\param nearNodePdfIn [in] 相邻格点的分布函数
     *\param nearNodeMacroIn [in] 相邻格点的宏观量
     *\param wVel [in] 边界给定的宏观速度
    */
    explicit NonEquilibriumOfPressure(
      const double *nearNodePdfIn,
      const double &wallPressureIn,
      const double latticeVelIn) :latticeVel(latticeVelIn), nearNodePdf(nearNodePdfIn), wallPressure(wallPressureIn) {}

    ~NonEquilibriumOfPressure() {}

    /** \brief 实现离散速度方向上的非平衡外推格式
     *\param pdfOut [out] 分布函数数组
     *\param pdfIn [in] 碰撞后的分布函数
     *\param macroVar [in] 宏观量
     *\return 返回bool型变量true
    */
    virtual bool operator ()(
      double *pdfOut,  // [out]分布函数数组
      const double *pdfIn,  // [in]碰撞后的分布函数
      const size_t &i
      ) const
    {
      nnwlbm::Equilibrium<LC> equilibrium;

      LC lc(latticeVel);
      double temp[LC::GetDim() + 1], boundaryVel[LC::GetDim() + 1];
      nnwlbm::cal_macro_var<LC>(temp, nearNodePdf, lc);

      double boundaryMacroVar[LC::GetDim() + 1];
      boundaryVel[0] = wallPressure;
      for (size_t i = 1; i <= LC::GetDim(); i++)
        boundaryVel[i] = temp[i];

      pdfOut[i] = equilibrium(boundaryVel, i, latticeVel) + nearNodePdf[i] - equilibrium(temp, i, latticeVel);

      return true;
    }
  protected:
    const double *nearNodePdf;   //存储相邻格点上的分布函数  
    double wallPressure;  //边界给定的宏观量
    const double latticeVel;
  };


  template<typename LC>
  class SingleNodeBoundaryCondition : virtual public SimpleBC {
  public:
    explicit SingleNodeBoundaryCondition() {}
    virtual ~SingleNodeBoundaryCondition() {}
    virtual bool operator ()(
      double *pdfOut,  // [out]分布函数数组
      const double *pdfIn,  // [in]碰撞后的分布函数
      double *macroVar // [in] 与边界相关的宏观量
      ) const {
      return false;
    }
  };

  /** \brief
    * PeriodicBoundary类：用于执行周期边界条件 \n
    * 数据包括：<周期边界对应格点上的分布函数、离散速度方向编号> \n
    * operator(): 调用以执行周期边界条件 \n
    * 示例，生成一个名为periodicBoundary的对象：
    ~~~~~{.cpp}
    nnwlbm::PeriodicBoundary<nnwlbm::D2Q9> periodicBoundary(i, nearNodePdf);
    ~~~~~
    * 示例，调用operator()执行周期边界条件：
    ~~~~~{.cpp}
     periodicBoundary(nodePdf, nodePdfPost, macro);
    ~~~~~
  */
  template<typename LC>
  class PeriodicBoundary : virtual public SingleNodeBoundaryCondition<LC>
  {
  public:
    /** \brief PeriodicBoundary 类的构造函数:
     *\param iIn [in] 离散速度方向编号
     *\param periodicNodePdfIn [in] 周期边界对应格点上的分布函数
    */
    explicit PeriodicBoundary(const size_t &iIn, const double *periodicNodePdfIn) :i(iIn), periodicNodePdf(periodicNodePdfIn)
    {}
    ~PeriodicBoundary() {}

    /** \brief 实现离散速度方向上的周期边界格式
     *\param pdfOut [out] 分布函数数组
     *\param pdfIn [in] 碰撞后的分布函数
     *\param macroVar [in] 格点上的宏观量
    */
    virtual bool operator ()(
      double *pdfOut,  // [out]分布函数数组
      const double *pdfIn,  // [in]碰撞后的分布函数
      double *macroVar // [in] 宏观量 
      ) const
    {
      pdfOut[i] = periodicNodePdf[i];
      return true;
    }

  protected:
    const double *periodicNodePdf; //周期边界对应格点上的分布函数
    const size_t i;  //离散速度方向编号
  };

  /** \brief
    * ConvectionBoundary类：用于执行充分对流边界条件 \n
    * 数据包括：<相邻格点上的分布函数、前一层格点上的最大速度、格子速度、离散速度方向编号> \n
    * operator(): 调用以执行充分对流边界条件 \n
    * 示例，生成一个名为convectionBoundary的对象：
    ~~~~~{.cpp}
    nnwlbm::ConvectionBoundary<nnwlbm::D2Q9> convectionBoundary(i,&velocityMaxx, nearNodePdf,lc.velocity());
    ~~~~~
    * 示例，调用operator()执行充分对流边界条件：
    ~~~~~{.cpp}
     convectionBoundary(nodePdf, nodePdfPost, macro);
    ~~~~~
  */
  template<typename LC>
  class ConvectionBoundary : virtual public SingleNodeBoundaryCondition<LC>
  {
  public:
    /** \brief ConvectionBoundary 类的构造函数:
     *\param iIn [in] 离散速度方向编号
     *\param velocityMaxxIn [in] 最大速度
     *\param nearNodePdfIn [in] 相邻格点的分布函数
     *\param latticeVelocityIn [in] 格子速度
    */
    explicit ConvectionBoundary(const size_t &iIn, const double *velocityMaxxIn, const double *nearNodePdfIn, const double &latticeVelocityIn) :i(iIn), velocityMaxx(velocityMaxxIn), nearNodePdf(nearNodePdfIn), latticeVelocity(latticeVelocityIn)
    {}

    ~ConvectionBoundary() {}

    /** \brief 实现离散速度方向上的充分对流边界条件
     *\param pdfOut [out] 分布函数数组
     *\param pdfIn [in] 碰撞后的分布函数
     *\param macroVar [in] 格点上的宏观量
     *\return 返回bool型变量true
    */
    virtual bool operator ()(
      double *pdfOut,  //[out]分布函数数组
      const double *pdfIn,  //[in]碰撞后的分布函数
      double *macroVar // [in] 宏观量
      ) const
    {
      double temp = (*velocityMaxx);
      pdfOut[i] = (pdfOut[i] + temp * nearNodePdf[i] / latticeVelocity) / (1.0 + temp / latticeVelocity);
      return true;
    }

  protected:
    const double *nearNodePdf; //相邻格点上的分布函数
    const double *velocityMaxx; //前一层格点上的最大速度
    const double latticeVelocity; //格子速度
    const size_t i;  //离散速度方向
  };

  /** \brief
    * NonEquilibriumExtrapolation1类：用于执行非平衡外推边界条件 \n
    * 数据包括：<相邻格点的分布函数指针、相邻格点的宏观量指针、给定的边界宏观量、离散速度方向编号> \n
    * operator(): 调用以执行非平衡外推边界条件 \n
    * 示例，生成一个名为nonEquilibriumExtrapolation1的对象：
    ~~~~~{.cpp}
    nnwlbm::NonEquilibriumExtrapolation1<nnwlbm::D2Q9> nonEquilibriumExtrapolation1(i, nearNodePdf, nearNodeMacro, nodeMacro);
    ~~~~~
    * 示例，调用operator()执行非平衡外推边界条件：
    ~~~~~{.cpp}
     nonEquilibriumExtrapolation1(nodePdf, nodePdfPost, macro);
    ~~~~~
  */
  template<typename LC>
  class NonEquilibriumExtrapolation1 : virtual public SingleNodeBoundaryCondition<LC>
  {
  public:
    /** \brief NonEquilibriumExtrapolation 类的构造函数:
     *\param iIn [in] 离散速度方向编号
     *\param nearNodePdfIn [in] 相邻格点的分布函数
     *\param nearNodeMacroIn [in] 相邻格点的宏观量
     *\param wVel [in] 边界给定的宏观速度
    */
    explicit NonEquilibriumExtrapolation1(const int &iIn, const double *nearNodePdfIn, const double *nearNodeMacroIn, const double *wVel) :i(iIn), nearNodePdf(nearNodePdfIn), nearNodeMacro(nearNodeMacroIn)
    {
      for (size_t i = 0; i < LC::GetDim() + 1; i++) boundaryVel[i] = wVel[i];
    }

    ~NonEquilibriumExtrapolation1() {}

    /** \brief 实现离散速度方向上的非平衡外推格式
     *\param pdfOut [out] 分布函数数组
     *\param pdfIn [in] 碰撞后的分布函数
     *\param macroVar [in] 宏观量
     *\return 返回bool型变量true
    */
    virtual bool operator ()(
      double *pdfOut,  // [out]分布函数数组
      const double *pdfIn,  // [in]碰撞后的分布函数
      double *macroVar // [in] 宏观量 
      ) const
    {
      for (size_t j = 0; j < LC::GetDim() + 1; j++) {
        macroVar[j] = boundaryVel[j];
      }
      pdfOut[i] = nnwlbm::Equilibrium<LC>()(macroVar, i, 1.0) + nearNodePdf[i] - nnwlbm::Equilibrium<LC>()(nearNodeMacro, i, 1.0);
      return true;
    }
  protected:
    const double *nearNodePdf;   //存储相邻格点上的分布函数 
    const double *nearNodeMacro;  //存储相邻格点上的宏观量 
    double boundaryVel[LC::GetDim() + 1];  //边界给定的宏观量
    const size_t i;  //离散速度方向
  };
}

#endif // !BOUNDARYCONDITION_H

