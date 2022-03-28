///////////////////////////////////////////////////////////////
// author : Wang,Junquan. email: as2697920806@gmail.com
///////////////////////////////////////////////////////////////
// participant 1 : Chen,Songze. email: 
///////////////////////////////////////////////////////////////
// tester 1 : 
///////////////////////////////////////////////////////////////

#pragma once
#ifndef BOUNDARYFACTORY_H
#define BOUNDARYFACTORY_H

#include "CommonAPI.h"
#include "Lattice.h"
#include "Range.h"
#include "BoundaryCondition.h"
#include "Geometry.h" 
#include "UniformGrid.h"

namespace nnwlbm
{
  template <typename LC>
  class BoundaryFactory
  {
  public:
    explicit BoundaryFactory() {}

    virtual std::shared_ptr<SimpleBC> CreateBoundary(
      const int *locIn,
      const int &i,
      const nnwlbm::LinkedMaskData<double, LC::GetDim()> &pdfData) const = 0;

    virtual ~BoundaryFactory() {}
  };

  /** \brief
    * BounceBackBc的特例化BoundaryFactory类：用于生成标准反弹类 \n
    * 数据包括：<格子类型、平衡态、计算区域内的宏观量、网格、分布函数、碰撞后的分布函数> \n
    * CreateBoundary()函数: 调用以生成标准反弹类 \n
    * 示例，生成一个名为bounceBackFactory的对象：
    ~~~~~{.cpp}
    auto bounceBackFactory = std::make_shared<nnwlbm::BoundaryFactory<Ltc,EQ, nnwlbm::SingleNodeBoundaryCondition<Ltc>, nnwlbm::BounceBackBc<Ltc>>>(lc, equilibrium,maskData,grid,pdfData, pdfDataPost);
    ~~~~~
    * 示例，调用CreateBoundary()生成标准反弹类：
    ~~~~~{.cpp}
     std::shared_ptr<nnwlbm::SingleNodeBoundaryCondition<Ltc>> bounceBack[Ltc::size()];
     for (size_t i = 0; i < Ltc::size(); i++)
      bounceBack[i] = bounceBackFactory->CreateBoundary(loc, i);
    ~~~~~
  */
  template<typename LC>
  class BounceBackFactory :virtual public BoundaryFactory<LC>
  {
  public:
    /** \brief BoundaryFactory类的构造函数:
     *\param lcIn [in] 格子类型
     *\param equilibriumIn [in] 平衡态
     *\param macroDataIn [in] 计算区域的宏观量
     *\param gridIn [in] 网格
     *\param pdfDataOut [in] 分布函数
     *\param pdfDataIn [in] 碰撞后的分布函数
    */
    explicit BounceBackFactory() {}

    ~BounceBackFactory() {}

  public:
    /** \brief 实现输入位置和离散速度方向编号，生成标准反弹类
     *\param locIn [out] 格点位置
     *\param i [in] 离散速度方向编号
     *\return 返回一个智能指针，指向生成的标准反弹类
    */
    virtual std::shared_ptr<SimpleBC> CreateBoundary(
      const int *locIn,
      const int &i,
      const nnwlbm::LinkedMaskData<double, LC::GetDim()> &pdfData
    ) const
    {
      return std::make_shared<nnwlbm::BounceBackBc<LC>>();
    }
  };

  /** \brief
    * MovingBounceBackBc的特例化BoundaryFactory类：用于生成MovingBounceBackBc类 \n
    * 数据包括：<格子类型、平衡态、计算区域内的宏观量、网格、分布函数、碰撞后的分布函数、移动边界的宏观速度> \n
    * CreateBoundary()函数: 调用以生成MovingBounceBackBc类 \n
    * 示例，生成一个名为movingBounceBackFactory的对象：
    ~~~~~{.cpp}
    auto movingBounceBackFactory = std::make_shared<nnwlbm::BoundaryFactory<Ltc, EQ, nnwlbm::SingleNodeBoundaryCondition<Ltc>, nnwlbm::MovingBounceBackBc<Ltc>>>(lc, equilibrium, maskData, grid, pdfData, pdfDataPost,macro);
    ~~~~~
    * 示例，调用CreateBoundary()生成MovingBounceBackBc类：
    ~~~~~{.cpp}
     std::shared_ptr<nnwlbm::SingleNodeBoundaryCondition<Ltc>> movingBounceBack[Ltc::size()];
    for (size_t i = 0; i < Ltc::size(); i++)
      movingBounceBack[i] = movingBounceBackFactory->CreateBoundary(loc, i);
    ~~~~~
  */
  template<typename LC>
  class MovingBounceBackFactory :virtual public BoundaryFactory<LC>
  {
    const double latticeVel;
    double macroVar[LC::GetDim()];

  public:
    /** \brief BoundaryFactory类的构造函数:
     *\param lcIn [in] 格子类型
     *\param equilibriumIn [in] 平衡态
     *\param macroDataIn [in] 计算区域的宏观量
     *\param gridIn [in] 网格
     *\param pdfDataOut [in] 分布函数
     *\param pdfDataIn [in] 碰撞后的分布函数
     *\param macroVarIn [in] 移动边界上的宏观速度
    */
    explicit MovingBounceBackFactory(
      const double &latticeVelIn,
      const double *macroVarIn
    ) :latticeVel(latticeVelIn) {
      for (int i = 0; i < LC::GetDim(); i++)
        macroVar[i] = macroVarIn[i];
    }

    ~MovingBounceBackFactory() {}

    /** \brief 实现输入位置和离散速度方向编号，生成标准反弹类
     *\param locIn [out] 格点位置
     *\param i [in] 离散速度方向编号
     *\return 返回生成的MovingBounceBackBc类的智能指针
    */
    virtual std::shared_ptr<SimpleBC> CreateBoundary(
      const int *locIn,
      const int &i,
      const nnwlbm::LinkedMaskData<double, LC::GetDim()> &pdfData
    ) const {
      int nodeId;
      nnwlbm::RangeIndex<LC::GetDim()> rIndex(pdfData.GetRange(), nnwlbm::RangeIndex<LC::GetDim()>::inputNodeRange);

      nodeId = rIndex.CalNodeID(locIn);
      const double *point = pdfData.GetVar(nodeId);

      double macroRho = 0.0;
      for (int j = 0; j < LC::size(); j++)
      {
        double temp = *(point + j);
        macroRho = macroRho + temp;
      }
      return std::make_shared<nnwlbm::MovingBounceBackBc<LC>>(macroRho, latticeVel, macroVar);
    }
  };

  /** \brief
    * NonEquilibriumExtrapolation的特例化BoundaryFactory类：用于生成NonEquilibriumExtrapolation类 \n
    * 数据包括：<格子类型、平衡态、计算区域内的宏观量、网格、分布函数、碰撞后的分布函数、边界给定的宏观速度> \n
    * CreateBoundary()函数: 调用以生成NonEquilibriumExtrapolation类 \n
    * 示例，生成一个名为nonEquilibriumExtrapolationFactory的对象：
    ~~~~~{.cpp}
    auto nonEquilibriumExtrapolationFactory = std::make_shared<nnwlbm::BoundaryFactory<Ltc, EQ, nnwlbm::SingleNodeBoundaryCondition<Ltc>, nnwlbm::NonEquilibriumExtrapolation<Ltc>>>(lc, equilibrium, maskData, grid, pdfData, pdfDataPost, macro);
    ~~~~~
    * 示例，调用CreateBoundary()生成NonEquilibriumExtrapolation类：
    ~~~~~{.cpp}
     std::shared_ptr<nnwlbm::SingleNodeBoundaryCondition<Ltc>> nonEquilibriumExtrapolation[Ltc::size()];
     for (size_t i = 0; i < Ltc::size(); i++)
       nonEquilibriumExtrapolation[i] = nonEquilibriumExtrapolationFactory->CreateBoundary(loc, i);
    ~~~~~
  */
  template<typename LC>
  class NonEquilibriumExtrapolationFactory :virtual public BoundaryFactory<LC>
  {
    const double &latticeVel;
    double macroVar[LC::GetDim()];
    int vector[LC::GetDim()];

  public:
    /** \brief BoundaryFactory类的构造函数:
     *\param lcIn [in] 格子类型
     *\param equilibriumIn [in] 平衡态
     *\param macroDataIn [in] 计算区域的宏观量
     *\param gridIn [in] 网格
     *\param pdfDataOut [in] 分布函数
     *\param pdfDataIn [in] 碰撞后的分布函数
     *\param macroVarIn [in] 边界给定的宏观速度
    */
    explicit NonEquilibriumExtrapolationFactory(
      const double &latticeVelIn,
      const double *macroVarIn,
      const int *vectorIn
    ) :latticeVel(latticeVelIn) {
      for (int i = 0; i < LC::GetDim(); i++) macroVar[i] = macroVarIn[i];
      for (int i = 0; i < LC::GetDim(); i++) vector[i] = vectorIn[i];
    }

    ~NonEquilibriumExtrapolationFactory() {}

    /** \brief 实现输入位置和离散速度方向编号，生成标准反弹类
     *\param locIn [out] 格点位置
     *\param i [in] 离散速度方向编号
     *\return 返回生成的NonEquilibriumExtrapolation类的智能指针
    */
    virtual std::shared_ptr<SimpleBC> CreateBoundary(
      const int *locIn,
      const int &i,
      const nnwlbm::LinkedMaskData<double, LC::GetDim()> &pdfData
    ) const {
      int locN[LC::GetDim()];
      int nodeIdN;

      nnwlbm::RangeIndex<LC::GetDim()> rIndex(pdfData.GetRange(), nnwlbm::RangeIndex<LC::GetDim()>::inputNodeRange);

      for (size_t j = 0; j < LC::GetDim(); j++)
      {
        locN[j] = locIn[j] + vector[j];
      }
      nodeIdN = rIndex.CalNodeID(locN);

      return std::make_shared<nnwlbm::NonEquilibriumExtrapolation<LC>>(pdfData.GetVar(nodeIdN), macroVar, latticeVel);
    }
  };

  /** \brief
    * NonEquilibriumOfPressure的特例化BoundaryFactory类：用于生成NonEquilibriumOfPressure类 \n
    * 数据包括：<格子类型、平衡态、计算区域内的宏观量、网格、分布函数、碰撞后的分布函数、边界给定的宏观速度> \n
    * CreateBoundary()函数: 调用以生成NonEquilibriumExtrapolation类 \n
    * 示例，生成一个名为nonEquilibriumExtrapolationFactory的对象：
    ~~~~~{.cpp}
    auto nonEquilibriumExtrapolationFactory = std::make_shared<nnwlbm::BoundaryFactory<Ltc, EQ, nnwlbm::SingleNodeBoundaryCondition<Ltc>, nnwlbm::NonEquilibriumExtrapolation<Ltc>>>(lc, equilibrium, maskData, grid, pdfData, pdfDataPost, macro);
    ~~~~~
    * 示例，调用CreateBoundary()生成NonEquilibriumExtrapolation类：
    ~~~~~{.cpp}
     std::shared_ptr<nnwlbm::SingleNodeBoundaryCondition<Ltc>> nonEquilibriumExtrapolation[Ltc::size()];
     for (size_t i = 0; i < Ltc::size(); i++)
       nonEquilibriumExtrapolation[i] = nonEquilibriumExtrapolationFactory->CreateBoundary(loc, i);
    ~~~~~
  */
  template<typename LC>
  class NonEquilibriumOfPressureFactory :virtual public BoundaryFactory<LC>
  {
    const double latticeVel;
    const double pressure;
    int vector[LC::GetDim()];

  public:
    /** \brief BoundaryFactory类的构造函数:
     *\param lcIn [in] 格子类型
     *\param equilibriumIn [in] 平衡态
     *\param macroDataIn [in] 计算区域的宏观量
     *\param gridIn [in] 网格
     *\param pdfDataOut [in] 分布函数
     *\param pdfDataIn [in] 碰撞后的分布函数
     *\param macroVarIn [in] 边界给定的宏观速度
    */
    explicit NonEquilibriumOfPressureFactory(
      const double &latticeVelIn,
      const double &pressureIn,
      const int *vectorIn
    ) :latticeVel(latticeVelIn), pressure(pressureIn) {
      for (int i = 0; i < LC::GetDim(); i++) vector[i] = vectorIn[i];
    }

    ~NonEquilibriumOfPressureFactory() {}

    /** \brief 实现输入位置和离散速度方向编号，生成标准反弹类
     *\param locIn [out] 格点位置
     *\param i [in] 离散速度方向编号
     *\return 返回生成的NonEquilibriumExtrapolation类的智能指针
    */
    virtual std::shared_ptr<SimpleBC> CreateBoundary(
      const int *locIn,
      const int &i,
      const nnwlbm::LinkedMaskData<double, LC::GetDim()> &pdfData
    ) const {
      int locN[LC::GetDim()];
      int nodeIdN;

      nnwlbm::RangeIndex<LC::GetDim()> rIndex(pdfData.GetRange(), nnwlbm::RangeIndex<LC::GetDim()>::inputNodeRange);

      for (size_t j = 0; j < LC::GetDim(); j++)
      {
        locN[j] = locIn[j] + vector[j];
      }
      nodeIdN = rIndex.CalNodeID(locN);

      return std::make_shared<nnwlbm::NonEquilibriumOfPressure<LC>>(pdfData.GetVar(nodeIdN), pressure, latticeVel);
    }
  };



  /** \brief
    * PeriodicBoundary的特例化BoundaryFactory类：用于生成PeriodicBoundary类 \n
    * 数据包括：<格子类型、平衡态、计算区域内的宏观量、网格、分布函数、碰撞后的分布函数、周期边界对应的格点位置> \n
    * CreateBoundary()函数: 调用以生成PeriodicBoundary类 \n
    * 示例，生成一个名为periodicBoundaryFactory的对象：
    ~~~~~{.cpp}
    auto periodicBoundaryFactory = std::make_shared<nnwlbm::BoundaryFactory<Ltc, EQ, nnwlbm::SingleNodeBoundaryCondition<Ltc>, nnwlbm::PeriodicBoundary<Ltc>>>(lc, equilibrium, maskData, grid, pdfData, pdfDataPost, loc);
    ~~~~~
    * 示例，调用CreateBoundary()生成PeriodicBoundary类：
    ~~~~~{.cpp}
     std::shared_ptr<nnwlbm::SingleNodeBoundaryCondition<Ltc>> periodicBoundary[Ltc::size()];
     for (size_t i = 0; i < Ltc::size(); i++)
      periodicBoundary[i] = periodicBoundaryFactory->CreateBoundary(loc, i);
    ~~~~~
  */
  template<typename LC, typename EQ>
  class PeriodicBoundaryFactory :virtual public BoundaryFactory<LC>
  {
    const LC &lc;
    const EQ &equilibrium;
    const nnwlbm::LinkedMaskData<double, LC::GetDim()> &macroData;
    const nnwlbm::UniformGrid<LC::GetDim()> &grid;
    nnwlbm::LinkedMaskData<double, LC::GetDim()> &pdfData;
    const nnwlbm::LinkedMaskData<double, LC::GetDim()> &pdfDataPost;
    const int *loc;

  public:
    /** \brief BoundaryFactory类的构造函数:
     *\param lcIn [in] 格子类型
     *\param equilibriumIn [in] 平衡态
     *\param macroDataIn [in] 计算区域的宏观量
     *\param gridIn [in] 网格
     *\param pdfDataOut [in] 分布函数
     *\param pdfDataIn [in] 碰撞后的分布函数
     *\param locIn [in] 周期边界对应的格点位置
    */
    explicit PeriodicBoundaryFactory(
      const LC &lcIn,
      const EQ &equilibriumIn,
      nnwlbm::LinkedMaskData<double, LC::GetDim()> &macroDataIn,
      const nnwlbm::UniformGrid<LC::GetDim()> &gridIn,
      nnwlbm::LinkedMaskData<double, LC::GetDim()> &pdfDataOut,
      const nnwlbm::LinkedMaskData<double, LC::GetDim()> &pdfDataIn,
      const int *locIn
    ) :lc(lcIn), equilibrium(equilibriumIn), macroData(macroDataIn), grid(gridIn), pdfData(pdfDataOut), pdfDataPost(pdfDataIn), loc(locIn) {}

    ~PeriodicBoundaryFactory() {}

    /** \brief 实现输入位置和离散速度方向编号，生成标准反弹类
     *\param locIn [out] 格点位置
     *\param i [in] 离散速度方向编号
     *\return 返回生成的PeriodicBoundary类的智能指针
    */
    virtual std::shared_ptr<SingleNodeBoundaryCondition<LC>> CreateBoundary(const int *locIn, const int &i) const {
      nnwlbm::RangeIndex<LC::GetDim()> rIndex(macroData.GetRange(), nnwlbm::RangeIndex<LC::GetDim()>::inputNodeRange);
      size_t nodeId = rIndex.CalNodeID(loc);

      return std::make_shared<nnwlbm::PeriodicBoundary<LC>>(i, pdfData.GetVar(nodeId));
    }
  };

  /** \brief
    * ConvectionBoundary的特例化BoundaryFactory类：用于生成ConvectionBoundary类 \n
    * 数据包括：<格子类型、平衡态、计算区域内的宏观量、网格、分布函数、碰撞后的分布函数、前一层格点上的最大速度> \n
    * CreateBoundary()函数: 调用以生成ConvectionBoundary类 \n
    * 示例，生成一个名为periodicBoundaryFactory的对象：
    ~~~~~{.cpp}
    auto convectionBoundaryFactory = std::make_shared<nnwlbm::BoundaryFactory<Ltc, EQ, nnwlbm::SingleNodeBoundaryCondition<Ltc>, nnwlbm::ConvectionBoundary<Ltc>>>(lc, equilibrium, maskData, grid, pdfData, pdfDataPost,&velocityMaxx);
    ~~~~~
    * 示例，调用CreateBoundary()生成ConvectionBoundary类：
    ~~~~~{.cpp}
     std::shared_ptr<nnwlbm::SingleNodeBoundaryCondition<Ltc>> convectionBoundary[Ltc::size()];
     for (size_t i = 0; i < Ltc::size(); i++)
      convectionBoundary[i] = convectionBoundaryFactory->CreateBoundary(loc, i);
    ~~~~~
  */
  template<typename LC, typename EQ>
  class ConvectionBoundaryFactory :virtual public BoundaryFactory<LC>
  {
    const LC &lc;
    const EQ &equilibrium;
    const nnwlbm::LinkedMaskData<double, LC::GetDim()> &macroData;
    const nnwlbm::UniformGrid<LC::GetDim()> &grid;
    nnwlbm::LinkedMaskData<double, LC::GetDim()> &pdfData;
    const nnwlbm::LinkedMaskData<double, LC::GetDim()> &pdfDataPost;
    const int *vector;

  public:
    /** \brief BoundaryFactory类的构造函数:
     *\param lcIn [in] 格子类型
     *\param equilibriumIn [in] 平衡态
     *\param macroDataIn [in] 计算区域的宏观量
     *\param gridIn [in] 网格
     *\param pdfDataOut [in] 分布函数
     *\param pdfDataIn [in] 碰撞后的分布函数
     *\param velocityMaxxIn [in] 最大速度
    */
    explicit ConvectionBoundaryFactory(
      const LC &lcIn,
      const EQ &equilibriumIn,
      nnwlbm::LinkedMaskData<double, LC::GetDim()> &macroDataIn,
      const nnwlbm::UniformGrid<LC::GetDim()> &gridIn,
      nnwlbm::LinkedMaskData<double, LC::GetDim()> &pdfDataOut,
      const nnwlbm::LinkedMaskData<double, LC::GetDim()> &pdfDataIn,
      const int *vectorIn
    ) :lc(lcIn), equilibrium(equilibriumIn), macroData(macroDataIn), grid(gridIn), pdfData(pdfDataOut), pdfDataPost(pdfDataIn), vector(vectorIn){}

    ~ConvectionBoundaryFactory() {}

    /** \brief 实现输入位置和离散速度方向编号，生成标准反弹类
     *\param locIn [out] 格点位置
     *\param i [in] 离散速度方向编号
     *\return 返回生成的ConvectionBoundary类的智能指针
    */
    virtual std::shared_ptr<SingleNodeBoundaryCondition<LC>> CreateBoundary(const int *locIn, const int &i) const {
      int locN[LC::GetDim()];
      int nodeIdN; int flagN;

      nnwlbm::RangeIndex<LC::GetDim()> rIndex(macroData.GetRange(), nnwlbm::RangeIndex<LC::GetDim()>::inputNodeRange);

      for (size_t j = 0; j < LC::GetDim(); j++)
      {
        locN[j] = locIn[j] + vector[j];
      }
      nodeIdN = rIndex.CalNodeID(locN);
      return std::make_shared<nnwlbm::ConvectionBoundary<LC>>(i, macroData.GetVar(nodeIdN) + 1, pdfData.GetVar(nodeIdN), lc.velocity());
    }
  };

  /** \brief
    * NonEquilibriumExtrapolation1的特例化BoundaryFactory类：用于生成NonEquilibriumExtrapolation1类 \n
    * 数据包括：<格子类型、平衡态、计算区域内的宏观量、网格、分布函数、碰撞后的分布函数、边界给定的宏观速度> \n
    * CreateBoundary()函数: 调用以生成NonEquilibriumExtrapolation1类 \n
    * 示例，生成一个名为nonEquilibriumExtrapolation1Factory的对象：
    ~~~~~{.cpp}
    auto nonEquilibriumExtrapolation1Factory = std::make_shared<nnwlbm::BoundaryFactory<Ltc, EQ, nnwlbm::SingleNodeBoundaryCondition<Ltc>, nnwlbm::NonEquilibriumExtrapolation1<Ltc>>>(lc, equilibrium, maskData, grid, pdfData, pdfDataPost, macro);
    ~~~~~
    * 示例，调用CreateBoundary()生成NonEquilibriumExtrapolation类：
    ~~~~~{.cpp}
     std::shared_ptr<nnwlbm::SingleNodeBoundaryCondition<Ltc>> nonEquilibriumExtrapolation1[Ltc::size()];
     for (size_t i = 0; i < Ltc::size(); i++)
       nonEquilibriumExtrapolation1[i] = nonEquilibriumExtrapolation1Factory->CreateBoundary(loc, i);
    ~~~~~
  */
  template<typename LC, typename EQ>
  class NonEquilibriumExtrapolationFactory1 :virtual public BoundaryFactory<LC>
  {
    const LC &lc;
    const EQ &equilibrium;
    const nnwlbm::LinkedMaskData<double, LC::GetDim()> &macroData;
    const nnwlbm::UniformGrid<LC::GetDim()> &grid;
    nnwlbm::LinkedMaskData<double, LC::GetDim()> &pdfData;
    const nnwlbm::LinkedMaskData<double, LC::GetDim()> &pdfDataPost;
    const double *macroVar;
    const int *vector;

  public:
    /** \brief BoundaryFactory类的构造函数:
     *\param lcIn [in] 格子类型
     *\param equilibriumIn [in] 平衡态
     *\param macroDataIn [in] 计算区域的宏观量
     *\param gridIn [in] 网格
     *\param pdfDataOut [in] 分布函数
     *\param pdfDataIn [in] 碰撞后的分布函数
     *\param macroVarIn [in] 边界给定的宏观速度
    */
    explicit NonEquilibriumExtrapolationFactory1(
      const LC &lcIn,
      const EQ &equilibriumIn,
      nnwlbm::LinkedMaskData<double, LC::GetDim()> &macroDataIn,
      const nnwlbm::UniformGrid<LC::GetDim()> &gridIn,
      nnwlbm::LinkedMaskData<double, LC::GetDim()> &pdfDataOut,
      const nnwlbm::LinkedMaskData<double, LC::GetDim()> &pdfDataIn,
      double *macroVarIn,
      const int *vectorIn
    ) :lc(lcIn), equilibrium(equilibriumIn), macroData(macroDataIn), grid(gridIn), pdfData(pdfDataOut), pdfDataPost(pdfDataIn), macroVar(macroVarIn),vector(vectorIn) {}

    ~NonEquilibriumExtrapolationFactory1() {}

    /** \brief 实现输入位置和离散速度方向编号，生成标准反弹类
     *\param locIn [out] 格点位置
     *\param i [in] 离散速度方向编号
     *\return 返回生成的NonEquilibriumExtrapolation1类的智能指针
    */
    virtual std::shared_ptr<SingleNodeBoundaryCondition<LC>> CreateBoundary(const int *locIn, const int &i) const {
      int locN[LC::GetDim()];
      int nodeIdN;

      nnwlbm::RangeIndex<LC::GetDim()> rIndex(macroData.GetRange(), nnwlbm::RangeIndex<LC::GetDim()>::inputNodeRange);

      for (size_t j = 0; j < LC::GetDim(); j++)
      {
        locN[j] = locIn[j] + vector[j];
      }
      nodeIdN = rIndex.CalNodeID(locN);
      return std::make_shared<nnwlbm::NonEquilibriumExtrapolation1<LC>>(i, pdfData.GetVar(nodeIdN), macroData.GetVar(nodeIdN), macroVar);
    }
  };
}
#endif