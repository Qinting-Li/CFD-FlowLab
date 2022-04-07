///////////////////////////////////////////////////////////////
// author : Wang,Junquan. email: as2697920806@gmail.com
///////////////////////////////////////////////////////////////
// tester 1 : Wang,Junquan
///////////////////////////////////////////////////////////////

#pragma once
#ifndef COLLISION_H
#define COLISSION_H

#include "CommonAPI.h"
#include "Range.h"
#include "UniformGrid.h"
#include "D2Q9.h"
#include "D3Q15.h"
#include "D3Q19.h"
//#include "SimpleBC.h"
#include <iostream>
#include <vector>


/** \brief 配合D2Q9BGK或者D3Q15BGK使用的计算碰撞的模板函数 
  *实现格点的碰撞模块，根据输入参数，执行碰撞，得到碰撞后的分布函数 
  * \param lattic [in] 格子构型，对应的是nnwlbm::D2Q9BGK 
  * \param macro [in] 存储宏观变量的数组 
  * \param pdf [in] 和格子离散速度相对应的分布函数
  * \param pdf_post [in,out] 储存碰撞后的分布函数
  * \param tau [in] 松弛时间，通过流体相关物理量求得
  */
template<class LC>
void bgk_collision (
  LC lattic, /**< 格子构型 */
  double *macro, /**< 宏观量，存在一个double macro[3]的数组里面，第一个元素表示u，第二个元素表示v，第三个元素表示rho */
  double *pdf, /**< 格点位置处的分布函数：D2Q9对应九个方向上的分布函数 */
  double *pdf_post, /**< 存储碰撞后的分布函数 */
  const double &tau /**< 松弛时间 */
)
{
  for (int i = 0; i < lattic.size(); i++)
    pdf_post[i] = pdf[i] - 1 / tau * (pdf[i] - lattic.CalEqui(macro, i));
}

/** \brief 配合D2Q9MRT或者D3Q15MRT使用的计算碰撞的模板函数
  *实现MRT模型的碰撞模块，根据输入参数，执行碰撞，得到碰撞后的分布函数
  * \param lattic [in] 格子构型，对应的是nnwlbm::D2Q9MRT
  * \param macro [in] 存储宏观变量的数组
  * \param pdf [in] 格子离散速度相对应的分布函数
  * \param pdf_post [in,out] 储存碰撞后的分布函数
  * \param s [in] 松弛矩阵，通过流体相关物理量求得，对应九个参数
  */
template<class LC>
void mrt_collision
(
  const LC &lattice, /**< [in]格子构型 */
  double *macro, /**< [in]宏观量，存在一个double macro[3]的数组里面，第一个元素表示u，第二个元素表示v，第三个元素表示rho */
  const double *pdf, /**< [in]格点位置处的分布函数：D2Q9对应九个方向上的分布函数 */
  double *pdfPost, /**< [in,out]存储碰撞后的分布函数 */
  const double *s /**< [in]MRT模型对应的松弛矩阵 */
)
{
  double moment[LC::size()];
  double momentEq[LC::size()];
  nnwlbm::pdf_to_moment<LC>(moment, pdf);

  double pdfEq[LC::size()];
  for(size_t i = 0; i < LC::size(); i++) pdfEq[i] = lattice.CalEqui(macro, i);
  nnwlbm::pdf_to_moment<LC>(momentEq, pdfEq);

  for (size_t i = 0; i < LC::size(); i++)
  {
    moment[i] = moment[i] - s[i] * (moment[i] - momentEq[i]);
  }

  nnwlbm::MomentToPdf<LC>(pdfPost, moment);
}


/** \brief 配合D2Q9格子构型使用的迁移函数
  *实现分布函数的迁移功能，输入格点(假设为A)的位置，然后将周围格点上的碰撞后的分布函数按离散速度迁移到格点A
  * \param lattice [in] 格子构型，对应的是D2Q9,D2Q9BGK,D2Q9MRT
  * \param md maskData [in] 对象，存储格点处的相关数据，如colormap
  * \param loc [in] 格点在网格上的坐标
  * \param pdf1 [in,out] 分布函数数组
  * \param pdf0 [in] 碰撞后的分布函数数组
  * 流体节点为0，边界为1，固体障碍物为2
  */
//template<class LC, const size_t dim>
//bool streaming
//(
//  LC lattice,/**< [in]格子构型 */ 
//  nnwlbm::MaskData<double, dim>& md, /**< [in]存储流场数据 */ 
//  const int *loc, /**< [in]格点位置坐标 */ 
//  double *pdf1, /**< [in,out]格点位置处的分布函数 */ 
//  const double *pdf0 /**< [in]碰撞后的分布函数 */ 
//)
//{
//  //迁移前格点的位置，序号，flag 
//  int locN[dim];
//  size_t nodeIdN;
//  int flagN;
//  //根据输入的maskData，创建一个rangindex对象，用来确定nodeId 
//  nnwlbm::RangeIndex<dim> rIndex(md.GetRange(), nnwlbm::RangeIndex<dim>::inputNodeRange);
//  size_t nodeId = rIndex.CalNodeID(loc);
//  nnwlbm::BounceBack<LC> bc;
//
//  int flag = md.GetColor(nodeId); 
//  double* cpdf= pdf1 + lattice.size() * rIndex.CalNodeID(loc);
//
//  switch (flag)
//  {
//    case 0:
//    {
//      for (size_t i = 0; i < lattice.size(); i++)
//      {
//        for (size_t j = 0; j < dim; j++)
//        {
//          locN[j] = loc[j] - *(lattice.direction(i) + j);
//        }
//        nodeIdN = rIndex.CalNodeID(locN);
//        flagN = md.GetColor(nodeIdN);
//        cpdf[i] = pdf0[i + lattice.size() * nodeIdN]; 
//      }
//      break;
//    }
//    case 1: 
//    {
//      for (int i = 0; i < lattice.size(); i++)
//      {
//        for (int j = 0; j < dim; j++)
//        {
//          locN[j] = loc[j] - *(lattice.direction(i) + j);
//        }
//        size_t nodeIdN = rIndex.CalNodeID(locN);
//        int flagN = md.GetColor(nodeIdN);
//        if (flagN == 2)
//        {
//          bc(pdf1 + lattice.size()*nodeId, pdf0 + lattice.size()*nodeId, i, nullptr);
//        }
//        else
//        {
//          cpdf[i] = pdf0[i + lattice.size() * nodeIdN];
//        }
//      }
//      break;
//    }
//    case 2:
//    { 
//      break; 
//    }
//    default:
//    {
//      std::cout << "No." << nodeId << ":flag error!" << std::endl;
//      break;
//    }
//  }
//  return true;
//}

namespace nnwlbm {
  template<typename LC>
  double CalTau(const LC &lc, const double &vis, const double &dt, const double &rho)
  {
    return 3.0 * vis/(rho * lc.velocity() * lc.velocity()) / dt + 0.5;
  }
  
  //！DnQm平衡态分布函数计算，目前模块包含的D2Q9,D3Q15,D3Q19的平衡态分布函数
  template<typename LC>
  class Equilibrium {
  public:
    template<typename LC1 = LC>
    typename std::enable_if<::std::is_base_of<Lattice<2>, LC1>::value, double>::type 
    operator() (
      const double *macro,
      const size_t i,
      const double &ltcVel = 1.0
    ) const {
      const double ev = macro[1] * LC::direction(i)[0] + macro[2] * LC::direction(i)[1];
      const double vv = macro[2] * macro[2] + macro[1] * macro[1];
      return LC::weight(i) * macro[0] * (1.0 + (3.0 * ev + 1.5 * (3.0 * ev * ev - vv) / ltcVel) / ltcVel);
    }

    template<typename LC1 = LC>
    typename std::enable_if<::std::is_base_of<Lattice<2>, LC1>::value, void>::type 
    operator() (
      double *pdf, 
      const double *macro,
      const double &ltcVel = 1.0
    ) const {
      const double vv = macro[1] * macro[1] + macro[2] * macro[2];
      for(size_t i = 0; i < LC::size(); i++) {
        const double ev = macro[1] * LC::direction(i)[0] + macro[2] * LC::direction(i)[1];
        pdf[i] = LC::weight(i) * macro[0] * (1.0 + (3.0 * ev + 1.5 * (3.0 * ev * ev - vv) / ltcVel) / ltcVel);
      }
    }

    template<typename LC1 = LC>
    typename std::enable_if<::std::is_base_of<Lattice<3>, LC1>::value, double>::type  operator() (
      const double *macro,
      const size_t i,
      const double &ltcVel = 1.0
    ) const {
      const double ev = macro[1] * (*(LC::direction(i))) + macro[2] * (*(LC::direction(i) + 1))
        + macro[3] * (*(LC::direction(i) + 2));
      const double vv = macro[3] * macro[3] + macro[1] * macro[1] + macro[2] * macro[2];
      return LC::weight(i) * macro[0] * (1.0 + (3.0 * ev + 1.5 * (3.0 * ev * ev - vv) / ltcVel) / ltcVel);
    }

    template<typename LC1 = LC>
    typename std::enable_if<::std::is_base_of<Lattice<3>, LC1>::value, void>::type  operator() (
      double *pdf, 
      const double *macro,
      const double &ltcVel = 1.0
    ) const {
      const double vv = macro[3] * macro[3] + macro[1] * macro[1] + macro[2] * macro[2];
      for(size_t i = 0; i < LC::size(); i++) {
        const double ev = macro[1] * (*(LC::direction(i))) + macro[2] * (*(LC::direction(i) + 1))
          + macro[3] * (*(LC::direction(i) + 2));
        pdf[i] = LC::weight(i)*macro[0] * (1.0 + (3.0 * ev + 1.5 * (3.0 * ev * ev - vv) / ltcVel) / ltcVel);
      }
    }
  };

  template<const size_t n>
  class TransMatrix
  {
  public:
    constexpr static inline size_t size() { return n; }
    static const double relaxation[n];
    static const double direct[n * n];
    static const double inverse[n * n];
  };

//	template<>
//  const double TransMatrix<9>::direct[9*9] =
//  {
//    1, 1, 1, 1, 1,1, 1, 1, 1, 
//    -4,-1,-1,-1,-1,2, 2, 2, 2,
//    4,-2,-2,-2,-2,1, 1, 1, 1,
//    0, 1, 0,-1, 0,1,-1,-1, 1,
//    0,-2, 0, 2, 0,1,-1,-1, 1,
//    0, 0, 1, 0,-1,1, 1,-1,-1,
//    0, 0,-2, 0, 2,1, 1,-1,-1,
//    0, 1,-1, 1,-1,0, 0, 0, 0,
//    0, 0, 0, 0, 0,1,-1, 1,-1
//  };
//
//	template<>
//  const double TransMatrix<9>::inverse[9*9] =
//  {
//    1.0/9, -1.0/9,   1.0/ 9,			0,			 0,      0,			  0,      0,      0,
//    1.0/9, -1.0/36, -1.0/18,  1.0/6,  -1.0/6,      0,			  0,  1.0/4,      0,
//    1.0/9, -1.0/36, -1.0/18,			0,			 0,  1.0/6, -1.0/ 6, -1.0/4,      0,
//    1.0/9, -1.0/36, -1.0/18, -1.0/6,   1.0/6,			 0,       0,  1.0/4,      0,
//    1.0/9, -1.0/36, -1.0/18,			0,			 0, -1.0/6,  1.0/ 6, -1.0/4,      0,
//    1.0/9,  1.0/18,  1.0/36,  1.0/6,  1.0/12,  1.0/6,  1.0/12,      0,  1.0/4,
//    1.0/9,  1.0/18,  1.0/36, -1.0/6, -1.0/12,  1.0/6,  1.0/12,      0, -1.0/4,
//    1.0/9,  1.0/18,  1.0/36, -1.0/6, -1.0/12, -1.0/6, -1.0/12,      0,  1.0/4,
//    1.0/9,  1.0/18,  1.0/36,  1.0/6,  1.0/12, -1.0/6, -1.0/12,      0, -1.0/4 
//  };
//
//	template<>
//  const double TransMatrix<9>::relaxation[9]=
//  {0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0};
//
//	template<>
//  const double TransMatrix<15>::direct[15*15] = {	
//    1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
//   -2, -1, -1, -1, -1, -1, -1,  1,  1,  1,  1,  1,  1,  1,  1,
//   16, -4, -4, -4, -4, -4, -4,  1,  1,  1,  1,  1,  1,  1,  1,
//    0,  1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1,
//    0, -4,  4,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1,
//    0,  0,  0,  1, -1,  0,  0,  1,  1, -1, -1,  1,  1, -1, -1,
//    0,  0,  0, -4,  4,  0,  0,  1,  1, -1, -1,  1,  1, -1, -1,
//    0,  0,  0,  0,  0,  1, -1,  1,  1,  1,  1, -1, -1, -1, -1,
//    0,  0,  0,  0,  0, -4,  4,  1,  1,  1,  1, -1, -1, -1, -1,
//    0,  2,  2, -1, -1, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0,
//    0,  0,  0,  1,  1, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0,
//    0,  0,  0,  0,  0,  0,  0,  1, -1, -1,  1,  1, -1, -1,  1,
//    0,  0,  0,  0,  0,  0,  0,  1,  1, -1, -1, -1, -1,  1,  1,
//    0,  0,  0,  0,  0,  0,  0,  1, -1,  1, -1, -1,  1, -1,  1,
//    0,  0,  0,  0,  0,  0,  0,  1, -1, -1,  1, -1,  1,  1, -1
//  }; 
//
//	template<>
//  const double TransMatrix<15>::inverse[15*15] = {
//    1./15., -1./ 9.,  2./ 45.,      0.,      0.,      0.,      0.,      0.,      0.,      0.,     0.,     0.,     0.,     0.,     0.,
//    1./15., -1./18., -1./ 90.,  1./10., -1./10.,      0.,      0.,      0.,      0.,  1./ 6.,     0.,     0.,     0.,     0.,     0.,
//    1./15., -1./18., -1./ 90., -1./10.,  1./10.,      0.,      0.,      0.,      0.,  1./ 6.,     0.,     0.,     0.,     0.,     0.,
//    1./15., -1./18., -1./ 90.,      0.,      0.,  1./10., -1./10.,      0.,      0., -1./12.,  1./4.,     0.,     0.,     0.,     0.,
//    1./15., -1./18., -1./ 90.,      0.,      0., -1./10.,  1./10.,      0.,      0., -1./12.,  1./4.,     0.,     0.,     0.,     0.,
//    1./15., -1./18., -1./ 90.,      0.,      0.,      0.,      0.,  1./10., -1./10., -1./12., -1./4.,     0.,     0.,     0.,     0.,
//    1./15., -1./18., -1./ 90.,      0.,      0.,      0.,      0., -1./10.,  1./10., -1./12., -1./4.,     0.,     0.,     0.,     0.,
//    1./15.,  1./18.,  1./360.,  1./10.,  1./40.,  1./10.,  1./40.,  1./10.,  1./40.,      0.,     0.,  1./8.,  1./8.,  1./8.,  1./8.,
//    1./15.,  1./18.,  1./360., -1./10., -1./40.,  1./10.,  1./40.,  1./10.,  1./40.,      0.,     0., -1./8.,  1./8., -1./8., -1./8.,
//    1./15.,  1./18.,  1./360.,  1./10.,  1./40., -1./10., -1./40.,  1./10.,  1./40.,      0.,     0., -1./8., -1./8.,  1./8., -1./8.,
//    1./15.,  1./18.,  1./360., -1./10., -1./40., -1./10., -1./40.,  1./10.,  1./40.,      0.,     0.,  1./8., -1./8., -1./8.,  1./8.,
//    1./15.,  1./18.,  1./360.,  1./10.,  1./40.,  1./10.,  1./40., -1./10., -1./40.,      0.,     0.,  1./8., -1./8., -1./8., -1./8.,
//    1./15.,  1./18.,  1./360., -1./10., -1./40.,  1./10.,  1./40., -1./10., -1./40.,      0.,     0., -1./8., -1./8.,  1./8.,  1./8.,
//    1./15.,  1./18.,  1./360.,  1./10.,  1./40., -1./10., -1./40., -1./10., -1./40.,      0.,     0., -1./8.,  1./8., -1./8.,  1./8.,
//    1./15.,  1./18.,  1./360., -1./10., -1./40., -1./10., -1./40., -1./10., -1./40.,      0.,     0.,  1./8.,  1./8.,  1./8., -1./8.
//  };
//
//	template<>
//  const double TransMatrix<15>::relaxation[15]= {
//    0.0,  // 0
//    0.5,  // 1
//    1.0,  // 2
//    0.0,  // 3
//    1.0,  // 4
//    0.0,  // 5
//    1.0,  // 6
//    0.0,  // 7
//    1.0,  // 8
//    1.0,  // 9
//    1.0,  //10
//    1.0,  //11
//    1.0,  //12
//    1.0,  //13
//    1.0   //14
//  };
//
//	template<>
//	const double TransMatrix<19>::direct[19 * 19] = {
//  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
//-30,-11,-11,-11,-11,-11,-11,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  8,
// 12, -4, -4, -4, -4, -4, -4,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
//  0,  1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0,
//  0, -4,  4,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0,
//  0,  0,  0,  1, -1,  0,  0,  1,  1, -1, -1,  0,  0,  0,  0,  1, -1,  1, -1,
//  0,  0,  0, -4,  4,  0,  0,  1,  1, -1, -1,  0,  0,  0,  0,  1, -1,  1, -1,
//  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  1,  1, -1, -1,  1,  1, -1, -1,
//  0,  0,  0,  0,  0, -4,  4,  0,  0,  0,  0,  1,  1, -1, -1,  1,  1, -1, -1,
//  0,  2,  2, -1, -1, -1, -1,  1,  1,  1,  1,  1,  1,  1,  1, -2, -2, -2, -2,
//  0, -4, -4,  2,  2,  2,  2,  1,  1,  1,  1,  1,  1,  1,  1, -2, -2, -2, -2,
//  0,  0,  0,  1,  1, -1, -1,  1,  1,  1,  1, -1, -1, -1, -1,  0,  0,  0,  0,
//  0,  0,  0, -2, -2,  2,  2,  1,  1,  1,  1, -1, -1, -1, -1,  0,  0,  0,  0,
//  0,  0,  0,  0,  0,  0,  0,  1, -1, -1,  1,  0,  0,  0,  0,  0,  0,  0,  0,
//  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1, -1, -1,  1,
//  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1, -1, -1,  1,  0,  0,  0,
//  0,  0,  0,  0,  0,  0,  0,  1, -1,  1, -1, -1,  1, -1,  1,  0,  0,  0,  0,
//  0,  0,  0,  0,  0,  0,  0, -1, -1,  1,  1,  0,  0,  0,  0,  1, -1,  1, -1,
//  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1, -1, -1, -1, -1,  1,  1
//	};
//
//	template<>
//	const double TransMatrix<19>::inverse[19 * 19] = {
//1./19.,   -5./399.,  1./21.,      0.,      0.,      0.,      0.,      0.,      0.,      0.,      0.,      0.,      0.,      0.,     0.,  0.,     0.,     0.,       0.,
//1./19., -11./2394., -1./63.,  1./10., -1./10.,      0.,      0.,      0.,      0.,  1./18., -1./18.,      0.,      0.,      0.,     0.,  0.,     0.,     0.,       0.,
//1./19., -11./2394., -1./63., -1./10.,  1./10.,      0.,      0.,      0.,      0.,  1./18., -1./18.,      0.,      0.,      0.,     0.,  0.,     0.,     0.,       0.,
//1./19., -11./2394., -1./63.,      0.,      0.,  1./10., -1./10.,      0.,      0., -1./36.,  1./36.,  1./12., -1./12.,      0.,     0.,  0.,     0.,     0.,       0.,
//1./19., -11./2394., -1./63.,      0.,      0., -1./10.,  1./10.,      0.,      0., -1./36.,  1./36.,  1./12., -1./12.,      0.,     0.,  0.,     0.,     0.,       0.,
//1./19., -11./2394., -1./63.,      0.,      0.,      0.,      0.,  1./10., -1./10., -1./36.,  1./36., -1./12.,  1./12.,      0.,     0.,  0.,     0.,     0.,       0.,
//1./19., -11./2394., -1./63.,      0.,      0.,      0.,      0., -1./10.,  1./10., -1./36.,  1./36., -1./12.,  1./12.,      0.,     0.,  0.,     0.,     0.,       0.,
//1./19.,   4./1197., 1./252.,  1./10.,  1./40.,  1./10.,  1./40.,      0.,      0.,  1./36.,  1./72.,  1./12.,  1./24.,   1./4.,     0.,  0.,  1./8., -1./8.,       0.,
//1./19.,   4./1197., 1./252., -1./10., -1./40.,  1./10.,  1./40.,      0.,      0.,  1./36.,  1./72.,  1./12.,  1./24.,  -1./4.,     0.,  0., -1./8., -1./8.,       0.,
//1./19.,   4./1197., 1./252.,  1./10.,  1./40., -1./10., -1./40.,      0.,      0.,  1./36.,  1./72.,  1./12.,  1./24.,  -1./4.,     0.,  0.,  1./8.,  1./8.,       0.,
//1./19.,   4./1197., 1./252., -1./10., -1./40., -1./10., -1./40.,      0.,      0.,  1./36.,  1./72.,  1./12.,  1./24.,   1./4.,     0.,  0., -1./8.,  1./8.,       0.,
//1./19.,   4./1197., 1./252.,      0.,      0.,  1./10.,  1./40.,   1./2.,   1./8., -1./18., -1./36.,      0.,      0.,      0.,  1./4., -1.,     0.,  1./8.,    3./8.,
//1./19.,   4./1197., 1./252.,      0.,      0., -1./10., -1./40., -3./10., -3./40.,   1./9.,  1./18.,  -1./6., -1./12.,      0., -1./4.,  1.,     0., -1./8.,   -1./8.,
//1./19.,   4./1197., 1./252.,   1./5.,  1./20., -1./10., -1./40., -1./2.,   -1./8.,   1./9.,  1./18.,  -1./6., -1./12.,      0., -1./4.,  1., -1./4., -1./8.,   -3./8.,
//1./19.,   4./1197., 1./252.,  -1./5., -1./20.,  1./10.,  1./40.,  3./10.,  3./40., -1./18., -1./36.,      0.,      0.,      0.,  1./4., -1.,  1./4.,  1./8.,    1./8.,
//1./19.,   4./1197., 1./252.,      0.,      0.,  1./10.,  1./40.,  1./10.,  1./40., -1./18., -1./36.,      0.,      0.,      0.,  1./4.,  0.,     0.,  1./8.,   -1./8.,
//1./19.,   4./1197., 1./252.,      0.,      0., -1./10., -1./40.,  1./10.,  1./40., -1./18., -1./36.,      0.,      0.,      0., -1./4.,  0.,     0., -1./8.,   -1./8.,
//1./19.,   4./1197., 1./252.,      0.,      0.,  1./10.,  1./40., -1./10., -1./40., -1./18., -1./36.,      0.,      0.,      0., -1./4.,  0.,     0.,  1./8.,    1./8.,
//1./19.,   4./1197., 1./252.,      0.,      0., -1./10., -1./40., -1./10., -1./40., -1./18., -1./36.,      0.,      0.,      0.,  1./4.,  0.,     0., -1./8.,    1./8.
//	};
//
//	template<>
//	const double TransMatrix<19>::relaxation[19] = {
//	
//	};

  /** 标签类，用于模板类和模板函数识别格子类型
  */
  template<
    typename LC, 
    typename EQ
    //template<typename LC1> typename EQ
  >
  class BGK {
  public:
    explicit BGK(const LC &lcIn, const double &tauIn) : lc(lcIn), tau(tauIn) {}
    void collide(
      double *pdf, /**< 格点位置处的分布函数：D2Q9对应九个方向上的分布函数 */
      double *macroVar
    ) const {
      cal_macro_var(macroVar, pdf, lc);
      for (int i = 0; i < LC::size(); i++)
        pdf[i] = pdf[i] - 1 / tau * (pdf[i] - equilibrium(macroVar, i, lc.velocity()));
    }
  private:
    const EQ equilibrium;
    const LC &lc;
    const double tau;
  };

  /** 标签类，用于模板类和模板函数识别格子类型
  */
  template< 
    typename LC, 
    typename TM, 
    typename EQ
    //template<typename LC1> typename EQ
  >
  class MRT {
  public:
    explicit MRT(
      const LC &lcIn, 
      const double &tau,
      const double *sIn = TM::relaxation
    ) : lc(lcIn), equilibrium() {
      for(size_t i = 0; i < LC::size(); i++) {
        s[i] = sIn[i] / tau;
      }
    }

    template<typename LC1 = LC, typename TM1 = TM>
    typename std::enable_if<LC1::size() == TM::size(), void>::type 
    pdf_to_moment(
      double *moment, 
      const double* pdf
    ) const {
      const double* row = TM::direct;
      for (size_t iM = 0; iM < LC::size(); iM++)
      {
        double& mom = moment[iM];
        mom = 0.0;
        for (size_t iVel = 0; iVel < LC::size(); iVel++)
        {
          mom += row[iVel] * pdf[iVel];
        }
        row += LC::size();
      }
    }

    template<typename LC1 = LC, typename TM1 = TM>
    typename std::enable_if<LC1::size() == TM::size(), void>::type 
    MomentToPdf(
      double *pdf, 
      const double *moment
    ) const {
      const double* row = TM::inverse;
      for (int iVel = 0; iVel < LC::size(); iVel++) {
        pdf[iVel] = 0.0;
        for (int iM = 0; iM < LC::size(); iM++)
        {
          pdf[iVel] += row[iM] * moment[iM];
        }
        row += LC::size();
      }
    }

    template<typename LC1 = LC, typename TM1 = TM>
    typename std::enable_if<LC1::size() == TM::size(), void>::type 
    collide(
      double* pdf, /**< 格点位置处的分布函数：D2Q9对应九个方向上的分布函数 */
      double* macroVar
    ) const {
      double moment[LC::size()];
      double momentEq[LC::size()];
      pdf_to_moment(moment, pdf);

      cal_macro_var(macroVar, pdf, lc);
      double pdfEq[LC::size()];
      equilibrium(pdfEq, macroVar, lc.velocity());
      pdf_to_moment(momentEq, pdfEq);

      for (size_t i = 0; i < LC::size(); i++)
      {
        moment[i] = moment[i] - s[i] * (moment[i] - momentEq[i]);
      }

      MomentToPdf(pdf, moment);
    }

  private:
    const EQ equilibrium;
    const LC &lc;
    double s[LC::size()];
  };

  extern COMMONAPI const double relaxationStandardD2Q9[D2Q9::size()];
  //extern COMMONAPI const double relaxationStandardD3Q15[D2Q9::size()];

  extern COMMONAPI const double transMatrixStandardD2Q9[D2Q9::size()*D2Q9::size()];
  extern COMMONAPI const double invTransMatrixStandardD2Q9[D2Q9::size()*D2Q9::size()];
}

#endif //!COLLISION_H
