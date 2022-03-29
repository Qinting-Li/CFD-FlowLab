////////////////////////////////////////////////////////////////////////////////////
// author: Chen, Songze. email:jacksongze@hust.edu.cn
////////////////////////////////////////////////////////////////////////////////////
// participant 1: Yang, Zhou. email:yz1249911657@gmail.com
////////////////////////////////////////////////////////////////////////////////////
// tester 1: Yang, Zhou. email: yz1249911657@gmail.com
////////////////////////////////////////////////////////////////////////////////////

// Lattice.h
#pragma once
#ifndef LATTICE_H
#define LATTICE_H
#include "CommonAPI.h"
#include "DimClass.h"
#include <iostream>

namespace nnwlbm {

  template<const size_t dim>
  class Lattice : public Dim<dim> { };

	//============================================求矩函数（微观->宏观）=====================================================
	/** 计算宏观量rho、u、v
	 * \param macroVar [out]储存宏观量的数组，按照rho、u、v的顺序存储
	 * \param pdf [in]储存平衡态分布函数值的数组
	 * \param lc [in]输入一个二维格子构型的对象
   * \return bool型，表示密度是否大于零
	 */
	template<typename T>
  typename std::enable_if<::std::is_base_of<Lattice<2>, T>::value, bool>::type 
    cal_macro_var(
      double *macroVar, 
      const double *pdf, 
      const T &lc
    ) 
  {
		double &rho = macroVar[0];
		double &u = macroVar[1];
		double &v = macroVar[2];
    rho = u = v = 0;
		for (size_t i = 0; i < lc.size(); i++)
		{
			const int *pointer;
			pointer = lc.direction(i);
			rho += pdf[i];
			u += *(pointer) * pdf[i];
			v += *(pointer + 1) * pdf[i];
		}
		u = u * lc.velocity() / rho;
		v = v * lc.velocity() / rho;
		return rho > 0.0;
	}

	//============================================求矩函数（微观->宏观）=====================================================
	/** \brief 计算宏观量rho、u、v
	 * \param macroVar [out]储存宏观量的数组，按照rho、u、v、w的顺序存储
	 * \param pdf [in]储存平衡态分布函数值的数组
	 * \param lc [in]表示3维的格子构型的对象
   * \return bool型，表示密度是否大于零
	*/
	template<typename T>
  typename std::enable_if<::std::is_base_of<Lattice<3>, T>::value, bool>::type
    cal_macro_var(
      double *macroVar, 
      const double *pdf, 
      const T &lc
    )
	{
		double &rho = macroVar[0];
		double &u = macroVar[1];
		double &v = macroVar[2];
		double &w = macroVar[3];
    rho = u = v = w = 0.0;
		for (int i = 0; i < lc.size(); i++)
		{
			const int *pointer;
			pointer = lc.direction(i);
			rho += pdf[i];
			u += *(pointer) * pdf[i];
			v += *(pointer + 1) * pdf[i];
			w += *(pointer + 2) * pdf[i];
		}
		u = u * lc.velocity() / rho;
		v = v * lc.velocity() / rho;
		w = w * lc.velocity() / rho;
		return rho > 0.0;
	}

  template<typename LC>
  void pdf_to_moment(double *moment, const double *pdf) {
    const double *row = LC::get_transMatrix();
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

  template<typename LC>
  void MomentToPdf(double *pdf, const double *moment) {
    const double *row = LC::get_invTransMatrix();
    for (int iVel = 0; iVel < LC::size(); iVel++) {
      pdf[iVel] = 0.0;
      for (int iM = 0; iM < LC::size(); iM++)
      {
        pdf[iVel] += row[iM] * moment[iM];
      }
      row += LC::size();
    }
  }

};

#endif
// end of Lattice.h