///////////////////////////////////////////////////////////////
// author : Wang,Junquan. email: as2697920806@gmail.com
///////////////////////////////////////////////////////////////
// participant 1 : Chen,Songze. email: 
///////////////////////////////////////////////////////////////
// tester 1 : Wang,Junquan
///////////////////////////////////////////////////////////////

#pragma once
#ifndef INITFIELD_H
#define INITFIELD_H

#include "CommonAPI.h"
#include "Range.h"
#include "Collision.h"
#include "UniformGrid.h"
#include "Function.h"
#include <iostream>

namespace nnwlbm {
	/** \brief 初始化流场
		*输入MaskData、initCond、Grid
		*首先判断初始化条件的参数和MaskData中的变量参数是否相同
		*从MaskData得到对应的nodeRange，得到Range内的location
		*location通过grid转化到Grid内的实际位置position
		*通过初始条件initCond，在对应的position对Maskdata中的数据进行初始化
		* \param md MaskData [in,out] 对象，存储相关数据
		* \param initCond [in] 初始化条件
		* \param grid [in] 实际的计算区域
		* \return 返回bool变量
		*/
	template<const size_t dim> //extern "C"
	bool ini_field(
		MaskData<double, dim> &md, /**< [in,out]存储需要处理的数据*/
		const Function &initCond, /**< [in]初始化条件*/
		const UniformGrid<dim> &grid /**< [in]实际的计算区域*/
	)
	{
		if (initCond.GetOutputDim() != md.GetVarNum()) {
			std::cout << "Dimension of the initial condition does not match. \n\t initCond->GetOutputDim() = "
				<< initCond.GetOutputDim() << "\n\t md.GetVarNum() = " << md.GetVarNum() << std::endl;
			return false;
		}

		nnwlbm::RangeIndex<dim> rIndex(md.GetRange(), nnwlbm::RangeIndex<dim>::inputNodeRange);
		auto nodeRange = rIndex.GetNodeRange();
		int loc[dim];
		double pos[dim];
		for (int carry = nodeRange.IniCarry(loc); carry == 0; carry = nodeRange.AddOne(loc)) {
			size_t nodeId = rIndex.CalNodeID(loc);
			grid.GetPos(pos, loc);
			if (!md.IsValid(nodeId)) continue;

			if (!initCond.cal_var(md.GetVar(nodeId), pos)) {
				std::cout << "************* Error : No initial condition is applied at node " << nodeId << ". *************\n Position: (";
				for (int i = 0; i < dim; i++)
					std::cout << pos[i] << ' ';
				std::cout << ')' << std::endl;
				getchar();
				exit(1);
			}
		}
		return true;
	}

	template<const size_t dim> //extern "C"
	bool ini_field(
		nnwlbm::LinkedMaskData<double, dim> &md, /**< [in,out]存储需要处理的数据*/
		const nnwlbm::Function &initCond, /**< [in]初始化条件*/
		const UniformGrid<dim> &grid /**< [in]实际的计算区域*/
	)
	{
		if (initCond.GetOutputDim() != md.GetVarNum()) {
			std::cout << "Dimension of the initial condition does not match. \n\t initCond->GetOutputDim() = "
				<< initCond.GetOutputDim() << "\n\t md.GetVarNum() = " << md.GetVarNum() << std::endl;
			return false;
		}

		nnwlbm::RangeIndex<dim> rIndex(md.GetRange(), nnwlbm::RangeIndex<dim>::inputNodeRange);
		auto nodeRange = rIndex.GetNodeRange();
		int loc[dim];
		double pos[dim];
		for (int carry = nodeRange.IniCarry(loc); carry == 0; carry = nodeRange.AddOne(loc)) {
			size_t nodeId = rIndex.CalNodeID(loc);
			grid.GetPos(pos, loc);
			if (!md.IsValid(nodeId)) continue;

			if (!initCond.cal_var(md.GetVar(nodeId), pos)) {
				std::cout << "************* Error : No initial condition is applied at node " << nodeId << ". *************\n Position: (";
				for (int i = 0; i < dim; i++)
					std::cout << pos[i] << ' ';
				std::cout << ')' << std::endl;
				getchar();
				exit(1);
			}
		}
		return true;
	}

  template<typename LC> 
  bool CalFieldPdf(
    nnwlbm::LinkedMaskData<double, LC::GetDim()> &pdfData, /**< [in,out]存储需要处理的数据*/
    const nnwlbm::LinkedMaskData<double, LC::GetDim()> &macroData, /**< [in,out]存储需要处理的数据*/
    const Range<LC::GetDim()> &loopRange /**< [in] 计算分布函数的范围 */
  ) {
    const size_t dim = LC::GetDim();
    nnwlbm::Equilibrium<LC> eqStat;
    int loc[dim];
    nnwlbm::RangeIndex<dim> rIndex(macroData.GetRange(), nnwlbm::RangeIndex<dim>::inputNodeRange);
    for (int carry = loopRange.IniCarry(loc); carry == 0; carry = loopRange.AddOne(loc))
    {
      size_t nodeId = rIndex.CalNodeID(loc);
      if(macroData.IsValid(nodeId)) {
        double* pdf = pdfData.GetVar(nodeId);
        for(size_t iVel = 0; iVel < LC::size(); iVel++){
          pdf[iVel] = eqStat(macroData.GetVar(nodeId), iVel);
        }
      }
    }
    return true;
  }

  template<typename LC> 
  bool CalFieldPdf_save(
    nnwlbm::LinkedMaskData<double, LC::GetDim()> &pdfData, /**< [in,out]存储需要处理的数据*/
    const nnwlbm::LinkedMaskData<double, LC::GetDim()> &macroData, /**< [in,out]存储需要处理的数据*/
    const Range<LC::GetDim()> &loopRange /**< [in] 计算分布函数的范围 */
  ) {
    const size_t dim = LC::GetDim();
    nnwlbm::Equilibrium<LC> eqStat;
    Range<dim> rg = Range<dim>::CalIntersection(pdfData.GetRange(), macroData.GetRange());
    rg = Range<dim>::CalIntersection(rg, loopRange);
    int loc[dim];
    nnwlbm::RangeIndex<dim> rIndex(macroData.GetRange(), nnwlbm::RangeIndex<dim>::inputNodeRange);
    for (int carry = loopRange.IniCarry(loc); carry == 0; carry = loopRange.AddOne(loc))
    {
      size_t nodeId = rIndex.CalNodeID(loc);
      if(macroData.IsValid(nodeId)) {
        double* pdf = pdfData.GetVar(nodeId);
        for(size_t iVel = 0; iVel < LC::size(); iVel++){
          pdf[iVel] = eqStat(macroData.GetVar(nodeId), iVel);
        }
      }
    }
  }

}

#endif // INIFIELD_H
