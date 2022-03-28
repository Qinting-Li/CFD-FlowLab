/////////////////////////////////////////////////
// author: Chen, Songze.  email: jacksongze@hust.edu.cn
////////////////////////////////////////////////
// collaborator 1: He, Shuangwu.  email: shuangwuh@gmail.com
///////////////////////////////////////////////
// test 1:He, Shuangwu.  email: shuangwuh@gmail.com
///////////////////////////////////////////////

#pragma once
#ifndef GENCOLORMAP_H
#define GENCOLORMAP_H
#include "Range.h"
#include "UniformGrid.h"
#include "LBMGeometry.h"


namespace nnwlbm {

	/**
		*@brief 函数ColorToBcBodyID()实现从color值至bcBody中Geo位置序号的转换。
		*@param cId: MaskData类对象中存储的color值。
		*@return Geo在vector<Geometry*>容器中的位置。
		*/
	inline int ColorToBcBodyID(int cId) { return -cId - 1; }

	/**
		*@brief 函数bcBodyId_to_color()实现从bcBody中Geo的位置序号至color值的转换。
		*@param bcId: Geo在vector<Geometry*>容器bcBody中的位置。
		*@return MaskData类对象中存储的color值。
		*/
	inline int bcBodyId_to_color(int bcId) { return -bcId - 1; }

	/**
  * \brief 函数GenColorMap用于标记计算区域：根据Range中各点到Geo的相对距离，分别记录流场内，边界和流场外的标记值color。
  * \note 调用函数后，会在MashData对象中存入每一个坐标点的color值，从而确定该点的计算位置(流场内，边界和流场外)
  * \note 标记值color < 0 表示计算区域之外(color==-ibc-1表示不同的流场外区域块)，color=0表示靠近边界的计算区域。
  * \param md: [in][out] MaskData类，存放每一个Range节点的color值。
  * \param grid: [in]UniformGrid类，笛卡尔直角网格。
  * \param bcBody: [in]容器vector<Geometry*>，存放流场区域外的所有几何体。
  * \return bool型
  * \retval 1.true 函数调用成功。
  * \retval 2.false 函数调用失败。
  */
	template<const size_t dim>
	bool GenColorMap(
		MaskData<double, dim> &md,/**< [in][out] MaskData类，存放每一个Range节点的color值。*/
		const UniformGrid<dim> &grid,/**< [in]UniformGrid类，笛卡尔直角网格。*/
		std::vector<Geometry*> &bcBody/**< [in]容器vector<Geometry*>，存放流场区域外的所有几何体。*/
	) 
	{
		const size_t &bcNum = bcBody.size();	// the number of bcBody,整型常量的引用
		RangeIndex<dim> rIndex(md.GetRange(), RangeIndex<dim>::inputCellRange);//GetRange()返回一个range类
		
		Range<dim> cellRange = rIndex.GetCellRange();		//取出MaskData对象中的range，cellrange
		auto faceRange = rIndex.GetFaceRange();			//取出MaskData对象中的facerange，facerange的数目等于维数

		double cellC[dim];
		int loc[dim];
		// gen colorMap
		for (int carry = cellRange.IniCarry(loc);
			carry == 0; carry = cellRange.AddOne(loc))		//遍历cellrange网格空间
		{
			size_t cellId = rIndex.CalCellID(loc);			//计算loc位置(Range内)的cellID
			// color < 0 represents the external of computational domain.
			// color == -ibc-1 represents the solid part.
			int color = geo::inside;
			grid.GetPos(cellC, loc);						//cellC为loc相对坐标原点的坐标

			//按照每一个cellC的相对位置设定color值，扫描整个vector容器
			for (int ibc = 0; ibc < bcNum; ibc++) {
				auto bcGeo = bcBody[ibc];
				if (bcGeo->CalRelativeLocation(cellC) == geo::inside) {
					color = bcBodyId_to_color(ibc);			//如果cellC相对Geo的位置为inside，则对color操作，每一个Geo附近的color值不同
					break;									
				}
			}
			md.SetColor(cellId, color);					//修改MaskData中每一个点的值		
		}
		return true;
	}

	template<const size_t dim, class Compare = Color>
	void mark_boundary(LinkedMaskRange<int, dim> &lmr, Compare comp) {
		Range<dim> rg = lmr.GetRange();
		int loc[dim];
		std::set<size_t> pending;
		std::set<size_t> completed;
		for (int carry = rg.IniCarry(loc); carry == 0; carry = rg.AddOne(loc))
		{
			size_t id = CAL_ID(rg, loc);
			if (comp(lmr.GetColor(id))) {
				Range<dim> neighborRange;
				neighborRange.SetStart(loc);
				neighborRange.SetEnd(loc);
				neighborRange.AddBuffer(1);
				neighborRange = Range<dim>::CalIntersection(neighborRange, rg);
        int loc1[dim];
				for (int c1 = neighborRange.IniCarry(loc1); c1 == 0; c1 = neighborRange.AddOne(loc1))
				{
					size_t neighborId = CAL_ID(rg, loc1);
					if ( !(comp(lmr.GetColor(neighborId))) ) {
            lmr.SetColor(neighborId, geo::boundary);
          }
				}
			}
		}
	}

	//template<const size_t dim, template<typename C, const unsigned int dim1>class MR>
	template<const size_t dim,	typename BC>
	bool GenColorMap(
		LinkedMaskRange<int, dim> &lmr,  /**< [in][out] LinkedMaskRange类，存放每一个Range节点的color值。*/
		const UniformGrid<dim> &grid, /**< [in]UniformGrid类，笛卡尔直角网格。*/
		std::map<int,std::pair<nnwlbm::Geometry*, BC*>> &bcBodyMap /**< [in]容器vector<Geometry*>，存放流场区域外的所有几何体。*/
	) 
	{
		const size_t &bcNum = bcBodyMap.size();	// the number of bcBody,整型常量的引用
		RangeIndex<dim> rIndex(lmr.GetRange(), RangeIndex<dim>::inputCellRange);//GetRange()返回一个range类
		
		Range<dim> cellRange = rIndex.GetCellRange();		//取出MaskData对象中的range，cellrange
		auto faceRange = rIndex.GetFaceRange();			//取出MaskData对象中的facerange，facerange的数目等于维数

		double cellC[dim];
		int loc[dim];
		// gen colorMap
		for (int carry = cellRange.IniCarry(loc);carry == 0; carry = cellRange.AddOne(loc))		//遍历cellrange网格空间
		{
			size_t cellId = rIndex.CalCellID(loc);			//计算loc位置(Range内)的cellID
			// color < 0 represents the external of computational domain.
			// color == -ibc-1 represents the solid part.
			int color = geo::inside;
			grid.GetPos(cellC, loc);						//cellC为loc相对坐标原点的坐标

			//按照每一个cellC的相对位置设定color值，扫描整个vector容器
      for(auto &elem : bcBodyMap) {
        auto bcGeo = elem.second.first;
				if (bcGeo->CalRelativeLocation(cellC) == geo::inside) {
					color = bcBodyId_to_color(elem.first);			//如果cellC相对Geo的位置为inside，则对color操作，每一个Geo附近的color值不同
					break;									
				}
      }
			lmr.SetColor(cellId, color);					//修改MaskData中每一个点的值		
		}
    // 标记计算区域内临近边界的格点
    mark_boundary<dim, nnwlbm::Comparator<int, std::less<int>>>(lmr, nnwlbm::Comparator<int, std::less<int>>(0));

		return true;
	}

}

#endif // !GENCOLORMAP

