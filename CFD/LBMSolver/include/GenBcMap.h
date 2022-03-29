/////////////////////////////////////////////////
// author: Chen, Songze.  email: jacksongze@hust.edu.cn
////////////////////////////////////////////////
// collaborator 1: He, Shuangwu.  email: shuangwuh@gmail.com
///////////////////////////////////////////////
// test 1:He, Shuangwu.  email: shuangwuh@gmail.com
///////////////////////////////////////////////

#pragma once
#ifndef GENBCMAP_H
#define GENBCMAP_H
#include "D2Q9.h"
#include "BoundaryCondition.h"
#include "UniformGrid.h"
#include "Range.h"
#include "GenColorMap.h"
#include "Geometry.h"
#include "BoundaryFactory.h"

namespace nnwlbm {
/**
* \brief BcData类的数据集为map容器，存放一个ID和ID对应的边界几何体的全部信息。
* 不同的格点可计算出不同的ID值，但不同格点可能具有相同的边界信息。
* 操作集包括构造函数，判断ID是否为map成员，添加map成员，删除map成员，返回ID对应的Geo成员
*/
template<typename BC, typename DC>
class BcData : public nnwlbm::DataClass
{
private:
	MemGen<DC, 64> memGen;
protected:
	/**整型存放ID值，pair型存放ID对应的边界Geo和Geo可能需要的数据*/
	std::map<size_t, std::pair<BC*, DC*>> bcAndVar;
public:
	/**
	* \brief构造函数
	* \param varNum: [in]memGen内存池开内存的数目 
	*/
	explicit BcData(const size_t varNum) : memGen(varNum) {}//explicit防止隐式类型转换
	/**析构函数，回收C* */
	~BcData() {
		for (auto &elem : bcAndVar)
		{
			memGen.recycle(elem.second.second);
		}
	}
	/**
	* \brief 判断ID是否为map容器中的成员。
	* \param id: [in]ID值，不同格点对应不同的ID值。
	* \return 返回bool值。
	* \return 1.true ID包含在map中。
	* \return 2.false ID未包含在map中。
	*/
	bool IsValid(
		const size_t &id
	) const{ 
		return bcAndVar.find(id) != bcAndVar.end(); 
	}
	/**
	* \brief 向map容器中添加几何体信息。
	* \param id:[in]ID值，不同格点对应不同的ID值。
	* \param bc:[in]Geometry*类，默认为空指针，指具体的几何体信息。
	* 如果ID已存在，则将bc推入，如果不存在，则新开辟一个ID并推入bc。
	*/
	void SetBC(
		const size_t &id, 
		Geometry *bc = nullptr
	) {
		if (IsValid(id)) bcAndVar[id] = std::make_pair(bc, bcAndVar.at(id).second);
		else bcAndVar.emplace(id, std::make_pair(bc, memGen.GetNew()));
	}
	/**
	* \brief 删除map容器中的ID值及其对应的几何体。
	* \param id:[in]ID值，不同的格点对应不同的ID值。
	* 如ID不属于map，则不操作，否则回收ID值对应的DC*所占的内存，并擦除ID值。
    */
	void del_bc(
		const int &id
	) {
		if (!IsValid(id)) return;
		memGen.recycle(bcAndVar[id].second);
		bcAndVar.erase(id);
	}
	/**
	* \brief 返回ID值对应Geo类相应的物理量的值。
	* \param id:[in]ID值，不同格点对应不同的ID值。
	*/
	DC* GetVar(
		const int& id
	) { 
		return bcAndVar.at(id).second; 
	}
	/**
	* \brief 返回ID值对应Geo类相应的物理量的值(常).
	* \param id:[in]ID值，不同格点对应不同的ID值。
	*/
	const DC *GetVar(
		const int &id
	) const { 
		return bcAndVar.at(id).second; 
	}
	/**
	* \brief 返回ID值对应的几何体
	* \param id:[in]ID值，不同格点对应不同的ID值。
	*/
	Geometry *GetBcFun(
		const int &id
	) { 
		return bcAndVar.at(id).first;
	}
	/**
	* \brief 返回ID值对应的几何体(常)
	* \param id:[in]ID值，不同格点对应不同的ID值。
	*/
	const Geometry *GetBcFun(
		const int &id
	) const { 
		return bcAndVar.at(id).first; 
	}
public:
	/** 返回类型名称. 通过name函数内部定义字符串实现一种类型一个名称 */
	virtual std::string name() const { return "BcData" + std::to_string(sizeof(DC)); } 
private:
	virtual bool OutputBinPrivate(
		std::ofstream &ostr,
		const std::string &filename,
		const std::string &path
	) const 
	{
		return false;
	}
};


	/**
	* \brief 函数gen_bcMap：根据color值设置边界条件。
	* \note 调用gen_bcMap后，所有的格点对应的边界条件(Geometry)信息均存在bcAndVar中。
	* \param bcData: [in][out]BcData类对象,开辟内存存放与每一个格点对应的几何体，每一个格点有唯一的ID值。
	* \param md：[in]MaskData类对象，提供color值及Range的范围。
	* \param grid: [in]UniformGrid类，笛卡尔直角网格。
	* \param bcBody: [in]容器vector<Geometry*>，包括计算区域内的所有几何体。
	* \param bcSurface: [in]容器vector<Geometry*>,当边界条件都相同时，BcSurface为空；否则BcSurface放置相应的边界条件。
	* \return bool型
	* \retval 1.true 函数调用成功。
	* \retval 2.false 函数调用失败。
	*/
	template<const size_t dim>
	bool gen_bcMap(
		BcData<Geometry, double>& bcData,/**< [in][out]开辟内存存放与每一个格点对应的几何体，每一个格点有唯一的ID值。 */
		MaskData<double, dim> &md,/**< [in]提供color值及Range的范围。*/
		const UniformGrid<dim> &grid,/**< [in]笛卡尔直角网格。*/
		std::vector<Geometry*> &bcBody,/**< [in]包括计算区域内的所有几何体。*/
		std::vector<Geometry*> &bcSurface/**< [in]补充的几何体*/
	)
	{
		int loc[dim], faceLoc[dim];

		RangeIndex<dim> rIndex(md.GetRange(), RangeIndex<dim>::inputCellRange);
		auto faceRange = rIndex.GetFaceRange();//从md中取得facerange(vector<Range<dim>>类型)

		for (int iDim = 0; iDim < dim; iDim++) {//维数循环,即facerange的个数,从而对每一个facerange操作
			Range<dim> loopRange = faceRange[iDim];
			loopRange.AddBuffer(-1);//Range中会多一层虚网格，去掉虚网格后为range的实际大小

			for (int carry = loopRange.IniCarry(faceLoc);
				carry == 0; carry = loopRange.AddOne(faceLoc))//loopRange循环
			{
				for (int i = 0; i < dim; i++)
					loc[i] = faceLoc[i];
				loc[iDim] = faceLoc[iDim] - 1;//loc在faceloc相应方向少1，相当于nodeRange

				size_t cellId = rIndex.CalCellID(faceLoc);
				size_t cellIdN = rIndex.CalCellID(loc);//减1后的faceloc

				if (md.IsValid(cellId) ^ md.IsValid(cellIdN))//按位异或运算，相同为0，相异为1.-1即为边界，face比node向前一位
				{
					int ibc = 0; 
          if(md.IsValid(cellId)) {
            ibc = ColorToBcBodyID(md.GetColor(cellIdN));
            md.SetColor(cellId, geo::boundary);
          }
          else {
            ibc = ColorToBcBodyID(md.GetColor(cellId));
            md.SetColor(cellIdN, geo::boundary);
          }
#ifdef _DEBUG
          if (ibc < 0) {  //ibc<0,即区域为流体，正确的程序应该不会发生
            std::cout << "bc error" << std::endl; 
#ifdef WIN32
            system("pause"); 
#endif
            exit(1); 
          }
#endif

					auto bcGeo = bcBody[ibc];
					if (bcGeo == nullptr) {
						std::cerr << "bcBody input is wrong" << std::endl;
						return false;
					}
					size_t faceId = rIndex.CalFaceID(faceLoc, iDim);
					bcData.SetBC(faceId, bcGeo);//将bcBody中相应的边界Geo插入BcData中，并给与唯一的faceID(包含方向)

					// set geometric property
					// bc point
					//处理边界处边界条件不同的情况
					double cellC[dim], cellN[dim], faceC[dim];
					grid.GetPos(cellC, loc);//物理空间
					grid.GetPos(cellN, faceLoc);//物理空间
					//if (!bcGeo->CalIntersection(faceC, cellC, cellN)) {//计算区域与周围Geo的交点，如果无法实现简单Geo的重载计算，则使用二分法
						//Geometry::CalIntersection(faceC, *bcGeo, cellC, cellN, dim);
					//}//faceC存放Geo交点的坐标
					for (unsigned ibc = 0; ibc < bcSurface.size(); ibc++) {
						auto bcGeo = bcSurface[ibc];
						if (bcGeo == nullptr) {
							std::cerr << "BcSurface input is wrong" << std::endl;
							return false;
						}
						if (bcGeo->CalRelativeLocation(faceC) == geo::inside) {
							bcData.SetBC(faceId, bcGeo);
							break;
						}
					}
					// bc normal direction
					//if(!bcGeo.cal_normal(current+2*varNum+dim, current+2*varNum)) {
					//	if(dim == 3) {
					//		if( !geo::d3::cal_normal_vector(current+2*varNum+dim, bcGeo, current+2*varNum, 0.5*minL) )
					//			cout << "cannot find the normal." <<endl;
					//	}
					//}
					//memcpy(current+2*varNum+dim, faceN+dim*dim*bcElem.first, dim*sizeof(double));
				}
			}
		}
		return true;
	}
}

namespace nnwlbm {
  typedef nnwlbm::Comparator<int, std::equal_to<int>> BcIdentifier;
  typedef nnwlbm::LinkedMaskData<std::shared_ptr<nnwlbm::SingleNodeBc<D2Q9>>, 2, BcIdentifier> 
    BcMapSingleNode2D;

  template<typename LC>
  class BcMapSingleNode 
    : public nnwlbm::MaskData<
    nnwlbm::SingleNodeBc<LC>*, 
    LC::GetDim(), 
    BcIdentifier
    > 
  {
  public:
    BcMapSingleNode(const size_t varN, const MaskRange<int, LC::GetDim()> &mask, const int i) 
      : nnwlbm::MaskData<
      nnwlbm::SingleNodeBc<LC>*,
      LC::GetDim(),
      BcIdentifier >(varN, mask, BcIdentifier(i)) {}

  };

	template<typename LC>
	bool gen_bcMap(
		BcMapSingleNode<LC> &bcMap,/**< [in][out]*/
		const UniformGrid<LC::GetDim()> &grid,/**< [in]笛卡尔直角网格。*/
    std::map<nnwlbm::Geometry*, nnwlbm::SingleNodeBc<LC>*> bcBodyMap,
    std::map<nnwlbm::Geometry*, nnwlbm::SingleNodeBc<LC>*> bcSurfaceMap
  )
	{
    const size_t dim = LC::GetDim();
		int loc[dim], locN[dim];
    double pos[dim], posN[dim];
		RangeIndex<dim> rIndex(bcMap.GetRange(), RangeIndex<dim>::inputNodeRange);
		auto nodeRange = rIndex.GetNodeRange();

		for (int carry = nodeRange.IniCarry(loc);
		  carry == 0; carry = nodeRange.AddOne(loc))
		{
      grid.GetPos(pos, loc);
			size_t nodeId = rIndex.CalNodeID(loc);
      auto directionalBc = bcMap.GetVar(nodeId);
      if(directionalBc == nullptr) continue;
      double bcPos[dim];
      for(size_t iVel = 0; iVel < LC::size(); iVel++)
      {
        auto dir = LC::direction(iVel);
        for(size_t iDim = 0; iDim < LC::GetDim(); iDim++) {
          locN[iDim] = loc[iDim] - dir[iDim]; 
        }
		    //int nodeIdN = rIndex.CalCellID(locN);
        grid.GetPos(posN, locN);
        for(auto& elem: bcBodyMap) {
          auto elemGeo = elem.first;
          auto elemBc = elem.second;
          if(geo::boundary <= elemGeo->CalRelativeLocation(posN)) {
            directionalBc[iVel] = dynamic_cast< nnwlbm::SingleNodeBc<LC> * > (elemBc);
            if (!elemGeo->CalIntersection(bcPos, pos, posN)) {
              Geometry::CalIntersection(bcPos, *elemGeo, pos, posN);
            }
            break;
          }
        }
        for(auto& elem :bcSurfaceMap) {
          auto elemGeo = elem.first;
          auto elemBc = elem.second;
          if(geo::boundary <= elemGeo->CalRelativeLocation(bcPos)) {
            directionalBc[iVel] = dynamic_cast< nnwlbm::SingleNodeBc<LC> * > (elemBc);
            break;
          }
        }
      }
    }

		return true;
	}

	template<typename LC>
	bool gen_bcMap(
    nnwlbm::LinkedMaskData<nnwlbm::SingleNodeBc<LC>*, LC::GetDim(), nnwlbm::BcIdentifier > &bcMap,/**< [in][out]*/
		const UniformGrid<LC::GetDim()> &grid,/**< [in]笛卡尔直角网格。*/
    std::map<nnwlbm::Geometry*, nnwlbm::SingleNodeBc<LC>*> bcBodyMap,
    std::map<nnwlbm::Geometry*, nnwlbm::SingleNodeBc<LC>*> bcSurfaceMap
  )
	{
    const size_t dim = LC::GetDim();
		int loc[dim], locN[dim];
    double pos[dim], posN[dim];
		RangeIndex<dim> rIndex(bcMap.GetRange(), RangeIndex<dim>::inputNodeRange);
		auto nodeRange = rIndex.GetNodeRange();

		for (int carry = nodeRange.IniCarry(loc);
		  carry == 0; carry = nodeRange.AddOne(loc))
		{
      grid.GetPos(pos, loc);
			size_t nodeId = rIndex.CalNodeID(loc);
      auto directionalBc = bcMap.GetVar(nodeId);
      if(directionalBc == nullptr) continue;
      double bcPos[dim];
      for(size_t iVel = 0; iVel < LC::size(); iVel++)
      {
        auto dir = LC::direction(iVel);
        for(size_t iDim = 0; iDim < LC::GetDim(); iDim++) {
          locN[iDim] = loc[iDim] - dir[iDim]; 
        }
		    //int nodeIdN = rIndex.CalCellID(locN);
        grid.GetPos(posN, locN);
        for(auto& elem: bcBodyMap) {
          auto elemGeo = elem.first;
          auto elemBc = elem.second;
          if(geo::boundary <= elemGeo->CalRelativeLocation(posN)) {
            directionalBc[iVel] = dynamic_cast< nnwlbm::SingleNodeBc<LC> * > (elemBc);
            if (!elemGeo->CalIntersection(bcPos, pos, posN)) {
              Geometry::CalIntersection(bcPos, *elemGeo, pos, posN);
            }
            break;
          }
        }
        for(auto& elem :bcSurfaceMap) {
          auto elemGeo = elem.first;
          auto elemBc = elem.second;
          if(geo::boundary <= elemGeo->CalRelativeLocation(bcPos)) {
            directionalBc[iVel] = dynamic_cast< nnwlbm::SingleNodeBc<LC> * > (elemBc);
            break;
          }
        }
      }
    }

		return true;
	}
}

namespace nnwlbm {

  template<typename LC>
  bool gen_bcMap(nnwlbm::LinkedMaskData<std::shared_ptr<nnwlbm::SimpleBC>, LC::GetDim(), nnwlbm::BcIdentifier > &bcMap,//< [in][out]
    const UniformGrid<LC::GetDim()> &grid,// [in]笛卡尔直角网格
    std::map<int,std::pair<nnwlbm::Geometry*, nnwlbm::BoundaryFactory<LC>*>> bcBodyMap,
    std::map<nnwlbm::Geometry*, nnwlbm::BoundaryFactory<LC>*> bcSurfaceMap,
    nnwlbm::LinkedMaskData<double, LC::GetDim()> &pdfData
  )
  {
    const size_t dim = LC::GetDim();
    int loc[dim], locN[dim];
    double pos[dim], posN[dim];
    RangeIndex<dim> rIndex(bcMap.GetRange(), RangeIndex<dim>::inputNodeRange);
    auto nodeRange = rIndex.GetNodeRange();

    for (int carry = nodeRange.IniCarry(loc); carry == 0; carry = nodeRange.AddOne(loc))
    {
      grid.GetPos(pos, loc);
      size_t nodeId = rIndex.CalNodeID(loc);
      auto directionalBc = bcMap.GetVar(nodeId);
      if (directionalBc == nullptr) continue;
      for (size_t iVel = 0; iVel < LC::size(); iVel++)
      {
        double bcPos[dim];
        auto dir = LC::direction(iVel);
        for (size_t iDim = 0; iDim < LC::GetDim(); iDim++) {
          locN[iDim] = loc[iDim] - dir[iDim];
        }
        //int nodeIdN = rIndex.CalCellID(locN);
        grid.GetPos(posN, locN);
        for (auto& elem : bcBodyMap) {
          auto elemGeo = elem.second.first;
          auto elemBc = elem.second.second;
          if (geo::boundary < elemGeo->CalRelativeLocation(posN)) {
            directionalBc[iVel] = elemBc->CreateBoundary(loc, iVel, pdfData);
            if (!elemGeo->CalIntersection(bcPos, pos, posN)) {
              Geometry::CalIntersection(bcPos, *elemGeo, pos, posN);
            }
            break;
          }
        }
        for (auto& elem : bcSurfaceMap) {
          auto elemGeo = elem.first;
          auto elemBc = elem.second;
          if (geo::boundary <= elemGeo->CalRelativeLocation(bcPos)) {
            directionalBc[iVel] = elemBc->CreateBoundary(loc, iVel, pdfData);
            //std::cout << loc[0] << "   " << loc[1] << std::endl;
            //bcPos[0] = -10.; bcPos[1] = -10.;
            break;
          }
        }
      }
    }
    return true;
  }
}

#endif // !BCMAP_H
