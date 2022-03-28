////////////////////////////////////////////////////////////////////////////////////
// author: Yang, Zhou. email:yz1249911657@gmail.com
////////////////////////////////////////////////////////////////////////////////////
// tester 1:: Yang, Zhou. email:yz1249911657@gmail.com
////////////////////////////////////////////////////////////////////////////////////
// tester 2: He, jiaming. email: 2457837217@qq.com
////////////////////////////////////////////////////////////////////////////////////

// D3Q15.h
#pragma once
#ifndef D3Q15_H
#define D3Q15_H

#include "CommonAPI.h"
#include "BinBase.h"
#include "Lattice.h"
#include <iostream>

// 当前代码都存在nnwlbm名字空间下，便于后期协同操作。
namespace nnwlbm
{
	class COMMONAPI D3Q15 :public nnwlbm::DataClass, public Lattice<3>
	{
	public:
		/**
			* \brief D3Q15类的构造函数：实现输入一个double类型的数据，将值赋给private中的格子速度c
			* \param cc [in] 外部用户给定的格子速度
			* 示例，D3Q15 dq(1.58),所以c=1.58
		 */
		explicit D3Q15(const double &cc) :c(cc) { };

		/**
			* \brief 参数输入检查函数
			* \param i [in] 离散速度空间编号
			* \return 0-14之内（包含0/14）返回true，其他返回false
		 */
		static inline bool check(const size_t &i)
		{
			if (i < 0 || i>14) { 
        std::cout << "error : index " << i << " is out of range!" << std::endl; 
        exit(1);
      }
			else return true;
		}

		/**
			* \brief 计算D3Q15类型对象的离散速度个数
			* \return sizeLc 15
		 */
		constexpr static inline size_t size() { return 15; }

		/**
			* \brief 计算D3Q15类型对象的格子速度
			* \return 外部用户给定的格子速度c
		 */
		inline const double &velocity() const
		{
			return c;
		}

		/**
			* \brief 计算D3Q15类型对象的迁移方向
			* 在debug状态下调用此函数会执行错误检查功能，即调用check函数来判断输入参数是否超出阈值
			* \param i [in] 离散速度空间编号
			* \return 返回一个int类型的指针p，(*p)表示x方向，*(p+1)表示y方向，*(p+2)表示z方向
		 */
		static inline const int *direction(const size_t &i)
		{
#ifdef _DEBUG
			check(i);
#endif // _DEBUG
			return(dir + 3 * i);
		}

		/**
			* \brief 计算D3Q15类型对象的权重系数
			* 在debug状态下调用此函数会执行错误检查功能，即调用check函数来判断输入参数是否超出阈值
			* \param i [in] 离散速度空间编号
			* \return 返回输入离散速度空间编号所对应的权系数
		 */
		static inline const double &weight(const size_t &i)
		{
#ifdef _DEBUG
			check(i);
#endif // _DEBUG
			return w[i];
		}

		/**
			* 计算D3Q15类型对象的碰撞后的离散速度空间编号
			* 在debug状态下调用此函数会执行错误检查功能，即调用check函数来判断输入参数是否超出阈值
			* \param i [in] 离散速度空间编号
			* \return 返回输入离散速度空间编号所对应的反向离散速度空间编号
		 */
		static inline const int InvertedDirection(const size_t &i)
		{
#ifdef _DEBUG
			check(i);
			#endif // _DEBUG
			return invDir[i];
		}

	private:
		/**
			* 格子速度c
			* 此参数由指定D3Q15对象的时候给定，随着对象的不同，c可以随之改变
		*/
		const double c;

		/**
			* 正向离散速度空间dir：
			* 此参数当对象类型指定为D3Q15时就已经确定
			* dir数组里面的数据按0-14的离散速度空间编号顺序存储
			* dir[45]={ 0,0,0,1,0,0,-1,0,0,0,1,0,0,-1,0,0,0,1,0,0,-1,1,1,1,-1,-1,-1,1,1,-1,-1,-1,1,1,-1,1,-1,1,-1,1,-1,-1,-1,1,1 }
		 */
		static const int dir[45];

		/**
			* 反向离散速度空间编号invDir：
			* 此参数当对象类型指定为D3Q15时就已经确定
			* invDir数组里面的数据按正向0-14的离散速度空间编号顺序存储其碰撞后的反向离散速度空间编号
			* invDir[15]={ 0,3,4,1,2,6,5,14,13,12,11,10,9,8,7 }
		*/
		static const int invDir[15];

		/**
			* 权重系数w：
			* 此参数当对象类型指定为D3Q15时就已经确定
			* w[15]={ 2.0 / 9.0,1.0 / 9.0,1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 72.0,1.0 / 72.0, 1.0 / 72.0, 1.0 / 72.0, 1.0 / 72.0, 1.0 / 72.0, 1.0 / 72.0, 1.0 / 72.0 }
		 */
		static const double w[15];

		/**
			* D3Q15模型的离散速度个数：
			* 此参数当对象类型指定为D3Q15时就已经确定
			* sizeLc=15
		 */
		static const int sizeLc;
	public:
		virtual std::string name() const;
	private:
		/**
			* \brief D3Q15类二进制输出
		*/
		virtual bool OutputBinPrivate
		(
			std::ofstream &ostr,
			const std::string &filename,
			const std::string &path
		) const;
	};

	//===================================计算平衡态分布函数(宏观->微观)==================================
	//																				D3Q15
	//																			/				\
	//																		 /					\
	//																D3Q15BGK			D3Q15MRT
	//=================================================================================================
	/**BGK模型时的计算平衡态分布函数的类，继承自D3Q15类。*/
	class COMMONAPI D3Q15BGK :public D3Q15
	{
	public:
		/**
			* \brief D3Q15BGK类的构造函数：实现输入一个double类型的数据，将值赋给D3Q15类private中的格子速度c
			*示例，D3Q15BGK dq(1.58),所以c=1.58
		 */
		explicit D3Q15BGK(const double &cc) :D3Q15(cc) { }

		/**
			* \brief 计算D3Q15 BGK模型的平衡态分布函数
			* \param macro [in] 储存宏观量的数组，按照rho、u、v、w顺序存储
			* \param i [in] 离散速度空间编号
			* \return bgk模型下的平衡态分布函数计算值
		*/
		double CalEqui(
      const double *macro, 
      int const i
    ) const
		{
			double bgkEq = 0.0;
			double ev;
			double vv;
			ev = macro[1] * (*(direction(i))) + macro[2] * (*(direction(i) + 1)) + macro[3] * (*(direction(i) + 2));
			vv = macro[3] * macro[3] + macro[1] * macro[1] + macro[2] * macro[2];
			bgkEq = weight(i) * macro[0] * (1.0 + (3.0 * ev + 1.5 * (3.0 * ev * ev - vv) / velocity()) / velocity());
			return bgkEq;
		}
	public:
		virtual std::string name() const { return "D3Q15BGK"; }
	};

	/** MRT模型时的求解平衡态分布函数的类，继承自D3Q15类。*/
	class COMMONAPI D3Q15MRT : public D3Q15
	{
	public:
		virtual std::string name() const { return "D3Q15MRT"; }
	public:
		/**
			* \brief D3Q15MRT类的构造函数：实现输入一个double类型的数据，将值赋给D3Q15类private中的格子速度c
			* 示例，D3Q15MRT dq(1.58),所以c=1.58
		 */
		explicit D3Q15MRT(const double &cc) :D3Q15(cc) { };

		/**
			* \brief 计算D3Q15格子构型下MRT模型的变换矩阵
			* 在debug状态下调用此函数会执行错误检查功能，即调用check函数来判断输入参数是否超出阈值
			* \param i [in] 离散速度空间编号
			* \return 返回一个int类型的指针p，（*p（i））表示拿到第i行的第一个值
		 */
		static inline const double *trans_matrix(const int &i)
		{
#ifdef _DEBUG
			check(i);
#endif // _DEBUG
			return (transMatrix + 15 * i);
		}

		/**
			* \brief 计算D3Q15 BGK模型的平衡态分布函数
			* \param macro [in] 储存宏观量的数组，按照rho、u、v、w顺序存储
			* \param i [in] 离散速度空间编号。
			* \return mrt模型下的平衡态分布函数计算值
		 */
		double CalEqui(
      const double *macro, 
      int const i
    ) const
		{
			D3Q15BGK a(velocity());
			double mrtEq = 0.0;
			for(int j=0;j<D3Q15::size();j++)
			{
				mrtEq = mrtEq + (*(trans_matrix(i) + j)) * a.CalEqui(macro, j);
			}
			return mrtEq;
		}
	private:
		/**
			* D3Q15模型的变换矩阵
			* 此参数当对象类型指定为D3Q15MRT时就已经确定
		 */
		static const double transMatrix[15 * 15];

		/**
		 * D3Q15格子构型MRT模型的逆矩阵；
		 * 此参数当对象类型指定为D3Q15时就已经确定。
		*/
		static const double invM[15 * 15];
	};
}
#endif // D3Q15_H


