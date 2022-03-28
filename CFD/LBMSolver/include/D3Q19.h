////////////////////////////////////////////////////////////////////////////////////
// author: Yang, Zhou. email:yz1249911657@gmail.com
////////////////////////////////////////////////////////////////////////////////////
// tester 1:: Yang, Zhou. email:yz1249911657@gmail.com
////////////////////////////////////////////////////////////////////////////////////
// tester 2: 
////////////////////////////////////////////////////////////////////////////////////

// D3Q19.h
#pragma once
#ifndef D3Q19_H
#define D3Q19_H

#include "CommonAPI.h"
#include "BinBase.h"
#include "Lattice.h"
#include <iostream>

namespace nnwlbm
{
	class COMMONAPI D3Q19 : public nnwlbm::DataClass, public Lattice<3>
	{
	public:
		/**
		 * \brief D3Q19类的构造函数：实现输入一个double类型的数据，将值赋给private中的格子速度c
		 * \param c [in] 外部用户给定格子速度
		 * 示例，D3Q19 dq(1.58),所以c=1.58
		 */
		explicit D3Q19(const double &cc) :c(cc) { }

		/**
		 * \brief 参数输入检查函数
		 * \param i [in] 离散速度空间编号
		 * \return 0-18之内（包含0/18）返回true，其他返回false
		 */
		static inline bool check(const size_t &i)
		{
			if (i > 18) {
				std::cout << "error : index " << i << " is out of range!" << std::endl;
				exit(1);
			}
			else return true;
		}

		/**
		 * \brief 返回D3Q19类型对象的离散速度个数
		 * \return 返回值 19
		 */
		constexpr static inline size_t size() { return 19; }

		/**
		 * \brief 计算D3Q19类型对象的格子速度
		 * \return 外部用户给定的格子速度c
		 */
		inline const double &velocity() const
		{
			return c;
		}

		/**
		 * \brief 计算D3Q19类型对象的迁移方向：
		 * 在debug状态下调用此函数会执行错误检查功能，即调用check函数来判断输入参数是否超出阈值
		 * \param i [in] 离散速度空间编号
		 * \return 返回一个int类型的指针p，(*p)表示x方向，*(p+1)表示y方向,*(p+2)表示z方向
		 */
		static inline const int *direction(const size_t &i)
		{
#ifdef _DEBUG
			check(i);
#endif // _DEBUG
			return (dir + 3 * i);
		}

		/**
		 * \brief 计算D3Q19类型对象的权重系数
		 * 在debug状态下调用此函数会执行错误检查功能，即调用check函数来判断输入参数是否超出阈值;
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
		 * \brief 计算D3Q19类型对象的碰撞后的离散速度空间编号
		 * 在debug状态下调用此函数会执行错误检查功能，即调用check函数来判断输入参数是否超出阈值;
		 * \param i [in] 离散速度空间编号
		 * \return 返回输入离散速度空间编号所对应的反向离散速度空间编号
		 */
		static inline int InvertedDirection(const size_t &i)
		{
#ifdef _DEBUG
			check(i);
#endif // _DEBUG
			return invDir[i];
		}
	public:
		virtual std::string name() const;
	private:
		/**
		 * 格子速度c
		 * 此参数由指定D3Q19对象的时候给定，随着对象的不同，c可以随之改变
		 */
		const double c;

		/**
		 * 正向离散速度空间dir：
		 * 此参数当对象类型指定为D3Q19时就已经确定
		 * dir数组里面的数据按0-18的离散速度空间编号顺序存储
		 * Ref. D. d’Humie`res, I. Ginzburg, M. Krafcyzk, P. Lallemand, and L.-S. Luo, Proc. R. Soc. London, Ser. A 360, 437 (2002).
		 * dir[3*19]={ 0, 0, 0, //  0
									 1, 0, 0, //  1
									-1, 0, 0, //  2
									 0, 1, 0, //  3
									 0,-1, 0, //  4
									 0, 0, 1, //  5
									 0, 0,-1, //  6
									 1, 1, 0, //  7
									-1,-1, 0, //  8
									 1,-1, 0, //  9
									-1, 1, 0, // 10
									 1, 0, 1, // 11
									-1, 0,-1, // 12
									 1, 0,-1, // 13
									-1, 0, 1, // 14
									 0, 1, 1, // 15
									 0,-1,-1, // 16
									 0, 1,-1, // 17
									 0,-1, 1  // 18}
		 */
		static const int dir[3 * 19];

		/**
		 * 反向离散速度空间编号invDir：
		 * 此参数当对象类型指定为D3Q19时就已经确定
		 * invDir数组里面的数据按正向0-8的离散速度空间编号顺序存储其碰撞后的反向离散速度空间编号
		 * invDir[19]={ 0,2,1,4,3,6,5,8,7,10,9,12,11,14,13,16,15,18,17 }
		 */
		static const int invDir[19];

		/**
		 * 权重系数w：
		 * 此参数当对象类型指定为D3Q19时就已经确定
		 * w[19]={1.0 / 3.0,  // 0
							1.0 / 18.0, // 1
							1.0 / 18.0, // 2
							1.0 / 18.0, // 3
							1.0 / 18.0, // 4
							1.0 / 18.0, // 5
							1.0 / 18.0, // 6
							1.0 / 36.0, // 7
							1.0 / 36.0, // 8
							1.0 / 36.0, // 9
							1.0 / 36.0, // 10
							1.0 / 36.0, // 11
							1.0 / 36.0, // 12
							1.0 / 36.0, // 13
							1.0 / 36.0, // 14
							1.0 / 36.0, // 15
							1.0 / 36.0, // 16
							1.0 / 36.0, // 17
							1.0 / 36.0  // 18}
		 */
		static const double w[19];

	private:
		/**
		* \brief D3Q19类二进制输出
		*/
		virtual bool OutputBinPrivate
		(
			std::ofstream &ostr,
			const std::string &filename,
			const std::string &path
		) const;
	};
}
#endif //D3Q19_H