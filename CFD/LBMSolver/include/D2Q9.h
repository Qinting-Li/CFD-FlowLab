////////////////////////////////////////////////////////////////////////////////////
// author: Yang, Zhou. email:yz1249911657@gmail.com
////////////////////////////////////////////////////////////////////////////////////
// participant 1: Wang, junquan. email: as2697920806@gmail.com
////////////////////////////////////////////////////////////////////////////////////
// tester 1:: Yang, Zhou. email:yz1249911657@gmail.com
////////////////////////////////////////////////////////////////////////////////////
// tester 2: He, jiaming. email: 2457837217@qq.com
////////////////////////////////////////////////////////////////////////////////////

// D2Q9.h
#pragma once
#ifndef D2Q9_H
#define D2Q9_H

#include "CommonAPI.h"
#include "BinBase.h"
#include "Lattice.h"
#include <iostream>

// 当前代码都存在nnwlbm名字空间下，便于后期协同操作。
namespace nnwlbm
{
	class COMMONAPI D2Q9 : public nnwlbm::DataClass, public Lattice<2>
	{
	public:
		/**
		 * \brief D2Q9类的构造函数：实现输入一个double类型的数据，将值赋给private中的格子速度c
		 * \param c [in] 外部用户给定格子速度
		 * 示例，D2Q9 dq(1.58),所以c=1.58
		 */
		explicit D2Q9(const double &cc) :c(cc) { }
    
		/**
		 * \brief 参数输入检查函数
		 * \param i [in] 离散速度空间编号
		 * \return 0-8之内（包含0/8）返回true，其他返回false
		 */
		static inline bool check(const size_t &i)
		{
			if (i > 8) { 
        std::cout << "error : index " << i << " is out of range!" << std::endl; 
        exit(1);
      }
			else return true;
		}

		/**
		 * \brief 返回D2Q9类型对象的离散速度个数
		 * \return 返回值 9
		 */
		constexpr static inline size_t size() { return 9; }

		/**
		 * \brief 计算D2Q9类型对象的格子速度
		 * \return 外部用户给定的格子速度c
		 */
		inline const double &velocity() const
		{
			return c;
		}
    
		/**
		 * \brief 计算D2Q9类型对象的迁移方向：
		 * 在debug状态下调用此函数会执行错误检查功能，即调用check函数来判断输入参数是否超出阈值
		 * \param i [in] 离散速度空间编号
		 * \return 返回一个int类型的指针p，(*p)表示x方向，*(p+1)表示y方向
		 */
		static inline const int *direction(const size_t &i)
		{
			#ifdef _DEBUG
			check(i);
			#endif // _DEBUG
			return (dir + 2 * i);
		}

		/**
		 * \brief 计算D2Q9类型对象的权重系数
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
		 * \brief 计算D2Q9类型对象的碰撞后的离散速度空间编号
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
		 * 此参数由指定D2Q9对象的时候给定，随着对象的不同，c可以随之改变
		 */
		const double c;

		/**
		 * 正向离散速度空间dir：
		 * 此参数当对象类型指定为D2Q9时就已经确定
		 * dir数组里面的数据按0-8的离散速度空间编号顺序存储
		 * dir[18]={ 0,0,1,0,0,1,-1,0,0,-1,1,1,-1,1,-1,-1,1,-1 }
		 */
		static const int dir[18];

		/**
		 * 反向离散速度空间编号invDir：
		 * 此参数当对象类型指定为D2Q9时就已经确定
		 * invDir数组里面的数据按正向0-8的离散速度空间编号顺序存储其碰撞后的反向离散速度空间编号
		 * invDir[9]={ 0,3,4,1,2,7,8,5,6 }
		 */
		static const int invDir[9];

		/**
		 * 权重系数w：
		 * 此参数当对象类型指定为D2Q9时就已经确定
		 * w[9]={ 4.0 / 9,1.0 / 9,1.0 / 9,1.0 / 9,1.0 / 9,1.0 / 36,1.0 / 36,1.0 / 36,1.0 / 36 }
		 */
		static const double w[9];

	private:
		/**
		* \brief D2Q9类二进制输出
		*/
		virtual bool OutputBinPrivate
		(
			std::ofstream &ostr,
			const std::string &filename,
			const std::string &path
		) const;
	};

	//===================================计算平衡态分布函数(宏观->微观)==================================
	//																				D2Q9
	//																			/				\
	//																		 /					\
	//																D2Q9BGK				D2Q9MRT
	//=================================================================================================
	/** BGK模型时的计算平衡态分布函数的类，继承自D2Q9类。*/
	class COMMONAPI D2Q9BGK : public D2Q9
	{
	public:
		/**
		 * \brief D2Q9BGK类的构造函数：实现输入一个double类型的数据，将值赋给D2Q9类private中的格子速度c
		 * 示例，D2Q9BGK dq(1.58),所以c=1.58
		 */
		explicit D2Q9BGK(const double &cc) :D2Q9(cc) {  }

		/**
		 * \brief 计算D2Q9 BGK模型的平衡态函数
		 * \param macro [in] 储存宏观量的数组，按照rho、u、v顺序存储
		 * \param i [in] 离散速度空间编号
		 * \return bgk模型下的平衡态分布函数计算值
		*/
		double CalEqui(
      const double *macro, 
      const size_t i
    ) const
		{
			double bgkEq = 0.0;
			double ev;
			double vv;
			ev = macro[1] * (*(direction(i))) + macro[2] * (*(direction(i) + 1));
			vv = macro[2] * macro[2] + macro[1] * macro[1];
			bgkEq = weight(i) * macro[0] * (1.0 + (3.0 * ev + 1.5 * (3.0 * ev * ev - vv) / velocity()) / velocity());
			return bgkEq;
		}
	};

	/** MRT模型时的求解平衡态分布函数的类，继承自D2Q9类。*/
	class COMMONAPI D2Q9MRT : public D2Q9
	{
	public:
		/**
		 * \brief D2Q9MRT类的构造函数：实现输入一个double类型的数据，将值赋给D2Q9类private中的格子速度c
		 * 示例，D2Q9MRT dq(1.58),所以c=1.58
		 */
		explicit D2Q9MRT(const double &cc) :D2Q9(cc) { }

		/**
		 * \brief 计算D2Q9格子构型下MRT模型的变换矩阵
		 * \param i [in] 离散速度空间编号
		 * \return 返回一个int类型的指针p,（*p（i））表示拿到第i行的第一个值
		 */
		static inline const double *trans_matrix(const size_t &i)
		{
			#ifdef _DEBUG
			check(i);
			#endif // _DEBUG
			return (transMatrix + size() * i);
		}

    static inline const double *get_transMatrix() { return transMatrix; }

		/**
		 * \brief 计算D2Q9格子构型下MRT模型的逆变换矩阵
		 * \param i [in] 离散速度空间编号
		 * \return 返回一个double类型的指针p,（*p（i））表示拿到第i行的第一个值
		 */
		static inline const double *inv_trans_matrix(const int &i)
		{
			#ifdef _DEBUG
			check(i);
			#endif // _DEBUG
			return (invTransMatrix + 9 * i);
		}

    static inline const double *get_invTransMatrix() { return invTransMatrix; }

		/**
		 * \brief 计算D2Q9 BGK模型的平衡态函数
		 * \param macro [in] 储存宏观量的数组，按照rho、u、v顺序存储
		 * \param i [in] 离散速度空间编号
		 * \return mrt模型下的平衡态分布函数计算值
		 */
		double CalEqui(
      const double *macro, 
      const size_t i
    ) const
		{
			D2Q9BGK a(velocity());
		  return a.CalEqui(macro, i);
		}

	private:
		/**
		 * D2Q9模型的变换矩阵
		 * 此参数当对象类型指定为D2Q9时就已经确定
		 * transMatrix[81] = { 1, 1, 1, 1, 1,1, 1, 1,1,
							 -4,-1,-1,-1,-1,2, 2, 2,2,
								4,-2,-2,-2,-2,1, 1, 1,1,
								0, 1, 0,-1, 0,1,-1,-1,1,
								0,-2, 0, 2, 0,1,-1,-1,1,
								0, 0, 1, 0,-1,1,1,-1,-1,
								0, 0,-2, 0, 2,1,1,-1,-1,
								0, 1,-1, 1,-1,0,0, 0, 0,
								0, 0, 0, 0, 0,1,-1,1,-1}

		 */
		static const double transMatrix[81];

		/**
		 * D2Q9格子构型MRT模型的逆矩阵；
		 * 此参数当对象类型指定为D2Q9时就已经确定。
		*/
		static const double invTransMatrix[81];
	};
}
#endif // !D2Q9NNWLBM_H

