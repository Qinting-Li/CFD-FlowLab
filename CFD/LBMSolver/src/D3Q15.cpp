////////////////////////////////////////////////////////////////////////////////////
// author: Yang, Zhou. email:yz1249911657@gmail.com
////////////////////////////////////////////////////////////////////////////////////

// D3Q15.cpp
#include "CommonLibDll.h"
#include "D3Q15.h"

namespace nnwlbm {

std::string D3Q15::name() const { return "D3Q15"; }

bool D3Q15::OutputBinPrivate
(
	std::ofstream &ostr,
	const std::string &filename,
	const std::string &path
) const
{
	ostr.write((const char*)(&c), 1 * sizeof(double));
	return true;
}

const int D3Q15::dir[45] = {
   0, 0, 0, //  0
   1, 0, 0, //  1
  -1, 0, 0, //  2
   0, 1, 0, //  3
   0,-1, 0, //  4
   0, 0, 1, //  5
   0, 0,-1, //  6
   1, 1, 1, //  7
  -1, 1, 1, //  8
   1,-1, 1, //  9
  -1,-1, 1, // 10
   1, 1,-1, // 11
  -1, 1,-1, // 12
   1,-1,-1, // 13
  -1,-1,-1  // 14
};

const int D3Q15::invDir[15] = {0,2,1,4,3,6,5,14,13,12,11,10,9,8,7};

const double D3Q15::w[15] = { 
  2.0/9.0, 
  1.0/9.0, 
  1.0/9.0, 
  1.0/9.0, 
  1.0/9.0, 
  1.0/9.0, 
  1.0/9.0, 
  1.0/72.0,
  1.0/72.0,
  1.0/72.0, 
  1.0/72.0,
  1.0/72.0, 
  1.0/72.0,
  1.0/72.0,
  1.0/72.0 
};

const double D3Q15MRT::transMatrix[15 * 15] = {	
   1,  1,  1,	 1,	 1,	 1,	 1,	1,	1,	1,	1,	1,	1,	1,	1,
  -2, -1, -1, -1, -1, -1, -1,	1,	1,	1,	1,	1,	1,	1,	1,
  16, -4, -4, -4, -4, -4, -4,	1,	1,	1,	1,	1,	1,	1,	1,
   0,  1, -1,  0,  0,  0,  0, 1, -1,	1, -1,	1, -1,	1, -1,
   0, -4,	 4,  0,  0,  0,  0, 1, -1,	1, -1,	1, -1,	1, -1,
   0,  0,  0,  1, -1,  0,  0,	1,	1, -1, -1,	1,	1, -1, -1,
   0,  0,  0, -4,  4,  0,  0,	1,	1, -1, -1,	1,	1, -1, -1,
   0,  0,  0,  0,  0,  1, -1,	1,	1,	1,	1, -1, -1, -1, -1,
   0,  0,  0,  0,  0, -4,  4,	1,	1,	1,	1, -1, -1, -1, -1,
   0,  2,  2, -1, -1, -1, -1,	0,	0,	0,	0,	0,	0,	0,	0,
   0,  0,  0,  1,  1, -1, -1,	0,	0,	0,	0,	0,	0,	0,	0,
   0,  0,  0,  0,  0,  0,  0,	1, -1, -1,	1,	1, -1, -1,	1,
   0,  0,  0,  0,  0,  0,  0,	1,	1, -1, -1, -1, -1,	1,	1,
   0,  0,  0,  0,  0,  0,  0,	1, -1,	1, -1, -1,	1, -1,	1,
   0,  0,  0,  0,  0,  0,  0,	1, -1, -1,	1, -1,	1,	1, -1
}; 

template<>
std::string BinIn<D3Q15>::name() const { return "D3Q15"; }

// D3Q15类二进制输入
template<>
std::shared_ptr<DataClass> BinIn<D3Q15>::InputBinPrivate(
	std::ifstream &istr, /**< [in,out]二进制文件输入流 */
	const std::string &filename, /**< [in]读入的文件名 */
	const std::string &path /**< [in]文件路径，默认路径为当前目录 */
) const
{
	// 临时内存空间
	double c;
	// 读入数据
	istr.read((char*)(&c), sizeof(double));
	// 使用合适的构造函数 构造智能指针
	return std::make_shared<D3Q15>(c);
}

template<>
std::string BinIn<D3Q15BGK>::name() const { return "D3Q15BGK"; }

// D3Q15BGK类二进制输入
template<>
std::shared_ptr<DataClass> BinIn<D3Q15BGK>::InputBinPrivate(
	std::ifstream &istr, /**< [in,out]二进制文件输入流 */
	const std::string &filename, /**< [in]读入的文件名 */
	const std::string &path /**< [in]文件路径，默认路径为当前目录 */
) const
{
	// 临时内存空间
	double c;
	// 读入数据
	istr.read((char*)(&c), sizeof(double));
	// 使用合适的构造函数 构造智能指针
	return std::make_shared<D3Q15>(c);
}

template<>
std::string BinIn<D3Q15MRT>::name() const { return "D3Q15MRT"; }

// D3Q15MRT类二进制输入
template<>
std::shared_ptr<DataClass> BinIn<D3Q15MRT>::InputBinPrivate(
	std::ifstream& istr, /**< [in,out]二进制文件输入流 */
	const std::string& filename, /**< [in]读入的文件名 */
	const std::string& path /**< [in]文件路径，默认路径为当前目录 */
) const
{
	// 临时内存空间
	double c;
	// 读入数据
	istr.read((char*)(&c), sizeof(double));
	// 使用合适的构造函数 构造智能指针
	return std::make_shared<D3Q15>(c);
}
}	// end of namespace nnwlbm
