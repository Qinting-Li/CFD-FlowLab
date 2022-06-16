////////////////////////////////////////////////////////////////////////////////////
// author: Yang, Zhou. email:yz1249911657@gmail.com
////////////////////////////////////////////////////////////////////////////////////

// D3Q19.cpp
#include "CommonLibDll.h"
#include "D3Q19.h"

namespace nnwlbm {

std::string D3Q19::name() const { return "D3Q19"; }

bool D3Q19::OutputBinPrivate
(
	std::ofstream &ostr,
	const std::string &filename,
	const std::string &path
) const
{
	ostr.write((const char*)(&c), 1 * sizeof(double));
	return true;
}

const int D3Q19::dir[3 * 19] = {
	 0, 0, 0, //  0
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
	 0,-1, 1  // 18
};

const int D3Q19::invDir[19] = { 0,2,1,4,3,6,5,8,7,10,9,12,11,14,13,16,15,18,17 };

const double D3Q19::w[19] = {
	1.0 / 3.0,  // 1
	1.0 / 18.0, // 2
	1.0 / 18.0, // 3
	1.0 / 18.0, // 4
	1.0 / 18.0, // 5
	1.0 / 18.0, // 6
	1.0 / 18.0, // 7
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
	1.0 / 36.0, // 18
	1.0 / 36.0  // 19
};

template<>
std::string BinIn<D3Q19>::name() const { return "D3Q19"; }

// D3Q19类二进制输入
template<>
std::shared_ptr<DataClass> BinIn<D3Q19>::InputBinPrivate(
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
	return std::make_shared<D3Q19>(c);
}
}	// end of namespace nnwlbm

