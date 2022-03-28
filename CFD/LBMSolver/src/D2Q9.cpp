////////////////////////////////////////////////////////////////////////////////////
// author: Yang, Zhou. email:yz1249911657@gmail.com
////////////////////////////////////////////////////////////////////////////////////

// D2Q9.cpp
#include "CommonLibDll.h"
#include "D2Q9.h"

namespace nnwlbm {

std::string D2Q9::name() const { return "D2Q9"; }

bool D2Q9::OutputBinPrivate
(
  std::ofstream &ostr,
  const std::string &filename,
  const std::string &path
) const
{
  ostr.write((const char*)(&c), 1 * sizeof(double));
  return true;
}

const int D2Q9::dir[18] = { 
   0, 0,  // 0
   1, 0,  // 1
   0, 1,  // 2
  -1, 0,  // 3
   0,-1,  // 4
   1, 1,  // 5
  -1, 1,  // 6
  -1,-1,  // 7
   1,-1   // 8
};

const int D2Q9::invDir[9] = { 0,3,4,1,2,7,8,5,6 };

const double D2Q9::w[9] = { 
  4.0 / 9,
  1.0 / 9,
  1.0 / 9,
  1.0 / 9,
  1.0 / 9,
  1.0 / 36,
  1.0 / 36,
  1.0 / 36,
  1.0 / 36 
};

const double D2Q9MRT::transMatrix[9 * 9]={ 
   1, 1, 1, 1, 1,1, 1, 1, 1,
  -4,-1,-1,-1,-1,2, 2, 2, 2,
   4,-2,-2,-2,-2,1, 1, 1, 1,
   0, 1, 0,-1, 0,1,-1,-1, 1,
   0,-2, 0, 2, 0,1,-1,-1, 1,
   0, 0, 1, 0,-1,1, 1,-1,-1,
   0, 0,-2, 0, 2,1, 1,-1,-1,
   0, 1,-1, 1,-1,0, 0, 0, 0,
   0, 0, 0, 0, 0,1,-1, 1,-1
};

const double D2Q9MRT::invTransMatrix[81]= {  
  1.0/9, -1.0/ 9,  1.0/ 9,      0,       0,      0,       0,      0,      0,
  1.0/9, -1.0/36, -1.0/18,  1.0/6, -1.0/ 6,      0,       0,  1.0/4,      0,
  1.0/9, -1.0/36, -1.0/18,      0,       0,  1.0/6, -1.0/ 6, -1.0/4,      0,
  1.0/9, -1.0/36, -1.0/18, -1.0/6,  1.0/ 6,      0,       0,  1.0/4,      0,
  1.0/9, -1.0/36, -1.0/18,      0,       0, -1.0/6,  1.0/ 6, -1.0/4,      0,
  1.0/9,  1.0/18,  1.0/36,  1.0/6,  1.0/12,  1.0/6,  1.0/12,      0,  1.0/4,
  1.0/9,  1.0/18,  1.0/36, -1.0/6, -1.0/12,  1.0/6,  1.0/12,      0, -1.0/4,
  1.0/9,  1.0/18,  1.0/36, -1.0/6, -1.0/12, -1.0/6, -1.0/12,      0,  1.0/4,
  1.0/9,  1.0/18,  1.0/36,  1.0/6,  1.0/12, -1.0/6, -1.0/12,      0, -1.0/4 
};

template<>
std::string BinIn<D2Q9>::name() const { return "D2Q9"; }

// D2Q9类二进制输入
template<>
std::shared_ptr<DataClass> BinIn<D2Q9>::InputBinPrivate(
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
	return std::make_shared<D2Q9>(c);
}

template<>
std::string BinIn<D2Q9BGK>::name() const { return "D2Q9BGK"; }

// D2Q9BGK类二进制输入
template<>
std::shared_ptr<DataClass> BinIn<D2Q9BGK>::InputBinPrivate(
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
	return std::make_shared<D2Q9>(c);
}

template<>
std::string BinIn<D2Q9MRT>::name() const { return "D2Q9MRT"; }

// D2Q9MRT类二进制输入
template<>
std::shared_ptr<DataClass> BinIn<D2Q9MRT>::InputBinPrivate(
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
	return std::make_shared<D2Q9>(c);
}
}	// end of namespace nnwlbm
