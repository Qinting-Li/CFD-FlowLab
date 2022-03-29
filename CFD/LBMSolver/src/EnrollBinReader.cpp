#include "CommonLibDll.h"
#include "BinBase.h"
#include "GlobalData.h"
#include "Range.h"
#include "UniformGrid.h"
#include "Function.h"
#include "LBMGeometry.h"
#include "D2Q9.h"
#include "D3Q15.h"
#include "D3Q19.h"

#ifdef _WIN32
extern "C" _declspec(dllexport)
#endif
void enroll_bin_reader(nnwlbm::BinReader& reader) {
  //! 模板参数既可以是DataClass的派生类T（要求定义了相应的BinIn<T>），也可以是BinInput的派生类
  nnwlbm::enroll_binInput<
    nnwlbm::Circle,
    nnwlbm::Rectangle,
    nnwlbm::Triangle,
    nnwlbm::Cube,
    nnwlbm::Sphere,
    nnwlbm::Ellipse,
    nnwlbm::EmptyGeo,
    nnwlbm::ConstFun,
    nnwlbm::LinearFun,
    nnwlbm::GlobalData,
    nnwlbm::RangeRead<1>,
    nnwlbm::RangeRead<2>,
    nnwlbm::RangeRead<3>,
    nnwlbm::RangeRead<4>,
    nnwlbm::MaskRangeRead<int, 1>,
    nnwlbm::MaskRangeRead<int, 2>,
    nnwlbm::MaskRangeRead<int, 3>,
    nnwlbm::MaskRangeRead<int, 4>,
    nnwlbm::MaskRangeRead<double, 1>,
    nnwlbm::MaskRangeRead<double, 2>,
    nnwlbm::MaskRangeRead<double, 3>,
    nnwlbm::MaskRangeRead<double, 4>,
    nnwlbm::UniformGridRead<1>,
    nnwlbm::UniformGridRead<2>,
    nnwlbm::UniformGridRead<3>,
    nnwlbm::UniformGridRead<4>,
		nnwlbm::D2Q9,
		nnwlbm::D2Q9BGK,
		nnwlbm::D2Q9MRT,
		nnwlbm::D3Q15,
		nnwlbm::D3Q15BGK,
		nnwlbm::D3Q15MRT,
		nnwlbm::D3Q19
  >(reader);
}

void nnwlbm::prepare_nnwlbm_bin_reader() { enroll_bin_reader(nnwlbm::binReader); }

#ifdef _WIN32
bool nnwlbm::enroll_bin_reader_from(nnwlbm::BinReader& r, const std::string& dllName)
{
  HINSTANCE tempDll = LoadLibraryA(dllName.c_str());
  if (tempDll == nullptr)
  {
    std::cout << "Cannot find " << dllName << "!" << std::endl;
    return false;
  }
  else {
    //  定义函数指针fp，并通过从动态链接库中取名为“fun”的函数地址赋值给fp
    auto fp = (nnwlbm::EnrollBinReader)GetProcAddress(tempDll, "enroll_bin_reader");
    if (fp != nullptr) {
      fp(r);
      return true;
    }
    else {
      std::cout << "Cannot enroll bin reader." << std::endl; 
      return false;
    }
  }
}
#else
#include <dlfcn.h>
bool nnwlbm::enroll_bin_reader_from(nnwlbm::BinReader& r, const std::string& dllName)
{
	void* tempLib = dlopen(dllName.c_str(), RTLD_NOW);
	if(tempLib == nullptr) {
    std::cout << "Cannot find " << dllName << "!" << std::endl;
    return false;
  }
  else {
    //  定义函数指针fp，并通过从动态链接库中取名为“fun”的函数地址赋值给fp
    auto fp = (nnwlbm::EnrollBinReader)dlsym(tempLib, "enroll_bin_reader");
    if (fp != nullptr) {
      fp(r);
      return true;
    }
    else {
      std::cout << "Cannot enroll bin reader." << std::endl; 
      return false;
    }
  }
}
#endif
