//////////////////////////////////////////////////////////
//Function.cpp
//////////////////////////////////////////////////////////
//author:Guan, Qingshan. email:qingsg@hust.edu.cn
//////////////////////////////////////////////////////////
//participant 1:He, Shuangwu. email:shuangwuh@gmail.com.
//participant 2:Wang, Junquan. email:1935405719@qq.com.
//participant 3:Chen, Songze. email: jacksongze@hust.edu.cn
//////////////////////////////////////////////////////////
//tester 1:He, Shuangwu. email:shuangwuh@gmail.com.
//tester 2:Wang, Junquan. email:1935405719@qq.com.

#include "CommonLibDll.h"
#include "Function.h"
#include <cstring>

namespace nnwlbm {
	::std::string Function::name() const { return "Function"; }

	bool Function::OutputBinPrivate(
		std::ofstream &ostr,  /**< [in,out]二进制文件输出流 */
		const std::string &filename, /**< [in]输出的文件名 */
		const std::string &path /**< [in]输出的路径，默认为当前目录 */
	) const { return true; }

	bool ConstFun::CalVarPrivate(double* va, const double* pos) const {
		memcpy(va, var.data(), sizeof(double) * var.size());
		return true;
	}

	::std::string ConstFun::name()const { return "ConstFun"; }
	bool ConstFun::OutputBinPrivate(
		std::ofstream& ostr,
		const std::string& filename,
		const std::string& path
	) const
	{
		ostr.write((const char*)(&outputDim), sizeof(size_t));
		geometry->output_bin(*geometry,filename+"_geo");
		size_t temp = var.size();
		ostr.write((const char*)(&temp), sizeof(size_t));
		ostr.write((const char*)(var.data()), temp * sizeof(double));
		return true;
	}

	template<>
	std::string BinIn<ConstFun>:: name() const { return "ConstFun"; }

	template<>
	std::shared_ptr<DataClass> BinIn<ConstFun>::InputBinPrivate(
		std::ifstream &istr,
		const std::string &filename,
		const std::string &path
	) const
	{
		size_t dim;
		istr.read((char*)(&dim), sizeof(size_t));
		std::shared_ptr<Geometry> geo = std::dynamic_pointer_cast<Geometry>(binReader.read(filename + "_geo"));
		size_t tempSize;
		istr.read((char*)(&tempSize), sizeof(size_t));
		std::vector<double> doubleVec(tempSize);
		istr.read((char*)(doubleVec.data()), tempSize * sizeof(double));
		return std::make_shared<ConstFun>(geo, doubleVec);
	}



	bool LinearFun::CalVarPrivate(double *va, const double *pos) const {
		va[0] = offset;
		for(size_t i = 0; i < get_input_dim(); i++)
		{
			va[0] += pos[i] * slope[i];
		}
		return true;
	}

	LinearFun::LinearFun(
		const std::shared_ptr<Geometry> &geo,
		const std::vector<double> &slopeIn,
		const double &offsetIn
	) : Function(1, geo), slope(slopeIn), offset(offsetIn)
	{
		if (geo->GetDim() != slope.size()) 
		{
			std::cout << "LinearFun : geo->GetDim() != slope.size()" << std::endl;
#ifdef WIN32
			system("pause");
#endif
			exit(1);
		}
	}

	::std::string LinearFun::name()const { return "LinearFun"; }

	bool LinearFun::OutputBinPrivate(
		std::ofstream &ostr,
		const std::string &filename,
		const std::string &path
	) const
	{
		ostr.write((const char*)(&outputDim), sizeof(size_t));
		geometry->output_bin(*geometry, filename + "_geo");
		size_t temp = slope.size();
		ostr.write((const char*)(&temp), sizeof(size_t));
		ostr.write((const char*)(slope.data()), temp*sizeof(double));
		ostr.write((const char*)(&offset), temp*sizeof(double));
		return true;
	}

	template<>
	std::string BinIn<LinearFun>::name() const { return "LinearFun"; }

	template<>
	std::shared_ptr<DataClass> BinIn<LinearFun>::InputBinPrivate(
		std::ifstream &istr,
		const std::string &filename,
		const std::string &path
	) const
	{
		size_t dim;
		istr.read((char*)(&dim), sizeof(size_t));
		std::shared_ptr<Geometry> geo = std::dynamic_pointer_cast<Geometry>(binReader.read(filename + "_geo"));
		size_t tempSize;
		istr.read((char*)(&tempSize), sizeof(size_t));
		std::vector<double> slope(tempSize);
		double offset;
		istr.read((char*)(slope.data()), tempSize*sizeof(double));
		istr.read((char*)(&offset), sizeof(double));
		return std::make_shared<LinearFun>(geo, slope, offset);
	}


	AssemblyFun::AssemblyFun(const std::vector<std::shared_ptr<Function>> &funVecIn) 
		: Function(**(funVecIn.begin())), funVec(funVecIn) 
	{
		outputDim = CalOutputDim();
		std::shared_ptr<const Geometry> geo = (*(funVec.begin()))->GetDomainofFunction();
		for(const auto& elemFun : funVec)
		{
			elemFun->SetDomainofFunction(geo);
		}
	}

	size_t AssemblyFun::CalOutputDim() const
	{
		size_t totoalOutputDim = 0;
		for(const auto& elemFun : funVec)
		{
			totoalOutputDim += elemFun->GetOutputDim();
		}
		return totoalOutputDim;
	}

	bool AssemblyFun::CalVarPrivate(
		double* var,
		const double* pos
	) const
	{
		for(const auto& elemFun : funVec)
		{
			elemFun->cal_var(var, pos);
			var += elemFun->GetOutputDim();
		}
		return true;
	}
}	// end of namespace nnwlbm
