/////////////////////////////////////////////////
// author: Chen, Songze.  email: jacksongze@hust.edu.cn
/////////////////////////////////////////////////
// collaborator 1 :
/////////////////////////////////////////////////
// test 1 :
///////////////////////////////////////////////

// UniformGrid.h 

#ifndef UNIFORMGRID_H 
#define UNIFORMGRID_H 

#include "CommonAPI.h"
#include "Range.h"
#include "PHIO.h"
//#include "IO.h"
#include <fstream>
#include <iomanip>
#include <sstream>
#include <array>

namespace nnwlbm {
  /**
  *@brief 直角均匀网格类
  * 数据集：步长dx(值保持不变)、原点坐标origin[dim]。
  * 操作集：初始化,设置或返回相对原点的坐标，设置或返回步长值，转换空间坐标，返回某个方向的坐标。
  */
  template<const size_t dim>
  class UniformGrid : public nnwlbm::DataClass {
  private:
    double dx;
    double origin[dim];
  public:
    /** \brief 初始化直角网格的步长dx和坐标原点
    */
    UniformGrid(const double &dxIn, const double *oIn) {
      SetDx(dxIn);
      SetOrigin(oIn);
    }
    /** \brief 通过复制初始化
    */
    UniformGrid(const UniformGrid<dim> &g) {
      SetDx(g.dx);
      SetOrigin(g.origin);
    }
    /** \brief 通过输入值设置相对原点的坐标，dim维
    * \param oIn:[in]预设置的坐标值数组指针
    */
    void SetOrigin(
      const double *oIn
    )
    { memcpy(origin, oIn, sizeof(double)*dim); }
    /** \brief 返回相对原点的坐标，dim维
    * \param oOut:[in][out]预输出的坐标值数组指针
    */
    void get_origin(
      double *oOut
    ) const
    { memcpy(oOut, origin, sizeof(double)*dim); }
    /** \brief 设置步长dx的值
    * \param dxIn:[in]预设置的步长值
    */
    void SetDx(
      const double &dxIn
    )
    { dx = dxIn; }

    /** \brief 返回步长dx的值，grid的常量
    */
    double get_dx() const
    { return dx; }

    /** \brief 取得当前网格的交错网格
    */
    UniformGrid<dim> GetStaggeredGrid() const
    {
      double originSG[dim];
      for (int i = 0; i < dim; i++) {
        originSG[i] = 0.5 * dx + origin[i];
      }
      return UniformGrid<dim>(dx, originSG);
    }
    /** \brief 根据loc的坐标，返回相对原点的坐标至pos中
    * param pos:[in][out]相对坐标原点origin下的坐标
    * param loc:[in]当地坐标，即非物理网格坐标
    */
    void GetPos(
      double *pos,
      const int *loc
    ) const {
      for (int i = 0; i < dim; i++) {
        pos[i] = origin[i] + loc[i] * dx;
      }
    }
    /** \brief 根据loc的坐标，返回iDim方向上相对坐标原点的值
    * \param loc:[in]当地坐标，即非物理网格坐标
    * \param iDim:[in]第iDim维，iDim<dim
    */
    inline double GetPosComp(
      const int *loc,
      const size_t iDim) const {
      return origin[iDim] + loc[iDim] * dx;
    }
    /** \brief 根据loc的坐标，返回iDim方向上相对坐标原点的值
    * \param iLoc:[in]iDim维上当地坐标值
    * \param iDim:[in]第iDim维，iDim<dim
    */
    inline double GetPosComp(
      const int &iLoc,
      const size_t iDim) const {
      return origin[iDim] + iLoc * dx;
    }
  public:
	  /** 返回类型名称. 通过name函数内部定义字符串实现一种类型一个名称 */
	  virtual std::string name() const { return "UniformGrid" + std::to_string(dim); }
  private:
	  /**UniformGrid类的二进制输出函数*/
	  virtual bool OutputBinPrivate(
		  std::ofstream &ostr,
		  const std::string &filename,
		  const std::string &path
	  ) const
	  {
		  ostr.write((const char*)(&dx), sizeof(double));
		  ostr.write((const char*)(origin), 2 * dim * sizeof(double));
		  return true;
	  }
  };

  /**UniformGrid类的二进制读入工具*/
  template<const size_t dim>
  class UniformGridRead : public nnwlbm::BinInput
  {
  public:
	  virtual std::string name() const { return "UniformGrid" + std::to_string(dim); }
  private:
	  virtual std::shared_ptr<nnwlbm::DataClass> InputBinPrivate(
		  std::ifstream &istr, /**< [in,out]二进制文件输入流 */
		  const std::string &filename, /**< [in]读入的文件名 */
		  const std::string &path = ".\\" /**< [in]文件路径，默认路径为当前目录 */
	  ) const {
		  double dx;
		  double origin[dim];
		  istr.read((char*)(&dx), sizeof(double));
		  istr.read((char*)(origin), dim * sizeof(double));
		  return std::make_shared<UniformGrid<dim>>(dx, origin);
	  }
  };

  template<
    typename Grid, 
    const size_t dim, 
    template<typename C, const size_t dim1, typename Compare = nnwlbm::Color, const int ndef = -99999>class MD = MaskData
  >
  bool output_tecplot_FeQuadrilateral(
    const MD<double, dim> &md, // the data to be outputed
    const Grid &grid, // uniform grid
    const std::string &fname = "", // information for title
    const std::string &suffix = "",
    const std::string &titleInfo = "" // information for title
  )
  {
    char str[75];
  #ifdef _WIN32
    sprintf_s(str, "./results/d%sfe%s.dat\0", fname.c_str(), suffix.c_str());
  #else
    sprintf(str, "./results/d%sfe%s.dat\0", fname.c_str(), suffix.c_str());
  #endif
    std::string filename(str);
    std::cout << "writing data... filename : " << filename << std::endl;

    const size_t counter = md.count_num(md.get_comp());
    RangeIndex<dim> rIndex(md.GetRange(), RangeIndex<dim>::inputCellRange);
    Range<dim> cRange = rIndex.GetCellRange();
    Range<dim> nRange = rIndex.GetNodeRange();
    const size_t varNum = md.GetVarNum();
    std::stringstream ss;
    ss << "TITLE = \"" << titleInfo << "\"\n";
    ss << "VARIABLES = ";
    ss << "X0";
    for (int i = 1; i < dim; i++) ss << ", X" << i;
    for (int i = 0; i < varNum; i++) ss << ", V" << i;
    ss << "\nZONE ";

    std::vector<size_t> nodeVec;
    std::map<size_t, size_t> nodeMap;
    std::map<size_t, std::array<size_t, 4>> elemMap;

    int loc[dim];
    for (int carry = cRange.IniCarry(loc); carry == 0; carry = cRange.AddOne(loc))
    {
      size_t cellId = rIndex.CalCellID(loc);

      if (md.IsValid(cellId)) {
        std::array<size_t, 4> nodeList;
        int tmpLoc[dim];
        auto add_node = [&](int id) {
          size_t nodeId = rIndex.CalNodeID(tmpLoc);
          auto iter = nodeMap.find(nodeId);
          if (iter == nodeMap.end())
          {
            nodeVec.push_back(nodeId);
            nodeMap.emplace(nodeId, nodeVec.size());
          }
          nodeList[id] = nodeId;
        };
        tmpLoc[0] = loc[0] + 0; tmpLoc[1] = loc[1] + 0;
        add_node(0);
        tmpLoc[0] = loc[0] + 1; tmpLoc[1] = loc[1] + 0;
        add_node(1);
        tmpLoc[0] = loc[0] + 1; tmpLoc[1] = loc[1] + 1;
        add_node(2);
        tmpLoc[0] = loc[0] + 0; tmpLoc[1] = loc[1] + 1;
        add_node(3);
        elemMap.emplace(cellId, nodeList);
      }
    }
    ss << "NODES = " << nodeVec.size() << ", ELEMENTS = " << elemMap.size();

    ss << "\nDATAPACKING=BLOCK";
    ss << "\nZONETYPE=FEQUADRILATERAL\n";
    ss << "VARLOCATION=([";
    for (int i = 0; i < varNum-1; i++) {
      ss << i + dim + 1 << ',';
    }
    ss << varNum + dim;
    ss << "]=CELLCENTERED)" << std::endl;
    std::string title = ss.str();

    std::ofstream fileOutput;
    fileOutput.open(filename);
    if (!(fileOutput.is_open())) return false;
    fileOutput << title << std::endl;

    fileOutput.setf(std::ios::scientific, std::ios::floatfield);
    fileOutput << std::setprecision(16);
    int breakCounter = 0;
    for (int iDim = 0; iDim < dim; iDim++) {
      breakCounter = 0;
      for (int iNode = 0; iNode < int(nodeVec.size()); iNode++)
      {
        rIndex.CalNodeIJK(loc, nodeVec[iNode]);
        fileOutput << grid.GetPosComp(loc, iDim) << '\t';
        breakCounter++;
        if (breakCounter % 5 == 0) fileOutput << '\n';
      }
      fileOutput << std::endl;
    }
    for (int iVar = 0; iVar < varNum; iVar++) {
      breakCounter = 0;
      for (const auto &elem : elemMap)
      {
        const size_t cellId = elem.first;
        const double *var = md.GetVar(cellId);
        fileOutput << var[iVar] << '\t';
        breakCounter++;
        if (breakCounter % 5 == 0) fileOutput << '\n';
      }
      fileOutput << std::endl;
    }
    fileOutput.unsetf(std::ios::scientific);
    for (const auto &elem : elemMap)
    {
      const auto &nodeList = elem.second;
      for (const auto &nodeId : nodeList) {
        fileOutput << nodeMap.find(nodeId)->second << ' ';
      }
      fileOutput << '\n';
    }
    fileOutput << std::endl;

    return true;
  }

  template<
    typename Grid, 
    const size_t dim, 
    template<typename C, const size_t dim1, typename Compare = nnwlbm::Color, const int ndef = -99999>class MD = MaskData
  >
  bool OutputCellCenteredcellData(
    const MD<double, dim> &md, // the data to be outputed
    const Grid &grid, // uniform grid
    const std::string &fname = "", // information for title
    const std::string &suffix = "",
    const std::string &titleInfo = "" // information for title
  )
  {
    char str[75];
  #ifdef _WIN32
    sprintf_s(str, "./results/d%sfe%s.dat\0", fname.c_str(), suffix.c_str());
  #else
    sprintf(str, "./results/d%sfe%s.dat\0", fname.c_str(), suffix.c_str());
  #endif
    std::string filename(str);
    /////////////////////////////////////////////////////
    std::cout << "writing data... filename : " << filename << '\n';
    /////////////////////////////////////////////////////

    const size_t counter = md.count_num(md.get_comp());
    RangeIndex<dim> rIndex(md.GetRange(), RangeIndex<dim>::inputCellRange);
    Range<dim> cRange = rIndex.GetCellRange();
    Range<dim> nRange = rIndex.GetNodeRange();
    const size_t varNum = md.GetVarNum();
    std::stringstream ss;
    ss << "TITLE = \"" << titleInfo << "\"\n";
    ss << "VARIABLES = ";
    ss << "X0";
    for (int i = 1; i < dim; i++) ss << ", X" << i;
    for (int i = 0; i < varNum; i++) ss << ", V" << i;
    ss << ", color \nZONE ";
    char ijk[] = "IJK";
    ss << "I = " << nRange.right(0) - nRange.left(0);
    for (int i = 1; i < dim; i++)
      ss << ", " << ijk[i] << " = " << nRange.right(i) - nRange.left(i);
    ss << "\nDATAPACKING=BLOCK\n";
    ss << "VARLOCATION=([";
    for (int i = 0; i < varNum; i++) {
      ss << i + dim + 1 << ',';
    }
    ss << varNum + dim + 1;
    ss << "]=CELLCENTERED)" << std::endl;
    std::string title = ss.str();

    std::ofstream fileOutput;
    fileOutput.open(filename);
    if (!(fileOutput.is_open())) return false;
    fileOutput << title << std::endl;

    fileOutput.setf(std::ios::scientific, std::ios::floatfield);
    fileOutput << std::setprecision(9);
    int loc[dim];
    int breakCounter = 0;
    for (int iDim = 0; iDim < dim; iDim++) {
      breakCounter = 0;
      for (int carry = nRange.IniCarry(loc); carry == 0; carry = nRange.AddOne(loc))
      {
        fileOutput << grid.GetPosComp(loc, iDim) << '\t';
        breakCounter++;
        if (breakCounter % 5 == 0) fileOutput << '\n';
      }
      fileOutput << std::endl;
    }
    for (int iVar = 0; iVar < varNum; iVar++) {
      breakCounter = 0;
      for (int carry = cRange.IniCarry(loc); carry == 0; carry = cRange.AddOne(loc))
      {
        size_t cellId = rIndex.CalCellID(loc);
        if (md.IsValid(cellId)) {
          const double *var = md.GetVar(cellId);
          fileOutput << var[iVar] << '\t';
        }
        else {
          fileOutput << 0.0 << '\t';
        }
        breakCounter++;
        if (breakCounter % 5 == 0) fileOutput << '\n';
      }
      fileOutput << std::endl;
    }
    breakCounter = 0;
    for (int carry = cRange.IniCarry(loc); carry == 0; carry = cRange.AddOne(loc))
    {
      size_t cellId = rIndex.CalCellID(loc);
      fileOutput << md.GetColor(cellId) << '\t';
      breakCounter++;
      if (breakCounter % 5 == 0) fileOutput << '\n';
    }
    fileOutput << std::endl;

    return true;
  }

  template<
    typename Grid, 
    const size_t dim, 
    template<typename C, const size_t dim1, typename Compare = nnwlbm::Color, const int ndef = -99999>class MD = MaskData
  >
  bool OutputCellCenteredDataHalfway(
    const MD<double, dim> &md, // the data to be outputed
    const Grid &grid, // uniform grid
    const std::string &fname = "", // information for title
    const std::string &suffix = "",
    const std::string &titleInfo = "" // information for title
  )
  {
    char str[75];
  #ifdef _WIN32
    sprintf_s(str, "./results/d%sfe%s.dat\0", fname.c_str(), suffix.c_str());
  #else
    sprintf(str, "./results/d%sfe%s.dat\0", fname.c_str(), suffix.c_str());
  #endif
    std::string filename(str);
    std::cout << "writing data... filename : " << filename << std::endl;

    const size_t counter = md.count_num(md.get_comp());
    Range<dim> range = md.GetRange();  // cell range
    RangeIndex<dim> rIndex(range, RangeIndex<dim>::inputCellRange);
    const size_t varNum = md.GetVarNum();
    std::stringstream ss;
    ss << "TITLE = \"" << titleInfo << "\"\n";
    ss << "VARIABLES = ";
    ss << "X0";
    for (int i = 1; i < dim; i++) ss << ", X" << i;
    for (int i = 0; i < varNum; i++) ss << ", V" << i;
    ss << ", color \nZONE ";
    char ijk[] = "IJK";
    ss << "I = " << range.right(0) - range.left(0);
    for (int i = 1; i < dim; i++)
      ss << ", " << ijk[i] << " = " << range.right(i) - range.left(i);
    std::string title = ss.str();

    std::ofstream fileOutput;
    fileOutput.open(filename);
    if (!(fileOutput.is_open())) return false;
    fileOutput << title << std::endl;

    fileOutput.setf(std::ios::scientific, std::ios::floatfield);
    fileOutput << std::setprecision(16);
    int loc[dim];
    for (int carry = range.IniCarry(loc); carry == 0; carry = range.AddOne(loc))
    {
      size_t cellId = rIndex.CalCellID(loc);
      double pos[dim];
      grid.GetPos(pos, loc);
      for (size_t i = 0; i < dim; i++) {
        fileOutput << pos[i] << '\t';
      }
      if (md.IsValid(cellId)) {
        const double *var = md.GetVar(cellId);
        for (int i = 0; i < varNum; i++) {
          fileOutput << var[i] << '\t';
        }
      }
      else {
        for (int i = 0; i < varNum; i++) {
          fileOutput << 0.0 << '\t';
        }
      }
      fileOutput << md.GetColor(cellId) << std::endl;
    }
    return true;
  }
};

#endif 
 
// end of CFDT.h
