//////////////////////////////////////////////////////////
//Geometry.cpp
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
#include "LBMGeometry.h"

namespace nnwlbm {
  /*点到圆心的距离等于半径则在边界上，大于半径则在圆外，小于半径则在圆内*/
  int Circle::CalRelativeLocation(const double* pos, const double& tolerance)const {
		
    double delta = (pos[0] - center[0]) * (pos[0] - center[0]) + (pos[1] - center[1]) * (pos[1] - center[1]) - radius * radius;
    if (delta > tolerance) return geo::outside;
    else if (delta < -tolerance) return geo::inside;
    else return geo::boundary;
  }

  bool Rectangle::CalIntersection(double *nod,const Rectangle &rec1, const Rectangle &rec2)const
  {
    double node1[4],node2[4];
	node1[0] = rec1.GetVertex1()[0];
	node1[1] = rec1.GetVertex1()[1];
	node1[2] = rec1.GetVertex2()[0];
	node1[3] = rec1.GetVertex2()[1];
	node2[0] = rec2.GetVertex1()[0];
	node2[1] = rec2.GetVertex1()[1];
	node2[2] = rec2.GetVertex2()[0];
	node2[3] = rec2.GetVertex2()[1];
	if (node1[0] > node1[2])
	{
		double temp[2];
		temp[0] = node1[0];temp[1] = node1[1];node1[0] = node1[2];
		node1[1] = node1[3];node1[2] = temp[0];node1[3] = temp[1];
	}
	if (node2[0] > node2[2])
	{
		double temp[2];
		temp[0] = node2[0];temp[1] = node2[1];node2[0] = node2[2];
		node2[1] = node2[3]; node2[2] = temp[0]; node2[3] = temp[1];
	}
	bool isInter1 = (node1[2] > node2[0]) && (node1[2] < node2[2]) && (node1[3] > node2[1]) && (node1[3] < node2[3]);
	bool isInter2 = (node1[0] > node2[0]) && (node1[0] < node2[2]) && (node1[1] > node2[1]) && (node1[1] < node2[3]);
	bool isInter = isInter1 || isInter2;
	if (isInter) {
		if (isInter1) 
			nod[0] = node2[0], nod[1] = node1[3], nod[2] = node1[2], nod[3] = node2[1];
		if (isInter2)
			nod[0] = node1[0], nod[1] = node2[3], nod[2] = node2[2], nod[3] = node1[1];
	}
	return true;
  }
  /*改变四边形的端点坐标，如果输入的两个端点恒左边或者总左边相同则不能组成四边形*/
  bool Rectangle::SetVertex(const double *d1, const double *d2) {
	  for (int i = 0; i < 2; i++)
	  {
		  vertex1[i] = d1[i];
		  vertex2[i] = d2[i];
	  }
	  return d1[0] < d2[0] && d1[1] < d2[1];
  }

  /*点小于长方形上边界并大于长方形下边界即界内，否则界外*/
  int Rectangle::CalRelativeLocation(const double *pos, const double &tolerance )const {
		
	bool isp2 = pos[0] < vertex2[0] - tolerance && pos[1] < vertex2[1] - tolerance;
	bool isp1 = pos[0] > vertex1[0] + tolerance && pos[1] > vertex1[1] + tolerance;
	bool isp3 = pos[0] < vertex2[0] + tolerance && pos[1] < vertex2[1] + tolerance;
	bool isp4 = pos[0] > vertex1[0] - tolerance && pos[1] > vertex1[1] - tolerance;
	bool ison = isp3 && isp4 && (!(isp1 && isp2));
	if (isp1&&isp2) return geo::inside;
	else if (ison) return geo::boundary;
	else return geo::outside;
  }


  /*设定三点坐标，如果三点在一条直线上则不能组成三角形*/
  bool Triangle::SetVertex(const double *d1, const double *d2, const double *d3) {
    if ((d1[1] - d2[1]) / (d1[0] - d2[0]) == (d2[1] - d3[1]) / (d2[0] - d3[0])) {
      std::cout << "错误！不能构成三角形！" << std::endl;
      return false;
    }
    else
    {
      vertex1[0] = d1[0]; vertex1[1] = d1[1]; vertex2[0] = d2[0]; vertex2[1] = d2[1];
      vertex3[0] = d3[0]; vertex3[1] = d3[1];
      return true;
    }
  }
  /*如果点与三角形三个顶点组成的三个三角形面积之合大于大三角形的面积则在三角形外，等于在三角形边界上，小于在三角形内*/
  int Triangle::CalRelativeLocation(const double *pos, const double &tolerance ) const {

    double abcArea, abdArea, bdcArea, acdArea;
    if ((vertex2[0] - vertex1[0]) * (vertex3[1] - vertex1[1]) - (vertex3[0] - vertex1[0]) * (vertex2[1] - vertex1[1]) < 0)
      abcArea = GetTriangleArea(vertex1, vertex3, vertex2);
    else abcArea = GetTriangleArea(vertex1, vertex2, vertex3);
    if ((vertex2[0] - vertex1[0]) * (pos[1] - vertex1[1]) - (pos[0] - vertex1[0]) * (vertex2[1] - vertex1[1]) < 0)
      abdArea = GetTriangleArea(vertex1, pos, vertex2);
    else abdArea = GetTriangleArea(vertex1, vertex2, pos);
    if ((vertex3[0] - vertex2[0]) * (pos[1] - vertex2[1]) - (pos[0] - vertex2[0]) * (vertex3[1] - vertex2[1]) < 0)
      bdcArea = GetTriangleArea(vertex2, pos, vertex3);
    else bdcArea = GetTriangleArea(vertex2, vertex3, pos);
    if ((vertex3[0] - vertex1[0]) * (pos[1] - vertex1[1]) - (pos[0] - vertex1[0]) * (vertex3[1] - vertex1[1]) < 0)
      acdArea = GetTriangleArea(vertex1, pos, vertex3);
    else acdArea = GetTriangleArea(vertex1, vertex3, pos);
    double sumSuqar = abdArea + bdcArea + acdArea;
    if ((fabs(bdcArea) < tolerance) || (fabs(acdArea) < tolerance) || (fabs(abdArea) < tolerance)) return geo::boundary;
    else if (sumSuqar - abcArea > tolerance) return geo::outside;
    else return geo::inside;
  }


  /*根据椭圆的定义该点到两焦点的距离小于长轴即在椭圆内部*/
  int Ellipse::CalRelativeLocation(const double *pos, const double &tolerance ) const {
		
    double delta1 = sqrt((pos[0] - centerA[0]) * (pos[0] - centerA[0]) + (pos[1] - centerA[1]) * (pos[1] - centerA[1]));
    double delta2 = sqrt((pos[0] - centerB[0]) * (pos[0] - centerB[0]) + (pos[1] - centerB[1]) * (pos[1] - centerB[1]));
    double delta = delta1 + delta2 - 2 * a;
    if (delta > tolerance) return geo::outside;
    else if (delta < -tolerance) return geo::inside;
    else return geo::boundary;
  }

  /*点到球心的距离等于半径则在边界上，大于半径则在球外，小于半径则在球内*/
  int Sphere::CalRelativeLocation(const double *pos, const double &tolerance )const {
    double delta = (pos[0] - center[0]) * (pos[0] - center[0]) + (pos[1] - center[1]) * (pos[1] - center[1])
      + (pos[2] - center[2]) * (pos[2] - center[2]) - radius * radius;
    if (delta > tolerance) return geo::outside;
    else if (delta < -tolerance) return geo::inside;
    else return geo::boundary;
  }


  /*改变长方体的端点坐标，如果输入的两点p1>p2则交换位置，并将p1值赋给vertex1，将p2值赋给vertex2*/
  bool Cube::SetVertex(const double *d1, const double *d2) {
		for (int i = 0; i < 3; i++)
		{
			vertex1[i] = d1[i];
			vertex2[i] = d2[i];
		}
		return d1[0] < d2[0] && d1[1] < d2[1] && d1[2] < d2[2];
  }
  /*点小于长方体上边界并大于长方体下边界即界内，否则界外*/
  int Cube::CalRelativeLocation(const double *pos, const double &tolerance )const {
		
    bool isp2 = pos[0] < vertex2[0] - tolerance && pos[1] < vertex2[1] - tolerance && pos[2] < vertex2[2] - tolerance;
    bool isp1 = pos[0] > vertex1[0] + tolerance && pos[1] > vertex1[1] + tolerance && pos[2] > vertex1[2] + tolerance;
    bool isp3 = pos[0] < vertex2[0] + tolerance && pos[1] < vertex2[1] + tolerance && pos[2] < vertex2[2] + tolerance;
    bool isp4 = pos[0] > vertex1[0] - tolerance && pos[1] > vertex1[1] - tolerance && pos[2] > vertex1[2] - tolerance;
    bool ison = isp3 && isp4 && (!(isp1 && isp2));
    if (isp1&&isp2) return geo::inside;
    else if (ison) return geo::boundary;
    else return geo::outside;
  }

  ::std::string EmptyGeo::name()const { return "EmptyGeo"; }
  bool EmptyGeo::OutputBinPrivate(
    std::ofstream &ostr,
    const std::string &filename,
    const std::string &path
  ) const
  {
    size_t dim = GetDim();
    ostr.write((const char *)(dim), sizeof(int));
    return false;
  }

	template<>
	std::string BinIn<EmptyGeo>::name() const { return "EmptyGeo"; }

	template<>
	std::shared_ptr<DataClass> BinIn<EmptyGeo>::InputBinPrivate(
		std::ifstream &istr,
		const std::string &filename,
		const std::string &path
	) const {
		int tempdim;
		istr.read((char*)(&tempdim), sizeof(int));
		return std::make_shared<EmptyGeo>(tempdim);
	};

	::std::string Circle::name()const { return "Circle"; }
	bool Circle::OutputBinPrivate(
		std::ofstream &ostr,
		const std::string &filename,
		const std::string &path
	) const
	{
		ostr.write((const char*)(center), 2 * sizeof(double));
		ostr.write((const char*)(&radius), sizeof(double));
		return true;
	}

	template<>
	std::string BinIn<Circle>::name() const { return "Circle"; }
	template<>
	std::shared_ptr<DataClass> BinIn<Circle>::InputBinPrivate(
		std::ifstream &istr,
		const std::string &filename,
		const std::string &path
	) const
	{
		double temp[2];
		double tempDouble;
		istr.read((char*)(temp), 2 * sizeof(double));
		istr.read((char*)(&tempDouble), sizeof(double));
		return std::make_shared<Circle>(temp[0], temp[1],tempDouble);
	}
	::std::string Rectangle::name()const { return "Rectangle"; }
	bool Rectangle::OutputBinPrivate(
		std::ofstream &ostr,
		const std::string &filename,
		const std::string &path
	) const
	{
		ostr.write((const char*)(vertex1), 2 * sizeof(double));
		ostr.write((const char*)(vertex2), 2 * sizeof(double));
		return true;
	}

	template<>
	std::string BinIn<Rectangle>::name() const { return "Rectangle"; }

	template<>
	std::shared_ptr<DataClass> BinIn<Rectangle>::InputBinPrivate(
		std::ifstream &istr,
		const std::string &filename,
		const std::string &path
	) const
	{
		double tempVertex1[2], tempCenter[2];
		istr.read((char*)(tempVertex1), 2 * sizeof(double));
		istr.read((char*)(tempCenter), 2 * sizeof(double));
		return std::make_shared<Rectangle>(tempVertex1, tempCenter);
	}

	::std::string Triangle::name()const { return "Triangle"; }
	bool Triangle::OutputBinPrivate(
		std::ofstream &ostr,
		const std::string &filename,
		const std::string &path
	) const
	{
		ostr.write((const char*)(vertex1), 2 * sizeof(double));
		ostr.write((const char*)(vertex2), 2 * sizeof(double));
		ostr.write((const char*)(vertex3), 2 * sizeof(double));
		return true;
	}

	template<>
	std::string BinIn<Triangle>::name() const { return "Triangle"; }

	template<>
	std::shared_ptr<DataClass> BinIn<Triangle>::InputBinPrivate(
		std::ifstream &istr,
		const std::string &filename,
		const std::string &path
	) const
	{
		double tempVertex1[2], tempVertex2[2], tempVertex3[2];
		istr.read((char*)(tempVertex1), 2 * sizeof(double));
		istr.read((char*)(tempVertex2), 2 * sizeof(double));
		istr.read((char*)(tempVertex3), 2 * sizeof(double));
		return std::make_shared<Triangle>(tempVertex1, tempVertex2, tempVertex3);
	}
	::std::string Ellipse::name()const { return "Ellipse"; }
	bool Ellipse::OutputBinPrivate(
		std::ofstream &ostr,
		const std::string &filename,
		const std::string &path
	) const
	{
		ostr.write((const char*)(centerA), 2 * sizeof(double));
		ostr.write((const char*)(centerB), 2 * sizeof(double));
		ostr.write((const char*)(&a), sizeof(double));
		ostr.write((const char*)(&b), sizeof(double));
		return true;
	}

	template<>
	std::string BinIn<Ellipse>::name() const { return "Ellipse"; }

	template<>
	std::shared_ptr<DataClass> BinIn<Ellipse>::InputBinPrivate(
		std::ifstream &istr,
		const std::string &filename,
		const std::string &path
	) const
	{
		double tempCenterA[2], tempCenterB[2], temp[2];
		istr.read((char*)(tempCenterA), 2 * sizeof(double));
		istr.read((char*)(tempCenterB), 2 * sizeof(double));
		istr.read((char*)(temp), 2 * sizeof(double));
		return std::make_shared<Ellipse>(tempCenterA, tempCenterB, temp);
	}

	::std::string Sphere::name()const { return "Sphere"; }
	bool Sphere::OutputBinPrivate(
		std::ofstream &ostr,
		const std::string &filename,
		const std::string &path
	) const
	{
		ostr.write((const char*)(center), 3 * sizeof(double));
		ostr.write((const char*)(&radius), sizeof(double));
		return true;
	}

	template<>
	std::string BinIn<Sphere>::name() const { return "Sphere"; }

	template<>
	std::shared_ptr<DataClass> BinIn<Sphere>::InputBinPrivate(
		std::ifstream &istr,
		const std::string &filename,
		const std::string &path
	) const
	{
		double temp[3], tempDouble;
		istr.read((char*)(temp), 3 * sizeof(double));
		istr.read((char*)(&tempDouble), sizeof(double));
		return std::make_shared<Sphere>(temp, tempDouble);
	}

	::std::string Cube::name()const { return "Cube"; }
	bool Cube::OutputBinPrivate(
		std::ofstream &ostr,
		const std::string &filename,
		const std::string &path
	) const
	{
		ostr.write((const char*)(vertex1), 3 * sizeof(double));
		ostr.write((const char*)(vertex2), 3 * sizeof(double));
		return true;
	}

	template<>
	std::string BinIn<Cube>::name() const { return "Cube"; }

	template<>
	std::shared_ptr<DataClass> BinIn<Cube>::InputBinPrivate(
		std::ifstream &istr,
		const std::string &filename,
		const std::string &path
	) const
	{
		double temp1[3], temp2[3];
		istr.read((char*)(temp1), 3 * sizeof(double));
		istr.read((char*)(temp2), 3 * sizeof(double));
		return std::make_shared<Cube>(temp1, temp2);
	}
		UnionGeo::UnionGeo(std::shared_ptr<Geometry> &g1, std::shared_ptr<Geometry> &g2) : Geometry(g1->GetDim()) {
			add_geo(g1);
			add_geo(g2);
		}

	int UnionGeo::CalRelativeLocation(const double *pos, const double &tolerance ) const {
		int ans = geo::outside;
		for (auto elemGeo : geoVec) {
			int rLoc = elemGeo->CalRelativeLocation(pos, tolerance);
			if (rLoc == geo::inside) {
				ans = geo::inside;
				return ans;
			}
			else if (rLoc == geo::boundary) {
				ans = geo::boundary;
			}
		}
		return ans;
	}

	bool UnionGeo::add_geo(std::shared_ptr<Geometry> &g1) {
		if (g1 == nullptr) return false;
		if (g1->GetDim() != GetDim()) return false;
		geoVec.emplace_back(g1);
		return true;
	}

	::std::string UnionGeo::name() const { return "UnionGeo"; }

	// 几何图形求补集类构造函数
	ComplementGeo::ComplementGeo(const std::shared_ptr<Geometry> &g1) :Geometry(g1->GetDim()), geometry(g1) {};
	// 几何图形求补集类计算相对位置函数
	int ComplementGeo::CalRelativeLocation(const double *pos, const double &tolerance)const
	{
		int type = geometry->CalRelativeLocation(pos, tolerance);
		switch (type)
		{
		case geo::inside: return geo::outside;
		case geo::outside: return geo::inside;
		default: return geo::boundary;
		}
	}

	// 几何图形求交集类构造函数
	IntersectionGeo::IntersectionGeo(std::shared_ptr<Geometry> g1, std::shared_ptr<Geometry> g2) 
		: Geometry(g1->GetDim()) 
	{
		add_geo(g1);
		add_geo(g2);
	}

	// 几何图形求补集类计算相对位置函数
	int IntersectionGeo::CalRelativeLocation(const double *pos, const double &tolerance)const
	{
		int ans = geo::inside;
		for (auto elemGeo : geoVec) {
			int rLoc = elemGeo->CalRelativeLocation(pos, tolerance);
			if (rLoc == geo::outside) {
				ans = geo::outside;
				return ans;
			}
			else if (rLoc == geo::boundary) {
				ans = geo::boundary;
			}
		}
		return ans;
	}

	bool IntersectionGeo::add_geo(std::shared_ptr<Geometry> &g1)
	{
		if (g1 == nullptr) return false;
		if (g1->GetDim() != GetDim()) return false;
		geoVec.emplace_back(g1);
		return true;
	}

	// 几何图形相减类构造函数
	SubtractGeo::SubtractGeo(std::shared_ptr<Geometry> g1, std::shared_ptr<Geometry> g2) :Geometry(g1->GetDim()), geo1(g1), geo2(g2)
	{
		if (g2->GetDim() != g1->GetDim())
		{
			geo2 = std::shared_ptr<UniversalGeo>(new UniversalGeo(g1->GetDim()));
		}
	}
	// 几何图形相减类计算相对位置函数
	int SubtractGeo::CalRelativeLocation(const double *pos, const double &tolerance)const
	{
		int ans = geo::inside;
		int type1 = geo1->CalRelativeLocation(pos, tolerance);
		int type2 = geo2->CalRelativeLocation(pos, tolerance);
		if (type1 == geo::inside&&type2 == geo::outside)
		{
			ans = geo::inside;
		}
		else if (type1 == geo::outside || type2 == geo::inside || type2 == geo::boundary)
		{
			ans = geo::outside;
		}
		else
		{
			ans = geo::boundary;
		}
		return ans;
	}

  // Box

  Box Box::FindBoundingBoxOf(      
    const double *nodes,
    const size_t dim,
    const size_t nodeNum
  ) {
    Box ans(dim);
    for(size_t i = 0; i < dim; i++) {
      ans[i] = ans[i + dim] = nodes[i + dim];
    }
    for(size_t i = 0; i < nodeNum; i++) {
      for(size_t iDim = 0; iDim < 3; iDim++) {
        if(nodes[iDim + 3 * i] < ans[iDim]) {
          ans[iDim] = nodes[iDim + 3 * i];
          continue;
        }
        if(nodes[iDim + 3 * i] > ans[iDim + 3]) {
          ans[iDim + 3] = nodes[iDim + 3 * i];
          continue;
        }
      }
    }
    return ans;
  }

  Box Box::CalIntersection(const Box &r1, const Box &r2) {
    Box ans;
    const size_t dim = r1.GetDim();
    for (int i = 0; i < dim; i++) {
      ans[i] = (r1[i] > r2[i]) ? r1[i] : r2[i];
      ans[i + dim] = (r1[i + dim] > r2[i + dim]) ? r2[i + dim] : r1[i + dim];
    }
    return std::move(ans);
  }

  bool Box::check() const {
    const size_t dim = GetDim();
    for (int i = 0; i < dim; i++) {
      if (box[i] >= box[i + dim]) return false;
    }
    return true;
  }

  int Box::CalRelativeLocation(const double *pos, const double &tolerance) const {
		const size_t dim = GetDim();
		bool isInside = true;
		for (int i = 0; i < dim && isInside; i++) {
			isInside = isInside && (pos[i] > box[i + dim * 0] + tolerance && pos[i] < box[i + dim * 1] - tolerance);
		}
		if (isInside) return geo::inside;
		else {
			bool isOnborder = true;
			for (int i = 0; i < dim && isOnborder; i++) {
				isOnborder = isOnborder && (pos[i] >= box[i + dim * 0] - tolerance && pos[i] <= box[i + dim * 1] + tolerance);
			}
			if (isOnborder) return geo::boundary;
			else return geo::outside;
		}
	}

	bool Geometry::cal_normal(
		double *norm,
		const double *pos,
		const double &tolerance
	) const
	{
		return false;
	}

	bool Geometry::CalIntersection(
		double *intersection,
		const double *p1,
		const double *p2,
		const double &tolerance
	) const {
		return false;
	}

	bool Geometry::CalIntersection(
		double *ans, 
		const Geometry &geoIn, 
		const double *p1, 
		const double *p2,
    const double &tolerance
	) {
		const size_t dim = geoIn.GetDim();
		std::vector<double> tempData(2 * dim);
		double* inPos = tempData.data();
		double* outPos = inPos + dim;
		int loc1, loc2;
		loc1 = geoIn.CalRelativeLocation(p1, tolerance);
		loc2 = geoIn.CalRelativeLocation(p2, tolerance);
		if (loc1 == geo::boundary) {
			for (int i = 0; i < dim; i++) ans[i] = p1[i];
			return true;
		}
		if (loc2 == geo::boundary) {
			for (int i = 0; i < dim; i++) ans[i] = p2[i];
			return true;
		}

		if ((loc1 == geo::inside && loc2 == geo::inside) || (loc1 == geo::outside && loc2 == geo::outside)) {
			return false;
		}

		if (loc1 == geo::inside) {
			for (int i = 0; i < dim; i++) {
				inPos[i] = p1[i]; outPos[i] = p2[i];
			}
		}
		else {
			for (int i = 0; i < dim; i++) {
				inPos[i] = p2[i]; outPos[i] = p1[i];
			}
		}
		for (int j = 0; (loc1 != geo::boundary) && (j < 50); j++) {
			for (int i = 0; i < dim; i++) {
				ans[i] = 0.5 * (inPos[i] + outPos[i]);
			}
			loc1 = geoIn.CalRelativeLocation(ans);
			if (loc1 == geo::inside) {
				for (int i = 0; i < dim; i++) inPos[i] = ans[i];
			}
			else if (loc1 == geo::outside) {
				for (int i = 0; i < dim; i++) outPos[i] = ans[i];
			}
			else {  // loc1 == geo::boundary
				break;
			}
		}
		if (loc1 != geo::boundary) {
			//cout << "cannot find pos on border." << endl;
			//loc1 = geo.CalRelativeLocation(inPos);
			//loc2 = geo.CalRelativeLocation(outPos);
			//return false;
		}
		return true;
	}
}





