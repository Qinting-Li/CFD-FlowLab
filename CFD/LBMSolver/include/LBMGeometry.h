//////////////////////////////////////////////////////////
//Geometry.h
//////////////////////////////////////////////////////////
//author:Guan, Qingshan. email:qingsg@hust.edu.cn
//////////////////////////////////////////////////////////
//participant 1:He, Shuangwu. email:shuangwuh@gmail.com.
//participant 2:Wang, Junquan. email:1935405719@qq.com.
//participant 3:Chen, Songze. email: jacksongze@hust.edu.cn
//////////////////////////////////////////////////////////
//tester 1:He, Shuangwu. email:shuangwuh@gmail.com.
//tester 2:Wang, Junquan. email:1935405719@qq.com.

#ifndef GEOMETRY_H
#define GEOMETRY_H
#include <math.h>
#include <iostream>
#include <vector>
#include "BinBase.h"
#include "CommonAPI.h"

namespace nnwlbm {
	/**
	* 枚举；
	* \param inside 值为0，表示内部；
	* \param boundary 值为1，表示边界；
	* \param outside 值为2，表示外部。
	*/
	enum geo { 
    ndef = -99999, 
    outside = -1, 
    boundary = 0, 
    inside = 1,
    neighbor = 2
  };
	/**
	  * Geometry基类，继承DataClass
	  * 其他二维图形和三维图形的基类
	  * 私有数据包括有维度dimension
	  */
	class COMMONAPI Geometry : public DataClass
	{
	public:
		/**
		 * Geometry类的构造函数：
		 * 输入一个int类型的数据，将值赋给private中的维度dimension；
		 * \param dimIn 外部用户给定的维度;
		 * \param geo 外部用户给定的图形，并将此图形的维度赋值到dimension；
		 * 示例，Geometry Geo(2),所以dimension=2。
		 */
		explicit Geometry(const size_t &dimIn) : dimension(dimIn) {}
		Geometry(const Geometry &geoIn) : dimension(geoIn.dimension) {}
		/**
		 * 参数输入计算纯虚函数，计算该点是否在图形内：
		 * 实现输入一个点数组，包含各维度的坐标值；
		 * \param tolerance为所允许误差，默认1e-9；
		 * return int值点在内部返回0，边界上返回1，外部返回2。
		 */
		virtual int CalRelativeLocation(
      const double *pos, 
      const double &tolerance = 1e-9
    ) const = 0;

    virtual bool cal_normal(
      double *norm, 
      const double *pos, 
      const double &tolerance = 1e-9
    ) const;

    virtual bool CalIntersection(
      double *intersection, 
      const double *p1, 
      const double *p2,
      const double &tolerance = 1e-9
    ) const;

    static bool CalIntersection(
      double *ans, 
      const Geometry &geoIn, 
      const double *pos1, 
      const double *pos2,
      const double &tolerance = 1e-9
    );


    /**
		 * 输出维度dimension：
		 * return int类型dimension值。
		 */
    const size_t GetDim() const { return dimension; }
	private:
		/**
		* 维度dimension；
		* 此参数由指定Geometry对象的时候给定，随着对象的不同，dimension可以随之改变；
		* int类型。
		*/
		const size_t dimension;
		/**
		* \brief 虚函数输出当前对象的二进制文件数据；
		* 实现自己的二进制输出；
		* \param ostr [in,out]二进制文件输出流；
		* \param filename [in]输出的文件名；
		* \param path [in]输出的路径，默认为当前目录；
		* \return 返回bool变量表示输出状态，正常为true，异常为false。
		*/
		virtual bool OutputBinPrivate(
			std::ofstream &ostr,
			const std::string &filename,
			const std::string &path
		) const
		{
			ostr.write((const char *)(&dimension), sizeof(int));
			return false;
		}
	};
	/**
	  * 空Geometry,即维度为0；
	  * 计算是否在内部时总是返回outside；
	  * 可用于检查测试，继承自Geometry。
	  */
	class COMMONAPI EmptyGeo : public Geometry {
	public:
		/**
		*返回类型名称. 通过name函数内部定义字符串实现一种类型一个名称。
		*/
		virtual ::std::string name()const;
		/**
		 * EmptyGeo类的构造函数：
		 * \param d 给定的维度int类型dimension，恒为0;
		 */
		EmptyGeo(const int d = 0) : Geometry(d) {}
		/**
		 * 参数输入检查纯虚函数，计算该点是否在图形内，由于空图形，所以恒为2，即外部：
		 * 实现输入一个点数组，包含各维度的坐标值；
		 * \param tolerance为所允许误差，默认1e-9；
		 * return 枚举中的outside，即为2。
		 */
		virtual int CalRelativeLocation(const double *, const double &tolerance = 1e-9) const {
      return geo::outside; 
    }
	private:
		/**
	   * \brief 虚函数输出当前对象的二进制文件数据；
	   * 实现自己的二进制输出；
	   * \param ostr [in,out]二进制文件输出流；
	   * \param filename [in]输出的文件名；
	   * \param path [in]输出的路径，默认为当前目录；
	   * \return 返回bool变量表示输出状态，正常为true，异常为false。
	   */
		virtual bool OutputBinPrivate(
			std::ofstream &ostr,
			const std::string &filename,
			const std::string &path
		) const;

	};

	  /**
	  * 圆形，继承自Geometry基类；
		* 通过一点以及一个半径构造一个圆形；
		* 私有成员包括圆心点的double类型的坐标数组center以及double类型的半径radius。
		*/
	class COMMONAPI Circle :public Geometry
	{
	public:
		/**
		* Circle的构造函数；
		* 输入一个double的数组赋值给圆心坐标，一个double值赋值给半径；
		* 默认维度为2维，不可更改；
		* 默认值半径为1.0，圆心坐标（0.0，0.0）可更改。
		*/
		explicit Circle() :Geometry(2), radius(1.0) { center[0] = 0.0; center[1] = 0.0; }
		explicit Circle(const double *c, const double &r) :Geometry(2)
		{
			SetCenter(c);
			SetRadius(r);
		}
		explicit Circle(const double &c1, const double &c2, const double &r) :Geometry(2)
		{
			SetCenter(c1, c2);
			SetRadius(r);
		}
		/**
    *返回类型名称. 通过name函数内部定义字符串实现一种类型一个名称。
    */
		virtual ::std::string name()const;
		/**
		* 圆心坐标设置内联函数；
		* 输入一个大小为2的double数组，实现给中心点坐标赋值；
		* 无返回值。
		*/
		inline void SetCenter(const double *c) { center[0] = c[0], center[1] = c[1]; }
		inline void SetCenter(const double &c1, const double &c2) { center[0] = c1, center[1] = c2; }
		/**
		* 半径设置内联函数;
		* 输入一个double值，赋值给半径;
		* 无返回值.
		*/
		inline void SetRadius(const double &r) { radius = r; }
		/**
		* 参数输入计算虚函数，计算该点是否在圆内：
		* 实现输入一个点坐标数组，以及误差tolerance；
		* \param * 要计算的点坐标数组；
		* \param tolerance为所允许误差，默认1e-9；
		* return int值点在内部返回0，边界上返回1，外部返回2。
		*/
		virtual int CalRelativeLocation(const double *, const double &tolerance = 1e-9)const;
		/**
		* 输出图形中的圆心坐标：
		* \return double值，圆心坐标数组，外部不可更改。
		*/
		inline const double *GetCenter()const { return center; }
		/**
		* 输出图形中的半径：
		* \return double值，圆心半径，外部不可更改。
		*/
		inline const double GetRadius()const { return radius; }
	private:
		/**
		* \param center 圆心点的坐标数组；
		* \param radius 半径；
		* 都为double类型。
		*/
		double center[2], radius;
		/**
		* \brief 虚函数输出当前对象的二进制文件数据；
		* 实现自己的二进制输出；
		* \param ostr 二进制文件输出流；
		* \param filename [in]输出的文件名；
		* \param path [in]输出的路径，默认为当前目录；
		* \return 返回bool变量表示输出状态，正常为true，异常为false。
		*/
		virtual bool OutputBinPrivate(
			std::ofstream &ostr,
			const std::string &filename,
			const std::string &path
		) const;

	
	};
  	/**
  	* 矩形，继承自Geometry；
		* 通过对角线上的两个顶点，构造一个矩形；
		* 私有数据包括两个顶点的坐标数组。
		*/
	class COMMONAPI Rectangle :public Geometry
	{
	public:
		/**
		* Rectangle构造函数，默认维度为2维，不可更改；
		* 默认坐标如下，可更改；
							1(-1.0,1.0)       4(1.0,1.0)
									------------
									|          |
									|     .    |
									|  center  |
									|   (0,0)  |
									------------
							 3(-1.0,-1.0)     2(1.0,-1.0)
		*/
		explicit Rectangle() :Geometry(2)
		{
			vertex1[0] = -1.0; vertex1[1] = 1.0; vertex2[0] = 1.0; vertex2[1] = -1.0;
		}
		explicit Rectangle(const double *p1, const double *p2) :Geometry(2)
		{
			SetVertex(p1, p2);
		}
		/**
		* 返回类型名称. 通过name函数内部定义字符串实现一种类型一个名称。
		*/
		virtual ::std::string name()const;
		/**
		* 更改坐标函数
		* 输入两个double数组，赋值给两个对角顶点坐标；
		* 两个顶点坐标必须严格满足第一个点的横坐标和纵坐标都小于第二个点；
									              2
									------------
									|          |
									|          |
									|          |
									|          |
									------------
							 1    
		* \param *d1,*d2 顶点坐标；
		* \return 返回bool变量表示输出状态，正常为true，异常为false。
		*/
		virtual bool SetVertex(const double *d1, const double *d2);
		/**
		* 参数输入计算虚函数，计算该点是否在矩形内：
		* 实现输入一个点坐标数组，以及误差tolerance；
		* \param * 要计算的点坐标数组；
		* \param tolerance为所允许误差，默认1e-9；
		* return int值点在内部返回0，边界上返回1，外部返回2。
		*/
		virtual int CalRelativeLocation(const double *, const double &tolerance = 1e-9)const;
		/**
		* 参数输入计算虚函数，计算两个矩形的相交点：
		* 实现输入两个矩形；
		* \param rec1，rec2[in] 两个矩形；
		* \param nod[in,out] 存放交点的数组，存放顺序为，x1,y1,x2,y2；
		* return 返回bool变量表示输出状态，正常为true，异常为false。
		*/
		virtual bool CalIntersection(double *nod, const Rectangle &rec1, const Rectangle &rec2)const;
		/**
		* 输出图形中的1顶点坐标：
		* \return 返回double值，1顶点坐标，外部不可更改。
		*/
		inline const double *GetVertex1()const { return vertex1; }
		/**
		* 输出图形中的2顶点坐标：
		* \return 返回double值，2顶点坐标，外部不可更改。
		*/
		inline const double *GetVertex2()const { return vertex2; }
	private:
		/**
		* 各点坐标；
		* \param vertex1-2 两个顶点坐标数组；
		* 都为double类型。
		*/
		double vertex1[2], vertex2[2];
		/**
		* \brief 虚函数输出当前对象的二进制文件数据；
		* 实现自己的二进制输出；
		* \param ostr 二进制文件输出流；
		* \param filename [in]输出的文件名；
		* \param path [in]输出的路径，默认为当前目录；
		* \return 返回bool变量表示输出状态，正常为true，异常为false。
		*/
		virtual bool OutputBinPrivate(
			std::ofstream &ostr,
			const std::string &filename,
			const std::string &path
		) const;
	
	};
    /**
    *三角形，继承自Geometry；
	  * 通过三个顶点，构造一个三角形；
	  * 私有数据包括三个顶点的坐标。
	  */
	class COMMONAPI Triangle : public Geometry
	{
	public:
		/**
		* Triangle构造函数；
		* 默认维度二维，不可更改；
		* 默认顶点坐标如下，可更改；

									3(0,1)
								  		、
								  	/  、
									 /     、
									/       、
								 /         、
								 ——————
							 1(-1,-1)      2(1,-1)
		*/
		explicit Triangle() :Geometry(2) {
			vertex1[0] = -1; vertex1[1] = -1; vertex2[0] = 1;
			vertex2[1] = -1; vertex3[0] = 0; vertex3[1] = 1;
		}
		explicit Triangle(const double *p1, const double *p2, const double *p3) :Geometry(2)
		{
			SetVertex(p1, p2, p3);
		}
		/**
		*返回类型名称. 通过name函数内部定义字符串实现一种类型一个名称。
		*/
		virtual ::std::string name()const;
		/**
		* 顶点坐标设定函数；
		* 输入三个double数组，分别给三个顶点赋值；
		* \param *d1-3[in] 三个顶点的坐标数组；
		* 判断是否能够组成三角形，如果不能则输出“错误！不能构成三角形！”。
		*/
		virtual bool SetVertex(const double *d1, const double *d2, const double *d3);
		/**
		* 计算三角形面积函数；
		* 输入三个double数组，分别给三个顶点赋值；
		* \param *p1-3[in] 三个顶点的坐标数组。
		*/
		inline double GetTriangleArea(const double *p1, const double *p2, const double *p3)const {
			return (p1[0] * p2[1] + p1[1] * p3[0] + p2[0] * p3[1] - p1[0] * p3[1] - p1[1] * p2[0] - p2[1] * p3[0]) / 2.0;
		}
		/**
		* 参数输入计算虚函数，计算该点是否在三角形内：
		* 实现输入一个点坐标数组，以及误差tolerance；
		* \param * 要计算的点坐标数组；
		* \param tolerance为所允许误差，默认1e-9；
		* return int值点在内部返回0，边界上返回1，外部返回2。
		*/
		virtual int CalRelativeLocation(const double *, const double &tolerance = 1e-9)const;
		/**
		* 输出图形中的1顶点坐标：
		* \return 返回double值，1顶点坐标，外部不可更改。
		*/
		inline const double *GetVertex1()const { return vertex1; }
		/**
		* 输出图形中的2顶点坐标：
		* \return 返回double值，2顶点坐标，外部不可更改。
		*/
		inline const double *GetVertex2()const { return vertex2; }
		/**
		* 输出图形中的3顶点坐标：
		* \return 返回double值，3顶点坐标，外部不可更改。
		*/
	inline const double *GetVertex3()const { return vertex3; }	private:
		/**
		* \param vertex1-3 三个顶点的坐标；
		* 都为double类型。
		*/
		double vertex1[2], vertex2[2], vertex3[2];
		/**
		* \brief 虚函数输出当前对象的二进制文件数据；
		* 实现自己的二进制输出；
		* \param ostr 二进制文件输出流；
		* \param filename [in]输出的文件名；
		* \param path [in]输出的路径，默认为当前目录；
		* \return 返回bool变量表示输出状态，正常为true，异常为false。
		*/
		virtual bool OutputBinPrivate(
			std::ofstream &ostr,
			const std::string &filename,
			const std::string &path
		) const;
	
	};
	  /**
	  *椭圆，继承自Geometry；
		* 通过设定两个焦点以及长轴短轴，构造一个椭圆；
		* 私有成员包括焦点坐标数组以及长轴和短轴值。
		*/
	class COMMONAPI Ellipse :public Geometry
	{
	public:
		/**
		* Ellipse构造函数；
		* 默认维度为二维，不可更改；
		* 默认两个焦点坐标（-1.0，0.0）和（1.0，0.0），默认长轴9.0，短轴4.0，可更改。
		*/
		explicit Ellipse() :Geometry(2) { centerA[0] = -1.0; centerA[1] = 0.0; centerB[0] = 1.0; centerB[1] = 0.0; a = 9.0; b = 4.0; }
		explicit Ellipse(const double *p1, const double *p2, const double *p3) :Geometry(2)
		{
			SetCenter(p1, p2);
			SetLength(p3[0], p3[1]);
		}
		/**
		*返回类型名称. 通过name函数内部定义字符串实现一种类型一个名称。
		*/
		virtual ::std::string name()const;
		/**
		* 更改焦点坐标内联函数；
		* \param *center1-2 焦点的坐标数组；
		* 都为double类型。
		*/
		inline void SetCenter(const double *center1, const double *center2) {
			centerA[0] = center1[0];
			centerA[1] = center1[1];
			centerB[0] = center2[0];
			centerB[1] = center2[1];
		}
		/**
		* 更改长轴短轴内联函数；
		* \param lengthA 长轴；
		* \param lengthB 短轴；
		* 都为double 类型。
		*/
		inline void SetLength(const double &lengthA, const double &lengthB) {
			a = (lengthA > lengthB) ? lengthA : lengthB;
			b = (lengthA < lengthB) ? lengthA : lengthB;
		}
		/**
		* 参数输入计算虚函数，计算该点是否在椭圆内：
		* 实现输入一个点坐标数组，以及误差tolerance；
		* \param * 要计算的点坐标数组；
		* \param tolerance为所允许误差，默认1e-9；
		* return int值点在内部返回0，边界上返回1，外部返回2。
		*/
		virtual int CalRelativeLocation(const double *, const double &tolerance = 1e-9) const;
		/**
		* 以指针形式输出图形中的1焦点坐标：
		* 输出顺序，依次输出1焦点的横坐标，纵坐标；
		* \return 返回double值，1焦点坐标，不可更改。
		*/
		inline const double *GetCenterA()const { return centerA; }
		/**
		* 以指针形式输出图形中的2焦点坐标：
		* 输出顺序，依次输出2焦点的横坐标，纵坐标；
		* \return 返回double值，2焦点坐标，不可更改。
		*/
		inline const double *GetCenterB()const { return centerB; }
		/**
		* 输出图形中的长轴值：
		* \return 返回double值，长轴值，不可更改。
		*/
		inline const double GetA()const { return a; }
		/**
		* 输出图形中的短轴值：
		* \return 返回double值，短轴值，不可更改。
		*/
		inline const double GetB()const { return b; }
	private:
		/**
		* \param centerA，centerB 两个焦点的坐标数组；
		* \param a,b 长轴短轴值；
		* \param 都为double类型。
		*/
		double centerA[2], centerB[2], a, b;
		/**
		* \brief 虚函数输出当前对象的二进制文件数据；
		* 实现自己的二进制输出；
		* \param ostr 二进制文件输出流；
		* \param filename [in]输出的文件名；
		* \param path [in]输出的路径，默认为当前目录；
		* \return 返回bool变量表示输出状态，正常为true，异常为false。
		*/
		virtual bool OutputBinPrivate(
			std::ofstream &ostr,
			const std::string &filename,
			const std::string &path
		) const;
	
	};

	  /**
	  * 球，继承自Geometry，三维；
		* 通过设定球心和半径，构造一个球体；
		* 私有成员包括球心点坐标数组以及半径。
		*/
	class COMMONAPI Sphere :public Geometry
	{
	public:
		/**
		* Sphere的构造函数；
		* 输入一个double的数组赋值给圆心坐标，一个double值赋值给半径；
		* 默认维度为3维，不可更改；
		* 默认值半径为1.0，圆心坐标（0.0，0.0，0.0）可更改。
		*/
		explicit Sphere() :Geometry(3), radius(1.0) { center[0] = 0.0; center[1] = 0.0, center[2] = 0.0; }
		explicit Sphere(const double *c, const double &r) :Geometry(3)
		{
			SetCenter(c);
			SetRadius(r);
		}
		explicit Sphere(const double &c1, const double &c2, const double &c3, const double &r) :Geometry(3)
		{
			SetCenter(c1, c2, c3);
			SetRadius(r);
		}
		/**
		*返回类型名称. 通过name函数内部定义字符串实现一种类型一个名称。
		*/
		virtual ::std::string name()const;
		/**
		* 球心坐标设置内联函数；
		* 输入double数组值，赋值给圆心坐标的x,y,z值；
		* 无返回值。
		*/
		inline void SetCenter(
      const double *c
    ) 
    { 
      center[0] = c[0], center[1] = c[1], center[2] = c[2]; 
    }
		inline void SetCenter(
      const double &c1, 
      const double &c2, 
      const double &c3
    ) 
    { 
      center[0] = c1, center[1] = c2, center[2] = c3; 
    }
		/**
		* 半径设置内联函数;
		* 输入一个double值，赋值给半径;
		* 无返回值.
		*/
		inline void SetRadius(const double &r) { radius = r; }
		/**
		* 参数输入计算虚函数，计算该点是否在球内：
		* 实现输入一个点坐标数组，以及误差tolerance；
		* \param * 要计算的点坐标数组；
		* \param tolerance为所允许误差，默认1e-9；
		* return int值点在内部返回0，边界上返回1，外部返回2。
		*/
		virtual int CalRelativeLocation(const double *, const double &tolerance = 1e-9)const;
		/**
		* 输出图形中的球心坐标：
		* \return 一个大小为2的double数组，球心坐标数组，外部不可更改。
		*/
		virtual const double *GetCenter()const { return center; }
		/**
		* 输出图形中的半径：
		* \return double值，球心半径，外部不可更改。
		*/
		virtual const double GetRadius()const { return radius; }
	private:
		/**
		* \param center，球心点的坐标数组；
		* \param radius 半径；
		*/
		double center[3], radius;
		/**
		* \brief 虚函数输出当前对象的二进制文件数据；
		* 实现自己的二进制输出；
		* \param ostr 二进制文件输出流；
		* \param filename [in]输出的文件名；
		* \param path [in]输出的路径，默认为当前目录；
		* \return 返回bool变量表示输出状态，正常为true，异常为false。
		*/
		virtual bool OutputBinPrivate(
			std::ofstream &ostr,
			const std::string &filename,
			const std::string &path
		) const;
	
	};

	  /**
	  * 长方体，继承自Geometry，三维
		* 通过设定对角的两个顶点，构造一个长方体；
		* 私有成员包括两个对角顶点的坐标。
		*/
	class COMMONAPI Cube :public Geometry
	{
	public:
		/**
		* Cube的构造函数；
		* 输入两个double的数组赋值给两个对角顶点；
		* 默认维度为3维，不可更改；
		* 默认值两点初始坐标（0，0，0）和（1，2，3）可更改。
		*/
		explicit Cube() :Geometry(3)
		{
			vertex1[0] = 0; vertex2[1] = 0; vertex2[2] = 0;
			vertex2[0] = 1.0; vertex2[1] = 2.0; vertex2[2] = 3.0;
		}
		explicit Cube(const double *p1, const double *p2) :Geometry(3)
		{
			SetVertex(p1, p2);
		}
		/**
		*返回类型名称. 通过name函数内部定义字符串实现一种类型一个名称。
		*/
		virtual ::std::string name()const;
		/**
		* 对角端点坐标设置内联函数；
    * 两个顶点坐标必须严格满足第一个点的x坐标和y坐标以及z坐标都小于第二个点；
		* \param p1,p2 两个顶点的坐标数组；
		* \return 返回bool变量表示输出状态，正常为true，异常为false。
		*/
		virtual bool SetVertex(const double *p1, const double *p2);
		/**
		* 参数输入计算虚函数，计算该点是否在长方体内：
		* 实现输入一个点坐标数组，以及误差tolerance；
		* \param * 要计算的点坐标数组；
		* \param tolerance为所允许误差，默认1e-9；
		* return int值点在内部返回0，边界上返回1，外部返回2。
		*/
		virtual int CalRelativeLocation(const double *, const double &tolerance = 1e-9)const;
		/**
		* 输出图形中的1顶点坐标：
		* \return 返回double值，1顶点坐标，外部不可更改。
		*/
		inline const double *GetVertex1()const { return vertex1; }
		/**
		* 输出图形中的2顶点坐标：
		* \return 返回double值，2顶点坐标，外部不可更改。
		*/
		inline const double *GetVertex2()const { return vertex2; }
	private:
		double vertex1[3], vertex2[3];
		/**
		* \brief 虚函数输出当前对象的二进制文件数据；
		* 实现自己的二进制输出；
		* \param ostr 二进制文件输出流；
		* \param filename [in]输出的文件名；
		* \param path [in]输出的路径，默认为当前目录；
		* \return 返回bool变量表示输出状态，正常为true，异常为false。
		*/
		virtual bool OutputBinPrivate(
			std::ofstream &ostr,
			const std::string &filename,
			const std::string &path
		) const;
	};

	/**
		* 并集Geometry，继承自Geometry；
		* 私有数据成员包括一个Geometry的智能指针容器。
		*/
	class COMMONAPI UnionGeo : public Geometry
	{
	public:
		/**
		* UnionGeo的构造函数，通过两个Geometry的智能指针，生成一个UnionGeo对象。
		*/
		explicit UnionGeo(std::shared_ptr<Geometry> &g1, std::shared_ptr<Geometry> &g2);
		explicit UnionGeo(const int &dim) : Geometry(dim) {}
		explicit UnionGeo(std::vector<std::shared_ptr<Geometry>> geoVecIn) 
      : Geometry(geoVecIn[0]->GetDim()), geoVec(geoVecIn) {}

		/**
		*返回类型名称. 通过name函数内部定义字符串实现一种类型一个名称。
		*/
		virtual ::std::string name()const;
		/**
		* 参数输入计算虚函数，计算该点是否在两个Geometry的并集区域内：
		* 实现输入一个点坐标数组，以及误差tolerance；
		* \param * 要计算的点坐标数组；
		* \param tolerance为所允许误差，默认1e-10；
		* return int值点在内部返回0，边界上返回1，外部返回2。
		*/
		virtual int CalRelativeLocation(const double *, const double &tolerance = 1e-10) const;
		/**
		* Geometry加入到私有vector中的函数：
		* 实现输入Geometry的智能指针；
		* \param g1 要加入的Geometry；
		* return int值点在内部返回0，边界上返回1，外部返回2。
		*/
		bool add_geo(std::shared_ptr<Geometry> &g1);
	private:
		/**
		* UnionGeo的私有数据成员vector；
		* 存有一系列Geometry的智能指针。
		*/
		std::vector<std::shared_ptr<Geometry>> geoVec;
	};

	/**/
	class COMMONAPI UniversalGeo : virtual public Geometry
	{
	public:
		virtual ::std::string name()const { return "UniversalGeo"; }
		explicit UniversalGeo(const size_t &dim) : Geometry(dim) {}
		virtual int CalRelativeLocation(const double* pos, const double &tolerance = 1e-9) const
		{
			return geo::inside;
		}
	};

	/* 几何图形求补集类*/
	class COMMONAPI ComplementGeo :virtual public Geometry
	{
	public:
		virtual ::std::string name()const { return "ComplementGeo"; }
		/*
		 *构造函数
		 */
		explicit ComplementGeo(const std::shared_ptr<Geometry> &g1);
		/*
		* \brief 计算相对位置函数
		* \param pos [in] 空间点位置,储存在数组中
		* \param tolerance [in] 设置容差为1e-9，用来判断空间点与几何图形的补集的相对位置
		* \return 返回geo的枚举数, inside=0, boundary=1, outside=2
		*/
		virtual int CalRelativeLocation(const double *pos, const double &tolerance = 1e-9)const;
	private:
		std::shared_ptr<Geometry> geometry;
	};

	/* 几何图形求交集类*/
	class COMMONAPI IntersectionGeo : virtual public Geometry
	{
	public:
		virtual ::std::string name() const { return "IntersectionGeo"; }
		/*
		 *构造函数
		 */
		explicit IntersectionGeo(std::shared_ptr<Geometry> g1, std::shared_ptr<Geometry> g2);
    explicit IntersectionGeo(const int &dim) : Geometry(dim) {}
    explicit IntersectionGeo(std::vector<std::shared_ptr<Geometry>> geoVecIn)
      : Geometry(geoVecIn[0]->GetDim()), geoVec(geoVecIn) {}

		/*
		* \brief 计算相对位置函数
		* \param pos [in] 空间点位置,储存在数组中
		* \param tolerance [in] 设置容差为1e-9，用来判断空间点与几何图形的交集的相对位置
		* \return 返回geo的枚举数, inside=0, boundary=1, outside=2
		*/
		virtual int CalRelativeLocation(const double *pos, const double &tolerance = 1e-9)const;

		bool add_geo(std::shared_ptr<Geometry> &g1);
	private:
		std::vector<std::shared_ptr<Geometry>> geoVec;
		std::shared_ptr<Geometry> geo1;
		std::shared_ptr<Geometry> geo2;
	};

	/* 几何图形相减类*/
	class COMMONAPI SubtractGeo :virtual public Geometry
	{
	public:
		virtual ::std::string name()const { return "SubtractGeo"; }
		/*
		 *构造函数
		 */
		explicit SubtractGeo(std::shared_ptr<Geometry> g1, std::shared_ptr<Geometry> g2);
		/*
		* \brief 计算相对位置函数
		* \param pos [in] 空间点位置,储存在数组中
		* \param tolerance [in] 设置容差为1e-9，用来判断空间点与两个几何图形相减之后的相对位置
		* \return 返回geo的枚举数, inside=0, boundary=1, outside=2
		*/
		virtual int CalRelativeLocation(const double *pos, const double &tolerance = 1e-9)const;
	private:
		std::shared_ptr<Geometry> geo1;
		std::shared_ptr<Geometry> geo2;
	};

  class COMMONAPI Box : public Geometry
  {
  public:
    explicit Box(const size_t &dim = 1) : Geometry(dim), box(2 * dim) {}
    explicit Box(const std::vector<double> &boxIn) : 
      Geometry(boxIn.size() / 2), box(boxIn) {}
    Box(const Box &b) : Geometry(b), box(b.box) {}
    virtual ~Box() {}
    virtual int CalRelativeLocation(const double *, const double &tolerance) const;
    const double *GetBox() const { return box.data(); }
    double &operator[](const size_t &i) { return box[i]; }
    const double &operator[](const size_t &i) const { return box[i]; }
    std::string name() const { return "Box"; }
    bool check() const;
  public:
    static Box FindBoundingBoxOf(
      const double *nodes,
      const size_t dim,
      const size_t nodeNum
    );
    static Box CalIntersection(
      const Box &r1, 
      const Box &r2
    );
  private:
    // usage: box[i+dim*0]; box[i+dim*1], dim is the dimension of the box. 0,1 represent two ends in a certain direction.
    std::vector<double> box;
  };
}

#endif 
