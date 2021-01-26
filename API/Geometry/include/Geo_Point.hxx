//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//                                                                             +
//       Hybrid Platform for Engineering and Research of Flows (HyperFlow)     +
//                        面向流体工程与研究的混合CFD平台                      +
//                        （C）Zhang Laiping and He Xin                        +
//                    State Key Laboratory of Aerodynamics                     +
//                    Computational Aerodynamics Institute                     +
//            China Aerodynamics Research and Development Center               +
//                                 Since 2007                                  +
//                                                                             +
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// 文件头部说明：                                                              +
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// FILE:       XXX.cpp (XXX.h)                                                 +
// AUTHOR(S):  He Xin                                                          +
// PURPOSE:    simulation object                                               +
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// 函数头部说明：                                                              +
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:   例如：Simulation::Simulation()                                  +
// AUTHOR(S):  He Xin                                                          +
// PURPOSE:    具体说明本模块的功能                                            +
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Geo_Point.h"

#include <iostream>
using namespace std;

namespace PHSPACE
{
template < typename T >
Geo_Point<T>::Geo_Point()
{
}

template < typename T >
Geo_Point<T>::Geo_Point(const T &x, const T &y, const T &z, int id)
{
	this->x  = x;
	this->y  = y;
	this->z  = z;
	this->id = id;
    tolerance = 0.0;
}

template < typename T >
Geo_Point<T>::Geo_Point(const point_type & rhs)
{
	this->x  = rhs.x;
	this->y  = rhs.y;
	this->z  = rhs.z;
	this->id = rhs.id;
    tolerance = 0.0;
}

template < typename T >
Geo_Point<T>::~Geo_Point(void)
{
}

template < typename T >
Geo_Point<T> & Geo_Point<T>::operator = (const Geo_Point<T> & rhs)
{
	if ( this == &rhs ) return *this;

	this->x  = rhs.x;
	this->y  = rhs.y;
	this->z  = rhs.z;
	this->id = rhs.id;

	return *this;
}

template < typename T >
bool Geo_Point<T>::operator < (const Geo_Point<T> & rhs) const
{
	T dx = x - rhs.x;
	T dy = y - rhs.y;
	T dz = z - rhs.z;

	T diff = 1.0e-8;

	if ( ABS( dx ) > diff ) return x < rhs.x;
	if ( ABS( dy ) > diff ) return y < rhs.y;
	if ( ABS( dz ) > diff ) return z < rhs.z;

	return false;
}

template < typename T >
int Geo_Point<T>::hash_key() const
{
	RDouble        dbl[3];
	int           *intg = (int *) dbl;

	dbl[0] = x;
	dbl[1] = y;
	dbl[2] = z;

	int           hashnum = 0;
	for (int i = 0; i < 3; i++) 
	{
		hashnum += (intg[i] >> i);
	}

	return (int) ABS(hashnum);   
}

}
