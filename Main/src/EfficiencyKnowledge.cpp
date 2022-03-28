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
#include <iostream>
#include "FYEfficiencyKnowledge.h"
#include "TK_Time.h"
using namespace std;

blitz::_bz_fortranTag blitz::fortranArray;

namespace PHSPACE
{

FYEfficiencyKnowledge::FYEfficiencyKnowledge()
{
	value1 = 101.123456789;
	value2 = 99.987654321;
	value3 = 3.123456789;
	value4 = 2.987654321;

	narray_for_test = 5;

	arrayef = new ArrayEfficiency * [narray_for_test];
	arrayef[0] = new FYArrayEfficiency();
	arrayef[1] = new BlitzArrayEfficiency();
	arrayef[2] = new CPPArrayEfficiency();
	arrayef[3] = new MatrixArrayEfficiency();
	arrayef[4] = new BoostArrayEfficiency();
}

FYEfficiencyKnowledge::~FYEfficiencyKnowledge()
{
	for ( int i = 0; i < narray_for_test; ++ i )
	{
		delete arrayef[i];
	}

	delete [] arrayef;
}

bool FYEfficiencyKnowledge::Run()
{
	const int n = 3;
	blitz::Array<int,2> A(n,n),B(n,n),C(n,n),D(n,n);
	A = 1;
	B = 2;
	C = 3;
    //D = 2 * A + B + C;
	D = 2 * A + 3.0 * B - C;
	D = 2 * A;
	D = 1 + A;
	D = A + 1;
	blitz::Range r1(0,2), r2(0,2);
	//D(r1,r2) = A(r1,r2) + B(r1,r2);
	D(r1,r2) = A(r1,r2) + A(r1,r2) + A(r1,r2);
	D(r1,r2) = abs(A(r1,r2));
	D(r1,r2) += -A(r1,r2);
	cout << D << "\n";

	//blitz::Array<int,1> A1(A(blitz::Range(0,n-1),0));
	blitz::Array<int,1> A1(n);
	A1 = A(blitz::Range(0,n-1),0);
	A(blitz::Range(0,n-1),0) = A1(blitz::Range(0,n-1));

	blitz::Array<int,3> AAA(n,n,n);
	AAA(blitz::Range(0,n-1),0,blitz::Range(0,n-1)) = A;

	blitz::Array<int,2> A2 = AAA(blitz::Range::all(), 2, blitz::Range::all());

	//FYArrayExprConstant<int> a(1),b(2);
	//FYArrayExprBinaryOp<int,int> c = 1 + 2;
	//FYArrayExpr<int> c = 1 + 2;
	FYArray<int, 1> c(Range(1,3));
	FYArray<int, 1> d(Range(1,3));

	blitz::Array<int,1> ABC(10);

	d = 10;
	c = 1 + d;
	//1 + d;

	FYArray<RFloat, 2> cc(Range(1,3),Range(1,2));
	FYArray<RFloat, 2> dd(Range(1,3),Range(1,2));
	FYArray<RFloat, 2> ee(Range(1,3),Range(1,2));
	FYArray<RFloat, 2> ff(Range(1,3),Range(1,2));
	ff = 2;
	ee = 10;
	cc = 0;
	//dd = 0;
	for ( int i = 1; i <= 3; ++ i )
	{
		for ( int j = 1; j <= 2; ++ j )
		{
			cc(i,j) = i * 100 + j;
		}
	}

	Range ri(1,2), rj(1,2);
	dd(ri,rj) = cc(ri,rj) + cc(ri,rj) + cc(ri,rj) + ee(ri,rj) + ff(ri,rj);
	//dd(ri,rj) = cc(ri,rj);
	dd(ri,rj) = MinMod(ee(ri,rj),ff(ri,rj));
	dd(ri,rj) = -cc(ri,rj);
	dd(ri,rj) = 0;
	dd(ri,rj) += -cc(ri,rj);

	Range I(1,2);
	Range J=I-1-2;
	Range K = J - 1;

	FYArray<double,2> fa(Range(1,4),Range(1,4),fortranArray);
	FYArray<double,2> fb(Range(1,4),Range(1,4),fortranArray);
	FYArray<double,2> fc(Range(1,4),Range(1,4),fortranArray);

	FYArray<double,3> fd(Range(1,4),Range(1,4),Range(1,4),fortranArray);

	Range all = Range::all();

	fd(Range(1,4),1,Range(1,4)) = 1;
	fd(Range(1,4),2,Range(1,4)) = 2;
	fd(Range(1,4),3,Range(1,4)) = 3;
	fd(Range(1,4),4,Range(1,4)) = 4;

	fd(all,1,all) = 100;
	fd(all,2,all) = 200;
	fd(all,3,all) = -300;
	fd(all,4,all) = 400;

	for ( int j = 1; j <= 4; ++ j )
	{
		for ( int i = 1; i <= 4; ++ i )
		{
			for ( int k = 1; k <= 4; ++ k )
			{
				cout << "fd["<< i << "]["<< j << "]["<< k << "]" << fd(i,j,k) << "\n";
			}
		}
	}
	cout << "\n";

	Range II(2,3),JJ(2,3);

	for ( int i = 1; i <= 4; ++ i )
	{
		for ( int j = 1; j <= 4; ++ j )
		{
			//if ( i <= 2 )
			//{
			//	fa(i,j) = -1;
			//}
			//else
			//{
			//	fa(i,j) = 1;
			//}
			fa(i,j) = ( i - 0 ) * 100 + j;
		}
	}

	for ( int i = 1; i <= 4; ++ i )
	{
		for ( int j = 1; j <= 4; ++ j )
		{
			for ( int k = 1; k <= 4; ++ k )
			{
				//if ( i <= 2 )
				//{
				//	fa(i,j) = -1;
				//}
				//else
				//{
				//	fa(i,j) = 1;
				//}

				fd(i,j,k) = (k-0) * 1000 + ( i - 0 ) * 100 + j;
			}
		}
	}

	//fb = fa;

	fb(II,JJ) = fa(II,JJ);
	//fc(II,JJ) = fa(II+1,JJ+1);
	//fc(II,JJ) = fd(II,JJ,Range(1,1));

	for ( int i = 1; i <= 4; ++ i )
	{
		for ( int j = 1; j <= 4; ++ j )
		{
			cout << "fa["<< i << "]["<< j << "]" << fa(i,j) << "\n";
		}
	}
	cout << "\n";

	for ( int i = 1; i <= 4; ++ i )
	{
		for ( int j = 1; j <= 4; ++ j )
		{
			cout << "fb["<< i << "]["<< j << "]" << fb(i,j) << "\n";
		}
	}
	cout << "\n";

	for ( int i = 1; i <= 4; ++ i )
	{
		for ( int j = 1; j <= 4; ++ j )
		{
			cout << "fc["<< i << "]["<< j << "]" << fc(i,j) << "\n";
		}
	}

	cout << "\n";

	blitz::Array<double,2> FA(blitz::Range(1,4),blitz::Range(1,4));
	blitz::Array<double,2> FB(blitz::Range(1,4),blitz::Range(1,4));
	blitz::Array<double,2> FC(blitz::Range(1,4),blitz::Range(1,4));
	blitz::Array<double,3> FD(blitz::Range(1,4),blitz::Range(1,4),blitz::Range(1,4));

	blitz::Range BII(2,3),BJJ(2,3);

	for ( int i = 1; i <= 4; ++ i )
	{
		for ( int j = 1; j <= 4; ++ j )
		{
			if ( i <= 2 )
			{
				FA(i,j) = -1;
			}
			else
			{
				FA(i,j) = 1;
			}
		}
	}

	for ( int i = 1; i <= 4; ++ i )
	{
		for ( int j = 1; j <= 4; ++ j )
		{
			for ( int k = 1; k <= 4; ++ k )
			{
				FD(i,j,k) = (k-0) * 1000 + ( i - 0 ) * 100 + j;
			}
		}
	}

	FB(BII,BJJ) = FA(BII  ,BJJ  );
	//FC(BII,BJJ) = FA(BII+1,BJJ+1);
	//FC(BII,BJJ) = FD(BII,BJJ,blitz::Range(1,1));


	for ( int i = 1; i <= 4; ++ i )
	{
		for ( int j = 1; j <= 4; ++ j )
		{
			cout << "FB["<< i << "]["<< j << "]" << FB(i,j) << "\n";
		}
	}
	cout << "\n";

	for ( int i = 1; i <= 4; ++ i )
	{
		for ( int j = 1; j <= 4; ++ j )
		{
			cout << "FC["<< i << "]["<< j << "]" << FC(i,j) << "\n";
		}
	}

	cout << "\n";

	//dd = 1 + cc;
	//dd = cc + 1;
	////dd = cc - 1;
	////dd = cc;
	//dd = cc * 2 + cc / 1.5;
	////dd = cc / 2;

	RFloat itest;

	for ( int i = 1; i <= 3; ++ i )
	{
		for ( int j = 1; j <= 2; ++ j )
		{
			itest = dd(i,j);
			cout << "dd["<< i << "]["<< j << "]" << dd(i,j) << "\n";
		}
	}

	int kkk = Add<int,int>::apply(1,2);


	//FYETBase<int> aa;
	//int cc = aa.unwrap();

	//blitz::Array<double, 3> aaaa(10,15,20,blitz::fortranArray);
	//aaaa = 0.0;
	//
	//blitz::Range I(1,10);
	//blitz::Range J(1,20);
	//blitz::Range K(1,30);
	//
	//blitzArrayDouble3D aaaa(I,J,K,blitz::fortranArray);
	//
	//
	//TestBasicOP();

	TestArrayOrder();

	return true;
}

void FYEfficiencyKnowledge::TestBasicOP()
{
	int iter_max = 1000;

	for ( int i = 0; i < narray_for_test; ++ i )
	{
		arrayef[i]->InitArray(iter_max);
	}

	for ( int i = 0; i < narray_for_test; ++ i )
	{
		arrayef[i]->AddConst(iter_max, value1);
	}

}

const int ni = 100;
const int nj = 100;
const int nk = 100;

double a[ni][nj][nk];

void FYEfficiencyKnowledge::TestArrayOrder()
{
	int nm = 5;
	Range I(0,ni-1);
	Range J(0,nj-1);
	Range K(0,nk-1);
	Range M(0,nm-1);

	int inum =  (I.last()-I.first())/I.stride();
	int jnum =  (J.last()-J.first())/J.stride();
	int knum =  (K.last()-K.first())/K.stride();

	RFloat4D a1(I,J,K,M,fortranArray);
	RFloat4D a2(M,I,J,K,fortranArray);

	int iter_max = 1000;

	int mst = M.first();
	int med = M.last();

	RFloat * prim = new RFloat [nm];

	for ( int m = 0; m < nm; ++ m )
	{
		prim[m] = m + 0.123456789;
	}

	int ist = I.first();
	int jst = J.first();
	int kst = K.first();

	int ied = I.last();
	int jed = J.last();
	int ked = K.last();

	TimeSpan *simu_ts = new TimeSpan();
	using namespace PHSPACE;

	for ( int it = 0; it < iter_max; ++ it )
	{
		for ( int m = mst; m <= med; ++ m )
		{
			for ( int k = kst; k <= ked; ++ k )
			{
				for ( int j = jst; j <= jed; ++ j )
				{
					for ( int i = ist; i <= ied; ++ i )
					{
						a1(i,j,k,m) = prim[m];
					}
				}
			}
		}
	}
	simu_ts->ShowTimeSpan("TestArrayOrder 1 ");

	simu_ts->ResetTime();

	for ( int it = 0; it < iter_max; ++ it )
	{
		for ( int k = kst; k <= ked; ++ k )
		{
			for ( int j = jst; j <= jed; ++ j )
			{
				for ( int i = ist; i <= ied; ++ i )
				{
					for ( int m = mst; m <= med; ++ m )
					{
						a1(i,j,k,m) = prim[m];
					}
				}
			}
		}
	}
	simu_ts->ShowTimeSpan("TestArrayOrder 2 ");

	simu_ts->ResetTime();

	for ( int it = 0; it < iter_max; ++ it )
	{
		for ( int m = mst; m <= med; ++ m )
		{
			a1(I,J,K,m) = prim[m];
		}
	}
	simu_ts->ShowTimeSpan("TestArrayOrder 3 ");
	delete simu_ts;

	delete [] prim;

}


void BoostArrayEfficiency::InitArray(int iter_max)
{
    typedef boost::multi_array<double, 3> array_type;
    typedef array_type::index index;
    array_type fa(boost::extents[ni][nj][nk]);

	TimeSpan *simu_ts = new TimeSpan();
	using namespace PHSPACE;

	for ( int i = 0; i < iter_max; ++ i )
	{
		SimpleAssign(fa);
		//SimpleAssign2(fa);
	}
	simu_ts->ShowTimeSpan("InitBoostArray ");
	delete simu_ts;
}

void BoostArrayEfficiency::AddConst(int iter_max, double &value)
{
    typedef boost::multi_array<double, 3> array_type;
    typedef array_type::index index;
    array_type fa(boost::extents[ni][nj][nk]);

	SimpleAssign(fa);

	TimeSpan *simu_ts = new TimeSpan();
	using namespace PHSPACE;

	for ( int i = 0; i < iter_max; ++ i )
	{
		AddConst(fa, value);
	}
	simu_ts->ShowTimeSpan("BoostArrayAddConst ");
	delete simu_ts;
}

void BoostArrayEfficiency::SimpleAssign(boostArrayDouble3D &a)
{
    typedef boostArrayDouble3D::index index;

    // Assign values to the elements
	for ( index k = 0; k < nk; ++k )
	{
		for ( index j = 0; j < nj; ++j )
		{
			for ( index i = 0; i < ni; ++ i )
			{
				a[i][j][k] = k;
			}
		}
	}
}

void BoostArrayEfficiency::SimpleAssign2(boostArrayDouble3D &a)
{
    typedef boostArrayDouble3D::index index;

    // Assign values to the elements
	for ( index k = 0; k < nk; ++k )
	{
		for ( index j = 0; j < nj; ++j )
		{
			for ( index i = 0; i < ni; ++ i )
			{
				a[i][j][k] = k;
			}
		}
	}
}

void BoostArrayEfficiency::AddConst(boostArrayDouble3D &a, double &value)
{
    typedef boostArrayDouble3D::index index;

    // Assign values to the elements
    for ( index i = 0; i < ni; ++ i )
	{
		for ( index j = 0; j < nj; ++j )
		{
			for ( index k = 0; k < nk; ++k )
			{
				a[i][j][k] += value;
			}
		}
	}
}

void BlitzArrayEfficiency::InitArray(int iter_max)
{
	blitz::Range I(0,ni-1);
	blitz::Range J(0,nj-1);
	blitz::Range K(0,nk-1);

	blitzArrayDouble3D fa(I,J,K,blitz::fortranArray);

	TimeSpan *simu_ts = new TimeSpan();
	using namespace PHSPACE;

	for ( int i = 0; i < iter_max; ++ i )
	{
		//SimpleAssign(fa);
		SimpleAssign2(fa);
	}
	simu_ts->ShowTimeSpan("InitBlitzArray ");
	delete simu_ts;
}

void BlitzArrayEfficiency::AddConst(int iter_max, double &value)
{
	blitz::Range I(0,ni-1);
	blitz::Range J(0,nj-1);
	blitz::Range K(0,nk-1);

	blitzArrayDouble3D fa(I,J,K,blitz::fortranArray);

	SimpleAssign2(fa);

	TimeSpan *simu_ts = new TimeSpan();
	using namespace PHSPACE;

	for ( int i = 0; i < iter_max; ++ i )
	{
		AddConst(fa, value);
	}
	simu_ts->ShowTimeSpan("BlitzArrayAddConst ");
	delete simu_ts;
}

void BlitzArrayEfficiency::SimpleAssign(blitzArrayDouble3D &a)
{
	for ( int i = 0; i < ni; ++ i )
	{
		for ( int j = 0; j < nj; ++ j )
		{
			for ( int k = 0; k < nk; ++ k )
			{
				a(i,j,k) = k;
			}
		}
	}
}

void BlitzArrayEfficiency::SimpleAssign2(blitzArrayDouble3D &a)
{
	for ( int k = 0; k < nk; ++ k )
	{
		for ( int j = 0; j < nj; ++ j )
		{
			for ( int i = 0; i < ni; ++ i )
			{
				a(i,j,k) = k;
			}
		}
	}
}

void BlitzArrayEfficiency::AddConst(blitzArrayDouble3D &a, double &value)
{
	for ( int k = 0; k < nk; ++ k )
	{
		for ( int j = 0; j < nj; ++ j )
		{
			for ( int i = 0; i < ni; ++ i )
			{
				a(i,j,k) += value;
			}
		}
	}
}

void FYArrayEfficiency::InitArray(int iter_max)
{
	Range I(0,ni-1);
	Range J(0,nj-1);
	Range K(0,nk-1);

	RFloat3D fa(I,J,K,fortranArray);

	TimeSpan *simu_ts = new TimeSpan();
	using namespace PHSPACE;

	for ( int i = 0; i < iter_max; ++ i )
	{
		//SimpleAssign(fa);
		SimpleAssign2(fa);
	}
	simu_ts->ShowTimeSpan("InitFYArray ");
	delete simu_ts;
}

void FYArrayEfficiency::AddConst(int iter_max, double &value)
{
	Range I(0,ni-1);
	Range J(0,nj-1);
	Range K(0,nk-1);
	RFloat3D fa(I,J,K,fortranArray);
	SimpleAssign2(fa);

	TimeSpan *simu_ts = new TimeSpan();
	using namespace PHSPACE;

	for ( int i = 0; i < iter_max; ++ i )
	{
		AddConst(fa, value);
	}
	simu_ts->ShowTimeSpan("FYArrayAddConst ");
	delete simu_ts;
}

void FYArrayEfficiency::SimpleAssign(RFloat3D &a)
{
	for ( int i = 0; i < ni; ++ i )
	{
		for ( int j = 0; j < nj; ++ j )
		{
			for ( int k = 0; k < nk; ++ k )
			{
				a(i,j,k) = k;
			}
		}
	}
}

void FYArrayEfficiency::SimpleAssign2(RFloat3D &a)
{
	for ( int k = 0; k < nk; ++ k )
	{
		for ( int j = 0; j < nj; ++ j )
		{
			for ( int i = 0; i < ni; ++ i )
			{
				a(i,j,k) = k;
			}
		}
	}
}

void FYArrayEfficiency::AddConst(RFloat3D &a, double &value)
{
	for ( int k = 0; k < nk; ++ k )
	{
		for ( int j = 0; j < nj; ++ j )
		{
			for ( int i = 0; i < ni; ++ i )
			{
				a(i,j,k) += value;
			}
		}
	}
}

void CPPArrayEfficiency::InitArray(int iter_max)
{
	TimeSpan *simu_ts = new TimeSpan();
	for ( int i = 0; i < iter_max; ++ i )
	{
		SimpleAssign();
		//SimpleAssign2();
	}
	simu_ts->ShowTimeSpan("InitCPPArray ");
	delete simu_ts;
}

void CPPArrayEfficiency::AddConst(int iter_max, double &value)
{
	SimpleAssign();

	TimeSpan *simu_ts = new TimeSpan();
	for ( int i = 0; i < iter_max; ++ i )
	{
		AddConst(value);
	}
	simu_ts->ShowTimeSpan("CPPArrayAddConst ");
	delete simu_ts;
}

void CPPArrayEfficiency::SimpleAssign()
{
	for ( int i = 0; i < ni; ++ i )
	{
		for ( int j = 0; j < nj; ++ j )
		{
			for ( int k = 0; k < nk; ++ k )
			{
				a[i][j][k] = k;
			}
		}
	}
}

void CPPArrayEfficiency::SimpleAssign2()
{
	for ( int k = 0; k < nk; ++ k )
	{
		for ( int j = 0; j < nj; ++ j )
		{
			for ( int i = 0; i < ni; ++ i )
			{
				a[i][j][k] = k;
			}
		}
	}
}

void CPPArrayEfficiency::AddConst(double &value)
{
	for ( int i = 0; i < ni; ++ i )
	{
		for ( int j = 0; j < nj; ++ j )
		{
			for ( int k = 0; k < nk; ++ k )
			{
				a[i][j][k] += value;
			}
		}
	}
}


void MatrixArrayEfficiency::InitArray(int iter_max)
{
	double *** fa = NewPointer3<double>(ni,nj,nk);
	TimeSpan *simu_ts = new TimeSpan();
	using namespace PHSPACE;

	for ( int i = 0; i < iter_max; ++ i )
	{
		SimpleAssign(fa);
		//SimpleAssign2(fa);
	}
	simu_ts->ShowTimeSpan("InitMatrixArray ");
	delete simu_ts;

	DelPointer3(fa);
}

void MatrixArrayEfficiency::AddConst(int iter_max, double &value)
{
	double *** fa = NewPointer3<double>(ni,nj,nk);
	SimpleAssign(fa);

	TimeSpan *simu_ts = new TimeSpan();
	using namespace PHSPACE;

	for ( int i = 0; i < iter_max; ++ i )
	{
		AddConst(fa, value);
	}
	simu_ts->ShowTimeSpan("MatrixArrayAddConst ");
	delete simu_ts;

	DelPointer3(fa);
}

void MatrixArrayEfficiency::SimpleAssign(double ***a)
{
	for ( int i = 0; i < ni; ++ i )
	{
		for ( int j = 0; j < nj; ++ j )
		{
			for ( int k = 0; k < nk; ++ k )
			{
				a[i][j][k] = k;
			}
		}
	}
}

void MatrixArrayEfficiency::SimpleAssign2(double ***a)
{
	for ( int k = 0; k < nk; ++ k )
	{
		for ( int j = 0; j < nj; ++ j )
		{
			for ( int i = 0; i < ni; ++ i )
			{
				a[i][j][k] = k;
			}
		}
	}
}

void MatrixArrayEfficiency::AddConst(double ***a, double &value)
{
	for ( int i = 0; i < ni; ++ i )
	{
		for ( int j = 0; j < nj; ++ j )
		{
			for ( int k = 0; k < nk; ++ k )
			{
				a[i][j][k] += value;
			}
		}
	}
}


}
