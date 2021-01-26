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
#pragma once
#include "FYSimpleArray.h"
#include "FYSmartMemory.h"
#include "FYTypeWrap.h"
#include <vector>
using namespace std;

namespace PHSPACE
{

template < typename T, int N >
class FYFastSmartArrayIterator;

template < typename T, int N >
class FYSmartArray;

template < typename T, int N >
class FYSmartArray : public FYTypeWrap< FYSmartArray< T, N > >
{
public:
	typedef T        TValueType;
	typedef T      * TypeP1; //TypeP1 T *
	typedef TypeP1 * TypeP2; //TypeP2 T **
	typedef TypeP2 * TypeP3; //TypeP3 T ***
	typedef TypeP3 * TypeP4; //TypeP4 T ****

	typedef FYSmartArray< T, N >        TArray;
	typedef FYFastSmartArrayIterator< T, N > TIterator;
public:
	T    * data ;
	TypeP1 data1, datap1; //TypeP1 T *
	TypeP2 data2, datap2; //TypeP2 T **
	TypeP3 data3, datap3; //TypeP3 T ***
	TypeP4 data4, datap4; //TypeP4 T ****

	SimpleArray< int, N > lengthArray;
	SimpleArray< int, N > strideArray;
	SimpleArray< int, N > baseArray;
	int zeroOffSet;
	int numberOfElements;
	bool ownedData;
	bool ownedRootData;
public:
	~FYSmartArray()
	{
		this->DeAllocateAll();
	}

    FYSmartArray( const FYSmartArray< T, N > & fatherArray )
        : FYTypeWrap< FYSmartArray<T, N> >( fatherArray )
    {
		this->lengthArray = fatherArray.lengthArray;
		this->baseArray   = fatherArray.baseArray;
		this->strideArray = fatherArray.strideArray;
		this->data        = const_cast< T * >( fatherArray.GetData() );
		this->zeroOffSet  = fatherArray.zeroOffSet;
		this->ownedData   = false;
    }

	FYSmartArray( const Range & r0 )
	{
		vector< Range > rangeArray;
		rangeArray.push_back( r0 );
		this->AllocateBasicData( rangeArray );
		this->ZeroAllPointer();

		FYArrayPointer< T > arrayPointer;
		PHSPACE::AllocateRange1< T >( r0, arrayPointer );

		data1  = arrayPointer.data1;
		datap1 = arrayPointer.data1New;

		data  = data1;
	}

	FYSmartArray( const Range & r0, const Range & r1 )
	{
		vector< Range > rangeArray;
		rangeArray.push_back( r0 );
		rangeArray.push_back( r1 );
		this->AllocateBasicData( rangeArray );
		this->ZeroAllPointer();

		FYArrayPointer< T > arrayPointer;
		PHSPACE::AllocateRange2< T >( r0, r1, arrayPointer );

		data1  = arrayPointer.data1;
		data2  = arrayPointer.data2;

		datap1 = arrayPointer.data1New;
		datap2 = arrayPointer.data2New;

		data  = data1;
	}

	FYSmartArray( const Range & r0, const Range & r1, const Range & r2 )
	{
		vector< Range > rangeArray;
		rangeArray.push_back( r0 );
		rangeArray.push_back( r1 );
		rangeArray.push_back( r2 );
		this->AllocateBasicData( rangeArray );
		this->ZeroAllPointer();

		FYArrayPointer< T > arrayPointer;
		PHSPACE::AllocateRange3< T >( r0, r1, r2, arrayPointer );

		data1  = arrayPointer.data1;
		data2  = arrayPointer.data2;
		data3  = arrayPointer.data3;

		datap1 = arrayPointer.data1New;
		datap2 = arrayPointer.data2New;
		datap3 = arrayPointer.data3New;

		data = data1;
	}

	FYSmartArray( const Range & r0, const Range & r1, const Range & r2, const Range & r3 )
	{
		vector< Range > rangeArray;
		rangeArray.push_back( r0 );
		rangeArray.push_back( r1 );
		rangeArray.push_back( r2 );
		rangeArray.push_back( r3 );
		this->AllocateBasicData( rangeArray );
		this->ZeroAllPointer();

		FYArrayPointer< T > arrayPointer;
		PHSPACE::AllocateRange4< T >( r0, r1, r2, r3, arrayPointer );

		data1  = arrayPointer.data1;
		data2  = arrayPointer.data2;
		data3  = arrayPointer.data3;
		data4  = arrayPointer.data4;

		datap1 = arrayPointer.data1New;
		datap2 = arrayPointer.data2New;
		datap3 = arrayPointer.data3New;
		datap4 = arrayPointer.data4New;

		data  = data1;
	}

	FYSmartArray( T * dataPointer, const Range & r0 )
	{
		vector< Range > rangeArray;
		rangeArray.push_back( r0 );
		this->AllocateBasicData( rangeArray );
		this->ZeroAllPointer();
		this->ownedRootData = false;

		FYArrayPointer< T > arrayPointer;
		PHSPACE::AllocateRangeMap1< T >( r0, dataPointer, arrayPointer );

		data1  = arrayPointer.data1;
		datap1 = arrayPointer.data1New;

		data  = data1;
	}

	FYSmartArray( T * dataPointer, const Range & r0, const Range & r1 )
	{
		vector< Range > rangeArray;
		rangeArray.push_back( r0 );
		rangeArray.push_back( r1 );
		this->AllocateBasicData( rangeArray );
		this->ZeroAllPointer();
		this->ownedRootData = false;

		FYArrayPointer< T > arrayPointer;
		PHSPACE::AllocateRangeMap2< T >( r0, r1, dataPointer, arrayPointer );

		data1  = arrayPointer.data1;
		data2  = arrayPointer.data2;

		datap1 = arrayPointer.data1New;
		datap2 = arrayPointer.data2New;

		data  = data1;
	}

	FYSmartArray( T * dataPointer, const Range & r0, const Range & r1, const Range & r2 )
	{
		vector< Range > rangeArray;
		rangeArray.push_back( r0 );
		rangeArray.push_back( r1 );
		rangeArray.push_back( r2 );
		this->AllocateBasicData( rangeArray );
		this->ZeroAllPointer();
		this->ownedRootData = false;

		FYArrayPointer< T > arrayPointer;
		PHSPACE::AllocateRangeMap3< T >( r0, r1, r2, dataPointer, arrayPointer );

		data1  = arrayPointer.data1;
		data2  = arrayPointer.data2;
		data3  = arrayPointer.data3;

		datap1 = arrayPointer.data1New;
		datap2 = arrayPointer.data2New;
		datap3 = arrayPointer.data3New;

		data = data1;
	}

	FYSmartArray( T * dataPointer, const Range & r0, const Range & r1, const Range & r2, const Range & r3 )
	{
		vector< Range > rangeArray;
		rangeArray.push_back( r0 );
		rangeArray.push_back( r1 );
		rangeArray.push_back( r2 );
		rangeArray.push_back( r3 );
		this->AllocateBasicData( rangeArray );
		this->ZeroAllPointer();
		this->ownedRootData = false;

		FYArrayPointer< T > arrayPointer;
		PHSPACE::AllocateRangeMap4< T >( r0, r1, r2, r3, dataPointer, arrayPointer );

		data1  = arrayPointer.data1;
		data2  = arrayPointer.data2;
		data3  = arrayPointer.data3;
		data4  = arrayPointer.data4;

		datap1 = arrayPointer.data1New;
		datap2 = arrayPointer.data2New;
		datap3 = arrayPointer.data3New;
		datap4 = arrayPointer.data4New;

		data = data1;
	}
public:
    TIterator GetFastIterator() const
    { 
		return TIterator( * this );
	}

    TArray & NonConst() const
    { 
		return const_cast< TArray & >( * this );
	}
public:
	FYSmartArray( TArray & fatherArray, const Range & r0 );
	FYSmartArray( TArray & fatherArray, const Range & r0, const Range & r1 );
	FYSmartArray( TArray & fatherArray, const Range & r0, const Range & r1, const Range & r2 );
	FYSmartArray( TArray & fatherArray, const Range & r0, const Range & r1, const Range & r2, const Range & r3 );

    TArray operator()( const Range & r0 ) const
    {
        return TArray( this->NonConst(), r0 );
    }

    TArray operator()( const Range & r0, const Range & r1 ) const
    {
        return TArray( this->NonConst(), r0, r1 );
    }

    TArray operator()( const Range & r0, const Range & r1, const Range & r2 ) const
    {
        return TArray( this->NonConst(), r0, r1, r2 );
    }

    TArray operator()( const Range & r0, const Range & r1, const Range & r2, const Range & r3 ) const
    {
        return TArray( this->NonConst(), r0, r1, r2, r3 );
    }

	void SliceArray ( int i, const Range & r );
	void SliceArray ( vector< Range > & rangeArray );
public:

	SimpleArray< int, N > & GetBaseArray(){ return baseArray; }

	void AllocateBasicData( vector< Range > & rangeArray );
	void SetBasicRange( vector< Range > & rangeArray );
	void ComputeStrides();
	void ComputeZeroOffset();
	void ZeroAllPointer();
	void DeAllocateAll();
	void ComputeNumberOfElements();
public:
    const T & operator()( int i0 ) const
    { 
        return datap1[ i0 ];
    }

    T & operator()( int i0 )
    {
		return datap1[ i0 ];
	}
	
	const T & operator()( int i0, int i1 ) const
    { 
		return datap2[ i1 ][ i0 ];
    }

    T & operator()( int i0, int i1 )
    {
		return datap2[ i1 ][ i0 ];
    }

    const T & operator()( int i0, int i1, int i2 ) const
    {
        return datap3[ i2 ][ i1 ][ i0 ];
    }

    T & operator()( int i0, int i1, int i2 ) 
    {
		return datap3[ i2 ][ i1 ][ i0 ];
    }

    const T & operator()( int i0, int i1, int i2, int i3 ) const
    {
		return datap4[ i3 ][ i2 ][ i1 ][ i0 ];
    }

    T & operator()( const int & i0, const int & i1, const int & i2, const int & i3 ) 
    {
		return datap4[ i3 ][ i2 ][ i1 ][ i0 ];
    }
public:
	T * GetData() { return data; }

    const T * GetData() const
    {
		return data;
	}

	T * GetDataPointer( int ist, int jst, int kst )
	{
		vector< int > startIndex;

		startIndex.push_back( ist );
		startIndex.push_back( jst );
		startIndex.push_back( kst );

		int offset = 0;

		for ( int i = 0; i < N; ++ i )
		{
			offset += strideArray[ i ] * ( startIndex[ i ] - baseArray[ i ] );
		}

		return data + offset; 
	}

	void Fill( const T & value )
	{
		// assign value to all elements
		for ( int i = 0; i < numberOfElements; ++ i )
		{
			data[ i ] = value;
		}
	}

	FYSmartArray< T, N > & operator = ( const T & value );
	FYSmartArray< T, N > & operator = ( const FYSmartArray< T, N > & rhs );

	template < typename TExpression, typename TUpdate >
	inline FYSmartArray< T, N > & Evaluate( TExpression tExpression, TUpdate );

	template < typename TExpression, typename TUpdate >
	inline FYSmartArray< T, N > & EvaluateWithStackTraversal1( TExpression tExpression, TUpdate );

	template < typename TExpression, typename TUpdate >
	inline FYSmartArray< T, N > & EvaluateWithStackTraversalN( TExpression tExpression, TUpdate );

	template < typename TExpression >
	inline FYSmartArray< T, N > & operator = ( const FYTypeWrap< TExpression > & tExpression );

    template < typename TExpression >
	inline TArray & operator += ( const TExpression & tExpression );

    template < typename TExpression >
	inline TArray & operator -= ( const TExpression & tExpression );

    template < typename TExpression >
	inline TArray & operator *= ( const TExpression & tExpression );

    template < typename TExpression >
	inline TArray & operator /= ( const TExpression & tExpression );

    int GetStride( int i ) const
    { 
		return strideArray[ i ];
	}
};

}

#include "FYFastSmartArrayIterator.h"
//#include "FYCfdMath.h"
#include "FYArrayExpression.h"
#include "FYArrayExpressionWrap.h"
#include "FYArrayExpressionTemplates.h"
#include "FYSmartArray.hxx"
