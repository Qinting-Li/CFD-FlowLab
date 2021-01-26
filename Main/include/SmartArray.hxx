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

namespace PHSPACE
{

template< typename X, typename Y >
class FYUpdateNewVersion
{
public:
    static inline void Update( X & x, Y y )
    {
		x = ( X )y; 
	}
};

#define FY_DEFINE_UPDATER_CLASS( FYClassName, Operator )    \
template< typename X, typename Y >                          \
class FYClassName                                           \
{                                                           \
public:                                                     \
    static inline void Update( X& x, Y y )                  \
    {                                                       \
		x Operator ( X )y;                                  \
	}                                                       \
};

template< typename T, int N >
void FYSmartArray< T, N >::AllocateBasicData( vector< Range > & rangeArray )
{
	this->SetBasicRange( rangeArray );
	this->ComputeStrides();
	this->ComputeNumberOfElements();
}

template< typename T, int N >
void FYSmartArray< T, N >::SetBasicRange( vector< Range > & rangeArray )
{
	for ( int i = 0; i < N; ++ i )
	{
		Range * range = & rangeArray[ i ];

		lengthArray[ i ] = range->length();
		baseArray  [ i ] = range->first();
	}
}

template< typename T, int N >
void FYSmartArray< T, N >::ComputeStrides()
{
	int stride = 1;
	for ( int i = 0; i < N; ++ i )
	{
		strideArray[ i ] = stride;
		stride *= lengthArray[ i ];
	}

	this->ComputeZeroOffset();
}

template< typename T, int N >
void FYSmartArray< T, N >::ComputeZeroOffset()
{
	zeroOffSet = 0;

	for ( int i = 0; i < N; ++ i )
	{
		zeroOffSet -= strideArray[ i ] * baseArray[ i ];
	}
}

template< typename T, int N >
void FYSmartArray< T, N >::ComputeNumberOfElements()
{
	numberOfElements = 1;

	for ( int i = 0; i < N; ++ i )
	{
		numberOfElements *= lengthArray[ i ];
	}
}

template< typename T, int N >
void FYSmartArray< T, N >::SliceArray( vector< Range > & rangeArray )
{
	for ( int i = 0; i < N; ++ i )
	{
		Range & range = rangeArray[ i ];

		this->SliceArray( i, range );
	}

	this->ownedData = false;
}

template< typename T, int N >
void FYSmartArray< T, N >::SliceArray( int i, const Range & r )
{
	int first  = r.first();
    int last   = r.last ();
    int stride = r.stride();

	this->lengthArray[ i ] = ( last - first ) / stride + 1;
    this->strideArray[ i ] *= stride;

    int offset = ( first - baseArray[ i ] * stride ) * strideArray[ i ];

    data       += offset;
    zeroOffSet += offset;
}

template< typename T, int N >
void FYSmartArray< T, N >::ZeroAllPointer()
{
	ownedData     = true;
	ownedRootData = true;
	data1 = 0;
	data2 = 0;
	data3 = 0;
	data4 = 0;
}

template< typename T, int N >
void FYSmartArray< T, N >::DeAllocateAll()
{
	if ( ! ownedData ) return;
	if ( ownedRootData )
	{
		//PHSPACE::Deallocate1( data1 );
		//PHSPACE::Deallocate2( data2 );
		//PHSPACE::Deallocate3( data3 );
		//PHSPACE::Deallocate4( data4 );
		delete [] data1;
		delete [] data2;
		delete [] data3;
		delete [] data4;
	}
	else
	{
		//PHSPACE::DeallocateMap1( data1 );
		//PHSPACE::DeallocateMap2( data2 );
		//PHSPACE::DeallocateMap3( data3 );
		//PHSPACE::DeallocateMap4( data4 );
		delete [] data2;
		delete [] data3;
		delete [] data4;
	}
}

template< typename T, int N >
FYSmartArray< T, N >::FYSmartArray( FYSmartArray< T, N > & fatherArray, const Range & r0 )
{
	this->lengthArray = fatherArray.lengthArray;
	this->baseArray   = fatherArray.baseArray;
	this->strideArray = fatherArray.strideArray;
	this->data        = fatherArray.GetData();
	this->zeroOffSet  = fatherArray.zeroOffSet;

	vector< Range > rangeArray;
	rangeArray.push_back( r0 );

	this->SliceArray( rangeArray );
}

template< typename T, int N >
FYSmartArray< T, N >::FYSmartArray( FYSmartArray< T, N > & fatherArray, const Range & r0, const Range & r1 )
{
	this->lengthArray = fatherArray.lengthArray;
	this->baseArray   = fatherArray.baseArray;
	this->strideArray = fatherArray.strideArray;
	this->data        = fatherArray.GetData();
	this->zeroOffSet  = fatherArray.zeroOffSet;

	vector< Range > rangeArray;
	rangeArray.push_back( r0 );
	rangeArray.push_back( r1 );

	this->SliceArray( rangeArray );
}

template< typename T, int N >
FYSmartArray< T, N >::FYSmartArray( FYSmartArray< T, N > & fatherArray, const Range & r0, const Range & r1, const Range & r2 )
{
	this->lengthArray = fatherArray.lengthArray;
	this->baseArray   = fatherArray.baseArray;
	this->strideArray = fatherArray.strideArray;
	this->data        = fatherArray.GetData();
	this->zeroOffSet  = fatherArray.zeroOffSet;

	vector< Range > rangeArray;
	rangeArray.push_back( r0 );
	rangeArray.push_back( r1 );
	rangeArray.push_back( r2 );

	this->SliceArray( rangeArray );
}

template< typename T, int N >
FYSmartArray< T, N >::FYSmartArray( FYSmartArray< T, N > & fatherArray, const Range & r0, const Range & r1, const Range & r2, const Range & r3 )
{
	this->lengthArray = fatherArray.lengthArray;
	this->baseArray   = fatherArray.baseArray;
	this->strideArray = fatherArray.strideArray;
	this->data        = fatherArray.GetData();
	this->zeroOffSet  = fatherArray.zeroOffSet;

	vector< Range > rangeArray;
	rangeArray.push_back( r0 );
	rangeArray.push_back( r1 );
	rangeArray.push_back( r2 );
	rangeArray.push_back( r3 );

	this->SliceArray( rangeArray );
}

template< typename T, int N >
FYSmartArray< T, N > & FYSmartArray< T, N >::operator = ( const T & value )
{
	this->Fill( value );
	return * this;
}

template< typename T, int N >
FYSmartArray< T, N > & FYSmartArray< T, N >::operator = ( const FYSmartArray< T, N > & rhs )
{
	( * this ) = FYArrayExpression< FYFastSmartArrayIterator< T, N > >( rhs.GetFastIterator() );
	return * this;
}

template < typename T, int N > template < typename TExpression >
inline FYSmartArray< T, N > & FYSmartArray< T, N >::operator = ( const FYTypeWrap< TExpression > & tExpression )
{
	this->Evaluate( tExpression.UnWrap(), FYUpdateNewVersion< TValueType, typename TExpression::TValueType >() );

	return * this;
}

template < typename T, int N > template < typename TExpression, typename TUpdate >
inline FYSmartArray< T, N > & FYSmartArray< T, N >::Evaluate( TExpression tExpression, TUpdate )
{
    if ( N == 1 )
	{
		return this->EvaluateWithStackTraversal1( tExpression, TUpdate() );
	}
    else
	{
		return this->EvaluateWithStackTraversalN( tExpression, TUpdate() );
	}
}

template < typename T, int N > template < typename TExpression, typename TUpdate >
inline FYSmartArray< T, N > & FYSmartArray< T, N >::EvaluateWithStackTraversal1( TExpression tExpression, TUpdate )
{
	FYFastSmartArrayIterator< T, N > tIterator( * this );

	int maxRank = 0;

    // Load the strides associated with the innermost loop.
    tIterator  .LoadStride( maxRank );
    tExpression.LoadStride( maxRank );

    const T * last = tIterator.GetData() + lengthArray[ maxRank ] * strideArray[ maxRank ];

  //  while ( tIterator.data() != last )
  //  {
		//TUpdate::Update( * const_cast< T * >( tIterator.GetData() ), * tExpression );
  //      tIterator.advance();
  //      tExpression.advance();
  //  }

	int lastLength = lengthArray[ maxRank ];

    T * data = const_cast< T * >( tIterator.GetData() );

    int commonStride = 1;
    bool useCommonStride = false;

    int upperBound = lastLength * commonStride;

    if ( commonStride == 1 )
    {
        for ( int i = 0; i < upperBound; ++ i )
		{
            TUpdate::Update( * data ++, tExpression.FastRead( i ) );
		}
    }
    else 
	{
        for ( int i = 0; i != upperBound; i += commonStride )
		{
            TUpdate::Update( data[ i ], tExpression.FastRead( i ) );
		}
    }

    return * this;
}


template < typename T, int N > template < typename TExpression, typename TUpdate >
inline FYSmartArray< T, N > & FYSmartArray< T, N >::EvaluateWithStackTraversalN( TExpression tExpression, TUpdate )
{
    const int maxRank = 0;

    // Create an iterator for the array receiving the result
    FYFastSmartArrayIterator< T, N > tIterator( * this );

    for ( int i = 1; i < N; ++ i )
    {
        tIterator  .Push( i );
        tExpression.Push( i );
    }

    // Load the strides associated with the innermost loop.
    tIterator  .LoadStride( maxRank );
    tExpression.LoadStride( maxRank );

	bool useUnitStride = true;

    int commonStride = 1;
    bool useCommonStride = false;

    const T * last[ N ];

    for ( int i = 1; i < N; ++ i )
	{
        last[ i ] = tIterator.GetData() + lengthArray[ i ] * strideArray[ i ];
	}

    int lastLength = lengthArray[ maxRank ];

    int firstNoncollapsedLoop = 1;

    while ( true )
	{
        if ( ( useUnitStride ) || ( useCommonStride ) )
        {
            // Compute the end of the innermost loop
            int upperBound = lastLength * commonStride;

            T * data = const_cast< T * >( tIterator.GetData() );

            if ( commonStride == 1 )
            {
                for ( int i = 0; i < upperBound; ++ i )
				{
                    TUpdate::Update( * data ++, tExpression.FastRead( i ) );
				}
            }
            else 
			{
                for ( int i = 0; i != upperBound; i += commonStride )
				{
                    TUpdate::Update( data[ i ], tExpression.FastRead( i ) );
				}
            }

            tIterator  .Advance( lastLength * commonStride );
            tExpression.Advance( lastLength * commonStride );
		}

        int j = firstNoncollapsedLoop;
        for ( ; j < N; ++ j )
        {
            // Get the next loop
            int r = j;

            // Pop-- this restores the data pointers to the first element
            // encountered in the loop.
            tIterator  .Pop( j );
            tExpression.Pop( j );

            // Load the stride associated with this loop, and increment
            // once.
            tIterator  .LoadStride( r );
            tExpression.LoadStride( r );
            tIterator  .Advance();
            tExpression.Advance();

            // If we aren't at the end of this loop, then stop popping.
            if ( tIterator.GetData() != last[ j ] )
			{
                break;
			}
        }

        if ( j == N )
		{
            break;
		}

        for ( ; j >= firstNoncollapsedLoop; -- j )
        {
            int r2 = j - 1;
            tIterator  .Push( j );
            tExpression.Push( j );
			last[ j - 1 ] = tIterator.GetData() + lengthArray[ r2 ] * strideArray[ r2 ];
        }

        // Load the stride for the innermost loop again.
        tIterator  .LoadStride( maxRank );
        tExpression.LoadStride( maxRank );
    }

    return * this;
}

#define FY_ARRAY_UPDATER_OPERATOR( UpdaterName, Operator )                                                \
template < typename T, int N > template < typename PExpression >                                          \
inline FYSmartArray< T, N > & FYSmartArray< T, N >::operator Operator( const PExpression & pExpression )  \
{                                                                                                         \
    this->Evaluate( typename FYExpressionWrap< PExpression >::TExpression( pExpression ),                 \
      UpdaterName< TValueType, typename FYExpressionWrap< PExpression >::TExpression::TValueType >() );   \
    return * this;                                                                                        \
}

FY_DEFINE_UPDATER_CLASS( FYUpdatePlus    , += );
FY_DEFINE_UPDATER_CLASS( FYUpdateMinus   , -= );
FY_DEFINE_UPDATER_CLASS( FYUpdateMultiply, *= );
FY_DEFINE_UPDATER_CLASS( FYUpdateDivide  , /= );

FY_ARRAY_UPDATER_OPERATOR( FYUpdatePlus    , += )
FY_ARRAY_UPDATER_OPERATOR( FYUpdateMinus   , -= )
FY_ARRAY_UPDATER_OPERATOR( FYUpdateMultiply, *= )
FY_ARRAY_UPDATER_OPERATOR( FYUpdateDivide  , /= )

template < typename T, int N >
void FYCopy( FYSmartArray< T, N > & a, FYSmartArray< T, N > & b, Range & r0, Range & r1, Range & r2 )
{
	int ist = r0.first();
	int ied = r0.last();
	int jst = r1.first();
	int jed = r1.last();
	int kst = r2.first();
	int ked = r2.last();
	for ( int k = kst; k <= ked; ++ k )
	{	
		for ( int j = jst; j <= jed; ++ j )
		{
			for ( int i = ist; i <= ied; ++ i )
			{
				a( i, j, k ) = b( i, j, k );
			}
		}
	}

}

class FYSmartOperation
{
public:
	FYSmartOperation()
	{
		ib = 0;
		jb = 0;
		kb = 0;

		ic = 0;
		jc = 0;
		kc = 0;

		coefb = 1.0;
		coefc = 1.0;
	}
public:
	Range I, J, K;
	int ib, jb, kb;
	int ic, jc, kc;
	double coefb, coefc;
};

template < typename T, int N >
void FYCopy( FYSmartArray< T, N > & a, FYSmartArray< T, N > & b, FYSmartArray< T, N > & c, FYSmartOperation & smartOperation )
{
	int ist = smartOperation.I.first();
	int ied = smartOperation.I.last();
	int jst = smartOperation.J.first();
	int jed = smartOperation.J.last();
	int kst = smartOperation.K.first();
	int ked = smartOperation.K.last();

	int ib = smartOperation.ib;
	int jb = smartOperation.jb;
	int kb = smartOperation.kb;

	int ic = smartOperation.ic;
	int jc = smartOperation.jc;
	int kc = smartOperation.kc;

	double coefb = smartOperation.coefb;
	double coefc = smartOperation.coefc;

	for ( int k = kst; k <= ked; ++ k )
	{	
		for ( int j = jst; j <= jed; ++ j )
		{
			for ( int i = ist; i <= ied; ++ i )
			{
				a( i, j, k ) = coefb * b( i + ib, j + jb, k + kb ) + coefc * c( i + ic, j + jc, k + kc );
			}
		}
	}
}

template < typename T, int N >
void FYCopyIJK( FYSmartArray< T, N > & a, FYSmartArray< T, N > & b, FYSmartArray< T, N > & c, FYSmartOperation & smartOperation )
{
	int ist = smartOperation.I.first();
	int ied = smartOperation.I.last();
	int jst = smartOperation.J.first();
	int jed = smartOperation.J.last();
	int kst = smartOperation.K.first();
	int ked = smartOperation.K.last();

	int ib = smartOperation.ib;
	int jb = smartOperation.jb;
	int kb = smartOperation.kb;

	int ic = smartOperation.ic;
	int jc = smartOperation.jc;
	int kc = smartOperation.kc;

	double coefb = smartOperation.coefb;
	double coefc = smartOperation.coefc;
	for ( int i = ist; i <= ied; ++ i )
	{
		for ( int j = jst; j <= jed; ++ j )
		{
			for ( int k = kst; k <= ked; ++ k )
			{	
				a( i, j, k ) = coefb * b( i + ib, j + jb, k + kb ) + coefc * c( i + ic, j + jc, k + kc );
			}
		}
	}
}

template < typename T, int N >
void FYCopyAdd( FYSmartArray< T, N > & a, FYSmartArray< T, N > & b, FYSmartArray< T, N > & c, FYSmartOperation & smartOperation )
{
	int ist = smartOperation.I.first();
	int ied = smartOperation.I.last();
	int jst = smartOperation.J.first();
	int jed = smartOperation.J.last();
	int kst = smartOperation.K.first();
	int ked = smartOperation.K.last();

	int ib = smartOperation.ib;
	int jb = smartOperation.jb;
	int kb = smartOperation.kb;

	int ic = smartOperation.ic;
	int jc = smartOperation.jc;
	int kc = smartOperation.kc;

	double coefb = smartOperation.coefb;
	double coefc = smartOperation.coefc;

	for ( int k = kst; k <= ked; ++ k )
	{	
		for ( int j = jst; j <= jed; ++ j )
		{
			for ( int i = ist; i <= ied; ++ i )
			{
				a( i, j, k ) += coefb * b( i + ib, j + jb, k + kb ) + coefc * c( i + ic, j + jc, k + kc );
			}
		}
	}
}


template < typename T, int N >
void FYCopy( FYSmartArray< T, N > & a, FYSmartArray< T, N > & b, FYSmartOperation & smartOperation )
{
	int ist = smartOperation.I.first();
	int ied = smartOperation.I.last();
	int jst = smartOperation.J.first();
	int jed = smartOperation.J.last();
	int kst = smartOperation.K.first();
	int ked = smartOperation.K.last();

	int ib = smartOperation.ib;
	int jb = smartOperation.jb;
	int kb = smartOperation.kb;

	double coefb = smartOperation.coefb;

	for ( int k = kst; k <= ked; ++ k )
	{	
		for ( int j = jst; j <= jed; ++ j )
		{
			for ( int i = ist; i <= ied; ++ i )
			{
				a( i, j, k ) = coefb * b( i + ib, j + jb, k + kb );
			}
		}
	}

}

template < typename T >
class FYTemplateFunctionClass
{
public:
	typedef T ( * FYFunction )( const T & x, const T & y );
};

template < typename T, int N >
void FYCopy( FYSmartArray< T, N > & a, FYSmartArray< T, N > & b, FYSmartArray< T, N > & c, FYSmartOperation & smartOperation, typename FYTemplateFunctionClass< T >::FYFunction function )
{
	int ist = smartOperation.I.first();
	int ied = smartOperation.I.last();
	int jst = smartOperation.J.first();
	int jed = smartOperation.J.last();
	int kst = smartOperation.K.first();
	int ked = smartOperation.K.last();

	int ib = smartOperation.ib;
	int jb = smartOperation.jb;
	int kb = smartOperation.kb;

	int ic = smartOperation.ic;
	int jc = smartOperation.jc;
	int kc = smartOperation.kc;

	double coefb = smartOperation.coefb;
	double coefc = smartOperation.coefc;

	for ( int k = kst; k <= ked; ++ k )
	{	
		for ( int j = jst; j <= jed; ++ j )
		{
			for ( int i = ist; i <= ied; ++ i )
			{
				a( i, j, k ) = function( coefb * b( i + ib, j + jb, k + kb ) , coefc * c( i + ic, j + jc, k + kc ) );
			}
		}
	}
}

}
