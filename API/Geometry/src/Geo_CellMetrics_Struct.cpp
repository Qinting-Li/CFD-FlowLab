#include "GlobalDataBase.h"
#include "Geo_CellMetrics_Struct.h"
namespace PHSPACE
{
Geo_CellMetrics_Struct::Geo_CellMetrics_Struct()
{
	xcc = 0;
	ycc = 0;
	zcc = 0;
	vol = 0;

	xlen = 0;
	ylen = 0;
	zlen = 0;

	//! add by Min Yaobing
	jacobian = 0;
	//! add by Min Yaobing

    largestLocalGridLength		= 0;
    subgridLength               = 0;
}
Geo_CellMetrics_Struct::~Geo_CellMetrics_Struct()
{
	DeallocateAll();

    delete largestLocalGridLength; largestLocalGridLength = 0;
    delete subgridLength         ; subgridLength          = 0;
}

void Geo_CellMetrics_Struct::AllocateMetrics(Range &I, Range &J, Range &K)
{
	DeallocateAll();
	xcc = new RDouble3D(I,J,K,fortranArray);
	ycc = new RDouble3D(I,J,K,fortranArray);
	zcc = new RDouble3D(I,J,K,fortranArray);
    vol = new RDouble3D(I,J,K,fortranArray);

	xlen = new RDouble3D(I,J,K,fortranArray);
	ylen = new RDouble3D(I,J,K,fortranArray);
	zlen = new RDouble3D(I,J,K,fortranArray);

	//! add by Min Yaobing
	int SolverStructOrder = PHSPACE::GlobalDataBase::GetIntParaFromDB( "SolverStructOrder" );
	if ( SolverStructOrder > 2 )
	{
		jacobian = new RDouble3D(I,J,K,fortranArray);
	}
	//! add by Min Yaobing

}
}