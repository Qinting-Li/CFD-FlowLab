#include "Geo_OversetGridTopo_Struct.h"
#include "OversetGridFactory.h"
namespace PHSPACE
{

Geo_OversetGridTopo_Struct::Geo_OversetGridTopo_Struct()
{
	numberOfCores = 0;
	zoneStartPointLabel = 0;
	iLinkPointLabel = 0;
	jLinkPointLabel = 0;
	kLinkPointLabel = 0;
	cellTypeContainer = 0;
	zoneStartCenterLabel = 0;
	hingedPointContainer = 0;
	xCoreContainer = 0;
	yCoreContainer = 0;
	zCoreContainer = 0;
}

Geo_OversetGridTopo_Struct::~Geo_OversetGridTopo_Struct()
{
	if( iLinkPointLabel ) delete [] iLinkPointLabel;
	if( jLinkPointLabel ) delete [] jLinkPointLabel;
	if( kLinkPointLabel ) delete [] kLinkPointLabel;
	if(cellTypeContainer) delete [] cellTypeContainer;
	if( hingedPointContainer ) delete [] hingedPointContainer;
	if( xCoreContainer ) delete [] xCoreContainer;
	if( yCoreContainer ) delete [] yCoreContainer;
	if( zCoreContainer ) delete [] zCoreContainer;
}

void Geo_OversetGridTopo_Struct::CreateLinkPointStructure(int nTotalNode)
{
	if( iLinkPointLabel ) delete [] iLinkPointLabel;
	if( jLinkPointLabel ) delete [] jLinkPointLabel;
	if( kLinkPointLabel ) delete [] kLinkPointLabel;
	iLinkPointLabel = new int [nTotalNode];
	jLinkPointLabel = new int [nTotalNode];
	kLinkPointLabel = new int [nTotalNode];
}

void Geo_OversetGridTopo_Struct::CreateCellContainer(int nTotalCell)
{
	if(cellTypeContainer) delete [] cellTypeContainer;
	cellTypeContainer = new int [nTotalCell];

	for( int iCell = 0; iCell < nTotalCell; ++iCell )
	{
		cellTypeContainer[iCell] = GENERIC_COLOR;
	}//very important part Guo Yongheng 20151202

}

void Geo_OversetGridTopo_Struct::CreateCoreParameterContainer(int numberOfCores)
{
	if( hingedPointContainer ) delete [] hingedPointContainer;
	if( xCoreContainer ) delete [] xCoreContainer;
	if( yCoreContainer ) delete [] yCoreContainer;
	if( zCoreContainer ) delete [] zCoreContainer;

	hingedPointContainer = new int [numberOfCores];
	xCoreContainer = new RDouble [numberOfCores];
	yCoreContainer = new RDouble [numberOfCores];
	zCoreContainer = new RDouble [numberOfCores];
}
}