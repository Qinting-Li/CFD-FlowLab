#include <iostream>
#include "PBSourceFlowDefault.h"

using namespace std;

namespace PHSPACE
{
	PBSourceFlowDefault::PBSourceFlowDefault(const string& name)
		:PBSource(name)
	{

	}

	PBSourceFlowDefault::~PBSourceFlowDefault()
	{

	}

	void PBSourceFlowDefault::Discrt(UnstructGrid* grid, string phiName)
	{
		DiscrtProxy(grid);
	}

	void PBSourceFlowDefault::DiscrtProxy(UnstructGrid* grid)
	{
		int totalCells = grid->GetNTotalCell();
		RDouble* cellVolume = grid->GetCellVolume();
		RDouble* dpdx = reinterpret_cast<RDouble*>(grid->GetDataPtr("dpdx"));
		RDouble* dpdy = reinterpret_cast<RDouble*>(grid->GetDataPtr("dpdy"));
		RDouble* dpdz = reinterpret_cast<RDouble*>(grid->GetDataPtr("dpdz"));

		RDouble* scU = reinterpret_cast<RDouble*>(grid->GetDataPtr("scU"));
		RDouble* scV = reinterpret_cast<RDouble*>(grid->GetDataPtr("scV"));
		RDouble* scW = reinterpret_cast<RDouble*>(grid->GetDataPtr("scW"));

		for (int iCell = 0; iCell < totalCells; ++iCell)
		{
			RDouble cv = cellVolume[iCell];
			scU[iCell] -= dpdx[iCell] * cv;
			scV[iCell] -= dpdy[iCell] * cv;
			scW[iCell] -= dpdz[iCell] * cv;
		}
	}

	DEFINE_PBSOLVER_REGISTER_FUNCTION(PBSource, PBSourceFlowDefault, SOURCE_FLOW_DEFAULT);
}