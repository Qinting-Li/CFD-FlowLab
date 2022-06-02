#include <iostream>
#include "PBSourceTurbKDefault.h"
#include "PBGrad.h"

using namespace std;

namespace PHSPACE
{
	PBSourceTurbKDefault::PBSourceTurbKDefault(const string& name)
		:PBSource(name)
	{
	}

	PBSourceTurbKDefault::~PBSourceTurbKDefault()
	{

	}

	void PBSourceTurbKDefault::InitMemory(UnstructGrid* grid,string phiName)
	{
		
	}

	void PBSourceTurbKDefault::Discrt(UnstructGrid* grid,string phiName)
	{
		RDouble* k = reinterpret_cast<RDouble*>(grid->GetDataPtr(phiName));
		RDouble* gen = reinterpret_cast<RDouble*>(grid->GetDataPtr("gen"));
		RDouble* ap = reinterpret_cast<RDouble*>(grid->GetDataPtr("ap"));
		RDouble* sc = reinterpret_cast<RDouble*>(grid->GetDataPtr("sc"));
		RDouble* rho = reinterpret_cast<RDouble*>(grid->GetDataPtr("rho"));
		RDouble* epsilon = reinterpret_cast<RDouble*>(grid->GetDataPtr("epsilon"));
		RDouble* vol = reinterpret_cast<RDouble*>(grid->GetCellVolume());

		for (int iCell = 0; iCell < grid->GetNTotalCell(); ++iCell)
		{
			sc[iCell] += gen[iCell] * vol[iCell];
			ap[iCell] += rho[iCell] * epsilon[iCell] / (k[iCell] + 1e-20) * vol[iCell];
		}
	}

	DEFINE_PBSOLVER_REGISTER_FUNCTION(PBSource, PBSourceTurbKDefault, SOURCE_TURB_K_DEFAULT);
}