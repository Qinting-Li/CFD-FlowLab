#include <iostream>

#include "PBGradLsq.h"
#include "Geo_UnstructGrid.h"

using namespace std;

namespace PHSPACE
{
	PBGradLsq::PBGradLsq(const string& name)
		:PBGrad(name)
	{

	}

	PBGradLsq::~PBGradLsq()
	{

	}

	void PBGradLsq::Calc(UnstructGrid* grid, string phiName, string dphidxName, string dphidyName, string dphidzName, bool isComm)
	{
		RDouble* phi = reinterpret_cast<RDouble*>(grid->GetDataPtr(phiName));
		RDouble* dphidx = reinterpret_cast<RDouble*>(grid->GetDataPtr(dphidxName));
		RDouble* dphidy = reinterpret_cast<RDouble*>(grid->GetDataPtr(dphidyName));
		RDouble* dphidz = reinterpret_cast<RDouble*>(grid->GetDataPtr(dphidzName));

		string gradientName = "lsq";
		grid->UpdateData("gradientName", &gradientName, PHSTRING, 1);
		grid->ComputeWeight();
		grid->CompGradientLSQ(phi, dphidx, dphidy, dphidz);
	}

	DEFINE_PBSOLVER_REGISTER_FUNCTION(PBGrad, PBGradLsq, GRAD_LSQ);
}