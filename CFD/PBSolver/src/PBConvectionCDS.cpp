#include <iostream>
#include <algorithm>

#include "PBConvectionCDS.h"
#include "Geo_UnstructGrid.h"

using namespace std;

namespace PHSPACE
{
	PBConvectionCDS::PBConvectionCDS(const string& name)
		:PBConvection(name)
	{

	}

	PBConvectionCDS::~PBConvectionCDS()
	{

	}

	void PBConvectionCDS::Discrt(UnstructGrid* grid, string phiName, string apName, string scName, RDouble* v,
		                         RDouble* apV, RDouble* scV, RDouble* w, RDouble* apW, RDouble* scW)
	{
		int faceBegin = grid->GetNBoundFace();
		int faceEnd = grid->GetNTotalFace();
		RDouble* flux    = reinterpret_cast<RDouble*>(grid->GetDataPtr("flux"));
		RDouble* ajUpper = reinterpret_cast<RDouble*>(grid->GetDataPtr("ajLeft"));
		RDouble* ajLower = reinterpret_cast<RDouble*>(grid->GetDataPtr("ajRight"));
		RDouble* phi     = reinterpret_cast<RDouble*>(grid->GetDataPtr(phiName));
		RDouble* ap      = reinterpret_cast<RDouble*>(grid->GetDataPtr("ap"));

		int total = faceEnd - faceBegin;
		for (int iFace = 0; iFace < total; ++iFace)
		{
			int iFacelocal = iFace - faceBegin;
			RDouble a12 = std::max(flux[iFace], 0.0);
			RDouble a21 = a12 - flux[iFace];
			ajUpper[iFacelocal] += a21;
			ajLower[iFacelocal] += a12;
		}

		RDouble* deltaL = reinterpret_cast<RDouble*>(grid->GetDataPtr("deltL"));
		RDouble* deltaR = reinterpret_cast<RDouble*>(grid->GetDataPtr("deltR"));

		if (scV == nullptr) 
		{
			RDouble* sc = reinterpret_cast<RDouble*>(grid->GetDataPtr("sc"));
			for (int iFace = 0; iFace < total; ++iFace)
			{
				int iFaceReal = iFace + grid->GetNBoundFace();
				int le = grid->GetLeftCellOfFace()[iFace];
				int re = grid->GetRightCellOfFace()[iFace];
				RDouble wfL = deltaL[iFace];
				RDouble wfR = deltaR[iFace];
				RDouble a12 = (-std::min(flux[iFace], 0.0) + flux[iFace] * wfR)*(phi[le] - phi[re]);

				sc[le] += a12;
				sc[re] -= a12;
			}
		}
		else
		{
			RDouble* sc = reinterpret_cast<RDouble*>(grid->GetDataPtr("scU"));
			for (int iFace = 0; iFace < total; ++iFace)
			{
				int iFaceReal = iFace + grid->GetNBoundFace();
				int le = grid->GetLeftCellOfFace()[iFace];
				int re = grid->GetRightCellOfFace()[iFace];
				RDouble wfL = deltaL[iFace];
				RDouble wfR = deltaR[iFace];
				RDouble a12 = ((-std::min(flux[iFace], 0.0) + flux[iFace] * wfR))* (phi[le] - phi[re]);
				sc[le] += a12;
				sc[re] -= a12;

				a12 = ((-std::min(flux[iFace], 0.0) + flux[iFace] * wfR)) * (v[le] - v[re]);
				scV[le] += a12;
				scV[re] -= a12;

				a12 = ((-std::min(flux[iFace], 0.0) + flux[iFace] * wfR)) * (w[le] - w[re]);
				scW[le] += a12;
				scW[re] -= a12;
			}
		}
	}

	DEFINE_PBSOLVER_REGISTER_FUNCTION(PBConvection, PBConvectionCDS, CONVECTION_CDS);
}