#include <iostream>
#include <algorithm>

#include "PBConvectionQUICK.h"
#include "Geo_UnstructGrid.h"

using namespace std;

namespace PHSPACE
{
	PBConvectionQUICK::PBConvectionQUICK(const string& name)
		:PBConvection(name)
	{

	}

	PBConvectionQUICK::~PBConvectionQUICK()
	{

	}

	void PBConvectionQUICK::Discrt(UnstructGrid* grid, string phiName, string apName, string scName, RDouble* v,
		                           RDouble* apV, RDouble* scV, RDouble* w, RDouble* apW, RDouble* scW)
	{
		int faceBegin = grid->GetNBoundFace();
		int faceEnd = grid->GetNTotalFace();
		RDouble* flux = reinterpret_cast<RDouble*>(grid->GetDataPtr("flux"));
		RDouble* ajUpper = reinterpret_cast<RDouble*>(grid->GetDataPtr("ajLeft"));
		RDouble* ajLower = reinterpret_cast<RDouble*>(grid->GetDataPtr("ajRight"));
		RDouble* phi = reinterpret_cast<RDouble*>(grid->GetDataPtr(phiName));
		RDouble* ap = reinterpret_cast<RDouble*>(grid->GetDataPtr("ap"));
		int total = faceEnd - faceBegin;
		for (int iFace = faceBegin; iFace < faceEnd; ++iFace)
		{
			int iFacelocal = iFace - faceBegin;
			RDouble a12 = std::max(flux[iFace], 0.0);
			RDouble a21 = a12 - flux[iFace];
			ajUpper[iFacelocal] += a21;
			ajLower[iFacelocal] += a12;
		}

		RDouble* ex = reinterpret_cast<RDouble*>(grid->GetCellCenterX());
		RDouble* ey = reinterpret_cast<RDouble*>(grid->GetCellCenterY());
		RDouble* ez = reinterpret_cast<RDouble*>(grid->GetCellCenterZ());

		if (scV == NULL)
		{
			RDouble* dphidx = reinterpret_cast<RDouble*>(grid->GetDataPtr("dudx"));
			RDouble* dphidy = reinterpret_cast<RDouble*>(grid->GetDataPtr("dudy"));
			RDouble* dphidz = reinterpret_cast<RDouble*>(grid->GetDataPtr("dudz"));
			RDouble* sc = reinterpret_cast<RDouble*>(grid->GetDataPtr("sc"));
			for (int iFace = faceBegin; iFace < faceEnd; ++iFace)
			{
				int le = grid->GetLeftCellOfFace()[iFace];
				int re = grid->GetRightCellOfFace()[iFace];
				RDouble dx = ex[re] - ex[le];
				RDouble dy = ey[re] - ey[le];
				RDouble dz = ez[re] - ez[le];
				RDouble cphi = phi[re] - 2 * (dphidx[le] * dx + dphidy[le] * dy + dphidz[le] * dz);
				RDouble cphi_r = (phi[le] - cphi) / (phi[re] - phi[le] + SMALL);

				RDouble a12 = (((3+ cphi_r)/8)*std::max(flux[iFace], 0.0) * (phi[le] - phi[re]) - ((3+ cphi_r)/8)*std::max(-flux[iFace], 0.0) * (phi[re] - phi[le]));

				sc[le] += a12;
				sc[re] -= a12;
			}
		}
		else
		{
			RDouble* dphidx = reinterpret_cast<RDouble*>(grid->GetDataPtr("dudx"));
			RDouble* dphidy = reinterpret_cast<RDouble*>(grid->GetDataPtr("dudy"));
			RDouble* dphidz = reinterpret_cast<RDouble*>(grid->GetDataPtr("dudz"));

			RDouble* dvdx = reinterpret_cast<RDouble*>(grid->GetDataPtr("dvdx"));
			RDouble* dvdy = reinterpret_cast<RDouble*>(grid->GetDataPtr("dvdy"));
			RDouble* dvdz = reinterpret_cast<RDouble*>(grid->GetDataPtr("dvdz"));

			RDouble* dwdx = reinterpret_cast<RDouble*>(grid->GetDataPtr("dwdx"));
			RDouble* dwdy = reinterpret_cast<RDouble*>(grid->GetDataPtr("dwdy"));
			RDouble* dwdz = reinterpret_cast<RDouble*>(grid->GetDataPtr("dwdz"));

			RDouble* sc = reinterpret_cast<RDouble*>(grid->GetDataPtr("scU"));
			for (int iFace = faceBegin; iFace < faceEnd; ++iFace)
			{
				int le = grid->GetLeftCellOfFace()[iFace];
				int re = grid->GetRightCellOfFace()[iFace];
				RDouble dx = ex[re] - ex[le];
				RDouble dy = ey[re] - ey[le];
				RDouble dz = ez[re] - ez[le];
				RDouble cphi = phi[re] - 2 * (dphidx[le] * dx + dphidy[le] * dy + dphidz[le] * dz);
				RDouble cphi_r = (phi[le] - cphi) / (phi[re] - phi[le] + SMALL);
				RDouble cv = v[re] - 2 * (dvdx[le] * dx + dvdy[le] * dy + dvdz[le] * dz);
				RDouble cv_r = (v[le] - cv) / (v[re] - v[le] + SMALL);
				RDouble cw = w[re] - 2 * (dwdx[le] * dx + dwdy[le] * dy + dwdz[le] * dz);
				RDouble cw_r = (w[le] - cw) / (w[re] - w[le] + SMALL);

				RDouble a12 = ((3 + cphi_r) / 8)*std::max(flux[iFace], 0.0) * (phi[le] - phi[re]) - ((3+ cphi_r) / 8)*std::max(-flux[iFace], 0.0)  * (phi[re] - phi[le]);
				sc[le] += a12;
				sc[re] -= a12;

				a12 = ((3 + cv_r) / 8) * std::max(flux[iFace], 0.0) * (v[le] - v[re]) - ((3 + cv_r) / 8) * std::max(-flux[iFace], 0.0)* (v[re] - v[le]);
				scV[le] += a12;
				scV[re] -= a12;

				a12 = ((3 + cw_r) / 8) * std::max(flux[iFace], 0.0) * (w[le] - w[re]) - ((3 + cw_r) / 8) * std::max(-flux[iFace], 0.0)* (w[re] - w[le]);
				scW[le] += a12;
				scW[re] -= a12;
			}
		}

	}

	DEFINE_PBSOLVER_REGISTER_FUNCTION(PBConvection, PBConvectionQUICK, CONVECTION_QUICK);
}