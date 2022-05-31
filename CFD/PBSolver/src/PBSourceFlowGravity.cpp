#include <iostream>
#include "PBSourceFlowGravity.h"
#include "PBGrad.h"

using namespace std;

namespace PHSPACE
{
	PBSourceFlowGravity::PBSourceFlowGravity(const string& name)
		:PBSource(name)
	{
		this->gx_ = GlobalDataBase::GetDoubleParaFromDB("gravityX");
		this->gy_ = GlobalDataBase::GetDoubleParaFromDB("gravityY");
		this->gz_ = GlobalDataBase::GetDoubleParaFromDB("gravityZ");
	}

	PBSourceFlowGravity::~PBSourceFlowGravity()
	{

	}

	void PBSourceFlowGravity::InitMemory(UnstructGrid* grid, string phiName)
	{
		this->InitProxy(grid);
	}

	void PBSourceFlowGravity::InitProxy(UnstructGrid* grid)
	{
		int totalCells = grid->GetNTotalCell();
		RDouble* cellVolume = grid->GetCellVolume();
		RDouble* rho = reinterpret_cast<RDouble*>(grid->GetDataPtr("rho"));

		RDouble* bfU = reinterpret_cast<RDouble*>(grid->GetDataPtr("bfU"));
		RDouble* bfV = reinterpret_cast<RDouble*>(grid->GetDataPtr("bfV"));
		RDouble* bfW = reinterpret_cast<RDouble*>(grid->GetDataPtr("bfW"));

		RDouble gx = this->gx_;
		RDouble gy = this->gy_;
		RDouble gz = this->gz_;

		//! For natural convection.
		int isBoussinesq = GlobalDataBase::GetIntParaFromDB("isBoussinesq");
		if (isBoussinesq)
		{
			RDouble cte = GlobalDataBase::GetDoubleParaFromDB("thermalExpansion");  //! coefficient of thermal expansion.
			RDouble refT = GlobalDataBase::GetDoubleParaFromDB("refT");
			RDouble* T = reinterpret_cast<RDouble*>(grid->GetDataPtr("T"));
			for (int iCell = 0; iCell < totalCells; ++iCell)
			{
				const RDouble rhoTmp = rho[iCell];
				RDouble bouCoe = -cte * (T[iCell] - refT);
				bfU[iCell] += gx * rhoTmp * bouCoe;
				bfV[iCell] += gy * rhoTmp * bouCoe;
				bfW[iCell] += gz * rhoTmp * bouCoe;
			}
		}
		else
		{
			for (int iCell = 0; iCell < totalCells; ++iCell)
			{
				const RDouble rhoTmp = rho[iCell];
				bfU[iCell] += gx * rhoTmp;
				bfV[iCell] += gy * rhoTmp;
				bfW[iCell] += gz * rhoTmp;
			}
		}
	}

	void PBSourceFlowGravity::Discrt(UnstructGrid* grid, string phiName)
	{
		DiscrtProxy(grid);
	}

	void PBSourceFlowGravity::DiscrtProxy(UnstructGrid* grid)
	{
		int totalCells = grid->GetNTotalCell();
		RDouble* cellVolume = grid->GetCellVolume();
		RDouble* rho = reinterpret_cast<RDouble*>(grid->GetDataPtr("rho"));

		RDouble* bfU = reinterpret_cast<RDouble*>(grid->GetDataPtr("bfU"));
		RDouble* bfV = reinterpret_cast<RDouble*>(grid->GetDataPtr("bfV"));
		RDouble* bfW = reinterpret_cast<RDouble*>(grid->GetDataPtr("bfW"));

		RDouble* scU = reinterpret_cast<RDouble*>(grid->GetDataPtr("scU"));
		RDouble* scV = reinterpret_cast<RDouble*>(grid->GetDataPtr("scV"));
		RDouble* scW = reinterpret_cast<RDouble*>(grid->GetDataPtr("scW"));

		RDouble gx = this->gx_;
		RDouble gy = this->gy_;
		RDouble gz = this->gz_;
		for (int iCell = 0; iCell < totalCells; ++iCell)
		{
			scU[iCell] += bfU[iCell] * cellVolume[iCell];
			scV[iCell] += bfV[iCell] * cellVolume[iCell];
			scW[iCell] += bfW[iCell] * cellVolume[iCell];
		}
	}

	DEFINE_PBSOLVER_REGISTER_FUNCTION(PBSource, PBSourceFlowGravity, SOURCE_FLOW_GRAVITY);
}
