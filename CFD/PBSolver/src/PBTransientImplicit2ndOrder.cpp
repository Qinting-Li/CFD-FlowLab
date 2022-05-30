#include <iostream>
#include "PBTransientImplicit2ndOrder.h"
#include "Geo_UnstructGrid.h"

using namespace std;

namespace PHSPACE
{
	PBTransientImplicit2ndOrder::PBTransientImplicit2ndOrder(const string& name)
		:PBTransient(name)
	{

	}

	PBTransientImplicit2ndOrder::~PBTransientImplicit2ndOrder()
	{

	}
	
	void PBTransientImplicit2ndOrder::AllocateMemory(UnstructGrid* grid, vector<string>& phiName, vector<RDouble>& initPhi)
	{
		int iVar = 0;
		vector<RDouble*> phiOld;
		vector<RDouble*> phiOldOld;
		RDouble nTotal = grid->GetNTotalCell() + grid->GetNBoundFace();

		for (vector<string>::iterator iter = phiName.begin(); iter != phiName.end(); ++iter)
		{			
			RDouble* phiOld_ = ObjAllocateMemory((*iter) + "Old", nTotal, initPhi[iVar], grid);
			phiOld.push_back(phiOld_);
			RDouble* phiOldOld_ = ObjAllocateMemory((*iter) + "OldOld", nTotal, initPhi[iVar], grid);
			phiOldOld.push_back(phiOldOld_);
		}

		vector<RDouble*>().swap(phiOld);
		vector<RDouble*>().swap(phiOldOld);
	}

	void PBTransientImplicit2ndOrder::Discrt(UnstructGrid* grid, string phiName, string phiOldName, string phiOldOldName, string rhoOldName, string rhoOldOldName)
	{
		RDouble* phi = reinterpret_cast<RDouble*>(grid->GetDataPtr(phiName));
		RDouble* cvOld = grid->GetCellVolume();
		RDouble* cvOldOld = grid->GetCellVolume();
		RDouble* phiOld = reinterpret_cast<RDouble*>(grid->GetDataPtr(phiOldName));
		RDouble* phiOldOld = reinterpret_cast<RDouble*>(grid->GetDataPtr(phiOldOldName));
		
		RDouble* rhoOld = reinterpret_cast<RDouble*>(grid->GetDataPtr(rhoOldName));
		RDouble* rhoOldOld = reinterpret_cast<RDouble*>(grid->GetDataPtr(rhoOldOldName));
		int cellBegin = 0;
		int cellEnd = grid->GetNTotalCell();
		RDouble dt = GlobalDataBase::GetDoubleParaFromDB("dt");
		if (phiName == "U")
		{
			RDouble* ap = reinterpret_cast<RDouble*>(grid->GetDataPtr("apU"));
			RDouble* sc = reinterpret_cast<RDouble*>(grid->GetDataPtr("scU"));
			RDouble* vOld = reinterpret_cast<RDouble*>(grid->GetDataPtr("VOld"));
			RDouble* wOld = reinterpret_cast<RDouble*>(grid->GetDataPtr("WOld"));
			RDouble* apV = reinterpret_cast<RDouble*>(grid->GetDataPtr("apV"));
			RDouble* apW = reinterpret_cast<RDouble*>(grid->GetDataPtr("apW"));
			RDouble* scV = reinterpret_cast<RDouble*>(grid->GetDataPtr("scV"));
			RDouble* scW = reinterpret_cast<RDouble*>(grid->GetDataPtr("scW"));

			for (int iCell = cellBegin; iCell < cellEnd; ++iCell)
			{
				const RDouble tmp = rhoOld[iCell] * cvOld[iCell] / dt;
				ap[iCell] += tmp;
				apV[iCell] += tmp;
				apW[iCell] += tmp;
				sc[iCell] += tmp * phiOld[iCell];
				scV[iCell] += tmp * vOld[iCell];
				scW[iCell] += tmp * wOld[iCell];
			}
		}
		else
		{
			RDouble* ap = reinterpret_cast<RDouble*>(grid->GetDataPtr("ap"));
			RDouble* sc = reinterpret_cast<RDouble*>(grid->GetDataPtr("sc"));
			if (rhoOld != nullptr && rhoOldOld != nullptr)
			{
				for (int iCell = cellBegin; iCell < cellEnd; ++iCell)
				{
					const RDouble tmpOld = rhoOld[iCell] * cvOld[iCell];
					const RDouble tmpOldOld = rhoOldOld[iCell] * cvOldOld[iCell];
					ap[iCell] += (RDouble(4.0) * tmpOld - tmpOldOld) / (2 * dt);
					sc[iCell] += (RDouble(4.0) * tmpOld * phiOld[iCell] - tmpOldOld * phiOldOld[iCell]) / (2 * dt);
				}
			}
		}
	}

	void PBTransientImplicit2ndOrder::SaveOldTimeValue(UnstructGrid* grid, vector<string>& phiName)
	{
		vector<RDouble*> phiNow;
		vector<RDouble*> phiOld;
		vector<RDouble*> phiOldOld;

		for (vector<string>::iterator iter = phiName.begin(); iter != phiName.end(); ++iter)
		{
			RDouble* phiNowPtr = reinterpret_cast<RDouble*>(grid->GetDataPtr((*iter)));
			RDouble* phiOldPtr = reinterpret_cast<RDouble*>(grid->GetDataPtr((*iter) + "Old"));
			RDouble* phiOldOldPtr = reinterpret_cast<RDouble*>(grid->GetDataPtr((*iter) + "OldOld"));
			phiNow.push_back(phiNowPtr);
			phiOld.push_back(phiOldPtr);
			phiOldOld.push_back(phiOldOldPtr);
		}

		for (int index = 0; index < phiName.size(); ++index)
		{
			for (int iCell = 0; iCell < grid->GetNTotalCell(); ++iCell)
			{
				phiOldOld[index][iCell] = phiOld[index][iCell];
				phiOld[index][iCell]    = phiNow[index][iCell];
			}
		}

		vector<RDouble*>().swap(phiNow);
		vector<RDouble*>().swap(phiOld);
		vector<RDouble*>().swap(phiOldOld);
	}

	DEFINE_PBSOLVER_REGISTER_FUNCTION(PBTransient, PBTransientImplicit2ndOrder, TRANSIENT_IMPLICIT_2ND_ORDER);
}