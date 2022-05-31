#include "PBTransientCrankNicolson.h"
#include "Geo_UnstructGrid.h"

namespace PHSPACE
{
    PBTransientCrankNicolson::PBTransientCrankNicolson(const string& name)
        :PBTransient(name)
    {

    }

    PBTransientCrankNicolson::~PBTransientCrankNicolson()
    {

    }

    void PBTransientCrankNicolson::AllocateMemory(UnstructGrid* grid, vector<string>& phiName, vector<RDouble>& initPhi)
    {
		int iVar = 0;
		vector<RDouble*> phi;
		RDouble nTotal = grid->GetNTotalCell() + grid->GetNBoundFace();
        for (vector<string>::iterator iter = phiName.begin(); iter != phiName.end(); ++iter)
		{
			RDouble* phi_ = ObjAllocateMemory((*iter) + "Old", nTotal, initPhi[iVar], grid);
			phi.push_back(phi_);
		}
		vector<RDouble*>().swap(phi);
    }

    void PBTransientCrankNicolson::Discrt(UnstructGrid* grid, string phiName, string phiOldName, string phiOldOldName, string rhoOldName, string rhoOldOldName)
    {
		RDouble* phi = reinterpret_cast<RDouble*>(grid->GetDataPtr(phiName));
		RDouble* phif = reinterpret_cast<RDouble*>(grid->GetDataPtr(phiName + "f"));
		RDouble* alphaFlux = reinterpret_cast<RDouble*>(grid->GetDataPtr("flux"));
		
		RDouble* phiOld = reinterpret_cast<RDouble*>(grid->GetDataPtr(phiOldName));
		RDouble* phiOldOld = reinterpret_cast<RDouble*>(grid->GetDataPtr(phiOldOldName));

		RDouble* cvOld = grid->GetCellVolume();
		RDouble* cvOldOld = grid->GetCellVolume();

		RDouble* rhoOld = reinterpret_cast<RDouble*>(grid->GetDataPtr(rhoOldName));
		RDouble* rhoOldOld = reinterpret_cast<RDouble*>(grid->GetDataPtr(rhoOldOldName));

		RDouble* ap = reinterpret_cast<RDouble*>(grid->GetDataPtr("ap"));
		RDouble* sc = reinterpret_cast<RDouble*>(grid->GetDataPtr("sc"));

		int cellBegin = 0;
		int cellEnd = grid->GetNTotalCell();

		RDouble dt = GlobalDataBase::GetDoubleParaFromDB("dt");

        if (rhoOld == NULL)
		{
			for (int iCell = cellBegin; iCell < cellEnd; ++iCell)
			{
				RDouble tmp = 2.0 * cvOld[iCell] / dt;
				ap[iCell] += tmp;
				sc[iCell] += tmp * phiOld[iCell];
			}
		}
		else
		{
			for (int iCell = cellBegin; iCell < cellEnd; ++iCell)
			{
				RDouble tmp = 2.0 * rhoOld[iCell] * cvOld[iCell] / dt;
				ap[iCell] += tmp;
				sc[iCell] += tmp * phiOld[iCell];
			}
		}

		RDouble* scOld = new RDouble[grid->GetNTotalCell()];
		for (int iCell = 0; iCell < grid->GetNTotalCell(); ++iCell)
		{
			scOld[iCell] = 0.0;
		}

		for (int iFace = 0; iFace < grid->GetNBoundFace(); ++iFace)
		{
			int le = grid->GetLeftCellOfFace()[iFace];
			int re = grid->GetRightCellOfFace()[iFace];
			scOld[le] += alphaFlux[iFace] * phif[iFace];
		}

		for (int iFace = grid->GetNBoundFace(); iFace < grid->GetNTotalFace(); ++iFace)
		{
			int le = grid->GetLeftCellOfFace()[iFace];
			int re = grid->GetRightCellOfFace()[iFace];
			RDouble tmp = alphaFlux[iFace] * phif[iFace];
			scOld[le] += tmp;
			scOld[re] -= tmp;
		}

		for (int iCell = 0; iCell < grid->GetNTotalCell(); ++iCell)
		{
			sc[iCell] -= scOld[iCell];
		}

		delete[] scOld;
    }

    void PBTransientCrankNicolson::SaveOldTimeValue(UnstructGrid* grid, vector<string>& phiName)
    {
        std::vector<RDouble*> phiNow;
        std::vector<RDouble*> phiOld;

		for (std::vector<std::string>::iterator iter = phiName.begin(); iter != phiName.end(); ++iter)
		{
			RDouble* phi = (RDouble*)grid->GetDataPtr((*iter));
			phiNow.push_back(phi);
			RDouble* phiold = (RDouble*)grid->GetDataPtr((*iter) + "Old");
			phiOld.push_back(phiold);
		}

		for (int index = 0; index < phiName.size(); ++index)
		{
			for (int iCell = 0; iCell < grid->GetNTotalCell(); ++iCell)
			{
				phiOld[index][iCell] = phiNow[index][iCell];
			}
		}

		vector<RDouble*>().swap(phiNow);
		vector<RDouble*>().swap(phiOld);
    }

	DEFINE_PBSOLVER_REGISTER_FUNCTION(PBTransient, PBTransientCrankNicolson, TRANSIENT_CRANK_NICOLSON);
}

