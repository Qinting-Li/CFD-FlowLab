#include <iostream>
#include "PBTransientImplicitEuler.h"
#include "Geo_UnstructGrid.h"

using namespace std;

namespace PHSPACE
{
    PBTransientImplicitEuler::PBTransientImplicitEuler(const string &name):PBTransient(name)
    {

    }

    PBTransientImplicitEuler::~PBTransientImplicitEuler()
    {

    }

    void PBTransientImplicitEuler::AllocateMemory(UnstructGrid *grid, vector<string> &phiName, vector<RDouble> &initPhi)
    {
        int iVar = 0;
        vector<RDouble *> phiOld;
        RDouble nTotal = grid->GetNTotalCell() + grid->GetNBoundFace();
        for (vector<string>::iterator iter = phiName.begin(); iter != phiName.end(); ++ iter)
        {
            RDouble *phi_ = ObjAllocateMemory((*iter) + "Old", nTotal, initPhi[iVar], grid);
            phiOld.push_back(phi_);
        }
        vector<RDouble *>().swap(phiOld);
    }

    void PBTransientImplicitEuler::Discrt(UnstructGrid *grid, string phiName, string phiOldName, string phiOldOldName, string rhoOldName, string rhoOldOldName)
    {
        RDouble *phi = reinterpret_cast<RDouble *>(grid->GetDataPtr(phiName));
        RDouble *cvOld = grid->GetCellVolume();
        RDouble *cvOldOld = grid->GetCellVolume();

        RDouble *phiOld = reinterpret_cast<RDouble *>(grid->GetDataPtr(phiOldName));
        RDouble *phiOldOld = reinterpret_cast<RDouble *>(grid->GetDataPtr(phiOldOldName));

        RDouble *rhoOld = reinterpret_cast<RDouble *>(grid->GetDataPtr(rhoOldName));
        RDouble *rhoOldOld = reinterpret_cast<RDouble *>(grid->GetDataPtr(rhoOldOldName));

        int cellBegin = 0;
        int cellEnd = grid->GetNTotalCell();

        RDouble dt = GlobalDataBase::GetDoubleParaFromDB("dt");
        if (phiName == "U")
        {
            RDouble* vOld = reinterpret_cast<RDouble*>(grid->GetDataPtr("VOld"));
            RDouble* wOld = reinterpret_cast<RDouble*>(grid->GetDataPtr("WOld"));
            RDouble* ap = reinterpret_cast<RDouble*>(grid->GetDataPtr("apU"));
            RDouble* sc = reinterpret_cast<RDouble*>(grid->GetDataPtr("scU"));
            RDouble* apV = reinterpret_cast<RDouble*>(grid->GetDataPtr("apV"));
            RDouble* apW = reinterpret_cast<RDouble*>(grid->GetDataPtr("apW"));
            RDouble* scV = reinterpret_cast<RDouble*>(grid->GetDataPtr("scV"));
            RDouble* scW = reinterpret_cast<RDouble*>(grid->GetDataPtr("scW"));

            for (int iCell = cellBegin; iCell < cellEnd; ++ iCell)
            {
                const RDouble tmp = rhoOld[iCell] * cvOld[iCell] / dt;
                ap[iCell]  += tmp;
                apV[iCell] += tmp;
                apW[iCell] += tmp;
                sc[iCell]  += tmp * phiOld[iCell];
                scV[iCell] += tmp * vOld[iCell];
                scW[iCell] += tmp * wOld[iCell];
            }
        }
        else
        {
            RDouble *ap = reinterpret_cast<RDouble *>(grid->GetDataPtr("ap"));
            RDouble *sc = reinterpret_cast<RDouble *>(grid->GetDataPtr("sc"));

            if (rhoOld == NULL)
            {
                for (int iCell = cellBegin; iCell < cellEnd; ++ iCell)
                {
                    const RDouble tmp = cvOld[iCell] / dt;
                    ap[iCell] += tmp;
                    sc[iCell] += tmp * phiOld[iCell];
                }
            }
            else
            {
                for (int iCell = cellBegin; iCell < cellEnd; ++ iCell)
                {
                    const RDouble tmp = rhoOld[iCell] * cvOld[iCell] / dt;
                    ap[iCell] += tmp;
                    sc[iCell] += tmp * phiOld[iCell];
                }
            }
        }
    }

    void PBTransientImplicitEuler::SaveOldTimeValue(UnstructGrid *grid, vector<string> &phiName)
    {
        vector<RDouble *> phiNow;
        vector<RDouble *> phiOld;
        for (vector<string>::iterator iter = phiName.begin(); iter != phiName.end(); ++ iter)
        {
            RDouble *phi = reinterpret_cast<RDouble *>(grid->GetDataPtr((*iter)));
            RDouble *phiOldPtr = reinterpret_cast<RDouble *>(grid->GetDataPtr((*iter) + "Old"));
            phiNow.push_back(phi);
            phiOld.push_back(phiOldPtr);
        }

        for (int index = 0; index < phiName.size(); ++ index)
        {
            for (int iCell = 0; iCell < grid->GetNTotalCell(); ++ iCell)
            {
                phiOld[index][iCell] = phiNow[index][iCell];

            }
        }

        vector<RDouble *>().swap(phiNow);
        vector<RDouble *>().swap(phiOld);
    }


    DEFINE_PBSOLVER_REGISTER_FUNCTION(PBTransient, PBTransientImplicitEuler, TRANSIENT_IMPLICIT_EULER);
}