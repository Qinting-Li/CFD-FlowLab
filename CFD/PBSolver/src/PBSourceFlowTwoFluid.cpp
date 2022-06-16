#include <iostream>
#include "PBSourceFlowTwoFluid.h"
#include "PBGrad.h"
#include "GlobalDataBase.h"

using namespace std;

namespace PHSPACE
{
    PBSourceFlowTwoFluid::PBSourceFlowTwoFluid(const std::string &name):PBSourceFlowDefault(name), PBPhaseBase(),
        number_(0), firstPhaseIndex_(nullptr), secondPhaseIndex_(nullptr), diameter_(nullptr), eog_(nullptr)
    {
        this->GetPhasesParams();

        this->number_ = GlobalDataBase::GetIntParaFromDB("twoFluidMomNumber");

        this->firstPhaseIndex_ = new int[this->number_];
        this->secondPhaseIndex_ = new int[this->number_];
        this->diameter_ = new RDouble[this->number_];
        this->eog_ = new RDouble[this->number_];

        GlobalDataBase::GetData("twoFluidMomFirstPhaseIndex", this->firstPhaseIndex_, PHINT, this->number_);
        GlobalDataBase::GetData("twoFluidMomSecondPhaseIndex", this->secondPhaseIndex_, PHINT, this->number_);
        GlobalDataBase::GetData("twoFluidMomDiameter", this->diameter_, PHDOUBLE, this->number_);

        for (int i = 0; i < this->number_; ++ i)
        {
            this->eog_[i] = 9.81 * (998.2 - 1.225) * std::pow(this->diameter_[i], 2) / 0.073;
        }
    }

    PBSourceFlowTwoFluid::~PBSourceFlowTwoFluid()
    {

    }

    void PBSourceFlowTwoFluid::Discrt(const string phiName, UnstructGrid &grid)
    {
        int ip = GlobalDataBase::GetIntParaFromDB("twoFluidPhaseIndexNow");
        const RDouble *dpdx = (RDouble*)grid.GetDataPtr("dpdx");
        const RDouble *dpdy = (RDouble*)grid.GetDataPtr("dpdy");
        const RDouble *dpdz = (RDouble*)grid.GetDataPtr("dpdz");

        this->DiscrtProxy(
            grid.GetNTotalCell(),
            grid.GetCellVolume(),
            dpdx,
            dpdy,
            dpdz,
            (RDouble *)grid.GetDataPtr("scU" + std::to_string(ip)),
            (RDouble *)grid.GetDataPtr("scV" + std::to_string(ip)),
            (RDouble *)grid.GetDataPtr("scW" + std::to_string(ip))
            );
    }

    void PBSourceFlowTwoFluid::DiscrtProxy
        (
        const int totalCells,
        const RDouble *cellVolume,
        const RDouble *dpdx,
        const RDouble *dpdy,
        const RDouble *dpdz,
        RDouble *scU,
        RDouble *scV,
        RDouble *scW
        )
    {
        int ip = GlobalDataBase::GetIntParaFromDB("twoFluidPhaseIndexNow");
        RDouble *vff = (RDouble *)(GetGrid(0)->GetDataPtr(this->phaseNames[ip]));
        for (int ic = 0; ic < totalCells; ++ic)
        {
            const RDouble cv = cellVolume[ic];
            scU[ic] -= vff[ic] * dpdx[ic] * cv;
            scV[ic] -= vff[ic] * dpdy[ic] * cv;
            scW[ic] -= vff[ic] * dpdz[ic] * cv;
        }
        if (ip == 0)
        {
            int nTotalCells = GetGrid(0)->GetNTotalCell();
            RDouble *src_mom = (RDouble*)(GetGrid(0)->GetDataPtr("srcMom"));
            RDouble *u1 = (RDouble *)(GetGrid(0)->GetDataPtr("U" + to_string(ip)));
            RDouble *u2 = (RDouble *)(GetGrid(0)->GetDataPtr("U" + to_string(ip + 1)));
            RDouble *v1 = (RDouble *)(GetGrid(0)->GetDataPtr("V" + to_string(ip)));
            RDouble *v2 = (RDouble *)(GetGrid(0)->GetDataPtr("V" + to_string(ip + 1)));
            RDouble *w1 = (RDouble *)(GetGrid(0)->GetDataPtr("W" + to_string(ip)));
            RDouble *w2 = (RDouble *)(GetGrid(0)->GetDataPtr("W" + to_string(ip + 1)));
            RDouble *apU1 = (RDouble *)(GetGrid(0)->GetDataPtr("apU" + to_string(ip)));
            RDouble *apV1 = (RDouble *)(GetGrid(0)->GetDataPtr("apV" + to_string(ip)));
            RDouble *apW1 = (RDouble *)(GetGrid(0)->GetDataPtr("apW" + to_string(ip)));
            RDouble *apU2 = (RDouble *)(GetGrid(0)->GetDataPtr("apU" + to_string(ip + 1)));
            RDouble *apV2 = (RDouble *)(GetGrid(0)->GetDataPtr("apV" + to_string(ip + 1)));
            RDouble *apW2 = (RDouble *)(GetGrid(0)->GetDataPtr("apW" + to_string(ip + 1)));
            RDouble *scU2 = (RDouble *)(GetGrid(0)->GetDataPtr("scU" + to_string(ip + 1)));
            RDouble *scV2 = (RDouble *)(GetGrid(0)->GetDataPtr("scV" + to_string(ip + 1)));
            RDouble *scW2 = (RDouble *)(GetGrid(0)->GetDataPtr("scW" + to_string(ip + 1)));

            for (int ic = 0; ic < totalCells; ++ ic)
            {
                apU1[ic] = apU1[ic] + src_mom[ic];
                apV1[ic] = apV1[ic] + src_mom[ic];
                apW1[ic] = apW1[ic] + src_mom[ic];
                apU2[ic] = apU2[ic] + src_mom[ic];
                apV2[ic] = apV2[ic] + src_mom[ic];
                apW2[ic] = apW2[ic] + src_mom[ic];

                scU[ic] = scU[ic] + src_mom[ic] * u2[ic];
                scV[ic] = scV[ic] + src_mom[ic] * v2[ic];
                scW[ic] = scW[ic] + src_mom[ic] * w2[ic];
                scU2[ic] = scU2[ic] + src_mom[ic] * u1[ic];
                scV2[ic] = scV2[ic] + src_mom[ic] * v1[ic];
                scW2[ic] = scW2[ic] + src_mom[ic] * w1[ic];
            }
        }
    }

    DEFINE_PBSOLVER_REGISTER_FUNCTION(PBSource, PBSourceFlowTwoFluid, SOURCE_TWO_FLUID_FLOW);
}