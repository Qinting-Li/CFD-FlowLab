#include <iostream>
#include <algorithm>

#include "PBBoundaryFlowVelocityInlet.h"
#include "Geo_UnstructGrid.h"

using namespace std;

namespace PHSPACE
{
	PBBoundaryFlowVelocityInlet::PBBoundaryFlowVelocityInlet(const string& name)
		:PBBoundaryFlow(name)
	{

	}

	PBBoundaryFlowVelocityInlet::~PBBoundaryFlowVelocityInlet()
	{

	}

	void PBBoundaryFlowVelocityInlet::CalcVelocityBoundaryCondition(UnstructGrid* grid, vector<int>* faceIndex)
	{
		int* leftCellOfFace  = grid->GetLeftCellOfFace();
		int* rightCellOfFace = grid->GetRightCellOfFace();
		RDouble* faceArea    = grid->GetFaceArea();
		RDouble* faceNormalX = grid->GetFaceNormalX();
		RDouble* faceNormalY = grid->GetFaceNormalY();
		RDouble* faceNormalZ = grid->GetFaceNormalZ();
		RDouble* faceCenterX = grid->GetFaceCenterX();
		RDouble* faceCenterY = grid->GetFaceCenterY();
		RDouble* faceCenterZ = grid->GetFaceCenterZ();
		RDouble* cellCenterX = grid->GetCellCenterX();
		RDouble* cellCenterY = grid->GetCellCenterY();
		RDouble* cellCenterZ = grid->GetCellCenterZ();

		RDouble* mu    = reinterpret_cast<RDouble*>(grid->GetDataPtr("mu"));
		RDouble* rho   = reinterpret_cast<RDouble*>(grid->GetDataPtr("rho"));
		RDouble* u     = reinterpret_cast<RDouble*>(grid->GetDataPtr("U"));
		RDouble* v     = reinterpret_cast<RDouble*>(grid->GetDataPtr("V"));
		RDouble* w     = reinterpret_cast<RDouble*>(grid->GetDataPtr("W"));
		RDouble* dudx  = reinterpret_cast<RDouble*>(grid->GetDataPtr("dudx"));
		RDouble* dudy  = reinterpret_cast<RDouble*>(grid->GetDataPtr("dudy"));
		RDouble* dudz  = reinterpret_cast<RDouble*>(grid->GetDataPtr("dudz"));
		RDouble* dvdx  = reinterpret_cast<RDouble*>(grid->GetDataPtr("dvdx"));
		RDouble* dvdy  = reinterpret_cast<RDouble*>(grid->GetDataPtr("dvdy"));
		RDouble* dvdz  = reinterpret_cast<RDouble*>(grid->GetDataPtr("dvdz"));
		RDouble* dwdx  = reinterpret_cast<RDouble*>(grid->GetDataPtr("dwdx"));
		RDouble* dwdy  = reinterpret_cast<RDouble*>(grid->GetDataPtr("dwdy"));
		RDouble* dwdz  = reinterpret_cast<RDouble*>(grid->GetDataPtr("dwdz"));
		RDouble* fluxb = reinterpret_cast<RDouble*>(grid->GetDataPtr("flux"));
		RDouble* apU   = reinterpret_cast<RDouble*>(grid->GetDataPtr("apU"));
		RDouble* apV   = reinterpret_cast<RDouble*>(grid->GetDataPtr("apV"));
		RDouble* apW   = reinterpret_cast<RDouble*>(grid->GetDataPtr("apW"));
		RDouble* scU   = reinterpret_cast<RDouble*>(grid->GetDataPtr("scU"));
		RDouble* scV   = reinterpret_cast<RDouble*>(grid->GetDataPtr("scV"));
		RDouble* scW   = reinterpret_cast<RDouble*>(grid->GetDataPtr("scW"));

		for (size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
		{
			int iFacelocal = (*faceIndex)[iFace];
			int le = leftCellOfFace[iFacelocal];
			int re = rightCellOfFace[iFacelocal];
			RDouble ub = u[re];
			RDouble vb = v[re];
			RDouble wb = w[re];
			RDouble muLocal = mu[le];
			RDouble area = faceArea[iFacelocal];

			RDouble fanx = faceNormalX[iFacelocal] * area;
			RDouble fany = faceNormalY[iFacelocal] * area;
			RDouble fanz = faceNormalZ[iFacelocal] * area;
			RDouble dx = faceCenterX[iFacelocal] - cellCenterX[le];
			RDouble dy = faceCenterY[iFacelocal] - cellCenterY[le];
			RDouble dz = faceCenterZ[iFacelocal] - cellCenterZ[le];

			RDouble fluxbValue = fluxb[iFacelocal];
			RDouble area2 = muLocal * (fanx * fanx + fany * fany + fanz * fanz) / (fanx * dx + fany * dy + fanz * dz);
			RDouble ax = muLocal * fanx - dx * area2;
			RDouble ay = muLocal * fany - dy * area2;
			RDouble az = muLocal * fanz - dz * area2;

			apU[le] = apU[le] + area2;
			apV[le] = apV[le] + area2;
			apW[le] = apW[le] + area2;
			scU[le] = scU[le] + area2 * ub + dudx[le] * ax + dudy[le] * ay + dudz[le] * az;
			scV[le] = scV[le] + area2 * vb + dvdx[le] * ax + dvdy[le] * ay + dvdz[le] * az;
			scW[le] = scW[le] + area2 * wb + dwdx[le] * ax + dwdy[le] * ay + dwdz[le] * az;

			if (fluxbValue < 0)
			{
				RDouble a12 = std::max(fluxbValue, 0.0);
				RDouble a21 = a12 - fluxbValue;
				apU[le] = apU[le] + a21;
				apV[le] = apV[le] + a21;
				apW[le] = apW[le] + a21;
				scU[le] = scU[le] + a21 * ub;
				scV[le] = scV[le] + a21 * vb;
				scW[le] = scW[le] + a21 * wb;
			}
		}
	}

	void PBBoundaryFlowVelocityInlet::UpdateVelocityBoundaryCondition(UnstructGrid* grid, vector<int>* faceIndex, Data_Param* bcData)
	{
		RDouble ub, vb, wb;
		int* leftCellOfFace = grid->GetLeftCellOfFace();
		int* rightCellOfFace = grid->GetRightCellOfFace();
		RDouble* u = reinterpret_cast<RDouble*>(grid->GetDataPtr("U"));
		RDouble* v = reinterpret_cast<RDouble*>(grid->GetDataPtr("V"));
		RDouble* w = reinterpret_cast<RDouble*>(grid->GetDataPtr("W"));
		for (size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
		{
			int iFacelocal = (*faceIndex)[iFace];

			int re = rightCellOfFace[iFacelocal];

			if (bcData)
			{
				if (bcData->IsExist("flowU", PHDOUBLE, 1))
				{
					bcData->GetData("flowU", &ub, PHDOUBLE, 1);
				}
				else
				{
					ub = 0.0;
				}

				if (bcData->IsExist("flowV", PHDOUBLE, 1))
				{
					bcData->GetData("flowV", &vb, PHDOUBLE, 1);
				}
				else
				{
					vb = 0.0;
				}

				if (bcData->IsExist("flowW", PHDOUBLE, 1))
				{
					bcData->GetData("flowW", &wb, PHDOUBLE, 1);
				}
				else
				{
					wb = 0.0;
				}

				u[re] = RDouble(ub);
				v[re] = RDouble(vb);
				w[re] = RDouble(wb);
			}
		}
	}


	void PBBoundaryFlowVelocityInlet::UpdatePressureBoundaryCondition(UnstructGrid* grid, vector<int>* faceIndex, Data_Param* bcData)
	{
		int* leftCellOfFace = grid->GetLeftCellOfFace();
		int* rightCellOfFace = grid->GetRightCellOfFace();
		RDouble* p = reinterpret_cast<RDouble*>(grid->GetDataPtr("P"));
		for (size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
		{
			int iFacelocal = (*faceIndex)[iFace];
			int le = leftCellOfFace[iFacelocal];
			int re = rightCellOfFace[iFacelocal];
			p[re] = p[le];
		}
	}

	void PBBoundaryFlowVelocityInlet::CalcPressureCorrectionBoundaryCondition(UnstructGrid* grid, vector<int>* faceIndex, Data_Param* bcData)
	{
		int* leftCellOfFace  = grid->GetLeftCellOfFace();
		int* rightCellOfFace = grid->GetRightCellOfFace();

		RDouble* faceArea    = grid->GetFaceArea();
		RDouble* faceNormalX = grid->GetFaceNormalX();
		RDouble* faceNormalY = grid->GetFaceNormalY();
		RDouble* faceNormalZ = grid->GetFaceNormalZ();
		RDouble* faceCenterX = grid->GetFaceCenterX();
		RDouble* faceCenterY = grid->GetFaceCenterY();
		RDouble* faceCenterZ = grid->GetFaceCenterZ();
		RDouble* cellCenterX = grid->GetCellCenterX();
		RDouble* cellCenterY = grid->GetCellCenterY();
		RDouble* cellCenterZ = grid->GetCellCenterZ();

		RDouble* mu   = reinterpret_cast<RDouble*>(grid->GetDataPtr("mu"));
		RDouble* rho  = reinterpret_cast<RDouble*>(grid->GetDataPtr("rho"));
		RDouble* dunb = reinterpret_cast<RDouble*>(grid->GetDataPtr("dun"));
		RDouble* dpdx = reinterpret_cast<RDouble*>(grid->GetDataPtr("dpdx"));
		RDouble* dpdy = reinterpret_cast<RDouble*>(grid->GetDataPtr("dpdy"));
		RDouble* dpdz = reinterpret_cast<RDouble*>(grid->GetDataPtr("dpdz"));
		RDouble* hu   = reinterpret_cast<RDouble*>(grid->GetDataPtr("hu"));
		RDouble* hv   = reinterpret_cast<RDouble*>(grid->GetDataPtr("hv"));
		RDouble* hw   = reinterpret_cast<RDouble*>(grid->GetDataPtr("hw"));
		RDouble* apU  = reinterpret_cast<RDouble*>(grid->GetDataPtr("apU"));
		RDouble* apV  = reinterpret_cast<RDouble*>(grid->GetDataPtr("apV"));
		RDouble* apW  = reinterpret_cast<RDouble*>(grid->GetDataPtr("apW"));
		RDouble urfFlux = ObjGetDoubleParaFromDB("urfFlux", 0.8);

		RDouble* u = reinterpret_cast<RDouble*>(grid->GetDataPtr("U"));
		RDouble* v = reinterpret_cast<RDouble*>(grid->GetDataPtr("V"));
		RDouble* w = reinterpret_cast<RDouble*>(grid->GetDataPtr("W"));
		RDouble* p = reinterpret_cast<RDouble*>(grid->GetDataPtr("P"));
		RDouble* fluxb = reinterpret_cast<RDouble*>(grid->GetDataPtr("flux"));
		RDouble* ap = reinterpret_cast<RDouble*>(grid->GetDataPtr("ap"));
		RDouble* sc = reinterpret_cast<RDouble*>(grid->GetDataPtr("sc"));

		for (size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
		{
			int iFacelocal = (*faceIndex)[iFace];
			int le = leftCellOfFace[iFacelocal];
			int re = rightCellOfFace[iFacelocal];

			RDouble ub = u[re];
			RDouble vb = v[re];
			RDouble wb = w[re];
			RDouble muLocal = mu[le];
			RDouble area = faceArea[iFacelocal];

			RDouble fanx = faceNormalX[iFacelocal] * area;
			RDouble fany = faceNormalY[iFacelocal] * area;
			RDouble fanz = faceNormalZ[iFacelocal] * area;

			RDouble fluxbValue = rho[le] * (ub * fanx + vb * fany + wb * fanz);
			fluxb[iFacelocal] = fluxbValue;
			sc[le] -= fluxbValue;
		}
	}

	void PBBoundaryFlowVelocityInlet::UpdatePressureCorrectionBoundaryCondition(UnstructGrid* grid, vector<int>* faceIndex, Data_Param* bcData)
	{
		int* leftCellOfFace  = grid->GetLeftCellOfFace();
		int* rightCellOfFace = grid->GetRightCellOfFace();
		RDouble* pc = reinterpret_cast<RDouble*>(grid->GetDataPtr("pc"));
		for (size_t i = 0; i < faceIndex->size(); ++i)
		{
			int iFacelocal = (*faceIndex)[i];
			int le = leftCellOfFace[iFacelocal];
			int re = rightCellOfFace[iFacelocal];
			pc[re] = pc[le];
		}
	}

	void PBBoundaryFlowVelocityInlet::UpdateVelocityCorrectionBoundary(UnstructGrid* grid, vector<int>* faceIndex, Data_Param* bcData)
	{

	}

	DEFINE_PBSOLVER_REGISTER_FUNCTION(PBBoundaryFlow, PBBoundaryFlowVelocityInlet, FLOW_VELOCITY_INLET);
}