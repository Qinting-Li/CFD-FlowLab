#include "PBBoundaryTurbEpsilonInlet.h"

using namespace std;

namespace PHSPACE
{
	PBBoundaryTurbEpsilonInlet::PBBoundaryTurbEpsilonInlet(const string& name)
		:PBBoundaryTurbDefault(name), cmu(0.09)
	{

	}

	PBBoundaryTurbEpsilonInlet::~PBBoundaryTurbEpsilonInlet()
	{

	}

	void PBBoundaryTurbEpsilonInlet::Update(UnstructGrid* grid, string phiName, vector<int>* faceIndex, Data_Param* bcData)
	{
		int nTotalCells = grid->GetNTotalCell();
		RDouble epsilonb = 0.0;
		RDouble* phi = reinterpret_cast<RDouble*>(grid->GetDataPtr("epsilon"));
		RDouble* phib = phi + nTotalCells;
		int* leftCellOfFace = grid->GetLeftCellOfFace();

		for (std::size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
		{
			const int iFacelocal = (*faceIndex)[iFace];

			if (bcData)
			{
				if (bcData->IsExist(phiName, PHDOUBLE, 1))
				{
					bcData->GetData(phiName, &epsilonb, PHDOUBLE, 1);
				}
			}
			else
			{
				std::cout << "Please assign a value of epsilon on inlet" << std::endl;;
			}
			phib[iFacelocal] = epsilonb;
		}
	}

	void PBBoundaryTurbEpsilonInlet::Discrt(UnstructGrid* grid, string phiName, vector<int>* faceIndex, Data_Param* bcData, RDouble* diffCoe)
	{
		int nTotalCells = grid->GetNTotalCell();
		RDouble* phi = reinterpret_cast<RDouble*>(grid->GetDataPtr(phiName));
		RDouble* phib = phi + nTotalCells;
		RDouble* diffCoeb = diffCoe + nTotalCells;
		RDouble* apPhi = reinterpret_cast<RDouble*>(grid->GetDataPtr("ap"));
		RDouble* scPhi = reinterpret_cast<RDouble*>(grid->GetDataPtr("sc"));
		RDouble* dphidx = reinterpret_cast<RDouble*>(grid->GetDataPtr("dphidx"));
		RDouble* dphidy = reinterpret_cast<RDouble*>(grid->GetDataPtr("dphidy"));
		RDouble* dphidz = reinterpret_cast<RDouble*>(grid->GetDataPtr("dphidz"));
		RDouble* fluxb = reinterpret_cast<RDouble*>(grid->GetDataPtr("flux"));

		for (std::size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
		{
			if (bcData)
			{
				//if (bcData->IsExist(phiName, PHDOUBLE, 1))
				//{
				//	int iFacelocal = (*faceIndex)[iFace]; //! The boundary element serial number is arranged in front of the element face.
				//	int le = grid->GetLeftCellOfFace()[iFacelocal];
				//	//! The boundary value is ranked behind the value on the cell. When it is passed in, the lower offset corresponds to the boundary sequence number.
				//	RDouble phibValue = phib[iFacelocal]; 
				//	RDouble diffCoebValue = diffCoeb[iFacelocal];
				//	RDouble area = grid->GetFaceArea()[iFacelocal];

				//	RDouble fanx = grid->GetFaceNormalX()[iFacelocal] * area;
				//	RDouble fany = grid->GetFaceNormalY()[iFacelocal] * area;
				//	RDouble fanz = grid->GetFaceNormalZ()[iFacelocal] * area;
				//	RDouble dx = grid->GetFaceCenterX()[iFacelocal] - grid->GetCellCenterX()[le];
				//	RDouble dy = grid->GetFaceCenterY()[iFacelocal] - grid->GetCellCenterY()[le];
				//	RDouble dz = grid->GetFaceCenterZ()[iFacelocal] - grid->GetCellCenterZ()[le];

				//	RDouble area2 = diffCoebValue * (fanx * fanx + fany * fany + fanz * fanz) / (fanx * dx + fany * dy + fanz * dz);
				//	RDouble ax = diffCoebValue * fanx - dx * area2;
				//	RDouble ay = diffCoebValue * fany - dy * area2;
				//	RDouble az = diffCoebValue * fanz - dz * area2;

				//	apPhi[le] += area2;
				//	scPhi[le] += area2 * phib[iFacelocal] + dphidx[le] * ax + dphidy[le] * ay + dphidz[le] * az;

				//	if (fluxb[iFacelocal] < 0)
				//	{
				//		apPhi[le] -= fluxb[iFacelocal];
				//		scPhi[le] -= fluxb[iFacelocal] * phib[iFacelocal];
				//	}
				//}
				//else
				//{
					int iFacelocal = (*faceIndex)[iFace]; //! The boundary element serial number is arranged in front of the element face.
					int le = grid->GetLeftCellOfFace()[iFacelocal];

					if (fluxb[iFacelocal] < 0)
					{
						apPhi[le] -= fluxb[iFacelocal];
						scPhi[le] -= fluxb[iFacelocal] * phib[iFacelocal];
					}
				/*}*/
			}
		}
	}

	DEFINE_PBSOLVER_REGISTER_FUNCTION(PBBoundary, PBBoundaryTurbEpsilonInlet, TURB_EPSILON_INLET);
} // PB_SPACE