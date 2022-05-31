#include <iostream>

#include "PBBoundaryScalarConvection.h"

namespace PHSPACE
{
	PBBoundaryScalarConvection::PBBoundaryScalarConvection(const string& name)
		:PBBoundaryScalar(name)
	{
	}

	PBBoundaryScalarConvection::~PBBoundaryScalarConvection()
	{

	}

	void PBBoundaryScalarConvection::Discrt(UnstructGrid* grid, string phiName, vector<int>* faceIndex, Data_Param* bcData, RDouble* diffCoe)
	{
		int nTotalCells = grid->GetNTotalCell();
		RDouble* phi = reinterpret_cast<RDouble*>(grid->GetDataPtr(phiName));
		RDouble* phib = phi + nTotalCells;
		RDouble* diffCoeb = diffCoe + nTotalCells;
		RDouble* fluxb = reinterpret_cast<RDouble*>(grid->GetDataPtr("flux"));
		RDouble* apPhi = reinterpret_cast<RDouble*>(grid->GetDataPtr("ap"));
		RDouble* scPhi = reinterpret_cast<RDouble*>(grid->GetDataPtr("sc"));

		for (std::size_t iFace = 0; iFace < faceIndex->size(); ++iFace)
		{
			int iFacelocal = (*faceIndex)[iFace];
			int le = grid->GetLeftCellOfFace()[iFacelocal];
			RDouble fluxbValue = fluxb[iFacelocal];
			if (fluxbValue < 0.0)
			{
				apPhi[le] -= fluxbValue;
				scPhi[le] -= fluxbValue * phib[iFacelocal];
			}
		}
	}

	DEFINE_PBSOLVER_REGISTER_FUNCTION(PBBoundary, PBBoundaryScalarConvection, SCALAR_CONVECTION);
}