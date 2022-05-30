#include <iostream>

#include "PBBoundaryScalarDefault.h"

namespace PHSPACE
{
	PBBoundaryScalarDefault::PBBoundaryScalarDefault(const string& name)
		:PBBoundaryScalar(name),
		convCalc(nullptr),
		diffCalc(nullptr)
	{
		convCalc = new PBBoundaryScalarConvection(name);
		diffCalc = new PBBoundaryScalarDiffusion(name);
	}

	PBBoundaryScalarDefault::~PBBoundaryScalarDefault()
	{
		if (convCalc != nullptr)
		{
			delete convCalc;
		}

		if (diffCalc != nullptr)
		{
			delete diffCalc;
		}
	}

	void PBBoundaryScalarDefault::Discrt(UnstructGrid* grid, string phiName, vector<int>* faceIndex, Data_Param* bcData, RDouble* diffCoe)
	{
		convCalc->Discrt(grid, phiName,  faceIndex, bcData, diffCoe);

		diffCalc->Discrt(grid, phiName,  faceIndex, bcData, diffCoe);
	}

	DEFINE_PBSOLVER_REGISTER_FUNCTION(PBBoundary, PBBoundaryScalarDefault, SCALAR_DEFAULT);
}