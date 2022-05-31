#include <iostream>

#include "PBBoundaryScalarDefault.h"

namespace PHSPACE
{
	PBBoundaryScalarDefault::PBBoundaryScalarDefault(const string& name)
		:PBBoundaryScalar(name),
		convCalc(NULL),
		diffCalc(NULL)
	{
		convCalc = new PBBoundaryScalarConvection(name);
		diffCalc = new PBBoundaryScalarDiffusion(name);
	}

	PBBoundaryScalarDefault::~PBBoundaryScalarDefault()
	{
		if (convCalc != NULL)
		{
			delete convCalc;
		}

		if (diffCalc != NULL)
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