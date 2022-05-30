#include "PBGrad.h"
#include "Solver.h"
#include "Geo_UnstructGrid.h"
#include "Task_ServerUpdateInterface.h"

namespace PHSPACE
{
	PBGrad::PBGrad(const string& name)
		:PBObject(name)
	{

	}

	PBGrad::~PBGrad()
	{

	}

	DEFINE_PBSOLVER_REGISTER_FUNCTION(PBGrad, PBGrad, GRAD_EMPTY);
}