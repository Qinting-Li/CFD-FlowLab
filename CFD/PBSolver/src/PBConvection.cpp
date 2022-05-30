#include "PBConvection.h"

namespace PHSPACE
{
	PBConvection::PBConvection(const string& name)
		:PBObject(name)
	{

	}

	PBConvection::~PBConvection()
	{

	}

	DEFINE_PBSOLVER_REGISTER_FUNCTION(PBConvection, PBConvection, CONVECTION_EMPTY);
}