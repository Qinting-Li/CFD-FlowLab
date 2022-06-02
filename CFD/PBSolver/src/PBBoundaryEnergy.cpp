#include "PBBoundaryEnergy.h"

using namespace std;

namespace PHSPACE
{
	PBBoundaryEnergy::PBBoundaryEnergy(const string& name)
		:PBBoundaryScalar(name)
	{
		
	}

	PBBoundaryEnergy::~PBBoundaryEnergy()
	{

	}

	DEFINE_PBSOLVER_REGISTER_FUNCTION(PBBoundary, PBBoundaryEnergy, ENERGY_EMPTY);
}