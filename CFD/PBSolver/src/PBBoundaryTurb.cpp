#include "PBBoundaryTurb.h"

using namespace std;

namespace PHSPACE
{
	PBBoundaryTurb::PBBoundaryTurb(const string& name)
		:PBBoundary(name)
	{

	}

	PBBoundaryTurb::~PBBoundaryTurb()
	{

	}

	DEFINE_PBSOLVER_REGISTER_FUNCTION(PBBoundary, PBBoundaryTurb, TURB_EMPTY);
} // PB_SPACE