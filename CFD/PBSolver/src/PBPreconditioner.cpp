#include "PBPreconditioner.h"
#include "Task_ServerUpdateInterface.h"

namespace PHSPACE
{
	PBPreconditioner::PBPreconditioner(const std::string& name)
		:PBObject(name)
	{

	}

	PBPreconditioner::~PBPreconditioner()
	{

	}

	void PBPreconditioner::Preconditioning(Region& region, vector<int>& zoneId, string& varStr, string& apStr, string& ajUpperStr,
		                                   string& ajLowerStr, string& aiiStr, string& scStr, string& rsStr)
	{

	}

	void PBPreconditioner::PreconditioningProxy(int nfaces, int ncells, int* leftCellOfFace, int* rightCellOfFace,
		                                        RDouble* ap, RDouble* ajLeft, RDouble* ajRight, RDouble* rs, string* phiName)
	{
		
	}

	void PBPreconditioner::ForwardAndBackward(int nfaces, int ncells, int* leftCellOfFace, int* rightCellOfFace, RDouble* ap0,
		                                      RDouble* ajLeft, RDouble* ajRight, RDouble* p, RDouble* rs, string* phiName)
	{
		
	}

	DEFINE_PBSOLVER_REGISTER_FUNCTION(PBPreconditioner, PBPreconditioner, ITERSOLVPRECOND_EMPTY);
}