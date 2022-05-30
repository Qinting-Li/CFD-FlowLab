#include "PBTransient.h"

namespace PHSPACE
{
	PBTransient::PBTransient(const string& name)
		:PBObject(name)
	{
		if (name == "EMPTY")
		{
			if (ObjGetIntParaFromDB("iunsteady", 0))
			{
				cout << "warning" << endl;
			}
		}
	}

	PBTransient::~PBTransient()
	{
		
	}

	DEFINE_PBSOLVER_REGISTER_FUNCTION(PBTransient, PBTransient, TRANSIENT_EMPTY);
}