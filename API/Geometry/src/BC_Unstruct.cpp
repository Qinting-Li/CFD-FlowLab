#include "BC_Unstruct.h"
#include "TK_Exit.h"

namespace PHSPACE
{
BoundaryCondition * BCRecord::GetBoundaryCondition() const 
{
	if(!this->boundaryCondition)
	{
		TK_Exit::ExceptionExit("Error: boundary file of .bc does not exit!");
		return this->boundaryCondition;
	}
	else
	{
		return this->boundaryCondition;
	}
}

}
