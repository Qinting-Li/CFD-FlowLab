#include <iostream>
#include "PBBoundaryFlow.h"

using namespace std;

namespace PHSPACE
{
PBBoundaryFlow::PBBoundaryFlow(const string &name):PBBoundary(name)
{

}

PBBoundaryFlow::~PBBoundaryFlow()
{

}
void PBBoundaryFlow::CalcVelocityBoundaryCondition(UnstructGrid *grid, vector<int> *faceIndex)
{

}

void PBBoundaryFlow::UpdateVelocityBoundaryCondition(UnstructGrid *grid, vector<int> *faceIndex, Data_Param *bcData)
{

}

void PBBoundaryFlow::UpdatePressureBoundaryCondition(UnstructGrid *grid, vector<int> *faceIndex, Data_Param *bcData)
{

}

void PBBoundaryFlow::CalcPressureCorrectionBoundaryCondition(UnstructGrid *grid, vector<int> *faceIndex, Data_Param *bcData)
{

}

void PBBoundaryFlow::UpdatePressureCorrectionBoundaryCondition(UnstructGrid *grid, vector<int> *faceIndex, Data_Param *bcData)
{

}

void PBBoundaryFlow::UpdateVelocityCorrectionBoundary(UnstructGrid *grid, vector<int> *faceIndex, Data_Param *bcData)
{

}

}