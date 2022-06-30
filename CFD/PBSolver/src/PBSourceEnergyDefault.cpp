#include <iostream>
#include "PBSourceEnergyDefault.h"

using namespace std;

namespace PHSPACE
{
PBSourceEnergyDefault::PBSourceEnergyDefault(const string &name):PBSource(name)
{

}

PBSourceEnergyDefault::~PBSourceEnergyDefault()
{

}

void PBSourceEnergyDefault::Discrt(UnstructGrid *grid,string &phiName)
{
    DiscrtProxy(grid);
}

void PBSourceEnergyDefault::DiscrtProxy(UnstructGrid *grid)
{

}

DEFINE_PBSOLVER_REGISTER_FUNCTION(PBSource, PBSourceEnergyDefault, SOURCE_ENERGY_DEFAULT);
}