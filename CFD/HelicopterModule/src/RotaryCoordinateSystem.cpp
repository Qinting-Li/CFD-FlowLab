#include "RotaryCoordinateSystem.h"
#include "MBConstTable.h"
#include "Geo_UnstructGrid.h"
using namespace std;

namespace PHSPACE
{
RotaryCoordinateSystem::RotaryCoordinateSystem()
{
    omiga.resize(MBDIM);

    this->grid = NULL;
    zone_st    = 0;
    zone_ed    = 0;
}

RotaryCoordinateSystem::~RotaryCoordinateSystem()
{

}

void RotaryCoordinateSystem::RecognizeGrid(UnstructGrid * grid)
{
    int zoneIndex = grid->GetGridID()->GetIndex();
    if (zone_st <= zoneIndex && zoneIndex <= zone_ed)
    { 
        this->grid = grid;
    }
    return;
}

}