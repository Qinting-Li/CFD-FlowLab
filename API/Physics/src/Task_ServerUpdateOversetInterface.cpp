#include "Task_ServerUpdateOversetInterface.h"
#include "Geo_Grid.h"
#include "OversetInformation.h"

namespace PHSPACE
{

Task_ServerUpdateOversetInterface::Task_ServerUpdateOversetInterface(){};

Task_ServerUpdateOversetInterface::~Task_ServerUpdateOversetInterface(){};

void Task_ServerUpdateOversetInterface::PreTask(ActionKey *actkey)
{

}

void Task_ServerUpdateOversetInterface::MainTask(ActionKey *actkey)
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        OversetTopologyManager * oversetTopologyManager = PHSPACE::GetOversetTopologyManager();

        uint_t numberOfOversetNeighbors = oversetTopologyManager->GetNumberOfOversetNeighbors(iZone);

        for (int iNeighbor = 0; iNeighbor < numberOfOversetNeighbors; ++ iNeighbor)
        {
            int oversetNeighborZoneIndex = oversetTopologyManager->GetNeighborOversetZoneIndex(iZone, iNeighbor);
            PH_Interface(actkey, iZone, oversetNeighborZoneIndex);
        }
    }
}

void Task_ServerUpdateOversetInterface::PostTask(ActionKey *actkey)
{

}
}