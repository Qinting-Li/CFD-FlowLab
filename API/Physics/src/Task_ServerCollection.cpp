#include "Task_ServerCollection.h"

namespace PHSPACE
{
Task_ServerCollection::Task_ServerCollection() {};

Task_ServerCollection::~Task_ServerCollection(){};

void Task_ServerCollection::PreTask(ActionKey *actkey)
{
    ActionKeyFUN pre_action = this->GetPreAction();
    pre_action(actkey);
}

void Task_ServerCollection::PostTask(ActionKey *actkey)
{
    using namespace PHMPI;
    int myid = GetCurrentProcessorID();
    if (myid == server)
    {
        ActionKeyFUN post_action = this->GetPostAction();
        post_action(actkey);
    }
}

void Task_ServerCollection::MainTask(ActionKey *actkey)
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();
    actkey->GetData()->moveto_begin();
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        LocalCollection(actkey, iZone);
    }

    this->GetMainAction()(actkey);
}

void Task_ServerCollection::LocalCollection(ActionKey *actkey, int iZone)
{
    using namespace PHMPI;

    int zone_proc = GetZoneProcessorID(iZone);
    int myid      = GetCurrentProcessorID();

    if (myid == zone_proc)
    {
        //! If the process zone i belongs to is the current process.
        GeneralAction(actkey, iZone);
    }
}
}