#include "Task_Client2Server.h"

namespace PHSPACE
{
Task_Client2Server::Task_Client2Server() {};

Task_Client2Server::~Task_Client2Server(){};

void Task_Client2Server::PreTask(ActionKey *actkey)
{
    using namespace PHMPI;

    int myid = GetCurrentProcessorID();

    int server = GetServerSepMode();
 
    if (myid == server)
    {
        ActionKeyFUN pre_action = this->GetPreAction();
        pre_action(actkey);
    }
}

void Task_Client2Server::PostTask(ActionKey *actkey)
{
    using namespace PHMPI;

    int myid = GetCurrentProcessorID();

    int server = GetServerSepMode();

    if (myid == server)
    {
        ActionKeyFUN post_action = this->GetPostAction();
        post_action(actkey);
    }
}

void Task_Client2Server::MainTask(ActionKey *actkey)
{
    using namespace PHMPI;

    int nZones = GetNumberofGlobalZones();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        Client2Server(actkey, iZone);
    }
}

void Task_Client2Server::Client2Server(ActionKey *actkey, int iZone)
{
    using namespace PHMPI;

    int send_proc = GetZoneProcessorIDSepMode(iZone);
    int recv_proc = GetServerSepMode();

    PH_WritePattern(actkey, GetMainAction(), send_proc, recv_proc, iZone);
}

}