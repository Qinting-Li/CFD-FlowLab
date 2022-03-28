#include "Task_Write.h"

#include "PHIO.h"
namespace PHSPACE
{
Task_Write::Task_Write(){}

Task_Write::~Task_Write(){}

void Task_Write::PreTask(ActionKey *actkey){}

void Task_Write::MainTask(ActionKey *actkey)
{
    //using namespace PHMPI;
    //! Server mode.
    if (! PHMPI::IsParallelRun())
    {
        ServerWrite(actkey);
    }
    else
    {
        ThreadWrite(actkey);
    }
}

void Task_Write::ServerWrite(ActionKey *actkey)
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    ParallelOpenFile(actkey);

    for (int izone = 0; izone < nZones; ++ izone)
    {
        int send_proc = GetZoneProcessorID(izone);
        int recv_proc = GetZoneFileID(izone);

        PH_WritePattern(actkey, GetMainAction(), send_proc, recv_proc, izone);
    }

    ParallelCloseFile(actkey);
}

void Task_Write::ThreadWrite(ActionKey *actkey)
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    ParallelOpenFile2(actkey);

    for (int izone = 0; izone < nZones; ++ izone)
    {
        int send_proc = GetZoneProcessorID(izone);
        int recv_proc = GetZoneFileID(izone);

        PH_WritePattern(actkey, GetMainAction(), send_proc, recv_proc, izone);
    }

    ParallelCloseFile(actkey);
}

void Task_Write::PostTask(ActionKey *actkey)
{
    using namespace PHMPI;
    int myid = GetCurrentProcessorID();
    if (myid == server)
    {
        ActionKeyFUN post_action = this->GetPostAction();
        post_action(actkey);
    }

}
}