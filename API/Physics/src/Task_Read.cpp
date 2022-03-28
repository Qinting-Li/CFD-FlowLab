#include "Task_Read.h"
#include "PHIO.h"

namespace PHSPACE
{
Task_Read::Task_Read() {};

Task_Read::~Task_Read(){};

void Task_Read::PreTask(ActionKey *actkey){}

void Task_Read::PostTask(ActionKey *actkey){}

void Task_Read::MainTask(ActionKey *actkey)
{
    using namespace PHMPI;
    //! Server mode.
    if (! PHMPI::IsParallelRun())
    {
        ServerRead(actkey);
    }
    else
    {
        ThreadRead(actkey);
    }
}


void Task_Read::ServerRead(ActionKey *actkey)
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    ParallelOpenFile(actkey);

    for (int izone = 0; izone < nZones; ++ izone)
    {
        int send_proc = server;
        int recv_proc = GetZoneProcessorID(izone);

        PH_ReadPattern(actkey, GetMainAction(), send_proc, recv_proc, izone);
    }

    ParallelCloseFile(actkey);
}

void Task_Read::ThreadRead(ActionKey *actkey)
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    ParallelOpenFile2(actkey);

    for (int izone = 0; izone < nZones; ++ izone)
    {
        //int send_proc = GetZoneProcessorID(izone);
        int send_proc = GetZoneFileID(izone);
        int recv_proc = GetZoneProcessorID(izone);

        PH_ReadPattern(actkey, GetMainAction(), send_proc, recv_proc, izone);
    }
    ParallelCloseFile(actkey);
}
}