#include "Post_ProbesWrite.h"
#include "IO_FileName.h"
#include "TK_Parse.h"
#include "PHIO.h"
#include "HyList.h"
#include "TK_Exit.h"
#include "Glb_Dimension.h"
#include "BC_Unstruct.h"
#include "BC_Struct.h"
#include "Post_Probes.h"
#include "Pointer.h"
#include "Solver.h"
#include "Geo_UnstructGrid.h"
#include "TK_Log.h"

namespace PHSPACE
{

LIB_EXPORT Post_ProbesWrite::Post_ProbesWrite(int flowType)
{
    this->flowType = flowType;
}

LIB_EXPORT Post_ProbesWrite::~Post_ProbesWrite()
{

}

LIB_EXPORT void Post_ProbesWrite::Run(bool dumpAllLevel)
{
    bool run_now = true;

    int nMGLevel = 1;
    if (dumpAllLevel)
    {
        GlobalDataBase::GetData("nMGLevel", &nMGLevel, PHINT, 1);
    }

    for (int iLevel = 0; iLevel < nMGLevel; ++ iLevel)
    {
        string probesflowfile = "sample.dat";
        GlobalDataBase::GetData("probesflowfile", &probesflowfile, PHSTRING, 1);

        string fname, fext;
        GetNameExt(probesflowfile, fname, fext, ".");

        ActionKey *actkey = new ActionKey();
        actkey->level    = iLevel;
        actkey->filename = AddDigitToFileName(fname, fext, ".", actkey->level);
        actkey->openmode = ios_base::out|ios_base::app;

        int isUnsteady = GlobalDataBase::GetIntParaFromDB("iunsteady");

        if (isUnsteady != 0)
        {
            int outnstep = PHSPACE::GlobalDataBase::GetIntParaFromDB("outnstep");
            actkey->filename = PHSPACE::AddSymbolToFileName(actkey->filename, '_', outnstep);
        }

        ServerProbesWrite(actkey);
        PostProbesWrite(actkey);
        delete actkey;
        this->Clear();
    }
}

void Post_ProbesWrite::ServerProbesWrite(ActionKey *actkey)
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    int recv_proc;
    int myid;

    for (int izone =  0; izone < nZones; ++ izone)
    {
        DataContainer *cdata = actkey->GetData();
        cdata->moveto_begin();

        int send_proc = GetZoneProcessorID(izone);
        recv_proc = GetServerProcessorID();

        myid = GetCurrentProcessorID();
        int tag = GetSendRecvTag(actkey, izone);


        if (myid == send_proc)
        {
            Grid *grid = GetGrid(izone, actkey->level);
            int nProbeVariables = GlobalDataBase::GetIntParaFromDB("nProbeVariables");
            int probeVariables[100];

            GlobalDataBase::GetData("probeVariables", probeVariables, PHINT, nProbeVariables);
            Post_Probes * postProbesVar = new Post_Probes(nProbeVariables, probeVariables, grid);
            int Nsolver = GlobalSolvers::GetNumberOfSolvers(izone);
            for (int isolver = 0; isolver < Nsolver; ++ isolver)
            {
                GlobalSolvers::GetSolver(izone, isolver)->ComputePostProbesVariables(postProbesVar);
            }

            StoreProbesVarariables(actkey, postProbesVar);
            FreePointer(postProbesVar);
        }

        PH_Trade(actkey, send_proc, recv_proc, tag);

        if (myid == recv_proc)
        {
            datalist.push_back(cdata);
            DataContainer * cdata1 = new DataContainer();
            actkey->SetData(cdata1);
        }
    }
}

void Post_ProbesWrite::StoreProbesVarariables(ActionKey *actkey, Post_Probes *postProbesVar)
{
    Grid *grid_in = postProbesVar->GetGrid();

    int nProbeVariables = postProbesVar->GetNumberofProbeVariables();

    RDouble ** q_probes = new RDouble *[nProbeVariables];
    postProbesVar->GetAllProbeVarPtr(q_probes);

    if (actkey->filename != "")
    {
        int zoneProbesNumber = grid_in->GetZoneProbesNumber();
        PHWrite(actkey->GetData(), &zoneProbesNumber, 1);

        vector<vector<RDouble> > probesCoordinates = grid_in->GetProbesCoordinates();
        int nDim;

        for (int iProbe = 0; iProbe < zoneProbesNumber; ++ iProbe)
        {
            nDim = static_cast<int>(probesCoordinates[iProbe].size());

            PHWrite(actkey->GetData(), &nDim, 1);

            for (int iDim = 0; iDim < nDim; ++ iDim)
            {
                PHWrite(actkey->GetData(), probesCoordinates[iProbe][iDim]);
            }
            for (int iVariable = 0; iVariable < nProbeVariables; ++ iVariable)
            {
                PHWrite(actkey->GetData(), q_probes[iVariable][iProbe]);
            }
        }
        delete [] q_probes;
    }
}

void Post_ProbesWrite::PostProbesWrite(ActionKey *actkey)
{
    using namespace PHMPI;

    int myid = GetCurrentProcessorID();

    if (datalist.size() > 0)
    {
        DumpToFileASCII dumpToFileASCII(datalist, actkey->filename);
        dumpToFileASCII.Run();
    }

    for (unsigned int izone = 0; izone < datalist.size(); ++ izone)
    {
        delete datalist[izone];
    }
}

void Post_ProbesWrite::Clear()
{
    this->datalist.clear();
}

}