#include "PostProcess.h"
#include "ForceProcessStruct.h"
#include "GridType.h"
#include "PHIO.h"
#include <cmath>

using namespace std;
namespace PHSPACE
{

PostProcess::PostProcess(Region *region_in)
{
    gridtype = 0;
    this->region = region_in;
    data = new PostProcessData();
}

PostProcess::~PostProcess()
{
    delete data;
}

void PostProcess::Run()
{
    ReadGridFile();
    ComputeMetrics();
    ReadFlowFile();
    InitGlobalValues();
    ReadBCFile();
    CompAirForceCoef();

}

void PostProcess::ReadGridFile()
{
    cout<<endl;
    cout<<"Begin to read grid finish ..."<<endl<<endl;
    region->ReadGrid();
    cout<<endl<<"Read grid finish ..."<<endl;

    Grid *grid = GetGrid(0);
    gridtype = grid->Type();
}

void PostProcess::ComputeMetrics()
{
    /*ActionKey *actkey = new ActionKey();
    actkey->action = COMPUTEMETRICS;
    actkey->kind   = GRID_BASED;
    actkey->level  = 0;

    Task *task = new Task_Client2Server();
    task->SetMainAction(WriteScreen);
    task->SetActionKey(actkey);
    PHSPACE::RunSimuTask(task, true);
    delete task;*/
}

void PostProcess::ReadFlowFile()
{
    IO_FlowReader *flowReader = new IO_FlowReader(this->region);
    flowReader->Run();

    delete flowReader;
}

void PostProcess::ReadBCFile()
{
    GlobalBoundaryCondition::ReadGlobalBoundaryCondition();

    GlobalBoundaryCondition::SetGlobalDataBaseByGlobalBC();
}

void PostProcess::InitGlobalValues()
{
    ComputeCoefficientOfStateEquation();
    ComputeReferenceTsuth();
    ComputeReferenceLengthDimensional();

    int nm    = GlobalDataBase::GetIntParaFromDB("nm");
    int nchem = GlobalDataBase::GetIntParaFromDB("nchem");
    int neqn = nm + nchem;

    RDouble *prim_inf = new RDouble[neqn];
    GlobalDataBase::UpdateDataPtr("prim_inf",prim_inf);

    RDouble attackd = GlobalDataBase::GetDoubleParaFromDB("attackd");
    RDouble angleSlide = GlobalDataBase::GetDoubleParaFromDB("angleSlide");
    RDouble coefficientOfStateEquation = GlobalDataBase::GetDoubleParaFromDB("coefficientOfStateEquation");

    RDouble attack   = attackd * PI / 180.0;
    RDouble sideslip = angleSlide * PI / 180.0;

    RDouble reference_density     = 1.0;
    RDouble reference_temperature = 1.0;

    RDouble reference_pressure, massoo_ref;
    massoo_ref = 1.0;

    reference_pressure = coefficientOfStateEquation * reference_density * reference_temperature / massoo_ref;

    using namespace IDX;

    prim_inf[IR] = reference_density;
    prim_inf[IU] = 1.0 * cos(attack) * cos(sideslip);
    prim_inf[IV] = 1.0 * sin(attack) * cos(sideslip);
    prim_inf[IW] = 1.0 *               sin(sideslip);
    prim_inf[IP] = reference_pressure;


}

void PostProcess::ComputeCoefficientOfStateEquation()
{
    RDouble refGama = GlobalDataBase::GetDoubleParaFromDB("refGama");

    RDouble refMachNumber = GlobalDataBase::GetDoubleParaFromDB("refMachNumber");

    RDouble coefficientOfStateEquation = 1.0 / (refGama * refMachNumber * refMachNumber);
    GlobalDataBase::UpdateData("coefficientOfStateEquation", &coefficientOfStateEquation, PHDOUBLE, 1);
}

void PostProcess::ComputeReferenceTsuth()
{
    RDouble tsi  = 110.4;

    RDouble refDimensionalTemperature = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");

    RDouble tsuth = tsi / refDimensionalTemperature;
    GlobalDataBase::UpdateData("tsuth", &tsuth, PHDOUBLE, 1);
}

void PostProcess::ComputeReferenceLengthDimensional()
{
	RDouble reynoldsReferenceLengthDimensional = 1.0;

	RDouble gridScaleFactor = GlobalDataBase::GetDoubleParaFromDB("gridScaleFactor");

	reynoldsReferenceLengthDimensional = gridScaleFactor;

	GlobalDataBase::UpdateData("reynoldsReferenceLengthDimensional", &reynoldsReferenceLengthDimensional, PHDOUBLE, 1);
}

void PostProcess::CompAirForceCoef()
{

    if ( gridtype == STRUCTGRID )
    {
        ForceProcessStruct *forceprocess = new ForceProcessStruct(this->data);
        forceprocess->Run();
    
        delete forceprocess;
    }
    else
    {




    }

    DumpAirForceCoef();


}

void PostProcess::DumpAirForceCoef()
{
    ostringstream oss;

    string header = "# iter		CL		Cd		Side		Cd_p		Cd_f\
    Cd-CL2par	Xcp		cfx		cfy		cfz		Cmx		Cmy		Cmz";

    int wordwidth = 16;
    oss << setiosflags(ios::left);
    oss << setiosflags(ios::scientific);
    oss << setprecision(4);

    int outnstep = GlobalDataBase::GetIntParaFromDB("outnstep");

    oss << header << endl;

    const int ALL = 0;

    RDouble *cfx            = data->GetCfx();
    RDouble *cfy            = data->GetCfy();
    RDouble *cfz            = data->GetCfz();
    RDouble *cpx            = data->GetCpx();
    RDouble *cpy            = data->GetCpy();
    RDouble *cpz            = data->GetCpz();
    RDouble *cmx            = data->GetCmx();
    RDouble *cmy            = data->GetCmy();
    RDouble *cmz            = data->GetCmz();
    RDouble *cl_tot         = data->GetClTot();
    RDouble *cd_tot         = data->GetCdTot();
    RDouble *cd_pr          = data->GetCdPr();
    RDouble *cd_sf          = data->GetCdSf();
    RDouble *cd_cl2pa       = data->GetCdCl2Pa();
    RDouble *xcp            = data->GetXcp();
    RDouble *side           = data->GetSide();
    RDouble *hingeMoment    = data->GetHingeMoment();

    oss << outnstep << "	";
    oss << cl_tot[ALL] << "	";
    oss << cd_tot[ALL] << "	";
    oss << side[ALL]   << "	";
    oss << cd_pr[ALL] << "	";
    oss << cd_sf[ALL] << "	";
    oss << cd_cl2pa[ALL] << "	";
    oss << xcp[ALL] << "	";
    oss << cfx[ALL] << "	";
    oss << cfy[ALL] << "	";
    oss << cfz[ALL] << "	";
    oss << cmx[ALL] << "	";
    oss << cmy[ALL] << "	";
    oss << cmz[ALL] << "	";

    oss << "\n";

    string filename = "results/AirForceCoef.dat";
    ios_base::openmode openmode = ios_base::out|ios_base::trunc;
    fstream file;
    ParallelOpenFile(file, filename, openmode);

    WriteASCIIFile(file, oss.str());

    ParallelCloseFile(file);
}

}