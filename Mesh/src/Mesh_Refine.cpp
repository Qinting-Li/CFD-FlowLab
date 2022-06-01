#include "Mesh_Refine.h"
#include "TK_Log.h"
#include "Pre_HDF5File.h"

namespace PHSPACE
{

Mesh_Refine::Mesh_Refine(const string &gridFileNameIn)
{
    this->gridFileName = gridFileNameIn;
    numberOfZones = 0;

    region = 0;
    refinedGrids = 0;
    refineParameter = 0;
}

Mesh_Refine::~Mesh_Refine()
{
    delete region;

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        delete refinedGrids[iZone];
    }
    delete [] refinedGrids;
}

void Mesh_Refine::Run()
{
    ReadGrid();

    AllocateMemory();

    ConstructGridTopo();

    BuildRefineProperty();

    FixAnisoRefineType();

    RefineGrid();

#ifdef USE_INHOUSE_VERSION
#ifdef USE_OCC
    ShapePreserving();
#endif
#endif

    GenerateAndDumpComputationalGrid();
}

void Mesh_Refine::SetRefineParameter(RefineParameter *refineParameterIn)
{
    this->refineParameter = refineParameterIn;
}

void Mesh_Refine::ReadGrid()
{
    GlobalDataBase::UpdateData("gridfile", &gridFileName, PHSTRING, 1);

    region = new Region();
    region->ReadGrid();
    region->SetupMultiBlockRelationship();
}

void Mesh_Refine::DumpComputationalGrid()
{
    string out_gfile;
    GlobalDataBase::GetData("out_gfile", &out_gfile, PHSTRING, 1);

    WriteLogFile(" Grid Writing ...");
    if (PHMPI::GetCurrentProcessorID() == 0) cout << "  Grid Writing ... " << endl;

    if (out_gfile == "") return;

    PrintToWindow("    Writing new HDF5 version grid to file ....\n");
    WriteLogFile("    Writing new HDF5 version grid to file ....\n");

    IO_HDF5Write writeGrid(out_gfile, refinedGrids, numberOfZones);
    writeGrid.Run();

    WriteLogFile(" End Writing Grid!");
}


void RefineParameter::SetAnisoRefineType(int anisoRefineIn)
{
    this->anisoRefine = anisoRefineIn;
}

int RefineParameter::GetAnisoRefineType()
{
    return this->anisoRefine;
}

#ifdef USE_INHOUSE_VERSION
#ifdef USE_OCC
void RefineParameter::SetIsProject(int isProjectIn)
{
    this->isProject = isProjectIn;
}

void RefineParameter::SetIsDeform(int isDeformIn)
{
    this->isDeform = isDeformIn;
}

void RefineParameter::SetExclusiveCase(int exclusiveCaseIn)
{
    this->exclusiveCase = exclusiveCaseIn;
}

void RefineParameter::SetProjectOrgPoint(int projectOrgPointIn)
{
    this->projectOrgPoint = projectOrgPointIn;
}

void RefineParameter::SetGeometryUnit(int geometryUnitIn)
{
    this->geometryUnit = geometryUnitIn;
}

void RefineParameter::SetGeometryFileName(int geometryFileNameIn)
{
    this->geometryFileName = geometryFileNameIn;
}

int RefineParameter::GetIsProject()
{
    return this->isProject;
}

int RefineParameter::GetIsDeform()
{
    return this->isDeform;
}

int RefineParameter::GetExclusiveCase()
{
    return this->exclusiveCase;
}

int RefineParameter::GetProjectOrgPoint()
{
    return this->projectOrgPoint;
}

int RefineParameter::GetGeometryUnit()
{
    return this->geometryUnit;
}

string RefineParameter::GetGeometryFileName()
{
    return this->geometryFileName;
}
#endif
#endif

}

