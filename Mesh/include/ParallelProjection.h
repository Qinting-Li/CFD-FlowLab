#pragma once
#ifdef USE_OCC
#include "AMRDef.h"
#include "ElementFactory.h"
#include "OCCInterFace.h"
#include "ZoneGridSolver.h"
#include "Geo_Point.h"
namespace PHSPACE
{
typedef Point3D SimplePoint;
typedef PHVector1D < Point3D > PointArray;
typedef vector< PointArray >  PointVector;
typedef PHVector1D < Point3D > ProjectVector;

class ParallelProjection
{
public:
    ParallelProjection(const string& geometryFileName, int dimension);
    ~ParallelProjection();

private:
    int numberOfZones;

    string globalGeometryFileName;

    //! Geometry unit.\n
    //!   1 -- meter;\n
    //!   2 -- millimeter;\n
    //!   3 -- inch.
    int geometryUnit;

    int dimension;
    
    ZoneGridSolver ** zoneGridSolverContainer;

    OCCInterFace ** occInterFaceList;

    PointVector * zoneProjectVector;

    int refineMethod;

    int symmetryPlane;

    //! ParallelProjection exclusive case.
    //! 0 - NON case;   \n
    //! 1 - JSM-C2-NPOFF case.
    int exclusiveCase;

    //! If project original wall points.
    int projectOrgPoint;

    //! Geometry file name list.
    set<string> geometryFileNameList;

    //! Map from geometry file name to OCC pointer.
    map< string, OCCInterFace * > fileGeometry;

    Region * region;
public:
    void Init();
    void Run();
    void SetSymmetryPlane(const int & symmetryPlane) { this->symmetryPlane = symmetryPlane; }
    void SetZoneGridSolverContainer(int number, ZoneGridSolver ** dataIn);
    void SetProjectVectorArray(PointVector * vectorArrayIn) { this->zoneProjectVector = vectorArrayIn; }
    void SetGeometryUnit(int unitIn) { this->geometryUnit = unitIn; }
    void SetExclusiveCase(int exclusiveCase) { this->exclusiveCase = exclusiveCase; }
    void SetIfProjectOriginalPoint(int projectOrgPoint) { this->projectOrgPoint = projectOrgPoint; }
    PointVector * GetProjectVecotrArray() const { return this->zoneProjectVector; }
    void ComputeGlobalMinWallEdgeLengh();
    void SetRigion(Region * region) { this->region = region; }

private:
    void ComputeWallPointTolerance(ZoneGridSolver * zoneGridSolver);
    void ProjectOnSolidWall();
    void ProjectOnSolidWall(int iZone, ZoneGridSolver * zoneGridSolver);
    void GetGeometryFileNameList();
    void LoadIGES();
    void CheckPointsAttribute(ZoneGridSolver * zoneGridSolver);
    void ProjectChildrenFaces(int iZone, ZoneGridSolver * zoneGridSolver);
    void ProjectAllSurfacePoints(int iZone, ZoneGridSolver * zoneGridSolver);
    bool IsParentFace(const int & parentFaceIndex)
    {
        return parentFaceIndex == PHSPACE::INVALID_INDEX;
    }
    bool IsExclusivePoint(const RDouble & x0, const RDouble & y0, const RDouble & z0);
    void ShowProjectionDistance();
    void ShowProjectionDistance(int iZone, const string & fileName);
    void GetFaceToNodeList(ZoneGridSolver * zoneGridSolver, int boundaryConditionType, PHVectorInt1D & linkMap, PHVectorInt2D & faceToNodeList);
    int  GetNumberofGeometryFile() const { return this->geometryFileNameList.size(); }

    void InitCellProjectType();
    void InitCellProjectType(ZoneGridSolver * zoneGridSolver);

    void FixCellProjectType();
    void FixCellProjectType(ZoneGridSolver * zoneGridSolver);

    void SwapInterfaceProjectType();

    int ResetCellProjectTypeByRightCell();
    int ResetCellProjectTypeByRightCell(ZoneGridSolver * zoneGridSolver);

    bool IfNeedProjection();
private:
    PHVectorRDouble1D pointTolerance;
    RFloat globalMinWallEdgeLength;
    PHVectorInt1D isOriginalWallFace;
    PHVectorInt1D isOriginalWallPoint;
    PHVectorInt1D isInterfacePoint;
    vector<bool>  isOnPatchInterface;
    vector<bool>  isNeedProject;

    RDouble minProject, maxProject;
};

void CompressCellProjectType(ZoneGridSolver * zoneGridSolver, DataContainer *& dataContainer, int neighborZoneIndex);
void DecompressCellProjectType(ZoneGridSolver * zoneGridSolver, DataContainer *& dataContainer, int neighborZoneIndex);

}

#endif