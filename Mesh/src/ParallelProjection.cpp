#ifdef USE_OCC
#include "ParallelProjection.h"
#include "Geo_PointFactory.h"
#include "MixGrid.h"
#include "Glb_Dimension.h"
#include "TK_Log.h"
#include "Region.h"
#include "Zone.h"
#include "cgnslib.h"
#include "PHHeader.h"
#include "PHIO.h"
#include "TK_Time.h"
#include "DataStruct_Sort.h"
#include "TK_Exit.h"
#include "Constants.h"

namespace PHSPACE
{
ParallelProjection::ParallelProjection(const string& geometryFileName, int dimension)
{
    this->globalGeometryFileName = geometryFileName;
    this->dimension = dimension;
    occInterFaceList = 0;
    zoneProjectVector = 0;

    exclusiveCase = 0;
    projectOrgPoint = 0;

    minProject = LARGE;
    maxProject = -LARGE;
}

ParallelProjection::~ParallelProjection()
{
    int nGeometry = this->GetNumberofGeometryFile();
    for(int iGeo = 0; iGeo < nGeometry; ++ iGeo)
    {
        delete occInterFaceList[iGeo];
    }
    delete [] occInterFaceList;
}

void ParallelProjection::SetZoneGridSolverContainer(int number, ZoneGridSolver ** dataIn)
{
    this->zoneGridSolverContainer = dataIn;
    this->numberOfZones = number;
}

void ParallelProjection::Init()
{
    //! refineMethod:
    //!   1 - Mesh deform for each refine(default). 
    //!   2 - Mesh deform after all times of mesh refine.
    refineMethod = 1;

    ComputeGlobalMinWallEdgeLengh();
}

void ParallelProjection::GetGeometryFileNameList()
{
    int nTBC = GlobalBoundaryCondition::GetNumberofBC();
    for (int iPart = 0; iPart < nTBC; ++ iPart)
    {
        if(!GlobalBoundaryCondition::IsSolidWallBC(iPart)) continue;

        SimpleBC * bc = GlobalBoundaryCondition::GetBC(iPart);
        
        string partGeometryFileName;
        bc->GetData("geometryFileName", &partGeometryFileName, PHSTRING, 1);

        if(partGeometryFileName != "")
        {
            //! If the partGeometryFileName is not been defined in boundary_condition.hypara, 
            //! it would be set to be "", meaning 'null'.
            geometryFileNameList.insert(partGeometryFileName);
        }
    }
    
}

void ParallelProjection::LoadIGES()
{
    int nGeometry = this->GetNumberofGeometryFile();
    occInterFaceList = new OCCInterFace * [nGeometry];

    int iGeo = 0;
    set<string>::iterator iter;
    for(iter = geometryFileNameList.begin(); iter != geometryFileNameList.end(); ++ iter)
    {
        OCCInterFace * partoccInterFace = new OCCInterFace();

        const string & geometryFileName = * iter;
        partoccInterFace->SetDimension(this->dimension);
        partoccInterFace->loadIGES(geometryFileName.c_str(), this->geometryUnit);

        PHSPACE::WriteLogFile(geometryFileName + "\n");
        PHSPACE::PrintToWindow(geometryFileName + "\n");
        partoccInterFace->PrintIGESBox();

        occInterFaceList[iGeo ++] = partoccInterFace;

        // Construct map from file name to geometry.
        fileGeometry.insert(pair< string, OCCInterFace * >(geometryFileName, partoccInterFace));
    }
}

void ParallelProjection::Run()
{
    if(IfNeedProjection())
    {
        GetGeometryFileNameList();

        LoadIGES();
    }

    FixCellProjectType();

    if(IfNeedProjection())
    {
        ProjectOnSolidWall();

        ShowProjectionDistance();
    }

    if (PHMPI::GetNumberOfProcessor() == PHMPI::GetNumberofGlobalZones())
    {
        WriteLogFile("    Start waiting reduce ...\n");

        RDouble globalMinDist, globalMaxDist;
        PH_Reduce(&minProject, &globalMinDist, 1, PH_MIN);
        PH_Reduce(&maxProject, &globalMaxDist, 1, PH_MAX);

        ostringstream oss;
        oss << "    Min && Max projection distance is: " << globalMinDist << ", " << globalMaxDist << endl;
        PrintToWindow(oss);
        WriteLogFile(oss);
    }
}

void ParallelProjection::ProjectOnSolidWall()
{
    PrintToWindow("Start Projecting ...", "\n");
    WriteLogFile("Start Projecting ...");

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        ZoneGridSolver * zoneGridSolver = zoneGridSolverContainer[ iZone ];

        //ComputeWallPointTolerance(zoneGridSolver);

        if(refineMethod == PHSPACE::OneRefineOneDeform)
        {
            this->ProjectOnSolidWall(iZone, zoneGridSolver);
        }
        else if(refineMethod == PHSPACE::DeformAfterAllRefine)
        {
            TK_Exit::ExceptionExit("refineMethod == DeformAfterAllRefine");
            //this->ProjectAllSurfacePoints(iZone, zoneGridSolver);
        }
    }

    PrintToWindow("End Projecting ...", "\n");
}

void ParallelProjection::ProjectOnSolidWall(int iZone, ZoneGridSolver * zoneGridSolver)
{
    CheckPointsAttribute(zoneGridSolver);

    ProjectChildrenFaces(iZone, zoneGridSolver);
}

// Check out the faces with SOLID_WALL BC.
void ParallelProjection::CheckPointsAttribute(ZoneGridSolver * zoneGridSolver)
{
    ElementFaceSolver *elementFaceSolver = zoneGridSolver->GetElementFaceSolver();
    Geo_PointFactory  *pointFactory = zoneGridSolver->GetPointFactory();

    int numberOfNodes  = pointFactory->GetNumberOfComputationalNodes();

    int numberOfTotalFaces = elementFaceSolver->GetNumberOfTotalFaces();
    PHVectorInt1D & faceBoundaryTypeContainer = elementFaceSolver->GetFaceBoundaryTypeContainer();
    PHVectorInt1D & parentFaceIndexContainer  = elementFaceSolver->GetParentFaceIndexContainer();
    PHVectorInt2D & faceNodeIndexesContainer  = elementFaceSolver->GetFaceNodeIndexesContainer();
    PHVectorInt1D & isOriginalWallFace        = this->isOriginalWallFace;
    PHVectorInt1D & isOriginalWallPoint       = this->isOriginalWallPoint;
    //PHVectorInt1D & isInterfacePoint        = this->isInterfacePoint;
    
    isOriginalWallFace.resize(numberOfTotalFaces, 0);
    isOriginalWallPoint.resize(numberOfNodes, 0);
    //isInterfacePoint.resize(numberOfNodes, 0);

    int parentFaceNum = 0;
    int childFaceNum  = 0;
    for (int iFace = 0; iFace < numberOfTotalFaces; ++ iFace)
    {
        int bcType = faceBoundaryTypeContainer[ iFace ];
        if(bcType != PHENGLEI::SOLID_SURFACE) continue;

        int parentFaceIndex = parentFaceIndexContainer[ iFace ];

        if (IsParentFace(parentFaceIndex))
        {
            isOriginalWallFace[ iFace ] = 1;

            for (unsigned int iNode = 0; iNode < faceNodeIndexesContainer[ iFace ].size(); ++ iNode)
            {
                int nodeIndex = faceNodeIndexesContainer[ iFace ][ iNode ];
                isOriginalWallPoint[ nodeIndex ] = 1;
            }
        }
    }

    //! Check if the point on the border of two BC: one need project, the other need not.
    vector<bool> & isNeedProject = this->isNeedProject;
    isNeedProject.resize(numberOfNodes);
    PHSPACE::SetField(isNeedProject, true, numberOfNodes);

    ElementPropertySolver * elementPropertySolver = zoneGridSolver->GetElementPropertySolver();
    PHVectorInt1D & projectType = * elementPropertySolver->GetProjectType();
    PHVectorInt1D & leftCell  = elementFaceSolver->GetLeftElementIndexContainer();
    vector<bool> & isSolidWallCell = elementPropertySolver->GetIsSolidWallCell();
    for (int iFace = 0; iFace < numberOfTotalFaces; ++ iFace)
    {
        int bcType = faceBoundaryTypeContainer[ iFace ];
        if(bcType != PHENGLEI::INTERFACE) continue;

        int parentFaceIndex = parentFaceIndexContainer[ iFace ];
        // Do not care of the parent faces.
        if (IsParentFace(parentFaceIndex)) continue;

        int leParent = leftCell[parentFaceIndex];

        if(isSolidWallCell[leParent] && !projectType[leParent])
        {
            for (unsigned iNode = 0; iNode < faceNodeIndexesContainer[ iFace ].size(); ++ iNode)
            {
                int nodeIndex = faceNodeIndexesContainer[ iFace ][ iNode ];
                //isInterfacePoint[ nodeIndex ] = 1;

                isNeedProject[nodeIndex] = false;
            }
        }
    }

    //! Check if the point on the border of different boundary conditions.
    //! The border points will be set to be non-project.
    vector<SimpleBC * > & faceBoundaryCondition = elementFaceSolver->GetFaceBoundaryCondition();
    map< string, OCCInterFace * >::iterator geoIter;
    vector< set<string> > pointBCName;
    pointBCName.resize(numberOfNodes);
    for (int iFace = 0; iFace < numberOfTotalFaces; ++ iFace)
    {
        int bcType = faceBoundaryTypeContainer[ iFace ];
        if(bcType != PHENGLEI::SOLID_SURFACE) continue;

        int parentFaceIndex = parentFaceIndexContainer[ iFace ];

        if (IsParentFace(parentFaceIndex)) continue;

        SimpleBC * boundaryCondition = faceBoundaryCondition[iFace];
        Data_Param * bcParamDataBase = boundaryCondition->GetBCParamDataBase();
        string geometryFileName;
        bcParamDataBase->GetData("geometryFileName", &geometryFileName, PHSTRING, 1);

        //! Children faces.
        for (unsigned iNode = 0; iNode < faceNodeIndexesContainer[ iFace ].size(); ++ iNode)
        {
            int nodeIndex = faceNodeIndexesContainer[ iFace ][ iNode ];
            if(geometryFileName != "")
            {
                pointBCName[nodeIndex].insert(geometryFileName);
            }
            else
            {
                pointBCName[nodeIndex].insert("null");
            }
        }
    }
    for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
    {
        if(pointBCName[iNode].size() > 1)
        {
            //! on the border.
            isNeedProject[iNode] = false;
        }
    }
}

void ParallelProjection::ComputeWallPointTolerance(ZoneGridSolver * zoneGridSolver)
{
    ElementFaceSolver *elementFaceSolver = zoneGridSolver->GetElementFaceSolver();
    Geo_PointFactory  *pointFactory = zoneGridSolver->GetPointFactory();

    PHVectorInt1D & parentFaceIndexContainer  = elementFaceSolver->GetParentFaceIndexContainer();
    PHVectorInt2D & faceNodeIndexesContainer  = elementFaceSolver->GetFaceNodeIndexesContainer();
    PHVectorInt1D & faceBoundaryTypeContainer = elementFaceSolver->GetFaceBoundaryTypeContainer();

    int numberOfNodes      = pointFactory->GetNumberOfComputationalNodes();
    int numberOfTotalFaces = elementFaceSolver->GetNumberOfTotalFaces();

    pointTolerance.clear();
    pointTolerance.resize(numberOfNodes, 1.0e-10);

    PHVectorRDouble1D minDist(numberOfNodes, PHSPACE::LARGE);

    for (int iFace = 0; iFace < numberOfTotalFaces; ++ iFace)
    {
        int bctype = faceBoundaryTypeContainer[iFace];
        if(bctype != PHENGLEI::SOLID_SURFACE) continue;

        int parentFaceIndex = parentFaceIndexContainer[ iFace ];
        if (!IsParentFace(parentFaceIndex))
        {
            int nPoint = faceNodeIndexesContainer[ iFace ].size();
            for (int iNode = 0; iNode < nPoint; ++ iNode)
            {
                int p1 = iNode; 
                int p2 = (iNode + 1) % nPoint;
                int nodeIndex1 = faceNodeIndexesContainer[ iFace ][ p1 ];
                int nodeIndex2 = faceNodeIndexesContainer[ iFace ][ p2 ];

                RDouble dx = pointFactory->GetPoint(nodeIndex1).X() - pointFactory->GetPoint(nodeIndex2).X();
                RDouble dy = pointFactory->GetPoint(nodeIndex1).Y() - pointFactory->GetPoint(nodeIndex2).Y();
                RDouble dz = pointFactory->GetPoint(nodeIndex1).Z() - pointFactory->GetPoint(nodeIndex2).Z();
                RDouble dist = dx * dx + dy * dy + dz * dz;
                
                minDist[nodeIndex1] = MIN(minDist[nodeIndex1], dist);
            }
        }
    }

    for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
    {
        pointTolerance[ iNode ] = sqrt(minDist[ iNode ]);
    }
}

void ParallelProjection::ProjectChildrenFaces(int iZone, ZoneGridSolver * zoneGridSolver)
{
    int readDist = GlobalDataBase::GetIntParaFromDB("readDist");
    ostringstream nameoss;
    nameoss << "./dist/project__" << PHMPI::GetNumberofGlobalZones() << "_" << PHMPI::GetCurrentProcessorID() << ".dist";
    string projectDistFileName = nameoss.str();
    fstream file;
    ios_base::openmode openMode;
    if(readDist)
    {
        openMode = ios_base::in | ios_base::binary;
    }else
    {
        openMode = ios_base::out | ios_base::binary;
    }
    PHSPACE::OpenFile(file, projectDistFileName, openMode);
    DataContainer * cdata = new DataContainer();
    if(readDist)
    {
        ReadFile(file, cdata);
        cdata->moveto_begin();
    }

    ElementFaceSolver *elementFaceSolver = zoneGridSolver->GetElementFaceSolver();
    Geo_PointFactory  *pointFactory = zoneGridSolver->GetPointFactory();

    int numberOfNodes  = pointFactory->GetNumberOfComputationalNodes();

    int numberOfTotalFaces = elementFaceSolver->GetNumberOfTotalFaces();
    vector<SimpleBC * > & faceBoundaryCondition = elementFaceSolver->GetFaceBoundaryCondition();
    PHVectorInt1D & faceBoundaryTypeContainer = elementFaceSolver->GetFaceBoundaryTypeContainer();
    PHVectorInt1D & parentFaceIndexContainer  = elementFaceSolver->GetParentFaceIndexContainer();
    PHVectorInt2D & faceNodeIndexesContainer  = elementFaceSolver->GetFaceNodeIndexesContainer();
    PHVectorInt1D & isOriginalWallFace        = this->isOriginalWallFace;
    PHVectorInt1D & isOriginalWallPoint       = this->isOriginalWallPoint;
    //PHVectorInt1D & isInterfacePoint        = this->isInterfacePoint;
    vector<bool>  & isNeedProject             = this->isNeedProject;

    PointVector * zoneProjectVector = this->GetProjectVecotrArray();
    (*zoneProjectVector)[ iZone ].resize(numberOfNodes);

    PHVectorInt1D isPointProjected;
    isPointProjected.resize(numberOfNodes, 0);

    if(!projectOrgPoint)
    {
        // If the original wall points need not to be projected, then set them to be projected.
        for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
        {
            if(isOriginalWallPoint[ iNode ])
            {
                isPointProjected[ iNode ] = 1;

                (*zoneProjectVector)[iZone][iNode].SetX(0.0);
                (*zoneProjectVector)[iZone][iNode].SetY(0.0);
                (*zoneProjectVector)[iZone][iNode].SetZ(0.0);
            }
        }
    }

    // Compute the number of wall faces.
    int nProjectWallFaces = 0;
    for (int iFace = 0; iFace < numberOfTotalFaces; ++ iFace)
    {
        int bcType = faceBoundaryTypeContainer[ iFace ];
        int parentFaceIndex = parentFaceIndexContainer[ iFace ];
        if(!IsParentFace(parentFaceIndex) && bcType == PHENGLEI::SOLID_SURFACE)
        {
            ++ nProjectWallFaces;
        }
    }

    map< string, OCCInterFace * >::iterator geoIter;
    Data_Param * bcParamDataBase = 0;

    RDouble border[4];

    int count = 0;
    int interval = GetProgressInterval(nProjectWallFaces, 20);
    for (int iFace = 0; iFace < numberOfTotalFaces; ++ iFace)
    {
        int bcType = faceBoundaryTypeContainer[ iFace ];
        // Deal with the solid face ONLY.
        if(bcType != PHENGLEI::SOLID_SURFACE) continue;

        int parentFaceIndex = parentFaceIndexContainer[ iFace ];

        // Do not care of the parent faces.
        if (IsParentFace(parentFaceIndex)) continue;

        ++ count;
        if (count > 0 && count % interval == 0)
        {
            ProgressMonitorToWindows(count, nProjectWallFaces, "Solid wall point projection");
            ProgressMonitorToLogFile(count, nProjectWallFaces, "Solid wall point projection");
        }

        SimpleBC * boundaryCondition = faceBoundaryCondition[parentFaceIndex];
        bcParamDataBase = boundaryCondition->GetBCParamDataBase();
        string geometryFileName;
        bcParamDataBase->GetData("geometryFileName", &geometryFileName, PHSTRING, 1);
        OCCInterFace * occInterFace = 0;
        geoIter = fileGeometry.find(geometryFileName);
        if(geoIter != fileGeometry.end())
        {
            occInterFace = geoIter->second;
        }else
        {
            //! If the geometry file not found, it is not need to be projected.!
            continue;
        }

        RDouble uBorderMin, uBorderMax, vBorderMin, vBorderMax;
        bcParamDataBase->GetData("uBorderMin", &uBorderMin, PHDOUBLE, 1);
        bcParamDataBase->GetData("uBorderMax", &uBorderMax, PHDOUBLE, 1);
        bcParamDataBase->GetData("vBorderMin", &vBorderMin, PHDOUBLE, 1);
        bcParamDataBase->GetData("vBorderMax", &vBorderMax, PHDOUBLE, 1);
        border[0] = uBorderMin; border[1] = uBorderMax; border[2] = vBorderMin; border[3] = vBorderMax;

        // For the children faces.
        for (unsigned iNode = 0; iNode < faceNodeIndexesContainer[ iFace ].size(); ++ iNode)
        {
            int nodeIndex = faceNodeIndexesContainer[ iFace ][ iNode ];

            //! isNeedProject should be deleted, because it's function has been instead of 'borderTolerancePercent'.
            if(!isNeedProject[nodeIndex]) continue;

            if(isPointProjected[ nodeIndex ]) continue;

            isPointProjected[ nodeIndex ] = 1;

            RDouble x0 = pointFactory->GetPoint(nodeIndex).X();
            RDouble y0 = pointFactory->GetPoint(nodeIndex).Y();
            RDouble z0 = pointFactory->GetPoint(nodeIndex).Z();

            if(IsExclusivePoint(x0, y0, z0)) continue;

            gp_Pnt point(x0, y0, z0);

            // for geometry entity.
            gp_Pnt pointProjected = point;
            RDouble projectDist = 1.0e10;
            bool isProjected;
            if (!readDist)
            {
                occInterFace->ComputeNearestPointOnGeometry(point, pointProjected, projectDist, isProjected, border);
                PHWrite(cdata, isProjected);
            }else
            {
                PHRead(cdata, isProjected);
            }

            if(isProjected)
            {
                RDouble xProjected;
                RDouble yProjected;
                RDouble zProjected;
                if(!readDist)
                {
                    xProjected = pointProjected.X();
                    yProjected = pointProjected.Y();
                    zProjected = pointProjected.Z();
                    PHWrite(cdata, xProjected);
                    PHWrite(cdata, yProjected);
                    PHWrite(cdata, zProjected);
                }
                else
                {
                    PHRead(cdata, xProjected);
                    PHRead(cdata, yProjected);
                    PHRead(cdata, zProjected);
                }

                (*zoneProjectVector)[ iZone ][ nodeIndex ].SetX(xProjected - x0);
                (*zoneProjectVector)[ iZone ][ nodeIndex ].SetY(yProjected - y0);
                (*zoneProjectVector)[ iZone ][ nodeIndex ].SetZ(zProjected - z0);

                // reset the coordinate.
                pointFactory->GetPoint(nodeIndex).SetX(xProjected);
                pointFactory->GetPoint(nodeIndex).SetY(yProjected);
                pointFactory->GetPoint(nodeIndex).SetZ(zProjected);

                RDouble distanceTemp = DISTANCE(xProjected - x0, yProjected - y0, zProjected - z0);
                minProject = MIN(minProject, distanceTemp);
                maxProject = MAX(maxProject, distanceTemp);
            }
            else
            {
                (*zoneProjectVector)[ iZone ][ nodeIndex ].SetX(0.0);
                (*zoneProjectVector)[ iZone ][ nodeIndex ].SetY(0.0);
                (*zoneProjectVector)[ iZone ][ nodeIndex ].SetZ(0.0);
            }
        }
        
    }

    // Reset the symmetry plan.
    for (int iFace = 0; iFace < numberOfTotalFaces; ++ iFace)
    {
        int bcType = faceBoundaryTypeContainer[ iFace ];
        if(bcType != PHENGLEI::SYMMETRY) continue;

        int parentFaceIndex = parentFaceIndexContainer[ iFace ];

        // The parent face dose not lie on the geometry.
        if (IsParentFace(parentFaceIndex)) continue;

        // For the children faces.
        for (unsigned iNode = 0; iNode < faceNodeIndexesContainer[ iFace ].size(); ++ iNode)
        {
            int nodeIndex = faceNodeIndexesContainer[ iFace ][ iNode ];

            if(symmetryPlane == 1)
            {
                // reset the X coordinate.
                pointFactory->GetPoint(nodeIndex).SetX(0.0);
                (*zoneProjectVector)[ iZone ][ nodeIndex ].SetX(0);
            }else if(symmetryPlane == 2)
            {
                // reset the Y coordinate.
                pointFactory->GetPoint(nodeIndex).SetY(0.0);
                (*zoneProjectVector)[ iZone ][ nodeIndex ].SetY(0);
            }else if(symmetryPlane == 3)
            {
                // reset the Z coordinate.
                pointFactory->GetPoint(nodeIndex).SetZ(0.0);
                (*zoneProjectVector)[ iZone ][ nodeIndex ].SetZ(0);
            }
        }
    }

    if(!readDist)
    {
        WriteFile(file, cdata);
    }
    
    PHSPACE::CloseFile(file);
    delete cdata;
}

void ParallelProjection::ProjectAllSurfacePoints(int iZone, ZoneGridSolver * zoneGridSolver)
{
    ElementFaceSolver *elementFaceSolver = zoneGridSolver->GetElementFaceSolver();
    Geo_PointFactory  *pointFactory = zoneGridSolver->GetPointFactory();
    vector<SimpleBC * > & faceBoundaryCondition = elementFaceSolver->GetFaceBoundaryCondition();
    PHVectorInt1D & faceBoundaryTypeContainer = elementFaceSolver->GetFaceBoundaryTypeContainer();
    PHVectorInt2D & faceNodeIndexesContainer  = elementFaceSolver->GetFaceNodeIndexesContainer();
    PHVectorInt1D & parentFaceIndexContainer  = elementFaceSolver->GetParentFaceIndexContainer();
    vector<bool>  & isNeedProject             = this->isNeedProject;

    int numberOfNodes = pointFactory->GetNumberOfComputationalNodes();
    int numberOfTotalFaces = elementFaceSolver->GetNumberOfTotalFaces();

    PointVector * zoneProjectVector = this->GetProjectVecotrArray();
    (*zoneProjectVector)[ iZone ].resize(numberOfNodes);

    PHVectorInt1D isPointProjected;
    isPointProjected.resize(numberOfNodes, 0);

    map< string, OCCInterFace * >::iterator geoIter;
    Data_Param * bcParamDataBase = 0;

    RDouble border[4];

    for (int iFace = 0; iFace < numberOfTotalFaces; ++ iFace)
    {
        int bcType = faceBoundaryTypeContainer[ iFace ];
        if(bcType != PHENGLEI::SOLID_SURFACE) continue;

        int parentFaceIndex = parentFaceIndexContainer[ iFace ];

        SimpleBC * boundaryCondition = faceBoundaryCondition[parentFaceIndex];
        bcParamDataBase = boundaryCondition->GetBCParamDataBase();
        string geometryFileName;
        bcParamDataBase->GetData("geometryFileName", &geometryFileName, PHSTRING, 1);
        OCCInterFace * occInterFace = 0;
        geoIter = fileGeometry.find(geometryFileName);
        if(geoIter != fileGeometry.end())
        {
            occInterFace = geoIter->second;
        }else
        {
            //! If the geometry file not found, it is not need to be projected.!
            continue;
        }

        RDouble uBorderMin, uBorderMax, vBorderMin, vBorderMax;
        bcParamDataBase->GetData("uBorderMin", &uBorderMin, PHDOUBLE, 1);
        bcParamDataBase->GetData("uBorderMax", &uBorderMax, PHDOUBLE, 1);
        bcParamDataBase->GetData("vBorderMin", &vBorderMin, PHDOUBLE, 1);
        bcParamDataBase->GetData("vBorderMax", &vBorderMax, PHDOUBLE, 1);
        border[0] = uBorderMin; border[1] = uBorderMax; border[2] = vBorderMin; border[3] = vBorderMax;

        // All NON-Wall faces.
        for (unsigned iNode = 0; iNode < faceNodeIndexesContainer[ iFace ].size(); ++ iNode)
        {
            int nodeIndex = faceNodeIndexesContainer[ iFace ][ iNode ];

            if(!isNeedProject[nodeIndex]) continue;

            if(isPointProjected[ nodeIndex ]) continue;
            isPointProjected[ nodeIndex ] = 1;

            RDouble x0 = pointFactory->GetPoint(nodeIndex).X();
            RDouble y0 = pointFactory->GetPoint(nodeIndex).Y();
            RDouble z0 = pointFactory->GetPoint(nodeIndex).Z();

            if(IsExclusivePoint(x0, y0, z0)) continue;

            gp_Pnt point(x0, y0, z0);

            // for geometry entity.
            gp_Pnt pointProjected = point;
            RDouble minDist = 1.0e10;
            bool isProjected;
            occInterFace->ComputeNearestPointOnGeometry(point, pointProjected, minDist, isProjected, border);

            RDouble xProjected = pointProjected.X();
            RDouble yProjected = pointProjected.Y();
            RDouble zProjected = pointProjected.Z();

            (*zoneProjectVector)[ iZone ][ nodeIndex ].SetX(xProjected - x0);
            (*zoneProjectVector)[ iZone ][ nodeIndex ].SetY(yProjected - y0);
            (*zoneProjectVector)[ iZone ][ nodeIndex ].SetZ(zProjected - z0);

            // reset the coordinate.
            pointFactory->GetPoint(nodeIndex).SetX(xProjected);
            pointFactory->GetPoint(nodeIndex).SetY(yProjected);
            pointFactory->GetPoint(nodeIndex).SetZ(zProjected);
        }

    }

    // Reset the symmetry plan.
    for (int iFace = 0; iFace < numberOfTotalFaces; ++ iFace)
    {
        int bcType = faceBoundaryTypeContainer[ iFace ];
        if(bcType != PHENGLEI::SYMMETRY) continue;

        for (unsigned iNode = 0; iNode < faceNodeIndexesContainer[ iFace ].size(); ++ iNode)
        {
            int nodeIndex = faceNodeIndexesContainer[ iFace ][ iNode ];

            // reset the Z coordinate.
            pointFactory->GetPoint(nodeIndex).SetZ(0.0);
        }
    }
}

void ParallelProjection::ComputeGlobalMinWallEdgeLengh()
{
    globalMinWallEdgeLength = PHSPACE::LARGE;

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        ZoneGridSolver * zoneGridSolver = zoneGridSolverContainer[ iZone ];

        ElementFaceSolver *elementFaceSolver = zoneGridSolver->GetElementFaceSolver();
        Geo_PointFactory  *pointFactory = zoneGridSolver->GetPointFactory();

        PHVectorInt1D & parentFaceIndexContainer  = elementFaceSolver->GetParentFaceIndexContainer();
        PHVectorInt2D & faceNodeIndexesContainer  = elementFaceSolver->GetFaceNodeIndexesContainer();
        PHVectorInt1D & faceBoundaryTypeContainer = elementFaceSolver->GetFaceBoundaryTypeContainer();

        int numberOfNodes      = pointFactory->GetNumberOfComputationalNodes();
        int numberOfTotalFaces = elementFaceSolver->GetNumberOfTotalFaces();

        for (int iFace = 0; iFace < numberOfTotalFaces; ++ iFace)
        {
            int bctype = faceBoundaryTypeContainer[iFace];
            if(bctype != PHENGLEI::SOLID_SURFACE) continue;

            int parentFaceIndex = parentFaceIndexContainer[ iFace ];
            //if (IsParentFace(parentFaceIndex)) continue;

            int nPoint = faceNodeIndexesContainer[ iFace ].size();
            for (int iNode = 0; iNode < nPoint; ++ iNode)
            {
                int p1 = iNode; 
                int p2 = (iNode + 1) % nPoint;
                int nodeIndex1 = faceNodeIndexesContainer[ iFace ][ p1 ];
                int nodeIndex2 = faceNodeIndexesContainer[ iFace ][ p2 ];

                RDouble dx = pointFactory->GetPoint(nodeIndex1).X() - pointFactory->GetPoint(nodeIndex2).X();
                RDouble dy = pointFactory->GetPoint(nodeIndex1).Y() - pointFactory->GetPoint(nodeIndex2).Y();
                RDouble dz = pointFactory->GetPoint(nodeIndex1).Z() - pointFactory->GetPoint(nodeIndex2).Z();
                RDouble dist = dx * dx + dy * dy + dz * dz;

                globalMinWallEdgeLength = MIN(globalMinWallEdgeLength, dist);
            }
        }
    }

    globalMinWallEdgeLength = sqrt(globalMinWallEdgeLength);

    PH_CompareMaxMin(globalMinWallEdgeLength, 2);

    PrintToWindow("  Global min wall edge length:", globalMinWallEdgeLength, "\n");
}

bool ParallelProjection::IsExclusivePoint(const RDouble & x0, const RDouble & y0, const RDouble & z0)
{
    if(exclusiveCase == 0) return false; // It needn't to check.

    const int JSM_SHAPE = 1;
    const int CHNT_SHAPE= 2;

    if(exclusiveCase == JSM_SHAPE)
    {
        // AIAA high-lift-3, JSM, C2 grid.
        // Left
        if(x0 > 1823.5000  && x0 < 1826.0000 &&
            y0 > -309.5000  && y0 < -305.5000 &&
            z0 > -131.0000  && z0 < -124.0000 )
        {
            return true;
        }

        // right
        if(x0 > 3084.0000   && x0 < 3094.5000  &&
            y0 > -2255.0000  && y0 < -2245.0000 &&
            z0 > -49.0000    && z0 < -37.0000  )
        {
            return true;
        }

        // Fuselage.
        if(x0 > 2300   && x0 < 2782  &&
            y0 < -100 &&
            z0 < -290)
        {
            return true;
        }

        if(x0 > 1970   && x0 < 1990  &&
            y0 > -194   && y0 < -172  &&
            z0 < -290)
        {
            return true;
        }


        if(x0 > 2556 && x0 < 2566 && 
            y0 > -863 && y0 < -264 &&
            z0 > -200 && z0 < -103)
        {
            return true;
        }

        if(x0 > 2353 && x0 < 2373 && 
            y0 > -266 && y0 < -264 &&
            z0 > -112 && z0 < -105)
        {
            return true;
        }

        if(x0 > 1689 && x0 < 1760 && 
            y0 > -20 && 
            z0 < -263)
        {
            return true;
        }

        if(x0 > 1957  && x0 <  1990 &&
            y0 > -285  && y0 < -260 &&
            z0 > -230  && z0 < -205  )
        {
            return true;
        }

        if(x0 > 1959  && x0 < 1980  &&
            y0 > -270  && y0 < -264 &&
            z0 > -90  && z0 < -80  )
        {
            return true;
        }

        if(x0 > 3025 && x0 < 3080  &&
            y0 > -10  &&
            z0 < -260  )
        {
            return true;
        }

        if(x0 > 3135 && x0 < 3145 && 
            y0 > -2220 && y0 < -2213 &&
            z0 > -35 && z0 < -34)
        {
            return true;
        }

        // most-inner-small
        if(x0 > 2690 && x0 < 2760 && 
            y0 > -321 && y0 < -300 &&
            z0 > -241 && z0 < -203)
        {
            return true;
        }

        // outer
        if(x0 > 2690 && x0 < 3102 && 
            y0 > -1555 && y0 < -1500 &&
            z0 > -180 && z0 < -90)
        {
            return true;
        }

        // middle
        if(x0 > 2460 && x0 < 2914 && 
            y0 > -1088 && y0 < -1035 &&
            z0 > -218 && z0 < -120)
        {
            return true;
        }

        // inner
        if(x0 > 2310 && x0 < 2828 && 
            y0 > -680 && y0 < -620 &&
            z0 > -300 && z0 < -145)
        {
            return true;
        }


    }
    else if(exclusiveCase == CHNT_SHAPE)
    {
        // vertical
        if(x0 > 21691 && y0 > 6815)
        {
            return true;
        }

        // horizontal
        if(x0 > 19263 && fabs(z0) > 4729)
        {
            return true;
        }
    }

    return false;
}

void ParallelProjection::ShowProjectionDistance()
{
    if(PHMPI::GetNumberOfProcessor() > 512)
    {
        PrintToWindow("Do not dump ParallelProjection distance ! \n");
        return;
    }
    PrintToWindow("ParallelProjection distance dumping ... \n");
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        ZoneGridSolver * zoneGridSolver = zoneGridSolverContainer[ iZone ];
        int localID  = zoneGridSolver->GetZoneIndex();
        int globalID = PHMPI::GetLocalZoneIDToGlobalZoneID(localID);

        string fileName = AddSymbolToFileName("./grid/projectDistance.dat", '_', globalID);
        fileName = ChangeExtensionOfFileName(fileName, "plt");

        ShowProjectionDistance(iZone, fileName);
    }
}

void ParallelProjection::ShowProjectionDistance(int iZone, const string & fileName)
{
    ZoneGridSolver * zoneGridSolver = zoneGridSolverContainer[ iZone ];

    ElementFaceSolver *elementFaceSolver = zoneGridSolver->GetElementFaceSolver();
    Geo_PointFactory  *pointFactory = zoneGridSolver->GetPointFactory();

    PHVectorInt1D & faceBoundaryTypeContainer = elementFaceSolver->GetFaceBoundaryTypeContainer();
    PHVectorInt2D & faceNodeIndexContainer = elementFaceSolver->GetFaceNodeIndexesContainer();

    PointVector * zoneProjectVector = this->GetProjectVecotrArray();

    int numberOfNodes         = pointFactory->GetNumberOfComputationalNodes();
    int numberOfCells         = zoneGridSolver->GetNumberOfElements();
    int numberOfBoundaryFaces = zoneGridSolver->GetNumberOfBoundaryFaces();
    int numberOfFaces         = zoneGridSolver->GetNumberOfTotalFaces();

    PHVectorInt1D & leftElementIndexContainer = elementFaceSolver->GetLeftElementIndexContainer();
    PHVectorInt1D & rightElementIndexContainer = elementFaceSolver->GetRightElementIndexContainer();

    PHVectorString1D titleOfTecplot;
    titleOfTecplot.push_back("title=\"Unstructured Grid\"");
    titleOfTecplot.push_back("variables=");
    titleOfTecplot.push_back("\"x\"");
    titleOfTecplot.push_back("\"y\"");
    titleOfTecplot.push_back("\"z\"");
    titleOfTecplot.push_back("\"projDist\"");

    int wordWidth = 20;
    //设置点之间的容许公差，小于此公差的两个点被认定为一个。
    RDouble diff = 1.0e-10;

    int iCount = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        iCount += faceNodeIndexContainer[ iFace ].size();
    }

    int numberOfWordsInEachLine = 5;

    using namespace PHENGLEI;
    PHVectorInt2D faceToNodeList(numberOfBoundaryFaces);
    PHVectorInt1D linkMap;
    set< int > bcList;

    bool doseNeedDump = false;
    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        int boundaryConditionType = faceBoundaryTypeContainer[ iFace ];
        if(boundaryConditionType != PHENGLEI::SOLID_SURFACE) continue;
        bcList.insert(boundaryConditionType);
        if(boundaryConditionType == PHENGLEI::SOLID_SURFACE)
        {
            doseNeedDump = true;
        }
    }


    if(!doseNeedDump)
    {
        return;
    }

    fstream file;
    PHSPACE::OpenFile(file, fileName, ios_base::out);

    for (size_t i = 0; i < titleOfTecplot.size(); ++ i)
    {
        file << titleOfTecplot[ i ] << "\n";
    }

    int localID  = zoneGridSolver->GetZoneIndex();
    int globalID = PHMPI::GetLocalZoneIDToGlobalZoneID(localID);

    for (set< int >::iterator iter = bcList.begin(); iter != bcList.end(); ++ iter)
    {
        if (* iter < 0 ) continue;
        if (* iter != PHENGLEI::SOLID_SURFACE && * iter != PHENGLEI::SYMMETRY) continue;

        faceToNodeList.resize(0);
        linkMap.resize(0);

        GetFaceToNodeList(zoneGridSolver, * iter, linkMap, faceToNodeList);
        ostringstream osZone;
        osZone << "\"" << globalID << "BC=" << * iter << "\" ";
        string zoneTitle = osZone.str();

        file << "zone T = " << zoneTitle << " N = " << linkMap.size() << "  E = " << faceToNodeList.size()  << " f = FEPOINT, ET = QUADRILATERAL\n";
        int wordWidth = 20;
        for (std::size_t i = 0; i < linkMap.size(); ++ i)
        {
            int it = linkMap[ i ];

            RDouble dx = (*zoneProjectVector)[ iZone ][ it ].X();
            RDouble dy = (*zoneProjectVector)[ iZone ][ it ].Y();
            RDouble dz = (*zoneProjectVector)[ iZone ][ it ].Z();
            RDouble dist = DISTANCE(dx, dy, dz);

            file << setiosflags(ios::left);
            file << setiosflags(ios::scientific);
            file << setprecision(10);
            file << setw(wordWidth) << pointFactory->GetPoint(it).X()
            << setw(wordWidth) << pointFactory->GetPoint(it).Y()
            << setw(wordWidth) << pointFactory->GetPoint(it).Z()
            << setw(wordWidth) << dist;
            file << "\n";
        }

        for (std::size_t i = 0; i < faceToNodeList.size(); ++ i)
        {
            //本单元面的顶点个数
            int np     = faceToNodeList[ i ].size();
            int index0 = faceToNodeList[ i ][ 0 ];
            for (int ip = 0; ip < np; ++ ip)
            {
                file << faceToNodeList[ i ][ ip ] + 1;
                if (ip != np - 1) file << " ";
            }

            for (int ip = np; ip < 4; ++ ip)
            {
                file << " " << index0 + 1;
            }
            file << "\n";
        }
    }

    PHSPACE::CloseFile(file);
}

void ParallelProjection::GetFaceToNodeList(ZoneGridSolver * zoneGridSolver, int boundaryConditionType, PHVectorInt1D & linkMap, PHVectorInt2D & faceToNodeList)
{
    ElementFaceSolver *elementFaceSolver = zoneGridSolver->GetElementFaceSolver();
    Geo_PointFactory  *pointFactory = zoneGridSolver->GetPointFactory();

    PHVectorInt1D & faceBoundaryTypeContainer = elementFaceSolver->GetFaceBoundaryTypeContainer();
    PHVectorInt2D & faceNodeIndexContainer = elementFaceSolver->GetFaceNodeIndexesContainer();

    PHVectorInt1D & leftElementIndexContainer = elementFaceSolver->GetLeftElementIndexContainer();
    PHVectorInt1D & rightElementIndexContainer = elementFaceSolver->GetRightElementIndexContainer();

    int numberOfNodes         = pointFactory->GetNumberOfComputationalNodes();
    int numberOfCells         = zoneGridSolver->GetNumberOfElements();
    int numberOfBoundaryFaces = zoneGridSolver->GetNumberOfBoundaryFaces();
    int numberOfFaces         = zoneGridSolver->GetNumberOfTotalFaces();

    //建立二叉树点集作为查询之用
    PHVectorInt1D locFaceToNode(4);
    set< DataStruct_Sort< int > > setForSorting;
    DataStruct_Sort< int > dataForSorting;

    faceToNodeList.resize(0);
    linkMap.resize(0);

    int nodePosition = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        int bcType = faceBoundaryTypeContainer[iFace];
        if(bcType != PHENGLEI::SOLID_SURFACE) continue;

        int faceNodeNumber = faceNodeIndexContainer[ iFace ].size();
        if (boundaryConditionType != bcType)
        {
            nodePosition += faceNodeNumber;
            continue;
        }
        locFaceToNode.resize(0);
        for (int iNodeInFace = 0; iNodeInFace < faceNodeNumber; ++ iNodeInFace)
        {
            int index = faceNodeIndexContainer[iFace][iNodeInFace];

            dataForSorting.value = index;
            set< DataStruct_Sort< int > >::iterator iter = setForSorting.find(dataForSorting);

            if (iter == setForSorting.end())
            {
                linkMap.push_back(index);
                dataForSorting.index = linkMap.size() - 1;

                setForSorting.insert(dataForSorting);
                locFaceToNode.push_back(dataForSorting.index);
            }
            else
            {
                locFaceToNode.push_back(iter->index);
            }

        }
        nodePosition += faceNodeNumber;
        faceToNodeList.push_back(locFaceToNode);
    }
}

void ParallelProjection::InitCellProjectType()
{
    int numberOfZones = PHMPI::GetNumberofLocalZones();
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        ZoneGridSolver * zoneGridSolver = zoneGridSolverContainer[ iZone ];

        InitCellProjectType(zoneGridSolver);
    }
}

void ParallelProjection::InitCellProjectType(ZoneGridSolver * zoneGridSolver)
{
    ElementPropertySolver * elementPropertySolver = zoneGridSolver->GetElementPropertySolver();
    ElementFaceSolver * elementFaceSolver = zoneGridSolver->GetElementFaceSolver();

    int nTCell = zoneGridSolver->GetNumberOfElements();
    int nBoundFace = elementFaceSolver->GetNumberOfBoundaryFaces();
    int nTotal = nTCell + nBoundFace;

    int numberOfTotalFaces = elementFaceSolver->GetNumberOfTotalFaces();
    vector<SimpleBC * > & faceBoundaryCondition = elementFaceSolver->GetFaceBoundaryCondition();
    PHVectorInt1D & faceBoundaryTypeContainer = elementFaceSolver->GetFaceBoundaryTypeContainer();
    PHVectorInt1D & parentFaceIndexContainer  = elementFaceSolver->GetParentFaceIndexContainer();
    PHVectorInt1D & leftCell  = elementFaceSolver->GetLeftElementIndexContainer();
    PHVectorInt1D & rightCell = elementFaceSolver->GetRightElementIndexContainer();

    PHVectorInt1D & elementType = elementPropertySolver->GetElementTypeContainer();
    PHVectorInt1D & projectType = * elementPropertySolver->GetProjectType();
    PHVectorInt1D & wallBCID    = * elementPropertySolver->GetWallBCID();
    projectType.resize(nTotal);
    wallBCID.resize(nTotal);
    SetField(projectType, 0, nTotal); // Set all cells to be non projected.
    SetField(wallBCID, -1, nTotal);

    map< string, OCCInterFace * >::iterator geoIter;
    set<string>::iterator geoFileNameIter;
    Data_Param * bcParamDataBase = 0;

    vector<bool> & isSolidWallCell = elementPropertySolver->GetIsSolidWallCell();
    isSolidWallCell.resize(nTCell);
    PHSPACE::SetField(isSolidWallCell, false, nTCell);
    for (int iFace = 0; iFace < numberOfTotalFaces; ++ iFace)
    {
        int bcType = faceBoundaryTypeContainer[ iFace ];
        if(bcType != PHENGLEI::SOLID_SURFACE) continue;

        int le = leftCell[iFace];
        isSolidWallCell[le] = true;
    }

    for (int iFace = 0; iFace < numberOfTotalFaces; ++ iFace)
    {
        int bcType = faceBoundaryTypeContainer[ iFace ];
        if(bcType != PHENGLEI::SOLID_SURFACE) continue;

        int parentFaceIndex = parentFaceIndexContainer[ iFace ];

        // Do not care of the children faces.
        if (!IsParentFace(parentFaceIndex)) continue;

        //! The boundary condition is same to it's parent face.
        //BoundaryCondition * boundaryCondition = faceBoundaryCondition[parentFaceIndex];
        SimpleBC * boundaryCondition = faceBoundaryCondition[iFace];

        if(boundaryCondition == 0)
        {
            TK_Exit::ExceptionExit("boundaryCondition == 0");
        }

        bcParamDataBase = boundaryCondition->GetBCParamDataBase();
        string geometryFileName;
        bcParamDataBase->GetData("geometryFileName", &geometryFileName, PHSTRING, 1);

        int bcID = 0;
        for(geoFileNameIter = geometryFileNameList.begin(); geoFileNameIter != geometryFileNameList.end(); ++ geoFileNameIter)
        {
            const string & name = * geoFileNameIter;
            if(geometryFileName == name)
            {
                int le = leftCell[iFace];
                if(le >= nTotal || le < 0)
                {
                    cout << le << endl;
                    TK_Exit::ExceptionExit("le >= nTotal || le < 0");
                }

                wallBCID[le]    = bcID;
                projectType[le] = 1;
                break;
            }

            ++ bcID;
        }
    }
}

void ParallelProjection::FixCellProjectType()
{
    WriteLogFile("Start fix cell project type!");

    if(GetDim() == PHSPACE::TWO_D) return;
    TimeSpan fixTime;

    InitCellProjectType();

//  int anisoRefine = GlobalDataBase::GetIntParaFromDB("anisoRefine");
//  if(anisoRefine == 0) return;

    int infectTime = 0;
    int isContinueInfect = 1;
    while (isContinueInfect)
    {
        isContinueInfect = 0;

        if(infectTime > 0)
        {
            PrintToWindow("\nInfect cell project type time: ", infectTime, "......");
            WriteLogFile("\nInfect cell project type time: ", infectTime, "......");
        }

        SwapInterfaceProjectType();

        isContinueInfect = ResetCellProjectTypeByRightCell();

        PH_CompareMaxMin(isContinueInfect, 1);

        ++ infectTime;
    }
    PrintToWindow("\n");
    

    fixTime.ShowTimeSpanToLogfile("Fix cell project type");
    WriteLogFile("End fix cell project type!");
}

void ParallelProjection::SwapInterfaceProjectType()
{
    CompressData compressData     = CompressCellProjectType;
    DeCompressData decompressData = DecompressCellProjectType;

    SwapInterfaceData(this->region, this->zoneGridSolverContainer, compressData, decompressData);
}

int ParallelProjection::ResetCellProjectTypeByRightCell()
{
    int isContinueInfect = 0;

    int numberOfZones = PHMPI::GetNumberofLocalZones();
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        ZoneGridSolver * zoneGridSolver = zoneGridSolverContainer[ iZone ];

        int temp = ResetCellProjectTypeByRightCell(zoneGridSolver);
        isContinueInfect = MAX(isContinueInfect, temp);
    }

    return isContinueInfect;
}

int ParallelProjection::ResetCellProjectTypeByRightCell(ZoneGridSolver * zoneGridSolver)
{
    ElementPropertySolver * elementPropertySolver = zoneGridSolver->GetElementPropertySolver();
    ElementFaceSolver * elementFaceSolver = zoneGridSolver->GetElementFaceSolver();

    int nTCell = zoneGridSolver->GetNumberOfElements();
    int nBoundFace = elementFaceSolver->GetNumberOfBoundaryFaces();
    int nTotalFace = elementFaceSolver->GetNumberOfTotalFaces();

    PHVectorInt1D & elementType = elementPropertySolver->GetElementTypeContainer();
    PHVectorInt1D & projectType = * elementPropertySolver->GetProjectType();
    PHVectorInt1D & wallBCID    = * elementPropertySolver->GetWallBCID();

    PHVectorInt1D & faceType = elementFaceSolver->GetFaceTypeContainer();
    PHVectorInt1D & faceBoundaryType = elementFaceSolver->GetFaceBoundaryTypeContainer();
    PHVectorInt1D & leftCell  = elementFaceSolver->GetLeftElementIndexContainer();
    PHVectorInt1D & rightCell = elementFaceSolver->GetRightElementIndexContainer();
    PHVectorInt1D & parentFaceIndexContainer  = elementFaceSolver->GetParentFaceIndexContainer();

    vector<bool> & isSolidWallCell = elementPropertySolver->GetIsSolidWallCell();

    int isContinueInfect = 0;
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        if(faceBoundaryType[iFace] != PHENGLEI::INTERFACE) continue;
        //if(faceType[iFace] != QUAD_4) continue;

        int parentFaceIndex = parentFaceIndexContainer[ iFace ];
        // Do not care of the children faces.
        if (!IsParentFace(parentFaceIndex)) continue;

        int le = leftCell[iFace];
        int re = rightCell[iFace];

        if(le >= nTCell)
        {
            SWAP(le, re);
        }
        if(re < 0)
        {
            re = -re;
        }

        if(isSolidWallCell[le])
        {
            if(projectType[re] < projectType[le])
            {
                //! re needn't project, so set le needn't project.
                projectType[le] = projectType[re];
                isContinueInfect = 1;
            }

            if(projectType[re] == 1 && projectType[le] == 1)
            {
                //! both le && re need project, judge the bc name then.
                if(wallBCID[re] != wallBCID[le])
                {
                    //! BC name of le is not equal to re, so set le needn't project.
                    projectType[le] = 0;
                    isContinueInfect = 1;
                }
            }
        }
    }

    return isContinueInfect;
}

bool ParallelProjection::IfNeedProjection()
{
    if(PHMPI::GetCurrentProcessorID() == PHMPI::GetServerProcessorID())
    {
        return true;
    }

    bool ifNeed = false;

    int nZones = PHMPI::GetNumberofLocalZones();
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        ZoneGridSolver * zoneGridSolver = zoneGridSolverContainer[ iZone ];
        ElementFaceSolver * elementFaceSolver = zoneGridSolver->GetElementFaceSolver();

        int numberOfTotalFaces = elementFaceSolver->GetNumberOfTotalFaces();
        PHVectorInt1D & faceBoundaryTypeContainer = elementFaceSolver->GetFaceBoundaryTypeContainer();
        for (int iFace = 0; iFace < numberOfTotalFaces; ++ iFace)
        {
            int bcType = faceBoundaryTypeContainer[ iFace ];

            if(bcType == PHENGLEI::SOLID_SURFACE)
            {
                ifNeed = true;
                break;
            }
        }
    }
    return ifNeed;
}

void CompressCellProjectType(ZoneGridSolver * zoneGridSolver, DataContainer *& dataContainer, int neighborZoneIndex)
{
    zoneGridSolver->CompressCellProjectType(dataContainer, neighborZoneIndex);
}

void DecompressCellProjectType(ZoneGridSolver * zoneGridSolver, DataContainer *& dataContainer, int neighborZoneIndex)
{
    zoneGridSolver->DecompressCellProjectType(dataContainer, neighborZoneIndex);
}


}

#endif