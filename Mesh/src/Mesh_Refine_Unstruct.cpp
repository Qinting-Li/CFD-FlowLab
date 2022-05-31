#include "Mesh_Refine_Unstruct.h"
#include "GlobalDataBase.h"
#include "TK_Log.h"
#include "TK_Time.h"
#include "Glb_Dimension.h"
#include "Geo_Element.h"
#include "GridType.h"
#include "PHHeader.h"
#include "ComputationalGrid.h"

using namespace std;
namespace PHSPACE
{
Mesh_Refine_Unstruct::Mesh_Refine_Unstruct(const string &gridFileNameIn) :
    Mesh_Refine(gridFileNameIn)
{
    gridRefineControllers = 0;
}

Mesh_Refine_Unstruct::~Mesh_Refine_Unstruct()
{
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        delete gridRefineControllers[iZone];
    }
    delete []gridRefineControllers;

}

void Mesh_Refine_Unstruct::GenerateAndDumpComputationalGrid()
{
    WriteLogFile("Start Generating Grid!");

    AllocateComputationalGrid();

    GenerateComputationalGrid();

    GenerateGlobalInterfaceLink();

    DumpComputationalGrid();
}

void Mesh_Refine_Unstruct::GenerateGlobalInterfaceLink()
{
    ComputationalGridProcess *computationalGridProcess = new ComputationalGridProcess(region);
    computationalGridProcess->Initialize(refinedGrids, numberOfZones);
    computationalGridProcess->PostProcessMultiZoneComputationalGrids();

    delete computationalGridProcess;
}

void Mesh_Refine_Unstruct::GenerateComputationalGrid()
{
    int nTCellLocal = 0;
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        gridRefineControllers[iZone]->GenerateComputationalGrid(refinedGrids[iZone]);

        nTCellLocal += refinedGrids[iZone]->GetNTotalCell();
    }

    RDouble nTCellGlobal = nTCellLocal * 1.0;
    RDouble ntCellLocalTemp = nTCellLocal * 1.0;
    PH_Reduce(&ntCellLocalTemp, &nTCellGlobal, 1, PH_SUM);

    int maxNTCellLocal, minNTCellLocal;
    PH_Reduce(&nTCellLocal, &maxNTCellLocal, 1, PH_MAX);
    PH_Reduce(&nTCellLocal, &minNTCellLocal, 1, PH_MIN);
    PrintToWindow("  Number of Global cell: ", nTCellGlobal, "\n");
    PrintToWindow("  Min & Max number of cell after adapt: ", minNTCellLocal, maxNTCellLocal, "\n");
}

void Mesh_Refine_Unstruct::AllocateComputationalGrid()
{
    refinedGrids = new Grid *[numberOfZones];
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        Grid *grid = CreateGridGeneral(UNSTRUCTGRID, new GridID(iZone), 0, GetDim());
        refinedGrids[iZone] = grid;
    }
}

void Mesh_Refine_Unstruct::AllocateMemory()
{
    numberOfZones = PHMPI::GetNumberofLocalZones();
    gridRefineControllers = new Mesh_RefineController *[numberOfZones];
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        Mesh_RefineController *gridRefineController = new Mesh_RefineController();
        gridRefineController->SetZoneIndex(iZone);
        gridRefineControllers[iZone] = gridRefineController;
    }
}

void Mesh_Refine_Unstruct::ConstructGridTopo()
{
    WriteLogFile("Start Constructing Grid Info!");

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int globalZoneIndex = PHMPI::GetLocalZoneIDToGlobalZoneID(iZone);

        Grid *grid = GetGrid(globalZoneIndex, 0);
        grid->ComputeMetrics();
        grid->SkewnessSummary();

        ConstructGridTopo(grid);
    }

    WriteLogFile("End Constructing Grid Info!");
}

void Mesh_Refine_Unstruct::ConstructGridTopo(Grid *grid)
{
    ConstructCellTopo(grid);

    ConstructNodeTopo(grid);

    ConstructFaceTopo(grid);

    ConstructMetrics(grid);

    ConstructInterFaceTopo(grid);
}

void Mesh_Refine_Unstruct::ConstructMetrics(Grid *grid)
{
    UnstructGrid *unstructuredGrid = UnstructGridCast(grid);

    RDouble *faceAreaContainer    = unstructuredGrid->GetFaceArea();
    RDouble *faceNormalXContainer = unstructuredGrid->GetFaceNormalX();
    RDouble *faceNormalYContainer = unstructuredGrid->GetFaceNormalY();
    RDouble *faceNormalZContainer = unstructuredGrid->GetFaceNormalZ();

    int localIndex = grid->GetZoneLocalID();
    Mesh_RefineController *currentController = gridRefineControllers[localIndex];
    Mesh_FaceTopo *faceTopo = currentController->GetFaceFactory();

    PHVectorRDouble1D &faceArea    = faceTopo->GetFaceArea();
    PHVectorRDouble1D &faceNormalX = faceTopo->GetFaceNormalX();
    PHVectorRDouble1D &faceNormalY = faceTopo->GetFaceNormalY();
    PHVectorRDouble1D &faceNormalZ = faceTopo->GetFaceNormalZ();

    int nTotalFace = unstructuredGrid->GetNTotalFace();
    faceArea.resize   (nTotalFace);
    faceNormalX.resize(nTotalFace);
    faceNormalY.resize(nTotalFace);
    faceNormalZ.resize(nTotalFace);

    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        faceArea   [iFace] = faceAreaContainer   [iFace];
        faceNormalX[iFace] = faceNormalXContainer[iFace];
        faceNormalY[iFace] = faceNormalYContainer[iFace];
        faceNormalZ[iFace] = faceNormalZContainer[iFace];
    }
}

void Mesh_Refine_Unstruct::ConstructCellTopo(Grid *grid)
{
    int dimension = GetDim();
    UnstructGrid *unstructuredGrid = UnstructGridCast(grid);

    int nTotalCell = unstructuredGrid->GetNTotalCell();
    int *nodeNumberOfEachCell = unstructuredGrid->GetNodeNumberOfEachCell();
    int *cell2Node  = unstructuredGrid->GetCell2Node();

    int localIndex = grid->GetZoneLocalID();
    Mesh_RefineController *currentController = gridRefineControllers[localIndex];
    Mesh_CellTopo *cellTopo = currentController->GetCellFactory();
    PHVectorInt1D &computationalToGlobalCellIndex = cellTopo->GetComputationalToGlobalCellIndex();
    PHVectorInt1D &globalToComputationalCellIndex = cellTopo->GetGlobalToComputationalCellIndex();

    computationalToGlobalCellIndex.resize(nTotalCell);
    globalToComputationalCellIndex.resize(nTotalCell);

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        computationalToGlobalCellIndex[iCell] = iCell;
        globalToComputationalCellIndex[iCell] = iCell;
    }

    PHVectorInt1D &cellType = cellTopo->GetCellType();
    cellType.resize(nTotalCell);

    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        int nodeNumber = nodeNumberOfEachCell[ iCell ];
        int iCellType = GetCellType(dimension, nodeNumber);
        cellType[iCell] = iCellType;
    }

    cellTopo->SetAllCellLevelToZero();
    cellTopo->SetAllParentIndexToRootState();

    PHVectorInt2D &cell2NodeArray = cellTopo->GetCell2Node();
    cell2NodeArray.resize(nTotalCell);
    int count = 0;
    for (int iCell = 0; iCell < nTotalCell; ++ iCell)
    {
        int nodeNumber = nodeNumberOfEachCell[iCell];
        cell2NodeArray[iCell].resize(nodeNumber);
        for (int iNode = 0; iNode < nodeNumber; ++ iNode)
        {
            cell2NodeArray[iCell][iNode] = cell2Node[count ++];
        }
    }

    int anisoRefine = this->refineParameter->GetAnisoRefineType();
    if(anisoRefine)
    {
        int *faceNumberOfEachCell = unstructuredGrid->GetFaceNumberOfEachCell();
        int **cell2Face = unstructuredGrid->GetCell2Face();

        PHVectorInt2D &cell2FaceArray = cellTopo->GetCell2Face();
        cell2FaceArray.resize(nTotalCell);

        for (int iCell = 0; iCell < nTotalCell; ++ iCell)
        {
            int faceNumber = faceNumberOfEachCell[iCell];
            cell2FaceArray[iCell].resize(faceNumber);
            for (int iFace = 0; iFace < faceNumber; ++ iFace)
            {
                cell2FaceArray[iCell][iFace] = cell2Face[iCell][iFace];
            }
        }
    }
}

void Mesh_Refine_Unstruct::ConstructNodeTopo(Grid *grid)
{
    UnstructGrid *unstructuredGrid = UnstructGridCast(grid);

    RDouble *x = unstructuredGrid->GetX();
    RDouble *y = unstructuredGrid->GetY();
    RDouble *z = unstructuredGrid->GetZ();

    int localIndex = grid->GetZoneLocalID();
    Mesh_RefineController *currentController = gridRefineControllers[localIndex];
    Geo_PointFactory *nodeTopo = currentController->GetNodeFactory();

    RDouble minDistance, maxDistance;
    unstructuredGrid->GetMinMaxDS(minDistance, maxDistance);
    nodeTopo->SetMinDistance(minDistance);

    int nTotalNode = unstructuredGrid->GetNTotalNode();
    for (int iPoint = 0; iPoint < nTotalNode; ++ iPoint)
    {
        Point3D *point = new Point3D(x[iPoint], y[iPoint], z[iPoint], iPoint);

        bool isNodeExist;
        int pointIndex = nodeTopo->AddPoint(point, isNodeExist);
    }

    int oldNumberOfComputationalNodes = nTotalNode;
    nodeTopo->SetOldNumberOfComputationalNodes(oldNumberOfComputationalNodes);
}

void Mesh_Refine_Unstruct::ConstructFaceTopo(Grid *grid)
{
    int dimension = GetDim();
    UnstructGrid *unstructuredGrid = UnstructGridCast(grid);

    int nTotalFace            = unstructuredGrid->GetNTotalFace();
    int nBoundFace            = unstructuredGrid->GetNBoundFace();
    int nTotalCell            = unstructuredGrid->GetNTotalCell();
    int *face2Node            = unstructuredGrid->GetFace2Node();
    int *leftCellOfFace       = unstructuredGrid->GetLeftCellOfFace();
    int *rightCellOfFace      = unstructuredGrid->GetRightCellOfFace();
    int *nodeNumberOfEachFace = unstructuredGrid->GetNodeNumberOfEachFace();

    int localIndex = grid->GetZoneLocalID();
    Mesh_RefineController *currentController = gridRefineControllers[localIndex];
    Mesh_FaceTopo *faceTopo = currentController->GetFaceFactory();

    PHVectorInt1D &computationalToGlobalFaceIndex = faceTopo->GetComputationalToGlobalFaceIndex();
    computationalToGlobalFaceIndex.resize(nTotalFace);

    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        computationalToGlobalFaceIndex[iFace] = iFace;
    }

    PHVectorInt1D &leftCellIndex  = faceTopo->GetLeftCellIndex();
    PHVectorInt1D &rightCellIndex = faceTopo->GetRightCellIndex();
    leftCellIndex.resize(nTotalFace);
    rightCellIndex.resize(nTotalFace);

    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        int leftCell  = leftCellOfFace [iFace];
        int rightCell = rightCellOfFace[iFace];
        if(rightCell >= nTotalCell)
        {
            rightCell = - rightCell;
        }

        leftCellIndex [iFace] = leftCell;
        rightCellIndex[iFace] = rightCell;
    }

    int nodeCount = 0;
    for (int iFace = 0; iFace < nTotalFace; ++iFace)
    {
        int nodeNumber = nodeNumberOfEachFace[iFace];
        PHVectorInt1D iFace2Node;
        iFace2Node.resize(nodeNumber);
        for (int iNode = 0; iNode < nodeNumber; ++ iNode)
        {
            iFace2Node[iNode] = face2Node[nodeCount ++];
        }

        sort(iFace2Node.begin(), iFace2Node.end());

        MultiIndex faceOfMultiIndex(nodeNumber, iFace);
        faceOfMultiIndex.SetData(&iFace2Node[0]);
        int id = faceTopo->FindFace(faceOfMultiIndex);
        if(id == PHSPACE::INVALID_INDEX)
        {
            faceTopo->InsertFace(faceOfMultiIndex);
        }
    }

    PHVectorInt2D &face2NodeArray = faceTopo->GetFace2Node();
    face2NodeArray.resize(nTotalFace);
    nodeCount = 0;
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        int nodeNumber = nodeNumberOfEachFace[iFace];
        face2NodeArray[iFace].resize(nodeNumber);

        for (int iNode = 0; iNode < nodeNumber; ++ iNode)
        {
            face2NodeArray[iFace][iNode] = face2Node[nodeCount ++];
        }
    }

    vector<SimpleBC * > &faceBoundaryCondition = faceTopo->GetFaceBoundaryCondition();
    faceBoundaryCondition.resize(nTotalFace);
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        faceBoundaryCondition[iFace] = 0;
    }
    UnstructBCSet *unstructBCSet = unstructuredGrid->GetUnstructBCSet();
    int *bcRegionIDofBCFace = unstructBCSet->GetBCFaceInfo();

    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        UnstructBC *bcRegion = unstructBCSet ->GetBCRegion(bcRegionIDofBCFace[iFace]);

        SimpleBC *iFaceBoundaryCondtion = new SimpleBC();
        int boundaryType = bcRegion->GetBCType();
        iFaceBoundaryCondtion->SetBCType(boundaryType);
        if(boundaryType != PHENGLEI::INTERFACE)
        {
            string bcName = bcRegion->GetBCName();
            iFaceBoundaryCondtion->SetBCName(bcName);
        }
        else
        {
            string bcName = "Interface";
            iFaceBoundaryCondtion->SetBCName(bcName);
        }
        faceBoundaryCondition[iFace] = iFaceBoundaryCondtion;
    }

    PHVectorInt1D &faceType = faceTopo->GetFaceType();
    faceType.resize(nTotalFace);
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        int nodeNumber  = nodeNumberOfEachFace[iFace];
        int iFaceType   = GetFaceType(dimension, nodeNumber);
        faceType[iFace] = iFaceType;
    }

    PHVectorInt1D &parentFaceIndex = faceTopo->GetParentFaceIndex();
    parentFaceIndex.resize(nTotalFace);
    SetField(parentFaceIndex, PHSPACE::INVALID_INDEX, nTotalFace);

    PHVectorInt2D &childrenFaceIndex = faceTopo->GetChildrenFaceIndex();
    childrenFaceIndex.resize(nTotalFace);

    PHVectorInt1D& childrenToCompositeFaceIndex = faceTopo->GetChildrenToCompositeFaceIndex();
    childrenToCompositeFaceIndex.resize(nTotalFace);
}

void Mesh_Refine_Unstruct::ConstructInterFaceTopo(Grid *grid)
{
    UnstructGrid *unstructuredGrid = UnstructGridCast(grid);

    int localIndex = grid->GetZoneLocalID();
    Mesh_RefineController *currentController = gridRefineControllers[localIndex];

    currentController->SetInterfaceInformation(unstructuredGrid->GetInterfaceInfo());
}

void Mesh_Refine_Unstruct::RefineGrid()
{
    WriteLogFile("Start Grid Adaptive!");
    PrintToWindow("Start Grid Adaptive!\n");

    PH_Barrier();
    TimeSpan timeSpan;

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        Mesh_RefineController *currentController = gridRefineControllers[iZone];
        currentController->RefineGrid();
    }

    PH_Barrier();
    WriteLogFile("Time consuming of grid refinement: ", timeSpan.GetTimeSpan());

    WriteLogFile("End Grid Adaptive!");
    PrintToWindow("End Grid Adaptive!\n");
}

#ifdef USE_INHOUSE_VERSION
#ifdef USE_OCC
void Mesh_Refine_Unstruct::ShapePreserving()
{
    string geometryFileName;
    GlobalDataBase::GetData("geometryFileName", &geometryFileName, PHSTRING, 1);
    int isProject = GlobalDataBase::GetIntParaFromDB("isProject");
    int isDeform = GlobalDataBase::GetIntParaFromDB("isDeform");
    int numberOfReferenceCP = GlobalDataBase::GetIntParaFromDB("numberOfReferenceCP");
    int symmetryPlane = GlobalDataBase::GetIntParaFromDB("symmetryPlane");
    RDouble influenceRadius = GlobalDataBase::GetDoubleParaFromDB("influenceRadius");

    PointVector *zoneProjectVector = new PointVector;
    zoneProjectVector->resize(numberOfZones);

    ParallelProjection *projection = new ParallelProjection(geometryFileName, GetDim());
    if (isProject)
    {
        PH_Barrier();
        TimeSpan timeSpan;

        int geometryUnit;
        GlobalDataBase::GetData("geometryUnit", &geometryUnit, PHINT, 1);
        int exclusiveCase;
        GlobalDataBase::GetData("exclusiveCase", &exclusiveCase, PHINT, 1);
        int projectOrgPoint;
        GlobalDataBase::GetData("projectOrgPoint", &projectOrgPoint, PHINT, 1);

        projection->SetZoneGridSolverContainer(numberOfZones, this->zoneGridSolverContainer);
        projection->SetSymmetryPlane(symmetryPlane);
        projection->SetGeometryUnit(geometryUnit);
        projection->SetExclusiveCase(exclusiveCase);
        projection->SetIfProjectOriginalPoint(projectOrgPoint);
        projection->SetProjectVectorArray(zoneProjectVector);
        projection->SetRigion(this->region);
        projection->Init();
        //if(IfNeedProjection() || PHMPI::GetCurrentProcessorID() == PHMPI::GetServerProcessorID())
        {
            projection->Run();
        }

        PH_Barrier();
        WriteLogFile("Time consuming of projection: ", timeSpan.GetTimeSpan());
    }

    ParallelRBF* parallelRBF = new ParallelRBF(GetDim(), influenceRadius);
    if (isProject && isDeform)
    {
        //! Preparing the data that RBF needs.
        int nZones = PHMPI::GetNumberofGlobalZones();
        vector< vector<bool> >* isSolidWallPoint = new vector< vector<bool> >;
        vector< vector<bool> >* isSymmetryPoint = new vector< vector<bool> >;
        vector< vector<bool> >* needDeform = new vector< vector<bool> >;
        PointVector* boundaryDeformedPoints = new PointVector;
        PointVector* deformedVector = new PointVector;
        PointVector* pointsCoord = new PointVector;
        PrepareCheckIsSymmetryPoint(isSymmetryPoint);
        PrepareCheckIsSolidWallPoint(isSolidWallPoint, needDeform);
        PrepareBoundaryDeformedPoint(isSolidWallPoint, zoneProjectVector, boundaryDeformedPoints, deformedVector);
        PreparePointsCoordinates(pointsCoord);

        //! RBF begin.
        PH_Barrier();
        TimeSpan timeSpan;

        parallelRBF->SetNumberOfZone(numberOfZones, nZones);
        parallelRBF->SetPointsCoordinates(pointsCoord);
        parallelRBF->SetBoundaryDeformedPoints(boundaryDeformedPoints);
        parallelRBF->SetDeformedVector(deformedVector);
        parallelRBF->SetIsSymmetryPoint(isSymmetryPoint);
        parallelRBF->SetIsPointsNeedDeform(needDeform);
        parallelRBF->SetSymmetryPlane(symmetryPlane);
        parallelRBF->SetNumberOfReferenceCP(numberOfReferenceCP);
        parallelRBF->SetIsParallelMeshGen(true);

        parallelRBF->Run();

        PH_Barrier();
        WriteLogFile("Time consuming of parallel RBF: ", timeSpan.GetTimeSpan());

        //! Update the points after RBF deformation.
        UpdatePointsCoordinates(pointsCoord, needDeform);

        delete isSymmetryPoint; delete isSolidWallPoint; delete needDeform;
        delete boundaryDeformedPoints; delete deformedVector;
        delete pointsCoord;
    }

    delete parallelRBF;
    delete projection;

    delete zoneProjectVector;
}

void Mesh_Refine_Unstruct::PrepareCheckIsSymmetryPoint(vector< vector<bool> > *isSymmetryPoint)
{
    int nLocalZones = PHMPI::GetNumberofLocalZones();
    isSymmetryPoint->resize(nLocalZones);
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        Mesh_RefineController *currentController = gridRefineControllers[iZone];
        Mesh_FaceTopo  *faceTopo = currentController->GetFaceFactory();
        Geo_PointFactory *nodeTopo = currentController->GetNodeFactory();

        int numberOfNodes = static_cast<int>(nodeTopo->GetNumberOfComputationalNodes());
        int numberOfTotalFaces = faceTopo->GetNTotalFaces();
        PHVectorInt2D      &face2Node             = faceTopo->GetFace2Node();
        vector<SimpleBC* > &faceBoundaryCondition = faceTopo->GetFaceBoundaryCondition();

        (*isSymmetryPoint)[iZone].clear();
        (*isSymmetryPoint)[iZone].resize(numberOfNodes, false);

        for (int iFace = 0; iFace < numberOfTotalFaces; ++ iFace)
        {
            int bcType = faceBoundaryCondition[iFace]->GetBCType();
            if (bcType != PHENGLEI::SYMMETRY) continue;

            for (unsigned iNode = 0; iNode < face2Node[iFace].size(); ++ iNode)
            {
                int nodeIndex = face2Node[iFace][iNode];
                (*isSymmetryPoint)[iZone][nodeIndex] = true;
            }
        }
    }
}

void Mesh_Refine_Unstruct::PrepareCheckIsSolidWallPoint(vector< vector<bool> > *isSolidWallPoint, vector< vector<bool> > *needDeform)
{
    int nLocalZones = PHMPI::GetNumberofLocalZones();
    isSolidWallPoint->resize(nLocalZones);
    needDeform->resize(nLocalZones);
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        Mesh_RefineController *currentController = gridRefineControllers[iZone];
        Mesh_FaceTopo  *faceTopo = currentController->GetFaceFactory();
        Geo_PointFactory *nodeTopo = currentController->GetNodeFactory();

        int numberOfNodes = static_cast<int>(nodeTopo->GetNumberOfComputationalNodes());
        int numberOfTotalFaces = faceTopo->GetNTotalFaces();
        PHVectorInt2D      &face2Node             = faceTopo->GetFace2Node();
        vector<SimpleBC* > &faceBoundaryCondition = faceTopo->GetFaceBoundaryCondition();

        (*isSolidWallPoint)[iZone].clear();
        (*isSolidWallPoint)[iZone].resize(numberOfNodes, false);

        (*needDeform)[iZone].clear();
        (*needDeform)[iZone].resize(numberOfNodes, true);

        for (int iFace = 0; iFace < numberOfTotalFaces; ++ iFace)
        {
            int bcType = faceBoundaryCondition[iFace]->GetBCType();
            if (bcType != PHENGLEI::SOLID_SURFACE) continue;

            for (unsigned iNode = 0; iNode < face2Node[iFace].size(); ++ iNode)
            {
                int nodeIndex = face2Node[iFace][iNode];
                (*isSolidWallPoint)[iZone][nodeIndex] = true;
            }
        }

        for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
        {
            if ((*isSolidWallPoint)[iZone][iNode])
            {
                (*needDeform)[iZone][iNode] = false;
            }
        }
    }
}

void Mesh_Refine_Unstruct::PrepareBoundaryDeformedPoint(vector< vector<bool> > *isSolidWallPoint, PointVector *zoneProjectVector, PointVector *boundaryDeformedPoints, PointVector *deformedVector)
{
    int nLocalZones = PHMPI::GetNumberofLocalZones();
    boundaryDeformedPoints->resize(nLocalZones);
    deformedVector->resize(nLocalZones);
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        Mesh_RefineController *currentController = gridRefineControllers[iZone];
        Geo_PointFactory *nodeTopo = currentController->GetNodeFactory();

        int numberOfNodes = static_cast<int>(nodeTopo->GetNumberOfComputationalNodes());
        int nBoundaryDeformedPoints = 0;
        for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
        {
            if ((*isSolidWallPoint)[iZone][iNode])
            {
                ++ nBoundaryDeformedPoints;
            }
        }
        (*boundaryDeformedPoints)[iZone].resize(nBoundaryDeformedPoints);
        (*deformedVector)[iZone].resize(nBoundaryDeformedPoints);

        nBoundaryDeformedPoints = 0;
        for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
        {
            if ((*isSolidWallPoint)[iZone][iNode])
            {
                double dx = (*zoneProjectVector)[iZone][iNode].X();
                double dy = (*zoneProjectVector)[iZone][iNode].Y();
                double dz = (*zoneProjectVector)[iZone][iNode].Z();

                RDouble x0 = nodeTopo->GetPoint(iNode).X();
                RDouble y0 = nodeTopo->GetPoint(iNode).Y();
                RDouble z0 = nodeTopo->GetPoint(iNode).Z();

                Point3D point1(x0, y0, z0);
                (*boundaryDeformedPoints)[iZone][nBoundaryDeformedPoints] = point1;

                Point3D point2(dx, dy, dz);
                (*deformedVector)[iZone][nBoundaryDeformedPoints] = point2;

                ++ nBoundaryDeformedPoints;
            }
        }
    }
}

void Mesh_Refine_Unstruct::PreparePointsCoordinates(PointVector* pointsCoord)
{
    int nLocalZones = PHMPI::GetNumberofLocalZones();
    pointsCoord->resize(nLocalZones);

    for (int iZone = 0; iZone < nLocalZones; ++iZone)
    {
        Mesh_RefineController *currentController = gridRefineControllers[iZone];
        Geo_PointFactory *nodeTopo = currentController->GetNodeFactory();

        int numberOfNodes = static_cast<int>(nodeTopo->GetNumberOfComputationalNodes());

        (*pointsCoord)[iZone].clear();
        (*pointsCoord)[iZone].resize(numberOfNodes);

        for (int iNode = 0; iNode < numberOfNodes; ++iNode)
        {
            double x0 = nodeTopo->GetPoint(iNode).X();
            double y0 = nodeTopo->GetPoint(iNode).Y();
            double z0 = nodeTopo->GetPoint(iNode).Z();

            (*pointsCoord)[iZone][iNode].SetX(x0);
            (*pointsCoord)[iZone][iNode].SetY(y0);
            (*pointsCoord)[iZone][iNode].SetZ(z0);
        }
    }
}

void Mesh_Refine_Unstruct::UpdatePointsCoordinates(PointVector* pointsCoord, vector< vector<bool> >* needDeform)
{
    int nLocalZones = PHMPI::GetNumberofLocalZones();
    pointsCoord->resize(nLocalZones);

    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        Mesh_RefineController *currentController = gridRefineControllers[iZone];
        Geo_PointFactory *nodeTopo = currentController->GetNodeFactory();

        int numberOfNodes = static_cast<int>(nodeTopo->GetNumberOfComputationalNodes());
        for (int iNode = 0; iNode < numberOfNodes; ++iNode)
        {
            if ((*needDeform)[iZone][iNode])
            {
                double xNew = (*pointsCoord)[iZone][iNode].X();
                double yNew = (*pointsCoord)[iZone][iNode].Y();
                double zNew = (*pointsCoord)[iZone][iNode].Z();

                nodeTopo->GetPoint(iNode).SetX(xNew);
                nodeTopo->GetPoint(iNode).SetY(yNew);
                nodeTopo->GetPoint(iNode).SetZ(zNew);
            }
        }
    }
}
#endif
#endif

void Mesh_Refine_Unstruct::FixAnisoRefineType()
{
    InitAnisoRefineType();

    if (GetDim() == PHSPACE::TWO_D) return;

    int anisoRefine = this->refineParameter->GetAnisoRefineType();
    if (anisoRefine != ANISOTROPICREFINE) return;

    WriteLogFile("Start to fix cell type!");
    TimeSpan fixTime;

    int infectTime = 0;
    int isContinueInfect = 1;
    while (isContinueInfect)
    {
        isContinueInfect = 0;

        if (infectTime > 0)
        {
            PrintToWindow("\nInfect isotropic refine type time: ", infectTime, "......");
            WriteLogFile("\nInfect isotropic refine type time: ", infectTime, "......");
        }

        for (int iZone = 0; iZone < numberOfZones; ++ iZone)
        {
            Mesh_RefineController *currentController = gridRefineControllers[iZone];
            currentController->FixAnisoRefineType();
        }

        SwapInterfaceRefineType();

        isContinueInfect = IsNeedInfectInterfaceIsotropicRefineType();

        PH_CompareMaxMin(isContinueInfect, 1);

        ++ infectTime;
    }

    PrintToWindow("\n");

    fixTime.ShowTimeSpanToLogFile("Fix cell refine type");
    WriteLogFile("End fix cell type!");
}

bool Mesh_Refine_Unstruct::IsNeedInfectInterfaceIsotropicRefineType()
{
    bool isContinue = false;

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        Mesh_RefineController *currentController = gridRefineControllers[iZone];

        if (currentController->IsNeedInfectInterfaceIsotropicRefineType())
        {
            isContinue = true;
            break;
        }
    }

    return isContinue;
}

void Mesh_Refine_Unstruct::SwapInterfaceRefineType()
{
    WriteLogFile("Start Swapping refine type!");

    int currentProcessor = PHMPI::GetCurrentProcessorID();
    int nZones = PHMPI::GetNumberofGlobalZones();
    vector< PH_Request > requestContainer;
    vector< vector< DataContainer* > > receivedDataBuffer;
    vector< vector< DataContainer* > > sendDataBuffer;
    receivedDataBuffer.resize(nZones);
    sendDataBuffer.resize(nZones);

    //! Compress the send information into the actkey.
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int sendProcessor = PHMPI::GetZoneProcessorID(iZone);

        ZoneNeighbor *globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
        std::size_t numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();
        if (currentProcessor == sendProcessor)
        {
            sendDataBuffer[iZone].resize(numberOfNeighbor);
        }

        for (std::size_t iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
        {
            int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);

            int receiveProcessor = PHMPI::GetZoneProcessorID(neighborZone);

            if (currentProcessor != sendProcessor && currentProcessor != receiveProcessor)
            {
                continue;
            }

            DataContainer *sendBuffer = 0;
            if (currentProcessor == sendProcessor)
            {
                sendBuffer = new DataContainer;
                sendDataBuffer[iZone][iNeighbor] = sendBuffer;

                int localZoneIndex = region->GetProcessGlobalZoneIndexToLocalZoneIndex(iZone);
                Mesh_RefineController *currentController = gridRefineControllers[localZoneIndex];
                currentController->CompressRefineType(sendBuffer, neighborZone);
            }
        }
    }

    //! Communication.
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        ZoneNeighbor *globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
        std::size_t numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();

        for (std::size_t iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
        {
            int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);

            int sendProcessor = PHMPI::GetZoneProcessorID(iZone);
            int receiveProcessor = PHMPI::GetZoneProcessorID(neighborZone);

            if (currentProcessor != sendProcessor && currentProcessor != receiveProcessor)
            {
                continue;
            }

            DataContainer *receivedBuffer = 0;
            DataContainer *sendBuffer = 0;
            if (currentProcessor == receiveProcessor)
            {
                receivedBuffer = new DataContainer;
                receivedDataBuffer[iZone].push_back(receivedBuffer);
            }
            if (currentProcessor == sendProcessor)
            {
                sendBuffer = sendDataBuffer[iZone][iNeighbor];
            }

            int tag = iZone;
            if (currentProcessor == sendProcessor)
            {
                if (sendProcessor != receiveProcessor)
                {
                    streamsize nlen = sendBuffer->size();
                    send(sendBuffer, receiveProcessor, nlen, requestContainer, tag);
                }
            }
            if (currentProcessor == receiveProcessor)
            {
                int neighborOrer = -1;
                ZoneNeighbor *globalNeighborZonesTemp = zoneConnectivity->GetZoneNeighbor(neighborZone);
                std::size_t numberOfNeighborTemp = globalNeighborZonesTemp->GetNumberOfNeighbor();

                for (std::size_t iNeighbor = 0; iNeighbor < numberOfNeighborTemp; ++ iNeighbor)
                {
                    uint_t neighborZoneTemp = globalNeighborZonesTemp->GetZoneIndexOfNeighbor(iNeighbor);
                    if (neighborZoneTemp == iZone)
                    {
                        neighborOrer = static_cast<int>(iNeighbor);
                        break;
                    }
                }

                ASSERT(neighborOrer != -1);

                if (sendProcessor != receiveProcessor)
                {
                    CharVecSizeType nlen = sendDataBuffer[neighborZone][neighborOrer]->size();
                    receive(receivedBuffer, sendProcessor, nlen, requestContainer, tag);
                }
                else
                {
                    uint_t lastID = receivedDataBuffer[iZone].size();
                    SWAP(receivedDataBuffer[iZone][lastID - 1], sendDataBuffer[neighborZone][neighborOrer]);
                }
            }
        }
    }

    if (PHMPI::GetNumberOfProcessor() > 1)
    {
        PH_Wait(static_cast<int>(requestContainer.size()), &(requestContainer[0]));
    }

    //! Decompress the interface data from datacontainer.
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        ZoneNeighbor *globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
        std::size_t numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();

        int count = 0;
        for (std::size_t iNeighbor = 0; iNeighbor < numberOfNeighbor; ++iNeighbor)
        {
            int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);

            int receiveProcessor = PHMPI::GetZoneProcessorID(neighborZone);
            if (currentProcessor == receiveProcessor)
            {
                DataContainer *receiveData = receivedDataBuffer[iZone][count];

                int localZoneIndex = region->GetProcessGlobalZoneIndexToLocalZoneIndex(neighborZone);
                Mesh_RefineController *currentController = gridRefineControllers[localZoneIndex];
                currentController->DecompressRefineType(receiveData, iZone);

                ++ count;
            }
        }
    }

    for (std::size_t iDim = 0; iDim < receivedDataBuffer.size(); ++ iDim)
    {
        for (std::size_t jDim = 0; jDim < receivedDataBuffer[iDim].size(); ++ jDim)
        {
            delete receivedDataBuffer[iDim][jDim];
        }
    }
    receivedDataBuffer.clear();
    for (std::size_t iDim = 0; iDim < sendDataBuffer.size(); ++ iDim)
    {
        for (std::size_t jDim = 0; jDim < sendDataBuffer[iDim].size(); ++ jDim)
        {
            delete sendDataBuffer[iDim][jDim];
        }
    }
    sendDataBuffer.clear();
    requestContainer.clear();

    WriteLogFile("End Swapping refine type!");
}

void Mesh_Refine_Unstruct::InitAnisoRefineType()
{
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        Mesh_RefineController *currentController = gridRefineControllers[iZone];
        currentController->InitAnisoRefineType();
    }
}

void Mesh_Refine_Unstruct::BuildRefineProperty()
{
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        SetDefaultRefineProperty(iZone);

        ComputeRefineProperty(iZone);

        ResetCellProperty(iZone);
    }
}

void Mesh_Refine_Unstruct::ResetCellProperty(int iZone)
{
    Mesh_RefineController *currentController = gridRefineControllers[iZone];
    Mesh_CellTopo *cellTopo = currentController->GetCellFactory();

    PHVectorInt1D &cellComputationalStatus = cellTopo->GetCellComputationalStatus();
    PHVectorInt1D &cellModifiedStatus = cellTopo->GetCellModifiedStatus();

    SetField(cellModifiedStatus, SPLIT, static_cast<int>(cellModifiedStatus.size()));
    SetField(cellComputationalStatus, HIDDEN, static_cast<int>(cellComputationalStatus.size()));
}

void Mesh_Refine_Unstruct::ComputeRefineProperty(int iZone)
{
    PrepareDataStructure(iZone);

    SetCellModifiedStatusToNoChange(iZone);

    SetCellComputationalStatusToHiddenState(iZone);

    SetCurrentCellModifiedStatusByOldCellModifiedStatus(iZone);

    SetCurrentCellComputationalStatus(iZone);
}

void Mesh_Refine_Unstruct::PrepareDataStructure(int iZone)
{
    Mesh_RefineController *currentController = gridRefineControllers[iZone];

    currentController->Initialize();
    currentController->ConstructNode2Cell();
    currentController->ConstructCell2Cell();
}

void Mesh_Refine_Unstruct::SetCellModifiedStatusToNoChange(int iZone)
{
    Mesh_RefineController *currentController = gridRefineControllers[iZone];
    Mesh_CellTopo* cellTopo = currentController->GetCellFactory();

    int numberOfTotal = cellTopo->GetNTotalCell() + currentController->GetNBoundaryFace();

    PHVectorInt1D &cellModifiedStatus = cellTopo->GetCellModifiedStatus();
    cellModifiedStatus.resize(numberOfTotal);
    SetField(cellModifiedStatus, NOCHANGE, numberOfTotal);
}

void Mesh_Refine_Unstruct::SetCellComputationalStatusToHiddenState(int iZone)
{
    Mesh_RefineController *currentController = gridRefineControllers[iZone];
    Mesh_CellTopo *cellTopo = currentController->GetCellFactory();

    int numberOfTotal = cellTopo->GetNTotalCell() + currentController->GetNBoundaryFace();

    PHVectorInt1D &cellComputationalStatus = cellTopo->GetCellComputationalStatus();
    cellComputationalStatus.resize(numberOfTotal);
    SetField(cellComputationalStatus, HIDDEN, numberOfTotal);
}

void Mesh_Refine_Unstruct::SetCurrentCellModifiedStatusByOldCellModifiedStatus(int iZone)
{
    Mesh_RefineController *currentController = gridRefineControllers[iZone];
    Mesh_CellTopo *cellTopo = currentController->GetCellFactory();

    PHVectorInt1D &computationalCellModifiedStatus = cellTopo->GetComputationalCellModifiedStatus();
    PHVectorInt1D &computationalToGlobalCellIndex = cellTopo->GetComputationalToGlobalCellIndex();

    int oldNumberOfCell = currentController->GetOldNumberOfCell();
    for (int iCell = 0; iCell < oldNumberOfCell; ++ iCell)
    {
        int cellIndex = computationalToGlobalCellIndex[iCell];
        int cellModifiedStatus = computationalCellModifiedStatus[iCell];
        cellTopo->SetCellModifiedStatus(cellIndex, cellModifiedStatus);
    }
}

void Mesh_Refine_Unstruct::SetCurrentCellComputationalStatus(int iZone)
{
    Mesh_RefineController *currentController = gridRefineControllers[iZone];
    Mesh_CellTopo *cellTopo = currentController->GetCellFactory();

    PHVectorInt1D &computationalCellModifiedStatus = cellTopo->GetComputationalCellModifiedStatus();
    PHVectorInt1D &computationalToGlobalCellIndex = cellTopo->GetComputationalToGlobalCellIndex();

    int oldNumberOfCell = currentController->GetOldNumberOfCell();
    for (int iCell = 0; iCell < oldNumberOfCell; ++ iCell)
    {
        int cellIndex = computationalToGlobalCellIndex[iCell];
        int cellModifiedStatus = computationalCellModifiedStatus[iCell];

        if (cellModifiedStatus == NOCHANGE)
        {
            cellTopo->SetComputationalStatusOfCell(cellIndex, ON);
        }
    }
}

void Mesh_Refine_Unstruct::SetDefaultRefineProperty(int iZone)
{
    Mesh_RefineController *currentController = gridRefineControllers[iZone];
    Mesh_CellTopo *cellTopo = currentController->GetCellFactory();

    int nTotalCell = cellTopo->GetNTotalCell();
    PHVectorInt1D &computationalCellModifiedStatus = cellTopo->GetComputationalCellModifiedStatus();
    computationalCellModifiedStatus.resize(nTotalCell);
    SetField(computationalCellModifiedStatus, NOCHANGE, nTotalCell);
}

int Mesh_Refine_Unstruct::GetCellType(int dimension, int nodeNumberOfCell)
{
    int cellType = -1;
    if (dimension == TWO_D)
    {
        switch (nodeNumberOfCell)
        {
        case 3:
            cellType = TRI_3;
            break;
        case 4:
            cellType = QUAD_4;
            break;
        default:
            break;
        }
    }
    else
    {
        switch (nodeNumberOfCell)
        {
        case 4:
            cellType = TETRA_4;
            break;
        case 5:
            cellType = PYRA_5;
            break;
        case 6:
            cellType = PENTA_6;
            break;
        case 8:
            cellType = HEXA_8;
            break;
        default:
            break;
        }
    }

    if (cellType == -1)
    {
        TK_Exit::PrintDebugInfoExit("This cell type is not supported!");
    }

    return cellType;
}

int Mesh_Refine_Unstruct::GetFaceType(int dimension, int nodeNumberOfFace)
{
    int faceType = -1;
    if (dimension == TWO_D)
    {
        switch (nodeNumberOfFace)
        {
        case 2:
            faceType = BAR_2;
            break;
        case 3:
            faceType = BAR_3;
            break;
        default:
            break;
        }
    }
    else
    {
        switch (nodeNumberOfFace)
        {
        case 3:
            faceType = TRI_3;
            break;
        case 4:
            faceType = QUAD_4;
            break;
        default:
            break;
        }
    }

    if (faceType == -1)
    {
        TK_Exit::PrintDebugInfoExit("This face type is not supported!");
    }

    return faceType;
}

}

