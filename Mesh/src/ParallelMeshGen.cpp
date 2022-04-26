#include "ParallelMeshGen.h"
#include "Zone.h"
#include "GridAdaptationCGNS.h"
#include "GridAdaptation.h"
#include "Glb_Dimension.h"
#include "Geometry.h"
#include "ComputationalGrid.h"
#include "PHIO.h"
#include "IO_FileName.h"
#include "Geo_Element.h"
#include "TK_Time.h"
#include "Pre_HDF5File.h"
#include "Pre_GridConversion.h"

using namespace std;
#pragma warning (disable:913)

namespace PHSPACE
{
void CompressVectorIntoDataContainer(DataContainer *& dataContainer, PHVector1D < int > & data)
{
    int size = static_cast<int>(data.size());
    PHWrite(dataContainer, &size, 1);

    PHWrite(dataContainer, data, static_cast<int>(data.size()));
}

void AnasysVectorFromDataContainer(DataContainer *& dataContainer, PHVector1D < int > & data)
{
    if(dataContainer->size() == 0) return;

    dataContainer->moveto_begin();

    int size;
    PHRead(dataContainer, &size, 1);
    data.resize(size);

    PHRead(dataContainer, data, size);
}

ParallelMeshGen::ParallelMeshGen()
{
    numberOfZones = 0;
    super_elements = 0;
    region   = 0;
    zoneGridSolverContainer = 0;
    currentZoneGridSolver = 0;
    fantasyGrids = 0;
    cgns_multizone = NULL;
    //computationalElementSolver = new ComputationalElementSolver();
}

ParallelMeshGen::~ParallelMeshGen()
{
    delete region;

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        delete zoneGridSolverContainer[iZone];
    }
    delete []zoneGridSolverContainer;
}

void ParallelMeshGen::ParallelGridGen()
{
    ReadZoneGridFromFTS();

    InitGlobalBoundaryCondition();

    AllocateZoneGridSolverContainer();

    ConstructZoneGridSolver();

    // 确定每个单元的自适应状态
    BuildElementAdaptProperty();

    FixCellRefineType();

    AdaptGrid();

    ShapePreserving();

    GenerateAndDumpComputationalGrid();
}

void ParallelMeshGen::ReadZoneGridFromFTS()
{
    string original_grid_file;
    GlobalDataBase::GetData("original_grid_file", &original_grid_file, PHSTRING, 1);
    GlobalDataBase::UpdateData("gridfile", &original_grid_file, PHSTRING, 1);

    string adapt_grid_file;
    GlobalDataBase::GetData("adapt_grid_file", &adapt_grid_file, PHSTRING, 1);

    region = new Region();

    // Step1: Read in the original grid.
    region->ReadGrid();

    region->SetupMultiBlockRelationship();
}

void ParallelMeshGen::AllocateZoneGridSolverContainer()
{
    numberOfZones = PHMPI::GetNumberofLocalZones();

    zoneGridSolverContainer = new ZoneGridSolver * [numberOfZones];

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        PHVector1D < CgnsZone * > cgnsZoneContainer;
        zoneGridSolverContainer[iZone] = new ZoneGridSolver(cgnsZoneContainer);
        zoneGridSolverContainer[iZone]->SetZoneHandler(this->region);
    }

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        zoneGridSolverContainer[iZone]->SetZoneIndex(iZone);
    }
}

void ParallelMeshGen::ConstructZoneGridSolver()
{
    WriteLogFile("Start Constructing ZoneGrid Solver!");

    VisulizationOriginalGrid();

    RecoverZoneGridSolverFromFTS();

    /*
    if(PHMPI::GetNumberOfProcessor() > 1)
        BroadcastNeighborZonesToEachProcessorNew();
    else
        BroadcastNeighborZonesToEachProcessor();
    //*/
    
    
    //region->SetupMultiBlockRelationship();

    /*
    InitialInterFaceInformation();
    
    SwapNeighborsSendContent();
    //*/
    
    //ComputeFaceRBFWeight();
    //*/
    
}

void ParallelMeshGen::BuildElementAdaptProperty()
{
    int nLocalZones = PHMPI::GetNumberofLocalZones();
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        ZoneGridSolver * zoneGridSolver = zoneGridSolverContainer[ iZone ];
        this->SetCurrentZoneGridSolver(zoneGridSolver);
        //this->SetCurrentAdaptiveFlowFieldSolver(adaptiveFlowFieldSolverVector[ iZone ]);

        this->DoSomethingIfResetElementPropertyBySelfDefine();
        //this->ReadGlobalOldQField();
        this->ReComputeAdaptionProperty();//根据判则重新判断哪些单元需要自适应
        this->ResetElementPropertyBySelfDefine();
    }
}

void ParallelMeshGen::ReComputeAdaptionProperty()
{
    AdaptivePropertySolver * adaptivePropertySolver = new AdaptivePropertySolver();

    adaptivePropertySolver->SetZoneGridSolver(currentZoneGridSolver);
    adaptivePropertySolver->ReComputeAdaptionProperty();

    delete adaptivePropertySolver;
}

void ParallelMeshGen::DoSomethingIfResetElementPropertyBySelfDefine()
{
    int adaptGridBySelfDefineMethod = 1;//GetAdaptGridBySelfDefineMethod();
    if(!adaptGridBySelfDefineMethod) return;

    //PHAdaptiveFlowFieldSolver * adaptiveFlowFieldSolver = this->GetCurrentAdaptiveFlowFieldSolver();
    //adaptiveFlowFieldSolver->SetNumberOfComputationalElements(numberOfComputationalElementsOfCurrentZone);

    //int numberOfTotalEquations = 5;
    //PHSPACE::SetNumberOfNavierStokesEquations(numberOfTotalEquations);

    int numberOfElements = currentZoneGridSolver->GetNumberOfElements();

    ElementPropertySolver * elementPropertySolver = currentZoneGridSolver->GetElementPropertySolver();
    PHVectorInt1D & computationalElementModifiedStatusContainer = elementPropertySolver->GetComputationalElementModifiedStatusContainer();
    computationalElementModifiedStatusContainer.resize(numberOfElements);
    SetField(computationalElementModifiedStatusContainer, PHSPACE::NOCHANGE, numberOfElements);
}

void ParallelMeshGen::ResetElementPropertyBySelfDefine()
{
    int adaptGridBySelfDefineMethod = 1;//GetAdaptGridBySelfDefineMethod();
    if(!adaptGridBySelfDefineMethod) return;

    //PHAdaptiveFlowFieldSolver * adaptiveFlowFieldSolver = currentAdaptiveFlowFieldSolver;
    //adaptiveFlowFieldSolver->ResetElementPropertyBySelfDefine();

    ElementPropertySolver * elementPropertySolver = currentZoneGridSolver->GetElementPropertySolver();

    PHVectorInt1D & elementComputationalStatusContainer = elementPropertySolver->GetElementComputationalStatusContainer();

    PHVectorInt1D & elementModifiedStatusContainer = elementPropertySolver->GetElementModifiedStatusContainer();

    if (adaptGridBySelfDefineMethod == 1)
    {
        // All Cells Split.
        SetField(elementModifiedStatusContainer, PHSPACE::SPLIT, static_cast<int>(elementModifiedStatusContainer.size()));
        SetField(elementComputationalStatusContainer, PHSPACE::HIDDEN, static_cast<int>(elementComputationalStatusContainer.size()));
    }
    else if (adaptGridBySelfDefineMethod == 2)
    {
        // All Cells No Change.
        SetField(elementModifiedStatusContainer, PHSPACE::NOCHANGE, static_cast<int>(elementModifiedStatusContainer.size()));
        SetField(elementComputationalStatusContainer, PHSPACE::ON, static_cast<int>(elementComputationalStatusContainer.size()));
    }
    else if (adaptGridBySelfDefineMethod == 3)
    {
        int zoneIndex = currentZoneGridSolver->GetZoneIndex();
        int globalZoneIndex = region->GetProcessLocalZoneIndexToGlobalZoneIndex(zoneIndex);
        if (globalZoneIndex % 2 == 1)
        {
            // All Cells Split.
            SetField(elementModifiedStatusContainer, PHSPACE::SPLIT, static_cast<int>(elementModifiedStatusContainer.size()));
            SetField(elementComputationalStatusContainer, PHSPACE::HIDDEN, static_cast<int>(elementComputationalStatusContainer.size()));
        }
    }
    else if (adaptGridBySelfDefineMethod == 4)
    {
        int zoneIndex = currentZoneGridSolver->GetZoneIndex();
        int globalZoneIndex = region->GetProcessLocalZoneIndexToGlobalZoneIndex(zoneIndex);
        if (globalZoneIndex % 2 == 0)
        {
            // All Cells Split.
            SetField(elementModifiedStatusContainer, PHSPACE::SPLIT, static_cast<int>(elementModifiedStatusContainer.size()));
            SetField(elementComputationalStatusContainer, PHSPACE::HIDDEN, static_cast<int>(elementComputationalStatusContainer.size()));
        }
    }
    else if (adaptGridBySelfDefineMethod == 5)
    {
        elementModifiedStatusContainer[ 0 ]      = PHSPACE::SPLIT;
        elementComputationalStatusContainer[ 0 ] = PHSPACE::HIDDEN;

        elementModifiedStatusContainer[ 1 ]      = PHSPACE::NOCHANGE;
        elementComputationalStatusContainer[ 1 ] = PHSPACE::ON;
    }
}


void ParallelMeshGen::AdaptGrid()
{
    WriteLogFile("Start Grid Adaptive!");
    PrintToWindow("Start Grid Adaptive!\n");

    PH_Barrier();
    TimeSpan timeSpan;

    int nLocalZones = PHMPI::GetNumberofLocalZones();
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        ZoneGridSolver * zoneGridSolver = zoneGridSolverContainer[ iZone ];
        this->SetCurrentZoneGridSolver(zoneGridSolver);
        //this->SetCurrentAdaptiveFlowFieldSolver(adaptiveFlowFieldSolverVector[ iZone ]);

        AdaptGrid(zoneGridSolver);
    }

    PH_Barrier();
    WriteLogFile("Time consuming of grid refinement: ", timeSpan.GetTimeSpan());

    WriteLogFile("End Grid Adaptive!");
    PrintToWindow("End Grid Adaptive!\n");
}

void ParallelMeshGen::AdaptGrid(ZoneGridSolver * zoneGridSolver)
{
    this->GridAdaptation(zoneGridSolver);

    this->ConstructBasicData();

    //this->DumpAdaptiveData();
}


void ParallelMeshGen::GridAdaptation(ZoneGridSolver * zoneGridSolver)
{
    WriteLogFile("Start Splitting And Merging Grid!");

    /*SetAdaptTropic();*/

    int numberOfGridSelfAdaptions = 1;//PHSPACE::GetNumberOfGridSelfAdaptions();

    for (int iGridSelfAdaption = 0; iGridSelfAdaption < numberOfGridSelfAdaptions; ++ iGridSelfAdaption)
    {
        this->SplitAllElements();

        //this->MergeElements();
    }

    currentZoneGridSolver->SetAdaptStatus(true);

    //cout << " Spliting and combining element complete \n";

}

void ParallelMeshGen::SplitAllElements()
{
    ElementPropertySolver * elementPropertySolver = currentZoneGridSolver->GetElementPropertySolver();

    Geo_PointFactory * pointFactory = currentZoneGridSolver->GetPointFactory();

    elementPropertySolver->SplitAllElements(pointFactory);
}

void ParallelMeshGen::ConstructBasicData()
{
    WriteLogFile("Start Constructing Basic Data!");

    ElementPropertySolver * elementPropertySolver = currentZoneGridSolver->GetElementPropertySolver();

    WriteLogFile("Start ScanAllElementsToReconstructFaceDataStructure!");
    elementPropertySolver->ScanAllElementsToReconstructFaceDataStructure(currentZoneGridSolver->GetAdaptStatus());

    WriteLogFile("Start GenerateElementIndexMapping!");
    this->GenerateElementIndexMapping();

    ElementFaceSolver     * elementFaceSolver     = currentZoneGridSolver->GetElementFaceSolver();

    WriteLogFile("Start GenerateGlobalToComputationalElementFaceIndexMapping!");
    elementFaceSolver->GenerateGlobalToComputationalElementFaceIndexMapping();

    WriteLogFile("Start GenerateComputationalElementFaceInformation!");
    elementFaceSolver->GenerateComputationalElementFaceInformation();

    WriteLogFile("Start GenerateComputationalFaceNodeIndex!");
    elementFaceSolver->GenerateComputationalFaceNodeIndex();

    //GenerateComputationalFaceNodeIndex();
    //在此没有交换面顺序（让前面为边界face,后面为内部face），否则会逻辑错乱
    //后面GenerateComputationalUnstructuredGrid，即转换为fts网格时再处理
}

void ParallelMeshGen::ShapePreserving()
{
#ifdef USE_OCC
    string geometryFileName;
    GlobalDataBase::GetData("geometryFileName", &geometryFileName, PHSTRING, 1);
    int isProject           = GlobalDataBase::GetIntParaFromDB("isProject");
    int isDeform            = GlobalDataBase::GetIntParaFromDB("isDeform");
    int numberOfReferenceCP = GlobalDataBase::GetIntParaFromDB("numberOfReferenceCP");
    int symmetryPlane       = GlobalDataBase::GetIntParaFromDB("symmetryPlane");
    RDouble influenceRadius = GlobalDataBase::GetDoubleParaFromDB("influenceRadius");

    PointVector *zoneProjectVector = new PointVector;
    zoneProjectVector->resize(numberOfZones);

    ParallelProjection * projection = new ParallelProjection(geometryFileName, GetDim());
    if(isProject)
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
    
    ParallelRBF * parallelRBF = new ParallelRBF(GetDim(), influenceRadius);
    if(isProject && isDeform)
    {
        //! Preparing the data that RBF needs.
        int nZones = PHMPI::GetNumberofGlobalZones();
        vector< vector<bool> > *isSolidWallPoint = new vector< vector<bool> >;
        vector< vector<bool> > *isSymmetryPoint  = new vector< vector<bool> >;
        vector< vector<bool> > *needDeform       = new vector< vector<bool> >;
        PointVector *boundaryDeformedPoints = new PointVector;
        PointVector *deformedVector = new PointVector;
        PointVector *pointsCoord = new PointVector; 
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
#endif
}


void ParallelMeshGen::GenerateElementIndexMapping()
{
    //cout << " GenerateElementIndexMapping begin \n";

    //当前单元数量,在加密或稀疏化的时候，网格单元数量可能会变化
    ElementPropertySolver * elementPropertySolver = currentZoneGridSolver->GetElementPropertySolver();
    int numberOfTotalElements =  elementPropertySolver->GetNumberOfTotalElements(); 

    PHVectorInt1D & computationalToGlobalElementIndexMapping = elementPropertySolver->GetComputationalToGlobalElementIndexMapping();
    PHVectorInt1D & globalToComputationalElementIndexMapping = elementPropertySolver->GetGlobalToComputationalElementIndexMapping();

    globalToComputationalElementIndexMapping.resize(numberOfTotalElements);

    int computationalElementIndex = 0;
    //从体单元集合中选出计算需要的简单体单元
    for (int elementIndex = 0; elementIndex < numberOfTotalElements; ++ elementIndex)
    {
        int computationalStatusOfElement = elementPropertySolver->GetComputationalStatusOfElement(elementIndex);
        if (computationalStatusOfElement == PHSPACE::ON)
        {
            //是计算需要的体单元则统计
            //生成网格编号并与与体单元编号对应
            ++ computationalElementIndex;
        }
    }

    int numberOfComputationalElements = computationalElementIndex;

    ElementFaceSolver     * elementFaceSolver     = currentZoneGridSolver->GetElementFaceSolver();

    ComputationalFaceSolver * computationalFaceSolver = elementFaceSolver->GetComputationalFaceSolver();

    computationalFaceSolver->SetNumberOfElements(numberOfComputationalElements);

    computationalToGlobalElementIndexMapping.resize(numberOfComputationalElements);

    computationalElementIndex = 0;
    //从体单元集合中选出计算需要的简单体单元
    for (int elementIndex = 0; elementIndex < numberOfTotalElements; ++ elementIndex)
    {
        int computationalStatusOfElement = elementPropertySolver->GetComputationalStatusOfElement(elementIndex);
        if (computationalStatusOfElement == PHSPACE::ON)
        {
            //是计算需要的体单元则统计
            //生成网格编号并与与体单元编号对应
            computationalToGlobalElementIndexMapping[ computationalElementIndex ] = elementIndex;
            globalToComputationalElementIndexMapping[ elementIndex              ] = computationalElementIndex;

            ++ computationalElementIndex;
        }
    }

    //cout << " GenerateElementIndexMapping complete \n";

}

void ParallelMeshGen::VisulizationOriginalGrid()
{
    if(PHMPI::GetNumberOfProcessor() >= 64) return;

    int nLocalZones = PHMPI::GetNumberofLocalZones();
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int globalZoneIndex = region->GetProcessLocalZoneIndexToGlobalZoneIndex(iZone);

        Zone * zone = region->GetZone(globalZoneIndex);
        UnstructGrid * unstructGrid = UnstructGridCast(zone->GetGeometry()->GetGrid());


        ostringstream fileName;
        fileName << "./grid/OriginalGrid" << globalZoneIndex << ".plt";

        if(PHSPACE::GetDim() == PHSPACE::TWO_D)
        {
            VisualizationMesh2D(unstructGrid, fileName.str());
        }
        else
        {
            VisualizationMesh3D(unstructGrid, fileName.str());
        }
    }
}

void ParallelMeshGen::RecoverZoneGridSolverFromFTS()
{
    WriteLogFile("Start Recovering Zone Grid Solver From .fts!");

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int globalZoneIndex = region->GetProcessLocalZoneIndexToGlobalZoneIndex(iZone);

        ZoneGridSolver * zoneGridSolver = zoneGridSolverContainer[ iZone ];
        this->SetCurrentZoneGridSolver(zoneGridSolver);

        Zone * zone = region->GetZone(globalZoneIndex);

        UnstructGrid * unstructGrid = UnstructGridCast(zone->GetGeometry()->GetGrid());
        unstructGrid->ComputeMetrics();
        unstructGrid->SkewnessSummary();

        this->RecoverZoneGridSolverFromFTS(zone);

        // Free some memory.
        // The following memory is retain, because the grids will be used latter.
        //delete zone->GetGeometry();
        //zone->SetGeometry(0);
    }
}

void ParallelMeshGen::RecoverZoneGridSolverFromFTS(Zone * zone)
{
    RecoverElements(zone);

    RecoverPoints(zone);

    RecoverFacess(zone);

    RecoverMetrics(zone);

    RecoverInterFaceInformation(zone);
}

void ParallelMeshGen::RecoverElements(Zone * zone)
{
    int dimension = GetDim();

    PHGeometry * geometry = zone->GetGeometry();
    UnstructGrid * unstructuredGrid = UnstructGridCast(geometry->GetGrid());
    
    int * cellNodeNumberContainer = unstructuredGrid->GetNodeNumberOfEachCell();
    int * cellNodeIndexContainer  = unstructuredGrid->GetCell2Node();
    
    int numberOfCells = unstructuredGrid->GetNTotalCell();

    ElementPropertySolver * elementPropertySolver = currentZoneGridSolver->GetElementPropertySolver();
    
    PHVectorInt1D & oldComputationalToGlobalElementIndexMapping = elementPropertySolver->GetOldComputationalToGlobalElementIndexMapping();
    PHVectorInt1D & oldGlobalToComputationalElementIndexMapping = elementPropertySolver->GetOldGlobalToComputationalElementIndexMapping();
    
    // 计算网格总数
    int numberOfComputationalElements = numberOfCells;

    // 读单元总数：计算单元总数与总的单元总数暂时是一样的
    int numberOfTotalElements = numberOfCells;

    oldComputationalToGlobalElementIndexMapping.resize(numberOfComputationalElements);
    oldGlobalToComputationalElementIndexMapping.resize(numberOfTotalElements);

    // 计算单元与整体单元对应关系：暂时是一样的
    for (int iCell = 0; iCell < numberOfCells; ++ iCell)
    {
        oldComputationalToGlobalElementIndexMapping[ iCell ] = iCell;
        oldGlobalToComputationalElementIndexMapping[ iCell ] = iCell;
    }
    
    // 单元类型
    PHVectorInt1D & elementTypeContainer = elementPropertySolver->GetElementTypeContainer();
    elementTypeContainer.resize(numberOfTotalElements);
    
    for (int iCell = 0; iCell < numberOfCells; ++ iCell)
    {
        int nodeNumberOfCell = cellNodeNumberContainer[ iCell ];

        int cellType      = GetElementTypeOfFTS(dimension, nodeNumberOfCell);

        elementTypeContainer[ iCell ] = cellType;
    }
    
    // 单元级别
    elementPropertySolver->SetAllElementsLevelToZero();

    // 单元parent
    elementPropertySolver->SetAllElementsParentIndexToRootState();

    // 单元点index
    PHVectorInt2D & elementNodeIndexContainer = elementPropertySolver->GetElementNodeIndexContainer(); 
    elementNodeIndexContainer.resize(numberOfCells);
    int count = 0;
    for (int iCell = 0; iCell < numberOfCells; ++ iCell)
    {
        int nodeNumberOfCell = cellNodeNumberContainer[ iCell ];
        elementNodeIndexContainer[ iCell ].resize(nodeNumberOfCell);
        for (int iNode = 0; iNode < nodeNumberOfCell; ++ iNode)
        {
            elementNodeIndexContainer[ iCell ][ iNode ] = cellNodeIndexContainer[ count ++ ];
        }
    }
    //*/

    int anisoRefine = GlobalDataBase::GetIntParaFromDB("anisoRefine");
    if(anisoRefine)
    {
        int * face_number_of_each_cell = unstructuredGrid->GetFaceNumberOfEachCell();
        int ** cell2face = unstructuredGrid->GetCell2Face();

        PHVectorInt2D * element2Face = new PHVectorInt2D;
        element2Face->resize(numberOfCells);
        for (int iCell = 0; iCell < numberOfCells; ++ iCell)
        {
            int faceNumberOfCell = face_number_of_each_cell[ iCell ];
            (*element2Face)[ iCell ].resize(faceNumberOfCell);
            for (int iFace = 0; iFace < faceNumberOfCell; ++ iFace)
            {
                (*element2Face)[ iCell ][ iFace ] = cell2face[iCell][iFace];
            }
        }

        elementPropertySolver->SetElement2Face(element2Face);
    }
    

    //* 以下部分还没有考虑.
    // parent、child关系
    PHVectorInt2D & childrenElementsIndexesContainer = elementPropertySolver->GetChildrenElementsIndexesContainer();

    // 实际parent单元总数
    int numberOfParentElements = 0;

    childrenElementsIndexesContainer.resize(numberOfParentElements);

    //读parent对应child单元
    for (int iParentElement = 0; iParentElement < numberOfParentElements; ++ iParentElement)
    {
        //读parent单元
        int parentElementIndex = 0;
        //PHRead(dataContainer, parentElementIndex);

        //读child单元
        PHVectorInt1D & childrenElementsIndexes = childrenElementsIndexesContainer[ parentElementIndex ]; 

        int  numberOfChildElements = 0;
        //PHRead(dataContainer, numberOfChildElements);

        if (numberOfChildElements != 0) //Green20140901
        {
            childrenElementsIndexes.resize(numberOfChildElements);
            //PHRead(dataContainer, childrenElementsIndexes, numberOfChildElements);
        }
    }
    //*/
}

void ParallelMeshGen::RecoverPoints(Zone * zone)
{
    PHGeometry * geometry = zone->GetGeometry();
    UnstructGrid * unstructuredGrid = UnstructGridCast(geometry->GetGrid());
    
    RDouble * x = unstructuredGrid->GetX();
    RDouble * y = unstructuredGrid->GetY();
    RDouble * z = unstructuredGrid->GetZ();

    Geo_PointFactory * pointFactory = currentZoneGridSolver->GetPointFactory();
    
    RDouble minDistance, maxDistance;
    unstructuredGrid->GetMinMaxDS(minDistance, maxDistance);
    //WriteLogFile("tolerance = ", minDistance/20);
    pointFactory->SetMinDistance(minDistance);


    // 整个单元点数
    int numberOfNodesInElements = unstructuredGrid->GetNTotalNode(); 

    // 点信息、恢复pointArray、pointList 
    for (int iPoint = 0; iPoint < numberOfNodesInElements; ++ iPoint)
    {
        Point3D * point = new Point3D(x[ iPoint ], y[ iPoint ], z[ iPoint ], iPoint);

        bool isNodeExist;
        int pointIndex = pointFactory->AddPoint(point, isNodeExist);
    }

    //读上一次计算单元或网格中的点数
    int oldNumberOfComputationalNodes = numberOfNodesInElements;
    pointFactory->SetOldNumberOfComputationalNodes(oldNumberOfComputationalNodes);
}

void ParallelMeshGen::RecoverFacess(Zone * zone)
{
    int dimension = GetDim();

    PHGeometry * geometry = zone->GetGeometry();
    UnstructGrid * unstructuredGrid = UnstructGridCast(geometry->GetGrid());

    int numberOfFaces         = unstructuredGrid->GetNTotalFace();
    int numberOfBoundaryFaces = unstructuredGrid->GetNBoundFace();
    int numberOfCells         = unstructuredGrid->GetNTotalCell();
    int * leftCellIndexContainer  = unstructuredGrid->GetLeftCellOfFace();
    int * rightCellIndexContainer = unstructuredGrid->GetRightCellOfFace();
    int * faceNodeNumberContainer = unstructuredGrid->GetNodeNumberOfEachFace();
    int * faceNodeIndexContainer  = unstructuredGrid->GetFace2Node();

    ElementFaceSolver * elementFaceSolver = currentZoneGridSolver->GetElementFaceSolver();
    ComputationalFaceSolver * computationalFaceSolver = elementFaceSolver->GetComputationalFaceSolver();
    
    // 上一次计算网格面总数
    PHVectorInt1D & oldComputationalToGlobalFaceIndexMapping = computationalFaceSolver->GetOldComputationalToGlobalFaceIndexMapping();

    int oldNumberOfTotalComputationalElementFaces = numberOfFaces;
    oldComputationalToGlobalFaceIndexMapping.resize(oldNumberOfTotalComputationalElementFaces);

    // 网格面与单元面对应关系
    for (int iFace = 0; iFace < oldNumberOfTotalComputationalElementFaces; ++ iFace)
    {
        oldComputationalToGlobalFaceIndexMapping[ iFace ] = iFace;
    }

    // 左右单元
    PHVectorInt1D & oldLeftElementIndexContainer  = computationalFaceSolver->GetOldLeftElementIndexContainer();
    PHVectorInt1D & oldRightElementIndexContainer = computationalFaceSolver->GetOldRightElementIndexContainer();

    oldLeftElementIndexContainer.resize(oldNumberOfTotalComputationalElementFaces);
    oldRightElementIndexContainer.resize(oldNumberOfTotalComputationalElementFaces);
    for (int iFace = 0; iFace < oldNumberOfTotalComputationalElementFaces; ++ iFace)
    {
        int leftCell  = leftCellIndexContainer[ iFace ];
        int rightCell = rightCellIndexContainer[ iFace ];
        //if(rightCell >= numberOfCells) rightCell = -1;
        if(rightCell >= numberOfCells) rightCell = - rightCell;

        oldLeftElementIndexContainer[ iFace ]  = leftCell;
        oldRightElementIndexContainer[ iFace ] = rightCell;
    }

    
    //读单元简单面总数
    int numberOfTotalFaces = numberOfFaces;

    // set< PHMultiIndex > 简单面集合
    int count = 0;
    for (int iFace = 0; iFace < numberOfTotalFaces; ++ iFace)
    {
        int  faceNodeNumber = faceNodeNumberContainer[ iFace ];
        PHVectorInt1D faceNodeContainer;
        faceNodeContainer.resize(faceNodeNumber);
        for (int iNode = 0; iNode < faceNodeNumber; ++ iNode)
        {
            faceNodeContainer[ iNode ] = faceNodeIndexContainer[ count ++ ];
        }

        sort(faceNodeContainer.begin(), faceNodeContainer.end());

        MultiIndex faceOfMultiIndex(faceNodeNumber, iFace);
        faceOfMultiIndex.SetData(& faceNodeContainer[ 0 ]);
        int id = elementFaceSolver->FindFace(faceOfMultiIndex);
        if(id == PHSPACE::INVALID_INDEX)
        {
            elementFaceSolver->InsertFace(faceOfMultiIndex);
        }
    }


    // 复合面总数
    int numberOfTotalCompositeFaces = 0;

    // set< PHMultiIndex > 复合面集合
    for (int iCompositeFace = 0; iCompositeFace < numberOfTotalCompositeFaces; ++ iCompositeFace)
    {
        TK_Exit::PrintDebugInfoExit("Error: numberOfTotalCompositeFaces != 0 ");
        /*
        int  numberOfNodesOfCompositeFace;
        int  compositeFaceID;

        PHVectorInt1D sortedCompositeFaceNodeIndex(numberOfNodesOfCompositeFace);

        MultiIndex compositeFaceOfMultiIndex(numberOfNodesOfCompositeFace, compositeFaceID);
        compositeFaceOfMultiIndex.SetData(& sortedCompositeFaceNodeIndex[ 0 ]);
        int id = elementFaceSolver->FindCompositeFace(compositeFaceOfMultiIndex);
        if (id == PHSPACE::INVALID_INDEX)
        {
            elementFaceSolver->InsertCompositeFace(compositeFaceOfMultiIndex);
        }
        else
        {
            StopProgram("compositeFaceOfMultiIndex read error \n");
        }*/
    }

    // 单元中face编号和Nodes编号对应关系
    PHVectorInt2D & faceNodeIndexesContainer = elementFaceSolver->GetFaceNodeIndexesContainer();
    faceNodeIndexesContainer.resize(numberOfTotalFaces);
    count = 0;
    for (int iFace = 0; iFace < numberOfTotalFaces; ++ iFace)
    {
        int  faceNodeNumber = faceNodeNumberContainer[ iFace ];
        faceNodeIndexesContainer[ iFace ].resize(faceNodeNumber);

        for (int iNode = 0; iNode < faceNodeNumber; ++ iNode)
        {
            faceNodeIndexesContainer[ iFace ][ iNode ] = faceNodeIndexContainer[ count ++ ];
        }
    }


    // 上一次计算单元或网格中的face编号和Nodes编号对应关系
    PHVectorInt2D & faceNodeIndexesContainerForComputationOld = computationalFaceSolver->GetOldFaceNodeIndexesContainer(); 
    faceNodeIndexesContainerForComputationOld.resize(oldNumberOfTotalComputationalElementFaces);
    count = 0;
    for (int iFace = 0; iFace < numberOfTotalFaces; ++ iFace)
    {
        int  faceNodeNumber = faceNodeNumberContainer[ iFace ];
        faceNodeIndexesContainerForComputationOld[ iFace ].resize(faceNodeNumber);

        for (int iNode = 0; iNode < faceNodeNumber; ++ iNode)
        {
            faceNodeIndexesContainerForComputationOld[ iFace ][ iNode ] = faceNodeIndexContainer[ count ++ ];
        }
    }

    // 单元中compositeFace编号和Nodes编号对应关系
    PHVectorInt2D & nodeIndexesOfCompositeFaceContainer = elementFaceSolver->GetCompositeFaceNodeIndexesContainer();
    nodeIndexesOfCompositeFaceContainer.resize(numberOfTotalCompositeFaces);
    for (int iCompositeFace = 0; iCompositeFace < numberOfTotalCompositeFaces; ++ iCompositeFace)
    {
        TK_Exit::PrintDebugInfoExit("Error: numberOfTotalCompositeFaces == 0 here !!!");
    }

    //读单元中face与左右单元对应关系
    PHVectorInt1D & leftElementIndexContainer = elementFaceSolver->GetLeftElementIndexContainer();
    PHVectorInt1D & rightElementIndexContainer = elementFaceSolver->GetRightElementIndexContainer();
    leftElementIndexContainer.resize(numberOfTotalFaces);
    rightElementIndexContainer.resize(numberOfTotalFaces);
    for (int iFace = 0; iFace < numberOfTotalFaces; ++ iFace)
    {
        int leftCell  = leftCellIndexContainer[ iFace ];
        int rightCell = rightCellIndexContainer[ iFace ];
        //if(rightCell >= numberOfCells) rightCell = -1;
        if(rightCell >= numberOfCells) rightCell = - rightCell;

        leftElementIndexContainer[ iFace ]  = leftCell;
        rightElementIndexContainer[ iFace ] = rightCell;
    }

    CompositeUnstructBC ** bcr = unstructuredGrid->GetBCRecord();

    //读单元中face对应边界条件类型
    PHVectorInt1D & faceBoundaryTypeContainer   = elementFaceSolver->GetFaceBoundaryTypeContainer();
    vector<SimpleBC * > & faceBoundaryCondition = elementFaceSolver->GetFaceBoundaryCondition();
    faceBoundaryTypeContainer.resize(numberOfTotalFaces);
    faceBoundaryCondition.resize(numberOfTotalFaces);
    SetField(faceBoundaryTypeContainer, PHENGLEI::NO_BOUNDARY_CONDITION, numberOfTotalFaces);
    for (int iFace = 0; iFace < numberOfTotalFaces; ++ iFace)
    {
        faceBoundaryCondition[iFace] = 0;
    }

    WriteLogFile("org grid, numberOfBoundaryFaces = ", numberOfBoundaryFaces);
    SimpleBC * interfaceBoundaryCondtion = new SimpleBC();
    Data_Param * bcParamDataBase = new Data_Param();
    interfaceBoundaryCondtion->SetBCName("Interface");
    interfaceBoundaryCondtion->SetBCParamDataBase(bcParamDataBase);
    interfaceBoundaryCondtion->SetBCType(PHENGLEI::INTERFACE);
    for (int iFace = 0; iFace < numberOfBoundaryFaces; ++ iFace)
    {
        int boundaryType = bcr[iFace]->GetKey();
        faceBoundaryTypeContainer[ iFace ] = boundaryType;
        if(boundaryType != PHENGLEI::INTERFACE)
        {
            faceBoundaryCondition[iFace] = bcr[iFace]->GetBoundaryCondition();
        }
        else
        {
            faceBoundaryCondition[iFace] = interfaceBoundaryCondtion;
        }
    }

    // 单元中face对应单元类型
    PHVectorInt1D & faceTypeContainer = elementFaceSolver->GetFaceTypeContainer();
    faceTypeContainer.resize(numberOfTotalFaces);
    for (int iFace = 0; iFace < numberOfTotalFaces; ++ iFace)
    {
        int  faceNodeNumber = faceNodeNumberContainer[ iFace ];
        int  faceType       = GetFaceTypeOfFTS(dimension, faceNodeNumber);
        faceTypeContainer[ iFace ] = faceType;
    }

    // 读parentFaceIndexContainer
    PHVectorInt1D & parentFaceIndexContainer = elementFaceSolver->GetParentFaceIndexContainer();
    parentFaceIndexContainer.resize(numberOfTotalFaces);
    SetField(parentFaceIndexContainer, PHSPACE::INVALID_INDEX, numberOfTotalFaces);

    // 读compositeFaceTypeContainer
    PHVectorInt1D & compositeFaceTypeContainer = elementFaceSolver->GetCompositeFaceTypeContainer();
    compositeFaceTypeContainer.resize(numberOfTotalCompositeFaces);

    // 读childFaceIndexContainer
    PHVectorInt2D & childFaceIndexContainer = elementFaceSolver->GetChildFaceIndexContainer();
    childFaceIndexContainer.resize(numberOfTotalFaces);


    //读简单面与复合面对应关系
    PHVectorInt1D & compositeToSimpleFaceIndexMapping = elementFaceSolver->GetCompositeToSimpleFaceIndexMapping();
    PHVectorInt1D & simpleToCompositeFaceIndexMapping = elementFaceSolver->GetSimpleToCompositeFaceIndexMapping();

    compositeToSimpleFaceIndexMapping.resize(numberOfTotalCompositeFaces);
    simpleToCompositeFaceIndexMapping.resize(numberOfTotalFaces);

    if (numberOfTotalCompositeFaces != 0)
    {
        TK_Exit::PrintDebugInfoExit("Error: numberOfTotalCompositeFaces != 0");
    }

    for (int iCompositeFace = 0; iCompositeFace < numberOfTotalCompositeFaces; ++ iCompositeFace)
    {
        TK_Exit::PrintDebugInfoExit("Error: Consider here!");
        int simpleFaceIndex = compositeToSimpleFaceIndexMapping[ iCompositeFace ];
        simpleToCompositeFaceIndexMapping[ simpleFaceIndex ] = iCompositeFace;
    }
}

void ParallelMeshGen::RecoverMetrics(Zone * zone)
{
    PHGeometry * geometry = zone->GetGeometry();
    UnstructGrid * unstructuredGrid = UnstructGridCast(geometry->GetGrid());
    
    unstructuredGrid->ComputeMetrics();
    
    RDouble * faceAreaContainer    = unstructuredGrid->GetFaceArea();
    RDouble * faceNormalXContainer = unstructuredGrid->GetFaceNormalX();
    RDouble * faceNormalYContainer = unstructuredGrid->GetFaceNormalY();
    RDouble * faceNormalZContainer = unstructuredGrid->GetFaceNormalZ();

    ElementFaceSolver * elementFaceSolver = currentZoneGridSolver->GetElementFaceSolver();
    PHVectorRDouble1D & faceArea = elementFaceSolver->GetFaceArea();
    PHVectorRDouble1D & faceNormalX = elementFaceSolver->GetFaceNormalX();
    PHVectorRDouble1D & faceNormalY = elementFaceSolver->GetFaceNormalY();
    PHVectorRDouble1D & faceNormalZ = elementFaceSolver->GetFaceNormalZ();

    int numberOfFaces         = unstructuredGrid->GetNTotalFace();
    faceArea.resize(numberOfFaces);
    faceNormalX.resize(numberOfFaces);
    faceNormalY.resize(numberOfFaces);
    faceNormalZ.resize(numberOfFaces);

    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        faceArea[ iFace ] = faceAreaContainer[ iFace ];
        faceNormalX[ iFace ] = faceNormalXContainer[ iFace ];
        faceNormalY[ iFace ] = faceNormalYContainer[ iFace ];
        faceNormalZ[ iFace ] = faceNormalZContainer[ iFace ];
    }
}

void ParallelMeshGen::RecoverInterFaceInformation(Zone * zone)
{
    PHGeometry * geometry = zone->GetGeometry();
    UnstructGrid * unstructuredGrid = UnstructGridCast(geometry->GetGrid());

    currentZoneGridSolver->SetInterfaceInformation(unstructuredGrid->GetInterfaceInfo());
}

void ParallelMeshGen::BroadcastNeighborZonesToEachProcessor()
{
    WriteLogFile("Start Broadcasting Neighbor Zones To Each Processor!");

    int nZones = PHMPI::GetNumberofGlobalZones();

    // 并行标准模式：遍历全局所有的zone
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        PHVector1D < int > neighborZones;

        if(ZoneBelongToCurrentProcessor(iZone))
        {
            int localZoneIndex = region->GetProcessGlobalZoneIndexToLocalZoneIndex(iZone);

            ZoneGridSolver * zoneGridSolver = zoneGridSolverContainer[ localZoneIndex ];
            this->SetCurrentZoneGridSolver(zoneGridSolver);
            //this->SetCurrentAdaptiveFlowFieldSolver(adaptiveFlowFieldSolverVector[ localZoneIndex ]);
            neighborZones = GetNeighborZones(zoneGridSolver);
        }
        
        this->BroadcastToEachProcessor(iZone, neighborZones);
    }

    //PH_Barrier();

    WriteLogFile("End Broadcasting Neighbor Zones To Each Processor!");


    /*
    PHVector2D < int > & globalNeighborZones = this->GetCurrentNeighborZones(0);
    for (int i = 0; i < globalNeighborZones.size(); ++ i)
    {
        cout << PHMPI::GetCurrentProcessorID() << ": " << i << " ";
        for (int j = 0; j < globalNeighborZones[ i ].size(); ++ j)
        {
            cout << globalNeighborZones[ i ][ j ] << " ";
        }
        cout << endl;
    }

    exit(0);
    //*/
}

void ParallelMeshGen::BroadcastNeighborZonesToEachProcessorNew()
{
    int tag = 0;
    int nLocalZones = PHMPI::GetNumberofLocalZones();
    for (int localZone = 0; localZone < nLocalZones; ++ localZone )
    {
        ZoneGridSolver * zoneGridSolver = this->GetZoneGridSolver(localZone);
        PHVector1D < int > neighborZones = GetNeighborZones(zoneGridSolver);
        int numberOfNeighborZones    = static_cast<int>(neighborZones.size());
        int currentZoneIndex         = zoneGridSolver->GetZoneIndex();
        int currentZoneglobalIndex   = region->GetProcessLocalZoneIndexToGlobalZoneIndex(currentZoneIndex);

        this->AddOneZoneNeighborToAllZonesNeighbor(neighborZones, currentZoneglobalIndex, currentZoneglobalIndex);

        // Find the mid-zone.
        int midN, iZone;
        for (iZone = 0; iZone < numberOfNeighborZones; ++ iZone)
        {
            if (neighborZones[ iZone ] >= currentZoneIndex)
            {
                break;
            }
        }
        midN = iZone;

        // Part1: receiving messages from neighbors with zone number less than myZone-1
        for (iZone = 0; iZone < midN; ++ iZone)
        {
            int neighborZoneIndex  = neighborZones[ iZone ];
            int sendProcessorIndex = PHMPI::GetZoneProcessorID(neighborZoneIndex);

            DataContainer * dataContainer = new DataContainer();

            receive(dataContainer, sendProcessorIndex, tag);

            // Decode.
            dataContainer->moveto_begin();
            int size;
            PHRead(dataContainer, &size, 1);
            PHVectorInt1D neighborZoneVec;
            neighborZoneVec.resize(size);
            PHRead(dataContainer, neighborZoneVec, size);

            this->AddOneZoneNeighborToAllZonesNeighbor(neighborZoneVec, neighborZoneIndex, currentZoneglobalIndex);

            delete dataContainer;
        }

        // Part2: sending messages to neighbors with zone number big than myZone-1
        if(midN < numberOfNeighborZones)
        {
            for (iZone = midN; iZone < numberOfNeighborZones; ++ iZone)
            {
                int neighborZoneIndex     = neighborZones[ iZone ];
                int receiveProcessorIndex = PHMPI::GetZoneProcessorID(neighborZoneIndex);

                // Encode.
                DataContainer * dataContainer = new DataContainer();

                int size = static_cast<int>(neighborZones.size());
                PHWrite(dataContainer, &size, 1);
                PHWrite(dataContainer, neighborZones, size);

                send(dataContainer, receiveProcessorIndex, tag);

                delete dataContainer;
            }
        }

        // Part3: receiving messages from neighbors with zone number big than myZone-1
        for (iZone = midN; iZone < numberOfNeighborZones; ++ iZone)
        {
            int neighborZoneIndex  = neighborZones[ iZone ];
            int sendProcessorIndex = PHMPI::GetZoneProcessorID(neighborZoneIndex);

            DataContainer * dataContainer = new DataContainer();

            receive(dataContainer, sendProcessorIndex, tag);

            // Decode.
            dataContainer->moveto_begin();
            int size;
            PHRead(dataContainer, &size, 1);
            PHVectorInt1D neighborZoneVec;
            neighborZoneVec.resize(size);
            PHRead(dataContainer, neighborZoneVec, size);

            this->AddOneZoneNeighborToAllZonesNeighbor(neighborZoneVec, neighborZoneIndex, currentZoneglobalIndex);

            delete dataContainer;
        }

        // Part4: sending messages to neighbors with zone number less than myZone-1  
        if (midN > 0)
        {
            for (iZone = 0; iZone < midN; ++ iZone)
            {
                int neighborZoneIndex     = neighborZones[ iZone ];
                int receiveProcessorIndex = PHMPI::GetZoneProcessorID(neighborZoneIndex);

                // Encode.
                DataContainer * dataContainer = new DataContainer();

                int size = static_cast<int>(neighborZones.size());
                PHWrite(dataContainer, &size, 1);
                PHWrite(dataContainer, neighborZones, size);

                send(dataContainer, receiveProcessorIndex, tag);

                delete dataContainer;
            }
        }
    }
}

void ParallelMeshGen::BroadcastToEachProcessor(const int zoneIndex, PHVector1D < int > & dataOfThisProcessor)
{
    PHVector1D < int > targetZones = GetAllOtherZones(zoneIndex);

    int numberOfTargetZones = static_cast<int>(targetZones.size());
    for (int iTarget = 0; iTarget < numberOfTargetZones; ++ iTarget)
    {
        int targetZoneIndex       = targetZones[ iTarget ];
        int sendProcessorIndex    = PHMPI::GetZoneProcessorID(zoneIndex);
        int receiveProcessorIndex = PHMPI::GetZoneProcessorID(targetZoneIndex);

        DataContainer * dataContainer = new DataContainer();

        if (PHMPI::GetCurrentProcessorID() == sendProcessorIndex)
        {
            CompressVectorIntoDataContainer(dataContainer, dataOfThisProcessor);
        }

        PHMPI::PH_Trade(dataContainer, sendProcessorIndex, receiveProcessorIndex);

        if (PHMPI::GetCurrentProcessorID() == receiveProcessorIndex)
        {
            PHVector1D < int > receivedData;
            AnasysVectorFromDataContainer(dataContainer, receivedData);

            this->AddOneZoneNeighborToAllZonesNeighbor(receivedData, zoneIndex, targetZoneIndex);
        }

        delete dataContainer;
    }
}

void ParallelMeshGen::InitialInterFaceInformation()
{
    int nZones = PHMPI::GetNumberofGlobalZones();

    ZoneConnectivity * zoneNeighborConnectivity = new ZoneConnectivity();
    zoneNeighborConnectivity->CreateAllZones(nZones);

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        if (!GetZone(iZone)) continue;

        int localZoneIndex = region->GetProcessGlobalZoneIndexToLocalZoneIndex(iZone);
        ZoneGridSolver * zoneGridSolver = this->GetZoneGridSolver(localZoneIndex);

        InterfaceInfo * interfaceInformation = zoneGridSolver->GetInterfaceInformation();
        if (!interfaceInformation) continue;

        interfaceInformation->InitNeighbors();
    }

    delete zoneNeighborConnectivity;
}

void ParallelMeshGen::SwapNeighborsSendContent()
{
    WriteLogFile("Start Swapping Neighbors Send Content!");

    int nZones = PHMPI::GetNumberofGlobalZones();
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    int tag = 0;

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int sendProcessorIndex = PHMPI::GetZoneProcessorID(iZone);

        if(GetZone(iZone))
        {
            int localZoneIndex = region->GetProcessGlobalZoneIndexToLocalZoneIndex(iZone);
            ZoneGridSolver * zoneGridSolver = this->GetZoneGridSolver(localZoneIndex);
            this->SetCurrentZoneGridSolver(zoneGridSolver);
            //this->SetCurrentAdaptiveFlowFieldSolver(adaptiveFlowFieldSolverVector[ localZoneIndex ]);
        }

        PHVector2D < int > & globalNeighborZones = this->GetCurrentNeighborZones(iZone);
        
        int numberOfNeighbors = static_cast<int>(globalNeighborZones[ iZone ].size());
        for (int neighborIndex = 0; neighborIndex < numberOfNeighbors; ++ neighborIndex)
        {
            int neighborZoneIndex = globalNeighborZones[ iZone ][ neighborIndex ];
            int receiveProcessorIndex = PHMPI::GetZoneProcessorID(neighborZoneIndex);
            int interfaceNumberBetweenTwoNeighboringZone = 0;
            int * interfaceIndexContainerForSend = 0;

            if (currentProcessorIndex == sendProcessorIndex)
            {
                InterfaceInfo * interfaceInformation = currentZoneGridSolver->GetInterfaceInformation();
                interfaceNumberBetweenTwoNeighboringZone = interfaceInformation->ComputeNIFaceOfNeighbor(neighborIndex);
                interfaceIndexContainerForSend           = interfaceInformation->ComputeFaceIndexForSend(neighborIndex);
            }
            
            PHMPI::PH_Trade(& interfaceNumberBetweenTwoNeighboringZone, sizeof(int), sendProcessorIndex, receiveProcessorIndex, tag);

            
            if (currentProcessorIndex == receiveProcessorIndex && 
                sendProcessorIndex    != receiveProcessorIndex)
            {
                interfaceIndexContainerForSend = new int [ interfaceNumberBetweenTwoNeighboringZone ];
            }

            PHMPI::PH_Trade(interfaceIndexContainerForSend, interfaceNumberBetweenTwoNeighboringZone * sizeof(int), sendProcessorIndex, receiveProcessorIndex, tag);

            if (currentProcessorIndex == receiveProcessorIndex)
            {
                int neighborLocalZoneIndex = region->GetProcessGlobalZoneIndexToLocalZoneIndex(neighborZoneIndex);
                InterfaceInfo * neighborZoneInterfaceInformation = this->GetZoneGridSolver(neighborLocalZoneIndex)->GetInterfaceInformation();
                neighborZoneInterfaceInformation->FillFaceIndexForSend(iZone, interfaceIndexContainerForSend);
            }

            delete [] interfaceIndexContainerForSend;    interfaceIndexContainerForSend = NULL;
        }
    }

    //PH_Barrier();
}

void ParallelMeshGen::SwapInterfaceRefineType()
{
    PHSPACE::WriteLogFile("Start Swapping refine type!");

    CompressData compressData     = CompressFixCellType;
    DeCompressData decompressData = DecompressFixCellType;

    SwapInterfaceData(this->region, this->zoneGridSolverContainer, compressData, decompressData);

    PHSPACE::WriteLogFile("End Swapping refine type!");
}

void ParallelMeshGen::CompressVector(DataContainer *& dataContainer, int globalZoneIndex, int neighborZoneIndex)
{
    int localZoneIndex = region->GetProcessGlobalZoneIndexToLocalZoneIndex(globalZoneIndex);
    ZoneGridSolver * zoneGridSolver = this->GetZoneGridSolver(localZoneIndex);
    zoneGridSolver->CompressRefineType(dataContainer, neighborZoneIndex);
}

void ParallelMeshGen::DecompressVector(DataContainer *& dataContainer, int globalZoneIndex, int neighborZoneIndex)
{
    // 注意这里的neighbor关系：
    // 此时是neighborZone解析数据，因此是站在neighborZone的角度，而globalZone反而成了他的邻居.
    int localZoneIndex = region->GetProcessGlobalZoneIndexToLocalZoneIndex(neighborZoneIndex);
    ZoneGridSolver * zoneGridSolver = this->GetZoneGridSolver(localZoneIndex);
    zoneGridSolver->DecompressRefineType(dataContainer, globalZoneIndex);
}

void ParallelMeshGen::AddOneZoneNeighborToAllZonesNeighbor(PHVector1D < int > & dataOfThisZone, int zoneIndex, int targetZoneIndex)
{
    if (dataOfThisZone.size() == 0)
    {
        return;
    }

    int nZones = PHMPI::GetNumberofGlobalZones();

    int localZoneIndexOfTargetZone = region->GetProcessGlobalZoneIndexToLocalZoneIndex(targetZoneIndex);
    ZoneGridSolver * zoneGridSolver = zoneGridSolverContainer[ localZoneIndexOfTargetZone ];
    PHVector2D < int > & globalNeighborZones       = zoneGridSolver->GetGlobalNeighborZones();
    if(globalNeighborZones.size() == 0) globalNeighborZones.resize(nZones);

    globalNeighborZones[ zoneIndex ].resize(dataOfThisZone.size());
    int dataOfThisZoneSize = static_cast<int>(dataOfThisZone.size());
    for (int iNeighbor = 0; iNeighbor < dataOfThisZoneSize; ++ iNeighbor)
    {
        globalNeighborZones[ zoneIndex ][ iNeighbor ] = dataOfThisZone[ iNeighbor ];
    }
}

PHVector1D < int > ParallelMeshGen::GetAllOtherZones(int zoneIndex)
{
    PHVector1D < int > allOtherZones;

    int nZones = PHMPI::GetNumberofGlobalZones();
    allOtherZones.resize(nZones);
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        allOtherZones[ iZone ] = iZone;
    }

    return allOtherZones;
}

int ParallelMeshGen::GetElementTypeOfFTS(const int & dimension, const int & nodeNumberOfCell)
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

int ParallelMeshGen::GetFaceTypeOfFTS(const int & dimension, const int & nodeNumberOfCell)
{
    int faceType = -1;
    if (dimension == TWO_D)
    {
        switch (nodeNumberOfCell)
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
        switch (nodeNumberOfCell)
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

Zone * ParallelMeshGen::GetZone(int izone)
{
    return region->GetZone(izone);
}

bool ParallelMeshGen::ZoneBelongToCurrentProcessor(int zoneIndex)
{
    if (GetZone(zoneIndex)) 
    {
        return true;
    }
    else
    {
        return false;
    }
}

PHVector1D < int > ParallelMeshGen::GetNeighborZones(ZoneGridSolver * zoneGridSolver)
{
    InterfaceInfo * interfaceInformation = zoneGridSolver->GetInterfaceInformation();
    int * interFace2ZoneID = interfaceInformation->GetInterFace2ZoneID();
    int numberOfInterfaces = interfaceInformation->GetNIFace();
    set< int > neighborZoneList;

    // Get neighbor zones information.
    for (int iFace = 0; iFace < numberOfInterfaces; ++ iFace)
    {
        int neighborZoneIndex = interFace2ZoneID[iFace];
        neighborZoneList.insert(neighborZoneIndex);
    }
    PHVector1D < int > neighborZones;
    for (set< int >::iterator iter = neighborZoneList.begin(); iter != neighborZoneList.end(); ++ iter)
    {
        neighborZones.push_back(* iter);
    }
    sort(neighborZones.begin(), neighborZones.end());

    return neighborZones;
}

PHVector2D < int > & ParallelMeshGen::GetCurrentNeighborZones(int zoneIndex)
{
    if(ZoneBelongToCurrentProcessor(zoneIndex))
    {
        int localZoneIndex = region->GetProcessGlobalZoneIndexToLocalZoneIndex(zoneIndex);
        ZoneGridSolver * zoneGridSolver = this->GetZoneGridSolver(localZoneIndex);
        return zoneGridSolver->GetGlobalNeighborZones();
    }

    int nZones = PHMPI::GetNumberofGlobalZones();
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        if(ZoneBelongToCurrentProcessor(iZone))
        {
            int localZoneIndex = region->GetProcessGlobalZoneIndexToLocalZoneIndex(iZone);
            ZoneGridSolver * zoneGridSolver = zoneGridSolverContainer[ localZoneIndex ];

            return zoneGridSolver->GetGlobalNeighborZones();
        }
    }

    // 这是权宜之计，因为总要返回一个东西。
    return this->ZERO_VECTOR2D;
}


void ParallelMeshGen::GenerateAndDumpComputationalGrid()
{
    WriteLogFile("Start Generating Grid!");

    AllocateComputationalGrid();

    GenerateComputationalGrid();

    DumpComputationalGrid();

    MeshQualityChecking();

    //DeAllocateComputationalGrid();
}

void ParallelMeshGen::AllocateComputationalGrid()
{
    int nLocalZones = PHMPI::GetNumberofLocalZones();
    fantasyGrids  = new Grid * [ nLocalZones ];
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        Grid *grid  = CreateBasicGrid(iZone);
        fantasyGrids[ iZone ] = grid;
    }
}

void ParallelMeshGen::GenerateComputationalGrid()
{
    int nLocalZones = PHMPI::GetNumberofLocalZones();
    int nTCellLocal = 0;
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        zoneGridSolverContainer[ iZone ]->GenerateComputationalGrid(fantasyGrids[ iZone ], GetGridType(iZone));
        
        nTCellLocal += zoneGridSolverContainer[iZone]->GetElementFaceSolver()->GetComputationalFaceSolver()->GetNumberOfElements();
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

void ParallelMeshGen::MeshQualityChecking()
{
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        Grid *grid = fantasyGrids[iZone];

        UnstructGrid * unstructuredGrid = UnstructGridCast(grid);

        unstructuredGrid->ComputeMetrics(0);

        unstructuredGrid->SkewnessSummary();
    }
}

void ParallelMeshGen::DumpComputationalGrid()
{
    string adapt_grid_file;
    GlobalDataBase::GetData("adapt_grid_file", &adapt_grid_file, PHSTRING, 1);

    int nLocalZones = PHMPI::GetNumberofLocalZones();

	ComputationalGridProcess* computationalGridProcess = new ComputationalGridProcess(region);
	computationalGridProcess->Initialize(fantasyGrids, nLocalZones);
	computationalGridProcess->PostProcessMultiZoneComputationalGrids();
	delete computationalGridProcess;

    DumpMultiZoneComputationalGrids(adapt_grid_file, fantasyGrids, nLocalZones);
}

void ParallelMeshGen::DumpMultiZoneComputationalGrids(const string & fileName, Grid ** gridContainer, int nZones)
{
    WriteLogFile(" Grid Writing ...");
    if(PHMPI::GetCurrentProcessorID() == 0) cout << "  Grid Writing ... " << endl;

    if(fileName =="") return;

    int dumpOldGrid = GlobalDataBase::GetIntParaFromDB("dumpOldGrid");
    if(!dumpOldGrid)
    {
        //! Dump new HDF5 version grid.
        PrintToWindow("    Writing new HDF5 version grid to file ....\n");
        WriteLogFile("    Writing  new HDF5 version grid to file ....\n");

        IO_HDF5Write writeGrid(fileName, gridContainer, nZones);
        writeGrid.Run();
    }
    else
    {
        //! Dump old version grid.
        PrintToWindow("    Writing old version grid to file ....\n");
        WriteLogFile("    Writing old version grid to file ....\n");
        int fileMode = P2P_MODE;
        if (fileMode == PHSPACE::P2P_MODE)
        {
            DumpMultiZoneComputationalGridsP2PMode(fileName, gridContainer, nZones);
        }
        else if (fileMode == PHSPACE::CS_MODE)
        {
            DumpMultiZoneComputationalGridsCSMode(fileName, gridContainer, nZones);
        }

        WriteFaceBC(fileName, gridContainer, nZones);

        WriteCellToNode(fileName);
    }

    WriteLogFile(" End Writing Grid.");
}

void ParallelMeshGen::DumpMultiZoneComputationalGridsP2PMode(const string &fileName, Grid **gridContainer, int nZones)
{
    using namespace PHMPI;

    //注意，dumpgrid直接写入网格文件，会影响网格块数，不可随意输出
    int numberOfSimulationZones = PHMPI::GetNumberofGlobalZones();

    int * zoneProcessorIndexContainer         = GetZoneProcessorID();
    int * zoneIndexContainer                  = GetZoneGridID();
    int * zoneTypeContainer                   = GetZoneGridType();
    int * zoneProcessorIndexContainerFromFile = GetZoneProcessorID_INP();

    if (!zoneProcessorIndexContainer)
    {
        CreateZoneProcessorID(numberOfSimulationZones);
        CreateZoneGridID(numberOfSimulationZones);
        CreateZoneGridType(numberOfSimulationZones);

        int * zoneProcessorIndexContainer         = GetZoneProcessorID();
        int * zoneIndexContainer                  = GetZoneGridID();
        int * zoneTypeContainer                   = GetZoneGridType();

        for (int iZone = 0; iZone < numberOfSimulationZones; ++ iZone)
        {
            zoneProcessorIndexContainer[ iZone ] = iZone;
            zoneIndexContainer         [ iZone ] = iZone;
            zoneTypeContainer          [ iZone ] = gridContainer[ iZone ]->Type();
        }
    }  

    if (zoneProcessorIndexContainerFromFile)
    {
        for (int iZone = 0; iZone < nZones; ++ iZone)
        {
            SetZoneProcessorID(iZone, zoneProcessorIndexContainerFromFile[ iZone ]);
        }
    }

    if (GetNumberOfProcessor() < numberOfSimulationZones)
    {
        zoneProcessorIndexContainer = new int [ numberOfSimulationZones ]();
        // Reset the zone processor.
        for (int iZone = 0; iZone < numberOfSimulationZones; ++ iZone)
        {
            zoneProcessorIndexContainer[ iZone ] = iZone;
        }
    }  

    // Self Mode.
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int globalZoneIndex = region->GetProcessLocalZoneIndexToGlobalZoneIndex(iZone);

        fstream file;
        string processorGridFileName = AddSymbolToFileName(fileName, '_', globalZoneIndex);
        processorGridFileName = ChangeExtensionOfFileName(processorGridFileName, "fts");
        PHSPACE::OpenFile(file, processorGridFileName, ios_base::out|ios_base::binary|ios_base::trunc);

        int localProcessor[1];
        int localZoneIndex[1]; 
        int localZoneType[1];

        int nZonesLocal = 1;
        localProcessor[0] = zoneProcessorIndexContainer[globalZoneIndex];
        localZoneIndex[0] = zoneIndexContainer[globalZoneIndex]         ;
        localZoneType[0]  = zoneTypeContainer[globalZoneIndex]          ;
        PHWrite(file, nZonesLocal);
        PHWrite(file, localProcessor, nZonesLocal);
        PHWrite(file, localZoneIndex, nZonesLocal);
        PHWrite(file, localZoneType , nZonesLocal);

        gridContainer[ iZone ]->WriteGrid(file);

        PHSPACE::CloseFile(file);
    }
    if (GetNumberOfProcessor() < numberOfSimulationZones)
    {
        delete [] zoneProcessorIndexContainer;    zoneProcessorIndexContainer = NULL;
    }
}

void ParallelMeshGen::DumpMultiZoneComputationalGridsCSMode(const string &fileName, Grid **gridContainer, int nZones)
{
    using namespace PHMPI;

    //注意，dumpgrid直接写入网格文件，会影响网格块数，不可随意输出
    int numberOfSimulationZones = PHMPI::GetNumberofGlobalZones();

    CreateZoneProcessorID(numberOfSimulationZones);
    CreateZoneGridID(numberOfSimulationZones);
    CreateZoneGridType(numberOfSimulationZones);

    int * zoneProcessorIndexContainer         = GetZoneProcessorID();
    int * zoneIndexContainer                  = GetZoneGridID();
    int * zoneTypeContainer                   = GetZoneGridType();
    int * zoneProcessorIndexContainerFromFile = GetZoneProcessorID_INP();

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        zoneProcessorIndexContainer[ iZone ] = iZone;
        zoneIndexContainer         [ iZone ] = iZone;
        zoneTypeContainer          [ iZone ] = gridContainer[ iZone ]->Type();
    }

    if (zoneProcessorIndexContainerFromFile)
    {
        for (int iZone = 0; iZone < nZones; ++ iZone)
        {
            SetZoneProcessorID(iZone, zoneProcessorIndexContainerFromFile[ iZone ]);
        }
    }

    fstream file;
    PHSPACE::OpenFile(file, fileName, ios_base::out|ios_base::binary|ios_base::trunc);

    PHWrite(file, nZones);
    PHWrite(file, zoneProcessorIndexContainer, nZones);
    PHWrite(file, zoneIndexContainer         , nZones);
    PHWrite(file, zoneTypeContainer          , nZones);

    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        gridContainer[ iZone ]->WriteGrid(file);
    }

    PHSPACE::CloseFile(file);

    WriteLogFile(" End Writing Grid.");
}

void ParallelMeshGen::WriteFaceBC(const string &gridfileName, Grid ** gridContainer, int nZones)
{
    int numberOfSimulationZones = PHMPI::GetNumberofGlobalZones();

    int fileMode = P2P_MODE;
    if (fileMode == PHSPACE::CS_MODE)
    {
        // serial model.
        string faceBCFileName = ChangeExtensionOfFileName(gridfileName, "bc");
        fstream file;
        PHSPACE::OpenFile(file, faceBCFileName, ios_base::out|ios_base::binary);

        PHWrite(file, numberOfSimulationZones);
        for (int iZone = 0; iZone < nZones; ++ iZone)
        {
            UnstructGrid * grid = UnstructGridCast(gridContainer[iZone]);
            grid->WriteFaceBC(file);
        }

        PHSPACE::CloseFile(file);
    }
    else if (fileMode == PHSPACE::P2P_MODE)
    {
        //! parallel model.
        for (int iZone = 0; iZone < nZones; ++ iZone)
        {
            int globalZoneIndex = region->GetProcessLocalZoneIndexToGlobalZoneIndex(iZone);

            fstream file;
            string processorGridFileName = AddSymbolToFileName(gridfileName, '_', globalZoneIndex);
            string faceBCFileName = ChangeExtensionOfFileName(processorGridFileName, "bc");
            PHSPACE::OpenFile(file, faceBCFileName, ios_base::out|ios_base::binary);

            PHWrite(file, numberOfSimulationZones);

            UnstructGrid * grid = UnstructGridCast(gridContainer[iZone]);
            grid->WriteFaceBC(file);

            PHSPACE::CloseFile(file);
        }
    }
}

void ParallelMeshGen::WriteCellToNode(const string &gridfileName)
{
    int fileMode = P2P_MODE;

    if(gridfileName =="") return;

    if (fileMode == PHSPACE::P2P_MODE)
    {
        int nZones = PHMPI::GetNumberofGlobalZones();

        for (int iZone = 0; iZone < numberOfZones; ++ iZone)
        {
            int globalZoneIndex = region->GetProcessLocalZoneIndexToGlobalZoneIndex(iZone);

            fstream file;
            string tempName = AddSymbolToFileName(gridfileName, '_', globalZoneIndex);
            string cellToNodeFileName = ChangeExtensionOfFileName(tempName, "c2n");
            PHSPACE::OpenFile(file, cellToNodeFileName, ios_base::out|ios_base::binary);

            PHWrite(file, nZones);

            VirtualFile * virtualFile = new VirtualFile(&file);
            virtualFile->BeginWriteWork();

            PHWrite(virtualFile, globalZoneIndex);

            ZoneGridSolver * zoneGridSolver = this->GetZoneGridSolver(iZone);

            ElementPropertySolver * elementPropertySolver = zoneGridSolver->GetElementPropertySolver();
            elementPropertySolver->DumpElementToNode(virtualFile);

            virtualFile->EndWriteWork();
            delete virtualFile;

            PHSPACE::CloseFile(file);
        }
    }
    else if (fileMode == PHSPACE::CS_MODE)
    {
        string cellToNodeFileName = ChangeExtensionOfFileName(gridfileName, "c2n");
        fstream file;
        PHSPACE::OpenFile(file, cellToNodeFileName, ios_base::out|ios_base::binary);

        PHWrite(file, numberOfZones);

        for (int iZone = 0; iZone < numberOfZones; ++ iZone)
        {
            VirtualFile * virtualFile = new VirtualFile(&file);
            virtualFile->BeginWriteWork();

            PHWrite(virtualFile, iZone);

            ZoneGridSolver * zoneGridSolver = this->GetZoneGridSolver(iZone);

            ElementPropertySolver * elementPropertySolver = zoneGridSolver->GetElementPropertySolver();
            elementPropertySolver->DumpElementToNode(virtualFile);

            virtualFile->EndWriteWork();
            delete virtualFile;
        }

        PHSPACE::CloseFile(file);
    }    
}

void ParallelMeshGen::ConstructGrid()
{
    ReadGrid();

    //从读入的标准网格信息（CGNS）解析出实际计算需要的网格信息
    ResolveComputationGrid();

    CreateGeneralMultiBlockGrid();
}

void ParallelMeshGen::ResolveComputationGrid()
{
    for (int izone = 0; izone < numberOfZones; ++ izone)
    {
        super_elements[izone]->ResolveComputationGrid();
    }
}

int ParallelMeshGen::GetGridType(int izone) const
{
    //return cgns_multizone->GetZone(izone)->GetZoneType();
    return (int)Unstructured;
}

Grid * ParallelMeshGen::CreateBasicGrid(int izone)
{
    int grid_type = GetGridType(izone);
    if (grid_type == Unstructured)
    {
        return CreateGridGeneral(UNSTRUCTGRID, new GridID(izone), 0, GetDim());
    }
    else
    {
        return CreateGridGeneral(STRUCTGRID, new GridID(izone), 0, GetDim());
    }
}

void ParallelMeshGen::CreateAllSingleGrid(Grid ** grids)
{
    for (int izone = 0; izone < numberOfZones; ++ izone)
    {
        Grid *grid  = CreateBasicGrid(izone);
        grids[izone] = grid;
        super_elements[izone]->CreateFantasyGrid(grid, GetGridType(izone));
    }
}

void ParallelMeshGen::CreateGeneralMultiBlockGrid()
{
    Grid ** grids  = new Grid * [numberOfZones];

    string out_gfile = "grid.fts";
    GlobalDataBase::GetData("out_gfile", &out_gfile, PHSTRING, 1);

    CreateAllSingleGrid(grids);

    ConstructGlobalInterfaceLink(grids, numberOfZones);

    DumpGrid(out_gfile, grids, numberOfZones);

    for (int izone = 0; izone < numberOfZones; ++ izone)
    {
        delete grids[izone];
    }
    delete [] grids;
}

void ParallelMeshGen::ReadGrid()
{
    string from_gfile = "grid.fts";
    GlobalDataBase::GetData("from_gfile", &from_gfile, PHSTRING, 1);

    int iadapt = 0;
    GlobalDataBase::GetData("iadapt", &iadapt, PHINT, 1);

    cgns_multizone = new CGNSMultiZone();
    cgns_multizone->ReadGrid(from_gfile);
    numberOfZones = cgns_multizone->GetNZones();

    super_elements = new SuperElement * [numberOfZones];
    for (int izone = 0; izone < numberOfZones; ++ izone)
    {
        super_elements[izone] = new SuperElement(cgns_multizone->GetZone(izone), iadapt);
    }
}

bool ParallelMeshGen::IfNeedProjection()
{
    bool ifNeed = false;

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
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

void ParallelMeshGen::SetAdaptTropic()
{
//  int anisoRefine = GlobalDataBase::GetIntParaFromDB("anisoRefine");
// 
//  int numberOfBasicElement = PHSPACE::GetNumberOfBasicElement();
// 
//  for (int iBasicElement = 0; iBasicElement < numberOfBasicElement; ++ iBasicElement)
//  {
//      GetBasicElement(iBasicElement)->SetAnisotropicAdapt(anisoRefine);
//  }
}

void ParallelMeshGen::InitCellRefineType()
{
    int nLocalZones = PHMPI::GetNumberofLocalZones();
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        ZoneGridSolver * zoneGridSolver = zoneGridSolverContainer[ iZone ];

        InitCellRefineType(zoneGridSolver);
    }
}

void ParallelMeshGen::InitCellRefineType(ZoneGridSolver * zoneGridSolver)
{
    ElementPropertySolver * elementPropertySolver = zoneGridSolver->GetElementPropertySolver();
    ElementFaceSolver * elementFaceSolver = zoneGridSolver->GetElementFaceSolver();

    int nTCell = zoneGridSolver->GetNumberOfElements();
    int nBoundFace = elementFaceSolver->GetNumberOfBoundaryFaces();
    int nTotal = nTCell + nBoundFace;

    PHVectorInt1D & elementType = elementPropertySolver->GetElementTypeContainer();
    PHVectorInt1D & refineType = * elementPropertySolver->GetRefineType();
    refineType.resize(nTotal);

    int anisoRefine = GlobalDataBase::GetIntParaFromDB("anisoRefine");

    if(anisoRefine == ISOTROPIC_REFINE)
    {
        // Set all cells to be anisotropic refine.
        SetField(refineType, ISOTROPIC_REFINE, nTotal);
    }
    else
    {
        // Set all cells to be anisotropic refine.
        SetField(refineType, ANISOTROPIC_REFINE, nTotal);

        // Set all the non-prism cells to be isotropic refine.
        for (int iCell = 0; iCell < nTCell; ++ iCell)
        {
            if(elementType[iCell] != PENTA_6)
            {
                refineType[iCell] = ISOTROPIC_REFINE;
            }
        }
    }
}

bool ParallelMeshGen::IsNeedInfectInterfaceIsotropicRefineType()
{
    bool isContinue = false;

    int nLocalZones = PHMPI::GetNumberofLocalZones();
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        ZoneGridSolver * zoneGridSolver = zoneGridSolverContainer[ iZone ];

        if(IsNeedInfectInterfaceIsotropicRefineType(zoneGridSolver))
        {
            isContinue = true;
            break;
        }
    }

    return isContinue;
}

bool ParallelMeshGen::IsNeedInfectInterfaceIsotropicRefineType(ZoneGridSolver * zoneGridSolver)
{
    ElementPropertySolver * elementPropertySolver = zoneGridSolver->GetElementPropertySolver();
    ElementFaceSolver * elementFaceSolver = zoneGridSolver->GetElementFaceSolver();

    int nTCell = zoneGridSolver->GetNumberOfElements();
    int nBoundFace = elementFaceSolver->GetNumberOfBoundaryFaces();

    PHVectorInt1D & refineType = * elementPropertySolver->GetRefineType();

    PHVectorInt1D & faceType = elementFaceSolver->GetFaceTypeContainer();
    PHVectorInt1D & faceBoundaryType = elementFaceSolver->GetFaceBoundaryTypeContainer();
    PHVectorInt1D & leftCell  = elementFaceSolver->GetLeftElementIndexContainer();
    PHVectorInt1D & rightCell = elementFaceSolver->GetRightElementIndexContainer();

    bool doseInterfaceIsotropicRefineTypeExist = false;

    // Boundary inter faces.
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        if(faceBoundaryType[iFace] != PHENGLEI::INTERFACE) continue;
        if(faceType[iFace] != QUAD_4) continue;

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

        if(refineType[le] != ISOTROPIC_REFINE && refineType[re] == ISOTROPIC_REFINE)
        {
            doseInterfaceIsotropicRefineTypeExist = true;
            break;
        }
    }

    return doseInterfaceIsotropicRefineTypeExist;
}

void ParallelMeshGen::FixCellRefineType()
{
    InitCellRefineType();

    if(GetDim() == PHSPACE::TWO_D) return;

    int anisoRefine = GlobalDataBase::GetIntParaFromDB("anisoRefine");
    if(anisoRefine != ANISOTROPIC_REFINE) return;

    WriteLogFile("Start to fix cell type!");
    TimeSpan fixTime;

    int infectTime = 0;

    int isContinueInfect = 1;
    while (isContinueInfect)
    {
        isContinueInfect = 0;

        if(infectTime > 0)
        {
            PrintToWindow("\nInfect isotropic refine type time: ", infectTime, "......");
            WriteLogFile("\nInfect isotropic refine type time: ", infectTime, "......");
        }

        int nLocalZones = PHMPI::GetNumberofLocalZones();
        for (int iZone = 0; iZone < nLocalZones; ++ iZone)
        {
            ZoneGridSolver * zoneGridSolver = zoneGridSolverContainer[ iZone ];

            FixCellRefineType(zoneGridSolver);
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

void ParallelMeshGen::FixCellRefineType(ZoneGridSolver * zoneGridSolver)
{
    ElementPropertySolver * elementPropertySolver = zoneGridSolver->GetElementPropertySolver();
    ElementFaceSolver * elementFaceSolver = zoneGridSolver->GetElementFaceSolver();

    int nTCell = zoneGridSolver->GetNumberOfElements();
    int nBoundFace = elementFaceSolver->GetNumberOfBoundaryFaces();
    int nTotalFace = elementFaceSolver->GetNumberOfTotalFaces();

    PHVectorInt2D & cell2Face = * elementPropertySolver->GetElement2Face();
    PHVectorInt1D & refineType = * elementPropertySolver->GetRefineType();
    
    PHVectorInt1D & faceType = elementFaceSolver->GetFaceTypeContainer();
    PHVectorInt1D & faceBoundaryType = elementFaceSolver->GetFaceBoundaryTypeContainer();
    PHVectorInt1D & leftCell  = elementFaceSolver->GetLeftElementIndexContainer();
    PHVectorInt1D & rightCell = elementFaceSolver->GetRightElementIndexContainer();

    deque<int> cellQueue;

    // Boundary inter faces.
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        if(faceType[iFace] != QUAD_4) continue;
        if(faceBoundaryType[iFace] != PHENGLEI::INTERFACE) continue;

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

        if(refineType[le] != ISOTROPIC_REFINE && refineType[re] == ISOTROPIC_REFINE)
        {
            cellQueue.push_back(le);
        }
    }

    // Interior faces.
    for (int iFace = nBoundFace; iFace < nTotalFace; ++ iFace)
    {
        if(faceType[iFace] != QUAD_4) continue;

        int le = leftCell[iFace];
        int re = rightCell[iFace];

        if(refineType[le] == ISOTROPIC_REFINE && refineType[re] != ISOTROPIC_REFINE)
        {
            cellQueue.push_back(re);
        }
        if(refineType[le] != ISOTROPIC_REFINE && refineType[re] == ISOTROPIC_REFINE)
        {
            cellQueue.push_back(le);
        }
    }

    while (!cellQueue.empty())
    {
        int currentCellID = cellQueue.front();
        cellQueue.pop_front();

        if(refineType[currentCellID] == ISOTROPIC_REFINE)
        {
            continue;
        }

        refineType[currentCellID] = ISOTROPIC_REFINE;

        int nFace = static_cast<int>(cell2Face[currentCellID].size());
        for (int iFace = 0; iFace < nFace; ++ iFace)
        {
            int faceID = cell2Face[currentCellID][iFace];

            if(faceID < nBoundFace) continue;  // only for interior face.
            if(faceType[faceID] != QUAD_4) continue;  // only for quad cell.

            int le = leftCell[faceID];
            int re = rightCell[faceID];

            int nextCell = re;
            if(re == currentCellID)
            {
                nextCell = le;
            }

            cellQueue.push_back(nextCell);
        }
    }
}

void ParallelMeshGen::InitGlobalBoundaryCondition()
{
    GlobalBoundaryCondition::ReadGlobalBoundaryCondition();

    this->region->ReSetBoundaryConditionByGlobalBC();
}

void ParallelMeshGen::PreparePointsCoordinates(PointVector *pointsCoord)
{
    int nLocalZones = PHMPI::GetNumberofLocalZones();
    pointsCoord->resize(nLocalZones);
    
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        ZoneGridSolver * zoneGridSolver = zoneGridSolverContainer[iZone];
        Geo_PointFactory * pointFactory = zoneGridSolver->GetPointFactory();

        int numberOfNodes  = static_cast<int>(pointFactory->GetNumberOfComputationalNodes());

        (*pointsCoord)[iZone].clear();
        (*pointsCoord)[iZone].resize(numberOfNodes);

        for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
        {
            double x0 = pointFactory->GetPoint(iNode).X();
            double y0 = pointFactory->GetPoint(iNode).Y();
            double z0 = pointFactory->GetPoint(iNode).Z();

            (*pointsCoord)[iZone][iNode].SetX(x0);
            (*pointsCoord)[iZone][iNode].SetY(y0);
            (*pointsCoord)[iZone][iNode].SetZ(z0);
        }
    }
}

void ParallelMeshGen::PrepareCheckIsSolidWallPoint(vector< vector<bool> > *isSolidWallPoint, vector< vector<bool> > *needDeform)
{
    int nLocalZones = PHMPI::GetNumberofLocalZones();
    isSolidWallPoint->resize(nLocalZones);
    needDeform->resize(nLocalZones);
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        ZoneGridSolver * zoneGridSolver = zoneGridSolverContainer[iZone];

        ElementFaceSolver *elementFaceSolver = zoneGridSolver->GetElementFaceSolver();
        Geo_PointFactory  *pointFactory = zoneGridSolver->GetPointFactory();

        int numberOfNodes  = static_cast<int>(pointFactory->GetNumberOfComputationalNodes());

        int numberOfTotalFaces = elementFaceSolver->GetNumberOfTotalFaces();
        PHVectorInt1D & faceBoundaryTypeContainer = elementFaceSolver->GetFaceBoundaryTypeContainer();
        PHVectorInt2D & faceNodeIndexesContainer  = elementFaceSolver->GetFaceNodeIndexesContainer();

        (*isSolidWallPoint)[iZone].clear();
        (*isSolidWallPoint)[iZone].resize(numberOfNodes, false);

        (*needDeform)[iZone].clear();
        (*needDeform)[iZone].resize(numberOfNodes, true);

        for (int iFace = 0; iFace < numberOfTotalFaces; ++ iFace)
        {
            int bcType = faceBoundaryTypeContainer[iFace];
            if(bcType != PHENGLEI::SOLID_SURFACE) continue;

            for (unsigned iNode = 0; iNode < faceNodeIndexesContainer[iFace].size(); ++ iNode)
            {
                int nodeIndex = faceNodeIndexesContainer[iFace ][ iNode];
                (*isSolidWallPoint)[iZone][nodeIndex] = true;
            }
        }
        
        for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
        {
            if((*isSolidWallPoint)[iZone][iNode])
            {
                (*needDeform)[iZone][iNode] = false;
            }
        }
    }
}

void ParallelMeshGen::PrepareCheckIsSymmetryPoint(vector< vector<bool> > *isSymmetryPoint)
{
    int nLocalZones = PHMPI::GetNumberofLocalZones();
    isSymmetryPoint->resize(nLocalZones);
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        ZoneGridSolver *zoneGridSolver = zoneGridSolverContainer[iZone];

        ElementFaceSolver *elementFaceSolver = zoneGridSolver->GetElementFaceSolver();
        Geo_PointFactory  *pointFactory = zoneGridSolver->GetPointFactory();

        int numberOfNodes  = static_cast<int>(pointFactory->GetNumberOfComputationalNodes());

        int numberOfTotalFaces = elementFaceSolver->GetNumberOfTotalFaces();
        PHVectorInt1D &faceBoundaryTypeContainer = elementFaceSolver->GetFaceBoundaryTypeContainer();
        PHVectorInt2D &faceNodeIndexesContainer  = elementFaceSolver->GetFaceNodeIndexesContainer();

        (*isSymmetryPoint)[iZone].clear();
        (*isSymmetryPoint)[iZone].resize(numberOfNodes, false);

        for (int iFace = 0; iFace < numberOfTotalFaces; ++ iFace)
        {
            int bcType = faceBoundaryTypeContainer[iFace];
            if(bcType != PHENGLEI::SYMMETRY) continue;

            for (unsigned iNode = 0; iNode < faceNodeIndexesContainer[iFace].size(); ++ iNode)
            {
                int nodeIndex = faceNodeIndexesContainer[iFace][iNode];
                (*isSymmetryPoint)[iZone][nodeIndex] = true;
            }
        }
    }
}

void ParallelMeshGen::PrepareBoundaryDeformedPoint(vector< vector<bool> > *isSolidWallPoint, PointVector *zoneProjectVector,
                                                   PointVector *boundaryDeformedPoints, PointVector *deformedVector)
{
    int nLocalZones = PHMPI::GetNumberofLocalZones();
    boundaryDeformedPoints->resize(nLocalZones);
    deformedVector->resize(nLocalZones);
    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        ZoneGridSolver *zoneGridSolver = zoneGridSolverContainer[iZone];

        Geo_PointFactory *pointFactory = zoneGridSolver->GetPointFactory();

        int numberOfNodes = static_cast<int>(pointFactory->GetNumberOfComputationalNodes());

        //! Compute number of wall points.
        int nBoundaryDeformedPoints = 0;
        for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
        {
            if((*isSolidWallPoint)[iZone][iNode])
            {
                ++ nBoundaryDeformedPoints;
            }
        }
        (*boundaryDeformedPoints)[iZone].resize(nBoundaryDeformedPoints);
        (*deformedVector)[iZone].resize(nBoundaryDeformedPoints);

        nBoundaryDeformedPoints = 0; 
        for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
        {
            if((*isSolidWallPoint)[iZone][iNode])
            {
                double dx = (*zoneProjectVector)[iZone][iNode].X();
                double dy = (*zoneProjectVector)[iZone][iNode].Y();
                double dz = (*zoneProjectVector)[iZone][iNode].Z();

                RDouble x0 = pointFactory->GetPoint(iNode).X();
                RDouble y0 = pointFactory->GetPoint(iNode).Y();
                RDouble z0 = pointFactory->GetPoint(iNode).Z();

                Point3D point1 (x0, y0, z0);
                (*boundaryDeformedPoints)[iZone][nBoundaryDeformedPoints] = point1;

                Point3D point2(dx, dy, dz);
                (*deformedVector)[iZone][nBoundaryDeformedPoints] = point2;

                ++ nBoundaryDeformedPoints;
            }
        }
    } 
}

void ParallelMeshGen::UpdatePointsCoordinates(PointVector *pointsCoord, vector< vector<bool> > *needDeform)
{
    int nLocalZones = PHMPI::GetNumberofLocalZones();
    pointsCoord->resize(nLocalZones);

    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        ZoneGridSolver *zoneGridSolver = zoneGridSolverContainer[iZone];
        Geo_PointFactory *pointFactory = zoneGridSolver->GetPointFactory();

        int numberOfNodes  = static_cast<int>(pointFactory->GetNumberOfComputationalNodes());
        for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
        {
            if((*needDeform)[iZone][iNode])
            {
                double xNew = (*pointsCoord)[iZone][iNode].X();
                double yNew = (*pointsCoord)[iZone][iNode].Y();
                double zNew = (*pointsCoord)[iZone][iNode].Z();

                pointFactory->GetPoint(iNode).SetX(xNew);
                pointFactory->GetPoint(iNode).SetY(yNew);
                pointFactory->GetPoint(iNode).SetZ(zNew);
            }
        }
    }
}

void CompressFixCellType(ZoneGridSolver * zoneGridSolver, DataContainer *& dataContainer, int neighborZoneIndex)
{
//  int localZoneIndex = region->GetProcessGlobalZoneIndexToLocalZoneIndex(globalZoneIndex);
//  ZoneGridSolver * zoneGridSolver = this->GetZoneGridSolver(localZoneIndex);
    zoneGridSolver->CompressRefineType(dataContainer, neighborZoneIndex);
}

void DecompressFixCellType(ZoneGridSolver * zoneGridSolver, DataContainer *& dataContainer, int neighborZoneIndex)
{
    // 注意这里的neighbor关系：
    // 此时是neighborZone解析数据，因此是站在neighborZone的角度，而globalZone反而成了他的邻居.
//  int localZoneIndex = region->GetProcessGlobalZoneIndexToLocalZoneIndex(neighborZoneIndex);
//  ZoneGridSolver * zoneGridSolver = this->GetZoneGridSolver(localZoneIndex);
    zoneGridSolver->DecompressRefineType(dataContainer, neighborZoneIndex);
}



}

