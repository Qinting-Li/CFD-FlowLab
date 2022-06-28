#include "UnstructuredOversetConfig.h"
#include "Pre_WalldistCompute.h"
#include "IO_FileName.h"
#include "Glb_Dimension.h"
#include "PHIO.h"
#include "IO_HDF5File.h"
#include "Region.h"
#include "Task.h"
#include "Task_ServerUpdateInterface.h"
#include "OversetKDTree.h"
using namespace std;

namespace PHSPACE
{

OversetConfigFactoryMulti *oversetConfigFactory = NULL;

Region *region;

void RunOversetGridConfig()
{
    if (region == NULL)
    {
        region = new Region;
    }

    region->InitGridForOversetConfig();

    RunUnstructuredOversetGridConfig();

    //delete region;
}

void RunUnstructuredOversetGridConfig()
{
    ReConfigUnstructOversetGrid();

    //OutUnstructOversetGridConfig();

    FreeUnstructOversetGridConfig();
}

void InitializeOversetGridConfig()
{
    if (oversetConfigFactory == NULL)
    {
        oversetConfigFactory = new OversetConfigFactoryMulti();
        oversetConfigFactory->Initialize();
    }
}

void ReConfigUnstructOversetGrid()
{
    InitializeOversetGridConfig();

    oversetConfigFactory->OversetConfigMulti();

    int outTecplotOverset = GlobalDataBase::GetIntParaFromDB("outTecplotOverset");

    if (outTecplotOverset == 1)
    {
        oversetConfigFactory->OutPutOversetFiles();
    }
}

void OutUnstructOversetGridConfig()
{
    oversetConfigFactory->OutPutOversetFiles();
}

void FreeUnstructOversetGridConfig()
{
    if (oversetConfigFactory != NULL)
    {
        delete oversetConfigFactory;
        oversetConfigFactory = NULL;
    }
}

OversetConfigFactoryMulti *GetOversetConfigFactoryMulti()
{
    return oversetConfigFactory;
}

OversetConfigFactoryMulti::OversetConfigFactoryMulti()
{
    numberOfInterfaceNodesInBlocks = NULL;
    keyOfInterfaceNodesInBlocks = NULL;

    numberOfBlocks = 0;
    numberOfZones = 0;

    oversetConfigMultiZone = NULL;
    zonesSearchNodes = NULL;
    zonesInterpolateCells = NULL;

    testTime = 0;
    geometricDimension = GetDim();

    pMinOfBlocks = NULL;
    pMaxOfBlocks = NULL;
}

OversetConfigFactoryMulti::~OversetConfigFactoryMulti()
{
    delete[] numberOfInterfaceNodesInBlocks;
    numberOfInterfaceNodesInBlocks = NULL;

    DelPointer2(keyOfInterfaceNodesInBlocks);
    keyOfInterfaceNodesInBlocks = NULL;

    int currentProcessor = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessor == processorIndex)
        {
            if (oversetConfigMultiZone[iZone]) delete oversetConfigMultiZone[iZone];
        }
        delete zonesInterpolateCells[iZone];
        zonesInterpolateCells[iZone] = NULL;
    }

    delete[] oversetConfigMultiZone;
    delete[] zonesSearchNodes;
    delete[] zonesInterpolateCells;

    for (unsigned iBlock = 0; iBlock < innerGrid.size(); ++iBlock)
    {
        if (innerGrid[iBlock]) delete innerGrid[iBlock];
    }
    innerGrid.resize(0);
}

void OversetConfigFactoryMulti::Initialize()
{
    //TimeSpan *timeSpan = new TimeSpan();
    if (GetIfOutPutOversetVisualization() == 1)
    {
        PHMakeDir("OversetOutput");
    }
    int currentProcessor = PHMPI::GetCurrentProcessorID();
    numberOfBlocks = PHSPACE::GlobalDataBase::GetIntParaFromDB("numberOfGridGroups");

    numberOfZones = PHMPI::GetNumberofGlobalZones();
    oversetConfigMultiZone = new OversetConfig * [numberOfZones];

    //vector< RDouble > *wallNodeCoordinateX = new vector< RDouble >[this->numberOfZones];
    //vector< RDouble > *wallNodeCoordinateY = new vector< RDouble >[this->numberOfZones];
    //vector< RDouble > *wallNodeCoordinateZ = new vector< RDouble >[this->numberOfZones];
    zonesSearchNodes = new ZoneSearchNodes * [numberOfZones];
    zonesInterpolateCells = new InterCells * [numberOfZones];
    for (int iZone = 0; iZone < this->numberOfZones; ++iZone)
    {
        oversetConfigMultiZone[iZone] = NULL;
        zonesSearchNodes[iZone] = NULL;
        zonesInterpolateCells[iZone] = NULL;
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);

        if (currentProcessor == processorIndex)
        {
            OversetConfig *oversetGrid = new OversetConfig();

            Grid *grid = GetGrid(iZone, 0);
            UnstructGrid *gridUnstr = PHSPACE::UnstructGridCast(grid);
            int iBlock = grid->GetIBlock();

            oversetGrid->Initialize(gridUnstr, this->numberOfZones, iZone, iBlock);
            oversetGrid->BuildBoundaryNodeList();
            this->oversetConfigMultiZone[iZone] = oversetGrid;
        }
    }

    //BuildNodesMinDistanceInThisBlocks(wallNodeCoordinateX, wallNodeCoordinateY, wallNodeCoordinateZ);               //! Minimum wall distance of nodes.
    ComputeNodeWalldist();                 //! Minimum wall distance of nodes.

    SetNodesMinDistanceInThisBlocksByWalldist();

    ComputeNodeAndCellTopology();          //! Calculate the topological relationship of point to point and volume to point.

    BuildLocalToGlobalInterfaceNodesMaps();

    ReadInInnerGrid();                     //! Read in the internal grid of the object to accurately estimate whether it is in the object.

    SetBufferRegionOnOversetBoundary();    //! Set buffer on the boundary of overset.
    //timeSpan->ShowTimeSpanToLogfile("OversetConfigFactoryMulti::Initialize");
    //delete timeSpan;
    OutPutOversetVisualizations("0_1_Initialize");
}

void OversetConfigFactoryMulti::BuildLocalToGlobalInterfaceNodesMaps()
{
    int *numberOfNodesInBlocksLocal = new int[numberOfBlocks];
    int *numberOfNodesInBlocks = new int[numberOfBlocks];
    SetField(numberOfNodesInBlocksLocal, 0, numberOfBlocks);
    int currentProcessor = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessor == processorIndex)
        {
            UnstructGrid *grid = oversetConfigMultiZone[iZone]->GetGrid();
            numberOfNodesInBlocksLocal[grid->GetIBlock()] += grid->GetNTotalNode();
        }
    }
    PH_AllReduce(numberOfNodesInBlocksLocal, numberOfNodesInBlocks, numberOfBlocks, MPI_SUM);
    delete[] numberOfNodesInBlocksLocal;
    numberOfNodesInBlocksLocal = NULL;

    numberOfInterfaceNodesInBlocks = new int[numberOfBlocks];
    SetField(numberOfInterfaceNodesInBlocks, 0, numberOfBlocks);
    for (int iBlock = 0; iBlock < numberOfBlocks; ++ iBlock)
    {
        BuildLocalToGlobalInterfaceNodesMap(iBlock, numberOfNodesInBlocks);
    }

    delete[] numberOfNodesInBlocks;
    numberOfNodesInBlocks = NULL;

    keyOfInterfaceNodesInBlocks = NewPointer2< int >(numberOfBlocks, numberOfInterfaceNodesInBlocks);
}

void OversetConfigFactoryMulti::BuildLocalToGlobalInterfaceNodesMap(int iBlockIn, int *numberOfNodesInBlocksIn)
{
    int &numberOfNodesInThisBlock = numberOfNodesInBlocksIn[iBlockIn];

    int *keyOfNodesInBlocks = new int[numberOfNodesInThisBlock];
    int *keyOfNodesInBlocksGlobal = new int[numberOfNodesInThisBlock];
    SetField(keyOfNodesInBlocks, 0, numberOfNodesInThisBlock);
    SetField(keyOfNodesInBlocksGlobal, 0, numberOfNodesInThisBlock);

    int currentProcessor = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessor == processorIndex)
        {
            UnstructGrid *grid = oversetConfigMultiZone[iZone]->GetGrid();
            int blockIndex = grid->GetIBlock();
            if (blockIndex == iBlockIn)
            {
                vector< int > &localToGlobalNodesMap = grid->GetLocalToGlobalNodesMap();
                vector< int > &interfaceBoundaryNodeList = oversetConfigMultiZone[iZone]->GetInterfaceBoundaryNodeList();
                int             numberOfinterfaceBoundaryNodes = interfaceBoundaryNodeList.size();
                for (int iNode = 0; iNode < numberOfinterfaceBoundaryNodes; ++ iNode)
                {
                    keyOfNodesInBlocks[localToGlobalNodesMap[interfaceBoundaryNodeList[iNode]]] = 1;
                }
            }
        }
    }

    PH_AllReduce(keyOfNodesInBlocks, keyOfNodesInBlocksGlobal, numberOfNodesInThisBlock, MPI_MAX);

    for (int iNode = 0; iNode < numberOfNodesInThisBlock; ++ iNode)
    {
        if (keyOfNodesInBlocksGlobal[iNode])
        {
            keyOfNodesInBlocks[iNode] = (numberOfInterfaceNodesInBlocks[iBlockIn]) ++;
        }
    }

    delete[] keyOfNodesInBlocksGlobal;
    keyOfNodesInBlocksGlobal = NULL;

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessor == processorIndex)
        {
            UnstructGrid *grid = oversetConfigMultiZone[iZone]->GetGrid();
            int blockIndex = grid->GetIBlock();
            if (blockIndex == iBlockIn)
            {
                vector< int > &localToGlobalNodesMap = grid->GetLocalToGlobalNodesMap();
                int *localToGlobalInterfaceNodesMap = oversetConfigMultiZone[iZone]->GetLocalToGlobalInterfaceNodesMap();
                vector< int > &interfaceBoundaryNodeList = oversetConfigMultiZone[iZone]->GetInterfaceBoundaryNodeList();
                int             numberOfinterfaceBoundaryNodes = interfaceBoundaryNodeList.size();

                for (int iNode = 0; iNode < numberOfinterfaceBoundaryNodes; ++ iNode)
                {
                    localToGlobalInterfaceNodesMap[iNode] = keyOfNodesInBlocks[localToGlobalNodesMap[interfaceBoundaryNodeList[iNode]]];
                }
            }
        }
    }

    delete[] keyOfNodesInBlocks;
    keyOfNodesInBlocks = NULL;
}

void OversetConfigFactoryMulti::BuildNodesMinDistanceInThisBlocks(vector< RDouble > *&wallNodeCoordinateXIn,
    vector< RDouble > *&wallNodeCoordinateYIn,
    vector< RDouble > *&wallNodeCoordinateZIn)
{
    vector< RDouble > *wallNodeCoordinateXInBlocks = new vector< RDouble >[numberOfBlocks];
    vector< RDouble > *wallNodeCoordinateYInBlocks = new vector< RDouble >[numberOfBlocks];
    vector< RDouble > *wallNodeCoordinateZInBlocks = new vector< RDouble >[numberOfBlocks];

    int blockIndexOfThisZone = 0;
    int numberOfWallNodesInThisZones = 0;
    int currentProcessor = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessor == processorIndex)
        {
            blockIndexOfThisZone = oversetConfigMultiZone[iZone]->GetGrid()->GetIBlock();
            numberOfWallNodesInThisZones = wallNodeCoordinateXIn[iZone].size();
        }
        PHMPI::PH_Bcast(&numberOfWallNodesInThisZones, sizeof(int), processorIndex);

        if (numberOfWallNodesInThisZones)
        {
            if (currentProcessor != processorIndex)
            {
                wallNodeCoordinateXIn[iZone].resize(numberOfWallNodesInThisZones);
                wallNodeCoordinateYIn[iZone].resize(numberOfWallNodesInThisZones);
                wallNodeCoordinateZIn[iZone].resize(numberOfWallNodesInThisZones);
            }

            PHMPI::PH_Bcast(& wallNodeCoordinateXIn[iZone][0], numberOfWallNodesInThisZones * sizeof(RDouble), processorIndex);
            PHMPI::PH_Bcast(& wallNodeCoordinateYIn[iZone][0], numberOfWallNodesInThisZones * sizeof(RDouble), processorIndex);
            PHMPI::PH_Bcast(& wallNodeCoordinateZIn[iZone][0], numberOfWallNodesInThisZones * sizeof(RDouble), processorIndex);

            PHMPI::PH_Bcast(& blockIndexOfThisZone, sizeof(int), processorIndex);
            for (int iNode = 0; iNode < numberOfWallNodesInThisZones; ++ iNode)
            {
                wallNodeCoordinateXInBlocks[blockIndexOfThisZone].push_back(wallNodeCoordinateXIn[iZone][iNode]);
                wallNodeCoordinateYInBlocks[blockIndexOfThisZone].push_back(wallNodeCoordinateYIn[iZone][iNode]);
                wallNodeCoordinateZInBlocks[blockIndexOfThisZone].push_back(wallNodeCoordinateZIn[iZone][iNode]);
            }
        }
    }

    delete[] wallNodeCoordinateXIn;
    delete[] wallNodeCoordinateYIn;
    delete[] wallNodeCoordinateZIn;
    wallNodeCoordinateXIn = NULL;
    wallNodeCoordinateYIn = NULL;
    wallNodeCoordinateZIn = NULL;

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessor == processorIndex)
        {
            oversetConfigMultiZone[iZone]->BuildNodesMinDistanceInThisBlock(wallNodeCoordinateXInBlocks,
                wallNodeCoordinateYInBlocks,
                wallNodeCoordinateZInBlocks);
        }
    }

    delete[] wallNodeCoordinateXInBlocks;
    delete[] wallNodeCoordinateYInBlocks;
    delete[] wallNodeCoordinateZInBlocks;
    wallNodeCoordinateXInBlocks = NULL;
    wallNodeCoordinateYInBlocks = NULL;
    wallNodeCoordinateZInBlocks = NULL;
}
void OversetConfigFactoryMulti::SetNodesMinDistanceInThisBlocksByWalldist()
{
    int currentProcessor = PHMPI::GetCurrentProcessorID();

    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessor == processorIndex)
        {
            oversetConfigMultiZone[iZone]->SetNodesMinDistanceInThisBlocksByWalldist();
        }
    }
}

void OversetConfigFactoryMulti::SetNodesMinDistanceInOtherBlocks()
{
    int currentProcessor = PHMPI::GetCurrentProcessorID();

    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessor == processorIndex)
        {
            oversetConfigMultiZone[iZone]->SetNodesMinDistanceInOtherBlocks();
        }
    }
}

bool OversetConfigFactoryMulti::ReadInInnerGridH5(const string &fileName, UnstructGrid *gridIn)
{
    bool HDF5FileRead = false;

    fstream file;
    OpenFile(file, fileName, ios_base::in | ios_base::binary);
    int nZone;
    int sizeInt = sizeof(int);
    file.read((char *)&nZone, sizeInt);
    CloseFile(file);

    if (nZone <= static_cast <int> ( 1e8 ))
    {
        return HDF5FileRead;
    }

    hid_t h5File, grpGrid, grpData;
    string grpName;
    h5File = OpenHDF5File(fileName);

    ostringstream oss;
    oss << "Grid-" << "0";
    grpName = oss.str();
    grpGrid = OpenGroup(h5File, grpName);

    int *iDimensions = ReadIntOneRow(grpGrid, "iDimensions");
    gridIn->SetNTotalNode(iDimensions[0]);
    gridIn->SetNTotalFace(iDimensions[1]);
    gridIn->SetNTotalCell(iDimensions[2]);
    delete[] iDimensions;

    //! Read Grid Coordinates.
    grpData = OpenGroup(grpGrid, "GridCoordinates");

    RDouble *CoordinateX = ReadDoubleOneRow(grpData, "CoordinateX");
    RDouble *CoordinateY = ReadDoubleOneRow(grpData, "CoordinateY");
    RDouble *CoordinateZ = ReadDoubleOneRow(grpData, "CoordinateZ");
    gridIn->SetX(CoordinateX);
    gridIn->SetY(CoordinateY);
    gridIn->SetZ(CoordinateZ);

    int axisup = 1;
    PHSPACE::GlobalDataBase::GetData("axisup", &axisup, PHINT, 1);
    if (axisup == 2)
    {
        gridIn->TurnZAxis2YAxis();
    }
    gridIn->ComputeMinMaxBox();

    H5Gclose(grpData);

    //! Read Face Topology.
    grpData = OpenGroup(grpGrid, "FaceTopology");

    int *nodeNumberOfEachFace = ReadIntOneRow(grpData, "nodeNumberOfEachFace");
    gridIn->SetNodeNumberOfEachFace(nodeNumberOfEachFace);

    int *face2Node = ReadIntOneRow(grpData, "face2Node");
    gridIn->SetFace2Node(face2Node);

    int *leftCellOfFace = ReadIntOneRow(grpData, "leftCellOfFace");
    int *rightCellOfFace = ReadIntOneRow(grpData, "rightCellOfFace");
    gridIn->SetLeftCellOfFace(leftCellOfFace);
    gridIn->SetRightCellOfFace(rightCellOfFace);

    int nTotalFace = gridIn->GetNTotalFace();
    int *nodePosi = new int[nTotalFace + 1];
    int nodeSum = 0;
    nodePosi[0] = 0;
    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        nodeSum += nodeNumberOfEachFace[iFace];
        nodePosi[iFace + 1] = nodeSum;
    }

    for (int iFace = 0; iFace < nTotalFace; ++ iFace)
    {
        if (leftCellOfFace[iFace] < 0)
        {
            std::reverse(face2Node + nodePosi[iFace], face2Node + nodePosi[iFace + 1]);
            SWAP(leftCellOfFace[iFace], rightCellOfFace[iFace]);
        }
    }
    delete[] nodePosi;

    int nBoundFace = 0;
    ReadData(grpData, &nBoundFace, "nBoundFace");
    gridIn->SetNBoundFace(nBoundFace);

    int *BcType = ReadIntOneRow(grpData, "BcType");
    UnstructBCSet **bcr = new UnstructBCSet * [nBoundFace];
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        bcr[iFace] = new UnstructBCSet();
        bcr[iFace]->SetKey(BcType[iFace]);
    }
    gridIn->SetBCRecord(bcr);
    delete[] BcType;

    gridIn->SpecifyRightCellofBC();

    H5Gclose(grpData);

    //! Read Interface Info.
    grpData = OpenGroup(grpGrid, "InterfaceInfo");

    InterfaceInfo *interfaceInfo = 0;
    gridIn->SetInterfaceInfo(interfaceInfo);

    int nIFace = 0;
    ReadData(grpData, &nIFace, "nIFace");

    if (nIFace > 0)
    {
        interfaceInfo = new InterfaceInfo(nIFace, gridIn);
        gridIn->SetInterfaceInfo(interfaceInfo);

        int *interFace2ZoneID = interfaceInfo->GetInterFace2ZoneID();
        int *interFace2InterFaceID = interfaceInfo->GetInterFace2InterFaceID();
        int *interFace2BoundaryFace = interfaceInfo->GetInterFace2BoundaryFace();

        int *interFace2ZoneIDTemp = ReadIntOneRow(grpData, "interFace2ZoneID");
        int *interFace2InterFaceIDTemp = ReadIntOneRow(grpData, "interFace2InterFaceID");
        int *interFace2BoundaryFaceTemp = ReadIntOneRow(grpData, "interFace2BoundaryFace");

        for (int iFace = 0; iFace < nIFace; ++ iFace)
        {
            interFace2ZoneID[iFace] = interFace2ZoneIDTemp[iFace];
            interFace2InterFaceID[iFace] = interFace2InterFaceIDTemp[iFace];
            interFace2BoundaryFace[iFace] = interFace2BoundaryFaceTemp[iFace];
        }

        delete[] interFace2ZoneIDTemp;
        delete[] interFace2InterFaceIDTemp;
        delete[] interFace2BoundaryFaceTemp;
    }

    H5Gclose(grpData);

    //! Read Cell Topology.
    grpData = OpenGroup(grpGrid, "CellTopology");

    int *cell2Node = ReadIntOneRow(grpData, "cell2Node");
    int *nodeNumberOfEachCell = ReadIntOneRow(grpData, "nodeNumberOfEachCell");

    gridIn->SetNodeNumberOfEachCell(nodeNumberOfEachCell);
    gridIn->SetCell2Node(cell2Node);

    H5Gclose(grpData);

    //! Read Unstr BCName.
    grpData = OpenGroup(grpGrid, "BCName");

    char *BCName = ReadCharData(grpData, "BCName");

    int nBFaces = gridIn->GetNBoundFace();
    string *BCNameList = new string[nBFaces];
    uint_long count = 0;
    for (int iFace = 0; iFace < nBFaces; ++ iFace)
    {
        while (BCName[count] != '\0')
        {
            BCNameList[iFace].push_back(BCName[count]);
            ++count;
        }
        ++count;
    }

    for (int iFace = 0; iFace < nBFaces; ++ iFace)
    {
        bcr[iFace]->SetBCName(BCNameList[iFace]);
    }

    delete[] BCNameList;
    delete[] BCName;
    H5Gclose(grpData);

    H5Gclose(grpGrid);
    H5Fclose(h5File);

    HDF5FileRead = true;
    return HDF5FileRead;
}

void OversetConfigFactoryMulti::ReadInInnerGrid()
{
    if (GetReadInAuxiliaryInnerGrid() != 1)
    {
        return;
    }

    int currentProcessor = PHMPI::GetCurrentProcessorID();
    int serverProcessor = PHMPI::GetServerSepMode();

    innerGrid.resize(numberOfBlocks);

    for (int iBlock = 0; iBlock < numberOfBlocks; ++iBlock)
    {
        ostringstream oss;

        oss << "auxiliaryInnerGrid" << iBlock;

        int fileID = PHMPI::GetFileIndexofCurrentProcessor();
        string fileName = AddSymbolToFileName(GlobalDataBase::GetStrParaFromDB(oss.str()), "_", fileID);

        if (!FileExist(fileName))
        {
            innerGrid[iBlock] = NULL;
            continue;
        }

        UnstructGrid *grid = new UnstructGrid();
        grid->InitGrid(new GridID(0, 0), 0, GetDim(), UNSTRUCTGRID);
        innerGrid[iBlock] = grid;
        innerGrid[iBlock]->SetIBlock(iBlock);

        bool fileIsH5File = false;
        fileIsH5File = ReadInInnerGridH5(fileName, grid);
        if (fileIsH5File)
        {
            //ostringstream oss;
            //oss << "./results/auxiliaryInnerGrid" << iBlock;
            //string fileName = oss.str();
            //
            //innerGrid[iBlock] ->Output(fileName);
            continue;
        }

        fstream file;

        DataContainer *dataContainer = new DataContainer;

        if (currentProcessor == serverProcessor)
        {
            cout << "Read in auxiliary inner grid" << fileName << endl;

            OpenFile(file, fileName, ios_base::in | ios_base::binary);

            int numberOfZones;

            vector < int > zoneProcessorIndexContainer;
            vector < int > zoneIndexContainer;
            vector < int > zoneTypeContainer;

            PHRead(file, numberOfZones);

            zoneProcessorIndexContainer.resize(numberOfZones);
            zoneIndexContainer.resize(numberOfZones);
            zoneTypeContainer.resize(numberOfZones);

            PHRead(file, zoneProcessorIndexContainer, numberOfZones);
            PHRead(file, zoneIndexContainer, numberOfZones);
            PHRead(file, zoneTypeContainer, numberOfZones);

            dataContainer->readfile(file);

            CloseFile(file);
        }

        int tag = 0;

        PHMPI::PH_BcastByServerSepMode(dataContainer, tag);

        innerGrid[iBlock]->Decode(dataContainer, 0);

        delete dataContainer;
    }
}

void OversetConfigFactoryMulti::ComputeNodeWalldist()
{
    int walldistMethod = 3;

    int level = 0;

    Pre_WalldistCompute *wallDistCompute = new Pre_WalldistCompute(walldistMethod, level);

    wallDistCompute->ComputeNodeDistance();

    delete wallDistCompute;
}

void OversetConfigFactoryMulti::ComputeNodeWalldistInOtherBlock()
{
    int walldistMethod = 3;

    int level = 0;

    Pre_WalldistCompute *wallDistCompute = new Pre_WalldistCompute(walldistMethod, level);

    wallDistCompute->ComputeNodeDistanceInOtherBlock();

    delete wallDistCompute;
}

void OversetConfigFactoryMulti::ComputeNodeAndCellTopology()
{
    int currentProcessor = PHMPI::GetCurrentProcessorID();

    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessor == processorIndex)
        {
            oversetConfigMultiZone[iZone]->ComputeNodeAndCellTopology();
        }
    }
}

//! Get the overall globalInterfaceNodes according to the localToGlobalNodeList information of each zone.
void OversetConfigFactoryMulti::GetGlobalInterfaceNodes()
{
    if (numberOfZones == numberOfBlocks)
    {
        return;
    }

    for (int iBlock = 0; iBlock < numberOfBlocks; ++ iBlock)
    {
        GetGlobalInterfaceNodes(iBlock);
    }
}

void OversetConfigFactoryMulti::GetGlobalInterfaceNodes(int iBlock)
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();

    //! step1: Distinguish globalInterNodes.
    vector< int > keyNodesGlobal;
    vector< int > keyNodesGlobalAllProcessor;

    int numberOfNodesInThisBlock = numberOfGlobalNodes[iBlock];

    keyNodesGlobal.resize(numberOfNodesInThisBlock);
    keyNodesGlobalAllProcessor.resize(numberOfNodesInThisBlock);

    SetField(keyNodesGlobal, 0, numberOfNodesInThisBlock);
    SetField(keyNodesGlobalAllProcessor, 0, numberOfNodesInThisBlock);

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int sendProcessorIndex = PHMPI::GetZoneProcessorIDSepMode(iZone);
        if (currentProcessorIndex == sendProcessorIndex)
        {
            oversetConfigMultiZone[iZone]->MarkGlobalInterfaceNodes(iBlock, keyNodesGlobal);
        }
    }

    PH_AllreduceSepModeForVector(keyNodesGlobal, keyNodesGlobalAllProcessor, numberOfNodesInThisBlock, MPI_INT, MPI_MAX);

    //! step2£ºkeyNodesGlobal becomes the number which points to the new sequence.
    int numberOfGlobalInterfaceNodes = 0;
    for (int iNode = 0; iNode < numberOfNodesInThisBlock; ++ iNode)
    {
        if (keyNodesGlobalAllProcessor[iNode] == 1)
        {
            keyNodesGlobal[iNode] = numberOfGlobalInterfaceNodes;
            numberOfGlobalInterfaceNodes++;
        }
    }

    //! step3£ºGet localToGlobalInterfaceNodes.
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int sendProcessorIndex = PHMPI::GetZoneProcessorIDSepMode(iZone);
        if (currentProcessorIndex == sendProcessorIndex)
        {
            oversetConfigMultiZone[iZone]->CalLocalToGlobalInterfaceNodes(iBlock, keyNodesGlobal);
        }
    }

    keyActiveOfGlobalInterNodes[iBlock].resize(numberOfGlobalInterfaceNodes);
    SetField(&keyActiveOfGlobalInterNodes[iBlock][0], 0, numberOfGlobalInterfaceNodes);
}

void OversetConfigFactoryMulti::OversetConfigMulti()
{
    SetInitialValue();    //! Set initial value.

    NodesOperator();

    CellsOperator();
}

void OversetConfigFactoryMulti::SpecifyCellsInBodyByAuxiliaryInnerGrid()
{
    //TimeSpan *timeSpan = new TimeSpan();
    vector< GridKDTree * > innerAuxiliaryGridKDTrees;

    innerAuxiliaryGridKDTrees.resize(numberOfBlocks);

    //! step1£ºInteral grid establishment ADT-Tree.
    for (int iBlock = 0; iBlock < numberOfBlocks; ++iBlock)
    {
        //! Set CurrentActionKey to 0, avoid unnecessary error in ComputMetrices.
        PHSPACE::SetCurrentActionKey(0);

        if (innerGrid[iBlock])
        {
            innerGrid[iBlock]->ComputeMetrics();
            innerAuxiliaryGridKDTrees[iBlock] = new GridKDTree(innerGrid[iBlock]);
        }
        else
        {
            innerAuxiliaryGridKDTrees[iBlock] = NULL;
        }
    }

    BuildInBlockCellsAndNodesByInnerAuxiliaryGridKDTrees(innerAuxiliaryGridKDTrees);

    for (int iBlock = 0; iBlock < numberOfBlocks; ++iBlock)
    {
        if (innerAuxiliaryGridKDTrees[iBlock] != 0)
        {
            delete innerAuxiliaryGridKDTrees[iBlock];
        }
    }
    //timeSpan->ShowTimeSpanToLogfile("SpecifyCellsInBodyByAuxiliaryInnerGrid");
    //delete timeSpan;
    OutPutOversetVisualizations("1_3_BuildInBlockCellsAndNodes");
}

void OversetConfigFactoryMulti::BuildInBlockCellsAndNodesByInnerAuxiliaryGridKDTrees(vector< GridKDTree * > &innerAuxiliaryGridKDTreesIn)
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->BuildInBlockCellsAndNodesByInnerAuxiliaryGridKDTree(innerAuxiliaryGridKDTreesIn);
        }
    }
}

//! Expand a row of active points, points in the innerBody and buffer are set to inactive.
void OversetConfigFactoryMulti::InterFaceBoundaryOptimize()
{
    //TimeSpan *timeSpan = new TimeSpan();
    int twoOrderInterpolationOrNot = GetTwoOrderInterpolationOrNot() + 1;

    for (int iter = 0; iter < twoOrderInterpolationOrNot; ++ iter)
    {
        CharacterNodeToCell();
        CharacterCellToNode();
        SetNegativeForInnerNodes();
    }

    //! Points in buffer are set to inactive.
    SetNegativeForBufferRegion();
    //timeSpan->ShowTimeSpanToLogfile("InterFaceBoundaryOptimize");
    //delete timeSpan;
    OutPutOversetVisualizations("1_5_InterFaceBoundaryOptimize");
}

void OversetConfigFactoryMulti::SetNegativeForInnerNodes()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->SetNegativeForInnerNodes();
        }
    }
}

void OversetConfigFactoryMulti::SetNegativeForInnerAndInitialInterNodes()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->SetNegativeForInnerAndInitialInterNodes();
        }
    }
}

void OversetConfigFactoryMulti::NodesOperator()
{
    BuildNodesMinDistanceInOtherBlocks();

    if (GetReadInAuxiliaryInnerGrid() == 1)
    {
        SpecifyCellsInBodyByAuxiliaryInnerGrid();
    }
    else
    {
        SpecifySpecicalCellsInteractWithWallBoundary();
    }

    SpecifyInitialInterNodesByDistance();

    FillTheActiveRegion();          //! Physical boundary advancement.

    InterFaceBoundaryOptimize();    //! Hole boundary optimization (expand active areas).
}

void OversetConfigFactoryMulti::ReComputeWallDistance()
{
    int numberOfMovingBodies = GlobalDataBase::GetIntParaFromDB("numberOfMovingBodies");

    int key = 0;

    for (int iBody = 0; iBody < numberOfMovingBodies; ++iBody)
    {
        ostringstream oss;

        oss << "morphing_" << iBody;

        int morphingManner = GlobalDataBase::GetIntParaFromDB(oss.str());

        if (morphingManner != 0)
        {
            key = 1;
            break;
        }
    }

    //! No need to solve the wall distance again if the overall objects don't deform.
    if (key == 0)
    {
        return;
    }

    ComputeNodeWalldist();
}

void OversetConfigFactoryMulti::BuildNodesMinDistanceInOtherBlocks()
{
    //TimeSpan *timeSpan = new TimeSpan();
    BuildKDTreesBySearchRegions();

    BuildNodesMinDistanceInOtherBlocksByKDTrees();
    //timeSpan->ShowTimeSpanToLogfile("BuildNodesMinDistanceInOtherBlocks");
    //delete timeSpan;
    OutPutOversetVisualizations("1_1_BuildNodesMinDistanceInOtherBlocks");
}

void OversetConfigFactoryMulti::BuildNodesMinDistanceInOtherBlocks2()
{
    //TimeSpan *timeSpan = new TimeSpan();
    //timeSpan->ShowTimeSpanToLogfile("BuildNodesMinDistanceInOtherBlocks");
    //delete timeSpan;
    BuildKDTreesBySearchRegions();

    BuildZonesSearchNodes();

    ComputeNodeWalldistInOtherBlock();

    SetNodesMinDistanceInOtherBlocks();

    OutPutOversetVisualizations("1_1_BuildNodesMinDistanceInOtherBlocks");
}

void OversetConfigFactoryMulti::BuildNodesMinDistanceInOtherBlocksByKDTrees()
{
    //TimeSpan *timeSpan = new TimeSpan();
    BuildZonesSearchNodes();
    //timeSpan->ShowTimeSpan("BuildZonesSearchNodes ");
    //delete timeSpan;

    //timeSpan = new TimeSpan();
    SearchZonesSearchNodesInKDTrees();
    //timeSpan->ShowTimeSpan("SearchZonesSearchNodesInKDTrees ");
    //delete timeSpan;

    AllReduceNodesMinDistanceInOtherBlock();

    SetNodesMinDistanceInOtherBlockByZonesNodes();
}

void OversetConfigFactoryMulti::SearchZonesSearchNodesInKDTrees()
{
    int currentProcessor = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessor == processorIndex)
        {
            oversetConfigMultiZone[iZone]->SearchZonesSearchNodesInKDTree();
        }
    }
}

void OversetConfigFactoryMulti::AllReduceNodesMinDistanceInOtherBlock()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        ZoneSearchNodes *&zoneSearchNodes = zonesSearchNodes[iZone];
        int numberOfNodes = zoneSearchNodes->GetNumberOfNodes();
        RDouble *nodesMinDistanceToOtherBlock = zoneSearchNodes->GetNodesMinDistanceToOtherBlock();

        if (numberOfNodes)
        {
            RDouble *nodesMinDistanceInOtherBlockGlobal = new RDouble[numberOfNodes];

            PH_AllReduce(nodesMinDistanceToOtherBlock, nodesMinDistanceInOtherBlockGlobal, numberOfNodes, PH_MIN);

            SetField(nodesMinDistanceToOtherBlock, nodesMinDistanceInOtherBlockGlobal, numberOfNodes);

            delete[] nodesMinDistanceInOtherBlockGlobal;
            nodesMinDistanceInOtherBlockGlobal = NULL;
        }
    }
}

void OversetConfigFactoryMulti::SetNodesMinDistanceInOtherBlockByZonesNodes()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        ZoneSearchNodes *&zoneSearchNodes = zonesSearchNodes[iZone];
        if (currentProcessorIndex == processorIndex)
        {
            RDouble *nodesMinDistanceToOtherBlock = zoneSearchNodes->GetNodesMinDistanceToOtherBlock();
            oversetConfigMultiZone[iZone]->CopyMinDist(nodesMinDistanceToOtherBlock);
        }
        delete zonesSearchNodes[iZone];
        zonesSearchNodes[iZone] = NULL;
    }
}

void OversetConfigFactoryMulti::CellsOperator()
{
    FindActiveCells();    //! Identify active cells.

    FindInterCells();     //! Identify interpolation type.

    BuildInterpolateRelationships();

    CheckInterpolateRelationships();

    ResetOversetInformation();    //! Rebuild interpolation information.
}

void OversetConfigFactoryMulti::SetInitialValue()
{
    int currentProcessor = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessor == processorIndex)
        {
            oversetConfigMultiZone[iZone]->SetInitialValue();
        }
        if (NULL != zonesInterpolateCells[iZone])
        {
            delete zonesInterpolateCells[iZone];
            zonesInterpolateCells[iZone] = NULL;
        }
    }
}

void OversetConfigFactoryMulti::CalMinBoxForSearch()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();

    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->BuildMinMaxBoxOfZoneBySearchRegion();
        }
    }
}

void OversetConfigFactoryMulti::SetBufferRegionOnOversetBoundary()
{
    //! Query physical boundary points and overset boundary points.
    int currentProcessor = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessor == processorIndex)
        {
            oversetConfigMultiZone[iZone]->SetInitialValueForBufferRigion();
        }
    }

    CharacterNodeToCell();    //! A row of points and cells are left in the buffer.
    //! Increase the interpolation area to improve the reconstruction accuracy of interpolation cells, two row of points and cells will be left in the buffer (There's no necessary actually).
    int keyEnlargeOfActiveNodes = GetKeyEnlargeOfNodes();
    for (int iter = 0; iter < keyEnlargeOfActiveNodes; ++ iter)
    {
        CharacterCellToNode();

        CharacterNodeToCell();
    }

    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessor == processorIndex)
        {
            oversetConfigMultiZone[iZone]->SetKeyBufferByKeyActive();
        }
    }
}

//! Identify searchRegion, starting from search point and expand twice.
//! Identify search point by judging whether it's in other blocks.
void OversetConfigFactoryMulti::BuildKDTreesBySearchRegions()
{
    //TimeSpan *timeSpan = new TimeSpan();
    BuildMinMaxBoxOfBlocks();

    SetActivePropertyForSearchRegions();

    //! Why the value is 2 (Ensure that searchRegion covers all blocks except this block).
    for (int i = 0; i < 2; ++i)
    {
        CharacterNodeToCell();

        CharacterCellToNode();
    }

    SetSearchRegionByActivePropertys();

    BuildMinMaxBoxOfZoneBySearchRegions();

    //TimeSpan *timeSpan = new TimeSpan();
    BuildKDTreeByMinMaxBoxOfZones();
    //timeSpan->ShowTimeSpan("BuildKDTreeByMinMaxBoxOfZones");
    //delete timeSpan;
}

void OversetConfigFactoryMulti::SetActivePropertyForSearchRegions()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->SetActivePropertyForSearchRegion();
        }
    }
    DelPointer2(pMinOfBlocks);
    DelPointer2(pMaxOfBlocks);
    pMinOfBlocks = NULL;
    pMaxOfBlocks = NULL;
}

void OversetConfigFactoryMulti::SetSearchRegionByActivePropertys()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->SetSearchRegionByActiveProperty();
        }
    }
}

void OversetConfigFactoryMulti::BuildMinMaxBoxOfZoneBySearchRegions()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->BuildMinMaxBoxOfZoneBySearchRegion();
        }
    }
}

void OversetConfigFactoryMulti::BuildKDTreeByMinMaxBoxOfZones()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->BuildKDTreeByMinMaxBoxOfZone();
        }
    }
}

void OversetConfigFactoryMulti::BuildMinMaxBoxOfBlocks()
{
    pMinOfBlocks = NewPointer2< RDouble >(numberOfBlocks, geometricDimension);
    pMaxOfBlocks = NewPointer2< RDouble >(numberOfBlocks, geometricDimension);
    RDouble **minBoxOfBlocksLocal = NewPointer2< RDouble >(numberOfBlocks, geometricDimension);
    RDouble **maxBoxOfBlocksLocal = NewPointer2< RDouble >(numberOfBlocks, geometricDimension);

    for (int iBlock = 0; iBlock < numberOfBlocks; ++iBlock)
    {
        SetField(pMinOfBlocks[iBlock], LARGE, geometricDimension);
        SetField(pMaxOfBlocks[iBlock], -LARGE, geometricDimension);
        SetField(minBoxOfBlocksLocal[iBlock], LARGE, geometricDimension);
        SetField(maxBoxOfBlocksLocal[iBlock], -LARGE, geometricDimension);
    }

    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->BuildMinMaxBoxOfZone();

            RDouble *minBoxOfZone = oversetConfigMultiZone[iZone]->GetPMin();
            RDouble *maxBoxOfZone = oversetConfigMultiZone[iZone]->GetPMax();
            int blockIndex = oversetConfigMultiZone[iZone]->GetGrid()->GetIBlock();
            for (int iDimension = 0; iDimension < geometricDimension; ++iDimension)
            {
                RDouble &minInThisDimension = minBoxOfBlocksLocal[blockIndex][iDimension];
                RDouble &maxInThisDimension = maxBoxOfBlocksLocal[blockIndex][iDimension];
                minInThisDimension = MIN(minInThisDimension, minBoxOfZone[iDimension]);
                maxInThisDimension = MAX(maxInThisDimension, maxBoxOfZone[iDimension]);
            }
        }
    }

    RDouble eps = GetToleranceForOversetBox();
    for (int iBlock = 0; iBlock < numberOfBlocks; ++iBlock)
    {
        PH_AllReduce(minBoxOfBlocksLocal[iBlock], pMinOfBlocks[iBlock], geometricDimension, PH_MIN);
        PH_AllReduce(maxBoxOfBlocksLocal[iBlock], pMaxOfBlocks[iBlock], geometricDimension, PH_MAX);

        for (int iDimension = 0; iDimension < geometricDimension; ++iDimension)
        {
            pMinOfBlocks[iBlock][iDimension] -= eps;
            pMaxOfBlocks[iBlock][iDimension] += eps;
        }
    }

    DelPointer2(minBoxOfBlocksLocal);
    DelPointer2(maxBoxOfBlocksLocal);
    minBoxOfBlocksLocal = NULL;
    maxBoxOfBlocksLocal = NULL;
}

void OversetConfigFactoryMulti::BuildZonesSearchNodes()
{
    int currentProcessor = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        int *keySearchOfNodes = NULL;
        int   numberOfNodes = 0;

        int   blockIndex = 0;
        int   numberOfSearchNodes = 0;

        if (currentProcessor == processorIndex)
        {
            keySearchOfNodes = oversetConfigMultiZone[iZone]->GetKeySearchOfNodes();
            numberOfNodes = oversetConfigMultiZone[iZone]->GetGrid()->GetNTotalNode();

            blockIndex = oversetConfigMultiZone[iZone]->GetGrid()->GetIBlock();

            for (int iNode = 0; iNode < numberOfNodes; ++iNode)
            {
                if (keySearchOfNodes[iNode])
                {
                    ++ numberOfSearchNodes;
                }
            }
        }

        PHMPI::PH_Bcast(&numberOfSearchNodes, sizeof(int), processorIndex, 1);

        PHMPI::PH_Bcast(&blockIndex, sizeof(int), processorIndex, 1);

        if(zonesSearchNodes[iZone])
        {
            delete zonesSearchNodes[iZone];
            zonesSearchNodes[iZone] = NULL;
        }
        zonesSearchNodes[iZone] = new ZoneSearchNodes(numberOfSearchNodes, blockIndex);

        RDouble *minBoxOfZoneSearchNodes = zonesSearchNodes[iZone]->GetPMin();
        RDouble *maxBoxOfZoneSearchNodes = zonesSearchNodes[iZone]->GetPMax();

        if (currentProcessor == processorIndex)
        {
            RDouble *minBoxOfZone = oversetConfigMultiZone[iZone]->GetPMin();
            RDouble *maxBoxOfZone = oversetConfigMultiZone[iZone]->GetPMax();

            SetField(minBoxOfZoneSearchNodes, minBoxOfZone, geometricDimension);
            SetField(maxBoxOfZoneSearchNodes, maxBoxOfZone, geometricDimension);
        }
        PHMPI::PH_Bcast(minBoxOfZoneSearchNodes, geometricDimension * sizeof(RDouble), processorIndex);
        PHMPI::PH_Bcast(maxBoxOfZoneSearchNodes, geometricDimension * sizeof(RDouble), processorIndex);

        if (numberOfSearchNodes)
        {
            RDouble *searchNodeCoordinateX = zonesSearchNodes[iZone]->GetX();
            RDouble *searchNodeCoordinateY = zonesSearchNodes[iZone]->GetY();
            RDouble *searchNodeCoordinateZ = zonesSearchNodes[iZone]->GetZ();

            if (currentProcessor == processorIndex)
            {
                RDouble *nodeCoordinateX = oversetConfigMultiZone[iZone]->GetGrid()->GetX();
                RDouble *nodeCoordinateY = oversetConfigMultiZone[iZone]->GetGrid()->GetY();
                RDouble *nodeCoordinateZ = oversetConfigMultiZone[iZone]->GetGrid()->GetZ();

                numberOfSearchNodes = 0;

                for (int iNode = 0; iNode < numberOfNodes; ++iNode)
                {
                    if (keySearchOfNodes[iNode])
                    {
                        searchNodeCoordinateX[numberOfSearchNodes] = nodeCoordinateX[iNode];
                        searchNodeCoordinateY[numberOfSearchNodes] = nodeCoordinateY[iNode];
                        searchNodeCoordinateZ[numberOfSearchNodes] = nodeCoordinateZ[iNode];

                        ++numberOfSearchNodes;
                    }
                }
            }
            PHMPI::PH_Bcast(searchNodeCoordinateX, numberOfSearchNodes * sizeof(RDouble), processorIndex);
            PHMPI::PH_Bcast(searchNodeCoordinateY, numberOfSearchNodes * sizeof(RDouble), processorIndex);
            PHMPI::PH_Bcast(searchNodeCoordinateZ, numberOfSearchNodes * sizeof(RDouble), processorIndex);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

void OversetConfigFactoryMulti::SpecifyInitialInterNodesByDistance()
{
    //TimeSpan *timeSpan = new TimeSpan();
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();

    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->SpecifyInitialInterNodesByDistance();
        }
    }
    //timeSpan->ShowTimeSpanToLogfile("SpecifyInitialInterNodesByDistance");
    //delete timeSpan;
    OutPutOversetVisualizations("1_2_BuildKeyBasicActiveNodesByDistances");
}

void OversetConfigFactoryMulti::SpecifySpecicalCellsInteractWithWallBoundary()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();

    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->SpecifySpecicalCellsInteractWithWallBoundary();
        }
    }
    OutPutOversetVisualizations("1_3_BuildInBlockCellsAndNodes");
}

void OversetConfigFactoryMulti::FillTheActiveRegion()
{
    //TimeSpan *timeSpan = new TimeSpan();
    int numberOfNewNodesPresentCycle = 0;
    int sum = 0;
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    int server = PHMPI::GetServerSepMode();

    SetNodeCharacterOnPhysicalBoundary();

    do
    {
        numberOfNewNodesPresentCycle = CharacterNodeToValidNodeByCell();

        //! mpi_reduce.
        PH_AllReduce(&numberOfNewNodesPresentCycle, &sum, 1, PH_SUM);
        numberOfNewNodesPresentCycle = sum;
        //PrintToWindow("numberOfNewNodesPresentCycle = ", numberOfNewNodesPresentCycle, "\n");

    } while (numberOfNewNodesPresentCycle > 0);
    //timeSpan->ShowTimeSpanToLogfile("FillTheActiveRegion");
    //delete timeSpan;
    FreeKeyBasicActiveNodes();

    OutPutOversetVisualizations("1_4_BuildActiveRegion");
}

void OversetConfigFactoryMulti::SetNodeCharacterOnPhysicalBoundary()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();

    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->SetNodeCharacterOnPhysicalBoundary();
        }
    }
}

void OversetConfigFactoryMulti::SetNegativeForBufferRegion()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->FreeKeyInBlockNodes();
            oversetConfigMultiZone[iZone]->SetNegativeForBufferRegion();
        }
    }
}
void OversetConfigFactoryMulti::BuildInterpolateRelationships()
{
    CreatGlobalInterCells();

    FindInterpolateRelationships();

    AllReduceZonesInterpolateCells();
}

void OversetConfigFactoryMulti::FindInterpolateRelationships()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->FindInterpolateRelationship();
        }
    }
}

void OversetConfig::FindInterpolateRelationship()
{
    RDouble *nodeCoordinate = new RDouble[THREE_D];
    int *keyActiveOfCells = this->GetKeyActiveOfCells();
    RDouble eps = GlobalDataBase::GetDoubleParaFromDB("toleranceForOversetBox");

    int blockIndex = grid->GetIBlock();
    int zoneIndex = grid->GetZoneID();
    InterCells **zonesInterpolateCells = GetOversetConfigFactoryMulti()->GetZonesInterpolateCells();

    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        InterCells *&zoneInterpolateCells = zonesInterpolateCells[iZone];

        if (blockIndex != zoneInterpolateCells->GetBlockIndex())
        {
            int numberOfInterpolateCells = zoneInterpolateCells->GetNumberOfInterCells();

            RDouble *interpolateCellCenterX = zoneInterpolateCells->GetX();
            RDouble *interpolateCellCenterY = zoneInterpolateCells->GetY();
            RDouble *interpolateCellCenterZ = zoneInterpolateCells->GetZ();

            int *donorZoneIndex = zoneInterpolateCells->GetDonorZoneIndex();
            int *donorCellIndex = zoneInterpolateCells->GetDonorCellIndex();

            int *donorLevel = zoneInterpolateCells->GetKey();
            RDouble *cellMinDistance = zoneInterpolateCells->GetMinDist();

            for (int iCell = 0; iCell < numberOfInterpolateCells; ++iCell)
            {
                nodeCoordinate[0] = interpolateCellCenterX[iCell];
                nodeCoordinate[1] = interpolateCellCenterY[iCell];
                nodeCoordinate[2] = interpolateCellCenterZ[iCell];

                int donorCellIndexOfThisCell = gridKDTree->SearchNodeInKDTree(nodeCoordinate);

                if (donorCellIndexOfThisCell != -1)
                {
                    int donorLevelNow = -3;

                    int &donorLevelOfThisCell = donorLevel[iCell];
                    RDouble &cellMinDistanceOfThisCell = cellMinDistance[iCell];

                    switch (keyActiveOfCells[donorCellIndexOfThisCell])
                    {
                    case -1:
                        donorLevelNow = -1;
                        break;
                    case   0:
                        donorLevelNow = -2;
                        break;
                    case   1:
                        donorLevelNow = 1;
                        break;
                    }

                    if (donorLevelNow > donorLevelOfThisCell)
                    {
                        donorLevelOfThisCell = donorLevelNow;
                        cellMinDistanceOfThisCell = ComputeDistanceInDonorZone(nodeCoordinate, donorCellIndexOfThisCell);
                        donorZoneIndex[iCell] = zoneIndex;
                        donorCellIndex[iCell] = donorCellIndexOfThisCell;
                    }
                    else if (donorLevelNow == donorLevelOfThisCell)
                    {
                        RDouble cellMinDistanceNow = ComputeDistanceInDonorZone(nodeCoordinate, donorCellIndexOfThisCell);
                        if (cellMinDistanceNow < cellMinDistanceOfThisCell)
                        {
                            cellMinDistanceOfThisCell = cellMinDistanceNow;
                            donorZoneIndex[iCell] = zoneIndex;
                            donorCellIndex[iCell] = donorCellIndexOfThisCell;
                        }
                    }
                }
            }
        }
    }

    delete gridKDTree;
    gridKDTree = NULL;

    delete[] minBoxOfZone;
    delete[] maxBoxOfZone;

    minBoxOfZone = NULL;
    maxBoxOfZone = NULL;

    delete[] nodeCoordinate;
    nodeCoordinate = NULL;
}

void OversetConfigFactoryMulti::CreatGlobalInterCells()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int   blockIndex = 0;
        int   numberOfInterpolateCells = 0;

        int zoneProcessorIndexOfThisZone = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessorIndex == zoneProcessorIndexOfThisZone)
        {
            blockIndex = oversetConfigMultiZone[iZone]->GetGrid()->GetIBlock();
            numberOfInterpolateCells = ( oversetConfigMultiZone[iZone]->GetInterpolateCellList() ).size();
        }

        PHMPI::PH_Bcast(&blockIndex, sizeof(int), zoneProcessorIndexOfThisZone);
        PHMPI::PH_Bcast(&numberOfInterpolateCells, sizeof(int), zoneProcessorIndexOfThisZone);

        if (zonesInterpolateCells[iZone])
        {
            delete zonesInterpolateCells[iZone];
            zonesInterpolateCells[iZone] = NULL;
        }

        zonesInterpolateCells[iZone] = new InterCells(numberOfInterpolateCells, blockIndex);

        if (numberOfInterpolateCells)
        {
            RDouble *interpolateCellCenterX = zonesInterpolateCells[iZone]->GetX();
            RDouble *interpolateCellCenterY = zonesInterpolateCells[iZone]->GetY();
            RDouble *interpolateCellCenterZ = zonesInterpolateCells[iZone]->GetZ();

            if (currentProcessorIndex == zoneProcessorIndexOfThisZone)
            {
                vector< int > &interpolateCellList = oversetConfigMultiZone[iZone]->GetInterpolateCellList();

                RDouble *cellCenterX = oversetConfigMultiZone[iZone]->GetGrid()->GetCellCenterX();
                RDouble *cellCenterY = oversetConfigMultiZone[iZone]->GetGrid()->GetCellCenterY();
                RDouble *cellCenterZ = oversetConfigMultiZone[iZone]->GetGrid()->GetCellCenterZ();

                for (int iCell = 0; iCell < numberOfInterpolateCells; ++iCell)
                {
                    int &cellIndex = interpolateCellList[iCell];

                    interpolateCellCenterX[iCell] = cellCenterX[cellIndex];
                    interpolateCellCenterY[iCell] = cellCenterY[cellIndex];
                    interpolateCellCenterZ[iCell] = cellCenterZ[cellIndex];
                }
            }

            PHMPI::PH_Bcast(interpolateCellCenterX, numberOfInterpolateCells * sizeof(RDouble), zoneProcessorIndexOfThisZone);
            PHMPI::PH_Bcast(interpolateCellCenterY, numberOfInterpolateCells * sizeof(RDouble), zoneProcessorIndexOfThisZone);
            PHMPI::PH_Bcast(interpolateCellCenterZ, numberOfInterpolateCells * sizeof(RDouble), zoneProcessorIndexOfThisZone);
        }
    }
}

void OversetConfigFactoryMulti::AllReduceZonesInterpolateCells()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        InterCells *zoneInterpolateCells = zonesInterpolateCells[iZone];

        int numberOfInterpolateCells = zoneInterpolateCells->GetNumberOfInterCells();

        if (numberOfInterpolateCells)
        {
            int *donorZoneIndex = zoneInterpolateCells->GetDonorZoneIndex();
            int *donorCellIndex = zoneInterpolateCells->GetDonorCellIndex();
            int *donorLevel = zoneInterpolateCells->GetKey();
            RDouble *cellMinDistance = zoneInterpolateCells->GetMinDist();

            int *donorZoneIndexGlobal = new int[numberOfInterpolateCells];
            int *donorCellIndexGlobal = new int[numberOfInterpolateCells];
            int *donorLevelGlobal = new int[numberOfInterpolateCells];
            RDouble *cellMinDistanceGlobal = new RDouble[numberOfInterpolateCells];

            PH_AllReduce(donorLevel, donorLevelGlobal, numberOfInterpolateCells, PH_MAX);
            for (int iCell = 0; iCell < numberOfInterpolateCells; ++iCell)
            {
                if (donorLevelGlobal[iCell] > donorLevel[iCell])
                {
                    cellMinDistance[iCell] = LARGE;
                    donorZoneIndex[iCell] = -1;
                }
            }

            PH_AllReduce(cellMinDistance, cellMinDistanceGlobal, numberOfInterpolateCells, PH_MIN);
            for (int iCell = 0; iCell < numberOfInterpolateCells; ++iCell)
            {
                if (cellMinDistanceGlobal[iCell] < cellMinDistance[iCell])
                {
                    donorZoneIndex[iCell] = -1;
                }
            }

            PH_AllReduce(donorZoneIndex, donorZoneIndexGlobal, numberOfInterpolateCells, PH_MAX);
            for (int iCell = 0; iCell < numberOfInterpolateCells; ++iCell)
            {
                if (donorZoneIndexGlobal[iCell] != donorZoneIndex[iCell])
                {
                    donorCellIndex[iCell] = -1;
                }
            }

            PH_AllReduce(donorCellIndex, donorCellIndexGlobal, numberOfInterpolateCells, PH_MAX);

            SetField(donorLevel, donorLevelGlobal, numberOfInterpolateCells);
            SetField(cellMinDistance, cellMinDistanceGlobal, numberOfInterpolateCells);
            SetField(donorZoneIndex, donorZoneIndexGlobal, numberOfInterpolateCells);
            SetField(donorCellIndex, donorCellIndexGlobal, numberOfInterpolateCells);

            delete[] donorLevelGlobal;
            delete[] cellMinDistanceGlobal;
            delete[] donorZoneIndexGlobal;
            delete[] donorCellIndexGlobal;

            donorLevelGlobal = NULL;
            cellMinDistanceGlobal = NULL;
            donorZoneIndexGlobal = NULL;
            donorCellIndexGlobal = NULL;
        }
        //int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        //if (currentProcessorIndex != processorIndex)
        //{
        //    delete zoneInterpolateCells;
        //    zoneInterpolateCells = NULL;
        //}
    }
}


//! Identify active cells.
void OversetConfigFactoryMulti::FindActiveCells()
{
    //TimeSpan *timeSpan = new TimeSpan();
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    //! Interpolation cells are determined.
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (processorIndex == currentProcessorIndex)
        {
            oversetConfigMultiZone[iZone]->FindActiveCells();
        }
    }
    //timeSpan->ShowTimeSpanToLogfile("FindActiveCells");
    //delete timeSpan;
    OutPutOversetVisualizations("2_1_FindActiveCells");
}

//! Identify interpolation cells.
void OversetConfigFactoryMulti::FindInterCells()
{
    //TimeSpan *timeSpan = new TimeSpan();
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();

    //! Interpolation cells are determined.

    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int num = 0;
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (processorIndex == currentProcessorIndex)
        {
            oversetConfigMultiZone[iZone]->FindInterCells();
        }
    }
    OutPutOversetVisualizations("2_2_FindInterpolateCells");
    //timeSpan->ShowTimeSpanToLogfile("FindInterCells");
    //delete timeSpan;
}

void OversetConfigFactoryMulti::CheckInterpolateRelationships()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    for (int iBlock = 0; iBlock < numberOfBlocks; ++iBlock)
    {
        int numberOfInterpolateCellsInThisBlock;

        int numberOfNegativeCellsInThisBlock_1;
        int numberOfNegativeCellsInThisBlock_2;
        int numberOfNegativeCellsInThisBlock_3;

        int numberOfInterpolateCellsInThisBlockLocal = 0;

        int numberOfNegativeCellsInThisBlock_1Local = 0;
        int numberOfNegativeCellsInThisBlock_2Local = 0;
        int numberOfNegativeCellsInThisBlock_3Local = 0;

        for (int iZone = 0; iZone < numberOfZones; ++iZone)
        {
            int zoneProcessorIndex = PHMPI::GetZoneProcessorID(iZone);
            if (currentProcessorIndex == zoneProcessorIndex)
            {
                oversetConfigMultiZone[iZone]->CheckInterpolateRelationship(iBlock, zonesInterpolateCells[iZone],
                    numberOfInterpolateCellsInThisBlockLocal,
                    numberOfNegativeCellsInThisBlock_1Local,
                    numberOfNegativeCellsInThisBlock_2Local,
                    numberOfNegativeCellsInThisBlock_3Local);
            }
        }

        PH_AllReduce(&numberOfInterpolateCellsInThisBlockLocal, &numberOfInterpolateCellsInThisBlock, 1, PH_SUM);
        PH_AllReduce(&numberOfNegativeCellsInThisBlock_1Local, &numberOfNegativeCellsInThisBlock_1, 1, PH_SUM);
        PH_AllReduce(&numberOfNegativeCellsInThisBlock_2Local, &numberOfNegativeCellsInThisBlock_2, 1, PH_SUM);
        PH_AllReduce(&numberOfNegativeCellsInThisBlock_3Local, &numberOfNegativeCellsInThisBlock_3, 1, PH_SUM);

        if (currentProcessorIndex == PHMPI::GetServerProcessorID())
        {
            //cout << "Block" << iBlock << "->InterpolateCellsNumber = " << numberOfInterpolateCellsInThisBlock << endl;
            //
            //if (numberOfNegativeCellsInThisBlock_1)
            //{
            //    cout << "Block" << iBlock << "->NegativeCell_1 = " << numberOfNegativeCellsInThisBlock_1 << endl;
            //}
            //if (numberOfNegativeCellsInThisBlock_2)
            //{
            //    cout << "Block" << iBlock << "->NegativeCell_2 = " << numberOfNegativeCellsInThisBlock_2 << endl;
            //}
            //if (numberOfNegativeCellsInThisBlock_3)
            //{
            //    cout << "Block" << iBlock << "->NegativeCell_3 = " << numberOfNegativeCellsInThisBlock_3 << endl;
            //    exit(- 1);
            //}
        }
    }

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int zoneProcessorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessorIndex == zoneProcessorIndex)
        {
            //zonesInterpolateCells[iZone]->FreeUselessMemery();
            //if (zonesInterpolateCells[iZone])
            //{
            //    delete zonesInterpolateCells[iZone];
            //    zonesInterpolateCells[iZone] = NULL;
            //}
        }
    }
}

void OversetConfigFactoryMulti::ResetOversetInformation()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();

    //! Initialize oversetInformation on each zone.
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->InitializeOversetInformation();
        }
    }

    //! Get oversetInterCells and oversetDonorCells at MainProcess, where oversetDonorCells should swap with neighbor zone.
    GetInterAndDonorCellManagersNew();

    //! PostProcessOversetInformation, get storageIndexContainer (Interpolation or contribution cell for a neighbor zone, its serial number in the whole interpolation or contribution cell).
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->PostProcessOversetInformation();
        }
    }
}

void OversetConfigFactoryMulti::GetInterAndDonorCellManagers()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();

    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int sendProcessorIndex = PHMPI::GetZoneProcessorIDSepMode(iZone);

        //! First get neighborZoneIndexContainer and broadcast to all process.
        PHVectorInt1D neighborZoneIndexContainer;
        if (currentProcessorIndex == sendProcessorIndex)
        {
            OversetInformationProxy *oversetInformationProxy = oversetConfigMultiZone[iZone]->GetOversetInformationProxy();
            oversetInformationProxy->GetNeighborZoneIndexContainer(neighborZoneIndexContainer);
        }

        int numberOfNeighbors = static_cast<int>( neighborZoneIndexContainer.size() );
        PHMPI::PH_Bcast(&numberOfNeighbors, sizeof(int), sendProcessorIndex);

        if (numberOfNeighbors != 0)
        {
            neighborZoneIndexContainer.resize(numberOfNeighbors);
            PHMPI::PH_Bcast(&neighborZoneIndexContainer[0], numberOfNeighbors * sizeof(int), sendProcessorIndex);
        }

        //! Loop query all neighbors one by one.
        for (int iNeighbor = 0; iNeighbor < numberOfNeighbors; ++iNeighbor)
        {
            int donorZoneIndex = neighborZoneIndexContainer[iNeighbor];
            int receiveProcessorIndex = PHMPI::GetZoneProcessorIDSepMode(donorZoneIndex);

            PHVectorInt1D     donorCellIndexContainer;
            PHVectorRDouble1D interCellCenterX;
            PHVectorRDouble1D interCellCenterY;
            PHVectorRDouble1D interCellCenterZ;

            if (currentProcessorIndex == sendProcessorIndex)
            {
                //! Get interCells ready.
                OversetInformationProxy *oversetInformationProxy = oversetConfigMultiZone[iZone]->GetOversetInformationProxy();
                OversetSwapManager *oversetSwapManager = oversetInformationProxy->GetOversetSwapManager();

                OversetInterplateCell *oversetInterplateCell = new OversetInterplateCell();
                oversetInterplateCell->SetZoneIndex(iZone);
                oversetInterplateCell->SetDonorZoneIndex(donorZoneIndex);
                PHVectorInt1D &cellIndexContainer = oversetInterplateCell->GetCellIndexContainer();
                oversetInformationProxy->GetCellIndexContainer(donorZoneIndex, cellIndexContainer);

                oversetSwapManager->AddOversetInterplateCell(oversetInterplateCell);

                //! Get necessary information ready for donorCells.
                oversetInformationProxy->GetDonorCellIndexContainer(donorZoneIndex, donorCellIndexContainer, interCellCenterX, interCellCenterY, interCellCenterZ);
            }

            int numberOfDonorCells = static_cast<int>( donorCellIndexContainer.size() );
            PHMPI::PH_Trade(&numberOfDonorCells, sizeof(int), sendProcessorIndex, receiveProcessorIndex, iZone);

            int *index = new int[numberOfDonorCells];
            RDouble *xCell = new RDouble[numberOfDonorCells];
            RDouble *yCell = new RDouble[numberOfDonorCells];
            RDouble *zCell = new RDouble[numberOfDonorCells];

            donorCellIndexContainer.resize(numberOfDonorCells);
            interCellCenterX.resize(numberOfDonorCells);
            interCellCenterY.resize(numberOfDonorCells);
            interCellCenterZ.resize(numberOfDonorCells);

            if (currentProcessorIndex == sendProcessorIndex)
            {
                for (int iCell = 0; iCell < numberOfDonorCells; ++iCell)
                {
                    index[iCell] = donorCellIndexContainer[iCell];
                    xCell[iCell] = interCellCenterX[iCell];
                    yCell[iCell] = interCellCenterY[iCell];
                    zCell[iCell] = interCellCenterZ[iCell];
                }
            }

            PHMPI::PH_Trade(index, numberOfDonorCells * sizeof(int), sendProcessorIndex, receiveProcessorIndex, iZone);
            PHMPI::PH_Trade(xCell, numberOfDonorCells * sizeof(RDouble), sendProcessorIndex, receiveProcessorIndex, iZone);
            PHMPI::PH_Trade(yCell, numberOfDonorCells * sizeof(RDouble), sendProcessorIndex, receiveProcessorIndex, iZone);
            PHMPI::PH_Trade(zCell, numberOfDonorCells * sizeof(RDouble), sendProcessorIndex, receiveProcessorIndex, iZone);

            if (currentProcessorIndex == receiveProcessorIndex)
            {
                for (int iCell = 0; iCell < numberOfDonorCells; ++iCell)
                {
                    donorCellIndexContainer[iCell] = index[iCell];
                    interCellCenterX[iCell] = xCell[iCell];
                    interCellCenterY[iCell] = yCell[iCell];
                    interCellCenterZ[iCell] = zCell[iCell];
                }

                OversetInformationProxy *oversetInformationProxy = oversetConfigMultiZone[donorZoneIndex]->GetOversetInformationProxy();
                OversetSwapManager *oversetSwapManager = oversetInformationProxy->GetOversetSwapManager();

                OversetDonorCell *oversetDonorCell = new OversetDonorCell();
                oversetDonorCell->SetZoneIndex(donorZoneIndex);
                oversetDonorCell->SetInterplateZoneIndex(iZone);
                oversetDonorCell->SetDonorCellIndexContainer(donorCellIndexContainer);
                oversetDonorCell->SetInterCellCenterX(interCellCenterX);
                oversetDonorCell->SetInterCellCenterY(interCellCenterY);
                oversetDonorCell->SetInterCellCenterZ(interCellCenterZ);

                oversetSwapManager->AddOversetDonorCell(oversetDonorCell);
            }

            delete[] index;
            delete[] xCell;
            delete[] yCell;
            delete[] zCell;
        }
    }
}

void OversetConfigFactoryMulti::GetInterAndDonorCellManagersNew()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();

    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int sendProcessorIndex = PHMPI::GetZoneProcessorID(iZone);

        //! First get neighborZoneIndexContainer and broadcast to all process.
        PHVectorInt1D neighborZoneIndexContainer;

        GetNeighborZoneIndexContainerForTargetZone(iZone, neighborZoneIndexContainer);

        int numberOfNeighbors = static_cast<int>( neighborZoneIndexContainer.size() );

        //! Loop query all neighbors one by one.
        for (int iNeighbor = 0; iNeighbor < numberOfNeighbors; ++iNeighbor)
        {
            DataContainer *dataContainer = new DataContainer;

            int donorZoneIndex = neighborZoneIndexContainer[iNeighbor];
            int receiveProcessorIndex = PHMPI::GetZoneProcessorID(donorZoneIndex);

            int numberOfDonorCell;
            PHVectorInt1D     donorCellIndexContainer;
            PHVectorRDouble1D interCellCenterX;
            PHVectorRDouble1D interCellCenterY;
            PHVectorRDouble1D interCellCenterZ;

            if (currentProcessorIndex != receiveProcessorIndex && currentProcessorIndex != sendProcessorIndex)
            {
                continue;
            }

            if (currentProcessorIndex == sendProcessorIndex)
            {
                //! Get interCells ready.
                OversetInformationProxy *oversetInformationProxy = oversetConfigMultiZone[iZone]->GetOversetInformationProxy();
                OversetSwapManager *oversetSwapManager = oversetInformationProxy->GetOversetSwapManager();

                OversetInterplateCell *oversetInterplateCell = new OversetInterplateCell();
                oversetInterplateCell->SetZoneIndex(iZone);
                oversetInterplateCell->SetDonorZoneIndex(donorZoneIndex);
                PHVectorInt1D &cellIndexContainer = oversetInterplateCell->GetCellIndexContainer();
                oversetInformationProxy->GetCellIndexContainer(donorZoneIndex, cellIndexContainer);

                oversetSwapManager->AddOversetInterplateCell(oversetInterplateCell);

                //! Get necessary information ready for donorCells.
                oversetInformationProxy->GetDonorCellIndexContainer(donorZoneIndex, donorCellIndexContainer, interCellCenterX, interCellCenterY, interCellCenterZ);

                numberOfDonorCell = static_cast<int>( donorCellIndexContainer.size() );

                PHSPACE::PHWrite(dataContainer, numberOfDonorCell);

                PHSPACE::PHWrite(dataContainer, donorCellIndexContainer, numberOfDonorCell);
                PHSPACE::PHWrite(dataContainer, interCellCenterX, numberOfDonorCell);
                PHSPACE::PHWrite(dataContainer, interCellCenterY, numberOfDonorCell);
                PHSPACE::PHWrite(dataContainer, interCellCenterZ, numberOfDonorCell);
            }

            PHMPI::PH_Trade(dataContainer, sendProcessorIndex, receiveProcessorIndex);

            if (currentProcessorIndex == receiveProcessorIndex)
            {
                dataContainer->moveto_begin();

                PHSPACE::PHRead(dataContainer, numberOfDonorCell);

                donorCellIndexContainer.resize(numberOfDonorCell);
                interCellCenterX.resize(numberOfDonorCell);
                interCellCenterY.resize(numberOfDonorCell);
                interCellCenterZ.resize(numberOfDonorCell);

                PHSPACE::PHRead(dataContainer, donorCellIndexContainer, numberOfDonorCell);
                PHSPACE::PHRead(dataContainer, interCellCenterX, numberOfDonorCell);
                PHSPACE::PHRead(dataContainer, interCellCenterY, numberOfDonorCell);
                PHSPACE::PHRead(dataContainer, interCellCenterZ, numberOfDonorCell);

                OversetInformationProxy *oversetInformationProxy = oversetConfigMultiZone[donorZoneIndex]->GetOversetInformationProxy();
                OversetSwapManager *oversetSwapManager = oversetInformationProxy->GetOversetSwapManager();

                OversetDonorCell *oversetDonorCell = new OversetDonorCell();
                oversetDonorCell->SetZoneIndex(donorZoneIndex);
                oversetDonorCell->SetInterplateZoneIndex(iZone);
                oversetDonorCell->SetDonorCellIndexContainer(donorCellIndexContainer);
                oversetDonorCell->SetInterCellCenterX(interCellCenterX);
                oversetDonorCell->SetInterCellCenterY(interCellCenterY);
                oversetDonorCell->SetInterCellCenterZ(interCellCenterZ);

                oversetSwapManager->AddOversetDonorCell(oversetDonorCell);
            }

            delete dataContainer;
        }
    }
}

void OversetConfigFactoryMulti::GetNeighborZoneIndexContainerForTargetZone(int iZone, PHVectorInt1D &neighborZoneList)
{
    InterCells *zoneInterpolateCells = zonesInterpolateCells[iZone];
    int numberOfInterCells = zoneInterpolateCells->GetNumberOfInterCells();

    int *donorZone = zoneInterpolateCells->GetDonorZoneIndex();

    set< int > neighborZoneIndexSet;

    for (int iCell = 0; iCell < numberOfInterCells; ++iCell)
    {
        neighborZoneIndexSet.insert(donorZone[iCell]);
    }

    for (set< int >::iterator iter = neighborZoneIndexSet.begin(); iter != neighborZoneIndexSet.end(); ++iter)
    {
        neighborZoneList.push_back(*iter);
    }
}

void OversetConfigFactoryMulti::OutPutOversetVisualizations(const string &name)
{
    if (GetIfOutPutOversetVisualization() == 0)
    {
        return;
    }

    string fileDir = "OversetOutput/" + name;
    PHMakeDir(fileDir);
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();

    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);

        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->OutPutOversetVisualization(name);
        }
    }
}

void OversetConfigFactoryMulti::FreeKeyBasicActiveNodes()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();

    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);

        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->FreeKeyBasicActiveNode();
        }
    }
}

void OversetConfigFactoryMulti::OutPutOversetFiles()
{
    OutPutBoundaryCells();

    PrintCellInformation();

    OutPutTecplot();

    OutPutOversetInformation();
}

void OversetConfigFactoryMulti::OutPutOversetInformation()
{
    int currentProcessor = PHMPI::GetCurrentProcessorID();
    int serverProcessor = PHMPI::GetServerSepMode();

    string oversetGridFileName = GlobalDataBase::GetStrParaFromDB("oversetGridFileName");
    fstream file;
    if (currentProcessor == serverProcessor)
    {
        PrintToWindow("Output Overset Configuration Information: ", oversetGridFileName, "\n");
        OpenFile(file, oversetGridFileName, ios_base::out | ios_base::binary);
    }

    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        DataContainer *dataContainer = new DataContainer;

        int sendProcessorIndex = PHMPI::GetZoneProcessorIDSepMode(iZone);

        if (currentProcessor == sendProcessorIndex)
        {
            oversetConfigMultiZone[iZone]->GetGrid()->Encode(dataContainer, 1);
        }

        PHMPI::PH_Trade(dataContainer, sendProcessorIndex, serverProcessor);

        if (currentProcessor == serverProcessor)
        {
            dataContainer->moveto_begin();
            dataContainer->writefile(file);
        }
        delete dataContainer;
    }

    if (currentProcessor == serverProcessor)
    {
        CloseFile(file);
    }
}

void OversetConfigFactoryMulti::OutPutTecplot()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();

    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);

        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->ShowGrid();
        }
    }
}

void OversetConfigFactoryMulti::OutPutBoundaryCells()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();

    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);

        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->ShowBoundaryCells();
        }
    }
}

void OversetConfigFactoryMulti::PrintCellInformation()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    int serverProcessorIndex = PHMPI::GetServerSepMode();

    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int numberOfActiveCells = 0;
        int numberOfNegativeCells = 0;
        int numberOfInterCells = 0;

        int sendProcessorIndex = PHMPI::GetZoneProcessorIDSepMode(iZone);

        DataContainer *dataContainer = new DataContainer;
        int tag = 0;

        if (currentProcessorIndex == sendProcessorIndex)
        {
            int *keyActiveOfCells = oversetConfigMultiZone[iZone]->GetKeyActiveOfCells();

            int numberOfCells = ( oversetConfigMultiZone[iZone]->grid )->GetNTotalCell();

            for (int iCell = 0; iCell < numberOfCells; ++iCell)
            {
                if (keyActiveOfCells[iCell] == 1)
                {
                    numberOfActiveCells++;
                }
                else if (keyActiveOfCells[iCell] == 0)
                {
                    numberOfNegativeCells++;
                }
                else if (keyActiveOfCells[iCell] == -1)
                {
                    numberOfInterCells++;
                }
            }

            if (currentProcessorIndex != serverProcessorIndex)
            {
                PHWrite(dataContainer, numberOfActiveCells);
                PHWrite(dataContainer, numberOfNegativeCells);
                PHWrite(dataContainer, numberOfInterCells);
                send(dataContainer, serverProcessorIndex, tag);
            }
        }

        if (currentProcessorIndex == serverProcessorIndex)
        {
            if (currentProcessorIndex != sendProcessorIndex)
            {
                receive(dataContainer, sendProcessorIndex, tag);
                PHRead(dataContainer, numberOfActiveCells);
                PHRead(dataContainer, numberOfNegativeCells);
                PHRead(dataContainer, numberOfInterCells);
            }

            cout << "iZone = " << iZone << endl;
            cout << "       numberOfActiveCells   = " << numberOfActiveCells << endl;
            cout << "       numberOfNegativeCells = " << numberOfNegativeCells << endl;
            cout << "       numberOfInterCells    = " << numberOfInterCells << endl;
        }

        delete dataContainer;
    }
}

void OversetConfigFactoryMulti::CharacterNodeToCell()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->CharacterNodeToCell();
        }
    }
}

//! For physical boundary advance, advanced to the interpolation boundary and stop, advanced to wallCell and stop too.
int OversetConfigFactoryMulti::CharacterNodeToValidNodeByCell()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();

    int numberOfNewNodesAllZones = 0;

    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);

        if (currentProcessorIndex == processorIndex)
        {
            int numberOfNewNodes = oversetConfigMultiZone[iZone]->CharacterNodeToValidNodeByCell();
            numberOfNewNodesAllZones += numberOfNewNodes;
        }
    }

    //! Transmission of point attributes through interface communication. InnerAndInitialInterNodes needs to be constrained.
    //CommunicateNodePropertyByKeyOfInterfaceNodesInBlocks();
    CommunicateNodeCharacterOnInterface();

    SetNegativeForInnerAndInitialInterNodes();

    return numberOfNewNodesAllZones;
}

void OversetConfigFactoryMulti::CharacterCellToNode()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();

    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->CharacterCellToNode();
        }
    }

    //! Transmission of point attributes through interface communication.
    //CommunicateNodePropertyByKeyOfInterfaceNodesInBlocks();
    CommunicateNodeCharacterOnInterface();
}

void OversetConfigFactoryMulti::CommunicateCellIBlank()
{
    using namespace PHMPI;
    ActionKey *actkey = new ActionKey();
    actkey->action = COMMCELLIBLANK;
    actkey->kind = GRID_BASED;
    actkey->level = 0;

    if (IsNeedNonBlockingCommunication(actkey))
    {
        MainTaskNonBlocking(actkey);
    }
    else
    {
        MainTaskBlocking(actkey);
    }

    delete actkey;
}

//! Through the communication of virtual cell information in interface, realize the indirect communication of point attributes.
void OversetConfigFactoryMulti::CommunicateNodeCharacterOnInterface()
{
    if (numberOfBlocks == numberOfZones)
    {
        return;
    }

    int currentProcessor = PHMPI::GetCurrentProcessorID();

    //! step1: Communication cell properties.
    CommunicateCellIBlank();

    //! step2: Judge the property of points through Interface.
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessor == processorIndex)
        {
            oversetConfigMultiZone[iZone]->CharacterInterfaceToNode();
        }
    }
}

void OversetConfigFactoryMulti::CommunicateNodePropertyByKeyOfInterfaceNodesInBlocks()
{
    if (numberOfBlocks == numberOfZones)
    {
        return;
    }

    for (int iBlock = 0; iBlock < numberOfBlocks; ++ iBlock)
    {
        SetField(keyOfInterfaceNodesInBlocks[iBlock], 0, numberOfInterfaceNodesInBlocks[iBlock]);
    }

    //! step1: Every zone update keyActiveOfGlobalNodes (Summary keyActiveOfNodes).
    int currentProcessor = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessor == processorIndex)
        {
            oversetConfigMultiZone[iZone]->UpdateInterfaceNodeCharacter();
        }
    }

    //! step2: mpiAllreduce.
    for (int iBlock = 0; iBlock < numberOfBlocks; ++ iBlock)
    {
        int &numberOfInterfaceNodesInThisBlock = numberOfInterfaceNodesInBlocks[iBlock];
        int *keyOfInterfaceNodesInThisBlockGlobal = new int[numberOfInterfaceNodesInThisBlock];

        PH_AllReduce(keyOfInterfaceNodesInBlocks[iBlock], keyOfInterfaceNodesInThisBlockGlobal, numberOfInterfaceNodesInThisBlock, MPI_MAX);
        SetField(keyOfInterfaceNodesInBlocks[iBlock], keyOfInterfaceNodesInThisBlockGlobal, numberOfInterfaceNodesInThisBlock);

        delete[] keyOfInterfaceNodesInThisBlockGlobal;
        keyOfInterfaceNodesInThisBlockGlobal = NULL;
    }

    //! step3: Assign the value of after mpiAllreduce to keyActiveOfNodes.
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessor == processorIndex)
        {
            oversetConfigMultiZone[iZone]->DownInterfaceNodeCharacter();
        }
    }
}

void OversetConfigFactoryMulti::OutInterCellInformation()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();

    string fileName = AddSymbolToFileName("interCells.dat", currentProcessorIndex);
    fstream file;
    PHSPACE::OpenFile(file, fileName, ios_base::out);

    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        InterCells *zoneInterpolateCells = zonesInterpolateCells[iZone];
        int numberOfInterCells = zoneInterpolateCells->GetNumberOfInterCells();

        RDouble *minDist = zoneInterpolateCells->GetMinDist();
        int *donorZone = zoneInterpolateCells->GetDonorZoneIndex();
        int *donorCell = zoneInterpolateCells->GetDonorCellIndex();

        file << "iZone = " << iZone << endl;
        for (int iCell = 0; iCell < numberOfInterCells; ++iCell)
        {
            file << "iCell = " << iCell << ":" << endl;
            file << "        " << minDist[iCell] << endl;
            file << "        " << donorZone[iCell] << endl;
            file << "        " << donorCell[iCell] << endl;
        }
    }
    PHSPACE::CloseFile(file);
}

void OversetConfigFactoryMulti::OutKeyActive()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorID(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->OutKeyActive();
        }
    }
}

void OversetConfigFactoryMulti::OutGlobalInterCells()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    string  fileName = AddSymbolToFileName("globalInterCell.dat", currentProcessorIndex);
    fstream file;

    PHSPACE::OpenFile(file, fileName, ios_base::out);

    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        InterCells *zoneInterpolateCells = zonesInterpolateCells[iZone];
        int blockIndex = zoneInterpolateCells->GetBlockIndex();
        int numberOfInterCells = zoneInterpolateCells->GetNumberOfInterCells();

        RDouble *x = zoneInterpolateCells->GetX();
        RDouble *y = zoneInterpolateCells->GetX();
        RDouble *z = zoneInterpolateCells->GetX();

        file << "blockIndex         = " << blockIndex << endl;
        file << "numberOfInterCells = " << numberOfInterCells << endl;

        for (int iCell = 0; iCell < numberOfInterCells; ++iCell)
        {
            file << x[iCell] << "  " << y[iCell] << "  " << z[iCell] << endl;
        }
    }

    PHSPACE::CloseFile(file);
}

void OversetConfigFactoryMulti::OutGlobalNodes()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    string fileName = AddSymbolToFileName("globalNodes.dat", currentProcessorIndex);
    fstream file;

    PHSPACE::OpenFile(file, fileName, ios_base::out);

    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        ZoneSearchNodes *zoneSearchNodes = zonesSearchNodes[iZone];
        int blockIndex = zoneSearchNodes->GetBlockIndex();
        int numberOfNodes = zoneSearchNodes->GetNumberOfNodes();

        RDouble *x = zoneSearchNodes->GetX();
        RDouble *y = zoneSearchNodes->GetX();
        RDouble *z = zoneSearchNodes->GetX();
        RDouble *pMin = zoneSearchNodes->GetPMin();
        RDouble *pMax = zoneSearchNodes->GetPMax();

        file << "blockIndex    = " << blockIndex << endl;
        file << "numberOfNodes = " << numberOfNodes << endl;

        file << pMin[0] << "  " << pMin[1] << "  " << pMin[2] << endl;

        file << pMax[0] << "  " << pMax[1] << "  " << pMax[2] << endl;

        for (int iNode = 0; iNode < numberOfNodes; ++iNode)
        {
            file << x[iNode] << "  " << y[iNode] << "  " << z[iNode] << endl;
        }
    }

    PHSPACE::CloseFile(file);
}

ZoneSearchNodes::ZoneSearchNodes(int numberOfSearchNodesIn,
    int blockIndexIn) :
    numberOfNodes(numberOfSearchNodesIn),
    blockIndex(blockIndexIn)
{
    int geometricDimension = GetDim();
    pMin = new RDouble[geometricDimension];
    pMax = new RDouble[geometricDimension];

    if (numberOfNodes)
    {
        x = new RDouble[numberOfNodes];
        y = new RDouble[numberOfNodes];
        z = new RDouble[numberOfNodes];

        nodesMinDistanceToOtherBlock = new RDouble[numberOfNodes];

        SetField(nodesMinDistanceToOtherBlock, LARGE, numberOfNodes);
    }
}

ZoneSearchNodes ::~ZoneSearchNodes()
{
    if (numberOfNodes)
    {
        delete[]x;
        delete[]y;
        delete[]z;
        delete[]nodesMinDistanceToOtherBlock;

        x = NULL;
        y = NULL;
        z = NULL;
        nodesMinDistanceToOtherBlock = NULL;
    }

    delete[]pMin;
    delete[]pMax;

    pMax = NULL;
    pMin = NULL;
}

InterCells::InterCells(int numberOfInterpolateCellsIn,
    int blockIndexIn) :
    numberOfInterCells(numberOfInterpolateCellsIn),
    blockIndex(blockIndexIn)
{
    if (numberOfInterCells)
    {
        donorZoneIndex = new int[numberOfInterCells];
        donorCellIndex = new int[numberOfInterCells];

        donorLevel = new int[numberOfInterCells];
        cellMinDistance = new RDouble[numberOfInterCells];

        x = new RDouble[numberOfInterCells];
        y = new RDouble[numberOfInterCells];
        z = new RDouble[numberOfInterCells];

        SetField(donorZoneIndex, -1, numberOfInterCells);
        SetField(donorCellIndex, -1, numberOfInterCells);

        SetField(donorLevel, -3, numberOfInterCells);
        SetField(cellMinDistance, LARGE, numberOfInterCells);
    }
}

InterCells::~InterCells()
{
    if (numberOfInterCells)
    {
        delete[] donorZoneIndex;
        delete[] donorCellIndex;

        donorZoneIndex = NULL;
        donorCellIndex = NULL;

        if (NULL != donorLevel)
        {
            delete[] donorLevel;
            delete[] cellMinDistance;

            delete[] x;
            delete[] y;
            delete[] z;

            donorLevel = NULL;
            cellMinDistance = NULL;

            x = NULL;
            y = NULL;
            z = NULL;
        }
    }
}

void InterCells::FreeUselessMemery()
{
    if (numberOfInterCells)
    {
        delete[] donorLevel;
        delete[] cellMinDistance;

        delete[] x;
        delete[] y;
        delete[] z;

        donorLevel = NULL;
        cellMinDistance = NULL;

        x = NULL;
        y = NULL;
        z = NULL;
    }
}

OversetConfig::OversetConfig()
{
    gridKDTree = NULL;

    minBoxOfZone = NULL;
    maxBoxOfZone = NULL;

    keyActiveOfNodes = NULL;

    keyBufferOfNodes = NULL;
    keyBufferOfCells = NULL;
    keySearchOfNodes = NULL;
    keySearchOfCells = NULL;
    keyInBlockNodes = NULL;
    keyBasicActiveNodes = NULL;

    nodesMinDistanceInThisBlock = NULL;
    nodesMinDistanceInOtherBlock = NULL;

    localToGlobalInterfaceNodesMap = NULL;

    zoneIndex = 0;
}

void OversetConfig::BuildBoundaryNodeList(vector< RDouble > &wallNodeCoordinateXIn,
    vector< RDouble > &wallNodeCoordinateYIn,
    vector< RDouble > &wallNodeCoordinateZIn)
{
    int   numberOfBoundaryElements = grid->GetNBoundFace();
    int *faceNodeNumber = grid->GetNodeNumberOfEachFace();
    int **faceNodeIndex = grid->GetFace2NodeArray();

    UnstructBCSet **unstructBCSet = grid->GetBCRecord();

    RDouble *nodeCoordinateX = grid->GetX();
    RDouble *nodeCoordinateY = grid->GetY();
    RDouble *nodeCoordinateZ = grid->GetZ();

    set< int > interfaceBoundaryNodeSet;
    set< int > wallBoundaryNodeSet;
    set< int > oversetBoundaryNodeSet;
    set< int > symmetryBoundaryNodeSet;
    set< int > physicalBoundaryNodeSet;

    for (int iFace = 0; iFace < numberOfBoundaryElements; ++iFace)
    {
        int boundaryConditionOfThisFace = unstructBCSet[iFace]->GetKey();
        int &numberOfNodesInThisFace = faceNodeNumber[iFace];

        switch (boundaryConditionOfThisFace)
        {
        case PHENGLEI::INTERFACE:
            for (int iNode = 0; iNode < numberOfNodesInThisFace; ++iNode)
            {
                interfaceBoundaryNodeSet.insert(faceNodeIndex[iFace][iNode]);
            }
            break;
        case PHENGLEI::SOLID_SURFACE:
            for (int iNode = 0; iNode < numberOfNodesInThisFace; ++iNode)
            {
                wallBoundaryNodeSet.insert(faceNodeIndex[iFace][iNode]);
            }
            break;
        case PHENGLEI::OVERSET:
            for (int iNode = 0; iNode < numberOfNodesInThisFace; ++iNode)
            {
                oversetBoundaryNodeSet.insert(faceNodeIndex[iFace][iNode]);
            }
            break;
        case PHENGLEI::SYMMETRY:
            for (int iNode = 0; iNode < numberOfNodesInThisFace; ++iNode)
            {
                symmetryBoundaryNodeSet.insert(faceNodeIndex[iFace][iNode]);
            }
            break;
        default:
            for (int iNode = 0; iNode < numberOfNodesInThisFace; ++iNode)
            {
                physicalBoundaryNodeSet.insert(faceNodeIndex[iFace][iNode]);
            }
        }
    }

    for (set< int >::iterator iter = interfaceBoundaryNodeSet.begin(); iter != interfaceBoundaryNodeSet.end(); ++iter)
    {
        interfaceBoundaryNodeList.push_back(*iter);
    }
    for (set< int >::iterator iter = oversetBoundaryNodeSet.begin(); iter != oversetBoundaryNodeSet.end(); ++iter)
    {
        innerBoundaryNodeList.push_back(*iter);
    }
    for (set< int >::iterator iter = physicalBoundaryNodeSet.begin(); iter != physicalBoundaryNodeSet.end(); ++iter)
    {
        physicalBoundaryNodeList.push_back(*iter);
    }

    for (set< int >::iterator iter = wallBoundaryNodeSet.begin(); iter != wallBoundaryNodeSet.end(); ++iter)
    {
        physicalBoundaryNodeList.push_back(*iter);
        wallNodeCoordinateXIn.push_back(nodeCoordinateX[*iter]);
        wallNodeCoordinateYIn.push_back(nodeCoordinateY[*iter]);
        wallNodeCoordinateZIn.push_back(nodeCoordinateZ[*iter]);
    }

    if (GetSymetryOrNot())
    {
        for (set< int >::iterator iter = symmetryBoundaryNodeSet.begin(); iter != symmetryBoundaryNodeSet.end(); ++iter)
        {
            physicalBoundaryNodeList.push_back(*iter);
        }
    }
}

void OversetConfig::BuildBoundaryNodeList()
{
    int   numberOfBoundaryElements = grid->GetNBoundFace();
    int *faceNodeNumber = grid->GetNodeNumberOfEachFace();
    int **faceNodeIndex = grid->GetFace2NodeArray();

    UnstructBCSet **unstructBCSet = grid->GetBCRecord();

    set< int > interfaceBoundaryNodeSet;
    set< int > wallBoundaryNodeSet;
    set< int > oversetBoundaryNodeSet;
    set< int > symmetryBoundaryNodeSet;
    set< int > physicalBoundaryNodeSet;

    for (int iFace = 0; iFace < numberOfBoundaryElements; ++iFace)
    {
        int boundaryConditionOfThisFace = unstructBCSet[iFace]->GetKey();
        int &numberOfNodesInThisFace = faceNodeNumber[iFace];

        switch (boundaryConditionOfThisFace)
        {
        case PHENGLEI::INTERFACE:
            for (int iNode = 0; iNode < numberOfNodesInThisFace; ++iNode)
            {
                interfaceBoundaryNodeSet.insert(faceNodeIndex[iFace][iNode]);
            }
            break;
        case PHENGLEI::SOLID_SURFACE:
            for (int iNode = 0; iNode < numberOfNodesInThisFace; ++iNode)
            {
                wallBoundaryNodeSet.insert(faceNodeIndex[iFace][iNode]);
            }
            break;
        case PHENGLEI::OVERSET:
            for (int iNode = 0; iNode < numberOfNodesInThisFace; ++iNode)
            {
                oversetBoundaryNodeSet.insert(faceNodeIndex[iFace][iNode]);
            }
            break;
        case PHENGLEI::SYMMETRY:
            for (int iNode = 0; iNode < numberOfNodesInThisFace; ++iNode)
            {
                symmetryBoundaryNodeSet.insert(faceNodeIndex[iFace][iNode]);
            }
            break;
        default:
            for (int iNode = 0; iNode < numberOfNodesInThisFace; ++iNode)
            {
                physicalBoundaryNodeSet.insert(faceNodeIndex[iFace][iNode]);
            }
        }
    }

    for (set< int >::iterator iter = interfaceBoundaryNodeSet.begin(); iter != interfaceBoundaryNodeSet.end(); ++iter)
    {
        interfaceBoundaryNodeList.push_back(*iter);
    }
    for (set< int >::iterator iter = oversetBoundaryNodeSet.begin(); iter != oversetBoundaryNodeSet.end(); ++iter)
    {
        innerBoundaryNodeList.push_back(*iter);
    }
    for (set< int >::iterator iter = physicalBoundaryNodeSet.begin(); iter != physicalBoundaryNodeSet.end(); ++iter)
    {
        physicalBoundaryNodeList.push_back(*iter);
    }

    for (set< int >::iterator iter = wallBoundaryNodeSet.begin(); iter != wallBoundaryNodeSet.end(); ++iter)
    {
        physicalBoundaryNodeList.push_back(*iter);
    }

    if (GetSymetryOrNot())
    {
        for (set< int >::iterator iter = symmetryBoundaryNodeSet.begin(); iter != symmetryBoundaryNodeSet.end(); ++iter)
        {
            physicalBoundaryNodeList.push_back(*iter);
        }
    }
    if (NULL == localToGlobalInterfaceNodesMap)
    {
        localToGlobalInterfaceNodesMap = new int[interfaceBoundaryNodeList.size()];
    }

}

void OversetConfig::BuildInBlockCellsAndNodesByInnerAuxiliaryGridKDTree(vector< GridKDTree * > &innerAuxiliaryGridKDTreesIn)
{
    RDouble *nodeCoordinate = new RDouble[THREE_D];

    int blockIndex = grid->GetIBlock();

    int numberOfNodes = grid->GetNTotalNode();
    int numberOfCells = grid->GetNTotalCell();
    int numberOfBlocks = GetOversetConfigFactoryMulti()->GetNumberOfBlocks();

    RDouble *cellCenterX = grid->GetCellCenterX();
    RDouble *cellCenterY = grid->GetCellCenterY();
    RDouble *cellCenterZ = grid->GetCellCenterZ();

    int **cellNodeIndexContainer = grid->GetCell2NodeArray();
    int *cellNodeNumberContainer = grid->GetNodeNumberOfEachCell();

    keyInBlockNodes = new int[numberOfNodes];
    SetField(keyInBlockNodes, 0, numberOfNodes);

    double tolerance = GetToleranceForOversetBox();
    int donorCellIndex;
    int symetryOrNot = GetSymetryOrNot();
    for (int iCell = 0; iCell < numberOfCells; ++iCell)
    {
        donorCellIndex = -2;

        if (!symetryOrNot || cellCenterZ[iCell] > -tolerance)
        {
            donorCellIndex = -1;

            nodeCoordinate[0] = cellCenterX[iCell];
            nodeCoordinate[1] = cellCenterY[iCell];
            nodeCoordinate[2] = cellCenterZ[iCell];

            for (int iBlock = 0; iBlock < numberOfBlocks; ++iBlock)
            {
                if (blockIndex != iBlock)
                {
                    if (innerAuxiliaryGridKDTreesIn[iBlock] != NULL)
                    {
                        donorCellIndex = innerAuxiliaryGridKDTreesIn[iBlock]->SearchNodeInKDTree(nodeCoordinate);
                    }
                }
            }
        }

        if (donorCellIndex != -1)
        {
            for (int iNode = 0; iNode < cellNodeNumberContainer[iCell]; ++iNode)
            {
                keyInBlockNodes[cellNodeIndexContainer[iCell][iNode]] = 1;
            }
        }
    }

    delete[] nodeCoordinate;
    nodeCoordinate = NULL;
}

OversetConfig::~OversetConfig()
{
    delete[] keyActiveOfNodes;
    keyActiveOfNodes = NULL;

    delete[] keyBufferOfNodes;
    keyBufferOfNodes = NULL;

    delete[] keyBufferOfCells;
    keyBufferOfCells = NULL;

    if (nodesMinDistanceInThisBlock != NULL)
    {
        delete[] nodesMinDistanceInThisBlock;
        nodesMinDistanceInThisBlock = NULL;
    }

    if (localToGlobalInterfaceNodesMap != NULL)
    {
        delete[] localToGlobalInterfaceNodesMap;
        localToGlobalInterfaceNodesMap = NULL;
    }
}

void OversetConfig::FreeKeyBasicActiveNode()
{
    delete[] keyBasicActiveNodes;
    keyBasicActiveNodes = NULL;
}

void OversetConfig::FreeKeyInBlockNodes()
{
    delete[] keyInBlockNodes;
    keyInBlockNodes = NULL;
}

void OversetConfig::Initialize(UnstructGrid *gridIn,
    int nZone, int iZone, int iBlock)
{
    geometricDimension = GetDim();
    numberOfZones = nZone;
    zoneIndex = iZone;
    blockIndex = iBlock;
    grid = gridIn;
    oversetInformationProxy = NULL;

    int numberOfNodes = grid->GetNTotalNode();
    int numberOfCells = grid->GetNTotalCell();
    int numberOfTotalCells = numberOfCells + grid->GetNBoundFace();

    keyActiveOfNodes = new int[numberOfNodes];
    SetField(keyActiveOfNodes, 0, numberOfNodes);

    nodesMinDistanceInThisBlock = new RDouble[numberOfNodes];
    SetField(nodesMinDistanceInThisBlock, LARGE, numberOfNodes);
}

void OversetConfig::SetInitialValue()
{
    interpolateCellList.resize(0);
}

//! Calculate the minimum box according to the cell.
void OversetConfig::BuildMinMaxBoxOfZoneBySearchRegion()
{
    int numberOfCells = grid->GetNTotalCell();

    int **cellNodeIndexContainer = grid->GetCell2NodeArray();
    int *cellNodeNumberContainer = grid->GetNodeNumberOfEachCell();

    vector< RDouble * > nodeCoordinate;
    RDouble *x = grid->GetX();
    nodeCoordinate.push_back(x);
    RDouble *y = grid->GetY();
    nodeCoordinate.push_back(y);
    RDouble *z = grid->GetZ();
    nodeCoordinate.push_back(z);

    SetField(minBoxOfZone, LARGE, geometricDimension);
    SetField(maxBoxOfZone, -LARGE, geometricDimension);

    RDouble eps = GetToleranceForOversetBox();
    for (int iDimension = 0; iDimension < geometricDimension; ++iDimension)
    {
        RDouble &minInThisDimension = minBoxOfZone[iDimension];
        RDouble &maxInThisDimension = maxBoxOfZone[iDimension];
        for (int iCell = 0; iCell < numberOfCells; ++iCell)
        {
            if (keySearchOfCells[iCell])
            {
                int cellNodeNumber = cellNodeNumberContainer[iCell];

                for (int iNode = 0; iNode < cellNodeNumber; ++iNode)
                {
                    int &nodeIndex = cellNodeIndexContainer[iCell][iNode];

                    RDouble &nodeCoordinateInThisDimension = nodeCoordinate[iDimension][nodeIndex];
                    minInThisDimension = MIN(minInThisDimension, nodeCoordinateInThisDimension);
                    maxInThisDimension = MAX(maxInThisDimension, nodeCoordinateInThisDimension);
                }
            }
        }
        minInThisDimension -= eps;
        maxInThisDimension += eps;
    }
}

void OversetConfig::SetInitialValueForBufferRigion()
{
    //! Because communicate_cell_iblank is the direct communication keyActiveOfCells.
    //! keyBuffer is determined by keyActiveOfCells array.
    int *keyActiveOfCells = this->GetKeyActiveOfCells();

    int numberOfNodes = grid->GetNTotalNode();
    int numberOfCells = grid->GetNTotalCell();

    SetField(keyActiveOfNodes, 0, numberOfNodes);
    SetField(keyActiveOfCells, 0, numberOfCells);

    int numberOfInnerBoundaryNodes = innerBoundaryNodeList.size();
    for (int iNode = 0; iNode < numberOfInnerBoundaryNodes; ++iNode)
    {
        keyActiveOfNodes[innerBoundaryNodeList[iNode]] = 1;
    }
}

void OversetConfig::BuildNodesMinDistanceInThisBlock(vector< RDouble > *wallNodeCoordinateXInBlocksIn,
    vector< RDouble > *wallNodeCoordinateYInBlocksIn,
    vector< RDouble > *wallNodeCoordinateZInBlocksIn)
{
    int blockIndex = grid->GetIBlock();
    int numberOfNodes = grid->GetNTotalNode();
    int numberOfWallNodesInBlock = wallNodeCoordinateXInBlocksIn[blockIndex].size();

    nodesMinDistanceInThisBlock = new RDouble[numberOfNodes];

    if (!numberOfWallNodesInBlock)
    {
        double wallDistanceOfBackground = GlobalDataBase::GetDoubleParaFromDB("walldistMainZone");
        SetField(nodesMinDistanceInThisBlock, wallDistanceOfBackground, numberOfNodes);
        return;
    }

    RDouble *nodeCoordinateX = grid->GetX();
    RDouble *nodeCoordinateY = grid->GetY();
    RDouble *nodeCoordinateZ = grid->GetZ();

    SetField(nodesMinDistanceInThisBlock, LARGE, numberOfNodes);
    for (int iNode = 0; iNode < numberOfNodes; ++iNode)
    {
        for (int iWallNode = 0; iWallNode < numberOfWallNodesInBlock; ++iWallNode)
        {
            RDouble dx = nodeCoordinateX[iNode] - wallNodeCoordinateXInBlocksIn[blockIndex][iWallNode];
            RDouble dy = nodeCoordinateY[iNode] - wallNodeCoordinateYInBlocksIn[blockIndex][iWallNode];
            RDouble dz = nodeCoordinateZ[iNode] - wallNodeCoordinateZInBlocksIn[blockIndex][iWallNode];
            RDouble tempDistance = DISTANCE(dx, dy, dz);

            nodesMinDistanceInThisBlock[iNode] = MIN(tempDistance, nodesMinDistanceInThisBlock[iNode]);
        }
    }
}

void OversetConfig::SetNodesMinDistanceInThisBlocksByWalldist()
{
    int numberOfNodes = grid->GetNTotalNode();
    vector<RDouble> &nodeWalldist = grid->GetNodeWallDist();

    SetField(nodesMinDistanceInThisBlock, &nodeWalldist[0], numberOfNodes);
}

void OversetConfig::SetNodesMinDistanceInOtherBlocks()
{
    int numberOfNodes = grid->GetNTotalNode();
    vector<RDouble> &nodeWalldistInOtherBlocks = grid->GetNodeWallDistInOtherBlock();
    nodesMinDistanceInOtherBlock = new RDouble[numberOfNodes];
    SetField(nodesMinDistanceInOtherBlock, &nodeWalldistInOtherBlocks[0], numberOfNodes);
}

void OversetConfig::SetKeyBufferByKeyActive()
{
    //! Because communicate_cell_iblank is the direct communication keyActiveOfCells.
    //! keyBuffer is determined by keyActiveOfCells array.
    int numberOfNodes = grid->GetNTotalNode();
    int numberOfCells = grid->GetNTotalCell();
    keyBufferOfNodes = new int[numberOfNodes];
    keyBufferOfCells = new int[numberOfCells];

    int *keyActiveOfCells = this->GetKeyActiveOfCells();

    SetField(keyBufferOfCells, keyActiveOfCells, numberOfCells);
    SetField(keyBufferOfNodes, keyActiveOfNodes, numberOfNodes);
}

void OversetConfig::ComputeNodeAndCellTopology()
{
    int **cellNodeIndexContainer = grid->GetCell2NodeArray();
    int **cellFaceIndexContainer = grid->GetCell2Face();
}

void OversetConfig::CopyMinDist(RDouble *nodesMinDistanceInOtherBlockIn)
{
    int numberOfNodes = grid->GetNTotalNode();
    nodesMinDistanceInOtherBlock = new RDouble[numberOfNodes];
    SetField(nodesMinDistanceInOtherBlock, LARGE, numberOfNodes);

    int iSearchNode = 0;

    for (int iNode = 0; iNode < numberOfNodes; ++iNode)
    {
        if (keySearchOfNodes[iNode])
        {
            nodesMinDistanceInOtherBlock[iNode] = nodesMinDistanceInOtherBlockIn[iSearchNode++];
        }
    }

    delete[] keySearchOfNodes;
    keySearchOfNodes = NULL;
}

void OversetConfig::BuildMinMaxBoxOfZone()
{
    minBoxOfZone = new RDouble[geometricDimension];
    maxBoxOfZone = new RDouble[geometricDimension];
    SetField(minBoxOfZone, LARGE, geometricDimension);
    SetField(maxBoxOfZone, -LARGE, geometricDimension);

    vector< RDouble * > nodeCoordinate;
    RDouble *x = grid->GetX();
    nodeCoordinate.push_back(x);
    RDouble *y = grid->GetY();
    nodeCoordinate.push_back(y);
    RDouble *z = grid->GetZ();
    nodeCoordinate.push_back(z);

    int numberOfNodes = grid->GetNTotalNode();
    for (int iDimension = 0; iDimension < geometricDimension; ++iDimension)
    {
        RDouble &minInThisDimension = minBoxOfZone[iDimension];
        RDouble &maxInThisDimension = maxBoxOfZone[iDimension];
        for (int iNode = 0; iNode < numberOfNodes; ++iNode)
        {
            RDouble &nodeCoordinateInThisDimension = nodeCoordinate[iDimension][iNode];
            minInThisDimension = MIN(minInThisDimension, nodeCoordinateInThisDimension);
            maxInThisDimension = MAX(maxInThisDimension, nodeCoordinateInThisDimension);
        }
    }
}

void OversetConfig::BuildKDTreeByMinMaxBoxOfZone()
{
    gridKDTree = new GridKDTree(grid, keySearchOfCells, minBoxOfZone, maxBoxOfZone);

    delete[] keySearchOfCells;
    keySearchOfCells = NULL;
}

RDouble OversetConfig::ComputeDistanceInDonorZone(RDouble *nodeCoordinateIn, int donorCell)
{
    if (donorCell < 0)
    {
        return LARGE;
    }

    int **cellNodeIndexContainer = grid->GetCell2NodeArray();
    int *cellNodeNumberContainer = grid->GetNodeNumberOfEachCell();

    RDouble *coordinateX = grid->GetX();
    RDouble *coordinateY = grid->GetY();
    RDouble *coordinateZ = grid->GetZ();

    int &numberOfNodeInCell = cellNodeNumberContainer[donorCell];

    RDouble weightDist = 0.0;
    RDouble weight = 0.0;

    for (int iNode = 0; iNode < numberOfNodeInCell; ++iNode)
    {
        int &nodeIndex = cellNodeIndexContainer[donorCell][iNode];

        RDouble dx = nodeCoordinateIn[0] - coordinateX[nodeIndex];
        RDouble dy = nodeCoordinateIn[1] - coordinateY[nodeIndex];
        RDouble dz = nodeCoordinateIn[2] - coordinateZ[nodeIndex];
        RDouble dist = DISTANCE(dx, dy, dz) + TINY;

        weight += 1.0 / dist;
        weightDist += 1.0 / dist * nodesMinDistanceInThisBlock[nodeIndex];
    }

    weightDist /= ( weight + TINY );

    return weightDist;
}

//! Interpolation boundary is determined according to the local distance and neighbor distance.
//! Overset boundary can be ignored.
void OversetConfig::SpecifyInitialInterNodesByDistance()
{
    int numberOfNodes = grid->GetNTotalNode();

    keyBasicActiveNodes = new int[numberOfNodes];
    SetField(keyBasicActiveNodes, 0, numberOfNodes);

    for (int iNode = 0; iNode < numberOfNodes; ++iNode)
    {
        if (nodesMinDistanceInOtherBlock[iNode] >= nodesMinDistanceInThisBlock[iNode])
        {
            keyBasicActiveNodes[iNode] = 1;
        }
    }

    delete[] nodesMinDistanceInOtherBlock;
    nodesMinDistanceInOtherBlock = NULL;
}

//! Any two points of an cell pass through the physical plane, the cell is a boundary cell.
//! This method isn't applicable in the presence of background grid.
void OversetConfig::SpecifySpecicalCellsInteractWithWallBoundary()
{
    int numberOfCells = grid->GetNTotalCell();
    int numberOfNodes = grid->GetNTotalNode();

    int **cellNodeIndexContainer = grid->GetCell2NodeArray();
    int *cellNodeNumberContainer = grid->GetNodeNumberOfEachCell();

    int *innerCellOrNot = new int[numberOfCells];
    keyInBlockNodes = new int[numberOfNodes];

    SetField(innerCellOrNot, 0, numberOfCells);
    SetField(keyInBlockNodes, 0, numberOfNodes);

    RDouble *xNode = grid->GetX();
    RDouble *yNode = grid->GetY();
    RDouble *zNode = grid->GetZ();

    for (int iCell = 0; iCell < numberOfCells; ++iCell)
    {
        int node1 = cellNodeIndexContainer[iCell][0];

        int numberOfNodeInCell = cellNodeNumberContainer[iCell];

        for (int iNode = 1; iNode < numberOfNodeInCell; ++iNode)
        {
            int node2 = cellNodeIndexContainer[iCell][iNode];

            RDouble dist1 = nodesMinDistanceInOtherBlock[node1];
            RDouble dist2 = nodesMinDistanceInOtherBlock[node2];

            RDouble dx = xNode[node1] - xNode[node2];
            RDouble dy = yNode[node1] - yNode[node2];
            RDouble dz = zNode[node1] - zNode[node2];

            RDouble length = DISTANCE(dx, dy, dz);

            if (dist1 < length || dist2 < length)
            {
                if (dist2 > 0.5 * LARGE || dist1 > 0.5 * LARGE)
                {
                    innerCellOrNot[iCell] = 1;
                    break;
                }
            }
        }
    }

    for (int iCell = 0; iCell < numberOfCells; ++iCell)
    {
        if (innerCellOrNot[iCell] == 1)
        {
            for (int iNode = 0; iNode < cellNodeNumberContainer[iCell]; ++iNode)
            {
                int nodeIndex = cellNodeIndexContainer[iCell][iNode];
                keyInBlockNodes[nodeIndex] = 1;
            }
        }
    }
    delete[]innerCellOrNot;
}

void OversetConfig::SetNegativeForBufferRegion()
{
    int numberOfNodes = grid->GetNTotalNode();
    for (int iNode = 0; iNode < numberOfNodes; ++iNode)
    {
        if (keyBufferOfNodes[iNode])
        {
            keyActiveOfNodes[iNode] = 0;
        }
    }
}

void OversetConfig::SetNegativeForInnerNodes()
{
    int numberOfNodes = grid->GetNTotalNode();

    for (int iNode = 0; iNode < numberOfNodes; ++iNode)
    {
        if (keyInBlockNodes[iNode])
        {
            keyActiveOfNodes[iNode] = 0;
        }
    }
}

void OversetConfig::SetNegativeForInnerAndInitialInterNodes()
{
    int numberOfNodes = grid->GetNTotalNode();

    for (int iNode = 0; iNode < numberOfNodes; ++iNode)
    {
        if (!keyBasicActiveNodes[iNode] || keyInBlockNodes[iNode])
        {
            keyActiveOfNodes[iNode] = 0;
        }
    }
}

void OversetConfig::SetNodeCharacterOnPhysicalBoundary()
{
    int *keyActiveOfCells = this->GetKeyActiveOfCells();

    int numberOfNodes = grid->GetNTotalNode();
    int numberOfCells = grid->GetNTotalCell();

    SetField(keyActiveOfNodes, 0, numberOfNodes);
    SetField(keyActiveOfCells, 0, numberOfCells);

    int numberOfPhysicalBoundaryNodes = physicalBoundaryNodeList.size();
    for (int iNode = 0; iNode < numberOfPhysicalBoundaryNodes; ++iNode)
    {
        keyActiveOfNodes[physicalBoundaryNodeList[iNode]] = 1;
    }
}

void OversetConfig::FillTheRemainingNodes()
{
    int numberOfNodes = grid->GetNTotalNode();

    for (int iNode = 0; iNode < numberOfNodes; ++iNode)
    {
        if (keyActiveOfNodes[iNode] == 0)
        {
            keyActiveOfNodes[iNode] = -1;
        }
    }
}

void OversetConfig::SetActivePropertyForSearchRegion()
{
    //! Because communicate_cell_iblank is the direct communication keyActiveOfCells.
    //! keyBuffer is determined by keyActiveOfCells array.
    OversetConfigFactoryMulti *overs = GetOversetConfigFactoryMulti();
    RDouble **pMinOfBlocks = overs->GetMinBoxOfBlocks();
    RDouble **pMaxOfBlocks = overs->GetMaxBoxOfBlocks();

    int numberOfNodes = grid->GetNTotalNode();
    int numberOfCells = grid->GetNTotalCell();
    int *keyActiveOfCells = GetKeyActiveOfCells();

    SetField(keyActiveOfNodes, 0, numberOfNodes);
    SetField(keyActiveOfCells, 0, numberOfCells);

    vector< RDouble * > nodeCoordinate;
    RDouble *x = grid->GetX();
    nodeCoordinate.push_back(x);
    RDouble *y = grid->GetY();
    nodeCoordinate.push_back(y);
    RDouble *z = grid->GetZ();
    nodeCoordinate.push_back(z);

    int blockIndex = grid->GetIBlock();
    int numberOfBlocks = overs->GetNumberOfBlocks();
    for (int iBlock = 0; iBlock < numberOfBlocks; ++iBlock)
    {
        if (iBlock != blockIndex)
        {
            for (int iNode = 0; iNode < numberOfNodes; ++iNode)
            {
                if (IfNodeIsInBox(nodeCoordinate, iNode, pMinOfBlocks[iBlock], pMaxOfBlocks[iBlock], geometricDimension))
                {
                    keyActiveOfNodes[iNode] = 1;
                }
            }
        }
    }
}

void OversetConfig::SetSearchRegionByActiveProperty()
{
    //! keySearch is determined by keyActive array.
    int numberOfNodes = grid->GetNTotalNode();
    int numberOfCells = grid->GetNTotalCell();

    keySearchOfNodes = new int[numberOfNodes];
    keySearchOfCells = new int[numberOfCells];

    int *keyActiveOfCells = this->GetKeyActiveOfCells();

    SetField(keySearchOfNodes, keyActiveOfNodes, numberOfNodes);
    SetField(keySearchOfCells, keyActiveOfCells, numberOfCells);
}

//! It is an active cell if the all cell nodes are 1.
void OversetConfig::FindActiveCells()
{
    int numberOfCells = grid->GetNTotalCell();

    int **cellNodeIndex = grid->GetCell2NodeArray();
    int *cellNodeNumber = grid->GetNodeNumberOfEachCell();

    int *keyActiveOfCells = this->GetKeyActiveOfCells();

    SetField(keyActiveOfCells, 1, numberOfCells);

    for (int iCell = 0; iCell < numberOfCells; ++iCell)
    {
        int &numberOfNodeInCell = cellNodeNumber[iCell];
        for (int iNode = 0; iNode < numberOfNodeInCell; ++iNode)
        {
            if (!keyActiveOfNodes[cellNodeIndex[iCell][iNode]])
            {
                keyActiveOfCells[iCell] = 0;
                break;
            }
        }
    }
}

void OversetConfig::SearchZonesSearchNodesInKDTree()
{
    RDouble *nodeCoordinate = new RDouble[THREE_D];
    int blockIndex = grid->GetIBlock();

    ZoneSearchNodes **zonesSearchNodes = GetOversetConfigFactoryMulti()->GetZonesSearchNodes();
    for (int iZone = 0; iZone < numberOfZones; ++iZone)
    {
        ZoneSearchNodes *&zoneSearchNodes = zonesSearchNodes[iZone];
        RDouble *pMin = zoneSearchNodes->GetPMin();
        RDouble *pMax = zoneSearchNodes->GetPMax();

        if (blockIndex != zoneSearchNodes->GetBlockIndex())
        {
            if (IfZoneSearchNodesIntersectionWithThisZone(pMin, pMax))
            {
                int numberOfSearchNodes = zoneSearchNodes->GetNumberOfNodes();

                RDouble *searchNodeCoordinateX = zoneSearchNodes->GetX();
                RDouble *searchNodeCoordinateY = zoneSearchNodes->GetY();
                RDouble *searchNodeCoordinateZ = zoneSearchNodes->GetZ();

                RDouble *nodesMinDistanceToOtherBlock = zoneSearchNodes->GetNodesMinDistanceToOtherBlock();

                for (int iNode = 0; iNode < numberOfSearchNodes; ++iNode)
                {
                    RDouble &nodesMinDistanceInOtherBlockOfThisNode = nodesMinDistanceToOtherBlock[iNode];

                    nodeCoordinate[0] = searchNodeCoordinateX[iNode];
                    nodeCoordinate[1] = searchNodeCoordinateY[iNode];
                    nodeCoordinate[2] = searchNodeCoordinateZ[iNode];

                    nodesMinDistanceInOtherBlockOfThisNode = MIN(nodesMinDistanceInOtherBlockOfThisNode,
                        ComputeDistanceInDonorZone(nodeCoordinate, FindDonorCell(nodeCoordinate)));
                }
            }
        }
    }

    delete[] nodeCoordinate;
    nodeCoordinate = NULL;
}

int OversetConfig::FindDonorCell(RDouble *nodeCoordinateIn)
{
    int donorCellIndex = gridKDTree->SearchNodeInKDTree(nodeCoordinateIn);
    if (donorCellIndex != -1)
    {
        if (keyBufferOfCells[donorCellIndex]) donorCellIndex = -1;
    }

    return donorCellIndex;
}

//! It is an interpolation cell if the cell node has both 1 and 0.
void OversetConfig::FindInterCells()
{
    int numberOfCells = grid->GetNTotalCell();

    int **cellNodeIndexContainer = grid->GetCell2NodeArray();
    int *cellNodeNumberContainer = grid->GetNodeNumberOfEachCell();

    int *keyActiveOfCells = this->GetKeyActiveOfCells();

    for (int iCell = 0; iCell < numberOfCells; ++iCell)
    {
        if (!keyActiveOfCells[iCell])
        {
            int numberOfNodeInCell = cellNodeNumberContainer[iCell];
            for (int iNode = 0; iNode < numberOfNodeInCell; ++iNode)
            {
                int nodeIndex = cellNodeIndexContainer[iCell][iNode];
                if (keyActiveOfNodes[nodeIndex])
                {
                    keyActiveOfCells[iCell] = -1;
                    interpolateCellList.push_back(iCell);
                    break;
                }
            }
        }
    }
}

void OversetConfig::InitializeOversetInformation()
{
    //! Movement of the mesh, because of the ressembly of the mesh, it is necessary to delete oversetInformationProxy and regenerate it.
    if (oversetInformationProxy != 0)
    {
        delete oversetInformationProxy;
        oversetInformationProxy = 0;
    }

    oversetInformationProxy = new OversetInformationProxy();
    oversetInformationProxy->SetGrid(grid);
    grid->SetOversetInformationProxy(oversetInformationProxy);

    //
    PHVector1D < OversetCell * > *oversetCellContainer = new PHVector1D < OversetCell * >();

    InterCells **zonesInterpolateCells = GetOversetConfigFactoryMulti()->GetZonesInterpolateCells();
    InterCells *zoneInterpolateCells = zonesInterpolateCells[zoneIndex];

    int numberOfInterCells = zoneInterpolateCells->GetNumberOfInterCells();

    int *donorCellIndex = zoneInterpolateCells->GetDonorCellIndex();
    int *donorZoneIndex = zoneInterpolateCells->GetDonorZoneIndex();

    for (int iCell = 0; iCell < numberOfInterCells; ++iCell)
    {
        int interCell = interpolateCellList[iCell];
        int donorCell = donorCellIndex[iCell];
        int donorZone = donorZoneIndex[iCell];

        OversetCell *oversetCell = new OversetCell();
        oversetCell->SetCellIndex(interCell);
        oversetCell->SetZoneIndex(zoneIndex);
        oversetCell->SetDonorCellIndex(donorCell);
        oversetCell->SetDonorZoneIndex(donorZone);

        oversetCellContainer->push_back(oversetCell);
    }

    oversetInformationProxy->InitializeReceivingInformation(oversetCellContainer);
}

void OversetConfig::PostProcessOversetInformation()
{
    oversetInformationProxy->PostProcessOversetInformation();
}

void OversetConfig::MarkGlobalInterfaceNodes(int iBlock, vector< int > &keyNodesGlobal)
{
    if (blockIndex != iBlock) return;

    int numberOfInterNodes = interfaceBoundaryNodeList.size();
    for (int iNode = 0; iNode < numberOfInterNodes; ++iNode)
    {
        int localNodeIndex = interfaceBoundaryNodeList[iNode];
        int globalNodeIndex = localToGlobalNodesIndex[localNodeIndex];
        keyNodesGlobal[globalNodeIndex] = 1;
    }
}

void OversetConfig::CalLocalToGlobalInterfaceNodes(int iBlock, vector< int > &keyNodesGlobal)
{
    if (blockIndex != iBlock) return;

    int numberOfInterNodes = interfaceBoundaryNodeList.size();
    localToGlobalInterNodesIndex.resize(numberOfInterNodes);

    for (int iNode = 0; iNode < numberOfInterNodes; ++iNode)
    {
        int localNodeIndex = interfaceBoundaryNodeList[iNode];
        int globalNodeIndex = localToGlobalNodesIndex[localNodeIndex];
        int globelInterNodeIndex = keyNodesGlobal[globalNodeIndex];

        //! Connection from local interNodes to global sequence.
        localToGlobalInterNodesIndex[iNode] = globelInterNodeIndex;
    }
}

//! If a cell has a node of 1, its attributes is 1, otherwise it is 0.
void OversetConfig::CharacterNodeToCell()
{
    int numberOfCells = grid->GetNTotalCell();

    int **cellNodeIndex = grid->GetCell2NodeArray();
    int *cellNodeNumber = grid->GetNodeNumberOfEachCell();

    int *keyActiveOfCells = this->GetKeyActiveOfCells();
    //SetField(keyActiveOfCells, 0, numberOfCells);

    for (int iCell = 0; iCell < numberOfCells; ++iCell)
    {
        int &numberOfNodeInCell = cellNodeNumber[iCell];
        for (int iNode = 0; iNode < numberOfNodeInCell; ++iNode)
        {
            if (keyActiveOfNodes[cellNodeIndex[iCell][iNode]])
            {
                keyActiveOfCells[iCell] = 1;
                break;
            }
        }
    }
}

//! Property of the point on interface is judged by the property of the virtual cell on interface.
void OversetConfig::CharacterInterfaceToNode()
{
    int numberOfBoundaryFaces = grid->GetNBoundFace();

    UnstructBCSet **unstructBCSet = grid->GetBCRecord();

    int *faceNodeNumber = grid->GetNodeNumberOfEachFace();
    int *faceNodeIndex = grid->GetFace2Node();

    int *leftCellIndex = grid->GetLeftCellOfFace();
    int *rightCellIndex = grid->GetRightCellOfFace();

    int *keyActiveOfCells = grid->GetBlankIndex();

    int iCount = 0;
    for (int iFace = 0; iFace < numberOfBoundaryFaces; ++iFace)
    {
        int boundaryType = unstructBCSet[iFace]->GetKey();

        int numberOfNodesInFace = faceNodeNumber[iFace];

        int leftCell = leftCellIndex[iFace];
        int rightCell = rightCellIndex[iFace];

        if (boundaryType == PHENGLEI::INTERFACE)
        {
            if (keyActiveOfCells[leftCell] != 1 && keyActiveOfCells[rightCell] != 1)
            {
                iCount += numberOfNodesInFace;
                continue;
            }

            for (int iNode = 0; iNode < numberOfNodesInFace; ++iNode)
            {
                int nodeIndex = faceNodeIndex[iCount];
                keyActiveOfNodes[nodeIndex] = 1;
                iCount++;
            }
        }
        else
        {
            iCount += numberOfNodesInFace;
        }
    }
}

void OversetConfig::CharacterCellToNode()
{
    int *keyActiveOfCells = this->GetKeyActiveOfCells();

    int numberOfCells = grid->GetNTotalCell();
    int numberOfNodes = grid->GetNTotalNode();

    int **cellNodeIndexContainer = grid->GetCell2NodeArray();
    int *cellNodeNumberContainer = grid->GetNodeNumberOfEachCell();

    for (int iCell = 0; iCell < numberOfCells; ++iCell)
    {
        if (keyActiveOfCells[iCell])
        {
            int &numberOfNodeInCell = cellNodeNumberContainer[iCell];
            for (int iNode = 0; iNode < numberOfNodeInCell; ++iNode)
            {
                int &nodeIndex = cellNodeIndexContainer[iCell][iNode];
                keyActiveOfNodes[nodeIndex] = 1;
            }
        }
    }
}

//! Property of the point on interface is judged by the property of the virtual cell on interface.
int OversetConfig::CharacterNodeToValidNodeByCell()
{
    int numberOfCells = grid->GetNTotalCell();

    int **cellNodeIndexContainer = grid->GetCell2NodeArray();
    int *cellNodeNumberContainer = grid->GetNodeNumberOfEachCell();

    int *keyActiveOfCells = this->GetKeyActiveOfCells();

    int numberOfNewNodes = 0;

    for (int iCell = 0; iCell < numberOfCells; ++iCell)
    {
        int cellNodeNumber = cellNodeNumberContainer[iCell];

        for (int iNode = 0; iNode < cellNodeNumber; ++iNode)
        {
            int nodeIndex = cellNodeIndexContainer[iCell][iNode];

            if (keyActiveOfNodes[nodeIndex])
            {
                keyActiveOfCells[iCell] = 1;
                break;
            }
        }

        if (keyActiveOfCells[iCell])
        {
            for (int iNode = 0; iNode < cellNodeNumber; ++iNode)
            {
                int nodeIndex = cellNodeIndexContainer[iCell][iNode];
                if (keyBasicActiveNodes[nodeIndex] && !keyInBlockNodes[nodeIndex] && !keyActiveOfNodes[nodeIndex])
                {
                    keyActiveOfNodes[nodeIndex] = 1;
                    numberOfNewNodes++;
                }
            }
        }
    }
    return numberOfNewNodes;
}

void OversetConfig::CheckInterpolateRelationship(int iBlockIn,
    InterCells *zoneInterpolateCellsIn,
    int &numberOfInterpolateCellsInThisBlockIn,
    int &numberOfNegativeCellsInThisBlock_1In,
    int &numberOfNegativeCellsInThisBlock_2In,
    int &numberOfNegativeCellsInThisBlock_3In)
{
    if (iBlockIn == zoneInterpolateCellsIn->GetBlockIndex())
    {
        int        zoneIndex = grid->GetZoneID();

        int        numberOfInterpolateCells = zoneInterpolateCellsIn->GetNumberOfInterCells();

        int *donorLevel = zoneInterpolateCellsIn->GetKey();

        int *donorZoneIndex = zoneInterpolateCellsIn->GetDonorZoneIndex();
        int *donorCellIndex = zoneInterpolateCellsIn->GetDonorCellIndex();

        RDouble *interpolateCellCenterX = zoneInterpolateCellsIn->GetX();
        RDouble *interpolateCellCenterY = zoneInterpolateCellsIn->GetY();
        RDouble *interpolateCellCenterZ = zoneInterpolateCellsIn->GetZ();

        for (int iCell = 0; iCell < numberOfInterpolateCells; ++iCell)
        {
            ++numberOfInterpolateCellsInThisBlockIn;
            int &donorLevelOfThisCell = donorLevel[iCell];
            switch (donorLevelOfThisCell)
            {
            case -3:
                ++numberOfNegativeCellsInThisBlock_3In;
                cout << "ERROR:" << endl;
                cout << "DonorLevel  = " << donorLevelOfThisCell << endl;
                cout << "Interpolate = " << zoneIndex << " -> " << interpolateCellList[iCell] << endl;
                cout << "Donor       = " << donorZoneIndex[iCell] << " -> " << donorCellIndex[iCell] << endl;
                cout << "Postion     = " << setprecision(8) << setw(12)
                    << interpolateCellCenterX[iCell] << " "
                    << interpolateCellCenterY[iCell] << " "
                    << interpolateCellCenterZ[iCell] << endl;
                break;
            case -2:
                ++numberOfNegativeCellsInThisBlock_2In;
                cout << "WARNING:" << endl;
                cout << "DonorLevel  = " << donorLevelOfThisCell << endl;
                cout << "Interpolate = " << zoneIndex << "->" << interpolateCellList[iCell] << endl;
                cout << "Donor       = " << donorZoneIndex[iCell] << "->" << donorCellIndex[iCell] << endl;
                cout << "Postion     = " << setprecision(8) << setw(12)
                    << interpolateCellCenterX[iCell] << " "
                    << interpolateCellCenterY[iCell] << " "
                    << interpolateCellCenterZ[iCell] << endl;
                break;
            case -1:
                ++numberOfNegativeCellsInThisBlock_1In;
                break;
            }
        }
    }
}

void OversetConfig::DownInterfaceNodeCharacter()
{
    int    blockIndex = grid->GetIBlock();
    int **keyOfInterfaceNodesInBlocks = GetOversetConfigFactoryMulti()->GetKeyOfInterfaceNodesInBlocks();

    int numberOfInterfaceBoundaryNodes = interfaceBoundaryNodeList.size();
    for (int iNode = 0; iNode < numberOfInterfaceBoundaryNodes; ++ iNode)
    {
        keyActiveOfNodes[interfaceBoundaryNodeList[iNode]] =
            keyOfInterfaceNodesInBlocks[blockIndex][localToGlobalInterfaceNodesMap[iNode]];
    }
}

void OversetConfig::UpdateInterfaceNodeCharacter()
{
    int blockIndex = grid->GetIBlock();
    int **keyOfInterfaceNodesInBlocks = GetOversetConfigFactoryMulti()->GetKeyOfInterfaceNodesInBlocks();

    int numberOfInterfaceBoundaryNodes = interfaceBoundaryNodeList.size();
    for (int iNode = 0; iNode < numberOfInterfaceBoundaryNodes; ++ iNode)
    {
        int &valueGlobal = keyOfInterfaceNodesInBlocks[blockIndex][localToGlobalInterfaceNodesMap[iNode]];
        valueGlobal = MAX(valueGlobal, keyActiveOfNodes[interfaceBoundaryNodeList[iNode]]);
    }
}

void OversetConfig::ShowBoundaryCells()
{
    string fileName = AddSymbolToFileName("./results/boundaryCells.dat", zoneIndex);

    fstream file;

    PHSPACE::OpenFile(file, fileName, ios_base::out);

    int numberOfNodes = grid->GetNTotalNode();
    int numberOfFaces = grid->GetNTotalFace();
    int numberOfCells = grid->GetNTotalCell();

    RDouble *x = grid->GetX();
    RDouble *y = grid->GetY();
    RDouble *z = grid->GetZ();

    int *faceNodeNumberContainer = grid->GetNodeNumberOfEachFace();
    int *faceNodeIndexContainer = grid->GetFace2Node();
    int *leftCellIndexContainer = grid->GetLeftCellOfFace();
    int *rightCellIndexContainer = grid->GetRightCellOfFace();

    int numberOfNodesInZone = 0;
    int numberOfFacesInZone = 0;
    int numberOfCellsInZone = 0;

    vector< int > keyNodes;
    vector< int > keyFaces;
    vector< int > nodeList;
    vector< int > cellList;

    keyNodes.resize(numberOfNodes);
    keyFaces.resize(numberOfFaces);
    nodeList.resize(numberOfNodes);
    cellList.resize(numberOfCells);

    SetField(keyNodes, 0, numberOfNodes);
    SetField(keyFaces, 0, numberOfFaces);
    SetField(nodeList, -1, numberOfNodes);
    SetField(cellList, 0, numberOfCells);

    int **cellNodeIndexContainer = grid->GetCell2NodeArray();
    int *cellNodeNumberContainer = grid->GetNodeNumberOfEachCell();

    int **cellFaceIndexContainer = grid->GetCell2Face();

    int *cellFaceNumberContainer = grid->GetFaceNumberOfEachCell();

    for (int iCell = 0; iCell < numberOfCells; ++iCell)
    {
        bool innerCellOrNot = false;
        for (int iNode = 0; iNode < cellNodeNumberContainer[iCell]; ++iNode)
        {
            int nodeIndex = cellNodeIndexContainer[iCell][iNode];
            if (keyInBlockNodes[nodeIndex])
            {
                innerCellOrNot = true;
                continue;
            }
        }
        if (innerCellOrNot)
        {
            numberOfCellsInZone++;
            cellList[iCell] = numberOfCellsInZone;

            for (int iNodes = 0; iNodes < cellNodeNumberContainer[iCell]; ++iNodes)
            {
                keyNodes[cellNodeIndexContainer[iCell][iNodes]] = 1;
            }

            for (int iFaces = 0; iFaces < cellFaceNumberContainer[iCell]; ++iFaces)
            {
                keyFaces[cellFaceIndexContainer[iCell][iFaces]] = 1;
            }
        }
    }

    for (int iNode = 0; iNode < numberOfNodes; ++iNode)
    {
        if (keyNodes[iNode] == 1)
        {
            numberOfNodesInZone++;
            nodeList[iNode] = numberOfNodesInZone;
        }
    }

    for (int iFace = 0; iFace < numberOfFaces; ++iFace)
    {
        if (keyFaces[iFace] == 1)
        {
            numberOfFacesInZone++;
        }
    }

    if (numberOfCellsInZone == 0)
    {
        return;
    }

    file << "title     = \"overset grid \"\n";
    file << "variables = \"x\" \"y\" \"z\"\n";

    file << "ZONE T=\"Boundary Cells\" \n";

    if (PHSPACE::GetDim() == THREE_D)
    {
        file << "ZoneType = FEPolyhedron\n";
    }
    else
    {
        file << "ZoneType = FEPolygon\n";
    }

    int totalNumFaceNodes = 0;

    for (int iFace = 0; iFace < numberOfFaces; ++iFace)
    {
        if (keyFaces[iFace] != 1) continue;
        totalNumFaceNodes += faceNodeNumberContainer[iFace];
    }

    file << "Nodes    = " << numberOfNodesInZone << "\n";
    file << "Faces    = " << numberOfFacesInZone << "\n";
    file << "Elements = " << numberOfCellsInZone << "\n";
    file << "TotalNumFaceNodes = " << totalNumFaceNodes << "\n";
    file << "NumConnectedBoundaryFaces = 0\n";
    file << "TotalNumBoundaryConnections = 0\n";

    int numberOfWordsInEachLine = 5;
    int iPrint = 0;

    //nodes coordinate x
    for (int iNode = 0; iNode < numberOfNodes; ++iNode)
    {
        if (keyNodes[iNode] != 1) continue;
        file << x[iNode] << " ";
        iPrint++;
        if (iPrint % numberOfWordsInEachLine == 0) file << "\n";
    }
    if (iPrint % numberOfWordsInEachLine != 0) file << "\n";

    //nodes coordinate y
    iPrint = 0;
    for (int iNode = 0; iNode < numberOfNodes; ++iNode)
    {
        if (keyNodes[iNode] != 1) continue;
        file << y[iNode] << " ";
        iPrint++;
        if (iPrint % numberOfWordsInEachLine == 0) file << "\n";
    }
    if (iPrint % numberOfWordsInEachLine != 0) file << "\n";

    //nodes coordinate z
    iPrint = 0;
    for (int iNode = 0; iNode < numberOfNodes; ++iNode)
    {
        if (keyNodes[iNode] != 1) continue;
        file << z[iNode] << " ";
        iPrint++;
        if (iPrint % numberOfWordsInEachLine == 0) file << "\n";
    }
    if (iPrint % numberOfWordsInEachLine != 0) file << "\n";

    //face nodes number
    iPrint = 0;
    if (PHSPACE::GetDim() == THREE_D)
    {
        for (int iFace = 0; iFace < numberOfFaces; ++iFace)
        {
            if (keyFaces[iFace] != 1) continue;
            file << faceNodeNumberContainer[iFace] << " ";
            iPrint++;
            if (iPrint % numberOfWordsInEachLine == 0) file << "\n";
        }
        if (iPrint % numberOfWordsInEachLine != 0) file << "\n";
    }

    //face nodes list
    int count = 0;
    iPrint = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++iFace)
    {
        int faceNodeNumber = faceNodeNumberContainer[iFace];

        if (keyFaces[iFace] == 0)
        {
            count += faceNodeNumber;
            continue;
        }

        for (int iNode = 0; iNode < faceNodeNumberContainer[iFace]; ++iNode)
        {
            int nodeIndex = nodeList[faceNodeIndexContainer[count++]];
            file << nodeIndex << " ";
            iPrint++;
            if (iPrint % numberOfWordsInEachLine == 0) file << "\n";
        }
    }
    if (iPrint % numberOfWordsInEachLine != 0) file << "\n";

    //left cell
    iPrint = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++iFace)
    {
        if (keyFaces[iFace] == 0) continue;

        int cellIndex = cellList[leftCellIndexContainer[iFace]];
        file << cellIndex << " ";
        iPrint++;
        if (iPrint % numberOfWordsInEachLine == 0) file << "\n";
    }
    if (iPrint % numberOfWordsInEachLine != 0) file << "\n";

    //right cell
    iPrint = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++iFace)
    {
        int cellIndex;
        if (keyFaces[iFace] == 0) continue;
        int rightCell = rightCellIndexContainer[iFace];

        if (rightCell >= numberOfCells || rightCell < 0)
        {
            cellIndex = 0;
        }
        else
        {
            cellIndex = cellList[rightCellIndexContainer[iFace]];
        }

        file << cellIndex << " ";
        iPrint++;
        if (iPrint % numberOfWordsInEachLine == 0) file << "\n";
    }
    if (iPrint % numberOfWordsInEachLine != 0) file << "\n";

    PHSPACE::CloseFile(file);
}

void OversetConfig::OutPutOversetVisualization(const string &name)
{
    ostringstream oss;
    oss << "OversetOutput/" << name << "/" << grid->GetZoneID() << ".dat";

    string fileName = oss.str();

    fstream file;

    OpenFile(file, fileName, ios_base::out | ios_base::trunc);

    OutPutOversetVisualizationHeader(file);
    OutPutOversetVisualizationCell(file);
    OutPutOversetVisualizationPoint(file);
    OutPutOversetVisualizationTail(file);

    CloseFile(file);
}

void OversetConfig::OutPutOversetVisualizationHeader(fstream &file)
{
    int numberOfNodes = grid->GetNTotalNode();
    int numberOfFaces = grid->GetNTotalFace();
    int numberOfCells = grid->GetNTotalCell();

    RDouble *nodeCoordinateX = grid->GetX();
    RDouble *nodeCoordinateY = grid->GetY();
    RDouble *nodeCoordinateZ = grid->GetZ();

    int *faceNodeNumber = grid->GetNodeNumberOfEachFace();

    file << "Title = \"OverSetVisualization\"" << endl << endl;

    file << "Variables = \"X\" \"Y\" \"Z\"" << endl;
    file << "\"keyActiveOfCells\" \"keyBufferOfCells\" \"keySearchOfCells\" \"interpolateCellList\"" << endl;
    file << "\"keyActiveOfNodes\" \"keyBufferOfNodes\" \"keySearchOfNodes\" \"keyInBlockNodes\" " << endl;
    file << "\"keyBasicActiveNodes\"  \"nodesMinDistanceInOtherBlock\" " << endl;
    file << endl;

    file << "Zone T    = \"Zone_" << grid->GetZoneID() << "\"" << endl << endl;
    file << "VARLOCATION = ( [1-3,8-13]=NODAL, [4-7]=CELLCENTERED )" << endl << endl;
    if (THREE_D == geometricDimension)
    {
        file << "ZoneType = FEPolyhedron" << endl << endl;
    }
    else
    {
        file << "ZoneType = FEPolygon" << endl << endl;
    }

    int totalNumberOfFaceNodes = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++iFace)
    {
        totalNumberOfFaceNodes += faceNodeNumber[iFace];
    }
    file << "Nodes                       = " << numberOfNodes << endl;
    file << "Faces                       = " << numberOfFaces << endl;
    file << "Elements                    = " << numberOfCells << endl;
    file << "TotalNumFaceNodes           = " << totalNumberOfFaceNodes << endl;
    file << "NumConnectedBoundaryFaces   = 0" << endl;
    file << "TotalNumBoundaryConnections = 0" << endl;
    file << endl;

    OutPutValueName(file, "nodeCoordinateX");
    OutPutPointerValue(file, numberOfNodes, nodeCoordinateX);

    OutPutValueName(file, "nodeCoordinateY");
    OutPutPointerValue(file, numberOfNodes, nodeCoordinateY);

    OutPutValueName(file, "nodeCoordinateZ");
    OutPutPointerValue(file, numberOfNodes, nodeCoordinateZ);
}

void OversetConfig::OutPutOversetVisualizationCell(fstream &file)
{
    int numberOfCells = grid->GetNTotalCell();

    OutPutValueName(file, "keyActiveOfCells");
    OutPutPointerValue(file, numberOfCells, GetKeyActiveOfCells());

    OutPutValueName(file, "keyBufferOfCells");
    OutPutPointerValue(file, numberOfCells, keyBufferOfCells);

    OutPutValueName(file, "keySearchOfCells");
    OutPutPointerValue(file, numberOfCells, keySearchOfCells);

    OutPutValueName(file, "interpolateCellList");
    OutPutVectorValue(file, numberOfCells, interpolateCellList);
}

void OversetConfig::OutPutOversetVisualizationPoint(fstream &file)
{
    int numberOfNodes = grid->GetNTotalNode();

    OutPutValueName(file, "keyActiveOfNodes");
    OutPutPointerValue(file, numberOfNodes, keyActiveOfNodes);

    OutPutValueName(file, "keyBufferOfNodes");
    OutPutPointerValue(file, numberOfNodes, keyBufferOfNodes);

    OutPutValueName(file, "keySearchOfNodes");
    OutPutPointerValue(file, numberOfNodes, keySearchOfNodes);

    OutPutValueName(file, "keyInBlockNodes");
    OutPutPointerValue(file, numberOfNodes, keyInBlockNodes);

    OutPutValueName(file, "keyBasicActiveNodes");
    OutPutPointerValue(file, numberOfNodes, keyBasicActiveNodes);

    OutPutValueName(file, "nodesMinDistanceInOtherBlock");
    OutPutPointerValue(file, numberOfNodes, nodesMinDistanceInOtherBlock);

}

void OversetConfig::OutPutOversetVisualizationTail(fstream &file)
{
    int    numberOfFaces = grid->GetNTotalFace();
    int    numberOfCells = grid->GetNTotalCell();

    int *faceNodeNumber = grid->GetNodeNumberOfEachFace();
    int **faceNodeIndex = grid->GetFace2NodeArray();
    int *leftCellIndex = grid->GetLeftCellOfFace();
    int *rightCellIndex = grid->GetRightCellOfFace();

    if (THREE_D == geometricDimension)
    {
        OutPutValueName(file, "faceNodeNumber");
        OutPutPointerValue(file, numberOfFaces, faceNodeNumber);
    }

    int numberOfWordsInEachLine = 5;
    int iPrint = 0;

    OutPutValueName(file, "faceNodeIndex");
    iPrint = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++iFace)
    {
        for (int iNode = 0; iNode < faceNodeNumber[iFace]; ++iNode)
        {
            int nodeIndex = faceNodeIndex[iFace][iNode] + 1;
            file << nodeIndex << " ";
            ++iPrint;
            if (iPrint % numberOfWordsInEachLine == 0) file << endl;
        }
    }
    if (iPrint % numberOfWordsInEachLine != 0) file << endl;
    file << endl;

    OutPutValueName(file, "leftCellIndex");
    iPrint = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++iFace)
    {
        int cellIndex = leftCellIndex[iFace] + 1;
        if (cellIndex > numberOfCells || cellIndex < 0) cellIndex = 0;
        file << cellIndex << " ";
        ++iPrint;
        if (iPrint % numberOfWordsInEachLine == 0) file << endl;
    }
    if (iPrint % numberOfWordsInEachLine != 0) file << endl;
    file << endl;

    OutPutValueName(file, "rightCellIndex");
    iPrint = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++iFace)
    {
        int cellIndex = rightCellIndex[iFace] + 1;
        if (cellIndex > numberOfCells || cellIndex < 0) cellIndex = 0;
        file << cellIndex << " ";
        ++iPrint;
        if (iPrint % numberOfWordsInEachLine == 0) file << endl;
    }
    if (iPrint % numberOfWordsInEachLine != 0) file << endl;
    file << endl;
}

void OversetConfig::ShowGrid()
{
    string fileName = AddSymbolToFileName("./results/oversetGrid.dat", zoneIndex);

    fstream file;

    PHSPACE::OpenFile(file, fileName, ios_base::out);

    PrintGridByOverSetConditionCell(file, -1);
    PrintGridByOverSetConditionCell(file, 1);

    PHSPACE::CloseFile(file);
}

void OversetConfig::PrintGridByOverSetConditionCell(fstream &file, int key)
{
    int numberOfNodes = grid->GetNTotalNode();
    int numberOfFaces = grid->GetNTotalFace();
    int numberOfCells = grid->GetNTotalCell();

    RDouble *x = grid->GetX();
    RDouble *y = grid->GetY();
    RDouble *z = grid->GetZ();

    int *faceNodeNumberContainer = grid->GetNodeNumberOfEachFace();
    int *faceNodeIndexContainer = grid->GetFace2Node();
    int *leftCellIndexContainer = grid->GetLeftCellOfFace();
    int *rightCellIndexContainer = grid->GetRightCellOfFace();

    int numberOfNodesInZone = 0;
    int numberOfFacesInZone = 0;
    int numberOfCellsInZone = 0;

    vector< int > keyNodes;
    vector< int > keyFaces;
    vector< int > nodeList;
    vector< int > cellList;

    keyNodes.resize(numberOfNodes);
    keyFaces.resize(numberOfFaces);
    nodeList.resize(numberOfNodes);
    cellList.resize(numberOfCells);

    SetField(keyNodes, 0, numberOfNodes);
    SetField(keyFaces, 0, numberOfFaces);
    SetField(nodeList, -1, numberOfNodes);
    SetField(cellList, 0, numberOfCells);

    int **cellNodeIndexContainer = grid->GetCell2NodeArray();
    int *cellNodeNumberContainer = grid->GetNodeNumberOfEachCell();

    int **cellFaceIndexContainer = grid->GetCell2Face();
    int *cellFaceNumberContainer = grid->GetFaceNumberOfEachCell();

    int *keyActiveOfCells = this->GetKeyActiveOfCells();

    for (int iCell = 0; iCell < numberOfCells; ++iCell)
    {
        if (keyActiveOfCells[iCell] == key)
        {
            numberOfCellsInZone++;
            cellList[iCell] = numberOfCellsInZone;

            for (int iNodes = 0; iNodes < cellNodeNumberContainer[iCell]; ++iNodes)
            {
                keyNodes[cellNodeIndexContainer[iCell][iNodes]] = 1;
            }

            for (int iFaces = 0; iFaces < cellFaceNumberContainer[iCell]; ++iFaces)
            {
                keyFaces[cellFaceIndexContainer[iCell][iFaces]] = 1;
            }
        }
    }

    for (int iNode = 0; iNode < numberOfNodes; ++iNode)
    {
        if (keyNodes[iNode] == 1)
        {
            numberOfNodesInZone++;
            nodeList[iNode] = numberOfNodesInZone;
        }
    }

    for (int iFace = 0; iFace < numberOfFaces; ++iFace)
    {
        if (keyFaces[iFace] == 1)
        {
            numberOfFacesInZone++;
        }
    }

    if (numberOfCellsInZone == 0)
    {
        return;
    }

    file << "title     = \"overset grid \"\n";
    file << "variables = \"x\" \"y\" \"z\"\n";

    if (key == 0)
    {
        file << "ZONE T=\"Negative Cells\" \n";
    }
    else if (key == 1)
    {
        file << "ZONE T=\"Active Cells\" \n";
    }
    else if (key == -1)
    {
        file << "ZONE T=\"Interpolation Cells\" \n";
    }

    if (PHSPACE::GetDim() == THREE_D)
    {
        file << "ZoneType = FEPolyhedron\n";
    }
    else
    {
        file << "ZoneType = FEPolygon\n";
    }

    int totalNumFaceNodes = 0;

    for (int iFace = 0; iFace < numberOfFaces; ++iFace)
    {
        if (keyFaces[iFace] != 1) continue;
        totalNumFaceNodes += faceNodeNumberContainer[iFace];
    }

    file << "Nodes    = " << numberOfNodesInZone << "\n";
    file << "Faces    = " << numberOfFacesInZone << "\n";
    file << "Elements = " << numberOfCellsInZone << "\n";
    file << "TotalNumFaceNodes = " << totalNumFaceNodes << "\n";
    file << "NumConnectedBoundaryFaces = 0\n";
    file << "TotalNumBoundaryConnections = 0\n";

    int numberOfWordsInEachLine = 5;
    int iPrint = 0;

    //nodes coordinate x
    for (int iNode = 0; iNode < numberOfNodes; ++iNode)
    {
        if (keyNodes[iNode] != 1) continue;
        file << x[iNode] << " ";
        iPrint++;
        if (iPrint % numberOfWordsInEachLine == 0) file << "\n";
    }
    if (iPrint % numberOfWordsInEachLine != 0) file << "\n";

    //nodes coordinate y
    iPrint = 0;
    for (int iNode = 0; iNode < numberOfNodes; ++iNode)
    {
        if (keyNodes[iNode] != 1) continue;
        file << y[iNode] << " ";
        iPrint++;
        if (iPrint % numberOfWordsInEachLine == 0) file << "\n";
    }
    if (iPrint % numberOfWordsInEachLine != 0) file << "\n";

    //nodes coordinate z
    iPrint = 0;
    for (int iNode = 0; iNode < numberOfNodes; ++iNode)
    {
        if (keyNodes[iNode] != 1) continue;
        file << z[iNode] << " ";
        iPrint++;
        if (iPrint % numberOfWordsInEachLine == 0) file << "\n";
    }
    if (iPrint % numberOfWordsInEachLine != 0) file << "\n";

    //face nodes number
    iPrint = 0;
    if (PHSPACE::GetDim() == THREE_D)
    {
        for (int iFace = 0; iFace < numberOfFaces; ++iFace)
        {
            if (keyFaces[iFace] != 1) continue;
            file << faceNodeNumberContainer[iFace] << " ";
            iPrint++;
            if (iPrint % numberOfWordsInEachLine == 0) file << "\n";
        }
        if (iPrint % numberOfWordsInEachLine != 0) file << "\n";
    }

    //face nodes list
    int count = 0;
    iPrint = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++iFace)
    {
        int faceNodeNumber = faceNodeNumberContainer[iFace];

        if (keyFaces[iFace] == 0)
        {
            count += faceNodeNumber;
            continue;
        }

        for (int iNode = 0; iNode < faceNodeNumberContainer[iFace]; ++iNode)
        {
            int nodeIndex = nodeList[faceNodeIndexContainer[count++]];
            file << nodeIndex << " ";
            iPrint++;
            if (iPrint % numberOfWordsInEachLine == 0) file << "\n";
        }
    }
    if (iPrint % numberOfWordsInEachLine != 0) file << "\n";

    //left cell
    iPrint = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++iFace)
    {
        if (keyFaces[iFace] == 0) continue;

        int cellIndex = cellList[leftCellIndexContainer[iFace]];
        file << cellIndex << " ";
        iPrint++;
        if (iPrint % numberOfWordsInEachLine == 0) file << "\n";
    }
    if (iPrint % numberOfWordsInEachLine != 0) file << "\n";

    //right cell
    iPrint = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++iFace)
    {
        int cellIndex;
        if (keyFaces[iFace] == 0) continue;
        int rightCell = rightCellIndexContainer[iFace];

        if (rightCell >= numberOfCells || rightCell < 0)
        {
            cellIndex = 0;
        }
        else
        {
            cellIndex = cellList[rightCellIndexContainer[iFace]];
        }

        file << cellIndex << " ";
        iPrint++;
        if (iPrint % numberOfWordsInEachLine == 0) file << "\n";
    }
    if (iPrint % numberOfWordsInEachLine != 0) file << "\n";
}

void OversetConfig::PrintGridByOverSetConditionPoint(fstream &file, int key)
{
    int numberOfNodes = grid->GetNTotalNode();

    RDouble *x = grid->GetX();
    RDouble *y = grid->GetY();
    RDouble *z = grid->GetZ();

    int numberOfSpecifyNodes = 0;

    vector< int > nodeList;

    for (int iNode = 0; iNode < numberOfNodes; ++iNode)
    {
        if (keyActiveOfNodes[iNode] == key)
        {
            numberOfSpecifyNodes++;
            nodeList.push_back(iNode);
        }
    }

    if (numberOfSpecifyNodes != 0)
    {
        file << "Title=\"THE SPECIFY NODES\" \n";
        file << "Variables=\"x\" \"y\" \"z\" \n";

        if (key == 1)
        {
            file << "ZONE T=\"Active Nodes\" \n";
        }
        else if (key == -1)
        {
            file << "ZONE T=\"Negative Nodes\" \n";
        }

        file << "N=" << numberOfSpecifyNodes << " E=" << numberOfSpecifyNodes << " F=FEPOINT, ET=quadrilateral\n";

        for (int iSpecifyNodes = 0; iSpecifyNodes < numberOfSpecifyNodes; ++iSpecifyNodes)
        {
            int activeNodeIndex = nodeList[iSpecifyNodes];

            file << x[activeNodeIndex] << " ";
            file << y[activeNodeIndex] << " ";
            file << z[activeNodeIndex] << " ";
            file << "\n";
        }

        for (int iSpecifyNodes = 0; iSpecifyNodes < numberOfSpecifyNodes; ++iSpecifyNodes)
        {
            file << iSpecifyNodes + 1 << " ";
            file << iSpecifyNodes + 1 << " ";
            file << iSpecifyNodes + 1 << " ";
            file << iSpecifyNodes + 1 << " ";
            file << "\n";
        }
    }
}

void OversetConfig::OutKeyBuffer()
{
    int numberOfNodes = grid->GetNTotalNode();
    int numberOfCells = grid->GetNTotalCell();

    string fileName = AddSymbolToFileName("keyBuffer.dat", zoneIndex);
    fstream file;

    PHSPACE::OpenFile(file, fileName, ios_base::out);

    for (int iCell = 0; iCell < numberOfCells; ++iCell)
    {
        file << keyBufferOfCells[iCell] << endl;
    }

    for (int iNode = 0; iNode < numberOfNodes; ++iNode)
    {
        file << keyBufferOfNodes[iNode] << endl;
    }

    CloseFile(file);
}

void OversetConfig::OutKeyActive()
{
    int numberOfNodes = grid->GetNTotalNode();
    string fileName = AddSymbolToFileName("keyActive.dat", zoneIndex);
    fstream file;

    PHSPACE::OpenFile(file, fileName, ios_base::out);

    file << "keyActiveOfNodes" << endl;
    for (int iNode = 0; iNode < numberOfNodes; ++iNode)
    {
        file << keyActiveOfNodes[iNode] << endl;
    }

    CloseFile(file);
}

void OversetConfig::OutputLocalNodeWalldistance()
{
    int numberOfNodes = grid->GetNTotalNode();

    string fileName = AddSymbolToFileName("localNodeWalldistance.dat", zoneIndex);
    fstream file;

    PHSPACE::OpenFile(file, fileName, ios_base::out);

    for (int iNode = 0; iNode < numberOfNodes; ++iNode)
    {
        file << nodeWalldistDistance[iNode] << endl;
    }

    PHSPACE::CloseFile(file);
}

int *OversetConfig::GetKeyActiveOfCells()
{
    int *keyActiveOfCells = this->grid->GetBlankIndex();
    return keyActiveOfCells;
}

RDouble GetToleranceForOversetBox()
{
    return GlobalDataBase::GetDoubleParaFromDB("toleranceForOversetBox");
}

int GetTwoOrderInterpolationOrNot()
{
    return GlobalDataBase::GetIntParaFromDB("twoOrderInterpolationOrNot");
}

int GetKeyEnlargeOfNodes()
{
    return GlobalDataBase::GetIntParaFromDB("keyEnlargeOfActiveNodes");
}

int GetReadInAuxiliaryInnerGrid()
{
    return GlobalDataBase::GetIntParaFromDB("readInAuxiliaryInnerGrid");
}

int GetReadInAuxiliaryOuterGrid()
{
    return GlobalDataBase::GetIntParaFromDB("readInAuxiliaryOuterGrid");
}

int GetReadInSklFileOrNot()
{
    return GlobalDataBase::GetIntParaFromDB("readInSklFileOrNot");
}

int GetSymetryOrNot()
{
    return GlobalDataBase::GetIntParaFromDB("symetryOrNot");
}

int GetIfOutPutOversetVisualization()
{
    return GlobalDataBase::GetIntParaFromDB("outPutOversetVisualization");
}

template< typename T >
void PH_AllreduceSepModeForVector(vector< T > &sendingBuffer, vector< T > &receivingBuffer, int numberOfElements, MPI_Datatype mpiDataType, PH_Op op)
{
    T *send = &sendingBuffer[0];
    T *recv = &receivingBuffer[0];

    PH_AllreduceSepMode(send, recv, numberOfElements, mpiDataType, op);
}

}
