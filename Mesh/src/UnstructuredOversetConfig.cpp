#include "UnstructuredOversetConfig.h"
#include "Pre_WalldistCompute.h"
#include "UnstructuredGridQuickSearch.h"
#include "IO_FileName.h"
#include "Glb_Dimension.h"
#include "PHIO.h"
#include "IO_HDF5File.h"
#include "Region.h"
#include "Task.h"
#include "Task_ServerUpdateInterface.h"
using namespace std;

namespace PHSPACE
{

OversetConfigFactoryMulti *oversetConfigFactory = NULL;

Region *region;

void RunOversetGridConfig()
{
    region = new Region;

    region->InitGridForOversetConfig();

    RunUnstructuredOversetGridConfig();

    //delete region;
}

void RunUnstructuredOversetGridConfig()
{
    ReConfigUnstructOversetGrid();

    OutUnstructOversetGridConfig();

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

OversetConfigFactoryMulti::OversetConfigFactoryMulti()
{
    numberOfBlocks = 0;
    numberOfZones = 0;

    oversetConfigMultiZone = 0;
    globalZoneNodes = 0;
    globalInterCells = 0;

    wallDistCompute = 0;
    testTime = 0;
    AdtMethodOrNod = 0;
}

OversetConfigFactoryMulti::~OversetConfigFactoryMulti()
{
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        if (oversetConfigMultiZone[iZone]) delete oversetConfigMultiZone[iZone];
        if (globalZoneNodes[iZone]) delete globalZoneNodes[iZone];
        if (globalInterCells[iZone]) delete globalInterCells[iZone];
    }

    delete [] oversetConfigMultiZone;
    delete [] globalZoneNodes;
    delete [] globalInterCells;

    for (unsigned iBlock = 0; iBlock < innerGrid.size(); ++ iBlock)
    {
        if (innerGrid[iBlock]) delete innerGrid[iBlock];
    }
    innerGrid.resize(0);

    for (unsigned iBlock = 0; iBlock < outerAuxilaryGrid.size(); ++ iBlock)
    {
        if (outerAuxilaryGrid[iBlock]) delete outerAuxilaryGrid[iBlock];
    }
    outerAuxilaryGrid.resize(0);

    delete wallDistCompute;
}

void OversetConfigFactoryMulti::Initialize()
{
    InitializeGlobalValue();    //! Initialize the necessary global variables.

    SearchPhysicalAndInterfaceBoundary();    //! Find and labal the points on physical boundary, inner boundary and interface boundary.

    if (GetReadInSklFileOrNot() == 1)
    {
        ReadInSklFiles();             //! Read in partition information file (used for communication node information).
        GetGlobalInterfaceNodes();    //! Interface boundary points on unpartitioned grids.
    }

    ComputeNodeWalldist();                 //! Minimum wall distance of nodes.
    
    ComputeNodeAndCellTopology();          //! Calculate the topological relationship of point to point and volume to point.
    
    ReadInInnerGrid();                     //! Read in the internal grid of the object to accurately estimate whether it is in the object.
    
    ReadInOuterGrid();                     //! Read in the external grid of the object for quick query.
    
    SetBufferRegionOnOversetBoundary();    //! Set buffer on the boundary of overset.
}

void OversetConfigFactoryMulti::InitializeGlobalValue()
{
    int currentProcessor = PHMPI::GetCurrentProcessorID();

    numberOfBlocks = PHSPACE::GlobalDataBase::GetIntParaFromDB("numberOfGridGroups");

    numberOfZones = PHMPI::GetNumberofGlobalZones();

    oversetConfigMultiZone = new OversetConfig * [numberOfZones];
    globalZoneNodes = new ZoneNodes * [numberOfZones];
    globalInterCells = new InterCells * [numberOfZones];

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorIDSepMode(iZone);

        ZoneNodes *zoneNodes = new ZoneNodes;
        InterCells *interCells = new InterCells;

        globalZoneNodes[iZone] = zoneNodes;
        globalInterCells[iZone] = interCells;

        if (currentProcessor == processorIndex)
        {
            OversetConfig *oversetGrid = new OversetConfig;
            Grid *grid = GetGrid(iZone, 0);
            UnstructGrid *gridUnstr = PHSPACE::UnstructGridCast(grid);
            int iBlock = grid->GetIBlock();

            oversetGrid->Initialize(this, gridUnstr, numberOfZones, iZone, iBlock, globalZoneNodes, globalInterCells);
            oversetConfigMultiZone[iZone] = oversetGrid;

            gridUnstr->SetOversetConfig(oversetGrid);
        }
        else
        {
            oversetConfigMultiZone[iZone] = 0;
        }
    }
}

void OversetConfigFactoryMulti::SearchPhysicalAndInterfaceBoundary()
{
    //! Check physical boundary points and overset boundary points.
    int currentProcessor = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorIDSepMode(iZone);
        if (currentProcessor == processorIndex)
        {
            oversetConfigMultiZone[iZone]->SearchPhysicalAndInterfaceBoundary();
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
    file.read((char *) &nZone, sizeInt);
    CloseFile(file);

    if (nZone <= static_cast <int> (1e8))
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
    delete [] iDimensions;

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
    int *nodePosi = new int [nTotalFace + 1];
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
    delete [] nodePosi;

    int nBoundFace = 0;
    ReadData(grpData, &nBoundFace, "nBoundFace");
    gridIn->SetNBoundFace(nBoundFace);

    int *BcType = ReadIntOneRow(grpData, "BcType");
    CompositeUnstructBC **bcr = new CompositeUnstructBC * [nBoundFace];
    for (int iFace = 0; iFace < nBoundFace; ++ iFace)
    {
        bcr[iFace] = new CompositeUnstructBC();
        bcr[iFace]->SetKey(BcType[iFace]);
    }
    gridIn->SetBCRecord(bcr);
    delete [] BcType;

    gridIn->SpecifyBC();

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

        for (int iface = 0; iface < nIFace; ++ iface)
        {
            interFace2ZoneID[iface] = interFace2ZoneIDTemp[iface];
            interFace2InterFaceID[iface] = interFace2InterFaceIDTemp[iface];
            interFace2BoundaryFace[iface] = interFace2BoundaryFaceTemp[iface];
        }

        delete [] interFace2ZoneIDTemp;
        delete [] interFace2InterFaceIDTemp;
        delete [] interFace2BoundaryFaceTemp;
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
    string *BCNameList = new string [nBFaces];
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

    delete [] BCNameList;
    delete [] BCName;
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

    for (int iBlock = 0; iBlock < numberOfBlocks; ++ iBlock)
    {
        ostringstream oss;

        oss << "auxiliaryInnerGrid" << iBlock;

        int fileID = PHMPI::GetFileIndexofCurrentProcessor();
        string fileName = AddSymbolToFileName(GlobalDataBase::GetStrParaFromDB(oss.str()), "_", fileID);

        if (!FileExist(fileName))
        {
            innerGrid[iBlock] = 0;
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
            continue;
        }

        fstream file;

        DataContainer *dataContainer = new DataContainer;

        if (currentProcessor == serverProcessor)
        {
            cout << "Read in auxiliary inner grid" << fileName << endl;

            OpenFile(file, fileName, ios_base::in | ios_base::binary);

            int nZones;

            vector < int > zoneProcessorIndexContainer;
            vector < int > zoneIndexContainer;
            vector < int > zoneTypeContainer;

            PHRead(file, nZones);

            zoneProcessorIndexContainer.resize(nZones);
            zoneIndexContainer.resize(nZones);
            zoneTypeContainer.resize(nZones);

            PHRead(file, zoneProcessorIndexContainer, nZones);
            PHRead(file, zoneIndexContainer, nZones);
            PHRead(file, zoneTypeContainer, nZones);

            dataContainer->readfile(file);

            CloseFile(file);
        }

        int tag = 0;

        PHMPI::PH_BcastByServerSepMode(dataContainer, tag);

        innerGrid[iBlock]->Decode(dataContainer, 0);

        delete dataContainer;
    }
}

void OversetConfigFactoryMulti::ReadInOuterGrid()
{
    if (GetReadInAuxiliaryOuterGrid() != 1)
    {
        return;
    }

    int currentProcessor = PHMPI::GetCurrentProcessorID();
    int serverProcessor = PHMPI::GetServerSepMode();

    outerAuxilaryGrid.resize(numberOfBlocks);

    for (int iBlock = 0; iBlock < numberOfBlocks; ++ iBlock)
    {
        ostringstream oss;

        oss << "auxiliaryOuterGrid" << iBlock;

        string fileName = GlobalDataBase::GetStrParaFromDB(oss.str());

        if (!FileExist(fileName))
        {
            outerAuxilaryGrid[iBlock] = 0;
            continue;
        }

        fstream file;

        DataContainer *dataContainer = new DataContainer;

        if (currentProcessor == serverProcessor)
        {
            cout << "Read in auxiliary outer grid " << fileName << endl;

            OpenFile(file, fileName, ios_base::in | ios_base::binary);

            int nZones;

            vector < int > zoneProcessorIndexContainer;
            vector < int > zoneIndexContainer;
            vector < int > zoneTypeContainer;

            PHRead(file, nZones);

            zoneProcessorIndexContainer.resize(nZones);
            zoneIndexContainer.resize(nZones);
            zoneTypeContainer.resize(nZones);

            PHRead(file, zoneProcessorIndexContainer, nZones);
            PHRead(file, zoneIndexContainer, nZones);
            PHRead(file, zoneTypeContainer, nZones);

            dataContainer->readfile(file);

            CloseFile(file);
        }

        int tag = 0;

        PHMPI::PH_BcastByServerSepMode(dataContainer, tag);

        UnstructGrid *grid = new UnstructGrid();

        grid->InitGrid(new GridID(0, 0), 0, GetDim(), UNSTRUCTGRID);

        grid->Decode(dataContainer, 0);

        grid->SetIBlock(iBlock);

        outerAuxilaryGrid[iBlock] = new OversetAuxilaryGrid;

        outerAuxilaryGrid[iBlock]->Initialize(grid, iBlock, this->GetWallDistCompute());

        delete dataContainer;
    }

    //! Calculate the minimum wall distance of external auxiliary grid nodes (Cell).
    for (int iBlock = 0; iBlock < numberOfBlocks; ++ iBlock)
    {
        PrintToWindow("\nCal Wall dist for outerAuxilaryGrid: ", iBlock, "\n");

        outerAuxilaryGrid[iBlock]->CalMetrics();
        outerAuxilaryGrid[iBlock]->CalCellWallDistance();
        outerAuxilaryGrid[iBlock]->CalGradientOfWallDistanceGGNode();
    }
}

void OversetConfigFactoryMulti::ComputeNodeWalldist()
{
    AdtMethodOrNod = true;

    int walldistMethod;

    GlobalDataBase::GetData("walldistMethod", &walldistMethod, PHINT, 1);

    int level = 0;

    wallDistCompute = new Pre_WalldistCompute(walldistMethod, level);

    wallDistCompute->ComputeNodeDistance();
}

void OversetConfigFactoryMulti::ComputeNodeAndCellTopology()
{
    int currentProcessor = PHMPI::GetCurrentProcessorID();

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorIDSepMode(iZone);
        if (currentProcessor == processorIndex)
        {
            oversetConfigMultiZone[iZone]->ComputeNodeAndCellTopology();
        }
    }
}

void OversetConfigFactoryMulti::ReadInSklFiles()
{
    if (numberOfZones == numberOfBlocks)
    {
        return;
    }

    int currentProcessor = PHMPI::GetCurrentProcessorID();
    int serverProcessor = PHMPI::GetServerSepMode();

    int numberOfGridGroups = PHSPACE::GlobalDataBase::GetIntParaFromDB("numberOfGridGroups");
    int zoneIndex = 0;

    numberOfGlobalNodes.resize(numberOfGridGroups);
    keyActiveOfGlobalInterNodes.resize(numberOfGridGroups);

    for (int iGridFile = 0; iGridFile < numberOfGridGroups; ++ iGridFile)
    {
        ostringstream oss;
        oss << "gridFileName" << iGridFile;
        string gridFileName = PHSPACE::GlobalDataBase::GetStrParaFromDB(oss.str());

        string mainName, extensionName;
        PHSPACE::GetFileNameExtension(gridFileName, mainName, extensionName, ".");

        ostringstream ossSkl;
        ossSkl << mainName << ".skl";
        string sklFileName = ossSkl.str();

        int nZones;
        int numberOfNodes;
        int numberOfFaces;
        int numberOfCells;

        fstream file;
        if (currentProcessor == serverProcessor)
        {
            cout << "sklFileName = " << sklFileName << endl;
            OpenFile(file, sklFileName, ios_base::in | ios_base::binary);

            PHSPACE::PHRead(file, nZones);
            PHSPACE::PHRead(file, numberOfNodes);
            PHSPACE::PHRead(file, numberOfFaces);
            PHSPACE::PHRead(file, numberOfCells);
        }

        PHMPI::PH_BcastSepMode(&nZones, sizeof(int), serverProcessor);
        PHMPI::PH_BcastSepMode(&numberOfNodes, sizeof(int), serverProcessor);

        numberOfGlobalNodes[iGridFile] = numberOfNodes;

        for (int iZone = 0; iZone < nZones; ++ iZone)
        {
            int numberOfLocalNodes;
            int numberOfLocalFaces;
            int numberOfLocalCells;

            int *listForRead = NULL;

            int receiveProcessor = PHMPI::GetZoneProcessorIDSepMode(zoneIndex);
            if (currentProcessor == serverProcessor)
            {
                PHSPACE::PHRead(file, numberOfLocalNodes);
                PHSPACE::PHRead(file, numberOfLocalFaces);
                PHSPACE::PHRead(file, numberOfLocalCells);

                listForRead = new int [numberOfLocalNodes];

                PHSPACE::PHRead(file, listForRead, numberOfLocalNodes);

                file.seekg((numberOfLocalFaces) * sizeof(int), ios::cur);
                file.seekg((numberOfLocalCells) * sizeof(int), ios::cur);
            }

            PHMPI::PH_Trade(&numberOfLocalNodes, 1, serverProcessor, receiveProcessor, iZone);

            if ((currentProcessor == receiveProcessor) && (receiveProcessor != serverProcessor))
            {
                listForRead = new int [numberOfLocalNodes];
            }

            PHMPI::PH_Trade(listForRead, numberOfLocalNodes, serverProcessor, receiveProcessor, zoneIndex);

            if ((currentProcessor == receiveProcessor))
            {
                vector< int > &localToGlobalNodesIndex = oversetConfigMultiZone[zoneIndex]->GetLocalToGlobalNodesIndex();
                localToGlobalNodesIndex.resize(numberOfLocalNodes);

                for (int iNode = 0; iNode < numberOfLocalNodes; ++ iNode)
                {
                    localToGlobalNodesIndex[iNode] = listForRead[iNode];
                }
            }

            delete [] listForRead;

            zoneIndex++;
        }

        if (currentProcessor == serverProcessor)
        {
            CloseFile(file);
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
    vector< MRADTSearch * > aDTSearchArray;

    aDTSearchArray.resize(numberOfBlocks);

    //! step1£ºInteral grid establishment ADT-Tree.
    for (int iBlock = 0; iBlock < numberOfBlocks; ++ iBlock)
    {
        if (innerGrid[iBlock] == 0)
        {
            aDTSearchArray[iBlock] = 0;
            continue;
        }

        //! Set CurrentActionKey to 0, avoid unnecessary error in ComputMetrices.
        PHSPACE::SetCurrentActionKey(0);

        innerGrid[iBlock]->ComputeMetrics();

        MRADTSearch *aDTSearch = new MRADTSearch(innerGrid[iBlock]);

        aDTSearch->BuildADTStructure();

        aDTSearchArray[iBlock] = aDTSearch;
    }

    //! step2£ºCheck the cell center in all innerGrid, sign it as cave cell if there is a host cell.
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorIDSepMode(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->SpecifyCellsInBody(aDTSearchArray);
        }
    }

    for (int iBlock = 0; iBlock < numberOfBlocks; ++ iBlock)
    {
        if (aDTSearchArray[iBlock] != 0)
        {
            delete aDTSearchArray[iBlock];
        }
    }

    //! Even if there's no communication here, may be some difference between partitioned and nonpartitioned interpolation boundaries, but it won't affect the results.
    //! The difference only lies in the local interpolation cell of the surface.
    //CommunicateNodeCharacterOnInterfaceBySKL(1);
}

void OversetConfigFactoryMulti::SpecifyNodesCharacterByDistanceAndInnerNodesCheck()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorIDSepMode(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->SpecifyNodesCharacterByDistanceAndInnerNodesCheck();
        }
    }
}

//! Expand a row of active points, points in the innerBody and buffer are set to inactive.
void OversetConfigFactoryMulti::InterFaceBoundaryOptimize()
{
    int keyEnlarge = GetKeyEnlargeOfNodes() + 1;

    for (int num = 0; num < keyEnlarge; ++ num)
    {
        CharacterNodeToCell();
        CharacterCellToNode();
        SetNegativeForInnerNodes();
    }

    //! Points in buffer are set to inactive.
    SetNegativeForBufferRegion();
}

void OversetConfigFactoryMulti::SetNegativeForInnerNodes()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorIDSepMode(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->SetNegativeForInnerNodes();
        }
    }
}

void OversetConfigFactoryMulti::SetNegativeForInnerAndInitialInterNodes()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorIDSepMode(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->SetNegativeForInnerAndInitialInterNodes();
        }
    }
}

void OversetConfigFactoryMulti::NodesOperator()
{
    if (GetReadInAuxiliaryOuterGrid() == 0)
    {
        CalMinDistByGlobalNodes();
    }
    else
    {
        CalMinDistByAuxiliaryOuterGrid();
    }

    SpecifyInitialInterNodesByDistance();

    if (GetReadInAuxiliaryInnerGrid() == 1)
    {
        SpecifyCellsInBodyByAuxiliaryInnerGrid();
    }
    else
    {
        SpecifySpecicalCellsInteractWithWallBoundary();
    }

    FillTheActiveRegion();          //! Physical boundary advancement.

    InterFaceBoundaryOptimize();    //! Hole boundary optimization (expand active areas).
}

void OversetConfigFactoryMulti::ReComputeWallDistance()
{
    int numberOfMovingBodies = GlobalDataBase::GetIntParaFromDB("numberOfMovingBodies");

    int key = 0;

    for (int iBody = 0; iBody < numberOfMovingBodies; ++ iBody)
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

    //! If the wall distance needs to be recalculated, the wall points are collected again, because the shape and displacement of the object have changed.
    delete wallDistCompute;

    ComputeNodeWalldist();
}

//! Collect globalNodes, check it in each area, and then communicate with MPI.
void OversetConfigFactoryMulti::CalMinDistByGlobalNodes()
{
    SetSearchRegionGlobal();

    BuildADTStructureAndMinBoxForSearch();

    CreatGlobalNodes();

    SpecifyTheDonorCellAndMinDistanceOfGlobalNodes();
}

//! Check in the auxiliary grid directly.
void OversetConfigFactoryMulti::CalMinDistByAuxiliaryOuterGrid()
{
    //SetSearchRegionGlobal();    //! Points in the non-overlapping area don't participate in any check or checked (Auxiliary grid, no needed for the moment).

    BuildADTStructureAndMinBoxForSearch();    //! Building a ADT tree, calculate the smallest box for each zone.

    //! step1. Calculating the wall distance of auxiliary grid points (This step is in initialization).
    for (int iBlock = 0; iBlock < numberOfBlocks; ++ iBlock)
    {
        //outerAuxilaryGrid[iBlock]->CalNodeWallDistance();
    }

    //! step2. Auxialiary grid to build ADT structure.
    for (int iBlock = 0; iBlock < numberOfBlocks; ++ iBlock)
    {
        outerAuxilaryGrid[iBlock]->CalMetrics();
        outerAuxilaryGrid[iBlock]->BuildADT();
    }

    //! step3. Calculate the minimum wall distance in the auxiliary grid.
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorIDSepMode(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->CalMinDistanceByAuxilaryGrid(outerAuxilaryGrid);
        }
    }
}

void OversetConfigFactoryMulti::CellsOperator()
{
    FindActiveCells();    //! Identify active cells.

    FindInterCells();     //! Identify interpolation type.

    CreatGlobalInterCellsServerCollectAndBcastNonBlocking();

    FindDonorZoneAndDonorCell();

    SelectDonorZoneAndDonorCell();

    CheckDonorZoneAndDonorCell();

    ResetOversetInformation();    //! Rebuild interpolation information.
}

void OversetConfigFactoryMulti::SetInitialValue()
{
    int currentProcessor = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorIDSepMode(iZone);
        if (currentProcessor == processorIndex)
        {
            oversetConfigMultiZone[iZone]->SetInitialValue();
        }
    }
}

void OversetConfigFactoryMulti::BuildADTStructureAndMinBoxForSearch()
{
    CalMinBoxForSearch();

    BuildADTStructure();
}

void OversetConfigFactoryMulti::CalMinBoxForSearch()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorIDSepMode(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->CalculateMinBox();
        }
    }
}

void OversetConfigFactoryMulti::BuildADTStructure()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorIDSepMode(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->BuildADTStructure();
        }
    }
}

void OversetConfigFactoryMulti::SetBufferRegionOnOversetBoundary()
{
    //! Query physical boundary points and overset boundary points.
    int currentProcessor = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorIDSepMode(iZone);
        if (currentProcessor == processorIndex)
        {
            oversetConfigMultiZone[iZone]->SetInitialValueForBufferRigion();
        }
    }

    CharacterNodeToCell();    //! A row of points and cells are left in the buffer.

    //! Increase the interpolation area to improve the reconstruction accuracy of interpolation cells, two row of points and cells will be left in the buffer (There's no necessary actually).
    int twoOrderInterpolationOrNot = GetTwoOrderInterpolationOrNot();
    if (twoOrderInterpolationOrNot == 1)
    {
        CharacterCellToNode();
        CharacterNodeToCell();
    }

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorIDSepMode(iZone);
        if (currentProcessor == processorIndex)
        {
            oversetConfigMultiZone[iZone]->SetKeyBufferByKeyActive();
        }
    }
}

//! Identify searchRegion, starting from search point and expand twice.
//! Identify search point by judging whether it's in other blocks.
void OversetConfigFactoryMulti::SetSearchRegionGlobal()
{
    CalMinMaxBoxOfBlocks();

    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorIDSepMode(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->SetSearchRegionGlobal(pMinOfBlocks, pMaxOfBlocks);
        }
    }

    //! Why the value is 2 (Ensure that searchRegion covers all blocks except this block).
    for (int i = 0; i < 2; ++ i)
    {
        CharacterNodeToCell();
        CharacterCellToNode();
    }

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorIDSepMode(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->SetKeySearchByKeyActive();
        }
    }
}

void OversetConfigFactoryMulti::CalMinMaxBoxOfBlocks()
{
    //! Calculate minBox of each block.
    vector< vector< RDouble > > pMinOfBlocksPresent;
    vector< vector< RDouble > > pMaxOfBlocksPresent;

    pMinOfBlocksPresent.resize(numberOfBlocks);
    pMaxOfBlocksPresent.resize(numberOfBlocks);

    pMinOfBlocks.resize(numberOfBlocks);
    pMaxOfBlocks.resize(numberOfBlocks);

    for (int iBlock = 0; iBlock < numberOfBlocks; ++ iBlock)
    {
        pMinOfBlocksPresent[iBlock].resize(3);
        pMaxOfBlocksPresent[iBlock].resize(3);
        pMinOfBlocks[iBlock].resize(3);
        pMaxOfBlocks[iBlock].resize(3);

        SetField(pMinOfBlocksPresent[iBlock],  LARGE, 3);
        SetField(pMaxOfBlocksPresent[iBlock], -LARGE, 3);
        SetField(pMinOfBlocks[iBlock],  LARGE, 3);
        SetField(pMaxOfBlocks[iBlock], -LARGE, 3);
    }

    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorIDSepMode(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            int blockIndex = oversetConfigMultiZone[iZone]->blockIndex;

            oversetConfigMultiZone[iZone]->GetGrid()->ComputeMinMaxBox();

            RDouble *pMin = (oversetConfigMultiZone[iZone]->grid)->GetMinBox();
            RDouble *pMax = (oversetConfigMultiZone[iZone]->grid)->GetMaxBox();

            for (int i = 0; i < 3; ++ i)
            {
                pMinOfBlocksPresent[blockIndex][i] = MIN(pMinOfBlocksPresent[blockIndex][i], pMin[i]);
                pMaxOfBlocksPresent[blockIndex][i] = MAX(pMaxOfBlocksPresent[blockIndex][i], pMax[i]);
            }
        }
    }

    for (int iBlock = 0; iBlock < numberOfBlocks; ++ iBlock)
    {
        PH_AllreduceSepModeForVector(pMinOfBlocksPresent[iBlock], pMinOfBlocks[iBlock], 3, MPI_DOUBLE, PH_MIN);
        PH_AllreduceSepModeForVector(pMaxOfBlocksPresent[iBlock], pMaxOfBlocks[iBlock], 3, MPI_DOUBLE, PH_MAX);

        RDouble eps = GetToleranceForOversetBox();

        for (int i = 0; i < 3; ++ i)
        {
            pMinOfBlocks[iBlock][i] -= eps;
            pMaxOfBlocks[iBlock][i] += eps;
        }
    }
}

void OversetConfigFactoryMulti::CreatGlobalNodesByFile()
{
    WriteGlobalNodes();

    ReadGlobalNodes();
}

void OversetConfigFactoryMulti::WriteGlobalNodes()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    int serverProcessorIndex = PHMPI::GetServerSepMode();

    fstream file;

    if (currentProcessorIndex == serverProcessorIndex)
    {
        OpenFile(file, "globalNodes.dat", ios_base::out | ios_base::binary);
    }

    //! step1. Use process server to collect GlobalNodes information.
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        DataContainer *dataContainer = new DataContainer();

        int sendProcessorIndex = PHMPI::GetZoneProcessorIDSepMode(iZone);

        if (currentProcessorIndex == sendProcessorIndex)
        {
            oversetConfigMultiZone[iZone]->CreatGlobalNodes(dataContainer);
        }

        PHMPI::PH_Trade(dataContainer, sendProcessorIndex, serverProcessorIndex);

        if (currentProcessorIndex == serverProcessorIndex)
        {
            dataContainer->moveto_begin();
            dataContainer->writefile(file);
        }

        delete dataContainer;
    }

    if (currentProcessorIndex == serverProcessorIndex)
    {
        CloseFile(file);
    }
}

void OversetConfigFactoryMulti::ReadGlobalNodes()
{
    fstream file;

    //! All processes read in file.
    OpenFile(file, "globalNodes.dat", ios_base::in | ios_base::binary);

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        DataContainer *dataContainer = new DataContainer();

        int blockIndex;
        int numberOfSearchNodes;

        dataContainer->readfile(file);

        PHRead(dataContainer, blockIndex);
        PHRead(dataContainer, numberOfSearchNodes);

        globalZoneNodes[iZone]->Initialize(iZone, blockIndex, numberOfSearchNodes);

        vector< RDouble > &x = globalZoneNodes[iZone]->GetX();
        vector< RDouble > &y = globalZoneNodes[iZone]->GetY();
        vector< RDouble > &z = globalZoneNodes[iZone]->GetZ();
        vector< RDouble > &pMin = globalZoneNodes[iZone]->GetPMin();
        vector< RDouble > &pMax = globalZoneNodes[iZone]->GetPMax();

        if (numberOfSearchNodes != 0)
        {
            PHRead(dataContainer, x, numberOfSearchNodes);
            PHRead(dataContainer, y, numberOfSearchNodes);
            PHRead(dataContainer, z, numberOfSearchNodes);
        }

        PHRead(dataContainer, pMin, 3);
        PHRead(dataContainer, pMax, 3);

        delete dataContainer;
    }

    CloseFile(file);
}

void OversetConfigFactoryMulti::CreatGlobalNodes()
{
    int currentProcessor = PHMPI::GetCurrentProcessorID();

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int numberOfNodes = 0;
        int numberOfSearchNodes = 0;
        int blockIndex = 0;

        ZoneNodes *zoneNodes = globalZoneNodes[iZone];
        OversetConfig *oversetConfig = oversetConfigMultiZone[iZone];

        vector< RDouble > &x = zoneNodes->GetX();
        vector< RDouble > &y = zoneNodes->GetY();
        vector< RDouble > &z = zoneNodes->GetZ();
        vector< RDouble > &pMin = zoneNodes->GetPMin();
        vector< RDouble > &pMax = zoneNodes->GetPMax();

        //! Broadcast the numberOfNodes and blockIndex first, initialization zoneNodes.
        int sendProcessor = PHMPI::GetZoneProcessorIDSepMode(iZone);

        if (currentProcessor == sendProcessor)
        {
            UnstructGrid *grid = PHSPACE::UnstructGridCast(GetGrid(iZone, 0));
            vector< int > &keySearchOfNodes = oversetConfig->GetKeySearchOfNodes();

            numberOfNodes = grid->GetNTotalNode();
            blockIndex = grid->GetIBlock();

            for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
            {
                if (keySearchOfNodes[iNode] == 1)
                {
                    numberOfSearchNodes++;
                }
            }
        }

        PHMPI::PH_BcastSepMode(&numberOfSearchNodes, sizeof(int), sendProcessor, 1);

        PHMPI::PH_BcastSepMode(&blockIndex, sizeof(int), sendProcessor, 1);

        zoneNodes->Initialize(iZone, blockIndex, numberOfSearchNodes);

        //! Then broadcast x/y/z.
        if (currentProcessor == sendProcessor)
        {
            UnstructGrid *grid = PHSPACE::UnstructGridCast(GetGrid(iZone, 0));
            RDouble *x_zone = grid->GetX();
            RDouble *y_zone = grid->GetY();
            RDouble *z_zone = grid->GetZ();

            vector< RDouble > &pMin_zone = oversetConfig->GetPMin();
            vector< RDouble > &pMax_zone = oversetConfig->GetPMax();

            numberOfSearchNodes = 0;

            vector< int > &keySearchOfNodes = oversetConfig->GetKeySearchOfNodes();

            for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
            {
                if (keySearchOfNodes[iNode] == 1)
                {
                    x[numberOfSearchNodes] = x_zone[iNode];
                    y[numberOfSearchNodes] = y_zone[iNode];
                    z[numberOfSearchNodes] = z_zone[iNode];

                    numberOfSearchNodes++;
                }
            }

            for (int i = 0; i < 3; ++ i)
            {
                pMin[i] = pMin_zone[i];
                pMax[i] = pMax_zone[i];
            }
        }

        if (numberOfSearchNodes != 0)
        {
            PHMPI::PH_BcastSepMode(&x[0], numberOfSearchNodes * sizeof(RDouble), sendProcessor);
            PHMPI::PH_BcastSepMode(&y[0], numberOfSearchNodes * sizeof(RDouble), sendProcessor);
            PHMPI::PH_BcastSepMode(&z[0], numberOfSearchNodes * sizeof(RDouble), sendProcessor);
        }

        PHMPI::PH_BcastSepMode(&pMin[0], 3 * sizeof(RDouble), sendProcessor);
        PHMPI::PH_BcastSepMode(&pMax[0], 3 * sizeof(RDouble), sendProcessor);
    }
}

//! Collect node coordinate information, form globalNodes in each process. Packaged bcast with low efficiency.
void OversetConfigFactoryMulti::CreatGlobalNodes_Old()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        DataContainer *dataContainer = new DataContainer();

        int processorIndex = PHMPI::GetZoneProcessorIDSepMode(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            UnstructGrid *grid = PHSPACE::UnstructGridCast(GetGrid(iZone, 0));

            int numberOfNodes = grid->GetNTotalNode();
            int blockIndex = grid->GetIBlock();

            RDouble *x = grid->GetX();
            RDouble *y = grid->GetY();
            RDouble *z = grid->GetZ();

            vector< RDouble > &pMin = oversetConfigMultiZone[iZone]->GetPMin();
            vector< RDouble > &pMax = oversetConfigMultiZone[iZone]->GetPMax();

            //! Send x\y\z to master node.
            PHSPACE::PHWrite(dataContainer, blockIndex);
            PHSPACE::PHWrite(dataContainer, numberOfNodes);
            PHSPACE::PHWrite(dataContainer, &x[0], numberOfNodes);
            PHSPACE::PHWrite(dataContainer, &y[0], numberOfNodes);
            PHSPACE::PHWrite(dataContainer, &z[0], numberOfNodes);
            PHSPACE::PHWrite(dataContainer, &pMin[0], 3);
            PHSPACE::PHWrite(dataContainer, &pMax[0], 3);
        }

        int tag = 0;

        //! Can we use more scientific and reasonable broadcasting or information communication mode for large data?
        PHMPI::PH_BcastSepMode(dataContainer, processorIndex, tag);

        dataContainer->moveto_begin();

        int numberOfNodes;

        int blockIndex;

        ZoneNodes *zoneNodes = new ZoneNodes;

        globalZoneNodes[iZone] = zoneNodes;

        PHSPACE::PHRead(dataContainer, blockIndex);

        PHSPACE::PHRead(dataContainer, numberOfNodes);

        zoneNodes->Initialize(iZone, blockIndex, numberOfNodes);

        vector< RDouble > &x = zoneNodes->GetX();
        vector< RDouble > &y = zoneNodes->GetY();
        vector< RDouble > &z = zoneNodes->GetZ();

        vector< RDouble > &pMin = zoneNodes->GetPMin();
        vector< RDouble > &pMax = zoneNodes->GetPMax();

        PHSPACE::PHRead(dataContainer, &x[0], numberOfNodes);
        PHSPACE::PHRead(dataContainer, &y[0], numberOfNodes);
        PHSPACE::PHRead(dataContainer, &z[0], numberOfNodes);

        PHSPACE::PHRead(dataContainer, &pMin[0], 3);
        PHSPACE::PHRead(dataContainer, &pMax[0], 3);

        delete dataContainer;
    }
}

void OversetConfigFactoryMulti::SpecifyTheDonorCellAndMinDistanceOfGlobalNodes()
{
    int serverProcessorIndex = PHMPI::GetServerSepMode();
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();

    //! Set initial value.
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        ZoneNodes *zoneNodes = globalZoneNodes[iZone];
        vector< RDouble > &minDist = zoneNodes->GetMinDist();
        int numberOfNodes = zoneNodes->GetNumberOfNodes();

        SetField(minDist, LARGE, numberOfNodes);
    }

    RDouble maximumSearchTimeForSingleZone = 0.0;
    //! Query in all localZone.
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int numberOfSearchNodes = 0;
        RDouble timeForSearch = 0.0;

        int processorIndex = PHMPI::GetZoneProcessorIDSepMode(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->SearchDonorCellAndComputeDistanceForGlobalNodes(numberOfSearchNodes, timeForSearch);
        }

        PHMPI::PH_Trade(&numberOfSearchNodes, sizeof(int), processorIndex, serverProcessorIndex, iZone);
        PHMPI::PH_Trade(&timeForSearch, sizeof(RDouble), processorIndex, serverProcessorIndex, iZone);

        maximumSearchTimeForSingleZone = MAX(maximumSearchTimeForSingleZone, timeForSearch);
    }

    //! Communication, globalNodes start to Reduce.
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        ZoneNodes *zoneNodes = globalZoneNodes[iZone];
        int numberOfNodes = zoneNodes->GetNumberOfNodes();

        if (numberOfNodes == 0)
        {
            continue;
        }

        vector< RDouble > &minDist = zoneNodes->GetMinDist();
        vector< RDouble > minDistAllProcessor;
        minDistAllProcessor.resize(numberOfNodes);

        PH_AllreduceSepModeForVector(minDist, minDistAllProcessor, numberOfNodes, MPI_DOUBLE, PH_MIN);

        SetField(&minDist[0], &minDistAllProcessor[0], numberOfNodes);
    }

    //! Assign the value of globalZone->minDist to overset->minDist.
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        if (currentProcessorIndex == PHMPI::GetZoneProcessorIDSepMode(iZone))
        {
            oversetConfigMultiZone[iZone]->CopyMinDist();
        }
    }
}

void OversetConfigFactoryMulti::SpecifyInitialInterNodesByDistance()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorIDSepMode(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->SpecifyInitialInterNodesByDistance();
        }
    }
}

void OversetConfigFactoryMulti::SpecifySpecicalCellsInteractWithWallBoundary()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorIDSepMode(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->SpecifySpecicalCellsInteractWithWallBoundary();
        }
    }
}

void OversetConfigFactoryMulti::FillTheActiveRegion()
{
    int numberOfNewNodesPresentCycle = 0;
    int sum = 0;

    SetNodeCharacterOnPhysicalBoundary();

    do 
    {
        numberOfNewNodesPresentCycle = CharacterNodeToValidNodeByCell();

        //! mpi_reduce.
        PH_AllreduceSepMode(&numberOfNewNodesPresentCycle, &sum, 1, MPI_INT, MPI_SUM);
        numberOfNewNodesPresentCycle = sum;

    } while (numberOfNewNodesPresentCycle != 0);
}

void OversetConfigFactoryMulti::SetNodeCharacterOnPhysicalBoundary()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorIDSepMode(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->SetNodeCharacterOnPhysicalBoundary();
        }
    }
}

void OversetConfigFactoryMulti::SetNegativeForBufferRegion()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorIDSepMode(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->SetNegativeForBufferRegion();
        }
    }
}

void OversetConfigFactoryMulti::CreatGlobalInterCells()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int blockIndex = 0;
        int numberOfInterCells = 0;

        InterCells *interCells = globalInterCells[iZone];

        vector< int > interCellIndex;
        vector< RDouble > x;
        vector< RDouble > y;
        vector< RDouble > z;

        int processorIndex = PHMPI::GetZoneProcessorIDSepMode(iZone);

        if (currentProcessorIndex == processorIndex)
        {
            blockIndex = oversetConfigMultiZone[iZone]->GetBlockIndex();
            numberOfInterCells = oversetConfigMultiZone[iZone]->GetNumberOfInterCells();
        }

        PHMPI::PH_BcastSepMode(&blockIndex, sizeof(int), processorIndex);
        PHMPI::PH_BcastSepMode(&numberOfInterCells, sizeof(int), processorIndex);
        interCells->Initialize(iZone, blockIndex, numberOfInterCells);

        if (numberOfInterCells == 0)
        {
            continue;
        }

        if (currentProcessorIndex == processorIndex)
        {
            interCellIndex = oversetConfigMultiZone[iZone]->GetInterCellIndex();
            oversetConfigMultiZone[iZone]->CalInterCellInformation(interCellIndex, x, y, z);
        }

        if (numberOfInterCells != 0)
        {
            PHMPI::PH_BcastSepMode(&interCellIndex[0], numberOfInterCells * sizeof(int), processorIndex);
            PHMPI::PH_BcastSepMode(&x[0], numberOfInterCells * sizeof(RDouble), processorIndex);
            PHMPI::PH_BcastSepMode(&y[0], numberOfInterCells * sizeof(RDouble), processorIndex);
            PHMPI::PH_BcastSepMode(&z[0], numberOfInterCells * sizeof(RDouble), processorIndex);
        }

        interCells->SetInterCellIndex(interCellIndex);
        interCells->SetCellCenterX(x);
        interCells->SetCellCenterY(y);
        interCells->SetCellCenterZ(z);
    }
}

void OversetConfigFactoryMulti::CreatGlobalInterCellsServerCollectAndBcastNonBlocking()
{
    int currentProcessor = PHMPI::GetCurrentProcessorID();

    int serverProcessor = PHMPI::GetServerSepMode();

    DataContainer *dataForSend = new DataContainer();

    int lengthForSend;

    int nZonesInThisProc = 0;

    //! Obtain the nZonesInThisProc number of grid blocks in each process.
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorIDSepMode(iZone);

        if (currentProcessor == processorIndex)
        {
            ++nZonesInThisProc;
        }
    }

    PHSPACE::PHWrite(dataForSend, nZonesInThisProc);

    //! Obtain the inerCell information of all zones and write the dataContainer for each process.
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorIDSepMode(iZone);

        int blockIndex = 0;
        int numberOfInterCells = 0;

        vector< int > interCellIndex;
        vector< RDouble > x;
        vector< RDouble > y;
        vector< RDouble > z;

        if (currentProcessor == processorIndex)
        {
            blockIndex = oversetConfigMultiZone[iZone]->GetBlockIndex();
            numberOfInterCells = oversetConfigMultiZone[iZone]->GetNumberOfInterCells();

            interCellIndex = oversetConfigMultiZone[iZone]->GetInterCellIndex();
            oversetConfigMultiZone[iZone]->CalInterCellInformation(interCellIndex, x, y, z);

            PHSPACE::PHWrite(dataForSend, iZone);
            PHSPACE::PHWrite(dataForSend, blockIndex);
            PHSPACE::PHWrite(dataForSend, numberOfInterCells);

            if (numberOfInterCells != 0)
            {
                PHSPACE::PHWrite(dataForSend, interCellIndex, numberOfInterCells);
                PHSPACE::PHWrite(dataForSend, x, numberOfInterCells);
                PHSPACE::PHWrite(dataForSend, y, numberOfInterCells);
                PHSPACE::PHWrite(dataForSend, z, numberOfInterCells);
            }
        }
    }

    lengthForSend = static_cast<int>(dataForSend->size());

    //! Main process collects information.
    vector< DataContainer * > dataForReceive;
    vector< int > lengthForReceive;

    vector< PH_Request > requestContainerForData;
    vector< PH_Request > requestContainerForLength;

    int numberOfOriginalProcessor = PHMPI::GetNumberOfProcessor();
    int numberOfGridProcessor = PHMPI::GetNumberOfGridProcessor();

    int startProcessor = 0;
    int endProcessor = numberOfOriginalProcessor;
    int numberOfProcessor = numberOfOriginalProcessor;

    if (PHMPI::CurrentProcessorIsGridProcessor())
    {
        startProcessor = numberOfOriginalProcessor;
        endProcessor = numberOfOriginalProcessor + numberOfGridProcessor;
        numberOfProcessor = numberOfGridProcessor;
    }

    //! Collect lengthForReceive.
    if (currentProcessor == serverProcessor)
    {
        lengthForReceive.resize(numberOfProcessor);
        lengthForReceive[serverProcessor - startProcessor] = lengthForSend;
    }

    for (int iProcessor = startProcessor; iProcessor < endProcessor; ++ iProcessor)
    {
        if (iProcessor == serverProcessor)
        {
            continue;
        }

        if (currentProcessor == iProcessor)
        {
            requestContainerForLength.push_back(MPI_REQUEST_NULL);
            PH_Send(&lengthForSend, sizeof(int), serverProcessor, iProcessor, &requestContainerForLength.back());
        }

        if (currentProcessor == serverProcessor)
        {
            requestContainerForLength.push_back(MPI_REQUEST_NULL);
            PH_Receive(&lengthForReceive[iProcessor - startProcessor], sizeof(int), iProcessor, iProcessor, &requestContainerForLength.back());
        }
    }

    if (numberOfProcessor > 1)
    {
        PH_Wait(static_cast<int>(requestContainerForLength.size()), &(requestContainerForLength[0]));
    }

    //! Collect dataForReceive.
    if (currentProcessor == serverProcessor)
    {
        dataForReceive.resize(numberOfProcessor);
        dataForReceive[serverProcessor - startProcessor] = dataForSend;
    }

    for (int iProcessor = startProcessor; iProcessor < endProcessor; ++ iProcessor)
    {
        if (iProcessor == serverProcessor)
        {
            continue;
        }

        if (currentProcessor == iProcessor)
        {
            send(dataForSend, serverProcessor, lengthForSend, requestContainerForData, iProcessor);
        }

        if (currentProcessor == serverProcessor)
        {
            int nlen = lengthForReceive[iProcessor - startProcessor];

            dataForReceive[iProcessor - startProcessor] = new DataContainer;

            receive(dataForReceive[iProcessor - startProcessor], iProcessor, nlen, requestContainerForData, iProcessor);
        }
    }

    if (numberOfProcessor > 1)
    {
        PH_Wait(static_cast<int>(requestContainerForData.size()), &(requestContainerForData[0]));
    }

    DataContainer *dataForBcast = new DataContainer;

    if (currentProcessor == serverProcessor)
    {
        for (unsigned iData = 0; iData < dataForReceive.size(); ++ iData)
        {
            dataForBcast->write(dataForReceive[iData]);
        }
    }

    int tag = 0;
    PHMPI::PH_BcastByServerSepMode(dataForBcast, tag);

    dataForBcast->moveto_begin();
    for (int iData = 0; iData < numberOfProcessor; ++ iData)
    {
        int nZones = 0;
        PHRead(dataForBcast, nZones);

        for (int iZone = 0; iZone < nZones; ++ iZone)
        {
            InitializeInterCellByDataContainer(dataForBcast);
        }
    }

    if (currentProcessor != serverProcessor)
    {
        delete dataForSend;
    }

    if (currentProcessor == serverProcessor)
    {
        for (int iProcessor = 0; iProcessor < numberOfProcessor; ++ iProcessor)
        {
            delete dataForReceive[iProcessor];
        }
    }

    delete dataForBcast;
}

void OversetConfigFactoryMulti::CreatGlobalInterCellsServerCollectAndBcast()
{
    CreatGlobalInterCellsServerCollect();

    CreatGlobalInterCellsServerBcast();
}

void OversetConfigFactoryMulti::CreatGlobalInterCellsServerCollect()
{
    int currentProcessor = PHMPI::GetCurrentProcessorID();

    int serverProcessor = PHMPI::GetServerSepMode();

    DataContainer *dataContainer = new DataContainer();

    int nZonesInThisProc = 0;

    //! Obtain the nZonesInThisProc number of grid blocks in each process.
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorIDSepMode(iZone);

        if (currentProcessor == processorIndex)
        {
            ++ nZonesInThisProc;
        }
    }

    PHSPACE::PHWrite(dataContainer, nZonesInThisProc);

    //! Obtain the inerCell information of all zones and write the dataContainer for each process.
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorIDSepMode(iZone);

        int blockIndex = 0;
        int numberOfInterCells = 0;

        vector< int > interCellIndex;
        vector< RDouble > x;
        vector< RDouble > y;
        vector< RDouble > z;

        if (currentProcessor == processorIndex)
        {
            blockIndex = oversetConfigMultiZone[iZone]->GetBlockIndex();
            numberOfInterCells = oversetConfigMultiZone[iZone]->GetNumberOfInterCells();

            interCellIndex = oversetConfigMultiZone[iZone]->GetInterCellIndex();
            oversetConfigMultiZone[iZone]->CalInterCellInformation(interCellIndex, x, y, z);

            PHSPACE::PHWrite(dataContainer, iZone);
            PHSPACE::PHWrite(dataContainer, blockIndex);
            PHSPACE::PHWrite(dataContainer, numberOfInterCells);

            if (numberOfInterCells != 0)
            {
                PHSPACE::PHWrite(dataContainer, interCellIndex, numberOfInterCells);
                PHSPACE::PHWrite(dataContainer, x, numberOfInterCells);
                PHSPACE::PHWrite(dataContainer, y, numberOfInterCells);
                PHSPACE::PHWrite(dataContainer, z, numberOfInterCells);
            }
        }
    }

    //! Main process collects information and interprets it.
    int numberOfOriginalprocessor = PHMPI::GetNumberOfProcessor();
    int numberOfGridprocessor = PHMPI::GetNumberOfGridProcessor();

    int startProcessor = 0;
    int endProcessor = numberOfOriginalprocessor;

    if (PHMPI::CurrentProcessorIsGridProcessor())
    {
        startProcessor = numberOfOriginalprocessor;
        endProcessor = numberOfOriginalprocessor + numberOfGridprocessor;
    }

    if (currentProcessor != serverProcessor)
    {
        send(dataContainer, serverProcessor, currentProcessor);
    }

    if (currentProcessor == serverProcessor)
    {
        dataContainer->moveto_begin();

        PHSPACE::PHRead(dataContainer, nZonesInThisProc);

        for (int iZone = 0; iZone < nZonesInThisProc; ++ iZone)
        {
            InitializeInterCellByDataContainer(dataContainer);
        }

        for (int iProcessor = startProcessor; iProcessor < endProcessor; ++ iProcessor)
        {
            if (iProcessor == serverProcessor)
            {
                continue;
            }

            delete dataContainer;

            dataContainer = new DataContainer();

            receive(dataContainer, iProcessor, iProcessor);

            dataContainer->moveto_begin();

            PHSPACE::PHRead(dataContainer, nZonesInThisProc);

            for (int iZone = 0; iZone < nZonesInThisProc; ++ iZone)
            {
                InitializeInterCellByDataContainer(dataContainer);
            }
        }
    }

    delete dataContainer;
}

void OversetConfigFactoryMulti::InitializeInterCellByDataContainer(DataContainer *dataContainer)
{
    int zoneIndex;
    int blockIndex;
    int numberOfInterCells;

    PHSPACE::PHRead(dataContainer, zoneIndex);
    PHSPACE::PHRead(dataContainer, blockIndex);
    PHSPACE::PHRead(dataContainer, numberOfInterCells);

    globalInterCells[zoneIndex]->Initialize(zoneIndex, blockIndex, numberOfInterCells);

    vector< int > interCellIndex;
    vector< RDouble > x;
    vector< RDouble > y;
    vector< RDouble > z;

    interCellIndex.resize(numberOfInterCells);
    x.resize(numberOfInterCells);
    y.resize(numberOfInterCells);
    z.resize(numberOfInterCells);

    if (numberOfInterCells != 0)
    {

        PHSPACE::PHRead(dataContainer, interCellIndex, numberOfInterCells);
        PHSPACE::PHRead(dataContainer, x, numberOfInterCells);
        PHSPACE::PHRead(dataContainer, y, numberOfInterCells);
        PHSPACE::PHRead(dataContainer, z, numberOfInterCells);
    }

    globalInterCells[zoneIndex]->SetInterCellIndex(interCellIndex);
    globalInterCells[zoneIndex]->SetCellCenterX(x);
    globalInterCells[zoneIndex]->SetCellCenterY(y);
    globalInterCells[zoneIndex]->SetCellCenterZ(z);
}

void OversetConfigFactoryMulti::CreatGlobalInterCellsServerBcast()
{
    int currentProcessor = PHMPI::GetCurrentProcessorID();

    int serverProcessor = PHMPI::GetServerSepMode();

    DataContainer *dataContainer = new DataContainer();

    if (currentProcessor == serverProcessor)
    {
        for (int iZone = 0; iZone < numberOfZones; ++ iZone)
        {
            int blockIndex = globalInterCells[iZone]->GetBlockIndex();
            int numberOfInterCells = globalInterCells[iZone]->GetNumberOfInterCells();

            PHSPACE::PHWrite(dataContainer, iZone);
            PHSPACE::PHWrite(dataContainer, blockIndex);
            PHSPACE::PHWrite(dataContainer, numberOfInterCells);

            vector< int > &interCellIndex = globalInterCells[iZone]->GetInterCellIndex();
            vector< RDouble > &x = globalInterCells[iZone]->GetX();
            vector< RDouble > &y = globalInterCells[iZone]->GetY();
            vector< RDouble > &z = globalInterCells[iZone]->GetZ();

            if (numberOfInterCells != 0)
            {
                PHSPACE::PHWrite(dataContainer, interCellIndex, numberOfInterCells);
                PHSPACE::PHWrite(dataContainer, x, numberOfInterCells);
                PHSPACE::PHWrite(dataContainer, y, numberOfInterCells);
                PHSPACE::PHWrite(dataContainer, z, numberOfInterCells);
            }
        }
    }

    PHMPI::PH_BcastByServerSepMode(dataContainer, 0);

    if (currentProcessor != serverProcessor)
    {
        dataContainer->moveto_begin();

        for (int iZone = 0; iZone < numberOfZones; ++ iZone)
        {
            InitializeInterCellByDataContainer(dataContainer);
        }
    }

    delete dataContainer;
}

void OversetConfigFactoryMulti::FindDonorZoneAndDonorCell()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();

    //! Set initial value.
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        InterCells *interCells = globalInterCells[iZone];

        vector< RDouble > &minDist = interCells->GetMinDist();

        vector< int > &key = interCells->GetKey();

        int numberOfInterCells = interCells->GetNumberOfInterCells();

        SetField(minDist, LARGE, numberOfInterCells);

        SetField(key, -3, numberOfInterCells);
    }

    //! Query in all localZone.
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorIDSepMode(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->FindDonorZoneAndDonorCell();
        }
    }
}

void OversetConfigFactoryMulti::SelectDonorZoneAndDonorCell()
{
    int numberOfTotalInterCells;

    vector< int >     numberOfInterCellsEachZone;
    vector< int >     keyLocal;
    vector< RDouble > minDistLocal;
    vector< int >     donorZoneLocal;
    vector< int >     donorCellLocal;

    //! Collect information.
    numberOfTotalInterCells = 0;

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        InterCells *interCells = globalInterCells[iZone];
        int numberOfInterCells = interCells->GetNumberOfInterCells();

        numberOfInterCellsEachZone.push_back(numberOfInterCells);

        if (numberOfInterCells == 0)
        {
            continue;
        }

        vector< int > &key = interCells->GetKey();
        vector< RDouble > &minDist = interCells->GetMinDist();
        vector< int > &donorZone = interCells->GetDonorZoneIndex();
        vector< int > &donorCell = interCells->GetDonorCellIndex();

        for (int iCell = 0; iCell < numberOfInterCells; ++ iCell)
        {
            keyLocal.push_back(key[iCell]);
            minDistLocal.push_back(minDist[iCell]);
            donorZoneLocal.push_back(donorZone[iCell]);
            donorCellLocal.push_back(donorCell[iCell]);
        }

        numberOfTotalInterCells += numberOfInterCells;
    }

    vector< int >     keyGlobal;
    vector< RDouble > minDistGlobal;
    vector< int >     donorZoneGlobal;
    vector< int >     donorCellGlobal;

    keyGlobal.resize(numberOfTotalInterCells);
    minDistGlobal.resize(numberOfTotalInterCells);
    donorZoneGlobal.resize(numberOfTotalInterCells);
    donorCellGlobal.resize(numberOfTotalInterCells);

    //! Judge by key first, if the key of this process is less than the global key, then the value of this process will be discarded.
    PH_AllreduceSepModeForVector(keyLocal, keyGlobal, numberOfTotalInterCells, MPI_INT, PH_MAX);

    for (int iCell = 0; iCell < numberOfTotalInterCells; ++ iCell)
    {
        if (keyLocal[iCell] < keyGlobal[iCell])
        {
            minDistLocal[iCell] = LARGE;
            donorZoneLocal[iCell] = -1;
        }
    }

    //! Then judge according to the nearest distance.
    PH_AllreduceSepModeForVector(minDistLocal, minDistGlobal, numberOfTotalInterCells, MPI_DOUBLE, PH_MIN);

    //! donorZone determined by the nearest distance, if there are more than one same nearest distance, the zone with large label is automatically selected by MPIAllreduce.
    for (int iCell = 0; iCell < numberOfTotalInterCells; ++ iCell)
    {
        if ((minDistLocal[iCell] - minDistGlobal[iCell]) > TINY)
        {
            donorZoneLocal[iCell] = -1;
        }
    }

    //! Go through the above two steps, The unique donorZone and donorCell are determined.
    //! Finally, donorCell is determined by donoeZone.
    PH_AllreduceSepModeForVector(donorZoneLocal, donorZoneGlobal, numberOfTotalInterCells, MPI_INT, PH_MAX);

    for (int iCell = 0; iCell < numberOfTotalInterCells; ++ iCell)
    {
        if (donorZoneLocal[iCell] != donorZoneGlobal[iCell])
        {
            if (donorCellLocal[iCell] >= 0)
            {
                donorCellLocal[iCell] = -1;
            }
        }
    }

    PH_AllreduceSepModeForVector(donorCellLocal, donorCellGlobal, numberOfTotalInterCells, MPI_INT, PH_MAX);

    numberOfTotalInterCells = 0;

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        InterCells *interCells = globalInterCells[iZone];
        int numberOfInterCells = interCells->GetNumberOfInterCells();

        if (numberOfInterCells == 0)
        {
            continue;
        }

        vector< int > &key = interCells->GetKey();
        vector< RDouble > &minDist = interCells->GetMinDist();
        vector< int > &donorZone = interCells->GetDonorZoneIndex();
        vector< int > &donorCell = interCells->GetDonorCellIndex();

        for (int iCell = 0; iCell < numberOfInterCells; ++ iCell)
        {
            key[iCell] = keyGlobal[numberOfTotalInterCells + iCell];
            minDist[iCell] = minDistGlobal[numberOfTotalInterCells + iCell];
            donorZone[iCell] = donorZoneGlobal[numberOfTotalInterCells + iCell];
            donorCell[iCell] = donorCellGlobal[numberOfTotalInterCells + iCell];
        }

        numberOfTotalInterCells += numberOfInterCells;
    }
}

void OversetConfigFactoryMulti::CreatGlobalInterCells_Old()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();

    globalInterCells = new InterCells * [numberOfZones];

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        DataContainer *dataContainer = new DataContainer();

        int processorIndex = PHMPI::GetZoneProcessorIDSepMode(iZone);

        if (currentProcessorIndex == processorIndex)
        {
            UnstructGrid *grid = PHSPACE::UnstructGridCast(GetGrid(iZone, 0));

            int* keyActiveOfCells = grid->GetBlankIndex();

            int blockIndex = grid->GetIBlock();
            int numberOfCells = grid->GetNTotalCell();

            RDouble *cellCenterXContainer = grid->GetCellCenterX();
            RDouble *cellCenterYContainer = grid->GetCellCenterY();
            RDouble *cellCenterZContainer = grid->GetCellCenterZ();

            int numberOfInterCells = 0;

            vector< int > interCellIndex;

            vector< RDouble > x;
            vector< RDouble > y;
            vector< RDouble > z;

            for (int iCell = 0; iCell < numberOfCells; ++ iCell)
            {
                if (keyActiveOfCells[iCell] == 0)
                {
                    interCellIndex.push_back(iCell);

                    x.push_back(cellCenterXContainer[iCell]);
                    y.push_back(cellCenterYContainer[iCell]);
                    z.push_back(cellCenterZContainer[iCell]);

                    numberOfInterCells++;
                }
            }

            PHSPACE::PHWrite(dataContainer, blockIndex);
            PHSPACE::PHWrite(dataContainer, numberOfInterCells);

            if (numberOfInterCells != 0)
            {
                PHSPACE::PHWrite(dataContainer, interCellIndex, numberOfInterCells);
                PHSPACE::PHWrite(dataContainer, x, numberOfInterCells);
                PHSPACE::PHWrite(dataContainer, y, numberOfInterCells);
                PHSPACE::PHWrite(dataContainer, z, numberOfInterCells);
            }
        }

        int tag = 0;

        //! Can we use more scientific and reasonable broadcasting or information communication mode for large data?
        PHMPI::PH_BcastSepMode(dataContainer, processorIndex, tag);

        dataContainer->moveto_begin();

        int blockIndex;
        int numberOfInterCells;
        PHSPACE::PHRead(dataContainer, blockIndex);
        PHSPACE::PHRead(dataContainer, numberOfInterCells);

        InterCells *interCells = new InterCells;
        globalInterCells[iZone] = interCells;
        interCells->Initialize(iZone, blockIndex, numberOfInterCells);

        if (numberOfInterCells == 0)
        {
            continue;
        }

        vector< int > &interCellIndex = interCells->GetInterCellIndex();
        vector< RDouble > &x = interCells->GetX();
        vector< RDouble > &y = interCells->GetY();
        vector< RDouble > &z = interCells->GetZ();

        interCellIndex.resize(numberOfInterCells);
        x.resize(numberOfInterCells);
        y.resize(numberOfInterCells);
        z.resize(numberOfInterCells);

        PHSPACE::PHRead(dataContainer, &interCellIndex[0], numberOfInterCells);
        PHSPACE::PHRead(dataContainer, &x[0], numberOfInterCells);
        PHSPACE::PHRead(dataContainer, &y[0], numberOfInterCells);
        PHSPACE::PHRead(dataContainer, &z[0], numberOfInterCells);

        delete dataContainer;
    }
}

//! Identify active cells.
void OversetConfigFactoryMulti::FindActiveCells()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    //! Interpolation cells are determined.
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorIDSepMode(iZone);
        if (processorIndex == currentProcessorIndex)
        {
            oversetConfigMultiZone[iZone]->FindActiveCells();
        }
    }
}

//! Identify interpolation cells.
void OversetConfigFactoryMulti::FindInterCells()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();

    //! Interpolation cells are determined.

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int num = 0;
        int processorIndex = PHMPI::GetZoneProcessorIDSepMode(iZone);
        if (processorIndex == currentProcessorIndex)
        {
            num = oversetConfigMultiZone[iZone]->FindInterCells();
        }
    }
}

bool OversetConfigFactoryMulti::CheckDonorZoneAndDonorCell()
{
    bool okOrNot = true;

    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();
    int serverProcessorIndex = PHMPI::GetServerSepMode();

    if (currentProcessorIndex == serverProcessorIndex)
    {
        okOrNot = CheckDonorZoneAndDonorCellServer();
    }

    PHMPI::PH_BcastSepMode(&okOrNot, sizeof(bool), serverProcessorIndex);

    return okOrNot;
}

bool OversetConfigFactoryMulti::CheckDonorZoneAndDonorCellServer()
{
    int numberOfNegativeCellsAll_1 = 0;
    int numberOfNegativeCellsAll_2 = 0;
    int numberOfNegativeCellsAll_3 = 0;

    for (int iBlock = 0; iBlock < numberOfBlocks; ++ iBlock)
    {
        int numberOfNegativeCellsBlock_1 = 0;
        int numberOfNegativeCellsBlock_2 = 0;
        int numberOfNegativeCellsBlock_3 = 0;
        int numberOfInterCellsBlock = 0;

        for (int iZone = 0; iZone < numberOfZones; ++ iZone)
        {
            InterCells *interCells = globalInterCells[iZone];

            if (interCells->GetBlockIndex() != iBlock)
            {
                continue;
            }

            int numberOfInterCells = interCells->GetNumberOfInterCells();

            vector< int > &key = interCells->GetKey();
            vector< int > &donorCell = interCells->GetDonorCellIndex();
            vector< RDouble > &xc = interCells->GetX();
            vector< RDouble > &yc = interCells->GetY();
            vector< RDouble > &zc = interCells->GetZ();

            for (int iCell = 0; iCell < numberOfInterCells; ++ iCell)
            {
                if (key[iCell] < 0 || donorCell[iCell] < 0)
                {
                    if (key[iCell] == -3)
                    {
                        numberOfNegativeCellsBlock_3++;
                        cout << "key:" << key[iCell] << "  iCell:" << iCell << "  donorCell:" << donorCell[iCell] << endl;
                        cout << "        " << xc[iCell] << " " << yc[iCell] << " " << zc[iCell] << endl;
                    }
                    else if (key[iCell] == -2)
                    {
                        numberOfNegativeCellsBlock_2++;
                    }
                    else if (key[iCell] == -1)
                    {
                        numberOfNegativeCellsBlock_1++;
                    }
                }
            }
            numberOfInterCellsBlock += numberOfInterCells;
        }

        numberOfNegativeCellsAll_1 += numberOfNegativeCellsBlock_1;
        numberOfNegativeCellsAll_2 += numberOfNegativeCellsBlock_2;
        numberOfNegativeCellsAll_3 += numberOfNegativeCellsBlock_3;

        if (numberOfNegativeCellsBlock_1 > 0)
        {
            cout << "       Number of negative cells(-1) = " << numberOfNegativeCellsBlock_1 << endl;
        }

        if (numberOfNegativeCellsBlock_2 > 0)
        {
            cout << "       Number of negative cells(-2) = " << numberOfNegativeCellsBlock_2 << endl;
        }

        if (numberOfNegativeCellsBlock_3 > 0)
        {
            cout << "       Number of negative cells(-3) = " << numberOfNegativeCellsBlock_3 << endl;
        }
    }

    if (numberOfNegativeCellsAll_2 > 0 || numberOfNegativeCellsAll_3 > 0)
    {
        return false;
    }
    else
    {
        return true;
    }
}

void OversetConfigFactoryMulti::ResetOversetInformation()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();

    //! Initialize oversetInformation on each zone.
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorIDSepMode(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->InitializeOversetInformation();
        }
    }

    //! Get oversetInterCells and oversetDonorCells at MainProcess, where oversetDonorCells should swap with neighbor zone.
    GetInterAndDonorCellManagersNew();

    //! PostProcessOversetInformation, get storageIndexContainer (Interpolation or contribution cell for a neighbor zone, its serial number in the whole interpolation or contribution cell).
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorIDSepMode(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->PostProcessOversetInformation();
        }
    }
}

void OversetConfigFactoryMulti::GetInterAndDonorCellManagers()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int sendProcessorIndex = PHMPI::GetZoneProcessorIDSepMode(iZone);

        //! First get neighborZoneIndexContainer and broadcast to all process.
        PHVectorInt1D neighborZoneIndexContainer;
        if (currentProcessorIndex == sendProcessorIndex)
        {
            OversetInformationProxy *oversetInformationProxy = oversetConfigMultiZone[iZone]->GetOversetInformationProxy();
            oversetInformationProxy->GetNeighborZoneIndexContainer(neighborZoneIndexContainer);
        }

        int numberOfNeighbors = static_cast<int>(neighborZoneIndexContainer.size());
        PHMPI::PH_BcastSepMode(&numberOfNeighbors, sizeof(int), sendProcessorIndex);

        if (numberOfNeighbors != 0)
        {
            neighborZoneIndexContainer.resize(numberOfNeighbors);
            PHMPI::PH_BcastSepMode(&neighborZoneIndexContainer[0], numberOfNeighbors * sizeof(int), sendProcessorIndex);
        }

        //! Loop query all neighbors one by one.
        for (int iNeighbor = 0; iNeighbor < numberOfNeighbors; ++ iNeighbor)
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

            int numberOfDonorCells = static_cast<int>(donorCellIndexContainer.size());
            PHMPI::PH_Trade(&numberOfDonorCells, sizeof(int), sendProcessorIndex, receiveProcessorIndex, iZone);

            int *index = new int [numberOfDonorCells];
            RDouble *xCell = new RDouble [numberOfDonorCells];
            RDouble *yCell = new RDouble [numberOfDonorCells];
            RDouble *zCell = new RDouble [numberOfDonorCells];

            donorCellIndexContainer.resize(numberOfDonorCells);
            interCellCenterX.resize(numberOfDonorCells);
            interCellCenterY.resize(numberOfDonorCells);
            interCellCenterZ.resize(numberOfDonorCells);

            if (currentProcessorIndex == sendProcessorIndex)
            {
                for (int iCell = 0; iCell < numberOfDonorCells; ++ iCell)
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
                for (int iCell = 0; iCell < numberOfDonorCells; ++ iCell)
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

            delete [] index;
            delete [] xCell;
            delete [] yCell;
            delete [] zCell;
        }
    }
}

void OversetConfigFactoryMulti::GetInterAndDonorCellManagersNew()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int sendProcessorIndex = PHMPI::GetZoneProcessorID(iZone);

        //! First get neighborZoneIndexContainer and broadcast to all process.
        PHVectorInt1D neighborZoneIndexContainer;

        GetNeighborZoneIndexContainerForTargetZone(iZone, neighborZoneIndexContainer);

        int numberOfNeighbors = static_cast<int>(neighborZoneIndexContainer.size());

        //! Loop query all neighbors one by one.
        for (int iNeighbor = 0; iNeighbor < numberOfNeighbors; ++ iNeighbor)
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

                numberOfDonorCell = static_cast<int>(donorCellIndexContainer.size());

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
    int numberOfInterCells = globalInterCells[iZone]->GetNumberOfInterCells();

    vector< int > &donorZone = globalInterCells[iZone]->GetDonorZoneIndex();

    set< int > neighborZoneIndexSet;

    for (int iCell = 0; iCell < numberOfInterCells; ++ iCell)
    {
        neighborZoneIndexSet.insert(donorZone[iCell]);
    }

    for (set< int >::iterator iter = neighborZoneIndexSet.begin(); iter != neighborZoneIndexSet.end(); ++ iter)
    {
        neighborZoneList.push_back(*iter);
    }
}

void OversetConfigFactoryMulti::OutPutOversetFiles()
{
    OutPutBoundaryCells();

    PrintCellInformation();

    int outTecplotOverset = GlobalDataBase::GetIntParaFromDB("outTecplotOverset");

    if (outTecplotOverset == 1)
    {
        OutPutTecplot();
    }

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

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
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

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorIDSepMode(iZone);

        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->ShowGrid();
        }
    }
}

void OversetConfigFactoryMulti::OutPutBoundaryCells()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorIDSepMode(iZone);

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

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
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

            int numberOfCells = (oversetConfigMultiZone[iZone]->grid)->GetNTotalCell();

            for (int iCell = 0; iCell < numberOfCells; ++ iCell)
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
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorIDSepMode(iZone);
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

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorIDSepMode(iZone);

        if (currentProcessorIndex == processorIndex)
        {
            int numberOfNewNodes = oversetConfigMultiZone[iZone]->CharacterNodeToValidNodeByCell();
            numberOfNewNodesAllZones += numberOfNewNodes;
        }
    }

    if (GetReadInSklFileOrNot() == 1)
    {
        //! Accurate communication of point attributes through skl file.
        CommunicateNodeCharacterOnInterfaceBySKL(0);
    }
    else
    {
        //! Transmission of point attributes through interface communication. InnerAndInitialInterNodes needs to be constrained.
        CommunicateNodeCharacterOnInterface();
        SetNegativeForInnerAndInitialInterNodes();
    }

    return numberOfNewNodesAllZones;
}

void OversetConfigFactoryMulti::CharacterCellToNode()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorIDSepMode(iZone);
        if (currentProcessorIndex == processorIndex)
        {
            oversetConfigMultiZone[iZone]->CharacterCellToNode();
        }
    }

    if (GetReadInSklFileOrNot() == 1)
    {
        //! Accurate communication of point attributes through skl file.
        CommunicateNodeCharacterOnInterfaceBySKL(0);
    }
    else
    {
        //! Transmission of point attributes through interface communication.
        CommunicateNodeCharacterOnInterface();
    }
}

void OversetConfigFactoryMulti::CommunicateNodeCharacterOnInterfaceBySKL(int key)
{
    if (numberOfBlocks == numberOfZones)
    {
        return;
    }

    for (int iBlock = 0; iBlock < numberOfBlocks; ++ iBlock)
    {
        int numberOfGlobalInterNodes = static_cast<int>(keyActiveOfGlobalInterNodes[iBlock].size());
        SetField(keyActiveOfGlobalInterNodes[iBlock], 0, numberOfGlobalInterNodes);
    }

    //! step1: Every zone update keyActiveOfGlobalNodes (Summary keyActiveOfNodes).
    int currentProcessor = PHMPI::GetCurrentProcessorID();
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorIDSepMode(iZone);
        if (currentProcessor == processorIndex)
        {
            oversetConfigMultiZone[iZone]->UpdateInterfaceNodeCharacter(keyActiveOfGlobalInterNodes, key);
        }
    }

    //! step2: mpiAllreduce.
    for (int iBlock = 0; iBlock < numberOfBlocks; ++ iBlock)
    {
        vector< int > num;

        int numberOfGlobalInterNodes = static_cast<int>(keyActiveOfGlobalInterNodes[iBlock].size());

        num.resize(numberOfGlobalInterNodes);

        PH_AllreduceSepMode(&keyActiveOfGlobalInterNodes[iBlock][0], &num[0], numberOfGlobalInterNodes, MPI_INT, PH_MAX);

        SetField(&keyActiveOfGlobalInterNodes[iBlock][0], &num[0], numberOfGlobalInterNodes);
    }

    //! step3: Assign the value of after mpiAllreduce to keyActiveOfNodes.
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorIDSepMode(iZone);
        if (currentProcessor == processorIndex)
        {
            oversetConfigMultiZone[iZone]->DownInterfaceNodeCharacter(keyActiveOfGlobalInterNodes, key);
        }
    }
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
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorIDSepMode(iZone);
        if (currentProcessor == processorIndex)
        {
            oversetConfigMultiZone[iZone]->CharacterInterfaceToNode();
        }
    }
}

void OversetConfigFactoryMulti::OutInterCellInformation()
{
    int currentProcessorIndex = PHMPI::GetCurrentProcessorID();

    string fileName = AddSymbolToFileName("interCells.dat", currentProcessorIndex);
    fstream file;
    PHSPACE::OpenFile(file, fileName, ios_base::out);

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        InterCells *interCells = globalInterCells[iZone];
        int numberOfInterCells = interCells->GetNumberOfInterCells();

        vector< RDouble > &minDist = interCells->GetMinDist();
        vector< int > &donorZone = interCells->GetDonorZoneIndex();
        vector< int > &donorCell = interCells->GetDonorCellIndex();

        file << "iZone = " << iZone << endl;
        for (int iCell = 0; iCell < numberOfInterCells; ++ iCell)
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
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int processorIndex = PHMPI::GetZoneProcessorIDSepMode(iZone);
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

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int blockIndex = globalInterCells[iZone]->GetBlockIndex();
        int numberOfInterCells = globalInterCells[iZone]->GetNumberOfInterCells();

        vector< int > &cellIndex = globalInterCells[iZone]->GetInterCellIndex();
        vector< RDouble > &x = globalInterCells[iZone]->GetX();
        vector< RDouble > &y = globalInterCells[iZone]->GetX();
        vector< RDouble > &z = globalInterCells[iZone]->GetX();

        file << "blockIndex         = " << blockIndex << endl;
        file << "numberOfInterCells = " << numberOfInterCells << endl;

        for (int iCell = 0; iCell < numberOfInterCells; ++ iCell)
        {
            file << cellIndex[iCell] << "  ";
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

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        int blockIndex = globalZoneNodes[iZone]->GetBlockIndex();
        int numberOfNodes = globalZoneNodes[iZone]->GetNumberOfNodes();

        vector< RDouble > &x = globalZoneNodes[iZone]->GetX();
        vector< RDouble > &y = globalZoneNodes[iZone]->GetX();
        vector< RDouble > &z = globalZoneNodes[iZone]->GetX();
        vector< RDouble > &pMin = globalZoneNodes[iZone]->GetPMin();
        vector< RDouble > &pMax = globalZoneNodes[iZone]->GetPMax();

        file << "blockIndex    = " << blockIndex << endl;
        file << "numberOfNodes = " << numberOfNodes << endl;

        file << pMin[0] << "  " << pMin[1] << "  " << pMin[2] << endl;

        file << pMax[0] << "  " << pMax[1] << "  " << pMax[2] << endl;

        for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
        {
            file << x[iNode] << "  " << y[iNode] << "  " << z[iNode] << endl;
        }
    }

    PHSPACE::CloseFile(file);
}

OversetAuxilaryGrid::OversetAuxilaryGrid()
{
    blockIndex = -1;
    grid = 0;
    adtSearch = 0;
    walldistCompute = 0;
}

OversetAuxilaryGrid::~OversetAuxilaryGrid()
{
    if (grid) delete grid;
    if (adtSearch) delete adtSearch;
}

void OversetAuxilaryGrid::Initialize(UnstructGrid *gridIn, int iBlock, Pre_WalldistCompute *walldistComputeIn)
{
    this->grid = gridIn;

    this->blockIndex = iBlock;

    this->walldistCompute = walldistComputeIn;
}

//! Because wallStrucPresentBlock isn't obtained from the auxiliary grid, it may lead to the object point distance not being 0.
//! Therefore, the object point distance should be set to 0.

void OversetAuxilaryGrid::BuildADT()
{
    if (adtSearch != 0)
    {
        delete adtSearch;
        adtSearch = 0;
    }

    adtSearch = new MRADTSearch(grid);

    adtSearch->BuildADTStructure();
}

void OversetAuxilaryGrid::CalMetrics()
{
    grid->ComputeMetrics();
}

void OversetAuxilaryGrid::CalGradientOfWallDistanceLSQ()
{

}

void OversetAuxilaryGrid::CalGradientOfWallDistanceGGNode()
{
    RDouble *faceNormalXContainer = grid->GetFaceNormalX();
    RDouble *faceNormalYContainer = grid->GetFaceNormalY();
    RDouble *faceNormalZContainer = grid->GetFaceNormalZ();

    RDouble *faceAreaContainer = grid->GetFaceArea();
    RDouble *cellVolumeContainer = grid->GetCellVolume();

    int *leftCellIndexContainer = grid->GetLeftCellOfFace();
    int *rightCellIndexContainer = grid->GetRightCellOfFace();

    int numberOfNodes = grid->GetNTotalNode();
    int numberOfCells = grid->GetNTotalCell();
    int numberOfFaces = grid->GetNTotalFace();
    int numberOfBoundaryFaces = grid->GetNBoundFace();
    int numberOfTotalCells = grid->GetNTotalCell();

    int *faceNodeNumberContainer = grid->GetNodeNumberOfEachFace();

    int *faceNodeIndexContainer = grid->GetFace2Node();

    RDouble *qNodeField = new RDouble [numberOfNodes];

    CompNodeVar(grid, qNodeField, &cellWallDistance[0]);

    wallDistanceGradientX.resize(numberOfTotalCells);
    wallDistanceGradientY.resize(numberOfTotalCells);
    wallDistanceGradientZ.resize(numberOfTotalCells);

    for (int iCell = 0; iCell < numberOfTotalCells; ++ iCell)
    {
        wallDistanceGradientX[iCell] = 0;
        wallDistanceGradientY[iCell] = 0;
        wallDistanceGradientZ[iCell] = 0;
    }

    int nodePosition = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        int leftCellIndex = leftCellIndexContainer[iFace];
        int rightCellIndex = rightCellIndexContainer[iFace];

        RDouble faceValue = 0.0;
        int faceNodeNumber = faceNodeNumberContainer[iFace];
        for (int iNodeInFace = 0; iNodeInFace < faceNodeNumber; ++ iNodeInFace)
        {
            int nodeIndex = faceNodeIndexContainer[nodePosition + iNodeInFace];
            faceValue += qNodeField[nodeIndex];
        }
        nodePosition += faceNodeNumber;
        faceValue /= faceNodeNumberContainer[iFace];

        RDouble faceAeraNormalX = faceNormalXContainer[iFace] * faceAreaContainer[iFace];
        RDouble faceAeraNormalY = faceNormalYContainer[iFace] * faceAreaContainer[iFace];
        RDouble faceAeraNormalZ = faceNormalZContainer[iFace] * faceAreaContainer[iFace];

        wallDistanceGradientX[leftCellIndex] += faceValue * faceAeraNormalX;
        wallDistanceGradientY[leftCellIndex] += faceValue * faceAeraNormalY;
        wallDistanceGradientZ[leftCellIndex] += faceValue * faceAeraNormalZ;

        if (iFace < numberOfBoundaryFaces) continue;
        wallDistanceGradientX[rightCellIndex] -= faceValue * faceAeraNormalX;
        wallDistanceGradientY[rightCellIndex] -= faceValue * faceAeraNormalY;
        wallDistanceGradientZ[rightCellIndex] -= faceValue * faceAeraNormalZ;
    }

    for (int iCell = 0; iCell < numberOfCells; ++ iCell)
    {
        wallDistanceGradientX[iCell] /= cellVolumeContainer[iCell];
        wallDistanceGradientY[iCell] /= cellVolumeContainer[iCell];
        wallDistanceGradientZ[iCell] /= cellVolumeContainer[iCell];
    }

    delete [] qNodeField;
}

void OversetAuxilaryGrid::CalGradientOfWallDistanceGGCell()
{
    RDouble *faceNormalXContainer = grid->GetFaceNormalX();
    RDouble *faceNormalYContainer = grid->GetFaceNormalY();
    RDouble *faceNormalZContainer = grid->GetFaceNormalZ();
            
    RDouble *faceCenterXContainer = grid->GetFaceCenterX();
    RDouble *faceCenterYContainer = grid->GetFaceCenterY();
    RDouble *faceCenterZContainer = grid->GetFaceCenterZ();
            
    RDouble *cellCenterXContainer = grid->GetCellCenterX();
    RDouble *cellCenterYContainer = grid->GetCellCenterY();
    RDouble *cellCenterZContainer = grid->GetCellCenterZ();
            
    RDouble *faceAreaContainer = grid->GetFaceArea();
    RDouble *cellVolumeContainer = grid->GetCellVolume();

    int *leftCellIndexContainer = grid->GetLeftCellOfFace();
    int *rightCellIndexContainer = grid->GetRightCellOfFace();

    int numberOfCells = grid->GetNTotalCell();
    int numberOfFaces = grid->GetNTotalFace();
    int numberOfBoundaryFaces = grid->GetNBoundFace();
    int numberOfTotalCells = grid->GetNTotalCell();

    wallDistanceGradientX.resize(numberOfTotalCells);
    wallDistanceGradientY.resize(numberOfTotalCells);
    wallDistanceGradientZ.resize(numberOfTotalCells);

    for (int iCell = 0; iCell < numberOfTotalCells; ++ iCell)
    {
        wallDistanceGradientX[iCell] = 0;
        wallDistanceGradientY[iCell] = 0;
        wallDistanceGradientZ[iCell] = 0;
    }

    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        int leftCellIndex = leftCellIndexContainer[iFace];
        int rightCellIndex = rightCellIndexContainer[iFace];

        RDouble dxl = faceCenterXContainer[iFace] - cellCenterXContainer[leftCellIndex];
        RDouble dyl = faceCenterYContainer[iFace] - cellCenterYContainer[leftCellIndex];
        RDouble dzl = faceCenterZContainer[iFace] - cellCenterZContainer[leftCellIndex];

        RDouble dxr = faceCenterXContainer[iFace] - cellCenterXContainer[rightCellIndex];
        RDouble dyr = faceCenterYContainer[iFace] - cellCenterYContainer[rightCellIndex];
        RDouble dzr = faceCenterZContainer[iFace] - cellCenterZContainer[rightCellIndex];

        RDouble delt1 = DISTANCE(dxl, dyl, dzl);
        RDouble delt2 = DISTANCE(dxr, dyr, dzr);
        RDouble delta = 1.0 / (delt1 + delt2 + SMALL);    //! distance weighted.
        RDouble cl = delt2 * delta;
        RDouble cr = delt1 * delta;

        RDouble faceValue = cl * cellWallDistance[leftCellIndex] + cr * cellWallDistance[rightCellIndex];

        RDouble faceAeraNormalX = faceNormalXContainer[iFace] * faceAreaContainer[iFace];
        RDouble faceAeraNormalY = faceNormalYContainer[iFace] * faceAreaContainer[iFace];
        RDouble faceAeraNormalZ = faceNormalZContainer[iFace] * faceAreaContainer[iFace];

        wallDistanceGradientX[leftCellIndex] += faceValue * faceAeraNormalX;
        wallDistanceGradientY[leftCellIndex] += faceValue * faceAeraNormalY;
        wallDistanceGradientZ[leftCellIndex] += faceValue * faceAeraNormalZ;

        if (iFace < numberOfBoundaryFaces) continue;
        wallDistanceGradientX[rightCellIndex] -= faceAeraNormalX * faceValue;
        wallDistanceGradientY[rightCellIndex] -= faceAeraNormalY * faceValue;
        wallDistanceGradientZ[rightCellIndex] -= faceAeraNormalZ * faceValue;
    }

    for (int iCell = 0; iCell < numberOfCells; ++ iCell)
    {
        RDouble reciprocalCellVolume = 1.0 / cellVolumeContainer[iCell];
        wallDistanceGradientX[iCell] *= reciprocalCellVolume;
        wallDistanceGradientY[iCell] *= reciprocalCellVolume;
        wallDistanceGradientZ[iCell] *= reciprocalCellVolume;
    }
}

void OversetAuxilaryGrid::CalCellWallDistance()
{
    int numberOfCells = grid->GetNTotalCell();
    int numberOfBoundaryFaces = grid->GetNBoundFace();
    int numberOfTotalCells = numberOfCells + numberOfBoundaryFaces;

    RDouble *x = grid->GetCellCenterX();
    RDouble *y = grid->GetCellCenterY();
    RDouble *z = grid->GetCellCenterZ();

    cellWallDistance.resize(numberOfTotalCells);

    SetField(&cellWallDistance[0], LARGE, numberOfTotalCells);

    walldistCompute->ComputeWallDistForNodesList(blockIndex, numberOfCells, x, y, z, &cellWallDistance[0]);

    CompositeUnstructBC **boundaryConditionRecord = grid->GetBCRecord();

    int *leftCellIndexContainer = grid->GetLeftCellOfFace();
    int *rightCellIndexContainer = grid->GetRightCellOfFace();

    for (int iFace = 0; iFace < numberOfBoundaryFaces; ++ iFace)
    {
        int boundaryConditionType = boundaryConditionRecord[iFace]->GetKey();

        int leftCellIndex = leftCellIndexContainer[iFace];
        int rightCellIndex = rightCellIndexContainer[iFace];

        if (boundaryConditionType == PHENGLEI::SOLID_SURFACE)
        {
            cellWallDistance[rightCellIndex] = -cellWallDistance[leftCellIndex];
        }
        else if (boundaryConditionType == PHENGLEI::SYMMETRY)
        {
            cellWallDistance[rightCellIndex] = cellWallDistance[leftCellIndex];
        }
        else
        {
            cellWallDistance[rightCellIndex] = cellWallDistance[leftCellIndex];
        }
    }
}

RDouble OversetAuxilaryGrid::SearchAndCalMinDistByVolumeCoordinate(RDouble xc, RDouble yc, RDouble zc)
{
    return 1.0;
}

RDouble OversetAuxilaryGrid::SearchAndCalMinDistByIDW(RDouble xc, RDouble yc, RDouble zc)
{
    return 1.0;
}

RDouble OversetAuxilaryGrid::SearchAndCalMinDistByGradient(RDouble xc, RDouble yc, RDouble zc)
{
    int donorCell = -1;

    adtSearch->SearchInAdt(xc, yc, zc, donorCell);

    if (donorCell < 0)
    {
        return LARGE;
    }

    RDouble *cellCenterXContainer = grid->GetCellCenterX();
    RDouble *cellCenterYContainer = grid->GetCellCenterY();
    RDouble *cellCenterZContainer = grid->GetCellCenterZ();

    RDouble dx = xc - cellCenterXContainer[donorCell];
    RDouble dy = yc - cellCenterYContainer[donorCell];
    RDouble dz = zc - cellCenterZContainer[donorCell];

    RDouble dqdx = wallDistanceGradientX[donorCell];
    RDouble dqdy = wallDistanceGradientY[donorCell];
    RDouble dqdz = wallDistanceGradientZ[donorCell];

    RDouble dist = cellWallDistance[donorCell] + dqdx * dx + dqdy * dy + dqdz * dz;

    return dist;
}

ZoneNodes::ZoneNodes()
{
    zoneIndex = 0;
    blockIndex = 0;
    numberOfNodes = 0;
}

ZoneNodes ::~ZoneNodes()
{

}

void ZoneNodes::Initialize(int iZone, int iBlock, int numberOfNodes)
{
    this->zoneIndex = iZone;
    this->blockIndex = iBlock;
    this->numberOfNodes = numberOfNodes;

    x.resize(numberOfNodes);
    y.resize(numberOfNodes);
    z.resize(numberOfNodes);

    minDist.resize(numberOfNodes);

    pMin.resize(3);
    pMax.resize(3);
}

InterCells::InterCells()
{

}

InterCells::~InterCells()
{
    zoneIndex = 0;
    blockIndex = 0;
    numberOfInterCells = 0;
}

void InterCells::Initialize(int iZone, int iBlock, int numberOfInterCells)
{
    this->zoneIndex = iZone;
    this->blockIndex = iBlock;
    this->numberOfInterCells = numberOfInterCells;

    donorZoneIndex.resize(numberOfInterCells);
    donorCellIndex.resize(numberOfInterCells);
    minDist.resize(numberOfInterCells);
    key.resize(numberOfInterCells);

    SetField(donorZoneIndex, -1, numberOfInterCells);
    SetField(donorCellIndex, -1, numberOfInterCells);
    SetField(key, -3, numberOfInterCells);
    SetField(minDist, static_cast <RDouble> (-1.0), numberOfInterCells);
}

OversetConfig::OversetConfig()
{
    numberOfZones = 0;
    blockIndex = 0;
    zoneIndex = 0;

    grid = 0;
    globalZoneNodes = 0;
    globalInterCells = 0;

    aDTSearch = 0;
    oversetInformationProxy = 0;

    pMin.resize(3);
    pMax.resize(3);

    initialInterNodeOrNot = 0;
    innerNodesOrNot = 0;
    innerCellsOrNot = 0;
    oversetConfigFactoryMulti = 0;
    numberOfInterCells = 0;
}

OversetConfig::~OversetConfig()
{
    if (aDTSearch != 0)
    {
        delete aDTSearch;
        aDTSearch = 0;
    }

    if (innerCellsOrNot != 0)
    {
        delete [] innerCellsOrNot;
    }

    if (innerNodesOrNot != 0)
    {
        delete [] innerNodesOrNot;
    }

    if (initialInterNodeOrNot != 0)
    {
        delete [] initialInterNodeOrNot;
    }
}

void OversetConfig::Initialize(OversetConfigFactoryMulti *oversetConfigFactoryMulti, UnstructGrid *gridIn, int nZone, int iZone, int iBlock, ZoneNodes **globalZoneNodesIn, InterCells **globalInterCellsIn)
{
    this->grid = gridIn;
    this->numberOfZones = nZone;
    this->zoneIndex = iZone;
    this->blockIndex = iBlock;
    this->globalZoneNodes = globalZoneNodesIn;
    this->globalInterCells = globalInterCellsIn;

    this->oversetConfigFactoryMulti = oversetConfigFactoryMulti;

    oversetInformationProxy = 0;

    Allocate();
}

void OversetConfig::Allocate()
{
    int numberOfNodes = grid->GetNTotalNode();
    int numberOfCells = grid->GetNTotalCell();
    int numberOfTotalCells = numberOfCells + grid->GetNBoundFace();

    keyActiveOfNodes.resize(numberOfNodes);
    SetField(keyActiveOfNodes, 0, numberOfNodes);

    int *keyActiveOfCells = this->GetKeyActiveOfCells();
    SetField(keyActiveOfCells, 0, numberOfTotalCells);

    keyBufferOfNodes.resize(numberOfNodes);
    SetField(keyBufferOfNodes, 0, numberOfNodes);

    keyBufferOfCells.resize(numberOfCells);
    SetField(keyBufferOfCells, 0, numberOfCells);

    keySearchOfNodes.resize(numberOfNodes);
    SetField(keySearchOfNodes, 0, numberOfNodes);

    keySearchOfCells.resize(numberOfCells);
    SetField(keySearchOfCells, 0, numberOfCells);

    minDist.resize(numberOfNodes);
    SetField(minDist, LARGE, numberOfNodes);

    initialInterNodeOrNot = new bool [numberOfNodes];
    SetField(initialInterNodeOrNot, false, numberOfNodes);

    innerCellsOrNot = new bool [numberOfCells];
    SetField(innerCellsOrNot, false, numberOfCells);

    innerNodesOrNot = new bool [numberOfNodes];
    SetField(innerNodesOrNot, false, numberOfNodes);
}

void OversetConfig::SetInitialValue()
{
    int numberOfNodes = grid->GetNTotalNode();
    int numberOfCells = grid->GetNTotalCell();
    int numberOfTotalCells = numberOfCells + grid->GetNBoundFace();

    SetField(keyActiveOfNodes, 0, numberOfNodes);

    int *keyActiveOfCells = this->GetKeyActiveOfCells();
    SetField(keyActiveOfCells, 0, numberOfTotalCells);

    SetField(keySearchOfNodes, 1, numberOfNodes);
    SetField(keySearchOfCells, 1, numberOfCells);

    SetField(minDist, LARGE, numberOfNodes);

    SetField(initialInterNodeOrNot, false, numberOfNodes);

    SetField(innerCellsOrNot, false, numberOfCells);

    SetField(innerNodesOrNot, false, numberOfNodes);
}

void OversetConfig::BuildADTStructure()
{
    //! ADTtree should be re-established after the mesh is moved, old ADTree needed to be delete firstly.
    if (aDTSearch != 0)
    {
        delete aDTSearch;
        aDTSearch = 0;
    }

    aDTSearch = new MRADTSearch(this->grid, this);

    aDTSearch->BuildADTStructure();
}

//! Calculate the minimum box according to the cell.
void OversetConfig::CalculateMinBox()
{
    int numberOfCells = grid->GetNTotalCell();

    int **cellNodeIndexContainer = grid->GetCell2NodeArray();
    int *cellNodeNumberContainer = grid->GetNodeNumberOfEachCell();

    RDouble *x = grid->GetX();
    RDouble *y = grid->GetY();
    RDouble *z = grid->GetZ();

    RDouble xmin, ymin, zmin, xmax, ymax, zmax;
    xmin = ymin = zmin = LARGE;
    xmax = ymax = zmax = -LARGE;

    for (int iCell = 0; iCell < numberOfCells; ++ iCell)
    {
        if (keySearchOfCells[iCell] != 1)
        {
            continue;
        }
        int cellNodeNumber = cellNodeNumberContainer[iCell];

        for (int iNode = 0; iNode < cellNodeNumber; ++ iNode)
        {
            int nodeIndex = cellNodeIndexContainer[iCell][iNode];

            xmin = PHSPACE::MIN(xmin, x[nodeIndex]);
            ymin = PHSPACE::MIN(ymin, y[nodeIndex]);
            zmin = PHSPACE::MIN(zmin, z[nodeIndex]);

            xmax = PHSPACE::MAX(xmax, x[nodeIndex]);
            ymax = PHSPACE::MAX(ymax, y[nodeIndex]);
            zmax = PHSPACE::MAX(zmax, z[nodeIndex]);
        }
    }

    RDouble eps = GetToleranceForOversetBox();

    pMin[0] = xmin - eps;
    pMin[1] = ymin - eps;
    pMin[2] = zmin - eps;

    pMax[0] = xmax + eps;
    pMax[1] = ymax + eps;
    pMax[2] = zmax + eps;
}

void OversetConfig::SearchPhysicalAndInterfaceBoundary()
{
    SearchBoundary();
}

void OversetConfig::SearchBoundary()
{
    int numberOfBoundaryFaces = grid->GetNBoundFace();

    CompositeUnstructBC **boundaryConditionRecord = grid->GetBCRecord();

    int *faceNodeNumber = grid->GetNodeNumberOfEachFace();
    int **faceNodeIndex = grid->GetFace2NodeArray();

    set< int > interfaceBoundaryNodeSet;
    set< int > wallBoundaryNodeSet;
    set< int > oversetBoundaryNodeSet;
    set< int > symmetryBoundaryNodeSet;
    set< int > physicalBoundaryNodeSet;

    int symetryOrNot = GetSymetryOrNot();

    for (int iFace = 0; iFace < numberOfBoundaryFaces; ++ iFace)
    {
        int boundaryType = boundaryConditionRecord[iFace]->GetKey();

        int numberOfNodesInFace = faceNodeNumber[iFace];

        switch (boundaryType)
        {
            case PHENGLEI::INTERFACE:
                for (int iNode = 0; iNode < numberOfNodesInFace; ++ iNode)
                {
                    interfaceBoundaryNodeSet.insert(faceNodeIndex[iFace][iNode]);
                }
                break;
            case PHENGLEI::SOLID_SURFACE:
                for (int iNode = 0; iNode < numberOfNodesInFace; ++ iNode)
                {
                    wallBoundaryNodeSet.insert(faceNodeIndex[iFace][iNode]);
                }
                break;
            case PHENGLEI::OVERSET:
                for (int iNode = 0; iNode < numberOfNodesInFace; ++ iNode)
                {
                    oversetBoundaryNodeSet.insert(faceNodeIndex[iFace][iNode]);
                }
                break;
            case PHENGLEI::SYMMETRY:
                for (int iNode = 0; iNode < numberOfNodesInFace; ++ iNode)
                {
                    symmetryBoundaryNodeSet.insert(faceNodeIndex[iFace][iNode]);
                }
                break;
            default:
                for (int iNode = 0; iNode < numberOfNodesInFace; ++ iNode)
                {
                    physicalBoundaryNodeSet.insert(faceNodeIndex[iFace][iNode]);
                }
        }
    }

    for (set< int >::iterator iter = interfaceBoundaryNodeSet.begin(); iter != interfaceBoundaryNodeSet.end(); ++ iter)
    {
        interfaceBoundaryNodesList.push_back(*iter);
    }
    for (set< int >::iterator iter = oversetBoundaryNodeSet.begin(); iter != oversetBoundaryNodeSet.end(); ++ iter)
    {
        innerBoundaryNodesList.push_back(*iter);
    }
    for (set< int >::iterator iter = physicalBoundaryNodeSet.begin(); iter != physicalBoundaryNodeSet.end(); ++ iter)
    {
        physicalBoundaryNodesList.push_back(*iter);
    }
    for (set< int >::iterator iter = wallBoundaryNodeSet.begin(); iter != wallBoundaryNodeSet.end(); ++ iter)
    {
        physicalBoundaryNodesList.push_back(*iter);
    }

    if (symetryOrNot)
    {
        for (set< int >::iterator iter = symmetryBoundaryNodeSet.begin(); iter != symmetryBoundaryNodeSet.end(); ++ iter)
        {
            physicalBoundaryNodesList.push_back(*iter);
        }
    }
}

void OversetConfig::SetInitialValueForBufferRigion()
{
    //! Because communicate_cell_iblank is the direct communication keyActiveOfCells.
    //! keyBuffer is determined by keyActiveOfCells array.
    int *keyActiveOfCells = this->GetKeyActiveOfCells();
    vector< int > &keyActiveOfNodes = this->GetKeyActiveOfNodes();

    int numberOfNodes = grid->GetNTotalNode();
    int numberOfCells = grid->GetNTotalCell();

    SetField(keyActiveOfNodes, 0, numberOfNodes);
    SetField(keyActiveOfCells, 0, numberOfCells);

    int numberOfInnerBoundaryNodes = innerBoundaryNodesList.size();
    for (int iNode = 0; iNode < numberOfInnerBoundaryNodes; ++ iNode)
    {
        int nodeIndex = innerBoundaryNodesList[iNode];

        keyActiveOfNodes[nodeIndex] = 1;
    }
}

void OversetConfig::SetKeyBufferByKeyActive()
{
    //! Because communicate_cell_iblank is the direct communication keyActiveOfCells.
    //! keyBuffer is determined by keyActiveOfCells array.
    int numberOfNodes = grid->GetNTotalNode();
    int numberOfCells = grid->GetNTotalCell();

    vector< int > &keyActiveOfNodes = this->GetKeyActiveOfNodes();
    int *keyActiveOfCells = this->GetKeyActiveOfCells();

    vector< int > &keyBufferOfCells = this->GetKeyBufferOfCells();
    vector< int > &keyBufferOfNodes = this->GetKeyBufferOfNodes();

    keyBufferOfCells.resize(numberOfCells);
    keyBufferOfNodes.resize(numberOfNodes);

    SetField(&keyBufferOfCells[0], &keyActiveOfCells[0], numberOfCells);
    SetField(&keyBufferOfNodes[0], &keyActiveOfNodes[0], numberOfNodes);

    SetField(keyActiveOfNodes, 0, numberOfNodes);
    SetField(keyActiveOfCells, 0, numberOfCells);
}

void OversetConfig::ComputeNodeAndCellTopology()
{
    int **cellNodeIndexContainer = grid->GetCell2NodeArray();
    int **cellFaceIndexContainer = grid->GetCell2Face();
}

void OversetConfig::CalMinDistanceByAuxilaryGrid(vector< OversetAuxilaryGrid* > &auxilaryGrid)
{
    int numberOfNodes = grid->GetNTotalNode();

    SetField(minDist, LARGE, numberOfNodes);

    RDouble *xNodes = grid->GetX();
    RDouble *yNodes = grid->GetY();
    RDouble *zNodes = grid->GetZ();

    for (unsigned iAuxilaryGrid = 0; iAuxilaryGrid < auxilaryGrid.size(); ++ iAuxilaryGrid)
    {
        if (iAuxilaryGrid == blockIndex)
        {
            continue;
        }

        for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
        {
            if (keySearchOfNodes[iNode] != 1)
            {
                continue;
            }

            RDouble xc = xNodes[iNode];
            RDouble yc = yNodes[iNode];
            RDouble zc = zNodes[iNode];
            RDouble dist;

            dist = auxilaryGrid[iAuxilaryGrid]->SearchAndCalMinDistByGradient(xc, yc, zc);

            if (minDist[iNode] > dist)
            {
                minDist[iNode] = dist;
            }
        }
    }
}

void OversetConfig::SearchDonorCellAndComputeDistanceForGlobalNodes(int &numberOfSearchNodes, RDouble &time)
{
    numberOfSearchNodes = 0;

    TimeSpan timeForSearch;

    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        ZoneNodes *zoneNodes = globalZoneNodes[iZone];
        vector< RDouble > &minDist = zoneNodes->GetMinDist();
        int numberOfNodes = zoneNodes->GetNumberOfNodes();

        //! Two zones belong to the same block, there is no need to query.
        if (blockIndex == zoneNodes->GetBlockIndex())
        {
            continue;
        }

        if (!TheTargetZoneOverLapWithMe(zoneNodes))
        {
            continue;
        }

        vector< RDouble > &xc = zoneNodes->GetX();
        vector< RDouble > &yc = zoneNodes->GetY();
        vector< RDouble > &zc = zoneNodes->GetZ();

        for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
        {
            RDouble x = xc[iNode];
            RDouble y = yc[iNode];
            RDouble z = zc[iNode];

            RDouble dist = SearchDonorCellAndComputeDistance(x, y, z);

            numberOfSearchNodes++;

            if (dist < minDist[iNode])
            {
                minDist[iNode] = dist;
            }
        }
    }

    time = timeForSearch.GetTimeSpan();
}

void OversetConfig::CopyMinDist()
{
    int numberOfNodes = grid->GetNTotalNode();

    vector< RDouble > &minDistInGlobalZoneNodes = globalZoneNodes[zoneIndex]->GetMinDist();

    int numberOfSearchNodes = 0;

    for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
    {
        if (keySearchOfNodes[iNode] == 1)
        {
            minDist[iNode] = minDistInGlobalZoneNodes[numberOfSearchNodes];
            numberOfSearchNodes++;
        }
        else
        {
            minDist[iNode] = LARGE;
        }
    }
}

RDouble OversetConfig::SearchDonorCellAndComputeDistance(RDouble xc, RDouble yc, RDouble zc)
{
    int donorCell = SearchDonorCell(xc, yc, zc);

    return ComputeDistanceInDonorZone(donorCell, xc, yc, zc);
}

bool OversetConfig::TheTargetZoneOverLapWithMe(ZoneNodes *zoneNodes)
{
    vector< RDouble > &pMin_targetZone = zoneNodes->GetPMin();
    vector< RDouble > &pMax_targetZone = zoneNodes->GetPMax();

    vector< RDouble > &pMin_presentZone = this->GetPMin();
    vector< RDouble > &pMax_presentZone = this->GetPMax();

    if (pMax_targetZone[0] < pMin_presentZone[0] ||
        pMin_targetZone[0] > pMax_presentZone[0] ||
        pMax_targetZone[1] < pMin_presentZone[1] ||
        pMin_targetZone[1] > pMax_presentZone[1] ||
        pMax_targetZone[2] < pMin_presentZone[2] ||
        pMin_targetZone[2] > pMax_presentZone[2])
    {
        return false;
    }

    return true;
}

int OversetConfig::SearchDonorCell(RDouble xc, RDouble yc, RDouble zc)
{
    int donorCell = -1;

    aDTSearch->SearchInAdt(xc, yc, zc, donorCell);

    //! If the contribution cell is a buffer cell, return -1.
    if (donorCell >= 0)
    {
        if (keyBufferOfCells[donorCell] == 1)
        {
            donorCell = -1;
        }
    }

    return donorCell;
}

void OversetConfig::SearchDonorCell(RDouble xc, RDouble yc, RDouble zc, int &donorCell)
{
    int numberOfCells = grid->GetNTotalCell();
    RDouble eps = GetToleranceForOversetSearch();

    int * leftCellIndexContainer = grid->GetLeftCellOfFace();
    int **cellFaceIndexContainer = grid->GetCell2Face();

    int *cellFaceNumberContainer = grid->GetFaceNumberOfEachCell();

    RDouble *faceNormalXContainer = grid->GetFaceNormalX();
    RDouble *faceNormalYContainer = grid->GetFaceNormalY();
    RDouble *faceNormalZContainer = grid->GetFaceNormalZ();

    RDouble *faceCenterXContainer = grid->GetFaceCenterX();
    RDouble *faceCenterYContainer = grid->GetFaceCenterY();
    RDouble *faceCenterZContainer = grid->GetFaceCenterZ();

    for (int iCell = 0; iCell < numberOfCells; ++ iCell)
    {
        int key = 1;

        int numberOfFacesInCell = cellFaceNumberContainer[iCell];

        for (int iFace = 0; iFace < numberOfFacesInCell; ++ iFace)
        {
            int faceIndex = cellFaceIndexContainer[iCell][iFace];

            RDouble xn = faceNormalXContainer[faceIndex];
            RDouble yn = faceNormalYContainer[faceIndex];
            RDouble zn = faceNormalZContainer[faceIndex];

            RDouble xf = faceCenterXContainer[faceIndex];
            RDouble yf = faceCenterYContainer[faceIndex];
            RDouble zf = faceCenterZContainer[faceIndex];

            int leftCellIndex = leftCellIndexContainer[faceIndex];

            RDouble volume = (xc - xf) * xn + (yc - yf) * yn + (zc - zf) * zn;

            if (leftCellIndex == iCell)
            {
                volume = -volume;
            }

            if (volume < -eps)
            {
                key = 0;
                break;
            }
        }

        if (key == 1)
        {
            donorCell = iCell;
        }
    }
}

RDouble OversetConfig::ComputeDistanceInDonorZone(int donorCell, RDouble xc, RDouble yc, RDouble zc)
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

    int numberOfNodeInCell = cellNodeNumberContainer[donorCell];

    RDouble weightDist = 0.0;
    RDouble weight = 0.0;

    vector< RDouble > &nodeWalldistDistance = grid->GetNodeWallDist();

    for (int iNode = 0; iNode < numberOfNodeInCell; ++ iNode)
    {
        int nodeIndex = cellNodeIndexContainer[donorCell][iNode];

        RDouble dx = xc - coordinateX[nodeIndex];
        RDouble dy = yc - coordinateY[nodeIndex];
        RDouble dz = zc - coordinateZ[nodeIndex];
        RDouble dist = nodeWalldistDistance[nodeIndex];

        weight += 1.0 / (DISTANCE(dx, dy, dz) + TINY);
        weightDist += 1.0 / (DISTANCE(dx, dy, dz) + TINY) * dist;
    }

    weightDist /= weight;

    return weightDist;
}

//! Internal point is set to 0.
//! Local distance is greater than neighbor distance, set to 0.
//! The others are set to 1.
void OversetConfig::SpecifyNodesCharacterByDistanceAndInnerNodesCheck()
{
    int numberOfNodes = grid->GetNTotalNode();

    ZoneNodes *zoneNodes = globalZoneNodes[zoneIndex];

    vector< RDouble > &minDistInNeighbourZone = zoneNodes->GetMinDist();

    SetField(keyActiveOfNodes, 1, numberOfNodes);

    int numberOfSearchNodes = 0;

    for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
    {
        if (keySearchOfNodes[iNode] == 1)
        {
            if (nodeWalldistDistance[iNode] >= minDistInNeighbourZone[numberOfSearchNodes])
            {
                keyActiveOfNodes[iNode] = 0;
            }
            numberOfSearchNodes++;
        }

        if (innerNodesOrNot[iNode])
        {
            keyActiveOfNodes[iNode] = 0;
        }
    }
}

//! Interpolation boundary is determined according to the local distance and neighbor distance.
//! Overset boundary can be ignored.
void OversetConfig::SpecifyInitialInterNodesByDistance()
{
    int numberOfNodes = grid->GetNTotalNode();

    vector< RDouble > &nodeWalldistDistance = grid->GetNodeWallDist();

    SetField(initialInterNodeOrNot, false, numberOfNodes);

    for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
    {
        if (nodeWalldistDistance[iNode] >= minDist[iNode])
        {
            initialInterNodeOrNot[iNode] = true;
        }
    }
}

void OversetConfig::SpecifyCellsInBody(vector< MRADTSearch* > &adtSearchArray)
{
    int numberOfNodes = grid->GetNTotalNode();
    int numberOfCells = grid->GetNTotalCell();
    int numberOfBlocks = oversetConfigFactoryMulti->GetNumberOfBlocks();

    int **cellNodeIndexContainer = grid->GetCell2NodeArray();
    int *cellNodeNumberContainer = grid->GetNodeNumberOfEachCell();

    RDouble *cellCenterX = grid->GetCellCenterX();
    RDouble *cellCenterY = grid->GetCellCenterY();
    RDouble *cellCenterZ = grid->GetCellCenterZ();

    SetField(innerCellsOrNot, false, numberOfCells);
    SetField(innerNodesOrNot, false, numberOfNodes);

    int symetryOrNot = GetSymetryOrNot();

    for (int iCell = 0; iCell < numberOfCells; ++ iCell)
    {
        bool negativeCell = false;

        RDouble xc = cellCenterX[iCell];
        RDouble yc = cellCenterY[iCell];
        RDouble zc = cellCenterZ[iCell];

        int donorCell = -1;

        for (int iBlock = 0; iBlock < numberOfBlocks; iBlock ++)
        {
            if (blockIndex == iBlock)
            {
                continue;
            }

            if (adtSearchArray[iBlock] != 0)
            {
                adtSearchArray[iBlock]->SearchInAdt(xc, yc, zc, donorCell);
            }
        }

        if ((symetryOrNot && zc < 0.0) || donorCell != -1)
        {
            negativeCell = true;
        }

        if (!negativeCell)
        {
            continue;
        }

        innerCellsOrNot[iCell] = true;

        for (int iNode = 0; iNode < cellNodeNumberContainer[iCell]; ++ iNode)
        {
            int nodeIndex = cellNodeIndexContainer[iCell][iNode];
            innerNodesOrNot[nodeIndex] = true;
        }
    }

    int numberOfInnerCells = 0;
    int numberOfInnerNodes = 0;

    for (int iCell = 0; iCell < numberOfCells; ++ iCell)
    {
        if (innerCellsOrNot[iCell])
        {
            numberOfInnerCells++;
        }
    }

    for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
    {
        if (innerNodesOrNot[iNode])
        {
            numberOfInnerNodes++;
        }
    }
}

//! Any two points of an cell pass through the physical plane, the cell is a boundary cell.
//! This method isn't applicable in the presence of background grid.
void OversetConfig::SpecifySpecicalCellsInteractWithWallBoundary()
{
    int numberOfCells = grid->GetNTotalCell();

    int **cellNodeIndexContainer = grid->GetCell2NodeArray();
    int *cellNodeNumberContainer = grid->GetNodeNumberOfEachCell();

    RDouble *xNode = grid->GetX();
    RDouble *yNode = grid->GetY();
    RDouble *zNode = grid->GetZ();

    for (int iCell = 0; iCell < numberOfCells; ++ iCell)
    {
        int node1 = cellNodeIndexContainer[iCell][0];

        int numberOfNodeInCell = cellNodeNumberContainer[iCell];

        for (int iNode = 1; iNode < numberOfNodeInCell; ++ iNode)
        {
            int node2 = cellNodeIndexContainer[iCell][iNode];

            RDouble dist1 = minDist[node1];
            RDouble dist2 = minDist[node2];

            RDouble dx = xNode[node1] - xNode[node2];
            RDouble dy = yNode[node1] - yNode[node2];
            RDouble dz = zNode[node1] - zNode[node2];

            RDouble length = DISTANCE(dx, dy, dz);

            if (dist1 < length || dist2 < length)
            {
                if (dist2 > 0.5 * LARGE || dist1 > 0.5 * LARGE)
                {
                    innerCellsOrNot[iCell] = true;
                    break;
                }
            }
        }
    }

    int numberOfInnerCells = 0;

    for (int iCell = 0; iCell < numberOfCells; ++ iCell)
    {
        if (innerCellsOrNot[iCell])
        {
            for (int iNode = 0; iNode < cellNodeNumberContainer[iCell]; ++ iNode)
            {
                int nodeIndex = cellNodeIndexContainer[iCell][iNode];
                innerNodesOrNot[nodeIndex] = true;
            }

            numberOfInnerCells++;
        }
    }
}

void OversetConfig::SetNegativeForBufferRegion()
{
    int numberOfNodes = grid->GetNTotalNode();
    for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
    {
        if (keyBufferOfNodes[iNode] == 1)
        {
            keyActiveOfNodes[iNode] = 0;
        }
    }
}

void OversetConfig::SetNegativeForInnerNodes()
{
    int numberOfNodes = grid->GetNTotalNode();

    for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
    {
        if (innerNodesOrNot[iNode])
        {
            keyActiveOfNodes[iNode] = 0;
        }
    }
}

void OversetConfig::SetNegativeForInnerAndInitialInterNodes()
{
    int numberOfNodes = grid->GetNTotalNode();

    for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
    {
        if (initialInterNodeOrNot[iNode] || innerNodesOrNot[iNode])
        {
            keyActiveOfNodes[iNode] = 0;
        }
    }
}

void OversetConfig::SetNodeCharacterOnPhysicalBoundary()
{
    vector< int > &keyActiveOfNodes = this->GetKeyActiveOfNodes();

    int numberOfNodes = grid->GetNTotalNode();

    SetField(keyActiveOfNodes, 0, numberOfNodes);

    int numberOfPhysicalBoundaryNodes = physicalBoundaryNodesList.size();
    for (int iNode = 0; iNode < numberOfPhysicalBoundaryNodes; ++ iNode)
    {
        keyActiveOfNodes[physicalBoundaryNodesList[iNode]] = 1;
    }
}

void OversetConfig::FillTheRemainingNodes()
{
    int numberOfNodes = grid->GetNTotalNode();

    for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
    {
        if (keyActiveOfNodes[iNode] == 0)
        {
            keyActiveOfNodes[iNode] = -1;
        }
    }
}

void OversetConfig::CreatGlobalNodes(DataContainer *dataContainer)
{
    RDouble *x_zone = grid->GetX();
    RDouble *y_zone = grid->GetY();
    RDouble *z_zone = grid->GetZ();

    vector< RDouble > x;
    vector< RDouble > y;
    vector< RDouble > z;

    int numberOfSearchNodes = 0;

    int numberOfNodes = grid->GetNTotalNode();

    for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
    {
        if (keySearchOfNodes[iNode] == 1)
        {
            x.push_back(x_zone[iNode]);
            y.push_back(y_zone[iNode]);
            z.push_back(z_zone[iNode]);

            numberOfSearchNodes++;
        }
    }

    PHSPACE::PHWrite(dataContainer, blockIndex);
    PHSPACE::PHWrite(dataContainer, numberOfSearchNodes);

    if (numberOfSearchNodes != 0)
    {
        PHSPACE::PHWrite(dataContainer, &x[0], numberOfSearchNodes);
        PHSPACE::PHWrite(dataContainer, &y[0], numberOfSearchNodes);
        PHSPACE::PHWrite(dataContainer, &z[0], numberOfSearchNodes);
    }

    PHSPACE::PHWrite(dataContainer, &pMin[0], 3);
    PHSPACE::PHWrite(dataContainer, &pMax[0], 3);
}

void OversetConfig::SetSearchRegionGlobal(vector< vector< RDouble > > &pMinOfBlocks, vector< vector< RDouble > > &pMaxOfBlocks)
{
    //! Because communicate_cell_iblank is the direct communication keyActiveOfCells.
    //! keyBuffer is determined by keyActiveOfCells array.
    int *keyActiveOfCells = this->GetKeyActiveOfCells();
    vector< int > &keyActiveOfNodes = this->GetKeyActiveOfNodes();

    int numberOfNodes = grid->GetNTotalNode();
    int numberOfCells = grid->GetNTotalCell();

    SetField(keyActiveOfNodes, 0, numberOfNodes);
    SetField(keyActiveOfCells, 0, numberOfCells);

    RDouble *x = grid->GetX();
    RDouble *y = grid->GetY();
    RDouble *z = grid->GetZ();

    for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
    {
        for (unsigned iBlock = 0; iBlock < pMinOfBlocks.size(); ++ iBlock)
        {
            if (iBlock == blockIndex)
            {
                continue;
            }

            if (x[iNode] > pMinOfBlocks[iBlock][0] && x[iNode] < pMaxOfBlocks[iBlock][0] &&
                y[iNode] > pMinOfBlocks[iBlock][1] && y[iNode] < pMaxOfBlocks[iBlock][1] &&
                z[iNode] > pMinOfBlocks[iBlock][2] && z[iNode] < pMaxOfBlocks[iBlock][2])
            {
                keyActiveOfNodes[iNode] = 1;
                break;
            }
        }
    }
}

void OversetConfig::SetKeySearchByKeyActive()
{
    //! Because communicate_cell_iblank is the direct communication keyActiveOfCells.
    //! keyBuffer is determined by keyActiveOfCells array.
    int numberOfNodes = grid->GetNTotalNode();
    int numberOfCells = grid->GetNTotalCell();

    vector< int > &keyActiveOfNodes = this->GetKeyActiveOfNodes();
    int* keyActiveOfCells = this->GetKeyActiveOfCells();

    vector< int > &keySearchOfCells = this->GetKeySearchOfCells();
    vector< int > &keySearchOfNodes = this->GetKeySearchOfNodes();

    SetField(&keySearchOfNodes[0], &keyActiveOfNodes[0], numberOfNodes);
    SetField(&keySearchOfCells[0], keyActiveOfCells, numberOfCells);

    SetField(keyActiveOfNodes, 0, numberOfNodes);
    SetField(keyActiveOfCells, 0, numberOfCells);
}

//! It is an active cell if the all cell nodes are 1.
void OversetConfig::FindActiveCells()
{
    int numberOfCells = grid->GetNTotalCell();

    int **cellNodeIndexContainer = grid->GetCell2NodeArray();
    int *cellNodeNumberContainer = grid->GetNodeNumberOfEachCell();

    int *keyActiveOfCells = this->GetKeyActiveOfCells();

    SetField(keyActiveOfCells, 1, numberOfCells);

    for (int iCell = 0; iCell < numberOfCells; ++ iCell)
    {
        int &numberOfNodeInCell = cellNodeNumberContainer[iCell];
        for (int iNode = 0; iNode < numberOfNodeInCell; ++ iNode)
        {
            int nodeIndex = cellNodeIndexContainer[iCell][iNode];
            if (keyActiveOfNodes[nodeIndex] != 1)
            {
                keyActiveOfCells[iCell] = 0;
                break;
            }
        }
    }
}

//! It is an interpolation cell if the cell node has both 1 and 0.
int OversetConfig::FindInterCells()
{
    int numberOfCells = grid->GetNTotalCell();

    int **cellNodeIndexContainer = grid->GetCell2NodeArray();
    int *cellNodeNumberContainer = grid->GetNodeNumberOfEachCell();

    int *keyActiveOfCells = this->GetKeyActiveOfCells();

    for (int iCell = 0; iCell < numberOfCells; ++ iCell)
    {
        if (keyActiveOfCells[iCell] == 1)
        {
            continue;
        }

        int numberOfNodeInCell = cellNodeNumberContainer[iCell];
        for (int iNode = 0; iNode < numberOfNodeInCell; ++ iNode)
        {
            int nodeIndex = cellNodeIndexContainer[iCell][iNode];
            if (keyActiveOfNodes[nodeIndex] == 1)
            {
                keyActiveOfCells[iCell] = -1;
                break;
            }
        }
    }

    numberOfInterCells = 0;
    interCellIndex.resize(0);

    for (int iCell = 0; iCell < numberOfCells; ++ iCell)
    {
        if (keyActiveOfCells[iCell] == -1)
        {
            interCellIndex.push_back(iCell);
            numberOfInterCells++;

            if (innerCellsOrNot[iCell])
            {
                cout << "Cell in body CANNOT be Interpolation cell!!! " << " iCell = " << iCell << endl;

                int numberOfNodeInCell = cellNodeNumberContainer[iCell];
                for (int iNode = 0; iNode < numberOfNodeInCell; ++ iNode)
                {
                    int nodeIndex = cellNodeIndexContainer[iCell][iNode];
                    cout << "   iNode = " << iNode << "innerOrNot =" << innerNodesOrNot[nodeIndex]  << endl;
                    cout << "   iNode = " << iNode << "keyActive  =" << keyActiveOfNodes[nodeIndex] << endl;
                }
                exit(0);
            }
        }
    }

    return numberOfInterCells;
}

void OversetConfig::CalInterCellInformation(vector< int > &index, vector< RDouble > &x, vector< RDouble > &y, vector< RDouble > &z)
{

    RDouble *cellCenterXContainer = grid->GetCellCenterX();
    RDouble *cellCenterYContainer = grid->GetCellCenterY();
    RDouble *cellCenterZContainer = grid->GetCellCenterZ();

    x.resize(numberOfInterCells);
    y.resize(numberOfInterCells);
    z.resize(numberOfInterCells);

    for (int iCell = 0; iCell < numberOfInterCells; ++ iCell)
    {
        int cellIndex = interCellIndex[iCell];

        x[iCell] = cellCenterXContainer[cellIndex];
        y[iCell] = cellCenterYContainer[cellIndex];
        z[iCell] = cellCenterZContainer[cellIndex];
    }
}

int OversetConfig::FindDonorZoneAndDonorCellOfInterCell(int &donorCell, RDouble xc, RDouble yc, RDouble zc)
{
    //! If there's no contribution cell, return -3.
    //! If contribution cell is inactive, return -2.
    //! If contribution cell is interpolation, return -1.

    donorCell = -1;
    int donorCell2 = -1;

    aDTSearch->SearchInAdt(xc, yc, zc, donorCell2);

    donorCell = donorCell2;

    int *keyActiveOfCells = GetKeyActiveOfCells();

    if (donorCell == -1)
    {
        return -3;
    }

    if (keyActiveOfCells[donorCell] == 0)
    {
        return -2;
    }

    if (keyActiveOfCells[donorCell] == -1)
    {
        return -1;
    }

    return 1;
}

void OversetConfig::FindDonorZoneAndDonorCell()
{
    for (int iZone = 0; iZone < numberOfZones; ++ iZone)
    {
        InterCells *interCells = globalInterCells[iZone];
        vector< RDouble > &minDist = interCells->GetMinDist();
        vector< int > &maxKey = interCells->GetKey();

        vector< int > &donorZone = interCells->GetDonorZoneIndex();
        vector< int > &donorCell = interCells->GetDonorCellIndex();

        int numberOfInterCells = interCells->GetNumberOfInterCells();

        //! Two zones belong to the same block, there is no need to query.
        if (blockIndex == interCells->GetBlockIndex())
        {
            continue;
        }

        vector< RDouble > &xc = interCells->GetX();
        vector< RDouble > &yc = interCells->GetY();
        vector< RDouble > &zc = interCells->GetZ();

        for (int iCell = 0; iCell < numberOfInterCells; ++ iCell)
        {
            RDouble x = xc[iCell];
            RDouble y = yc[iCell];
            RDouble z = zc[iCell];

            int targetCell = -1;

            int key = FindDonorZoneAndDonorCellOfInterCell(targetCell, x, y, z);

            if (key == -3)
            {
                continue;
            }

            RDouble dist = ComputeDistanceInDonorZone(targetCell, x, y, z);

            if (key > maxKey[iCell])
            {
                minDist[iCell] = dist;
                donorZone[iCell] = zoneIndex;
                donorCell[iCell] = targetCell;
                maxKey[iCell] = key;
            }
            else if (key == maxKey[iCell])
            {
                if (dist < minDist[iCell])
                {
                    minDist[iCell] = dist;
                    donorZone[iCell] = zoneIndex;
                    donorCell[iCell] = targetCell;
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
    PHVector1D < OversetCell* > *oversetCellContainer = new PHVector1D < OversetCell* >();

    int numberOfInterCells = globalInterCells[zoneIndex]->GetNumberOfInterCells();

    vector< int > &interCellIndex = globalInterCells[zoneIndex]->GetInterCellIndex();
    vector< int > &donorCellIndex = globalInterCells[zoneIndex]->GetDonorCellIndex();
    vector< int > &donorZoneIndex = globalInterCells[zoneIndex]->GetDonorZoneIndex();

    for (int iCell = 0; iCell < numberOfInterCells; ++ iCell)
    {
        int interCell = interCellIndex[iCell];
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

    int numberOfInterNodes = interfaceBoundaryNodesList.size();
    for (int iNode = 0; iNode < numberOfInterNodes; ++ iNode)
    {
        int localNodeIndex = interfaceBoundaryNodesList[iNode];
        int globalNodeIndex = localToGlobalNodesIndex[localNodeIndex];
        keyNodesGlobal[globalNodeIndex] = 1;
    }
}

void OversetConfig::CalLocalToGlobalInterfaceNodes(int iBlock, vector< int > &keyNodesGlobal)
{
    if (blockIndex != iBlock) return;

    int numberOfInterNodes = interfaceBoundaryNodesList.size();
    localToGlobalInterNodesIndex.resize(numberOfInterNodes);

    for (int iNode = 0; iNode < numberOfInterNodes; ++ iNode)
    {
        int localNodeIndex = interfaceBoundaryNodesList[iNode];
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

    int **cellNodeIndexContainer = grid->GetCell2NodeArray();
    int *cellNodeNumberContainer = grid->GetNodeNumberOfEachCell();

    int *keyActiveOfCells = this->GetKeyActiveOfCells();
    SetField(keyActiveOfCells, 0, numberOfCells);

    for (int iCell = 0; iCell < numberOfCells; ++ iCell)
    {
        int numberOfNodeInCell = cellNodeNumberContainer[iCell];
        for (int iNode = 0; iNode < numberOfNodeInCell; ++ iNode)
        {
            int nodeIndex = cellNodeIndexContainer[iCell][iNode];
            if (keyActiveOfNodes[nodeIndex] == 1)
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

    CompositeUnstructBC **boundaryConditionRecord = grid->GetBCRecord();

    int *faceNodeNumber = grid->GetNodeNumberOfEachFace();
    int *faceNodeIndex = grid->GetFace2Node();

    int *leftCellIndex = grid->GetLeftCellOfFace();
    int *rightCellIndex = grid->GetRightCellOfFace();

    int *keyActiveOfCells = grid->GetBlankIndex();

    int iCount = 0;
    for (int iFace = 0; iFace < numberOfBoundaryFaces; ++ iFace)
    {
        int boundaryType = boundaryConditionRecord[iFace]->GetKey();

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

            for (int iNode = 0; iNode < numberOfNodesInFace; ++ iNode)
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

    SetField(keyActiveOfNodes, 0, numberOfNodes);

    for (int iCell = 0; iCell < numberOfCells; ++ iCell)
    {
        if (keyActiveOfCells[iCell] == 1)
        {
            int numberOfNodeInCell = cellNodeNumberContainer[iCell];
            for (int iNode = 0; iNode < numberOfNodeInCell; ++ iNode)
            {
                int nodeIndex = cellNodeIndexContainer[iCell][iNode];
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

    SetField(keyActiveOfCells, 0, numberOfCells);

    int numberOfNewNodes = 0;

    for (int iCell = 0; iCell < numberOfCells; ++ iCell)
    {
        int cellNodeNumber = cellNodeNumberContainer[iCell];

        for (int iNode = 0; iNode < cellNodeNumber; ++ iNode)
        {
            int nodeIndex = cellNodeIndexContainer[iCell][iNode];

            if (keyActiveOfNodes[nodeIndex] == 1)
            {
                keyActiveOfCells[iCell] = 1;
                break;
            }
        }

        if (keyActiveOfCells[iCell] != 1)
        {
            continue;
        }

        for (int iNode = 0; iNode < cellNodeNumber; ++ iNode)
        {
            int nodeIndex = cellNodeIndexContainer[iCell][iNode];

            if (initialInterNodeOrNot[nodeIndex] || innerNodesOrNot[nodeIndex])
            {
                continue;
            }

            if (keyActiveOfNodes[nodeIndex] == 0)
            {
                keyActiveOfNodes[nodeIndex] = 1;
                numberOfNewNodes++;
            }
        }
    }

    return numberOfNewNodes;
}

void OversetConfig::DownInterfaceNodeCharacter(vector< vector< int > > &keyActiveOfGlobalInterNodes, int key)
{
    vector< int > &localToGlobalInterNodesIndex = this->localToGlobalInterNodesIndex;
    vector< int > &localInterNodesList = this->interfaceBoundaryNodesList;
    vector< int > &keyActiveOfNodes = this->GetKeyActiveOfNodes();
    int blockIndex = this->blockIndex;

    for (unsigned iNode = 0; iNode < localInterNodesList.size(); ++ iNode)
    {
        int localNodesIndex = localInterNodesList[iNode];
        int globalNodesIndex = localToGlobalInterNodesIndex[iNode];

        int numGlobal = keyActiveOfGlobalInterNodes[blockIndex][globalNodesIndex];

        if (key == 1)
        {
            innerNodesOrNot[localNodesIndex] = numGlobal;
        }
        else if (key == 0)
        {
            keyActiveOfNodes[localNodesIndex] = numGlobal;
        }
    }
}

void OversetConfig::UpdateInterfaceNodeCharacter(vector< vector< int > > &keyActiveOfGlobalInterNodes, int key)
{
    for (unsigned iNode = 0; iNode < interfaceBoundaryNodesList.size(); ++ iNode)
    {
        int localNodesIndex = interfaceBoundaryNodesList[iNode];
        int globalNodesIndex = localToGlobalInterNodesIndex[iNode];

        int numGlobal = keyActiveOfGlobalInterNodes[blockIndex][globalNodesIndex];

        int numLocalZone;

        if (key == 1)
        {
            numLocalZone = innerNodesOrNot[localNodesIndex];
        }
        else
        {
            numLocalZone = keyActiveOfNodes[localNodesIndex];
        }

        keyActiveOfGlobalInterNodes[blockIndex][globalNodesIndex] = MAX(numGlobal, numLocalZone);
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
    int *faceNodeIndexContainer  = grid->GetFace2Node();
    int *leftCellIndexContainer  = grid->GetLeftCellOfFace();
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

    SetField(keyNodes,  0, numberOfNodes);
    SetField(keyFaces,  0, numberOfFaces);
    SetField(nodeList, -1, numberOfNodes);
    SetField(cellList,  0, numberOfCells);

    int **cellNodeIndexContainer = grid->GetCell2NodeArray();
    int *cellNodeNumberContainer = grid->GetNodeNumberOfEachCell();

    int **cellFaceIndexContainer = grid->GetCell2Face();

    int *cellFaceNumberContainer = grid->GetFaceNumberOfEachCell();

    for (int iCell = 0; iCell < numberOfCells; ++ iCell)
    {
        if (innerCellsOrNot[iCell])
        {
            numberOfCellsInZone ++;
            cellList[iCell] = numberOfCellsInZone;

            for (int iNodes = 0; iNodes < cellNodeNumberContainer[iCell]; ++ iNodes)
            {
                keyNodes[cellNodeIndexContainer[iCell][iNodes]] = 1;
            }

            for (int iFaces = 0; iFaces < cellFaceNumberContainer[iCell]; ++ iFaces)
            {
                keyFaces[cellFaceIndexContainer[iCell][iFaces]] = 1;
            }
        }
    }

    for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
    {
        if (keyNodes[iNode] == 1)
        {
            numberOfNodesInZone ++;
            nodeList[iNode] = numberOfNodesInZone;
        }
    }

    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        if (keyFaces[iFace] == 1)
        {
            numberOfFacesInZone ++;
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

    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
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
    for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
    {
        if (keyNodes[iNode] != 1) continue;
        file << x[iNode] << " ";
        iPrint ++;
        if (iPrint % numberOfWordsInEachLine == 0) file << "\n";
    }
    if (iPrint % numberOfWordsInEachLine != 0) file << "\n";

    //nodes coordinate y
    iPrint = 0;
    for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
    {
        if (keyNodes[iNode] != 1) continue;
        file << y[iNode] << " ";
        iPrint ++;
        if (iPrint % numberOfWordsInEachLine == 0) file << "\n";
    }
    if (iPrint % numberOfWordsInEachLine != 0) file << "\n";

    //nodes coordinate z
    iPrint = 0;
    for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
    {
        if (keyNodes[iNode] != 1) continue;
        file << z[iNode] << " ";
        iPrint ++;
        if (iPrint % numberOfWordsInEachLine == 0) file << "\n";
    }
    if (iPrint % numberOfWordsInEachLine != 0) file << "\n";

    //face nodes number
    iPrint = 0;
    if (PHSPACE::GetDim() == THREE_D)
    {
        for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
        {
            if (keyFaces[iFace] != 1) continue;
            file << faceNodeNumberContainer[iFace] << " ";
            iPrint ++;
            if (iPrint % numberOfWordsInEachLine == 0) file << "\n";
        }
        if (iPrint % numberOfWordsInEachLine != 0) file << "\n";
    }

    //face nodes list
    int count = 0;
    iPrint = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        int faceNodeNumber = faceNodeNumberContainer[iFace];

        if (keyFaces[iFace] == 0)
        {
            count += faceNodeNumber;
            continue;
        }

        for (int iNode = 0; iNode < faceNodeNumberContainer[iFace]; ++ iNode)
        {
            int nodeIndex = nodeList[faceNodeIndexContainer[count++]];
            file << nodeIndex << " ";
            iPrint ++;
            if (iPrint % numberOfWordsInEachLine == 0) file << "\n";
        }
    }
    if (iPrint % numberOfWordsInEachLine != 0) file << "\n";

    //left cell
    iPrint = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        if (keyFaces[iFace] == 0) continue;

        int cellIndex = cellList[leftCellIndexContainer[iFace]];
        file << cellIndex << " ";
        iPrint ++;
        if (iPrint % numberOfWordsInEachLine == 0) file << "\n";
    }
    if (iPrint % numberOfWordsInEachLine != 0) file << "\n";

    //right cell
    iPrint = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
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
        iPrint ++;
        if (iPrint % numberOfWordsInEachLine == 0) file << "\n";
    }
    if (iPrint % numberOfWordsInEachLine != 0) file << "\n";

    PHSPACE::CloseFile(file);
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

    SetField(keyNodes,  0, numberOfNodes);
    SetField(keyFaces,  0, numberOfFaces);
    SetField(nodeList, -1, numberOfNodes);
    SetField(cellList,  0, numberOfCells);

    int **cellNodeIndexContainer = grid->GetCell2NodeArray();
    int *cellNodeNumberContainer = grid->GetNodeNumberOfEachCell();

    int **cellFaceIndexContainer = grid->GetCell2Face();
    int *cellFaceNumberContainer = grid->GetFaceNumberOfEachCell();

    int *keyActiveOfCells = this->GetKeyActiveOfCells();

    for (int iCell = 0; iCell < numberOfCells; ++ iCell)
    {
        if (keyActiveOfCells[iCell] == key)
        {
            numberOfCellsInZone ++;
            cellList[iCell] = numberOfCellsInZone;

            for (int iNodes = 0; iNodes < cellNodeNumberContainer[iCell]; ++ iNodes)
            {
                keyNodes[cellNodeIndexContainer[iCell][iNodes]] = 1;
            }

            for (int iFaces = 0; iFaces < cellFaceNumberContainer[iCell]; ++ iFaces)
            {
                keyFaces[cellFaceIndexContainer[iCell][iFaces]] = 1;
            }
        }
    }

    for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
    {
        if (keyNodes[iNode] == 1)
        {
            numberOfNodesInZone ++;
            nodeList[iNode] = numberOfNodesInZone;
        }
    }

    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        if (keyFaces[iFace] == 1)
        {
            numberOfFacesInZone ++;
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

    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
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
    for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
    {
        if (keyNodes[iNode] != 1) continue;
        file << x[iNode] << " ";
        iPrint ++;
        if (iPrint % numberOfWordsInEachLine == 0) file << "\n";
    }
    if (iPrint % numberOfWordsInEachLine != 0) file << "\n";

    //nodes coordinate y
    iPrint = 0;
    for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
    {
        if (keyNodes[iNode] != 1) continue;
        file << y[iNode] << " ";
        iPrint ++;
        if (iPrint % numberOfWordsInEachLine == 0) file << "\n";
    }
    if (iPrint % numberOfWordsInEachLine != 0) file << "\n";

    //nodes coordinate z
    iPrint = 0;
    for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
    {
        if (keyNodes[iNode] != 1) continue;
        file << z[iNode] << " ";
        iPrint ++;
        if (iPrint % numberOfWordsInEachLine == 0) file << "\n";
    }
    if (iPrint % numberOfWordsInEachLine != 0) file << "\n";

    //face nodes number
    iPrint = 0;
    if (PHSPACE::GetDim() == THREE_D)
    {
        for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
        {
            if (keyFaces[iFace] != 1) continue;
            file << faceNodeNumberContainer[iFace] << " ";
            iPrint ++;
            if (iPrint % numberOfWordsInEachLine == 0) file << "\n";
        }
        if (iPrint % numberOfWordsInEachLine != 0) file << "\n";
    }

    //face nodes list
    int count = 0;
    iPrint = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        int faceNodeNumber = faceNodeNumberContainer[iFace];

        if (keyFaces[iFace] == 0)
        {
            count += faceNodeNumber;
            continue;
        }

        for (int iNode = 0; iNode < faceNodeNumberContainer[iFace]; ++ iNode)
        {
            int nodeIndex = nodeList[faceNodeIndexContainer[count ++]];
            file << nodeIndex << " ";
            iPrint ++;
            if (iPrint % numberOfWordsInEachLine == 0) file << "\n";
        }
    }
    if (iPrint % numberOfWordsInEachLine != 0) file << "\n";

    //left cell
    iPrint = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
    {
        if (keyFaces[iFace] == 0) continue;

        int cellIndex = cellList[leftCellIndexContainer[iFace]];
        file << cellIndex << " ";
        iPrint ++;
        if (iPrint % numberOfWordsInEachLine == 0) file << "\n";
    }
    if (iPrint % numberOfWordsInEachLine != 0) file << "\n";

    //right cell
    iPrint = 0;
    for (int iFace = 0; iFace < numberOfFaces; ++ iFace)
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
        iPrint ++;
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

    for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
    {
        if (keyActiveOfNodes[iNode] == key)
        {
            numberOfSpecifyNodes ++;
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

        for (int iSpecifyNodes = 0; iSpecifyNodes < numberOfSpecifyNodes; ++ iSpecifyNodes)
        {
            int activeNodeIndex = nodeList[iSpecifyNodes];

            file << x[activeNodeIndex] << " ";
            file << y[activeNodeIndex] << " ";
            file << z[activeNodeIndex] << " ";
            file << "\n";
        }

        for (int iSpecifyNodes = 0; iSpecifyNodes < numberOfSpecifyNodes; ++ iSpecifyNodes)
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

    for (int iCell = 0; iCell < numberOfCells; ++ iCell)
    {
        file << keyBufferOfCells[iCell] << endl;
    }

    for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
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
    for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
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

    for (int iNode = 0; iNode < numberOfNodes; ++ iNode)
    {
        file << nodeWalldistDistance[iNode] << endl;
    }

    PHSPACE::CloseFile(file);
}

int * OversetConfig::GetKeyActiveOfCells()
{
    int *keyActiveOfCells = this->grid->GetBlankIndex();
    return keyActiveOfCells;
}

RDouble GetToleranceForOversetSearch()
{
    return GlobalDataBase::GetDoubleParaFromDB("toleranceForOversetSearch");
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

template< typename T >
void PH_AllreduceSepModeForVector(vector< T > &sendingBuffer, vector< T > &receivingBuffer, int numberOfElements, MPI_Datatype mpiDataType, PH_Op op)
{
    T *send = &sendingBuffer[0];
    T *recv = &receivingBuffer[0];

    PH_AllreduceSepMode(send, recv, numberOfElements, mpiDataType, op);
}

}
