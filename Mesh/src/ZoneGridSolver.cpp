#include "ZoneGridSolver.h"
#include "GridAdaptationCGNS.h"
#include "TK_Log.h"
#include "Region.h"
#include "Geo_Element.h"
using namespace std;

namespace PHSPACE
{
ZoneGridSolver::ZoneGridSolver(PHVector1D < CgnsZone * > &cgnsZoneContainer)
{
    minimumLength           = 0.0;
    maximumLength           = 0.0;
    zoneIndex               = 0;
    hasBeenAdapted          = false;
    interfaceInformation    = 0;
    this->cgnsZoneContainer = cgnsZoneContainer;
    pointFactory            = new Geo_PointFactory();
    elementFaceSolver       = new ElementFaceSolver(pointFactory);
    elementPropertySolver   = new ElementPropertySolver();

    elementFaceSolver->SetElementPropertySolver(elementPropertySolver);
    elementPropertySolver->SetElementFaceSolver(elementFaceSolver);

    region = 0;
}

ZoneGridSolver::~ZoneGridSolver()
{
    delete pointFactory;
    delete elementFaceSolver;
    delete elementPropertySolver;
}

void ZoneGridSolver::GenerateComputationalGrid(Grid *gridIn, int gridType)
{
    UnstructGrid * grid = UnstructGridCast(gridIn);

    int parallelAdapt = 1;
    if (parallelAdapt)
    {
        int localZoneIndex  = this->GetZoneIndex();
        int globalZoneIndex = region->GetProcessLocalZoneIndexToGlobalZoneIndex(localZoneIndex);
        grid->SetZoneID(globalZoneIndex);
    }

    ElementFaceSolver *elementFaceSolver = this->GetElementFaceSolver();

    ComputationalFaceSolver *computationalFaceSolver = elementFaceSolver->GetComputationalFaceSolver();

    Geo_PointFactory *pointFactory = GetPointFactory();

    int numberOfComputationalElements = computationalFaceSolver->GetNumberOfElements();
    int numberOfNodes = static_cast<int>(pointFactory->GetNumberOfComputationalNodes());

    PrintToWindow("  nTCell of Computational: ", numberOfComputationalElements, "\n");
    grid->SetNTotalCell(numberOfComputationalElements);
    grid->SetNTotalNode(numberOfNodes);

    RDouble *x = new RDouble [numberOfNodes];
    RDouble *y = new RDouble [numberOfNodes];
    RDouble *z = new RDouble [numberOfNodes];

    grid->SetX(x);
    grid->SetY(y);
    grid->SetZ(z);

    PHVectorInt1D & computationalToGlobalNodeIndexMapping = pointFactory->GetComputationalToGlobalNodeIndexMapping();
    uint_t numberOfNodesForComputation = computationalToGlobalNodeIndexMapping.size();
    for (int iNode = 0; iNode < numberOfNodesForComputation; ++ iNode)
    {
        int nodeIndex = computationalToGlobalNodeIndexMapping[ iNode ];

        x[ iNode ] = pointFactory->GetPoint(nodeIndex).X();
        y[ iNode ] = pointFactory->GetPoint(nodeIndex).Y();
        z[ iNode ] = pointFactory->GetPoint(nodeIndex).Z();
    }

    elementFaceSolver->GenerateComputationalUnstructuredGrid(grid);

    //! Generate cell to node from elementPropertySolver.
    ElementPropertySolver *elementPropertySolver = this->GetElementPropertySolver();
    elementPropertySolver->GenerateCellToNode(grid);
}

void ZoneGridSolver::ComputeMinimumAndMaximumElementEdgeLength()
{
}

void ZoneGridSolver::InitializeAllCgnsZoneElements()
{
}

void ZoneGridSolver::ScanElementBoundaryFace()
{
}

void ZoneGridSolver::PrepareUnstructuredComputationalGrid()
{
    this->ComputeMinimumAndMaximumElementEdgeLength();
    this->InitializeAllCgnsZoneElements();
    
    ElementPropertySolver * elementPropertySolver = this->GetElementPropertySolver();
    elementPropertySolver->InitializeElementIndexAndState();

    //解析所有单元面
    elementPropertySolver->ScanAllElementsToReconstructFaceDataStructure(this->hasBeenAdapted); //Green

    //解析单元边界面
    this->ScanElementBoundaryFace();

    //继续解析？有空看一下
    elementPropertySolver->ScanAllElementsToReconstructFaceDataStructure(this->hasBeenAdapted);

    GenerateInformationForComputationalElement();
}

void ZoneGridSolver::GenerateInformationForComputationalElement()
{
    ElementPropertySolver * elementPropertySolver = this->GetElementPropertySolver();
    int numberOfTotalElements =  elementPropertySolver->GetNumberOfTotalElements(); 

    ElementFaceSolver     * elementFaceSolver = this->GetElementFaceSolver();

    ComputationalFaceSolver * computationalFaceSolver = elementFaceSolver->GetComputationalFaceSolver();

    computationalFaceSolver->SetNumberOfElements(numberOfTotalElements);

    int numberOfTotalFaces = elementFaceSolver->GetNumberOfTotalFaces();

    computationalFaceSolver->SetNumberOfFaces(numberOfTotalFaces);

    PHVectorInt1D & leftElementIndexContainer  = elementFaceSolver->GetLeftElementIndexContainer();
    PHVectorInt1D & rightElementIndexContainer = elementFaceSolver->GetRightElementIndexContainer();

    PHVectorInt1D & leftElementIndexContainerForComputation  = computationalFaceSolver->GetLeftElementIndexContainer();   
    PHVectorInt1D & rightElementIndexContainerForComputation = computationalFaceSolver->GetRightElementIndexContainer();

    PHVectorInt1D & faceBoundaryTypeContainerForComputation = computationalFaceSolver->GetFaceBoundaryTypeContainer();
    PHVectorInt1D & faceBoundaryTypeContainer               = elementFaceSolver->GetFaceBoundaryTypeContainer();
    
    PHVectorInt2D & faceNodeIndexesContainer              = elementFaceSolver->GetFaceNodeIndexesContainer();
    PHVectorInt2D & computationalFaceNodeIndexesContainer = computationalFaceSolver->GetFaceNodeIndexesContainer();

    int iBoundaryFace = 0;

    for (int iFace = 0; iFace < numberOfTotalFaces; ++ iFace)
    {
        int rightElementIndex = rightElementIndexContainer[ iFace ];

        //没有转换关系时两者等价
        leftElementIndexContainerForComputation .push_back(leftElementIndexContainer [ iFace ]);
        rightElementIndexContainerForComputation.push_back(rightElementIndexContainer[ iFace ]);
        computationalFaceNodeIndexesContainer.push_back(faceNodeIndexesContainer[ iFace ]);

        if (rightElementIndex == PHSPACE::INVALID_INDEX) //边界
        {
            faceBoundaryTypeContainerForComputation.push_back(faceBoundaryTypeContainer[ iFace ]);
            ++ iBoundaryFace;
        }
    }

    computationalFaceSolver->SetNumberOfBoundaryFaces(iBoundaryFace);

    Geo_PointFactory * pointFactory = this->GetPointFactory();

    PHVectorInt1D & computationalToGlobalNodeIndexMapping = pointFactory->GetComputationalToGlobalNodeIndexMapping(); 
    uint_t numberOfTotalPoints = pointFactory->GetNumberOfPoints();
    computationalToGlobalNodeIndexMapping.resize(numberOfTotalPoints);
    for (int iNode = 0; iNode < numberOfTotalPoints; ++ iNode)
    {
        computationalToGlobalNodeIndexMapping[ iNode ] = iNode;
    }
}

void ZoneGridSolver::CompressStatusIntoDataContainer(DataContainer *& dataContainer, int neighborZoneIndex)
{

}

void ZoneGridSolver::CompressFaceRBFWeightIntoDataContainer(DataContainer *& dataContainer, int neighborZoneIndex)
{

}

void ZoneGridSolver::AnasysStatusFromDataContainer(DataContainer *& dataContainer, int neighborZoneIndex)
{

}

void ZoneGridSolver::AnasysFaceRBFWeightFromDataContainer(DataContainer *& dataContainer, int neighborZoneIndex)
{

}

void ZoneGridSolver::GetSourceCellIndex(int iFace, int iPosition, int & sourceIndex)
{
    InterfaceInfo * interfaceInformation = this->GetInterfaceInformation();
    int * localInterfaceToBoundaryFaceIndexesMapping = interfaceInformation->GetInterFace2BoundaryFace();

    PHVectorInt1D & leftCellIndexContainer = this->GetElementFaceSolver()->GetLeftElementIndexContainer();
    int boundaryFaceIndex = localInterfaceToBoundaryFaceIndexesMapping[ iFace ];
    sourceIndex = leftCellIndexContainer[ boundaryFaceIndex ];
}

void ZoneGridSolver::GetTargetIndex(int iFace, int iPosition, int & targetIndex)
{
    InterfaceInfo * interfaceInformation = this->GetInterfaceInformation();
    int * localInterfaceToBoundaryFaceIndexesMapping = interfaceInformation->GetInterFace2BoundaryFace();
    
    PHVectorInt1D & rightCellIndexContainer = this->GetElementFaceSolver()->GetRightElementIndexContainer();
    int boundaryFaceIndex = localInterfaceToBoundaryFaceIndexesMapping[ iFace ];
    targetIndex = rightCellIndexContainer[ boundaryFaceIndex ];
    targetIndex = - targetIndex;
}

void ZoneGridSolver::GetSourceFaceIndex(int iFace, int & sourceIndex)
{
    InterfaceInfo * interfaceInformation = this->GetInterfaceInformation();
    int * localInterfaceToBoundaryFaceIndexesMapping = interfaceInformation->GetInterFace2BoundaryFace();

    int boundaryFaceIndex = localInterfaceToBoundaryFaceIndexesMapping[ iFace ];
    sourceIndex = boundaryFaceIndex;
}

void ZoneGridSolver::GetTargetFaceIndex(int iFace, int & targetIndex)
{
    InterfaceInfo * interfaceInformation = this->GetInterfaceInformation();
    int * localInterfaceToBoundaryFaceIndexesMapping = interfaceInformation->GetInterFace2BoundaryFace();

    int boundaryFaceIndex = localInterfaceToBoundaryFaceIndexesMapping[ iFace ];
    targetIndex = boundaryFaceIndex;
}

void ZoneGridSolver::ResetInterfaceElementProperty()
{
    PHVectorInt1D & faceBoundaryTypeContainer = this->GetElementFaceSolver()->GetFaceBoundaryTypeContainer();

    int numberOfBoundaryFaces = this->GetNumberOfBoundaryFaces();

    PHVectorInt1D & elementModifiedStatusContainer      = this->GetElementPropertySolver()->GetElementModifiedStatusContainer();
    PHVectorInt1D & elementComputationalStatusContainer = this->GetElementPropertySolver()->GetElementComputationalStatusContainer();

    PHVectorInt1D & leftElementIndexContainer = this->GetElementFaceSolver()->GetLeftElementIndexContainer();
    PHVectorInt1D & rightCellIndexContainer   = this->GetElementFaceSolver()->GetRightElementIndexContainer();
    for (int iFace = 0; iFace < numberOfBoundaryFaces; ++ iFace)
    {
        int boundaryType = faceBoundaryTypeContainer[ iFace ];
        if (boundaryType == PHENGLEI::INTERFACE)
        {
            int leftCellIndex  = leftElementIndexContainer[ iFace ];
            int rightCellIndex = - rightCellIndexContainer[ iFace ];

            int modifiedStatusStatus = elementModifiedStatusContainer[ rightCellIndex ];
            if (modifiedStatusStatus == PHSPACE::SPLIT || modifiedStatusStatus == PHSPACE::MERGE_ELEMENT)
            {
                elementModifiedStatusContainer[ leftCellIndex ] = modifiedStatusStatus;
            }

            if (modifiedStatusStatus == PHSPACE::SPLIT)
            {
                elementComputationalStatusContainer[ leftCellIndex ] = PHSPACE::HIDDEN;
            }
            else if (modifiedStatusStatus == PHSPACE::MERGE_ELEMENT)
            {
                elementComputationalStatusContainer[ leftCellIndex ] = PHSPACE::ON;
            }
        }
    }
}

void ZoneGridSolver::CompressRefineType(DataContainer *& dataContainer, int neighborZoneIndex)
{
    InterfaceInfo * interfaceInformation = this->GetInterfaceInformation();
    if (!interfaceInformation) return;

    ElementPropertySolver * elementPropertySolver = this->GetElementPropertySolver();
    PHVectorInt1D & refineType = * elementPropertySolver->GetRefineType();

    int iNeighborZone                            = interfaceInformation->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = interfaceInformation->GetNIFaceOfNeighbor(iNeighborZone);
    int * interfaceIndexContainerForSend         = interfaceInformation->GetFaceIndexForSend(iNeighborZone);
    dataContainer->moveto_begin();

    int * dataTemp = new int [interfaceNumberBetweenTwoNeighboringZone];
    for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
    {
        int sourceCell;
        int iFace = interfaceIndexContainerForSend[ iLocalFace ];
        this->GetSourceCellIndex(iFace, 1, sourceCell);

        dataTemp[iLocalFace] = refineType[sourceCell];
    }
    PHWrite(dataContainer, dataTemp, interfaceNumberBetweenTwoNeighboringZone);
    delete [] dataTemp;

 /*
    for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
    {
        int sourceCell;
        int iFace = interfaceIndexContainerForSend[ iLocalFace ];
        this->GetSourceCellIndex(iFace, 1, sourceCell);

        PHWrite(dataContainer, refineType[sourceCell]);
    }
    //*/
}

void ZoneGridSolver::DecompressRefineType(DataContainer *& dataContainer, int neighborZoneIndex)
{
    InterfaceInfo * interfaceInformation = this->GetInterfaceInformation();
    if (!interfaceInformation) return;

    int iNeighborZone                            = interfaceInformation->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = interfaceInformation->GetNIFaceOfNeighbor(iNeighborZone);
    int * interfaceIndexContainerForReceive      = interfaceInformation->GetFaceIndexForRecv(iNeighborZone);

    ElementPropertySolver * elementPropertySolver = this->GetElementPropertySolver();
    PHVectorInt1D & refineType = * elementPropertySolver->GetRefineType();

    dataContainer->moveto_begin();
    int * dataTemp = new int [interfaceNumberBetweenTwoNeighboringZone];
    PHRead(dataContainer, dataTemp, interfaceNumberBetweenTwoNeighboringZone);

    for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
    {
        int iFace = interfaceIndexContainerForReceive[ iLocalFace ];
        int targetCell;
        this->GetTargetIndex(iFace, 1, targetCell);

        refineType[targetCell] = dataTemp[iLocalFace];
    }

    delete [] dataTemp;
}

void ZoneGridSolver::CompressCellProjectType(DataContainer *& dataContainer, int neighborZoneIndex)
{
    InterfaceInfo * interfaceInformation = this->GetInterfaceInformation();
    if (!interfaceInformation) return;

    ElementPropertySolver * elementPropertySolver = this->GetElementPropertySolver();
    PHVectorInt1D & projectType = * elementPropertySolver->GetProjectType();
    PHVectorInt1D & wallBCID    = * elementPropertySolver->GetWallBCID();

    int iNeighborZone                            = interfaceInformation->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = interfaceInformation->GetNIFaceOfNeighbor(iNeighborZone);
    int * interfaceIndexContainerForSend         = interfaceInformation->GetFaceIndexForSend(iNeighborZone);
    dataContainer->moveto_begin();

    int * dataTemp = new int [interfaceNumberBetweenTwoNeighboringZone * 2];
    for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
    {
        int sourceCell;
        int iFace = interfaceIndexContainerForSend[ iLocalFace ];
        this->GetSourceCellIndex(iFace, 1, sourceCell);

        dataTemp[iLocalFace] = projectType[sourceCell];
    }
    for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
    {
        int sourceCell;
        int iFace = interfaceIndexContainerForSend[ iLocalFace ];
        this->GetSourceCellIndex(iFace, 1, sourceCell);

        dataTemp[interfaceNumberBetweenTwoNeighboringZone + iLocalFace] = wallBCID[sourceCell];
    }
    PHWrite(dataContainer, dataTemp, interfaceNumberBetweenTwoNeighboringZone * 2);
    delete [] dataTemp;
}

void ZoneGridSolver::DecompressCellProjectType(DataContainer *& dataContainer, int neighborZoneIndex)
{
    InterfaceInfo * interfaceInformation = this->GetInterfaceInformation();
    if (!interfaceInformation) return;

    int iNeighborZone                            = interfaceInformation->FindIthNeighbor(neighborZoneIndex);
    int interfaceNumberBetweenTwoNeighboringZone = interfaceInformation->GetNIFaceOfNeighbor(iNeighborZone);
    int * interfaceIndexContainerForReceive      = interfaceInformation->GetFaceIndexForRecv(iNeighborZone);

    ElementPropertySolver * elementPropertySolver = this->GetElementPropertySolver();
    PHVectorInt1D & projectType = * elementPropertySolver->GetProjectType();
    PHVectorInt1D & wallBCID    = * elementPropertySolver->GetWallBCID();

    dataContainer->moveto_begin();
    int * dataTemp = new int [interfaceNumberBetweenTwoNeighboringZone * 2];
    PHRead(dataContainer, dataTemp, interfaceNumberBetweenTwoNeighboringZone * 2);

    for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
    {
        int iFace = interfaceIndexContainerForReceive[ iLocalFace ];
        int targetCell;
        this->GetTargetIndex(iFace, 1, targetCell);

        projectType[targetCell] = dataTemp[iLocalFace];
    }
    for (int iLocalFace = 0; iLocalFace < interfaceNumberBetweenTwoNeighboringZone; ++ iLocalFace)
    {
        int iFace = interfaceIndexContainerForReceive[ iLocalFace ];
        int targetCell;
        this->GetTargetIndex(iFace, 1, targetCell);

        wallBCID[targetCell] = dataTemp[interfaceNumberBetweenTwoNeighboringZone + iLocalFace];
    }

    delete [] dataTemp;
}

void SwapInterfaceData(Region * region, ZoneGridSolver ** zoneGridSolverList, CompressData compressData, DeCompressData decompressData)
{
    int currentProcessor = PHMPI::GetCurrentProcessorID();
    int nZones = PHMPI::GetNumberofGlobalZones();
    vector< PH_Request > requestContainer;
    vector< vector< DataContainer * > > receivedDataBuffer;
    vector< vector< DataContainer * > > sendDataBuffer;
    receivedDataBuffer.resize(nZones);
    sendDataBuffer.resize(nZones);

    //! Very Important!!!!!
    //! Warning: this function should not be used before the interfaces of newly refined grid
    //!          have been reconstructed! Because the connection is old.

    // Step 0: Compressing data firstly.
    WriteLogFile("        Compressing data ...");
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        int sendProcessor = PHMPI::GetZoneProcessorID(iZone);

        ZoneNeighbor * globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
        std::size_t numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();
        if(currentProcessor == sendProcessor)
        {
            sendDataBuffer[iZone].resize(numberOfNeighbor);
        }

        for (std::size_t iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
        {
            int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);

            int receiveProcessor = PHMPI::GetZoneProcessorID(neighborZone);

            if(currentProcessor != sendProcessor && currentProcessor != receiveProcessor)
            {
                continue;
            }

            // Allocate and compress the buffers for sending.
            DataContainer * sendBuffer = 0;
            if(currentProcessor == sendProcessor)
            {
                sendBuffer = new DataContainer;
                sendDataBuffer[ iZone ][iNeighbor] = sendBuffer;

                // Compress the send information into the actkey.
                int localZoneIndex = region->GetProcessGlobalZoneIndexToLocalZoneIndex(iZone);
                ZoneGridSolver * zoneGridSolver = zoneGridSolverList[localZoneIndex];
                compressData(zoneGridSolver, sendBuffer, neighborZone);
            }
        }
    }


    // Step 1: Communication.
    WriteLogFile("        Communication ...");
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        ZoneNeighbor * globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
        std::size_t numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();

        // Communicating.
        for (std::size_t iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
        {
            int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);

            int sendProcessor    = PHMPI::GetZoneProcessorID(iZone);
            int receiveProcessor = PHMPI::GetZoneProcessorID(neighborZone);

            if(currentProcessor != sendProcessor && currentProcessor != receiveProcessor)
            {
                continue;
            }

            // Allocate the buffers for sending and receiving.
            DataContainer * receivedBuffer = 0;
            DataContainer * sendBuffer = 0;
            if(currentProcessor == receiveProcessor)
            {
                receivedBuffer = new DataContainer;
                receivedDataBuffer[ iZone ].push_back(receivedBuffer);
            }
            if(currentProcessor == sendProcessor)
            {
                sendBuffer = sendDataBuffer[iZone][iNeighbor];
            }

            int tag = iZone;

            //if (sendProcessor == receiveProcessor) continue;

            // Communication.
            if (currentProcessor == sendProcessor)
            {
                if (sendProcessor != receiveProcessor)
                {
                    streamsize nlen = sendBuffer->size();

                    // Send the data to neighbors.
                    send(sendBuffer, receiveProcessor, nlen, requestContainer, tag);
                }
            }
            if (currentProcessor == receiveProcessor)
            {
                // Get the neighbor order: which order that 'iZone' on neighbors of neighborZone.
                int neighborOrer = -1;
                ZoneNeighbor * globalNeighborZonesTemp = zoneConnectivity->GetZoneNeighbor(neighborZone);
                std::size_t numberOfNeighborTemp = globalNeighborZonesTemp->GetNumberOfNeighbor();

                for (std::size_t iNeighbor = 0; iNeighbor < numberOfNeighborTemp; ++ iNeighbor)
                {
                    uint_t neighborZoneTemp = globalNeighborZonesTemp->GetZoneIndexOfNeighbor(iNeighbor);
                    if(neighborZoneTemp == iZone)
                    {
                        neighborOrer = static_cast<int>(iNeighbor);
                        break;
                    }
                }

                ASSERT(neighborOrer != -1);

                if (sendProcessor != receiveProcessor)
                {
                    // the data length of the received data is same to the length that send to the neighbor.
                    CharVecSizeType nlen = sendDataBuffer[neighborZone][neighborOrer]->size();

                    // Receive data from neighbors.
                    receive(receivedBuffer, sendProcessor, nlen, requestContainer, tag);
                }
                else
                {
                    uint_t lastID = receivedDataBuffer[iZone].size();
                    SWAP(receivedDataBuffer[ iZone ][lastID - 1], sendDataBuffer[neighborZone][neighborOrer]);
                }
            }
        }
    }

    // Step2: MPI waiting for Non-blocking communications.
    if(PHMPI::GetNumberOfProcessor() > 1)
    {
        PH_Wait(static_cast<int>(requestContainer.size()), &(requestContainer[0]));
    }

    // Step3: Translate the data container.
    WriteLogFile("        Translating data ...");
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        ZoneNeighbor * globalNeighborZones = zoneConnectivity->GetZoneNeighbor(iZone);
        std::size_t numberOfNeighbor = globalNeighborZones->GetNumberOfNeighbor();

        int count = 0;
        for (std::size_t iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
        {
            int neighborZone = globalNeighborZones->GetZoneIndexOfNeighbor(iNeighbor);

            int receiveProcessor = PHMPI::GetZoneProcessorID(neighborZone);

            if (currentProcessor == receiveProcessor)
            {
                DataContainer * receiveData = receivedDataBuffer[ iZone ][ count ];

                //! Because get into here when currentProcessor==receiveProcessor,
                //! so, use zoneGridSolver of neighborZone, this zone is on the current processor!
                int localZoneIndex = region->GetProcessGlobalZoneIndexToLocalZoneIndex(neighborZone);
                ZoneGridSolver * zoneGridSolver = zoneGridSolverList[localZoneIndex];

                // Decompress the interface data from datacontainer.
                decompressData(zoneGridSolver, receiveData, iZone);

                ++ count;
            }
        }
    }

    WriteLogFile("        Free buffer data ...");
    // Step4: Free the buffers.
    for (std::size_t iDim = 0; iDim < receivedDataBuffer.size(); ++ iDim)
    {
        for (std::size_t jDim = 0; jDim < receivedDataBuffer[ iDim ].size(); ++ jDim)
        {
            delete receivedDataBuffer[ iDim ][ jDim ];
        }
    }
    receivedDataBuffer.clear();
    for (std::size_t iDim = 0; iDim < sendDataBuffer.size(); ++ iDim)
    {
        for (std::size_t jDim = 0; jDim < sendDataBuffer[ iDim ].size(); ++ jDim)
        {
            delete sendDataBuffer[ iDim ][ jDim ];
        }
    }
    sendDataBuffer.clear();
    requestContainer.clear();
}

}