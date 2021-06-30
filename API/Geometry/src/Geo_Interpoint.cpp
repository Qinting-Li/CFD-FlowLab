#include "Geo_Interpoint.h"
#include "PHMatrix.h"
#include "PHHeader.h"
#include "TK_Exit.h"
#include "Geo_SimpleGrid.h"
#include "Geo_Interface.h"
#include "Constants.h"

namespace PHSPACE
{
InterpointPatch::InterpointPatch()
{
    pointIndexForReceive = 0;
    pointIndexForSend = 0;
    dataStorageSend = new Data_ParamFieldSuite();
    dataStorageReceive = new Data_ParamFieldSuite();
    zoneID = 0;
    zoneIndexOfNeighbor = 0;
    pointNumberOfPatch = 0;
}

InterpointPatch::InterpointPatch(const InterpointPatch &interPatchIn)
{
    zoneID = interPatchIn.GetZoneID();
    zoneIndexOfNeighbor = interPatchIn.GetZoneIndexOfNeighbor();
    pointNumberOfPatch = interPatchIn.GetNumberOfPoint();

    if (interPatchIn.GetPointIndexForReceive())
    {
        pointIndexForReceive = new int[pointNumberOfPatch];
        PHSPACE::CopyField(pointIndexForReceive, interPatchIn.GetPointIndexForReceive(), pointNumberOfPatch);
    }

    if (interPatchIn.GetPointIndexForSend())
    {
        pointIndexForSend = new int[pointNumberOfPatch];
        PHSPACE::CopyField(pointIndexForSend, interPatchIn.GetPointIndexForSend(), pointNumberOfPatch);
    }

    dataStorageSend = new Data_ParamFieldSuite();
    dataStorageReceive = new Data_ParamFieldSuite();
}

InterpointPatch::~InterpointPatch()
{
    delete dataStorageSend;
    delete dataStorageReceive;
    if (pointIndexForReceive)
    {
        delete []pointIndexForReceive;
    }
    if (pointIndexForSend)
    {
        delete []pointIndexForSend;
    }
}

void InterpointPatch::AllocateData()
{
    int pointNumberOfPatch = this->GetNumberOfPoint();

    int *receivePoint = new int[pointNumberOfPatch];
    int *sendPoint = new int[pointNumberOfPatch];

    this->SetPointIndexForReceive(receivePoint);
    this->SetPointIndexForSend(sendPoint);
}

InterpointInformation::InterpointInformation()
{
    numberOfInterpoints = 0;
    numberOfNeighbor = 0;
    cellNumberOfInterPoint = 0;
    totalZonesOfInterPoint = 0;
    labelOfInterPoint = 0;
    parent = 0;
    interPoint2ZoneID = 0;
    interPoint2InterPointID = 0;
    interPoint2GlobalPoint = 0;

    dataForSend = new Data_ParamFieldSuite[GetNumberOfGhostCellLayers()];
    dataForReceive = new Data_ParamFieldSuite[GetNumberOfGhostCellLayers()];

    interPointPatch = 0;
    isNeighborInterpointFound = 0;
}

InterpointInformation::InterpointInformation(const InterpointInformation &interpointInformationIn)
{
    numberOfInterpoints = interpointInformationIn.GetNumberOfInterpoints();
    numberOfNeighbor = interpointInformationIn.GetNumberOfNeighbor();
    if (interpointInformationIn.GetInterPoint2ZoneID())
    {
        interPoint2ZoneID = new int[numberOfInterpoints];
        PHSPACE::CopyField(interPoint2ZoneID, interpointInformationIn.GetInterPoint2ZoneID(), numberOfInterpoints);
    }

    if (interpointInformationIn.GetInterPoint2InterPointID())
    {
        interPoint2InterPointID = new int[numberOfInterpoints];
        PHSPACE::CopyField(interPoint2InterPointID, interpointInformationIn.GetInterPoint2InterPointID(), numberOfInterpoints);
    }

    if (interpointInformationIn.GetInterPoint2GlobalPoint())
    {
        interPoint2GlobalPoint = new int[numberOfInterpoints];
        PHSPACE::CopyField(interPoint2GlobalPoint, interpointInformationIn.GetInterPoint2GlobalPoint(), numberOfInterpoints);
    }

    if (interpointInformationIn.GetCellNumberOfInterPoint())
    {
        cellNumberOfInterPoint = new int[numberOfInterpoints];
        PHSPACE::CopyField(cellNumberOfInterPoint, interpointInformationIn.GetCellNumberOfInterPoint(), numberOfInterpoints);
    }

    if (interpointInformationIn.GetTotalZonesOfInterPoint())
    {
        totalZonesOfInterPoint = new int[numberOfInterpoints];
        PHSPACE::CopyField(totalZonesOfInterPoint, interpointInformationIn.GetTotalZonesOfInterPoint(), numberOfInterpoints);
    }

    if (interpointInformationIn.GetLabelOfInterPoint())
    {
        labelOfInterPoint = new int[numberOfInterpoints];
        PHSPACE::CopyField(labelOfInterPoint, interpointInformationIn.GetLabelOfInterPoint(), numberOfInterpoints);
    }
    
    if (interpointInformationIn.GetInterpointPatch())
    {
        int numberOfNeighbor = this->GetNumberOfNeighbor();
        interPointPatch = new InterpointPatch *[numberOfNeighbor];
        for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
        {
            interPointPatch[iNeighbor] = new InterpointPatch(*interpointInformationIn.GetInterpointPatch(iNeighbor));
        }
    }

    dataForSend = new Data_ParamFieldSuite[GetNumberOfGhostCellLayers()];
    dataForReceive = new Data_ParamFieldSuite[GetNumberOfGhostCellLayers()];

    parent = 0;

    isNeighborInterpointFound = 0;
}

InterpointInformation::InterpointInformation(int numberOfInterpoints, Grid *parent)
{
    this->numberOfInterpoints = numberOfInterpoints;
    interPoint2ZoneID = new int[numberOfInterpoints];
    interPoint2InterPointID = new int[numberOfInterpoints];
    interPoint2GlobalPoint = new int[numberOfInterpoints];
    cellNumberOfInterPoint = new int[numberOfInterpoints];
    totalZonesOfInterPoint = new int[numberOfInterpoints];
    labelOfInterPoint = new int[numberOfInterpoints];

    dataForSend = new Data_ParamFieldSuite[GetNumberOfGhostCellLayers()];
    dataForReceive = new Data_ParamFieldSuite[GetNumberOfGhostCellLayers()];
    this->parent = parent;

    interPointPatch = 0;

    isNeighborInterpointFound = 0;
    numberOfNeighbor = 0;
}

InterpointInformation::~InterpointInformation()
{
    delete []interPoint2ZoneID;
    delete []interPoint2InterPointID;
    delete []interPoint2GlobalPoint;
    delete []cellNumberOfInterPoint;
    delete []totalZonesOfInterPoint;
    delete []labelOfInterPoint;
    delete []dataForSend;
    delete []dataForReceive;

    if (interPointPatch)
    {
        for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
        {
            delete interPointPatch[iNeighbor];
        }
        delete []interPointPatch;
    }

    delete []isNeighborInterpointFound;
}

void InterpointInformation::AllocateNeighborInfo()
{
    int numberOfNeighbor = this->GetNumberOfNeighbor();
    interPointPatch = new InterpointPatch *[numberOfNeighbor];
    for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
    {
        interPointPatch[iNeighbor] = new InterpointPatch();
    }
}

int InterpointInformation::ComputeNumberOfNeighbor(int *zoneFlag)
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    int *interPoint2ZoneID = this->GetInterPoint2ZoneID();
    int numberOfInterpoints = this->GetNumberOfInterpoints();

    //! Set zoneflag to zero,assuming that the current zone has no neighbor.
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        zoneFlag[iZone] = 0;
    }

    //! Set zoneflag to zero,assuming that the current zone has no neighbor.
    for (int iPoint = 0; iPoint < numberOfInterpoints; ++ iPoint)
    {
        zoneFlag[interPoint2ZoneID[iPoint]] = 1;
    }

    //! Compute the number of neighbor zones of the current zone,including itself.
    int numberOfNeighbor = 0;
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        numberOfNeighbor += zoneFlag[iZone];
    }

    return numberOfNeighbor;
}

void InterpointInformation::FindNeighbors()
{
    using namespace PHMPI;
    int nZones = GetNumberofGlobalZones();

    int *zoneFlag = new int[nZones];

    int numberOfNeighbor = ComputeNumberOfNeighbor(zoneFlag);

    this->SetNumberOfNeighbor(numberOfNeighbor);
    this->AllocateNeighborInfo();

    int iNeighbor = 0;
    for (int iZone = 0; iZone < nZones; ++ iZone)
    {
        if (zoneFlag[iZone])
        {
            InterpointPatch *interPointPatch = this->GetInterpointPatch(iNeighbor ++);
            interPointPatch->SetZoneIndexOfNeighbor(iZone);
        }
    }

    delete []zoneFlag;
}

void InterpointInformation::ComputePointIndexForReceive(int iNeighbor)
{
    int numberOfInterpoints = this->GetNumberOfInterpoints();

    int *interPoint2ZoneID = this->GetInterPoint2ZoneID();
    int zoneIDOfNeighbor = this->GetZoneIndexOfNeighbor(iNeighbor);

    InterpointPatch *interPointPatch = this->GetInterpointPatch(iNeighbor);
    int *pointIndexForReceive = interPointPatch->GetPointIndexForReceive();
    int count = 0;
    for (int iPoint = 0; iPoint < numberOfInterpoints; ++ iPoint)
    {
        if (interPoint2ZoneID[iPoint] == zoneIDOfNeighbor)
        {
            //! It means that pointIndexForReceive is count locally  based on the interfaces of current zone.
            pointIndexForReceive[count] = iPoint;
            ++ count;
        }
    }
}

int InterpointInformation::ComputeNumberOfInterpointsForNeighbor(int iNeighbor)
{
    int numberOfInterpoints = this->GetNumberOfInterpoints();
    int *interPoint2ZoneID = this->GetInterPoint2ZoneID();
    int zoneIDOfNeighbor = this->GetZoneIndexOfNeighbor(iNeighbor);
    int NumberOfInterpointsForNeighbor = 0;

    for (int iPoint = 0; iPoint < numberOfInterpoints; ++ iPoint)
    {
        if (interPoint2ZoneID[iPoint] == zoneIDOfNeighbor)
        {
            ++ NumberOfInterpointsForNeighbor;
        }
    }

    return NumberOfInterpointsForNeighbor;
}

void InterpointInformation::AllocateSendReceive(int iNeighbor)
{
    int NumberOfInterpointsForNeighbor = this->ComputeNumberOfInterpointsForNeighbor(iNeighbor);

    InterpointPatch *interPointPatch = this->GetInterpointPatch(iNeighbor);
    interPointPatch->SetNumberOfPoint(NumberOfInterpointsForNeighbor);
    interPointPatch->AllocateData();
}

void InterpointInformation::InitNeighbors(int iNeighbor)
{
    this->AllocateSendReceive(iNeighbor);
    this->ComputePointIndexForReceive(iNeighbor);
}

void InterpointInformation::InitNeighbors()
{
    FindNeighbors();

    int numberOfNeighbor = this->GetNumberOfNeighbor();
    for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
    {
        InitNeighbors(iNeighbor);
    }
}

int InterpointInformation::FindIthNeighbor(int zone)
{
    int idx = -1;
    for (int iNeighbor = 0; iNeighbor < numberOfNeighbor; ++ iNeighbor)
    {
        if (zone == this->GetZoneIndexOfNeighbor(iNeighbor))
        {
            idx = iNeighbor;
            break;
        }
    }
    return idx;
}

int *InterpointInformation::ComputePointIndexForSend(int iNeighbor)
{
    int *pointIndexForSend = new int[this->GetNumberOfInterpointsForNeighbor(iNeighbor)];
    int *interPoint2InterPointID = this->GetInterPoint2InterPointID();
    int zoneIDOfNeighbor = this->GetZoneIndexOfNeighbor(iNeighbor);

    int count = 0;
    for (int iPoint = 0; iPoint < numberOfInterpoints; ++ iPoint)
    {
        if (interPoint2ZoneID[iPoint] == zoneIDOfNeighbor)
        {
            pointIndexForSend[count] = interPoint2InterPointID[iPoint];
            ++ count;
        }
    }
    return pointIndexForSend;
}

void InterpointInformation::FillPointIndexForSend(int zone, int *pointIndexForSendIn)
{
    int iNeighbor = this->FindIthNeighbor(zone);
    InterpointPatch *interPointPatch = this->GetInterpointPatch(iNeighbor);
    int *pointIndexForSend = interPointPatch->GetPointIndexForSend();
    int pointNumberOfPatch = interPointPatch->GetNumberOfPoint();
    for (int iPoint = 0; iPoint < pointNumberOfPatch; ++ iPoint)
    {
		//! Here we can make out that pointIndexForSend is corresponds to the neighbor zone.
		//! That is to say,it comes according to the corresponding sequence of corresponding neighbor zone.
        pointIndexForSend[iPoint] = pointIndexForSendIn[iPoint];
    }
}

void InterpointInformation::ReSize(int numberOfNewInterpoints)
{
    if (this->numberOfInterpoints == numberOfNewInterpoints)
    {
        return;
    }

    delete []interPoint2ZoneID;
    delete []interPoint2InterPointID;
    delete []interPoint2GlobalPoint;
    delete []cellNumberOfInterPoint;
    delete []totalZonesOfInterPoint;
    delete []labelOfInterPoint;

    this->numberOfInterpoints = numberOfNewInterpoints;

    interPoint2ZoneID = new int[numberOfInterpoints];
    interPoint2InterPointID = new int[numberOfInterpoints];
    interPoint2GlobalPoint = new int[numberOfInterpoints];
    cellNumberOfInterPoint = new int[numberOfInterpoints];
    totalZonesOfInterPoint = new int[numberOfInterpoints];
    labelOfInterPoint = new int[numberOfInterpoints];
}

InterpointFields::InterpointFields(InterpointInformation *interpointInformationIn)
{
    this->interpointInformation = interpointInformationIn;
}

InterpointFields::~InterpointFields()
{
    //! The send &&receive data are allocated by InterpointFields,
    //! so, they should be free by InterpointFields.
    FreeInterpointVariable();
}

void InterpointFields::RegisterField(const string &name, const int type, const int dimesion, const int solverID)
{
    variableNames.push_back(name);
    variableTypes.push_back(type);
    variableDimensions.push_back(dimesion);
    solverIndex.push_back(solverID);

    for (int IDOfGhostLayer = GetNumberOfGhostCellLayers() - 1; IDOfGhostLayer >= 0; -- IDOfGhostLayer)
    {
        AllocateInterpointVariable(name, type, dimesion, interpointInformation->GetSendDataStorage(IDOfGhostLayer));
        AllocateInterpointVariable(name, type, dimesion, interpointInformation->GetReceiveDataStorage(IDOfGhostLayer));
    }
}

void InterpointFields::AllocateInterpointVariable(const string &name, const int type, const int dimesion, Data_ParamFieldSuite *datastore)
{
    int numberOfInterpoints = interpointInformation->GetNumberOfInterpoints();

    switch (type)
    {
    case PHINT:
        {
            int **qInt = NewPointer2<int>(dimesion, numberOfInterpoints);
            for (int iDimension = 0; iDimension < dimesion; ++ iDimension)
            {
                for (int jDimension = 0; jDimension < numberOfInterpoints; ++ jDimension)
                {
                    qInt[iDimension][jDimension] = 0;
                }
            }
            datastore->UpdateDataPtr(name, qInt);
            break;
        }
    case PHFLOAT:
        {
            float **qFloat = NewPointer2<float>(dimesion, numberOfInterpoints);
            for (int iDimension = 0; iDimension < dimesion; ++ iDimension)
            {
                for (int jDimension = 0; jDimension < numberOfInterpoints; ++ jDimension)
                {
                    qFloat[iDimension][jDimension] = 0.0;
                }
            }
            datastore->UpdateDataPtr(name, qFloat);
            break;
        }
    case PHDOUBLE:
        {
            RDouble **qDouble = NewPointer2<RDouble>(dimesion, numberOfInterpoints);
            for (int iDimension = 0; iDimension < dimesion; ++ iDimension)
            {
                for (int jDimension = 0; jDimension < numberOfInterpoints; ++ jDimension)
                {
                    qDouble[iDimension][jDimension] = 0.0;
                }
            }
            datastore->UpdateDataPtr(name, qDouble);
            break;
        }
    case PHBOOL:
        {
            bool **qBool = NewPointer2<bool>(dimesion, numberOfInterpoints);
            for (int iDimension = 0; iDimension < dimesion; ++ iDimension)
            {
                for (int jDimension = 0; jDimension < numberOfInterpoints; ++ jDimension)
                {
                    qBool[iDimension][jDimension] = 0;
                }
            }
            datastore->UpdateDataPtr(name, qBool);
            break;
        }
    default:
        {
            TK_Exit::UnexpectedVarValue("type", type);
            break;
        }
    }
}

void InterpointFields::FreeInterpointVariable()
{
    for (int IDOfGhostLayer = GetNumberOfGhostCellLayers() - 1; IDOfGhostLayer >= 0; -- IDOfGhostLayer)
    {
        FreeInterpointVariable(interpointInformation->GetSendDataStorage(IDOfGhostLayer));
        FreeInterpointVariable(interpointInformation->GetReceiveDataStorage(IDOfGhostLayer));
    }
}

void InterpointFields::FreeInterpointVariable(Data_ParamFieldSuite *dataStore)
{
    for (int iVariable = 0; iVariable < Size(); ++ iVariable)
    {
        string variableName = variableNames[iVariable];
        void **q = reinterpret_cast<void **> (dataStore->GetDataPtr(variableName));
        DelPointer2(q);
    }
}

template<typename T>
T **InterpointFields::GetSendVariable(const string &name, const int IDOfGhostLayer)
{
     return reinterpret_cast<T **> (interpointInformation->GetSendDataStorage(IDOfGhostLayer)->GetDataPtr(name));
}

template<typename T>
T **InterpointFields::GetReceiveVariable(const string &name, const int IDOfGhostLayer)
{
    return reinterpret_cast<T **> (interpointInformation->GetReceiveDataStorage(IDOfGhostLayer)->GetDataPtr(name));
}

InterpointDataProxy::InterpointDataProxy()
{
    ;
}

InterpointDataProxy::~InterpointDataProxy()
{
}

vector<RDouble **> &InterpointDataProxy::GetVectorData()
{
    return vdata;
}


vector<int> &InterpointDataProxy::GetVectorDimension()
{
    return vdim;
}

}