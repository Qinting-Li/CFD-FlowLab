inline int * InterfacePatch::GetFaceIndexForRecv() const
{
    return faceIndexForRecv;
}

inline void InterfacePatch::SetFaceIndexForRecv(int *faceIndexForRecv)
{
    this->faceIndexForRecv = faceIndexForRecv;
}

inline int * InterfacePatch::GetFaceIndexForSend() const
{
    return faceIndexForSend;
}

inline void InterfacePatch::SetFaceIndexForSend(int *faceIndexForSend)
{
    this->faceIndexForSend = faceIndexForSend;
}

inline int InterfacePatch::GetNumberOfFace() const
{
    return faceNumberOfPatch;
}

inline void InterfacePatch::SetNumberOfFace(int faceNumberOfPatch)
{
    this->faceNumberOfPatch = faceNumberOfPatch;
}

inline int InterfacePatch::GetZoneIndexOfNeighbor() const
{
    return zoneIndexOfNeighbor;
}

inline void InterfacePatch::SetZoneIndexOfNeighbor(int zoneIndexOfNeighbor)
{
    this->zoneIndexOfNeighbor = zoneIndexOfNeighbor;
}

inline int InterfacePatch::GetZoneID() const
{
    return zoneid;
}

inline void InterfaceInfo::SetNIFace(int nIFace)
{
    this->nIFace = nIFace;
}

inline int InterfaceInfo::GetNIFace() const
{
    return nIFace;
}

inline void InterfaceInfo::SetNumberOfNeighbor(int numberOfNeighbor)
{
    this->numberOfNeighbor = numberOfNeighbor;
}

inline int InterfaceInfo::GetNumberOfNeighbor() const
{
    return numberOfNeighbor;
}

inline vector <InterfacePatch *> InterfaceInfo::GetInterfacePatch() const
{
    return interFacePatch;
}

inline InterfacePatch * InterfaceInfo::GetInterfacePatch(int ineighbor) const
{
    return interFacePatch[ineighbor];
}

inline int * InterfaceInfo::GetInterFace2ZoneID() const
{
    return interFace2ZoneID;
}

inline int * InterfaceInfo::GetInterFace2CellID() const 
{ 
    return interFace2CellID; 
}

inline int * InterfaceInfo::GetInterFace2InterFaceID() const
{
    return interFace2InterFaceID;
}

inline int * InterfaceInfo::GetInterFace2BoundaryFace() const
{
    return interFace2BoundaryFace;
}

inline int * InterfaceInfo::GetInterFaceDirection() const
{
    return interFaceDirection;
}

inline void InterfaceInfo::SetIsNeighborInterfaceFound(bool *data)
{
    this->isNeighborInterfaceFound = data;
}

inline bool * InterfaceInfo::GetIsNeighborInterfaceFound()
{
    return isNeighborInterfaceFound;
}

inline int InterfaceInfo::GetZoneIndexOfNeighbor(int ineighbor) const
{
    InterfacePatch *interFacePatch = this->GetInterfacePatch(ineighbor);
    return interFacePatch->GetZoneIndexOfNeighbor();
}

inline int * InterfaceInfo::GetFaceIndexForRecv(int ineighbor) const
{
    InterfacePatch *interFacePatch = this->GetInterfacePatch(ineighbor);
    return interFacePatch->GetFaceIndexForRecv();
}

inline int * InterfaceInfo::GetFaceIndexForSend(int ineighbor) const
{
    InterfacePatch *interFacePatch = this->GetInterfacePatch(ineighbor);
    return interFacePatch->GetFaceIndexForSend();
}

inline Data_ParamFieldSuite * InterfaceInfo::GetSendDataStorage(int istore) const
{
    return &dsend[istore];
}

inline Data_ParamFieldSuite * InterfaceInfo::GetRecvDataStorage(int istore) const
{
    return &drecv[istore];
}

inline int InterfaceInfo::GetNIFaceOfNeighbor(int ineighbor) const
{
    InterfacePatch *interFacePatch = this->GetInterfacePatch(ineighbor);
    return interFacePatch->GetNumberOfFace();
}

inline uint_t InterfaceFields::Size()
{
    return variableNames.size();
}

inline int InterfaceFields::GetDim(const int iData)
{
    return variableDimensions[iData];
}

inline string & InterfaceFields::GetName(const int iData)
{
    return variableNames[iData];
}

inline int InterfaceFields::GetSolverIndex(const int iData)
{
    return solverIndex[iData];
}