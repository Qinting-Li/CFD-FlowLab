#include "Geo_Sample.h"
#include "IO_FileName.h"
#include "GridType.h"
#include "TK_Parse.h"
#include "PHIO.h"
#include "TK_Log.h"
#include "TK_Exit.h"
#include "Glb_Dimension.h"
#include "Geo_MultiGridInfo_Struct.h"
#ifdef PH_PARALLEL
#endif

#include "PHMpi.h"

using namespace std;
using namespace PHSPACE;
using namespace PHMPI;

namespace PHSPACE
{

ProbeData::ProbeData()
{
    probeToCellDistance =1.0e10;
    probeCellType = -1;
    probeCellZoneID = -1;
    probeCellID = -1;
    probeCellNI = -1;
    probeCellNJ = -1;
    probeCellNK = -1;
    probeGlobalID = -1;
    probeSurfaceID = -1;
    probeLineID = -1;
    isProbeValid = false;
    probeCoordinates.resize(0);
}

ProbeData::~ProbeData()
{
    probeCoordinates.clear();
}

void ProbeData::SetProbeCoordinates(const vector <RDouble> &probeCoordinates)
{
    this->probeCoordinates = probeCoordinates;
}

void ProbeData::SetProbeGlobalID(const int &probeGlobalID)
{
    this->probeGlobalID = probeGlobalID;
}

void ProbeData::SetProbeLineID(const int &probeLineID)
{
    this->probeLineID = probeLineID;
}

void ProbeData::SetProbeSurfaceID(const int &probeSurfaceID)
{
    this->probeSurfaceID = probeSurfaceID;
}

void ProbeData::SetProbeValidity(const bool &isProbeValid)
{
    this->isProbeValid = isProbeValid;
}

void ProbeData::SetProbeToCellDistance(const RDouble &probeToCellDistance)
{
    this->probeToCellDistance = probeToCellDistance;
}

void ProbeData::SetProbeCellType(const int &probeCellType)
{
    this->probeCellType = probeCellType;
}

void ProbeData::SetProbeCellZoneID(const int &probeCellZoneID)
{
    this->probeCellZoneID = probeCellZoneID;
}

void ProbeData::SetProbeCellID(const int &probeCellID)
{
    this->probeCellID = probeCellID;
}

void ProbeData::SetProbeCellNI(const int &probeCellNI)
{
    this->probeCellNI = probeCellNI;
}

void ProbeData::SetProbeCellNJ(const int &probeCellNJ)
{
    this->probeCellNJ = probeCellNJ;
}

void ProbeData::SetProbeCellNK(const int &probeCellNK)
{
    this->probeCellNK = probeCellNK;
}

vector <RDouble> ProbeData::GetProbeCoordinates()
{
    return this->probeCoordinates;
}

int ProbeData::GetProbeGlobalID()
{
    return this->probeGlobalID;
}

int ProbeData::GetProbeLineID()
{
    return this->probeLineID;
}

int ProbeData::GetProbeSurfaceID()
{
    return this->probeSurfaceID;
}

bool ProbeData::GetProbeValidity()
{
    return this->isProbeValid;
}

RDouble ProbeData::GetProbeToCellDistance()
{
    return this->probeToCellDistance;
}

int ProbeData::GetProbeCellType()
{
    return this->probeCellType;
}

int ProbeData::GetProbeCellZoneID()
{
    return this->probeCellZoneID;
}

int ProbeData::GetProbeCellID()
{
    return this->probeCellID;
}

int ProbeData::GetProbeCellNI()
{
    return this->probeCellNI;
}

int ProbeData::GetProbeCellNJ()
{
    return this->probeCellNJ;
}

int ProbeData::GetProbeCellNK()
{
    return this->probeCellNK;
}

void ProbeData::SearchNearestCell(Grid *grid)
{
    if (grid->Type() == PHSPACE::STRUCTGRID)
    {
        StructGrid * structGrid = StructGridCast(grid);
        SearchNearestCell(structGrid);
    }
    else
    {
        UnstructGrid * unstructGrid = UnstructGridCast(grid);
        SearchNearestCell(unstructGrid);
    }
}

void ProbeData::SearchNearestCell(StructGrid * grid)
{
    int gridtype = grid->Type();
    int zoneID = grid->GetGridID()->GetIndex();

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, 0, -1, I, J, K);

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    RDouble3D &xcc = *(grid->GetCellCenterX());
    RDouble3D &ycc = *(grid->GetCellCenterY());
    RDouble3D &zcc = *(grid->GetCellCenterZ());

    RDouble maxDistance = grid->CalMaxEdgeLength();

    vector <RDouble> coordinates = GetProbeCoordinates();
    bool isValid = GetProbeValidity();
    RDouble nearestCellDistance = GetProbeToCellDistance();
    int nearestCellZoneID = GetProbeCellZoneID();
    int nearestCellNI = GetProbeCellNI();
    int nearestCellNJ = GetProbeCellNJ();
    int nearestCellNK = GetProbeCellNK();

    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                RDouble CurrentDistance = DISTANCE(xcc(i, j, k) - coordinates[0], ycc(i, j, k) - coordinates[1], zcc(i, j, k) - coordinates[2]);
                if (CurrentDistance < nearestCellDistance)
                {
                    nearestCellDistance = CurrentDistance;
                    nearestCellNI = i;
                    nearestCellNJ = j;
                    nearestCellNK = k;
                    nearestCellZoneID = zoneID;
                    if (nearestCellDistance <= maxDistance) isValid = true;
                }
            }
        }
    }

    SetProbeCellType(gridtype);
    SetProbeValidity(isValid);
    SetProbeToCellDistance(nearestCellDistance);
    SetProbeCellZoneID(nearestCellZoneID);
    SetProbeCellNI(nearestCellNI);
    SetProbeCellNJ(nearestCellNJ);
    SetProbeCellNK(nearestCellNK);
}

void ProbeData::SearchNearestCell(UnstructGrid * grid)
{
    int gridtype = grid->Type();
    int zoneID = grid->GetGridID()->GetIndex();

    int nTotalCell = grid->GetNTotalCell();

    RDouble * xcc = grid->GetCellCenterX();
    RDouble * ycc = grid->GetCellCenterY();
    RDouble * zcc = grid->GetCellCenterZ();

    RDouble maxDistance = grid->CalMaxEdgeLength();

    vector <RDouble> coordinates = GetProbeCoordinates();
    bool isValid = GetProbeValidity();
    RDouble nearestCellDistance = GetProbeToCellDistance();
    int nearestCellID = GetProbeCellID();
    int nearestCellZoneID = GetProbeCellZoneID();

    for (int icell = 0; icell < nTotalCell; ++ icell)
    {
        RDouble CurrentDistance = DISTANCE(xcc[icell]-coordinates[0], ycc[icell]-coordinates[1], zcc[icell]-coordinates[2]);
        if (CurrentDistance < nearestCellDistance)
        {
            nearestCellDistance = CurrentDistance;
            nearestCellID = icell;
            nearestCellZoneID = zoneID;
            if (nearestCellDistance <= maxDistance) isValid = true;
        }
    }

    SetProbeCellType(gridtype);
    SetProbeValidity(isValid);
    SetProbeToCellDistance(nearestCellDistance);
    SetProbeCellZoneID(nearestCellZoneID);
    SetProbeCellID(nearestCellID);

}

void ProbeData::SearchRealCell(Grid *grid, bool &flag)
{
    if (grid->Type() == PHSPACE::STRUCTGRID)
    {
        StructGrid * structGrid = StructGridCast(grid);
        SearchRealCell(structGrid, flag);
    }
    else
    {
        UnstructGrid * unstructGrid = UnstructGridCast(grid);
        SearchRealCell(unstructGrid, flag);
    }
}

void ProbeData::SearchRealCell(StructGrid * grid, bool &flag)
{
    int gridtype = grid->Type();
    int zoneID = grid->GetGridID()->GetIndex();

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange(ni, nj, nk, 0, -1, I, J, K);

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    RDouble3D &xcc = *(grid->GetCellCenterX());
    RDouble3D &ycc = *(grid->GetCellCenterY());
    RDouble3D &zcc = *(grid->GetCellCenterZ());

    RDouble4D &xfn = *(grid->GetFaceNormalX());
    RDouble4D &yfn = *(grid->GetFaceNormalY());
    RDouble4D &zfn = *(grid->GetFaceNormalZ());

    vector <RDouble> coordinates = GetProbeCoordinates();
    RDouble realCellDistance = GetProbeToCellDistance();
    bool isValid = GetProbeValidity();
    int realCellZoneID = GetProbeCellZoneID();
    int realCellNI = GetProbeCellNI();
    int realCellNJ = GetProbeCellNJ();
    int realCellNK = GetProbeCellNK();

    int nDim = GetDim();
    int il, jl, kl;
    RDouble small = 1.0e-8;
    for (int k = kst; k <= ked; ++ k)
    {
        for (int j = jst; j <= jed; ++ j)
        {
            for (int i = ist; i <= ied; ++ i)
            {
                for (int iSurface = 1; iSurface <= nDim; ++ iSurface)
                { 
                    int il1, jl1, kl1;
                    GetNsurfIndex(iSurface, il1, jl1, kl1);

                    il = i + il1;
                    jl = j + jl1;
                    kl = k + kl1;


                    RDouble xfn1 = xfn(i, j, k, iSurface);
                    RDouble yfn1 = yfn(i, j, k, iSurface);
                    RDouble zfn1 = zfn(i, j, k, iSurface);

                    RDouble xfn2 = xfn(il, jl, kl, iSurface);
                    RDouble yfn2 = yfn(il, jl, kl, iSurface);
                    RDouble zfn2 = zfn(il, jl, kl, iSurface);

                    RDouble xfc1, yfc1, zfc1, xfc2, yfc2, zfc2;
                    grid->FaceCoor(i, j, k, iSurface, xfc1, yfc1, zfc1);

                    RDouble probeToFaceVecX1 = coordinates[0] - xfc1;
                    RDouble probeToFaceVecY1 = coordinates[1] - yfc1;
                    RDouble probeToFaceVecZ1 = coordinates[2] - zfc1;
                    RDouble probeToFaceDistance1 = DISTANCE(probeToFaceVecX1, probeToFaceVecY1, probeToFaceVecZ1);
                    RDouble dotProduct11 = probeToFaceVecX1 * xfn1 + probeToFaceVecY1 * yfn1 + probeToFaceVecZ1 * zfn1;

                    RDouble centerToFaceVecX1 = xcc(i, j, k) - xfc1;
                    RDouble centerToFaceVecY1 = ycc(i, j, k) - yfc1;
                    RDouble centerToFaceVecZ1 = zcc(i, j, k) - zfc1;
                    RDouble dotProduct21 = centerToFaceVecX1 * xfn1 + centerToFaceVecY1 * yfn1 + centerToFaceVecZ1 * zfn1;

                    grid->FaceCoor(il ,jl ,kl , iSurface, xfc2, yfc2, zfc2);

                    RDouble probeToFaceVecX2 = coordinates[0] - xfc2;
                    RDouble probeToFaceVecY2 = coordinates[1] - yfc2;
                    RDouble probeToFaceVecZ2 = coordinates[2] - zfc2;
                    RDouble probeToFaceDistance2 = DISTANCE(probeToFaceVecX2, probeToFaceVecY2, probeToFaceVecZ2);
                    RDouble dotProduct12 = probeToFaceVecX2 * xfn2 + probeToFaceVecY2 * yfn2 + probeToFaceVecZ2 * zfn2;


                    RDouble centerToFaceVecX2 = xcc(i, j, k) - xfc2;
                    RDouble centerToFaceVecY2 = ycc(i, j, k) - yfc2;
                    RDouble centerToFaceVecZ2 = zcc(i, j, k) - zfc2;
                    RDouble dotProduct22 = centerToFaceVecX2 * xfn2 + centerToFaceVecY2 * yfn2 + centerToFaceVecZ2 * zfn2;

                    // Judge if the probe and cell center at the same side of the current face.
                    // First face.

                    if (ABS(dotProduct11 / (probeToFaceDistance1 + TINY)) < small)
                    {
                        flag = true;
                    }
                    else if (dotProduct11 * dotProduct21 > zero)
                    {
                        flag = true;
                    }
                    else
                    {
                        flag = false;
                        break;
                    }

                    // Second face.
                    if (ABS(dotProduct12 / (probeToFaceDistance2 + TINY)) < small)
                    {
                        flag = true;
                    }
                    else if (dotProduct12 * dotProduct22 > zero)
                    {
                        flag = true;
                    }
                    else
                    {
                        flag = false;
                        break;
                    }
                }

                if (flag)
                {
                    isValid = true;
                    realCellNI = i;
                    realCellNJ = j;
                    realCellNK = k;
                    realCellZoneID = zoneID;
                    realCellDistance = DISTANCE(xcc(i, j, k) - coordinates[0], ycc(i, j, k) - coordinates[1], zcc(i, j, k) - coordinates[2]); 
                    break;
                }

            }

            if (flag) break;

        }

        if (flag) break;

    }

    SetProbeCellType(gridtype);
    SetProbeToCellDistance(realCellDistance);
    SetProbeCellZoneID(realCellZoneID);
    SetProbeCellNI(realCellNI);
    SetProbeCellNJ(realCellNJ);
    SetProbeCellNK(realCellNK);
    SetProbeValidity(isValid);
}

void ProbeData::SearchRealCell(UnstructGrid * grid, bool &flag)
{
    int gridtype = grid->Type();
    int zoneID = grid->GetGridID()->GetIndex();

    int nTotalCell = grid->GetNTotalCell();

    RDouble *xcc = grid->GetCellCenterX();
    RDouble *ycc = grid->GetCellCenterY();
    RDouble *zcc = grid->GetCellCenterZ();

    RDouble *xfn = grid->GetFaceNormalX();
    RDouble *yfn = grid->GetFaceNormalY();
    RDouble *zfn = grid->GetFaceNormalZ();

    RDouble *xfc = grid->GetFaceCenterX();
    RDouble *yfc = grid->GetFaceCenterY();
    RDouble *zfc = grid->GetFaceCenterZ();

    int ** cell2face = grid->GetCell2Face();
    int * faceNumberOfEachCell = grid->GetFaceNumberOfEachCell();

    vector <RDouble> coordinates = GetProbeCoordinates();
    RDouble realCellDistance = GetProbeToCellDistance();
    bool isValid = GetProbeValidity();
    int realCellID = GetProbeCellID();
    int realCellZoneID = GetProbeCellZoneID();

    RDouble small = 1.0e-8;
    for (int icell = 0; icell < nTotalCell; ++ icell)
    {
        for (int iface = 0; iface < faceNumberOfEachCell[icell]; ++ iface)
        {
            int faceID = cell2face[icell][iface];

            RDouble probeToFaceVecX = coordinates[0] - xfc[faceID];
            RDouble probeToFaceVecY = coordinates[1] - yfc[faceID];
            RDouble probeToFaceVecZ = coordinates[2] - zfc[faceID];
            RDouble probeToFaceDistance = DISTANCE(probeToFaceVecX, probeToFaceVecY, probeToFaceVecZ);
            RDouble dotProduct1 = probeToFaceVecX * xfn[faceID] + probeToFaceVecY * yfn[faceID] + probeToFaceVecZ * zfn[faceID];

            RDouble centerToFaceVecX = xcc[icell] - xfc[faceID];
            RDouble centerToFaceVecY = ycc[icell] - yfc[faceID];
            RDouble centerToFaceVecZ = zcc[icell] - zfc[faceID];
            RDouble dotProduct2 = centerToFaceVecX * xfn[faceID] + centerToFaceVecY * yfn[faceID] + centerToFaceVecZ * zfn[faceID];

            // Judge if the probe and cell center at the same side of the current face.
            if (ABS(dotProduct1 / (probeToFaceDistance + TINY)) < small)
            {
                flag = true;
            }
            else if (dotProduct1 * dotProduct2 > zero)
            {
                flag = true;
            }
            else
            {
                flag = false;
                break;
            }
        }

        if (flag) 
        {
            isValid = true;
            realCellID = icell;
            realCellZoneID = zoneID;
            realCellDistance = DISTANCE(xcc[realCellID] - coordinates[0], ycc[realCellID] - coordinates[1], zcc[realCellID] - coordinates[2]);
            break;
        }

    }

    SetProbeCellType(gridtype);
    SetProbeToCellDistance(realCellDistance);
    SetProbeCellZoneID(realCellZoneID);
    SetProbeCellID(realCellID);
    SetProbeValidity(isValid);

}

void ProbeData::WriteStructProbeCellInfo(fstream &file)
{
    int dataMonitorType = GlobalDataBase::GetIntParaFromDB("dataMonitorType");
    ostringstream oss;

    oss << setiosflags(ios::left);
    oss << setiosflags(ios::scientific);
    oss << setprecision(10);

    vector <RDouble> coordinates = GetProbeCoordinates();
    int probeID = GetProbeGlobalID();
    int lineID = GetProbeLineID();
    int surfaceID = GetProbeSurfaceID();
    //int probeNumber = probeID + 1;

    RDouble distance = GetProbeToCellDistance();
    int cellNI = GetProbeCellNI();
    int cellNJ = GetProbeCellNJ();
    int cellNK = GetProbeCellNK();

    int cellZoneID   = GetProbeCellZoneID();
    int cellType = GetProbeCellType(); 

    int wordWidth = 7;

    oss << setw(wordWidth) << probeID << "	";
    if (dataMonitorType == LINESMONITOR)
    {
        oss << setw(wordWidth) << lineID << "	";
    }
    if (dataMonitorType == SURFACESMONITOR)
    {
        oss << setw(wordWidth) << surfaceID << "\t";
    }
    oss << setw(wordWidth) << distance << "\t";
    oss << " ( " << cellNI << " , " << cellNJ << " , " << cellNK << " ) " << "\t";
    oss << setw(wordWidth) << cellType << "\t";
    oss << setw(wordWidth) << cellZoneID << "\t";
    oss << "\n";

    WriteASCIIFile(file, oss.str());
}

void ProbeData::WriteUnstructProbeCellInfo(fstream &file)
{
    int dataMonitorType = GlobalDataBase::GetIntParaFromDB("dataMonitorType");
    ostringstream oss;

    oss << setiosflags(ios::left);
    oss << setiosflags(ios::scientific);
    oss << setprecision(10);

    
    vector <RDouble> coordinates = GetProbeCoordinates();
    int probeID = GetProbeGlobalID();
    int lineID = GetProbeLineID();
    int surfaceID = GetProbeSurfaceID();

    RDouble distance = GetProbeToCellDistance();
    int cellID       = GetProbeCellID();
    int cellZoneID   = GetProbeCellZoneID();
    int cellType = GetProbeCellType();
    
    int wordWidth = 7;
    oss << setw(wordWidth) << probeID << "	";
    if (dataMonitorType == LINESMONITOR)
    {
        oss << setw(wordWidth) << lineID << "	";
    }
    if (dataMonitorType == SURFACESMONITOR)
    {
        oss << setw(wordWidth) << surfaceID << "	";
    }
    oss << setw(wordWidth) << distance << "	";
    oss << setw(wordWidth) << cellID << "	";
    oss << setw(wordWidth) << cellType << "	";
    oss << setw(wordWidth) << cellZoneID << "	";
    oss << "\n";

    WriteASCIIFile(file, oss.str());
}

void ProbeData::CompressUnstructProbeCellInfo(DataContainer *cdata)
{
    vector <RDouble> coordinates = GetProbeCoordinates();
    int probeID = GetProbeGlobalID();
    int lineID = GetProbeLineID();
    int surfaceID = GetProbeSurfaceID();
    bool isValid = GetProbeValidity();
    int cellType = GetProbeCellType();
    int cellZoneID = GetProbeCellZoneID();
    int cellID = GetProbeCellID();
    RDouble distance = GetProbeToCellDistance();

    int nDim =static_cast<int>(coordinates.size());
    PHWrite(cdata, nDim);
    for (int iDim = 0; iDim < nDim; ++ iDim)
    {
        PHWrite(cdata,coordinates[iDim]);
    }

    PHWrite(cdata, probeID);
    PHWrite(cdata, lineID);
    PHWrite(cdata, surfaceID);
    PHWrite(cdata, isValid);
    PHWrite(cdata, cellType);
    PHWrite(cdata, cellZoneID);
    PHWrite(cdata, cellID);
    PHWrite(cdata, distance);
}

void ProbeData::CompressStructProbeCellInfo(DataContainer *cdata)
{
    vector <RDouble> coordinates = GetProbeCoordinates();
    int probeID = GetProbeGlobalID();
    int lineID = GetProbeLineID();
    int surfaceID = GetProbeSurfaceID();
    bool isValid = GetProbeValidity();
    int cellType = GetProbeCellType();
    int cellZoneID = GetProbeCellZoneID();
    int cellNI = GetProbeCellNI();
    int cellNJ = GetProbeCellNJ();
    int cellNK = GetProbeCellNK();
    RDouble distance = GetProbeToCellDistance();

    int nDim =static_cast<int>(coordinates.size());
    PHWrite(cdata, nDim);

    for (int iDim = 0; iDim < nDim; ++ iDim)
    {
        PHWrite(cdata,coordinates[iDim]);
    }
    PHWrite(cdata, probeID);
    PHWrite(cdata, lineID);
    PHWrite(cdata, surfaceID);
    PHWrite(cdata, isValid);
    PHWrite(cdata, cellType);
    PHWrite(cdata, cellZoneID);
    PHWrite(cdata, cellNI);
    PHWrite(cdata, cellNJ);
    PHWrite(cdata, cellNK);
    PHWrite(cdata, distance);

}

void ProbeData::DecompressStructProbeCellInfo(DataContainer *cdata)
{
    RDouble distanceRecv;
    int cellZoneIDRecv, cellTypeRecv, cellNIRecv, cellNJRecv, cellNKRecv;
    vector <RDouble> coordinatesRecv;
    int nDimRecv, probeIDRecv, lineIDRecv, surfaceIDRecv;
    bool isValidRecv;

    PHRead(cdata, nDimRecv);
    coordinatesRecv.resize(nDimRecv);

    for (int iDim = 0; iDim < nDimRecv; ++ iDim)
    {
        PHRead(cdata,coordinatesRecv[iDim]);
    }

    PHRead(cdata, probeIDRecv);
    PHRead(cdata, lineIDRecv);
    PHRead(cdata, surfaceIDRecv);
    PHRead(cdata, isValidRecv);
    PHRead(cdata, cellTypeRecv);
    PHRead(cdata, cellZoneIDRecv);
    PHRead(cdata, cellNIRecv);
    PHRead(cdata, cellNJRecv);
    PHRead(cdata, cellNKRecv);
    PHRead(cdata, distanceRecv);

    SetProbeCoordinates(coordinatesRecv);
    SetProbeGlobalID(probeIDRecv);
    SetProbeLineID(lineIDRecv);
    SetProbeSurfaceID(surfaceIDRecv);
    SetProbeValidity(isValidRecv);
    SetProbeCellType(cellTypeRecv);
    SetProbeCellZoneID(cellZoneIDRecv);
    SetProbeCellNI(cellNIRecv);
    SetProbeCellNJ(cellNJRecv);
    SetProbeCellNK(cellNKRecv);
    SetProbeToCellDistance(distanceRecv);
}

void ProbeData::DecompressUnstructProbeCellInfo(DataContainer *cdata)
{
    RDouble distanceRecv;
    int cellZoneIDRecv, cellTypeRecv, cellIDRecv;
    vector <RDouble> coordinatesRecv;
    int nDimRecv, probeIDRecv, lineIDRecv, surfaceIDRecv;
    bool isValidRecv;

    PHRead(cdata, nDimRecv);
    coordinatesRecv.resize(nDimRecv);
    for (int iDim = 0; iDim < nDimRecv; ++ iDim)
    {
        PHRead(cdata,coordinatesRecv[iDim]);
    }

    PHRead(cdata, probeIDRecv);
    PHRead(cdata, lineIDRecv);
    PHRead(cdata, surfaceIDRecv);
    PHRead(cdata, isValidRecv);
    PHRead(cdata, cellTypeRecv);
    PHRead(cdata, cellZoneIDRecv);
    PHRead(cdata, cellIDRecv);
    PHRead(cdata, distanceRecv);


    SetProbeCoordinates(coordinatesRecv);
    SetProbeGlobalID(probeIDRecv);
    SetProbeLineID(lineIDRecv);
    SetProbeSurfaceID(surfaceIDRecv);
    SetProbeValidity(isValidRecv);
    SetProbeCellType(cellTypeRecv);
    SetProbeCellZoneID(cellZoneIDRecv);
    SetProbeCellID(cellIDRecv);
    SetProbeToCellDistance(distanceRecv);
}

SampleLocationInfo::SampleLocationInfo()
{

    numberToMonitor = 0;
    numberOfTotalProbes    = 0;
    numberOfProbes.resize(0);
    surfaceProbes_ni.resize(0);
    surfaceProbes_nj.resize(0);
    typeOfTotalProbes.resize(0);
    sampleLocationInfo.resize(0);
    defineFileExist = false;
}

SampleLocationInfo::~SampleLocationInfo()
{
    for (vector<ProbeData *>::iterator iter = sampleLocationInfo.begin(); iter != sampleLocationInfo.end(); ++ iter)
    {
        if (NULL != *iter)
        {
            delete *iter;
            *iter = NULL;
        }
    }
    sampleLocationInfo.clear();
}

void SampleLocationInfo::Run()
{
    int dataMonitorType = GlobalDataBase::GetIntParaFromDB("dataMonitorType");

    if (dataMonitorType == PROBESMONITOR)
    {
        string probesDefineFile = GlobalDataBase::GetStrParaFromDB("probesDefineFile");
        SetDefineFileName(probesDefineFile);

        CheckIfDefineFileExist();

        if (IfDefineFileExist())
        {
            ReadProbesData();
        }
        else
        {
            TK_Exit::ExceptionExit("Error: the file used to define probes location information does not exist !\n");
            return;
        }
    }
    else if (dataMonitorType == LINESMONITOR)
    {
        string linesDefineFile = GlobalDataBase::GetStrParaFromDB("linesDefineFile");
        SetDefineFileName(linesDefineFile);

        CheckIfDefineFileExist();

        if (IfDefineFileExist())
        {
            int nLines = GlobalDataBase::GetIntParaFromDB("nLines");
            SetNumberToMonitor(nLines);

            ReadLinesData();
        }
        else
        {
            TK_Exit::ExceptionExit("Error: the file used to define lines location information does not exist !\n");
            return;
        }
    }
    else if (dataMonitorType == SURFACESMONITOR)
    {
        string surfacesDefineFile = GlobalDataBase::GetStrParaFromDB("surfacesDefineFile");
        SetDefineFileName(surfacesDefineFile);

        CheckIfDefineFileExist();
        if (IfDefineFileExist())
        {
            int nSurfaces = GlobalDataBase::GetIntParaFromDB("nSurfaces");
            SetNumberToMonitor(nSurfaces);

            ReadSurfacesData();
        }
        else
        {
            TK_Exit::ExceptionExit("Error: the file used to define surfaces location information does not exist !\n");
            return;
        }
    }
    else
    {
        TK_Exit::UnexpectedVarValue("dataMonitorType", dataMonitorType);
    }

    UpdateAllProbesData();

    UpdateProbesGlobalPara();

    ServerCollectionProbesData();

    int currentProcessorID = GetCurrentProcessorID();
    if (currentProcessorID == server)
    {
        CheckInvalidProbesData();

        string gridfile = GlobalDataBase::GetStrParaFromDB("gridfile");
        string probesCellFileName = PHSPACE::AddSymbolToFileName(gridfile, '_', currentProcessorID);
        WriteProbesCellInfo(probesCellFileName);
    }

    BroadCastProbesData();

    SampleLocationInfoToGrid();
}

void SampleLocationInfo::SetNumberToMonitor(const int &numberToMonitor)
{
    this->numberToMonitor = numberToMonitor;
}

int SampleLocationInfo::GetNumberToMonitor()
{
    return this->numberToMonitor;
}

void SampleLocationInfo::SetDefineFileName(const string &defineFileName)
{
    this->defineFileName = defineFileName;
}
void SampleLocationInfo::SetTotalProbesNumber(const int &numberOfTotalProbes)
{
    this->numberOfTotalProbes = numberOfTotalProbes;
}

void SampleLocationInfo::SetProbesNumber(const vector <int> &numberOfProbes)
{
    this->numberOfProbes = numberOfProbes;
}

void SampleLocationInfo::SetSurfaceProbesNI(const vector <int> &surfaceProbes_ni)
{
    this->surfaceProbes_ni = surfaceProbes_ni;
}

void SampleLocationInfo::SetSurfaceProbesNJ(const vector <int> &surfaceProbes_nj)
{
    this->surfaceProbes_nj = surfaceProbes_nj;
}

void SampleLocationInfo::SetTotalProbesType(const vector <int> &typeOfTotalProbes)
{
    this->typeOfTotalProbes = typeOfTotalProbes;
}

string SampleLocationInfo::GetDefineFileName()
{
    return this->defineFileName;
}

int SampleLocationInfo::GetTotalProbesNumber()
{
    return this->numberOfTotalProbes;
}

vector <int> SampleLocationInfo::GetProbesNumber()
{
    return this->numberOfProbes;
}

vector <int> SampleLocationInfo::GetSurfaceProbesNI()
{
    return this->surfaceProbes_ni;
}

vector <int> SampleLocationInfo::GetSurfaceProbesNJ()
{
    return this->surfaceProbes_nj;
}

vector <int> SampleLocationInfo::GetTotalProbesType()
{
    return this->typeOfTotalProbes;
}

vector <ProbeData*> SampleLocationInfo::GetSampleLocationInfo()
{
    return this->sampleLocationInfo;
}

bool SampleLocationInfo::IfDefineFileExist()
{
    return defineFileExist;
}

void SampleLocationInfo::Nondimensional()
{
    RDouble gridScaleFactor = GlobalDataBase::GetDoubleParaFromDB("gridScaleFactor");
    for (int iProbe = 0; iProbe < numberOfTotalProbes; ++ iProbe)
    {
        vector <RDouble> coordinates = sampleLocationInfo[iProbe]->GetProbeCoordinates();

        for (int iDim = 0; iDim < coordinates.size(); ++ iDim)
        {
            coordinates[iDim] *= gridScaleFactor;
        }

        sampleLocationInfo[iProbe]->SetProbeCoordinates(coordinates);
    }
}

void SampleLocationInfo::CheckIfDefineFileExist()
{
    defineFileExist = false;

    fstream file;
    string keyFile = GetDefineFileName();
    file.open(keyFile.c_str(), ios_base::in);

    if (!file)
    {
        return;
    }

    defineFileExist = true;
}

void SampleLocationInfo::ReadProbesData()
{
    fstream file;
    string keyFile = GetDefineFileName();
    file.open(keyFile.c_str(), ios_base::in);

    string line, word;
    string separator  = " =\t\r\n#$,;\"'";
    string sep1 = "\r\n";

    while (!file.eof())
    {
        getline(file, line);

        if (line == "") continue;
        FindNextWord(line, word, sep1);
        if (word.substr(0, 1) == "#" || word.substr(0, 2) == "//") continue;
        line = FindNextWord(line, word, separator);

        int nProbesIn;
        from_string <int> (nProbesIn, word, std::dec);
        numberOfTotalProbes = nProbesIn;
        break;
    }

    const int nDim = 3;
    vector<vector<RDouble> > coordinatesList;

    while (!file.eof())
    {
        getline(file, line);

        if (line == "") continue;
        FindNextWord(line, word, sep1);
        if (word.substr(0, 1) == "#" || word.substr(0, 2) == "//") continue;

        vector <RDouble> coordinates;
        for (int iDim = 0; iDim < nDim; ++ iDim)
        {
            line = FindNextWord(line, word, separator);
            RDouble coordinatesIn;
            from_string <RDouble> (coordinatesIn, word, std::dec);
            coordinates.push_back(coordinatesIn);
        }
        coordinatesList.push_back(coordinates);
    }

    if (coordinatesList.size() == numberOfTotalProbes)
    {
        for (int iProbe = 0; iProbe < numberOfTotalProbes; ++ iProbe)
        {
            ProbeData *probe = new ProbeData();

            probe->SetProbeCoordinates(coordinatesList[iProbe]);
            probe->SetProbeGlobalID(iProbe);

            sampleLocationInfo.push_back(probe);
        }
    }
    else
    {
        ostringstream oss;
        oss<<"Error: the total number of monitored probes does not match the number of coordinates !" <<endl;
        TK_Exit::ExceptionExit(oss);
    }

    file.close();
    file.clear();
}

void SampleLocationInfo::ReadLinesData()
{
    fstream file;
    string keyFile = GetDefineFileName();
    file.open(keyFile.c_str(), ios_base::in);

    string line, word;
    string separator  = " =\t\r\n#$,;\"'";

    const int nDim = 3;

    vector<vector<RDouble> > coordinatesList;
    vector <int> probeIDList;
    vector <int> lineIDList;
    for (int iLine = 0; iLine < numberToMonitor; ++ iLine)
    {
        getline(file, line);

        line = FindNextWord(line, word, separator);
        int nProbesIn;
        from_string <int> (nProbesIn, word, std::dec);
        numberOfProbes.push_back(nProbesIn);

        for (int iProbe = 0; iProbe < nProbesIn; ++ iProbe)
        {
            getline(file, line);

            vector <RDouble> coordinates;
            for (int iDim = 0; iDim < nDim; ++ iDim)
            {
                line = FindNextWord(line, word, separator);
                RDouble coordinatesIn;
                from_string <RDouble> (coordinatesIn, word, std::dec);
                coordinates.push_back(coordinatesIn);
            }
            coordinatesList.push_back(coordinates);
            probeIDList.push_back(iProbe);
            lineIDList.push_back(iLine);
        }

        numberOfTotalProbes += nProbesIn;
    }

    if (coordinatesList.size() == numberOfTotalProbes)
    {
        for (int iProbe = 0; iProbe < numberOfTotalProbes; ++ iProbe)
        {
            ProbeData *probe = new ProbeData();

            probe->SetProbeCoordinates(coordinatesList[iProbe]);
            probe->SetProbeGlobalID(probeIDList[iProbe]);
            probe->SetProbeLineID(lineIDList[iProbe]);

            sampleLocationInfo.push_back(probe);
        }
    }
    else
    {
        ostringstream oss;
        oss<<"Error: the total number of monitored probes does not match the number of coordinates !" <<endl;
        TK_Exit::ExceptionExit(oss);
    }

    file.close();
    file.clear();
}

void SampleLocationInfo::ReadSurfacesData()
{
    fstream file;
    string keyFile = GetDefineFileName();
    file.open(keyFile.c_str(), ios_base::in);

    string line, word;
    string separator  = " =\t\r\n#$,;\"'";
    string sep1 = "\r\n";

    const int nDim = 3;

    vector<vector<RDouble> > coordinatesList;
    vector <int> probeIDList;
    vector <int> surfaceIDList;

    for (int iSurface = 0; iSurface < numberToMonitor; ++ iSurface)
    {
        getline(file, line);
        line = FindNextWord(line, word, separator);
        int nProbesIn_i;
        from_string <int> (nProbesIn_i, word, std::dec);
        surfaceProbes_ni.push_back(nProbesIn_i);

        int nProbesIn_j;
        line = FindNextWord(line, word, separator);
        from_string <int> (nProbesIn_j, word, std::dec);
        surfaceProbes_nj.push_back(nProbesIn_j);

        int nProbesIn = nProbesIn_i * nProbesIn_j;
        numberOfProbes.push_back(nProbesIn);

        for (int iProbe = 0; iProbe < nProbesIn; ++ iProbe)
        {
            getline(file, line);
            vector <RDouble> coordinates;
            for (int iDim = 0; iDim < nDim; ++ iDim)
            {
                line = FindNextWord(line, word, separator);
                RDouble coordinatesIn;
                from_string <RDouble> (coordinatesIn, word, std::dec);
                coordinates.push_back(coordinatesIn);
            }
            coordinatesList.push_back(coordinates);
            probeIDList.push_back(iProbe);
            surfaceIDList.push_back(iSurface);
        }

        numberOfTotalProbes += nProbesIn;
    }

    if (coordinatesList.size() == numberOfTotalProbes)
    {
        for (int iProbe = 0; iProbe < numberOfTotalProbes; ++ iProbe)
        {
            ProbeData *probe = new ProbeData();

            probe->SetProbeCoordinates(coordinatesList[iProbe]);
            probe->SetProbeGlobalID(probeIDList[iProbe]);
            probe->SetProbeSurfaceID(surfaceIDList[iProbe]);

            sampleLocationInfo.push_back(probe);
        }
    }
    else
    {
        ostringstream oss;
        oss<<"Error: the total number of monitored probes does not match the number of coordinates !" <<endl;
        TK_Exit::ExceptionExit(oss);
    }

    file.close();
    file.clear();
}

void SampleLocationInfo::UpdateAllProbesData()
{
    int searchCellsMethod = GlobalDataBase::GetIntParaFromDB("searchCellsMethod");

    if (searchCellsMethod == NEARESTCELLDATA)
    {
        UpdateNearestCellsData();
    }
    else if (searchCellsMethod == REALCELLDATA)
    {
        UpdateRealCellsData();
    }
    else
    {
        TK_Exit::UnexpectedVarValue("searchCellsMethod", searchCellsMethod);
    }

}

void SampleLocationInfo::UpdateNearestCellsData()
{
    int nLocalZones = GetNumberofLocalZones();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (int iProbe = 0; iProbe < numberOfTotalProbes; ++ iProbe)
    {
        vector <RDouble> coordinates = sampleLocationInfo[iProbe]->GetProbeCoordinates();
        for (int izone = 0; izone < nLocalZones; ++ izone)
        {
            int zoneID = GetLocalZoneIDToGlobalZoneID(izone);

            Grid *grid = GetGrid(zoneID, 0);

            RDouble *pmin = grid->GetMinBox();
            RDouble *pmax = grid->GetMaxBox();

            bool ProbeInMinMaxBox = CheckIfProbeInMinMaxBox(coordinates[0],coordinates[1],coordinates[2],pmin,pmax);
            if (!ProbeInMinMaxBox) continue;
            sampleLocationInfo[iProbe]->SearchNearestCell(grid);
        }
    }
}

void SampleLocationInfo::UpdateRealCellsData()
{
    int nLocalZones = GetNumberofLocalZones();

#ifdef USE_OMP
#pragma omp parallel for
#endif

    for (int iProbe = 0; iProbe < numberOfTotalProbes; ++ iProbe)
    {
        bool flag = false;
        vector <RDouble> coordinates = sampleLocationInfo[iProbe]->GetProbeCoordinates();

        for (int izone = 0; izone < nLocalZones; ++ izone)
        {
            int zoneID = GetLocalZoneIDToGlobalZoneID(izone);

            Grid *grid = GetGrid(zoneID, 0);

            RDouble *pmin = grid->GetMinBox();
            RDouble *pmax = grid->GetMaxBox();

            bool ProbeInMinMaxBox = CheckIfProbeInMinMaxBox(coordinates[0],coordinates[1],coordinates[2],pmin,pmax);
            if (!ProbeInMinMaxBox) continue;
            sampleLocationInfo[iProbe]->SearchRealCell(grid, flag);
            if (flag) break;
        }

    }
}

void SampleLocationInfo::UpdateProbesGlobalPara()
{
    int dataMonitorType = GlobalDataBase::GetIntParaFromDB("dataMonitorType");
    int nTotalProbes = GetTotalProbesNumber();
    vector <int> np = GetProbesNumber();
    GlobalDataBase::UpdateData("nTotalProbes", &nTotalProbes, PHINT, 1);

    if (dataMonitorType == LINESMONITOR)
    {
        int nLines = GetNumberToMonitor();
        GlobalDataBase::UpdateData("nLines", &nLines, PHINT, 1);

        int *nProbesOfLine = new int [nLines];
        for (int iLine = 0; iLine < nLines; ++ iLine)
        {
            nProbesOfLine[iLine] = np[iLine];
        }
        GlobalDataBase::UpdateDataPtr("nProbesOfLine", nProbesOfLine);
    }

    if (dataMonitorType == SURFACESMONITOR)
    {
        int nSurface = GetNumberToMonitor();
        GlobalDataBase::UpdateData("nSurface", &nSurface, PHINT, 1);

        vector <int> np = GetProbesNumber();

        int *nProbesOfSurface = new int [nSurface];

        for (int iSurface = 0; iSurface < nSurface; ++ iSurface)
        {
            nProbesOfSurface[iSurface] = np[iSurface];
        }

        GlobalDataBase::UpdateDataPtr("nProbesOfSurface", nProbesOfSurface);
    }


}

void SampleLocationInfo::ServerCollectionProbesData()
{
    int searchCellsMethod = GlobalDataBase::GetIntParaFromDB("searchCellsMethod");
    int server = GetServerProcessorID();
    int myid   = GetCurrentProcessorID();
    int numberOfProcessor = GetNumberOfProcessor();
    DataContainer *cdata = new DataContainer();
    int nTotalProbes = GetTotalProbesNumber();

    cdata->moveto_begin();
    RDouble distance;
    int cellZoneID, cellType, cellID, cellNI, cellNJ, cellNK;
    vector <RDouble> coordinates;
    int nDim, probeID, lineID, surfaceID;
    bool isValid;

    for (int iProbe = 0; iProbe < nTotalProbes; ++ iProbe)
    {
        coordinates = sampleLocationInfo[iProbe]->GetProbeCoordinates();
        probeID = sampleLocationInfo[iProbe]->GetProbeGlobalID();
        lineID = sampleLocationInfo[iProbe]->GetProbeLineID();
        surfaceID = sampleLocationInfo[iProbe]->GetProbeSurfaceID();
        isValid = sampleLocationInfo[iProbe]->GetProbeValidity();
        cellType = sampleLocationInfo[iProbe]->GetProbeCellType();
        cellZoneID = sampleLocationInfo[iProbe]->GetProbeCellZoneID();
        cellID = sampleLocationInfo[iProbe]->GetProbeCellID();
        cellNI = sampleLocationInfo[iProbe]->GetProbeCellNI();
        cellNJ = sampleLocationInfo[iProbe]->GetProbeCellNJ();
        cellNK = sampleLocationInfo[iProbe]->GetProbeCellNK();
        distance = sampleLocationInfo[iProbe]->GetProbeToCellDistance();

        nDim =static_cast<int>(coordinates.size());

        PHWrite(cdata, myid);
        PHWrite(cdata, nDim);
        for (int iDim = 0; iDim < nDim; ++ iDim)
        {
            PHWrite(cdata,coordinates[iDim]);
        }
        PHWrite(cdata, probeID);
        PHWrite(cdata, lineID);
        PHWrite(cdata, surfaceID);
        PHWrite(cdata, isValid);
        PHWrite(cdata, cellType);
        PHWrite(cdata, cellZoneID);
        PHWrite(cdata, cellID);
        PHWrite(cdata, cellNI);
        PHWrite(cdata, cellNJ);
        PHWrite(cdata, cellNK);
        PHWrite(cdata, distance);
    }

    if (myid != server)
    {
        int tag = myid;
        send(cdata, server, tag);
        delete cdata;
    }
    else
    {
        for (int iProcessor = 0; iProcessor < numberOfProcessor; ++ iProcessor)
        {
            int tag = iProcessor;

            if (server != iProcessor)
            {
                if (cdata != NULL)
                {
                    delete cdata;
                }
                cdata = new DataContainer();
                receive(cdata, iProcessor, tag);
            }
            else
            {
                continue;
            }
            int nTotalProbesOfProcessor = GetTotalProbesNumber();
            cdata->moveto_begin();

            int processorOfminDisResRecv = 0;
            RDouble distanceRecv;
            int cellZoneIDRecv, cellTypeRecv, cellIDRecv, cellNIRecv, cellNJRecv, cellNKRecv;
            vector <RDouble> coordinatesRecv;
            int nDimRecv,probeIDRecv, lineIDRecv,surfaceIDRecv;
            bool isValidRecv;
            for (int iProbe = 0; iProbe < nTotalProbesOfProcessor; ++ iProbe)
            {
                PHRead(cdata, processorOfminDisResRecv);
                PHRead(cdata,nDimRecv);
                coordinatesRecv.resize(nDimRecv);
                for (int iDim = 0; iDim < nDimRecv; ++ iDim)
                {
                    PHRead(cdata,coordinatesRecv[iDim]);
                }
                PHRead(cdata, probeIDRecv);
                PHRead(cdata, lineIDRecv);
                PHRead(cdata, surfaceIDRecv);
                PHRead(cdata, isValidRecv);
                PHRead(cdata, cellTypeRecv);
                PHRead(cdata, cellZoneIDRecv);
                PHRead(cdata, cellIDRecv);
                PHRead(cdata, cellNIRecv);
                PHRead(cdata, cellNJRecv);
                PHRead(cdata, cellNKRecv);
                PHRead(cdata, distanceRecv);

                double distanceTemp = sampleLocationInfo[iProbe]->GetProbeToCellDistance();
                int cellZoneIDTemp = sampleLocationInfo[iProbe]->GetProbeCellZoneID();

                if(searchCellsMethod == 0)
                {
                    if (distanceRecv < distanceTemp)
                    {
                        sampleLocationInfo[iProbe]->SetProbeCoordinates(coordinatesRecv);
                        sampleLocationInfo[iProbe]->SetProbeGlobalID(probeIDRecv);
                        sampleLocationInfo[iProbe]->SetProbeLineID(lineIDRecv);
                        sampleLocationInfo[iProbe]->SetProbeSurfaceID(surfaceIDRecv);
                        sampleLocationInfo[iProbe]->SetProbeValidity(isValidRecv);
                        sampleLocationInfo[iProbe]->SetProbeCellType(cellTypeRecv);
                        sampleLocationInfo[iProbe]->SetProbeCellZoneID(cellZoneIDRecv);
                        sampleLocationInfo[iProbe]->SetProbeCellID(cellIDRecv);
                        sampleLocationInfo[iProbe]->SetProbeCellNI(cellNIRecv);
                        sampleLocationInfo[iProbe]->SetProbeCellNJ(cellNJRecv);
                        sampleLocationInfo[iProbe]->SetProbeCellNK(cellNKRecv);
                        sampleLocationInfo[iProbe]->SetProbeToCellDistance(distanceRecv);
                    }
                }
                else
                {
                    if (cellZoneIDTemp < 0)
                    {
                        if(cellZoneIDRecv > 0)
                        {
                            sampleLocationInfo[iProbe]->SetProbeCoordinates(coordinatesRecv);
                            sampleLocationInfo[iProbe]->SetProbeGlobalID(probeIDRecv);
                            sampleLocationInfo[iProbe]->SetProbeLineID(lineIDRecv);
                            sampleLocationInfo[iProbe]->SetProbeSurfaceID(surfaceIDRecv);
                            sampleLocationInfo[iProbe]->SetProbeValidity(isValidRecv);
                            sampleLocationInfo[iProbe]->SetProbeCellType(cellTypeRecv);
                            sampleLocationInfo[iProbe]->SetProbeCellZoneID(cellZoneIDRecv);
                            sampleLocationInfo[iProbe]->SetProbeCellID(cellIDRecv);
                            sampleLocationInfo[iProbe]->SetProbeCellNI(cellNIRecv);
                            sampleLocationInfo[iProbe]->SetProbeCellNJ(cellNJRecv);
                            sampleLocationInfo[iProbe]->SetProbeCellNK(cellNKRecv);
                            sampleLocationInfo[iProbe]->SetProbeToCellDistance(distanceRecv);
                        }
                    }
                    else
                    {
                        if(cellZoneIDRecv > 0 && cellZoneIDRecv < cellZoneIDTemp)
                        {
                            sampleLocationInfo[iProbe]->SetProbeCoordinates(coordinatesRecv);
                            sampleLocationInfo[iProbe]->SetProbeGlobalID(probeIDRecv);
                            sampleLocationInfo[iProbe]->SetProbeLineID(lineIDRecv);
                            sampleLocationInfo[iProbe]->SetProbeSurfaceID(surfaceIDRecv);
                            sampleLocationInfo[iProbe]->SetProbeValidity(isValidRecv);
                            sampleLocationInfo[iProbe]->SetProbeCellType(cellTypeRecv);
                            sampleLocationInfo[iProbe]->SetProbeCellZoneID(cellZoneIDRecv);
                            sampleLocationInfo[iProbe]->SetProbeCellID(cellIDRecv);
                            sampleLocationInfo[iProbe]->SetProbeCellNI(cellNIRecv);
                            sampleLocationInfo[iProbe]->SetProbeCellNJ(cellNJRecv);
                            sampleLocationInfo[iProbe]->SetProbeCellNK(cellNKRecv);
                            sampleLocationInfo[iProbe]->SetProbeToCellDistance(distanceRecv);
                        }
                    }
                }
            }
        }
        delete cdata;
        cdata = NULL;
    }
}

void SampleLocationInfo::CheckInvalidProbesData()
{
    int dataMonitorType = GlobalDataBase::GetIntParaFromDB("dataMonitorType");
    int nTotalProbes = GetTotalProbesNumber();
    vector<ProbeData *> probesLocationInfo = GetSampleLocationInfo();

    int count = 0;
    for (int iProbe = 0; iProbe < nTotalProbes; ++ iProbe)
    {
        vector <RDouble> coordinates = probesLocationInfo[iProbe]->GetProbeCoordinates();
        int lineID = probesLocationInfo[iProbe]->GetProbeLineID();
        int surfaceID = probesLocationInfo[iProbe]->GetProbeSurfaceID();
        bool isValid = probesLocationInfo[iProbe]->GetProbeValidity();
        if (!isValid)
        {
            std::ostringstream oss;
            if (dataMonitorType == PROBESMONITOR)
            {
                oss << " Warning: the monitored probe " << "(x, y, z) = " << coordinates[0] << " " << coordinates[1] << " " << coordinates[2] << " " << " is not in the computational domain! "<< endl;
            }
            else if (dataMonitorType == LINESMONITOR)
            {
                oss << " Warning: the monitored probe " << "(x, y, z) = " << coordinates[0] << " " << coordinates[1] << " " << coordinates[2] << " " << " of line "<< lineID <<" is not in the computational domain! "<< endl;
            }
            else if (dataMonitorType == SURFACESMONITOR)
            {
                oss << " Warning: the monitored probe " << "(x, y, z) = " << coordinates[0] << " " << coordinates[1] << " " << coordinates[2] << " " << " of surface "<< surfaceID <<" is not in the computational domain! "<< endl;
            }
            else
            {
                TK_Exit::UnexpectedVarValue("dataMonitorType", dataMonitorType);
            }

            PrintToWindow(oss);
            count ++;
        }
    }

    if (count > 0)
    {
        std::ostringstream oss;
        if (dataMonitorType == PROBESMONITOR)
        {
            oss << " Error: there are " << count << " invalid probes, please modify the probes define file!"<< endl;
        }
        else if (dataMonitorType == LINESMONITOR)
        {
            oss << " Error: there are " << count << " invalid probes, please modify the lines define file!"<< endl;
        }
        else if (dataMonitorType == SURFACESMONITOR)
        {
            oss << " Error: there are " << count << " invalid probes, please modify the surfaces define file!"<< endl;
        }
        else
        {
            TK_Exit::UnexpectedVarValue("dataMonitorType", dataMonitorType);
        }

        TK_Exit::ExceptionExit(oss);
    }
}

void SampleLocationInfo::BroadCastProbesData()
{

    int server = GetServerProcessorID();
    int myid   = GetCurrentProcessorID();

    DataContainer *cdata = new DataContainer();

    if (myid == server)
    {
        CompressProbesData(cdata);
    }

    PH_Bcast(cdata, server);

    if (myid != server)
    {
        DecompressProbesData(cdata);
    }

    delete cdata;
}

void SampleLocationInfo::CompressProbesData(DataContainer *cdata)
{
    cdata->moveto_begin();

    int nTotalProbes = GetTotalProbesNumber();
    PHWrite(cdata, nTotalProbes);

    for (int iProbe = 0; iProbe < nTotalProbes; ++ iProbe)
    {
        int gridType = sampleLocationInfo[iProbe]->GetProbeCellType();

        typeOfTotalProbes.push_back(gridType);
        PHWrite(cdata, typeOfTotalProbes[iProbe]);

        if (gridType == PHSPACE::STRUCTGRID)
        {
            sampleLocationInfo[iProbe]->CompressStructProbeCellInfo(cdata);
        }
        else
        {
            sampleLocationInfo[iProbe]->CompressUnstructProbeCellInfo(cdata);
        }
    }
}

void SampleLocationInfo::DecompressProbesData(DataContainer *cdata)
{
    cdata->moveto_begin();
    sampleLocationInfo.resize(0);

    int nTotalProbesRecv;
    PHRead(cdata, nTotalProbesRecv);

    vector <int> probesTypeRecv;
    probesTypeRecv.resize(nTotalProbesRecv);

    for (int iProbe = 0; iProbe < nTotalProbesRecv; ++ iProbe)
    {
        ProbeData *sampleLocationOnProcessor = new ProbeData();

        PHRead(cdata, probesTypeRecv[iProbe]);

        int gridType = probesTypeRecv[iProbe];

        if (gridType == PHSPACE::STRUCTGRID)
        {
            sampleLocationOnProcessor->DecompressStructProbeCellInfo(cdata);
        }
        else
        {
            sampleLocationOnProcessor->DecompressUnstructProbeCellInfo(cdata);
        }

        sampleLocationInfo.push_back(sampleLocationOnProcessor);
    }

}

void SampleLocationInfo::SampleLocationInfoToGrid()
{
    int nLocalZones = GetNumberofLocalZones();

    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);
        Grid *grid = GetGrid(zoneID, 0);

        int nTotalProbes= GetTotalProbesNumber();
        vector<ProbeData *> probesSampleLocation = GetSampleLocationInfo();

        if (grid->Type() == PHSPACE::UNSTRUCTGRID)
        {
            UnstructGrid * unstructGrid = UnstructGridCast(grid);

            vector <int> probesCellID;
            vector <int> probesGlobalID;
            vector <int> probesLineID;
            vector <int> probesSurfaceID;
            vector <vector<RDouble> > probesCoordinates;
            int count = 0;
            for(int iProbe = 0; iProbe < nTotalProbes; ++ iProbe)
            {
                int cellZoneID = probesSampleLocation[iProbe]->GetProbeCellZoneID();
                bool isValid = probesSampleLocation[iProbe]->GetProbeValidity();
                if (cellZoneID != zoneID || !isValid) continue;
                int cellID = probesSampleLocation[iProbe]->GetProbeCellID();
                vector<RDouble> coordinates = probesSampleLocation[iProbe]->GetProbeCoordinates();
                int probeID = probesSampleLocation[iProbe]->GetProbeGlobalID();
                int lineID = probesSampleLocation[iProbe]->GetProbeLineID();
                int surfaceID = probesSampleLocation[iProbe]->GetProbeSurfaceID();
                probesCellID.push_back(cellID);
                probesGlobalID.push_back(probeID);
                probesLineID.push_back(lineID);
                probesSurfaceID.push_back(surfaceID);
                probesCoordinates.push_back(coordinates);
                count++;
            }

            unstructGrid->SetZoneProbesNumber(count);
            unstructGrid->SetZoneProbesCellID(probesCellID);
            unstructGrid->SetZoneProbesCoordinates(probesCoordinates);
            unstructGrid->SetZoneProbesGlobalID(probesGlobalID);
            unstructGrid->SetZoneProbesLineID(probesLineID);
            unstructGrid->SetZoneProbesSurfaceID(probesSurfaceID);

            probesCellID.clear();
            probesGlobalID.clear();
            probesLineID.clear();
            probesSurfaceID.clear();
            probesCoordinates.clear();
        }
        else
        {
            StructGrid * structGrid = StructGridCast(grid);

            vector <int> probesCellNI;
            vector <int> probesCellNJ;
            vector <int> probesCellNK;
            vector <int> probesGlobalID;
            vector <int> probesLineID;
            vector <int> probesSurfaceID;
            vector <vector<RDouble> > probesCoordinates;

            int count = 0;
            for (int iProbe = 0; iProbe < nTotalProbes; ++ iProbe)
            {
                int cellZoneID = probesSampleLocation[iProbe]->GetProbeCellZoneID();
                bool isValid = probesSampleLocation[iProbe]->GetProbeValidity();
                if (cellZoneID != zoneID || !isValid) continue;
                int cellNI = probesSampleLocation[iProbe]->GetProbeCellNI();
                int cellNJ = probesSampleLocation[iProbe]->GetProbeCellNJ();
                int cellNK = probesSampleLocation[iProbe]->GetProbeCellNK();
                vector<RDouble> coordinates = probesSampleLocation[iProbe]->GetProbeCoordinates();
                int probeID = probesSampleLocation[iProbe]->GetProbeGlobalID();
                int lineID = probesSampleLocation[iProbe]->GetProbeLineID();
                int surfaceID = probesSampleLocation[iProbe]->GetProbeSurfaceID();

                probesCellNI.push_back(cellNI);
                probesCellNJ.push_back(cellNJ);
                probesCellNK.push_back(cellNK);

                probesGlobalID.push_back(probeID);
                probesLineID.push_back(lineID);
                probesSurfaceID.push_back(surfaceID);
                probesCoordinates.push_back(coordinates); 
                count++;
            }
            structGrid->SetZoneProbesNumber(count);
            structGrid->SetZoneProbesCellNI(probesCellNI);
            structGrid->SetZoneProbesCellNJ(probesCellNJ);
            structGrid->SetZoneProbesCellNK(probesCellNK);
            structGrid->SetZoneProbesCoordinates(probesCoordinates);
            structGrid->SetZoneProbesGlobalID(probesGlobalID);
            structGrid->SetZoneProbesLineID(probesLineID);
            structGrid->SetZoneProbesSurfaceID(probesSurfaceID);

            probesCellNI.clear();
            probesCellNJ.clear();
            probesCellNK.clear();
            probesCoordinates.clear();
            probesGlobalID.clear();
            probesLineID.clear();
            probesSurfaceID.clear();
        }
    }
}

void SampleLocationInfo::WriteProbesCellInfo(const string &filename)
{
    string nearestCellInfoFileName = ChangeExtensionOfFileName(filename, "p2c");
    ios_base::openmode openmode = ios_base::out|ios_base::trunc;
    fstream file;
    ParallelOpenFile(file, nearestCellInfoFileName, openmode);

    int nTotalProbes = GetTotalProbesNumber();

    vector<ProbeData *> probesSampleLocation = GetSampleLocationInfo();

    WriteProbesCellInfoHeader(file);

    for (int iProbe = 0; iProbe < nTotalProbes; ++ iProbe)
    {
        int gridType = probesSampleLocation[iProbe]->GetProbeCellType();

        if (gridType == PHSPACE::STRUCTGRID)
        {
            probesSampleLocation[iProbe]->WriteStructProbeCellInfo(file);
        }
        else
        {
            probesSampleLocation[iProbe]->WriteUnstructProbeCellInfo(file);
        }

    }

    ParallelCloseFile(file);

    PrintToWindow("Finally, ", nTotalProbes," probes will be successfully monitored!\n");
}

bool CheckIfProbeInMinMaxBox(RDouble coordX,RDouble coordY,RDouble coordZ,RDouble *pmin,RDouble *pmax)
{
    if (coordX < pmin[0] || coordX > pmax[0]) return false;
    if (coordY < pmin[1] || coordY > pmax[1]) return false;
    if (coordZ < pmin[2] || coordZ > pmax[2]) return false;
    return true;
}

void WriteProbesCellInfoHeader(fstream &file)
{
    int dataMonitorType = GlobalDataBase::GetIntParaFromDB("dataMonitorType");

    ostringstream oss;

    if (IfFileEmpty(file))
    {
        oss << "Title=\"THE CELL INFORMATION OF PROBE\"" << endl;
        oss << "Notice : ProbeCellID means (i, j, k) for structured grid." << endl;
        oss << "Variables=" << endl;
        oss << "ProbeID" << endl;
        if (dataMonitorType == LINESMONITOR)
        {
            oss << "ProbeLineID" << endl;
        }
        if (dataMonitorType == SURFACESMONITOR)
        {
            oss << "ProbeSurfaceID" << endl;
        }

        oss << "ProbeToCellDistance" << endl;
        oss << "ProbeCellID" << endl;
        oss << "ZoneGridType" << endl;
        oss << "ZoneID" << endl;

    }

    WriteASCIIFile(file, oss.str());

}

}
