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
    probeCellZoneType = -1;
    probeCellZoneID = -1;
    probeCellID = -1;
    probeCellNI = -1;
    probeCellNJ = -1;
    probeCellNK = -1;
    surfaceID = -1;
    lineID = -1;
    isValid = false;
    coordinates.resize(0);
}

ProbeData::~ProbeData()
{
    coordinates.clear();
}

void ProbeData::SetCoordinates(const vector <RDouble> &coordinates)
{
    this->coordinates = coordinates;
}

void ProbeData::SetLineID(const int &lineID)
{
    this->lineID = lineID;
}

void ProbeData::SetSurfaceID(const int &surfaceID)
{
    this->surfaceID = surfaceID;
}

void ProbeData::SetProbeValidity(const bool &isValid)
{
    this->isValid = isValid;
}

void ProbeData::SetProbeToCellDistance(const RDouble &probeToCellDistance)
{
    this->probeToCellDistance = probeToCellDistance;
}

void ProbeData::SetProbeCellZoneType(const int &probeCellZoneType)
{
    this->probeCellZoneType = probeCellZoneType;
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

vector <RDouble> ProbeData::GetCoordinates()
{
    return this->coordinates;
}

int ProbeData::GetLineID()
{
    return this->lineID;
}

int ProbeData::GetSurfaceID()
{
    return this->surfaceID;
}

bool ProbeData::GetProbeValidity()
{
    return this->isValid;
}

RDouble ProbeData::GetProbeToCellDistance()
{
    return this->probeToCellDistance;
}

int ProbeData::GetProbeCellZoneType()
{
    return this->probeCellZoneType;
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
    if(grid->Type() == PHSPACE::STRUCTGRID)
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
    using namespace PHMPI;

    int gridtype = grid->Type();
    int zoneID = grid->GetGridID()->GetIndex();

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange( ni, nj, nk, 0, -1, I, J, K );

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    RDouble3D &xcc = *(grid->GetCellCenterX());
    RDouble3D &ycc = *(grid->GetCellCenterY());
    RDouble3D &zcc = *(grid->GetCellCenterZ());

    RDouble maxDistance = grid->CalMaxEdgeLength();

    vector <RDouble> coordinates = GetCoordinates();
    RDouble nearestCellDistance = GetProbeToCellDistance();
    int nearestCellZoneID = GetProbeCellZoneID();
    int nearestCellNI = GetProbeCellNI();
    int nearestCellNJ = GetProbeCellNJ();
    int nearestCellNK = GetProbeCellNK();

    for ( int k = kst; k <= ked; ++ k )
    {
        for ( int j = jst; j <= jed; ++ j )
        {
            for ( int i = ist; i <= ied; ++ i )
            {
                RDouble CurrentDistance = DISTANCE(xcc(i, j, k) - coordinates[0], ycc(i, j, k) - coordinates[1], zcc(i, j, k) - coordinates[2]);
                if(CurrentDistance < nearestCellDistance)
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

    SetProbeCellZoneType(gridtype);
    SetProbeValidity(isValid);
    SetProbeToCellDistance(nearestCellDistance);
    SetProbeCellZoneID(nearestCellZoneID);
    SetProbeCellNI(nearestCellNI);
    SetProbeCellNJ(nearestCellNJ);
    SetProbeCellNK(nearestCellNK);
}

void ProbeData::SearchNearestCell(UnstructGrid * grid)
{
    using namespace PHMPI;

    int gridtype = grid->Type();
    int zoneID = grid->GetGridID()->GetIndex();

    int nTotalCell = grid->GetNTotalCell();
    int nBoundFace = grid->GetNBoundFace();

    RDouble * xcc = grid->GetCellCenterX();
    RDouble * ycc = grid->GetCellCenterY();
    RDouble * zcc = grid->GetCellCenterZ();

    RDouble maxDistance = grid->CalMaxEdgeLength();

    vector <RDouble> coordinates = GetCoordinates();
    bool isValid = GetProbeValidity();
    RDouble nearestCellDistance = GetProbeToCellDistance();
    int nearestCellID = GetProbeCellID();
    int nearestCellZoneID = GetProbeCellZoneID();

    for ( int icell = 0; icell < nTotalCell; ++ icell )
    {
        RDouble CurrentDistance = DISTANCE( xcc[icell]-coordinates[0], ycc[icell]-coordinates[1], zcc[icell]-coordinates[2] );
        if(CurrentDistance < nearestCellDistance)
        {
            nearestCellDistance = CurrentDistance;
            nearestCellID = icell;
            nearestCellZoneID = zoneID;
            if (nearestCellDistance <= maxDistance) isValid = true;
        }
    }

    SetProbeCellZoneType(gridtype);
    SetProbeValidity(isValid);
    SetProbeToCellDistance(nearestCellDistance);
    SetProbeCellZoneID(nearestCellZoneID);
    SetProbeCellID(nearestCellID);

}

void ProbeData::SearchRealCell(Grid *grid, bool &flag)
{
    if(grid->Type() == PHSPACE::STRUCTGRID)
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
    using namespace PHMPI;

    int gridtype = grid->Type();
    int zoneID = grid->GetGridID()->GetIndex();

    int ni = grid->GetNI();
    int nj = grid->GetNJ();
    int nk = grid->GetNK();

    Range I, J, K;
    GetRange( ni, nj, nk, 0, -1, I, J, K );

    int ist, ied, jst, jed, kst, ked;
    grid->GetCellIterationIndex(ist, ied, jst, jed, kst, ked);

    RDouble3D &xcc = *(grid->GetCellCenterX());
    RDouble3D &ycc = *(grid->GetCellCenterY());
    RDouble3D &zcc = *(grid->GetCellCenterZ());

    RDouble4D &xfn = *(grid->GetFaceNormalX());
    RDouble4D &yfn = *(grid->GetFaceNormalY());
    RDouble4D &zfn = *(grid->GetFaceNormalZ());

    vector <RDouble> coordinates = GetCoordinates();
    RDouble realCellDistance = GetProbeToCellDistance();
    bool isValid = GetProbeValidity();
    int realCellZoneID = GetProbeCellZoneID();
    int realCellNI = GetProbeCellNI();
    int realCellNJ = GetProbeCellNJ();
    int realCellNK = GetProbeCellNK();

    int nDim = GetDim();
    int il, jl, kl;
    RDouble small = 1.0e-8;
    for ( int k = kst; k <= ked; ++ k )
    {
        for ( int j = jst; j <= jed; ++ j )
        {
            for ( int i = ist; i <= ied; ++ i )
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

                    if (ABS(dotProduct11 / probeToFaceDistance1) < small)
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
                    if (ABS(dotProduct12 / probeToFaceDistance2) < small)
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

    SetProbeCellZoneType(gridtype);
    SetProbeToCellDistance(realCellDistance);
    SetProbeCellZoneID(realCellZoneID);
    SetProbeCellNI(realCellNI);
    SetProbeCellNJ(realCellNJ);
    SetProbeCellNK(realCellNK);
    SetProbeValidity(isValid);
}

void ProbeData::SearchRealCell(UnstructGrid * grid, bool &flag)
{
    using namespace PHMPI;

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

    vector <RDouble> coordinates = GetCoordinates();
    RDouble realCellDistance = GetProbeToCellDistance();
    bool isValid = GetProbeValidity();
    int realCellID = GetProbeCellID();
    int realCellZoneID = GetProbeCellZoneID();

    RDouble small = 1.0e-8;
    for ( int icell = 0; icell < nTotalCell; ++ icell )
    {
        for ( int iface = 0; iface < faceNumberOfEachCell[icell]; ++ iface )
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
            if (ABS(dotProduct1 / probeToFaceDistance) < small)
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

    SetProbeCellZoneType(gridtype);
    SetProbeToCellDistance(realCellDistance);
    SetProbeCellZoneID(realCellZoneID);
    SetProbeCellID(realCellID);
    SetProbeValidity(isValid);

}

void ProbeData::WriteStructProbeCellInfo(fstream &file,int &iprobe)
{
    int dataMonitorType = GlobalDataBase::GetIntParaFromDB("dataMonitorType");
    ostringstream oss;

    int wordwidth = 16;
    oss << setiosflags(ios::left);
    oss << setiosflags(ios::scientific);
    oss << setprecision(10);

    int ProbeNumber = iprobe + 1;
    vector <RDouble> coordinates = GetCoordinates();
    int ProbeLineID = GetLineID();
    int ProbeSurfaceID = GetSurfaceID();

    RDouble ProbeToCellDistance = GetProbeToCellDistance();
    int ProbeCellNI = GetProbeCellNI();
    int ProbeCellNJ = GetProbeCellNJ();
    int ProbeCellNK = GetProbeCellNK();

    int ProbeCellZoneID   = GetProbeCellZoneID();
    int ProbeCellZoneType = GetProbeCellZoneType(); 

    int wordWidth = 7;

    oss << setw(wordWidth) << ProbeNumber << "	";
    if (dataMonitorType == LINESMONITOR)
    {
        oss << setw(wordWidth) << ProbeLineID << "	";
        oss << " ( " << coordinates[0] << " , " << coordinates[1] << " , " << coordinates[2] << " ) " << "	";
    }
    if (dataMonitorType == SURFACESMONITOR)
    {
        oss << setw(wordWidth) << ProbeSurfaceID << "	";
        oss << " ( " << coordinates[0] << " , " << coordinates[1] << " , " << coordinates[2] << " ) " << "	";
    }
    oss << setw(wordWidth) << ProbeToCellDistance << "	";
    oss << " ( " << ProbeCellNI << " , " << ProbeCellNJ << " , " << ProbeCellNK << " ) " << "	";
    oss << setw(wordWidth) << ProbeCellZoneType << "	";
    oss << setw(wordWidth) << ProbeCellZoneID << "	";
    oss << "\n";

    WriteASCIIFile(file, oss.str());
}

void ProbeData::WriteUnstructProbeCellInfo(fstream &file,int &iprobe)
{
    int dataMonitorType = GlobalDataBase::GetIntParaFromDB("dataMonitorType");
    ostringstream oss;

    int wordwidth = 16;
    oss << setiosflags(ios::left);
    oss << setiosflags(ios::scientific);
    oss << setprecision(10);

    int ProbeNumber = iprobe + 1;
    vector <RDouble> coordinates = GetCoordinates();
    int ProbeLineID = GetLineID();
    int ProbeSurfaceID = GetSurfaceID();

    RDouble ProbeToCellDistance = GetProbeToCellDistance();
    int ProbeCellID       = GetProbeCellID();
    int ProbeCellZoneID   = GetProbeCellZoneID();
    int ProbeCellZoneType = GetProbeCellZoneType();
    
    int wordWidth = 7;
    oss << setw(wordWidth) << ProbeNumber << "	";
    if (dataMonitorType == LINESMONITOR)
    {
        oss << setw(wordWidth) << ProbeLineID << "	";
        oss << " ( " << coordinates[0] << " , " << coordinates[1] << " , " << coordinates[2] << " ) " << "	";
    }
    if (dataMonitorType == SURFACESMONITOR)
    {
        oss << setw(wordWidth) << ProbeSurfaceID << "	";
        oss << " ( " << coordinates[0] << " , " << coordinates[1] << " , " << coordinates[2] << " ) " << "	";
    }
    oss << setw(wordWidth) << ProbeToCellDistance << "	";
    oss << setw(wordWidth) << ProbeCellID << "	";
    oss << setw(wordWidth) << ProbeCellZoneType << "	";
    oss << setw(wordWidth) << ProbeCellZoneID << "	";
    oss << "\n";

    WriteASCIIFile(file, oss.str());
}

void ProbeData::CompressUnstructProbeCellInfo(DataContainer *cdata)
{
    vector <RDouble> coordinates = GetCoordinates();
    int lineID = GetLineID();
    int surfaceID = GetSurfaceID();
    bool isValid = GetProbeValidity();
    int ProbeCellZoneType = GetProbeCellZoneType();
    int ProbeCellZoneID = GetProbeCellZoneID();
    int ProbeCellID = GetProbeCellID();
    RDouble ProbeToCellDistance = GetProbeToCellDistance();

    int nDim =static_cast<int>( coordinates.size());
    PHWrite(cdata, nDim);
    for (int iDim = 0; iDim < nDim; ++ iDim)
    {
        PHWrite(cdata,coordinates[iDim]);
    }

    PHWrite(cdata, lineID);
    PHWrite(cdata, surfaceID);
    PHWrite(cdata, isValid);
    PHWrite(cdata, ProbeCellZoneType);
    PHWrite(cdata, ProbeCellZoneID);
    PHWrite(cdata, ProbeCellID);
    PHWrite(cdata, ProbeToCellDistance);
}

void ProbeData::CompressStructProbeCellInfo(DataContainer *cdata)
{
    vector <RDouble> coordinates = GetCoordinates();
    int lineID = GetLineID();
    int surfaceID = GetSurfaceID();
    bool isValid = GetProbeValidity();
    int ProbeCellZoneType = GetProbeCellZoneType();
    int ProbeCellZoneID = GetProbeCellZoneID();
    int ProbeCellNI = GetProbeCellNI();
    int ProbeCellNJ = GetProbeCellNJ();
    int ProbeCellNK = GetProbeCellNK();
    RDouble ProbeToCellDistance = GetProbeToCellDistance();

    int nDim =static_cast<int>( coordinates.size());
    PHWrite(cdata, nDim);

    for (int iDim = 0; iDim < nDim; ++ iDim)
    {
        PHWrite(cdata,coordinates[iDim]);
    }
    PHWrite(cdata, lineID);
    PHWrite(cdata, surfaceID);
    PHWrite(cdata, isValid);
    PHWrite(cdata, ProbeCellZoneType);
    PHWrite(cdata, ProbeCellZoneID);
    PHWrite(cdata, ProbeCellNI);
    PHWrite(cdata, ProbeCellNJ);
    PHWrite(cdata, ProbeCellNK);
    PHWrite(cdata, ProbeToCellDistance);

}

void ProbeData::DecompressStructProbeCellInfo(DataContainer *cdata)
{
    RDouble ProbeToCellDistanceRecv;
    int ProbeCellZoneIDRecv, ProbeCellZoneTypeRecv, ProbeCellNIRecv, ProbeCellNJRecv, ProbeCellNKRecv;
    vector <RDouble> coordinatesRecv;
    int nDimRecv, lineIDRecv, surfaceIDRecv;
    bool isValidRecv;

    PHRead(cdata, nDimRecv);
    coordinatesRecv.resize(nDimRecv);

    for (int iDim = 0; iDim < nDimRecv; ++ iDim)
    {
        PHRead(cdata,coordinatesRecv[iDim]);
    }

    PHRead(cdata, lineIDRecv);
    PHRead(cdata, surfaceIDRecv);
    PHRead(cdata, isValidRecv);
    PHRead(cdata, ProbeCellZoneTypeRecv);
    PHRead(cdata, ProbeCellZoneIDRecv);
    PHRead(cdata, ProbeCellNIRecv);
    PHRead(cdata, ProbeCellNJRecv);
    PHRead(cdata, ProbeCellNKRecv);
    PHRead(cdata, ProbeToCellDistanceRecv);

    SetCoordinates(coordinatesRecv);
    SetLineID(lineIDRecv);
    SetSurfaceID(surfaceIDRecv);
    SetProbeValidity(isValidRecv);
    SetProbeCellZoneType(ProbeCellZoneTypeRecv);
    SetProbeCellZoneID(ProbeCellZoneIDRecv);
    SetProbeCellNI(ProbeCellNIRecv);
    SetProbeCellNJ(ProbeCellNJRecv);
    SetProbeCellNK(ProbeCellNKRecv);
    SetProbeToCellDistance(ProbeToCellDistanceRecv);
}

void ProbeData::DecompressUnstructProbeCellInfo(DataContainer *cdata)
{
    RDouble ProbeToCellDistanceRecv;
    int ProbeCellZoneIDRecv, ProbeCellZoneTypeRecv, ProbeCellIDRecv;
    vector <RDouble> coordinatesRecv;
    int nDimRecv, lineIDRecv, surfaceIDRecv;
    bool isValidRecv;

    PHRead(cdata, nDimRecv);
    coordinatesRecv.resize(nDimRecv);
    for (int iDim = 0; iDim < nDimRecv; ++ iDim)
    {
        PHRead(cdata,coordinatesRecv[iDim]);
    }

    PHRead(cdata, lineIDRecv);
    PHRead(cdata, surfaceIDRecv);
    PHRead(cdata, isValidRecv);
    PHRead(cdata, ProbeCellZoneTypeRecv);
    PHRead(cdata, ProbeCellZoneIDRecv);
    PHRead(cdata, ProbeCellIDRecv);
    PHRead(cdata, ProbeToCellDistanceRecv);


    SetCoordinates(coordinatesRecv);
    SetLineID(lineIDRecv);
    SetSurfaceID(surfaceIDRecv);
    SetProbeValidity(isValidRecv);
    SetProbeCellZoneType(ProbeCellZoneTypeRecv);
    SetProbeCellZoneID(ProbeCellZoneIDRecv);
    SetProbeCellID(ProbeCellIDRecv);
    SetProbeToCellDistance(ProbeToCellDistanceRecv);
}

SampleLocationInfo::SampleLocationInfo()
{

    numberToMonitor = 0;
    nTotalProbes    = 0;
    nProbes.resize(0);
    ProbesType.resize(0);
    SampleLocation.resize(0);
    defineFileExist = false;
}

SampleLocationInfo::~SampleLocationInfo()
{
    for (vector<ProbeData *>::iterator iter = SampleLocation.begin(); iter != SampleLocation.end(); ++ iter)
    {
        if (NULL != *iter)
        {
            delete *iter;
            *iter = NULL;
        }
    }
    SampleLocation.clear();
}

void SampleLocationInfo::Run()
{
    using namespace PHMPI;
    int dataMonitorType = GlobalDataBase::GetIntParaFromDB("dataMonitorType");

    if (dataMonitorType == PROBESMONITOR)
    {
        string probesDefineFile = GlobalDataBase::GetStrParaFromDB("probesDefineFile");
        SetDefineFileName(probesDefineFile);

        CheckIfDefineFileExist();

        if ( IfDefineFileExist())
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

        if ( IfDefineFileExist())
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
        if ( IfDefineFileExist())
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
        TK_Exit::UnexpectedVarValue( "dataMonitorType", dataMonitorType );
    }

    UpdateAllProbesData();

    ServerCollectionProbesData();

    int currentProcessorID = GetCurrentProcessorID();
    if (currentProcessorID == server)
    {
        CheckInvalidProbesData();

        string gridfile = GlobalDataBase::GetStrParaFromDB( "gridfile" );
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
void SampleLocationInfo::SetTotalProbesNumber(const int &nTotalProbes)
{
    this->nTotalProbes = nTotalProbes;
}

void SampleLocationInfo::SetProbesNumber(const vector <int> &nProbes)
{
    this->nProbes = nProbes;
}

void SampleLocationInfo::SetProbesType(const vector <int> &ProbesType)
{
    this->ProbesType = ProbesType;
}

string SampleLocationInfo::GetDefineFileName()
{
    return this->defineFileName;
}

int SampleLocationInfo::GetTotalProbesNumber()
{
    return this->nTotalProbes;
}

vector <int> SampleLocationInfo::GetProbesNumber()
{
    return this->nProbes;
}

vector <int> SampleLocationInfo::GetProbesType()
{
    return this->ProbesType;
}

vector <ProbeData*> SampleLocationInfo::GetSampleLocationInfo()
{
    return this->SampleLocation;
}

bool SampleLocationInfo::IfDefineFileExist()
{
    return defineFileExist;
}

void SampleLocationInfo::Nondimensional()
{
    RDouble gridScaleFactor = GlobalDataBase::GetDoubleParaFromDB("gridScaleFactor");
    for (int iProbe = 0; iProbe < nTotalProbes; ++ iProbe)
    {
        vector <RDouble> coordinates = SampleLocation[iProbe]->GetCoordinates();

        for (int iDim = 0; iDim < coordinates.size(); ++ iDim)
        {
            coordinates[iDim] *= gridScaleFactor;
        }

        SampleLocation[iProbe]->SetCoordinates(coordinates);
    }
}

void SampleLocationInfo::CheckIfDefineFileExist()
{
    using namespace PHMPI;
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
    string::size_type npos = string::npos;

    while (!file.eof())
    {
        getline(file, line);

        if (line == "") continue;
        FindNextWord(line, word, sep1);
        if (word.substr(0, 1) == "#" || word.substr(0, 2) == "//") continue;
        line = FindNextWord(line, word, separator);

        int nProbesIn;
        from_string <int> (nProbesIn, word, std::dec);
        nTotalProbes = nProbesIn;
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

    if (coordinatesList.size() == nTotalProbes)
    {
        for (int iProbe = 0; iProbe < nTotalProbes; ++ iProbe)
        {
            ProbeData *probe = new ProbeData();

            probe->SetCoordinates(coordinatesList[iProbe]);

            SampleLocation.push_back(probe);
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
    string::size_type npos = string::npos;

    const int nDim = 3;

    vector<vector<RDouble> > coordinatesList;
    vector <int> lineIDList;

    for (int iLine = 1; iLine <= numberToMonitor; ++ iLine)
    {
        getline(file, line);

        line = FindNextWord(line, word, separator);
        int nProbesIn;
        from_string <int> (nProbesIn, word, std::dec);
        nProbes.push_back(nProbesIn);

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
            lineIDList.push_back(iLine);
        }

        nTotalProbes += nProbesIn;
    }

    if (coordinatesList.size() == nTotalProbes)
    {
        for (int iProbe = 0; iProbe < nTotalProbes; ++ iProbe)
        {
            ProbeData *probe = new ProbeData();

            probe->SetCoordinates(coordinatesList[iProbe]);
            probe->SetLineID(lineIDList[iProbe]);

            SampleLocation.push_back(probe);
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
    string::size_type npos = string::npos;

    const int nDim = 3;

    vector<vector<RDouble> > coordinatesList;
    vector <int> surfaceIDList;

    for (int iSurface = 1; iSurface <= numberToMonitor; ++ iSurface)
    {
        getline(file, line);
        line = FindNextWord(line, word, separator);
        int nProbesIn_i;
        from_string <int> (nProbesIn_i, word, std::dec);

        int nProbesIn_j;
        line = FindNextWord(line, word, separator);
        from_string <int> (nProbesIn_j, word, std::dec);

        int nProbesIn = nProbesIn_i *nProbesIn_j;

        nProbes.push_back(nProbesIn);

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
            surfaceIDList.push_back(iSurface);
        }

        nTotalProbes += nProbesIn;
    }

    if (coordinatesList.size() == nTotalProbes)
    {
        for (int iProbe = 0; iProbe < nTotalProbes; ++ iProbe)
        {
            ProbeData *probe = new ProbeData();

            probe->SetCoordinates(coordinatesList[iProbe]);
            probe->SetSurfaceID(surfaceIDList[iProbe]);

            SampleLocation.push_back(probe);
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
    using namespace PHMPI;
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
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();
    int nTotalProbes = GetTotalProbesNumber();
    vector<ProbeData *> SampleLocation = GetSampleLocationInfo();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for ( int iProbe = 0; iProbe < nTotalProbes; ++ iProbe )
    {
        vector <RDouble> coordinates = SampleLocation[iProbe]->GetCoordinates();
        for ( int izone = 0; izone < nLocalZones; ++ izone )
        {
            int zoneID = GetLocalZoneIDToGlobalZoneID(izone);

            Grid *grid = GetGrid(zoneID, 0);

            RDouble *pmin = grid->GetMinBox();
            RDouble *pmax = grid->GetMaxBox();

            bool ProbeInMinMaxBox = CheckIfProbeInMinMaxBox(coordinates[0],coordinates[1],coordinates[2],pmin,pmax);
            if (!ProbeInMinMaxBox) continue;
            SampleLocation[iProbe]->SearchNearestCell(grid);
        }
    }
}

void SampleLocationInfo::CheckInvalidProbesData()
{
    int dataMonitorType = GlobalDataBase::GetIntParaFromDB("dataMonitorType");
    int nTotalProbes = GetTotalProbesNumber();
    vector<ProbeData *> SampleLocation = GetSampleLocationInfo();

    int count = 0;
    for ( int iProbe = 0; iProbe < nTotalProbes; ++ iProbe )
    {
        vector <RDouble> coordinates = SampleLocation[iProbe]->GetCoordinates();
        RDouble probeToCellDistance = SampleLocation[iProbe]->GetProbeToCellDistance();
        int probeLineID = SampleLocation[iProbe]->GetLineID();
        int probeSurfaceID = SampleLocation[iProbe]->GetSurfaceID();
        bool isProbeValid = SampleLocation[iProbe]->GetProbeValidity();
        if (!isProbeValid)
        {
            std::ostringstream oss;
            if (dataMonitorType == PROBESMONITOR)
            {
                oss << " Warning: the monitored probe " << "(x, y, z) = " << coordinates[0] << " " << coordinates[1] << " " << coordinates[2] << " " << " is not in the computational domain! "<< endl;
            }
            else if (dataMonitorType == LINESMONITOR)
            {
                oss << " Warning: the monitored probe " << "(x, y, z) = " << coordinates[0] << " " << coordinates[1] << " " << coordinates[2] << " " << " of line "<< probeLineID <<" is not in the computational domain! "<< endl;
            }
            else if (dataMonitorType == SURFACESMONITOR)
            {
                oss << " Warning: the monitored probe " << "(x, y, z) = " << coordinates[0] << " " << coordinates[1] << " " << coordinates[2] << " " << " of surface "<< probeSurfaceID <<" is not in the computational domain! "<< endl;
            }
            else
            {
                TK_Exit::UnexpectedVarValue( "dataMonitorType", dataMonitorType );
            }

            //WriteLogFile(oss);
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
            TK_Exit::UnexpectedVarValue( "dataMonitorType", dataMonitorType );
        }

        TK_Exit::ExceptionExit(oss);
    }
}

void SampleLocationInfo::UpdateRealCellsData()
{
    using namespace PHMPI;
    int nLocalZones = GetNumberofLocalZones();
    int nTotalProbes = GetTotalProbesNumber();
    vector<ProbeData *> SampleLocation = GetSampleLocationInfo();

#ifdef USE_OMP
#pragma omp parallel for
#endif

    for ( int iProbe = 0; iProbe < nTotalProbes; ++ iProbe )
    {
        bool flag = false;
        vector <RDouble> coordinates = SampleLocation[iProbe]->GetCoordinates();

        for ( int izone = 0; izone < nLocalZones; ++ izone )
        {
            int zoneID = GetLocalZoneIDToGlobalZoneID(izone);

            Grid *grid = GetGrid(zoneID, 0);

            RDouble *pmin = grid->GetMinBox();
            RDouble *pmax = grid->GetMaxBox();

            bool ProbeInMinMaxBox = CheckIfProbeInMinMaxBox(coordinates[0],coordinates[1],coordinates[2],pmin,pmax);
            if (!ProbeInMinMaxBox) continue;
            SampleLocation[iProbe]->SearchRealCell(grid, flag);
            if (flag) break;
        }

    }
}

void SampleLocationInfo::ServerCollectionProbesData()
{
#ifndef PH_PARALLEL
    return;
#endif

    int server = PHMPI::GetServerProcessorID();
    int myid   = PHMPI::GetCurrentProcessorID();
    int numberOfProcessor = PHMPI::GetNumberOfProcessor();
    DataContainer *cdata = new DataContainer();
    int nTotalProbes = GetTotalProbesNumber();

    cdata->moveto_begin();
    RDouble ProbeToCellDistance;
    int ProbeCellZoneID, ProbeCellZoneType, ProbeCellID, ProbeCellNI, ProbeCellNJ, ProbeCellNK;
    vector <RDouble> coordinates;
    int nDim, lineID, surfaceID;
    bool isValid;

    for (int iProbe = 0; iProbe < nTotalProbes; ++ iProbe)
    {
        coordinates = SampleLocation[iProbe]->GetCoordinates();
        lineID = SampleLocation[iProbe]->GetLineID();
        surfaceID = SampleLocation[iProbe]->GetSurfaceID();
        isValid = SampleLocation[iProbe]->GetProbeValidity();
        ProbeCellZoneType = SampleLocation[iProbe]->GetProbeCellZoneType();
        ProbeCellZoneID = SampleLocation[iProbe]->GetProbeCellZoneID();
        ProbeCellID = SampleLocation[iProbe]->GetProbeCellID();
        ProbeCellNI = SampleLocation[iProbe]->GetProbeCellNI();
        ProbeCellNJ = SampleLocation[iProbe]->GetProbeCellNJ();
        ProbeCellNK = SampleLocation[iProbe]->GetProbeCellNK();
        ProbeToCellDistance = SampleLocation[iProbe]->GetProbeToCellDistance();

        nDim =static_cast<int>( coordinates.size());

        PHWrite(cdata, myid);
        PHWrite(cdata, nDim);
        for (int iDim = 0; iDim < nDim; ++ iDim)
        {
            PHWrite(cdata,coordinates[iDim]);
        }
        PHWrite(cdata, lineID);
        PHWrite(cdata, surfaceID);
        PHWrite(cdata, isValid);
        PHWrite(cdata, ProbeCellZoneType);
        PHWrite(cdata, ProbeCellZoneID);
        PHWrite(cdata, ProbeCellID);
        PHWrite(cdata, ProbeCellNI);
        PHWrite(cdata, ProbeCellNJ);
        PHWrite(cdata, ProbeCellNK);
        PHWrite(cdata, ProbeToCellDistance);
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

            cdata->moveto_begin();

            int processorOfminDisResRecv = 0;
            RDouble ProbeToCellDistanceRecv;
            int ProbeCellZoneIDRecv, ProbeCellZoneTypeRecv, ProbeCellIDRecv, ProbeCellNIRecv, ProbeCellNJRecv, ProbeCellNKRecv;
            vector <RDouble> coordinatesRecv;
            int nDimRecv,lineIDRecv,surfaceIDRecv;
            bool isValidRecv;
            for (int iProbe = 0; iProbe < nTotalProbes; ++ iProbe)
            {
                PHRead(cdata, processorOfminDisResRecv);
                PHRead(cdata,nDimRecv);
                coordinatesRecv.resize(nDimRecv);
                for (int iDim = 0; iDim < nDimRecv; ++ iDim)
                {
                    PHRead(cdata,coordinatesRecv[iDim]);
                }
                PHRead(cdata, lineIDRecv);
                PHRead(cdata, surfaceIDRecv);
                PHRead(cdata, isValidRecv);
                PHRead(cdata, ProbeCellZoneTypeRecv);
                PHRead(cdata, ProbeCellZoneIDRecv);
                PHRead(cdata, ProbeCellIDRecv);
                PHRead(cdata, ProbeCellNIRecv);
                PHRead(cdata, ProbeCellNJRecv);
                PHRead(cdata, ProbeCellNKRecv);
                PHRead(cdata, ProbeToCellDistanceRecv);


                if (ProbeToCellDistanceRecv < SampleLocation[iProbe]->GetProbeToCellDistance())
                {
                    SampleLocation[iProbe]->SetCoordinates(coordinatesRecv);
                    SampleLocation[iProbe]->SetLineID(lineIDRecv);
                    SampleLocation[iProbe]->SetSurfaceID(surfaceIDRecv);
                    SampleLocation[iProbe]->SetProbeValidity(isValidRecv);
                    SampleLocation[iProbe]->SetProbeCellZoneType(ProbeCellZoneTypeRecv);
                    SampleLocation[iProbe]->SetProbeCellZoneID(ProbeCellZoneIDRecv);
                    SampleLocation[iProbe]->SetProbeCellID(ProbeCellIDRecv);
                    SampleLocation[iProbe]->SetProbeCellNI(ProbeCellNIRecv);
                    SampleLocation[iProbe]->SetProbeCellNJ(ProbeCellNJRecv);
                    SampleLocation[iProbe]->SetProbeCellNK(ProbeCellNKRecv);
                    SampleLocation[iProbe]->SetProbeToCellDistance(ProbeToCellDistanceRecv);
                }

            }
            delete cdata;
            cdata = NULL;
        }
    }
}

void SampleLocationInfo::BroadCastProbesData()
{
    int searchCellsMethod = GlobalDataBase::GetIntParaFromDB("searchCellsMethod");
    if (searchCellsMethod == 1) return;

#ifndef PH_PARALLEL
    return;
#endif

    using namespace PHMPI;
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

    nTotalProbes = GetTotalProbesNumber();
    PHWrite(cdata, nTotalProbes);

    for (int iProbe = 0; iProbe < nTotalProbes; ++ iProbe)
    {
        int gridType = SampleLocation[iProbe]->GetProbeCellZoneType();

        ProbesType.push_back(gridType);
        PHWrite(cdata, ProbesType[iProbe]);

        if(gridType == PHSPACE::STRUCTGRID)
        {
            SampleLocation[iProbe]->CompressStructProbeCellInfo(cdata);
        }
        else
        {
            SampleLocation[iProbe]->CompressUnstructProbeCellInfo(cdata);
        }
    }
}

void SampleLocationInfo::DecompressProbesData(DataContainer *cdata)
{
    cdata->moveto_begin();
    SampleLocation.resize(0);

    int nProbesRecv;
    PHRead(cdata, nProbesRecv);

    vector <int> ProbesTypeRecv;
    ProbesTypeRecv.resize(nProbesRecv);

    for (int iProbe = 0; iProbe < nProbesRecv; ++ iProbe)
    {
        ProbeData *SampleLocationOnProcessor = new ProbeData();

        PHRead(cdata, ProbesTypeRecv[iProbe]);

        int gridType = ProbesTypeRecv[iProbe];

        if(gridType == PHSPACE::STRUCTGRID)
        {
            SampleLocationOnProcessor->DecompressStructProbeCellInfo(cdata);
        }
        else
        {
            SampleLocationOnProcessor->DecompressUnstructProbeCellInfo(cdata);
        }

        SampleLocation.push_back(SampleLocationOnProcessor);
    }

}

void SampleLocationInfo::SampleLocationInfoToGrid()
{
    using namespace PHMPI;

    int nLocalZones = GetNumberofLocalZones();

    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);
        Grid *grid = GetGrid(zoneID, 0);

        int nTotalProbes= GetTotalProbesNumber();
        vector<ProbeData *> SampleLocation = GetSampleLocationInfo();

        if(grid->Type() == PHSPACE::UNSTRUCTGRID)
        {
            UnstructGrid * unstructGrid = UnstructGridCast(grid);

            vector <int> ProbesCellID;
            vector <int> ProbesLineID;
            vector <int> ProbesSurfaceID;
            vector <vector<RDouble> > ProbesCoordinates;
            int count = 0;
            for(int iProbe = 0; iProbe < nTotalProbes; ++ iProbe)
            {
                int ProbeCellZoneID = SampleLocation[iProbe]->GetProbeCellZoneID();
                bool isProbeValid = SampleLocation[iProbe]->GetProbeValidity();
                if (ProbeCellZoneID != zoneID || !isProbeValid) continue;
                int ProbeCellID = SampleLocation[iProbe]->GetProbeCellID();
                vector<RDouble> coordinates = SampleLocation[iProbe]->GetCoordinates();
                int lineID = SampleLocation[iProbe]->GetLineID();
                int surfaceID = SampleLocation[iProbe]->GetSurfaceID();
                ProbesCellID.push_back(ProbeCellID);
                ProbesLineID.push_back(lineID);
                ProbesSurfaceID.push_back(surfaceID);
                ProbesCoordinates.push_back(coordinates);
                count++;
            }

            unstructGrid->SetZoneProbesNumber(count);
            unstructGrid->SetProbesCellID(ProbesCellID);
            unstructGrid->SetProbesCoordinates(ProbesCoordinates);
            unstructGrid->SetProbesLineID(ProbesLineID);
            unstructGrid->SetProbesSurfaceID(ProbesSurfaceID);

            ProbesCellID.clear();
            ProbesLineID.clear();
            ProbesSurfaceID.clear();
            ProbesCoordinates.clear();
        }
        else
        {
            StructGrid * structGrid = StructGridCast(grid);

            vector <int> ProbesCellNI;
            vector <int> ProbesCellNJ;
            vector <int> ProbesCellNK;
            vector <int> ProbesLineID;
            vector <int> ProbesSurfaceID;
            vector <vector<RDouble> > ProbesCoordinates;

            int count = 0;
            for(int iProbe = 0; iProbe < nTotalProbes; ++ iProbe)
            {
                int ProbeCellZoneID = SampleLocation[iProbe]->GetProbeCellZoneID();
                bool isProbeValid = SampleLocation[iProbe]->GetProbeValidity();
                if (ProbeCellZoneID != zoneID || !isProbeValid) continue;
                int ProbeCellNI = SampleLocation[iProbe]->GetProbeCellNI();
                int ProbeCellNJ = SampleLocation[iProbe]->GetProbeCellNJ();
                int ProbeCellNK = SampleLocation[iProbe]->GetProbeCellNK();
                vector<RDouble> coordinates = SampleLocation[iProbe]->GetCoordinates();
                int lineID = SampleLocation[iProbe]->GetLineID();
                int surfaceID = SampleLocation[iProbe]->GetSurfaceID();

                ProbesCellNI.push_back(ProbeCellNI);
                ProbesCellNJ.push_back(ProbeCellNJ);
                ProbesCellNK.push_back(ProbeCellNK);

                ProbesLineID.push_back(lineID);
                ProbesSurfaceID.push_back(surfaceID);
                ProbesCoordinates.push_back(coordinates); 
                count++;
            }
            structGrid->SetZoneProbesNumber(count);
            structGrid->SetProbesCellNI(ProbesCellNI);
            structGrid->SetProbesCellNJ(ProbesCellNJ);
            structGrid->SetProbesCellNK(ProbesCellNK);
            structGrid->SetProbesCoordinates(ProbesCoordinates);
            structGrid->SetProbesLineID(ProbesLineID);
            structGrid->SetProbesSurfaceID(ProbesSurfaceID);

            ProbesCellNI.clear();
            ProbesCellNJ.clear();
            ProbesCellNK.clear();
            ProbesCoordinates.clear();
            ProbesLineID.clear();
            ProbesSurfaceID.clear();
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

    vector<ProbeData *> SampleLocation = GetSampleLocationInfo();

    WriteProbesCellInfoHeader(file);

    for ( int iProbe = 0; iProbe < nTotalProbes; ++ iProbe )
    {
        int gridType = SampleLocation[iProbe]->GetProbeCellZoneType();

        if(gridType == PHSPACE::STRUCTGRID)
        {
            SampleLocation[iProbe]->WriteStructProbeCellInfo(file, iProbe);
        }
        else
        {
            SampleLocation[iProbe]->WriteUnstructProbeCellInfo(file, iProbe);
        }

    }

    ParallelCloseFile(file);

    PrintToWindow("Finally, ", nTotalProbes," probes will be successfully monitored!\n");
}

bool CheckIfProbeInMinMaxBox(RDouble coordX,RDouble coordY,RDouble coordZ,RDouble *pmin,RDouble *pmax)
{
    if ( coordX < pmin[0] || coordX > pmax[0] ) return false;
    if ( coordY < pmin[1] || coordY > pmax[1] ) return false;
    if ( coordZ < pmin[2] || coordZ > pmax[2] ) return false;
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
        oss << "Number" << endl;
        if (dataMonitorType == LINESMONITOR)
        {
            oss << "ProbeLineID" << endl;
            oss << "ProbeCoordinates" << endl;
        }
        if (dataMonitorType == SURFACESMONITOR)
        {
            oss << "ProbeSurfaceID" << endl;
            oss << "ProbeCoordinates" << endl;
        }
        oss << "ProbeToCellDistance" << endl;
        oss << "ProbeCellID" << endl;
        oss << "ZoneGridType" << endl;
        oss << "ZoneID" << endl;


    }

    WriteASCIIFile(file, oss.str());

}

}
