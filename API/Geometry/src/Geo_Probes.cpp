#include "Geo_Probes.h"
#include "IO_FileName.h"
#include "GridType.h"
#include "Geo_MultiGridInfo_Struct.h"
#include "TK_Parse.h"
#include "PHIO.h"
#include "TK_Log.h"
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
    NearestDistance =1.0e10;
    NearestCellZoneType = -1;
    NearestCellZoneID = -1;
    NearestCellID = -1;
    NearestNI = -1;
    NearestNJ = -1;
    NearestNK = -1;
    coordinates.resize(0);
}

ProbeData::~ProbeData()
{
    coordinates.clear();
}

void ProbeData::SetProbeDataName(const string &dataName)
{
    this->dataName = dataName;
}

void ProbeData::SetCoordinates(const vector <RDouble> &coordinates)
{
    this->coordinates = coordinates;
}

void ProbeData::SetNearestDistance(const RDouble &NearestDistance)
{
    this->NearestDistance = NearestDistance;
}

void ProbeData::SetNearestCellZoneType(const RDouble &NearestCellZoneType)
{
    this->NearestCellZoneType = NearestCellZoneType;
}

void ProbeData::SetNearestCellZoneID(const int &NearestCellZoneID)
{
    this->NearestCellZoneID = NearestCellZoneID;
}

void ProbeData::SetNearestCellID(const int &NearestCellID)
{
    this->NearestCellID = NearestCellID;
}

void ProbeData::SetNearestNI(const int &NearestNI)
{
    this->NearestNI = NearestNI;
}

void ProbeData::SetNearestNJ(const int &NearestNJ)
{
    this->NearestNJ = NearestNJ;
}

void ProbeData::SetNearestNK(const int &NearestNK)
{
    this->NearestNK = NearestNK;
}

string ProbeData::GetProbeDataName()
{
    return this->dataName;
}

vector <RDouble> ProbeData::GetCoordinates()
{
    return this->coordinates;
}

RDouble ProbeData::GetNearestDistance()
{
    return this->NearestDistance;
}

int ProbeData::GetNearestCellZoneType()
{
    return this->NearestCellZoneType;
}

int ProbeData::GetNearestCellZoneID()
{
    return this->NearestCellZoneID;
}

int ProbeData::GetNearestCellID()
{
    return this->NearestCellID;
}

int ProbeData::GetNearestNI()
{
    return this->NearestNI;
}

int ProbeData::GetNearestNJ()
{
    return this->NearestNJ;
}

int ProbeData::GetNearestNK()
{
    return this->NearestNK;
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

    RDouble3D * xcc = grid->GetCellCenterX();
    RDouble3D * ycc = grid->GetCellCenterY();
    RDouble3D * zcc = grid->GetCellCenterZ();

    RDouble3D * largestLocalGridLength = new RDouble3D(I,J,K,fortranArray);
    PHSPACE::SetField(largestLocalGridLength, -LARGE, ist,ied,jst,jed,kst,ked);

    vector <RDouble> coordinates = GetCoordinates();
    RDouble NearestDistance = GetNearestDistance();
    int NearestCellZoneID = GetNearestCellZoneID();
    int NearestNI = GetNearestNI();
    int NearestNJ = GetNearestNJ();
    int NearestNK = GetNearestNK();

    for ( int k = kst; k <= ked; ++ k )
    {
        for ( int j = jst; j <= jed; ++ j )
        {
            for ( int i = ist; i <= ied; ++ i )
            {
                RDouble CurrentDistance = DISTANCE( (*xcc)(i,j,k)-coordinates[0], (*ycc)(i,j,k)-coordinates[1],(*zcc)(i,j,k)-coordinates[2] );
                if(CurrentDistance < NearestDistance)
                {
                    NearestDistance = CurrentDistance;
                    NearestNI = i;
                    NearestNJ = j;
                    NearestNK = k;
                    NearestCellZoneID = zoneID;
                }
            }
        }
    }

    SetNearestCellZoneType(gridtype);
    SetNearestDistance(NearestDistance);
    SetNearestCellZoneID(NearestCellZoneID);
    SetNearestNI(NearestNI);
    SetNearestNJ(NearestNJ);
    SetNearestNK(NearestNK);
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

    vector <RDouble> coordinates = GetCoordinates();
    RDouble NearestDistance = GetNearestDistance();
    int NearestCellID = GetNearestCellID();
    int NearestCellZoneID = GetNearestCellZoneID();

    for ( int icell = 0; icell < nTotalCell; ++ icell )
    {
        RDouble CurrentDistance = DISTANCE( xcc[icell]-coordinates[0], ycc[icell]-coordinates[1], zcc[icell]-coordinates[2] );
        if(CurrentDistance < NearestDistance)
        {
            NearestDistance = CurrentDistance;
            NearestCellID = icell;
            NearestCellZoneID = zoneID;
        }
    }

    SetNearestCellZoneType(gridtype);
    SetNearestDistance(NearestDistance);
    SetNearestCellZoneID(NearestCellZoneID);
    SetNearestCellID(NearestCellID);

}

void ProbeData::WriteStructNearestCellInfo(fstream &file,int &iprobe)
{
    ostringstream oss;

    int wordwidth = 16;
    oss << setiosflags(ios::left);
    oss << setiosflags(ios::scientific);
    oss << setprecision(4);

    string iProbeName       = GetProbeDataName();
    RDouble NearestDistance = GetNearestDistance();

    int NearestNI = GetNearestNI();
    int NearestNJ = GetNearestNJ();
    int NearestNK = GetNearestNK();

    int NearestCellZoneID   = GetNearestCellZoneID();
    int NearestCellZoneType = GetNearestCellZoneType();


    oss << "	" << iprobe << "	";
    oss << "	" << iProbeName << "	";
    oss << "	" << NearestDistance << "	";
    oss << "	" << " ( " << NearestNI << " , " << NearestNJ << " , " << NearestNK << " ) " << "	";
    oss << "	" << NearestCellZoneType << "	";
    oss << "	" << NearestCellZoneID << "	";
    oss << "\n";

    WriteASCIIFile(file, oss.str());
}

void ProbeData::WriteUnstructNearestCellInfo(fstream &file,int &iprobe)
{
    ostringstream oss;

    int wordwidth = 16;
    oss << setiosflags(ios::left);
    oss << setiosflags(ios::scientific);
    oss << setprecision(4);

    string iProbeName       = GetProbeDataName();
    RDouble NearestDistance = GetNearestDistance();
    int NearestCellID       = GetNearestCellID();
    int NearestCellZoneID   = GetNearestCellZoneID();
    int NearestCellZoneType = GetNearestCellZoneType();


    oss << "	" << iprobe << "	";
    oss << "	" << iProbeName << "	";
    oss << "	" << NearestDistance << "	";
    oss << "	" << NearestCellID << "	";
    oss << "	" << NearestCellZoneType << "	";
    oss << "	" << NearestCellZoneID << "	";
    oss << "\n";

    WriteASCIIFile(file, oss.str());
}

void ProbeData::CompressUnstructNearestCellInfo(DataContainer *cdata)
{
    vector <RDouble> coordinates = GetCoordinates();
    int nearestCellZoneType = GetNearestCellZoneType();
    int nearestCellZoneID = GetNearestCellZoneID();
    int nearestCellID = GetNearestCellID();
    RDouble nearestDistance = GetNearestDistance();

    int nDim =static_cast<int>( coordinates.size());
    PHWrite(cdata, nDim);
    for (int iDim = 0; iDim < nDim; ++ iDim)
    {
        PHWrite(cdata,coordinates[iDim]);
    }

    PHWrite(cdata, nearestCellZoneType);
    PHWrite(cdata, nearestCellZoneID);
    PHWrite(cdata, nearestCellID);
    PHWrite(cdata, nearestDistance);
}

void ProbeData::CompressStructNearestCellInfo(DataContainer *cdata)
{
    vector <RDouble> coordinates = GetCoordinates();
    int nearestCellZoneType = GetNearestCellZoneType();
    int nearestCellZoneID = GetNearestCellZoneID();
    int nearestNI = GetNearestNI();
    int nearestNJ = GetNearestNJ();
    int nearestNK = GetNearestNK();
    RDouble nearestDistance = GetNearestDistance();

    int nDim =static_cast<int>( coordinates.size());
    PHWrite(cdata, nDim);

    for (int iDim = 0; iDim < nDim; ++ iDim)
    {
        PHWrite(cdata,coordinates[iDim]);
    }

    PHWrite(cdata, nearestCellZoneType);
    PHWrite(cdata, nearestCellZoneID);
    PHWrite(cdata, nearestNI);
    PHWrite(cdata, nearestNJ);
    PHWrite(cdata, nearestNK);
    PHWrite(cdata, nearestDistance);

}

void ProbeData::DecompressStructNearestCellInfo(DataContainer *cdata)
{
    RDouble nearestDistanceRecv;
    int nearestCellZoneIDRecv, nearestCellZoneTypeRecv, nearestNIRecv, nearestNJRecv, nearestNKRecv;
    vector <RDouble> coordinatesRecv;
    int nDimRecv;

    PHRead(cdata, nDimRecv);
    coordinatesRecv.resize(nDimRecv);

    for (int iDim = 0; iDim < nDimRecv; ++ iDim)
    {
        PHRead(cdata,coordinatesRecv[iDim]);
    }

    PHRead(cdata, nearestCellZoneTypeRecv);
    PHRead(cdata, nearestCellZoneIDRecv);
    PHRead(cdata, nearestNIRecv);
    PHRead(cdata, nearestNJRecv);
    PHRead(cdata, nearestNKRecv);
    PHRead(cdata, nearestDistanceRecv);

    SetCoordinates(coordinatesRecv);
    SetNearestCellZoneType(nearestCellZoneTypeRecv);
    SetNearestCellZoneID(nearestCellZoneIDRecv);
    SetNearestNI(nearestNIRecv);
    SetNearestNJ(nearestNJRecv);
    SetNearestNK(nearestNKRecv);
    SetNearestDistance(nearestDistanceRecv);
}

void ProbeData::DecompressUnstructNearestCellInfo(DataContainer *cdata)
{
    RDouble nearestDistanceRecv;
    int nearestCellZoneIDRecv, nearestCellZoneTypeRecv, nearestCellIDRecv;
    vector <RDouble> coordinatesRecv;
    int nDimRecv;

    PHRead(cdata, nDimRecv);
    coordinatesRecv.resize(nDimRecv);
    for (int iDim = 0; iDim < nDimRecv; ++ iDim)
    {
        PHRead(cdata,coordinatesRecv[iDim]);
    }

    PHRead(cdata, nearestCellZoneTypeRecv);
    PHRead(cdata, nearestCellZoneIDRecv);
    PHRead(cdata, nearestCellIDRecv);
    PHRead(cdata, nearestDistanceRecv);


    SetCoordinates(coordinatesRecv);
    SetNearestCellZoneType(nearestCellZoneTypeRecv);
    SetNearestCellZoneID(nearestCellZoneIDRecv);
    SetNearestCellID(nearestCellIDRecv);
    SetNearestDistance(nearestDistanceRecv);
}

ProbesLocationInfo::ProbesLocationInfo()
{

    nProbes   = 0;
    ProbesType.resize(0);
    ProbesLocation.resize(0);
    ProbesFileExist = false;
}

ProbesLocationInfo::~ProbesLocationInfo()
{
    for (vector<ProbeData *>::iterator iter = ProbesLocation.begin(); iter != ProbesLocation.end(); ++ iter)
    {
        if (NULL != *iter)
        {
            delete *iter;
            *iter = NULL;
        }
    }
    ProbesLocation.clear();
}

void ProbesLocationInfo::Run()
{
    using namespace PHMPI;
    string probesDefineFile = GlobalDataBase::GetStrParaFromDB("probesDefineFile");
    SetProbesFileName(probesDefineFile);

    CheckIfProbesFileExist();
    if ( IfProbesFileExist())
    {
        ReadProbesData();
    }
    else
    {
        PrintToWindow("Warning : The file used to define probes does not exist! \n");
        return;
    }

    UpdateAllProbesData();

    ServerCollectionProbesData();

    int currentProcessorID = GetCurrentProcessorID();
    if (currentProcessorID == server)
    {
        string gridfile = GlobalDataBase::GetStrParaFromDB( "gridfile" );

        string nearestCellFileName = PHSPACE::AddSymbolToFileName(gridfile, '_', currentProcessorID);

        WriteNearestCellInfo(nearestCellFileName);
    }

    BroadCastProbesData();

    ProbesLocationInfoToGrid();
}

void ProbesLocationInfo::SetNumberOfProbes(const int &nProbes)
{
    this->nProbes = nProbes;
}

int ProbesLocationInfo::GetNumberOfProbes()
{
    return this->nProbes;
}

void ProbesLocationInfo::SetProbesFileName(const string &ProbesFileName)
{
    this->ProbesFileName = ProbesFileName;
}

void ProbesLocationInfo::SetProbesType(const vector <int> &ProbesType)
{
    this->ProbesType = ProbesType;
}

string ProbesLocationInfo::GetProbesFileName()
{
    return this->ProbesFileName;
}

vector <int> ProbesLocationInfo::GetProbesType()
{
    return this->ProbesType;
}

vector <ProbeData*> ProbesLocationInfo::GetProbesLocationInfo()
{
    return this->ProbesLocation;
}

bool ProbesLocationInfo::IfProbesFileExist()
{
    return ProbesFileExist;
}

void ProbesLocationInfo::Nondimensional()
{
    RDouble gridScaleFactor = GlobalDataBase::GetDoubleParaFromDB("gridScaleFactor");
    for (int iProbe = 0; iProbe < GetNumberOfProbes(); ++ iProbe)
    {
        vector <RDouble> coordinates = ProbesLocation[iProbe]->GetCoordinates();

        for (int iDim = 0; iDim < coordinates.size(); ++ iDim)
        {
            coordinates[iDim] *= gridScaleFactor;
        }

        ProbesLocation[iProbe]->SetCoordinates(coordinates);
    }
}

void ProbesLocationInfo::CheckIfProbesFileExist()
{
    using namespace PHMPI;
    ProbesFileExist = false;

    fstream file;
    string keyFile = GetProbesFileName();
    file.open(keyFile.c_str(), ios_base::in);

    if (!file)
    {
        return;
    }

    ProbesFileExist = true;
}

void ProbesLocationInfo::ReadProbesData()
{
    fstream file;
    string keyFile = GetProbesFileName();
    file.open(keyFile.c_str(), ios_base::in);

    string line, word;
    string separator  = " =\t\r\n#$,;\"'";
    string::size_type npos = string::npos;

    string dataName;

    getline(file, line);
    line = FindNextWord(line, word, separator);
    from_string <int> (nProbes, word, std::dec);

    for (int iProbe = 0; iProbe < nProbes; ++ iProbe)
    {
        ProbeData *probe = new ProbeData;

        getline(file, line);
        line = FindNextWord(line, word, separator);
        from_string <string> (dataName, word, std::dec);

        const int nDim = 3;
        vector <RDouble> coordinates;
        getline(file, line);
        for (int iDim = 0; iDim < nDim; ++ iDim)
        {
            line = FindNextWord(line, word, separator);
            RDouble coordinatesIn;
            from_string <RDouble> (coordinatesIn, word, std::dec);
            coordinates.push_back(coordinatesIn);
        }

        probe->SetProbeDataName(dataName);
        probe->SetCoordinates(coordinates);

        ProbesLocation.push_back(probe);
    }

    file.close();
    file.clear();
}

void ProbesLocationInfo::UpdateAllProbesData()
{
    int nLocalZones = GetNumberofLocalZones();
    int nProbes = GetNumberOfProbes();
    vector<ProbeData *> ProbesLocation = GetProbesLocationInfo();

#ifdef USE_OMP
#pragma omp parallel for
#endif
    for ( int izone = 0; izone < nLocalZones; ++ izone )
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(izone);

        Grid *grid = GetGrid(zoneID, 0);

        RDouble *pmin = grid->GetMinBox();
        RDouble *pmax = grid->GetMaxBox();

        for ( int iProbe = 0; iProbe < nProbes; ++ iProbe )
        {
            vector <RDouble> coordinates = ProbesLocation[iProbe]->GetCoordinates();
            bool ProbeInMinMaxBox = CheckIfProbeInMinMaxBox(coordinates[0],coordinates[1],coordinates[2],pmin,pmax);
            if (!ProbeInMinMaxBox) continue;
            ProbesLocation[iProbe]->SearchNearestCell(grid);
        }
    }

}

void ProbesLocationInfo::ServerCollectionProbesData()
{
#ifndef PH_PARALLEL
    return;
#endif

    int server = PHMPI::GetServerProcessorID();
    int myid   = PHMPI::GetCurrentProcessorID();
    int numberOfProcessor = PHMPI::GetNumberOfProcessor();
    DataContainer *cdata = new DataContainer();
    int nProbes = GetNumberOfProbes();

    cdata->moveto_begin();
    RDouble nearestDistance;
    int nearestCellZoneID, nearestCellZoneType, nearestCellID, nearestNI, nearestNJ, nearestNK;
    vector <RDouble> coordinates;
    int nDim;

    for (int iProbe = 0; iProbe < nProbes; ++ iProbe)
    {
        coordinates = ProbesLocation[iProbe]->GetCoordinates();
        nearestCellZoneType = ProbesLocation[iProbe]->GetNearestCellZoneType();
        nearestCellZoneID = ProbesLocation[iProbe]->GetNearestCellZoneID();
        nearestCellID = ProbesLocation[iProbe]->GetNearestCellID();
        nearestNI = ProbesLocation[iProbe]->GetNearestNI();
        nearestNJ = ProbesLocation[iProbe]->GetNearestNJ();
        nearestNK = ProbesLocation[iProbe]->GetNearestNK();
        nearestDistance = ProbesLocation[iProbe]->GetNearestDistance();

        nDim =static_cast<int>( coordinates.size());

        PHWrite(cdata, myid);
        PHWrite(cdata, nDim);
        for (int iDim = 0; iDim < nDim; ++ iDim)
        {
            PHWrite(cdata,coordinates[iDim]);
        }
        PHWrite(cdata, nearestCellZoneType);
        PHWrite(cdata, nearestCellZoneID);
        PHWrite(cdata, nearestCellID);
        PHWrite(cdata, nearestNI);
        PHWrite(cdata, nearestNJ);
        PHWrite(cdata, nearestNK);
        PHWrite(cdata, nearestDistance);
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
            RDouble nearestDistanceRecv;
            int nearestCellZoneIDRecv, nearestCellZoneTypeRecv, nearestCellIDRecv, nearestNIRecv, nearestNJRecv, nearestNKRecv;
            vector <RDouble> coordinatesRecv;
            int nDimRecv;

            for (int iProbe = 0; iProbe < nProbes; ++ iProbe)
            {
                PHRead(cdata, processorOfminDisResRecv);
                PHRead(cdata,nDimRecv);
                coordinatesRecv.resize(nDimRecv);
                for (int iDim = 0; iDim < nDimRecv; ++ iDim)
                {
                    PHRead(cdata,coordinatesRecv[iDim]);
                }
                PHRead(cdata, nearestCellZoneTypeRecv);
                PHRead(cdata, nearestCellZoneIDRecv);
                PHRead(cdata, nearestCellIDRecv);
                PHRead(cdata, nearestNIRecv);
                PHRead(cdata, nearestNJRecv);
                PHRead(cdata, nearestNKRecv);
                PHRead(cdata, nearestDistanceRecv);


                if (nearestDistanceRecv < ProbesLocation[iProbe]->GetNearestDistance())
                {
                    ProbesLocation[iProbe]->SetCoordinates(coordinatesRecv);
                    ProbesLocation[iProbe]->SetNearestCellZoneType(nearestCellZoneTypeRecv);
                    ProbesLocation[iProbe]->SetNearestCellZoneID(nearestCellZoneIDRecv);
                    ProbesLocation[iProbe]->SetNearestCellID(nearestCellIDRecv);
                    ProbesLocation[iProbe]->SetNearestNI(nearestNIRecv);
                    ProbesLocation[iProbe]->SetNearestNJ(nearestNJRecv);
                    ProbesLocation[iProbe]->SetNearestNK(nearestNKRecv);
                    ProbesLocation[iProbe]->SetNearestDistance(nearestDistanceRecv);
                }

            }
            delete cdata;
            cdata = NULL;
        }
    }

}

void ProbesLocationInfo::BroadCastProbesData()
{
#ifndef PH_PARALLEL
    return;
#endif

    using namespace PHMPI;
    int server = GetServerProcessorID();
    int myid   = GetCurrentProcessorID();

    DataContainer *cdata = new DataContainer();

    if (myid == server)
    {
        CompressProbeData(cdata);
    }

    PH_Bcast(cdata, server);

    if (myid != server)
    {
        DecompressProbeData(cdata);
    }

    delete cdata;
}

void ProbesLocationInfo::CompressProbeData(DataContainer *cdata)
{
    cdata->moveto_begin();

    nProbes = GetNumberOfProbes();
    PHWrite(cdata, nProbes);

    for (int iProbe = 0; iProbe < nProbes; ++ iProbe)
    {
        int gridType = ProbesLocation[iProbe]->GetNearestCellZoneType();

        ProbesType.push_back(gridType);
        PHWrite(cdata, ProbesType[iProbe]);

        if(gridType == PHSPACE::STRUCTGRID)
        {
            ProbesLocation[iProbe]->CompressStructNearestCellInfo(cdata);
        }
        else
        {
            ProbesLocation[iProbe]->CompressUnstructNearestCellInfo(cdata);
        }
    }
}

void ProbesLocationInfo::DecompressProbeData(DataContainer *cdata)
{
    cdata->moveto_begin();
    ProbesLocation.resize(0);

    int nProbesRecv;
    PHRead(cdata, nProbesRecv);

    vector <int> ProbesTypeRecv;
    ProbesTypeRecv.resize(nProbesRecv);

    for (int iProbe = 0; iProbe < nProbesRecv; ++ iProbe)
    {
        ProbeData *ProbesLocationOnProcessor = new ProbeData();

        PHRead(cdata, ProbesTypeRecv[iProbe]);

        int gridType = ProbesTypeRecv[iProbe];

        if(gridType == PHSPACE::STRUCTGRID)
        {
            ProbesLocationOnProcessor->DecompressStructNearestCellInfo(cdata);
        }
        else
        {
            ProbesLocationOnProcessor->DecompressUnstructNearestCellInfo(cdata);
        }

        ProbesLocation.push_back(ProbesLocationOnProcessor);
    }

}

void ProbesLocationInfo::ProbesLocationInfoToGrid()
{
    using namespace PHMPI;

    int nLocalZones = GetNumberofLocalZones();

    for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int zoneID = GetLocalZoneIDToGlobalZoneID(iZone);
        Grid *grid = GetGrid(zoneID, 0);

        int nProbes= GetNumberOfProbes();
        vector<ProbeData *> ProbesLocation = GetProbesLocationInfo();

        if(grid->Type() == PHSPACE::UNSTRUCTGRID)
        {
            UnstructGrid * unstructGrid = UnstructGridCast(grid);

            vector <int> probesNearestCellID;
            vector <vector<RDouble> > probesCoordinates;
            int count = 0;
            for(int iProbe = 0; iProbe < nProbes; ++ iProbe)
            {
                int nearestCellZoneID = ProbesLocation[iProbe]->GetNearestCellZoneID();
                if (nearestCellZoneID != zoneID) continue;
                int nearestCellID = ProbesLocation[iProbe]->GetNearestCellID();
                vector<RDouble> coordinates = ProbesLocation[iProbe]->GetCoordinates();

                probesNearestCellID.push_back(nearestCellID);
                probesCoordinates.push_back(coordinates);
                count++;
            }

            unstructGrid->SetZoneProbesNumber(count);
            unstructGrid->SetProbesNearestCellID(probesNearestCellID);
            unstructGrid->SetProbesCoordinates(probesCoordinates);

            probesNearestCellID.clear();
            probesCoordinates.clear();

        }
        else
        {
            StructGrid * structGrid = StructGridCast(grid);

            vector <int> probesNearestNI;
            vector <int> probesNearestNJ;
            vector <int> probesNearestNK;
            vector <vector<RDouble> > probesCoordinates;

            int count = 0;
            for(int iProbe = 0; iProbe < nProbes; ++ iProbe)
            {
                int nearestCellZoneID = ProbesLocation[iProbe]->GetNearestCellZoneID();
                if (nearestCellZoneID != zoneID) continue;
                int nearestNI = ProbesLocation[iProbe]->GetNearestNI();
                int nearestNJ = ProbesLocation[iProbe]->GetNearestNJ();
                int nearestNK = ProbesLocation[iProbe]->GetNearestNK();
                vector<RDouble> coordinates = ProbesLocation[iProbe]->GetCoordinates();

                probesNearestNI.push_back(nearestNI);
                probesNearestNJ.push_back(nearestNJ);
                probesNearestNK.push_back(nearestNK);
                probesCoordinates.push_back(coordinates); 
                count++;
            }
            structGrid->SetZoneProbesNumber(count);
            structGrid->SetProbesNearestNI(probesNearestNI);
            structGrid->SetProbesNearestNJ(probesNearestNJ);
            structGrid->SetProbesNearestNK(probesNearestNK);
            structGrid->SetProbesCoordinates(probesCoordinates);

            probesNearestNI.clear();
            probesNearestNJ.clear();
            probesNearestNK.clear();
            probesCoordinates.clear();
        }

    }
}

void ProbesLocationInfo::WriteNearestCellInfo(const string &filename)
{
    string nearestCellInfoFileName = ChangeExtensionOfFileName(filename, "p2c");
    ios_base::openmode openmode = ios_base::out|ios_base::trunc;
    fstream file;
    ParallelOpenFile(file, nearestCellInfoFileName, openmode);

    int nProbes = GetNumberOfProbes();

    vector<ProbeData *> ProbesLocation = GetProbesLocationInfo();

    WriteNearestCellInfoHeader(file);

    for ( int iProbe = 0; iProbe < nProbes; ++ iProbe )
    {
        int gridType = ProbesLocation[iProbe]->GetNearestCellZoneType();

        if(gridType == PHSPACE::STRUCTGRID)
        {
            ProbesLocation[iProbe]->WriteStructNearestCellInfo(file, iProbe);
        }
        else
        {
            ProbesLocation[iProbe]->WriteUnstructNearestCellInfo(file, iProbe);
        }

    }

    ParallelCloseFile(file);

    PrintToWindow("Finally, ", nProbes," probes will be successfully monitored!\n");
}

bool CheckIfProbeInMinMaxBox(RDouble coordX,RDouble coordY,RDouble coordZ,RDouble *pmin,RDouble *pmax)
{
    if ( coordX < pmin[0] || coordX > pmax[0] ) return false;
    if ( coordY < pmin[1] || coordY > pmax[1] ) return false;
    if ( coordZ < pmin[2] || coordZ > pmax[2] ) return false;
    return true;
}

void WriteNearestCellInfoHeader(fstream &file)
{
    ostringstream oss;

    string note = " Notice : NearestCellID means (i, j, k) for structured grid.";
    string header = "# No		ProbeName		NearestDistance		NearestCellID		ZoneGridType		ZoneID";

    oss << note << endl;
    oss << header << endl;

    WriteASCIIFile(file, oss.str());

}

}
