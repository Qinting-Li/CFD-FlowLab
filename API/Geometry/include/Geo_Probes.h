//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//          PPPPP  H   H  EEEEE  N    N  GGGGG  L      EEEEE  III         +
//          P   P  H   H  E      NN   N  G      L      E       I          +
//          PPPPP  HHHHH  EEEEE  N N  N  G  GG  L      EEEEE   I          +
//          P      H   H  E      N  N N  G   G  L      E       I          +
//          P      H   H  EEEEE  N    N  GGGGG  LLLLL  EEEEE  III         +
//------------------------------------------------------------------------+
//          Platform for Hybrid Engineering Simulation of Flows           +
//          China Aerodynamics Research and Development Center            +
//                     (C) Copyright, Since 2010                          +
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      Geo_Probes.h
//! @brief     To set probes information and data structure for Sampling.
//! @author    MengLiyuan.

#pragma once
#include "Math_BasisFunction.h"
#include "Geo_Grid.h"
#include "Geo_UnstructGrid.h"
#include "Geo_StructGrid.h"
using namespace std;

namespace PHSPACE
{
class ProbeData
{
private:
    string dataName;
    vector <RDouble> coordinates;

    RDouble NearestDistance;
    int NearestCellZoneType;
    int NearestCellZoneID;
    int NearestCellID;
    int NearestNI;
    int NearestNJ;
    int NearestNK;

public:
    ProbeData();
    ~ProbeData();

public:

    //! To set the probe data name.
    void SetProbeDataName(const string &dataName);
    //! To set coordinates by in-parameter.
    void SetCoordinates(const vector <RDouble> &coordinates);
    //! To set distance from the probe to the nearest cell.
    void SetNearestDistance(const RDouble &NearestDistance);
    //! To set zone grid type of the nearest cell.
    void SetNearestCellZoneType(const RDouble &NearestCellZoneType);
    //! To set zone index of the the nearest cell.
    void SetNearestCellZoneID(const int &NearestCellZoneID);
    //! To set index of the the nearest cell(unstructured grid).
    void SetNearestCellID(const int &NearestCellID);

    //! To set the nearest cell node dimension of i direction (structured grid).
    void SetNearestNI(const int &NearestNI);
    //! To set the nearest cell node dimension of j direction (structured grid).
    void SetNearestNJ(const int &NearestNJ);
    //! To set the nearest cell node dimension of k direction (structured grid).
    void SetNearestNK(const int &NearestNK);

    //! To get the probe data name.
    string GetProbeDataName();
    //! To get coordinates of the probe.
    vector <RDouble> GetCoordinates();
    //! To get distance from the probe to the nearest cell.
    RDouble GetNearestDistance();
    //! To get zone grid type of the nearest cell.
    int GetNearestCellZoneType();
    //! To get zone index of the the nearest cell.
    int GetNearestCellZoneID();
    //! To get index of the the nearest cell(unstructured grid).
    int GetNearestCellID();

    //! To get the nearest cell node dimension of i direction (structured grid).
    int GetNearestNI();
    //! To get the nearest cell node dimension of j direction (structured grid).
    int GetNearestNJ();
    //! To get the nearest cell node dimension of k direction (structured grid).
    int GetNearestNK();

    //! To search the cell nearest to the probe.
    void SearchNearestCell(Grid *grid);
    //! To search the cell nearest to the probe for structured grid.
    void SearchNearestCell(StructGrid * grid);
    //! To search the cell nearest to the probe for unstructured grid.
    void SearchNearestCell(UnstructGrid * grid);

    //! To dump the nearest cell information of the each probe for structured grid.
    void WriteStructNearestCellInfo(fstream &file,int &iprobe);
    //! To dump the nearest cell information of the each probe for unstructured grid.
    void WriteUnstructNearestCellInfo(fstream &file,int &iprobe);

    //! Compress each ProbeData nearest cell information to DataContainer for structured grid.
    void CompressStructNearestCellInfo(DataContainer *cdata);
    //! Compress each ProbeData nearest cell information to DataContainer for unstructured grid.
    void CompressUnstructNearestCellInfo(DataContainer *cdata);
   
    //! Decompress DataContainer's nearest cell information into each ProbeData for structured grid.
    void DecompressStructNearestCellInfo(DataContainer *cdata);
    //! Decompress DataContainer's nearest cell information into each ProbeData for unstructured grid.
    void DecompressUnstructNearestCellInfo(DataContainer *cdata);

};

class ProbesLocationInfo
{
private:
    string ProbesFileName;
    int nProbes;
    bool ProbesFileExist;
    vector <int> ProbesType;
    vector <ProbeData *> ProbesLocation;

public:
    ProbesLocationInfo();
    ~ProbesLocationInfo();

public:
    void Run();
    //! To set the file name of probes parameter.
    void SetProbesFileName(const string &ProbesFileName);
    //!To set the grid type of the probes in current zone.
    void SetProbesType(const vector <int> &ProbesType);
    //! To set the number of probes.
    void SetNumberOfProbes(const int &nProbes);
    //! Check if the probes file exist.
    void CheckIfProbesFileExist();

    //! To get the number of probes.
    int GetNumberOfProbes();
    //! To get the file name of probes parameter.
    string GetProbesFileName();

    //!To set the grid type of the probes in current zone.
    vector <int> GetProbesType();

    //! To get all of the probes location.
    vector<ProbeData *> GetProbesLocationInfo();
    //! If the probes file exist or not.
    bool IfProbesFileExist();
    //! No dimension the parameter.
    void Nondimensional();

private:
    //! To read probes parameter in probesDefineFile.
    void ReadProbesData();
    //! Add or update probes' nearest cell information to ProbData.
    //! If there was a data with same name, replace the old with the new.
    void UpdateAllProbesData();
    //! Communicate with each processor,then determine nearest cell of all zones with minimum distance.
    void ServerCollectionProbesData();
    //! Bcast probes data to all processor.
    void BroadCastProbesData();
    //! Compress each ProbeData nearest cell information to DataContainer.
    void CompressProbeData(DataContainer *cdata);
    //! Decompress DataContainer's nearest cell information into each ProbeData
    void DecompressProbeData(DataContainer *cdata);
    //! To save the probes nearest cells data to class grid.
    void ProbesLocationInfoToGrid();
    //! To dump the information of the probes nearest cells.
    void WriteNearestCellInfo(const string &filename);
};

//! Check if the probes in Min-Max box.
bool CheckIfProbeInMinMaxBox(RDouble coordX,RDouble coordY,RDouble coordZ,RDouble *pmin,RDouble *pmax);

//! To dump the nearest cells information file title.
void WriteNearestCellInfoHeader(fstream &file);

}
