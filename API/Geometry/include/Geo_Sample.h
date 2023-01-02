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
//! @file      Geo_Sample.h
//! @brief     To set location information and data structure for Sampling.
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
    vector <RDouble> probeCoordinates;
    int probeGlobalID;
    int probeLineID;
    int probeSurfaceID;
    bool isProbeValid;
    RDouble probeToCellDistance;
    int probeCellType;
    int probeCellZoneID;
    int probeCellID;
    int probeCellNI;
    int probeCellNJ;
    int probeCellNK;

public:
    ProbeData();
    ~ProbeData();

public:
    //! To set coordinates by in-parameter.
    void SetProbeCoordinates(const vector <RDouble> &probeCoordinates);
    //! To set global index of the probe.
    void SetProbeGlobalID(const int &probeGlobalID);
    //! To set line index of the probe.
    void SetProbeLineID(const int &probeLineID);
    //! To set surface index of the probe.
    void SetProbeSurfaceID(const int &probeSurfaceID);
    //! To set validity of the probe.
    void SetProbeValidity(const bool &isProbeValid);

    //! To set distance from the probe to its cell.
    void SetProbeToCellDistance(const RDouble &probeToCellDistance);
    //! To set grid type of the probe cell.
    void SetProbeCellType(const int &probeCellType);
    //! To set zone index of the probe cell.
    void SetProbeCellZoneID(const int &probeCellZoneID);
    //! To set index of the probe cell(unstructured grid).
    void SetProbeCellID(const int &probeCellID);

    //! To set the probe cell dimension of i direction (structured grid).
    void SetProbeCellNI(const int &probeCellNI);
    //! To set the probe cell dimension of j direction (structured grid).
    void SetProbeCellNJ(const int &probeCellNJ);
    //! To set the probe cell dimension of k direction (structured grid).
    void SetProbeCellNK(const int &probeCellNK);

    //! To get coordinates of the probe.
    vector <RDouble> GetProbeCoordinates();
    //! To get validity of the probe.
    bool GetProbeValidity();

    //! To get distance from the probe to its cell.
    RDouble GetProbeToCellDistance();
    //! To get grid type of the probe cell.
    int GetProbeCellType();
    //! To get zone index of the probe cell.
    int GetProbeCellZoneID();
    //! To get index of the probe cell(unstructured grid).
    int GetProbeCellID();
    //! To Get global index of the probe.
    int GetProbeGlobalID();
    //! To Get line index of the probe.
    int GetProbeLineID();
    //! To Get surface index of the probe.
    int GetProbeSurfaceID();
    //! To get the probe cell dimension of i direction (structured grid).
    int GetProbeCellNI();
    //! To get the probe cell node dimension of j direction (structured grid).
    int GetProbeCellNJ();
    //! To get the probe cell node dimension of k direction (structured grid).
    int GetProbeCellNK();

    //! To search the cell nearest to the probe.
    void SearchNearestCell(Grid *grid);
    //! To search the cell nearest to the probe for structured grid.
    void SearchNearestCell(StructGrid *grid);
    //! To search the cell nearest to the probe for unstructured grid.
    void SearchNearestCell(UnstructGrid *grid);

    //! To search the real cell where the probe is located.
    void SearchRealCell(Grid *grid, bool &flag);
    //! To search the real cell where the probe is located for structured grid.
    void SearchRealCell(StructGrid *grid, bool &flag);
    //! To search the real cell where the probe is located for unstructured grid.
    void SearchRealCell(UnstructGrid *grid, bool &flag);

    //! To dump the probe cell information of the each probe for structured grid.
    void WriteStructProbeCellInfo(fstream &file);
    //! To dump the probe cell information of the each probe for unstructured grid.
    void WriteUnstructProbeCellInfo(fstream &file);

    //! Compress each ProbeData probe cell information to DataContainer for structured grid.
    void CompressStructProbeCellInfo(DataContainer *cdata);
    //! Compress each ProbeData probe cell information to DataContainer for unstructured grid.
    void CompressUnstructProbeCellInfo(DataContainer *cdata);
   
    //! Decompress DataContainer's probe cell information into each ProbeData for structured grid.
    void DecompressStructProbeCellInfo(DataContainer *cdata);
    //! Decompress DataContainer's probe cell information into each ProbeData for unstructured grid.
    void DecompressUnstructProbeCellInfo(DataContainer *cdata);
};

class SampleLocationInfo
{
private:

    string defineFileName;
    int numberToMonitor;
    bool defineFileExist;
    int numberOfTotalProbes;
    vector <int> numberOfProbes;
    vector <int> surfaceProbes_ni;
    vector <int> surfaceProbes_nj;
    vector <int> typeOfTotalProbes;
    vector <ProbeData *> sampleLocationInfo;

public:
    SampleLocationInfo();
    ~SampleLocationInfo();

public:
    void Run();
    //! To set the file name of data monitor parameter.
    void SetDefineFileName(const string &defineFileName);
    //! To set the number of lines or surfaces.
    void SetNumberToMonitor(const int &numberToMonitor);
    //! To set the total number of probes.
    void SetTotalProbesNumber(const int &numberOfTotalProbes);
    //! To set the probes number of each line or surface.
    void SetProbesNumber(const vector <int> &numberOfProbes);
    //! To set the probe dimension of I direction of each surface.
    void SetSurfaceProbesNI(const vector <int> &surfaceProbes_ni);
    //! To set the probe dimension of J direction of each surface..
    void SetSurfaceProbesNJ(const vector <int> &surfaceProbes_nj);
    //!To set the grid type of all probes.
    void SetTotalProbesType(const vector <int> &typeOfTotalProbes);

    //! Check if the defined file exist(xxx_XYZ.dat).
    void CheckIfDefineFileExist();
    //! To get the number of lines or surfaces.
    int GetNumberToMonitor();

    //! To get the file name of data monitor parameter.
    string GetDefineFileName();
    //!To get the total number of probes.
    int GetTotalProbesNumber();
    //!To get the probes number of each line or surface.
    vector <int> GetProbesNumber();
    //! To get the probe dimension of I direction of each surface.
    vector <int> GetSurfaceProbesNI();
    //! To get the probe dimension of J direction of each surface..
    vector <int> GetSurfaceProbesNJ();
    //!To get the grid type of all probes.
    vector <int> GetTotalProbesType();

    //! To get the location information of all probes.
    vector<ProbeData *> GetSampleLocationInfo();
    //! If the defined file exist or not(xxx_XYZ.dat).
    bool IfDefineFileExist();
    //! No dimension the parameter.
    void Nondimensional();

private:
    //! To read coordination parameter in probesDefineFile.
    void ReadProbesData();
    //! To read coordination parameter in linesDefineFile.
    void ReadLinesData();
    //! To read coordination parameter in surfacesDefineFile.
    void ReadSurfacesData();
    //! Check the invalid probes data in define file.
    void CheckInvalidProbesData();
    //! Add or update all the probes' cell information to ProbData.
    void UpdateAllProbesData();
    //! Add or update all the probes' nearest cell information to ProbData.
    void UpdateNearestCellsData();
    //! Add or update all the probes' real cell information to ProbData.
    void UpdateRealCellsData();
    //! Add or update some global parameters about probes information.
    void UpdateProbesGlobalPara();
    //! Communicate with each processor,then determine probes' cell of all zones with minimum distance.
    void ServerCollectionProbesData();
    //! Broadcast all probes data to all processor.
    void BroadcastProbesData();
    //! Compress each ProbeData probe cell information to DataContainer.
    void CompressProbesData(DataContainer *cdata);
    //! Decompress DataContainer's probe cell information into each ProbeData
    void DecompressProbesData(DataContainer *cdata);
    //! To save the all probes cells data to class grid.
    void SampleLocationInfoToGrid();
    //! To dump the information of the all probes cells.
    void WriteProbesCellInfo(const string &filename);
};

//! Check if the probe in Min-Max box.
bool CheckIfProbeInMinMaxBox(RDouble coordX, RDouble coordY, RDouble coordZ, RDouble *pmin, RDouble *pmax);

//! To dump the probes cells information file title.
void WriteProbesCellInfoHeader(fstream &file);

}
