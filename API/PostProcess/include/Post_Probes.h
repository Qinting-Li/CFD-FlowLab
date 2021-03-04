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
//! @file      Post_Probes.h
//! @brief     Post_Probes defines the post variables of probes that need to be monitored.
//! @author    Meng Liyuan.

#pragma once
#include "GridType.h"
#include "Geo_UnstructGrid.h"
#include "Geo_StructGrid.h"

using namespace std;

namespace PHSPACE
{

class Grid;
class ProbesLocationInfo;
//! @brief Post_Probes class is used to manage post variables data of probes.
class Post_Probes
{
public:
    //! @param[in] nVariables_in    Number of probe variables, parsed from parameter file usually.
    //! @param[in] Variables_in     Type of probe variables, parsed from parameter file usually.
    LIB_EXPORT Post_Probes(int nProbeVariables_in, int *probeVariables_in, Grid *grid_in = NULL);

    LIB_EXPORT ~Post_Probes();

public:
    //! Return the ASCII type file title header.
    LIB_EXPORT void GetASCIIFileHeader(vector<string> &title_tecplot);

    //! Return number of probe variables.
    int GetNumberofProbeVariables();

    //! Return if the variableType should be Monitored.
    bool IsNeedMonitoring(int variableType);

    //! Compute variables data of probe through variables on nearest cell center.
    //! Unstructured grid data type.
    //! @param[in] name              Variable data's name.
    //! @param[in] cellCenterData    Pointer the nearest cell center data.
    LIB_EXPORT void UpdateProbeVarPtr(const string &name, RDouble *NearestCellCenterData);

    //! Compute variables data of probe through variables on nearest cell center.
    //! Structured grid data type.
    //! @param[in] name              Variable data's name.
    //! @param[in] cellCenterData    Pointer the nearest cell center data.
    //! @param[in] index             The array index from RDouble4D array,
    //!                              i.e we only need the cellCenterData[...][...][...][index].
    LIB_EXPORT void UpdateProbeVarPtr(const string &name, RDouble4D &NearestCellCenterData, int index);

    //! Obtain probe variable data's pointer(address) from database by name.
    //! If the dataSet doesn't has a same name data, returns null.
    //! @param[in] name    Variable data's name.
    void * GetProbeVarPtr(const string &name);

    //! Get all probe variables' pointer.
    LIB_EXPORT void GetAllProbeVarPtr(RDouble **q_probes);

    //! Get probe variables' name through variable's index.
    //! @param[in] variableNumber    variable's index.
    string GetProbeVariableName(int variableNumber) const;

    //! Get set of variables.
    const set<int> & GetProbeVariables() const;

    map<int, string> probeVariablesMap;

    Grid * GetGrid() const;

private:
    //! Geometry information
    Grid *grid;

    //! Grid type;
    int gridType;

    //! Number of probe variables.
    int nProbeVariables;

    //! Type of probe variables, each number represents ONE variable.
    set<int> probeVariables;

    //! Store probe variables.
    Data_Field * probeVariablesPtr;
private:
    //! Free all allocated memory by probe variables' pointer.
    void FreeAllProbeVarPtr();
    //! Free probe variables memory.
    //! @param[in] name    variables' name.
    void DeleteProbeVarPtr(const string &name);
};

const int PROBE_DENSITY     = 0;
const int PROBE_U           = 1;
const int PROBE_V           = 2;
const int PROBE_W           = 3;
const int PROBE_PRESSURE    = 4;
const int PROBE_TEMPERATURE = 5;

#include "Post_Probes.hxx"

class DumpToFileASCII
{
private:
    vector <DataContainer *> dataList;
    string fileName;

public:
    DumpToFileASCII(vector<DataContainer *> datalist_in, string filename_in);
    ~DumpToFileASCII();

public:
    void Run();

    //! Set the variables.
    void SetVariables(string Variables_in) { this->Variables = Variables_in; }

    //! Get the variables.
    string GetVariables() { return Variables; }

    //! Set number of variables in file.
    void SetnVariablesToDump(int nVariable_in) { this->nVariable = nVariable_in; }

    //! Return number of variables in file.
    int GetnVariablesToDump() { return nVariable; }

    //! Get the file name.
    string GetFilename() { return fileName; }

private:
    //! Init the probe variables need to dump.
    void InitProbeVariables();

    //! Write the probe variables data into file.
    void WriteProbeFlowFile();

private:
    int nVariable;
    string Variables;
};

}