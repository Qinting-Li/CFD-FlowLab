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
//! @file      Pre_HDF5File.h
//! @brief     It defines HDF5 file read and write method.
//! @author    Baka.

#pragma once
#include "IO_HDF5File.h"
#include "Region.h"

namespace PHSPACE
{

class IO_HDF5Write
{
private:
    Grid **grids;
    int nBlocks;
    string gridFileName;

public:
    //! Dump nBolcks grids to file.
    //! @param[in] gridFileName    Grid file name.
    //! @param[in] grids           The grids that need to be dumped.
    //! @param[in] nBlocks         The number of the grid zones.
    IO_HDF5Write(const string &gridFileName, Grid **grids, int nBlocks);
    ~IO_HDF5Write();

public:
    void Run();

private:
    hid_t CreateFile(const string &filename);
    void WriteVersionInfo(hid_t loc_id);
    void CreateTotalInfo();
    void WriteTotalInfo(hid_t loc_id);
    void WriteGridCoordinates(hid_t loc_id, Grid *gridIn);
    void WriteInterfaceInfo(hid_t loc_id, Grid *gridIn);

    void WriteStructGrid(hid_t loc_id, Grid *gridIn);
    void WriteBCRegion(hid_t loc_id, Grid *gridIn);
    void WriteStrBCName(hid_t loc_id, Grid *gridIn);

    void WriteUnstructGrid(hid_t loc_id, Grid *gridIn, int iZone);
    void WriteFaceTopology(hid_t loc_id, Grid *gridIn);
    void WriteCellTopology(hid_t loc_id, Grid *gridIn);
    void WriteUnstrBCName(hid_t loc_id, Grid *gridIn);
    void WriteInterpointInfo(hid_t loc_id, Grid *gridIn);

private:
    //! The grids int current processor
    Grid **gridsLocal;
    int *nPartOfEachProcessor;

    int *block_proc;
    int *block_idx;
    int *block_type;
    int *block_group;
};

class IO_HDF5Read
{
private:
    PHVectorString1D gridNameList;
    double version;

public:
    //! Read grid from file.
    //! @param[in] gridFileName    Grid file name.
    //! @param[in] region_in       The class to manage the geometry.
    IO_HDF5Read(PHVectorString1D &gridNameListIn, Region *region_in);
    ~IO_HDF5Read();

public:
    bool Run();

private:
    void CheckNumberofFiles(int iGrid);
    void InitGridFileNameListAndFileProcessor();
    hid_t OpenH5File(const string &filename);

    void ReadVersionInfo(hid_t loc_id);
    void ReadTotalInfoOld(hid_t loc_id);
    void ReadNumberOfZonesOld(hid_t loc_id);
    void ReadEachGridOld(hid_t loc_id);

    void ReadTotalInfo();
    void ReadNumberOfZones();
    void ReadEachGrid();

    void ReadGridCoordinates(hid_t loc_id, Grid *gridIn);
    void ReadInterfaceInfo(hid_t loc_id, Grid *gridIn);

    void ReadStructGrid(hid_t loc_id, int iZone);
    void ReadBCRegion(hid_t loc_id, Grid *gridIn);
    void ReadStrBCName(hid_t loc_id, Grid *gridIn);

    void ReadUnstructGrid(hid_t loc_id, int iZone);
    void ReadFaceTopology(hid_t loc_id, Grid *gridIn);
    void ReadCellTopology(hid_t loc_id, Grid *gridIn);
    void ReadUnstrBCName(hid_t loc_id, Grid *gridIn);
    void ReadInterpointInfo(hid_t loc_id, Grid *gridIn);

    void SetNZones(int nZonesIn) { this->nZones = nZonesIn; }
    int GetNZones() const { return nZones; }

    void RedistributeZonesIfNeed();
    bool IsZoneLayoutInitialized();
    void SetGlobalZoneLayout();
    void ClassifyGridSystem();
    void AddGrid(Grid *grid);
    void CreateZones();
    void ReadWallDist();

    void SetBlockProcDump(int *blockProcDumpIn) { this->blockProcDump = blockProcDumpIn; }
    void SetBlockProcGrid(int *blockProcGridIn) { this->blockProcGrid = blockProcGridIn; }
    void SetBlockProc(int *blockProcIn) { this->blockProc = blockProcIn; }
    void SetBlockType(int *blockTypeIn) { this->blockType = blockTypeIn; }
    void SetBlockIndex(int *blockIdxIn) { this->blockIdx = blockIdxIn; }
    int *GetBlockProcDump() { return this->blockProcDump; }
    int *GetBlockProc()     { return this->blockProc; }
    int *GetBlockProcGrid() { return this->blockProcGrid; }
    int *GetBlockType()     { return this->blockType; }
    int *GetBlockIndex()    { return this->blockIdx; }

private:
    //! The file index which the current processor charged.
    int fileID;
    bool isFileProcessor;
    bool isP2PMode;

    //! The processor index of the current one.
    int myProcID;

    //! The processor index of the serve.
    int serverProcID;

    //! Number of global grid files.
    int globalFileNumber;
    vector < vector<string> > localGridNameList;

    //! The total number of global zones.
    int nZones;

    //! Number of zones in each .fts file.
    VVInt nZonesList;

    uint_t numberOfGridGroups;

    int zoneStart;
    int presentGridGroup;
    int *blockProc;
    int *blockType;
    int *blockIdx ;
    int *blockProcDump;
    int *blockProcGrid;
    int *blockFileProc;
	int *blockGroup;
    int taskType;

    vector<Grid *> grids;
    Region *region;
};

LIB_EXPORT void DumpHDF5Grid(const string &gridFileName, Grid **grids, int nBlocks);
LIB_EXPORT bool ReadHDF5File(Region *region_in, PHVectorString1D &gridGroupNameList);

}