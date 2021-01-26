#pragma once
#include "Pointer.h"
#include "BC_Struct.h"
#include "TK_Exit.h"
#include "MultigridFactory.h"
#include "MixGrid.h"

using namespace std;

namespace PHSPACE
{
class Region;

class GridGroup
{
public:
	GridGroup(int zoneStart = 0);
	~GridGroup();
private:
	int   nzones;
	int * block_proc_dump;		// Bell 20131124 add
	int * block_proc;
	int * block_proc_grid;
	int * block_type;
	int * block_idx ;
	int * block_fileProc;
	vector< Grid *> grids;
	int zoneStart;
	Region *region;
public:
	void AddGrid(Grid *grid);

public:
	bool IsZoneLayoutInitialized();
	void SetZoneHandler(Region *region) { this->region = region; };
	Region * GetRegion() { return region; };
	void SetNZones(int nZonesIn) { this->nzones = nZonesIn; }
	int GetNZones() const { return nzones; }	
	int GetNumberofGrid() const;
	Grid *GetGrid(int izone) const { return grids[izone]; }
	
	//void ReadSingleGrid(fstream &file, int izone);
	void InitZoneLayoutInformation(const string &filename);

	void SetBlockProcDump(int * block_proc_dump) { this->block_proc_dump = block_proc_dump; }
	void SetBlockProc(int * block_proc) { this->block_proc = block_proc; }
	void SetBlockProcGrid(int * block_proc) { this->block_proc_grid = block_proc; }
	void SetBlockType(int * block_type) { this->block_type = block_type; }
	void SetBlockIndex(int * block_idx) { this->block_idx = block_idx; }
	void SetBockFileProc(int * block_fileProc) { this->block_fileProc = block_fileProc; }

	int * GetBlockProcDump() { return this->block_proc_dump; }
	int * GetBlockProc()     { return this->block_proc; }
	int * GetBlockProcGrid() { return this->block_proc_grid; }
	int * GetBlockType()     { return this->block_type; }
	int * GetBlockIndex()    { return this->block_idx; }
	int * GetBlockFileProc() { return this->block_fileProc; }
	int   GetZoneStart()     { return this->zoneStart; }

	void SetGlobalZoneLayout();
private:
	void InitZoneLayout(fstream &file);
	
};
}