#pragma once
#include "Region.h"

using namespace std;

namespace PHSPACE
{

class GridTranslator
{
private:
    int numberOfZones;
    int codeOfOversetGrid;
    vector<StructGrid *> *gridGroup;
    vector<StructFace *> *structFaceGroup;
public:
     GridTranslator();
     GridTranslator(Grid **&strgrids, int &nBlocks);
    ~GridTranslator();
public:
    void Run();
private:
    void ReadGridgenFiles();
    void ReadGridgenCoordinateFile();
    void ReadGridgenBoundaryFile();
    void ProcessBoundaryCondition();
    void GenerateStructFaceGroup();
    void ColorateStructFaceGroup();
    void CollectInterfaceInformation();
    void OutputMultiZoneComputationalGrids();
    void OutputGrid();
    void CheckGrid();
};

}