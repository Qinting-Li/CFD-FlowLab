#pragma once
#include "ElementFactory.h"
using namespace std;

namespace PHSPACE
{
class Region;
class CgnsZone;

class ZoneGridSolver
{
public:
    ZoneGridSolver(PHVector1D < CgnsZone * > & cgnsZoneContainer);
    ~ZoneGridSolver();
protected:
    Geo_PointFactory      * pointFactory;
    ElementPropertySolver * elementPropertySolver;
    ElementFaceSolver     * elementFaceSolver;
    PHVector1D < CgnsZone * > cgnsZoneContainer;
    Region * region;

    InterfaceInfo *interfaceInformation;
    PHVector2D < int > globalNeighborZones;

    RDouble minimumLength, maximumLength;
    int     zoneIndex;
    bool    hasBeenAdapted;
public:
    int  GetNumberOfElements()      { return this->GetElementPropertySolver()->GetNumberOfTotalElements(); }
    int  GetNumberOfFaces()         { return this->GetElementFaceSolver()->GetNumberOfFaces(); }
    int  GetNumberOfBoundaryFaces() { return this->GetElementFaceSolver()->GetNumberOfBoundaryFaces(); }
    int  GetNumberOfTotalFaces()    { return this->GetElementFaceSolver()->GetNumberOfTotalFaces(); }
    int  GetNumberOfTotalElements() { return this->GetNumberOfElements() + this->GetNumberOfBoundaryFaces(); }

    void SetZoneHandler(Region * region) { this->region = region; }
    void SetAdaptStatus(const bool & statusIn) { this->hasBeenAdapted = statusIn; }
    bool GetAdaptStatus() const { return hasBeenAdapted; }
    void SetInterfaceInformation(InterfaceInfo * interfaceInformationIn) { this->interfaceInformation = interfaceInformationIn; }
    PHVector2D < int > & GetGlobalNeighborZones() { return this->globalNeighborZones; }
    InterfaceInfo * GetInterfaceInformation() { return interfaceInformation; }
    void    SetZoneIndex(const int & indexIn) { this->zoneIndex = indexIn; }
    int     GetZoneIndex() const     { return zoneIndex; }
    RDouble GetMinimumLength() const { return minimumLength; }
    RDouble GetMaximumLength() const { return maximumLength; }
    void ComputeMinimumAndMaximumElementEdgeLength();

    Geo_PointFactory          * GetPointFactory         () { return pointFactory;          };
    ElementFaceSolver     * GetElementFaceSolver    () { return elementFaceSolver;     };
    ElementPropertySolver * GetElementPropertySolver() { return elementPropertySolver; };

    void GenerateComputationalGrid            (Grid *grid, int gridType);

    void InitializeAllCgnsZoneElements ();
    void ScanElementBoundaryFace();
    void ScanNewElementBoundaryFace();

    void GenerateInformationForComputationalElement();
    void PrepareUnstructuredComputationalGrid();
    void CompressStatusIntoDataContainer(DataContainer *& dataContainer, int neighborZoneIndex);
    void CompressFaceRBFWeightIntoDataContainer(DataContainer *& dataContainer, int neighborZoneIndex);
    void AnasysStatusFromDataContainer(DataContainer *& dataContainer, int neighborZoneIndex);
    void AnasysFaceRBFWeightFromDataContainer(DataContainer *& dataContainer, int neighborZoneIndex);
    void GetSourceCellIndex(int iFace, int iPosition, int & sourceIndex);
    void GetTargetIndex(int iFace, int iPosition, int & targetIndex);
    void GetSourceFaceIndex(int iFace, int & sourceIndex);
    void GetTargetFaceIndex(int iFace, int & targetIndex);
    void ResetInterfaceElementProperty();

    void CompressRefineType(DataContainer *& dataContainer, int neighborZoneIndex);
    void DecompressRefineType(DataContainer *& dataContainer, int neighborZoneIndex);

    void CompressCellProjectType(DataContainer *& dataContainer, int neighborZoneIndex);
    void DecompressCellProjectType(DataContainer *& dataContainer, int neighborZoneIndex);
};

typedef void (*CompressData)(ZoneGridSolver * zoneGridSolver, DataContainer *& dataContainer, int neighborZoneIndex);
typedef void (*DeCompressData)(ZoneGridSolver * zoneGridSolver, DataContainer *& dataContainer, int neighborZoneIndex);

void SwapInterfaceData(Region * region, ZoneGridSolver ** zoneGridSolverList, CompressData compressData, DeCompressData decompressData);

}