//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//          PPPPP  H   H  EEEEE  N    N  GGGGG  L      EEEEE  III         +
//          P   P  H   H  E      NN   N  G      L      E       I          +
//          PPPPP  HHHHH  EEEEE  N N  N  G  GG  L      EEEEE   I          +
//          P      H   H  E      N  N N  G   G  L      E       I          +
//          P      H   H  EEEEE  N    N  GGGGG  LLLLL  EEEEE  III         +
//------------------------------------------------------------------------+
//          Platform for Hybrid Engineering Simulation of Flows           +
//          China Aerodynamics Research and Development Center             +
//                     (C) Copyright, Since 2010                          +
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      ParticlePointGroup.h
//! @brief     It is the base class of Particle Point Group.
//! @author    Lei Yinghaonan (Lanzhou University).

#pragma once
#include "OnePointVariable.h"
#include "IndexCell.h"
#include "DataContainer.h"
namespace PHSPACE
{

//! Typedef for particle ID's iter.
typedef map<int, OnePointVariable*>::iterator ParticleIDIter;

class ParticlePointGroup
{
private:
    //! The particle map
    //! The first varibal is the global id for particle for all zones.
    //! The second varibal is each particle variable on this zones.
    map<int, OnePointVariable*> *particlePointVariable;

    //! The mapping relationship of cell local id and particle global id .
    IndexCell* cellParticleID;

public:
    //! Constructor function, new map and new IndexCell.
    ParticlePointGroup();

    //! Destructor function , delete the map and IndexCell.
    ~ParticlePointGroup();

public:
    //! Init number of Cell and new map of cellParticleID.
    void InitCellIndex(Grid* grid);

    //! Add particle ID whitin search particle.
    //! Return the cell ID of particel
    //! 1. Search the cell of particle.
    //! 2. Get the ID of cell including particle.
    //! 3. Add the particeID ot cellParticleID(IndexCell).
    //! 3. Set the celllID of this particle on particlePointVariable.
    //! 4. Add this particle to particlePointVariable(Map).
    int AddParticle(Grid* grid, int idParticle, OnePointVariable* onePointVariable, bool &inZone);
    void QuickAddParticle(int idParticle, int IDCellOfParticle, OnePointVariable* onePointVariable, bool &inZone);

    void RemoveParticle(int particleID);

    //! Get particle variable.
    OnePointVariable* GetOnePointVariableByParticleID(int particleID);

    //! Get particle coordinate.
    SPDouble* GetParticleCoordinate(int particleID);

    //! Get particle velocity.
    SPDouble* GetParticleVelocity(int particleID);

    //! Get particle coordinate on RDouble**.
    //! This is used by dump HDF5.
    RDouble* GetParticleDiameter(RDouble refLength);
    RDouble* GetParticleDensity(RDouble refDensity);

    RDouble* GetParticleReNum();
    RDouble* GetParticleMaNum();
    RDouble* GetParticleSpecificHeatCapacity(RDouble refVelocityDimensitional, RDouble refTemperatureDimensitional);
    RDouble* GetParticleConvectiveHeatTransfer(RDouble refDensityDimensitional, RDouble refVelocityDimensitional);
    RDouble* GetParticleRadiativeHeatTransfer(RDouble refDensityDimensitional, RDouble refVelocityDimensitional);
    RDouble* GetParticleNuseltNum();
    RDouble* GetParticleWork(RDouble refDensityDimensitional, RDouble refVelocityDimensitional);
    RDouble* GetParticleDissipation(RDouble refDensityDimensitional, RDouble refVelocityDimensitional);

    RDouble* GetParticleCoordinateDim(RDouble refLength,int iDim);
    RDouble* GetParticleVelocityDim(RDouble refVelocity,int iDim);
    RDouble* GetParticleRelativeVelocityDim(RDouble refVelocity, int iDim);
    RDouble* GetParticleAngularVelocityDim(RDouble refAngularVel,int iDim);
    RDouble* GetParticleAccelerationDim(RDouble refAcceleration,int iDim);
    RDouble* GetParticleTemperatureDim(RDouble refTemperatuer,int iDim);

    RDouble* GetParticleTotalForceSt();
    RDouble* GetParticleStokesSt();
    RDouble* GetParticleModifiedSt();

    RDouble* GetParticleForceDim(RDouble refForce,int forceType, int iDim);

    //! Get the global particle ID list of cell by local cellID.
    int* GetIDOfParticleOfCell(int cellID);

    //! Get global particleID list on current local zone.
    int* GetGlobalParticleIDOnLocalZone();

    //! Get IndexCell.
    IndexCell* GetIndexCell();

    //! Get the cell ID of particle
    int GetCellIDOfParticle(int particleID);

    //! Get the global particle ID by local particleID
    int GetGlobalParticleIDByLocalIParticle(int iParticle);

    //! Get number of particle on current zone.
    int GetNumOfLocalParticle();

    //! Get number of cell on current zone.
    int GetNumOfLocalCell();

    //! Get the number of particle on cell by local cellID.
    int GetNumOfParticleOfCell(int cellID);

    //! Print particle variable.
    void PrintParticleVariable2Window(int particleID);

    //! Get the begin of iter.
    ParticleIDIter GetIterBegin();

    //! Get the end of iter.
    ParticleIDIter GetIterEnd();

    bool isParticleInZone(int particleID);

    void UpdateParticle(Grid*grid);

    //! Add particle variable to map.
    void AddParticle2Map(int particleID, OnePointVariable* onePointVariable);

    void CompressData(StructGrid* grid, DataContainer* &cdata, int &particleID, int &zoneIDGlobalRecv, int &size);

    //! Get particlePointVariable.
    map<int, OnePointVariable*>* GetMap();

    OnePointVariable* GetOnePoint(int particleID);

private:
    //! Delete the pointer for map.
    void DeleteMap();
};

}