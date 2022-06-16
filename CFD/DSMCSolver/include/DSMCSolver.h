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
//! @file      DSMCSolver.h
//! @brief     Solver of DSMC.
//! @author    Li Jin.

#pragma once
#include "Precision.h"
#include <vector>
#include <string>
using namespace std;

namespace PHSPACE
{

class DSMC 
{
private:
    class Region * region;
    class TimeSpan * testTime;

    class DSMCUnstructGrid ** dsmcgrid;
    class Molecule ** molecule;
    class Sample ** sample;
    class Index ** index;
    class Collision ** collision;

    int ZoneNum,ZoneLeft,ZonePerCPU,ZoneStart,ZonePerCPU1,ZoneEnd,WallSurfNum;
    int MyRank,numCPU,ierr;
    RDouble timeStep;
    vector< int > ZoneNo,Zone2CPU;
    vector< RDouble > FNUMI,DTMI,TIME,TIMI;
    int NSMP;
    RDouble *** SPM;
    vector< int > ISPV;

    struct Species 
    {          // info on each particle species
        string id;
        RDouble percent;
        RDouble molwt;
        RDouble mass;
        RDouble rotrel;
        RDouble vibrel;
        RDouble vibtemp;
        RDouble specwt;
        RDouble charge;
        int rotdof,vibdof;
        int internaldof;
        RDouble diam;
        RDouble omega;
        RDouble tref;
        RDouble alpha;
        RDouble rotc1;
        RDouble rotc2;
        RDouble rotc3;
        RDouble vibc1;
        RDouble vibc2;
    };
    Species * species;         // list of particle species info
    int nSpecies;            // # of defined species
    int MMVM;
    int nweight;

    RDouble preCollisionEtran,preCollisionErot,preCollisionEvib,preCollisionEtotal,preCollisionEint;
    RDouble postCollisionEtran,postCollisionErot,postCollisionEvib,postCollisionEtotal,postCollisionEint;

public:
    DSMC();
    ~DSMC();
public:     
    void Run();
    void PreProcess();
    void InitSample();
    void MoleculeMove(){}
    void MoleculeIndex();
    void MoleculeCollision();
    void FieldSample();
    void PostProcess();
    void WriteResumeFiles(){}

    void AllocateMemory();
    void DeleteMemory();
    void ReadGeometry();
    void UpdateAllZonesPara();
    void InitGeometry();
    void ZoneToCPU();
    
    void ReadRestart(){}

    void ReadSpecies();
    void CalTimeStep();
    void CalRealtoSimRatio();
    void InitData();
    void SelectMolecule(int,int,int & ,int &);
    void NonReactingEDisposal(int,int,int);
    void Larsen_Borgenakke(RDouble,RDouble,RDouble &);
    void TwobodyScattering(int,int,int,RDouble,int &);

    RDouble *** GetSPM() const { return SPM;}
    int GetnSpecies() const { return nSpecies;}
    int GetMMVM() const { return MMVM;}   
    int Getnweight() const { return nweight;} 
    RDouble GetTimeStep() const { return timeStep;}

    DSMCUnstructGrid ** Getdsmcgrid() const {return dsmcgrid;}
    Molecule ** GetMolecule() const {return molecule;}

    void SetZoneStart(int ZoneStart){this->ZoneStart=ZoneStart;}
    void SetZoneEnd(int ZoneEnd){this->ZoneEnd=ZoneEnd;}
    int GetZoneStart()const{return ZoneStart;}
    int GetZoneEnd()const{return ZoneEnd;}

    void CPUtoZone();
    const vector< int > & GetZoneNo() { return ZoneNo; }
    const vector< RDouble > & GetFNUMI() { return FNUMI; }
    const vector< RDouble > & GetDTMI() { return DTMI; }
    const vector< RDouble > & GetTIME() { return TIME; }
    const vector< RDouble > & GetTIMI() { return TIMI; }
    const vector< int > & GetZone2CPU() { return Zone2CPU; }

};
    
    RDouble gamma(RDouble);
    RDouble gamma(int);
    RDouble myrand();

}






