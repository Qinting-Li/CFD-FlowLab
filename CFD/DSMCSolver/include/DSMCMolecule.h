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
//! @file      DSMCMolecule.h
//! @brief     Molecule related operations of DSMC.
//! @author    Li Jin.

#pragma once
#include "Precision.h"

using namespace std;

namespace PHSPACE
{

class Molecule
{
private:
    int * moleculeLocation, * moleculeSpecies, * IPCP;
    int * moleculeVibrationLevel;
    RDouble ** moleculePosition, ** moleculeVelocity;
    RDouble * moleculeRotationalEnergy, * moleculeVibrationalEnergy;
public:
    Molecule();
    ~Molecule();
public:
    int * GetMoleculeLocation()const{return moleculeLocation;}

    int * GetMoleculeSpecies()const{return moleculeSpecies;}

    int * GetIPCP()const{return IPCP;}

    RDouble ** GetMoleculePosition() const {return moleculePosition;}

    RDouble ** GetMoleculeVelocity() const {return moleculeVelocity;}

    RDouble CalMoleculeRotationalEnergy(int rotationalDegreeOfFreedom, RDouble temperature);
    RDouble * GetMoleculeRotationalEnergy() const { return moleculeRotationalEnergy;}

    RDouble CalMoleculeVibrationalEnergy(int vibrationalDegreeOfFreedom, RDouble temperature, RDouble vibtemperature);
    RDouble * GetMoleculeVibrationalEnergy() const { return moleculeVibrationalEnergy;}

    int * GetMoleculeVibrationLevel() {return moleculeVibrationLevel;}

    void init();
};

}
