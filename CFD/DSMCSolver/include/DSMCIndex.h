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
//! @file      DSMCIndex.h
//! @brief     Index operator for DSMC.
//! @author    Li Jin.

#pragma once

namespace PHSPACE
{

class Index 
{
private:   
    int * moleculeNoInOrderOfCells;
    int * startAddressOfMoleculeInCells;
    int * numberOfMoleculeInCells;

public:
    Index();
    ~Index();
    int * GetStartAddressOfMoleculeInCells() const {return this->startAddressOfMoleculeInCells;}
    void SetStartAddressOfMoleculeInCells(int size) {startAddressOfMoleculeInCells = new int [size];}
    int * GetNumberOfMoleculeInCells() const {return this->numberOfMoleculeInCells;}
    void SetNumberOfMoleculeInCells(int size) {numberOfMoleculeInCells = new int [size];}
    int * GetMoleculeNoInOrderOfCells() const {return this->moleculeNoInOrderOfCells;}
    void SetMoleculeNoInOrderOfCells(int size) {moleculeNoInOrderOfCells = new int [size];}
};

}



