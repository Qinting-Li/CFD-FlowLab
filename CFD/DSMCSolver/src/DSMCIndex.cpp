#include "DSMCIndex.h"
#include <stdio.h>
using namespace std;

namespace PHSPACE
{

Index::Index()
{
    moleculeNoInOrderOfCells = NULL;
    startAddressOfMoleculeInCells = NULL;
    numberOfMoleculeInCells = NULL;
}

Index::~Index()
{
    delete [] startAddressOfMoleculeInCells ;
    delete [] numberOfMoleculeInCells;
    delete [] moleculeNoInOrderOfCells;

    moleculeNoInOrderOfCells = NULL;
    startAddressOfMoleculeInCells = NULL;
    numberOfMoleculeInCells = NULL;
}

}



