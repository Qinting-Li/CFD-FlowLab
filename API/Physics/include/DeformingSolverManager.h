#pragma once

#include "AleManager.h"
#include "Force.h"
#include "SixDofManager.h"

namespace PHSPACE
{

namespace MORPHINGMANNER
{
    const int NOMPRPHING  = 1;
    const int FISH        = 2;
    const int INSECT      = 3;
    const int BIRD        = 4;
    const int MISSILE     = 5;
    const int OSCILLATE   = 6;
    const int MORPHINGAIRCRAFTFOLD = 7;
    const int MORPHINGAIRCRAFTTAILEDGEUP = 8;
}

class MassCharacter;

//该参数仅用来确定网格的运动
class DeformingParameter
{
public:
    DeformingParameter();
    ~DeformingParameter();
protected:
    int numberOfParameters;
    vector< RDouble > parameters;
public:
    DeformingParameter & operator = ( const DeformingParameter & rhs );
    bool operator == ( const DeformingParameter & rhs );
public:
    void SetNumberOfParameters( int numberOfParameters ) { this->numberOfParameters = numberOfParameters; };
    int  GetNumberOfParameters() { return numberOfParameters; };

    vector< RDouble > & GetParameters()   { return parameters; };
    void SetParameters(PHVectorRDouble1D& parameterVector );
};

DeformingParameter * GetDeformingParameter( int bodyIndex );

class DeformingSolverManager
{
public:
    DeformingSolverManager();
    ~DeformingSolverManager();

public:
    void Initialize();
    void Restart();
    void Run();
    void Post();
    void Dump();
    void DumpRestart( fstream & file );
};

}