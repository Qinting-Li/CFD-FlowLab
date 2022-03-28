#pragma once
#include "Precision.h"
#include "Force.h"
#include "Pointer.h"
#include "AMRDef.h"
#include "Data_ParamFieldSuite.h"

using namespace std;

namespace PHSPACE
{
class SixDofSolverManager;
class MovingMeshManager;
class DeformingSolverManager;
class AleForceManager;
class Grid;
class Data_ParamFieldSuite;

class AleManager
{
public:
    AleManager();
    ~AleManager();

private:
    AleForceManager     *aleForceManager;
    SixDofSolverManager* sixDofSolverManager;
    DeformingSolverManager* deformingManager;
    MovingMeshManager   *movingMeshManager;

public:
    void Initialize();
    void Allocate();
    void Restart();
    void Run();
    void UpdateAleForce();
    void Post();
    void Dump();
    void DumpRestart();

protected:
    void InitializeDimensionalReferenceTime();
    fstream* currentAleFile;

public:
    fstream & GetCurrentAleFile() { return * currentAleFile; };
    AleForceManager        * GetAleForceManager()        { return aleForceManager;        };
    SixDofSolverManager    * GetSixDofSolverManager()    { return sixDofSolverManager;    };
    DeformingSolverManager* GetDeformingSolverManager()  { return deformingManager; };
    void ComputeMetrics();
    void CommunicateCellCenterData();
    void GridVerticeVelocityGeneral(Grid* grid, RDouble* xt, RDouble* yt, RDouble* zt);
    void ComputeGridFaceVelocity();
};

//整体上的调用
void CreateAleManager();
void FreeAleManager();
void SolveSingleInnerIterationStepForAleSolvers();
void PostprocessAfterInnerIterationStepForAleSolvers();

AleManager * GetAleManager();

void BackupAllZonesOldGrid();

//ALE备份文件
bool IsReadAleFile();
void SetIsReadAleFile(bool isReadAleFileIn);
fstream & GetCurrentAleFile();

int      GetIntegerParameterFromDataBase(int bodyIndex, const string& parameterName);
RDouble  GetRDoubleParameterFromDataBase(int bodyIndex, const string& parameterName);
string   GetStringParameterFromDataBase(int bodyIndex, const string& parameterName);

int*     GetIntegerArrayFromDataBase(const string& nameOfIntegerArray, int numberOfElements);
void     GetRDoubleVectorFromDataBase(vector< RDouble >& parameter, int bodyIndex, const string& parameterName, int numberOfElements);

RDouble   GetAleTimeStep();

bool     IsImplicitAleMethod();
bool     IsExplicitAleMethod();

template< typename T >
class SuperFieldNameProperty;

class AllocateParameter
{
public:
    AllocateParameter(): numberOfEquations(0),dimensions(0),numberOfBoundaryFaces(0),allocateBoundaryFaces(true){};
    ~AllocateParameter() {};
protected:
    int numberOfEquations;
    int dimensions;
    int numberOfBoundaryFaces;
    bool allocateBoundaryFaces;
public:
    int GetNumberOfEquations() { return numberOfEquations; };
    int GetDimensions() { return dimensions; };
    int GetNumberOfBoundaryFaces() { return numberOfBoundaryFaces; };

    void SetNumberOfEquations(int numberOfEquations) { this->numberOfEquations = numberOfEquations; };
    void SetDimensions(int dimensions) { this->dimensions = dimensions; };
    void SetNumberOfBoundaryFaces(int numberOfBoundaryFaces) { this->numberOfBoundaryFaces = numberOfBoundaryFaces; };

    void SetAllocateBoundaryFaces(bool allocateBoundaryFaces) { this->allocateBoundaryFaces = allocateBoundaryFaces; };
    bool GetAllocateBoundaryFaces() { return this->allocateBoundaryFaces; };
};

class FieldNameProperty
{
public:
    FieldNameProperty() {};
    ~FieldNameProperty() {};
protected:
    static std::map< string, int > fieldNamePropertyData;
public:
    static void AddFieldNameProperty(const string& fieldName, int numberOfEquations)
    {
        fieldNamePropertyData[fieldName] = numberOfEquations;
    };

    static int GetNumberOfEquations(const string& fieldName)
    {
        std::map< string, int >::iterator iter;
        iter = fieldNamePropertyData.find(fieldName);
        if (iter != fieldNamePropertyData.end())
        {
            return iter->second;
        }
        return -1;
    };
};

template< typename T >
class InterfaceFieldNameProperty
{
public:
    InterfaceFieldNameProperty();
    ~InterfaceFieldNameProperty();
protected:
    static std::map< string, AllocateParameter > fieldNameParameterData;
public:
    static void AddFieldNameProperty(const string& fieldName, int numberOfEquations);
    static void AddFieldNameProperty(const string& fieldName, int numberOfEquations, int dimensions);
    static int GetNumberOfEquations(const string& fieldName);
    static int GetNumberOfDimensions(const string& fieldName);
public:
    static void AllocateInterfaceVariables(int numberOfInterfaces, Data_ParamFieldSuite* dataStorage);
    static void DeAllocateInterfaceVariables(Data_ParamFieldSuite* dataStorage);
};

template< typename T >
void GreatAllocateOversetInterfaceVariables(Grid* gridIn);

template< typename T >
void GreatDeAllocateOversetInterfaceVariables(Grid* gridIn);

void DeAllocateOversetStorage();

void AllocateOversetStorage();

}