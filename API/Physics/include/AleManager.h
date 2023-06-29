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
    void SychonizeAleTime();
    void InitializeDimensionalReferenceTime();
    fstream* currentAleFile;

public:
    fstream & GetCurrentAleFile() { return * currentAleFile; };
    AleForceManager        * GetAleForceManager()        { return aleForceManager;        };
    SixDofSolverManager    * GetSixDofSolverManager()    { return sixDofSolverManager;    };
    DeformingSolverManager* GetDeformingSolverManager()  { return deformingManager; };
    void ComputeMetrics();
    void CommunicateCellCenterData();
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
int      GetIntegerParameterFromDataBaseIfExist(int bodyIndex, const string& parameterName, int defaultValve);
RDouble  GetRDoubleParameterFromDataBase(int bodyIndex, const string& parameterName);
string   GetStringParameterFromDataBase(int bodyIndex, const string& parameterName);

RDouble  GetRDoubleParameterFromDataBaseIfExist(int bodyIndex, const string &parameterName, RDouble defaultValue);

int*     GetIntegerArrayFromDataBase(const string& nameOfIntegerArray, int numberOfElements);
void     GetRDoubleVectorFromDataBase(vector< RDouble >& parameter, int bodyIndex, const string& parameterName, int numberOfElements);

RDouble   GetAleTimeStep();

bool     IsImplicitAleMethod();
bool     IsExplicitAleMethod();

void DeAllocateOversetStorage();

void AllocateOversetStorage();

}