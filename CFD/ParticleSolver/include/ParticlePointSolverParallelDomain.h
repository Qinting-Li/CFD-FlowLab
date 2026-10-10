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
//! @file      ParticlePointSolverParallelDomain.h
//! @brief     The Particle point solver.
//! @author    Lei Yinghaonan (Lanzhou University).

#pragma once
#include "ParticleSolver.h"
#include "CFDSolver.h"
#include "Flux_Inviscid.h"
#include "LIB_Macro.h"
#include "ParticlePointGroup.h"
#include "LagrangianVariable.h"
#include "RungeKuttaMethod.h"
#include <iostream>
#include <time.h>
#include <iomanip>
#include <math.h>
#include <fstream>
#include<cstdio> 
#include<stdio.h>
#include "DataContainer.h"
#include "ctime"
#include "cstdlib"
#include <cmath>

#pragma warning(disable:4996)
using namespace std;

namespace PHSPACE
{
class Param_ParticlePointSolverParallelDomain;

class ParticlePointSolverParallelDomain : public ParticleSolver
{
public:
    //! Nothing to do.
    //! The variable init by new is not in this Constructor function.
    ParticlePointSolverParallelDomain();

    //! Destructor function.
    //! 1.free the variable on grid.
    //! 2.free the parameter on grid.
    ~ParticlePointSolverParallelDomain();
public:
    LIB_EXPORT Param_ParticlePointSolverParallelDomain *GetControlParameters() const;

    //! This is a virtual function by PHSolver.
    //! Init the memory for Param and ParticlePointGroup by new and init.
    void InitMemory();

    //! This is a virtual function by PHSolver.
    //! Initialize particle field.
    void InitFlow();

    //! This is a virtual function by CFDSolver.
    void DumpResultFile(Grid *grid, int level);

    //! Solving the equation of particle motion.
    void SolverMultiphaseFlow(FieldProxy *rhsProxy);

    //! Communicate particle variable.
    void CommunicateMultiphaseFlow(FieldProxy *rhsProxy);

private:
    //! =================================================================
    //! ===                  Used by InitMemory().                    ===
    //! =================================================================
    //! New and init Param_PartcilePointSolvers.
    void InitControlParameters();
    //! Read particle point hypara
    void ReadParameter();
    //! New particle point group.
    void AllocateGlobalVariables();
    //! Get ParticlePointGroup
    ParticlePointGroup *GetParticlePointGroup(string varibaleName);
    //! Get Lagrangian variable
    LagrangianVariable<RDouble> *GetLagVariableDouble(string varibaleName);
    //! Get number of particle on currentlocal zone
    int GetNumOfLocalParticle();
    int GetNumOfTotalParticle();

    //! =================================================================
    //! ===                  Used by InitFlow().                      ===
    //! =================================================================
    //! Judge if need to read restart particle file.
    bool JudgeIfRestart();
    //! Judge if need to write restart particle file.
    bool JudgeIfWriteRestart();
    //! Judge if need to read init particle file.
    bool JudgeIfInit();
    //! Read init particle file on hdf5.
    void ReadParticleH5(ActionKey *actkey);
    //! Fill action key for task in this solver.
    void FillActionKey(ActionKey *actkey, int action, int level);
    //! Get the init particle file name.
    const string GetInitFileName();
    void InitLagrangianVariable();
    //! Check sum number of particle for each zone.
    void CheckNumOfParticle();
    void UpdateNumOfParticle();
    //! Check particle cell.
    void CheckParticleCell();
    int QuickSearchParticleCell(Grid *grid, SPDouble &onePointCoordinate, bool &inZone);
    int SearchParticelForOrthogonalStrGrid(Grid *gridIn, SPDouble &onePointCoordinate, bool &inZone);
    void ReadParticleHDF5File(ActionKey *actkey);

    //! =================================================================
    //! ===       Used by DumpResultFile(Grid *grid, int level).      ===
    //! =================================================================
    //! This is a virtual function by CFDSolver.
    void DumpRestartData(ActionKey *actkey);
    //! Note here,CFDSolver::CreateH5RestartFile(ActionKey *actkey)
    //! is not a virtual function,which will be used by
    //! CFDSolver::DumpRestartData,which is a virtual function.
    //! So here, we creat a new function.
    void CreateH5ParticleFile(ActionKey *actkey);
    //! This is a virtual function by CFDSolver.
    void DumpRestartH5(ActionKey *actkey);

    //! =================================================================
    //! ===                Particle Interpolation                     ===
    //! =================================================================
    //! New version code.
    int *GetCellIndex(Grid *grid, int cellID);
    //! Calculate the fluid variables needed for particle interpolation.
    void CalculateVariableForParticleOnEular(ActionKey *actkey);
    void CalculateVariableForParticleOnEular(StructGrid *grid, ActionKey *actkey);
    void CalculateVariableForParticleOnEular(UnstructGrid *grid, ActionKey *actkey);

    //! Calc the flow at particle location after CalculateVariableForParticleOnEular.
    void CalcFlowOnLocationOfParticle(ActionKey *actkey);

    void ParticleInterpolation(ActionKey *actkey, Grid *grid, OnePointVariable *onePointVariable);

    //! actkey-subact : 
    //!   0 -- : init particle as flow.
    //!   1 -- : Interpolation in RK.

    //! Gimenez2019JCP
    //! Grad interpolation.
    void InterpolationGimenez2019JCP(ActionKey *actkey, Grid *grid, OnePointVariable *onePointVariable);
    void InterpolationGimenez2019JCP(ActionKey *actkey, StructGrid *grid,OnePointVariable *onePointVariable,
        RDouble4D &flowOnParticleInEular, RDouble4D &flowOnParticleInEularGrad, RDouble4D  &flowOnParticleInEularGradSecond);
    void InterpolationGimenez2019JCP(ActionKey *actkey, UnstructGrid *grid, OnePointVariable *onePointVariable,
        RDouble **flowOnParticleInEular, RDouble **flowOnParticleInEularGrad, RDouble **flowOnParticleInEularGradSecond);
    RDouble GradDot(SPDouble &relativeCoord, RDouble4D &gradient,int &iP,int &jP,int&kP,int index);
    RDouble GradDoubleDot(SPDouble &relativeCoord, RDouble4D &gradient2ord, int &iP, int &jP, int &kP, int index);

    //! TrilinearInterpolation
    //! Simple linear interpolation.
    RDouble LinearInterpolation(RDouble &x0, RDouble &f0, RDouble &x1, RDouble &f1, RDouble &x);
    RDouble CheckInterpolation(StructGrid *grid, SPDouble &particleCoordinate, int cellID, RDouble4D &q, int mVar);

    void CheckCornerAndGhost();

    //! =================================================================
    //! ===                   Particle Boundary                       ===
    //! =================================================================
    //! Particle boundary condition.
    bool ParticleBoundary(Grid *gridIn,OnePointVariable *onePointVariable, IndexCell *indexCell);

    //! This is only a temporary treatment, 
    //! and the particle boundary conditions need to be deleted and optimized later
    void CheckParticleBCForOne(SPDouble &oneParticleCoordinate);

    bool FullyElasticCollisionWallStrcut(StructGrid *grid, OnePointVariable *onePointVariable, IndexParticleOfCell *indexParticleOfCell, int &iBCFace);

    void ResetParticleWhenOut(Grid *grid, OnePointVariable *onePointVariable, IndexParticleOfCell *indexParticleOfCell);
    void AddParticleAsInitBC(Grid *grid);
    RDouble GetRandomNumber(RDouble ValueOfLeftOrMean, RDouble ValueOfRightOrStandardDeviation, int mode);
    RDouble GetMaxParIDLocal();
    //! =================================================================
    //! ===              Lagrangian Particle Track                   ===
    //! =================================================================
    //! Particle Track
    void ParticleTrack(Grid *grid, OnePointVariable *onePointVariable);
    void ParticleTrack(StructGrid *grid, OnePointVariable *onePointVariable);
    void ParticleTrack(UnstructGrid *grid, OnePointVariable *onePointVariable);

    //! Used by ParticleTrack(StructGrid).
    void GetIndexIJKStructOnTrackSurface(int &indexSurface, int &iM, int &jM, int &kM);

    //! =================================================================
    //! ===                 Parallel Communication                    ===
    //! =================================================================
    void CommunicationParticle(Grid *grid, FieldProxy *rhsProxy);
    void CommunicationParticle(StructGrid *grid, ParticlePointGroup *particlePointGroup, FieldProxy *rhsProxy);
    void CommunicationParticle(UnstructGrid *grid, ParticlePointGroup *particlePointGroup, FieldProxy *rhsProxy);
    void CalcNumOfParticleToSendStruct(FieldProxy *rhsProxy, ActionKey *action, DataContainer *cdataIn, int &nParticleIn, int &size_in, DataContainer *cdataOut, int &nParticleToSend, int &size_send, set<int> &sendIDSet);
    void DecompressParticleDataContainerStruct(DataContainer *cdata, int &iGlobalzone, int &particleID, int &iP, int &jP, int &kP, OnePointVariable *onePointVariable);
    int CreatDataContanier(int &nParticle, DataContainer *cdata);

    void AddParticle(FieldProxy *rhsProxy, StructGrid *grid, int &particleID, int &iP, int &jP, int &kP, OnePointVariable *onePointVariable);

    void CommunicationFlowOnParticleInEular(ActionKey *actkey, Grid *grid);
    //! The main function here is mpi communication exchange for variables.
    //! (note that the ghost and corner is included).
    //! The main purpose of the actionKey is to control basic task information.
    //! actkey->subact : neighbour or corner to communicate in mpi.
    //!   -- 0 : only communicate neighbor.
    //!   -- 1 : only communicate corner.
    //! actkey->ipos : the zone global index to recv.
    //! actkey->format : the variable type.
    //!   -- 0 : RDouble3D
    //!   -- 1 : RDouble4D
    //!   -- 2 : OnePointVariable
    //! actkey->filename : the variable name.
    //! actkey->kind : the plan.
    //!   -- 0 : var corner only.
    //!   -- 1 : var ghost and corner.
    //!   -- 2 : grad and grad 2 ghost and corner.
    void CommunicateSpecificVariable(ActionKey *actkey);
    
    //! CompressSpecificArray : compress variable to DataContainer.
    //! actkey->
    //! varType : variable type for mpi.
    //!   -- 0 : RDouble3D
    //!   -- 1 : RDouble4D
    //!   -- 2 : OnePointVariable
    void CompressSpecificArray(DataContainer *&dataContainer, ActionKey *actkey);
    void CompressSpecificArray(UnstructGrid *grid, DataContainer *&dataContainer, ActionKey *actkey);
    void CompressSpecificArray(StructGrid *grid,DataContainer *&dataContainer, ActionKey *actkey);

    void DecompressSpecificArray(DataContainer *&dataContainer, ActionKey *actkey);
    void DecompressSpecificArray(UnstructGrid *grid, DataContainer *&dataContainer, ActionKey *actkey);
    void DecompressSpecificArray(StructGrid *grid, DataContainer *&dataContainer, ActionKey *actkey);

    //! =================================================================
    //! ===                        Solver                             ===
    //! =================================================================
    //! Stores the start values of RK iteration for particle coord,velocity
    void StoresParticleStartValue(OnePointVariable *onePointVariable);
    //! Update RK.
    void UpdateRK4(Grid *grid,OnePointVariable *onePointVariable, RungeKuttaIndex *indexRK, Data_Param *particleParam);
    void UpdateRK1(Grid *grid, OnePointVariable *onePointVariable, RungeKuttaIndex *indexRK, Data_Param *particleParam);

    //! particle CFL number
    void CalcCFLNumber();
    RDouble CalcCFLNumber(Grid *grid,OnePointVariable *onePointVariable, RDouble &particleCFLThreshold);
    RDouble CalcCFLNumber(StructGrid *grid, OnePointVariable *onePointVariable, RDouble &particleCFLThreshold);
    RDouble CalcCFLNumber(UnstructGrid *grid, OnePointVariable *onePointVariable, RDouble &particleCFLThreshold);

    void InitSource(Grid *grid);

    //! Calculate Two-way coupling.
    void GetSourceForRestart();
    void GetSourceFromParticle2Flow(Grid *grid, OnePointVariable *onePointVariable, Data_Param *particleParam);
    void GetSourceFromParticle2Flow(StructGrid *grid, OnePointVariable *onePointVariable, Data_Param *particleParam);
    void GetSourceFromParticle2Flow(UnstructGrid *grid, OnePointVariable *onePointVariable, Data_Param *particleParam);
    void KernelBoundaryTransformation(RDouble fn[3], RDouble particleVelocity[3],  RDouble parVelNew[3]);
    void ConvertVectorToHomogeneous(RDouble originVector[3], RDouble homoVector[4]);
    void MirrorCoordinate(RDouble originPoint[3], RDouble wallDistance, RDouble Fn[3], RDouble mirrorPoint[3]);
    RDouble DistanceBetweenPointAndSurface(RDouble point[3], RDouble surfacePoint[3], RDouble fn[3]);


};

}