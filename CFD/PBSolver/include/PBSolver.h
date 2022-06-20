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
//! @file      PBSolver.h
//! @brief     Pressure Based solver.
//! @author    Wan Yunbo, Ming Pingjian, Meng Liyuan, Sun Huawen, Bell.

#pragma once

#include "PBObject.h"
#include "CFDSolver.h"
#include "Parallel.h"
#include "Post_Visual.h"
#include "Post_Probes.h"
#include "PBIterSolv.h"

namespace PHSPACE
{
class Grid;

class PBIterSolv;

class PBSolver : public CFDSolver, public PBObject
{
public:
    PBSolver(const string &name);
    ~PBSolver();

public:
    void AllocateGlobalVar(Grid *grid);
    void DeAllocateGlobalVar(Grid *grid);
    void ZeroResiduals(Grid *grid);
    void TimeStep(Grid *grid);
    void ComputeResidual(Grid *grid, FieldProxy *rhs_proxy);
    void SourceFlux(Grid *grid);
    void InviscidFlux(Grid *grid);
    void ViscousFlux(Grid *grid);
    void RungeKutta(Grid *grid, FieldProxy *rhs_proxy);
    void Post();
    void PostSolve(Grid *grid, int stage, int level);
    RDouble UnsteadyConvergence(Grid *grid);


    virtual void InitMemory();
    virtual void AllocateMemory() = 0;
    virtual void InitField() = 0;
    virtual void ReinitField() = 0;
    virtual void UpdateBoundary() = 0;
    //virtual bool SolutionIsConverged() = 0;
    virtual void OutputResidual(ofstream &file, int flag)=0;

    inline void SetControllerName(string name) { this->controllerName = name; }
    inline string GetControllerName() { return this->controllerName; }
    inline void SetNSolvers(int nSolvers) { this->nSolvers = nSolvers; }
    inline int GetNSolvers() { return this->nSolvers; }
    inline void SetPrintName(string name) { this->printName = name; }
    inline string GetPrintName() { return this->printName; }

public: //! Extended parallel communication
    void CommunicateInterfaceData(vector<string> &nameList);
    void UploadInterfaceData(vector<string> &nameList);
    void DownloadInterfaceData(vector<string> &nameList);
    void UploadInterfaceValue(vector<string> &nameList);
    void DownloadInterfaceValue(vector<string> &nameList);
    void UpdateInterfaceData();
    void CommunicationInterfaceData();
    void FillActionKey(ActionKey *actkey, int action, int level);

private:
    string printName;            //! Screen output variable names
    string controllerName;       //! Name of controllers
    int nSolvers;                     //! Number of solvers
};

DEFINE_PBOBJECT_FACTORY(PBSolver);
}