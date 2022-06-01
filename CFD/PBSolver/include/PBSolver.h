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
//! @author    Ming Pingjian, Sun Huawen.

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
    virtual void AllocateGlobalVar(Grid * grid);
    virtual void DeAllocateGlobalVar(Grid *grid);
    virtual void ZeroResiduals(Grid * grid);
    virtual void TimeStep(Grid * grid);
    virtual void ComputeResidual(Grid *grid, FieldProxy *rhs_proxy);
    virtual void SourceFlux(Grid *grid);
    virtual void InviscidFlux(Grid *grid);
    virtual void ViscousFlux(Grid *grid);
    virtual void RungeKutta(Grid *grid, FieldProxy *rhs_proxy);
    virtual void Post();
    virtual void PostSolve(Grid * grid, int stage, int level);
    virtual RDouble UnsteadyConvergence(Grid* grid);


    virtual void InitMemory() {};
    virtual void AllocateMemory() = 0;
    virtual void InitField() = 0;
    virtual void ReinitField() = 0;
    virtual void UpdateBoundary() = 0;
    //virtual bool SolutionIsConverged() = 0;
	virtual void OutputResidual(ofstream& file, int flag)=0;
	
    inline void SetControllerName(string name) { this->controllerName = name; }
    inline std::string GetControllerName() { return this->controllerName; }
    inline void SetNSolvers(int nSolvers) { this->nSolvers = nSolvers; }
    inline int GetNSolvers() { return this->nSolvers; }
    inline void SetPrintName(string name) { this->printName = name; }
    inline std::string GetPrintName() { return this->printName; }

public: //! Extended parallel communication
    virtual void CommunicateInterfaceData(vector<string>& nameList);
    virtual void UploadInterfaceData(vector<string>& nameList);
    virtual void DownloadInterfaceData(vector<string>& nameList);
    virtual void UploadInterfaceValue(vector<string>& nameList);
    virtual void DownloadInterfaceValue(vector<string>& nameList);
    virtual void UpdateInterfaceData();
    virtual void CommunicationInterfaceData();
    virtual void FillActionKey(ActionKey* actkey, int action, int level);

private:
    std::string printName;            //! Screen output variable names
    std::string controllerName;       //! Name of controllers
    int nSolvers;                     //! Number of solvers
};

DEFINE_PBOBJECT_FACTORY(PBSolver);
}