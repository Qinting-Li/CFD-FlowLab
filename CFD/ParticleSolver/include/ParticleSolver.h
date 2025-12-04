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
//! @file      ParticleSolver.h
//! @brief     The base particle solver.
//! @author    Lei Yinghaonan (Lanzhou University).

#pragma once
#include "CFDSolver.h"
#include "LIB_Macro.h"

namespace PHSPACE
{
class Param_ParticleSolver;

class ParticleSolver : public CFDSolver
{
public:
    ParticleSolver();
    ~ParticleSolver();

public:
    LIB_EXPORT Param_ParticleSolver *GetControlParameters() const;

    //! Pure virtual function from PHSolver,
    //! but given the function body by CFDSolver.
    virtual void InitMemory();
    virtual void InitFlow();

    //! Virtual function from PHSolver,
    //! but given the function body by CFDSolver.
    virtual void CommunicationInterfaceData();
    virtual void CommunicationInterpointData();
    virtual void Action(ActionKey *actkey);
    virtual void TranslateAction(ActionKey *actkey);
    virtual streamsize TranslateActionLength(ActionKey *actkey);
    virtual void PrepareInterfaceData(Grid *grid, Data_ParamFieldSuite *datastore, InterfaceDataProxy *interfaceDataProxy);

    //! Pure virtual function from PHSolver,
    //! with function body by CFDSolver but empty.
    virtual void Solve();
    virtual void ComputePostVisualVariables(Post_Visual *postVisualization);

    //! Pure virtual function from PHSolver,
    //! with function body by CFDSolver but empty.
    virtual void InitCoarseGridsFlow();
    virtual void ActionReflect(ActionTag *acttag);
    virtual void ComputePostProbesVariables(Post_Probes *postProbesVar);

    //! Pure virtual function from PHSolver,
    //! but without function body in CFDSolver.
    virtual int GetNumberOfEquations();

    //! Pure virtual function from PHSolver,
    //! but didn't appear in CFDSolver.
    virtual void Post();

    //! Virtual function from CFDSolver
    //! for output resualtfile
    virtual void DumpResultFile(Grid *grid, int level);
    //! for RK
    virtual void RungeKuttaResidual(Grid *grid, FieldProxy *dqProxy, RDouble coef);

    //! For Controller solving.
    void ZeroResiduals(Grid* grid);
};

}