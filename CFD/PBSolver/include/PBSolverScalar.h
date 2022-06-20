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
//! @file      PBSolverScalar.h
//! @brief     Scalar.
//! @author    Wan Yunbo, Ming Pingjian, Meng Liyuan, Sun Huawen, Bell.

#pragma once

#include "PBSolver.h"

namespace PHSPACE
{
class PBGrad;
class PBBoundary;
class PBTransient;
class PBConvection;
class PBDiffusion;
class PBSource;
class UnstructGrid;

#define SOLVER_SCALAR SCALAR

class PBSolverScalar : public PBSolver
{
protected:
    //! discrete scheme
    map<string, PBGrad *> gradCalc; //! gradient calculation
    vector<PBBoundary *> boundCalc; //! Boundary condition calculation
    PBTransient *tranCalc; //! Transient calculation
    PBConvection *convCalc; //! Convection calculation
    PBDiffusion *diffCalc; //! Diffusion term calculation
    vector<PBSource *> sourceCalc; //! Source term calculation

    UnstructGrid *grid;

    string scalarName;
    string gradientName;
    string *scalarNameList;

protected:
    string gradCalcConfigNameSuffix;
    string boundaryCalcConfigName;
    string convCalcConfigName;
    string diffCalcConfigName;
    string sourceCalcConfigName;
    string tranCalcConfigName;
    string gradCalcConfigName;
    string setFieldCalcConfigName;

protected: //! From configuration parameters
    RDouble initPhi; //! initial value of scalar variable
    RDouble initRho; //! initial value of density
    RDouble initMu; //! initial value of viscosity
    RDouble urf; //! relaxation factor
    RDouble res; //! residual of scalar equation

protected: //! Global register
    RDouble *phi; //! The scalar value to be solved, the number of elements + the number of boundary element faces.
    RDouble *phib; //! The boundary scalar value start pointer.
    RDouble *phic; //! Scalar correction
    RDouble *rho; //! density.
    RDouble *rhob; //! Boundary density start pointer.
    RDouble *rhof; //!  density on face.
    RDouble *mu; //!  viscosity on face.
    RDouble *mub; //!  Boundary viscosity start pointer.
    RDouble *diffCoe; //! Diffusion coefficient.
    RDouble *diffCoeb; //!Boundary diffusion coefficient start pointer.

protected:
    RDouble *ap; //! diagonal coefficient of general transport equation.
    RDouble *ajLeft; //!Upper triangular coefficient.
    RDouble *ajRight; //!Lower triangular coefficient.
    RDouble *sc; //! right-hand side of general transport equation.

protected: //Variables are used in calculations
    RDouble *phiOld; //! scalar at last time step.
    RDouble *phiOldOld;
    RDouble *flux; //! flux on face.
    RDouble *faceFlux; //! Pointer to initial position of flux on internal cell surface.
    RDouble *dphidx; //! gradient in X direction.
    RDouble *dphidy;
    RDouble *dphidz;
    RDouble resNow; //! residual of scalar

    int maxIterStep;
    int minIterStep;
    int innerIter;

protected:
    //! Solver iter Pointer
    PBIterSolv *iterSolv; //! for scalar
    int maxSweep; //! max iter steps 
    RDouble iterSolvTol;//! tolerance of iteration.

protected:
    RDouble *deltL; //! weighting coefficient of face.
    RDouble *deltR; //! weighting coefficient of face.
    int nTotalInFaces; //! number of interior faces.
    inline const int &NTotalInFace() 
    {
        return this->nTotalInFaces;
    }

public:
    inline const string &ScalarName() const
    {
        return this->scalarName;
    }

    inline void SetScalarName(const string  &name)
    {
        this->scalarName = name;
    }

    inline const string &GradientPhiName() const
    {
        return this->gradientName;
    }

    inline void SetGradientPhiName(const string &name)
    {
        this->gradientName = name;
    }

public:
    PBSolverScalar(const string &name);
    ~PBSolverScalar();

public:
    void InitMemory();
    void AllocateMemory();
    void InitField();
    void ReinitField();
    void UpdateUnsteadyVariable(Grid *grid);
    bool SolutionIsConverged();
    virtual void OutputResidual(ofstream &file, int flag);
    void InitIncompressibleResultFile(Grid *grid, string fileName);
    void DumpIncompressibleResultFile(Grid *grid, string fileName);

public:
    virtual void Discrt();
    virtual void Update();
    virtual void UpdateBoundary();
    virtual void UpdateProperties() {};
    virtual void UpdateAfterIterloop() {};
    virtual void SetOutputs(vector<string> &phiNameList);

public:
    void InitMatrix();
    void CalcGrad();
    void CalcBoundaryConditions();
    void CalcTransientTerm();
    void CalcConvectionTerm();
    void CalcDiffusionTerm();
    void CalcSourceTerm();
    void PrimeEquations();
    void CalcRes();
    void Relax();
    void SolveScalarEquation(Region &region);

public:
    void UpdatePhi(string phiName, int cellBegin, int cellEnd, RDouble *phic, RDouble *phi);
    void SolutionIsConverged(int iZone, bool *flag);

};

DECLARE_PBSOLVER_REGISTER_FUNCTION(PBSolver, PBSolverScalar, SOLVER_SCALAR);
}
