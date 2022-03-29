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
//! @file      UGKSSimulation.h
//! @brief     It defines UGKS's detail function of flow simulation.
//! @author    

#pragma once
#include "Post_Visual.h"
#include "Solver.h"
#include "Geometry.h"
#include "Param_UGKSSolver.h"
namespace PHSPACE
{
    class StructBC;
    class FluxParaProxy;
    class TimeSpan;
    //! @brief UGKSSolver defines UGKS solver's detail function of flow simulation.
    class UGKSSolver : public PHSolver
    {
    private:
        Param_UGKSSolver * controlParameters;
        Post_Visual * postVisualization;

    public:
        //! Constructor
        UGKSSolver();

        //! Destructor
        ~UGKSSolver();
    public:

        //! Allocate global variable.
        void InitMemory();

        void ReleaseMemory();

         //! Initialize inflow variables.
        void InitFlow();

        //! Interface data communications in intial state.
        void CommunicationInterfaceDataInit();

        //! Compute minimum time step.
        void ComputeMinTimeStep(Grid * grid_in, RDouble &minDt);

        //! Solve UGKS algorithm in one time step.
        void UGKSSolve(Grid *grid_in);

        //! Solve Post-processing algorithm in one time step.
        void PostSolve(Grid * grid);

        //! Solve Post-processing algorithm when flow converges to an ideal state.
        void Post();
        
        void VisulizeGhostCell(Grid * grid);

    private:
        void Solve();

        //! Initialize control parameters.
        void InitControlParameters();

        //! Free control parameters.
        void FreeControlParameters();

        //! Read parameters about limiter.
        void ReadParameter();

        //! Get ID of limiter in UGKSSolver.
        //! @param[in] limiter_name        name of limiter
        int GetUGKSLimiterID(const string &limiter_name);

        //! Allocate variables of visualization.
        void AllocateVisualizationVariables();

        //! Allocate variables of global variable.
        void AllocateGlobalVar();

        //! Get discrete velocity variable.
        //! @param[in] controlParameters    parameter of UGKSSolver 
        //! @param[out] us                  discrete velocity of u
        //! @param[out] vs                  discrete velocity of v
        //! @param[out] ws                  discrete velocity of w
        //! @param[out] wx                  Velocity integrate coeffient of u in x direction
        //! @param[out] wy                  Velocity integrate coeffient of v in y direction
        //! @param[out] wz                  Velocity integrate coeffient of w in z direction
        void GetDiscreteVelocityVariable(Param_UGKSSolver * controlParameters, RDouble * us, RDouble * vs, RDouble * ws, RDouble * wx, RDouble * wy, RDouble * wz);

        //! Get discrete velocity variable's coordinate range.
        //! @param[in] controlParameters    parameter of UGKSSolver 
        //! @param[out] uBegin              minimum discrete velocity of u 
        //! @param[out] vBegin              minimum discrete velocity of v
        //! @param[out] wBegin              minimum discrete velocity of w
        //! @param[out] uEnd                maximum discrete velocity of u
        //! @param[out] vEnd                maximum discrete velocity of v
        //! @param[out] wEnd                maximum discrete velocity of w
        void GetDiscreteVelocityRange(Param_UGKSSolver * controlParameters, int &uBegin, int &vBegin, int &wBegin, int &uEnd, int &vEnd, int &wEnd);

        //! Deallocate global variable.
        void DeallocateGlobalVar(); 

        //! Register an interface field (variable) into buffer which is loaded on grid, used in parallel communication.
        //! The interface field variable buffer is set to be ZERO.
        //! @param[in] name    variable name of the interface field that need to be communicated.
        //! @param[in] type    variable type of the interface field
        //! -# FYINT   : Integer type.
        //! -# FYFLOAT : Float type.
        //! -# FYDouble: Double type.
        //! @param[in] dimension  variable dimension
        //! -# 1 :  1-D array.
        //! -# 2 :  2-D array.
        //! -# n :  n-D array.
        void RegisterInterfaceField(const string & name, const int type, const int dimesion);

         //! Register an interface field (variable) into buffer in CFDSolver.
        void RegisterCFDSolverInterfaceField();
        
        //! Set grid derivative.
        void SetGridDerivativeNT();

        //! Set cell's thickness.
        void SetThickness();

        //! Set ghost cell's thickness.
        //! param[in] w     thickness of cell
        void GhostCellThickness3D(RDouble4D &w);

        //! InitFlow
        //! Judge as restart UGKSSolver.
        bool JudgeIfRestartUGKS();

        //! Initialize flow as restart.
        void InitFlowAsRestart();

        //! Intialize flow as reading restart.
        void InitFlowAsReadingRestart();
        
        //! InitFlowAsRestart
        //! Calculate heat flux of inflow.
        void GetHeatFlux();

        //! Compute heat flux. 
        //! @param[in] prim                     primitive variable
        //! @param[out] heatFlux                heat flux of primitive variable
        void ComputeHeatFlux(RDouble *prim, RDouble * heatFlux);

        //! Get equilibrium state by caculating a Maxwellian distribution.
        //! @param[in] cx                       temporary variable: discrete velocity of u minus macroscopic velocity of u
        //! @param[in] cx                       temporary variable: discrete velocity of v minus macroscopic velocity of v
        //! @param[in] cx                       temporary variable: discrete velocity of w minus macroscopic velocity of w
        //! @param[in] eig                      represent m/2kT, m is molecular mass, k is the Boltzmann constant, T is the temperature
        //! @param[in] density                  macroscopic density
        //! @param[out] uv2                     temporary variable: represent cx * cx + cy * cy + cz * cz
        //! @param[out] f0                      initial distribution function
        void GetGDistribution(RDouble cx, RDouble cy, RDouble cz, RDouble eig,
                              RDouble density, RDouble &uv2, RDouble &f0);
        
        //! Convert primitive variable to conservative variable.
        //! @param[in] prim                     primitive variable
        //! @param[in] gamma                    fluid constant
        //! @param[out] cons                    conservative variable
        void PrimitiveToConservative(RDouble *prim, RDouble *cons, RDouble gamma);
        
        //! Convert conservative variable to primitive variable.
        //! @param[in] cons                     conservative variable
        //! @param[in] gamma                    fluid constant
        //! @param[out] prim                    primitive variable
        void ConservativeToPrimitive(RDouble *cons, RDouble *prim, RDouble gamma);

        //! Initialize inflow coeffcients 
        void InitInfowCoefs();

        //! Initialize flow variable
        void InitFlowVariable();

        //! Get discrete heat flux modified equilibrium state in BGK-Shakhov.
        //! @param[in] cx                       temporary variable: discrete velocity of u minus macroscopic velocity of u
        //! @param[in] cx                       temporary variable: discrete velocity of v minus macroscopic velocity of v
        //! @param[in] cx                       temporary variable: discrete velocity of w minus macroscopic velocity of w
        //! @param[in] qx                       represent heat flux of x
        //! @param[in] qy                       represent heat flux of y
        //! @param[in] qz                       represent heat flux of z
        //! @param[in] eig                      represent m/2kT, m is molecular mass, k is the Boltzmann constant,  T is the temperature
        //! @param[in] density                  macroscopic density
        //! @param[out] prandtl                 prandtl number
        //! @param[out] f0                      initial distribution function
        void GetDiscreteShakhovFlux(RDouble cx, RDouble cy, RDouble cz, RDouble qx, RDouble qy, 
        RDouble qz, RDouble eig, RDouble density, RDouble prandtl, RDouble &f0);

        //! Interface data communications.
        void CommunicationInterfaceData();

        //! Intialize task and add task to simu.
        //! param[in] acttag                    tag of action        
        void ActionReflect(ActionTag *acttag);

        //! Fill actionkey.
        //! @param[in] actkey                   action key variable
        //! @param[in] action                   simulation action
        //! @param[in] level                    multigrid level
        void FillActionKey(ActionKey *actkey, int action, int level);

        //! Call action function according to action.
        //! @param[in] actkey                   action key variable 
        void Action(ActionKey *actkey);

        //! Prepare interface data for GetInterfaceData/TranslateInterfaceData function.
        //! @param[in] grid                     store information of grid
        //! @param[in] datastore                data storage for sending/receving
        //! @param[in] interfaceDataProxy       proxy for interface data 
        void PrepareInterfaceData(Grid * grid, Data_ParamFieldSuite *datastore, InterfaceDataProxy *interfaceDataProxy);
        
        //! Get interface data by interface topology  information.
        //! @param[in] interfaceInfo                    Interface topology information of interface faces
        //! @param[in] actkey                   represent action key variable
        //! @param[in] interfaceDataProxy       proxy for interface data 
        void GetInterfaceData  (InterfaceInfo *interfaceInfo, ActionKey *actkey, InterfaceDataProxy * interfaceDataProxy);
        
        //! Get interface data by ActionKey class.
        //! @param[in] actkey                   represent action key variable
        void GetInterfaceData  (ActionKey *actkey);

        //! Upload interface value of 3D array.
        //! @param[in] grid                     store information of structured grid.
        //! @param[in] f                        3D array variable pointer of which need to be uploaded
        //! @param[in] name                     3D array variable name of which need to be uploaded 
        void UploadInterfaceValue3DArray(StructGrid *grid, RDouble3D *f, const string &name);

        //! Download interface value of 3D array.
        //! @param[in] grid                     store information of structured grid.
        //! @param[in] f                        3D array variable pointer of which need to be downloaded
        //! @param[in] name                     3D array variable name of which need to be downloaded 
        void DownloadInterfaceValue3DArray(StructGrid *grid, RDouble3D *f, const string &name);

        //! Upload interface value of 3D array.
        //! @param[in] grid                     store information of structured grid.
        //! @param[in] f                        6D array variable pointer of which need to be uploaded
        //! @param[in] name                     6D array variable name of which need to be uploaded 
        //! @param[in] neqn                     6D array variable dimension of which need to be uploaded 
        void UploadInterfaceValue6DArray(StructGrid *grid, RDouble6D &f, const string &name, int neqn);
        
         //! Download interface value of 3D array.
        //! @param[in] grid                     store information of structured grid.
        //! @param[in] f                        6D array variable pointer of which need to be downloaded
        //! @param[in] name                     6D array variable name of which need to be downloaded 
        //! @param[in] neqn                     6D array variable dimension of which need to be downloaded 
        void DownloadInterfaceValue6DArray(StructGrid *grid, RDouble6D &f, const string &name, int neqn);

        //! Upload interface data by ActionKey class.
        //! @param[in] actkey                   represent action key variable
        void UploadInterfaceData  (ActionKey *actkey);

        //! Download interface data ActionKey class.
        //! @param[in] actkey                   represent action key variable
        void DownloadInterfaceData(ActionKey *actkey);

        //! Translate interface data when action is UPDATE_INTERFACE_DATA.
        //! @param[in] actkey        represent action key variable 
        void TranslateAction(ActionKey *actkey);

        //! Get interface data by ActionKey class. 
        //! @param[in] actkey                   represent action key variable
        void TranslateInterfaceData(ActionKey *actkey);

        //! Translate interface data by interface topology  information.
        //! @param[in] interfaceInfo                    Interface topology information of interface faces
        //! @param[in] actkey                   represent action key variable
        //! @param[in] interfaceDataProxy       proxy for interface data 
        void TranslateInterfaceData(InterfaceInfo *interfaceInfo, ActionKey *actkey, InterfaceDataProxy * interfaceDataProxy);
        
        void UploadInterfaceDataAction(bool isInit);
        void UpdateInterfaceDataAction();
        void DownloadInterfaceDataAction(bool isInit);
        
        //! Upload interface data of intial state by ActionKey class.
        //! @param[in] actkey                   represent action key variable
        void UploadInterfaceDataInit(ActionKey *actkey);

        //! Download interface data by ActionKey class.
        //! @param[in] actkey                   represent action key variable
        void DownloadInterfaceDataInit(ActionKey *actkey);

        //! Translate length of interface data when action is UPDATE_INTERFACE_DATA.
        //! @param[in] actkey                   represent action key variable 
        streamsize TranslateActionLength(ActionKey *actkey);

        //! Translate length of interface data by ActionKey class.
        //! @param[in] actkey                   represent action key variable
        streamsize TranslateInterfaceDataLength(ActionKey *actkey);

        //! Translate length of interface data by interface topology  information.
        //! @param[in] interfaceInfo                    Interface topology information of interface faces
        //! @param[in] actkey                   represent action key variable
        //! @param[in] interfaceDataProxy       proxy for interface data 
        streamsize TranslateInterfaceDataLength(InterfaceInfo *interfaceInfo, ActionKey *actkey, InterfaceDataProxy * interfaceDataProxy);

        //!UGKSSolveSimu--boundary.
        //! Boundary condition.
        //! @param[in] grid                     store information of grid
        void Boundary(Grid *grid);

        //! Restrict index of UGKSSolver.
        //! @param[in] i                        index of x direction
        //! @param[in] j                        index of y direction
        //! @param[in] k                        index of z direction
        //! @param[in] maxi                     maximum index of x direction
        //! @param[in] maxj                     maximum index of y direction
        //! @param[in] maxk                     maximum index of z direction
        void UGKSRestrIndex(int &i, int &j, int &k, int maxi, int maxj, int maxk);

        //! Inflow boundary condition. 
        //! @param[in] grid                     store information of grid
        //! @param[in] bcr                      store information of boundary condition
        void InFlowBC3D(Grid * grid_in, StructBC * bcr);

        //! Symmetry boundary condition. 
        //! @param[in] grid                     store information of grid
        //! @param[in] bcr                      store information of boundary condition
        void SymmetryBC3D(Grid * grid_in, StructBC * bcr);

        //! Outflow boundary condition. 
        //! @param[in] grid                     store information of grid
        //! @param[in] bcr                      store information of boundary condition
        void OutFlowBC3D(Grid * grid_in, StructBC * bcr);

        //! Mark solid boundary condition.
        //! @param[in] grid                     store information of grid
        //! @param[in] bcr                      store information of boundary condition
        void MarkSolidBC3D(Grid * grid_in, StructBC * bcr);

        //! Mark moving solid boundary condition.
        //! @param[in] grid                     store information of grid
        //! @param[in] bcr                      store information of boundary condition
        void MarkMovingSolidBC3D(Grid * grid_in, StructBC * bcr);

        //! Calculate ghost cell solid boundary.
        void CalculateGhostCellSolidBoundary(Grid * grid_in);

        //! Mark solid boundary condition.
        //! @param[in] grid                     store information of grid
        //! @param[in] bcr                      store information of boundary condition
        void SolidBC3D(Grid * grid_in, StructBC * bcr);

        //! Mark moving solid boundary condition.
        //! @param[in] grid                     store information of grid
        //! @param[in] bcr                      store information of boundary condition
        void MovingSolidBC3D(Grid * grid_in, StructBC * bcr);
        
        //! Corner point of ghost cell.
        //! @param[in] grid                     store information of grid
        void CornerPoint(Grid * grid_in);

        //! Fill K direction in 2D for output
        //! @param[in] q                        4D variable
        //! @param[in] ni                       num of x direction
        //! @param[in] nj                       num of y direction
        //! @param[in] nk                       num of z direction
        //! @param[in] neqn                     num of equation
        void FillKDirection2D(RDouble4D & q, int ni, int nj, int nk, int neqn);

        //! Calculate viscidity and timestep. 
        //! @param[in] grid                     store information of grid
        void GetCTAndVisl(StructGrid *grid);

        //! Convert primitive variable to conservative variable.
        //! @param[in] prim                     primitive variable
        //! @param[in] gamma                    fluid constant
        //! @param[out] cons                    conservative variable
        void PrimitiveToConservativeBoundary(RDouble *prim, RDouble *cons,RDouble gamma);

        //! Get face normal vector.
        //! @param[in] grid                     store information of grid
        //! @param[in] i                        index of x direction
        //! @param[in] j                        index of y direction
        //! @param[in] k                        index of z direction
        //! @param[in] nx                       face normal vector of x direction 
        //! @param[in] ny                       face normal vector of y direction
        //! @param[in] nz                       face normal vector of z direction
        //! @param[in] nt                       face grid derivative 
        //! @param[in] nsurf                    face direction
        //! -# 0:  default option, for unstructured grid.
        //! -# 1:  face of x direction
        //! -# 2:  face of y direction
        //! -# 3:  face of z direction
        void GetFaceNormal(StructGrid * grid, int i, int j, int k, 
            RDouble &nx, RDouble &ny, RDouble &nz, RDouble &nt, int nsurf);

        //! Caculate UGKS's flux.
        //! @param[in] grid                     store information of grid
        void UGKSFlux(StructGrid *grid);

        //! Calculate UGKS's flux of certain direction.
        //! @param[in] xfn                      face normal vector of x direction
        //! @param[in] yfn                      face normal vector of y direction
        //! @param[in] zfn                      face normal vector of z direction
        //! @param[in] area                     face area
        //! @param[in] qcvs                     conservative variable 
        //! @param[in] dq                       change of conservative variable in one time step
        //! @param[in] thickness                cell thickness
        //! @param[in] f                        face distribution function
        //! @param[in] df                       change of face distribution function in one time step
        //! @param[in] st                       minimum index in certain direction 
        //! @param[in] ed                       maximum index in certain direction
        //! @param[in] i                        index of x direction
        //! @param[in] j                        index of y direction
        //! @param[in] k                        index of z direction
        //! @param[in] nsurf                    face direction
        //! -# 0:  default option, for unstructured grid.
        //! -# 1:  face of x direction
        //! -# 2:  face of y direction
        //! -# 3:  face of z direction
        void UGKSFluxCertainDirection(FluxParaProxy * fluxParaProxy, RDouble4D &xfn, RDouble4D &yfn, RDouble4D &zfn, RDouble4D &area, RDouble4D &qcvs, 
                                      RDouble4D &dq, RDouble4D &thickness, RDouble6D &f, RDouble6D &df,
                                      int st, int ed, int i, int j, int k, int nsurf);

        //! Set dq and df to zero.
        //! @param[in] grid                     store information of grid
        //! @param[in] dq                       change of conservative variable in one time step
        //! @param[in] df                       change of face distribution function in one time step
        void zeroDqDf(StructGrid * grid,RDouble4D &dq, RDouble6D &df);

        //! Calculate flux of certain direction in line.
        //! @param[in] grid                     store information of grid
        //! @param[in] dq                       change of conservative variable in one time step
        //! @param[in] df                       change of face distribution function in one time step
        void CalFluxLine(FluxParaProxy * fluxParaProxy, int nMax);

        //! Set faceFlux and fFaceFlux to zero.
        //! @param[in] grid                     store information of grid
        //! @param[in] faceFlux                 flux of conservative variable in one time step
        //! @param[in] fFaceFlux                flux of face distribution function in one time step
        void ZeroFaceFlux(RDouble *faceFlux, RDouble3D &fFaceFlux);

        //! Get face tanget vectors by face normal vectors.
        //! It can be used to inverse global coordinates to local coordinates.
        //! @param[in] normalVector             3D face normal vectors
        //! @param[out] tangetVector1           3D face tanget vectors 
        //! @param[out] tangetVector2           3D face tanget vectors
        void CalculateCoordinateTranformationMatrix(RDouble *normalVector, RDouble *tangetVector1, RDouble *tangetVector2);

        //! Calculate tanget vector.
        //! @param[in] dim                      dimension
        //! @param[in] normalVector             3D face normal vectors
        //! @param[out] tangetVector1           3D face tanget vectors
        void CalculateTangetVector(int dim, RDouble *normalVector, RDouble *tangetVector1);

        //! Calculate tanget vector.
        //! @param[in] controlParameters        parameters of UGKSSolver 
        //! @param[in] normalVector             3D face normal vectors
        //! @param[in] tangetVector1            3D face tanget vectors 
        //! @param[in] tangetVector2            3D face tanget vectors
        //! @param[out] discreteVelocity        discrete velocity 
        //! @param[out] wxyz                    minimum discrete velocity of v
        void ComputeDiscreteVariable(Param_UGKSSolver *controlParameters,  RDouble *normalVector, RDouble *tangetVector1, 
                                     RDouble *tangetVector2, RDouble4D &discreteVelocity, RDouble3D &wxyz);

        //! Reconstruction.
        //! @param[in] fLine                    f distribution function in line 
        //! @param[in] dxll                     dx of the left cell of the left cell   
        //! @param[in] dxl                      dx of the left cell  
        //! @param[in] dxr                      dx of the right cell   
        //! @param[in] dxrr                     dx of the right cell of the right cell
        //! @param[in] i                        index 
        //! @param[out] gfl                     face's left gradient 
        //! @param[out] gfr                     face's right gradient
        //! @param[out] ffl                     face's left distribution function f
        //! @param[out] ffr                     face's right distribution function f 
        void Reconstruction(RDouble4D &fLine, RDouble3D &gfl, RDouble3D &gfr, RDouble3D &ffl, 
            RDouble3D &ffr, RDouble dxll, RDouble dxl, RDouble dxr, RDouble dxrr, int i);
        
        //! Flux of kinetic wall boundary condition.
        //! @param[in] fluxParaProxy             proxy of flux parameters
        //! @param[in] controlParameters         parameters of UGKSSolver   
        //! @param[in] uBInf                     variables of wall
        //! @param[in] ff                        face's left/right distribution function f 
        //! @param[in] gf                        face's left/right gradient 
        //! @param[in] wb                        conservative variable of boundary face 
        //! @param[in] wc                        conservative variable of boundary cell
        //! @param[in] dx                        dx of boundary cell
        //! @param[in] direction                 direction of boundary face
        //! -# 0:  left boundary face
        //! -# 1:  right boundary face
        //! @param[out] cf                       frictional drag  
        //! @param[out] qWall                    conservative variable  of wall  
        void FluxBGKUnifiedKineticSlipBC2D(FluxParaProxy * fluxParaProxy, Param_UGKSSolver * controlParameters, 
            RDoubleDBase &uBInf, RDouble3D &ff, RDouble3D &gf, RDouble *wb, RDouble *wc,  RDouble dx, int direction, RDouble *cf, RDouble *qWall);

        //! Judge if calculate gas distribution of boundary cell.
        //! @param[in] direction                 direction of boundary face
        //! -# 0:       left boundary face
        //! -# 1:       right boundary face
        //! @param[out] un                       discrete normal velocity u
        bool JudgeIfCalculate(int direction, RDouble un);

        //! Judge if calculate gas distribution of boundary cell except zero.
        //! @param[in] direction                 direction of boundary face
        //! -# 0:   left boundary face
        //! -# 1:   right boundary face
        //! @param[out] un                       discrete normal velocity u
        bool JudgeIfCalculateExceptZero(int direction, RDouble un);

        //! Compute lambda. 
        //! @param[in] cvs                 conservative variable
        //! @retval     (return value) lamdba
        RDouble ComputeEIG(RDouble * cvs);

        //! Compute collision time. 
        //! @param[in] density             density
        //! @param[in] eig                 represent m/2kT,m is molecular mass,k is the Boltzmann constant, T is the temperature
        //! @param[in] gamma               fluid constant
        //! @param[in] moo                 mach number of inflow
        //! @param[in] refReNumber         reynolds
        //! @retval     (return value) collision time
        RDouble GetCollisionTime(RDouble density,  RDouble eig, RDouble gamma, RDouble moo, RDouble refReNumber);

        //! Compute coefficient related to dt when calculating flux.
        //! @param[in] tau                 collision time
        //! @param[out] rt                  store coefficient of dt
        void GetFluxDtCoefficient(RDouble tau, RDouble * rt);

        //! Solve linearized equation when calculating flux.
        //! @param[in] eig                 represent m/2kT, m is molecular mass, k is the Boltzmann constant, T is the temperature
        //! @param[in] u                   velocity u
        //! @param[in] v                   velocity v
        //! @param[in] w                   velocity w
        //! @param[in] g                   array of right side
        //! @param[in] coefficent          coefficient related to dt when calculating flux
        void SolveEquation(RDouble eig, RDouble u, RDouble v, RDouble w,RDouble *g, RDouble *coefficent);
        
        //! Flux of cell face except solid boundary.
        //! @param[in] fluxParaProxy             proxy of flux parameters
        //! @param[in] controlParameters         parameters of UGKSSolver   
        void FluxBGKUnified3D(FluxParaProxy * fluxParaProxy, Param_UGKSSolver * controlParameters);

        //! Calculate interface conservation.
        //! @param[in] controlParameters         parameters of UGKSSolver   
        //! @param[in] fluxParaProxy             proxy of flux parameters   
        //! @param[in] wn0                       conservative variables of face  
        //! @param[in] eig0                      eig of face  
        void InderfaceVar(Param_UGKSSolver * controlParameters, FluxParaProxy * fluxParaProxy, RDouble *wn0, RDouble &eig0);

        //! Heviside function.
        //! @param[out] un                       discrete normal velocity u
        //! @retval     (returan value) 0 un < 0
        //! @retval     (returan value) 1 un >= 0
        RDouble HevisideFunction(RDouble un);

        //! Calculate flux of equilibrium term.
        //! @param[in] controlParameters         parameters of UGKSSolver   
        //! @param[in] fluxParaProxy             proxy of flux parameters   
        //! @param[in] prim0                     primitive variable of face 
        //! @param[in] eig0                      eig of face  
        //! @param[in] rt                        store coefficient of dt
        //! @param[out] faceFlux                 flux of conservative variable in one time step
        //! @param[out] fFaceFlux                flux of face distribution function in one time step
        void FluxEquilibriumTerm(Param_UGKSSolver * controlParameters,  FluxParaProxy * fluxParaProxy, RDouble * prim0, RDouble * wn0, 
                                 RDouble &eig0, RDouble * rt, RDouble3D &fFaceFlux, RDouble * faceFlux);
        
        //! Calculate flux of non-equilibrium term.
        //! @param[in] controlParameters         parameters of UGKSSolver   
        //! @param[in] fluxParaProxy             proxy of flux parameters   
        //! @param[in] rt                        store coefficient of dt
        //! @param[out] faceFlux                 flux of conservative variable in one time step
        //! @param[out] fFaceFlux                flux of face distribution function in one time step
        void FluxNonEquilibriumTerm(Param_UGKSSolver * controlParameters, FluxParaProxy * fluxParaProxy, 
                                    RDouble * rt, RDouble3D &fFaceFlux, RDouble * faceFlux);
        
        //! Calculate flux of shakhov modified term.
        //! @param[in] controlParameters         parameters of UGKSSolver   
        //! @param[in] fluxParaProxy             proxy of flux parameters   
        //! @param[in] prim0                     primitive variable of face 
        //! @param[in] eig0                      eig of face  
        //! @param[in] rt                        store coefficient of dt
        //! @param[out] faceFlux                 flux of conservative variable in one time step
        //! @param[out] fFaceFlux                flux of face distribution function in one time step
        void FluxShakhov(Param_UGKSSolver * controlParameters, FluxParaProxy * fluxParaProxy, RDouble * prim0, 
                          RDouble &eig0, RDouble * rt, RDouble3D &fFaceFlux, RDouble * faceFlux);

        //! Update fluid variables of UGKS solver.
        //! @param[in] grid                     store information of grid
        void UpdateUGKS(StructGrid *grid);

        //! Update distribution function f of UGKS solver.
        //! @param[in] i                        index of x direction
        //! @param[in] j                        index of y direction
        //! @param[in] k                        index of z direction
        //! @param[in] consNew                  conservation after updating
        void UpdateMicroF(int i, int j, int k, RDouble *consNew);

        //! Update macro fluid variables of UGKS solver by distribution function f.
        //! @param[in] i                        index of x direction
        //! @param[in] j                        index of y direction
        //! @param[in] k                        index of z direction
        void UpdateMicroVariable3D(int iIndex, int jIndex, int kIndex);

        //! Compute post-processing visual variables.
        //! @param[in] postVisualization        store post visualization variables data
        void ComputePostVisualVariables(Post_Visual * postVisualization);

        //! compute Mach number.
        //! @param[in] grid_in                   store information of grid
        RDouble4D * ComputeMachNumber(Grid *grid_in);

        //! compute Global Kn number.
        //! @param[in] grid_in                   store information of grid
        RDouble4D * CompKnGlobal(Grid *grid_in);

        //! Compute air force coefficient.
        //! @param[in] actkey                   action key variable
        void AirForceCoef(ActionKey *actkey);

        //! Compute part of air force coefficient.
        //! @param[in] actkey                   action key variable
        void AirForceCoefParts(ActionKey *actkey, int partID);

        //! Compute residual.
        //! @param[in] actkey                   action key variable
        void GetResidual(ActionKey * actkey);

        //! Init residual.
        //! @param[in] grid_in                  store information of grid
        //! @param[in] res                      residual array
        void InitResidual(StructGrid *grid,RDouble4D &res);

        //! Read restart file.
        //! @param[in] actkey                   action key variable
        void ReadRestart(ActionKey *actkey);

        //! Write restart file.
        //! @param[in] actkey                   action key variable
        void DumpRestart(ActionKey *actkey);

        //! Read restart file of NS.
        //! @param[in] actkey                   action key variable
        void InitFlowFromNS(ActionKey *actkey);

        //int GetNumberOfEquations(); ///4为二维，5为三维

    private:
        //! Store post-processing visual variables.
        //! @param[in] postVisualization_in        store post visualization variables data
        void SetPostVisualization(Post_Visual * postVisualization_in){ this ->postVisualization = postVisualization; }

        //! Get grid information.
        //! @param[in] level        multigrid level
        Grid * GetGrid(int level = 0) { return this->GetGeometry()->GetGrid(level); }

        //! Store constrol parameters.
        Param_UGKSSolver * GetControlParameters(){ return this->controlParameters; }
    
        int GetNumberOfEquations() {return this->controlParameters->GetNumberOfEquations();}
}; 
}