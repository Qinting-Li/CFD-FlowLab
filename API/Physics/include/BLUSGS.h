#include "Geo_UnstructGrid.h"
#include "Geo_StructGrid.h"
#include "FieldProxy.h"
#include "FaceProxy.h"
#include "PHMatrix.h"

namespace PHSPACE
{
// this is for Block Lower-Upper Symmetric Gauss¨CSeidel Scheme according wzj reference in 2000 year;
//! forward sweep in BLUSGS time methods.
//! @param[in] grid    computational mesh.
//! @param[in] dqProxy        a proxy to store dq.
//! @param[in] LUplusDQ       a proxy to store dq.
//! @param[in] sweepNormal    a value depends the number of inner iterations.
void SolveBLUSGSForwardInHouse(UnstructGrid *grid, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal);

//! backward sweep in BLUSGS time methods.
void SolveBLUSGSBackwardInHouse(UnstructGrid *grid, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal);

//! To execute the forward sweep of the block LU-SGS method.
//! @param[in ]: gridIn denotes the current computational regions of the grids.
//! @param[out]: dqProxy is an array that stores the difference values between the conservative variables at n and n+1 time step.
void SolveBLUSGSForwardInHouse(StructGrid *gridIn, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal);

//! To execute the backward sweep of the block LU-SGS method.
//! @param[in ]: gridIn denotes the current computational regions of the grids.
//! @param[out]: dqProxy is an array that stores the difference values between the conservative variables at n and n+1 time step.
void SolveBLUSGSBackwardInHouse(StructGrid *gridIn, FieldProxy *dqProxy, FieldProxy *LUplusDQ, RDouble &sweepNormal);

//! it computes (dF/dQ)*dQ in cell j. for the upper and lower sweeps.
void FluxDiffFace3D(RDouble *q, RDouble *dq, RDouble *nxyz, RDouble gama, RDouble *df, int nm);

//! it computes |Cij|*dQ in cell j. for the upper and lower sweeps for consitive variables.
void RoeDisDiff3D(RDouble *ql, RDouble *qr, RDouble *dq_loc, RDouble *nxyz, RDouble gama, 
                  int RoeEntropyFixMethod, RDouble RoeEntropyFixCoef1, RDouble RoeEntropyFixCoef2, RDouble *flux, int nm);

//! it computes Jv*dQ in cell j. for the upper and lower sweeps for consitive variables.
//! Jv is viscous flux jacobi matrix for consitive variables.
void VisFluxDiffOneFace3D(RDouble *q, RDouble *dq, RDouble *flux, RDouble vis, RDouble tm, RDouble kcp, RDouble delta, RDouble *nxyz);

//! Notes: This function is same to the SetEntropyFixCoefficients, the two functions should be merged future.
//! Compute the entropy fix coefficients by the inflow mach number.
//! @param[in ]: type, '1' denotes the diagonal matrix. '2' denote the lower/upper sweep.
void GetEntropyFixMethodAndCoefficients(int &RoeEntropyFixMethod, RDouble &RoeEntropyFixCoef1, RDouble &RoeEntropyFixCoef2);

void SetGhostDQLUSGSInHouse(UnstructGrid *gridIn, RDouble **dq, RDouble **q);

//! Compute and save the diagonal matrice on grids that are used for block LU-SGS method.
//! @param[in ]: gridIn denotes the current computational regions of the grids.
void ComputeDiagonalMatrice(StructGrid *gridIn);

//! To compute the spectrum radius of the inviscid flux Jacobian matrix in current control volume.
//! @param[in ]: primitiveVars is the array that stores the primitive variables of the current cell.
//! @param[in ]: xfn, yfn and zfn are components of the average normal vector located in the center of the current cell.
//! @param[in ]: area denotes the avarage area located in the center of the current cell.
//! @param[in ]: vgn denotes the moment velocity located in the center of the current cell, and it is used for dynamic mesh.
//! @param[in ]: gama denotes the specific heat ratio.
//! @param[out]: return the spectrum radius of the inviscid flux Jacobian matrix computed by the variables of the current cell.
RDouble ObtainInviscidSpectrumRadius(RDouble *primitiveVars, RDouble xfn, RDouble yfn, RDouble zfn, RDouble area, RDouble vgn, RDouble gama);

//! To compute the spectrum radius of the viscous flux Jacobian matrix in current control volume.
//! @param[in ]: primitiveVars is the array that stores the primitive variables of the current cell.
//! @param[in ]: area denotes the avarage area located in the center of the current cell.
//! @param[in ]: vol denotes the avarage volume located in the center of the current cell.
//! @param[in ]: gama denotes the specific heat ratio.
//! @param[in ]: laminarViscosity and turbulentViscosity denote the laminar viscosity and turbulent viscosity, respectively.
//! @param[out]: return the spectrum radius of the viscous flux Jacobian matrix computed by the variables of the current cell.
RDouble ObtainViscousSpectrumRadius(RDouble *primitiveVars, RDouble area, RDouble vol, RDouble gama, RDouble laminarViscosity, RDouble turbulentViscosity);

//! Construct the Jacobian matrix of the inviscid flux item.
//! @param[in ]: primitiveVars is the array that stores the primitive variables of the current cell.
//! @param[in ]: temperature is the array that stores the temperature values of the current cell.
//! @param[in ]: xfn, yfn and zfn are components of the average normal vector located in the center of the current cell.
//! @param[in ]: area denotes the avarage area located in the center of the current cell.
//! @param[in ]: vgn denotes the moment velocity located in the center of the current cell, and it is used for dynamic mesh.
//! @param[in ]: gama denotes the specific heat ratio.
//! @param[out]: tempMatrix is an matrix that stores the Jacobian matrix of the inviscid flux item.
void ObtainInviscidMatrix(RDouble *primitiveVars, RDouble *temperature, RDouble xfn, RDouble yfn, RDouble zfn, RDouble area, RDouble vgn, 
                          RDouble gama, RDouble2D *tempMatrix);

//! Construct the Jacobian matrix of the viscous flux item.
//! @param[in ]: primitiveVars is the array that stores the primitive variables of the current cell.
//! @param[in ]: temperature is the array that stores the temperature values of the current cell.
//! @param[in ]: xfn, yfn and zfn are components of the average normal vector located in the center of the current cell.
//! @param[in ]: area denotes the avarage area located in the center of the current cell.
//! @param[in ]: vol denotes the avarage volume located in the center of the current cell.
//! @param[in ]: vgn denotes the moment velocity located in the center of the current cell, and it is used for dynamic mesh.
//! @param[in ]: gama denotes the specific heat ratio.
//! @param[in ]: laminarViscosity and turbulentViscosity denote the laminar viscosity and turbulent viscosity, respectively.
//! @param[out]: tempMatrix is an matrix that stores the Jacobian matrix of the viscous flux item.
void ObtainViscousMatrix(RDouble *primitiveVars, RDouble *temperature, RDouble xfn, RDouble yfn, RDouble zfn, RDouble area, RDouble vol, 
                         RDouble vgn, RDouble gama, RDouble laminarViscosity, RDouble turbulentViscosity, RDouble2D *tempMatrix);

void FillField(FieldProxy *field1Proxy, FieldProxy *field2Proxy);

FieldProxy * GetFieldProxy(StructGrid *gridIn, const string &field_name);

//! Compute the Matrixs of left hand side for time integration on unstruct grid:\n
//! @param[in] gridIn    computational mesh.
void LHSMatrixsInHouse(UnstructGrid *gridIn);

// create the array for LHS Matrixs and assign to zero.
//! @param[in] gridIn    computational mesh.
void InitLHSMatrix(UnstructGrid *gridIn);

// add the diagonal part to the LHS Matrixs.
//! @param[in] gridIn    computational mesh.
void LHSDiagMatrix(UnstructGrid *gridIn);

//! calculate the inviscid matrixs for LHS matrixs in BLUSGS
void InviscidMatrixs(UnstructGrid *gridIn);

//! calculate the viscous matrixs for LHS matrixs in BLUSGS.
void ViscousMatrixs (UnstructGrid *gridIn);

//! calculate the viscous matrix from face localStart to localEnd in BLUSGS.
void ViscousMatrix(UnstructGrid *gridIn, FaceProxy *faceProxy, MAT3D< RDouble > *dfdq, int localStart, int localEnd);

//! calculate the roe matrix from face localStart to localEnd in BLUSGS.
void RoeDispMatrix(UnstructGrid *gridIn, FaceProxy *faceProxy, MAT3D< RDouble > *dfdq, int localStart, int localEnd);

//! calculate Jacobian matrix from face localStart to localEnd BLUSGS.
void GetJacobianFull(UnstructGrid *gridIn, RDouble **q, MAT3D< RDouble > *dfdq, int localStart, int localEnd);

//! add dfdq to the left cell of the face .
void AddImpLHSLeft(UnstructGrid *gridIn, MAT3D< RDouble > *dfdq, int localStart, int localEnd);

//! add dfdq to the right cell of the face .
void AddImpLHSRight(UnstructGrid *gridIn, MAT3D< RDouble > *dfdq, int localStart, int localEnd);

FaceProxy * CreateFaceProxyInHouse(UnstructGrid *gridIn);

void GetQlQrInHouse(UnstructGrid *gridIn, FaceProxy *faceProxy, int localStart, int localEnd);

}