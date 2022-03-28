#include "Geo_UnstructGrid.h"
#include "PHHeapSort.h"
#include "Constants.h"
#include "PHMatrix.h"
#include "Flux_Inviscid.h"

namespace PHSPACE
{
//! calculate the three Dialog, Lower, Upper Matrixs in all Lines.
void CalDLUMatrixLineLUSGS(Grid* gridIn);

//! get the three Matrixs pointer.
void GetDLUMatrixByLine(Grid* gridIn);

//! calculate the the face is shock face or not.
void CalIsShockFace(Grid* gridIn, int* IsShockFace);

//! calculate viscous jacobian matrix of left hand side.
void CalVisJacobianMatrix(RDouble** matrix, RDouble* xyzfn, RDouble* q, RDouble refGama, RDouble refMachNumber, RDouble miu, RDouble k_mod, RDouble sd1);

//! calculate roe matrix of left hand side.
void CalRoeMatrix(RDouble** matrix, RDouble* qL, RDouble* qR, RDouble* xyzfn, RDouble refGama, RDouble area);

//! calculate conv matrix of left hand side.
void CalConvJacobianMatrix(RDouble** Matrix, RDouble* xzyfn, RDouble* q, RDouble refGama, RDouble area);

//! calculate conv flux of left hand side.
void FluxLUSGS3D(RDouble* flux, RDouble* q, RDouble* DQ, RDouble* xyzfn, RDouble refGama, RDouble lhs_omga);

//! calculate conv flux without eigen of left hand side.
void FluxLUSGS3DNoEig(RDouble* flux, RDouble* q, RDouble* DQ, RDouble* xyzfn, RDouble refGama);

//! calculate eigen value flux of left hand side.
void RoeDeltaFlux(RDouble* flux, RDouble* qL, RDouble* qR, RDouble* DQ, RDouble* xyzfn, RDouble refGama, RDouble RoeEntropyFixCoef1, RDouble RoeEntropyFixCoef2);

//! calculate vicous flux of left hand side.
void ViscousDeltaFlux(RDouble* flux, RDouble* q, RDouble* dQ, RDouble* xyzfn, RDouble refGama,
    RDouble vis, RDouble kcp, RDouble T);

//! Forward sweep in line implicit LUSGS time methods.
void SolveLineLUSGSForwardInHouse(Grid* grid, FieldProxy* dqProxy, FieldProxy* LUplusDQ, RDouble& sweepNormal);

//! Backward sweep in line implicit LUSGS time methods.
void SolveLineLUSGSBackwardInHouse(Grid* grid, FieldProxy* dqProxy, FieldProxy* LUplusDQ, RDouble& sweepNormal);

RDouble GetRadInHouse(RDouble *prim, RDouble *nxyz, RDouble gama);

}