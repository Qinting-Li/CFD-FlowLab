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
//! @file      Jacobian.h
//! @brief     Explain this file briefly.
//! @author    He Xin, Bell, He Kun, He Xianyao, Liu Jian, Xu Qingxin, etc.

#pragma once
#include "Precision.h"
namespace PHSPACE
{
//!The fuction computes the formula (M +- R) * dQ, M and R denote the Jacobian matrix and spectrum radius of inviscid flux,respectively,\n
//!dQ is the increment of conservative vector.
void MXDQ_STD(RDouble *prim, RDouble nxs, RDouble nys, RDouble nzs, RDouble area, RDouble vgn, RDouble gama, RDouble *dq, RDouble *f, 
              int nm, int nl, RDouble radius, int ipn);

void MXDQ_STD(RDouble *prim, RDouble *gykb, RDouble *dq, RDouble *df, RDouble gama, int nm, int nl, RDouble radius, int ipn);

void ComputePreconMatrix(int sign, RDouble **GG, RDouble *prim, RDouble tm, RDouble gama);
void PreconMatrixDF(RDouble *df, RDouble **GG, int nm);
void PreconMatrixRes(RDouble **dq, RDouble **GG, int nm, int icell);

//! The fuction computes the formula (M +- R) * dQ, where M and R denote the Jacobian matrix and spectrum radius of inviscid flux,respectively.
//! @param[in ]: primitiveVars is the array that stores the primitive variables of the current cell.
//! @param[in ]: temperature is the array that stores the temperature values of the current cell.
//! @param[in ]: nx, ny and nz are components of the average normal vector located in the center of the current cell.
//! @param[in ]: faceArea denotes the avarage area located in the center of the current cell.
//! @param[in ]: moveVelocity denotes the movement velocity located in the center of the current cell, and it is used for dynamic mesh.
//! @param[in ]: deltaQ is the increment of conservative vector.
//! @param[out]: flux is the inviscid flux that equals to the Jacobian matrix multiplies the difference values of the conservative variables.
//! @param[in ]: radius denotes the spectrum radius of invsicid Jacobian matrix.
//! @param[in ]: iSign denotes the computation symbol, -1 denotes M substracts R, 1 denotes M plus R.
void ChemicalMXDQ(RDouble *primitiveVars, RDouble *temperature, RDouble nx, RDouble ny, RDouble nz, RDouble faceArea, RDouble moveVelocity,
                  RDouble *deltaQ, RDouble *flux, int nNSEquation, int nLaminar, int nTemperatureModel, RDouble radius, int iSign,int nElectronIndex);
void ChemicalMXDQR(RDouble *primitiveVars, RDouble *temperature, RDouble nx, RDouble ny, RDouble nz, RDouble faceArea, RDouble moveVelocity,
                  RDouble *deltaQ, RDouble *flux, int nNSEquation, int nLaminar, int nTemperatureModel, RDouble radius, int iSign, int nElectronIndex,
                  RDouble *speciesCvs, RDouble *speciesEtrs, RDouble totalH,  RDouble totalCv, RDouble totalCvtr, RDouble totalCvv, RDouble totalCve);

//!The fuction computes the formula (Mv +- Rv) * dQ, Mv and Rv denote the Jacobian matrix and spectrum radius of viscous flux,respectively,\n
//!dQ is the increment of conservative vector.
void MVXDQ_STD_LP(RDouble *prim, RDouble nxs, RDouble nys, RDouble nzs, RDouble area, RDouble vol, RDouble visl, RDouble vist, RDouble gama, 
                  RDouble *dq, RDouble *f, int nm, int nl, RDouble radius, int ipn);

//! this is for unstructgrid solver. 
void PreconMatrixFreeDF(RDouble *q, RDouble *dq, RDouble *nxyz, RDouble **MG, RDouble gama, RDouble *df, int nm);

//! this is for structgrid solver. 
void PreconMatrixFreeDF(RDouble *q, RDouble *dq, RDouble *gykb, RDouble **MG, RDouble tm, RDouble gama, RDouble *df, int nm, int ipn);
void PreconMatrixFreeDF(RDouble *q, RDouble *dq, RDouble nxs, RDouble nys, RDouble nzs, RDouble area, RDouble vgn, RDouble **MG, RDouble gama, RDouble *df, int nm, int ipn);

}
