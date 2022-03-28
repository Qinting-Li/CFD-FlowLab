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
//! @brief     It defines parameters when calculating flux.
//! @author    Wang Pei, original by Jiang Dingwu.

#pragma once
#include "TypeDefine.h"

namespace PHSPACE
{
    //! @brief Param_UGKSSolver defines parameters when calculating flux.
    class FluxParaProxy
    {
    public:
        FluxParaProxy();
        ~FluxParaProxy();
    private:
        int neqn;
        //RDouble * gamal;
       // RDouble * gamar;
        RDouble * wlc;
        RDouble * wrc;
        RDouble3D * gfl;
        RDouble3D * gfr;
        RDouble3D * ffl;
        RDouble3D * ffr;
        RDouble  dxl;
        RDouble  dxr;
        RDouble * faceFlux;
        RDouble * varWall;
        RDouble4D * fFluxLine;
        RDouble2D * fluxLine;
        RDouble3D * fFaceFlux;
        RDouble4D * fLine;
        RDouble2D * qLine;
        RDouble4D * discreteVelocity;
        RDouble3D * wxyz;

        RDouble1D * xfnc;
        RDouble1D * yfnc;
        RDouble1D * zfnc;
        RDouble1D * areac;
        RDouble1D * thicknessLine;
        RDouble * aInverse;
        RDouble * bInverse;
        RDouble * cInverse;
    public:
        void CreateFluxParaProxyLine(Range IS, Range NU,Range NV,Range NW,Range M);
        void CreateFluxParaProxyFace(Range NU,Range NV,Range NW,Range M);

        void SetNumberofTotalEquation(const int &neqn) { this->neqn = neqn; }
        //void SetLeftAndRightGama(RDouble * gamal, RDouble * gamar) { this->gamal = gamal; this->gamar = gamar; }
        void SetLeftAndRightW(RDouble * wlc, RDouble * wrc) { this->wlc = wlc; this->wrc = wrc; }
        void SetLeftAndRightGrad(RDouble3D * gfl, RDouble3D * gfr) { this->gfl = gfl; this->gfr = gfr; }
        void SetLeftAndRightFDistribution(RDouble3D * ffl, RDouble3D * ffr) { this->ffl = ffl; this->ffr = ffr; }
        void SetLeftAndRightDs(RDouble  dxl, RDouble  dxr) { this->dxl = dxl; this->dxr = dxr; }

        void SetFFaceFlux(RDouble3D *fFaceFlux){this->fFaceFlux = fFaceFlux;}
        void SetFaceFlux(RDouble * faceFlux) {this->faceFlux = faceFlux;}
        void SetVarWall(RDouble * varWall) {this->varWall = varWall;}
        void SetFaceNormal(RDouble1D * xfn, RDouble1D *yfn, RDouble1D *zfn) { this->xfnc = xfn; this->yfnc = yfn; this->zfnc = zfn; }
        void SetThicknessLine(RDouble1D * thicknessLine) {this->thicknessLine = thicknessLine;}
        //void SetFLineBoundCell(RDouble3D *fLineBoundCell){ this->fLineBoundCell = fLineBoundCell; }
        void SetDiscreteVelocity(RDouble4D *discreteVelocity){ this->discreteVelocity = discreteVelocity; }
        void SetWxyz(RDouble3D * wxyz){ this->wxyz = wxyz; }

        void SetFaceArea(RDouble1D * areac) {this->areac = areac; }
        void SetFaceConservation(RDouble2D * qLine) {this->qLine = qLine; }
        void SetFaceF(RDouble4D * fLine) {this->fLine = fLine; }

        int GetNumberofTotalEquation() { return this->neqn; }

        //RDouble * GetLeftGama() { return this->gamal; }
        //RDouble * GetRightGama() { return this->gamar; }
        RDouble * GetLeftW() { return this->wlc; }
        RDouble * GetRightW() { return this->wrc; }
        RDouble3D * GetLeftGrad() { return this->gfl; }
        RDouble3D * GetRightGrad() { return this->gfr; }
        RDouble3D * GetLeftFDistribution() { return this->ffl; }
        RDouble3D * GetRightFDistribution() { return this->ffr; }
        RDouble  GetLeftDs() { return this->dxl; }
        RDouble  GetRightDs() { return this->dxr; }

        RDouble4D *GetFFluxLine() { return this->fFluxLine;}
        RDouble2D *GetFluxLine() { return this->fluxLine;}
        RDouble3D *GetFFaceFlux() { return this->fFaceFlux;}
        RDouble *GetFaceFlux() { return this->faceFlux;}
        RDouble *GetVarWall() { return this->varWall;}

        RDouble1D * GetFaceNormalX() { return this->xfnc; }
        RDouble1D * GetFaceNormalY() { return this->yfnc; }
        RDouble1D * GetFaceNormalZ() { return this->zfnc; }
        RDouble1D * GetThicknessLine() { return this->thicknessLine;}
        RDouble4D * GetDiscreteVelocity(){ return this->discreteVelocity; }
        RDouble3D * GetWxyz(){ return this->wxyz; }
        RDouble2D * GetFaceConservation() { return this->qLine; }
        RDouble4D * GetFaceF() { return this->fLine; }
        RDouble1D * GetFaceArea() { return this->areac; }
        RDouble * GetTransformCoeffUNormal() { return this->aInverse; }
        RDouble * GetTransformCoeffVTangential() { return this->bInverse; }
        RDouble * GetTransformCoeffWTangential() { return this->cInverse; }
    };
}