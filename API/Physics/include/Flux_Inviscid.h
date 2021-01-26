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
//! @file      Flux_Inviscid.h
//! @brief     Inviscid (convection) flux computing, including number of inviscid scheme,
//!            such as Roe, Steger-Warmming, Van Leer, HLLC, AUSM, etc.
//! @author    Bell, He Xin.
#pragma once
#include "FaceProxy.h"
namespace PHSPACE
{
class InviscidSchemeParameter;

typedef void (*INVScheme) (FaceProxy *face_proxy, InviscidSchemeParameter *invSchemePara);

int GetSchemeID(const string &scheme_name);

INVScheme GetINVScheme(int scheme_id);

//! This Roe scheme is coding by the conservative variable way, in which the DF is expressed by dH.
void Roe_Scheme_ConservativeForm(FaceProxy *face_proxy, InviscidSchemeParameter *invSchemePara);

//! This Roe scheme is coding by the primitive variable way, in which the DF is split into three parts.
void Roe_Scheme_PrimitiveForm(FaceProxy *face_proxy, InviscidSchemeParameter *invSchemePara);

void Steger_Scheme (FaceProxy *face_proxy, InviscidSchemeParameter *invSchemePara);

void Ausmdv_Scheme (FaceProxy *face_proxy, InviscidSchemeParameter *invSchemePara);

void Ausmpw_Scheme (FaceProxy *face_proxy, InviscidSchemeParameter *invSchemePara);

void Hlle_Scheme   (FaceProxy *face_proxy, InviscidSchemeParameter *invSchemePara);

void Vanleer_Scheme(FaceProxy *face_proxy, InviscidSchemeParameter *invSchemePara);

RDouble FMSplit1(const RDouble &mach, const RDouble &ipn);
RDouble FMSplit2(const RDouble &mach, const RDouble &ipn);
RDouble FMSplit4(const RDouble &mach, const RDouble &beta, const RDouble &ipn);
RDouble FPSplit5(const RDouble &mach, const RDouble &alpha, const RDouble &ipn);

void RoeAveragedVariables(RDouble rhoLeft, RDouble uLeft, RDouble vLeft, RDouble wLeft, RDouble pLeft, RDouble gamaLeft, 
                          RDouble rhoRight, RDouble uRight, RDouble vRight, RDouble wRight, RDouble pRight, RDouble gamaRight,
                          RDouble &rhoRoe, RDouble &uRoe, RDouble &vRoe, RDouble &wRoe, RDouble &pRoe);

//! @brief The parameters that inviscid flux computing need.
//!        These data are compacted into this class to fill the need of OpenMP parallel.
//!        If not, would lead to data completion problems.
class InviscidSchemeParameter
{
public:
    InviscidSchemeParameter();
    ~InviscidSchemeParameter();

private:	
    int nLen;
    int nm;
    int nLaminar;
    int nEquation;
    int nChemical;
    int nPrecondition;
    int nTemperatureModel;
    int nEntrFix;

    RDouble refMachNumber;

    RDouble entrFixCoefLeft;
    RDouble entrFixCoefRight;

    RDouble *gamal; RDouble *gamar;
    RDouble *lrtl; RDouble *lrtr;
    RDouble **ql; RDouble **qr;
    RDouble **tl; RDouble **tr; //left and right temperature.
    RDouble **flux;
    RDouble *xfn; RDouble *yfn; RDouble *zfn;
    RDouble *area;
    RDouble *vgn ;	

public:
    void SetLength(const int nLen) {this->nLen = nLen;}
    void SetNumberOfLaminar(const int nLaminar) {this->nLaminar = nLaminar;}
    void Setnm(const int nm) {this->nm = nm;}
    void SetIfPrecondition(const bool nPrecondition) {this->nPrecondition = nPrecondition;}
    void SetNumberOfChemical(const int nChemical) {this->nChemical = nChemical;}
    void SetNumberOfTemperatureModel(const int nTemperatureModel) {this->nTemperatureModel = nTemperatureModel;}
    void SetNumberOfTotalEquation(const int nEquation) {this->nEquation = nEquation;}	
    void SetEntropyFixMethod(const int nEntrFix) {this->nEntrFix = nEntrFix;}
    void SetEntropyFixCoefficients(const RDouble entrFixCoefLeft, const RDouble entrFixCoefRight)
    {this->entrFixCoefLeft = entrFixCoefLeft; this->entrFixCoefRight = entrFixCoefRight;}
    void SetMachNumber(const RDouble refMachNumber) {this->refMachNumber = refMachNumber;}
    void SetLeftAndRightGama(RDouble *gamal, RDouble *gamar) {this->gamal = gamal; this->gamar = gamar;}
    void SetLeftAndRightPressureCoefficient(RDouble *lrtl, RDouble *lrtr) {this->lrtl = lrtl; this->lrtr = lrtr;}
    void SetLeftAndRightQ(RDouble **ql, RDouble **qr) {this->ql = ql; this->qr = qr;}
    void SetLeftAndRightTemperature(RDouble **tl, RDouble **tr) {this->tl = tl; this->tr = tr;}
    void SetFlux(RDouble **flux) {this->flux = flux;}
    void SetFaceNormal(RDouble *xfn, RDouble *yfn, RDouble *zfn) {this->xfn = xfn; this->yfn = yfn; this->zfn = zfn;}
    void SetFaceArea(RDouble *area) {this->area = area;}
    void SetFaceVelocity(RDouble *vgn) {this->vgn = vgn;}

    int GetLength() {return this->nLen;}
    int Getnm() {return this->nm;}
    int GetNumberOfLaminar() {return this->nLaminar;}
    int GetNumberOfTotalEquation() {return this->nEquation;}
    int GetNumberOfChemical() {return this->nChemical;}
    int GetNumberOfTemperatureModel() {return this->nTemperatureModel;}
    int GetEntropyFixMethod() {return this->nEntrFix;}
    int GetIfPrecondition() {return this->nPrecondition;}
    RDouble GetLeftEntropyFixCoefficient() {return this->entrFixCoefLeft;}
    RDouble GetRightEntropyFixCoefficient() {return this->entrFixCoefRight;}
    RDouble GetMachNumber() {return this->refMachNumber;}
    RDouble *GetLeftGama() {return this->gamal;}
    RDouble *GetRightGama() {return this->gamar;}

    RDouble *GetLeftPressureCoefficient() {return this->lrtl;}
    RDouble *GetRightPressureCoefficient() {return this->lrtr;}

    RDouble **GetLeftQ() {return this->ql;}
    RDouble **GetRightQ() {return this->qr;}

    RDouble** GetLeftTemperature() {return this->tl;}
    RDouble** GetRightTemperature() {return this->tr;}

    RDouble **GetFlux() {return this->flux;}
    RDouble *GetFaceNormalX() {return this->xfn;}
    RDouble *GetFaceNormalY() {return this->yfn;}
    RDouble *GetFaceNormalZ() {return this->zfn;}
    RDouble *GetFaceArea() {return this->area;}
    RDouble *GetFaceVelocity() {return this->vgn;}
};

}
