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
//! @file      DSMCZone_related_vardef.h
//! @brief     Zone_related var definition for DSMC.
//! @author    Jia Chuan.
#pragma once

namespace PHSPACE
{

	using namespace std;
	class Zone_related 
{
	private:
		int PntNum,SurfNum,EleNum,BCSurfNum,MoleNum,MoleStart,MoleNum2;          
		int SubEleNum,EleStart,EleEnd;            
		int SurfNum_Wall,SurfNum_In,SurfNum_Out,SurfNum_Sym,SurfNum_Inner;            
		int *nSurfBndryMap,*nSurfWallMap,*nSurfFreeMap,*nSurfVacuMap,*nSurfSymmMap,*nSurfInneMap;     
		int *nSurfBndryControl,*nOutVect_FR,*nOutVect_SM;
		int *ISC,*SurfPntNum;  
		int **ElePnt,**nSurfVolMap,**EleSurf,**SurfPnt,**nRatioSub;           
		int **MoleMap;   
		double *fX,*fY,*fZ,*EleWeight,*EleVol,*SurfArea,*BME,*BMR;  
		double **EleCenter,**SurfCenter,**SurfNormal;  
		double ***SurfTan_Wall,***SurfTan_In; 
	public:
		Zone_related();
		~Zone_related();


	public:
		void setPntNum(int PntNum){this->PntNum=PntNum;}; 
		int getPntNum(){return this->PntNum;};

		void setSurfNum(int SurfNum){this->SurfNum=SurfNum;}; 
		int getSurfNum(){return this->SurfNum;};

		void setEleNum(int EleNum){this->EleNum=EleNum;}; 
		int getEleNum(){return this->EleNum;};

		void setBCSurfNum(int BCSurfNum){this->BCSurfNum=BCSurfNum;}; 
		int getBCSurfNum(){return this->BCSurfNum;};

		void setMoleNum(int MoleNum){this->MoleNum=MoleNum;}; 
		int getMoleNum(){return this->EleNum;};

		void setMoleStart(int MoleStart){this->MoleNum=MoleNum;}; 
		int getMoleStart(){return this->MoleStart;};

		void setMoleNum2(int MoleNum2){this->MoleNum2=MoleNum2;}; 
		int getMoleNum2(){return this->MoleNum2;};

		void setSubEleNum(int SubEleNum){this->SubEleNum=SubEleNum;}; 
		int getSubEleNum(){return this->SubEleNum;};

		void setEleStart(int EleStart){this->EleStart=EleStart;}; 
		int getEleStart(){return this->EleStart;};

		void setEleEnd(int EleEnd){this->EleEnd=EleEnd;}; 
		int getEleEnd(){return this->EleEnd;};

		void setSurfNum_Wall(int SurfNum_Wall){this->SurfNum_Wall=SurfNum_Wall;}; 
		int getSurfNum_Wall(){return this->SurfNum_Wall;};

		void setSurfNum_In(int SurfNum_In){this->SurfNum_In=SurfNum_In;}; 
		int getSurfNum_In(){return this->SurfNum_In;};

		void setSurfNum_Out(int SurfNum_Out){this->SurfNum_Out=SurfNum_Out;}; 
		int getSurfNum_Out(){return this->SurfNum_Out;};

		void setSurfNum_Sym(int SurfNum_Sym){this->SurfNum_Sym=SurfNum_Sym;}; 
		int getSurfNum_Sym(){return this->SurfNum_Sym;};

		void setSurfNum_Inner(int SurfNum_Inner){this->SurfNum_Inner=SurfNum_Inner;}; 
		int getSurfNum_Inner(){return this->SurfNum_Inner;};

	public:

		void setnSurfBndryMap(int *nSurfBndryMap){this->nSurfBndryMap=nSurfBndryMap;}; 
		int *getnSurfBndryMap(){return this->nSurfBndryMap;};

		void setnSurfWallMap(int *nSurfWallMap){this->nSurfWallMap=nSurfWallMap;}; 
		int *getnSurfWallMap(){return this->nSurfWallMap;};

		void setnSurfFreeMap(int *nSurfFreeMap){this->nSurfFreeMap=nSurfFreeMap;}; 
		int *getnSurfFreeMap(){return this->nSurfFreeMap;};

		void setnSurfVacuMap(int *nSurfVacuMap){this->nSurfVacuMap=nSurfVacuMap;}; 
		int *getnSurfVacuMap(){return this->nSurfVacuMap;};

		void setnSurfSymmMap(int *nSurfSymmMap){this->nSurfSymmMap=nSurfSymmMap;}; 
		int *getnSurfSymmMap(){return this->nSurfSymmMap;};

		void setnSurfInneMap(int *nSurfInneMap){this->nSurfInneMap=nSurfInneMap;}; 
		int *getnSurfInneMap(){return this->nSurfInneMap;};

		void setnSurfBndryControl(int *nSurfBndryControl){this->nSurfBndryControl=nSurfBndryControl;}; 
		int *getnSurfBndryControl(){return this->nSurfBndryControl;};

		void setnOutVect_FR(int *nOutVect_FR){this->nOutVect_FR=nOutVect_FR;}; 
		int *getnOutVect_FR(){return this->nOutVect_FR;};

		void setnOutVect_SM(int *nOutVect_SM){this->nOutVect_SM=nOutVect_SM;}; 
		int *getnOutVect_SM(){return this->nOutVect_SM;};

		void setISC(int *ISC){this->ISC=ISC;}; 
		int *getISC(){return this->ISC;};

		void setSurfPntNum(int *SurfPntNum){this->SurfPntNum=SurfPntNum;}; 
		int *getSurfPntNum(){return this->SurfPntNum;};

	public:   

		void SetElePnt(int n, int value);
		int ** GetElePnt(){return ElePnt;};

		void SetnSurfVolMap(int n, int value);
		int ** GetnSurfVolMap(){return nSurfVolMap;};

		void SetEleSurf(int n, int value);
		int ** GetEleSurf(){return EleSurf;};

		void SetSurfPnt(int n, int value);
		int ** GetSurfPnt(){return SurfPnt;};

		void SetnRatioSub(int n, int value);
		int ** GetnRatioSub(){return nRatioSub;};

		void SetMoleMap(int n, int value);
		int ** GetMoleMap(){return MoleMap;};


	public:

		void setfX(double *fX){this->fX=fX;}; 
		double *getfX(){return this->fX;};

		void setfY(double *fY){this->fY=fY;}; 
		double *getfY(){return this->fY;};
	
		void setfZ(double *fZ){this->fZ=fZ;}; 
		double *getfZ(){return this->fZ;};

		void setEleWeight(double *EleWeight){this->EleWeight=EleWeight;}; 
		double *getEleWeight(){return this->EleWeight;};

		void setEleVol(double *EleVol){this->EleVol=EleVol;}; 
		double *getEleVol(){return this->EleVol;};

		void setSurfArea(double *SurfArea){this->SurfArea=SurfArea;}; 
		double *getSurfArea(){return this->SurfArea;};

		void setBME(double *BME){this->BME=BME;}; 
		double *getBME(){return this->BME;};

		void setBMR(double *BMR){this->BMR=BMR;}; 
		double *getBMR(){return this->BMR;};

	public:

		void SetEleCenter(double n, double value);
		double ** GetEleCenter(){return EleCenter;};

		void SetSurfCenter(double n, double value);
		double ** GetSurfCenter(){return SurfCenter;};

		void SetSurfNormal(double n, double value);
		double ** GetSurfNormal(){return SurfNormal;};

		void SetSurfTan_Wall(double n, double value,double value_1);
		double *** GetSurfTan_Wall(){return SurfTan_Wall;};

		void SetSurfTan_In(double n, double value);
		double *** GetSurfTan_In(){return SurfTan_In;};


};
}