#include "DSMCMolecule.h"
#include "DSMCSolver.h"
#include "DSMCConsts.h"
#include "GlobalDataBase.h"
#include "Pointer.h"
#include <cmath>
using namespace std;

namespace PHSPACE
{

Molecule::Molecule()
{
    moleculeLocation = new int[MaxMoleculeNumber];
    moleculeSpecies = new int[MaxMoleculeNumber];
    IPCP = new int[MaxMoleculeNumber];

    moleculeVibrationLevel = new int [MaxMoleculeNumber];
    moleculeRotationalEnergy = new RDouble [MaxMoleculeNumber];
    moleculeVibrationalEnergy = new RDouble [MaxMoleculeNumber];
    moleculePosition = NewPointer2<RDouble>(MaxMoleculeNumber,3);
    moleculeVelocity = NewPointer2<RDouble>(MaxMoleculeNumber,3);

}
Molecule::~Molecule()
{
    delete [] moleculeLocation; moleculeLocation = NULL;
    delete [] moleculeSpecies; moleculeSpecies = NULL;
    delete [] IPCP; IPCP = NULL; 
    
    delete [] moleculeVibrationLevel ; moleculeVibrationLevel = NULL;
    
    delete [] moleculeRotationalEnergy; moleculeRotationalEnergy = NULL;
    delete [] moleculeVibrationalEnergy; moleculeVibrationalEnergy = NULL;
    DelPointer2 (moleculePosition) ; moleculePosition = NULL;
    DelPointer2 (moleculeVelocity) ; moleculeVelocity = NULL;

    
}

//void Molecule::init()
//{
//	int nZones = 1;
//
//	
//	molecule = new Molecule * [nZones]; 
//
//	for(int i = 0; i < nZones;i++)
//	{
//		molecule[i] = new Molecule [10];
//	}
//
//	for(int i=0;i<nZones;i++)
//	{
//		Molecule * molecule = molecule[i];
//		for (int n=0;n<10;++n)
//		{
//			RDouble value[]= {n+1.0,n+1.0,n+1.0};
//			molecule->SetPP(n,value);
//		}
//	}
//
//	for(int i=0;i<nZones;i++)
//	{
//		Molecule * molecule = molecules[i];
//
//		RDouble ** PX = molecule->GetPP();
//		//for(int j = 0;j<10;j++)
//		//{
//		//	cout<<"moleculeLocation"<<":"<<i<<" "<<j<<" "<<PX[0][j]<<endl;
//		//}
//	}
//
//	for(int i=0;i<nZones;i++)
//	{
//		delete [] molecules[i];
//	}
//	delete [] molecules;
//}
//
//void Molecule::SetMoleculePosition(int n, RDouble value[])
//{
//	this->moleculePosition[0][n] = value[0];
//	this->moleculePosition[1][n] = value[1];
//	this->moleculePosition[2][n] = value[2];
//
//}
//
//void Molecule::SetMoleculeVelocity(int n, RDouble value[])
//{
//	this->moleculeVelocity[0][n] = value[0];
//	this->moleculeVelocity[1][n] = value[1];
//	this->moleculeVelocity[2][n] = value[2];
//}

RDouble Molecule::CalMoleculeRotationalEnergy(int rotationalDegreeOfFreedom, RDouble temperature)
{

    RDouble eng,a,erm,b;
    eng = 0.0;
    if(rotationalDegreeOfFreedom == 2)
    {
        eng = -log(myrand())*BOLTZ*temperature;
    
    }
    else
    {
        a = 0.5*rotationalDegreeOfFreedom - 1.0;
        while(1)
        {
            erm = 10.0*myrand();
            b = pow(erm/a,a)*exp(a-erm);
            if(b > myrand()) break;
        }
    }
    return eng;
}

RDouble Molecule::CalMoleculeVibrationalEnergy(int vibrationalDegreeOfFreedom, RDouble temperature, RDouble vibtemperature)
{

    RDouble eng;
    eng = 0.0;
    if(vibrationalDegreeOfFreedom == 2)
    {
        int ivib = static_cast <int> (-log(myrand())*temperature/vibtemperature);
        eng = ivib *BOLTZ* vibtemperature;
    }
    return eng;
}

}