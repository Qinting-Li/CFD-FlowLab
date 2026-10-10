#include "ParticlePointGroup.h"
#include "GlobalDataBase.h"
#include "Pointer.h"
#include "TK_Exit.h"
#include "TK_Parse.h"
#include "Geo_StructGrid.h"
#include "Glb_Dimension.h"

namespace PHSPACE
{

ParticlePointGroup::ParticlePointGroup()
{
    this->particlePointVariable = new map<int,OnePointVariable*>;
    this->cellParticleID = new IndexCell();
}

ParticlePointGroup::~ParticlePointGroup()
{
    //! Destructors can cause many problems of deletion by mistake.
    //! We must be careful.
    this->DeleteMap();
    delete cellParticleID;
}

void ParticlePointGroup::InitCellIndex(Grid *grid)
{
    this->cellParticleID->InitNumOfCell(grid);
    this->cellParticleID->InitCellID(grid);
    this->cellParticleID->SetBoundaryCellIndex(grid);
}

int ParticlePointGroup::AddParticle(Grid *grid, int idParticle, OnePointVariable *onePointVariable,bool &inZone)
{
    //! Judge the particle and add to IndexCell.
    inZone = false;
    int IDCellOfParticle = this->cellParticleID->InitIndexOfParticleInCell(grid, idParticle, onePointVariable, inZone);
    onePointVariable->SetIsParticleInZone(inZone);

    if (inZone)
    {
        onePointVariable->SetCellID(IDCellOfParticle);
        onePointVariable->SetIsParticleInZone(inZone);
        this->AddParticle2Map(idParticle, onePointVariable);
    }
    
    return IDCellOfParticle;
}

void ParticlePointGroup::QuickAddParticle(int idParticle,int IDCellOfParticle, OnePointVariable *onePointVariable, bool &inZone)
{
    onePointVariable->SetIsParticleInZone(inZone);

    if (inZone)
    {
        onePointVariable->SetCellID(IDCellOfParticle);
        this->AddParticle2Map(idParticle, onePointVariable);
    }
}

void ParticlePointGroup::RemoveParticle(int particleID)
{
    ParticleIDIter iter;
    iter = this->particlePointVariable->find(particleID);
    int IDCellOfParticle = iter->second->GetCellID();
    if (iter != particlePointVariable->end())
    {
        delete iter->second;
        iter->second = nullptr;
        particlePointVariable->erase(iter);
    }
    else
    {
        TK_Exit::UnexpectedVarValue("no particleID can be remove :", particleID);
    }

    this->cellParticleID->RemoveParticleID(IDCellOfParticle, particleID);
}

void ParticlePointGroup::UpdateParticle(Grid *grid)
{
    for (map<int, OnePointVariable*>::iterator iter = GetIterBegin(); iter != GetIterEnd(); ++iter)
    {
        int particleID = iter->first;
        int cellID = iter->second->GetCellID();
        SPDouble *particleCoordinate = (iter->second->GetOnePointCoordinate());
        bool inZone = false;
        int IDCellOfParticle = this->cellParticleID->SearchParticleCell(grid, particleCoordinate, inZone);
        if (IDCellOfParticle != cellID)
        {
            if (-1 == IDCellOfParticle)
            {
                //! This particle is not in current zone
                continue;
            }
            else
            {
                //! Add particle ID.
                cellParticleID->AddParticleID(IDCellOfParticle, particleID);
                cellParticleID->RemoveParticleID(cellID, particleID);
                iter->second->SetCellID(IDCellOfParticle);
            }
        }
    }
}

bool ParticlePointGroup::isParticleInZone(int particleID)
{
    bool inZone = false;

    ParticleIDIter iter;
    iter = this->particlePointVariable->find(particleID);
    if (iter != particlePointVariable->end())
    {
        inZone = iter->second->IsParticleInZone();
    }

    return inZone;
}

void ParticlePointGroup::AddParticle2Map(int particleID, OnePointVariable *onePointVariable)
{
    bool check = this->particlePointVariable->insert(make_pair(particleID, onePointVariable)).second;
    if (!check)
    {
        ostringstream oss;
        oss << "Error: insert particle: "<<"particleID = "<< particleID << endl;
        TK_Exit::ExceptionExit(oss);
    }
}

OnePointVariable *ParticlePointGroup::GetOnePointVariableByParticleID(int particleID)
{
    ParticleIDIter iter;
    iter = this->particlePointVariable->find(particleID);
    if (iter != particlePointVariable->end())
    {
        return iter->second;
    }
    else
    {
        return 0;
    }
}

void ParticlePointGroup::PrintParticleVariable2Window(int particleID)
{
    ParticleIDIter iter;
    iter = this->particlePointVariable->find(particleID);
    if (iter != particlePointVariable->end())
    {
        cout << " check particleID  = " << particleID << endl;
        iter->second->Print2Window();
    }
    else
    {
        cout << " Un exist particleID  = " << particleID << endl;
    }
}

map<int, OnePointVariable*> *ParticlePointGroup::GetMap()
{
    return this->particlePointVariable;
}

OnePointVariable *ParticlePointGroup::GetOnePoint(int particleID)
{
    OnePointVariable *onePointVariable = 0;
    ParticleIDIter iter;
    iter = this->particlePointVariable->find(particleID);
    if (iter != particlePointVariable->end())
    {
        onePointVariable = iter->second;
    }
    else
    {
        cout << " Un exist particleID  = " << particleID << endl;
    }

    return onePointVariable;
}

map<int, OnePointVariable*>::iterator ParticlePointGroup::GetIterBegin()
{
    return this->particlePointVariable->begin();
}

map<int, OnePointVariable*>::iterator ParticlePointGroup::GetIterEnd()
{
    return this->particlePointVariable->end();
}

SPDouble *ParticlePointGroup::GetParticleCoordinate(int particleID)
{
    ParticleIDIter iter;
    iter = this->particlePointVariable->find(particleID);
    if (iter != particlePointVariable->end())
    {
        return iter->second->GetOnePointCoordinate();
    }
    else
    {
        cout << "Un exist particleID  = " << particleID << endl;
        return nullptr;
    }
}

SPDouble *ParticlePointGroup::GetParticleVelocity(int particleID)
{
    ParticleIDIter iter;
    iter = this->particlePointVariable->find(particleID);
    if (iter != particlePointVariable->end())
    {
        return iter->second->GetOnePointVelocity();
    }
    else
    {
        cout << "Un exist particleID  = " << particleID << endl;
        return nullptr;
    }
}

RDouble *ParticlePointGroup::GetParticleDiameter(RDouble refLength)
{
    RDouble *particleDiameter = new RDouble[GetNumOfLocalParticle()];
    int countParticle = 0;
    ParticleIDIter iter;
    for (iter = particlePointVariable->begin(); iter != particlePointVariable->end(); ++iter)
    {
        particleDiameter[countParticle] = iter->second->GetOnePointDiameter();
        particleDiameter[countParticle] *= refLength;
        countParticle++;
    }
    return particleDiameter;
}

RDouble *ParticlePointGroup::GetParticleSpecificHeatCapacity(RDouble refVelocityDimensitional, RDouble refTemperatureDimensional)
{
    RDouble *particleSpecificHeatCapacity = new RDouble[GetNumOfLocalParticle()];
    int countParticle = 0;
    ParticleIDIter iter;
    for (iter = particlePointVariable->begin(); iter != particlePointVariable->end(); ++iter)
    {
        particleSpecificHeatCapacity[countParticle] = iter->second->GetOnePointSpecificHeatCapacity();
        particleSpecificHeatCapacity[countParticle] *= (refVelocityDimensitional* refVelocityDimensitional/ refTemperatureDimensional);
        countParticle++;
    }
    return particleSpecificHeatCapacity;
}

RDouble *ParticlePointGroup::GetParticleConvectiveHeatTransfer(RDouble refDensityDimensitional,RDouble refVelocityDimensitional)
{
    RDouble *particleConvectiveHeatTransfer = new RDouble[GetNumOfLocalParticle()];
    int countParticle = 0;
    ParticleIDIter iter;
    for (iter = particlePointVariable->begin(); iter != particlePointVariable->end(); ++iter)
    {
        particleConvectiveHeatTransfer[countParticle] = iter->second->GetOnePointParticleConvectiveHeatTransfer();
        particleConvectiveHeatTransfer[countParticle] *= (refDensityDimensitional*refVelocityDimensitional * refVelocityDimensitional * refVelocityDimensitional);
        countParticle++;
    }
    return particleConvectiveHeatTransfer;
}

RDouble *ParticlePointGroup::GetParticleRadiativeHeatTransfer(RDouble refDensityDimensitional, RDouble refVelocityDimensitional)
{
    RDouble *particleRadiativeHeatTransfer = new RDouble[GetNumOfLocalParticle()];
    int countParticle = 0;
    ParticleIDIter iter;
    for (iter = particlePointVariable->begin(); iter != particlePointVariable->end(); ++iter)
    {
        particleRadiativeHeatTransfer[countParticle] = iter->second->GetOnePointParticleRadiativeHeatTransfer();
        particleRadiativeHeatTransfer[countParticle] *= (refDensityDimensitional * refVelocityDimensitional * refVelocityDimensitional * refVelocityDimensitional);
        countParticle++;
    }
    return particleRadiativeHeatTransfer;
}

RDouble *ParticlePointGroup::GetParticleNuseltNum()
{
    RDouble *particleNuseltNum = new RDouble[GetNumOfLocalParticle()];
    int countParticle = 0;
    ParticleIDIter iter;
    for (iter = particlePointVariable->begin(); iter != particlePointVariable->end(); ++iter)
    {
        particleNuseltNum[countParticle] = iter->second->GetOnePointParticleNuseltNum();
        countParticle++;
    }
    return particleNuseltNum;
}

RDouble *ParticlePointGroup::GetParticleWork(RDouble refDensityDimensitional, RDouble refVelocityDimensitional)
{
    RDouble *particleWork = new RDouble[GetNumOfLocalParticle()];
    int countParticle = 0;
    ParticleIDIter iter;
    for (iter = particlePointVariable->begin(); iter != particlePointVariable->end(); ++iter)
    {
        particleWork[countParticle] = iter->second->GetOnePointParticleWork();
        particleWork[countParticle] *= (refDensityDimensitional * refVelocityDimensitional * refVelocityDimensitional * refVelocityDimensitional);
        countParticle++;
    }
    return particleWork;
}

RDouble *ParticlePointGroup::GetParticleDissipation(RDouble refDensityDimensitional, RDouble refVelocityDimensitional)
{
    RDouble *particleDissipation = new RDouble[GetNumOfLocalParticle()];
    int countParticle = 0;
    ParticleIDIter iter;
    for (iter = particlePointVariable->begin(); iter != particlePointVariable->end(); ++iter)
    {
        particleDissipation[countParticle] = iter->second->GetOnePointParticleDissipation();
        particleDissipation[countParticle] *= (refDensityDimensitional * refVelocityDimensitional * refVelocityDimensitional * refVelocityDimensitional);
        countParticle++;
    }
    return particleDissipation;
}

RDouble *ParticlePointGroup::GetParticleDensity(RDouble refDensity)
{
    RDouble *particleDensity = new RDouble[GetNumOfLocalParticle()];
    int countParticle = 0;
    ParticleIDIter iter;
    for (iter = particlePointVariable->begin(); iter != particlePointVariable->end(); ++iter)
    {
        particleDensity[countParticle] = iter->second->GetOnePointDensity();
        particleDensity[countParticle] *= refDensity;
        countParticle++;
    }
    return particleDensity;
}

RDouble *ParticlePointGroup::GetParticleReNum()
{
    RDouble *particleReNum = new RDouble[GetNumOfLocalParticle()];
    int countParticle = 0;
    ParticleIDIter iter;
    for (iter = particlePointVariable->begin(); iter != particlePointVariable->end(); ++iter)
    {
        particleReNum[countParticle] = iter->second->GetOnePointParticleReNum();
        countParticle++;
    }
    return particleReNum;
}

RDouble *ParticlePointGroup::GetParticleMaNum()
{
    RDouble *particleMaNum = new RDouble[GetNumOfLocalParticle()];
    int countParticle = 0;
    ParticleIDIter iter;
    for (iter = particlePointVariable->begin(); iter != particlePointVariable->end(); ++iter)
    {
        particleMaNum[countParticle] = iter->second->GetOnePointParticleMaNum();
        countParticle++;
    }
    return particleMaNum;
}


RDouble *ParticlePointGroup::GetParticleCoordinateDim(RDouble refLength,int iDim)
{
    int nDim = 3;
    RDouble *particleCoordinate = new RDouble[GetNumOfLocalParticle()];
    //! Loop for local particle.
    int countParticle = 0;
    ParticleIDIter iter;
    for (iter = particlePointVariable->begin(); iter != particlePointVariable->end(); ++iter)
    {
        particleCoordinate[countParticle] = iter->second->GetOnePointCoordinateOnDirction(iDim);
        particleCoordinate[countParticle] *= refLength;
        ++countParticle;
    }
    return particleCoordinate;
}

RDouble *ParticlePointGroup::GetParticleVelocityDim(RDouble refVelocity,int iDim)
{
    int nDim = 3;
    RDouble *particleVelocity = new RDouble[GetNumOfLocalParticle()];
    //! Loop for local particle.
    int countParticle = 0;
    ParticleIDIter iter;
    for (iter = particlePointVariable->begin(); iter != particlePointVariable->end(); ++iter)
    {
        particleVelocity[countParticle] = iter->second->GetOnePointVelocityOnDirction(iDim);
        particleVelocity[countParticle] *= refVelocity;
        ++countParticle;
    }
    return particleVelocity;
}


RDouble *ParticlePointGroup::GetParticleRelativeVelocityDim(RDouble refVelocity, int iDim)
{
    int nDim = 3;
    RDouble *particleRelativeVelocity = new RDouble[GetNumOfLocalParticle()];
    //! Loop for local particle.
    int countParticle = 0;
    ParticleIDIter iter;
    for (iter = particlePointVariable->begin(); iter != particlePointVariable->end(); ++iter)
    {
        particleRelativeVelocity[countParticle] = iter->second->GetOnePointRelativeVelocityOnDirction(iDim);
        particleRelativeVelocity[countParticle] *= refVelocity;
        ++countParticle;
    }
    return particleRelativeVelocity;
}


RDouble *ParticlePointGroup::GetParticleTotalForceSt()
{
    int nDim = 3;
    RDouble *particleTotalForceSt = new RDouble[GetNumOfLocalParticle()];
    //! Loop for local particle.
    int countParticle = 0;
    ParticleIDIter iter;
    for (iter = particlePointVariable->begin(); iter != particlePointVariable->end(); ++iter)
    {
        particleTotalForceSt[countParticle] = iter->second->GetOnePointTotalForceSt();
        ++countParticle;
    }
    return particleTotalForceSt;
}

RDouble *ParticlePointGroup::GetParticleStokesSt()
{
    int nDim = 3;
    RDouble *particleStokesSt = new RDouble[GetNumOfLocalParticle()];
    //! Loop for local particle.
    int countParticle = 0;
    ParticleIDIter iter;
    for (iter = particlePointVariable->begin(); iter != particlePointVariable->end(); ++iter)
    {
        particleStokesSt[countParticle] = iter->second->GetOnePointStokesSt();
        ++countParticle;
    }
    return particleStokesSt;
}


RDouble *ParticlePointGroup::GetParticleModifiedSt()
{
    int nDim = 3;
    RDouble *particleModifiedSt = new RDouble[GetNumOfLocalParticle()];
    //! Loop for local particle.
    int countParticle = 0;
    ParticleIDIter iter;
    for (iter = particlePointVariable->begin(); iter != particlePointVariable->end(); ++iter)
    {
        particleModifiedSt[countParticle] = iter->second->GetOnePointModifiedSt();
        ++countParticle;
    }
    return particleModifiedSt;
}

RDouble *ParticlePointGroup::GetParticleAngularVelocityDim(RDouble refAngularVel,int iDim)
{
    int nDim = 3;
    RDouble *particleAngularVelocity = new RDouble[GetNumOfLocalParticle()];
    //! Loop for local particle.
    int countParticle = 0;
    ParticleIDIter iter;
    for (iter = particlePointVariable->begin(); iter != particlePointVariable->end(); ++iter)
    {
        particleAngularVelocity[countParticle] = iter->second->GetOnePointAngularVelocityOnDirction(iDim);
        particleAngularVelocity[countParticle] *= refAngularVel;
        ++countParticle;
    }
    return particleAngularVelocity;
}

RDouble *ParticlePointGroup::GetParticleAccelerationDim(RDouble refAcceleration,int iDim)
{
    int nDim = 3;
    RDouble *particleAcceleration = new RDouble[GetNumOfLocalParticle()];
    //! Loop for local particle.
    int countParticle = 0;
    ParticleIDIter iter;
    for (iter = particlePointVariable->begin(); iter != particlePointVariable->end(); ++iter)
    {

        particleAcceleration[countParticle] = iter->second->GetOnePointAccelerationOnDirction(iDim);
        particleAcceleration[countParticle] *= refAcceleration;
        ++countParticle;
    }

    return particleAcceleration;
}

RDouble *ParticlePointGroup::GetParticleTemperatureDim(RDouble refTemperatuer,int iDim)
{
    int nDim = 1;
    RDouble *particleTemperature = new RDouble[GetNumOfLocalParticle()];
    //! Loop for local particle.
    int countParticle = 0;
    ParticleIDIter iter;
    for ( iter = particlePointVariable->begin(); iter != particlePointVariable->end(); ++iter)
    {
        particleTemperature[countParticle] = iter->second->GetOnePointTemperatureOnDirction(iDim);
        particleTemperature[countParticle] *= refTemperatuer;
        ++countParticle;
    }
    return particleTemperature;
}

RDouble *ParticlePointGroup::GetParticleForceDim(RDouble refForce,int forceType, int iDim)
{
    int nDim = 3;
    RDouble *particleForce = new RDouble[GetNumOfLocalParticle()];
    //! Loop for local particle.
    int countParticle = 0;
    ParticleIDIter iter;
    for (iter = particlePointVariable->begin(); iter != particlePointVariable->end(); ++iter)
    {
        particleForce[countParticle] = iter->second->GetOnePointForceOnDirction(iDim, forceType);
        particleForce[countParticle] *= refForce;
        ++countParticle;
    }
    return particleForce;
}

IndexCell *ParticlePointGroup::GetIndexCell()
{
    return this->cellParticleID;
}

int ParticlePointGroup::GetCellIDOfParticle(int particleID)
{
    ParticleIDIter iter;
    iter = this->particlePointVariable->find(particleID);
    if (iter != particlePointVariable->end())
    {
        return iter->second->GetCellID();
    }
    else
    {
        cout << "Un exist particleID  = " << particleID <<"###Can not find particle ID" << endl;
        return NULL;
    }
}

int ParticlePointGroup::GetGlobalParticleIDByLocalIParticle(int iParticle)
{
    int particleID = -1;
    ParticleIDIter iter;
    int count = 0;
    for (iter = particlePointVariable->begin(); iter != particlePointVariable->end(); ++iter)
    {
        if (iParticle == count)
        {
            particleID = count;
            break;
        }
        else
        {
            count++;
        }
    }
    return particleID;
}

int ParticlePointGroup::GetNumOfLocalParticle()
{
    return static_cast<int>(this->particlePointVariable->size());
}

int ParticlePointGroup::GetNumOfLocalCell()
{
    return static_cast<int>(this->cellParticleID->GetNumOfCell());
}

int ParticlePointGroup::GetNumOfParticleOfCell(int cellID)
{
    return this->cellParticleID->GetNumOfParticleOnCell(cellID);
}

int *ParticlePointGroup::GetIDOfParticleOfCell(int cellID)
{
    return this->cellParticleID->GetIDOfParticleForCell(cellID);
}

int *ParticlePointGroup::GetGlobalParticleIDOnLocalZone()
{
    int nLocalParticle = particlePointVariable->size();
    int *particleIDList = new int[nLocalParticle];
    int count = 0;
    for (map<int, OnePointVariable*>::iterator iter = particlePointVariable->begin(); iter != particlePointVariable->end(); ++iter)
    {
        particleIDList[count] = iter->first;
        count++;
    }
    return particleIDList;
}

void ParticlePointGroup::DeleteMap()
{
    if (0 == this)
    {
        return;
    }

    ParticleIDIter iter = particlePointVariable->begin();
    while (iter != particlePointVariable->end())
    {
        delete iter->second;
        iter->second = NULL;
        particlePointVariable->erase(iter++);
    }

    delete particlePointVariable;
    particlePointVariable = NULL;
}

void ParticlePointGroup::CompressData(StructGrid *grid,DataContainer *&cdata,int &particleID,int &zoneIDGlobalRecv, int &size)
{
    OnePointVariable *onePointVariable = this->GetOnePointVariableByParticleID(particleID);

    size = 0;

    //! int data
    cdata->Write(&particleID, sizeof(int));

    size += sizeof(int);

    //! int data
    cdata->Write(&zoneIDGlobalRecv, sizeof(int));

    size += sizeof(int);

    //! Get index IJK of particle position in current zone.
    int IDCellOfParticle = onePointVariable->GetCellID();
    int nI, nJ, nK;
    int iP, jP, kP;
    grid->GetND(nI, nJ, nK);
    nI -= 1;
    nJ -= 1;
    nK -= 1;
    if (TWO_D == GetDim())
    {
        nK = 1;
    }
    GetIndexCellIJKStruct(nI, nJ, nK, iP, jP, kP, GetNumberOfGhostCellLayers(), IDCellOfParticle);

    cout << "send " << " ip " << iP << " " << jP << endl;
    if (iP > nI)
    {
        iP = 1;
    }
    else if(iP < 1)
    {
        iP = -1;
    }

    if (jP > nJ)
    {
        jP = 1;
    }
    else if(jP < 1)
    {
        jP = -1;
    }

    if (kP > nK)
    {
        kP = 1;
    }
    else if(kP < 1)
    {
        kP = -1;
    }

    cdata->Write(&iP, sizeof(int));
    cdata->Write(&jP, sizeof(int));
    cdata->Write(&kP, sizeof(int));

    size += sizeof(int);
    size += sizeof(int);
    size += sizeof(int);

    //! double data
    SPDouble &onePointCoordinate = *(onePointVariable->GetOnePointCoordinate());
    SPDouble &onePointVelocity = *(onePointVariable->GetOnePointVelocity());
    SPDouble &onePointAcceleration = *(onePointVariable->GetOnePointAcceleration());
    SPDouble &onePointAngularVelocity = *(onePointVariable->GetOnePointAngularVelocity());

    int nDimVar = 12;
    RDouble *data = new RDouble[nDimVar];

    int index = -1;
    int indexDim = -1;
    data[++index] = onePointCoordinate[++indexDim];
    data[++index] = onePointCoordinate[++indexDim];
    data[++index] = onePointCoordinate[++indexDim];

    indexDim = -1;
    data[++index] = onePointVelocity[++indexDim];
    data[++index] = onePointVelocity[++indexDim];
    data[++index] = onePointVelocity[++indexDim];

    indexDim = -1;
    data[++index] = onePointAcceleration[++indexDim];
    data[++index] = onePointAcceleration[++indexDim];
    data[++index] = onePointAcceleration[++indexDim];

    indexDim = -1;
    data[++index] = onePointAngularVelocity[++indexDim];
    data[++index] = onePointAngularVelocity[++indexDim];
    data[++index] = onePointAngularVelocity[++indexDim];

    for (int i = 0; i < nDimVar; ++i)
    {
        cdata->Write(&data[i], sizeof(RDouble));
        size += sizeof(RDouble);
    }

    delete [] data;    data = nullptr;
}

};