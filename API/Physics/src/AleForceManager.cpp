#include "AleForceManager.h"
#include "Math_BasisFunction.h"
#include "Constants.h"
#include "PHIO.h"
#include "IO_FileName.h"
#include "Geo_SimpleBC.h"
#include "TK_Time.h"
#include "TK_Parse.h"
#include "TK_Exit.h"
#include "Glb_Dimension.h"
#include "GlobalDataBase.h"
#include "PHHeader.h"
#include "PHMpi.h"
#include "AleManager.h"

using namespace std;

namespace PHSPACE
{
bool AerodynamicForceBody::initialize = false;
vector< std::set< string > >  AerodynamicForceBody::bodys;

AerodynamicForceBody::AerodynamicForceBody()
{
}

AerodynamicForceBody::~AerodynamicForceBody()
{
}

void AerodynamicForceBody::Initialize()
{
    if (!AerodynamicForceBody::initialize)
    {
        int numberOfMovingBodies = GlobalDataBase::GetIntParaFromDB("numberOfMovingBodies");
        AerodynamicForceBody::bodys.resize(numberOfMovingBodies);
        SimpleBC *boundaryCondition = 0;
        int iMovingBody = 0;
        int nTBC = GlobalBoundaryCondition::GetNumberOfBC();
        for (int iBC = 0; iBC < nTBC; ++ iBC)
        {
            if (GlobalBoundaryCondition::IsSolidWallBC(iBC))
            {
                boundaryCondition = GlobalBoundaryCondition::GetBC(iBC);
                string bodyName = boundaryCondition->GetBodyName();
                AerodynamicForceBody::bodys[iMovingBody++].insert(bodyName);
            }
        }

        AerodynamicForceBody::initialize = true;
    }
}

bool AerodynamicForceBody::FaceInAerodynamicForceBodys(string bodyName, int aerodynamicForceBodyIndex)
{
    AerodynamicForceBody::Initialize();
    set< string >::iterator iter = AerodynamicForceBody::bodys[aerodynamicForceBodyIndex].find(bodyName);
    return iter != AerodynamicForceBody::bodys[aerodynamicForceBodyIndex].end();
}

bool FaceInAerodynamicForceBodys(string bodyName, int aerodynamicForceBodyIndex)
{
    return AerodynamicForceBody::FaceInAerodynamicForceBodys(bodyName, aerodynamicForceBodyIndex);
}

AerodynamicForceManager aerodynamicForceManager;

void InitializeAerodynamicForce()
{
    AerodynamicForceManager *aerodynamicForceManager = PHSPACE::GetAerodynamicForceManager();
    aerodynamicForceManager->InitializeAerodynamicForceManager();
}

void AddAerodynamicForce(ActionKey *actionKey)
{
    AerodynamicForceManager *aerodynamicForceManager = PHSPACE::GetAerodynamicForceManager();
    aerodynamicForceManager->AddAerodynamicForce(actionKey);
}

void DumpAerodynamicForceCoefficient(ActionKey *actionKey)
{
    AerodynamicForceManager *aerodynamicForceManager = PHSPACE::GetAerodynamicForceManager();
    aerodynamicForceManager->DumpAerodynamicForceCoefficient(actionKey);
}

void CollectAerodynamicForceOfAllProcessors()
{
    AerodynamicForceManager *aerodynamicForceManager = PHSPACE::GetAerodynamicForceManager();
    aerodynamicForceManager->CollectAerodynamicForceOfAllProcessors();
}

BasicAerodynamicForce **CreateBasicAerodynamicForceContainer()
{
    int numberOfAerodynamicForceBodys = GlobalDataBase::GetIntParaFromDB("numberOfMovingBodies");

    BasicAerodynamicForce **basicAerodynamicForceContainer = new BasicAerodynamicForce * [numberOfAerodynamicForceBodys + 1];

    for (int iMovingBody = 0; iMovingBody <= numberOfAerodynamicForceBodys; ++ iMovingBody)
    {
        basicAerodynamicForceContainer[iMovingBody] = new BasicAerodynamicForce();
    }

    return basicAerodynamicForceContainer;
}

vector< BasicAerodynamicForce * > *CreateBasicAerodynamicForceVector()
{
    int numberOfAerodynamicForceBodys = GlobalDataBase::GetIntParaFromDB("numberOfMovingBodies");

    vector< BasicAerodynamicForce * > *basicAerodynamicForceVector = new vector< BasicAerodynamicForce * >(numberOfAerodynamicForceBodys + 1);

    for (int iMovingBody = 0; iMovingBody <= numberOfAerodynamicForceBodys; ++ iMovingBody)
    {
        (*basicAerodynamicForceVector)[iMovingBody] = new BasicAerodynamicForce();
    }

    return basicAerodynamicForceVector;
}

void FreeBasicAerodynamicForceVector(vector<BasicAerodynamicForce *> *basicAerodynamicForceVector)
{
    int numberOfAerodynamicForceBodys = basicAerodynamicForceVector->size();

    for (int iMovingBody = 0; iMovingBody < numberOfAerodynamicForceBodys; ++ iMovingBody)
    {
        delete (*basicAerodynamicForceVector)[iMovingBody];
    }

    delete  basicAerodynamicForceVector;
}

void FreeBasicAerodynamicForceContainer(BasicAerodynamicForce **basicAerodynamicForceContainer)
{
    int numberOfAerodynamicForceBodys = GlobalDataBase::GetIntParaFromDB("numberOfMovingBodies");

    for (int iMovingBody = 0; iMovingBody <= numberOfAerodynamicForceBodys; ++ iMovingBody)
    {
        delete basicAerodynamicForceContainer[iMovingBody];
    }

    delete[] basicAerodynamicForceContainer;
}

void WriteAerodynamicForceToActionKey(ActionKey *actionKey, BasicAerodynamicForce **basicAerodynamicForceContainer)
{
    int numberOfAerodynamicForceBodys = GlobalDataBase::GetIntParaFromDB("numberOfMovingBodies");

    DataContainer *dataContainer = actionKey->GetData();
    dataContainer->moveto_begin();

    for (int iMovingBody = 0; iMovingBody <= numberOfAerodynamicForceBodys; ++ iMovingBody)
    {
        basicAerodynamicForceContainer[iMovingBody]->WriteAerodynamicForceToDataContainer(dataContainer);
    }
}

void WriteAerodynamicForceToActionKey(ActionKey *actionKey, vector< BasicAerodynamicForce * > *basicAerodynamicForceVector)
{
    int numberOfAerodynamicForceBodys = GlobalDataBase::GetIntParaFromDB("numberOfMovingBodies");

    DataContainer *dataContainer = actionKey->GetData();
    dataContainer->moveto_begin();

    for (int iMovingBody = 0; iMovingBody <= numberOfAerodynamicForceBodys; ++ iMovingBody)
    {
        (*basicAerodynamicForceVector)[iMovingBody]->WriteAerodynamicForceToDataContainer(dataContainer);
    }
}

BasicAerodynamicForce::BasicAerodynamicForce()
{
    this->ZeroBasicAerodynamicForce();
}

BasicAerodynamicForce::~BasicAerodynamicForce()
{
}

void BasicAerodynamicForce::ZeroBasicAerodynamicForce()
{
    forceX = 0.0;
    forceY = 0.0;
    forceZ = 0.0;

    forceXByPressure = 0.0;
    forceYByPressure = 0.0;
    forceZByPressure = 0.0;

    forceXByViscosity = 0.0;
    forceYByViscosity = 0.0;
    forceZByViscosity = 0.0;

    momentX = 0.0;
    momentY = 0.0;
    momentZ = 0.0;

    momentXByPressure = 0.0;
    momentYByPressure = 0.0;
    momentZByPressure = 0.0;

    momentXByViscosity = 0.0;
    momentYByViscosity = 0.0;
    momentZByViscosity = 0.0;
}

void BasicAerodynamicForce::SetMomentReferencePoint(RDouble momentReferenceX, RDouble momentReferenceY, RDouble momentReferenceZ)
{
    this->momentReferenceX = momentReferenceX;
    this->momentReferenceY = momentReferenceY;
    this->momentReferenceZ = momentReferenceZ;
}

void BasicAerodynamicForce::SetPressureAerodynamicForce(RDouble forceXByPressure, RDouble forceYByPressure, RDouble forceZByPressure)
{
    this->forceXByPressure = forceXByPressure;
    this->forceYByPressure = forceYByPressure;
    this->forceZByPressure = forceZByPressure;
}

void BasicAerodynamicForce::SetViscousAerodynamicForce(RDouble forceXByViscosity, RDouble forceYByViscosity, RDouble forceZByViscosity)
{
    this->forceXByViscosity = forceXByViscosity;
    this->forceYByViscosity = forceYByViscosity;
    this->forceZByViscosity = forceZByViscosity;
}

void BasicAerodynamicForce::ComputeResultantAerodynamicForce()
{
    forceX = forceXByPressure + forceXByViscosity;
    forceY = forceYByPressure + forceYByViscosity;
    forceZ = forceZByPressure + forceZByViscosity;
}

void BasicAerodynamicForce::ComputeResultantAerodynamicMoment()
{
    momentX = momentXByPressure + momentXByViscosity;
    momentY = momentYByPressure + momentYByViscosity;
    momentZ = momentZByPressure + momentZByViscosity;
}

void BasicAerodynamicForce::SetPressureAerodynamicMoment(RDouble faceCenterX, RDouble faceCenterY, RDouble faceCenterZ)
{
    momentXByPressure = (faceCenterY - momentReferenceY) * forceZByPressure - (faceCenterZ - momentReferenceZ) * forceYByPressure;
    momentYByPressure = (faceCenterZ - momentReferenceZ) * forceXByPressure - (faceCenterX - momentReferenceX) * forceZByPressure;
    momentZByPressure = (faceCenterX - momentReferenceX) * forceYByPressure - (faceCenterY - momentReferenceY) * forceXByPressure;
}

void BasicAerodynamicForce::SetViscousAerodynamicMoment(RDouble faceCenterX, RDouble faceCenterY, RDouble faceCenterZ)
{
    momentXByViscosity = (faceCenterY - momentReferenceY) * forceZByViscosity - (faceCenterZ - momentReferenceZ) * forceYByViscosity;
    momentYByViscosity = (faceCenterZ - momentReferenceZ) * forceXByViscosity - (faceCenterX - momentReferenceX) * forceZByViscosity;
    momentZByViscosity = (faceCenterX - momentReferenceX) * forceYByViscosity - (faceCenterY - momentReferenceY) * forceXByViscosity;
}

void BasicAerodynamicForce::ComputeAerodynamicMoment(RDouble faceCenterX, RDouble faceCenterY, RDouble faceCenterZ)
{
    momentX = (faceCenterY - momentReferenceY) * forceZ - (faceCenterZ - momentReferenceZ) * forceY;
    momentY = (faceCenterZ - momentReferenceZ) * forceX - (faceCenterX - momentReferenceX) * forceZ;
    momentZ = (faceCenterX - momentReferenceX) * forceY - (faceCenterY - momentReferenceY) * forceX;
}

void BasicAerodynamicForce::AddAerodynamicForce(BasicAerodynamicForce *basicAerodynamicForce)
{
    this->forceX += basicAerodynamicForce->forceX;
    this->forceY += basicAerodynamicForce->forceY;
    this->forceZ += basicAerodynamicForce->forceZ;

    this->forceXByPressure += basicAerodynamicForce->forceXByPressure;
    this->forceYByPressure += basicAerodynamicForce->forceYByPressure;
    this->forceZByPressure += basicAerodynamicForce->forceZByPressure;

    this->forceXByViscosity += basicAerodynamicForce->forceXByViscosity;
    this->forceYByViscosity += basicAerodynamicForce->forceYByViscosity;
    this->forceZByViscosity += basicAerodynamicForce->forceZByViscosity;

    this->momentX += basicAerodynamicForce->momentX;
    this->momentY += basicAerodynamicForce->momentY;
    this->momentZ += basicAerodynamicForce->momentZ;

    this->momentXByPressure += basicAerodynamicForce->momentXByPressure;
    this->momentYByPressure += basicAerodynamicForce->momentYByPressure;
    this->momentZByPressure += basicAerodynamicForce->momentZByPressure;

    this->momentXByViscosity += basicAerodynamicForce->momentXByViscosity;
    this->momentYByViscosity += basicAerodynamicForce->momentYByViscosity;
    this->momentZByViscosity += basicAerodynamicForce->momentZByViscosity;
}

void BasicAerodynamicForce::WriteAerodynamicForceToDataContainer(DataContainer *dataContainer)
{
    PHWrite(dataContainer, forceX);
    PHWrite(dataContainer, forceY);
    PHWrite(dataContainer, forceZ);

    PHWrite(dataContainer, forceXByPressure);
    PHWrite(dataContainer, forceYByPressure);
    PHWrite(dataContainer, forceZByPressure);

    PHWrite(dataContainer, forceXByViscosity);
    PHWrite(dataContainer, forceYByViscosity);
    PHWrite(dataContainer, forceZByViscosity);

    PHWrite(dataContainer, momentX);
    PHWrite(dataContainer, momentY);
    PHWrite(dataContainer, momentZ);

    PHWrite(dataContainer, momentXByPressure);
    PHWrite(dataContainer, momentYByPressure);
    PHWrite(dataContainer, momentZByPressure);

    PHWrite(dataContainer, momentXByViscosity);
    PHWrite(dataContainer, momentYByViscosity);
    PHWrite(dataContainer, momentZByViscosity);
}

void BasicAerodynamicForce::ReadAerodynamicForceFromDataContainer(DataContainer *dataContainer)
{
    PHRead(dataContainer, forceX);
    PHRead(dataContainer, forceY);
    PHRead(dataContainer, forceZ);

    PHRead(dataContainer, forceXByPressure);
    PHRead(dataContainer, forceYByPressure);
    PHRead(dataContainer, forceZByPressure);

    PHRead(dataContainer, forceXByViscosity);
    PHRead(dataContainer, forceYByViscosity);
    PHRead(dataContainer, forceZByViscosity);

    PHRead(dataContainer, momentX);
    PHRead(dataContainer, momentY);
    PHRead(dataContainer, momentZ);

    PHRead(dataContainer, momentXByPressure);
    PHRead(dataContainer, momentYByPressure);
    PHRead(dataContainer, momentZByPressure);

    PHRead(dataContainer, momentXByViscosity);
    PHRead(dataContainer, momentYByViscosity);
    PHRead(dataContainer, momentZByViscosity);
}

BasicAleForce::BasicAleForce()
{
    this->ZeroBasicAleForce();
}

BasicAleForce::~BasicAleForce()
{
}

void BasicAleForce::ZeroBasicAleForce()
{
    forceX = 0.0;
    forceY = 0.0;
    forceZ = 0.0;

    momentX = 0.0;
    momentY = 0.0;
    momentZ = 0.0;
}

void BasicAleForce::SetBasicAleForceByAerodynamicForce(AerodynamicForce *aerodynamicForce)
{
    forceX = aerodynamicForce->GetAerodynamicForceX();
    forceY = aerodynamicForce->GetAerodynamicForceY();
    forceZ = aerodynamicForce->GetAerodynamicForceZ();

    aerodynamicForce->GetMomentByOriginOfCoordinate(momentX, momentY, momentZ);
}

void BasicAleForce::ReadBasicAleForceFromBasicAleForce(BasicAleForce *basicAleForce)
{
    (*this) = (*basicAleForce);
}

void BasicAleForce::ReadBasicAleForceFromFile(fstream &file)
{
    PHSPACE::PHRead(file, forceX);
    PHSPACE::PHRead(file, forceY);
    PHSPACE::PHRead(file, forceZ);

    PHSPACE::PHRead(file, momentX);
    PHSPACE::PHRead(file, momentY);
    PHSPACE::PHRead(file, momentZ);

    PHSPACE::PHRead(file, powerConsumption);
}

void BasicAleForce::WriteBasicAleForceToFile(fstream &file)
{
    PHSPACE::PHWrite(file, forceX);
    PHSPACE::PHWrite(file, forceY);
    PHSPACE::PHWrite(file, forceZ);

    PHSPACE::PHWrite(file, momentX);
    PHSPACE::PHWrite(file, momentY);
    PHSPACE::PHWrite(file, momentZ);

    PHSPACE::PHWrite(file, powerConsumption);
}

AleForceOfBody::AleForceOfBody()
{
    bodyIndex = 0;
    numberOfTimeLayers = 0;
}

AleForceOfBody::~AleForceOfBody()
{
    this->DeallocateMemory();
}

void AleForceOfBody::AllocateMemory()
{
    basicAleForceContainer.resize(numberOfTimeLayers);
    for (int iTimeLayer = 0; iTimeLayer < numberOfTimeLayers; ++iTimeLayer)
    {
        basicAleForceContainer[iTimeLayer] = new BasicAleForce();
    }
}

void AleForceOfBody::DeallocateMemory()
{
    for (int iTimeLayer = 0; iTimeLayer < numberOfTimeLayers; ++iTimeLayer)
    {
        delete basicAleForceContainer[iTimeLayer];
    }
}

BasicAleForce *AleForceOfBody::GetBasicAleForce(int timeLayerIndex)
{
    return basicAleForceContainer[timeLayerIndex];
}

void AleForceOfBody::DumpRestartAleForceOfBodyContainer(fstream &file)
{
    int numberOfTimeLayers = this->GetNumberOfTimeLayers();
    for (int iTimeLayer = 0; iTimeLayer < numberOfTimeLayers; ++iTimeLayer)
    {
        BasicAleForce *basicAleForce = this->GetBasicAleForce(iTimeLayer);
        basicAleForce->WriteBasicAleForceToFile(file);
    }
}

void AleForceOfBody::InitializeAleForceOfBodyContainer()
{
    int numberOfTimeLayers = this->GetNumberOfTimeLayers();
    int bodyIndex = this->GetBodyIndex();
    for (int iTimeLayer = 0; iTimeLayer < numberOfTimeLayers; ++iTimeLayer)
    {
        BasicAleForce *basicAleForce = this->GetBasicAleForce(iTimeLayer);
        AerodynamicForce *aerodynamicForce = PHSPACE::GetCurrentAerodynamicForce(bodyIndex);
        basicAleForce->SetBasicAleForceByAerodynamicForce(aerodynamicForce);
    }
}

void AleForceOfBody::ReadAleForceOfBodyContainer()
{
    fstream &aleFile = PHSPACE::GetCurrentAleFile();

    int numberOfTimeLayers = this->GetNumberOfTimeLayers();
    for (int iTimeLayer = 0; iTimeLayer < numberOfTimeLayers; ++iTimeLayer)
    {
        BasicAleForce *basicAleForce = this->GetBasicAleForce(iTimeLayer);
        fstream &file = aleFile;
        basicAleForce->ReadBasicAleForceFromFile(file);
    }
}

void AleForceOfBody::UpdateAleForceOfBodyOfNewTimeLayer()
{
    int iTimeLayer = 0;
    int bodyIndex = this->GetBodyIndex();

    AerodynamicForce *aerodynamicForce = PHSPACE::GetCurrentAerodynamicForce(bodyIndex);

    BasicAleForce *basicAleForce = this->GetBasicAleForce(iTimeLayer);
    basicAleForce->SetBasicAleForceByAerodynamicForce(aerodynamicForce);
}

void AleForceOfBody::UpdateUnsteadyAleForceTerm()
{
    int iTimeLevel0 = 0;
    int iTimeLevel1 = 1;
    int iTimeLevel2 = 2;

    BasicAleForce *basicAleForce = this->GetBasicAleForce(iTimeLevel0);
    BasicAleForce *basicAleForceNTimeLevel = this->GetBasicAleForce(iTimeLevel1);
    BasicAleForce *basicAleForceN1TimeLevel = this->GetBasicAleForce(iTimeLevel2);

    PHSPACE::CopyBasicAleForce(basicAleForceN1TimeLevel, basicAleForceNTimeLevel);
    PHSPACE::CopyBasicAleForce(basicAleForceNTimeLevel, basicAleForce);
}

AleForceManager::AleForceManager()
{
    numberOfMovingBodies = 0;

    dimensionalForceCoefficient = 1.0;
    dimensionalMomentCoefficient = 1.0;
    dimensionalPowerCoefficient = 1.0;
    numberOfTimeLayers = 0;
}

AleForceManager::~AleForceManager()
{
    this->DeallocateAleForceOfBodyContainer();
}

void AleForceManager::Initialize()
{
    this->AllocateAleForceOfBodyContainer();
    this->InitializeAleForceOfBodyContainer();
}

void AleForceManager::Restart()
{
    this->ReadAleForceOfBodyContainer();
}

void AleForceManager::InitializeDimensionalReferenceVariables()
{
    RDouble referenceDensity = GlobalDataBase::GetDoubleParaFromDB("referenceDensity");
    RDouble referenceVelocity = GlobalDataBase::GetDoubleParaFromDB("referenceVelocity");
    RDouble referenceLength = GlobalDataBase::GetDoubleParaFromDB("referenceLength");

    int geometricDimension = GetDim();

    RDouble dynamicPressure = referenceDensity * PHSPACE::SQR(referenceVelocity);

    dimensionalForceCoefficient = half * dynamicPressure * PHSPACE::SQR(referenceLength);
    dimensionalMomentCoefficient = half * dynamicPressure * PHSPACE::POWER3(referenceLength);
    dimensionalPowerCoefficient = dimensionalForceCoefficient * referenceVelocity;

    if (geometricDimension == TWO_D)
    {
        dimensionalForceCoefficient = half * dynamicPressure * referenceLength;
        dimensionalMomentCoefficient = half * dynamicPressure * PHSPACE::SQR(referenceLength);
    }
}

void AleForceManager::AllocateAleForceOfBodyContainer()
{
    numberOfTimeLayers = 3;
    numberOfMovingBodies = GlobalDataBase::GetIntParaFromDB("numberOfMovingBodies");
    this->InitializeDimensionalReferenceVariables();

    aleForceOfBodyContainer.resize(numberOfMovingBodies);

    for (int iMovingBody = 0; iMovingBody < numberOfMovingBodies; ++ iMovingBody)
    {
        AleForceOfBody *aleForceOfBody = new AleForceOfBody();
        aleForceOfBodyContainer[iMovingBody] = aleForceOfBody;

        aleForceOfBody->SetBodyIndex(iMovingBody);
        aleForceOfBody->SetNumberOfTimeLayers(numberOfTimeLayers);
        aleForceOfBody->AllocateMemory();
    }
}

void AleForceManager::DeallocateAleForceOfBodyContainer()
{
    for (int iMovingBody = 0; iMovingBody < numberOfMovingBodies; ++ iMovingBody)
    {
        delete aleForceOfBodyContainer[iMovingBody];
    }
}

void AleForceManager::InitializeAleForceOfBodyContainer()
{
    for (int iMovingBody = 0; iMovingBody < numberOfMovingBodies; ++ iMovingBody)
    {
        AleForceOfBody *aleForceOfBody = this->GetAleForceOfBody(iMovingBody);

        aleForceOfBody->InitializeAleForceOfBodyContainer();
    }
}

void AleForceManager::ReadAleForceOfBodyContainer()
{
    for (int iMovingBody = 0; iMovingBody < numberOfMovingBodies; ++ iMovingBody)
    {
        AleForceOfBody *aleForceOfBody = this->GetAleForceOfBody(iMovingBody);

        aleForceOfBody->ReadAleForceOfBodyContainer();
    }
}

void AleForceManager::Run()
{
    for (int iMovingBody = 0; iMovingBody < numberOfMovingBodies; ++ iMovingBody)
    {
        aleForceOfBodyContainer[iMovingBody]->UpdateAleForceOfBodyOfNewTimeLayer();
    }
}

void AleForceManager::Post()
{
    int numberOfMovingBodies = GlobalDataBase::GetIntParaFromDB("numberOfMovingBodies");

    for (int iMovingBody = 0; iMovingBody < numberOfMovingBodies; ++ iMovingBody)
    {
        AleForceOfBody *aleForceOfBody = this->GetAleForceOfBody(iMovingBody);
        aleForceOfBody->UpdateUnsteadyAleForceTerm();
    }
}

void AleForceManager::DumpRestartAleForceOfBodyContainer(fstream &file)
{
    for (int iMovingBody = 0; iMovingBody < numberOfMovingBodies; ++ iMovingBody)
    {
        AleForceOfBody *aleForceOfBody = this->GetAleForceOfBody(iMovingBody);
        aleForceOfBody->DumpRestartAleForceOfBodyContainer(file);
    }
}

BasicAleForce *AleForceManager::GetBasicAleForce(int bodyIndex, int iTimeLayer)
{
    AleForceOfBody *aleForceOfBody = this->GetAleForceOfBody(bodyIndex);
    return aleForceOfBody->GetBasicAleForce(iTimeLayer);
}

void CopyBasicAleForce(BasicAleForce *basicAleForce, BasicAleForce *basicAleForceIn)
{
    *basicAleForce = (*basicAleForceIn);
}

AerodynamicForce::AerodynamicForce()
{
    bodyName = "AllWall";
    this->Zero();
}

AerodynamicForce::~AerodynamicForce()
{
}

void AerodynamicForce::Zero()
{
    forceX = 0.0;
    forceY = 0.0;
    forceZ = 0.0;

    forceXByPressure = 0.0;
    forceYByPressure = 0.0;
    forceZByPressure = 0.0;

    forceXByViscosity = 0.0;
    forceYByViscosity = 0.0;
    forceZByViscosity = 0.0;

    momentX = 0.0;
    momentY = 0.0;
    momentZ = 0.0;

    momentXByPressure = 0.0;
    momentYByPressure = 0.0;
    momentZByPressure = 0.0;

    momentXByViscosity = 0.0;
    momentYByViscosity = 0.0;
    momentZByViscosity = 0.0;

    coefForceX = 0.0;
    coefForceY = 0.0;
    coefForceZ = 0.0;

    coefForceXByPressure = 0.0;
    coefForceYByPressure = 0.0;
    coefForceZByPressure = 0.0;

    coefForceXByViscosity = 0.0;
    coefForceYByViscosity = 0.0;
    coefForceZByViscosity = 0.0;

    coefMomentX = 0.0;
    coefMomentY = 0.0;
    coefMomentZ = 0.0;

    coefMomentXByPressure = 0.0;
    coefMomentYByPressure = 0.0;
    coefMomentZByPressure = 0.0;

    coefMomentXByViscosity = 0.0;
    coefMomentYByViscosity = 0.0;
    coefMomentZByViscosity = 0.0;
}

void AerodynamicForce::Initialize(RDouble angleOfAttackDenotedByRadian, RDouble angleOfSideslipDenotedByRadian, RDouble referenceAreaOfAerodynamicForce, RDouble referenceLengthOfAerodynamicMoment, RDouble referenceSpainLengthOfAerodynamicMoment)
{
    this->angleOfAttackDenotedByRadian = angleOfAttackDenotedByRadian;
    this->angleOfSideslipDenotedByRadian = angleOfSideslipDenotedByRadian;
    this->referenceArea = referenceAreaOfAerodynamicForce;
    this->referenceLength = referenceLengthOfAerodynamicMoment;
    this->referenceSpainLength = referenceSpainLengthOfAerodynamicMoment;

    this->Zero();
}

void AerodynamicForce::SetTorqueRef(RDouble torqueRefX, RDouble torqueRefY, RDouble torqueRefZ)
{
    this->torqueRefX = torqueRefX;
    this->torqueRefY = torqueRefY;
    this->torqueRefZ = torqueRefZ;
};

void AerodynamicForce::SetHingeMomentRef(RDouble *localCoordAxis0, RDouble *localCoordAxis1)
{
    for (int m = 0; m < 3; ++m)
    {
        this->localCoordAxis0[m] = localCoordAxis0[m];
        this->localCoordAxis1[m] = localCoordAxis1[m];
    }
};

void AerodynamicForce::Add(DataContainer *dataContainer)
{
    const int numberOfForceVariables = 18;
    RDouble forceMomentContainer[numberOfForceVariables];

    dataContainer->moveto_begin();

    PHSPACE::PHRead(dataContainer, forceMomentContainer, numberOfForceVariables);

    forceX += forceMomentContainer[0];
    forceY += forceMomentContainer[1];
    forceZ += forceMomentContainer[2];

    forceXByPressure += forceMomentContainer[3];
    forceYByPressure += forceMomentContainer[4];
    forceZByPressure += forceMomentContainer[5];

    forceXByViscosity += forceMomentContainer[6];
    forceYByViscosity += forceMomentContainer[7];
    forceZByViscosity += forceMomentContainer[8];

    momentX += forceMomentContainer[9];
    momentY += forceMomentContainer[10];
    momentZ += forceMomentContainer[11];

    momentXByPressure += forceMomentContainer[12];
    momentYByPressure += forceMomentContainer[13];
    momentZByPressure += forceMomentContainer[14];

    momentXByViscosity += forceMomentContainer[15];
    momentYByViscosity += forceMomentContainer[16];
    momentZByViscosity += forceMomentContainer[17];
}

void AerodynamicForce::GetMomentByOriginOfCoordinate(RDouble &aerodynamicMomentX0, RDouble &aerodynamicMomentY0, RDouble &aerodynamicMomentZ0)
{
    RDouble momentReferencePoint[3];

    momentReferencePoint[0] = GlobalDataBase::GetDoubleParaFromDB("TorqueRefX");
    momentReferencePoint[1] = GlobalDataBase::GetDoubleParaFromDB("TorqueRefY");
    momentReferencePoint[2] = GlobalDataBase::GetDoubleParaFromDB("TorqueRefZ");

    RDouble aerodynamicForce[3];

    aerodynamicForce[0] = forceX;
    aerodynamicForce[1] = forceY;
    aerodynamicForce[2] = forceZ;

    RDouble changeOfAerodynamicMoment[3];
    PHSPACE::CrossProduct(momentReferencePoint, aerodynamicForce, changeOfAerodynamicMoment);

    aerodynamicMomentX0 = momentX + changeOfAerodynamicMoment[0];
    aerodynamicMomentY0 = momentY + changeOfAerodynamicMoment[1];
    aerodynamicMomentZ0 = momentZ + changeOfAerodynamicMoment[2];
}

void AerodynamicForce::DumpCoefficient(ActionKey *actionKey)
{
    using namespace PHMPI;
    int myID = GetCurrentProcessorID();
    int server = GetServerProcessorID();
    if (myID != server) return;

    PHSPACE::ParallelOpenFile(actionKey);

    fstream &file = *actionKey->GetFile();

    this->DumpCoefficient(file);

    PHSPACE::ParallelCloseFile(file);
}

void AerodynamicForce::DumpComponentCoefficient(ActionKey *actionKey, int iComponent)
{
    using namespace PHMPI;
    int myID = GetCurrentProcessorID();
    int server = GetServerProcessorID();
    if (myID != server) return;

    int nTotalBC = GlobalBoundaryCondition::GetNumberOfBC();
    for (int iPart = 0; iPart < nTotalBC; ++iPart)
    {
        SimpleBC *boundaryCondition = GlobalBoundaryCondition::GetBC(iPart);
        string partName = boundaryCondition->GetBodyName();
        if (FaceInAerodynamicForceBodys(partName, iComponent))
        {
            bodyName = partName;
        }
    }

    string totalForceFileName = actionKey->GetFileName();
    string partAircoefFileName = PHSPACE::AddSymbolToFileName(totalForceFileName, "_", bodyName);

    fstream file;
    file.open(partAircoefFileName.c_str(), ios_base::out | ios_base::app);

    this->DumpCoefficient(file);

    file.close();
    file.clear();
}

void AerodynamicForce::DumpCoefficient(fstream &file)
{
    this->ComputeCoefficient();

    ostringstream oss;

    int outerIterationSteps = GlobalDataBase::GetIntParaFromDB("outnstep");

    oss << setiosflags(ios::right);
    oss << setiosflags(ios::scientific);

    if (PHSPACE::IfFileEmpty(file))
    {
        vector<string> titleOfTecplot;
        titleOfTecplot.push_back("Title=\"THE AIRCOEF\"");
        titleOfTecplot.push_back("Variables=");
        titleOfTecplot.push_back("iter");
        titleOfTecplot.push_back("<i>C<sub>L</sub></i>");
        titleOfTecplot.push_back("<i>C<sub>D</sub></i>");
        titleOfTecplot.push_back("<i>C<sub>Z</sub></i>");
        titleOfTecplot.push_back("<i>C<sub>D_p</sub></i>");
        titleOfTecplot.push_back("<i>C<sub>D_f</sub></i>");
        titleOfTecplot.push_back("<i>C<sub>D</sub></i>-<i>C<sub>L</sub></i><sup>2</sup>/(<greek>p</greek><math>4</math>A<sub>r</sub>)</i>");
        titleOfTecplot.push_back("<i>X<sub>cp</sub></i>");
        titleOfTecplot.push_back("<i>C<sub>A</sub></i>");
        titleOfTecplot.push_back("<i>C<sub>N</sub></i>");
        titleOfTecplot.push_back("<i>C<sub>Z1</sub></i>");
        titleOfTecplot.push_back("<i>C<sub>ml</sub></i>");
        titleOfTecplot.push_back("<i>C<sub>mn</sub></i>");
        titleOfTecplot.push_back("<i>C<sub>m</sub></i>");
        if (bodyName == "AllWall")
        {
            titleOfTecplot.push_back("\"WallTime\"");
        }

        for (std::size_t iTitle = 0; iTitle < titleOfTecplot.size(); ++iTitle)
        {
            oss << titleOfTecplot[iTitle] << "\n";
        }
    }

    oss << outerIterationSteps << "\t";

    oss << setprecision(10);

    oss << totalLiftCoefficient << "\t";
    oss << totalDragCoefficient << "\t";
    oss << totalCrossCoefficient << "\t";
    oss << totalPressureDragCoefficient << "\t";
    oss << totalViscosityDragCoefficient << "\t";
    oss << DragCoefLiftCoef2piar << "\t";
    oss << pressureCenter << "\t";

    oss << coefForceX << "\t";
    oss << coefForceY << "\t";
    oss << coefForceZ << "\t";
    oss << coefMomentX << "\t";
    oss << coefMomentY << "\t";
    oss << coefMomentZ << "\t";

    if (bodyName == "AllWall")
    {
        RDouble wallTime = PHSPACE::TIME_SPACE::GetWallTime();
        oss << wallTime;
    }

    oss << "\n";
    PHSPACE::WriteASCIIFile(file, oss.str());
}

void AerodynamicForce::ComputeCoefficient()
{
    RDouble sina = sin(angleOfAttackDenotedByRadian);
    RDouble cosa = cos(angleOfAttackDenotedByRadian);
    RDouble sinb = sin(angleOfSideslipDenotedByRadian);
    RDouble cosb = cos(angleOfSideslipDenotedByRadian);

    const RDouble AR = 9.5;

    RDouble referenceCoefficientOfForce = 1.0 / referenceArea;
    RDouble referenceCoefficientOfMoment = 1.0 / (referenceArea * referenceLength);
    RDouble referenceCoefficientOfMomentSpain = 1.0 / (referenceArea * referenceSpainLength);

    coefForceX = forceX * referenceCoefficientOfForce;
    coefForceY = forceY * referenceCoefficientOfForce;
    coefForceZ = forceZ * referenceCoefficientOfForce;

    coefForceXByPressure = forceXByPressure * referenceCoefficientOfForce;
    coefForceYByPressure = forceYByPressure * referenceCoefficientOfForce;
    coefForceZByPressure = forceZByPressure * referenceCoefficientOfForce;

    coefForceXByViscosity = forceXByViscosity * referenceCoefficientOfForce;
    coefForceYByViscosity = forceYByViscosity * referenceCoefficientOfForce;
    coefForceZByViscosity = forceZByViscosity * referenceCoefficientOfForce;

    coefMomentX = momentX * referenceCoefficientOfMomentSpain;
    coefMomentY = momentY * referenceCoefficientOfMomentSpain;
    coefMomentZ = momentZ * referenceCoefficientOfMoment;

    coefMomentXByPressure = momentXByPressure * referenceCoefficientOfMomentSpain;
    coefMomentYByPressure = momentYByPressure * referenceCoefficientOfMomentSpain;
    coefMomentZByPressure = momentZByPressure * referenceCoefficientOfMoment;

    coefMomentXByViscosity = momentXByViscosity * referenceCoefficientOfMomentSpain;
    coefMomentYByViscosity = momentYByViscosity * referenceCoefficientOfMomentSpain;
    coefMomentZByViscosity = momentZByViscosity * referenceCoefficientOfMoment;

    totalLiftCoefficient = -sina * coefForceX + cosa * coefForceY;
    totalDragCoefficient = cosa * cosb * coefForceX + sina * cosb * coefForceY + sinb * coefForceZ;
    totalCrossCoefficient = -cosa * sinb * coefForceX - sina * sinb * coefForceY + cosb * coefForceZ;
    totalPressureDragCoefficient = cosa * cosb * coefForceXByPressure + sina * cosb * coefForceYByPressure + sinb * coefForceZByPressure;
    totalViscosityDragCoefficient = cosa * cosb * coefForceXByViscosity + sina * cosb * coefForceYByViscosity + sinb * coefForceZByViscosity;
    DragCoefLiftCoef2piar = totalDragCoefficient - totalLiftCoefficient * totalLiftCoefficient / (PI * AR);
    pressureCenter = PHSPACE::SIGN(1.0, coefForceY) * coefMomentZ / (PHSPACE::ABS(coefForceY) + SMALL);
}

void AerodynamicForce::CollectOfAllProcessors()
{
    using namespace PHMPI;
    RDouble forceXSum = forceX;
    RDouble forceYSum = forceY;
    RDouble forceZSum = forceZ;

    RDouble forceXByPressureSum = forceXByPressure;
    RDouble forceYByPressureSum = forceYByPressure;
    RDouble forceZByPressureSum = forceZByPressure;

    RDouble forceXByViscositySum = forceXByViscosity;
    RDouble forceYByViscositySum = forceYByViscosity;
    RDouble forceZByViscositySum = forceZByViscosity;

    RDouble momentXSum = momentX;
    RDouble momentYSum = momentY;
    RDouble momentZSum = momentZ;

    RDouble momentXByPressureSum = momentXByPressure;
    RDouble momentYByPressureSum = momentYByPressure;
    RDouble momentZByPressureSum = momentZByPressure;

    RDouble momentXByViscositySum = momentXByViscosity;
    RDouble momentYByViscositySum = momentYByViscosity;
    RDouble momentZByViscositySum = momentZByViscosity;

    PH_Reduce(&forceX, &forceXSum, 1, PH_SUM);
    PH_Reduce(&forceY, &forceYSum, 1, PH_SUM);
    PH_Reduce(&forceZ, &forceZSum, 1, PH_SUM);

    PH_Reduce(&forceXByPressure, &forceXByPressureSum, 1, PH_SUM);
    PH_Reduce(&forceYByPressure, &forceYByPressureSum, 1, PH_SUM);
    PH_Reduce(&forceZByPressure, &forceZByPressureSum, 1, PH_SUM);

    PH_Reduce(&forceXByViscosity, &forceXByViscositySum, 1, PH_SUM);
    PH_Reduce(&forceYByViscosity, &forceYByViscositySum, 1, PH_SUM);
    PH_Reduce(&forceZByViscosity, &forceZByViscositySum, 1, PH_SUM);

    PH_Reduce(&momentX, &momentXSum, 1, PH_SUM);
    PH_Reduce(&momentY, &momentYSum, 1, PH_SUM);
    PH_Reduce(&momentZ, &momentZSum, 1, PH_SUM);

    PH_Reduce(&momentXByPressure, &momentXByPressureSum, 1, PH_SUM);
    PH_Reduce(&momentYByPressure, &momentYByPressureSum, 1, PH_SUM);
    PH_Reduce(&momentZByPressure, &momentZByPressureSum, 1, PH_SUM);

    PH_Reduce(&momentXByViscosity, &momentXByViscositySum, 1, PH_SUM);
    PH_Reduce(&momentYByViscosity, &momentYByViscositySum, 1, PH_SUM);
    PH_Reduce(&momentZByViscosity, &momentZByViscositySum, 1, PH_SUM);

    forceX = forceXSum;
    forceY = forceYSum;
    forceZ = forceZSum;

    forceXByPressure = forceXByPressureSum;
    forceYByPressure = forceYByPressureSum;
    forceZByPressure = forceZByPressureSum;

    forceXByViscosity = forceXByViscositySum;
    forceYByViscosity = forceYByViscositySum;
    forceZByViscosity = forceZByViscositySum;

    momentX = momentXSum;
    momentY = momentYSum;
    momentZ = momentZSum;

    momentXByPressure = momentXByPressureSum;
    momentYByPressure = momentYByPressureSum;
    momentZByPressure = momentZByPressureSum;

    momentXByViscosity = momentXByViscositySum;
    momentYByViscosity = momentYByViscositySum;
    momentZByViscosity = momentZByViscositySum;
}

AerodynamicForceManager::AerodynamicForceManager()
{
    isAllocateMemory = false;
    numberOfAerodynamicForceBodys = 0;
    bodyID = 0;
}

AerodynamicForceManager::~AerodynamicForceManager()
{
    this->DeallocateMemory();
}

void AerodynamicForceManager::AllocateMemory()
{
    if (isAllocateMemory) return;

    isAllocateMemory = true;

    aerodynamicForceContainer.resize(numberOfAerodynamicForceBodys + 1);

    for (int iMovingBody = 0; iMovingBody <= numberOfAerodynamicForceBodys; ++ iMovingBody)
    {
        aerodynamicForceContainer[iMovingBody] = new AerodynamicForce();
    }
}

void AerodynamicForceManager::DeallocateMemory()
{
    if (!isAllocateMemory) return;

    for (int iMovingBody = 0; iMovingBody <= numberOfAerodynamicForceBodys; ++ iMovingBody)
    {
        delete aerodynamicForceContainer[iMovingBody];
    }
}

void AerodynamicForceManager::DumpAerodynamicForceCoefficient(ActionKey *actionKey)
{
    aerodynamicForceContainer[numberOfAerodynamicForceBodys]->DumpCoefficient(actionKey);

    for (int iMovingBody = 0; iMovingBody < numberOfAerodynamicForceBodys; ++ iMovingBody)
    {
        aerodynamicForceContainer[iMovingBody]->DumpComponentCoefficient(actionKey, iMovingBody);
    }

    DumpAllForceData(actionKey);
}

void AerodynamicForceManager::DumpAllForceData(ActionKey *actkey)
{
    string fileName = actkey->filename;

    string mainName, fileDir, extensionName;
    string globalFileName, localFileName;
    GetFileNameExtension(fileName, mainName, extensionName, ".");
    GetFileNameExtension(mainName, fileDir, extensionName, "/");

    ostringstream ossGlobalFileName, ossLocalFileName, ossGlobalData, ossLocalData;
    ossGlobalFileName << fileDir << "/force_global.dat";
    ossLocalFileName << fileDir << "/force_local.dat";
    globalFileName = ossGlobalFileName.str();
    localFileName = ossLocalFileName.str();

    string gridFile = GlobalDataBase::GetStrParaFromDB("gridfile");
    RDouble refMachNumber = GlobalDataBase::GetDoubleParaFromDB("refMachNumber");
    RDouble attackd = GlobalDataBase::GetDoubleParaFromDB("attackd");
    RDouble angleSlide = GlobalDataBase::GetDoubleParaFromDB("angleSlide");
    RDouble refReNumber = GlobalDataBase::GetDoubleParaFromDB("refReNumber");
    RDouble referenceArea = GlobalDataBase::GetDoubleParaFromDB("forceReferenceArea");
    RDouble referenceLength = GlobalDataBase::GetDoubleParaFromDB("forceReferenceLength");
    RDouble referenceSpainLength = GlobalDataBase::GetDoubleParaFromDB("forceReferenceLengthSpanWise");
    RDouble TorqueRefX = GlobalDataBase::GetDoubleParaFromDB("TorqueRefX");
    RDouble TorqueRefY = GlobalDataBase::GetDoubleParaFromDB("TorqueRefY");
    RDouble TorqueRefZ = GlobalDataBase::GetDoubleParaFromDB("TorqueRefZ");
    RDouble gridScaleFactor = GlobalDataBase::GetDoubleParaFromDB("gridScaleFactor");
    RDouble wallTemperature = GlobalDataBase::GetDoubleParaFromDB("wallTemperature");
    RDouble reynoldsReferenceLengthDimensional = GlobalDataBase::GetDoubleParaFromDB("reynoldsReferenceLengthDimensional");
    int outerIterationSteps = GlobalDataBase::GetIntParaFromDB("outnstep");

    RDouble height = 0.0;
    int inflowParaType = GlobalDataBase::GetIntParaFromDB("inflowParaType");
    if (inflowParaType == 1)
    {
        height = GlobalDataBase::GetDoubleParaFromDB("height");
    }

    //! Write global force data file title
    ossGlobalData << "#     MeshFile                     :       " << gridFile << endl;
    if (inflowParaType == 1)
    {
        ossGlobalData << "#     Height                       :       " << height << "  km" << endl;
    }

    ossGlobalData << setiosflags(ios::scientific);
    ossGlobalData << setprecision(6);
    ossGlobalData << "#     Reynolds                     :       " << refReNumber / reynoldsReferenceLengthDimensional << "  1/m" << endl;
    ossGlobalData << resetiosflags(ios::scientific);

    ossGlobalData << "#     Mach                         :       " << refMachNumber << endl;
    ossGlobalData << "#     Alpha                        :       " << attackd << "  degree" << endl;
    ossGlobalData << "#     Beta                         :       " << angleSlide << "  degree" << endl;

    if (wallTemperature < 0.0)
    {
        ossGlobalData << "#     Wall Temperature             :       " << "Adiabatic Wall" << endl;
    }
    else
    {
        ossGlobalData << "#     Wall Temperature             :       " << wallTemperature << "  K" << endl;
    }

    ossGlobalData << "#     gridScaleFactor              :       " << gridScaleFactor << endl;
    ossGlobalData << "#     forceReferenceLength         :       " << referenceLength * reynoldsReferenceLengthDimensional << "  m" << endl;
    ossGlobalData << "#     forceReferenceLengthSpanWise :       " << referenceSpainLength * reynoldsReferenceLengthDimensional << "  m" << endl;

    if (GetDim() == THREE_D)
    {
        ossGlobalData << "#     forceReferenceArea           :       " << referenceArea * reynoldsReferenceLengthDimensional * reynoldsReferenceLengthDimensional << "  m2" << endl;
    }
    else
    {
        ossGlobalData << "#     forceReferenceArea           :       " << referenceArea * reynoldsReferenceLengthDimensional << "  m^2" << endl;
    }

    ossGlobalData << "#     TorqueRefX                   :       " << TorqueRefX * reynoldsReferenceLengthDimensional << endl;
    ossGlobalData << "#     TorqueRefY                   :       " << TorqueRefY * reynoldsReferenceLengthDimensional << endl;
    ossGlobalData << "#     TorqueRefZ                   :       " << TorqueRefZ * reynoldsReferenceLengthDimensional << endl;
    ossGlobalData << "#     iter                         :       " << outerIterationSteps << endl << endl;
    ossGlobalData << "#        CA             CN            CZ1           Cml             Cmn            Cm" << endl;
    ossGlobalData << "#        CA_p           CN_p          CZ1_p         Cml_p           Cmn_p          Cm_p" << endl;
    ossGlobalData << "#        CA_f           CN_f          CZ1_f         Cml_f           Cmn_f          Cm_f" << endl;

    //! Write local force data file title
    ossLocalData << "#     MeshFile             :       " << gridFile << endl;
    if (inflowParaType == 1)
    {
        ossLocalData << "#     Height               :       " << height << "  km" << endl;
    }

    ossLocalData << setiosflags(ios::scientific);
    ossLocalData << setprecision(6);
    ossLocalData << "#     Reynolds             :       " << refReNumber / reynoldsReferenceLengthDimensional << "  1/m" << endl;
    ossLocalData << resetiosflags(ios::scientific);

    ossLocalData << "#     Mach                 :       " << refMachNumber << endl;
    ossLocalData << "#     Alpha                :       " << attackd << "  degree" << endl;
    ossLocalData << "#     Beta                 :       " << angleSlide << "  degree" << endl;

    if (wallTemperature < 0.0)
    {
        ossLocalData << "#     Wall Temperature     :       " << "Adiabatic Wall" << endl;
    }
    else
    {
        ossLocalData << "#     Wall Temperature     :       " << wallTemperature << "  K" << endl;
    }

    ossLocalData << "#     gridScaleFactor      :       " << gridScaleFactor << endl;
    ossLocalData << "#     iter                 :       " << outerIterationSteps << endl << endl;
    ossLocalData << "#        CA             CN            CZ1           Cml             Cmn            Cm" << endl;
    ossLocalData << "#        CA_p           CN_p          CZ1_p         Cml_p           Cmn_p          Cm_p" << endl;
    ossLocalData << "#        CA_f           CN_f          CZ1_f         Cml_f           Cmn_f          Cm_f" << endl;
    ossLocalData << "#        CMh" << endl;

    string partName = "whole";
    aerodynamicForceContainer[numberOfAerodynamicForceBodys]->DumpAllGlobalForceData(partName, numberOfAerodynamicForceBodys, ossGlobalData);

    //! Write force data file
    for (int iMovingBody = 0; iMovingBody <= numberOfAerodynamicForceBodys; ++ iMovingBody)
    {

        const string &bodyName = aerodynamicForceContainer[numberOfAerodynamicForceBodys]->GetBodyName();
        partName = bodyName;

        aerodynamicForceContainer[iMovingBody]->DumpAllGlobalForceData(partName, iMovingBody, ossGlobalData);
        aerodynamicForceContainer[iMovingBody]->DumpAllLocalForceData(partName, iMovingBody, ossLocalData);

        //! Write global force data file
        fstream globalFile;
        ParallelOpenFile(globalFile, globalFileName, ios_base::out | ios_base::trunc);
        WriteASCIIFile(globalFile, ossGlobalData.str());
        ParallelCloseFile(globalFile);

        //! Write local force data file
        fstream localFile;
        ParallelOpenFile(localFile, localFileName, ios_base::out | ios_base::trunc);
        WriteASCIIFile(localFile, ossLocalData.str());
        ParallelCloseFile(localFile);
    }
}

void AerodynamicForce::DumpAllLocalForceData(const string &partName, int partID, std::ostringstream &oss)
{
    RDouble reynoldsReferenceLengthDimensional = GlobalDataBase::GetDoubleParaFromDB("reynoldsReferenceLengthDimensional");

    oss << resetiosflags(ios::scientific);

    oss << endl;
    oss << "*" << partName << endl;
    oss << "#forceReferenceLength         :    " << referenceLength * reynoldsReferenceLengthDimensional << "  m" << endl;
    oss << "#forceReferenceLengthSpanWise :    " << referenceSpainLength * reynoldsReferenceLengthDimensional << "  m" << endl;

    if (GetDim() == THREE_D)
    {
        oss << "#forceReferenceArea           :    " << referenceArea * reynoldsReferenceLengthDimensional * reynoldsReferenceLengthDimensional << "  m2" << endl;
    }
    else
    {
        oss << "#forceReferenceArea           :    " << referenceArea * reynoldsReferenceLengthDimensional << "  m^2" << endl;
    }

    oss << "#TorqueRefX                   :    " << torqueRefX * reynoldsReferenceLengthDimensional << endl;
    oss << "#TorqueRefY                   :    " << torqueRefY * reynoldsReferenceLengthDimensional << endl;
    oss << "#TorqueRefZ                   :    " << torqueRefZ * reynoldsReferenceLengthDimensional << endl;

    oss << setiosflags(ios::right);
    oss << setiosflags(ios::scientific);
    oss << setprecision(10);

    oss << coefForceX << "	";
    oss << coefForceY << "	";
    oss << coefForceZ << "	";
    oss << coefMomentX << "	";
    oss << coefMomentY << "	";
    oss << coefMomentZ << endl;

    oss << coefForceXByPressure << "	";
    oss << coefForceYByPressure << "	";
    oss << coefForceZByPressure << "	";
    oss << coefMomentXByPressure << "	";
    oss << coefMomentYByPressure << "	";
    oss << coefMomentZByPressure << endl;

    oss << coefForceXByViscosity << "	";
    oss << coefForceYByViscosity << "	";
    oss << coefForceZByViscosity << "	";
    oss << coefMomentXByViscosity << "	";
    oss << coefMomentYByViscosity << "	";
    oss << coefMomentZByViscosity << endl;
}

void AerodynamicForce::DumpAllGlobalForceData(const string &partName, int partID, std::ostringstream &oss)
{
    oss << setiosflags(ios::right);
    oss << setiosflags(ios::scientific);
    oss << setprecision(10);

    oss << endl;
    oss << "*" << partName << endl;

    oss << coefForceX << "	";
    oss << coefForceY << "	";
    oss << coefForceZ << "	";
    oss << coefMomentX << "	";
    oss << coefMomentY << "	";
    oss << coefMomentZ << endl;

    oss << coefForceXByPressure << "	";
    oss << coefForceYByPressure << "	";
    oss << coefForceZByPressure << "	";
    oss << coefMomentXByPressure << "	";
    oss << coefMomentYByPressure << "	";
    oss << coefMomentZByPressure << endl;

    oss << coefForceXByViscosity << "	";
    oss << coefForceYByViscosity << "	";
    oss << coefForceZByViscosity << "	";
    oss << coefMomentXByViscosity << "	";
    oss << coefMomentYByViscosity << "	";
    oss << coefMomentZByViscosity << endl;

}

void AerodynamicForceManager::InitializeAerodynamicForceManager()
{
    numberOfAerodynamicForceBodys = GlobalDataBase::GetIntParaFromDB("numberOfMovingBodies");
    this->AllocateMemory();

    RDouble angleOfAttackDenotedByRadian = GlobalDataBase::GetDoubleParaFromDB("attack");
    RDouble angleOfSideslipDenotedByRadian = GlobalDataBase::GetDoubleParaFromDB("sideslip");

    RDouble referenceAreaOfAerodynamicForce = GlobalDataBase::GetDoubleParaFromDB("forceReferenceArea");
    RDouble referenceLengthOfAerodynamicMoment = GlobalDataBase::GetDoubleParaFromDB("forceReferenceLength");
    RDouble referenceSpainLengthOfAerodynamicMoment = GlobalDataBase::GetDoubleParaFromDB("forceReferenceLengthSpanWise");

    for (int iMovingBody = 0; iMovingBody <= numberOfAerodynamicForceBodys; ++ iMovingBody)
    {
        AerodynamicForce *aerodynamicForce = this->GetAerodynamicForce(iMovingBody);
        aerodynamicForce->Initialize(angleOfAttackDenotedByRadian, angleOfSideslipDenotedByRadian, referenceAreaOfAerodynamicForce,
            referenceLengthOfAerodynamicMoment, referenceSpainLengthOfAerodynamicMoment);
    }
}

void AerodynamicForceManager::ExtractForceDataFromActionKey(ActionKey *actionKey, vector< DataContainer * > *dataContainerArray)
{
    DataContainer *dataContainer = actionKey->GetData();
    dataContainer->moveto_begin();

    const int numberOfData = 18;
    RDouble bufferForSwap[numberOfData];

    for (int iMovingBody = 0; iMovingBody <= numberOfAerodynamicForceBodys; ++ iMovingBody)
    {
        PHSPACE::PHRead(dataContainer, bufferForSwap, numberOfData);
        PHSPACE::PHWrite((*dataContainerArray)[iMovingBody], bufferForSwap, numberOfData);
    }
}

vector< DataContainer * > *AerodynamicForceManager::CreateDataContainerArray()
{
    vector< DataContainer * > *dataContainerArray = new vector< DataContainer * >(numberOfAerodynamicForceBodys + 1);
    for (int iMovingBody = 0; iMovingBody <= numberOfAerodynamicForceBodys; ++ iMovingBody)
    {
        (*dataContainerArray)[iMovingBody] = new DataContainer();
    }
    return dataContainerArray;
}

void AerodynamicForceManager::FreeDataContainerArray(vector< DataContainer * > *dataContainerArray)
{
    for (int iMovingBody = 0; iMovingBody <= numberOfAerodynamicForceBodys; ++ iMovingBody)
    {
        delete (*dataContainerArray)[iMovingBody];
    }
    delete dataContainerArray;
}

void AerodynamicForceManager::AddAerodynamicForce(ActionKey *actionKey)
{
    vector< DataContainer * > *dataContainerArray = this->CreateDataContainerArray();

    this->ExtractForceDataFromActionKey(actionKey, dataContainerArray);

    for (int iMovingBody = 0; iMovingBody < numberOfAerodynamicForceBodys; ++ iMovingBody)
    {
        DataContainer *dataContainer = (*dataContainerArray)[iMovingBody];
        aerodynamicForceContainer[iMovingBody]->Add(dataContainer);
        aerodynamicForceContainer[numberOfAerodynamicForceBodys]->Add(dataContainer);
    }
    this->FreeDataContainerArray(dataContainerArray);
}

void AerodynamicForceManager::CollectAerodynamicForceOfAllProcessors()
{
    for (int iMovingBody = 0; iMovingBody <= numberOfAerodynamicForceBodys; ++ iMovingBody)
    {
        aerodynamicForceContainer[iMovingBody]->CollectOfAllProcessors();
    }
}

AerodynamicForceManager *GetAerodynamicForceManager()
{
    return &aerodynamicForceManager;
}

AerodynamicForce *GetCurrentAerodynamicForce(int bodyIndex)
{
    return GetAerodynamicForceManager()->GetAerodynamicForce(bodyIndex);
}

BasicAerodynamicForce *GetBasicAerodynamicForce(int iBoundaryConditionRegion, BasicAerodynamicForce **basicAerodynamicForceContainer)
{
    return basicAerodynamicForceContainer[iBoundaryConditionRegion];
}

BasicAerodynamicForce *GetBasicAerodynamicForce(int iBoundaryConditionRegion, vector<BasicAerodynamicForce *> *basicAerodynamicForceVector)
{
    return (*basicAerodynamicForceVector)[iBoundaryConditionRegion];
}

bool OutPutSubIterationForceHeader()
{
    int outerIterationSteps = GlobalDataBase::GetIntParaFromDB("outnstep");
    int innerIterationSteps = GlobalDataBase::GetIntParaFromDB("innstep");
    int numberOfForceSaveSteps = GlobalDataBase::GetIntParaFromDB("intervalStepForce");

    if (outerIterationSteps / numberOfForceSaveSteps == 1 && innerIterationSteps / numberOfForceSaveSteps == 1)
    {
        return true;
    }

    return false;
}

}

