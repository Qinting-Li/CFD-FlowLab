#include "ParticleBoundaryCondition.h"
#include "TK_Exit.h"
#include "TK_Parse.h"

using namespace std;
namespace PHSPACE
{
vector< vector<string>* > *GlobalParticleBoundaryLink::globalLinkName = new vector< vector<string>* >;

GlobalParticleBoundaryLink::GlobalParticleBoundaryLink()
{

}

GlobalParticleBoundaryLink::~GlobalParticleBoundaryLink()
{

}

void GlobalParticleBoundaryLink::InitGlobalParticleBoundaryLink(int nBC)
{
    for (int iBC = 0; iBC < nBC; iBC++)
    {
        vector<string> *bcLinkName = new vector<string>;
        globalLinkName->push_back(bcLinkName);
    }
}

void GlobalParticleBoundaryLink::AddBCLinkName(int iBC, string bcLinkName)
{
    (*globalLinkName)[iBC]->push_back(bcLinkName);
}

int GlobalParticleBoundaryLink::GetNumOfGlobalParticleBoundaryLink()
{
    return globalLinkName->size();
}

int GlobalParticleBoundaryLink::GetNumOfLinkBCByID(int iBC)
{
    return (*globalLinkName)[iBC]->size();
}

bool GlobalParticleBoundaryLink::IsBCLink(int iBC)
{
    if ((*globalLinkName)[iBC] == 0)
    {
        return false;
    }
    else
    {
        return true;
    }
}

vector< vector<string>* > *GlobalParticleBoundaryLink::GetGlobalLinkName()
{
    return globalLinkName;
}

ParticleBoundaryCondition::ParticleBoundaryCondition()
{
    this->particleBCType = -1;
    this->particleBCParamDataBase = 0;
    this->particleBCLink = new vector<ParticleBoundaryCondition* >;
}

ParticleBoundaryCondition::~ParticleBoundaryCondition()
{
    //! no new here,no delete.
    delete particleBCParamDataBase;
    delete this->particleBCLink;
}

void ParticleBoundaryCondition::SetParticleBCName(string &particleBCName)
{
    this->particleBCName = particleBCName;
}

void ParticleBoundaryCondition::SetParticleBCTye(int &particleBCType)
{
    this->particleBCType = particleBCType;
}

void ParticleBoundaryCondition::SetParticleBCParam(Data_Param *particleBCParamDataBase)
{
    this->particleBCParamDataBase = particleBCParamDataBase;
}

const string &ParticleBoundaryCondition::GetParticleBCName()
{
    return this->particleBCName;
}

int ParticleBoundaryCondition::GetParticleBCType()
{
    return this->particleBCType;
}

Data_Param *ParticleBoundaryCondition::GetParticleBCParam()
{
    return this->particleBCParamDataBase;
}

bool ParticleBoundaryCondition::IsPeriodicBC()
{
    if (this->particleBCType == PARTICLEBCSPACE::PERIODIC)
    {
        return true;
    }
    else
    {
        return false;
    }
}

}