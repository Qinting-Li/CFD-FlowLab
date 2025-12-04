#include "IndexParticleOfCell.h"
#include "TK_Exit.h"
#include "ParticleBoundaryCondition.h"

namespace PHSPACE
{
IndexParticleOfCell::IndexParticleOfCell()
{
    this->particleCellID = new set<int>;
}

IndexParticleOfCell::IndexParticleOfCell(const IndexParticleOfCell &rhs)
{
    this->particleCellID = new set<int>;
    set<int>::iterator iter;
    for (iter = rhs.particleCellID->begin(); iter != rhs.particleCellID->end(); ++iter)
    {
        this->particleCellID->insert(*iter);
    }
    this->bcIndex.clear();
}

IndexParticleOfCell::~IndexParticleOfCell()
{
    //! Destructors can cause many problems of deletion by mistake.
    //! We must be careful. 
    this->DeleteSet();
}

bool IndexParticleOfCell::IsExistParticleInCell(int particleID)
{
    //! Is similar to Data_Param::IsExist.
    set<int>::iterator iter;
    iter = particleCellID->find(particleID);
    if (iter != particleCellID->end())
    {
        //! Particle is exist on current cell.
        return true;
    }
    else
    {
        //! Particle is not exist on current cell.
        return false;
    }
}

bool IndexParticleOfCell::IsEmpty()
{
    return this->particleCellID->empty();
}

void IndexParticleOfCell::PrintParticleIDInCell()
{
    if (IsEmpty())
    {
        cout << "This cell has no particles" << endl;
    }
    else
    {
        set<int>::iterator iter;
        for (iter = particleCellID->begin(); iter != particleCellID->end(); ++iter)
        {
            cout << "Particle ID on Cell :" << *iter << endl;
        }
    }
}

void IndexParticleOfCell::AddParticleID(int particleID)
{
    set<int>::iterator iter;
    iter = particleCellID->find(particleID);
    if (iter != particleCellID->end())
    {

    }
    else
    {
        this->particleCellID->insert(particleID);
    }
}

void IndexParticleOfCell::RemoveParticleID(int particleID)
{
    set<int>::iterator iter;
    iter = particleCellID->find(particleID);
    if (iter != particleCellID->end())
    {
        this->particleCellID->erase(iter);
    }
    else
    {
        TK_Exit::UnexpectedVarValue("Un exist Particle ID = ", particleID);
    }
}

int IndexParticleOfCell::GetNumOfParticleOfCell()
{
    if (IsEmpty())
    {
        return NULL;
    }
    else
    {
        //! Note the different of size_t and int.
        //! size_t is unsigned integer,which mat be 4 or 8 byte.
        //! int is signed integer,which must be 4 byte.
        return static_cast<int>(this->particleCellID->size());
    }
}

int *IndexParticleOfCell::GetParticleIDOfCell()
{
    if (IsEmpty())
    {
        return NULL;
    }
    else
    {
        int count = 0;
        set<int>::iterator iter;
        int nParticle = GetNumOfParticleOfCell();
        int *particeIDOfCell = new int[nParticle];
        for (iter = particleCellID->begin(); iter != particleCellID->end(); ++iter)
        {
            if (count >= nParticle)
            {
                TK_Exit::UnexpectedVarValue("error on iter, count = ", count);
                break;
            }
            particeIDOfCell[count] = *iter;
            ++count;
        }
        return particeIDOfCell;
    }
}

void IndexParticleOfCell::ClearParticle()
{
    this->particleCellID->clear();
}

void IndexParticleOfCell::DeleteSet()
{
    //! In STL, the vector container is often used to empty and store data for many times. 
    //! Using clear() only empties elements without releasing memory. 
    //! You can use swap() to empty elements and release memory.
    //! But for set/map associative containers, use clear() is ok.
    if (0 == this)
    {
        return;
    }

    if (0 == particleCellID)
    {
        return;
    }
    
    this->particleCellID->clear();
    delete this->particleCellID;
    this->particleCellID = nullptr;
}

void IndexParticleOfCell::InitBCIndex(int nFace)
{
    for (int iFace = 0; iFace < nFace; ++iFace)
    {
        using namespace PARTICLEBCSPACE;
        int initBC = NO_PARTICLE_BC;
        bool check = this->bcIndex.insert(make_pair(iFace, initBC)).second;

        if (!check)
        {
            ostringstream ossError;
            ossError << "Error: insert init cell bc in IndexParticleOfCell: " << endl;
            TK_Exit::ExceptionExit(ossError);
        }
    }
}

void IndexParticleOfCell::SetBCType(int iFace, int iBC)
{
    auto iterBCIndex = bcIndex.find(iFace);
    if (iterBCIndex != bcIndex.end())
    {
        iterBCIndex->second = iBC;
    }
    else
    {
        TK_Exit::UnexpectedVarValue("Un exist iFace ID = ", iFace);
    }
}

void IndexParticleOfCell::SetCellType(int cellType)
{
    this->cellType = cellType;
}

int IndexParticleOfCell::GetCellType()
{
    return this->cellType;
}

map<int, int> *IndexParticleOfCell::GetBCType()
{
    return &(this->bcIndex);
}

int IndexParticleOfCell::GetBCType(int iFace)
{
    using namespace PARTICLEBCSPACE;
    int bcType = NO_PARTICLE_BC;
    map<int, int>::iterator iterBCIndex = bcIndex.find(iFace);
    if (iterBCIndex != bcIndex.end())
    {
        bcType = iterBCIndex->second;
    }
    else
    {
        TK_Exit::UnexpectedVarValue("Un exist iFace ID = ", iFace);
    }
    return bcType;
}

}