#include "RungeKuttaMethod.h"

namespace PHSPACE
{
RungeKuttaIndex::RungeKuttaIndex(int nRK)
{
    this->nRK = nRK;
    this->kindex = new RDouble[nRK];
    this->findex = new RDouble[nRK];

    if (4 == nRK)
    {
        kindex[0] = 1;
        kindex[1] = 2;
        kindex[2] = 2;
        kindex[3] = 1;

        findex[0] = 0.5;
        findex[1] = 0.5;
        findex[2] = 1;
        findex[3] = 1.0/6.0;

    }
    this->index = 0;
}

RungeKuttaIndex::~RungeKuttaIndex()
{
    delete [] kindex;    kindex = nullptr;
    delete [] findex;    findex = nullptr;
    this->index = -1;
}

void RungeKuttaIndex::AddIndex()
{
    this->index++;
}

int RungeKuttaIndex::GetIndex()
{
    return this->index;
}

int RungeKuttaIndex::GetNumRK()
{
    return this->nRK;
}

RDouble RungeKuttaIndex::GetFindex()
{
    return this->findex[index];
}

RDouble RungeKuttaIndex::GetKindex()
{
    return this->kindex[index];
}

void RungeKuttaIndex::ResetIndex()
{
    this->index = 0;
}

}