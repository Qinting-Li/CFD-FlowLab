#include "PBTransient.h"

namespace PHSPACE
{
PBTransient::PBTransient(const string &name):PBObject(name)
{
    if (name == "EMPTY")
    {
        if (ObjGetIntParaFromDB("iunsteady", 0))
        {
            cout << "warning" << endl;
        }
    }
}

PBTransient::~PBTransient()
{

}

void PBTransient::AllocateMemory(UnstructGrid *grid, vector<string> &phiName, vector<RDouble> &initPhi)
{

}

void PBTransient::Discrt(UnstructGrid *grid, string phiName, string phiOldName, string phiOldOldName, string rhoOldName,string rhoOldOldName)
{

}

void PBTransient::SaveOldTimeValue(UnstructGrid *grid, vector<string> &phiName)
{

}

DEFINE_PBSOLVER_REGISTER_FUNCTION(PBTransient, PBTransient, TRANSIENT_EMPTY);
}