#include "ParticleAnalysis.h"
#include "GlobalDataBase.h"
#include "TK_Exit.h"
#include "TK_Parse.h"
#include "Constants.h"

namespace PHSPACE
{
void ParticleAnalysis()
{
    using namespace PARTICLE_ANALYSIS_SPACE;

    int particleTask = GlobalDataBase::GetIntParaFromDB("particleTask");

    switch (particleTask)
    {
    case ADDPARTICLEBCTOGRID:
        AddParticleBC();
        break;
    default:
        TK_Exit::UnexpectedVarValue("particleTask = ", particleTask);
        break;
    }
}

void AddParticleBC()
{

}

namespace PARTICLE_PARAM
{
void ReadParticleParamFile(Data_Param *parameter, const string &paraFileName)
{
    fstream file;
    file.open(paraFileName.c_str(), ios_base::in);

    if (!file)
    {
        TK_Exit::FileOpenErrorExit(paraFileName);
    }

    TK_Parse::ReadBasicData(file, parameter);

    file.close();
    file.clear();
}

string GetParticleStrPara(Data_Param *parameter, const string &name)
{
    string data;
    parameter->GetData(name, &data, PHSTRING, 1);
    return data;
}

int GetParticleIntPara(Data_Param *parameter, const string &name)
{
    int data;
    parameter->GetData(name, &data, PHINT, 1);
    return data;
}

RDouble GetParticleDoublePara(Data_Param *parameter, const string &name)
{
    RDouble data;
    parameter->GetData(name, &data, PHDOUBLE, 1);
    return data;
}

}

}