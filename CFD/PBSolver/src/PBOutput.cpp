#include "PBOutput.h"

namespace PHSPACE
{
PBOutput::PBOutput(const string &name):PBObject(name)
{

}

PBOutput::~PBOutput()
{

}

void PBOutput::Run(Zone *zone, vector<string> &varNameList)
{

}

DEFINE_PBSOLVER_REGISTER_FUNCTION(PBOutput, PBOutput, OUTPUT_EMPTY);
}