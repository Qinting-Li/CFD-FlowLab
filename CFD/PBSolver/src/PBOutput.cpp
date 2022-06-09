#include "PBOutput.h"

namespace PHSPACE
{
    PBOutput::PBOutput(const string &name):PBObject(name)
    {

    }

    PBOutput::~PBOutput()
    {

    }

    DEFINE_PBSOLVER_REGISTER_FUNCTION(PBOutput, PBOutput, OUTPUT_EMPTY);
}