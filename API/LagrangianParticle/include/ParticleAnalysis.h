//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//          PPPPP  H   H  EEEEE  N    N  GGGGG  L      EEEEE  III         +
//          P   P  H   H  E      NN   N  G      L      E       I          +
//          PPPPP  HHHHH  EEEEE  N N  N  G  GG  L      EEEEE   I          +
//          P      H   H  E      N  N N  G   G  L      E       I          +
//          P      H   H  EEEEE  N    N  GGGGG  LLLLL  EEEEE  III         +
//------------------------------------------------------------------------+
//          Platform for Hybrid Engineering Simulation of Flows           +
//          China Aerodynamics Research and Development Center             +
//                     (C) Copyright, Since 2010                          +
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      ParticleAnalysis.h
//! @brief     It is the controller of ParticleAnalysis, which include the
//!            tool function of particle,as add particle boundary condtion 
//!            to *.fts, the function of Data_Param.
//! @author    Lei Yinghaonan (Lanzhou University).

#include "Data_Param.h"
#include "Precision.h"

namespace PHSPACE
{

namespace PARTICLE_ANALYSIS_SPACE
{
const int ADDPARTICLEBCTOGRID = 1;
}

void ParticleAnalysis();

void AddParticleBC();

namespace PARTICLE_PARAM
{
    //! DataParam IO.

    //! Read Data_Param file.
    void ReadParticleParamFile(Data_Param *parameter, const string &paraFileName);

    //! Get param from Data_Param.
    string GetParticleStrPara(Data_Param *parameter, const string &name);
    int GetParticleIntPara(Data_Param *parameter, const string &name);
    RDouble GetParticleDoublePara(Data_Param *parameter, const string &name);

}

}