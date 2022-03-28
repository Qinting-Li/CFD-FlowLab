//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//          PPPPP  H   H  EEEEE  N    N  GGGGG  L      EEEEE  III         +
//          P   P  H   H  E      NN   N  G      L      E       I          +
//          PPPPP  HHHHH  EEEEE  N N  N  G  GG  L      EEEEE   I          +
//          P      H   H  E      N  N N  G   G  L      E       I          +
//          P      H   H  EEEEE  N    N  GGGGG  LLLLL  EEEEE  III         +
//------------------------------------------------------------------------+
//          Platform for Hybrid Engineering Simulation of Flows           +
//          China Aerodynamics Research and Development Center            +
//                     (C) Copyright, Since 2010                          +
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      DSMCCollision.h
//! @brief     Collision operator for DSMC.
//! @author    Li Jin.

#pragma once

namespace PHSPACE
{

class Collision 
{
private:   
    RDouble VRR,VR,CVR;
    vector<RDouble> VRC,VCM;
    RDouble ** CCG;
public:
    Collision(int);
    ~Collision();
    RDouble ** GetCCG() const {return CCG;}
    vector<RDouble> GetVRC() const {return VRC;}
    void SetVRC(vector<RDouble> var) {VRC = var;}
    vector<RDouble> GetVCM() const {return VCM;}
    void SetVRR(RDouble var){VRR = var;}
    void SetVR(RDouble var){VR = var;}
    void SetCVR(RDouble var){CVR = var;}

    RDouble GetCVR() const{ return CVR;}
    RDouble GetVRR() const{ return VRR;}
    RDouble GetVR() const{ return VR;}
};

}



