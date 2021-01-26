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
//! @file      JetBC.h
//! @brief     To set boundary information and data structure for Jet boundary conditions.
//! @author    He Kun, Chen Liangzhong.

#pragma once
#include "Math_BasisFunction.h"
#include "Geo_Grid.h"
#include "BC_Struct.h"
#include "BC_Unstruct.h"

using namespace std;

namespace PHSPACE
{
class JetData
{
private:
    string nozzleName;
    int BCID;
    bool isJetOn;
    //! Basic five variables of jet parameters.
    vector <RDouble> velocity;
    RDouble pressure, temperature, density;

    //! Advanced parameters of jet, which will useful later.
    RDouble specificHeatRatio, molecularWeight, machNumber;
    int nChemicalComponent;
    vector <RDouble> chemicalConstituents;

public:
    JetData();
    ~JetData();

public:
    //! To set the nozzle name to a jet.
    void SetJetName(const string &nozzleName);
    //! To set jet on/off flag by in-parameter isJetOn. \n
    //! true - on \n
    //! false - off
    void SetJetOnJudgement(const bool isJetOn);
    //! To set the id number to BCID.
    void SetBCID(const int BCID);
    //! To set velocity by in-parameter.
    void SetVelocity(const vector <RDouble> &velocity);
    //! To set pressure by in-parameter.
    void SetPressure(const RDouble pressure);
    //! To set density by in-parameter.
    void SetDensity(const RDouble density);
    //! To set temperature by in-parameter.
    void SetTemperature(const RDouble temperature);
    //! To set specific heat ratio by in-parameter.
    void SetSpecificHeatRatio(const RDouble specificHeatRatio);
    //! To set molecular weight by in-parameter.
    void SetMolecularWeight(const RDouble molecularWeight);
    //! To set mach number by in-parameter.
    void SetmachNumber(const RDouble machNumber);
    //! To set number of chemical constituents by in-parameter.
    void SetNumberOfChemicalComponents(const int nChemicalComponent);
    //! To set chemical constituents by in-parameter.
    void SetChemicalConstituents(const vector <RDouble> &chemicalConstituents);

    //! To Get the nozzle name.
    string GetJetName();
    //! To Get the nozzle ID.
    int GetBCID();
    //! To Get the nozzle on/off flag.
    bool GetJetOnJudgement();
    //! To Get air velocity on jet boundary.
    vector <RDouble> GetVelocity();
    //! To Get pressure on jet boundary.
    RDouble GetPressure();
    //! To Get density on jet boundary.
    RDouble GetDensity();
    //! To Get temperature on jet boundary.
    RDouble GetTemperature();
    //! To Get specific heat ratio on jet boundary.
    RDouble GetSpecificHeatRatio();
    //! To Get molecular weight on jet boundary.
    RDouble GetMolecularWeight();
    //! To Get mach number on jet boundary.
    RDouble GetMachNumber();
    //! To Get number of chemical constituents by in-parameter.
    int getNumberOfChemicalComponents();
    //! To Get chemical constituents on jet boundary.
    vector <RDouble> GetChemicalConstituents();

};

class JetBCInfo
{
private:
    string jetFileName;
    int BCIDMax, BCIDMin;
    int nJets;
    int hasSimilarBC;
    bool jetFileExist;
    vector <JetData *> jetBCs;

public:
    JetBCInfo();
    ~JetBCInfo();

public:
    //! To Set the flag if the jet boundary condition exist or not.
    void SetSimilarBCExistLabel();
    //! To set the range of jet boundary condition ID according to the jet.hypara file.
    void SetBCIDRange();
    //! To set the file name of jet parameter.
    void SetBCDataFileName(const string &fileName);
    //! To set the file name of jet parameter.
    void CheckIfBCDataFileExist();

    //! To get the Lowest range of jet boundary condition ID.
    int GetBCIDRangLower();
    //! To get the uppest range of jet boundary condition ID.
    int GetBCIDRangUpper();
    //! To get the number of jets.
    int GetNumberOfJets();
    //! To get the file name of jet parameter.
    string GetBCDataFileName();
    //! To get all of the jet boundary condition.
    vector<JetData *> GetBCInfo();
    //! To check if the BCType exist in BCIDs of jetDefineFile.
    void CheckSpecificBCExist();
    //! To check if the jet parameter define file exist.
    bool IfDefineFileExist();
    //! Nodimension the jet parameter
    void Nondimensional();
    //! To get the flag if Jet exist or not.
    int GetJetBCExistLabel();
    //! To get the jet which is in use.
    JetData *GetInUseJets(const int BCType);
    //! To read jet parameter in jetDefineFile.
    void ReadBCData();
    //! Bcast jet parameter to all processor.
    void BcastBCData();
    //! Compress jet parameter to DataContainer.
    void CompressJetData(DataContainer *Jetdata);
    //! Decompress DataContainer into jet parameter.
    void DecompressJetData(DataContainer *Jetdata);

};

//! The controller to set jet boundary informations.
void GetJetBCData();

//! To return true if jet Boundary condition define file exist. /n
//! Return false elsewise.
//void CheckifJetBCFileExist();

//! To check if the input BCType is jet type and if this jet is open./n
//! true - Is jet and open.
//! false - Is not jet or jet is closed.
bool ifIsJetBC(const int BCType);

//! Function to do jet boundary condition.
//! TO DO:Turbulence is not considered here now!
//! Presently, this is only useful where the jet conditions are prescribed, 
//! Set the primitive variables' value on ghost cell according to the jet parameters od r, u, p. \n

//! 1. Srtucture.
//! @param[out]   grid         The grid and flow data which is to  modified or set
//! @param[in]    BCRegion     The boundary which is to specified.
void JetBC3D(Grid *grid, BCRegion *BCRegion); //! Structure solver

//! 2. Unsrtucture.
//! @param[out]    primTarget    The primitive data which is to  modified or set
//! @param[in]     BCType        The boundary type which is to specified.
void JetBC3D(RDouble *primTarget, int BCType); //! Unstructure solver

//! 2. Unsrtucture.
//! JetBC boundary condition for one Region of a bounday.
//! @param[in] gridIn              an unstruct grid pointer.
//! @param[in] bcRegionUnstruct    an object of some BCRegionUnstruct class, contain the information of this BC.
void JetBC3DRegion(Grid *grid, BCRegionUnstruct *bcRegionUnstruct); //! Unstructure solver
}
