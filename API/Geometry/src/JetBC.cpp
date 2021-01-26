#include "JetBC.h"
#include "GridType.h"
#include "BC_Struct.h"
#include "BC_Unstruct.h"
#include "Geo_StructGrid.h"
#include "Geo_MultiGridInfo_Struct.h"
#include "Geo_UnstructGrid.h"
#include "TK_Parse.h"
#include "Geo_Interface.h"
#include "Geo_Interpoint.h"
#include "Geo_FaceTopo_Unstruct.h"
#include "Geo_NodeTopo_Struct.h"
#include "Geo_FaceMetrics_Struct.h"
#include "Geo_CellMetrics_Struct.h"
#include "Geo_DynamicGridMetrics_Struct.h"
#include "Geo_OversetGridTopo_Struct.h"
#include "Geo_FaceMetrics_Unstruct.h"
#include "Geo_CellMetrics_Unstruct.h"
#include "Geo_DynamicGridMetrics_Unstruct.h"
#include "Geo_LSQWeight_Unstruct.h"
#ifdef PH_PARALLEL
#endif

#include "PHMpi.h"

using namespace std;
using namespace PHSPACE;
using namespace PHMPI;

namespace PHSPACE
{

JetData::JetData()
{
    velocity.resize(0);

    pressure    = 0.0;
    temperature = 0.0;
    density     = 0.0;

    specificHeatRatio = 1.4;
    molecularWeight   = 0.0;
    machNumber        = 0.0;

    SetNumberOfChemicalComponents(0);
    chemicalConstituents.resize(0);
}

JetData::~JetData()
{
    velocity.clear();
    chemicalConstituents.clear();

}

void JetData::SetJetName(const string &nozzleName)
{
    this->nozzleName = nozzleName;
}

void JetData::SetJetOnJudgement(const bool isJetOn)
{
    this->isJetOn = isJetOn;
}

void JetData::SetNumberOfChemicalComponents(const int nChemicalComponent)
{
    this->nChemicalComponent = nChemicalComponent;
}

void JetData::SetBCID(const int BCID)
{
    this->BCID = BCID;
}

void JetData::SetVelocity(const vector <RDouble> &velocity)
{
    this->velocity = velocity;
}

void JetData::SetPressure(const RDouble pressure)
{
    this->pressure = pressure;
}

void JetData::SetDensity(const RDouble density)
{
    this->density = density;
}

void JetData::SetTemperature(const RDouble temperature)
{
    this->temperature = temperature;
}

void JetData::SetSpecificHeatRatio(const RDouble specificHeatRatio)
{
    this->specificHeatRatio = specificHeatRatio;
}

void JetData::SetMolecularWeight(const RDouble molecularWeight)
{
    this->molecularWeight = molecularWeight;
}

void JetData::SetmachNumber(const RDouble machNumber)
{
    this->machNumber = machNumber;
}

void JetData::SetChemicalConstituents(const vector <RDouble> &chemicalConstituents)
{
    this->chemicalConstituents = chemicalConstituents;
}

string JetData::GetJetName()
{
    return this->nozzleName;
}

int JetData::GetBCID()
{
    return this->BCID;
}

bool JetData::GetJetOnJudgement()
{
    return this->isJetOn;
}

vector <RDouble> JetData::GetVelocity()
{
    return this->velocity;
}

RDouble JetData::GetPressure()
{
    return this->pressure;
}

RDouble JetData::GetDensity()
{
    return this->density;
}

RDouble JetData::GetTemperature()
{
    return this->temperature;
}

RDouble JetData::GetSpecificHeatRatio()
{
    return this->specificHeatRatio;
}

RDouble JetData::GetMolecularWeight()
{
    return this->molecularWeight;
}

RDouble JetData::GetMachNumber()
{
    return this->machNumber;
}

int JetData::getNumberOfChemicalComponents()
{
    return this->nChemicalComponent;
}

vector <RDouble> JetData::GetChemicalConstituents()
{
    return this->chemicalConstituents;
}

JetBCInfo::JetBCInfo()
{
    BCIDMax = -100000;
    BCIDMin = +100000;
    nJets   = 0;
    hasSimilarBC = 0;
    jetBCs.resize(0);
    jetFileExist = false;
}

JetBCInfo::~JetBCInfo()
{
    for (vector<JetData *>::iterator iter = jetBCs.begin(); iter != jetBCs.end(); ++ iter)
    {
        if (NULL != *iter)
        {
            delete *iter;
            *iter = NULL;
        }
    }
    jetBCs.clear();
}

void JetBCInfo::SetBCIDRange()
{
    for (int iJet = 0; iJet < nJets; ++ iJet)
    {
        BCIDMax = MAX(BCIDMax, jetBCs[iJet]->GetBCID());
        BCIDMin = MIN(BCIDMin, jetBCs[iJet]->GetBCID());
    }
}

void JetBCInfo::SetSimilarBCExistLabel()
{
    this->hasSimilarBC = 1;
}

int JetBCInfo::GetBCIDRangLower()
{
    return this->BCIDMin;
}

int JetBCInfo::GetBCIDRangUpper()
{
    return this->BCIDMax;
}

int JetBCInfo::GetNumberOfJets()
{
    return this->nJets;
}

void JetBCInfo::SetBCDataFileName(const string &fileName)
{
    this->jetFileName = fileName;
}

string JetBCInfo::GetBCDataFileName()
{
    return this->jetFileName;
}

vector <JetData*> JetBCInfo::GetBCInfo()
{
    return this->jetBCs;
}

void JetBCInfo::CheckSpecificBCExist()
{
    int hasSpecificBCGlobalProcessor = 0;
    int hasSpecificBCLocalProcessor  = 0;

    int nZones = GetNumberofGlobalZones();

    for (int izone = 0; izone < nZones; ++ izone)
    {
        Grid *gridIn = GetGrid(izone, 0);

        if (!gridIn) continue;

        int systemGridType = STRUCTGRID;
        systemGridType = GlobalDataBase::GetIntParaFromDB("sys_gridtype");

        if (systemGridType == PHSPACE::STRUCTGRID)
        {
            StructGrid *grid = StructGridCast(gridIn);

            CompositeBCRegion *compositeBCRegion = grid->GetCompositeBCRegion();
            int nBCRegion = compositeBCRegion->GetnBCRegion();

            for (int iBCRegion = 0; iBCRegion < nBCRegion; ++ iBCRegion)
            {
                BCRegion *BCRegion = compositeBCRegion->GetBCRegion(iBCRegion);
                int BCType = BCRegion->GetBCType();
                for (int iJet = 0; iJet < GetNumberOfJets(); ++ iJet)
                {
                    if (BCType == jetBCs[iJet]->GetBCID() && !jetBCs[iJet]->GetJetOnJudgement())
                    {
                        BCType = PHENGLEI::SOLID_SURFACE;
                        BCRegion->SetBCType(BCType);
                    }
                }
                if (BCType >= BCIDMin && BCType <= BCIDMax)
                {
                    hasSpecificBCLocalProcessor = 1;
                    break;
                }
            }
        }
        else if (systemGridType == PHSPACE::UNSTRUCTGRID)
        {
            UnstructGrid *grid = UnstructGridCast(gridIn);
            BCRecord **BCRecord = grid->GetBCRecord();
            for (int iFace = 0; iFace < grid->GetNBoundFace(); ++ iFace)
            {
                int BCType = BCRecord[iFace]->GetKey();
                for (int iJet = 0; iJet < GetNumberOfJets(); ++ iJet)
                {
                    if (BCType == jetBCs[iJet]->GetBCID() && !jetBCs[iJet]->GetJetOnJudgement())
                    {
                        BCType = PHENGLEI::SOLID_SURFACE;
                        BCRecord[iFace]->SetKey(BCType);
                    }
                }
                if (BCType >= BCIDMin && BCType <= BCIDMax)
                {
                    hasSpecificBCLocalProcessor = 1;
                    break;
                }
            }
        }
    }

    PH_Reduce(&hasSpecificBCLocalProcessor, &hasSpecificBCGlobalProcessor, 1, PH_MAX);

    if (hasSpecificBCGlobalProcessor != 0)
    {
        SetSimilarBCExistLabel();
    }
}

bool JetBCInfo::IfDefineFileExist()
{
    return jetFileExist;
}

int JetBCInfo::GetJetBCExistLabel()
{
    return this->hasSimilarBC;
}

void JetBCInfo::Nondimensional()
{
//! The jet conditions are passed in through a common in non-dimensional form as:
//!      UJET   = u(jet) / Vref
//!      VJET   = v(jet) / Vref
//!      WJET   = w(jet) / Vref
//!      RJET   = r(jet) / rhoref
//!      TJET   = T(jet) / Tref
//!      PJET   = p(jet) / (rhoref * Vref**2)

    for (int iJet = 0; iJet < GetNumberOfJets(); ++ iJet)
    {
        RDouble density           = jetBCs[iJet]->GetDensity();
        RDouble pressure          = jetBCs[iJet]->GetPressure();
        RDouble temperature       = jetBCs[iJet]->GetTemperature();
        vector <RDouble> velocity = jetBCs[iJet]->GetVelocity();

        RDouble refDimensionalTemperature = GlobalDataBase::GetDoubleParaFromDB("refDimensionalTemperature");
        RDouble refDimensionalSonicSpeed  = GlobalDataBase::GetDoubleParaFromDB("refDimensionalSonicSpeed");
        RDouble refDimensionalDensity     = GlobalDataBase::GetDoubleParaFromDB("refDimensionalDensity");
        RDouble refMachNumber             = GlobalDataBase::GetDoubleParaFromDB("refMachNumber");

        RDouble refVelocityDimensional = refMachNumber * refDimensionalSonicSpeed;

        density /= refDimensionalDensity;
        pressure /= (refDimensionalDensity * refVelocityDimensional * refVelocityDimensional);
        temperature /= refDimensionalTemperature;

        for (int iDim = 0; iDim < velocity.size(); ++ iDim)
        {
            velocity[iDim] /= refVelocityDimensional;
        }

        jetBCs[iJet]->SetDensity(density);
        jetBCs[iJet]->SetPressure(pressure);
        jetBCs[iJet]->SetTemperature(temperature);
        jetBCs[iJet]->SetVelocity(velocity);
    }
}

void JetBCInfo::CheckIfBCDataFileExist()
{
    using namespace PHMPI;
    jetFileExist = false;

    fstream file;
    string keyFile = GetBCDataFileName();
    file.open(keyFile.c_str(), ios_base::in);

    if (!file)
    {
        return;
    }

    jetFileExist = true;
}

void JetBCInfo::ReadBCData()
{
    fstream file;
    string keyFile = GetBCDataFileName();
    file.open(keyFile.c_str(), ios_base::in);

    string line, word;
    string separator  = " =\t\r\n#$,;\"'";
    string::size_type npos = string::npos;

    string nozzleName;
    int BCID;
    bool isJetOn;
    RDouble pressure, temperature, density;
    RDouble specificHeatRatio, molecularWeight, machNumber;

    getline(file, line);
    line = FindNextWord(line, word, separator);
    from_string <int> (nJets, word, std::dec);

    for (int iJet = 0; iJet < nJets; ++ iJet)
    {
        JetData *jet = new JetData;
        SkipLines(file, 1);

        getline(file, line);
        line = FindNextWord(line, word, separator);
        from_string <string> (nozzleName, word, std::dec);

        line = FindNextWord(line, word, separator);
        from_string <int> (BCID, word, std::dec);

        line = FindNextWord(line, word, separator);
        from_string <bool> (isJetOn, word, std::dec);

        const int nDim = 3;
        vector <RDouble> velocity;
        getline(file, line);
        for (int iDim = 0; iDim < nDim; ++ iDim)
        {
            line = FindNextWord(line, word, separator);
            RDouble velocityIn;
            from_string <RDouble> (velocityIn, word, std::dec);
            velocity.push_back(velocityIn);
        }

        getline(file, line);

        line = FindNextWord(line, word, separator);
        from_string <RDouble> (pressure, word, std::dec);

        line = FindNextWord(line, word, separator);
        from_string <RDouble> (temperature, word, std::dec);

        line = FindNextWord(line, word, separator);
        from_string <RDouble> (density, word, std::dec);

        getline(file, line);

        line = FindNextWord(line, word, separator);
        from_string <RDouble> (specificHeatRatio, word, std::dec);

        line = FindNextWord(line, word, separator);
        from_string <RDouble> (molecularWeight, word, std::dec);

        line = FindNextWord(line, word, separator);
        from_string <RDouble> (machNumber, word, std::dec);

        getline(file, line);
        int nChemicalComponent = 5;
        vector <RDouble> chemicalConstituents;
        for (int iChemComponent = 0; iChemComponent < nChemicalComponent; ++ iChemComponent)
        {
            line = FindNextWord(line, word, separator);
            RDouble chemicalConstituentsIn;
            from_string <RDouble> (chemicalConstituentsIn, word, std::dec);
            chemicalConstituents.push_back(chemicalConstituentsIn);
        }

        jet->SetBCID(BCID);
        jet->SetJetName(nozzleName);
        jet->SetJetOnJudgement(isJetOn);
        jet->SetVelocity(velocity);
        jet->SetPressure(pressure);
        jet->SetTemperature(temperature);
        jet->SetDensity(density);
        jet->SetSpecificHeatRatio(specificHeatRatio);
        jet->SetMolecularWeight(molecularWeight);
        jet->SetmachNumber(machNumber);
        jet->SetNumberOfChemicalComponents(nChemicalComponent);
        jet->SetChemicalConstituents(chemicalConstituents);

        jetBCs.push_back(jet);
    }

    file.close();
    file.clear();
}

void JetBCInfo::BcastBCData()
{
#ifndef PH_PARALLEL
    return;
#endif

    int serverProcessorID = GetServerProcessorID();
    int myProcessorID   = GetCurrentProcessorID();
    
    DataContainer *Jetdata = new DataContainer();
    
    if (myProcessorID == serverProcessorID)
    {
        CompressJetData(Jetdata);
    }
    
    PH_Bcast(Jetdata, server);
    
    if (myProcessorID != serverProcessorID)
    {
        DecompressJetData(Jetdata);
    }

    delete Jetdata;
}

void JetBCInfo::CompressJetData(DataContainer *jetcData)
{
    jetcData->moveto_begin();

    int nJets = GetNumberOfJets();
    jetcData->append(&nJets, sizeof(int));

    for (int iJet = 0; iJet < nJets; ++ iJet)
    {
        string nozzleName = jetBCs[iJet]->GetJetName();
        int BCID = jetBCs[iJet]->GetBCID();
        bool isJetOn = jetBCs[iJet]->GetJetOnJudgement();
        vector <RDouble> velocity = jetBCs[iJet]->GetVelocity();
        RDouble pressure          = jetBCs[iJet]->GetPressure();
        RDouble temperature       = jetBCs[iJet]->GetTemperature();
        RDouble density           = jetBCs[iJet]->GetDensity();
        RDouble specificHeatRatio = jetBCs[iJet]->GetSpecificHeatRatio();
        RDouble molecularWeight   = jetBCs[iJet]->GetMolecularWeight();
        RDouble machNumber        = jetBCs[iJet]->GetMachNumber();
        vector <RDouble> chemicalConstituents = jetBCs[iJet]->GetChemicalConstituents();

        int nDim =static_cast<int>( velocity.size());
        int nChemicalComponent = static_cast<int>(chemicalConstituents.size());

        jetcData->append_string(nozzleName);
        jetcData->append(&BCID, sizeof(int));
        jetcData->append(&isJetOn, sizeof(bool));
        jetcData->append(&nDim, sizeof(int));
        for (int iDim = 0; iDim < nDim; ++ iDim)
        {
            jetcData->append(&velocity[iDim], sizeof(RDouble));
        }
        jetcData->append(&pressure, sizeof(RDouble));
        jetcData->append(&temperature, sizeof(RDouble));
        jetcData->append(&density, sizeof(RDouble));
        jetcData->append(&specificHeatRatio, sizeof(RDouble));
        jetcData->append(&molecularWeight, sizeof(RDouble));
        jetcData->append(&machNumber, sizeof(RDouble));
        jetcData->append(&nChemicalComponent, sizeof(int));
        for (int iChemicalComponent = 0; iChemicalComponent < nChemicalComponent; ++ iChemicalComponent)
        {
            jetcData->append(&chemicalConstituents[iChemicalComponent], sizeof(RDouble));
        }
    }
}

void JetBCInfo::DecompressJetData(DataContainer *jetcData)
{
    jetcData->moveto_begin();

    nJets = GetNumberOfJets();
    jetcData->read(&nJets, sizeof(int));

    jetBCs.resize(0);

    for (int iJet = 0; iJet < nJets; ++ iJet)
    {
        JetData *jetOnProcessor = new JetData;
        int BCID;
        string nozzleName;
        bool isJetOn;
        RDouble pressure, temperature, density;
        RDouble specificHeatRatio, molecularWeight, machNumber;

        int nDim, nChemicalComponent;

        jetcData->read_string(nozzleName);
        jetcData->read(&BCID, sizeof(int));
        jetcData->read(&isJetOn, sizeof(bool));

        vector <RDouble> velocity;
        jetcData->read(&nDim, sizeof(int));
        velocity.resize(nDim);
        for (int iDim = 0; iDim < nDim; ++ iDim)
        {
            jetcData->read(&velocity[iDim], sizeof(RDouble));
        }

        jetcData->read(&pressure, sizeof(RDouble));
        jetcData->read(&temperature, sizeof(RDouble));
        jetcData->read(&density, sizeof(RDouble));
        jetcData->read(&specificHeatRatio, sizeof(RDouble));
        jetcData->read(&molecularWeight, sizeof(RDouble));
        jetcData->read(&machNumber, sizeof(RDouble));

        vector <RDouble> chemicalConstituents;
        jetcData->read(&nChemicalComponent, sizeof(int));
        chemicalConstituents.resize(nChemicalComponent);
        for (int iChemicalComponent = 0; iChemicalComponent < nChemicalComponent; ++ iChemicalComponent)
        {
            jetcData->read(&chemicalConstituents[iChemicalComponent], sizeof(RDouble));
        }

        jetOnProcessor->SetBCID(BCID);
        jetOnProcessor->SetJetName(nozzleName);
        jetOnProcessor->SetJetOnJudgement(isJetOn);
        jetOnProcessor->SetVelocity(velocity);
        jetOnProcessor->SetPressure(pressure);
        jetOnProcessor->SetTemperature(temperature);
        jetOnProcessor->SetDensity(density);
        jetOnProcessor->SetSpecificHeatRatio(specificHeatRatio);
        jetOnProcessor->SetMolecularWeight(molecularWeight);
        jetOnProcessor->SetmachNumber(machNumber);
        jetOnProcessor->SetChemicalConstituents(chemicalConstituents);

        jetBCs.push_back(jetOnProcessor);
    }
}

JetData *JetBCInfo::GetInUseJets(const int BCType)
{
    JetData *jetNull = NULL;
    for (int iJet = 0; iJet < this->GetNumberOfJets(); ++ iJet)
    {
        JetData *jet = this->GetBCInfo()[iJet];

        if (BCType == jet->GetBCID() && jet->GetJetOnJudgement())
        {
            return jet;
        }
    }
    return jetNull;
}

JetBCInfo globalJetBC;

void GetJetBCData()
{
    string jetDefineFile = GlobalDataBase::GetStrParaFromDB("jetDefineFile");
    globalJetBC.SetBCDataFileName(jetDefineFile);
    globalJetBC.CheckIfBCDataFileExist();

    using namespace PHMPI;

    if (GetCurrentProcessorID() == server && globalJetBC.IfDefineFileExist())
    {
        globalJetBC.ReadBCData();
    }

    if (!globalJetBC.IfDefineFileExist())
    {
        return;
    }

    globalJetBC.BcastBCData();

    globalJetBC.SetBCIDRange();

    globalJetBC.CheckSpecificBCExist();

    globalJetBC.Nondimensional();
}

bool ifIsJetBC(const int BCType)
{
    if (globalJetBC.GetBCIDRangLower() <= BCType && globalJetBC.GetBCIDRangUpper() >= BCType)
    {
        for (int iJet = 0; iJet < globalJetBC.GetNumberOfJets(); ++ iJet)
        {
            JetData *jet = globalJetBC.GetBCInfo()[iJet];

            if (BCType == jet->GetBCID() && jet->GetJetOnJudgement())
            {
                return true;
            }
        }
    }
    return false;
}

void JetBC3D(Grid *grid, BCRegion *BCRegion)
{
    if (!globalJetBC.IfDefineFileExist())
    {
        return;
    }

    StructGrid *gridStructure = StructGridCast(grid);
    int ni = gridStructure->GetNI();
    int nj = gridStructure->GetNJ();
    int nk = gridStructure->GetNK();

    RDouble4D &q = *reinterpret_cast <RDouble4D *> (gridStructure->GetDataPtr("q"));

    int iStart, iEnd, jStart, jEnd, kStart, kEnd;
    BCRegion->GetIJKRegion(iStart, iEnd, jStart, jEnd, kStart, kEnd);

    int *faceDirection3D = BCRegion->GetFaceDirectionIndex();

    JetData *jetInUse = globalJetBC.GetInUseJets(BCRegion->GetBCType());

    vector <RDouble> velocity = jetInUse->GetVelocity();
    RDouble pressure          = jetInUse->GetPressure();
    RDouble density           = jetInUse->GetDensity();
    //vector <RDouble> chemicalConstituents = jet.GetChemicalConstituents();

    for (int k = kStart; k <= kEnd; ++ k)
    {
        for (int j = jStart; j <= jEnd; ++ j)
        {
            for (int i = iStart; i <= iEnd; ++ i)
            {
                for (int iGhostlayer = 1; iGhostlayer <= 2; ++ iGhostlayer)
                {
                    int iIndexGhost = i + faceDirection3D[0] * iGhostlayer;
                    int jIndexGhost = j + faceDirection3D[1] * iGhostlayer;
                    int kIndexGhost = k + faceDirection3D[2] * iGhostlayer;

                    q(iIndexGhost, jIndexGhost, kIndexGhost, 0) = density;
                    q(iIndexGhost, jIndexGhost, kIndexGhost, 1) = velocity[0];
                    q(iIndexGhost, jIndexGhost, kIndexGhost, 2) = velocity[1];
                    q(iIndexGhost, jIndexGhost, kIndexGhost, 3) = velocity[2];
                    q(iIndexGhost, jIndexGhost, kIndexGhost, 4) = pressure;
                }
            }
        }
    }
}

void JetBC3D(RDouble *primTarget, int BCType)
{
    if (!globalJetBC.IfDefineFileExist())
    {
        return;
    }

    JetData *jetInUse = globalJetBC.GetInUseJets(BCType);

    vector <RDouble> velocity = jetInUse->GetVelocity();
    RDouble pressure          = jetInUse->GetPressure();
    RDouble density           = jetInUse->GetDensity();
    //RDouble *chemicalConstituents = jet.GetChemicalConstituents();

    primTarget[0] = density;
    primTarget[1] = velocity[0];
    primTarget[2] = velocity[1];
    primTarget[3] = velocity[2];
    primTarget[4] = pressure;
}

void JetBC3DRegion(Grid *gridIn, BCRegionUnstruct *bcRegionUnstruct)
{
    if (!globalJetBC.IfDefineFileExist())
    {
        return;
    }

    // get the basic variables for boundary calculating.
    UnstructGrid *grid   = UnstructGridCast(gridIn);
    int *rightCellofFace = grid->GetRightCellOfFace();
    RDouble   *xfn       = grid->GetFaceNormalX();
    RDouble   *yfn       = grid->GetFaceNormalY();
    RDouble   *zfn       = grid->GetFaceNormalZ();
    RDouble   *vgn       = grid->GetFaceNormalVelocity();
    RDouble **q = reinterpret_cast< RDouble ** > (grid->GetDataPtr("q"));

    int BCType = bcRegionUnstruct->GetBCType();
    JetData *jetInUse = globalJetBC.GetInUseJets(BCType);

    vector <RDouble> velocity = jetInUse->GetVelocity();
    RDouble pressure          = jetInUse->GetPressure();
    RDouble density           = jetInUse->GetDensity();

    int iFace, re;	
    //! get the face indexes of this boundary condition region.
    vector<int> *faceIndex = bcRegionUnstruct->GetFaceIndex();
    using namespace IDX;
    for (vector<int>::iterator iter = faceIndex->begin(); iter != faceIndex->end(); iter++)
    {
        iFace = *iter;
        re = rightCellofFace[iFace];		
        q[IR][re] = density;
        q[IU][re] = velocity[0];
        q[IV][re] = velocity[1];
        q[IW][re] = velocity[2];
        q[IP][re] = pressure;
    }   
}

}
