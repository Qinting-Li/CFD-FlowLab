#include "Force.h"
#include "Math_BasisFunction.h"
#include "Constants.h"
#include "PHIO.h"
#include "IO_FileName.h"
#include "Geo_SimpleBC.h"
#include "TK_Time.h"
#include "TK_Parse.h"
#include "TK_Exit.h"
#include "Glb_Dimension.h"
#include "GlobalDataBase.h"
#include "PHHeader.h"
#include "PHMpi.h"
#include "AleManager.h"

using namespace std;

namespace PHSPACE
{

Force globalAerodynamicForce;

void InitForce (ActionKey *actkey)
{

    RDouble attack,sideslip;
    RDouble attackd = GlobalDataBase::GetDoubleParaFromDB("attackd");
    RDouble angleSlide = GlobalDataBase::GetDoubleParaFromDB("angleSlide");

    attack   = attackd * PI / 180.0;
    sideslip = angleSlide * PI / 180.0;

    RDouble forceReferenceArea   = GlobalDataBase::GetDoubleParaFromDB("forceReferenceArea");
    RDouble forceReferenceLength = GlobalDataBase::GetDoubleParaFromDB("forceReferenceLength");
    RDouble forceReferenceLengthSpanWise = GlobalDataBase::GetDoubleParaFromDB("forceReferenceLengthSpanWise");

    globalAerodynamicForce.Init(attack, sideslip, forceReferenceArea, forceReferenceLength, forceReferenceLengthSpanWise);
}


void CollectionForce(ActionKey *actkey)
{
    globalAerodynamicForce.CollectionForce(actkey->GetData());
}

void DumpForce(ActionKey *actkey)
{
    int iterationStep;
    GlobalDataBase::GetData("outnstep", &iterationStep, PHINT, 1);

    globalAerodynamicForce.DumpForce(actkey, iterationStep);
}

Force::Force()
{
    partBCID = 0;

    cpx = 0;
    cpy = 0;
    cpz = 0;

    cl_tot   = 0;
    cd_tot   = 0;
    cd_pr    = 0;
    cd_sf    = 0;
    cd_cl2pa = 0;
    xcp      = 0;
    ycp      = 0;
    zcp      = 0;
    side     = 0;

    CA_f = 0;
    CA_p = 0;
    CN_f = 0;
    CN_p = 0;
    CZ_f = 0;
    CZ_p = 0;
    Cl_f = 0;
    Cl_p = 0;
    Cn_f = 0;
    Cn_p = 0;
    Cm_f = 0;
    Cm_p = 0;

    CA_f_body  = 0;
    CN_f_body  = 0;
    CZ_f_body = 0;
    Cl_f_body = 0;
    Cn_f_body = 0;
    Cm_f_body  = 0;
    CA_p_body  = 0;
    CN_p_body  = 0;
    CZ_p_body = 0;
    Cl_p_body = 0;
    Cn_p_body = 0;
    Cm_p_body  = 0;

    numberOfBody = 0;
    wallNameOfEachBody.resize(0);
}

Force::~Force()
{
    delete [] partBCID;    partBCID = nullptr;

    delete [] cpx;    cpx = nullptr;
    delete [] cpy;    cpy = nullptr;
    delete [] cpz;    cpz = nullptr;

    delete [] cl_tot;    cl_tot = nullptr;
    delete [] cd_tot;    cd_tot = nullptr;
    delete [] cd_pr;    cd_pr = nullptr;
    delete [] cd_sf;    cd_sf = nullptr;
    delete [] cd_cl2pa;    cd_cl2pa = nullptr;
    delete [] xcp;    xcp = nullptr;
    delete [] ycp;    ycp = nullptr;
    delete [] zcp;    zcp = nullptr;
    delete [] side;    side = nullptr;

    delete [] CA_f;    CA_f = nullptr;
    delete [] CA_p;    CA_p = nullptr;
    delete [] CN_f;    CN_f = nullptr;
    delete [] CN_p;    CN_p = nullptr;
    delete [] CZ_f;    CZ_f = nullptr;
    delete [] CZ_p;    CZ_p = nullptr;
    delete [] Cl_f;    Cl_f = nullptr;
    delete [] Cl_p;    Cl_p = nullptr;
    delete [] Cn_f;    Cn_f = nullptr;
    delete [] Cn_p;    Cn_p = nullptr;
    delete [] Cm_f;    Cm_f = nullptr;
    delete [] Cm_p;    Cm_p = nullptr;

    delete [] CA_f_body;    CA_f_body  = nullptr;
    delete [] CN_f_body;    CN_f_body  = nullptr;
    delete [] CZ_f_body;    CZ_f_body  = nullptr;
    delete [] Cl_f_body;    Cl_f_body  = nullptr;
    delete [] Cn_f_body;    Cn_f_body  = nullptr;
    delete [] Cm_f_body;    Cm_f_body  = nullptr;
    delete [] CA_p_body;    CA_p_body  = nullptr;
    delete [] CN_p_body;    CN_p_body  = nullptr;
    delete [] CZ_p_body;    CZ_p_body  = nullptr;
    delete [] Cl_p_body;    Cl_p_body  = nullptr;
    delete [] Cn_p_body;    Cn_p_body  = nullptr;
    delete [] Cm_p_body;    Cm_p_body  = nullptr;
}

void Force::Init(RDouble attack, RDouble sideslip, RDouble forceReferenceArea, RDouble forceReferenceLength, RDouble forceReferenceLengthSpanWise)
{
    this->attack   = attack;
    this->sideslip = sideslip;
    this->forceReferenceArea   = forceReferenceArea;
    this->forceReferenceLength = forceReferenceLength;
    this->forceReferenceLengthSpanWise = forceReferenceLengthSpanWise;

    isPartExist = (GlobalBoundaryCondition::GetGlobalBoundaryConditionList() != 0);
    numberofSolidWallPart = 0;
    if(isPartExist)
    {
        numberofSolidWallPart = GlobalBoundaryCondition::GetNumberofSolidWallPart();
    }
    //! Because of global and local part force, numberofSolidWallPart * 2
    numberofTotalPart = numberofSolidWallPart * 2 + 1;

    int myid   = PHMPI::GetCurrentProcessorID();
    int server = PHMPI::GetServerProcessorID();
    if((myid == server) && (!cpx))
    {
        cpx = new RDouble [numberofTotalPart];
        cpy = new RDouble [numberofTotalPart];
        cpz = new RDouble [numberofTotalPart];
        hingeMoment = new RDouble[numberofTotalPart];

        CA_f = new RDouble [numberofTotalPart];
        CA_p = new RDouble [numberofTotalPart];
        CN_f = new RDouble [numberofTotalPart];
        CN_p = new RDouble [numberofTotalPart];
        CZ_f = new RDouble [numberofTotalPart];
        CZ_p = new RDouble [numberofTotalPart];
        Cl_f = new RDouble [numberofTotalPart];
        Cl_p = new RDouble [numberofTotalPart];
        Cn_f = new RDouble [numberofTotalPart];
        Cn_p = new RDouble [numberofTotalPart];
        Cm_f = new RDouble [numberofTotalPart];
        Cm_p = new RDouble [numberofTotalPart];

        PHSPACE::SetField(cpx, 0.0, numberofTotalPart);
        PHSPACE::SetField(cpy, 0.0, numberofTotalPart);
        PHSPACE::SetField(cpz, 0.0, numberofTotalPart);
        PHSPACE::SetField(hingeMoment, 0.0, numberofTotalPart);

        PHSPACE::SetField(CA_f, 0.0, numberofTotalPart);
        PHSPACE::SetField(CA_p, 0.0, numberofTotalPart);
        PHSPACE::SetField(CN_f, 0.0, numberofTotalPart);
        PHSPACE::SetField(CN_p, 0.0, numberofTotalPart);
        PHSPACE::SetField(CZ_f, 0.0, numberofTotalPart);
        PHSPACE::SetField(CZ_p, 0.0, numberofTotalPart);
        PHSPACE::SetField(Cl_f, 0.0, numberofTotalPart);
        PHSPACE::SetField(Cl_p, 0.0, numberofTotalPart);
        PHSPACE::SetField(Cn_f, 0.0, numberofTotalPart);
        PHSPACE::SetField(Cn_p, 0.0, numberofTotalPart);
        PHSPACE::SetField(Cm_f, 0.0, numberofTotalPart);
        PHSPACE::SetField(Cm_p, 0.0, numberofTotalPart);
    }

    numberOfBody = GlobalBoundaryCondition::GetNumberOfBody();
    if ((myid == server) && (!CA_f_body) && (numberOfBody))
    {
        CA_f_body = new RDouble[numberOfBody];
        CN_f_body = new RDouble[numberOfBody];
        CZ_f_body = new RDouble[numberOfBody];
        Cl_f_body = new RDouble[numberOfBody];
        Cn_f_body = new RDouble[numberOfBody];
        Cm_f_body = new RDouble[numberOfBody];
        CA_p_body = new RDouble[numberOfBody];
        CN_p_body = new RDouble[numberOfBody];
        CZ_p_body = new RDouble[numberOfBody];
        Cl_p_body = new RDouble[numberOfBody];
        Cn_p_body = new RDouble[numberOfBody];
        Cm_p_body = new RDouble[numberOfBody];

        set < string > *bodyList = GlobalBoundaryCondition::GetBodyNameList();
        set < string >::iterator bodyListIter;
        wallNameOfEachBody.resize(numberOfBody);
        int bodyIndex = 0;
        for (bodyListIter = bodyList->begin(); bodyListIter != bodyList->end(); ++ bodyListIter)
        {
            string bodyNameInList = *bodyListIter;
            wallNameOfEachBody[bodyIndex].push_back(bodyNameInList);
            bodyIndex ++;
        }
    }

    if ((myid == server) && (numberOfBody))
    {
        PHSPACE::SetField(CA_f_body, 0.0, numberOfBody);
        PHSPACE::SetField(CN_f_body, 0.0, numberOfBody);
        PHSPACE::SetField(CZ_f_body, 0.0, numberOfBody);
        PHSPACE::SetField(Cl_f_body, 0.0, numberOfBody);
        PHSPACE::SetField(Cn_f_body, 0.0, numberOfBody);
        PHSPACE::SetField(Cm_f_body, 0.0, numberOfBody);
        PHSPACE::SetField(CA_p_body, 0.0, numberOfBody);
        PHSPACE::SetField(CN_p_body, 0.0, numberOfBody);
        PHSPACE::SetField(CZ_p_body, 0.0, numberOfBody);
        PHSPACE::SetField(Cl_p_body, 0.0, numberOfBody);
        PHSPACE::SetField(Cn_p_body, 0.0, numberOfBody);
        PHSPACE::SetField(Cm_p_body, 0.0, numberOfBody);
    }

    if ((!partBCID) && (myid == server))
    {
        partBCID = new int [numberofTotalPart];

        const int ALL_SOLID_WALL = -1;
        PHSPACE::SetField(partBCID, ALL_SOLID_WALL, numberofTotalPart);

        int partCount = 1;
        uint_t nTBC = GlobalBoundaryCondition::GetNumberofBC();
        for (int iBC = 0; iBC < nTBC; ++ iBC)
        {
            if((!GlobalBoundaryCondition::IsSolidWallBC(iBC)) && ((GlobalBoundaryCondition::GetBC(iBC)->GetBCType() != PHENGLEI::ABLATION_SURFACE))) continue;

            //! Because of global and local part force, partBCID set twice
            partBCID[partCount ++] = iBC;
            partBCID[partCount ++] = iBC;

            if (numberOfBody == 0)
            {
                continue;
            }
            string wallName = GlobalBoundaryCondition::GetBoundaryName(iBC);
            int bodyIndex = GlobalBoundaryCondition::GetBodyIndex(iBC);
            wallNameOfEachBody[bodyIndex].push_back(wallName);
        }
    }
}

void Force::CollectionForce(DataContainer *data_in)
{
    data_in->moveto_begin();
    int nLocalZones = PHMPI::GetNumberofLocalZones();

    int nTotalVarDim = nForceVar * numberofTotalPart;
    RDouble *localForceVar = new RDouble [nTotalVarDim];
    PHSPACE::SetField(localForceVar, 0.0, nTotalVarDim);

    int nTBC = static_cast<int>( GlobalBoundaryCondition::GetNumberofBC());
   for (int iZone = 0; iZone < nLocalZones; ++ iZone)
    {
        int varCount = 0;

        //! '-1'     , represent total wall force;\n
        //! 0 - nWALL, represent each wall part force.
        for (int iPart = ALL_SOLID_WALL; iPart < nTBC; ++ iPart)
        {
            bool ifRead = (iPart == ALL_SOLID_WALL); // total wall force.
            if(iPart >= 0)
            {
                ifRead = ((GlobalBoundaryCondition::IsSolidWallBC(iPart)) || ((GlobalBoundaryCondition::GetBC(iPart)->GetBCType() == PHENGLEI::ABLATION_SURFACE)));
            }
            if(!ifRead) continue;
            
            RDouble dCA_f, dCA_p, dCN_f, dCN_p, dCZ_f, dCZ_p, dcpx, dcpy, dcpz, dCl_f, dCl_p, dCn_f, dCn_p, dCm_f, dCm_p, hingeMoment;

            data_in->read(&dCA_f, sizeof(RDouble));
            data_in->read(&dCA_p, sizeof(RDouble));
            data_in->read(&dCN_f, sizeof(RDouble));
            data_in->read(&dCN_p, sizeof(RDouble));
            data_in->read(&dCZ_f, sizeof(RDouble));
            data_in->read(&dCZ_p, sizeof(RDouble));

            data_in->read(&dcpx,sizeof(RDouble));
            data_in->read(&dcpy,sizeof(RDouble));
            data_in->read(&dcpz,sizeof(RDouble));

            data_in->read(&dCl_f, sizeof(RDouble));
            data_in->read(&dCl_p, sizeof(RDouble));
            data_in->read(&dCn_f, sizeof(RDouble));
            data_in->read(&dCn_p, sizeof(RDouble));
            data_in->read(&dCm_f, sizeof(RDouble));
            data_in->read(&dCm_p, sizeof(RDouble));

            data_in->read(&hingeMoment, sizeof(RDouble));

            localForceVar[varCount ++] += dCA_f;
            localForceVar[varCount ++] += dCA_p;
            localForceVar[varCount ++] += dCN_f;
            localForceVar[varCount ++] += dCN_p;
            localForceVar[varCount ++] += dCZ_f;
            localForceVar[varCount ++] += dCZ_p;

            localForceVar[varCount ++] += dcpx;
            localForceVar[varCount ++] += dcpy;
            localForceVar[varCount ++] += dcpz;

            localForceVar[varCount ++] += dCl_f;
            localForceVar[varCount ++] += dCl_p;
            localForceVar[varCount ++] += dCn_f;
            localForceVar[varCount ++] += dCn_p;
            localForceVar[varCount ++] += dCm_f;
            localForceVar[varCount ++] += dCm_p;

            localForceVar[varCount ++] += hingeMoment;

            //! Because of global and local part force, when (iPart > ALL_SOLID_WALL), read second.
            if (iPart > ALL_SOLID_WALL)
            {
                data_in->read(&dCA_f, sizeof(RDouble));
                data_in->read(&dCA_p, sizeof(RDouble));
                data_in->read(&dCN_f, sizeof(RDouble));
                data_in->read(&dCN_p, sizeof(RDouble));
                data_in->read(&dCZ_f, sizeof(RDouble));
                data_in->read(&dCZ_p, sizeof(RDouble));

                data_in->read(&dcpx,sizeof(RDouble));
                data_in->read(&dcpy,sizeof(RDouble));
                data_in->read(&dcpz,sizeof(RDouble));

                data_in->read(&dCl_f, sizeof(RDouble));
                data_in->read(&dCl_p, sizeof(RDouble));
                data_in->read(&dCn_f, sizeof(RDouble));
                data_in->read(&dCn_p, sizeof(RDouble));
                data_in->read(&dCm_f, sizeof(RDouble));
                data_in->read(&dCm_p, sizeof(RDouble));

                data_in->read(&hingeMoment, sizeof(RDouble));

                localForceVar[varCount ++] += dCA_f;
                localForceVar[varCount ++] += dCA_p;
                localForceVar[varCount ++] += dCN_f;
                localForceVar[varCount ++] += dCN_p;
                localForceVar[varCount ++] += dCZ_f;
                localForceVar[varCount ++] += dCZ_p;

                localForceVar[varCount ++] += dcpx;
                localForceVar[varCount ++] += dcpy;
                localForceVar[varCount ++] += dcpz;

                localForceVar[varCount ++] += dCl_f;
                localForceVar[varCount ++] += dCl_p;
                localForceVar[varCount ++] += dCn_f;
                localForceVar[varCount ++] += dCn_p;
                localForceVar[varCount ++] += dCm_f;
                localForceVar[varCount ++] += dCm_p;

                localForceVar[varCount ++] += hingeMoment;
            }
        }
    }

    //! Collection.
    int myid   = PHMPI::GetCurrentProcessorID();
    int server = PHMPI::GetServerProcessorID();

    RDouble * globalForceVar = new RDouble[nTotalVarDim];
    if(PHMPI::GetNumberOfProcessor() > 1)
    {
        PH_AllreduceSepMode(localForceVar, globalForceVar, nTotalVarDim, MPI_DOUBLE, PH_SUM);
    }
    else
    {
        for(int iVar = 0; iVar < nTotalVarDim; ++ iVar)
        {
            globalForceVar[iVar] = localForceVar[iVar];
        }
    }

    if(myid == server)
    {
        int varCount = 0;
        for (int iPart = 0; iPart < numberofTotalPart; ++ iPart)
        {
            CA_f[iPart] = globalForceVar[varCount ++];
            CA_p[iPart] = globalForceVar[varCount ++];
            CN_f[iPart] = globalForceVar[varCount ++];
            CN_p[iPart] = globalForceVar[varCount ++];
            CZ_f[iPart] = globalForceVar[varCount ++];
            CZ_p[iPart] = globalForceVar[varCount ++];

            cpx[iPart] = globalForceVar[varCount ++];
            cpy[iPart] = globalForceVar[varCount ++];
            cpz[iPart] = globalForceVar[varCount ++];

            Cl_f[iPart] = globalForceVar[varCount ++];
            Cl_p[iPart] = globalForceVar[varCount ++];
            Cn_f[iPart] = globalForceVar[varCount ++];
            Cn_p[iPart] = globalForceVar[varCount ++];
            Cm_f[iPart] = globalForceVar[varCount ++];
            Cm_p[iPart] = globalForceVar[varCount ++];

            hingeMoment[iPart] = globalForceVar[varCount ++];
        }
    }

    delete []localForceVar;
    delete []globalForceVar;
}

void Force::CompForce()
{
    RDouble sina,cosa,sinb,cosb;

    sina = sin(attack);
    cosa = cos(attack);
    sinb = sin(sideslip);
    cosb = cos(sideslip);

    bool isGlobalCoordinate = false;
    for (int iPart = 0; iPart < numberofTotalPart; ++ iPart)
    {
        RDouble scf,vcm, Sb, forceReferenceAreaPart, forceReferenceLengthPart, forceReferenceLengthPart_B;
        int bcID = partBCID[iPart];
        if((bcID == ALL_SOLID_WALL) || (isGlobalCoordinate == true))
        {
            forceReferenceAreaPart   = this->forceReferenceArea;
            forceReferenceLengthPart = this->forceReferenceLength;
            forceReferenceLengthPart_B = this->forceReferenceLengthSpanWise;
            isGlobalCoordinate = false;
        }
        else
        {
            SimpleBC *bc = GlobalBoundaryCondition::GetBC(bcID);
            Data_Param *bcParamDB = bc->GetBCParamDataBase();
            bcParamDB->GetData("forceReferenceArea",   &forceReferenceAreaPart,   PHDOUBLE, 1);
            bcParamDB->GetData("forceReferenceLength", &forceReferenceLengthPart, PHDOUBLE, 1);
            bcParamDB->GetData("forceReferenceLengthSpanWise", &forceReferenceLengthPart_B, PHDOUBLE, 1);
            isGlobalCoordinate = true;
        }

        scf = 1.0 / forceReferenceAreaPart;
        vcm = 1.0 / (forceReferenceAreaPart * forceReferenceLengthPart);
        Sb = 1.0 / (forceReferenceAreaPart * forceReferenceLengthPart_B);

        CA_f[iPart] *= scf;
        CA_p[iPart] *= scf;
        CN_f[iPart] *= scf;
        CN_p[iPart] *= scf;
        CZ_f[iPart] *= scf;
        CZ_p[iPart] *= scf;

        cpx[iPart] *= scf;
        cpy[iPart] *= scf;
        cpz[iPart] *= scf;

        Cl_f[iPart] *= Sb;
        Cl_p[iPart] *= Sb;
        Cn_f[iPart] *= Sb;
        Cn_p[iPart] *= Sb;
        Cm_f[iPart] *= vcm;
        Cm_p[iPart] *= vcm;

        hingeMoment[iPart] *= vcm;
    }

    //* GRID_SIZE (F10.0) Number of grid nodes or cells
    //* MACH      (F10.2) Mach Number.
    //* ALPHA     (F10.3) Angle of Attack in Degrees.
    //* CL_TOT    (F10.3) Total Lift Coefficient.
    //* CD_TOT    (F10.5) Total Drag Coefficient.
    //* CD_PR     (F10.5) Surface-Pressure Integrated Drag Coefficient.
    //* CD_SF     (F10.5) Skin-Friction Integrated Drag Coefficient.
    //* CD-CL2/PA (F10.5) CD_TOT - (CL_TOT*CL_TOT)/(PI*AR)£¬Note: AR=9.5
    //* CM_TOT    (F10.4) Total Moment Coefficient.

    if(!cl_tot || !cd_tot || !cd_pr || !cd_sf || !cd_cl2pa || !xcp || !side)
    {
        cl_tot    = new RDouble[numberofTotalPart];
        cd_tot    = new RDouble[numberofTotalPart];
        cd_pr     = new RDouble[numberofTotalPart];
        cd_sf     = new RDouble[numberofTotalPart];
        cd_cl2pa  = new RDouble[numberofTotalPart];
        xcp       = new RDouble[numberofTotalPart];
        ycp       = new RDouble[numberofTotalPart];
        zcp       = new RDouble[numberofTotalPart];
        side      = new RDouble[numberofTotalPart];
    }

    const RDouble AR = 9.5; 
    for (int iPart = 0; iPart < numberofTotalPart; ++ iPart)
    {
        RDouble cfx, cfy, cfz;
        cfx = CA_f[iPart] + CA_p[iPart];
        cfy = CN_f[iPart] + CN_p[iPart];
        cfz = CZ_f[iPart] + CZ_p[iPart];

        RDouble cmx, cmy, cmz;
        cmx = Cl_f[iPart] + Cl_p[iPart];
        cmy = Cn_f[iPart] + Cn_p[iPart];
        cmz = Cm_f[iPart] + Cm_p[iPart];

        cl_tot[iPart]   = - sina *        cfx + cosa        * cfy;
        cd_tot[iPart]   = cosa   * cosb * cfx + sina * cosb * cfy + sinb * cfz;
        side[iPart]     = -cosa  * sinb * cfx - sina * sinb * cfy + cosb * cfz;
        cd_pr[iPart]    = cosa   * cosb * cpx[iPart] + sina * cosb * cpy[iPart] + sinb * cpz[iPart];
        cd_sf[iPart]    = cd_tot[iPart] - cd_pr[iPart];
        cd_cl2pa[iPart] = cd_tot[iPart] - ( cl_tot[iPart] * cl_tot[iPart] ) / ( PI * AR );

        const RDouble one = 1.0;
        xcp[iPart] = SIGN(one, cfy) * cmz/( ABS(cfy) + SMALL );
    }

    if (numberOfBody == 0)
    {
        return;
    }

    isGlobalCoordinate = false;
    for (int iPart = 0; iPart < numberofTotalPart; ++ iPart)
    {
        int bcID = partBCID[iPart];
        if (bcID == ALL_SOLID_WALL)
        {
            continue;
        }

        if (isGlobalCoordinate == false)
        {
            isGlobalCoordinate = true;
            continue;
        }

        int bodyIndex = GlobalBoundaryCondition::GetBodyIndex(bcID);

        CA_f_body[bodyIndex] += CA_f[iPart];
        CN_f_body[bodyIndex] += CN_f[iPart];
        CZ_f_body[bodyIndex] += CZ_f[iPart];
        Cl_f_body[bodyIndex] += Cl_f[iPart];
        Cn_f_body[bodyIndex] += Cn_f[iPart];
        Cm_f_body[bodyIndex] += Cm_f[iPart];
        CA_p_body[bodyIndex] += CA_p[iPart];
        CN_p_body[bodyIndex] += CN_p[iPart];
        CZ_p_body[bodyIndex] += CZ_p[iPart];
        Cl_p_body[bodyIndex] += Cl_p[iPart];
        Cn_p_body[bodyIndex] += Cn_p[iPart];
        Cm_p_body[bodyIndex] += Cm_p[iPart];

        isGlobalCoordinate = false;
    }
}

//! Bell modify 20120605.
void Force::DumpForce(int outnstep,std::ostringstream &oss)
{
    oss << setiosflags(ios::left);
    oss << setiosflags(ios::scientific);
    oss << setprecision(10);

    const int ALL = 0;
    oss << outnstep << "	";
    oss << cl_tot[ALL] << "	";
    oss << cd_tot[ALL] << "	";
    oss << side[ALL] << "	";
    oss << cd_pr[ALL] << "	";
    oss << cd_sf[ALL] << "	";
    oss << cd_cl2pa[ALL] << "	";
    oss << xcp[ALL] << "	";
    oss << CA_f[ALL] + CA_p[ALL] << "	";
    oss << CN_f[ALL] + CN_p[ALL] << "	";
    oss << CZ_f[ALL] + CZ_p[ALL] << "	";
    oss << Cl_f[ALL] + Cl_p[ALL] << "	";
    oss << Cn_f[ALL] + Cn_p[ALL] << "	";
    oss << Cm_f[ALL] + Cm_p[ALL] << "	";

    RDouble wallTime = PHSPACE::TIME_SPACE::GetWallTime();
    oss << wallTime;

    oss << "\n";
}

void Force::DumpPartForce(int outnstep, int partID, const string & partName, std::ostringstream &oss)
{
    oss << setiosflags(ios::left);
    oss << setiosflags(ios::scientific);
    oss << setprecision(10);

    oss << outnstep << "	";

    oss << cl_tot[partID] << "	";
    oss << cd_tot[partID] << "	";
    oss << side[partID] << "	";
    oss << cd_pr[partID] << "	";
    oss << cd_sf[partID] << "	";
    oss << cd_cl2pa[partID] << "	";
    oss << xcp[partID] << "	";
    oss << CA_f[partID] + CA_p[partID] << "	";
    oss << CN_f[partID] + CN_p[partID] << "	";
    oss << CZ_f[partID] + CZ_p[partID] << "	";
    oss << Cl_f[partID] + Cl_p[partID] << "	";
    oss << Cn_f[partID] + Cn_p[partID] << "	";
    oss << Cm_f[partID] + Cm_p[partID] << "	";

    oss << "\n";
}


void Force::DumpForce(ActionKey *actkey, int outnstep)
{
    CompForce();

    if (!actkey->IsNeedOpenFile())
    {
        return;
    }

    //! For total wall force.
    ParallelOpenFile(actkey);

    std::ostringstream oss;
    fstream &file = *(actkey->file);

    if (IfFileEmpty(file))
    {
        oss << "Title=\"THE AIRCOEF\"" << endl;
        oss << "Variables=" << endl;
        oss << "iter" << endl;
        oss << "<i>C<sub>L</sub></i>" << endl;
        oss << "<i>C<sub>D</sub></i>" << endl;
        oss << "<i>C<sub>Z</sub></i>" << endl;
        oss << "<i>C<sub>D_p</sub></i>" << endl;
        oss << "<i>C<sub>D_f</sub></i>" << endl;
        oss << "<i>C<sub>D</sub></i>-<i>C<sub>L</sub></i><sup>2</sup>/(<greek>p</greek><math>4</math>A<sub>r</sub>)</i>" << endl;
        oss << "<i>X<sub>cp</sub></i>" << endl;
        oss << "<i>C<sub>A</sub></i>" << endl;
        oss << "<i>C<sub>N</sub></i>" << endl;
        oss << "<i>C<sub>Z1</sub></i>" << endl;
        oss << "<i>C<sub>ml</sub></i>" << endl;
        oss << "<i>C<sub>mn</sub></i>" << endl;
        oss << "<i>C<sub>m</sub></i>" << endl;

        oss << "\"Walltime\"" << endl;
    }

    DumpForce(outnstep,oss);

    WriteASCIIFile(file, oss.str());

    ParallelCloseFile(actkey);

    //! For part force.
    bool isGlobalCoordinate = false;
    for (int iPart = 1; iPart < numberofTotalPart; ++ iPart)
    {
        int bcID = partBCID[iPart];

        SimpleBC *bc = GlobalBoundaryCondition::GetBC(bcID);
        const string &bcName = bc->GetBCName();

        string totalForceFileName  = actkey->filename;
        string partAircoefFileName = PHSPACE::AddSymbolToFileName(totalForceFileName, "_", bcName);

        if (isGlobalCoordinate == false)
        {
            isGlobalCoordinate = true;
            continue;
        }
        else
        {
            isGlobalCoordinate = false;
        }

        fstream file;
        ParallelOpenFile(file, partAircoefFileName, ios_base::out|ios_base::app);

        ostringstream oss;
        if (IfFileEmpty(file))
        {
            oss << "Title=\"Aerodynamic Force of " << bcName << "\"" << endl;
            oss << "Variables=" << endl;
            oss << "iter" << endl;

            oss << "<i>C<sub>L</sub></i>" << endl;
            oss << "<i>C<sub>D</sub></i>" << endl;
            oss << "<i>C<sub>Z</sub></i>" << endl;
            oss << "<i>C<sub>D_p</sub></i>" << endl;
            oss << "<i>C<sub>D_f</sub></i>" << endl;
            oss << "<i>C<sub>D</sub></i>-<i>C<sub>L</sub></i><sup>2</sup>/(<greek>p</greek><math>4</math>A<sub>r</sub>)</i>" << endl;
            oss << "<i>X<sub>cp</sub></i>" << endl;
            oss << "<i>C<sub>A</sub></i>" << endl;
            oss << "<i>C<sub>N</sub></i>" << endl;
            oss << "<i>C<sub>Z1</sub></i>" << endl;
            oss << "<i>C<sub>ml</sub></i>" << endl;
            oss << "<i>C<sub>mn</sub></i>" << endl;
            oss << "<i>C<sub>m</sub></i>" << endl;
        }
        DumpPartForce(outnstep, iPart, bcName, oss);
        WriteASCIIFile(file, oss.str());

        ParallelCloseFile(file);
    }

    DumpAllForceData(actkey, outnstep);
}

void Force::DumpAllForceData(ActionKey *actkey, int outnstep)
{
    string fileName = actkey->filename;

    string mainName, fileDir, extensionName;
    string globalFileName, localFileName;
    GetFileNameExtension(fileName, mainName, extensionName, ".");
    GetFileNameExtension(mainName, fileDir, extensionName, "/");

    ostringstream ossGlobalFileName, ossLocalFileName, ossGlobalData, ossLocalData;
    ossGlobalFileName << fileDir <<"/force_global.dat";
    ossLocalFileName << fileDir <<"/force_local.dat";
    globalFileName = ossGlobalFileName.str();
    localFileName = ossLocalFileName.str();

    string gridFile = GlobalDataBase::GetStrParaFromDB("gridfile");
    RDouble refMachNumber = GlobalDataBase::GetDoubleParaFromDB("refMachNumber");
    RDouble attackd = GlobalDataBase::GetDoubleParaFromDB("attackd");
    RDouble angleSlide = GlobalDataBase::GetDoubleParaFromDB("angleSlide");
    RDouble refReNumber = GlobalDataBase::GetDoubleParaFromDB("refReNumber");
    RDouble TorqueRefX = GlobalDataBase::GetDoubleParaFromDB("TorqueRefX");
    RDouble TorqueRefY = GlobalDataBase::GetDoubleParaFromDB("TorqueRefY");
    RDouble TorqueRefZ = GlobalDataBase::GetDoubleParaFromDB("TorqueRefZ");
    RDouble gridScaleFactor = GlobalDataBase::GetDoubleParaFromDB("gridScaleFactor");
    RDouble wallTemperature = GlobalDataBase::GetDoubleParaFromDB("wallTemperature");
    RDouble reynoldsReferenceLengthDimensional = GlobalDataBase::GetDoubleParaFromDB( "reynoldsReferenceLengthDimensional" );
    RDouble height = 0.0;
    int inflowParaType = GlobalDataBase::GetIntParaFromDB("inflowParaType");
    if (inflowParaType == 1)
    {
        height = GlobalDataBase::GetDoubleParaFromDB("height");
    }


    //! Write global force data file title
    ossGlobalData << "#     MeshFile                     :       " << gridFile << endl;
    if (inflowParaType == 1)
    {
        ossGlobalData << "#     Height                       :       " << height << "  km" << endl;
    }

    ossGlobalData << setiosflags(ios::scientific);
    ossGlobalData << setprecision(6);
    ossGlobalData << "#     Reynolds                     :       " << refReNumber / reynoldsReferenceLengthDimensional << "  1/m" << endl;
    ossGlobalData << resetiosflags(ios::scientific);

    ossGlobalData << "#     Mach                         :       " << refMachNumber << endl;
    ossGlobalData << "#     Alpha                        :       " << attackd << "  degree" << endl;
    ossGlobalData << "#     Beta                         :       " << angleSlide << "  degree" << endl;

    if (wallTemperature < 0.0)
    {
        ossGlobalData << "#     Wall Temperature             :       " << "Adiabatic Wall" << endl;
    }
    else
    {
        ossGlobalData << "#     Wall Temperature             :       " << wallTemperature << "  K" << endl;
    }

    ossGlobalData << "#     gridScaleFactor              :       " << gridScaleFactor << endl;
    ossGlobalData << "#     forceReferenceLength         :       " << this->forceReferenceLength * reynoldsReferenceLengthDimensional << "  m" << endl;
    ossGlobalData << "#     forceReferenceLengthSpanWise :       " << this->forceReferenceLengthSpanWise * reynoldsReferenceLengthDimensional << "  m" << endl;

    if (GetDim() == THREE_D)
    {
        ossGlobalData << "#     forceReferenceArea           :       " << this->forceReferenceArea * reynoldsReferenceLengthDimensional * reynoldsReferenceLengthDimensional << "  m2" << endl;
    }
    else
    {
        ossGlobalData << "#     forceReferenceArea           :       " << this->forceReferenceArea * reynoldsReferenceLengthDimensional << "  m^2" << endl;
    }

    ossGlobalData << "#     TorqueRefX                   :       " << TorqueRefX * reynoldsReferenceLengthDimensional << endl;
    ossGlobalData << "#     TorqueRefY                   :       " << TorqueRefY * reynoldsReferenceLengthDimensional << endl;
    ossGlobalData << "#     TorqueRefZ                   :       " << TorqueRefZ * reynoldsReferenceLengthDimensional << endl;
    ossGlobalData << "#     iter                         :       " << outnstep << endl << endl;
    ossGlobalData << "#        CA             CN            CZ1           Cml             Cmn            Cm"   << endl;
    ossGlobalData << "#        CA_p           CN_p          CZ1_p         Cml_p           Cmn_p          Cm_p" << endl;
    ossGlobalData << "#        CA_f           CN_f          CZ1_f         Cml_f           Cmn_f          Cm_f" << endl;


    //! Write local force data file title
    ossLocalData << "#     MeshFile             :       " << gridFile << endl;
    if (inflowParaType == 1)
    {
        ossLocalData << "#     Height               :       " << height << "  km" << endl;
    }

    ossLocalData << setiosflags(ios::scientific);
    ossLocalData << setprecision(6);
    ossLocalData << "#     Reynolds             :       " << refReNumber / reynoldsReferenceLengthDimensional << "  1/m" << endl;
    ossLocalData << resetiosflags(ios::scientific);

    ossLocalData << "#     Mach                 :       " << refMachNumber << endl;
    ossLocalData << "#     Alpha                :       " << attackd << "  degree" << endl;
    ossLocalData << "#     Beta                 :       " << angleSlide << "  degree" << endl;

    if (wallTemperature < 0.0)
    {
        ossLocalData << "#     Wall Temperature     :       " << "Adiabatic Wall" << endl;
    }
    else
    {
        ossLocalData << "#     Wall Temperature     :       " << wallTemperature << "  K" << endl;
    }
    

    ossLocalData << "#     gridScaleFactor      :       " << gridScaleFactor << endl;
    ossLocalData << "#     iter                 :       " << outnstep << endl << endl;
    ossLocalData << "#        CA             CN            CZ1           Cml             Cmn            Cm"   << endl;
    ossLocalData << "#        CA_p           CN_p          CZ1_p         Cml_p           Cmn_p          Cm_p" << endl;
    ossLocalData << "#        CA_f           CN_f          CZ1_f         Cml_f           Cmn_f          Cm_f" << endl;
    ossLocalData << "#        CMh" << endl;


    bool isGlobalData = true;
    //! Write force data file
    for (int iPart = 0; iPart < numberofTotalPart; ++ iPart)
    {
        int bcID = partBCID[iPart];

        string partName;
        if (bcID == -1)
        {
            partName = "whole";
        }
        else
        {
            SimpleBC *bc = GlobalBoundaryCondition::GetBC(bcID);
            const string &bcName = bc->GetBCName();
            partName = bcName;
        }

        if (isGlobalData == true)
        {
            DumpAllGlobalForceData(partName, iPart, ossGlobalData);
            isGlobalData = false;
        }
        else
        {
            DumpAllLocalForceData(partName, iPart, ossLocalData);
            isGlobalData = true;
        }
    }

    for (int iBody = 0; iBody < numberOfBody; ++ iBody)
    {
        DumpBodyForceData(iBody, ossGlobalData);
    }

    //! Write global force data file
    fstream globalFile;
    ParallelOpenFile(globalFile, globalFileName, ios_base::out|ios_base::trunc);
    WriteASCIIFile(globalFile, ossGlobalData.str());
    ParallelCloseFile(globalFile);

    //! Write local force data file
    fstream localFile;
    ParallelOpenFile(localFile, localFileName, ios_base::out|ios_base::trunc);
    WriteASCIIFile(localFile, ossLocalData.str());
    ParallelCloseFile(localFile);
}

void Force::DumpAllLocalForceData(const string &partName, int partID, std::ostringstream &oss)
{
    RDouble TorqueRefX, TorqueRefY, TorqueRefZ, forceReferenceAreaPart, forceReferenceLengthPart, forceReferenceLengthSpanWisePart;
    int bcID = partBCID[partID];
    SimpleBC *bc = GlobalBoundaryCondition::GetBC(bcID);
    Data_Param *bcParamDB = bc->GetBCParamDataBase();
    bcParamDB->GetData("forceReferenceArea",   &forceReferenceAreaPart,   PHDOUBLE, 1);
    bcParamDB->GetData("forceReferenceLength", &forceReferenceLengthPart, PHDOUBLE, 1);
    bcParamDB->GetData("forceReferenceLengthSpanWise", &forceReferenceLengthSpanWisePart, PHDOUBLE, 1);
    bcParamDB->GetData("TorqueRefX", &TorqueRefX, PHDOUBLE, 1);
    bcParamDB->GetData("TorqueRefY", &TorqueRefY, PHDOUBLE, 1);
    bcParamDB->GetData("TorqueRefZ", &TorqueRefZ, PHDOUBLE, 1);

    RDouble reynoldsReferenceLengthDimensional = GlobalDataBase::GetDoubleParaFromDB("reynoldsReferenceLengthDimensional");

    oss << resetiosflags(ios::scientific);

    oss <<endl;
    oss << "*"<<partName <<endl;
    oss << "#forceReferenceLength         :    " << forceReferenceLengthPart * reynoldsReferenceLengthDimensional << "  m" << endl;
    oss << "#forceReferenceLengthSpanWise :    " << forceReferenceLengthSpanWisePart * reynoldsReferenceLengthDimensional << "  m" <<endl;

    if (GetDim() == THREE_D)
    {
        oss << "#forceReferenceArea           :    " << forceReferenceAreaPart * reynoldsReferenceLengthDimensional * reynoldsReferenceLengthDimensional << "  m2" <<endl;
    }
    else
    {
        oss << "#forceReferenceArea           :    " << forceReferenceAreaPart * reynoldsReferenceLengthDimensional << "  m^2" <<endl;
    }

    oss << "#TorqueRefX                   :    " << TorqueRefX * reynoldsReferenceLengthDimensional << endl;
    oss << "#TorqueRefY                   :    " << TorqueRefY * reynoldsReferenceLengthDimensional << endl;
    oss << "#TorqueRefZ                   :    " << TorqueRefZ * reynoldsReferenceLengthDimensional << endl;

    oss << setiosflags(ios::right);
    oss << setiosflags(ios::scientific);
    oss << setprecision(10);

    oss << CA_f[partID] + CA_p[partID] << "	";
    oss << CN_f[partID] + CN_p[partID] << "	";
    oss << CZ_f[partID] + CZ_p[partID] << "	";
    oss << Cl_f[partID] + Cl_p[partID] << "	";
    oss << Cn_f[partID] + Cn_p[partID] << "	";
    oss << Cm_f[partID] + Cm_p[partID] << endl;

    oss << CA_p[partID] << "	";
    oss << CN_p[partID] << "	";
    oss << CZ_p[partID] << "	";
    oss << Cl_p[partID] << "	";
    oss << Cn_p[partID] << "	";
    oss << Cm_p[partID] << endl;

    oss << CA_f[partID] << "	";
    oss << CN_f[partID] << "	";
    oss << CZ_f[partID] << "	";
    oss << Cl_f[partID] << "	";
    oss << Cn_f[partID] << "	";
    oss << Cm_f[partID] << endl;

    oss << hingeMoment[partID] << endl;
}

void Force::DumpAllGlobalForceData(const string &partName, int partID, std::ostringstream &oss)
{
    oss << setiosflags(ios::right);
    oss << setiosflags(ios::scientific);
    oss << setprecision(10);

    oss <<endl;
    oss << "*"<<partName <<endl;

    oss << CA_f[partID] + CA_p[partID] << "	";
    oss << CN_f[partID] + CN_p[partID] << "	";
    oss << CZ_f[partID] + CZ_p[partID] << "	";
    oss << Cl_f[partID] + Cl_p[partID] << "	";
    oss << Cn_f[partID] + Cn_p[partID] << "	";
    oss << Cm_f[partID] + Cm_p[partID] << endl;

    oss << CA_p[partID] << "	";
    oss << CN_p[partID] << "	";
    oss << CZ_p[partID] << "	";
    oss << Cl_p[partID] << "	";
    oss << Cn_p[partID] << "	";
    oss << Cm_p[partID] << endl;

    oss << CA_f[partID] << "	";
    oss << CN_f[partID] << "	";
    oss << CZ_f[partID] << "	";
    oss << Cl_f[partID] << "	";
    oss << Cn_f[partID] << "	";
    oss << Cm_f[partID] << endl;
}

void Force::DumpBodyForceData(int bodyID, std::ostringstream &oss)
{
    oss << setiosflags(ios::right);
    oss << setiosflags(ios::scientific);
    oss << setprecision(10);

    oss << endl;
    oss << "*" << wallNameOfEachBody[bodyID][0] << "(";
    oss << wallNameOfEachBody[bodyID][1];

    for (int ibc = 2; ibc < wallNameOfEachBody[bodyID].size(); ++ ibc)
    {
        oss << " + "<<wallNameOfEachBody[bodyID][ibc];
    }
    oss << ")"<<endl;

    oss << CA_f_body[bodyID] + CA_p_body[bodyID] << "	";
    oss << CN_f_body[bodyID] + CN_p_body[bodyID] << "	";
    oss << CZ_f_body[bodyID] + CZ_p_body[bodyID] << "	";
    oss << Cl_f_body[bodyID] + Cl_p_body[bodyID] << "	";
    oss << Cn_f_body[bodyID] + Cn_p_body[bodyID] << "	";
    oss << Cm_f_body[bodyID] + Cm_p_body[bodyID] << endl;

    oss << CA_p_body[bodyID] << "	";
    oss << CN_p_body[bodyID] << "	";
    oss << CZ_p_body[bodyID] << "	";
    oss << Cl_p_body[bodyID] << "	";
    oss << Cn_p_body[bodyID] << "	";
    oss << Cm_p_body[bodyID] << endl;

    oss << CA_f_body[bodyID] << "	";
    oss << CN_f_body[bodyID] << "	";
    oss << CZ_f_body[bodyID] << "	";
    oss << Cl_f_body[bodyID] << "	";
    oss << Cn_f_body[bodyID] << "	";
    oss << Cm_f_body[bodyID] << endl;
}

void Force::DumpData(DataContainer * data)
{
    data->moveto_begin();

    PHWrite(data, numberofSolidWallPart);
    PHWrite(data, numberofTotalPart);

    PHWrite(data, partBCID, numberofTotalPart);

    PHWrite(data, cpx, numberofTotalPart);
    PHWrite(data, cpy, numberofTotalPart);
    PHWrite(data, cpz, numberofTotalPart);

    PHWrite(data, cl_tot, numberofTotalPart);
    PHWrite(data, cd_tot, numberofTotalPart);
    PHWrite(data, cd_pr, numberofTotalPart);
    PHWrite(data, cd_sf, numberofTotalPart);
    PHWrite(data, cd_cl2pa, numberofTotalPart);
    PHWrite(data, xcp, numberofTotalPart);
    PHWrite(data, side, numberofTotalPart);
    PHWrite(data, hingeMoment, numberofTotalPart);

}

    bool AerodynamicForceComponent::initialize = false;
    vector< std::set< string > >  AerodynamicForceComponent::components;

    AerodynamicForceComponent::AerodynamicForceComponent()
    {
    }

    AerodynamicForceComponent::~AerodynamicForceComponent()
    {
    }

    void AerodynamicForceComponent::Initialize()
    {
        if (!AerodynamicForceComponent::initialize)
        {
            int numberOfAerodynamicForceComponents = GlobalDataBase::GetIntParaFromDB("numberOfAerodynamicForceComponents");
            AerodynamicForceComponent::components.resize(numberOfAerodynamicForceComponents);
            int nTBCWall = GlobalBoundaryCondition::GetNumberofSolidWallPart();
            SimpleBC* boundaryCondition = 0;
            int iComponent = 0;
            int nTBC = GlobalBoundaryCondition::GetNumberofBC();
            if (numberOfAerodynamicForceComponents != nTBCWall)
            {
                TK_Exit::ExceptionExit("Error in AerodynamicForceComponent::initialize : numberOfAerodynamicForceComponents != nTBCWall");

            }
            for (int iBC = 0; iBC < nTBC; ++iBC)
            {
                if (GlobalBoundaryCondition::IsSolidWallBC(iBC))
                {
                    boundaryCondition = GlobalBoundaryCondition::GetBC(iBC);
                    string bcName = boundaryCondition->GetBCName();
                    AerodynamicForceComponent::components[iComponent++].insert(bcName);
                }
            }

            AerodynamicForceComponent::initialize = true;
        }
    }

    bool AerodynamicForceComponent::FaceInAerodynamicForceComponents(string componentName, int aerodynamicForceComponentIndex)
    {
        AerodynamicForceComponent::Initialize();
        set< string >::iterator iter = AerodynamicForceComponent::components[aerodynamicForceComponentIndex].find(componentName);
        return iter != AerodynamicForceComponent::components[aerodynamicForceComponentIndex].end();
    }

    bool FaceInAerodynamicForceComponents(string componentName, int aerodynamicForceComponentIndex)
    {
        return AerodynamicForceComponent::FaceInAerodynamicForceComponents(componentName, aerodynamicForceComponentIndex);
    }

    AerodynamicForceManager aerodynamicForceManager;

    void InitializeAerodynamicForce()
    {
        AerodynamicForceManager* aerodynamicForceManager = PHSPACE::GetAerodynamicForceManager();
        aerodynamicForceManager->InitializeAerodynamicForceManager();
    }

    void AddAerodynamicForce(ActionKey* actionKey)
    {
        AerodynamicForceManager* aerodynamicForceManager = PHSPACE::GetAerodynamicForceManager();
        aerodynamicForceManager->AddAerodynamicForce(actionKey);
    }

    void DumpAerodynamicForceCoefficient(ActionKey* actionKey)
    {
        AerodynamicForceManager* aerodynamicForceManager = PHSPACE::GetAerodynamicForceManager();
        aerodynamicForceManager->DumpAerodynamicForceCoefficient(actionKey);
    }

    void ComputeAerodynamicForceCoefficient()
    {
        AerodynamicForceManager* aerodynamicForceManager = PHSPACE::GetAerodynamicForceManager();
        aerodynamicForceManager->ComputeAerodynamicForceCoefficient();
    }

    void CollectAerodynamicForceOfAllProcessors()
    {
        AerodynamicForceManager* aerodynamicForceManager = PHSPACE::GetAerodynamicForceManager();
        aerodynamicForceManager->CollectAerodynamicForceOfAllProcessors();
    }

    BasicAerodynamicForce** CreateBasicAerodynamicForceContainer()
    {
        int numberOfAerodynamicForceComponents = GlobalDataBase::GetIntParaFromDB("numberOfAerodynamicForceComponents");
        int numberOfAerodynamicForces = numberOfAerodynamicForceComponents + 1;

        BasicAerodynamicForce** basicAerodynamicForceContainer = new BasicAerodynamicForce * [numberOfAerodynamicForces];

        for (int iAerodynamicForces = 0; iAerodynamicForces < numberOfAerodynamicForces; ++iAerodynamicForces)
        {
            basicAerodynamicForceContainer[iAerodynamicForces] = new BasicAerodynamicForce();
        }

        return basicAerodynamicForceContainer;
    }

    vector< BasicAerodynamicForce* >* CreateBasicAerodynamicForceVector()
    {
        int numberOfAerodynamicForceComponents = GlobalDataBase::GetIntParaFromDB("numberOfAerodynamicForceComponents");
        int numberOfAerodynamicForces = numberOfAerodynamicForceComponents + 1;

        vector< BasicAerodynamicForce* >* basicAerodynamicForceVector = new vector< BasicAerodynamicForce* >(numberOfAerodynamicForces);

        for (int iAerodynamicForces = 0; iAerodynamicForces < numberOfAerodynamicForces; ++iAerodynamicForces)
        {
            (*basicAerodynamicForceVector)[iAerodynamicForces] = new BasicAerodynamicForce();
        }

        return basicAerodynamicForceVector;
    }

    void FreeBasicAerodynamicForceVector(vector< BasicAerodynamicForce* >* basicAerodynamicForceVector)
    {
        int numberOfAerodynamicForces = basicAerodynamicForceVector->size();

        for (int iAerodynamicForces = 0; iAerodynamicForces < numberOfAerodynamicForces; ++iAerodynamicForces)
        {
            delete (*basicAerodynamicForceVector)[iAerodynamicForces];
        }

        delete  basicAerodynamicForceVector;
    }

    void FreeBasicAerodynamicForceContainer(BasicAerodynamicForce** basicAerodynamicForceContainer)
    {
        int numberOfAerodynamicForceComponents = GlobalDataBase::GetIntParaFromDB("numberOfAerodynamicForceComponents");
        int numberOfAerodynamicForces = numberOfAerodynamicForceComponents + 1;

        for (int iAerodynamicForces = 0; iAerodynamicForces < numberOfAerodynamicForces; ++iAerodynamicForces)
        {
            delete basicAerodynamicForceContainer[iAerodynamicForces];
        }

        delete [] basicAerodynamicForceContainer;
    }

    void WriteAerodynamicForceToActionKey(ActionKey* actionKey, BasicAerodynamicForce** basicAerodynamicForceContainer)
    {
        int numberOfAerodynamicForceComponents = GlobalDataBase::GetIntParaFromDB("numberOfAerodynamicForceComponents");
        int numberOfAerodynamicForces = numberOfAerodynamicForceComponents + 1;

        DataContainer* dataContainer = actionKey->GetData();
        dataContainer->moveto_begin();

        for (int iAerodynamicForces = 0; iAerodynamicForces < numberOfAerodynamicForces; ++iAerodynamicForces)
        {
            basicAerodynamicForceContainer[iAerodynamicForces]->WriteAerodynamicForceToDataContainer(dataContainer);
        }
    }

    void WriteAerodynamicForceToActionKey(ActionKey* actionKey, vector< BasicAerodynamicForce* >* basicAerodynamicForceVector)
    {
        int numberOfAerodynamicForceComponents = GlobalDataBase::GetIntParaFromDB("numberOfAerodynamicForceComponents");
        int numberOfAerodynamicForces = numberOfAerodynamicForceComponents + 1;

        DataContainer* dataContainer = actionKey->GetData();
        dataContainer->moveto_begin();

        for (int iAerodynamicForces = 0; iAerodynamicForces < numberOfAerodynamicForces; ++iAerodynamicForces)
        {
            (*basicAerodynamicForceVector)[iAerodynamicForces]->WriteAerodynamicForceToDataContainer(dataContainer);
        }
    }

    BasicAerodynamicForce::BasicAerodynamicForce()
    {
        this->ZeroBasicAerodynamicForce();
    }

    BasicAerodynamicForce::~BasicAerodynamicForce()
    {
    }

    void BasicAerodynamicForce::ZeroBasicAerodynamicForce()
    {
        forceX = 0.0;
        forceY = 0.0;
        forceZ = 0.0;

        forceXByPressure = 0.0;
        forceYByPressure = 0.0;
        forceZByPressure = 0.0;

        forceXByViscosity = 0.0;
        forceYByViscosity = 0.0;
        forceZByViscosity = 0.0;

        momentX = 0.0;
        momentY = 0.0;
        momentZ = 0.0;

        momentXByPressure = 0.0;
        momentYByPressure = 0.0;
        momentZByPressure = 0.0;

        momentXByViscosity = 0.0;
        momentYByViscosity = 0.0;
        momentZByViscosity = 0.0;
    }

    void BasicAerodynamicForce::SetMomentReferencePoint(RDouble momentReferenceX, RDouble momentReferenceY, RDouble momentReferenceZ)
    {
        this->momentReferenceX = momentReferenceX;
        this->momentReferenceY = momentReferenceY;
        this->momentReferenceZ = momentReferenceZ;
    }

    void BasicAerodynamicForce::SetPressureAerodynamicForce(RDouble forceXByPressure, RDouble forceYByPressure, RDouble forceZByPressure)
    {
        this->forceXByPressure = forceXByPressure;
        this->forceYByPressure = forceYByPressure;
        this->forceZByPressure = forceZByPressure;
    }

    void BasicAerodynamicForce::SetViscousAerodynamicForce(RDouble forceXByViscosity, RDouble forceYByViscosity, RDouble forceZByViscosity)
    {
        this->forceXByViscosity = forceXByViscosity;
        this->forceYByViscosity = forceYByViscosity;
        this->forceZByViscosity = forceZByViscosity;
    }

    void BasicAerodynamicForce::ComputeResultantAerodynamicForce()
    {
        forceX = forceXByPressure + forceXByViscosity;
        forceY = forceYByPressure + forceYByViscosity;
        forceZ = forceZByPressure + forceZByViscosity;
    }

    void BasicAerodynamicForce::ComputeResultantAerodynamicMoment()
    {
        momentX = momentXByPressure + momentXByViscosity;
        momentY = momentYByPressure + momentYByViscosity;
        momentZ = momentZByPressure + momentZByViscosity;
    }

    void BasicAerodynamicForce::SetPressureAerodynamicMoment(RDouble faceCenterX, RDouble faceCenterY, RDouble faceCenterZ)
    {
        momentXByPressure = (faceCenterY - momentReferenceY) * forceZByPressure - (faceCenterZ - momentReferenceZ) * forceYByPressure;
        momentYByPressure = (faceCenterZ - momentReferenceZ) * forceXByPressure - (faceCenterX - momentReferenceX) * forceZByPressure;
        momentZByPressure = (faceCenterX - momentReferenceX) * forceYByPressure - (faceCenterY - momentReferenceY) * forceXByPressure;
    }

    void BasicAerodynamicForce::SetViscousAerodynamicMoment(RDouble faceCenterX, RDouble faceCenterY, RDouble faceCenterZ)
    {
        momentXByViscosity = (faceCenterY - momentReferenceY) * forceZByViscosity - (faceCenterZ - momentReferenceZ) * forceYByViscosity;
        momentYByViscosity = (faceCenterZ - momentReferenceZ) * forceXByViscosity - (faceCenterX - momentReferenceX) * forceZByViscosity;
        momentZByViscosity = (faceCenterX - momentReferenceX) * forceYByViscosity - (faceCenterY - momentReferenceY) * forceXByViscosity;
    }

    void BasicAerodynamicForce::ComputeAerodynamicMoment(RDouble faceCenterX, RDouble faceCenterY, RDouble faceCenterZ)
    {
        momentX = (faceCenterY - momentReferenceY) * forceZ - (faceCenterZ - momentReferenceZ) * forceY;
        momentY = (faceCenterZ - momentReferenceZ) * forceX - (faceCenterX - momentReferenceX) * forceZ;
        momentZ = (faceCenterX - momentReferenceX) * forceY - (faceCenterY - momentReferenceY) * forceX;
    }

    void BasicAerodynamicForce::AddAerodynamicForce(BasicAerodynamicForce* basicAerodynamicForce)
    {
        this->forceX += basicAerodynamicForce->forceX;
        this->forceY += basicAerodynamicForce->forceY;
        this->forceZ += basicAerodynamicForce->forceZ;

        this->forceXByPressure += basicAerodynamicForce->forceXByPressure;
        this->forceYByPressure += basicAerodynamicForce->forceYByPressure;
        this->forceZByPressure += basicAerodynamicForce->forceZByPressure;

        this->forceXByViscosity += basicAerodynamicForce->forceXByViscosity;
        this->forceYByViscosity += basicAerodynamicForce->forceYByViscosity;
        this->forceZByViscosity += basicAerodynamicForce->forceZByViscosity;

        this->momentX += basicAerodynamicForce->momentX;
        this->momentY += basicAerodynamicForce->momentY;
        this->momentZ += basicAerodynamicForce->momentZ;

        this->momentXByPressure += basicAerodynamicForce->momentXByPressure;
        this->momentYByPressure += basicAerodynamicForce->momentYByPressure;
        this->momentZByPressure += basicAerodynamicForce->momentZByPressure;

        this->momentXByViscosity += basicAerodynamicForce->momentXByViscosity;
        this->momentYByViscosity += basicAerodynamicForce->momentYByViscosity;
        this->momentZByViscosity += basicAerodynamicForce->momentZByViscosity;
    }

    void BasicAerodynamicForce::WriteAerodynamicForceToDataContainer(DataContainer* dataContainer)
    {
        PHWrite(dataContainer, forceX);
        PHWrite(dataContainer, forceY);
        PHWrite(dataContainer, forceZ);

        PHWrite(dataContainer, forceXByPressure);
        PHWrite(dataContainer, forceYByPressure);
        PHWrite(dataContainer, forceZByPressure);

        PHWrite(dataContainer, forceXByViscosity);
        PHWrite(dataContainer, forceYByViscosity);
        PHWrite(dataContainer, forceZByViscosity);

        PHWrite(dataContainer, momentX);
        PHWrite(dataContainer, momentY);
        PHWrite(dataContainer, momentZ);

        PHWrite(dataContainer, momentXByPressure);
        PHWrite(dataContainer, momentYByPressure);
        PHWrite(dataContainer, momentZByPressure);

        PHWrite(dataContainer, momentXByViscosity);
        PHWrite(dataContainer, momentYByViscosity);
        PHWrite(dataContainer, momentZByViscosity);
    }

    void BasicAerodynamicForce::ReadAerodynamicForceFromDataContainer(DataContainer* dataContainer)
    {
        PHRead(dataContainer, forceX);
        PHRead(dataContainer, forceY);
        PHRead(dataContainer, forceZ);

        PHRead(dataContainer, forceXByPressure);
        PHRead(dataContainer, forceYByPressure);
        PHRead(dataContainer, forceZByPressure);

        PHRead(dataContainer, forceXByViscosity);
        PHRead(dataContainer, forceYByViscosity);
        PHRead(dataContainer, forceZByViscosity);

        PHRead(dataContainer, momentX);
        PHRead(dataContainer, momentY);
        PHRead(dataContainer, momentZ);

        PHRead(dataContainer, momentXByPressure);
        PHRead(dataContainer, momentYByPressure);
        PHRead(dataContainer, momentZByPressure);

        PHRead(dataContainer, momentXByViscosity);
        PHRead(dataContainer, momentYByViscosity);
        PHRead(dataContainer, momentZByViscosity);
    }

    BasicAleForce::BasicAleForce()
    {
        this->ZeroBasicAleForce();
    }

    BasicAleForce::~BasicAleForce()
    {
    }

    void BasicAleForce::ZeroBasicAleForce()
    {
        forceX = 0.0;
        forceY = 0.0;
        forceZ = 0.0;

        momentX = 0.0;
        momentY = 0.0;
        momentZ = 0.0;
    }

    void BasicAleForce::SetBasicAleForceByAerodynamicForce(AerodynamicForce* aerodynamicForce)
    {
        forceX = aerodynamicForce->GetAerodynamicForceX();
        forceY = aerodynamicForce->GetAerodynamicForceY();
        forceZ = aerodynamicForce->GetAerodynamicForceZ();

        aerodynamicForce->GetMomentByOriginOfCoordinate(momentX, momentY, momentZ);
    }

    void BasicAleForce::ReadBasicAleForceFromBasicAleForce(BasicAleForce* basicAleForce)
    {
        (*this) = (*basicAleForce);
    }

    void BasicAleForce::ReadBasicAleForceFromFile(fstream& file)
    {
        PHSPACE::PHRead(file, forceX);
        PHSPACE::PHRead(file, forceY);
        PHSPACE::PHRead(file, forceZ);

        PHSPACE::PHRead(file, momentX);
        PHSPACE::PHRead(file, momentY);
        PHSPACE::PHRead(file, momentZ);

        PHSPACE::PHRead(file, powerConsumption);
    }

    void BasicAleForce::WriteBasicAleForceToFile(fstream& file)
    {
        PHSPACE::PHWrite(file, forceX);
        PHSPACE::PHWrite(file, forceY);
        PHSPACE::PHWrite(file, forceZ);

        PHSPACE::PHWrite(file, momentX);
        PHSPACE::PHWrite(file, momentY);
        PHSPACE::PHWrite(file, momentZ);

        PHSPACE::PHWrite(file, powerConsumption);
    }

    AleForceOfBody::AleForceOfBody()
    {
        bodyIndex = 0;
        numberOfTimeLayers = 0;
    }

    AleForceOfBody::~AleForceOfBody()
    {
        this->DeallocateMemory();
    }

    void AleForceOfBody::AllocateMemory()
    {
        basicAleForceContainer.resize(numberOfTimeLayers);
        for (int iTimeLayer = 0; iTimeLayer < numberOfTimeLayers; ++iTimeLayer)
        {
            basicAleForceContainer[iTimeLayer] = new BasicAleForce();
        }
    }

    void AleForceOfBody::DeallocateMemory()
    {
        for (int iTimeLayer = 0; iTimeLayer < numberOfTimeLayers; ++iTimeLayer)
        {
            delete basicAleForceContainer[iTimeLayer];
        }
    }

    BasicAleForce* AleForceOfBody::GetBasicAleForce(int timeLayerIndex)
    {
        return basicAleForceContainer[timeLayerIndex];
    }

    void AleForceOfBody::DumpRestartAleForceOfBodyContainer(fstream& file)
    {
        int numberOfTimeLayers = this->GetNumberOfTimeLayers();
        for (int iTimeLayer = 0; iTimeLayer < numberOfTimeLayers; ++iTimeLayer)
        {
            BasicAleForce* basicAleForce = this->GetBasicAleForce(iTimeLayer);
            basicAleForce->WriteBasicAleForceToFile(file);
        }
    }

    void AleForceOfBody::InitializeAleForceOfBodyContainer()
    {
        int numberOfTimeLayers = this->GetNumberOfTimeLayers();
        int bodyIndex = this->GetBodyIndex();
        for (int iTimeLayer = 0; iTimeLayer < numberOfTimeLayers; ++iTimeLayer)
        {
            BasicAleForce* basicAleForce = this->GetBasicAleForce(iTimeLayer);
            AerodynamicForce* aerodynamicForce = PHSPACE::GetCurrentAerodynamicForce(bodyIndex);
            basicAleForce->SetBasicAleForceByAerodynamicForce(aerodynamicForce);
        }
    }

    void AleForceOfBody::ReadAleForceOfBodyContainer()
    {
        fstream& aleFile = PHSPACE::GetCurrentAleFile();

        int numberOfTimeLayers = this->GetNumberOfTimeLayers();
        for (int iTimeLayer = 0; iTimeLayer < numberOfTimeLayers; ++iTimeLayer)
        {
            BasicAleForce* basicAleForce = this->GetBasicAleForce(iTimeLayer);
            fstream& file = aleFile;
            basicAleForce->ReadBasicAleForceFromFile(file);
        }
    }

    void AleForceOfBody::UpdateAleForceOfBodyOfNewTimeLayer()
    {
        int iTimeLayer = 0;
        int bodyIndex = this->GetBodyIndex();

        AerodynamicForce* aerodynamicForce = PHSPACE::GetCurrentAerodynamicForce(bodyIndex);

        BasicAleForce* basicAleForce = this->GetBasicAleForce(iTimeLayer);
        basicAleForce->SetBasicAleForceByAerodynamicForce(aerodynamicForce);
    }

    void AleForceOfBody::UpdateUnsteadyAleForceTerm()
    {
        int iTimeLevel0 = 0;
        int iTimeLevel1 = 1;
        int iTimeLevel2 = 2;

        BasicAleForce* basicAleForce = this->GetBasicAleForce(iTimeLevel0);
        BasicAleForce* basicAleForceNTimeLevel = this->GetBasicAleForce(iTimeLevel1);
        BasicAleForce* basicAleForceN1TimeLevel = this->GetBasicAleForce(iTimeLevel2);

        PHSPACE::CopyBasicAleForce(basicAleForceN1TimeLevel, basicAleForceNTimeLevel);
        PHSPACE::CopyBasicAleForce(basicAleForceNTimeLevel, basicAleForce);
    }

    AleForceManager::AleForceManager()
    {
        numberOfMovingBodies = 0;

        dimensionalForceCoefficient = 1.0;
        dimensionalMomentCoefficient = 1.0;
        dimensionalPowerCoefficient = 1.0;
        numberOfTimeLayers = 0;
    }

    AleForceManager::~AleForceManager()
    {
        this->DeallocateAleForceOfBodyContainer();
    }

    void AleForceManager::Initialize()
    {
        this->AllocateAleForceOfBodyContainer();
        this->InitializeAleForceOfBodyContainer();
    }

    void AleForceManager::Restart()
    {
        this->ReadAleForceOfBodyContainer();
    }

    void AleForceManager::InitializeDimensionalReferenceVariables()
    {
        RDouble referenceDensity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalDensity");
        RDouble referenceVelocity = GlobalDataBase::GetDoubleParaFromDB("refDimensionalVelocity");
        RDouble referenceLength = GlobalDataBase::GetDoubleParaFromDB("forceReferenceLength");

        int geometricDimension = GetDim();

        RDouble dynamicPressure = referenceDensity * PHSPACE::SQR(referenceVelocity);

        dimensionalForceCoefficient = half * dynamicPressure * PHSPACE::SQR(referenceLength);
        dimensionalMomentCoefficient = half * dynamicPressure * PHSPACE::POWER3(referenceLength);
        dimensionalPowerCoefficient = dimensionalForceCoefficient * referenceVelocity;

        if (geometricDimension == TWO_D)
        {
            dimensionalForceCoefficient = half * dynamicPressure * referenceLength;
            dimensionalMomentCoefficient = half * dynamicPressure * PHSPACE::SQR(referenceLength);
        }
    }

    void AleForceManager::AllocateAleForceOfBodyContainer()
    {
        numberOfTimeLayers = 3;
        numberOfMovingBodies = GlobalDataBase::GetIntParaFromDB("numberOfMovingBodies");
        this->InitializeDimensionalReferenceVariables();

        aleForceOfBodyContainer.resize(numberOfMovingBodies);

        for (int iMovingBody = 0; iMovingBody < numberOfMovingBodies; ++iMovingBody)
        {
            AleForceOfBody* aleForceOfBody = new AleForceOfBody();
            aleForceOfBodyContainer[iMovingBody] = aleForceOfBody;

            aleForceOfBody->SetBodyIndex(iMovingBody);
            aleForceOfBody->SetNumberOfTimeLayers(numberOfTimeLayers);
            aleForceOfBody->AllocateMemory();
        }
    }

    void AleForceManager::DeallocateAleForceOfBodyContainer()
    {
        for (int iMovingBody = 0; iMovingBody < numberOfMovingBodies; ++iMovingBody)
        {
            delete aleForceOfBodyContainer[iMovingBody];
        }
    }

    void AleForceManager::InitializeAleForceOfBodyContainer()
    {
        for (int iMovingBody = 0; iMovingBody < numberOfMovingBodies; ++iMovingBody)
        {
            AleForceOfBody* aleForceOfBody = this->GetAleForceOfBody(iMovingBody);

            aleForceOfBody->InitializeAleForceOfBodyContainer();
        }
    }

    void AleForceManager::ReadAleForceOfBodyContainer()
    {
        for (int iMovingBody = 0; iMovingBody < numberOfMovingBodies; ++iMovingBody)
        {
            AleForceOfBody* aleForceOfBody = this->GetAleForceOfBody(iMovingBody);

            aleForceOfBody->ReadAleForceOfBodyContainer();
        }
    }

    void AleForceManager::Run()
    {
        for (int iMovingBody = 0; iMovingBody < numberOfMovingBodies; ++iMovingBody)
        {
            aleForceOfBodyContainer[iMovingBody]->UpdateAleForceOfBodyOfNewTimeLayer();
        }
    }

    void AleForceManager::Post()
    {
        int numberOfMovingBodies = GlobalDataBase::GetIntParaFromDB("numberOfMovingBodies");

        for (int iMovingBody = 0; iMovingBody < numberOfMovingBodies; ++iMovingBody)
        {
            AleForceOfBody* aleForceOfBody = this->GetAleForceOfBody(iMovingBody);
            aleForceOfBody->UpdateUnsteadyAleForceTerm();
        }
    }

    void AleForceManager::DumpRestartAleForceOfBodyContainer(fstream& file)
    {
        for (int iMovingBody = 0; iMovingBody < numberOfMovingBodies; ++iMovingBody)
        {
            AleForceOfBody* aleForceOfBody = this->GetAleForceOfBody(iMovingBody);
            aleForceOfBody->DumpRestartAleForceOfBodyContainer(file);
        }
    }

    BasicAleForce* AleForceManager::GetBasicAleForce(int bodyIndex, int iTimeLayer)
    {
        AleForceOfBody* aleForceOfBody = this->GetAleForceOfBody(bodyIndex);
        return aleForceOfBody->GetBasicAleForce(iTimeLayer);
    }

    void CopyBasicAleForce(BasicAleForce* basicAleForce, BasicAleForce* basicAleForceIn)
    {
        *basicAleForce = (*basicAleForceIn);
    }

    AerodynamicForce::AerodynamicForce()
    {
        componentName = "AllWall";
        this->Zero();
    }

    AerodynamicForce::~AerodynamicForce()
    {
    }

    void AerodynamicForce::Zero()
    {
        forceX = 0.0;
        forceY = 0.0;
        forceZ = 0.0;

        forceXByPressure = 0.0;
        forceYByPressure = 0.0;
        forceZByPressure = 0.0;

        forceXByViscosity = 0.0;
        forceYByViscosity = 0.0;
        forceZByViscosity = 0.0;

        momentX = 0.0;
        momentY = 0.0;
        momentZ = 0.0;

        momentXByPressure = 0.0;
        momentYByPressure = 0.0;
        momentZByPressure = 0.0;

        momentXByViscosity = 0.0;
        momentYByViscosity = 0.0;
        momentZByViscosity = 0.0;

        coefForceX = 0.0;
        coefForceY = 0.0;
        coefForceZ = 0.0;

        coefForceXByPressure = 0.0;
        coefForceYByPressure = 0.0;
        coefForceZByPressure = 0.0;

        coefForceXByViscosity = 0.0;
        coefForceYByViscosity = 0.0;
        coefForceZByViscosity = 0.0;

        coefMomentX = 0.0;
        coefMomentY = 0.0;
        coefMomentZ = 0.0;

        coefMomentXByPressure = 0.0;
        coefMomentYByPressure = 0.0;
        coefMomentZByPressure = 0.0;

        coefMomentXByViscosity = 0.0;
        coefMomentYByViscosity = 0.0;
        coefMomentZByViscosity = 0.0;
    }

    void AerodynamicForce::Initialize(RDouble angleOfAttackDenotedByRadian, RDouble angleOfSideslipDenotedByRadian, RDouble referenceAreaOfAerodynamicForce, RDouble referenceLengthOfAerodynamicMoment, RDouble referenceSpainLengthOfAerodynamicMoment)
    {
        this->angleOfAttackDenotedByRadian = angleOfAttackDenotedByRadian;
        this->angleOfSideslipDenotedByRadian = angleOfSideslipDenotedByRadian;
        this->referenceArea = referenceAreaOfAerodynamicForce;
        this->referenceLength = referenceLengthOfAerodynamicMoment;
        this->referenceSpainLength = referenceSpainLengthOfAerodynamicMoment;

        this->Zero();
    }

    void AerodynamicForce::SetTorqueRef(RDouble torqueRefX, RDouble torqueRefY, RDouble torqueRefZ)
    {
        this->torqueRefX = torqueRefX;
        this->torqueRefY = torqueRefY;
        this->torqueRefZ = torqueRefZ;
    };

    void AerodynamicForce::SetHingeMomentRef(RDouble *localCoordAxis0, RDouble *localCoordAxis1)
    {
        for (int m = 0; m < 3; ++m)
        {
            this->localCoordAxis0[m] = localCoordAxis0[m];
            this->localCoordAxis1[m] = localCoordAxis1[m];
        }
    };

    void AerodynamicForce::Add(DataContainer* dataContainer)
    {
        const int numberOfForceVariables = 18;
        RDouble forceMomentContainer[numberOfForceVariables];

        dataContainer->moveto_begin();

        PHSPACE::PHRead(dataContainer, forceMomentContainer, numberOfForceVariables);

        forceX += forceMomentContainer[0];
        forceY += forceMomentContainer[1];
        forceZ += forceMomentContainer[2];

        forceXByPressure += forceMomentContainer[3];
        forceYByPressure += forceMomentContainer[4];
        forceZByPressure += forceMomentContainer[5];

        forceXByViscosity += forceMomentContainer[6];
        forceYByViscosity += forceMomentContainer[7];
        forceZByViscosity += forceMomentContainer[8];

        momentX += forceMomentContainer[9];
        momentY += forceMomentContainer[10];
        momentZ += forceMomentContainer[11];

        momentXByPressure += forceMomentContainer[12];
        momentYByPressure += forceMomentContainer[13];
        momentZByPressure += forceMomentContainer[14];

        momentXByViscosity += forceMomentContainer[15];
        momentYByViscosity += forceMomentContainer[16];
        momentZByViscosity += forceMomentContainer[17];
    }

    void AerodynamicForce::GetMomentByOriginOfCoordinate(RDouble& aerodynamicMomentX0, RDouble& aerodynamicMomentY0, RDouble& aerodynamicMomentZ0)
    {
        RDouble momentReferencePoint[3];

        momentReferencePoint[0] = GlobalDataBase::GetDoubleParaFromDB("TorqueRefX");
        momentReferencePoint[1] = GlobalDataBase::GetDoubleParaFromDB("TorqueRefY");
        momentReferencePoint[2] = GlobalDataBase::GetDoubleParaFromDB("TorqueRefZ");

        RDouble aerodynamicForce[3];

        aerodynamicForce[0] = forceX;
        aerodynamicForce[1] = forceY;
        aerodynamicForce[2] = forceZ;

        RDouble changeOfAerodynamicMoment[3];
        PHSPACE::CrossProduct(momentReferencePoint, aerodynamicForce, changeOfAerodynamicMoment);

        aerodynamicMomentX0 = momentX + changeOfAerodynamicMoment[0];
        aerodynamicMomentY0 = momentY + changeOfAerodynamicMoment[1];
        aerodynamicMomentZ0 = momentZ + changeOfAerodynamicMoment[2];
    }

    void AerodynamicForce::DumpCoefficient(ActionKey* actionKey)
    {
        using namespace PHMPI;
        int myID = GetCurrentProcessorID();
        int server = GetServerProcessorID();
        if (myID != server) return;

        PHSPACE::ParallelOpenFile(actionKey);

        fstream& file = *actionKey->GetFile();

        this->DumpCoefficient(file);

        PHSPACE::ParallelCloseFile(file);
    }

    void AerodynamicForce::DumpComponentCoefficient(ActionKey* actionKey, int iComponent)
    {
        using namespace PHMPI;
        int myID = GetCurrentProcessorID();
        int server = GetServerProcessorID();
        if (myID != server) return;

        int nTotalBC = GlobalBoundaryCondition::GetNumberofBC();
        for (int iPart = 0; iPart < nTotalBC; ++iPart)
        {
            SimpleBC* boundaryCondition = GlobalBoundaryCondition::GetBC(iPart);
            string partName = boundaryCondition->GetBCName();
            if (FaceInAerodynamicForceComponents(partName, iComponent))
            {
                componentName = partName;
            }
        }

        string totalForceFileName = actionKey->GetFileName();
        string partAircoefFileName = PHSPACE::AddSymbolToFileName(totalForceFileName, "_", componentName);

        fstream file;
        file.open(partAircoefFileName.c_str(), ios_base::out | ios_base::app);

        this->DumpCoefficient(file);

        file.close();
        file.clear();
    }

    void AerodynamicForce::DumpCoefficient(fstream& file)
    {
        this->ComputeCoefficient();

        ostringstream oss;

        int outerIterationSteps = GlobalDataBase::GetIntParaFromDB("outnstep");

        oss << setiosflags(ios::right);
        oss << setiosflags(ios::scientific);

        if (PHSPACE::IfFileEmpty(file))
        {
            vector<string> titleOfTecplot;
            titleOfTecplot.push_back("Title=\"THE AIRCOEF\"");
            titleOfTecplot.push_back("Variables=");
            titleOfTecplot.push_back("iter");
            titleOfTecplot.push_back("<i>C<sub>L</sub></i>");
            titleOfTecplot.push_back("<i>C<sub>D</sub></i>");
            titleOfTecplot.push_back("<i>C<sub>Z</sub></i>");
            titleOfTecplot.push_back("<i>C<sub>D_p</sub></i>");
            titleOfTecplot.push_back("<i>C<sub>D_f</sub></i>");
            titleOfTecplot.push_back("<i>C<sub>D</sub></i>-<i>C<sub>L</sub></i><sup>2</sup>/(<greek>p</greek><math>4</math>A<sub>r</sub>)</i>");
            titleOfTecplot.push_back("<i>X<sub>cp</sub></i>");
            titleOfTecplot.push_back("<i>C<sub>A</sub></i>");
            titleOfTecplot.push_back("<i>C<sub>N</sub></i>");
            titleOfTecplot.push_back("<i>C<sub>Z1</sub></i>");
            titleOfTecplot.push_back("<i>C<sub>ml</sub></i>");
            titleOfTecplot.push_back("<i>C<sub>mn</sub></i>");
            titleOfTecplot.push_back("<i>C<sub>m</sub></i>");
            if (componentName == "AllWall")
            {
                titleOfTecplot.push_back("\"WallTime\"");
            }

            for (std::size_t iTitle = 0; iTitle < titleOfTecplot.size(); ++iTitle)
            {
                oss << titleOfTecplot[iTitle] << "\n";
            }
        }

        oss << outerIterationSteps << "\t";

        oss << setprecision(10);

        oss << totalLiftCoefficient << "\t";
        oss << totalDragCoefficient << "\t";
        oss << totalCrossCoefficient << "\t";
        oss << totalPressureDragCoefficient << "\t";
        oss << totalViscosityDragCoefficient << "\t";
        oss << DragCoefLiftCoef2piar << "\t";
        oss << pressureCenter << "\t";

        oss << coefForceX << "\t";
        oss << coefForceY << "\t";
        oss << coefForceZ << "\t";
        oss << coefMomentX << "\t";
        oss << coefMomentY << "\t";
        oss << coefMomentZ << "\t";

       if (componentName == "AllWall")
       {
            RDouble wallTime = PHSPACE::TIME_SPACE::GetWallTime();
            oss << wallTime;
       }

        oss << "\n";
        PHSPACE::WriteASCIIFile(file, oss.str());
    }

    void AerodynamicForce::ComputeCoefficient()
    {
        RDouble sina = sin(angleOfAttackDenotedByRadian);
        RDouble cosa = cos(angleOfAttackDenotedByRadian);
        RDouble sinb = sin(angleOfSideslipDenotedByRadian);
        RDouble cosb = cos(angleOfSideslipDenotedByRadian);

        const RDouble AR = 9.5;

        RDouble referenceCoefficientOfForce = 1.0 / referenceArea;
        RDouble referenceCoefficientOfMoment = 1.0 / (referenceArea * referenceLength);
        RDouble referenceCoefficientOfMomentSpain = 1.0 / (referenceArea * referenceSpainLength);

        coefForceX = forceX * referenceCoefficientOfForce;
        coefForceY = forceY * referenceCoefficientOfForce;
        coefForceZ = forceZ * referenceCoefficientOfForce;

        coefForceXByPressure = forceXByPressure * referenceCoefficientOfForce;
        coefForceYByPressure = forceYByPressure * referenceCoefficientOfForce;
        coefForceZByPressure = forceZByPressure * referenceCoefficientOfForce;

        coefForceXByViscosity = forceXByViscosity * referenceCoefficientOfForce;
        coefForceYByViscosity = forceYByViscosity * referenceCoefficientOfForce;
        coefForceZByViscosity = forceZByViscosity * referenceCoefficientOfForce;

        coefMomentX = momentX * referenceCoefficientOfMomentSpain;
        coefMomentY = momentY * referenceCoefficientOfMomentSpain;
        coefMomentZ = momentZ * referenceCoefficientOfMoment;

        coefMomentXByPressure = momentXByPressure * referenceCoefficientOfMomentSpain;
        coefMomentYByPressure = momentYByPressure * referenceCoefficientOfMomentSpain;
        coefMomentZByPressure = momentZByPressure * referenceCoefficientOfMoment;

        coefMomentXByViscosity = momentXByViscosity * referenceCoefficientOfMomentSpain;
        coefMomentYByViscosity = momentYByViscosity * referenceCoefficientOfMomentSpain;
        coefMomentZByViscosity = momentZByViscosity * referenceCoefficientOfMoment;

        totalLiftCoefficient = -sina * coefForceX + cosa * coefForceY;
        totalDragCoefficient = cosa * cosb * coefForceX + sina * cosb * coefForceY + sinb * coefForceZ;
        totalCrossCoefficient = -cosa * sinb * coefForceX - sina * sinb * coefForceY + cosb * coefForceZ;
        totalPressureDragCoefficient = cosa * cosb * coefForceXByPressure + sina * cosb * coefForceYByPressure + sinb * coefForceZByPressure;
        totalViscosityDragCoefficient = cosa * cosb * coefForceXByViscosity + sina * cosb * coefForceYByViscosity + sinb * coefForceZByViscosity;
        DragCoefLiftCoef2piar = totalDragCoefficient - totalLiftCoefficient * totalLiftCoefficient / (PI * AR);
        pressureCenter = PHSPACE::SIGN(1.0, coefForceY) * coefMomentZ / (PHSPACE::ABS(coefForceY) + SMALL);
    }

    void AerodynamicForce::CollectOfAllProcessors()
    {
        using namespace PHMPI;
        RDouble forceXSum = forceX;
        RDouble forceYSum = forceY;
        RDouble forceZSum = forceZ;

        RDouble forceXByPressureSum = forceXByPressure;
        RDouble forceYByPressureSum = forceYByPressure;
        RDouble forceZByPressureSum = forceZByPressure;

        RDouble forceXByViscositySum = forceXByViscosity;
        RDouble forceYByViscositySum = forceYByViscosity;
        RDouble forceZByViscositySum = forceZByViscosity;

        RDouble momentXSum = momentX;
        RDouble momentYSum = momentY;
        RDouble momentZSum = momentZ;

        RDouble momentXByPressureSum = momentXByPressure;
        RDouble momentYByPressureSum = momentYByPressure;
        RDouble momentZByPressureSum = momentZByPressure;

        RDouble momentXByViscositySum = momentXByViscosity;
        RDouble momentYByViscositySum = momentYByViscosity;
        RDouble momentZByViscositySum = momentZByViscosity;

        PH_Reduce(&forceX, &forceXSum, 1, PH_SUM);
        PH_Reduce(&forceY, &forceYSum, 1, PH_SUM);
        PH_Reduce(&forceZ, &forceZSum, 1, PH_SUM);

        PH_Reduce(&forceXByPressure, &forceXByPressureSum, 1, PH_SUM);
        PH_Reduce(&forceYByPressure, &forceYByPressureSum, 1, PH_SUM);
        PH_Reduce(&forceZByPressure, &forceZByPressureSum, 1, PH_SUM);

        PH_Reduce(&forceXByViscosity, &forceXByViscositySum, 1, PH_SUM);
        PH_Reduce(&forceYByViscosity, &forceYByViscositySum, 1, PH_SUM);
        PH_Reduce(&forceZByViscosity, &forceZByViscositySum, 1, PH_SUM);

        PH_Reduce(&momentX, &momentXSum, 1, PH_SUM);
        PH_Reduce(&momentY, &momentYSum, 1, PH_SUM);
        PH_Reduce(&momentZ, &momentZSum, 1, PH_SUM);

        PH_Reduce(&momentXByPressure, &momentXByPressureSum, 1, PH_SUM);
        PH_Reduce(&momentYByPressure, &momentYByPressureSum, 1, PH_SUM);
        PH_Reduce(&momentZByPressure, &momentZByPressureSum, 1, PH_SUM);

        PH_Reduce(&momentXByViscosity, &momentXByViscositySum, 1, PH_SUM);
        PH_Reduce(&momentYByViscosity, &momentYByViscositySum, 1, PH_SUM);
        PH_Reduce(&momentZByViscosity, &momentZByViscositySum, 1, PH_SUM);

        forceX = forceXSum;
        forceY = forceYSum;
        forceZ = forceZSum;

        forceXByPressure = forceXByPressureSum;
        forceYByPressure = forceYByPressureSum;
        forceZByPressure = forceZByPressureSum;

        forceXByViscosity = forceXByViscositySum;
        forceYByViscosity = forceYByViscositySum;
        forceZByViscosity = forceZByViscositySum;

        momentX = momentXSum;
        momentY = momentYSum;
        momentZ = momentZSum;

        momentXByPressure = momentXByPressureSum;
        momentYByPressure = momentYByPressureSum;
        momentZByPressure = momentZByPressureSum;

        momentXByViscosity = momentXByViscositySum;
        momentYByViscosity = momentYByViscositySum;
        momentZByViscosity = momentZByViscositySum;
    }

    AerodynamicForceManager::AerodynamicForceManager()
    {
        isAllocateMemory  = false;
        isComponentsExist = false;
        numberOfAerodynamicForceComponents = 0;
        partBCID = NULL;
    }

    AerodynamicForceManager::~AerodynamicForceManager()
    {
        this->DeallocateMemory();
    }

    void AerodynamicForceManager::AllocateMemory()
    {
        if (isAllocateMemory) return;

        isAllocateMemory = true;

        aerodynamicForceContainer.resize(numberOfAerodynamicForceComponents+1);

        for (int iAerodynamicForces = 0; iAerodynamicForces <= numberOfAerodynamicForceComponents; ++iAerodynamicForces)
        {
            aerodynamicForceContainer[iAerodynamicForces] = new AerodynamicForce();
        }
    }

    void AerodynamicForceManager::DeallocateMemory()
    {
        if (!isAllocateMemory) return;

        for (int iAerodynamicForces = 0; iAerodynamicForces <= numberOfAerodynamicForceComponents; ++iAerodynamicForces)
        {
            delete aerodynamicForceContainer[iAerodynamicForces];
        }
    }

    void AerodynamicForceManager::InitializeAerodynamicForceManager()
    {
        isComponentsExist = (GlobalBoundaryCondition::GetGlobalBoundaryConditionList() != 0);
        numberOfAerodynamicForceComponents = 0;
        if(isComponentsExist)
        {
            numberOfAerodynamicForceComponents = GlobalBoundaryCondition::GetNumberofSolidWallPart();
        }
        GlobalDataBase::UpdateData("numberOfAerodynamicForceComponents", &numberOfAerodynamicForceComponents, PHINT, 1);

        this->AllocateMemory();

        if (!partBCID)
        {
            partBCID = new int[numberOfAerodynamicForceComponents+1];
        }
        const int ALL_SOLID_WALL = -1;
        PHSPACE::SetField(partBCID, ALL_SOLID_WALL, numberOfAerodynamicForceComponents+1);

        RDouble forceReferenceArea, forceReferenceLength, forceReferenceLengthSpanWise;
        RDouble TorqueRefX, TorqueRefY, TorqueRefZ;
        RDouble angleOfAttackDenotedByRadian = GlobalDataBase::GetDoubleParaFromDB("attack");
        RDouble angleOfSideslipDenotedByRadian = GlobalDataBase::GetDoubleParaFromDB("sideslip");

        int partCount = 0;
        uint_t nTBC = GlobalBoundaryCondition::GetNumberofBC();
        for (int iBC = 0; iBC < nTBC; ++iBC)
        {
            if (!GlobalBoundaryCondition::IsSolidWallBC(iBC)) continue;

            partBCID[partCount] = iBC;

            SimpleBC *bc = GlobalBoundaryCondition::GetBC(iBC);
            Data_Param *bcParamDB = bc->GetBCParamDataBase();

            AerodynamicForce* aerodynamicForce = this->GetAerodynamicForce(partCount);

            bcParamDB->GetData("forceReferenceArea",           &forceReferenceArea,           PHDOUBLE, 1);
            bcParamDB->GetData("forceReferenceLength",         &forceReferenceLength,         PHDOUBLE, 1);
            bcParamDB->GetData("forceReferenceLengthSpanWise", &forceReferenceLengthSpanWise, PHDOUBLE, 1);
            aerodynamicForce->Initialize(angleOfAttackDenotedByRadian, angleOfSideslipDenotedByRadian, forceReferenceArea,
                forceReferenceLength, forceReferenceLengthSpanWise);

            bcParamDB->GetData("TorqueRefX", &TorqueRefX, PHDOUBLE, 1);
            bcParamDB->GetData("TorqueRefY", &TorqueRefY, PHDOUBLE, 1);
            bcParamDB->GetData("TorqueRefZ", &TorqueRefZ, PHDOUBLE, 1);
            aerodynamicForce->SetTorqueRef(TorqueRefX, TorqueRefY, TorqueRefZ);
            partCount ++;
        }

        forceReferenceArea = GlobalDataBase::GetDoubleParaFromDB("forceReferenceArea");
        forceReferenceLength = GlobalDataBase::GetDoubleParaFromDB("forceReferenceLength");
        forceReferenceLengthSpanWise = GlobalDataBase::GetDoubleParaFromDB("forceReferenceLengthSpanWise");
        TorqueRefX = GlobalDataBase::GetDoubleParaFromDB("TorqueRefX");
        TorqueRefY = GlobalDataBase::GetDoubleParaFromDB("TorqueRefY");
        TorqueRefZ = GlobalDataBase::GetDoubleParaFromDB("TorqueRefZ");

        //! Total aerodynamic force reference
        int iAerodynamicForces = numberOfAerodynamicForceComponents;
        AerodynamicForce* aerodynamicForce = this->GetAerodynamicForce(iAerodynamicForces);
        aerodynamicForce->Initialize(angleOfAttackDenotedByRadian, angleOfSideslipDenotedByRadian, forceReferenceArea,
            forceReferenceLength, forceReferenceLengthSpanWise);
        aerodynamicForce->SetTorqueRef(TorqueRefX, TorqueRefY, TorqueRefZ);
    }

    void AerodynamicForceManager::ExtractForceDataFromActionKey(ActionKey* actionKey, vector< DataContainer* >* dataContainerArray)
    {
        DataContainer* dataContainer = actionKey->GetData();
        dataContainer->moveto_begin();

        const int numberOfData = 18;
        RDouble bufferForSwap[numberOfData];

        for (int iAerodynamicForces = 0; iAerodynamicForces <= numberOfAerodynamicForceComponents; ++iAerodynamicForces)
        {
            PHSPACE::PHRead(dataContainer, bufferForSwap, numberOfData);
            PHSPACE::PHWrite((*dataContainerArray)[iAerodynamicForces], bufferForSwap, numberOfData);
        }
    }

    vector< DataContainer* >* AerodynamicForceManager::CreateDataContainerArray()
    {
        vector< DataContainer* >* dataContainerArray = new vector< DataContainer* >(numberOfAerodynamicForceComponents+1);
        for (int iAerodynamicForces = 0; iAerodynamicForces <= numberOfAerodynamicForceComponents; ++iAerodynamicForces)
        {
            (*dataContainerArray)[iAerodynamicForces] = new DataContainer();
        }
        return dataContainerArray;
    }

    void AerodynamicForceManager::FreeDataContainerArray(vector< DataContainer* >* dataContainerArray)
    {
        for (int iAerodynamicForces = 0; iAerodynamicForces <= numberOfAerodynamicForceComponents; ++iAerodynamicForces)
        {
            delete (*dataContainerArray)[iAerodynamicForces];
        }
        delete dataContainerArray;
    }

    void AerodynamicForceManager::AddAerodynamicForce(ActionKey* actionKey)
    {
        vector< DataContainer* >* dataContainerArray = this->CreateDataContainerArray();

        this->ExtractForceDataFromActionKey(actionKey, dataContainerArray);

        for (int iAerodynamicForces = 0; iAerodynamicForces < numberOfAerodynamicForceComponents; ++iAerodynamicForces)
        {
            DataContainer* dataContainer = (*dataContainerArray)[iAerodynamicForces];
            aerodynamicForceContainer[iAerodynamicForces]->Add(dataContainer);
            aerodynamicForceContainer[numberOfAerodynamicForceComponents]->Add(dataContainer);
        }
        this->FreeDataContainerArray(dataContainerArray);
    }

    void AerodynamicForceManager::DumpAerodynamicForceCoefficient(ActionKey* actionKey)
    {
        aerodynamicForceContainer[numberOfAerodynamicForceComponents]->DumpCoefficient(actionKey);

        for (int iAerodynamicForces = 0; iAerodynamicForces < numberOfAerodynamicForceComponents; ++iAerodynamicForces)
        {
            aerodynamicForceContainer[iAerodynamicForces]->DumpComponentCoefficient(actionKey, iAerodynamicForces);
        }

        DumpAllForceData(actionKey);
    }

    void AerodynamicForceManager::DumpAllForceData(ActionKey *actkey)
    {
        string fileName = actkey->filename;

        string mainName, fileDir, extensionName;
        string globalFileName, localFileName;
        GetFileNameExtension(fileName, mainName, extensionName, ".");
        GetFileNameExtension(mainName, fileDir, extensionName, "/");

        ostringstream ossGlobalFileName, ossLocalFileName, ossGlobalData, ossLocalData;
        ossGlobalFileName << fileDir <<"/force_global.dat";
        ossLocalFileName << fileDir <<"/force_local.dat";
        globalFileName = ossGlobalFileName.str();
        localFileName = ossLocalFileName.str();

        string gridFile = GlobalDataBase::GetStrParaFromDB("gridfile");
        RDouble refMachNumber = GlobalDataBase::GetDoubleParaFromDB("refMachNumber");
        RDouble attackd = GlobalDataBase::GetDoubleParaFromDB("attackd");
        RDouble angleSlide = GlobalDataBase::GetDoubleParaFromDB("angleSlide");
        RDouble refReNumber = GlobalDataBase::GetDoubleParaFromDB("refReNumber");
        RDouble referenceArea = GlobalDataBase::GetDoubleParaFromDB("forceReferenceArea");
        RDouble referenceLength = GlobalDataBase::GetDoubleParaFromDB("forceReferenceLength");
        RDouble referenceSpainLength = GlobalDataBase::GetDoubleParaFromDB("forceReferenceLengthSpanWise");
        RDouble TorqueRefX = GlobalDataBase::GetDoubleParaFromDB("TorqueRefX");
        RDouble TorqueRefY = GlobalDataBase::GetDoubleParaFromDB("TorqueRefY");
        RDouble TorqueRefZ = GlobalDataBase::GetDoubleParaFromDB("TorqueRefZ");
        RDouble gridScaleFactor = GlobalDataBase::GetDoubleParaFromDB("gridScaleFactor");
        RDouble wallTemperature = GlobalDataBase::GetDoubleParaFromDB("wallTemperature");
        RDouble reynoldsReferenceLengthDimensional = GlobalDataBase::GetDoubleParaFromDB( "reynoldsReferenceLengthDimensional" );
        int outerIterationSteps = GlobalDataBase::GetIntParaFromDB("outnstep");

        RDouble height = 0.0;
        int inflowParaType = GlobalDataBase::GetIntParaFromDB("inflowParaType");
        if (inflowParaType == 1)
        {
            height = GlobalDataBase::GetDoubleParaFromDB("height");
        }

        //! Write global force data file title
        ossGlobalData << "#     MeshFile                     :       " << gridFile << endl;
        if (inflowParaType == 1)
        {
            ossGlobalData << "#     Height                       :       " << height << "  km" << endl;
        }

        ossGlobalData << setiosflags(ios::scientific);
        ossGlobalData << setprecision(6);
        ossGlobalData << "#     Reynolds                     :       " << refReNumber / reynoldsReferenceLengthDimensional << "  1/m" << endl;
        ossGlobalData << resetiosflags(ios::scientific);

        ossGlobalData << "#     Mach                         :       " << refMachNumber << endl;
        ossGlobalData << "#     Alpha                        :       " << attackd << "  degree" << endl;
        ossGlobalData << "#     Beta                         :       " << angleSlide << "  degree" << endl;

        if (wallTemperature < 0.0)
        {
            ossGlobalData << "#     Wall Temperature             :       " << "Adiabatic Wall" << endl;
        }
        else
        {
            ossGlobalData << "#     Wall Temperature             :       " << wallTemperature << "  K" << endl;
        }

        ossGlobalData << "#     gridScaleFactor              :       " << gridScaleFactor << endl;
        ossGlobalData << "#     forceReferenceLength         :       " << referenceLength * reynoldsReferenceLengthDimensional << "  m" << endl;
        ossGlobalData << "#     forceReferenceLengthSpanWise :       " << referenceSpainLength * reynoldsReferenceLengthDimensional << "  m" << endl;

        if (GetDim() == THREE_D)
        {
            ossGlobalData << "#     forceReferenceArea           :       " << referenceArea * reynoldsReferenceLengthDimensional * reynoldsReferenceLengthDimensional << "  m2" << endl;
        }
        else
        {
            ossGlobalData << "#     forceReferenceArea           :       " << referenceArea * reynoldsReferenceLengthDimensional << "  m^2" << endl;
        }

        ossGlobalData << "#     TorqueRefX                   :       " << TorqueRefX * reynoldsReferenceLengthDimensional << endl;
        ossGlobalData << "#     TorqueRefY                   :       " << TorqueRefY * reynoldsReferenceLengthDimensional << endl;
        ossGlobalData << "#     TorqueRefZ                   :       " << TorqueRefZ * reynoldsReferenceLengthDimensional << endl;
        ossGlobalData << "#     iter                         :       " << outerIterationSteps << endl << endl;
        ossGlobalData << "#        CA             CN            CZ1           Cml             Cmn            Cm"   << endl;
        ossGlobalData << "#        CA_p           CN_p          CZ1_p         Cml_p           Cmn_p          Cm_p" << endl;
        ossGlobalData << "#        CA_f           CN_f          CZ1_f         Cml_f           Cmn_f          Cm_f" << endl;

        //! Write local force data file title
        ossLocalData << "#     MeshFile             :       " << gridFile << endl;
        if (inflowParaType == 1)
        {
            ossLocalData << "#     Height               :       " << height << "  km" << endl;
        }

        ossLocalData << setiosflags(ios::scientific);
        ossLocalData << setprecision(6);
        ossLocalData << "#     Reynolds             :       " << refReNumber / reynoldsReferenceLengthDimensional << "  1/m" << endl;
        ossLocalData << resetiosflags(ios::scientific);

        ossLocalData << "#     Mach                 :       " << refMachNumber << endl;
        ossLocalData << "#     Alpha                :       " << attackd << "  degree" << endl;
        ossLocalData << "#     Beta                 :       " << angleSlide << "  degree" << endl;

        if (wallTemperature < 0.0)
        {
            ossLocalData << "#     Wall Temperature     :       " << "Adiabatic Wall" << endl;
        }
        else
        {
            ossLocalData << "#     Wall Temperature     :       " << wallTemperature << "  K" << endl;
        }

        ossLocalData << "#     gridScaleFactor      :       " << gridScaleFactor << endl;
        ossLocalData << "#     iter                 :       " << outerIterationSteps << endl << endl;
        ossLocalData << "#        CA             CN            CZ1           Cml             Cmn            Cm"   << endl;
        ossLocalData << "#        CA_p           CN_p          CZ1_p         Cml_p           Cmn_p          Cm_p" << endl;
        ossLocalData << "#        CA_f           CN_f          CZ1_f         Cml_f           Cmn_f          Cm_f" << endl;
        ossLocalData << "#        CMh" << endl;

        string partName = "whole";;
        aerodynamicForceContainer[numberOfAerodynamicForceComponents]->DumpAllGlobalForceData(partName, numberOfAerodynamicForceComponents, ossGlobalData);

        //! Write force data file
        for (int iComponent = 0; iComponent < numberOfAerodynamicForceComponents; ++iComponent)
        {
            int bcID = partBCID[iComponent];

            SimpleBC* bc = GlobalBoundaryCondition::GetBC(bcID);
            const string& bcName = bc->GetBCName();
            partName = bcName;

            aerodynamicForceContainer[iComponent]->DumpAllGlobalForceData(partName, iComponent, ossGlobalData);
            aerodynamicForceContainer[iComponent]->DumpAllLocalForceData(partName, iComponent, ossLocalData);

            //! Write global force data file
            fstream globalFile;
            ParallelOpenFile(globalFile, globalFileName, ios_base::out | ios_base::trunc);
            WriteASCIIFile(globalFile, ossGlobalData.str());
            ParallelCloseFile(globalFile);

            //! Write local force data file
            fstream localFile;
            ParallelOpenFile(localFile, localFileName, ios_base::out | ios_base::trunc);
            WriteASCIIFile(localFile, ossLocalData.str());
            ParallelCloseFile(localFile);
        }
    }

    void AerodynamicForce::DumpAllLocalForceData(const string &partName, int partID, std::ostringstream &oss)
    {
        RDouble reynoldsReferenceLengthDimensional = GlobalDataBase::GetDoubleParaFromDB("reynoldsReferenceLengthDimensional");

        oss << resetiosflags(ios::scientific);

        oss <<endl;
        oss << "*"<<partName <<endl;
        oss << "#forceReferenceLength         :    " << referenceLength * reynoldsReferenceLengthDimensional << "  m" << endl;
        oss << "#forceReferenceLengthSpanWise :    " << referenceSpainLength * reynoldsReferenceLengthDimensional << "  m" <<endl;

        if (GetDim() == THREE_D)
        {
            oss << "#forceReferenceArea           :    " << referenceArea * reynoldsReferenceLengthDimensional * reynoldsReferenceLengthDimensional << "  m2" <<endl;
        }
        else
        {
            oss << "#forceReferenceArea           :    " << referenceArea * reynoldsReferenceLengthDimensional << "  m^2" <<endl;
        }

        oss << "#TorqueRefX                   :    " << torqueRefX * reynoldsReferenceLengthDimensional << endl;
        oss << "#TorqueRefY                   :    " << torqueRefY * reynoldsReferenceLengthDimensional << endl;
        oss << "#TorqueRefZ                   :    " << torqueRefZ * reynoldsReferenceLengthDimensional << endl;

        oss << setiosflags(ios::right);
        oss << setiosflags(ios::scientific);
        oss << setprecision(10);

        oss << coefForceX << "	";
        oss << coefForceY << "	";
        oss << coefForceZ << "	";
        oss << coefMomentX << "	";
        oss << coefMomentY << "	";
        oss << coefMomentZ << endl;

        oss << coefForceXByPressure << "	";
        oss << coefForceYByPressure << "	";
        oss << coefForceZByPressure << "	";
        oss << coefMomentXByPressure << "	";
        oss << coefMomentYByPressure << "	";
        oss << coefMomentZByPressure << endl;

        oss << coefForceXByViscosity << "	";
        oss << coefForceYByViscosity << "	";
        oss << coefForceZByViscosity << "	";
        oss << coefMomentXByViscosity << "	";
        oss << coefMomentYByViscosity << "	";
        oss << coefMomentZByViscosity << endl;
    }

    void AerodynamicForce::DumpAllGlobalForceData(const string &partName, int partID, std::ostringstream &oss)
    {
        oss << setiosflags(ios::right);
        oss << setiosflags(ios::scientific);
        oss << setprecision(10);

        oss <<endl;
        oss << "*"<<partName <<endl;

        oss << coefForceX << "	";
        oss << coefForceY << "	";
        oss << coefForceZ << "	";
        oss << coefMomentX << "	";
        oss << coefMomentY << "	";
        oss << coefMomentZ << endl;

        oss << coefForceXByPressure << "	";
        oss << coefForceYByPressure << "	";
        oss << coefForceZByPressure << "	";
        oss << coefMomentXByPressure << "	";
        oss << coefMomentYByPressure << "	";
        oss << coefMomentZByPressure << endl;

        oss << coefForceXByViscosity << "	";
        oss << coefForceYByViscosity << "	";
        oss << coefForceZByViscosity << "	";
        oss << coefMomentXByViscosity << "	";
        oss << coefMomentYByViscosity << "	";
        oss << coefMomentZByViscosity << endl;

    }

    void AerodynamicForceManager::ComputeAerodynamicForceCoefficient()
    {
        aerodynamicForceContainer[numberOfAerodynamicForceComponents]->ComputeCoefficient();
    }

    void AerodynamicForceManager::CollectAerodynamicForceOfAllProcessors()
    {
        for (int iAerodynamicForces = 0; iAerodynamicForces <= numberOfAerodynamicForceComponents; ++iAerodynamicForces)
        {
            aerodynamicForceContainer[iAerodynamicForces]->CollectOfAllProcessors();
        }
    }

    AerodynamicForceManager* GetAerodynamicForceManager()
    {
        return &aerodynamicForceManager;
    }

    AerodynamicForce* GetCurrentAerodynamicForce(int bodyIndex)
    {
        return GetAerodynamicForceManager()->GetAerodynamicForce(bodyIndex);
    }

    BasicAerodynamicForce* GetBasicAerodynamicForce(int iBoundaryConditionRegion, BasicAerodynamicForce** basicAerodynamicForceContainer)
    {
        return basicAerodynamicForceContainer[iBoundaryConditionRegion];
    }

    BasicAerodynamicForce* GetBasicAerodynamicForce(int iBoundaryConditionRegion, vector< BasicAerodynamicForce* >* basicAerodynamicForceVector)
    {
        return (*basicAerodynamicForceVector)[iBoundaryConditionRegion];
    }

    bool OutPutSubIterationForceHeader()
    {
        int outerIterationSteps = GlobalDataBase::GetIntParaFromDB("outnstep");
        int innerIterationSteps = GlobalDataBase::GetIntParaFromDB("innstep");
        int numberOfForceSaveSteps = GlobalDataBase::GetIntParaFromDB("intervalStepForce");

        if (outerIterationSteps / numberOfForceSaveSteps == 1 && innerIterationSteps / numberOfForceSaveSteps == 1)
        {
            return true;
        }

        return false;
    }

}

