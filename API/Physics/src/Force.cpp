#include "Force.h"
#include "Math_BasisFunction.h"
#include "Constants.h"
#include "PHIO.h"
#include "IO_FileName.h"
#include "BC.h"
#include "TK_Time.h"
#include "Glb_Dimension.h"
#include "GlobalDataBase.h"
#include "PHHeader.h"

using namespace std;

namespace PHSPACE
{

Force globalAerodynamicForce;

void InitForce (ActionKey *actkey)
{
    int myid   = PHMPI::GetCurrentProcessorID();
    int server = PHMPI::GetServerProcessorID();
    //if ( myid != server ) return;

    RDouble attack,sideslip;
    RDouble attackd = GlobalDataBase::GetDoubleParaFromDB("attackd");
    RDouble angleSlide = GlobalDataBase::GetDoubleParaFromDB("angleSlide");

    attack   = attackd * PI / 180.0;
    sideslip = angleSlide * PI / 180.0;

    RDouble forceRefenenceArea   = GlobalDataBase::GetDoubleParaFromDB("forceRefenenceArea");
    RDouble forceRefenenceLength = GlobalDataBase::GetDoubleParaFromDB("forceRefenenceLength");
    RDouble forceRefenenceLengthSpanWise = GlobalDataBase::GetDoubleParaFromDB("forceRefenenceLengthSpanWise");

    globalAerodynamicForce.Init(attack, sideslip, forceRefenenceArea, forceRefenenceLength, forceRefenenceLengthSpanWise);
}


void CollectionForce(ActionKey *actkey)
{
    globalAerodynamicForce.CollectionForce(actkey->GetData());
}

void DumpForce(ActionKey *actkey)
{
	int iterationStep;

	//int nMGLevel = GlobalDataBase::GetIntParaFromDB("nMGLevel");
	////! Steady
	//if (nMGLevel <= 1)
	//{
	//	//! Finest grid.
		GlobalDataBase::GetData("outnstep",&iterationStep, PHINT, 1);
	//}
	//else
	//{
	//	//! Coarse grid.
	//	GlobalDataBase::GetData("newnstep", &iterationStep, PHINT, 1);
	//}

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
}

Force::~Force()
{
    delete []partBCID;

    delete []cpx;
    delete []cpy;
    delete []cpz;

    delete []cl_tot  ;
    delete []cd_tot  ;
    delete []cd_pr   ;
    delete []cd_sf   ;
    delete []cd_cl2pa;
    delete []xcp     ;
    delete []ycp     ;
    delete []zcp     ;
    delete []side    ;

    delete []CA_f;
    delete []CA_p;
    delete []CN_f;
    delete []CN_p;
    delete []CZ_f;
    delete []CZ_p;
    delete []Cl_f;
    delete []Cl_p;
    delete []Cn_f;
    delete []Cn_p;
    delete []Cm_f;
    delete []Cm_p;
}

void Force::Init(RDouble attack, RDouble sideslip, RDouble forceRefenenceArea, RDouble forceRefenenceLength, RDouble forceRefenenceLengthSpanWise)
{
    this->attack   = attack;
    this->sideslip = sideslip;
    this->forceRefenenceArea   = forceRefenenceArea;
    this->forceRefenenceLength = forceRefenenceLength;
    this->forceRefenenceLengthSpanWise = forceRefenenceLengthSpanWise;

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

    if(!partBCID)
    {
        partBCID = new int [numberofTotalPart];
    }
    const int ALL_SOLID_WALL = -1;
    PHSPACE::SetField(partBCID, ALL_SOLID_WALL, numberofTotalPart);

    int partCount = 1;
    uint_t nTBC = GlobalBoundaryCondition::GetNumberofBC();
    for (int iBC = 0; iBC < nTBC; ++ iBC)
    {
        if(!GlobalBoundaryCondition::IsSolidWallBC(iBC)) continue;

        //! Because of global and local part force, partBCID set twice
        partBCID[partCount ++] = iBC;
        partBCID[partCount ++] = iBC;
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
                ifRead = GlobalBoundaryCondition::IsSolidWallBC(iPart);
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
        PH_AllreduceSepMode(localForceVar, globalForceVar, nTotalVarDim, MPI_DOUBLE ,PH_SUM);
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
        RDouble scf,vcm, Sb, forceRefenenceAreaPart, forceRefenenceLengthPart, forceRefenenceLengthPart_B;
        int bcID = partBCID[iPart];
        if((bcID == ALL_SOLID_WALL) || (isGlobalCoordinate == true))
        {
            forceRefenenceAreaPart   = this->forceRefenenceArea;
            forceRefenenceLengthPart = this->forceRefenenceLength;
            forceRefenenceLengthPart_B = this->forceRefenenceLengthSpanWise;
            isGlobalCoordinate = false;
        }
        else
        {
            BoundaryCondition *bc = GlobalBoundaryCondition::GetBC(bcID);
            Data_Param *bcParamDB = bc->GetBCParamDataBase();
            bcParamDB->GetData("forceRefenenceArea",   &forceRefenenceAreaPart,   PHDOUBLE, 1);
            bcParamDB->GetData("forceRefenenceLength", &forceRefenenceLengthPart, PHDOUBLE, 1);
            bcParamDB->GetData("forceRefenenceLengthSpanWise", &forceRefenenceLengthPart_B, PHDOUBLE, 1);
            isGlobalCoordinate = true;
        }

        scf = 1.0 / forceRefenenceAreaPart;
        vcm = 1.0 / (forceRefenenceAreaPart * forceRefenenceLengthPart);
        Sb = 1.0 / (forceRefenenceAreaPart * forceRefenenceLengthPart_B);

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

}

//! Bell modify 20120605.
void Force::DumpForce(int outnstep,std::ostringstream &oss)
{
    int intervalStepForce = GlobalDataBase::GetIntParaFromDB("intervalStepForce");

    //string header = "# iter		CL		Cd		Side		Cd_p		Cd_f\
    //Cd-CL2par	Xcp		cfx		cfy		cfz		Cmx		Cmy		Cmz		Walltime";

    int wordwidth = 16;
    oss << setiosflags(ios::left);
    oss << setiosflags(ios::scientific);
    oss << setprecision(10);

    //if (outnstep / intervalStepForce == 1)
    //{
    //	oss << header << endl;
    //}
    //int nMGLevel = GlobalDataBase::GetIntParaFromDB("nMGLevel");
    //int mgInitStep = PHSPACE::GlobalDataBase::GetIntParaFromDB( "mgInitStep" );
    //int totalMGInitStep = 0;
    //if (nMGLevel > 1)
    //{
    //	for (int iLevel = 1; iLevel < nMGLevel; iLevel ++)
    //	{
    //		totalMGInitStep += iLevel * mgInitStep;
    //	}
    //}
    const int ALL = 0;

    //oss << outnstep + totalMGInitStep << "	";
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
    int intervalStepForce = GlobalDataBase::GetIntParaFromDB("intervalStepForce");

    //string header = "# iter		CL		Cd		Side		Cd_p		Cd_f		Cd-CL2par		Xcp		cfx		cfy		cfz		Cmx		Cmy		Cmz		Hinge";

    int wordwidth = 16;
    oss << setiosflags(ios::left);
    oss << setiosflags(ios::scientific);
    oss << setprecision(10);

	//if (outnstep / intervalStepForce == 1)
	//{
	//	oss << "# Aerodynamic Force of " << partName << endl;
	//	oss << header << endl;
	//}
	//int nMGLevel = GlobalDataBase::GetIntParaFromDB("nMGLevel");
	//int mgInitStep = PHSPACE::GlobalDataBase::GetIntParaFromDB( "mgInitStep" );
	//int totalMGInitStep = 0;
	//if (nMGLevel > 1)
	//{
	//	for (int iLevel = 1; iLevel < nMGLevel; iLevel ++)
	//	{
	//		totalMGInitStep += iLevel * mgInitStep;
	//	}
	//}

	//oss << outnstep + totalMGInitStep << "	";
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

        /*string header = "# iter             CL                  Cd               CL/Cd                Side                Cd_p                Cd_f\
           Cd-CL2par                 Xcp                 Ycp                 Zcp                 cfx                 cfy                 cfz                 Cmx                 Cmy                 Cmz            Walltime";
        oss << header << endl;*/
    }

    DumpForce(outnstep,oss);

    WriteASCIIFile(file, oss.str());

    ParallelCloseFile(actkey);

    //! For part force.
    bool isGlobalCoordinate = false;
    for (int iPart = 1; iPart < numberofTotalPart; ++ iPart)
    {
        int bcID = partBCID[iPart];

        BoundaryCondition *bc = GlobalBoundaryCondition::GetBC(bcID);
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
            /*string header = "# iter             CL                  Cd               CL/Cd                Side                Cd_p                Cd_f\
            Cd-CL2par                 Xcp                 Ycp                 Zcp                 cfx                 cfy                 cfz                 Cmx                 Cmy                 Cmz               Hinge";
            oss << "# Aerodynamic Force of " << bcName << endl;
            oss << header << endl;*/

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
    ossGlobalData << "#     forceRefenenceLength         :       " << this->forceRefenenceLength * reynoldsReferenceLengthDimensional << "  m" << endl;
    ossGlobalData << "#     forceRefenenceLengthSpanWise :       " << this->forceRefenenceLengthSpanWise * reynoldsReferenceLengthDimensional << "  m" << endl;

    if (GetDim() == THREE_D)
    {
        ossGlobalData << "#     forceRefenenceArea           :       " << this->forceRefenenceArea * reynoldsReferenceLengthDimensional * reynoldsReferenceLengthDimensional << "  m2" << endl;
    }
    else
    {
        ossGlobalData << "#     forceRefenenceArea           :       " << this->forceRefenenceArea * reynoldsReferenceLengthDimensional << "  m^2" << endl;
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
            BoundaryCondition *bc = GlobalBoundaryCondition::GetBC(bcID);
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
    RDouble TorqueRefX, TorqueRefY, TorqueRefZ, forceRefenenceAreaPart, forceRefenenceLengthPart, forceRefenenceLengthSpanWisePart;
    int bcID = partBCID[partID];
    BoundaryCondition *bc = GlobalBoundaryCondition::GetBC(bcID);
    Data_Param *bcParamDB = bc->GetBCParamDataBase();
    bcParamDB->GetData("forceRefenenceArea",   &forceRefenenceAreaPart,   PHDOUBLE, 1);
    bcParamDB->GetData("forceRefenenceLength", &forceRefenenceLengthPart, PHDOUBLE, 1);
    bcParamDB->GetData("forceRefenenceLengthSpanWise", &forceRefenenceLengthSpanWisePart, PHDOUBLE, 1);
    bcParamDB->GetData("TorqueRefX", &TorqueRefX, PHDOUBLE, 1);
    bcParamDB->GetData("TorqueRefY", &TorqueRefY, PHDOUBLE, 1);
    bcParamDB->GetData("TorqueRefZ", &TorqueRefZ, PHDOUBLE, 1);

    RDouble reynoldsReferenceLengthDimensional = GlobalDataBase::GetDoubleParaFromDB("reynoldsReferenceLengthDimensional");

    oss << resetiosflags(ios::scientific);

    oss <<endl;
    oss << "*"<<partName <<endl;
    oss << "#forceRefenenceLength         :    " << forceRefenenceLengthPart * reynoldsReferenceLengthDimensional << "  m" << endl;
    oss << "#forceRefenenceLengthSpanWise :    " << forceRefenenceLengthSpanWisePart * reynoldsReferenceLengthDimensional << "  m" <<endl;

    if (GetDim() == THREE_D)
    {
        oss << "#forceRefenenceArea           :    " << forceRefenenceAreaPart * reynoldsReferenceLengthDimensional * reynoldsReferenceLengthDimensional << "  m2" <<endl;
    }
    else
    {
        oss << "#forceRefenenceArea           :    " << forceRefenenceAreaPart * reynoldsReferenceLengthDimensional << "  m^2" <<endl;
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

void Force::DumpData(DataContainer * data)
{
    data->moveto_begin();

    PHWrite(data, numberofSolidWallPart);
    PHWrite(data, numberofTotalPart);

    PHWrite(data, partBCID, numberofTotalPart);

    PHWrite(data, cpx,			numberofTotalPart);
    PHWrite(data, cpy,			numberofTotalPart);
    PHWrite(data, cpz,			numberofTotalPart);

    PHWrite(data, cl_tot,		numberofTotalPart);
    PHWrite(data, cd_tot,		numberofTotalPart);
    PHWrite(data, cd_pr,		numberofTotalPart);
    PHWrite(data, cd_sf,		numberofTotalPart);
    PHWrite(data, cd_cl2pa,		numberofTotalPart);
    PHWrite(data, xcp,			numberofTotalPart);
    PHWrite(data, side,			numberofTotalPart);
    PHWrite(data, hingeMoment,	numberofTotalPart);

}

}
